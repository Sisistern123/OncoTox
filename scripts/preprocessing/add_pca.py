import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA

from scripts.preprocessing.expression import kinker_transform
from scripts.preprocessing.layout import PipelinePaths, add_data_args


DEFAULT_N_COMPS = 512
# Project-wide seed. scanpy's own default is random_state=0, so this must be passed
# explicitly to sc.pp.pca -- leaving it out silently uses 0.
DEFAULT_SEED = 42
# z-score clipping cap, following Seurat's ScaleData(scale.max = 10) default. Shared by both
# the all-cells and the train-fitted path so the two cannot drift apart.
SCALE_MAX_VALUE = 10.0

# Fixed-split columns (written by create_splits.py) that get their own train-fitted PCA.
# The 5-fold CV is deliberately not covered: its folds are drawn at training time, so a
# train-only fit there cannot be a single stored matrix. See docs/steps/02.
TRAIN_SPLIT_COLS = ("split_ctrp", "split_paclitaxel")


def train_pca_key(split_col: str) -> str:
    """``obsm`` key holding the PCA fitted on that split's training cells only."""
    return "X_pca_train_" + split_col.removeprefix("split_")


def _pca_fitted_on_train(
    X: np.ndarray,
    train_mask: np.ndarray,
    n_comps: int,
    seed: int,
    max_value: float = SCALE_MAX_VALUE,
) -> np.ndarray:
    """Standardize and project, with every statistic estimated on training cells only.

    ``X`` is the transformed expression matrix (``log2(1 + CPM/10)``) for *all* cells,
    before per-gene scaling; ``train_mask`` selects the
    training ones. Per-gene mean and standard deviation **and** the principal-component
    rotation are fitted on ``X[train_mask]`` alone and then applied to every cell, so a
    held-out cell's coordinates depend on no held-out cell.

    This is the hand-rolled equivalent of ``sc.pp.scale`` + ``sc.pp.pca`` because neither
    separates fitting from transforming, which is precisely what a train-only fit needs.
    Returns the projection for all cells, in the row order of ``X``.
    """
    train = X[train_mask]
    mean = train.mean(axis=0)
    std = train.std(axis=0)
    # Genes with no variance across training cells carry no information and would divide
    # by zero. Leaving std at 1 keeps them mean-centred but unscaled, so a held-out cell
    # that does express such a gene is not silently amplified.
    std[std == 0] = 1.0

    Z = (X - mean) / std
    np.clip(Z, -max_value, max_value, out=Z)

    pca = PCA(n_components=n_comps, random_state=seed)
    pca.fit(Z[train_mask])
    return pca.transform(Z).astype(np.float32)


def run(
    h5ad_path: str,
    force: bool = False,
    counts_h5ad: str | None = None,
    n_comps: int = DEFAULT_N_COMPS,
    seed: int = DEFAULT_SEED,
):
    """Compute the PCA baseline (``X_pca``) and store it in the targets AnnData.

    ``n_comps`` PCA components are kept (default 512, matching the scGPT embedding
    width so the PCA-vs-scGPT comparison uses the same input dimensionality).

    PCA is computed on the **HVG-filtered counts** taken from ``counts_h5ad`` (the
    ``convert`` output ``SCP542_CCLE.h5ad``), *not* on the targets file's own ``.X``.

    Why: the ``scgpt`` step drops scGPT-out-of-vocabulary genes from ``.X`` (e.g.
    5,000 HVG -> 4,576). Running PCA on that matrix would silently couple the PCA
    baseline to scGPT's vocabulary and shrink its gene set. Sourcing the counts from
    the convert output keeps PCA on the single HVG filter, so PCA and scGPT are a
    clean, like-for-like comparison. The targets file's ``.X`` is left untouched.

    The transform applied before PCA is ``log2(1 + CPM/10)`` then per-gene scaling (see the
    comments at the call site). The input is already CPM, so no library-size normalization
    is applied here -- that step belongs on the full gene set and has already happened.

    **Two kinds of key are written.** ``X_pca`` is fitted on all cells and is the
    descriptive representation -- UMAPs and latent-space validation, where holding cells
    out would be wrong. In addition, for every fixed-split column in ``TRAIN_SPLIT_COLS``
    present in ``obs``, a ``X_pca_train_<split>`` key is fitted on that split's training
    cells only and is what a model evaluated on that split should read.

    If ``counts_h5ad`` is None, PCA falls back to a copy of the targets ``.X``
    (legacy behaviour). Keys already present are skipped unless ``force=True``.
    """
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    split_cols = []
    for col in TRAIN_SPLIT_COLS:
        if col in adata.obs:
            split_cols.append(col)
        else:
            print(f"  obs['{col}'] not present -- skipping {train_pca_key(col)}.")

    wanted = ["X_pca"] + [train_pca_key(col) for col in split_cols]
    todo = wanted if force else [key for key in wanted if key not in adata.obsm]
    if not todo:
        print(f"Already present: {', '.join(wanted)}. Pass force=True (or --force) to recompute.")
        return adata

    for key in todo:
        adata.obsm.pop(key, None)

    if counts_h5ad is not None:
        print(f"Computing PCA on HVG-filtered counts from {counts_h5ad}...")
        src = sc.read_h5ad(counts_h5ad)
        if not np.array_equal(np.asarray(adata.obs_names), np.asarray(src.obs_names)):
            raise ValueError(
                "Cell order/identity mismatch between targets and counts h5ad; "
                "cannot align X_pca. Re-run convert/scgpt for this variant."
            )
    else:
        print("No counts file given; computing PCA on a copy of the targets .X (legacy).")
        src = adata.copy()

    # SCP542 ships as CPM (Kinker et al. distribute `CPM_data.txt`), so library-size
    # normalization has already been applied -- once, to the full gene set, which is where the
    # standard recipe puts it. What remains is the log transform, and we use the one the
    # dataset's own authors use: log2(1 + CPM/10). Same call as the HVG step, so genes are
    # selected and projected under an identical transform. See expression.py for the citation.
    #
    # Until 05.08.2026 a `normalize_total(target_sum=1e4)` sat here. It was not a rescale to a
    # different target but a *second* library-size normalization computed over the HVG subset
    # only: each cell divided by its own retained-gene sum, so a cell whose expression sits
    # mostly outside the HVG set got inflated relative to one whose expression sits inside it,
    # and the same cell was scaled differently in hvg1000 than in hvg5000.
    kinker_transform(src)
    max_comps = min(src.n_obs, src.n_vars)
    if n_comps > max_comps:
        raise ValueError(
            f"n_comps={n_comps} exceeds min(n_obs, n_vars)={max_comps} for {counts_h5ad}."
        )

    # The train-fitted keys are computed FIRST, and read src.X directly: sc.pp.scale below
    # standardizes src.X in place, so afterwards it is no longer the transformed matrix.
    X_log = src.X.toarray() if sparse.issparse(src.X) else np.asarray(src.X)
    for col in split_cols:
        key = train_pca_key(col)
        if key not in todo:
            continue
        train_mask = adata.obs[col].to_numpy() == "train"
        n_train = int(train_mask.sum())
        if n_train == 0:
            raise ValueError(f"obs['{col}'] labels no cell 'train'; cannot fit {key}.")
        if n_comps > min(n_train, src.n_vars):
            raise ValueError(
                f"n_comps={n_comps} exceeds min(n_train={n_train}, n_genes={src.n_vars}) "
                f"for obs['{col}']."
            )
        adata.obsm[key] = _pca_fitted_on_train(X_log, train_mask, n_comps, seed)
        print(f"  {key}: fitted on {n_train} train cells -> shape {adata.obsm[key].shape}.")

    if "X_pca" in todo:
        # Per-gene standardization across cells -- the axis CPM does not touch. Without it the
        # leading PCs track absolute expression level rather than variation.
        #
        # Every statistic here is fitted on ALL cells, held-out lines included, and that is
        # intended: X_pca is the descriptive representation (UMAPs, latent-space validation),
        # where excluding cells would be wrong. Models trained on a fixed split read the
        # train-fitted key above instead. The 5-fold CV still reads X_pca -- see docs/steps/02.
        sc.pp.scale(src, max_value=SCALE_MAX_VALUE)
        sc.pp.pca(src, n_comps=n_comps, random_state=seed)
        adata.obsm["X_pca"] = src.obsm["X_pca"]
        print(f"  X_pca (all cells) computed on {src.n_vars} genes -> shape {adata.obsm['X_pca'].shape}.")

    print("Saving updated AnnData with X_pca (targets .X unchanged)...")
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(h5ad_path, convert_strings_to_categoricals=False)
    print("Done! You can now run baseline training.")
    return adata


def _parse_args():
    parser = argparse.ArgumentParser(description="Add PCA baseline embedding to AnnData file.")
    add_data_args(parser)
    parser.add_argument(
        "--path",
        type=Path,
        default=None,
        help="Targets h5ad to write X_pca into (default: <variant>/..._with_targets.h5ad).",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=None,
        help="Counts h5ad to compute PCA from (default: <variant>/SCP542_CCLE.h5ad, the HVG-filtered convert output).",
    )
    parser.add_argument("--force", action="store_true", help="Recompute X_pca even if it exists.")
    parser.add_argument(
        "--n-comps",
        type=int,
        default=DEFAULT_N_COMPS,
        help=f"Number of PCA components to keep (default: {DEFAULT_N_COMPS}, matches scGPT width).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"random_state for sc.pp.pca (default: {DEFAULT_SEED}; scanpy's own default is 0).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant, args.score)
    counts = args.counts or paths.raw_h5ad
    run(
        str(args.path or paths.targets_h5ad),
        args.force,
        counts_h5ad=str(counts),
        n_comps=args.n_comps,
        seed=args.seed,
    )
