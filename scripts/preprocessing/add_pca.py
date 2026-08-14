import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA

from scripts.preprocessing.expression import kinker_transform
from scripts.layout import PipelinePaths, add_data_args


# 512 to MATCH scGPT's embedding width, so the two arms are compared at equal dimensionality.
# It is a comparability choice, NOT a statement about how many components this atlas needs --
# no scree or variance-explained criterion selected it, and it has never been varied
# (docs/TODO.md item 4A; `uns['pca_fits']['variance_ratio']`, written by _pca_record below, is
# what would answer the question). Noted here 14.08.2026: the justification existed only in
# report/results_numbers.tex, where nobody reading this module would find it.
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
# `split_paclitaxel` was removed on 12.08.2026 with the rest of the single-drug chain: it cost a
# second 512-component train-only fit on every run, for a column whose only reader
# (`ScGPTDrugDataset`) had no caller. h5ads written before that date still carry
# `X_pca_train_paclitaxel`; nothing reads it.
TRAIN_SPLIT_COLS = ("split_ctrp",)


def train_pca_key(split_col: str) -> str:
    """``obsm`` key holding the PCA fitted on that split's training cells only."""
    return "X_pca_train_" + split_col.removeprefix("split_")


def _pca_record(
    genes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    center: np.ndarray,
    components: np.ndarray,
    variance: np.ndarray,
    variance_ratio: np.ndarray,
    n_comps: int,
    seed: int,
    fitted_on: str,
) -> dict:
    """Everything needed to reproduce a stored projection, or apply it to new cells.

    ``obsm["X_pca*"]`` holds only where each cell landed. That is the *output* of the
    projection and does not contain the projection itself: it says nothing about which
    genes built each component, nor how much of the total variance the kept components
    account for. Both are computed on the way and were discarded until 10.08.2026
    (audit item 04b) -- so neither "what fraction of the variance does PCA(512) retain?"
    nor "which genes dominate PC1?" could be answered without re-running the fit, and a
    new cell could not be placed in the same space at all.

    That last point is the one that matters for ``X_pca_train_*``: those are fitted on
    training cells only, so re-running reproduces them only while every input is
    bit-identical (same cells, same HVG list, same split, same seed). Across a
    re-preprocessing sweep that stops holding, and the recipe behind the stored
    coordinates would be gone.

    Stored in ``uns`` rather than ``varm`` because ``varm`` is indexed by the file's own
    ``var``, and the two axes do not match: the ``scgpt`` step drops out-of-vocabulary genes
    from the targets file's ``.X`` while PCA runs on the full HVG set from the convert output.
    How many genes that leaves changes with the gene-symbol repair, which is exactly why the
    record carries ``genes`` rather than relying on the file's own axis -- it labels the gene
    axis of ``mean``/``std``/``center`` and the columns of ``components``.

    The transform this reproduces, for a cell's ``log2(1 + CPM/10)`` gene vector ``x``::

        z = clip((x - mean) / std, -scale_max_value, +scale_max_value)
        coords = (z - center) @ components.T

    ``mean``/``std`` are the per-gene standardization; ``center`` is the separate mean
    that the PCA itself subtracts, which is near but not exactly zero because clipping
    happens after standardization.
    """
    return {
        "genes": np.asarray(genes, dtype=object).astype(str),
        "mean": np.asarray(mean, dtype=np.float32),
        "std": np.asarray(std, dtype=np.float32),
        "center": np.asarray(center, dtype=np.float32),
        "scale_max_value": float(SCALE_MAX_VALUE),
        "components": np.asarray(components, dtype=np.float32),
        "variance": np.asarray(variance, dtype=np.float32),
        "variance_ratio": np.asarray(variance_ratio, dtype=np.float32),
        "n_comps": int(n_comps),
        "seed": int(seed),
        "fitted_on": str(fitted_on),
    }


def _pca_fitted_on_train(
    X: np.ndarray,
    train_mask: np.ndarray,
    n_comps: int,
    seed: int,
    genes: np.ndarray,
    fitted_on: str,
    max_value: float = SCALE_MAX_VALUE,
) -> tuple[np.ndarray, dict]:
    """Standardize and project, with every statistic estimated on training cells only.

    ``X`` is the transformed expression matrix (``log2(1 + CPM/10)``) for *all* cells,
    before per-gene scaling; ``train_mask`` selects the
    training ones. Per-gene mean and standard deviation **and** the principal-component
    rotation are fitted on ``X[train_mask]`` alone and then applied to every cell, so a
    held-out cell's coordinates depend on no held-out cell.

    This is the hand-rolled equivalent of ``sc.pp.scale`` + ``sc.pp.pca`` because neither
    separates fitting from transforming, which is precisely what a train-only fit needs.
    Returns the projection for all cells in the row order of ``X``, and the ``_pca_record``
    describing the fit that produced it.
    """
    train = X[train_mask]
    mean = train.mean(axis=0)
    # ddof=1, matching sc.pp.scale on the all-cells path (scanpy applies the Bessel
    # correction; numpy's .std() defaults to ddof=0). Harmonized 10.08.2026, audit item 04b:
    # the two PCA fits are meant to differ only in which cells they see, and on this atlas
    # the correction itself is worth well under 0.01%.
    std = train.std(axis=0, ddof=1)
    # Genes with no variance across training cells carry no information and would divide
    # by zero. Leaving std at 1 keeps them mean-centred but unscaled, so a held-out cell
    # that does express such a gene is not silently amplified.
    std[std == 0] = 1.0

    Z = (X - mean) / std
    np.clip(Z, -max_value, max_value, out=Z)

    pca = PCA(n_components=n_comps, random_state=seed)
    pca.fit(Z[train_mask])
    record = _pca_record(
        genes=genes,
        mean=mean,
        std=std,
        # pca.mean_ is the column mean of Z[train_mask] -- the centring the PCA does itself,
        # on top of the standardization above. Not zero, because the clip lands after the
        # divide by std.
        center=pca.mean_,
        components=pca.components_,
        variance=pca.explained_variance_,
        # Relative to the variance of the TRAINING cells, since that is the only thing this
        # fit saw. Not comparable as a percentage with the all-cells key below.
        variance_ratio=pca.explained_variance_ratio_,
        n_comps=n_comps,
        seed=seed,
        fitted_on=fitted_on,
    )
    return pca.transform(Z).astype(np.float32), record


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

    Each key recomputed here also gets an entry in ``uns["pca_fits"]`` holding the fit
    itself -- loadings, variances and standardization statistics -- see ``_pca_record``.

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
    # Records are keyed by the obsm key they describe, and only ever written for keys
    # actually recomputed in this call: a key skipped because it already exists keeps
    # whatever record it came with, rather than acquiring one from a different fit.
    pca_fits = dict(adata.uns.get("pca_fits", {}))
    for key in todo:
        pca_fits.pop(key, None)

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
    # float32 (Selin, 12.08.2026). Chosen for the per-fold CV fits, where the counts matrix is
    # 53,513 x 5,000 and float64 doubles peak memory to ~4.4 GB, and extended here so the
    # descriptive all-cells fit and the per-fold fits stay identical in kind -- the property
    # she established on 10.08.2026 when she harmonized `ddof` to 1, that the two differ only in
    # which cells they see. It is a preprocessing change: `X_pca` moves in its last digits.
    X_log = (src.X.toarray() if sparse.issparse(src.X) else np.asarray(src.X)).astype(
        np.float32, copy=False
    )
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
        adata.obsm[key], pca_fits[key] = _pca_fitted_on_train(
            X_log,
            train_mask,
            n_comps,
            seed,
            genes=np.asarray(src.var_names),
            fitted_on=f"{n_train} cells labelled 'train' in obs['{col}']",
        )
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
        # sc.pp.scale writes the statistics it used into var["mean"]/var["std"]; taking them
        # from there rather than recomputing means the record cannot disagree with the
        # transform it describes. (It uses ddof=1, where numpy's .std() defaults to ddof=0 --
        # recomputing by hand got this wrong and the synthetic check in the docstring caught it.)
        Z_all = src.X.toarray() if sparse.issparse(src.X) else np.asarray(src.X)
        sc.pp.pca(src, n_comps=n_comps, random_state=seed)
        adata.obsm["X_pca"] = src.obsm["X_pca"]
        pca_fits["X_pca"] = _pca_record(
            genes=np.asarray(src.var_names),
            mean=src.var["mean"].to_numpy(),
            std=src.var["std"].to_numpy(),
            center=Z_all.mean(axis=0),
            # scanpy stores loadings as (n_genes, n_comps); sklearn uses the transpose. Both
            # records are written in sklearn's orientation so one formula reprojects either.
            components=np.asarray(src.varm["PCs"]).T,
            variance=src.uns["pca"]["variance"],
            variance_ratio=src.uns["pca"]["variance_ratio"],
            n_comps=n_comps,
            seed=seed,
            fitted_on=f"all {src.n_obs} cells",
        )
        print(f"  X_pca (all cells) computed on {src.n_vars} genes -> shape {adata.obsm['X_pca'].shape}.")
        print(f"    {100 * float(pca_fits['X_pca']['variance_ratio'].sum()):.1f}% of variance retained.")

    adata.uns["pca_fits"] = pca_fits

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
