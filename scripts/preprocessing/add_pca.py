import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

from scripts.preprocessing.layout import PipelinePaths, add_data_args


def run(h5ad_path: str, force: bool = False, counts_h5ad: str | None = None):
    """Compute the PCA baseline (``X_pca``) and store it in the targets AnnData.

    PCA is computed on the **HVG-filtered counts** taken from ``counts_h5ad`` (the
    ``convert`` output ``SCP542_CCLE.h5ad``), *not* on the targets file's own ``.X``.

    Why: the ``scgpt`` step drops scGPT-out-of-vocabulary genes from ``.X`` (e.g.
    5,000 HVG -> 4,576). Running PCA on that matrix would silently couple the PCA
    baseline to scGPT's vocabulary and shrink its gene set. Sourcing the counts from
    the convert output keeps PCA on the single HVG filter, so PCA and scGPT are a
    clean, like-for-like comparison. The targets file's ``.X`` is left untouched.

    If ``counts_h5ad`` is None, PCA falls back to a copy of the targets ``.X``
    (legacy behaviour). If ``X_pca`` is already present, it is skipped unless
    ``force=True``.
    """
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    if "X_pca" in adata.obsm and not force:
        print("X_pca already exists! Pass force=True (or --force) to recompute.")
        return adata

    adata.obsm.pop("X_pca", None)

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

    sc.pp.normalize_total(src, target_sum=1e4)
    sc.pp.log1p(src)
    sc.pp.pca(src)
    adata.obsm["X_pca"] = src.obsm["X_pca"]
    print(f"  X_pca computed on {src.n_vars} genes -> shape {adata.obsm['X_pca'].shape}.")

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
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant)
    counts = args.counts or paths.raw_h5ad
    run(str(args.path or paths.targets_h5ad), args.force, counts_h5ad=str(counts))
