import argparse
from pathlib import Path

import anndata as ad
import scanpy as sc

from scripts.preprocessing.layout import PipelinePaths, add_data_args


def run(h5ad_path: str, force: bool = False):
    """Permanently save a PCA baseline (X_pca) into the AnnData file.

    If X_pca is already present, the calculation is skipped unless `force=True`.
    """
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    if "X_pca" in adata.obsm and not force:
        print("X_pca already exists! Pass force=True (or --force) to recompute.")
        return adata

    if "X_pca" in adata.obsm and force:
        print("Force flag set; deleting existing X_pca before recompute...")
        del adata.obsm["X_pca"]

    print("Calculating PCA baseline...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)

    print("Saving updated AnnData with X_pca...")
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
        help="Targets h5ad (default: <variant>/..._with_targets.h5ad).",
    )
    parser.add_argument("--force", action="store_true", help="Recompute X_pca even if it exists.")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant)
    run(str(args.path or paths.targets_h5ad), args.force)
