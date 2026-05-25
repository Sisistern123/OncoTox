import argparse

import anndata as ad
import scanpy as sc

DEFAULT_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"


def run(h5ad_path: str = DEFAULT_PATH, force: bool = False):
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
    # Replicates the exact prep used for the standard UMAPs.
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
    parser.add_argument("--path", default=DEFAULT_PATH)
    parser.add_argument("--force", action="store_true", help="Recompute X_pca even if it exists.")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run(args.path, args.force)
