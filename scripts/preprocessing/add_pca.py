import scanpy as sc
import anndata as ad
# permanently saves PCA to .h5ad file (for UMAP visualization and training)

file_path = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"

print(f"Loading {file_path}...")
adata = sc.read_h5ad(file_path)

if "X_pca" not in adata.obsm:
    print("Calculating PCA baseline...")
    # Replicating the exact prep you used for your standard UMAPs
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)

    print("Saving updated AnnData with X_pca...")
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(file_path, convert_strings_to_categoricals=False)
    print("Done! You can now run baseline training.")
else:
    print("X_pca already exists! You are good to go.")