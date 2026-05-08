import pandas as pd
import anndata as ad
import scanpy as sc

# ---------------------------------------------------------
# 1. Load the Expression Matrix
# ---------------------------------------------------------
print("Loading expression matrix... (this may take a few minutes and require high RAM)")
# Note: Swap 'CPM_data.txt' with 'UMIcount_data.txt' if you want raw pre-QC counts
df_expr = pd.read_csv('/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/expression/CPM_data.txt', sep='\t', index_col=0)

# Transpose the dataframe so Cells are rows and Genes are columns
adata = ad.AnnData(X=df_expr.T)

# ---------------------------------------------------------
# 2. Load the Metadata
# ---------------------------------------------------------
print("Loading metadata...")
df_meta = pd.read_csv('/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/Metadata.txt', sep='\t', low_memory=False)

# Drop the "TYPE" row (which is the first row of data at index 0)
df_meta = df_meta.drop(0)

# Set the cell ID column ('NAME') as the index so it aligns with the expression data
df_meta = df_meta.set_index('NAME')

# ---------------------------------------------------------
# 3. Combine and Align
# ---------------------------------------------------------
print("Aligning metadata with expression data...")
# Only keep metadata for cells that actually exist in the expression matrix
# and align them in the exact same order
adata.obs = df_meta.loc[adata.obs_names]

# (Optional) If you want to load the tSNE coordinates for a specific cancer type:
# tsne_df = pd.read_csv('tSNE_Breast_Cancer.txt', sep='\t').drop(0).set_index('NAME')
# Note: You would only apply this if your AnnData object was subsetted to just breast cancer cells.

# ---------------------------------------------------------
# 4. Save to .h5ad
# ---------------------------------------------------------
print("Saving to h5ad...")
ad.settings.allow_write_nullable_strings = True
adata.write('/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE.h5ad')

print(f"Success! Created AnnData object: {adata}")