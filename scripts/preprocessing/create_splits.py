import anndata as ad
import scanpy as sc
import pandas as pd
from sklearn.model_selection import train_test_split

# --- 1. Load the Data ---
file_path = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
print(f"Loading {file_path}...")
adata = sc.read_h5ad(file_path)

TARGET_DRUG = "paclitaxel"
mask_col = f'train_mask_{TARGET_DRUG}'
split_col = f'split_{TARGET_DRUG}'

# Initialize the split column with 'unassigned'
adata.obs[split_col] = 'unassigned'

# --- 2. Isolate the Valid Cells & Find@ Unique Cell Lines ---
# Get a dataframe of only the cells that actually have drug labels
valid_cells_df = adata.obs[adata.obs[mask_col] == True]

# The 'Cell_line' column contains the cell line names (e.g., 'MCF7_BREAST')
# Let's extract the unique cell lines present in our valid data
unique_cell_lines = valid_cells_df['Cell_line'].unique()
print(f"Found {len(unique_cell_lines)} unique cell lines with {TARGET_DRUG} labels.")

# --- 3. Perform the Grouped 70/15/15 Split on CELL LINES ---
SEED = 42

# Step A: Split off 30% of cell lines for Val/Test, keep 70% for Train
train_lines, temp_lines = train_test_split(unique_cell_lines, test_size=0.30, random_state=SEED)

# Step B: Split that 30% perfectly in half to get 15% Val and 15% Test
val_lines, test_lines = train_test_split(temp_lines, test_size=0.50, random_state=SEED)

print(f"Cell Line Split -> Train: {len(train_lines)}, Val: {len(val_lines)}, Test: {len(test_lines)}")

# --- 4. Assign the Labels Back to the Individual Cells ---
# Find the indices (cell barcodes) that belong to the cell lines in each split
train_idx = valid_cells_df[valid_cells_df['Cell_line'].isin(train_lines)].index
val_idx = valid_cells_df[valid_cells_df['Cell_line'].isin(val_lines)].index
test_idx = valid_cells_df[valid_cells_df['Cell_line'].isin(test_lines)].index

# Safely assign the split labels back to the main AnnData object
adata.obs.loc[train_idx, split_col] = 'train'
adata.obs.loc[val_idx, split_col] = 'val'
adata.obs.loc[test_idx, split_col] = 'test'

print(f"\nFinal Cell Split distribution for {TARGET_DRUG}:")
print(adata.obs[split_col].value_counts())

# --- 5. Save the Updated AnnData ---
print("\nSaving updated AnnData...")
adata.obs.index = adata.obs.index.astype(str).astype(object) # standard index sanitization
ad.settings.allow_write_nullable_strings = True
adata.write_h5ad(file_path, convert_strings_to_categoricals=False)
print("Done! Leakage-free grouped splits are permanently saved.")