import anndata as ad
import scanpy as sc
import numpy as np
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

# --- 2. Isolate the Valid Cells ---
# We only want to split the cells that actually have drug labels
valid_indices = adata.obs[adata.obs[mask_col] == True].index.tolist()

# --- 3. Perform the 70/15/15 Split ---
SEED = 42

# Step A: Split off 30% for Val/Test, keep 70% for Train
train_idx, temp_idx = train_test_split(valid_indices, test_size=0.30, random_state=SEED)

# Step B: Split that 30% perfectly in half to get 15% Val and 15% Test
val_idx, test_idx = train_test_split(temp_idx, test_size=0.50, random_state=SEED)

# --- 4. Assign the Labels ---
# .loc is the safest way to update pandas/anndata dataframe columns
adata.obs.loc[train_idx, split_col] = 'train'
adata.obs.loc[val_idx, split_col] = 'val'
adata.obs.loc[test_idx, split_col] = 'test'

print(f"Split distribution for {TARGET_DRUG}:")
print(adata.obs[split_col].value_counts())

# --- 5. Save the Updated AnnData ---
print("Saving updated AnnData...")
adata.obs.index = adata.obs.index.astype(str).astype(object) # standard index sanitization

# Add this line to bypass the string writing restriction
ad.settings.allow_write_nullable_strings = True

adata.write_h5ad(file_path, convert_strings_to_categoricals=False)
print("Done! Splits are permanently saved.")