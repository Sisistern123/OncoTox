import scanpy as sc
import pandas as pd
import numpy as np

# --- 1. Load your single-cell data ---
print("Loading AnnData...")
adata = sc.read_h5ad('/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings.h5ad')

# --- 2. Load CTRPv2 Data ---
print("Loading CTRPv2 metadata...")
# 2a. Load the actual viability values (cpd_avg_pv)
ctrp_values = pd.read_csv(
    "/Users/selin/Desktop/OncoTox/data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.per_cpd_post_qc.txt",
    sep="\t",
    usecols=["experiment_id", "master_cpd_id", "cpd_avg_pv"]
)

# 2b. Load experiment metadata (links experiment to cell line)
ctrp_exp_meta = pd.read_csv(
    "/Users/selin/Desktop/OncoTox/data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt",
    sep="\t",
    usecols=["experiment_id", "master_ccl_id"]
)

# 2c. Load cell line metadata (gets the actual cell line name)
ctrp_cell_meta = pd.read_csv(
    "/Users/selin/Desktop/OncoTox/data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt",
    sep="\t",
    usecols=["master_ccl_id", "ccl_name"]
)

# 2d. Load compound metadata (gets the actual drug name)
ctrp_cpd_meta = pd.read_csv(
    "/Users/selin/Desktop/OncoTox/data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt",
    sep="\t",
    usecols=["master_cpd_id", "cpd_name"]
)

# --- 3. Merge CTRPv2 into one clean table ---
print("Merging CTRPv2 tables...")
ctrp_full = (
    ctrp_values.merge(ctrp_exp_meta, on="experiment_id", how="inner")
            .merge(ctrp_cell_meta, on="master_ccl_id", how="inner")
            .merge(ctrp_cpd_meta, on="master_cpd_id", how="inner")
)

# Normalize cell line names (lowercase, no spaces/hyphens) to ensure perfect matching
ctrp_full['ccl_name_norm'] = ctrp_full['ccl_name'].astype(str).str.strip().str.lower().str.replace("-", "")

# --- 4. Choose your target drug ---
TARGET_DRUG = "paclitaxel"  # Make sure this matches CTRP's naming (lowercase)
ctrp_full['cpd_name_norm'] = ctrp_full['cpd_name'].astype(str).str.strip().str.lower()

# Filter for only the drug we care about
drug_data = ctrp_full[ctrp_full['cpd_name_norm'] == TARGET_DRUG]

# If there are multiple experiments for the same cell line & drug, take the mean viability
cell_line_to_viability = drug_data.groupby('ccl_name_norm')['cpd_avg_pv'].mean().to_dict()

# --- 5. Map scores to the single cells in AnnData ---
print(f"Mapping {TARGET_DRUG} scores to single cells...")

# Normalize the AnnData cell line names the exact same way
adata.obs['Cell_line_norm'] = adata.obs['Cell_line'].astype(str).str.split("_").str[0].str.strip().str.lower().str.replace("-", "")

# Map the dictionary to a new column in adata.obs
target_col_name = f'viability_{TARGET_DRUG}'
adata.obs[target_col_name] = adata.obs['Cell_line_norm'].map(cell_line_to_viability)

# --- 6. Create a Training Mask ---
# Not all cell lines were tested with Paclitaxel. Cells with NaN should be ignored during training.
mask_col_name = f'train_mask_{TARGET_DRUG}'
adata.obs[mask_col_name] = adata.obs[target_col_name].notna()

print(f"Summary for {TARGET_DRUG}:")
print(f"Total cells: {adata.n_obs}")
print(f"Cells with a valid viability score: {adata.obs[mask_col_name].sum()}")

# Clean up the temporary normalization column
adata.obs = adata.obs.drop(columns=['Cell_line_norm'])

# --- 7. Save the updated AnnData ---
output_path = '/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad'

print("Sanitizing metadata to fix H5AD string/index compatibility...")

# 1. Sanitize the cell metadata (.obs)
adata.obs.index = adata.obs.index.astype(str).astype(object)
for col in adata.obs.columns:
    if pd.api.types.is_string_dtype(adata.obs[col]) or pd.api.types.is_object_dtype(adata.obs[col]):
        adata.obs[col] = adata.obs[col].astype(str).astype(object)

# 2. Sanitize the gene metadata (.var)
adata.var.index = adata.var.index.astype(str).astype(object)
for col in adata.var.columns:
    if pd.api.types.is_string_dtype(adata.var[col]) or pd.api.types.is_object_dtype(adata.var[col]):
        adata.var[col] = adata.var[col].astype(str).astype(object)

print(f"Saving updated AnnData to {output_path}...")
adata.write_h5ad(output_path, convert_strings_to_categoricals=False)
print("Done!")