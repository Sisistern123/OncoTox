import torch
from torch.utils.data import Dataset, DataLoader
import scanpy as sc

class ScGPTDrugDataset(Dataset):
    def __init__(self, h5ad_path, target_drug="paclitaxel", use_rep="X_scGPT"):
        """
        Loads the h5ad file, filters for cells with valid drug labels,
        and extracts the specified embeddings and targets.
        """
        print(f"Loading AnnData from {h5ad_path}...")
        adata = sc.read_h5ad(h5ad_path)

        # 1. Apply the mask to keep only cells with valid target scores
        mask_col = f'train_mask_{target_drug}'
        valid_adata = adata[adata.obs[mask_col] == True].copy()
        print(f"Filtered down to {valid_adata.n_obs} cells with valid {target_drug} labels.")

        # 2. Extract Features (X)
        if use_rep in valid_adata.obsm.keys():
            self.X = torch.tensor(valid_adata.obsm[use_rep], dtype=torch.float32)
        else:
            raise ValueError(f"Representation '{use_rep}' not found in adata.obsm!")

        # 3. Extract Targets (y)
        target_col = f'viability_{target_drug}'
        self.y = torch.tensor(valid_adata.obs[target_col].values, dtype=torch.float32).unsqueeze(1)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

# --- Execution ---
file_path = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"

# Instantiate the dataset using your scGPT embeddings
scgpt_dataset = ScGPTDrugDataset(h5ad_path=file_path, target_drug="paclitaxel", use_rep="X_scGPT")

# Create the DataLoader (using a reasonable batch size for single-cell data)
train_loader = DataLoader(scgpt_dataset, batch_size=128, shuffle=True)

# Quick check to ensure the tensors are shaped correctly
features, targets = next(iter(train_loader))
print(f"Batch X shape: {features.shape}")
print(f"Batch y shape: {targets.shape}")