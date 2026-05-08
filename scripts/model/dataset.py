import torch
from torch.utils.data import Dataset
import scanpy as sc

class ScGPTDrugDataset(Dataset):
    def __init__(self, h5ad_path, target_drug="paclitaxel", use_rep="X_scGPT", split="train"):
        """
        split: Should be 'train', 'val', or 'test'
        """
        self.split = split
        print(f"Loading {split} split from {h5ad_path}...")
        adata = sc.read_h5ad(h5ad_path)

        # 1. Filter using the exact split column we just created
        split_col = f'split_{target_drug}'
        if split_col not in adata.obs.columns:
            raise ValueError(f"Split column '{split_col}' not found! Run the split generation script first.")

        valid_adata = adata[adata.obs[split_col] == split].copy()
        print(f"Loaded {valid_adata.n_obs} cells for the '{split}' set.")

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