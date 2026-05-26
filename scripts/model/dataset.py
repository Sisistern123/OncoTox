import numpy as np
import scanpy as sc
import torch
from torch.utils.data import Dataset


class ScGPTDrugDataset(Dataset):
    """Single-drug dataset (used by ``train_baseline.py`` / ``train_scGPT.py``)."""

    def __init__(self, h5ad_path, target_drug="paclitaxel", use_rep="X_scGPT", split="train"):
        """
        split: Should be 'train', 'val', or 'test'
        """
        self.split = split
        print(f"Loading {split} split from {h5ad_path}...")
        adata = sc.read_h5ad(h5ad_path)

        split_col = f'split_{target_drug}'
        if split_col not in adata.obs.columns:
            raise ValueError(f"Split column '{split_col}' not found! Run the split generation script first.")

        split_indices = adata.obs[adata.obs[split_col] == split].index

        if len(split_indices) == 0:
            raise ValueError(f"No cells found for split '{split}' in {h5ad_path}")

        valid_adata = adata[split_indices].copy()
        print(f"Loaded {valid_adata.n_obs} cells for the '{split}' set.")

        if use_rep in valid_adata.obsm.keys():
            self.X = torch.tensor(valid_adata.obsm[use_rep], dtype=torch.float32)
        else:
            raise ValueError(f"Representation '{use_rep}' not found in adata.obsm! Please verify your embeddings.")

        target_col = f'viability_{target_drug}'
        self.y = torch.tensor(valid_adata.obs[target_col].values, dtype=torch.float32).unsqueeze(1)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


class MultiDrugDataset(Dataset):
    """Multi-task dataset over all CTRP drugs persisted by ``ctrp_to_h5ad.py``.

    Returns ``(x, y_vec, mask_vec)`` triples:
      * ``x``        : (D,)   float32 cell embedding (PCA or scGPT).
      * ``y_vec``    : (K,)   float32 viability per drug (zeros where missing).
      * ``mask_vec`` : (K,)   float32 1.0 where observed, 0.0 where missing.

    ``self.drug_names`` exposes the column order of ``y_vec`` / ``mask_vec`` so
    the training loop can report per-drug metrics.
    """

    def __init__(
        self,
        h5ad_path,
        use_rep: str = "X_scGPT",
        split: str = "train",
        split_col: str = "split_ctrp",
        y_obsm_key: str = "Y_ctrp",
        mask_obsm_key: str = "M_ctrp",
        drugs_uns_key: str = "ctrp_drugs",
        drugs: "list[str] | None" = None,
    ):
        self.split = split
        print(f"Loading multi-drug {split} split from {h5ad_path}...")
        adata = sc.read_h5ad(h5ad_path)

        if split_col not in adata.obs.columns:
            raise ValueError(
                f"Split column '{split_col}' not found in adata.obs. "
                f"Run create_splits.run_multi first."
            )
        for key in (y_obsm_key, mask_obsm_key):
            if key not in adata.obsm:
                raise ValueError(
                    f"obsm['{key}'] not found. Run ctrp_to_h5ad first."
                )
        if drugs_uns_key not in adata.uns:
            raise ValueError(
                f"uns['{drugs_uns_key}'] not found. Run ctrp_to_h5ad first."
            )

        split_mask = (adata.obs[split_col] == split).to_numpy()
        n_split = int(split_mask.sum())
        if n_split == 0:
            raise ValueError(f"No cells found for split '{split}' in {h5ad_path}")

        if use_rep not in adata.obsm:
            raise ValueError(
                f"Representation '{use_rep}' not found in adata.obsm! "
                f"Please verify your embeddings."
            )

        all_drugs: list[str] = list(adata.uns[drugs_uns_key])
        if drugs is None:
            keep_idx = np.arange(len(all_drugs))
            drug_names = all_drugs
        else:
            requested = [d.strip().lower() for d in drugs]
            available = {d: i for i, d in enumerate(all_drugs)}
            missing = [d for d in requested if d not in available]
            if missing:
                raise ValueError(
                    f"Requested drugs not present in obsm['{y_obsm_key}']: {missing}. "
                    f"Re-run ctrp_to_h5ad with a lower --min-cell-lines or pass the "
                    f"drug name as it appears in uns['{drugs_uns_key}']."
                )
            keep_idx = np.array([available[d] for d in requested], dtype=int)
            drug_names = requested

        X = np.asarray(adata.obsm[use_rep], dtype=np.float32)[split_mask]
        Y = np.asarray(adata.obsm[y_obsm_key], dtype=np.float32)[split_mask][:, keep_idx]
        M = np.asarray(adata.obsm[mask_obsm_key], dtype=bool)[split_mask][:, keep_idx]
        # Replace NaNs in masked-out entries with 0.0 so they're safe to feed
        # through PyTorch even though the loss will mask them out.
        Y = np.where(M, Y, 0.0).astype(np.float32, copy=False)

        self.X = torch.from_numpy(X)
        self.y = torch.from_numpy(Y)
        self.mask = torch.from_numpy(M.astype(np.float32, copy=False))
        self.drug_names: list[str] = drug_names

        n_per_cell = M.sum(axis=1)
        n_per_drug = M.sum(axis=0)
        print(
            f"Loaded {n_split} cells for the '{split}' set. "
            f"K={len(self.drug_names)} drugs | "
            f"mean drugs/cell={n_per_cell.mean():.1f} | "
            f"min drug coverage={int(n_per_drug.min())} cells | "
            f"max drug coverage={int(n_per_drug.max())} cells."
        )

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.mask[idx]
