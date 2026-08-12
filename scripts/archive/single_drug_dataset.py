"""``ScGPTDrugDataset`` -- the single-drug dataset, archived 12.08.2026 (Selin).

**Not importable from here and not runnable.** It was lifted out of ``scripts/model/dataset.py``
unchanged; see ``scripts/archive/README.md`` for the rule this directory follows.

**Why it went.** Nothing called it. Its docstring named ``train_baseline.py`` and ``train_scGPT.py``
as its callers, and both were deleted on 26.05.2026 when ``train_multitask.py --drugs <one>``
replaced them -- K=1 reduces exactly to plain MSE
(``docs/steps/corrections-and-dead-ends.md``). It outlived them by two and a half months as the
last link of a chain nothing consumed: ``ctrp_to_h5ad`` wrote ``viability_<drug>`` /
``train_mask_<drug>``, ``create_splits.run`` turned those into ``split_<drug>``, and ``add_pca``
fitted a second 512-component train-only decomposition for that column on every pipeline run --
all of it to feed this class, which no notebook or script instantiated.

**What went with it**, on the same decision: ``targets()`` no longer writes the legacy per-drug
columns, ``splits()`` writes only ``split_ctrp``, and ``add_pca.TRAIN_SPLIT_COLS`` is down to
``("split_ctrp",)``.

**It cannot be run as it stands**, independently of the above: it reads ``obs["viability_<drug>"]``
and ``obs["split_<drug>"]``, neither of which the pipeline writes any more, and the h5ads that do
carry them were built on ``mean_pv`` -- a target removed on 11.08.2026. The Step 04 results it
produced are void and explicitly unregenerable
(``docs/steps/04-single-task-results.md``).

Restoring single-drug work does not mean restoring this file: the multi-task path already covers it
with ``MultiDrugDataset(..., drugs=["<one>"])``, which reads ``Y_ctrp``/``M_ctrp`` and needs no
per-drug columns at all.
"""

import scanpy as sc
import torch
from torch.utils.data import Dataset


class ScGPTDrugDataset(Dataset):
    """Single-drug dataset (used by ``train_baseline.py`` / ``train_scGPT.py``, both deleted)."""

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
