"""Cell-line-grouped train/val/test splits.

Two entry points:

* ``run``        : per-drug split. Reads ``train_mask_<drug>`` and writes
  ``split_<drug>`` -- preserved for back-compat with the original single-drug
  pipeline (paclitaxel).
* ``run_multi``  : drug-agnostic split for the multi-task setting. Uses
  ``adata.obsm["M_ctrp"]`` (any-drug-observed) to decide which cell lines are
  eligible, and writes a single ``split_ctrp`` column. Because the split is
  grouped by cell line, it's leakage-free for every drug head simultaneously.
"""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import scanpy as sc
from sklearn.model_selection import train_test_split

DEFAULT_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
DEFAULT_DRUG = "paclitaxel"
DEFAULT_MULTI_SPLIT_COL = "split_ctrp"
DEFAULT_MASK_OBSM_KEY = "M_ctrp"


def _split_cell_lines(
    cell_lines: np.ndarray, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (train_lines, val_lines, test_lines) for a 70/15/15 split."""
    train_lines, temp_lines = train_test_split(cell_lines, test_size=0.30, random_state=seed)
    val_lines, test_lines = train_test_split(temp_lines, test_size=0.50, random_state=seed)
    return train_lines, val_lines, test_lines


def run(h5ad_path: str = DEFAULT_PATH, target_drug: str = DEFAULT_DRUG, seed: int = 42):
    """Per-drug cell-line-grouped 70/15/15 train/val/test split."""
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    mask_col = f"train_mask_{target_drug}"
    split_col = f"split_{target_drug}"

    adata.obs[split_col] = "unassigned"

    valid_cells_df = adata.obs[adata.obs[mask_col] == True]  # noqa: E712 - nullable bool col
    unique_cell_lines = valid_cells_df["Cell_line"].unique()
    print(f"Found {len(unique_cell_lines)} unique cell lines with {target_drug} labels.")

    train_lines, val_lines, test_lines = _split_cell_lines(unique_cell_lines, seed=seed)

    print(
        f"Cell Line Split -> Train: {len(train_lines)}, "
        f"Val: {len(val_lines)}, Test: {len(test_lines)}"
    )

    train_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(train_lines)].index
    val_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(val_lines)].index
    test_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(test_lines)].index

    adata.obs.loc[train_idx, split_col] = "train"
    adata.obs.loc[val_idx, split_col] = "val"
    adata.obs.loc[test_idx, split_col] = "test"

    print(f"\nFinal Cell Split distribution for {target_drug}:")
    print(adata.obs[split_col].value_counts())

    _save(adata, h5ad_path)
    return adata


def run_multi(
    h5ad_path: str = DEFAULT_PATH,
    seed: int = 42,
    split_col: str = DEFAULT_MULTI_SPLIT_COL,
    mask_obsm_key: str = DEFAULT_MASK_OBSM_KEY,
):
    """Drug-agnostic cell-line-grouped split using the CTRP mask matrix.

    A cell line is eligible if any of its cells has at least one observed drug
    label in ``adata.obsm[mask_obsm_key]``. The split is identical across drug
    heads, which keeps val/test untouched as new heads are added later.
    """
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    if mask_obsm_key not in adata.obsm:
        raise ValueError(
            f"obsm['{mask_obsm_key}'] not found. Run ctrp_to_h5ad first so the "
            f"multi-drug mask is available."
        )

    M = np.asarray(adata.obsm[mask_obsm_key], dtype=bool)
    has_any_label = M.any(axis=1)

    adata.obs[split_col] = "unassigned"
    valid_cells_df = adata.obs.loc[has_any_label]
    unique_cell_lines = valid_cells_df["Cell_line"].unique()
    print(
        f"Found {len(unique_cell_lines)} unique cell lines with at least one CTRP "
        f"drug label across {valid_cells_df.shape[0]} cells."
    )

    train_lines, val_lines, test_lines = _split_cell_lines(unique_cell_lines, seed=seed)

    print(
        f"Cell Line Split -> Train: {len(train_lines)}, "
        f"Val: {len(val_lines)}, Test: {len(test_lines)}"
    )

    train_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(train_lines)].index
    val_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(val_lines)].index
    test_idx = valid_cells_df[valid_cells_df["Cell_line"].isin(test_lines)].index

    adata.obs.loc[train_idx, split_col] = "train"
    adata.obs.loc[val_idx, split_col] = "val"
    adata.obs.loc[test_idx, split_col] = "test"

    print(f"\nFinal Cell Split distribution for multi-drug ({split_col}):")
    print(adata.obs[split_col].value_counts())

    _save(adata, h5ad_path)
    return adata


def _save(adata, h5ad_path: str) -> None:
    print("\nSaving updated AnnData...")
    adata.obs.index = adata.obs.index.astype(str).astype(object)
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(h5ad_path, convert_strings_to_categoricals=False)
    print("Done! Leakage-free grouped splits are permanently saved.")


def _parse_args():
    parser = argparse.ArgumentParser(description="Create cell-line-grouped train/val/test splits.")
    parser.add_argument("--path", default=DEFAULT_PATH)
    parser.add_argument(
        "--mode",
        choices=("single", "multi"),
        default="single",
        help="single: per-drug split (split_<drug>); multi: drug-agnostic split (split_ctrp).",
    )
    parser.add_argument("--drug", default=DEFAULT_DRUG, help="Used only when --mode single.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--split-col", default=DEFAULT_MULTI_SPLIT_COL, help="Used only when --mode multi.")
    parser.add_argument("--mask-obsm-key", default=DEFAULT_MASK_OBSM_KEY, help="Used only when --mode multi.")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    if args.mode == "multi":
        run_multi(
            h5ad_path=args.path,
            seed=args.seed,
            split_col=args.split_col,
            mask_obsm_key=args.mask_obsm_key,
        )
    else:
        run(args.path, args.drug, args.seed)
