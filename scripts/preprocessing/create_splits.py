import argparse

import anndata as ad
import scanpy as sc
from sklearn.model_selection import train_test_split

DEFAULT_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
DEFAULT_DRUG = "paclitaxel"


def run(h5ad_path: str = DEFAULT_PATH, target_drug: str = DEFAULT_DRUG, seed: int = 42):
    """Create a cell-line-grouped 70/15/15 train/val/test split.

    Grouping by cell line prevents leakage from cells of the same line appearing
    in multiple splits.
    """
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    mask_col = f"train_mask_{target_drug}"
    split_col = f"split_{target_drug}"

    adata.obs[split_col] = "unassigned"

    valid_cells_df = adata.obs[adata.obs[mask_col] == True]  # noqa: E712 - obs col is nullable boolean
    unique_cell_lines = valid_cells_df["Cell_line"].unique()
    print(f"Found {len(unique_cell_lines)} unique cell lines with {target_drug} labels.")

    # 70/30, then split the 30 in half -> 70/15/15.
    train_lines, temp_lines = train_test_split(unique_cell_lines, test_size=0.30, random_state=seed)
    val_lines, test_lines = train_test_split(temp_lines, test_size=0.50, random_state=seed)

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

    print("\nSaving updated AnnData...")
    adata.obs.index = adata.obs.index.astype(str).astype(object)
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(h5ad_path, convert_strings_to_categoricals=False)
    print("Done! Leakage-free grouped splits are permanently saved.")
    return adata


def _parse_args():
    parser = argparse.ArgumentParser(description="Create cell-line-grouped train/val/test splits.")
    parser.add_argument("--path", default=DEFAULT_PATH)
    parser.add_argument("--drug", default=DEFAULT_DRUG)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run(args.path, args.drug, args.seed)
