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
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split

from scripts.layout import PipelinePaths, add_data_args

DEFAULT_DRUG = "paclitaxel"
DEFAULT_MULTI_SPLIT_COL = "split_ctrp"
DEFAULT_MASK_OBSM_KEY = "M_ctrp"

# Frozen split assignments live in the repository, not under the (gitignored) data root:
# they are small, they must be versioned with the results they define, and a split that
# cannot be cited cannot be defended.
SPLIT_DIR = Path(__file__).resolve().parents[2] / "splits"


def frozen_split(
    cell_lines: np.ndarray,
    split_file: Path,
    seed: int,
    regenerate: bool = False,
) -> dict[str, str]:
    """Return a ``cell line -> {train, val, test}`` assignment, frozen on disk.

    The split is stored rather than recomputed because its input is not stable. Eligible
    lines are those carrying at least one CTRP label, so a change to the drug panel or to
    ``ctrp_to_h5ad``'s filters silently moves lines between train, val and test -- and runs
    from either side of that change look comparable when they are not.

    On first use, or with ``regenerate=True``, the 70/15/15 split is drawn with ``seed``
    and written to ``split_file``. After that the file is authoritative. If the data
    contains an eligible line the file does not cover, this raises rather than guessing:
    a quietly reassigned line is precisely the failure being prevented.
    """
    if split_file.exists() and not regenerate:
        stored = pd.read_csv(split_file)
        assignment = dict(zip(stored["Cell_line"], stored["split"]))
        unknown = sorted(set(cell_lines) - assignment.keys())
        if unknown:
            shown = ", ".join(unknown[:10]) + (" ..." if len(unknown) > 10 else "")
            raise ValueError(
                f"{len(unknown)} eligible cell line(s) are absent from {split_file.name}: "
                f"{shown}\nThe frozen split cannot place them. Re-freeze deliberately "
                f"(--regenerate-split), and treat every earlier result as having been "
                f"scored on a different split."
            )
        absent = sorted(assignment.keys() - set(cell_lines))
        if absent:
            shown = ", ".join(absent[:10]) + (" ..." if len(absent) > 10 else "")
            print(f"  {len(absent)} frozen line(s) not eligible in this data: {shown}")
        print(f"  Using frozen split from {split_file} ({len(assignment)} lines).")
        return assignment

    train_lines, val_lines, test_lines = _split_cell_lines(np.asarray(cell_lines), seed=seed)
    assignment = {line: "train" for line in train_lines}
    assignment.update({line: "val" for line in val_lines})
    assignment.update({line: "test" for line in test_lines})

    split_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sorted(assignment.items()), columns=["Cell_line", "split"]).to_csv(
        split_file, index=False
    )
    print(f"  Wrote frozen split to {split_file} ({len(assignment)} lines, seed={seed}).")
    return assignment


def _split_cell_lines(
    cell_lines: np.ndarray, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (train_lines, val_lines, test_lines) for a 70/15/15 split."""
    train_lines, temp_lines = train_test_split(cell_lines, test_size=0.30, random_state=seed)
    val_lines, test_lines = train_test_split(temp_lines, test_size=0.50, random_state=seed)
    return train_lines, val_lines, test_lines


def _apply_assignment(adata, valid_cells_df, split_col: str, assignment: dict[str, str]) -> None:
    """Write a frozen line -> split assignment onto the cells of ``valid_cells_df``."""
    adata.obs.loc[valid_cells_df.index, split_col] = (
        valid_cells_df["Cell_line"].map(assignment).to_numpy()
    )
    present = pd.Series(assignment)
    present = present[present.index.isin(valid_cells_df["Cell_line"].unique())]
    counts = present.value_counts().to_dict()
    print(
        "Cell Line Split -> "
        + ", ".join(f"{name}: {counts.get(name, 0)}" for name in ("train", "val", "test"))
    )


def run(
    h5ad_path: str,
    target_drug: str = DEFAULT_DRUG,
    seed: int = 42,
    regenerate: bool = False,
):
    """Per-drug cell-line-grouped 70/15/15 train/val/test split."""
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    mask_col = f"train_mask_{target_drug}"
    split_col = f"split_{target_drug}"

    adata.obs[split_col] = "unassigned"

    valid_cells_df = adata.obs[adata.obs[mask_col] == True]  # noqa: E712 - nullable bool col
    unique_cell_lines = valid_cells_df["Cell_line"].unique()
    print(f"Found {len(unique_cell_lines)} unique cell lines with {target_drug} labels.")

    assignment = frozen_split(
        unique_cell_lines, SPLIT_DIR / f"{split_col}.csv", seed=seed, regenerate=regenerate
    )
    _apply_assignment(adata, valid_cells_df, split_col, assignment)

    print(f"\nFinal Cell Split distribution for {target_drug}:")
    print(adata.obs[split_col].value_counts())

    _save(adata, h5ad_path)
    return adata


def run_multi(
    h5ad_path: str,
    seed: int = 42,
    split_col: str = DEFAULT_MULTI_SPLIT_COL,
    mask_obsm_key: str = DEFAULT_MASK_OBSM_KEY,
    regenerate: bool = False,
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

    assignment = frozen_split(
        unique_cell_lines, SPLIT_DIR / f"{split_col}.csv", seed=seed, regenerate=regenerate
    )
    _apply_assignment(adata, valid_cells_df, split_col, assignment)

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
    add_data_args(parser)
    parser.add_argument(
        "--path",
        type=Path,
        default=None,
        help="Targets h5ad (default: <variant>/..._with_targets.h5ad).",
    )
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
    parser.add_argument(
        "--regenerate-split",
        action="store_true",
        help=(
            f"Redraw the split and overwrite the frozen file in {SPLIT_DIR.name}/. "
            "Every result produced before it was regenerated was scored on different "
            "held-out lines and is no longer comparable."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant, args.score)
    h5ad_path = str(args.path or paths.targets_h5ad)
    if args.mode == "multi":
        run_multi(
            h5ad_path=h5ad_path,
            seed=args.seed,
            split_col=args.split_col,
            mask_obsm_key=args.mask_obsm_key,
            regenerate=args.regenerate_split,
        )
    else:
        run(h5ad_path, args.drug, args.seed, regenerate=args.regenerate_split)
