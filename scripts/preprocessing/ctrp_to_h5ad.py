"""Map CTRPv2 viability targets onto cells in the embedded SCP542 AnnData.

Two outputs are written into the AnnData (both happen by default):

1. Multi-drug target matrix (used by the multi-task training):
    * ``adata.obsm["Y_ctrp"]``  : float32 (n_cells, K), NaN where missing.
    * ``adata.obsm["M_ctrp"]``  : bool    (n_cells, K), True where Y_ctrp is observed.
    * ``adata.uns["ctrp_drugs"]``: list[str] of length K (normalized drug names),
      giving the column order of Y_ctrp / M_ctrp.

   Drugs are kept only if at least ``min_cell_lines`` distinct SCP542-overlapping
   cell lines were screened against them (default 50) so we don't add heads with
   too little support.

2. Per-drug flat columns (back-compat with the original single-drug pipeline):
    * ``adata.obs["viability_<drug>"]``    : per-cell viability, NaN when missing.
    * ``adata.obs["train_mask_<drug>"]``   : per-cell bool, True when present.

   Controlled by ``extra_single_drug_cols``. Defaults to ("paclitaxel",) so the
   25.05.2026 baseline remains reproducible without any code changes.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from scripts.preprocessing.layout import PipelinePaths, add_data_args

DEFAULT_MIN_CELL_LINES = 50
DEFAULT_EXTRA_SINGLE_DRUG_COLS: tuple[str, ...] = ("paclitaxel",)


def _normalize_cell_line(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().str.replace("-", "")


def _normalize_drug(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower()


def _load_ctrp_long(ctrp_dir: Path) -> pd.DataFrame:
    """Return the merged CTRPv2 long table with normalized name columns."""
    ctrp_values = pd.read_csv(
        ctrp_dir / "v20.data.per_cpd_post_qc.txt",
        sep="\t",
        usecols=["experiment_id", "master_cpd_id", "cpd_avg_pv"],
    )
    ctrp_exp_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_experiment.txt",
        sep="\t",
        usecols=["experiment_id", "master_ccl_id"],
    )
    ctrp_cell_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_cell_line.txt",
        sep="\t",
        usecols=["master_ccl_id", "ccl_name"],
    )
    ctrp_cpd_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_compound.txt",
        sep="\t",
        usecols=["master_cpd_id", "cpd_name"],
    )

    ctrp_full = (
        ctrp_values.merge(ctrp_exp_meta, on="experiment_id", how="inner")
        .merge(ctrp_cell_meta, on="master_ccl_id", how="inner")
        .merge(ctrp_cpd_meta, on="master_cpd_id", how="inner")
    )
    ctrp_full["ccl_name_norm"] = _normalize_cell_line(ctrp_full["ccl_name"])
    ctrp_full["cpd_name_norm"] = _normalize_drug(ctrp_full["cpd_name"])
    return ctrp_full


def _build_drug_table(
    ctrp_full: pd.DataFrame,
    overlap_cell_lines_norm: set[str],
    min_cell_lines: int,
    target_drugs: Sequence[str] | None,
) -> tuple[pd.DataFrame, list[str]]:
    """Aggregate to one row per (cell line, drug) and decide which drugs to keep.

    Returns
    -------
    long_overlap : DataFrame with columns ``ccl_name_norm``, ``cpd_name_norm``,
        ``cpd_avg_pv`` (one row per (cell line, drug) inside the SCP542 overlap).
    kept_drugs : ordered list of drug names (normalized) to use as Y_ctrp columns.
    """
    long_overlap = (
        ctrp_full[ctrp_full["ccl_name_norm"].isin(overlap_cell_lines_norm)]
        .groupby(["ccl_name_norm", "cpd_name_norm"], as_index=False)["cpd_avg_pv"]
        .mean()
    )

    if target_drugs is not None:
        target_drugs_norm = [d.strip().lower() for d in target_drugs]
        kept_drugs = [d for d in target_drugs_norm if d in set(long_overlap["cpd_name_norm"])]
        missing = sorted(set(target_drugs_norm) - set(kept_drugs))
        if missing:
            print(
                f"  Warning: {len(missing)} requested drugs have no overlap data "
                f"and will be skipped (e.g. {missing[:5]})."
            )
    else:
        coverage = long_overlap.groupby("cpd_name_norm")["ccl_name_norm"].nunique()
        kept = coverage[coverage >= min_cell_lines].sort_values(ascending=False)
        kept_drugs = kept.index.tolist()
        print(
            f"  Drug filter: {len(kept_drugs)} / {coverage.shape[0]} drugs kept "
            f"(>= {min_cell_lines} overlapping cell lines)."
        )

    long_overlap = long_overlap[long_overlap["cpd_name_norm"].isin(set(kept_drugs))]
    return long_overlap, kept_drugs


def run(
    input_h5ad: str,
    output_h5ad: str,
    ctrp_dir: str,
    min_cell_lines: int = DEFAULT_MIN_CELL_LINES,
    target_drugs: Sequence[str] | None = None,
    extra_single_drug_cols: Sequence[str] = DEFAULT_EXTRA_SINGLE_DRUG_COLS,
):
    """Map CTRPv2 viability scores onto cells in the embedded AnnData.

    Parameters
    ----------
    target_drugs:
        If provided, restrict the multi-drug matrix to these drug names (after
        lower-casing). When ``None`` (default), include every CTRPv2 drug that
        passes ``min_cell_lines``.
    min_cell_lines:
        Drug filter. A drug must have been screened against at least this many
        SCP542-overlapping cell lines. Ignored if ``target_drugs`` is given.
    extra_single_drug_cols:
        For each drug name in this iterable also write the legacy flat columns
        ``viability_<drug>`` / ``train_mask_<drug>`` so the original single-drug
        training scripts continue to work. Pass ``()`` to disable.
    """
    print("Loading AnnData...")
    adata = sc.read_h5ad(input_h5ad)

    print("Loading CTRPv2 metadata...")
    ctrp_full = _load_ctrp_long(Path(ctrp_dir))

    cell_line_norm = (
        adata.obs["Cell_line"]
        .astype(str)
        .str.split("_")
        .str[0]
        .pipe(_normalize_cell_line)
    )
    overlap_cell_lines_norm = set(cell_line_norm.unique()) & set(ctrp_full["ccl_name_norm"])
    print(
        f"  Overlap with SCP542: {len(overlap_cell_lines_norm)} cell lines "
        f"out of {cell_line_norm.nunique()} in AnnData."
    )

    long_overlap, kept_drugs = _build_drug_table(
        ctrp_full,
        overlap_cell_lines_norm=overlap_cell_lines_norm,
        min_cell_lines=min_cell_lines,
        target_drugs=target_drugs,
    )

    print("Building (cell line x drug) viability matrix...")
    cl_drug_matrix = long_overlap.pivot(
        index="ccl_name_norm", columns="cpd_name_norm", values="cpd_avg_pv"
    )
    # Reindex columns so the ordering matches uns["ctrp_drugs"] exactly.
    cl_drug_matrix = cl_drug_matrix.reindex(columns=kept_drugs)

    print(f"Mapping {len(kept_drugs)} drugs to {adata.n_obs} single cells...")
    Y_full = cl_drug_matrix.reindex(cell_line_norm.values)
    Y = Y_full.to_numpy(dtype=np.float32)
    M = ~np.isnan(Y)

    adata.obsm["Y_ctrp"] = Y
    adata.obsm["M_ctrp"] = M.astype(bool)
    adata.uns["ctrp_drugs"] = list(kept_drugs)

    has_any_label = M.any(axis=1)
    print(
        f"  Multi-drug summary: {has_any_label.sum()} / {adata.n_obs} cells have "
        f"at least one CTRP label; mean drugs/cell = {M.sum(axis=1).mean():.1f}."
    )

    if extra_single_drug_cols:
        print(f"Writing legacy per-drug columns for: {list(extra_single_drug_cols)}")
        drug_to_idx = {d: i for i, d in enumerate(kept_drugs)}
        for raw_drug in extra_single_drug_cols:
            drug = raw_drug.strip().lower()
            target_col = f"viability_{drug}"
            mask_col = f"train_mask_{drug}"
            if drug not in drug_to_idx:
                print(
                    f"  Warning: '{drug}' is not in the kept drug list "
                    f"(min_cell_lines={min_cell_lines}); skipping legacy cols."
                )
                # Write empty columns so downstream code can still run without KeyError.
                adata.obs[target_col] = np.nan
                adata.obs[mask_col] = False
                continue
            col_idx = drug_to_idx[drug]
            adata.obs[target_col] = Y[:, col_idx]
            adata.obs[mask_col] = M[:, col_idx]
            n_present = int(M[:, col_idx].sum())
            print(f"  {drug}: {n_present} / {adata.n_obs} cells with viability.")

    print("Sanitizing metadata to fix H5AD string/index compatibility...")
    adata.obs.index = adata.obs.index.astype(str).astype(object)
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col]) or pd.api.types.is_object_dtype(
            adata.obs[col]
        ):
            adata.obs[col] = adata.obs[col].astype(str).astype(object)

    adata.var.index = adata.var.index.astype(str).astype(object)
    for col in adata.var.columns:
        if pd.api.types.is_string_dtype(adata.var[col]) or pd.api.types.is_object_dtype(
            adata.var[col]
        ):
            adata.var[col] = adata.var[col].astype(str).astype(object)

    print(f"Saving updated AnnData to {output_h5ad}...")
    adata.write_h5ad(output_h5ad, convert_strings_to_categoricals=False)
    print("Done!")
    return adata


def _parse_args():
    parser = argparse.ArgumentParser(description="Merge CTRPv2 viability targets into embedded AnnData.")
    add_data_args(parser)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Embedding h5ad (default: <variant>/SCP542_CCLE_scGPT_human_embeddings.h5ad).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output h5ad (default: <variant>/..._with_targets.h5ad).",
    )
    parser.add_argument(
        "--ctrp-dir",
        type=Path,
        default=None,
        help="CTRPv2 directory (default: <data-root>/metadata/CTRPv2...).",
    )
    parser.add_argument(
        "--min-cell-lines",
        type=int,
        default=DEFAULT_MIN_CELL_LINES,
        help="Drug filter: minimum number of SCP542-overlapping cell lines required.",
    )
    parser.add_argument(
        "--all-drugs",
        action="store_true",
        help="Shortcut for --min-cell-lines 0 (keep every CTRPv2 drug with any overlap).",
    )
    parser.add_argument(
        "--drugs",
        nargs="+",
        default=None,
        help="Optional explicit drug list (overrides --min-cell-lines / --all-drugs).",
    )
    parser.add_argument(
        "--single-drug-cols",
        nargs="+",
        default=list(DEFAULT_EXTRA_SINGLE_DRUG_COLS),
        help="Drugs to also expose as legacy viability_<drug>/train_mask_<drug> columns.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant)
    min_cell_lines = 0 if args.all_drugs else args.min_cell_lines
    run(
        input_h5ad=str(args.input or paths.embed_h5ad),
        output_h5ad=str(args.output or paths.targets_h5ad),
        ctrp_dir=str(args.ctrp_dir or paths.ctrp_dir),
        min_cell_lines=min_cell_lines,
        target_drugs=args.drugs,
        extra_single_drug_cols=tuple(args.single_drug_cols),
    )
