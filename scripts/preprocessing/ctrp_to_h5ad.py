"""Map CTRPv2 drug-response targets onto cells in the embedded SCP542 AnnData.

The target score is selected with ``score`` (CLI: ``--score``). In all cases a
*higher* value means a *less* sensitive (more resistant) cell line:

* ``auc_z``   -- per-drug z-scored ``auc`` (**retired as the default 27.07.2026**). Standardizing within each
  drug puts every multi-task head on the same scale, so the MSE of a
  well-covered drug does not dominate the loss, and it removes the per-drug
  potency offset the model cannot infer from expression anyway.
* ``auc``     -- ``area_under_curve / conc_pts_fit`` from the post-QC sigmoid
  fits (``v20.data.curves_post_qc.txt``). CTRP's raw AUC is an *integrated*
  area, so it scales with how many concentration points were fitted (8-29,
  usually 16); dividing by ``conc_pts_fit`` puts drugs on a common axis.
* ``mean_pv`` -- legacy score: the unweighted mean of ``cpd_avg_pv`` over the
  dose grid (``v20.data.per_cpd_post_qc.txt``). Kept so the pre-AUC results
  stay reproducible; correlates ~0.97 with raw AUC but ignores the curve fit.

Two outputs are written into the AnnData (both happen by default):

1. Multi-drug target matrix (used by the multi-task training):
    * ``adata.obsm["Y_ctrp"]``  : float32 (n_cells, K), NaN where missing.
    * ``adata.obsm["M_ctrp"]``  : bool    (n_cells, K), True where Y_ctrp is observed.
    * ``adata.uns["ctrp_drugs"]``: list[str] of length K (normalized drug names),
      giving the column order of Y_ctrp / M_ctrp.
    * ``adata.uns["ctrp_score"]``: which score Y_ctrp holds.
    * ``adata.uns["ctrp_score_center"]`` / ``["ctrp_score_scale"]``: per-drug mean
      and std of the pre-z-score ``auc`` (length K, aligned with ``ctrp_drugs``),
      so predictions can be mapped back to AUC units. Both are 0/1 filler when
      ``score != "auc_z"``.

   Drugs are kept only if at least ``min_cell_lines`` distinct SCP542-overlapping
   cell lines were screened against them (default 50) so we don't add heads with
   too little support.

2. Per-drug flat columns (back-compat with the original single-drug pipeline):
    * ``adata.obs["viability_<drug>"]``    : per-cell target, NaN when missing.
    * ``adata.obs["train_mask_<drug>"]``   : per-cell bool, True when present.

   The ``viability_`` prefix is historical -- the column holds whichever score
   was selected, not necessarily a percent viability. Controlled by
   ``extra_single_drug_cols``. Defaults to ("paclitaxel",) so the 25.05.2026
   baseline remains reproducible without any code changes.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from scripts.preprocessing.layout import (
    CTRP_SCORES,
    DEFAULT_CTRP_SCORE,
    PipelinePaths,
    add_data_args,
)

DEFAULT_MIN_CELL_LINES = 50
DEFAULT_EXTRA_SINGLE_DRUG_COLS: tuple[str, ...] = ("paclitaxel",)


def _normalize_cell_line(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().str.replace("-", "")


def _normalize_drug(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower()


def _load_score_values(ctrp_dir: Path, score: str) -> pd.DataFrame:
    """Return one ``score`` value per (experiment_id, master_cpd_id).

    ``auc``/``auc_z`` read the post-QC sigmoid fits and normalize the integrated
    area by the number of fitted concentration points; ``mean_pv`` averages the
    raw dose grid (legacy behaviour).
    """
    if score == "mean_pv":
        per_cpd = pd.read_csv(
            ctrp_dir / "v20.data.per_cpd_post_qc.txt",
            sep="\t",
            usecols=["experiment_id", "master_cpd_id", "cpd_avg_pv"],
        )
        return (
            per_cpd.groupby(["experiment_id", "master_cpd_id"], as_index=False)["cpd_avg_pv"]
            .mean()
            .rename(columns={"cpd_avg_pv": "score"})
        )

    curves = pd.read_csv(
        ctrp_dir / "v20.data.curves_post_qc.txt",
        sep="\t",
        usecols=["experiment_id", "master_cpd_id", "conc_pts_fit", "area_under_curve"],
    )
    n_bad = int(curves["area_under_curve"].isna().sum())
    if n_bad:
        print(f"  Dropping {n_bad} curve fits with no area_under_curve.")
        curves = curves.dropna(subset=["area_under_curve"])
    # CTRP's area_under_curve is integrated, not averaged, so it grows with the
    # size of the concentration grid (conc_pts_fit is 8-29, usually 16).
    curves["score"] = curves["area_under_curve"] / curves["conc_pts_fit"]
    return curves[["experiment_id", "master_cpd_id", "score"]]


def _load_ctrp_long(ctrp_dir: Path, score: str) -> pd.DataFrame:
    """Return the merged CTRPv2 long table with normalized name columns."""
    ctrp_values = _load_score_values(ctrp_dir, score)
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
        ``score`` (one row per (cell line, drug) inside the SCP542 overlap;
        replicate experiments are averaged).
    kept_drugs : ordered list of drug names (normalized) to use as Y_ctrp columns.
    """
    long_overlap = (
        ctrp_full[ctrp_full["ccl_name_norm"].isin(overlap_cell_lines_norm)]
        .groupby(["ccl_name_norm", "cpd_name_norm"], as_index=False)["score"]
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


def _zscore_per_drug(
    long_overlap: pd.DataFrame, kept_drugs: list[str]
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Standardize ``score`` within each drug, over the overlapping cell lines.

    Statistics are computed per cell line (not per cell), so a cell line
    contributing many cells does not pull the mean toward itself. A drug whose
    scores have zero spread would blow up, so its scale is left at 1.0 --
    ``min_cell_lines`` plus the learnability filter should remove those anyway.

    Returns the standardized table plus the per-drug ``center`` and ``scale``
    arrays (aligned with ``kept_drugs``) needed to map predictions back to AUC.
    """
    grouped = long_overlap.groupby("cpd_name_norm")["score"]
    center = grouped.mean().reindex(kept_drugs)
    scale = grouped.std(ddof=0).reindex(kept_drugs)

    degenerate = scale[(scale.isna()) | (scale <= 0)].index.tolist()
    if degenerate:
        print(
            f"  Warning: {len(degenerate)} drugs have zero AUC spread across cell lines; "
            f"leaving their scale at 1.0 (e.g. {degenerate[:5]})."
        )
    scale = scale.where(scale > 0, 1.0)

    out = long_overlap.copy()
    out["score"] = (
        out["score"] - out["cpd_name_norm"].map(center)
    ) / out["cpd_name_norm"].map(scale)
    return out, center.to_numpy(dtype=np.float32), scale.to_numpy(dtype=np.float32)


def run(
    input_h5ad: str,
    output_h5ad: str,
    ctrp_dir: str,
    min_cell_lines: int = DEFAULT_MIN_CELL_LINES,
    target_drugs: Sequence[str] | None = None,
    extra_single_drug_cols: Sequence[str] = DEFAULT_EXTRA_SINGLE_DRUG_COLS,
    score: str = DEFAULT_CTRP_SCORE,
):
    """Map CTRPv2 drug-response scores onto cells in the embedded AnnData.

    Parameters
    ----------
    score:
        Which CTRPv2 response score to use as the target: ``auc`` (default since 27.07.2026),
        ``auc``, or the legacy ``mean_pv``. See the module docstring.
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
    if score not in CTRP_SCORES:
        raise ValueError(f"score must be one of {CTRP_SCORES}, got {score!r}")

    print("Loading AnnData...")
    adata = sc.read_h5ad(input_h5ad)

    print(f"Loading CTRPv2 metadata (score={score})...")
    ctrp_full = _load_ctrp_long(Path(ctrp_dir), score)

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

    if score == "auc_z":
        print("Z-scoring AUC within each drug (over overlapping cell lines)...")
        long_overlap, center, scale = _zscore_per_drug(long_overlap, kept_drugs)
    else:
        center = np.zeros(len(kept_drugs), dtype=np.float32)
        scale = np.ones(len(kept_drugs), dtype=np.float32)

    print(f"Building (cell line x drug) {score} matrix...")
    cl_drug_matrix = long_overlap.pivot(
        index="ccl_name_norm", columns="cpd_name_norm", values="score"
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
    adata.uns["ctrp_score"] = score
    adata.uns["ctrp_score_center"] = center
    adata.uns["ctrp_score_scale"] = scale

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
    parser = argparse.ArgumentParser(description="Merge CTRPv2 response targets into embedded AnnData.")
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
    paths = PipelinePaths.build(args.data_root, args.variant, args.score)
    min_cell_lines = 0 if args.all_drugs else args.min_cell_lines
    run(
        input_h5ad=str(args.input or paths.embed_h5ad),
        output_h5ad=str(args.output or paths.targets_h5ad),
        ctrp_dir=str(args.ctrp_dir or paths.ctrp_dir),
        min_cell_lines=min_cell_lines,
        target_drugs=args.drugs,
        extra_single_drug_cols=tuple(args.single_drug_cols),
        score=args.score,
    )
