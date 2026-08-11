"""Map CTRPv2 drug-response targets onto cells in the embedded SCP542 AnnData.

**Source (since 11.08.2026): DrEval's reprocessed CTRPv2**, not CTRPv2's own 2015 distribution.
DrEval take CTRPv2's raw dose-response measurements, normalise each replicate against its no-drug
control, and re-fit every curve with CurveCurator -- including replicate variability in the fit rather
than averaging replicates before fitting, "as is the standard practice", which they argue "leads to
inaccurate or misleading drug response measures in the case of large discrepancies between
replicates" (``papers/DrEval_s41467-026-72903-w.pdf``, Methods, "Benchmark data"). The table is
pinned to one Zenodo record and fetched by ``fetch_ctrp_response``.

The measure is selected with ``score`` (CLI: ``--score``):

* ``auc_cc``     -- ``AUC_curvecurator``: area under the fitted curve, viability normalised against
  the vehicle control. A *higher* value means a *less* sensitive (more resistant) cell line, and 1.0
  is the no-effect level, because CTRP's fit pins the low-concentration asymptote to the DMSO value
  (Seashore-Ludlow et al., *Cancer Discov* 5(11):1210-1223, 2015, Methods).
* ``ln_ic50_cc`` -- ``LN_IC50_curvecurator``: natural log of the IC50 from the same fit. Here a
  *lower* value means *more* sensitive -- the opposite direction. Missing for ~40 % of curves by
  construction; CTRPv2 itself publishes no IC50.

**What this replaced, and why it is a replacement rather than a repair.** Until 11.08.2026 the target
was CTRP's published ``area_under_curve`` divided by ``conc_pts_fit``. That divisor counts the
concentration points which survived outlier censoring during CTRP's curve fit, not the width of the
integral, so the target was inflated for cell lines whose measurements were noisy enough to lose
points. Every number computed on it is void -- see
``docs/steps/corrections-and-dead-ends.md#the-auc-target-was-divided-by-the-wrong-quantity``.

Two outputs are written into the AnnData (both happen by default):

1. Multi-drug target matrix (used by the multi-task training):
    * ``adata.obsm["Y_ctrp"]``  : float32 (n_cells, K), NaN where missing.
    * ``adata.obsm["M_ctrp"]``  : bool    (n_cells, K), True where Y_ctrp is observed.
    * ``adata.uns["ctrp_drugs"]``: list[str] of length K (normalized drug names),
      giving the column order of Y_ctrp / M_ctrp.
    * ``adata.uns["ctrp_score"]``: which measure Y_ctrp holds.
    * ``adata.obs["cellosaurus_id"]``: the cell line's Cellosaurus accession -- **recorded, never
      joined on**; the join key is still the normalized name. See
      ``notebooks/data_and_harmonization/cell_line_join_verification.ipynb``.

   Drugs are kept only if at least ``min_cell_lines`` distinct SCP542-overlapping
   cell lines were screened against them (default 50) so we don't add heads with
   too little support.

2. Per-drug flat columns (back-compat with the original single-drug pipeline):
    * ``adata.obs["viability_<drug>"]``    : per-cell target, NaN when missing.
    * ``adata.obs["train_mask_<drug>"]``   : per-cell bool, True when present.

   The ``viability_`` prefix is historical -- the column holds whichever measure
   was selected, not necessarily a percent viability. Controlled by
   ``extra_single_drug_cols``. Defaults to ("paclitaxel",) so the 25.05.2026
   baseline remains reproducible without any code changes.

**Measures that were removed, not kept.** ``auc`` and the legacy ``mean_pv`` went with the CTRP-file
readers on 11.08.2026. ``auc_z`` (per-drug z-scored ``auc``) was the default from 13.07.2026 until
retired on 27.07.2026, and its code was removed on 11.08.2026: any normalization across drugs belongs
in the loss, where it can be estimated per fold, rather than baked into the stored target. The ``uns``
keys that existed only to invert it -- ``ctrp_score_center`` / ``ctrp_score_scale`` -- went with it.
h5ads written before that date still carry them, and nothing reads them.
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

# CTRPv2 and SCP542 share no persistent identifier, so the cell-line NAME is the only
# join key (see _normalize_cell_line). Where CTRP spells a line differently from CCLE --
# which is what SCP542's Cell_line column holds -- the line is silently dropped even
# though it was screened. This table maps a normalized CTRP name onto its normalized
# CCLE/SCP542 name. Every entry carries the evidence for the identification; do not add
# one without it.
#
# h292 -> ncih292  (accepted 10.08.2026, audit item 02)
#   Cellosaurus CVCL_0455 is "NCI-H292", listing both "H292" and "NCIH292" among its
#   synonyms, disease "lung mucoepidermoid carcinoma", CCLE member, DepMap ACH-001075.
#   CTRP's row (master_ccl_id 290) agrees on every field it carries:
#   ccl_availability=ccle;public, ccle_primary_site=lung,
#   ccle_hist_subtype_1=mucoepidermoid_carcinoma -- the only lung mucoepidermoid
#   carcinoma among CTRP's 1,107 lines. CTRP writes 106 other lines as `NCIH...`; this is
#   the single place it drops the prefix. Recovers 213 SCP542 cells and 454 drug labels,
#   taking the trainable overlap from 180 to 181 cell lines.
CTRP_CELL_LINE_ALIASES: dict[str, str] = {"h292": "ncih292"}


def _normalize_cell_line(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().str.replace("-", "")


def _normalize_drug(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower()


#: Which column of ``CTRPv2.csv`` each response measure reads. Both come from the same CurveCurator
#: fit, so selecting between them changes only what is summarised from that fit -- not which curves
#: were fitted, nor how.
SCORE_COLUMNS: dict[str, str] = {
    "auc_cc": "AUC_curvecurator",
    "ln_ic50_cc": "LN_IC50_curvecurator",
}


def _load_drevalpy_long(response_csv: Path, score: str) -> pd.DataFrame:
    """Read DrEval's reprocessed CTRPv2 into a long table with normalized name columns.

    One row per measurement as published. A (cell line, drug) pair can appear more than once --
    on record 21807175 always as an exact duplicate row -- so :func:`_deduplicate_measurements`
    reduces it afterwards.

    Columns returned: ``ccl_name_norm``, ``cpd_name_norm``, ``score``, ``r2``, ``cellosaurus_id``.

    Two columns are joined on rather than CTRP's own name fields:

    * ``ccl_name`` is CTRP's spelling, **not** DrEval's display ``cell_line_name``. The display name
      loses 19 SCP542 lines to punctuation (``MDA-MB-361`` vs ``MDAMB361``); CTRP's spelling is what
      ``_normalize_cell_line`` was built against and matches 180 of 198.
    * ``cellosaurus_id`` is carried through but **never joined on** -- resolving accessions is
      strictly worse as a key (172 of 198). It rides along so the target carries a persistent
      identifier, and so the join can be verified against one
      (``notebooks/data_and_harmonization/cell_line_join_verification.ipynb``).

    Rows with no value for the requested measure are dropped here rather than silently becoming
    unobserved entries. That matters for ``ln_ic50_cc``, which is absent for ~40 % of curves by
    construction -- see :data:`scripts.preprocessing.layout.CTRP_SCORES`.
    """
    column = SCORE_COLUMNS[score]
    df = pd.read_csv(
        response_csv,
        usecols=["ccl_name", "drug_name", "cellosaurus_id", "R2", column],
    )
    n_all = len(df)
    df = df.dropna(subset=[column])
    if len(df) < n_all:
        print(
            f"  {n_all - len(df):,} of {n_all:,} rows have no {column} "
            f"({100 * (n_all - len(df)) / n_all:.1f} %) and are dropped."
        )

    df = df.rename(columns={column: "score", "R2": "r2"})
    # Aliases are applied on the CTRP side only: they map CTRP's spelling onto the CCLE name SCP542
    # uses, so both sides end up in the same key space.
    df["ccl_name_norm"] = _normalize_cell_line(df["ccl_name"]).replace(CTRP_CELL_LINE_ALIASES)
    df["cpd_name_norm"] = _normalize_drug(df["drug_name"])
    print(
        f"  {len(df):,} measurements | {df.ccl_name_norm.nunique():,} cell lines | "
        f"{df.cpd_name_norm.nunique():,} drugs"
    )
    return df[["ccl_name_norm", "cpd_name_norm", "score", "r2", "cellosaurus_id"]]


def _deduplicate_measurements(long: pd.DataFrame) -> pd.DataFrame:
    """One row per (cell line, drug), by dropping **exact** duplicates -- and refusing anything else.

    ``CTRPv2.csv`` repeats rows: on record 21807175, 15,946 of its 395,024 rows are duplicates of
    another row across **all 46 columns** -- same compound identifiers, same fitted values, same
    ``R2`` -- forming 7,331 pairs and 428 triples. They are row duplication in the published table,
    not repeat experiments of the same pair, so removing them changes no value.

    **Decided 11.08.2026 (Selin): de-duplicate exact rows, and raise on anything else.** A rule that
    silently reconciles disagreeing rows -- averaging them, or keeping the better ``R2`` -- would keep
    working if a future record started carrying genuinely different repeats, and the target would
    change with nothing to show for it. Failing loudly turns that into an error someone has to decide
    about, which is the only thing that makes the decision revisitable.

    :raises ValueError: if a (cell line, drug) pair survives with rows that are not identical.
    """
    keys = ["ccl_name_norm", "cpd_name_norm"]
    n_before = len(long)
    out = long.drop_duplicates().sort_values(keys, kind="mergesort").reset_index(drop=True)
    if n_before > len(out):
        print(
            f"  {n_before - len(out):,} of {n_before:,} rows are exact duplicates of another row "
            f"({100 * (n_before - len(out)) / n_before:.1f} %) and are dropped."
        )

    conflicting = out[out.duplicated(keys, keep=False)]
    if not conflicting.empty:
        sample = conflicting.sort_values(keys).head(6).to_string(index=False)
        raise ValueError(
            f"{conflicting.drop_duplicates(keys).shape[0]:,} (cell line, drug) pairs have more than "
            f"one non-identical measurement. The pipeline has no rule for reconciling them, "
            f"deliberately -- decide how to aggregate before proceeding.\n{sample}"
        )
    return out


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
    # Already one row per (cell line, drug) -- :func:`_deduplicate_measurements` did that, and did it on
    # the full table so the choice of surviving row cannot depend on which cell lines overlap.
    long_overlap = ctrp_full[ctrp_full["ccl_name_norm"].isin(overlap_cell_lines_norm)]

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
    response_csv: str,
    min_cell_lines: int = DEFAULT_MIN_CELL_LINES,
    target_drugs: Sequence[str] | None = None,
    extra_single_drug_cols: Sequence[str] = DEFAULT_EXTRA_SINGLE_DRUG_COLS,
    score: str = DEFAULT_CTRP_SCORE,
):
    """Map CTRPv2 drug-response measures onto cells in the embedded AnnData.

    Parameters
    ----------
    response_csv:
        DrEval's reprocessed CTRPv2 table (``CTRPv2.csv``), fetched and MD5-verified by
        ``fetch_ctrp_response`` from a pinned Zenodo record.
    score:
        Which response measure to use as the target: ``auc_cc`` (default) or ``ln_ic50_cc``.
        See :data:`scripts.preprocessing.layout.CTRP_SCORES`.
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

    print(f"Loading CTRPv2 responses from {Path(response_csv).name} (score={score})...")
    ctrp_full = _deduplicate_measurements(_load_drevalpy_long(Path(response_csv), score))

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

    # Persistent identifier per cell line, taken from the response table. **Recorded, never joined
    # on** -- accessions resolve fewer SCP542 lines than names do (172 vs 180 of 198), so they make a
    # worse key; they are here so the target carries a citable identifier and so the name join can be
    # checked against an external authority. Verification and the ambiguity rules:
    # notebooks/data_and_harmonization/cell_line_join_verification.ipynb.
    line_to_cvcl = ctrp_full.drop_duplicates("ccl_name_norm").set_index("ccl_name_norm")[
        "cellosaurus_id"
    ]
    adata.obs["cellosaurus_id"] = pd.array(
        cell_line_norm.map(line_to_cvcl).to_numpy(), dtype="string"
    )
    n_lines_with_id = int(adata.obs.loc[:, ["Cell_line", "cellosaurus_id"]]
                          .drop_duplicates()["cellosaurus_id"].notna().sum())
    print(f"  Cellosaurus accessions attached for {n_lines_with_id} cell lines "
          f"(cells: {int(adata.obs['cellosaurus_id'].notna().sum()):,} / {adata.n_obs:,}).")

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
        "--response-csv",
        type=Path,
        default=None,
        help="DrEval CTRPv2 response table (default: the pinned cache under <data-root>/metadata/).",
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
        response_csv=str(args.response_csv or paths.ctrp_response_csv),
        min_cell_lines=min_cell_lines,
        target_drugs=args.drugs,
        extra_single_drug_cols=tuple(args.single_drug_cols),
        score=args.score,
    )
