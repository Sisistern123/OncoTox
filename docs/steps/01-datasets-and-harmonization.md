# Step 01 — Datasets & harmonization

*Part of [OncoTox project progress](../project_progress.md). Covers: the raw datasets
collected, the cross-dataset overlap/coverage audit, the drug catalog, and the harmonization
strategy that aligns cell lines and compounds across sources.*

Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition**.

---

## Data collection (26.03–30.03.2026)

| Dataset | Role | Key numbers | Used? |
|---|---|---|---|
| **SCP542** scRNA-seq (PERCEPTION paper, Kinker et al. 2020) | single-cell input | **53,513 cells × 22,722 genes**, **198 unique cell lines** | ✅ primary |
| **CTRPv2** | viability labels | 1,107 cell lines, **545 compounds**, target `cpd_avg_pv` | ✅ primary |
| **PRISM** Repurposing (Public 24Q2) | single-dose LFC | 915 cell lines, 6,575 compounds | downloaded, not used |
| **GDSC2** | IC50 / AUC | 967 cell lines, 295 drugs | downloaded, not used |

- SCP542 source: SCP542 / Broad Single Cell Portal; `.X` stored as **CPM**.
- CTRPv2 raw tables used: `v20.data.per_cpd_post_qc.txt` (`cpd_avg_pv`),
  `v20.meta.per_experiment.txt`, `v20.meta.per_cell_line.txt`, `v20.meta.per_compound.txt`.

✅ On-plan: SCP542 + CTRPv2 are the designated primary pair; PRISM/GDSC reserved for later.

---

## Overlap & coverage audit (03.04.2026) — `notebooks/compare_GDSC_CTRP.ipynb`

The work behind the plan's **Fig. 1 / Fig. 2**.

**Cell-line overlap with SCP542** (normalized: trim + lowercase; SCP542/PRISM split on `_`):

| Dataset | Total lines | Overlap w/ SCP542 | Missing |
|---|---|---|---|
| GDSC | 967 | 133 | 65 |
| CTRPv2 | 1,107 | **190** | 8 |
| PRISM | 915 | 182 | 16 |

**Drug / compound harmonization** — unified catalog `data/drug/all_sources_drug_catalog.csv`
(7,415 source rows; GDSC 295 / CTRPv2 545 / PRISM 6,575; union 7,040):

- Name overlap (normalized): CTRPv2↔GDSC 66, CTRPv2↔PRISM 218, GDSC↔PRISM 144.
- **BRD-ID overlap CTRPv2↔PRISM: 243** (higher-confidence link).
- Exports: `data/drug/drug_overlap_candidates.csv`.

**DrugBank match** (from `full database.xml`, normalized names + synonyms):
GDSC 118/295 (40.0 %), CTRPv2 173/545 (31.7 %), PRISM 3,483/6,575 (53.0 %),
overall 3,774/7,415 (50.9 %). Exports: `drugbank_overlap_matches.csv`, `..._unmatched.csv`.

**Applicable (non-null) response coverage within SCP542 overlap:**

| Dataset | Metric | Non-null / total | % | % within overlap subset |
|---|---|---|---|---|
| GDSC | `LN_IC50` | 8,007 / 242,036 | 3.31 % | 100 % |
| CTRPv2 | `cpd_avg_pv` | 1,521,028 / 7,227,951 | 21.04 % | **100 %** |
| PRISM | extended primary matrix | 1,210,432 / 4,213,048 | 28.73 % | 97.95 % |

✅ On-plan: satisfies sub-goal 1 (harmonization incl. BRD + DrugBank) and supplies the
Fig. 1/2 numbers the plan rests sub-goal 3 on.

> ⚠️ **Number to reconcile:** this audit reports **190** overlapping cell lines
> (case-insensitive). The training pipeline (`ctrp_to_h5ad.py`) normalizes more strictly
> (also strips `-`) and reports **180** at run time (see
> [Step 05](05-multitask-results.md)). Same data, different normalization — pick one for the
> thesis figure and the pipeline.

---

## Harmonization strategy (reference)

How the sources are aligned — gathered for the writeup:

- **Cell lines:** name normalization (trim + lowercase; SCP542/PRISM split on `_`; the pipeline
  additionally strips `-`) → SCP542∩CTRPv2 = **190** by the audit normalization, **180** by the
  stricter pipeline one (the number to reconcile above).
- **Drugs:** matched three ways — normalized **name**, **BRD-ID** (CTRPv2↔PRISM 243, the
  higher-confidence link), and **DrugBank** name+synonym match (CTRPv2 173/545) — see the audit
  above for all pairwise counts and exports (`drug_overlap_candidates.csv`,
  `drugbank_overlap_matches.csv`).
- **Expression:** SCP542 `.X` kept as **CPM**; HVG selection done on a `log1p` copy but the saved
  matrix stays CPM (see [Step 02](02-preprocessing-and-embeddings.md)). Two harmonized variants on
  disk: `hvg5000` (default) and `all_genes`.
- **PRISM / GDSC are harmonized but not yet wired into training** (cross-database heads = the open
  Phase-3b — see the scorecard in the [index](../project_progress.md)).

## Drugs — scope and catalog (reference)

- **CTRPv2 compound set: K = 545** drugs (the `--all-drugs` run). Optional coverage filter keeps
  a drug only if screened on ≥ `--min-cell-lines` overlapping lines (default 50) — the headline
  run used min 0 = all 545 (see [Step 05](05-multitask-results.md)).
- Each drug = one **column/head** of the shared output layer; identity is *not* an input feature.
- Single-task work uses **paclitaxel** as the reference drug (see
  [Step 04](04-single-task-results.md)).
- Unified cross-source catalog `data/drug/all_sources_drug_catalog.csv` (GDSC 295 / CTRPv2 545 /
  PRISM 6,575; union 7,040) — built for the eventual cross-database heads.

---

## Code, notebooks & key variables

- **Notebook — the whole audit:** `notebooks/compare_GDSC_CTRP.ipynb` — cell-line overlap,
  name/BRD/DrugBank drug matching, the Fig. 1 / Fig. 2 coverage numbers above.
- **Catalog + export artifacts** (`data/drug/`): `all_sources_drug_catalog.csv` (unified catalog),
  `drug_overlap_candidates.csv`, `drugbank_overlap_matches.csv`, `drugbank_overlap_unmatched.csv`.
- **Raw CTRPv2 source tables** (`metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/`):
  `v20.data.per_cpd_post_qc.txt` (`cpd_avg_pv`), `v20.meta.per_experiment.txt`,
  `v20.meta.per_cell_line.txt` (`ccl_name`), `v20.meta.per_compound.txt` (`cpd_name`).
- **Normalization used by the pipeline** (the stricter 180-overlap rule): `_normalize_cell_line`
  and `_normalize_drug` in `scripts/preprocessing/ctrp_to_h5ad.py` (trim + lowercase + strip `-`),
  producing the `ccl_name_norm` / `cpd_name_norm` keys everything downstream joins on.
