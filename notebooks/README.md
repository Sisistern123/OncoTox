# OncoTox notebooks

Numbered in **pipeline order**. All figures/tables are written to `notebooks/outputs/`; full
per-run model artifacts live under `runs/` (git-ignored) and are indexed in `runs/runs_index.csv`.

## TL;DR — what you actually need

Only **two** notebooks are on the critical path for the pipeline:

1. **`05_preprocessing.ipynb`** — builds / refreshes the trainable data.
2. **`07_training.ipynb`** — produces every model result (matrix, CV, correlation, HVG sweep).

Everything else (`01`, `02`, `03`, `04`, `06`) is **exploration / harmonization / QC**: it informed
design decisions and helps interpret the results, but is **not required** to (re)produce the headline
numbers. Run order for a clean reproduction is simply **05 → 07**.

| # | Notebook | Role | Essential for Pipeline? |
|---|---|---|---|
| 01 | `01_scDAExploration.ipynb` | Initial single-cell (SCP542) data exploration | No — exploration |
| 02 | `02_compare_GDSC_CTRP.ipynb` | Drug-catalog harmonization (CTRP/GDSC/DrugBank); writes `data/drug/*` catalogs | No — one-off harmonization |
| 03 | `03_analysis.ipynb` | CTRP→PRISM drug-repurposing / clinical-phase mapping | No — metadata enrichment |
| 04 | `04_drug_coverage.ipynb` | Per-drug coverage & response variance ("which drugs are learnable") | No — informs interpretation/thresholds |
| **05** | **`05_preprocessing.ipynb`** | **Build the trainable h5ad (incl. 512-d PCA, HVG variants)** | **Yes — data** |
| 06 | `06_verify_variants.ipynb` | QC audit of preprocessing outputs + PCA-vs-scGPT UMAPs | No — validation/QC |
| **07** | **`07_training.ipynb`** | **All training + evaluation** | **Yes — results** |

The scripts these notebooks call (`scripts/preprocessing/run_preprocessing.py`,
`scripts/training/train_multitask.py`) do **not** read any output of `02/03/04` — those are analysis
side-products, not pipeline inputs.

---

## The two essential notebooks

### `05_preprocessing.ipynb` — data
A documented front-end to `scripts/preprocessing/run_preprocessing.py` (the 5-step pipeline:
`convert → scgpt → targets → splits → pca`). It does not reimplement anything — it calls the script
so the notebook and CLI can't drift.

- **§A — recompute the 512-d PCA baseline** for the two built variants (`hvg5000`, `all_genes`).
  Idempotent (`--start-at pca --force-pca --pca-n-comps 512`); this is the step the 512-d switch needed.
- **§B — HVG-count sweep data-gen**: build `hvg1000/2000/3000` (full pipeline incl. **scGPT
  re-embedding**, hours + GBs). Gated behind `RUN_HVG_SWEEP` so it doesn't run by accident.

Output per variant: the trainable `…_with_targets.h5ad` carrying `X_scGPT`, `X_pca`, `Y_ctrp`,
`M_ctrp`, `split_ctrp`. (The core `hvg5000`/`all_genes` data was first built from the CLI; this
notebook documents and refreshes it.)

### `07_training.ipynb` — results
Everything model-side, all on the matched setup (same `(128,64)` trunk, same 512-d input; only the
representation differs). Every fit calls `train_multitask.train_rep` / `cv_evaluate` (the same code the
CLI uses).

- **§1 — the 8-run matrix** `{hvg5000, all_genes} × {X_pca, X_scGPT} × {all-drugs, single-paclitaxel}`,
  **load-or-train**: re-running loads the saved `runs/` instead of retraining (`RETRAIN_MATRIX=True` to
  force). Produces the all-drugs / single-task tables + the per-drug scatter.
- **§2 — cross-validation**: 5-fold GroupKFold over `Cell_line`, **test held out** (resamples the 153
  train+val lines) → mean ± std for heads-beating, **Δmse** (continuous model−baseline), val MSE, and
  the overfitting gap. The train/val/test split (70/15/15) is documented here.
- **§3 — per-drug correlation**: Spearman/Pearson of predicted vs true viability across held-out cell
  lines, for drugs with real response variance.
- **§4 — HVG sweet spot**: heads-beating vs HVG count (1k/2k/3k/5k) under CV, all drugs.

**Caching flags** (default `False` = load saved results, fast re-run): `RETRAIN_MATRIX` (§1),
`RECOMPUTE_CV` (§2, uses `outputs/cv_folds.csv`), `RECOMPUTE_SWEEP` (§4, uses `outputs/hvg_sweep.csv`).
Set a flag `True` to recompute that section. Metric definitions are in
[`docs/steps/05`](../docs/steps/05-multitask-results.md#metrics--what-each-number-means).

---

## Supporting notebooks (understanding / data exploration)

- **`01_scDAExploration.ipynb`** — first look at the SCP542 single-cell data.
- **`02_compare_GDSC_CTRP.ipynb`** — cross-database drug-name/BRD/DrugBank harmonization; writes the
  catalogs under `data/drug/`. A one-off audit ([Step 01](../docs/steps/01-datasets-and-harmonization.md)).
- **`03_analysis.ipynb`** — maps CTRP compounds to PRISM repurposing metadata + clinical phase.
- **`04_drug_coverage.ipynb`** — per-drug coverage and response-distribution analysis; this is where
  "which drugs are even learnable" (and the real-variance threshold used in `07` §3) comes from.
- **`06_verify_variants.ipynb`** — re-runnable QC of the preprocessing outputs (gene counts, `X_pca`
  source, cell alignment) and the PCA-vs-scGPT UMAPs (Fig. 3/4).
