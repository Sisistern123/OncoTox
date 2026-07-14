# OncoTox notebooks

Numbered in **pipeline order**. All figures/tables are written to `notebooks/outputs/`; full
per-run model artifacts live under `runs/` (git-ignored) and are indexed in `runs/runs_index.csv`.

## TL;DR — what you actually need

> ⚠️ **`07`'s conclusions are superseded (13.07.2026).** They were produced at K=545 on the legacy
> `mean_pv` target, whose **unstandardized per-drug variance was destroying the signal** — `11`
> reproduces that failure on demand. Read **`08 → 09 → 10 → 11`** for the current picture, and treat
> `07`'s 8-run matrix / "ρ ≈ 0" results as pending a re-run on `--score auc_z`.

Critical path for the pipeline:

1. **`05_preprocessing.ipynb`** — builds / refreshes the trainable data (per `--score`).
2. **`07_training.ipynb`** — the 8-run matrix, CV, correlation, HVG sweep (**on `mean_pv` — stale**).
3. **`08` → `09`** — learnability filter → the current headline result (PCA vs scGPT on the learnable
   drugs).
4. **`10`, `11`** — why it works: what is *not* the bottleneck (model) and what *is* (the target).

`01`, `02`, `03`, `04`, `06` are **exploration / harmonization / QC** — not required to reproduce the
numbers.

| # | Notebook | Role | Essential for Pipeline? |
|---|---|---|---|
| 01 | `01_scDAExploration.ipynb` | Initial single-cell (SCP542) data exploration | No — exploration |
| 02 | `02_compare_GDSC_CTRP.ipynb` | Drug-catalog harmonization (CTRP/GDSC/DrugBank); writes `data/drug/*` catalogs | No — one-off harmonization |
| 03 | `03_analysis.ipynb` | CTRP→PRISM drug-repurposing / clinical-phase mapping | No — metadata enrichment |
| 04 | `04_drug_coverage.ipynb` | Per-drug coverage & response variance ("which drugs are learnable") | No — superseded by `08` on the AUC target |
| **05** | **`05_preprocessing.ipynb`** | **Build the trainable h5ad (incl. 512-d PCA, HVG variants, `--score`)** | **Yes — data** |
| 06 | `06_verify_variants.ipynb` | QC audit of preprocessing outputs + PCA-vs-scGPT UMAPs | No — validation/QC |
| 07 | `07_training.ipynb` | 8-run matrix, CV, per-drug correlation, HVG sweep | ⚠️ **Stale** (`mean_pv`, K=545) |
| **08** | **`08_learnability_filter.ipynb`** | **Harsh learnability filter on the AUC target → 5 drugs** | **Yes — drug scope** |
| **09** | **`09_learnable5_training.ipynb`** | **PCA vs scGPT on the 5 drugs — the current headline** | **Yes — results** |
| 10 | `10_ablations.ipynb` | Is it the model? (regularization / capacity / batch / weighting + **ridge control**) | Yes — closes model tuning |
| 11 | `11_auc_vs_aucz.ipynb` | Does the per-drug z-scoring help? (**yes, decisively, at K=545**) | Yes — justifies the target |

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

Output per variant **and target score** (`--score {auc_z,auc,mean_pv}`, default `auc_z`): the
trainable `…_with_targets[_<score>].h5ad` carrying `X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`,
`split_ctrp`. (The core `hvg5000`/`all_genes` data was first built from the CLI; this notebook
documents and refreshes it.)

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
`RECOMPUTE_CV` (§2, uses `outputs/03_training_545/cv_folds.csv`), `RECOMPUTE_SWEEP` (§4, uses `outputs/03_training_545/hvg_sweep.csv`).
Set a flag `True` to recompute that section. Metric definitions are in
[`docs/steps/05`](../docs/steps/05-multitask-results.md#metrics--what-each-number-means).

### `08_learnability_filter.ipynb` → `09_learnable5_training.ipynb` — the learnable-drug pair (13.07.2026)

Run in order; `09` reads `08`'s CSV, so re-running `08` with different gates changes what `09` trains.

- **`08`** — harsh learnability filter on the `auc_z` target. The `04` score (`resp_std × cov_frac`) is
  **degenerate** here (z-scoring makes every drug's std 1.0), so spread is read off the raw `auc` scale
  via `uns["ctrp_score_scale"]`. Adds the condition the loose filter lacked — a drug must **kill** a real
  population of lines *and* **leave one alive**. **5 / 545 survive** → `outputs/04_learnability/ctrp_drug_learnability_auc.csv`,
  `outputs/04_learnability/learnability_filter_auc.png`.
- **`09`** — PCA vs scGPT trained on those 5 (via `train_rep`, matched trunk). Headline metric is per-drug
  Spearman on **cross-validated out-of-fold predictions** (5-fold GroupKFold over the 153 train+val lines,
  ~150 lines/drug — the fixed val split only has 27). **Mean Spearman 0.43 (PCA) / 0.49 (scGPT)** vs ≈ 0
  over all 545 drugs → the 545-drug null result was a **drug-selection artifact**
  ([Step 05](../docs/steps/05-multitask-results.md)). Outputs: `04_learnability/learnable5_per_drug_correlation.csv`,
  `04_learnability/learnable5_fixed_split_mse.csv`, `04_learnability/learnable5_pred_vs_true.png`, `04_learnability/learnable5_pca_vs_scgpt.png`.

⚠️ Both are a **best-case diagnostic**: the 5 drugs are selected using all 180 lines (val/test included).
Fine for "does any signal exist?", not a generalization estimate — see [TODO](../docs/TODO.md).

### `10_ablations.ipynb` — is it the model? (no)

Four model-side knobs on the 5 learnable drugs, same out-of-fold Spearman metric as `09`: **regularization**
(none→heavy), **capacity** (74,629→2,565 params, down to a linear head), **batch size** (32/128/512) and
**sample reweighting** (line-balanced, focus-extremes). **All flat** — the defaults are at/near the best on
every axis, so **model-side tuning is closed**. Includes the control that now sets the bar: **`RidgeCV` on
the 150 cell-line mean embeddings** (no cells, no network) scores **0.428**, *tying* the PCA MLP and within
0.06 of scGPT's — because the label is per cell line, so there are only ~150 independent examples.
Outputs: `06_ablations/ablation_regularization.csv`, `06_ablations/ablation_capacity.csv`, `06_ablations/ablation_batch_weighting.csv`,
`06_ablations/ablation_reg_capacity.png`. See [Step 03](../docs/steps/03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026).

### `11_auc_vs_aucz.ipynb` — which target? (`mean_pv` vs `auc` vs `auc_z`)

All three CTRPv2 scores, identical model, K=5 and K=545, scored **out-of-fold** on one common yardstick
(the curve-fit AUC ranking — `mean_pv` is *not* an affine map of `auc`, so scoring each against its own
target would make the columns answer different questions). Reports **Spearman + Pearson**, a **95%
bootstrap CI over the ~150 held-out lines**, per-drug dots, and a **3-seed stability check**.

| | `mean_pv` | `auc` | `auc_z` |
|---|---|---|---|
| K=5 (PCA / scGPT) | 0.450 / 0.481 | 0.439 / 0.482 | 0.424 / **0.488** |
| **K=545** (PCA / scGPT) | **+0.027 / −0.070** | **+0.016 / −0.087** | **+0.378 / +0.430** |

- **At K=5 all three tie** (CIs overlap) — nothing to equalize, so on a spread-homogeneous subset
  `--score auc` is equally good *and* keeps native AUC units.
- **At K=545 both unstandardized targets collapse** while `auc_z` holds → per-drug standardization is what
  makes a 545-head masked MSE trainable, and this **reproduces the `07` §3 null result on demand**.
- ⚠️ **The curve fit buys no accuracy** (`mean_pv` ≈ `auc` everywhere) — keep it for GDSC comparability,
  not performance.
- **scGPT − PCA = +0.075 ± 0.038**, sign-consistent over 3 seeds (`05_target/seed_stability.csv`).

Outputs: `05_target/target_comparison.csv`, `05_target/target_comparison_ci.csv`, `05_target/seed_stability.csv`, `05_target/target_comparison.png`.

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
