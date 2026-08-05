# OncoTox notebooks

**A number means pipeline.** `1_` → `2_` → `3_` is the run order that rebuilds the data and the results.
Everything else is analysis and lives in a named directory — no numbers, because there is no order to run
them in. Figures and tables go to [`outputs/`](outputs/); model artifacts to `runs/` (git-ignored),
indexed by `runs/runs_index.csv`.

> ⛔ **The drug panel is void and a pipeline review is in progress** ([`docs/TODO.md`](../docs/TODO.md)).
> Nothing computed on that panel is quotable until it is rebuilt — which affects `3_panel_training` and
> everything under `drug_selection/`.

## The pipeline

| # | Notebook | What it does | Drives |
|---|---|---|---|
| **1** | `1_preprocessing.ipynb` | Builds the trainable data: recomputes the 512-d PCA baseline for both variants (§A) and builds the HVG-sweep variants including the scGPT re-embed (§B) | `scripts/preprocessing/run_preprocessing.py` |
| **2** | `2_training.ipynb` | The PCA-vs-scGPT harness: 8-run matrix (load-or-train), 5-fold GroupKFold CV with test held out, per-drug correlation | `train_multitask.train_rep`, `cv_evaluate` |
| **3** | `3_panel_training.ipynb` | The current training run: raw `auc`, per-fold density weighting, out-of-fold scoring against the ridge control | `scripts/training/cv.py`, `density_weighting.py` |

`1_` and `2_` call the **same entry points the CLI uses**, so the notebooks and the command line cannot
drift — they are documentation *and* a re-run, not a fork.

> ⚠️ **`2_training`'s conclusions are superseded; its machinery is not.** The 8-run matrix and the
> "ρ ≈ 0, the model cannot rank cell lines" reading were produced at K=545 on the legacy `mean_pv`
> target, whose unstandardised per-drug variance was destroying the signal
> ([why](../docs/steps/corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result)).
> It stays in the pipeline because re-running the matrix on the current target is an open task and this is
> the only place that harness exists.

## Analysis

### `data_and_harmonization/` — what is in the data, and does the join hold?

| Notebook | Question it answers |
|---|---|
| `drug_catalog.ipynb` | Cross-database compound harmonization (CTRP / GDSC / PRISM / DrugBank) → writes `data/drug/*`. The **only** analysis notebook whose outputs another step consumes |
| `drug_coverage.ipynb` | Per-drug coverage and response spread; the label distribution behind "why the task is hard". ⚠️ Its *learnability* section was built on `mean_pv` and is superseded; the target-distribution figures still stand |
| `verify_variants.ipynb` | QC of `hvg5000` vs `all_genes`, the PCA-vs-scGPT UMAP latent validation, and (§9) the gene-set sweep — heads-beating vs gene count under CV, moved here from `2_training` §4 on 03.08.2026 and re-targeted to `auc` |

### `result_evaluation/` — is the number real?

| Notebook | Question it answers |
|---|---|
| `target_comparison.ipynb` | **Which target?** `mean_pv` vs `auc` vs `auc_z` at K=10 and K=545, with bootstrap CIs, Pearson alongside Spearman, and seed stability |
| `ablations_and_rescue.ipynb` | **Why did the 545-head model fail, and what fixes it?** The implicit σ²-weighting of the loss, the causal rescue test on the broken setting, the model-knob ablations on the corrected one, and the ridge control |
| `dreval_benchmark.ipynb` | **How strong is this by the field's standard?** Our data and model through the real **DrEval** package (`drevalpy` 1.5.1): their LCO splits, their baselines, their metrics |
| `diagnostics.ipynb` | The drug-selection gate defect, the proliferation test, and result dispersion across folds *and* drugs |

### `drug_selection/` — which compounds enter the model

The project's main open deliverable. All three currently rest on a criterion that
[measured the wrong quantity](../docs/steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability),
and all three get re-run after the rebuild.

| Notebook | What it was |
|---|---|
| `learnability_filter.ipynb` | The kill/spare gate that took 545 drugs to 10 — **discredited criterion** |
| `learnable_subset_training.ipynb` | PCA vs scGPT on that subset — a best-case diagnostic, never a generalization number |
| `panel_distributions.ipynb` | Response distributions and the weighting design on the 8-drug panel — **void panel** |

### `archive/` — nothing here is load-bearing

**The dividing line: everything outside `archive/` is cited by a step file and something documented
depends on it. Nothing in `archive/` is.** That is the test, applied by checking which notebooks the step
files actually reference — not by judging how interesting they look.

| Notebook | Why it is here |
|---|---|
| `scdrugatlas_exploration.ipynb` | Explores **scDrugAtlas**, a data source that was evaluated and [rejected](../docs/steps/corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources). Kept as the record of that decision. *(Long mislabelled in the docs as SCP542 exploration — "scDA" is scDrugAtlas.)* |
| `ctrp_prism_repurposing.ipynb` | CTRP→PRISM repurposing and clinical-phase mapping. Read-only, writes no artifact, and no step depends on it. Worth knowing it exists: it is the only notebook that loads `GDSC2_fitted_dose_response_27Oct23.xlsx`, which the "externalize the spread requirement" task will need |

Two things are **not** grounds for archiving. **Superseded conclusions** — `2_training`'s results are void
but its harness is the only way to re-run the matrix, so it stays in the pipeline. **A discredited
criterion** — everything in `drug_selection/` is built on one, but it is the documented record of how
selection was done and it re-runs after the rebuild.

## Re-running

Every training notebook has a **`RETRAIN` flag, default `False`** — it then loads the saved CSVs from
`outputs/` and only redraws the figures, in seconds. Set it `True` to refit.

Fits use `TrainConfig(epochs=25)`: across 36 recorded runs the best epoch was **median 6, max 11**, and
early stopping (patience 10) never came close to 25 — the earlier cap of 50 only cost wall-clock.

Data is built by the CLI, not by the notebooks (`--score` defaults to `auc`):

```bash
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --score auc --start-at targets --skip-scgpt      # one targets h5ad per --score
```

All results in this repo use the **`hvg5000`** variant: scGPT embeds only the subset of those genes in
its vocabulary, PCA is computed on all of them, and both representations come out **512-d**. Gene counts
per stage: [Step 02](../docs/steps/02-preprocessing-and-embeddings.md#hvg-5000-pipeline-outputs).
⛔ The OOV drop is **not** clean — most of it is a symbol-matching defect, not scGPT's vocabulary:
[Corrections](../docs/steps/corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols).

The directories under [`outputs/`](outputs/) are deliberately **not** named after notebooks, because they
do not map one-to-one: `data/` is written by `drug_coverage` *and* `ablations_and_rescue`, `target/` by
`target_comparison` *and* `ablations_and_rescue`. They are named for what they contain.

## Where the written record lives

`docs/steps/01`–`05` hold the numbers and own every claim ·
[`corrections-and-dead-ends.md`](../docs/steps/corrections-and-dead-ends.md) holds everything superseded,
retracted, refuted or abandoned · [`docs/TODO.md`](../docs/TODO.md) is the action list.
