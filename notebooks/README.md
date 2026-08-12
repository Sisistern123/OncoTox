# OncoTox notebooks

**A number means pipeline.** `1_` → `5_` is one path from raw data to a scored result, and running them
in order rebuilds everything. Everything else is analysis and lives under [`analysis/`](analysis/) — no
numbers, because there is no order to run those in. Figures and tables go to [`outputs/`](outputs/);
model artifacts to `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

> ⛔ **The drug panel is void and a pipeline review is in progress** ([`docs/TODO.md`](../docs/TODO.md)).
> Nothing computed on that panel is quotable until a run exists on the rebuilt one — which affects
> `4a_percell_training` and everything derived from it.

## The pipeline

Renumbered 12.08.2026 (Selin) from three stages to five, so that the numbered chain is a complete
end-to-end pipeline rather than a partial one with the first build done by hand.

| # | Notebook | What it does | Drives |
|---|---|---|---|
| **1** | `1_data.ipynb` | Fetches the pinned CTRPv2 response table and converts SCP542 into the raw h5ad, HVG-filtered | `pipeline.fetch`, `pipeline.convert` |
| **2** | `2_drug_selection.ipynb` | **Which drugs does the model predict, and on whose authority?** FDA-approved list → what CTRPv2 screened → coverage over the overlap → the 11 with a verified published claim. Writes `outputs/panel/panel.csv` | `scripts/annotation/drug_annotation.py`, `scripts/sources/pubchem.py` |
| **3** | `3_representations.ipynb` | scGPT embedding, CTRP targets, the frozen cell-line-grouped split, and the PCA baseline | `pipeline.scgpt/targets/splits/pca` |
| **4** | `4a_percell_training.ipynb` | The training run: `auc_cc`, per-fold density weighting, out-of-fold scoring against the ridge control. ⚠️ Its stored outputs are cleared and its panel is the void 8-drug one; it re-runs on `outputs/panel/panel.csv` at R4 of the sweep | `scripts/training/cv.py`, `density_weighting.py` |
| **5** | `5_evaluation.ipynb` | ⛔ **Stub.** External benchmark + diagnostics. Neither notebook it should absorb can run yet — both hardcode the removed `'auc'`, and `dreval_benchmark` also imports a deleted module. Blocked on **review item 11** | — |

**Why drug selection is stage 2 and not an analysis notebook.** It reads exactly two things — the
cell-line roster stage 1 writes and the response table stage 1 fetches — and it writes `panel.csv`,
which stage 4 consumes. It is a step in the chain, and it sits before the representations because it
needs none of them.

**The panel does not feed stage 3.** `ctrp_to_h5ad` keeps every drug with enough coverage and the panel
is applied as a column selection at training time. Split eligibility is derived from *"this line has at
least one observed label"*, so restricting the target matrix to 11 drugs could drop lines and silently
re-freeze the split (Selin, 12.08.2026 — see `3_representations` §B).

**The notebooks are the orchestrator.** Each stage calls `scripts/preprocessing/pipeline.py`, where every
step owns its own guard and preconditions, so running them out of order fails rather than producing
something subtly wrong. `run_preprocessing.py`, the CLI that used to hold the step order, was
[archived](../scripts/archive/README.md) that day: the order is the numbering of these notebooks, and a
second copy of it in a CLI was a second thing to keep in step.

> ⚠️ **`4a_percell_training` has two sections, and they answer different questions.** §A is the panel run. §B is
> the PCA-vs-scGPT 8-run matrix, folded in from `2_training.ipynb` on 12.08.2026, which is retired by
> that merge. §B's conclusions are superseded — the matrix and the "ρ ≈ 0, the model cannot rank cell
> lines" reading were produced at K=545 on the legacy `mean_pv` target, whose unstandardised per-drug
> variance was destroying the signal
> ([why](../docs/steps/corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result))
> — but its harness is the only place that matrix exists, so it survives the notebook that held it.
>
> §B keeps `cv_evaluate` rather than moving to `oof_predictions`: the two became consistent on
> 12.08.2026 (both initialise the head bias and share `cv.grouped_folds`), and only their *return*
> differs — §A needs predictions, §B needs per-fold metrics such as `gap`, which `oof_predictions` does
> not return. Outputs go to `outputs/archive/training_545_mean_pv/`, not §A's `outputs/panel/`.
>
> §B's outputs staying under `archive/` (named `legacy/` when this was decided) is accepted, not a loose end (Selin, 12.08.2026).
>
> *(Corrected 12.08.2026: this read `outputs/matrix/`, a directory that has never existed. `OUT_MATRIX`
> named `outputs`/`matrix` while the call sites supplied `legacy`/`training_545_mean_pv` themselves, so
> anything reading it resolved into a directory that was not there — which is why the CV guard fell
> through and recomputed while reporting that it had loaded the committed folds. Fixed in `f6cbef4`,
> then consolidated so the variable means the directory it is named for; this sentence was the
> documentation that still described the broken value.)*

## Analysis

Nothing here is on the pipeline path; each answers a question about it.

### `analysis/qc/` — is the input what we think it is?

| Notebook | Question it answers |
|---|---|
| `verify_variants.ipynb` | QC of `hvg5000` vs `all_genes`, the PCA-vs-scGPT UMAP latent validation, and (§9) the gene-set sweep — heads-beating vs gene count under CV |
| `hvg_sweep_build.ipynb` | Builds the `hvg1000/2000/3000` variants §9 compares, through the same `pipeline.*` steps the numbered stages use. Was `1_preprocessing` §B; moved out of the pipeline 12.08.2026 because it exists for one analysis. ⚠️ Hours and gigabytes — gated behind `RUN_HVG_SWEEP` |
| `gene_symbol_rescue.ipynb` | How many genes scGPT discarded that are in its vocabulary under a **current** HGNC symbol — the symbol-matching defect, quantified |

### `analysis/harmonization/` — does the join hold?

| Notebook | Question it answers |
|---|---|
| `cell_line_join_verification.ipynb` | **Does the name join pair the right two lines?** Resolves SCP542's names against the pinned Cellosaurus release and checks them against the accessions DrEval ship, with an independent tissue cross-check. Runs `scripts/sources/cellosaurus.py`, so the rules it documents are the rules the pipeline uses |
| `drug_coverage.ipynb` | Per-drug coverage and response spread; the label distribution behind "why the task is hard". ⚠️ Its *learnability* section was built on `mean_pv` and is superseded; the target-distribution figures still stand |
| `drug_catalog.ipynb` | Cross-database compound harmonization (CTRP / GDSC / PRISM / DrugBank). ⛔ **Needs a rewrite, not a re-run:** it hardcodes absolute `/Users/...` paths, never uses `PipelinePaths`, reads CTRPv2's retired `v20.*` tables, and its `../data/*` inputs are not in the repository |

### `analysis/evaluation/` — is the number real?

Both are **temporarily** here rather than in `5_evaluation`, because both raise on their first cell.
They move up when they run.

| Notebook | Question it answers |
|---|---|
| `dreval_benchmark.ipynb` | **How strong is this by the field's standard?** Our data and model through the real **DrEval** package (`drevalpy` 1.5.1): their LCO splits, their baselines, their metrics. ⛔ Broken twice over — see `5_evaluation` |
| `diagnostics.ipynb` | The drug-selection gate defect, the proliferation test, and result dispersion across folds *and* drugs. ⛔ Hardcodes the removed `'auc'` |

### `archive/` — nothing here is load-bearing

**Two grounds for being here, and only two.**

1. **Nothing documented depends on it.** No step file cites it. This is the original test, applied by
   checking which notebooks the step files actually reference — not by judging how interesting they
   look.
2. **It cannot be re-run**, because the data or the target it reads no longer exists (added
   11.08.2026). A notebook can land here under this ground *while still being cited*, since archiving
   does not retract the numbers it produced. Where that happens, the citing step file must say the
   notebook is archived and why, so a reader who follows the path is not left thinking it is runnable.

**Their imports were deliberately not updated (12.08.2026).** `scripts/` was reorganized that day —
`layout`, the fetchers, `pubchem`, `cellosaurus`, `drug_annotation` and `gene_symbols` left
`scripts/preprocessing/` — and every live notebook, script and document was repointed. Six notebooks
here still import the old paths, and they keep them: an archived notebook is a record of a run, and
rewriting its imports would make it look runnable when ground 2 says it is not. The current location of
any module is in [CLAUDE.md](../CLAUDE.md) under *Where things live*.

| Notebook | Why it is here |
|---|---|
| `1_preprocessing.ipynb` | **Ground 1**, archived 12.08.2026. The old stage 1: §A recomputed the 512-d PCA for both built variants, §B built the HVG-sweep variants. Superseded by [`1_data`](1_data.ipynb) + [`3_representations`](3_representations.ipynb), which split preprocessing at the point [`2_drug_selection`](2_drug_selection.ipynb) needs. It drove `run_preprocessing.py`, archived the same day, so it cannot run as written. Its §B is the source for `analysis/qc/hvg_sweep_build.ipynb` |
| `scdrugatlas_exploration.ipynb` | Explores **scDrugAtlas**, a data source that was evaluated and [rejected](../docs/steps/corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources). Kept as the record of that decision. *(Long mislabelled in the docs as SCP542 exploration — "scDA" is scDrugAtlas.)* |
| `learnability_filter.ipynb` | The kill/spare gate that took 545 drugs to 10. Archived 12.08.2026: the criterion [measured potency, not rankability](../docs/steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability), and the [rebuilt panel](../docs/steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026) uses no statistic of our labels at all |
| `learnable_subset_training.ipynb` | PCA vs scGPT on that subset — a best-case diagnostic, never a generalization number. Archived 12.08.2026 with the gate that produced the subset |
| `panel_distributions.ipynb` | Response distributions and the density-weighting design on the [void 8-drug panel](../docs/steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it). Archived 12.08.2026. ⚠️ It also held the only justification for `alpha=0.5` and `cap=3` in `scripts/training/density_weighting.py`, so **audit 09 must re-derive those or drop the weighting** |
| `ctrp_prism_repurposing.ipynb` | CTRP→PRISM repurposing and clinical-phase mapping. Read-only, writes no artifact, and no step depends on it. Worth knowing it exists: it is the only notebook that loads `GDSC2_fitted_dose_response_27Oct23.xlsx`, which the "externalize the spread requirement" task will need |
| `target_comparison.ipynb` | **Ground 2**, archived 11.08.2026. `mean_pv` vs `auc` vs `auc_z` at K=10 and K=545. All three targets were removed that day, so `layout.CTRP_SCORES` rejects every one it asks for and it cannot run. Its conclusion — retire `auc_z` — is in [Corrections](../docs/steps/corrections-and-dead-ends.md#auc_z-as-the-training-target) |
| `replicate_variation.ipynb` | **Ground 2**, archived 11.08.2026. How far apart two screenings of the same (cell line, drug) fall — it existed because the pipeline *averaged* them, and it no longer does. Its measurements are kept in [Step 01](../docs/steps/01-datasets-and-harmonization.md#genuine-repeats-are-averaged-and-they-disagree-more-than-the-targets-own-spread-10082026), reframed as the evidence for the target switch |
| `ablations_and_rescue.ipynb` | **Ground 2**, archived 11.08.2026. *Why did the 545-head model fail, and what fixes it?* — the implicit σ²-weighting of the loss, the causal rescue test, the model-knob ablations, the ridge control. The whole argument is a head-to-head between `auc` and `auc_z`, both removed, so it cannot run. Re-wiring it to `auc_cc` / `ln_ic50_cc` would not preserve the argument; it would ask a different question. ⛔ **Its ablation and ridge tables are void (12.08.2026):** its `oof()` early-stopped on the fold it scored, so the MLP rows are flattered and the ridge row is not — [Corrections](../docs/steps/corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning). The within-notebook *rescue* ranking survives, because every arm there carries the same leak |

Two things are **not** grounds for archiving. **Superseded conclusions** — the 8-run matrix's results are
void but its harness is the only way to re-run it, so it survives as `4a_percell_training` §B rather than going
here. **A discredited criterion** — `2_drug_selection`'s predecessors were built on one, but it is the
documented record of how selection was done, and the rebuilt panel replaced it rather than deleting the
record.

## Re-running

Every training notebook has a **`RETRAIN` flag, default `False`** — it then loads the saved CSVs from
`outputs/` and only redraws the figures, in seconds. Set it `True` to refit.

Fits use `TrainConfig(epochs=25)`: across 36 recorded runs the best epoch was **median 6, max 11**, and
early stopping (patience 10) never came close to 25 — the earlier cap of 50 only cost wall-clock.

**The notebooks build the data** — stages 1 and 3, by calling `scripts/preprocessing/pipeline.py`. There
is no CLI to run instead; `run_preprocessing.py` was archived on 12.08.2026. To rebuild without a
browser, execute the stages headless:

```bash
uv run jupyter nbconvert --execute --inplace notebooks/1_data.ipynb
uv run jupyter nbconvert --execute --inplace notebooks/3_representations.ipynb
```

`--score` defaults to `auc_cc`; each score writes its own targets h5ad, so `ln_ic50_cc` does not
overwrite it. Set `SCORE` in the stage-3 bootstrap cell to switch.

> ⚠️ Until 12.08.2026 this section read *"Data is built by the CLI, not by the notebooks"*, which
> contradicted three statements above it and the header cell of the notebook that has built the data
> since it was written. What was true underneath it is narrower and is now in stage 1: the `hvg5000`
> targets h5ad was originally built by a direct CLI call, because no notebook covered `fetch → targets`
> at the time. One does now.

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
