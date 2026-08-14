# `scripts/archive/` — superseded code, kept readable, not runnable

Same rule as [`notebooks/archive/`](../../notebooks/README.md): **nothing here is load-bearing.** A file
lands here when the thing it computed is no longer part of the analysis, and it stays because the
reasoning inside it is part of the record — deleting it would leave the write-ups in
[`docs/steps/corrections-and-dead-ends.md`](../../docs/steps/corrections-and-dead-ends.md) pointing at
nothing.

Assume every file here is **broken**: the target moved to `auc_cc` on 11.08.2026 and the artifacts these
scripts read were never regenerated. Do not run them; read them.

⚠️ **One exception, so the blanket claim above is not read as covering it.** `run_preprocessing.py` was
archived on 12.08.2026 because it was **superseded, not because it broke** — it ran correctly on the day
it was moved. The distinction matters: a broken file is evidence of what went wrong, whereas this one is
a design that was replaced, and reading it as defective would misattribute the reason.

## `run_preprocessing.py` — the CLI orchestrator

**Archived 12.08.2026 (Selin). Replaced by [`scripts/preprocessing/pipeline.py`](../preprocessing/pipeline.py).**

**What it was.** A single `argparse` entry point that ran the six preprocessing steps in order —
`fetch → convert → scgpt → targets → splits → pca` — from a `STEP_ORDER` list, with `--start-at` to
resume partway and per-step guards against clobbering expensive artifacts.

**Why it went.** The notebooks were renumbered into five end-to-end stages on 12.08.2026, which split
preprocessing across two of them: `1_data` runs `fetch` and `convert`, `3_representations` runs the
remaining four. The step *order* therefore became the notebooks' property, and a second copy of it inside
a CLI is a second thing that has to be kept in step with them — the failure this repository has already
had with paths and targets. What was not plumbing was kept: every guard, every between-step precondition
and the scGPT subprocess bridge are now functions in `pipeline.py`, one per step, each owning its own
checks so a step cannot run out of order regardless of who calls it.

**What changed in the move, rather than being transcribed:**

- `guard_output` now raises `FileExistsError` instead of `SystemExit`, and its message names
  `overwrite=True` rather than the `--overwrite` flag that no longer exists (`scripts/layout.py`).
  `resolve_data_root` likewise raises `NotADirectoryError`.
- The response-table check moved from the `targets` step into `fetch`, because the two are no longer
  consecutive lines in one script but a whole stage apart.
- `scgpt()` has **no interactive fallback**. The original blocked on `input()` when given no
  interpreter; run headless via `jupyter nbconvert --execute` that hangs the kernel on an invisible
  prompt instead of failing.
- The single-drug steps are gone — see below.

**There is deliberately no CLI replacement.** The notebooks are the entry point; to run without a
browser, execute them with `jupyter nbconvert --execute`.

## `single_drug_dataset.py` — `ScGPTDrugDataset`

**Archived 12.08.2026 (Selin), with the whole chain that fed it.** Nothing called it. The two scripts its
docstring named as callers were deleted on 26.05.2026. `ctrp_to_h5ad` no longer writes
`viability_<drug>` / `train_mask_<drug>`, `create_splits.run_multi` is the only split writer, and
`add_pca.TRAIN_SPLIT_COLS` is down to `("split_ctrp",)` — which had been costing a second 512-component
train-only decomposition on every pipeline run for a column no code read. Full reasoning in the file's
own docstring. For a single drug, `MultiDrugDataset(..., drugs=["<one>"])` covers it and needs no
per-drug column at all.

## Removed rather than archived — `dreval_normalize.py`'s cell-line-effect diagnostic

**Deleted 12.08.2026 (Selin), not archived** — the one exception to the rule above, and deliberately so.

**What it was.** A **locally invented** metric that removed the **cell-line effect** from our own
out-of-fold predictions and re-scored them, reporting `rho_raw`, `rho_normalized` and
`rho_naive_baseline` per (heads × representation × drug). It asked whether a per-drug correlation was
drug-specific biology or merely *"this cell line is fragile"*, and it read held-out labels to do it, so
it was always a diagnostic rather than a predictor. A second local addition scored every model against
one common `auc` ranking regardless of which score it trained on.

**Why it was deleted rather than kept.** It has **no counterpart in DrEval's paper**, and it lived in a
file named after DrEval — an arrangement in which its output gets read as DrEval's metric. The paper
describes subtracting the `NaiveMeanEffectsPredictor` from truth and prediction and nothing more. The
standing instruction for this strand is that it *"should be as contained as the paper itself in
functionality, and nothing new"*, and an archived copy of a home-grown metric is still an invitation to
revive it without re-deciding it.

**Where it is if it is ever wanted:** commit `bf93084`, path `scripts/evaluation/dreval_normalize.py`
(`git show bf93084:scripts/evaluation/dreval_normalize.py`). It was introduced on 27.07.2026 to
reconstruct `notebooks/outputs/dreval/dreval_normalized.csv`, whose own producing code had already been
lost — so the metric is older than the file that implemented it.

**What is live instead.** `scripts/evaluation/dreval_normalize.py`, rewritten 12.08.2026 to apply
DrEval's normalization and nothing else: their `NaiveMeanEffectsPredictor` and their
`drevalpy.evaluation.evaluate`, no re-implementation, defaulting to `notebooks/outputs/panel/panel.csv`
on `auc_cc`.

**The open question the deletion does not settle**, for **audit 11 (Evaluation)**: under our
leave-cell-line-out splits, DrEval's normalization removes only the **drug** effect, because a held-out
line's effect is unseen and therefore zero. A synthetic predictor emitting nothing but
`mean + line effect + drug effect` — no drug-specific signal at all — would still score highly.
⚠️ *A figure of **0.98** stood here until 14.08.2026 and is removed: the 12.08.2026 audit cleared that
claim because no code produces it and no artifact records it, but reached only the report. The
mechanism stands; the number was never measured.* So *"is this drug-specific signal or general
fragility?"* is a real question that the paper's
metric does not answer under this split design. Audit 11 has to decide how to answer it, and reviving a
metric that reads held-out labels is only one of the options.

**Already broken when it went**, independently of all of the above: it built `PipelinePaths(..., "auc")`,
and `auc` stopped being a valid score on 11.08.2026, so it raised on construction. Its committed outputs
under `notebooks/outputs/dreval/` were computed on the retired target and the voided 8-drug panel, and
are void with them. `notebooks/analysis/evaluation/dreval_benchmark.ipynb` imports the removed module and
also hardcodes `'auc'`; it is untouched pending audit 11.

## `merge_gate.sh` — retired 14.08.2026 (Selin)

**Verified a branch against `main` before a merge. There are no branches.** `git branch -a` shows
`main` only, the workflow is commit-direct-to-main, and the last merge commit (`99637f4`, 13.08.2026)
is **140 commits** behind HEAD. A gate that cannot be run as intended is not a safety net; it is a
file that looks like one.

**What it did that `verify_main.sh` does not:** it diffed a branch tip against `main` — the
`coderefs.py` before/after comparison at its checks 4-6, which catches a reference that a *branch*
breaks. That capability goes with it, and nothing replaces it, because nothing produces branches to
compare.

**What survives:** every content check it called is still called by `scripts/gate/verify_main.sh`,
which runs against `main` as it stands and is the gate that is actually used. Recoverable from git
history if branching ever returns.

