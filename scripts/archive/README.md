# `scripts/archive/` — superseded code, kept readable, not runnable

Same rule as [`notebooks/archive/`](../../notebooks/README.md): **nothing here is load-bearing.** A file
lands here when the thing it computed is no longer part of the analysis, and it stays because the
reasoning inside it is part of the record — deleting it would leave the write-ups in
[`docs/steps/corrections-and-dead-ends.md`](../../docs/steps/corrections-and-dead-ends.md) pointing at
nothing.

Assume every file here is **broken**: the target moved to `auc_cc` on 11.08.2026 and the artifacts these
scripts read were never regenerated. Do not run them; read them.

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
`mean + line effect + drug effect` — no drug-specific signal at all — still scores normalized Spearman
**0.98**. So *"is this drug-specific signal or general fragility?"* is a real question that the paper's
metric does not answer under this split design. Audit 11 has to decide how to answer it, and reviving a
metric that reads held-out labels is only one of the options.

**Already broken when it went**, independently of all of the above: it built `PipelinePaths(..., "auc")`,
and `auc` stopped being a valid score on 11.08.2026, so it raised on construction. Its committed outputs
under `notebooks/outputs/dreval/` were computed on the retired target and the voided 8-drug panel, and
are void with them. `notebooks/result_evaluation/dreval_benchmark.ipynb` imports the removed module and
also hardcodes `'auc'`; it is untouched pending audit 11.
