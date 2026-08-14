# Open decisions

Choices that are Selin's and are **not taken**. Opened 13.08.2026, during the wrap-up run.

Each entry states the choice as a choice, what each option implies, the assumption underneath, and —
where there is one — a reading marked as a reading. Nothing here is a recommendation that has been
acted on. When one is settled, move it into the file that owns the topic (`docs/TODO.md` for review
items, `docs/steps/` for the scientific record) and delete it here, so this file never becomes a
second place where a decision lives.

Ranked by what it blocks for the 14.08.2026 lab meeting. ⚠️ *That date has passed; the one
remaining entry (§7) blocks nothing, so the ranking no longer orders anything.*

---

## 1 · ~~`worktree-docs-report-post-rerun`~~ — SETTLED 14.08.2026

**Discarded (Selin).** The worktree was removed and the branch deleted. What it contained, why it was
discarded rather than merged, and how to get it back live in
[Corrections § Dead ends](./steps/corrections-and-dead-ends.md#worktree-docs-report-post-rerun--a-results-section-written-against-the-previous-run).
Numbering is kept so earlier references still resolve.

---

## 2 · ~~The aggregation convention~~ — SETTLED 14.08.2026

**Line-level, pooled (Selin).** Predictions are averaged to the cell line, then one Spearman per drug
across all held-out lines, then the mean over drugs.

**Why, recorded because the alternatives were measured and are defensible.** The label is per
(cell line, drug), so one point per line is one point per label — cell-level scoring adds ~300 points
per label that carry no independent information, costs resolution to the within-line scatter, and
weights lines by cell count over a 36× range (56–1,990). Pooling estimates from all ~153 held-out
lines at once rather than averaging five correlations of ~30, where Spearman's small-sample upward
bias is visible in the data: **every arm scores higher under per-fold**, both representations, all six
arms.

**What it costs: nothing.** `5_evaluation` §1.8's `order` was already computed this way — verified
arm by arm against `panel_aggregation_comparison.csv`, identical to four decimals for all twelve.
No number moves and no code changes.

> ⚠️ **A note about the artifact — and the correction is the point (14.08.2026, 19:20 JST).**
> `panel_aggregation_comparison.csv` went dirty in the working tree during the consolidation sweep:
> one row, `X_pca / mlp / alpha=0 / mse / seed 42`, moved in the fourth decimal
> (`0.24210 → 0.24397` on the first column). **It has been restored to `HEAD` and the tree is clean,
> so the verification above stands exactly as committed.**
>
> ⛔ **I first recorded here that a live Jupyter kernel had written it. That was wrong, and I caused
> it.** `scripts/gate/verify_main.sh`'s module-import check imports every module under `scripts/`
> except `archive` and `gate`, and nine files in `scripts/evaluation/` are straight-line scripts with
> no `if __name__ == "__main__"` guard — so **importing them runs them**.
> `aggregation_comparison.py` writes this very file at top level. Running the post-merge gate
> therefore rewrote it. Proven twice: the gate's own log carries that script's
> *"skip … no matching rows in panel_oof_predictions.csv"* output, and after the file was restored it
> **stayed** restored while the kernel kept running.
>
> **Nothing about the decision changes.** The row that moved is the single arm §3 already disqualifies
> as unstable, so even had it been a genuine re-run it would not bear on the convention. What changed
> is the gate: `evaluation` is now excluded from its import check, and the cost of that exclusion is
> stated at the check itself.

**What it does and does not settle.** It fixes which arm has the highest `order` on `X_pca`
(**α=0.5/mae, 0.2824**). It does **not** settle item 9A, which *blocks* that arm on `values` and
`spread_slope` — raising the density exponent buys ranking and pays for it in calibration and
absolute error. Convention and decision rule are separate questions.

**And it does not touch Q1**, which leads under all four conventions in all 24 cells — the whole
point of having measured them.

Full table: [Step 05](./steps/05-multitask-results.md#the-aggregation-convention--q1-survives-it-the-loss-ranking-does-not).

---

## 3 · ~~The alpha axis / item 9A~~ — SETTLED 14.08.2026, and re-anchored

**The rule cannot select a winner, and the blocker is identified.** Full record in
[Step 05 §Item 9A](./steps/05-multitask-results.md#item-9a--settled-14082026-the-rule-cannot-select-a-winner-and-why).

Applied under the line-level pooled convention (§2), the rule returned **two different answers** on
two runs — "no challenger wins" at an incumbent `order` of 0.2541, and "`α=0`/`mae` wins" at 0.2473.
**The verdict flips inside the incumbent's own measured band**, 0.2450–0.2541 over eight executions.

The reference point is the problem, and it is a poor one on two measured grounds:

1. `α=0`/`mse`/`X_pca` is the **only** unstable configuration in the sweep — twelve of fourteen
   arm × rep rows are identical to six decimals across all eight runs, including the ridge control.
2. It is **last under every scoring convention** — all four rank `α=0` bottom on both losses.

It is also the arm with the earliest median `best_epoch` of all twelve (**1.0**; 8 of its 15 fits stop
at epoch 1). It barely trains, so its score sits near the head-bias initialisation.

**A device warm-up was tried and does not fix it** — two runs gave 0.2525 and 0.2480. Reverted
(`e6c087d`); reapplicable from `664f3e8` if the separate, smaller first-fit effect is ever worth
chasing.

### ✅ Closed 14.08.2026 — re-anchored on `α=0`/`mae`/`X_pca`

Selin took the reading offered here. The incumbent is now one of the twelve arms identical across all
eight executions, so the verdict is stable, and it is written into `5_evaluation` §1.8 as an explicit
`INCUMBENT_ARM` rather than falling out of iteration order.

**Result: no challenger wins — all thirteen blocked.** The unweighted MAE arm on `X_pca` stands;
density weighting does not displace it. `α=0.5`/`mae` has the sweep's highest `order` (0.2824 against
0.2617) and is blocked on `values` and `spread_slope`.

⚠️ It changes what the rule asks — challengers are now judged against `mae` — and that is stated
where it is used rather than left implicit.

---

## 4 · ~~The three section flags in `4a`~~ — SETTLED 14.08.2026

**All three stay `True` (Selin).** A full top-to-bottom run of `4a` is therefore 330 fits, about
1 h 24 min, and it reproduces every artifact the notebook's own markdown discusses. The reasoning is
recorded where it applies — in the §C, §D and §E banners of
`notebooks/4a_percell_training.ipynb` — rather than here.

---

## 5 · ~~`5_evaluation` sections 2 and 3~~ — SETTLED 14.08.2026

**Written and executed (Selin: promote everything each notebook computes).** §2 absorbs
`analysis/evaluation/diagnostics` — all five CSVs, the three figures named rather than re-rendered.
§3 absorbs `analysis/evaluation/dreval_benchmark` — the per-fold table, DrEval's own baseline
leaderboard, and the normalized read on the pipeline's own predictions.

Both **read**; neither recomputes. Every number keeps the single owner that wrote it, so nothing is
stated in two places. Nothing was selected on your behalf: all artifacts each notebook produces are
promoted, in the order those notebooks compute them.

Two things the absorption surfaced that were not visible before:

- **DrEval does not separate the arms, and now says so numerically.** The spread *within* `X_pca`
  across folds is **0.1245**; the largest gap *between* the arms in any fold is **0.0918**.
- **Under DrEval's own normalization applied to our predictions, `X_pca` leads `X_scGPT` at every
  α** (0.3247/0.3310/0.3169 against 0.2910/0.2743/0.2559). The external check does not contradict
  Q1 — the earlier reading that it "pointed the other way" rested on one fold of one run.

---

---

## 6 · ~~The two arms are not optimised comparably~~ — LARGELY DISSOLVED 14.08.2026

**Opened this morning on the reading that the 104× input-scale gap confounded Q1. That reading was
tested and is wrong.**

`X_scGPT` was rescaled by 103.4× to `X_pca`'s exact magnitude and the α=0/mse arm re-run over three
seeds and five folds. Mean `best_epoch` went **6.73 → 7.13** — unchanged, with several folds
identical fold-for-fold. **AdamW normalises each step by a running second moment of the gradient, so
it is approximately invariant to a uniform input rescale**, which makes review item 4A's premise
(*"one learning rate is not one setting"*) a statement about SGD-like updates rather than about this
optimizer.

**What that dissolves.** The scale gap is **not** a confound on Q1, and §C's capacity comparison is
not compromised by it. The three options originally listed here — standardise the inputs, per-arm
learning rates, report the confound — were answers to a problem that does not exist in the form
stated.

**What remains, and it is genuinely smaller.** The 104× difference is real and worth stating as a
difference between the arms; it simply does not drive the training-regime crossover or the score gap.
The crossover itself is an interaction between **capacity and representation**
([Step 05](./steps/05-multitask-results.md#training-dynamics-do-not-explain-the-gap--best_epoch-does-not-track-score)),
and `best_epoch` does not track score in any case.

**Still open, if you want it:** whether to standardise the representations anyway, on the general
principle that two arms being compared should reach the optimizer alike. That is now a tidiness
argument rather than a correctness one, and it would still change every number, so my reading is: not
before the talk, and not urgent after it.

**My earlier reading — that this was "the strongest remaining threat to Q1" — is withdrawn.** It was
stated with a hypothesis marked as a hypothesis, and the test refuted it within the day.

---

### ⚠️ Below: the body item 6 had BEFORE it dissolved — superseded, kept as the record (marked 14.08.2026)

**None of what follows is live.** It is the reasoning item 6 was opened on, and the paragraph above
withdraws its conclusion by name. It was left in the file unheaded, beneath a table whose header row
had been lost, so the section ran straight into item 7 and read as current — including *"the strongest
remaining threat to Q1"*, which is the exact sentence line 150 retracts. Preserved rather than deleted
because the reasoning is what the test refuted, and a retraction without its claim cannot be checked.
**The header row is restored below; the three options and the closing reading are superseded by the
dissolution above.**

Median `best_epoch`, by architecture — re-derived 14.08.2026 from
`notebooks/outputs/panel/panel_arch_folds.csv` (60 fits), which is the artifact that carries
`best_epoch`; `panel_arch_summary.csv` does not:

| arch | `X_pca` | `X_scGPT` |
|---|---|---|
| trunk (128,64) | **1** | 8 |
| linear | **12** | 2 |

A crossover, so it is not "PCA trains fast". `X_pca` reaches the optimizer at **104×** the magnitude of
`X_scGPT` (median per-dimension sd 1.1062 against 0.0107,
`notebooks/outputs/diagnostics/input_scale.csv`) under **one shared learning rate**. Large inputs with
high capacity overfit inside one epoch; tiny inputs with low capacity cannot move off the head-bias
initialisation. Both are optimisation artefacts, not statements about the representations.

**What it costs.** The α=0/mse `X_pca` arm peaks at **epoch 1** — barely past its initialisation —
which is why it is the one arm that does not reproduce (§3). And **§C's "capacity does not carry Q1"
was measured correctly but rests on an uncontrolled difference**: the two arms are not receiving
comparable optimisation, so a capacity effect and a scale effect cannot be separated in that design.

**This is review item 4A's second half**, which has never been tested: *"if one arm reaches the
optimizer with values ~78× larger than the other, one learning rate is not one setting, and an arm can
look worse for a reason that has nothing to do with the representation."* The measured factor is now
104×, not 78×.

**The options.**

- **Standardise each representation** (z-score to unit variance) before the model. Makes the learning
  rate mean the same thing for both. Changes every number and requires a full re-run.
- **Per-arm learning rate**, tuned so each peaks in a comparable epoch range. Keeps the inputs as they
  are, but introduces a second per-arm setting to justify.
- **Leave it and report the confound.** Costs nothing, and the Q1 margins stay as measured with a
  stated qualification.

**Label supply is excluded as the cause, tested 14.08.2026.** §E varies only the label budget with
architecture, α, loss and input held fixed. `Spearman(n_label_lines, best_epoch)` is **−0.274**
(p = 0.034) for `X_pca` and **+0.463** (p = 0.0002) for `X_scGPT` — opposite signs, both significant.
Label supply moves where training peaks, but in opposite directions per representation, so it cannot
be the common cause. Full table in
[Step 05](./steps/05-multitask-results.md#why-training-peaks-so-early--label-count-excluded-capacity--scale-implicated-14082026).

**The assumption underneath.** That the Q1 ordering would survive equalisation. Untested — and it is
the one open confound that could plausibly move it, since it is the only known difference between the
arms that is not the representation itself.

**My reading, as a reading.** Report the confound now (option 3) and standardise afterwards, because
option 1 changes every number in the project and there is no time to re-verify them before the talk.
But say it out loud in the talk: it is the strongest remaining threat to Q1, and it is better said by
you than found by someone in the room.

---

## 7 · ~~`input_dropout`~~ — MOVED TO OUTLOOK 14.08.2026 (Selin)

**Not settled by choosing a value, and not left open either: demoted.** Selin's ruling — *"if it is
unnecessary, keep it as an outlook."*

**The sweep she remembered wanting has already been run**: six values — 0, 0.02, 0.05, 0.10, 0.20,
0.30 — three seeds each, both arms, 36 fits, committed as
`notebooks/outputs/diagnostics/input_dropout_test.csv` and `input_dropout_sweep_extra.csv`. So there
is no experiment outstanding. What the sweep shows: **no interior optimum** (`X_pca` peaks at 0.02,
0.2753, against the shipped 0.10's 0.2608; `X_scGPT` is flat at 0.0020 across the whole range, inside
one seed sd).

**Nothing changes.** `input_dropout` stays at **0.1** on both arms, which is what every reported
number was measured at and is the **conservative** setting — at `X_pca`'s best rate the Q1 margin
would be ≈0.046 rather than the reported 0.0317.

**Where it now lives, as future work rather than an open decision:**
`report/sections/06_limitations_and_outlook.tex` and `docs/final_presentation.md` §7. The finding
itself is owned by [Step 05](./steps/05-multitask-results.md).

**Why it is not a decision.** It blocks nothing, the interim course is implemented (0.1, named
wherever a margin is quoted — `report/sections/04_results.tex`), changing it re-runs every number,
and Q1's direction survives every rate in the sweep.

---

## Nothing is currently open

Every entry above is settled, dissolved or demoted. **This file is empty of open decisions as of
14.08.2026** — which is the state it should be left in, not a sign it has been forgotten. Add the
next one as §8.
