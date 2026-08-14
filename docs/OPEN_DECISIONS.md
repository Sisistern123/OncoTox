# Open decisions

Choices that are Selin's and are **not taken**. Opened 13.08.2026, during the wrap-up run.

Each entry states the choice as a choice, what each option implies, the assumption underneath, and —
where there is one — a reading marked as a reading. Nothing here is a recommendation that has been
acted on. When one is settled, move it into the file that owns the topic (`docs/TODO.md` for review
items, `docs/steps/` for the scientific record) and delete it here, so this file never becomes a
second place where a decision lives.

Ranked by what it blocks for the 14.08.2026 lab meeting.

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

**What it does and does not settle.** It fixes which arm has the highest `order` on `X_pca`
(**α=0.5/mae, 0.2824**). It does **not** settle item 9A, which *blocks* that arm on `values` and
`spread_slope` — raising the density exponent buys ranking and pays for it in calibration and
absolute error. Convention and decision rule are separate questions.

**And it does not touch Q1**, which leads under all four conventions in all 24 cells — the whole
point of having measured them.

Full table: [Step 05](./steps/05-multitask-results.md#the-aggregation-convention--q1-survives-it-the-loss-ranking-does-not).

---

## 3 · ~~The alpha axis / item 9A~~ — SETTLED 14.08.2026, with one narrower question left

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

### The one question left, and it is yours

**Which arm should the rule be anchored to instead?** Any of the twelve stable rows makes it
evaluable. That is a decision about what the comparison is anchored to, not a measurement, so it has
not been taken.

**My reading, as a reading.** Anchor on `α=0`/`mae`/`X_pca` — it is the same α level, so the axis
still reads as "does weighting help", but on the stable loss; and `mae` was already measured to be
the better-behaved objective on this data. But note this changes what the rule *asks*, so it is not a
free substitution.

⚠️ **Until then, do not report "α=0 wins" or "MAE wins".** Both are readings of one arm's noise.

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

## 6 · The two arms are not optimised comparably — input scale against one learning rate

**Opened 14.08.2026.** Found while asking why `best_epoch` is so early.

**The observation.** Median best epoch, out of 50 available, patience 10:

| | `X_pca` | `X_scGPT` |
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

