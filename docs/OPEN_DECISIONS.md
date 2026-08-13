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

