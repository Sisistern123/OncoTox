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

## 3 · The alpha axis — the decision is not booked as taken

**Blocks:** whether the density-weighting sweep can be described as closed.

**The choice.** Whether the axis (`alpha` in {0, 0.5, 1}) is closed, and at which level.

**Why it is open.** The earlier report that *"alpha=0 won"* was wrong: it read the guard column (delta
against the null model) as the item-9A criterion. On mean per-drug Spearman, from
`notebooks/outputs/panel/panel_leaderboard.csv`:

| arm | X_pca | X_scGPT |
|---|---|---|
| MLP mse alpha=0 | 0.2541 → 0.2473 ⚠️ unstable | 0.2009 |
| MLP mse alpha=0.5 | **0.2754** | 0.1927 |
| MLP mae alpha=0 | 0.2617 | 0.2403 |
| MLP mae alpha=0.5 | **0.2824** | 0.2395 |

alpha=0.5 beats alpha=0 on `X_pca` by **+0.0213** (mse) and **+0.0207** (mae) — in both losses. On
`X_scGPT` it does nothing or slight harm. The closure decision rested on the wrong summary and is
**not recorded as taken**.

**What is settled and does not need re-deciding:** alpha=1 puts `X_scGPT` below the per-drug
constant, so that level is out (`4a` cell 38's comment; `docs/steps/corrections-and-dead-ends.md`,
the inverse-density entry).

**The assumption underneath.** That the leaderboard's point estimates are comparable without a band.
`panel_leaderboard.csv` reports one number per arm and carries no seed band, although the underlying
`panel_oof_predictions.csv` has three seeds — so a band **could** be computed and has not been.

**⚠️ UPDATED 13.08.2026, later the same evening — the rule has now been run.**
`5_evaluation` §1.8 was executed for the first time and applied item 9A's rule
(`notebooks/outputs/panel/panel_metrics.csv`, three seeds). Against the `alpha=0` / `mse` / `X_pca`
incumbent, **all thirteen challengers are blocked**. `alpha=0.5` / `mae` / `X_pca` has the highest
`order` of any arm in the sweep — **0.2824**, above the ridge control's 0.2767 — and still fails, on
`values` and `spread_slope`. Raising the density exponent buys ranking and pays for it in calibration
and absolute error.

**What that leaves you to decide.** The rule returns "no challenger displaces the incumbent", which
is an answer to item 9A. Whether that *closes* the axis is still yours, because it depends on
something the rule cannot settle: **whether a guard failure should veto a real gain in the primary
quantity, or merely be reported alongside it.** As written, the rule vetoes. If you would rather rank
on `order` and report the guards, `alpha=0.5`/`mae` wins on `X_pca` and the answer reverses.

**⛔ CORRECTED AGAIN, 23:35 the same evening — the rule gave a different answer on a re-run.**
`4a` was then executed top to bottom in a fresh kernel. The incumbent arm's `order` fell from 0.2541
to **0.2473** — same code, same seeds, same inputs — and on the re-evaluated rule **`alpha=0` / `mae`
/ `X_pca` now WINS**, where four hours earlier all thirteen challengers were blocked.

**So the axis is not answered, and the reason is worse than an open choice.** The gap the verdict
turns on is +0.0076 in one run and +0.0144 in the other, and the run-to-run instability that moved it
is of the same size as the effect. A decision rule cannot resolve a difference smaller than the
variation in its own inputs.

**✅ The band was measured, 14.08.2026 — five executions.** The instability is far narrower than two
runs suggested: **eleven of twelve arms are identical to six decimals across all five**, and within
the twelfth it is a single fit — fold 1, seed 42, the **first fit executed in the process**. 179 of
§A's 180 fits reproduce exactly. Full record in
[Step 05 §Reproducibility](./steps/05-multitask-results.md#reproducibility--measured-over-five-executions-and-it-is-one-fit-not-the-pipeline).

**That changes what is at stake here.** The pipeline is not unreliable; one fit is. But that fit is
the α=0/mse `X_pca` arm, which is item 9A's **incumbent**, so its 0.0091 range still decides the
verdict: at 0.2541 nothing clears the bar, at 0.2450–0.2490 `alpha=0`/`mae` does. Four of five runs
fall below the flip point.

**Three ways out, and the cheapest is now the most attractive.**

1. **Fix the first-fit instability.** Predicted mechanism: something warm on the second fit is cold on
   the first. The test is one §A run with the arm order reversed (~38 min) — if the wobble follows
   the first *position*, a throwaway warm-up fit before the grid removes it and the comparison becomes
   exactly reproducible. **Not run**, and the mechanism is a hypothesis.
2. **Report a band over executions.** Now cheap for the honest version: only the one arm needs it.
3. **State the loss comparison unresolved** and stop.

**My reading, as a reading.** State it as unresolved and say why. Both verdicts are defensible
readings of the same design, which is exactly what "unresolved" means, and presenting either one as
*the* answer would be presenting a coin-flip as a finding. The guard-veto question below is still
real, but it is downstream of this.

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

