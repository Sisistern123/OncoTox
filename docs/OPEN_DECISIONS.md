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

## 2 · The aggregation rule

**Blocks:** the interpretation of every Spearman in the deck.

**The choice.** How predictions are aggregated before scoring: **per cell or per cell line**, and
**pooled across folds or computed per fold and averaged**. Recorded as open and yours in
`docs/TODO.md` items 9A and 11.

**What is in force in the code today**, so the gap between the open decision and the running code is
visible: `scripts/training/cv.py::line_level` aggregates to the **cell line**, and `4a`'s §C/§D score
**per drug, pooled over the out-of-fold predictions**, then take the mean over drugs. That is what
produced `panel_arch_summary.csv` and `panel_heads_summary.csv`, and therefore every §C/§D number.

**What each option implies.** Per-cell scoring weights cell lines by how many cells they contribute
(56–1,990 cells per line, an ~35x range), so large lines dominate. Per-line scoring gives each of the
~153 held-out lines one vote, which matches the fact that the label is per line. Pooling across folds
mixes folds of different difficulty into one correlation; per-fold-then-average separates them but
has ~30 lines per fold to compute a correlation from.

**The assumption underneath.** That the choice does not change the *ordering* of the arms, only the
magnitude. That has not been tested — no run has scored the same predictions both ways.

**✅ MEASURED 14.08.2026 — the choice is now informed, and it splits in two.**
All four conventions computed on the same predictions:
`notebooks/outputs/panel/panel_aggregation_comparison.csv`, written up in
[Step 05](./steps/05-multitask-results.md#the-aggregation-convention--q1-survives-it-the-loss-ranking-does-not).

- **For Q1 the choice does not matter.** `X_pca` leads `X_scGPT` in all 24 cells (6 arms × 4
  conventions). Only the size moves, by up to 2×.
- **For the loss comparison it decides the answer.** The best arm is α=0.5/mae, α=0.5/mse,
  α=0.5/mae or α=1/mse depending on which convention you pick. α=0 is last under all four.

**So what is actually still yours** is narrower than the original entry suggested: not "which
convention is right in general", but **which convention the loss comparison is judged under** — and
whether that comparison should be made at all, given it also fails to survive a re-run (§3).

**My reading, as a reading.** Line-level, pooled — the convention already in force. The label is per
cell line, so a line is the natural unit and cell-level scoring spends resolution on within-line
scatter that carries no label information; and pooling keeps one correlation over ~153 lines rather
than six over ~30, where the per-fold inflation above comes from. But this is a preference between
defensible options, and it is the one that decides the loss arm, so I have not acted on it.

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

**What I would put to you, and it is not the α question any more.** Before α is decided, the
comparison needs an error bar that includes **re-execution**, not only seeds. Everything in this
project is currently reported with a seed band; `4a` §A now demonstrates that a fresh run of the same
seeds moves one arm by more than a third of that band. Options: report a band over N repeated
executions (costly — `4a` is 1 h 23 min a run); find and fix the non-determinism (review item 10
recorded it as *not* reproducing under current code, which tonight refutes); or state the loss
comparison as unresolved at this noise level and stop there.

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

## 5 · `5_evaluation` sections 2 and 3 do not exist

**Blocks:** whether `5_evaluation` can be called complete.

**The choice.** Whether to write them before the meeting.

Gate 5's execution log, row 15, records them as *"not yet written"*: they are meant to absorb
`analysis/evaluation/diagnostics` and `analysis/evaluation/dreval_benchmark`, both of which have run.
Writing them is new code, and the absorption involves choosing which of each notebook's results are
promoted into the evaluation chain — a selection, and therefore yours.

**My reading, as a reading.** Do not write them tonight. Both notebooks stand on their own and their
outputs are committed; folding them in is a structural improvement that changes no number, and it is
the kind of work that goes wrong when done against a deadline without the person who owns the
selection.
