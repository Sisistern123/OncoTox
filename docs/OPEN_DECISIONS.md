# Open decisions

Choices that are Selin's and are **not taken**. Opened 13.08.2026, during the wrap-up run.

Each entry states the choice as a choice, what each option implies, the assumption underneath, and —
where there is one — a reading marked as a reading. Nothing here is a recommendation that has been
acted on. When one is settled, move it into the file that owns the topic (`docs/TODO.md` for review
items, `docs/steps/` for the scientific record) and delete it here, so this file never becomes a
second place where a decision lives.

Ranked by what it blocks for the 14.08.2026 lab meeting.

---

## 1 · `worktree-docs-report-post-rerun` — merge, discard, or push

**Blocks:** nothing tonight. It blocks the repository being in a final state.

**The choice.** The branch holds two commits that are not on `main`:

| commit | subject |
|---|---|
| `7077096` | docs: retire the two stop-the-world banners — the rerun happened and Q1 reversed |
| `49611ef` | report: write the results from the rerun — both questions, with what each does not settle |

They describe the **previous** re-run — the one-seed 4a that Gate 5 recorded
(`docs/gate5-rerun-report.md`, execution log row 5, *"⚠️ One seed"*). The artifacts they were written
against were replaced on 13.08.2026 at 17:41 by `02f0fe6`, which committed a three-seed 4a with a
corrected early-stopping rule.

**Options.**

- **Merge** — pulls prose describing superseded numbers into `main`. The banner-retirement half of
  `7077096` is work that still needs doing, but `49611ef`'s results section is written against
  numbers that no longer hold.
- **Discard** (`git branch -D worktree-docs-report-post-rerun`) — loses both, including the banner
  wording, which would then be rewritten from scratch.
- **Push the branch without merging** — preserves it on the remote without putting it on `main`.
  Costs a remote branch you would later clean up, and puts superseded conclusions somewhere they can
  be read as current.
- **Leave it local, unmerged** — what is in force now. The commits survive in the shared `.git`
  (branch refs are not stored in the worktree, so they outlive the worktree directory).

**The assumption underneath.** That `49611ef`'s results section cannot be salvaged by editing. It
would have to be checked paragraph by paragraph against the three-seed artifacts; whether that is
cheaper than rewriting is a judgement about the prose, which has not been made.

**My reading, as a reading.** Do not push and do not merge. Discard `49611ef` and rewrite the
banner retirement directly on `main`, because the banners have to be corrected tonight anyway and a
merge would conflict with that work. But the branch is costing nothing where it sits, so leaving it
is equally defensible and is what I have done.

**Not done by me:** nothing on this branch was merged, pushed, deleted, or edited.

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

**No reading offered.** This is the choice CLAUDE.md names first, and I have no measurement that
would inform it.

---

## 3 · The alpha axis — the decision is not booked as taken

**Blocks:** whether the density-weighting sweep can be described as closed.

**The choice.** Whether the axis (`alpha` in {0, 0.5, 1}) is closed, and at which level.

**Why it is open.** The earlier report that *"alpha=0 won"* was wrong: it read the guard column (delta
against the null model) as the item-9A criterion. On mean per-drug Spearman, from
`notebooks/outputs/panel/panel_leaderboard.csv`:

| arm | X_pca | X_scGPT |
|---|---|---|
| MLP mse alpha=0 | 0.2541 | 0.2009 |
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

**My reading, as a reading.** The axis is not closed. +0.021 on one representation with no band
attached is not enough to close it. Computing the seed band from `panel_oof_predictions.csv` would
settle it and is cheap, but it needs the aggregation rule (decision 2) first.

---

## 4 · `RUN_SECTION_C`, `RUN_SECTION_D`, `RUN_SECTION_E` are all committed as `True`

**Blocks:** nothing scientific. It sets the cost of reproducing `4a`.

**The choice.** Whether the three section flags in `notebooks/4a_percell_training.ipynb` stay `True`
in the committed notebook. Each section's own banner says to set it back to `False` once the question
is answered; all three are `True`, because that is what produced the stored outputs.

**What each option implies.** Left at `True`, a top-to-bottom run of `4a` retrains
**90 (A) + 60 (C) + 60 (D) + 120 (E) = 330 fits**, of which 240 answer questions already answered.
Set to `False`, a top-to-bottom run reproduces section A only and the other three print their skip
message — so the notebook stops regenerating the artifacts its own markdown discusses.

**The assumption underneath.** That "reproducible" means a fresh clone regenerates the artifacts,
rather than that the notebook is cheap to re-run.

**My reading, as a reading.** Leave all three `True`. A notebook whose committed flags do not
reproduce its committed outputs is the defect class this review has been finding all week; the cost
is wall-clock, which is recoverable, and the alternative costs reproducibility, which is not.

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
