# Step 04 — Single-task results (paclitaxel) & the data-leak fix

*Part of [OncoTox project progress](../project_progress.md). Covers: the single-task paclitaxel
baseline, and the random-split data leak that the grouped split fixed — the methodological result this
page exists for.*

This is plan-Phase-2 (single-task continuous regression). Model/training design is in
[Step 03](03-model-and-training-design.md).

> ⛔ **Every number this page used to carry is void and cannot be regenerated.** All of them were
> trained on **`mean_pv`**, which was removed with its reader code on 11.08.2026 when the target moved
> to DrEval's reprocessed CTRPv2
> ([Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)).
> `--score mean_pv` now raises, so this page's earlier instruction to reproduce it that way was not
> merely stale but false. The measurements are kept as a record in
> [Corrections](corrections-and-dead-ends.md#the-steps-0405-numbers-as-a-comparable-baseline);
> what remains here is the design and the one finding that does not depend on the target.

---

## The setup

The first predictor regresses per-cell **paclitaxel response** — the column
`obs["viability_paclitaxel"]`, i.e. the bulk per-line value broadcast to every cell of the matching
line. It is loaded by `ScGPTDrugDataset` (`scripts/model/dataset.py`, `target_drug="paclitaxel"`) and
trained with `train_multitask.py --use-rep {X_scGPT|X_pca} --drugs paclitaxel` (`output_dim = 1`).
Built **smallest-first** (plan §Prototyping) as a methodological probe before scaling out.

The split is `obs["split_paclitaxel"]`, its own dedicated split over paclitaxel-labelled lines,
written by `create_splits.py` `run()`: sklearn `train_test_split` over **whole cell lines**
(`random_state=42`, group = `Cell_line`, 70/15/15 as `test_size=0.30` then `0.50`).

> **Which single-task is this?** This dedicated-`split_paclitaxel` baseline, **not** the 8-run
> experiment matrix's single-task cells, which use the *shared* `split_ctrp` and live in
> [Step 05](05-multitask-results.md). The two are on different splits and are not comparable — an
> apples-to-apples *"does multi-task help paclitaxel?"* comparison still needs a single-task re-run on
> `split_ctrp`.

> **Scope — 1 database, 1 score.** One CTRPv2 metric, one compound: the narrowest slice of the
> project. The widening happens in [Step 05](05-multitask-results.md) (still CTRPv2, across drugs) and
> ultimately in [Step 06](06-planned-work.md#a-cross-database-integration).

---

## The data leak, and why the grouped split was forced (08.05.2026)

**This is the result the page exists for, and it survives the target change** — it follows from the
*structure* of the label, not from which score was loaded.

**Step 1 — random 70/15/15 split over cells (the deliberate mistake).** PCA's validation MSE came out
implausibly low — far better than scGPT's, and better than any honest generalization number should be.

**Why.** With cells split at random, the same cell line lands in both train and validation. The label
is **constant within a line** ([Step 03](03-model-and-training-design.md#every-cell-of-a-line-carries-the-identical-label)),
and PCA isolates each line as its own tissue "island" in the embedding, so the model reduces to a
nearest-neighbour lookup of a label it has already memorized. The validation score measured
memorization, not generalization; the implausibly low PCA number was the tell.

**Step 2 — cell-line-grouped split**, the correct design for this label structure: whole cell lines
are partitioned, so no line appears in two splits. PCA's validation MSE **collapsed**, confirming the
cheating; scGPT's barely moved.

**Step 3 — aggressive regularization** (hidden 256→64, dropout 0.3→0.5, weight decay 1e-5→1e-3).
Validation MSE ended up near-equal for the two representations, but **PCA's train/val gap was roughly
twice scGPT's** — the core hypothesis, and the plan's Fig. 4 prediction: PCA cheats by classifying cell
line, scGPT overfits far less.

✅ On-plan, and the random-split → leak-diagnosis → grouped-split arc is not in the plan but is the
"find failures cheaply, document even suboptimal versions" discipline it asks for. Worth keeping as a
result in its own right.

**What has to be re-measured at the sweep:** every MSE, every gap, and the split sizes — the overlap
moves from 180 to 181 cell lines and the labelled-cell count moves with it.
