# Step 06 — Cross-database integration (CTRPv2 + PRISM + GDSC, efficacy + toxicity)

> **Status: ❌ NOT STARTED — placeholder.** This is the true target of the project's
> multi-task goal and the step where "combine all" actually happens. Data is downloaded and
> harmonized ([Step 01](01-datasets-and-harmonization.md)); no integrated training has run yet.
> This file documents the intended scope so the whole project structure is visible.

*Part of [OncoTox project progress](../project_progress.md). Covers (future): extending the
masked-loss multi-task model from one database/one metric to **multiple databases and multiple
response metrics**.*

---

## Why this is a distinct step (and why 04/05 are not it)

The work so far is narrow on **two axes** that this step widens:

| Axis | Steps 04–05 (done) | Step 06 (this, planned) |
|---|---|---|
| **Database** | CTRPv2 only | CTRPv2 **+ PRISM + GDSC** |
| **Response metric** | `cpd_avg_pv` viability only | viability **+ LN_IC50 / AUC + toxicity** (efficacy *and* toxicity) |
| **What "multi-task" means** | multi-**drug** (545 heads, one metric) | multi-**database / multi-metric** heads |

So the 545-head run in [Step 05](05-multitask-results.md) is multi-task only across drugs of a
single source. Step 06 is the plan's actual **Goal Option B** — "simultaneously predict multiple
response metrics (efficacy and toxicity) across pan-cancer cell lines using a multi-task learning
setup (via masked losses) to handle sparse or missing labels."

## Planned approach (from the project plan)

1. **Add PRISM first** (the plan's explicit next move): the much larger, far sparser dataset
   (915 lines, 6,575 compounds, ~29 % non-null). This needs a PRISM/GDSC analog of
   `scripts/preprocessing/ctrp_to_h5ad.py` that emits additional `Y_*` / `M_*` blocks (or extends
   `Y_ctrp`/`M_ctrp` into a unified matrix), using the harmonized drug catalog
   ([Step 01](01-datasets-and-harmonization.md): name + BRD-ID + DrugBank links, already in
   `data/drug/`) as the join keys to align PRISM compounds onto existing CTRPv2 heads and to add
   PRISM-only heads. `scripts/preprocessing/layout.py` gains the new source files.
2. **Then GDSC** (`LN_IC50` / `AUC`) — a *different metric type*, so this is where heads stop being
   homogeneous and the model becomes genuinely multi-metric.
3. **Masked loss across the union:** generalize the `MultiDrugDataset` mask machinery
   (`scripts/model/dataset.py`) and the masked loss ([Step 03](03-model-and-training-design.md)) to
   a block-sparse label matrix spanning all sources — each (cell line × drug × metric) entry is
   observed in only some databases — with per-database/per-metric head grouping and weighting in
   `scripts/training/train_multitask.py`.
4. **Cross-database splits:** keep the cell-line-grouped, leakage-free split discipline of
   `create_splits.py` ([Step 04](04-single-task-results.md)) across the unified cell-line set.

## Open design questions (to resolve when starting)

- One shared trunk with per-database/per-metric head groups, vs. per-metric normalization before a
  shared head? (Metrics live on different scales: viability ≈ [0,1], `LN_IC50` unbounded.)
- How to weight databases/metrics in the loss given wildly different coverage (CTRPv2 21 % vs PRISM
  29 % vs GDSC 3 % non-null)? Ties into the [Step 05](05-multitask-results.md) open question on
  per-head / uncertainty weighting.
- Whether to harmonize metrics onto a common response scale, or keep them as separate heads.
- Reconcile the 190-vs-180 cell-line normalization ([Step 01](01-datasets-and-harmonization.md))
  before unioning sources.

## Definition of done

- A single model with masked heads spanning ≥2 databases trains leakage-free.
- Per-database / per-metric "heads-beating-baseline" reported (same honest metric as
  [Step 05](05-multitask-results.md)).
- Result documented here with run IDs in the versioning ledger
  ([Step 05](05-multitask-results.md)).
