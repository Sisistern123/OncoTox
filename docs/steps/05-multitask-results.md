# Step 05 — Multi-task masked loss (545 CTRPv2 drugs) & run versioning

*Part of [OncoTox project progress](../project_progress.md). Covers: the multi-task masked-loss
model over all 545 CTRPv2 drugs, its results vs. the per-drug-mean baseline, and the run-versioning
ledger that records every training run.*

This moves from plan-Phase-2 (single-task) into plan-Phase-3 (masked-loss multi-task). Masked-loss
mechanics are in [Step 03](03-model-and-training-design.md). These are the **multi-task (all-drugs,
K=545) rows of the 8-run experiment matrix**
([index](../project_progress.md#experiment-matrix--pca-vs-scgpt)).

> **Scope — still 1 database, 1 score; "multi-task" here = multi-*drug* only.** Every one of the
> K=545 heads predicts the **same** metric (`cpd_avg_pv` viability) from the **same** database
> (CTRPv2). This validates the masked-loss machinery on intra-CTRPv2 sparsity, but it is **not**
> the plan's ultimate multi-task goal, which is **cross-database** (CTRPv2 + PRISM + GDSC) **and
> multi-metric** (efficacy *and* toxicity). That integration — the real "combine all" — is
> [Step 06](06-cross-database-integration.md) and is **not yet started**. Do not read the 545-head
> run as "multi-task complete."

---

## Multi-task masked loss over all 545 CTRPv2 drugs (26.05.2026)

**New target artifacts written by `ctrp_to_h5ad.py`:**

- `obsm["Y_ctrp"]` — float32 (n_cells, K), per-cell viability, NaN where missing.
- `obsm["M_ctrp"]` — bool (n_cells, K), True where observed.
- `uns["ctrp_drugs"]` — ordered length-K drug-name list (column order of Y/M).
- `obs["split_ctrp"]` — **one drug-agnostic, cell-line-grouped 70/15/15 split** written by
  `create_splits.py` `run_multi()`, shared across all heads (leakage-free for every drug at once;
  a single shared split is only possible *because* the leakage control is at the cell-line level).
- Legacy flat `viability_<drug>` / `train_mask_<drug>` / `split_<drug>` kept for back-compat.

**Drug-scope filter:** keep a drug only if screened on ≥ `--min-cell-lines` overlapping
cell lines (default 50). This run used **`--all-drugs` (= min 0) → K = 545 drugs**.

**Run-time overlap reported by the pipeline:** **180 / 198** SCP542 cell lines overlap
CTRPv2 (the stricter-normalization 180; cf. the audit's 190 in
[Step 01](01-datasets-and-harmonization.md)).

**`split_ctrp` distribution (shared by all four runs below):**

- Cells: **train 34,126 / val 7,121 / test 5,980 / unassigned 6,286**
- Cell lines: **126 train / 27 val / 27 test**

**Model & training:** a single `OncoMLP` with `output_dim = K`, fed by `MultiDrugDataset`
(`scripts/model/dataset.py`) whose 3-tuple `(x, y, mask)` batches `train_model` auto-detects to
switch into **masked MSE** (mean over observed entries only). Up front,
`train_multitask._per_drug_constant_mse` computes a **per-drug-mean sanity baseline** — the proper
null model here: for each drug it predicts the constant train-set mean viability over that drug's
observed cells. Because labels cluster near 1.0, that constant is already a strong predictor, so a
head only counts as having *learned* response if it **beats its own drug's constant**.

**Shared hyperparameters** (from `config.json` / `run_meta.json`): batch 128, epochs 50
(early-stopped), lr 1e-3, weight_decay 1e-3, dropout 0.5, input_dropout 0.1, grad_clip 1.0,
scheduler patience 3, early-stop patience 10, seed 42, loss MSE, norm LayerNorm.
scGPT input_dim **512** / hidden (128,64); PCA input_dim per `X_pca` / hidden (64,32).

**The 8-run matrix (refreshed 13.06.2026; all share `split_ctrp`, n_train 34,126 / n_val 7,121).**
Baseline mean MSE over drugs is the per-drug constant: **0.0434** for K=1 paclitaxel, **0.0097** for
K=545. Run dirs `runs/20260613_1648xx–1651xx_*` (see `runs/runs_index.csv`).

| Gene set | Task | Rep | Best val MSE | Best ep | Model mean MSE | Heads beat baseline |
|---|---|---|---|---|---|---|
| `hvg5000` | single (K=1) | scGPT | 0.0406 | 3 | 0.0406 | 1 / 1 |
| `hvg5000` | single (K=1) | PCA | **0.0372** | 5 | 0.0372 | 1 / 1 |
| `hvg5000` | all (K=545) | scGPT | 0.0107 | 6 | 0.0106 | **135 / 545** |
| `hvg5000` | all (K=545) | PCA | 0.0110 | 8 | 0.0112 | 103 / 545 |
| `all_genes` | single (K=1) | scGPT | 0.0442 | 9 | 0.0442 | **0 / 1** |
| `all_genes` | single (K=1) | PCA | **0.0334** | 9 | 0.0334 | 1 / 1 |
| `all_genes` | all (K=545) | scGPT | 0.0105 | 7 | 0.0104 | **141 / 545** |
| `all_genes` | all (K=545) | PCA | 0.0114 | 10 | 0.0116 | 80 / 545 |

**Reading the results:**

- **All-drugs (K=545), heads-beating-baseline — the honest metric:** **scGPT beats PCA in both gene
  sets** — `hvg5000` **135 vs 103**, `all_genes` **141 vs 80**. scGPT's margin is *larger* on the
  full transcriptome. This is the core-hypothesis signal: scGPT learns real response on more drugs
  than the PCA baseline. (Absolute MSE ≈ 0.011 stays misleadingly low because the baseline is already
  0.0097 — read heads, not raw MSE.)
- **Single paclitaxel (K=1) flips it:** **PCA edges scGPT** in both gene sets (`hvg5000` 0.0372 <
  0.0406; `all_genes` 0.0334 < 0.0442), and `all_genes`/scGPT even fails to beat the constant
  baseline (**0/1**). On one low-variance drug a 50-d PCA is competitive; scGPT's advantage emerges
  **across many drugs**, not a single one.
- ⚠️ **Capacity caveat:** PCA is 50-d into hidden `(64,32)`, scGPT 512-d into `(128,64)` — the reps
  differ in dimensionality and head size, so the single-drug PCA edge may partly reflect capacity,
  not representation quality. Matching capacity is a [TODO](../TODO.md) item.

✅ On-plan: masked-loss multi-task, correctly gated behind a working single-task baseline,
with the cheap sanity baseline the plan's prototyping section calls for.

> ⚠️ **Key deviation — what "multi-task" means today:** the plan frames multi-task as
> **cross-database** (CTRPv2 + PRISM + GDSC heads). What's built is multi-task **across the
> 545 drugs of one database (CTRPv2)**. A legitimate *intermediate* step that validates the
> masked-loss machinery — but PRISM/GDSC are **not yet integrated**, so plan-Phase-3 is only
> half done. Don't read the 545-head run as "the multi-task goal is complete."

> ⚠️ **Split note — these are the matrix single-task cells, not Step 04's.** The K=1 rows above use
> `--drugs paclitaxel` on the **shared `split_ctrp`** (27 held-out lines), the same split as the
> K=545 runs — so within this table every comparison is apples-to-apples. They are **not** comparable
> to [Step 04](04-single-task-results.md)'s progression, which uses the separate `split_paclitaxel`
> (25 held-out lines). Different splits = different held-out cell lines.

---

## Run versioning (26.05.2026)

**Run versioning** (`training_utils.create_run_dir` / `save_run`): every
`train_multitask.py` run writes a self-contained `runs/<timestamp>_<tag>/`:

- `config.json` — exact `TrainConfig`.
- `run_meta.json` — scope, rep, dataset sizes, hidden_dims, host/python/torch info, drug list.
- `history.csv` — epoch, train_mse, val_mse, lr.
- `summary.json` — best_val_mse, best_epoch, baseline-vs-model mean MSE, heads-beating count.
- `best_model.pt` — best-val-MSE state_dict.
- `per_drug_results.csv` — drug, model_val_mse, baseline_val_mse, delta, n_val.

Plus one row per run in `runs/runs_index.csv` (columns: run_id, tag, scope, rep, K,
n_train_cells, n_val_cells, best_epoch, best_val_mse, baseline_mean_mse, model_mean_mse,
n_beats_baseline, n_total_heads, started_at, finished_at). `runs/` is gitignored.

✅ On-plan: satisfies "retain every working version + data to re-run + results, even
suboptimal ones."

The full 545-head run is reproduced with `train_multitask.py --use-rep {X_scGPT|X_pca}` (omitting
`--drugs` selects all K).
