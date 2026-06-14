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
CTRPv2 (180 = lines with actual post-QC measurements; the audit's 190 counts roster name-matches — see
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
**Matched trunk (14.06.2026):** both reps now use the **same** hidden layers `(128,64)`, so only the
input representation (scGPT 512-d / PCA ~50-d, and its first projection) differs — a fair
PCA-vs-scGPT comparison. (Earlier runs used `(64,32)` for PCA vs `(128,64)` for scGPT, which
handicapped PCA; see the capacity note below.)

**The 8-run matrix (matched trunk, 14.06.2026; all share `split_ctrp`, n_train 34,126 / n_val 7,121).**
Per-drug-mean baseline: **0.0434** (K=1 paclitaxel), **0.0097** (K=545). Run dirs
`runs/20260614_2056xx_*` (see `runs/runs_index.csv`).

**Single-task (K=1 paclitaxel) — the overfitting story:**

| Gene set | Rep | Train MSE | Val MSE | Train/val gap |
|---|---|---|---|---|
| `hvg5000` | scGPT | 0.023 | 0.041 | **0.018** |
| `hvg5000` | PCA | 0.010 | 0.039 | 0.028 |
| `all_genes` | scGPT | 0.033 | 0.045 | **0.012** |
| `all_genes` | PCA | 0.010 | 0.038 | 0.029 |

**All-drugs (K=545) — heads beating baseline:**

| Gene set | Rep | Val MSE | Heads beat baseline |
|---|---|---|---|
| `hvg5000` | scGPT | 0.0105 | 158 / 545 |
| `hvg5000` | PCA | 0.0105 | 156 / 545 |
| `all_genes` | scGPT | 0.0106 | 137 / 545 |
| `all_genes` | PCA | 0.0104 | **196 / 545** |

**Reading the results (matched trunk):**

- **Core hypothesis — supported (single-task):** scGPT **overfits far less** — train/val gap
  **0.012–0.018** vs PCA **0.028–0.029** — even though PCA's *raw* val MSE is slightly lower. scGPT
  trades a little fit for much better generalization, exactly the denoised-prior claim.
- **All-drugs, with capacity matched:** PCA is now **competitive/better** — `hvg5000` 158 vs 156
  (≈ tie), `all_genes` **PCA 196 vs scGPT 137**. The earlier scGPT advantage (135 vs 103; 141 vs 80)
  was **largely a capacity artifact**: PCA had been handicapped by the smaller `(64,32)` trunk.
- **Net:** scGPT's clear, robust win is **lower overfitting**, not higher absolute accuracy. Once PCA
  isn't handicapped, the two are close on raw predictive metrics (PCA even ahead on `all_genes`).
- **Which heads are even learnable** is driven by coverage + response variance — see
  `notebooks/drug_coverage.ipynb`: the ≈16-line drugs (n_val 221) are the unreliable/hardest heads,
  while high-coverage high-variance drugs (docetaxel, gemcitabine, oligomycin a) are the easiest.

> ⚠️ **Remaining caveat:** input dimensionality still differs (PCA ~50 vs scGPT 512), so the *first*
> projection isn't matched — for a fully controlled test, raise PCA `n_comps` toward 512
> ([TODO](../TODO.md)). Also note the all-drugs "gap" is not a clean overfit measure (train MSE is
> logged with dropout active, so it can exceed the masked val MSE).

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
