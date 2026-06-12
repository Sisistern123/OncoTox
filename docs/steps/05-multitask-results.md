# Step 05 — Multi-task masked loss (545 CTRPv2 drugs) & run versioning

*Part of [OncoTox project progress](../project_progress.md). Covers: the multi-task masked-loss
model over all 545 CTRPv2 drugs, its results vs. the per-drug-mean baseline, and the run-versioning
ledger that records every training run.*

This moves from plan-Phase-2 (single-task) into plan-Phase-3 (masked-loss multi-task). Masked-loss
mechanics are in [Step 03](03-model-and-training-design.md).

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

**The four runs (all share `split_ctrp`; n_train 34,126 / n_val 7,121):**

| Run id | Rep | K | Best epoch | Best val MSE | Baseline mean MSE | Model mean MSE | Heads beat baseline |
|---|---|---|---|---|---|---|---|
| `20260526_132914_multitask_X_scGPT_subset_K1` | X_scGPT | 1 (paclitaxel) | 11 | 0.0412 | 0.0434 | 0.0412 | 1 / 1 |
| `20260526_132952_multitask_X_pca_subset_K1` | X_pca | 1 (paclitaxel) | 5 | 0.0393 | 0.0434 | 0.0393 | 1 / 1 |
| `20260526_133012_multitask_X_scGPT_all_drugs` | X_scGPT | 545 | 7 | 0.0105 | 0.0097 | 0.0103 | **142 / 545** |
| `20260526_133112_multitask_X_pca_all_drugs` | X_pca | 545 | 6 | 0.0112 | 0.0097 | 0.0114 | 97 / 545 |

**Reading the results:**

- The K=545 ~0.0105 looks good **only because most viability values sit near 1.0**, so the
  per-drug-mean baseline is already 0.0097. The honest metric is **heads-beating-baseline**:
  **scGPT 142/545 vs PCA 97/545** — scGPT wins on ~47 % more heads at the same K and split.
- Worst heads (model < baseline) are the lowest-coverage ones (n_val = 221): `brd-k30748066`,
  `vx-680`, `brd-k33514849`, `brd9876:mk-1775 (4:1 mol/mol)`, `bafilomycin a1` — candidates
  to drop or down-weight.
- Largest single win in both reps: `gsk-j4` (model ≈ 0.000 vs baseline 0.011, n = 221) —
  sanity check that a head can fit a low-variance drug-line combination.

✅ On-plan: masked-loss multi-task, correctly gated behind a working single-task baseline,
with the cheap sanity baseline the plan's prototyping section calls for.

> ⚠️ **Key deviation — what "multi-task" means today:** the plan frames multi-task as
> **cross-database** (CTRPv2 + PRISM + GDSC heads). What's built is multi-task **across the
> 545 drugs of one database (CTRPv2)**. A legitimate *intermediate* step that validates the
> masked-loss machinery — but PRISM/GDSC are **not yet integrated**, so plan-Phase-3 is only
> half done. Don't read the 545-head run as "the multi-task goal is complete."

> ⚠️ **Not comparable:** the K=1 paclitaxel numbers here (val 0.0412 scGPT / 0.0393 PCA on
> `split_ctrp`, 27 held-out lines) are **not** comparable to the single-task numbers in
> [Step 04](04-single-task-results.md) (0.0336 / 0.0351 on `split_paclitaxel`, different
> held-out lines). An apples-to-apples "does multi-task help paclitaxel?" comparison still
> needs a single-task re-run on `split_ctrp`.

> ⚠️ **Provenance:** these four `run_meta.json` files record the targets h5ad at the **old
> flat path** `data/scRNAseq_SCP542/metadata/…_with_targets.h5ad`, i.e. they predate the
> variant-based `processed/<variant>/` layout refactor (commit `900abe6`). Re-running under
> `processed/hvg5000/` should reproduce them but hasn't been done.

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
