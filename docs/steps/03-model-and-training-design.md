# Step 03 — Model & training design

*Part of [OncoTox project progress](../project_progress.md). Covers: what a single training
example is (input, output, target, mask), the supervised learning paradigm and the weak-supervision
it rests on, the masked-loss formulation, the exact MSE definitions, and the model architecture —
each tied to the code that implements it (`scripts/model/`, `scripts/training/`).*

---

## The learning problem — weakly-supervised, fully-supervised regression

The downstream task is a **continuous regression**: map one cell's transcriptomic representation to
a drug-response scalar. It is **fully supervised** — `train_model` in
`scripts/training/training_utils.py` optimizes a (masked) **MSE** or **Huber** loss directly against
observed labels, with no classification, pseudo-labeling, or consistency/contrastive objective.

The *weak supervision* is in the **labels, not the algorithm**: the response value is a **bulk**
(cell-line-level) measurement broadcast onto every single cell of that line (see the target section
below). Mapping one average bulk score onto thousands of heterogeneous single cells deliberately
injects label noise — this is the project's central modeling assumption (plan §Strategy): by forcing
the network to map a heterogeneous single-cell input to the *average* bulk response, it must learn
the transcriptomic signatures of sensitivity rather than any per-cell idiosyncrasy.

Two points that are easy to misstate:

- **It is not semi-supervised.** Cells lacking a label for a given (cell line × drug) are **dropped
  from the loss** by the mask `M` (below), not exploited as unlabeled data. The mask handles label
  **sparsity**; it does not add an unsupervised objective.
- **The only self-supervised component is upstream and frozen.** scGPT is a foundation model
  pretrained self-supervised on ~33 M human cells; here it is a **fixed feature extractor** — the
  512-dim cell embedding is read from `obsm["X_scGPT"]` and never fine-tuned. The pipeline is
  therefore *"self-supervised representation (frozen) → supervised regression head."* The PCA
  baseline `X_pca` is an **unsupervised** linear transform feeding the identical head, and exists
  only to test the core hypothesis (that scGPT's denoised manifold overfits less than PCA's
  tissue-clustered one — see [Step 04](04-single-task-results.md)).

---

## A training example = one single cell

The unit fed to the network is **one cell**, never a cell line / cancer type / drug aggregate.
Cell-line, cancer-type and drug identity are **not** input features — they only determine the label
and the train/val/test grouping.

- **Input `X` (per cell):** a single embedding vector, shape `(D,)`, selected with `--use-rep` and
  read from `adata.obsm[use_rep]` by the dataset classes in `scripts/model/dataset.py`:
  `X_scGPT` → **512-dim** scGPT embedding, or `X_pca` → the (HVG-5000 / all-genes) PCA baseline
  (≈50-dim). The genes themselves are never seen by the network — only these representations.
- **Output (per cell):** the final layer is a single `Linear(prev_dim → output_dim)`
  (`OncoMLP.py`). For **single-task** `output_dim = 1` (one drug's viability); for **multi-task**
  `output_dim = K`, so the "`K` drug heads" are literally the **K rows of that one output matrix**
  over a shared trunk — there are no separate per-drug sub-networks. The default catalog is
  **K = 545** CTRPv2 drugs.

---

## Target `y` — what "viability" is and at what resolution it is defined

The label is CTRPv2 **`cpd_avg_pv`** (compound average percent viability): the dose-averaged
fraction of cells surviving relative to vehicle controls, so **most values sit near 1.0** (a fully
resistant line ≈ 1, a sensitive one < 1). It is a **bulk, per-(cell line × drug)** quantity — *not*
measured per single cell.

`scripts/preprocessing/ctrp_to_h5ad.py` turns this into per-cell labels in three steps, all visible
in the code:

1. **Aggregate** every CTRPv2 measurement to one value per (cell line, drug) by mean —
   `groupby(["ccl_name_norm","cpd_name_norm"]).mean()` in `_build_drug_table`.
2. **Pivot** to a (cell line × drug) matrix, column order pinned to `uns["ctrp_drugs"]`.
3. **Broadcast** each bulk value to every cell of the matching line —
   `Y_full = cl_drug_matrix.reindex(cell_line_norm.values)`.

⇒ **every cell of a given line carries the identical label vector** — the structural fact that makes
grouped splitting mandatory (below) and absolute MSE small (next). The result is stored as
`obsm["Y_ctrp"]` `(n_cells, K)` float32 (NaN where unscreened) with the length-K column→drug map in
`uns["ctrp_drugs"]`. A drug column is kept only if screened on ≥ `--min-cell-lines` overlapping
lines (default 50; the K=545 run used `--all-drugs`, i.e. min 0). Cancer type is never a label or a
feature — it only colors the UMAPs in [Step 02](02-preprocessing-and-embeddings.md).

---

## Mask `M` — the sparsity-handling mechanism (plan sub-goal 2)

Most (cell line × drug) pairs were never assayed, so the label matrix is sparse. `ctrp_to_h5ad.py`
records `obsm["M_ctrp"]` `(n_cells, K)` bool, **True iff that cell's line was actually screened
against that drug**. `MultiDrugDataset` (`scripts/model/dataset.py`) fills missing `Y` with 0.0 only
so the tensor is finite, then carries `M` alongside so the loss can ignore those zeros.

This is the masked-loss core that the plan's sub-goal 2 calls for, and it is what will generalize to
the cross-database block-sparse matrix in [Step 06](06-cross-database-integration.md).

---

## Exact loss & MSE definitions (`training_utils.py`)

`train_model` auto-detects multi-task batches by peeking for a 3-tuple `(x, y, mask)`
(`_is_multitask_loader`) and switches loss accordingly (`_make_loss_fn`):

- **Per-element error** `sq = (pred − y)²` (MSE), or `smooth_l1_loss(beta=0.05)` for `--loss huber`
  (robust to the occasional outlier viability while staying quadratic near 0).
- **Masked batch loss** = `(sq · M).sum() / M.sum()` (`_masked_mean`, denominator clamped ≥ 1) —
  the mean over **only observed (cell, drug) entries**. Gradients therefore flow *only* through
  observed entries; unscreened pairs contribute nothing to loss, gradient, or metric.
- **Epoch MSE** accumulates `Σ(sq·M)` and `ΣM` across batches, then divides — i.e. it is
  **entry-pooled**: every observed (cell, drug) pair is weighted equally, so high-coverage drugs and
  cell lines with more cells count proportionally more. This entry-pooled val MSE is `best_val_mse`
  in `summary.json` and `val_mse` in `history.csv`, and is what early-stopping/scheduler watch.
- **Single-task** uses plain `((pred − y)²).mean()` (no mask).

**Two aggregations are reported — do not conflate them:**

| Name | Where | Definition | K=545 scGPT |
|---|---|---|---|
| **Entry-pooled MSE** | `best_val_mse`, `history.csv` | `Σ sq·M / Σ M` over all observed entries | **0.0105** |
| **Macro per-drug MSE** | `model_mean_mse` / `baseline_mean_mse` | per-drug `Σ_cells sq_k / n_k`, then `np.nanmean` over drugs (equal weight per drug) | **0.0103** |

Only the **macro per-drug** numbers feed the **per-drug-mean baseline** comparison
(`train_multitask._per_drug_constant_mse`): a null model predicting, for each drug, the constant
train-set mean viability over its observed cells. "**Heads beating baseline**" counts drugs whose
model per-drug val MSE < that constant's (scGPT 142/545, PCA 97/545 — [Step 05](05-multitask-results.md)).
**This is the honest metric**, because absolute MSE ≈ 0.01 is misleadingly tiny: with `cpd_avg_pv`
clustered near 1.0 the constant baseline already achieves ≈ 0.0097, so beating it — not the raw MSE —
is what signals a head actually learned response.

---

## Why splits must be cell-line-grouped

Because the label is **constant within a cell line**, a random cell-level split would place cells of
the same line (hence the same label and near-identical tissue signature) in both train and val,
letting the model memorize the per-line label instead of learning response. `create_splits.py`
therefore partitions **whole cell lines** 70/15/15 (`split_ctrp`, drug-agnostic so it is leakage-free
for all heads at once; `split_paclitaxel` for the single-task case). This is not a detail — it is the
methodological control that exposes the PCA-vs-scGPT overfitting gap in
[Step 04](04-single-task-results.md).

---

## Model architecture & regularization (`OncoMLP.py`, 25.05.2026)

`OncoMLP` is a deliberately **small** MLP (the plan asks for the smallest functional model):
input dropout → `[Linear → LayerNorm → GELU → Dropout]` per hidden layer → `Linear(→ output_dim)`.
Defaults encode specific choices for this regime:

- **`hidden_dims` = (64,32) for PCA, (128,64) for scGPT** (`train_multitask.DEFAULT_HIDDEN_DIMS`) —
  scGPT's 512-d input warrants a slightly wider trunk than PCA's ≈50-d.
- **LayerNorm** (not BatchNorm): batches are cell-line-grouped and small, so BatchNorm running
  statistics are noisy; LayerNorm normalizes per-sample and is stable here.
- **GELU** rather than ReLU — a smoother activation for continuous-valued targets.
- **Heavy regularization** — `input_dropout=0.1` on the raw embedding plus `dropout=0.5` in the
  trunk, and Adam **weight decay 1e-3** (L2) — all aimed at the same failure mode: suppressing
  cell-line memorization given the broadcast labels.

Training (`training_utils.train_model`, all in `TrainConfig`) is seeded (42) and uses **Adam**
(lr 1e-3), **ReduceLROnPlateau** (factor 0.5, patience 3) on the val MSE, **gradient clipping**
(max-norm 1.0), and **early stopping** (patience 10); the **best-val-MSE checkpoint is restored** at
the end rather than the last-epoch weights. The single entrypoint `train_multitask.py` exposes these
as flags (`--use-rep`, `--drugs`, `--batch-size 128`, `--epochs 50`, `--lr`, `--weight-decay`,
`--dropout`, `--input-dropout`, `--loss {mse,huber}`, `--hidden-dims`, `--seed`); run artifacts are
written by `create_run_dir`/`save_run` ([Step 05](05-multitask-results.md)).
