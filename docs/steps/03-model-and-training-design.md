# Step 03 — Model & training design

*Part of [OncoTox project progress](../project_progress.md). Covers: exactly what a training
example is (input, output, target, mask), how the reported MSE is computed, and the high-level
design choices — target, masking strategy, and the supervised training paradigm.*

Everything below is read straight out of the code (`scripts/model/`, `scripts/training/`,
`scripts/preprocessing/`).

---

## Model mechanics — exact input, output, targets, and MSE

*The single most-asked detail: what does a training example actually look like, what is
the label defined over, and how is the reported MSE computed?*

### Unit of a training example = **one single cell**

The model is fed **per single cell**, never per cell line / per cancer type / per drug. Cell
line, cancer type, and drug identity are *not* input features — they only determine the label
and the train/val/test grouping.

- **Input `X` (per cell):** one embedding vector, shape `(D,)`.
  - `X_scGPT` → **512-dim** scGPT embedding, or `X_pca` → PCA of the (HVG-5000 / all-genes)
    expression. Selected via `--use-rep`; read from `adata.obsm[use_rep]` (`dataset.py`).
- **Output (per cell):**
  - **Single-task:** 1 scalar = predicted viability for the one drug.
  - **Multi-task:** **`K`-dim vector**, one scalar per CTRPv2 drug "head". The net is one shared
    trunk + a single final `Linear(prev_dim → K)` — so the "K heads" are the K **rows of that
    one output layer**; all hidden layers are shared across drugs (`OncoMLP.py:56`).
    The default catalog is **K = 545** drugs.
- Architecture: input-dropout → `[Linear → LayerNorm → GELU → Dropout]×hidden_dims` → `Linear→K`.
  hidden_dims **(64,32) for PCA, (128,64) for scGPT** (`train_multitask.py:DEFAULT_HIDDEN_DIMS`).

### Target `y` — what "viability" is, and at what granularity it is defined

- Viability = CTRPv2 **`cpd_avg_pv`** (Compound Average Percent Viability — weighted average of
  surviving cells across all tested doses; most values sit **near 1.0**).
- It is defined **per (cell line × drug)** — a *bulk* number, **not** measured per single cell.
  `ctrp_to_h5ad.py` aggregates every CTRPv2 measurement to one value per (cell line, drug) by
  **mean** (`_build_drug_table`: `groupby(["ccl_name_norm","cpd_name_norm"]).mean()`), pivots to
  a (cell line × drug) matrix, then **broadcasts each bulk value to every single cell of the
  matching cell line** (`Y_full = cl_drug_matrix.reindex(cell_line_norm.values)`).
  ⇒ **All cells of a given cell line carry the identical label vector.**
- Stored as `obsm["Y_ctrp"]` `(n_cells, K)` float32 (NaN where missing) + `obsm["M_ctrp"]`
  `(n_cells, K)` bool mask; `uns["ctrp_drugs"]` is the length-K column→drug ordering.
- **Not** per cancer type. Cancer type is used only to color the UMAPs (see
  [Step 02](02-preprocessing-and-embeddings.md)), never as label/feature.
- Drug column kept only if screened on ≥ `--min-cell-lines` overlapping cell lines (default 50;
  the K=545 run used `--all-drugs` = min 0).

### Mask `M` — handles the sparsity

`M[cell, k] = True` iff that cell's **cell line** was actually screened against drug `k`. Missing
(cell-line-never-tested-on-drug-k) entries are **excluded from the loss and from all metrics** —
this is the masked-loss machinery (`MultiDrugDataset` fills missing `Y` with 0.0 only so it's
safe to pass through PyTorch; the mask zeroes them out before they touch the loss).

### Exact MSE computation (`training_utils.py`)

**Multi-task (masked) — per batch:**
1. Per-element squared error `sq = (preds − y)²`, shape `(batch, K)`.
2. Mask it: `sq * M` zeroes unobserved (cell, drug) pairs.
3. Batch loss `= (sq*M).sum() / M.sum()` (`_masked_mean`) — mean over **only observed
   (cell, drug) entries**; returns 0 if a batch has no observed entries.

**Epoch MSE (train and val)** accumulate across batches:
`running_sq_sum += (sq*M).sum()`, `running_n += M.sum()`, then `MSE = running_sq_sum / running_n`.
This is **entry-pooled**: every observed *(cell, drug)* pair is weighted equally, so high-coverage
drugs — and cell lines contributing more cells — count proportionally more. It is **not** a
per-drug average. `best_val_mse` in `summary.json` / `history.csv` `val_mse` is this number.

**Single-task:** plain `((preds − y)²).mean()` over the batch (no mask); epoch = `sq.sum()/numel`.

**Gradients** flow only through observed entries (the mask zeroes the rest). Optimizer: **Adam**
(lr 1e-3, weight_decay 1e-3) on masked **MSE** (or masked-Huber, `beta=0.05`, via `--loss huber`);
grad-clip max-norm 1.0; ReduceLROnPlateau (factor 0.5, patience 3); early-stop (patience 10) on
the entry-pooled val MSE; best-val-MSE checkpoint is restored at the end.

### ⚠️ Two different aggregations are reported — do not confuse them

| Name | Where | Definition | K=545 scGPT |
|---|---|---|---|
| **Entry-pooled MSE** | `best_val_mse`, `history.csv` | `Σ sq·M / Σ M` over all observed entries | **0.0105** |
| **Macro per-drug MSE** | `model_mean_mse` / `baseline_mean_mse` | per-drug `Σ_cells sq_k / n_k`, then `np.nanmean` over drugs (each drug equal weight) | **0.0103** |

Only the **macro per-drug** numbers are what the **per-drug-mean baseline** is compared against
(`train_multitask._per_drug_constant_mse`): the baseline predicts, for each drug, the constant
**train-set mean viability over observed cells**; "**heads beating baseline**" counts drugs whose
model per-drug val MSE < that constant's per-drug val MSE (scGPT 142/545, PCA 97/545 — see
[Step 05](05-multitask-results.md)).

### Why this matters for reading every MSE in this doc

- MSE ≈ 0.01 looks tiny **because `cpd_avg_pv` clusters near 1.0** — the constant baseline already
  sits at 0.0097, so absolute MSE is nearly meaningless; **heads-beating-baseline** is the honest
  signal.
- Because the label is constant within a cell line, splits **must** be grouped by `Cell_line`
  (70/15/15, `create_splits.py`) — otherwise the model memorizes the per-cell-line label
  (this is exactly the leakage story in [Step 04](04-single-task-results.md)).

---

## Design choices — target · masking · training paradigm (reference)

*Consolidated answers to the recurring "how did you set the problem up?" questions, so the
thesis writeup can cite design choices in one place.*

### Training paradigm — **supervised regression** (with a self-supervised feature prior)

- The downstream predictor is trained **fully supervised**: continuous-valued regression of
  observed viability labels, optimized with (masked) **MSE / Huber** (`training_utils.py`).
  There is **no** classification, no pseudo-labeling, no consistency/contrastive objective on
  the downstream task.
- **Not semi-supervised.** Cells without a label for a given (cell line × drug) are simply
  **excluded from the loss** via the mask — they are *dropped*, not used as unlabeled data. The
  masked loss handles label **sparsity**, it does **not** make the method semi-supervised.
- The **only** self-/unsupervised component is **upstream and frozen**: scGPT is a pretrained
  (self-supervised) foundation model used purely as a **fixed feature extractor** — its 512-dim
  embedding is read from `obsm` and never fine-tuned here. So the pipeline is
  **"self-supervised representation (frozen) → supervised regression head."**
- PCA (the `X_pca` baseline rep) is likewise an **unsupervised** feature transform feeding the
  same supervised head — it exists only to test the core hypothesis (scGPT overfits less).

### Target — what is regressed

- **CTRPv2 `cpd_avg_pv`** (Compound Average Percent Viability), continuous, mostly near 1.0.
- Defined **per (cell line × drug)** as a bulk number, then **broadcast to every single cell**
  of that cell line ⇒ all cells of a line share one label vector (full mechanics above).
- **Granularity that is *not* used:** cancer type (UMAP coloring only), and the label is never
  per-single-cell despite the per-cell input.

### Masking strategy — how label sparsity is handled

- `M[cell, k] = True` iff that cell's **cell line was actually screened against drug k**;
  stored as `obsm["M_ctrp"]` `(n_cells, K)` bool alongside `obsm["Y_ctrp"]` (NaN→0.0 filled).
- **Loss:** per-element `sq = (pred−y)²` is multiplied by `M`; batch loss `= (sq·M).sum()/M.sum()`
  (`_masked_mean`), so **unobserved (cell, drug) pairs contribute nothing to loss, gradients, or
  metrics**. Epoch MSE is **entry-pooled** (every observed pair weighted equally) — see the
  entry-pooled vs. macro-per-drug distinction above.
- **Splitting interacts with masking:** because the label is constant within a cell line, the
  split is **cell-line-grouped** (`split_ctrp`, drug-agnostic, 70/15/15) so it is leakage-free for
  every head at once (leak story in [Step 04](04-single-task-results.md), distribution in
  [Step 05](05-multitask-results.md)).
- This is the plan's **sub-goal 2 (masked-loss sparsity handling)** — done within CTRPv2.

### Model & training upgrade (25.05.2026)

- **Model** (`scripts/model/OncoMLP.py`): default **LayerNorm + GELU**, input dropout 0.1,
  configurable `hidden_dims` — **(64,32) for PCA, (128,64) for scGPT**.
- **Training** (`scripts/training/training_utils.py`): seeded (seed 42), Adam,
  **ReduceLROnPlateau** (factor 0.5, patience 3), **gradient clipping** (max-norm 1.0),
  **early stopping** (patience 10), best-val checkpoint restore, masked-loss support.

---

## Code & key variables

**Model — `scripts/model/OncoMLP.py`:**
`OncoMLP(input_dim, hidden_dims=(64,32), dropout_rate=0.5, input_dropout=0.1, norm="layer",
output_dim=1)`. `output_dim=K` makes the multi-task heads; `norm ∈ {layer, batch, none}`.

**Datasets — `scripts/model/dataset.py`:**
- `MultiDrugDataset` (multi-task) reads, by default, `use_rep`, `split_col="split_ctrp"`,
  `y_obsm_key="Y_ctrp"`, `mask_obsm_key="M_ctrp"`, `drugs_uns_key="ctrp_drugs"`; returns
  `(x, y, mask)` 3-tuples that `train_model` auto-detects as masked.
- `ScGPTDrugDataset` (single-task) reads `use_rep`, `target_drug` (→ `viability_<drug>` /
  `split_<drug>` columns); returns `(x, y)` pairs.

**Training utils — `scripts/training/training_utils.py`:** `train_model` (the loop), `_masked_mean`
(the masked MSE), plus run-versioning `create_run_dir` / `save_run` (see
[Step 05](05-multitask-results.md)).

**Entrypoint — `scripts/training/train_multitask.py`:** `DEFAULT_HIDDEN_DIMS` = `{X_pca:(64,32),
X_scGPT:(128,64)}`. CLI flags: `--use-rep {X_scGPT,X_pca}`, `--drugs <names…>` (omit = all K),
`--path`, `--batch-size 128`, `--epochs 50`, `--lr 1e-3`, `--weight-decay 1e-3`, `--dropout 0.5`,
`--input-dropout 0.1`, `--loss {mse,huber}`, `--hidden-dims …`, `--seed 42`, `--tag`,
`--baseline-topk 5`.

**Key variables consumed:** `obsm["X_scGPT"]` / `obsm["X_pca"]` (input `X`), `obsm["Y_ctrp"]`
(target), `obsm["M_ctrp"]` (mask), `uns["ctrp_drugs"]` (column→drug order), `obs["split_ctrp"]` /
`obs["split_paclitaxel"]` (grouping). All written in [Step 02](02-preprocessing-and-embeddings.md).
