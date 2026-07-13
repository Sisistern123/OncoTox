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

The *weak supervision* lives in the **labels, not the algorithm**. The response value is a **bulk**
(cell-line-level) measurement, broadcast onto every single cell of that line (target section below).
Mapping one average score onto thousands of heterogeneous cells deliberately injects label noise.
That is the project's central assumption (plan §Strategy): forced to map a noisy single-cell input
to the *average* bulk response, the network must learn the transcriptomic signature of sensitivity
rather than any per-cell quirk.

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
  (**512-dim**, matched to the scGPT width via `add_pca.DEFAULT_N_COMPS`). The genes themselves are
  never seen by the network — only these representations.
- **Output (per cell):** the final layer is a single `Linear(prev_dim → output_dim)`
  (`OncoMLP.py`). For **single-task** `output_dim = 1` (one drug's response score); for **multi-task**
  `output_dim = K`, so the "`K` drug heads" are literally the **K rows of that one output matrix**
  over a shared trunk — there are no separate per-drug sub-networks. The default catalog is
  **K = 545** CTRPv2 drugs.

---

## Target `y` — the response score, and at what resolution it is defined

**This section is the canonical definition of the target; the other steps refer back to it.**

The label is a CTRPv2 **drug-response score**, always a **bulk, per-(cell line × drug)** quantity —
*not* measured per single cell. Which score is selected with `--score`
(`layout.CTRP_SCORES`), and every score is written to its **own** targets h5ad, so two scores can be
trained head-to-head without rebuilding anything. In **all** of them a *higher* value = a *more
resistant* line.

| `--score` | Definition | Source table |
|---|---|---|
| **`auc_z`** (default, 13.07.2026) | per-drug **z-score** of `auc`, over the 180 overlapping cell lines | ⟵ derived |
| `auc` | `area_under_curve / conc_pts_fit` — the sigmoid-fit AUC, **normalized** by the size of the concentration grid | `v20.data.curves_post_qc.txt` |
| `mean_pv` | *legacy:* unweighted mean of `cpd_avg_pv` (percent viability) over the dose grid; clusters near 1.0 | `v20.data.per_cpd_post_qc.txt` |

**Decision (13.07.2026): train on `auc_z`.** Three reasons, each fixing a defect of `mean_pv`:

1. **Use the curve fit, not the raw dose grid.** `mean_pv` averaged the measured viability points, so
   it was dominated by however many concentrations happened to sit in the flat, pre-response part of
   the curve, and it inherited every noisy well. `area_under_curve` comes from CTRP's post-QC sigmoid
   fit — the standard target in the field (Rees et al. 2016; DepMap/PharmacoGx), and the same
   quantity GDSC2 reports, which is what makes the [Step 06](06-cross-database-integration.md) merge
   possible at all.
2. **Normalize by the concentration grid.** CTRP's `area_under_curve` is an *integrated* area, so it
   grows with `conc_pts_fit` (8–29 points, usually 16) and is not comparable across compounds. `auc`
   divides it out (`_load_score_values`).
3. **Z-score within each drug** (`_zscore_per_drug`). Standardizing per drug removes the potency
   offset the model cannot infer from expression anyway, and puts every one of the K heads on the
   **same scale**, so no single well-covered drug dominates the shared masked loss. The per-drug mean
   and std are kept in `uns["ctrp_score_center"]` / `["ctrp_score_scale"]`, so predictions can be
   mapped back to AUC units.

> The two scores are **not** interchangeable. Globally they correlate at ρ ≈ 0.97, but that is
> inflated by between-drug potency differences. *Within* a drug, across cell lines — the only
> variation the model must actually predict — the median Spearman between `auc_z` and `mean_pv` is
> only **0.72** (min 0.42). Steps 04–05 were trained on `mean_pv` and are **not** directly comparable
> to `auc_z` numbers.
>
> **Known caveat:** the z-score statistics are computed over *all* overlapping lines, including those
> that later land in val/test. This is standard practice, but it is technically mild leakage (the
> held-out labels' mean/spread inform the normalization). Making it train-only requires computing the
> splits before the targets step.

`scripts/preprocessing/ctrp_to_h5ad.py` turns the chosen score into per-cell labels:

1. **Aggregate** to one value per (cell line, drug), averaging replicate experiments —
   `groupby(["ccl_name_norm","cpd_name_norm"]).mean()` in `_build_drug_table`.
2. **Standardize** per drug if `--score auc_z` (statistics computed per *cell line*, so a line with
   many cells does not pull its own mean).
3. **Pivot** to a (cell line × drug) matrix, column order pinned to `uns["ctrp_drugs"]`.
4. **Broadcast** each bulk value to every cell of the matching line —
   `Y_full = cl_drug_matrix.reindex(cell_line_norm.values)`.

⇒ **every cell of a line carries the identical label vector.** Two consequences follow: grouped
splitting becomes mandatory (below), and per-cell MSE is not the honest metric (next section).

The result is stored as `obsm["Y_ctrp"]` `(n_cells, K)` float32, NaN where unscreened, with the
length-K column→drug map in `uns["ctrp_drugs"]` and the score name in `uns["ctrp_score"]`. A drug
column is kept only if screened on ≥ `--min-cell-lines` overlapping lines (default 50; the K=545 runs
used `--all-drugs`, i.e. min 0). The legacy flat columns `obs["viability_<drug>"]` /
`obs["train_mask_<drug>"]` hold **whichever score was selected** — the `viability_` prefix is
historical. Cancer type is never a label or a feature — it only colors the UMAPs in
[Step 02](02-preprocessing-and-embeddings.md).

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

| Name | Where | Definition | K=545 scGPT (`mean_pv`) |
|---|---|---|---|
| **Entry-pooled MSE** | `best_val_mse`, `history.csv` | `Σ sq·M / Σ M` over all observed entries | **0.0105** |
| **Macro per-drug MSE** | `model_mean_mse` / `baseline_mean_mse` | per-drug `Σ_cells sq_k / n_k`, then `np.nanmean` over drugs (equal weight per drug) | **0.0103** |

Only the **macro per-drug** numbers feed the **per-drug-mean baseline** comparison
(`train_multitask._per_drug_constant_mse`). That baseline is a null model: for each drug it predicts
the constant train-set mean score over that drug's observed cells. "**Heads beating baseline**"
then counts drugs whose model per-drug val MSE beats that constant (scGPT 142/545, PCA 97/545 —
[Step 05](05-multitask-results.md)).

**Beating that baseline, not the raw MSE, is the honest metric** — and the target scale decides how
obvious that is:

- On the legacy `mean_pv`, absolute MSE ≈ 0.01 looked impressively tiny but was meaningless: labels
  cluster near 1.0, so the constant baseline already reached ≈ 0.0097.
- On **`auc_z` the two coincide by construction**: a z-scored target has unit variance, so predicting
  a drug's mean scores **MSE ≈ 1.0** and any MSE below 1 is real per-drug signal. (Observed: 0.92 for
  the val baseline, since val holds only 27 lines.) This readability is a direct benefit of the
  z-scoring.

---

## Why splits must be cell-line-grouped

Because the label is **constant within a cell line**, a random cell-level split would put cells of
the same line in both train and val. Those cells share the same label and a near-identical tissue
signature, so the model could memorize the per-line label instead of learning response.
`create_splits.py` avoids this by partitioning **whole cell lines** 70/15/15: `split_ctrp` is
drug-agnostic (one split is leakage-free for all heads at once), and `split_paclitaxel` is the
single-task version. This is not a detail — it is the control that exposes the PCA-vs-scGPT
overfitting gap in [Step 04](04-single-task-results.md).

---

## Model architecture & regularization (`OncoMLP.py`, 25.05.2026)

`OncoMLP` is a deliberately **small** MLP (the plan asks for the smallest functional model):
input dropout → `[Linear → LayerNorm → GELU → Dropout]` per hidden layer → `Linear(→ output_dim)`.
Defaults encode specific choices for this regime:

- **`hidden_dims` = (128,64) for both reps** (`train_multitask.DEFAULT_HIDDEN_DIMS`, matched
  14.06.2026) — a **matched trunk** so only the input representation (and its first projection)
  differs, making PCA vs scGPT a fair comparison. (Earlier runs used (64,32) for PCA, which
  handicapped it — see [Step 05](05-multitask-results.md).)
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

Both the CLI and `notebooks/07_training.ipynb` drive one training run through the same
`train_multitask.train_rep(...)` function (datasets → per-drug-mean baseline → `OncoMLP` → `train_model`
→ `save_run`), returning the run dir, history, and per-drug MSE arrays. The notebook is the
**reproducible PCA-vs-scGPT comparison**: it trains both reps at the matched 512-d width and writes
the comparison figures/tables to `notebooks/outputs/`. Because both paths call `train_rep`, the
notebook and command line cannot diverge.
