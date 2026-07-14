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

### How `auc_z` is computed, step by step

CTRPv2 screens each (cell line × drug) over a **concentration series** and fits a sigmoid to it. The
pipeline (`ctrp_to_h5ad.py`) turns that into one number per (cell line × drug) in four steps.

**Step 1 — take the curve fit, not the raw dose points.** Read `area_under_curve` from
`v20.data.curves_post_qc.txt` (`_load_score_values`), one row per (`experiment_id`, `master_cpd_id`).
This is the integrated area under the **fitted** sigmoid, in log2-concentration space. It is not a
percentage and **not bounded by 1** — for a 16-point grid it runs roughly 0–16.

**Step 2 — divide out the concentration grid → `auc`.** CTRP integrates rather than averages, so the
raw area grows with how many concentration points were fitted. `conc_pts_fit` varies **8–29** (usually
16), so a raw AUC of 13 means different things for different compounds. Dividing removes it:

```
auc = area_under_curve / conc_pts_fit        # ~0-1, viability-like: the mean height of the fitted curve
```

The result reads like an average survival fraction across the tested range: **low `auc` = the drug
kills that line, high `auc` = the line survives it.** (Observed: median 0.877, p05 0.536, p95 1.071 —
values slightly above 1 occur where a line grows faster than the vehicle control.)

**Step 3 — one value per (cell line, drug).** Replicate experiments are averaged
(`_build_drug_table`: `groupby(["ccl_name_norm","cpd_name_norm"]).mean()`), giving the (180 × 545)
matrix the rest of the pipeline works from.

**Step 4 — z-score *within each drug*, across cell lines** (`_zscore_per_drug`):

```
center[d] = mean of auc[:, d] over the cell lines screened against drug d   # -> uns["ctrp_score_center"]
scale[d]  = std  of auc[:, d] over those same lines                         # -> uns["ctrp_score_scale"]

auc_z[l, d] = (auc[l, d] - center[d]) / scale[d]
```

Note **which axis is standardized**: the statistics are taken **down the cell-line axis, separately
for every drug** — *not* over the whole matrix, and *not* per cell line. Each drug ends up with mean 0
and std 1 across the panel, so `auc_z` answers **"how sensitive is this line to this drug, relative to
the typical line for this drug?"** −1 means one standard deviation more sensitive than that drug's
average line; +1 means one more resistant.

The statistics are computed **per cell line, not per cell** — a line with 500 cells must not pull the
mean toward itself more than a line with 50, since both are one measurement. Degenerate drugs (zero
spread) keep `scale = 1.0` rather than dividing by zero.

*Worked example — `dasatinib`* (`center = 0.631`, `scale = 0.155`): a resistant line with `auc = 0.85`
becomes `(0.85 − 0.631) / 0.155 = +1.41`; a sensitive line with `auc = 0.40` becomes **−1.49**. Both
are stored in `Y_ctrp` and broadcast to every cell of the line.

**The map is exactly invertible**, which is why `center`/`scale` are saved: `auc = auc_z * scale +
center` recovers the original units for any prediction (needed for cross-drug comparisons, which are
meaningless on the z-scale — see below).

### Why train on `auc_z` rather than `auc`

1. **It equalizes the heads in the shared loss** — the real reason. Per-drug `auc` spread ranges from
   **0.034 to 0.302** across the 545 drugs (a 9× span, ~80× in squared-error terms). Under a masked
   MSE on raw `auc`, the wide-spread drugs would supply nearly all the gradient to the shared trunk
   and the narrow ones almost none — purely because of their *units*, not their learnability.
   Z-scoring gives every head equal weight.
2. **It makes the metric readable.** With unit variance, the per-drug-mean null model scores **MSE =
   1.0** exactly, so any value below 1 is real per-drug signal. On raw `auc` every drug has its own
   null value and no number is interpretable on its own. (This is presentation, not learning.)
3. **It removes the per-drug potency offset** — the *weakest* of the three, and worth stating honestly:
   each drug is its own output row **with its own bias term**, so the head absorbs that offset either
   way. This argument is close to vacuous; do not lean on it.

### Measured: all three targets head-to-head (`notebooks/11_auc_vs_aucz.ipynb`, 13.07.2026)

The argument above was tested, not assumed. All three scores were trained with the identical model and
scored on **one common yardstick — the curve-fit `auc` ranking of cell lines** (`mean_pv` is *not* an
affine map of `auc` — within a drug they agree only ρ ≈ 0.72 — so scoring each model against its own
target would make the columns answer different questions; `auc`/`auc_z` *do* share a ranking exactly).
Out-of-fold Spearman on the 5 learnable drugs, ±95% bootstrap CI over the ~150 held-out lines:

| | `mean_pv` (legacy) | `auc` (curve fit) | `auc_z` (z-scored) |
|---|---|---|---|
| **K=5** · PCA | 0.450 [0.39, 0.50] | 0.439 [0.38, 0.49] | 0.424 [0.36, 0.48] |
| **K=5** · scGPT | 0.481 [0.42, 0.53] | 0.482 [0.42, 0.53] | **0.488** [0.43, 0.54] |
| **K=545** · PCA | **+0.027** [−0.04, 0.10] | **+0.016** [−0.06, 0.09] | **+0.378** [0.31, 0.44] |
| **K=545** · scGPT | **−0.070** [−0.14, 0.00] | **−0.087** [−0.15, −0.02] | **+0.430** [0.37, 0.48] |

(Pearson agrees throughout, within ±0.02 — neither metric is doing anything special. Per-drug spread is
wide, 0.10–0.65, so the dots in the notebook figure matter as much as the bars.)

**Two results, one of which corrects an earlier claim in these docs:**

1. **The z-scoring is load-bearing — and it is the *whole* effect.** At K=545 both unstandardized targets
   collapse to zero (or below) while `auc_z` holds at ~0.4. With per-drug spreads spanning 9× (0.034 →
   0.302, ~80× in squared error), the wide-spread heads monopolize the shared trunk's gradient **because
   of their units, not their learnability**, and nothing transferable is learned.
2. ⚠️ **The curve fit buys no measurable accuracy.** An earlier version of this section credited
   `area_under_curve` over the dose-averaged `mean_pv`. **That is falsified:** trained head-to-head,
   `mean_pv` and raw `auc` are statistically identical at *both* K (0.450 vs 0.439 at K=5; +0.027 vs
   +0.016 at K=545 — CIs fully overlapping). Reason 1 in the list above is **principled, not empirical**:
   keep the curve fit for the post-QC sigmoid, the confidence intervals, and because it is the metric
   family GDSC2 reports (which [Step 06](06-cross-database-integration.md) needs) — but do **not** claim
   it improves accuracy. It does not.

> **Decision: keep `auc_z` as the default** (`layout.DEFAULT_CTRP_SCORE`). At the full catalog it is the
> difference between a model that ranks cell lines and one that does not. **At K=5 all three targets tie**
> — on a small, spread-homogeneous subset `--score auc` is equally good and keeps predictions in native,
> interpretable AUC units.

**And it retro-diagnoses the project's central null result.** [Step 05](05-multitask-results.md) §3
("neither rep ranks cell lines", ρ ≈ 0 across 545 drugs) ran at **K=545 on `mean_pv`** — the grey column
above, which sits at **−0.070 / +0.027**. The table **reproduces that null on demand**, and shows it
vanishes the moment the heads are standardized. It was never clean evidence about scGPT vs PCA; it was
substantially an artifact of an **unstandardized multi-task loss**.

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

---

## These hyperparameters are not worth tuning (ablated 13.07.2026)

`notebooks/10_diagnosis.ipynb` sweeps four model-side knobs on the 5 learnable drugs, scored with
out-of-fold per-drug Spearman (the metric of [Step 05](05-multitask-results.md)). **Every axis is flat,
and the defaults above are at or within noise of the best setting on all of them:**

| Knob | Range tested | PCA | scGPT |
|---|---|---|---|
| Regularization (`dropout`/`input_dropout`/`weight_decay`) | none → (0.7, 0.2, 1e-2) | 0.42–0.44 | 0.44–**0.49** |
| Capacity (`hidden_dims`) | 74,629 → 2,565 params | 0.41–0.43 | 0.44–**0.49** |
| `batch_size` | 32 / 128 / 512 | 0.43–0.44 | 0.46–**0.49** |
| Sample reweighting | line-balanced, focus-extremes | 0.41–0.43 | 0.48–0.49 |

> ⚠️ **Scope correction (14.07.2026).** The four ablations above were run on the **K=5** setup. They show
> the knobs do not *improve* the corrected model — they do **not** show that the knobs could not have
> *fixed* the K=545 collapse. That claim was tested separately, and **one of them partially does**
> (`notebooks/outputs/ablations/rescue_k545.csv`). Applied to the **broken** setting (K=545, raw `auc`,
> scGPT, ρ = −0.063):
>
> | Intervention | ρ |
> |---|---|
> | heavy regularization | −0.091 |
> | line-balanced **sample** reweighting | −0.078 |
> | smaller model (74,629 → 16,645 params) | −0.053 |
> | batch size 32 | +0.027 |
> | **no regularization** | **+0.234** |
> | **per-drug (task) reweighting = `auc_z`** | **+0.433** |
>
> **Removing regularization recovers ~70% of the collapse.** This is mechanistically consistent: the
> failure is a **capacity competition** between heads — with dropout 0.5 the trunk cannot serve both the
> loud noisy drugs and the learnable ones, so the loud ones win. Adding capacity (dropping the
> regularizer) lets it fit both.
>
> **But it is a symptom fix, not a cause fix**, and the interaction proves it: on the **corrected**
> (`auc_z`) setting the *same* regularization is **optimal** (K=5: current 0.488 vs none 0.456). The model
> was never over-regularized in absolute terms — it was **over-regularized relative to a mis-weighted
> loss**. And the price is steep: without regularization the network memorizes the training lines (train
> MSE ≈ 0.01) and still reaches only half of what the weighting fix delivers.

**Decision: stop tuning the model** *(on the corrected loss)*. At ~150 independent labels, architecture
search cannot buy signal. Three findings are worth carrying forward, because each kills a
plausible-sounding "fix":

1. **Not over-regularized — the opposite.** With regularization *off*, PCA drives **train** MSE to ≈ 0.01
   (near-perfect memorization of the training lines) and still reaches only 0.42 out-of-fold. A model
   that can overfit that hard is not being suppressed by its regularizer. Heavy regularization is the
   only setting that hurts, and it does so via over-shrinkage (`pred_std` 0.33), not lost ranking.
2. **The prediction shrinkage is correct, not a defect.** `pred_std ≈ ρ × true_std` (scGPT: 0.47 vs
   ρ = 0.48) is exactly what an MSE-optimal predictor must do — the conditional mean shrinks toward the
   prior in proportion to how little signal exists. Loosening dropout to "fix" it would *raise* MSE. To
   report in AUC units, divide by ρ (Spearman is unchanged).
3. **Line-balanced reweighting is principled but empty.** The entry-pooled loss lets a 500-cell line pull
   10× harder than a 50-cell line, though both carry exactly **one** independent label — worth fixing on
   principle, but it changes nothing (scGPT 0.485 vs 0.488). In hindsight this is forced: the ridge
   control below **is** the fully line-balanced limit.

### The baseline that actually binds: ridge on 150 line-mean embeddings

| Model | PCA | scGPT |
|---|---|---|
| `OncoMLP` (128,64) | 0.428 | **0.487** |
| linear head (`hidden_dims=()`) | 0.412 | 0.438 |
| **RidgeCV on the 150 cell-line mean embeddings** | **0.428** | **0.428** |

A ridge regression on one mean vector per cell line — **no single cells, no network, seconds to fit** —
**ties the PCA MLP** and comes within 0.06 of the scGPT MLP. This is the direct consequence of the target
resolution: the label is per (cell line × drug), so there are **~150 independent training examples**, and
the 34k cells are an illusion of sample size.

> **Decision: `RidgeCV` on line-mean embeddings is the baseline to beat from now on** — not the
> per-drug-mean null, which is far too weak a bar. A result that does not clear ridge is not a statement
> about single-cell modelling. The one thing that currently does clear it is **scGPT + a hidden layer**
> (+0.06): scGPT's own linear head drops to 0.438, so it genuinely *needs* the nonlinearity, while PCA
> gains nothing from one. That asymmetry is the strongest evidence to date that the scGPT embedding
> encodes something PCA's does not.
