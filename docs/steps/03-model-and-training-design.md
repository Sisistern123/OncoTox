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

**What the target *is* — its source, the two measures, the counts and the provenance — is stated once
in [Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026).
This section covers only what that choice forces on the *model*.**

The label is a **bulk, per-(cell line × drug)** quantity, never measured per single cell. Selected with
`--score`; one targets h5ad per measure, so two can be trained head-to-head without rebuilding anything.

| `--score` | Direction | Scale |
|---|---|---|
| **`auc_cc`** (default) | higher = **more resistant** | ~0.02–1.83, median 0.925. **1.0 is the no-effect level** — the curve fit pins the low-concentration asymptote to the vehicle |
| `ln_ic50_cc` | lower = **more sensitive** — the opposite | log concentration, −11.4 to 8.6 |

⚠️ **The two run in opposite directions.** Anything that reads a sign — a loss term, a plot axis, a
"most sensitive lines" table — is wrong for one of them unless it is written to ask which is loaded.

### Every cell of a line carries the identical label

`ctrp_to_h5ad.py` pivots to a (cell line × drug) matrix and broadcasts each bulk value to every cell of
the matching line (`Y_full = cl_drug_matrix.reindex(cell_line_norm.values)`). Two consequences run
through the rest of this file:

- **Grouped splitting is mandatory** ([below](#why-splits-must-be-cell-line-grouped)) — a random cell
  split would put the same label on both sides.
- **Per-cell MSE is not the honest metric** (next section) — a line with 1,990 sequenced cells would
  count 35× a line with 56, for one measurement.

It is also what makes **research question 2 structurally unanswerable under this architecture**: the
objective penalises exactly the within-line variation Q2 asks about. MIL is the instrument for that
question, not a capacity lever — see [TODO](../TODO.md).

### Two mechanics forced by a target centred near 0.9 rather than 0

Both are silent if missed — nothing errors, the model simply trains against a handicap — so they are
recorded with the code that causes them.

1. **The output layer must be excluded from weight decay.** `optim.Adam(..., weight_decay=1e-3)`
   (`scripts/training/training_utils.py:179`) decays **every** parameter, head biases included. Each
   head's bias must sit near its drug's mean — around 0.9 on `auc_cc` — and decay pulls it toward 0.
   Biases and LayerNorm parameters therefore go in a `weight_decay=0` group:
   `TrainConfig.exclude_output_from_decay`, default **off** so older runs reproduce unchanged.
2. **Head biases must be initialized at the train-fold per-drug means.** `nn.Linear` initializes them
   within ±0.125, so the model otherwise starts far from the null predictor rather than at it.

⚠️ Neither carries over to `ln_ic50_cc` unchanged: its per-drug means are spread across a
log-concentration scale instead of clustering near 0.9. That has to be checked before that measure is
trained, not assumed.

### Per-drug Spearman is affine-invariant — and what that implies

For a *fixed* model, any within-drug affine transform of the target scores **identically**. So a
per-drug rescaling can never show up in the metric: it is a **loss-weighting scheme in disguise**, which
is what retired `auc_z` ([Corrections](corrections-and-dead-ends.md#auc_z-as-the-training-target)).

It also means `auc_cc` vs `ln_ic50_cc` is a real comparison rather than a relabelling — the two are
*not* affine images of one another, so they can rank cell lines differently.

### Is AUC the right target at all? — answered 11.08.2026

Raised 27.07.2026 by Selin's supervisor and by DrEval itself, on two grounds. Both now have an answer:

- **AUC conflates potency with efficacy.** True of any AUC: the area mixes *how little drug is needed*
  with *how much killing is achievable*, so two pharmacologically opposite compounds can share one. The
  pipeline now carries `ln_ic50_cc` beside `auc_cc` from the **same** CurveCurator fit, so the two are
  compared rather than argued about. CTRPv2 publishes no IC50 of its own — it exists only in the re-fit.
- **The tested concentration range is not standardized.** Across the 545 compounds the top test
  concentration spans **0.13 µM to 600 µM**. This does not break comparisons *between* drugs, since the
  metric is within-drug — but it bites *within* one: a range badly matched to a compound's potency
  compresses every line toward one end of the curve, shrinking the spread the model is asked to predict
  toward the noise floor. **Unresolved, and recorded as a limitation rather than a defect**;
  re-integrating over a per-drug common window was considered and not adopted ([Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)).

The retired measures — `auc`, `auc_z`, `mean_pv` — their definitions, the head-to-head comparison that
once chose between them, and the normalisation defect that voided `auc`, are all in
[Corrections](corrections-and-dead-ends.md#the-auc-target-was-divided-by-the-wrong-quantity).

---

## Mask `M` — the sparsity-handling mechanism (plan sub-goal 2)

Most (cell line × drug) pairs were never assayed, so the label matrix is sparse. `ctrp_to_h5ad.py`
records `obsm["M_ctrp"]` `(n_cells, K)` bool, **True iff that cell's line was actually screened
against that drug**. `MultiDrugDataset` (`scripts/model/dataset.py`) fills missing `Y` with 0.0 only
so the tensor is finite, then carries `M` alongside so the loss can ignore those zeros.

This is the masked-loss core that the plan's sub-goal 2 calls for, and it is what will generalize to
the cross-database block-sparse matrix in [Step 06](06-planned-work.md#a-cross-database-integration).

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

### ⚠️ Open defect — the loss weights cell lines by how deeply they were sequenced

The loss is a mean over observed **cells** (`scripts/training/training_utils.py:90`), but the label is
per **cell line**. Two distributions therefore multiply into the weight each line carries:

| | min | median | max |
|---|---|---|---|
| cells per line | 56 | 226 | 1,990 |
| observed drugs per line | 38 | 465 | 505 |

Combined, `NCIH2110_LUNG` carries **4.37 %** of the total loss and `PANC1_PANCREAS` **0.054 %** — a factor
of **82**. Equal shares would be 0.56 % each; the top 10 lines carry **18.3 %** instead of 5.6 %. With
~150 independent examples, letting one line count 82× another is indefensible, and most of the imbalance
is the cell count — a sequencing artifact carrying no information about the label. (The z-scoring code was
careful about exactly this — "statistics are computed per cell line, not per cell" — the loss was not.)

**Reweighting the existing objective does not fix it:** line-balanced reweighting was tested and is empty
([Corrections](corrections-and-dead-ends.md#line-balanced-reweighting-will-help)), which in hindsight is
forced, since the ridge-on-line-means control *is* the fully line-balanced limit and it ties the PCA MLP.
MIL removes the defect structurally — one bag is one line is one example — which is one of the reasons it
is next ([TODO](../TODO.md)).

**Sample weighting rides in the mask (27.07.2026).** `M` is a float tensor, so substituting a weight
at each observed entry for the 0/1 flag turns `Σ(sq·M)/ΣM` into `Σ(w·sq)/Σw` exactly — the weighted
objective, with **no change to the training loop**. `scripts/training/cv.py` does this when
`density_weighting=True`: it fits one weight function per drug inside each fold, on that fold's
training *lines* only (`density_weighting.fit_weight_fns`), and writes
`weight_matrix(...)` into `dataset.mask`. Unobserved entries stay 0 and remain excluded.
Consequences worth stating: the per-epoch val MSE printed during a weighted run is the **weighted**
MSE — the right quantity for early stopping, but not comparable to an unweighted run, so the
comparable numbers are recomputed afterwards from the out-of-fold predictions.
Anatomy of the objective: `docs/figures/loss_01_objective.png`; the weight curve:
`loss_02_weights.png`.

**Two aggregations are reported — do not conflate them:**

| Name | Where | Definition |
|---|---|---|
| **Entry-pooled MSE** | `best_val_mse`, `history.csv` | `Σ sq·M / Σ M` over all observed entries |
| **Macro per-drug MSE** | `model_mean_mse` / `baseline_mean_mse` | per-drug `Σ_cells sq_k / n_k`, then `np.nanmean` over drugs (equal weight per drug) |

Only the **macro per-drug** numbers feed the **per-drug-mean baseline** comparison
(`train_multitask._per_drug_constant_mse`). That baseline is a null model: for each drug it predicts
the constant train-set mean score over that drug's observed cells. "**Heads beating baseline**"
then counts drugs whose model per-drug val MSE beats that constant (scGPT 142/545, PCA 97/545 —
[Step 05](05-multitask-results.md)).

**Beating that baseline, not the raw MSE, is the honest metric**, and on `auc_cc` the raw number is
actively misleading: the labels cluster near 0.9 with a small spread, so an absolute MSE around 0.01
looks impressively tiny while the constant baseline already reaches nearly the same value. A
standardized target would make the two coincide by construction — unit variance means the baseline
scores 1.0 and anything below is signal — but that readability was the one real benefit of `auc_z`, and
it did not survive its costs ([Corrections](corrections-and-dead-ends.md#auc_z-as-the-training-target)).
Read the baseline comparison, never the MSE alone.

---

## Why splits must be cell-line-grouped

Because the label is **constant within a cell line**, a random cell-level split would put cells of
the same line in both train and val. Those cells share the same label and a near-identical tissue
signature, so the model could memorize the per-line label instead of learning response.
`create_splits.py` avoids this by partitioning **whole cell lines** 70/15/15: `split_ctrp` is
drug-agnostic (one split is leakage-free for all heads at once), and `split_paclitaxel` is the
single-task version. This is not a detail — it is the control that exposes the PCA-vs-scGPT
overfitting gap in [Step 04](04-single-task-results.md).

### Out-of-fold predictions carry their fold (12.08.2026)

The cross-validation is 5-fold `GroupKFold` over the train+val lines
(`scripts/training/cv.py::oof_predictions`), and each cell line is held out by **exactly one** fold —
which is what makes a fold label well defined at the line level at all.

`oof_predictions` records each fold's held-out lines in its log, and `line_level_predictions(folds=...)`
stamps a `fold` column onto every row of `outputs/panel/panel_oof_predictions.csv`. It raises if a line
appears in two folds, or if a predicted line is claimed by none.

**Why the column exists.** Any baseline fitted *after* the fact — DrEval's mean-effects predictor above
all — has to be fitted on the folds a prediction did **not** come from. Without the fold label the only
option is to fit it on the same out-of-fold rows it will then be subtracted from, which lets held-out
labels define the baseline they are scored against. `scripts/evaluation/dreval_normalize.py` therefore
requires the column and refuses to run without it.

⚠️ The committed `panel_oof_predictions.csv` predates this and has no `fold` column, so that script
correctly raises on it. It becomes runnable when `3_panel_training.ipynb` re-runs at R4 of the sweep.

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

Both the CLI and `notebooks/2_training.ipynb` drive one training run through the same
`train_multitask.train_rep(...)` function (datasets → per-drug-mean baseline → `OncoMLP` → `train_model`
→ `save_run`), returning the run dir, history, and per-drug MSE arrays. The notebook is the
**reproducible PCA-vs-scGPT comparison**: it trains both reps at the matched 512-d width and writes
the comparison figures/tables to `notebooks/outputs/`. Because both paths call `train_rep`, the
notebook and command line cannot diverge.

---

## These hyperparameters are not worth tuning (ablated 13.07.2026)

`notebooks/archive/ablations_and_rescue.ipynb` sweeps four model-side knobs on the 5 learnable drugs, scored with
out-of-fold per-drug Spearman (the metric of [Step 05](05-multitask-results.md)). **Every axis is flat,
and the defaults above are at or within noise of the best setting on all of them:**

| Knob | Range tested | PCA | scGPT |
|---|---|---|---|
| Regularization (`dropout`/`input_dropout`/`weight_decay`) | none → (0.7, 0.2, 1e-2) | 0.42–0.44 | 0.44–**0.49** |
| Capacity (`hidden_dims`) | 74,629 → 2,565 params | 0.41–0.43 | 0.44–**0.49** |
| `batch_size` | 32 / 128 / 512 | 0.43–0.44 | 0.46–**0.49** |
| Sample reweighting | line-balanced, focus-extremes | 0.41–0.43 | 0.48–0.49 |

> ⚠️ **These ablations ran on the *corrected* setup, so they show the knobs do not *improve* it — not
> that the knobs could not have *fixed* the K=545 collapse.** That was tested separately, and removing
> regularization recovers ~70 % of it. Why that is a symptom fix rather than a cause fix, with the full
> rescue table:
> [Corrections](corrections-and-dead-ends.md#the-model-is-over-regularized-or-too-small).

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
