# Step 03 — Model & training design

*Part of [OncoTox project progress](../project_progress.md). Covers: what a single training
example is (input, output, target, mask), the supervised learning paradigm and the weak-supervision
it rests on, the masked-loss formulation, the exact MSE definitions, and the model architecture —
each tied to the code that implements it (`scripts/model/`, `scripts/training/`).*

---

## The learning problem — weakly-supervised, fully-supervised regression

The downstream task is a **continuous regression**: map one cell's transcriptomic representation to
a drug-response scalar. It is **fully supervised** — `train_model` in
`scripts/training/training_utils.py` optimizes a (masked) **MSE** or **MAE** loss directly against
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

### The drug panel is a training-time choice, not a property of the target file (12.08.2026)

`ctrp_to_h5ad` *can* restrict the target matrix to a named compound list — `target_drugs`, exposed as
`--drugs` — so which drugs reach the h5ad is a real choice rather than a default.

**Decision (Selin, 12.08.2026): it is not used. The panel is applied at training and nowhere earlier.**
The reason is that the panel's entire effect is on the model: it sets `output_dim`, that is, how many
heads share one trunk. It is a statement about what the network is asked to predict, not about what the
data contains. So `Y_ctrp` / `M_ctrp` keep the full screened catalog, and the panel is a column selection
made by `MultiDrugDataset(drugs=…)` when a run is configured.

Two consequences worth having written down:

- **Changing the panel requires no preprocessing re-run.** One targets h5ad serves any panel, which is
  why `literature_panel.ipynb` can rebuild the panel under the freeze while the h5ads cannot be rebuilt.
- **Filtering upstream would have moved the splits.** A cell line is eligible for splitting if any of its
  cells carries at least one observed label — `has_any_label = M.any(axis=1)`, `create_splits.py:164`.
  That test is taken over the *width* of `M`, so narrowing `M` from ~545 columns to 11 would re-evaluate
  eligibility against the panel and could drop lines that are screened but not against a panel compound,
  silently redrawing `split_ctrp`. Found by the code-quality session and verified here against
  `create_splits.py`; it is a consequence of the decision, not its motivation.

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

### The uncentred target is handled the same way in every training path

Two mechanics follow from a target sitting near 0.9 rather than 0. Both are silent if missed — nothing
errors, the model simply trains against a handicap — so they are recorded with the code that causes them.

1. **Weight decay is applied to the weight matrices, not to the biases and normalization parameters.**
   `optim.Adam(..., weight_decay=1e-3)` decays **every** parameter it is given, head biases included,
   and each head's bias must sit near its drug's mean. The grouping is the standard one —
   HuggingFace `transformers`' `Trainer.create_optimizer` (`no_decay = ["bias", "LayerNorm.weight"]`),
   inherited from the BERT reference implementation — and it is on because it is the convention, not
   because this target invented a need for it: `TrainConfig.no_decay_bias_and_norm`, default **on**,
   implemented in `training_utils._decay_param_groups`.
2. **Head biases are initialized at the fitting-fold per-drug means** (`OncoMLP.init_head_bias_`), so
   the model does not spend its first epochs climbing from `nn.Linear`'s ±0.125 to ~0.9. Initializing a
   final layer at the base rate is standard for exactly this reason — Lin et al., *Focal Loss for Dense
   Object Detection*, ICCV 2017 §4.1, done there for a class prior rather than a regression mean. The
   means are taken **per cell line** (`cv.per_drug_line_mean`), not per cell, because the label is
   broadcast across a line's cells.

⚠️ **It starts the model near the null predictor's *level*, not at the null predictor.** The head's
weight rows are still randomly initialized, so at initialization predictions already scatter with a
standard deviation of ≈0.31 across cells, against a true across-line spread of order 0.17 on `auc_cc`;
the mean prediction lands at 0.76 and 0.96 for a requested 0.90 at seeds 42 and 0. Measured on synthetic
unit-variance input, `init_spread.py` under audit 08 — an architecture property, no data read. Starting
genuinely *at* the null would also require shrinking the output layer's weight initialization, which
Lin et al. do (σ = 0.01) and **this project has not decided**.

⚠️ Neither mechanic carries over to `ln_ic50_cc` unchanged: its per-drug means are spread across a
log-concentration scale instead of clustering near 0.9. That has to be checked before that measure is
trained, not assumed.

**Until 12.08.2026 only one of the three training paths did either of these** — `cv.oof_predictions`,
which `4a_percell_training` drives. `train_multitask.cv_evaluate` (the CV behind the 8-run matrix) and
`train_multitask.train_rep` (the fixed-split path) initialized no head bias and ran with
`exclude_output_from_decay` at its `False` default, so on an uncentred target the matrix trained against
an offset the panel run did not, and the two were never the same experiment. Both now take
`init_head_bias=True` by default. The record, including what the superseded flag actually did:
[Corrections](corrections-and-dead-ends.md#the-two-uncentred-target-mechanics-ran-in-one-training-path-of-three).
**Takes effect at R4.**

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

- **Per-element error** `sq = (pred − y)²` (MSE), or `l1_loss` for `--loss mae`
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

**Confirmed as standing until MIL (Selin, audit 09, 12.08.2026)** rather than fixed in the loss now —
MIL dissolves it structurally, one bag is one line is one example, so a line-balanced weighting would be
machinery built to be discarded.

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

### The loss is plain masked MSE, and stays that way until MIL (audit 09, 12.08.2026)

**Decided by Selin, 12.08.2026.** Three things were on the table — a spread term, a ranking term, and
tuning their weights — and none enters the objective:

- **Ranking cannot go in this loss.** Spearman is piecewise constant in the predictions, so its gradient
  is zero almost everywhere. A differentiable surrogate exists and is standard (RankNet; LambdaRank,
  Burges et al.), but it needs *one score per ranked item*, and under per-cell training a batch is
  cells, so it would rank cells against each other on broadcast labels — the wrong object. It becomes
  well-posed only under MIL, where one bag is one cell line, and it is deferred to there
  (`notebooks/4b_mil_training.ipynb`).
- **Spread does not need the loss.** Predictions are shrunk because an MSE-optimal predictor must
  shrink in proportion to how little it knows (`pred_std ≈ ρ × true_std`); that is correct behaviour,
  not a defect. Where calibrated *values* are reported, a per-drug linear recalibration fitted on the
  training folds fixes it exactly, costs nothing, and leaves every ranking untouched because it is
  monotone. Putting a spread penalty in the objective would buy the same thing approximately, with a
  hyperparameter.
- **Tuning three loss weights** on ~153 independent labels is selection on very little data, and it
  would land at the same time as MIL, making neither attributable.

What replaces them is **measurement**: order, order-at-the-top, values and spread are four things to
*look at*, computed in `notebooks/5_evaluation.ipynb`, not four terms to optimize
([Step 05](05-multitask-results.md)).

**Which losses get compared, and when.** **MSE / MAE**, each with the density weighting at
`alpha ∈ {off, 0.5, 1.0}` — **six arms** — on the per-cell architecture; ranking losses only under MIL.
Two conditions carry over from the last comparison, which failed because it lacked them: the decision
rule is fixed *before* the run, and ≥3 seeds are needed for any difference to be readable against the
±0.04 seed band.

> ⚠️ **Huber was dropped from this grid (Selin, 12.08.2026).** It used to be the third loss here. Its
> position on the robustness axis is set entirely by `beta`: at the current `huber_beta = 0.05` it is
> linear above a residual far below the typical one, so it behaves close to MAE and the grid would carry
> **two near-duplicate columns**. Choosing a principled `beta` instead is worse, not better — it imports
> a new unsourced constant into a comparison whose entire purpose is to *be* the justification for the
> loss. **MSE and MAE already bracket the axis** (pure quadratic against pure absolute), so Huber adds
> cost without adding a corner. This retires review item 9C's `huber_beta` question with it.
>
> ✅ **Huber has since left the code as well** (`f16b3ec`, merged 13.08.2026): `--loss` is `{mse, mae}`,
> and `TrainConfig.huber_beta`, the name check and the `smooth_l1_loss` branch are all gone. The
> mechanics described elsewhere on this page were rewritten in the same change.
> *(This block previously read "the code still exposes `--loss huber`", which was accurate when the
> comparison dropped Huber but the code had not. It became false the moment that branch merged, and the
> sentences it protected were exactly the ones that then needed rewriting — an inversion worth keeping
> visible: text left alone **because** it was accurate is the text that goes stale first.)*

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
stamps a `fold` column onto every row of `outputs/legacy/panel_void_8drug/panel_oof_predictions.csv`. It raises if a line
appears in two folds, or if a predicted line is claimed by none.

### The early-stopping set is nested inside the training lines (12.08.2026)

⛔ **Until 12.08.2026 the fold's held-out lines were the early-stopping set as well as the scored set.**
`train_model` restores the checkpoint with the lowest validation MSE, and `oof_predictions` passed it
the held-out fold and then predicted that same fold — so every out-of-fold prediction came from the
epoch that best fit the lines it was about to be scored on, and `cv_evaluate`'s `best_val_mse` is a
minimum over epochs on its own scored fold. Both are optimistically biased. **This is the same defect
found in the DrEval benchmark on 14.07.2026 and fixed there** (`ee07b00`, by switching to DrEval's own
`sp['validation']`) — the fix was never carried into the code that produces the headline numbers.
Recorded in [Corrections](corrections-and-dead-ends.md).

The bias is not uniform across the two arms, which is what makes it more than cosmetic here. In the
last run on record the selected epoch was `[1,1,3,1,1]` for PCA against `[10,11,2,21,4]` for scGPT: the
PCA arm's scored checkpoint was picked, on the scored lines, from among near-tied epoch-1 states, while
scGPT's came from a real trajectory. The PCA-vs-scGPT margin therefore compared two differently-selected
quantities. (Those epochs are from the void 8-drug panel and are cited as evidence about the mechanism,
not as a result.)

**The fix (Selin, 12.08.2026):** `scripts/training/cv.py::inner_holdout` withholds **15 % of each
fold's training lines** as the early-stopping set, drawn by `GroupShuffleSplit(n_splits=1,
test_size=0.15)` — whose `test_size` is a proportion of *groups*, so it is 15 % of lines and the cell
counts fall where they fall. The scored fold is then never seen by the checkpoint that predicts it.
It costs training data: ~104 fitting lines per fold where there were 122.

- **0.15 is arbitrary** and is documented as such. It is the conventional size of a validation slice;
  nothing in this data sets it. The alternative considered and rejected was reusing a whole
  neighbouring `GroupKFold` fold as the early-stopping set — no new fraction and no new seed, but 20 %
  out instead of 15 %, on a panel where the label side is already the binding constraint.
- **The inner split has its own seed** (`INNER_SPLIT_SEED = 42`), deliberately *not* `TrainConfig.seed`.
  Repeat runs at different training seeds then early-stop on the same lines, so the four-identical-runs
  check that established the `mps` nondeterminism keeps measuring one thing
  ([Corrections](corrections-and-dead-ends.md)); tying it to the run seed would mix kernel
  nondeterminism with inner-split variability, and would also stop the two representations from being
  compared on identical inner splits unless both were run at the same seed.

**Takes effect at the sweep.** Every number on record predates it, and all of them will move down: the
optimistic selection is removed and the training set shrinks at the same time.

### What is and is not standard about the cross-validation

After the fix, the *model* fit is a textbook grouped 5-fold CV — every line held out exactly once, the
checkpoint chosen on lines the scored fold does not contain. Four things around it are not, and are
stated here so no reader has to infer them:

1. **It covers 153 of the 180 labelled lines.** The 27 fixed `test` lines are outside CV entirely
   (`eligible_splits=("train","val")`), and the 18 lines with no CTRPv2 label are outside everything.
2. **The folds are unshuffled and unseeded.** `GroupKFold` assigns whole lines greedily to balance
   *cell* counts, so the folds hold out 29/31/31/31/31 lines rather than equal numbers
   (`outputs/legacy/panel_void_8drug/panel_training_folds.csv`). Deterministic, but not the shuffled `KFold` the phrase
   usually implies.
3. **It is one partition, not repeated CV.** The fold-to-fold spread quoted in `diagnostics.ipynb`
   comes from a single draw.
4. **The representation is fitted outside the loop.** HVG selection and the PCA rotation are computed
   once over all 53,513 cells, before any fold exists; a fully nested CV would refit both per fold.
   Open under review item 7 — [Step 02](02-preprocessing-and-embeddings.md#what-transform-pca-sees--corrected-05082026).

A fifth, which is not a property of the CV but bears on how independent it is: the trunk width, dropout
and epoch budget were settled on the fixed `val` lines ([Step 05](05-multitask-results.md)), and CV then
pools train+val. Lines that informed the architecture are inside the folds it is scored on.

**Why the column exists.** Any baseline fitted *after* the fact — DrEval's mean-effects predictor above
all — has to be fitted on the folds a prediction did **not** come from. Without the fold label the only
option is to fit it on the same out-of-fold rows it will then be subtracted from, which lets held-out
labels define the baseline they are scored against. `scripts/evaluation/dreval_normalize.py` therefore
requires the column and refuses to run without it.

⚠️ The committed `panel_oof_predictions.csv` predates this and has no `fold` column, so that script
correctly raises on it. It becomes runnable when `4a_percell_training.ipynb` re-runs at R4 of the sweep.

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
the end rather than the last-epoch weights. Which lines that validation set is drawn from is the
subject of *The early-stopping set is nested inside the training lines* above — under cross-validation
it is no longer the scored fold. The single entrypoint `train_multitask.py` exposes these
as flags (`--use-rep`, `--drugs`, `--batch-size 128`, `--epochs 50`, `--lr`, `--weight-decay`,
`--dropout`, `--input-dropout`, `--loss {mse,mae}`, `--hidden-dims`, `--seed`); run artifacts are
written by `create_run_dir`/`save_run` ([Step 05](05-multitask-results.md)).

Both the CLI and `notebooks/4a_percell_training.ipynb (§B)` drive one training run through the same
`train_multitask.train_rep(...)` function (datasets → per-drug-mean baseline → `OncoMLP` → `train_model`
→ `save_run`), returning the run dir, history, and per-drug MSE arrays. The notebook is the
**reproducible PCA-vs-scGPT comparison**: it trains both reps at the matched 512-d width and writes
the comparison figures/tables to `notebooks/outputs/`. Because both paths call `train_rep`, the
notebook and command line cannot diverge.

---

## These hyperparameters are not worth tuning (ablated 13.07.2026)

> ⛔ **Every MLP number in this section is void (found 12.08.2026, audit 08), and the ridge row is not.**
> The notebook behind both tables chose each model's checkpoint on the fold it then scored — the same
> early-stopping leak audit 07 fixed in the pipeline — on the retired `auc_z` target, over five drugs
> chosen by the discredited learnability gate, and it can no longer be run because those target files
> were deleted. `RidgeCV` has no early stopping, so **the baseline was never flattered and the MLP rows
> were**, which is the one direction that matters when the two are tied. What has to be re-derived and
> what does not:
> [Corrections](corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning).

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

⚠️ **The +0.06 is a margin between a flattered number and a clean one (12.08.2026).** Only the MLP rows
carry the early-stopping leak; ridge selects its penalty by generalized cross-validation inside the
training fold and was always honest. The MLP-minus-ridge differences in the table are therefore upper
bounds, and the PCA row — where the two are *equal* — can go negative once the leak is removed. Which
would say the per-cell model is beaten by a ridge on 153 line-mean vectors, a different claim from the
one written above. **The decision on the row that matters is recorded in [TODO](../TODO.md) item 8C:
re-derive trunk vs bare linear head, both representations, against ridge on the same folds, at R4.**
