# Step 05 — Multi-task masked loss (545 CTRPv2 drugs) & run versioning

*Part of [OncoTox project progress](../project_progress.md). Covers: the multi-task masked-loss
model over all 545 CTRPv2 drugs, its results vs. the per-drug-mean baseline, and the run-versioning
ledger that records every training run.*

This moves from plan-Phase-2 (single-task) into plan-Phase-3 (masked-loss multi-task). Masked-loss
mechanics are in [Step 03](03-model-and-training-design.md). These are the **multi-task (all-drugs,
K=545) rows of the 8-run experiment matrix**
([index](../project_progress.md#experiment-matrix--pca-vs-scgpt)).

> **Scope — still 1 database, 1 score; "multi-task" here = multi-*drug* only.** Every one of the
> K=545 heads predicts the **same** metric from the **same** database
> (CTRPv2). This validates the masked-loss machinery on intra-CTRPv2 sparsity, but it is **not**
> the plan's ultimate multi-task goal, which is **cross-database** (CTRPv2 + PRISM + GDSC) **and
> multi-metric** (efficacy *and* toxicity). That integration — the real "combine all" — is
> [Step 06](06-cross-database-integration.md) and is **not yet started**. Do not read the 545-head
> run as "multi-task complete."

> ⚠️ **Legacy target score.** Every number on this page was trained on **`mean_pv`**, the only target
> until 13.07.2026; the default is now raw **`auc`** ([Step 03](03-model-and-training-design.md)),
> via `auc_z`, which was the default 13.07–27.07 and is retired.
> Absolute MSEs across the two scores are **not comparable** (a z-scored target has unit variance, so
> its baseline sits near 1.0 rather than 0.0097) — only *heads beating baseline* and the per-drug
> correlations transfer. Reproduce this page exactly with `--score mean_pv`.

---

## Multi-task masked loss over all 545 CTRPv2 drugs (26.05.2026)

The target artifacts (`obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`, the legacy flat
`viability_<drug>` columns) are defined once in [Step 03](03-model-and-training-design.md). What this
run adds is the split and the drug scope:

- `obs["split_ctrp"]` — **one drug-agnostic, cell-line-grouped 70/15/15 split** written by
  `create_splits.py` `run_multi()`, shared across all heads (leakage-free for every drug at once;
  a single shared split is only possible *because* the leakage control is at the cell-line level).
- **Drug-scope filter:** keep a drug only if screened on ≥ `--min-cell-lines` overlapping cell lines
  (default 50). This run used **`--all-drugs` (= min 0) → K = 545 drugs**.

**Run-time overlap reported by the pipeline:** **180 / 198** SCP542 cell lines overlap
CTRPv2 (180 = lines with actual post-QC measurements; the audit's 190 counts roster name-matches — see
[Step 01](01-datasets-and-harmonization.md)).

**`split_ctrp` distribution (one cell-line-grouped **70/15/15** split, shared by all heads):**

| split | lines | % of lines | cells | % of measured cells |
|---|---|---|---|---|
| train | 126 | 70.0% | 34,126 | 72.3% |
| val   | 27  | 15.0% | 7,121  | 15.1% |
| test  | 27  | 15.0% | 5,980  | 12.7% |

- 70/15/15 is the design target at the **cell-line** level (`create_splits._split_cell_lines`); the
  **cell** percentages differ slightly because lines carry different cell counts.
- `unassigned` = **18 lines / 6,286 cells** (SCP542 lines with no CTRP measurement; 198 → 180 measured).
- **Cross-validation** (`notebooks/07_training.ipynb` §2) **holds `test` out** and resamples only the
  153 train+val lines via 5-fold GroupKFold, so test is never seen in CV.

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
**Matched trunk + matched width.** Both reps use the **same** hidden layers `(128,64)`
(`DEFAULT_HIDDEN_DIMS`, set 14.06.2026) **and** the same **512-d** input (`X_pca` raised from scanpy's
~50 default to `add_pca.DEFAULT_N_COMPS = 512` on 27.06.2026), so the entire network — including the
first projection's parameter count — is identical and **only the representation differs**. This closes
the last comparison confound. (History: the original matrix used a `(64,32)` PCA trunk and a ~50-d PCA,
both of which handicapped PCA; the numbers below supersede those.)

### Metrics — what each number means

Every result below is one of these. They are reported on the **val** split (single fixed split) or as
**5-fold CV mean ± std** (test held out); read them together — MSE alone is misleading near a viability
of 1.0.

- **Masked (val) MSE** — the training objective. Per-cell squared error `(pred − viability)²`, averaged
  **only over observed `(cell × drug)` entries** (`mask = 1`); missing labels contribute nothing
  (`_masked_mean` in `training_utils.py`). For a single drug it's plain MSE. **Why ≈ 0.01 is
  misleading:** viability is per (cell line × drug), broadcast to all the line's cells, and clusters
  near 1.0 with tiny variance — so even predicting a constant scores ~0.01. Absolute MSE therefore says
  little; what matters is whether it beats the constant and whether it *ranks* lines (below). *Train*
  MSE is logged with dropout **active**, so it can sit below or above the dropout-free val MSE.

- **Per-drug-mean baseline** — the null model (`_per_drug_constant_mse`). For each drug it predicts that
  drug's **train-set mean viability** for every cell, then is scored on val. Because labels are near
  constant this is already a *strong* predictor, so it's the bar every head must clear; a head only
  counts as having *learned* response if it beats its own drug's constant.

- **Heads beating baseline (`heads_beat`)** — the **count** of the K = 545 drugs whose model val MSE is
  below their per-drug-mean baseline. Intuitive, but a **thresholded count of near-ties**: most heads
  have model ≈ constant (labels ≈ 1.0), so they sit on the decision boundary, and the per-fold baseline
  is recomputed from that fold's train lines. If a fold's held-out lines are collectively a little
  above/below the train mean, **hundreds of heads flip together** (common-mode), so the CV std is huge
  (±73–94; cf. √(K·p(1−p)) ≈ 11 if heads were independent — observed is ~8× that). **Treat as
  directional, not precise.**

- **Δmse (model − baseline)** — the **continuous** counterpart of heads-beating: the mean over drugs of
  `model_mse − baseline_mse`. **Negative ⇒ model better** than the constant on average; it is not
  thresholded, so it doesn't suffer the count's instability (its CV std is small relative to the mean).
  Reported as CV mean ± std; the per-fold `legacy/training_545_mean_pv/cv_folds.csv` also carries `median_delta` and `frac_beat`
  (= `heads_beat / n_total`).

- **Overfitting gap** — `val_mse − train_mse` at the best epoch (single-task). Larger ⇒ more
  memorization; the core hypothesis predicts scGPT < PCA. Same dropout-in-train caveat as above, so it
  is indicative, not exact.

- **Per-drug correlation (Spearman / Pearson)** — the metric that actually asks *does the model rank
  cell lines?* For each drug, predictions are averaged to one value **per held-out cell line** and
  correlated with the true per-line viability across lines (Spearman = rank, Pearson = linear).
  Restricted to drugs with **real response variance** (per-line true std ≥ 0.05) and ≥ 5 val lines —
  otherwise there is nothing to rank. Insensitive to the near-1.0 offset that dominates MSE.

- **5-fold GroupKFold CV (test held out)** — robustness wrapper: `GroupKFold(5)` over `Cell_line`
  resamples the 153 train+val lines into 5 train/val folds (no line on both sides), each retrained from
  scratch; we report **mean ± std**. The fixed `test` set is never touched, so **CV numbers are a
  stability check, not a test-set estimate**.

**The 8-run matrix (512-d, 27.06.2026; all share `split_ctrp`, n_train 34,126 / n_val 7,121).**
Per-drug-mean baseline: **~0.043** (K=1 paclitaxel, data-derived, rep-independent), **0.0097** (K=545).
Reproducible in `notebooks/07_training.ipynb`; run dirs `runs/20260627_1913xx_*` (see
`runs/runs_index.csv`).

**Single-task (K=1 paclitaxel) — the overfitting story** (gap = val − train, at the best epoch):

| Gene set | Rep | Train MSE | Val MSE | Gap (val−train) |
|---|---|---|---|---|
| `hvg5000` | scGPT | 0.037 | 0.041 | **0.004** |
| `hvg5000` | PCA | 0.011 | 0.045 | 0.033 |
| `all_genes` | scGPT | 0.032 | 0.045 | 0.013 |
| `all_genes` | PCA | 0.042 | 0.039 | −0.003 |

**All-drugs (K=545) — heads beating the per-drug-mean baseline:**

| Gene set | Rep | Val MSE | Heads beat baseline |
|---|---|---|---|
| `hvg5000` | scGPT | 0.0105 | 147 / 545 |
| `hvg5000` | PCA | 0.0103 | **169 / 545** |
| `all_genes` | scGPT | 0.0106 | 131 / 545 |
| `all_genes` | PCA | 0.0106 | **138 / 545** |

**Reading the results (matched trunk + matched 512-d width):**

- **Core hypothesis — supported (single-task, `hvg5000`):** scGPT's train/val gap is **0.004** vs
  PCA's **0.033** — scGPT overfits far less. Matching PCA to 512-d *sharpened* this: PCA's extra
  first-layer capacity lets it fit the train set harder (train 0.011) while val stays high (0.045),
  exactly the memorization the denoised scGPT prior is meant to avoid.
- **All-drugs — PCA competitive/better on raw accuracy:** heads-beating `hvg5000` **PCA 169 vs scGPT
  147**, `all_genes` **PCA 138 vs scGPT 131**; val MSEs are within 0.0003. scGPT does **not** win on
  absolute predictive metrics.
- **Net:** scGPT's robust, reproducible win is **lower overfitting**, not higher accuracy — and this
  now holds with input dimensionality matched, so it can no longer be dismissed as a capacity artifact.
- **Which heads are even learnable** is driven by coverage + response variance — see
  `notebooks/04_drug_coverage.ipynb`: the ≈16-line drugs (n_val 221) are the unreliable/hardest heads,
  while high-coverage high-variance drugs (docetaxel, gemcitabine, oligomycin a) are the easiest.

> ⚠️ **Gap-metric caveat.** Train MSE is logged with dropout (0.5) + input-dropout (0.1) **active**, so
> it can sit *below or above* the (dropout-free) masked val MSE; the gap is indicative, not exact. The
> `all_genes` rows early-stop very fast (best epoch 1–4), so their gaps are noisy — `all_genes`·PCA's
> **−0.003** reflects near-no learning + the dropout offset, not genuine negative generalization. The
> clean comparison is `hvg5000` single-task (scGPT 0.004 vs PCA 0.033).

### Is the difference real? — 5-fold cross-validation (27.06.2026)

The single-split numbers above rest on **27 val lines**, so they are point estimates. To test
robustness, `cv_evaluate` (`notebooks/07_training.ipynb` §2) runs **5-fold GroupKFold over `Cell_line`,
holding the fixed `test` set out** and resampling only the 153 train+val lines (~122 train / ~31 val
per fold). On `hvg5000`:

| Rep | Heads beating baseline (mean ± std) | Δmse model−baseline (mean ± std) | All-drugs val MSE | Paclitaxel gap (val − train) |
|---|---|---|---|---|
| `X_pca` | **207 ± 73** / 545 | **+0.00058 ± 0.00040** | 0.0106 ± 0.0008 | **+0.011 ± 0.020** |
| `X_scGPT` | **191 ± 94** / 545 | **+0.00072 ± 0.00047** | 0.0107 ± 0.0009 | **−0.002 ± 0.014** |

- **The continuous metric is the honest one — and it's negative news:** Δmse is **positive for both
  reps** (4 of 5 folds), i.e. on average the model is **marginally *worse* than the per-drug-mean
  constant**. The heads-beating count (~190–207 of 545, i.e. < 40% of heads) said the same thing all
  along; the continuous Δ just makes it unambiguous and stable (std ≪ the count's). PCA's Δ (+0.00058)
  is slightly *less bad* than scGPT's (+0.00072) — same direction as heads-beating.
- **The heads-beating count itself is *not* robust:** the fold std (±73–94) **dwarfs** the PCA−scGPT
  difference (~16). The single-split "169 vs 147" is within fold noise — don't read it as a real PCA
  advantage. (See *Metrics* above for why the count swings so hard.)
- **The overfitting direction survives, weakly:** mean paclitaxel gap is lower for scGPT (−0.002) than
  PCA (+0.011), consistent with the denoised-prior claim, but the spreads overlap.

### Better metric — per-drug correlation (27.06.2026)

Because viability clusters near 1.0, beating the per-drug-mean on MSE is a weak bar. §3 of the notebook
instead correlates **predicted vs true viability across held-out cell lines**, per drug (Spearman +
Pearson), restricted to the 461 drugs with real per-line variance (std ≥ 0.05, ≥ 5 val lines):

| Rep | mean Spearman | median Spearman | frac. drugs ρ > 0.3 |
|---|---|---|---|
| `X_pca` | −0.02 | −0.01 | 4.3% |
| `X_scGPT` | −0.05 | −0.05 | 3.9% |

- **Sobering:** per-drug rank correlation is **≈ 0 for both reps** — the models do **not** rank cell
  lines by drug response. The marginal MSE "wins" over the per-drug-mean reflect shrinking toward the
  constant, **not** real per-line predictive power. At this resolution (per-line viability broadcast to
  cells, values ≈ 1.0) the task is barely learnable beyond the mean — for *either* representation.
- This reframes the whole comparison: the scGPT-vs-PCA question is secondary to the fact that **neither
  rep yet predicts response variation across lines**. Motivates the better-target / better-metric work
  in [TODO.md](../TODO.md) (correlation-based selection, drugs with real variance).

> ⚠️ **Superseded (13.07.2026) — and partly an artifact.** Two independent problems with this verdict:
>
> 1. **It is an average over 545 drugs**, and averaging destroys it: on the 5 drugs that carry real
>    signal, the same model reaches Spearman **0.43–0.49** (next sections).
> 2. **The multi-task loss was unstandardized.** These runs used `mean_pv`, whose per-drug variance is
>    wildly heterogeneous, so a minority of wide-spread heads monopolized the shared trunk's gradient.
>    `notebooks/11_auc_vs_aucz.ipynb` **reproduces this failure on demand**: training K=545 on raw `auc`
>    (also unstandardized) drives per-drug Spearman to **−0.087 (scGPT) / +0.016 (PCA)**, while the
>    z-scored `auc_z` on the *same* drugs, model and split reaches **+0.430 / +0.378**.
>
> ⇒ The "neither rep ranks cell lines" conclusion was **never clean evidence about scGPT vs PCA**. It is
> substantially an artifact of a variance-dominated loss, and it does not survive per-drug
> standardization ([Step 03](03-model-and-training-design.md#measured-auc-vs-auc_z-notebooks11_auc_vs_aucz-ipynb-13072026)).

### Learnability-filtered subset — the signal was there all along (13.07.2026)

`notebooks/08_learnability_filter.ipynb` → `notebooks/09_learnable5_training.ipynb`. The §3 null result
above pooled a few learnable heads with hundreds of flat, inert ones. **Filter first, then ask.**

**The filter (`08`).** The learnability score of [`04`](../../notebooks/04_drug_coverage.ipynb)
(`resp_std × cov_frac`) is **degenerate on `auc_z`** — the target is z-scored per drug, so every drug
has std exactly 1.0 and all 545 tie. Spread is therefore measured on the **raw `auc` scale**, recovered
exactly via `uns["ctrp_score_scale"]`/`["ctrp_score_center"]` (and that `scale` vector *is* the per-drug
std of `auc`). The loose `04` gates (`cov ≥ 100 & std ≥ 0.05`) kept 439/545 and so never bit; the
missing condition is **differential response** — a drug must both **kill** a real population of lines
(`n_sens`: `auc ≤ 0.5`) and **leave one alive** (`n_res`: `auc ≥ 0.8`). A uniformly inert or uniformly
toxic drug has no cross-line ranking to learn, however well covered it is. **6 / 545 pass; the top 5 by
learnability are trained.**

**The result (`09`).** Both reps trained on those 5 heads (matched trunk, `--score auc_z`); the honest
metric is per-drug Spearman on **cross-validated out-of-fold predictions** — 5-fold GroupKFold over the
153 train+val lines, so every line is ranked by a model that never saw it (~150 lines per drug, versus
the 27 the fixed val split would allow):

| Rep | mean Spearman | mean Pearson | heads beating baseline | best val MSE |
|---|---|---|---|---|
| `X_pca` | **0.432** | 0.416 | 3 / 5 | 0.925 |
| `X_scGPT` | **0.488** | 0.482 | 4 / 5 | 0.777 |

(On `auc_z` the per-drug-mean null model scores MSE = **1.0** by construction, so these MSEs are
readable directly.) Per drug: `ml162` 0.59/**0.65**, `1s,3r-rsl-3` 0.58/**0.59**, `dasatinib`
0.52/**0.56**, `cay10618` **0.36**/0.35, `kx2-391` 0.11/**0.28** (PCA / scGPT).

- **Signal exists.** Against −0.02 / −0.05 over 545 drugs, the same architecture reaches ~0.45 here.
  **The 545-drug null result was a drug-selection artifact, not a representation failure.** The
  standing conclusion "the ceiling is the label, no gene representation can help" is *true on average
  and false on the drugs that matter* — **drug selection is a first-class lever**, and a cheap one.
- **The biology checks out.** The two strongest drugs are the **GPX4 inhibitors** (`ml162`,
  `1s,3r-rsl-3`): ferroptosis sensitivity tracks a cell's lipid-peroxidation/redox state, which is a
  *transcriptional* state. `dasatinib` (SRC/ABL) follows target addiction. The filter selected drugs
  whose variance has a transcriptional cause, not merely high-variance drugs.
- **First non-tie between the reps — and it survives a seed check (13.07.2026).** scGPT leads on every
  aggregate and on 4/5 drugs, most clearly where PCA collapses (`kx2-391`, 0.28 vs 0.11). Repeating the
  **K=545 `auc_z`** configuration over **3 seeds** (`notebooks/11`, `outputs/target/seed_stability.csv`):

  | seed | PCA | scGPT | gap |
  |---|---|---|---|
  | 42 | 0.388 | 0.430 | +0.043 |
  | 1 | 0.367 | 0.434 | +0.066 |
  | 7 | 0.355 | 0.472 | +0.117 |

  **Gap = +0.075 ± 0.038, sign-consistent across all three seeds.** No longer a one-seed accident — but
  3 seeds × 5 evaluation drugs is **consistent evidence, not a proven margin**. Do not upgrade it to a
  headline claim without more seeds and a wider drug set.
- **Ranking ≫ calibration.** `pred_std` is 0.53 (PCA) / 0.47 (scGPT) against a true spread of 1.0 — both
  models hedge toward each drug's mean. This is **not** an over-regularization artifact: `pred_std ≈ ρ ×
  true_std` is exactly what an MSE-optimal predictor must do (see the ablations below). Fine for ranking;
  to report in AUC units, divide by ρ.

**Is it the model? No — four knobs, all flat (`notebooks/10_diagnosis.ipynb`, 13.07.2026).**
Regularization (none → heavy), capacity (74,629 → 2,565 params), batch size (32/128/512) and sample
reweighting (line-balanced, focus-extremes) all leave out-of-fold Spearman within noise of the defaults
(PCA 0.41–0.44, scGPT 0.44–0.49). With regularization *off*, PCA memorizes the training lines (train MSE
≈ 0.01) and still reaches only 0.42 out-of-fold — the model is not being suppressed, it is out of signal.
Full table and reasoning in [Step 03](03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026).

**And the baseline that actually binds:** `RidgeCV` on the **150 cell-line mean embeddings** — no single
cells, no network — scores **0.428**, *tying* the PCA MLP (0.428) and within 0.06 of the scGPT MLP
(0.487). The whole deep single-cell apparatus currently buys **+0.06 Spearman, and only for scGPT**
(whose linear head drops to 0.438 — it *needs* the hidden layer; PCA does not). **Ridge on line means is
the baseline to beat from now on.** The cause is structural: the label is per cell line, so there are
~150 independent examples and the 34k cells are an illusion of sample size — which is why the remaining
levers are **label-side** (more lines, bulk pretraining, denoising), not model-side.

### The learnability gate measured the wrong quantity (27.07.2026)

The gate defined above — a drug must **kill** (`auc ≤ 0.5`) and **spare** (`auc ≥ 0.8`) a real population
— is not the right criterion for this project, and the 10-drug panel inherits the error.

**The mismatch.** `auc ≤ 0.5` asks *does the line die?*, i.e. it filters on absolute potency, which is
essentially the per-drug mean. But the target is `auc_z`, which **subtracts that mean**
([Step 03](03-model-and-training-design.md#known-problems-with-auc_z--the-scaling-is-not-yet-right-27072026)),
and the evaluation metric is **Spearman**, which only reads the *ordering* of lines. Whether a drug's
values sit around 0.4 or around 0.9 is irrelevant to both. The gate optimizes for a property the model
is neither given nor scored on.

**What it costs — `nutlin-3` as the clean example.** Raw-scale spread across the panel:

| drug | `auc_mean` | `auc_std` | kill (`≤0.5`) | gate verdict |
|---|---|---|---|---|
| `dasatinib` | 0.631 | **0.155** | 35 | selected |
| `nutlin-3` | 0.874 | **0.147** | **0** | rejected |

**Nutlin-3 has essentially the same spread as dasatinib** — its lines differ just as much, the whole
distribution simply sits higher (~0.6–1.0 instead of ~0.3–0.9). The reason is pharmacological, not
technical: nutlin-3 is **cytostatic**, not cytotoxic. p53 activation drives arrest and senescence, so
viability never falls below 50 % however sensitive the line is. Any threshold phrased as "does it kill"
is structurally blind to every cytostatic agent — which is a large fraction of targeted therapy.

This is not one unlucky drug: **116 of the 545** have zero kills yet `auc_std ≥ 0.10` and coverage
≥ 90 % (`oxaliplatin` among them). The gate discarded all of them silently.

**Why it matters beyond drug choice.** `nutlin-3`/TP53 is the single strongest association in the GDSC
pharmacogenomic screen, and it is *expression*-readable, not only genomic: `MDM2`, `CDKN1A` and
`RPL22L1` — p53 target genes — are selected in ~90–100 % of published gene sets predicting nutlin-3a
sensitivity. It is close to a best case for this model, and the filter threw it out.

**Correction to adopt:** replace kill/spare with **spread on the raw AUC scale** (`auc_std`, recoverable
exactly via `uns["ctrp_score_scale"]`, which *is* the per-drug std) plus coverage. This is the same
quantity that governs the noise-amplification problem in Step 03, so one criterion fixes both: high
`auc_std` = real signal to rank *and* a safe denominator for the z-score. Not yet re-run — the 10-drug
results above and the 8-drug panel below both still rest on the old gate.

### Literature-anchored drug panel — selecting without looking at our labels (25.07.2026)

**Why replace the filter.** The 10-drug set above is a *best-case* subset: the `08` gates were computed
on all 180 lines, val/test included, so selection saw held-out labels — the blocking limitation of the
15.07 progress report. Worse, the gates are not stable: shifting the kill/spare thresholds from
0.5/0.8 to 0.7/0.8 yields a completely **different** ten drugs of the same quality, so the filter
enriches reliably but *which* drugs it names is arbitrary. A panel anchored in **published** evidence
fixes the second problem outright — the drugs are named by citation, not by a threshold we chose — and
substantially reduces the first. It does **not** fully eliminate the first; see "Residual selection
effect" below, which is why the train-only check stays on the list.

**Selection rule.** Restrict `data/drug/all_sources_drug_catalog.csv` (`dataset == "CTRPv2"`) to single
agents with `compound_status ∈ {FDA, clinical}` — 173 of the 545, i.e. real therapeutics rather than
Broad screening probes — and keep those with an **independently published sensitivity determinant in
cancer cell lines**. The scientific rationale: where sensitivity is a documented, mechanistically
understood function of cell state, a transcriptome-based model *ought* to work; a failure there is a
model result, not a label artifact. This inverts the previous logic — the old filter asked "which drugs
have spread in *our* labels", which is unanswerable without peeking.

The catalog itself is built in `notebooks/02_compare_GDSC_CTRP.ipynb` from CTRP's official
`v20.meta.per_compound.txt`, with `gene_symbol_of_protein_target` → `target`,
`target_or_activity_of_compound` → `moa_or_pathway`, `cpd_status` → `compound_status`
([Step 01](01-datasets-and-harmonization.md)). Note this promotes the catalog from the "exploratory,
not consumed by any model" status recorded there to a **selection input**.

**Decision (25.07.2026) — the 8-drug panel.** `kill` = lines with raw `auc ≤ 0.5`, `spare` = `auc ≥ 0.8`,
`cov` = fraction of the 180 trainable lines, rank = position in
`outputs/learnability/ctrp_drug_learnability_auc.csv`:

| drug | target | kill / spare | cov | `08` rank | published determinant |
|---|---|---|---|---|---|
| `methotrexate` | DHFR | 52 / 31 | 0.94 | 6 | **`SLC19A1`** (reduced folate carrier) governs uptake; its loss is a classical resistance mechanism in cell lines — [Zhao & Goldman 2014](https://pubmed.ncbi.nlm.nih.gov/24396145/), [Wright et al., *Nature* 2022](https://www.nature.com/articles/s41586-022-05168-0) |
| `dasatinib` | SRC/ABL, EPHA2, KIT | 35 / 27 | 0.98 | 9 | six-gene **expression** model predicts sensitivity in 92 % of held-out breast and 83 % of lung lines — [Huang et al., *Cancer Res* 2007](https://aacrjournals.org/cancerres/article/67/5/2226/534297/Identification-of-Candidate-Molecular-Markers); `LYN` in lung ADC — [*Oncotarget* 2016](https://www.oncotarget.com/article/12657/text/) |
| `paclitaxel` | tubulin | 66 / 25 | 0.94 | 12 | **`ABCB1`** efflux + **`TUBB3`** — [*Oncotarget* 2016](https://www.oncotarget.com/article/9118/text/), [*Br J Cancer* 2016](https://www.nature.com/articles/bjc2016203) |
| `vincristine` | tubulin | 61 / 19 | 0.99 | 17 | same `ABCB1`/`TUBB3` axis (shared microtubule-disruptor resistance mechanism) |
| `afatinib` | EGFR, ERBB2 | 19 / 26 | 0.91 | 19 | `EGFR`+`ERBB2` co-amplification — [*Cancer Discov* 2019](https://aacrjournals.org/cancerdiscovery/article/9/2/199/10771/EGFR-and-MET-Amplifications-Determine-Response-to). ⚠ receptor *expression* alone did **not** correlate in pancreatic lines — [*Br J Cancer* 2011](https://www.nature.com/articles/bjc2011396) |
| `topotecan` | TOP1 | 37 / 18 | 0.97 | 21 | **`SLFN11`** expression, the canonical topoisomerase-inhibitor marker — Zoppoli et al., *PNAS* 2012 (NCI-60 + CCLE), pan-cancer replication [*PLOS One* 2019](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0224267), [review](https://www.sciencedirect.com/science/article/pii/S1359644625002922) |
| `tanespimycin` (17-AAG) | HSP90 | 14 / 44 | 0.96 | 30 | **`NQO1`** expression bioactivates the benzoquinone to its potent hydroquinone form; correlation confirmed in CCLE *and* GDSC across 7 cancer types — [*PLOS One* 2016](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0153181), [*Br J Cancer* 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4032580/) |
| `selumetinib` | MAP2K1/2 (MEK) | 12 / 80 | 0.97 | 43 | `BRAF` / `RAS` mutation — [*Mol Cancer Ther* 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2939826/) |

**All eight pass the `08` gate unchanged** (coverage ≥ 90 %, a real killed *and* a real spared
population), so this is a **re-ranking inside the gate-passing set, not a relaxation of it**. Only
`dasatinib` and `methotrexate` overlap the old ten — six of the eight are new drugs the previous filter
never named.

**Residual selection effect — what this panel does and does not fix.** As executed, the candidate list
was produced by ranking the 173 clinical/FDA compounds by `min(kill, spare)` and *then* applying the
literature criterion. Those counts come from our own `auc` values over all 180 lines, val/test included.
Concretely, several compounds with a published determinant dropped out because of **our label
statistics**, not because of the literature: `sirolimus` (6 kill / 63 spare), `neratinib` (12 / 75),
`clofarabine` (15 / 68), `cytarabine hydrochloride` (5 / 109), `gdc-0941` (5 / 45). So:

- **Fixed:** the arbitrariness of *which* drugs (they are now named by citation, and the choice is
  reproducible by someone who never sees our AUCs) and the threshold instability (0.5/0.8 → 0.7/0.8 no
  longer changes the panel).
- **Not fixed:** the panel is still enriched for drugs that happen to separate *our* 180 lines, so
  per-drug ρ measured on it retains an optimistic component. Honest description: **literature-anchored,
  spread-verified** — not label-blind.
- **Consequence:** the train-only-selection check ([TODO](../TODO.md)) stays **blocking** for any
  headline number computed on this panel.

**⚠️ The panel inherited the defect it was meant to escape (found 27.07.2026, after the run).** The
candidate list was ranked by `min(kill, spare)` before the literature criterion was applied — the very
kill-based quantity the section above shows to be the wrong one. The consequence is not marginal:
**32 of the 116 wrongly-discarded drugs are approved or in clinical trials**, among them `oxaliplatin`,
`bortezomib`, `ruxolitinib`, `regorafenib`, `entinostat` — and **`nutlin-3` itself**, the drug used to
demonstrate the defect (spread 0.147, coverage 0.96, status `clinical`, but zero kills, so balance 0 and
never a candidate).

So the eight are defensible as compounds, and every number computed on them stands, but the *pool they
were drawn from* was silently pre-filtered by the discredited criterion. Stating it plainly: **the panel
is literature-anchored, spread-verified, and drawn from a kill-filtered pool.**

**The fix, and why it is not applied yet.** Rebuild the pool on coverage plus `auc_std` — no kill counts
anywhere — then apply the literature criterion to *that*. It would very likely admit `nutlin-3` and
several of the other 32, giving a larger and better-justified panel. It also invalidates every number in
[Step 05](#step-1-executed--raw-auc--density-weighting-on-the-panel-notebooks14_panel_trainingipynb-27072026)
and requires re-running notebooks 13-15, so it is the next data step rather than a same-day correction.
 A cleaner variant — measure the kill/spare requirement on
  **GDSC2/PRISM** instead of on our CTRP labels — is recorded there as the follow-up.

**The panel is a hypothesis test, not a grab bag.** The determinants split by *data modality*, and our
input is expression only:

- **Expression-determined** — `methotrexate` (`SLC19A1`), `paclitaxel`/`vincristine` (`ABCB1`/`TUBB3`),
  `topotecan` (`SLFN11`), `tanespimycin` (`NQO1`), `dasatinib` (six-gene signature). The causal variable
  is *in* `X`, so these should be learnable.
- **Mutation-determined** — `selumetinib` (`BRAF`/`RAS` point mutations), `afatinib` (amplification;
  and expression explicitly failed to predict it in the pancreatic panel). The causal variable is
  **not** in `X` except through downstream expression, so weak per-drug ρ here is *expected* and is not
  evidence against the representation.

**Prediction to check:** ρ should be systematically higher in the first group. If it is, that is a
mechanistic validation of the whole approach; if `selumetinib`/`afatinib` also score well, the model is
picking up lineage rather than the stated mechanism and needs scrutiny.

**Caveat on the annotation gate.** CTRP's `target` column is **empty** for `paclitaxel` and
`vincristine` (and for `ml210`/`ml162`), which are annotated only via `moa_or_pathway`. A naive "keep
rows with a validated target" filter silently drops four of the most informative compounds — selection
must read `moa_or_pathway` as well. Not yet trained on this panel.

### Step 1 executed — raw AUC + density weighting on the panel (`notebooks/14_panel_training.ipynb`, 27.07.2026)

First run of the retired-`auc_z` setup: target raw `auc` winsorized at 1.1, the 8-drug literature panel,
per-sample inverse-density weights fitted **per fold on training lines only**, output layer excluded from
weight decay, head biases initialized to the train-fold per-drug means. Architecture, splits, optimizer
and batching unchanged, so the change is attributable. 5-fold GroupKFold over the 153 train+val lines,
**one seed (42)**.

| model | ρ `X_pca` | ρ `X_scGPT` | MSE `X_pca` | MSE `X_scGPT` |
|---|---|---|---|---|
| MLP, unweighted | 0.316 ± 0.003 | **0.377** | 0.0265 | 0.0254 |
| MLP, density-weighted | 0.308 | 0.369 | 0.0274 | 0.0254 |
| `RidgeCV` on line means | 0.306 | 0.299 | 0.0270 | 0.0268 |

Null (per-drug mean) MSE is 0.030, so the numbers are readable directly: the scGPT model explains ~15 %
of the variance **in AUC units**, RMSE ≈ 0.16 viability.

**Dispersion — the 5-fold numbers, which the table above omits.** Computed from the stored out-of-fold
predictions without retraining (`notebooks/15_diagnostics.ipynb` §5,
`outputs/diagnostics/result_dispersion.csv`):

| | pooled ρ | sd across the 5 folds | sd across the 8 drugs | per-drug range |
|---|---|---|---|---|
| PCA, unweighted | 0.315 | ±0.028 | 0.111 | 0.19 – 0.53 |
| scGPT, unweighted | 0.377 | ±0.043 | 0.091 | 0.30 – 0.55 |

Two things follow, and the second one qualifies the headline. The two dispersions answer different
questions — fold spread says how much the result depends on *which lines were held out*, drug spread says
how unevenly the model performs *across compounds* — and they are not interchangeable. And the
**scGPT−PCA gap of +0.062 is about the size of one fold standard deviation**, so it is consistent
evidence rather than an established margin, on top of being a single seed.

Note the estimators differ: pooling gives one correlation over ~150 lines, while the fold-wise mean
averages five correlations over ~30 lines each and sits slightly higher (0.341 / 0.387). Report the
pooled value as the point estimate and the fold spread as the dispersion; do not fuse them into a single
`mean ± sd` as though they came from one calculation.


**Confirmed — the June collapse was a K=545 effect, not a property of the target.** `auc_z` was adopted
because raw `auc` at 545 heads scored **−0.069** (scGPT). The same raw target on 8 comparable heads scores
**+0.377**. So the standardization was never fixing the *target*; it was compensating for pooling drugs
whose variances differ by 81×. Removing the cause (the panel) works at least as well as compensating for
it, without the side effect of amplifying noise-dominated drugs.

**Confirmed — the ridge tie for PCA, and the scGPT margin over it, both replicate.** PCA MLP 0.316 vs its
ridge 0.306 is the third independent panel on which averaging a line's cells into one vector loses
nothing. scGPT MLP 0.377 vs its ridge 0.299 is **+0.077**, against **+0.082** on the 14.07 10-drug panel
(0.402 vs 0.320, `ablations/ablation_capacity.csv`). This is a **replication on an independently chosen
drug set**, not a new result — which is the stronger claim of the two, since the earlier panel was
selected on our own labels and this one was not.

**Refuted — inverse-density loss weighting is not a lever here.** −0.006 (PCA) / −0.008 (scGPT) mean
Spearman; per drug a wash (`selumetinib` +0.09/+0.06, `tanespimycin` −0.06/−0.06). The mechanism did fire:
predicted spread rose 0.062 → 0.082 (PCA) and 0.062 → 0.080 (scGPT), the reduced shrinkage the method is
designed to produce. So the objective was reweighted as intended and the ranking did not follow. Coherent
with [`13`](../../notebooks/13_panel_distributions.ipynb): after winsorizing the artifacts above `auc` 1.1
every drug has |skew| ≤ 0.47, so there was little imbalance left for an imbalance correction to act on.
**Decision: do not carry the weighting into Step 2.** The pre-registered expectation ("MSE worse, Spearman
better") also failed — both stayed flat — which is what a null intervention looks like.

Two points worth stating explicitly, because they are what make this a usable result rather than a shrug:

- **The rising `pred_std` is the evidence that this is a null, not a bug.** A weighting that never reached
  the loss — wrong sign, weights dropped by the mask, an indexing slip — would leave the predictions
  numerically identical to the unweighted run. They are not: the model demonstrably hedges less. So the
  objective was changed as designed and the ranking still did not move. "Did not work" and "was broken"
  are different claims, and only the first is supported here.
- **It removes a candidate explanation for the shrinkage.** Predictions span 0.08 against a true spread of
  0.171, and one standing hypothesis was that the objective is dominated by the crowded middle of each
  drug's response range. That hypothesis is now tested and rejected: pointing the loss at the sparse
  extremes does not close the gap. What remains is that there is too little signal for ~150 independent
  cell lines to support — which is a label-side problem, not a loss-side one, and it is the direct
  argument for Step 2 and for more lines rather than for further objective engineering.

**Absolute numbers are lower than the 14.07 panel (0.356 / 0.402), and that is the expected direction.**
Those drugs were selected using all 180 lines including val/test; these were not chosen for spread on our
labels at all. Lower and more defensible is the trade that was made deliberately.

**Not addressed by this run.** Research question 2 — whether heterogeneity is learned implicitly — is
untouched and remains structurally untestable under a constant-within-line label. That is Step 2.

### Benchmarked with the real DrEval package (`notebooks/12_dreval_benchmark.ipynb`, 14.07.2026)

Not a re-implementation: `pip install drevalpy` (v1.5.1, <https://github.com/daisybio/drevalpy>), and our
data/model run through **their** `DrugResponseDataset`, **their** `split_dataset(mode="LCO")`, **their**
baselines and **their** `evaluate()`. `OncoMLP` is trained on the single cells of each fold's train lines
and predicts the held-out lines' cells, averaged back per cell line — scored on exactly the same pairs.

**LCO, 5-fold CV, the 5 learnable drugs, native `auc` units** (mean over folds; *normalized* = their
recipe: subtract the `NaiveMeanEffects` prediction from `y_true` **and** `y_pred`, then re-evaluate):

| Model | Spearman (raw) | **Spearman (norm.)** | **R² (norm.)** |
|---|---|---|---|
| `NaivePredictor` | 0.000 | 0.020 | −0.052 |
| `NaiveDrugMeanPredictor` | 0.197 | 0.000 | −0.002 |
| `NaiveCellLineMeanPredictor` | 0.000 | 0.020 | −0.052 |
| **`NaiveMeanEffectsPredictor`** | 0.197 | **0.000** | −0.002 |
| their `SingleDrugElasticNet` (scGPT) | 0.197 | **0.000** ❌ | −0.002 |
| their `SingleDrugRandomForest` (PCA) | 0.245 | 0.148 | 0.022 |
| their `SingleDrugElasticNet` (PCA) | 0.320 | 0.300 | 0.056 |
| their `SingleDrugRandomForest` (scGPT) | 0.468 | 0.438 | 0.178 |
| **`OncoMLP` (X_pca)** | 0.481 | **0.442 ± 0.071** | 0.178 |
| **`OncoMLP` (X_scGPT)** | **0.549** | **0.511 ± 0.085** | **0.224** |

1. **OncoMLP clears `NaiveMeanEffects` decisively** (normalized ρ = 0.511 vs 0.000) — the bar **half the
   published models in the DrEval paper fail**. Our normalized **R² = 0.224** is directly comparable to
   the numbers they report for their best models in LCO (**DIPK 11%**, **Random Forest 19%**).
2. **scGPT > PCA is confirmed externally** (+0.07 normalized, on *their* splits with *their* metrics) —
   an independent replication of the +0.075 ± 0.038 we measured ourselves.
3. **The single-cell MLP beats their line-level reference models on the same embeddings** (0.511 vs 0.438
   for `SingleDrugRandomForest` on scGPT; 0.442 vs 0.148 on PCA). This **qualifies the ridge result**
   above: against a *stronger* per-drug regressor on line-mean embeddings, the per-cell model does add
   something — small (+0.07) but consistent across both representations.

> **Note on LCO:** a held-out line is unseen, so `NaiveMeanEffectsPredictor` sets its cell-line effect to
> **0** and reduces to *global mean + drug effect* (hence it ties `NaiveDrugMean` here). Their normalized
> metric therefore removes the **drug** mean — the fix for the Simpson's-paradox artifact they describe.
> It does *not* remove the cell-line effect, because in LCO no honest predictor can know it.

⚠️ Still a **best-case subset**: the 5 drugs were selected by a filter that saw all 180 lines
(`notebooks/08`). DrEval fixes the *evaluation*, not our *selection*.

### Own-implementation check: what if the cell-line effect is also removed? (14.07.2026)

**Reference:** Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response
prediction models with DrEval*, **Nat. Commun. (2026)**. Their headline: *"deep learning models barely
outperform a naive model that predicts only the mean drug and cell line effects"* — about **half** of
published models fail to beat their `NaiveMeanEffectsPredictor`. Our setting **is** their **LCO**
(leave-cell-line-out) with per-drug evaluation, i.e. the split they recommend.

**The gap this exposed in our metric.** `auc_z` removes the **drug** mean but *not* the **cell-line**
mean. Some lines are simply sensitive to everything (σ of the line effect = **0.40**), so a model can
score a good per-drug correlation by learning *"this line is fragile"* — with zero drug-specific biology.
DrEval's normalized metric subtracts the mean-effects predictor from **prediction and truth**, then
correlates; what remains is **differential sensitivity only**. Mean effects are fit on **train lines
only**, inside each fold (`notebooks/outputs/dreval/dreval_normalized.csv`):

| | raw ρ | **normalized ρ** | naive baseline |
|---|---|---|---|
| K=5 · PCA | 0.427 | **0.368** | 0.291 |
| K=5 · scGPT | 0.488 | **0.396** | 0.291 |
| K=545 · PCA | 0.378 | **0.297** | 0.291 |
| K=545 · scGPT | 0.430 | **0.323** | 0.291 |

✅ **~80% of the effect survives** — it is genuine drug-specific signal, not the cell-line effect. For
scale, DrEval report their best models (DIPK, Random Forest) explaining **11–19%** of differential
sensitivity in LCO; ours is ρ² ≈ **0.16**, the same ballpark, on the filtered subset.

**Per drug (K=5, scGPT) — one of our five is an artifact:**

| drug | raw | **normalized** | naive baseline |
|---|---|---|---|
| `ml162` | 0.655 | **0.587** | 0.196 |
| `1s,3r-rsl-3` | 0.591 | **0.530** | 0.178 |
| `dasatinib` | 0.563 | **0.548** | 0.269 |
| `cay10618` | 0.347 | **0.306** | 0.226 |
| `kx2-391` | 0.283 | **0.006** ⚠️ | **0.584** |

⚠️ **`kx2-391` collapses to zero**: its entire apparent signal *was* the cell-line effect — exactly the
artifact class DrEval describes, found in our own results. The other four (both GPX4/ferroptosis
inducers, dasatinib, CAY10618) are real.

> **Decision (14.07.2026): adopt `NaiveMeanEffects` (drug mean + cell-line mean) as the standard
> baseline, and report raw *and* normalized correlations.** The per-drug-mean null is too weak; even
> ridge-on-line-means does not control for the cell-line effect. Any future claim must clear the
> normalized bar.

---

> **Decision — this is a best-case diagnostic, not a headline number.** The 5 drugs were selected using
> **all 180 lines, val/test included**, so the selection saw the held-out labels. It answers the
> question it was built to answer ("does *any* cross-line signal exist here?" — yes). Turning it into a
> reportable result requires selecting the subset **train-only inside each CV fold** and repeating over
> seeds; that is now the top [TODO](../TODO.md) item.

### Which change actually produced the gain? (decomposition, 13.07.2026)

Three things changed at once between the §3 null result and the 5-drug result: the **target**
(`mean_pv` → `auc_z`), the **training scope** (K=545 → K=5) and the **measurement** (27 fixed-val lines
→ ~150 out-of-fold lines). A Spearman over 27 points has a standard error of roughly ±0.2, so the
measurement change alone could have manufactured the effect. Isolating them — mean per-drug Spearman
**on the same 5 drugs** throughout:

| Config | PCA | scGPT |
|---|---|---|
| **Old** — K=545, `mean_pv`, 27 val lines (`07` §3) | −0.036 | **−0.286** |
| K=545, **`auc_z`**, 27 val lines | +0.254 | **+0.349** |
| K=545, `auc_z`, **150 OOF lines** | +0.378 | **+0.430** |
| **New** — **K=5**, `auc_z`, 150 OOF lines (`09`) | +0.434 | **+0.488** |

Read down the column — each row changes exactly one thing:

- **The target switch is the dominant term: +0.29 (PCA) and +0.64 (scGPT)**, with head count and the
  27-line measurement held fixed. This is a genuine improvement in the *predictions*, not in the
  metric: on `mean_pv` the model's ranking of these drugs was **negative** (scGPT −0.29 — worse than a
  coin flip), and on `auc_z` the same model, same drugs, same 27 lines ranks them at +0.35.
- **Honest measurement adds ~+0.1** (27 → 150 held-out lines). Real, but it is a *precision* gain, not
  a model gain — and it is why the old −0.47-type numbers should never have been read as findings.
- **Drug filtering adds only ~+0.06.** The learnability filter is the *smallest* of the three effects.
  Its value is that it identifies *where* the signal lives — at K=545 these drugs already reach 0.430
  (scGPT), so the filter is not what created the signal.

> ⚠️ **Which part of the target switch did the work? The z-scoring — *not* the curve fit.** An earlier
> version of this section credited the curve-fit AUC ("the dose-averaged viability target was destroying
> the signal the curve fit preserves"). `notebooks/11` **falsifies that**: trained head-to-head, the
> legacy `mean_pv` and the raw curve-fit `auc` behave *identically* — at K=5 they tie
> (0.450 / 0.481 vs 0.439 / 0.482, CIs fully overlapping), and at K=545 they **both collapse**
> (+0.027 / −0.070 vs +0.016 / −0.087). The curve fit is worth keeping for principled reasons (post-QC
> fitting, and the same metric family GDSC2 reports — [Step 06](06-cross-database-integration.md)), but
> it buys **no measurable accuracy**. **Every bit of the +0.64 comes from standardizing the per-drug
> variance**, i.e. from fixing the multi-task loss.

⇒ **The model is genuinely better than it was this morning, and the credit goes to one thing: per-drug
standardization of a 545-head shared loss.**

### Gene-set sweep — heads-beating vs gene count (incl. all_genes, 28.06.2026)

Does either rep have a preferred gene-set size? `notebooks/07_training.ipynb` §4 builds each variant
(1k/2k/3k/5k **plus `all_genes`**, full pipeline incl. scGPT re-embed; `05` §B) and runs the same
**5-fold GroupKFold, test held out, all 545 drugs** — so the HVG-vs-all-genes comparison is
apples-to-apples under identical CV:

| Gene set | genes | `X_pca` heads-beat | `X_scGPT` heads-beat | Δmse (PCA / scGPT) |
|---|---|---|---|---|
| `hvg1000` | 1,000 | 207 ± 75 | 193 ± 83 | +0.00058 / +0.00060 |
| `hvg2000` | 2,000 | 203 ± 78 | 185 ± 84 | +0.00062 / +0.00064 |
| `hvg3000` | 3,000 | 216 ± 85 | 190 ± 83 | +0.00053 / +0.00063 |
| `hvg5000` | 5,000 | 210 ± 73 | 189 ± 94 | +0.00055 / +0.00074 |
| `all_genes` | 22,722 | 204 ± 86 | 184 ± 90 | +0.00058 / +0.00069 |

- **No sweet spot, and no all-genes advantage.** Both reps are **flat across the whole axis** (PCA
  ~203–216, scGPT ~184–193) — filtering does not help scGPT (contrary to the earlier hunch), and
  `all_genes` is **no better than HVG** for either rep (PCA's `all_genes` 204 sits mid-band, below
  hvg3000's 216; the earlier "PCA prefers all genes" is not reproduced). Val MSE ~constant
  (0.0105–0.0107) throughout.
- PCA is marginally higher than scGPT at every gene count, but the ±73–94 fold spread overlaps
  completely — within noise at all sizes, consistent with the CV finding above.
- **Δmse > 0 at every gene-set size** for both reps: the model stays marginally *worse* than the
  per-drug-mean baseline regardless of how many genes feed it.

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
