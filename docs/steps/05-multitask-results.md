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
  Reported as CV mean ± std; the per-fold `cv_folds.csv` also carries `median_delta` and `frac_beat`
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
