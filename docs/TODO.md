# OncoTox — TODO

My running task list. Scientific open questions live in
[project_progress.md](./project_progress.md#open-questions-carried-forward); this is the action list.

## Done

- [x] **8-run matrix trained** (13.06): `{all_genes, hvg5000}` × `{X_pca, X_scGPT}` ×
      `{single paclitaxel, all-drugs K=545}`, all on the shared `split_ctrp`
      ([Step 05](./steps/05-multitask-results.md)).
- [x] **Matched trunk** (14.06): both reps now use `(128,64)` so head capacity is equal — only the
      input representation differs. Re-ran the matrix. Result: scGPT's win is *lower overfitting*, not
      higher accuracy; with matched capacity PCA is competitive/better on all-drugs
      (`hvg5000` 158 vs 156; `all_genes` PCA 196 vs scGPT 137).
- [x] Verified both variants' gene counts, cell alignment, `.X` (CPM), and `X_pca`/`X_scGPT` dims;
      embeddings confirmed on filtered (4,576) and non-filtered (20,570) data.

## Next steps (from 15.06 progress report)

**Guiding question:** is scGPT's lower-overfitting edge a *real, useful* advantage, or is PCA
sufficient for this task? The steps below are designed to answer it.

**1. Is the difference real?**
- [x] Grouped **5-fold cross-validation** (`GroupKFold` over cell lines, **test held out**, CV over the
      153 train+val lines) → `notebooks/07_training.ipynb` §2. **Answer: largely no.** Heads-beating
      `hvg5000` PCA **207 ± 73** vs scGPT **191 ± 94** — the fold std dwarfs the ~16 difference, so the
      single-split "169 vs 147" is within noise. Overfitting direction survives weakly (paclitaxel gap
      PCA +0.011 ± 0.020 vs scGPT −0.002 ± 0.014). See [Step 05](./steps/05-multitask-results.md).

**2. Fair comparison & better metric**
- [x] **Match input dimensionality**: PCA now uses **512 components** (`add_pca.DEFAULT_N_COMPS`,
      overridable with `--pca-n-comps`) so PCA and scGPT share the same input width. The **full 8-run
      matrix was re-run at 512-d** in `notebooks/07_training.ipynb` (run dirs `runs/20260627_1913xx_*`),
      superseding the ~50-d matrix; results in [Step 05](./steps/05-multitask-results.md).
- [x] Add a **per-drug correlation** metric (Spearman/Pearson, predicted vs true across held-out
      lines), restricted to drugs with real response variance → `notebooks/07_training.ipynb` §3.
      **Sobering result:** per-drug rank correlation ≈ 0 for *both* reps (mean Spearman PCA −0.02,
      scGPT −0.05; ~4% of 461 drugs ρ > 0.3) — neither rep ranks cell lines; MSE wins are shrinkage to
      the mean, not real predictive power. See [Step 05](./steps/05-multitask-results.md).
- [x] **HVG-count sweep** (1k/2k/3k/5k) under CV (all drugs, test held out) → `07` §4, variants built in
      `05` §B. **No sweet spot:** scGPT heads-beating is flat (~185–193) across 1k–5k (filtering doesn't
      help it); PCA marginally higher everywhere but within fold noise. See [Step 05](./steps/05-multitask-results.md).

**3. Understand the result**
- [x] Per-drug **coverage & response-distribution** analysis → `notebooks/04_drug_coverage.ipynb`
      (coverage per drug, variance per drug, and coverage/variance vs beats-baseline). Finding:
      no drug covers all 180 lines (max 179, median 171); 80 drugs < 50% coverage; the low-coverage
      drugs (≈16 lines, n_val 221) are the unreliable/hardest heads.
- [ ] **Predicted-vs-true** diagnostics: is the model just predicting the per-cell-line mean? Does
      averaging per-cell predictions back to the line help? Compare single-task paclitaxel (K=1) vs
      the paclitaxel head inside the K=545 run.
- [ ] Confirm scGPT input preprocessing in `gen_embeds.py` (raw counts vs CPM) so scGPT isn't handicapped.

**4. Levers most likely to move the needle**
- [ ] **Bulk RNA-seq pretraining / scDEAL-style denoising + domain adaptation** — attacks the
      noisy-label bottleneck (the real ceiling, and where scGPT-style representations should pay off).
- [ ] Cross-database **PRISM + GDSC** (masked multi-task) — see Roadmap / [Step 06](./steps/06-cross-database-integration.md).
- [ ] **XAI** — feature importance → resistance drivers — see Roadmap / [Step 07](./steps/07-xai-feature-interpretability.md).

## Housekeeping

- [ ] **Presentation slides** (15.06) — outstanding fixes: add a Core-Hypothesis slide before the
      first UMAP; de-duplicate the two UMAP slides; Results slide needs "heads beating baseline
      (out of 545)" + the conclusion line; fix "5.000" → "5,000"; reword target "(AUC)".
- [ ] Run the **UMAP cells** of `notebooks/06_verify_variants.ipynb` (compute-heavy; not yet run).
- [ ] *(Optional)* Regenerate scGPT embeddings from scratch — identical output; only for a
      reproducibility pass.
- [ ] *(Optional)* Re-run the `split_paclitaxel` single-task to fill
      [Step 04](./steps/04-single-task-results.md)'s PCA column, or retire that progression.

## Roadmap (project plan)

- [ ] **Cross-database integration** — PRISM then GDSC, efficacy + toxicity heads
      ([Step 06](./steps/06-cross-database-integration.md)).
- [ ] **XAI / feature interpretability** ([Step 07](./steps/07-xai-feature-interpretability.md)).
- [ ] **Foundation model + clinical fine-tuning** ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)).
- [x] **190 vs 180 resolved** (14.06): not normalization — 190 = CTRPv2 roster name-matches, 180 =
      lines with actual post-QC measurements (10 listed-but-unscreened drop out). Use 180.
