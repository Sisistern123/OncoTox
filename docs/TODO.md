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
- [ ] Grouped **5-fold cross-validation** (`GroupKFold` over cell lines); report **mean ± std** for
      the overfitting gap and heads-beating (turns "158 vs 156" into "155 ± 9 vs 152 ± 11"). Only
      27 val lines in one split → current numbers are point estimates.

**2. Fair comparison & better metric**
- [ ] **Match input dimensionality**: sweep PCA `n_comps` up to 512 so only the representation
      differs (currently hard-coded to scanpy's default 50; max = `min(n_cells, n_genes)−1`).
- [ ] Add a **per-drug correlation** metric (Spearman/Pearson, predicted vs true across held-out
      lines), restricted to drugs with real response variance — more informative than beating a constant.
- [ ] **HVG-count sweep** (1k/2k/3k/5k) under CV — scGPT improved with filtering while PCA preferred
      all genes; find scGPT's sweet spot.

**3. Understand the result**
- [ ] Per-drug **response-distribution** analysis (how many drugs are learnable, given viability ≈ 1.0?).
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
- [ ] Run the **UMAP cells** of `notebooks/verify_variants.ipynb` (compute-heavy; not yet run).
- [ ] *(Optional)* Regenerate scGPT embeddings from scratch — identical output; only for a
      reproducibility pass.
- [ ] *(Optional)* Re-run the `split_paclitaxel` single-task to fill
      [Step 04](./steps/04-single-task-results.md)'s PCA column, or retire that progression.

## Roadmap (project plan)

- [ ] **Cross-database integration** — PRISM then GDSC, efficacy + toxicity heads
      ([Step 06](./steps/06-cross-database-integration.md)).
- [ ] **XAI / feature interpretability** ([Step 07](./steps/07-xai-feature-interpretability.md)).
- [ ] **Foundation model + clinical fine-tuning** ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)).
- [ ] Reconcile the 190 (audit) vs 180 (pipeline) cell-line overlap and pick one for the thesis.
