# OncoTox — TODO

My running task list. Scientific open questions live in
[project_progress.md](./project_progress.md#open-questions-carried-forward); this is the action list.

## Next up (rerun day — 2026-06-13)

- [ ] **Rerun the whole pipeline from scratch, both variants** (`hvg5000`, `all_genes`):
      `convert → scGPT embeddings → targets → splits → pca`. **Includes regenerating the scGPT
      embeddings.**
- [ ] **Re-train all baselines — the 8-run matrix:** `{all_genes, hvg5000}` × `{X_pca, X_scGPT}` ×
      `{single-task paclitaxel, all-drugs K=545}`. Refresh the tables in
      [Step 04](./steps/04-single-task-results.md) / [Step 05](./steps/05-multitask-results.md).
      Per condition, PCA uses the full filtered set (5,000 / 22,722) and scGPT its in-vocab subset
      (4,576 / 20,570) — **the OOV gap is intentionally part of the comparison** (scGPT's vocabulary
      coverage is a property of the model; *shouganai*).
- [ ] Sanity-check the rerun with `notebooks/verify_variants.ipynb` (gene counts, cell alignment,
      `.X` state, PCA-vs-scGPT UMAPs).
- [ ] **Presentation slides.**

## Make the PCA-vs-scGPT comparison properly fair

- [ ] Address the dimensionality / capacity mismatch: PCA `X_pca` is ≈50-d with hidden `(64,32)`,
      scGPT is 512-d with hidden `(128,64)`. Match input dim and/or head capacity so the comparison
      measures representation quality, not capacity (e.g. sweep PCA `n_comps`, or give both reps the
      same trunk).
- [ ] Confirm the scGPT input preprocessing in `gen_embeds.py` (raw counts vs CPM) so scGPT is fed
      what it expects and isn't handicapped.
- [ ] Add multiple seeds / CV folds for error bars (only 27 val cell lines → noisy point estimates).
- [ ] Single-task re-run on `split_ctrp` for an apples-to-apples "does multi-task help paclitaxel?".

## Roadmap (project plan)

- [ ] **Cross-database integration** — PRISM then GDSC, efficacy + toxicity heads
      ([Step 06](./steps/06-cross-database-integration.md)).
- [ ] **XAI / feature interpretability** ([Step 07](./steps/07-xai-feature-interpretability.md)).
- [ ] **Foundation model + clinical fine-tuning** ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)).
- [ ] Evaluate HVG-5000 vs `all_genes` (`notebooks/hvg_vs_all_genes_umap.ipynb`).
- [ ] Reconcile the 190 (audit) vs 180 (pipeline) cell-line overlap and pick one for the thesis.
