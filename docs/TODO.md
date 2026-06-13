# OncoTox — TODO

My running task list. Scientific open questions live in
[project_progress.md](./project_progress.md#open-questions-carried-forward); this is the action list.

## Done — 2026-06-13

- [x] **8-run matrix trained** on the fixed `X_pca`: `{all_genes, hvg5000}` × `{X_pca, X_scGPT}` ×
      `{single paclitaxel, all-drugs K=545}`, all on the shared `split_ctrp`. Results in
      [Step 05](./steps/05-multitask-results.md). Headline: all-drugs heads-beating **scGPT > PCA in
      both gene sets** (135/103, 141/80); single-paclitaxel **PCA > scGPT**.
- [x] Verified both variants' gene counts, cell alignment, `.X` (CPM), and `X_pca`/`X_scGPT` dims
      (the `verify_variants.ipynb` non-UMAP checks). Embeddings confirmed on filtered (4,576) and
      non-filtered (20,570) data.

## Still to do

- [ ] **Presentation slides.**
- [ ] Run the **UMAP cells** of `notebooks/verify_variants.ipynb` for the visual PCA-vs-scGPT
      comparison (compute-heavy; not run on 13.06).
- [ ] *(Skipped 13.06, optional)* Regenerate scGPT embeddings from scratch — the existing embeddings
      are verified correct, so re-embedding both variants (~hours, 32 GB, identical output) was a
      no-op; do it only if a from-scratch reproducibility pass is wanted.
- [ ] *(Optional)* Re-run the `split_paclitaxel` single-task to fill [Step 04](./steps/04-single-task-results.md)'s
      PCA column, or retire that progression in favour of the matrix's `split_ctrp` single-task.

## Make the PCA-vs-scGPT comparison properly fair

- [ ] Address the dimensionality / capacity mismatch: PCA `X_pca` is ≈50-d with hidden `(64,32)`,
      scGPT is 512-d with hidden `(128,64)`. Match input dim and/or head capacity so the comparison
      measures representation quality, not capacity (e.g. sweep PCA `n_comps`, or give both reps the
      same trunk).
- [ ] Confirm the scGPT input preprocessing in `gen_embeds.py` (raw counts vs CPM) so scGPT is fed
      what it expects and isn't handicapped.
- [ ] Add multiple seeds / CV folds for error bars (only 27 val cell lines → noisy point estimates).
- [ ] Compare the single-task paclitaxel (K=1, now on `split_ctrp`) against the paclitaxel **head**
      inside the K=545 run — "does multi-task help paclitaxel?" (the single-task numbers now exist).

## Roadmap (project plan)

- [ ] **Cross-database integration** — PRISM then GDSC, efficacy + toxicity heads
      ([Step 06](./steps/06-cross-database-integration.md)).
- [ ] **XAI / feature interpretability** ([Step 07](./steps/07-xai-feature-interpretability.md)).
- [ ] **Foundation model + clinical fine-tuning** ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)).
- [ ] Evaluate HVG-5000 vs `all_genes` (`notebooks/verify_variants.ipynb`).
- [ ] Reconcile the 190 (audit) vs 180 (pipeline) cell-line overlap and pick one for the thesis.
