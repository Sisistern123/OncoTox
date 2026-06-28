# OncoTox — TODO

Action list. Scientific narrative + full numbers live in
[project_progress.md](./project_progress.md) and [`docs/steps/`](./steps/); this is the running tasks.

## Done

- [x] **8-run matrix** `{hvg5000, all_genes} × {X_pca, X_scGPT} × {single-paclitaxel, all-drugs K=545}`,
      shared cell-line-grouped `split_ctrp`, matched `(128,64)` trunk → [Step 05](./steps/05-multitask-results.md).
- [x] **Matched input dim**: PCA raised to **512 components** (`add_pca.DEFAULT_N_COMPS`, `--pca-n-comps`)
      so PCA and scGPT share input width; full matrix re-run at 512-d (supersedes the ~50-d one).
- [x] **5-fold GroupKFold CV** (test held out, 153 train+val lines) → `07` §2. **Difference not robust:**
      heads-beating `hvg5000` PCA **207 ± 73** vs scGPT **191 ± 94** (fold std ≫ the ~16 gap); Δmse > 0 for
      both (model marginally *worse* than the per-drug-mean). scGPT's only edge = slightly lower overfitting.
- [x] **Per-drug correlation** (Spearman/Pearson, pred vs true across lines, 461 real-variance drugs) → `07` §3.
      **≈ 0 for both** (mean Spearman PCA −0.02, scGPT −0.05; ~4% ρ > 0.3) — neither ranks cell lines.
- [x] **Gene-set sweep** 1k/2k/3k/5k **+ all_genes** under CV → `07` §4, variants built in `05` §B.
      **No sweet spot, no all-genes advantage** — flat across the axis, within noise.
- [x] **Target distribution** (data only) → `04` (`outputs/target_distribution.png`): viability clusters
      near 1.0 (median 0.91; 75% ≥ 0.8); per-drug std median 0.088, only 3% truly flat; a loose
      cov ≥ 100 & std ≥ 0.05 filter keeps **439/545** → coverage+std alone removes few.
- [x] **Coverage & learnability** analysis → `04` (`outputs/*_drug_learnability.csv`, `drug_coverage.png`).
- [x] **Cancer-type UMAPs** → `06` §8: 2-panel PCA-vs-scGPT (`outputs/umap_cancertype_pca_vs_scgpt.png`,
      dpi 300) + full gene-set sweep grid (`outputs/umap_sweep_cancertype.png`, dpi 200). Tissue islands
      (PCA) vs continuous manifold (scGPT) at every gene count.
- [x] **Initial informative-drug list** (CTRPv2) from `04` shared with Hashimoto-san (known not-final;
      GDSC version was for her only, not the modelling work).
- [x] **190 vs 180 resolved**: 190 = CTRPv2 roster name-matches, 180 = lines with post-QC measurements.

**Net read:** with a fair 512-d + cross-validated comparison **PCA ≈ scGPT**, and the model learns the
per-drug mean, not cross-line sensitivity. The ceiling is the **label** (bulk value broadcast to a
line's cells → ~126 independent lines; viability compressed near 1.0), not the gene representation.

## Next focus — drug-learnability filtering (28.06.2026)

Filter drugs by learnability *before* training/eval, and pin down whether **any** real signal exists.
- [ ] **Define learnability properly** — stricter than coverage+std (that keeps 439/545). Candidate:
      coverage ≥ threshold **AND genuine differential response** (e.g. ≥ N lines with viability < 0.7),
      optionally CV per-drug correlation > 0. Drop near-flat / low-coverage drugs.
- [ ] **Best-case diagnostic** — per-drug correlation + predicted-vs-true scatter on the *most-responsive,
      well-covered* drugs. Decides whether the ceiling is the data or the model.
- [ ] **Predicted-vs-true diagnostic** — is the model just predicting the per-cell-line mean? Does
      averaging per-cell predictions back to the line help?
- [ ] **Filter → re-run** matrix / CV / sweep on the learnable subset; do PCA-vs-scGPT and the absolute
      metrics separate once dead heads are removed?
- [ ] *(Stretch)* cluster cell lines by response and **stratify train/val/test** (high/med/low) for
      lower-variance evaluation.

## Levers / later

- [ ] **Bulk RNA-seq pretraining / scDEAL-style denoising + domain adaptation** — attacks the
      noisy-label bottleneck (the real ceiling).
- [ ] **Cross-database PRISM** (masked multi-task) — [Step 06](./steps/06-cross-database-integration.md).
      (GDSC not a modelling priority; was only for Hashimoto-san's list.)
- [ ] **XAI** — feature importance → resistance drivers — [Step 07](./steps/07-xai-feature-interpretability.md).
- [ ] Confirm scGPT input preprocessing in `gen_embeds.py` (raw counts vs CPM) so scGPT isn't handicapped.
- [ ] *(Optional)* regenerate scGPT embeddings from scratch (reproducibility pass; identical output).
- [ ] *(Optional)* re-run `split_paclitaxel` single-task to fill [Step 04](./steps/04-single-task-results.md)'s
      PCA column, or retire that progression.

## Roadmap (project plan)

- [ ] Cross-database integration — PRISM then GDSC, efficacy + toxicity ([Step 06](./steps/06-cross-database-integration.md)).
- [ ] XAI / feature interpretability ([Step 07](./steps/07-xai-feature-interpretability.md)).
- [ ] Foundation model + clinical fine-tuning ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)).
