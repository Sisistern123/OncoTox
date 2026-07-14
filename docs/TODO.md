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
- [x] **Target score → `auc_z`** (13.07.2026): curve-fit `area_under_curve / conc_pts_fit`, z-scored per
      drug; `--score {auc_z,auc,mean_pv}` on every script, one targets h5ad per score →
      [Step 03](./steps/03-model-and-training-design.md). Within-drug Spearman vs the old `mean_pv` is
      only **0.72** (median), so this is *not* cosmetic — but nothing is re-trained on it yet.

- [x] **Learnability filter + best-case diagnostic** (13.07.2026) → `08` + `09`. Harsh gates (coverage,
      AUC spread, **and** a real killed *and* surviving population) keep **5/545**. Trained on those:
      out-of-fold per-drug Spearman **0.43 (PCA) / 0.49 (scGPT)** vs ≈ 0 over all 545 →
      [Step 05](./steps/05-multitask-results.md).

**Net read (revised 13.07.2026):** the old net read — *"the ceiling is the label; the model learns the
per-drug mean, not cross-line sensitivity; PCA ≈ scGPT"* — **does not survive today's work.** Two separate
defects produced it, and the **target was the bigger one**:

1. **An unstandardized multi-task loss.** With 545 heads of wildly different variance, a few wide-spread
   drugs monopolized the shared trunk. `11` reproduces the failure on demand (raw `auc` at K=545 → scGPT
   ρ = **−0.087**) and shows `auc_z` fixes it (**+0.430**, same drugs/model/split).
2. **A drug set that was mostly unlearnable**, which dragged the *average* to zero even where signal
   existed.

Decomposed on the 5 learnable drugs (`09`), the gain from ≈0 → 0.49 splits as: **target ≈ +0.64 (scGPT)**,
honest 150-line measurement ≈ +0.1, drug filtering ≈ +0.06. **The target change is the single biggest
improvement this project has made.** The label ceiling is still real (~150 independent lines; ridge on
line means ties the MLP — see below), but "no gene representation can help" was never established.

## Next focus — make the 5-drug result honest (13.07.2026)

The `08`/`09` numbers are a **best-case diagnostic**: the 5 drugs were selected using all 180 lines,
val/test included, so the selection saw held-out labels. Turning it into a reportable result:

- [ ] **Train-only selection** — run the learnability gates *inside each CV fold* (train lines only) and
      re-measure. If the effect survives, it is real; this is the blocking item.
- [x] **Multiple seeds** (13.07.2026) — done for K=545 `auc_z`: gap **+0.075 ± 0.038**, sign-consistent
      over 3 seeds. Still thin; see the "more seeds + wider drug set" item above.
- [ ] **Loosen to ~20–50 drugs** — 5 is a diagnostic, not a model. Where does the signal die as the gates
      relax? (`ctrp_drug_learnability_auc.csv` is already ranked for this.)
- [ ] **Re-run the full 8-run matrix + CV on `--score auc_z`** for a like-for-like against the `mean_pv`
      Steps 04–05 numbers. **Expect this to overturn them:** at K=545 the old unstandardized target was
      destroying the signal (below), so the 8-run matrix's conclusions are suspect, not just stale.
- [x] **All 3 targets measured head-to-head** (13.07.2026) → `11`, with bootstrap CIs, Pearson, and per-drug
      dots. **z-scoring is essential at K=545** (`mean_pv` −0.070 / `auc` −0.087 / **`auc_z` +0.430**,
      scGPT) and **irrelevant at K=5** (all three tie, ρ ≈ 0.42–0.49). ⚠️ **The curve fit buys no
      accuracy** — `mean_pv` ≈ `auc` everywhere; keep it for GDSC comparability, not performance. Keep
      `auc_z` as the default; `--score auc` is fine on a spread-homogeneous subset (native units).
- [x] **Seed check on scGPT vs PCA** (13.07.2026) → `11`. K=545 `auc_z` gap **+0.075 ± 0.038**, sign-
      consistent over seeds 42/1/7. Consistent evidence, **not** a proven margin.
- [ ] **More seeds + a wider drug set** before scGPT > PCA becomes a headline claim (3 seeds, 5 drugs is
      thin). Pair it with the train-only selection below.
- [x] ~~Fix the calibration shrinkage with lighter regularization~~ — **withdrawn 13.07.2026.** `pred_std ≈
      ρ × true_std` is what an MSE-optimal predictor *must* do; the shrinkage is correct calibration, not
      over-regularization, and loosening dropout raises MSE (`10` §1). To report in AUC units, divide by ρ.

## Model-side tuning is closed (13.07.2026)

`notebooks/10_ablations.ipynb`: regularization (none→heavy), capacity (74,629→2,565 params), batch size
(32/128/512) and reweighting (line-balanced, focus-extremes) **all leave Spearman flat** (PCA 0.41–0.44,
scGPT 0.44–0.49); the current defaults are already at/near the best on every axis. With regularization
off, PCA memorizes the train lines (train MSE ≈ 0.01) and still only reaches 0.42 out-of-fold — it is out
of *signal*, not out of *capacity*. **Don't spend more time on architecture or hyperparameters.**

**New baseline to beat: `RidgeCV` on the 150 cell-line mean embeddings** (ρ = **0.428**) — it *ties* the
PCA MLP and comes within 0.06 of the scGPT MLP, with no single cells and no network. The per-drug-mean
null is too weak a bar; report ridge alongside it from now on. scGPT + hidden layer (0.487) is the only
configuration that clears it — and scGPT's *linear* head drops to 0.438, so it genuinely needs the
nonlinearity while PCA does not.

- [ ] **Make the single-cell dimension earn itself** — averaging a line's cells into one vector currently
      loses nothing. Test MIL / attention pooling over a line's cells (predict the line label from a *bag*
      of cells), which at least matches the true label resolution. If that doesn't beat ridge either, the
      per-cell framing needs a different justification.
- [ ] **Add ridge (line-level) to `07`'s comparison tables** so every future claim is scored against it.
- [ ] *(Optional)* **z-score train-only.** The per-drug mean/std currently use all 180 lines, val/test
      included — mild leakage. Fixing it means computing splits before the targets step.
- [ ] *(Stretch)* cluster cell lines by response and **stratify train/val/test** (high/med/low) for
      lower-variance evaluation.

## DrEval alignment (14.07.2026)

Paper: Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response
prediction models with DrEval*, Nat. Commun. 2026. Half of published models don't beat a naive
drug-mean + cell-line-mean predictor. **Our ridge ≈ MLP finding is the field's norm, independently
reproduced.** Our split *is* their LCO.

- [x] **Normalized evaluation** → `outputs/dreval_normalized.csv`. After subtracting mean effects from
      prediction *and* truth: **scGPT K=5 ρ = 0.396** (raw 0.488), naive baseline 0.291 → **~80% of the
      signal is genuine differential sensitivity.** ⚠️ **`kx2-391` collapses to 0.006** — its signal was
      entirely the cell-line effect. The other 4 drugs survive.
- [ ] **Make `NaiveMeanEffects` the default baseline** in `train_multitask.py` (currently: per-drug mean,
      too weak — it does not control for the cell-line effect at all).
- [ ] **Report raw + normalized** correlations everywhere from now on.
- [ ] *(Consider)* their other splits — LTO (leave-tissue-out) and LDO (leave-drug-out). LDO would test
      whether anything generalizes across chemical space; DrEval found **no model** beats naive there.
- [ ] *(Consider)* **CurveCurator** for standardized dose-response fitting, as they recommend.

## Levers / later

- [ ] **Bulk RNA-seq pretraining / scDEAL-style denoising + domain adaptation** — attacks the
      noisy-label bottleneck (the real ceiling). **Promoted by the `10` ablations:** with model-side
      tuning closed and ridge-on-150-lines matching the MLP, the only remaining levers are label-side —
      above all **more independent cell lines** (SCP542×CTRPv2 caps at 180; CTRPv2 itself has ~1,100).
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
