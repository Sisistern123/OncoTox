# OncoTox — TODO

Action list. Scientific narrative + full numbers live in
[project_progress.md](./project_progress.md) and [`docs/steps/`](./steps/); this is the running tasks.
A standalone write-up of the current state is `../report/` (LaTeX → `main.pdf`).

## Next up — prioritized (15.07.2026, from the progress-report feedback)

> **Framing that governs this list:** *more performance ≠ a bigger model.* Model-side tuning is
> demonstrably closed (see "Model-side tuning is closed" below); the levers are label-/data-side.
> The audience's "bigger MLP / more capacity" suggestion was already tested and is flat — only MIL
> (S2) is a genuinely untested capacity lever.
>
> **Update 25.07.2026:** the report's first next-step — *drug selection from the literature instead of my
> filter* — is **done as a definition** (8-drug panel, see "Next focus" below); the training run on it is
> pending and should be paired with S1. The panel is literature-anchored but **not yet label-blind**, so
> train-only selection remains blocking for any headline number.

1. **S1 — DrEval-aligned target (the top performance lever).** Train on the double-normalized residual
   `resid[i,j] = auc[i,j] − (μ_drug[j] + μ_line[i] − μ_global)`, means computed **train-only per fold**,
   so the objective matches DrEval's normalized metric. Motivated by `dreval/dreval_normalized.csv`
   (~20% of the signal is pure cell-line effect; `kx2-391` is entirely that artifact). New `--score`
   option in `ctrp_to_h5ad.py` (pattern: `_zscore_per_drug`, `DEFAULT_CTRP_SCORE`); compare in `12`.
   Success = normalized DrEval ρ rises above 0.357 (scGPT) without inflating the raw correlation.
2. **S2 — MIL / attention pooling over a line's cells** — the only untested capacity lever (bag of cells
   → line label). Must beat the ridge baseline (0.342 PCA). *(Detailed item under "Model-side tuning".)*
3. **S3 — More independent cell lines** — SCP542×CTRPv2 caps at 180; CTRPv2 has ~1,100. Attacks the
   real ceiling. *(Overlaps the scDEAL/label-side lever under "Levers / later".)*
4. **S4 — Diagnostic explainability (now, low-risk)** — where errors concentrate (drugs/tissues/lines
   with ρ<0) and how much residual error is the cell-line effect. Uses existing per-drug-ρ CSVs.
5. **Later — Discovery explainability (gated)** — gene-level XAI only once per-drug ρ is substantially
   higher and stable, else it interprets noise. *(Project-plan stretch goal; see Step 07.)*

**Communication fixes:** rebuild the single overview image (data → rep → shared trunk → K heads →
out-of-fold eval); make the rescue figure show the **ceiling** (ridge = MLP, no-reg memorizes), because
the "model-side is closed" message did not land in the talk.

**Working agreement:** for any analysis beyond the explicitly agreed fix — *especially how a plot is
computed/displayed* — confirm first, don't decide silently. (Unasked decisions produced process bugs and
undefendable slides last round.)

## Agreed plan — order of work (27.07.2026)

**Governing rule: never change the target and the architecture in the same run.** The June result took
weeks to unpick precisely because two changes landed together; if MIL and a new target arrive at once and
the number moves, the cause is unattributable.

**Step 1 — target + loss weighting, on the existing per-cell MLP.** Everything else unchanged
(architecture, splits, optimizer, batching), so the change is attributable.

**Decision (27.07.2026): `auc_z` is retired as the target.** Its centering is inert (the per-drug head
bias absorbs it) and its scaling is the defect. Target becomes raw `auc`; the scaling moves into the loss
as an explicit per-drug weight, where it can be estimated per fold and combined with other factors.

- **1.0 — prerequisite:** estimate `σ_noise` pooled from the 7,708 replicated curve fits in
  `v20.data.curves_post_qc.txt` (2.0 % of 387,130). Per-drug is not feasible; the pooled estimate assumes
  homoscedastic assay noise — state that, and sanity-check whether replicate discrepancy tracks a drug's
  mean AUC.
- **1.1 — refactor, no scientific change.** Target = raw `auc`; per-drug weight `w_j = 1/σ_j²` with σ
  estimated **per fold on training lines only** (this also retires the standing z-score leak, since the
  statistics are no longer baked into the targets h5ad). Three traps, all verified in the code and all
  silent if missed:
  1. `optim.Adam(model.parameters(), ..., weight_decay=1e-3)` (`training_utils.py:179`) decays **all**
     parameters including head biases. On `auc_z` the optimal bias was 0 so decay was free; on raw `auc`
     the bias must sit near 0.7 and is actively pulled to 0. Put biases + LayerNorm params in a
     `weight_decay=0` group.
  2. `nn.Linear` initializes head biases at ±0.125 against a target near 0.7 — initialize them to the
     train-fold per-drug means instead.
  3. Per-drug Spearman is invariant to within-drug affine transforms, so for a *fixed* model this setup
     and the old `auc_z` score identically.
  > **Pass/fail:** must reproduce the current `auc_z` numbers. Any gap is either an implementation bug or
  > the size of the leak — both worth knowing before measuring anything new.
- **1.2 — the experiment:** change the weight from `1/σ_j²` to `r_j/σ_j²`, `r_j = (σ_j² − σ_noise²)/σ_j²`.
  Exactly one change, against a verified baseline. Rationale: the metric averages per-drug Spearman with
  equal weight per drug, so the loss should too — discounted by the share of each drug's variance that is
  real rather than assay noise.

> **Superseded 27.07.2026 — 1.0/1.1/1.2 above were never needed.** Per-drug variance weighting is a
> *K=545* problem; on the 8-drug panel the variance ratio is 2.5× and the raw target works unweighted
> (confirmed empirically: raw `auc` scores −0.069 at 545 heads and **+0.377** at 8). `σ_noise` was
> therefore never estimated. Kept above only so the reasoning is not re-derived from scratch.

### Step 1 — DONE (27.07.2026), `notebooks/14_panel_training.ipynb`

Raw `auc` winsorized at 1.1, 8-drug panel, per-sample inverse-density weights fitted per fold on training
lines only, output layer excluded from weight decay (`TrainConfig.exclude_output_from_decay`, default off
so old runs are unchanged), head biases initialized to train-fold per-drug means. One seed.

- [x] **Confirmed:** the June collapse was a K=545 effect, not a target property → retiring `auc_z` is free.
- [x] **Confirmed (replication):** PCA MLP ties its ridge (0.316 vs 0.306); scGPT MLP clears its ridge by
      **+0.077**, against +0.082 on the 14.07 panel — now on a drug set chosen without our labels.
- [x] **Refuted:** inverse-density loss weighting (−0.006 / −0.008). Mechanism fired (pred_std 0.062 →
      0.08) but ranking did not follow — after winsorizing, |skew| ≤ 0.47, so there was no imbalance left
      to correct. **Do not carry into Step 2.**
- [ ] **Reproducibility:** the PCA-unweighted arm is not bit-reproducible on `mps` — three identical
      runs gave 0.313 / 0.315 / 0.317 / 0.320 (four draws), every other arm reproduced exactly. Cause: PCA peaks at epoch 1
      (`[1,1,3,1,1]` vs scGPT `[10,11,2,21,4]`), so its checkpoint is chosen among near-tied states. The
      weighting deltas lie inside that band, so **do not report their sign**.
- [ ] **Seeds.** Everything above is one seed against ±0.04 documented seed variation. Repeat over ≥3
      seeds before scGPT − PCA (+0.061) or scGPT − ridge (+0.077) is quoted as a margin. **Blocking for
      any headline number.**
- [ ] **Report raw + normalized** (DrEval) on this panel, once seeds are in.

**Step 2 — MIL / attention pooling, against the target fixed in Step 1.** Controls that must both be
beaten: the per-cell MLP and RidgeCV on cell-line mean embeddings. If MIL beats neither, the single-cell
resolution has again failed to justify itself and that is the reportable result.

- MIL makes the per-line weighting problem disappear structurally (one bag = one example), so the
  82× cells×labels imbalance needs no separate fix under it.
- Open design decisions, to settle before building: fixed bag size with per-epoch subsampling (acts as
  augmentation) vs. padding + mask, given 56–1,990 cells per line; and the optimizer regime, since an
  epoch becomes ~120 bags instead of ~34,000 cells so the current epoch/LR/early-stop settings do not
  carry over.
- Attention weights are the readout for *which subpopulation drives response* — the link to the relapse
  motivation and to the annotated heterogeneity programs.

**Not in scope for either step:** the base-quantity question (AUC vs EC50/Emax, T3) and learned/adaptive
task weights. Both are deferred deliberately; adaptive weights estimate *residual* variance, which mixes
label noise with model error and risks a self-reinforcing loop, especially at 545 tasks over ~120 bags.

## Target & drug-selection defects found 27.07.2026 (do these before any new headline number)

Both found by asking why `nutlin-3` was rejected by the filter. Write-ups:
[Step 03](./steps/03-model-and-training-design.md#known-problems-with-auc_z--the-scaling-is-not-yet-right-27072026),
[Step 05](./steps/05-multitask-results.md#the-learnability-gate-measured-the-wrong-quantity-27072026).

- [ ] **T1 — Replace the kill/spare gate with `auc_std` + coverage.** The gate filters on absolute
      potency; `auc_z` subtracts the per-drug mean and Spearman reads only the ordering, so we selected
      on a quantity the model never sees. `nutlin-3` σ = 0.147 vs `dasatinib` σ = 0.155 — same spread,
      rejected only because it is cytostatic. **116/545** drugs have zero kills but σ ≥ 0.10 and
      coverage ≥ 90 %. Everything downstream (10-drug panel, 8-drug literature panel, all K=10 numbers)
      rests on the old gate and must be re-derived.
- [ ] **T2 — Fix the `auc_z` denominator.** Dividing by `auc_std` forces noise-floor drugs to variance 1
      and hands them full weight in the shared loss — the mirror of the June σ² bug. Use
      `sqrt(auc_std² + σ_noise²)`, or weight each drug by its reliable variance fraction.
      `σ_noise` is estimable **pooled** from the 7,708 replicated (line × compound) fits (2.0 % of
      387,130) in `v20.data.curves_post_qc.txt`; per-drug is not feasible.
- [ ] **T3 — Reconsider AUC as the target** (raised by Selin's supervisor, DrEval co-author; DrEval lists
      inconsistent viability data as an obstacle and recommends **CurveCurator**). AUC conflates potency
      with efficacy, and CTRP's own fit already separates them in the file we parse:
      `apparent_ec50_umol`, `pred_pv_high_conc` (≈ Emax), `p3_total_decline`. Top test concentration
      spans **0.13–600 µM** across the 545 — harmless across drugs (z-scoring is within-drug) but it
      compresses spread *within* a drug, so it is a **cause of T2**, not a separate issue.
      Interacts with S1 — decide the target once, for both.

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
- [x] **Target distribution** (data only) → `04` (`outputs/data/target_distribution.png`): viability clusters
      near 1.0 (median 0.91; 75% ≥ 0.8); per-drug std median 0.088, only 3% truly flat; a loose
      cov ≥ 100 & std ≥ 0.05 filter keeps **439/545** → coverage+std alone removes few.
- [x] **Coverage & learnability** analysis → `04` (`outputs/*_drug_learnability.csv`, `data/drug_coverage.png`).
- [x] **Cancer-type UMAPs** → `06` §8: 2-panel PCA-vs-scGPT (`outputs/embeddings/umap_cancertype_pca_vs_scgpt.png`,
      dpi 300) + full gene-set sweep grid (`outputs/embeddings/umap_sweep_cancertype.png`, dpi 200). Tissue islands
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
      *(14.07.2026: filter widened to **10 drugs** (coverage ≥ 90 %); re-run gives **0.360 / 0.396** —
      see "Next up" above.)*

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

- [x] **Literature-anchored panel defined** (25.07.2026) — 8 drugs chosen from *published* cell-line
      sensitivity determinants instead of my gates, so selection never sees our labels:
      `methotrexate`, `dasatinib`, `paclitaxel`, `vincristine`, `afatinib`, `topotecan`, `tanespimycin`,
      `selumetinib`. Rule = CTRPv2 ∩ `compound_status ∈ {FDA, clinical}` (173/545) ∩ published
      determinant. **All eight already pass the `08` gate**, so it is a re-ranking inside the
      gate-passing set, not a loosening; only 2 overlap the old ten. Table + citations in
      [Step 05](./steps/05-multitask-results.md#literature-anchored-drug-panel--selecting-without-looking-at-our-labels-25072026),
      decision log in [project_notes](./project_notes.md). **Not trained yet — that is the next run.**
      Six of the eight are expression-determined and two (`selumetinib`, `afatinib`) are
      mutation-determined, so the panel doubles as a hypothesis test on the representation.
      ⚠ **Literature-anchored, spread-verified — not label-blind:** the candidate list was ranked by
      `min(kill, spare)` on our own AUCs before the literature criterion was applied, so drugs with
      published determinants but little spread here dropped out (`sirolimus`, `neratinib`,
      `clofarabine`, `cytarabine hydrochloride`, `gdc-0941`). Fixed: arbitrariness of *which* drugs +
      threshold instability. Not fixed: the optimistic component in per-drug ρ.
- [ ] **Train-only selection** — run the learnability gates *inside each CV fold* (train lines only) and
      re-measure. If the effect survives, it is real; **this is still the blocking item**, including for
      the literature panel (see the caveat above).
- [ ] **Externalize the spread requirement** — re-derive the panel with the kill/spare condition measured
      on **GDSC2** (`data/GDSC2_fitted_dose_response_27Oct23.xlsx`) or PRISM instead of on the CTRP labels
      we train on. Keeps "the drug must actually kill something" while removing the leak; cheaper than
      the fold-internal selection above and would make the panel genuinely label-blind.
- [x] **Multiple seeds** (13.07.2026) — done for K=545 `auc_z`: gap **+0.075 ± 0.038**, sign-consistent
      over 3 seeds. Still thin; see the "more seeds + wider drug set" item above.
- [ ] **Loosen to ~20–50 drugs** — 5 is a diagnostic, not a model. Where does the signal die as the gates
      relax? (`learnability/ctrp_drug_learnability_auc.csv` is already ranked for this.)
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

`notebooks/10_diagnosis.ipynb`: regularization (none→heavy), capacity (74,629→2,565 params), batch size
(32/128/512) and reweighting (line-balanced, focus-extremes) **all leave Spearman flat** (PCA 0.41–0.44,
scGPT 0.44–0.49); the current defaults are already at/near the best on every axis. With regularization
off, PCA memorizes the train lines (train MSE ≈ 0.01) and still only reaches 0.42 out-of-fold — it is out
of *signal*, not out of *capacity*. **Don't spend more time on architecture or hyperparameters.**

**New baseline to beat: `RidgeCV` on the 150 cell-line mean embeddings** (ρ = **0.428**) — it *ties* the
PCA MLP and comes within 0.06 of the scGPT MLP, with no single cells and no network. The per-drug-mean
null is too weak a bar; report ridge alongside it from now on. scGPT + hidden layer (0.487) is the only
configuration that clears it — and scGPT's *linear* head drops to 0.438, so it genuinely needs the
nonlinearity while PCA does not.
*(14.07.2026, 10-drug re-run — same conclusion, smaller numbers: ridge **0.343** ties the PCA MLP
**0.356**; scGPT MLP **0.402**, scGPT linear **0.292**; `ablations/ablation_capacity.csv`.)*

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

- [x] **Benchmarked with the REAL package** (`pip install drevalpy` 1.5.1) → `notebooks/12_dreval_benchmark.ipynb`,
      `outputs/dreval/dreval_lco_results.csv`. Their `DrugResponseDataset` + LCO splits + baselines + `evaluate()`.
      **OncoMLP (scGPT) clears `NaiveMeanEffects`: normalized ρ = 0.511 ± 0.085, normalized R² = 0.224** —
      vs the 11% (DIPK) / 19% (RF) the paper reports for its best LCO models. **scGPT > PCA confirmed
      externally** (+0.07). And our per-cell MLP **beats their `SingleDrugRandomForest` on the same
      embeddings** (0.511 vs 0.438) → qualifies the ridge≈MLP result.
      *(⚠️ 14.07.2026 — the 0.511/0.224 above was the 5-drug best-case **with** a val-split leak. After
      fixing the leak (`ee07b00`, test fold was the val loader) and re-running on the 10-drug set, the
      mean over 5 folds is **scGPT ρ 0.357 / R² 0.114**, PCA 0.340 / 0.086, their `SingleDrugRF` (scgpt)
      0.339 — still above naive, but only faintly above PCA and their RF. Use these.)*
- [x] **Own-implementation check** (`outputs/dreval/dreval_normalized.csv`): additionally removing the *cell-line*
      effect (which LCO's naive predictor cannot know) still leaves scGPT at ρ = 0.396. ⚠️ **`kx2-391`
      collapses to 0.006** — its signal was entirely the cell-line effect.
- [ ] **Make `NaiveMeanEffects` the default baseline** in `train_multitask.py` (currently: per-drug mean,
      too weak).
- [ ] **Report raw + normalized** correlations everywhere from now on.
- [ ] **Run DrEval on all 545 drugs**, not just the best-case 5 — and with their LTO / LDO settings.
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
