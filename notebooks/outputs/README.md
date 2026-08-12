# `notebooks/outputs/` — figures & tables

Grouped by pipeline stage, so a file's provenance is never in doubt. Model weights and per-run configs
are not here — they live in `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

**The dividing line, applied 12.08.2026 (Selin): can a standard pipeline run recreate this?** If the
notebook that wrote it is archived, or its criterion was retracted, the answer is permanently no and
the artifact lives under `legacy/`. Everything outside `legacy/` is reproducible by running the
numbered notebooks — or, for `diagnostics/` and `dreval/`, will be once a one-line fix to the
producing notebook lands. Those two are **blocked, not dead**, which is why they did not move.

## Current

| Directory | Written by | Recreatable? | Contents |
|---|---|---|---|
| **`data/`** | `analysis/harmonization/drug_coverage` | **yes** | Target distribution and per-drug coverage on the raw CTRPv2 labels. |
| **`embeddings/`** | `analysis/qc/verify_variants`, `gene_symbol_rescue` | **yes** | Latent-space validation: PCA-vs-scGPT UMAPs, gene-set variants, and what the gene set actually delivered to scGPT (`gene_symbol_rescue.csv`, below). |
| **`panel/`** | `2_drug_selection` | **yes** | `panel.csv` and `literature_panel_candidates.csv` — the [rebuilt 11-drug panel](../../docs/steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026) and the 57 candidates it was drawn from, each with its reference and what that reference establishes. The stage reads only the response CSV, so it runs under the freeze. |
| **`dreval/`** | `analysis/evaluation/dreval_benchmark` | ⛔ **blocked** | External benchmark against **DrEval** (`drevalpy` 1.5.1): their LCO splits, baselines, metrics. The notebook imports a module deleted on 12.08.2026 and hardcodes the removed `'auc'`, so it raises on its first cell. Review item 11. |
| **`diagnostics/`** | `analysis/evaluation/diagnostics` | ⛔ **blocked** | The drug-selection gate defect (`gate_per_drug.csv`, `gate_potency_vs_spread.png`), the proliferation test (`line_effect_vs_programs.csv`, `line_effect_vs_proliferation.png`), the input-scale asymmetry (`input_scale.csv`), and result dispersion (`result_dispersion.csv`). Hardcodes the removed `'auc'` and raises. Its §5 dispersion figures were computed on the void panel. |

### The figures that carry the argument

⚠️ **Retitled 12.08.2026 — three of these four are now under `legacy/`.** This was headed *"the five
figures that carry the current story"*, which stopped being true when the artifacts a standard run
cannot recreate moved. They still carry the argument that got the project here; they are no longer
current results, and only the last is at a live path.

| Figure | What it shows |
|---|---|
| `legacy/target/loss_weighting_bug.png` | An unweighted MSE weights each drug by σ². The widest 10% of drugs carry **30%** of the loss — and the very widest (`ifosfamide`, `ciclopirox`) kill **zero** cell lines. |
| `legacy/target/target_comparison.png` | `mean_pv` / `auc` / `auc_z` at K=5 and K=545, with bootstrap CIs. Both unstandardized targets collapse at K=545; `auc_z` holds. ⛔ All three targets have since been removed. |
| `legacy/ablations/rescue_k545.png` | Every June hypothesis applied to the **broken** setting. Only *task* reweighting fixes it (+0.43); removing regularization partially rescues (+0.23) — the symptom, not the cause. |
| `dreval/dreval_lco.png` | OncoMLP vs the DrEval baselines under their LCO protocol — it clears `NaiveMeanEffects`, which half the published field does not. Still at its original path: `dreval_benchmark` re-runs at R5. |

### Two DrEval tables that look alike but answer different questions

- `dreval/dreval_lco_results.csv` — the **real package** (`drevalpy`), their LCO splits and baselines.
  In LCO the naive predictor cannot know a held-out line's effect, so it reduces to *global mean + drug
  effect*; its normalized metric therefore removes the **drug** mean.
- `dreval/dreval_normalized.csv` — ⛔ **void, and no longer reproducible as written.** It was our
  **own, stricter** variant, which *additionally* removed the **cell-line** effect using that line's own
  labels — not what DrEval does — answering *"how much of our signal is merely 'this line is sensitive
  to everything'?"* That metric was **deleted** on 12.08.2026 as a local invention, and the producing
  script cut back to DrEval's recipe
  ([why](../../scripts/archive/README.md)), so the current `dreval_normalize.py` writes to this path but
  computes the **paper's** metric, not the stricter one. Whether the stricter diagnostic returns is for
  audit 11.

  ⚠️ Worth knowing before that decision: under our leave-cell-line-out splits DrEval's normalization
  removes only the **drug** effect, because a held-out line's effect is unseen and therefore zero. A
  synthetic predictor emitting nothing but mean + line effect + drug effect scores normalized Spearman
  **0.98**. A high normalized score is not evidence of drug-specific signal here.

## `legacy/` — what a standard run cannot recreate

Two grounds, and both are permanent: the notebook that wrote it is **archived**, or its **criterion
was retracted**. Nothing here is regenerated by re-running the pipeline, and nothing here is a live
result. Five directories moved in on 12.08.2026 (Selin) — before that they sat alongside the current
outputs with only a prose warning, which is not a structure a reader can rely on.

| Path | Why it is legacy |
|---|---|
| `legacy/learnability/` | **Moved 12.08.2026.** The drug filter (545 → 10, `learnability = min(#killed, #spared)`) and the PCA-vs-scGPT result on that subset. The criterion was [retracted](../../docs/steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability) — it measured potency, not rankability — and both producing notebooks are archived. |
| `legacy/target/` | **Moved 12.08.2026.** Which target to train on; the per-drug **loss-weighting bug**; seed stability. From `target_comparison` and `ablations_and_rescue` §2, both archived: the comparison is between `mean_pv`, `auc` and `auc_z`, all three removed on 11.08.2026. |
| `legacy/ablations/` | **Moved 12.08.2026.** The **causal rescue test** on the broken K=545 setting, the model-knob ablations, and the ridge control. From `ablations_and_rescue`, archived. ⛔ Its ablation and ridge tables are void on a second ground: the `oof()` it used early-stopped on the fold it scored. |
| `legacy/panel_void_8drug/` | **Moved 12.08.2026.** Everything `4_training` wrote on the [voided 8-drug panel](../../docs/steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it) — response distributions, weight curves, out-of-fold predictions, the ridge baseline and the per-drug correlations. The rebuilt panel shares only three compounds with it, so these are not earlier versions of what R4 will produce. |
| `legacy/replicate_variation.{csv,png}` | **Moved 12.08.2026.** How far CTRPv2's repeated measurements diverge — 2,637 pairs screened twice, median disagreement 0.49× the drug's spread across cell lines. Its notebook is archived: it reads CTRPv2's own `v20.*` tables, which stopped being the target source on 11.08.2026. The measurement stands as the evidence for that switch ([Step 01](../../docs/steps/01-datasets-and-harmonization.md#genuine-repeats-are-averaged-and-they-disagree-more-than-the-targets-own-spread-10082026)). |
| `legacy/training_545_mean_pv/` | The 8-run matrix, 5-fold CV, HVG sweep and per-drug correlations from what is now `4_training` §B. Produced at **K=545 on `mean_pv`**, i.e. with the unstandardized loss that `ablations_and_rescue` shows was destroying the signal. Its "ρ ≈ 0" conclusion does not survive; the numbers are kept only as evidence of the failure mode. |
| `legacy/ctrp_drug_learnability_mean_pv.csv` | The old learnability table (`resp_std × coverage`) on `mean_pv`. Its score is **degenerate** on the z-scored target and its gates kept 439/545, so it never bit. Superseded by `legacy/learnability/ctrp_drug_learnability_auc.csv`, itself retracted. |
| `legacy/gdsc_drug_learnability.csv` | GDSC2 learnability list, produced once and shared on request outside this project. Never consumed by the modelling pipeline. |
