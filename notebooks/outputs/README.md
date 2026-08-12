# `notebooks/outputs/` — figures & tables

Grouped by pipeline stage, so a file's provenance is never in doubt. Superseded artifacts live under
`legacy/` and are **not** cited by any current claim. Model weights and per-run configs are not here —
they live in `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

## Current

| Directory | Written by | Contents |
|---|---|---|
| **`data/`** | `drug_coverage`, `replicate_variation` | Target distribution and per-drug coverage on the raw CTRPv2 labels, plus how far CTRPv2's repeated measurements diverge (`replicate_variation.png/.csv` — 2,637 pairs screened twice, median disagreement 0.49× the drug's spread across cell lines). |
| **`embeddings/`** | `verify_variants`, `gene_symbol_rescue` | Latent-space validation: PCA-vs-scGPT UMAPs, gene-set variants, and what the gene set actually delivered to scGPT (`gene_symbol_rescue.csv`, below). |
| **`learnability/`** | *(archived notebooks)* | ⛔ The drug filter (545 → 10, `learnability = min(#killed, #spared)`) and the PCA-vs-scGPT result on that subset. **The criterion was [retracted](../../docs/steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability)** and both producing notebooks archived 12.08.2026; nothing regenerates these. |
| **`target/`** | `target_comparison`, `ablations_and_rescue` §2 | Which target to train on; the per-drug **loss-weighting bug**; seed stability. |
| **`ablations/`** | `ablations_and_rescue` §3–§5 | The **causal rescue test** on the broken K=545 setting, the model-knob ablations on the corrected one, and the ridge control. |
| **`dreval/`** | `dreval_benchmark` | External benchmark against **DrEval** (`drevalpy` 1.5.1): their LCO splits, baselines, metrics. |
| **`diagnostics/`** | `diagnostics` | The drug-selection gate defect (`gate_per_drug.csv`, `gate_potency_vs_spread.png`), the proliferation test (`line_effect_vs_programs.csv`, `line_effect_vs_proliferation.png`), the input-scale asymmetry (`input_scale.csv`), and result dispersion (`result_dispersion.csv`). |
| **`panel/`** | `literature_panel`, `3_panel_training` | **`panel.csv` and `literature_panel_candidates.csv` are current** — the [rebuilt 11-drug panel](../../docs/steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026) and the 57 candidates it was drawn from, each with its reference and what that reference establishes. ⛔ **Everything else in the directory is computed on the [voided 8-drug panel](../../docs/steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it)** — the response distributions, weight curves and training result. Do not quote those; they re-run at R4. |

### The five figures that carry the current story

| Figure | What it shows |
|---|---|
| `target/loss_weighting_bug.png` | An unweighted MSE weights each drug by σ². The widest 10% of drugs carry **30%** of the loss — and the very widest (`ifosfamide`, `ciclopirox`) kill **zero** cell lines. |
| `target/target_comparison.png` | `mean_pv` / `auc` / `auc_z` at K=5 and K=545, with bootstrap CIs. Both unstandardized targets collapse at K=545; `auc_z` holds. |
| `ablations/rescue_k545.png` | Every June hypothesis applied to the **broken** setting. Only *task* reweighting fixes it (+0.43); removing regularization partially rescues (+0.23) — the symptom, not the cause. |
| `dreval/dreval_lco.png` | OncoMLP vs the DrEval baselines under their LCO protocol — it clears `NaiveMeanEffects`, which half the published field does not. |

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

## `legacy/` — superseded, kept as the record of what was overturned

| Path | Why it is legacy |
|---|---|
| `legacy/training_545_mean_pv/` | The 8-run matrix, 5-fold CV, HVG sweep and per-drug correlations from `2_training`. Produced at **K=545 on `mean_pv`**, i.e. with the unstandardized loss that `ablations_and_rescue` shows was destroying the signal. Its "ρ ≈ 0" conclusion does not survive; the numbers are kept only as evidence of the failure mode. |
| `legacy/ctrp_drug_learnability_mean_pv.csv` | The old learnability table (`resp_std × coverage`) on `mean_pv`. Its score is **degenerate** on the z-scored target and its gates kept 439/545, so it never bit. Replaced by `learnability/ctrp_drug_learnability_auc.csv`. |
| `legacy/gdsc_drug_learnability.csv` | GDSC2 learnability list, produced once and shared on request outside this project. Never consumed by the modelling pipeline. |
