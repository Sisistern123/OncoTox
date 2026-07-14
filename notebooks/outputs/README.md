# `notebooks/outputs/` — figures & tables

Grouped by pipeline stage, so a file's provenance is never in doubt. Superseded artifacts live under
`legacy/` and are **not** cited by any current claim. Model weights and per-run configs are not here —
they live in `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

## Current

| Directory | Written by | Contents |
|---|---|---|
| **`data/`** | `04_drug_coverage`, `10_diagnosis` §1 | Target distribution and per-drug coverage on the raw CTRPv2 labels; `target_biology.png` = what the targets measure on a real dose–response curve. |
| **`embeddings/`** | `06_verify_variants` | Latent-space validation: PCA-vs-scGPT UMAPs, gene-set variants. |
| **`learnability/`** | `08_learnability_filter`, `09_learnable5_training` | The drug filter (545 → 5) and the PCA-vs-scGPT result on that subset. |
| **`target/`** | `11_auc_vs_aucz`, `10_diagnosis` §2 | Which target to train on; the per-drug **loss-weighting bug**; seed stability. |
| **`ablations/`** | `10_diagnosis` §3–§5 | The **causal rescue test** on the broken K=545 setting, the model-knob ablations on the corrected one, and the ridge control. |
| **`dreval/`** | `12_dreval_benchmark` | External benchmark against **DrEval** (`drevalpy` 1.5.1): their LCO splits, baselines, metrics. |

### The five figures that carry the current story

| Figure | What it shows |
|---|---|
| `data/target_biology.png` | On one real CTRPv2 curve: `mean_pv` averages the measured points, `auc` integrates the fitted sigmoid. They land within ~0.03 of each other — *why* the curve fit alone buys no accuracy. |
| `target/loss_weighting_bug.png` | An unweighted MSE weights each drug by σ². The widest 10% of drugs carry **30%** of the loss — and the very widest (`ifosfamide`, `ciclopirox`) kill **zero** cell lines. |
| `target/target_comparison.png` | `mean_pv` / `auc` / `auc_z` at K=5 and K=545, with bootstrap CIs. Both unstandardized targets collapse at K=545; `auc_z` holds. |
| `ablations/rescue_k545.png` | Every June hypothesis applied to the **broken** setting. Only *task* reweighting fixes it (+0.43); removing regularization partially rescues (+0.23) — the symptom, not the cause. |
| `dreval/dreval_lco.png` | OncoMLP vs the DrEval baselines under their LCO protocol — it clears `NaiveMeanEffects`, which half the published field does not. |

### Two DrEval tables that look alike but answer different questions

- `dreval/dreval_lco_results.csv` — the **real package** (`drevalpy`), their LCO splits and baselines.
  In LCO the naive predictor cannot know a held-out line's effect, so it reduces to *global mean + drug
  effect*; its normalized metric therefore removes the **drug** mean.
- `dreval/dreval_normalized.csv` — our **own, stricter** variant, which *additionally* removes the
  **cell-line** effect (using that line's own labels). Not what DrEval does. It answers a separate
  question: *"how much of our signal is merely 'this line is sensitive to everything'?"* — answer: ~20%,
  and `kx2-391` is **entirely** that.

## `legacy/` — superseded, kept as the record of what was overturned

| Path | Why it is legacy |
|---|---|
| `legacy/training_545_mean_pv/` | The 8-run matrix, 5-fold CV, HVG sweep and per-drug correlations from `07_training`. Produced at **K=545 on `mean_pv`**, i.e. with the unstandardized loss that `10_diagnosis` shows was destroying the signal. Its "ρ ≈ 0" conclusion does not survive; the numbers are kept only as evidence of the failure mode. |
| `legacy/ctrp_drug_learnability_mean_pv.csv` | The old learnability table (`resp_std × coverage`) on `mean_pv`. Its score is **degenerate** on the z-scored target and its gates kept 439/545, so it never bit. Replaced by `learnability/ctrp_drug_learnability_auc.csv`. |
| `legacy/gdsc_drug_learnability.csv` | GDSC2 learnability list, produced once for Hashimoto-san. Never consumed by the modelling pipeline. |
