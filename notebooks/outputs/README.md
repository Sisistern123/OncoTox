# `notebooks/outputs/` — figures & tables

Grouped by pipeline stage, so a file's provenance is never in doubt. Model artifacts (weights, per-run
configs) are **not** here — they live in `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

| Directory | Written by | Contents |
|---|---|---|
| **`01_data/`** | `04_drug_coverage`, `10_diagnosis` (§1) | Target distribution and per-drug coverage on the raw CTRPv2 labels; `target_biology.png` = what the targets measure on a real dose–response curve. |
| **`02_embeddings/`** | `06_verify_variants` | Latent-space validation: PCA-vs-scGPT UMAPs, gene-set variants. |
| **`03_training_545/`** | `07_training` | The original 8-run matrix, 5-fold CV, HVG sweep, per-drug correlations. ⚠️ **Legacy `mean_pv` at K=545 — superseded.** |
| **`04_learnability/`** | `08_learnability_filter`, `09_learnable5_training` | The drug filter (545 → 5) and the PCA-vs-scGPT result on that subset. |
| **`05_target/`** | `11_auc_vs_aucz`, `10_diagnosis` (§2) | Which target to train on; the per-drug **loss-weighting bug**; seed stability. |
| **`06_ablations/`** | `10_diagnosis` (§3–§5) | The **causal rescue test** on the broken K=545 setting, the model-knob ablations on the corrected one, and the ridge-on-line-means control. |
| **`07_dreval/`** | `12_dreval_benchmark` | External benchmark against **DrEval** (`drevalpy` 1.5.1): their LCO splits, baselines, metrics. |

## The five figures that carry the current story

| Figure | What it shows |
|---|---|
| `01_data/target_biology.png` | On one real CTRPv2 curve: `mean_pv` averages the measured points, `auc` integrates the fitted sigmoid. They land within ~0.03 of each other — *why* the curve fit alone buys no accuracy. |
| `05_target/loss_weighting_bug.png` | An unweighted MSE weights each drug by σ². The widest 10% of drugs carry **30%** of the loss — and the very widest (`ifosfamide`, `ciclopirox`) kill **zero** cell lines. |
| `05_target/target_comparison.png` | `mean_pv` / `auc` / `auc_z` at K=5 and K=545, with bootstrap CIs. Both unstandardized targets collapse at K=545; `auc_z` holds. |
| `06_ablations/rescue_k545.png` | Every June hypothesis applied to the **broken** setting. Only *task* reweighting fixes it (+0.43); removing regularization partially rescues (+0.23) — the symptom, not the cause. |
| `07_dreval/dreval_lco.png` | OncoMLP vs the DrEval baselines under their LCO protocol — it clears `NaiveMeanEffects`, which half the published field does not. |

## Two tables that look similar but answer different questions

- `07_dreval/dreval_lco_results.csv` — the **real package** (`drevalpy`), their LCO splits and baselines.
  In LCO the naive predictor cannot know a held-out line's effect, so it reduces to *global mean + drug
  effect*; the normalized metric therefore removes the **drug** mean.
- `07_dreval/dreval_normalized.csv` — our **own** stricter variant, which *additionally* removes the
  **cell-line** effect (using that line's own labels). Not what DrEval does — it answers the separate
  question *"how much of our signal is just 'this line is sensitive to everything'?"* (Answer: ~20%; and
  `kx2-391` is entirely that.)
