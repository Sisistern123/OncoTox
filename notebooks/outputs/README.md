# `notebooks/outputs/` — figures & tables

Grouped by pipeline stage. Each subdirectory says which notebook writes it, so a file's provenance is
never in doubt. Model artifacts themselves (weights, per-run configs) are **not** here — they live in
`runs/` (git-ignored) and are indexed in `runs/runs_index.csv`.

| Directory | Written by | Contents |
|---|---|---|
| **`01_data/`** | `04_drug_coverage.ipynb` | Target distribution & per-drug coverage on the raw CTRPv2 labels. `target_biology.png` shows what `mean_pv` and `AUC` actually measure on a real dose–response curve. |
| **`02_embeddings/`** | `06_verify_variants.ipynb` | Latent-space validation: PCA-vs-scGPT UMAPs, gene-set variants. |
| **`03_training_545/`** | `07_training.ipynb` | The original 8-run matrix, 5-fold CV, HVG sweep, per-drug correlations. ⚠️ **All on the legacy `mean_pv` target at K=545 — superseded** (see `05_target/`). |
| **`04_learnability/`** | `08_learnability_filter.ipynb`, `09_learnable5_training.ipynb` | The harsh drug filter (545 → 5) and the PCA-vs-scGPT result on the learnable subset. |
| **`05_target/`** | `11_auc_vs_aucz.ipynb` | Which target to train on (`mean_pv` / `auc` / `auc_z`), the per-drug **loss-weighting bug**, seed stability. |
| **`06_ablations/`** | `10_ablations.ipynb` | What is **not** the bottleneck: regularization, capacity, batch size, sample reweighting — plus the ridge-on-line-means control and the K=545 rescue attempts. |
| **`07_dreval/`** | `12_dreval_benchmark.ipynb` | External benchmark against **DrEval** (`drevalpy` 1.5.1): their LCO splits, their baselines, their metrics. |

## The four figures that carry the current story

| Figure | Shows |
|---|---|
| `01_data/target_biology.png` | On one real CTRPv2 dose–response curve: `mean_pv` averages the measured points, AUC integrates the fitted sigmoid. They land within ~0.03 of each other — which is *why* the curve fit alone buys no accuracy. |
| `05_target/loss_weighting_bug.png` | An unweighted MSE weights each drug by σ². The widest 10% of drugs carry 30% of the loss — and the very widest kill no cell lines at all. |
| `06_ablations/what_did_not_work.png` | Regularization, model size, batch size and sample reweighting are flat; only re-weighting the drug heads recovers the signal. |
| `07_dreval/dreval_lco.png` | OncoMLP vs the DrEval baselines under their LCO protocol. |
