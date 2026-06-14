
# OncoTox Project Notes
## 10–14.06.2026

### Documentation restructure
Split the source-of-truth into a thin index `project_progress.md` + eight step files under
`docs/steps/` (01 datasets/harmonization … 08 foundation model). Added `docs/TODO.md` and a pipeline
overview figure.

### PCA fix — compute on the HVG counts, not the OOV-dropped `.X`
Found `add_pca` was computing `X_pca` on the targets `.X`, which the `scgpt` step had already reduced
to scGPT's in-vocab genes (5,000→4,576; 22,722→20,570). Changed it to compute PCA from the **convert
counts** (full HVG-5000 / 22,722 genes) and to leave `.X` as CPM. **Decision:** the scGPT OOV drop
is kept as part of the comparison (it's a real property of using scGPT), so PCA uses the full
filtered set and scGPT its in-vocab subset.

### Matched trunk + 8-run matrix rerun
Set `DEFAULT_HIDDEN_DIMS` to `(128,64)` for **both** reps (was (64,32) for PCA) so head capacity is
equal — only the input representation differs. Re-ran the full 8-run matrix
`{hvg5000, all_genes} × {X_pca, X_scGPT} × {single, all-drugs}` on `split_ctrp`.
**Result:** scGPT overfits less (single-task gap 0.012–0.018 vs PCA 0.028–0.029), but with matched
capacity PCA is competitive/better on raw accuracy (all-drugs heads-beating: `hvg5000` 158 vs 156;
`all_genes` PCA 196 vs scGPT 137). The earlier scGPT lead (135 vs 103; 141 vs 80) was a capacity
artifact. Net: scGPT's edge is generalization, not predictive power.

### Drug coverage / learnability (`notebooks/drug_coverage.ipynb`)
Quantified per-drug coverage and response variance. No drug covers all 180 lines (max 179, median
171); 382 drugs ≥90% coverage, 80 drugs <50%; 14 drugs have std<0.05 (unlearnable). Low-coverage
drugs (~16 lines, n_val 221) are exactly the hardest/unreliable heads.

### 190 vs 180 cell-line overlap — resolved
**Not** a normalization difference (both rules give 190). **190** = SCP542 names in CTRPv2's
cell-line roster; **180** = those with actual post-QC viability measurements. 10 roster-listed but
unscreened lines drop out: `abc1, hs939t, jhh7, mdamb436, mfe280, ncih1048, ncih2073, ncih2347,
rerflckj, ten`. Use 180 (the trainable set).

## 26.05.2026

### Multi-drug refactor (v2-plan iterative step: single-task -> masked multi-task)
Replaced the paclitaxel-only training with a single multi-task entry point that
runs masked MSE over any subset of CTRPv2 drugs. K=1 reduces exactly to plain
MSE (no missing entries), so paclitaxel single-task is now just one flag combo
on `train_multitask.py` — the dedicated `train_baseline.py` / `train_scGPT.py`
scripts were removed.

**New artifacts written into the h5ad (by `ctrp_to_h5ad.py`):**
* `adata.obsm["Y_ctrp"]`  : float32 (n_cells, K), per-cell viability for each drug (NaN where missing).
* `adata.obsm["M_ctrp"]`  : bool    (n_cells, K), True where a label is observed.
* `adata.uns["ctrp_drugs"]`: ordered list of K drug names = column order of Y/M.
* `adata.obs["split_ctrp"]`: cell-line-grouped 70/15/15 split shared across every head (drug-agnostic, leakage-free).
* Per-drug legacy columns (`viability_<drug>`, `train_mask_<drug>`, `split_<drug>`) still written for downstream tools that read flat columns.

**Drug-scope knobs:**
* Default: `--min-cell-lines 50` -> keeps drugs screened on >=50 SCP542-overlapping cell lines.
* All CTRPv2 drugs: `--all-drugs` (equivalent to `--min-cell-lines 0`).
* Intermediate K: `train_multitask.py --drugs paclitaxel docetaxel gemcitabine ...` (no preprocessing rerun needed as long as those drugs passed the prep-time filter).

**Multi-task model & training:**
* `OncoMLP` takes `output_dim=K` (default 1).
* `training_utils.train_model` auto-switches to masked MSE / masked Huber when it sees 3-tuple `(x, y, mask)` batches; logs top-k best/worst per-drug val MSE.
* Only entry point: `scripts/training/train_multitask.py`.

**Sanity baseline (cheap failure detection, per v2 plan):**
* Per-drug-mean predictor (predicts the train-set mean viability per head) is computed up-front and compared to the trained model's per-drug val MSE at the end. Heads where the model fails to beat this floor have not learned anything useful regardless of absolute MSE.

### Data layout and path resolution

Defaults live in one place: `scripts/preprocessing/layout.py` (`DEFAULT_DATA_ROOT`, scGPT script/model paths). Override with `--data-root` / `--scgpt-script` if needed. No shell env vars required.

Directory structure under `DEFAULT_DATA_ROOT` (`/Users/selin/Desktop/OncoTox/data`):

```
data/
  scRNAseq_SCP542/expression/CPM_data.txt
  scRNAseq_SCP542/metadata/Metadata.txt
  metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/
  processed/scRNAseq_SCP542/hvg5000/     # default training (--variant hvg5000)
  processed/scRNAseq_SCP542/all_genes/   # full transcriptome (--variant all_genes)
```

`scripts/preprocessing/layout.py` is the **only** module that maps `(data_root, variant)` → file paths. The orchestrator passes those paths into each step. Rerunning `convert` or scGPT on the same variant **fails** unless `--overwrite` is passed; `hvg5000` and `all_genes` never share an output folder.

### Pipeline order

**A. Preprocessing (`run_preprocessing.py` orchestrates 5 steps):**
1. `scp542_conversion`       -> `processed/.../<variant>/SCP542_CCLE.h5ad`
2. `gen_embeds.py` (external) -> `..._scGPT_human_embeddings.h5ad`
3. `ctrp_to_h5ad`            -> `..._with_targets.h5ad`
4. `create_splits`           -> `split_paclitaxel` + `split_ctrp`
5. `add_pca`                 -> `X_pca`

Use `--variant hvg5000` (default) or `--variant all_genes` to pick the output subdirectory.

**B. Training (all via `train_multitask.py`):**

| Scenario | Command |
|---|---|
| scGPT, all CTRPv2 drugs that survived the prep filter | `train_multitask.py --use-rep X_scGPT` |
| PCA baseline, all CTRPv2 drugs | `train_multitask.py --use-rep X_pca` |
| scGPT, few-drug intermediate (validate masked loss on K=3 before going wide) | `train_multitask.py --use-rep X_scGPT --drugs paclitaxel docetaxel gemcitabine` |
| scGPT single-task paclitaxel (replaces old `train_scGPT.py`) | `train_multitask.py --use-rep X_scGPT --drugs paclitaxel` |
| PCA single-task paclitaxel (replaces old `train_baseline.py`) | `train_multitask.py --use-rep X_pca --drugs paclitaxel` |

Each training run writes a versioned directory under `runs/` (see below).

**C. End-to-end commands:**
```bash
uv run scripts/preprocessing/run_preprocessing.py \
    --variant hvg5000 --start-at targets --skip-scgpt --all-drugs

uv run scripts/preprocessing/run_preprocessing.py \
    --variant all_genes --all-drugs --scgpt-python "$SCGPT_PYTHON"

uv run scripts/training/train_multitask.py --use-rep X_scGPT
uv run scripts/training/train_multitask.py --variant hvg5000 --use-rep X_pca --drugs paclitaxel
```

### Run versioning + artifact saving
Every `train_multitask.py` run writes a self-contained directory under `runs/` (gitignored) and appends one row to `runs/runs_index.csv` for cross-run comparison. All versioning lives in `scripts/training/training_utils.py` (`create_run_dir`, `save_run`).

**Per-run files:**
* `config.json`           - the exact TrainConfig used.
* `run_meta.json`         - drug scope (`single_drug` / `subset` / `all_drugs`), rep, dataset sizes, hidden_dims, host info.
* `history.csv`           - epoch, train_mse, val_mse, lr.
* `summary.json`          - best_val_mse, best_epoch, mean baseline vs model MSE, heads-beating-baseline count.
* `best_model.pt`         - state_dict of the best-val-MSE checkpoint.
* `per_drug_results.csv`  - drug, model_val_mse, baseline_val_mse, delta, n_val.

**Ledger columns (`runs/runs_index.csv`):** `run_id, tag, scope, rep, K, n_train_cells, n_val_cells, best_epoch, best_val_mse, baseline_mean_mse, model_mean_mse, n_beats_baseline, n_total_heads, started_at, finished_at`.

Historical paclitaxel results (best val MSE 0.0375 / 0.0371 on 08.05, then 0.0351 / 0.0336 on 25.05 HVG-5000) are documented below in this file and serve as the single-task reference points for any future multi-task comparison.

### First multi-task results (HVG-5000, `--all-drugs`, `split_ctrp`)

Preprocessing was re-run with `--start-at targets --skip-scgpt --all-drugs` so that `Y_ctrp` / `M_ctrp` cover all 545 CTRPv2 drugs (no `--min-cell-lines` filter; 180 / 198 cell-line overlap). New `split_ctrp` distribution: train 34,126 / val 7,121 / test 5,980 / unassigned 6,286 cells over 126 / 27 / 27 cell lines. All four runs share these splits and use the per-drug-mean sanity baseline.

| Run id | Rep | K | Best epoch | Best val MSE | Baseline mean MSE | Model mean MSE | Heads beating baseline |
|---|---|---|---|---|---|---|---|
| `20260526_132914_multitask_X_scGPT_subset_K1` | X_scGPT | 1 (paclitaxel) | 11 | **0.0412** | 0.0434 | 0.0412 | 1 / 1 |
| `20260526_132952_multitask_X_pca_subset_K1`   | X_pca   | 1 (paclitaxel) |  5 | **0.0393** | 0.0434 | 0.0393 | 1 / 1 |
| `20260526_133012_multitask_X_scGPT_all_drugs` | X_scGPT | 545           |  7 | **0.0105** | 0.0097 | 0.0103 | **142 / 545** |
| `20260526_133112_multitask_X_pca_all_drugs`   | X_pca   | 545           |  6 | **0.0112** | 0.0097 | 0.0114 | 97 / 545 |

Observations:
* **Paclitaxel K=1 on `split_ctrp` (6,497 val labels) is not directly comparable to the historical 25.05 numbers (`split_paclitaxel`, 5,035 val labels):** the multi-task split holds out a different set of 27 cell lines, so the absolute MSE level shifts. Within this new split, PCA (0.0393) beats scGPT (0.0412) on paclitaxel alone, mirroring the previous-day trend.
* **K=545 numerically looks better (~0.0105) only because most viability values are near 1.0** — the per-drug-mean baseline drops to 0.0097, so the multi-task model only beats baseline on **142 / 545** heads (scGPT) vs **97 / 545** (PCA). scGPT beats baseline on ~47% more heads than PCA at the same K and split.
* The same low-coverage heads (n_val=221) dominate the "worse than baseline" list for both reps: `brd-k30748066`, `vx-680`, `brd-k33514849`, `brd9876:mk-1775 (4:1 mol/mol)`, `bafilomycin a1`. These are candidates for either dropping or weighting down in the next iteration.
* `gsk-j4` is the largest single win in both reps (model MSE ≈ 0.000 vs baseline 0.011 at n=221) — a sanity check that the multi-task head can fit a low-variance drug-line combination.

Saved artifacts:
* `runs/20260526_132914_multitask_X_scGPT_subset_K1/`
* `runs/20260526_132952_multitask_X_pca_subset_K1/`
* `runs/20260526_133012_multitask_X_scGPT_all_drugs/`
* `runs/20260526_133112_multitask_X_pca_all_drugs/`
* Aggregated ledger: `runs/runs_index.csv`.

### New notebook: HVG-5000 vs all-genes UMAP comparison
`notebooks/hvg_vs_all_genes_umap.ipynb` — uses `PipelinePaths.build` for `hvg5000` vs `all_genes`.

### Open questions for the next iteration
* Does multi-task hurt or help paclitaxel's val MSE vs the 0.0336 / 0.0351 single-task numbers? (Need to re-train paclitaxel single-task on `split_paclitaxel` to apples-to-apples compare; `split_ctrp` paclitaxel numbers above are on a different val set.)
* Which heads consistently fail to beat the per-drug-mean baseline? (Concrete starting list: small-n heads above — drop them or down-weight.)
* Loss weighting: currently uniform per observed entry; per-head weighting or per-drug uncertainty is a future tweak — especially relevant given the 142/545 vs 97/545 split.
* Does HVG-5000 leave noticeable cancer-type / drug-response signal on the table compared to the full transcriptome? (Visual answer from the new notebook once the all-genes side is regenerated.)

## 25.05.2026

### Preprocessing orchestrator + HVG-5000
Added `scripts/preprocessing/run_preprocessing.py` to run preprocessing end-to-end with tunable HVG count.

**Pipeline order:**
1. `scp542_conversion` → `SCP542_CCLE.h5ad`
2. `gen_embeds.py` (external: `/Users/selin/PycharmProjects/scGPT/gen_embeds.py`) → `..._scGPT_human_embeddings.h5ad`
3. `ctrp_to_h5ad` → `..._with_targets.h5ad`
4. `create_splits` (cell-line-grouped 70/15/15, writes `split_paclitaxel`)
5. `add_pca` (writes `X_pca`, optional `--force`)

**Command used:**
```bash
uv run scripts/preprocessing/run_preprocessing.py \
  --n-top-genes 5000 \
  --scgpt-python /Users/selin/PycharmProjects/scGPT/.venv/bin/python
```

**How HVG-5000 is filtered (`scp542_conversion.py`):**
1. Start from full CPM matrix (`53,513 × 22,722`).
2. Create copy, run `log1p`, then `sc.pp.highly_variable_genes(..., n_top_genes=5000, flavor="seurat")`.
3. Subset original CPM matrix to selected genes (saved `.X` remains CPM).
4. Save chosen count in `adata.uns["hvg_n_top_genes"]`.

**Pipeline outputs (HVG-5000):**
* Gene count: `22,722 → 5,000`
* scGPT vocab match: `4,576 / 5,000` (424 OOV)
* Embedded AnnData: `53,513 × 5,000`
* Paclitaxel labels: `44,367 / 53,513` cells
* Split counts: train `31,824` | val `5,035` | test `7,508` | unassigned `9,146`

### Training outputs kept from terminal runs

#### A) HVG-5000 with old model (BatchNorm+ReLU, hidden 64→32, fixed LR)
**PCA baseline (X_pca):**
* Epoch [01/50] | Train MSE: 0.1731 | Val MSE: 0.0547
* Epoch [05/50] | Train MSE: 0.0192 | Val MSE: 0.0362
* Epoch [10/50] | Train MSE: 0.0114 | Val MSE: 0.0373
* Epoch [15/50] | Train MSE: 0.0093 | Val MSE: 0.0377
* Epoch [20/50] | Train MSE: 0.0091 | Val MSE: 0.0388
* Epoch [25/50] | Train MSE: 0.0092 | Val MSE: 0.0379
* Epoch [30/50] | Train MSE: 0.0090 | Val MSE: 0.0381
* Epoch [35/50] | Train MSE: 0.0088 | Val MSE: 0.0395
* Epoch [40/50] | Train MSE: 0.0091 | Val MSE: 0.0368
* Epoch [45/50] | Train MSE: 0.0089 | Val MSE: 0.0358
* Epoch [50/50] | Train MSE: 0.0089 | Val MSE: 0.0393

**scGPT (X_scGPT):**
* Epoch [01/50] | Train MSE: 0.1187 | Val MSE: 0.0492
* Epoch [05/50] | Train MSE: 0.0232 | Val MSE: 0.0393
* Epoch [10/50] | Train MSE: 0.0179 | Val MSE: 0.0377
* Epoch [15/50] | Train MSE: 0.0169 | Val MSE: 0.0360
* Epoch [20/50] | Train MSE: 0.0173 | Val MSE: 0.0423
* Epoch [25/50] | Train MSE: 0.0175 | Val MSE: 0.0367
* Epoch [30/50] | Train MSE: 0.0176 | Val MSE: 0.0366
* Epoch [35/50] | Train MSE: 0.0177 | Val MSE: 0.0366
* Epoch [40/50] | Train MSE: 0.0177 | Val MSE: 0.0391
* Epoch [45/50] | Train MSE: 0.0175 | Val MSE: 0.0359
* Epoch [50/50] | Train MSE: 0.0177 | Val MSE: 0.0354

#### B) HVG-5000 with upgraded model/training (LayerNorm+GELU+scheduler+early stopping)
**PCA baseline (X_pca):**
* [baseline] Epoch [01/50] | Train MSE: 0.0895 | Val MSE: 0.0355 | LR: 1.0e-03  <- best
* [baseline] Epoch [05/50] | Train MSE: 0.0169 | Val MSE: 0.0363 | LR: 1.0e-03
* [baseline] Epoch [08/50] | Train MSE: 0.0136 | Val MSE: 0.0351 | LR: 5.0e-04  <- best
* [baseline] Epoch [10/50] | Train MSE: 0.0127 | Val MSE: 0.0374 | LR: 5.0e-04
* [baseline] Epoch [15/50] | Train MSE: 0.0115 | Val MSE: 0.0361 | LR: 2.5e-04
* [baseline] Early stopping at epoch 18 (no val improvement for 10 epochs).
* [baseline] Training complete. Best Val MSE: 0.0351 at epoch 8.

**scGPT (X_scGPT):**
* [scGPT] Epoch [01/50] | Train MSE: 0.0532 | Val MSE: 0.0363 | LR: 1.0e-03  <- best
* [scGPT] Epoch [03/50] | Train MSE: 0.0359 | Val MSE: 0.0350 | LR: 1.0e-03  <- best
* [scGPT] Epoch [05/50] | Train MSE: 0.0314 | Val MSE: 0.0343 | LR: 1.0e-03  <- best
* [scGPT] Epoch [10/50] | Train MSE: 0.0254 | Val MSE: 0.0350 | LR: 5.0e-04
* [scGPT] Epoch [12/50] | Train MSE: 0.0248 | Val MSE: 0.0341 | LR: 5.0e-04  <- best
* [scGPT] Epoch [14/50] | Train MSE: 0.0241 | Val MSE: 0.0336 | LR: 5.0e-04  <- best
* [scGPT] Epoch [15/50] | Train MSE: 0.0240 | Val MSE: 0.0401 | LR: 5.0e-04
* [scGPT] Epoch [20/50] | Train MSE: 0.0218 | Val MSE: 0.0357 | LR: 2.5e-04
* [scGPT] Early stopping at epoch 24 (no val improvement for 10 epochs).
* [scGPT] Training complete. Best Val MSE: 0.0336 at epoch 14.

### Model/training upgrades (concise)
* `scripts/model/OncoMLP.py`: default `LayerNorm`, `GELU`, input dropout (`0.1`), configurable `hidden_dims` (`(64,32)` PCA, `(128,64)` scGPT).
* `scripts/training/training_utils.py` (new): seeded training, ReduceLROnPlateau, early stopping, best-val checkpoint restore, gradient clipping.
* `train_baseline.py` / `train_scGPT.py`: migrated to shared `TrainConfig` + `train_model(...)`. (Both scripts were later removed on 26.05 -- subsumed by `train_multitask.py --drugs paclitaxel --use-rep X_pca|X_scGPT`.)

### Quick comparison (best val MSE)
| Setup | PCA best val | scGPT best val |
|-------|--------------|----------------|
| No HVG, regularized run from 08.05 | ~0.0375 (epoch 10) | ~0.0371 (epoch 10) |
| HVG-5000, old model | 0.0362 (epoch 5) | 0.0354 (epoch 50) |
| HVG-5000, upgraded model | **0.0351 (epoch 8)** | **0.0336 (epoch 14)** |

## 08.05.2026

### 1. Initial Setup & Random Splitting (The Data Leak)
* **Total cells:** 53,513
* **Cells with a valid paclitaxel viability score:** 44,367
* Added `ad.settings.allow_write_nullable_strings = True` to the split script to bypass AnnData string writing restrictions.
* Saved a **random 70/15/15 data split** directly into the `.h5ad` object (`split_paclitaxel` column).
    * **Train:** 31,056 cells | **Val:** 6,655 cells | **Test:** 6,656 cells | **Unassigned:** 9,146 cells

#### Training on Random Split (Data Leakage Identified)
* **scGPT Run:** Train MSE dropped to 0.0132, Val MSE to 0.0137.
* **PCA Baseline Run:** Train MSE dropped to 0.0022, Val MSE to 0.0011.
* *Conclusion:* The PCA baseline artificially outperformed scGPT because of data leakage. Since cells from the same cell line were randomly distributed across Train and Val, the PCA model simply memorized the tissue-of-origin "islands" rather than learning true biological resistance.

---

### 2. Fixing Data Leakage (Cell Line Grouped Split)
To prevent the model from cheating, the split logic was rewritten to group by cell line. If a cell line is in the training set, none of its cells appear in the validation set.
* Found 170 unique cell lines with paclitaxel labels.
* **Cell Line Split:** Train: 119 | Val: 25 | Test: 26
* **Final Cell Split Distribution:**
    * **Train:** 31,824 cells | **Val:** 5,035 cells | **Test:** 7,508 cells | **Unassigned:** 9,146 cells

#### Training on Grouped Split (Unregularized)
* **scGPT Run:** Train MSE 0.0110 | Val MSE: 0.0437
* **PCA Baseline Run:** Train MSE 0.0018 | Val MSE: 0.0390
* *Conclusion:* Data leak fixed. The PCA model completely failed to generalize to unseen cell lines (Val MSE stuck at ~0.04), proving it was previously cheating. However, scGPT also plateaued, indicating the highly parameterized MLP (256 hidden dims) was simply memorizing the noisy training labels.

---

### 3. Applying Aggressive Regularization (Proving the Hypothesis)
To force the MLP to learn generalized pathways instead of memorizing noise, aggressive regularization was applied:
* Reduced `hidden_dim` from 256 to 64.
* Increased `dropout_rate` from 0.3 to 0.5.
* Increased Adam optimizer `weight_decay` (L2 regularization) from 1e-5 to 1e-3.

#### Final Regularized Runs
**scGPT (Regularized)**
* Epoch [01/50] | Train MSE: 0.1014 | Val MSE: 0.0490
* Epoch [10/50] | Train MSE: 0.0260 | Val MSE: 0.0371
* Epoch [50/50] | Train MSE: 0.0260 | Val MSE: 0.0391
* *Gap:* ~0.013

**PCA Baseline (Regularized)**
* Epoch [01/50] | Train MSE: 0.0790 | Val MSE: 0.0480
* Epoch [10/50] | Train MSE: 0.0090 | Val MSE: 0.0375
* Epoch [50/50] | Train MSE: 0.0082 | Val MSE: 0.0380
* *Gap:* ~0.029

**Final Conclusion:** The foundation model successfully demonstrated a mathematically superior prior. While absolute validation performance hit a "weak supervision wall" around ~0.037 for both, the scGPT representation vastly reduced overfitting (Train/Val gap of 0.013 vs the baseline's 0.029). This proves standard PCA relies heavily on memorization, whereas scGPT forces the network to learn generalized, cross-tissue biological signatures.


## 21.04.2026 - 07.05.2026
### Preprocessing Timeline
scp542_conversion.py: massive, raw CPM_data.txt and Metadata.txt files and compiled them into foundational SCP542_CCLE.h5ad object.

scGPT Execution (The Embedding Phase): The scGPT foundation model was run using that raw .h5ad file to generate the continuous biological prior, outputting SCP542_CCLE_scGPT_human_embeddings.h5ad.

ctrp_to_h5ad.py (The Target Mapping): loaded those embeddings, parsed the CTRPv2 databases, merged the metadata, and cleanly mapped the paclitaxel viability scores to specific cell lines. This generated the final, all-in-one file: SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad.

The UMAP Script (Latent Space Validation): loaded that final file into memory, calculated the standard PCA baselines, ran the UMAPs, and generated the comparative visual proof

## 20.04.2026
### thoughts about available scores
#### GDSC2 (Genomics of Drug Sensitivity in Cancer)
Focus: Continuous dose-response profiling across a wide concentration range.
- LN_IC50 (Natural Log of Half-Maximal Inhibitory Concentration): The concentration required to kill 50% of the cells. Lower values mean higher toxicity to the cancer cell line.
- AUC (Area Under the Curve): The integrated area under the dose-response curve. It captures both the potency (how much drug is needed) and efficacy (the maximum kill rate). Lower values indicate strong sensitivity.
Best Pick for Multi-task: AUC is generally preferred over IC50 for machine learning because it captures the entire behavior of the curve, whereas IC50 can be noisy if the drug never actually reaches a 50% kill rate (resulting in extrapolated, artificial values).

#### CTRPv2 (Cancer Therapeutics Response Portal)
Focus: Similar to GDSC, this focuses on multi-dose response curves, but uses different curve-fitting algorithms and raw viability measurements.
- area_under_curve: As highlighted in your data dictionary, this is the integrated area under their sigmoid-fit concentration-response curve. Lower values equal higher sensitivity.
- cpd_avg_pv (Compound Average Percent Viability): The weighted average of surviving cells across all tested doses.
Best Pick for Multi-task: area_under_curve. Using AUC here allows your CTRPv2 head to learn a conceptually similar target to your GDSC head, even though the raw scales differ.

#### PRISM Repurposing (Public 24Q2)
Focus: High-throughput, single-dose screening. According to your provided Readme, all compounds here were screened at a single dose of 2.5 μM. Because there is no dose curve, there is no AUC or IC50.
- LMFI (Log2 Median Fluorescence Intensity): The raw optical readout of the barcodes.
- LFC (Log2 Fold Change): The median collapsed log-ratio of treated cells versus negative control (DMSO) cells. A negative LFC means the cancer cell line was depleted (killed) by the drug.
Best Pick for Multi-task: LFC (specifically from the Extended_Primary_Data_Matrix.csv). This is your only true efficacy metric for this dataset.

Multi-Task Recommendations
multi-task neural network with three output heads, your best target matrix configuration would be:

Head 1 (GDSC2): AUC
Head 2 (CTRPv2): area_under_curve
Head 3 (PRISM): LFC



hard to actually find toxicity definition and annotations, since toxicity is wanted for cancer treatments but excessive toxicity results in withdrawing from trials.

## 06.04.2026
### Ideas
- train on scRNA-seq data first and fine-tune on specific cancer types?
  - depends on available datasets
- train on scRNA-seq data first and fine-tune on clinical data?
  - depends on labels chosen for training?
  - depends on finding sufficient and properly labeled clinical data
- combine the two above ideas
  - find clinical scRNA-seq data sets for specific cancer types
  - maybe do bulk or pseudo bulk RNA-seq training as well, since scRNA-seq data is sparse?
- make a scRNA-seq transformer to use as a pretrained model? VAE? something else entirely?
- use scGPT embeddings (https://www.nature.com/articles/s41592-024-02201-0)
- look how (sc) DeepInsight could work with this


## 03.04.2026
### Advisor updates and alignment
- Clarification from Artem:
  - toxicity/response definition should be driven by available labels in selected datasets
  - if possible, include multiple response types instead of a single endpoint
  - prefer continuous outputs (IC50/viability-like) and binarize only for specific evaluations
  - prepare overlap and applicable sample counts to judge feasibility before modeling
- Methodology guidance:
  - DrugBank can be used for FDA/drug annotation support: https://go.drugbank.com/
  - multi-task setup can handle missing labels via masked losses (no need to force full intersection)
  - output/task weighting may be needed during training

### Notebook work completed (`notebooks/compare_GDSC_CTRP.ipynb`)
#### 1) Cell-line overlap with SCP542 (case-insensitive, unique names)
- SCP542 unique cell lines: `198`
- GDSC unique cell lines: `967`; overlap with SCP542: `133`; missing from GDSC: `65`
- CTRPv2 unique cell lines: `1107`; overlap with SCP542: `190`; missing from CTRPv2: `8`
- PRISM unique cell lines: `915`; overlap with SCP542: `182`; missing from PRISM: `16`
- Matching used normalized names (trim + lowercase; SCP542/PRISM names split on `_` where applicable).

#### 2) Drug/compound harmonization and cross-dataset overlap
- Built unified catalog: `data/drug/all_sources_drug_catalog.csv`
  - source rows: GDSC `295`, CTRPv2 `545`, PRISM `6575`
  - total rows: `7415`
- Name-based unique compounds:
  - GDSC `286`, CTRPv2 `545`, PRISM `6575`
  - overall union: `7040`
- Pairwise name overlap (normalized name):
  - CTRPv2 vs GDSC: `66`
  - CTRPv2 vs PRISM: `218`
  - GDSC vs PRISM: `144`
- Added CTRPv2<->PRISM BRD-based matching (from canonical `BRD-*` IDs):
  - BRD overlap: `243` (higher-confidence link)
- Exported overlap candidates:
  - `data/drug/drug_overlap_candidates.csv`

#### 3) DrugBank overlap
- Loaded DrugBank XML (`/Users/selin/Desktop/OncoTox/data/full database.xml`) and matched to catalog via normalized names (+ synonym expansion where available).
- Dataset-level matches to DrugBank:
  - GDSC: `118 / 295` (`40.00%`)
  - CTRPv2: `173 / 545` (`31.74%`)
  - PRISM: `3483 / 6575` (`52.97%`)
  - overall: `3774 / 7415` (`50.9%`)
- Exported review files:
  - `data/drug/drugbank_overlap_matches.csv`
  - `data/drug/drugbank_overlap_unmatched.csv`

#### 4) Applicable sample numbers (non-null response values)
- Computed for SCP542-overlapping cell lines vs each dataset's own full response space:
  - GDSC (`LN_IC50`): `8,007 / 242,036` (`3.31%`)
  - CTRPv2 (`cpd_avg_pv` viability): `1,521,028 / 7,227,951` (`21.04%`)
  - PRISM (Extended Primary matrix values): `1,210,432 / 4,213,048` (`28.73%`)
- Added second completeness metric in notebook:
  - `pct_non_null_within_overlap_subset` = non-null overlap rows / total overlap rows
  - Current values: GDSC `100%`, CTRPv2 `100%`, PRISM `97.95%`

### Interpretation notes
- "Applicable sample numbers" are comparable as response coverage indicators, but raw units differ:
  - GDSC/CTRPv2: long-table response rows
  - PRISM: matrix entries (compound x cell line)
- Therefore percentages are best interpreted as coverage estimates, not strict one-to-one sample equivalence.

### Cleanup/refactor performed
- Refactored repetitive SCP542 comparison code into reusable notebook helpers.
- Removed trailing empty notebook cells and tiny redundant summary CSV files.

### Pending decisions / next actions
- Decide primary training endpoint strategy:
  - option A: single dataset first (higher consistency)
  - option B: multi-task with masked missing labels (higher coverage)
- Decide whether to harmonize drug identity using:
  - conservative set (name + BRD agreement), or
  - broader set (name/synonym candidates for manual curation).
- Prepare a compact slide/table for in-person discussion with Artem next week:
  - cell-line overlap, response coverage, and recommended modeling path.

## 31.03.2026
### Advisor response to project-definition questions
- Project framing:
  - define toxicity according to available experimental labels in selected datasets
  - if multiple toxicity/response types exist and effort is reasonable, try to model all
- Target label format:
  - prefer continuous outputs (e.g., IC50/viability-like values) when available
  - convert model outputs to binary labels only when needed for downstream comparison
- Data source guidance:
  - Sanger site was temporarily unavailable at the time
  - scDrugAtlas can be used if all data are obtained successfully; note Harmony processing and avoid direct cross-dataset merging
  - can follow up in person at the lab

## 30.03.2026
- downloaded raw scRNA-seq cell line data from the PERCEPTION paper
  - original publication: https://doi.org/10.1038/s41588-020-00726-6
  - download link: https://singlecell.broadinstitute.org/single_cell/study/SCP542/%20pan-cancer-cell-line-heterogeneity#/
  - 53513 total cells, 22722 genes

-  CTRP^2 downloaded viability scores and annotations for chemical compounds
  - https://depmap.org/portal/data_page/?tab=allData
  - Cell Line ID: v20.meta.per_cell_line (master_ccl_id)
  - Experiments: v20.meta.per_experiment (master_ccl_id is inside)

- PRISM Repurposing for also failed drugs
  - https://depmap.org/portal/download/all/?release=PRISM+Primary+Repurposing+DepMap+Public+24Q2&file=Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv

- GDSC IC50 scores
  - https://cellmodelpassports.sanger.ac.uk/downloads

- clinical phase
- https://repo-hub.broadinstitute.org/repurposing#download-data

## 27.03.2026
### Clearing up project definition
Sent message to Artem to align on project direction before processing data. Waiting on replies for the following core questions:
- Project Focus (Toxicity Definition): Are we targeting Cytotoxicity/Efficacy (predicting if the drug successfully kills heterogeneous cancer cells) OR Adverse Patient Toxicity (predicting dangerous side effects to healthy tissue)?
- Target Variable Labeling: Should the model predict continuous values (e.g., IC50 scores) or binary categorical values (toxic/nontoxic, sensitive/resistant)?
- Dataset Priority: Should scDrugAtlas be the primary focus right now, or is there a better alternative?

## 26.03.2026
### Data Collection
#### Databases

- GDSC2
  - Drug Sensitivity Data -> GDSC2 IC50 Data & GDSC2 Raw Data
  - https://cellmodelpassports.sanger.ac.uk/downloads
  - Contacted DepMap/GDSC team
    - dead documentation link: https://depmap.sanger.ac.uk/documentation/gdsc/
    - asked specifically how the functional dataset “Drug Sensitivity Data” -> “GDSC2 IC50 Data” was processed
  - Status: waiting for response

- ClinTox
  - Binary toxicity prediction based on SMILES
  - https://tdcommons.ai/single_pred_tasks/tox/#clintox
  - maybe not well suited, only binary values

- scDrugAtlas
  - http://drug.hliulab.tech/scDrugAtlas/
  - Contacted Prof. Liu about:
    - original cell line IDs in the consolidated downloadable files
    - source publications for each dataset
    - how the datasets were integrated/combined
    - whether documentation exists for interpreting the consolidated h5ad files
    - whether there is a data dictionary / schema for the main variables and annotations
    - whether additional annotations exist (e.g. IC50 values or sample/cell line IDs) for matching to GDSC
  - Observation:
    - in `breast_cancer_palbociclib`, drug response appears to be encoded as binary values
  - Current issue:
    - hard to identify which cell line is which
    - consolidated file format / variable definitions unclear
  - Status: waiting for response
