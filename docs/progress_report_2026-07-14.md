# Progress Report — Prediction of Anti-Cancer Drug Efficacy/Toxicity Scores
### 14 July 2026 · Selin Tuerkoglu

Slide text + figure references. Same structure as the previous decks.
Backing facts, parameters and CSV paths: `docs/progress_report_2026-07-14_notes.md`.

---

## Slide 1 — Title

Progress Report
Prediction of Anti-Cancer Drug Efficacy/Toxicity Scores
14th July 2026 · Selin Tuerkoglu

---

## Slide 2 — Our Approach

- **The Goal:** a pan-cancer single-cell foundation model capable of predicting pharmacological response
- **Current Phase:** June gave ρ ≈ 0 for both representations — this phase isolates the cause
  (target / evaluation / drug selection / model)

---

## Slide 3 — Core Hypothesis

**Overcoming Tissue Bias:**
- Standard feature extraction (like PCA) artificially clusters cells by tissue-of-origin (memorizing the
  cell line) rather than functional state
- Foundation model embeddings (like scGPT) project cells into a continuous, shared pan-cancer manifold,
  stripping away artificial tissue bias

---

## Slide 4 — Pipeline

**Figure:** `docs/pipeline_overview.png` (full width, no text)

---

## Slide 5 — Data

- Single-cell input: **SCP542** (Broad Single Cell Portal) — 53,513 cells × 22,722 genes, ~198 pan-cancer
  cell lines
- Labels: **CTRPv2** — 545 compounds, defined per **(cell line × drug)**
- Working set: 190 lines match CTRPv2's roster, **180 have post-QC viability** → the trainable set
- Bulk cell-line labels are broadcast to all of a line's cells
- **New:** we now read CTRPv2's **post-QC curve fits**, not only the raw dose measurements
  → three selectable targets (slide 9)

---

## Slide 6 — Preprocessing (unchanged since June)

- Gene-set variant used throughout: **5,000 HVG** (Seurat-style, on log1p; CPM kept)
- scGPT embeds in-vocab genes only; **OOV dropped: 5,000 → 4,576**
- PCA baseline at **512-d**, matched to the scGPT embedding width
- Fair-comparison choice: PCA on the full filtered set, scGPT on its in-vocab subset

---

## Slide 7 — Model Architecture & Task

**Figure:** `docs/model_architecture.png`

- 2 MLPs: **PCA vs scGPT**, both 512-d in, **same trunk** (128 → 64), one head per drug
- Output: **K = 545** heads (all drugs), or **K = 10** (learnable subset, slide 10)
- Masked MSE: `Σ(sq·M)/ΣM` — only observed (cell, drug) pairs contribute
- 74,954 params at K=10 · 109,729 at K=545

---

## Slide 8 — Training & Evaluation (what changed)

| | June | Now |
|---|---|---|
| Target | `mean_pv` (dose-averaged viability) | **`auc_z`** (per-drug z-scored AUC) |
| Evaluation | fixed val split, **27 cell lines** | **5-fold GroupKFold over 153 lines**, out-of-fold |
| Metric | per-drug Spearman/Pearson | same, but **averaged cell → cell line first** |
| Epochs | 50 | 25 (best epoch over 36 runs: median 6) |

- **Why the evaluation changed:** at n = 27, SE(ρ) ≈ **±0.2** — anything below ρ ≈ 0.4 was
  indistinguishable from zero. The June null was partly a *measurement* null.
- Label lives on the **line**, not the cell: 53,513 cells but **~150 independent examples per drug**
  (June: "effective n ≈ 126 lines"). Cells are pseudo-replicates → aggregate before correlating.

---

## Slide 9 — Action item: *"try out AUC"*

**Figure:** `outputs/data/target_biology.png` — all three targets on one CTRPv2 dose–response curve

| Target | Definition |
|---|---|
| `mean_pv` | unweighted mean of the 16 measured viability points — **the June target** |
| `auc` | area under the **fitted** sigmoid, normalized by the dose grid |
| `auc_z` | **`auc`, z-scored per drug** across cell lines |

- June expectation: *"AUC → more variance"*. In fact `mean_pv` and `auc` land within **~0.03** of each
  other on real curves.
- `auc_z` changes the **question**: not *"how potent is this drug?"* (= drug identity) but
  **"is this cell line more sensitive than average to this drug?"**

---

## Slide 10 — Action item: *"define learnability properly"*

**Figure:** `outputs/learnability/learnability_filter_auc.png`

June proposal: *"coverage ≥ threshold AND genuine response (≥ N lines with viability < 0.7), not just
std."* → implemented as one score:

```
learnability = min(#killed, #spared)     killed: auc ≤ 0.5   spared: auc ≥ 0.8
gate:          coverage ≥ 90 % of the 180 labelled lines      → top 10 drugs
```

- A drug is only rankable across lines if it **separates** them — kills a real set, spares a real set
- **This is a simplification for diagnosis, not a result.** It is *not* the claim that only 10 drugs are
  learnable (slide 14)
- ⚠️ The filter sees all 180 lines (incl. val/test) → **best-case subset**

---

## Slide 11 — Results: per-drug variance mis-scales the loss

**Figure:** `outputs/target/loss_weighting_bug.png`

- The masked MSE weights every (cell × drug) entry equally — but a drug's **squared error scales with
  σ²**
- σ across the 545 drugs: **0.034 – 0.302** → **≈ 80×** in squared error
- The widest **10 %** of drugs carry **31 %** of the loss; all 545 heads share **one trunk**

| | σ | min `auc` | max `auc` | lines killed | loss share |
|---|---|---|---|---|---|
| `fqi-2` | 0.296 | **0.59** | 2.11 | **0** | **1.18 %** |
| `dasatinib` | 0.155 | **0.07** | 1.09 | **35** | ~0.4 % |

- `fqi-2` **never drops below `auc = 0.59`** — it kills no cell line. Its variance is on the *other* side
  (up to 2.11 = lines growing faster under treatment): proliferation and assay noise.
- The **three largest loss carriers kill zero cell lines** between them. The 10 selected drugs carry ~2 %.

⇒ The loss weights each head by **σ² — by scale, not by learnability.**

---

## Slide 12 — Results: the target comparison

**Figure:** `outputs/target/target_comparison.png` (3 targets × K=10 / K=545, bootstrap CIs)

Out-of-fold ρ, **evaluated on the same 10 drugs in every row**:

| | `mean_pv` | `auc` | `auc_z` |
|---|---|---|---|
| K = 10, PCA | 0.370 | 0.360 | 0.360 |
| K = 10, scGPT | 0.396 | 0.405 | 0.396 |
| K = 545, PCA | 0.073 | −0.012 | **0.316** |
| K = 545, scGPT | −0.078 | −0.069 | **0.328** |

1. At **K = 10 all three targets are equivalent** — the curve fit alone buys nothing
2. At **K = 545 both unstandardized targets collapse**; scGPT even goes negative
3. **`auc_z` holds at K = 545**

- Per-drug z-scoring ≡ **weighting each head by 1/σ²** — the regression analogue of class-imbalance
  reweighting (Kendall et al. 2018). Changes **no per-drug metric**, only the loss.
- ⇒ This is the *"try reweighting to focus the model"* item — **and it is the one that works**
- ⇒ **The head count was never the problem.** Same 545 heads, same ~150 lines/drug, same model.

---

## Slide 13 — Results: what did *not* fix it

**Figure:** `outputs/ablations/rescue_k545.png`

Every remaining June item, applied to the **same failing config** (K=545, raw `auc`, scGPT;
baseline ρ = **−0.084**), one change at a time:

| June item | Change | ρ |
|---|---|---|
| *too much regularization?* | heavy (dropout 0.7, wd 1e-2) | −0.116 |
| *try reweighting* — as **sample** weights | line-balanced sampling | −0.116 |
| *shrink model?* | 74,954 → 16,810 params | −0.051 |
| *reduce batch size* | batch 32 | −0.008 |
| *too much regularization?* | **no regularization** | **+0.210** |
| *try reweighting* — as **task** weights | **`auc_z`** | **+0.333** |

- **Removing regularization recovers ~60 %** — the failure is a **capacity competition between heads**:
  under dropout 0.5 the trunk cannot serve both the noisy high-σ drugs and the learnable ones, and the
  noisy ones own the loss.
- But it treats the **symptom**: it memorizes the training lines (train MSE ≈ 0.01), reaches **half** of
  the weighting fix, and on the corrected loss the same regularization is **optimal again**
  (0.396 with / 0.372 without).
- On the corrected setting **all model knobs are flat** (ρ 0.32–0.40): regularization, capacity
  (74,954 → 2,560 params), batch size, sample reweighting.
  **Figure:** `outputs/ablations/ablation_reg_capacity.png`

---

## Slide 14 — Results: how many drugs are learnable?

Out-of-fold ρ for **all 537 scorable drugs**, K=545, `auc_z`
*(8 of 545 have < 20 labelled lines → no correlation computable)*

| | drugs | share | **June** |
|---|---|---|---|
| ρ > 0 | 408 | **76 %** | — |
| ρ > 0.2 | 170 | 32 % | — |
| ρ > 0.3 | 56 | **10 %** | **~4 %** |
| ρ > 0.4 | 12 | 2 % | — |
| median | | **0.12** | ≈ 0 |

The filter, validated against what the model actually achieves:
- it **enriches** — selected ρ = 0.40 vs rejected 0.12
- but it is a **weak ranker** — Spearman(learnability, achieved ρ) = **+0.36**
- **9 of the 12 drugs with ρ > 0.4 were rejected.** `ml210` (ρ = 0.516, best in the panel) failed the old
  coverage gate at 0.94 vs 0.95 → gate lowered to 0.90

⇒ **"10 / 545 pass the gates" ≠ "only 10 drugs are learnable."** The filter *understates* the model.

**What ρ means:** at ρ ≈ 0.45, picking the 15 lines the model calls most sensitive returns **≈ 4 truly
sensitive lines instead of 1.5 by chance** (≈ 3× enrichment).

---

## Slide 15 — Results: PCA vs scGPT

**Figure:** `outputs/learnability/pca_vs_scgpt.png`

| | ρ, out-of-fold |
|---|---|
| PCA | 0.360 |
| **scGPT** | **0.396** |

- Seed stability (K=545): PCA **0.332 / 0.295 / 0.299** · scGPT **0.328 / 0.344 / 0.385** — at seed 42 PCA is *ahead*
- June's *"within noise"* no longer holds — but the margin (+0.04…+0.08) is the size of the seed spread
  (±0.03) ⇒ **consistent, not large**

**Control — `RidgeCV` on the 150 cell-line mean embeddings (no deep learning, no single cells):**

| | ρ |
|---|---|
| Ridge (line-level), PCA | **0.342** |
| MLP, PCA | **0.356** |
| MLP, scGPT | **0.402** |

⇒ **Relative to a classical line-level regressor, the benefit comes from scGPT — not from the
single-cell resolution.**

---

## Slide 16 — External benchmark: DrEval

**Figure:** `outputs/dreval/dreval_lco.png` · notebook `12_dreval_benchmark.ipynb`

> Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response prediction
> models with DrEval*, **Nat. Commun. (2026)** · package `drevalpy` v1.5.1

Run through **their** package: **their** LCO splits (leave-cell-line-out = personalized medicine),
**their** baselines, **their** metrics. Nothing re-implemented.

**How the paper says results must be read:** normalized metrics **remove the mean drug and cell line
effects** from y_true *and* y_pred, *"because most of the explainable variation is driven by the drug
identity"*. In their study raw Pearson **0.91 → 0.56** per drug (Simpson's paradox).
**A model beats the baseline only in the normalized metrics.**

LCO, 5-fold, **normalized** (mean ± std), on our 10 drugs / 179 lines:

| | Spearman | R² |
|---|---|---|
| `NaiveMeanEffectsPredictor` — *the bar* | 0.000 | 0.000 |
| their `SingleDrugRandomForest` (our scGPT features) | 0.340 ± 0.065 | 0.099 |
| our `OncoMLP` (PCA) | 0.345 ± 0.050 | 0.087 |
| our `OncoMLP` (scGPT) | **0.357 ± 0.070** | **0.114** |

**Supported:** we clear `NaiveMeanEffects` — the bar **~half the published field fails**.
**Not supported:** we are at R² = 11 % vs **19 %** for their best LCO model; the scGPT–PCA gap
(0.357 vs 0.345) is **within fold spread**; we do **not** significantly beat their Random Forest.

---

## Slide 17 — Discussion & Limitations

- **The June null was substantially an artifact of a mis-scaled multi-task loss**, not a property of the
  data. June's *"the bottleneck is the label"* was right — but the mechanism was the label's **variance**,
  not its compression near 1.0. **No new data was needed.**
- **Model-side tuning is closed** — regularization, capacity, batch size, sample reweighting: all flat.
- **scGPT > PCA** in our metric, but **not** under DrEval's normalized metric ⇒ the core hypothesis is
  **supported, not proven**.
- **Ridge on line-means ties the PCA MLP** ⇒ the single-cell resolution pays off *only* with scGPT.
  This is the same result DrEval reports for the whole field.

**Limitations:**
1. **Drug selection saw val/test lines** → best-case subset. **Blocking.**
2. `learnability` correlates only **+0.36** with achieved ρ; under bootstrap the passing set is unstable
   (median 9 drugs, range 2–17)
3. DrEval numbers before 14.07 were inflated (test fold used for early stopping) — **corrected**
4. ~20 % of the DrEval signal is a pure **cell-line effect**; for `kx2-391` it is *all* of it

---

## Slide 18 — Next Steps

**Immediate**
- **Train-only drug selection inside each CV fold** — turns the best-case diagnosis into a
  generalization estimate *(blocking)*
- Re-run DrEval on **all 545 drugs**, and on **LTO / LDO**

**Model**
- **scDEAL** — still not attempted *(June: "if nothing works" — something now works, so it becomes a
  genuine comparison, not a fallback)*
- Learned task weights (Kendall) instead of fixed z-scoring
- Attention pooling per cell line *(June "Further Ideas")* — uses tumor heterogeneity as a feature

**Data**
- Cross-database PRISM + GDSC (masked multi-task)
