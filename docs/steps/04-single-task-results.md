# Step 04 — Single-task results (paclitaxel) & the data-leak fix

*Part of [OncoTox project progress](../project_progress.md). Covers: the single-task paclitaxel
baseline, the random-split data leak and its grouped-split fix, regularization, and the
progression of best single-task val MSE through the HVG-5000 + model upgrade.*

This is plan-Phase-2 (single-task continuous `cpd_avg_pv` regression). Model/training design
is in [Step 03](03-model-and-training-design.md).

> **Scope — 1 database, 1 score.** Everything here predicts the **single** CTRPv2 metric
> `cpd_avg_pv` (viability) for the **single** drug paclitaxel. This is the narrowest slice of
> the project: one dataset, one response type, one compound. The widening happens in
> [Step 05](05-multitask-results.md) (still CTRPv2 only, across drugs) and ultimately in
> [Step 06](06-cross-database-integration.md) (the true "combine all databases + metrics" goal).

---

## Single-task paclitaxel baseline + data-leak fix (08.05.2026)

The first predictor regresses per-cell **paclitaxel viability** — the column
`obs["viability_paclitaxel"]`, i.e. the bulk `cpd_avg_pv` broadcast to every cell of the matching
line. It is loaded by `ScGPTDrugDataset` (`scripts/model/dataset.py`, `target_drug="paclitaxel"`)
and trained with `train_multitask.py --use-rep {X_scGPT|X_pca} --drugs paclitaxel`
(`output_dim = 1`). This was deliberately built **smallest-first** (plan §Prototyping) and used as
a methodological probe before scaling out.

- Total cells **53,513**; cells with a valid paclitaxel label **44,367**.

**Step 1 — random 70/15/15 split (the deliberate mistake):** `split_paclitaxel` →
train 31,056 / val 6,655 / test 6,656 / unassigned 9,146.

- scGPT: train MSE 0.0132 / val 0.0137
- PCA: train MSE 0.0022 / **val 0.0011** ← implausibly good

→ **Data leakage**: with cells split randomly, the same cell line lands in both train and val.
Since the label is constant within a line and PCA isolates each line as a tissue "island", the model
reduces to a nearest-neighbour lookup of the memorized per-line label — the val score measures
memorization, not generalization. The implausibly low PCA val MSE is the tell.

**Step 2 — cell-line-grouped split** — the correct cross-validation design for this label
structure: `create_splits.py` `run()` partitions **whole cell lines** with sklearn
`train_test_split` (`random_state=42`, group = `Cell_line`; 70/15/15 = test_size 0.30 then 0.50)
into `obs["split_paclitaxel"]`, so no line appears in two splits:

- **170 cell lines with paclitaxel labels → 119 train / 25 val / 26 test**
- Cells: **train 31,824 / val 5,035 / test 7,508 / unassigned 9,146**
- Re-trained unregularized (old 256-dim MLP): scGPT val 0.0437, PCA val 0.0390 → PCA's
  generalization collapsed, confirming prior cheating.

**Step 3 — aggressive regularization** (hidden 256→64, dropout 0.3→0.5, weight_decay
1e-5→1e-3):

| Model | Train MSE (ep 50) | Val MSE | Train/val gap |
|---|---|---|---|
| scGPT | 0.0260 | ~0.0371 (ep 10) | **≈ 0.013** |
| PCA | 0.0082 | ~0.0380 (ep 10) | **≈ 0.029** |

✅ On-plan + **core hypothesis confirmed**: near-equal val MSE, but scGPT overfits far less
— exactly the Fig. 4 prediction that PCA cheats by classifying cell line.

> ⚠️ **Addition (good practice):** the random-split → leak-diagnosis → grouped-split arc
> isn't written in the plan, but it's the "find failures cheaply, document even suboptimal
> versions" discipline the plan asks for. Worth keeping as a result.

---

## Best single-task val MSE — progression (with HVG-5000 + model upgrade, 25.05.2026)

After the HVG-5000 variant ([Step 02](02-preprocessing-and-embeddings.md)) and the
LayerNorm/GELU model + scheduler/early-stop training upgrade
([Step 03](03-model-and-training-design.md)):

**Paclitaxel single-task results — progression of best val MSE:**

| Setup | PCA best val | scGPT best val |
|---|---|---|
| No-HVG, regularized (08.05) | ~0.0375 (ep 10) | ~0.0371 (ep 10) |
| HVG-5000, old model (BatchNorm/ReLU) | 0.0362 (ep 5) | 0.0354 (ep 50) |
| **HVG-5000, upgraded model** | **0.0351 (ep 8)** | **0.0336 (ep 14)** |

These are the **single-task reference points** (`split_paclitaxel`, 5,035 val cells), loaded by
`ScGPTDrugDataset` over `obsm["X_scGPT"]` / `obsm["X_pca"]`.

✅ On-plan (still single-task CTRPv2 viability on the overlap; best result to date).

> ⚠️ **Not directly comparable to the multi-task K=1 numbers:** these use `split_paclitaxel`
> (25/26 held-out lines); the multi-task paclitaxel head in [Step 05](05-multitask-results.md)
> uses `split_ctrp` (27 held-out lines). An apples-to-apples "does multi-task help paclitaxel?"
> comparison still needs a single-task re-run on `split_ctrp`.
