# OncoTox notebooks

Figures and tables go to [`outputs/`](outputs/) (grouped by pipeline stage — see its README).
Model artifacts live in `runs/` (git-ignored), indexed by `runs/runs_index.csv`.

> ⚠️ **`07_training.ipynb` is superseded (13.07.2026).** Its 8-run matrix and its "ρ ≈ 0, the model
> cannot rank cell lines" conclusion were produced at K=545 on the legacy `mean_pv` target, whose
> **unstandardized per-drug variance was destroying the signal**. `10_diagnosis` reproduces that failure
> on demand and fixes it. Treat `07`'s conclusions as **pending a re-run**, not merely stale.

## Reading order

**`08 → 09 → 10 → 11 → 12`** is the current story, in the order it should be read:

| # | Notebook | Question it answers |
|---|---|---|
| **08** | `08_learnability_filter` | **Which drugs *can* be learned?** `learnability = min(#killed, #spared)`, coverage ≥ 90 % → **top 10 of 545**. *(A diagnostic simplification, not a result — see `10` §on validity: the score correlates only +0.36 with the ρ the model actually reaches.)* |
| **09** | `09_learnable5_training` | **Does the pipeline learn anything at all?** PCA vs scGPT on those 10, out-of-fold on ~150 held-out lines: **ρ = 0.360 / 0.396** — against ≈ 0 in June. |
| **10** | `10_diagnosis` | **Why did the 545-head model fail, and what fixes it?** The targets' biology, the implicit σ²-weighting of the loss, **the causal rescue test on the broken setting**, the model-knob ablations on the corrected one, and the ridge control. |
| **11** | `11_auc_vs_aucz` | **Which target?** `mean_pv` vs `auc` vs `auc_z`, at K=10 and K=545, with bootstrap CIs, Pearson, and seed stability. |
| **12** | `12_dreval_benchmark` | **How strong is this by the field's standard?** Our data + model through the real **DrEval** package (`drevalpy` 1.5.1): their LCO splits, their baselines, their metrics. |

**Supporting (not on the critical path):**

| # | Notebook | Role |
|---|---|---|
| 02 | `02_compare_GDSC_CTRP` | Drug-catalog harmonization (CTRP/GDSC/DrugBank) → writes `data/drug/*`. One-off. |
| 04 | `04_drug_coverage` | Coverage & response variance on the raw labels. Its *learnability* section is superseded by `08` (it was built on `mean_pv`); its target-distribution figures still stand. |
| 06 | `06_verify_variants` | Preprocessing QC + the PCA-vs-scGPT UMAPs. |
| 07 | `07_training` | The 8-run matrix / CV / HVG sweep. ⚠️ **Superseded** (see banner). |

### Why the numbering has gaps

**01, 03, 05 are missing on purpose** — they were moved to `archive/` and their numbers are *not*
reused: `07` and `04` are cited by name in the docs, the commit history and the progress reports, so
renumbering would break every reference to buy nothing. Treat the numbers as **stable identifiers, not
a sequence**.

| Archived | Why |
|---|---|
| `01_scDAExploration` | first look at the SCP542 data; superseded by everything after it |
| `03_analysis` | CTRP→PRISM repurposing map; never consumed by the pipeline |
| `05_preprocessing` | a front-end to the preprocessing CLI — duplicated by [Step 02](../docs/steps/02-preprocessing-and-embeddings.md#reproduce) |

`10_diagnosis` also absorbed a `13` that existed only briefly (it asked the same question as `10` on a
different setting); its number is retired.

The directories under [`outputs/`](outputs/) are **not** numbered, precisely because they do not map
one-to-one onto notebook numbers (`data/` is written by `04` *and* `10`; `target/` by `11` *and* `10`).
They are named for what they contain.

## Re-running

Every training notebook has a **`RETRAIN` flag, default `False`** — it then loads the saved CSVs from
`outputs/` and only redraws the figures (seconds). Set it to `True` to refit everything.

Fits use `TrainConfig(epochs=25)`: across 36 recorded runs the best epoch was **median 6, max 11**, and
early stopping (patience 10) never approached 25 — the previous cap of 50 only cost wall-clock.

Data is built by the CLI, not the notebooks:

```bash
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --score auc_z --start-at targets --skip-scgpt      # one targets h5ad per --score
```

All results in this repo use the **`hvg5000`** variant: 5,000 HVGs; scGPT embeds the **4,576** of them in
its vocabulary (424 OOV); PCA is computed on all 5,000 → both representations **512-d**.
