# OncoTox — Project Progress (index)

*Top-level index. The detailed, self-contained record is split into thematic step files under
[`docs/steps/`](./steps/) — this page holds the pipeline overview, the full project arc, the
current-status scorecard, and the doc-maintenance conventions. `project_notes.md` is a
complementary dated thought/decision log.*

Reference plan: `~/Desktop/OncoTox/project_plan/project_planning_v2.pdf`.
Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition** inside each step file.

> **Scope reality check — read this first.** Everything trained so far (Steps 04–05) uses **one
> database and one response score**: CTRPv2 `cpd_avg_pv` (viability). The 545-head "multi-task"
> run is multi-**drug**, *not* multi-database or multi-metric. The project's **ultimate goal is to
> combine all** — CTRPv2 + PRISM + GDSC, **efficacy and toxicity** — via cross-database masked
> multi-task ([Step 06](./steps/06-cross-database-integration.md)), then turn the result into a
> reusable foundation model fine-tunable on clinical outcomes
> ([Step 08](./steps/08-foundation-model-and-clinical-finetuning.md)). Don't read the current
> results as the finished goal.

---

## Pipeline overview (at a glance)

![OncoTox pipeline status overview](./pipeline_overview.png)

Green = done / on-plan · amber = addition or partial · red (dashed) = still missing.
Stages 1–6 are complete; the red boxes (cross-database PRISM/GDSC heads and the XAI stretch goal)
are the remaining work. Regenerate with `uv run docs/make_pipeline_overview.py`
(source: `docs/make_pipeline_overview.py`).

---

## The full project arc — document map

Each step is a self-contained file. **Steps 01–05 are done; 06–08 are placeholders** for planned
work, kept here so the entire project structure is visible end-to-end.

| Step | Status | What it covers |
|---|---|---|
| **[01 — Datasets & harmonization](./steps/01-datasets-and-harmonization.md)** | ✅ Done | Raw datasets (SCP542, CTRPv2, PRISM, GDSC), overlap/coverage audit, drug catalog, cell-line & compound harmonization. |
| **[02 — Preprocessing & embeddings](./steps/02-preprocessing-and-embeddings.md)** | ✅ Done | AnnData build, scGPT embeddings, UMAP latent validation, HVG-5000, `all_genes` variant, on-disk layout, reproduce commands. |
| **[03 — Model & training design](./steps/03-model-and-training-design.md)** | ✅ Done | Exact input/output/target/mask of a training example, MSE definition, **supervised** training paradigm. |
| **[04 — Single-task results](./steps/04-single-task-results.md)** | ✅ Done | Paclitaxel baseline + data-leak fix. **1 database, 1 score, 1 drug.** |
| **[05 — Multi-task results & versioning](./steps/05-multitask-results.md)** | ✅ Done | Masked-loss across 545 CTRPv2 drugs + run ledger. **Still 1 database, 1 score; multi-*drug* only.** |
| **[06 — Cross-database integration](./steps/06-cross-database-integration.md)** | ❌ Not started | **The "combine all" goal:** CTRPv2 + PRISM + GDSC, efficacy + toxicity, cross-database masked multi-task. |
| **[07 — XAI / feature interpretability](./steps/07-xai-feature-interpretability.md)** | ❌ Not started | Stretch goal: feature importance → transcriptomic drivers of resistance. |
| **[08 — Foundation model & clinical fine-tuning](./steps/08-foundation-model-and-clinical-finetuning.md)** | ❌ Not started | Overarching goal: reusable pan-cancer foundation model, fine-tunable on clinical (binary) outcomes. |

**Where this is going (the two axes that widen):**

```
Step 04   1 database · 1 score · 1 drug        (CTRPv2 cpd_avg_pv, paclitaxel)
Step 05   1 database · 1 score · K=545 drugs   (CTRPv2 cpd_avg_pv, all drugs)   ← here now
Step 06   3 databases · 2 metric types         (CTRPv2+PRISM+GDSC, efficacy+toxicity)
Step 08   + clinical fine-tuning               (continuous pre-train → binary clinical head)
```

**Fast facts you'll want regardless of which step you open** (full detail in
[Step 03](./steps/03-model-and-training-design.md)):

- A training example = **one single cell**; input = a 512-dim scGPT embedding (`X_scGPT`) or PCA
  (`X_pca`). Cell line / cancer type / drug are **not** input features.
- Target = CTRPv2 `cpd_avg_pv` viability, defined per **(cell line × drug)** and broadcast to every
  cell of that line — so MSE ≈ 0.01 is misleadingly tiny (values cluster near 1.0).
- Training is **fully supervised regression** (masked MSE/Huber). scGPT is a **frozen** self-supervised
  feature prior; the mask handles label sparsity but does **not** make it semi-supervised.

---

## Experiment matrix — PCA vs scGPT

The central comparison (the plan's core hypothesis) is run as a **2 × 2 × 2 = 8-run matrix**.
All eight runs share the cell-line-grouped split and the same MLP head + training protocol; only
the input representation (and its gene set) changes.

| Axis | Values |
|---|---|
| **Gene set** | `all_genes` (full transcriptome) · `hvg5000` (top-5,000 HVG from raw) |
| **Representation** | `X_pca` (standard single-cell PCA baseline) · `X_scGPT` (512-d embedding) |
| **Task** | single-task (paclitaxel) · multi-task (all drugs, K = 545) |

**Genes per condition** — PCA uses the full filtered set; scGPT uses only its in-vocabulary subset.
This OOV gap is **intentional** (scGPT's vocabulary coverage is part of the model — see
[Step 02](./steps/02-preprocessing-and-embeddings.md)):

| Gene set | PCA genes | scGPT genes (in-vocab) |
|---|---|---|
| `all_genes` | 22,722 | 20,570 |
| `hvg5000` | 5,000 | 4,576 |

Results live in [Step 04](./steps/04-single-task-results.md) (single-task) and
[Step 05](./steps/05-multitask-results.md) (multi-task). Action list for running the full matrix:
[TODO.md](./TODO.md).

---

## Where everything is saved (file map)

Data root: `DEFAULT_DATA_ROOT = /Users/selin/Desktop/OncoTox/data` (override with `--data-root`).
Each gene-set variant has its own folder `processed/scRNAseq_SCP542/<variant>/` holding three h5ad
files in pipeline order. **There is no separate file per representation, per drug, or per task** —
one trainable file per variant bundles everything, and the representation / drug / task are
*selected at training time*, not stored as separate files.

**Raw inputs** (shared, not per-variant):
- scRNA-seq counts → `data/scRNAseq_SCP542/expression/CPM_data.txt`
- cell metadata → `data/scRNAseq_SCP542/metadata/Metadata.txt`
- CTRPv2 tables → `data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.*`

**Exact location of each artifact** (paths under `processed/scRNAseq_SCP542/`):

| Artifact | HVG filtered? | File | Stored as | Shape / genes |
|---|---|---|---|---|
| Counts (CPM) | **filtered** | `hvg5000/SCP542_CCLE.h5ad` | `.X` | 53,513 × 5,000 |
| Counts (CPM) | **non-filtered** | `all_genes/SCP542_CCLE.h5ad` | `.X` | 53,513 × 22,722 |
| **scGPT embeddings — filtered HVG** | **filtered** | `hvg5000/SCP542_CCLE_scGPT_human_embeddings.h5ad` | `obsm["X_scGPT"]` | 53,513 × 512 (from 4,576 in-vocab genes) |
| **scGPT embeddings — non-filtered** | **non-filtered** | `all_genes/SCP542_CCLE_scGPT_human_embeddings.h5ad` | `obsm["X_scGPT"]` | 53,513 × 512 (from 20,570 in-vocab genes) |
| **PCA — filtered HVG** | **filtered** | `hvg5000/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad` | `obsm["X_pca"]` | 53,513 × 50 (computed on the 5,000 HVG) |
| **PCA — non-filtered** | **non-filtered** | `all_genes/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad` | `obsm["X_pca"]` | 53,513 × 50 (computed on all 22,722) |
| Drug labels — **all 545 drugs** | both | `<variant>/…_with_targets.h5ad` | `obsm["Y_ctrp"]` (+ `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`) | 53,513 × 545 |
| Drug labels — **one drug (paclitaxel)** | both | same `…_with_targets.h5ad` | one column of `Y_ctrp` selected via `--drugs paclitaxel`; legacy `obs["viability_paclitaxel"]` | 53,513 × 1 |
| Split — shared, cell-line-grouped | both | same `…_with_targets.h5ad` | `obs["split_ctrp"]` | per-cell |
| Split — paclitaxel-only (legacy) | both | same `…_with_targets.h5ad` | `obs["split_paclitaxel"]` | per-cell |

**The trainable file** is `<variant>/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad` — the only
file passed to training. It contains, together: `X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`, `split_ctrp`,
`split_paclitaxel`, `viability_paclitaxel`, and `uns["ctrp_drugs"]`.

**How the 8 matrix runs select from these files** (no new files are written for a run's inputs):
- gene set → `--variant {hvg5000, all_genes}` (which folder)
- representation → `--use-rep {X_scGPT, X_pca}` (which `obsm` key)
- task → `--drugs paclitaxel` (one drug) vs omitted (all 545)

**Training outputs:** each run writes `runs/<timestamp>_<tag>/` (gitignored) with `best_model.pt`,
`config.json`, `run_meta.json` (records the variant via the targets path), `history.csv`,
`summary.json`, `per_drug_results.csv`; one index row per run in `runs/runs_index.csv`.

---

## The plan (for reference)

A staged prototype (from the plan PDF):

1. **Latent-space validation** — generate scGPT embeddings, compare to full-transcriptome
   PCA via UMAP (Fig. 3 by cancer type, Fig. 4 by paclitaxel viability); confirm scGPT
   removes tissue-of-origin bias.
2. **Single-task baseline** — regress the continuous CTRPv2 `cpd_avg_pv` (viability)
   score from the embeddings on the **highest-confidence intersection** SCP542×CTRPv2
   (**190 cell lines, 545 compounds, 100 % non-null in overlap**). *Do not start
   multi-task / PRISM / GDSC until this works.*
3. **Iterate outward** — add masked-loss multi-task and integrate the larger, sparser
   PRISM (and GDSC) datasets — efficacy **and** toxicity.
4. **Stretch goal** — XAI / feature importance.

Overarching main goal: a reusable pan-cancer single-cell **foundation model** fine-tunable for
specific cancer types / clinical (binary) datasets.

**Core hypothesis:** scGPT embeddings are a denoised biological prior that forces the
regressor to learn real resistance signatures instead of memorizing cell line / tissue
identity → should show as **less overfitting (smaller train/val gap) for scGPT than PCA**.

---

## Current status — plan vs. reality

| Plan item | Status | Evidence |
|---|---|---|
| Sub-goal 1: compound harmonization (names + BRD + DrugBank) | ✅ Done | [Step 01](./steps/01-datasets-and-harmonization.md) |
| Sub-goal 2: masked-loss sparsity handling | ✅ Done (intra-CTRPv2) | [Step 05](./steps/05-multitask-results.md) |
| Sub-goal 3: baseline on SCP542×CTRPv2 highest-confidence intersection | ✅ Done | [Step 04](./steps/04-single-task-results.md)–[05](./steps/05-multitask-results.md) |
| Phase 1: scGPT embeddings + UMAP latent validation | ✅ Done | [Step 02](./steps/02-preprocessing-and-embeddings.md); Fig. 3/4 |
| Phase 2: single-task continuous `cpd_avg_pv` regression | ✅ Done | best scGPT val **0.0336** ([Step 04](./steps/04-single-task-results.md)) |
| Core hypothesis: scGPT overfits less than PCA | ✅ Confirmed | gap 0.013 vs 0.029 ([Step 04](./steps/04-single-task-results.md)) |
| Phase 3a: multi-task masked loss | ✅ Done **within CTRPv2 only** | [Step 05](./steps/05-multitask-results.md) |
| Phase 3b: integrate PRISM / GDSC (cross-database, efficacy+toxicity) | ❌ Not started | data downloaded + harmonized only ([Step 06](./steps/06-cross-database-integration.md)) |
| Stretch: XAI / feature importance | ❌ Not started | [Step 07](./steps/07-xai-feature-interpretability.md) |
| Main goal: foundation model + clinical fine-tuning | ❌ Not started (horizon) | [Step 08](./steps/08-foundation-model-and-clinical-finetuning.md) |

**Additions beyond the written plan (all defensible — document them):** random→leak→grouped
split ([Step 04](./steps/04-single-task-results.md)); HVG-5000 + all-genes comparison
([Step 02](./steps/02-preprocessing-and-embeddings.md)); per-drug-mean sanity baseline +
run-versioning ledger ([Step 05](./steps/05-multitask-results.md)).

**Two things to flag clearly in the writeup:**

1. **Multi-task today = 545 CTRPv2 drugs, not CTRPv2+PRISM+GDSC** — plan-Phase-3 half done
   (the real "combine all" is [Step 06](./steps/06-cross-database-integration.md)).
2. **Cell-line overlap is quoted as 190 (audit/Fig. 1) vs 180 (pipeline)** — same data,
   different name normalization; pick one and use it consistently.

---

## Open questions carried forward

*Action items live in [TODO.md](./TODO.md); this list is the scientific open questions.*

- Does multi-task help or hurt paclitaxel? Single-task on `split_ctrp` now exists (scGPT 0.0406,
  PCA 0.0372 on `hvg5000`); compare against the paclitaxel **head** inside the K=545 run
  ([Step 05](./steps/05-multitask-results.md)).
- Which low-coverage heads (n_val = 221) to drop or down-weight?
- Move loss from uniform-per-entry to per-head / uncertainty weighting?
- Does HVG-5000 lose signal vs the full transcriptome? (Pending the all-genes side of
  `hvg_vs_all_genes_umap.ipynb`.)
- When to integrate PRISM/GDSC as additional masked heads (the true Phase-3,
  [Step 06](./steps/06-cross-database-integration.md))?

---

## Maintaining these docs (conventions)

- **Source of truth = this index + [`steps/`](./steps/).** Together they must hold **all**
  important steps, numbers, parameters, and results, so everything is derivable from these files
  alone. `project_notes.md` is a dated thought/decision log — an *addition*, not the primary
  record; mine it for context but put authoritative numbers here.
- **When new work lands:** update the relevant **step file** with the hard details (cell/gene
  counts, split distributions, hyperparameters, run IDs, MSEs, deviations), then refresh the
  **scorecard + arc** above if the plan-vs-reality picture changed. Each step carries its own
  ✅ on-plan / ⚠️ deviation callouts against the plan PDF.
- **Keep numbers consistent** between this index, the step files, and `project_notes.md`.
- **Known cross-doc inconsistency to keep flagging:** SCP542×CTRPv2 cell-line overlap is **190**
  (audit notebook, case-insensitive) vs **180** (pipeline `ctrp_to_h5ad.py`, which also strips
  `-`). Same data, different normalization — pick one per context
  ([Step 01](./steps/01-datasets-and-harmonization.md)).
