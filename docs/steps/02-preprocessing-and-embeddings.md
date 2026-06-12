# Step 02 — Preprocessing pipeline, embeddings & latent validation

*Part of [OncoTox project progress](../project_progress.md). Covers: building the AnnData,
generating scGPT embeddings, the UMAP latent-space validation, HVG-5000 filtering, the
orchestrator, the `all_genes` variant, and the on-disk data layout / reproduce commands.*

Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition**.

---

## Build the AnnData + latent-space validation (21.04–07.05.2026)

Pipeline scripts, in order:

1. **`scp542_conversion.py`** — read raw `CPM_data.txt` (genes × cells) + `Metadata.txt`,
   transpose to cells × genes, align metadata → **`SCP542_CCLE.h5ad`** (53,513 × 22,722,
   `.X` = CPM, no gene filter at this stage).
2. **scGPT embedding** — external `gen_embeds.py` (separate scGPT venv) using the
   `scGPT_human` weights → `SCP542_CCLE_scGPT_human_embeddings.h5ad` with
   **`obsm["X_scGPT"]` = 512-dim** per cell.
3. **`ctrp_to_h5ad.py`** — merge the 4 CTRPv2 tables, aggregate `cpd_avg_pv` to one value
   per (cell line, drug), map onto matching cells → `..._with_targets.h5ad`.
4. **UMAP validation** — `notebooks/scgpt_umap.ipynb`: standard PCA vs scGPT UMAP, colored
   by cancer type (**Fig. 3**) and by paclitaxel viability (**Fig. 4**).

✅ On-plan and in the right order: embeddings + comparative UMAP came before any predictor,
and visually confirmed the hypothesis (PCA = discrete tissue "islands", scGPT = continuous
shared manifold; paclitaxel sensitivity mixed on the scGPT manifold).

---

## Orchestrator + HVG-5000 filtering (25.05.2026)

**Orchestrator** `scripts/preprocessing/run_preprocessing.py` runs 5 steps in order:
`convert → scgpt → targets → splits → pca`. Paths derived once from `(data_root, variant)`
in `layout.py`; outputs live under `processed/scRNAseq_SCP542/<variant>/`. Expensive steps
refuse to overwrite without `--overwrite`; `hvg5000` and `all_genes` never share a folder.

**HVG-5000 filtering** (`scp542_conversion.py`): on a `log1p` **copy**, run
`sc.pp.highly_variable_genes(n_top_genes=5000, flavor="seurat")`, subset the **original CPM**
matrix to the selected genes (saved `.X` stays CPM), record `uns["hvg_n_top_genes"]=5000`.

**HVG-5000 pipeline outputs:**

- Genes: **22,722 → 5,000**
- scGPT vocab match: **4,576 / 5,000** (424 OOV)
- Embedded AnnData: 53,513 × 5,000
- Paclitaxel labels: 44,367 / 53,513 cells
- `split_paclitaxel`: train **31,824** / val **5,035** / test **7,508** / unassigned **9,146**

The model/training upgrade that landed alongside this work is documented in
[Step 03](03-model-and-training-design.md); the resulting single-task numbers are in
[Step 04](04-single-task-results.md).

> ⚠️ **Addition:** the plan only mentions full-transcriptome PCA; HVG-5000 (5,000-gene
> reduction) is a new variant — fewer scGPT OOV genes, smaller files. Justify it against
> the full-transcriptome path (the `all_genes` variant below exists for this comparison).

---

## `all_genes` (full-transcriptome) variant (26.05.2026)

The whole pipeline (convert → scGPT → targets) was regenerated without HVG filtering and now
exists on disk under `processed/scRNAseq_SCP542/all_genes/`.
`notebooks/hvg_vs_all_genes_umap.ipynb` compares HVG-5000 vs all-genes UMAPs. Evaluation of the
all-genes side is still pending.

✅ On-plan / closes part of the HVG deviation by enabling the full-transcriptome comparison.

---

## Current data layout (on disk)

`DEFAULT_DATA_ROOT = /Users/selin/Desktop/OncoTox/data`

```
data/
  scRNAseq_SCP542/expression/CPM_data.txt
  scRNAseq_SCP542/metadata/Metadata.txt
  metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/
  drug/                                  # harmonization catalogs + DrugBank exports
  processed/scRNAseq_SCP542/hvg5000/     # default training variant
  processed/scRNAseq_SCP542/all_genes/   # full transcriptome variant
```

Per variant, three h5ad files: `SCP542_CCLE.h5ad` → `..._scGPT_human_embeddings.h5ad`
→ `..._with_targets.h5ad` (the trainable file: `X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`,
`split_ctrp`, `split_paclitaxel`).

**Reproduce end-to-end:**
```bash
# preprocessing (HVG-5000, all CTRPv2 drugs, skip the external scGPT step if embeddings exist)
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --start-at targets --skip-scgpt --all-drugs
# training
uv run scripts/training/train_multitask.py --use-rep X_scGPT            # all 545 drugs
uv run scripts/training/train_multitask.py --use-rep X_pca --drugs paclitaxel  # single-task PCA
```

---

## Code, notebooks & key variables

**Orchestrator:** `scripts/preprocessing/run_preprocessing.py` — `STEP_ORDER = [convert, scgpt,
targets, splits, pca]`. Key flags: `--variant {hvg5000,all_genes}`, `--start-at <step>`,
`--skip-scgpt`, `--all-drugs`, `--min-cell-lines` (default 50), `--target-drug` (default
paclitaxel), `--force-pca`, `--overwrite`, `--n-top-genes`, `--seed`.

**Per-step scripts (in `STEP_ORDER`):**

| Step | Script | Does | Writes |
|---|---|---|---|
| convert | `scripts/preprocessing/scp542_conversion.py` | raw CPM + metadata → AnnData; optional HVG via `sc.pp.highly_variable_genes(n_top_genes, flavor="seurat")` on a `log1p` copy | `SCP542_CCLE.h5ad`, `uns["hvg_n_top_genes"]` |
| scgpt | external `gen_embeds.py` (separate scGPT venv) | scGPT_human weights → embeddings | `obsm["X_scGPT"]` (512-dim) |
| targets | `scripts/preprocessing/ctrp_to_h5ad.py` | map CTRPv2 viability onto cells ([Step 03](03-model-and-training-design.md) for mechanics) | `obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]` |
| splits | `scripts/preprocessing/create_splits.py` | cell-line-grouped splits | `obs["split_ctrp"]`, `obs["split_paclitaxel"]` |
| pca | `scripts/preprocessing/add_pca.py` | `sc.pp.normalize_total(1e4)` → `log1p` → `sc.pp.pca` (default 50 comps) | `obsm["X_pca"]` |

**Path layout:** `scripts/preprocessing/layout.py` — `PipelinePaths.build(data_root, variant)`,
`VARIANTS=("hvg5000","all_genes")`, filenames `SCP542_CCLE.h5ad` → `..._scGPT_human_embeddings.h5ad`
→ `..._scGPT_human_embeddings_with_targets.h5ad`.

**Notebooks:**
- `notebooks/scgpt_umap.ipynb` — **Fig. 3 / Fig. 4**: builds `X_umap_standard` (from `X_pca`) vs
  `X_umap_scGPT` (from `X_scGPT`) via `sc.pp.neighbors` + UMAP, colored by `Cancer_type` and
  `viability_paclitaxel`.
- `notebooks/hvg_vs_all_genes_umap.ipynb` — HVG-5000 vs `all_genes` UMAP comparison (`X_scGPT` /
  `X_pca` across both variants).

**Key variables produced:** `obsm["X_scGPT"]` (512), `obsm["X_pca"]` (≈50), `uns["hvg_n_top_genes"]`,
plus the target arrays consumed in [Step 03](03-model-and-training-design.md).
