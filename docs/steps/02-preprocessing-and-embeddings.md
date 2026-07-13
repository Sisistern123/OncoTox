# Step 02 — Preprocessing pipeline, embeddings & latent validation

*Part of [OncoTox project progress](../project_progress.md). Covers: the one orchestrator that
runs the whole pipeline, what each step reads/writes, exactly what HVG filtering removes and how it
couples to the embeddings, the UMAP latent-space validation, the `all_genes` variant, and the
on-disk layout / reproduce commands.*

Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition**.

---

## The orchestrator runs everything (`run_preprocessing.py`)

One entry point builds a complete, trainable h5ad from raw files:
`scripts/preprocessing/run_preprocessing.py`. It derives all paths once from
`(--data-root, --variant, --score)` via `layout.py`, then runs **five steps in a fixed order** —
`STEP_ORDER = [convert → scgpt → targets → splits → pca]` — writing only under
`processed/scRNAseq_SCP542/<variant>/`. You normally never call the individual scripts by hand:
the gene-set variant is chosen with `--variant {hvg5000,all_genes}`, the **response score** with
`--score {auc_z,auc,mean_pv}` ([Step 03](03-model-and-training-design.md)), the drug scope with
`--all-drugs` / `--min-cell-lines`, `--start-at <step>` resumes mid-pipeline, `--skip-scgpt` reuses
existing embeddings, and `--overwrite` is required to replace the guarded `convert`/`scgpt` outputs
(everything is seeded via `--seed`, default 42).

`--variant` and `--score` are the **two axes of the output layout**: the variant picks the folder,
the score picks the targets filename. Every script that touches a targets file takes both flags
(`layout.add_data_args`), so a `mean_pv` run and an `auc_z` run can coexist and be compared without
rebuilding the expensive `convert`/`scgpt` outputs, which they **share**.

### What each step reads and writes (in order)

| # | Step / script | Reads | Writes (added to the h5ad) |
|---|---|---|---|
| 1 | **convert** — `scp542_conversion.py` | `expression/CPM_data.txt` (genes×cells) + `metadata/Metadata.txt` | `SCP542_CCLE.h5ad`: cells×genes, `.X` = **CPM**. **HVG filtering happens here** (see below); records `uns["hvg_n_top_genes"]`. |
| 2 | **scgpt** — external `gen_embeds.py` (separate scGPT venv) | the **convert output** `SCP542_CCLE.h5ad` | `..._scGPT_human_embeddings.h5ad`: adds `obsm["X_scGPT"]` (**512-dim**) **and drops scGPT-OOV genes from `.X`** (hvg5000: 5,000→4,576). |
| 3 | **targets** — `ctrp_to_h5ad.py` | the embeddings h5ad + the CTRPv2 tables (curve fits for `auc`/`auc_z`, dose grid for `mean_pv`) | `..._with_targets[_<score>].h5ad`: adds `obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`, `uns["ctrp_score"]` + the de-standardization stats ([Step 03](03-model-and-training-design.md) for mechanics). |
| 4 | **splits** — `create_splits.py` | the targets h5ad (in place) | `obs["split_paclitaxel"]` (`run`) + `obs["split_ctrp"]` (`run_multi`) — cell-line-grouped. |
| 5 | **pca** — `add_pca.py` | the targets h5ad + the **convert counts** `SCP542_CCLE.h5ad` | `obsm["X_pca"]`: `normalize_total(1e4)` → `log1p` → `sc.pp.pca` (**512 comps**, matching the scGPT width) computed on the **HVG-filtered convert counts** (5,000 genes), *not* the targets `.X`. Targets `.X` left unchanged. |

✅ On-plan order: embeddings + comparative UMAP (below) come **before** any predictor — the plan's
Phase-1 latent validation gates the regression work.

---

## What HVG filtering removes, and what `.X` holds at each stage

HVG filtering happens **inside step 1 (`convert`)**, never as a later pass. Whether it runs is set
by `--variant` (`layout.VARIANT_N_TOP_GENES`: `hvg5000 → 5000`, `all_genes → None`).

**How the selection works** (`scp542_conversion.py`), starting from `.X` = CPM
(counts-per-million; 22,722 genes × 53,513 cells):

1. Copy the matrix and `log1p` the **copy only**. (`sc.pp.highly_variable_genes(flavor="seurat")`
   ranks genes by normalized dispersion — a statistic defined on the log scale — so it needs
   log-transformed input.)
2. On that copy, keep the **top 5,000 genes** by dispersion.
3. Subset the **original CPM** matrix to those 5,000 genes. Discard the log1p copy.

So at the `convert` step, two things are true:

- **Only genes are filtered** — `22,722 → 5,000`. All 53,513 cells are kept.
- **The values are not transformed.** `log1p` only *ranked* the genes; it never touched the saved
  numbers. The kept genes keep their CPM values, and `convert` records `uns["hvg_n_top_genes"]`.

**A second, scGPT-specific reduction happens at the `scgpt` step — and it does *not* propagate to
PCA.** scGPT can only embed genes in its own fixed vocabulary, so `gen_embeds.py` drops the
out-of-vocabulary (OOV) genes from the embeddings file's `.X`: `5,000 → 4,576` for `hvg5000`
(424 OOV) and `22,722 → 20,570` for `all_genes` (2,152 OOV). This shrinks **only the gene set scGPT
embeds**. The HVG filter is applied **once**, and PCA uses that full filtered set (below).

**What `.X` holds along the pipeline** (`hvg5000` gene counts shown):

| After step | `.X` holds | genes |
|---|---|---|
| convert | CPM, subset to HVG | 5,000 |
| scgpt | CPM, scGPT-OOV genes dropped (+ `obsm["X_scGPT"]`) | 4,576 |
| targets, splits | CPM, unchanged | 4,576 |
| pca | CPM, **unchanged** (+ `obsm["X_pca"]`, computed from the convert file) | 4,576 |

So the trainable file's `.X` stays CPM throughout (the `pca` step no longer rewrites it). The model
never reads `.X` anyway — only `obsm["X_scGPT"]` / `obsm["X_pca"]`.

**Filter once — where each representation's genes come from.**

```
convert : 22,722 → 5,000 genes (HVG)  — the single filter        [.X = CPM]
   ├─ scgpt : embeds the 4,576 of those in scGPT's vocabulary ──► X_scGPT (512-d)
   └─ pca   : PCA of all 5,000 HVG genes (read from convert) ───► X_pca   (512-d)
```

`add_pca.py` reads the **convert counts** `SCP542_CCLE.h5ad` (the full HVG set) to compute `X_pca`,
*not* the targets `.X` (which lost the OOV genes). So `X_pca` is a genuine HVG-5000 (or, for
`all_genes`, full-transcriptome) PCA — a standard single-cell PCA baseline — while scGPT uses the
vocabulary subset it is able to.

**This gene-count asymmetry (PCA on the full filtered set, scGPT on its in-vocab subset) is
intentional.** Dropping out-of-vocabulary genes is a real property of *using* scGPT, so each method
is compared exactly as it would be applied in practice: PCA on the genes the HVG step selects, scGPT
on the genes it can actually embed. The comparison therefore includes scGPT's vocabulary coverage as
part of the model, not as a confound to be normalized away.

Changing the gene set means re-running `convert`, which forces a re-embed and a re-PCA; that is why
`hvg5000` and `all_genes` live in **separate folders that never share files** (`guard_output`
enforces it). `notebooks/06_verify_variants.ipynb` checks these gene counts and the `X_pca` source at
any time.

### The two representations — what they are scientifically

- **`X_scGPT` (512-dim, the prior).** scGPT (`gen_embeds.py`, `scGPT_human` weights) is a transformer
  foundation model, pretrained self-supervised on ~33 M human cells. For each cell it reads the
  expressed genes and their binned expression values and outputs one fixed-length **cell embedding**.
  Genes outside scGPT's vocabulary are dropped as **out-of-vocabulary (OOV)**, so only 4,576 / 5,000
  HVGs contribute (424 OOV). This is the hypothesized *denoised biological prior* — it aligns
  functional cell states across tissues.
- **`X_pca` (512-dim, the baseline).** The standard single-cell linear baseline, sized to **512
  components to match the scGPT embedding width** so the two reps differ only in *how* the genes are
  encoded, not in input dimensionality (`add_pca.DEFAULT_N_COMPS = 512`, overridable with
  `--pca-n-comps`). `add_pca.py` runs
  `normalize_total(1e4)` → `log1p` → `sc.pp.pca` on the **full HVG-5000 convert counts** (`all_genes`:
  all 22,722), keeping the directions of greatest variance. That variance is dominated by
  tissue-of-origin markers, so PCA clusters cells into discrete lineage "islands" — the failure mode
  the scGPT prior is meant to overcome (Fig. 3/4 below).

### HVG-5000 pipeline outputs

- Genes after HVG (convert): **22,722 → 5,000** — the single filter
- scGPT embeds the **4,576** of those in its vocabulary (424 OOV); PCA uses all **5,000**
- Trainable AnnData: `.X` = CPM, **53,513 × 4,576** (OOV-dropped), carrying `X_scGPT` (512-d, from the
  4,576 vocab genes) and `X_pca` (512-d, from the 5,000 HVG genes) in `obsm`
- Paclitaxel labels: 44,367 / 53,513 cells
- `split_paclitaxel`: train **31,824** / val **5,035** / test **7,508** / unassigned **9,146**

The model/training upgrade that landed alongside this work is in
[Step 03](03-model-and-training-design.md); the single-task numbers are in
[Step 04](04-single-task-results.md).

> ⚠️ **Addition + history:** the first build (21.04–07.05.2026) used the **full transcriptome**
> (53,513 × 22,722, no HVG). HVG-5000 was added **inside `convert`** on 25.05.2026 — fewer scGPT
> OOV genes, smaller files. The plan only mentions full-transcriptome PCA, so HVG-5000 is a
> deviation justified against the full path via the `all_genes` variant below.

---

## `all_genes` (full-transcriptome) variant (26.05.2026)

Re-running the **whole** orchestrator with `--variant all_genes` (HVG off) regenerates an independent
gene set under `processed/scRNAseq_SCP542/all_genes/`. `convert` keeps all 22,722 genes; the `scgpt`
OOV-drop then leaves **20,570** in `.X` (what scGPT embeds), while **`X_pca` is computed on the full
22,722 convert counts** — a genuine full-transcriptome PCA. So the trainable file is **53,513 ×
20,570** in `.X`, carrying `X_scGPT` (from the 20,570 in-vocab genes) and `X_pca` (from all 22,722).
`notebooks/06_verify_variants.ipynb` checks the gene counts directly and plots the two variants'
UMAPs side by side. Evaluation of the all-genes side is still pending.

✅ On-plan / closes part of the HVG deviation by enabling the full-transcriptome comparison.

---

## Latent-space validation (UMAP, Fig. 3 / Fig. 4)

`notebooks/06_verify_variants.ipynb` (§7) is the **standalone validation** (not part of the
orchestrator): it builds PCA-vs-scGPT UMAPs for **both** variants via `sc.pp.neighbors` + UMAP,
colored by `Cancer_type` (**Fig. 3**) and `viability_paclitaxel` (**Fig. 4**). It visually confirmed
the hypothesis: PCA = discrete tissue "islands", scGPT = continuous shared manifold; paclitaxel
sensitivity mixed across the scGPT manifold.

---

## Current data layout (on disk)

`DEFAULT_DATA_ROOT = /Users/selin/Desktop/OncoTox/data`

```
data/
  scRNAseq_SCP542/expression/CPM_data.txt
  scRNAseq_SCP542/metadata/Metadata.txt
  metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/
  drug/                                  # harmonization catalogs + DrugBank exports (Step 01)
  processed/scRNAseq_SCP542/hvg5000/     # default training variant
  processed/scRNAseq_SCP542/all_genes/   # full transcriptome variant
```

Per variant, the two expensive files are built **once** and shared by every score:
`SCP542_CCLE.h5ad` → `..._scGPT_human_embeddings.h5ad`. The targets step then forks per score into
the trainable file (`X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`, `split_ctrp`, `split_paclitaxel`):

- `..._with_targets_auc_z.h5ad` — **default** (`--score auc_z`)
- `..._with_targets_auc.h5ad` — `--score auc`
- `..._with_targets.h5ad` — `--score mean_pv` (legacy name kept, so the Step 04–05 runs still resolve)

**Reproduce** (a documented, runnable walk-through of these commands lives in
`notebooks/05_preprocessing.ipynb`). `--score` defaults to `auc_z`; passing it explicitly is what
makes a comparison run unambiguous:
```bash
# From scratch (runs convert+HVG → embeddings → targets → splits → pca).
# The scgpt step needs the separate scGPT env, hence --scgpt-python.
# PCA width defaults to 512 (--pca-n-comps) to match the scGPT embedding.
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --score auc_z --scgpt-python /path/to/scgpt-venv/bin/python

# Add a second target score on top of existing convert+embeddings (this is the
# score-comparison build; convert/scgpt are untouched and reused).
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --score mean_pv --start-at targets --skip-scgpt

# Recompute only the 512-d PCA baseline in-place (what the 512-d switch needed).
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 \
    --start-at pca --skip-scgpt --force-pca --pca-n-comps 512

# training (--score selects which targets file to train on)
uv run scripts/training/train_multitask.py --use-rep X_scGPT --score auc_z   # all 545 drugs
uv run scripts/training/train_multitask.py --use-rep X_pca --score auc_z --drugs paclitaxel
```
