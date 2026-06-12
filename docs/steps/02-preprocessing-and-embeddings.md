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
`(--data-root, --variant)` via `layout.py`, then runs **five steps in a fixed order** —
`STEP_ORDER = [convert → scgpt → targets → splits → pca]` — writing only under
`processed/scRNAseq_SCP542/<variant>/`. You normally never call the individual scripts by hand:
the gene-set variant is chosen with `--variant {hvg5000,all_genes}`, the drug scope with
`--all-drugs` / `--min-cell-lines`, `--start-at <step>` resumes mid-pipeline, `--skip-scgpt` reuses
existing embeddings, and `--overwrite` is required to replace the guarded `convert`/`scgpt` outputs
(everything is seeded via `--seed`, default 42).

### What each step reads and writes (in order)

| # | Step / script | Reads | Writes (added to the h5ad) |
|---|---|---|---|
| 1 | **convert** — `scp542_conversion.py` | `expression/CPM_data.txt` (genes×cells) + `metadata/Metadata.txt` | `SCP542_CCLE.h5ad`: cells×genes, `.X` = **CPM**. **HVG filtering happens here** (see below); records `uns["hvg_n_top_genes"]`. |
| 2 | **scgpt** — external `gen_embeds.py` (separate scGPT venv) | the **convert output** `SCP542_CCLE.h5ad` | `..._scGPT_human_embeddings.h5ad`: adds `obsm["X_scGPT"]` (**512-dim**). |
| 3 | **targets** — `ctrp_to_h5ad.py` | the embeddings h5ad + the 4 CTRPv2 tables | `..._with_targets.h5ad`: adds `obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]` ([Step 03](03-model-and-training-design.md) for mechanics). |
| 4 | **splits** — `create_splits.py` | the targets h5ad (in place) | `obs["split_paclitaxel"]` (`run`) + `obs["split_ctrp"]` (`run_multi`) — cell-line-grouped. |
| 5 | **pca** — `add_pca.py` | the targets h5ad (in place) | `obsm["X_pca"]`: `normalize_total(1e4)` → `log1p` → `sc.pp.pca` (≈50 comps) over the **same gene set**. |

✅ On-plan order: embeddings + comparative UMAP (below) come **before** any predictor — the plan's
Phase-1 latent validation gates the regression work.

---

## What HVG filtering removes — and why it forces re-embedding

**HVG filtering is part of step 1 (`convert`), not a later pass.** The gene-set choice is set by
`--variant` (`layout.VARIANT_N_TOP_GENES`: `hvg5000 → 5000`, `all_genes → None`). The raw `.X` is
**CPM** (counts-per-million: library-size-normalized but not log-transformed), so in
`scp542_conversion.py`:

- Highly variable genes are selected on a **`log1p` copy** via
  `sc.pp.highly_variable_genes(n_top_genes=5000, flavor="seurat")` — the Seurat flavor is a
  **dispersion-based** selector (it bins genes by mean expression and ranks each by its normalized
  dispersion = variance/mean within its bin), and the dispersion statistic assumes a log scale,
  which is why selection runs on the `log1p` copy. The **original CPM** matrix is then subset to the
  chosen genes (saved `.X` stays CPM, recording `uns["hvg_n_top_genes"]`); CPM is kept because scGPT
  applies its own value binning and PCA re-normalizes downstream, so neither wants a pre-logged `.X`.
- **What is filtered:** the **genes** — `22,722 → 5,000` for `hvg5000` (the most informative,
  cell-to-cell-variable genes). Cells (53,513) are untouched.
- **What is *not* the input to the model:** the genes / `.X` themselves are **never fed to the
  network**. They exist only to be turned into the two **representations** the model actually uses
  (next).

**The coupling that makes step order matter** — both downstream representations are derived from
whatever gene set step 1 kept:

```
convert (HVG choice fixes the gene set in .X)
        │
        ├─► scgpt  : X_scGPT  (512-d)  ← embeds exactly those genes
        └─► pca    : X_pca    (≈50-d)  ← PCA of exactly those genes
```

So **re-filtering = re-running `convert`, which invalidates both `X_scGPT` (must re-embed) and
`X_pca` (must recompute)**. Because `scgpt` reads the `convert` output, the fixed `convert → scgpt`
order already guarantees the embeddings reflect the current gene set — you cannot filter "after"
embedding. This is exactly why `hvg5000` and `all_genes` are **separate folders that never share
files** (`guard_output` enforces it): each is one self-consistent (gene set, X_scGPT, X_pca) triple.

### The two representations — what they are scientifically

- **`X_scGPT` (512-dim, the prior).** scGPT (`gen_embeds.py`, `scGPT_human` weights) is a
  transformer foundation model pretrained self-supervised on ~33 M human cells via masked
  gene-expression modeling. For each cell it tokenizes the expressed genes together with their
  value-binned expression and attends over them, emitting a single fixed-length **cell embedding**.
  Genes outside scGPT's gene vocabulary are **out-of-vocabulary (OOV)** and ignored — hence only
  4,576 / 5,000 HVGs contribute (424 OOV). This is the hypothesized *denoised biological prior* that
  aligns functional cell states across tissues.
- **`X_pca` (≈50-dim, the baseline).** A standard linear baseline: `add_pca.py` runs
  `normalize_total(target_sum=1e4)` → `log1p` → `sc.pp.pca`, capturing the directions of greatest
  linear variance. Because that variance is dominated by tissue-of-origin markers, PCA clusters
  cells into discrete lineage "islands" (the failure mode the scGPT prior is meant to overcome —
  Fig. 3/4 below).

### HVG-5000 pipeline outputs

- Genes: **22,722 → 5,000**
- scGPT vocab match: **4,576 / 5,000** (424 OOV — genes absent from the scGPT vocabulary)
- Embedded AnnData: 53,513 × 5,000 (with `X_scGPT` 512-d in `obsm`)
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

Re-running the **whole** orchestrator with `--variant all_genes` (HVG off) regenerates an
independent `(gene set = all 22,722, X_scGPT, X_pca)` triple under
`processed/scRNAseq_SCP542/all_genes/`. `notebooks/hvg_vs_all_genes_umap.ipynb` compares the two
variants' UMAPs. Evaluation of the all-genes side is still pending.

✅ On-plan / closes part of the HVG deviation by enabling the full-transcriptome comparison.

---

## Latent-space validation (UMAP, Fig. 3 / Fig. 4)

`notebooks/scgpt_umap.ipynb` is a **standalone validation notebook** (not part of the orchestrator):
it builds `X_umap_standard` (from `X_pca`) vs `X_umap_scGPT` (from `X_scGPT`) via `sc.pp.neighbors`
+ UMAP, colored by `Cancer_type` (**Fig. 3**) and `viability_paclitaxel` (**Fig. 4**). It visually
confirmed the hypothesis: PCA = discrete tissue "islands", scGPT = continuous shared manifold;
paclitaxel sensitivity mixed across the scGPT manifold.

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

Per variant, three h5ad files in pipeline order: `SCP542_CCLE.h5ad` →
`..._scGPT_human_embeddings.h5ad` → `..._scGPT_human_embeddings_with_targets.h5ad` (the trainable
file: `X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`, `split_ctrp`, `split_paclitaxel`).

**Reproduce:**
```bash
# From scratch (runs convert+HVG → embeddings → targets → splits → pca).
# The scgpt step needs the separate scGPT env, hence --scgpt-python.
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --scgpt-python /path/to/scgpt-venv/bin/python

# Re-derive only targets/splits/pca when convert + embeddings already exist.
# (convert/scgpt refuse to overwrite without --overwrite, so resume past them.)
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --start-at targets --skip-scgpt

# training
uv run scripts/training/train_multitask.py --use-rep X_scGPT            # all 545 drugs
uv run scripts/training/train_multitask.py --use-rep X_pca --drugs paclitaxel  # single-task PCA
```
