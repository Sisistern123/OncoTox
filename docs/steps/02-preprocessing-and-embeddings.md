# Step 02 — Preprocessing pipeline, embeddings & latent validation

*Part of [OncoTox project progress](../project_progress.md). Covers: the one orchestrator that
runs the whole pipeline, what each step reads/writes, exactly what HVG filtering removes and how it
couples to the embeddings, the UMAP latent-space validation, the `all_genes` variant, and the
on-disk layout / reproduce commands.*

Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition**.

---

## The steps, and who runs them (`pipeline.py`)

**Six steps in a fixed order** — `fetch → convert → scgpt → targets → splits → pca` — building a
complete, trainable h5ad from raw files. Everything is written under
`processed/scRNAseq_SCP542/<variant>/`, except `fetch`, which caches third-party data under
`metadata/`. All paths derive once from `(data_root, variant, score)` via `layout.py`, and everything
is seeded at 42.

Each step is one function in `scripts/preprocessing/pipeline.py`, and each owns the guard and the
preconditions belonging to it — `convert` and `scgpt` refuse to overwrite their outputs without
`overwrite=True`, `targets` refuses to run without the embedding and the response table, `pca`
refuses without `split_ctrp`. So a step cannot run out of order regardless of who calls it.

**The order itself is not in the code.** It is the numbering of the notebooks:
[`1_data`](../../notebooks/1_data.ipynb) runs `fetch` and `convert`;
[`3_representations`](../../notebooks/3_representations.ipynb) runs the other four, with
[`2_drug_selection`](../../notebooks/2_drug_selection.ipynb) between them because it needs the roster
and the response table and nothing else.

> ⚠️ **Rewritten 12.08.2026.** This section previously described `scripts/preprocessing/run_preprocessing.py`,
> a CLI holding the same six steps behind `STEP_ORDER` and `--start-at`, and said *"you normally never
> call the individual scripts by hand"*. That script was
> [archived](../../scripts/archive/README.md) when the notebooks were renumbered into five end-to-end
> stages: the ordering half became theirs, and a second copy of the order in a CLI is a second thing to
> keep in step with them. **There is no CLI replacement**, and the flags this section used to document
> — `--start-at`, `--skip-scgpt`, `--all-drugs`, `--overwrite` — no longer exist. Anything selective is
> now a direct call to the step you want, which is what one-function-per-step buys; see
> [Reproduce](#reproduce). *(It also said "five steps" while listing six.)*

`variant` and `score` are the **two axes of the output layout**: the variant picks the folder, the
score picks the targets filename. Both are carried on `PipelinePaths`, and the scripts that still have
CLIs take them as flags (`layout.add_data_args`), so an `auc_cc` run and an `ln_ic50_cc` run coexist
and can be compared without rebuilding the expensive `convert` / `scgpt` outputs, which they **share**.

Since 11.08.2026 the measure names carry a `_cc` suffix, for CurveCurator. That is deliberate rather
than cosmetic: the previous `auc` was a *different quantity* computed from CTRP's own distribution and
was defective, so reusing the name would have left old and new targets indistinguishable on disk under
identical filenames. `auc`, `auc_z` and `mean_pv` were all removed with their reader code
([Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026),
[corrections](corrections-and-dead-ends.md#the-auc-target-was-divided-by-the-wrong-quantity)).

### What each step reads and writes (in order)

| # | Step / script | Reads | Writes (added to the h5ad) |
|---|---|---|---|
| 0 | **fetch** — `fetch_ctrp_response.py` | Zenodo record **`21807175`**, pinned in `layout.ZENODO_RESPONSE_RECORD` | `metadata/drevalpy_CTRPv2_zenodo_21807175/`: the CTRPv2 response table and Cellosaurus 52.0, each **MD5-verified against the record's published checksums**, plus a `provenance.json` recording record, DOI and retrieval date. Idempotent — a cached archive matching its MD5 is neither re-downloaded nor re-extracted, so this costs a checksum pass. The **only** step that touches the network; `--start-at convert` skips it. |
| 1 | **convert** — `scp542_conversion.py` | `expression/CPM_data.txt` (genes×cells) + `metadata/Metadata.txt` | `SCP542_CCLE.h5ad`: cells×genes, `.X` = **CPM**. **HVG filtering happens here** (see below); records `uns["hvg_n_top_genes"]`. Since 05.08.2026 also `var["hgnc_symbol"]` — the current HGNC symbol per row, read only by step 2. |
| 2 | **scgpt** — `scripts/preprocessing/gen_embeds.py`, run under the separate scGPT venv via `--scgpt-python` (vendored into the repo 03.08.2026; it was previously an untracked file outside it) | the **convert output** `SCP542_CCLE.h5ad` | `..._scGPT_human_embeddings.h5ad`: adds `obsm["X_scGPT"]` (**512-dim**) **and drops scGPT-OOV genes from `.X`** (hvg5000: 5,000→4,576). |
| 3 | **targets** — `ctrp_to_h5ad.py` | the embeddings h5ad + the fetched response table `CTRPv2.csv` (one column per measure, both from the same CurveCurator fit) | `..._with_targets_<score>.h5ad`: adds `obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`, `uns["ctrp_score"]`, and `obs["cellosaurus_id"]` — the persistent accession per cell line, **recorded but never joined on**. Exact duplicate rows are dropped and a disagreeing pair **raises**. Resulting matrices, counts and every discrepancy against CTRPv2's own distribution: [Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026). |
| 4 | **splits** — `create_splits.py` | the targets h5ad (in place) | `obs["split_paclitaxel"]` (`run`) + `obs["split_ctrp"]` (`run_multi`) — cell-line-grouped. |
| 5 | **pca** — `add_pca.py` | the targets h5ad + the **convert counts** `SCP542_CCLE.h5ad` | `obsm["X_pca"]`: `log2(1 + CPM/10)` → `sc.pp.scale(max_value=10)` → `sc.pp.pca` (**512 comps**, matching the scGPT width) computed on the **HVG-filtered convert counts** (5,000 genes), *not* the targets `.X` — plus one train-fitted `X_pca_train_<split>` per fixed split, and a `uns["pca_fits"]` record per key (both below). Targets `.X` left unchanged. |

✅ On-plan order: embeddings + comparative UMAP (below) come **before** any predictor — the plan's
Phase-1 latent validation gates the regression work.

---

## What HVG filtering removes, and what `.X` holds at each stage

HVG filtering happens **inside step 1 (`convert`)**, never as a later pass. Whether it runs is set
by `--variant` (`layout.VARIANT_N_TOP_GENES`: `hvg5000 → 5000`, `all_genes → None`).

**How the selection works** (`scp542_conversion.py`), starting from `.X` = CPM
(counts-per-million; 22,722 genes × 53,513 cells):

1. Copy the matrix and transform the **copy only** to `log2(1 + CPM/10)`
   (`expression.kinker_transform`). `sc.pp.highly_variable_genes(flavor="seurat")` ranks genes by
   normalized dispersion — a statistic defined on the log scale — so it needs log-transformed
   input, and it reads `uns["log1p"]["base"]` to invert the transform before computing
   dispersions. *(Was a plain `log1p(CPM)` until 05.08.2026; see
   [the transform](#the-expression-transform-is-the-datasets-own-05082026).)*
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
vocabulary subset it is able to. *(The **transform** applied before that PCA was not standard until
05.08.2026 — see [What transform PCA sees](#what-transform-pca-sees--corrected-05082026) below.)*

**This gene-count asymmetry (PCA on the full filtered set, scGPT on its in-vocab subset) was intended
as part of the model.** Dropping genes scGPT cannot embed is a real property of *using* it, so each
method would be compared as it is applied in practice: PCA on the genes the HVG step selects, scGPT on
the genes it can embed.

> ⛔ **That justification does not hold — corrected 05.08.2026.** Most of the drop is **not** scGPT's
> vocabulary coverage; it is a symbol-matching defect. `gen_embeds.py` matches SCP542's symbols exactly,
> against an older HGNC annotation than the vocabulary uses, so **775 genes carrying 3.6 % of every
> cell's transcriptome** are discarded despite being in the vocabulary under their current names —
> RACK1, ATP5F1E, H2AZ1, the H3-3A/B and H4C3 histones among them. Only 0.21 % of the loss is genuine
> vocabulary coverage. Full write-up, scale and direction of the bias:
> [Corrections](corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols).
>
> **Repair applied to the code 05.08.2026 — `scripts/annotation/gene_symbols.py`.**
> `scp542_conversion.py` adds a **`var['hgnc_symbol']`** column carrying the current HGNC symbol —
> `var_names` keeps SCP542's identifiers exactly as distributed — and
> `gen_embeds.py::resolve_gene_names` resolves the vocabulary through it, **consulting it only where a
> row's own symbol is not a token**, so no gene embedded today can be lost. The **12 collisions**, where
> a recovered symbol already exists as its own row, are **left unmapped**: they carry 0.0021 % of the
> transcriptome, and merging them would change expression values and risk folding an antisense lncRNA
> into its sense gene. Renames are also refused where the old symbol is itself a currently approved
> symbol of another gene, since HGNC records reassignments (`OSR1`→`OXSR1`) alongside true renames.
> Every choice runs in the same direction — **no expression value changes**, so `X_pca` and the HVG
> selection stay bit-identical and the re-embed is attributable to the recovered genes alone. The
> three decisions and what each rejects are in
> [Corrections](corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols).
>
> ⚠️ **Code only — nothing is recomputed.** Every embeddings file on disk, and therefore every scGPT
> number in these docs and in the report, still carries the defect in full. The repair takes effect at
> the [clean sweep](../TODO.md); it landed before the sweep so the sweep is not paid for twice. The
> gene counts stated throughout this file (4,576 / 20,570) are what is **on disk**; after the sweep
> they become **4,704 / 21,332**.

Changing the gene set means re-running `convert`, which forces a re-embed and a re-PCA; that is why
`hvg5000` and `all_genes` live in **separate folders that never share files** (`guard_output`
enforces it). `notebooks/analysis/qc/verify_variants.ipynb` checks these gene counts and the `X_pca` source at
any time.

### The expression transform is the dataset's own (05.08.2026)

**`.X` really is raw CPM.** Verified against the distributed file: `CPM_data.txt` carries values such
as `31.54` and `42.24`, far outside the 0–15 range a log-transformed matrix would give. So the portal
ships the CPM matrix, not the `E` matrix the paper analyses.

**But Kinker et al. do not analyse raw CPM.** Their Methods ("Processing of scRNA-seq data") quantify
expression as

> `E[i,j] = log2(1 + CPM[i,j]/10)`

and state the reason for the divisor: the average number of UMIs detected per cell is **under
100,000**, so without dividing by 10 the difference between detected (`E > 0`) and undetected
(`E = 0`) genes is inflated. — Kinker et al., *Pan-cancer single-cell RNA-seq identifies recurring
programs of cellular heterogeneity*, **Nature Genetics 52, 1208–1218 (2020)**,
doi:10.1038/s41588-020-00726-6. Cited in `references.bib` as `scp542`.

**Why transform at all, rather than feeding CPM straight in.** PCA — and every method built on
Euclidean distance — works on *squared differences*, so each gene contributes in proportion to its
absolute spread. On raw CPM a housekeeping gene sitting at 5,000 CPM that varies by a mild 20% moves
cells 1,000 units apart, while a transcription factor at 50 CPM that varies by a full 100% moves them
50. The leading components then describe which cells expressed a lot of ribosomal and mitochondrial
RNA, not which cells differ biologically. The log fixes the mismatch at its source: expression
differences are *multiplicative*, and a log turns a fold-change into a fixed additive distance
regardless of where on the scale it happens, so a doubling counts the same in a rare gene as in an
abundant one.

Per-gene scaling (below) equalizes variance too, but it does not replace the log and cannot: gene
selection happens *before* scaling and its dispersion statistic is defined on the log scale, and
z-scoring a heavily right-skewed raw-CPM distribution leaves a handful of extreme cells dominating the
gene even after standardization. Log first, then scale.

**We now use their transform**, in both places that need it, through the single entry point
`scripts/preprocessing/expression.py::kinker_transform` so the two cannot drift apart. Until
05.08.2026 both used a plain `log1p(CPM)` — natural log, no divisor — which has no justification
specific to this data, while the authors' choice is argued from a measured property of it. The base
is set through `sc.pp.log1p(base=2)` rather than by hand, because that records
`uns["log1p"]["base"]`, which `highly_variable_genes(flavor="seurat")` reads to invert the transform
correctly.

⚠️ **This can change which genes are selected.** The earlier claim that HVG selection is invariant to
a global constant relied on gene means being ≫ 1, where `log1p(mean)` behaves like a shift and the
equal-width mean-bins move with it. Dividing by 10 pushes many means below 1, where that no longer
holds, so bin membership — and therefore the HVG set — may genuinely differ. Not measurable until the
pipeline is re-run.

### What transform PCA sees — corrected 05.08.2026

`add_pca.py::run` reads the convert file and applies this before `sc.pp.pca`:

```python
kinker_transform(src)            # log2(1 + CPM/10); .X is CPM, already library-size normalized
sc.pp.scale(src, max_value=10)   # per-gene z-score across cells, clipped
sc.pp.pca(src, n_comps=512, random_state=42)
```

**The standard this follows.** Scanpy's and Seurat's recipes both normalize library size **once, on
the full gene matrix, before gene selection**, then log-transform, then standardize genes, then take
components. SCP542 arrives as CPM — Kinker et al. distribute `CPM_data.txt` — so that normalization
has already happened and must not be repeated. The log step is the dataset authors' (above) rather
than a generic `log1p`. `max_value=10` follows Seurat's `ScaleData(scale.max = 10)` default and
scanpy's own tutorial setting; with 53,513 cells a single outlier can otherwise hand one gene an
entire component.

**Two defects this replaces**, both found in the 05.08.2026 audit of every line in the repository
that transforms expression values:

1. **Over-normalization.** `sc.pp.normalize_total(src, target_sum=1e4)` ran *after* HVG subsetting.
   That is not a rescale to a different target but a **second** library-size normalization computed
   over the retained genes only: each cell divided by its own HVG sum, so a cell whose expression
   sits largely outside the selected set was inflated relative to one whose expression sits inside
   it — a biological property turned into a scale factor. It also made variants mutually
   incomparable (the same cell scaled differently in `hvg1000` than in `hvg5000`) and the pipeline
   internally inconsistent: genes *ranked* on `log1p(CPM)`, then *projected* on
   `log1p(CP10K-over-HVGs)`.
2. **Under-normalization.** No per-gene standardization ran before PCA. CPM normalizes *within a
   cell, across genes*; it says nothing about a gene's variance *across cells*, so the leading
   components tracked absolute expression level rather than variation between cells. That is the
   axis `sc.pp.scale` addresses, and it is why Seurat and scanpy both scale before PCA whatever
   normalization preceded it.

**What it costs — and what was already being paid.** `sc.pp.scale` fits each gene's mean and standard
deviation over **all** 53,513 cells, held-out lines included, so a test cell is centred using
statistics it helped compute. That is leakage, but it is *not* new: `sc.pp.pca` already fits the
entire 512-dimensional rotation on all cells, and always has, which is a far larger all-cells fit
than two numbers per gene. Adding scaling makes an existing leak marginally wider; it does not open
one. HVG selection (`convert`) is the third, and applies to both arms.

**The asymmetry this creates between the arms matters more than the leak itself.** All three fits are
unsupervised — none sees a response label — which is why standard pipelines accept them. But `X_pca`
is fitted *on our cells* (HVG + scaling + rotation), while `X_scGPT` comes from frozen pretrained
weights and a binning that uses quantiles of each cell's **own** values, so a cell's scGPT embedding
draws on no other cell at all. Only the HVG gene set is shared. **The PCA baseline therefore gets two
all-cells fits that scGPT structurally cannot have, and the bias runs toward the control**: any
scGPT-over-PCA margin measured this way is conservative. This qualifies every PCA-vs-scGPT number in
the project and belongs next to review item 4's input-scale asymmetry in `docs/TODO.md`.

**What was done about it (05.08.2026).** For the **fixed splits** the leak is removed. `add_pca.py`
writes, alongside the all-cells `X_pca`, one train-fitted key per fixed-split column:

| key | fitted on | used for |
|---|---|---|
| `X_pca` | all 53,513 cells | UMAPs and latent-space validation — **and every CV run**, see below |
| `X_pca_train_ctrp` | cells of `split_ctrp == "train"` | models scored on `split_ctrp` |
| `X_pca_train_paclitaxel` | cells of `split_paclitaxel == "train"` | models scored on `split_paclitaxel` |

Both the per-gene mean/sd and the rotation are fitted on training cells and then applied to every cell
(`add_pca.py::_pca_fitted_on_train`, sklearn rather than scanpy because neither `sc.pp.scale` nor
`sc.pp.pca` separates fitting from transforming). `scripts/model/dataset.py::resolve_rep` maps
`--use-rep X_pca` to the train-fitted key whenever the run uses a fixed split, and **raises** if that
key is absent rather than falling back to the leaky matrix. The key actually read is recorded in each
run's `run_meta.json` as `rep_key`.

**Still not fixed: cross-validation.** CV folds are drawn at training time (`cv.py`, GroupKFold), so
`resolve_rep` leaves them on the all-cells `X_pca` — five fold-specific matrices cannot be one stored
array. HVG selection also remains an all-cells step, for both arms. So the asymmetry above is reduced,
not eliminated, and every CV number still carries it.

**"All cells" includes cells that never train (noted 10.08.2026).** `convert` runs before the CTRP
join, so the gene selection, `sc.pp.scale`'s per-gene mean/sd and the all-cells rotation are all fitted
on the **full** 53,513 cells — including those of the SCP542 cell lines that carry no CTRPv2 label at
all and are dropped from every split by `create_splits.py::has_any_label`. That is **18 lines /
6,286 cells** in the files on disk, and **17 lines / 6,073 cells** after the `H292` fix takes effect
([Step 01](01-datasets-and-harmonization.md#the-join-dropped-a-screened-cell-line-h292-10082026)) —
11.4 % of the atlas. It is not a test leak: those cells appear in no split, so no held-out label is
touched. But the representation is partly shaped by data the model never sees, which is a separate
question from the leak above and is decided with it under review item 7.

**Splits are frozen to a file (05.08.2026) — but the file does not exist yet (found 12.08.2026).**
`create_splits.py::frozen_split` reads `splits/split_ctrp.csv` (repo-versioned by design, not under the
gitignored data root) and only redraws under `--regenerate-split`. ⛔ **That file was never committed and
is not on disk** (`git log --all -- splits/` is empty; `splits/` is not gitignored), so `frozen_split`
has taken its redraw branch on every run and the guard has never once been in force. This paragraph
previously stated the file was versioned in the repo, which was untrue when written.

The consequence is bounded and was accepted rather than patched (Selin, 12.08.2026): the assignment on
disk is reproducible from the targets h5ad at seed 42, and it is in any case recoverable from a
committed artifact — `outputs/panel/panel_oof_predictions.csv` names all **153** train+val lines, so the
test set is the labelled lines it omits. Nothing is lost by letting the sweep redraw, and every number
scored on the current split is void on target and panel grounds anyway. **`frozen_split` writes the file
itself when it is missing, so R2 creates it; it is committed then**, which is the point at which the
guard starts protecting anything. The reason is that
eligibility is *not* stable — a line qualifies if it carries at least one CTRP label, so changing the
drug panel or `ctrp_to_h5ad`'s filters silently moves lines between train, val and test, and runs from
either side of the change look comparable when they are not. A line present in the data but missing from
the frozen file raises rather than being assigned: regenerating is a deliberate act that invalidates
comparability with everything run before it.

> ⚠️ **Corrected 12.08.2026.** This paragraph used to end *"**This will hard-fail the panel rebuild**,
> which is the intended behaviour"*. That is no longer true, and the reason is a decision taken the same
> day: **the panel does not feed `ctrp_to_h5ad`**
> ([Step 01](01-datasets-and-harmonization.md#the-panel-does-not-enter-the-target-build-decided-12082026-selin)).
> `M_ctrp` is built over every drug clearing `--min-cell-lines`, so a panel revision does not change
> which lines carry at least one label, does not change eligibility, and therefore cannot trip the
> guard. The coupling this sentence described was real while the panel *could* reach the target matrix;
> severing it is what removed the hard-fail. Changing `--min-cell-lines` or the response source still
> moves eligibility and still trips it — those remain genuine re-freeze events.

**What it invalidates.** Every `X_pca` on disk — and therefore every PCA number in these docs and in
the report — predates this change. ⛔ Nothing is recomputed until the review is finished; see the
banner in [`docs/TODO.md`](../TODO.md). Expect the **~78× input-scale asymmetry** between `X_pca`
and `X_scGPT` (review item 4) to move too, since standardizing genes changes component magnitudes.

**By-product — the HVG ranking scale.** `flavor="seurat"` `expm1`s the matrix back
(`scanpy/preprocessing/_highly_variable_genes.py:373`, scanpy 1.12), bins genes by `log1p(mean)`, and
z-scores log-dispersion within bin. A global constant on the input therefore shifts log-dispersion by
a constant that the within-bin z-score removes exactly — *provided bin membership survives*, which
holds while gene means are ≫ 1 but not once they are pushed below 1. Ranking and projection now read
the same matrix either way, so there is no CPM-vs-CP10K mismatch left; what remains is that the
`/10` divisor adopted above can itself move bin membership. See the warning in the previous section.

### What scGPT is fed, and why its scale does not matter

**scGPT needs none of this.** Its embedding path applies no normalization and no `log1p`
(`scgpt/tasks/cell_emb.py` builds the model directly; `DataCollator(do_binning=True)` bins in the
loader, `scgpt/data_collator.py:87-90`, `n_bins=51`).

**Why the input scale does not matter — an argument from the source, not a measurement
(restated 10.08.2026).** `binning` (`scgpt/preprocess.py:273`) takes the bin edges as quantiles of
**the cell's own non-zero values** (`np.quantile(non_zero_row, np.linspace(0, 1, n_bins - 1))`) and
digitizes that same row against them. Any strictly monotone transform of the row moves the values and
the edges identically, so `np.digitize` returns the same integers — which covers CPM (a per-cell
linear rescale) and `log1p` (a global monotone map) alike. The one thing that would break it is a
transform that collapses distinct values onto equal floats, since that changes which values are tied.

This paragraph replaces an earlier one that claimed the invariance had been *measured*; the claim was
retracted 10.08.2026 and the reasons are in
[Retracted claims](corrections-and-dead-ends.md#the-scgpt-binning-invariance-was-verified-on-200-cells).
The conclusion is unchanged — only its status, from measurement to argument.

**Ties are broken at random, which is why the embeddings are seeded.** `_digitize`
(`scgpt/preprocess.py:239`) resolves values falling on a repeated bin edge with `np.random.rand`, so
binning is **stochastic** wherever a cell has tied expression values — pervasive in single-cell data.
`gen_embeds.py:243-250` seeds `np.random.seed(42)` alongside `torch.manual_seed(42)` for exactly this
reason (and for the 1,199-gene subsampling); without it, embedding the same cell twice would give
different bins.

### The fit is stored, not only the coordinates (10.08.2026)

`obsm["X_pca*"]` records **where each cell landed**. That is the output of the projection and does not
contain the projection: nothing in it says which genes build each component, nor what fraction of the
total variance the 512 kept components account for. Both were computed on the way and discarded, so
until 10.08.2026 the pipeline could answer neither *how much variance does PCA(512) retain* nor *which
genes dominate PC1*, and could not place a new cell in an existing space at all. `add_pca.py::run` now
writes a record per recomputed key into `uns["pca_fits"]` (`add_pca.py::_pca_record`):

| field | what it is |
|---|---|
| `genes` | gene names labelling the axis of `mean`/`std`/`center` and the columns of `components` |
| `mean`, `std` | the per-gene standardization |
| `center` | the mean the PCA subtracts on top of it — not zero, because the clip lands after the divide |
| `scale_max_value` | the ±10 clip (Seurat's `ScaleData(scale.max = 10)` default) |
| `components` | the loadings, `(n_comps, n_genes)` |
| `variance`, `variance_ratio` | per component; for a train-fitted key the ratio is of the **training** cells' variance, so it is not comparable as a percentage with the all-cells key |
| `n_comps`, `seed`, `fitted_on` | provenance of the fit |

The transform it reproduces, for a cell's `log2(1 + CPM/10)` gene vector `x`:

```
z      = clip((x - mean) / std, -scale_max_value, +scale_max_value)
coords = (z - center) @ components.T
```

**Why `uns` and not `varm`.** `varm` is indexed by the file's own `var`, and the targets file keeps
4,576 genes after the scGPT OOV drop while PCA runs on the 5,000 from the convert output — the axes do
not match. Hence the explicit `genes` vector.

**Which key this matters most for.** `X_pca` can be recovered by re-running, since it depends only on
the convert output. `X_pca_train_*` cannot: it is fitted on training cells alone, so re-running
reproduces it only while every input is bit-identical — same cells, same HVG list, same split, same
seed. That stops holding across a re-preprocessing sweep, and the recipe behind the stored coordinates
would be gone. Cost is ~10 MB per fit at `hvg5000` and ~47 MB at `all_genes`, against targets files of
2.2 GB and 8.6 GB.

**Checked while it was written, not by a standing test.** During development the stored record was
applied back to the input matrix on synthetic data to see whether it reproduced the stored coordinates.
It did not, and the reason was a real defect: the record's `mean`/`std` were recomputed by hand with
numpy's `ddof=0` while `sc.pp.scale` uses `ddof=1`. The record now reads scanpy's own
`var["mean"]`/`var["std"]` instead of recomputing them, so the recipe and the transform it describes
come from one source and cannot drift apart that way again. That throwaway check is not in the
repository and no figure from it is quoted here — this is a note on how the code reached its present
form, not a result that can be re-run.

**Consequence — the two fits are now standardized identically (10.08.2026).** That defect exposed a
pre-existing mismatch: `_pca_fitted_on_train` used `ddof=0`, `sc.pp.scale` uses `ddof=1`. The gap is a
factor of √(n/(n−1)), under 0.01 % on this atlas, so no conclusion depends on it — but the two **PCA
fits** were not produced by identical code. `_pca_fitted_on_train` now passes
`ddof=1`, leaving *which cells are seen* as the only difference between the two fits. This changes
`X_pca_train_*` in its last digits and takes effect with the sweep, like everything else here.

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

### Why HVG-5000 is the default (03.08.2026)

Three reasons, in descending order of weight.

**1 — More genes buy nothing measurable.** The gene-set sweep
([Step 05](05-multitask-results.md#gene-set-sweep--heads-beating-vs-gene-count-incl-all_genes-28062026),
`notebooks/analysis/qc/verify_variants.ipynb` §9, 28.06.2026) puts 1k/2k/3k/5k **and**
`all_genes` through the same
5-fold GroupKFold over all 545 drugs. Heads-beating-baseline is **flat across the whole axis** for both
representations, and `all_genes` (PCA 204 ± 86, scGPT 184 ± 90) is no better than `hvg5000`
(210 ± 73 / 189 ± 94) — it sits mid-band for PCA and lowest of all for scGPT. Val MSE is constant at
0.0105–0.0107 throughout.

**2 — At the input length we run, `all_genes` is randomly subsampled and `hvg5000` is not.** We embed
with `max_length=1200` (`gen_embeds.py`), which is `embed_data`'s default and matches the `scGPT_human`
checkpoint's own pretraining configuration — `"max_seq_len": 1200`, `"trunc_by_sample": true`
(`scGPT_human/args.json`). When a cell carries more expressed genes than the cap, the data collator draws
a **random** subset rather than truncating (`scgpt/data_collator.py:143-169`). Measured over all 53,513
cells of every variant's embeddings `.X` (post-OOV, i.e. the matrix actually tokenized) in
`notebooks/analysis/qc/verify_variants.ipynb` §10b–§10c, 05.08.2026; per-cell counts cached
to `outputs/embeddings/scgpt_nonzero_per_cell.npz`:

| Variant | in scGPT vocab | expressed genes / cell (median · mean · max) | cells at or above the cap | **genes scGPT actually saw** (mean) |
|---|---|---|---|---|
| `hvg1000` | 939 | 119 · 122 · 277 | 0 | 122 |
| `hvg2000` | 1,846 | 253 · 257 · 555 | 0 | 257 |
| `hvg3000` | 2,771 | 367 · 373 · 790 | 0 | 373 |
| `hvg5000` | 4,576 | 580 · 589 · 1,208 | **1** (0.002 %) | 589 |
| `all_genes` | 20,570 | 3,461 · 3,550 · 7,797 | **53,513** (100 %) | **1,199** |

The last column is the one that governs how the gene-set axis reads: it is `min(expressed, 1199)`
averaged over cells — what the model received, as against the size the variant is named for. In capped
cells `hvg5000` keeps 99.3 % of a cell's genes; `all_genes` keeps **34.6 %**.

**The axis is real for scGPT, but compressed at the top.** Across the five variants scGPT's actual input
grows 122 → 257 → 373 → 589 → 1,199, so more genes *do* reach it at every step — but the nominal
22,722-gene condition delivers only about **twice** what `hvg5000` does, not twenty times, and delivers
it as a random draw rather than a dispersion-selected set.

> ⚠️ **Supersedes the 03.08.2026 figures.** This table previously carried `hvg5000` and `all_genes` only,
> from an ad-hoc `h5py` pass that existed as a shell command rather than as code. The re-runnable
> measurement reproduces both rows exactly and adds the three smaller variants.

So at HVG-5000 scGPT sees every gene surviving the filter in all but one cell, whereas at `all_genes` it
sees a random ~34 % of each cell — a different third on each run. Under this configuration, feeding the
full transcriptome does not deliver more information to scGPT; it only randomises which fraction reaches
the model.

**Confirmed against the publication** — Cui, H., Wang, C., Maan, H. *et al.* "scGPT: toward building a
foundation model for single-cell multi-omics using generative AI." *Nature Methods* **21**, 1470–1480
(2024), doi:10.1038/s41592-024-02201-0, Methods, *Implementation details*:

> "Note that, in pretraining, only genes with non-zero expression are input to the model. We set a
> maximum input length of 1,200. For cells with a number of non-zero genes larger than the maximum input
> length, 1,200 input genes would be randomly sampled at each iteration."

So the random subsampling our pipeline hits at `all_genes` is exactly the operation used in pretraining,
not an artifact of the inference wrapper.

> **1,200 is a configuration choice, not a limit of the model.** The paper defines M as "a predefined
> maximum input length" and states that "the input dimension M **can reach tens of thousands of genes**,
> substantially exceeding the input length of conventional transformers commonly used in NLG", handled
> via FlashAttention (Methods, *Input embeddings*). The code agrees: `TransformerModel` builds gene
> inputs through `GeneEncoder` alone and never instantiates `PositionalEncoding`
> (`scgpt/model/model.py:86` vs `:742`), so genes are an unordered set and nothing architecturally caps
> the sequence, and the authors' own tutorials run the same pretrained weights at `max_seq_len` 1536,
> 3001 and 4001 (`max_seq_len = n_hvg + 1`). The argument above therefore justifies HVG-5000 **given**
> our 1,200 setting, and would have to be redone if that setting changed. Note that the FlashAttention
> the paper relies on to reach those lengths is **not available to us** — see
> [Compute environment and its limits](#compute-environment-and-its-limits-03082026).

**The paper's HVG counts bracket ours.** Depending on task the authors select **1,200** HVGs (multiomic
integration), **3,000** (cell annotation) and **5,000** (perturbation prediction; Methods), using Scanpy
for "normalization, log transformation and HVG selection". Our 5,000 sits at the top of that range, not
outside it.

---

## Compute environment and its limits (03.08.2026)

Everything scGPT-related runs on one Apple Silicon machine (arm64, macOS 26.6, `torch` 2.3.1 in the
separate scGPT venv). Two hardware limits shape what the embedding step can and cannot do, and both were
verified by running them, not by reading docs.

### FlashAttention is unavailable, and cannot be installed

The paper leans on FlashAttention to reach long inputs. We cannot use it. `pip install flash-attn` fails
at **metadata generation**, before any compilation:

```
OSError: CUDA_HOME environment variable is not set.
  torch/utils/cpp_extension.py -> CUDAExtension -> library_paths(cuda=True)
```

flash-attn ships CUDA kernels and requires an NVIDIA GPU; there is no Apple Silicon build, and
`torch.cuda.is_available()` is `False` here. **This is not a configuration problem and has no workaround
on this machine.**

The consequence is easy to miss: `embed_data` requests `use_fast_transformer=True` by default, and when
flash-attn is absent the model **silently downgrades** to the standard PyTorch transformer with a
`UserWarning`, not an error (`scgpt/model/model.py:75-82`). So attention is the ordinary quadratic
implementation, and any increase in `max_length` costs O(n²) with none of the mitigation the paper
assumes.

### MPS works — 4× faster, same embeddings — but needs a fallback flag

`gen_embeds.py` **used to** hardcode `device = "cpu"` with the comment "MPS disabled" (changed
03.08.2026). The reason was a single missing operator:

```
NotImplementedError: The operator 'aten::_nested_tensor_from_mask_left_aligned'
is not currently implemented for the MPS device.
```

Setting **`PYTORCH_ENABLE_MPS_FALLBACK=1`** runs that one op on CPU and the rest natively on MPS.
Benchmarked on 256 real cells (`hvg5000`, `max_length=1200`, `batch_size=64`, seeded identically for both
devices; scratchpad `mps_smoke2.py`, 03.08.2026):

| Device | 256 cells | extrapolated to 53,513 cells |
|---|---|---|
| CPU | 26.9 s | ~94 min |
| **MPS** (with fallback) | **6.6 s** | **~23 min** |

**The two devices agree numerically:** max absolute difference **2.7 × 10⁻⁷**, cosine similarity
**1.000000** (min and mean) across all 256 cells. Cell embeddings are L2-normalised by scGPT
(`cell_emb.py`), so that residual is float32 noise, not a difference in result.

> ⚠️ The full-run figures are **extrapolated linearly from 256 cells** — no complete run has been timed.

**`gen_embeds.py` now selects MPS and seeds itself** (03.08.2026): it sets
`PYTORCH_ENABLE_MPS_FALLBACK=1` *before* the scgpt import — PyTorch reads that variable when the MPS
backend registers its dispatch keys, i.e. at `torch` import time, so a later assignment is silently
ignored — then picks `mps` when `torch.backends.mps.is_available()` and falls back to `cpu` otherwise,
and seeds `torch` and `numpy` with **42** before embedding. That covers both random draws in the path:
the collator's gene subsampling (`data_collator.py:169`) and the tie-breaking inside value binning
(`preprocess.py:_digitize`); the DataLoader uses `num_workers=0`, so there is no per-worker reseeding to
handle.

Verified by running `gen_embeds.py` **twice** over a 256-cell `all_genes` subset — the case where every
cell is subsampled, so the seed is actually exercised — and comparing `obsm["X_scGPT"]`: **bitwise
identical, max absolute difference 0.0**.

> ⚠️ **The embeddings currently on disk predate this change.** Every `..._scGPT_human_embeddings.h5ad`
> under `processed/scRNAseq_SCP542/` was generated on CPU without a seed and has **not** been
> regenerated. Until it is, the files do not match what the script now produces.

### What this rules out

Raising `max_length` to cover every gene in `all_genes` would need ~7,798 (the largest expressed-gene
count in any cell). Without FlashAttention that is roughly a 40× increase in attention cost over 1,200
and forces a smaller batch size for memory, on top of putting the input far outside the 1,200-token
regime the checkpoint was pretrained in. **Full-coverage single-pass embedding is therefore expensive and
untested here.**

### Decision — one seeded draw at 1,200; `all_genes` is a sanity check (03.08.2026)

Four routes were considered for making the `all_genes` scGPT embedding something other than one
unreproducible random draw: (**A**) keep `max_length=1200` and seed it; (**B**) keep 1,200 and average
*k* independently seeded draws per cell; (**C**) raise `max_length` to ~7,798 for full single-pass
coverage; (**D**) drop `all_genes` as a scGPT comparator entirely.

**Chosen: A — a single draw at `max_length=1200`, seeded with 42.** Decision by Selin, 03.08.2026.
`all_genes` is hereby **a sanity check, not a full-transcriptome comparator**; `hvg5000` is the primary
result.

Why A over B, the closest alternative:

1. **Interpretability.** Every A embedding is an actual model output produced at the input length the
   checkpoint was pretrained on. B's averaged vector is a **constructed object the model never emitted**:
   scGPT L2-normalises each cell embedding (`cell_emb.py`), so the mean of *k* unit vectors falls inside
   the sphere and has to be renormalised by hand. For a representation that is fed to a downstream
   predictor and then interpreted, "this is what scGPT returned" is a materially easier thing to defend
   than "this is the renormalised mean of five draws".
2. **Cost.** B multiplies the embedding step by *k* and obliges a re-run of anything computed on the
   embeddings. Measured: ~23 min per pass on MPS (extrapolated from 256 cells, above), so *k* = 5 is
   roughly 2 h of embedding plus downstream re-runs — **not prohibitive, and this is the weaker of the
   two reasons.** Reason 1 is the load-bearing one.
3. **It keeps the two variants comparable.** Only **1 cell of 53,513** is subsampled in `hvg5000`, so B
   would build `all_genes` by a procedure `hvg5000` never uses — trading a known sampling artifact for a
   new asymmetry between the very variants being compared. For a sanity check, that defeats the purpose.

> ⚠️ **The consequence, stated plainly.** Under A the `all_genes` scGPT embeddings still represent a
> **random ~34 % of each cell's expressed genes**. They become *reproducible*, not *complete*. No claim
> of the form "scGPT does not benefit from the full transcriptome" can rest on them, because scGPT never
> received the full transcriptome. What `all_genes` can support is the narrower statement that
> HVG-5000 loses nothing detectable relative to a same-length draw from the unfiltered gene set.
> Route **C** remains the only one that would license the stronger claim, and it is not being taken.

**3 — Storage cost (but *not* embedding time).** `all_genes` occupies **26 GB** on disk against **11 GB**
for `hvg5000` (`du -sh data/processed/scRNAseq_SCP542/*/`, 03.08.2026).

> ⚠️ **Correction, 03.08.2026.** An earlier version of this section claimed re-embedding time "scales
> similarly" with the gene set. **It does not.** Because `max_length` caps the sequence at 1,200, scGPT
> reads at most 1,199 genes per cell in *either* variant, so the embedding step costs essentially the
> same for `all_genes` as for `hvg5000`. The variant cost difference is disk plus the convert/PCA/IO
> steps, not the transformer. Reason 3 is therefore a **storage** argument only, and is the weakest of
> the three. See [Compute environment and its limits](#compute-environment-and-its-limits-03082026) for
> the measured numbers.

> ⛔ **The sweep has no live numbers.** Those above were trained on `mean_pv`, removed with its reader
> code on 11.08.2026, so they cannot be regenerated; §9 now pins `SCORE = 'auc_cc'` and must be re-run.
> The `all_genes` column additionally predates the seeding fix and the gene-symbol repair, so it was
> never exactly reproducible ([below](#decision--one-seeded-draw-at-1200-all_genes-is-a-sanity-check-03082026)).
> **The conclusion is expected to carry over but is not currently evidenced.**

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
`notebooks/analysis/qc/verify_variants.ipynb` checks the gene counts directly and plots the two variants'
UMAPs side by side. Evaluation of the all-genes side is still pending.

✅ On-plan / closes part of the HVG deviation by enabling the full-transcriptome comparison.

---

## Latent-space validation (UMAP, Fig. 3 / Fig. 4)

`notebooks/analysis/qc/verify_variants.ipynb` (§7) is the **standalone validation** (not part of the
orchestrator): it builds PCA-vs-scGPT UMAPs for **both** variants via `sc.pp.neighbors` + UMAP,
colored by `Cancer_type` (**Fig. 3**) and `viability_paclitaxel` (**Fig. 4**). It visually confirmed
the hypothesis: PCA = discrete tissue "islands", scGPT = continuous shared manifold; paclitaxel
sensitivity mixed across the scGPT manifold.

**Outputs** (all under `notebooks/outputs/embeddings/`, written by `analysis/qc/verify_variants.ipynb` §8):

| File | What it shows |
|---|---|
| `umap_cancertype_pca_vs_scgpt.png` (dpi 300) | the 2-panel PCA-vs-scGPT comparison by cancer type — the headline latent-validation figure |
| `umap_sweep_cancertype.png` (dpi 200) | the same contrast across the full gene-set sweep grid; **the tissue-islands vs continuous-manifold split holds at every gene count**, so it is not an artifact of the HVG choice |
| `variants.png` | QC check that the `hvg5000` and `all_genes` outputs agree |

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

**Raw inputs** (shared, not per-variant):

- scRNA-seq counts → `data/scRNAseq_SCP542/expression/CPM_data.txt`
- cell metadata → `data/scRNAseq_SCP542/metadata/Metadata.txt`
- CTRPv2 tables → `data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.*`

**There is no separate file per representation, per drug, or per task.** One trainable file per variant
bundles everything, and the representation / drug / task are *selected at training time*. Exact location
of each artifact, under `processed/scRNAseq_SCP542/` (`<score>` = the `--score` suffix, `auc` by default):

| Artifact | HVG filtered? | File | Stored as | Shape / genes |
|---|---|---|---|---|
| Counts (CPM) | **filtered** | `hvg5000/SCP542_CCLE.h5ad` | `.X` | 53,513 × 5,000 |
| Counts (CPM) | **non-filtered** | `all_genes/SCP542_CCLE.h5ad` | `.X` | 53,513 × 22,722 |
| scGPT embeddings | **filtered** | `hvg5000/SCP542_CCLE_scGPT_human_embeddings.h5ad` | `obsm["X_scGPT"]` | 53,513 × 512 (from 4,576 in-vocab genes) |
| scGPT embeddings | **non-filtered** | `all_genes/SCP542_CCLE_scGPT_human_embeddings.h5ad` | `obsm["X_scGPT"]` | 53,513 × 512 (from 20,570 in-vocab genes) |
| PCA | **filtered** | `hvg5000/…_with_targets_<score>.h5ad` | `obsm["X_pca"]` | 53,513 × **512** (computed on the 5,000 HVG) |
| PCA | **non-filtered** | `all_genes/…_with_targets_<score>.h5ad` | `obsm["X_pca"]` | 53,513 × **512** (computed on all 22,722) |
| Drug labels — all 545 | both | `<variant>/…_with_targets_<score>.h5ad` | `obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`, `uns["ctrp_score"]` | 53,513 × 545 |
| Drug labels — one drug | both | same targets file | one column of `Y_ctrp` via `--drugs paclitaxel`; legacy `obs["viability_paclitaxel"]` | 53,513 × 1 |
| Split — shared, cell-line-grouped | both | same targets file | `obs["split_ctrp"]` | per-cell |
| Split — paclitaxel-only (legacy) | both | same targets file | `obs["split_paclitaxel"]` | per-cell |

**How a run selects from these** — no new input files are written per run:

| Choice | Flag | Selects |
|---|---|---|
| gene set | `--variant {hvg5000, all_genes}` | which folder |
| target score | `--score {auc, auc_z, mean_pv}` | which targets file in that folder |
| representation | `--use-rep {X_scGPT, X_pca}` | which `obsm` key |
| task | `--drugs paclitaxel` vs omitted | K = 1 vs K = 545 |

**Training outputs:** each run writes `runs/<timestamp>_<tag>/` (gitignored) holding `best_model.pt`,
`config.json`, `run_meta.json` (records the variant via the targets path), `history.csv`,
`summary.json`, `per_drug_results.csv`, plus one index row in `runs/runs_index.csv`. Column-level detail
in [Step 05](05-multitask-results.md#run-versioning-26052026).

Per variant, the two expensive files are built **once** and shared by every score:
`SCP542_CCLE.h5ad` → `..._scGPT_human_embeddings.h5ad`. The targets step then forks per score into
the trainable file (`X_scGPT`, `X_pca`, `Y_ctrp`, `M_ctrp`, `split_ctrp`, `split_paclitaxel`):

- `..._with_targets_auc_cc.h5ad` — **the default** (`--score auc_cc`), CurveCurator AUC
- `..._with_targets_ln_ic50_cc.h5ad` — `--score ln_ic50_cc`, the same fit's IC50 (incomplete by
  construction; see [Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026))

⛔ Files named `..._with_targets_auc.h5ad`, `..._auc_z.h5ad` and `..._with_targets.h5ad` may still sit
on disk from before 11.08.2026. They hold **different quantities under similar names** and no code
produces them any more — `auc`, `auc_z` and `mean_pv` were removed with their readers
([Corrections](corrections-and-dead-ends.md#the-auc-target-was-divided-by-the-wrong-quantity)). They
are the reason the new measures carry the `_cc` suffix instead of reusing `auc`.

### Reproduce

The runnable walk-through is the numbered notebooks themselves: `notebooks/1_data.ipynb` (`fetch` →
`convert`) and `notebooks/3_representations.ipynb` (`scgpt` → `targets` → `splits` → `pca`). Both drive
`scripts/preprocessing/pipeline.py`, one function per step, each owning its own guard and preconditions.

⚠️ **`run_preprocessing.py` was archived on 12.08.2026** and there is no CLI replacement — the step order
is the notebook numbering, and a second copy of it in a CLI was a second thing to keep in step
([`scripts/archive/README.md`](../../scripts/archive/README.md)). The commands this section used to give
therefore no longer run. To rebuild without a browser:

```bash
# From scratch. Set SCGPT_PYTHON in the stage-3 bootstrap cell first: the scgpt
# step needs the separate scGPT venv and raises without it (no interactive fallback).
uv run jupyter nbconvert --execute --inplace notebooks/1_data.ipynb
uv run jupyter nbconvert --execute --inplace notebooks/3_representations.ipynb
```

Anything more selective is a direct call, which is what one-function-per-step buys — there is no
`--start-at` to reconstruct:

```python
from scripts.layout import PipelinePaths
from scripts.preprocessing import pipeline

# Add the second measure on top of the existing convert + embeddings. The embedding
# is score-independent, so scgpt is simply not called; nothing is skipped or reused
# by flag, it is just not invoked.
paths = PipelinePaths.build(None, "hvg5000", "ln_ic50_cc")
pipeline.targets(paths)
pipeline.splits(paths)
pipeline.pca(paths)

# Recompute only the 512-d PCA baseline in place (what the 512-d switch needed).
pipeline.pca(PipelinePaths.build(None, "hvg5000", "auc_cc"), n_comps=512, force=True)
```

`--score` defaults to `auc_cc`; each measure writes its own targets h5ad, so the two never overwrite
each other. Training still has a CLI:

```bash
uv run scripts/training/train_multitask.py --use-rep X_scGPT --score auc_cc   # every kept drug
uv run scripts/training/train_multitask.py --use-rep X_pca --score ln_ic50_cc --drugs paclitaxel
```
