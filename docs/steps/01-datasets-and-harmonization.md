# Step 01 — Datasets & harmonization

*Part of [OncoTox project progress](../project_progress.md). Covers: the raw datasets, what each
response assay actually measures, the cross-dataset overlap/coverage audit, and how cell lines and
compounds are harmonized across sources — distinguishing the one-off exploratory audit from the
normalization that actually feeds training.*

Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition**.

---

## Data collection (26.03–30.03.2026)

| Dataset | Role | Key numbers | Used? |
|---|---|---|---|
| **SCP542** scRNA-seq (Kinker et al. 2020; used in PERCEPTION) | single-cell input | **53,513 cells × 22,722 genes**, **198 unique cell lines** | ✅ primary |
| **CTRPv2** (Cancer Therapeutics Response Portal v2) | dose-response labels | 1,107 cell lines, **545 compounds**; training target `auc_cc`, with `ln_ic50_cc` as the alternative — both from DrEval's CurveCurator re-fit ([below](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)) | ✅ primary |
| **PRISM** Repurposing (Public 24Q2) | multiplexed viability (LFC) | 915 cell lines, 6,575 compounds | downloaded, not used |
| **GDSC2** | `LN_IC50` / AUC | 967 cell lines, 295 drugs | downloaded, not used |

**What the assays measure (and why the resolution mismatch matters).** SCP542 (from the Broad
Single Cell Portal; `.X` stored as **CPM** — counts-per-million, library-size-normalized) is a
pan-cancer cell-line scRNA-seq atlas capturing **single-cell** heterogeneity. The three response
datasets are all **bulk, cell-line-level** drug screens:

- **CTRPv2** screens each (cell line × compound) over a **concentration series** (usually 16 points)
  with a CellTiter-Glo-type readout, and ships it at two levels: the raw per-concentration percent
  viability `cpd_avg_pv` (fraction surviving vs vehicle/DMSO controls) in
  `v20.data.per_cpd_post_qc.txt`, and the **post-QC sigmoid fit** of that series —
  `area_under_curve`, `apparent_ec50_umol`, slope, and per-parameter confidence intervals — in
  `v20.data.curves_post_qc.txt`. Either is joined to names via `v20.meta.per_experiment.txt` →
  `v20.meta.per_cell_line.txt` (`ccl_name`) and `v20.meta.per_compound.txt` (`cpd_name`). It is an
  **efficacy** metric. ⚠️ Since 11.08.2026 the label comes from **neither** of these directly: the
  target is CurveCurator's re-fit of the same raw measurements
  ([below](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)). What that choice forces on the
  model is [Step 03](03-model-and-training-design.md).
- **GDSC2** reports `LN_IC50` (natural-log half-maximal inhibitory concentration) and AUC — the same
  curve-fit family as CTRP's `area_under_curve`, which is what makes
  [Step 06](06-planned-work.md#a-cross-database-integration) tractable.
- **PRISM** is a barcoded, multiplexed viability assay reporting log fold-change vs control — very
  large and very sparse.

The whole project exists to **bridge this bulk-to-single-cell gap** (plan §Understanding): the
bulk labels are the only high-volume continuous signal available, so they are mapped onto single
cells as weak supervision ([Step 03](03-model-and-training-design.md)).

✅ On-plan: SCP542 + CTRPv2 are the designated primary pair; PRISM/GDSC reserved for later.

### Provenance — what was retrieved, from where, when

*Moved here 30.07.2026 from the dated log, which is where it had been decaying. This is the record
[TODO](../TODO.md) review item 1 and FAIRER **A**/**E** ask for. Release identifiers, filenames and
roots re-checked against the filesystem 10.08.2026 (audit item 01).*

**The sources sit in two different roots.** Nothing else in the project said so, and paths quoted
elsewhere in this file are relative to whichever root the source landed in:

| Root | Sources it holds |
|---|---|
| `~/Desktop/OncoTox/data/` — the pipeline's `--data-root`, hard-coded as `DEFAULT_DATA_ROOT` in `scripts/layout.py` | SCP542 (`scRNAseq_SCP542/`), CTRPv2 (`metadata/CTRPv2.0_2015_ctd2_ExpandedDataset/`), PRISM (`metadata/PRISM_REPURPOSED/`), Repurposing Hub annotations (`metadata/repo-drug-annotation-20200324.txt`), DrugBank (`full database.xml`) |
| `<repo>/data/` — untracked, see [Licences](#licences-and-terms-of-use-of-the-source-data-checked-28072026) | GDSC2, the compound catalog the [overlap audit](#overlap--coverage-audit-03042026) writes (`drug/`), and the abandoned scDrugAtlas download (`scDrugAtlas/`, identified below) |

Only the first root is reachable from code; the second is referenced by the audit notebook alone.

**Two downloaded files belong to no pipeline and had no record until 10.08.2026.** Both are identified
here so neither is mistaken later for something the analysis depends on:

- `scRNAseq_SCP542/other/CCLE_scRNAseq_github/` (~3.9 GB, files dated 29.08.2020) is the **authors' own
  reproduction bundle**, from `https://github.com/gabrielakinker/CCLE_heterogeneity` — the repository
  named in the Kinker et al. Code availability statement (p. 14). It holds `CCLE_heterogeneity_Rfiles/`
  (their CPM matrix and metadata as `.RDS`, copy number by gene, gene loci, the tumour comparison
  matrix, the literature metaprogram table) and `Expected_results/` (3.2 GB of module 1–6 reference
  outputs). Nothing in this project reads it. It is the reference implementation of *their*
  heterogeneity analysis, which is the natural comparison point for research question 2, not for the
  response prediction.
- `data/scDrugAtlas/ea472fa4aec64cabaef5194af0ba5ba0.h5ad` (45 MB) is **not SCP542**: 1,761 cells ×
  18,919 genes, `obs` = `drug` (control 1,061 / `Pal` 693 / `palbociclib` 7), binary `response`
  (1: 1,064, 0: 697), two batches, genes as Ensembl IDs, already carrying `X_pca` and a neighbour
  graph. It is a palbociclib treatment-versus-control perturbation experiment with sensitivity labels,
  evidently two sources stitched together — batch 2 is seven cells under a different spelling of the
  drug. ⚠️ **Which study it came from cannot be recovered:** the hash filename carries no provenance and
  nothing accompanies the file. It arrived with the scDrugAtlas line of work, which was abandoned
  ([Corrections](corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources)), and no result
  in this project rests on it.

| Source | Retrieved | From |
|---|---|---|
| **SCP542** scRNA-seq (53,513 cells × 22,722 genes) | 30.03.2026 | Broad Single Cell Portal, study **SCP542** *"pan-cancer cell line heterogeneity"* — <https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity>. Source publication Kinker et al., *Nat Genet* 2020, <https://doi.org/10.1038/s41588-020-00726-6>. This is the dataset the PERCEPTION paper used. **Also deposited at GEO under accession `GSE157220`** (from the paper's Data availability statement, p. 13) — the only persistent accession among our four primary sources, and the durable route if the portal study is ever revised or withdrawn. ⚠️ The portal exposes **no version or revision identifier**, so the retrieval date is the only version reference this source has; recorded as sufficient (Selin, 10.08.2026). Files taken: `expression/CPM_data.txt` (5.4 GB, the matrix everything downstream uses), `metadata/Metadata.txt`, `other/UMIcount_data.txt` (3.5 GB raw UMI counts — **downloaded but never used**; see [Step 02](02-preprocessing-and-embeddings.md#the-expression-transform-is-the-datasets-own-05082026)), plus the per-cancer-type tSNE cluster files. `file_supplemental_info.tsv` lists the study's full file set. |
| **CTRPv2** viability + compound annotations — **no longer the target source** (see the row below); retained for compound metadata (`broad_cpd_id`) and read by several analysis notebooks | 30.03.2026 | DepMap portal, all-data tab — <https://depmap.org/portal/data_page/?tab=allData>. Release `CTRPv2.0_2015_ctd2_ExpandedDataset`. Cell lines join via `v20.meta.per_cell_line` (`master_ccl_id`) and experiments via `v20.meta.per_experiment`. The 15 files carry their upstream 2015 dates, and the release ships **`MANIFEST.txt` with an MD5 per file** — the only integrity anchor any source provided. ✅ **All 15 verified against it 10.08.2026: every checksum matches**, so this copy is byte-identical to the release the Broad published. The other three sources shipped no checksums and no comparable check is possible for them. |
| **CTRPv2 responses, reprocessed by DrEval** — ⭐ **the training target since 11.08.2026** | 11.08.2026 | Zenodo record **`21807175`** (*"Dataset for drevalpy"*, published 2026-08-05), DOI <https://doi.org/10.5281/zenodo.21807175>. Files `CTRPv2.zip` and `meta.zip`, fetched and **MD5-verified against the record's published checksums** by `scripts/sources/fetch_ctrp_response.py`, which writes `provenance.json` beside the data. Not obtained via `drevalpy.datasets.loader.load_ctrpv2()`: that resolves the *concept* DOI to whatever release is current and re-downloads unconditionally, so the version it returns changes silently — a target whose version cannot be named is not a citable source. The record is pinned in `layout.ZENODO_RESPONSE_RECORD`; bumping it is a target change. **What it contains:** CTRPv2's own raw dose-response measurements, normalised per replicate against the no-drug control and re-fitted with CurveCurator — see [what changed and why](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026). `meta.zip` also carries **Cellosaurus release 52.0 (10 April 2025, CC BY 4.0)**, pinned by this record because Cellosaurus publishes no stable per-release download URL of its own. |
| **PRISM** Repurposing, Public 24Q2 | 30.03.2026 | DepMap — `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`. Taken specifically because it covers **failed** drugs as well as approved ones. The full 24Q2 file set was downloaded (LFC, LMFI, QC, metadata); only the extended primary matrix is read. |
| **GDSC2** IC50 + raw data | 26.03.2026 | Sanger Cell Model Passports — <https://cellmodelpassports.sanger.ac.uk/downloads>. Release **27 Oct 2023**: `GDSC2_fitted_dose_response_27Oct23.xlsx` ("GDSC2 IC50 Data") and `GDSC2_public_raw_data_27Oct23/` ("GDSC2 Raw Data", 2.1 GB + description PDF). |
| **Clinical-phase annotations** | 30.03.2026 | Broad Drug Repurposing Hub — <https://repo-hub.broadinstitute.org/repurposing#download-data>. File `repo-drug-annotation-20200324.txt`, i.e. the **24.03.2020** annotation release. |
| **DrugBank** full database (XML) | 06.03.2026 | <https://go.drugbank.com/>, academic download; account required. Licence terms above. ⚠️ The date is the file's mtime, not a recorded download — no acquisition record exists for this source, and the XML carries no version in its filename. |
| **FDA-approved anticancer drug list** — ⭐ the criterion the [drug panel](#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026) is selected on | 11.08.2026 | Sun, J. *et al.*, *BMC Systems Biology* **11**(Suppl 5), 87 (2017), [doi:10.1186/s12918-017-0464-7](https://doi.org/10.1186/s12918-017-0464-7), **Table 1** — 150 drugs approved 1949–2014, 61 cytotoxic and 89 targeted. Retrieved as JATS XML from NCBI E-utilities (`efetch`, `db=pmc`, `PMC5629554`) by `scripts/sources/fetch_sun2017_drugs.py` → `reference/sun2017_fda_anticancer_drugs.csv` (committed; CC BY 4.0). **Table 1 is the whole dataset** — the paper ships no supplementary file. ⚠️ PMC publishes no checksum for a rendered article, so the parse is verified against the counts the paper's own Results state (150 / 61 / 89) and **raises** if they do not reproduce. ⚠️ **The list stops at 2014**: nothing approved since is in it — no CDK4/6 inhibitor, no third-generation EGFR inhibitor, and no checkpoint inhibitor beyond pembrolizumab and nivolumab. Harmless against a screen run in 2012–2015, but it is not a current statement of what is approved. |
| **PubChem** compound identifiers and parent-compound relations | 11.08.2026 | PubChem PUG-REST (<https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>) — `compound/name/{name}/cids` and `compound/cid/{cid}/cids?cids_type=parent`, via `scripts/sources/pubchem.py`. Needed because Sun names drugs by INN and salt form while CTRPv2 uses development codes and free bases, so no string rule relates them. **Cached and committed** as `reference/pubchem_parent_cids.csv` (566 CIDs), so the panel reproduces on a machine with no network and cannot drift when PubChem is re-curated. Citation: Kim, S. *et al.* PubChem 2023 update. *Nucleic Acids Research* **51**, D1373–D1380 (2023), [doi:10.1093/nar/gkac956](https://doi.org/10.1093/nar/gkac956). |

⚠️ **How GDSC2's `LN_IC50` was processed is not documented.** The Sanger documentation link
(`depmap.sanger.ac.uk/documentation/gdsc/`) was dead, and the DepMap/GDSC team was asked directly with no
response ([Corrections](corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources)). This is an
open gap for any future GDSC head.

### Upstream QC — what SCP542 already had done to it (05.08.2026)

*The data arrives quality-controlled. Dying cells, low-quality cells and suspected doublets **have
been removed** — by Kinker et al., before publication — so the pipeline applies no second filtering
pass, and needs none. Recording their decisions is the point: the `53,513 × 22,722` figures quoted
throughout this project are the output of someone else's quality control, not a raw measurement, and a
choice we cannot defend is a choice we should not be making silently. Source: Kinker et al.,
**Nature Genetics 52, 1208–1218 (2020)**,
doi:10.1038/s41588-020-00726-6 — Methods, "Processing of scRNA-seq data", and Results, "Pan-cancer
scRNA-seq of human cell lines". `references.bib` key `scp542`.*

**Design.** 198 CCLE lines profiled in nine multiplexed pools — eight CCLE pools of 24–27 lines each,
grouped by doubling time, plus one custom head-and-neck pool — on 10x Genomics Chromium 3′ v2, for an
average of **~280 cells per line**.

| Stage | What they did |
|---|---|
| Barcode filtering, alignment, UMI counting | **Cell Ranger 3.0.1** (10x Genomics) |
| Assigning cells to cell lines | Consensus of two approaches, genetic (SNP) and expression-based. Inconsistent assignments — mostly cells with low SNP coverage — were **excluded** |
| Cell quality | Detected genes used as the quality proxy; cells **"conservatively" retained at 2,000–9,000 detected genes**. Low-quality cells and suspected doublets excluded |
| Cell-line quality | Lines with **fewer than 50 assigned cells** excluded |
| Result | **53,513 cells, 198 cell lines, 22 cancer types**; mean **19,264 UMIs** and **3,802 genes** detected per cell |

**What the 2,000–9,000 window actually is.** For each cell, count how many of the 22,722 genes have a
non-zero value — its *detected genes*, averaging 3,802 across this dataset. Cells were kept only when
that count fell inside the window. The two bounds catch different failures:

- **Floor (2,000)** — too few genes detected means almost no RNA was captured: an empty droplet
  containing only ambient RNA, or a dying cell whose RNA has leaked out. The profile is mostly noise.
- **Ceiling (9,000)** — too many usually means the droplet held **two** cells rather than one. Two
  transcriptomes together show more distinct genes than either alone, so an unusually high count is the
  signature of a doublet.

Together: *keep droplets that look like exactly one intact cell.* Every cell in our 53,513 therefore
has between 2,000 and 9,000 detected genes; nothing outside that range exists anywhere in our data.

⚠️ **What the measure cannot distinguish.** A broken cell and a real cell that genuinely contains
little RNA look identical by this test — a quiescent or senescence-like cell is small and
transcriptionally quiet, so it has few detected genes and falls near the floor. Both are cut. That is
standard and usually harmless; it matters here because those quiet states are plausibly among the rare
survivors that **research question 2** is about.

**Quantification.** `CPM[i,j] = 10⁶ × UMI[i,j] / Σ UMI[·,j]`, then `E[i,j] = log2(1 + CPM[i,j]/10)`.
The portal distributes the **CPM** matrix (`CPM_data.txt`), not `E` — verified from the values, which
run to tens and hundreds. How and why we apply their `E` transform ourselves is
[Step 02](02-preprocessing-and-embeddings.md#the-expression-transform-is-the-datasets-own-05082026).

**Why we do not re-filter.** Their doublet call rests on a signal the distributed matrix does not
carry: cells were assigned to lines **by genotype**, so a droplet holding two different cell lines
shows mixed SNPs and is caught directly. That is a far stronger discriminator than the
expression-based inference a tool such as Scrublet performs, and it cannot be reconstructed from
expression alone. Running scanpy's standard QC on top would filter an already-filtered population a
second time, cutting real cells at the margins for no gain. **Cell-level QC is therefore inherited and
cited, not repeated** — which is ordinary practice for a published, quality-controlled dataset.

**Their gene filters are *not* in the file we hold.** For their own heterogeneity analyses they
filtered genes two ways — per cell line, `E > 3.5` in at least 2% of cells (≈6,758 genes per line);
across lines, the 7,000 most highly expressed (minimum average expression 12 CPM). Our matrix has
**22,722** genes, so neither was applied to the distributed data. Those filters sit downstream of what
we received, and our HVG selection is independent of them.

**Their centering (`Er`, `ER`) — one is already covered, the other would break the task.**

- **`ER`**, global centering (`E − mean(E)` over all 53,513 cells), is a subset of what we already do:
  `sc.pp.scale` in [Step 02](02-preprocessing-and-embeddings.md) centers *and* scales each gene to unit
  variance. Adopting `ER` explicitly would mean centering without scaling — strictly less. No action.
- **`Er`**, per-cell-line centering (`E − mean(E)` over that line's cells), removes between-line
  variation and leaves only within-line variation. **The response label is constant within a cell
  line**, so all of its signal is between-line; centering per line would delete exactly the variance
  the target depends on and the task would collapse. Not adopted, and not a close call.
  It is, however, the natural representation for **research question 2** — Kinker et al. use `Er`
  precisely because they study within-line heterogeneity. Worth revisiting when the objective stops
  being "predict a line-constant scalar from each cell", i.e. at the MIL / attention-pooling step
  ([TODO](../TODO.md) S2).

**Two things this pins down.**

1. **The detected-gene window is an adopted analysis decision, not an inherent property of the data.**
   It is the only cell-level QC the project stands on, it was made by someone else, and it bounds every
   heterogeneity claim we can make (see the window and its blind spot above).
2. **198 is where the cell-line funnel starts.** The 190 roster name matches and 180 screened lines in
   the [overlap audit](#overlap--coverage-audit-03042026) below are counted down from this 198, not
   from CCLE at large.

✅ **Closed 10.08.2026 (FAIRER E): SCP542 carries no study-level licence beyond the portal terms of
use.** Checked in the source publication itself — Kinker et al., *Nature Genetics* 52, 1208–1218 (2020),
**Data availability** statement, p. 13: *"Raw and processed scRNA-seq data are available through the
Broad Institute's single-cell portal (SCP542) and at the Gene Expression Omnibus (GEO) (accession number
GSE157220)."* No licence is named and no restriction is declared; the Reporting Summary (p. 28) repeats
the statement under a heading that expressly asks for "a description of any restrictions on data
availability" and lists none. The portal ToS is therefore the operative term, as recorded under
[Licences](#licences-and-terms-of-use-of-the-source-data-checked-28072026).

### Design decisions taken with the advisor (27.03–03.04.2026)

The project's framing was settled by asking, not assumed. The three questions put to Artem on 27.03 were:
is the target **cytotoxicity/efficacy** (does the drug kill heterogeneous cancer cells) or **adverse
patient toxicity** (dangerous side effects in healthy tissue); should the model predict **continuous** or
**binary** values; and which dataset should be primary. The answers, which several current choices rest
on:

- **Define toxicity by the labels the chosen datasets actually carry** — do not import an external
  definition. This is why the project predicts efficacy-type response and not patient-level toxicity.
- **If several response types exist and the effort is reasonable, model them all** — the origin of the
  multi-metric goal in [Step 06 · A](06-planned-work.md#a-cross-database-integration).
- **Prefer continuous outputs** (IC50 / viability-like); binarize only where a specific downstream
  evaluation needs it. This is why the target is a continuous regression
  ([Step 03](03-model-and-training-design.md)) and why binary clinical outcomes are a *fine-tuning* step
  rather than the training objective.
- **Multi-task can absorb missing labels via masked losses** — no need to force a full intersection. The
  direct origin of the mask `M` ([Step 03](03-model-and-training-design.md#mask-m--the-sparsity-handling-mechanism-plan-sub-goal-2)).
- **Output/task weighting may be needed during training** — flagged this early; still an open question.
- **Report overlap and applicable sample counts before modelling**, to judge feasibility — which is what
  the audit below was built to do.
- **DrugBank** is acceptable for FDA/drug annotation support.
- **scDrugAtlas** is usable if the data can be obtained, but it is Harmony-processed and cross-dataset
  merging should be avoided — a caution that contributed to abandoning it
  ([Corrections](corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources)).

> **On the word "toxicity".** Finding usable toxicity definitions and annotations is genuinely hard, and
> the reason is structural: toxicity toward tumour cells is precisely what is wanted, while *excessive*
> toxicity is what withdraws a compound from trials — so the quantity is recorded inconsistently and
> mostly as an absence. That is why this project operationalises toxicity as cell-line efficacy, and why
> PRISM's coverage of failed drugs was attractive.

> **One GDSC-derived artifact exists and is not part of the modelling work.** An initial
> informative-drug list was produced from `notebooks/analysis/harmonization/drug_coverage.ipynb` and shared
> on request outside this project — a CTRPv2 version and a GDSC version
> (`notebooks/outputs/legacy/gdsc_drug_learnability.csv`, alongside
> `ctrp_drug_learnability_mean_pv.csv`; both moved to `legacy/` when the outputs were reorganized). Both
> were shared as **explicitly not final**, and the GDSC one for that use only; no GDSC drug list has
> ever fed this project's drug selection or training. GDSC remains downloaded-but-unused, and is not a
> modelling priority.

### Licences and terms of use of the source data (checked 28.07.2026)

The repository's own MIT licence covers **code and documentation only** — it grants nothing over the
data analysed here, each source of which carries its own terms. Checked against the providers:

| Source | Terms | Redistribution | Commercial use |
|---|---|---|---|
| **CTRPv2** (CTD² / Broad, via DepMap) | **CC BY 4.0** | permitted with attribution | permitted with attribution |
| **PRISM** Repurposing (DepMap) | **CC BY 4.0** | permitted with attribution | permitted with attribution |
| **SCP542** (Broad Single Cell Portal) | **no named licence, confirmed at both levels.** The portal's Terms of Service state that data in *public* studies is available for "unrestricted public view, redistribution and reuse"; the portal does not own the data, the contributing study does — and the study names no licence and declares no restriction either (checked 10.08.2026, see above). The same data is deposited at GEO under **GSE157220** | permitted per the ToS | not addressed |
| **GDSC2 / Cell Model Passports** (Sanger) | ⚠️ **no open licence.** A bespoke policy grants "a non-exclusive, non-transferable right to use data files for **internal** proprietary research and educational purposes", and explicitly excludes resale, combination with other data or product offerings, and provision of commercial services | **not granted** | **excluded** |
| **DrugBank** (compound harmonization only) | **CC BY-NC 4.0** for the academic full download; account required. Only the separate *DrugBank Open Data* identifier set is CC0 | permitted, non-commercial | **requires a separate agreement** |

**The operative asymmetry: GDSC is the only source that does not grant redistribution, and DrugBank is
the only one restricted to non-commercial use.** Both feed
`data/drug/all_sources_drug_catalog.csv` — 295 GDSC rows with targets and pathways, plus three DrugBank
match files.

**No violation exists today** — verified 10.08.2026: `.gitignore` line 1 excludes `data/`, and
`git ls-files data/` returns nothing. Note this is the **repo-local** root, which is where both
restricted sources happen to live (GDSC2 and the DrugBank-derived catalog); the Desktop root is outside
the working tree and cannot be committed at all. But that file is exactly the one that would be
committed by accident during a cleanup, and the
repository's MIT licence would then appear to grant rights over GDSC and DrugBank content that neither
provider allows. **Keep `data/` untracked**, and if a harmonized compound table ever needs to be shared,
share the CTRPv2 and PRISM columns only, or regenerate it from the providers.

Attribution obligations that follow: cite Kinker et al. 2020 for SCP542, the CTRP publications for
CTRPv2 ([Step 05](05-multitask-results.md) carries them), and acknowledge the CTD² Data Portal URL where
the funding is acknowledged. Sources for the table above are recorded in `references.bib`.

✅ *Checked 10.08.2026:* SCP542 carries **no** study-level licence of its own. The portal delegates to
the contributing study, and the study's own Data availability statement (Kinker et al. 2020, p. 13)
names no licence and declares no restriction — see the
[upstream QC section](#upstream-qc--what-scp542-already-had-done-to-it-05082026) for the quoted text.
The portal ToS stands as the operative term, and the FAIRER **F** item is closed.

---

## Overlap & coverage audit (03.04.2026)

The audit lives in the standalone notebook `notebooks/analysis/harmonization/drug_catalog.ipynb` — a one-off
exploratory analysis (not part of the training pipeline) that produces the plan's **Fig. 1 / Fig. 2**
and writes the drug-catalog CSVs. Its purpose is to pick the **highest-confidence intersection** to
start from before any modeling (plan sub-goal 3).

**Cell-line overlap with SCP542** (the notebook's normalization: trim + lowercase; SCP542/PRISM
split on `_`):

| Dataset | Total lines | Overlap w/ SCP542 | Missing |
|---|---|---|---|
| GDSC | 967 | 133 | 65 |
| CTRPv2 | 1,107 | **190** | 8 |
| PRISM | 915 | 182 | 16 |

**Drug / compound harmonization** — the notebook builds a unified catalog
`data/drug/all_sources_drug_catalog.csv` (7,415 source rows; GDSC 295 / CTRPv2 545 / PRISM 6,575;
union 7,040) by matching compounds three ways, in increasing confidence:

- **Normalized name** overlap: CTRPv2↔GDSC 66, CTRPv2↔PRISM 218, GDSC↔PRISM 144.
- **BRD-ID** overlap CTRPv2↔PRISM: **243** — the higher-confidence link, because Broad **BRD
  identifiers** are canonical per compound (stable across name synonyms and salt forms), unlike
  free-text names. Candidate pairs exported to `data/drug/drug_overlap_candidates.csv`.
- **DrugBank** name+synonym match (from `full database.xml`), which additionally enables future
  FDA-status filtering (plan sub-goal 1): GDSC 118/295 (40.0 %), CTRPv2 173/545 (31.7 %),
  PRISM 3,483/6,575 (53.0 %), overall 3,774/7,415 (50.9 %). Exports:
  `drugbank_overlap_matches.csv`, `drugbank_overlap_unmatched.csv`.

**Applicable (non-null) response coverage within the SCP542 overlap** — this is what motivates
choosing CTRPv2 as the starting database:

| Dataset | Metric | Non-null / total | % | % within overlap subset |
|---|---|---|---|---|
| GDSC | `LN_IC50` | 8,007 / 242,036 | 3.31 % | 100 % |
| CTRPv2 | `cpd_avg_pv` | 1,521,028 / 7,227,951 | 21.04 % | **100 %** |
| PRISM | extended primary matrix | 1,210,432 / 4,213,048 | 28.73 % | 97.95 % |

The decisive number is CTRPv2's **100 % non-null within the SCP542 overlap** — by far the densest of
the three, which is what motivated starting there.

⚠️ **But that 100 % is `cpd_avg_pv`, the raw dose grid, and it is not the density of the target.** Every
pair that was screened has dose points; a **curve-fit** measure additionally drops pairs whose fit
failed QC. The trainable matrix is **84.7 %** dense for `auc_cc` and **64.5 %** for `ln_ic50_cc`
([below](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)). This section previously called it a
"complete target matrix" — that has been wrong since the target became a curve fit on 27.07.2026, and
is corrected here 11.08.2026. The comparison between the three databases stands; the word *complete*
does not.

✅ On-plan: satisfies sub-goal 1 (harmonization incl. BRD + DrugBank) and supplies the Fig. 1/2
numbers sub-goal 3 rests on.

---

## What actually feeds training vs. what is forward-looking

The audit above is **exploratory** — its catalogs (`all_sources_drug_catalog.csv`,
`drug_overlap_candidates.csv`, the DrugBank exports) are **not yet consumed by any model**. They are
built for the cross-database join in [Step 06](06-planned-work.md#a-cross-database-integration), where PRISM/GDSC
heads finally need a unified compound vocabulary.

Today the trained model depends on Step 01 through exactly **one** thing: the cell-line and drug
**name normalization inside the pipeline**, `_normalize_cell_line` / `_normalize_drug` in
`scripts/preprocessing/ctrp_to_h5ad.py` (trim + lowercase + **strip `-`**). These produce the
`ccl_name_norm` / `cpd_name_norm` join keys that map CTRPv2 response scores onto SCP542 cells during
the **targets** step ([Step 02](02-preprocessing-and-embeddings.md)). At pipeline run time the overlap
is **180**, not the audit's 190 — and the reason is **data availability, not normalization**:

> ✅ **190 vs 180 vs 181 — the whole funnel (resolved 14.06.2026, extended 11.08.2026).** Both
> normalizations (with/without stripping `-`) give **190** SCP542 names present in CTRPv2's cell-line
> roster, but only **180** of those were actually **screened**. The 10 roster-listed-but-unscreened
> lines are `abc1, hs939t, jhh7, mdamb436, mfe280, ncih1048, ncih2073, ncih2347, rerflckj, ten`.
> The pipeline then recovers one more via the `h292 → ncih292` alias
> ([below](#the-join-dropped-a-screened-cell-line-h292-10082026)), giving the **181** it produces today.
> **190 is only a name-match count and is never the trainable set.**
>
> ⚠️ Every committed artifact, and `\NLines` in `report/results_numbers.tex`, still rests on **180**.
> Code and artifacts agree again only after the [clean sweep](../TODO.md).

### The target moved to DrEval's reprocessed CTRPv2 (11.08.2026)

**Decided 11.08.2026 (Selin).** The training target is no longer derived from CTRPv2's own 2015
distribution. It now comes from **DrEval's reprocessing of the same underlying screen**, pinned to
Zenodo record [`21807175`](https://doi.org/10.5281/zenodo.21807175) and fetched by
`scripts/sources/fetch_ctrp_response.py` (the `fetch` step, [Step 02](02-preprocessing-and-embeddings.md)).

**What DrEval do differently.** They download CTRPv2's *raw* dose-response measurements, normalise each
replicate against its own no-drug control, and fit one curve across replicates with CurveCurator —
rather than averaging replicates first, which is what CTRP did and what the field does generally:

> "Instead of aggregating replicates prior to normalization and curve fitting, **as is the standard
> practice**, DrEval includes replicate variability into quality control measures. This source of
> experimental variability is often overlooked when aggregating replicates prior to fitting, which
> leads to **inaccurate or misleading drug response measures** in the case of large discrepancies
> between replicates."
>
> — DrEval, Methods, "Benchmark data" (`papers/DrEval_s41467-026-72903-w.pdf`). Curve fitting via
> CurveCurator; datasets on Zenodo, concept DOI `10.5281/zenodo.12633909`.

They do **not** criticise CTRPv2's normalisation, and CTRP's published `area_under_curve` is not
mentioned in their paper — their stated motivation is cross-dataset consistency. The reason we
switched is separate and is recorded as an error, not an improvement: our own use of that column was
defective ([corrections](corrections-and-dead-ends.md#the-auc-target-was-divided-by-the-wrong-quantity)).

**The two measures now available**, both columns of the same CurveCurator fit
(`layout.CTRP_SCORES`, `ctrp_to_h5ad.SCORE_COLUMNS`):

| | column | direction | completeness |
|---|---|---|---|
| **`auc_cc`** (default) | `AUC_curvecurator` | **higher = more resistant**; 1.0 is the no-effect level, because the fit pins the low-concentration asymptote to the vehicle value | every curve |
| `ln_ic50_cc` | `LN_IC50_curvecurator` | **lower = more sensitive** — the opposite direction | 59.7 % of curves |

`ln_ic50_cc` is incomplete by construction, not by accident: DrEval discard an IC50 falling more than
an order of magnitude outside the measured dose range, which is most compounds that never reach
half-killing. **CTRPv2 itself publishes no IC50 at all** — its curve table carries `apparent_ec50_umol`
(50 % of the *fitted decline*, not of absolute viability), so `ln_ic50_cc` exists only in the
reprocessing.

**The sparsity is not spread evenly, and what it tracks is potency.** A curve that never crosses 50 %
viability has no IC50 to report, so a compound's IC50 count is a proxy for how hard it kills. Across the
181 overlapping lines the range is extreme — `crizotinib` 174 and `dasatinib` 172, against `decitabine`
5, `temozolomide` 5, `procarbazine` 4, `thalidomide` 3, `ifosfamide` 3 and `fulvestrant` 0 — and the
low-count compounds are cytostatic or slow-acting rather than badly measured. **This is why IC50
coverage is not a selection criterion** (see the panel below): filtering on it would re-create the
potency filter that voided the
[learnability gate](corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability),
under a name that sounds like a data-quality rule. Restricting an `auc_cc`-vs-`ln_ic50_cc` comparison to
drugs with enough IC50s is legitimate; doing it at selection time is not.

**What the target build produces.** Counts from
`notebooks/analysis/harmonization/cell_line_join_verification.ipynb` §3, which calls the pipeline's own
loader, de-duplication and drug filter rather than reimplementing them:

| measure | lines | drugs | observed | density | min | median | max |
|---|---|---|---|---|---|---|---|
| `auc_cc` | 181 | 534 | 81,906 | 84.7 % | 0.020 | 0.925 | 1.830 |
| `ln_ic50_cc` | 181 | 365 | 42,584 | 64.5 % | −11.385 | 2.487 | 8.634 |

**Why the counts differ from CTRPv2's own**, same notebook section:

| observation | why |
|---|---|
| 886 cell lines, not the roster's 1,107 | 887 lines have at least one curve fit; DrEval's re-fit keeps 886 — only `KRIJ` is absent |
| 395,024 rows against CTRP's 395,263 post-QC fits | **not a subset.** 2,668 (line, drug) pairs exist only in CTRP's post-QC set and 2,375 only in DrEval's, with 384,462 in common. The difference is *symmetric* because DrEval re-fit under CurveCurator's quality control instead of adopting CTRP's — so the two are not the same curve set, which is a consequence of the switch rather than a fault in it |
| 8,187 duplicate rows for `auc_cc`, 4,346 for `ln_ic50_cc` | `CTRPv2.csv` repeats rows: 15,946 of them duplicate another row across **all 46 columns**, in 7,331 pairs and 428 triples. Dropping the rows with no IC50 first removes many partners, hence the smaller count. The pipeline drops exact duplicates and **raises** on any pair that disagrees — none does |
| 534 of 545 drugs (`auc_cc`) | 11 drugs reach fewer than 50 of the 181 overlapping lines |
| 365 of **539** (`ln_ic50_cc`) | the denominator is not 545: six drugs have no valid IC50 anywhere in the overlap |
| 181 overlapping lines, where the raw name join gives 180 | the pipeline applies the `h292 → ncih292` alias ([below](#the-join-dropped-a-screened-cell-line-h292-10082026)) |
| `auc_cc` tops out at 1.830 where the old `auc` reached 2.310 | the old upper tail was substantially the divisor artifact |

**No winsorization and no quality filter (decided 11.08.2026, Selin).** The target is used as
published: nothing is clipped, and no row is dropped on fit quality, although `CTRPv2.csv` ships
`R2`, `RMSE`, `pValue` and `SignalQuality` per curve. This follows the benchmark rather than inventing
a rule:

> "The datasets can be filtered for quality using the statistical measures provided by CurveCurator
> … If not stated differently, **we did not apply any quality filter in our benchmark experiments to
> maintain comparability to previous studies and avoid data loss**."
>
> — DrEval, Methods, "Benchmark data"

**What this retires.** The pipeline previously carried `DEFAULT_WINSOR = 1.1` — clip the response at
1.1 — which had no source beyond a code comment, and whose stated purpose was to stop inverse-density
loss weighting from handing the sparse upper tail the largest weights. Three things ended it: the
threshold was never sourced; that weighting is itself a
[refuted hypothesis](corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking);
and most of the old upper tail was the divisor artifact rather than data. On `auc_cc` the tail is
thin — **3.65 %** of measurements exceed 1.1, 0.82 % exceed 1.2, p99 = 1.185
(`cell_line_join_verification.ipynb` §3 builds the distribution). If audit 09 keeps density weighting,
whether to clip the weighting's *input* — never the target — is reopened there.

**Also gained.** Every row carries a **Cellosaurus accession**, and `meta.zip` ships Cellosaurus 52.0,
so the name join can now be checked against an external authority — see the join audit below.

### The join audit — what was checked and what held (10.08.2026)

Walked as review item 2. The name join is the only thing linking CTRPv2 to SCP542, and until this date
it had never been checked for the failure mode that matters most: matching two lines that are *not* the
same. It does not.

| Check | Result |
|---|---|
| Normalized-name collisions, SCP542 | none — 198 names → 198 distinct keys |
| Normalized-name collisions, CTRPv2 roster | none — 1,107 → 1,107 |
| Normalized-name collisions, CTRPv2 compounds | none — 545 → 545 |
| `Cell_line.split("_")[0]` truncating a line name | never — 58 SCP542 names carry more than one `_`, all of it in the CCLE tissue suffix |
| Curve-fit rows lost to the three inner merges | none |
| **False matches** | **none** — 189 of the 190 name matches also agree on CCLE primary site (SCP542's own name suffix vs CTRPv2 `ccle_primary_site`); the single exception is `EKVX`, where CTRPv2 records no site at all |

Reproduced from the shipped files by
`scripts/preprocessing/ctrp_to_h5ad.py::_load_ctrp_long` plus the checks in review item 2; the two
defects it did find are below.

### The join dropped a screened cell line: `H292` (10.08.2026)

CTRPv2 spells one line differently from CCLE, and since the **name is the only join key**, the line was
invisible to the pipeline although it had been screened:

| | Name | Normalized key |
|---|---|---|
| SCP542 (CCLE naming) | `NCIH292_LUNG` | `ncih292` |
| CTRPv2 (`master_ccl_id` 290) | `H292` | `h292` |

The identification rests on three independent sources, all agreeing:

- **Cellosaurus `CVCL_0455`** is *NCI-H292*, listing **both** `H292` and `NCIH292` among its synonyms,
  disease *lung mucoepidermoid carcinoma*, a CCLE member, DepMap `ACH-001075`
  (<https://www.cellosaurus.org/CVCL_0455>, retrieved 10.08.2026).
- **CTRPv2's own row** agrees on every field it carries: `ccl_availability=ccle;public`,
  `ccle_primary_site=lung`, `ccle_hist_subtype_1=mucoepidermoid_carcinoma` — the **only** lung
  mucoepidermoid carcinoma among its 1,107 lines.
- **CTRPv2's naming is otherwise consistent**: 106 roster entries are written `NCIH…`, and this is the
  one place the prefix is dropped. (The only other bare `H<digit>` name, `H4`, is a genuinely different
  line — a CNS glioma, not an NCI-H line — which is why no general "strip a leading `NCI`" rule was
  adopted.)

**Fix (accepted by Selin, 10.08.2026):** an explicit alias table, `CTRP_CELL_LINE_ALIASES` in
`scripts/preprocessing/ctrp_to_h5ad.py`, applied to the CTRPv2 side inside `_load_ctrp_long`. The
evidence above is recorded in the table itself, so a future entry cannot be added without its own.
**Effect at the next sweep:** trainable overlap **180 → 181**, recovering **213 cells** and **454 drug
labels**. `\NLines` in `report/results_numbers.tex` still reads 180 and must not be changed until an
artifact supports 181.

Eight SCP542 lines matched no CTRPv2 name at all. `ncih292` was this defect; the other **seven** are
genuinely absent from CTRPv2 and stay unlabelled: `93vu, jhu006, jhu011, jhu029, ncih2077, scc47,
scc90` — the JHU and SCC head-and-neck lines were never in the CTRPv2 panel. With
the ten roster-listed-but-unscreened lines that leaves **17 of 198** SCP542 lines unlabelled — **6,073
cells**, which `scripts/preprocessing/create_splits.py` excludes from every split via `has_any_label`,
but which remain in the h5ad and therefore still enter HVG selection, scaling and the PCA fit
([Step 02](02-preprocessing-and-embeddings.md)).

### The drug-name join, and the compounds it hid (12.08.2026)

Walked as review item 6. The **cell-line** join was audited on 10.08.2026; the **drug** join never had
been, and it is worse. Two files describe the same 545 compounds — DrEval's `CTRPv2.csv` holds the
response values, and `data/drug/all_sources_drug_catalog.csv` (built from CTRP's own
`v20.meta.per_compound.txt` in `notebooks/analysis/harmonization/drug_catalog.ipynb`) holds approval
status, protein target and mechanism. Every consumer joined them **by name**.

**102 of 545 do not match by name**, because DrEval renamed to preferred names, changed separators
(`:` → `-`) and altered hyphenation. **15 of the unmatched are single-agent and FDA-approved or in
clinical trials** — exactly what a clinically-motivated selection looks for:

| CTRP name | DrEval name | | CTRP name | DrEval name |
|---|---|---|---|---|
| `abt-199` | Venetoclax | | `byl-719` | Alpelisib |
| `sirolimus` | Rapamycin | | `cal-101` | Idelalisib |
| `fluorouracil` | 5-Fluorouracil | | `gdc-0941` | Pictilisib |
| `mitomycin` | Mitomycin-C | | `lbh-589` | Panobinostat |
| `alvocidib` | Flavopiridol | | `nvp-bez235` | Dactolisib |

(plus `ex-527`, `tg-101348`, `tipifarnib-p1`, `vx-680`, `ym-155`.) Two of them, `sirolimus` and
`gdc-0941`, are compounds the [voided literature panel](corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it)
had already lost once.

**Fix: join on `master_cpd_id`.** CTRP's own compound identifier is present in both files — the
catalog's `identifier`, with `source_id_type == "master_cpd_id"` — and matches **545 of 545 in both
directions**. `scripts/annotation/drug_annotation.py` is now the single implementation and **raises**
if the join is ever inexact, rather than dropping rows.

Note this is the **opposite conclusion to the cell-line join**, where names resolve 180 SCP542 lines and
Cellosaurus accessions only 172, so the accession rides along as an attribute and the name stays the key
([above](#the-join-audit--what-was-checked-and-what-held-10082026)). Neither "always prefer the
identifier" nor "always prefer the name" would have got both right; each join was decided on its own
measured coverage.

#### Cisplatin was invisible: a wrong structure on a nonstandard name

CTRPv2 screens cisplatin under the name **`Platin`** (`master_cpd_id` 375395), and DrEval's
`pubchem_id` for it is **23939 — elemental platinum**, the element rather than the drug. With a
nonstandard name *and* a wrong structure, the most widely used platinum agent in oncology could not be
found by any key an external drug list can use.

The identification rests on CTRP's own compound record: SMILES `N[Pt](N)(Cl)Cl`, i.e.
*cis*-diamminedichloroplatinum(II), and vendor entry **Selleck S1166**, which is cisplatin.

**Fix:** `CTRP_PUBCHEM_OVERRIDES` in `drug_annotation.py` maps 375395 → PubChem CID **5702198**,
explicitly and one compound at a time, with the evidence in the table — the same pattern as
`CTRP_CELL_LINE_ALIASES`. The panel refers to the drug by the key the data uses, `platin`, so no
downstream translation table exists to fall out of date.

### The drug panel — FDA-approved compounds this screen covers (12.08.2026)

**Decided 12.08.2026 (Selin), rebuilding review item 6.** Both previous panels were selected using our
own response values and both were voided for it — the
[learnability gate](corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability),
which filtered on absolute potency and so discarded every cytostatic compound, and the
[8-drug literature panel](corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it),
whose candidate list was ranked on our AUCs before any citation was consulted. The rebuild therefore
takes the criterion **entirely outside our labels**.

Produced end to end by `notebooks/2_drug_selection.ipynb` →
`notebooks/outputs/panel/panel.csv`. It reads no pipeline artifact, only the response CSV, so it runs
before the sweep and does not go stale when the h5ads do.

**The panel does not enter the target build.** `Y_ctrp` / `M_ctrp` keep the full screened catalog and the
panel is applied at training time, because what it determines is the number of heads rather than what the
data contains (Selin, 12.08.2026):
[Step 03](03-model-and-training-design.md#the-drug-panel-is-a-training-time-choice-not-a-property-of-the-target-file-12082026).

**The criterion, in three conditions.**

1. **FDA-approved for a cancer indication**, from Sun, J., Wei, Q., Zhou, Y., Wang, J., Liu, Q. & Xu, H.
   *A systematic analysis of FDA-approved anticancer drugs.* **BMC Systems Biology 11**(Suppl 5), 87
   (2017), [doi:10.1186/s12918-017-0464-7](https://doi.org/10.1186/s12918-017-0464-7), Table 1 — 150
   drugs approved 1949–2014 with approval year, therapeutic class, target gene and delivery type.
   Retrieved by `scripts/sources/fetch_sun2017_drugs.py`; see [provenance](#provenance--what-was-retrieved-from-where-when).
2. **Screened by CTRPv2 against ≥ 90 % of the 181 overlapping cell lines.** The cut is where the
   coverage distribution breaks, and the break is sharp: 45 candidates sit at or above it, the lowest
   being `afatinib` at 91.2 %, and the next compound down is `omacetaxine mepesuccinate` at 69.1 % — a
   22-point gap — with the tail reaching `fulvestrant` at 28.7 %.
3. **Carries a published claim about it**, recorded per drug with its reference and what that
   reference actually establishes.

Conditions 1 and 2 are derived in code. Condition 3 is a literature judgement and cannot be, so it is
carried as a table inside the notebook — the pattern the voided panel's `DETERMINANTS` established — and
the code asserts every entry lies inside the set that conditions 1 and 2 produce.

**The funnel.** 150 FDA-approved drugs → **120** with a PubChem structure (30 are biologics — antibodies,
enzymes, a cell therapy, a radiopharmaceutical — which cannot appear in a small-molecule screen) →
**57** also screened by CTRPv2 → **45** at ≥ 90 % coverage → **11** with a verified published claim.

**Matching needed four keys, and each was necessary.** Sun names drugs by INN; CTRPv2 uses development
codes and salt forms. `drug_annotation.match_external_list` records per drug which key succeeded:
DrEval's spelling alone finds `idelalisib`; CTRP's alone finds `fluorouracil` and `mitomycin`; the
PubChem structure connects `Vemurafenib` to `plx-4032`; and **PubChem's parent-compound relation**,
applied to *both* sides, connects `Imatinib mesylate` to `imatinib` and, in the other direction,
`Cytarabine` to `cytarabine hydrochloride`. That last key alone accounts for 13 matches including
imatinib, doxorubicin, vincristine and topotecan — dropped, without it, on a suffix.

**The panel (11).**

| drug | key in the data | class | evidence |
|---|---|---|---|
| Gemcitabine | `gemcitabine` | cytotoxic | strong |
| Doxorubicin | `doxorubicin` | cytotoxic | medium |
| Etoposide | `etoposide` | cytotoxic | strong |
| Paclitaxel | `paclitaxel` | cytotoxic | medium |
| Cisplatin | `platin` | cytotoxic | **contested** |
| Imatinib | `imatinib` | targeted | strong |
| Sorafenib | `sorafenib` | targeted | strong |
| Dasatinib | `dasatinib` | targeted | medium |
| Erlotinib | `erlotinib` | targeted | medium |
| Crizotinib | `crizotinib` | targeted | medium |
| Afatinib | `afatinib` | targeted | medium |

Coverage runs 91.2 %–98.3 % of the 181 lines. Each drug's reference, the claim that reference makes and
the setting it was established in are columns of `panel.csv`; they are not restated here.

**Ten of the eleven are named in the four papers this project is built on** — Kinker 2020 (the
expression atlas), Seashore-Ludlow 2015 and Rees 2016 (the CTRPv2 papers) and scDEAL 2022 (the
single-cell response benchmark). Every mention was **read in the PDF**, not counted by a string search,
and two compounds a search had flagged were dropped on inspection: doxorubicin's Kinker "mentions" are
the names of published **senescence gene programs** ("lung cancer doxorubicin (n = 414)"), signatures
borrowed from other studies rather than a compound Kinker screened, and carboplatin's only appearance is
a title in scDEAL's bibliography.

**Only `doxorubicin` rests on outside evidence alone**, admitted to represent the anthracycline class,
which the others lack; the `setting` column in `panel.csv` records that. Cisplatin does **not** belong in
that category, though an earlier version of this page said so: it is one of the five drugs the scDEAL
benchmark is built on, and a perturbation in Kinker Fig. 6 tested for association with expression
heterogeneity alongside etoposide. What is contested for cisplatin is its *determinant*, not its
presence in the literature this project builds on.

**One platinum, deliberately.** CTRPv2 screens cisplatin, carboplatin and oxaliplatin, all at ~176 of
181 lines, so coverage does not choose between them. Cisplatin and carboplatin share the *cis*-diammine
carrier and form the same 1,2-intrastrand GpG adduct; they are largely cross-resistant and driven by the
same repair biology, so a second one would spend a panel slot on the same resistance mechanism.
Oxaliplatin's DACH–Pt adduct evades mismatch repair and is genuinely complementary, but its expression
determinant is the least established of the three. The class is represented **once**, by cisplatin.

**Cisplatin's determinant is contested, and this is recorded rather than smoothed over.** ERCC1 is the
standard marker of platinum resistance, and the two directly comparable cell-line studies disagree:

- Britten, R. A. *et al.* *Int J Cancer* **89**, 453–457 (2000) — ERCC1 **mRNA** correlates with
  cisplatin resistance in cervical carcinoma lines (*p* ≤ 0.011), while ERCC1 **protein** does not.
- Shimizu, J. *et al.* *Respirology* **13**, 510–517 (2008) — **no** correlation between ERCC1 mRNA and
  cisplatin or carboplatin sensitivity across 20 lung cancer lines.

Neither supersedes the other on the usual tie-breaks: the negative study is both the more recent (2008
vs 2000) and the far less cited (28 vs 155 citations, Semantic Scholar, 12.08.2026), and it is a
different lineage and assay. What explains the disagreement is Friboulet, L. *et al.*, *N Engl J Med*
**368**, 1101–1110 (2013), on 494 patients across two phase 3 trials: **ERCC1 has four splice isoforms
and only ERCC1-202 repairs DNA**, and available antibodies cannot distinguish them.

That lands directly on this project: our features are **gene-level CPM**, so the functional isoform is
not representable in them even in principle. If the model ranks cell lines poorly on cisplatin, this is
the first place to look — and it was known before the run rather than after it. A panel in which every
drug had a settled strong expression marker would be the stacked deck this rebuild exists to avoid.

**What is deliberately not in the criterion.** IC50 coverage — see the
[target section](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026): it tracks potency, so
selecting on it would rebuild the discredited gate. The `auc_cc`-versus-`ln_ic50_cc` comparison is
restricted at **evaluation** instead, on whatever subset of the panel has enough IC50s, with per-drug
counts reported alongside. Within the panel those run from cisplatin's 14 to crizotinib's 174.

*(A second write-up of the panel/target-build decision stood here from 12.08.2026 until the same day.
It gave split eligibility as the deciding reason, which is not the ruling: the panel stays out of the
target build because it determines how many heads the model has. Removed rather than corrected in
place, because the decision is already recorded once, in
[Step 03](03-model-and-training-design.md#the-drug-panel-is-a-training-time-choice-not-a-property-of-the-target-file-12082026),
and this file points at it above.)*

### Replicate experiments were double-counted (10.08.2026)

> ⛔ **Moot since 11.08.2026. The code this describes no longer exists.** `_load_ctrp_long` and the `v20.*` readers were removed
> on 11.08.2026 when the target moved to DrEval's reprocessing ([above](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)),
> which reads a single response table and never merges an experiment roster. The bug cannot recur, and
> the fix it describes is no longer in the tree. Kept because the numbers below quantify how much a
> silent duplication can move a *ranking* — the point that made it worth finding.

`v20.meta.per_experiment.txt` carries **one row per (`experiment_id`, `experiment_date`)**: 153 of its
907 experiments ran across two calendar days and therefore appear twice, 1,061 rows in total. They are
exact duplicates in every field the pipeline uses — `master_ccl_id` is constant within an
`experiment_id` — but `_load_ctrp_long` merged the curve fits on `experiment_id` **without dropping
them**, inflating 395,263 curve rows to 462,784, and the per-(line, drug) mean that follows is taken
over *rows*. A duplicated experiment therefore counted twice.

It bites only where a line has one duplicated experiment **and** a second, unduplicated one — exactly
one line, `NCIH1299`, on **460 of its 545 drugs**. The values moved by a median of 0.013 and at most
0.112, i.e. up to **0.95×** the mean per-drug spread across cell lines (0.117). Because the evaluation
metric is a per-drug Spearman *across* lines, the consequential number is not the 0.57 % of targets
touched but the ranking: `NCIH1299`'s rank moved on **427 of its 469** drugs, a median of **8** places
out of 180. Per-drug means and SDs shifted by ≤ 0.0014, so `auc_z`'s centering was never affected.

**Fixed 10.08.2026** by deduplicating the experiment table before the merge
(`ctrp_to_h5ad.py::_load_ctrp_long`); the loader now returns 395,263 rows, exactly the number of
post-QC curve fits. Applies identically to `auc`, `auc_z` and `mean_pv`. Like the alias above, it
changes no committed artifact until the sweep.

### Genuine repeats are averaged, and they disagree more than the target's own spread (10.08.2026)

> ✅ **Resolved 11.08.2026 — not superseded, and it is the strongest evidence for the target switch.** This section
> measured how far apart two screenings of the same (cell line, drug) fall, and found a typical
> disagreement of **0.49×** the spread the model is asked to predict. The old pipeline averaged them,
> which discards exactly that information. DrEval's reprocessing instead fits **one curve across the
> replicates** and folds their disagreement into the fit's quality measures — which is precisely the
> practice their Methods argue for, quoted [above](#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026).
>
> Confirmed on the data 11.08.2026: the six re-screened lines named below are among the **19** cell
> lines whose rows appear more than once in DrEval's `CTRPv2.csv`, and every such group is an **exact
> duplicate across all 46 columns**. A repeat therefore now arrives as one fitted value repeated per
> source experiment, not two values to be averaged. The pipeline drops the exact duplicates and
> **raises** if any pair ever disagrees, so the averaging this section warned about cannot silently
> return (`ctrp_to_h5ad._deduplicate_measurements`).
>
> The numbers below stand as the measurement of the problem; they are no longer a description of what
> the pipeline does.

Separate from the double-counting above: some (cell line, drug) combinations really were screened
**twice**, in two different experiments. `ctrp_to_h5ad.py::_build_drug_table` averages them into the
single value that becomes the target. How far apart those two measurements are was never examined
until now — quantified in `notebooks/archive/replicate_variation.ipynb`, artifacts
`notebooks/outputs/legacy/replicate_variation.{png,csv}`. **That notebook was archived on 11.08.2026**:
it reads CTRPv2's own `v20.*` tables, which are no longer the target source, so it can no longer be
re-run against what the pipeline uses. The numbers below are its output and stand as measured.

**2,637** of the 81,626 (cell line, drug) pairs were screened twice — never three times, so a median
would be identical to the mean. They are **not spread over the panel**: they come from just **six**
cell lines (`a375, aspc1, ccfsttg1, ncih1299, ncih460, oaw28`), each re-screened against 534 of the
545 drugs.

A raw difference in normalized AUC is not interpretable on its own — drugs differ by an order of
magnitude in how much they vary across cell lines — so each difference is also expressed **relative to
that drug's spread across the cell lines**, which is precisely the quantity the model is asked to
predict. A relative difference of 1.0 means the same line measured twice differs as much as two
randomly chosen lines typically do.

| | median | p90 | max |
|---|---|---|---|
| absolute \|rep1 − rep2\| | 0.053 | 0.205 | 0.903 |
| **relative to the drug's spread across lines** | **0.49×** | 1.75× | 7.01× |

Order statistics, not averages — a handful of extreme pairs would drag a mean upward and misrepresent
the typical case. **Median** = the middle pair: half the repeated pairs disagree by less than this and
half by more, so the headline **0.49×** reads as *"a typical repeat differs by about half the spread
the model is asked to predict"*. **p90** = only one repeated pair in ten disagrees by more than this.
**max** = the single worst pair, shown so the tail is visible rather than hidden.

**720 of the 2,637 repeats (27.3 %) differ by more than the full spread the model is trained to
predict.** It is not uniform across the six lines: `ccfsttg1` disagrees by a median 0.19×, `a375` by
0.77× and `ncih1299` by 0.74×. In the scatter (`replicate_variation.png`) the cloud is visibly loose
around the identity line, and loosest in the dense region at AUC ≈ 0.8–1.0 where most drugs barely
kill anything; the strong killers below 0.5 track the diagonal better.

> ✅ **Decision (Selin, 10.08.2026): keep averaging.** The repeats are averaged as before, with no
> weighting by replicate count and no exclusion of discordant pairs. What changes is that the
> disagreement is now measured and on the record rather than invisible.

**What this does and does not license.** It is an estimate of the screen's reproducibility on **six of
181** cell lines; extending it to the other 175 is an assumption, not a result — CTRPv2 chose which
lines to repeat, and that choice was not random. Read as an order of magnitude, it says a substantial
part of what a model is asked to predict is screening noise, so a modest per-drug ρ is partly a
property of the labels rather than of the representation. It does **not** support a numeric ceiling on
achievable ρ; that would need the repeats to be a random sample of the panel, and they are not.

From this overlap, a drug becomes a model **head** (one column of `obsm["Y_ctrp"]`, one row of the
output layer — never an input feature) only if it was screened on ≥ `--min-cell-lines` overlapping
lines; the headline run used `--all-drugs` (min 0) → **K = 545**
([Step 05](05-multitask-results.md)). The single-task work uses **paclitaxel** as its reference drug
([Step 04](04-single-task-results.md)). PRISM and GDSC are harmonized but **not wired into training**
— their integration is the open Phase-3b in the [scorecard](../project_progress.md).
