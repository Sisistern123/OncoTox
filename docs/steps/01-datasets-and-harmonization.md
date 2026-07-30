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
| **CTRPv2** (Cancer Therapeutics Response Portal v2) | dose-response labels | 1,107 cell lines, **545 compounds**; training target `auc_z` ([Step 03](03-model-and-training-design.md)) | ✅ primary |
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
  **efficacy** metric. *Which* of these becomes the label — and why the AUC fit wins — is
  [Step 03](03-model-and-training-design.md).
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
[TODO](../TODO.md) review item 1 and FAIRER **A**/**E** ask for.*

| Source | Retrieved | From |
|---|---|---|
| **SCP542** scRNA-seq (53,513 cells × 22,722 genes) | 30.03.2026 | Broad Single Cell Portal, study **SCP542** *"pan-cancer cell line heterogeneity"* — <https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity>. Source publication Kinker et al., *Nat Genet* 2020, <https://doi.org/10.1038/s41588-020-00726-6>. This is the dataset the PERCEPTION paper used. |
| **CTRPv2** viability + compound annotations | 30.03.2026 | DepMap portal, all-data tab — <https://depmap.org/portal/data_page/?tab=allData>. Release `CTRPv2.0_2015_ctd2_ExpandedDataset`. Cell lines join via `v20.meta.per_cell_line` (`master_ccl_id`) and experiments via `v20.meta.per_experiment`. |
| **PRISM** Repurposing, Public 24Q2 | 30.03.2026 | DepMap — `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`. Taken specifically because it covers **failed** drugs as well as approved ones. |
| **GDSC2** IC50 + raw data | 26.03–30.03.2026 | Sanger Cell Model Passports — <https://cellmodelpassports.sanger.ac.uk/downloads> ("GDSC2 IC50 Data", "GDSC2 Raw Data"). |
| **Clinical-phase annotations** | 30.03.2026 | Broad Drug Repurposing Hub — <https://repo-hub.broadinstitute.org/repurposing#download-data> |
| **DrugBank** full database (XML) | — | <https://go.drugbank.com/>, academic download; account required. Licence terms above. |

⚠️ **How GDSC2's `LN_IC50` was processed is not documented.** The Sanger documentation link
(`depmap.sanger.ac.uk/documentation/gdsc/`) was dead, and the DepMap/GDSC team was asked directly with no
response ([Corrections](corrections-and-dead-ends.md#scdrugatlas-and-clintox-as-data-sources)). This is an
open gap for any future GDSC head.

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
> informative-drug list was produced from `notebooks/data_and_harmonization/drug_coverage.ipynb` and shared with
> Hashimoto-san — a CTRPv2 version and a GDSC version (`outputs/data/gdsc_drug_learnability.csv`). Both
> were shared as **explicitly not final**, and the GDSC one was for her use only; no GDSC drug list has
> ever fed this project's drug selection or training. GDSC remains downloaded-but-unused, and is not a
> modelling priority.

### Licences and terms of use of the source data (checked 28.07.2026)

The repository's own MIT licence covers **code and documentation only** — it grants nothing over the
data analysed here, each source of which carries its own terms. Checked against the providers:

| Source | Terms | Redistribution | Commercial use |
|---|---|---|---|
| **CTRPv2** (CTD² / Broad, via DepMap) | **CC BY 4.0** | permitted with attribution | permitted with attribution |
| **PRISM** Repurposing (DepMap) | **CC BY 4.0** | permitted with attribution | permitted with attribution |
| **SCP542** (Broad Single Cell Portal) | **no named licence.** The portal's Terms of Service state that data in *public* studies is available for "unrestricted public view, redistribution and reuse"; the portal does not own the data, the contributing study does | permitted per the ToS | not addressed |
| **GDSC2 / Cell Model Passports** (Sanger) | ⚠️ **no open licence.** A bespoke policy grants "a non-exclusive, non-transferable right to use data files for **internal** proprietary research and educational purposes", and explicitly excludes resale, combination with other data or product offerings, and provision of commercial services | **not granted** | **excluded** |
| **DrugBank** (compound harmonization only) | **CC BY-NC 4.0** for the academic full download; account required. Only the separate *DrugBank Open Data* identifier set is CC0 | permitted, non-commercial | **requires a separate agreement** |

**The operative asymmetry: GDSC is the only source that does not grant redistribution, and DrugBank is
the only one restricted to non-commercial use.** Both feed
`data/drug/all_sources_drug_catalog.csv` — 295 GDSC rows with targets and pathways, plus three DrugBank
match files.

**No violation exists today**, because `data/` is excluded in `.gitignore` and nothing under it is
tracked. But that file is exactly the one that would be committed by accident during a cleanup, and the
repository's MIT licence would then appear to grant rights over GDSC and DrugBank content that neither
provider allows. **Keep `data/` untracked**, and if a harmonized compound table ever needs to be shared,
share the CTRPv2 and PRISM columns only, or regenerate it from the providers.

Attribution obligations that follow: cite Kinker et al. 2020 for SCP542, the CTRP publications for
CTRPv2 ([Step 05](05-multitask-results.md) carries them), and acknowledge the CTD² Data Portal URL where
the funding is acknowledged. Sources for the table above are recorded in `references.bib`.

*Not yet checked:* whether SCP542 carries a study-level licence of its own beyond the portal ToS — the
portal delegates that to the contributing study, so the Kinker et al. data-availability statement is the
place to look.

---

## Overlap & coverage audit (03.04.2026)

The audit lives in the standalone notebook `notebooks/data_and_harmonization/drug_catalog.ipynb` — a one-off
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

The decisive number is CTRPv2's **100 % non-null within the SCP542 overlap**: the 190-line × 545-drug
block is a **complete target matrix**, exactly the dense, highest-confidence intersection the plan
wants for the initial baseline.

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

> ✅ **190 vs 180 — resolved (14.06.2026).** Both normalizations (with/without stripping `-`) give
> **190** SCP542 names that appear in CTRPv2's cell-line roster (`v20.meta.per_cell_line.txt`). But
> only **180** of those were actually **screened post-QC** — identical whether counted from the raw
> dose grid or the curve fits, so the 13.07.2026 switch to `auc_z` did not change the trainable set.
> The **10**
> roster-listed-but-unscreened lines are `abc1, hs939t, jhh7, mdamb436, mfe280, ncih1048, ncih2073,
> ncih2347, rerflckj, ten`. Use **180** (the trainable set); 190 is just the name-match count.

From this overlap, a drug becomes a model **head** (one column of `obsm["Y_ctrp"]`, one row of the
output layer — never an input feature) only if it was screened on ≥ `--min-cell-lines` overlapping
lines; the headline run used `--all-drugs` (min 0) → **K = 545**
([Step 05](05-multitask-results.md)). The single-task work uses **paclitaxel** as its reference drug
([Step 04](04-single-task-results.md)). PRISM and GDSC are harmonized but **not wired into training**
— their integration is the open Phase-3b in the [scorecard](../project_progress.md).
