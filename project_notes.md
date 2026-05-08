
# OncoTox Project Notes
## 08.05.2026

### 1. Initial Setup & Random Splitting (The Data Leak)
* **Total cells:** 53,513
* **Cells with a valid paclitaxel viability score:** 44,367
* Added `ad.settings.allow_write_nullable_strings = True` to the split script to bypass AnnData string writing restrictions.
* Saved a **random 70/15/15 data split** directly into the `.h5ad` object (`split_paclitaxel` column).
    * **Train:** 31,056 cells | **Val:** 6,655 cells | **Test:** 6,656 cells | **Unassigned:** 9,146 cells

#### Training on Random Split (Data Leakage Identified)
* **scGPT Run:** Train MSE dropped to 0.0132, Val MSE to 0.0137.
* **PCA Baseline Run:** Train MSE dropped to 0.0022, Val MSE to 0.0011.
* *Conclusion:* The PCA baseline artificially outperformed scGPT because of data leakage. Since cells from the same cell line were randomly distributed across Train and Val, the PCA model simply memorized the tissue-of-origin "islands" rather than learning true biological resistance.

---

### 2. Fixing Data Leakage (Cell Line Grouped Split)
To prevent the model from cheating, the split logic was rewritten to group by cell line. If a cell line is in the training set, none of its cells appear in the validation set.
* Found 170 unique cell lines with paclitaxel labels.
* **Cell Line Split:** Train: 119 | Val: 25 | Test: 26
* **Final Cell Split Distribution:**
    * **Train:** 31,824 cells | **Val:** 5,035 cells | **Test:** 7,508 cells | **Unassigned:** 9,146 cells

#### Training on Grouped Split (Unregularized)
* **scGPT Run:** Train MSE 0.0110 | Val MSE: 0.0437
* **PCA Baseline Run:** Train MSE 0.0018 | Val MSE: 0.0390
* *Conclusion:* Data leak fixed. The PCA model completely failed to generalize to unseen cell lines (Val MSE stuck at ~0.04), proving it was previously cheating. However, scGPT also plateaued, indicating the highly parameterized MLP (256 hidden dims) was simply memorizing the noisy training labels.

---

### 3. Applying Aggressive Regularization (Proving the Hypothesis)
To force the MLP to learn generalized pathways instead of memorizing noise, aggressive regularization was applied:
* Reduced `hidden_dim` from 256 to 64.
* Increased `dropout_rate` from 0.3 to 0.5.
* Increased Adam optimizer `weight_decay` (L2 regularization) from 1e-5 to 1e-3.

#### Final Regularized Runs
**scGPT (Regularized)**
* Epoch [01/50] | Train MSE: 0.1014 | Val MSE: 0.0490
* Epoch [10/50] | Train MSE: 0.0260 | Val MSE: 0.0371
* Epoch [50/50] | Train MSE: 0.0260 | Val MSE: 0.0391
* *Gap:* ~0.013

**PCA Baseline (Regularized)**
* Epoch [01/50] | Train MSE: 0.0790 | Val MSE: 0.0480
* Epoch [10/50] | Train MSE: 0.0090 | Val MSE: 0.0375
* Epoch [50/50] | Train MSE: 0.0082 | Val MSE: 0.0380
* *Gap:* ~0.029

**Final Conclusion:** The foundation model successfully demonstrated a mathematically superior prior. While absolute validation performance hit a "weak supervision wall" around ~0.037 for both, the scGPT representation vastly reduced overfitting (Train/Val gap of 0.013 vs the baseline's 0.029). This proves standard PCA relies heavily on memorization, whereas scGPT forces the network to learn generalized, cross-tissue biological signatures.


## 21.04.2026 - 07.05.2026
### Preprocessing Timeline
scp542_conversion.py: massive, raw CPM_data.txt and Metadata.txt files and compiled them into foundational SCP542_CCLE.h5ad object.

scGPT Execution (The Embedding Phase): The scGPT foundation model was run using that raw .h5ad file to generate the continuous biological prior, outputting SCP542_CCLE_scGPT_human_embeddings.h5ad.

ctrp_to_h5ad.py (The Target Mapping): loaded those embeddings, parsed the CTRPv2 databases, merged the metadata, and cleanly mapped the paclitaxel viability scores to specific cell lines. This generated the final, all-in-one file: SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad.

The UMAP Script (Latent Space Validation): loaded that final file into memory, calculated the standard PCA baselines, ran the UMAPs, and generated the comparative visual proof

## 20.04.2026
### thoughts about available scores
#### GDSC2 (Genomics of Drug Sensitivity in Cancer)
Focus: Continuous dose-response profiling across a wide concentration range.
- LN_IC50 (Natural Log of Half-Maximal Inhibitory Concentration): The concentration required to kill 50% of the cells. Lower values mean higher toxicity to the cancer cell line.
- AUC (Area Under the Curve): The integrated area under the dose-response curve. It captures both the potency (how much drug is needed) and efficacy (the maximum kill rate). Lower values indicate strong sensitivity.
Best Pick for Multi-task: AUC is generally preferred over IC50 for machine learning because it captures the entire behavior of the curve, whereas IC50 can be noisy if the drug never actually reaches a 50% kill rate (resulting in extrapolated, artificial values).

#### CTRPv2 (Cancer Therapeutics Response Portal)
Focus: Similar to GDSC, this focuses on multi-dose response curves, but uses different curve-fitting algorithms and raw viability measurements.
- area_under_curve: As highlighted in your data dictionary, this is the integrated area under their sigmoid-fit concentration-response curve. Lower values equal higher sensitivity.
- cpd_avg_pv (Compound Average Percent Viability): The weighted average of surviving cells across all tested doses.
Best Pick for Multi-task: area_under_curve. Using AUC here allows your CTRPv2 head to learn a conceptually similar target to your GDSC head, even though the raw scales differ.

#### PRISM Repurposing (Public 24Q2)
Focus: High-throughput, single-dose screening. According to your provided Readme, all compounds here were screened at a single dose of 2.5 μM. Because there is no dose curve, there is no AUC or IC50.
- LMFI (Log2 Median Fluorescence Intensity): The raw optical readout of the barcodes.
- LFC (Log2 Fold Change): The median collapsed log-ratio of treated cells versus negative control (DMSO) cells. A negative LFC means the cancer cell line was depleted (killed) by the drug.
Best Pick for Multi-task: LFC (specifically from the Extended_Primary_Data_Matrix.csv). This is your only true efficacy metric for this dataset.

Multi-Task Recommendations
multi-task neural network with three output heads, your best target matrix configuration would be:

Head 1 (GDSC2): AUC
Head 2 (CTRPv2): area_under_curve
Head 3 (PRISM): LFC



hard to actually find toxicity definition and annotations, since toxicity is wanted for cancer treatments but excessive toxicity results in withdrawing from trials.

## 06.04.2026
### Ideas
- train on scRNA-seq data first and fine-tune on specific cancer types?
  - depends on available datasets
- train on scRNA-seq data first and fine-tune on clinical data?
  - depends on labels chosen for training?
  - depends on finding sufficient and properly labeled clinical data
- combine the two above ideas
  - find clinical scRNA-seq data sets for specific cancer types
  - maybe do bulk or pseudo bulk RNA-seq training as well, since scRNA-seq data is sparse?
- make a scRNA-seq transformer to use as a pretrained model? VAE? something else entirely?
- use scGPT embeddings (https://www.nature.com/articles/s41592-024-02201-0)
- look how (sc) DeepInsight could work with this


## 03.04.2026
### Advisor updates and alignment
- Clarification from Artem:
  - toxicity/response definition should be driven by available labels in selected datasets
  - if possible, include multiple response types instead of a single endpoint
  - prefer continuous outputs (IC50/viability-like) and binarize only for specific evaluations
  - prepare overlap and applicable sample counts to judge feasibility before modeling
- Methodology guidance:
  - DrugBank can be used for FDA/drug annotation support: https://go.drugbank.com/
  - multi-task setup can handle missing labels via masked losses (no need to force full intersection)
  - output/task weighting may be needed during training

### Notebook work completed (`notebooks/compare_GDSC_CTRP.ipynb`)
#### 1) Cell-line overlap with SCP542 (case-insensitive, unique names)
- SCP542 unique cell lines: `198`
- GDSC unique cell lines: `967`; overlap with SCP542: `133`; missing from GDSC: `65`
- CTRPv2 unique cell lines: `1107`; overlap with SCP542: `190`; missing from CTRPv2: `8`
- PRISM unique cell lines: `915`; overlap with SCP542: `182`; missing from PRISM: `16`
- Matching used normalized names (trim + lowercase; SCP542/PRISM names split on `_` where applicable).

#### 2) Drug/compound harmonization and cross-dataset overlap
- Built unified catalog: `data/drug/all_sources_drug_catalog.csv`
  - source rows: GDSC `295`, CTRPv2 `545`, PRISM `6575`
  - total rows: `7415`
- Name-based unique compounds:
  - GDSC `286`, CTRPv2 `545`, PRISM `6575`
  - overall union: `7040`
- Pairwise name overlap (normalized name):
  - CTRPv2 vs GDSC: `66`
  - CTRPv2 vs PRISM: `218`
  - GDSC vs PRISM: `144`
- Added CTRPv2<->PRISM BRD-based matching (from canonical `BRD-*` IDs):
  - BRD overlap: `243` (higher-confidence link)
- Exported overlap candidates:
  - `data/drug/drug_overlap_candidates.csv`

#### 3) DrugBank overlap
- Loaded DrugBank XML (`/Users/selin/Desktop/OncoTox/full database.xml`) and matched to catalog via normalized names (+ synonym expansion where available).
- Dataset-level matches to DrugBank:
  - GDSC: `118 / 295` (`40.00%`)
  - CTRPv2: `173 / 545` (`31.74%`)
  - PRISM: `3483 / 6575` (`52.97%`)
  - overall: `3774 / 7415` (`50.9%`)
- Exported review files:
  - `data/drug/drugbank_overlap_matches.csv`
  - `data/drug/drugbank_overlap_unmatched.csv`

#### 4) Applicable sample numbers (non-null response values)
- Computed for SCP542-overlapping cell lines vs each dataset's own full response space:
  - GDSC (`LN_IC50`): `8,007 / 242,036` (`3.31%`)
  - CTRPv2 (`cpd_avg_pv` viability): `1,521,028 / 7,227,951` (`21.04%`)
  - PRISM (Extended Primary matrix values): `1,210,432 / 4,213,048` (`28.73%`)
- Added second completeness metric in notebook:
  - `pct_non_null_within_overlap_subset` = non-null overlap rows / total overlap rows
  - Current values: GDSC `100%`, CTRPv2 `100%`, PRISM `97.95%`

### Interpretation notes
- "Applicable sample numbers" are comparable as response coverage indicators, but raw units differ:
  - GDSC/CTRPv2: long-table response rows
  - PRISM: matrix entries (compound x cell line)
- Therefore percentages are best interpreted as coverage estimates, not strict one-to-one sample equivalence.

### Cleanup/refactor performed
- Refactored repetitive SCP542 comparison code into reusable notebook helpers.
- Removed trailing empty notebook cells and tiny redundant summary CSV files.

### Pending decisions / next actions
- Decide primary training endpoint strategy:
  - option A: single dataset first (higher consistency)
  - option B: multi-task with masked missing labels (higher coverage)
- Decide whether to harmonize drug identity using:
  - conservative set (name + BRD agreement), or
  - broader set (name/synonym candidates for manual curation).
- Prepare a compact slide/table for in-person discussion with Artem next week:
  - cell-line overlap, response coverage, and recommended modeling path.

## 31.03.2026
### Advisor response to project-definition questions
- Project framing:
  - define toxicity according to available experimental labels in selected datasets
  - if multiple toxicity/response types exist and effort is reasonable, try to model all
- Target label format:
  - prefer continuous outputs (e.g., IC50/viability-like values) when available
  - convert model outputs to binary labels only when needed for downstream comparison
- Data source guidance:
  - Sanger site was temporarily unavailable at the time
  - scDrugAtlas can be used if all data are obtained successfully; note Harmony processing and avoid direct cross-dataset merging
  - can follow up in person at the lab

## 30.03.2026
- downloaded raw scRNA-seq cell line data from the PERCEPTION paper
  - original publication: https://doi.org/10.1038/s41588-020-00726-6
  - download link: https://singlecell.broadinstitute.org/single_cell/study/SCP542/%20pan-cancer-cell-line-heterogeneity#/
  - 53513 total cells, 22722 genes

-  CTRP^2 downloaded viability scores and annotations for chemical compounds
  - https://depmap.org/portal/data_page/?tab=allData
  - Cell Line ID: v20.meta.per_cell_line (master_ccl_id)
  - Experiments: v20.meta.per_experiment (master_ccl_id is inside)

- PRISM Repurposing for also failed drugs
  - https://depmap.org/portal/download/all/?release=PRISM+Primary+Repurposing+DepMap+Public+24Q2&file=Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv

- GDSC IC50 scores
  - https://cellmodelpassports.sanger.ac.uk/downloads

- clinical phase
- https://repo-hub.broadinstitute.org/repurposing#download-data

## 27.03.2026
### Clearing up project definition
Sent message to Artem to align on project direction before processing data. Waiting on replies for the following core questions:
- Project Focus (Toxicity Definition): Are we targeting Cytotoxicity/Efficacy (predicting if the drug successfully kills heterogeneous cancer cells) OR Adverse Patient Toxicity (predicting dangerous side effects to healthy tissue)?
- Target Variable Labeling: Should the model predict continuous values (e.g., IC50 scores) or binary categorical values (toxic/nontoxic, sensitive/resistant)?
- Dataset Priority: Should scDrugAtlas be the primary focus right now, or is there a better alternative?

## 26.03.2026
### Data Collection
#### Databases

- GDSC2
  - Drug Sensitivity Data -> GDSC2 IC50 Data & GDSC2 Raw Data
  - https://cellmodelpassports.sanger.ac.uk/downloads
  - Contacted DepMap/GDSC team
    - dead documentation link: https://depmap.sanger.ac.uk/documentation/gdsc/
    - asked specifically how the functional dataset “Drug Sensitivity Data” -> “GDSC2 IC50 Data” was processed
  - Status: waiting for response

- ClinTox
  - Binary toxicity prediction based on SMILES
  - https://tdcommons.ai/single_pred_tasks/tox/#clintox
  - maybe not well suited, only binary values

- scDrugAtlas
  - http://drug.hliulab.tech/scDrugAtlas/
  - Contacted Prof. Liu about:
    - original cell line IDs in the consolidated downloadable files
    - source publications for each dataset
    - how the datasets were integrated/combined
    - whether documentation exists for interpreting the consolidated h5ad files
    - whether there is a data dictionary / schema for the main variables and annotations
    - whether additional annotations exist (e.g. IC50 values or sample/cell line IDs) for matching to GDSC
  - Observation:
    - in `breast_cancer_palbociclib`, drug response appears to be encoded as binary values
  - Current issue:
    - hard to identify which cell line is which
    - consolidated file format / variable definitions unclear
  - Status: waiting for response
