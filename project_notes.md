
# OncoTox Project Notes

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

-
