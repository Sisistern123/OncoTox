# `reference/` — external mapping tables, versioned with the code

Provenance-critical lookups that a result depends on. They live here rather than under the gitignored
data root because a mapping that cannot be reproduced makes every number derived from it unverifiable.

## `hgnc_complete_set.txt`

The HGNC approved-symbol set, used to test whether genes that `gen_embeds.py` discarded as
out-of-vocabulary are in fact present in scGPT's vocabulary under a **current** symbol. See
[the write-up](../docs/steps/corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols)
and `notebooks/data_and_harmonization/gene_symbol_rescue.ipynb`.

### Exact source

| | |
|---|---|
| Publisher | HUGO Gene Nomenclature Committee (HGNC), EMBL-EBI — genenames.org |
| Dataset | `hgnc_complete_set`, the complete set of HGNC-approved gene symbols, TSV |
| URL retrieved from | `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt` |
| Retrieved | **05.08.2026**, 16:45 local (`curl`, no redirect) |
| Server `last-modified` | **04.08.2026 13:01:56 GMT** — the release snapshot this file is |
| SHA-256 | `2106d1f237d6c542a85a4c399225e011ee9c5822199e7a400b7b9f842c8c8ca0` |
| Size | 16,930,015 bytes |
| Rows | 45,031, all `status == Approved` |
| Rows with a `prev_symbol` | 12,700 |

### Which release this is, and why it is identified this way

HGNC does not stamp a version into the file. The release is therefore pinned by three independent
markers, all derived from the artifact itself:

- the server's `last-modified`, **04.08.2026** — when this snapshot was published;
- the latest `date_modified` across all rows, **03.08.2026** — the newest record it contains;
- the latest `date_symbol_changed`, **15.07.2026** — the most recent rename it knows about.

**There is no permanent URL for this release.** HGNC publishes monthly and **overwrites the URL above in
place**, so re-downloading it later yields a *different* file under the same address. Three dated
archive paths were probed on 05.08.2026 and all returned HTTP 404:

```
https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-08-01.txt
https://storage.googleapis.com/public-download-files/hgnc/archive/monthly/tsv/hgnc_complete_set_2026-08-01.txt
https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2026-08-01.txt
```

**So the committed file is the citation.** Any result derived from it must be reproduced against *this*
copy, verified by the checksum above (`shasum -a 256 reference/hgnc_complete_set.txt`) — not against
whatever the URL serves on the day someone re-runs the notebook. That is the whole reason a 17 MB file
is versioned with the code rather than fetched on demand.

**Columns this project uses.** `symbol` (current approved), `prev_symbol` (pipe-separated former
approved symbols — a rename, authoritative), `alias_symbol` (pipe-separated synonyms — looser, may map
two genes onto one name), `status`, `ensembl_gene_id`.

**Citation.** Seal, R. L. *et al.* Genenames.org: the HGNC resources in 2023. *Nucleic Acids Research*
**51**, D1003–D1009 (2023). doi:10.1093/nar/gkac888

## `sun2017_fda_anticancer_drugs.csv`

The FDA-approved anticancer drug list the [drug panel](../docs/steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026)
is selected on. Retrieved and parsed by `scripts/preprocessing/fetch_sun2017_drugs.py`; used in
`notebooks/drug_selection/literature_panel.ipynb` §1.

### Exact source

| | |
|---|---|
| Publication | Sun, J., Wei, Q., Zhou, Y., Wang, J., Liu, Q. & Xu, H. A systematic analysis of FDA-approved anticancer drugs. *BMC Systems Biology* **11**(Suppl 5), 87 (2017) |
| DOI | `10.1186/s12918-017-0464-7` — PMC5629554 |
| Licence | **CC BY 4.0**, which is why the parsed table may be committed here |
| What was taken | **Table 1**, *"Summary of FDA-approved anticancer drugs from 1949 to 2014"* |
| Retrieved | **11.08.2026**, NCBI E-utilities `efetch.fcgi?db=pmc&id=5629554&rettype=xml` |
| Rows | 150 drugs — 61 cytotoxic, 89 targeted |

**Table 1 is the entire dataset.** The paper ships no supplementary file: *"All data generated or
analysed during this study are included in this published article."* It is fetched as JATS XML rather
than scraped from the article HTML, because the XML marks table structure explicitly.

**How the version is pinned.** PMC publishes no checksum for a rendered article, so there is nothing to
verify bytes against. Instead the parse must reproduce the three counts the paper states in its own
Results — 150 drugs, 61 cytotoxic, 89 targeted — and **raises** otherwise. If NCBI changes its markup or
the article record is revised, that fails loudly rather than handing a different table to the panel.

**Columns.** `drug`, `approval_year`, `therapeutic_class`, `target_gene`, `delivery_type`, `drug_class`
(cytotoxic or targeted — carried in the source as two section rows interleaved with the data, and the
only machine-readable place that distinction exists), and `pubchem_cids`.

⚠️ **`pubchem_cids` is not from the paper.** Sun et al. give names only; the CIDs were resolved
separately so the list can be matched to a screen that uses development codes. The provenance file
records that separately from the article.

⚠️ **The list stops at 2014.** Nothing approved since is in it. The concrete cost, recorded so it is not
rediscovered: `selumetinib` (approved 2020) cannot be a panel candidate under this criterion.

## `pubchem_parent_cids.csv`

PubChem compound identifiers mapped to their **parent** compound — the neutral form carrying the active
moiety — for every CID either the drug list or CTRPv2 uses. 566 rows, written by
`scripts/preprocessing/pubchem.py::parent_cids`.

**Why it exists.** An approval list names `Imatinib mesylate` because that is what the FDA approved; a
screen names `imatinib` because that is what it dissolved. These are different molecules with different
PubChem records, so neither the name nor the compound identifier relates them — only the parent relation
does. Applied to **both** sides it recovered 13 panel candidates that name and structure matching
missed, including imatinib, doxorubicin, vincristine and topotecan.

| | |
|---|---|
| Source | PubChem PUG-REST, `compound/cid/{cid}/cids/JSON?cids_type=parent` |
| Retrieved | **11.08.2026** |
| Citation | Kim, S. *et al.* PubChem 2023 update. *Nucleic Acids Research* **51**, D1373–D1380 (2023). doi:10.1093/nar/gkac956 |

**Committed rather than fetched on demand**, for the same reason as the HGNC file above: PubChem is
curated continuously, so a panel that depends on a live lookup is not reproducible. `parent_cids` fetches
only CIDs the cache does not already hold, so re-running costs nothing and a machine with no network
reproduces the same 11 drugs.
