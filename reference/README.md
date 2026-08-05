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
