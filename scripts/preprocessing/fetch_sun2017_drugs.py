"""Fetch the FDA-approved anticancer drug list the panel is selected from.

**What this downloads.** Table 1 of a systematic review that enumerates every drug the FDA approved
with a cancer indication between 1949 and 2014, with its approval year, therapeutic class, target
gene and whether it is given single-agent or in combination:

    Sun, J., Wei, Q., Zhou, Y., Wang, J., Liu, Q. & Xu, H. A systematic analysis of FDA-approved
    anticancer drugs. *BMC Systems Biology* **11**(Suppl 5), 87 (2017).
    doi:10.1186/s12918-017-0464-7 -- PMC5629554, CC BY 4.0.

The authors compiled the list from NCI drug information, the MediLexicon cancer drug list and
NavigatingCancer, then verified each entry against Drugs@FDA, DailyMed and DrugBank.

**Why this list defines the drug panel.** The panel has to be justified by something external to our
own labels; selecting on response statistics is what voided the previous two panels
(``docs/steps/corrections-and-dead-ends.md``). "Approved by the FDA for a cancer indication" is
external, dated and reproducible by someone who never sees our AUCs, and this paper is the citable
form of it with mechanism already attached -- which a bare approval list is not.

**Its one limitation, stated because it bounds the panel.** The cut-off is end-2014, so nothing
approved since is in it: no checkpoint inhibitor beyond ``Pembrolizumab``/``Nivolumab``, no CDK4/6
inhibitor, no third-generation EGFR inhibitor. For an intersection with CTRPv2 -- screened 2012-2015 --
that costs little, but the list is not a current statement of what is approved.

**Table 1 is the dataset.** The paper carries no supplementary file: *"All data generated or analysed
during this study are included in this published article."* It is fetched from NCBI E-utilities as
JATS XML rather than scraped from the article HTML, because the XML marks up table structure
explicitly and is the archival form.

**How the version is pinned.** PMC publishes no checksum for a rendered article, so the retrieval is
verified against the counts the paper itself states in its Results -- 150 drugs, 61 cytotoxic and 89
targeted. If NCBI's markup or the article record ever changes shape, that assertion fails rather than
a silently different table reaching the panel.
"""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import pandas as pd
import requests

from scripts.preprocessing.pubchem import resolve_names_to_cids

#: PubMed Central identifier of Sun et al. 2017. The article is the citation; this is how it is
#: retrieved.
PMC_ID = "PMC5629554"
DOI = "10.1186/s12918-017-0464-7"
#: E-utilities returns JATS XML for ``db=pmc``, which marks up ``<table-wrap>`` structure explicitly.
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

#: What the paper's Results state, and what the parse is checked against. Not a checksum -- a
#: statement made in the text that the extracted table has to reproduce.
EXPECTED_TOTAL = 150
EXPECTED_BY_CLASS = {"Cytotoxic": 61, "Targeted": 89}

#: Table 1's own column labels, in order, mapped to the names used downstream.
_COLUMNS = {
    "Drug": "drug",
    "Approval year": "approval_year",
    "Therapeutic class": "therapeutic_class",
    "Target gene": "target_gene",
    "Delivery type": "delivery_type",
}


def parse_table1(xml_text: str) -> pd.DataFrame:
    """Extract Table 1 from the article's JATS XML as one row per drug.

    The table interleaves data rows with single-cell rows that label the two sections -- ``Cytotoxic``
    and ``Targeted``. Those are not drugs; they carry the classification for every row beneath them,
    so they are consumed into a ``drug_class`` column rather than dropped, which is the only place
    that distinction exists in the source.

    :raises ValueError: if Table 1 is absent, its header is not the expected five columns, or the
        result does not reproduce the counts the paper states (:data:`EXPECTED_TOTAL`,
        :data:`EXPECTED_BY_CLASS`).
    """
    root = ET.fromstring(xml_text)
    wraps = [w for w in root.iter("table-wrap") if (w.findtext("label") or "").strip() == "Table 1"]
    if not wraps:
        raise ValueError(f"{PMC_ID}: no <table-wrap> labelled 'Table 1' in the fetched XML.")

    rows = [[" ".join("".join(cell.itertext()).split()) for cell in tr] for tr in wraps[0].iter("tr")]
    header, body = rows[0], rows[1:]
    if header != list(_COLUMNS):
        raise ValueError(f"Table 1 header changed: expected {list(_COLUMNS)}, got {header}.")

    drug_class, records = None, []
    for row in body:
        if len(row) == 1:                      # a section label, not a drug
            drug_class = row[0]
            continue
        if len(row) != len(_COLUMNS):
            raise ValueError(f"Table 1 row has {len(row)} cells, expected {len(_COLUMNS)}: {row}")
        records.append({**dict(zip(_COLUMNS.values(), row)), "drug_class": drug_class})

    table = pd.DataFrame(records)
    counts = table.drug_class.value_counts().to_dict()
    if len(table) != EXPECTED_TOTAL or counts != EXPECTED_BY_CLASS:
        raise ValueError(
            f"Parsed {len(table)} drugs {counts}, but the paper states "
            f"{EXPECTED_TOTAL} {EXPECTED_BY_CLASS}. Do not select a panel from this table."
        )
    return table


def fetch_sun2017_drugs(reference_dir: str | Path, *, force: bool = False) -> Path:
    """Retrieve Table 1 and write it to ``reference_dir`` as CSV, with its provenance beside it.

    Idempotent: an existing CSV is re-read and re-verified rather than re-fetched, so this is safe to
    call from a notebook on every run. The CSV is committed to the repository -- it is 150 rows, and a
    selection criterion that depends on a live network call is not reproducible on a machine that
    cannot reach NCBI.

    Returns the path of the CSV. Writes ``sun2017_fda_anticancer_drugs.provenance.json`` next to it
    recording the article, its DOI, the E-utilities request and the retrieval date, so the version is
    readable from the data rather than only from this file.
    """
    dest = Path(reference_dir)
    dest.mkdir(parents=True, exist_ok=True)
    csv_path = dest / "sun2017_fda_anticancer_drugs.csv"

    if csv_path.exists() and not force:
        table = pd.read_csv(csv_path)
        counts = table.drug_class.value_counts().to_dict()
        if len(table) == EXPECTED_TOTAL and counts == EXPECTED_BY_CLASS:
            print(f"  {csv_path.name}: cached, {len(table)} drugs {counts} -- skipping fetch")
            return csv_path
        print(f"  {csv_path.name}: cached copy fails verification ({len(table)} {counts}) -- refetching")

    params = {"db": "pmc", "id": PMC_ID.removeprefix("PMC"), "rettype": "xml"}
    print(f"Sun et al. 2017 Table 1 -- {PMC_ID}, doi:{DOI}")
    resp = requests.get(_EFETCH, params=params, timeout=120)
    resp.raise_for_status()

    table = parse_table1(resp.text)
    print(f"  parsed {len(table)} drugs {table.drug_class.value_counts().to_dict()}")

    table["pubchem_cids"] = resolve_names_to_cids(table.drug)
    n_resolved = int((table.pubchem_cids != "").sum())
    print(f"  PubChem: {n_resolved} of {len(table)} names resolved to a CID "
          f"({len(table) - n_resolved} biologics with no small-molecule record)")
    table.to_csv(csv_path, index=False)
    print(f"  written to {csv_path}")

    csv_path.with_suffix(".provenance.json").write_text(json.dumps({
        "citation": ("Sun J, Wei Q, Zhou Y, Wang J, Liu Q, Xu H. A systematic analysis of "
                     "FDA-approved anticancer drugs. BMC Systems Biology 11(Suppl 5), 87 (2017)."),
        "doi": DOI,
        "pmc_id": PMC_ID,
        "licence": "CC BY 4.0",
        "source": "Table 1, 'Summary of FDA-approved anticancer drugs from 1949 to 2014'",
        "retrieved_from": f"{_EFETCH}?db=pmc&id={params['id']}&rettype=xml",
        "retrieved": date.today().isoformat(),
        "pubchem_cids": {
            "note": ("The pubchem_cids column is NOT from the paper. Sun et al. list names only; "
                     "CIDs were resolved separately so the list can be matched to CTRPv2, which "
                     "screened many of these compounds under development codes."),
            "source": "PubChem PUG-REST, /compound/name/{name}/cids/JSON",
            "resolved": date.today().isoformat(),
            "n_resolved": n_resolved,
            "unresolved_are": "biologics (antibodies, enzymes, cell therapy, radiopharmaceutical)",
        },
        "verified_against": {"total": EXPECTED_TOTAL, "by_class": EXPECTED_BY_CLASS},
        "coverage_cutoff": "1949-2014; nothing approved after 2014 is listed",
        "fetched_by": "scripts/preprocessing/fetch_sun2017_drugs.py",
    }, indent=2) + "\n")
    return csv_path
