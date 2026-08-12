"""PubChem lookups that let two drug lists be compared by structure rather than by name.

Drug lists name the same molecule differently, and no amount of string normalization relates the
spellings. Three kinds of difference occur here, and only the third has a chemical answer:

* **development code vs INN** -- ``plx-4032`` is ``Vemurafenib``, ``abt-199`` is ``Venetoclax``;
* **preferred-name drift** between two curations of the same source, handled for CTRPv2 by joining on
  ``master_cpd_id`` (:mod:`scripts.preprocessing.drug_annotation`);
* **salt form vs free base** -- a regulatory list names ``Imatinib mesylate`` because that is what the
  FDA approved, a screen names ``imatinib`` because that is what it dissolved. These are different
  molecules with different PubChem records, so matching on the compound identifier alone still misses
  them.

PubChem answers the third by publishing, for a salt, the CID of its **parent** -- the neutral form
carrying the active moiety. Comparing each side's CIDs *together with their parents* matched 12 more
FDA-approved drugs to CTRPv2 than names and raw CIDs did, ``Imatinib mesylate``, ``Doxorubicin
hydrochloride``, ``Vincristine sulfate`` and ``Topotecan hydrochloride`` among them, and lost none.

**Retrieved, not derived.** Both lookups are network calls against a database that is curated
continuously, so a result is only reproducible if the answer is stored. Parent lookups are therefore
cached in a committed CSV: re-running fetches nothing it already has, and a machine with no network
reproduces the same panel from the cache.

    PubChem PUG-REST, https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest -- Kim, S. *et al.* PubChem 2023
    update. *Nucleic Acids Research* **51**, D1373-D1380 (2023). doi:10.1093/nar/gkac956
"""

from __future__ import annotations

import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

#: PubChem asks for no more than 5 requests/second.
_DELAY_S = 0.22
_NAME_TO_CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
_CID_TO_PARENT = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/cids/JSON?cids_type=parent"


def _get_cids(url: str) -> list[int]:
    """One PUG-REST call returning its CID list; an unknown identifier answers 404, not an empty list."""
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return list(resp.json()["IdentifierList"]["CID"])
    except (requests.RequestException, KeyError, ValueError):
        return []


def resolve_names_to_cids(names: pd.Series) -> pd.Series:
    """Resolve drug names to semicolon-joined CIDs; empty where the name has no PubChem record.

    All CIDs a name returns are kept rather than only the first, because a query can answer with a
    free base and its salts and either may be the form that was screened.

    An empty result is informative rather than a failure: for Sun et al. 2017 the 30 names that
    resolve to nothing are all biologics -- antibodies, enzymes, a cell therapy, a radiopharmaceutical
    -- which have no small-molecule record to find and equally cannot appear in a small-molecule
    screen.
    """
    out = []
    for name in names:
        out.append(";".join(str(c) for c in _get_cids(_NAME_TO_CID.format(name=quote(str(name))))))
        time.sleep(_DELAY_S)
    return pd.Series(out, index=names.index, name="pubchem_cids")


def parent_cids(cids: list[int], cache_csv: str | Path) -> dict[int, list[int]]:
    """Parent CID(s) of each input CID, fetching only what ``cache_csv`` does not already hold.

    The cache is a committed two-column CSV (``cid``, ``parent_cids``) shared by every list that needs
    resolving, so the ~600 lookups this project makes happen once. A CID that is already a parent
    answers with itself; the entry is stored either way, so a cached empty is distinguishable from a
    CID that was never asked about.
    """
    cache_path = Path(cache_csv)
    if cache_path.exists():
        cached = pd.read_csv(cache_path, dtype={"parent_cids": str}).fillna({"parent_cids": ""})
        known = {int(r.cid): r.parent_cids for r in cached.itertuples()}
    else:
        known = {}

    missing = [c for c in dict.fromkeys(cids) if c not in known]
    if missing:
        print(f"  PubChem parents: {len(missing)} new CIDs to resolve "
              f"({len(known)} already cached in {cache_path.name})")
        for cid in missing:
            known[cid] = ";".join(str(p) for p in _get_cids(_CID_TO_PARENT.format(cid=cid)))
            time.sleep(_DELAY_S)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        (pd.DataFrame({"cid": list(known), "parent_cids": list(known.values())})
           .sort_values("cid").to_csv(cache_path, index=False))

    return {c: [int(p) for p in known[c].split(";") if p] for c in dict.fromkeys(cids)}


def with_parents(cids: list[int], cache_csv: str | Path) -> dict[int, list[int]]:
    """Each CID mapped to itself **and** its parents -- the set two lists should be compared on.

    Applied to both sides of a comparison, not one: CTRPv2 screened ``cytarabine hydrochloride`` while
    the FDA list names ``Cytarabine``, so resolving only the external list would still miss it.
    """
    parents = parent_cids(cids, cache_csv)
    return {c: list(dict.fromkeys([c, *parents[c]])) for c in parents}
