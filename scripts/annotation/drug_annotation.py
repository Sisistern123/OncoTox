"""Join CTRPv2's compound annotation to the response data on ``master_cpd_id``, never on names.

The response values come from DrEval's reprocessed ``CTRPv2.csv``; the annotation a drug panel is
selected with -- approval status, protein target, mechanism -- comes from
``data/drug/all_sources_drug_catalog.csv``, built from CTRP's own ``v20.meta.per_compound.txt``
(``notebooks/data_and_harmonization/drug_catalog.ipynb``). The two name the same 545 compounds
differently, so anything that joins them by name loses compounds silently.

**How badly.** 102 of 545 do not match on lower-cased names, because DrEval renamed to preferred
names (``abt-199`` -> ``Venetoclax``, ``byl-719`` -> ``Alpelisib``, ``sirolimus`` -> ``Rapamycin``),
changed separators (``:`` -> ``-``) and altered hyphenation (``azd7545`` -> ``azd-7545``). **15 of
the unmatched are single-agent, FDA-approved or clinical** -- exactly the compounds a
clinically-motivated panel selects, and two of them (``sirolimus``, ``gdc-0941``) were already lost
once by the voided literature panel.

**The key that works.** ``CTRPv2.csv`` carries ``master_cpd_id``, CTRP's own compound identifier,
which is also the catalog's ``identifier`` (``source_id_type == "master_cpd_id"``): **545 of 545 in
both directions, nothing unmatched.** This is the opposite conclusion to the cell-line join, where
names beat Cellosaurus accessions 180 to 172 (``scripts/preprocessing/ctrp_to_h5ad.py``) -- there the
identifier was incomplete, here it is exact, so each join is decided on its own evidence.

Nothing here is written to disk. The annotation is derived on demand from the two files that already
exist, so there is no third copy of it to drift out of date.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from scripts.sources.pubchem import with_parents

#: A CTRPv2 combination is named for its two agents and their molar ratio -- ``"alisertib-navitoclax
#: (2-1 mol-mol)"`` in DrEval's spelling, ``"alisertib:navitoclax (2:1 mol/mol)"`` in CTRP's. The two
#: patterns agree on the same **49** compounds, which is what makes either safe to use as the flag.
_COMBINATION = re.compile(r"\(\d+(?:\.\d+)?[-:]\d+(?:\.\d+)?\s*mol[-/]mol\)")

#: Sourced corrections to the ``pubchem_id`` DrEval publishes, keyed by ``master_cpd_id``. Applied in
#: :func:`ctrp_compounds`, explicitly and one at a time -- never a pattern -- so that every correction
#: names the evidence that justifies it, in the manner of ``CTRP_CELL_LINE_ALIASES`` in
#: :mod:`scripts.preprocessing.ctrp_to_h5ad`.
CTRP_PUBCHEM_OVERRIDES: dict[int, int] = {
    # 375395, CTRP compound name "Platin". DrEval records CID 23939, which is *elemental platinum* --
    # the element, not a drug. CTRP's own `v20.meta.per_compound.txt` gives this compound the SMILES
    # `N[Pt](N)(Cl)Cl`, i.e. cis-diamminedichloroplatinum(II) = cisplatin, and the vendor catalogue
    # entry Selleck S1166, which is cisplatin. Corrected to CID 5702198.
    #
    # This matters beyond tidiness: with a nonstandard name *and* a wrong structure, cisplatin is
    # unreachable by every key an external drug list can use, so the most widely used platinum agent
    # silently fails to appear in any literature-based selection.
    375395: 5702198,
}


def ctrp_compounds(response_csv: str | Path) -> pd.DataFrame:
    """The compound axis of DrEval's ``CTRPv2.csv`` -- one row per compound, indexed by identifier.

    Columns: ``drug`` (DrEval's name, lower-cased -- the spelling ``Y_ctrp``'s columns carry, see
    ``ctrp_to_h5ad._normalize_drug``), ``drug_dreval`` (as published), ``pubchem_cid`` and
    ``is_combination``.

    ``pubchem_cid`` is ``NA`` for 56 compounds: the 49 combinations, which have no single structure,
    and 7 single-agent BRD probes with no PubChem record. DrEval fills the field with the compound
    name in those cases rather than leaving it empty, so it is coerced to a number here and the
    non-numeric entries become missing -- which is what they mean.
    """
    df = pd.read_csv(
        response_csv,
        usecols=["master_cpd_id", "drug_name", "pubchem_id"],
        dtype={"pubchem_id": str},
    ).drop_duplicates("master_cpd_id")

    out = pd.DataFrame({
        "drug": df.drug_name.str.strip().str.lower(),
        "drug_dreval": df.drug_name,
        "pubchem_cid": pd.to_numeric(df.pubchem_id, errors="coerce").astype("Int64"),
        "is_combination": df.drug_name.str.contains(_COMBINATION),
    })
    out.index = pd.Index(df.master_cpd_id.astype(int), name="master_cpd_id")

    missing = set(CTRP_PUBCHEM_OVERRIDES) - set(out.index)
    if missing:
        raise ValueError(f"pubchem_id overrides for compounds not in the response file: {missing}")
    for cpd_id, cid in CTRP_PUBCHEM_OVERRIDES.items():
        out.loc[cpd_id, "pubchem_cid"] = cid

    return out.sort_index()


#: Catalog columns carried through, and what each is for. ``compound_status`` is CTRP's own
#: ``cpd_status`` (``FDA`` / ``clinical`` / ``probe`` / ``GE-active``); ``top_test_conc_umol`` is the
#: highest concentration that compound was screened at, which spans 0.13-600 uM across the 545 and so
#: bounds what any response measure can mean for it.
_CATALOG_COLUMNS = {
    "compound_name_norm": "drug_ctrp",
    "compound_status": "compound_status",
    "target": "target",
    "moa_or_pathway": "moa_or_pathway",
    "top_test_conc_umol": "top_test_conc_umol",
}


def annotate_compounds(compounds: pd.DataFrame, catalog_csv: str | Path) -> pd.DataFrame:
    """Attach CTRP's compound annotation to :func:`ctrp_compounds`, joined on ``master_cpd_id``.

    Adds the columns in :data:`_CATALOG_COLUMNS`. ``drug_ctrp`` is CTRP's own spelling and is kept
    alongside DrEval's ``drug`` because **neither name space is sufficient on its own** when matching
    an external drug list: against Sun et al. 2017, DrEval's spelling is the only one that finds
    ``idelalisib`` and CTRP's is the only one that finds ``fluorouracil`` and ``mitomycin``.

    :raises ValueError: if the join is not exact in both directions. It is exact today (545/545), and
        an inexact one means the catalog and the response file have drifted apart -- which is the
        failure this module exists to prevent, so it must stop the run rather than drop rows.
    """
    catalog = pd.read_csv(catalog_csv)
    catalog = catalog[catalog.dataset == "CTRPv2"]
    if not (catalog.source_id_type == "master_cpd_id").all():
        raise ValueError("Catalog CTRPv2 rows are not all keyed by master_cpd_id -- cannot join.")
    catalog = catalog.set_index(catalog.identifier.astype(int).rename("master_cpd_id"))

    only_response = compounds.index.difference(catalog.index)
    only_catalog = catalog.index.difference(compounds.index)
    if len(only_response) or len(only_catalog):
        raise ValueError(
            f"master_cpd_id join is not exact: {len(only_response)} compounds in the response file "
            f"are absent from the catalog and {len(only_catalog)} the reverse. "
            f"Response-only: {sorted(only_response)[:10]}; catalog-only: {sorted(only_catalog)[:10]}."
        )

    annotated = compounds.join(catalog[list(_CATALOG_COLUMNS)].rename(columns=_CATALOG_COLUMNS))
    annotated["drug_ctrp"] = annotated.drug_ctrp.str.strip().str.lower()
    return annotated


def match_external_list(
    annotated: pd.DataFrame,
    external: pd.DataFrame,
    *,
    name_col: str,
    cids_col: str | None = None,
    parent_cache: str | Path | None = None,
) -> pd.DataFrame:
    """Match an external drug list onto CTRPv2, returning it with ``master_cpd_id`` and ``matched_by``.

    Every row of ``external`` is returned, matched or not -- a selection criterion has to be able to
    show what it *failed* to find as well as what it found.

    Four keys are tried, and ``matched_by`` records every one that succeeded, because they are not
    equally strong and the panel should not hide which compound rests on which:

    * ``name_dreval`` -- the external name equals DrEval's spelling;
    * ``name_ctrp`` -- it equals CTRP's own spelling. Both are needed: against Sun et al. 2017,
      DrEval's finds ``idelalisib`` and CTRP's finds ``fluorouracil`` and ``mitomycin``.
    * ``pubchem_cid`` -- the structures are identical though no name matches, which is how
      ``Vemurafenib`` reaches the ``plx-4032`` CTRPv2 screened it as.
    * ``pubchem_parent`` -- the structures agree only after both sides are resolved to PubChem's
      parent compound, which is how ``Imatinib mesylate`` reaches ``imatinib`` and, in the other
      direction, ``Cytarabine`` reaches ``cytarabine hydrochloride``. Requires ``parent_cache``.
      This is the weakest of the four and the one to inspect first if a match looks wrong.

    ``cids_col`` names a semicolon-joined CID column on ``external``; ``parent_cache`` is the CSV
    :func:`scripts.sources.pubchem.with_parents` caches lookups in.

    :raises ValueError: if a row matches two different compounds, which would make the panel depend on
        the order the keys were tried rather than on the data.
    """
    def _cids(value: object) -> list[int]:
        # Read back from CSV, a column of mostly-single CIDs infers inconsistent dtypes across
        # pandas' read chunks, so a CID can arrive as "4033.0". Accept that, but only when it is
        # exactly an integer -- a fractional value would mean the column is not CIDs at all.
        out = []
        for part in str(value).split(";"):
            if not part.strip() or part.strip().lower() == "nan":
                continue
            number = float(part)
            if not number.is_integer():
                raise ValueError(f"{part!r} is not a PubChem CID.")
            out.append(int(number))
        return out

    by_dreval = {n: i for i, n in annotated.drug.items()}
    by_ctrp = {n: i for i, n in annotated.drug_ctrp.items()}
    ctrp_cids = {int(c): i for i, c in annotated.pubchem_cid.dropna().items()}

    external_cids = {c for _, row in external.iterrows() for c in _cids(row[cids_col])} if cids_col \
        else set()
    if parent_cache is not None:
        # Both sides are expanded, not just the external list: CTRPv2 screened `cytarabine
        # hydrochloride` while the FDA list names `Cytarabine`, so a one-sided expansion still misses it.
        expanded = with_parents(sorted(set(ctrp_cids) | external_cids), parent_cache)
        by_parent: dict[int, set[int]] = {}
        for cid, index in ctrp_cids.items():
            for equivalent in expanded[cid]:
                by_parent.setdefault(equivalent, set()).add(index)
    else:
        expanded, by_parent = {}, {}

    def _lookup(cids: list[int]) -> int | None:
        """First compound reachable from ``cids`` through the parent relation, if unambiguous."""
        for cid in cids:
            for equivalent in expanded.get(cid, [cid]):
                hits = by_parent.get(equivalent, set())
                if len(hits) > 1:
                    raise ValueError(
                        f"PubChem CID {equivalent} reaches {len(hits)} CTRPv2 compounds "
                        f"({sorted(hits)}); the parent relation is ambiguous here."
                    )
                if hits:
                    return next(iter(hits))
        return None

    ids, how = [], []
    for _, row in external.iterrows():
        name = str(row[name_col]).strip().lower()
        cids = _cids(row[cids_col]) if cids_col else []
        exact = next((ctrp_cids[c] for c in cids if c in ctrp_cids), None)
        hits = {
            "name_dreval": by_dreval.get(name),
            "name_ctrp": by_ctrp.get(name),
            "pubchem_cid": exact,
            # Reported only when the parent relation is what found it, so the column stays readable.
            "pubchem_parent": None if exact is not None else _lookup(cids),
        }
        found = {k: v for k, v in hits.items() if v is not None}
        if len(set(found.values())) > 1:
            raise ValueError(f"{row[name_col]!r} matches more than one CTRPv2 compound: {found}")
        ids.append(next(iter(found.values()), pd.NA))
        how.append("+".join(found) if found else "")

    return external.assign(
        master_cpd_id=pd.array(ids, dtype="Int64"),
        matched_by=how,
    )
