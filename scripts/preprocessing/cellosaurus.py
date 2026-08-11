"""Resolve cell-line names to Cellosaurus accessions (``CVCL_``).

The pipeline joins SCP542 to CTRPv2 on a **normalised name**, not on an identifier -- see
``ctrp_to_h5ad._normalize_cell_line``. That join is verified, not replaced, by this module: it gives
every cell line a persistent accession so the join can be checked against an external authority and so
the target carries a citable identifier. Nothing here is a join key.

**Source.** Cellosaurus flat file (``cellosaurus.txt``), release recorded in the file's own header and
reported by :func:`load_cellosaurus`. The copy the pipeline reads ships inside Zenodo record
``21807175`` (*Dataset for drevalpy*), so it is pinned by the same accession as the response data --
Cellosaurus itself publishes releases without a stable per-release download URL.

    Bairoch A. The Cellosaurus, a cell-line knowledge resource. *Journal of Biomolecular Techniques*
    29(2):25-38 (2018). https://doi.org/10.7171/jbt.18-2902-002 -- CC BY 4.0.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

#: Line codes we keep. Cellosaurus defines many more; these are the ones any resolution rule needs.
#:   ID = the entry's own name   AC = accession   SY = synonyms
#:   OX = species                CA = category (e.g. "Cancer cell line")
_KEEP = ("ID", "AC", "SY", "OX", "CA")

_SPECIES = re.compile(r"!\s*(?P<name>[^(]+)")


@dataclass(frozen=True)
class CellosaurusRelease:
    """What the flat file says about itself -- the version string a citation needs."""

    version: str
    last_update: str
    path: Path

    def __str__(self) -> str:
        return f"Cellosaurus {self.version} ({self.last_update})"


def load_cellosaurus(path: str | Path) -> tuple[pd.DataFrame, CellosaurusRelease]:
    """Parse the flat file into one row per **name**, plus the release it declares.

    Returns ``(names, release)`` where ``names`` has one row per way an entry can be referred to::

        accession   name                kind          species        category
        CVCL_0035   PC-3                identifier    Homo sapiens   Cancer cell line
        CVCL_0035   PC3                 synonym       Homo sapiens   Cancer cell line

    Long form on purpose: ``kind`` is what lets a caller prefer an entry's own name over some other
    entry's synonym, and ``species``/``category`` are what let it drop candidates that cannot be the
    line in question. Resolution rules belong to the caller; this function only reads the file.
    """
    path = Path(path)
    version = last_update = "unknown"
    rows: list[tuple[str, str, str, str | None, str | None]] = []
    acc = idn = species = category = None
    syns: list[str] = []

    def flush() -> None:
        nonlocal acc, idn, species, category, syns
        if acc and idn:
            rows.append((acc, idn, "identifier", species, category))
            rows.extend((acc, s, "synonym", species, category) for s in syns)
        acc = idn = species = category = None
        syns = []

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            code, _, rest = line.partition("   ")
            rest = rest.strip()
            if version == "unknown" and line.startswith(" Version:"):
                version = line.split(":", 1)[1].strip()
            elif last_update == "unknown" and line.startswith(" Last update:"):
                last_update = line.split(":", 1)[1].strip()
            if code not in _KEEP:
                continue
            if code == "ID":
                flush()
                idn = rest
            elif code == "AC":
                acc = rest
            elif code == "SY":
                syns = [s.strip() for s in rest.split(";") if s.strip()]
            elif code == "OX":
                m = _SPECIES.search(rest)
                species = m.group("name").strip() if m else rest
            elif code == "CA":
                category = rest
        flush()

    names = pd.DataFrame(rows, columns=["accession", "name", "kind", "species", "category"])
    return names, CellosaurusRelease(version=version, last_update=last_update, path=path)


def lookup_key(s: str) -> str:
    """Normalise a name for *looking it up in Cellosaurus* -- never for joining two datasets.

    Deliberately looser than ``ctrp_to_h5ad._normalize_cell_line``: Cellosaurus writes names with
    spaces and dots (``PC.3``, ``C32 [Human melanoma]``) that neither SCP542 nor CTRP produce, so the
    production key never has to handle them. Keeping the two functions separate means widening this
    one cannot silently widen the join.
    """
    return str(s).strip().lower().replace("-", "").replace(" ", "").replace(".", "")


#: Tie-break rules, applied **in this order and only while more than one candidate remains**.
#: Each drops candidates on a property Cellosaurus records, so no rule expresses a view about which
#: line was *meant*. They are tie-breakers, never filters: a name matching exactly one entry keeps
#: that entry even if the rule would have rejected it, so a wrong-species match surfaces as a
#: mismatch instead of vanishing.
TIE_BREAKS: tuple[tuple[str, str, str], ...] = (
    ("species", "Homo sapiens", "human only"),
    ("category", "Cancer cell line", "cancer cell lines only"),
    ("kind", "identifier", "a primary identifier beats a synonym"),
)


def resolve_accessions(names_to_resolve, cellosaurus_names: pd.DataFrame) -> pd.DataFrame:
    """Resolve each name to at most one accession, recording which rule did it.

    One row per input name::

        name    key     accession   kind        n_candidates  resolved_by   status
        PC3     pc3     CVCL_0035   identifier  2             human only    resolved
        FOO     foo     <NA>        <NA>        0             <NA>          absent

    ``status`` is ``resolved``, ``absent`` (no candidate) or ``ambiguous`` (candidates survived every
    rule). **Ambiguous is reported, never guessed** -- picking a winner there would be a judgement
    about cell-line identity, which is not a lookup and does not belong in a preprocessing step.
    """
    cand = cellosaurus_names.assign(key=lambda d: d.name.map(lookup_key))
    by_key = {k: g for k, g in cand.groupby("key", sort=False)}

    out = []
    for nm in pd.unique(pd.Series(list(names_to_resolve), dtype=object)):
        key = lookup_key(nm)
        g = by_key.get(key)
        if g is None or g.empty:
            out.append((nm, key, pd.NA, pd.NA, 0, pd.NA, "absent"))
            continue

        n_cand = g.accession.nunique()
        resolved_by = pd.NA
        for col, wanted, label in TIE_BREAKS:
            if g.accession.nunique() <= 1:
                break
            narrowed = g[g[col] == wanted]
            if not narrowed.empty and narrowed.accession.nunique() < g.accession.nunique():
                g, resolved_by = narrowed, label

        if g.accession.nunique() == 1:
            row = g.sort_values("kind").iloc[0]  # 'identifier' sorts before 'synonym'
            out.append((nm, key, row.accession, row.kind, n_cand, resolved_by, "resolved"))
        else:
            out.append((nm, key, pd.NA, pd.NA, n_cand, resolved_by, "ambiguous"))

    return pd.DataFrame(out, columns=["name", "key", "accession", "kind",
                                      "n_candidates", "resolved_by", "status"])
