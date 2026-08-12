"""Map SCP542's gene symbols onto the current HGNC approved set.

SCP542 was annotated against an older gene nomenclature than scGPT's vocabulary, so a gene
renamed between the two releases fails ``gen_embeds.py``'s exact string match and is thrown
away although the vocabulary holds it under its current name -- 775 genes carrying 3.6% of
every cell in the ``all_genes`` variant. The defect and its measurement are written up in
``docs/steps/corrections-and-dead-ends.md``; the counts come from
``notebooks/outputs/embeddings/gene_symbol_rescue.csv``.

This module only *annotates*. It writes ``var['hgnc_symbol']`` and leaves ``var_names``,
``.X`` and every expression value untouched, so the HVG set and ``X_pca`` stay bit-identical
and any change in a scGPT result is attributable to the recovered genes alone. Deciding which
symbol scGPT is actually given is left to ``gen_embeds.py``, which is the only consumer that
knows the vocabulary.

Source table: ``reference/hgnc_complete_set.txt``, pinned by checksum in
``reference/README.md`` (HGNC overwrites the file in place at its URL, so the committed copy
is the citation). Seal et al., *Genenames.org: the HGNC resources in 2023*, Nucleic Acids
Research 51, D1003-D1009 (2023), doi:10.1093/nar/gkac888.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HGNC_FILE = Path(__file__).resolve().parents[2] / "reference" / "hgnc_complete_set.txt"

# alias_symbol is deliberately not used: an alias is a synonym rather than an official
# rename, and one alias can name two distinct genes. Decided 05.08.2026 (Selin) together
# with dropping the GENCODE/Ensembl-ID route -- the conservative rescue over the larger one.
RENAME_FIELD = "prev_symbol"
CURRENT_FIELD = "symbol"


def load_rename_map(hgnc_file: Path | str = HGNC_FILE) -> tuple[dict[str, str], set[str]]:
    """Return ``{former symbol -> current symbol}`` and the set of symbols refused as unsafe.

    HGNC stores former symbols pipe-separated in one field, so it is exploded to one row per
    ``(former -> current)`` pair. Two classes of pair are then refused, because ``prev_symbol``
    records **reassignments** as well as renames and a wrong lookup routes one gene's counts to
    another gene's token, which is worse than the OOV drop it repairs:

    1. **Former symbols HGNC attributes to more than one current gene** (86 of them). There is
       no way to pick a winner without a judgement about gene identity.
    2. **Former symbols that are themselves a current approved symbol of a different gene.**
       ``OSR1`` is the approved symbol of odd-skipped related 1 *and* a former symbol of
       ``OXSR1``; likewise ``NTNG1`` -> ``NTNG2``, ``ADCY3`` -> ``ADCY8``, ``SRGAP2`` ->
       ``SRGAP3``. Where SCP542 uses a symbol that is approved today, it is read as meaning the
       gene that holds it today. *Added 05.08.2026 (Selin) after the guard was found missing
       from the rescue notebook; it costs two rescues in ``all_genes`` (``RNU12``,
       ``EPB41L4A-AS2``, 0.003% of the rescued expression) and none in ``hvg5000``, so the
       775/129 counts already reported stand at the precision they are quoted to.*

    The assumption in guard 2 is that SCP542's annotation postdates each reassignment it is
    applied to. It cannot be checked -- SCP542 ships no annotation version -- which is the
    reason for resolving the doubt toward leaving the gene alone.
    """
    hgnc = pd.read_csv(hgnc_file, sep="\t", dtype=str, low_memory=False)
    approved = set(hgnc[CURRENT_FIELD].dropna())

    pairs = (
        hgnc[[CURRENT_FIELD, RENAME_FIELD]]
        .dropna(subset=[RENAME_FIELD])
        .assign(**{RENAME_FIELD: lambda d: d[RENAME_FIELD].str.split("|")})
        .explode(RENAME_FIELD)
    )
    pairs[RENAME_FIELD] = pairs[RENAME_FIELD].str.strip()
    pairs = pairs[pairs[RENAME_FIELD] != ""]

    counts = pairs[RENAME_FIELD].value_counts()
    refused = set(counts[counts > 1].index) | (set(pairs[RENAME_FIELD]) & approved)

    usable = pairs[~pairs[RENAME_FIELD].isin(refused)]
    return usable.set_index(RENAME_FIELD)[CURRENT_FIELD].to_dict(), refused


def annotate_hgnc_symbols(adata, hgnc_file: Path | str = HGNC_FILE) -> None:
    """Add ``var['hgnc_symbol']`` in place: the current HGNC symbol for each row.

    Rows that were never renamed, or whose rename is refused by ``load_rename_map``, carry
    their own ``var_name``, so the column is always a usable symbol and never empty.

    **Collisions are left unmapped, not merged.** A rename is dropped if its target is already
    a row of this matrix (SCP542 carried both names), or if two rows would be renamed onto the
    same target. Resolving either means combining two rows' expression values, which is a
    decision about the data rather than a lookup, and is deliberately not taken here
    (Selin, 05.08.2026). The affected rows keep their original symbol and stay
    out-of-vocabulary exactly as they are today.

    Collisions are matrix-dependent, so a variant with fewer genes can have fewer of them; the
    count is printed rather than assumed.
    """
    rename_map, _ = load_rename_map(hgnc_file)
    existing = set(adata.var_names)

    proposed = {g: rename_map[g] for g in adata.var_names if g in rename_map}
    target_counts = pd.Series(list(proposed.values())).value_counts() if proposed else pd.Series(dtype=int)
    blocked = {
        g
        for g, target in proposed.items()
        if target in existing or target_counts.get(target, 0) > 1
    }

    applied = {g: t for g, t in proposed.items() if g not in blocked}
    adata.var["hgnc_symbol"] = [applied.get(g, g) for g in adata.var_names]

    print(
        f"  hgnc_symbol: {len(applied):,} of {adata.n_vars:,} rows renamed to their current "
        f"symbol, {len(blocked):,} renames withheld as collisions"
    )
