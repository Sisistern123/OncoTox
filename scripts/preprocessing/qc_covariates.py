"""Per-cell sequencing-depth covariates, recovered from the raw UMI counts.

SCP542 reaches this project as **CPM**, and every matrix downstream of that is normalised: the
library size a cell was sequenced to is divided out before anything here sees it. Two of the four
covariates the Q2 confound veto is defined on -- total counts and mitochondrial fraction
(`notebooks/4b_mil_training.ipynb` §2.5) -- are therefore not recoverable from any processed file.
They are recoverable from the study's own `UMIcount_data.txt`, which is what this module reads.

**Why the confound veto needs them at all.** Stages 1 and 2 establish that per-cell predictions vary
within a cell line and that independent seeds agree on *which* cells. Sequencing depth produces
exactly that pattern -- it varies within a line, and it reproduces across seeds perfectly, because it
is a property of the cell rather than of the model. Without these two columns the strongest rival
explanation for a positive Q2 result goes untested. `Genes_expressed`, already in `obs`, is partial
cover only: gene detection saturates, so two cells at 3x different depth can carry similar gene
counts.

**Over the full gene set, not the highly variable subset (Selin, 13.08.2026).** The decision is not
about precision:

* A total over the HVG genes is not depth. It is depth times the fraction of that cell's counts
  falling in the HVG set, and **that fraction is biological** -- a cell running an HVG-heavy program
  has a higher one by construction. A confound regressor carrying the signal under test can veto a
  true positive, which is worse than a noisy one.
* Of the 13 `MT-` genes in the raw matrix, **4** survive HVG selection (`ND1`, `ND3`, `ND6`, `CYB`);
  the nine dropped include `CO1`, `CO2`, `CO3`, `ATP6` and `ND4`, the high expressers. An HVG-
  restricted mitochondrial fraction is a truncated numerator over a denominator missing most of the
  transcriptome -- a different quantity, not a weaker version of the same one.

Both are read in one streaming pass, because the file is 3.5 GB and only two column-wise reductions
are wanted from it. The existing `scp542_conversion.run` loads its matrix whole; that is right when
the matrix becomes the `AnnData` and wrong here.

⚠️ **The layout is not the CPM file's.** `CPM_data.txt` begins with a `GENE` header and then gene
rows. `UMIcount_data.txt` has an unnamed first header field and carries **`Cell_line` and `Pool_ID`
as its first two data rows**, which are strings and would poison the accumulators. They are skipped
explicitly rather than coerced.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

#: Rows of ``UMIcount_data.txt`` that are metadata rather than genes, 0-indexed from the header.
UMI_METADATA_ROWS = [1, 2]

#: HGNC prefix for the mitochondrially encoded genes. 13 of them are in the raw matrix.
MT_PREFIX = "MT-"


def umi_depth_and_mito(
    umi_txt: str | Path,
    *,
    chunk_genes: int = 512,
    verbose: bool = True,
) -> pd.DataFrame:
    """Stream the raw UMI matrix once; return per-cell depth and mitochondrial fraction.

    The matrix is genes x cells, so a chunk of rows is a slab of genes across **all** cells and the
    reductions are column sums accumulated over chunks. Nothing larger than
    ``chunk_genes x n_cells`` is ever resident.

    :returns: DataFrame indexed by cell barcode with
        ``total_counts`` (the cell's library size over every gene in the file),
        ``mt_counts`` and ``pct_counts_mt`` (percent, the usual scale for this quantity), and
        ``n_genes_umi`` (genes with a non-zero count). The last is not a veto covariate -- it is
        carried because it is free in this pass and gives an independent check on the barcode join:
        it must track ``obs['Genes_expressed']``, which came from the study's own metadata by a
        different route.
    :raises ValueError: if no mitochondrial gene is found, which would silently yield a
        mitochondrial fraction of exactly zero for every cell.
    """
    umi_txt = Path(umi_txt)
    # No `dtype=`: it is applied to the index column as well, which holds gene symbols, so asking
    # for float32 fails on the first row with "could not convert string to float: 'A1BG'". The
    # counts are inferred as int64 and cast per chunk instead, which costs one copy of a slab and
    # nothing else.
    reader = pd.read_csv(
        umi_txt, sep="\t", index_col=0, skiprows=UMI_METADATA_ROWS, chunksize=chunk_genes,
    )

    cells: pd.Index | None = None
    total = mt = n_genes = None
    n_rows = n_mt_genes = 0
    for i, chunk in enumerate(reader, start=1):
        if cells is None:
            cells = chunk.columns
            total = np.zeros(len(cells), dtype=np.float64)
            mt = np.zeros(len(cells), dtype=np.float64)
            n_genes = np.zeros(len(cells), dtype=np.int64)
        elif not chunk.columns.equals(cells):
            raise ValueError(f"chunk {i} of {umi_txt.name} has a different cell order than chunk 1.")

        values = chunk.to_numpy(dtype=np.float32)
        is_mt = np.asarray(chunk.index.str.upper().str.startswith(MT_PREFIX))
        total += values.sum(axis=0, dtype=np.float64)
        n_genes += (values > 0).sum(axis=0)
        if is_mt.any():
            mt += values[is_mt].sum(axis=0, dtype=np.float64)
            n_mt_genes += int(is_mt.sum())
        n_rows += len(chunk)
        if verbose and i % 10 == 0:
            print(f"  {n_rows:6d} genes read, {n_mt_genes} mitochondrial")

    if n_mt_genes == 0:
        raise ValueError(
            f"no gene in {umi_txt.name} starts with {MT_PREFIX!r}, so the mitochondrial fraction "
            f"would be exactly zero for every cell -- which is a silent failure, not a measurement.")

    out = pd.DataFrame(
        {"total_counts": total, "mt_counts": mt,
         "pct_counts_mt": 100.0 * mt / np.where(total > 0, total, np.nan),
         "n_genes_umi": n_genes},
        index=pd.Index(cells, name="cell"),
    )
    if verbose:
        print(f"{n_rows} genes x {len(out)} cells | {n_mt_genes} mitochondrial genes | "
              f"median depth {out.total_counts.median():,.0f} | "
              f"median mito {out.pct_counts_mt.median():.2f}%")
    return out


def load_or_compute(
    umi_txt: str | Path,
    cache_csv: str | Path,
    *,
    recompute: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """:func:`umi_depth_and_mito`, cached to ``cache_csv``.

    **The cache is deliberately variant-independent**, and that is the whole reason it exists. These
    covariates are computed over the *full* gene set, so ``hvg5000`` and ``all_genes`` produce
    byte-identical values; caching per variant would scan 3.5 GB twice for one answer. It sits beside
    the variant directories rather than inside one, so nothing suggests it belongs to either.

    Unlike a cached PCA projection -- which ``cv.fold_pca_projections`` refuses to persist, because it
    is meaningful only alongside the fold assignment that produced it -- this is a **function of the
    raw file alone**. It cannot silently describe a different run, which is what makes caching it
    safe where caching that was not.
    """
    cache_csv = Path(cache_csv)
    if cache_csv.exists() and not recompute:
        if verbose:
            print(f"[qc] reading cached UMI covariates from {cache_csv.name}")
        return pd.read_csv(cache_csv, index_col=0)

    qc = umi_depth_and_mito(umi_txt, verbose=verbose)
    cache_csv.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(cache_csv)
    if verbose:
        print(f"[qc] wrote {cache_csv}")
    return qc


#: The obs columns this module adds. Named for what they are, and matching the convention scanpy's
#: own ``calculate_qc_metrics`` uses, so a reader who knows that function knows these.
QC_OBS_COLUMNS = ("total_counts", "pct_counts_mt")


def attach(adata, qc: pd.DataFrame, *, verbose: bool = True) -> None:
    """Join ``qc`` onto ``adata.obs`` by cell barcode, in place.

    **By barcode, never by position.** The UMI matrix holds 56,982 cells against the processed
    object's 53,513 -- preprocessing dropped 3,469 -- and its column order differs from ``obs_names``.
    A positional join would therefore align almost every cell to the wrong one while raising nothing.

    The join is checked in both directions, because either failure is silent: every cell must be
    found, and the recovered gene count must track the ``Genes_expressed`` the study's own metadata
    supplied. That second check is what distinguishes "the join worked" from "the file parsed" --
    the two quantities come from different files by different routes and cannot agree by accident.

    :raises ValueError: if any cell is missing from ``qc``, or if the two gene counts disagree
        badly enough that the join cannot be the intended one.
    """
    missing = adata.obs_names.difference(qc.index)
    if len(missing):
        raise ValueError(
            f"{len(missing)} of {adata.n_obs} cells are absent from the UMI covariates "
            f"(e.g. {list(missing[:3])}). The join is by barcode, so a missing cell means the two "
            f"files describe different cell sets, not that a value is unavailable.")

    aligned = qc.loc[adata.obs_names]
    if "Genes_expressed" in adata.obs.columns:
        theirs = pd.to_numeric(adata.obs["Genes_expressed"], errors="coerce").to_numpy(float)
        ours = aligned["n_genes_umi"].to_numpy(float)
        ok = np.isfinite(theirs) & np.isfinite(ours)
        r = float(np.corrcoef(theirs[ok], ours[ok])[0, 1]) if ok.sum() > 2 else float("nan")
        if not (r > 0.9):
            raise ValueError(
                f"genes detected from the UMI matrix correlate r={r:.3f} with obs['Genes_expressed'], "
                f"which came from the study's own metadata. Below 0.9 the barcode join is not the "
                f"intended one -- these are two measurements of the same quantity by different "
                f"routes and should agree almost exactly.")
        if verbose:
            print(f"[qc] join check: genes detected vs obs['Genes_expressed'] r={r:.4f}")

    for col in QC_OBS_COLUMNS:
        adata.obs[col] = aligned[col].to_numpy()
    if verbose:
        print(f"[qc] added {list(QC_OBS_COLUMNS)} | median depth "
              f"{adata.obs['total_counts'].median():,.0f} | median mito "
              f"{adata.obs['pct_counts_mt'].median():.2f}%")
