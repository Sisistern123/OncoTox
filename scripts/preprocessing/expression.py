"""The expression transform, in one place so the HVG and PCA steps cannot diverge."""

from __future__ import annotations

import scanpy as sc

# Kinker et al. quantify expression as E[i,j] = log2(1 + CPM[i,j]/10), dividing CPM by 10
# because the average number of UMIs detected per cell is under 100,000; without it the
# difference between detected (E > 0) and undetected (E = 0) genes is inflated.
CPM_DIVISOR = 10.0
LOG_BASE = 2


def kinker_transform(adata) -> None:
    """Apply the SCP542 authors' expression transform in place: ``log2(1 + CPM/10)``.

    ``adata.X`` must hold raw CPM, which is what the portal distributes -- values in
    ``CPM_data.txt`` run into the tens and hundreds, so it is not the already-transformed
    ``E`` matrix.

    Source: Kinker et al., *Pan-cancer single-cell RNA-seq identifies recurring programs of
    cellular heterogeneity*, Nature Genetics 52, 1208-1218 (2020),
    doi:10.1038/s41588-020-00726-6, Methods, "Processing of scRNA-seq data". Following the
    dataset's own authors is preferred here over a generic ``log1p(CPM)``: their divisor is
    argued from a property of this data (UMIs per cell), and ``log1p(CPM)`` has no
    comparable justification for it.

    Why transform at all rather than feed CPM to PCA: Euclidean distance weighs each gene by
    its *absolute* spread, so a housekeeping gene at 5,000 CPM varying by 20% moves cells
    1,000 units apart while a transcription factor at 50 CPM varying by 100% moves them 50 --
    the leading components then read out ribosomal and mitochondrial content rather than
    biology. Expression differences are multiplicative; a log turns a fold change into a fixed
    additive distance wherever on the scale it happens. Per-gene scaling does not replace this:
    gene selection runs before it, and z-scoring a right-skewed distribution still leaves a few
    extreme cells dominating the gene.

    ``sc.pp.log1p(base=2)`` records ``uns['log1p']['base']``, which
    ``sc.pp.highly_variable_genes(flavor="seurat")`` reads in order to invert the transform
    before computing dispersions -- so the base must be set through scanpy, not by hand.
    """
    adata.X = adata.X / CPM_DIVISOR
    sc.pp.log1p(adata, base=LOG_BASE)
