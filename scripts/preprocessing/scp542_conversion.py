import argparse

import anndata as ad
import pandas as pd
import scanpy as sc

from scripts.preprocessing import qc_covariates
from scripts.preprocessing.expression import kinker_transform
from scripts.annotation.gene_symbols import annotate_hgnc_symbols
from scripts.layout import (
    PipelinePaths,
    VARIANT_N_TOP_GENES,
    add_data_args,
    guard_output,
)


def run(
    input_expr: str,
    input_meta: str,
    output_path: str,
    n_top_genes: int | None = None,
    umi_txt: str | None = None,
    umi_qc_csv: str | None = None,
):
    """Build the foundational SCP542_CCLE.h5ad object from raw CPM + metadata.

    If `n_top_genes` is given, the AnnData is subset to that many highly variable
    genes. The ranking runs on a transformed *copy* (``log2(1 + CPM/10)``, the SCP542
    authors' own quantification -- see ``expression.py``), so the saved ``.X`` keeps the
    original CPM values.

    ``var`` gains one column, ``hgnc_symbol`` -- SCP542's own identifiers stay in
    ``var_names`` and no expression value is touched, so the HVG set and ``X_pca`` are
    unaffected by it. Only ``gen_embeds.py`` reads it; see ``gene_symbols.py``.

    ``obs`` gains two columns, ``total_counts`` and ``pct_counts_mt``, when ``umi_txt`` is given
    (13.08.2026). They are read from SCP542's raw UMI matrix because ``input_expr`` is CPM, so the
    library size is already divided out of it, and they are the two covariates the Q2 confound veto
    is defined on that no processed file can supply -- see ``qc_covariates.py``. **No expression
    value is touched and no gene is added or removed**, so the HVG set, ``X_pca`` and every
    downstream number are unchanged by this; it is additive metadata.

    :param umi_txt: SCP542's ``UMIcount_data.txt``. Omitting it produces an object on which the
        confound veto cannot be evaluated, and says so loudly rather than failing.
    :param umi_qc_csv: where the per-cell covariates are cached. One scan of 3.5 GB serves every
        variant, because the covariates are computed over the full gene set.
    """
    print("Loading expression matrix... (this may take a few minutes and require high RAM)")
    df_expr = pd.read_csv(input_expr, sep="\t", index_col=0)

    adata = ad.AnnData(X=df_expr.T)

    print("Loading metadata...")
    df_meta = pd.read_csv(input_meta, sep="\t", low_memory=False)
    df_meta = df_meta.drop(0)
    df_meta = df_meta.set_index("NAME")

    print("Aligning metadata with expression data...")
    adata.obs = df_meta.loc[adata.obs_names]

    # Sequencing depth and mitochondrial fraction, from the raw UMI matrix (13.08.2026). They enter
    # here, with the rest of obs, and BEFORE the HVG filter below -- which is the point rather than a
    # convenience. Both are defined over the full gene set: a total restricted to the highly variable
    # genes is not depth but depth times the biological fraction of a cell's counts falling in that
    # set, and only 4 of the 13 MT- genes survive HVG selection. Computed once and cached, because
    # they are variant-independent for exactly that reason.
    if umi_txt is not None:
        qc = qc_covariates.load_or_compute(umi_txt, umi_qc_csv)
        qc_covariates.attach(adata, qc)
    else:
        print(
            "⚠️  No UMI matrix given: obs will carry NEITHER total_counts NOR pct_counts_mt, and the "
            "Q2 confound veto (4b_mil_training §2.5) cannot be evaluated on the result. Pass "
            "PipelinePaths.umi_file unless that is intended."
        )

    # Runs before HVG selection deliberately: the collision check then sees the full
    # transcriptome, so a rename is never applied to a symbol whose current holder exists in
    # the data but was dropped by HVG. It costs one recovered gene in hvg5000 against checking
    # each variant's own gene set, and none in the smaller variants (Selin, 05.08.2026).
    print("Annotating current HGNC symbols...")
    annotate_hgnc_symbols(adata)

    if n_top_genes is not None and n_top_genes > 0:
        n_before = adata.n_vars
        if n_top_genes >= n_before:
            print(
                f"Requested n_top_genes={n_top_genes} >= total genes ({n_before}); "
                "skipping HVG filtering."
            )
        else:
            print(f"Selecting top {n_top_genes} highly variable genes...")
            adata_hvg = adata.copy()
            # log2(1 + CPM/10), the transform the SCP542 authors use -- see expression.py.
            # flavor="seurat" needs log-scale input and reads uns["log1p"]["base"] to invert
            # it before computing dispersions.
            kinker_transform(adata_hvg)
            sc.pp.highly_variable_genes(
                adata_hvg,
                n_top_genes=n_top_genes,
                flavor="seurat",
            )
            hvg_mask = adata_hvg.var["highly_variable"].to_numpy()
            adata = adata[:, hvg_mask].copy()
            adata.uns["hvg_n_top_genes"] = int(n_top_genes)
            print(f"  Gene count: {n_before} -> {adata.n_vars}")

    print(f"Saving to {output_path}...")
    ad.settings.allow_write_nullable_strings = True
    adata.write(output_path)

    print(f"Success! Created AnnData object: {adata}")
    return adata


def _parse_args():
    parser = argparse.ArgumentParser(description="Build SCP542_CCLE.h5ad from raw CPM + metadata.")
    add_data_args(parser)
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=None,
        help="Override HVG count (default follows --variant).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace output h5ad if it already exists.",
    )
    args = parser.parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant, args.score)
    n_top = args.n_top_genes
    if n_top is None:
        n_top = VARIANT_N_TOP_GENES[args.variant]
    hvg = n_top if n_top and n_top > 0 else None
    guard_output(paths.raw_h5ad, overwrite=args.overwrite, step="scp542_conversion")
    paths.processed_dir.mkdir(parents=True, exist_ok=True)
    run(
        str(paths.expr_file),
        str(paths.meta_file),
        str(paths.raw_h5ad),
        hvg,
        umi_txt=str(paths.umi_file),
        umi_qc_csv=str(paths.umi_qc_csv),
    )


if __name__ == "__main__":
    _parse_args()
