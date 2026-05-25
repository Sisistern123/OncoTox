import argparse

import anndata as ad
import pandas as pd
import scanpy as sc

DEFAULT_EXPR = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/expression/CPM_data.txt"
DEFAULT_META = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/Metadata.txt"
DEFAULT_OUTPUT = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE.h5ad"


def run(
    input_expr: str = DEFAULT_EXPR,
    input_meta: str = DEFAULT_META,
    output_path: str = DEFAULT_OUTPUT,
    n_top_genes: int | None = None,
):
    """Build the foundational SCP542_CCLE.h5ad object from raw CPM + metadata.

    If `n_top_genes` is given, the AnnData is subset to that many highly variable
    genes (computed on a log1p copy so the saved .X keeps the original CPM values).
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

    if n_top_genes is not None and n_top_genes > 0:
        n_before = adata.n_vars
        if n_top_genes >= n_before:
            print(
                f"Requested n_top_genes={n_top_genes} >= total genes ({n_before}); "
                "skipping HVG filtering."
            )
        else:
            print(f"Selecting top {n_top_genes} highly variable genes...")
            # Compute HVGs on a log1p copy so the saved .X stays as CPM (what scGPT expects).
            adata_hvg = adata.copy()
            sc.pp.log1p(adata_hvg)
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
    parser.add_argument("--input-expr", default=DEFAULT_EXPR)
    parser.add_argument("--input-meta", default=DEFAULT_META)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=None,
        help="If set, filter to top N highly variable genes before saving.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run(args.input_expr, args.input_meta, args.output, args.n_top_genes)
