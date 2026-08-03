"""Generate scGPT cell embeddings for SCP542 (step 2 of ``run_preprocessing.py``).

**This script does not run under the OncoTox environment.** It imports ``scgpt``, which
lives in a separate checkout and virtualenv, so ``run_preprocessing.py`` invokes it as a
subprocess via ``--scgpt-python /path/to/scgpt-venv/bin/python``. It is vendored here
anyway so the embedding step -- the one part of the pipeline that was previously an
untracked file outside the repository -- is versioned with the results it produces.

Reads the ``convert`` output (``SCP542_CCLE.h5ad``), writes an h5ad carrying
``obsm["X_scGPT"]`` (512-d) with scGPT-out-of-vocabulary genes dropped from ``.X``, plus
an OOV gene table and summary alongside it.

Two behaviours worth knowing, both documented in
``docs/steps/02-preprocessing-and-embeddings.md``:

* ``max_length=1200`` matches the ``scGPT_human`` checkpoint's pretraining configuration.
  Cells with more non-zero genes have 1,199 of them **randomly sampled** -- the same
  operation used during pretraining (Cui et al., *Nature Methods* 21, 1470-1480, 2024).
* Embedding runs on MPS when available. That needs ``PYTORCH_ENABLE_MPS_FALLBACK=1`` set
  before ``torch`` is imported, which is why the ``os.environ`` line below precedes the
  imports.
"""

import argparse
import json
import os
from pathlib import Path

# PyTorch reads PYTORCH_ENABLE_MPS_FALLBACK when the MPS backend registers its dispatch keys,
# i.e. at torch import time -- so this must run before the scgpt import below, which pulls
# torch in. Without it, embedding on MPS dies on the one operator MPS does not implement:
#   NotImplementedError: aten::_nested_tensor_from_mask_left_aligned
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import torch  # noqa: E402
from scipy import sparse  # noqa: E402
from scgpt.tasks.cell_emb import embed_data  # noqa: E402

# Fixed so embeddings are reproducible; see the seeding block in main().
SEED = 42


def parse_args():
    parser = argparse.ArgumentParser(description="Generate scGPT cell embeddings for SCP542.")
    parser.add_argument("--input", required=True, help="Input SCP542_CCLE.h5ad path.")
    parser.add_argument("--output", required=True, help="Output embeddings h5ad path.")
    parser.add_argument("--model-dir", required=True, help="scGPT_human model directory.")
    parser.add_argument("--gene-col", default="index")
    return parser.parse_args()


def export_oov_genes(adata, model_dir, gene_col, input_file, oov_genes_file, oov_summary_file):
    vocab_file = Path(model_dir) / "vocab.json"
    with open(vocab_file, "r") as f:
        vocab_payload = json.load(f)

    vocab_tokens = set(vocab_payload.keys())
    gene_names = adata.var.index if gene_col == "index" else adata.var[gene_col]
    gene_names = pd.Series(gene_names, index=adata.var_names)
    in_vocab_mask = gene_names.isin(vocab_tokens).to_numpy()
    oov_mask = ~in_vocab_mask

    oov_var = adata.var.loc[oov_mask].copy()
    if oov_var.shape[0] == 0:
        print("All genes are in vocabulary. No OOV gene file written.")
        return

    oov_var["gene_name"] = gene_names.loc[oov_var.index].astype(str).values
    oov_var["in_vocab"] = False

    oov_X = adata[:, oov_mask].X
    if sparse.issparse(oov_X):
        n_cells_expressed = np.asarray((oov_X > 0).sum(axis=0)).ravel()
        total_expr = np.asarray(oov_X.sum(axis=0)).ravel()
    else:
        n_cells_expressed = np.asarray((oov_X > 0).sum(axis=0)).ravel()
        total_expr = np.asarray(oov_X.sum(axis=0)).ravel()

    n_cells = adata.n_obs
    pct_cells_expressed = (n_cells_expressed / max(n_cells, 1)) * 100.0
    mean_expr_if_expressed = np.divide(
        total_expr,
        n_cells_expressed,
        out=np.zeros_like(total_expr, dtype=float),
        where=n_cells_expressed > 0,
    )

    oov_var["n_cells_expressed"] = n_cells_expressed.astype(int)
    oov_var["pct_cells_expressed"] = pct_cells_expressed
    oov_var["total_expression"] = total_expr
    oov_var["mean_expression_if_expressed"] = mean_expr_if_expressed
    oov_var["source_input_file"] = input_file
    oov_var["model_vocab_file"] = str(vocab_file)
    oov_var["gene_column_used"] = gene_col

    print(f"Writing OOV gene metadata to {oov_genes_file}...")
    oov_var.to_csv(oov_genes_file, index=True)

    oov_summary = {
        "input_file": input_file,
        "model_vocab_file": str(vocab_file),
        "gene_column_used": gene_col,
        "n_total_genes": int(adata.n_vars),
        "n_oov_genes": int(oov_mask.sum()),
        "pct_oov_genes": float((oov_mask.sum() / max(adata.n_vars, 1)) * 100.0),
        "n_cells": int(adata.n_obs),
    }
    with open(oov_summary_file, "w") as f:
        json.dump(oov_summary, f, indent=2)
    print(f"Writing OOV summary to {oov_summary_file}...")


def run_embedding(adata, model_dir, gene_col, target_device):
    print(f"Running scGPT embedding on {target_device}...")
    return embed_data(
        adata_or_file=adata,
        model_dir=model_dir,
        gene_col=gene_col,
        batch_size=64,
        max_length=1200,
        device=target_device,
        return_new_adata=False,
    )


def sanitize_for_h5ad(adata):
    adata = adata.copy()

    def to_py_str(value):
        if pd.isna(value):
            return ""
        return str(value)

    def fix_df(df):
        clean_cols = {}
        for col in df.columns:
            series = df[col]
            if pd.api.types.is_numeric_dtype(series.dtype) or pd.api.types.is_bool_dtype(
                series.dtype
            ):
                clean_cols[col] = series.to_numpy()
            else:
                clean_cols[col] = np.array([to_py_str(v) for v in series.tolist()], dtype=object)

        clean_df = pd.DataFrame(clean_cols, index=df.index)
        for col in clean_df.columns:
            if not (
                pd.api.types.is_numeric_dtype(clean_df[col].dtype)
                or pd.api.types.is_bool_dtype(clean_df[col].dtype)
            ):
                clean_df[col] = clean_df[col].astype(object)
        clean_df.index = pd.Index([to_py_str(v) for v in clean_df.index.tolist()], dtype=object)
        return clean_df

    print("Sanitizing metadata for HDF5 compatibility...")
    adata.obs = fix_df(adata.obs)
    adata.var = fix_df(adata.var)

    if "X_scGPT" in adata.obsm:
        adata.obsm["X_scGPT"] = np.array(adata.obsm["X_scGPT"]).astype(np.float32)

    return adata


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    model_dir = args.model_dir
    gene_col = args.gene_col

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    oov_genes_file = str(Path(output_file).with_name(f"{Path(output_file).stem}_oov_genes.csv"))
    oov_summary_file = str(Path(output_file).with_name(f"{Path(output_file).stem}_oov_summary.json"))

    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"Using {device.upper()} for embedding.")
    if device == "mps":
        print(
            "  PYTORCH_ENABLE_MPS_FALLBACK=1: aten::_nested_tensor_from_mask_left_aligned "
            "is not implemented for MPS and runs on CPU; the rest runs natively on MPS."
        )

    # Two sources of randomness sit inside the embedding path, both seeded here:
    #   1. cells with more non-zero genes than max_length have their genes randomly
    #      subsampled (torch.randperm, scgpt/data_collator.py:169) -- this is the same
    #      operation used in pretraining, not a defect;
    #   2. value binning breaks ties with np.random.rand (scgpt/preprocess.py:_digitize).
    # The DataLoader runs with num_workers=0, so no per-worker reseeding is needed.
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    print(f"Seeded torch and numpy with {SEED}.")

    print(f"Loading AnnData from {input_file}...")
    adata = sc.read_h5ad(input_file)

    export_oov_genes(adata, model_dir, gene_col, input_file, oov_genes_file, oov_summary_file)
    adata = run_embedding(adata, model_dir, gene_col, device)
    adata = sanitize_for_h5ad(adata)

    print(f"Saving to {output_file}...")
    adata.write_h5ad(output_file, convert_strings_to_categoricals=False)
    print("Success.")


if __name__ == "__main__":
    main()
