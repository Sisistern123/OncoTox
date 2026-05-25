import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

DEFAULT_INPUT = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings.h5ad"
DEFAULT_OUTPUT = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
DEFAULT_CTRP_DIR = "/Users/selin/Desktop/OncoTox/data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset"
DEFAULT_DRUG = "paclitaxel"


def run(
    input_h5ad: str = DEFAULT_INPUT,
    output_h5ad: str = DEFAULT_OUTPUT,
    ctrp_dir: str = DEFAULT_CTRP_DIR,
    target_drug: str = DEFAULT_DRUG,
):
    """Map CTRPv2 viability scores onto cells in the embedded AnnData."""
    print("Loading AnnData...")
    adata = sc.read_h5ad(input_h5ad)

    print("Loading CTRPv2 metadata...")
    ctrp_dir = Path(ctrp_dir)
    ctrp_values = pd.read_csv(
        ctrp_dir / "v20.data.per_cpd_post_qc.txt",
        sep="\t",
        usecols=["experiment_id", "master_cpd_id", "cpd_avg_pv"],
    )
    ctrp_exp_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_experiment.txt",
        sep="\t",
        usecols=["experiment_id", "master_ccl_id"],
    )
    ctrp_cell_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_cell_line.txt",
        sep="\t",
        usecols=["master_ccl_id", "ccl_name"],
    )
    ctrp_cpd_meta = pd.read_csv(
        ctrp_dir / "v20.meta.per_compound.txt",
        sep="\t",
        usecols=["master_cpd_id", "cpd_name"],
    )

    print("Merging CTRPv2 tables...")
    ctrp_full = (
        ctrp_values.merge(ctrp_exp_meta, on="experiment_id", how="inner")
        .merge(ctrp_cell_meta, on="master_ccl_id", how="inner")
        .merge(ctrp_cpd_meta, on="master_cpd_id", how="inner")
    )

    # Normalize cell line names (lowercase, no spaces/hyphens) to ensure perfect matching.
    ctrp_full["ccl_name_norm"] = (
        ctrp_full["ccl_name"].astype(str).str.strip().str.lower().str.replace("-", "")
    )
    ctrp_full["cpd_name_norm"] = ctrp_full["cpd_name"].astype(str).str.strip().str.lower()

    drug_data = ctrp_full[ctrp_full["cpd_name_norm"] == target_drug]
    # If there are multiple experiments for the same cell line and drug, take the mean viability.
    cell_line_to_viability = drug_data.groupby("ccl_name_norm")["cpd_avg_pv"].mean().to_dict()

    print(f"Mapping {target_drug} scores to single cells...")
    adata.obs["Cell_line_norm"] = (
        adata.obs["Cell_line"]
        .astype(str)
        .str.split("_")
        .str[0]
        .str.strip()
        .str.lower()
        .str.replace("-", "")
    )

    target_col = f"viability_{target_drug}"
    adata.obs[target_col] = adata.obs["Cell_line_norm"].map(cell_line_to_viability)

    # Cells with NaN should be ignored during training.
    mask_col = f"train_mask_{target_drug}"
    adata.obs[mask_col] = adata.obs[target_col].notna()

    print(f"Summary for {target_drug}:")
    print(f"  Total cells: {adata.n_obs}")
    print(f"  Cells with a valid viability score: {adata.obs[mask_col].sum()}")

    adata.obs = adata.obs.drop(columns=["Cell_line_norm"])

    print("Sanitizing metadata to fix H5AD string/index compatibility...")
    adata.obs.index = adata.obs.index.astype(str).astype(object)
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col]) or pd.api.types.is_object_dtype(
            adata.obs[col]
        ):
            adata.obs[col] = adata.obs[col].astype(str).astype(object)

    adata.var.index = adata.var.index.astype(str).astype(object)
    for col in adata.var.columns:
        if pd.api.types.is_string_dtype(adata.var[col]) or pd.api.types.is_object_dtype(
            adata.var[col]
        ):
            adata.var[col] = adata.var[col].astype(str).astype(object)

    print(f"Saving updated AnnData to {output_h5ad}...")
    adata.write_h5ad(output_h5ad, convert_strings_to_categoricals=False)
    print("Done!")
    return adata


def _parse_args():
    parser = argparse.ArgumentParser(description="Merge CTRPv2 viability targets into embedded AnnData.")
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--ctrp-dir", default=DEFAULT_CTRP_DIR)
    parser.add_argument("--drug", default=DEFAULT_DRUG)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run(args.input, args.output, args.ctrp_dir, args.drug)
