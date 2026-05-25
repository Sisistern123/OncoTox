"""Orchestrate the full OncoTox preprocessing pipeline in the correct order.

Steps (training is NOT included here):
    1. scp542_conversion : raw CPM + Metadata  ->  SCP542_CCLE.h5ad
                           (optionally filtered to top-N highly variable genes)
    2. gen_embeds.py     : SCP542_CCLE.h5ad    ->  SCP542_CCLE_scGPT_human_embeddings.h5ad
                           (lives in a separate repo, invoked via subprocess)
    3. ctrp_to_h5ad      : map CTRPv2 viability targets onto the embedded AnnData
                           -> SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad
    4. create_splits     : cell-line-grouped 70/15/15 train/val/test split
    5. add_pca           : PCA baseline embedding (X_pca)

Examples:
    # Full pipeline, top-5000 HVGs, with scGPT auto-run via the given interpreter
    python scripts/preprocessing/run_preprocessing.py \
        --n-top-genes 5000 \
        --scgpt-python /Users/selin/PycharmProjects/scGPT/.venv/bin/python

    # Same, but you'll run gen_embeds.py manually when prompted
    python scripts/preprocessing/run_preprocessing.py --n-top-genes 5000

    # Resume the pipeline starting at the splits step (e.g. after manually
    # re-running an earlier step)
    python scripts/preprocessing/run_preprocessing.py --start-at splits
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.preprocessing import (  # noqa: E402 - sys.path manipulation above
    add_pca,
    create_splits,
    ctrp_to_h5ad,
    scp542_conversion,
)

DATA_ROOT = Path("/Users/selin/Desktop/OncoTox")
EXPR_FILE = DATA_ROOT / "data/scRNAseq_SCP542/expression/CPM_data.txt"
META_FILE = DATA_ROOT / "data/scRNAseq_SCP542/metadata/Metadata.txt"
RAW_H5AD = DATA_ROOT / "data/scRNAseq_SCP542/metadata/SCP542_CCLE.h5ad"
EMBED_H5AD = DATA_ROOT / "data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings.h5ad"
TARGETS_H5AD = (
    DATA_ROOT
    / "data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
)
CTRP_DIR = DATA_ROOT / "data/metadata/CTRPv2.0_2015_ctd2_ExpandedDataset"

DEFAULT_SCGPT_SCRIPT = Path("/Users/selin/PycharmProjects/scGPT/gen_embeds.py")

STEP_ORDER = ["convert", "scgpt", "targets", "splits", "pca"]


def _print_step(idx: int, total: int, label: str) -> None:
    bar = "=" * 70
    print(f"\n{bar}\n[{idx}/{total}] {label}\n{bar}")


def _run_scgpt(scgpt_python: str | None, scgpt_script: Path) -> None:
    if not scgpt_script.exists():
        raise FileNotFoundError(
            f"scGPT embedding script not found at: {scgpt_script}\n"
            f"Pass --scgpt-script to point to the right path."
        )

    if scgpt_python:
        cmd = [scgpt_python, str(scgpt_script)]
        print(f"Running scGPT via subprocess:\n  {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    else:
        print()
        print("-" * 70)
        print("MANUAL STEP: generate scGPT embeddings.")
        print("Re-run this orchestrator with --scgpt-python <path> to automate, or")
        print("run the following in your scGPT env and then come back here:")
        print(f"  python {scgpt_script}")
        print()
        print(f"  expected input  : {RAW_H5AD}")
        print(f"  expected output : {EMBED_H5AD}")
        print("-" * 70)
        input("Press Enter once gen_embeds.py has completed successfully... ")

    if not EMBED_H5AD.exists():
        raise RuntimeError(
            f"Expected scGPT output not found after this step:\n  {EMBED_H5AD}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Run the full OncoTox preprocessing pipeline in order.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=5000,
        help="Number of highly variable genes to retain (0 disables HVG filtering).",
    )
    parser.add_argument("--target-drug", default="paclitaxel")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--scgpt-python",
        default=None,
        help="Python interpreter that has the scgpt package installed. "
        "If omitted, you will be prompted to run gen_embeds.py manually.",
    )
    parser.add_argument(
        "--scgpt-script",
        default=str(DEFAULT_SCGPT_SCRIPT),
        help="Path to gen_embeds.py in the scGPT repo.",
    )
    parser.add_argument(
        "--skip-scgpt",
        action="store_true",
        help="Skip scGPT embedding generation (assume embeddings already exist).",
    )
    parser.add_argument(
        "--force-pca",
        action="store_true",
        help="Recompute X_pca even if it already exists in the h5ad file.",
    )
    parser.add_argument(
        "--start-at",
        choices=STEP_ORDER,
        default="convert",
        help="Resume the pipeline starting at the given step.",
    )

    args = parser.parse_args()

    hvg = args.n_top_genes if args.n_top_genes and args.n_top_genes > 0 else None
    start_idx = STEP_ORDER.index(args.start_at)
    total = len(STEP_ORDER)

    if start_idx <= STEP_ORDER.index("convert"):
        hvg_label = f"with top-{hvg} HVGs" if hvg else "no HVG filter"
        _print_step(1, total, f"scp542_conversion ({hvg_label})")
        scp542_conversion.run(
            input_expr=str(EXPR_FILE),
            input_meta=str(META_FILE),
            output_path=str(RAW_H5AD),
            n_top_genes=hvg,
        )

    if start_idx <= STEP_ORDER.index("scgpt"):
        if args.skip_scgpt:
            _print_step(2, total, "scGPT embedding generation (SKIPPED via --skip-scgpt)")
            if not EMBED_H5AD.exists():
                raise RuntimeError(
                    f"--skip-scgpt was set but expected embeddings file is missing:\n"
                    f"  {EMBED_H5AD}"
                )
        else:
            _print_step(2, total, "scGPT embedding generation (gen_embeds.py)")
            _run_scgpt(args.scgpt_python, Path(args.scgpt_script))

    if start_idx <= STEP_ORDER.index("targets"):
        _print_step(3, total, f"ctrp_to_h5ad (drug='{args.target_drug}')")
        ctrp_to_h5ad.run(
            input_h5ad=str(EMBED_H5AD),
            output_h5ad=str(TARGETS_H5AD),
            ctrp_dir=str(CTRP_DIR),
            target_drug=args.target_drug,
        )

    if start_idx <= STEP_ORDER.index("splits"):
        _print_step(4, total, f"create_splits (drug='{args.target_drug}', seed={args.seed})")
        create_splits.run(
            h5ad_path=str(TARGETS_H5AD),
            target_drug=args.target_drug,
            seed=args.seed,
        )

    if start_idx <= STEP_ORDER.index("pca"):
        _print_step(5, total, f"add_pca (force={args.force_pca})")
        add_pca.run(h5ad_path=str(TARGETS_H5AD), force=args.force_pca)

    print("\nPreprocessing pipeline complete.")
    print(f"Final file: {TARGETS_H5AD}")


if __name__ == "__main__":
    main()
