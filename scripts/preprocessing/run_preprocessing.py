"""Orchestrate the full OncoTox preprocessing pipeline in the correct order.

Paths are derived once from ``--data-root`` (default in ``layout.py``) and
``--variant``, then passed explicitly to each step. Outputs live only under::

    <data-root>/processed/scRNAseq_SCP542/<variant>/

Expensive steps (convert, scGPT) refuse to overwrite existing files unless
``--overwrite`` is passed. ``hvg5000`` and ``all_genes`` never share a folder.

Examples::

    uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \\
        --start-at targets --skip-scgpt --scgpt-python ...

    uv run scripts/preprocessing/run_preprocessing.py --variant all_genes --all-drugs \\
        --scgpt-python ...
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.preprocessing import (  # noqa: E402
    add_pca,
    create_splits,
    ctrp_to_h5ad,
    scp542_conversion,
)
from scripts.preprocessing.layout import (  # noqa: E402
    DEFAULT_SCGPT_MODEL_DIR,
    DEFAULT_SCGPT_SCRIPT,
    PipelinePaths,
    VARIANT_N_TOP_GENES,
    add_data_args,
    guard_output,
)

STEP_ORDER = ["convert", "scgpt", "targets", "splits", "pca"]


def _print_step(idx: int, total: int, label: str) -> None:
    bar = "=" * 70
    print(f"\n{bar}\n[{idx}/{total}] {label}\n{bar}")


def _run_scgpt(
    scgpt_python: str | None,
    scgpt_script: Path,
    model_dir: Path,
    input_path: Path,
    output_path: Path,
    overwrite: bool,
) -> None:
    guard_output(output_path, overwrite=overwrite, step="scGPT embeddings")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not scgpt_script.exists():
        raise FileNotFoundError(f"scGPT script not found: {scgpt_script}")

    if scgpt_python:
        cmd = [
            scgpt_python,
            str(scgpt_script),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--model-dir",
            str(model_dir),
        ]
        print(f"Running scGPT via subprocess:\n  {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    else:
        print("-" * 70)
        print("MANUAL STEP: run gen_embeds.py, then press Enter.")
        print(
            f"  python {scgpt_script} --input {input_path} --output {output_path} "
            f"--model-dir {model_dir}"
        )
        print("-" * 70)
        input("Press Enter once embeddings exist... ")

    if not output_path.exists():
        raise RuntimeError(f"Expected scGPT output missing:\n  {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run the full OncoTox preprocessing pipeline in order.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_data_args(parser)
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=None,
        help="Override HVG count for convert (default follows --variant).",
    )
    parser.add_argument("--min-cell-lines", type=int, default=50)
    parser.add_argument("--all-drugs", action="store_true")
    parser.add_argument("--target-drug", default="paclitaxel")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-multi-split", action="store_true")
    parser.add_argument(
        "--scgpt-python",
        default=None,
        help="Python with scgpt installed (only needed when running the scGPT step).",
    )
    parser.add_argument(
        "--scgpt-script",
        type=Path,
        default=DEFAULT_SCGPT_SCRIPT,
        help="gen_embeds.py path.",
    )
    parser.add_argument(
        "--scgpt-model-dir",
        type=Path,
        default=DEFAULT_SCGPT_MODEL_DIR,
        help="scGPT_human weights directory.",
    )
    parser.add_argument("--skip-scgpt", action="store_true")
    parser.add_argument(
        "--force-pca",
        action="store_true",
        help="Recompute X_pca inside the targets h5ad (does not rebuild embeddings).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow convert/scGPT to replace existing raw or embedding h5ad files.",
    )
    parser.add_argument("--start-at", choices=STEP_ORDER, default="convert")

    args = parser.parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant)
    paths.processed_dir.mkdir(parents=True, exist_ok=True)

    default_hvg = VARIANT_N_TOP_GENES[args.variant]
    n_top = args.n_top_genes if args.n_top_genes is not None else default_hvg
    hvg = n_top if n_top and n_top > 0 else None
    min_cell_lines = 0 if args.all_drugs else args.min_cell_lines
    start_idx = STEP_ORDER.index(args.start_at)
    total = len(STEP_ORDER)

    print(f"data_root : {paths.data_root}")
    print(f"variant   : {paths.variant} -> {paths.processed_dir}")

    if start_idx <= STEP_ORDER.index("convert"):
        hvg_label = f"top-{hvg} HVGs" if hvg else "no HVG filter"
        _print_step(1, total, f"scp542_conversion ({hvg_label})")
        guard_output(paths.raw_h5ad, overwrite=args.overwrite, step="scp542_conversion")
        scp542_conversion.run(
            str(paths.expr_file),
            str(paths.meta_file),
            str(paths.raw_h5ad),
            hvg,
        )

    if start_idx <= STEP_ORDER.index("scgpt"):
        if args.skip_scgpt:
            _print_step(2, total, "scGPT (skipped)")
            if not paths.embed_h5ad.exists():
                raise RuntimeError(f"--skip-scgpt but missing:\n  {paths.embed_h5ad}")
        else:
            _print_step(2, total, "scGPT embeddings")
            _run_scgpt(
                args.scgpt_python,
                Path(args.scgpt_script),
                Path(args.scgpt_model_dir),
                paths.raw_h5ad,
                paths.embed_h5ad,
                args.overwrite,
            )

    if start_idx <= STEP_ORDER.index("targets"):
        _print_step(3, total, f"ctrp_to_h5ad (min_cell_lines={min_cell_lines})")
        if not paths.embed_h5ad.exists():
            raise RuntimeError(f"Missing embeddings input:\n  {paths.embed_h5ad}")
        ctrp_to_h5ad.run(
            str(paths.embed_h5ad),
            str(paths.targets_h5ad),
            str(paths.ctrp_dir),
            min_cell_lines=min_cell_lines,
            extra_single_drug_cols=(args.target_drug,) if args.target_drug else (),
        )

    if start_idx <= STEP_ORDER.index("splits"):
        _print_step(4, total, "create_splits")
        if not paths.targets_h5ad.exists():
            raise RuntimeError(f"Missing targets h5ad:\n  {paths.targets_h5ad}")
        if args.target_drug:
            create_splits.run(str(paths.targets_h5ad), args.target_drug, args.seed)
        if not args.skip_multi_split:
            create_splits.run_multi(str(paths.targets_h5ad), seed=args.seed)

    if start_idx <= STEP_ORDER.index("pca"):
        _print_step(5, total, f"add_pca (force={args.force_pca})")
        # PCA baseline is computed on the HVG-filtered convert counts (paths.raw_h5ad),
        # NOT the targets .X (which has had scGPT-OOV genes dropped) -- keeps PCA on the
        # single HVG filter for a clean PCA-vs-scGPT comparison.
        add_pca.run(
            str(paths.targets_h5ad),
            force=args.force_pca,
            counts_h5ad=str(paths.raw_h5ad),
        )

    print("\nPreprocessing pipeline complete.")
    print(f"Final file: {paths.targets_h5ad}")


if __name__ == "__main__":
    main()
