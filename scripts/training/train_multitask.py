"""Train one OncoMLP with K heads (one per CTRPv2 drug) using masked MSE.

Switch between the PCA baseline and the scGPT embedding with ``--use-rep``:

    # scGPT multi-task (matches the project_planning_v2.pdf next step)
    uv run scripts/training/train_multitask.py --use-rep X_scGPT

    # PCA baseline multi-task
    uv run scripts/training/train_multitask.py --use-rep X_pca

    # Few-drug intermediate (validates the masked-loss machinery on a small K
    # before scaling out to the full catalog -- recommended by the v2 plan).
    uv run scripts/training/train_multitask.py --use-rep X_scGPT \\
        --drugs paclitaxel docetaxel gemcitabine

    # All CTRPv2 drugs (requires preprocessing with --min-cell-lines 0).
    uv run scripts/training/train_multitask.py --use-rep X_scGPT

Requires ``ctrp_to_h5ad`` + ``create_splits --mode multi`` to have run, so the
h5ad file has Y_ctrp / M_ctrp obsm matrices, the ctrp_drugs uns list, and a
``split_ctrp`` obs column.

A per-drug-mean predictor (predicts the train-set mean viability per head) is
always evaluated alongside the model. This is the cheapest possible sanity
baseline: any head where the model fails to beat it has not learned anything
useful, regardless of the absolute MSE.
"""

from __future__ import annotations

import argparse

import numpy as np
import torch
from torch.utils.data import DataLoader

from scripts.model.OncoMLP import OncoMLP
from scripts.model.dataset import MultiDrugDataset
from scripts.training.training_utils import (
    TrainConfig,
    pick_device,
    train_model,
)

FILE_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"

DEFAULT_HIDDEN_DIMS = {
    "X_pca": (64, 32),
    "X_scGPT": (128, 64),
}


def _parse_args():
    parser = argparse.ArgumentParser(description="Multi-task CTRP drug-response training.")
    parser.add_argument("--path", default=FILE_PATH)
    parser.add_argument(
        "--use-rep",
        default="X_scGPT",
        choices=("X_pca", "X_scGPT"),
        help="Cell representation to feed the MLP.",
    )
    parser.add_argument(
        "--drugs",
        nargs="+",
        default=None,
        help="Restrict training to this drug subset (must be present in uns['ctrp_drugs']). "
        "Defaults to all drugs persisted by ctrp_to_h5ad.",
    )
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--dropout", type=float, default=0.5)
    parser.add_argument("--input-dropout", type=float, default=0.1)
    parser.add_argument("--loss", default="mse", choices=("mse", "huber"))
    parser.add_argument(
        "--hidden-dims",
        nargs="+",
        type=int,
        default=None,
        help="Override hidden dims (default: 64,32 for X_pca; 128,64 for X_scGPT).",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tag", default=None, help="Logging tag (default: --use-rep value).")
    parser.add_argument(
        "--baseline-topk",
        type=int,
        default=5,
        help="How many best/worst (model vs per-drug-mean) deltas to print.",
    )
    return parser.parse_args()


def _per_drug_train_mean(train_dataset: MultiDrugDataset) -> np.ndarray:
    """Per-drug train mean viability over observed (mask=1) entries only."""
    y = train_dataset.y.numpy()
    m = train_dataset.mask.numpy()
    sums = (y * m).sum(axis=0)
    counts = m.sum(axis=0)
    means = np.full(sums.shape, np.nan, dtype=np.float64)
    means[counts > 0] = sums[counts > 0] / counts[counts > 0]
    return means


def _per_drug_constant_mse(constants: np.ndarray, dataset: MultiDrugDataset) -> tuple[np.ndarray, np.ndarray]:
    """Per-drug MSE if we predict `constants[k]` for every cell of head k.

    Returns (mse_per_drug, counts_per_drug). MSE is NaN for heads with no
    observed val entries.
    """
    y = dataset.y.numpy()
    m = dataset.mask.numpy()
    safe_const = np.where(np.isnan(constants), 0.0, constants)
    preds = np.broadcast_to(safe_const[None, :], y.shape)
    sq = (preds - y) ** 2 * m
    sums = sq.sum(axis=0)
    counts = m.sum(axis=0)
    mse = np.full(sums.shape, np.nan, dtype=np.float64)
    mse[counts > 0] = sums[counts > 0] / counts[counts > 0]
    # Heads with no train support => no baseline prediction => NaN.
    mse[np.isnan(constants)] = np.nan
    return mse, counts


def _print_baseline_comparison(
    drug_names: list[str],
    baseline_mse: np.ndarray,
    model_mse: np.ndarray,
    counts: np.ndarray,
    topk: int,
    tag: str,
) -> None:
    """Print which heads the model beats the per-drug-mean baseline on."""
    finite = np.isfinite(baseline_mse) & np.isfinite(model_mse)
    if finite.sum() == 0:
        print(f"[{tag}] No head has both baseline and model val MSE; skipping comparison.")
        return

    deltas = model_mse - baseline_mse  # negative => model beats baseline.
    n_beats = int(np.logical_and(deltas < 0, finite).sum())
    n_total = int(finite.sum())
    mean_baseline = float(np.nanmean(baseline_mse))
    mean_model = float(np.nanmean(model_mse))

    print(f"\n[{tag}] Per-drug-mean baseline vs model (val):")
    print(f"  mean MSE over drugs : baseline={mean_baseline:.4f} | model={mean_model:.4f}")
    print(f"  heads beating baseline: {n_beats} / {n_total}")

    finite_idx = np.flatnonzero(finite)
    ranked_in_finite = finite_idx[np.argsort(deltas[finite_idx])]
    best = ranked_in_finite[: min(topk, n_total)]
    worst = ranked_in_finite[::-1][: min(topk, n_total)]

    def fmt(idx):
        return (
            f"{drug_names[idx]}: model={model_mse[idx]:.3f} "
            f"baseline={baseline_mse[idx]:.3f} d={deltas[idx]:+.3f} "
            f"(n={int(counts[idx])})"
        )

    print("  best  (largest model gains):")
    for i in best:
        print(f"    {fmt(i)}")
    print("  worst (model worse than baseline):")
    for i in worst:
        print(f"    {fmt(i)}")


def _evaluate_model_per_drug_mse(
    model: torch.nn.Module,
    val_loader: DataLoader,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-drug val MSE for the (already-best-state) model."""
    model.eval()
    sums = None
    counts = None
    with torch.no_grad():
        for batch_x, batch_y, batch_mask in val_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            batch_mask = batch_mask.to(device)
            preds = model(batch_x)
            sq = ((preds - batch_y) ** 2) * batch_mask
            s = sq.sum(dim=0).cpu().numpy()
            n = batch_mask.sum(dim=0).cpu().numpy()
            sums = s if sums is None else sums + s
            counts = n if counts is None else counts + n

    if sums is None or counts is None:
        return np.array([]), np.array([])
    mse = np.full(sums.shape, np.nan, dtype=np.float64)
    mse[counts > 0] = sums[counts > 0] / counts[counts > 0]
    return mse, counts


def main():
    args = _parse_args()

    hidden_dims = tuple(args.hidden_dims) if args.hidden_dims else DEFAULT_HIDDEN_DIMS[args.use_rep]
    tag = args.tag or args.use_rep

    train_dataset = MultiDrugDataset(
        h5ad_path=args.path, use_rep=args.use_rep, split="train", drugs=args.drugs
    )
    val_dataset = MultiDrugDataset(
        h5ad_path=args.path, use_rep=args.use_rep, split="val", drugs=args.drugs
    )

    if train_dataset.drug_names != val_dataset.drug_names:
        raise RuntimeError(
            "Train and val splits disagree on the drug column ordering; "
            "rerun ctrp_to_h5ad to regenerate Y_ctrp consistently."
        )

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False)

    sample_x, _, _ = train_dataset[0]
    input_dim = sample_x.shape[0]
    output_dim = len(train_dataset.drug_names)

    baseline_const = _per_drug_train_mean(train_dataset)
    baseline_mse, val_counts = _per_drug_constant_mse(baseline_const, val_dataset)
    print(
        f"[{tag}] Per-drug-mean sanity baseline ready: "
        f"mean baseline val MSE = {np.nanmean(baseline_mse):.4f} "
        f"over {int(np.isfinite(baseline_mse).sum())} / {output_dim} heads."
    )

    model = OncoMLP(
        input_dim=input_dim,
        hidden_dims=hidden_dims,
        dropout_rate=args.dropout,
        input_dropout=args.input_dropout,
        norm="layer",
        output_dim=output_dim,
    )

    config = TrainConfig(
        epochs=args.epochs,
        lr=args.lr,
        weight_decay=args.weight_decay,
        grad_clip=1.0,
        scheduler_patience=3,
        early_stop_patience=10,
        log_every=5,
        seed=args.seed,
        loss=args.loss,
    )

    print(
        f"Starting multi-task training: rep={args.use_rep}, K={output_dim} drugs, "
        f"input_dim={input_dim}, hidden_dims={hidden_dims}."
    )
    best_model, _history = train_model(
        model,
        train_loader,
        val_loader,
        config=config,
        tag=tag,
        drug_names=train_dataset.drug_names,
    )

    device = pick_device()
    model_mse, _ = _evaluate_model_per_drug_mse(best_model, val_loader, device)
    _print_baseline_comparison(
        drug_names=train_dataset.drug_names,
        baseline_mse=baseline_mse,
        model_mse=model_mse,
        counts=val_counts,
        topk=args.baseline_topk,
        tag=tag,
    )


if __name__ == "__main__":
    main()
