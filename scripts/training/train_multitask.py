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

import scanpy as sc
from sklearn.model_selection import GroupKFold

from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS as OncoMLP_DEFAULT_HIDDEN_DIMS
from scripts.model.OncoMLP import OncoMLP, init_head_bias_
from scripts.model.dataset import MultiDrugDataset
from scripts.training.cv import grouped_folds, inner_holdout, per_drug_line_mean
from scripts.training.density_weighting import line_level
from scripts.training.training_utils import (
    TrainConfig,
    create_run_dir,
    pick_device,
    save_run,
    train_model,
    utc_now_iso,
)

from scripts.layout import PipelinePaths, add_data_args

# DEFAULT_HIDDEN_DIMS now lives in scripts/model/OncoMLP.py (it is a property of the architecture, and
# putting it there lets scripts.training.cv read it without an import cycle). Re-exported so the many
# notebooks that import it from here keep working.
DEFAULT_HIDDEN_DIMS = OncoMLP_DEFAULT_HIDDEN_DIMS


def _parse_args():
    parser = argparse.ArgumentParser(description="Multi-task CTRP drug-response training.")
    add_data_args(parser)
    parser.add_argument(
        "--path",
        type=str,
        default=None,
        help="Override targets h5ad (default: derived from --data-root and --variant).",
    )
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
    parser.add_argument("--loss", default="mse", choices=("mse", "mae", "huber"))
    parser.add_argument(
        "--hidden-dims",
        nargs="+",
        type=int,
        default=None,
        help="Override hidden dims (default: 128,64 for both X_pca and X_scGPT -- the trunk is "
             "matched so that only the representation differs).",
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


def _to_lines(
    dataset: MultiDrugDataset, preds: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Collapse a dataset from cells to cell lines: ``(y_lines, obs_lines, pred_lines)``.

    **Everything this module scores is scored per cell line, not per cell (audit 11, 12.08.2026).**
    The label is one measurement per (cell line, drug), broadcast onto that line's cells, so a
    cell-level mean weighs a 1,990-cell line 35x a 56-cell line for the same single measurement --
    the imbalance docs/steps/03 records as a factor of 82 across the loss. Until 12.08.2026 these
    metrics were per cell, which meant the "heads beating baseline" counts were reported at the
    resolution that same document calls dishonest.

    Labels collapse by lookup (:func:`density_weighting.line_level`), since they are constant within
    a line; predictions collapse by **averaging** the line's cells, which is what
    :func:`cv.line_level_predictions` already did for the out-of-fold path. The two now agree.
    """
    if dataset.groups is None:
        raise ValueError(
            "dataset has no cell-line labels, so it cannot be scored per line. Check that "
            "obs['Cell_line'] exists -- scoring per cell is deliberately not offered as a fallback, "
            "because it is the defect this function exists to remove."
        )
    groups = dataset.groups
    lines = np.unique(groups)
    y_lines, obs_lines = line_level(
        dataset.y.numpy(), dataset.mask.numpy().astype(bool), groups, lines
    )
    pred_lines = None
    if preds is not None:
        pred_lines = np.vstack([preds[groups == ln].mean(axis=0) for ln in lines])
    return y_lines, obs_lines, pred_lines


def _per_drug_train_mean(train_dataset: MultiDrugDataset) -> np.ndarray:
    """Per-drug train mean response, averaged over **cell lines**, over observed entries only."""
    y_lines, obs_lines, _ = _to_lines(train_dataset)
    return per_drug_line_mean(y_lines, obs_lines).astype(np.float64)


def _per_drug_constant_mse(constants: np.ndarray, dataset: MultiDrugDataset) -> tuple[np.ndarray, np.ndarray]:
    """Per-drug MSE over **cell lines** if we predict ``constants[k]`` for head k.

    Returns ``(mse_per_drug, n_lines_per_drug)`` -- the count is now cell lines, not cells. NaN for
    heads with no observed line here, and for heads whose constant is NaN for want of train support.
    """
    y_lines, obs_lines, _ = _to_lines(dataset)
    safe_const = np.where(np.isnan(constants), 0.0, constants)
    sq = (safe_const[None, :] - y_lines) ** 2
    sums = np.where(obs_lines, sq, 0.0).sum(axis=0)
    counts = obs_lines.sum(axis=0)
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
    val_dataset: MultiDrugDataset,
    device: torch.device,
    batch_size: int = 256,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-drug val MSE for the (already-best-state) model, **over cell lines**.

    Predictions are made per cell and then averaged to one value per line before being scored, which
    is the same collapse :func:`cv.line_level_predictions` performs. Returns
    ``(mse_per_drug, n_lines_per_drug)``.

    Takes the dataset rather than a ``DataLoader`` (changed 12.08.2026, audit 11): the loader has no
    idea which line each cell came from, and scoring per line needs that. The loader it builds is
    unshuffled, so rows stay aligned with ``dataset.groups``.
    """
    model.eval()
    loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    chunks = []
    with torch.no_grad():
        for batch_x, _, _ in loader:
            chunks.append(model(batch_x.to(device)).cpu().numpy())
    if not chunks:
        return np.array([]), np.array([])

    y_lines, obs_lines, pred_lines = _to_lines(val_dataset, np.vstack(chunks))
    sq = (pred_lines - y_lines) ** 2
    sums = np.where(obs_lines, sq, 0.0).sum(axis=0)
    counts = obs_lines.sum(axis=0)
    mse = np.full(sums.shape, np.nan, dtype=np.float64)
    mse[counts > 0] = sums[counts > 0] / counts[counts > 0]
    return mse, counts


def train_rep(
    *,
    use_rep: str,
    h5ad_path: str,
    config: TrainConfig,
    drugs: list[str] | None = None,
    hidden_dims: tuple[int, ...] | None = None,
    batch_size: int = 128,
    dropout: float = 0.5,
    input_dropout: float = 0.1,
    init_head_bias: bool = True,
    data_root: str | None = None,
    variant: str | None = None,
    tag: str | None = None,
    baseline_topk: int = 5,
    print_comparison: bool = True,
) -> dict:
    """Train one multi-task OncoMLP for ``use_rep`` and persist a run dir.

    This is the single source of truth for a training run; both the CLI
    (``main``) and ``notebooks/4a_percell_training.ipynb`` call it so they cannot drift.

    ``init_head_bias`` starts each head at that drug's mean over the **train split's cell
    lines** (``OncoMLP.init_head_bias_``). Added 12.08.2026: this path had never done it, so on a
    target centred near 0.9 it trained against an offset ``cv.oof_predictions`` did not.

    Returns a results dict with ``run_dir``, ``summary``, ``history``, the
    per-drug val MSE arrays (model + per-drug-mean baseline), ``drug_names``,
    and ``input_dim`` / ``output_dim`` for in-notebook plotting.
    """
    if hidden_dims is None:
        hidden_dims = DEFAULT_HIDDEN_DIMS[use_rep]
    hidden_dims = tuple(hidden_dims)
    tag = tag or use_rep

    train_dataset = MultiDrugDataset(
        h5ad_path=h5ad_path, use_rep=use_rep, split="train", drugs=drugs
    )
    val_dataset = MultiDrugDataset(
        h5ad_path=h5ad_path, use_rep=use_rep, split="val", drugs=drugs
    )
    # What was actually read: on a fixed split X_pca resolves to the train-fitted key.
    # Recorded separately from `rep` so run tags stay stable while provenance is exact.
    rep_key = train_dataset.use_rep
    if rep_key != use_rep:
        print(f"  rep '{use_rep}' resolved to obsm['{rep_key}'] (fitted on train cells only).")

    if train_dataset.drug_names != val_dataset.drug_names:
        raise RuntimeError(
            "Train and val splits disagree on the drug column ordering; "
            "rerun ctrp_to_h5ad to regenerate Y_ctrp consistently."
        )

    # Explicit generator so the shuffle order does not depend on train_model happening to
    # call set_seed before the loader is first iterated.
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        generator=torch.Generator().manual_seed(config.seed),
    )
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

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
        dropout_rate=dropout,
        input_dropout=input_dropout,
        norm="layer",
        output_dim=output_dim,
    )
    if init_head_bias:
        # Train split only. Per line, not per cell -- see cv.per_drug_line_mean.
        if train_dataset.groups is None:
            raise ValueError(
                "train_dataset has no cell-line labels, so the per-drug means cannot be taken per "
                "line. Check that obs['Cell_line'] exists, or pass init_head_bias=False."
            )
        y_lines, obs_lines = line_level(
            train_dataset.y.numpy(),
            train_dataset.mask.numpy().astype(bool),
            train_dataset.groups,
            np.unique(train_dataset.groups),
        )
        init_head_bias_(model, per_drug_line_mean(y_lines, obs_lines))

    if drugs:
        scope = "subset"
    elif output_dim <= 1:
        scope = "single_drug"
    else:
        scope = "all_drugs"
    run_tag = f"multitask_{use_rep}_{scope}"
    if scope == "subset":
        run_tag += f"_K{output_dim}"
    run_dir = create_run_dir(run_tag)
    started_at = utc_now_iso()

    print(
        f"Starting multi-task training: rep={use_rep}, K={output_dim} drugs, "
        f"input_dim={input_dim}, hidden_dims={hidden_dims}."
    )
    best_model, history = train_model(
        model,
        train_loader,
        val_loader,
        config=config,
        tag=tag,
        drug_names=train_dataset.drug_names,
    )

    device = pick_device()
    model_mse, _ = _evaluate_model_per_drug_mse(best_model, val_dataset, device)
    if print_comparison:
        _print_baseline_comparison(
            drug_names=train_dataset.drug_names,
            baseline_mse=baseline_mse,
            model_mse=model_mse,
            counts=val_counts,
            topk=baseline_topk,
            tag=tag,
        )

    summary = save_run(
        run_dir=run_dir,
        tag=run_tag,
        config=config,
        history=history,
        model=best_model,
        run_meta={
            "scope": scope,
            "drug_scope_kind": "multi_drug",
            "drugs_requested": drugs,
            "rep": use_rep,
            "rep_key": rep_key,
            "data_root": str(data_root) if data_root is not None else None,
            "variant": variant,
            "h5ad_path": h5ad_path,
            "input_dim": input_dim,
            "output_dim": output_dim,
            "hidden_dims": list(hidden_dims),
            "dropout_rate": dropout,
            "input_dropout": input_dropout,
            "norm": "layer",
            "init_head_bias": init_head_bias,
            "no_decay_bias_and_norm": config.no_decay_bias_and_norm,
            "batch_size": batch_size,
            "loss": config.loss,
            "n_train_cells": len(train_dataset),
            "n_val_cells": len(val_dataset),
            "script": "scripts/training/train_multitask.py",
        },
        started_at=started_at,
        drug_names=train_dataset.drug_names,
        model_per_drug_val_mse=model_mse,
        baseline_per_drug_val_mse=baseline_mse,
        n_val_per_drug=val_counts,
    )

    return {
        "run_dir": run_dir,
        "summary": summary,
        "history": history,
        "model_per_drug_val_mse": model_mse,
        "baseline_per_drug_val_mse": baseline_mse,
        "n_val_per_drug": val_counts,
        "drug_names": train_dataset.drug_names,
        "input_dim": input_dim,
        "output_dim": output_dim,
        "rep": use_rep,
        "rep_key": rep_key,
    }


def cv_evaluate(
    *,
    use_rep: str,
    config: TrainConfig,
    h5ad_path: str | None = None,
    adata=None,
    n_splits: int = 5,
    drugs: list[str] | None = None,
    hidden_dims: tuple[int, ...] | None = None,
    batch_size: int = 128,
    dropout: float = 0.5,
    input_dropout: float = 0.1,
    init_head_bias: bool = True,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
) -> list[dict]:
    """K-fold **GroupKFold** cross-validation over cell lines for one rep.

    Answers "is the PCA-vs-scGPT difference real, or a one-split artifact?" by
    re-fitting on `n_splits` cell-line-grouped folds and returning per-fold
    metrics (so the caller can report mean ± std). Grouping is by ``group_col``
    (cell line), so no line appears in both train and val of a fold — the same
    leakage control as the fixed `split_ctrp`, but resampled `n_splits` ways.

    Only cells whose `split_ctrp` is in ``eligible_splits`` are used. The default
    ``("train", "val")`` **holds the fixed test set out entirely** (CV resamples
    only the 153 train+val lines); pass ``("train", "val", "test")`` to pool all
    180 measured lines. The h5ad is read once and sliced per fold via ``cell_mask``.

    Within a fold, 15 % of the training lines are withheld for early stopping
    (``cv.inner_holdout``) so the scored fold decides nothing about the model that
    predicts it. The per-drug-mean baseline is fitted on the fitting lines for the
    same reason.

    ``init_head_bias`` starts each head at that drug's mean over the fold's fitting lines
    (``OncoMLP.init_head_bias_``), which ``cv.oof_predictions`` has always done and this
    function did not until 12.08.2026 — so the two CV paths were not running the same
    experiment on a target centred near 0.9.

    Returns a list of per-fold dicts: best_val_mse, train_mse_at_best, gap
    (early-stopping − fitting MSE at the best epoch), n_beats / n_total (heads
    beating the per-drug-mean baseline **on the scored fold**), model_mean_mse,
    baseline_mean_mse, and fold line counts.
    """
    if hidden_dims is None:
        hidden_dims = DEFAULT_HIDDEN_DIMS[use_rep]
    hidden_dims = tuple(hidden_dims)
    if adata is None:
        if h5ad_path is None:
            raise ValueError("Provide either h5ad_path or a preloaded adata.")
        adata = sc.read_h5ad(h5ad_path)

    if group_col not in adata.obs.columns:
        raise ValueError(f"group_col '{group_col}' not in adata.obs.")
    groups_all = adata.obs[group_col].astype(str).to_numpy()
    device = pick_device()

    # Same partition helper as scripts.training.cv.oof_predictions, so metrics computed here and
    # predictions produced there refer to identical folds.
    idx, fold_split = grouped_folds(
        adata, n_splits=n_splits, group_col=group_col, eligible_splits=eligible_splits
    )
    folds: list[dict] = []
    for fold, (tr, va) in enumerate(fold_split, start=1):
        train_cells = np.zeros(adata.n_obs, dtype=bool)
        val_cells = np.zeros(adata.n_obs, dtype=bool)
        train_cells[idx[tr]] = True
        val_cells[idx[va]] = True

        # The training lines split once more: the model is fitted on `fit_cells`, early stopping
        # watches `stop_cells`, and the fold `val_cells` is only ever scored. Until 12.08.2026 the
        # scored fold was also the early-stopping set, so every metric below was a minimum over
        # epochs on its own scored data (docs/steps/03, "The early-stopping set is nested").
        fit_cells, stop_cells = inner_holdout(groups_all, train_cells)

        fit_ds = MultiDrugDataset(adata=adata, use_rep=use_rep, cell_mask=fit_cells, drugs=drugs)
        stop_ds = MultiDrugDataset(adata=adata, use_rep=use_rep, cell_mask=stop_cells, drugs=drugs)
        val_ds = MultiDrugDataset(adata=adata, use_rep=use_rep, cell_mask=val_cells, drugs=drugs)
        fit_loader = DataLoader(
            fit_ds,
            batch_size=batch_size,
            shuffle=True,
            generator=torch.Generator().manual_seed(config.seed),
        )
        stop_loader = DataLoader(stop_ds, batch_size=batch_size, shuffle=False)
        # No val loader: the scored fold is read through `val_ds` so that predictions can be
        # averaged to one value per cell line before scoring (audit 11). A loader cannot say which
        # line a row came from.

        baseline_const = _per_drug_train_mean(fit_ds)
        baseline_mse, _ = _per_drug_constant_mse(baseline_const, val_ds)

        model = OncoMLP(
            input_dim=fit_ds.X.shape[1],
            hidden_dims=hidden_dims,
            dropout_rate=dropout,
            input_dropout=input_dropout,
            norm="layer",
            output_dim=len(fit_ds.drug_names),
        )
        if init_head_bias:
            # Fitting lines only -- the early-stopping slice and the scored fold are both excluded,
            # for the same reason the baseline constant above is fitted here.
            g_fit = groups_all[fit_cells]
            y_lines, obs_lines = line_level(
                fit_ds.y.numpy(), fit_ds.mask.numpy().astype(bool), g_fit, np.unique(g_fit)
            )
            init_head_bias_(model, per_drug_line_mean(y_lines, obs_lines))
        tag = f"cv{fold}/{n_splits}_{use_rep}"
        best_model, history = train_model(
            model, fit_loader, stop_loader, config=config, tag=tag,
            drug_names=fit_ds.drug_names,
        )
        model_mse, _ = _evaluate_model_per_drug_mse(best_model, val_ds, device)

        finite = np.isfinite(baseline_mse) & np.isfinite(model_mse)
        n_total = int(finite.sum())
        delta = model_mse[finite] - baseline_mse[finite]  # per-drug model − baseline MSE
        n_beats = int((delta < 0).sum())                  # heads where model beats the constant
        be = history.best_epoch
        folds.append({
            "fold": fold,
            "rep": use_rep,
            "n_train_lines": int(np.unique(groups_all[idx][tr]).size),
            "n_fit_lines": int(np.unique(groups_all[fit_cells]).size),
            "n_stop_lines": int(np.unique(groups_all[stop_cells]).size),
            "n_val_lines": int(np.unique(groups_all[idx][va]).size),
            # best_val_mse, train_mse_at_best and gap all come from the training curve, so they are
            # measured on the fitting and early-stopping slices -- NOT on the scored fold. The
            # numbers that describe the fold are model_mean_mse / baseline_mean_mse and the
            # heads-beating counts below, which are evaluated on val_ds.
            "best_val_mse": float(history.best_val_mse),
            "train_mse_at_best": float(history.train_mse[be - 1]),
            "gap": float(history.val_mse[be - 1] - history.train_mse[be - 1]),
            "n_beats": n_beats,
            "n_total": n_total,
            # Continuous counterpart of heads-beating: mean/median per-drug delta.
            # Negative => model better than the per-drug-mean baseline on average.
            "mean_delta": float(delta.mean()) if delta.size else float("nan"),
            "median_delta": float(np.median(delta)) if delta.size else float("nan"),
            "frac_beat": float((delta < 0).mean()) if delta.size else float("nan"),
            "model_mean_mse": float(np.nanmean(model_mse)),
            "baseline_mean_mse": float(np.nanmean(baseline_mse)),
        })
    return folds


def main():
    args = _parse_args()
    paths = PipelinePaths.build(args.data_root, args.variant, args.score)
    h5ad_path = args.path or str(paths.targets_h5ad)

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

    train_rep(
        use_rep=args.use_rep,
        h5ad_path=h5ad_path,
        config=config,
        drugs=args.drugs,
        hidden_dims=tuple(args.hidden_dims) if args.hidden_dims else None,
        batch_size=args.batch_size,
        dropout=args.dropout,
        input_dropout=args.input_dropout,
        data_root=paths.data_root,
        variant=paths.variant,
        tag=args.tag,
        baseline_topk=args.baseline_topk,
    )


if __name__ == "__main__":
    main()
