"""Shared training utilities for OncoTox MLP heads.

Contains both the training loop and the run-versioning layer:
  * ``TrainConfig`` / ``TrainHistory`` dataclasses.
  * ``train_model(...)`` -- single-task + multi-task masked-loss in one fn.
  * ``set_seed``, ``pick_device``.
  * ``create_run_dir`` / ``save_run`` -- writes ``runs/<ts>_<tag>/`` with
    config.json, run_meta.json, history.csv, summary.json, best_model.pt
    (+ per_drug_results.csv for multi-task), and appends to
    ``runs/runs_index.csv``.
"""

from __future__ import annotations

import copy
import csv
import json
import platform
import random
import socket
import sys
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader


@dataclass
class TrainConfig:
    # 25 is a safe cap: over 36 recorded runs the best epoch was median 6, max 11, and
    # early stopping (patience 10) has never reached 25. 50 only cost wall-clock.
    #
    # The default is NOT what the training notebook passes. `4a_percell_training` sets
    # epochs=50 for both of its sections (Selin, 12.08.2026, decided when the 8-run matrix
    # was folded in) so that one
    # config covers both rather than 25 for one and 50 for the other. Left at 25 here deliberately:
    # whether the default itself should move is review item 10, which owns the epoch question, and
    # the evidence above is a fixed-split measurement that the inner-holdout change may shift --
    # each fold now fits on ~104 lines rather than 122, so where training peaks can move.
    epochs: int = 25
    lr: float = 1e-3
    # 0.0, and this is a decision rather than a default (Selin, 12.08.2026). **The model trains
    # without weight decay**, with dropout 0.5 and input dropout 0.1 as the only regularizers.
    #
    # It follows from the AdamW switch below. Measured before deciding: AdamW at the inherited 1e-3
    # gives a weight-matrix L2 norm of 8.971 against 8.975 with decay off -- 0.04 %, i.e. nothing --
    # where Adam at the same nominal value gave 5.597, a 38 % shrinkage. Matching Adam under AdamW
    # needs wd ~= 1.0, three orders of magnitude up. So 1e-3 could not be carried across: it would
    # have read as a deliberate setting while doing nothing at all.
    #
    # The rejected alternative was wd ~= 1.0, which reverse-engineers a value to reproduce the
    # shrinkage of runs that are themselves void on target and panel grounds.
    #
    # ⚠️ **Whether this model needs weight decay is OPEN, not settled.** Do not justify 0.0 by "the
    # model is not regularization-limited": that refutation rested on the weight-decay axis of an
    # ablation which is itself void, so the evidence for it is gone. The honest justification is
    # narrower -- no sourced value exists for either optimizer, dropout is already substantial, and
    # a decay nobody can justify is worse than none.
    weight_decay: float = 0.0
    grad_clip: float | None = 1.0
    scheduler_patience: int = 3
    scheduler_factor: float = 0.5
    early_stop_patience: int = 10
    log_every: int = 5
    seed: int = 42
    loss: str = "mse"  # "mse" or "mae" -- the two arms of review item 9A
    # Decay the weight matrices, exempt the biases and the normalization parameters. This is the
    # standard grouping -- HuggingFace transformers' `Trainer.create_optimizer`
    # (`no_decay = ["bias", "LayerNorm.weight"]`), inherited from the BERT reference
    # implementation -- and it is on by default because it is the convention, not because this
    # target needs it. It happens to matter more than usual here: `auc_cc` is centred near 0.9,
    # so each head bias must sit at its drug's mean and decay pulls it toward 0.
    #
    # Replaced `exclude_output_from_decay` (12.08.2026), which exempted the whole output
    # `Linear` -- its weight matrix included -- and left the LayerNorm parameters and the hidden
    # biases decayed. That is neither the convention nor what docs/steps/03 described.
    no_decay_bias_and_norm: bool = True
    log_per_drug_topk: int = 5  # how many best/worst per-drug val MSEs to print


@dataclass
class TrainHistory:
    train_mse: list[float] = field(default_factory=list)
    val_mse: list[float] = field(default_factory=list)
    lr: list[float] = field(default_factory=list)
    best_val_mse: float = float("inf")
    best_epoch: int = -1
    # Multi-task only; remains empty in the single-target case.
    per_drug_val_mse: list[np.ndarray] = field(default_factory=list)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)


def pick_device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def _is_multitask_loader(loader: DataLoader) -> bool:
    """Detect masked multi-task loaders by peeking at one batch."""
    try:
        sample = next(iter(loader))
    except StopIteration:
        return False
    return isinstance(sample, (list, tuple)) and len(sample) == 3


def _masked_mean(per_elem: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Mean of ``per_elem`` over entries where ``mask > 0``.

    Returns 0.0 when no entries are observed in the batch, to avoid NaNs.
    """
    denom = mask.sum().clamp_min(1.0)
    return (per_elem * mask).sum() / denom


def _make_loss_fn(config: TrainConfig, multitask: bool):
    """Return a callable ``loss_fn(preds, targets, mask=None) -> scalar tensor``.

    Single-task callers can pass ``mask=None`` and get the plain loss.
    Multi-task callers must pass ``mask`` (same shape as ``targets``).

    **Both losses go through the identical masked reduction**, ``sum(err * M) / sum(M)``, with
    the density weights substituting for the 0/1 mask when weighting is on. That is what makes them
    comparable as arms of one experiment: an MAE that bypassed the mask would be a different
    estimator rather than a different loss, and the comparison in review item 9A would not mean what
    it says. Only ``per_elem`` differs between them.
    """
    name = config.loss.lower()
    if name not in {"mse", "mae"}:
        raise ValueError(f"Unknown loss: {config.loss!r} (expected 'mse' or 'mae')")

    if name == "mse":
        per_elem = lambda preds, targets: (preds - targets) ** 2
    else:
        # Added 12.08.2026 (Selin) as the second arm of item 9A's loss comparison, and since
        # 12.08.2026 the only alternative to MSE: Huber was removed rather than kept, because its
        # position in the grid is set entirely by `beta` -- so keeping it meant deciding a parameter
        # -- and at the 0.05 it carried it behaves close to MAE, so the grid would have run two
        # near-duplicate columns. MAE carries no parameter at all, which is half its appeal.
        # Less pulled by the sparse extremes than squared error, which is the same axis the
        # density weighting acts on -- so the two are not independent and the grid is read as a
        # grid, not as two separate one-dimensional sweeps.
        per_elem = lambda preds, targets: (preds - targets).abs()

    if multitask:
        def loss_fn(preds, targets, mask):
            return _masked_mean(per_elem(preds, targets), mask)
    else:
        def loss_fn(preds, targets, mask=None):
            return per_elem(preds, targets).mean()

    return loss_fn


def _decay_param_groups(model: nn.Module, weight_decay: float) -> list[dict]:
    """Split parameters into a decayed and a non-decayed group.

    Weight matrices are decayed; every bias and every normalization parameter is not. Modules are
    walked with ``recurse=False`` so each parameter is classified exactly once by the module that
    owns it, which is what lets a LayerNorm's ``weight`` be recognized as a norm parameter rather
    than as a weight matrix.
    """
    decay: list[nn.Parameter] = []
    no_decay: list[nn.Parameter] = []
    norm_types = (nn.LayerNorm, nn.BatchNorm1d, nn.GroupNorm)
    for module in model.modules():
        is_norm = isinstance(module, norm_types)
        for name, p in module.named_parameters(recurse=False):
            if not p.requires_grad:
                continue
            (no_decay if (is_norm or name.endswith("bias")) else decay).append(p)
    return [
        {"params": decay, "weight_decay": weight_decay},
        {"params": no_decay, "weight_decay": 0.0},
    ]


def _format_per_drug_block(
    per_drug_mse: np.ndarray,
    per_drug_n: np.ndarray,
    drug_names: list[str] | None,
    topk: int,
) -> str:
    if drug_names is None or topk <= 0:
        return ""
    # Heads with zero val support are NaN; surface them but don't rank by them.
    finite = np.isfinite(per_drug_mse)
    if finite.sum() == 0:
        return ""
    order = np.argsort(np.where(finite, per_drug_mse, np.inf))
    best = order[: min(topk, finite.sum())]
    worst = order[::-1][: min(topk, finite.sum())]
    # `per_drug_n` is the column sum of the mask, which is a count only while the mask is 0/1.
    # Under density weighting the mask carries weights, so it is the weight sum -- printing it as
    # `n=` read as a sample size (found 12.08.2026, audit 09). Labelled `w=` when it is not an
    # integer count. The MSE beside it is likewise the weighted per-drug MSE in that case, which is
    # the objective being optimized but is not comparable to an unweighted run.
    def fmt(idx):
        n = per_drug_n[idx]
        amount = f"n={int(n)}" if float(n).is_integer() else f"w={n:.1f}"
        return f"{drug_names[idx]}:{per_drug_mse[idx]:.3f}({amount})"
    return (
        "  best : " + ", ".join(fmt(i) for i in best) + "\n"
        "  worst: " + ", ".join(fmt(i) for i in worst)
    )


def train_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader,
    config: TrainConfig | None = None,
    tag: str = "model",
    drug_names: list[str] | None = None,
) -> tuple[nn.Module, TrainHistory]:
    """Train ``model`` and return (best_model, history).

    The returned model is loaded with the state_dict that achieved the lowest
    val MSE seen during training, not the final-epoch weights.

    Multi-task batches (3-tuples ``(x, y, mask)``) are detected automatically
    from ``train_loader``. In that case the loss is masked-MSE / masked-MAE
    and the printed metric is the overall masked val MSE; the top-k best/worst
    per-drug val MSEs are also printed when ``drug_names`` is provided.
    """
    config = config or TrainConfig()
    set_seed(config.seed)

    device = pick_device()
    print(f"[{tag}] Training on device: {device}")
    model.to(device)

    multitask = _is_multitask_loader(train_loader)
    if multitask:
        print(f"[{tag}] Multi-task mode detected (masked {config.loss.upper()}).")

    loss_fn = _make_loss_fn(config, multitask=multitask)
    # AdamW (Selin, 12.08.2026, review item 10), and `weight_decay` set to 0.0 with it -- see the
    # field's own comment for why those two are one decision rather than two.
    #
    # The argument for the switch: Loshchilov & Hutter, *Decoupled Weight Decay
    # Regularization*, ICLR 2019, show that passing `weight_decay` to Adam adds an L2 term to the
    # gradient which Adam's adaptive scaling then divides through, so decay strength is entangled
    # with gradient history, and that the decoupled form is what "weight decay" normally means. That
    # so the decoupled form is what "weight decay" normally means, and it is now what we use.
    #
    # The migration carried a price, measured before switching
    # (8 epochs, real OncoMLP, real grouping, synthetic data at workload scale):
    #
    #     Adam   wd=1e-3 : weight-matrix L2 norm 5.597
    #     AdamW  wd=1e-3 : weight-matrix L2 norm 8.971
    #     AdamW  wd=0    : weight-matrix L2 norm 8.975   <- no-decay reference
    #
    # AdamW at the SAME nominal value is indistinguishable from no weight decay at all -- 0.04 % from
    # the zero-decay run -- while Adam at 1e-3 shrinks the weight matrices by 38 %. AdamW's step is
    # `theta -= lr*wd*theta`, a factor of 1 - 1e-6 per step at these settings; Adam's L2 enters the
    # gradient and is amplified by 1/sqrt(v). Reproducing Adam's shrinkage under AdamW needs
    # wd ~= 1.0, three orders of magnitude up.
    #
    # So migrating meant re-deriving `weight_decay`, and no sourced value exists for either
    # optimizer. Rather than carry 1e-3 -- which under AdamW does nothing while looking deliberate --
    # the decay is set to 0.0 and the model trains without it.
    #
    # Audit 08's grouping is untouched by any of this: *which* parameters are decayed was settled
    # there, and only *how* decay is applied was ever in question.
    if config.no_decay_bias_and_norm:
        groups = _decay_param_groups(model, config.weight_decay)
        print(f"[{tag}] Weight decay {config.weight_decay:g} on "
              f"{sum(p.numel() for p in groups[0]['params'])} weight parameters; "
              f"{sum(p.numel() for p in groups[1]['params'])} bias/norm parameters exempt.")
        optimizer = optim.AdamW(groups, lr=config.lr)
    else:
        optimizer = optim.AdamW(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=config.scheduler_factor,
        patience=config.scheduler_patience,
    )

    history = TrainHistory()
    best_state = copy.deepcopy(model.state_dict())
    epochs_without_improvement = 0

    for epoch in range(1, config.epochs + 1):
        # -- Training phase --
        model.train()
        running_train_sq_sum = 0.0
        running_train_n = 0.0
        for batch in train_loader:
            if multitask:
                batch_x, batch_y, batch_mask = batch
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                batch_mask = batch_mask.to(device)
            else:
                batch_x, batch_y = batch
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                batch_mask = None

            optimizer.zero_grad()
            preds = model(batch_x)
            loss = loss_fn(preds, batch_y, batch_mask) if multitask else loss_fn(preds, batch_y)
            loss.backward()
            if config.grad_clip is not None:
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=config.grad_clip)
            optimizer.step()

            with torch.no_grad():
                sq = (preds.detach() - batch_y) ** 2
                if multitask:
                    running_train_sq_sum += (sq * batch_mask).sum().item()
                    running_train_n += batch_mask.sum().item()
                else:
                    running_train_sq_sum += sq.sum().item()
                    running_train_n += float(batch_y.numel())

        train_mse = running_train_sq_sum / max(running_train_n, 1.0)

        # -- Validation phase --
        model.eval()
        running_val_sq_sum = 0.0
        running_val_n = 0.0
        per_drug_sq_sum: np.ndarray | None = None
        per_drug_n: np.ndarray | None = None

        with torch.no_grad():
            for batch in val_loader:
                if multitask:
                    batch_x, batch_y, batch_mask = batch
                    batch_x = batch_x.to(device)
                    batch_y = batch_y.to(device)
                    batch_mask = batch_mask.to(device)
                else:
                    batch_x, batch_y = batch
                    batch_x = batch_x.to(device)
                    batch_y = batch_y.to(device)
                    batch_mask = None

                preds = model(batch_x)
                sq = (preds - batch_y) ** 2
                if multitask:
                    masked_sq = sq * batch_mask
                    running_val_sq_sum += masked_sq.sum().item()
                    running_val_n += batch_mask.sum().item()
                    sq_per_drug = masked_sq.sum(dim=0).cpu().numpy()
                    n_per_drug = batch_mask.sum(dim=0).cpu().numpy()
                    if per_drug_sq_sum is None:
                        per_drug_sq_sum = sq_per_drug
                        per_drug_n = n_per_drug
                    else:
                        per_drug_sq_sum += sq_per_drug
                        per_drug_n += n_per_drug
                else:
                    running_val_sq_sum += sq.sum().item()
                    running_val_n += float(batch_y.numel())

        val_mse = running_val_sq_sum / max(running_val_n, 1.0)

        current_lr = optimizer.param_groups[0]["lr"]
        history.train_mse.append(train_mse)
        history.val_mse.append(val_mse)
        history.lr.append(current_lr)

        per_drug_mse_arr = None
        if multitask and per_drug_sq_sum is not None and per_drug_n is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                per_drug_mse_arr = per_drug_sq_sum / np.where(per_drug_n > 0, per_drug_n, np.nan)
            history.per_drug_val_mse.append(per_drug_mse_arr)

        scheduler.step(val_mse)

        improved = val_mse < history.best_val_mse - 1e-6
        if improved:
            history.best_val_mse = val_mse
            history.best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1

        should_log = (
            epoch == 1
            or epoch == config.epochs
            or epoch % config.log_every == 0
            or improved
        )
        if should_log:
            marker = "  <- best" if improved else ""
            print(
                f"[{tag}] Epoch [{epoch:02d}/{config.epochs}] "
                f"| Train MSE: {train_mse:.4f} | Val MSE: {val_mse:.4f} "
                f"| LR: {current_lr:.1e}{marker}"
            )
            if multitask and per_drug_mse_arr is not None and drug_names is not None:
                block = _format_per_drug_block(
                    per_drug_mse_arr,
                    per_drug_n if per_drug_n is not None else np.zeros_like(per_drug_mse_arr),
                    drug_names,
                    config.log_per_drug_topk,
                )
                if block:
                    print(block)

        if epochs_without_improvement >= config.early_stop_patience:
            print(
                f"[{tag}] Early stopping at epoch {epoch} "
                f"(no val improvement for {config.early_stop_patience} epochs)."
            )
            break

    model.load_state_dict(best_state)
    print(
        f"[{tag}] Training complete. "
        f"Best Val MSE: {history.best_val_mse:.4f} at epoch {history.best_epoch}."
    )
    return model, history


# ----------------------------------------------------------------------------
# Run versioning + artifact saving
# ----------------------------------------------------------------------------

DEFAULT_RUNS_ROOT = Path("runs")
INDEX_FILENAME = "runs_index.csv"
INDEX_COLUMNS = [
    "run_id",
    "tag",
    "scope",
    "rep",
    "K",
    "n_train_cells",
    "n_val_cells",
    "best_epoch",
    "best_val_mse",
    "baseline_mean_mse",
    "model_mean_mse",
    "n_beats_baseline",
    "n_total_heads",
    "started_at",
    "finished_at",
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _now_local_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _serialize(obj: Any) -> Any:
    """JSON-serialize configs / metadata, coercing NaN/Inf to None."""
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        return obj if (obj == obj and obj not in (float("inf"), float("-inf"))) else None
    if isinstance(obj, Path):
        return str(obj)
    if is_dataclass(obj):
        return _serialize(asdict(obj))
    if isinstance(obj, Mapping):
        return {str(k): _serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _serialize(obj.tolist())
    if isinstance(obj, np.floating):
        v = obj.item()
        return v if (v == v and v not in (float("inf"), float("-inf"))) else None
    if isinstance(obj, (np.integer, np.bool_)):
        return obj.item()
    return str(obj)


def _maybe_float(x: Any, digits: int = 6) -> Any:
    if x is None:
        return ""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ""
    if np.isnan(v):
        return ""
    return round(v, digits)


def create_run_dir(tag: str, root: Path | str = DEFAULT_RUNS_ROOT) -> Path:
    """Create ``runs/<timestamp>_<tag>/`` and return the path."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    safe_tag = "".join(c if c.isalnum() or c in {"-", "_"} else "_" for c in tag)
    run_id = f"{_now_local_id()}_{safe_tag}"
    run_dir = root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def _write_history_csv(run_dir: Path, history: Any) -> None:
    train_mse = list(getattr(history, "train_mse", []))
    val_mse = list(getattr(history, "val_mse", []))
    lr = list(getattr(history, "lr", []))
    n = max(len(train_mse), len(val_mse), len(lr))
    if n == 0:
        return
    with open(run_dir / "history.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["epoch", "train_mse", "val_mse", "lr"])
        for i in range(n):
            w.writerow(
                [
                    i + 1,
                    _maybe_float(train_mse[i]) if i < len(train_mse) else "",
                    _maybe_float(val_mse[i]) if i < len(val_mse) else "",
                    _maybe_float(lr[i]) if i < len(lr) else "",
                ]
            )


def _write_per_drug_csv(
    run_dir: Path,
    drug_names: Sequence[str],
    model_mse: np.ndarray,
    baseline_mse: np.ndarray | None,
    n_val_per_drug: np.ndarray | None,
) -> None:
    has_baseline = baseline_mse is not None
    has_counts = n_val_per_drug is not None
    with open(run_dir / "per_drug_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        header = ["drug", "model_val_mse"]
        if has_baseline:
            header += ["baseline_val_mse", "delta_model_minus_baseline"]
        if has_counts:
            header += ["n_val"]
        w.writerow(header)
        for i, name in enumerate(drug_names):
            row: list[Any] = [name, _maybe_float(model_mse[i])]
            if has_baseline:
                b = float(baseline_mse[i])  # type: ignore[index]
                delta = (
                    float("nan")
                    if (np.isnan(b) or np.isnan(model_mse[i]))
                    else float(model_mse[i]) - b
                )
                row += [_maybe_float(b), _maybe_float(delta)]
            if has_counts:
                row += [int(n_val_per_drug[i])]  # type: ignore[index]
            w.writerow(row)


def _append_index_row(root: Path, row: dict[str, Any]) -> None:
    path = root / INDEX_FILENAME
    new_file = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=INDEX_COLUMNS)
        if new_file:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in INDEX_COLUMNS})


def save_run(
    *,
    run_dir: Path,
    tag: str,
    config: Any,
    history: Any,
    model: torch.nn.Module,
    run_meta: Mapping[str, Any],
    started_at: str,
    drug_names: Sequence[str] | None = None,
    model_per_drug_val_mse: np.ndarray | None = None,
    baseline_per_drug_val_mse: np.ndarray | None = None,
    n_val_per_drug: np.ndarray | None = None,
) -> dict[str, Any]:
    """Persist all artifacts for one training run. Returns the summary dict."""
    finished_at = utc_now_iso()
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    with open(run_dir / "config.json", "w") as f:
        json.dump(_serialize(config), f, indent=2)

    full_meta = dict(run_meta)
    full_meta.setdefault("started_at", started_at)
    full_meta["finished_at"] = finished_at
    full_meta.setdefault("python_version", sys.version.split()[0])
    full_meta.setdefault("platform", platform.platform())
    full_meta.setdefault("hostname", socket.gethostname())
    full_meta.setdefault("torch_version", torch.__version__)
    if drug_names is not None:
        full_meta["drug_names"] = list(drug_names)
        full_meta.setdefault("K", len(drug_names))
    with open(run_dir / "run_meta.json", "w") as f:
        json.dump(_serialize(full_meta), f, indent=2)

    _write_history_csv(run_dir, history)
    torch.save(model.state_dict(), run_dir / "best_model.pt")

    summary: dict[str, Any] = {
        "tag": tag,
        "run_dir": str(run_dir),
        "best_val_mse": float(getattr(history, "best_val_mse", float("nan"))),
        "best_epoch": int(getattr(history, "best_epoch", -1)),
        "n_epochs_trained": len(getattr(history, "val_mse", [])),
        "final_train_mse": float(history.train_mse[-1]) if getattr(history, "train_mse", []) else None,
        "final_val_mse": float(history.val_mse[-1]) if getattr(history, "val_mse", []) else None,
        "started_at": started_at,
        "finished_at": finished_at,
        "scope": full_meta.get("scope"),
        "rep": full_meta.get("rep"),
        "K": full_meta.get("K", 1),
        "n_train_cells": full_meta.get("n_train_cells"),
        "n_val_cells": full_meta.get("n_val_cells"),
    }

    n_beats: int | None = None
    n_total: int | None = None
    baseline_mean: float | None = None
    model_mean: float | None = None
    if model_per_drug_val_mse is not None and drug_names is not None:
        _write_per_drug_csv(
            run_dir=run_dir,
            drug_names=drug_names,
            model_mse=model_per_drug_val_mse,
            baseline_mse=baseline_per_drug_val_mse,
            n_val_per_drug=n_val_per_drug,
        )
        finite_model = np.isfinite(model_per_drug_val_mse)
        if finite_model.any():
            model_mean = float(np.nanmean(model_per_drug_val_mse))
        if baseline_per_drug_val_mse is not None:
            finite_pair = finite_model & np.isfinite(baseline_per_drug_val_mse)
            n_total = int(finite_pair.sum())
            if n_total > 0:
                n_beats = int(
                    (model_per_drug_val_mse[finite_pair] < baseline_per_drug_val_mse[finite_pair]).sum()
                )
                baseline_mean = float(np.nanmean(baseline_per_drug_val_mse))

    summary["baseline_mean_mse"] = baseline_mean
    summary["model_mean_mse"] = model_mean
    summary["n_beats_baseline"] = n_beats
    summary["n_total_heads"] = n_total

    with open(run_dir / "summary.json", "w") as f:
        json.dump(_serialize(summary), f, indent=2)

    _append_index_row(
        root=run_dir.parent,
        row={
            "run_id": run_dir.name,
            "tag": tag,
            "scope": full_meta.get("scope", ""),
            "rep": full_meta.get("rep", ""),
            "K": full_meta.get("K", 1),
            "n_train_cells": full_meta.get("n_train_cells", ""),
            "n_val_cells": full_meta.get("n_val_cells", ""),
            "best_epoch": summary["best_epoch"],
            "best_val_mse": summary["best_val_mse"],
            "baseline_mean_mse": baseline_mean if baseline_mean is not None else "",
            "model_mean_mse": model_mean if model_mean is not None else "",
            "n_beats_baseline": n_beats if n_beats is not None else "",
            "n_total_heads": n_total if n_total is not None else "",
            "started_at": started_at,
            "finished_at": finished_at,
        },
    )

    print(f"\n[{tag}] Saved run artifacts to: {run_dir}")
    return summary
