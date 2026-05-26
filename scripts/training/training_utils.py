"""Shared training utilities for OncoTox MLP heads.

  * ``TrainConfig`` / ``TrainHistory`` dataclasses.
  * ``train_model(...)`` -- single-task + multi-task masked-loss in one fn.
  * ``set_seed``, ``pick_device`` helpers.

Single-task training keeps the original (preds, targets) -> MSE/Huber path.
When the data loader yields 3-tensor batches ``(x, y, mask)``, training
auto-switches to masked MSE / masked Huber so missing (cell, drug) entries
are skipped; per-drug val MSE is tracked alongside.
"""

from __future__ import annotations

import copy
import random
from dataclasses import dataclass, field

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader


@dataclass
class TrainConfig:
    epochs: int = 50
    lr: float = 1e-3
    weight_decay: float = 1e-3
    grad_clip: float | None = 1.0
    scheduler_patience: int = 3
    scheduler_factor: float = 0.5
    early_stop_patience: int = 10
    log_every: int = 5
    seed: int = 42
    loss: str = "mse"  # "mse" or "huber"
    huber_beta: float = 0.05
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

    Single-task callers can pass ``mask=None`` and get plain MSE / Huber.
    Multi-task callers must pass ``mask`` (same shape as ``targets``).
    """
    name = config.loss.lower()
    if name not in {"mse", "huber"}:
        raise ValueError(f"Unknown loss: {config.loss!r}")

    if name == "mse":
        per_elem = lambda preds, targets: (preds - targets) ** 2
    else:
        per_elem = lambda preds, targets: nn.functional.smooth_l1_loss(
            preds, targets, beta=config.huber_beta, reduction="none"
        )

    if multitask:
        def loss_fn(preds, targets, mask):
            return _masked_mean(per_elem(preds, targets), mask)
    else:
        def loss_fn(preds, targets, mask=None):
            return per_elem(preds, targets).mean()

    return loss_fn


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
    def fmt(idx):
        return f"{drug_names[idx]}:{per_drug_mse[idx]:.3f}(n={int(per_drug_n[idx])})"
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
    from ``train_loader``. In that case the loss is masked-MSE / masked-Huber
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
    optimizer = optim.Adam(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
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
