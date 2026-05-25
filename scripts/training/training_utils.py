"""Shared training utilities for OncoTox MLP heads.

Adds the pieces missing from the original training scripts:
  * Reproducible seeding (Python / NumPy / PyTorch / CUDA / MPS).
  * LR scheduling via ReduceLROnPlateau.
  * Early stopping with patience.
  * Best-val-MSE checkpointing in-memory (so the returned model is the best
    seen, not the noisy final-epoch one).
  * Gradient clipping to dampen the val-MSE oscillations we saw on scGPT
    (val MSE bounced between 0.0354 and 0.0423 across epochs with the
    original setup).
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


@dataclass
class TrainHistory:
    train_mse: list[float] = field(default_factory=list)
    val_mse: list[float] = field(default_factory=list)
    lr: list[float] = field(default_factory=list)
    best_val_mse: float = float("inf")
    best_epoch: int = -1


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


def _make_loss(name: str) -> nn.Module:
    name = name.lower()
    if name == "mse":
        return nn.MSELoss()
    if name == "huber":
        # Huber/SmoothL1 reduces sensitivity to viability outliers near 0 / >1.
        return nn.SmoothL1Loss(beta=0.05)
    raise ValueError(f"Unknown loss: {name!r}")


def train_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader,
    config: TrainConfig | None = None,
    tag: str = "model",
) -> tuple[nn.Module, TrainHistory]:
    """Train ``model`` and return (best_model, history).

    The returned model is loaded with the state_dict that achieved the lowest
    val MSE seen during training, not the final-epoch weights.
    """
    config = config or TrainConfig()
    set_seed(config.seed)

    device = pick_device()
    print(f"[{tag}] Training on device: {device}")
    model.to(device)

    criterion = _make_loss(config.loss)
    optimizer = optim.Adam(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=config.scheduler_factor,
        patience=config.scheduler_patience,
    )

    # Always report MSE (even when training with Huber loss) so runs are comparable.
    mse_metric = nn.MSELoss(reduction="sum")

    history = TrainHistory()
    best_state = copy.deepcopy(model.state_dict())
    epochs_without_improvement = 0

    for epoch in range(1, config.epochs + 1):
        # -- Training phase --
        model.train()
        running_train_mse_sum = 0.0
        n_train = 0
        for batch_x, batch_y in train_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)

            optimizer.zero_grad()
            preds = model(batch_x)
            loss = criterion(preds, batch_y)
            loss.backward()
            if config.grad_clip is not None:
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=config.grad_clip)
            optimizer.step()

            with torch.no_grad():
                running_train_mse_sum += mse_metric(preds, batch_y).item()
            n_train += batch_x.size(0)

        train_mse = running_train_mse_sum / max(n_train, 1)

        # -- Validation phase --
        model.eval()
        running_val_mse_sum = 0.0
        n_val = 0
        with torch.no_grad():
            for batch_x, batch_y in val_loader:
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                preds = model(batch_x)
                running_val_mse_sum += mse_metric(preds, batch_y).item()
                n_val += batch_x.size(0)
        val_mse = running_val_mse_sum / max(n_val, 1)

        current_lr = optimizer.param_groups[0]["lr"]
        history.train_mse.append(train_mse)
        history.val_mse.append(val_mse)
        history.lr.append(current_lr)

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
