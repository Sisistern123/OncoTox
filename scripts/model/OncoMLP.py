from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import torch
import torch.nn as nn

# Matched trunk for a fair PCA-vs-scGPT comparison: both reps use the same hidden layers, so only the
# input representation (and its first projection) differs. Defined here rather than in the training
# module so that scripts.training.cv can read it without importing train_multitask, which imports cv.
DEFAULT_HIDDEN_DIMS = {
    "X_pca": (128, 64),
    "X_scGPT": (128, 64),
}


class OncoMLP(nn.Module):
    """Small regression MLP for cell-level viability prediction.

    Designed to be robust across input regimes:
      * PCA baseline (~50 dims) and scGPT embeddings (512 dims)
      * Small cell-line-grouped batches where BatchNorm running stats can be noisy

    Differences from the original 256->64 redesign:
      * LayerNorm by default (more stable than BatchNorm for embedding regression
        across heterogeneous cell-line batches).
      * GELU activation (smoother than ReLU for continuous targets).
      * Optional input dropout to regularize the raw embedding directly.
      * Configurable hidden_dims so PCA / scGPT can be tuned independently if needed.
    """

    def __init__(
        self,
        input_dim: int,
        hidden_dims: Sequence[int] = (64, 32),
        dropout_rate: float = 0.5,
        input_dropout: float = 0.1,
        norm: str = "layer",
        output_dim: int = 1,
    ):
        super().__init__()

        # "batch" and "none" look dead -- every construction in scripts/ passes "layer" -- but they
        # are not. `norm` is persisted to run_meta.json by save_run and read back verbatim when
        # `4_training` reconstructs a saved model (`OncoMLP(..., norm=m['norm'])`), so any run that
        # recorded another value still needs its branch to reload. Checked 12.08.2026 (Selin):
        # keep them, and say why, rather than narrow the constructor to what today's callers happen
        # to pass. LayerNorm remains the only value any current path uses, for the reason in the
        # class docstring.
        if norm not in {"layer", "batch", "none"}:
            raise ValueError(f"norm must be 'layer', 'batch', or 'none' (got {norm!r})")
        if output_dim < 1:
            raise ValueError(f"output_dim must be >= 1 (got {output_dim})")

        self.output_dim = output_dim
        layers: list[nn.Module] = []

        if input_dropout and input_dropout > 0:
            layers.append(nn.Dropout(input_dropout))

        prev_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev_dim, h))
            if norm == "layer":
                layers.append(nn.LayerNorm(h))
            elif norm == "batch":
                layers.append(nn.BatchNorm1d(h))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout_rate))
            prev_dim = h

        layers.append(nn.Linear(prev_dim, output_dim))
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)


def init_head_bias_(model: nn.Module, means: np.ndarray) -> None:
    """Start each head at its drug's mean instead of at ~0. In place.

    ``nn.Linear`` initializes biases uniformly in ±1/sqrt(fan_in), so an untouched model predicts
    near 0 while ``auc_cc`` sits near 0.9 -- the first epochs then go on climbing to the drug means
    rather than on the differences between cell lines, and early stopping here regularly selects
    epoch 1-3. Initializing the final bias at the base rate is standard practice for exactly this
    reason (Lin et al., *Focal Loss for Dense Object Detection*, ICCV 2017, §4.1 "prior" init, done
    there for a class prior rather than a regression mean).

    ``means`` must be one value per head, computed on the *fitting* lines of the fold and nowhere
    else -- it is a statistic of the labels, so a mean over held-out lines would inform training.

    ⚠️ **This does not make the model start at the null predictor, only near its level.** The head's
    weight rows are still randomly initialized (mean row norm ~0.58 over a LayerNorm'd 64-d hidden
    vector), so at initialization the predictions already scatter with a standard deviation of ~0.31
    across cells -- against a true across-line spread of order 0.17 on ``auc_cc``. Measured on
    synthetic unit-variance input at seeds 42 and 0; the mean prediction lands at 0.76 and 0.96 for
    a requested 0.90. Starting genuinely *at* the null would also require shrinking the output
    layer's weight initialization, which Lin et al. do (sigma=0.01 on the final layer) and this
    project has not decided. Open, audit 08.

    Applied 12.08.2026 to all three training paths. Before that only ``cv.oof_predictions`` did it,
    so the fixed-split and 8-run-matrix paths trained against an offset the panel run did not:
    docs/steps/03, *The uncentred target is handled the same way in every training path*.
    """
    head = [m for m in model.modules() if isinstance(m, nn.Linear)][-1]
    means = np.asarray(means, dtype=np.float32)
    if means.shape != head.bias.shape:
        raise ValueError(
            f"means has shape {means.shape}, but the output layer has {tuple(head.bias.shape)} "
            f"heads -- one mean per drug is required, in the drug column order."
        )
    with torch.no_grad():
        head.bias.copy_(torch.from_numpy(means))
