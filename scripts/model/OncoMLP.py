from __future__ import annotations

from collections.abc import Sequence

import torch.nn as nn


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
