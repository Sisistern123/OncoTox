"""Cross-validated out-of-fold predictions -- the one implementation every caller should use.

``cv_evaluate`` in :mod:`scripts.training.train_multitask` answers "is the difference real?" and
returns per-fold *metrics*. This module returns the *predictions*, which is what any analysis
downstream of training needs: per-drug correlations, the cell-line-effect normalization, calibration
diagnostics. Those were previously re-implemented in each notebook that wanted them, with small
divergences (line-level versus cell-level aggregation, different epoch counts, bias initialized or
not), which is exactly how two runs stop being comparable.

The protocol is the project standard and does not vary between callers: K-fold ``GroupKFold`` over
cell lines, restricted to the ``split_ctrp`` train+val lines so the fixed test set is never touched,
and every held-out line predicted by a model that has not seen it. Per-cell predictions are averaged
back to one value per line by :func:`line_level_predictions`, because the label is per line and
scoring per cell would score pseudo-replicates.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import GroupKFold
from torch.utils.data import DataLoader

from scripts.model.dataset import MultiDrugDataset
from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS, OncoMLP
from scripts.training.density_weighting import (
    DEFAULT_ALPHA,
    DEFAULT_CAP,
    fit_weight_fns,
    line_level,
    weight_matrix,
)
from scripts.training.training_utils import TrainConfig, train_model


def grouped_folds(
    adata,
    *,
    n_splits: int = 5,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
) -> tuple[np.ndarray, list[tuple[np.ndarray, np.ndarray]]]:
    """The project's fold partition: ``(idx, folds)`` over the eligible cells.

    Exposed separately so that anything scored alongside the MLP -- the ridge control above all --
    uses the *same* partition rather than an independently constructed one that happens to share a
    seed. ``idx`` indexes into the full cell array; each fold is a pair of positions into ``idx``.
    """
    groups = adata.obs[group_col].astype(str).to_numpy()
    eligible = adata.obs["split_ctrp"].isin(eligible_splits).to_numpy()
    idx = np.flatnonzero(eligible)
    return idx, list(GroupKFold(n_splits=n_splits).split(idx, groups=groups[idx]))


def oof_predictions(
    adata,
    rep: str,
    drugs: list[str],
    *,
    config: TrainConfig | None = None,
    hidden_dims: tuple[int, ...] | None = None,
    n_splits: int = 5,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
    batch_size: int = 128,
    dropout: float = 0.5,
    input_dropout: float = 0.1,
    density_weighting: bool = False,
    alpha: float = DEFAULT_ALPHA,
    cap: float = DEFAULT_CAP,
    init_head_bias: bool = True,
    tag: str = "oof",
) -> tuple[np.ndarray, list[dict]]:
    """Fit ``n_splits`` cell-line-grouped folds and return out-of-fold per-cell predictions.

    :param density_weighting: weight each observation by inverse label density
        (:mod:`scripts.training.density_weighting`). The density is fitted **inside each fold, on the
        training lines only** -- it is a function of the labels, so fitting it once over all lines
        would let held-out labels inform training. The weights are carried in the mask tensor, which
        turns the masked mean into a weighted one without touching the loss code.
    :param init_head_bias: initialize each head's bias to that drug's mean over the fold's training
        *lines*, so the model starts at the null predictor. Necessary on an uncentred target such as
        raw AUC, where the bias must reach ~0.7 from a default near 0; harmless on a centred one.
    :returns: ``(pred, folds)`` -- ``pred`` is (n_cells, len(drugs)) with NaN where a cell was never
        held out, ``folds`` is a per-fold log.
    """
    config = config or TrainConfig()
    hidden_dims = tuple(hidden_dims or DEFAULT_HIDDEN_DIMS[rep])

    groups = adata.obs[group_col].astype(str).to_numpy()
    idx, fold_split = grouped_folds(
        adata, n_splits=n_splits, group_col=group_col, eligible_splits=eligible_splits
    )

    all_drugs = list(adata.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in drugs]
    Y = np.asarray(adata.obsm["Y_ctrp"], dtype=np.float32)[:, kcol]
    M = np.asarray(adata.obsm["M_ctrp"], dtype=bool)[:, kcol]

    pred = np.full((adata.n_obs, len(drugs)), np.nan, dtype=float)
    folds: list[dict] = []
    for fold, (tr, va) in enumerate(fold_split, start=1):
        trc = np.zeros(adata.n_obs, dtype=bool)
        trc[idx[tr]] = True
        vac = np.zeros(adata.n_obs, dtype=bool)
        vac[idx[va]] = True

        tr_ds = MultiDrugDataset(adata=adata, use_rep=rep, cell_mask=trc, drugs=drugs)
        va_ds = MultiDrugDataset(adata=adata, use_rep=rep, cell_mask=vac, drugs=drugs)

        # per-drug statistics from this fold's training LINES (never cells, never held-out lines)
        tr_lines = np.unique(groups[trc])
        y_lines, obs_lines = line_level(Y, M, groups, tr_lines)

        if density_weighting:
            fns = fit_weight_fns(y_lines, obs_lines, alpha=alpha, cap=cap)
            tr_ds.mask = torch.from_numpy(weight_matrix(fns, Y[trc], M[trc]))
            va_ds.mask = torch.from_numpy(weight_matrix(fns, Y[vac], M[vac]))

        model = OncoMLP(
            input_dim=tr_ds.X.shape[1],
            hidden_dims=hidden_dims,
            dropout_rate=dropout,
            input_dropout=input_dropout,
            norm="layer",
            output_dim=len(drugs),
        )
        if init_head_bias:
            means = np.array(
                [y_lines[obs_lines[:, j], j].mean() if obs_lines[:, j].any() else 0.0
                 for j in range(len(drugs))], dtype=np.float32)
            head = [m for m in model.modules() if isinstance(m, torch.nn.Linear)][-1]
            with torch.no_grad():
                head.bias.copy_(torch.from_numpy(means))

        # Explicit generator: without one, RandomSampler seeds itself from the global torch
        # RNG at iteration time, which happens to be seeded because train_model calls
        # set_seed(config.seed) first -- an ordering dependency that breaks silently.
        best, hist = train_model(
            model,
            DataLoader(
                tr_ds,
                batch_size=batch_size,
                shuffle=True,
                generator=torch.Generator().manual_seed(config.seed),
            ),
            DataLoader(va_ds, batch_size=batch_size, shuffle=False),
            config=config,
            tag=f"{tag}_f{fold}",
            drug_names=drugs,
        )
        best = best.to("cpu").eval()  # train_model leaves it on mps/cuda
        with torch.no_grad():
            x = torch.from_numpy(np.asarray(adata.obsm[rep], dtype=np.float32)[vac])
            pred[vac] = best(x).numpy()  # held-out lines only

        folds.append({
            "fold": fold, "rep": rep, "n_train_lines": int(len(tr_lines)),
            "n_val_lines": int(np.unique(groups[vac]).size),
            "best_epoch": hist.best_epoch, "best_val_obj": float(hist.best_val_mse),
        })
    return pred, folds


def line_level_predictions(
    pred: np.ndarray,
    adata,
    drugs: list[str],
    *,
    truth: np.ndarray | None = None,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
    **extra,
) -> pd.DataFrame:
    """Average per-cell predictions to one row per (drug, cell line) and attach the truth.

    :param truth: (n_cells, len(drugs)) labels to score against; defaults to the ``Y_ctrp`` columns
        of ``drugs``. Pass an explicit array to score every model against one common yardstick even
        when they trained on different scores.
    :param extra: constant columns to stamp on every row (e.g. ``rep=..., weighted=...``), so
        several configurations can be concatenated into one tidy table.
    :returns: DataFrame(drug, cell_line, y_true, y_pred, **extra) -- the shape every downstream
        analysis expects, including the cell-line-effect normalization.
    """
    groups = adata.obs[group_col].astype(str).to_numpy()
    eligible = adata.obs["split_ctrp"].isin(eligible_splits).to_numpy()
    all_drugs = list(adata.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in drugs]
    M = np.asarray(adata.obsm["M_ctrp"], dtype=bool)[:, kcol]
    Y = np.asarray(adata.obsm["Y_ctrp"], dtype=float)[:, kcol] if truth is None else np.asarray(truth)

    rows = []
    for j, d in enumerate(drugs):
        for ln in np.unique(groups[eligible]):
            ci = np.flatnonzero((groups == ln) & eligible)
            obs = ci[M[ci, j] & np.isfinite(pred[ci, j])]
            if obs.size:
                rows.append({"drug": d, "cell_line": ln,
                             "y_true": float(Y[obs, j].mean()),
                             "y_pred": float(pred[obs, j].mean()), **extra})
    return pd.DataFrame(rows)
