"""Multiple-instance training: one bag is one cell line, and the bag is complete.

The per-cell path in :mod:`scripts.training.cv` draws batches of *cells* across lines and asks each
cell for its line's label. This module draws batches of *bags*: every sequenced cell of one line goes
through the encoder, the per-cell predictions are averaged, and that mean is compared against the
line's single label. The model, the trunk, the optimizer settings and the fold partition are the ones
``cv.py`` uses -- **only the unit of the loss changes**, which is what makes the two runs comparable
(docs/TODO.md, the governing rule).

Why that change is the whole point, in one line. For a line with per-cell predictions ``p_c`` and
label ``y``::

    (1/n) Σ_c (p_c − y)²  =  (p̄ − y)²  +  Var_c(p)

an identity, not an approximation. The left side is what ``cv.oof_predictions`` minimizes, regrouped
by line; the first term on the right is what this module minimizes. The two objectives differ by
exactly the within-line variance of the predictions, so the per-cell model charges for the quantity
Q2 asks about and this one does not
(docs/steps/03-model-and-training-design.md, *The penalty on within-line variation is exact*).

**Bags are complete, and that is not a batching convenience.** For a sub-bag of ``B`` cells drawn from
a line of ``n`` the expected loss is ``(p̄_n − y)² + (σ²/B)(1 − (B−1)/(n−1))``: at ``B = 1`` it is
exactly the per-cell loss, at ``B = n`` the variance term is gone. Bag size therefore dials
continuously between the two architectures and sets how hard the objective charges for the very
quantity under study. Fixed-size sub-bags would also restore the depth weighting that full-line bags
remove, since a deeply sequenced line would yield proportionally more bags. Decided by Selin,
13.08.2026; the criterion is ``notebooks/4b_mil_training.ipynb`` §2.

The costs, taken deliberately: one gradient step per line, so an epoch here is ~122 steps against the
per-cell path's ~330 at ``batch_size=128``, and peak memory scales with the largest line (1,990
cells).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np
import torch
import torch.optim as optim

from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS, OncoMLP, init_head_bias_
from scripts.training.cv import (
    CV_FOLD_PCA_KEY,
    fold_pca_projections_for,
    grouped_folds,
    inner_holdout,
    per_drug_line_mean,
)
from scripts.training.density_weighting import (
    DEFAULT_ALPHA,
    DEFAULT_CAP,
    WeightFn,
    fit_weight_fns,
    line_level,
    weight_matrix,
)
# `_decay_param_groups` and `_make_loss_fn` are private to training_utils and are imported rather
# than reimplemented, for the reason cv.py gives for importing `_pca_fitted_on_train`: which
# parameters decay and what the per-element error is were *decided* (audit 08; item 9A), and a second
# implementation is how the bag path and the per-cell path silently stop being arms of one
# experiment. Only the loop below is new -- the model, the optimizer, the grouping, the loss and the
# seeding are the same objects 4a uses.
from scripts.training.training_utils import (
    TrainConfig,
    TrainHistory,
    _decay_param_groups,
    _make_loss_fn,
    pick_device,
    set_seed,
)


@dataclass
class Bags:
    """One fold's cells, grouped into complete per-line bags.

    Cells are stored contiguously per bag so a bag is a slice rather than a gather -- the loop below
    touches every bag once per epoch, and fancy-indexing 1,990 rows out of a 40,000-row tensor each
    time is the one avoidable cost in it.

    :ivar X: (n_cells, D) representation, in bag order.
    :ivar y: (n_bags, K) one label row per bag. The label is per (cell line, drug), so this is a
        lookup, not an average.
    :ivar w: (n_bags, K) the masked loss's denominator -- 0/1 observation mask, or the density
        weights standing in for it, exactly as the per-cell path carries them in the mask tensor.
    :ivar bounds: (n_bags + 1,) start/stop offsets into ``X``.
    :ivar lines: (n_bags,) cell-line name per bag, in bag order.
    :ivar cells: (n_cells,) positions into the *full* cell array, in bag order, so a prediction can
        be written back to the row it came from.
    """

    X: torch.Tensor
    y: torch.Tensor
    w: torch.Tensor
    bounds: np.ndarray
    lines: np.ndarray
    cells: np.ndarray

    def __len__(self) -> int:
        return len(self.lines)

    def bag(self, i: int) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """``(X_bag, y_bag, w_bag)`` for bag ``i`` -- a view on ``X``, not a copy."""
        lo, hi = self.bounds[i], self.bounds[i + 1]
        return self.X[lo:hi], self.y[i], self.w[i]


def build_bags(
    X_all: np.ndarray,
    Y: np.ndarray,
    M: np.ndarray,
    groups: np.ndarray,
    cell_mask: np.ndarray,
    *,
    universe: np.ndarray,
    weight_fns: list[WeightFn] | None = None,
) -> Bags:
    """Group the cells selected by ``cell_mask`` into one complete bag per cell line.

    ``cell_mask`` comes from the same ``cv.inner_holdout`` split the per-cell path uses, which splits
    *lines* rather than cells -- so a line is wholly inside the mask or wholly outside it, and every
    bag built here is complete with respect to the sequenced cells of its line.

    **Both properties the bag loss rests on are checked here rather than trusted**, because each
    fails silently and produces plausible numbers. A mask that cut a line in half would build
    sub-bags, and the variance term in the module docstring would return at a weight nobody chose. A
    label that varied within a line would break the regrouping the identity performs.

    :param universe: the cells a bag could be drawn from -- the run's eligible set, i.e. the
        ``split_ctrp`` train+val cells. Completeness is defined against it and cannot be defined
        without it: a line absent from ``cell_mask`` is a held-out line, which is correct, while a
        line *partly* present is a sub-bag, which is not.
    :param weight_fns: per-drug density weights fitted on this fold's fitting lines
        (:func:`density_weighting.fit_weight_fns`). Applied to the **line-level** labels, which is
        where they were always fitted; ``None`` leaves the plain 0/1 mask.
    """
    cell_mask = np.asarray(cell_mask, dtype=bool)
    pos = np.flatnonzero(cell_mask)
    g = groups[pos]
    order = np.argsort(g, kind="stable")  # stable so cells keep their original order inside a bag
    cells = pos[order]
    g_sorted = g[order]

    lines, starts = np.unique(g_sorted, return_index=True)
    bounds = np.append(starts, len(cells))

    y_lines, obs_lines = line_level(Y, M, groups, lines)
    universe = np.asarray(universe, dtype=bool)
    for i, ln in enumerate(lines):
        sel = cells[bounds[i]:bounds[i + 1]]
        n_available = int(((groups == ln) & universe).sum())
        if sel.size != n_available:
            raise ValueError(
                f"cell line {str(ln)!r} contributes {sel.size} of its {n_available} eligible cells, so "
                f"this bag is a sub-bag. Full-line bags are what delete the within-line variance "
                f"term from the objective (Selin, 13.08.2026); a partial bag restores it at a "
                f"weight nobody chose. Split by line, not by cell."
            )
        obs = obs_lines[i]
        if obs.any() and not np.allclose(Y[sel][:, obs], y_lines[i, obs], equal_nan=True):
            raise ValueError(
                f"cell line {str(ln)!r} carries more than one label value across its cells; the bag loss "
                f"and the per-cell loss are only comparable while the label is per line."
            )

    W = (weight_matrix(weight_fns, y_lines, obs_lines) if weight_fns is not None
         else obs_lines.astype(np.float32))
    return Bags(
        X=torch.from_numpy(np.asarray(X_all, dtype=np.float32)[cells]),
        y=torch.from_numpy(np.nan_to_num(y_lines, nan=0.0).astype(np.float32)),
        w=torch.from_numpy(np.asarray(W, dtype=np.float32)),
        bounds=bounds,
        lines=lines,
        cells=cells,
    )


def bag_predictions(model, bags: Bags, device) -> tuple[torch.Tensor, torch.Tensor]:
    """``(per_cell, per_bag)`` -- (n_cells, K) predictions and their (n_bags, K) bag means.

    The bag prediction is the **mean** of its cells' predicted responses, matching the definition of
    the label it is compared against (Selin, 13.08.2026: mean pooling, fixed before the model
    existed and not revisited). Instance-level, so the per-cell tensor is the model's native output
    and the bag tensor is derived from it -- not the other way round, which is what makes every
    downstream stage of the criterion computable at all.
    """
    per_cell = model(bags.X.to(device))
    sums = torch.zeros(len(bags), per_cell.shape[1], device=device, dtype=per_cell.dtype)
    idx = torch.from_numpy(
        np.repeat(np.arange(len(bags)), np.diff(bags.bounds))
    ).to(device)
    sums.index_add_(0, idx, per_cell)
    counts = torch.from_numpy(np.diff(bags.bounds).astype(np.float32)).to(device).unsqueeze(1)
    return per_cell, sums / counts


def train_bag_model(
    model,
    train_bags: Bags,
    stop_bags: Bags,
    config: TrainConfig | None = None,
    tag: str = "mil",
) -> tuple[object, TrainHistory]:
    """Train on complete per-line bags; return the best checkpoint and its history.

    One bag is one example and **one gradient step**, so lines are weighted equally regardless of how
    deeply each was sequenced -- which is the resolution of the per-cell path's depth-weighting defect
    (docs/steps/03, *the loss weights cell lines by how deeply they were sequenced*) rather than an
    incidental property of this loop. It costs an epoch of ~120 steps against the per-cell path's
    ~330, at the same learning rate and epoch cap, because the configuration is 4a's and the
    architecture is the only thing permitted to move.

    Early stopping watches the bag MSE on ``stop_bags``, which under ``cv.inner_holdout`` are lines
    the model never fits on. The per-cell path's early-stopping quantity is a per-*cell* MSE; here it
    is the objective actually being minimized, which is the bag MSE. That difference is forced by the
    architecture, not chosen.

    ⚠️ **Dropout does not act the same way here.** It is applied per cell and independently, so
    averaging a 1,990-cell bag cancels far more of it than averaging a 56-cell one -- the effective
    regularization on the bag prediction falls with bag size. The rates are 4a's and are deliberately
    not retuned (one change at a time), so this is recorded as a known asymmetry rather than fixed.
    """
    config = config or TrainConfig()
    set_seed(config.seed)
    device = pick_device()
    print(f"[{tag}] Training on device: {device} | {len(train_bags)} train bags, "
          f"{len(stop_bags)} early-stopping bags")
    model.to(device)

    loss_fn = _make_loss_fn(config, multitask=True)
    groups_ = _decay_param_groups(model, config.weight_decay)
    optimizer = optim.AdamW(groups_, lr=config.lr)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=config.scheduler_factor, patience=config.scheduler_patience,
    )

    history = TrainHistory()
    best_state = copy.deepcopy(model.state_dict())
    stalled = 0
    # Explicit generator, for the reason cv.py gives: without one the shuffle seeds itself from the
    # global torch RNG at iteration time, which happens to be seeded because set_seed ran above -- an
    # ordering dependency that breaks silently.
    gen = torch.Generator().manual_seed(config.seed)

    for epoch in range(1, config.epochs + 1):
        model.train()
        sq_sum = w_sum = 0.0
        for i in torch.randperm(len(train_bags), generator=gen).tolist():
            x, y, w = train_bags.bag(i)
            x, y, w = x.to(device), y.to(device).unsqueeze(0), w.to(device).unsqueeze(0)
            optimizer.zero_grad()
            pooled = model(x).mean(dim=0, keepdim=True)
            loss = loss_fn(pooled, y, w)
            loss.backward()
            if config.grad_clip is not None:
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=config.grad_clip)
            optimizer.step()
            with torch.no_grad():
                sq = (pooled.detach() - y) ** 2
                sq_sum += float((sq * w).sum())
                w_sum += float(w.sum())
        train_mse = sq_sum / max(w_sum, 1.0)

        model.eval()
        with torch.no_grad():
            _, pooled = bag_predictions(model, stop_bags, device)
            w = stop_bags.w.to(device)
            sq = (pooled - stop_bags.y.to(device)) ** 2
            val_mse = float((sq * w).sum() / w.sum().clamp_min(1.0))
            per_drug_w = w.sum(0)
            history.per_drug_val_mse.append(
                torch.where(per_drug_w > 0, (sq * w).sum(0) / per_drug_w.clamp_min(1.0),
                            torch.full_like(per_drug_w, float("nan"))).cpu().numpy()
            )

        history.train_mse.append(train_mse)
        history.val_mse.append(val_mse)
        history.lr.append(optimizer.param_groups[0]["lr"])
        scheduler.step(val_mse)

        if val_mse < history.best_val_mse:
            history.best_val_mse, history.best_epoch = val_mse, epoch
            best_state = copy.deepcopy(model.state_dict())
            stalled = 0
        else:
            stalled += 1
        if epoch % config.log_every == 0 or epoch == 1:
            print(f"[{tag}] epoch {epoch:3d} | train bag MSE {train_mse:.5f} | "
                  f"val bag MSE {val_mse:.5f} | lr {history.lr[-1]:.2e}")
        if stalled >= config.early_stop_patience:
            print(f"[{tag}] early stop at epoch {epoch} (best {history.best_epoch}, "
                  f"{history.best_val_mse:.5f})")
            break

    model.load_state_dict(best_state)
    return model, history


def bag_oof_predictions(
    adata,
    rep: str,
    drugs: list[str],
    *,
    config: TrainConfig | None = None,
    hidden_dims: tuple[int, ...] | None = None,
    n_splits: int = 5,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
    dropout: float = 0.5,
    input_dropout: float = 0.1,
    density_weighting: bool = False,
    alpha: float = DEFAULT_ALPHA,
    cap: float = DEFAULT_CAP,
    init_head_bias: bool = True,
    counts_h5ad=None,
    pca_seed: int | None = None,
    tag: str = "mil",
) -> tuple[np.ndarray, list[dict]]:
    """The bag-model twin of :func:`cv.oof_predictions`, returning the same object it does.

    ``(pred, folds)`` with ``pred`` of shape ``(n_cells, len(drugs))`` and NaN where a cell was never
    held out -- **the identical contract**, so ``cv.line_level_predictions`` and 4a's within-line
    spread table read this array unchanged. That is what "same scorer as the per-cell model"
    (4b §1) means in code rather than in convention.

    **The fold partition is not re-derived here; it is imported.** ``grouped_folds`` and
    ``inner_holdout`` are called with the same arguments 4a passes and are both deterministic, so
    fold *f* holds out the same cell lines in both notebooks and every seed. Stage 1 compares this
    model's within-line spread against 4a's on the same lines, drugs and folds, and a second fold
    implementation here is precisely how that premise would quietly stop holding.

    Everything fitted per fold -- the density weights, the head bias, the PCA rotation -- is fitted on
    the **fitting** lines alone, as in the per-cell path: on ``fitc``, not on ``trc``, so the
    early-stopping lines shape nothing they then judge.

    :param pca_seed: the seed for the per-fold PCA fit, held **separate from the model seed**.
        ``cv.oof_predictions`` passes ``config.seed`` to both, which is harmless there because it
        trains one seed. Here stage 2 asks whether independent initializations rank the same cells
        alike, and if the representation moved with the seed as well, the comparison would mix
        model-initialization variation with PCA-fit variation and could not tell them apart. Defaults
        to ``config.seed`` -- the per-cell path's behaviour, so nothing changes silently for a caller
        that does not set it. **4b passes 42 for all three model seeds (Selin, 13.08.2026)**, so its
        stage 2 measures the model given a fixed representation rather than the pipeline as a whole;
        the cost is that the seeds are then independent draws of the model and not of the PCA fit.
    """
    config = config or TrainConfig()
    hidden_dims = tuple(hidden_dims or DEFAULT_HIDDEN_DIMS[rep])

    groups = adata.obs[group_col].astype(str).to_numpy()
    idx, fold_split = grouped_folds(
        adata, n_splits=n_splits, group_col=group_col, eligible_splits=eligible_splits
    )
    eligible = np.zeros(adata.n_obs, dtype=bool)
    eligible[idx] = True

    all_drugs = list(adata.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in drugs]
    Y = np.asarray(adata.obsm["Y_ctrp"], dtype=np.float32)[:, kcol]
    M = np.asarray(adata.obsm["M_ctrp"], dtype=bool)[:, kcol]

    masks = []
    for tr, va in fold_split:
        trc = np.zeros(adata.n_obs, dtype=bool)
        trc[idx[tr]] = True
        vac = np.zeros(adata.n_obs, dtype=bool)
        vac[idx[va]] = True
        fitc, stopc = inner_holdout(groups, trc)
        masks.append((trc, fitc, stopc, vac))

    projections = fold_pca_projections_for(
        rep, counts_h5ad, [m[1] for m in masks],
        seed=config.seed if pca_seed is None else pca_seed,
    )

    device = pick_device()
    pred = np.full((adata.n_obs, len(drugs)), np.nan, dtype=float)
    folds: list[dict] = []
    for fold, (trc, fitc, stopc, vac) in enumerate(masks, start=1):
        rep_key = rep
        if projections is not None:
            rep_key = CV_FOLD_PCA_KEY
            adata.obsm[rep_key] = projections[fold - 1]
        X_all = np.asarray(adata.obsm[rep_key], dtype=np.float32)

        fit_lines = np.unique(groups[fitc])
        y_lines, obs_lines = line_level(Y, M, groups, fit_lines)
        fns = (fit_weight_fns(y_lines, obs_lines, alpha=alpha, cap=cap)
               if density_weighting else None)

        fit_bags = build_bags(X_all, Y, M, groups, fitc, universe=eligible, weight_fns=fns)
        stop_bags = build_bags(X_all, Y, M, groups, stopc, universe=eligible, weight_fns=fns)
        val_bags = build_bags(X_all, Y, M, groups, vac, universe=eligible)

        model = OncoMLP(
            input_dim=X_all.shape[1],
            hidden_dims=hidden_dims,
            dropout_rate=dropout,
            input_dropout=input_dropout,
            norm="layer",
            output_dim=len(drugs),
        )
        if init_head_bias:
            init_head_bias_(model, per_drug_line_mean(y_lines, obs_lines))

        best, hist = train_bag_model(
            model, fit_bags, stop_bags, config=config, tag=f"{tag}_f{fold}"
        )
        best = best.to(device).eval()
        with torch.no_grad():
            per_cell, _ = bag_predictions(best, val_bags, device)
        # Written back through `val_bags.cells`, never through `vac`: build_bags sorts cells by line
        # so that a bag is a contiguous slice, and the two orders coincide only by accident.
        pred[val_bags.cells] = per_cell.cpu().numpy()

        val_lines = np.unique(groups[vac])
        folds.append({
            "fold": fold, "rep": rep, "n_train_lines": int(np.unique(groups[trc]).size),
            "n_fit_lines": int(fit_lines.size),
            "n_stop_lines": int(np.unique(groups[stopc]).size),
            "n_val_lines": int(val_lines.size),
            # One bag is one line, so these are also the numbers of training examples and of gradient
            # steps per epoch -- which is the point of recording them separately from the cell counts.
            "n_fit_bags": len(fit_bags), "n_fit_cells": int(fit_bags.X.shape[0]),
            "largest_bag": int(np.diff(fit_bags.bounds).max()),
            "val_lines": val_lines.tolist(),
            "best_epoch": hist.best_epoch, "best_val_obj": float(hist.best_val_mse),
        })
    return pred, folds
