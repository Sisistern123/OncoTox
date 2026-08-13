"""Cross-validated out-of-fold predictions -- the one implementation every caller should use.

``cv_evaluate`` in :mod:`scripts.training.train_multitask` answers "is the difference real?" and
returns per-fold *metrics*. This module returns the *predictions*, which is what any analysis
downstream of training needs: per-drug correlations, the cell-line-effect normalization, calibration
diagnostics. Those were previously re-implemented in each notebook that wanted them, with small
divergences (line-level versus cell-level aggregation, different epoch counts, bias initialized or
not), which is exactly how two runs stop being comparable.

The protocol is the project standard and does not vary between callers: K-fold ``GroupKFold`` over
cell lines, restricted to the ``split_ctrp`` train+val lines so the fixed test set is never touched,
and every held-out line predicted by a model that has not seen it -- neither in training nor in
early stopping, which watches a slice of the training lines instead (:func:`inner_holdout`, 12.08.2026;
before that it watched the scored fold). Per-cell predictions are averaged
back to one value per line by :func:`line_level_predictions`, because the label is per line and
scoring per cell would score pseudo-replicates.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy import sparse
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from torch.utils.data import DataLoader

from scripts.model.dataset import MultiDrugDataset
# `_pca_fitted_on_train` is private to add_pca, and is imported here deliberately rather than
# reimplemented: the per-fold fit must be the *same* fit as the descriptive one, differing only in
# which cells it sees. A second implementation is how those two silently drift apart.
from scripts.preprocessing.add_pca import DEFAULT_N_COMPS, _pca_fitted_on_train
from scripts.preprocessing.expression import kinker_transform
from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS, OncoMLP, init_head_bias_
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


INNER_VAL_FRACTION = 0.15
INNER_SPLIT_SEED = 42


def inner_holdout(
    groups: np.ndarray,
    train_cells: np.ndarray,
    *,
    frac: float = INNER_VAL_FRACTION,
    seed: int = INNER_SPLIT_SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Split one fold's training cells into a fitting set and an early-stopping set.

    Until 12.08.2026 a fold's held-out lines served as the early-stopping set **and** as the
    scored set, so every out-of-fold prediction came from the checkpoint that best fit the lines
    it was about to be scored on. The same defect was found in the DrEval benchmark on
    14.07.2026 and fixed there (``ee07b00``) without being fixed here
    ([corrections](../../docs/steps/corrections-and-dead-ends.md)). Nesting the early-stopping
    set inside the training lines is the fix (Selin, 12.08.2026).

    ``frac`` is a proportion of the fold's training **lines**, not of its cells --
    ``GroupShuffleSplit`` splits groups -- so no line contributes to both sides. **0.15 is
    arbitrary:** it is the conventional size of a validation slice and nothing in this data sets
    it; it was chosen over reusing a whole neighbouring fold (20 %) to spend fewer training
    lines.

    :returns: ``(fit_cells, stop_cells)`` -- boolean arrays over all cells, disjoint, and
        together equal to ``train_cells``.
    :raises ValueError: if either side comes out empty, which means the fold is too small for
        this fraction rather than that the split is merely unlucky.
    """
    pos = np.flatnonzero(train_cells)
    splitter = GroupShuffleSplit(n_splits=1, test_size=frac, random_state=seed)
    fit_i, stop_i = next(splitter.split(pos, groups=groups[pos]))
    if fit_i.size == 0 or stop_i.size == 0:
        raise ValueError(
            f"inner split of {np.unique(groups[pos]).size} training lines at frac={frac} left "
            f"{np.unique(groups[pos[fit_i]]).size} fitting and "
            f"{np.unique(groups[pos[stop_i]]).size} early-stopping lines."
        )

    fit_cells = np.zeros_like(train_cells, dtype=bool)
    stop_cells = np.zeros_like(train_cells, dtype=bool)
    fit_cells[pos[fit_i]] = True
    stop_cells[pos[stop_i]] = True
    return fit_cells, stop_cells


_FOLD_PCA_CACHE: dict[str, list[np.ndarray]] = {}


def _fold_key(counts_h5ad: Path, fit_masks: list[np.ndarray], n_comps: int, seed: int) -> str:
    """Digest of the fold assignment itself, not of fold indices.

    Two runs with different eligible line sets must never share an entry, so the key is built from
    the masks' contents. Fold *index* would collide across different partitions of the same size.
    """
    h = hashlib.sha256()
    h.update(f"{counts_h5ad}|{n_comps}|{seed}".encode())
    for m in fit_masks:
        h.update(np.packbits(np.asarray(m, dtype=bool)).tobytes())
    return h.hexdigest()


def fold_pca_projections(
    counts_h5ad: str | Path,
    fit_masks: list[np.ndarray],
    *,
    n_comps: int = DEFAULT_N_COMPS,
    seed: int = 42,
) -> list[np.ndarray]:
    """One PCA per fold, each fitted on that fold's own cells: ``[(n_cells, n_comps)] * n_folds``.

    **Decision 2 of the three fits (Selin, 12.08.2026): the cross-validated PCA is fitted per fold,
    at training time, and this closes the leak rather than documenting it.** Until then CV read the
    all-cells ``X_pca``, so a held-out line's coordinates depended on held-out lines --
    ``resolve_rep`` returns ``use_rep`` unchanged when ``split_col is None``, which is the CV case.

    **It must read the counts file, not the targets ``.X``.** The two carry different gene sets: PCA
    is fitted on the full HVG set (5,000 for ``hvg5000``) while the targets ``.X`` has scGPT's
    out-of-vocabulary genes dropped (4,704 after the symbol repair). Refitting from the targets file
    would silently lose 296 genes relative to the descriptive fit and produce coordinates that are
    not comparable with it. Measured 12.08.2026; this is the reason the extra 2.15 GB read exists,
    and the reason not to "simplify" it to whatever matrix is already in memory.

    **Load once, fit every fold, release.** The counts matrix is dense 53,513 x 5,000, cast to
    float32 on load -- 1.07 GB rather than the 2.14 GB it occupies on disk as float64. The five
    projections it produces total 0.55 GB. Reading the file per fold would cost 10.75 GB of I/O for
    an identical result. Peak is ~2.2 GB during the fits and ~0.55 GB afterwards, because the matrix
    is dropped before the caller starts training.

    Results are cached **in process** under a digest of the fold assignment. Not on disk: a cached
    projection under ``runs/`` would be an artifact outliving the code that produced it -- the exact
    class the 03.08 freeze exists to police -- and gitignored, so invisible in review, and failing by
    producing plausible wrong numbers rather than an error.

    **float32 throughout, here and in ``add_pca`` (Selin, 12.08.2026).** Both fits were extended
    together so they stay identical in kind, preserving the property she established on 10.08.2026
    when she harmonized ``ddof`` to 1 -- that the descriptive fit and the train-only fits differ
    only in which cells they see. Extending it makes ``X_pca`` move in its last digits, so this is a
    preprocessing change and R2 runs after it.

    **Fitted on the fold's ``fit_cells``, not its whole training side (Selin, 12.08.2026).** The two
    arms already receive the same cells, the same per-cell CPM and log transform, and the same
    all-cells HVG gene set. The single asymmetry between them is that **PCA needs a fit estimated
    across cells -- mean, std, rotation -- and scGPT needs none**: its weights are pretrained and
    frozen, and its value binning digitizes each cell against that cell's own distribution. Because
    that fit is the only asymmetry, the narrower it is, the more of any measured difference is
    attributable to the representation rather than to how much the control's fit was allowed to see.
    Using the whole training side would hand the control ~18 % more cells that scGPT structurally
    cannot receive.

    It also removes an exception: PCA is now fitted on the same set as ``fit_weight_fns`` and the
    head-bias initialization, so no per-fold statistic is different for a reason someone has to
    remember.

    ⚠️ **The cost, which is not hidden:** ~15 % fewer cells for a 512-component fit, so a scGPT win
    invites the question whether PCA was undersold. The answer is that the fit still uses ~85 % of
    the fold's training cells, and the alternative would give the control an advantage the other arm
    cannot have.

    :param fit_masks: one boolean mask over all cells per fold, marking the cells that fold's PCA may
        be fitted on. Whatever set the caller passes is the set the fit sees.
    """
    counts_h5ad = Path(counts_h5ad)
    key = _fold_key(counts_h5ad, fit_masks, n_comps, seed)
    if key in _FOLD_PCA_CACHE:
        return _FOLD_PCA_CACHE[key]

    src = sc.read_h5ad(counts_h5ad)
    kinker_transform(src)  # log2(1 + CPM/10), the dataset authors' own transform -- as add_pca does
    X_log = (src.X.toarray() if sparse.issparse(src.X) else np.asarray(src.X)).astype(
        np.float32, copy=False
    )
    genes = np.asarray(src.var_names)
    projections = []
    for fold, mask in enumerate(fit_masks, start=1):
        n_fit = int(np.asarray(mask, dtype=bool).sum())
        if n_comps > min(n_fit, X_log.shape[1]):
            raise ValueError(
                f"n_comps={n_comps} exceeds min(n_fit={n_fit}, n_genes={X_log.shape[1]}) for fold "
                f"{fold}. A fold cannot be projected onto more components than it has cells."
            )
        proj, _record = _pca_fitted_on_train(
            X_log, np.asarray(mask, dtype=bool), n_comps, seed,
            genes=genes, fitted_on=f"{n_fit} cells in the fitting set of CV fold {fold}",
        )
        projections.append(proj)
        print(f"  fold {fold}: PCA fitted on {n_fit} cells -> {proj.shape}")
    del X_log, src  # the caller trains next; it needs the projections, not the gene matrix
    _FOLD_PCA_CACHE[key] = projections
    return projections


#: Scratch ``obsm`` key holding the current fold's PCA projection. Overwritten per fold; never
#: persisted, because it is meaningful only alongside the fold assignment that produced it.
CV_FOLD_PCA_KEY = "_X_pca_cv_fold"


def fold_pca_projections_for(
    rep: str,
    counts_h5ad: str | Path | None,
    fit_masks: list[np.ndarray],
    *,
    seed: int = 42,
) -> list[np.ndarray] | None:
    """Per-fold projections when ``rep`` is a PCA representation, ``None`` when it is not.

    The one place that decides whether cross-validation needs a per-fold fit, so that
    :func:`oof_predictions` and :func:`cv_evaluate` cannot answer it differently. scGPT needs
    nothing: its weights are pretrained and frozen, and its value binning digitizes each cell
    against that cell's own distribution, so no statistic is estimated across cells at all.

    :raises ValueError: if a PCA rep is requested without ``counts_h5ad``. Falling back to the
        stored all-cells ``X_pca`` is exactly the leak decision 2 closes, and a silent fallback is
        how it would come back.
    """
    if not rep.startswith("X_pca"):
        return None
    if counts_h5ad is None:
        raise ValueError(
            f"rep={rep!r} under cross-validation needs `counts_h5ad` so the PCA can be fitted per "
            f"fold (Selin, 12.08.2026). Without it the only available projection is the all-cells "
            f"X_pca, whose coordinates depend on the held-out lines -- the leak this replaces. Pass "
            f"PipelinePaths.raw_h5ad, or use a representation that needs no fit."
        )
    return fold_pca_projections(counts_h5ad, fit_masks, seed=seed)


def per_drug_line_mean(y_lines: np.ndarray, obs_lines: np.ndarray) -> np.ndarray:
    """Each drug's mean over the cell lines that were screened against it.

    Per **line**, not per cell: the label is broadcast to every cell of a line, so a cell-level mean
    would weigh a 1,990-cell line 35x a 56-cell line for the same single measurement. Drugs with no
    observation in the input get 0.0 -- a head with nothing to fit has no mean to start from, and 0
    is the untouched-``nn.Linear`` neighbourhood rather than a value invented for it.

    Take ``y_lines`` / ``obs_lines`` from :func:`density_weighting.line_level` over the fitting lines
    of a fold. The same array serves the head-bias initialization and any per-drug null.
    """
    return np.array(
        [y_lines[obs_lines[:, j], j].mean() if obs_lines[:, j].any() else 0.0
         for j in range(y_lines.shape[1])],
        dtype=np.float32,
    )


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
    counts_h5ad: str | Path | None = None,
    pca_seed: int | None = None,
    n_label_lines: int | None = None,
    tag: str = "oof",
) -> tuple[np.ndarray, list[dict]]:
    """Fit ``n_splits`` cell-line-grouped folds and return out-of-fold per-cell predictions.

    Each fold trains on its training lines minus a 15 % early-stopping slice
    (:func:`inner_holdout`), so the checkpoint that predicts a held-out line was chosen without
    looking at it. Everything fitted per fold -- the density weighting, the head bias -- is fitted on
    the fitting lines alone.

    :param density_weighting: weight each observation by inverse label density
        (:mod:`scripts.training.density_weighting`). The density is fitted **inside each fold, on the
        fitting lines only** -- it is a function of the labels, so fitting it once over all lines
        would let held-out labels inform training. The weights are carried in the mask tensor, which
        turns the masked mean into a weighted one without touching the loss code.
    :param init_head_bias: initialize each head's bias to that drug's mean over the fold's *fitting*
        lines, so the model starts at the null predictor. Necessary on an uncentred target such as
        raw AUC, where the bias must reach ~0.7 from a default near 0; harmless on a centred one.
    :param pca_seed: the seed for the per-fold PCA fit, held **separate from the model seed**. They
        were the same value until 13.08.2026, which made "seed" mean two different things here and in
        :func:`mil.bag_oof_predictions`: a run over three model seeds also refitted the PCA three
        times, so the representation moved with the initialisation. That breaks the premise of `4b`'s
        stage 1, which compares the two architectures **seed by seed** and claims the architecture is
        the only thing differing between the columns. Defaults to ``config.seed``, the old behaviour,
        so no existing caller changes silently.

        It is also why an ``X_pca`` sweep was slow: the projection cache keys on this seed
        (:func:`_fold_key`), so three model seeds meant three full 512-component fits per alpha
        instead of one.
    :returns: ``(pred, folds)`` -- ``pred`` is (n_cells, len(drugs)) with NaN where a cell was never
        held out, ``folds`` is a per-fold log.
    """
    config = config or TrainConfig()
    # `is None`, NOT `or` (13.08.2026). `()` is a legitimate value -- it means "no hidden layer",
    # which `OncoMLP` resolves to a bare `Linear(input_dim, K)` and which the capacity control in
    # `4a` section C passes. Under `or` it is falsy, so that call would have silently received
    # DEFAULT_HIDDEN_DIMS and compared the trunk against itself, reporting a capacity comparison
    # between two identical networks.
    hidden_dims = tuple(DEFAULT_HIDDEN_DIMS[rep] if hidden_dims is None else hidden_dims)

    groups = adata.obs[group_col].astype(str).to_numpy()
    idx, fold_split = grouped_folds(
        adata, n_splits=n_splits, group_col=group_col, eligible_splits=eligible_splits
    )

    all_drugs = list(adata.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in drugs]
    Y = np.asarray(adata.obsm["Y_ctrp"], dtype=np.float32)[:, kcol]
    M = np.asarray(adata.obsm["M_ctrp"], dtype=bool)[:, kcol]

    # Every fold's masks are computed before any training, because the per-fold PCA below fits all
    # folds from one load of the counts matrix. The fold's training lines split once more: the model
    # is fitted on `fitc` and early stopping watches `stopc`. The scored fold `vac` enters neither,
    # so the checkpoint that predicts it was not chosen by it -- see inner_holdout above.
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

    pred = np.full((adata.n_obs, len(drugs)), np.nan, dtype=float)
    folds: list[dict] = []
    for fold, (trc, fitc, stopc, vac) in enumerate(masks, start=1):
        # Under a per-fold PCA the representation differs per fold, so it is written to a scratch
        # obsm key the datasets then read. Overwritten each iteration: MultiDrugDataset copies the
        # array at construction, so no fold can see another's projection.
        rep_key = rep
        if projections is not None:
            rep_key = CV_FOLD_PCA_KEY
            adata.obsm[rep_key] = projections[fold - 1]

        fit_ds = MultiDrugDataset(adata=adata, use_rep=rep_key, cell_mask=fitc, drugs=drugs)
        stop_ds = MultiDrugDataset(adata=adata, use_rep=rep_key, cell_mask=stopc, drugs=drugs)

        # per-drug statistics from the FITTING lines only (never cells, never held-out lines, and
        # not the early-stopping lines either: a set that decides when to stop must not also have
        # shaped what it is judging).
        fit_lines = np.unique(groups[fitc])

        # `n_label_lines` thins the LABEL supply without touching the input (13.08.2026, for `4a`
        # section E's learning curve). The cells of the dropped lines stay in the fold, stay in the
        # batches and stay in the per-fold PCA -- only their labels go. That is what separates "how
        # many labelled cell lines" from "how many cells the representation was fitted on": both
        # arms keep the identical input at every point on the curve, so what moves along it is the
        # label supply and nothing else.
        #
        # `label_lines` is a SEPARATE name on purpose. `fit_lines` is logged as `n_fit_lines` and
        # means "lines in this fold's fitting set"; reassigning it would silently redefine an
        # existing column to mean something else in some runs and not in others.
        #
        # The thinning reaches all THREE label-derived quantities below, not only the loss. Leaving
        # the head-bias init or the density fit on the full set would hand the model the dropped
        # lines' labels through the back door, and the curve would flatten for a reason that has
        # nothing to do with the representations.
        label_lines, unlabelled = fit_lines, None
        if n_label_lines is not None and n_label_lines < fit_lines.size:
            rng = np.random.default_rng(config.seed)
            label_lines = np.sort(rng.choice(fit_lines, size=n_label_lines, replace=False))
            unlabelled = ~np.isin(groups[fitc], label_lines)

        y_lines, obs_lines = line_level(Y, M, groups, label_lines)

        if density_weighting:
            fns = fit_weight_fns(y_lines, obs_lines, alpha=alpha, cap=cap)
            fit_ds.mask = torch.from_numpy(weight_matrix(fns, Y[fitc], M[fitc]))
            stop_ds.mask = torch.from_numpy(weight_matrix(fns, Y[stopc], M[stopc]))

        # After the density weighting, which rewrites the mask wholesale.
        if unlabelled is not None:
            fit_ds.mask[torch.from_numpy(unlabelled)] = 0

        model = OncoMLP(
            input_dim=fit_ds.X.shape[1],
            hidden_dims=hidden_dims,
            dropout_rate=dropout,
            input_dropout=input_dropout,
            norm="layer",
            output_dim=len(drugs),
        )
        if init_head_bias:
            init_head_bias_(model, per_drug_line_mean(y_lines, obs_lines))

        # Explicit generator: without one, RandomSampler seeds itself from the global torch
        # RNG at iteration time, which happens to be seeded because train_model calls
        # set_seed(config.seed) first -- an ordering dependency that breaks silently.
        best, hist = train_model(
            model,
            DataLoader(
                fit_ds,
                batch_size=batch_size,
                shuffle=True,
                generator=torch.Generator().manual_seed(config.seed),
            ),
            DataLoader(stop_ds, batch_size=batch_size, shuffle=False),
            config=config,
            tag=f"{tag}_f{fold}",
            drug_names=drugs,
        )
        best = best.to("cpu").eval()  # train_model leaves it on mps/cuda
        with torch.no_grad():
            # rep_key, not rep: under a per-fold PCA the held-out cells must be projected
            # through THIS fold's rotation, the one fitted without them.
            x = torch.from_numpy(np.asarray(adata.obsm[rep_key], dtype=np.float32)[vac])
            pred[vac] = best(x).numpy()  # held-out lines only

        val_lines = np.unique(groups[vac])
        folds.append({
            "fold": fold, "rep": rep, "n_train_lines": int(np.unique(groups[trc]).size),
            # n_train_lines is the fold's training side as a whole; the model only ever saw
            # n_fit_lines of it, the rest being on loan to early stopping.
            "n_fit_lines": int(fit_lines.size),
            "n_label_lines": int(label_lines.size),
            "n_stop_lines": int(np.unique(groups[stopc]).size),
            "n_val_lines": int(val_lines.size),
            # Which lines this fold held out, not just how many. The split is grouped by cell line, so
            # a line belongs to exactly one fold and this assignment is total and unambiguous. It is
            # recorded because anything fitting a baseline on the *training* folds needs it -- DrEval's
            # normalization above all (scripts/evaluation/dreval_normalize.py) -- and re-deriving it
            # downstream would mean a second copy of the fold rule that can drift from this one.
            "val_lines": val_lines.tolist(),
            # best_val_obj is measured on the early-stopping slice, not on the scored fold: it is the
            # quantity that chose the checkpoint, and it is not a generalization estimate. The
            # out-of-fold predictions are.
            "best_epoch": hist.best_epoch, "best_val_obj": float(hist.best_val_mse),
        })
    return pred, folds


def line_level_predictions(
    pred: np.ndarray,
    adata,
    drugs: list[str],
    *,
    truth: np.ndarray | None = None,
    folds: list[dict] | None = None,
    group_col: str = "Cell_line",
    eligible_splits: tuple[str, ...] = ("train", "val"),
    **extra,
) -> pd.DataFrame:
    """Average per-cell predictions to one row per (drug, cell line) and attach the truth.

    :param truth: (n_cells, len(drugs)) labels to score against; defaults to the ``Y_ctrp`` columns
        of ``drugs``. Pass an explicit array to score every model against one common yardstick even
        when they trained on different scores.
    :param folds: the fold log :func:`oof_predictions` returns. When given, each row gains a ``fold``
        column saying which fold held that cell line out. **Pass it** unless there is a reason not to:
        without it a prediction cannot be traced to the split that produced it, and any downstream
        baseline has to be fitted on the same rows it will be subtracted from --- which is what
        ``scripts/evaluation/dreval_normalize.py`` refuses to do.
    :param extra: constant columns to stamp on every row (e.g. ``rep=..., weighted=...``), so
        several configurations can be concatenated into one tidy table.
    :returns: DataFrame(drug, cell_line, y_true, y_pred[, fold], **extra) -- the shape every
        downstream analysis expects.
    :raises ValueError: if ``folds`` does not assign every predicted cell line to exactly one fold.
    """
    line_fold: dict[str, int] = {}
    if folds is not None:
        for entry in folds:
            for line in entry["val_lines"]:
                if line in line_fold:
                    raise ValueError(
                        f"cell line {line!r} is held out by folds {line_fold[line]} and "
                        f"{entry['fold']}; the split is not grouped by cell line."
                    )
                line_fold[line] = entry["fold"]

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
                row = {"drug": d, "cell_line": ln,
                       "y_true": float(Y[obs, j].mean()),
                       "y_pred": float(pred[obs, j].mean())}
                if folds is not None:
                    if ln not in line_fold:
                        raise ValueError(
                            f"cell line {ln!r} has predictions but no fold assigns it; the fold log "
                            f"does not describe the run these predictions came from."
                        )
                    row["fold"] = line_fold[ln]
                rows.append({**row, **extra})
    return pd.DataFrame(rows)
