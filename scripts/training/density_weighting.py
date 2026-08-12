"""Inverse label-density sample weights for imbalanced regression (one implementation).

The response values of a compound are not spread evenly over the cell lines: most sit in a narrow
band and the pharmacologically interesting extremes are sparse. Under unweighted squared error the
crowded middle owns the gradient and the optimum shrinks toward the compound's mean. The remedy is to
weight each observation by the inverse of a *smoothed* estimate of the label density:

    Yang, Zha, Chen, Wang & Katabi. Delving into Deep Imbalanced Regression. ICML 2021 (LDS).
    Steininger, Kobs, Davidson, Krause & Hotho. Density-based weighting for imbalanced regression.
    Machine Learning 110 (2021) -- DenseWeight, source of the compression exponent ``alpha``.

Gaussian KDE (Scott's rule) stands in for LDS's binned-and-smoothed density; the bandwidth is the one
real knob and is currently unexamined.

**Fit per drug, on cell lines.** The evaluation metric is within-drug rank correlation, so a density
pooled over compounds would conflate potency differences between them with sensitivity differences
between lines. And the density must be fitted on line-level values: the label is constant within a
line, so fitting on cells would let a line with 1,990 sequenced cells bend the density toward its own
value while a line with 56 barely registers.

**Parameter values — settled by audit 09 (Selin, 12.08.2026): ``alpha`` becomes an arm of the loss
comparison rather than a constant needing a prior justification.** The evidence for ``alpha=0.5`` had
stopped being live: it was recorded in ``panel_distributions.ipynb``, archived 12.08.2026 with the void
8-drug panel it was computed on. Re-deriving it would have meant justifying a number in advance in order
to test a hypothesis the same run was about to test anyway, so instead the loss comparison sweeps
``alpha`` over **{off, 0.5, 1.0}** and the comparison *is* the justification. Nothing here has to be
argued before it runs; what wins is reported with the setting that won.

The original reasoning is kept because it explains what the sweep's endpoints mean: at ``alpha=1`` (full
inverse density) with a loose cap the weight curve saturates the cap across wide stretches of the
response range, so the cap -- an arbitrary safety limit -- sets the weights rather than the density
does. ``alpha=0.5`` compresses the range so the density still orders the samples without dictating the
magnitude.

``cap=3`` is **held fixed and is arbitrary**, documented as such rather than derived. It exists to stop
a handful of extreme lines driving a fit over ~150 points, and it is a safety limit, not a modelling
choice. Sweeping it as well would confound the alpha result, since alpha and cap trade against each
other by construction.

⚠️ Note that inverse-density weighting remains a
`refuted hypothesis <../../docs/steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking>`_
on the void panel and the retired target. It is re-tested rather than assumed, and it is re-tested
because the explanation given for that null rested on a winsorization the pipeline no longer applies.

Note that clipping and renormalizing interact: renormalizing after a clip can push values back over
the cap (an earlier version reported a maximum of 13.1 under a cap of 10). The normalizing constant is
therefore found by iterating to a fixed point, and the result is asserted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import gaussian_kde

DEFAULT_ALPHA = 0.5
DEFAULT_CAP = 3.0
# RETIRED CONSTANT `DEFAULT_WINSOR`. Kept as the record of why it went; it documents nothing below
# this line. Was a `#:` doc-comment until 12.08.2026, which made Sphinx attach it to WeightFn, the
# next declaration -- a paragraph about a deleted constant appearing as that class's documentation.
# Retired 11.08.2026 (Selin) -- do not reintroduce without a source. It clipped the response at
# 1.1 on the theory that anything above it was assay artifact, and existed so that inverse-density
# weighting would not hand the sparse upper tail the largest weights. It goes because: the threshold
# was never sourced beyond this comment; the benchmark we follow applies **no** quality filter or
# clipping at all, "to maintain comparability to previous studies and avoid data loss" (DrEval,
# Methods); the weighting it protected is itself a refuted hypothesis; and most of the old upper tail
# was an artifact of a defective normalisation rather than data. On the current ``auc_cc`` target only
# 3.65 % of measurements exceed 1.1 and p99 is 1.185. See
# docs/steps/01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026.
# If audit 09 keeps density weighting, clipping the weighting's *input* -- never the target -- may be
# reopened there.


@dataclass(frozen=True)
class WeightFn:
    """Weight as a function of label value, fitted on one drug's cell-line values.

    Callable on any array of label values -- the fitting values themselves, held-out cells, or a grid
    for plotting the curve -- and always with the normalization frozen at fit time, so weights on
    unseen values are on the same scale as the ones the model trained with.
    """

    kde: gaussian_kde
    norm: float
    alpha: float = DEFAULT_ALPHA
    cap: float = DEFAULT_CAP

    def __call__(self, y: np.ndarray) -> np.ndarray:
        raw = np.clip(self.kde(np.asarray(y, dtype=float)), 1e-12, None) ** (-self.alpha)
        return np.clip(raw / self.norm, 1.0 / self.cap, self.cap)


def fit_weight_fn(
    values: np.ndarray,
    *,
    alpha: float = DEFAULT_ALPHA,
    cap: float = DEFAULT_CAP,
    max_iter: int = 50,
    tol: float = 1e-9,
) -> WeightFn:
    """Fit the weight function for one drug, normalized to mean weight 1 over ``values``.

    Mean 1 means the drug's total contribution to the loss is unchanged, so the weighting acts
    strictly *within* the drug and cannot reintroduce an imbalance between drugs.
    """
    v = np.asarray(values, dtype=float)
    kde = gaussian_kde(v)
    raw = np.clip(kde(v), 1e-12, None) ** (-alpha)
    norm = float(raw.mean())
    for _ in range(max_iter):
        w = np.clip(raw / norm, 1.0 / cap, cap)
        norm_new = norm * float(w.mean())
        if abs(norm_new - norm) < tol:
            norm = norm_new
            break
        norm = norm_new

    fn = WeightFn(kde=kde, norm=norm, alpha=alpha, cap=cap)
    w = fn(v)
    assert abs(w.mean() - 1) < 1e-6, "weights must average to 1 over the fitting values"
    assert w.max() / w.min() <= cap**2 + 1e-6, "weight spread must respect the cap"
    return fn


def fit_weight_fns(
    y_lines: np.ndarray,
    obs_lines: np.ndarray,
    *,
    alpha: float = DEFAULT_ALPHA,
    cap: float = DEFAULT_CAP,
) -> list[WeightFn]:
    """One :class:`WeightFn` per drug from a (cell line x drug) label matrix and its mask.

    Pass only the *training* lines of a fold: the density is a function of the labels, so fitting it
    on all lines would let held-out labels inform training.
    """
    return [fit_weight_fn(y_lines[obs_lines[:, j], j], alpha=alpha, cap=cap)
            for j in range(y_lines.shape[1])]


def weight_matrix(fns: list[WeightFn], y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """(n, K) weights at observed entries, 0 elsewhere -- drop-in for a 0/1 mask.

    The masked loss computes ``sum(err * mask) / sum(mask)``, so substituting this for the mask turns
    it into ``sum(w * err) / sum(w)`` exactly, with no change to the training code. Unobserved entries
    stay 0 and remain excluded.
    """
    W = np.zeros_like(np.asarray(y), dtype=np.float32)
    for j, fn in enumerate(fns):
        sel = np.asarray(mask)[:, j]
        if sel.any():
            W[sel, j] = fn(np.asarray(y)[sel, j]).astype(np.float32)
    return W


def line_level(values: np.ndarray, mask: np.ndarray, groups: np.ndarray,
               lines: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Collapse cell-level labels to (n_lines, K) values and observation mask.

    The label is constant within a cell line, so this is a lookup rather than an average -- but it has
    to happen before any per-drug statistic, otherwise cell-rich lines are counted many times over.
    """
    out = np.full((len(lines), values.shape[1]), np.nan)
    obs = np.zeros_like(out, dtype=bool)
    for i, ln in enumerate(lines):
        ci = np.flatnonzero(groups == ln)
        m = mask[ci].any(0)
        out[i, m] = values[ci][:, m][0]
        obs[i] = m
    return out, obs
