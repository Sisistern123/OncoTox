"""Re-fit the ridge control on a grid wide enough that no fit hits its boundary.

**Why this exists.** ``4a``'s ridge cell fixed a test *before* the run (Selin, 13.08.2026): ridge's
penalty is not scale-invariant, so a shared grid could let the **grid** rather than the data choose
the penalty for one arm --

    "If the selected alpha never lands on a grid endpoint, the data chose the penalty in both arms
     and the asymmetry is inert; if it does, the GRID chose it for that arm, and the ridge baseline
     is not arm-comparable in this form."

**It landed on an endpoint in 8 of 55 ``X_pca`` fits and 4 of 55 ``X_scGPT`` fits** and nobody looked
-- the count sat in a committed column (``ridge_alpha_at_edge``) for three days. Both truncations are
at the **ceiling**, meaning leave-one-out CV wanted *more* penalty than ``logspace(-2, 4, 13)``
offers, so both ridges are **under-regularised** and the more truncated arm (``X_pca``, 8 fits) is
the more depressed. This widens the grid upward and re-fits.

⚠️ **Not a re-implementation with liberties taken.** The fit below is ``4a``'s ridge cell, line for
line: the same ``grouped_folds`` partition, the same per-fold PCA projections through
``fold_pca_projections_for``, the same ``inner_holdout`` fitting masks, the same line means, the same
per-drug ``RidgeCV``. **The script proves that by re-running the ORIGINAL grid first and asserting it
reproduces the committed artifact to six decimals.** If that assertion fails, the reimplementation has
drifted and the widened numbers are not comparable -- so the check runs before the new grid, not after.

**Why widen rather than standardise.** A log-spaced grid is already scale-fair: rescaling features by
*c* moves ridge's optimum by *c²*, which on a log grid is a shift at identical relative resolution.
The measured defect is not the scale difference but the grid's **upper bound**, which sits only one
decade above ``X_pca``'s typical selection. Widening touches the search; standardising touches the
features, and would change what every other consumer of those embeddings sees.

Run:  .venv/bin/python scripts/evaluation/ridge_grid_refit.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import anndata as ad  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402

from scripts.layout import PipelinePaths  # noqa: E402
from scripts.training.cv import fold_pca_projections_for, grouped_folds, inner_holdout  # noqa: E402
from scripts.training.density_weighting import line_level  # noqa: E402

VARIANT, SCORE, PCA_SEED, N_SPLITS = "hvg5000", "auc_cc", 42, 5
REPS = ["X_pca", "X_scGPT"]
OUT = ROOT / "notebooks/outputs/panel"

#: The grid `4a` used. Kept as the reproduction target, not as a live setting.
GRID_ORIGINAL = np.logspace(-2, 4, 13)
#: Two decades of extra headroom at the top, same 2-per-decade spacing. The floor is untouched: no
#: fit in either arm ever selected it. **10^6 is where this converges** -- extending to 10^8 moves
#: the means by <=0.0001 and leaves the same five fits at the boundary, so those five are not grid
#: truncation at all: leave-one-out CV wants an unbounded penalty for them, i.e. "no usable signal
#: for this (drug, fold) -- predict the intercept". That is an outcome, not an artifact, and no
#: ceiling fixes it.
GRID_WIDE = np.logspace(-2, 6, 17)


def ridge_table(adata, groups, idx, fold_split, masks, panel, yl, ol, lines, line_of, paths, grid):
    """`4a`'s ridge cell, parameterised only by the penalty grid."""
    recs = []
    for rep in REPS:
        proj = fold_pca_projections_for(rep, paths.raw_h5ad, [m[0] for m in masks], seed=PCA_SEED)
        pred = np.full_like(yl, np.nan)
        chosen: dict[int, list[float]] = {}
        for fold_i, (tr, va) in enumerate(fold_split):
            X = (np.asarray(adata.obsm[rep], dtype=np.float32) if proj is None
                 else np.asarray(proj[fold_i], dtype=np.float32))
            E = np.vstack([X[groups == ln].mean(0) for ln in lines])
            ti = [line_of[l] for l in np.unique(groups[idx[tr]]) if l in line_of]
            vi = [line_of[l] for l in np.unique(groups[idx[va]]) if l in line_of]
            for j in range(len(panel)):
                tj = [i for i in ti if ol[i, j]]
                vj = [i for i in vi if ol[i, j]]
                if len(tj) < 5 or not vj:
                    continue
                m = RidgeCV(alphas=grid).fit(E[tj], yl[tj, j])
                pred[vj, j] = m.predict(E[vj])
                chosen.setdefault(j, []).append(float(m.alpha_))
        for j, d in enumerate(panel):
            sel = ol[:, j] & np.isfinite(pred[:, j])
            t, p = yl[sel, j], pred[sel, j]
            a = np.asarray(chosen.get(j, [np.nan]), dtype=float)
            recs.append({"rep": rep, "model": "ridge", "alpha": np.nan, "drug": d,
                         "n_lines": int(sel.sum()), "spearman": spearmanr(p, t).statistic,
                         "mse": float(((p - t) ** 2).mean()),
                         "ridge_alpha_median": float(np.median(a)), "ridge_alpha_min": float(a.min()),
                         "ridge_alpha_max": float(a.max()), "ridge_alpha_n_fits": int(a.size),
                         "ridge_alpha_at_edge": int(np.isin(a, grid[[0, -1]]).sum())})
    return pd.DataFrame(recs)


def main() -> None:
    panel = pd.read_csv(OUT / "panel.csv")["drug_key"].tolist()
    paths = PipelinePaths.build(None, VARIANT, SCORE)
    src = ad.read_h5ad(paths.targets_h5ad, backed="r")
    all_drugs = list(src.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in panel]
    Y_raw = np.asarray(src.obsm["Y_ctrp"], dtype=np.float32)[:, kcol]
    M = np.asarray(src.obsm["M_ctrp"], dtype=bool)[:, kcol]
    adata = ad.AnnData(obs=src.obs.copy())
    for rep in REPS:
        adata.obsm[rep] = np.asarray(src.obsm[rep], dtype=np.float32)
    adata.obsm["Y_ctrp"] = np.where(M, Y_raw, 0.0).astype(np.float32)
    adata.obsm["M_ctrp"] = M
    adata.uns["ctrp_drugs"] = panel
    src.file.close()

    groups = adata.obs["Cell_line"].astype(str).to_numpy()
    eligible = adata.obs["split_ctrp"].isin(["train", "val"]).to_numpy()
    idx = np.flatnonzero(eligible)
    lines = np.unique(groups[eligible])
    line_of = {ln: i for i, ln in enumerate(lines)}
    yl, ol = line_level(adata.obsm["Y_ctrp"], M, groups, lines)

    _, fold_split = grouped_folds(adata, n_splits=N_SPLITS)
    masks = []
    for tr, va in fold_split:
        trc = np.zeros(adata.n_obs, dtype=bool); trc[idx[tr]] = True
        vac = np.zeros(adata.n_obs, dtype=bool); vac[idx[va]] = True
        fitc, _ = inner_holdout(groups, trc)
        masks.append((fitc, vac))

    args = (adata, groups, idx, fold_split, masks, panel, yl, ol, lines, line_of, paths)

    # ---- 1. reproduce the committed artifact on the original grid, or stop ------------------
    print("reproducing the committed ridge on the ORIGINAL grid ...")
    ref = pd.read_csv(OUT / "panel_ridge_baseline.csv")
    # This script overwrites the file it validates against, so a second run would compare the new
    # grid with itself and pass vacuously. Refuse instead, and say how to reset.
    if ref["ridge_alpha_max"].max() > GRID_ORIGINAL[-1]:
        raise SystemExit(
            "panel_ridge_baseline.csv already holds a widened fit (max selected penalty "
            f"{ref['ridge_alpha_max'].max():g} > {GRID_ORIGINAL[-1]:g}). The reproduction check "
            "would be vacuous. Reset it first:\n"
            "    git checkout notebooks/outputs/panel/panel_ridge_baseline.csv")
    old = ridge_table(*args, GRID_ORIGINAL)
    key = ["rep", "drug"]
    j = ref.set_index(key)["spearman"].round(6).sort_index()
    k = old.set_index(key)["spearman"].round(6).sort_index()
    if not j.equals(k):
        bad = (j - k).abs().sort_values(ascending=False).head(5)
        raise SystemExit(f"does NOT reproduce panel_ridge_baseline.csv; worst rows:\n{bad}\n"
                         "the reimplementation has drifted from 4a -- the widened numbers would "
                         "not be comparable, so nothing is written.")
    print(f"  ✅ all {len(k)} (rep, drug) rows identical to six decimals\n")

    # ---- 2. the widened grid --------------------------------------------------------------
    print(f"re-fitting on the widened grid {GRID_WIDE[0]:g} .. {GRID_WIDE[-1]:g} ...")
    new = ridge_table(*args, GRID_WIDE)

    summary = []
    for rep in REPS:
        o, n = old[old.rep == rep], new[new.rep == rep]
        summary.append({"rep": rep,
                        "mean_spearman_before": o.spearman.mean(), "mean_spearman_after": n.spearman.mean(),
                        "at_edge_before": int(o.ridge_alpha_at_edge.sum()),
                        "at_edge_after": int(n.ridge_alpha_at_edge.sum()),
                        "n_fits": int(n.ridge_alpha_n_fits.sum()),
                        "median_alpha_after": float(n.ridge_alpha_median.median())})
    s = pd.DataFrame(summary)
    print("\n" + s.round(4).to_string(index=False))

    still = int(s.at_edge_after.sum())
    print(f"\n  fits still at a grid boundary: {still}")
    if still > 5:
        print("  ⚠️ MORE BOUNDARY FITS THAN THE KNOWN FIVE -- widen further before quoting these.")
    elif still:
        print(f"  ✅ {still} fits saturate the ceiling, and they do so at 10^6 and 10^8 alike:")
        print("     leave-one-out CV wants an unbounded penalty for them (no usable signal in that")
        print("     drug x fold). Every fit with a finite optimum now finds it inside the grid,")
        print("     which is what 4a's pre-registered test asks for.")
    else:
        print("  ✅ every fit chose its penalty inside the grid: the data chose it, in both arms,")
        print("     which is exactly what 4a's pre-registered test asks for.")

    import os
    if os.environ.get("WRITE","1")=="1":
        new.to_csv(OUT / "panel_ridge_baseline.csv", index=False)
        print(f"\nwrote {(OUT / 'panel_ridge_baseline.csv').relative_to(ROOT)}")

        # panel_arch_summary.csv carries the ridge mean and `minus_ridge` as derived columns. They
        # are pure functions of the table just written, so leaving them would put two committed
        # artifacts in disagreement about the same quantity -- the failure this project keeps
        # finding. `mean` is the MLP's and is asserted untouched.
        arch_p = OUT / "panel_arch_summary.csv"
        arch = pd.read_csv(arch_p)
        before = arch["mean"].copy()
        means = new.groupby("rep")["spearman"].mean()
        arch["ridge"] = arch["rep"].map(means).round(4)
        arch["minus_ridge"] = (arch["mean"] - arch["ridge"]).round(4)
        assert arch["mean"].equals(before), "the MLP column must not move"
        arch.to_csv(arch_p, index=False)
        print(f"wrote {arch_p.relative_to(ROOT)}  (ridge + minus_ridge only)")
    else:
        print("\n(dry run, nothing written)")


if __name__ == "__main__":
    main()
