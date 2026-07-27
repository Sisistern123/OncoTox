"""Cell-line-effect normalization of our own out-of-fold predictions (DrEval recipe, extended).

Reconstructs ``notebooks/outputs/dreval/dreval_normalized.csv``, whose producing code was lost
(the committed CSV predates the val-split fix ``ee07b00`` and no notebook or script wrote it).

**What this answers, and why DrEval itself cannot.** DrEval's normalized metric subtracts the
``NaiveMeanEffectsPredictor`` from truth *and* prediction, leaving only what a model adds beyond the
average drug and cell-line effect. But in the LCO setting a held-out cell line was never seen, so its
cell-line effect is 0 and the correction collapses to *removing the drug mean only*
(``notebooks/12_dreval_benchmark.ipynb``). That is the honest choice for a prediction benchmark --
no model can know an unseen line's effect -- but it leaves a diagnostic question open:

    is our per-drug correlation genuine drug-specific biology, or just "this cell line is fragile"?

Some lines are sensitive to everything (sigma of the line effect ~ 0.40 on the raw AUC scale). A model
can score a good per-drug Spearman by learning general fragility, with zero drug-specific signal. This
script removes the cell-line effect *as well*. It is deliberately **not** a legitimate predictor -- it
reads held-out labels -- it is a diagnostic that asks what survives when fragility is taken away.

Three numbers per (heads x rep x drug), Spearman over the ~150 out-of-fold cell lines:

    rho_raw            corr(prediction, truth)
    rho_normalized     corr(prediction - c, truth - c)    <- what the model adds beyond fragility
    rho_naive_baseline corr(c, truth)                     <- what fragility alone explains

where ``c`` is the **fragility index**: the cell line's mean response across all 545 drugs, in per-drug
z units. A drug with a high ``rho_naive_baseline`` and a ``rho_normalized`` near 0 had no drug-specific
signal at all -- its apparent performance *was* the cell-line effect. That is the artifact class DrEval
describes, measured on our own results.

**Why ``c`` is estimated from all 545 drugs and not from the evaluation panel.** Subtracting the same
estimated quantity from truth and prediction injects that estimate's noise into *both* residuals, so a
noisy ``c`` manufactures correlation out of nothing. Estimated from a panel of 3 drugs, a synthetic
model that has learned *only* fragility scores rho_normalized = 0.94 instead of 0. Averaged over ~545
drugs the estimation noise is negligible and the statistic behaves. The evaluated drug is additionally
left out of its own ``c`` (leave-one-drug-out), so no drug contributes to the quantity it is scored
against. Do not change this to a panel-based estimate.

**Usage -- prefer the route that trains nothing.** The normalization is post-processing of
out-of-fold predictions, and the fragility index is computed from labels, so model fitting only ever
served to *produce* the predictions. Notebook 14 already writes them, so:

    # 0 model fits: score the predictions notebook 14 produced
    uv run scripts/evaluation/dreval_normalize.py \\
        --oof-csv notebooks/outputs/panel/panel_oof_predictions.csv --weighted false

Refitting is the fallback for a configuration no run has covered yet. It costs 5 folds per
(heads x rep), so restrict the heads: the K=545 arm is the retired configuration and its contrast
value is already recorded in docs/steps/05 (raw auc scores -0.069 there against +0.377 on the panel).

    # 10 fits: panel heads only, both representations
    uv run scripts/evaluation/dreval_normalize.py --heads panel --drugs methotrexate dasatinib \\
        paclitaxel vincristine afatinib topotecan tanespimycin selumetinib

The training defaults follow the current setup (27.07.2026): winsorized at 1.1, output layer out of
weight decay, head bias at the train-fold per-drug means, 25 epochs. Otherwise a refit would
normalize a model that is no longer in use.

The cross-validation itself is :func:`scripts.training.cv.oof_predictions`, the project's single
implementation, so this script, ``notebooks/14_panel_training.ipynb`` and anything else that needs
out-of-fold predictions cannot drift apart: 5-fold GroupKFold over the train+val cell lines (the 27
test lines are never touched), each fold predicting only the lines it did not see, per-cell
predictions averaged to one value per cell line. Scoring against the common curve-fit ``auc`` ranking,
whatever score a model trained on, is this module's own addition.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

from scripts.preprocessing.layout import PipelinePaths, add_data_args
from scripts.training.cv import line_level_predictions, oof_predictions
from scripts.training.training_utils import TrainConfig

DEFAULT_OUT = Path("notebooks/outputs/dreval/dreval_normalized.csv")
DEFAULT_PANEL_CSV = Path("notebooks/outputs/learnability/ctrp_drug_learnability_auc.csv")
REPS = ("X_pca", "X_scGPT")


def load_panel(panel_csv: Path) -> list[str]:
    """Read the selected drug panel from a learnability CSV (the ``selected`` column)."""
    lrn = pd.read_csv(panel_csv)
    if "selected" not in lrn.columns:
        raise ValueError(f"{panel_csv} has no 'selected' column -- pass --drugs explicitly.")
    return lrn[lrn["selected"]].sort_values("learnability", ascending=False)["drug"].tolist()


def line_by_drug_truth(variant: str, data_root: Path | None) -> pd.DataFrame:
    """(cell line x all drugs) matrix of curve-fit ``auc``, averaged over each line's cells.

    Restricted to the train+val lines, i.e. the same population the out-of-fold predictions cover.
    NaN where a (line, drug) pair was never measured.
    """
    adata = sc.read_h5ad(PipelinePaths.build(data_root, variant, "auc").targets_h5ad)
    eligible = adata.obs["split_ctrp"].isin(["train", "val"]).to_numpy()
    lines = adata.obs["Cell_line"].astype(str).to_numpy()
    y = np.asarray(adata.obsm["Y_ctrp"], dtype=float)
    m = np.asarray(adata.obsm["M_ctrp"], dtype=bool)
    drugs = list(adata.uns["ctrp_drugs"])

    uniq = np.unique(lines[eligible])
    out = np.full((uniq.size, len(drugs)), np.nan)
    for i, line in enumerate(uniq):
        ci = np.flatnonzero((lines == line) & eligible)
        obs = m[ci]
        counts = obs.sum(axis=0)
        sums = np.where(obs, y[ci], 0.0).sum(axis=0)
        out[i] = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)
    return pd.DataFrame(out, index=pd.Index(uniq, name="cell_line"), columns=drugs)


def fragility_index(truth_all: pd.DataFrame, *, leave_out: str | None = None) -> pd.Series:
    """Per-line mean response across all drugs, in per-drug z units ("this line is fragile").

    Each drug column is standardized across cell lines first, so drugs with a wide raw AUC spread do
    not dominate the average. Negative = sensitive to everything, positive = resistant to everything.

    :param leave_out: drug to exclude, so a drug never contributes to the index it is scored against.
    """
    cols = [c for c in truth_all.columns if c != leave_out]
    sub = truth_all[cols]
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0)
    return z.mean(axis=1, skipna=True)


def load_oof_csv(path: Path, *, rep: str, weighted: bool | None = None) -> dict[str, pd.DataFrame]:
    """Read out-of-fold predictions produced elsewhere, in the shape the normalization wants.

    **This is the preferred route — it costs no model fits.** The normalization is pure
    post-processing of predictions, and the fragility index comes from labels alone, so refitting
    only ever served to *produce* the predictions. ``notebooks/14_panel_training.ipynb`` already
    writes them, keyed by cell line, to ``outputs/panel/panel_oof_predictions.csv``.

    Reusing them is also the more correct option, not merely the cheaper one: a refit here would
    normalize a model configuration that is no longer in use unless the arguments below are matched
    to it exactly.

    :param path: tidy CSV with columns rep, weighted, drug, cell_line, y_true, y_pred.
    :param weighted: pick one weighting arm; ``None`` requires the file to hold only one.
    """
    df = pd.read_csv(path)
    df = df[df["rep"] == rep]
    if weighted is not None and "weighted" in df.columns:
        df = df[df["weighted"].astype(bool) == weighted]
    if "weighted" in df.columns and df["weighted"].nunique() > 1:
        raise ValueError("CSV holds several weighting arms; pass weighted=True/False to choose one.")
    if df.empty:
        raise ValueError(f"no rows for rep={rep!r} in {path}")
    return {d: g[["cell_line", "y_true", "y_pred"]].reset_index(drop=True)
            for d, g in df.groupby("drug", sort=False)}


def oof_line_predictions(
    score: str,
    drugs: list[str],
    rep: str,
    eval_drugs: list[str],
    *,
    variant: str,
    data_root: Path | None,
    seed: int = 42,
    # 25, not 50: over 36 recorded runs the best epoch was median 6 / max 11 and early stopping
    # (patience 10) has never reached 25, so the extra epochs were pure wall-clock.
    epochs: int = 25,
    winsor: float | None = 1.1,
    exclude_output_from_decay: bool = True,
    init_head_bias: bool = True,
) -> dict[str, pd.DataFrame]:
    """Train K=len(drugs) heads and return out-of-fold per-cell-line predictions.

    Prefer :func:`load_oof_csv` when predictions already exist — this function is the expensive path.

    Defaults match the current training setup (27.07.2026, ``notebooks/14_panel_training.ipynb``) so
    that a refit normalizes the model actually in use. On a raw-AUC target that means: responses
    winsorized at ``winsor`` (above ~1.1 the compound apparently improved growth over control, i.e.
    assay artifact), the output layer kept out of weight decay (its bias must sit near each drug's
    mean, and decay pulls it to 0), and that bias initialized to the train-fold per-drug means so the
    model starts at the null predictor rather than at zero. Pass ``winsor=None,
    exclude_output_from_decay=False, init_head_bias=False, epochs=50`` to reproduce the pre-27.07
    protocol of ``notebooks/11_auc_vs_aucz.ipynb`` instead.

    :returns: drug -> DataFrame(cell_line, y_true, y_pred); ``y_true`` is always on the curve-fit
        ``auc`` scale so models trained on different scores share one yardstick.
    """
    truth = sc.read_h5ad(PipelinePaths.build(data_root, variant, "auc").targets_h5ad, backed="r")
    truth_y = np.asarray(truth.obsm["Y_ctrp"], dtype=float)
    truth_drugs = list(truth.uns["ctrp_drugs"])

    adata = sc.read_h5ad(PipelinePaths.build(data_root, variant, score).targets_h5ad)
    mask = np.asarray(adata.obsm["M_ctrp"], dtype=bool)
    if winsor is not None:
        y = np.asarray(adata.obsm["Y_ctrp"], dtype=np.float32)
        adata.obsm["Y_ctrp"] = np.where(mask, np.clip(y, None, winsor), 0.0).astype(np.float32)

    cfg = TrainConfig(
        epochs=epochs, seed=seed, log_every=1000,
        exclude_output_from_decay=exclude_output_from_decay,
    )
    # One CV implementation for the whole project (scripts/training/cv.py): grouped folds, per-fold
    # line-level statistics, optional density weighting, head-bias initialization.
    pred, _ = oof_predictions(
        adata, rep, drugs,
        config=cfg,
        init_head_bias=init_head_bias,
        tag=f"norm_{score}_{rep}_K{len(drugs)}_s{seed}",
    )

    # Score every model against the common curve-fit `auc` ranking, whatever it trained on.
    jcol = [drugs.index(d) for d in eval_drugs]
    tcol = [truth_drugs.index(d) for d in eval_drugs]
    tidy = line_level_predictions(
        pred[:, jcol], adata, eval_drugs, truth=truth_y[:, tcol]
    )
    return {d: g[["cell_line", "y_true", "y_pred"]].reset_index(drop=True)
            for d, g in tidy.groupby("drug", sort=False)}


def _partial_spearman(p: np.ndarray, t: np.ndarray, c: np.ndarray) -> float:
    """Spearman between ``p`` and ``t`` controlling for ``c``.

    DrEval subtracts its naive prediction from truth and prediction with a coefficient of 1, which is
    valid there because the naive predictor is fitted in the response's own units. The fragility index
    is not: it lives in per-drug z units, so a unit subtraction leaves a residual proportional to
    fragility in *both* series and the correlation survives almost intact (0.62 instead of 0 on the
    synthetic control). Regressing ``c`` out of both rank vectors is the scale-free generalization --
    the residuals are orthogonal to fragility by construction.
    """
    from scipy.stats import rankdata

    pr, tr, cr = (rankdata(v).astype(float) for v in (p, t, c))
    design = np.c_[np.ones(len(cr)), cr]
    resid_p = pr - design @ np.linalg.lstsq(design, pr, rcond=None)[0]
    resid_t = tr - design @ np.linalg.lstsq(design, tr, rcond=None)[0]
    return float(spearmanr(resid_p, resid_t).statistic)


def normalized_correlations(
    oof: dict[str, pd.DataFrame],
    truth_all: pd.DataFrame,
    *,
    heads: str,
    rep: str,
) -> list[dict]:
    """Per-drug raw / normalized / naive-baseline Spearman for one (heads, rep) configuration."""
    rows = []
    for drug, df in oof.items():
        c = fragility_index(truth_all, leave_out=drug)
        d = df.set_index("cell_line").join(c.rename("c"), how="inner").dropna()
        p, t, cc = d["y_pred"].to_numpy(), d["y_true"].to_numpy(), d["c"].to_numpy()
        rows.append(
            {
                "heads": heads,
                "rep": rep,
                "drug": drug,
                "n_lines": len(d),
                "rho_raw": spearmanr(p, t).statistic,
                "rho_normalized": _partial_spearman(p, t, cc),
                "rho_naive_baseline": spearmanr(cc, t).statistic,
            }
        )
    return rows


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    add_data_args(parser)
    parser.add_argument(
        "--drugs",
        nargs="+",
        default=None,
        help="Evaluation panel; default: the `selected` drugs in --panel-csv.",
    )
    parser.add_argument("--panel-csv", type=Path, default=DEFAULT_PANEL_CSV)
    parser.add_argument(
        "--heads",
        choices=["panel", "all", "both"],
        default="both",
        help="Train K=len(panel) heads, K=545, or both (default).",
    )
    parser.add_argument("--reps", nargs="+", default=list(REPS), choices=list(REPS))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=25,
                        help="Best epoch is median 6 over 36 recorded runs; 25 is never reached.")
    parser.add_argument(
        "--oof-csv",
        type=Path,
        default=None,
        help="Score existing out-of-fold predictions instead of refitting (no model fits at all). "
             "Expects rep, drug, cell_line, y_true, y_pred -- e.g. "
             "notebooks/outputs/panel/panel_oof_predictions.csv from notebook 14.",
    )
    parser.add_argument("--weighted", choices=["true", "false"], default=None,
                        help="With --oof-csv: pick one weighting arm if the file holds both.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    panel = args.drugs or load_panel(args.panel_csv)
    score = args.score  # from add_data_args; the training target

    truth_all = line_by_drug_truth(args.variant, args.data_root)
    all_drugs = list(truth_all.columns)
    missing = [d for d in panel if d not in all_drugs]
    if missing:
        raise ValueError(f"panel drugs not in the targets file: {missing}")

    print(f"fragility index over {len(all_drugs)} drugs, {len(truth_all)} lines")
    if args.oof_csv is None:
        print(f"panel ({len(panel)}): {panel}")

    if args.oof_csv is not None:
        # No training: the normalization is post-processing of predictions, and the fragility
        # index comes from labels alone.
        weighted = None if args.weighted is None else args.weighted == "true"
        configs = [(None, f"oof:{args.oof_csv.name}")]
        print(f"scoring existing predictions from {args.oof_csv} (0 model fits)\n")
    else:
        configs = []
        if args.heads in ("panel", "both"):
            configs.append((panel, f"K={len(panel)}"))
        if args.heads in ("all", "both"):
            configs.append((all_drugs, f"K={len(all_drugs)}"))
        print(f"score={score} variant={args.variant} seed={args.seed}")
        print(f"{len(configs) * len(args.reps)} runs x 5 folds "
              f"-- consider --oof-csv instead, which needs none\n")

    rows = []
    for drugs, heads in configs:
        for rep in args.reps:
            if args.oof_csv is not None:
                oof = load_oof_csv(args.oof_csv, rep=rep, weighted=weighted)
            else:
                oof = oof_line_predictions(
                    score,
                    drugs,
                    rep,
                    panel,
                    variant=args.variant,
                    data_root=args.data_root,
                    seed=args.seed,
                    epochs=args.epochs,
                )
            res = normalized_correlations(oof, truth_all, heads=heads, rep=rep)
            rows.extend(res)
            print(
                f'{heads:8s} {rep:8s} raw={np.mean([r["rho_raw"] for r in res]):+.3f}  '
                f'normalized={np.mean([r["rho_normalized"] for r in res]):+.3f}  '
                f'naive={np.mean([r["rho_naive_baseline"] for r in res]):+.3f}',
                flush=True,
            )

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
