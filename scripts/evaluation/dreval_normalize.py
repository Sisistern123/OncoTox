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

Usage (the training is the expensive part -- 5 folds per heads x rep combination):

    # panel = whatever `selected` marks in the learnability CSV
    uv run scripts/evaluation/dreval_normalize.py

    # explicit panel, e.g. the literature-anchored one (docs/steps/05)
    uv run scripts/evaluation/dreval_normalize.py --drugs methotrexate dasatinib paclitaxel \\
        vincristine afatinib topotecan tanespimycin selumetinib

    # panel only, skip the K=545 run (halves the compute)
    uv run scripts/evaluation/dreval_normalize.py --heads panel

Reproduces the evaluation protocol of ``notebooks/11_auc_vs_aucz.ipynb`` exactly: 5-fold GroupKFold
over the train+val cell lines (the 27 test lines are never touched), each fold predicts only the lines
it did not see, per-cell predictions averaged to one value per cell line, and every model scored
against the common curve-fit ``auc`` ranking regardless of which score it trained on.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold
from torch.utils.data import DataLoader

from scripts.model.dataset import MultiDrugDataset
from scripts.model.OncoMLP import OncoMLP
from scripts.preprocessing.layout import PipelinePaths, add_data_args
from scripts.training.train_multitask import DEFAULT_HIDDEN_DIMS
from scripts.training.training_utils import TrainConfig, train_model

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


def oof_line_predictions(
    score: str,
    drugs: list[str],
    rep: str,
    eval_drugs: list[str],
    *,
    variant: str,
    data_root: Path | None,
    seed: int = 42,
    epochs: int = 50,
) -> dict[str, pd.DataFrame]:
    """Train K=len(drugs) heads and return out-of-fold per-cell-line predictions.

    Mirrors ``oof_run`` in ``notebooks/11_auc_vs_aucz.ipynb``, with one addition: the cell-line
    identity is kept, because the normalization needs to align drugs and the fragility index by line.

    :returns: drug -> DataFrame(cell_line, y_true, y_pred); ``y_true`` is always on the curve-fit
        ``auc`` scale so models trained on different scores share one yardstick.
    """
    truth = sc.read_h5ad(PipelinePaths.build(data_root, variant, "auc").targets_h5ad, backed="r")
    truth_y = np.asarray(truth.obsm["Y_ctrp"], dtype=float)
    truth_m = np.asarray(truth.obsm["M_ctrp"], dtype=bool)
    truth_drugs = list(truth.uns["ctrp_drugs"])

    adata = sc.read_h5ad(PipelinePaths.build(data_root, variant, score).targets_h5ad)
    eligible = adata.obs["split_ctrp"].isin(["train", "val"]).to_numpy()  # test never touched
    groups = adata.obs["Cell_line"].astype(str).to_numpy()
    idx = np.flatnonzero(eligible)
    mask = np.asarray(adata.obsm["M_ctrp"], dtype=bool)
    all_drugs = list(adata.uns["ctrp_drugs"])

    cfg = TrainConfig(epochs=epochs, seed=seed, log_every=1000)
    pred = np.full((adata.n_obs, len(drugs)), np.nan)
    for fold, (tr, va) in enumerate(GroupKFold(n_splits=5).split(idx, groups=groups[idx]), 1):
        trc = np.zeros(adata.n_obs, bool)
        trc[idx[tr]] = True
        vac = np.zeros(adata.n_obs, bool)
        vac[idx[va]] = True
        tr_ds = MultiDrugDataset(adata=adata, use_rep=rep, cell_mask=trc, drugs=drugs)
        va_ds = MultiDrugDataset(adata=adata, use_rep=rep, cell_mask=vac, drugs=drugs)
        model = OncoMLP(
            input_dim=tr_ds.X.shape[1],
            hidden_dims=DEFAULT_HIDDEN_DIMS[rep],
            dropout_rate=0.5,
            input_dropout=0.1,
            norm="layer",
            output_dim=len(tr_ds.drug_names),
        )
        best, _ = train_model(
            model,
            DataLoader(tr_ds, batch_size=128, shuffle=True),
            DataLoader(va_ds, batch_size=256, shuffle=False),
            config=cfg,
            tag=f"norm_{score}_{rep}_K{len(drugs)}_s{seed}_f{fold}",
            drug_names=tr_ds.drug_names,
        )
        best = best.to("cpu").eval()  # train_model leaves it on mps/cuda
        with torch.no_grad():
            x = torch.tensor(np.asarray(adata.obsm[rep], dtype=np.float32)[vac])
            pred[vac] = best(x).numpy()  # held-out lines only

    out: dict[str, pd.DataFrame] = {}
    for drug in eval_drugs:
        j, k, kt = drugs.index(drug), all_drugs.index(drug), truth_drugs.index(drug)
        lines, y_true, y_pred = [], [], []
        for line in np.unique(groups[eligible]):
            ci = np.flatnonzero((groups == line) & eligible)
            obs = ci[mask[ci, k] & truth_m[ci, kt] & np.isfinite(pred[ci, j])]
            if obs.size:
                lines.append(line)
                y_pred.append(pred[obs, j].mean())  # cells -> one value per cell line
                y_true.append(truth_y[obs, kt].mean())
        out[drug] = pd.DataFrame({"cell_line": lines, "y_true": y_true, "y_pred": y_pred})
    return out


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
    parser.add_argument("--epochs", type=int, default=50)
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

    configs = []
    if args.heads in ("panel", "both"):
        configs.append((panel, f"K={len(panel)}"))
    if args.heads in ("all", "both"):
        configs.append((all_drugs, f"K={len(all_drugs)}"))

    print(f"panel ({len(panel)}): {panel}")
    print(f"score={score} variant={args.variant} seed={args.seed}")
    print(f"fragility index over {len(all_drugs)} drugs, {len(truth_all)} lines")
    print(f"{len(configs) * len(args.reps)} runs x 5 folds\n")

    rows = []
    for drugs, heads in configs:
        for rep in args.reps:
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
