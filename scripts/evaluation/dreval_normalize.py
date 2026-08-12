"""DrEval's normalized evaluation of our own out-of-fold predictions — their recipe, nothing added.

**What the normalization is.** A model can score a good correlation by learning only that some drugs are
potent and some cell lines are fragile, without predicting anything drug-specific. DrEval controls for
that by subtracting a naive baseline — ``overall mean + cell-line effect + drug effect`` — from the truth
*and* from the prediction, and re-scoring what is left:

    y_true_norm = y_true - y_naive        y_pred_norm = y_pred - y_naive

    — DrEval, ``drevalpy/visualization/utils.py::_normalize_metrics_by_mean_effects``, applied per
    algorithm and per CV split. Correlation metrics only; DrEval exclude MSE/RMSE/MAE from the
    normalized set, and so does this module.

**Nothing here re-implements it.** The baseline is drevalpy's own
``MODEL_FACTORY["NaiveMeanEffectsPredictor"]`` and the metrics are ``drevalpy.evaluation.evaluate``, so
this module is wiring: it reads our out-of-fold predictions, hands them to DrEval's code in the shape
that code expects, and writes the result out. If DrEval change the metric, we inherit the change.

**Scope, deliberately.** This file was **archived on 12.08.2026 and restored paper-only** the same day
(Selin): *"it should be as contained as the paper itself in functionality, and nothing new."* Its
previous version additionally removed the **cell-line effect** using held-out labels — a diagnostic that
asked how much of a per-drug correlation was mere line fragility — and scored every model against one
common ``auc`` ranking regardless of what it trained on. Both were **local inventions with no
counterpart in the paper**, and both were deleted rather than archived, so that a home-grown metric
cannot be revived from inside a file named after DrEval without being re-decided
(``scripts/archive/README.md``; recoverable at ``git show bf93084:scripts/evaluation/dreval_normalize.py``).
Whether that question gets answered again, and how, is for **audit 11 (Evaluation)**.

**⚠️ In our split design this metric does not control for cell-line fragility, and that is not a defect
in it.** DrEval evaluate in several test modes; ours is leave-cell-line-out, where every held-out line
is unseen by the baseline, so ``cell_line_effects.get(cl, 0)`` returns **0** and the subtraction removes
the *drug* effect only. Demonstrated, not assumed: a synthetic predictor emitting nothing but
``overall mean + line effect + drug effect`` — zero drug-specific signal — scores **normalized Spearman
0.98** here, pooled and per drug. The line effect survives in truth and prediction alike and correlates
with itself.

So a high normalized score is **not** evidence of drug-specific biology under leave-cell-line-out. It is
the correct behaviour for a prediction benchmark, because no model can know an unseen line's effect —
but the question *"is this drug-specific signal or general fragility?"* is left open by it, and reading
this output as though it were answered is the mistake to avoid. Answering it needs held-out labels and
is therefore a diagnostic rather than a metric, not something to fold back into this file.
**Audit 11** decides how to answer it.

**The fold column is required, and this is not pedantry.** The naive baseline must be fitted on the
folds a prediction did *not* come from; fitting it on the same out-of-fold rows it will be subtracted
from would let held-out labels define the baseline they are scored against. ⚠️ The current writer,
``notebooks/3_panel_training.ipynb``, does **not** emit a fold column
(``drug, cell_line, y_true, y_pred, rep, weighted``), so this script cannot run on the predictions in
the tree. That is a change owed to notebook 3 at R5 of the sweep, not something to work around here.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from drevalpy.datasets.dataset import DrugResponseDataset, FeatureDataset
from drevalpy.datasets.utils import CELL_LINE_IDENTIFIER, DRUG_IDENTIFIER
from drevalpy.evaluation import evaluate
from drevalpy.models import MODEL_FACTORY

#: The drug panel, rebuilt 12.08.2026 on FDA approval and published determinants
#: (``notebooks/drug_selection/literature_panel.ipynb``). Until then this pointed at
#: ``outputs/learnability/ctrp_drug_learnability_auc.csv``, the artifact of the discredited kill/spare
#: gate -- so the default panel was one whose selection criterion had been retracted.
DEFAULT_PANEL_CSV = Path("notebooks/outputs/panel/panel.csv")
DEFAULT_OUT = Path("notebooks/outputs/dreval/dreval_normalized.csv")

#: DrEval compute the normalized metrics for correlations only, because subtracting a fitted baseline
#: from both sides changes the scale of the residuals and makes an absolute error uninterpretable.
NORMALIZED_METRICS = ["Pearson", "Spearman", "Kendall", "R^2"]
RAW_METRICS = ["Pearson", "Spearman", "Kendall", "R^2", "MSE", "RMSE", "MAE"]

REQUIRED_COLUMNS = ("drug", "cell_line", "fold", "y_true", "y_pred")


def load_panel(panel_csv: Path = DEFAULT_PANEL_CSV) -> list[str]:
    """Read the drug panel, in the spelling the response data uses.

    ``panel.csv`` carries both ``drug`` (the FDA list's name, e.g. *Cisplatin*) and ``drug_key`` (what
    CTRPv2 calls it, e.g. ``platin``). Predictions are keyed by the latter, so that is what is returned.
    """
    panel = pd.read_csv(panel_csv)
    if "drug_key" not in panel.columns:
        raise ValueError(
            f"{panel_csv} has no 'drug_key' column. Expected the output of "
            f"notebooks/drug_selection/literature_panel.ipynb -- pass --drugs explicitly instead."
        )
    return panel["drug_key"].tolist()


def load_oof(path: Path, *, rep: str, weighted: bool | None = None) -> pd.DataFrame:
    """Read one arm of the out-of-fold predictions and check it can carry the normalization.

    :param rep: representation to keep, ``X_pca`` or ``X_scGPT``.
    :param weighted: keep only rows with this ``weighted`` flag; ``None`` keeps all.
    :raises ValueError: if the ``fold`` column is missing, which makes the baseline unfittable without
        leaking held-out labels into it (see the module docstring).
    """
    oof = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in oof.columns]
    if missing:
        raise ValueError(
            f"{path} is missing {missing}. Present: {list(oof.columns)}. "
            f"A 'fold' column is required: the naive baseline has to be fitted on the folds a "
            f"prediction did not come from."
        )
    oof = oof[oof["rep"] == rep]
    if weighted is not None:
        oof = oof[oof["weighted"] == weighted]
    if oof.empty:
        raise ValueError(f"No rows in {path} for rep={rep!r}, weighted={weighted!r}.")
    return oof.reset_index(drop=True)


def _identity_features(ids: np.ndarray, view: str) -> FeatureDataset:
    """A ``FeatureDataset`` mapping each identifier to itself.

    ``NaiveMeanEffectsPredictor.train`` takes feature datasets but reads only the identifiers out of
    them -- its prediction is ``mean + cell_line_effect + drug_effect`` and uses no features at all. An
    identity view is therefore the honest way to call their code without inventing features it would
    ignore.
    """
    return FeatureDataset(features={i: {view: np.array([i])} for i in np.unique(ids)})


def naive_predictions(oof: pd.DataFrame) -> pd.Series:
    """DrEval's ``NaiveMeanEffectsPredictor``, fitted per fold on the *other* folds.

    Returns one prediction per row of ``oof``, aligned to its index. A cell line or drug unseen in the
    training folds gets effect 0, which is drevalpy's own behaviour, not a choice made here.
    """
    out = pd.Series(np.nan, index=oof.index, dtype=float)
    for fold in oof["fold"].unique():
        held_out = oof["fold"] == fold
        train, test = oof[~held_out], oof[held_out]

        model = MODEL_FACTORY["NaiveMeanEffectsPredictor"]()
        model.train(
            output=DrugResponseDataset(
                response=train["y_true"].to_numpy(),
                cell_line_ids=train["cell_line"].to_numpy(),
                drug_ids=train["drug"].to_numpy(),
            ),
            cell_line_input=_identity_features(train["cell_line"].to_numpy(), CELL_LINE_IDENTIFIER),
            drug_input=_identity_features(train["drug"].to_numpy(), DRUG_IDENTIFIER),
        )
        out.loc[test.index] = model.predict(
            cell_line_ids=test["cell_line"].to_numpy(),
            drug_ids=test["drug"].to_numpy(),
            cell_line_input=_identity_features(test["cell_line"].to_numpy(), CELL_LINE_IDENTIFIER),
            drug_input=_identity_features(test["drug"].to_numpy(), DRUG_IDENTIFIER),
        )
    return out


def _score(frame: pd.DataFrame, metrics: list[str], *, normalized: bool) -> dict[str, float]:
    """Metrics for one group, through drevalpy's own :func:`evaluate`."""
    truth, pred = frame["y_true"], frame["y_pred"]
    if normalized:
        truth, pred = truth - frame["y_naive"], pred - frame["y_naive"]
    return evaluate(
        dataset=DrugResponseDataset(
            response=truth.to_numpy(),
            predictions=pred.to_numpy(),
            cell_line_ids=frame["cell_line"].to_numpy(),
            drug_ids=frame["drug"].to_numpy(),
        ),
        metric=metrics,
    )


def normalized_evaluation(oof: pd.DataFrame) -> pd.DataFrame:
    """Raw and normalized metrics, pooled over all pairs and per drug.

    Both groupings are reported because they answer different questions and DrEval report both. Pooled,
    the normalization removes the fact that drugs differ in potency and lines in fragility, so what is
    left is genuine (line × drug) specificity. Per drug, the drug effect is already constant and only
    the cell-line effect is removed, so it asks whether the ranking within a drug is more than a
    restatement of which lines are broadly sensitive.
    """
    scored = oof.assign(y_naive=naive_predictions(oof))
    rows = [{
        "grouping": "pooled", "drug": "", "n": len(scored),
        **{f"{k}": v for k, v in _score(scored, RAW_METRICS, normalized=False).items()},
        **{f"{k}: normalized": v
           for k, v in _score(scored, NORMALIZED_METRICS, normalized=True).items()},
    }]
    for drug, group in scored.groupby("drug", sort=False):
        rows.append({
            "grouping": "per_drug", "drug": drug, "n": len(group),
            **{f"{k}": v for k, v in _score(group, RAW_METRICS, normalized=False).items()},
            **{f"{k}: normalized": v
               for k, v in _score(group, NORMALIZED_METRICS, normalized=True).items()},
        })
    return pd.DataFrame(rows)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--oof-csv", type=Path, required=True,
                        help="Out-of-fold predictions; needs drug, cell_line, fold, y_true, y_pred.")
    parser.add_argument("--rep", default="X_scGPT", choices=["X_pca", "X_scGPT"])
    parser.add_argument("--weighted", type=lambda s: s.lower() == "true", default=None,
                        help="Keep only rows with this weighted flag; omit to keep all.")
    parser.add_argument("--panel-csv", type=Path, default=DEFAULT_PANEL_CSV)
    parser.add_argument("--drugs", nargs="+", default=None,
                        help=f"Drugs to score; default: the panel in --panel-csv.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    drugs = args.drugs or load_panel(args.panel_csv)
    oof = load_oof(args.oof_csv, rep=args.rep, weighted=args.weighted)

    absent = sorted(set(drugs) - set(oof["drug"]))
    if absent:
        raise ValueError(f"{len(absent)} panel drugs have no predictions in {args.oof_csv}: {absent}")
    oof = oof[oof["drug"].isin(drugs)]

    results = normalized_evaluation(oof)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(f"{len(drugs)} drugs | {oof.fold.nunique()} folds | rep={args.rep} -> {args.out}")
    print(results.set_index(["grouping", "drug"])[["Spearman", "Spearman: normalized"]].round(3))


if __name__ == "__main__":
    main()
