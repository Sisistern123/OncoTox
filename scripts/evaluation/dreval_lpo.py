"""DrEval's **leave-pairs-out** benchmark — the mode in which their bias baselines are informative.

**The question this answers.** Is the model learning drug-specific biology, or is it exploiting the
fact that some cell lines are fragile to everything and some compounds are potent against everything?

**Why the existing LCO run cannot answer it, and this one can.** The project's benchmark runs
leave-**cell-line**-out. There a held-out line has never been seen, so ``NaiveCellLineMeanPredictor``
has no line mean to predict with: it emits a constant and scores a Spearman of exactly **0.0000**, and
DrEval's normalization has no line effect to subtract either. The line-bias channel is therefore
invisible under LCO **by construction**, not by oversight.

Under **LPO** the held-out pairs come from lines that *are* in training. ``NaiveCellLineMeanPredictor``
can estimate each line's mean, and its score **is** the line-fragility share — measured as an honest
predictor evaluated on held-out data, with no post-hoc subtraction and no diagnostic that reads
held-out labels. The same argument makes ``NaiveTissueMeanPredictor`` the direct test of this
project's motivating hypothesis: if lineage is what the fitted baseline encodes, a tissue-mean
predictor should already capture much of the achievable signal.

⚠️ **LPO IS AN EASIER PROTOCOL AND ITS NUMBERS ARE NOT COMPARABLE TO THE LCO ONES.** Every held-out
pair here belongs to a cell line the model has seen in other pairs. Nothing produced by this script
may be quoted beside a leave-cell-line-out figure, and nothing here replaces the project's headline
protocol — LCO remains the design the results are reported under, because predicting an *unseen* line
is the task the project is about. This run exists to decompose where the signal comes from, not to
report a better score.

⚠️ **Our own model is not run here.** Its out-of-fold predictions are cross-validated by cell line, so
it has no LPO predictions and scoring it would require retraining under a different split. The
comparison drawn is therefore between DrEval's own baselines: what a bias-only predictor achieves
against what a feature-using predictor achieves, on identical folds.

Written 14.08.2026 (Selin: *"incorporate it and run"*), after establishing that ``drevalpy`` ships
four test modes and six naive baselines while this project had only ever run one mode.

**Data wiring.** The ``DrugResponseDataset`` is built exactly as
``notebooks/analysis/evaluation/dreval_benchmark.ipynb`` §1 builds it — per (line, drug) mean of
``auc_cc``, cell-line features being the mean of that line's cell embeddings — with one addition, the
``tissue`` view the tissue baselines require, taken from ``obs['Cancer_type']``. ⚠️ That construction
is currently written twice, here and in the notebook; the notebook should import
:func:`build_response_dataset` rather than repeat it, and until it does the two can drift.

Run:  .venv/bin/python scripts/evaluation/dreval_lpo.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import anndata as ad  # noqa: E402
from drevalpy.datasets.dataset import DrugResponseDataset, FeatureDataset  # noqa: E402
from drevalpy.evaluation import evaluate  # noqa: E402
from drevalpy.models import MODEL_FACTORY  # noqa: E402

from scripts.evaluation.dreval_normalize import load_panel  # noqa: E402
from scripts.layout import DEFAULT_CTRP_SCORE, PipelinePaths  # noqa: E402

N_SPLITS = 5
SEED = 42
OUT = ROOT / "notebooks/outputs/dreval/dreval_lpo_results.csv"

#: The six naive baselines drevalpy ships, plus two feature-using models as a reference ceiling.
#: The naive ones are the point: each isolates one channel a model could exploit without learning
#: anything drug-specific about a cell line's biology.
MODELS = [
    "NaivePredictor",                 # global mean — the floor
    "NaiveDrugMeanPredictor",         # compound potency only
    "NaiveCellLineMeanPredictor",     # LINE FRAGILITY only — degenerate under LCO, informative here
    "NaiveTissueMeanPredictor",       # lineage only — the motivating hypothesis, as a number
    "NaiveTissueDrugMeanPredictor",   # lineage x compound
    "NaiveMeanEffectsPredictor",      # line + drug effects, no interaction
]

#: Their per-drug reference models, fitted on OUR embedding as the cell-line view -- the same trick
#: the LCO notebook uses (``hp['cell_line_views'] = view``), because drevalpy's defaults expect a
#: view named ``gene_expression`` and ours are named ``pca`` / ``scgpt``.
SINGLE_DRUG = [("SingleDrugElasticNet", "pca"), ("SingleDrugRandomForest", "pca"),
               ("SingleDrugElasticNet", "scgpt"), ("SingleDrugRandomForest", "scgpt")]


def per_drug_spearman(y_true, y_pred, drug_ids) -> float:
    """Mean Spearman **within each drug**, across the held-out lines — *this project's* metric.

    DrEval reports a pooled correlation over all (line, drug) pairs, in which compound potency
    dominates: a predictor that knows only each drug's mean already scores ~0.74 there. Scoring
    within a drug removes that by construction, which is why the project's headline numbers are on
    this scale and not on the pooled one. Reported alongside so a bias baseline can be compared with
    the model's own 0.2824 rather than against a number on a different scale.
    """
    from scipy.stats import spearmanr
    out = []
    for d in np.unique(drug_ids):
        m = drug_ids == d
        if np.unique(y_true[m]).size > 2 and np.unique(y_pred[m]).size > 1:
            out.append(spearmanr(y_true[m], y_pred[m]).statistic)
    return float(np.nanmean(out)) if out else float("nan")


def build_response_dataset():
    """``(response, cell_features, drug_features)`` over the panel — the notebook's §1, plus tissue."""
    panel = load_panel()
    paths = PipelinePaths.build(None, "hvg5000", DEFAULT_CTRP_SCORE)
    a = ad.read_h5ad(paths.targets_h5ad)
    Y = np.asarray(a.obsm["Y_ctrp"], float)
    M = np.asarray(a.obsm["M_ctrp"], bool)
    all_drugs = list(a.uns["ctrp_drugs"])
    lines = a.obs["Cell_line"].astype(str).to_numpy()
    tissue = a.obs["Cancer_type"].astype(str).to_numpy()
    emb = {"scgpt": np.asarray(a.obsm["X_scGPT"], float),
           "pca": np.asarray(a.obsm["X_pca"], float)}

    missing = [d for d in panel if d not in all_drugs]
    if missing:
        raise ValueError(f"panel drugs absent from the targets file: {missing}")

    pairs, feats = [], {}
    for line in np.unique(lines):
        ci = np.flatnonzero(lines == line)
        feats[line] = {
            "cell_line_name": np.array([line]),
            "tissue": np.array([tissue[ci][0]]),
            "scgpt": emb["scgpt"][ci].mean(0),
            "pca": emb["pca"][ci].mean(0),
        }
        for d in panel:
            k = all_drugs.index(d)
            obs = ci[M[ci, k]]
            if obs.size:
                pairs.append((line, d, float(Y[obs, k].mean())))

    response = DrugResponseDataset(
        response=np.array([p[2] for p in pairs], dtype=float),
        cell_line_ids=np.array([p[0] for p in pairs]),
        drug_ids=np.array([p[1] for p in pairs]),
        dataset_name="OncoTox_SCP542_CTRPv2",
    )
    return (response,
            FeatureDataset(features=feats),
            FeatureDataset(features={d: {"pubchem_id": np.array([d])} for d in panel}),
            panel)


def main() -> None:
    response, cell_feats, drug_feats, panel = build_response_dataset()
    print(f"{len(response.response)} (cell line x drug) pairs | "
          f"{len(set(response.cell_line_ids))} cell lines | {len(panel)} drugs")

    response.split_dataset(n_cv_splits=N_SPLITS, mode="LPO", split_validation=True,
                           random_state=SEED)
    folds = response.cv_splits
    print(f"{len(folds)} LPO folds\n")

    rows = []

    def record(name, i, te, pred):
        pred = np.asarray(pred, dtype=float)
        scored = DrugResponseDataset(response=te.response, cell_line_ids=te.cell_line_ids,
                                     drug_ids=te.drug_ids, predictions=pred)
        met = evaluate(scored, metric=["Pearson", "Spearman", "R^2", "MSE"])
        # n_scored is not decoration. The naive baselines are scored on every test pair; the
        # per-drug models are scored only where prediction succeeded (the ``ok`` mask below), so
        # without this column a comparison between the two groups cannot be checked for being
        # like-for-like. Added 14.08.2026 before the LPO ordering was quoted anywhere.
        rows.append({"algorithm": name, "fold": i, "n_scored": int(len(pred)), **met,
                     "per_drug_Spearman": per_drug_spearman(te.response, pred, te.drug_ids)})

    for name in MODELS:
        if name not in MODEL_FACTORY:
            print(f"  skip {name} — not in this drevalpy build")
            continue
        for i, fold in enumerate(folds, start=1):
            tr, te = fold["train"], fold["test"]
            model = MODEL_FACTORY[name]()
            hp = MODEL_FACTORY[name].get_hyperparameter_set()[0].copy()
            try:
                model.build_model(hp)
                model.train(output=tr, cell_line_input=cell_feats, drug_input=drug_feats)
                pred = model.predict(cell_line_ids=te.cell_line_ids, drug_ids=te.drug_ids,
                                     cell_line_input=cell_feats, drug_input=drug_feats)
            except Exception as exc:                      # a model this data cannot feed
                print(f"  {name} fold {i}: FAILED — {type(exc).__name__}: {exc}"[:150])
                continue
            record(name, i, te, pred)
        done = [r for r in rows if r["algorithm"] == name]
        if done:
            print(f"  {name:30s} Spearman {np.mean([r['Spearman'] for r in done]):+.4f}  "
                  f"over {len(done)} fold(s)")

    # Their per-drug references, on our views. Fitted one drug at a time, as the model expects.
    for name, view in SINGLE_DRUG:
        label = f"{name} ({view})"
        for i, fold in enumerate(folds, start=1):
            tr, te = fold["train"], fold["test"]
            pred = np.full(len(te.response), np.nan)
            for d in np.unique(te.drug_ids):
                trm, tem = tr.drug_ids == d, te.drug_ids == d
                if trm.sum() < 5 or tem.sum() == 0:
                    continue
                model = MODEL_FACTORY[name]()
                hp = MODEL_FACTORY[name].get_hyperparameter_set()[0].copy()
                hp["cell_line_views"] = view
                try:
                    model.build_model(hp)
                    model.train(
                        output=DrugResponseDataset(response=tr.response[trm],
                                                   cell_line_ids=tr.cell_line_ids[trm],
                                                   drug_ids=tr.drug_ids[trm]),
                        cell_line_input=cell_feats, drug_input=drug_feats)
                    pred[tem] = model.predict(cell_line_ids=te.cell_line_ids[tem],
                                              drug_ids=te.drug_ids[tem],
                                              cell_line_input=cell_feats, drug_input=drug_feats)
                except Exception as exc:
                    print(f"  {label} fold {i} drug {d}: {type(exc).__name__}: {exc}"[:120])
            ok = ~np.isnan(pred)
            if ok.sum() == 0:
                continue
            sub = DrugResponseDataset(response=te.response[ok], cell_line_ids=te.cell_line_ids[ok],
                                      drug_ids=te.drug_ids[ok])
            record(label, i, sub, pred[ok])
        done = [r for r in rows if r["algorithm"] == label]
        if done:
            print(f"  {label:30s} Spearman {np.mean([r['Spearman'] for r in done]):+.4f}  "
                  f"per-drug {np.mean([r['per_drug_Spearman'] for r in done]):+.4f}")

    if not rows:
        raise SystemExit("no model produced a score — refusing to write an empty artifact")
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(df)} rows)\n")

    summary = df.groupby("algorithm")[["Spearman", "per_drug_Spearman", "Pearson", "R^2"]].mean().round(4)
    summary = summary.sort_values("Spearman", ascending=False)
    print("=== DrEval LPO, mean over folds ===")
    print("  Spearman = pooled over all pairs (drug potency dominates it);"
          " per_drug_Spearman = within-drug, this project's metric\n")
    print(summary.to_string())


if __name__ == "__main__":
    main()
