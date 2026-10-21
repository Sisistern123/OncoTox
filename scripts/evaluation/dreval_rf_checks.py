"""Why DrEval's random forest fails on PCA: two checks on the LCO benchmark (Selin, 21.10.2026).

**The observation.** On DrEval's LCO protocol, with the 153 CV lines and the per-fold PCA,
``SingleDrugRF (pca)`` scored 0.0458 on the normalised metric against 0.2988 for ``SingleDrugRF
(scgpt)`` and 0.3304 for the MLP on PCA (``dreval/dreval_lco_results.csv``). Linear models on the
same PCA line means do not fail: DrEval's elastic net scored 0.2218 there, the ridge 0.2725 in the
pipeline's CV.

**Two readings, each tested by changing one setting** (both readings are hypotheses, not results):

1. *Starved trees.* drevalpy's first hyperparameter set gives ``max_samples=0.2`` and ``max_depth=5``,
   and sklearn's default considers all features at every split -- so each tree sees ~23 of ~115
   training lines while choosing among 512 features. Check (ii): ``max_samples=1.0``, everything else
   unchanged. Run on both representations.
2. *Spread signal.* PCA components are uncorrelated by construction, so a response signal that is a
   weighted sum of many components is carried by none of them alone, which axis-aligned splits on few
   lines cannot find. Check (i): the forest on the first k components only, k in {10, 50, 512}
   (Selin, 21.10.2026; the values are arbitrary: few, moderate, all). PCA only: scGPT's dimensions
   have no order, so "the first k" would mean something else there (Selin, 21.10.2026, option a).

**Everything else is the reported setup**, so the k=512 / default row must reproduce ``SingleDrugRF
(pca)`` and the scGPT default row ``SingleDrugRF (scgpt)`` in ``dreval_lco_results.csv`` exactly:
the same response dataset (:func:`dreval_lpo.build_response_dataset`, 153 CV lines), drevalpy's LCO
split with the notebook's arguments, the PCA fitted per fold on each fold's ``train`` lines, line
means as features, forests seeded with 42, and DrEval's normalisation against
``NaiveMeanEffectsPredictor`` fitted on the fold's training pairs.

Run:  .venv/bin/python scripts/evaluation/dreval_rf_checks.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from drevalpy.datasets.dataset import DrugResponseDataset, FeatureDataset  # noqa: E402
from drevalpy.evaluation import evaluate  # noqa: E402
from drevalpy.models import MODEL_FACTORY  # noqa: E402

import anndata as ad  # noqa: E402
from scripts.evaluation.dreval_lpo import build_response_dataset  # noqa: E402
from scripts.layout import DEFAULT_CTRP_SCORE, PipelinePaths  # noqa: E402
from scripts.training.cv import fold_pca_projections  # noqa: E402

SEED = 42
OUT = ROOT / "notebooks/outputs/dreval/dreval_rf_checks.csv"

#: (view, k components or None for all, max_samples or None for drevalpy's default)
ARMS = [("pca", 10, None), ("pca", 50, None), ("pca", 512, None), ("pca", 512, 1.0),
        ("scgpt", None, None), ("scgpt", None, 1.0)]


def run_forest(tr, te, feats, drug_feats, view, max_samples):
    """drevalpy's SingleDrugRandomForest, one model per drug, as the LCO notebook fits it."""
    name = "SingleDrugRandomForest"
    hp = MODEL_FACTORY[name].get_hyperparameter_set()[0].copy()
    hp["cell_line_views"] = view
    preds = np.full(len(te.response), np.nan)
    for d in np.unique(te.drug_ids):
        trm, tem = tr.drug_ids == d, te.drug_ids == d
        m = MODEL_FACTORY[name]()
        m.build_model(hp)
        m.model.set_params(random_state=SEED)
        if max_samples is not None:
            m.model.set_params(max_samples=max_samples)
        m.train(output=DrugResponseDataset(response=tr.response[trm], cell_line_ids=tr.cell_line_ids[trm],
                                           drug_ids=tr.drug_ids[trm]),
                cell_line_input=feats, drug_input=drug_feats)
        preds[tem] = m.predict(cell_line_ids=te.cell_line_ids[tem], drug_ids=te.drug_ids[tem],
                               cell_line_input=feats, drug_input=drug_feats)
    return preds


def main() -> None:
    response, base_feats, drug_feats, _panel = build_response_dataset()
    response.split_dataset(n_cv_splits=5, mode="LCO", split_validation=True,
                           validation_ratio=0.1, random_state=SEED)
    folds = response.cv_splits

    paths = PipelinePaths.build(None, "hvg5000", DEFAULT_CTRP_SCORE)
    lines = ad.read_h5ad(paths.targets_h5ad, backed="r").obs["Cell_line"].astype(str).to_numpy()
    fit_masks = []
    for sp in folds:
        fit_lines = set(sp["train"].cell_line_ids)
        assert not fit_lines & (set(sp["validation"].cell_line_ids) | set(sp["test"].cell_line_ids))
        fit_masks.append(np.isin(lines, list(fit_lines)))
    proj = fold_pca_projections(paths.raw_h5ad, fit_masks, seed=SEED)

    rows = []
    for fold, sp in enumerate(folds, start=1):
        tr, te = sp["train"], sp["test"]
        feats = FeatureDataset(features={
            l: {**base_feats.features[l], "pca": proj[fold - 1][lines == l].mean(0)}
            for l in base_feats.features})
        naive = MODEL_FACTORY["NaiveMeanEffectsPredictor"]()
        naive.build_model({})
        naive.train(output=tr, cell_line_input=feats, drug_input=drug_feats)
        ref = naive.predict(cell_line_ids=te.cell_line_ids, drug_ids=te.drug_ids,
                            cell_line_input=feats, drug_input=drug_feats)
        for view, k, max_samples in ARMS:
            f = feats
            if k is not None and k < proj[fold - 1].shape[1]:
                f = FeatureDataset(features={l: {**feats.features[l], "pca": feats.features[l]["pca"][:k]}
                                             for l in feats.features})
            p = run_forest(tr, te, f, drug_feats, view, max_samples)
            norm = evaluate(DrugResponseDataset(response=te.response - ref, cell_line_ids=te.cell_line_ids,
                                                drug_ids=te.drug_ids, predictions=p - ref), ["Spearman"])
            rows.append({"fold": fold, "view": view, "k": k or 512,
                         "max_samples": max_samples or 0.2, "Spearman: normalized": norm["Spearman"]})
        print(f"fold {fold} done", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(OUT, index=False)
    print(res.groupby(["view", "k", "max_samples"])["Spearman: normalized"].agg(["mean", "std"]).round(4))


if __name__ == "__main__":
    main()
