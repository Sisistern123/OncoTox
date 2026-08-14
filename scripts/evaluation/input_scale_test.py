"""Does the 104x input-scale gap CAUSE the training-regime difference?

X_pca's input is refitted per fold inside cv.fold_pca_projections, so it cannot be rescaled
without changing the training path. X_scGPT's comes straight from obsm, so it can. The test
is therefore run in the other direction: scale X_scGPT UP to X_pca's magnitude and see
whether it adopts X_pca's training behaviour.

Prediction if scale is the driver: best_epoch collapses toward 1 (X_scGPT natively peaks at
7-8 in section A's alpha=0/mse arm) and the arm stops reproducing.
Prediction if it is not: nothing moves, because a uniform rescale changes no relative
structure whatsoever -- the same directions, the same ordering, the same everything but units.

Writes only to the job scratch directory; touches no committed artifact.
"""
import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path('/Users/selin/PycharmProjects/OncoTox')
sys.path.insert(0, str(ROOT))
from scripts.layout import PipelinePaths
from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS
from scripts.training.cv import oof_predictions
from scripts.training.training_utils import TrainConfig


if __name__ == "__main__":   # guard added 14.08.2026 (Selin): importing this file used to RUN it

    PANEL = pd.read_csv(ROOT / 'notebooks/outputs/panel/panel.csv')['drug_key'].tolist()
    paths = PipelinePaths.build(None, 'hvg5000', 'auc_cc')
    src = ad.read_h5ad(paths.targets_h5ad, backed='r')
    drugs = list(src.uns['ctrp_drugs'])
    k = [drugs.index(d) for d in PANEL]
    Y = np.asarray(src.obsm['Y_ctrp'], dtype=np.float32)[:, k]
    M = np.asarray(src.obsm['M_ctrp'], dtype=bool)[:, k]

    adata = ad.AnnData(obs=src.obs.copy())
    X = np.asarray(src.obsm['X_scGPT'], dtype=np.float32)
    adata.obsm['X_scGPT'] = X
    adata.obsm['Y_ctrp'] = np.where(M, Y, 0.0).astype(np.float32)
    adata.obsm['M_ctrp'] = M
    adata.uns['ctrp_drugs'] = PANEL
    src.file.close()

    # The measured gap, from diagnostics/input_scale.csv: 1.1062 / 0.0107.
    FACTOR = 1.1062 / 0.0107
    adata.obsm['X_scGPT_scaled'] = (X * FACTOR).astype(np.float32)
    print(f'scaling factor {FACTOR:.1f}x')
    for r in ('X_scGPT', 'X_scGPT_scaled'):
        print(f'  {r:16s} median per-dim sd = {np.median(adata.obsm[r].std(axis=0)):.4f}')

    cfg = TrainConfig(epochs=50, seed=42, loss='mse')
    rows = []
    for rep in ('X_scGPT', 'X_scGPT_scaled'):
        for seed in (42, 43, 44):
            pred, folds = oof_predictions(
                adata, rep, PANEL,
                config=replace(cfg, seed=seed), hidden_dims=DEFAULT_HIDDEN_DIMS['X_scGPT'],
                n_splits=5, density_weighting=False, alpha=0.0, init_head_bias=True,
                tag=f'scale_{rep}_s{seed}')
            for f in folds:
                rows.append({'rep': rep, 'seed': seed, 'fold': f['fold'],
                             'best_epoch': f['best_epoch'], 'best_val_obj': f['best_val_obj']})
            print(f'  done {rep} seed={seed}', flush=True)
            del pred

    d = pd.DataFrame(rows)
    d.to_csv('/Users/selin/.claude/jobs/ce8d4fe5/tmp/scale_test.csv', index=False)
    print('\n=== median best_epoch ===')
    print(d.groupby('rep')['best_epoch'].agg(['median', 'mean', 'min', 'max']).round(2).to_string())
    print('\n=== per fold/seed ===')
    print(d.pivot_table(index=['seed', 'fold'], columns='rep', values='best_epoch').to_string())
