"""Item 8B: does the input_dropout asymmetry between the arms affect the Q1 margin?

input_dropout=0.1 zeroes each input coordinate independently. scGPT's dimensions are entangled and
comparable in magnitude; PCA's are variance-ordered, so the same rate removes a heavier-tailed share
for PCA -- measured 14.08.2026 as a 9.4x per-dimension sd spread against scGPT's 5.6x, with PC1 at
5.9% of retained variance. The item asks whether that ASYMMETRY matters, which was never tested.

It matters more now than when the item was written: weight_decay is 0.0, so dropout is the ONLY
regularizer.

Design: section C's linear row (alpha=0, mse, linear head, 3 seeds x 5 folds), both representations,
at input_dropout 0.1 (as shipped) and 0.0 (off). 60 fits. If the margin is unchanged the asymmetry is
harmless; if it moves, every Q1 margin carries an uncontrolled regularizer difference.
"""
import sys
from dataclasses import replace
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path('/Users/selin/PycharmProjects/OncoTox'); sys.path.insert(0, str(ROOT))
from scripts.layout import PipelinePaths
from scripts.training.cv import oof_predictions, line_level_predictions
from scripts.training.training_utils import TrainConfig

PANEL = pd.read_csv(ROOT/'notebooks/outputs/panel/panel.csv')['drug_key'].tolist()
paths = PipelinePaths.build(None, 'hvg5000', 'auc_cc')
src = ad.read_h5ad(paths.targets_h5ad, backed='r')
drugs = list(src.uns['ctrp_drugs']); k = [drugs.index(d) for d in PANEL]
Y = np.asarray(src.obsm['Y_ctrp'], dtype=np.float32)[:, k]
M = np.asarray(src.obsm['M_ctrp'], dtype=bool)[:, k]
adata = ad.AnnData(obs=src.obs.copy())
for r in ('X_pca', 'X_scGPT'):
    adata.obsm[r] = np.asarray(src.obsm[r], dtype=np.float32)
adata.obsm['Y_ctrp'], adata.obsm['M_ctrp'] = np.where(M, Y, 0.0).astype(np.float32), M
adata.uns['ctrp_drugs'] = PANEL
src.file.close()

cfg = TrainConfig(epochs=50, seed=42, loss='mse')
rows = []
for idrop in (0.1, 0.0):
    for rep in ('X_pca', 'X_scGPT'):
        for seed in (42, 43, 44):
            pred, folds = oof_predictions(
                adata, rep, PANEL, config=replace(cfg, seed=seed), hidden_dims=(),
                n_splits=5, density_weighting=False, alpha=0.0, init_head_bias=True,
                input_dropout=idrop, counts_h5ad=paths.raw_h5ad, pca_seed=42,
                tag=f'drop{idrop}_{rep}_s{seed}')
            g = line_level_predictions(pred, adata, PANEL, folds=folds)
            rho = [spearmanr(d.y_true, d.y_pred).statistic
                   for _, d in g.groupby('drug') if d.y_true.nunique() > 2]
            rows.append({'input_dropout': idrop, 'rep': rep, 'seed': seed,
                         'order': float(np.nanmean(rho))})
            print(f'  {idrop} {rep} s{seed}: {rows[-1]["order"]:.4f}', flush=True)
            del pred

d = pd.DataFrame(rows)
d.to_csv(ROOT/'notebooks/outputs/diagnostics/input_dropout_test.csv', index=False)
m = d.groupby(['input_dropout', 'rep'])['order'].mean().unstack()
m['margin'] = m.X_pca - m.X_scGPT
print('\n=== mean per-drug Spearman, linear head, alpha=0, mse, 3 seeds ===')
print(m.round(4).to_string())
print(f"\n  Q1 margin at input_dropout=0.1 (shipped): {m.loc[0.1,'margin']:+.4f}")
print(f"  Q1 margin at input_dropout=0.0 (off)    : {m.loc[0.0,'margin']:+.4f}")
print(f"  change: {m.loc[0.0,'margin']-m.loc[0.1,'margin']:+.4f}")
