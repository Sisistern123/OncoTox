"""Section E2: the 'smaller study' scenario -- shrink the cell lines AND refit PCA on only those cells.

Section E thinned the LABEL supply while keeping every cell in the fold, in the batches and in the
per-fold PCA, so both arms saw an identical input at every point. That refuted the small-data reading
of X_pca's lead. The explanation it left standing is the INPUT side: X_pca is fitted on this atlas and
X_scGPT is not, so PCA's directions adapt to exactly these cells while scGPT's are frozen.

E2 tests that directly. Cell lines are dropped ENTIRELY -- their cells leave the study -- so the
per-fold PCA is refitted on a smaller atlas while X_scGPT's frozen embedding is unaffected in kind.

  If PCA's advantage comes from adaptation, its margin should SHRINK as the atlas shrinks.
  If it comes from what variance-maximisation captures regardless, the margin should hold.

The x-axis is read from the fold log, as in section E, rather than from the constants below.

⚠️ **Corrected 14.08.2026 — this paragraph described an implementation the file does not have.** It
read *"No code change: subsetting the AnnData is the scenario."* Subsetting **was** the first
approach and it does not work: the per-fold PCA refits from the full counts h5ad, so the masks handed
to it must stay full-length, and a subset AnnData breaks that correspondence with an ``IndexError``.
The dropped lines are instead marked **ineligible** by overwriting ``split_ctrp`` (see the loop
below), which removes them from the folds, from the batches and from the PCA's fitting set alike.
**The scenario is unchanged** — those cells leave the study either way — but it is achieved by
relabelling rather than by subsetting, and *there is* a code change. The string written is
``'test'``, used purely as an ineligibility marker because ``eligible_splits`` is
``("train", "val")``; **the project's real held-out test lines are not read here or anywhere else.**
``base_split`` is snapshotted with ``.copy()`` before the loop precisely because the overwrite is
in place, so each iteration rebuilds from the original assignment rather than from the previous
iteration's.
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
full = ad.AnnData(obs=src.obs.copy())
for r in ('X_pca', 'X_scGPT'):
    full.obsm[r] = np.asarray(src.obsm[r], dtype=np.float32)
full.obsm['Y_ctrp'], full.obsm['M_ctrp'] = np.where(M, Y, 0.0).astype(np.float32), M
full.uns['ctrp_drugs'] = PANEL
src.file.close()

base_split = full.obs['split_ctrp'].astype(str).to_numpy().copy()
elig = full.obs['split_ctrp'].isin(['train', 'val']).to_numpy()
lines = np.unique(full.obs['Cell_line'].astype(str).to_numpy()[elig])
print(f'eligible cell lines: {len(lines)}')

SIZES = [31, 62, 94, None]     # None = the full eligible set
cfg = TrainConfig(epochs=50, seed=42, loss='mse')
rows, folds_log = [], []
for n_lines in SIZES:
    for seed in (42, 43, 44):
        rng = np.random.default_rng(seed)
        keep = lines if n_lines is None else np.sort(rng.choice(lines, size=n_lines, replace=False))
        # The per-fold PCA refits from the FULL counts h5ad, so its masks must stay full-length --
        # subsetting the AnnData breaks that correspondence (IndexError, 14.08.2026). Instead the
        # dropped lines are marked INELIGIBLE, which removes them from the folds, from the batches
        # and from the PCA's fitting set alike. That is the E2 scenario: the study shrinks.
        sub = full
        sub.obs['split_ctrp'] = np.where(
            np.isin(full.obs['Cell_line'].astype(str).to_numpy(), keep),
            base_split, 'test')
        for rep in ('X_pca', 'X_scGPT'):
            pred, folds = oof_predictions(
                sub, rep, PANEL, config=replace(cfg, seed=seed), hidden_dims=(),
                n_splits=5, density_weighting=False, alpha=0.0, init_head_bias=True,
                counts_h5ad=paths.raw_h5ad, pca_seed=42,
                tag=f'E2_{n_lines or "all"}_{rep}_s{seed}')
            g = line_level_predictions(pred, sub, PANEL, folds=folds)
            rho = [spearmanr(d.y_true, d.y_pred).statistic
                   for _, d in g.groupby('drug') if d.y_true.nunique() > 2]
            rows.append({'n_lines_kept': len(keep), 'rep': rep, 'seed': seed,
                         'n_fit_lines': int(np.median([f['n_fit_lines'] for f in folds])),
                         'order': float(np.nanmean(rho))})
            print(f"  {len(keep):3d} lines {rep:8s} s{seed}: {rows[-1]['order']:.4f}", flush=True)
            del pred, g

d = pd.DataFrame(rows)
d.to_csv(ROOT/'notebooks/outputs/panel/panel_curve_e2.csv', index=False)
m = d.groupby(['n_lines_kept', 'rep'])['order'].mean().unstack()
m['margin'] = m.X_pca - m.X_scGPT
m['median_n_fit_lines'] = d.groupby('n_lines_kept')['n_fit_lines'].median()
print('\n=== E2: study shrunk, PCA refitted on the smaller atlas ===')
print(m.round(4).to_string())
print('\n  section E (labels thinned, atlas intact): +0.0036 / +0.0090 / +0.0505 / +0.0317')
