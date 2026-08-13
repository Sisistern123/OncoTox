"""Score the SAME out-of-fold predictions under all four aggregation conventions.

The open decision (docs/OPEN_DECISIONS.md 2) is two independent binary choices:

  level : per CELL, or per CELL LINE (cells averaged to their line first)
  folds : POOLED into one correlation, or per fold then averaged

This computes all four and changes nothing else -- same predictions, same drugs, same
held-out lines, same Spearman. It picks no winner; which convention becomes canonical is
Selin's, and the point of the table is to show how much the choice moves the answer.

Per-cell predictions come from runs/percell/ (gitignored, 110 MB, written by 4a section A
cell 14). Line-level truth and the fold map come from panel_oof_predictions.csv, joined on
(drug, cell_line): a pair is scored only if it appears there, so the observed-label mask is
inherited rather than re-derived.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PERCELL = Path('runs/percell')
OUT = Path('notebooks/outputs/panel')

cell_index = pd.read_csv(PERCELL / 'cell_index.csv')
drugs = json.load(open(PERCELL / 'drug_order.json'))
oof = pd.read_csv(OUT / 'panel_oof_predictions.csv')

# cell -> row position, and cell_line per row, in the .npy row order
line_of_cell = cell_index['cell_line'].to_numpy()

def rho(x, y):
    """Spearman, or NaN when it is undefined rather than 0."""
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan
    return spearmanr(x, y).statistic

def score(arm_oof, P):
    """Return the four aggregations for one arm-seed, each a mean over the 11 drugs."""
    out = {}
    line_rows, cell_rows = [], []
    for j, drug in enumerate(drugs):
        sub = arm_oof[arm_oof.drug == drug]
        if sub.empty:
            continue
        truth = dict(zip(sub.cell_line, sub.y_true))
        foldof = dict(zip(sub.cell_line, sub.fold))

        # --- line level: the cells of a line averaged, then one point per line ---
        keep = sub.cell_line.to_numpy()
        pred_line = {}
        for cl in keep:
            m = line_of_cell == cl
            v = P[m, j]
            v = v[~np.isnan(v)]
            if v.size:
                pred_line[cl] = float(v.mean())
        cls = [c for c in keep if c in pred_line]
        if cls:
            line_rows.append(dict(
                drug=drug,
                pooled=rho(np.array([truth[c] for c in cls]),
                           np.array([pred_line[c] for c in cls])),
                per_fold=np.nanmean([
                    rho(np.array([truth[c] for c in cls if foldof[c] == f]),
                        np.array([pred_line[c] for c in cls if foldof[c] == f]))
                    for f in sorted({foldof[c] for c in cls})]),
            ))

        # --- cell level: every cell of a held-out line is its own point ---
        sel = np.isin(line_of_cell, keep)
        yc = np.array([truth.get(cl, np.nan) for cl in line_of_cell[sel]])
        pc = P[sel, j]
        fc = np.array([foldof.get(cl, -1) for cl in line_of_cell[sel]])
        ok = ~np.isnan(yc) & ~np.isnan(pc)
        if ok.sum() >= 3:
            cell_rows.append(dict(
                drug=drug,
                pooled=rho(yc[ok], pc[ok]),
                per_fold=np.nanmean([rho(yc[ok & (fc == f)], pc[ok & (fc == f)])
                                     for f in sorted(set(fc[ok]))]),
            ))

    for name, rows in (('line', line_rows), ('cell', cell_rows)):
        d = pd.DataFrame(rows)
        out[f'{name}_pooled'] = d.pooled.mean() if len(d) else np.nan
        out[f'{name}_per_fold'] = d.per_fold.mean() if len(d) else np.nan
    return out

PAT = re.compile(r'^percell_(oof|mil)_(X_pca|X_scGPT)_a([\d.]+)(?:_(mae|mse))?_s(\d+)\.npy$')
records = []
for f in sorted(PERCELL.glob('*.npy')):
    m = PAT.match(f.name)
    if not m:
        continue                      # the two seedless legacy files
    kind, rep, alpha, loss, seed = m.groups()
    model = 'mil' if kind == 'mil' else 'mlp'
    loss = loss or 'mse'              # mil files carry no loss token; they are mse
    alpha, seed = float(alpha), int(seed)
    sel = oof[(oof.rep == rep) & (oof.model == model) & (oof.alpha == alpha)
              & (oof.loss == loss) & (oof.seed == seed)]
    if sel.empty:
        print(f'  skip {f.name} -- no matching rows in panel_oof_predictions.csv')
        continue
    P = np.load(f)
    records.append({'rep': rep, 'model': model, 'alpha': alpha, 'loss': loss,
                    'seed': seed, **score(sel, P)})
    print(f'  {f.name}', flush=True)

assert records, 'no per-cell arrays matched the pattern -- refusing to write an empty comparison'
df = pd.DataFrame(records)
df.to_csv(OUT / 'panel_aggregation_comparison.csv', index=False)
print(f'\nwrote panel_aggregation_comparison.csv -- {len(df)} arm-seeds')

mean = df.groupby(['model', 'rep', 'alpha', 'loss'])[
    ['line_pooled', 'line_per_fold', 'cell_pooled', 'cell_per_fold']].mean().round(4)
print('\nmean over the three seeds:')
print(mean.to_string())
