"""Test the reading offered for why the bag objective helps X_scGPT and hurts X_pca.

The account: within-line variation is noise under a per-cell objective (every cell of a line is
asked to predict that line's single label), X_scGPT carries more of it, and pooling into a bag
removes exactly that. PREDICTION: the bag objective's advantage should be larger where within-line
variation is larger.

Test, on committed artifacts only. Within each drug, split its held-out cell lines at the median of
their within-line prediction spread, then compare the MIL-minus-per-cell gap in the high half against
the low half. If the account is right the gap is larger in the high half.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


if __name__ == "__main__":   # guard added 14.08.2026 (Selin): importing this file used to RUN it

    ROOT = Path('/Users/selin/PycharmProjects/OncoTox')
    P = ROOT / 'notebooks/outputs'
    A, L = 0.5, 'mse'                       # the only arm MIL was run at

    spread = pd.read_csv(P / 'panel/panel_within_line_spread.csv')
    spread = spread[(spread.alpha == A) & (spread.loss == L)]
    sd = spread.groupby(['rep', 'drug', 'cell_line'])['within_line_sd'].mean().reset_index()

    mlp = pd.read_csv(P / 'panel/panel_oof_predictions.csv')
    mlp = mlp[(mlp.alpha == A) & (mlp.loss == L) & (mlp.model == 'mlp')]
    mil = pd.read_csv(P / 'mil/mil_oof_predictions.csv')
    mil = mil[(mil.alpha == A) & (mil.loss == L)]

    def per_drug_rho(df, keep):
        """Mean per-drug Spearman over the (drug, line) pairs in `keep`."""
        d = df.merge(keep, on=['rep', 'drug', 'cell_line'], how='inner')
        out = {}
        for (rep, drug), g in d.groupby(['rep', 'drug']):
            g = g.groupby('cell_line')[['y_true', 'y_pred']].mean()
            if len(g) >= 5 and g.y_true.nunique() > 2:
                out[(rep, drug)] = spearmanr(g.y_true, g.y_pred).statistic
        return pd.Series(out)

    rows = []
    for rep in ['X_pca', 'X_scGPT']:
        s = sd[sd.rep == rep]
        for drug, g in s.groupby('drug'):
            med = g.within_line_sd.median()
            for half, sel in [('low', g[g.within_line_sd <= med]), ('high', g[g.within_line_sd > med])]:
                rows.append(sel.assign(half=half))
    halves = pd.concat(rows, ignore_index=True)[['rep', 'drug', 'cell_line', 'half']]

    print(f'{halves.groupby(["rep","half"]).size().to_dict()}\n')
    res = []
    for half in ['low', 'high']:
        keep = halves[halves.half == half][['rep', 'drug', 'cell_line']]
        a, b = per_drug_rho(mlp, keep), per_drug_rho(mil, keep)
        common = a.index.intersection(b.index)
        gap = (b[common] - a[common])
        for rep in ['X_pca', 'X_scGPT']:
            g = gap[[i for i in gap.index if i[0] == rep]]
            res.append({'half': half, 'rep': rep, 'n_drugs': len(g),
                        'per_cell': a[[i for i in common if i[0] == rep]].mean(),
                        'bag': b[[i for i in common if i[0] == rep]].mean(),
                        'bag_minus_percell': g.mean()})
    r = pd.DataFrame(res)
    print('=== mean per-drug Spearman, by within-line-spread half ===')
    print(r.round(4).to_string(index=False))
    print()
    for rep in ['X_pca', 'X_scGPT']:
        lo = r[(r.rep == rep) & (r.half == 'low')].bag_minus_percell.iloc[0]
        hi = r[(r.rep == rep) & (r.half == 'high')].bag_minus_percell.iloc[0]
        print(f'  {rep:8s} bag advantage: low-spread half {lo:+.4f}   high-spread half {hi:+.4f}   '
              f'change {hi-lo:+.4f}')
    print('\nPREDICTION: the bag advantage should be LARGER in the high-spread half.')
    r.to_csv(P / 'diagnostics/mil_by_within_line_spread.csv', index=False)
    print(f'wrote diagnostics/mil_by_within_line_spread.csv')
