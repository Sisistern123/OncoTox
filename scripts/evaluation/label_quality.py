"""Does label quality explain per-drug performance?

The 14.08 limitation said the label-noise ceiling was 'not currently measured'. That is true of
replicate DISAGREEMENT -- the live target folds it into one CurveCurator fit -- but the fit itself
carries quality columns (R2, RMSE, pValue, conc_pts_fit) that were never read. This reads them, on the
exact (cell line, drug) pairs this project trains on.
"""
import sys
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path('/Users/selin/PycharmProjects/OncoTox')
sys.path.insert(0, str(ROOT))
from scripts.preprocessing.ctrp_to_h5ad import _normalize_cell_line, _normalize_drug
from scripts.layout import PipelinePaths


if __name__ == "__main__":   # guard added 14.08.2026 (Selin): importing this file used to RUN it

    p = PipelinePaths.build(None, 'hvg5000', 'auc_cc')
    resp = pd.read_csv(p.ctrp_response_csv,
                       usecols=['ccl_name', 'drug_name', 'AUC_curvecurator', 'R2', 'RMSE',
                                'pValue', 'conc_pts_fit'])
    resp['line'] = _normalize_cell_line(resp.ccl_name)
    resp['drug'] = _normalize_drug(resp.drug_name)

    panel = pd.read_csv(ROOT / 'notebooks/outputs/panel/panel.csv')
    drugs = set(_normalize_drug(panel.drug_ctrp))
    split = pd.read_csv(ROOT / 'splits/split_ctrp.csv')
    lines = set(_normalize_cell_line(split.Cell_line.str.split('_').str[0]))

    r = resp[resp.drug.isin(drugs) & resp.line.isin(lines)].dropna(subset=['AUC_curvecurator'])
    print(f'curves behind this project\'s labels: {len(r):,}  '
          f'({r.drug.nunique()} drugs x {r.line.nunique()} lines)\n')

    print('=== 1 · how good are the fits the labels come from? ===')
    for c in ['R2', 'RMSE', 'pValue', 'conc_pts_fit']:
        q = r[c].quantile([.05, .25, .5, .75, .95])
        print(f'  {c:13s} median {q[.5]:8.4f}   IQR {q[.25]:8.4f}-{q[.75]:8.4f}   5-95% {q[.05]:.4f}-{q[.95]:.4f}')
    print(f"  R2 < 0.5      : {100*(r.R2 < 0.5).mean():5.1f} % of curves")
    print(f"  R2 < 0.25     : {100*(r.R2 < 0.25).mean():5.1f} %")
    print(f"  pValue > 0.05 : {100*(r.pValue > 0.05).mean():5.1f} %  (fit not significant)")

    print('\n=== 2 · do worse-fitting drugs get predicted worse? ===')
    q = (r.groupby('drug').agg(median_R2=('R2', 'median'),
                               frac_R2_lt_half=('R2', lambda v: (v < 0.5).mean()),
                               frac_ns=('pValue', lambda v: (v > 0.05).mean()),
                               n=('R2', 'size')))
    corr = pd.read_csv(ROOT / 'notebooks/outputs/panel/panel_per_drug_correlation.csv')
    perf = (corr[(corr.loss == 'mse') & (corr.alpha == 0.0)]
            .groupby(['rep', 'drug'])['spearman'].mean().unstack(0))
    j = q.join(perf, how='inner').sort_values('median_R2')
    print(j.round(4).to_string())
    print()
    for rep in ['X_pca', 'X_scGPT']:
        for metric, sense in [('median_R2', 'higher = better fit'), ('frac_ns', 'higher = worse fit')]:
            s = spearmanr(j[metric], j[rep])
            print(f'  {rep:8s} vs {metric:16s} ({sense:19s}): rho={s.statistic:+.3f}  p={s.pvalue:.3f}  n={len(j)}')
    j.to_csv(ROOT / 'notebooks/outputs/diagnostics/label_quality_vs_performance.csv')
    print(f'\nwrote diagnostics/label_quality_vs_performance.csv')
