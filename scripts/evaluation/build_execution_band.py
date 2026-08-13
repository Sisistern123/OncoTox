"""Assemble the re-execution band from five independent runs of 4a section A."""
import io, subprocess
from pathlib import Path
import pandas as pd

TMP = Path('/Users/selin/.claude/jobs/ce8d4fe5/tmp')

def tidy(df, label):
    """panel_leaderboard.csv is a two-row-header pivot; flatten it to long form."""
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    first = df.columns[0]
    reps = df.iloc[0]                      # row 0 holds the rep per column
    body = df.iloc[2:]                     # rows 0,1 are the header rows
    out = []
    for col in df.columns[1:]:
        if not col.startswith('spearman'):
            continue
        rep = reps[col]
        for _, r in body.iterrows():
            arm = r[first]
            if not isinstance(arm, str) or not arm.startswith('MLP'):
                continue
            out.append({'execution': label, 'arm': arm, 'rep': rep,
                        'spearman': float(r[col])})
    return pd.DataFrame(out)

def from_git(rev):
    txt = subprocess.run(['git', 'show', f'{rev}:notebooks/outputs/panel/panel_leaderboard.csv'],
                         capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(txt))

frames = [
    tidy(from_git('9732b6f^'), 'exec_1_committed_13.08'),
    tidy(from_git('9732b6f'),  'exec_2_full_rerun_13.08'),
]
for k in (1, 2, 3):
    d = pd.read_csv(TMP / f'leaderboard_band_{k}.csv').drop(columns=['execution'])
    frames.append(tidy(d, f'exec_{k + 2}_band_14.08'))

band = pd.concat(frames, ignore_index=True)
wide = band.pivot_table(index=['arm', 'rep'], columns='execution', values='spearman')
wide['min'], wide['max'] = wide.min(axis=1), wide.max(axis=1)
wide['range'] = wide['max'] - wide['min']
wide = wide.sort_values('range', ascending=False)

out = Path('notebooks/outputs/panel/panel_execution_band.csv')
wide.round(6).to_csv(out)
print(f'wrote {out.name} -- {len(wide)} arm x rep rows over {band.execution.nunique()} executions\n')
print(wide.round(4).to_string())
print(f'\nlargest range over any arm : {wide["range"].max():.4f}')
print(f'median range               : {wide["range"].median():.4f}')
print(f'arms whose range exceeds 0.01: {(wide["range"] > 0.01).sum()} of {len(wide)}')
