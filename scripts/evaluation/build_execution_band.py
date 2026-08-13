"""Full execution band: every arm x rep, over all eight independent executions of 4a section A."""
import io, subprocess
from pathlib import Path
import pandas as pd

TMP = Path('/Users/selin/.claude/jobs/ce8d4fe5/tmp')

def tidy(df, label, cond):
    df = df.copy(); df.columns = [str(c) for c in df.columns]
    first = df.columns[0]; reps = df.iloc[0]; body = df.iloc[2:]
    rows = []
    for col in df.columns[1:]:
        if not col.startswith('spearman'):
            continue
        for _, r in body.iterrows():
            arm = r[first]
            if isinstance(arm, str) and (arm.startswith('MLP') or arm.startswith('ridge')):
                rows.append({'execution': label, 'condition': cond,
                             'arm': arm, 'rep': reps[col], 'spearman': float(r[col])})
    return pd.DataFrame(rows)

def from_git(rev):
    return pd.read_csv(io.StringIO(subprocess.run(
        ['git','show',f'{rev}:notebooks/outputs/panel/panel_leaderboard.csv'],
        capture_output=True, text=True).stdout))

frames = [
    tidy(from_git('9732b6f^'), 'exec_1_committed_13.08',  'normal order'),
    tidy(from_git('9732b6f'),  'exec_2_full_rerun_13.08', 'normal order'),
]
for i, k in enumerate([1, 2, 3], start=3):
    frames.append(tidy(pd.read_csv(TMP/f'pre_warmup_leaderboard_band_{k}.csv').drop(columns=['execution']),
                       f'exec_{i}_band_14.08', 'normal order'))
frames.append(tidy(pd.read_csv(TMP/'leaderboard_reversed.csv'), 'exec_6_reversed_14.08', 'REPS reversed'))
for i, k in enumerate([1, 2], start=7):
    frames.append(tidy(pd.read_csv(TMP/f'leaderboard_band_{k}.csv').drop(columns=['execution']),
                       f'exec_{i}_warmup_14.08', 'warm-up active'))

band = pd.concat(frames, ignore_index=True)
wide = band.pivot_table(index=['arm','rep'], columns='execution', values='spearman')
wide['min'], wide['max'] = wide.min(axis=1), wide.max(axis=1)
wide['range'] = wide['max'] - wide['min']
wide = wide.sort_values('range', ascending=False)

out = Path('notebooks/outputs/panel/panel_execution_band.csv')
wide.round(6).to_csv(out)
print(f'{out.name}: {len(wide)} arm x rep rows over {band.execution.nunique()} executions\n')
print(wide[['min','max','range']].round(4).to_string())
print(f"\narms with range > 0     : {(wide['range'] > 1e-9).sum()} of {len(wide)}")
print(f"arms with range > 0.001 : {(wide['range'] > 0.001).sum()} of {len(wide)}")
