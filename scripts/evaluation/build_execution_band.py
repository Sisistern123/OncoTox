"""Execution band: every arm x rep, over the independent executions of 4a section A that can be REPLAYED.

⚠️ **Rebuilt 14.08.2026 (Selin) from git-replayable executions only, and it used to be eight.**
Six of the original eight were read from ``/Users/selin/.claude/jobs/ce8d4fe5/tmp/`` -- an agent
session's scratch directory, never committed, one machine, transient by construction. So six eighths
of a band that ``docs/steps/05`` cites and that ``docs/OPEN_DECISIONS.md`` §3 leans on could not be
re-derived by anyone, including us once that directory is cleared.

**What the rebuild costs and what it buys.** The band narrows from 0.2450-0.2541 (width 0.0091) to
0.2473-0.2541 (width 0.0068), and "twelve of fourteen arm x rep rows identical" becomes thirteen of
fourteen -- fewer executions, so less opportunity to differ. **What it does not cost is §3's
argument**: the two values that argument turns on, 0.2541 ("no challenger wins") and 0.2473
("alpha=0/mae wins"), are exactly the two executions that remain. The verdict flip it demonstrates is
entirely inside the reproducible pair.

The eight-execution measurement is not deleted -- it happened -- it is recorded as superseded in
``docs/steps/corrections-and-dead-ends.md`` with its provenance limitation, and this file is now the
only thing that writes the live artifact.

⚠️ **Adding more executions is possible and is NOT done here.** Seven commits have touched
``panel_leaderboard.csv``, but only some are the same pipeline: the earlier ones predate the target
and panel corrections, so including them would measure pipeline changes rather than execution noise.
**Which commits count as comparable executions is an analysis decision** and was not taken --
parked as ``docs/OPEN_DECISIONS.md`` §8. The candidate set is one commit: ``ffe13be`` is the only
other execution on the corrected pipeline, and it carries 4 arms against these 7.

**21.10.2026 (Selin, options a and i): a third execution, and the band is now MLP-only.**
After the methods/results audit, 4a was re-executed in full. The audit changed only the ridge
control (fitted on each fold's fitting lines, grid ``logspace(-4, 6, 21)``), not the MLP path, so the
MLP arms of that run are a third execution of the same pipeline. The evidence that they are: 11 of
the 12 MLP arm x rep rows equal both 13.08 executions to four decimals. The ridge rows are dropped,
not carried. The 13.08 rows describe the pre-audit ridge the report no longer uses, and the new ridge
has only one recorded execution (``ridge_grid_refit.py`` and 4a both wrote the same file), which
would give a range of zero without measuring anything. Being a closed-form fit, it has no run-to-run
nondeterminism to measure.

The third execution was read from the working tree until Selin committed it as ``157a385``
(21.10.2026); it is now replayed from that commit like the other two.
"""
import io, subprocess
from pathlib import Path
import pandas as pd


if __name__ == "__main__":   # guard added 14.08.2026 (Selin): importing this file used to RUN it

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

    # All three replayed from committed history, so anyone with the repository can rebuild this file.
    frames = [
        tidy(from_git('9732b6f^'), 'exec_1_committed_13.08',  'normal order'),
        tidy(from_git('9732b6f'),  'exec_2_full_rerun_13.08', 'normal order'),
        tidy(from_git('157a385'),  'exec_3_audit_rerun_21.10', 'normal order'),
    ]

    band = pd.concat(frames, ignore_index=True)
    # MLP arms only since 21.10.2026 (Selin, option i): see the module docstring.
    band = band[band.arm.str.startswith('MLP')]
    assert band.groupby(['arm', 'rep']).execution.nunique().eq(3).all(), 'an MLP arm is missing an execution'

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
