"""Decisive single-run test of the first-fit hypothesis.

Section A trains `for rep in REPS: for alpha: for loss: for seed:`, so the first fit of the
process is X_pca / alpha=0 / mse / seed 42 -- precisely the one arm that moves across
executions (0.2473-0.2541 over the two replayable runs, everything else identical to six decimals;
   was 0.2450-0.2541 over eight until the band was rebuilt on 14.08.2026 -- six of those eight lived
   only in an agent scratch directory).

Reverse REPS and X_scGPT / alpha=0 / mse becomes first instead. That arm has read exactly
0.2009 in all five runs while it was NOT first. If it moves now, position causes the
instability, not the arm -- and a throwaway warm-up fit before the grid fixes it. If it
stays at 0.2009 and X_pca still wobbles, the hypothesis is wrong and the arm itself is
the problem.

Runs from a scratch copy; section A's artifacts are restored by the caller afterwards.
"""
import shutil
from pathlib import Path
import nbformat
from nbclient import NotebookClient
import pandas as pd


if __name__ == "__main__":   # guard added 14.08.2026 (Selin): importing this file used to RUN it

    REPO = Path('/Users/selin/PycharmProjects/OncoTox')
    TMP = Path('/Users/selin/.claude/jobs/ce8d4fe5/tmp')
    SECTION_A = [1, 3, 5, 7, 9, 11, 12, 14, 16, 17, 18]

    scratch = REPO / 'notebooks' / '_reverse_scratch.ipynb'
    shutil.copy2(REPO / 'notebooks/4a_percell_training.ipynb', scratch)
    nb = nbformat.read(scratch, as_version=4)

    OLD = "REPS = ['X_pca', 'X_scGPT']"
    NEW = "REPS = ['X_scGPT', 'X_pca']   # REVERSED for the first-fit test -- scratch copy only"
    assert nb.cells[1].source.count(OLD) == 1, 'REPS literal not found in cell 1'
    nb.cells[1].source = nb.cells[1].source.replace(OLD, NEW)
    print('REPS reversed in the scratch copy', flush=True)

    client = NotebookClient(nb, kernel_name='python3', timeout=None, record_timing=True,
                            allow_errors=False,
                            resources={'metadata': {'path': str(REPO / 'notebooks')}})
    with client.setup_kernel():
        for i in SECTION_A:
            print(f'  cell[{i}]', flush=True)
            client.execute_cell(nb.cells[i], i)

    lb = pd.read_csv(REPO / 'notebooks/outputs/panel/panel_leaderboard.csv')
    lb.to_csv(TMP / 'leaderboard_reversed.csv', index=False)
    scratch.unlink(missing_ok=True)
    print('\nREVERSED-ORDER RUN COMPLETE')
    print(lb.to_string())
