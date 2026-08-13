"""Run 4a section A N times and record the spread ACROSS EXECUTIONS.

Section A is cells 1,3,5,7,9,11,12,14,16,17,18 (the code cells before section B). It is run
from a SCRATCH COPY of the committed notebook so the committed executed version is not
overwritten by a partial run, with the working directory still notebooks/ so every path
resolves identically.

Each run overwrites notebooks/outputs/panel/*. The leaderboard is captured after each, and
the committed artifacts are restored with `git checkout --` at the end -- so the repository
ends where it started plus one new file, the band itself.
"""
import shutil, sys, json
from pathlib import Path
import nbformat
from nbclient import NotebookClient
import pandas as pd

REPO = Path('/Users/selin/PycharmProjects/OncoTox')
TMP = Path('/Users/selin/.claude/jobs/ce8d4fe5/tmp')
SECTION_A = [1, 3, 5, 7, 9, 11, 12, 14, 16, 17, 18]
N = int(sys.argv[1])

captured = []
for run in range(1, N + 1):
    scratch = REPO / 'notebooks' / f'_band_scratch.ipynb'
    shutil.copy2(REPO / 'notebooks/4a_percell_training.ipynb', scratch)
    nb = nbformat.read(scratch, as_version=4)
    client = NotebookClient(nb, kernel_name='python3', timeout=None, record_timing=True,
                            allow_errors=False,
                            resources={'metadata': {'path': str(REPO / 'notebooks')}})
    print(f'=== run {run}/{N} start', flush=True)
    with client.setup_kernel():
        for i in SECTION_A:
            client.execute_cell(nb.cells[i], i)
    lb = pd.read_csv(REPO / 'notebooks/outputs/panel/panel_leaderboard.csv')
    lb.insert(0, 'execution', f'band_{run}')
    captured.append(lb)
    (TMP / f'leaderboard_band_{run}.csv').write_text(lb.to_csv(index=False))
    scratch.unlink(missing_ok=True)
    print(f'=== run {run}/{N} done', flush=True)

pd.concat(captured, ignore_index=True).to_csv(TMP / 'band_raw.csv', index=False)
print('BAND RUNS COMPLETE')
