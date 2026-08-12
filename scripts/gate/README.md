# `scripts/gate/` — pre-merge and post-merge checks

**Nothing here is pipeline code and nothing imports it**, the same exception
`check_resolved_paths.py` and `check_unbound_names.py` occupy one directory up. These run by hand
(or from a hook) around a merge.

| file | question it answers |
|---|---|
| `merge_gate.sh <branch>` | is this branch safe to merge into `main`? |
| `verify_main.sh` | is `main` **as merged** sound? A three-way merge can produce a tree neither parent had, so the gate cannot answer this. |
| `links.py` | does every relative markdown link and heading anchor resolve — in `.md` **and** in notebook markdown cells? |
| `artifacts.py` | does every backticked reference to `outputs/` or `figures/` point at something that exists? |
| `cmdpaths.py` | do the paths inside fenced command blocks and LaTeX `\code{}` resolve? |
| `coderefs.py` | which code references (`module.function`) named in prose no longer exist? Compared branch-vs-main, so only **new** breakage blocks. |

## The one rule these all share

**Every check prints how many candidates it examined, and has a floor that blocks when that count
collapses.** A check that passes because it examined nothing is the most common defect this project
has found in its own tooling — five instances, including one inside `check_unbound_names.py`, which
reported `checked 0 file(s) / no unbound names` and exited 0 when run from outside the repository.
`0 problems` and `0 examined` print identically unless the denominator is on screen.

## Known limits

Recorded here because a check is only as trustworthy as its documented ways of lying.

- **Two floors have never been fault-injected**, only reasoned about: `merge_gate.sh` check 1's
  merge-tree output size and check 6's diff file count. Every other floor has been verified by
  breaking the thing it guards.
- **The floors are round numbers**, chosen to sit obviously below the true count. They catch
  collapse to zero or near-zero. They do not catch a check that examines 80% of what it should.
- **`check_resolved_paths.py` sees `Path.glob` only** — a module-level `glob.glob` is invisible to
  it. That gap let `4a` call `glob.glob` unimported for weeks (`c351851`); `check_unbound_names.py`
  now covers it from the other side, by catching the unbound name rather than the dead pattern.
- **Branch-specific checks are not carried forward.** `merge_gate.sh` briefly held checks that
  pinned particular strings in `4a` and particular output directories against a rename. They were
  right for one merge and noise on every other. Add one for the merge that needs it, delete it
  after.
- **`links.py` already strips fenced blocks and inline code**, so cisplatin's SMILES
  `` `N[Pt](N)(Cl)Cl` `` is not read as a link. An ad-hoc regex written during a review *did* trip
  on it and report a false broken link — the lesson is about ad-hoc checks, not about this one.

## What they do not do

They are read-only. Nothing trains, reads the `h5ad` data, or writes outside `report/` and `/tmp`.
They check that the repository is internally consistent — never that a result is correct.
