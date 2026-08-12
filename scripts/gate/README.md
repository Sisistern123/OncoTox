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

## A defect class no check here can catch

Every limitation above is about a check examining the wrong thing, or too little of it. This one is
different, and it is recorded because the reflex after reading the rest of this file is to believe
that more checks would have caught more.

**13.08.2026 — cell-level arrays passed into a line-level contract.** While fixing the head-bias
defect in `dreval_benchmark` (audit 08's fix had reached three training paths of four), the first
attempt passed the dataset's per-**cell** `y` and `mask` to `cv.per_drug_line_mean`, whose contract is
per-**line**. On a fixture with three lines of true value 0.2 / 0.9 / 0.9 and cell counts
1990 / 56 / 100:

    line-level (correct)   0.667
    cell-level (the bug)   0.251     — dominated by the one cell-rich line

**No test would have failed.** Both are finite floats in the right range; both look like a plausible
head bias; nothing raises, nothing is empty, no denominator collapses. A unit test asserting
"returns a float per drug" passes on the wrong one. It was caught by *reading the callee's docstring
before wiring the call* — the docstring says, in as many words, that a cell-level mean "would weigh a
1,990-cell line 35x a 56-cell line for the same single measurement".

The generalisation, such as it is: **checks catch the absent and the malformed; only reading catches
the plausible-but-wrong.** Where a function's correctness depends on which *level* its inputs are
aggregated at — cell versus line, fold versus pooled, drug versus panel — the contract has to be read,
because every level produces a number of the right type and magnitude.

(`check_unbound_names.py` did catch the *intermediate* state of that same edit — two names used before
being imported — in its first hour in the repository. That is the class it exists for, and it is a
different class from this one.)

## What they do not do

They are read-only. Nothing trains, reads the `h5ad` data, or writes outside `report/` and `/tmp`.
They check that the repository is internally consistent — never that a result is correct.
