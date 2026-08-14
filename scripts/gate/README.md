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

**But the denominator is only half of it, and the sixth instance is what shows why (13.08.2026).** An
ad-hoc sweep for hardcoded absolute paths annotated each candidate with
`os.path.exists(path.rstrip('/*'))`, intending "strip a trailing glob before testing". `rstrip` takes
a *character set*, not a suffix, and strips only from the end — so a path ending in `.csv` lost
nothing, the `*` stayed mid-string where `os.path.exists` has no glob semantics, and **every**
glob-containing path was reported dead. It named `drug_coverage` c13's `runs/*_all_drugs/` glob as
matching nothing. It matches 17 files.

Its denominator was correct and on screen throughout: 24 examined, 23 resolved, one dead. The first
five instances were silent zeros in the **denominator**; this one is a silent false **positive in a
single row**, which no floor can see. Three things kept it alive, and they generalise better than the
mechanism does: it was a **minority verdict** (a check that fails everything gets debugged; one that
fails a single thing gets believed), it was **corroborated by a true fact** (`runs/` really is
gitignored and untracked, so the wrong conclusion had real evidence beside it), and it **agreed with
what the sweep was looking for**, so it drew less scrutiny than a surprising answer would have. The
class does not come from checks that are wrong everywhere; it comes from checks that are wrong only
where nobody looks, and it survives because the false verdict is the one you were expecting.

So the mitigation is a **pair**, not one habit: *report the denominator, and keep a case the check
must pass.* A check needs an input it is required to accept, not only inputs it is required to
reject — the fault-injected floors below have both, and the two that are not fault-injected have
neither. Second: use the matching operation rather than a proxy (a path containing `*` is a glob;
`glob.glob` answers the question, `os.path.exists` answers a different one), and print the operation,
not the verdict — had the output shown the literal string tested, the `*` would have been visible.

*How it was found is a rule about routing, not about authors.* It surfaced only because the finding
was sent back to its author **with a question attached** — investigate this, do not repoint it —
which forced a look at the artifact instead of a re-read of the sweep's own output. Had the finding
been accepted, it would have stood. Route a finding back with a question; never with an acceptance.

## `check_calls.py` — three checks for defect classes nothing here was catching (13.08.2026)

Added at Gate 4 of the pipeline review. Each generalises a defect found **by hand** during it, on the
principle that a defect found once by reading will be reintroduced unless something looks for it.

| check | generalises | what it does |
|---|---|---|
| **signatures** | — | every call into `scripts.*` against the real signature: unexpected keywords, too many positionals, missing required arguments |
| **preconditions** | `cv.oof_predictions`'s `counts_h5ad`, which **no caller passed**, so R4 would have died on its first arm | enumerates parameters that *have* a default but whose value makes the function `raise`. A signature check cannot see these — the call is well-formed and the failure is hours into a run |
| **producers** | `drug_catalog` §5, which read a path §3 never wrote, so it raised on any machine | every file a notebook reads must be written by something, or be a known raw source |

**The preconditions check enumerates; it does not pass or fail.** Whether a caller satisfies a
precondition depends on the *values* it passes, which are not static. The output is a list to read.
Treating it as a gate would be claiming an assurance it cannot give.

**The producers check reports orphans rather than failing on them.** It follows string literals only,
so a path built through a variable — `out = DIR / name` then `df.to_csv(out)` — reads as an orphan
when it is not. It did exactly that on its first run, and the finding was verified by hand and left
in as a reminder that the check is a prompt, not a verdict.

## Two checks added 14.08.2026, against the defects that actually recurred

`retired_values.py` — **a registry of values that have been corrected, flagged wherever they reappear
as live text.** The most repeated defect here is not a wrong number but a *corrected* number surviving
where the correction did not reach: `Q2_CONTROL_THRESHOLD` in five places after it was closed, the
`0.98` synthetic-predictor figure cleared from one of three named locations, `4,704` gene counts, `545`
standing in for the head count, `180` trainable lines. Each was found by hand, late. On its first run
this check found **28 live occurrences in 15 files**, seven of them genuine — four being the same
`4,704` a docs-only sweep had "fixed" days earlier without touching `scripts/` or `notebooks/`.
Registering a value is meant to be part of retiring it.

`shell_safety.py` — **backticks inside a `$( ... )` substitution in these shell scripts.** There a
backtick opens a nested substitution, so a markdown-style code span in a comment makes the script fail
at runtime; `bash -n` does not catch it and the script still prints most of its output. It reached
committed `main` on 14.08.2026. ⚠️ Its first two versions were themselves broken — one summed
`grep -c` with `awk -F:` and reported 0 against a true 11 (the silent-zero failure this directory
exists to prevent, inside the check written to prevent it), the other used a regex that skipped
exactly the multi-line substitutions the defect lives in. It is now a bracket scanner **verified by
fault injection**.

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
- ⛔ **`verify_main.sh` was not read-only until 14.08.2026, and its own limits said it was.** Its
  module-import check imported every module under `scripts/` except `archive` and `gate`. **Nine
  files in `scripts/evaluation/` have no `if __name__ == "__main__"` guard**, so importing them ran
  them: `aggregation_comparison.py` and `build_execution_band.py` write committed artifacts at top
  level, and `section_e2_smaller_study.py` and `input_dropout_test.py` call `oof_predictions` — they
  **train**. Running the post-merge gate therefore retrained models and left the tree dirty. Fixed by
  excluding `evaluation`, exactly as `gate` was already excluded for the same reason.
  **Residual limitation, open:** those nine modules are no longer import-checked, so a syntax error
  or missing import in `scripts/evaluation/` is invisible to this gate. The real fix is a `__main__`
  guard on each, after which the exclusion can be dropped — nine files, not done in a consolidation
  pass. **What that costs meanwhile:** every other `scripts/` module is still covered (24 imported,
  0 failed, 1.7 s), and the evaluation scripts are exercised whenever a notebook runs them, so the
  gap is narrow but real.
- ⚠️ **Both gates invoke every check through `uv run`, which blocks while another `uv` process holds
  the project environment (found 14.08.2026).** With a `uv run jupyter lab` server and a live kernel
  running, `verify_main.sh` sat at 0 % CPU for over twenty minutes without producing a line, and had
  to be killed. It is not a hang in the checks themselves: run directly against `.venv/bin/python`,
  the same checks complete in seconds. **Consequence for anyone verifying main:** if the gate appears
  to stall, it is waiting on the environment lock, not working — reproduce the individual checks with
  `.venv/bin/python` rather than assuming a slow check. Whether the gates should call the interpreter
  directly instead of through `uv` is a change to two scripts and is not made here.
- ✅ **`merge_gate.sh` was retired on 14.08.2026 (Selin) and moved to `scripts/archive/`.** It gated a
  branch before merging; there are no branches, and the last merge is 140 commits back. The capability
  lost with it is the branch-vs-main `coderefs.py` comparison, which nothing replaces because nothing
  produces branches to compare. Every content check it called is still called by `verify_main.sh`.
  **The floors below that refer to `merge_gate.sh` describe a retired script and are kept as the
  record of what was verified, not as live checks.**
- **Branch-specific checks are not carried forward.** `merge_gate.sh` briefly held checks that
  pinned particular strings in `4a` and particular output directories against a rename. They were
  right for one merge and noise on every other. Add one for the merge that needs it, delete it
  after.
- **Three of the seven components had no documented limits until 14.08.2026** — `artifacts.py`,
  `coderefs.py` and `cmdpaths.py`. Each documents its own *scope* in its module docstring, but this
  section, which is where a reviewer looks, named only two checkers and the `merge_gate.sh` floors.
  What the four path-checkers leave between them, stated once: `coderefs.py` reads only
  the `.py`, `.ipynb`, `.png` and `.tex` extensions, so a dangling `.csv` or `.h5ad` in inline code is invisible to it and
  is caught only if it lies under `outputs/` or `figures/`, which is `artifacts.py`'s scope;
  `cmdpaths.py` requires a token to contain `/` **and** start with a tracked top-level directory or
  end with a known source extension, so a bare filename in a shell block is skipped by design. **No
  checker covers a path built by string concatenation in prose**, and none validates a path inside a
  notebook *output* rather than its source.
- ⚠️ **A stale justification lived inside two of them until 14.08.2026.** `artifacts.py` and
  `coderefs.py` both explained skipping csv references partly by "splits/split_ctrp.csv is documented
  as not existing until R2 creates it". **R2 has run**: that file is tracked and 181 rows long. The
  scope restrictions are unchanged — narrowing or widening them is a judgement about noise, not a bug
  fix — but their stated reason named a file that now exists, which is the same defect class the
  checkers themselves exist to catch.
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
