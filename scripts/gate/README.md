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
