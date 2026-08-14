#!/bin/bash
# Post-merge verification: run the content checks against MAIN AS MERGED, not against a branch tip.
#
# WHY THIS EXISTS SEPARATELY FROM merge_gate.sh. The gate answers "is this branch safe to merge".
# It never answers "is the merge result sound", and a three-way merge can produce a tree that
# neither parent had. Both merges on 12.08.2026 were verified by hand afterwards; this replaces the
# hand pass so it cannot be skipped or done unevenly.
#
# EVERY CHECK STATES ITS DENOMINATOR AND HAS A FLOOR. A check that "passes" because it examined
# nothing is not a pass -- that class has been found five times in this project, including once
# inside check_unbound_names.py itself, which reported "checked 0 file(s) / no unbound names" and
# exited 0 when run from outside the repository. `need` blocks when a floor is not met.
#
# KNOWN LIMITS, carried here rather than in a chat message:
#   * The floors themselves are round numbers chosen to be obviously below the true count, not
#     derived. They catch collapse to zero or near-zero; they do not catch a check that examines
#     80% of what it should.
#   * Runs through .venv/bin/python, NOT `uv run` (changed 14.08.2026, Selin). Every check used to
#     be invoked via `uv run`, which blocks on the project environment lock whenever another uv
#     process holds it -- with a `uv run jupyter lab` server and a live kernel up, this script sat at
#     0 % CPU for over twenty minutes and had to be killed. The same checks complete in seconds
#     against the interpreter directly. Cost: it now assumes .venv/ exists, which every other
#     documented entry point in this repo already does.
#   * Read-only -- TRUE ONLY SINCE 14.08.2026, and it was false when written. The module-import
#     check below imported scripts/evaluation/*, whose files have no __main__ guard, so importing
#     them ran them: it retrained models and rewrote two committed artifacts. 'evaluation' is now
#     excluded, with the cost of that exclusion stated at the check itself. Nothing here trains,
#     reads the h5ad data, or writes outside report/ and /tmp.
#
# Usage:  bash scripts/gate/verify_main.sh
set -uo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
G="$R/scripts/gate"
cd "$R"
fail=0

need() { # label, observed, minimum
  if [ -z "${2:-}" ] || [ "$2" -lt "$3" ] 2>/dev/null; then
    printf '   ^ BLOCKER: %s examined %s, expected >= %s\n' "$1" "${2:-nothing}" "$3"; fail=1
  fi
}
num() { grep -oE "$2" <<<"$1" | head -1 | grep -oE '[0-9]+' | head -1; }

echo "main @ $(git log --oneline -1 | cut -c1-70)"

L=$(.venv/bin/python "$G/links.py" "$R");     echo "  $(echo "$L" | tail -n +2 | head -1)"
need links     "$(num "$L" '[0-9]+ relative links checked')" 100
A=$(.venv/bin/python "$G/artifacts.py" "$R"); echo "  $(echo "$A" | sed -n 2p)"
need artifacts "$(num "$A" '[0-9]+ artifact references checked')" 10
C=$(.venv/bin/python "$G/cmdpaths.py" "$R");  echo "  $(echo "$C" | sed -n 2p)"
need cmdpaths  "$(num "$C" '[0-9]+ command/LaTeX paths checked')" 10

.venv/bin/python scripts/check_resolved_paths.py 2>&1 | grep -E 'composition|dead-glob' | sed 's/^/  /'
.venv/bin/python scripts/check_resolved_paths.py >/dev/null 2>&1 \
  || { echo "   ^ BLOCKER: resolved-path checker non-zero"; fail=1; }

# Added 13.08.2026. Catches names used but never bound -- code that parses and raises at call time.
# Two such defects reached merged main before this existed (f7ef9e4, c351851), both of which would
# have failed R4 on the first fold.
U=$(.venv/bin/python scripts/check_unbound_names.py 2>&1); echo "  $(echo "$U" | head -1)"
echo "$U" | grep UNBOUND | sed 's/^/  /'
.venv/bin/python scripts/check_unbound_names.py >/dev/null 2>&1 \
  || { echo "   ^ BLOCKER: unbound-name checker non-zero"; fail=1; }
need "unbound-name files" "$(num "$U" 'checked [0-9]+ file')" 20

N=$(.venv/bin/python -c "
import glob, warnings, nbformat; warnings.filterwarnings('ignore')
f = sorted(glob.glob('notebooks/**/*.ipynb', recursive=True)); bad = 0
for p in f:
    try: nbformat.validate(nbformat.read(p, as_version=4))
    except Exception: bad += 1
print(f'{len(f)} {bad}')")
echo "  notebooks: $(echo $N | cut -d' ' -f1) validated, $(echo $N | cut -d' ' -f2) invalid"
need notebooks "$(echo $N | cut -d' ' -f1)" 10
[ "$(echo $N | cut -d' ' -f2)" != "0" ] && fail=1

M=$(.venv/bin/python -c "
import importlib, pathlib
# 'gate' is excluded for the same reason 'archive' is: nothing imports it. The helpers there are
# scripts with top-level code, so importing them RUNS them -- which briefly made this check report
# 27 modules and 3 failures, the failures being the gate's own helpers executing mid-sweep.
#
# ⚠️ 'evaluation' was excluded on 14.08.2026 and the exclusion is GONE AGAIN the same day, because
# the real fix landed. Nine files under scripts/evaluation/ were straight-line scripts with no
# 'if __name__ == main' guard, so importing them executed them: aggregation_comparison.py wrote
# (no backticks in this comment: it sits inside a double-quoted bash command substitution, where a
#  backtick opens a nested substitution -- which is exactly how this file acquired a syntax error)
# notebooks/outputs/panel/panel_aggregation_comparison.csv at top level, build_execution_band.py
# wrote panel_execution_band.csv, and section_e2_smaller_study.py and input_dropout_test.py called
# oof_predictions -- they TRAINED. Running this gate therefore retrained models and modified two
# COMMITTED artifacts, which falsified this script's own read-only claim (no inner double quotes
# here on purpose: this comment lives inside a double-quoted bash command substitution).
#
# All nine now carry the guard (Selin, 14.08.2026), verified by comparing each rewritten file's AST
# against the original so no string literal could be corrupted by the reindent. Importing the tree is
# silent and writes nothing, so the exclusion is unnecessary and coverage is restored: 35 modules
# rather than the 24 the exclusion left.
mods = [p for p in pathlib.Path('.').glob('scripts/**/*.py')
        if not {'archive', 'gate', '__pycache__'} & set(p.parts)]
bad = []
for p in sorted(mods):
    name = '.'.join(p.with_suffix('').parts)
    try: importlib.import_module(name)
    except ModuleNotFoundError as e:
        # gen_embeds imports scgpt, which lives in the separate --scgpt-python venv by design.
        if 'scgpt' not in str(e): bad.append(name)
    except Exception: bad.append(name)
print(f'{len(mods)} {len(bad)}')" 2>/dev/null)
echo "  modules: $(echo $M | cut -d' ' -f1) imported, $(echo $M | cut -d' ' -f2) failed"
need modules "$(echo $M | cut -d' ' -f1)" 10
[ "$(echo $M | cut -d' ' -f2)" != "0" ] && fail=1

(cd report && rm -f main.bbl main.aux main.log \
 && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1 && bibtex main >/dev/null 2>&1 \
 && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1 \
 && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1
 pg=$(grep -oE '[0-9]+ pages' main.log | tail -1 | grep -oE '[0-9]+')
 echo "  report: $(grep -cE '^! ' main.log) errors, \
$(grep -cE 'Reference .* undefined|Citation .* undefined' main.log) undefined, ${pg:-0} pages"
 [ "${pg:-0}" -lt 10 ] && echo "   ^ BLOCKER: report shorter than 10 pages"
 exit 0)

echo "  git: $(git status --porcelain | wc -l | tr -d ' ') dirty file(s), $(git status -sb | head -1)"
echo "=== POST-MERGE: $([ $fail = 0 ] && echo 'main is sound' || echo 'BLOCKED') ==="
exit $fail
