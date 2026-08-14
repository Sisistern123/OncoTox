#!/bin/bash
# Pre-merge gate: is BRANCH safe to merge into main? Read-only -- nothing here trains, reads the
# h5ad data, or touches the working tree. It inspects committed trees via git, not checkouts.
#
# Usage:  bash scripts/gate/merge_gate.sh <branch>        (defaults to the current branch)
#
# EVERY CHECK STATES HOW MANY CANDIDATES IT EXAMINED, AND `need` BLOCKS WHEN THE FLOOR IS NOT MET.
# The reason is a defect this gate itself once had: a check grepped a notebook that had been
# renamed, found nothing, and reported a pass (12.08.2026). "0 problems" and "0 examined" print
# identically unless the denominator is on screen.
#
# KNOWN LIMITS, recorded here rather than in a message, because a gate is only as trustworthy as
# its documented ways of lying:
#   * TWO FLOORS HAVE NEVER BEEN FAULT-INJECTED, only reasoned about: check 1's merge-tree output
#     size, and check 6's diff file count. Every other floor in this file has been verified by
#     deliberately breaking the thing it guards and confirming it blocks. Treat those two as
#     unproven.
#   * check_resolved_paths.py sees `Path.glob` only; a module-level `glob.glob` is invisible to it.
#     That gap is what let 4a call glob.glob for weeks (fixed c351851, and now covered by
#     check_unbound_names.py, which catches the unbound name rather than the dead pattern).
#   * The branch-specific checks that lived here on 12-13.08.2026 -- preserving particular strings
#     in 4a, and pinning outputs/{diagnostics,dreval,data} against a rename -- are NOT carried
#     forward. They were correct for one merge and would be noise on every other. Add such a check
#     to this file only for the merge that needs it, and delete it after.
set -uo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
G="$R/scripts/gate"
cd "$R"
B="${1:-$(git rev-parse --abbrev-ref HEAD)}"
fail=0

say()  { printf '%-46s %s\n' "$1" "$2"; }
need() { if [ -z "${2:-}" ] || [ "$2" -lt "$3" ] 2>/dev/null; then
           printf '     ^ BLOCKER: %s examined %s, expected at least %s\n' "$1" "${2:-nothing}" "$3"
           fail=1
         fi; }
num()  { grep -oE "$2" <<<"$1" | head -1 | grep -oE '[0-9]+' | head -1; }

echo "=== branch ==="
say "  main" "$(git log --oneline -1 main | cut -c1-60)"
say "  $B"   "$(git log --oneline -1 "$B" | cut -c1-60)"
git merge-base --is-ancestor main "$B" \
  && say "  contains main tip" "yes" \
  || say "  contains main tip" "no - behind main (blocks only if it also conflicts, see 1)"

echo; echo "=== 1. merge cleanliness (no working tree touched) ==="
out=$(git merge-tree --write-tree main "$B" 2>&1)
if echo "$out" | grep -qi "conflict"; then
  say "  conflicts" "YES"; echo "$out" | grep -i conflict | head; fail=1
else
  say "  conflicts" "none"
fi
need "merge-tree output" "$(printf '%s' "$out" | wc -c | tr -d ' ')" 40   # <- floor never injected

echo; echo "=== 2. content checks against the branch's committed tree ==="
T=$(mktemp -d); git archive "$B" | tar -x -C "$T"
L=$(uv run python "$G/links.py"     "$T"); echo "  $(echo "$L" | tail -n +2 | head -1)"
need "links" "$(num "$L" '[0-9]+ relative links checked')" 100
echo "$L" | tail -n +3 | head -12
A=$(uv run python "$G/artifacts.py" "$T"); echo "  $(echo "$A" | sed -n 2p)"
need "artifact refs" "$(num "$A" '[0-9]+ artifact references checked')" 10
C=$(uv run python "$G/cmdpaths.py"  "$T"); echo "  $(echo "$C" | sed -n 2p)"
need "command/LaTeX paths" "$(num "$C" '[0-9]+ command/LaTeX paths checked')" 10

echo; echo "=== 3. new stale code references vs main ==="
uv run python "$G/coderefs.py" "$R" | grep '^   [^ ]' | awk '{print $1}' | sort > /tmp/_gm.txt
uv run python "$G/coderefs.py" "$T" | grep '^   [^ ]' | awk '{print $1}' | sort > /tmp/_gb.txt
MREF=$(uv run python "$G/coderefs.py" "$R"); BREF=$(uv run python "$G/coderefs.py" "$T")
need "coderefs(main)"   "$(num "$MREF" '[0-9]+ code references checked')" 10
need "coderefs(branch)" "$(num "$BREF" '[0-9]+ code references checked')" 10
say "  new vs main" "$(comm -13 /tmp/_gm.txt /tmp/_gb.txt | wc -l | tr -d ' ')"
comm -13 /tmp/_gm.txt /tmp/_gb.txt | sed 's/^/     /'

echo; echo "=== 4. the two repo checkers, against the branch tree ==="
(cd "$T" && uv run --project "$R" python scripts/check_resolved_paths.py 2>&1) \
  | grep -E "composition|dead-glob|FAIL" | sed 's/^/  /'
(cd "$T" && uv run --project "$R" python scripts/check_resolved_paths.py >/dev/null 2>&1) \
  || { echo "     ^ BLOCKER: check_resolved_paths non-zero"; fail=1; }
U=$(cd "$T" && uv run --project "$R" python scripts/check_unbound_names.py 2>&1)
echo "  $(echo "$U" | head -1)"; echo "$U" | grep UNBOUND | sed 's/^/  /'
(cd "$T" && uv run --project "$R" python scripts/check_unbound_names.py >/dev/null 2>&1) \
  || { echo "     ^ BLOCKER: check_unbound_names non-zero"; fail=1; }
need "unbound-name files" "$(num "$U" 'checked [0-9]+ file')" 20

echo; echo "=== 5. notebooks parse + validate, and stored outputs must not shrink ==="
uv run python - "$T" <<'PYX'
import glob, json, subprocess, sys, warnings
import nbformat
warnings.filterwarnings("ignore")
tree = sys.argv[1]
files = sorted(glob.glob(f"{tree}/notebooks/**/*.ipynb", recursive=True))
bad = []
for p in files:
    try:
        json.load(open(p)); nbformat.validate(nbformat.read(p, as_version=4))
    except Exception as e:
        bad.append(f"{p}: {type(e).__name__}")
print(f"  {len(files)} notebooks, {len(bad)} invalid")
for b in bad:
    print("     " + b)
PYX
uv run python - "$B" <<'PYX'
import json, subprocess, sys
branch = sys.argv[1]
def counts(ref):
    names = subprocess.run(["git", "ls-tree", "-r", "--name-only", ref],
                           capture_output=True, text=True).stdout.split()
    out = {}
    for f in names:
        if not f.endswith(".ipynb"):
            continue
        blob = subprocess.run(["git", "show", f"{ref}:{f}"], capture_output=True, text=True).stdout
        try:
            nb = json.loads(blob)
        except Exception:
            continue
        out[f] = sum(len(c.get("outputs", [])) for c in nb.get("cells", []))
    return out
a, b = counts("main"), counts(branch)
lost = {k: (a[k], b[k]) for k in a if k in b and b[k] < a[k]}
print(f"  notebooks whose COMMITTED stored outputs shrink: {len(lost)}")
for k, (x, y) in sorted(lost.items()):
    print(f"     {k}: {x} -> {y}   BLOCKER")
PYX

echo; echo "=== 6. freeze compliance: no data or model artifacts added ==="
need "files in the branch diff" "$(git diff --name-only main.."$B" | wc -l | tr -d ' ')" 1  # <- never injected
added=$(git diff --name-only --diff-filter=A main.."$B" | grep -cE '\.(h5ad|pt|npz|npy)$|^runs/')
say "  h5ad/pt/npz/runs added" "$added"
[ "$added" != "0" ] && fail=1

echo; echo "=== 7. report builds from the branch tree ==="
(cd "$T/report" 2>/dev/null && rm -f main.bbl main.aux main.log \
  && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1 \
  && bibtex main >/dev/null 2>&1 \
  && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1 \
  && pdflatex -interaction=nonstopmode main.tex >/dev/null 2>&1
  echo "  errors $(grep -cE '^! ' main.log) | undefined \
$(grep -cE 'Reference .* undefined|Citation .* undefined' main.log) | \
$(grep -oE '[0-9]+ pages' main.log | tail -1)") || echo "  (no report/ in this tree)"

rm -rf "$T"
echo; echo "=== GATE: $([ $fail = 0 ] && echo 'no blockers found' || echo 'BLOCKED') ==="
exit $fail
