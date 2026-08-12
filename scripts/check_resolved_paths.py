"""Check paths that are BUILT FROM VARIABLES, which no string search can see.

Run it:  python <this file>          -> prints findings, exits 1 if any
Exit 0 = clean.

WHY THIS EXISTS. Three defects in one week lived in path *composition* rather than in any path
string, so grep, link checkers and code-reference checkers all passed them:

    OUT_MATRIX = NB_DIR / 'outputs' / 'matrix'          # definition
    OUT_MATRIX / 'legacy' / 'training_545_mean_pv'      # call site, one level too deep

Neither line contains a wrong string. The defect is the composition. One of the three would have
silently recomputed a cross-validation while printing that it had loaded committed folds.

TWO CHECKS, AND NEITHER SUBSUMES THE OTHER. Do not "simplify" them into one:

  (a) COMPOSITION   VAR / 'a' / 'b'      -> the parent directory must exist
  (b) DEAD GLOB     VAR.glob('*.csv')    -> must match at least one file

The proof that both are needed is the defect that prompted this file. With
``OUT_MATRIX = notebooks/outputs``, check (a) PASSES -- that directory exists -- and the notebook
still reports nothing, because the files it globs for live two levels deeper. Only (b) catches it.
(b) is also the more dangerous class: its failure mode is silence, not an exception. A reporting
cell that finds zero files looks exactly like a run with nothing to report.

TWO BLIND SPOTS. Read them before treating a clean run as proof.

  1. Only string-literal chains from known roots resolve. A path built through an f-string, a loop
     variable, a function return or a dict lookup is invisible here. This is a FLOOR, not a
     ceiling: it would have caught all three of the defects it was written for, and it would miss a
     fourth built differently.

  2. The dead-glob check is only meaningful where the target artifacts are actually present, so a
     glob into a gitignored tree is a false positive, not a defect. ``runs/``, ``data/`` and
     ``splits/`` are skipped for that reason -- their contents legitimately do not exist in a fresh
     clone. Widening the root table without widening that skip list will produce noise, and a noisy
     checker trains people to skim its output, which is worse than not having one.

  3. ``#`` comments are stripped before parsing, so a comment quoting a path expression is not
     read as live code. Without this the checker fires on its own explanatory comments and on every
     correction anyone documents -- penalising the habit of recording what changed.

EXPLICITLY OUT OF SCOPE: STALE COMMENTS. A comment that misdescribes the code around it -- "the
call sites append that suffix" after they stopped, or "that question is open" after it was decided --
is a real defect and this checker will never catch one. Do not widen it to try. The reason is that
there is no bad path to find: every path named in such a comment usually resolves perfectly, and the
falsehood is a claim about the code's *structure* or a status, not about a location. Parsing comments
would instead produce false positives on deliberate historical references, which this repository uses
by convention -- `cf3ad3f:notebooks/2_training.ipynb cell 1` is *meant* not to resolve against the
working tree, and the archived-notebook citations are the same. A checker that cries wolf on
legitimate history gets ignored, and the real findings go with it.

The durable fix for that class is a convention, not a parser: **a change that supersedes a comment
retires that comment in the same commit.** Stated here so the next person does not automate it and
conclude the checker is broken when it cannot.

WHY IT PRINTS ITS DENOMINATOR. Each check reports how many candidates it examined, not only how
many failed, because "0 failed" and "0 examined" are otherwise indistinguishable -- and a checker
that silently parses nothing reports a clean pass. That is the same silent-success failure this
file exists to catch, turned on the file itself. If a notebook schema changes, or the root table
stops matching how paths are written, the candidate count collapses and THAT is the signal. A pass
that cannot state its denominator is not a pass.
"""

from __future__ import annotations

import glob as globmod
import json
import os
import re
import sys
from pathlib import Path

# Resolve the repository root from this file, not from the working directory: the gate runs it from
# wherever it happens to be, and every glob below is repo-relative.
REPO_ROOT = Path(__file__).resolve().parent.parent

#: Path variables whose value is known statically. Extend deliberately -- see blind spot 2.
KNOWN_ROOTS = {
    "ROOT": ".",
    "NB_DIR": "notebooks",
    "HERE": "docs",
    "PANEL_OUT": "notebooks/outputs/panel",
    "LEGACY_PANEL": "notebooks/outputs/legacy/panel_void_8drug",
}

#: Gitignored or externally-supplied trees: absent by design, so a miss there is not a defect.
SKIP_PREFIXES = ("runs", "data", "splits")

_ASSIGN = re.compile(
    r"^\s*([A-Z][A-Z0-9_]*)\s*=\s*([A-Za-z_]\w*(?:\s*/\s*'[^']+')+)\s*(?:;|$)", re.M)
_USE = re.compile(r"\b([A-Z][A-Z0-9_]*)((?:\s*/\s*'[^']+')+)")
_GLOB = re.compile(r"\b([A-Z][A-Z0-9_]*)((?:\s*/\s*'[^']+')*)\s*\.glob\(\s*'([^']+)'\s*\)")


def _literals(expr: str) -> list[str]:
    return re.findall(r"'([^']+)'", expr)


def _skip(path: str) -> bool:
    head = os.path.normpath(path).split(os.sep)[0]
    return head in SKIP_PREFIXES


def _strip_comments(text: str) -> str:
    """Remove ``#`` comments so prose *describing* a path is not read as one.

    Necessary rather than cosmetic: comments in this repository routinely quote the path expression
    they are explaining ("was OUT_MATRIX / 'legacy' / ...", "fixed in ..."), and a checker that
    reads those as live code reports a defect for every correction anyone documents -- punishing
    exactly the habit the project wants. Quote-aware, so a ``#`` inside a string literal survives.
    """
    out = []
    for line in text.splitlines():
        quote = None
        for i, ch in enumerate(line):
            if quote:
                if ch == quote:
                    quote = None
            elif ch in "'\"":
                quote = ch
            elif ch == "#":
                line = line[:i]
                break
        out.append(line)
    return "\n".join(out)


def _sources() -> dict[str, str]:
    """Every live source file, as text. Archived trees are records, not live code."""
    out: dict[str, str] = {}
    paths = sorted(globmod.glob("notebooks/**/*.ipynb", recursive=True))
    paths += sorted(globmod.glob("scripts/**/*.py", recursive=True))
    paths += sorted(globmod.glob("docs/**/*.py", recursive=True))
    for f in paths:
        if "/archive/" in f or ".ipynb_checkpoints" in f:
            continue
        # This file quotes broken path expressions in its own docstring, as the worked example that
        # explains why both checks are needed. Analysing itself would report those as defects --
        # documentation read as code. Excluded by name rather than by stripping docstrings, because
        # docstring-stripping is reliable for .py and not for notebook cells, and a rule that works
        # in one place and not the other is worse than an explicit exception.
        if os.path.basename(f) == os.path.basename(__file__):
            continue
        if f.endswith(".ipynb"):
            nb = json.load(open(f))
            out[f] = _strip_comments("\n".join(
                c["source"] if isinstance(c["source"], str) else "".join(c["source"])
                for c in nb["cells"] if c["cell_type"] == "code"))
        else:
            out[f] = _strip_comments(open(f, errors="ignore").read())
    return out


def check() -> tuple[list[str], dict[str, int]]:
    findings: list[str] = []
    seen = {"composition": 0, "dead-glob": 0}
    for f, text in _sources().items():
        table: dict[str, str] = {}
        for var, expr in _ASSIGN.findall(text):
            base = expr.split("/")[0].strip()
            root = table.get(base, KNOWN_ROOTS.get(base))
            if root is not None:
                table[var] = os.path.normpath(os.path.join(root, *_literals(expr)))

        for var, chain in _USE.findall(text):
            if var not in table:
                continue
            seen["composition"] += 1
            target = os.path.normpath(os.path.join(table[var], *_literals(chain)))
            parent = os.path.dirname(target)
            if parent and not _skip(parent) and not os.path.isdir(parent):
                findings.append(
                    f"COMPOSITION  {f}\n"
                    f"    {var} = {table[var]}\n"
                    f"    {var}{chain.strip()} -> {target}\n"
                    f"    parent directory does not exist: {parent}")

        for var, chain, pattern in _GLOB.findall(text):
            if var not in table:
                continue
            d = os.path.normpath(os.path.join(table[var], *_literals(chain))) if chain else table[var]
            if _skip(d):
                continue
            seen["dead-glob"] += 1
            if not globmod.glob(os.path.join(d, pattern)):
                findings.append(
                    f"DEAD GLOB    {f}\n"
                    f"    {var} = {table[var]}\n"
                    f"    {var}{chain}.glob('{pattern}') matches nothing in {d}\n"
                    f"    (a reporting cell that finds nothing looks like a run with nothing to report)")
    return findings, seen


def main() -> int:
    os.chdir(REPO_ROOT)
    findings, seen = check()
    for x in findings:
        print(x)
    n_comp = sum(1 for f in findings if f.startswith("COMPOSITION"))
    n_glob = sum(1 for f in findings if f.startswith("DEAD GLOB"))
    print(f"\ncomposition: {seen['composition']} expressions checked, {n_comp} failed")
    print(f"dead-glob  : {seen['dead-glob']} globs checked, {n_glob} empty")
    if not any(seen.values()):
        print("WARNING: nothing was examined -- the parser matched no candidates at all, which is "
              "itself a defect, not a pass.")
        return 1
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
