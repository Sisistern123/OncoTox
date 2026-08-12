"""Find names a notebook or module uses but never binds — the defect that parses and then raises.

Like ``check_resolved_paths.py`` beside it, this is **not pipeline code**: nothing imports it, and it
runs as a pre-merge check. It exists for the same reason — it earned its place by catching defects
that had already bitten, twice, on merged ``main``:

1. ``train_multitask.cv_evaluate`` used ``CV_FOLD_PCA_KEY`` and ``fold_pca_projections_for`` and
   imported neither (fixed ``f7ef9e4``). The module imported cleanly, so the ``NameError`` waited
   until call time — section B's CV harness, first fold, both representation arms.
2. ``4a_percell_training`` called ``glob.glob`` and ``json.load`` and imported neither module
   (fixed ``c351851``). Same shape: every cell parsed, section B died on its first call.

Both were missed by a hardcoded-value sweep, a panel sweep and a link check, because none of those
ask the question this asks. Syntax checks do not catch it either — the code is perfectly well-formed.
Only executing the notebook or walking its bindings will find it, and executing it costs a training
run.

**How it decides.** For each file it walks the AST accumulating *bindings* — imports, assignments,
``def``/``class``, function arguments, ``global``/``nonlocal``, ``except ... as`` — and then flags
every ``Load`` of a name that is in neither the accumulated bindings nor ``builtins``. A notebook is
treated as one namespace across its code cells in order, which is how it actually executes.

**Bindings for a cell are collected BEFORE that cell's loads are checked**, deliberately. A function
defined in a cell may legitimately reference a global defined lower down in the same cell, and more
importantly this is what stops the check flooding.

⚠️ **Its failure mode is over-reporting, and that is only safe while the over-report is absurd.**
The first version of this check collected bindings *after* checking each cell's loads, so every name
defined and used within one cell was flagged. It reported **870 findings across 14 notebooks**. That
was caught only because 870 is an implausible number — at 4 findings it would have been believed, and
the real defect would have been buried in three false ones. If this check ever returns a large
number, suspect the check before the code. A checker that can flood is only safe while somebody
notices the flood.

**Known false positives, excluded:** ``__file__`` and the IPython injections (``get_ipython``,
``display``, ``In``, ``Out``, ``_``) are module- or kernel-provided rather than bound in the source.

Usage::

    python scripts/check_unbound_names.py                 # scripts/ + notebooks/, skipping archive/
    python scripts/check_unbound_names.py path/to/x.ipynb # one file

Exits non-zero if anything is found, so it can gate a merge.
"""

from __future__ import annotations

import ast
import builtins
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Names present at runtime but never bound in the source. Kept short on purpose: every entry here is
#: a hole in the check, so an addition needs a reason beyond "it was noisy".
PROVIDED = {"__file__", "get_ipython", "display", "In", "Out", "_"}

SKIP_DIRS = {"archive", ".ipynb_checkpoints", ".venv", "outputs"}


def _bindings(tree: ast.AST) -> set[str]:
    """Every name this tree binds, by any mechanism Python has for binding one."""
    bound: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                bound.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            bound.add(node.id)
        elif isinstance(node, ast.arg):
            bound.add(node.arg)
        elif isinstance(node, (ast.Global, ast.Nonlocal)):
            bound.update(node.names)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)
    return bound


def _unbound(sources: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Walk ``(label, source)`` chunks sharing one namespace; return ``(label, name)`` for each miss.

    Chunks are a module's single source, or a notebook's code cells in execution order.
    """
    known = set(dir(builtins)) | PROVIDED
    found, seen = [], set()
    for label, src in sources:
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue                      # a syntax error is a different check's finding
        known |= _bindings(tree)          # BEFORE the loads below -- see the flooding note above
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
                if node.id not in known and (label, node.id) not in seen:
                    seen.add((label, node.id))
                    found.append((label, node.id))
    return found


def check_file(path: Path) -> list[str]:
    """Return human-readable findings for one ``.py`` or ``.ipynb``."""
    if path.suffix == ".ipynb":
        cells = json.loads(path.read_text())["cells"]
        sources = [
            (f"code-cell {i}", c["source"] if isinstance(c["source"], str) else "".join(c["source"]))
            for i, c in enumerate(cells) if c["cell_type"] == "code"
        ]
    else:
        sources = [("module", path.read_text())]
    rel = path.relative_to(ROOT) if path.is_absolute() else path
    return [f"{rel} {label}: {name}" for label, name in _unbound(sources)]


def iter_targets(roots: list[Path]):
    for root in roots:
        if root.is_file():
            yield root
            continue
        for pattern in ("*.py", "*.ipynb"):
            for p in sorted(root.rglob(pattern)):
                if not SKIP_DIRS & set(p.parts):
                    yield p


def main(argv: list[str]) -> int:
    roots = [Path(a) for a in argv[1:]] or [ROOT / "scripts", ROOT / "notebooks"]
    targets = list(iter_targets(roots))
    findings = [f for p in targets for f in check_file(p)]

    # The denominator is printed always, not only on failure: "0 findings" is not interpretable
    # without knowing whether it examined 2 files or 200. That distinction has cost this project
    # real time -- a module sweep once covered 2 of 22 and reported success.
    print(f"checked {len(targets)} file(s)")
    for f in findings:
        print(f"  UNBOUND  {f}")
    if findings:
        print(f"\n{len(findings)} unbound name(s). If this number is large, suspect this checker "
              f"before the code -- see the flooding note in the module docstring.")
        return 1
    print("no unbound names")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
