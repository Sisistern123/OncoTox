"""Three static checks over the notebooks, for defect classes nothing here was catching.

Written 13.08.2026 (Gate 4). Each generalises a defect found by hand during the pipeline review,
because a defect found once by reading will be reintroduced unless something looks for it:

1. **Signature conformance** — every call into ``scripts.*`` checked against the real signature.
   Unexpected keywords, too many positionals, missing required arguments.
2. **Preconditions enforced by ``raise``** — parameters that HAVE a default but whose omission or
   value makes the function raise. A signature check cannot see these: the call is well-formed and
   the failure is three hours into a run. This is the class of ``cv.oof_predictions``'s
   ``counts_h5ad``, which no caller passed, so R4 would have died on its first arm.
3. **Reads without a producer** — every file a notebook reads must be written by something or be a
   known raw source. This is the class of ``drug_catalog`` §5, which read a path that §3 never
   wrote, so the notebook raised there on any machine.

**Every check prints how many candidates it examined and fails when that count collapses to zero**,
the rule the rest of ``scripts/gate/`` follows: "0 failed" and "0 examined" are otherwise
indistinguishable, and a check that silently stops checking is worse than no check.

Known limits, stated rather than discovered later:

* Check 1 resolves only calls whose callee is imported from ``scripts.*`` by name, or reached as
  ``module.function``. A callable passed through a variable is invisible to it.
* Check 2 **enumerates** preconditions; it does not prove the callers satisfy them. The output is a
  list to read, not a pass. Deciding whether a caller satisfies a precondition needs the argument
  values, which are not static.
* Check 3 follows string literals only. A path built through a variable -- ``out = DIR / name`` then
  ``df.to_csv(out)`` -- is not matched, so it can report a false orphan. It did exactly that on its
  first run and the finding was verified by hand.
"""
from __future__ import annotations

import ast
import glob
import importlib
import inspect
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

NOTEBOOKS = [f for f in sorted(glob.glob(str(ROOT / 'notebooks/**/*.ipynb'), recursive=True))
             if '/archive/' not in f and '.ipynb_checkpoints' not in f]
SCRIPTS = [f for f in sorted(glob.glob(str(ROOT / 'scripts/**/*.py'), recursive=True))
           if '/archive/' not in f]


def notebook_code(path: str) -> str:
    nb = json.load(open(path))
    src = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    return '\n'.join(l for l in src.split('\n') if not l.lstrip().startswith(('!', '%')))


def check_signatures() -> tuple[int, list[str]]:
    examined, problems = 0, []
    for path in NOTEBOOKS:
        try:
            tree = ast.parse(notebook_code(path))
        except SyntaxError as e:
            problems.append(f'SYNTAX  {path}: {e}')
            continue
        resolved, modules = {}, {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith('scripts'):
                try:
                    mod = importlib.import_module(node.module)
                except Exception as e:
                    problems.append(f'IMPORT  {path}: {node.module}: {e}')
                    continue
                for a in node.names:
                    obj = getattr(mod, a.name, None)
                    if callable(obj):
                        resolved[a.asname or a.name] = obj
                    try:
                        modules[a.asname or a.name] = importlib.import_module(
                            f'{node.module}.{a.name}')
                    except Exception:
                        pass
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = label = None
            if isinstance(node.func, ast.Name) and node.func.id in resolved:
                fn, label = resolved[node.func.id], node.func.id
            elif (isinstance(node.func, ast.Attribute)
                  and isinstance(node.func.value, ast.Name)
                  and node.func.value.id in modules):
                cand = getattr(modules[node.func.value.id], node.func.attr, None)
                if callable(cand):
                    fn, label = cand, f'{node.func.value.id}.{node.func.attr}'
            if fn is None:
                continue
            examined += 1
            try:
                sig = inspect.signature(fn)
            except (ValueError, TypeError):
                continue
            params = list(sig.parameters.values())
            given = {k.arg for k in node.keywords if k.arg}
            splat = any(k.arg is None for k in node.keywords)
            n_pos = len(node.args)
            if not any(p.kind == p.VAR_KEYWORD for p in params):
                for k in given - {p.name for p in params}:
                    problems.append(f'UNEXPECTED KWARG  {Path(path).name}: {label}(..., {k}=...)')
            positional = [p for p in params
                          if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
            if not any(p.kind == p.VAR_POSITIONAL for p in params) and n_pos > len(positional):
                problems.append(
                    f'TOO MANY POSITIONAL  {Path(path).name}: {label}: {n_pos} > {len(positional)}')
            if not splat:
                for p in params:
                    if p.default is not inspect.Parameter.empty:
                        continue
                    if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
                        continue
                    by_pos = p in positional and positional.index(p) < n_pos
                    if p.name not in given and not by_pos:
                        problems.append(
                            f'MISSING REQUIRED  {Path(path).name}: {label}(...) needs {p.name!r}')
    return examined, problems


def list_preconditions() -> tuple[int, list[str]]:
    """Enumerate raises conditioned on a DEFAULTED parameter. A list to read, not a pass."""
    examined, found = 0, []
    for path in SCRIPTS:
        tree = ast.parse(open(path).read())
        for fn in [n for n in ast.walk(tree)
                   if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
            examined += 1
            params = {a.arg for a in fn.args.args + fn.args.kwonlyargs}
            defaulted = {a.arg for a, _ in zip(fn.args.args[::-1], (fn.args.defaults or [])[::-1])}
            defaulted |= {a.arg for a, d in zip(fn.args.kwonlyargs, fn.args.kw_defaults or [])
                          if d is not None}
            for node in ast.walk(fn):
                if not isinstance(node, ast.If):
                    continue
                if not any(isinstance(n, ast.Raise) for n in ast.walk(node)):
                    continue
                used = {n.id for n in ast.walk(node.test) if isinstance(n, ast.Name)} & params
                if used & defaulted:
                    rel = Path(path).relative_to(ROOT)
                    found.append(f'{rel}:{node.lineno}  {fn.name}()  <- '
                                 f'{", ".join(sorted(used & defaulted))}')
    return examined, found


_LIT = re.compile(r"""['"]([^'"]*\.(?:csv|json|npy|npz|png|h5ad|txt|xlsx|pt))['"]""")
_READ = re.compile(r'(read_csv|read_excel|read_h5ad|read_xml|np\.load|load\()')
_WRITE = re.compile(r'(to_csv|savefig|np\.save|np\.savez_compressed|write_text|to_json|\.write)')
RAW_SOURCES = (
    'CPM_data', 'Metadata.txt', 'UMIcount', 'v20.', 'CTRPv2.csv', 'cellosaurus', 'Repurposing_',
    'full database', 'GDSC2', 'hgnc', 'provenance.json', 'sun2017', 'pubchem_parent',
    'cell_line_names.csv', 'vocab.json', 'best_model.pt', 'run_meta.json', 'summary.json',
    'history.csv', 'per_drug_results.csv', 'args.json', 'oov_genes.csv', 'oov_summary.json',
)


def check_producers() -> tuple[int, list[str]]:
    reads: dict[str, list[str]] = {}
    writes: set[str] = set()
    for path in NOTEBOOKS:
        nb = json.load(open(path))
        for i, c in enumerate(nb['cells']):
            if c['cell_type'] != 'code':
                continue
            for l in ''.join(c['source']).split('\n'):
                st = l.strip()
                if not st or st.startswith('#'):
                    continue
                for m in _LIT.finditer(st):
                    name = m.group(1).split('/')[-1]
                    if _READ.search(st):
                        reads.setdefault(name, []).append(f'{Path(path).name} [{i}]')
                    if _WRITE.search(st):
                        writes.add(name)
    for path in SCRIPTS:
        for m in _LIT.finditer(open(path).read()):
            writes.add(m.group(1).split('/')[-1])
    orphans = [f'{n}  <- {", ".join(sites)}' for n, sites in sorted(reads.items())
               if n not in writes and not any(r in n for r in RAW_SOURCES)]
    return len(reads), orphans


def main() -> int:
    failed = False

    examined, problems = check_signatures()
    print(f'signatures : {examined} call(s) into scripts.* examined, {len(problems)} problem(s)')
    if examined == 0:
        print('  FAIL: resolved zero calls -- the check is not checking anything')
        failed = True
    for p in problems:
        print(f'  {p}')
    failed |= bool(problems)

    n_fns, preconds = list_preconditions()
    print(f'preconditions: {n_fns} function(s) examined, {len(preconds)} enforced by raise on a '
          f'defaulted parameter')
    if n_fns == 0:
        print('  FAIL: examined zero functions')
        failed = True
    print('  (enumerated for reading -- satisfying them needs argument values, which are not static)')

    n_reads, orphans = check_producers()
    print(f'producers  : {n_reads} distinct read filename(s) examined, {len(orphans)} without a '
          f'producer')
    if n_reads == 0:
        print('  FAIL: found no reads -- the check is not checking anything')
        failed = True
    for o in orphans:
        print(f'  ORPHAN READ  {o}')
    # Orphans are reported, not failed: check 3 follows string literals only, so a path built
    # through a variable reads as an orphan when it is not. Judge each one.

    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
