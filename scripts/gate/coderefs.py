"""Find inline-code references to tracked repo paths that no longer exist.

Restricted to .py/.ipynb/.png/.tex on purpose. Extending it to .csv/.h5ad only produces
noise: runs/, data/ and the outputs artifacts are gitignored, and splits/split_ctrp.csv is
documented as not existing until R2 creates it -- all legitimately absent from the repo.

The link checker deliberately strips code spans (so a SMILES string in backticks is not read as a
link), which means `notebooks/2_training.ipynb` in prose is invisible to it. After a rename that is
exactly where the breakage hides: the reference still reads correctly and points at nothing.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
SKIP = {".git", ".claude", "node_modules", ".venv", ".ipynb_checkpoints"}
PAT = re.compile(r"`([^`\n]*?\.(?:ipynb|py|png|tex))`")

def text_of(p):
    if p.suffix == ".ipynb":
        nb = json.loads(p.read_text())
        return "\n".join(
            (c.get("source", "") if isinstance(c.get("source", ""), str) else "".join(c["source"]))
            for c in nb.get("cells", []) if c.get("cell_type") == "markdown")
    return p.read_text()

missing = {}
n_checked = 0
for p in sorted(ROOT.rglob("*")):
    if p.suffix not in {".md", ".ipynb", ".tex"} or not p.is_file():
        continue
    if any(x in SKIP for x in p.relative_to(ROOT).parts):
        continue
    try:
        body = text_of(p)
    except Exception:
        continue
    for ref in set(PAT.findall(body)):
        ref = ref.replace("\\_", "_").strip()
        ref = ref.split()[-1] if " " in ref else ref   # `uv run docs/make_figures.py`
        if "/" not in ref:          # bare filenames are usually prose, too noisy
            continue
        # Same multi-root resolution as artifacts.py: docs cite these in short form relative to
        # the directory a reader is assumed to be in, so `archive/target/x.png` in a docs/ file
        # means notebooks/outputs/archive/target/x.png. Resolving against ROOT and the citing
        # file's parent alone reports the whole convention as broken.
        n_checked += 1
        roots = ('', 'notebooks/', 'docs/', 'notebooks/outputs/', 'docs/figures/', 'report/')
        if (p.parent / ref).resolve().exists() or any((ROOT / r / ref).exists() for r in roots):
            continue
        missing.setdefault(ref, []).append(str(p.relative_to(ROOT)))

print(f"{ROOT}\n  {n_checked} code references checked, {len(missing)} do not exist")
for ref, where in sorted(missing.items(), key=lambda kv: -len(kv[1])):
    print(f"   {ref}   <- {len(where)} file(s)")
    for w in sorted(where):
        print(f"        {w}")
