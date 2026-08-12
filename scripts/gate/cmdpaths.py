"""Repo paths written OUTSIDE inline backticks -- in fenced shell blocks and LaTeX \\code{}/\\texttt{}.

The fourth match class. links.py sees markdown links, coderefs.py and artifacts.py see inline-code
spans; none sees `cp ../notebooks/outputs/ablations/x.png` inside a ``` block or
\\code{notebooks/outputs/data/x.csv} in the report. Those are the repo's own quickstart and
figure-build commands, so when they dangle the documented way to reproduce something is broken.

Deliberately narrow, because a shell block is mostly not paths: a token qualifies only if it contains
a "/" and either starts with a known tracked top-level directory or ends with a known source
extension. Tokens with a shell variable, a glob, a URL, or a gitignored root are skipped -- data/,
runs/ and the processed h5ads are absent by design, not by defect.

Usage: cmdpaths.py <tree_root>
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
SKIP_DIRS = {".git", ".claude", "node_modules", ".venv", ".ipynb_checkpoints"}
TRACKED_TOPS = ("scripts/", "notebooks/", "docs/", "report/", "reference/", "outputs/", "figures/")
EXTS = (".py", ".ipynb", ".png", ".tex", ".md", ".csv", ".json", ".bib", ".cff", ".txt")
# Gitignored roots, plus third-party packages quoted by import path (scgpt, torch) and the data
# root's own layout -- none of these are repo files, so absence is not a defect.
IGNORED_ROOTS = ("data/", "runs/", "~/", "/Users/", "splits/", "scgpt/", "torch/",
                 "scRNAseq_SCP542/", "//")
PREFIXES = ["", "notebooks/", "docs/", "notebooks/outputs/", "docs/figures/", "report/"]

FENCE = re.compile(r"```[^\n]*\n(.*?)```", re.S)
TEXSPAN = re.compile(r"\\(?:code|texttt|path)\{([^}]*)\}")
TOKEN = re.compile(r"[A-Za-z0-9_./~-]+")


def markdown_of(p: Path) -> str:
    if p.suffix != ".ipynb":
        return p.read_text()
    nb = json.loads(p.read_text())
    out = []
    for c in nb.get("cells", []):
        if c.get("cell_type") != "markdown":
            continue
        s = c.get("source", "")
        out.append(s if isinstance(s, str) else "".join(s))
    return "\n".join(out)


def candidates(path: Path):
    if path.suffix == ".tex":
        for m in TEXSPAN.finditer(path.read_text()):
            yield m.group(1).replace("\\_", "_").replace("\\", "").split("::")[0]
        return
    for block in FENCE.findall(markdown_of(path)):
        for tok in TOKEN.findall(block):
            yield tok.split("::")[0]


def interesting(ref: str) -> bool:
    if "/" not in ref or "*" in ref or "$" in ref or "://" in ref:
        return False
    if not ref.endswith(EXTS):        # `notebooks/07` is prose, not a path
        return False
    if ref.startswith(IGNORED_ROOTS):
        return False
    core = ref.lstrip("./")
    return core.startswith(TRACKED_TOPS) or ref.endswith(EXTS)


def resolves(ref: str, near: Path) -> bool:
    clean = ref.lstrip("./") if ref.startswith("./") else ref
    if (near.parent / ref).exists() or (near.parent / clean).exists():
        return True
    return any((ROOT / pre / clean).exists() for pre in PREFIXES)


missing: dict[str, list[str]] = {}
n = 0
for path in sorted(ROOT.rglob("*")):
    if path.suffix not in {".md", ".ipynb", ".tex"} or not path.is_file():
        continue
    if any(x in SKIP_DIRS for x in path.relative_to(ROOT).parts):
        continue
    try:
        seen = {c for c in candidates(path) if interesting(c)}
    except Exception:
        continue
    for ref in sorted(seen):
        n += 1
        if not resolves(ref, path):
            missing.setdefault(ref, []).append(str(path.relative_to(ROOT)))

print(f"{ROOT}\n  {n} command/LaTeX paths checked, {len(missing)} unresolved")
for ref, where in sorted(missing.items()):
    print(f"   {ref}")
    for w in sorted(where):
        print(f"        {w}")
