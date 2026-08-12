"""Backticked artifact references that no longer resolve.

A third match class, distinct from markdown links and from code references: an output file cited in
prose. `outputs/panel/panel.csv` in a table cell is neither a link nor a module path, so neither of
the other two checkers sees it, and it is exactly what a directory move leaves dangling.

Scoped to outputs/ and figures/ deliberately. A blanket check over every backticked .csv is noise,
because runs/, data/ and the processed h5ads are gitignored and splits/split_ctrp.csv is documented
as not existing until R2 creates it -- all legitimately absent. Under outputs/ and figures/ the files
are tracked, so absence is a defect.

Both prefixes are written in docs in their short form, relative to the directory a reader is assumed
to be in: `outputs/...` means notebooks/outputs/..., `figures/...` means docs/figures/... . Each
candidate is resolved against every plausible root before being called missing.

Usage: artifacts.py <tree_root>
"""

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
SKIP = {".git", ".claude", "node_modules", ".venv", ".ipynb_checkpoints"}
# Gitignored files are not repo content -- docs/progress_report_*.md is Selin's personal working
# record, untracked by design, and a stale artifact reference inside it is hers, not the repo's.
_ignored = set(subprocess.run(["git", "-C", str(ROOT), "ls-files", "--others", "--ignored",
                               "--exclude-standard", "--directory"],
                              capture_output=True, text=True).stdout.split())
PAT = re.compile(r"`((?:outputs|figures)/[^`\n]*?\.(?:csv|png|npz|json|h5ad|txt))`")
PREFIXES = ["", "notebooks/", "docs/", "notebooks/outputs/", "docs/figures/"]


def text_of(p: Path) -> str:
    if p.suffix == ".ipynb":
        nb = json.loads(p.read_text())
        out = []
        for c in nb.get("cells", []):
            if c.get("cell_type") != "markdown":
                continue
            s = c.get("source", "")
            out.append(s if isinstance(s, str) else "".join(s))
        return "\n".join(out)
    return p.read_text()


def resolves(ref: str, near: Path) -> bool:
    if (near.parent / ref).exists():
        return True
    return any((ROOT / pre / ref).exists() for pre in PREFIXES)


missing: dict[str, list[str]] = {}
n = 0
for path in sorted(ROOT.rglob("*")):
    if path.suffix not in {".md", ".ipynb", ".tex"} or not path.is_file():
        continue
    if any(x in SKIP for x in path.relative_to(ROOT).parts):
        continue
    if str(path.relative_to(ROOT)) in _ignored:
        continue
    try:
        body = text_of(path)
    except Exception:
        continue
    for ref in sorted(set(PAT.findall(body))):
        if '*' in ref:      # a glob is a description, not a path
            continue
        n += 1
        if not resolves(ref, path):
            missing.setdefault(ref, []).append(str(path.relative_to(ROOT)))

print(f"{ROOT}\n  {n} artifact references checked, {len(missing)} unresolved")
for ref, where in sorted(missing.items()):
    print(f"   {ref}")
    for w in sorted(where):
        print(f"        {w}")
