"""Resolve every relative markdown link in a tree against the filesystem.

Covers .md files AND the markdown cells of .ipynb notebooks, which the first version missed
entirely -- notebook prose carries deep-links into docs/steps just as the docs do.

Two failure modes this is written to avoid, both of which a pattern-only checker passes:
  * a notebook cell's ``source`` is sometimes a str and sometimes a list of lines; iterating the
    str yields characters, so those cells get silently skipped;
  * link targets are resolved against the filesystem and heading anchors against the real
    headings, so a rewrite that merely looks plausible ("../.." -> "../../..") is still caught.

Anchor rule is GitHub's: lowercase, drop inline code/emphasis/link syntax, strip every character
that is not word/space/hyphen (underscores are word characters and survive), then replace EACH
space with one hyphen -- runs of spaces are not collapsed.

Usage: links.py [tree_root]
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "/Users/selin/PycharmProjects/OncoTox").resolve()
SKIP_DIRS = {".git", ".claude", "node_modules", ".venv", ".ipynb_checkpoints"}

_anchor_cache: dict[Path, set[str]] = {}


def _cells(nb_path):
    try:
        nb = json.loads(nb_path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "markdown":
            continue
        src = cell.get("source", "")
        # str or list-of-lines -- both occur in this repo; joining a str would iterate characters
        yield src if isinstance(src, str) else "".join(src)


def text_of(path: Path) -> str:
    if path.suffix == ".ipynb":
        return "\n".join(_cells(path))
    try:
        return path.read_text()
    except UnicodeDecodeError:
        return ""


def anchors(path: Path) -> set[str]:
    if path in _anchor_cache:
        return _anchor_cache[path]
    out = set()
    for line in text_of(path).splitlines():
        m = re.match(r"^#{1,6}\s+(.*)", line)
        if not m:
            continue
        t = m.group(1).strip().lower()
        t = re.sub(r"`([^`]*)`", r"\1", t)
        t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", t)
        t = re.sub(r"\*", "", t)
        t = re.sub(r"[^\w\s-]", "", t)
        out.add(t.replace(" ", "-"))
    _anchor_cache[path] = out
    return out


def strip_code(text: str) -> str:
    """Drop fenced blocks and inline code, so `N[Pt](N)(Cl)Cl` is not read as a link."""
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    return re.sub(r"`[^`\n]*`", "", text)


bad, checked = [], 0
for path in sorted(ROOT.rglob("*")):
    if path.suffix not in {".md", ".ipynb"} or not path.is_file():
        continue
    if any(p in SKIP_DIRS for p in path.relative_to(ROOT).parts):
        continue
    for m in re.finditer(r"\[[^\]]*\]\(([^)\s]+)\)", strip_code(text_of(path))):
        tgt = m.group(1)
        if tgt.startswith(("http", "mailto:", "#!")):
            continue
        checked += 1
        rel, _, frag = tgt.partition("#")
        dest = (path.parent / rel).resolve() if rel else path
        if not dest.exists():
            bad.append(f"{path.relative_to(ROOT)}\n      -> {tgt}   MISSING FILE")
        elif frag and dest.suffix in {".md", ".ipynb"} and frag not in anchors(dest):
            bad.append(f"{path.relative_to(ROOT)}\n      -> {tgt}   MISSING ANCHOR")

print(f"{ROOT}\n  {checked} relative links checked, {len(bad)} broken")
for b in bad:
    print("   " + b)
