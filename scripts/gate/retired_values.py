"""Values and phrases that have been RETIRED — flagged wherever they reappear as live text.

**Why this exists.** The single most repeated defect in this repository is not a wrong number. It is a
*corrected* number that survives somewhere the correction did not reach. Every instance was found by
hand, weeks or hours late, and each time the fix was applied to the places someone happened to look:

* ``Q2_CONTROL_THRESHOLD`` was recorded as outstanding in **four** places after it was closed, and a
  **fifth** was found in the report a day later (``b8b3131``, then ``ce80f30``).
* The synthetic mean-effects figure **0.98** was cleared on 12.08.2026 from three named locations. It
  reached **one**; the other two kept asserting it, one under the word *Demonstrated* (``f75e9bb``).
* The gene counts ``4,576`` / ``20,570`` stayed live in fourteen places in ``docs/steps/02`` after the
  sweep replaced them, alongside a prediction of ``4,704`` that the measurement superseded (``9fafcff``).
* ``545`` stood in for the model's head count long after the count became ``534`` (``ef45d75``).
* ``180`` trainable lines outlived the ``H292`` alias that made it ``181`` (``914bb50``).

A grep finds these in seconds **once you know to look**. This file is the "know to look" part, made
mechanical: a registry of what has been retired, where it is still legitimately allowed to appear, and
a check that fails when it appears anywhere else.

**What "allowed" means.** A retired value is not deleted from the repository — that would destroy the
record. It stays in the files whose job is to hold superseded content (``corrections-and-dead-ends``),
in dated correction markers, and in the entry that retired it. ``allowed`` lists those files. Anywhere
else, the value is being asserted as live and the check fails.

**Adding an entry is part of retiring a value.** When a number is corrected, add it here with the
commit that retired it. That is the step every instance above skipped.

⚠️ **Known limits, because a check is only as good as its documented ways of lying.**

* It matches **text**, so it cannot tell an assertion from a mention. A file that discusses a retired
  value legitimately must be listed in ``allowed`` — which means a careless addition to ``allowed``
  silences the check for that file entirely. Keep ``allowed`` as narrow as possible.
* It only knows what is registered. It cannot find a retired value nobody wrote down, which is the
  same failure it exists to prevent, one level up.
* Notebooks are read as **source cells only**. A retired value printed into a stored output is a
  record of a past run and is not flagged — which also means a stale number that exists *only*
  in an output is invisible here.
* ``docs/final_presentation*.md`` and ``docs/progress_report_*.md`` are **not scanned**: they are
  untracked by design and one is a dated historical record.
* Numbers too generic to grep (``0.31``, ``0.7``) are deliberately absent — matching them would drown
  the check in false positives, which is how a check gets ignored.

Usage:  retired_values.py <tree_root>
"""
import re
import sys
from pathlib import Path

#: (pattern, what it was, what replaced it, commit that retired it, files allowed to still contain it)
#: Paths in ``allowed`` are prefixes, matched against the repo-relative path.
RETIRED = [
    (r"Q2_CONTROL_THRESHOLD", "a threshold said to be outstanding", "no such constant ever existed; "
     "§2.4 fixes the veto as a comparison", "b8b3131",
     ["docs/steps/corrections-and-dead-ends.md", "docs/TODO.md", "docs/project_progress.md",
      "docs/steps/05-multitask-results.md", "scripts/gate/retired_values.py", "scripts/gate/README.md"]),

    (r"normalized Spearman\s*\n?\s*\*{0,2}0\.98|Spearman\s+0\.98|ρ\s*=\s*\*{0,2}0\.98",
     "a synthetic mean-effects predictor's normalized score", "removed: no code produces it and no "
     "artifact records it", "f75e9bb",
     ["docs/TODO.md", "docs/steps/corrections-and-dead-ends.md", "report/sections/03_methods.tex",
      # dated removal notes that quote the figure they removed
      "scripts/evaluation/dreval_normalize.py", "scripts/archive/README.md",
      "notebooks/outputs/README.md", "notebooks/analysis/evaluation/dreval_benchmark.ipynb",
      "scripts/gate/retired_values.py", "scripts/gate/README.md"]),

    (r"0\.2450", "the lower bound of the eight-execution band",
     "0.2473, from the two replayable executions", "a2558d5",
     ["docs/steps/corrections-and-dead-ends.md", "docs/OPEN_DECISIONS.md",
      "docs/steps/05-multitask-results.md", "report/sections/06_limitations_and_outlook.tex",
      "scripts/evaluation/build_execution_band.py", "scripts/evaluation/first_fit_order_test.py",
      "scripts/gate/retired_values.py", "scripts/gate/README.md"]),

    (r"4,704", "the PREDICTED post-repair in-vocabulary gene count",
     "4,765, measured on the live h5ad", "9fafcff",
     ["docs/TODO.md", "docs/steps/corrections-and-dead-ends.md",
      "docs/steps/02-preprocessing-and-embeddings.md", "report/results_numbers.tex",
      # dated corrections that name 4,704 as the superseded prediction
      "scripts/preprocessing/gen_embeds.py", "scripts/training/cv.py",
      "notebooks/analysis/qc/verify_variants.ipynb", "scripts/gate/retired_values.py", "scripts/gate/README.md"]),

    (r"K\s*=\s*545|545\s+heads|545-head", "the catalogue standing in for the model's head count",
     "534 — the drugs that clear the 50-overlapping-line cut", "ef45d75",
     ["docs/steps/corrections-and-dead-ends.md", "docs/TODO.md", "docs/project_progress.md",
      "docs/steps/05-multitask-results.md", "docs/steps/03-model-and-training-design.md",
      # descriptions of the archived K=545 era: the 26.05 --all-drugs runs really did have 545 heads
      "docs/steps/01-datasets-and-harmonization.md", "docs/make_figures.py",
      "notebooks/README.md", "notebooks/outputs/README.md", "notebooks/4a_percell_training.ipynb",
      "notebooks/analysis/evaluation/dreval_benchmark.ipynb", "scripts/gate/retired_values.py", "scripts/gate/README.md"]),

    (r"180 trainable|180 labelled lines|153 of the 180",
     "the trainable line count before the H292 alias", "181", "914bb50",
     ["docs/steps/corrections-and-dead-ends.md", "docs/TODO.md",
      "docs/steps/01-datasets-and-harmonization.md", "docs/steps/03-model-and-training-design.md",
      "report/results_numbers.tex", "docs/steps/05-multitask-results.md", "docs/make_figures.py",
      "scripts/gate/retired_values.py", "scripts/gate/README.md"]),
]

SCAN_SUFFIXES = {".md", ".tex", ".py", ".ipynb"}
SKIP_DIRS = {".git", ".venv", ".claude", "node_modules", ".ipynb_checkpoints", "__pycache__",
             "archive"}
#: Untracked by design, or a dated historical record: not part of the live claim surface.
SKIP_FILES = {"docs/final_presentation.md", "docs/final_presentation_slides.md",
              "docs/final_presentation_supplement.md", "docs/gate5-rerun-report.md"}


def _text_of(p: Path) -> str:
    """Readable source of a file. For a notebook: **source cells only**, never stored outputs.

    A number printed into a committed cell output is the record of a past run, not a live claim, and
    matching it would make this check noisy — which is how a check stops being read. Line numbers are
    therefore reported against the concatenated cell sources, which is what a reader edits.
    """
    if p.suffix != ".ipynb":
        return p.read_text(errors="ignore")
    import json
    try:
        nb = json.loads(p.read_text(errors="ignore"))
    except Exception:
        return ""
    parts = []
    for c in nb.get("cells", []):
        src = c.get("source", "")
        parts.append(src if isinstance(src, str) else "".join(src))
    return "\n".join(parts)


def _files(root: Path):
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.suffix not in SCAN_SUFFIXES:
            continue
        if SKIP_DIRS & set(p.parts):
            continue
        rel = str(p.relative_to(root))
        if rel in SKIP_FILES or rel.startswith("docs/progress_report_"):
            continue
        yield rel, p


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    files = list(_files(root))
    findings, checked = [], 0
    for rel, p in files:
        text = _text_of(p)
        if not text:
            continue
        for pattern, was, replaced_by, commit, allowed in RETIRED:
            checked += 1
            if any(rel.startswith(a) for a in allowed):
                continue
            for m in re.finditer(pattern, text):
                line = text[:m.start()].count("\n") + 1
                findings.append(
                    f"RETIRED VALUE  {rel}:{line}\n"
                    f"    matched : {m.group(0)!r}\n"
                    f"    was     : {was}\n"
                    f"    now     : {replaced_by}   (retired in {commit})\n"
                    f"    -> either this is a live assertion and must be corrected, or this file\n"
                    f"       legitimately holds the record and belongs in that entry's `allowed` list")
    print(f"  {len(RETIRED)} retired values x {len(files)} files = {checked} checks, "
          f"{len(findings)} live occurrence(s)")
    for f in findings:
        print(f)
    if not files or not RETIRED:
        print("   ^ BLOCKER: nothing examined — '0 found' would be indistinguishable from '0 looked at'")
        return 2
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
