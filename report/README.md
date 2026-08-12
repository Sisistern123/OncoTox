# Preliminary Report (LaTeX) — Prediction of Anti-Cancer Drug Efficacy Scores

Modular LaTeX source for the preliminary project report: one `main.tex`, sections `\input`
from `sections/`, and **plain `bibtex`** — `sections/07_bibliography.tex` is
`\bibliographystyle` + `\bibliography{../references}`, generated from the project's single
`references.bib`, not a hand-written `thebibliography`. See *Render the PDF* for the build
command; `pdflatex` alone is not enough.

The structure follows a standard scientific manuscript (Abstract, Introduction, Methods,
Results, Discussion, Limitations & Future Work), matching the ToxFam manuscript layout.
`project_planning_v2.pdf` is kept alongside as the baseline reference.

## Layout

```
main.tex               # preamble + \input of results_numbers and each section
results_numbers.tex    # single source of truth for the headline numbers (see below)
sections/              # 01_abstract, 02_introduction, 03_methods, 04_results,
                       #   05_discussion, 06_limitations_and_outlook, 07_bibliography
figures/               # PNGs copied from notebooks/outputs/ (self-contained)
project_planning_v2.pdf# the project plan, kept as the baseline to orient by
```

## The numbers: `results_numbers.tex`

Every headline figure cited in the prose lives once, as a `\newcommand` macro in
`results_numbers.tex`, and the section text references the macro (e.g. `\NLines`, `\NDrugs`)
instead of a literal digit. `main.tex` does `\input{results_numbers}` in the
preamble, so a number appears in exactly one place and editing it there updates every
mention. This mirrors the ToxFam manuscript's `results_numbers`, with one difference: here
the file is **maintained by hand** (there is no extraction script yet), so each macro carries
a comment naming the source CSV under `notebooks/outputs/` and, for aggregates, how the value
is reduced (e.g. mean over per-drug or per-fold rows). To update after re-running a notebook:
read the value from the refreshed CSV, edit the macro, recompile.

⛔ **As of 12.08.2026 `results_numbers.tex` defines DATA ONLY.** Every macro carrying a model
result was removed when the results were withdrawn, and the section text with it — so the
`\rhoFullScgpt` / `\ridgePca` this paragraph used to name as examples no longer exist, and the
"two results tables carrying their per-cell values inline" it used to describe are gone with
`04_results.tex`, which is now a withdrawal note. Results macros are **not** re-added by hand:
TODO item 12 replaces this file with an extraction script, and hand-transcription is the defect
that item exists to remove. Read the withdrawal block at the foot of `results_numbers.tex`
before adding anything to it.

## Render the PDF

Needs a TeX distribution with `pdflatex` (e.g. BasicTeX, as in the CV repo). The only
non-default package is `booktabs`; install once if missing:

```sh
sudo tlmgr update --self
sudo tlmgr install booktabs
```

Then, from this `report/` directory:

```sh
pdflatex main.tex   # writes main.aux with the \cite keys
bibtex   main       # resolves them against ../references.bib -> main.bbl
pdflatex main.tex   # pulls main.bbl in
pdflatex main.tex   # settles \ref / page numbers
```

Output: `main.pdf`. (On Overleaf/ShareLaTeX it renders directly, like the CV.)

⚠️ **`bibtex` is not optional and its omission is invisible in a tree that has been built
before.** `report/*.bbl` is gitignored, so in a fresh clone or a git worktree the
`pdflatex`-twice recipe this file used to document yields **33 undefined citations and a
bibliography-free PDF**; where a stale `main.bbl` is lying around it looks correct. Corrected
in `.gitignore` on 12.08.2026 (`cf3ad3f`) and here on the same day. Verified by building this
worktree, which had no `.bbl`: `pdflatex` alone reported `No file main.bbl`, and the four-step
sequence above gives 0 undefined references over 17 pages.

**`main.pdf` is not in git** (untracked since 10.08.2026): it is a build product, and while it was
committed it went stale every time a `.tex` changed without someone remembering to rebuild. Build it
from the source above. `project_planning_v2.pdf` *is* tracked — it is a source document, not something
this repository produces.

## Updating figures

The figures are copies from the notebook outputs. **Only one is used by the report** —
`fig_umap.png`, in `02_introduction.tex`; it is the sole `\includegraphics` in the whole
document. Refresh it with:

```sh
cp ../notebooks/outputs/embeddings/umap_cancertype_pca_vs_scgpt.png figures/fig_umap.png
```

⚠️ **`figures/fig_rescue.png` and `figures/fig_dreval.png` are orphans** — referenced by no
`.tex` file, and both left over from the results withdrawn on 12.08.2026. Their sources are
dead too: `fig_rescue.png` came from `outputs/legacy/ablations/rescue_k545.png`, whose notebook
is archived and cannot be re-run, and `fig_dreval.png` from `outputs/dreval/dreval_lco.png`,
computed on the retired target and the voided panel. Whether they return, and from which
artifact, is decided when §Results is rewritten at **R6**; refresh commands for them are not
documented here in the meantime, because a `cp` from a stale source is how a withdrawn figure
walks back into a document.

Provenance of every number: `../docs/project_progress.md`, `../docs/steps/`, `../docs/TODO.md`,
and `../notebooks/outputs/*.csv`.
