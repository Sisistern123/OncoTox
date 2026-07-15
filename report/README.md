# Preliminary Report (LaTeX) — Prediction of Anti-Cancer Drug Efficacy Scores

Modular LaTeX source for the preliminary project report. Same build model as the CV repo:
plain `pdflatex`, one `main.tex`, sections `\input` from `sections/`. No `bibtex`/`biber`
run needed — references are a manual `thebibliography` in `sections/07_bibliography.tex`.

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
`results_numbers.tex`, and the section text references the macro (e.g. `\rhoFullScgpt`,
`\ridgePca`) instead of a literal digit. `main.tex` does `\input{results_numbers}` in the
preamble, so a number appears in exactly one place and editing it there updates every
mention. This mirrors the ToxFam manuscript's `results_numbers`, with one difference: here
the file is **maintained by hand** (there is no extraction script yet), so each macro carries
a comment naming the source CSV under `notebooks/outputs/` and, for aggregates, how the value
is reduced (e.g. mean over per-drug or per-fold rows). To update after re-running a notebook:
read the value from the refreshed CSV, edit the macro, recompile. The two results tables carry
their full per-cell values inline, since a table is itself a single location.

## Render the PDF

Needs a TeX distribution with `pdflatex` (e.g. BasicTeX, as in the CV repo). The only
non-default package is `booktabs`; install once if missing:

```sh
sudo tlmgr update --self
sudo tlmgr install booktabs
```

Then, from this `report/` directory:

```sh
pdflatex main.tex   # run twice so \ref / \cite cross-references resolve
pdflatex main.tex
```

Output: `main.pdf`. (On Overleaf/ShareLaTeX it renders directly, like the CV.)

## Updating figures

The figures are copies from the notebook outputs — refresh them with:

```sh
cp ../notebooks/outputs/embeddings/umap_cancertype_pca_vs_scgpt.png figures/fig_umap.png
cp ../notebooks/outputs/ablations/rescue_k545.png                   figures/fig_rescue.png
cp ../notebooks/outputs/dreval/dreval_lco.png                       figures/fig_dreval.png
```

Provenance of every number: `../docs/project_progress.md`, `../docs/steps/`, `../docs/TODO.md`,
and `../notebooks/outputs/*.csv`.
