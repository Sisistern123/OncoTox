# OncoTox

Predicting pharmacological response (drug viability/toxicity) from single-cell RNA-seq by
mapping bulk CTRPv2 viability labels onto SCP542 cells, using **scGPT** foundation-model
embeddings as a denoised biological prior and comparing against a PCA baseline.

## Documentation

All project documentation lives in [`docs/`](docs/):

- **[docs/project_progress.md](docs/project_progress.md)** — the source of truth: an index
  (pipeline overview, full project arc, current-status scorecard, doc conventions) linking to
  thematic step files under [`docs/steps/`](docs/steps/) that hold the complete record of every
  step, number, parameter, and result, with alignment notes against the project plan. Start here.
- **[docs/steps/](docs/steps/)** — the eight self-contained step files (01–05 done, 06–08
  planned placeholders) covering the whole project structure end-to-end.
- **[docs/project_notes.md](docs/project_notes.md)** — dated thought/decision log
  (reasoning, advisor updates, ideas, open questions).

### Pipeline status at a glance

![OncoTox pipeline status overview](docs/pipeline_overview.png)

Regenerate the figure with `uv run docs/make_pipeline_overview.py`.

## Layout

```
scripts/preprocessing/   # SCP542 conversion, CTRP target mapping, splits, PCA, orchestrator
scripts/model/           # OncoMLP + datasets
scripts/training/        # train_multitask.py + shared training utils / run versioning
notebooks/               # overlap audit, scGPT vs PCA UMAPs, HVG-vs-all-genes comparison
runs/                    # per-run artifacts + runs_index.csv (gitignored)
docs/                    # project_progress.md (index) + steps/ (01-08), project_notes.md, figures
```

## Quickstart

```bash
# Preprocess (HVG-5000 variant, all CTRPv2 drugs; skips the external scGPT step if embeddings exist)
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --start-at targets --skip-scgpt --all-drugs

# Train multi-task (all 545 CTRPv2 drugs)
uv run scripts/training/train_multitask.py --use-rep X_scGPT     # scGPT embeddings
uv run scripts/training/train_multitask.py --use-rep X_pca       # PCA baseline
```

See [docs/project_progress.md](docs/project_progress.md) for full commands, data layout, and results.
