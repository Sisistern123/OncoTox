# Step 07 — XAI / feature interpretability (stretch goal)

> **Status: ❌ NOT STARTED — placeholder.** This is the plan's explicit **stretch goal**
> (time permitting), gated behind a stable, working predictor. Documented here so the full
> project structure is visible.

*Part of [OncoTox project progress](../project_progress.md). Covers (future): Explainable AI
methods to extract feature importance from a trained predictor and bridge predictions back to
biological mechanism.*

---

## Goal (from the project plan)

> "If the baseline regression model is successfully established and time allows, we will employ
> Explainable AI (XAI) methods to extract feature importance. This will allow us to bridge the
> predictive model back to underlying biological mechanisms by highlighting the key transcriptomic
> drivers of drug resistance."

## What this needs first

A stable, trained predictor worth interpreting — a checkpoint `runs/<…>/best_model.pt` loaded into
`OncoMLP` (`scripts/model/OncoMLP.py`), with cells fed via `scripts/model/dataset.py`. Ideally the
cross-database model ([Step 06](06-cross-database-integration.md)), but the existing CTRPv2
multi-task model ([Step 05](05-multitask-results.md)) is already a valid target. The analysis itself
would most naturally live in a new `notebooks/` notebook, alongside the existing
`notebooks/scgpt_umap.ipynb` / `notebooks/analysis.ipynb`.

## Open design questions (to resolve when starting)

- **Attribution target.** The model consumes **embeddings** (`X_scGPT` 512-dim or `X_pca`), not raw
  genes, so importance is over embedding dimensions. The hard part is mapping embedding-dim
  importance back to **transcriptomic drivers** — e.g. propagating attributions through to input
  genes, or correlating salient embedding dims with known gene programs.
- **Method:** gradient-based attribution (Integrated Gradients / saliency), SHAP, or per-head
  permutation importance.
- **Per-drug-head interpretation:** which transcriptomic signatures drive resistance for a specific
  compound.

## Definition of done

- Feature-importance attributions produced for ≥1 drug head and linked to plausible biology.
- Method + findings documented here.
