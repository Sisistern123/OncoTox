# Step 06 — Planned work: cross-database integration, XAI, foundation model

> **Status: ❌ NONE OF THIS HAS STARTED.** This file documents the three stages that remain, so the whole
> project structure is visible end-to-end. It holds intended scope, not results. It replaces the former
> `06-cross-database-integration.md`, `07-xai-feature-interpretability.md` and
> `08-foundation-model-and-clinical-finetuning.md`, which were three separate placeholders with the same
> skeleton.

*Part of [OncoTox project progress](../project_progress.md). Current work is in
[Steps 01–05](../project_progress.md); what went wrong is in
[Corrections](corrections-and-dead-ends.md); the action list is [TODO](../TODO.md).*

The three stages are strictly ordered: **A** widens the label space, **B** interprets whatever **A**
produces, **C** turns it into something reusable. Each depends on the one before, which is why none of
them is a candidate for "next" while the panel and target are still being settled.

---

## A: Cross-database integration

*CTRPv2 + PRISM + GDSC, efficacy and toxicity.* This is the plan's actual **Goal Option B** — *"simultaneously predict multiple response metrics
(efficacy and toxicity) across pan-cancer cell lines using a multi-task learning setup (via masked
losses) to handle sparse or missing labels"* — and the step where "combine all" actually happens. Data
is downloaded and harmonized ([Step 01](01-datasets-and-harmonization.md)); no integrated training has
run.

**Why Steps 04–05 are not already this.** The work so far is narrow on two axes that this stage widens:

| Axis | Steps 04–05 (done) | Stage A (planned) |
|---|---|---|
| **Database** | CTRPv2 only | CTRPv2 **+ PRISM + GDSC** |
| **Response metric** | one CTRPv2 score | viability **+ LN_IC50 / AUC + toxicity** |
| **What "multi-task" means** | multi-**drug** (545 heads, one metric) | multi-**database / multi-metric** heads |

So the 545-head run in [Step 05](05-multitask-results.md) is multi-task across the drugs of a single
source. It validates the masked-loss machinery on intra-CTRPv2 sparsity and nothing more.

**Planned approach.**

1. **PRISM first** — the plan's explicit next move, and the larger, far sparser dataset (915 lines,
   6,575 compounds, ~29 % non-null). Needs a PRISM/GDSC analogue of
   `scripts/preprocessing/ctrp_to_h5ad.py` emitting additional `Y_*` / `M_*` blocks (or extending
   `Y_ctrp`/`M_ctrp` into a unified matrix), joined via the harmonized drug catalog in `data/drug/`
   (name + BRD-ID + DrugBank — [Step 01](01-datasets-and-harmonization.md)) to map PRISM compounds onto
   existing CTRPv2 heads and add PRISM-only heads. `scripts/layout.py` gains the new
   source files.
2. **Then GDSC** (`LN_IC50` / AUC) — a *different metric type*, so this is where heads stop being
   homogeneous and the model becomes genuinely multi-metric.
3. **Masked loss across the union** — generalize the `MultiDrugDataset` mask machinery
   (`scripts/model/dataset.py`) and the masked loss ([Step 03](03-model-and-training-design.md)) to a
   block-sparse label matrix spanning all sources, where each (cell line × drug × metric) entry is
   observed in only some databases, with per-database/per-metric head grouping and weighting in
   `scripts/training/train_multitask.py`.
4. **Cross-database splits** — keep the cell-line-grouped, leakage-free discipline of
   `scripts/preprocessing/create_splits.py` across the unified cell-line set.

**Open design questions.**

- One shared trunk with per-database/per-metric head groups, or per-metric normalization before a shared
  head? The metrics live on different scales — viability ≈ [0,1], `LN_IC50` unbounded.
- How to weight databases and metrics given very different coverage (CTRPv2 21 % vs PRISM 29 % vs GDSC
  3 % non-null). Ties into the open per-head / uncertainty weighting question.
- Harmonize metrics onto a common response scale, or keep them as separate heads?
- Use the cell lines with actual measurements, not the 190 roster name-matches, before unioning
  sources — see [Step 01](01-datasets-and-harmonization.md). **That count is 181**, not the 180 this
  line carried until 13.08.2026: `ctrp_to_h5ad` applies a sourced `h292 → ncih292` alias
  (Cellosaurus `CVCL_0455`, audit 02) which recovers one screened line. Confirmed by running the
  `targets` step — 198 lines in the atlas, 181 with CTRPv2 labels, 17 without.

**Definition of done.** A single model with masked heads spanning ≥ 2 databases trains leakage-free;
per-database / per-metric "heads beating baseline" reported on the same honest metric as
[Step 05](05-multitask-results.md); result documented with run IDs in the versioning ledger.

---

## B: XAI and feature interpretability

*A stretch goal, gated behind a stable predictor.* The plan, time permitting: *"If the baseline regression model is successfully
established and time allows, we will employ Explainable AI (XAI) methods to extract feature importance.
This will allow us to bridge the predictive model back to underlying biological mechanisms by
highlighting the key transcriptomic drivers of drug resistance."*

**Gated, and the gate is quantitative.** Gene-level XAI runs only once per-drug ρ is substantially
higher and stable — below that it interprets noise. A cheaper *diagnostic* form is available now and is
tracked separately in [TODO](../TODO.md): where errors concentrate (drugs, tissues, lines with ρ < 0)
and how much residual error is the cell-line effect, using the existing per-drug-ρ CSVs.

**What it needs first.** A stable, trained predictor worth interpreting — a checkpoint
`runs/<…>/best_model.pt` loaded into `scripts/model/OncoMLP.py`, with cells fed via
`scripts/model/dataset.py`. Ideally the Stage A model, though the existing CTRPv2 multi-task model
([Step 05](05-multitask-results.md)) is already a valid target. The analysis would live in a new
notebook.

**Open design questions.**

- **Attribution target.** The model consumes **embeddings** (`X_scGPT` 512-d or `X_pca`), not raw genes,
  so importance is over embedding dimensions. The hard part is mapping that back to **transcriptomic
  drivers** — propagating attributions through to input genes, or correlating salient dimensions with
  known gene programs.
- **Method** — gradient-based attribution (Integrated Gradients / saliency), SHAP, or per-head
  permutation importance.
- **Per-drug-head interpretation** — which transcriptomic signatures drive resistance for a specific
  compound.

**Definition of done.** Feature-importance attributions produced for ≥ 1 drug head and linked to
plausible biology, with method and findings documented.

---

## C: Foundation model and clinical fine-tuning

The project's **main goal** as stated in the plan, and the reason every earlier step is built the way it
is: *"To develop a domain-specific pan-cancer single-cell foundation model capable of predicting
pharmacological response (efficacy/toxicity), which can subsequently be fine-tuned for specific cancer
types and/or clinical datasets (binary clinical outcomes)."*

**Why the continuous bulk labels are the right pre-training signal.** CTRPv2/PRISM/GDSC responses
capture a richer distribution than binary cutoffs, so a model pre-trained on them should yield a
representation that fine-tunes well onto scarce clinical data with binary responder / non-responder
outcomes.

**Why it sits on top of Stages A–B.** The whole pipeline — frozen scGPT prior → supervised regression
head, masked cross-database multi-task — exists to produce a **transferable representation**, not just a
CTRPv2 predictor. The inference payoff the plan describes: run the trained model on individual cells to
predict a **distribution of sensitivities within one sample**, flagging rare naturally-resistant
sub-clones computationally before selection pressure acts on them. That is the direct link to how
baseline tumour heterogeneity drives resistance, and it is the same motivation as research question 2.

**What it needs first.**

- A stable cross-database multi-task model (Stage A) — its `runs/<…>/best_model.pt` becomes the
  pre-trained **trunk**, with `scripts/model/OncoMLP.py` given a swappable clinical (binary) head.
- A fine-tune entrypoint alongside `scripts/training/train_multitask.py`, plus a clinical dataset loader
  analogous to `scripts/model/dataset.py`.
- **Access to a clinical dataset with binary outcomes** — not in the project, and the plan itself notes
  large-scale clinical data as the standing bottleneck.

**Open questions (long horizon).** Freeze trunk and retrain the head, or full fine-tune on clinical
data? How to bridge in-vitro cell-line training to in-vivo clinical tumours (domain shift)? Single-cell
inference plus sub-clone distribution analysis as a downstream evaluation.

**Definition of done.** Pre-trained checkpoint packaged as reusable starting weights, and at least one
fine-tuning / transfer demonstration documented.
