# OncoTox — Project Progress (index)

*Top-level index. The record itself lives in the step files under [`docs/steps/`](./steps/); this page
holds the pipeline overview, the project arc, the status scorecard and the doc conventions — and links
everywhere else. What went wrong is collected in
[Corrections](./steps/corrections-and-dead-ends.md); open work is in [TODO](./TODO.md).*

Reference plan: `~/Desktop/OncoTox/project_plan/project_planning_v2.pdf`.
Plan-alignment is marked **✅ on-plan** or **⚠️ deviation/addition** inside each step file.
A standalone LaTeX write-up of the current state lives in [`../report/`](../report/) (→ `main.pdf`).

> ## ⛔ 28.07.2026 — no number on this page is currently quotable
>
> **The drug panel was rebuilt on 12.08.2026, and every number computed on the old one is still void.**
> The previous panel's candidate list had been ranked on our own response values before the literature
> criterion was applied
> ([Corrections](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it)).
> The replacement is 11 drugs selected on FDA approval and verified published determinants
> ([Step 01](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026)),
> but **nothing has been re-run on it**, so no figure on this page is quotable until the sweep produces
> one. The review that must finish first is the [TODO](./TODO.md) checklist.
>
> **The report's numbers were withdrawn rather than caveated (12.08.2026).** `04_results.tex` is now a
> withdrawal note and `results_numbers.tex` defines data only; the values live in
> [Corrections](./steps/corrections-and-dead-ends.md).
>
> **Two earlier "current numbers" banners lived here and have been removed**, because both were built on
> since-retired foundations: the 13.07 box on `auc_z`
> ([retired](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target)) and the 14.07 box on
> the 10-drug panel ([superseded](./steps/corrections-and-dead-ends.md#the-1307-five-drug-numbers)).
> Their numbers are preserved in
> [Corrections](./steps/corrections-and-dead-ends.md) rather than deleted.
>
> **Scope reality check, which has not changed.** Everything trained so far is **one database, one
> response score** (CTRPv2). The 545-head "multi-task" run is multi-**drug**, *not* multi-database or
> multi-metric. The goal of combining CTRPv2 + PRISM + GDSC across efficacy *and* toxicity is
> [Step 06 · A](./steps/06-planned-work.md#a-cross-database-integration), and the reusable foundation model
> is [Step 06 · C](./steps/06-planned-work.md#c-foundation-model-and-clinical-fine-tuning). Neither has
> started. Do not read the current results as the finished goal.

---

## Pipeline overview (at a glance)

> ⛔ **The pipeline diagram is [archived](./figures/archive/), not shown here (12.08.2026).** It is
> derived from the targets h5ad and from a training run, and neither exists on the current target:
> preprocessing has not re-run under the [freeze](TODO.md), and the last training run used the void
> 8-drug panel. It is deliberately **not** rebuilt from the retired `auc` h5ad still on disk — a
> figure reproducible only from a target the pipeline no longer writes is not reproducible by a
> standard run. The same applies to `model_architecture.png`, `loss_02_weights.png` and
> `loss_03_effect.png`. All four return at R4; `make_figures.py` skips each with a printed reason
> until then.

What the pipeline does, stage by stage — the sparse (cell line × drug) response matrix, the drug panel
funnel, the cell-line-grouped folds, the two representations that are compared, the per-cell MLP,
the weighted loss, and the out-of-fold scoring — is in
[Step 03](steps/03-model-and-training-design.md), and the notebooks that run it are
[`notebooks/`](../notebooks/README.md) stages 1–5.

![OncoTox pipeline status overview](./figures/pipeline_overview.png)

Green = done / on-plan · amber = addition beyond plan · grey = results withdrawn · red (dashed) = still
missing. Stages 1–3 are complete; **stages 4 and 5 have had every model result withdrawn** (12.08.2026 —
the target was replaced, the panel rebuilt, and the representations predate the preprocessing
corrections), and they are re-measured at R4. The red boxes — cross-database PRISM/GDSC heads and the
XAI stretch goal — are the unstarted work. This figure is a pure drawing and carries no data, so it was
corrected and re-rendered on 12.08.2026 without touching a frozen artifact.

### The figure set

**Current — pure drawings, always reproducible:**

| figure | shows | argument lives in |
|---|---|---|
| `figures/pipeline_overview.png` | status against the written plan | this file |
| `figures/loss_01_objective.png` | what the objective is made of | working report §4 |

**[Archived](./figures/archive/) 12.08.2026 — derived from data, and not reproducible today:**

| figure | shows | why it is archived |
|---|---|---|
| `figures/archive/pipeline.png` | the pipeline, stage by stage | needs the targets cache *and* the panel CSVs |
| `figures/archive/model_architecture.png` | one cell in, one AUC per panel drug out | needs the targets cache |
| `figures/archive/loss_02_weights.png` | one drug's label density and the weight curve it produces | needs the targets cache |
| `figures/archive/loss_03_effect.png` | what the weighting did: spread up, ranking flat | needs a **training output**, so no h5ad can rebuild it |

The first three read `figure_data.npz`, rebuilt from the **`auc_cc`** targets h5ad — which has never
been written, because preprocessing has not re-run under the [freeze](TODO.md). They are deliberately
**not** built from the retired `auc` h5ad still on disk: a figure reproducible only from a target the
pipeline no longer writes is not reproducible by a standard run. `uv run docs/make_figures.py` builds
the two current figures and skips these four with a printed reason; every skip clears by itself at R4.
The drug panel is read from `notebooks/outputs/panel/panel.csv` rather than duplicated in
`make_figures.py`, which is how the old 8-drug copy went stale.

---

## The full project arc — document map

Each step is a self-contained file. **Steps 01–05 are done; 06–08 are placeholders** for planned
work, kept here so the entire project structure is visible end-to-end.

| Step | Status | What it covers |
|---|---|---|
| **[01 — Datasets & harmonization](./steps/01-datasets-and-harmonization.md)** | ✅ Done | Raw datasets (SCP542, CTRPv2, PRISM, GDSC), overlap/coverage audit, drug catalog, cell-line & compound harmonization. |
| **[02 — Preprocessing & embeddings](./steps/02-preprocessing-and-embeddings.md)** | ✅ Done | AnnData build, scGPT embeddings, UMAP latent validation, HVG-5000, `all_genes` variant, on-disk layout, reproduce commands. |
| **[03 — Model & training design](./steps/03-model-and-training-design.md)** | ✅ Done | Exact input/output/target/mask of a training example, MSE definition, **supervised** training paradigm. |
| **[04 — Single-task results](./steps/04-single-task-results.md)** | ✅ Done | Paclitaxel baseline + data-leak fix. **1 database, 1 score, 1 drug.** |
| **[05 — Multi-task results & versioning](./steps/05-multitask-results.md)** | ✅ Done | Masked-loss across 545 CTRPv2 drugs + run ledger. **Still 1 database, 1 score; multi-*drug* only.** The current Q1 result — capacity, head count, label supply and objective — is [§Q1 on the rebuilt panel](./steps/05-multitask-results.md#q1-on-the-rebuilt-panel--what-carries-the-pca-lead-13082026), which is **not** covered by that page's void banner. |
| **[06 — Cross-database integration](./steps/06-planned-work.md#a-cross-database-integration)** | ❌ Not started | **The "combine all" goal:** CTRPv2 + PRISM + GDSC, efficacy + toxicity, cross-database masked multi-task. |
| **[07 — XAI / feature interpretability](./steps/06-planned-work.md#b-xai-and-feature-interpretability)** | ❌ Not started | Stretch goal: feature importance → transcriptomic drivers of resistance. |
| **[08 — Foundation model & clinical fine-tuning](./steps/06-planned-work.md#c-foundation-model-and-clinical-fine-tuning)** | ❌ Not started | Overarching goal: reusable pan-cancer foundation model, fine-tunable on clinical (binary) outcomes. |

**Where this is going (the two axes that widen):**

```
Step 04   1 database · 1 score · 1 drug        (CTRPv2, paclitaxel)
Step 05   1 database · 1 score · K=545 drugs   (CTRPv2, all drugs)              ← here now
Step 06   3 databases · 2 metric types         (CTRPv2+PRISM+GDSC, efficacy+toxicity)
Step 08   + clinical fine-tuning               (continuous pre-train → binary clinical head)
```

**Fast facts you'll want regardless of which step you open** (full detail in
[Step 03](./steps/03-model-and-training-design.md)):

- A training example = **one single cell**; input = a 512-dim scGPT embedding (`X_scGPT`) or PCA
  (`X_pca`). Cell line / cancer type / drug are **not** input features.
- Target = a CTRPv2 response score chosen with `--score`, defined per **(cell line × drug)** and
  broadcast to every cell of that line. **Default since 11.08.2026: `auc_cc`** — the area under
  DrEval's CurveCurator re-fit of CTRPv2's raw dose-response data, in native viability units, with
  **no winsorization and no quality filter**
  ([Step 01](./steps/01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)).
  `ln_ic50_cc` is the alternative measure from the same fit and is not the default: it is undefined
  for ~40 % of curves by construction. Three earlier targets are
  [retired](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target) and none is comparable
  to `auc_cc` — `auc` (divided by the wrong quantity), `auc_z` (per-drug z-scored `auc`), and the
  legacy `mean_pv` that still backs Steps 04–05.
- Training is **fully supervised regression** (masked MSE/MAE). scGPT is a **frozen** self-supervised
  feature prior; the mask handles label sparsity but does **not** make it semi-supervised.

---

## Experiment matrix — PCA vs scGPT

The central comparison (the plan's core hypothesis) is run as a **2 × 2 × 2 = 8-run matrix**.
All eight runs share the cell-line-grouped split and a **matched trunk** `(128,64)` (set 14.06.2026)
+ identical training protocol; only the input representation (and its gene set) changes.

| Axis | Values |
|---|---|
| **Gene set** | `all_genes` (full transcriptome) · `hvg5000` (top-5,000 HVG from raw) |
| **Representation** | `X_pca` (standard single-cell PCA baseline, **512-d** to match scGPT) · `X_scGPT` (512-d embedding) |
| **Task** | single-task (paclitaxel) · multi-task (all drugs, K = 545) |

**Genes per condition** — PCA uses the full filtered set, scGPT only its in-vocabulary subset. Counts
per variant: [Step 02](./steps/02-preprocessing-and-embeddings.md#hvg-5000-pipeline-outputs).

> ⛔ **Corrected 05.08.2026.** This page previously called that gap **intentional** — "scGPT's
> vocabulary coverage is part of the model". It is not. **775 of the discarded genes are present in
> scGPT's vocabulary under their current symbols** and were thrown away by an exact match against an
> older annotation, costing 3.6 % of every cell's transcriptome:
> [Corrections](./steps/corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols).

> ⛔ **The 27.06.2026 result that stood here is superseded — twice over.** This page stated it without
> markers, which is exactly what an index must not do. Its conclusions were **overturned 13.07.2026**
> (they rest on the K=545 null the unstandardized loss produced) and the `all_genes` half was **struck
> 05.08.2026** (at `max_length=1200` scGPT never received the full transcriptome, so those rows compared
> two different gene sets). What survives, what does not, and the current numbers:
> [Step 05](./steps/05-multitask-results.md#multi-task-masked-loss-over-all-545-ctrpv2-drugs-26052026)
> and [Corrections](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions). **Do not quote a
> PCA-vs-scGPT result from this page** — it is an index and holds none of its own.

Results: [Step 04](./steps/04-single-task-results.md) (single-task), [Step 05](./steps/05-multitask-results.md)
(multi-task); per-drug coverage & learnability in `notebooks/analysis/harmonization/drug_coverage.ipynb`. Action list:
[TODO.md](./TODO.md).

> **PCA width matched to scGPT (27.06.2026).** The original matrix used a **~50-d** PCA (scanpy
> default) + a smaller `(64,32)` PCA trunk, both of which handicapped PCA. PCA now keeps **512
> components** (`add_pca.DEFAULT_N_COMPS`, override `--pca-n-comps`) on the matched `(128,64)` trunk,
> so PCA and scGPT share input width *and* parameter count — the last comparison confound is closed.
> The **full 8-run matrix was re-run at 512-d** (reproducible in `notebooks/4a_percell_training.ipynb (§B)`; run
> dirs `runs/20260627_1913xx_*`); the numbers above and in [Step 05](./steps/05-multitask-results.md)
> are these 512-d results and supersede the 14.06 (~50-d) matrix.

---

## Where everything is saved

The on-disk layout, the exact artifact for every representation / score / split, and how a run selects
among them: [Step 02](./steps/02-preprocessing-and-embeddings.md#current-data-layout-on-disk).
Run artifacts and the ledger: [Step 05](./steps/05-multitask-results.md#run-versioning-26052026).

## Notebooks and outputs

- **What each notebook does, and the order to read them:**
  [`notebooks/README.md`](../notebooks/README.md).
- **What each output directory holds and which notebook wrote it:**
  [`notebooks/outputs/README.md`](../notebooks/outputs/README.md).
- **The numbers themselves** are owned by the step files, not by any catalog.

Artifacts that current claims rest on, with the step that owns each number:

Paths updated 12.08.2026: everything a standard pipeline run cannot recreate moved under
`archive/`, and the rows below follow it. A `archive/` prefix means the producing notebook is archived
or its criterion was retracted — nothing regenerates those.

| Path (under `notebooks/outputs/`) | Owns the discussion |
|---|---|
| `data/target_distribution.png` | [Step 05](./steps/05-multitask-results.md) — re-run at R5, the line count moves 180 → 181 |
| `archive/drug_coverage.png` | [Step 05](./steps/05-multitask-results.md) — model output on a retired target; its cells were dropped from `drug_coverage` 13.08.2026, so nothing regenerates it |
| `embeddings/umap_cancertype_pca_vs_scgpt.png`, `umap_sweep_cancertype.png`, `variants.png` | [Step 02](./steps/02-preprocessing-and-embeddings.md#latent-space-validation-umap-fig-3--fig-4) |
| `dreval/dreval_lco*.{png,csv}`, `dreval/dreval_normalized*.csv` | [Step 05](./steps/05-multitask-results.md) — `dreval_benchmark` waits for R4's `outputs/panel/` files, then re-runs at R5 (*corrected 13.08.2026: this read "blocked on review item 11"; item 11's code fixes have landed and the notebook's imports all resolve*) |
| `diagnostics/*` | [Step 05](./steps/05-multitask-results.md), [Corrections](./steps/corrections-and-dead-ends.md) — waits for R4's `outputs/panel/` files, then re-runs at R5 |
| `archive/target/target_comparison.*`, `archive/target/loss_weighting_bug.png`, `archive/target/seed_stability.csv` | [Step 03](./steps/03-model-and-training-design.md), [Corrections](./steps/corrections-and-dead-ends.md) |
| `archive/ablations/rescue_k545.*`, `archive/ablations/ablation_*` | [Step 03](./steps/03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026) |
| `archive/learnability/*` | [Step 05](./steps/05-multitask-results.md) — criterion retracted |
| `archive/ctrp_drug_learnability_mean_pv.csv`, `archive/gdsc_drug_learnability.csv` | [Step 01](./steps/01-datasets-and-harmonization.md) — the GDSC list was never part of the modelling work |
| `panel/panel.csv`, `panel/literature_panel_candidates.csv` | **current** — the [rebuilt 11-drug panel](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026) and its 57 candidates |
| `panel/*` (everything else) | ⛔ computed on the [voided panel](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it) |
| `archive/training_545_mean_pv/*` | superseded — the `mean_pv` 8-run matrix, CV, per-drug ρ and gene-set sweep; numbers in [Step 05](./steps/05-multitask-results.md), status in [Corrections](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions) |

---

## The plan (for reference)

A staged prototype (from the plan PDF):

1. **Latent-space validation** — generate scGPT embeddings, compare to full-transcriptome
   PCA via UMAP (Fig. 3 by cancer type, Fig. 4 by paclitaxel viability); confirm scGPT
   removes tissue-of-origin bias.
2. **Single-task baseline** — regress the continuous CTRPv2 response
   score from the embeddings on the **highest-confidence intersection** SCP542×CTRPv2
   (**190 cell lines, 545 compounds, 100 % non-null in overlap**). *Do not start
   multi-task / PRISM / GDSC until this works.*
3. **Iterate outward** — add masked-loss multi-task and integrate the larger, sparser
   PRISM (and GDSC) datasets — efficacy **and** toxicity.
4. **Stretch goal** — XAI / feature importance.

Overarching main goal: a reusable pan-cancer single-cell **foundation model** fine-tunable for
specific cancer types / clinical (binary) datasets.

**Core hypothesis:** scGPT embeddings are a denoised biological prior that forces the
regressor to learn real resistance signatures instead of memorizing cell line / tissue
identity → should show as **less overfitting (smaller train/val gap) for scGPT than PCA**.

---

## Current status — plan vs. reality

> ⚠️ **Every verdict in this table rests on measurements withdrawn on 12.08.2026** and is pending
> re-measurement. The target was replaced, the drug panel rebuilt, and the representations on disk
> predate the preprocessing corrections, so a ✅ here means *"this was the finding on the pipeline as it
> then stood"*, not *"this currently holds"*. The report withdrew the corresponding numbers rather than
> caveating them (`report/sections/04_results.tex`); the verdicts are kept because a scorecard that
> silently emptied would lose the record of what was believed and why. Each is re-decided at **R6** of
> the sweep, from the regenerated artifacts.

| Plan item | Status | Where it is written up |
|---|---|---|
| Sub-goal 1: compound harmonization (names + BRD + DrugBank) | ✅ Done | [Step 01](./steps/01-datasets-and-harmonization.md) |
| Sub-goal 2: masked-loss sparsity handling | ✅ Done (intra-CTRPv2 only) | [Step 03](./steps/03-model-and-training-design.md#mask-m--the-sparsity-handling-mechanism-plan-sub-goal-2) |
| Sub-goal 3: baseline on the SCP542×CTRPv2 intersection | ✅ Done | [Step 04](./steps/04-single-task-results.md), [Step 05](./steps/05-multitask-results.md) |
| Phase 1: scGPT embeddings + UMAP latent validation | ✅ Done | [Step 02](./steps/02-preprocessing-and-embeddings.md#latent-space-validation-umap-fig-3--fig-4) |
| Phase 2: single-task continuous regression | ✅ Done, on legacy `mean_pv` | [Step 04](./steps/04-single-task-results.md) |
| Phase 3a: multi-task masked loss | ✅ Done, **CTRPv2 only** | [Step 05](./steps/05-multitask-results.md) |
| Phase 3b: PRISM / GDSC cross-database | ❌ Not started — data harmonized only | [Step 06 · A](./steps/06-planned-work.md#a-cross-database-integration) |
| Stretch: XAI / feature importance | ❌ Not started | [Step 06 · B](./steps/06-planned-work.md#b-xai-and-feature-interpretability) |
| Main goal: foundation model + clinical fine-tuning | ❌ Not started (horizon) | [Step 06 · C](./steps/06-planned-work.md#c-foundation-model-and-clinical-fine-tuning) |
| Core hypothesis: scGPT overfits less than PCA | ✅ Confirmed — **generalization only**, not accuracy | [Step 05](./steps/05-multitask-results.md) |
| Does the model rank cell lines at all? | ✅ Yes, on drugs that carry signal | [Step 05](./steps/05-multitask-results.md); the earlier "no" is [superseded](./steps/corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result) |
| Is scGPT's lead over PCA real? | 🟡 Sign-consistent, **not a proven margin** — single seed, on the since-voided panel | [Step 05](./steps/05-multitask-results.md), [TODO](./TODO.md) (seeds are blocking) |
| Does the deep single-cell apparatus beat ridge on line means? | 🟡 Open — the measured margin is void, and only the network arm was flattered | [Step 03](./steps/03-model-and-training-design.md#the-baseline-that-actually-binds-ridge-on-150-line-mean-embeddings) |
| Is the model over-regularized / too big? | ❌ No — model-side tuning is **closed** | [Step 03](./steps/03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026), [Corrections](./steps/corrections-and-dead-ends.md#the-model-is-over-regularized-or-too-small) |
| External benchmark (DrEval LCO, normalized) | ✅ Above naive, below best-in-class | [Step 05](./steps/05-multitask-results.md); the first run's leak in [Corrections](./steps/corrections-and-dead-ends.md#the-first-dreval-benchmark--a-val-split-leak) |
| Does the signal survive removing the cell-line effect? | 🟡 Open — the measurement is void, and the metric that produced it does not test this under LCO | [Corrections](./steps/corrections-and-dead-ends.md#the-step-1-training-run-on-the-voided-panel) (measured on the voided panel) |
| Is the cell-line effect proliferation? | ❌ No — tested and refuted | [Corrections](./steps/corrections-and-dead-ends.md#the-cell-line-effect-is-largely-proliferation) |
| Does inverse-density loss weighting help? | ❌ No — clean negative | [Corrections](./steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking) |
| Was `auc_z` the right target? | ❌ No — retired 27.07.2026 | [Corrections](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target) |
| Is the drug-selection criterion sound? | ❌ No — it measured potency, not rankability | [Corrections](./steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability) |
| Drug selection | ✅ **Rebuilt 12.08.2026 — 11 drugs, FDA approval + published determinants**; nothing has been run on it yet | [Step 01](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026); the voided predecessor is in [Corrections](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it) |
| 8-run matrix conclusions | ⚠️ Suspect — re-run required, expect them to change | [Corrections](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions) |
| Can the setup test research question 2 (implicit heterogeneity)? | ❌ **Not as built** — the objective penalizes it | [Step 03](./steps/03-model-and-training-design.md), and *Where this goes next* below |

**Additions beyond the written plan (all defensible — document them):** random→leak→grouped
split ([Step 04](./steps/04-single-task-results.md)); HVG-5000 + all-genes comparison
([Step 02](./steps/02-preprocessing-and-embeddings.md)); per-drug-mean sanity baseline +
run-versioning ledger ([Step 05](./steps/05-multitask-results.md)); **512-d PCA** (input width
matched to scGPT, removing the dimensionality confound — [Step 05](./steps/05-multitask-results.md));
**reproducible preprocessing/training notebooks** (see [`notebooks/README.md`](../notebooks/README.md)).

**Two things to flag clearly in the writeup:**

1. **Multi-task today = 545 CTRPv2 drugs, not CTRPv2+PRISM+GDSC** — plan-Phase-3 half done
   (the real "combine all" is [Step 06](./steps/06-planned-work.md#a-cross-database-integration)).
2. **Cell-line overlap: 190 vs 180** — 190 = name matches in CTRPv2's roster; 180 = lines with
   actual post-QC measurements (10 listed-but-unscreened lines drop out). It's **data availability,
   not normalization** (verified 14.06). Use 180 (the trainable set) — **181 after the next sweep**,
   since the join audit found one screened line that name matching had been dropping
   ([Step 01](./steps/01-datasets-and-harmonization.md#the-join-dropped-a-screened-cell-line-h292-10082026),
   10.08.2026).

---

## Where this goes next, and why in that order

*The governing rule: never change the target and the architecture in the same run. That is what made the
June result take weeks to unpick, and it is why the 27.07 step moved only the target and the loss.*

**1. MIL — next, and it is `4b_mil_training`.** Not one more architecture to try, but the **minimal
change that makes research question 2 askable at all**: today the constant-within-line label means the
objective penalizes any difference the model predicts between two cells of one line. A bag constrains
only the aggregate. It removes the line-weighting artifact structurally (one bag = one line = one
example), and it makes the per-cell prediction the clinically interesting readout — which subpopulation
drives the response. It needs no new data.

> ⛔ **Two claims were removed from the paragraph above (12.08.2026), both already retracted in
> [TODO](./TODO.md) on 11.08.2026 and left standing here.** MIL was called *"the only untested capacity
> lever"* — it is not a capacity lever at all, and there is no untested one; that phrasing survives from
> the Q1 framing and implies a reserve of unexplored performance that does not exist. And the success
> criterion was given as *"it has to beat the per-cell MLP and ridge on line means, and failing both is
> itself a reportable result"*, which scores a Q2 experiment on a Q1 criterion. The controls are a
> **floor, not the criterion**. The 82× figure was also dropped from this sentence — it is owned by
> [Step 03](./steps/03-model-and-training-design.md), and this page must not restate it.
>
> ⚠️ **This heading read *"MIL / attention pooling"*, and that is now the wrong branch (`58fadd7`).**
> The design is **instance-level** MIL: every cell carries its own predicted response rather than an
> attention weight over a pooled embedding — readable at the individual cell, at an expected cost in
> predictive accuracy, taken deliberately because Q2 is a question about readability. **And the success
> criterion is no longer open**, contrary to what this page said an hour earlier: it is pre-registered
> in [`4b_mil_training.ipynb`](../notebooks/4b_mil_training.ipynb) §2 as a synthetic positive control
> (precondition), within-line spread (necessary condition), cross-seed reproducibility against a
> shuffled-cell control (**the test**) and confound regression (**veto**). One number is outstanding —
> `Q2_CONTROL_THRESHOLD` — and it is the single blank keeping 4b a stub.

**2. scDEAL-style bulk pretraining + more cell lines — after MIL.** The remaining levers are argued to
be label-side. ⛔ **The three measurements behind that argument are void (12.08.2026), for two different
reasons.** *Model tuning is closed* and *ridge ties the MLP* come from runs whose checkpoint was chosen
on the fold it was scored on, on a retired target and a voided panel — re-derived minimally at R4
([TODO](./TODO.md) item 8C). *The density weighting was a null* is **not** one of those runs and fails
differently: audit 09 found the metrics it was judged on could not have seen the effect it was used to
rule out, so it is **re-tested rather than retired**, with `alpha` swept over {off, 0.5, 1.0} as an arm
of the loss comparison ([TODO](./TODO.md) item 9A). What survives is the
structural argument, which needs no run: the label is per cell line, so the independent sample size is
the line count and not the cell count. The screens themselves are much larger —
CTRPv2 ~1,100 lines, GDSC ~970, PRISM ~900 — so the labels exist and single-cell expression is what is
missing. scDEAL pretrains a denoising autoencoder on bulk and aligns the bulk and single-cell latent
spaces by domain adaptation, which attacks that gap directly instead of copying one bulk value onto ~300
cells. It comes second because it changes *where the representation comes from* while MIL changes *how
cells map to a line-level prediction* — landing both at once makes the result unattributable — and
because MIL's outcome decides whether the single-cell framing is worth building on.

**Deferred with reasons, not as a backlog:** the base quantity (EC50/Emax instead of AUC — a second
target change, would collide with MIL); seeds (preliminary results with a fixed seed are acceptable until
a margin is quoted as a number); the input-scale asymmetry between the two arms under one learning rate
(cheap, but it qualifies an old conclusion rather than advancing a new one — and its size has to be
re-measured before it can be argued from, since `sc.pp.scale` entered the PCA path on 05.08.2026 and
standardizing genes moves component magnitudes; the ~78× this line used to quote is
[stale](./TODO.md#after-the-sweep--the-one-review-item-that-needs-new-runs), and the naive
"different effective step size" reading of it does not survive LayerNorm plus Adam); learned task
weights (they estimate residual variance, mixing label noise with model error).

---

## Maintaining these docs (conventions)

- **Source of truth = [`steps/`](./steps/).** The step files hold every number, parameter and result;
  this page is an **index and scorecard only** and must not become a second copy of them. Anything
  superseded, retracted, refuted or abandoned goes in
  [Corrections](./steps/corrections-and-dead-ends.md) — the steps carry a pointer, not the content.
- **When new work lands:** update the relevant **step file** with the hard details (cell/gene
  counts, split distributions, hyperparameters, run IDs, MSEs, deviations), then refresh the
  **scorecard + arc** above if the plan-vs-reality picture changed. Each step carries its own
  ✅ on-plan / ⚠️ deviation callouts against the plan PDF.
- **Never state a number in two places.** If it belongs in a step file, this page links to it. The one
  exception is a superseded number, which lives only in [Corrections](./steps/corrections-and-dead-ends.md).
- **Every claim names the code that produced it** — the script and function, or the notebook and section,
  plus the `outputs/` artifact a number was read from.
- **The 190-vs-180 cell-line overlap** is resolved and written up where it belongs, in
  [Step 01](./steps/01-datasets-and-harmonization.md) — with the roster-vs-screened distinction and the
  ten unscreened lines. It used to be restated here in full, which contradicted the two rules above.
