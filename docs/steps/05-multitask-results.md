# Step 05 — Multi-task masked loss (545 CTRPv2 drugs) & run versioning

*Part of [OncoTox project progress](../project_progress.md). Covers: the multi-task masked-loss
model over all 545 CTRPv2 drugs, its results vs. the per-drug-mean baseline, and the run-versioning
ledger that records every training run.*

> # ⛔ 13.08.2026 — NO NUMBER ON THIS PAGE IS CURRENTLY SUPPORTED
>
> Every result here predates the pipeline review and rests on at least one of: a **retired target**
> (`mean_pv`, `auc`, `auc_z` — the target has been `auc_cc` since 11.08.2026), the **void 8-drug panel**
> (rebuilt to 11 drugs on 12.08.2026), or the **early-stopping leak** fixed 12.08.2026, under which
> every out-of-fold prediction and CV metric on record is a minimum over epochs on its own scored data.
> Do not quote from this page, and do not compare a regenerated number to one here without saying which
> of those three the old one carried.
>
> ⚠️ **The DrEval benchmark section is named separately because it reads as a standing external
> validation and carries no marker of its own.** Its normalized ρ **0.511 ± 0.085** and R² **0.224**,
> and the claims built on them — *"clears `NaiveMeanEffects` decisively"*, *"directly comparable to
> DIPK 11 % and Random Forest 19 %"*, *"scGPT > PCA is confirmed externally"* — are void on two further
> grounds. They were measured on the retired `auc` target; and **since `66442d2` the benchmark runs on
> the rebuilt 11-drug panel, which shares exactly one drug — `dasatinib` — with the ten it used.** That
> is a different benchmark, not a refreshed one, so a future DrEval number is not comparable to this one
> at all. The *"+0.075 ± 0.038 we measured ourselves"* that the section calls an independent replication
> was itself **withdrawn** (`1971990`).
>
> **This lifts per section as R4–R6 regenerate what each rests on — it is not one event.** The sections
> already carrying their own dated banners keep them; this one covers what they do not.

This moves from plan-Phase-2 (single-task) into plan-Phase-3 (masked-loss multi-task). Masked-loss
mechanics are in [Step 03](03-model-and-training-design.md). These are the **multi-task (all-drugs,
K=545) rows of the 8-run experiment matrix**
([index](../project_progress.md#experiment-matrix--pca-vs-scgpt)).

> **Scope — still 1 database, 1 score; "multi-task" here = multi-*drug* only.** Every one of the
> K=545 heads predicts the **same** metric from the **same** database
> (CTRPv2). This validates the masked-loss machinery on intra-CTRPv2 sparsity, but it is **not**
> the plan's ultimate multi-task goal, which is **cross-database** (CTRPv2 + PRISM + GDSC) **and
> multi-metric** (efficacy *and* toxicity). That integration — the real "combine all" — is
> [Step 06](06-planned-work.md#a-cross-database-integration) and is **not yet started**. Do not read the 545-head
> run as "multi-task complete."

> ⛔ **Every number on this page is void, and none of it can be regenerated.** All of it was trained on
> **`mean_pv`**, the only target until 13.07.2026. That measure — along with `auc` and `auc_z` — was
> removed with its reader code on 11.08.2026 when the target moved to DrEval's reprocessed CTRPv2
> ([Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)),
> so `--score mean_pv` now raises. This page's earlier instruction to reproduce it that way was not
> merely stale but false.
>
> **The content is kept until the sweep regenerates it** (decided 11.08.2026, Selin): the run-versioning
> ledger, the split and drug-scope design, the panel and DrEval sections carry structure that is not
> just numbers, and gutting it now means rebuilding it in a fortnight. Read every figure on this page as
> *what was measured then*, never as a current result. Why the two scores are not interchangeable, and
> what does and does not transfer across them:
> [Corrections](corrections-and-dead-ends.md#the-steps-0405-numbers-as-a-comparable-baseline).

---

## Q1 on the rebuilt panel — what carries the PCA lead (13.08.2026)

> ✅ **This section is not covered by the page banner above.** It was measured after the pipeline
> review, on the rebuilt 11-drug panel, the `auc_cc` target and the corrected early stopping — the
> three grounds the banner names. Every figure below is read from a committed artifact and names it.

**What was asked.** Section A of `notebooks/4a_percell_training.ipynb` put `X_pca` ahead of
`X_scGPT`. Two explanations had to be excluded before that can be called a property of the
representations: that the shared trunk was doing the work, and that the eleven-drug panel is a
scoring set that happens to favour PCA. Sections C and D test one each, changing nothing else.

### C · Capacity does not carry it

`4a` §C, `RUN_SECTION_C = True` → `notebooks/outputs/panel/panel_arch_summary.csv`.
Grid: two architectures x two representations x three seeds x five folds = 60 fits, at `alpha=0`,
`loss=mse`. Bands are the seed half-range; "mean" is the mean per-drug out-of-fold Spearman over the
eleven panel drugs.

| arch | rep | mean | band | vs ridge |
|---|---|---|---|---|
| linear | `X_pca` | 0.2608 | ±0.0072 | **−0.0158** |
| linear | `X_scGPT` | 0.2291 | ±0.0019 | +0.0378 |
| trunk (128,64) | `X_pca` | 0.2429 | ±0.0173 | **−0.0337** |
| trunk (128,64) | `X_scGPT` | 0.2047 | ±0.0011 | +0.0133 |

**The Q1 margin survives the capacity change and moves little:** `X_pca` − `X_scGPT` is **+0.0317**
with a linear head and **+0.0383** with the trunk, both outside the wider of the two seed bands. The
ordering does not depend on the trunk.

**The smaller head is the better one on both arms** — +0.0179 on `X_pca`, +0.0245 on `X_scGPT`.
⚠️ On `X_pca` that gain is **1.03× the seed band** (+0.0179 against ±0.0173), i.e. it clears the band
and no more; on `X_scGPT` it is 12.9× (+0.0245 against ±0.0019). So "the linear head is better" is
load-bearing for scGPT and barely established for PCA.
*(Corrected 13.08.2026: an earlier summary put this ratio at 1.6×. The artifact gives 1.03×, which
makes the PCA-side claim weaker, not stronger.)*

**The training regimes swap with capacity**, which is why "PCA cannot train here" was a statement
about the trunk rather than the representation — median best epoch, from the same cell:

| arch | `X_pca` | `X_scGPT` |
|---|---|---|
| linear | 12 | 2 |
| trunk | 1 | 8 |

**⚠️ Neither representation clears the ridge control on `X_pca`.** `RidgeCV` on line-mean embeddings
scores 0.2767, above both PCA arms. A linear model on averaged embeddings beats the per-cell network
on the representation that wins Q1 — recorded here because it bounds what the per-cell model has
been shown to buy, and it is visible in `panel_alpha_response.png` as the dashed ridge line sitting
above the PCA median at every `alpha`.

### D · The evaluation set moves the answer, not the head count

`4a` §D, `RUN_SECTION_D = True` → `notebooks/outputs/panel/panel_heads_summary.csv`. The same design
with **534 heads** instead of 11, scored twice: on the eleven panel drugs, and on all 534.

| scored on | arch | Q1 margin (`X_pca` − `X_scGPT`) | verdict against the band |
|---|---|---|---|
| 11 panel drugs | linear | **+0.0479** | outside |
| 11 panel drugs | trunk | **+0.0252** | outside |
| all 534 drugs | linear | **+0.0231** | outside |
| all 534 drugs | trunk | **−0.0077** | **inside — a tie, and the sign flips** |

The 534-head model still puts `X_pca` ahead on the same eleven drugs, so **head count does not flip
the ordering**. The same model scored on all 534 halves the margin, and with a trunk tips it into a
tie with the sign reversed. What the panel and the gene-set sweep disagree about is **which drugs
are measured on**, not the model — and an earlier test excluded their response spread as the cause.

**⚠️ The caveat that travels with §D wherever it is quoted.** Across all 534 drugs the model barely
beats the per-drug null: the `vs_null` column runs **+0.00017 to +0.00043**, and for linear/`X_pca`
it is **−0.00005** — worse than the constant. That is the noise floor, and it is why the sweep's
ordering was fragile.
*(Corrected 13.08.2026: an earlier summary gave the lower bound as +0.00002. The artifact's smallest
positive value is +0.00017.)*

### E · The lead is not a small-data effect — the margin grows with labels

`4a` §E, `RUN_SECTION_E = True` → `notebooks/outputs/panel/panel_curve_summary.csv`.
120 fits (four label budgets x two representations x three seeds x five folds), linear head,
`alpha=0`, `loss=mse`, executed 13.08.2026 in 14 min 21 s.

**Only the label supply moves.** The dropped lines' cells stay in the fold, in the batches and in the
per-fold PCA, so both arms keep an identical input at every point on the curve.
`scripts/training/cv.py::oof_predictions(n_label_lines=…)` does the thinning, and it reaches the loss
mask, the density fit and the head-bias init together — leaving any of them on the full set would
feed the dropped lines' labels back in through the side door.

| labelled lines per fold | `X_pca` | `X_scGPT` | Q1 margin | against the wider band |
|---|---|---|---|---|
| 25 | 0.1374 | 0.1338 | **+0.0036** | ±0.0230 — **inside; the arms are indistinguishable** |
| 50 | 0.1947 | 0.1857 | **+0.0090** | ±0.0096 — **inside; still indistinguishable** |
| 75 | 0.2469 | 0.1964 | **+0.0505** | ±0.0181 — outside |
| 103 (all) | 0.2608 | 0.2291 | **+0.0317** | ±0.0072 — outside |

**The small-data reading is refuted.** It was the live alternative: that `X_scGPT` holds more but
needs more supervision before a head can use it, and would overtake given enough labels. That
predicts a margin which *shrinks* as labels are added. The margin instead runs **+0.0036 → +0.0317**,
and at the two smallest budgets the two representations cannot be told apart at all. Nor is the
margin *flat*, which would have placed the advantage purely on the input side — that `X_pca` is
fitted on this atlas and `X_scGPT` is not.

⚠️ **What this does not license.** Four points against bands this wide do not support a fitted slope,
and the notebook computes no trend statistic on purpose. The margin is also not monotone — it peaks
at 75 lines (+0.0505) and is smaller at the full 103 (+0.0317). Read the table as *"indistinguishable
at 25–50, separated at 75–103"*, not as a growth rate.

**Internal check.** §E's full-budget point (0.2608 / 0.2291) reproduces §C's linear arm to four
decimals — the same design reached by a different code path.

**Consequence.** The follow-up this curve was to motivate — **E2**, which shrinks the lines *and*
refits the PCA on only those cells — is now the more interesting of the two, because the input-side
explanation is the one still standing. So is a PCA fitted on a *foreign* dataset.

### The aggregation convention — Q1 survives it, the loss ranking does not

The scoring convention was an open decision (`docs/OPEN_DECISIONS.md` §2): predictions can be scored
**per cell or per cell line**, and folds **pooled into one correlation or averaged over folds**. Those
are two independent binary choices, so four conventions. All four were computed on the **same**
out-of-fold predictions, changing nothing else —
`notebooks/outputs/panel/panel_aggregation_comparison.csv`, 36 arm-seeds, per-cell predictions from
`runs/percell/` (gitignored, written by `4a` §A cell 14).

**Q1 margin (`X_pca` − `X_scGPT`), mean over three seeds:**

| arm | line/pooled | line/per-fold | cell/pooled | cell/per-fold |
|---|---|---|---|---|
| α=0 mae | 0.0214 | 0.0192 | 0.0371 | 0.0298 |
| α=0 mse | 0.0464 | 0.0364 | 0.0536 | 0.0319 |
| α=0.5 mae | 0.0428 | 0.0320 | 0.0608 | 0.0460 |
| α=0.5 mse | 0.0827 | 0.0811 | 0.0939 | 0.0691 |
| α=1 mae | 0.0550 | 0.0397 | 0.0772 | 0.0553 |
| α=1 mse | 0.0648 | 0.0566 | 0.0859 | 0.0644 |

✅ **`X_pca` is ahead in all 24 cells.** The Q1 ordering does not depend on the convention. Its
*size* does — the same arm's margin varies by up to a factor of two across conventions (α=0/mae runs
0.0192 to 0.0371) — so a margin may be quoted only with its convention named, but the **direction is
not at risk**.

⛔ **The loss-arm ranking is a different story: the best arm changes with the convention.**

| convention | best arm on `X_pca` | full ranking |
|---|---|---|
| line / pooled | **α=0.5 mae** (0.2824) | 0.5/mae > 1/mae > 0.5/mse > 1/mse > 0/mae > 0/mse |
| line / per-fold | **α=0.5 mse** (0.2985) | 0.5/mse > 0.5/mae > 1/mae > 1/mse > 0/mae > 0/mse |
| cell / pooled | **α=0.5 mae** (0.2500) | 0.5/mae > 1/mae > 1/mse > 0.5/mse > 0/mae > 0/mse |
| cell / per-fold | **α=1 mse** (0.2572) | 1/mse > 1/mae > 0.5/mae > 0.5/mse > 0/mae > 0/mse |

Three different winners across four defensible conventions. Only one thing is stable: **α=0 is last
under every convention**, on both losses.

**Two systematic effects, worth knowing before any of these numbers is compared to another.**
Per-fold scores are uniformly *higher* than pooled — each fold holds ~30 lines and a correlation over
fewer points runs higher. Cell-level scores are uniformly *lower* than line-level — a line's cells
share one label, so the within-line scatter is noise the correlation cannot use.

> ### ✅ DECIDED 14.08.2026 (Selin) — line-level, pooled
>
> The convention in force is the one already used: aggregate to the cell line, pool the folds. The
> label is per (line, drug), so one point per line is one point per label; cell-level scoring adds
> points with no independent content and weights lines by cell count over a 36× range, and per-fold
> averaging runs on ~30 lines where Spearman is upward-biased — visible above as *every* arm scoring
> higher in that column. **`5_evaluation` §1.8's `order` was already computed this way**, verified arm
> by arm against `panel_aggregation_comparison.csv`, so no number moves.
>
> It fixes the highest-`order` arm (α=0.5/mae). It does **not** settle item 9A, which blocks that arm
> on its guards, and it does not touch Q1.

**Taken with the re-execution instability below, the loss comparison is unresolved twice over:** it
does not survive a re-run, and it does not survive a change of scoring convention. Q1 survives both.
⚠️ The convention half of that is now closed by decision; the re-execution half is being fixed.

⚠️ The MIL arms are absent from this table. Their predictions live in
`notebooks/outputs/mil/mil_oof_predictions.csv`, a separate file this comparison does not read, so
the bag-objective axis is not re-scored under the four conventions.

### Reproducibility — measured over five executions, and it is one fit, not the pipeline

`4a` §A was executed **five independent times** with identical code, seeds and inputs: the run
committed in `02f0fe6`, the full top-to-bottom re-run in `9732b6f`, and three more on 14.08.2026
(~38 min each). Every arm's mean per-drug Spearman was captured each time —
`notebooks/outputs/panel/panel_execution_band.csv`, built by
`scripts/evaluation/aggregation_comparison.py`'s sibling procedure and reproducible from the
committed notebook.

| arm | rep | exec 1 | exec 2 | exec 3 | exec 4 | exec 5 | range |
|---|---|---|---|---|---|---|---|
| MLP mse α=0 | `X_pca` | 0.2541 | 0.2473 | 0.2459 | 0.2490 | 0.2450 | **0.0091** |
| *the other eleven arms* | both | — | — | — | — | — | **0.0000** |

**Eleven of twelve arms are identical to six decimal places across all five runs.** The pipeline is
deterministic. Exactly one arm is not.

**And it is one *fit*, not one arm.** Within that arm, the fold log's `best_epoch` is identical in
all fifteen (fold × seed) combinations across runs, and `best_val_obj` differs in exactly **one**:
fold 1, seed 42 — the arm's first fold, of the first seed, of the first arm trained. **That is the
first fit executed in the process.** 179 of §A's 180 fits reproduce exactly; fit number 1 does not.

**Tested 14.08.2026 — and confirmed, with a complication.** §A was run once with `REPS` reversed, so
`X_scGPT`/α=0/mse became the first fit instead of `X_pca`/α=0/mse:

| arm | rep | normal order | reversed | moved? |
|---|---|---|---|---|
| MLP mse α=0 | `X_scGPT` | 0.2009 | **0.2024** | **yes** |
| MLP mse α=0 | `X_pca` | 0.2473 | 0.2508 | yes |
| the other ten | both | — | — | no |

**`X_scGPT`/α=0/mse had read exactly 0.2009 in all five previous executions** while it sat at
position ~19. Placed first, it moved. That is the confirmation, and it is clean: position causes
instability.

⚠️ **But it is not the only cause.** `X_pca`/α=0/mse moved as well, to 0.2508, **despite no longer
being first** — a value still inside its five-run band. So there appear to be **two** first-use
effects: the first fit of the process, and the first fit that exercises the **PCA path**, which
computes per-fold projections and is always `X_pca`/α=0/mse whatever order `REPS` is in.

A device warm-up addresses the first. Whether it addresses the second is what the verification runs
below measure rather than assume.

**What it costs, stated at its real size.** The instability is **0.0091 on one arm**, and it is not
where the Q1 result lives — §C, §D, §E and the other eleven arms reproduce exactly. But that one arm
is the item-9A **incumbent**, so a 0.0091 wobble decides the loss comparison's verdict
(see below). Any margin under ~0.01 that involves the α=0/mse `X_pca` arm is not resolvable by a
single run; every other comparison in this section is.

⚠️ **Review item 10's finding that "the `mps` nondeterminism does not reproduce under current code"
is wrong, but only just** — it reproduces in one fit of 180, which is why a smaller check missed it.

### What this settles about Q1

The margin moves along four axes, all measured:

| axis | from | to | source |
|---|---|---|---|
| **Objective** | +0.0827 per-cell | **+0.0236 bag** | `panel_metrics.csv`, `alpha=0.5`, `mse` |
| **Capacity** | +0.0383 trunk | +0.0317 linear | `panel_arch_summary.csv` |
| **Scoring set** | +0.0479 on 11 drugs | +0.0231 linear / −0.0077 trunk on 534 | `panel_heads_summary.csv` |
| **Label supply** | +0.0036 at 25 lines (a tie) | +0.0317 at 103 lines | `panel_curve_summary.csv` |

**The objective is the largest of the four axes.** `5_evaluation` §1.8 was executed for the first
time on 13.08.2026 and wrote `notebooks/outputs/panel/panel_metrics.csv`, which scores the per-cell
and MIL arms through one scorer. At `alpha=0.5`, `mse`, over three seeds:

| objective | `X_pca` | `X_scGPT` | Q1 margin |
|---|---|---|---|
| per-cell (`mlp`) | 0.2754 | 0.1927 | **+0.0827** |
| bag (`mil`) | 0.2412 | 0.2177 | **+0.0236** |

**71.5 % of the per-cell margin does not survive the move to a bag objective** — `X_scGPT` gains
(0.1927 → 0.2177) while `X_pca` loses (0.2754 → 0.2412). This is what makes *"PCA beats scGPT"*
dependent on the objective: the statement is about a per-cell loss, not about the representations.

### The item-9A rule, evaluated — no arm wins

The rule (`5_evaluation` §1.3, Selin's: win on `order`, non-inferior on `top_of_order`, `values` and
`spread_slope`, each margin the quantity's own seed band, half-range over three seeds) was applied
for the first time on 13.08.2026. **It was applied twice, four hours apart, and it gave two different
answers.**

| run | incumbent `order` | verdict |
|---|---|---|
| `4a` as committed at `02f0fe6` | 0.2541 | **no challenger wins** — all thirteen blocked |
| `4a` re-executed top to bottom, 23:30 | **0.2473** | **`alpha=0` / `mae` / `X_pca` WINS** |

Nothing changed between them but the run. Same code, same seeds, same inputs; the incumbent arm's
`order` fell by **0.0068** and that was enough to let a challenger through. Every other arm in the
leaderboard is identical across both.

**Five executions later, the incumbent's range is 0.2450–0.2541** (§*Reproducibility*), so the
verdict depends on which run you ask: at 0.2541 no challenger clears the bar, and at 0.2450–0.2490
`alpha=0` / `mae` does. Four of the five runs put it below the flip point, but "four of five" is not
a basis on which to declare a winner.

⛔ **What this means for item 9A, stated plainly: the sweep's answer is not stable, and the
instability is the same size as the effect.** The gap it turns on (incumbent vs `mae` at
`alpha=0`) is +0.0076 in one run and +0.0144 in the other. A decision rule cannot resolve a
difference smaller than the run-to-run variation of its own inputs, and this one is not.

⚠️ **This refutes review item 10's finding that "the `mps` nondeterminism does not reproduce under
current code."** It reproduces. See §*Reproducibility* below.

### Counter-evidence — the external benchmark does not separate the arms

`notebooks/analysis/evaluation/dreval_benchmark.ipynb` →
`notebooks/outputs/dreval/dreval_lco_results.csv`, five leave-cell-line-out folds. Normalised
Spearman for the two OncoMLP arms:

Re-executed 13.08.2026 against the current artifacts (1 min 57 s, clean):

| fold | `X_pca` | `X_scGPT` | ahead |
|---|---|---|---|
| 1 | **0.2331** ⚠️ unstable | 0.2439 | scGPT |
| 2 | 0.2764 | 0.2903 | scGPT |
| 3 | **0.3576** | 0.2657 | PCA |
| 4 | 0.2625 | 0.2638 | scGPT, by 0.0013 |
| 5 | 0.2581 | 0.2962 | scGPT |

The fold-to-fold spread on `X_pca` alone (0.2331–0.3576) is larger than any between-arm gap in the
table, and **fold 1 is not stable across runs at all** (0.2025 → 0.2776 → 0.2331 over three
executions; folds 2–5 are byte-identical every time). Under this protocol the two representations are
**not distinguished**, and the fold that most affects the tally is the one that will not sit still.

> ⛔ **Correction 13.08.2026 — a fold-1 value had been reported as the result.** The Gate 5 execution
> log recorded *"scGPT edges PCA (0.747 vs 0.738 raw; 0.244 vs 0.203 normalised)"*. Those figures are
> **fold 1 only**, taken from a 50-row file of five folds × ten algorithms, and the direction they
> assert did not survive the other four folds even then. On the re-run, fold 1 itself reverses
> (`X_pca` 0.2025 → 0.2776, the only fold that moved), so the sentence is now wrong about its own
> fold as well. Do not repeat it.

**Ruled out by measurement, not by argument:** the metric, drug-level response spread, a misuse of
scGPT (the call matches `Tutorial_Reference_Mapping` from the scGPT repo), and broken embeddings.
CPM as scGPT's input is correct because its per-cell binning takes quantiles of the cell's own
non-zero values ([Step 02](02-preprocessing-and-embeddings.md), §*What scGPT is fed*).

**Consequence for how Q1 is stated.** *"PCA beats scGPT"* is not defensible without naming the
objective and the scoring set: on the full 534-drug catalogue with a trunk the margin is a tie with
the sign reversed, and near the per-drug noise floor throughout.

---

## Q2 on the rebuilt panel — does a per-cell model learn heterogeneity implicitly? (14.08.2026)

> ✅ **Not covered by the page banner above.** Measured after the pipeline review, on the rebuilt
> 11-drug panel and the `auc_cc` target. **This is the first entry Q2 has had in the scientific
> record** — until 14.08.2026 its results existed only in `notebooks/4b_mil_training.ipynb` and
> `notebooks/outputs/mil/`, with nothing in `docs/steps/` or the report.

**What is being asked, and why it is hard.** There are no per-cell labels: one bulk response value
per (cell line × drug) is broadcast onto ~300 cells. So any within-line structure in the predictions
is something the model **imposed**, not something it was shown. `4b` therefore does not ask "is the
model right about individual cells" — it cannot — but "does the model produce within-line structure
that is more than noise, and is that structure anything other than technical artefact".

**The instrument.** Seven staged tests, each a null, a comparison or a collapse test, run over three
seeds on both representations. Source: `notebooks/4b_mil_training.ipynb` →
`notebooks/outputs/mil/q2_verdict.csv` and `stage0`–`stage7`.

| stage | question | `X_pca` | `X_scGPT` |
|---|---|---|---|
| 0 | does the **input** carry within-line variation at all? | 0.416 share, no collapse | 0.461, no collapse |
| 1 | does MIL predict **more** within-line spread than the per-cell model? | 100 % of (drug, line) pairs, all 3 seeds | 82.3 %, all 3 seeds |
| 2 | does that structure **reproduce across seeds**? | median ρ 0.256, 83.0 % beat a shuffled-cell null | ρ 0.861, 99.2 % |
| 6 | do **technical confounds** explain it? | R²adj 0.0533 — **82 % as much as the signal reproduces** | 0.2656 — 36 % |
| 7 | can the instrument detect a **known** gap? | AUROC **0.521**, 46.0 % of pairs significant | **0.537**, 48.4 % |

**The verdict, in the notebook's own words: `Q2(a) POSITIVE` for both representations** — but stated
there with the qualifier that matters, *"at a measured instrument sensitivity of AUROC 0.521"*.

### Why that positive is weak, and stage 7 is the reason

Stage 7 is a **positive control**: it asks whether the instrument can see a difference it already
knows is there. The answer is *barely* — 0.521 and 0.537 against a 0.5 null, with **under half** the
pairs reaching significance. A positive from stages 1 and 2 read through an instrument that
insensitive does not support a claim about biology.

For `X_pca` the confound result compounds it: sequencing depth, genes detected, mitochondrial
fraction and cell-cycle scores together explain **82 % as much variance as the signal reproduces
across seeds**. The structure is real in the sense that it is not noise; it is not established to be
*biological*.

⚠️ **Q2(b) and Q2(c) are not addressed and cannot be with these measurements** — that is `4b` §1's
own scoping, not a caveat added here. Whether the within-line variation is real **drug-response**
heterogeneity (b), and whether it predicts **which cells survive** (c), both need post-treatment
single-cell data, which this project does not have
([Step 01](01-datasets-and-harmonization.md#post-treatment-single-cell-data--what-would-be-needed-for-q2-b-and-what-exists)).

### What is verified, and what is not

**Verified 14.08.2026:** all eleven `outputs/mil/` artifacts match `HEAD` unmodified; `q2_verdict.csv`
re-derives from its own stage tables under the notebook's aggregation (median across seeds, and the
median of per-seed medians for stage 7); three seeds are present in every stage table; and stage 1's
`median_sd_4a` matches the committed `panel_within_line_spread.csv` at α=0.5/mse exactly for all six
representation × seed combinations.

⚠️ **Not verified by execution.** `4b` has not been *run* since 13.08.2026 09:47. The above is a
consistency check against current artifacts, not a fresh run.

⛔ **Q2 has no figure.** `notebooks/outputs/mil/` holds eleven CSVs and no images. Every Q2 claim in
this project is a table or a number.

---

## Multi-task masked loss over all 545 CTRPv2 drugs (26.05.2026)

The target artifacts (`obsm["Y_ctrp"]`, `obsm["M_ctrp"]`, `uns["ctrp_drugs"]`, the legacy flat
`viability_<drug>` columns) are defined once in [Step 03](03-model-and-training-design.md). What this
run adds is the split and the drug scope:

- `obs["split_ctrp"]` — **one drug-agnostic, cell-line-grouped 70/15/15 split** written by
  `create_splits.py` `run_multi()`, shared across all heads (leakage-free for every drug at once;
  a single shared split is only possible *because* the leakage control is at the cell-line level).
- **Drug-scope filter:** keep a drug only if screened on ≥ `--min-cell-lines` overlapping cell lines
  (default 50). This run used **`--all-drugs` (= min 0) → K = 545 drugs**.

**Run-time overlap reported by the pipeline:** **180 / 198** SCP542 cell lines overlap
CTRPv2 (180 = lines with actual post-QC measurements; the audit's 190 counts roster name-matches — see
[Step 01](01-datasets-and-harmonization.md)).

> ⚠️ **10.08.2026 — this becomes 181 at the next sweep, and every number on this page with it.** The
> name join was dropping `NCIH292`, which CTRPv2 spells `H292`, and was double-counting experiments
> listed once per calendar day. Both are fixed in `ctrp_to_h5ad.py`; no artifact reflects them yet, so
> everything below still describes the 180-line matrix. Evidence and effect sizes:
> [Step 01](01-datasets-and-harmonization.md#the-join-dropped-a-screened-cell-line-h292-10082026).

**`split_ctrp` distribution (one cell-line-grouped **70/15/15** split, shared by all heads):**

| split | lines | % of lines | cells | % of measured cells |
|---|---|---|---|---|
| train | 126 | 70.0% | 34,126 | 72.3% |
| val   | 27  | 15.0% | 7,121  | 15.1% |
| test  | 27  | 15.0% | 5,980  | 12.7% |

- 70/15/15 is the design target at the **cell-line** level (`create_splits._split_cell_lines`); the
  **cell** percentages differ slightly because lines carry different cell counts.
- `unassigned` = **18 lines / 6,286 cells** (SCP542 lines with no CTRP measurement; 198 → 180 measured).
- **Cross-validation** (`notebooks/4a_percell_training.ipynb` §B2) **holds `test` out** and resamples only the
  153 train+val lines via 5-fold GroupKFold, so test is never seen in CV.

**Model & training:** a single `OncoMLP` with `output_dim = K`, fed by `MultiDrugDataset`
(`scripts/model/dataset.py`) whose 3-tuple `(x, y, mask)` batches `train_model` auto-detects to
switch into **masked MSE** (mean over observed entries only). Up front,
`train_multitask._per_drug_constant_mse` computes a **per-drug-mean sanity baseline** — the proper
null model here: for each drug it predicts the constant train-set mean viability over that drug's
observed cells. Because labels cluster near 1.0, that constant is already a strong predictor, so a
head only counts as having *learned* response if it **beats its own drug's constant**.

**Shared hyperparameters** (from `config.json` / `run_meta.json`): batch 128, epochs 50
(early-stopped), lr 1e-3, weight_decay 1e-3, dropout 0.5, input_dropout 0.1, grad_clip 1.0,
scheduler patience 3, early-stop patience 10, seed 42, loss MSE, norm LayerNorm.
**Matched trunk + matched width.** Both reps use the **same** hidden layers `(128,64)`
(`DEFAULT_HIDDEN_DIMS`, set 14.06.2026) **and** the same **512-d** input (`X_pca` raised from scanpy's
~50 default to `add_pca.DEFAULT_N_COMPS = 512` on 27.06.2026), so the entire network — including the
first projection's parameter count — is identical and **only the representation differs**. This closes
the last comparison confound. (History: the original matrix used a `(64,32)` PCA trunk and a ~50-d PCA,
both of which handicapped PCA; the numbers below supersede those.)

### Metrics — what each number means

Every result below is one of these. They are reported on the **val** split (single fixed split) or as
**5-fold CV mean ± std** (test held out); read them together — MSE alone is misleading near a viability
of 1.0.

- **Masked (val) MSE** — the training objective. Per-cell squared error `(pred − viability)²`, averaged
  **only over observed `(cell × drug)` entries** (`mask = 1`); missing labels contribute nothing
  (`_masked_mean` in `training_utils.py`). For a single drug it's plain MSE. **Why ≈ 0.01 is
  misleading:** viability is per (cell line × drug), broadcast to all the line's cells, and clusters
  near 1.0 with tiny variance — so even predicting a constant scores ~0.01. Absolute MSE therefore says
  little; what matters is whether it beats the constant and whether it *ranks* lines (below). *Train*
  MSE is logged with dropout **active**, so it can sit below or above the dropout-free val MSE.

- **Per-drug-mean baseline** — the null model (`_per_drug_constant_mse`). For each drug it predicts that
  drug's **train-set mean viability** for every cell, then is scored on val. Because labels are near
  constant this is already a *strong* predictor, so it's the bar every head must clear; a head only
  counts as having *learned* response if it beats its own drug's constant.

- **Heads beating baseline (`heads_beat`)** — the **count** of the K = 545 drugs whose model val MSE is
  below their per-drug-mean baseline. Intuitive, but a **thresholded count of near-ties**: most heads
  have model ≈ constant (labels ≈ 1.0), so they sit on the decision boundary, and the per-fold baseline
  is recomputed from that fold's train lines. If a fold's held-out lines are collectively a little
  above/below the train mean, **hundreds of heads flip together** (common-mode), so the CV std is huge
  (±73–94; cf. √(K·p(1−p)) ≈ 11 if heads were independent — observed is ~8× that). **Treat as
  directional, not precise.**

- **Δmse (model − baseline)** — the **continuous** counterpart of heads-beating: the mean over drugs of
  `model_mse − baseline_mse`. **Negative ⇒ model better** than the constant on average; it is not
  thresholded, so it doesn't suffer the count's instability (its CV std is small relative to the mean).
  Reported as CV mean ± std; the per-fold `archive/training_545_mean_pv/cv_folds.csv` also carries `median_delta` and `frac_beat`
  (= `heads_beat / n_total`).

- **Overfitting gap** — `val_mse − train_mse` at the best epoch (single-task). Larger ⇒ more
  memorization; the core hypothesis predicts scGPT < PCA. Same dropout-in-train caveat as above, so it
  is indicative, not exact.

- **Per-drug correlation (Spearman / Pearson)** — the metric that actually asks *does the model rank
  cell lines?* For each drug, predictions are averaged to one value **per held-out cell line** and
  correlated with the true per-line viability across lines (Spearman = rank, Pearson = linear).
  Restricted to drugs with **real response variance** (per-line true std ≥ 0.05) and ≥ 5 val lines —
  otherwise there is nothing to rank. Insensitive to the near-1.0 offset that dominates MSE.

- **5-fold GroupKFold CV (test held out)** — robustness wrapper: `GroupKFold(5)` over `Cell_line`
  resamples the 153 train+val lines into 5 train/val folds (no line on both sides), each retrained from
  scratch; we report **mean ± std**. The fixed `test` set is never touched, so **CV numbers are a
  stability check, not a test-set estimate**.

**The 8-run matrix (512-d, 27.06.2026; all share `split_ctrp`, n_train 34,126 / n_val 7,121).**
Per-drug-mean baseline: **~0.043** (K=1 paclitaxel, data-derived, rep-independent), **0.0097** (K=545).
Reproducible in `notebooks/4a_percell_training.ipynb (§B)`; run dirs `runs/20260627_1913xx_*` (see
`runs/runs_index.csv`).

**Single-task (K=1 paclitaxel) — the overfitting story** (gap = val − train, at the best epoch):

| Gene set | Rep | Train MSE | Val MSE | Gap (val−train) |
|---|---|---|---|---|
| `hvg5000` | scGPT | 0.037 | 0.041 | **0.004** |
| `hvg5000` | PCA | 0.011 | 0.045 | 0.033 |
| `all_genes` | scGPT | 0.032 | 0.045 | 0.013 |
| `all_genes` | PCA | 0.042 | 0.039 | −0.003 |

**All-drugs (K=545) — heads beating the per-drug-mean baseline:**

| Gene set | Rep | Val MSE | Heads beat baseline |
|---|---|---|---|
| `hvg5000` | scGPT | 0.0105 | 147 / 545 |
| `hvg5000` | PCA | 0.0103 | **169 / 545** |
| `all_genes` | scGPT | 0.0106 | 131 / 545 |
| `all_genes` | PCA | 0.0106 | **138 / 545** |

> ⛔ **05.08.2026 — the `all_genes` rows are not a full-transcriptome comparison.** At `max_length=1200`
> every cell in `all_genes` exceeds the cap, so scGPT received a random fraction of each cell's expressed
> genes while PCA received the whole gene set — counts in
> [Step 02](02-preprocessing-and-embeddings.md#why-hvg-5000-is-the-default-03082026). The two arms of an
> `all_genes` row therefore differ in *gene set* as well as in encoding, so **no PCA-vs-scGPT contrast
> may be drawn within those rows**, and none of them supports a statement about scGPT and the full
> transcriptome — scGPT never received it. The `hvg5000` rows are unaffected; the cap binds in a single
> cell there. The decision to keep 1,200 is
> [here](02-preprocessing-and-embeddings.md#decision--one-seeded-draw-at-1200-all_genes-is-a-sanity-check-03082026).
> These embeddings were additionally generated **unseeded** — and that part is *not* confined to
> `all_genes` (widened 10.08.2026, review item 4). Two things in the embedding path draw on the RNG: the
> gene subsample above, and the tie-breaking inside scGPT's value binning, where `_digitize`
> (`scgpt/preprocess.py:239`) resolves values landing on a repeated bin edge with `np.random.rand`. The
> second touches **every cell that has tied expression values, in both variants**. So no embedding on
> disk is exactly reproducible, `hvg5000` included; only the *truncation* caveat above is
> `all_genes`-specific. The seed fix — `np.random.seed(42)` beside `torch.manual_seed(42)`,
> `gen_embeds.py:243-250` — postdates every embedding on disk and takes effect at the sweep.

**Reading the results (matched trunk + matched 512-d width):**

- **Core hypothesis — supported (single-task, `hvg5000`):** scGPT's train/val gap is **0.004** vs
  PCA's **0.033** — scGPT overfits far less. Matching PCA to 512-d *sharpened* this: PCA's extra
  first-layer capacity lets it fit the train set harder (train 0.011) while val stays high (0.045),
  exactly the memorization the denoised scGPT prior is meant to avoid.
- **All-drugs — PCA competitive/better on raw accuracy:** heads-beating `hvg5000` **PCA 169 vs scGPT
  147**, ~~`all_genes` **PCA 138 vs scGPT 131**~~ (struck 05.08.2026 — not like-for-like, see the block
  above); val MSEs are within 0.0003. scGPT does **not** win on absolute predictive metrics.
- **Net:** scGPT's robust, reproducible win is **lower overfitting**, not higher accuracy — and this
  now holds with input dimensionality matched, so it can no longer be dismissed as a capacity artifact.
- **Which heads are even learnable** is driven by coverage + response variance — see
  `notebooks/analysis/harmonization/drug_coverage.ipynb`: the ≈16-line drugs (n_val 221) are the unreliable/hardest heads,
  while high-coverage high-variance drugs (docetaxel, gemcitabine, oligomycin a) are the easiest.

> ⚠️ **Gap-metric caveat.** Train MSE is logged with dropout (0.5) + input-dropout (0.1) **active**, so
> it can sit *below or above* the (dropout-free) masked val MSE; the gap is indicative, not exact. The
> `all_genes` rows early-stop very fast (best epoch 1–4), so their gaps are noisy — `all_genes`·PCA's
> **−0.003** reflects near-no learning + the dropout offset, not genuine negative generalization. The
> clean comparison is `hvg5000` single-task (scGPT 0.004 vs PCA 0.033).

### Is the difference real? — 5-fold cross-validation (27.06.2026)

The single-split numbers above rest on **27 val lines**, so they are point estimates. To test
robustness, `cv_evaluate` (`notebooks/4a_percell_training.ipynb` §B2) runs **5-fold GroupKFold over `Cell_line`,
holding the fixed `test` set out** and resampling only the 153 train+val lines (~122 train / ~31 val
per fold). On `hvg5000`:

| Rep | Heads beating baseline (mean ± std) | Δmse model−baseline (mean ± std) | All-drugs val MSE | Paclitaxel gap (val − train) |
|---|---|---|---|---|
| `X_pca` | **207 ± 73** / 545 | **+0.00058 ± 0.00040** | 0.0106 ± 0.0008 | **+0.011 ± 0.020** |
| `X_scGPT` | **191 ± 94** / 545 | **+0.00072 ± 0.00047** | 0.0107 ± 0.0009 | **−0.002 ± 0.014** |

- **The continuous metric is the honest one — and it's negative news:** Δmse is **positive for both
  reps** (4 of 5 folds), i.e. on average the model is **marginally *worse* than the per-drug-mean
  constant**. The heads-beating count (~190–207 of 545, i.e. < 40% of heads) said the same thing all
  along; the continuous Δ just makes it unambiguous and stable (std ≪ the count's). PCA's Δ (+0.00058)
  is slightly *less bad* than scGPT's (+0.00072) — same direction as heads-beating.
- **The heads-beating count itself is *not* robust:** the fold std (±73–94) **dwarfs** the PCA−scGPT
  difference (~16). The single-split "169 vs 147" is within fold noise — don't read it as a real PCA
  advantage. (See *Metrics* above for why the count swings so hard.)
- **The overfitting direction survives, weakly:** mean paclitaxel gap is lower for scGPT (−0.002) than
  PCA (+0.011), consistent with the denoised-prior claim, but the spreads overlap.

### Better metric — per-drug correlation (27.06.2026)

Because viability clusters near 1.0, beating the per-drug-mean on MSE is a weak bar. §3 of the notebook
instead correlates **predicted vs true viability across held-out cell lines**, per drug (Spearman +
Pearson), restricted to the 461 drugs with real per-line variance (std ≥ 0.05, ≥ 5 val lines):

| Rep | mean Spearman | median Spearman | frac. drugs ρ > 0.3 |
|---|---|---|---|
| `X_pca` | −0.02 | −0.01 | 4.3% |
| `X_scGPT` | −0.05 | −0.05 | 3.9% |

- **Sobering:** per-drug rank correlation is **≈ 0 for both reps** — the models do **not** rank cell
  lines by drug response. The marginal MSE "wins" over the per-drug-mean reflect shrinking toward the
  constant, **not** real per-line predictive power. At this resolution (per-line viability broadcast to
  cells, values ≈ 1.0) the task is barely learnable beyond the mean — for *either* representation.
- This reframes the whole comparison: the scGPT-vs-PCA question is secondary to the fact that **neither
  rep yet predicts response variation across lines**. Motivates the better-target / better-metric work
  in [TODO.md](../TODO.md) (correlation-based selection, drugs with real variance).

> ⚠️ **Superseded (13.07.2026), and partly an artifact.** This verdict does not survive: it averages
> over 545 drugs, and the multi-task loss it was measured under was unstandardized.
> `notebooks/archive/target_comparison.ipynb` reproduces the failure on demand. Full account, including the
> decomposition of what actually produced the later gain:
> [Corrections](corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result).

### Learnability-filtered subset — the signal was there all along (13.07.2026)

`notebooks/archive/learnability_filter.ipynb` → `notebooks/archive/learnable_subset_training.ipynb`. The §3 null result
above pooled a few learnable heads with hundreds of flat, inert ones. **Filter first, then ask.**

**The filter (`learnability_filter`).** The learnability score of [`drug_coverage`](../../notebooks/analysis/harmonization/drug_coverage.ipynb)
(`resp_std × cov_frac`) is **degenerate on `auc_z`** — the target is z-scored per drug, so every drug
has std exactly 1.0 and all 545 tie. Spread is therefore measured on the **raw `auc` scale**, recovered
exactly via `uns["ctrp_score_scale"]`/`["ctrp_score_center"]` (and that `scale` vector *is* the per-drug
std of `auc`). The loose `drug_coverage` gates (`cov ≥ 100 & std ≥ 0.05`) kept 439/545 and so never bit; the
missing condition is **differential response** — a drug must both **kill** a real population of lines
(`n_sens`: `auc ≤ 0.5`) and **leave one alive** (`n_res`: `auc ≥ 0.8`). A uniformly inert or uniformly
toxic drug has no cross-line ranking to learn, however well covered it is. **6 / 545 pass; the top 5 by
learnability are trained.**

**What the raw label distribution looks like, and why a spread filter alone cannot bite.**
`notebooks/analysis/harmonization/drug_coverage.ipynb` → `outputs/data/target_distribution.png`, four panels: (A) the
viability histogram, (B) the per-drug response-std histogram, (C) the coverage-vs-std filter scatter,
(D) per-drug response bands.

| Quantity | Value |
|---|---|
| viability across all (line × drug) pairs | clusters near 1.0 — median **0.91**, 75 % ≥ 0.8, bands squeezed into ~0.8–1.0 |
| per-drug response std | median **0.088**; only **3 %** of drugs are truly flat |
| what `cov ≥ 100 & std ≥ 0.05` keeps | **439 / 545** |

So coverage and spread together remove barely a fifth of the catalog — which is exactly why the gate
above needed a differential-response condition, and why the near-1.0 label distribution makes absolute
MSE uninformative (see *Metrics* above).

**Per-drug coverage** (`analysis/harmonization/drug_coverage.ipynb` §2–§3, printed output): **no drug
covers all 180 lines** — max 179, median 171 — **382 drugs clear 90 %** coverage, 80 drugs fall below
50 %, and 14 have std < 0.05.

> ⚠️ **Dated marker, 13.08.2026 (Selin).** Two further claims stood in this paragraph and have been
> removed. Recorded rather than silently dropped, because both were stated here as fact.
>
> *"The ~16-line drugs (`n_val` 221) are the unreliable heads that dominate the worse-than-baseline
> lists"* — **retracted, not awaiting refresh.** Both quantities came from `drug_coverage` §4–§5, which
> merged the newest `runs/*_all_drugs/per_drug_results.csv` **by modification time**. Of the 17 such runs,
> 16 are on the void `auc` target and one on the retired `auc_z`, and the one that selection picked ran
> **one epoch in four seconds** (`runs/20260713_141229_multitask_X_pca_all_drugs/run_meta.json`,
> `05:12:29`→`05:12:33`). R4 re-runs §B of `4a_percell_training` and will produce a K=545 run on `auc_cc`,
> but the cells that turned such a run into this claim are archived, and archived notebooks are not
> re-run — so **no rerun regenerates it**. Restoring the claim means restoring those cells.
>
> *"Per-drug values are in `outputs/*_drug_learnability.csv`"* — **removed permanently.** The
> learnability criterion is [retracted](./corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability)
> and its tables sit in `outputs/archive/`, which is by definition what a standard run cannot recreate.
>
> The coverage figures above are **unaffected** — they come from §2–§3, which read no model output. They
> move 180 → 181 lines at R5. The parenthetical formerly credited the drug-coverage figure with these
> numbers; that figure contains none of them (both its panels are model output), so the attribution has
> been corrected to the notebook sections that print them. The figure itself now sits at
> `outputs/archive/drug_coverage.png`, having been under `data/` until its producing cells were dropped
> on the same day.

**The result (`learnable_subset_training`).** Both reps trained on those 5 heads (matched trunk, on the then-current `auc_z`); the honest
metric is per-drug Spearman on **cross-validated out-of-fold predictions** — 5-fold GroupKFold over the
153 train+val lines, so every line is ranked by a model that never saw it (~150 lines per drug, versus
the 27 the fixed val split would allow):

| Rep | mean Spearman | mean Pearson | heads beating baseline | best val MSE |
|---|---|---|---|---|
| `X_pca` | **0.432** | 0.416 | 3 / 5 | 0.925 |
| `X_scGPT` | **0.488** | 0.482 | 4 / 5 | 0.777 |

(On `auc_z` the per-drug-mean null model scores MSE = **1.0** by construction, so these MSEs are
readable directly.) Per drug: `ml162` 0.59/**0.65**, `1s,3r-rsl-3` 0.58/**0.59**, `dasatinib`
0.52/**0.56**, `cay10618` **0.36**/0.35, `kx2-391` 0.11/**0.28** (PCA / scGPT).

- **Signal exists.** Against −0.02 / −0.05 over 545 drugs, the same architecture reaches ~0.45 here.
  **The 545-drug null result was a drug-selection artifact, not a representation failure.** The
  standing conclusion "the ceiling is the label, no gene representation can help" is *true on average
  and false on the drugs that matter* — **drug selection is a first-class lever**, and a cheap one.
- **The biology checks out.** The two strongest drugs are the **GPX4 inhibitors** (`ml162`,
  `1s,3r-rsl-3`): ferroptosis sensitivity tracks a cell's lipid-peroxidation/redox state, which is a
  *transcriptional* state. `dasatinib` (SRC/ABL) follows target addiction. The filter selected drugs
  whose variance has a transcriptional cause, not merely high-variance drugs.
- **First non-tie between the reps — and it survives a seed check (13.07.2026).** scGPT leads on every
  aggregate and on 4/5 drugs, most clearly where PCA collapses (`kx2-391`, 0.28 vs 0.11). Repeating the
  **K=545 `auc_z`** configuration over **3 seeds** (`target_comparison`, `outputs/archive/target/seed_stability.csv`):

  | seed | PCA | scGPT | gap |
  |---|---|---|---|
  | 42 | 0.388 | 0.430 | +0.043 |
  | 1 | 0.367 | 0.434 | +0.066 |
  | 7 | 0.355 | 0.472 | +0.117 |

  **Gap = +0.075 ± 0.038, sign-consistent across all three seeds.** No longer a one-seed accident — but
  3 seeds × 5 evaluation drugs is **consistent evidence, not a proven margin**. Do not upgrade it to a
  headline claim without more seeds and a wider drug set.
- **Ranking ≫ calibration.** `pred_std` is 0.53 (PCA) / 0.47 (scGPT) against a true spread of 1.0 — both
  models hedge toward each drug's mean. This is **not** an over-regularization artifact: `pred_std ≈ ρ ×
  true_std` is exactly what an MSE-optimal predictor must do (see the ablations below). Fine for ranking;
  to report in AUC units, divide by ρ.

**Is it the model? No — four knobs, all flat (`notebooks/archive/ablations_and_rescue.ipynb`, 13.07.2026).**
Regularization (none → heavy), capacity (74,629 → 2,565 params), batch size (32/128/512) and sample
reweighting (line-balanced, focus-extremes) all leave out-of-fold Spearman within noise of the defaults
(PCA 0.41–0.44, scGPT 0.44–0.49). With regularization *off*, PCA memorizes the training lines (train MSE
≈ 0.01) and still reaches only 0.42 out-of-fold — the model is not being suppressed, it is out of signal.
Full table and reasoning in [Step 03](03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026).

**And the baseline that actually binds:** `RidgeCV` on the **150 cell-line mean embeddings** — no single
cells, no network — scores **0.428**, *tying* the PCA MLP (0.428) and within 0.06 of the scGPT MLP
(0.487). The whole deep single-cell apparatus currently buys **+0.06 Spearman, and only for scGPT**
(whose linear head drops to 0.438 — it *needs* the hidden layer; PCA does not). **Ridge on line means is
the baseline to beat from now on.** The cause is structural: the label is per cell line, so there are
~150 independent examples and the 34k cells are an illusion of sample size — which is why the remaining
levers are **label-side** (more lines, bulk pretraining, denoising), not model-side.

### The learnability gate measured the wrong quantity (27.07.2026)

The gate defined above — a drug must **kill** (`auc ≤ 0.5`) and **spare** (`auc ≥ 0.8`) a real population
— is not the right criterion for this project, and the 10-drug panel inherits the error.

**The mismatch.** `auc ≤ 0.5` asks *does the line die?*, i.e. it filters on absolute potency, which is
essentially the per-drug mean. But the target is `auc_z`, which **subtracts that mean**
([Corrections](corrections-and-dead-ends.md#auc_z-as-the-training-target)),
and the evaluation metric is **Spearman**, which only reads the *ordering* of lines. Whether a drug's
values sit around 0.4 or around 0.9 is irrelevant to both. The gate optimizes for a property the model
is neither given nor scored on.

**What it costs — `nutlin-3` as the clean example.** Raw-scale spread across the panel:

| drug | `auc_mean` | `auc_std` | kill (`≤0.5`) | gate verdict |
|---|---|---|---|---|
| `dasatinib` | 0.631 | **0.155** | 35 | selected |
| `nutlin-3` | 0.874 | **0.147** | **0** | rejected |

**Nutlin-3 has essentially the same spread as dasatinib** — its lines differ just as much, the whole
distribution simply sits higher (~0.6–1.0 instead of ~0.3–0.9). The reason is pharmacological, not
technical: nutlin-3 is **cytostatic**, not cytotoxic. p53 activation drives arrest and senescence, so
viability never falls below 50 % however sensitive the line is. Any threshold phrased as "does it kill"
is structurally blind to every cytostatic agent — which is a large fraction of targeted therapy.

This is not one unlucky drug: **116 of the 545** have zero kills yet `auc_std ≥ 0.10` and coverage
≥ 90 % (`oxaliplatin` among them). The gate discarded all of them silently.

**Why it matters beyond drug choice.** `nutlin-3`/TP53 is the single strongest association in the GDSC
pharmacogenomic screen, and it is *expression*-readable, not only genomic: `MDM2`, `CDKN1A` and
`RPL22L1` — p53 target genes — are selected in ~90–100 % of published gene sets predicting nutlin-3a
sensitivity. It is close to a best case for this model, and the filter threw it out.

**Correction to adopt:** replace kill/spare with **spread on the raw AUC scale** (`auc_std`, recoverable
exactly via `uns["ctrp_score_scale"]`, which *is* the per-drug std) plus coverage. This is the same
quantity that governs the noise-amplification problem in Step 03, so one criterion fixes both: high
`auc_std` = real signal to rank *and* a safe denominator for the z-score. Not yet re-run — the 10-drug
results above and the 8-drug panel below both still rest on the old gate.

### Published sensitivity determinants for CTRPv2 compounds

*Evidence base, not a panel.* These are compounds whose cell-line sensitivity has an **independently
published determinant**, collected 25.07.2026 to select drugs by citation rather than by our own label
statistics. The scientific rationale, which the rebuild keeps: where sensitivity is a documented,
mechanistically understood function of cell state, a transcriptome-based model *ought* to work — so a
failure there is a model result rather than a label artifact.

> ⛔ **The 8-drug panel assembled from this table is VOID** — the candidate list it was drawn from had
> been pre-filtered on our own response values. The panel decision, why it failed, and the 32 approved or
> clinical compounds it wrongly excluded are in
> [Corrections](corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it). **The citations below are
> unaffected.**
>
> ⚠️ **The rebuild happened on 12.08.2026 and did not use the criterion this passage promised.** It said
> the same criterion would be re-applied to "a pool built on coverage and `auc_std` only"; `auc_std` was
> rejected, because spread is still our own label statistic and selecting on it keeps the selection
> label-dependent. The
> [panel](01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026)
> is selected on FDA approval and published determinants instead, and of the eight compounds below only
> `paclitaxel`, `dasatinib` and `afatinib` are in it.

Coverage (`cov`) is the fraction of the 180 trainable lines. `kill`/`spare` counts are shown only because
they document how the void panel was ranked — **the rebuild must not use them**
([why](corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability)).

| drug | target | kill / spare | cov | published determinant |
|---|---|---|---|---|
| `methotrexate` | DHFR | 52 / 31 | 0.94 | **`SLC19A1`** (reduced folate carrier) governs uptake; its loss is a classical resistance mechanism in cell lines — [Zhao & Goldman 2014](https://pubmed.ncbi.nlm.nih.gov/24396145/), [Wright et al., *Nature* 2022](https://www.nature.com/articles/s41586-022-05168-0) |
| `dasatinib` | SRC/ABL, EPHA2, KIT | 35 / 27 | 0.98 | six-gene **expression** model predicts sensitivity in 92 % of held-out breast and 83 % of lung lines — [Huang et al., *Cancer Res* 2007](https://aacrjournals.org/cancerres/article/67/5/2226/534297/Identification-of-Candidate-Molecular-Markers); `LYN` in lung ADC — [*Oncotarget* 2016](https://www.oncotarget.com/article/12657/text/) |
| `paclitaxel` | tubulin | 66 / 25 | 0.94 | **`ABCB1`** efflux + **`TUBB3`** — [*Oncotarget* 2016](https://www.oncotarget.com/article/9118/text/), [*Br J Cancer* 2016](https://www.nature.com/articles/bjc2016203) |
| `vincristine` | tubulin | 61 / 19 | 0.99 | same `ABCB1`/`TUBB3` axis (shared microtubule-disruptor resistance mechanism) |
| `afatinib` | EGFR, ERBB2 | 19 / 26 | 0.91 | `EGFR`+`ERBB2` co-amplification — [*Cancer Discov* 2019](https://aacrjournals.org/cancerdiscovery/article/9/2/199/10771/EGFR-and-MET-Amplifications-Determine-Response-to). ⚠ receptor *expression* alone did **not** correlate in pancreatic lines — [*Br J Cancer* 2011](https://www.nature.com/articles/bjc2011396) |
| `topotecan` | TOP1 | 37 / 18 | 0.97 | **`SLFN11`** expression, the canonical topoisomerase-inhibitor marker — Zoppoli et al., *PNAS* 2012 (NCI-60 + CCLE), pan-cancer replication [*PLOS One* 2019](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0224267), [review](https://www.sciencedirect.com/science/article/pii/S1359644625002922) |
| `tanespimycin` (17-AAG) | HSP90 | 14 / 44 | 0.96 | **`NQO1`** expression bioactivates the benzoquinone to its potent hydroquinone form; correlation confirmed in CCLE *and* GDSC across 7 cancer types — [*PLOS One* 2016](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0153181), [*Br J Cancer* 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4032580/) |
| `selumetinib` | MAP2K1/2 (MEK) | 12 / 80 | 0.97 | `BRAF` / `RAS` mutation — [*Mol Cancer Ther* 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2939826/) |

Compounds with a published determinant that the void panel dropped **on our label statistics rather than
on the literature**, and which the rebuild should reconsider: `sirolimus`, `neratinib`, `clofarabine`,
`cytarabine hydrochloride`, `gdc-0941`.

**Considered and set aside on coverage**, recorded so they are not re-proposed without checking it first:
`trametinib` and `at13387` (coverage only **0.46** of the 180 trainable lines). `kx2-391` was also
excluded, for a different and stronger reason: its signal was entirely the cell-line effect
([Corrections](corrections-and-dead-ends.md#kx2-391-carries-drug-specific-signal)).

> ⚠️ **`gemcitabine` was listed here and is not set aside — corrected 12.08.2026.** This passage
> recorded it at coverage **0.86**, below the 0.9 threshold, and told the rebuild to revisit it only if
> the threshold moved. On the current response source its coverage is **0.983**, the joint highest in
> the panel, and it is *in* the rebuilt panel
> ([Step 01](01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026),
> `outputs/panel/panel.csv`). The 0.86 was measured on CTRPv2's own 2015 distribution, which stopped
> being the target source on 11.08.2026; the threshold never moved. Left visible rather than deleted
> because this paragraph is written as instructions to the rebuild, and a reader following it would
> otherwise exclude a compound the rebuild had already selected.

**The determinants split by data modality, which makes any panel drawn from them a hypothesis test.**
Our input is expression only:

- **Expression-determined** — `methotrexate` (`SLC19A1`), `paclitaxel`/`vincristine` (`ABCB1`/`TUBB3`),
  `topotecan` (`SLFN11`), `tanespimycin` (`NQO1`), `dasatinib` (six-gene signature). The causal variable
  is *in* `X`, so these should be learnable.
- **Mutation-determined** — `selumetinib` (`BRAF`/`RAS` point mutations) and `afatinib` (amplification;
  expression explicitly failed to predict it in the pancreatic panel). The causal variable is **not** in
  `X` except through downstream expression, so weak per-drug ρ here is *expected* and is not evidence
  against the representation.

**Prediction to check on whatever panel is rebuilt:** ρ should be systematically higher in the first
group. If it is, that is a mechanistic validation of the approach; if the mutation-determined compounds
also score well, the model is picking up lineage rather than the stated mechanism and needs scrutiny.

> ⚠️ **Trap for the rebuild — read `moa_or_pathway`, not just `target`.** CTRP's `target` column is
> **empty** for `paclitaxel` and `vincristine` (and for `ml210` / `ml162`); they are annotated only via
> `moa_or_pathway`. A naive "keep rows with a validated target" filter silently drops four of the most
> informative compounds.

### Reporting convention — pooled estimate, fold spread, never fused

Established on the Step-1 run and independent of which drugs are in the panel: the two estimators differ.
**Pooling** gives one correlation over ~150 held-out lines; the **fold-wise mean** averages five
correlations over ~30 lines each and sits slightly higher. Report the **pooled value as the point
estimate** and the **fold spread as the dispersion** — do not fuse them into a single `mean ± sd` as
though they came from one calculation. And keep fold spread distinct from drug spread: the first says how
much the result depends on which lines were held out, the second how unevenly the model performs across
compounds.

The Step-1 run itself, with its numbers and dispersion, is in
[Corrections](corrections-and-dead-ends.md#the-step-1-training-run-on-the-voided-panel) — it was computed on the voided panel.

### Benchmarked with the real DrEval package (`notebooks/analysis/evaluation/dreval_benchmark.ipynb`, 14.07.2026)

Not a re-implementation: `pip install drevalpy` (v1.5.1, <https://github.com/daisybio/drevalpy>), and our
data/model run through **their** `DrugResponseDataset`, **their** `split_dataset(mode="LCO")`, **their**
baselines and **their** `evaluate()`. `OncoMLP` is trained on the single cells of each fold's train lines
and predicts the held-out lines' cells, averaged back per cell line — scored on exactly the same pairs.

**LCO, 5-fold CV, the 5 learnable drugs, native `auc` units** (mean over folds; *normalized* = their
recipe: subtract the `NaiveMeanEffects` prediction from `y_true` **and** `y_pred`, then re-evaluate):

| Model | Spearman (raw) | **Spearman (norm.)** | **R² (norm.)** |
|---|---|---|---|
| `NaivePredictor` | 0.000 | 0.020 | −0.052 |
| `NaiveDrugMeanPredictor` | 0.197 | 0.000 | −0.002 |
| `NaiveCellLineMeanPredictor` | 0.000 | 0.020 | −0.052 |
| **`NaiveMeanEffectsPredictor`** | 0.197 | **0.000** | −0.002 |
| their `SingleDrugElasticNet` (scGPT) | 0.197 | **0.000** ❌ | −0.002 |
| their `SingleDrugRandomForest` (PCA) | 0.245 | 0.148 | 0.022 |
| their `SingleDrugElasticNet` (PCA) | 0.320 | 0.300 | 0.056 |
| their `SingleDrugRandomForest` (scGPT) | 0.468 | 0.438 | 0.178 |
| **`OncoMLP` (X_pca)** | 0.481 | **0.442 ± 0.071** | 0.178 |
| **`OncoMLP` (X_scGPT)** | **0.549** | **0.511 ± 0.085** | **0.224** |

1. **OncoMLP clears `NaiveMeanEffects` decisively** (normalized ρ = 0.511 vs 0.000) — the bar **half the
   published models in the DrEval paper fail**. Our normalized **R² = 0.224** is directly comparable to
   the numbers they report for their best models in LCO (**DIPK 11%**, **Random Forest 19%**).
2. **scGPT > PCA is confirmed externally** (+0.07 normalized, on *their* splits with *their* metrics) —
   an independent replication of the +0.075 ± 0.038 we measured ourselves.
3. **The single-cell MLP beats their line-level reference models on the same embeddings** (0.511 vs 0.438
   for `SingleDrugRandomForest` on scGPT; 0.442 vs 0.148 on PCA). This **qualifies the ridge result**
   above: against a *stronger* per-drug regressor on line-mean embeddings, the per-cell model does add
   something — small (+0.07) but consistent across both representations.

> **Note on LCO:** a held-out line is unseen, so `NaiveMeanEffectsPredictor` sets its cell-line effect to
> **0** and reduces to *global mean + drug effect* (hence it ties `NaiveDrugMean` here). Their normalized
> metric therefore removes the **drug** mean — the fix for the Simpson's-paradox artifact they describe.
> It does *not* remove the cell-line effect, because in LCO no honest predictor can know it.

⚠️ Still a **best-case subset**: the 5 drugs were selected by a filter that saw all 180 lines
(`learnability_filter`). DrEval fixes the *evaluation*, not our *selection*.

### Own-implementation check: what if the cell-line effect is also removed? (14.07.2026)

**Reference:** Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response
prediction models with DrEval*, **Nat. Commun. (2026)**. Their headline: *"deep learning models barely
outperform a naive model that predicts only the mean drug and cell line effects"* — about **half** of
published models fail to beat their `NaiveMeanEffectsPredictor`. Our setting **is** their **LCO**
(leave-cell-line-out) with per-drug evaluation, i.e. the split they recommend.

**The gap this exposed in our metric.** `auc_z` removes the **drug** mean but *not* the **cell-line**
mean. Some lines are simply sensitive to everything (σ of the line effect = **0.40**), so a model can
score a good per-drug correlation by learning *"this line is fragile"* — with zero drug-specific biology.
DrEval's normalized metric subtracts the mean-effects predictor from **prediction and truth**, then
correlates; what remains is **differential sensitivity only**. Mean effects are fit on **train lines
only**, inside each fold (`notebooks/outputs/dreval/dreval_normalized.csv`).

> ⛔ **The numbers below cannot be reproduced by the current code (12.08.2026).** They came from a
> **stricter, locally invented** variant that removed the **cell-line** effect as well, using held-out
> labels. That has no counterpart in DrEval's paper and was deleted
> ([why](../../scripts/archive/README.md)). `scripts/evaluation/dreval_normalize.py` still exists but
> now applies the paper's normalization only. ⚠️ And under **leave-cell-line-out that normalization
> removes only the drug effect**, because a held-out line's effect is unseen and therefore zero — so
> the claim below that ~80 % of the effect survives *the cell-line effect* is exactly what the paper
> metric does **not** test. Re-deciding this is review item 11.

| | raw ρ | **normalized ρ** | naive baseline |
|---|---|---|---|
| K=5 · PCA | 0.427 | **0.368** | 0.291 |
| K=5 · scGPT | 0.488 | **0.396** | 0.291 |
| K=545 · PCA | 0.378 | **0.297** | 0.291 |
| K=545 · scGPT | 0.430 | **0.323** | 0.291 |

✅ **~80% of the effect survives** — it is genuine drug-specific signal, not the cell-line effect. For
scale, DrEval report their best models (DIPK, Random Forest) explaining **11–19%** of differential
sensitivity in LCO; ours is ρ² ≈ **0.16**, the same ballpark, on the filtered subset.

**Per drug (K=5, scGPT) — one of our five is an artifact:**

| drug | raw | **normalized** | naive baseline |
|---|---|---|---|
| `ml162` | 0.655 | **0.587** | 0.196 |
| `1s,3r-rsl-3` | 0.591 | **0.530** | 0.178 |
| `dasatinib` | 0.563 | **0.548** | 0.269 |
| `cay10618` | 0.347 | **0.306** | 0.226 |
| `kx2-391` | 0.283 | **0.006** ⚠️ | **0.584** |

⚠️ **`kx2-391` collapses to zero**: its entire apparent signal *was* the cell-line effect — exactly the
artifact class DrEval describes, found in our own results. The other four (both GPX4/ferroptosis
inducers, dasatinib, CAY10618) are real.

> **Decision (14.07.2026): adopt `NaiveMeanEffects` (drug mean + cell-line mean) as the standard
> baseline, and report raw *and* normalized correlations.** The per-drug-mean null is too weak; even
> ridge-on-line-means does not control for the cell-line effect. Any future claim must clear the
> normalized bar.

---

> **These are a best-case diagnostic, not headline numbers** — the drug subset was selected using all
> 180 lines, val and test included. Train-only selection inside each fold is what would make them
> reportable, and it remains blocking:
> [Corrections](corrections-and-dead-ends.md#the-1307-five-drug-numbers).

### Which change actually produced the gain?

Three things changed at once between the null result and the working one — the target, the head count and
the measurement. Isolated on the same drugs throughout, the target switch dominates (+0.29 PCA / +0.64
scGPT), honest out-of-fold measurement adds ~+0.1, and drug filtering only ~+0.06. The full table, and the
correction to an earlier claim that credited the curve fit rather than the standardization:
[Corrections](corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result) and
[Corrections](corrections-and-dead-ends.md#the-curve-fit-preserves-signal-the-dose-average-destroys).

### Gene-set sweep — heads-beating vs gene count (incl. all_genes, 28.06.2026)

> ⚠️ **12.08.2026 — after the sweep this comparison stops being like-for-like, and that is a third and
> separate problem from the two banners below.** Those say the numbers are superseded and that the
> `all_genes` point is mislabelled for scGPT. This one is about the *re-run*: **R1 re-embeds `hvg5000`
> and `all_genes` only** (decided 12.08.2026, Selin — the middle option, covering every number the
> report quotes; scGPT embedding is the expensive step, which is why the scope had to be fixed before
> R2). This table spans `hvg1000/2000/3000/5000`, so once R2 lands it will **mix one re-embedded variant
> with three embedded by the older code** — before the gene-symbol repair (4,576 → 4,704 in-vocab genes
> at `hvg5000`), before `gen_embeds.py` was seeded, and before the `ddof=1` harmonization. The gene-set
> axis would then vary the embedding code alongside the gene count, which is precisely the confound the
> sweep exists to exclude.
>
> **So any conclusion drawn across these points needs this stated, or the three remaining variants
> re-embedded as a top-up first.** Rejected alternatives, recorded so the cost is not rediscovered:
> re-embedding **all five** keeps the sweep like-for-like but is the longest run; **`hvg5000` alone** is
> cheapest but voids the sweep entirely *and* leaves every `all_genes` number in the report stale. Full
> decision: [TODO](../TODO.md), R1.

> ⛔ **03.08.2026 — the numbers in this table are superseded.** They were produced on the retired
> **`mean_pv`** target and cached at `outputs/archive/training_545_mean_pv/hvg_sweep.csv`. The sweep
> moved to `notebooks/analysis/qc/verify_variants.ipynb` §9 and was re-targeted to **`auc`**,
> which no longer reads that cache — so the sweep currently has **no live numbers**. The table is kept
> as the record of what was believed on 28.06.2026; do not quote it as current. Two further caveats
> on it are in [Step 02](02-preprocessing-and-embeddings.md#decision--one-seeded-draw-at-1200-all_genes-is-a-sanity-check-03082026):
> the `all_genes`/scGPT column came from unseeded embeddings, and the PCA column will move again once
> the pending `add_pca.py` changes land.

> ⛔ **05.08.2026 — the `all_genes` point does not mean what its label says, for scGPT only.** The
> `max_length=1200` cap binds in every cell at `all_genes` and in a single cell at `hvg5000`; per-variant
> counts are in [Step 02](02-preprocessing-and-embeddings.md#why-hvg-5000-is-the-default-03082026). So
> the four HVG points are genuine, and the PCA column is unaffected throughout — PCA reads every gene it
> is given at every point.
>
> **What the `all_genes` scGPT point therefore does support (B2, 05.08.2026).** It is not "fewer genes":
> at `all_genes` scGPT got roughly **twice** as many genes as at `hvg5000`, drawn at random from the
> in-vocab set instead of selected by dispersion. The flat result across the two is a real finding,
> stated narrowly: **doubling the gene count while randomising which genes are chosen buys nothing over
> half as many dispersion-selected genes.** It is *not* evidence about scGPT and the full transcriptome,
> which was never fed to it.
>
> ⚠️ `hvg1000`–`hvg3000` have no measured expressed-gene counts, but they cannot reach the cap. The HVG
> sets are **strictly nested** — `hvg1000 ⊂ hvg2000 ⊂ hvg3000 ⊂ hvg5000 ⊂ all_genes`, zero genes outside
> the larger set at every step, verified 05.08.2026 in
> `notebooks/analysis/qc/verify_variants.ipynb` §10a — so their per-cell counts are bounded by
> `hvg5000`'s, whose own maximum sits below the cap. `hvg1000` is settled independently and needs no
> check: 939 in-vocab genes cannot fill a 1,200-token sequence.

Does either rep have a preferred gene-set size?
`notebooks/analysis/qc/verify_variants.ipynb` §9 builds each variant
(1k/2k/3k/5k **plus `all_genes`**, full pipeline incl. scGPT re-embed; built by `1_preprocessing` §B until
12.08.2026, when the sweep build moved out of the numbered pipeline to `analysis/qc/`) and runs the same
**5-fold GroupKFold, test held out, all 545 drugs** — so the HVG-vs-all-genes comparison is
apples-to-apples under identical CV:

| Gene set | genes | `X_pca` heads-beat | `X_scGPT` heads-beat | Δmse (PCA / scGPT) |
|---|---|---|---|---|
| `hvg1000` | 1,000 | 207 ± 75 | 193 ± 83 | +0.00058 / +0.00060 |
| `hvg2000` | 2,000 | 203 ± 78 | 185 ± 84 | +0.00062 / +0.00064 |
| `hvg3000` | 3,000 | 216 ± 85 | 190 ± 83 | +0.00053 / +0.00063 |
| `hvg5000` | 5,000 | 210 ± 73 | 189 ± 94 | +0.00055 / +0.00074 |
| `all_genes` | 22,722 | 204 ± 86 | 184 ± 90 | +0.00058 / +0.00069 |

- **No sweet spot, and no all-genes advantage.** Both reps are **flat across the whole axis** (PCA
  ~203–216, scGPT ~184–193) — filtering does not help scGPT (contrary to the earlier hunch), and
  `all_genes` is **no better than HVG** ~~for either rep~~ **for PCA** (PCA's `all_genes` 204 sits
  mid-band, below hvg3000's 216; the earlier "PCA prefers all genes" is not reproduced). Val MSE
  ~constant (0.0105–0.0107) throughout. *(Amended 05.08.2026: for scGPT the `all_genes` point is a
  capped random draw, so it supports only the narrower claim in the block above, not a
  no-all-genes-advantage statement.)*
- PCA is marginally higher than scGPT at every gene count, but the ±73–94 fold spread overlaps
  completely — within noise at all sizes, consistent with the CV finding above.
- **Δmse > 0 at every gene-set size** for both reps: the model stays marginally *worse* than the
  per-drug-mean baseline regardless of how many genes feed it.

✅ On-plan: masked-loss multi-task, correctly gated behind a working single-task baseline,
with the cheap sanity baseline the plan's prototyping section calls for.

> ⚠️ **Key deviation — what "multi-task" means today:** the plan frames multi-task as
> **cross-database** (CTRPv2 + PRISM + GDSC heads). What's built is multi-task **across the
> 545 drugs of one database (CTRPv2)**. A legitimate *intermediate* step that validates the
> masked-loss machinery — but PRISM/GDSC are **not yet integrated**, so plan-Phase-3 is only
> half done. Don't read the 545-head run as "the multi-task goal is complete."

> ⚠️ **Split note — these are the matrix single-task cells, not Step 04's.** The K=1 rows above use
> `--drugs paclitaxel` on the **shared `split_ctrp`** (27 held-out lines), the same split as the
> K=545 runs — so within this table every comparison is apples-to-apples. They are **not** comparable
> to [Step 04](04-single-task-results.md)'s progression, which uses the separate `split_paclitaxel`
> (25 held-out lines). Different splits = different held-out cell lines.

---

### The first multi-task runs (26.05.2026) — which heads learned and which did not

The four runs that first established the masked-loss machinery, on `mean_pv` (run dirs are gitignored, so
these IDs are the only surviving trace; ledger row per run in `runs/runs_index.csv`):

| Run id | Rep | K | Best epoch | Best val MSE | Baseline mean MSE | Heads beating baseline |
|---|---|---|---|---|---|---|
| `20260526_132914_multitask_X_scGPT_subset_K1` | scGPT | 1 (paclitaxel) | 11 | 0.0412 | 0.0434 | 1 / 1 |
| `20260526_132952_multitask_X_pca_subset_K1` | PCA | 1 (paclitaxel) | 5 | 0.0393 | 0.0434 | 1 / 1 |
| `20260526_133012_multitask_X_scGPT_all_drugs` | scGPT | 545 | 7 | 0.0105 | 0.0097 | **142 / 545** |
| `20260526_133112_multitask_X_pca_all_drugs` | PCA | 545 | 6 | 0.0112 | 0.0097 | 97 / 545 |

Two things from these runs are still useful, because they name *specific* heads:

- **The heads that consistently fail** are the same for both representations, and they are the
  low-coverage ones (n_val = 221): `brd-k30748066`, `vx-680`, `brd-k33514849`,
  `brd9876:mk-1775 (4:1 mol/mol)`, `bafilomycin a1`. These are the concrete candidates for dropping or
  down-weighting — the open question under *Levers* in [TODO](../TODO.md).
- **The largest single win in both representations is `gsk-j4`** — model MSE ≈ 0.000 against a baseline of
  0.011 at n = 221. Worth keeping as an existence proof that a multi-task head *can* fit a low-variance
  (cell line × drug) combination, so a failing head is not evidence that the architecture cannot fit
  small-n drugs.

⚠️ The paclitaxel K=1 rows here use the shared `split_ctrp` (6,497 val labels) and are **not** comparable
to [Step 04](04-single-task-results.md)'s progression on `split_paclitaxel` (5,035 val labels) — different
held-out lines. Within this table PCA (0.0393) beats scGPT (0.0412) on paclitaxel alone.

## Run versioning (26.05.2026)

**Run versioning** (`training_utils.create_run_dir` / `save_run`): every
`train_multitask.py` run writes a self-contained `runs/<timestamp>_<tag>/`:

- `config.json` — exact `TrainConfig`.
- `run_meta.json` — scope, rep, dataset sizes, hidden_dims, host/python/torch info, drug list.
- `history.csv` — epoch, train_mse, val_mse, lr.
- `summary.json` — best_val_mse, best_epoch, baseline-vs-model mean MSE, heads-beating count.
- `best_model.pt` — best-val-MSE state_dict.
- `per_drug_results.csv` — drug, model_val_mse, baseline_val_mse, delta, n_val.

Plus one row per run in `runs/runs_index.csv` (columns: run_id, tag, scope, rep, K,
n_train_cells, n_val_cells, best_epoch, best_val_mse, baseline_mean_mse, model_mean_mse,
n_beats_baseline, n_total_heads, started_at, finished_at). `runs/` is gitignored.

✅ On-plan: satisfies "retain every working version + data to re-run + results, even
suboptimal ones."

The full 545-head run is reproduced with `train_multitask.py --use-rep {X_scGPT|X_pca}` (omitting
`--drugs` selects all K).
