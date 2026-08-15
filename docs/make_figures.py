"""Generate the OncoTox slide/doc graphics into ``docs/figures/``.

**Pure drawings — always reproducible, always built:**

  pipeline_overview.png    status of the whole project against the plan (steps 01-08)
  loss_01_objective.png    what the objective is made of

**Plots of committed results — built from an artifact, skipped if it is absent:**

  q2_instrument.png        stage 2 and stage 7 of the Q2 instrument, as distributions
  lpo_bias.png             mean per-drug Spearman under DrEval leave-pairs-out, by predictor

"Pure" means no expression matrix and no plotting of results: the layout is drawn, not measured.
Both still read the small committed CSVs the module derives its *counts* from (``panel.csv``,
``literature_panel_candidates.csv``, ``splits/split_ctrp.csv``, ``panel_heads_summary.csv``), which
is rule 4 doing its job rather than an exception to it -- a status figure carrying a typed count is
the failure this module has now had four times.

FIGURE CONVENTIONS — READ BEFORE WRITING ANY TITLE, CAPTION OR ANNOTATION
=========================================================================

Set by Selin 12.08.2026: *"stop making LLM like figures and titles. this is a scientific report."*
Recorded here rather than in a message so it does not have to be restated to whoever comes next.
Three rules, and every one of them was broken somewhere in this file when they were written.

1. **A title names what is plotted. It never states what the result means.**
   Not ``"The weighting fired — and the ranking did not move"``.
   But  ``"Prediction spread and per-drug Spearman across density-weighting levels"``.
   The test is not "is it true" — a conclusion-shaped title is the wrong register even when it is
   correct, and a title that would still hold after the rerun is still wrong. A title that *cannot*
   survive a rerun is worse: it is a claim in a PNG, where it cannot be qualified, cited or
   retracted. See ``loss_01_objective.png``, which asserted a squared loss for weeks and could not
   fail because it reads no data.

2. **Annotations label; they do not interpret.**
   ``"high AUC = resistant"`` is a legend and belongs. ``"above the line: the model hedges less"``
   is an interpretation and does not — interpretation goes in the report text, where it can be
   argued, sourced and disagreed with. The same applies to an argument set as a caption: if a
   sentence would need a citation in prose, it needs one here too, and therefore belongs there.

3. **Every mark earns its place.** Axes, data, and the minimum labelling needed to read them.
   Circles, arrows, heat strips, rounded colour-filled callout boxes, italic grey asides and
   explanatory sub-captions are infographic register, not scientific figures. If a mark's job is to
   look explanatory rather than to carry information the reader needs, it goes.

A fourth, from an unrelated defect class found the same day: **never hardcode a count in a label.**
Derive it — ``len(PANEL)``, a row count, a config value. ``"the 8 heads"`` was already wrong when it
was found (the rebuilt panel has 11); ``180 trainable`` became 181 at the sweep and is now
derived from the committed split (:func:`_n_trainable_lines`). A number typed
into a string cannot be checked by anything.

AUDIT AGAINST THOSE RULES — 12.08.2026, not yet applied
-------------------------------------------------------
Every title, caption, annotation and axis label in this file, judged twice: *does it assert a
result?* and *is it a scientific caption or a headline?* **Nothing below has been rewritten.** The
figures are held until R4, and the ⛔ entries cannot be rewritten before then anyway — a title naming
what the run showed has to wait for the run. A checklist, so the pass is not a rediscovery.

**Keyed by function and by the opening words of the string — never by line number.** The first
version of this audit was keyed by line, and the commit that recorded it is the commit that broke
it: inserting this docstring pushed every cited site down by about 76 lines. Line numbers into a
file under active edit decay silently, and a checklist that sends the reader to the wrong place is
worse than none, because it will be trusted. ``grep -n`` on the quoted words finds each site in any
ref.

⛔ **Asserts a result — rewrite from what the run shows, not before:**
  ✅ ``build_pipeline_flow``   ~~"one seed — only scGPT clears the ridge control"~~
        **Closed 14.08.2026 (Selin): the claim is dropped and the caption names what is plotted** —
        "mean per-drug Spearman over N drugs, N seeds — each MLP against its own ridge". Both halves
        were wrong, in different ways. *"one seed"* was false about the bars, which average three
        (``panel_corr``); stage 7's scatter is one seed on purpose, stage 8 never was. And the claim
        half, while true of the plotted bars, put a **+0.0096** margin — against a seed sd of 0.0053
        — into a PNG, where it can be neither qualified nor retracted; the tallest bar on the panel
        is in fact the *PCA ridge* (0.2767), so "only scGPT clears" reads as "scGPT wins" and is not
        that. The finding lives in ``docs/steps/05``, sourced and disputable.
  ✅ ``build_loss_effect``  ~~"The weighting fired — and the ranking did not move"~~ ·
     ~~"above the line: the model hedges less"~~ · ~~"on the line: no gain, no loss"~~
        **All three closed in `8fd1c85`; this block was not updated and still listed them as
        outstanding until 14.08.2026.** The title is now descriptive — *"Prediction spread and
        per-drug Spearman, unweighted against density-weighted"* — and the two interpretive
        annotations are gone, leaving only the ``dashed: y = x`` legend. Verified twice over, by
        rendering ``loss_03_effect.png`` and looking at it, and by grepping the source: the only
        surviving occurrences of all three strings are in this audit block and in two comments at
        the call site recording what was removed.
        ⚠️ **This is the third audit entry here to outlive its defect** — after the head-count
        entries and the ``"Correct as written; do not sweep in"`` line about the PLOTTED disclaimer.
        The pattern is not that the entries are wrong when written; it is that **nothing re-reads
        this block when a figure changes**, so it decays in exactly one direction: toward claiming
        more is broken than is.

⛔ **Asserts a swept arm as a fact:**
  ✅ ``draw_architecture``   ~~"per-cell MLP · trained with the density-weighted masked MSE"~~
        **Fixed 14.08.2026.** It stated two arms as facts at once, and the bars beneath it are
        ``alpha = 0`` — so the figure asserted a weighting it does not show. The subtitle now names
        only what does not vary (the masking); *which* arm the bars are is stated once, with them,
        from ``EXAMPLE_ARM``.
  ✅ ``build_pipeline_flow`` ~~"unscreened pairs dropped, rare response values weighted up"~~
        **Fixed 14.08.2026 (Selin).** Same class as the one above: *"weighted up"* is true only for
        ``alpha > 0``, and every panel on that figure plots ``alpha = 0``, where ``W`` is identically
        1 and nothing is weighted at all. The caption now separates what the objective **always**
        does (masking) from what is **swept** (the weight), and reads ``alpha`` off ``EXAMPLE_ARM``
        so it follows the arm the figure draws. This matches the stance ``LOSS_TEX_MACROS`` already
        takes deliberately — ``\\ell`` and ``\\alpha`` stay symbolic there so the figure cannot
        pre-empt the R4 comparison — which the caption had been quietly contradicting.

⚠️ **Headline register — true, but the wrong form:**
  ``draw_architecture``     "Model architecture — one cell in, one AUC per panel drug out"
        The clause after the dash is a strapline; "Model architecture" is the caption.
  ``draw_loss_objective``   "The objective — a weighted, masked mean error"     (written 12.08, mine)
  ``draw_loss_objective``   "per-drug scaling belongs here, in the loss, rather than in the labels"
  ✅ ``draw_architecture``  ~~"The target is uncentred: each head's bias is initialized to …"~~
        **MOVED TO METHODS 14.08.2026 (Selin), and gone from the figure.** It was four lines of
        reasoning set as a caption, which is rule 2 exactly: a sentence needing a citation in prose
        needs one here too, and therefore belongs there. It now sits in
        ``report/sections/03_methods.tex`` §Representation and model with the three citations it
        always needed and a PNG could not carry — Lin et al. 2017 §4.1 (prior init), Loshchilov &
        Hutter 2019 (why an Adam-era decay coefficient does not carry to AdamW), and the
        ``transformers`` library (the bias/LayerNorm exemption). The figure keeps the PLOTTED
        provenance line, which is a *different passage* ten lines away in the same function and has
        been confused with this one once already.
        ⛔ **Before it moved, this entry called the passage "Accurate", and it was not — both defects
        were repaired first so the corrected text is what travelled (moving a passage is not a way
        to stop owning what it says):**
          · **"(~0.7)"**, a single hardcoded value for a per-drug quantity spanning **0.58–0.99**.
            It was below the panel's own mean (0.790), well below the full catalogue's (0.893), and
            it contradicted ``OncoMLP.init_head_bias_``'s own docstring — *"``auc_cc`` sits near
            0.9"*. A leftover from the retired winsorized ``auc``. Now derived, and as a range,
            because the *shape* was wrong too: the caption argued for a **per-head** initialization
            using a number that implied one head would do (:func:`_panel_drug_mean_range`).
          · **the causal clause** — *"excluded from weight decay, **because** … decay would pull it
            to 0"*. Both mechanisms are real, but ``TrainConfig.weight_decay`` is **0.0** and so is
            the CLI default, so the decayed group is decayed by zero and the exemption exempts the
            biases from nothing. The figure named a force no run applies and presented it as what
            holds the bias in place. Nothing in it was false; **the register was** — an inert
            grouping read as a load-bearing one. Now stated as conditional rather than deleted,
            because the condition is live: ``TrainConfig``'s own comment records that whether this
            model needs weight decay is **open**, the evidence that closed it being void.
        **The lesson is the one two entries down:** "true, but the wrong form" and "correct as
        written" are both claims about code that moves, and both expired here on the same day.

⚠️ **Infographic register (rule 3).** ``draw_loss_objective``'s three rounded colour-filled callout
  boxes and grey italic asides; ``draw_architecture``'s circles, arrows and heat strips. Both mine
  to redraw at R4.

⚠️ **Hardcoded counts (rule 4).** Mostly discharged 14.08.2026; what is left is listed as left.
  ✅ ``draw_architecture``   ~~"the 8 heads are the 8 rows of one Linear(64 → 8)"~~ and the head
        layer's own ``"heads / 8 drugs"`` caption — **both derived from ``len(PANEL)``**. The far
        worse defect was underneath them and this checklist had not spotted it: ``counts = [6,5,4,8]``
        and ``ys = … range(8)`` drew **eight** heads, so ``zip`` truncated the output silently and
        ``dasatinib``, ``crizotinib`` and ``afatinib`` were missing from a figure titled *"one AUC
        per panel drug out"*. Spacing, radius and bar height now derive too, with the head *band*
        held fixed so the drawing around them does not move.
  ✅ ``build_pipeline``      ~~"K=545 · out-of-fold over 153 lines"~~ — **K is 534**, derived via
        :func:`_n_heads`; 545 is CTRPv2's catalogue, not the model's head count. ``180 trainable``
        was already derived (:func:`_n_trainable_lines`) and reads 181.
  ✅ ``build_pipeline_flow`` ~~"545 drugs  →"~~ — that grid is ``line_mask``, which is
        ``(lines × 534)``; taken from the array's own shape. ~~"8  the panel"~~ was already derived.
  ⚠️ **Still hardcoded, and left:** ``build_pipeline``'s "CTRPv2 545 drugs" and "* 190 = name-matches
        …" (both **correct** — the catalogue and the name-match count), "out-of-fold CV over 153
        lines"; ``build_pipeline_flow``'s "one cell line = 56–1,990 cells"; ``draw_architecture``'s
        "Dropout 0.5 · input dropout 0.1 · AdamW, early stopping (patience 10)", which restates
        ``TrainConfig`` fields — deriving those means importing torch into the figure build, which is
        a heavier change than this pass is scoped for.

✅ **Correct as written; do not sweep in:**
  ⚠️ ``draw_architecture``   ~~"PLOTTED: the superseded run … CURRENT PIPELINE: auc_cc, …"~~
        **This entry went stale on 14.08.2026 and is the cautionary one.** The sentence was true, and
        recorded here as true, when the bars were eight hardcoded values from the void run. Deriving
        them from ``panel_oof_predictions.csv`` made *every clause of it false at once* — the bars
        now **are** the current pipeline, so a disclaimer contrasting them with it says the opposite
        of the truth, and this list was still vouching for it. Replaced by a provenance line naming
        the plotted arm, every field read from ``EXAMPLE_ARM`` and ``PANEL``. **A "correct as
        written" entry is a claim about code that changes; it expires like any other.**
  ``draw_architecture``     "high AUC = resistant" and its sensitive counterpart — legends.
  ``build_pipeline``        "Results withdrawn" — a plan-status key.

**And one artifact that is not a figure**, written to ``report/`` rather than ``docs/figures/``:

  report/loss_objective.tex   the objective's equations, as LaTeX macros the report ``\\input``s

It is generated from ``LOSS_TEX_MACROS`` — the same strings ``loss_01_objective.png`` renders — so
that the maths in the report and the maths in the figure would have one source and could not drift
apart.

✅ **Wired in 14.08.2026 (Selin), and it had never been.** Until then ``report/main.tex`` carried no
``\\input{loss_objective}`` and no section cited any of the five macros — the file was regenerated on
every figure build, committed, and read by nothing, while both this docstring and the generated file's
own header asserted that the report used it. §Representation and model now sets ``\\LossObjective`` and
``\\LossWeightMatrix`` as equations, so the guarantee is real for the first time.

⚠️ **Only two of the five are cited, deliberately.** ``\\LossWeightFn`` and ``\\LossDensity`` describe
the density weighting, which is **off in every arm the report reports** (``alpha = 0``, so ``w_j``
is identically 1). Typesetting the weighting apparatus in Methods would give formal prominence to a
mechanism that never fires in the results — the same register error this file's own audit block
catalogues in its captions. They belong wherever the ``alpha`` sweep is reported. ⚠️ **Still missing:
the regenerate-and-diff pre-merge check** named below as owned by the gate session. It does not exist —
``loss_objective`` appears in no script under ``scripts/gate/`` — so nothing yet catches a hand-edit of
the generated file. Decided
by Selin 12.08.2026 over rendering the formula into the PNG (raster, fonts would not match the
report body, not referenceable by LaTeX) and over a PGF/vector figure (fonts match, but it puts a
LaTeX installation in the figure build path). It is committed so a clean checkout still compiles
the report; a regenerate-and-diff pre-merge check that fails on drift is owned by the gate session.

**Derived from data — all four BUILD again as of 14.08.2026** (they were skipped from 12.08.2026,
and the archived copies under ``docs/figures/archive/`` are the superseded ones):

  pipeline.png             the pipeline as a picture, stage by stage
  model_architecture.png   one cell in, one AUC per panel drug out
  loss_02_weights.png      one drug's label density and the weight curve fitted to it
  loss_03_effect.png       what the weighting did, per drug: spread up, ranking flat

The four read either ``figure_data.npz`` -- rebuilt from the ``auc_cc`` targets h5ad -- or a
training output under ``notebooks/outputs/panel/``. Both now exist on the current target: the
13.08.2026 re-run produced the panel artifacts, and the last two skips (``pipeline.png`` and
``model_architecture.png``, held by ``_example_matches_panel``) cleared on 14.08.2026 when
:func:`_example_predictions` began deriving the example vector from ``panel_oof_predictions.csv``.

**They are still not built from the retired ``auc`` h5ad that is on disk**, and the skip machinery
stays: a figure that can only be produced from a target the pipeline no longer writes is not
reproducible by a standard run -- and one that renders anyway is worse than one that is absent,
since nothing on its face says which target it came from. Each skip prints its reason.

The drug panel is read from ``notebooks/outputs/panel/panel.csv`` rather than duplicated here. It
was a hardcoded 8-compound list until 12.08.2026, which is how it went stale when the panel was
rebuilt -- the two share only three compounds.

The figures carry a caption at most; the argument behind them lives in the prose. Note that
``docs/progress_report_*.md`` is **untracked by design**, so a fresh clone will not have it --
the tracked source for every claim is ``docs/steps/`` (§4 the weighting, §9 the Step-1 result
correspond to Step 03 and Step 05) and
``docs/steps/03-model-and-training-design.md`` (architecture, loss, target).

Run:  uv run docs/make_figures.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG = HERE / "figures"
PANEL_OUT = ROOT / "notebooks" / "outputs" / "panel"
MIL_OUT = ROOT / "notebooks" / "outputs" / "mil"
# The 27.07.2026 training run's outputs, archived 12.08.2026: they were produced on the void 8-drug
# panel and cannot be recreated by a standard run, so they moved out of outputs/panel/. Only the
# archived figures read them, and each guards on existence first.
LEGACY_PANEL = ROOT / "notebooks" / "outputs" / "archive" / "panel_void_8drug"

GREEN = "#2e7d32"; GREEN_FILL = "#c8e6c9"
AMBER = "#b8860b"; AMBER_FILL = "#ffe9b3"
RED = "#c62828"; RED_FILL = "#ffcdd2"
BLUE = "#1f6fb2"; BLUE_FILL = "#dbe7f3"
GREY = "#777777"; GREY_FILL = "#e8e8e8"
INK = "#1a1a1a"
MUTED = "#52514e"

def _panel() -> list[str]:
    """The drug panel, read from ``outputs/panel/panel.csv`` — the file that defines it.

    Hardcoded as eight compounds until 12.08.2026, which is how it went stale: the panel was
    rebuilt on FDA approval and published determinants and this copy was not. The rebuilt panel
    shares only ``paclitaxel``, ``dasatinib`` and ``afatinib`` with the old one, so it was never an
    extension of it. ``drug_key`` is the spelling the response data uses (*Cisplatin* is
    ``platin``), which is what the h5ad is indexed by.
    """
    import pandas as pd

    return pd.read_csv(PANEL_OUT / "panel.csv")["drug_key"].tolist()


def _panel_display() -> dict[str, str]:
    """``{drug_key: name a reader recognises}`` — for axis labels, never for indexing data.

    ``panel.csv`` carries both spellings: ``drug_key`` is what CTRPv2 and every artifact are keyed
    by, ``drug`` is the compound's name. They differ enough to matter on a slide -- ``platin`` is
    *Cisplatin* -- and an audience reading an internal key off an axis is a defect, not a detail.
    The salt/formulation suffix is dropped (*Imatinib mesylate* -> *Imatinib*): it is accurate and it
    is noise at this size.
    """
    import pandas as pd

    p = pd.read_csv(PANEL_OUT / "panel.csv")
    drop = (" hydrochloride", " mesylate", " tosylate")
    out = {}
    for key, name in zip(p["drug_key"], p["drug"]):
        for suffix in drop:
            name = name.replace(suffix, "")
        out[key] = name
    return out


def _n_candidates() -> int:
    """How many FDA-approved anticancer compounds CTRPv2 actually screened.

    Read from ``outputs/panel/literature_panel_candidates.csv`` — the set the panel was drawn from,
    written by ``2_drug_selection`` §3, which states the funnel itself: *"panel: 11 drugs from 57
    candidates from 150 FDA-approved drugs"* (cell 11) and *"57 of 150 FDA-approved anticancer drugs
    were screened by CTRPv2 (57 of 120 small molecules)"* (cell 7).

    The tier was hardcoded ``173`` and labelled "FDA / clinical" until 13.08.2026. **173 is not a
    count of compounds at all** — it is a value of ``n_auc_cc``, the number of *cell lines* with an
    ``auc_cc`` measurement for one drug, visible in that notebook's cell 9 output. A per-drug
    cell-line count had been placed in a funnel of drugs.

    Derived rather than typed for the same reason as :func:`_panel`: ``n`` also sets the tier's drawn
    width, so a stale literal made the funnel the wrong *shape*, not only the wrong number. Both
    tiers now fail together if either moves.
    """
    import pandas as pd

    return len(pd.read_csv(PANEL_OUT / "literature_panel_candidates.csv"))


def _n_heads() -> int:
    """How many drugs the multi-task model actually has a head for — **534, not 545**.

    The two are different quantities and the figures conflated them. **545** is CTRPv2's catalogue,
    and it is right wherever the catalogue is meant — the drug funnel narrows *from* it. **534** is
    what survives ``ctrp_to_h5ad.DEFAULT_MIN_CELL_LINES`` (:file:`scripts/preprocessing/ctrp_to_h5ad.py:78`,
    applied at ``:246``): a compound is kept only if it was screened against at least 50 distinct
    SCP542-overlapping cell lines, and eleven compounds miss that. So the targets h5ad is
    53,513 x 534, ``uns['ctrp_drugs']`` has 534 entries, and the model has 534 heads.

    ⚠️ **The threshold of 50 carries no source** — see
    :file:`docs/steps/01-datasets-and-harmonization.md` §*The 50-cell-line drug cut*, which records
    it as arbitrary rather than leaving it to read as principled. Note it is **not** the other 50 in
    the pipeline (lines with fewer than 50 assigned *cells*); the two are unrelated cuts that happen
    to share a number.

    Read from ``outputs/panel/panel_heads_summary.csv``, the artifact behind Step 05 §D's 534-head
    run, rather than typed — rule 4, and the reason the count was wrong here in the first place.
    """
    import pandas as pd

    heads = pd.read_csv(PANEL_OUT / "panel_heads_summary.csv")["heads"].unique()
    if len(heads) != 1:                  # never quietly pick one
        raise ValueError(f"panel_heads_summary.csv mixes head counts: {sorted(heads)}")
    return int(heads[0])


def _n_trainable_lines() -> int:
    """Cell lines with post-QC response measurements, read from the committed split.

    Derived rather than typed for the reason the other two tiers are: it was hardcoded ``180`` and
    became **181** on 13.08.2026, when ``ctrp_to_h5ad`` deduplicated the experiment table and resolved
    the ``H292`` alias. A status figure carrying a stale count is worse than one carrying none,
    because it is read as the current state of the project.

    ``splits/split_ctrp.csv`` is the committed frozen split — one row per trainable line, 126 train /
    27 val / 28 test — so it is the artifact that defines this number rather than a copy of it.
    """
    import pandas as pd

    return len(pd.read_csv(ROOT / "splits" / "split_ctrp.csv"))


#: The drug panel — the order every figure uses. Read from panel.csv, not maintained here.
PANEL = _panel()

#: The middle tier of the drug funnel — see :func:`_n_candidates` for why it is not 173.
N_CANDIDATES = _n_candidates()

#: Trainable cell lines — see :func:`_n_trainable_lines` for why it is not 180.
N_LINES = _n_trainable_lines()

#: Model output heads — see :func:`_n_heads` for why it is 534 and not CTRPv2's 545.
N_HEADS = _n_heads()


def box(ax, x, y, w, h, title, lines, edge, fill, title_color=None, dashed=False,
        title_size=10.5, body_size=8.0):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=1.4",
        linewidth=2.2, edgecolor=edge, facecolor=fill,
        linestyle="--" if dashed else "-", mutation_aspect=1.0, zorder=2))
    ax.text(x + w / 2, y + h - 3.0, title, ha="center", va="top",
            fontsize=title_size, fontweight="bold", color=title_color or INK, zorder=3)
    if lines:
        ax.text(x + w / 2, y + h - 7.6, "\n".join(lines), ha="center", va="top",
                fontsize=body_size, color=INK, zorder=3)


def arrow(ax, x1, y1, x2, y2, color=INK, dashed=False):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
                 linewidth=2.0, color=color, linestyle="--" if dashed else "-", zorder=1))


# ============================================================ 1) pipeline overview
def build_pipeline():
    fig, ax = plt.subplots(figsize=(17.0, 9.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    W, H = 21.5, 20
    XS = [3.5, 28.0, 52.5, 77.0]
    ROW_A, ROW_B = 58, 28

    ax.text(50, 98, "OncoTox Pipeline — Status Overview", ha="center", va="top",
            fontsize=17, fontweight="bold", color=INK)
    ax.text(50, 93.5, "as of 2026-08-14   ·   reference: project_planning_v2.pdf   ·   steps: docs/steps/",
            ha="center", va="top", fontsize=9.5, color=GREY)
    handles = [
        mpatches.Patch(facecolor=GREEN_FILL, edgecolor=GREEN, label="Done / on-plan"),
        mpatches.Patch(facecolor=AMBER_FILL, edgecolor=AMBER, label="Addition beyond plan"),
        mpatches.Patch(facecolor=GREY_FILL, edgecolor=GREY, label="Results withdrawn"),
        mpatches.Patch(facecolor=RED_FILL, edgecolor=RED, label="Not started (planned)"),
    ]
    ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.885),
              ncol=4, fontsize=9.5, frameon=True, framealpha=0.9)

    box(ax, XS[0], ROW_A, W, H, "01 · Datasets & harmonization",
        ["SCP542 53,513 cells x 22,722 g", "CTRPv2 545 drugs · target: auc_cc",
         f"overlap 190* lines · {N_LINES} trainable"], GREEN, GREEN_FILL)
    box(ax, XS[1], ROW_A, W, H, "02 · Preprocessing & embeddings",
        ["scGPT X_scGPT = 512-d", "gene-set sweep 1k-5k + all_genes",
         "X_pca = 512-d · cancer-type UMAPs"], GREEN, GREEN_FILL)
    box(ax, XS[2], ROW_A, W, H, "03 · Model & training design",
        ["per-cell input -> viability", "masked MSE · matched (128,64) MLP",
         "PCA & scGPT both 512-d"], GREEN, GREEN_FILL)
    # Results withdrawn 12.08.2026 -- the target was replaced, the panel rebuilt, and the
    # representations predate the preprocessing corrections. These boxes named specific numbers
    # ("best scGPT val MSE 0.0336"; "target fix (auc_z): rho ~0 -> ~0.4") until then. They say so
    # rather than quietly dropping them: a status figure with the number silently removed reads as
    # a stage that was never measured, which is not what happened.
    box(ax, XS[3], ROW_A, W, H, "04 · Single-task baseline",
        ["paclitaxel, leak-free split", "results WITHDRAWN 12.08.2026",
         "1 DB · 1 score · 1 drug"], GREY, GREY_FILL, title_color=GREY)

    # Re-measured 13.-14.08.2026 and no longer withdrawn: 4a sections A/C/D/E and 4b ran on the
    # rebuilt 11-drug panel, the auc_cc target and the corrected early stopping, over three seeds.
    # The box states SCOPE and STATUS only -- what the run showed belongs in prose that can be cited
    # and disputed, not in a stage label of a pipeline diagram (this file's own audit rule).
    # K corrected 545 -> 534 on 14.08.2026, and derived with it. 545 is CTRPv2's catalogue; the
    # model's head count is what clears the 50-overlapping-lines cut in ctrp_to_h5ad, which is 534
    # and is what docs/steps/05 §D reports. See :func:`_n_heads`.
    box(ax, XS[0], ROW_B, W, H, "05 · Multi-task + fair eval",
        [f"K={N_HEADS} · out-of-fold over 153 lines", "re-run 13.08.2026 · 3 seeds",
         "results: docs/steps/05"], GREEN, GREEN_FILL)
    box(ax, XS[1], ROW_B, W, H, "06 · Cross-database  (MISSING)",
        ["CTRPv2 + PRISM + GDSC", "efficacy + toxicity heads",
         "the 'combine all' goal"], RED, RED_FILL, title_color=RED, dashed=True)
    box(ax, XS[2], ROW_B, W, H, "07 · XAI / interpretability  (MISSING)",
        ["feature importance -> drivers", "stretch goal", "not started"],
        RED, RED_FILL, title_color=RED, dashed=True)
    box(ax, XS[3], ROW_B, W, H, "08 · Foundation model  (HORIZON)",
        ["reusable pan-cancer FM", "fine-tune on clinical (binary)",
         "overarching main goal"], RED, RED_FILL, title_color=RED, dashed=True)

    # The band listed the retired auc_z target and the retracted learnability filter among the
    # additions. Both were removed on 12.08.2026 (Selin) rather than relabelled: this figure states
    # what the project has, and a retired target is not an addition it still carries. The record of
    # both is in docs/steps/corrections-and-dead-ends.md.
    BAND_Y, BAND_H = 5, 13
    box(ax, XS[0], BAND_Y, 94, BAND_H, "Additions beyond the written plan",
        ["512-d PCA (matched to scGPT)  ·  out-of-fold CV over 153 lines  ·  per-drug correlation metric  ·  "
         "gene-set sweep  ·  cancer-type UMAPs  ·  cell-line-grouped split (leak fix)  ·  run versioning\n"
         "ridge line-level control  ·  external benchmark against DrEval (drevalpy)  ·  "
         "nested early-stopping split"],
        AMBER, AMBER_FILL, title_color=AMBER)

    for i in range(3):
        arrow(ax, XS[i] + W, ROW_A + H / 2, XS[i + 1], ROW_A + H / 2)
    yb = ROW_A - 3.0
    arrow(ax, XS[3] + W / 2, ROW_A, XS[3] + W / 2, yb)
    ax.add_patch(FancyArrowPatch((XS[3] + W / 2, yb), (XS[0] + W / 2, yb),
                 arrowstyle="-", linewidth=2.0, color=INK, zorder=1))
    arrow(ax, XS[0] + W / 2, yb, XS[0] + W / 2, ROW_B + H)
    for i in range(3):
        arrow(ax, XS[i] + W, ROW_B + H / 2, XS[i + 1], ROW_B + H / 2, color=RED, dashed=True)

    ax.text(99.5, 1.2, f"* 190 = name-matches in CTRPv2's roster; {N_LINES} = lines with actual post-QC measurements",
            ha="right", va="bottom", fontsize=8, color=GREY, style="italic")

    out = FIG / "pipeline_overview.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


# ============================================================ 2) input + model + task
def _heat_strip(ax, xc, y0, y1, vals, cmap, w=2.4):
    """Vertical heatmap strip (a 'vector') of len(vals) cells."""
    ys = np.linspace(y0, y1, len(vals) + 1)
    cm = plt.colormaps[cmap]
    for i, v in enumerate(vals):
        ax.add_patch(Rectangle((xc - w / 2, ys[i]), w, ys[i + 1] - ys[i],
                     facecolor=cm(v), edgecolor="white", lw=0.4, zorder=3))
    ax.add_patch(Rectangle((xc - w / 2, ys[0]), w, ys[-1] - ys[0], fill=False,
                 edgecolor=INK, lw=1.3, zorder=4))


#: One real held-out cell line, scGPT, unweighted run (notebooks/outputs/archive/panel_void_8drug/panel_oof_predictions.csv,
#: fold in which SKES1_BONE was held out). Predicted vs measured AUC for the eight panel drugs. Used
#: instead of an invented vector so the figure shows the actual output scale -- including the visible
#: shrinkage (predictions span 0.37-0.75 against a measured 0.04-0.91).
#: The example cell line the architecture drawing shows. Unchanged from the original figure
#: (Selin, 27.07.2026) so that repointing it onto the rebuilt panel updates the DATA and not the
#: editorial choice; ``SKES1_BONE`` is still held out and still carries all eleven panel drugs.
EXAMPLE_LINE = "SKES1_BONE"

#: The arm those predictions come from: the transcript-level representation, unweighted, squared
#: error, first seed. One real run rather than a mean across seeds, because the figure exists to show
#: the actual output scale -- including the visible shrinkage -- and an average is not a run's output.
#:
#: **Chosen by Selin, 14.08.2026** -- recorded here because it is a display decision and this is where
#: it lives. The docs had said the repair was hers, an agent proposed the arm, and she confirmed it as
#: her choice while keeping ``EXAMPLE_LINE`` at her original 27.07.2026 pick. ⚠️ Note what it is: not
#: the project's best arm (that is ``X_pca``/``mae``/alpha=0.5 at 0.2824) but ``X_scGPT`` unweighted --
#: so the bars show the transcript-level arm's own output, and a reader must not take them for the
#: headline result. The PLOTTED line under the figure names all four fields for exactly that reason.
EXAMPLE_ARM = dict(rep="X_scGPT", alpha=0.0, loss="mse", seed=42)


def _example_predictions():
    """Real out-of-fold predictions for one held-out cell line, in ``PANEL`` order.

    **Derived, not typed (14.08.2026).** These were eight hardcoded values from the **void 8-drug**
    run, while the drawing labelled bar ``j`` with ``PANEL[j]`` -- so every bar carried another
    compound's name and four named compounds that were never in that run. The figure was skipped
    rather than corrected because repair needed out-of-fold predictions for the rebuilt panel, which
    did not exist until the 13.08.2026 re-run.

    Reading them from the artifact is what stops the defect recurring: the vector and the panel are
    now the same file's two views, so they cannot disagree. It is the third quantity in this module
    to be derived for that reason, after :func:`_panel` and :func:`_n_candidates`.
    """
    import pandas as pd

    oof = pd.read_csv(PANEL_OUT / "panel_oof_predictions.csv")
    sel = oof[(oof.rep == EXAMPLE_ARM["rep"]) & (oof.alpha == EXAMPLE_ARM["alpha"])
              & (oof.loss == EXAMPLE_ARM["loss"]) & (oof.seed == EXAMPLE_ARM["seed"])
              & (oof.cell_line == EXAMPLE_LINE)].set_index("drug")
    missing = [d for d in PANEL if d not in sel.index]
    if missing:                      # never silently draw a short vector
        raise ValueError(f"{EXAMPLE_LINE} has no out-of-fold row for {missing}")
    return ([round(float(sel.loc[d, "y_pred"]), 3) for d in PANEL],
            [round(float(sel.loc[d, "y_true"]), 3) for d in PANEL])


def _panel_drug_mean_range() -> tuple[float, float]:
    """``(min, max)`` of the **per-drug** mean ``auc_cc`` over the panel's cross-validated lines.

    What ``OncoMLP.init_head_bias_`` actually targets, and the point is that it is a *vector*: each
    head is started at **its own** drug's mean over the fold's fitting lines, which is the entire
    reason the initialization is per head rather than one constant.

    **Replaces a hardcoded "~0.7" (14.08.2026, Selin).** That number was wrong three ways at once —
    below the panel's own mean of 0.790, well below the full catalogue's 0.893
    (``panel_heads_oof.csv``), and in direct contradiction with ``init_head_bias_``'s own docstring,
    which says ``auc_cc`` *"sits near 0.9"*. It reads as a leftover from the retired winsorized
    ``auc``. Worse than the value, a single figure stood in for a quantity spanning 0.58–0.99, so the
    caption argued for a per-head initialization using a number that implied one head would do.

    Taken over the out-of-fold rows of the plotted arm, so the range describes the same lines the
    bars beside it are drawn from.
    """
    import pandas as pd

    oof = pd.read_csv(PANEL_OUT / "panel_oof_predictions.csv")
    sel = oof[(oof.rep == EXAMPLE_ARM["rep"]) & (oof.alpha == EXAMPLE_ARM["alpha"])
              & (oof.loss == EXAMPLE_ARM["loss"]) & (oof.seed == EXAMPLE_ARM["seed"])
              & (oof.drug.isin(PANEL))]
    means = sel.groupby("drug").y_true.mean()
    if len(means) != len(PANEL):     # never describe a partial panel as the panel
        raise ValueError(f"per-drug means cover {len(means)} of {len(PANEL)} panel drugs")
    return float(means.min()), float(means.max())


EXAMPLE_DRUGS = list(PANEL)
EXAMPLE_PRED, EXAMPLE_TRUE = _example_predictions()


def _example_matches_panel(name: str) -> bool:
    """True (and prints why) if the architecture example vector does not describe ``PANEL``.

    The bars are real out-of-fold predictions, so each one is a measurement and carries the name of
    the drug it was measured on. When the example vector and the panel disagree, drawing them
    together publishes each measurement under some other compound's name -- which no amount of
    caption disclaims. Skip-with-a-reason, the same convention the data-derived figures use.
    """
    if EXAMPLE_DRUGS == list(PANEL):
        return False
    print(f"  {name}: SKIPPED — the example prediction vector is the void 8-drug panel "
          f"({', '.join(EXAMPLE_DRUGS[:3])}…) and the current panel is {len(PANEL)} drugs "
          f"({', '.join(PANEL[:3])}…). Drawing them together labels each measurement with "
          f"another drug's name. Needs out-of-fold predictions for the rebuilt panel (R4).")
    return True


def draw_architecture(ax, *, compact: bool = False):
    """Draw the architecture into ``ax``.

    One drawing, two sizes: the standalone ``model_architecture.png`` uses ``compact=False``
    (with the mechanics spelled out underneath), stage 5 of ``pipeline.png`` uses ``compact=True``
    (the same picture, without the paragraphs that would be unreadable at panel size). Nothing is
    redrawn by hand for the small version, so the two cannot diverge.
    """
    # compact = a genuinely tighter drawing (layers closer, bars shorter, type scaled down),
    # not the same geometry shrunk -- otherwise the panel eats half the pipeline figure.
    fs = 0.70 if compact else 1.0
    xc, xe, lx, x0, span = ((4.5, 9.5, [16, 23, 29, 35], 49.0, 10.5) if compact else
                            (6.0, 14.0, [26, 39, 51, 62], 76.0, 17.0))
    # Non-compact lower bound: 14 -> 11 on 12.08.2026 to fit the "plotted" note under the
    # uncentred-target paragraph, 11 -> 10 earlier on 14.08.2026 when that paragraph grew to three
    # lines, and back to 14 once the paragraph moved to Methods the same day -- the note is now the
    # only thing below the legends and sits where the paragraph began. Compact is untouched: it is
    # embedded in pipeline.png, which has its own layout, and it draws neither note.
    ax.set_xlim(0, 62 if compact else 100); ax.set_ylim(23 if compact else 14, 47 if compact else 52)
    ax.set_aspect("equal"); ax.axis("off")

    if not compact:
        ax.text(0, 51, "Model architecture — one cell in, one AUC per panel drug out",
                ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
        # Was "trained with the density-weighted masked MSE (loss_01-03)" until 14.08.2026, which
        # stated two swept arms as facts: `alpha` is in {0.0, 0.5, 1.0} and the loss in {mse, mae},
        # and the arm actually drawn below is alpha = 0 -- every observed pair weighs exactly 1, so
        # the figure asserted a weighting it does not show. The subtitle now names only what does
        # not vary (the masking); which arm the bars come from is stated once, with the bars.
        ax.text(0, 48.0, "per-cell MLP · a masked objective — only screened (cell line, drug) pairs "
                         "contribute (loss_01–03)",
                ha="left", va="top", fontsize=8.5, color=GREY)

    # ---------- INPUT: one cell -> embedding vector ----------
    ax.add_patch(Circle((xc, 37), 2.4 * (0.85 if compact else 1), facecolor="#fde0c5", edgecolor="#d2691e", lw=1.8, zorder=3))
    for dx, dy in [(-0.8, 0.45), (0.6, -0.35), (0.15, 0.9), (-0.25, -0.8)]:
        ax.add_patch(Circle((xc + dx, 37 + dy), 0.5, facecolor="#d2691e", lw=0, zorder=4))
    ax.text(xc, 32.6, "single cell\n(scRNA-seq)", ha="center", va="top", fontsize=9 * fs, color=INK)
    arrow(ax, xc + 2.8, 37, xe - 2.0, 37, color=INK)

    _heat_strip(ax, xe, 30.5, 43.5, np.linspace(0.05, 0.95, 14), "viridis", w=2.4 * (0.8 if compact else 1))
    ax.text(xe, 29.6, "512-d embedding", ha="center", va="top", fontsize=9.5 * fs,
            fontweight="bold", color=BLUE)
    ax.text(xe, 26.9, "PCA  or  scGPT" + ("" if compact else "\n(frozen — never fine-tuned)"),
            ha="center", va="top", fontsize=8.2 * fs, color=INK)
    arrow(ax, xe + 1.6, 37, lx[0] - 3.5, 37, color=INK)

    # ---------- MODEL: MLP drawn as neurons ----------
    layers_x = lx
    # The first three layers are drawn truncated (a '⋮' says so); the head layer is drawn IN FULL,
    # one circle per panel drug, so its count is len(PANEL) and never a literal.
    counts = [6, 5, 4, len(PANEL)]
    cy, r = 37, 1.25 * (0.72 if compact else 1)
    sp = 3.0 if not compact else 2.3
    # Spacing, radius and bar height follow the panel size. All three were fixed at eight heads'
    # worth (2.05 / 1.62 apart, radius 0.85 / 0.55) until 14.08.2026, and `zip` then truncated the
    # drawing to the first eight of eleven drugs without a word -- so a figure titled "one AUC per
    # panel drug out" silently dropped dasatinib, crizotinib and afatinib. What is held constant is
    # the vertical BAND the heads occupy (7 x 2.05 = 14.35 units, 7 x 1.62 = 11.34 compact), because
    # it is the band the dashed trunk box and the layer captions below were laid out around: keeping
    # it fixed lets the head count change without moving anything else on the drawing.
    head_band = 14.35 if not compact else 11.34
    head_sp = head_band / max(len(PANEL) - 1, 1)
    head_r = head_sp * (0.85 / 2.05 if not compact else 0.55 / 1.62)
    pos = []
    for li, (lxx, n) in enumerate(zip(layers_x, counts)):
        s = head_sp if li == 3 else sp
        pos.append([(lxx, cy + (i - (n - 1) / 2) * s) for i in range(n)])
    for a, b in zip(pos[:-1], pos[1:]):
        for (x1, y1) in a:
            for (x2, y2) in b:
                ax.plot([x1, x2], [y1, y2], color="#bcd0e6", lw=0.5, zorder=1)
    for li, layer in enumerate(pos):
        rr = head_r if li == 3 else r
        for (x, y) in layer:
            ax.add_patch(Circle((x, y), rr, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.5, zorder=3))
    for lxx, n in zip(layers_x[:3], counts[:3]):  # '...' only where neurons are omitted
        ax.text(lxx, cy - ((n - 1) / 2) * sp - r - 0.5, "⋮", ha="center", va="top",
                fontsize=12 * fs, color=GREY)
    for lxx, t in zip(layers_x, ["input\n512", "hidden\n128", "hidden\n64",
                                 f"heads\n{len(PANEL)} drugs"]):
        ax.text(lxx, 25.8, t, ha="center", va="top", fontsize=9 * fs, color=INK)
    ax.add_patch(FancyBboxPatch((layers_x[0] - 3.5, 27.5), layers_x[3] - layers_x[0] + 7, 19,
                 boxstyle="round,pad=0.3,rounding_size=1.2",
                 fill=False, edgecolor=BLUE, lw=1.2, linestyle="--", zorder=0))
    if not compact:
        # Epoch count moved to the "plotted" line above: it is a property of the run shown, not of
        # the architecture, and the cap is itself under review (item 10 owns whether TrainConfig's
        # default of 25 moves to the 50 that 4a_percell_training passes). Naming it here made the
        # architecture caption go stale every time the cap moved.
        # "Adam" corrected to AdamW 14.08.2026: the optimizer has been `optim.AdamW`
        # (`training_utils.py::train_model`) since Selin's 12.08.2026 switch under review item 10,
        # and the label had not followed. It names the optimizer only -- `TrainConfig.weight_decay`
        # is 0.0, which is a *training* setting and belongs with the arm, not in the architecture
        # caption.
        ax.text(43.2, 22.6, "hidden block  =  Linear → LayerNorm → GELU → Dropout 0.5      ·      "
                            "input dropout 0.1      ·      AdamW, early stopping (patience 10)",
                ha="center", va="top", fontsize=8.2, color=INK)
        ax.text(43.2, 19.8, f"the {len(PANEL)} heads are the {len(PANEL)} rows of one "
                            f"Linear(64 → {len(PANEL)}) over a shared trunk — "
                            "there is no per-drug sub-network",
                ha="center", va="top", fontsize=8.2, color=GREY, style="italic")

    # ---------- OUTPUT: predicted raw AUC per drug ----------
    amax = 1.15                                # AUC 0 .. 1.15 maps to x0 .. x0+span
    cm = plt.colormaps["coolwarm"]             # low AUC = sensitive (blue) .. high = resistant (red)
    shade = lambda v: cm(np.clip(0.5 + (v - 0.5) / 1.2, 0, 1))  # white anchored at AUC 0.5
    # One row per head, in the same order and at the same heights as the head circles: PANEL[j],
    # EXAMPLE_PRED[j] and EXAMPLE_TRUE[j] are three views of one row and are all len(PANEL) long.
    bar_h = head_sp * 0.75
    ys = [cy + (i - (len(PANEL) - 1) / 2) * head_sp for i in range(len(PANEL))][::-1]
    for j, (yy, p, t) in enumerate(zip(ys, EXAMPLE_PRED, EXAMPLE_TRUE)):
        arrow(ax, layers_x[3] + 1.0, yy, x0 - (5.5 if compact else 6.5), yy, color="#9db6cf")
        ax.plot([x0, x0 + span], [yy, yy], color="#e4e4e0", lw=0.8, zorder=1)
        ax.add_patch(Rectangle((x0, yy - bar_h / 2), p / amax * span, bar_h, facecolor=shade(p),
                               edgecolor="white", lw=0.6, zorder=3))
        ax.plot([x0 + t / amax * span], [yy], marker="D", ms=3.4, color=INK, zorder=4)
        ax.text(x0 - 0.8, yy, PANEL[j], ha="right", va="center", fontsize=8 * fs, color=INK)
    for v in (0.0, 0.5, 1.0):
        ax.plot([x0 + v / amax * span] * 2, [ys[-1] - 0.9, ys[0] + 0.9], color=GREY, lw=0.6,
                linestyle=":", zorder=1)
        ax.text(x0 + v / amax * span, ys[-1] - 1.4, f"{v:.1f}", ha="center", va="top",
                fontsize=7.5 * fs, color=GREY)
    ax.text(x0 + span / 2, 26.6, "predicted AUC per drug", ha="center", va="top",
            fontsize=9.5 * fs, fontweight="bold", color=INK)
    if compact:
        ax.text(x0 + span / 2, 24.3, "low = sensitive · high = resistant",
                ha="center", va="top", fontsize=6.6, color=GREY)
        return
    ax.text(x0 - 4.5, 24.0, "low AUC = sensitive (the drug kills this line)", ha="left", va="top",
            fontsize=8, color=BLUE)
    ax.text(x0 - 4.5, 21.8, "high AUC = resistant", ha="left", va="top", fontsize=8, color=RED)
    ax.text(x0 - 4.5, 19.6, f"bars: out-of-fold prediction for {EXAMPLE_LINE}\n"
                            "◆ = the measured CTRPv2 value (never seen in training)",
            ha="left", va="top", fontsize=7.6, color=GREY)

    # Two corrections, 12.08.2026. The retired auc_z target is no longer named -- the point stands
    # on its own, and Methods no longer describes that target at all. And "the output layer is
    # excluded from weight decay" described exclude_output_from_decay, which exempted the whole
    # output Linear including its weight matrix while still decaying LayerNorm; it was replaced by
    # no_decay_bias_and_norm the same day.
    #
    # Two further corrections, 14.08.2026 (Selin), both inside the second clause.
    #
    # 1. THE NUMBER. "(~0.7)" is derived now -- see :func:`_panel_drug_mean_range` for why a single
    #    value was the wrong SHAPE of fact here, not merely the wrong one.
    #
    # 2. THE CAUSAL REGISTER, which is the substantive fix. The clause read "biases and LayerNorm
    #    are excluded from weight decay, BECAUSE ... decay would pull it to 0". Both mechanisms are
    #    real -- `init_head_bias_` runs in all three training paths, and `no_decay_bias_and_norm` is
    #    True by default (training_utils.py:88, grouping at :304-309, sourced to HuggingFace's
    #    Trainer.create_optimizer). But `TrainConfig.weight_decay` is 0.0 and so is the CLI default,
    #    so the "decayed" group is decayed by zero and the exemption exempts the biases from
    #    nothing. The sentence named a force that is not applied in any run this pipeline produces,
    #    and presented it as what keeps the bias in place -- so a reader took an inert grouping for
    #    a load-bearing one. Nothing in it was false; the register was.
    #
    #    It is stated as conditional rather than deleted because the condition is live: TrainConfig's
    #    own comment records that "whether this model needs weight decay is OPEN, not settled", the
    #    evidence that used to close it being void. The grouping is insurance against a decision not
    #    yet taken, and saying exactly that is both accurate today and still accurate if it is.
    # ---- the uncentred-target paragraph lived here until 14.08.2026 ----
    #
    # MOVED TO report/sections/03_methods.tex, §Representation and model (Selin's decision). It was
    # four lines of reasoning set as a caption, which is this module's rule 2: a sentence that would
    # need a citation in prose needs one here too, and therefore belongs in prose. In Methods it now
    # carries the three it always needed and could not have -- Lin et al. 2017 §4.1 for prior init,
    # Loshchilov & Hutter 2019 for why an Adam-era decay coefficient does not carry to AdamW, and
    # the transformers library for the bias/LayerNorm exemption it follows.
    #
    # Both defects the passage was carrying were repaired first and travelled with it, so what moved
    # is the corrected text: the "~0.7" became a derived RANGE (\PanelMeanLo..\PanelMeanHi in
    # report/results_numbers.tex, from :func:`_panel_drug_mean_range`), and the inert "decay would
    # pull it to 0" became the honest statement that the coefficient is zero and the grouping
    # therefore binds on nothing. Moving a passage is not a way to stop owning what it says.
    #
    # :func:`_panel_drug_mean_range` is deliberately KEPT even though no figure draws it now. It is
    # what re-derives the two macros when the run changes, and the macro file says so; deleting it
    # would leave the report's only copy of the number hand-maintained with nothing to check it.

    # Rewritten 14.08.2026, and it is a REPLACEMENT rather than an edit, because every clause of the
    # old sentence had gone false at once. It read "PLOTTED: the superseded run -- auc winsorized at
    # 1.1, the void 8-drug panel, 25 epochs. CURRENT PIPELINE: auc_cc, no winsorization, the rebuilt
    # 11-drug panel, 50 epochs." That was accurate when written and stopped being so the moment
    # `_example_predictions` repointed the bars onto `panel_oof_predictions.csv`: the bars ARE the
    # current pipeline now, so a warning contrasting them with it says the opposite of the truth.
    # (The module docstring listed this string under "Correct as written; do not sweep in" -- that
    # entry went stale with the repointing and is corrected there too.)
    #
    # What replaces it is not a warning but provenance: which of the swept arms these eleven bars
    # are, since `alpha` and `loss` are both arms and one figure can only show one of them. Every
    # field is read from EXAMPLE_ARM and PANEL, so it cannot drift from the vector it describes --
    # the same reason the vector itself is derived. Colour drops from RED to MUTED with the change
    # of job. Plain ASCII: matplotlib's default font renders the warning emoji as a tofu box.
    ax.text(0, 16.0,
            f"PLOTTED: one arm of the current sweep — target auc_cc, the {len(PANEL)}-drug panel; "
            f"{EXAMPLE_ARM['rep']}, alpha = {EXAMPLE_ARM['alpha']:g} (unweighted: every observed "
            f"pair weighs 1), {EXAMPLE_ARM['loss'].upper()}, seed {EXAMPLE_ARM['seed']}"
            "    ·    source: notebooks/outputs/panel/panel_oof_predictions.csv",
            ha="left", va="top", fontsize=8.0, color=MUTED)


def build_architecture():
    if _needs_data("model_architecture.png") or _example_matches_panel("model_architecture.png"):
        return
    fig, ax = plt.subplots(figsize=(17.0, 6.6))
    draw_architecture(ax)
    out = FIG / "model_architecture.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")



# ============================================================ shared data for the drawn-from-data panels
#: Line-level quantities the figures need, cached so nothing here depends on the 2.4 GB targets h5ad
#: (which lives outside the repo).  Everything is per *cell line* over the train+val split, which is
#: the population every fold's statistics are estimated from.
CACHE = FIG / "figure_data.npz"


def figure_data() -> dict:
    """``{line_mask, panel_auc, n_cells, fold, n_by_split, pca2}`` — from the targets h5ad once."""
    if CACHE.exists():
        with np.load(CACHE) as z:
            return {k: z[k] for k in z.files}

    import sys

    sys.path.insert(0, str(ROOT))
    import anndata as ad
    from sklearn.model_selection import GroupKFold

    from scripts.layout import DEFAULT_CTRP_SCORE, PipelinePaths

    # auc_cc ONLY, and deliberately no fallback (Selin, 12.08.2026). Until 11.08.2026 this read the
    # retired `auc` target and clipped at DEFAULT_WINSOR; both were removed from the pipeline, so
    # the call raised ValueError and the import above it ImportError. Neither ever fired, because
    # the committed cache short-circuits this function -- the breakage was invisible.
    #
    # The retired `auc` h5ad is still on disk and these figures COULD be built from it. They are
    # not: a figure that can only be produced from a target the pipeline no longer writes is not
    # reproducible by a standard run, and one that renders anyway is worse than one that is absent,
    # because nothing on it says which target it came from. If the auc_cc h5ad is missing, this
    # raises and the figures that need it are skipped by their callers.
    target_h5ad = PipelinePaths.build(None, "hvg5000", DEFAULT_CTRP_SCORE).targets_h5ad
    if not target_h5ad.exists():
        raise FileNotFoundError(
            f"{target_h5ad.name} does not exist, so the figure cache cannot be built.\n"
            f"  Run stages 1 and 3 (1_data, 3_representations) first.\n"
            f"  The retired `auc` h5ad is NOT used as a substitute -- see the comment above."
        )
    src = ad.read_h5ad(target_h5ad, backed="r")
    all_drugs = list(src.uns["ctrp_drugs"])
    kcol = [all_drugs.index(d) for d in PANEL]
    Y = np.asarray(src.obsm["Y_ctrp"], dtype=float)
    M = np.asarray(src.obsm["M_ctrp"], dtype=bool)
    groups = src.obs["Cell_line"].astype(str).to_numpy()
    split = src.obs["split_ctrp"].astype(str).to_numpy()
    eligible = np.isin(split, ["train", "val"])
    # how many *cell lines* sit in each split -- the numbers the split panel has to be honest about
    first = {ln: split[groups == ln][0] for ln in np.unique(groups)}
    n_by_split = np.array([sum(v == s for v in first.values())
                           for s in ("train", "val", "test", "unassigned")], dtype=int)
    # a 2-d view of the PCA space, subsampled -- the icon for "linear, unsupervised" in stage 4
    rs = np.random.default_rng(0)
    take = np.sort(rs.choice(src.n_obs, size=min(4000, src.n_obs), replace=False))
    pca2 = np.asarray(src.obsm["X_pca"], dtype=float)[take][:, :2]
    src.file.close()

    lines = np.unique(groups[eligible])
    line_mask = np.zeros((len(lines), M.shape[1]), dtype=bool)
    panel_auc = np.full((len(lines), len(PANEL)), np.nan)
    n_cells = np.zeros(len(lines), dtype=int)
    for i, ln in enumerate(lines):
        ci = np.flatnonzero((groups == ln) & eligible)
        n_cells[i] = ci.size
        line_mask[i] = M[ci].any(0)
        obs = line_mask[i][kcol]
        panel_auc[i, obs] = Y[ci][:, np.array(kcol)[obs]][0]

    # the project's fold partition (scripts.training.cv.grouped_folds), collapsed to one id per line
    idx = np.flatnonzero(eligible)
    pos = {ln: i for i, ln in enumerate(lines)}
    fold = np.zeros(len(lines), dtype=int)
    for f, (_, va) in enumerate(GroupKFold(n_splits=5).split(idx, groups=groups[idx]), start=1):
        for ln in np.unique(groups[idx[va]]):
            fold[pos[ln]] = f

    FIG.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE, n_by_split=n_by_split, pca2=pca2, line_mask=line_mask, panel_auc=panel_auc,
                        n_cells=n_cells, fold=fold)
    return figure_data()


#: Which loss the drawings hold fixed when they need one. `mse` is the level sections C and D hold
#: fixed, so the figures inherit the sweep's reference rather than introducing a fourth convention.
FIG_LOSS = "mse"

#: Which alpha level stands for "the weighting is on". `cv.py`'s own note: "today's weighted=True is
#: exactly alpha=0.5 because DEFAULT_ALPHA is 0.5", so this is a documented identity, not a choice
#: made here. `weighted` (bool) was replaced by `alpha` (numeric) in 059d548; this reader was missed.
FIG_ALPHA_ON = 0.5


def panel_corr():
    """``notebooks/outputs/panel/panel_per_drug_correlation.csv``, averaged over seeds.

    **Repointed 14.08.2026 (Selin) from ``outputs/archive/panel_void_8drug/``.** It read the void
    8-drug panel, which shares three compounds with the rebuilt eleven, so every consumer either
    skipped itself or would have plotted three points as though they were a complete comparison.

    The current file carries three seeds and two losses per (rep, drug) where the legacy one carried
    a single row, so it is reduced here rather than at each call site: the loss is held at
    :data:`FIG_LOSS` and the three seeds are averaged, which is the convention every number in
    ``docs/steps/05`` uses. Callers therefore still get one row per (rep, alpha, drug).
    """
    import pandas as pd

    d = pd.read_csv(PANEL_OUT / "panel_per_drug_correlation.csv")
    d = d[d.loss == FIG_LOSS]
    return (d.groupby(["rep", "alpha", "drug"], as_index=False)
             .agg({"spearman": "mean", "mse": "mean", "pred_std": "mean",
                   "true_std": "mean", "n_lines": "first"}))


def _n_seeds() -> int:
    """How many seeds :func:`panel_corr` averages — read from the file it averages.

    Rule 4, in its other direction: stage 8 of ``pipeline.png`` was captioned *"one seed"* while its
    bars averaged three, and nothing could catch it because the count was typed rather than taken
    from the data it described. Any caption stating the aggregation reads it from here.
    """
    import pandas as pd

    d = pd.read_csv(PANEL_OUT / "panel_per_drug_correlation.csv")
    return int(d[d.loss == FIG_LOSS].seed.nunique())


def _needs_data(name: str) -> bool:
    """True if ``name`` cannot be built right now; prints why.

    Added 12.08.2026. The figures split cleanly in two: some are pure drawings and always build,
    the rest are derived from the targets h5ad or from a training output. The derived ones are only
    reproducible by a standard pipeline run on ``auc_cc``, and that has not happened -- so rather
    than render them from whatever happens to be on disk, they are skipped and their superseded
    PNGs are archived under ``docs/figures/archive/``.
    """
    try:
        figure_data()
    except FileNotFoundError as exc:
        print(f"  {name}: SKIPPED — {str(exc).splitlines()[0]}")
        return True
    return False


def _tidy(ax, *, grid="", ticks=True):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#bdbdb8")
    ax.tick_params(colors=MUTED, labelsize=7.5, length=3)
    if grid:
        ax.grid(axis=grid, color="#ecece8", lw=0.6)
        ax.set_axisbelow(True)
    if not ticks:
        ax.set_xticks([]); ax.set_yticks([])


# ============================================================ 3) the pipeline, drawn
FOLD_COLORS = ["#1f6fb2", "#4b9cd3", "#8fbfe0", "#c3dcef", "#e3eef7"]


def build_pipeline_flow():
    """Eight panels on two rows — a picture per stage, a caption of at most two short lines.

    Stages 5 and 6 call the same drawing functions as ``model_architecture.png`` and
    ``loss_01_objective.png`` with ``compact=True``, so the pipeline cannot drift away from the
    standalone figures. Everything that needs a sentence lives in the docs, not on the slide.
    """
    if _needs_data("pipeline.png") or _example_matches_panel("pipeline.png"):
        return
    d = figure_data()
    corr = panel_corr()
    import pandas as pd

    fig = plt.figure(figsize=(16.0, 8.2))
    bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, 100); bg.set_ylim(0, 100); bg.axis("off")
    bg.text(1.5, 99.0, "OncoTox — the pipeline", ha="left", va="top",
            fontsize=17, fontweight="bold", color=INK)
    bg.text(1.5, 94.0, "one cell in, one AUC per drug out  ·  everything grouped by cell line",
            ha="left", va="top", fontsize=9.5, color=GREY)

    ROW1_TITLE, ROW1_CAP = 86.5, 58.5
    ROW2_TITLE, ROW2_CAP = 46.0, 8.5

    def stage(x, y, n, title, caption, cap_y):
        bg.text(x, y, f"{n}", ha="left", va="bottom", fontsize=12, fontweight="bold", color="#c3c3bd")
        bg.text(x + 2.4, y, title, ha="left", va="bottom", fontsize=10.5, fontweight="bold", color=INK)
        if caption:
            bg.text(x, cap_y, caption, ha="left", va="top", fontsize=8.2, color=MUTED)

    # ================================================== 1 · data
    stage(1.5, ROW1_TITLE, "1", "Data",
          "one AUC per (cell line, drug) — and that one\nvalue labels every cell of the line",
          ROW1_CAP)
    ax = fig.add_axes([0.015, 0.625, 0.20, 0.215]); ax.set_xlim(0, 34); ax.set_ylim(0, 14)
    ax.axis("off")

    rs = np.random.default_rng(3)
    for _ in range(22):                                  # a bounded blob, so nothing leaves the panel
        r, th = 2.0 * np.sqrt(rs.random()), 2 * np.pi * rs.random()
        ax.add_patch(Circle((5.0 + r * np.cos(th) * 1.15, 9.8 + r * np.sin(th)), 0.40,
                            facecolor="#fde0c5", edgecolor="#d2691e", lw=0.7))
    ax.text(5.0, 6.6, "one cell line\n= 56–1,990 cells", ha="center", va="top", fontsize=7.5,
            color=INK)
    arrow(ax, 8.2, 9.8, 11.0, 9.8, color=MUTED)

    x0, y0, cw, ch = 12.2, 3.2, 1.75, 1.15
    sub = d["line_mask"][:6, :10]
    for i in range(6):
        for j in range(10):
            ax.add_patch(Rectangle((x0 + j * cw, y0 + (5 - i) * ch), cw * 0.9, ch * 0.85,
                                   facecolor=BLUE if bool(sub[i, j]) else "#efefeb", lw=0))
    ax.add_patch(Rectangle((x0 - 0.35, y0 + 5 * ch - 0.2), 10 * cw + 0.4, ch * 1.25, fill=False,
                           edgecolor="#d2691e", lw=1.4))
    ax.text(x0 + 10 * cw + 0.9, y0 + 5 * ch + 0.45, "= this line's row", ha="left", va="center",
            fontsize=7.6, color="#d2691e")
    # The grid is a 6x10 corner of `line_mask`, which is (lines x drugs kept by the 50-line cut) --
    # so its width is 534, not CTRPv2's catalogue of 545. Labelled "545 drugs" until 14.08.2026,
    # which described a matrix this figure has never drawn. Taken from the array itself, so the
    # label and the thing it labels are one object. The funnel in stage 2 keeps 545: that tier IS
    # the catalogue, and narrowing away from it is the funnel's whole point.
    ax.text(x0 + 5 * cw, y0 + 6 * ch + 0.5, f"{d['line_mask'].shape[1]} drugs  →", ha="center",
            va="bottom", fontsize=7.6, color=MUTED)
    ax.text(x0 - 0.7, y0 + 3 * ch, "cell lines", ha="center", va="center", fontsize=7.4,
            color=MUTED, rotation=90)
    for dx, (col, lab) in enumerate([(BLUE, "screened"), ("#efefeb", "not screened")]):
        ax.add_patch(Rectangle((12.2 + dx * 9.5, 0.9), 1.4, 0.9, facecolor=col, lw=0))
        ax.text(14.1 + dx * 9.5, 1.35, lab, ha="left", va="center", fontsize=7.4, color=MUTED)

    # ================================================== 2 · drug panel
    stage(26.0, ROW1_TITLE, "2", "Drug panel",
          "only compounds with a published\nsensitivity determinant", ROW1_CAP)
    ax = fig.add_axes([0.255, 0.615, 0.125, 0.225]); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")
    # Both lower tiers are derived. They were hardcoded 173 and 8: the first is a per-drug cell-line
    # count that was never a compound count, the second the void 8-drug panel against a rebuilt one
    # of 11. `n` also sets the tier's drawn width, so each literal made the funnel the wrong *shape*
    # as well as the wrong number.
    for i, (n, lab, col) in enumerate([(545, "545  CTRPv2 compounds", "#c3dcef"),
                                       (N_CANDIDATES,
                                        f"{N_CANDIDATES}  FDA-approved, screened by CTRPv2",
                                        "#6ba7d6"),
                                       (len(PANEL), f"{len(PANEL)}  the panel", BLUE)]):
        half = 4.6 * (0.55 + 0.45 * (n / 545) ** 0.35)
        yb = 7.4 - i * 2.9
        ax.add_patch(mpatches.Polygon(
            [(5 - half, yb + 2.2), (5 + half, yb + 2.2),
             (5 + half * 0.62, yb), (5 - half * 0.62, yb)],
            closed=True, facecolor=col, edgecolor="white", lw=1.2))
        ax.text(5, yb + 1.05, lab, ha="center", va="center", fontsize=7.8,
                color="white" if i == 2 else INK, fontweight="bold" if i == 2 else "normal")
        if i < 2:
            ax.annotate("", xy=(5, yb - 0.15), xytext=(5, yb - 0.6),
                        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.3))

    # ================================================== 3 · split
    n_tr, n_va, n_te = (int(d["n_by_split"][i]) for i in (0, 1, 2))
    n_el = n_tr + n_va
    stage(44.0, ROW1_TITLE, "3", "Split",
          f"the {n_te} test lines are locked away once;\n"
          f"5 folds over the other {n_el}, grouped by line", ROW1_CAP)
    ax = fig.add_axes([0.435, 0.625, 0.225, 0.215]); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.axis("off")

    tot = n_el + n_te
    xw = 96.0 / tot
    xc = 2.0
    for n, col, lab in [(n_el, "#cfe0ef", "train + val"), (n_te, "#f3c9c9", "test")]:
        ax.add_patch(Rectangle((xc, 76), n * xw, 13, facecolor=col, edgecolor="white", lw=1.0))
        ax.text(xc + n * xw / 2, 82.5, lab, ha="center", va="center", fontsize=8, color=INK)
        ax.text(xc + n * xw / 2, 91.5, str(n), ha="center", va="bottom", fontsize=8, color=MUTED)
        xc += n * xw
    ax.text(2.0, 71, f"{tot} cell lines with CTRPv2 labels", ha="left", va="top",
            fontsize=7.6, color=MUTED)
    ax.add_patch(Rectangle((2.0, 76), n_el * xw, 13, fill=False, edgecolor=BLUE, lw=1.6))
    ax.annotate("", xy=(2 + n_el * xw * 0.78, 60), xytext=(2 + n_el * xw * 0.78, 74),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.6))

    order = np.argsort(d["fold"])
    fw = (n_el * xw) / len(order)
    for f in range(1, 6):
        yb = 46 - (f - 1) * 9.5
        held = (d["fold"][order] == f)
        for i, h in enumerate(held):
            ax.add_patch(Rectangle((2.0 + i * fw, yb), fw, 7.4, lw=0,
                                   facecolor=BLUE if h else "#e6e6e2"))
        ax.add_patch(Rectangle((2.0 + n_el * xw + 1.5, yb), n_te * xw, 7.4, lw=0,
                               facecolor="#f6e2e2", hatch="///", edgecolor="#e0b4b4"))
        ax.text(0.5, yb + 3.7, f"fold {f}", ha="right", va="center", fontsize=7.4, color=MUTED)
    for dx, (col, lab) in enumerate([(BLUE, "held out"), ("#e6e6e2", "train")]):
        ax.add_patch(Rectangle((2.0 + dx * 21, -1), 5, 5, facecolor=col, lw=0))
        ax.text(8.0 + dx * 21, 1.5, lab, ha="left", va="center", fontsize=7.4, color=MUTED)
    ax.add_patch(Rectangle((44.0, -1), 5, 5, facecolor="#f6e2e2", hatch="///",
                           edgecolor="#e0b4b4", lw=0))
    ax.text(50.0, 1.5, "test — out of every fold", ha="left", va="center", fontsize=7.4, color=MUTED)

    # ================================================== 4 · representation
    stage(71.0, ROW1_TITLE, "4", "Representation",
          "same trunk, same folds —\nonly the 512-d input differs", ROW1_CAP)
    ax = fig.add_axes([0.715, 0.655, 0.098, 0.165])
    ax.scatter(d["pca2"][:, 0], d["pca2"][:, 1], s=2.2, color="#b9c9d8", lw=0)
    lim = np.abs(d["pca2"]).max() * 1.05
    ax.annotate("", xy=(lim * 0.85, 0), xytext=(-lim * 0.85, 0),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    ax.annotate("", xy=(0, lim * 0.6), xytext=(0, -lim * 0.6),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim * 0.72, lim * 0.72)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_color("#dededa")
    ax.set_title("X_pca", fontsize=9, fontweight="bold", color=GREY, pad=3)
    ax.set_xlabel("linear, unsupervised", fontsize=7.4, color=MUTED, labelpad=2)

    ax = fig.add_axes([0.845, 0.655, 0.098, 0.165]); ax.set_xlim(0, 10); ax.set_ylim(0, 12)
    ax.axis("off")
    ax.set_title("X_scGPT", fontsize=9, fontweight="bold", color=AMBER, pad=3)
    for j in range(6):
        ax.add_patch(Rectangle((1.4 + j * 1.2, 10.2), 0.9, 0.85,
                               facecolor="#f3e2c4", edgecolor=AMBER, lw=0.8))
    ax.text(5.0, 9.7, "the genes of one cell", ha="center", va="top", fontsize=7.2, color=MUTED)
    ax.annotate("", xy=(5, 7.9), xytext=(5, 8.6), arrowprops=dict(arrowstyle="-|>", color=AMBER, lw=1.2))
    for k in range(3):
        ax.add_patch(FancyBboxPatch((1.2, 4.5 + k * 1.15), 7.6, 0.8,
                     boxstyle="round,pad=0.06,rounding_size=0.2",
                     facecolor="#fbf1de", edgecolor=AMBER, lw=1.0))
    ax.text(5.0, 6.05, "transformer (frozen)", ha="center", va="center", fontsize=6.8, color=INK)
    ax.annotate("", xy=(5, 3.2), xytext=(5, 4.1), arrowprops=dict(arrowstyle="-|>", color=AMBER, lw=1.2))
    for j in range(8):
        ax.add_patch(Rectangle((0.7 + j * 1.1, 2.0), 0.95, 0.95,
                               facecolor=plt.colormaps["YlOrBr"](0.25 + 0.55 * ((j * 3) % 7) / 7),
                               edgecolor="white", lw=0.5))
    ax.text(5.0, 1.4, "pretrained on ~33 M cells", ha="center", va="top", fontsize=7.2, color=MUTED)

    # ================================================== 5 · model  (the standalone drawing)
    stage(1.5, ROW2_TITLE, "5", "Model",
          "one cell in, one raw AUC per panel drug out", ROW2_CAP)
    ax = fig.add_axes([0.015, 0.135, 0.375, 0.30])
    draw_architecture(ax, compact=True)

    # ================================================== 6 · loss  (the standalone drawing)
    # Rewritten 14.08.2026 (Selin). The second clause -- "rare response values weighted up" --
    # described alpha > 0, while every panel on this figure is alpha = 0, where W is identically 1
    # and nothing is weighted up at all. Same defect as the architecture subtitle: a swept arm stated
    # as a fact. The caption now separates what the objective always does (masking) from what is
    # swept (the weight), and reads alpha off EXAMPLE_ARM so it tracks the arm the figure plots.
    stage(41.0, ROW2_TITLE, "6", "Loss",
          f"unscreened pairs dropped; the weight W is\n"
          f"a swept arm — α = {EXAMPLE_ARM['alpha']:g} here, so every pair weighs 1", ROW2_CAP)
    ax = fig.add_axes([0.405, 0.215, 0.215, 0.16])
    draw_loss_objective(ax, compact=True)

    # ================================================== 7 · evaluation
    stage(64.5, ROW2_TITLE, "7", "Evaluation",
          "cells → one value per line,\nSpearman within each drug", ROW2_CAP)
    ax = fig.add_axes([0.655, 0.155, 0.125, 0.255])
    # Repointed 14.08.2026 from the void panel, and reduced to ONE seed on purpose: this is a
    # scatter of individual held-out lines, so averaging predictions across seeds would draw a
    # cloud no single run produced. Seed 42 is the first of the three.
    oof = pd.read_csv(PANEL_OUT / "panel_oof_predictions.csv")
    g8 = oof[(oof.rep == "X_scGPT") & (oof.alpha == 0.0) & (oof.loss == FIG_LOSS)
             & (oof.seed == 42) & (oof.drug == "dasatinib")]
    ax.scatter(g8.y_true, g8.y_pred, s=11, color=BLUE, alpha=0.55, edgecolor="white", lw=0.35)
    rho = corr[(corr.rep == "X_scGPT") & (corr.alpha == 0.0)
               & (corr.drug == "dasatinib")].spearman.iloc[0]
    ax.set_xlabel("measured AUC", fontsize=7.2, labelpad=1)
    ax.set_ylabel("predicted", fontsize=7.2, labelpad=1)
    ax.text(0.04, 0.96, f"dasatinib\nρ = {rho:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.6, color=INK)
    _tidy(ax, grid="both")
    ax.tick_params(labelsize=6.5)

    # ================================================== 8 · result
    # Rewritten 14.08.2026 (Selin). Both halves were wrong in different ways. "one seed" was simply
    # false: the bars come from `panel_corr`, which averages seeds 42/43/44 -- stage 7's scatter is
    # one seed on purpose, this is not. And "only scGPT clears the ridge control" asserted a result
    # on a pipeline diagram: true of the plotted bars (scGPT MLP 0.2009 vs its ridge 0.1914,
    # sign-consistent 3/3 seeds; PCA MLP 0.2473 vs its ridge 0.2767) but a +0.0096 margin against a
    # seed sd of 0.0053, and the tallest bar on the panel is the PCA ridge -- so a reader takes
    # "only scGPT clears" for "scGPT wins". The caption now names the axes and the aggregation;
    # the finding lives in docs/steps/05, where it can be sourced and disputed.
    stage(82.0, ROW2_TITLE, "8", "Result",
          f"mean per-drug Spearman over {len(PANEL)} drugs,\n"
          f"{_n_seeds()} seeds — each MLP against its own ridge", ROW2_CAP)
    ax = fig.add_axes([0.868, 0.155, 0.115, 0.255])
    ridge = pd.read_csv(PANEL_OUT / "panel_ridge_baseline.csv")
    bars = [
        ("scGPT MLP", corr[(corr.rep == "X_scGPT") & (corr.alpha == 0.0)].spearman.mean(), AMBER),
        ("PCA MLP", corr[(corr.rep == "X_pca") & (corr.alpha == 0.0)].spearman.mean(), GREY),
        ("scGPT ridge", ridge[ridge.rep == "X_scGPT"].spearman.mean(), "#cfcfc9"),
        ("PCA ridge", ridge[ridge.rep == "X_pca"].spearman.mean(), "#cfcfc9"),
    ]
    y = np.arange(len(bars))[::-1]
    ax.barh(y, [b[1] for b in bars], height=0.6, color=[b[2] for b in bars])
    for yy, (lab, v, _) in zip(y, bars):
        ax.text(v + 0.01, yy, f"{v:.3f}", va="center", ha="left", fontsize=6.8, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([b[0] for b in bars], fontsize=6.8)
    ax.set_xlim(0, 0.50); ax.set_xticks([0, 0.2, 0.4])
    ax.set_xlabel("mean per-drug Spearman", fontsize=7.2, labelpad=1)
    _tidy(ax, grid="x")
    ax.tick_params(labelsize=6.5)

    # ================================================== flow arrows
    for x0_, x1_ in [(22.5, 25.0), (39.0, 41.5), (67.5, 70.0)]:
        arrow(bg, x0_, 73.0, x1_, 73.0, color="#b6b6b0")
    bg.add_patch(FancyArrowPatch((97.0, 66.0), (97.0, 53.0), arrowstyle="-",
                 linewidth=1.8, color="#b6b6b0"))
    bg.add_patch(FancyArrowPatch((97.0, 53.0), (3.0, 53.0), arrowstyle="-",
                 linewidth=1.8, color="#b6b6b0"))
    arrow(bg, 3.0, 53.0, 3.0, 49.0, color="#b6b6b0")
    for x0_, x1_ in [(38.0, 40.5), (60.0, 62.5), (78.5, 81.0)]:
        arrow(bg, x0_, 28.0, x1_, 28.0, color="#b6b6b0")

    out = FIG / "pipeline.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")



# ============================================================ 4) the loss, in three figures

#: Where the generated formula macros land. `report/main.tex` picks them up with `\input`, the same
#: way it already picks up `results_numbers.tex` -- but unlike that file, this one is generated and
#: must never be hand-edited. It is committed so a clean checkout still compiles the report, and a
#: pre-merge check regenerates it and fails on any difference (owned by the gate session, not here).
LOSS_TEX = ROOT / "report" / "loss_objective.tex"

#: The objective, as LaTeX macro bodies -- **the single source for both destinations**. The figure
#: renders these same strings through mathtext, so the equation in `docs/figures/` and the equation
#: in `main.pdf` cannot disagree: there is one place to change and both follow.
#:
#: Notation, and why each symbol is the one it is (settled by Selin 12.08.2026, see the module
#: docstring of `scripts/training/density_weighting.py` for the mechanism):
#:   n in N   -- CELLS. The loss is summed over cells, because training is per cell.
#:   i in I_j -- CELL LINES of the training fold observed for drug j. The density is fitted here,
#:               *not* over cells: a line with 1,990 sequenced cells would otherwise bend the density
#:               toward its own response while a line with 56 barely registered. The two index sets
#:               are shown separately because that asymmetry IS the design decision -- an equation
#:               summing only over cells would be equally true of the naive variant this rejects.
#:   j        -- drug, 1..K.  c -- the weight cap (held fixed, documented as arbitrary).
#:
#: Everything the loss comparison (item 9A) varies is left SYMBOLIC: the pointwise loss `l`, the
#: density exponent `alpha`. Rendering a concrete value here would pre-empt a decision that is
#: Selin's and that the R4 run exists to make.
#:
#: **Huber was dropped from the comparison on 12.08.2026** (`8bda87a`): at `beta=0.05` it sits close
#: to MAE, so the grid would carry two near-duplicate columns, and a principled `beta` would import a
#: new unsourced constant into the very comparison meant to supply the justification. MSE and MAE
#: bracket the axis. `beta` is therefore gone from these macros -- it parameterised only Huber -- and
#: item 9C's `huber_beta` question is retired with it. The code still exposes `--loss huber`; whether
#: the option leaves the code is the model session's call, not this file's.
#:
#: The normalizer is stated as a CONSTRAINT rather than as a constant. `Z_j` is the solution of a
#: fixed point (`fit_weight_fn`, max_iter=50, tol=1e-9) because the clip and the mean-1 condition
#: interact -- writing `/Z_j` would assert a closed form that does not exist, which makes showing it
#: less accurate than omitting it. The `1e-12` density floor is a numerical guard and is deliberately
#: absent; it belongs in the docstring *provided it never binds*, which is measurable on `auc_cc`
#: once R2 produces real values and is not yet established.
#:
#: ⚠️ These strings must parse under BOTH LaTeX and matplotlib's mathtext, which is a subset. Two
#: constructs are therefore avoided, and both were verified failing rather than assumed:
#:   `\tfrac`            -> use `\frac`  (mathtext: "Unknown symbol: \tfrac")
#:   `\lvert` / `\rvert` -> use `\left|` / `\right|`  (mathtext: "Unknown symbol: \lvert")
#: `\text{...}` is supported and is fine. Substituting either back would still compile the report and
#: would break only the figure -- the asymmetry is the reason this warning is here rather than left
#: to be rediscovered.
#:
#: For the same reason there is no `\!` after a *subscripted* symbol (`w_j\left(`, not `w_j\!\left(`).
#: LaTeX tightens the gap; mathtext pulls the parenthesis back over the subscript, so `w_j(y)` renders
#: as an unreadable `w/y`. `\ell\!\left(` is kept -- `\ell` carries no subscript to collide with.
LOSS_TEX_MACROS: dict[str, str] = {
    # The `\;` between the two sums is load-bearing for the figure, not typographic fussiness:
    # mathtext sets both limits *below* their sigma, so without it "n \in \mathcal{N}" and "j = 1"
    # abut and read as one nonsense subscript. In LaTeX it is an ordinary thin space.
    "LossObjective": (
        r"\mathcal{L}\;=\;"
        r"\frac{\sum_{n \in \mathcal{N}}\;\sum_{j=1}^{K} W_{nj}\,"
        r"\ell\!\left(\hat{y}_{nj},\, y_{nj}\right)}"
        r"{\sum_{n \in \mathcal{N}}\;\sum_{j=1}^{K} W_{nj}}"
    ),
    "LossWeightMatrix": r"W_{nj}\;=\;M_{nj}\;w_j\left(y_{nj}\right)",
    "LossWeightFn": (
        r"w_j(y)\;\propto\;\hat{p}_j(y)^{-\alpha},"
        r"\qquad w_j \in \left[\frac{1}{c},\, c\right],"
        r"\qquad \frac{1}{\left| \mathcal{I}_j \right|}"
        r"\sum_{i \in \mathcal{I}_j} w_j\left(y_{ij}\right)\;=\;1"
    ),
    "LossDensity": (
        r"\hat{p}_j\;=\;\text{Gaussian KDE fitted on }"
        r"\left\{\, y_{ij} \;:\; i \in \mathcal{I}_j \,\right\}"
    ),
    "LossArms": (
        r"\ell \in \{\text{squared},\, \text{absolute}\},"
        r"\qquad \alpha \in \{\text{off},\, \frac{1}{2},\, 1\}"
    ),
}


def build_loss_formula_tex():
    """Write the objective's macros to ``report/loss_objective.tex`` for the report to ``\\input``.

    Emits definitions only -- no display environment, no numbering. The section file decides whether
    a macro lands in ``equation``, ``align`` or inline, so the generator never dictates typesetting
    it cannot see the context for.
    """
    lines = [
        "% " + "=" * 76,
        "% loss_objective.tex -- the training objective, as LaTeX macros.",
        "% " + "=" * 76,
        "%",
        "% GENERATED FILE -- DO NOT EDIT.",
        "%   Written by docs/make_figures.py :: build_loss_formula_tex() from LOSS_TEX_MACROS,",
        "%   which is also what docs/figures/loss_01_objective.png renders -- so that the equation",
        "%   in the report and the equation in the figure would have one source and could not",
        "%   drift apart. To change a formula, edit LOSS_TEX_MACROS and re-run the generator.",
        "%",
        "% HOW THEY ARE USED",
        "%   INTENDED USE: main.tex does \\input{loss_objective} in the preamble; a section then",
        "%   writes e.g. \\begin{equation}\\LossObjective\\end{equation}. Definitions only -- the",
        "%   section chooses the environment.",
        "%",
        "%   WIRED IN 14.08.2026: main.tex \\input{loss_objective}s in the preamble and 03_methods",
        "%   sets \\LossObjective and \\LossWeightMatrix as equations (eq:objective, eq:weight).",
        "%   \\LossWeightFn / \\LossDensity are deliberately NOT cited: they describe the density",
        "%   weighting, which is off (alpha = 0) in every arm the report reports.",
        "%",
        "% WHAT IS DELIBERATELY SYMBOLIC",
        "%   \\ell and \\alpha are left as symbols: the loss comparison (docs/TODO.md item 9A)",
        "%   sweeps MSE / MAE against alpha in {off, 0.5, 1.0} -- six arms. A concrete value here",
        "%   would pre-empt a decision the R4 run exists to make. Huber was dropped from the",
        "%   comparison on 12.08.2026, and \\beta went with it: it parameterised only Huber.",
        "%",
        "% NOTATION",
        "%   n in N    cells -- the loss is summed over cells, training being per cell.",
        "%   i in I_j  cell lines of the training fold observed for drug j -- where the density is",
        "%             fitted. Kept a separate index on purpose: fitting on cells would let a line",
        "%             with 1,990 sequenced cells bend the density and one with 56 barely register.",
        "%   j         drug, 1..K.        c   the weight cap (fixed; arbitrary, and documented so).",
        "",
    ]
    lines += [f"\\newcommand{{\\{name}}}{{{body}}}" for name, body in LOSS_TEX_MACROS.items()]
    lines.append("")

    LOSS_TEX.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {LOSS_TEX}")


def draw_loss_objective(ax, *, compact: bool = False):
    """Draw the objective into ``ax`` — same drawing for the standalone figure and for stage 6.

    Every equation here is rendered from :data:`LOSS_TEX_MACROS`, the same strings
    :func:`build_loss_formula_tex` writes into the report. Nothing is retyped, so the figure and
    ``main.pdf`` cannot disagree about what is being minimised.

    Until 12.08.2026 they *did* disagree, and only the figure was wrong. It hardcoded its own
    ``(\\hat{y}-y)^2``, titled the objective "a weighted, masked mean squared error", and labelled a
    single index as ``i = cell`` — asserting a fixed squared loss when the loss comparison
    (``docs/TODO.md`` item 9A) sweeps MSE / MAE, and hiding the cell/cell-line split that is
    the whole design of the weighting. It went stale silently: it reads no data, so it rendered
    happily every run with nothing to fail.
    """
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    tex = LOSS_TEX_MACROS

    if not compact:
        ax.text(1, 99, "The objective — a weighted, masked mean error",
                ha="left", va="top", fontsize=14, fontweight="bold", color=INK)
        ax.text(1, 92, "target = auc_cc, the curve-fit AUC  ·  per-drug scaling belongs here, in the "
                       "loss, rather than in the labels",
                ha="left", va="top", fontsize=9, color=GREY)
        # The three swept quantities are named as swept, so the figure cannot be read as claiming a
        # loss the run has not chosen yet.
        ax.text(1, 86, r"$\ell$ and $\alpha$ are what the loss comparison sweeps "
                       "(item 9A) — the figure fixes neither",
                ha="left", va="top", fontsize=9, color=GREY)

    ax.text(50, 80 if not compact else 92, f"${tex['LossObjective']}$",
            ha="center", va="top", fontsize=20 if not compact else 15, color=INK)
    if compact:
        return

    # The weight definition is not decoration: the objective above sums over cells (n), while the
    # density is fitted on the training fold's cell LINES (i). Dropping these two lines would leave
    # an equation equally true of the naive per-cell variant this deliberately rejects.
    ax.text(50, 47, f"${tex['LossWeightMatrix']}$",
            ha="center", va="top", fontsize=13, color=INK)
    ax.text(50, 39, f"${tex['LossWeightFn']}$",
            ha="center", va="top", fontsize=13, color=INK)

    for x, edge, fill, head, body in [
        (1.5, GREY, GREY_FILL, r"$M_{nj}$   mask",
         "1 if the line was screened\nagainst drug $j$, else 0"),
        (35.0, BLUE, BLUE_FILL, r"$w_j(y_{nj})$   sample weight",
         "inverse label density — fitted\non cell lines, mean 1"),
        (68.5, RED, RED_FILL, r"$\ell(\hat{y}_{nj},\,y_{nj})$   error",
         "squared or absolute\n$n$ = cell,   $i$ = cell line,   $j$ = drug"),
    ]:
        ax.add_patch(FancyBboxPatch((x, 1), 30, 23, boxstyle="round,pad=0.5,rounding_size=2.0",
                     linewidth=1.8, edgecolor=edge, facecolor=fill, zorder=2))
        ax.text(x + 15, 21.0, head, ha="center", va="top", fontsize=11.5,
                fontweight="bold", color=edge, zorder=3)
        ax.text(x + 15, 13.0, body, ha="center", va="top", fontsize=8.5, color=INK,
                zorder=3, linespacing=1.5)


def build_loss_objective():
    """What the objective is made of — the formula and the three factors, nothing else."""
    fig, ax = plt.subplots(figsize=(12.0, 5.6))
    draw_loss_objective(ax)
    out = FIG / "loss_01_objective.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def draw_weight_curve(ax, *, drug: str = "dasatinib", compact: bool = False):
    """The fitted weight curve of one drug — stage 6 of the pipeline and the lower half of loss_02."""
    import sys

    sys.path.insert(0, str(ROOT))
    from scripts.training.density_weighting import DEFAULT_CAP, fit_weight_fn

    vals = figure_data()["panel_auc"][:, PANEL.index(drug)]
    vals = vals[np.isfinite(vals)]
    fn = fit_weight_fn(vals)
    grid = np.linspace(vals.min() - 0.03, 1.1, 300)
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
    ax.axhline(DEFAULT_CAP, color=RED, lw=0.9, ls="--")
    ax.plot(grid, fn(grid), color=BLUE, lw=2.2)
    ax.plot(vals, fn(vals), "o", ms=2.6, color=BLUE, alpha=0.3, markeredgecolor="white",
            markeredgewidth=0.4)
    ax.set_ylim(0, DEFAULT_CAP * 1.15)
    fs = 6.6 if compact else 7.5
    ax.set_xlabel("AUC", fontsize=fs, labelpad=1)
    ax.set_ylabel("weight" if compact else "sample weight", fontsize=fs, labelpad=1)
    if compact:   # keep the label off the cap line, which the curve saturates at both ends
        ax.text(0.98, 0.04, f"{drug} · cap {DEFAULT_CAP:g}×", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=fs, color=MUTED)
    else:
        ax.text(0.98, 0.94, f"{drug}\ncap {DEFAULT_CAP:g}×", transform=ax.transAxes,
                ha="right", va="top", fontsize=fs, color=MUTED)
    _tidy(ax)
    if compact:
        ax.tick_params(labelsize=6.2)
        ax.set_yticks([0, 1, 3])
    return fn, vals


def build_loss_weights():
    """Where the weights come from: one drug's label density, and the curve it produces."""
    if _needs_data("loss_02_weights.png"):
        return
    import sys

    sys.path.insert(0, str(ROOT))
    from scripts.training.density_weighting import DEFAULT_ALPHA, DEFAULT_CAP, fit_weight_fn

    drug = "dasatinib"
    vals = figure_data()["panel_auc"][:, PANEL.index(drug)]
    vals = vals[np.isfinite(vals)]
    fn = fit_weight_fn(vals)
    grid = np.linspace(vals.min() - 0.03, 1.1, 400)
    dens, wcurve = fn.kde(grid), fn(grid)

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.0, 6.8), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 1.15], hspace=0.12))
    fig.subplots_adjust(top=0.845)
    fig.suptitle("Rare response values get more weight", x=0.02, ha="left",
                 fontsize=13, fontweight="bold", color=INK, y=0.995)
    fig.text(0.02, 0.945, f"{drug}, {len(vals)} cell lines — the density is fitted per drug, "
                          "inside each fold, on training lines only",
             ha="left", va="top", fontsize=8.8, color=GREY)
    fig.text(0.02, 0.905, r"$w(y)\;\propto\;\hat{p}(y)^{-\alpha}$,   "
                          rf"$\alpha={DEFAULT_ALPHA:g}$,   clipped to "
                          rf"$[1/{DEFAULT_CAP:g},\,{DEFAULT_CAP:g}]$,   normalized to mean 1",
             ha="left", va="top", fontsize=8.8, color=INK)

    top = float(dens.max())
    a1.hist(vals, bins=22, color="#dededa", edgecolor="white", lw=0.7, density=True)
    a1.plot(grid, dens, color=INK, lw=2.0)
    a1.set_ylim(0, top * 1.55)
    a1.set_ylabel("label density", fontsize=9)
    a1.annotate("crowded middle", xy=(float(grid[int(np.argmax(dens))]), top * 1.02),
                xytext=(0.66, top * 1.42), fontsize=9, color=MUTED, ha="left", va="top",
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
    a1.annotate("sparse extremes", xy=(0.20, 0.10), xytext=(0.06, top * 1.42),
                fontsize=9, color=MUTED, ha="left", va="top",
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))

    a2.axhline(1.0, color=MUTED, lw=0.9, ls=":")
    for lvl in (DEFAULT_CAP, 1.0 / DEFAULT_CAP):
        a2.axhline(lvl, color=RED, lw=1.0, ls="--")
    a2.plot(grid, wcurve, color=BLUE, lw=2.4, zorder=3)
    a2.plot(vals, fn(vals), "o", ms=4.5, color=BLUE, alpha=0.35, markeredgecolor="white",
            markeredgewidth=0.5, zorder=4)
    a2.set_ylim(0, DEFAULT_CAP * 1.2)
    a2.set_ylabel("sample weight  $w(y)$", fontsize=9)
    a2.set_xlabel("AUC   (winsorized at 1.1)", fontsize=9)
    a2.text(0.995, 1.06, "unweighted (w = 1)", transform=a2.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=8, color=MUTED)
    a2.text(0.5, DEFAULT_CAP + 0.06, f"cap {DEFAULT_CAP:g}×", transform=a2.get_yaxis_transform(),
            ha="center", va="bottom", fontsize=8, color=RED)
    a2.text(0.02, 1 / DEFAULT_CAP + 0.06, f"floor 1/{DEFAULT_CAP:g}",
            transform=a2.get_yaxis_transform(), ha="left", va="bottom", fontsize=8, color=RED)
    for a in (a1, a2):
        _tidy(a)
        a.tick_params(labelsize=8)

    out = FIG / "loss_02_weights.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def build_loss_effect():
    """What the weighting did: the spread it was designed to raise, and the ranking it did not move.

    ⚠️ **The drawing below is superseded and must not be rebuilt as it stands.** It is a paired
    scatter — unweighted on x, density-weighted on y, one point per drug, with a diagonal — and that
    form encodes an assumption that has stopped holding: that the comparison has exactly **two**
    arms. From R4 the density-weighting arm is three levels, ``alpha`` in {0.0, 0.5, 1.0}, replacing
    the boolean ``weighted``. There is no single delta to put on two axes.

    **DESIGN DECIDED 12.08.2026 by the visualization agent** — figure decisions were delegated by
    Selin on that date; she keeps the analysis decisions. Recorded here rather than in a message so
    it cannot drift from the code it governs. *Not implemented*: the figure is skipped until R4 (its
    input is a training output), and this project does not commit a drawing nobody has rendered and
    looked at. Implement it when the R4 outputs exist, then render it and read it.

    **Replace the paired scatter with a response curve over the sweep — put ``alpha`` on the x-axis.**

      * A **2x2 grid**: rows are the two quantities, spread of the predictions (``pred_std``) and
        per-drug Spearman; columns are the two representations, scGPT and PCA. Faceting the
        representation rather than colouring it is what buys room for the seed band below without
        four series colliding in one panel.
      * Rows are the two quantities and not one, because that pairing *is* audit 09's finding — the
        weighting acts on spread, and the previous null was read off rank correlation and error,
        both structurally blind to it. Keeping both is what lets the re-test distinguish "no effect"
        from "an effect the metrics cannot see".
      * x = ``alpha`` as ordered ticks; y = the quantity.
      * One faint line per drug across the ticks, at that drug's median over seeds, so *within-drug*
        movement stays visible — the evaluation metric is within-drug, and a mean over drugs would
        hide a reversal.
      * One bold line: the median across drugs. **Around it, a band spanning the seeds.**
      * A horizontal reference line for the ridge baseline, which has no ``alpha``.

    **The seed band is not decoration and must not be dropped for tidiness.** The decision rule's
    margin is *defined* as the spread across seeds (>=3, item 9A), so a figure without it invites
    exactly the reading the rule exists to prevent: an alpha-to-alpha difference smaller than the
    noise, read off the plot as an effect. Whatever else is cut, this stays.

    ``loss`` (the MSE / MAE arm of item 9A) is **held fixed** here and named in the caption. This
    figure answers "what does ``alpha`` do"; crossing it with the loss arm as well would put four
    dimensions in one drawing to answer a one-dimensional question. Which loss it is fixed at is an
    analysis decision and therefore Selin's, not a display choice.

    **Why this form and not the alternatives.**

      * *0.5 and 1.0 each against 0.0*, or *best-alpha against 0.0*, keep the scatter but need a new
        panel per pair and grow as O(N) panels — and "best-alpha" hides the axis whose legibility is
        the entire reason the sweep exists.
      * A curve **scales to any number of arms**: adding a level adds a tick, not a panel. The arm
        list is explicitly provisional and expected to grow, so a design that only works for exactly
        three points would need redoing on the next change.
      * It **puts the ridge baseline where it belongs**. Selin resolved the schema question on
        12.08.2026: a separate ``model`` column in {mlp, ridge, and mil when 4b lands}, with
        ``alpha`` numeric and meaningful for ``model='mlp'`` rows only. Ridge therefore has no
        ``alpha`` and is drawn as a horizontal reference line — not a workaround for a sentinel, but
        the honest picture of a baseline that does not vary along this axis. The design was chosen
        before that decision landed and needed no change when it did.
      * It **degenerates gracefully**: at two ticks it is a slopegraph, so this is a generalisation
        of the current figure rather than a different figure.

    What is lost, stated plainly: the diagonal, and with it the "above the line: the model hedges
    less" annotation. Both only work for exactly two arms. The replacement reading is that a **flat**
    line means ``alpha`` does nothing, and a **rising** line in the spread panel means the weighting
    reduced the shrinkage it was designed to reduce.

    The same decision settles ``4a``'s summary rows: one row per level (``MLP alpha=0.0`` /
    ``0.5`` / ``1.0``) instead of ``MLP weighted`` / ``MLP unweighted``, with ridge as its own row
    rather than a value in the column.

    ⚠️ **The title is a conclusion and must be rewritten after the run, not before.** "The weighting
    fired — and the ranking did not move" is the 27.07 result, on the void panel and the retired
    target, and it is exactly the claim R4 re-tests. Carrying it into a figure drawn on new outputs
    would assert the answer the run exists to find.
    """
    corr = panel_corr()

    # This figure is the one thing here that cannot follow the panel. Its input is a TRAINING
    # output -- per-drug correlations from the 27.07 run on the void 8-drug panel -- so unlike the
    # label-derived panels it cannot be rebuilt by re-reading the h5ad; it needs a re-run, which the
    # freeze holds until R4. With the rebuilt panel only 3 of 11 drugs have a row, so it would
    # quietly plot three points and read as a complete comparison. Skipped instead, leaving the
    # committed PNG (captioned as void) in place.
    have = set(corr["drug"])
    missing = [d for d in PANEL if d not in have]
    if missing:
        print(f"  loss_03_effect.png: SKIPPED -- panel_per_drug_correlation.csv has no rows for "
              f"{len(missing)} of {len(PANEL)} panel drugs ({', '.join(missing[:4])}"
              f"{'...' if len(missing) > 4 else ''}). It is a training output; re-runs at R4.")
        return

    reps = [("X_scGPT", AMBER, "scGPT"), ("X_pca", GREY, "PCA")]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 5.4))
    fig.subplots_adjust(top=0.82)
    fig.suptitle("Prediction spread and per-drug Spearman, unweighted against density-weighted",
                 x=0.02, ha="left", fontsize=13, fontweight="bold", color=INK, y=0.99)
    fig.text(0.02, 0.925, f"one point per drug, out of fold, mean of three seeds, loss={FIG_LOSS} — "
                          f"alpha=0 against alpha={FIG_ALPHA_ON:g}",
             ha="left", va="top", fontsize=8.8, color=GREY)

    for ax, col, label, lim in [(a1, "pred_std", "spread of the predictions", (0.03, 0.11)),
                                (a2, "spearman", "per-drug Spearman", (0.0, 0.65))]:
        ax.plot(lim, lim, color="#c9c9c4", lw=1.0, ls="--", zorder=1)
        for rep, c, name in reps:
            u = corr[(corr.rep == rep) & (corr.alpha == 0.0)].set_index("drug").reindex(PANEL)[col]
            w = corr[(corr.rep == rep) & (corr.alpha == FIG_ALPHA_ON)].set_index("drug").reindex(PANEL)[col]
            ax.scatter(u, w, s=46, color=c, alpha=0.85, edgecolor="white", lw=0.8,
                       label=name, zorder=3)
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.set_xlabel(f"unweighted — {label}", fontsize=9)
        ax.set_ylabel(f"density-weighted — {label}", fontsize=9)
        ax.set_aspect("equal")
        _tidy(ax, grid="both")
        ax.tick_params(labelsize=8)
    a1.legend(frameon=False, fontsize=8.5, loc="lower right")
    # Were "above the line: the model hedges less" and "on the line: no gain, no loss". The second
    # was a reading of the plot AND is now false: at alpha=0.5 the per-drug Spearman moves +0.0281
    # on X_pca and -0.0082 on X_scGPT (panel_metrics.csv). The diagonal needs identifying, not
    # interpreting; what a departure from it means belongs in prose that can be cited.
    for _ax in (a1, a2):
        _ax.text(0.04, 0.96, "dashed: y = x", transform=_ax.transAxes,
                 ha="left", va="top", fontsize=8.5, color=MUTED)

    out = FIG / "loss_03_effect.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")



def build_q2_instrument():
    """Q2's two halves in one figure: the structure reproduces, the instrument is near chance.

    **Why this figure exists.** Until 14.08.2026 Q2 had no figure at all -- eleven CSVs and no
    image -- while being one of the two questions the project is organised around.

    **What is on the axes, and why these two stages of the seven.** Q2's verdict is ``POSITIVE`` and
    its weakness is measured, and both facts live in different stage tables, so quoting either alone
    misrepresents it. Stage 2 asks whether the within-line structure the model imposes *reproduces
    across seeds* -- if it did not, there would be nothing to discuss. Stage 7 is the **positive
    control**: it asks whether the instrument can detect a between-line gap it already knows is
    there. Putting them side by side is the honest statement, because the second bounds what the
    first is worth.

    Distributions rather than the medians alone: the medians are in
    ``notebooks/outputs/mil/q2_verdict.csv`` and a bar chart of four numbers would hide that stage
    7's mass sits against its null.

    The title says what is plotted. What it means belongs in prose that can be cited and disputed --
    this file's own rule against asserting results inside a drawing.
    """
    import pandas as pd

    s2 = pd.read_csv(MIL_OUT / "stage2_cross_seed_agreement.csv")
    s7 = pd.read_csv(MIL_OUT / "stage7_positive_control.csv")

    reps = [("X_pca", BLUE, "PCA"), ("X_scGPT", AMBER, "scGPT")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    fig.subplots_adjust(top=0.80, wspace=0.28)
    fig.suptitle("Q2 instrument: cross-seed agreement, and the positive control",
                 x=0.02, ha="left", fontsize=13, fontweight="bold", color=INK, y=0.99)
    fig.text(0.02, 0.90,
             "left: stage 2, one point per (drug, cell line, seed pair).   "
             "right: stage 7, one point per (drug, line pair, seed).   "
             "dashed = the null each is measured against",
             ha="left", va="top", fontsize=8.6, color=GREY)

    for rep, c, name in reps:
        a1.hist(s2[s2.rep == rep]["rho"].dropna(), bins=60, range=(-1, 1), histtype="step",
                lw=1.8, color=c, label=name, density=True)
    a1.axvline(0.0, color="#9a9a95", lw=1.1, ls="--", zorder=1)
    a1.set_xlabel("cross-seed Spearman of per-cell predictions", fontsize=9)
    a1.set_ylabel("density", fontsize=9)
    a1.set_xlim(-1, 1)

    for rep, c, name in reps:
        a2.hist(s7[s7.rep == rep]["auroc"].dropna(), bins=60, range=(0, 1), histtype="step",
                lw=1.8, color=c, label=name, density=True)
    a2.axvline(0.5, color="#9a9a95", lw=1.1, ls="--", zorder=1)
    a2.set_xlabel("within-bag AUROC, known between-line gap", fontsize=9)
    a2.set_ylabel("density", fontsize=9)
    a2.set_xlim(0, 1)

    # Medians read from the same rows that are plotted, so the annotation cannot drift from the bars.
    for ax, df, col in ((a1, s2, "rho"), (a2, s7, "auroc")):
        lines = []
        for rep, _, name in reps:
            med = df[df.rep == rep][col].median()
            lines.append(f"{name} median {med:.3f}")
        # Labelled "pooled" on purpose. q2_verdict.csv's headline medians aggregate PER SEED
        # (and, for stage 7, take the median of per-seed medians), so a pooled median over every
        # point drawn here differs in the third decimal. Saying which statistic this is keeps the
        # figure self-consistent with what it plots instead of appearing to contradict the verdict.
        ax.text(0.03, 0.97, "pooled median of points shown\n" + "\n".join(lines),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=8.5, color=MUTED, linespacing=1.5)
        ax.legend(frameon=False, fontsize=8.5, loc="upper right")
        _tidy(ax, grid="y")
        ax.tick_params(labelsize=8)

    out = FIG / "q2_instrument.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def build_panel_response():
    """The target itself, per panel compound: how `auc_cc` is distributed across the cell lines.

    **Why it exists.** The panel slide named eleven compounds and showed nothing about them. What a
    viewer needs in order to read every later result is *what there is to predict for each* — the
    label's location and its spread. A compound whose lines all sit at the same response cannot be
    ranked by any model, and one of the eleven is very close to that.

    **Ordered by interquartile range**, widest at the top. The ordering is a display choice and is
    stated in the subtitle; nothing downstream depends on it.

    ⚠️ **This figure must not be read as a selection criterion, and the panel is its own proof.**
    Response spread is a statistic of *our own labels*; selecting on it is what voided two earlier
    panels ([Step 05](../docs/steps/05-multitask-results.md) and the corrections file). The panel was
    chosen on approval status, a published resistance determinant, and screen coverage — never on
    this. **The compound with the narrowest spread of the eleven is in the panel**, which is what
    selecting on spread would have prevented.

    Population: the train+val cell lines, i.e. the ones every fold is fitted and scored on, not all
    181 in the overlap. Read from ``figure_data()['panel_auc']``, which is built from the ``auc_cc``
    targets h5ad.

    The title names what is plotted. What it means is prose, per this file's own rule.
    """
    d = figure_data()
    auc, panel = d["panel_auc"], _panel()

    order = sorted(range(len(panel)),
                   key=lambda k: np.nanpercentile(auc[:, k], 75) - np.nanpercentile(auc[:, k], 25))
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    fig.subplots_adjust(top=0.86, left=0.24)
    fig.suptitle("Response (auc_cc) across cell lines, per panel compound",
                 x=0.02, ha="left", fontsize=13, fontweight="bold", color=INK, y=1.0)
    # One line, and only what cannot be inferred from the plot: n, the whisker convention (5-95 is
    # NOT matplotlib's default 1.5xIQR) and the row ordering. Direction and the meaning of 1.0 were
    # here too and are gone -- this is an expert audience and the axis carries them.
    #
    # Stated as PERCENTILES rather than as "quartiles / whiskers / IQR" (Selin, 15.08.2026, who
    # asked what they meant). The percentile form is no longer and strictly more precise -- the
    # jargon does not say *which* whisker convention, which is the one thing here that is
    # non-default -- and the presenter should not have to define a term to answer a question.
    fig.text(0.02, 0.905,
             f"{int(np.isfinite(auc).any(1).sum())} train+val cell lines  ·  "
             "box = 25th–75th percentile, line = median, whiskers = 5th–95th  ·  "
             "ordered by box width",
             ha="left", va="top", fontsize=8.6, color=GREY)

    rng = np.random.default_rng(0)                # jitter only, never the values
    for row, k in enumerate(order):
        v = auc[:, k][np.isfinite(auc[:, k])]
        ax.scatter(v, np.full(v.size, row) + rng.uniform(-0.17, 0.17, v.size),
                   s=5, color=GREY, alpha=0.35, linewidths=0, zorder=2)
        ax.boxplot(v, positions=[row], vert=False, widths=0.62, whis=(5, 95), showfliers=False,
                   patch_artist=True, zorder=3,
                   boxprops=dict(facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.2),
                   medianprops=dict(color=INK, lw=1.6),
                   whiskerprops=dict(color=BLUE, lw=1.1), capprops=dict(color=BLUE, lw=1.1))

    ax.axvline(1.0, color=RED, lw=1.2, ls="--", zorder=4)
    ax.text(1.01, len(panel) - 0.4, "no effect", fontsize=8.4, color=RED, va="center")

    names = _panel_display()
    ax.set_yticks(range(len(panel)))
    ax.set_yticklabels([names[panel[k]] for k in order], fontsize=9.2)
    ax.set_xlabel("auc_cc", fontsize=9.5)
    ax.set_ylim(-0.7, len(panel) - 0.3)
    _tidy(ax, grid="x")
    ax.tick_params(labelsize=8.5)

    out = FIG / "panel_response.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def build_lpo_bias():
    """Where the within-drug ranking comes from: lineage, cell-line identity, or the embeddings.

    **The question.** The project's motivating hypothesis is that a model can score well by
    recognising tissue of origin rather than by learning drug response. Under the project's own
    leave-**cell-line**-out protocol that cannot be measured at all -- a held-out line is unseen, so a
    line-mean baseline emits a constant and scores exactly 0.0000. Under leave-**pairs**-out the
    held-out pairs come from lines that are in training, so each bias channel can be run as an honest
    predictor and read off directly.

    **Why these six bars and not the ten the run produced.** Four of the ten are unplottable or
    uninformative here. ``NaivePredictor`` and ``NaiveDrugMeanPredictor`` are **constant within a
    drug by definition**, so a within-drug rank correlation is undefined for them -- an absence by
    construction, not a result. ``NaiveMeanEffectsPredictor`` scores **identically** to
    ``NaiveCellLineMeanPredictor`` (0.3220 both), and must: within one drug the drug effect is an
    additive constant, so the two induce the same ranking of cell lines. Drawing it again would
    suggest two measurements where there is one. ``SingleDrugElasticNet (scgpt)`` **is** drawn, as a
    zero-height bar labelled *collapsed*, because unlike the first two it is a fitted model that
    found no within-drug signal and shrank to the intercept -- that is a result and hiding it would
    flatter the scGPT column.

    **The axis is this project's own metric**, mean Spearman within each drug across held-out cell
    lines, so the bars are commensurable with the numbers reported everywhere else -- but they are
    **LPO and the project's headline is LCO**, which is why no OncoTox bar appears here and why the
    subtitle says so. The model has never been run under this protocol; inventing a bar for it by
    quoting its LCO score is exactly the comparison the record forbids.

    ⚠️ **The two random-forest bars do not reproduce across executions** and are hatched to say so.
    ``drevalpy``'s default hyperparameter set pins no ``random_state``, so re-running the script moved
    ``SingleDrugRandomForest (scgpt)`` 0.2401 -> 0.2141 and ``(pca)`` 0.0466 -> 0.0271 while every
    other bar was bit-identical. **No claim the figure supports rests on them:** both sit below the
    line-mean baseline in either execution.

    The title states what is plotted. What it means is prose, per this file's own rule.
    """
    import pandas as pd

    src = ROOT / "notebooks" / "outputs" / "dreval" / "dreval_lpo_results.csv"
    d = pd.read_csv(src)
    g = d.groupby("algorithm").agg(m=("per_drug_Spearman", "mean"),
                                   s=("per_drug_Spearman", "std"),
                                   n=("n_scored", "sum"))

    # Every bar is scored on the same pairs, or the comparison is not one. Checked rather than
    # assumed: the per-drug models are scored only where prediction succeeded, so this could
    # silently have been false.
    if g["n"].nunique() != 1:
        raise SystemExit(f"algorithms scored on different numbers of pairs: {g['n'].to_dict()}")

    #: (key, label, group, colour, hatched-because-nondeterministic)
    BARS = [
        ("NaiveTissueMeanPredictor", "lineage only", 0, GREY, False),
        ("NaiveCellLineMeanPredictor", "cell-line identity", 0, INK, False),
        ("SingleDrugElasticNet (pca)", "elastic net · PCA", 1, BLUE, False),
        ("SingleDrugRandomForest (pca)", "random forest · PCA", 1, BLUE, True),
        ("SingleDrugElasticNet (scgpt)", "elastic net · scGPT", 1, AMBER, False),
        ("SingleDrugRandomForest (scgpt)", "random forest · scGPT", 1, AMBER, True),
    ]

    #: Group headers occupy their own row rather than a footnote, so the two blocks are named where
    #: they are read. A caption strip under the axes collided with the last bar.
    HEADERS = {0: "NO EXPRESSION FEATURES", 1: "OUR EMBEDDINGS · DrEval's per-drug models"}

    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    fig.subplots_adjust(top=0.84, left=0.28)
    # Title names what is plotted, per rule 1. It read "Where the within-drug ranking of cell lines
    # comes from" in the first draft -- a conclusion-shaped title, and exactly what rule 1 forbids.
    fig.suptitle("Mean per-drug Spearman under DrEval leave-pairs-out",
                 x=0.02, ha="left", fontsize=13, fontweight="bold", color=INK, y=1.0)
    # Trimmed to one line on the same rule as build_panel_response(): only what cannot be read off
    # the plot. The long "held-out pairs come from lines that are in training" clause is gone --
    # "leave-pairs-out" already says it to anyone who knows the protocol, and the warning that
    # matters (not comparable with LCO) is kept as a single clause rather than a sentence.
    fig.text(0.02, 0.905,
             f"5 folds, {int(g['n'].iloc[0] / 5):,} held-out pairs per fold  ·  "
             "whiskers = sd over folds  ·  hatched = unseeded, does not reproduce  ·  "
             "not comparable with leave-cell-line-out numbers",
             ha="left", va="top", fontsize=8.6, color=GREY)

    ref = float(g.loc["NaiveCellLineMeanPredictor", "m"])
    xmax = max(0.46, float((g["m"] + g["s"]).max()) + 0.055)
    label_x = xmax - 0.004                       # one right-aligned column, so values cannot drift

    rows, y = [], 0.0                            # laid out top-down, then flipped by set_ylim
    for i, (key, label, grp, colour, hatch) in enumerate(BARS):
        if i == 0 or BARS[i - 1][2] != grp:      # a header row opening each group
            ax.text(0.003, y, HEADERS[grp], va="center", ha="left", fontsize=8.2,
                    color=MUTED, fontweight="bold", zorder=5)
            y -= 1.0
        rows.append((y, label))
        m, s = g.loc[key, "m"], g.loc[key, "s"]
        if pd.isna(m):                           # fitted, then collapsed to the intercept
            # Descriptive, not interpretive (rule 2): it read "no within-drug signal", which is the
            # reading. What is observable is that the fitted model emits one value per drug.
            ax.text(0.003, y, "shrank to the intercept — constant within each drug",
                    va="center", ha="left", fontsize=8.4, color=AMBER, style="italic", zorder=5)
            ax.text(label_x, y, "n/a", va="center", ha="right", fontsize=9,
                    color=AMBER, fontweight="bold", zorder=5)
        else:
            ax.barh(y, m, height=0.62, color=colour, alpha=0.32 if hatch else 0.88,
                    edgecolor=colour, lw=1.2, hatch="///" if hatch else None, zorder=3)
            ax.errorbar(m, y, xerr=s, fmt="none", ecolor=INK, elinewidth=1.1, capsize=3, zorder=4)
            ax.text(label_x, y, f"{m:.3f}", va="center", ha="right", fontsize=9,
                    color=INK, fontweight="bold", zorder=5)
        y -= 1.0

    ax.axvline(ref, color=RED, lw=1.2, ls="--", zorder=2)
    ax.text(ref - 0.006, 0.62, "knowing only which cell line it is",
            fontsize=8.4, color=RED, va="bottom", ha="right")

    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=9.2)
    ax.set_xlabel("mean per-drug Spearman  (this project's metric)", fontsize=9.5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(y + 0.45, 1.15)
    _tidy(ax, grid="x")
    ax.tick_params(labelsize=8.5)

    out = FIG / "lpo_bias.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    build_pipeline()
    build_pipeline_flow()
    build_architecture()
    build_loss_objective()
    # Written from the same macros build_loss_objective() just rendered, and immediately after it, so
    # the report's equations and the figure's cannot be regenerated apart.
    build_loss_formula_tex()
    build_loss_weights()
    build_loss_effect()
    build_q2_instrument()
    build_panel_response()
    build_lpo_bias()
