"""Generate the OncoTox slide/doc graphics into ``docs/figures/``.

**Pure drawings — always reproducible, always built:**

  pipeline_overview.png    status of the whole project against the plan (steps 01-08)
  loss_01_objective.png    what the objective is made of

**And one artifact that is not a figure**, written to ``report/`` rather than ``docs/figures/``:

  report/loss_objective.tex   the objective's equations, as LaTeX macros the report ``\\input``s

It is generated from ``LOSS_TEX_MACROS`` — the same strings ``loss_01_objective.png`` renders — so
the maths in the report and the maths in the figure have one source and cannot drift apart. Decided
by Selin 12.08.2026 over rendering the formula into the PNG (raster, fonts would not match the
report body, not referenceable by LaTeX) and over a PGF/vector figure (fonts match, but it puts a
LaTeX installation in the figure build path). It is committed so a clean checkout still compiles
the report; a regenerate-and-diff pre-merge check that fails on drift is owned by the gate session.

**Derived from data — currently SKIPPED, and archived (12.08.2026):**

  pipeline.png             the pipeline as a picture, stage by stage
  model_architecture.png   one cell in, one AUC per panel drug out
  loss_02_weights.png      one drug's label density and the weight curve fitted to it
  loss_03_effect.png       what the weighting did, per drug: spread up, ranking flat

The four read either ``figure_data.npz`` -- rebuilt from the ``auc_cc`` targets h5ad -- or a
training output under ``notebooks/outputs/panel/``. Neither exists on the current target:
preprocessing has not re-run under the freeze, and the last training run was on the void 8-drug
panel. **They are not built from the retired ``auc`` h5ad that is still on disk**, because a figure
that can only be produced from a target the pipeline no longer writes is not reproducible by a
standard run -- and one that renders anyway is worse than one that is absent, since nothing on its
face says which target it came from. Each is skipped with a printed reason; the superseded PNGs are
in ``docs/figures/archive/`` with a README. Every skip clears by itself at R4.

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
# The 27.07.2026 training run's outputs, archived 12.08.2026: they were produced on the void 8-drug
# panel and cannot be recreated by a standard run, so they moved out of outputs/panel/. Only the
# archived figures read them, and each guards on existence first.
LEGACY_PANEL = ROOT / "notebooks" / "outputs" / "legacy" / "panel_void_8drug"

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


#: The drug panel — the order every figure uses. Read from panel.csv, not maintained here.
PANEL = _panel()


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
    ax.text(50, 93.5, "as of 2026-08-12   ·   reference: project_planning_v2.pdf   ·   steps: docs/steps/",
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
         "overlap 190* lines · 180 trainable"], GREEN, GREEN_FILL)
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

    box(ax, XS[0], ROW_B, W, H, "05 · Multi-task + fair eval",
        ["K=545 · out-of-fold over 153 lines", "results WITHDRAWN 12.08.2026",
         "re-measured at R4 of the sweep"], GREY, GREY_FILL, title_color=GREY)
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

    ax.text(99.5, 1.2, "* 190 = name-matches in CTRPv2's roster; 180 = lines with actual post-QC measurements",
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


#: One real held-out cell line, scGPT, unweighted run (notebooks/outputs/legacy/panel_void_8drug/panel_oof_predictions.csv,
#: fold in which SKES1_BONE was held out). Predicted vs measured AUC for the eight panel drugs. Used
#: instead of an invented vector so the figure shows the actual output scale -- including the visible
#: shrinkage (predictions span 0.37-0.75 against a measured 0.04-0.91).
EXAMPLE_LINE = "SKES1_BONE"
EXAMPLE_PRED = [0.438, 0.663, 0.384, 0.372, 0.748, 0.545, 0.666, 0.739]
EXAMPLE_TRUE = [0.511, 0.797, 0.122, 0.036, 0.692, 0.177, 0.655, 0.908]


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
    # Non-compact lower bound dropped 14 -> 11 on 12.08.2026 to make room for the "plotted" note
    # below the uncentred-target paragraph. Compact is untouched: it is embedded in pipeline.png,
    # which has its own layout, and it draws neither of those notes.
    ax.set_xlim(0, 62 if compact else 100); ax.set_ylim(23 if compact else 11, 47 if compact else 52)
    ax.set_aspect("equal"); ax.axis("off")

    if not compact:
        ax.text(0, 51, "Model architecture — one cell in, one AUC per panel drug out",
                ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
        # The bars come from figure_data.npz, built on the superseded run, so "winsorized at 1.1",
        # "8-drug" and "25 epochs" are accurate about what is PLOTTED and wrong only if read as the
        # current pipeline. Split into two lines on 12.08.2026 (Selin) so the distinction is on the
        # figure rather than only in the README caption.
        ax.text(0, 48.0, "per-cell MLP · trained with the density-weighted masked MSE (loss_01–03)",
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
    counts = [6, 5, 4, 8]
    cy, r = 37, 1.25 * (0.72 if compact else 1)
    sp = 3.0 if not compact else 2.3
    head_sp = 2.05 if not compact else 1.62  # 8 heads drawn in full need tighter spacing
    pos = []
    for lx, n in zip(layers_x, counts):
        s = head_sp if n == 8 else sp
        pos.append([(lx, cy + (i - (n - 1) / 2) * s) for i in range(n)])
    for a, b in zip(pos[:-1], pos[1:]):
        for (x1, y1) in a:
            for (x2, y2) in b:
                ax.plot([x1, x2], [y1, y2], color="#bcd0e6", lw=0.5, zorder=1)
    for li, layer in enumerate(pos):
        rr = (0.85 if not compact else 0.55) if li == 3 else r
        for (x, y) in layer:
            ax.add_patch(Circle((x, y), rr, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.5, zorder=3))
    for lx, n in zip(layers_x[:3], counts[:3]):  # '...' only where neurons are omitted
        ax.text(lx, cy - ((n - 1) / 2) * sp - r - 0.5, "⋮", ha="center", va="top",
                fontsize=12 * fs, color=GREY)
    for lxx, t in zip(layers_x, ["input\n512", "hidden\n128", "hidden\n64", "heads\n8 drugs"]):
        ax.text(lxx, 25.8, t, ha="center", va="top", fontsize=9 * fs, color=INK)
    ax.add_patch(FancyBboxPatch((layers_x[0] - 3.5, 27.5), layers_x[3] - layers_x[0] + 7, 19,
                 boxstyle="round,pad=0.3,rounding_size=1.2",
                 fill=False, edgecolor=BLUE, lw=1.2, linestyle="--", zorder=0))
    if not compact:
        # Epoch count moved to the "plotted" line above: it is a property of the run shown, not of
        # the architecture, and the cap is itself under review (item 10 owns whether TrainConfig's
        # default of 25 moves to the 50 that 4a_percell_training passes). Naming it here made the
        # architecture caption go stale every time the cap moved.
        ax.text(43.2, 22.6, "hidden block  =  Linear → LayerNorm → GELU → Dropout 0.5      ·      "
                            "input dropout 0.1      ·      Adam, early stopping (patience 10)",
                ha="center", va="top", fontsize=8.2, color=INK)
        ax.text(43.2, 19.8, "the 8 heads are the 8 rows of one Linear(64 → 8) over a shared trunk — "
                            "there is no per-drug sub-network",
                ha="center", va="top", fontsize=8.2, color=GREY, style="italic")

    # ---------- OUTPUT: predicted raw AUC per drug ----------
    amax = 1.15                                # AUC 0 .. 1.15 maps to x0 .. x0+span
    cm = plt.colormaps["coolwarm"]             # low AUC = sensitive (blue) .. high = resistant (red)
    shade = lambda v: cm(np.clip(0.5 + (v - 0.5) / 1.2, 0, 1))  # white anchored at AUC 0.5
    ys = [cy + (i - 3.5) * head_sp for i in range(8)][::-1]
    for j, (yy, p, t) in enumerate(zip(ys, EXAMPLE_PRED, EXAMPLE_TRUE)):
        arrow(ax, layers_x[3] + 1.0, yy, x0 - (5.5 if compact else 6.5), yy, color="#9db6cf")
        ax.plot([x0, x0 + span], [yy, yy], color="#e4e4e0", lw=0.8, zorder=1)
        ax.add_patch(Rectangle((x0, yy - 0.62), p / amax * span, 1.24, facecolor=shade(p),
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
    ax.text(0, 16.0,
            "The target is uncentred:  the head bias is initialized to the fold's per-drug mean AUC, "
            "and biases and LayerNorm are excluded from weight decay,\nbecause the bias must sit near "
            "the drug's mean (~0.7) and decay would pull it to 0.  Any per-drug scaling belongs in the "
            "loss, not in the target.",
            ha="left", va="top", fontsize=8.2, color=INK)

    # The bars come from figure_data.npz, built on the superseded run, so "winsorized at 1.1",
    # "8-drug" and "25 epochs" are accurate about what is PLOTTED and wrong only if read as the
    # current pipeline. Stated on the figure rather than only in the README caption. Plain ASCII:
    # matplotlib's default font has no glyph for the warning emoji and renders it as a tofu box.
    ax.text(0, 12.6,
            "PLOTTED: the superseded run — auc winsorized at 1.1, the void 8-drug panel, 25 epochs."
            "    CURRENT PIPELINE: auc_cc, no winsorization, the rebuilt 11-drug panel, 50 epochs.",
            ha="left", va="top", fontsize=8.0, color=RED)


def build_architecture():
    if _needs_data("model_architecture.png"):
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


def panel_corr():
    """``notebooks/outputs/legacy/panel_void_8drug/panel_per_drug_correlation.csv`` as a DataFrame."""
    import pandas as pd

    return pd.read_csv(LEGACY_PANEL / "panel_per_drug_correlation.csv")


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
    if _needs_data("pipeline.png"):
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
    ax.text(x0 + 5 * cw, y0 + 6 * ch + 0.5, "545 drugs  →", ha="center", va="bottom",
            fontsize=7.6, color=MUTED)
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
    for i, (n, lab, col) in enumerate([(545, "545  CTRPv2 compounds", "#c3dcef"),
                                       (173, "173  FDA / clinical", "#6ba7d6"),
                                       (8, "8  the panel", BLUE)]):
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
    stage(41.0, ROW2_TITLE, "6", "Loss",
          "unscreened pairs dropped, rare\nresponse values weighted up", ROW2_CAP)
    ax = fig.add_axes([0.405, 0.215, 0.215, 0.16])
    draw_loss_objective(ax, compact=True)

    # ================================================== 7 · evaluation
    stage(64.5, ROW2_TITLE, "7", "Evaluation",
          "cells → one value per line,\nSpearman within each drug", ROW2_CAP)
    ax = fig.add_axes([0.655, 0.155, 0.125, 0.255])
    oof = pd.read_csv(LEGACY_PANEL / "panel_oof_predictions.csv")
    g8 = oof[(oof.rep == "X_scGPT") & (~oof.weighted) & (oof.drug == "dasatinib")]
    ax.scatter(g8.y_true, g8.y_pred, s=11, color=BLUE, alpha=0.55, edgecolor="white", lw=0.35)
    rho = corr[(corr.rep == "X_scGPT") & (~corr.weighted) & (corr.drug == "dasatinib")].spearman.iloc[0]
    ax.set_xlabel("measured AUC", fontsize=7.2, labelpad=1)
    ax.set_ylabel("predicted", fontsize=7.2, labelpad=1)
    ax.text(0.04, 0.96, f"dasatinib\nρ = {rho:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.6, color=INK)
    _tidy(ax, grid="both")
    ax.tick_params(labelsize=6.5)

    # ================================================== 8 · result
    stage(82.0, ROW2_TITLE, "8", "Result",
          "one seed — only scGPT clears\nthe ridge control", ROW2_CAP)
    ax = fig.add_axes([0.868, 0.155, 0.115, 0.255])
    ridge = pd.read_csv(LEGACY_PANEL / "panel_ridge_baseline.csv")
    bars = [
        ("scGPT MLP", corr[(corr.rep == "X_scGPT") & (~corr.weighted)].spearman.mean(), AMBER),
        ("PCA MLP", corr[(corr.rep == "X_pca") & (~corr.weighted)].spearman.mean(), GREY),
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
#: density exponent `alpha`, Huber's `beta`. Rendering a concrete value here would pre-empt a
#: decision that is Selin's and that the R4 run exists to make.
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
        r"\ell \in \{\text{squared},\, \text{absolute},\, \text{Huber}_\beta\},"
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
        "%   which is also what docs/figures/loss_01_objective.png renders. One source, both",
        "%   destinations, so the equation in the report and the equation in the figure cannot",
        "%   drift apart. To change a formula, edit LOSS_TEX_MACROS and re-run the generator.",
        "%",
        "% HOW THEY ARE USED",
        "%   main.tex does \\input{loss_objective} in the preamble; a section then writes e.g.",
        "%     \\begin{equation}\\LossObjective\\end{equation}",
        "%   Definitions only -- the section chooses the environment.",
        "%",
        "% WHAT IS DELIBERATELY SYMBOLIC",
        "%   \\ell, \\alpha and \\beta are left as symbols: the loss comparison (docs/TODO.md item",
        "%   9A) sweeps MSE / MAE / Huber against alpha in {off, 0.5, 1.0}, and Huber's beta must",
        "%   be derived from the residual scale (item 9C) rather than inherited. A concrete value",
        "%   here would pre-empt a decision the R4 run exists to make.",
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
    (``docs/TODO.md`` item 9A) sweeps MSE / MAE / Huber, and hiding the cell/cell-line split that is
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
        # loss the run has not chosen yet. beta is Selin's to derive (item 9C), not TrainConfig's.
        ax.text(1, 86, r"$\ell$, $\alpha$ and $\beta$ are what the loss comparison sweeps "
                       "(item 9A) — the figure fixes none of them",
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
         "squared / absolute / Huber$_\\beta$\n$n$ = cell,   $i$ = cell line,   $j$ = drug"),
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
    fig.suptitle("The weighting fired — and the ranking did not move", x=0.02, ha="left",
                 fontsize=13, fontweight="bold", color=INK, y=0.99)
    fig.text(0.02, 0.925, "one point per drug, out of fold, one seed — unweighted against "
                          "density-weighted",
             ha="left", va="top", fontsize=8.8, color=GREY)

    for ax, col, label, lim in [(a1, "pred_std", "spread of the predictions", (0.03, 0.11)),
                                (a2, "spearman", "per-drug Spearman", (0.0, 0.65))]:
        ax.plot(lim, lim, color="#c9c9c4", lw=1.0, ls="--", zorder=1)
        for rep, c, name in reps:
            u = corr[(corr.rep == rep) & (~corr.weighted)].set_index("drug").reindex(PANEL)[col]
            w = corr[(corr.rep == rep) & (corr.weighted)].set_index("drug").reindex(PANEL)[col]
            ax.scatter(u, w, s=46, color=c, alpha=0.85, edgecolor="white", lw=0.8,
                       label=name, zorder=3)
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.set_xlabel(f"unweighted — {label}", fontsize=9)
        ax.set_ylabel(f"density-weighted — {label}", fontsize=9)
        ax.set_aspect("equal")
        _tidy(ax, grid="both")
        ax.tick_params(labelsize=8)
    a1.legend(frameon=False, fontsize=8.5, loc="lower right")
    a1.text(0.04, 0.96, "above the line:\nthe model hedges less", transform=a1.transAxes,
            ha="left", va="top", fontsize=8.5, color=BLUE)
    a2.text(0.04, 0.96, "on the line:\nno gain, no loss", transform=a2.transAxes,
            ha="left", va="top", fontsize=8.5, color=MUTED)

    out = FIG / "loss_03_effect.png"
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
