"""Generate the OncoTox slide/doc graphics:

  1. docs/pipeline_overview.png   — status overview of the whole pipeline (steps 01-08)
  2. docs/model_architecture.png  — input + model + task on one figure (to merge slides)
  3. docs/pipeline_full.png       — the pipeline as a data flow: data -> split -> representation
                                    (PCA vs scGPT) -> drug panel -> model+loss -> evaluation
  4. docs/loss_function.png       — anatomy of the density-weighted masked MSE (reweighted
                                    regression), with the real weight curve of one panel drug

Figures 2-4 describe the setup as of 27.07.2026: target = raw AUC (winsorized at 1.1), no per-drug
z-score, the per-drug scaling moved into the loss as an inverse-label-density sample weight.

Green = done / on-plan · Amber = addition beyond plan · Red (dashed) = not started.

Run:  uv run docs/make_pipeline_overview.py
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
PANEL_OUT = ROOT / "notebooks" / "outputs" / "panel"

GREEN = "#2e7d32"; GREEN_FILL = "#c8e6c9"
AMBER = "#b8860b"; AMBER_FILL = "#ffe9b3"
RED = "#c62828"; RED_FILL = "#ffcdd2"
BLUE = "#1f6fb2"; BLUE_FILL = "#dbe7f3"
GREY = "#777777"; GREY_FILL = "#e8e8e8"
INK = "#1a1a1a"
MUTED = "#52514e"

#: The literature panel (docs/steps/01, progress report §2) — the order every figure uses.
PANEL = ["methotrexate", "dasatinib", "paclitaxel", "vincristine",
         "afatinib", "topotecan", "tanespimycin", "selumetinib"]


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
    ax.text(50, 93.5, "as of 2026-07-14   ·   reference: project_planning_v2.pdf   ·   steps: docs/steps/",
            ha="center", va="top", fontsize=9.5, color=GREY)
    handles = [
        mpatches.Patch(facecolor=GREEN_FILL, edgecolor=GREEN, label="Done / on-plan"),
        mpatches.Patch(facecolor=AMBER_FILL, edgecolor=AMBER, label="Addition beyond plan"),
        mpatches.Patch(facecolor=RED_FILL, edgecolor=RED, label="Not started (planned)"),
    ]
    ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.885),
              ncol=3, fontsize=9.5, frameon=True, framealpha=0.9)

    box(ax, XS[0], ROW_A, W, H, "01 · Datasets & harmonization",
        ["SCP542 53,513 cells x 22,722 g", "CTRPv2 545 drugs · target: auc_z",
         "overlap 190* lines · 180 trainable"], GREEN, GREEN_FILL)
    box(ax, XS[1], ROW_A, W, H, "02 · Preprocessing & embeddings",
        ["scGPT X_scGPT = 512-d", "gene-set sweep 1k-5k + all_genes",
         "X_pca = 512-d · cancer-type UMAPs"], GREEN, GREEN_FILL)
    box(ax, XS[2], ROW_A, W, H, "03 · Model & training design",
        ["per-cell input -> viability", "masked MSE · matched (128,64) MLP",
         "PCA & scGPT both 512-d"], GREEN, GREEN_FILL)
    box(ax, XS[3], ROW_A, W, H, "04 · Single-task baseline",
        ["paclitaxel, leak-free split", "best scGPT val MSE 0.0336",
         "1 DB · 1 score · 1 drug"], GREEN, GREEN_FILL)

    box(ax, XS[0], ROW_B, W, H, "05 · Multi-task + fair eval",
        ["K=545 · out-of-fold over 153 lines", "target fix (auc_z): rho ~0 -> ~0.4",
         "scGPT >= PCA · benchmarked (DrEval)"], GREEN, GREEN_FILL)
    box(ax, XS[1], ROW_B, W, H, "06 · Cross-database  (MISSING)",
        ["CTRPv2 + PRISM + GDSC", "efficacy + toxicity heads",
         "the 'combine all' goal"], RED, RED_FILL, title_color=RED, dashed=True)
    box(ax, XS[2], ROW_B, W, H, "07 · XAI / interpretability  (MISSING)",
        ["feature importance -> drivers", "stretch goal", "not started"],
        RED, RED_FILL, title_color=RED, dashed=True)
    box(ax, XS[3], ROW_B, W, H, "08 · Foundation model  (HORIZON)",
        ["reusable pan-cancer FM", "fine-tune on clinical (binary)",
         "overarching main goal"], RED, RED_FILL, title_color=RED, dashed=True)

    BAND_Y, BAND_H = 5, 13
    box(ax, XS[0], BAND_Y, 94, BAND_H, "Additions beyond the written plan",
        ["512-d PCA (matched to scGPT)  ·  out-of-fold CV over 153 lines  ·  per-drug correlation metric  ·  "
         "gene-set sweep  ·  cancer-type UMAPs  ·  cell-line-grouped split (leak fix)  ·  run versioning\n"
         "per-drug z-scored target (auc_z = 1/sigma^2 head weighting)  ·  learnability filter  ·  "
         "ridge line-level control  ·  external benchmark against DrEval (drevalpy)"],
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

    out = HERE / "pipeline_overview.png"
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


#: One real held-out cell line, scGPT, unweighted run (notebooks/outputs/panel/panel_oof_predictions.csv,
#: fold in which SKES1_BONE was held out). Predicted vs measured AUC for the eight panel drugs. Used
#: instead of an invented vector so the figure shows the actual output scale -- including the visible
#: shrinkage (predictions span 0.37-0.75 against a measured 0.04-0.91).
EXAMPLE_LINE = "SKES1_BONE"
EXAMPLE_PRED = [0.438, 0.663, 0.384, 0.372, 0.748, 0.545, 0.666, 0.739]
EXAMPLE_TRUE = [0.511, 0.797, 0.122, 0.036, 0.692, 0.177, 0.655, 0.908]


def build_architecture():
    fig, ax = plt.subplots(figsize=(17.0, 6.6))
    ax.set_xlim(0, 100); ax.set_ylim(14, 52); ax.set_aspect("equal"); ax.axis("off")

    ax.text(0, 51, "Model architecture — one cell in, one AUC per panel drug out",
            ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
    ax.text(0, 48.0, "per-cell MLP · target = raw AUC (winsorized at 1.1) · 8-drug literature panel · "
                     "trained with the density-weighted masked MSE (see loss_function.png)",
            ha="left", va="top", fontsize=8.5, color=GREY)

    # ---------- INPUT: one cell -> embedding vector ----------
    ax.add_patch(Circle((6, 37), 2.7, facecolor="#fde0c5", edgecolor="#d2691e", lw=1.8, zorder=3))
    for dx, dy in [(-0.9, 0.5), (0.7, -0.4), (0.2, 1.0), (-0.3, -0.9)]:
        ax.add_patch(Circle((6 + dx, 37 + dy), 0.55, facecolor="#d2691e", lw=0, zorder=4))
    ax.text(6, 32.3, "single cell\n(scRNA-seq)", ha="center", va="top", fontsize=9, color=INK)
    arrow(ax, 9, 37, 11.3, 37, color=INK)

    _heat_strip(ax, 14, 30.5, 43.5, np.linspace(0.05, 0.95, 14), "viridis")
    ax.text(14, 29.6, "512-d embedding", ha="center", va="top", fontsize=9.5, fontweight="bold", color=BLUE)
    ax.text(14, 26.9, "PCA  or  scGPT\n(frozen — never fine-tuned)", ha="center", va="top",
            fontsize=8.2, color=INK)
    arrow(ax, 15.8, 37, 22.5, 37, color=INK)

    # ---------- MODEL: MLP drawn as neurons ----------
    layers_x = [26, 39, 51, 62]
    counts = [6, 5, 4, 8]
    cy, sp, r = 37, 3.0, 1.25
    head_sp = 2.05  # 8 heads are drawn in full, so they need tighter spacing than the trunk
    pos = []
    for lx, n in zip(layers_x, counts):
        s = head_sp if n == 8 else sp
        pos.append([(lx, cy + (i - (n - 1) / 2) * s) for i in range(n)])
    for a, b in zip(pos[:-1], pos[1:]):
        for (x1, y1) in a:
            for (x2, y2) in b:
                ax.plot([x1, x2], [y1, y2], color="#bcd0e6", lw=0.5, zorder=1)
    for li, layer in enumerate(pos):
        rr = 0.85 if li == 3 else r
        for (x, y) in layer:
            ax.add_patch(Circle((x, y), rr, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.5, zorder=3))
    for lx, n in zip(layers_x[:3], counts[:3]):  # '...' only where neurons are omitted
        ax.text(lx, cy - (n / 2) * sp - 0.6, "⋮", ha="center", va="top", fontsize=12, color=GREY)
    for lx, t in zip(layers_x, ["input\n512", "hidden\n128", "hidden\n64", "heads\n8 drugs"]):
        ax.text(lx, 25.5, t, ha="center", va="top", fontsize=9, color=INK)
    ax.add_patch(FancyBboxPatch((22.5, 27.5), 42.5, 19, boxstyle="round,pad=0.3,rounding_size=1.2",
                 fill=False, edgecolor=BLUE, lw=1.2, linestyle="--", zorder=0))
    ax.text(43.7, 21.4, "hidden block  =  Linear → LayerNorm → GELU → Dropout 0.5      ·      "
                        "input dropout 0.1      ·      Adam, 25 epochs, early stopping",
            ha="center", va="top", fontsize=8.2, color=INK)
    ax.text(43.7, 18.3, "the 8 heads are the 8 rows of one Linear(64 → 8) over a shared trunk — "
                        "there is no per-drug sub-network",
            ha="center", va="top", fontsize=8.2, color=GREY, style="italic")

    # ---------- OUTPUT: predicted raw AUC per drug ----------
    x0, span, amax = 71.5, 20.0, 1.15          # AUC 0 .. 1.15 maps to x0 .. x0+span
    cm = plt.colormaps["coolwarm_r"]           # low AUC = sensitive (blue) .. high = resistant (red)
    ys = [cy + (i - 3.5) * head_sp for i in range(8)][::-1]
    for j, (yy, p, t) in enumerate(zip(ys, EXAMPLE_PRED, EXAMPLE_TRUE)):
        arrow(ax, 63.2, cy + (j - 3.5) * head_sp, 66.0, cy + (j - 3.5) * head_sp, color="#9db6cf")
        ax.plot([x0, x0 + span], [yy, yy], color="#e4e4e0", lw=0.8, zorder=1)
        ax.add_patch(Rectangle((x0, yy - 0.62), p / amax * span, 1.24, facecolor=cm(p / amax),
                               edgecolor="white", lw=0.6, zorder=3))
        ax.plot([x0 + t / amax * span], [yy], marker="D", ms=3.4, color=INK, zorder=4)
        ax.text(x0 - 0.8, yy, PANEL[j], ha="right", va="center", fontsize=8, color=INK)
    for v in (0.0, 0.5, 1.0):
        ax.plot([x0 + v / amax * span] * 2, [ys[-1] + 1.2, ys[0] - 1.2], color=GREY, lw=0.6,
                linestyle=":", zorder=1)
        ax.text(x0 + v / amax * span, ys[0] - 1.8, f"{v:.1f}", ha="center", va="top",
                fontsize=7.5, color=GREY)
    ax.text(x0 + span / 2, 29.6, "predicted AUC per drug", ha="center", va="top",
            fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x0, 27.0, "low AUC = sensitive (the drug kills this line)", ha="left", va="top",
            fontsize=8, color=BLUE)
    ax.text(x0, 24.6, "high AUC = resistant", ha="left", va="top", fontsize=8, color=RED)
    ax.text(x0, 22.2, f"bars: out-of-fold prediction for {EXAMPLE_LINE}\n"
                      "◆ = the measured CTRPv2 value (never seen in training)",
            ha="left", va="top", fontsize=7.6, color=GREY)

    ax.text(0, 16.0,
            "Raw AUC, not the old per-drug z-score:  the head bias is initialized to the fold's per-drug "
            "mean AUC and the output layer is excluded from weight decay,\nbecause on an uncentred target "
            "the bias must sit near 0.7 and decay would pull it to 0.  The per-drug scaling that auc_z "
            "used to apply now lives in the loss.",
            ha="left", va="top", fontsize=8.2, color=INK)

    out = HERE / "model_architecture.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    build_pipeline()
    build_architecture()
