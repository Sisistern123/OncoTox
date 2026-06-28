"""Generate the OncoTox slide/doc graphics:

  1. docs/pipeline_overview.png   — status overview of the whole pipeline (steps 01-08)
  2. docs/model_architecture.png  — input + model + task on one figure (to merge slides)

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

GREEN = "#2e7d32"; GREEN_FILL = "#c8e6c9"
AMBER = "#b8860b"; AMBER_FILL = "#ffe9b3"
RED = "#c62828"; RED_FILL = "#ffcdd2"
BLUE = "#1f6fb2"; BLUE_FILL = "#dbe7f3"
GREY = "#777777"; GREY_FILL = "#e8e8e8"
INK = "#1a1a1a"


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
    ax.text(50, 93.5, "as of 2026-06-28   ·   reference: project_planning_v2.pdf   ·   steps: docs/steps/",
            ha="center", va="top", fontsize=9.5, color=GREY)
    handles = [
        mpatches.Patch(facecolor=GREEN_FILL, edgecolor=GREEN, label="Done / on-plan"),
        mpatches.Patch(facecolor=AMBER_FILL, edgecolor=AMBER, label="Addition beyond plan"),
        mpatches.Patch(facecolor=RED_FILL, edgecolor=RED, label="Not started (planned)"),
    ]
    ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.885),
              ncol=3, fontsize=9.5, frameon=True, framealpha=0.9)

    box(ax, XS[0], ROW_A, W, H, "01 · Datasets & harmonization",
        ["SCP542 53,513 cells x 22,722 g", "CTRPv2 545 drugs (cpd_avg_pv)",
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
        ["K=545 drugs · 5-fold GroupKFold CV", "PCA ~ scGPT (within noise)",
         "per-drug rho ~ 0 · 1 DB·1 score"], GREEN, GREEN_FILL)
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
        ["512-d PCA (matched to scGPT)  ·  5-fold GroupKFold CV  ·  per-drug correlation metric  ·  "
         "gene-set sweep (1k-5k + all_genes)  ·  cancer-type UMAPs  ·  cell-line-grouped split (leak fix)  ·  "
         "per-drug-mean baseline + run versioning"],
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
    """Vertical heatmap strip (the predicted 545-vector)."""
    ys = np.linspace(y0, y1, len(vals) + 1)
    cm = plt.colormaps[cmap]
    for i, v in enumerate(vals):
        ax.add_patch(Rectangle((xc - w / 2, ys[i]), w, ys[i + 1] - ys[i],
                     facecolor=cm(v), edgecolor="white", lw=0.4, zorder=3))
    ax.add_patch(Rectangle((xc - w / 2, ys[0]), w, ys[-1] - ys[0], fill=False,
                 edgecolor=INK, lw=1.3, zorder=4))


def build_architecture():
    # equal aspect so the neurons are round
    fig, ax = plt.subplots(figsize=(16.0, 8.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 50); ax.set_aspect("equal"); ax.axis("off")

    ax.text(50, 49, "OncoTox — Input, Model & Task", ha="center", va="top",
            fontsize=15, fontweight="bold", color=INK)

    # ---------- one cell feeds the input layer ----------
    ax.add_patch(Circle((6, 35), 2.6, facecolor="#fde0c5", edgecolor="#d2691e", lw=1.8, zorder=3))
    for dx, dy in [(-0.8, 0.4), (0.6, -0.4), (0.1, 0.9)]:
        ax.add_patch(Circle((6 + dx, 35 + dy), 0.5, facecolor="#d2691e", lw=0, zorder=4))
    ax.text(6, 30.4, "single cell\n(scRNA-seq)", ha="center", va="top", fontsize=9, color=INK)
    arrow(ax, 8.7, 35, 16.0, 35, color=INK)

    # ---------- MLP drawn as neurons (size labelled under each layer) ----------
    layers_x = [19, 34, 49, 63]
    counts = [6, 5, 4, 6]
    cy, sp, r = 35, 3.0, 1.25
    pos = [[(lx, cy + (i - (n - 1) / 2) * sp) for i in range(n)] for lx, n in zip(layers_x, counts)]
    for a, b in zip(pos[:-1], pos[1:]):
        for (x1, y1) in a:
            for (x2, y2) in b:
                ax.plot([x1, x2], [y1, y2], color="#bcd0e6", lw=0.5, zorder=1)
    for layer in pos:
        for (x, y) in layer:
            ax.add_patch(Circle((x, y), r, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.6, zorder=3))
    for lx, n in zip(layers_x, counts):
        ax.text(lx, cy - (n / 2) * sp - 0.4, "⋮", ha="center", va="top", fontsize=11, color=GREY)
    sizelabels = ["input\n512", "hidden\n128", "hidden\n64", "output heads\n545 drugs"]
    for lx, t in zip(layers_x, sizelabels):
        ax.text(lx, 23.5, t, ha="center", va="top", fontsize=9, color=INK)
    ax.text(19, 19.6, "PCA  or  scGPT", ha="center", va="top", fontsize=8.5,
            fontweight="bold", color=BLUE)
    ax.text(41, 46.6, "OncoMLP  ·  LayerNorm + dropout 0.5", ha="center", va="top",
            fontsize=10.5, fontweight="bold", color=BLUE)
    ax.add_patch(FancyBboxPatch((15.5, 25.0), 51.0, 19.5, boxstyle="round,pad=0.4,rounding_size=1.2",
                 fill=False, edgecolor=BLUE, lw=1.2, linestyle="--", zorder=0))
    arrow(ax, 65, 35, 73, 35, color=INK)

    # ---------- output: predicted viability per drug (a real 545-vector) ----------
    out_vals = np.array([0.95, 0.4, 0.88, 0.2, 0.7, 0.55, 0.97, 0.35, 0.8, 0.6, 0.9, 0.45, 0.75, 0.3])
    _heat_strip(ax, 77, 28, 42, out_vals, "RdYlGn")
    ax.text(78.9, 41.7, "1 = survives", ha="left", va="center", fontsize=7.5, color=GREY)
    ax.text(78.9, 28.3, "0 = killed", ha="left", va="center", fontsize=7.5, color=GREY)
    ax.text(77, 27.0, "predicted viability\nper drug (545-vector)", ha="center", va="top",
            fontsize=9, fontweight="bold", color=INK)

    # ---------- target + loss as diagram nodes (info woven in, no prose box) ----------
    box(ax, 12, 3.0, 34, 13.5, "Target — CTRPv2 viability (cpd_avg_pv)",
        ["1 BULK value per (cell line × drug)", "broadcast to all the line's cells → noisy labels",
         "1.0 = no effect · 0 = killed · ~126 train lines"],
        AMBER, AMBER_FILL, title_color=AMBER, title_size=9.5, body_size=8.3)
    box(ax, 64, 4.5, 21, 10.5, "Masked MSE",
        ["loss on observed", "(cell × drug) pairs only"],
        GREY, GREY_FILL, title_color=INK, title_size=10.0, body_size=8.5)

    arrow(ax, 77, 28, 74.5, 15.0, color=INK)        # prediction ↓ into loss
    arrow(ax, 46, 9.8, 64, 9.8, color=INK)          # target → loss

    ax.text(50, 0.8,
            "leak-free cell-line-grouped 70/15/15 split  ·  eval vs per-drug-mean baseline + per-drug "
            "correlation (5-fold CV)", ha="center", va="bottom", fontsize=8.5, color=GREY, style="italic")

    out = HERE / "model_architecture.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    build_pipeline()
    build_architecture()
