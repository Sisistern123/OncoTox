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
    """Vertical heatmap strip (a 'vector') of len(vals) cells."""
    ys = np.linspace(y0, y1, len(vals) + 1)
    cm = plt.colormaps[cmap]
    for i, v in enumerate(vals):
        ax.add_patch(Rectangle((xc - w / 2, ys[i]), w, ys[i + 1] - ys[i],
                     facecolor=cm(v), edgecolor="white", lw=0.4, zorder=3))
    ax.add_patch(Rectangle((xc - w / 2, ys[0]), w, ys[-1] - ys[0], fill=False,
                 edgecolor=INK, lw=1.3, zorder=4))


def build_architecture():
    fig, ax = plt.subplots(figsize=(16.0, 8.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 50); ax.set_aspect("equal"); ax.axis("off")

    ax.text(50, 49, "OncoTox — Input, Model & Task", ha="center", va="top",
            fontsize=16, fontweight="bold", color=INK)

    # ---------- INPUT: one cell -> embedding vector ----------
    ax.add_patch(Circle((6, 37), 2.7, facecolor="#fde0c5", edgecolor="#d2691e", lw=1.8, zorder=3))
    for dx, dy in [(-0.9, 0.5), (0.7, -0.4), (0.2, 1.0), (-0.3, -0.9)]:
        ax.add_patch(Circle((6 + dx, 37 + dy), 0.55, facecolor="#d2691e", lw=0, zorder=4))
    ax.text(6, 32.3, "single cell\n(scRNA-seq)", ha="center", va="top", fontsize=9, color=INK)
    arrow(ax, 9, 37, 11.3, 37, color=INK)

    _heat_strip(ax, 14, 30.5, 43.5, np.linspace(0.05, 0.95, 14), "viridis")
    ax.text(14, 29.6, "512-d embedding", ha="center", va="top", fontsize=9.5, fontweight="bold", color=BLUE)
    ax.text(14, 26.7, "PCA  or  scGPT", ha="center", va="top", fontsize=9, color=INK)
    arrow(ax, 15.8, 37, 22.5, 37, color=INK)

    # ---------- MODEL: MLP drawn as neurons ----------
    layers_x = [26, 40, 54, 66]
    counts = [6, 5, 4, 6]
    cy, sp, r = 37, 3.0, 1.25
    pos = [[(lx, cy + (i - (n - 1) / 2) * sp) for i in range(n)] for lx, n in zip(layers_x, counts)]
    for a, b in zip(pos[:-1], pos[1:]):
        for (x1, y1) in a:
            for (x2, y2) in b:
                ax.plot([x1, x2], [y1, y2], color="#bcd0e6", lw=0.5, zorder=1)
    for layer in pos:
        for (x, y) in layer:
            ax.add_patch(Circle((x, y), r, facecolor=BLUE_FILL, edgecolor=BLUE, lw=1.6, zorder=3))
    for lx, n in zip(layers_x, counts):
        ax.text(lx, cy - (n / 2) * sp - 0.6, "⋮", ha="center", va="top", fontsize=12, color=GREY)
    for lx, t in zip(layers_x, ["input\n512", "hidden\n128", "hidden\n64", "heads\n545 drugs"]):
        ax.text(lx, 25.5, t, ha="center", va="top", fontsize=9, color=INK)
    ax.text(46, 47.0, "OncoMLP  ·  LayerNorm + dropout 0.5", ha="center", va="top",
            fontsize=10.5, fontweight="bold", color=BLUE)
    ax.add_patch(FancyBboxPatch((22.5, 27.5), 47, 19, boxstyle="round,pad=0.3,rounding_size=1.2",
                 fill=False, edgecolor=BLUE, lw=1.2, linestyle="--", zorder=0))
    arrow(ax, 68, 37, 75, 37, color=INK)

    # ---------- OUTPUT: predicted viability per drug ----------
    out_vals = np.array([0.95, 0.4, 0.88, 0.2, 0.7, 0.55, 0.97, 0.35, 0.8, 0.6, 0.9, 0.45, 0.75, 0.3])
    _heat_strip(ax, 79, 30.5, 43.5, out_vals, "RdYlGn")
    ax.text(80.8, 43.2, "green ≈ 1: no effect", ha="left", va="center", fontsize=7.5, color=GREEN)
    ax.text(80.8, 30.8, "red < 1: killed", ha="left", va="center", fontsize=7.5, color=RED)
    ax.text(79, 29.6, "viability per drug", ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(79, 26.7, "≈1 = no effect (relative\nto control; can exceed 1)", ha="center", va="top",
            fontsize=8.0, color=INK)

    # ---------- TASK: why it's hard (bulk label broadcast) ----------
    ax.add_patch(FancyBboxPatch((3, 2), 94, 17, boxstyle="round,pad=0.3,rounding_size=1.5",
                 facecolor="#fbf4e6", edgecolor=AMBER, lw=1.6, zorder=0))
    ax.text(6, 17.2, "The task — and why it is hard", ha="left", va="top",
            fontsize=11, fontweight="bold", color=AMBER)

    # bulk value broadcast to a line's cells
    ax.add_patch(FancyBboxPatch((9, 12.0), 12, 3.2, boxstyle="round,pad=0.2,rounding_size=0.8",
                 facecolor="white", edgecolor=AMBER, lw=1.4, zorder=3))
    ax.text(15, 13.6, "bulk viability 0.8", ha="center", va="center", fontsize=8.5, color=INK)
    cell_x = [10.5, 15, 19.5]
    for cxp in cell_x:
        ax.add_patch(Circle((cxp, 7.2), 1.3, facecolor="#fde0c5", edgecolor="#d2691e", lw=1.3, zorder=3))
        arrow(ax, 15, 11.9, cxp, 8.6, color="#d2691e")
    ax.text(15, 4.6, "1 cell line (~300 cells)", ha="center", va="top", fontsize=8, color=INK)
    ax.text(26, 13.2,
            "One BULK value per (cell line × drug) is copied to every cell of that line\n"
            "→ noisy per-cell labels; the model can only learn per-LINE signal (~126 train lines).",
            ha="left", va="top", fontsize=9, color=INK)
    ax.text(26, 7.4,
            "Masked MSE: loss only on measured (cell, drug) pairs   ·   split grouped by cell line\n"
            "(leak-free) — evaluate on unseen lines vs a per-drug-mean baseline + per-drug correlation.",
            ha="left", va="top", fontsize=9, color=INK)

    out = HERE / "model_architecture.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    build_pipeline()
    build_architecture()
