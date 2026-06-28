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
def build_architecture():
    fig, ax = plt.subplots(figsize=(16.0, 8.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    ax.text(50, 97, "OncoTox — Input, Model & Task", ha="center", va="top",
            fontsize=17, fontweight="bold", color=INK)

    # ---------- one cell feeds the model ----------
    ax.add_patch(Circle((9.5, 90), 2.4, facecolor="#fde0c5", edgecolor="#d2691e", lw=1.8, zorder=3))
    for dx, dy in [(-0.8, 0.4), (0.6, -0.4), (0.1, 0.9)]:
        ax.add_patch(Circle((9.5 + dx, 90 + dy), 0.5, facecolor="#d2691e", lw=0, zorder=4))
    ax.text(13.2, 90, "one single cell (scRNA-seq) · one prediction per cell",
            ha="left", va="center", fontsize=9.5, color=INK)
    arrow(ax, 9.5, 87.5, 9.5, 78.5, color=INK)

    # ---------- forward block diagram (academic: a box per layer, sizes on arrows) ----------
    yB, hB = 58, 20
    blocks = [
        (2.0, 17.0, "Input", ["cell embedding", "512-d   (PCA or scGPT)"], BLUE, BLUE_FILL, BLUE),
        (24.0, 16.5, "Hidden layer 1", ["Linear 512 → 128", "LayerNorm · ReLU · dropout"], BLUE, BLUE_FILL, BLUE),
        (44.5, 15.0, "Hidden layer 2", ["Linear 128 → 64", "LayerNorm · ReLU · dropout"], BLUE, BLUE_FILL, BLUE),
        (63.0, 15.0, "Output layer", ["Linear 64 → 545", "one head per drug"], BLUE, BLUE_FILL, BLUE),
        (81.5, 16.5, "Prediction", ["viability / drug", "545 values ∈ [0, 1]"], GREEN, GREEN_FILL, GREEN),
    ]
    spans = []
    for x, w, title, lines, edge, fill, tc in blocks:
        box(ax, x, yB, w, hB, title, lines, edge, fill, title_color=tc, title_size=10.5, body_size=9.0)
        spans.append((x, w))
    dims = ["512", "128", "64", "545"]
    for (x, w), nx, d in zip(spans[:-1], [s[0] for s in spans[1:]], dims):
        arrow(ax, x + w, yB + hB / 2, nx, yB + hB / 2, color=INK)
        ax.text((x + w + nx) / 2, yB + hB / 2 + 2.0, d, ha="center", va="bottom",
                fontsize=8.5, color=GREY, style="italic")

    # OncoMLP bracket around the three trainable layers
    ax.add_patch(FancyBboxPatch((23.0, 54.5), 56.5, 27.5, boxstyle="round,pad=0.4,rounding_size=1.4",
                 fill=False, edgecolor=BLUE, lw=1.3, linestyle="--", zorder=0))
    ax.text(51.0, 81.0, "OncoMLP  —  MLP, 2 hidden layers", ha="center", va="bottom",
            fontsize=10.5, fontweight="bold", color=BLUE)

    # ---------- target + loss woven into the diagram ----------
    yT, hT = 22, 22
    box(ax, 60.0, yT, 23.5, hT, "Target — CTRPv2 viability",
        ["cpd_avg_pv  (1.0 = no effect, 0 = killed)", "ONE bulk value per (cell line × drug)",
         "broadcast to all the line's cells", "→ noisy per-cell labels · ~126 train lines"],
        AMBER, AMBER_FILL, title_color=AMBER, title_size=10.0, body_size=8.5)
    box(ax, 86.0, yT + 4, 12.5, hT - 8, "Masked MSE",
        ["loss on observed", "(cell × drug)", "pairs only"],
        GREY, GREY_FILL, title_color=INK, title_size=10.0, body_size=8.5)

    # prediction (top) and target (left) feed the loss node
    arrow(ax, 89.75, yB, 92.25, yT + (hT - 8) + 4, color=INK)        # prediction ↓ into loss
    arrow(ax, 83.5, yT + hT / 2, 86.0, yT + hT / 2, color=INK)        # target → loss

    ax.text(50, 9.0,
            "Trained on a leak-free, cell-line-grouped 70/15/15 split  ·  evaluated on unseen cell lines "
            "vs a per-drug-mean baseline + per-drug correlation (5-fold CV)",
            ha="center", va="center", fontsize=9.5, color=GREY, style="italic")

    out = HERE / "model_architecture.png"
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    build_pipeline()
    build_architecture()
