"""Generate docs/pipeline_overview.png — a status overview of the OncoTox pipeline.

Green  = done / on-plan
Amber  = done but an addition or only partially covering the plan
Red    = not started (still missing)

Run:  uv run docs/make_pipeline_overview.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "pipeline_overview.png"

GREEN = "#2e7d32"
GREEN_FILL = "#c8e6c9"
AMBER = "#b8860b"
AMBER_FILL = "#ffe9b3"
RED = "#c62828"
RED_FILL = "#ffcdd2"
GREY = "#777777"
INK = "#1a1a1a"

fig, ax = plt.subplots(figsize=(15.5, 9.0))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# Grid
W, H = 27, 17
XS = [4, 36.5, 69]
ROW_A, ROW_B, ROW_C = 67, 41, 9


def box(x, y, w, h, title, lines, edge, fill, title_color=None, dashed=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=1.6",
        linewidth=2.2, edgecolor=edge, facecolor=fill,
        linestyle="--" if dashed else "-", mutation_aspect=1.0, zorder=2,
    ))
    ax.text(x + w / 2, y + h - 3.2, title, ha="center", va="top",
            fontsize=11.5, fontweight="bold", color=title_color or INK, zorder=3)
    ax.text(x + w / 2, y + h - 8.0, "\n".join(lines), ha="center", va="top",
            fontsize=8.7, color=INK, zorder=3)


def arrow(x1, y1, x2, y2, color=INK, dashed=False):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
        linewidth=2.0, color=color, linestyle="--" if dashed else "-", zorder=1,
    ))


# Title + legend
ax.text(50, 98, "OncoTox Pipeline — Status Overview", ha="center", va="top",
        fontsize=17, fontweight="bold", color=INK)
ax.text(50, 93.5, "as of 2026-06-10   ·   reference: project_planning_v2.pdf",
        ha="center", va="top", fontsize=9.5, color=GREY)
handles = [
    mpatches.Patch(facecolor=GREEN_FILL, edgecolor=GREEN, label="Done / on-plan"),
    mpatches.Patch(facecolor=AMBER_FILL, edgecolor=AMBER, label="Addition / partial"),
    mpatches.Patch(facecolor=RED_FILL, edgecolor=RED, label="Missing (not started)"),
]
ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.875),
          ncol=3, fontsize=9.5, frameon=True, framealpha=0.9)

# Row A — done
box(XS[0], ROW_A, W, H, "1 · Data collection",
    ["SCP542 scRNA-seq", "53,513 cells x 22,722 genes", "198 cell lines",
     "CTRPv2: 545 drugs (cpd_avg_pv)"], GREEN, GREEN_FILL)
box(XS[1], ROW_A, W, H, "2 · Overlap & harmonization",
    ["SCP542 x CTRPv2: 190 lines*", "545 compounds, 100% non-null",
     "BRD-ID match 243 · DrugBank"], GREEN, GREEN_FILL)
box(XS[2], ROW_A, W, H, "3 · Embeddings + AnnData",
    ["scGPT  X_scGPT = 512-d", "HVG-5000 (4,576/5,000 vocab)",
     "+ X_pca baseline"], GREEN, GREEN_FILL)

# Row B — done
box(XS[0], ROW_B, W, H, "4 · Latent validation",
    ["UMAP scGPT vs PCA (Fig 3/4)", "scGPT = continuous manifold",
     "PCA = tissue 'islands'"], GREEN, GREEN_FILL)
box(XS[1], ROW_B, W, H, "5 · Single-task baseline",
    ["paclitaxel viability regression", "cell-line-grouped split (leak-free)",
     "best scGPT val MSE 0.0336"], GREEN, GREEN_FILL)
box(XS[2], ROW_B, W, H, "6 · Multi-task (masked MSE)",
    ["K=545 CTRP drugs, split_ctrp", "scGPT beats baseline 142/545",
     "(PCA 97/545)"], GREEN, GREEN_FILL)

# Row C — additions (amber) + missing (red)
box(XS[0], ROW_C, W, H, "Additions beyond plan",
    ["HVG-5000 + all_genes variant", "leak -> grouped-split fix",
     "per-drug-mean sanity baseline", "run versioning (runs/ + ledger)"],
    AMBER, AMBER_FILL, title_color=AMBER)
box(XS[1], ROW_C, W, H, "PRISM + GDSC heads  (MISSING)",
    ["true cross-database multi-task", "downloaded only, not integrated",
     "-> plan Phase 3 half done"], RED, RED_FILL, title_color=RED, dashed=True)
box(XS[2], ROW_C, W, H, "XAI / interpretability  (MISSING)",
    ["feature importance on drivers", "stretch goal, not started"],
    RED, RED_FILL, title_color=RED, dashed=True)

# Arrows — Row A
arrow(XS[0] + W, ROW_A + H / 2, XS[1], ROW_A + H / 2)
arrow(XS[1] + W, ROW_A + H / 2, XS[2], ROW_A + H / 2)
# wrap A -> B (down from box3, across, down into box4)
yb = ROW_A - 3.5
arrow(XS[2] + W / 2, ROW_A, XS[2] + W / 2, yb)
ax.add_patch(FancyArrowPatch((XS[2] + W / 2, yb), (XS[0] + W / 2, yb),
             arrowstyle="-", linewidth=2.0, color=INK, zorder=1))
arrow(XS[0] + W / 2, yb, XS[0] + W / 2, ROW_B + H)
# Arrows — Row B
arrow(XS[0] + W, ROW_B + H / 2, XS[1], ROW_B + H / 2)
arrow(XS[1] + W, ROW_B + H / 2, XS[2], ROW_B + H / 2)
# Missing: box6 -> PRISM/GDSC (dashed red)
arrow(XS[1] + W / 2, ROW_B, XS[1] + W / 2, ROW_C + H, color=RED, dashed=True)
arrow(XS[2] + W / 2, ROW_B, XS[2] + W / 2, ROW_C + H, color=RED, dashed=True)

ax.text(99.5, 2.0, "* 190 case-insensitive (audit) vs 180 stricter pipeline normalization",
        ha="right", va="bottom", fontsize=8, color=GREY, style="italic")

fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
print(f"wrote {OUT}")
