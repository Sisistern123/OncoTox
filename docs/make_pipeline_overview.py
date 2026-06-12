"""Generate docs/pipeline_overview.png — a status overview of the OncoTox pipeline.

Layout mirrors the eight thematic step files in docs/steps/:
  01-05 = done (green)   ·   06-08 = planned placeholders (red, dashed)

Green  = done / on-plan
Amber  = done but an addition beyond the plan
Red    = not started (still missing / horizon)

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

fig, ax = plt.subplots(figsize=(17.0, 9.0))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# Grid — 4 columns x 2 rows of step boxes, plus a full-width amber band at the bottom.
W, H = 21.5, 20
XS = [3.5, 28.0, 52.5, 77.0]
ROW_A, ROW_B = 58, 28


def box(x, y, w, h, title, lines, edge, fill, title_color=None, dashed=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=1.4",
        linewidth=2.2, edgecolor=edge, facecolor=fill,
        linestyle="--" if dashed else "-", mutation_aspect=1.0, zorder=2,
    ))
    ax.text(x + w / 2, y + h - 3.0, title, ha="center", va="top",
            fontsize=10.5, fontweight="bold", color=title_color or INK, zorder=3)
    ax.text(x + w / 2, y + h - 7.6, "\n".join(lines), ha="center", va="top",
            fontsize=8.0, color=INK, zorder=3)


def arrow(x1, y1, x2, y2, color=INK, dashed=False):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
        linewidth=2.0, color=color, linestyle="--" if dashed else "-", zorder=1,
    ))


# Title + legend
ax.text(50, 98, "OncoTox Pipeline — Status Overview", ha="center", va="top",
        fontsize=17, fontweight="bold", color=INK)
ax.text(50, 93.5, "as of 2026-06-12   ·   reference: project_planning_v2.pdf   ·   steps: docs/steps/",
        ha="center", va="top", fontsize=9.5, color=GREY)
handles = [
    mpatches.Patch(facecolor=GREEN_FILL, edgecolor=GREEN, label="Done / on-plan"),
    mpatches.Patch(facecolor=AMBER_FILL, edgecolor=AMBER, label="Addition beyond plan"),
    mpatches.Patch(facecolor=RED_FILL, edgecolor=RED, label="Not started (planned)"),
]
ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.885),
          ncol=3, fontsize=9.5, frameon=True, framealpha=0.9)

# Row A — steps 01-04 (done)
box(XS[0], ROW_A, W, H, "01 · Datasets & harmonization",
    ["SCP542 53,513 cells x 22,722 g", "CTRPv2 545 drugs (cpd_avg_pv)",
     "overlap 190* lines · BRD/DrugBank"], GREEN, GREEN_FILL)
box(XS[1], ROW_A, W, H, "02 · Preprocessing & embeddings",
    ["scGPT  X_scGPT = 512-d", "HVG-5000 (4,576/5,000 vocab)",
     "+ X_pca · UMAP Fig 3/4"], GREEN, GREEN_FILL)
box(XS[2], ROW_A, W, H, "03 · Model & training design",
    ["per-cell input -> viability", "supervised regression (masked)",
     "frozen scGPT prior · small MLP"], GREEN, GREEN_FILL)
box(XS[3], ROW_A, W, H, "04 · Single-task baseline",
    ["paclitaxel, leak-free split", "best scGPT val MSE 0.0336",
     "1 DB · 1 score · 1 drug"], GREEN, GREEN_FILL)

# Row B — step 05 (done) + steps 06-08 (missing)
box(XS[0], ROW_B, W, H, "05 · Multi-task (masked MSE)",
    ["K=545 CTRP drugs, split_ctrp", "scGPT beats baseline 142/545",
     "still 1 DB · 1 score (multi-drug)"], GREEN, GREEN_FILL)
box(XS[1], ROW_B, W, H, "06 · Cross-database  (MISSING)",
    ["CTRPv2 + PRISM + GDSC", "efficacy + toxicity heads",
     "the 'combine all' goal"], RED, RED_FILL, title_color=RED, dashed=True)
box(XS[2], ROW_B, W, H, "07 · XAI / interpretability  (MISSING)",
    ["feature importance -> drivers", "stretch goal", "not started"],
    RED, RED_FILL, title_color=RED, dashed=True)
box(XS[3], ROW_B, W, H, "08 · Foundation model  (HORIZON)",
    ["reusable pan-cancer FM", "fine-tune on clinical (binary)",
     "overarching main goal"], RED, RED_FILL, title_color=RED, dashed=True)

# Bottom amber band — additions beyond the written plan
BAND_Y, BAND_H = 5, 13
box(XS[0], BAND_Y, 94, BAND_H, "Additions beyond the written plan",
    ["HVG-5000 + all_genes variant   ·   random -> cell-line-grouped split (leak fix)   ·   "
     "per-drug-mean sanity baseline   ·   run versioning (runs/ + ledger)"],
    AMBER, AMBER_FILL, title_color=AMBER)

# Arrows — Row A (01 -> 02 -> 03 -> 04)
for i in range(3):
    arrow(XS[i] + W, ROW_A + H / 2, XS[i + 1], ROW_A + H / 2)
# wrap 04 -> 05 (down from box04, across, down into box05)
yb = ROW_A - 3.0
arrow(XS[3] + W / 2, ROW_A, XS[3] + W / 2, yb)
ax.add_patch(FancyArrowPatch((XS[3] + W / 2, yb), (XS[0] + W / 2, yb),
             arrowstyle="-", linewidth=2.0, color=INK, zorder=1))
arrow(XS[0] + W / 2, yb, XS[0] + W / 2, ROW_B + H)
# Arrows — Row B (05 -> 06 -> 07 -> 08), dashed red into the not-yet-built chain
for i in range(3):
    arrow(XS[i] + W, ROW_B + H / 2, XS[i + 1], ROW_B + H / 2, color=RED, dashed=True)

ax.text(99.5, 1.2, "* 190 case-insensitive (audit) vs 180 stricter pipeline normalization",
        ha="right", va="bottom", fontsize=8, color=GREY, style="italic")

fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
print(f"wrote {OUT}")
