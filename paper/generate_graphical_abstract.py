"""
Graphical abstract — J Ethnopharmacology submission
Output: paper/figures/graphical_abstract.{png,pdf}
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size":       9,
    "figure.dpi":      300,
    "savefig.dpi":     300,
    "savefig.bbox":    "tight",
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})

OUT = os.path.join(os.path.dirname(__file__), "figures")

C_SIG    = "#2C4A7C"
C_TREND  = "#B87318"
C_NULL   = "#8B9BB4"
C_STEP   = "#1E3A5F"
C_GSEA   = "#2C7A50"
C_BG     = "#F7F9FC"
C_TEXT   = "#1A1A2E"
C_GRID   = "#DDE3ED"
C_RAND   = "#9E9E9E"
C_ACCENT = "#C0392B"

fig = plt.figure(figsize=(9.5, 4.4))
fig.patch.set_facecolor(C_BG)

fig.text(0.5, 0.972,
         "Rank-based enrichment analysis of Korean traditional medicine "
         "clinical practice guidelines",
         ha="center", va="top", fontsize=10, fontweight="bold",
         color=C_TEXT, fontstyle="italic")

gs = gridspec.GridSpec(
    1, 3, figure=fig,
    left=0.03, right=0.97, top=0.87, bottom=0.08,
    wspace=0.12
)

# ════════════════════════════════════════════════════════════════════
# Panel A — Pipeline
# ════════════════════════════════════════════════════════════════════
axA = fig.add_subplot(gs[0])
axA.set_xlim(0, 1); axA.set_ylim(0, 1)
axA.axis("off"); axA.set_facecolor(C_BG)

axA.text(0.05, 0.97, "A  Pipeline", fontsize=9, fontweight="bold",
         color=C_TEXT, va="top")

steps = [
    (0.81, ["Disease gene targets", "(OpenTargets)"],           "#DDEEFF"),
    (0.57, ["Herb NP scoring", "(TCMSP · 500 herbs)"],         "#DDEEFF"),
    (0.33, ["CPG formula pool", "(Korean TM CPG 2021)"],       "#DDEEFF"),
    (0.09, ["AUROC enrichment test", "(permutation n=1,000)"], "#D6EFE0"),
]

box_w, box_h = 0.74, 0.13
for cy, lines, fc in steps:
    rect = FancyBboxPatch(
        (0.5 - box_w/2, cy - box_h/2), box_w, box_h,
        boxstyle="round,pad=0.015", linewidth=0.8,
        edgecolor=C_STEP, facecolor=fc,
        transform=axA.transAxes, zorder=3)
    axA.add_patch(rect)
    axA.text(0.50, cy, "\n".join(lines), ha="center", va="center",
             fontsize=7.5, color=C_TEXT, linespacing=1.35,
             transform=axA.transAxes, zorder=4,
             fontweight="bold" if cy == 0.09 else "normal")

for y_from, y_to in [(0.81-box_h/2, 0.57+box_h/2),
                     (0.57-box_h/2, 0.33+box_h/2),
                     (0.33-box_h/2, 0.09+box_h/2)]:
    axA.annotate("", xy=(0.50, y_to), xytext=(0.50, y_from),
                 arrowprops=dict(arrowstyle="-|>", color=C_STEP,
                                 lw=1.0, mutation_scale=9),
                 xycoords="axes fraction", textcoords="axes fraction")

for cy, num in [(0.81,"①"),(0.57,"②"),(0.33,"③"),(0.09,"④")]:
    axA.text(0.04, cy, num, ha="center", va="center",
             fontsize=8.5, color=C_STEP, fontweight="bold",
             transform=axA.transAxes)

# ════════════════════════════════════════════════════════════════════
# Panel B — GSEA analogy (waterfall bars kept narrow; label inside)
# ════════════════════════════════════════════════════════════════════
axB = fig.add_subplot(gs[1])
axB.set_facecolor(C_BG); axB.axis("off")
axB.set_xlim(0, 1); axB.set_ylim(0, 1)

axB.text(0.05, 0.97, "B  GSEA-style enrichment", fontsize=9,
         fontweight="bold", color=C_TEXT, va="top")

n_total = 38
np.random.seed(42)
formula_ranks = sorted(np.random.choice(range(0, 12), 8, replace=False))

bar_h    = 0.017
bar_top  = 0.86
bar_left = 0.22
bar_wmax = 0.38   # narrower — keeps bracket inside panel
y_gap    = 0.002
scores   = np.exp(-np.linspace(0, 3.5, n_total)) * 0.88 + 0.02

for i in range(n_total):
    y = bar_top - i * (bar_h + y_gap)
    if y < 0.10:
        break
    w = scores[i] * bar_wmax
    color = C_SIG if i in formula_ranks else C_GRID
    alpha = 1.0 if i in formula_ranks else 0.55
    axB.add_patch(mpatches.Rectangle((bar_left, y), w, bar_h,
                                     color=color, alpha=alpha, zorder=3))

# Top / bottom labels — short, inside panel
axB.text(bar_left, bar_top + bar_h + 0.01, "← High NP relevance",
         ha="left", va="bottom", fontsize=6.2, color=C_TEXT, style="italic")
axB.text(bar_left, 0.13, "← Low NP relevance",
         ha="left", va="top", fontsize=6.2, color=C_TEXT, style="italic")

# Bracket: positioned to stay inside panel (max x ≈ 0.75)
fr_y_top = bar_top - formula_ranks[0]*(bar_h+y_gap) + bar_h
fr_y_bot = bar_top - formula_ranks[-1]*(bar_h+y_gap)
bx = bar_left + bar_wmax + 0.03   # ≈ 0.63

# Draw bracket manually (two horizontal ticks + vertical bar)
mid_y = (fr_y_top + fr_y_bot) / 2
for tick_y in [fr_y_top, fr_y_bot]:
    axB.plot([bx, bx + 0.04], [tick_y, tick_y],
             color=C_SIG, lw=1.2, transform=axB.transAxes, clip_on=False)
axB.plot([bx + 0.04, bx + 0.04], [fr_y_bot, fr_y_top],
         color=C_SIG, lw=1.2, transform=axB.transAxes, clip_on=False)

# Label inside bracket area — max x ≈ 0.92
axB.text(bx + 0.07, mid_y,
         "CPG\nformula\nherbs",
         ha="left", va="center", fontsize=6.5, color=C_SIG,
         fontweight="bold", transform=axB.transAxes, clip_on=False)

# AUROC callout at bottom
axB.text(0.50, 0.04,
         "AUROC  =  enrichment statistic\n(analogous to GSEA enrichment score)",
         ha="center", va="bottom", fontsize=6.8, color=C_GSEA,
         fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.28", facecolor="#D6EFE0",
                   edgecolor=C_GSEA, linewidth=0.8),
         transform=axB.transAxes)

# ════════════════════════════════════════════════════════════════════
# Panel C — Results
# ════════════════════════════════════════════════════════════════════
axC = fig.add_subplot(gs[2])
axC.set_facecolor(C_BG)

# Panel label inside, top-left
axC.text(0.03, 0.97, "C  Results", fontsize=9, fontweight="bold",
         color=C_TEXT, va="top", transform=axC.transAxes)

diseases  = ["Essential\nhypertension", "Insomnia\ndisorder", "Dementia"]
aurocs    = [0.655, 0.596, 0.571]
colors    = [C_SIG, C_TREND, C_NULL]
pvals     = ["p < 0.001  ***", "p = 0.054", "p = 0.146"]
sig_flags = [True, False, False]
y_pos     = [2, 1, 0]

for y, auroc, col, pv, sig in zip(y_pos, aurocs, colors, pvals, sig_flags):
    axC.barh(y, auroc - 0.5, left=0.5, height=0.38,
             color=col, alpha=0.88, zorder=3)
    axC.text(auroc + 0.005, y + 0.11, f"{auroc:.3f}",
             va="center", ha="left", fontsize=9,
             fontweight="bold", color=col)
    axC.text(auroc + 0.005, y - 0.18, pv,
             va="center", ha="left", fontsize=7,
             color=C_ACCENT if sig else C_NULL,
             fontweight="bold" if sig else "normal")

# Random chance line — label below the line at left to avoid overlap
axC.axvline(0.5, color=C_RAND, lw=1.0, ls="--", zorder=2)
axC.text(0.497, 1.62, "Random\nchance",    # between INS and DEM bars
         ha="right", va="center", fontsize=6.3, color=C_RAND,
         style="italic")

axC.set_yticks(y_pos)
axC.set_yticklabels(diseases, fontsize=8, color=C_TEXT)
axC.set_xlim(0.44, 0.80)
axC.set_ylim(-0.65, 2.78)
axC.set_xlabel("Pool AUROC", fontsize=8.5, color=C_TEXT)
axC.spines[["top", "right", "left"]].set_visible(False)
axC.spines["bottom"].set_color(C_GRID)
axC.tick_params(axis="x", colors=C_TEXT, labelsize=7.5)
axC.tick_params(axis="y", left=False)
axC.xaxis.grid(True, color=C_GRID, lw=0.6, zorder=1)
axC.set_axisbelow(True)
axC.set_facecolor(C_BG)

# Sensitivity footnote at very bottom
axC.text(0.50, -0.60,
         "HTN: significant in 5 of 7 independent\n"
         "sensitivity scenarios (AUROC 0.600–0.659)",
         ha="left", va="bottom", fontsize=6.8, color=C_SIG,
         style="italic", transform=axC.get_yaxis_transform())

# ── Save ─────────────────────────────────────────────────────────────
for ext in ("png", "pdf"):
    path = os.path.join(OUT, f"graphical_abstract.{ext}")
    fig.savefig(path, dpi=300, facecolor=C_BG)
    print(f"Saved: {path}")

plt.close(fig)
