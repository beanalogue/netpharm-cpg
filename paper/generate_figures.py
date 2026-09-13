"""
Journal figure generation for:
  'Reverse Network Pharmacology as a Molecular Evidence Layer for
   Traditional Medicine CPGs' — J Ethnopharmacology submission

Generates:
  Figure 1 — pipeline flowchart
  Figure 2 — disease pool AUROC bar chart
  Figure 3 — sensitivity analysis line plot

All figures saved as 300 DPI PNG + PDF (vector) in paper/figures/
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size":         8,
    "axes.labelsize":    9,
    "axes.titlesize":    9,
    "xtick.labelsize":   8,
    "ytick.labelsize":   8,
    "legend.fontsize":   8,
    "figure.dpi":        300,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth":   1.4,
    "pdf.fonttype":      42,   # embed fonts in PDF
    "ps.fonttype":       42,
})

OUT = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT, exist_ok=True)

# ── Colour palette (matches graphical abstract) ─────────────────────
C_SIG    = "#2C4A7C"   # navy  — HTN (significant)
C_TREND  = "#B87318"   # amber — INS (trend)
C_NULL   = "#8B9BB4"   # slate — DEM (null)
C_RANDOM = "#BDBDBD"   # light grey — chance line
C_TEXT   = "#1A1A2E"   # near-black
C_STEP   = "#1E3A5F"   # step box fill (dark navy)
C_STEP2  = "#2C7A50"   # step 5 (result box, green)
C_ARROW  = "#4A4A6A"

# ════════════════════════════════════════════════════════════════════
# Figure 1 — Pipeline flowchart
# ════════════════════════════════════════════════════════════════════
def make_figure1():
    fig, ax = plt.subplots(figsize=(7.0, 7.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    steps = [
        {
            "step": "Step 1",
            "title": "Disease Gene Targets",
            "detail": "OpenTargets Platform (MONDO ontology)\nscore ≥ 0.10  ·  MONDO IDs curated per indication",
            "y": 8.7,
            "color": C_STEP,
        },
        {
            "step": "Step 2",
            "title": "Herb Compound–Target Retrieval",
            "detail": "TCMSP database  ·  500 indexed herbs\nOB ≥ 30%  ·  DL ≥ 0.18  (Lipinski-based filters)",
            "y": 6.9,
            "color": C_STEP,
        },
        {
            "step": "Step 3",
            "title": "Hypergeometric Herb Scoring",
            "detail": "p_h = P(X ≥ k),  X ~ Hypergeom(N=20 000, |T_h|, |D|)\nAll 500 herbs scored → ranked list",
            "y": 5.1,
            "color": C_STEP,
        },
        {
            "step": "Step 4",
            "title": "CPG Formula Herb Pool",
            "detail": "Korean TM CPG 2021 editions\nChinese Pharmacopoeia 2020 standard compositions",
            "y": 3.3,
            "color": C_STEP,
        },
        {
            "step": "Step 5",
            "title": "Pool AUROC + Permutation Test",
            "detail": "AUROC via Mann–Whitney U  ·  n = 1,000 permutations, seed = 42\nBonferroni correction  α = 0.05 / 3 = 0.0167",
            "y": 1.5,
            "color": C_STEP2,
        },
    ]

    box_h = 1.25
    box_x = 0.8
    box_w = 8.4

    for i, s in enumerate(steps):
        # main box
        box = FancyBboxPatch(
            (box_x, s["y"] - box_h / 2),
            box_w, box_h,
            boxstyle="round,pad=0.05",
            linewidth=1.0,
            edgecolor=s["color"],
            facecolor=s["color"] + "20",   # 12% opacity
        )
        ax.add_patch(box)

        # step number badge
        badge = FancyBboxPatch(
            (box_x, s["y"] - box_h / 2),
            1.05, box_h,
            boxstyle="round,pad=0.05",
            linewidth=0,
            facecolor=s["color"],
        )
        ax.add_patch(badge)
        ax.text(
            box_x + 0.525, s["y"],
            s["step"],
            ha="center", va="center",
            color="white", fontsize=8, fontweight="bold",
        )

        # title
        ax.text(
            box_x + 1.35, s["y"] + 0.17,
            s["title"],
            ha="left", va="center",
            color=s["color"], fontsize=9, fontweight="bold",
        )
        # detail
        ax.text(
            box_x + 1.35, s["y"] - 0.25,
            s["detail"],
            ha="left", va="center",
            color=C_TEXT, fontsize=7.5, linespacing=1.4,
        )

        # arrow to next step
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(box_x + box_w / 2, steps[i + 1]["y"] + box_h / 2 + 0.01),
                xytext=(box_x + box_w / 2, s["y"] - box_h / 2 - 0.01),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=C_ARROW,
                    lw=1.2,
                    mutation_scale=12,
                ),
            )

    # Side note: exclusions
    ax.annotate(
        "Exclusions documented\nat each step →\nconservative bias",
        xy=(box_x + box_w + 0.1, 5.1),
        xytext=(box_x + box_w + 0.15, 5.1),
        ha="left", va="center",
        fontsize=6.5, color="#666666",
        annotation_clip=False,
    )

    ax.set_title(
        "Figure 1. Reverse network pharmacology pipeline for molecular evidence evaluation of traditional medicine CPGs",
        fontsize=8.5, pad=10, loc="left", color=C_TEXT,
    )

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"figure1.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("  Figure 1 saved.")


# ════════════════════════════════════════════════════════════════════
# Figure 2 — Disease pool AUROC bar chart
# ════════════════════════════════════════════════════════════════════
def make_figure2():
    diseases  = ["Essential\nhypertension", "Insomnia\ndisorder", "Dementia /\ncognitive impairment"]
    aurocs    = [0.655, 0.596, 0.571]
    perm_ps   = [0.001, 0.054, 0.146]   # HTN perm_p < 0.001 (0/1000); plotted at 0.001 for annotation only
    colors    = [C_SIG, C_TREND, C_NULL]
    labels    = ["Significant\n(p < 0.001)", "Trend\n(p = 0.054)", "Null\n(p = 0.146)"]

    # Approximate 95% CI of null distribution (analytical, Mann-Whitney)
    # n_formula / n_total from paper: HTN 27/450, INS 21/412, DEM 20/456
    # SE = sqrt(n1*n2*(n1+n2+1)/12) / (n1*n2)
    n1s = [27, 21, 20]
    n2s = [423, 391, 436]
    null_cis = []
    for n1, n2 in zip(n1s, n2s):
        se = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12) / (n1 * n2)
        null_cis.append(1.96 * se)   # half-width of 95% CI

    fig, ax = plt.subplots(figsize=(3.5, 4.0))

    x = np.arange(len(diseases))
    bars = ax.bar(
        x, aurocs,
        color=colors, width=0.55,
        linewidth=0.8, edgecolor="white",
        zorder=3,
    )

    # Null distribution 95% CI as shaded band
    ax.axhspan(
        0.5 - null_cis[0], 0.5 + null_cis[0],   # use HTN CI (similar across diseases)
        color=C_RANDOM, alpha=0.25, zorder=1, label="Null 95% CI (approx.)",
    )

    # Reference line at 0.5
    ax.axhline(0.5, color=C_RANDOM, lw=1.2, ls="--", zorder=2)
    ax.text(2.35, 0.5 + 0.003, "chance", va="bottom", ha="right",
            fontsize=7, color="#888888")

    # Observed AUROC dots (on top of bars)
    ax.scatter(x, aurocs, color="white", s=22, zorder=5, linewidths=0.8,
               edgecolors=[c for c in colors])

    # Significance marker
    sig_y = aurocs[0] + 0.018
    ax.text(x[0], sig_y, "*", ha="center", va="bottom",
            fontsize=14, color=C_SIG, fontweight="bold")
    ax.text(x[0], sig_y + 0.030, "p < 0.001", ha="center", va="bottom",
            fontsize=6.5, color=C_SIG)

    # p-value annotations for other diseases
    for i in [1, 2]:
        ax.text(x[i], aurocs[i] + 0.010, f"p = {perm_ps[i]:.3f}",
                ha="center", va="bottom", fontsize=6.5, color=colors[i])

    # AUROC value inside each bar
    for i, (xi, auroc) in enumerate(zip(x, aurocs)):
        ax.text(xi, auroc - 0.015, f"{auroc:.3f}",
                ha="center", va="top", fontsize=8, color="white",
                fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(diseases, fontsize=8)
    ax.set_ylabel("Pool AUROC", fontsize=9)
    ax.set_ylim(0.44, 0.74)
    ax.set_yticks([0.5, 0.55, 0.60, 0.65, 0.70])
    ax.yaxis.grid(True, lw=0.5, color="#E0E0E0", zorder=0)
    ax.set_axisbelow(True)

    # Legend for significance categories
    legend_elements = [
        mpatches.Patch(facecolor=C_SIG,   label="Significant (Bonferroni p < 0.0167)"),
        mpatches.Patch(facecolor=C_TREND, label="Positive trend (p = 0.054)"),
        mpatches.Patch(facecolor=C_NULL,  label="Null (p = 0.146)"),
        mpatches.Patch(facecolor=C_RANDOM, alpha=0.4, label="Null 95% CI (approx.)"),
    ]
    ax.legend(handles=legend_elements, fontsize=6.5, loc="upper right",
              framealpha=0.9, edgecolor="#CCCCCC")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title(
        "Figure 2. Disease pool AUROC across three Korean medicine CPG indications",
        fontsize=8, loc="left", pad=8, color=C_TEXT,
    )

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"figure2.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("  Figure 2 saved.")


# ════════════════════════════════════════════════════════════════════
# Figure 3 — Sensitivity analysis line plot
# ════════════════════════════════════════════════════════════════════
def make_figure3():
    scenarios = list(range(1, 12))
    labels_x = [
        "Primary\n(OT≥0.10,\nk≥1)",
        "OT≥0.05†\nk≥1",
        "OT≥0.20\nk≥1",
        "OT≥0.30\nk≥1",
        "OB≥20%‡\nk≥1",
        "OB≥40%‡\nk≥1",
        "OT≥0.10\nk≥2",
        "OT≥0.10\nk≥3",
        "OT≥0.20\nk≥2",
        "OT≥0.30\nk≥2",
        "SVM≥0.7§\nk≥1",
    ]

    # Scenarios 1-10: OT/OB/overlap parameter sweep
    # Scenario 11: TCMSP SVM/RF confidence threshold >= 0.7
    htn_auroc = [0.655, 0.655, 0.659, 0.600, 0.655, 0.655, 0.640, 0.623, 0.644, 0.600, 0.662]
    htn_sigp  = [0.0005, 0.0005, 0.0005, 0.030, 0.0005, 0.0005, 0.006, 0.014, 0.003, 0.046, 0.003]
    ins_auroc = [0.596, 0.596, 0.587, 0.599, 0.596, 0.596, 0.607, 0.636, 0.622, 0.608, 0.618]
    ins_sigp  = [0.054, 0.054, 0.081, 0.059, 0.054, 0.054, 0.048, 0.027, 0.026, 0.041, 0.039]
    dem_auroc = [0.571, 0.571, 0.571, 0.571, 0.571, 0.571, 0.557, 0.574, 0.557, 0.557, 0.506]
    dem_sigp  = [0.146, 0.146, 0.146, 0.146, 0.146, 0.146, 0.188, 0.136, 0.188, 0.188, 0.451]

    BONF = 0.0167
    x = np.arange(1, 12)

    fig, ax = plt.subplots(figsize=(8.0, 4.0))

    htn_sig_mask = [p < BONF for p in htn_sigp]

    # Plot lines
    ax.plot(x, htn_auroc, color=C_SIG, marker="o", ms=6, lw=1.5,
            label="Essential hypertension", zorder=4)
    ax.plot(x, ins_auroc, color=C_TREND, marker="^", ms=6, lw=1.5,
            label="Insomnia disorder", zorder=4)
    ax.plot(x, dem_auroc, color=C_NULL, marker="s", ms=6, lw=1.5,
            label="Dementia", zorder=4)

    # Mark duplicate scenarios (2, 5, 6) with open markers
    for idx in [1, 4, 5]:
        ax.plot(x[idx], htn_auroc[idx], "o", ms=8, mfc="white",
                mec=C_SIG, mew=1.5, zorder=5)
        ax.plot(x[idx], ins_auroc[idx], "^", ms=8, mfc="white",
                mec=C_TREND, mew=1.5, zorder=5)
        ax.plot(x[idx], dem_auroc[idx], "s", ms=8, mfc="white",
                mec=C_NULL, mew=1.5, zorder=5)

    # Scenario 11: diamond marker to distinguish from main scenarios
    ax.plot(x[10], htn_auroc[10], "D", ms=7, color=C_SIG,   zorder=6)
    ax.plot(x[10], ins_auroc[10], "D", ms=7, color=C_TREND, zorder=6)
    ax.plot(x[10], dem_auroc[10], "D", ms=7, color=C_NULL,  zorder=6)

    # Vertical separator before Scenario 11
    ax.axvline(10.5, color="#CCCCCC", lw=0.8, ls=":", zorder=1)

    # Mark significant HTN scenarios
    for i, (auroc, sig) in enumerate(zip(htn_auroc, htn_sig_mask)):
        if sig:
            ax.text(x[i], auroc + 0.012, "*", ha="center", va="bottom",
                    fontsize=10, color=C_SIG, fontweight="bold")

    # Reference lines
    ax.axhline(0.5,   color=C_RANDOM, lw=1.0, ls="--", zorder=1)
    ax.axhline(0.655, color=C_SIG,    lw=0.8, ls=":",  alpha=0.5, zorder=1)

    ax.text(11.35, 0.502, "chance\n(0.5)", va="bottom", ha="left",
            fontsize=6.5, color="#888888", clip_on=False)
    ax.text(11.35, 0.656, "primary\n(0.655)", va="bottom", ha="left",
            fontsize=6.5, color=C_SIG, alpha=0.7, clip_on=False)

    # Annotations
    ax.annotate("‡ OB cached\nat 30%",
        xy=(5.5, 0.542), xytext=(5.2, 0.519), fontsize=6.5, color="#888888",
        arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=0.8), ha="center")
    ax.annotate("† OT page-\ncapped",
        xy=(2, 0.596), xytext=(1.55, 0.519), fontsize=6.5, color="#888888",
        arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=0.8), ha="center")
    ax.annotate("§ SVM/RF\nconfidence",
        xy=(11, 0.506), xytext=(10.6, 0.520), fontsize=6.5, color="#888888",
        arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=0.8), ha="center")

    ax.set_xticks(x)
    ax.set_xticklabels(labels_x, fontsize=6.5, ha="center")
    ax.set_ylabel("Pool AUROC", fontsize=9)
    ax.set_xlim(0.3, 11.7)
    ax.set_ylim(0.47, 0.70)
    ax.set_yticks([0.50, 0.55, 0.60, 0.65])
    ax.yaxis.grid(True, lw=0.5, color="#E0E0E0", zorder=0)
    ax.set_axisbelow(True)

    # Shaded columns for HTN non-significant scenarios (4, 10)
    for idx in [3, 9]:
        ax.axvspan(x[idx] - 0.35, x[idx] + 0.35, color="#F0F0F0", zorder=0, alpha=0.9)

    legend_elements = [
        Line2D([0], [0], color=C_SIG,   marker="o", ms=5, label="Essential hypertension"),
        Line2D([0], [0], color=C_TREND, marker="^", ms=5, label="Insomnia disorder"),
        Line2D([0], [0], color=C_NULL,  marker="s", ms=5, label="Dementia"),
        Line2D([0], [0], color="none",  marker="o", ms=6, mfc="white",
               mec="#666666", mew=1.2, label="Duplicate of primary (not independent)"),
        Line2D([0], [0], color="none",  marker="D", ms=6, mfc="#666666",
               label="Scenario 11: SVM/RF ≥ 0.7 confidence filter"),
        Line2D([0], [0], color="none",  marker="*", ms=8, mfc=C_SIG, linestyle="None",
               label="Significant (Bonferroni p < 0.0167)"),
        mpatches.Patch(facecolor="#F0F0F0", label="HTN non-significant (shaded)"),
    ]
    ax.legend(handles=legend_elements, fontsize=6.2, loc="upper left",
              framealpha=0.9, edgecolor="#CCCCCC", ncol=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title(
        "Figure 3. Sensitivity analysis — pool AUROC across 11 scenarios (7 independent + Scenario 11)",
        fontsize=8, loc="left", pad=8, color=C_TEXT,
    )

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"figure3.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("  Figure 3 saved.")


# ── Main ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating figures...")
    make_figure1()
    make_figure2()
    make_figure3()
    print(f"Done. Figures saved to: {OUT}/")
