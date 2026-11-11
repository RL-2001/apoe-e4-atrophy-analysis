#!/usr/bin/env python3
"""
Top-N regions by APOE e4/e4 excess longitudinal atrophy.

Reads results/BRAINMAP_effect_sizes.csv with columns:
    Atlas_ID, Region, Slope_0alleles, Slope_1allele, Slope_2alleles,
    Excess_atrophy_SDyr, FDR_unadj, FDR_amyloidAdj

Horizontal bar chart of the N regions with the largest excess atrophy rate
(e4/e4 minus non-carrier slope, SD/year). Bars are shaded by whether the
effect survives FDR after amyloid adjustment.

Output: figure_top_regions.(png|pdf)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

TOP_N = 50
HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE.parent / "results" / "BRAINMAP_effect_sizes.csv")

df = df.dropna(subset=["Excess_atrophy_SDyr"])
top = df.sort_values("Excess_atrophy_SDyr", ascending=False).head(TOP_N).iloc[::-1]

# colour: survives amyloid-adjusted FDR vs. unadjusted only
def shade(row):
    if pd.notna(row.get("FDR_amyloidAdj")) and row["FDR_amyloidAdj"] < 0.05:
        return "#C0392B"      # robust to amyloid adjustment
    if pd.notna(row.get("FDR_unadj")) and row["FDR_unadj"] < 0.05:
        return "#E79A86"      # significant unadjusted only
    return "#BBBBBB"          # n.s.

colors = [shade(r) for _, r in top.iterrows()]

fig, ax = plt.subplots(figsize=(8.5, max(5, 0.26 * len(top) + 1.5)))
ax.barh(np.arange(len(top)), top["Excess_atrophy_SDyr"], color=colors)
ax.set_yticks(np.arange(len(top)))
ax.set_yticklabels(top["Region"], fontsize=7)
ax.set_xlabel("e4/e4 excess atrophy (SD/year)", fontsize=11)
ax.set_title(f"Top {TOP_N} regions by APOE ε4-dose–accelerated atrophy",
             fontsize=13, fontweight="bold")

legend = [
    plt.Rectangle((0, 0), 1, 1, color="#C0392B", label="FDR-sig, amyloid-adjusted"),
    plt.Rectangle((0, 0), 1, 1, color="#E79A86", label="FDR-sig, unadjusted only"),
    plt.Rectangle((0, 0), 1, 1, color="#BBBBBB", label="n.s."),
]
ax.legend(handles=legend, loc="lower right", frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "figure_top_regions.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(HERE / "figure_top_regions.pdf", bbox_inches="tight")
print(f"Saved figure_top_regions.png / .pdf (top {TOP_N})")
