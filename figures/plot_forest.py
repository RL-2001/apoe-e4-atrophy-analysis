#!/usr/bin/env python3
"""
Amyloid-stratified forest plot of the ε4/ε4-vs-non-carrier baseline difference.

Reads results/forest_amyloid.csv with columns:
    Region, Stratum, diff, lo, hi, n0, n2
where `diff` is the ε4/ε4 minus non-carrier difference (SD units) with
95% CI [lo, hi], within each amyloid Stratum (Amyloid+ / Amyloid-).

Output: figure_forest_amyloid.(png|pdf)
"""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE.parent / "results" / "forest_amyloid.csv")

strata = ["Amyloid+", "Amyloid-"]
colors = {"Amyloid+": "#C0392B", "Amyloid-": "#2C6FB0"}

plot_items, labels, y = [], [], []
ypos = 0
for s in strata:
    block = df[df["Stratum"] == s].sort_values("diff")
    for _, r in block.iterrows():
        plot_items.append((ypos, r, s)); labels.append(r["Region"]); y.append(ypos); ypos += 1
    ypos += 1  # gap between strata

fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(plot_items) + 2)))
for yp, r, s in plot_items:
    ax.plot([r["lo"], r["hi"]], [yp, yp], color=colors[s], lw=2, zorder=2)
    ax.plot(r["diff"], yp, "o", color=colors[s], ms=6, zorder=3)
ax.axvline(0, color="#555", lw=1, ls="--", zorder=1)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8); ax.invert_yaxis()
ax.set_xlabel("ε4/ε4 − non-carrier baseline difference (SD)", fontsize=11)
ax.set_title("APOE ε4/ε4 effect by amyloid status", fontsize=13, fontweight="bold")
handles = [plt.Line2D([0], [0], color=colors[s], marker="o", lw=2, label=s) for s in strata]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "figure_forest_amyloid.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(HERE / "figure_forest_amyloid.pdf", bbox_inches="tight")
print("Saved figure_forest_amyloid.png / .pdf")
