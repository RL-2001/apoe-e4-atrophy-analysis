#!/usr/bin/env python3
"""
Subcortical panel: APOE e4/e4 excess atrophy for subcortical structures.

Reads results/subcortical_panel_data.csv with columns:
    structure, hemi, effect, FDR, sig

Grouped horizontal bars (left/right hemisphere) per structure, shaded by
FDR significance.

Output: figure_subcortical.(png|pdf)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE.parent / "results" / "subcortical_panel_data.csv")

structures = list(dict.fromkeys(df["structure"]))  # preserve order
hemis = ["Left", "Right"]
hcol = {"Left": "#2C6FB0", "Right": "#C0392B"}

y = np.arange(len(structures))
h = 0.38
fig, ax = plt.subplots(figsize=(7.5, max(4, 0.5 * len(structures) + 1.5)))
for i, hemi in enumerate(hemis):
    sub = df[df["hemi"] == hemi].set_index("structure")
    vals = [sub["effect"].get(s, np.nan) for s in structures]
    sig = [bool(sub["sig"].get(s, False)) for s in structures]
    offset = (i - 0.5) * h
    bars = ax.barh(y + offset, vals, height=h, label=hemi, color=hcol[hemi])
    # per-bar transparency encodes FDR significance (faded = not significant)
    for b, s in zip(bars, sig):
        b.set_alpha(1.0 if s else 0.4)

ax.set_yticks(y)
ax.set_yticklabels(structures, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("e4/e4 excess atrophy (SD/year)", fontsize=11)
ax.set_title("Subcortical APOE ε4-dose effects", fontsize=13, fontweight="bold")
ax.legend(loc="lower right", frameon=False, fontsize=9,
          title="(faded = not FDR-significant)", title_fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "figure_subcortical.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(HERE / "figure_subcortical.pdf", bbox_inches="tight")
print("Saved figure_subcortical.png / .pdf")
