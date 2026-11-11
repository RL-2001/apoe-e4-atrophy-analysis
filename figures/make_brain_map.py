#!/usr/bin/env python3
###############################################################################
#  Statistical brain map of the APOE e4-dose atrophy effect  (Destrieux / nilearn)
#
#  Paints each Destrieux cortical region onto an inflated brain surface, coloured
#  by how much faster e4/e4 homozygotes lose it (Excess atrophy, SD/year).
#
#  KEY ADVANTAGE: nilearn DOWNLOADS the Destrieux surface atlas itself, so you do
#  NOT need FreeSurfer or any local atlas files - just an internet connection once.
#
#  Run in your own Python. Needs: nilearn, matplotlib, numpy, pandas.
#     pip install nilearn matplotlib numpy pandas
###############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nilearn import datasets, surface, plotting

# ---- 1. Load your effect sizes -------------------------------------------
eff = pd.read_csv("nilearn_cortical_effects.csv")   # hemi, destrieux_name, effect, FDR_unadj, ...

# OPTIONAL: mask non-significant regions (set their effect to NaN -> shown grey)
MASK_NONSIG = True
if MASK_NONSIG:
    eff.loc[eff["FDR_unadj"] >= 0.05, "effect"] = np.nan

# ---- 2. Fetch the Destrieux surface atlas + a surface mesh ---------------
# nilearn downloads these the first time (cached afterwards).
destrieux = datasets.fetch_atlas_surf_destrieux()
fsaverage = datasets.fetch_surf_fsaverage("fsaverage")   # inflated meshes + sulcal maps

# The atlas gives, per hemisphere, an array 'map_left'/'map_right' assigning each
# vertex an integer index, and a 'labels' list (byte strings) for those indices.
labels = [l.decode() if isinstance(l, bytes) else l for l in destrieux["labels"]]
# labels[0] = 'Unknown', labels[1] = 'Medial_wall', then the 74 regions.

def build_vertex_values(hemi, vertex_map):
    """Return a per-vertex array of effect sizes for one hemisphere."""
    sub = eff[eff["hemi"] == hemi].set_index("destrieux_name")["effect"].to_dict()
    vals = np.full(vertex_map.shape, np.nan, dtype=float)
    for idx, lab in enumerate(labels):
        if lab in ("Unknown", "Medial_wall"):
            continue
        if lab in sub and not pd.isna(sub[lab]):
            vals[vertex_map == idx] = sub[lab]
    return vals

vals_left  = build_vertex_values("left",  destrieux["map_left"])
vals_right = build_vertex_values("right", destrieux["map_right"])

# ---- 3. Sanity check the matching (IMPORTANT) ----------------------------
atlas_names = set(l for l in labels if l not in ("Unknown", "Medial_wall"))
our_names   = set(eff["destrieux_name"])
missing = our_names - atlas_names
if missing:
    print("WARNING - these names did not match the atlas (fix before trusting the map):")
    for m in sorted(missing): print("   ", m)
else:
    print("All region names matched the atlas cleanly.")

# ---- 4. Plot: lateral + medial views, both hemispheres -------------------
vmax = np.nanmax(np.abs(np.concatenate([vals_left, vals_right])))
fig, axes = plt.subplots(2, 2, figsize=(11, 9), subplot_kw={"projection": "3d"})
views = [("left","lateral",vals_left,0,0), ("left","medial",vals_left,0,1),
         ("right","lateral",vals_right,1,0), ("right","medial",vals_right,1,1)]
for hemi, view, vals, r, c in views:
    mesh   = fsaverage[f"infl_{hemi}"]
    bg     = fsaverage[f"sulc_{hemi}"]
    plotting.plot_surf_roi(
        mesh, roi_map=vals, hemi=hemi, view=view,
        bg_map=bg, bg_on_data=True, darkness=0.5,
        cmap="hot_r", vmin=0, vmax=vmax,
        axes=axes[r, c], colorbar=False)
    axes[r, c].set_title(f"{hemi} {view}", fontsize=11)

# one shared colourbar
sm = plt.cm.ScalarMappable(cmap="hot_r",
        norm=plt.Normalize(vmin=0, vmax=vmax))
cbar = fig.colorbar(sm, ax=axes, shrink=0.5, pad=0.04)
cbar.set_label("e4/e4 excess atrophy (SD/year)", fontsize=11)
fig.suptitle("APOE \u03b54-dose acceleration of cortical atrophy", fontsize=14, fontweight="bold")

# ---- 5. Save at publication resolution -----------------------------------
fig.savefig("figure_brain_map.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig("figure_brain_map.pdf", bbox_inches="tight")   # vector, best for journals
print("Saved figure_brain_map.png and .pdf")

###############################################################################
#  If step 3 prints WARNINGs about unmatched names, copy them to me and I'll
#  give you the exact rename to fix them. Otherwise you're done.
###############################################################################
