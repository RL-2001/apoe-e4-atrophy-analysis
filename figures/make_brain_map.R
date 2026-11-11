###############################################################################
#  Statistical brain map of the APOE e4-dose atrophy effect  (Destrieux / ggseg)
#  Colours each of the 148 cortical Destrieux regions by how much faster
#  e4/e4 homozygotes lose it (Excess atrophy, SD/year).
#
#  You run this in R in your own environment. It needs the internet once to
#  install packages. It does NOT need your MRI files - ggseg has the atlas
#  geometry built in; we only supply the per-region numbers.
###############################################################################

## ---- 1. Install / load packages (first time only) -------------------------
# install.packages("ggseg")          # base ggseg + built-in atlases
# install.packages("ggplot2")
# install.packages("dplyr")
# If the Destrieux atlas is not bundled in your ggseg version, also:
# install.packages("ggseg3d")        # sometimes ships extra atlas data
# (Destrieux is included in recent ggseg as `ggseg::destrieux`.)

library(ggseg)
library(ggplot2)
library(dplyr)

## ---- 2. Load your effect sizes -------------------------------------------
# Point this at the file I prepared for you.
dat <- read.csv("ggseg_cortical_effects.csv", stringsAsFactors = FALSE)
# columns: Region, ggseg_label, hemi, region_noheim, Atlas_ID, effect, FDR_unadj, FDR_amyloidAdj

## ---- 3. Inspect the atlas's own region labels (CRITICAL matching step) ----
# Pull the label vocabulary ggseg expects for Destrieux, so we can match to it.
atlas_labels <- as.data.frame(destrieux)      # the built-in atlas data frame
# Look at how ggseg names its regions/labels:
print(head(unique(atlas_labels$label), 20))
print(head(unique(atlas_labels$region), 20))
# -> compare these to dat$ggseg_label. They should line up; if the separators
#    differ (e.g. ggseg uses "_and_" vs "&", or "." vs "-"), adjust in step 4.

## ---- 4. Join your data to the atlas ---------------------------------------
# ggseg matches on a column called `label` (full, hemisphere-prefixed) OR
# `region` (name without hemisphere). We try `label` first.
plotdat <- dat %>%
  rename(label = ggseg_label) %>%
  mutate(label = tolower(label))    # ggseg labels are lower-case; safe to match

# OPTIONAL: mask non-significant regions to grey by setting their effect to NA
plotdat_sig <- plotdat %>%
  mutate(effect = ifelse(FDR_unadj < 0.05, effect, NA))

## ---- 5. Draw the map ------------------------------------------------------
# Diverging/sequential fill; here a white-to-red scale since all effects are
# "faster loss" (positive excess). Reverse or recentre if you prefer.
p <- ggseg(.data = plotdat_sig,
           atlas   = destrieux,
           mapping = aes(fill = effect),
           position = "stacked",       # lateral + medial stacked; or "dispersed"
           colour  = "grey30", size = 0.2) +
  scale_fill_gradient(low = "#FFE9E4", high = "#8B1A12",
                      na.value = "grey88",
                      name = "Excess atrophy\n(SD / year)") +
  theme_void(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "APOE \u03b54-dose acceleration of cortical atrophy")

print(p)

## ---- 6. Save at publication resolution ------------------------------------
ggsave("figure_brain_map.png", p, width = 9, height = 6, dpi = 400, bg = "white")
ggsave("figure_brain_map.pdf", p, width = 9, height = 6)          # vector, best for journals

###############################################################################
#  If step 3 shows the labels DON'T match (common), the fix is almost always
#  a string-format difference. Send me the printed atlas labels from step 3
#  and I'll give you the exact transform to line them up with your data.
###############################################################################
