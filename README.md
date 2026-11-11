# APOE ε4 dose and longitudinal brain atrophy — analysis & figures

Post-processing statistical analysis and figure-generation code for the study
of APOE ε4 allele dose and regional brain atrophy in the ADNI cohort
(deep-learning MRI parcellation; Lispi, Banti & Chen).

This repository covers the **downstream** half of the project: starting from the
ComBat-harmonized regional volumes, it produces the cross-sectional and
longitudinal statistics and all the figures reported in the paper. Effects are
reported in **SD/year** (regional volumes residualized in the cognitively-normal
reference sample and z-scored).

The **upstream** half : MRI preprocessing, DL+DiReCt parcellation, volumetric
extraction and ComBat harmonization (across scanner manufacturer), lives in a
separate repository:
<https://github.com/Mappo23/mri-preprocess-inference-pipeline>

## Validation

Running this pipeline on the upstream harmonized files reproduces the published
results:

- **Cohort matches Table 1:** 2,304 participants — ε4 dose 0 = 1,219, dose 1 =
  816, dose 2 = 269.
- **Longitudinal:** 169 / 179 regions show a dose-dependent difference in
  atrophy rate (one-way ANOVA across dose, Benjamini-Hochberg FDR < 0.05).
- **Distributed topography:** the largest ε4/ε4 excess-atrophy effects are in
  amygdala and hippocampus **and** lingual/occipitotemporal, bilateral
  middle-temporal, orbital and posterior-cingulate cortex — i.e. well beyond the
  medial temporal lobe.
- **Stage specificity:** baseline dose effects concentrate at the MCI stage
  (97 / 179 regions) and are far fewer at the AD stage (13 / 179).
- **Amyloid dependence:** longitudinally, 115 / 179 regions are significant in
  amyloid-positive participants versus only 14 / 179 in amyloid-negative;
  amyloid-positive rates rise with dose (31 % / 70 % / 90 %), matching Table 1.
- Per-region effect sizes correlate with the previously reported values at
  **r = 0.95** across all 179 regions.

**Head-size term.** The size covariate is **total brain volume = the sum of the
179 regional volumes**. (The harmonized files carry a `brain_volume_mm3` column
that is a broken constant, so it is not used; because ComBat is applied upstream
and volumes are residualized/standardized here, no separate eTIV is required to
reproduce the published SD/year effects.)

## Repository layout

```
analysis/
  build_panel.py          # load the 6 harmonized genotype files -> residualized
                          #   SD-unit panel (time from scan date; within-visit
                          #   averaging; CN-referenced residualization)
  baseline_analyses.py    # cross-sectional ANOVA + BH-FDR, by stage (MCI, AD)
  longitudinal_analyses.py# two-stage per-subject slopes; dose-trend ANOVA +
                          #   BH-FDR; per-region excess-atrophy table
figures/
  make_brain_map.py       # Destrieux cortical surface map (nilearn; auto-atlas)
  make_brain_map.R        # R/ggseg equivalent
  plot_top_regions_bar.py # top-N regions by excess atrophy
  plot_subcortical.py     # subcortical structures panel
  plot_forest.py          # amyloid-stratified forest plot
results/                  # region-level summary tables (publishable; Suppl. S1–S6)
data/                     # <- place ADNI-derived harmonized inputs here (gitignored)
requirements.txt
LICENSE
```

## Inputs (not included — ADNI-derived)

Place the six upstream harmonized files in `data/`, named exactly:

```
data/Genotype_2_2_harmonized_volumes_residuals.csv
data/Genotype_2_3_harmonized_volumes_residuals.csv
data/Genotype_3_3_harmonized_volumes_residuals.csv
data/Genotype_2_4_harmonized_volumes_residuals.csv
data/Genotype_3_4_harmonized_volumes_residuals.csv
data/Genotype_4_4_harmonized_volumes_residuals.csv
```

Each is a per-scan table with the 179 DL+DiReCt region columns
(`Left-Cerebral-White-Matter` … `rh-S_temporal_transverse`) plus metadata
(`subject_id`, `image_date`, `Diagnosis_Label`, `Age`, `Sex_Code`, …). These are
**derived from ADNI and must not be committed or redistributed** under the ADNI
Data Use Agreement; `.gitignore` blocks them. Qualified investigators regenerate
them from ADNI raw data with the upstream pipeline.

**Amyloid status** comes from ADNI's Berkeley amyloid-PET summary table
(per-participant amyloid-positive flag and Centiloids), merged to each subject by
the amyloid scan nearest their MRI baseline. To rebuild from scratch, put that
table in `data/` and set `AMYLOID_FILE` in `build_panel.py`; `baseline_analyses.py`
and `longitudinal_analyses.py` then emit the amyloid-adjusted and amyloid-stratified
results (included in `results/`) automatically.

## How to run

```bash
pip install -r requirements.txt

python analysis/build_panel.py           # data/*.csv -> panel_residualized.csv
python analysis/longitudinal_analyses.py # -> results/RESULT_longitudinal_ALL.csv, BRAINMAP_effect_sizes.csv
python analysis/baseline_analyses.py     # -> results/RESULT_baseline_{MCI,AD}.csv

python figures/plot_top_regions_bar.py   # reads results/, writes PNG+PDF
python figures/plot_subcortical.py
python figures/plot_forest.py            # amyloid-stratified forest
python figures/make_brain_map.py         # needs internet once to fetch the atlas
```

The `figures/plot_*.py` scripts run from the region-level tables in `results/`,
so figures regenerate without ADNI access.

## Data availability

Raw MRI and clinical data are available from ADNI (<https://adni.loni.usc.edu>)
to qualified investigators under its Data Use Agreement. Region-level summary
statistics are provided in `results/`. Subject-level derived regional volumes are
available to ADNI-approved investigators on request, consistent with ADNI policy.

## Notes

- Portions of this analysis code were drafted with the assistance of generative
  AI tools and subsequently reviewed, verified and edited by the authors,
  consistent with the manuscript's declaration of generative-AI use.
- Released under the MIT License (see `LICENSE`); update the copyright holders if
  the author list changes.
- Citation for the associated paper to be added on acceptance.
