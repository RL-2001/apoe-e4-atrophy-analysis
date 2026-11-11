#!/usr/bin/env python3
"""
Build the residualized analysis panel from the upstream ComBat-harmonized
volume files (one file per APOE genotype) produced by the preprocessing/
inference pipeline:  https://github.com/Mappo23/mri-preprocess-inference-pipeline

Each input is a per-scan table with the 179 DL+DiReCt regional volumes
(columns 'Left-Cerebral-White-Matter' ... 'rh-S_temporal_transverse') plus
metadata (subject_id, image_date, Diagnosis_Label, Age, Sex_Code, Batch, ...).
ComBat harmonization across scanner manufacturer has already been applied
upstream, so this script only:

  1. tags each file with its ε4 allele dose (0/1/2),
  2. derives years-from-baseline from image_date (robust to ADNI's many visit
     codes; do NOT rely on visit labels),
  3. averages within-visit replicate scans,
  4. forms the head-size term as total brain volume = sum of the 179 regions
     (NB: the files' `brain_volume_mm3` column is a broken constant, so it is
     not used),
  5. residualizes every region in the cognitively-normal (CN) baseline sample
     against total brain volume, age and sex, then z-scores to SD units,
  6. (optional) merges amyloid-PET status if an amyloid table is supplied.

Output: panel_residualized.csv  (one row per subject-visit; region columns in
SD units — this is what the paper reports as "SD/year" effects).

Place the harmonized inputs in ./data/ (they are ADNI-derived and must NOT be
committed; see README and .gitignore).
"""
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

DATA = Path(__file__).resolve().parent.parent / "data"

# genotype file -> ε4 allele dose
FILES = {
    "Genotype_2_2_harmonized_volumes_residuals.csv": 0,
    "Genotype_2_3_harmonized_volumes_residuals.csv": 0,
    "Genotype_3_3_harmonized_volumes_residuals.csv": 0,   # = harmonized_volumes_residuals_33.csv
    "Genotype_2_4_harmonized_volumes_residuals.csv": 1,
    "Genotype_3_4_harmonized_volumes_residuals.csv": 1,   # = harmonized_volumes_residuals_34.csv
    "Genotype_4_4_harmonized_volumes_residuals.csv": 2,
}
# optional amyloid-PET table (Berkeley), with columns PTID, SCANDATE,
# AMYLOID_STATUS (0/1), CENTILOIDS. Set to None to skip amyloid merge.
AMYLOID_FILE = None   # e.g. "All_Subjects_UCBERKELEY_AMY_6MM.csv"


def load_one(path, dose):
    df = pd.read_csv(path, low_memory=False)
    cols = list(df.columns)
    regions = cols[cols.index("Left-Cerebral-White-Matter"):cols.index("brain_volume_mm3")]
    df["e4_dose"] = dose
    df["subj"] = df["subject_id"].astype(str)
    df["date"] = pd.to_datetime(df["image_date"], errors="coerce")
    return df, regions


def main():
    frames, regions = [], None
    for fname, dose in FILES.items():
        p = DATA / fname
        if not p.exists():
            raise FileNotFoundError(f"missing harmonized input: {p}")
        df, R = load_one(p, dose)
        regions = regions or R
        assert R == regions, f"region columns differ in {fname}"
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    # 2) years from baseline (per subject earliest scan)
    df["time_yr"] = (df["date"] - df.groupby("subj")["date"].transform("min")).dt.days / 365.25

    # 3) average within-visit replicate scans
    agg = {r: "mean" for r in regions}
    agg.update({"e4_dose": "first", "Age": "first", "Sex_Code": "first",
                "Diagnosis_Label": "first", "time_yr": "first"})
    panel = df.groupby(["subj", "date"]).agg(agg).reset_index()

    # 4) head-size term = total brain volume (sum of regions)
    panel["TBV"] = panel[regions].sum(axis=1)

    # 5) CN-baseline residualization -> z (SD units)
    cn = panel[(panel["Diagnosis_Label"] == "CN") & (panel["time_yr"] < 0.01)].copy()
    Z = pd.DataFrame(index=panel.index)
    for r in regions:
        m = smf.ols(f'Q("{r}") ~ TBV + Age + C(Sex_Code)',
                    data=cn[[r, "TBV", "Age", "Sex_Code"]].dropna()).fit()
        e = panel[r] - m.predict(panel[["TBV", "Age", "Sex_Code"]])
        mu, sd = e.loc[cn.index].mean(), e.loc[cn.index].std()
        Z[r] = (e - mu) / sd

    out = pd.concat([panel[["subj", "time_yr", "e4_dose", "Age", "Sex_Code",
                            "Diagnosis_Label"]], Z], axis=1)

    # 6) optional amyloid merge (one call per subject, scan nearest MRI baseline)
    if AMYLOID_FILE:
        amy = pd.read_csv(DATA / AMYLOID_FILE, low_memory=False)
        amy["subj"] = amy["PTID"].astype(str)
        amy["adate"] = pd.to_datetime(amy["SCANDATE"], errors="coerce")
        amy = amy[amy["AMYLOID_STATUS"].isin([0, 1])]
        base = panel.sort_values("time_yr").groupby("subj")["date"].first()
        amy = amy.merge(base.rename("mri_base").reset_index(), on="subj")
        amy["dd"] = (amy["adate"] - amy["mri_base"]).abs().dt.days
        best = amy.sort_values("dd").groupby("subj").first()
        out["amyloid_pos"] = out["subj"].map(best["AMYLOID_STATUS"])
        out["centiloids"] = out["subj"].map(best["CENTILOIDS"])
    else:
        out["amyloid_pos"] = np.nan

    out.to_csv(Path(__file__).resolve().parent.parent / "panel_residualized.csv", index=False)
    print(f"panel_residualized.csv: {out['subj'].nunique()} subjects, "
          f"{len(out)} subject-visits, {len(regions)} regions")
    print("by dose:", panel.groupby("e4_dose")["subj"].nunique().to_dict())


if __name__ == "__main__":
    main()
