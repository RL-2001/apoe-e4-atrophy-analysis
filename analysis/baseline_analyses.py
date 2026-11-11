#!/usr/bin/env python3
"""
Cross-sectional (baseline) analysis of APOE ε4 dose and regional volume.

Reads panel_residualized.csv (SD-unit regional residuals, from build_panel.py),
restricts to each subject's baseline visit, and tests whether baseline regional
volume differs by ε4 dose (0/1/2) within each disease stage (MCI, AD), with a
one-way ANOVA and Benjamini-Hochberg FDR across the 179 regions. If an
`amyloid_pos` column is present, it also reports the amyloid-adjusted dose
effect and the amyloid-positive-only stratum.

Outputs (written to ../results/):
    RESULT_baseline_MCI.csv, RESULT_baseline_AD.csv
    RESULT_baseline_MCI_amyloidPOS.csv, RESULT_baseline_AD_amyloidPOS.csv (if amyloid)
"""
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
META = {"subj", "time_yr", "e4_dose", "Age", "Sex_Code", "Diagnosis_Label",
        "amyloid_pos", "centiloids"}

d = pd.read_csv(ROOT / "panel_residualized.csv")
REGIONS = [c for c in d.columns if c not in META]
HAS_AMY = "amyloid_pos" in d.columns and d["amyloid_pos"].notna().any()
MCI = {"MCI", "EMCI", "LMCI"}


def block(sub, tag):
    sub = sub[sub["time_yr"] < 0.01].copy()          # baseline visit
    rows = []
    for r in REGIONS:
        g = [sub[sub.e4_dose == k][r].dropna() for k in (0, 1, 2)]
        if min(len(g[0]), len(g[2])) < 5:
            rows.append({"Region": r}); continue
        F, p = stats.f_oneway(*g)
        rec = {"Region": r, "F": F, "p_unadj": p,
               "mean0": g[0].mean(), "mean2": g[2].mean(),
               "e4e4_worse": g[2].mean() < g[0].mean()}
        if HAS_AMY:
            xa = sub[[r, "e4_dose", "amyloid_pos"]].dropna()
            if xa["amyloid_pos"].nunique() > 1 and len(xa) > 30:
                m = smf.ols(f'Q("{r}") ~ e4_dose + C(amyloid_pos)', data=xa).fit()
                rec["p_amyAdj"] = m.pvalues.get("e4_dose", np.nan)
        rows.append(rec)
    res = pd.DataFrame(rows)
    for col, name in [("p_unadj", "FDR_unadj"), ("p_amyAdj", "FDR_amyloidAdj")]:
        if col in res:
            ok = res[col].notna()
            res[name] = np.nan
            if ok.sum() > 1:
                res.loc[ok, name] = multipletests(res.loc[ok, col], method="fdr_bh")[1]
    res.to_csv(ROOT / "results" / f"RESULT_baseline_{tag}.csv", index=False)
    nsig = int((res.get("FDR_unadj", pd.Series(dtype=float)) < 0.05).sum())
    print(f"BASELINE {tag}: n={sub.subj.nunique()} | FDR-sig(unadj)={nsig}/{len(REGIONS)}")
    return res


mci = d[d["Diagnosis_Label"].isin(MCI)]
ad = d[d["Diagnosis_Label"] == "AD"]
block(mci, "MCI")
block(ad, "AD")
if HAS_AMY:
    block(mci[mci.amyloid_pos == 1], "MCI_amyloidPOS")
    block(ad[ad.amyloid_pos == 1], "AD_amyloidPOS")
