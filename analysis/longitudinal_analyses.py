#!/usr/bin/env python3
"""
Longitudinal (two-stage) analysis of APOE ε4 dose and atrophy rate.

Reads panel_residualized.csv (SD-unit regional residuals, from build_panel.py),
estimates a per-subject linear slope for every region (subjects with >=3
timepoints), then tests whether the slope differs by ε4 dose (0/1/2) with a
one-way ANOVA and Benjamini-Hochberg FDR across the 179 regions. If an
`amyloid_pos` column is present, it also reports the amyloid-adjusted dose
effect and the amyloid-positive-only stratum.

Effects are in SD/year; "excess atrophy" = (dose-0 slope) - (dose-2 slope),
so positive = ε4/ε4 atrophies faster.

Outputs (written to ../results/):
    RESULT_longitudinal_ALL.csv
    BRAINMAP_effect_sizes.csv            (Region, per-dose slopes, excess, FDR)
    RESULT_longitudinal_amyloidPOS.csv   (if amyloid present)
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


def subject_slopes(data):
    """per-subject OLS slope for each region; subjects with >=3 timepoints."""
    nt = data.groupby("subj")["time_yr"].nunique()
    data = data[data["subj"].isin(nt[nt >= 3].index)]

    def _fit(g):
        t = g["time_yr"].values
        if len(t) < 3 or np.var(t) == 0:
            return None
        return pd.Series({r: np.polyfit(t, g[r].values, 1)[0] for r in REGIONS})

    S = data.groupby("subj").apply(_fit)
    S["e4_dose"] = data.groupby("subj")["e4_dose"].first()
    if HAS_AMY:
        S["amyloid_pos"] = data.groupby("subj")["amyloid_pos"].first()
    return S


def analyse(S, tag):
    rows = []
    for r in REGIONS:
        a = S[S.e4_dose == 0][r].dropna()
        b = S[S.e4_dose == 1][r].dropna()
        c = S[S.e4_dose == 2][r].dropna()
        if min(len(a), len(c)) < 5:
            rows.append({"Region": r}); continue
        F, p = stats.f_oneway(a, b, c)
        rec = {"Region": r, "Slope_0alleles": a.mean(), "Slope_1allele": b.mean(),
               "Slope_2alleles": c.mean(), "Excess_atrophy_SDyr": a.mean() - c.mean(),
               "faster_e4e4": c.mean() < a.mean(), "p_unadj": p}
        if HAS_AMY:
            xa = S[[r, "e4_dose", "amyloid_pos"]].dropna()
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
    res.to_csv(ROOT / "results" / f"RESULT_longitudinal_{tag}.csv", index=False)
    nsig = int((res.get("FDR_unadj", pd.Series(dtype=float)) < 0.05).sum())
    print(f"LONGITUDINAL {tag}: subjects={len(S)} | FDR-sig(unadj)={nsig}/{len(REGIONS)}")
    return res


S = subject_slopes(d)
res = analyse(S, "ALL")

# brain-map effect table (used by figures/)
cols = ["Region", "Slope_0alleles", "Slope_1allele", "Slope_2alleles",
        "Excess_atrophy_SDyr", "FDR_unadj"] + (["FDR_amyloidAdj"] if HAS_AMY else [])
res[cols].to_csv(ROOT / "results" / "BRAINMAP_effect_sizes.csv", index=False)

if HAS_AMY:
    analyse(S[S.amyloid_pos == 1], "amyloidPOS")
    analyse(S[S.amyloid_pos == 0], "amyloidNEG")

# headline
top = res.dropna(subset=["Excess_atrophy_SDyr"]).sort_values("Excess_atrophy_SDyr", ascending=False)
print("\nTop 8 regions by ε4/ε4 excess atrophy (SD/yr):")
for _, x in top.head(8).iterrows():
    print(f"  {x['Region']:32s} {x['Excess_atrophy_SDyr']:+.3f}  FDR={x['FDR_unadj']:.1e}")
