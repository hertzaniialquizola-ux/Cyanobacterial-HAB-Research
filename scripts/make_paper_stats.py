"""make_paper_stats.py -- recomputes outputs/paper_stats.json, _county_summary.csv and _hist_data.csv from
outputs/county_year_panel_v2.csv for the manuscript tables/figures. The analytic sample here is the NO_LST sample
(every lake-matched county-year with exposure, covariates and rurality; the MODIS lake-temperature covariate is dropped).
Run from the project root:  python3 scripts/make_paper_stats.py"""
import json, shutil
from pathlib import Path
import numpy as np, pandas as pd
R = Path(__file__).resolve().parent.parent
p = pd.read_csv(R / "outputs/county_year_panel_v2.csv", dtype={"county_fips": str})
p = p[~p.county_fips.str[:2].isin(["02","15","72","78","66","60","69"])]   # CONUS + DC only, as in the models
old = R / "outputs/paper_stats.json"
if old.exists(): shutil.copy(old, R / "outputs/paper_stats_prev.json")
p["deaths_suppressed"] = p["deaths_suppressed"].astype(bool)
cov = ["well_water_fraction", "surface_supply_fraction", "poverty_rate", "median_household_income", "pct_65plus", "pct_black", "pct_hispanic", "pct_aian", "rucc"]
base = p["population"].gt(0) & p[cov].notna().all(axis=1)
an = p[p["cei"].notna() & p["px_share_20k"].notna() & base].copy()
ex = p[p["cei"].notna()]
S = dict(panel_rows=int(len(p)), panel_counties=int(p.county_fips.nunique()), panel_supp_pct=float(p.deaths_suppressed.mean() * 100),
         n_cei_rows=int(len(ex)), n_cei_counties=int(ex.county_fips.nunique()), n_lst=int(len(an)),
         allcounty_rows=int(base.sum()), allcounty_counties=int(p[base].county_fips.nunique()),
         analytic_rows=int(len(an)), analytic_counties=int(an.county_fips.nunique()), analytic_supp_pct=float(an.deaths_suppressed.mean() * 100),
         analytic_zero_rows=int(((an.deaths == 0) & ~an.deaths_suppressed).sum()), analytic_zero_pct=float(((an.deaths == 0) & ~an.deaths_suppressed).mean() * 100),
         pop_filled=int(p["pop_was_missing"].sum()), new_cei_ge99_pct=float((ex.cei >= 0.99).mean() * 100), old_cei_ge99_pct=float((ex.cei_old_buggy >= 0.99).mean() * 100),
         analytic_years=sorted(an.year.unique().tolist()), px_iqr=float(an.px_share_20k.quantile(.75) - an.px_share_20k.quantile(.25)), px_median=float(an.px_share_20k.median()))
an["tert"] = pd.qcut(an.px_share_20k, 3, labels=["Low", "Middle", "High"])
def summ(d):
    return dict(rows=int(len(d)), counties=int(d.county_fips.nunique()), supp_pct=float(d.deaths_suppressed.mean() * 100), pop_med=float(d.population.median()),
                age65=float(d.pct_65plus.mean() * 100), black=float(d.pct_black.mean() * 100), hisp=float(d.pct_hispanic.mean() * 100), pov=float(d.poverty_rate.mean() * 100),
                inc=float(d.median_household_income.median()), ww=float(d.well_water_fraction.mean() * 100), sw=float(d.surface_supply_fraction.mean() * 100),
                expo_med=float(d.px_share_20k.median()), expo_lo=float(d.px_share_20k.min()), expo_hi=float(d.px_share_20k.max()),
                expo_q1=float(d.px_share_20k.quantile(.25)), expo_q3=float(d.px_share_20k.quantile(.75)),
                deaths_med_nonsupp=float(d.loc[~d.deaths_suppressed, "deaths"].median()), rucc_med=float(d.rucc.median()),
                age65_sd=float(d.pct_65plus.std() * 100), pov_sd=float(d.poverty_rate.std() * 100),
                inc_iqr=float(d.median_household_income.quantile(.75) - d.median_household_income.quantile(.25)),
                ww_sd=float(d.well_water_fraction.std() * 100), sw_sd=float(d.surface_supply_fraction.std() * 100))
T = {"All": summ(an)}
for t in ["Low", "Middle", "High"]: T[t] = summ(an[an.tert == t])
S["table1"] = T
S["pct_counties_multi_tertile"] = float((an.groupby("county_fips")["tert"].nunique() > 1).mean() * 100)
lake_rows = p["cei"].notna()
S["county_supp_lake"] = float(p[lake_rows].deaths_suppressed.mean() * 100); S["county_supp_nolake"] = float(p[~lake_rows].deaths_suppressed.mean() * 100)
S["pop_med_lake"] = float(p[lake_rows].population.median()); S["pop_med_nolake"] = float(p[~lake_rows].population.median())
S["n_lake_counties_any"] = S["n_cei_counties"]; S["n_nolake"] = S["panel_counties"] - S["n_cei_counties"]
json.dump(S, open(R / "outputs/paper_stats.json", "w"), indent=1)
cs = p.groupby("county_fips").agg(supp=("deaths_suppressed", "mean"), lake=("cei", lambda x: bool(x.notna().any())), px=("px_share_20k", "mean"), cei=("cei", "mean")).reset_index()
cs.to_csv(R / "outputs/_county_summary.csv", index=False)
p[["county_fips", "year", "cei", "cei_old_buggy", "px_share_20k"]].to_csv(R / "outputs/_hist_data.csv", index=False)
print({k: v for k, v in S.items() if k != "table1"})
for t, d in T.items(): print(t, d["rows"], d["counties"], round(d["expo_med"], 3), round(d["supp_pct"], 1), round(d["sw"], 1))
