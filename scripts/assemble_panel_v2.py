"""
assemble_panel_v2.py -- builds outputs/county_year_panel_v2.csv: the old panel
(mortality, population, well water, poverty, income, LST) with
  * the CORRECTED CyAN exposure (outputs/alt_cei_county_year.csv) replacing the buggy CEI
    (old value kept as cei_old_buggy)
  * ACS age/sex/race-ethnicity, USGS surface-supply fraction, optional USDA RUCC.
Prints a sample-attrition table (use it for the flow diagram).

Run order:  build_alt_cei.py, get_census_demographics.py, get_usgs_surface_supply.py,
            [optional: put USDA RUCC csv at data/raw/usda/rucc_2023.csv], then this script.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent.parent
D = {"county_fips": str}
old = pd.read_csv(ROOT / "outputs/county_year_panel_full.csv", dtype=D)
old = old.rename(columns={"cei": "cei_old_buggy", "n_periods_observed": "n_periods_old"})
alt = pd.read_csv(ROOT / "outputs/alt_cei_county_year.csv", dtype=D)
p = old.merge(alt, on=["county_fips", "year"], how="left")

dem = json.load(open(ROOT / "data/raw/census/demographics_by_county_2022.json"))
dm = pd.DataFrame(dem[1:], columns=dem[0])
dm["county_fips"] = dm["state"] + dm["county"]
num = lambda c: pd.to_numeric(dm[c], errors="coerce")
tot = num("B01001_001E")
old65 = sum(num(f"B01001_0{n}E") for n in list(range(20, 26)) + list(range(44, 50)))
dm["pct_65plus"] = old65 / tot
dm["pct_male"] = num("B01001_002E") / tot
dm["pct_black"] = num("B03002_004E") / num("B03002_001E")
dm["pct_hispanic"] = num("B03002_012E") / num("B03002_001E")
dm["pct_aian"] = num("B03002_005E") / num("B03002_001E")      # needs re-run of get_census_demographics.py
dm["acs_total_pop"] = tot
p = p.merge(dm[["county_fips", "pct_65plus", "pct_male", "pct_black", "pct_hispanic", "pct_aian", "acs_total_pop"]], on="county_fips", how="left")

# ---- POPULATION FIX (2026-10-03) -------------------------------------------------
# CDC WONDER was exported with "Show Zero Values: False": county-years with ZERO deaths have no
# row, hence no population. They were silently dropped by the model's population > 0 filter =
# selection on the outcome. Fill their denominator from the same county's other CDC years
# (median), falling back to the ACS 2022 total population.
p["population"] = pd.to_numeric(p["population"], errors="coerce")
p["pop_was_missing"] = p["population"].isna()
cty_med = p.groupby("county_fips")["population"].transform("median")
p.loc[p["pop_was_missing"], "population"] = cty_med[p["pop_was_missing"]]
still = p["population"].isna()
p.loc[still, "population"] = p.loc[still, "acs_total_pop"]
print(f"Population filled for {int(p['pop_was_missing'].sum())} county-years "
      f"({int(still.sum())} needed the ACS fallback); of the filled rows, share with deaths == 0: "
      f"{(p.loc[p['pop_was_missing'], 'deaths'] == 0).mean():.3f}; still missing: {int(p['population'].isna().sum())}")

ss = pd.read_csv(ROOT / "data/raw/usgs/surface_supply_by_county.csv", dtype=D)
p = p.merge(ss, on="county_fips", how="left")

rucc_path = ROOT / "data/raw/usda/rucc_2023.csv"
if rucc_path.exists():
    r = pd.read_csv(rucc_path, dtype=str, encoding="latin1")
    if "Attribute" in r.columns:            # USDA long format
        r = r[r["Attribute"] == "RUCC_2023"].rename(columns={"Value": "rucc"})
    else:
        r = r.rename(columns={"RUCC_2023": "rucc"})
    r["county_fips"] = r["FIPS"].str.zfill(5); r["rucc"] = pd.to_numeric(r["rucc"], errors="coerce")
    p = p.merge(r[["county_fips", "rucc"]], on="county_fips", how="left")
    print("RUCC merged.")
else:
    print("No RUCC file at data/raw/usda/rucc_2023.csv -- rurality not included (download from USDA ERS if you want it).")

# ---- optional CDC PLACES (diabetes, obesity, binge drinking, uninsured) ----------------
places_dir = ROOT / "data/raw/places"
pf = sorted(places_dir.glob("*.csv")) if places_dir.exists() else []
if pf:
    pl = pd.read_csv(pf[0], dtype=str, low_memory=False)
    pl = pl[pl["MeasureId"].isin(["DIABETES", "OBESITY", "BINGE", "ACCESS2"])]
    if "DataValueTypeID" in pl.columns and (pl["DataValueTypeID"] == "AgeAdjPrv").any():
        pl = pl[pl["DataValueTypeID"] == "AgeAdjPrv"]
    pl["county_fips"] = pl["LocationID"].str.zfill(5)
    pl["v"] = pd.to_numeric(pl["Data_Value"], errors="coerce")
    wide = pl.pivot_table(index="county_fips", columns="MeasureId", values="v", aggfunc="first").reset_index()
    wide = wide.rename(columns={"DIABETES": "diabetes", "OBESITY": "obesity", "BINGE": "binge", "ACCESS2": "uninsured"})
    p = p.merge(wide, on="county_fips", how="left")
    print("PLACES merged from", pf[0].name, "->", [c for c in ["diabetes", "obesity", "binge", "uninsured"] if c in p.columns])
else:
    print("No PLACES csv in data/raw/places/ -- diabetes/obesity/alcohol/insurance covariates not included.")

# ---- lagged / lead / within-between versions of the primary exposure -----------------
p = p.sort_values(["county_fips", "year"]).reset_index(drop=True)
g = p.groupby("county_fips")["px_share_20k"]
p["px_lag1"] = g.shift(1)
p["px_lead1"] = g.shift(-1)
p["px_cum3"] = pd.concat([g.shift(1), g.shift(2), g.shift(3)], axis=1).mean(axis=1, skipna=True).where(
    pd.concat([g.shift(1), g.shift(2), g.shift(3)], axis=1).notna().sum(axis=1) >= 2)
p["px_cty_mean"] = g.transform("mean")
p["px_within"] = p["px_share_20k"] - p["px_cty_mean"]

p.to_csv(ROOT / "outputs/county_year_panel_v2.csv", index=False)
print("\n=== SAMPLE ATTRITION (county-year rows) ===")
steps = [("all county-year rows", np.ones(len(p), bool)),
         ("has corrected CyAN exposure (cei)", p["cei"].notna()),
         ("+ has LST", p["cei"].notna() & p["lst_day_c"].notna())]
cov = ["well_water_fraction", "surface_supply_fraction", "poverty_rate", "median_household_income", "pct_65plus", "pct_black", "pct_hispanic", "pct_aian", "population"]
steps.append(("+ has all other covariates", steps[-1][1] & p[cov].notna().all(axis=1)))
for name, m in steps:
    print(f"{name:42s} {int(m.sum()):7d} rows, {p.loc[m,'county_fips'].nunique():5d} counties, {p.loc[m,'deaths_suppressed'].mean():.1%} suppressed")
print("\nExposure distribution (corrected):")
print(p[["cei", "px_share_20k", "px_share_100k", "wk_frac_5pct_20k", "log10_cells_mean"]].describe().round(3).to_string())
print("\nOLD (buggy) vs corrected CEI correlation:", round(p[["cei_old_buggy", "cei"]].corr().iloc[0, 1], 3))
print("share of county-years with corrected cei >= 0.99:", round((p["cei"] >= 0.99).mean() / p["cei"].notna().mean(), 3))
