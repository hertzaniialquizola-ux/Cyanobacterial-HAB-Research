"""Build the 7-year POOLED, age-standardized county dataset.

Inputs (CDC WONDER exports, tab-delimited .txt, saved in data/raw/cdc_pooled/ -- see the settings in
the instructions):
  county_kidney.txt     county-level, 2018-2024 pooled, UCD N18
  county_liver.txt      county-level, 2018-2024 pooled, UCD K74 + K76
  county_k70.txt        (optional) alcohol-related liver, UCD K70
  county_negctrl.txt    (optional) negative-control cause (unintentional injury V01-X59,Y85-Y86)
  national_age_kidney.txt / national_age_liver.txt / national_age_k70.txt / national_age_negctrl.txt
                        same causes, NOT by county: grouped by Ten-Year Age Groups only, 2018-2024
Also needs: data/raw/census/acs_age_groups_2022.csv (get_acs_age_groups.py), outputs/county_year_panel_v2.csv
Writes outputs/county_pooled_panel.csv.

Method: indirect standardization.  expected_deaths(county) = person_years x sum_g share_g(county) * rate_g(nation),
where share_g is the county's 2022 ACS age distribution in 9 broad groups and rate_g the national 2018-2024
age-specific death rate for the same cause.  The model then uses log(expected) as the offset, so the exposure
coefficient is a rate ratio adjusted for age structure (not just % aged 65+)."""
import csv, os, re, sys
import numpy as np, pandas as pd
D = "data/raw/cdc_pooled"
def find_file(stem):
    """accept .txt, .tsv or .csv (WONDER TSV export) for a given stem"""
    for ext in (".txt", ".tsv", ".csv"):
        if os.path.exists(stem + ext): return stem + ext
    return None
def read_wonder(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = f.read().splitlines()
    hdr = None; rows = []
    for ln in lines:
        if ln.startswith('"---"') or ln.startswith("---"): break
        cells = next(csv.reader([ln], delimiter="\t", quotechar='"'))
        if hdr is None:
            if "Deaths" in cells: hdr = cells
            continue
        if len(cells) == len(hdr): rows.append(cells)
    df = pd.DataFrame(rows, columns=hdr)
    return df
def num(s):
    s = str(s).strip().replace(",", "")
    try: return float(s)
    except ValueError: return np.nan            # "Suppressed", "Not Available", "Not Applicable", "Unreliable", ...
CAUSES = {"kidney": "kidney", "liver": "liver", "k70": "k70", "negctrl": "negctrl"}
age = pd.read_csv("data/raw/census/acs_age_groups_2022.csv", dtype={"county_fips": str}).set_index("county_fips")
GROUPS = list(age.columns)                                   # ag_0_14 ... ag_85p
lo_map = {0: "ag_0_14", 15: "ag_15_24", 25: "ag_25_34", 35: "ag_35_44", 45: "ag_45_54", 55: "ag_55_64",
          65: "ag_65_74", 75: "ag_75_84", 85: "ag_85p"}
def grp_of(code, label):
    for s in (code, label):
        m = re.search(r"\d+", str(s))
        if m:
            lo = int(m.group())
            if lo < 15: return "ag_0_14"
            return lo_map.get(lo)
    return None
share = age.div(age.sum(axis=1), axis=0)
out = None
for key in CAUSES:
    cpath, npath = find_file(f"{D}/county_{key}"), find_file(f"{D}/national_age_{key}")
    if not (cpath and npath):
        print(f"[skip] {key}: need county_{key} and national_age_{key} (.txt/.tsv) in {D}"); continue
    nat = read_wonder(npath)
    cod = [c for c in nat.columns if "Age Groups" in c and "Code" in c]
    lab = [c for c in nat.columns if "Age Groups" in c and "Code" not in c]
    nat["g"] = [grp_of(nat[cod[0]][i] if cod else "", nat[lab[0]][i] if lab else "") for i in range(len(nat))]
    nat["Deaths"] = nat["Deaths"].map(num); nat["Population"] = nat["Population"].map(num)
    ng = nat.dropna(subset=["g"]).groupby("g")[["Deaths", "Population"]].sum()
    rate = (ng["Deaths"] / ng["Population"]).reindex(GROUPS)
    print(f"{key}: national age-specific rates per 100k:\n{(rate*1e5).round(2).to_string()}")
    if rate.isna().any(): sys.exit(f"{key}: missing national rate for some age group -- check the age-group export")
    cty = read_wonder(cpath)
    cc = [c for c in cty.columns if c == "County Code"][0]
    c = pd.DataFrame({"county_fips": cty[cc].str.zfill(5), "deaths_raw": cty["Deaths"].astype(str),
                      "py": cty["Population"].map(num)})
    c[f"{key}_suppressed"] = c.deaths_raw.str.contains("Suppressed")
    c[f"{key}_deaths"] = c.deaths_raw.map(num)
    c = c[c.county_fips.str.fullmatch(r"\d{5}")]            # drops WONDER "Total" rows
    c = c.drop(columns="deaths_raw").set_index("county_fips")
    c = c[~c.index.duplicated()]
    n_fill = int(c.py.isna().sum())
    c["py"] = c["py"].fillna(age.sum(axis=1).reindex(c.index) * 7)   # e.g. Connecticut counties: WONDER has no population after the 2022 planning-region change
    if n_fill: print(f"  {key}: person-years filled from ACS total x 7 for {n_fill} counties (CT)")
    exp_rate = (share.reindex(c.index) * rate.values).sum(axis=1)       # expected deaths per person-year
    c[f"{key}_expected"] = c["py"] * exp_rate
    c = c.rename(columns={"py": f"py_{key}"})
    out = c if out is None else out.join(c.drop(columns=[]), how="outer", rsuffix="_dup")
if out is None: sys.exit("No WONDER files found.")
out["person_years"] = out[[c for c in out.columns if c.startswith("py_")]].bfill(axis=1).iloc[:, 0]
# sanity: pooled population should be ~7x the annual population
acs_tot = age.sum(axis=1).reindex(out.index)
print("median person_years / ACS total (expect ~7):", round(float((out.person_years / acs_tot).median()), 2))
# exposure + covariate means from the annual panel
p = pd.read_csv("outputs/county_year_panel_v2.csv", dtype={"county_fips": str})
cols = {"px_share_20k": "px_mean", "cei": "cei_mean", "wk_frac_5pct_20k": "wk5_mean", "log10_cells_mean": "logc_mean",
        "well_water_fraction": "well", "surface_supply_fraction": "surf", "poverty_rate": "pov",
        "median_household_income": "inc", "pct_black": "blk", "pct_hispanic": "hsp", "pct_aian": "aian",
        "rucc": "rucc", "lst_day_c": "lst", "diabetes": "diabetes", "obesity": "obesity", "binge": "binge",
        "uninsured": "uninsured"}
g = p.groupby("county_fips")[list(cols)].mean().rename(columns=cols)
g["n_years_exposure"] = p.groupby("county_fips").px_share_20k.count()
out = out.join(g, how="left")
os.makedirs("outputs", exist_ok=True)
out.reset_index().rename(columns={"index": "county_fips"}).to_csv("outputs/county_pooled_panel.csv", index=False)
print(f"{len(out)} counties; with exposure: {int(out.px_mean.notna().sum())} -> outputs/county_pooled_panel.csv")
for key in CAUSES:
    if f"{key}_suppressed" in out: print(f"  {key}: {100*out[f'{key}_suppressed'].mean():.1f}% of counties still suppressed (pooled 7 yrs; annual was ~50%)")
