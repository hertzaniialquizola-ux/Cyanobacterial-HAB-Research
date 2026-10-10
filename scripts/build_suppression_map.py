import sys
import sqlite3
import glob
import os
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import matplotlib.patches as mpatches

sys.path.insert(0, "/home/claude")
from shp_dbf import read_dbf, read_shp_polygons, point_in_rings, ring_bbox

BASE = "/mnt/user-data/uploads/Cyanobacterial-HAB-Research/data/raw"
SHP = f"{BASE}/census/county_boundaries/cb_2022_us_county_20m.shp"
DBF = f"{BASE}/census/county_boundaries/cb_2022_us_county_20m.dbf"
WONDER_CSV = f"{BASE}/cdc-wonder/Multiple Cause of Death, 2018-2024, Single Race.csv"
NHD_DIR = f"{BASE}/nhdplus/by_state"

YEARS = [2018, 2019, 2020, 2021, 2022, 2023, 2024]

STATE_FIPS_TO_ABBR = {
    "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA", "08": "CO", "09": "CT",
    "10": "DE", "11": "DC", "12": "FL", "13": "GA", "15": "HI", "16": "ID", "17": "IL",
    "18": "IN", "19": "IA", "20": "KS", "21": "KY", "22": "LA", "23": "ME", "24": "MD",
    "25": "MA", "26": "MI", "27": "MN", "28": "MS", "29": "MO", "30": "MT", "31": "NE",
    "32": "NV", "33": "NH", "34": "NJ", "35": "NM", "36": "NY", "37": "NC", "38": "ND",
    "39": "OH", "40": "OK", "41": "OR", "42": "PA", "44": "RI", "45": "SC", "46": "SD",
    "47": "TN", "48": "TX", "49": "UT", "50": "VT", "51": "VA", "53": "WA", "54": "WV",
    "55": "WI", "56": "WY", "72": "PR", "60": "AS", "66": "GU", "69": "MP", "78": "VI",
}
# CONUS + DC only, matching CyAN/NHDPlus scope in the project log
CONUS_ABBR = set(STATE_FIPS_TO_ABBR.values()) - {"AK", "HI", "PR", "AS", "GU", "MP", "VI"}

print("=== 1. Loading county boundaries ===", flush=True)
dbf_fields, dbf_rows = read_dbf(DBF)
shp_features = read_shp_polygons(SHP)
assert len(dbf_rows) == len(shp_features), (len(dbf_rows), len(shp_features))

counties = {}  # GEOID -> dict(name, statefp, rings, bbox)
for row, rings in zip(dbf_rows, shp_features):
    geoid = row["GEOID"]
    counties[geoid] = {
        "name": row["NAME"],
        "state_name": row.get("STATE_NAME", ""),
        "statefp": row["STATEFP"],
        "rings": rings,
        "bbox": ring_bbox(rings) if rings else None,
    }
print(f"Loaded {len(counties)} counties", flush=True)

print("=== 2. Parsing CDC WONDER suppression data ===", flush=True)
wonder = pd.read_csv(WONDER_CSV, dtype=str)
wonder = wonder.dropna(subset=["Year", "County Code"])
wonder["Year"] = wonder["Year"].str.strip()
wonder = wonder[wonder["Year"].str.isdigit()]
wonder["Year"] = wonder["Year"].astype(int)
wonder["County Code"] = wonder["County Code"].str.strip().str.zfill(5)
wonder = wonder[wonder["Year"].isin(YEARS)]

# status per (county, year): 'suppressed' | 'numeric'
status = {}
for _, r in wonder.iterrows():
    key = (r["County Code"], r["Year"])
    status[key] = "suppressed" if str(r["Deaths"]).strip() == "Suppressed" else "numeric"

suppression_summary = {}
for geoid in counties:
    n_supp = 0
    n_numeric = 0
    n_absent_zero = 0
    for y in YEARS:
        s = status.get((geoid, y))
        if s == "suppressed":
            n_supp += 1
        elif s == "numeric":
            n_numeric += 1
        else:
            n_absent_zero += 1  # row absent in WONDER export == true zero deaths
    suppression_summary[geoid] = {
        "n_suppressed": n_supp,
        "n_numeric": n_numeric,
        "n_absent_zero": n_absent_zero,
        "suppressed_frac": n_supp / len(YEARS),
        "ever_suppressed": n_supp > 0,
    }

total_cells = len(counties) * len(YEARS)
total_supp = sum(v["n_suppressed"] for v in suppression_summary.values())
print(f"County-year cells: {total_cells}, suppressed: {total_supp} "
      f"({100*total_supp/total_cells:.1f}%)", flush=True)

print("=== 3. Building per-state county spatial index (bbox pre-filter) ===", flush=True)
counties_by_state = {}
for geoid, c in counties.items():
    counties_by_state.setdefault(c["statefp"], []).append((geoid, c))

print("=== 4. Reading NHDPlus lake bounding-box centers per state, assigning to counties ===", flush=True)
# Two definitions of "lake-adjacent", since NHDPlus's raw waterbody layer
# includes every mapped pond/marsh/reservoir, most far too small for CyAN's
# 300m-pixel product to ever resolve a bloom in:
#  - lake_adjacent_any: any NHDPlus waterbody polygon at all (matches what the
#    project's join pipeline currently does with no size/type filter)
#  - lake_adjacent_big: LakePond/Reservoir features >= 1.0 sq km, roughly the
#    scale at which CyAN's CI_cyano product can actually resolve a bloom away
#    from shoreline mixed pixels
BIG_LAKE_MIN_AREASQKM = 1.0
BIG_LAKE_FTYPES = {"LakePond", "Reservoir"}

lake_adjacent_any = set()
lake_adjacent_big = set()
state_files = sorted(glob.glob(os.path.join(NHD_DIR, "waterbodies_*.gpkg")))
state_lake_counts = {}

abbr_to_fips = {v: k for k, v in STATE_FIPS_TO_ABBR.items()}

for path in state_files:
    abbr = os.path.basename(path).replace("waterbodies_", "").replace(".gpkg", "")
    fips = abbr_to_fips.get(abbr)
    if fips is None or fips not in counties_by_state:
        continue
    con = sqlite3.connect(path)
    cur = con.cursor()
    table = f"waterbodies_{abbr}"
    rtree_table = f"rtree_{table}_geom"
    cur.execute(
        f"SELECT t.ftype, t.areasqkm, r.minx, r.maxx, r.miny, r.maxy "
        f"FROM {table} t JOIN {rtree_table} r ON t.fid = r.id"
    )
    rows = cur.fetchall()
    con.close()

    candidates = counties_by_state[fips]
    matched_any = 0
    matched_big = 0
    for ftype, areasqkm, minx, maxx, miny, maxy in rows:
        cx = (minx + maxx) / 2.0
        cy = (miny + maxy) / 2.0
        is_big = (ftype in BIG_LAKE_FTYPES) and (areasqkm is not None) and (areasqkm >= BIG_LAKE_MIN_AREASQKM)
        for geoid, c in candidates:
            bb = c["bbox"]
            if bb is None:
                continue
            if not (bb[0] <= cx <= bb[2] and bb[1] <= cy <= bb[3]):
                continue
            if point_in_rings(cx, cy, c["rings"]):
                lake_adjacent_any.add(geoid)
                matched_any += 1
                if is_big:
                    lake_adjacent_big.add(geoid)
                    matched_big += 1
                break
    state_lake_counts[abbr] = (len(rows), matched_any, matched_big)
    print(f"  {abbr}: {len(rows)} waterbodies -> {matched_any} matched (any), {matched_big} matched (>=1 km2 lake/reservoir)", flush=True)

lake_adjacent = lake_adjacent_any  # backward-compat name used below for the "any" overlay

covered_abbrs = {os.path.basename(p).replace("waterbodies_", "").replace(".gpkg", "") for p in state_files}
missing_states = sorted(CONUS_ABBR - covered_abbrs)
print(f"States/DC covered: {len(covered_abbrs)} of {len(CONUS_ABBR)}. Missing: {missing_states}", flush=True)
print(f"Total counties flagged lake-adjacent: {len(lake_adjacent)}", flush=True)

print("=== 5. Computing suppression rate among lake-adjacent counties (national, non-SC) ===", flush=True)
conus_geoids = [g for g, c in counties.items() if STATE_FIPS_TO_ABBR.get(c["statefp"]) in CONUS_ABBR]

def group_rate(geoids):
    supp = sum(suppression_summary[g]["n_suppressed"] for g in geoids)
    total = len(geoids) * len(YEARS)
    return supp, total, (100 * supp / total if total else float("nan"))

national_supp, national_total, national_rate = group_rate(conus_geoids)

summary_txt = []
summary_txt.append(f"CONUS+DC counties: {len(conus_geoids)}")
summary_txt.append(f"National county-year suppression rate: {national_supp}/{national_total} = {national_rate:.1f}%")
summary_txt.append("")
summary_txt.append(f"-- Definition A: 'lake-adjacent' = any NHDPlus waterbody polygon (ponds/marshes/reservoirs included, no size filter) --")
lake_any = [g for g in conus_geoids if g in lake_adjacent_any]
nonlake_any = [g for g in conus_geoids if g not in lake_adjacent_any]
s, t, r = group_rate(lake_any)
summary_txt.append(f"  Lake-adjacent counties: {len(lake_any)} ({100*len(lake_any)/len(conus_geoids):.0f}% of CONUS). Suppression rate: {s}/{t} = {r:.1f}%")
s, t, r = group_rate(nonlake_any)
summary_txt.append(f"  Non-lake-adjacent counties: {len(nonlake_any)}. Suppression rate: {s}/{t} = {r:.1f}%")
summary_txt.append(f"  --> essentially no difference from the national rate: this definition is too inclusive to be useful (almost every county has *some* mapped waterbody).")
summary_txt.append("")
summary_txt.append(f"-- Definition B: 'lake-adjacent' = LakePond/Reservoir >= {BIG_LAKE_MIN_AREASQKM} sq km (roughly CyAN's 300m-pixel detection floor) --")
lake_big = [g for g in conus_geoids if g in lake_adjacent_big]
nonlake_big = [g for g in conus_geoids if g not in lake_adjacent_big]
s, t, r = group_rate(lake_big)
summary_txt.append(f"  Lake-adjacent counties: {len(lake_big)} ({100*len(lake_big)/len(conus_geoids):.0f}% of CONUS). Suppression rate: {s}/{t} = {r:.1f}%")
s, t, r = group_rate(nonlake_big)
summary_txt.append(f"  Non-lake-adjacent counties: {len(nonlake_big)}. Suppression rate: {s}/{t} = {r:.1f}%")
diff = r - group_rate(lake_big)[2]
summary_txt.append(f"  --> difference vs non-lake counties: {group_rate(lake_big)[2]-r:+.1f} pts")
print("\n".join(summary_txt), flush=True)

with open("/home/claude/suppression_map_summary.txt", "w") as f:
    f.write("\n".join(summary_txt) + "\n\n")
    f.write("Per-state waterbody counts (found -> matched any, matched >=1km2 lake/reservoir):\n")
    for abbr, (n, ma, mb) in sorted(state_lake_counts.items()):
        f.write(f"  {abbr}: {n} -> {ma}, {mb}\n")

print("=== 6. Rendering map ===", flush=True)

# simple equirectangular-ish projection for CONUS, cos-latitude corrected
REF_LAT = 38.0
def project(x, y):
    return x * math.cos(math.radians(REF_LAT)), y

fig, ax = plt.subplots(figsize=(16, 10), dpi=150)

def frac_to_color(frac, lake_flag):
    # never suppressed -> light blue-gray; partially -> orange gradient; fully -> dark red
    if frac == 0:
        base = (0.85, 0.88, 0.90)
    else:
        # interpolate light yellow -> dark red as frac goes 0->1
        base = (1.0 - 0.15 * frac, 0.75 - 0.65 * frac, 0.55 - 0.55 * frac)
    return base

for geoid, c in counties.items():
    abbr = STATE_FIPS_TO_ABBR.get(c["statefp"])
    if abbr not in CONUS_ABBR:
        continue
    frac = suppression_summary[geoid]["suppressed_frac"]
    color = frac_to_color(frac, geoid in lake_adjacent)
    for ring in c["rings"]:
        if len(ring) < 3:
            continue
        pts = [project(x, y) for x, y in ring]
        poly = mpatches.Polygon(pts, closed=True, facecolor=color,
                                 edgecolor="white", linewidth=0.05, zorder=1)
        ax.add_patch(poly)

# overlay lake-adjacent counties (Definition B: >=1 sq km LakePond/Reservoir,
# the size CyAN can plausibly resolve) with a bold outline
for geoid in lake_adjacent_big:
    c = counties.get(geoid)
    if c is None:
        continue
    abbr = STATE_FIPS_TO_ABBR.get(c["statefp"])
    if abbr not in CONUS_ABBR:
        continue
    for ring in c["rings"]:
        if len(ring) < 3:
            continue
        pts = [project(x, y) for x, y in ring]
        poly = mpatches.Polygon(pts, closed=True, facecolor="none",
                                 edgecolor="#1a5fb4", linewidth=0.7, zorder=2)
        ax.add_patch(poly)

ax.set_xlim(*[project(x, 0)[0] for x in (-125.5, -66.5)])
ax.set_ylim(23.5, 49.5)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title(
    "CDC WONDER suppression (N18/K74/K76, 2018-2024) by county, vs. NHDPlus lake-adjacency\n"
    f"Color = fraction of the 7 years suppressed  |  Blue outline = county has a NHDPlus lake/reservoir "
    f">= {BIG_LAKE_MIN_AREASQKM} sq km ({len(covered_abbrs)}/{len(CONUS_ABBR)} states; SC pending a geometry fix)",
    fontsize=11,
)

# legend
legend_handles = [
    mpatches.Patch(facecolor=frac_to_color(0, False), edgecolor="gray", label="0/7 years suppressed"),
    mpatches.Patch(facecolor=frac_to_color(0.5, False), edgecolor="gray", label="~3-4/7 years suppressed"),
    mpatches.Patch(facecolor=frac_to_color(1.0, False), edgecolor="gray", label="7/7 years suppressed"),
    mpatches.Patch(facecolor="none", edgecolor="#1a5fb4", linewidth=1.5, label=f"Lake/reservoir >= {BIG_LAKE_MIN_AREASQKM} km2 in county"),
]
ax.legend(handles=legend_handles, loc="lower right", fontsize=9, frameon=True)

footer_lines = [
    f"National suppression rate (all CONUS+DC counties, 2018-2024): {national_rate:.1f}%",
    f"Among counties with a >= {BIG_LAKE_MIN_AREASQKM} sq km lake/reservoir (n={len(lake_big)}): {group_rate(lake_big)[2]:.1f}%   |   Among counties without one (n={len(nonlake_big)}): {group_rate(nonlake_big)[2]:.1f}%",
]
fig.text(0.01, 0.01, "\n".join(footer_lines),
          fontsize=8, va="bottom", ha="left", family="monospace")

out_path = "/home/claude/suppression_lake_map.png"
fig.savefig(out_path, bbox_inches="tight", facecolor="white")
print(f"Saved map to {out_path}", flush=True)
