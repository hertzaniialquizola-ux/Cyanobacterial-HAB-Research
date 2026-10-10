#!/usr/bin/env python3
"""
sdwis_feasibility.py -- FEASIBILITY TEST of a drinking-water-weighted bloom exposure.

Question: can we replace "share of ALL lake pixels in the county >= 20,000 cells/mL" with an
exposure that reflects who actually drinks surface water? Test on a few states and two years.

Stage 1  per-LAKE bloom metric (same CyAN rules as build_alt_cei_v2.py: DN 254/255 removed,
         >=50% valid pixels, DN>=42 ~ 20,000 cells/mL), kept per lake instead of per county.
Stage 2  link each active community water system (EPA service-area polygon) to lakes
         (a) lakes that overlap the polygon, else (b) lakes within RADIUS_KM of it.
Stage 3  county exposure = population-weighted over surface-water systems; compare with the
         current county exposure (outputs/alt_cei_county_year.csv).

NOTHING here changes your existing outputs. Everything is written to outputs/sdwis_test/.

Run from the project root (Ohio+Iowa, 2021+2022 by default; re-run is resumable):
    python3 scripts/sdwis_feasibility.py
    STATES=39,19 YEARS=2021,2022 RADIUS_KM=25 python3 scripts/sdwis_feasibility.py

Honest limits (also in the summary): the polygon is where water is DELIVERED, not where the
intake is; "purchased surface water" (SWP) systems draw from a wholesaler's lake, so they are
reported separately; multi-county systems split population evenly across their counties.
"""
import os, re, sys, time, pickle
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
from shapely.geometry import box as shapely_box

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_county_year_panel as b

ROOT = b.ROOT
STATES = os.environ.get("STATES", "39,19").split(",")
YEARS = [int(y) for y in os.environ.get("YEARS", "2021,2022").split(",")]
RADIUS_KM = float(os.environ.get("RADIUS_KM", "25"))
OUT = ROOT / "outputs/sdwis_test"; OUT.mkdir(parents=True, exist_ok=True)
DN_20K, MIN_VALID_FRAC, MIN_WEEKS = 42, 0.5, 10
ACRS = b.ASSIGNMENT_CRS
SDWIS = ROOT / "data/raw/epa_sdwis"
tag = "_".join(STATES) + "_" + "_".join(map(str, YEARS))
log = []
def say(s=""):
    print(s, flush=True); log.append(s)

# ------------------------------------------------------------------ Stage 1: per-lake metrics
counties = b.load_counties()
cs = counties[counties.STATEFP.isin(STATES)].copy()
if cs.empty: sys.exit(f"No counties for STATES={STATES}")
reach = gpd.GeoDataFrame(geometry=cs.to_crs(ACRS).dissolve().geometry.buffer(RADIUS_KM * 1000), crs=ACRS)
nhd_path, _ = b.resolve_nhdplus_path()
lakes = gpd.read_file(nhd_path, mask=reach)
lakes = lakes[lakes["ftype"].isin({"LakePond", "Reservoir"}) & (lakes["areasqkm"] >= 1.0)].reset_index(drop=True)
lakes["lake_id"] = lakes.index
lakes_a = lakes.to_crs(ACRS)
lakes_a["geometry"] = lakes_a.geometry.buffer(0)
lakes_a["lake_km2"] = lakes_a.geometry.area / 1e6
say(f"STATES={STATES} YEARS={YEARS} RADIUS_KM={RADIUS_KM}")
say(f"Lakes (>=1 km2, LakePond/Reservoir) within the states + {RADIUS_KM:.0f} km: {len(lakes_a)}")

groups = {}
for tif in sorted(b.CYAN_DIR.glob("*_CYAN_CONUS_300m*.tif")):
    try: label, comp = b._parse_cyan_date(tif.name)
    except ValueError: continue
    if int(label[:4]) in YEARS:
        groups.setdefault((label, comp), []).append(tif)
say(f"CyAN composite groups to read: {len(groups)}")

CACHE = OUT / f"lake_extract_cache_{tag}.pkl"
recs, done_groups = [], set()
if CACHE.exists():
    recs, done_groups = pickle.load(open(CACHE, "rb"))
    say(f"Resuming: {len(done_groups)}/{len(groups)} groups already done")
proj, t0 = {}, time.time()
for i, ((label, comp), files) in enumerate(groups.items(), 1):
    if (label, comp) in done_groups: continue
    print(f"  group {i}/{len(groups)} ({label}) {time.time()-t0:.0f}s", flush=True)
    non_tiles = [f for f in files if not b._TILE_SUFFIX_RE.search(f.name)]
    tiles = [f for f in files if f not in non_tiles]
    done = set()
    for tif in non_tiles + tiles:
        with rasterio.open(tif) as src:
            key = str(src.crs)
            if key not in proj: proj[key] = lakes_a.to_crs(src.crs)
            lk = proj[key]
            cand = [p for p in lk.sindex.query(shapely_box(*src.bounds), predicate="intersects") if p not in done]
            if not cand: continue
            full = src.read(1); H, W = full.shape
            for pos in cand:
                geom = lk.geometry.iloc[pos]
                try:
                    win = from_bounds(*geom.bounds, transform=src.transform).round_offsets().round_lengths()
                    win = win.intersection(Window(0, 0, W, H))
                except Exception:
                    continue
                if win.width < 1 or win.height < 1: continue
                sub = full[win.toslices()]
                inside = geometry_mask([geom], out_shape=sub.shape, transform=src.window_transform(win), invert=True)
                px = sub[inside].astype(np.int16)
                if px.size == 0: continue
                valid = px[(px != 254) & (px != 255)]
                if valid.size / px.size < MIN_VALID_FRAC or valid.size == 0: continue
                recs.append((int(lk["lake_id"].iloc[pos]), label, comp, int(valid.size), int((valid >= DN_20K).sum())))
                done.add(pos)
            del full
    done_groups.add((label, comp))
    if len(done_groups) % 5 == 0: pickle.dump((recs, done_groups), open(CACHE, "wb"))
pickle.dump((recs, done_groups), open(CACHE, "wb"))

lw = pd.DataFrame(recs, columns=["lake_id", "period", "comp", "n_valid", "n_ge20k"])
lw["share"] = lw.n_ge20k / lw.n_valid
lw["year"] = lw.period.str.slice(0, 4).astype(int)
ly = lw.groupby(["lake_id", "year"]).agg(n_weeks=("share", "size"), lake_share_20k=("share", "mean")).reset_index()
ly = ly[ly.n_weeks >= MIN_WEEKS]
keep_cols = [c for c in ["lake_id", "permanent_identifier", "comid", "gnis_name", "areasqkm"] if c in lakes.columns]
ly = ly.merge(lakes[keep_cols], on="lake_id", how="left")
ly.to_csv(OUT / f"lake_year_metrics_{tag}.csv", index=False)
say(f"Lake-years with a usable metric (>= {MIN_WEEKS} valid weeks): {len(ly)} over {ly.lake_id.nunique()} lakes")
lake_metric = {(r.lake_id, r.year): r.lake_share_20k for r in ly.itertuples()}
lake_km2 = lakes_a.set_index("lake_id")["lake_km2"].to_dict()

# ------------------------------------------------------------------ Stage 2: systems -> lakes
sysd = pd.read_csv(SDWIS / "system_search/SDWA_SYSTEM_SEARCH.csv", dtype=str,
                   usecols=["PWSID", "PWS_NAME", "PWS_TYPE_CODE", "PWS_ACTIVITY_CODE", "PRIMARY_SOURCE_CODE",
                            "POPULATION_SERVED_COUNT", "FIPS_CODES"])
sysd = sysd[(sysd.PWS_TYPE_CODE == "CWS") & (sysd.PWS_ACTIVITY_CODE == "A")].copy()
sysd["pop"] = pd.to_numeric(sysd.POPULATION_SERVED_COUNT, errors="coerce").fillna(0)
sysd["fips"] = sysd.FIPS_CODES.fillna("").apply(lambda s: re.findall(r"\d{5}", s))
sysd = sysd[sysd.fips.apply(lambda f: any(x[:2] in STATES for x in f))].copy()
say(f"Active community water systems in these states: {len(sysd)}; "
    f"source mix (people, millions): " + str((sysd.groupby('PRIMARY_SOURCE_CODE')['pop'].sum() / 1e6).round(2).to_dict()))

gp = SDWIS / "3_0/Service_Areas_V_3_0.gpkg"
poly = gpd.read_file(gp, layer="CWS", mask=cs.dissolve()[["geometry"]])
pid_col = next((c for c in poly.columns if c.upper() == "PWSID"), None)
if pid_col is None: sys.exit(f"No PWSID column in CWS layer; columns = {list(poly.columns)}")
poly = poly[poly[pid_col].isin(sysd.PWSID)].dissolve(by=pid_col).reset_index().rename(columns={pid_col: "PWSID"})
poly = poly.to_crs(ACRS); poly["geometry"] = poly.geometry.buffer(0)
say(f"Service-area polygons matched to systems: {len(poly)} of {len(sysd)} "
    f"({poly.PWSID.isin(sysd[sysd.PRIMARY_SOURCE_CODE.isin(['SW','SWP'])].PWSID).sum()} surface-water)")

lake_geoms = lakes_a.geometry; sidx = lake_geoms.sindex
link = {}
sw = sysd[sysd.PRIMARY_SOURCE_CODE.isin(["SW", "SWP", "GU"])]
poly_sw = poly[poly.PWSID.isin(sw.PWSID)]
for r in poly_sw.itertuples():
    g = r.geometry
    cand = list(sidx.query(g.buffer(RADIUS_KM * 1000), predicate="intersects"))
    ins = [int(lakes_a.lake_id.iloc[i]) for i in cand if lake_geoms.iloc[i].intersects(g)]
    near = [int(lakes_a.lake_id.iloc[i]) for i in cand]
    link[r.PWSID] = (ins, near)

rows = []
for p, (ins, near) in link.items():
    for y in YEARS:
        def wmean(ids):
            v = [(lake_metric[(i, y)], lake_km2[i]) for i in ids if (i, y) in lake_metric]
            return np.average([a for a, _ in v], weights=[w for _, w in v]) if v else np.nan
        e_in = wmean(ins); e_near = wmean(near)
        rows.append(dict(PWSID=p, year=y, n_lakes_in=len(ins), n_lakes_near=len(near), E_in=e_in,
                         E_best=e_in if not np.isnan(e_in) else e_near,
                         link_type="in" if not np.isnan(e_in) else ("near" if not np.isnan(e_near) else "none")))
se = pd.DataFrame(rows).merge(sysd[["PWSID", "PWS_NAME", "PRIMARY_SOURCE_CODE", "pop"]], on="PWSID")
se.to_csv(OUT / f"system_exposure_{tag}.csv", index=False)
for s in ["SW", "SWP", "GU"]:
    t = se[(se.PRIMARY_SOURCE_CODE == s) & (se.year == YEARS[0])]
    if len(t):
        say(f"  {s}: {len(t)} systems, {t['pop'].sum()/1e6:.2f}M people | linked to a lake: "
            f"overlapping {(t.link_type=='in').mean():.0%}, within {RADIUS_KM:.0f} km {(t.link_type=='near').mean():.0%}, "
            f"none {(t.link_type=='none').mean():.0%} (people-weighted linked: {t.loc[t.link_type!='none','pop'].sum()/t['pop'].sum():.0%})")

# ------------------------------------------------------------------ Stage 3: county exposure
cty_rows = []
for r in sysd.itertuples():
    f = [x for x in r.fips if x[:2] in STATES]
    for c in f: cty_rows.append((c, r.PWSID, r.PRIMARY_SOURCE_CODE, r.pop / len(r.fips)))
cp = pd.DataFrame(cty_rows, columns=["county_fips", "PWSID", "src", "pop_c"])
tot = cp.groupby("county_fips").pop_c.sum().rename("pop_total")
def county_table(src_set, label):
    d = cp[cp.src.isin(src_set)].merge(se[["PWSID", "year", "E_best"]], on="PWSID")
    d = d.dropna(subset=["E_best"])
    d["w"] = d.E_best * d.pop_c
    o = d.groupby(["county_fips", "year"]).agg(sw=("w", "sum"), pl=("pop_c", "sum")).reset_index()
    o[f"E_linked_{label}"] = np.where(o.pl > 0, o.sw / o.pl, np.nan)
    o = o.rename(columns={"pl": f"pop_linked_{label}"}).drop(columns="sw")
    o = o.merge(tot.reset_index(), on="county_fips")
    o[f"E_dose_{label}"] = o[f"E_linked_{label}"] * o[f"pop_linked_{label}"] / o.pop_total
    return o
c1 = county_table(["SW"], "SW"); c2 = county_table(["SW", "SWP"], "SWSWP")
ct = c1.merge(c2.drop(columns="pop_total"), on=["county_fips", "year"], how="outer")
cur = pd.read_csv(ROOT / "outputs/alt_cei_county_year.csv", dtype={"county_fips": str})[["county_fips", "year", "px_share_20k"]]
cmp = cur[cur.county_fips.str[:2].isin(STATES) & cur.year.isin(YEARS)].merge(ct, on=["county_fips", "year"], how="outer")
swfrac = cp[cp.src.isin(["SW", "SWP"])].groupby("county_fips").pop_c.sum() / tot
cmp["surface_water_pop_frac"] = cmp.county_fips.map(swfrac)
cmp.to_csv(OUT / f"county_exposure_comparison_{tag}.csv", index=False)

say("")
say("=== COUNTY-LEVEL COMPARISON (current vs drinking-water-weighted) ===")
say(f"counties in these states: {cs.county_fips.nunique()} | with CURRENT exposure: {cur[cur.county_fips.str[:2].isin(STATES)].county_fips.nunique()} "
    f"| with NEW exposure (SW only): {c1.county_fips.nunique()} | (SW+SWP): {c2.county_fips.nunique()}")
for y in YEARS:
    m = cmp[(cmp.year == y)].dropna(subset=["px_share_20k", "E_linked_SW"])
    if len(m) > 5:
        say(f"{y}: n={len(m)} counties with both | Spearman current vs new(linked SW) = {m.px_share_20k.corr(m.E_linked_SW, method='spearman'):.2f}"
            f" | vs dose version = {m.px_share_20k.corr(m.E_dose_SW, method='spearman'):.2f}")
hot = cmp[(cmp.px_share_20k > 0.05) & cmp.surface_water_pop_frac.notna()]
if len(hot):
    say(f"Counties-years where CURRENT exposure shows blooms (>5% of lake pixels): {len(hot)}; "
        f"of those, {(hot.surface_water_pop_frac < 0.25).mean():.0%} have <25% of residents on surface water "
        f"(blooming lakes the population mostly does not drink from)")
say(f"Counties with current exposure but NO linked surface-water system: {cmp.px_share_20k.notna().sum() - cmp.dropna(subset=['px_share_20k','E_linked_SW']).shape[0]} county-years")
say("")
say("Files: lake_year_metrics_*.csv, system_exposure_*.csv, county_exposure_comparison_*.csv in outputs/sdwis_test/")
(OUT / f"feasibility_summary_{tag}.txt").write_text("\n".join(log))
