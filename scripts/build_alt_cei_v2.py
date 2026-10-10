"""
build_alt_cei_v2.py -- same corrected CyAN extraction as build_alt_cei.py, but with a
FIXED lake-to-county assignment.

WHY (critique 2026-10-03): the old code gave each lake wholly to the county containing the
lake's CENTROID. That (a) put the whole of a lake that straddles county lines into one county,
and (b) dropped any lake whose centroid is not inside a county polygon -- e.g. Lake Erie, whose
centroid is in open water, which is why Lucas County, OH (Toledo, the textbook bloom case)
had no exposure data.

NEW RULE: every lake is cut by the county boundaries (overlay), and each lake-in-county piece
is extracted separately, so a county gets exactly the lake pixels that lie inside it. Water that
lies outside every county polygon (Great Lakes, coastal edge, boundary generalization slivers)
is assigned to the county whose shoreline is nearest, but only within NEARSHORE_KM of that
county (default 5 km): a nearshore band, not the whole open lake. Where two counties' bands
overlap on the water, the strip goes to the lower county FIPS (a deliberate, documented
simplification; the strips are narrow).

Everything else (DN decoding, 254/255 removal, >=50% valid pixels, >=10 weeks/year) is unchanged.

Run from the project root (slow: it re-reads all rasters, like the original; resumable):
    python3 scripts/build_alt_cei_v2.py            # NEARSHORE_KM=5 by default
    NEARSHORE_KM=0 python3 scripts/build_alt_cei_v2.py     # overlay only, no nearshore band
Backs up the old outputs as *_centroid_assignment.csv, then overwrites
outputs/county_cyan_weekly_fixed.csv, outputs/alt_cei_county_year.csv, outputs/alt_cei_summary.txt
and writes outputs/lake_county_pieces_summary.csv (diagnostics).
Then re-run: assemble_panel_v2.py, then the R model scripts.
"""
import os, sys, shutil, pickle, time
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
from shapely.geometry import box as shapely_box

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_county_year_panel as b

ROOT = b.ROOT
DN_20K, DN_100K = 42, 100
MIN_VALID_FRAC, MIN_WEEKS = 0.5, 10
NEARSHORE_KM = float(os.environ.get("NEARSHORE_KM", "5"))
MIN_PIECE_KM2 = 0.01
ACRS = b.ASSIGNMENT_CRS            # EPSG:5070 (metres)

def dn_to_cells(dn):
    return 1e8 * 10 ** (3.0 / 250.0 * dn - 4.2)
assert abs(dn_to_cells(0) - 6310) < 50
assert dn_to_cells(DN_20K) >= 20000 > dn_to_cells(DN_20K - 1)


def build_lake_pieces(counties, nhd_path):
    lakes = gpd.read_file(nhd_path)
    n0 = len(lakes)
    lakes = lakes[lakes["ftype"].isin(b.BIG_LAKE_FTYPES if hasattr(b, "BIG_LAKE_FTYPES") else {"LakePond", "Reservoir"})
                  & (lakes["areasqkm"] >= 1.0)].copy()
    print(f"Lakes: kept {len(lakes)} of {n0} NHDPlus waterbodies (LakePond/Reservoir, >= 1 km2)", flush=True)
    lakes = lakes.to_crs(ACRS)
    lakes["geometry"] = lakes.geometry.buffer(0)
    lakes = lakes.reset_index(drop=True)
    lakes["lake_id"] = lakes.index
    lakes["lake_km2"] = lakes.geometry.area / 1e6
    cty = counties.to_crs(ACRS)[["county_fips", "geometry"]].copy()
    cty["geometry"] = cty.geometry.buffer(0)

    # 1. inside-county pieces
    inside = gpd.overlay(lakes[["lake_id", "geometry"]], cty, how="intersection", keep_geom_type=True)
    inside["km2"] = inside.geometry.area / 1e6
    inside = inside[inside["km2"] >= MIN_PIECE_KM2]
    inside["kind"] = "inside"
    covered = inside.groupby("lake_id")["km2"].sum()
    lakes["covered_frac"] = lakes["lake_id"].map(covered).fillna(0) / lakes["lake_km2"]
    print(f"  {int((lakes.covered_frac < 0.995).sum())} lakes have > 0.5% of their area outside every county polygon "
          f"(Great Lakes, coast, generalization slivers)", flush=True)

    # 2. nearshore band for the part of each lake outside every county
    pieces = [inside[["lake_id", "county_fips", "geometry", "kind"]]]
    if NEARSHORE_KM > 0:
        buf = cty.copy(); buf["geometry"] = cty.geometry.buffer(NEARSHORE_KM * 1000)
        csi, bsi = cty.sindex, buf.sindex
        rows = []
        todo = lakes[lakes.covered_frac < 0.995]
        for k, lk in enumerate(todo.itertuples(), 1):
            geom = lk.geometry
            near_land = cty.geometry.iloc[csi.query(geom.buffer(NEARSHORE_KM * 1000), predicate="intersects")]
            left = geom.difference(near_land.union_all() if hasattr(near_land, "union_all") else near_land.unary_union) if len(near_land) else geom
            if left.is_empty:
                continue
            cand = sorted(bsi.query(left, predicate="intersects"), key=lambda i: cty.county_fips.iloc[i])
            for ci in cand:
                p = left.intersection(buf.geometry.iloc[ci])
                if p.is_empty or p.area / 1e6 < MIN_PIECE_KM2:
                    continue
                rows.append((lk.lake_id, cty.county_fips.iloc[ci], p, "nearshore"))
                left = left.difference(p)
                if left.is_empty:
                    break
            if k % 200 == 0:
                print(f"    nearshore pieces: {k}/{len(todo)} lakes", flush=True)
        if rows:
            pieces.append(gpd.GeoDataFrame(rows, columns=["lake_id", "county_fips", "geometry", "kind"], crs=ACRS))
    pc = pd.concat(pieces, ignore_index=True)
    pc = gpd.GeoDataFrame(pc, geometry="geometry", crs=ACRS).explode(index_parts=False).reset_index(drop=True)
    pc["km2"] = pc.geometry.area / 1e6
    pc = pc[pc["km2"] >= MIN_PIECE_KM2].reset_index(drop=True)
    return pc, lakes


counties = b.load_counties()
nhd_path, _ = b.resolve_nhdplus_path()
lakes, lake_tbl = build_lake_pieces(counties, nhd_path)
lakes.index = range(len(lakes))

# ---- diagnostics: how many counties gain/lose lake area vs the centroid rule -------------
summ = lakes.groupby(["county_fips", "kind"])["km2"].sum().unstack(fill_value=0).reset_index()
summ.to_csv(ROOT / "outputs/lake_county_pieces_summary.csv", index=False)
print(f"{len(lakes)} lake-in-county pieces in {lakes.county_fips.nunique()} counties "
      f"({int((lakes.kind == 'nearshore').sum())} nearshore pieces)", flush=True)
for fips, nm in [("39095", "Lucas OH (Toledo)"), ("13177", "Lee GA")]:
    s = lakes[lakes.county_fips == fips]
    print(f"  {nm}: {len(s)} pieces, {s.km2.sum():.1f} km2 of lake", flush=True)

all_tifs = sorted(b.CYAN_DIR.glob("*_CYAN_CONUS_300m*.tif"))
groups = {}
for tif in all_tifs:
    try:
        label, comp = b._parse_cyan_date(tif.name)
    except ValueError:
        continue
    groups.setdefault((label, comp), []).append(tif)

CACHE = ROOT / "outputs/_cyan_extract_cache_v2.pkl"       # NEW cache: the old one used centroid assignment
recs, done_groups = [], set()
if CACHE.exists():
    recs, done_groups = pickle.load(open(CACHE, "rb"))
    print(f"Resuming: {len(done_groups)} of {len(groups)} groups already done.", flush=True)
t0 = time.time()
for i, ((label, comp), files) in enumerate(groups.items(), 1):
    if (label, comp) in done_groups:
        continue
    print(f"  group {i}/{len(groups)} ({label}) - {time.time()-t0:.0f}s elapsed", flush=True)
    non_tiles = [f for f in files if not b._TILE_SUFFIX_RE.search(f.name)]
    tiles = [f for f in files if f not in non_tiles]
    done = set()
    for tif in non_tiles + tiles:
        with rasterio.open(tif) as src:
            lk = lakes.to_crs(src.crs) if lakes.crs != src.crs else lakes
            cand = list(lk.sindex.query(shapely_box(*src.bounds), predicate="intersects"))
            cand = [p for p in cand if lk.index[p] not in done]
            if not cand:
                continue
            full = src.read(1)
            H, W = full.shape
            for pos in cand:
                pid = lk.index[pos]
                geom = lk.geometry.iloc[pos]
                try:
                    win = from_bounds(*geom.bounds, transform=src.transform).round_offsets().round_lengths()
                    win = win.intersection(Window(0, 0, W, H))
                except Exception:
                    continue
                if win.width < 1 or win.height < 1:
                    continue
                sub = full[win.toslices()]
                inside = geometry_mask([geom], out_shape=sub.shape, transform=src.window_transform(win), invert=True)
                px = sub[inside].astype(np.int16)
                if px.size == 0:
                    continue
                valid = px[(px != 254) & (px != 255)]
                if valid.size / px.size < MIN_VALID_FRAC or valid.size == 0:
                    continue
                recs.append((lk["county_fips"].iloc[pos], label, comp, int(valid.size), int(px.size),
                             int((valid >= DN_20K).sum()), int((valid >= DN_100K).sum()),
                             float(np.log10(dn_to_cells(valid.astype(float))).mean()),
                             float(dn_to_cells(valid.astype(float)).mean())))
                done.add(pid)
            del full
    done_groups.add((label, comp))
    if len(done_groups) % 5 == 0:
        pickle.dump((recs, done_groups), open(CACHE, "wb"))
pickle.dump((recs, done_groups), open(CACHE, "wb"))

lw = pd.DataFrame(recs, columns=["county_fips", "period", "composite", "n_valid", "n_total",
                                 "n_ge20k", "n_ge100k", "sum_log10_mean", "mean_cells"])
lw["w_log10"] = lw["sum_log10_mean"] * lw["n_valid"]
lw["w_cells"] = lw["mean_cells"] * lw["n_valid"]
cw = lw.groupby(["county_fips", "period", "composite"], as_index=False).agg(
    n_valid=("n_valid", "sum"), n_ge20k=("n_ge20k", "sum"), n_ge100k=("n_ge100k", "sum"),
    w_log10=("w_log10", "sum"), w_cells=("w_cells", "sum"))
cw["share_20k"] = cw["n_ge20k"] / cw["n_valid"]
cw["share_100k"] = cw["n_ge100k"] / cw["n_valid"]
cw["log10_cells"] = cw["w_log10"] / cw["n_valid"]
cw["mean_cells"] = cw["w_cells"] / cw["n_valid"]
cw["wk_above_who"] = (cw["mean_cells"] >= 20000).astype(float)
cw["year"] = cw["period"].str.slice(0, 4).astype(int)

for f in ["county_cyan_weekly_fixed.csv", "alt_cei_county_year.csv"]:      # keep the old outputs
    src_f = ROOT / "outputs" / f
    if src_f.exists():
        shutil.copy(src_f, ROOT / "outputs" / f.replace(".csv", "_centroid_assignment.csv"))
cw.to_csv(ROOT / "outputs/county_cyan_weekly_fixed.csv", index=False)

g = cw.groupby(["county_fips", "year"])
out = g.agg(n_weeks_valid=("share_20k", "size"), cei=("wk_above_who", "mean"),
            px_share_20k=("share_20k", "mean"), px_share_100k=("share_100k", "mean"),
            log10_cells_mean=("log10_cells", "mean"),
            wk_frac_5pct_20k=("share_20k", lambda x: (x >= 0.05).mean())).reset_index()
out = out[out["n_weeks_valid"] >= MIN_WEEKS]
out.to_csv(ROOT / "outputs/alt_cei_county_year.csv", index=False)

old = ROOT / "outputs/alt_cei_county_year_centroid_assignment.csv"
lines = [f"NEARSHORE_KM={NEARSHORE_KM}", f"county-years kept (>= {MIN_WEEKS} valid weeks): {len(out)}, counties: {out.county_fips.nunique()}"]
if old.exists():
    o = pd.read_csv(old, dtype={"county_fips": str})
    lines.append(f"previous (centroid rule): {len(o)} county-years, {o.county_fips.nunique()} counties; "
                 f"counties gained: {len(set(out.county_fips) - set(o.county_fips))}, lost: {len(set(o.county_fips) - set(out.county_fips))}")
lines.append("Lucas OH (39095) county-years: %d; Lee GA (13177): %d" % ((out.county_fips == "39095").sum(), (out.county_fips == "13177").sum()))
for c in ["cei", "px_share_20k", "px_share_100k", "wk_frac_5pct_20k"]:
    lines.append(f"{c:18s} ==0: {(out[c]==0).mean():6.1%}  >=0.99: {(out[c]>=0.99).mean():6.1%}  "
                 f"median {out[c].median():.3f}  IQR {out[c].quantile(.25):.3f}-{out[c].quantile(.75):.3f}")
txt = "\n".join(lines)
(ROOT / "outputs/alt_cei_summary.txt").write_text(txt)
print(txt)
