"""
build_alt_cei.py -- CORRECTED CyAN extraction + continuous exposure metrics.

WHY THIS EXISTS (found 2026-10-02): the CyAN GeoTIFFs are uint8 digital numbers
(DN), with NO nodata tag set. DN 255 = no data/cloud (~252M pixels in the test
file), DN 254 = land (~133M), DN 0 = below detection, DN 1-253 = detected.
build_county_year_panel.py masks with src.nodata (=None), so 254/255 pixels
inside lake polygons were averaged in as if they were real values, and the
pixels outside the polygon were filled with 0. It then compared that mean DN
to 10 -- but DN 10 is only ~8,300 cells/mL, and WHO's 20,000 cells/mL (the
10 ug/L chl-a-equivalent tier) is DN ~42.
Decoding (EPA/NASA CyAN):  CI = 10^(3/250 * DN - 4.2);  cells/mL = CI * 1e8.
(DN 0 -> 6,310 cells/mL is the documented detection limit, a useful check.)

This script re-extracts properly:
  * only pixels inside the lake polygon are used (masked array, filled=False)
  * DN 254/255 are dropped; DN 0 counts as a valid "not detected" pixel (0 cells/mL)
  * county-week values are pixel-count-weighted across lakes (= area weighted)
  * lake-weeks with < MIN_VALID_FRAC valid pixels (cloud) are discarded
  * metrics are computed only over weeks with valid data

Run from project root:   python3 scripts/build_alt_cei.py     (slow: reads all rasters)
Writes: outputs/county_cyan_weekly_fixed.csv, outputs/alt_cei_county_year.csv,
        outputs/alt_cei_summary.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
import rasterio
from shapely.geometry import box as shapely_box

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_county_year_panel as b

ROOT = b.ROOT
DN_20K = 42        # >= 20,000 cells/mL (WHO moderate/10 ug/L chl-a tier)
DN_100K = 100      # >= 100,000 cells/mL (WHO high tier)
MIN_VALID_FRAC = 0.5
MIN_WEEKS = 10

def dn_to_cells(dn):
    return 1e8 * 10 ** (3.0 / 250.0 * dn - 4.2)

assert abs(dn_to_cells(0) - 6310) < 50, "DN->cells/mL decode sanity check failed"
assert dn_to_cells(DN_20K) >= 20000 > dn_to_cells(DN_20K - 1)

counties = b.load_counties()
nhd_path, _ = b.resolve_nhdplus_path()
lakes = b.load_lakes_with_county(counties, nhd_path)

all_tifs = sorted(b.CYAN_DIR.glob("*_CYAN_CONUS_300m*.tif"))
groups = {}
for tif in all_tifs:
    try:
        label, comp = b._parse_cyan_date(tif.name)
    except ValueError:
        continue
    groups.setdefault((label, comp), []).append(tif)

import pickle, time
CACHE = ROOT / "outputs/_cyan_extract_cache.pkl"
recs, done_groups = [], set()
if CACHE.exists():                       # resume after Ctrl-C / crash
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
            # SPEED FIX: these GeoTIFFs are stored as ONE strip, so every windowed read
            # (rio_mask per lake) decompresses the entire 396-megapixel raster (~30 s each
            # time). Read the band into memory ONCE per file, then cut lake windows from it.
            full = src.read(1)
            H, W = full.shape
            for pos in cand:
                lid = lk.index[pos]
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
                px = sub[inside].astype(np.int16)          # only pixels inside the polygon
                if px.size == 0:
                    continue
                valid = px[(px != 254) & (px != 255)]
                if valid.size / px.size < MIN_VALID_FRAC or valid.size == 0:
                    continue
                recs.append((lk["county_fips"].iloc[pos], label, comp, int(valid.size), int(px.size),
                             int((valid >= DN_20K).sum()), int((valid >= DN_100K).sum()),
                             float(np.log10(dn_to_cells(valid.astype(float))).mean()),
                             float(dn_to_cells(valid.astype(float)).mean())))
                done.add(lid)
            del full
    done_groups.add((label, comp))
    if len(done_groups) % 5 == 0:
        pickle.dump((recs, done_groups), open(CACHE, "wb"))

lw = pd.DataFrame(recs, columns=["county_fips", "period", "composite", "n_valid", "n_total",
                                 "n_ge20k", "n_ge100k", "sum_log10_mean", "mean_cells"])
# n_valid-weighted county-week aggregation
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
cw.to_csv(ROOT / "outputs/county_cyan_weekly_fixed.csv", index=False)

g = cw.groupby(["county_fips", "year"])
out = g.agg(n_weeks_valid=("share_20k", "size"),
            cei=("wk_above_who", "mean"),      # corrected CEI: fraction of valid weeks with county-mean >= 20,000 cells/mL
            px_share_20k=("share_20k", "mean"),
            px_share_100k=("share_100k", "mean"),
            log10_cells_mean=("log10_cells", "mean"),
            wk_frac_5pct_20k=("share_20k", lambda x: (x >= 0.05).mean())).reset_index()
out = out[out["n_weeks_valid"] >= MIN_WEEKS]
out.to_csv(ROOT / "outputs/alt_cei_county_year.csv", index=False)

lines = [f"county-years kept (>= {MIN_WEEKS} valid weeks): {len(out)}"]
for c in ["cei", "px_share_20k", "px_share_100k", "wk_frac_5pct_20k"]:
    lines.append(f"{c:18s} ==0: {(out[c]==0).mean():6.1%}  >=0.99: {(out[c]>=0.99).mean():6.1%}  "
                 f"median {out[c].median():.3f}  IQR {out[c].quantile(.25):.3f}-{out[c].quantile(.75):.3f}")
lines.append("log10_cells_mean:\n" + out["log10_cells_mean"].describe().to_string())
txt = "\n".join(lines)
(ROOT / "outputs/alt_cei_summary.txt").write_text(txt)
print(txt)
