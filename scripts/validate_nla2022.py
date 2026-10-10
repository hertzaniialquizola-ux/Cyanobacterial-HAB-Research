"""
validate_nla2022.py -- does the satellite cyanobacteria index track measured microcystin?

Matches each EPA NLA 2022 lake visit (microcystin, ug/L) to the CyAN 7-day composite that covers the
sampling date (and, as a sensitivity, the composites one week either side), using the pixels within
RADIUS_M of the NLA sample point. NLA 2017 cannot be matched here because the extraction covers 2018-2024.

Outputs (outputs/validation/):
  nla2022_matched.csv       one row per NLA visit with satellite metrics and microcystin
  nla2022_validation.txt    summary: n, detection rates, Spearman, AUC, and a county-level check
Run from the project root:  python scripts/validate_nla2022.py
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.windows import Window
from pyproj import Transformer
sys.path.insert(0, "scripts")
import build_county_year_panel as b

ROOT = b.ROOT; OUT = ROOT / "outputs/validation"; OUT.mkdir(parents=True, exist_ok=True)
NLA = ROOT / "data/raw/nla/nla2022"; RADIUS_M = float(__import__("os").environ.get("RADIUS_M", "1000")); MIN_PIX = 3
dn_cells = lambda dn: 1e8 * 10 ** (3.0 / 250.0 * dn - 4.2)

tox = pd.read_csv(NLA / "nla22_algalToxins.csv"); tox = tox[(tox.ANALYTE == "MICX") & (tox.VISIT_NO == 1)].copy()
tox["micx"] = pd.to_numeric(tox.RESULT, errors="coerce").fillna(0.0)          # ND = below detection -> 0
tox["detect"] = (tox.micx > 0).astype(int)
site = pd.read_csv(NLA / "nla22_siteinfo.csv"); site = site[site.VISIT_NO == 1]
site["lat"] = site.INDEX_LAT_DD.fillna(site.LAT_DD83); site["lon"] = site.INDEX_LON_DD.fillna(site.LON_DD83)
d = tox.merge(site[["UID", "lat", "lon", "AREA_HA", "CNTYNAME", "PSTL_CODE"]].drop_duplicates("UID"), on="UID", how="inner",
              suffixes=("", "_s"))
d["date"] = pd.to_datetime(d.DATE_COL); d["doy"] = d.date.dt.dayofyear
d = d.dropna(subset=["lat", "lon"]).reset_index(drop=True)
print(f"{len(d)} NLA 2022 index visits with microcystin and coordinates; {int(d.detect.sum())} detects, {(d.micx>=8).sum()} >= 8 ug/L")

rx = re.compile(r"L(\d{4})(\d{3})(\d{4})(\d{3})\.")
files = {}
for tif in sorted(b.CYAN_DIR.glob("L2022*_CYAN_CONUS_300m*.tif")):
    m = rx.match(tif.name)
    if m: files.setdefault((int(m.group(2)), int(m.group(4))), []).append(tif)
print(f"{len(files)} CyAN 2022 composites")
starts = sorted(files)

def metrics(px):
    v = px[(px != 254) & (px != 255)]
    if v.size < MIN_PIX: return None
    c = dn_cells(v.astype(float))
    return dict(n_valid=int(v.size), share20k=float((v >= 42).mean()), mean_cells=float(c.mean()),
                log10_cells=float(np.log10(c).mean()), dn_max=int(v.max()), dn_p90=float(np.percentile(v, 90)))

def read_point(tifs, lon, lat):
    for tif in tifs:                                    # merged mosaic first, tiles as fallback
        with rasterio.open(tif) as src:
            x, y = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True).transform(lon, lat)
            if not (src.bounds.left <= x <= src.bounds.right and src.bounds.bottom <= y <= src.bounds.top): continue
            r, c = src.index(x, y); h = int(np.ceil(RADIUS_M / abs(src.res[0])))
            win = Window(c - h, r - h, 2 * h + 1, 2 * h + 1).intersection(Window(0, 0, src.width, src.height))
            if win.width < 1 or win.height < 1: continue
            arr = src.read(1, window=win)
            rr, cc = np.mgrid[win.row_off:win.row_off + arr.shape[0], win.col_off:win.col_off + arr.shape[1]]
            circ = np.hypot((rr - r) * abs(src.res[1]), (cc - c) * abs(src.res[0])) <= RADIUS_M
            m = metrics(arr[circ])
            if m: return m
    return None

rows = []
for i, rw in d.iterrows():
    k = [j for j, (s, e) in enumerate(starts) if s <= rw.doy <= e]
    if not k: continue
    j = k[0]; res = {}
    for lab, jj in [("same", [j]), ("pm1wk", [j - 1, j, j + 1])]:
        ms = [read_point(files[starts[q]], rw.lon, rw.lat) for q in jj if 0 <= q < len(starts)]
        ms = [m for m in ms if m]
        if ms:
            res[lab] = dict(share20k=np.mean([m["share20k"] for m in ms]), log10_cells=np.mean([m["log10_cells"] for m in ms]),
                            dn_max=max(m["dn_max"] for m in ms), mean_cells=np.mean([m["mean_cells"] for m in ms]), n_valid=sum(m["n_valid"] for m in ms))
    rec = dict(UID=rw.UID, SITE_ID=rw.SITE_ID, date=rw.date.date(), lat=rw.lat, lon=rw.lon, area_ha=rw.AREA_HA, micx=rw.micx, detect=rw.detect)
    for lab, m in res.items():
        for kk, vv in m.items(): rec[f"{lab}_{kk}"] = vv
    rows.append(rec)
    if (i + 1) % 100 == 0: print(f"  {i+1}/{len(d)}", flush=True)
M = pd.DataFrame(rows); M.to_csv(OUT / "nla2022_matched.csv", index=False)

from scipy.stats import spearmanr
def auc(y, s):
    y = np.asarray(y); s = np.asarray(s); pos, neg = s[y == 1], s[y == 0]
    if len(pos) < 3 or len(neg) < 3: return np.nan
    return float((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean())
L = []
def say(s): print(s); L.append(s)
say(f"NLA 2022 visits: {len(d)}; matched to a CyAN composite with >= {MIN_PIX} valid pixels within {RADIUS_M:.0f} m: same-week {M.get('same_share20k', pd.Series(dtype=float)).notna().sum()}, +/-1 week {M.get('pm1wk_share20k', pd.Series(dtype=float)).notna().sum()}")
for lab in ["same", "pm1wk"]:
    c = f"{lab}_share20k"
    if c not in M: continue
    g = M.dropna(subset=[c]).copy()
    say(f"\n== {lab}: n={len(g)} lakes; microcystin detected {g.detect.mean()*100:.0f}%, >=1 ug/L {(g.micx>=1).mean()*100:.0f}%, >=8 ug/L {(g.micx>=8).sum()} lakes")
    for var in [f"{lab}_share20k", f"{lab}_log10_cells", f"{lab}_dn_max"]:
        rho = spearmanr(g[var], g.micx)
        say(f"   {var:22s} Spearman vs microcystin: {rho.correlation:+.2f} (p={rho.pvalue:.3g}); AUC detect {auc(g.detect, g[var]):.2f}; AUC >=1 ug/L {auc((g.micx>=1).astype(int), g[var]):.2f}; AUC >=8 ug/L {auc((g.micx>=8).astype(int), g[var]):.2f}")
    r = np.corrcoef(g[f"{lab}_log10_cells"], np.log10(g.micx.where(g.micx > 0, 0.05)))[0, 1]
    say(f"   Pearson r (log10 cells vs log10 microcystin, ND=0.05): {r:.2f}  (r^2 = {r*r:.2f})")
    hi = g[g[f"{lab}_share20k"] >= 0.2]; lo = g[g[f"{lab}_share20k"] < 0.02]
    say(f"   lakes with share>=20k >=20% of pixels: n={len(hi)}, median microcystin {hi.micx.median():.2f}, >=1 ug/L {(hi.micx>=1).mean()*100:.0f}%;  share<2%: n={len(lo)}, median {lo.micx.median():.2f}, >=1 ug/L {(lo.micx>=1).mean()*100:.0f}%")

# county-level check: the paper's exposure is a county-year mean, so compare the 2022 county exposure with the mean
# microcystin of NLA lakes sampled in that county
try:
    cty = b.load_counties().to_crs(4326)[["county_fips", "geometry"]]
    cty["geometry"] = cty.geometry.buffer(0)               # repair invalid rings before the point-in-polygon join
    pts = gpd.GeoDataFrame(M, geometry=gpd.points_from_xy(M.lon, M.lat), crs=4326)
    pts = gpd.sjoin(pts, cty, how="left", predicate="within")
    ex = pd.read_csv(ROOT / "outputs/county_year_panel_v2.csv", dtype={"county_fips": str})
    ex = ex[ex.year == 2022][["county_fips", "px_share_20k"]].dropna()
    g = pts.groupby("county_fips").agg(n=("micx", "size"), micx_mean=("micx", "mean"), micx_max=("micx", "max"), detect=("detect", "mean")).reset_index().merge(ex, on="county_fips")
    say(f"\n== county level: {pts.county_fips.nunique()} counties had NLA 2022 lakes; {len(g)} of them are in the analytic sample with a 2022 exposure")
    if len(g) >= 8:
        for v in ["micx_mean", "micx_max", "detect"]:
            rho = spearmanr(g.px_share_20k, g[v]); say(f"   county 2022 px_share_20k vs NLA {v}: Spearman {rho.correlation:+.2f} (p={rho.pvalue:.3g}, n={len(g)})")
    g.to_csv(OUT / "nla2022_county_check.csv", index=False)
except Exception as e:
    say(f"\n(county-level check skipped: {e})")
say("\nCaveats: single-date grab samples vs weekly satellite composites; sample point may be near shore where pixels are invalid; "
    "this measures agreement at the lake-week level, not the reliability of a county-year mean, so treat it as a guide for choosing reliability values in the simulation.")
(OUT / "nla2022_validation.txt").write_text("\n".join(L))
