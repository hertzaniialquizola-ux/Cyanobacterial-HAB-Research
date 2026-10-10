"""
get_usgs_surface_supply.py -- county measure of how much drinking water comes from
SURFACE water (lake/reservoir/river intakes), the pathway cyanotoxins use to reach people.

USGS 2015 reports population-by-source (PS-SWPop) for only ~1,900 counties; 15 states
(AL, CA, ID, IN, IA, KS, KY, MI, NE, NV, ND, SD, TX, UT, WA) have '--' (not reported).
Public-supply WITHDRAWALS by source (PS-WSWFr surface, PS-WGWFr ground) ARE reported
for every county, so the main variable is withdrawal-based:

  surface_supply_fraction = (PS-TOPop / TP-TotPop) * PS-WSWFr / (PS-WSWFr + PS-WGWFr)
     = (share of residents on public supply) x (share of public-supply water from surface sources)

Caveats (state in the paper): withdrawals are attributed to the county where the intake
sits, not where the water is delivered; this is a share of water, not measured exposure.
A population-based version (PS-SWPop/TP-TotPop) is also written for the ~1,900 counties
that have it, and the script prints how well the two agree.
Usage: python3 scripts/get_usgs_surface_supply.py
Writes: data/raw/usgs/surface_supply_by_county.csv
"""
import csv
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
rd = csv.reader(open(ROOT / "data/raw/usgs/usco2015_raw.txt", encoding="latin1", newline=""))
next(rd); h = next(rd); ix = {c: i for i, c in enumerate(h)}
need = ["FIPS", "PS-SWPop", "PS-TOPop", "PS-WSWFr", "PS-WGWFr", "TP-TotPop"]
assert all(c in ix for c in need), f"missing columns; have {h[:20]}"
def num(r, c):
    try: return float(r[ix[c]])
    except ValueError: return None
rows, a, b = [], [], []
for r in rd:
    if len(r) < len(h) or not r[ix["FIPS"]].strip():
        continue
    tp, top, wsw, wgw, sw = (num(r, "TP-TotPop"), num(r, "PS-TOPop"), num(r, "PS-WSWFr"), num(r, "PS-WGWFr"), num(r, "PS-SWPop"))
    if not tp or tp <= 0 or top is None or wsw is None or wgw is None:
        continue
    frac_ps = min(top / tp, 1.0)
    share_sw = (wsw / (wsw + wgw)) if (wsw + wgw) > 0 else 0.0
    f_w = frac_ps * share_sw
    f_p = min(sw / tp, 1.0) if sw is not None else ""
    rows.append((r[ix["FIPS"]].zfill(5), f_w, f_p))
    if f_p != "":
        a.append(f_w); b.append(f_p)
out = ROOT / "data/raw/usgs/surface_supply_by_county.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["county_fips", "surface_supply_fraction", "surface_supply_fraction_pop"]); w.writerows(rows)
print(f"{len(rows)} counties written -> {out}")
print("mean surface_supply_fraction (withdrawal-based) = %.3f" % np.mean([r[1] for r in rows]))
print(f"agreement with population-based version on {len(a)} counties: r = {np.corrcoef(a, b)[0,1]:.3f}")
