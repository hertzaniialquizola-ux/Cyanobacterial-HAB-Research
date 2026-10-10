#!/usr/bin/env python3
"""
sensitivity_reliability.py -- exposure-measurement-error SENSITIVITY table (not a correction).

Takes the observed estimates and rescales their log-RR and CI by the attenuation factors
already measured in sim_power_v2.R (outputs/model_final_v3_nolst/sim_mde.csv), giving the
range of TRUE per-IQR effects compatible with each result under each assumed error scenario.

    true log RR = observed log RR / attenuation        (same for the CI limits)

No new model fits and no new data. Needs only pandas + numpy. Run from the project root:
    python3 scripts/sensitivity_reliability.py
Writes outputs/model_final_v3_nolst/sensitivity_reliability.csv and .md (paper-ready table).

Caveats (also written into the output notes):
  * Primary (county-year) rows use the simulation's attenuation, which came from a NON-spatial
    fit; the observed estimate is the spatial BYM2 one. Treat as approximate.
  * Pooled kidney/liver rows were NOT simulated. They are cross-sectional (no county random
    effect), so attenuation is approximated by the reliability itself (classical-error rule of
    thumb). Persistent vs classical makes no difference there, so one row per reliability.
  * Reliability values (1.0, 0.5, 0.3) are the ones simulated; none is an estimate of the true
    reliability, which is unknown.
"""
import os, sys
import numpy as np, pandas as pd

D = "outputs/model_final_v3_nolst"
P = "outputs/model_pooled_nolst"
for f in (f"{D}/sim_mde.csv", f"{D}/results_v3.csv", f"{P}/results_pooled.csv"):
    if not os.path.exists(f):
        sys.exit(f"Missing {f}. Run from the project root.")

mde = pd.read_csv(f"{D}/sim_mde.csv")
v3 = pd.read_csv(f"{D}/results_v3.csv")
pooled = pd.read_csv(f"{P}/results_pooled.csv")

prim = v3[v3["model"] == "Baseline BYM2 primary"].iloc[0]
obs = [("Primary (county-year, BYM2)", prim.RR, prim.lo95, prim.hi95, "sim")]
for _, r in pooled.iterrows():
    if r.outcome in ("kidney", "liver"):
        obs.append((f"Pooled 2018-2024, {r.outcome}", r.RR_per_IQR, r.lo95, r.hi95, "approx"))

def rescale(rr, a):
    return float(np.exp(np.log(rr) / a))

rows = []
for name, rr, lo, hi, kind in obs:
    if kind == "sim":
        scen = [("none", 1.0, 1.0)]
        for _, m in mde.iterrows():
            if m.error_type != "none":
                scen.append((m.error_type, float(m.reliability), float(m.attenuation)))
    else:
        scen = [("none", 1.0, 1.0), ("any", 0.5, 0.5), ("any", 0.3, 0.3)]
    for etype, rel, a in scen:
        rows.append(dict(analysis=name, error_type=etype, reliability=rel,
                         attenuation=round(a, 3), observed_RR=round(rr, 3),
                         true_RR=round(rescale(rr, a), 3),
                         true_lo95=round(rescale(lo, a), 3), true_hi95=round(rescale(hi, a), 3),
                         basis="simulated attenuation (non-spatial sim)" if kind == "sim"
                               else "approximation: attenuation = reliability"))
out = pd.DataFrame(rows)
out["_o"] = out["error_type"].map({"none": 0, "classical": 1, "persistent": 2, "any": 1})
out = out.sort_values(["analysis", "_o", "reliability"], ascending=[True, True, False]).drop(columns="_o")
out.to_csv(f"{D}/sensitivity_reliability.csv", index=False)

# paper-ready markdown table
lines = ["| Analysis | Assumed error | Reliability | Attenuation | True RR per IQR (95% CI) |",
         "|---|---|---|---|---|"]
for _, r in out.iterrows():
    et = {"none": "none", "classical": "year-to-year (classical)", "persistent": "persistent within county",
          "any": "cross-sectional"}[r.error_type]
    lines.append(f"| {r.analysis} | {et} | {r.reliability:.1f} | {r.attenuation:.2f} | "
                 f"{r.true_RR:.2f} ({r.true_lo95:.2f}-{r.true_hi95:.2f}) |")
md = "\n".join(lines)
open(f"{D}/sensitivity_reliability.md", "w").write(md + "\n")
print(md)
print("\nUpper confidence limits under each scenario (largest true effect still compatible with the data):")
print(out[["analysis", "error_type", "reliability", "true_hi95"]].to_string(index=False))
print(f"\nWrote {D}/sensitivity_reliability.csv and .md")
