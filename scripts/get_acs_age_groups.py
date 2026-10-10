"""Pull full ACS 2022 5-year age distribution (B01001) by county -> 9 broad age groups for
indirect age standardization.  Needs: export CENSUS_API_KEY="<your key>"
Writes data/raw/census/acs_age_groups_2022.csv  (county_fips, ag_0_14, ag_15_24, ... ag_85p)"""
import os, sys, json, urllib.request, csv
key = os.environ.get("CENSUS_API_KEY")
if not key: sys.exit('Set CENSUS_API_KEY first: export CENSUS_API_KEY="<key>"')
male = list(range(3, 26)); female = list(range(27, 50))        # B01001_003..025 (male), 027..049 (female)
vars_ = [f"B01001_{i:03d}E" for i in male + female]
url = ("https://api.census.gov/data/2022/acs/acs5?get=" + ",".join(vars_) +
       "&for=county:*&key=" + key)
with urllib.request.urlopen(url, timeout=120) as r: rows = json.load(r)
hdr, data = rows[0], rows[1:]
ix = {h: i for i, h in enumerate(hdr)}
# groups in terms of male variable numbers (female = +24)
groups = {"ag_0_14": [3, 4, 5], "ag_15_24": [6, 7, 8, 9, 10], "ag_25_34": [11, 12], "ag_35_44": [13, 14],
          "ag_45_54": [15, 16], "ag_55_64": [17, 18, 19], "ag_65_74": [20, 21, 22], "ag_75_84": [23, 24], "ag_85p": [25]}
os.makedirs("data/raw/census", exist_ok=True)
out = "data/raw/census/acs_age_groups_2022.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["county_fips"] + list(groups))
    for r in data:
        vals = []
        for g, nums in groups.items():
            tot = 0
            for n in nums:
                for k in (n, n + 24):
                    v = r[ix[f"B01001_{k:03d}E"]]; tot += int(v) if v not in (None, "", "-666666666") else 0
            vals.append(tot)
        w.writerow([r[ix["state"]] + r[ix["county"]]] + vals)
print(f"{len(data)} counties -> {out}")
