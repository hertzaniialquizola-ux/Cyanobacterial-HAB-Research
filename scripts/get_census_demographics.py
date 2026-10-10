"""
get_census_demographics.py -- ACS 5-year (2022) county age/sex/race-ethnicity
covariates for the final model (age is the dominant driver of CKD/cirrhosis
mortality; this is a covariate-level adjustment, NOT true age standardization).

  B01001_001E total pop; _002E male; _020E-_025E male 65+; _044E-_049E female 65+
  B03002_001E total; _003E NH White; _004E NH Black; _005E NH American Indian/Alaska Native; _012E Hispanic

Usage:  export CENSUS_API_KEY="<key>";  python3 scripts/get_census_demographics.py
Writes: data/raw/census/demographics_by_county_2022.json
"""
import json, os, sys, urllib.request, urllib.error
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/raw/census/demographics_by_county_2022.json"
key = os.environ.get("CENSUS_API_KEY") or sys.exit('Set CENSUS_API_KEY first: export CENSUS_API_KEY="<key>"')
vars_ = (["B01001_001E", "B01001_002E"] + [f"B01001_0{n}E" for n in range(20, 26)] +
         [f"B01001_0{n}E" for n in range(44, 50)] + ["B03002_001E", "B03002_003E", "B03002_004E", "B03002_005E", "B03002_012E"])
url = f"https://api.census.gov/data/2022/acs/acs5?get=NAME,{','.join(vars_)}&for=county:*&key={key}"
try:
    data = json.loads(urllib.request.urlopen(url, timeout=90).read().decode())
except urllib.error.HTTPError as e:
    sys.exit(f"Census API failed: HTTP {e.code} {e.read().decode()[:400]}")
OUT.parent.mkdir(parents=True, exist_ok=True)
json.dump(data, open(OUT, "w"))
print(f"{len(data)-1} counties -> {OUT}\n{data[0]}\n{data[1]}")
