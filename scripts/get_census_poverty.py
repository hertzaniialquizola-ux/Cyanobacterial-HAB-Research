"""
get_census_poverty.py

Pulls county-level poverty rate and median household income from the ACS
5-year estimates -- the SES covariates needed for the suppression multiple-
imputation model (scripts/fit_bym2_model_mi.R). Same API, same pattern as
the well-water pull (data/raw/census/well_water_by_county_2022.json), just
a different table.

Variables:
  B17001_002E  -- population with income in the past 12 months below poverty level
  B17001_001E  -- population for whom poverty status is determined (denominator)
  B19013_001E  -- median household income (inflation-adjusted dollars)

Usage:
  export CENSUS_API_KEY="<key>"   # same key used for the well-water pull
  python3 scripts/get_census_poverty.py

Writes: data/raw/census/poverty_income_by_county_2022.json
"""
import json
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = ROOT / "data/raw/census/poverty_income_by_county_2022.json"

API_KEY = os.environ.get("CENSUS_API_KEY")
if not API_KEY:
    sys.exit(
        "CENSUS_API_KEY is not set. Export the same key used for the "
        "well-water pull:\n  export CENSUS_API_KEY=\"<key>\"\n"
        "Rotate it at https://api.census.gov/data/key_signup.html if it's "
        "ever been shared or pasted in plaintext."
    )

URL = (
    "https://api.census.gov/data/2022/acs/acs5"
    "?get=NAME,B17001_002E,B17001_001E,B19013_001E"
    "&for=county:*"
    f"&key={API_KEY}"
)

print("Requesting ACS poverty + median household income by county...")
try:
    with urllib.request.urlopen(URL, timeout=60) as resp:
        data = json.loads(resp.read().decode())
except urllib.error.HTTPError as e:
    sys.exit(f"Census API request failed: HTTP {e.code} -- {e.read().decode()[:500]}")
except urllib.error.URLError as e:
    sys.exit(f"Census API request failed: {e.reason}")

n_rows = len(data) - 1  # first row is the header
print(f"Received {n_rows} county rows.")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(data, f)

print(f"Wrote {OUT_PATH}")
print("Schema:", data[0])
print("Example row:", data[1])
