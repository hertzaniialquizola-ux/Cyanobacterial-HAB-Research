"""
Builds data/raw/usgs/well_water_by_county.csv -- the REAL well-water
dependency proxy, replacing the mislabeled Census B25040 ("House Heating
Fuel") variable that was wrongly used for this in every model run before
2026-09-30.

USGS does not offer this as a simple scriptable API pull (ScienceBase's
landing pages block automated fetching), so this is a two-step, mostly-
manual process:

STEP 1 (you, in your own browser -- takes ~2 minutes):
    1. Open: https://doi.org/10.5066/F7TB15V5
       ("Estimated Use of Water in the United States County-Level Data
       for 2015" -- the most recent USGS release with this variable in
       simple county-level tabular form.)
    2. Look for a download link/button for the county-level data file
       (it may be a single .txt/.csv, or a .zip containing one). USGS
       has called this file "usco2015.txt" in past releases, but the
       exact name on the current ScienceBase page may differ -- take
       whatever the "county" data file is called.
    3. Save it into this project at:
           data/raw/usgs/usco2015_raw.txt
       (create the data/raw/usgs/ folder if it doesn't exist -- any
       extension is fine, this script sniffs the delimiter.)

STEP 2 (this script, run by you in Terminal):
       python scripts/get_well_water_usgs.py --diagnose
   This prints every column name found in the raw file, so we can
   confirm the self-supplied-population and total-population columns
   before trusting them (this is the exact check that was skipped last
   time, which is how the B25040 mislabeling happened). Paste that
   output back before running the real conversion.

   Once columns are confirmed (or if the --diagnose output already
   matches the DEFAULT_COLUMN_CANDIDATES guesses below, which it prints
   clearly either way), run:
       python scripts/get_well_water_usgs.py
   to write the clean data/raw/usgs/well_water_by_county.csv that
   build_county_year_panel.py now expects.
"""
import argparse
import csv
import sys
from pathlib import Path

RAW_PATH = Path("data/raw/usgs/usco2015_raw.txt")
OUT_PATH = Path("data/raw/usgs/well_water_by_county.csv")

# Known USGS county water-use naming convention (used consistently across
# the 1995-2015 vintages per USGS's own data dictionaries). We match
# case-insensitively and fuzzily rather than hard-coding one exact name,
# specifically so a naming mismatch is CAUGHT, not silently wrong.
SELF_SUPPLIED_POP_HINTS = ["do-sspop", "dosspop", "do_sspop", "self_supplied_pop", "selfsuppliedpop"]
TOTAL_POP_HINTS = ["tp-totpop", "tptotpop", "tp_totpop", "total_pop", "totpop"]
STATE_FIPS_HINTS = ["statefips", "state_fips", "state_cd", "fipsst"]
COUNTY_FIPS_HINTS = ["countyfips", "county_fips", "fipscty", "fipscounty"]
COMBINED_FIPS_HINTS = ["fips"]  # fallback: a single 5-digit combined FIPS column


def sniff_and_load(path: Path):
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        lines = f.readlines()

    # USGS's usco2015v2.0.csv leads with a one-line citation string before
    # the real header row (confirmed 2026-09-30 against the actual
    # downloaded file). Skip any leading line(s) that don't look like the
    # real header (no FIPS/STATE token) rather than assuming a fixed
    # number of lines to skip, so this still works if a future vintage
    # doesn't have the citation line at all.
    start = 0
    for i, line in enumerate(lines[:5]):
        if "FIPS" in line.upper() or "STATEFIPS" in line.upper():
            start = i
            break

    text = "".join(lines[start:])
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
    except csv.Error:
        dialect = csv.get_dialect("excel-tab") if "\t" in sample else csv.get_dialect("excel")
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    rows = list(reader)
    return rows


def find_col(columns, hints):
    lower = {c.lower().replace(" ", "").replace("-", "").replace("_", ""): c for c in columns}
    for hint in hints:
        key = hint.replace(" ", "").replace("-", "").replace("_", "")
        for lc, orig in lower.items():
            if key in lc:
                return orig
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnose", action="store_true", help="Print columns and exit without writing output.")
    args = ap.parse_args()

    if not RAW_PATH.exists():
        print(f"ERROR: {RAW_PATH} not found.")
        print("See this script's docstring (top of file) for the manual download step.")
        sys.exit(1)

    rows = sniff_and_load(RAW_PATH)
    if not rows:
        print(f"ERROR: {RAW_PATH} parsed to zero rows -- wrong delimiter or empty file.")
        sys.exit(1)

    columns = list(rows[0].keys())
    print(f"Loaded {len(rows)} rows, {len(columns)} columns from {RAW_PATH}")
    print("Columns found:")
    for c in columns:
        print(f"  - {c!r}")

    sspop_col = find_col(columns, SELF_SUPPLIED_POP_HINTS)
    totpop_col = find_col(columns, TOTAL_POP_HINTS)
    state_fips_col = find_col(columns, STATE_FIPS_HINTS)
    county_fips_col = find_col(columns, COUNTY_FIPS_HINTS)
    combined_fips_col = find_col(columns, COMBINED_FIPS_HINTS) if not (state_fips_col and county_fips_col) else None

    print("\nAuto-matched columns (VERIFY these before trusting them):")
    print(f"  self-supplied population -> {sspop_col!r}")
    print(f"  total population         -> {totpop_col!r}")
    print(f"  state FIPS               -> {state_fips_col!r}")
    print(f"  county FIPS              -> {county_fips_col!r}")
    print(f"  combined FIPS (fallback) -> {combined_fips_col!r}")
    print(f"\nFirst data row, for a sanity check:\n  {rows[0]}")

    if args.diagnose:
        print("\n--diagnose mode: not writing output. Confirm the matches above are correct.")
        return

    missing = []
    if not sspop_col:
        missing.append("self-supplied population")
    if not totpop_col:
        missing.append("total population")
    if not (combined_fips_col or (state_fips_col and county_fips_col)):
        missing.append("county FIPS (state+county or combined)")
    if missing:
        print(f"\nERROR: could not confidently auto-match: {missing}.")
        print("Re-run with --diagnose, check the raw column list above, and edit the")
        print("*_HINTS lists at the top of this script to match the real names, or")
        print("tell Claude the exact column names so the script can be corrected.")
        sys.exit(1)

    out_rows = []
    for r in rows:
        if combined_fips_col:
            fips = str(r[combined_fips_col]).strip().zfill(5)
        else:
            fips = str(r[state_fips_col]).strip().zfill(2) + str(r[county_fips_col]).strip().zfill(3)
        try:
            sspop = float(str(r[sspop_col]).replace(",", "").strip())
            totpop = float(str(r[totpop_col]).replace(",", "").strip())
        except ValueError:
            continue
        out_rows.append({"county_fips": fips, "self_supplied_population": sspop, "total_population": totpop})

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["county_fips", "self_supplied_population", "total_population"])
        w.writeheader()
        w.writerows(out_rows)

    print(f"\nWrote {len(out_rows)} county rows to {OUT_PATH}")


if __name__ == "__main__":
    main()
