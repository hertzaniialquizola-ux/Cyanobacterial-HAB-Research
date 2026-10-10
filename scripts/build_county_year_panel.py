"""
build_county_year_panel.py
Cyanobacterial HAB Research — join pipeline

Joins the core data pulls into one county-year panel ready to hand to
R/INLA for the BYM2 model:

    NHDPlus lake polygons  -->  CyAN CI_cyano raster extraction
        --> aggregate to county  --> merge with CDC WONDER mortality
        --> merge with Census well-water proxy
        --> (optional) merge with MODIS/Landsat LST

This script now runs against EITHER scale of data, automatically:
  - Smoke test: data/raw/nhdplus/lake_erie_toledo_smoketest.gpkg (Toledo
    lakes only) + whatever CyAN weeks happen to be in data/raw/cyan/.
  - Full scale: data/raw/nhdplus/lake_national.gpkg (all US lakes, built by
    get_nhdplus_national.R + combine_nhdplus_states.R) + the full
    2018-2024 bloom-season CyAN pull from get_cyan_data_full.py.
The national lake file is preferred when present; otherwise this falls back
to the smoke-test file with a warning. Output filename reflects which one
ran (county_year_panel_smoketest.csv vs county_year_panel_full.csv), so a
full run never silently overwrites the smoke-test result.

2026-09-XX correctness fix before the first full-scale run: the original
version of load_cyan_layer() picked ONE file per (date, composite) grid key,
which was harmless at smoke-test scale (only 2 tiles + a merged mosaic, all
covering the same tiny area) but would have silently thrown away most of the
country's data at national scale, where a given week has many tiles
covering DIFFERENT regions and often no single merged mosaic at all.
load_cyan_layer() now processes every tile for a given (date, composite),
using a merged mosaic (when one exists) as the first choice per lake and
falling back to tiles only for lakes it didn't cover -- so nothing is
double-counted and nothing is silently dropped. It also pre-filters lakes
against each raster's bounding box with a spatial index before attempting
per-lake extraction, which matters once there are thousands of lakes and
hundreds of raster files instead of a few hundred lakes and 6 files.

Prerequisites (run in the `hab-sts` conda env, NOT the Cowork sandbox):
    conda activate hab-sts
    python scripts/build_county_year_panel.py

Required data (see docs/ and the project's progress log for how each was
pulled):
    data/raw/census/county_boundaries/cb_2022_us_county_20m.shp
    data/raw/nhdplus/lake_erie_toledo_smoketest.gpkg (smoke test) and/or
        data/raw/nhdplus/lake_national.gpkg (full scale)
    data/raw/cyan/*.tif
    data/raw/cdc-wonder/Multiple Cause of Death, 2018-2024, Single Race.csv
    data/raw/census/well_water_by_county_2022.json
    data/raw/modis-lst/*.tif (optional -- not pulled yet as of 2026-09-07)
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask as rio_mask
from shapely.geometry import box as shapely_box

# ---------------------------------------------------------------------------
# Paths (relative to project root — run this script from the project root)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
NHDPLUS_NATIONAL_GPKG = ROOT / "data/raw/nhdplus/lake_national.gpkg"
NHDPLUS_SMOKETEST_GPKG = ROOT / "data/raw/nhdplus/lake_erie_toledo_smoketest.gpkg"
CYAN_DIR = ROOT / "data/raw/cyan"
COUNTY_SHP = ROOT / "data/raw/census/county_boundaries/cb_2022_us_county_20m.shp"
CDC_CSV = ROOT / "data/raw/cdc-wonder/Multiple Cause of Death, 2018-2024, Single Race.csv"
WELLWATER_JSON = ROOT / "data/raw/census/well_water_by_county_2022.json"
POVERTY_JSON = ROOT / "data/raw/census/poverty_income_by_county_2022.json"
MODIS_LST_DIR = ROOT / "data/raw/modis-lst"

# WHO (2003) recreational-water "low probability of acute health effects"
# guideline value for cyanobacteria, chlorophyll-a tier: 10 ug/L
# (paired cell-density tier: 20,000 cells/mL). Source: WHO 1999/2003
# Guidelines for Safe Recreational Water Environments Vol. 1, as adopted
# by EPA's 2019 cyanotoxin recreational criteria recommendations.
# CORRECTED 2026-09-30: this constant was previously set to 12.0 with no
# citation, flagged in-code as "provisional pending a cited sensitivity
# check" that was never done. Replaced with the actual cited WHO value.
WHO_THRESHOLD = 10.0  # ug/L chlorophyll-a equivalent (WHO low-risk tier)


# ---------------------------------------------------------------------------
# 1. County boundaries
# ---------------------------------------------------------------------------
def load_counties() -> gpd.GeoDataFrame:
    if not COUNTY_SHP.exists():
        raise FileNotFoundError(
            f"Missing county boundary file: {COUNTY_SHP}\n"
            "See docs/ for the one-time curl command to fetch it (run it on "
            "your own machine, not in a sandboxed shell)."
        )
    counties = gpd.read_file(COUNTY_SHP)[["GEOID", "NAME", "STATEFP", "geometry"]]
    counties = counties.rename(columns={"GEOID": "county_fips"})
    return counties


# ---------------------------------------------------------------------------
# 2. NHDPlus lakes -> assign each lake to a county
# ---------------------------------------------------------------------------
# CONUS Albers Equal Area -- centroid + point-in-polygon math needs a
# projected (metric) CRS, not the geographic lat/lon CRS NHDPlus and the
# Census cartographic boundary file both ship in.
ASSIGNMENT_CRS = "EPSG:5070"


def resolve_nhdplus_path() -> tuple[Path, bool]:
    """Returns (path, is_full_scale). Prefers the national file when present."""
    if NHDPLUS_NATIONAL_GPKG.exists():
        return NHDPLUS_NATIONAL_GPKG, True
    if NHDPLUS_SMOKETEST_GPKG.exists():
        warnings.warn(
            f"{NHDPLUS_NATIONAL_GPKG.name} not found -- falling back to the "
            "Toledo smoke-test lake file. Run get_nhdplus_national.R then "
            "combine_nhdplus_states.R (both in scripts/) for the full-scale run."
        )
        return NHDPLUS_SMOKETEST_GPKG, False
    raise FileNotFoundError(
        "No NHDPlus lake file found -- run get_nhdplus_national.R (full scale) "
        "or the smoke-test AOI query documented in the project's progress log."
    )


def load_lakes_with_county(counties: gpd.GeoDataFrame, nhdplus_path: Path) -> gpd.GeoDataFrame:
    lakes = gpd.read_file(nhdplus_path)

    # Size/type filter (added 2026-09-09): NHDPlus's waterbody layer includes
    # every mapped pond and marsh, most far too small for CyAN's 300m-pixel
    # CI_cyano product to ever resolve a bloom in (median LakePond feature is
    # ~0.02 sq km). Restricting to LakePond/Reservoir >= 1.0 sq km matches
    # "Definition B" in the suppression-map writeup (progress log section 11),
    # which found the unfiltered version carries no suppression signal at all.
    BIG_LAKE_MIN_AREASQKM = 1.0
    BIG_LAKE_FTYPES = {"LakePond", "Reservoir"}
    n_before = len(lakes)
    lakes = lakes[lakes["ftype"].isin(BIG_LAKE_FTYPES) & (lakes["areasqkm"] >= BIG_LAKE_MIN_AREASQKM)]
    print(f"Lake size/type filter: kept {len(lakes)} of {n_before} NHDPlus waterbodies "
          f"(ftype in {BIG_LAKE_FTYPES}, areasqkm >= {BIG_LAKE_MIN_AREASQKM})")

    # Simplification worth flagging to a mentor (Waller/Wakefield, small-area
    # estimation): a lake straddling a county line is assigned wholly to the
    # county containing its centroid, not split/area-weighted across both.
    # Revisit if county-straddling lakes turn out to matter at full scale.
    lakes_proj = lakes.to_crs(ASSIGNMENT_CRS)
    counties_proj = counties.to_crs(ASSIGNMENT_CRS)

    lake_centroids = lakes_proj.copy()
    lake_centroids["geometry"] = lakes_proj.geometry.centroid
    joined = gpd.sjoin(
        lake_centroids, counties_proj[["county_fips", "NAME", "geometry"]],
        how="left", predicate="within",
    ).drop(columns="geometry")

    lakes = lakes.merge(joined[["county_fips", "NAME"]], left_index=True, right_index=True)
    n_unmatched = lakes["county_fips"].isna().sum()
    if n_unmatched:
        warnings.warn(f"{n_unmatched} lake(s) did not fall inside any county polygon "
                       "(likely open water / coastal edge falling outside every county "
                       "polygon) — dropped.")
    return lakes.dropna(subset=["county_fips"])


# ---------------------------------------------------------------------------
# 3. CyAN CI_cyano rasters -> per-lake mean, per date
# ---------------------------------------------------------------------------
_CYAN_FNAME_RE = re.compile(
    r"L(?P<start>\d{7})(?:(?P<end>\d{7}))?\.L3m_(?P<comp>DAY|7D)_CYAN_CI_cyano_CYAN_CONUS_300m(?:_\d+_\d+)?\.tif$"
)
_TILE_SUFFIX_RE = re.compile(r"_\d+_\d+\.tif$")


def _parse_cyan_date(fname: str) -> tuple[str, str]:
    """Returns (period_label, composite_type) e.g. ('2026-day213', 'DAY')."""
    m = _CYAN_FNAME_RE.search(fname)
    if not m:
        raise ValueError(f"Unrecognized CyAN filename pattern: {fname}")
    comp = m.group("comp")
    doy7 = m.group("start")
    year, doy = doy7[:4], doy7[4:]
    label = f"{year}-day{doy}"
    return label, comp


def load_cyan_layer(lakes: gpd.GeoDataFrame) -> pd.DataFrame:
    """
    Extracts mean CI_cyano within each lake polygon for every CyAN raster
    found, grouped by (date, composite).

    A given (date, composite) may have a single merged/non-tiled mosaic, or
    it may have several tiled files each covering a different part of the
    country (the common case at national scale, where no single mosaic
    exists) -- both cases are handled: the merged file (when present) is
    tried first for every lake, then any tiles are tried only for lakes not
    already matched, so a lake is never extracted twice and no region's
    tiles are silently skipped.

    Each raster is also spatially pre-filtered against the lakes with a
    bounding-box index query before attempting per-lake masking -- with
    thousands of lakes and hundreds of raster files at national scale,
    trying every lake against every file would be very slow for almost no
    payoff (most lake/file pairs don't overlap at all).
    """
    # CONUS only, matching the CONUS-only scope decision: the full-scale
    # pull also returns a "CYANAK" (Alaska) product variant with a different
    # filename shape (..._CYANAK_CI_cyano_CYAN_AK_300m...) that
    # _CYAN_FNAME_RE doesn't match. get_cyan_data_full.py now excludes these
    # at download time, but ~394 already landed in data/raw/cyan from the
    # first full-scale run -- filter them out here too rather than relying
    # solely on the downloader, so this script is safe even if AK files
    # (or anything else unexpected) end up in the directory.
    all_tifs = sorted(CYAN_DIR.glob("*_CYAN_CONUS_300m*.tif"))
    if not all_tifs:
        raise FileNotFoundError(f"No CyAN rasters found in {CYAN_DIR}")

    groups: dict[tuple[str, str], list[Path]] = {}
    n_skipped_unparsed = 0
    for tif in all_tifs:
        try:
            label, comp = _parse_cyan_date(tif.name)
        except ValueError as e:
            warnings.warn(f"{e} -- skipping this file.")
            n_skipped_unparsed += 1
            continue
        groups.setdefault((label, comp), []).append(tif)
    if n_skipped_unparsed:
        print(f"  skipped {n_skipped_unparsed} CyAN file(s) with an unrecognized name")

    records = []
    n_groups = len(groups)
    for i, ((label, comp), files) in enumerate(groups.items(), start=1):
        if i == 1 or i % 25 == 0 or i == n_groups:
            print(f"  extracting CyAN group {i}/{n_groups} ({label}, {comp})...")

        non_tiles = [f for f in files if not _TILE_SUFFIX_RE.search(f.name)]
        tiles = [f for f in files if f not in non_tiles]
        ordered_sources = non_tiles + tiles  # merged mosaic first, tiles as fallback

        already_matched: set = set()  # lake index labels already extracted this period
        for tif_path in ordered_sources:
            with rasterio.open(tif_path) as src:
                lakes_r = lakes.to_crs(src.crs) if lakes.crs != src.crs else lakes
                raster_extent = shapely_box(*src.bounds)
                candidate_positions = list(
                    lakes_r.sindex.query(raster_extent, predicate="intersects")
                )
                for pos in candidate_positions:
                    lake_idx = lakes_r.index[pos]
                    if lake_idx in already_matched:
                        continue
                    geom = lakes_r.geometry.iloc[pos]
                    try:
                        out_img, _ = rio_mask(src, [geom], crop=True, nodata=src.nodata)
                    except ValueError:
                        continue  # lake polygon doesn't overlap this raster after all
                    band = out_img[0].astype("float64")
                    nodataval = src.nodata if src.nodata is not None else np.nan
                    valid = band[band != nodataval]
                    valid = valid[~np.isnan(valid)]
                    if valid.size == 0:
                        continue
                    records.append({
                        "county_fips": lakes_r["county_fips"].iloc[pos],
                        "period": label,
                        "composite": comp,
                        "ci_cyano_mean": float(valid.mean()),
                        "ci_cyano_max": float(valid.max()),
                        "n_pixels": int(valid.size),
                    })
                    already_matched.add(lake_idx)

    cyan_df = pd.DataFrame.from_records(records)
    if cyan_df.empty:
        raise RuntimeError("CyAN extraction produced zero records — check CRS/extent overlap.")

    # county-level: pixel-count-weighted mean across lakes in a county
    county_cyan = (
        cyan_df.groupby(["county_fips", "period", "composite"])
        .apply(lambda g: np.average(g["ci_cyano_mean"], weights=g["n_pixels"]), include_groups=False)
        .rename("ci_cyano_mean")
        .reset_index()
    )
    county_cyan["year"] = county_cyan["period"].str.slice(0, 4).astype(int)
    return county_cyan


def build_cei(county_cyan: pd.DataFrame) -> pd.DataFrame:
    """
    Cumulative Exposure Intensity per county-year: fraction of observed
    periods that year where mean CI_cyano exceeded the WHO threshold (a
    rate, not a raw count, so it's comparable across counties/years with
    different numbers of observed periods).
    """
    df = county_cyan.copy()
    df["exceeds"] = df["ci_cyano_mean"] > WHO_THRESHOLD
    cei = (
        df.groupby(["county_fips", "year"])
        .agg(cei=("exceeds", "mean"), n_periods_observed=("exceeds", "size"))
        .reset_index()
    )
    return cei


# ---------------------------------------------------------------------------
# 4. CDC WONDER mortality (handles suppression + hidden zero-death rows)
# ---------------------------------------------------------------------------
def load_mortality() -> pd.DataFrame:
    raw = pd.read_csv(CDC_CSV, dtype=str)
    # data rows have a 5-digit FIPS in "County Code"; footer/notes rows don't
    raw = raw[raw["County Code"].str.match(r"^\d{5}$", na=False)].copy()

    raw["year"] = raw["Year"].str.strip().astype(int)
    raw["county_fips"] = raw["County Code"]
    raw["deaths_suppressed"] = raw["Deaths"].eq("Suppressed")
    raw["deaths"] = pd.to_numeric(raw["Deaths"], errors="coerce")  # NaN where suppressed
    raw["population"] = pd.to_numeric(raw["Population"], errors="coerce")

    mort = raw[["county_fips", "year", "deaths", "population", "deaths_suppressed"]]

    # WONDER was exported with "Show Zero Values: False" -- county-years
    # with zero N18/K74/K76 deaths simply don't appear as rows at all. Any
    # (county, year) combo present in `population` data but absent here
    # is a true zero, NOT a missing value -- reconstruct the full grid so
    # the model doesn't accidentally drop or NaN-out real zero counties.
    all_fips = mort["county_fips"].unique()
    all_years = range(2018, 2025)
    full_grid = pd.MultiIndex.from_product([all_fips, all_years], names=["county_fips", "year"]).to_frame(index=False)
    mort_full = full_grid.merge(mort, on=["county_fips", "year"], how="left")
    mort_full["deaths_suppressed"] = mort_full["deaths_suppressed"].fillna(False).infer_objects(copy=False).astype(bool)
    mort_full.loc[mort_full["deaths"].isna() & ~mort_full["deaths_suppressed"], "deaths"] = 0
    return mort_full


# ---------------------------------------------------------------------------
# 5. Well-water dependency proxy (USGS domestic self-supplied population)
# ---------------------------------------------------------------------------
# CORRECTED 2026-09-30: the previous version of this function read
# B25040_010E/B25040_001E from Census ACS table B25040. That table is
# "House Heating Fuel" -- B25040_010E is households reporting "No fuel
# used" for HEATING, not anything about water source. It had been
# mislabeled as a well-water variable since the original data pull and
# was silently wrong in every model run to date. See
# claude/hab-well-water-covariate-fix-2026-09-30.md for the full writeup.
#
# Replaced with USGS's actual domestic self-supplied population dataset
# (the standard public proxy for private-well dependency), built by
# scripts/get_well_water_usgs.py into WELLWATER_USGS_CSV. That script
# prints the raw source file's column headers before mapping them, and
# this loader re-validates the mapped columns at runtime -- on purpose,
# so a future column-name mismatch fails loudly instead of silently
# producing a plausible-looking wrong number again.
WELLWATER_USGS_CSV = ROOT / "data/raw/usgs/well_water_by_county.csv"


def load_wellwater() -> pd.DataFrame:
    if not WELLWATER_USGS_CSV.exists():
        raise FileNotFoundError(
            f"Missing {WELLWATER_USGS_CSV}. Run scripts/get_well_water_usgs.py "
            "first (see that script's docstring for the USGS source file "
            "you need to download manually first)."
        )
    df = pd.read_csv(WELLWATER_USGS_CSV, dtype={"county_fips": str})
    required = {"county_fips", "self_supplied_population", "total_population"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{WELLWATER_USGS_CSV.name} is missing expected column(s) {missing}. "
            f"Actual columns: {list(df.columns)}. Re-check "
            "scripts/get_well_water_usgs.py's column mapping against the raw "
            "USGS file rather than assuming the old names still apply."
        )
    df["county_fips"] = df["county_fips"].str.zfill(5)
    df["well_water_fraction"] = (
        df["self_supplied_population"] / df["total_population"]
    )
    return df[["county_fips", "well_water_fraction"]]


def load_poverty() -> pd.DataFrame | None:
    """SES covariates for the suppression multiple-imputation model
    (scripts/fit_bym2_model_mi.R). Run scripts/get_census_poverty.py first
    if POVERTY_JSON doesn't exist yet -- this returns None and the join is
    skipped (with a warning) rather than failing, so the rest of the panel
    still builds without it."""
    if not POVERTY_JSON.exists():
        warnings.warn(
            f"{POVERTY_JSON.name} not found -- poverty_rate/median_income "
            "will be absent from the panel. Run scripts/get_census_poverty.py "
            "first if you need them (e.g. for the suppression MI model)."
        )
        return None
    with open(POVERTY_JSON) as f:
        rows = json.load(f)
    header, data = rows[0], rows[1:]
    df = pd.DataFrame(data, columns=header)
    df["county_fips"] = df["state"] + df["county"]
    below_poverty = pd.to_numeric(df["B17001_002E"], errors="coerce")
    poverty_denom = pd.to_numeric(df["B17001_001E"], errors="coerce")
    df["poverty_rate"] = below_poverty / poverty_denom
    df["median_household_income"] = pd.to_numeric(df["B19013_001E"], errors="coerce")
    return df[["county_fips", "poverty_rate", "median_household_income"]]


# ---------------------------------------------------------------------------
# 6. MODIS LST (MOD11A2, 8-day 1km) -- pull script written 2026-09-19
#    (scripts/get_modis_lst.py); this loader is a FIRST DRAFT, UNTESTED
#    against a real downloaded file as of 2026-09-19. Run the smoke test in
#    get_modis_lst.py's docstring first and confirm the subdataset name /
#    HDF4 driver work before trusting this at scale.
# ---------------------------------------------------------------------------
_MODIS_LST_FNAME_RE = re.compile(
    r"MOD11A2\.A(?P<year>\d{4})(?P<doy>\d{3})\.h\d{2}v\d{2}\.061\.\d+\.hdf$"
)
_MODIS_LST_SCALE = 0.02   # LST_Day_1km / LST_Night_1km scale factor -> Kelvin
_MODIS_LST_FILL = 0       # fill value for both subdatasets
_MODIS_SUBDATASETS = {
    "lst_day_c": "LST_Day_1km",
    "lst_night_c": "LST_Night_1km",
}


def _parse_modis_lst_date(fname: str) -> str:
    """Returns a period label like '2020-day213' -- same label shape as
    _parse_cyan_date's output, for eyeballing/debugging consistency (the
    two sources are NOT joined on this label directly; both ultimately
    roll up to county-year before merging into the panel)."""
    m = _MODIS_LST_FNAME_RE.search(fname)
    if not m:
        raise ValueError(f"Unrecognized MODIS LST filename pattern: {fname}")
    return f"{m.group('year')}-day{m.group('doy')}"


def load_lst(lakes: gpd.GeoDataFrame) -> pd.DataFrame | None:
    """
    Per-lake MOD11A2 LST extraction, mirroring load_cyan_layer()'s pattern:
    spatial-index-prefiltered mask extraction per raster, per-lake mean,
    aggregated to county-year. Two things differ from CyAN and are the
    parts most worth double-checking against a real file before trusting
    this at scale (see the module docstring in get_modis_lst.py for the
    smoke-test snippet):
      1. MODIS ships as HDF-EOS sinusoidal tiles, not a plain GeoTIFF --
         opened here via GDAL's `HDF4_EOS:EOS_GRID:"<file>":<grid>:<var>`
         connection-string syntax, which requires GDAL to have been built
         with HDF4 support.
      2. LST_Day_1km / LST_Night_1km are stored as raw uint16 with a 0.02
         scale factor (-> Kelvin) and fill value 0 -- converted to Celsius
         here (raw * 0.02 - 273.15).
    """
    if not MODIS_LST_DIR.exists() or not any(MODIS_LST_DIR.iterdir()):
        warnings.warn(
            "data/raw/modis-lst is empty -- run scripts/get_modis_lst.py "
            "(smoke test first, per its docstring) before this will do "
            "anything. Continuing the join WITHOUT LST."
        )
        return None

    all_hdf = sorted(MODIS_LST_DIR.glob("MOD11A2.*.hdf"))
    if not all_hdf:
        warnings.warn(
            f"{MODIS_LST_DIR} has files but none match the expected "
            "MOD11A2.A<year><doy>.h##v##.061.*.hdf pattern -- skipping LST."
        )
        return None

    records = []
    n_skipped_unparsed = 0
    n_skipped_open_error = 0
    for i, hdf_path in enumerate(all_hdf, start=1):
        if i == 1 or i % 50 == 0 or i == len(all_hdf):
            print(f"  extracting MODIS LST file {i}/{len(all_hdf)}...")
        try:
            label = _parse_modis_lst_date(hdf_path.name)
        except ValueError as e:
            warnings.warn(f"{e} -- skipping.")
            n_skipped_unparsed += 1
            continue

        for var_name, subdataset in _MODIS_SUBDATASETS.items():
            gdal_path = f'HDF4_EOS:EOS_GRID:"{hdf_path}":MODIS_Grid_8Day_1km_LST:{subdataset}'
            try:
                with rasterio.open(gdal_path) as src:
                    lakes_r = lakes.to_crs(src.crs) if lakes.crs != src.crs else lakes
                    raster_extent = shapely_box(*src.bounds)
                    candidate_positions = list(
                        lakes_r.sindex.query(raster_extent, predicate="intersects")
                    )
                    for pos in candidate_positions:
                        geom = lakes_r.geometry.iloc[pos]
                        try:
                            out_img, _ = rio_mask(src, [geom], crop=True, nodata=_MODIS_LST_FILL)
                        except ValueError:
                            continue
                        band = out_img[0].astype("float64")
                        valid = band[band != _MODIS_LST_FILL]
                        if valid.size == 0:
                            continue
                        celsius = valid * _MODIS_LST_SCALE - 273.15
                        records.append({
                            "county_fips": lakes_r["county_fips"].iloc[pos],
                            "period": label,
                            "variable": var_name,
                            "value_c": float(celsius.mean()),
                            "n_pixels": int(valid.size),
                        })
            except Exception as e:
                n_skipped_open_error += 1
                if n_skipped_open_error <= 3:
                    warnings.warn(
                        f"Could not open {subdataset} in {hdf_path.name}: {e} -- "
                        "if this happens for every file, GDAL's HDF4 driver is "
                        "probably missing or the subdataset name is wrong; run "
                        "the smoke-test snippet in get_modis_lst.py's docstring "
                        "directly to debug before re-running this."
                    )
                continue

    if n_skipped_unparsed:
        print(f"  skipped {n_skipped_unparsed} file(s) with an unrecognized name")
    if n_skipped_open_error:
        print(f"  {n_skipped_open_error} file/subdataset open error(s) -- see warnings above")

    if not records:
        warnings.warn("MODIS LST extraction produced zero records -- continuing WITHOUT LST.")
        return None

    lst_df = pd.DataFrame.from_records(records)
    lst_df["year"] = lst_df["period"].str.slice(0, 4).astype(int)

    wide = (
        lst_df.groupby(["county_fips", "year", "variable"])
        .apply(lambda g: np.average(g["value_c"], weights=g["n_pixels"]), include_groups=False)
        .rename("value_c")
        .reset_index()
        .pivot(index=["county_fips", "year"], columns="variable", values="value_c")
        .reset_index()
    )
    wide.columns.name = None
    return wide


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    nhdplus_path, is_full_scale = resolve_nhdplus_path()
    out_csv = ROOT / (
        "outputs/county_year_panel_full.csv" if is_full_scale
        else "outputs/county_year_panel_smoketest.csv"
    )
    print(f"NHDPlus source: {nhdplus_path.name} ({'full scale' if is_full_scale else 'smoke test'})")

    counties = load_counties()
    lakes = load_lakes_with_county(counties, nhdplus_path)
    print(f"{len(lakes)} lakes assigned to counties; extracting CyAN data from {CYAN_DIR}...")
    county_cyan = load_cyan_layer(lakes)
    cei = build_cei(county_cyan)
    mortality = load_mortality()
    wellwater = load_wellwater()
    poverty = load_poverty()
    lst = load_lst(lakes)

    panel = mortality.merge(cei, on=["county_fips", "year"], how="left")
    panel = panel.merge(wellwater, on="county_fips", how="left")
    if poverty is not None:
        panel = panel.merge(poverty, on="county_fips", how="left")
    if lst is not None:
        panel = panel.merge(lst, on=["county_fips", "year"], how="left")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(out_csv, index=False)

    n_cei_matched = panel["cei"].notna().sum()
    n_counties_matched = panel.loc[panel["cei"].notna(), "county_fips"].nunique()

    print(f"Wrote {len(panel)} county-year rows to {out_csv}")
    print(f"  county-year rows with CyAN exposure data: {n_cei_matched}")
    print(f"  distinct counties with any CyAN exposure data: {n_counties_matched}")
    print(f"  suppressed mortality rows: {panel['deaths_suppressed'].sum()}")
    print(f"  MODIS LST merged: {lst is not None}")
    print(f"  Poverty/income covariates merged: {poverty is not None}")

    if n_cei_matched == 0:
        cyan_years = sorted(county_cyan["year"].unique())
        print(
            f"\n  NOTE: zero rows matched CyAN exposure to mortality by (county, year).\n"
            f"  CyAN data covers year(s) {cyan_years}; CDC WONDER mortality only covers "
            f"2018-2024.\n"
            f"  If {cyan_years} falls outside 2018-2024, that is the reason -- not a join bug."
        )


if __name__ == "__main__":
    main()
