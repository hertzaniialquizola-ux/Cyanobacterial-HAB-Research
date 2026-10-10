"""
get_modis_lst.py
Full-scale MODIS Land Surface Temperature pull -- the one core data source
in the research design (project doc section 5) that has never been started.

Product: MOD11A2 v061 (Terra, 8-Day, 1km, Land Surface Temperature/Emissivity,
L3 Global, MODIS Sinusoidal Grid), served through NASA's LP DAAC / LPCLOUD
provider -- same earthaccess-based pull mechanism already proven out for
CyAN (get_cyan_data_full.py), collection resolved dynamically by
short_name+version below rather than a hardcoded concept_id (a concept_id
could not be independently verified from this end; short_name/version is
the documented, stable way to ask earthaccess/CMR for "give me this product,
whatever its current concept_id is").

Why MOD11A2 (Terra) and not also MYD11A2 (Aqua): starting with one satellite
keeps the first real run simple and matches the CyAN precedent (smoke test
narrow, then decide whether to expand). Add MYD11A2 later as a second
SHORT_NAME/VERSION pair in the same loop if daily-scale (Terra+Aqua
combined) coverage turns out to matter for the lake-thermal-regime signal.

Scope decisions baked in below (same pattern as get_cyan_data_full.py --
edit the constants if you want something different):
  - Bloom season only, May 1 - Oct 31, matching CyAN's window (this is a
    lake-thermal-regime covariate for the exposure model, so it only needs
    to cover the same season CEI is computed over).
  - Full CONUS bounding box (no Alaska/Hawaii/territories -- matches CyAN
    and the rest of the pipeline's CONUS-only scope).
  - 2018-2024, matching CDC WONDER's mortality years.
  - 8-day composite only (that's what MOD11A2 is -- there is no "daily"
    variant of this exact short_name; MOD11A1 is the daily product, out of
    scope here).

IMPORTANT -- read before running at full scale:
  1. RUN THE SMOKE TEST FIRST (see bottom of this docstring). MODIS
     granules download as HDF-EOS (.hdf) files, not the same GeoTIFF format
     CyAN used -- county_year_panel_full.csv's load_lst() extraction step
     (build_county_year_panel.py) needs GDAL's HDF4 driver to read the
     LST_Day_1km/LST_Night_1km subdatasets, which is common but not
     guaranteed to be present in every GDAL build. Confirm one file opens
     and extracts correctly (see the smoke-test snippet below) BEFORE
     kicking off a multi-hour, many-GB full pull that might turn out to be
     unusable.
  2. earthaccess must be installed and logged in the same way as the CyAN
     pull (`pip install earthaccess` if needed; `earthaccess.login()` below
     uses the same .netrc/credentials as get_cyan_data_full.py already
     relies on).
  3. This is a MODIS SINUSOIDAL-GRID product (h/v tile numbers, not a
     CONUS-mosaic file like CyAN) -- CONUS is covered by roughly a dozen
     tiles (h07-h13, v04-v06), all of which download as part of a normal
     bounding-box search; there is no single "CONUS mosaic" file to prefer
     the way load_cyan_layer() does. The loader (build_county_year_panel.py)
     handles this by extracting per-lake means from whichever tile(s)
     spatially overlap that lake, same sindex-prefilter idea as CyAN.

Safe to interrupt and re-run: each period's granules are checked against
what's already on disk before downloading (same already_on_disk() pattern
as get_cyan_data_full.py), so a killed/crashed run picks back up close to
where it left off.

SMOKE TEST before the full run -- confirm the mechanism works on ONE
period/tile before trusting the full loop:
    python3 -c "
import earthaccess
earthaccess.login()
r = earthaccess.search_data(short_name='MOD11A2', version='061',
    temporal=('2020-08-01','2020-08-07'),
    bounding_box=(-83.7, 41.5, -82.6, 42.3))
print(len(r), 'granule(s) found')
files = earthaccess.download(r, 'data/raw/modis-lst')
print(files)
"
Then confirm the downloaded .hdf file actually opens with rasterio/GDAL:
    python3 -c "
import rasterio
path = 'HDF4_EOS:EOS_GRID:\"data/raw/modis-lst/<the .hdf filename>\":MODIS_Grid_8Day_1km_LST:LST_Day_1km'
with rasterio.open(path) as src:
    print(src.crs, src.bounds, src.read(1).shape)
"
If either of those fails (module not found, HDF4 driver missing, wrong
subdataset name), STOP and fix that before running main() at scale -- the
subdataset name in particular is exactly the kind of thing worth checking
against the real file rather than assuming (same lesson as the NHDPlus
`get_waterbodies()` API-guessing cost in section 1 of the project doc).
"""

import datetime as dt
import threading
import time
from pathlib import Path

import earthaccess

SHORT_NAME = "MOD11A2"
VERSION = "061"
OUT_DIR = Path("data/raw/modis-lst")
YEARS = range(2018, 2025)          # 2018-2024 inclusive, matches CDC WONDER
SEASON_START = (5, 1)              # May 1
SEASON_END = (10, 31)              # Oct 31
PERIOD_STEP_DAYS = 8                # MOD11A2 is an 8-day composite
# CONUS bounding box (west, south, east, north) -- same idea as CyAN's
# full-CONUS scope, expressed as a literal bbox here since earthaccess's
# bounding_box search takes one directly (CyAN's collection-level search
# had no bbox filter at all; this product does, and CONUS-only keeps the
# tile count -- and download volume -- down to what the pipeline needs).
CONUS_BBOX = (-125.0, 24.5, -66.9, 49.4)

# Same rationale as get_cyan_data_full.py: a single period's worth of tiles
# is small, so if a call hasn't returned in this long, something is stuck,
# not just slow.
SEARCH_TIMEOUT_S = 120
DOWNLOAD_TIMEOUT_S = 600


def _run_with_timeout(fn, timeout_s, *args, **kwargs):
    """Runs fn(*args, **kwargs) on a fresh, one-off DAEMON thread and
    enforces a real wall-clock timeout.

    NOTE on why this is a fresh thread per call, not a shared thread pool:
    get_cyan_data_full.py's first timeout fix reused one shared
    single-worker ThreadPoolExecutor across every call. That looked fine
    but had a real bug (project doc section 4, Round 2): Python cannot
    actually kill a thread stuck in a blocking call, so once one call
    genuinely hung, that single shared worker was occupied forever and
    every subsequent call silently queued behind it -- they weren't each
    hanging independently, they just never got a turn, which showed up as
    a run of back-to-back identical-duration "timeouts." A fresh daemon
    thread per call means a stuck call only strands its own one-off
    thread; every later call still gets to run immediately, and daemon
    threads don't block process exit either.
    """
    result: list = []
    error: list = []

    def target():
        try:
            result.append(fn(*args, **kwargs))
        except Exception as e:  # noqa: BLE001 -- re-raised on the calling thread below
            error.append(e)

    t = threading.Thread(target=target, daemon=True)
    t.start()
    t.join(timeout_s)
    if t.is_alive():
        raise TimeoutError(f"{fn.__name__} did not return within {timeout_s}s")
    if error:
        raise error[0]
    return result[0]


def period_windows_for_year(year: int):
    start = dt.date(year, *SEASON_START)
    end = dt.date(year, *SEASON_END)
    cur = start
    while cur <= end:
        window_end = min(cur + dt.timedelta(days=PERIOD_STEP_DAYS - 1), end)
        yield cur, window_end
        cur = window_end + dt.timedelta(days=1)


def already_on_disk(granule) -> bool:
    fnames = [Path(link).name for link in granule.data_links()]
    return bool(fnames) and all((OUT_DIR / f).exists() for f in fnames)


def main() -> None:
    earthaccess.login()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    total_periods = 0
    total_downloaded_files = 0
    total_skipped_granules = 0
    periods_with_zero_granules = []
    periods_with_download_errors = []
    periods_with_timeouts = []

    for year in YEARS:
        for start, end in period_windows_for_year(year):
            total_periods += 1
            temporal = (start.isoformat(), end.isoformat())
            t0 = time.monotonic()
            print(f"[{year}] {temporal[0]} .. {temporal[1]} ...", end=" ", flush=True)

            try:
                results = _run_with_timeout(
                    earthaccess.search_data, SEARCH_TIMEOUT_S,
                    short_name=SHORT_NAME, version=VERSION,
                    temporal=temporal, bounding_box=CONUS_BBOX,
                )
            except TimeoutError:
                elapsed = time.monotonic() - t0
                print(f"SEARCH TIMED OUT after {elapsed:.0f}s")
                periods_with_timeouts.append((temporal, "search", elapsed))
                continue
            except Exception as e:
                print(f"SEARCH FAILED: {e}")
                continue

            if not results:
                print("0 granules found")
                periods_with_zero_granules.append(temporal)
                continue

            to_fetch = [g for g in results if not already_on_disk(g)]
            total_skipped_granules += len(results) - len(to_fetch)

            if not to_fetch:
                print(f"{len(results)} granule(s) already on disk, skipped "
                      f"({time.monotonic() - t0:.0f}s)")
                continue

            try:
                files = _run_with_timeout(
                    earthaccess.download, DOWNLOAD_TIMEOUT_S,
                    to_fetch, str(OUT_DIR),
                )
                total_downloaded_files += len(files)
                print(f"downloaded {len(files)} file(s) ({time.monotonic() - t0:.0f}s)")
            except TimeoutError:
                elapsed = time.monotonic() - t0
                print(f"DOWNLOAD TIMED OUT after {elapsed:.0f}s")
                periods_with_timeouts.append((temporal, "download", elapsed))
                continue
            except Exception as e:
                periods_with_download_errors.append((temporal, str(e)))
                print(f"DOWNLOAD FAILED: {e}")
                continue

    print()
    print(f"Done. {total_periods} periods processed, {total_downloaded_files} files "
          f"downloaded, {total_skipped_granules} granule(s) already present.")
    if periods_with_zero_granules:
        print(f"{len(periods_with_zero_granules)} period(s) had zero granules found. "
              f"First few: {periods_with_zero_granules[:5]}")
    if periods_with_timeouts:
        print(f"\n{len(periods_with_timeouts)} period(s) hit a hard timeout and were skipped:")
        for temporal, stage, elapsed in periods_with_timeouts[:10]:
            print(f"  {temporal[0]}..{temporal[1]}: {stage} still running after {elapsed:.0f}s")
    if periods_with_download_errors:
        print(f"\n{len(periods_with_download_errors)} period(s) hit a download error:")
        for temporal, err in periods_with_download_errors[:10]:
            print(f"  {temporal[0]}..{temporal[1]}: {err[:150]}")
    if periods_with_timeouts or periods_with_download_errors or periods_with_zero_granules:
        print("\nJust re-run this script again -- already-downloaded files are "
              "skipped, so it will only retry what's still missing.")


if __name__ == "__main__":
    main()
