"""
get_cyan_data_full.py
Full-scale CyAN CI_cyano pull: bloom season, 2018-2024, full continental US,
7-day composite product only.

This is a MUCH larger pull than the original smoke test (get_cyan_data.py,
one week/one small bounding box). Scope decisions baked in below -- edit the
constants if you want something different:

  - Bloom season only: May 1 - Oct 31 each year. Cyanobacteria blooms are
    overwhelmingly a warm-season phenomenon; this cuts total volume roughly
    in half vs a full 12-month pull, with minimal loss of real bloom signal.
  - Full CONUS: no bounding_box restriction, so this pulls every region the
    product covers, not just the Toledo smoke-test area -- EXCEPT Alaska
    (see EXCLUDE_ALASKA below).
  - 7-day composites (L3m_7D) only, not daily (L3m_DAY) -- matches the
    "Temporal Resolution: 7-day composite" the research design specifies,
    and roughly halves volume again vs also pulling daily granules. (The
    smoke test proved the daily product downloads fine too, if a future
    need for it comes up.)

These are defaults I (Claude) picked without checking with you first --
change YEARS / SEASON_START / SEASON_END / DAILY_INSTEAD_OF_7D below if you
want a different scope, then just re-run (already-downloaded files are
skipped, so re-running after a change only fetches what's new).

Safe to interrupt and re-run: each week's granules are checked against what's
already on disk before downloading, so a killed/crashed run just picks back
up close to where it left off next time.

Prerequisites: `pip install earthaccess` in the hab-sts env if not already
done (it's not part of the conda env's original package list -- this bit
the first full-scale attempt on 2026-09-07).

2026-09-07 fixes (after the first full-scale run):
  - The collection also serves a "CYANAK" (Alaska) product variant mixed in
    with the CONUS one -- same "L3m_7D" substring, different filename shape
    (..._CYANAK_CI_cyano_CYAN_AK_300m...). Alaska is out of scope (matches
    no county in the CDC WONDER / Census join, and build_county_year_panel.py's
    filename parser doesn't recognize it), so it's now excluded up front
    instead of being downloaded and later choking the join script. ~394 AK
    files from the first run are sitting in data/raw/cyan already; they're
    harmless to leave (the join script's file glob was also fixed to ignore
    them) but you can delete them to reclaim disk space if you want.
  - search_data()/download() calls are now wrapped with a hard timeout
    (SEARCH_TIMEOUT_S / DOWNLOAD_TIMEOUT_S below) running on a background
    thread. The first full-scale run appeared to silently stop mid-2019 --
    with no exception, no traceback, nothing -- which is the signature of a
    network call blocking forever inside a library with no timeout of its
    own (unlike the ConnectTimeoutErrors seen elsewhere, which DID eventually
    raise and get caught). A plain try/except can't catch a hang that never
    raises; only a timeout that runs the call on a separate thread and gives
    up on waiting for it can. On a timeout the week is now logged and
    skipped, same as any other failure, and the run keeps going.
  - Each week now prints how long it took, so a future stall is visible in
    real time instead of only discoverable after the fact from file mtimes.
"""

import concurrent.futures
import datetime as dt
import time
from pathlib import Path

import earthaccess

CONCEPT_ID = "C3416412382-OB_CLOUD"
OUT_DIR = Path("data/raw/cyan")
YEARS = range(2018, 2025)        # 2018-2024 inclusive
SEASON_START = (5, 1)            # May 1
SEASON_END = (10, 31)            # Oct 31
WEEK_STEP_DAYS = 7
EXCLUDE_ALASKA = True            # CONUS-only scope decision -- see note above

# Hard wall-clock caps on the two network calls. A single lake-week's worth
# of files is small (tens of MB); if either call hasn't returned in this
# long, something is stuck, not just slow -- give up on that week and move
# on rather than blocking the whole multi-hour run indefinitely.
SEARCH_TIMEOUT_S = 120
DOWNLOAD_TIMEOUT_S = 600

_executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)


def _run_with_timeout(fn, timeout_s, *args, **kwargs):
    """Runs fn(*args, **kwargs) on a worker thread and enforces a real
    wall-clock timeout. A plain try/except around a blocking network call
    cannot do this -- if the call never raises (a true hang), the only way
    to move on is to stop waiting on a separate thread. Note the underlying
    call keeps running in the background after a timeout (Python can't kill
    a thread), so a timed-out request may still finish writing files later;
    already_on_disk() will just see it as done on a future run."""
    future = _executor.submit(fn, *args, **kwargs)
    return future.result(timeout=timeout_s)


def week_windows_for_year(year: int):
    start = dt.date(year, *SEASON_START)
    end = dt.date(year, *SEASON_END)
    cur = start
    while cur <= end:
        window_end = min(cur + dt.timedelta(days=WEEK_STEP_DAYS - 1), end)
        yield cur, window_end
        cur = window_end + dt.timedelta(days=1)


def is_alaska_granule(granule) -> bool:
    return any(
        "CYANAK" in link or "_AK_" in link
        for link in granule.data_links()
    )


def already_on_disk(granule) -> bool:
    fnames = [Path(link).name for link in granule.data_links()]
    return bool(fnames) and all((OUT_DIR / f).exists() for f in fnames)


def main() -> None:
    earthaccess.login()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    total_weeks = 0
    total_downloaded_files = 0
    total_skipped_granules = 0
    total_alaska_skipped = 0
    weeks_with_zero_granules = []
    weeks_with_download_errors = []
    weeks_with_timeouts = []

    for year in YEARS:
        for start, end in week_windows_for_year(year):
            total_weeks += 1
            temporal = (start.isoformat(), end.isoformat())
            t0 = time.monotonic()
            print(f"[{year}] {temporal[0]} .. {temporal[1]} ...", end=" ", flush=True)

            try:
                results = _run_with_timeout(
                    earthaccess.search_data, SEARCH_TIMEOUT_S,
                    concept_id=CONCEPT_ID, temporal=temporal,
                )
            except concurrent.futures.TimeoutError:
                elapsed = time.monotonic() - t0
                print(f"SEARCH TIMED OUT after {elapsed:.0f}s")
                weeks_with_timeouts.append((temporal, "search", elapsed))
                continue
            except Exception as e:
                print(f"SEARCH FAILED: {e}")
                continue

            # the collection mixes L3m_DAY and L3m_7D granules together --
            # keep only the 7-day composites per the scope decision above
            seven_day = [
                g for g in results
                if any("L3m_7D" in link for link in g.data_links())
            ]

            if EXCLUDE_ALASKA:
                before = len(seven_day)
                seven_day = [g for g in seven_day if not is_alaska_granule(g)]
                total_alaska_skipped += before - len(seven_day)

            if not seven_day:
                print("0 seven-day granules found")
                weeks_with_zero_granules.append(temporal)
                continue

            to_fetch = [g for g in seven_day if not already_on_disk(g)]
            total_skipped_granules += len(seven_day) - len(to_fetch)

            if not to_fetch:
                print(f"{len(seven_day)} granule(s) already on disk, skipped "
                      f"({time.monotonic() - t0:.0f}s)")
                continue

            # NASA's on-prem GSFC server (which even OB_CLOUD-provider
            # granules seem to route through for actual file bytes, not just
            # an AWS/cloud-native path) is known to be flaky -- the original
            # smoke test hit a full outage there once already. A single
            # connection timeout on ONE file used to kill this entire
            # multi-hour run; now it's caught, logged, and the loop moves on.
            # A hard timeout (not just try/except) also guards against a
            # request that hangs forever instead of raising -- see
            # _run_with_timeout's docstring.
            try:
                files = _run_with_timeout(
                    earthaccess.download, DOWNLOAD_TIMEOUT_S,
                    to_fetch, str(OUT_DIR),
                )
                total_downloaded_files += len(files)
                print(f"downloaded {len(files)} file(s) ({time.monotonic() - t0:.0f}s)")
            except concurrent.futures.TimeoutError:
                elapsed = time.monotonic() - t0
                print(f"DOWNLOAD TIMED OUT after {elapsed:.0f}s")
                weeks_with_timeouts.append((temporal, "download", elapsed))
                continue
            except Exception as e:
                weeks_with_download_errors.append((temporal, str(e)))
                print(f"DOWNLOAD FAILED: {e}")
                continue

    print()
    print(f"Done. {total_weeks} weeks processed, {total_downloaded_files} files "
          f"downloaded, {total_skipped_granules} granule(s) already present, "
          f"{total_alaska_skipped} Alaska granule(s) excluded.")
    if weeks_with_zero_granules:
        print(f"{len(weeks_with_zero_granules)} week(s) had zero 7-day granules "
              f"found (CyAN coverage isn't perfectly uniform week to week -- this "
              f"can be normal). First few: {weeks_with_zero_granules[:5]}")
    if weeks_with_timeouts:
        print(f"\n{len(weeks_with_timeouts)} week(s) hit a hard timeout (call never "
              f"returned) and were skipped:")
        for temporal, stage, elapsed in weeks_with_timeouts[:10]:
            print(f"  {temporal[0]}..{temporal[1]}: {stage} still running after {elapsed:.0f}s")
    if weeks_with_download_errors:
        print(f"\n{len(weeks_with_download_errors)} week(s) hit a download error "
              f"(often a transient NASA server timeout) and were skipped:")
        for temporal, err in weeks_with_download_errors[:10]:
            print(f"  {temporal[0]}..{temporal[1]}: {err[:150]}")
    if weeks_with_timeouts or weeks_with_download_errors or weeks_with_zero_granules:
        print("\nJust re-run this script again -- already-downloaded files are "
              "skipped, so it will only retry what's still missing.")


if __name__ == "__main__":
    main()
