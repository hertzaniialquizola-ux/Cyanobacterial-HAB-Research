import earthaccess

earthaccess.login()

# 2020-08-01..07 (not a live/current date) so this smoke-test pull actually
# falls inside CDC WONDER's 2018-2024 mortality range -- the original
# 2026-08-01..07 window pulled "this week" as a live test, which downloaded
# fine but produced zero county-year matches when joined against mortality
# (WONDER simply has no 2026 data). Same bounding box as the NHDPlus smoke
# test (Toledo/western Lake Erie). These files land alongside the original
# 2026 pull in data/raw/cyan/ -- no collision, since granule filenames
# encode their own date.
results = earthaccess.search_data(
    concept_id="C3416412382-OB_CLOUD",
    temporal=("2020-08-01", "2020-08-07"),
    bounding_box=(-83.7, 41.5, -82.6, 42.3),
)

print(f"Found {len(results)} granules")
files = earthaccess.download(results, "data/raw/cyan")
print(files)
