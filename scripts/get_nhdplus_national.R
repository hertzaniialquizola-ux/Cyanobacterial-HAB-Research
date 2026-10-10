# get_nhdplus_national.R
# Scale-up of the Toledo smoke-test NHDPlus pull to nationwide lake coverage.
#
# download_nhdplusv2() doesn't support a regional/bbox filter (confirmed
# during the smoke test -- it only offers one multi-GB national 7z file), so
# this loops the same get_waterbodies(AOI=...) call the smoke test validated
# over each US state instead of trying one nationwide AOI in a single call,
# since a single AOI covering the whole country risks a timeout or memory
# blowup on the live query service the smoke test used.
#
# Safe to interrupt and re-run: any state whose output file already exists
# is skipped, so a killed/crashed run just resumes with the next state.
#
# After this finishes, run combine_nhdplus_states.R to merge the per-state
# files into the single data/raw/nhdplus/lake_national.gpkg that
# build_county_year_panel.py looks for.

library(nhdplusTools)
library(sf)
sf::sf_use_s2(FALSE)  # NHDPlus SC returns a spherically-invalid ring; GEOS planar mode tolerates it
sf::sf_use_s2(FALSE)  # NHDPlus SC returns a spherically-invalid ring; GEOS planar mode tolerates it
library(dplyr)

if (!requireNamespace("tigris", quietly = TRUE)) install.packages("tigris")
library(tigris)
options(tigris_use_cache = TRUE)

out_dir <- "data/raw/nhdplus/by_state"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# CONUS only -- matches CyAN's CONUS-only coverage, so pulling AK/HI/territory
# lakes would just be lakes with no possible exposure data to join against.
us_states <- tigris::states(cb = TRUE, resolution = "20m") |>
  sf::st_transform(4326) |>
  dplyr::filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "AS", "MP"))

message("Pulling waterbodies for ", nrow(us_states), " states/DC...")

for (i in seq_len(nrow(us_states))) {
  state_abbr <- us_states$STUSPS[i]
  out_file <- file.path(out_dir, paste0("waterbodies_", state_abbr, ".gpkg"))

  if (file.exists(out_file)) {
    message(state_abbr, ": already pulled, skipping")
    next
  }

  message(state_abbr, ": querying...")
  aoi <- sf::st_geometry(sf::st_make_valid(us_states[i, ]))

  wb <- tryCatch(
    nhdplusTools::get_waterbodies(AOI = aoi),
    error = function(e) {
      message(state_abbr, ": FAILED - ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(wb) || nrow(wb) == 0) {
    message(state_abbr, ": 0 waterbodies returned")
    next
  }

  sf::st_write(wb, out_file, delete_dsn = TRUE, quiet = TRUE)
  message(state_abbr, ": wrote ", nrow(wb), " waterbodies")
}

message("Done. Run combine_nhdplus_states.R next to merge everything into ",
        "data/raw/nhdplus/lake_national.gpkg.")
