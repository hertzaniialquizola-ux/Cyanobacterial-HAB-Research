# combine_nhdplus_states.R
# Run this AFTER get_nhdplus_national.R finishes pulling every state.
#
# Merges all per-state waterbody files into one national file and
# de-duplicates lakes that straddle a state line (a lake near, say, the
# Ohio/Michigan border can get returned by both neighboring states' AOI
# queries) -- COMID is NHDPlus's unique waterbody identifier, so a lake
# appearing under more than one state collapses to a single row.

library(sf)
library(dplyr)

files <- list.files("data/raw/nhdplus/by_state", pattern = "\\.gpkg$", full.names = TRUE)

if (length(files) == 0) {
  stop("No per-state files found in data/raw/nhdplus/by_state -- run ",
       "get_nhdplus_national.R first.")
}

message("Combining ", length(files), " state file(s)...")
all_wb <- do.call(rbind, lapply(files, sf::st_read, quiet = TRUE))

n_before <- nrow(all_wb)
all_wb <- all_wb |> dplyr::distinct(comid, .keep_all = TRUE)
n_after <- nrow(all_wb)

sf::st_write(all_wb, "data/raw/nhdplus/lake_national.gpkg", delete_dsn = TRUE)
message("Wrote ", n_after, " unique national waterbodies to ",
        "data/raw/nhdplus/lake_national.gpkg (", n_before - n_after,
        " cross-state duplicate(s) removed).")
