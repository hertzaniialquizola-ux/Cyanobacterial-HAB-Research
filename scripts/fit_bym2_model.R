# ============================================================================
# fit_bym2_model.R
#
# BYM2 (Besag-York-Mollié 2) spatial disease-mapping model for the HAB
# CyanoHAB exposure -> hepatic/GI mortality study, fit via R-INLA.
#
# Reads:  outputs/county_year_panel_full.csv (from build_county_year_panel.py)
#         data/raw/census/county_boundaries/cb_2022_us_county_20m.shp
# Writes: outputs/model/bym2_baseline_summary.txt
#         outputs/model/bym2_baseline_fit.rds       (full INLA object, for later reuse)
#         outputs/model/county_rr_estimates.csv     (posterior county RRs)
#         outputs/model/morans_i_diagnostics.txt
#
# Run from the project root:
#   cd ~/Projects/Cyanobacterial-HAB-Research
#   Rscript scripts/fit_bym2_model.R
#
# NOTE ON SUPPRESSION HANDLING (read before you trust these numbers):
# 10,930 of 21,826 county-year rows (50.1%) have deaths_suppressed == TRUE and
# blank deaths/cei. Dr. Waller's explicit recommendation (project doc, section 10)
# was a covariate-informed (poverty/SES) approach to those cells, NOT dropping
# them and NOT flat substitution. This baseline script does the simplest
# defensible thing instead -- complete-case (drops suppressed rows) -- because
# no SES/poverty covariate is joined into the panel yet. That is a real
# limitation, not a final modeling choice: see the "SUPPRESSION TODO" block
# below for what to change once you've decided how to incorporate SES data or
# a censored likelihood (cenpoisson2). Do not present this baseline run's
# numbers in the manuscript without addressing that first.
# ============================================================================

suppressed_pkgs <- character(0)
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) suppressed_pkgs <<- c(suppressed_pkgs, pkg)
}
for (p in c("sf", "spdep", "dplyr", "INLA")) check_pkg(p)
if (length(suppressed_pkgs) > 0) {
  if ("INLA" %in% suppressed_pkgs) {
    message("INLA is not installed. Install it with:")
    message('  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)')
    suppressed_pkgs <- setdiff(suppressed_pkgs, "INLA")
  }
  if (length(suppressed_pkgs) > 0) {
    message("Missing CRAN packages: ", paste(suppressed_pkgs, collapse = ", "))
    message('  install.packages(c(', paste(sprintf('"%s"', suppressed_pkgs), collapse = ", "), '))')
  }
  stop("Install the packages above, then re-run this script.")
}

suppressWarnings({
  library(sf)
  library(spdep)
  library(dplyr)
  library(INLA)
})

# Work around a known crash on macOS (especially Apple Silicon) where INLA's
# internal GMRFLib sorting routine hits a threading race condition under the
# default multi-threaded execution, specifically with scale.model=TRUE (BYM2)
# -- it aborts with "Assertion failed: (0 == 1), function GMRFLib_idxval_nsort_x_core"
# rather than a normal R error. Forcing single-threaded execution avoids the
# race condition; it's slower but correct. If this run succeeds, you can try
# raising it later (e.g. inla.setOption(num.threads = "4:1")) for speed once
# you've confirmed the model itself is not the problem.
inla.setOption(num.threads = "1:1")

dir.create("outputs/model", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load panel data and county boundaries
# ---------------------------------------------------------------------------
panel <- read.csv("outputs/county_year_panel_full.csv", colClasses = c(county_fips = "character"))
panel$county_fips <- sprintf("%05s", panel$county_fips)
panel$county_fips <- gsub(" ", "0", panel$county_fips)

counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
# Restrict to CONUS + DC to match the rest of the pipeline (CyAN is CONUS-only)
excluded_state_fips <- c("02", "15", "72", "78", "66", "60", "69") # AK, HI, PR, VI, GU, AS, MP
counties <- counties[!(counties$STATEFP %in% excluded_state_fips), ]
counties$GEOID <- as.character(counties$GEOID)
# Reproject to CONUS Albers (meters) -- the raw shapefile is in geographic
# degrees, where poly2nb's snap tolerance is nearly meaningless and small
# generalization gaps between adjacent counties (this is the 20m/1:20M
# cartographic boundary file, not full-resolution TIGER) are more likely to
# register as non-touching. Matches the projection used elsewhere in the
# pipeline (build_county_year_panel.py uses EPSG:5070 for the same reason).
counties <- st_transform(counties, 5070)

# Keep only counties present in both the panel and the boundary file, in the
# SAME row order in both objects -- this order defines the adjacency graph
# index, and a mismatch here silently misassigns every spatial effect.
county_ids <- sort(intersect(unique(panel$county_fips), counties$GEOID))
counties <- counties[match(county_ids, counties$GEOID), ]
stopifnot(all(counties$GEOID == county_ids))
cat(sprintf("Counties in model: %d\n", length(county_ids)))

# ---------------------------------------------------------------------------
# 2. Build the adjacency graph (queen contiguity, first-order), then force
#    full connectivity.
#
# WHY: BYM2's ICAR component (scale.model=TRUE) requires the graph to be a
# single connected component -- INLA's own scaling routine crashes outright
# (GMRFLib_idxval_nsort_x_core assertion) on a disconnected graph rather
# than failing gracefully. Two real things cause disconnection here: (1)
# generalization gaps in the 20m cartographic boundary file between
# counties that truly share a border (partly addressed by the snap= below
# and the reprojection above), and (2) genuine geographic islands (e.g.
# Nantucket MA, San Juan County WA) with no land border to any other
# county. Both need to end up connected for the model to run at all, so
# rather than trying to tell the two apart, every remaining disconnected
# component is stitched to its nearest neighbor by centroid distance.
# ---------------------------------------------------------------------------
nb <- poly2nb(counties, queen = TRUE, snap = 1000)  # 1km snap, in meters post-reprojection

force_connect_nb <- function(nb, coords) {
  # coords: n x 2 matrix of projected (planar, meters) centroid coordinates,
  # same row order as nb. Iteratively connects the nearest pair of counties
  # across the two largest remaining disconnected components until the
  # whole graph is one component. Returns the patched nb plus a data.frame
  # log of every edge it had to add (for the connectivity report below).
  edge_log <- data.frame(from_idx = integer(0), to_idx = integer(0), distance_km = numeric(0))
  repeat {
    comp <- n.comp.nb(nb)
    if (comp$nc <= 1) break
    comp_id <- comp$comp.id
    idx1 <- which(comp_id == 1)
    idx_other <- which(comp_id != 1)
    dx <- outer(coords[idx1, 1], coords[idx_other, 1], "-")
    dy <- outer(coords[idx1, 2], coords[idx_other, 2], "-")
    d <- sqrt(dx^2 + dy^2)
    best <- which(d == min(d), arr.ind = TRUE)[1, ]
    i <- idx1[best[1]]
    j <- idx_other[best[2]]
    nb[[i]] <- if (identical(as.integer(nb[[i]]), 0L)) as.integer(j) else sort(unique(c(nb[[i]], j)))
    nb[[j]] <- if (identical(as.integer(nb[[j]]), 0L)) as.integer(i) else sort(unique(c(nb[[j]], i)))
    edge_log <- rbind(edge_log, data.frame(from_idx = i, to_idx = j, distance_km = d[best[1], best[2]] / 1000))
  }
  attr(nb, "edge_log") <- edge_log
  nb
}

centroid_coords <- st_coordinates(st_centroid(counties))
n_isolated_before <- sum(card(nb) == 0)
n_components_before <- n.comp.nb(nb)$nc
nb <- force_connect_nb(nb, centroid_coords)
edge_log <- attr(nb, "edge_log")

cat(sprintf("Adjacency graph: %d counties had zero neighbors, %d disconnected components before fixing.\n",
            n_isolated_before, n_components_before))
if (nrow(edge_log) > 0) {
  cat(sprintf("Added %d edge(s) to force full connectivity (nearest-centroid stitching):\n", nrow(edge_log)))
  connectivity_report <- data.frame(
    county_a = county_ids[edge_log$from_idx],
    county_b = county_ids[edge_log$to_idx],
    distance_km = round(edge_log$distance_km, 1)
  )
  print(connectivity_report)
  write.csv(connectivity_report, "outputs/model/connectivity_fixes.csv", row.names = FALSE)
  cat("  (full list written to outputs/model/connectivity_fixes.csv -- worth a one-line\n")
  cat("   mention in the manuscript's spatial-structure-robustness discussion, Section 6b/6c)\n")
}
stopifnot(n.comp.nb(nb)$nc == 1)  # must be fully connected before nb2INLA/INLA will run

nb2INLA("outputs/model/county_graph.adj", nb)
g <- inla.read.graph("outputs/model/county_graph.adj")

county_idx_map <- data.frame(county_fips = county_ids, idx = seq_along(county_ids))

# ---------------------------------------------------------------------------
# 3. SUPPRESSION TODO -- replace this block once the suppression-handling
#    decision (manuscript doc, Section 3 checklist item, project doc section
#    10) is made. Current behavior: complete-case (drop suppressed rows).
# ---------------------------------------------------------------------------
panel_model <- panel %>%
  filter(deaths_suppressed == "False" | deaths_suppressed == FALSE) %>%
  mutate(
    deaths = as.numeric(deaths),
    cei = as.numeric(cei),
    population = as.numeric(population)
  ) %>%
  filter(!is.na(deaths), !is.na(cei), !is.na(population), population > 0) %>%
  inner_join(county_idx_map, by = "county_fips")

cat(sprintf("Rows used (complete-case, non-suppressed): %d of %d (%.1f%%)\n",
            nrow(panel_model), nrow(panel), 100 * nrow(panel_model) / nrow(panel)))

# Year as a simple random effect (RW1 or IID) alongside the spatial term;
# start with IID -- swap to model="rw1" if you want to allow smooth drift.
panel_model$year_idx <- panel_model$year - min(panel_model$year) + 1
# BYM2 needs one spatial index per row (INLA replicates internally); expected
# deaths (offset) uses population so the fitted term is a relative risk.
panel_model$log_pop <- log(panel_model$population)

# ---------------------------------------------------------------------------
# 4. PC priors (see manuscript doc Section 3c) -- P(sigma>1)=0.01, P(phi>0.5)=0.5
# ---------------------------------------------------------------------------
pc_prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

# ---------------------------------------------------------------------------
# 5. Fit the baseline BYM2 model
#    deaths_ij ~ Poisson(population_ij * exp(beta0 + beta1*cei_ij + beta2*well_water_ij + u_i + year_j))
# ---------------------------------------------------------------------------
formula_bym2 <- deaths ~ 1 + cei + well_water_fraction +
  f(idx, model = "bym2", graph = g, scale.model = TRUE, hyper = pc_prior_bym2) +
  f(year_idx, model = "iid")

fit_bym2 <- inla(
  formula_bym2,
  family = "poisson",
  data = panel_model,
  offset = log_pop,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)
)

saveRDS(fit_bym2, "outputs/model/bym2_baseline_fit.rds")

# ---------------------------------------------------------------------------
# 6. Non-spatial IID comparator (for the DIC/WAIC comparison in Section 3d)
# ---------------------------------------------------------------------------
formula_iid <- deaths ~ 1 + cei + well_water_fraction +
  f(idx, model = "iid") +
  f(year_idx, model = "iid")

fit_iid <- inla(
  formula_iid,
  family = "poisson",
  data = panel_model,
  offset = log_pop,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# ---------------------------------------------------------------------------
# 7. Write a plain-text summary (fixed effects, hyperparameters, fit stats,
#    CPO flags) -- paste these numbers into the Results template in the
#    manuscript doc (Section 3d)
# ---------------------------------------------------------------------------
sink("outputs/model/bym2_baseline_summary.txt")
cat("=== BYM2 baseline model summary ===\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("Rows used: %d / %d total panel rows (suppressed rows dropped -- see SUPPRESSION TODO in script)\n\n",
            nrow(panel_model), nrow(panel)))

cat("--- Fixed effects (log-RR scale) ---\n")
print(fit_bym2$summary.fixed)
cat("\n--- Fixed effects (RR scale, exp of above) ---\n")
fe <- fit_bym2$summary.fixed
fe_rr <- data.frame(
  term = rownames(fe),
  RR = exp(fe$mean),
  RR_lower95 = exp(fe$`0.025quant`),
  RR_upper95 = exp(fe$`0.975quant`)
)
print(fe_rr)

cat("\n--- BYM2 hyperparameters ---\n")
print(fit_bym2$summary.hyperpar)

cat("\n--- Model fit: BYM2 vs. non-spatial IID comparator ---\n")
cat(sprintf("BYM2    DIC = %.1f   WAIC = %.1f\n", fit_bym2$dic$dic, fit_bym2$waic$waic))
cat(sprintf("IID     DIC = %.1f   WAIC = %.1f\n", fit_iid$dic$dic, fit_iid$waic$waic))

cat("\n--- CPO diagnostics ---\n")
mean_log_cpo <- mean(log(fit_bym2$cpo$cpo), na.rm = TRUE)
cat(sprintf("Mean log(CPO): %.4f\n", mean_log_cpo))
low_cpo_cutoff <- quantile(fit_bym2$cpo$cpo, 0.01, na.rm = TRUE)
n_low_cpo <- sum(fit_bym2$cpo$cpo < low_cpo_cutoff, na.rm = TRUE)
cat(sprintf("Rows below 1st percentile CPO (%.5f): %d -- inspect these as potential influential observations\n",
            low_cpo_cutoff, n_low_cpo))
sink()

# ---------------------------------------------------------------------------
# 8. County-level posterior RR estimates -> CSV (for the choropleth map)
# ---------------------------------------------------------------------------
bym2_re <- fit_bym2$summary.random$idx[1:length(county_ids), ]
county_rr <- data.frame(
  county_fips = county_ids,
  posterior_mean_log_rr = bym2_re$mean,
  rr = exp(bym2_re$mean),
  rr_lower95 = exp(bym2_re$`0.025quant`),
  rr_upper95 = exp(bym2_re$`0.975quant`)
)
write.csv(county_rr, "outputs/model/county_rr_estimates.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. Moran's I diagnostics: pre-model (raw SMR) and post-model (BYM2 residual)
# ---------------------------------------------------------------------------
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

smr_by_county <- panel_model %>%
  group_by(county_fips) %>%
  summarise(observed = sum(deaths), expected = sum(population) * (sum(panel_model$deaths) / sum(panel_model$population))) %>%
  mutate(smr = observed / expected) %>%
  right_join(county_idx_map, by = "county_fips") %>%
  arrange(idx)

moran_raw <- moran.test(smr_by_county$smr, lw, zero.policy = TRUE, na.action = na.exclude)

fitted_vals <- fit_bym2$summary.fitted.values$mean
resid_bym2 <- panel_model$deaths - fitted_vals
resid_by_county <- tapply(resid_bym2, panel_model$idx, mean, na.rm = TRUE)
resid_by_county_full <- rep(NA_real_, length(county_ids))
resid_by_county_full[as.integer(names(resid_by_county))] <- resid_by_county
moran_resid <- moran.test(resid_by_county_full, lw, zero.policy = TRUE, na.action = na.exclude)

sink("outputs/model/morans_i_diagnostics.txt")
cat("=== Pre-model: Moran's I on raw county SMRs ===\n")
print(moran_raw)
cat("\n=== Post-model: Moran's I on BYM2 residuals ===\n")
print(moran_resid)
sink()

cat("\nDone. Outputs written to outputs/model/:\n")
cat("  bym2_baseline_summary.txt   -- fixed effects, hyperparameters, DIC/WAIC, CPO\n")
cat("  morans_i_diagnostics.txt    -- pre/post-model spatial autocorrelation\n")
cat("  county_rr_estimates.csv     -- posterior county-level RR (for the choropleth)\n")
cat("  bym2_baseline_fit.rds       -- full INLA object (reload with readRDS() for further work)\n")
