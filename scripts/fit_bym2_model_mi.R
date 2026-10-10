# ============================================================================
# fit_bym2_model_mi.R
#
# BYM2 model with suppression handled via covariate-informed multiple
# imputation (Option B from the suppression-handling discussion): CDC
# WONDER's suppressed cells have a true count known to be an integer in
# [1,9] (that's WONDER's own suppression threshold -- true zeros already
# appear as explicit 0s in the panel, per build_county_year_panel.py, so
# every suppressed row really is a nonzero count between 1 and 9).
#
# For each suppressed county-year, we:
#   1. Fit a Poisson count model on the NON-suppressed rows, using
#      poverty_rate, median_household_income, well_water_fraction, and cei
#      as predictors (population as an offset) -- this is the "similar
#      counties" part: a poorer, higher-exposure county's predicted rate
#      informs where in [1,9] its suppressed count most likely sits.
#   2. Convert that county's fitted rate into a set of probabilities over
#      the 9 possible true values (1..9), via a Poisson pmf normalized to
#      sum to 1 over just that range (since we know the true value is in
#      that range, not the full Poisson support).
#   3. Draw M = 20 completed datasets, each with a random draw per
#      suppressed cell from those cell-specific probabilities.
#   4. Fit BYM2 on each of the M completed datasets.
#   5. Pool the M fits with Rubin's rules for the fixed effects, and report
#      the range of DIC/WAIC/hyperparameters across imputations (Rubin's
#      rules don't strictly apply to DIC/WAIC/hyperparameters the way they
#      do to fixed effects -- report those as a range, not a pooled point
#      estimate, and say so in the manuscript).
#
# This is the "most likely value, informed by similar counties, without
# pretending we're certain" version of suppression handling -- more
# defensible than single point-imputation (Rubin's rules propagate the
# extra imputation uncertainty into the final CIs) and simpler to implement
# and explain than the interval-censored-likelihood (cenpoisson2) approach.
#
# PREREQUISITE: run scripts/get_census_poverty.py and re-run
# scripts/build_county_year_panel.py FIRST so poverty_rate and
# median_household_income are actually columns in county_year_panel_full.csv.
# This script checks for those columns and stops with a clear message if
# they're missing rather than silently running without them.
#
# Run from the project root:
#   cd ~/Projects/Cyanobacterial-HAB-Research
#   Rscript scripts/fit_bym2_model_mi.R
#
# Writes: outputs/model_mi/pooled_fixed_effects.txt
#         outputs/model_mi/per_imputation_hyperpar_and_fit.csv
#         outputs/model_mi/imputed_datasets.rds   (all M completed datasets, for reuse/audit)
#         outputs/model_mi/bym2_mi_fits.rds       (all M fitted BYM2 INLA objects)
#         outputs/model_mi/iid_mi_fits.rds        (all M fitted non-spatial IID comparator objects, added 2026-09-29)
#
# 2026-09-29: added lst_day_c (daytime lake-surface temperature, rescaled to
# lst_day_c_10) as a fixed effect, and added a non-spatial IID comparator
# fit on every imputation (mirroring fit_bym2_model.R's baseline comparator,
# but now inside the MI loop) so DIC/WAIC has a spatial-vs-non-spatial test
# under the same suppression handling as the headline results.
# ============================================================================

set.seed(20260917)  # reproducible imputation draws -- report this seed in the manuscript
M <- 20              # number of imputations; 20 is a reasonable default given ~50% missingness

suppressed_pkgs <- character(0)
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) suppressed_pkgs <<- c(suppressed_pkgs, pkg)
}
for (p in c("sf", "spdep", "dplyr", "INLA")) check_pkg(p)
if (length(suppressed_pkgs) > 0) {
  if ("INLA" %in% suppressed_pkgs) {
    message('Install INLA: install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)')
    suppressed_pkgs <- setdiff(suppressed_pkgs, "INLA")
  }
  if (length(suppressed_pkgs) > 0) message("Install: install.packages(c(", paste(sprintf('"%s"', suppressed_pkgs), collapse = ", "), "))")
  stop("Install the packages above, then re-run this script.")
}

suppressWarnings({
  library(sf)
  library(spdep)
  library(dplyr)
  library(INLA)
})

# NOTE: this was previously (incorrectly) diagnosed as a macOS/Apple Silicon
# threading race condition in INLA's internal GMRFLib sorting routine, and
# forced to single-threaded execution as a workaround. That did NOT fix the
# crash -- it reproduced identically single-threaded, which a genuine race
# condition would not do. The actual cause was an unscaled covariate
# (median_household_income, raw dollars, ~$18k-$170k) creating a huge scale
# mismatch against the other [0,1]-range covariates in the same formula --
# see the median_household_income_10k fix below. Left single-threaded for
# this confirmation run; safe to try inla.setOption(num.threads = "4:1")
# again for speed once a run completes cleanly with the scaling fix.
inla.setOption(num.threads = "1:1")

dir.create("outputs/model_mi", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load panel + boundaries, build the adjacency graph (same as the
#    complete-case baseline script -- kept independent/self-contained here
#    on purpose so this script can be read and run without fit_bym2_model.R)
# ---------------------------------------------------------------------------
panel <- read.csv("outputs/county_year_panel_full.csv", colClasses = c(county_fips = "character"))

required_cols <- c("poverty_rate", "median_household_income")
missing_cols <- setdiff(required_cols, names(panel))
if (length(missing_cols) > 0) {
  stop(
    "county_year_panel_full.csv is missing: ", paste(missing_cols, collapse = ", "), "\n",
    "Run these two steps first, then re-run this script:\n",
    "  1. export CENSUS_API_KEY=\"<your key>\"; python3 scripts/get_census_poverty.py\n",
    "  2. python3 scripts/build_county_year_panel.py"
  )
}

counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
excluded_state_fips <- c("02", "15", "72", "78", "66", "60", "69")
counties <- counties[!(counties$STATEFP %in% excluded_state_fips), ]
counties$GEOID <- as.character(counties$GEOID)
# Reproject to CONUS Albers (meters) -- see fit_bym2_model.R for why (poly2nb
# snap tolerance and centroid distances are meaningless in geographic degrees).
counties <- st_transform(counties, 5070)

county_ids <- sort(intersect(unique(panel$county_fips), counties$GEOID))
counties <- counties[match(county_ids, counties$GEOID), ]
stopifnot(all(counties$GEOID == county_ids))

# Build the adjacency graph and force full connectivity -- see the detailed
# comment in fit_bym2_model.R (same block, kept self-contained here). BYM2's
# scale.model=TRUE crashes INLA outright on a disconnected graph rather than
# failing gracefully, so every disconnected component gets stitched to its
# nearest neighbor by centroid distance before the graph is written out.
nb <- poly2nb(counties, queen = TRUE, snap = 1000)  # 1km snap, in meters post-reprojection

force_connect_nb <- function(nb, coords) {
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
  cat(sprintf("Added %d edge(s) to force full connectivity (nearest-centroid stitching) -- see outputs/model_mi/connectivity_fixes.csv\n", nrow(edge_log)))
  connectivity_report <- data.frame(
    county_a = county_ids[edge_log$from_idx],
    county_b = county_ids[edge_log$to_idx],
    distance_km = round(edge_log$distance_km, 1)
  )
  write.csv(connectivity_report, "outputs/model_mi/connectivity_fixes.csv", row.names = FALSE)
}
stopifnot(n.comp.nb(nb)$nc == 1)

nb2INLA("outputs/model_mi/county_graph.adj", nb)
g <- inla.read.graph("outputs/model_mi/county_graph.adj")
county_idx_map <- data.frame(county_fips = county_ids, idx = seq_along(county_ids))

panel <- panel %>%
  filter(county_fips %in% county_ids) %>%
  mutate(
    population = as.numeric(population),
    cei = as.numeric(cei),
    lst_day_c = as.numeric(lst_day_c),
    deaths = as.numeric(deaths),  # NA for suppressed rows -- that's expected
    deaths_suppressed = as.logical(deaths_suppressed)
  ) %>%
  filter(!is.na(population), population > 0) %>%
  # Restrict to county-years with a real CyAN lake match. A missing cei means
  # "no qualifying lake mapped nearby," not "confirmed zero bloom exposure" --
  # treating it as a hard zero would assume the lake-match process captures
  # every real exposure source, which isn't defensible (unmapped/sub-threshold
  # water bodies could still matter). This also matches fit_bym2_model.R's
  # baseline scope, which drops these rows too (its complete-case filter
  # requires !is.na(cei)) -- keeping both scripts on the same analytic sample
  # is what makes their DIC/WAIC actually comparable. Chosen 2026-09-19,
  # after the alternative (keep all CONUS counties, cei=0 for no-lake)
  # produced a well_water_fraction estimate that didn't replicate once
  # fit_bym2_model.R's smaller lake-only sample was checked against it.
  filter(!is.na(cei))

# LST (lake thermal regime) added as a covariate 2026-09-29. MODIS's coverage
# doesn't perfectly overlap CyAN's (project doc section 15: 7,966 of 8,155
# cei-matched rows have a lst_day_c value) -- same "missing means no data,
# not a confirmed value" logic as the cei filter above, so these rows are
# dropped rather than imputed or zero-filled.
n_before_lst <- nrow(panel)
panel <- panel %>% filter(!is.na(lst_day_c))
cat(sprintf(
  "LST covariate: dropped %d of %d cei-matched rows with no lst_day_c match (MODIS/CyAN coverage mismatch) -- %d rows remain.\n",
  n_before_lst - nrow(panel), n_before_lst, nrow(panel)
))

panel <- panel %>%
  inner_join(county_idx_map, by = "county_fips") %>%
  mutate(
    year_idx = year - min(year) + 1,
    log_pop = log(population),
    # Rescaled the same way median_household_income was below: lst_day_c's
    # raw range (~14-50 degC) is close enough to the other covariates' scale
    # that it's very unlikely to trigger the same GMRFLib crash on its own,
    # but dividing by 10 brings it in line with them anyway and makes the
    # fitted coefficient read as "RR per 10 degC of daytime lake-surface
    # temperature" -- consistent with how median_household_income_10k is
    # reported.
    lst_day_c_10 = lst_day_c / 10
  )

# Numerical stability fix: median_household_income spans ~$18k-$170k while
# the other fixed effects (cei, well_water_fraction, poverty_rate) are all
# confined to [0,1] -- passing raw dollar amounts straight into INLA's
# linear predictor alongside those creates a ~5-6 order-of-magnitude scale
# mismatch between covariates in the same formula. That mismatch is a
# documented cause of internal numerical overflow in GMRFLib's sparse
# index/value sorting routine (R-INLA discussion group: "Assertion failed
# ... GMRFLib_idxval_nsort_x_core", root-caused to unscaled covariates and
# fixed there by scaling). This is almost certainly the actual cause of the
# identical crash seen here on every attempt -- NOT a threading race
# condition (a real race condition would not reproduce identically under
# single-threaded execution, which it did). Rescaling to tens of thousands
# of dollars brings this covariate to the same order of magnitude as the
# others and keeps its coefficient directly interpretable as "per $10,000
# of county median household income".
panel$median_household_income_10k <- panel$median_household_income / 10000

cat(sprintf(
  "Panel: %d rows total, %d suppressed (%.1f%%), %d usable for the imputation predictor model.\n",
  nrow(panel), sum(panel$deaths_suppressed), 100 * mean(panel$deaths_suppressed),
  sum(!panel$deaths_suppressed)
))

# ---------------------------------------------------------------------------
# 2. Fit the imputation (predictor) model on non-suppressed rows only.
#    This is a plain glm, not INLA -- it just needs to produce a fitted
#    rate per suppressed county-year to build the truncated draw
#    probabilities in step 3. Keep it simple and inspectable.
# ---------------------------------------------------------------------------
impute_model <- glm(
  deaths ~ poverty_rate + median_household_income + well_water_fraction + cei + offset(log_pop),
  family = poisson(link = "log"),
  data = panel %>% filter(!deaths_suppressed)
)
cat("\n--- Imputation predictor model (fit on non-suppressed rows) ---\n")
print(summary(impute_model))

suppressed_rows <- panel %>% filter(deaths_suppressed)
if (any(is.na(suppressed_rows$poverty_rate) | is.na(suppressed_rows$median_household_income))) {
  n_bad <- sum(is.na(suppressed_rows$poverty_rate) | is.na(suppressed_rows$median_household_income))
  warning(sprintf(
    "%d suppressed rows have missing poverty_rate/median_household_income (ACS nonresponse for that county) -- these fall back to the marginal [1,9] distribution below, not a covariate-informed one. Inspect these counties before trusting their imputed values.",
    n_bad
  ))
}
fitted_rate <- predict(impute_model, newdata = suppressed_rows, type = "response")
# Guard against a completely missing prediction (NA covariates): fall back
# to the empirical suppressed-cell rate implied by the model's average
# suppressed-cell population, rather than crashing.
fallback_rate <- mean(fitted_rate, na.rm = TRUE)
fitted_rate[is.na(fitted_rate)] <- fallback_rate

# ---------------------------------------------------------------------------
# 3. Build per-cell probabilities over the known censoring range {1..9} and
#    draw M completed datasets.
# ---------------------------------------------------------------------------
censor_range <- 1:9
draw_probs <- t(sapply(fitted_rate, function(lam) {
  p <- dpois(censor_range, lambda = lam)
  p / sum(p)
}))

imputed_datasets <- vector("list", M)
for (m in seq_len(M)) {
  draws <- apply(draw_probs, 1, function(p) sample(censor_range, size = 1, prob = p))
  panel_m <- panel
  panel_m$deaths[panel_m$deaths_suppressed] <- draws
  imputed_datasets[[m]] <- panel_m
}
saveRDS(imputed_datasets, "outputs/model_mi/imputed_datasets.rds")
cat(sprintf("\nBuilt %d imputed datasets.\n", M))

# ---------------------------------------------------------------------------
# 4. Fit BYM2 on each imputed dataset
# ---------------------------------------------------------------------------
pc_prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)
formula_bym2 <- deaths ~ 1 + cei + well_water_fraction + poverty_rate + median_household_income_10k + lst_day_c_10 +
  f(idx, model = "bym2", graph = g, scale.model = TRUE, hyper = pc_prior_bym2) +
  f(year_idx, model = "iid")

# Non-spatial IID comparator, added 2026-09-29 to mirror fit_bym2_model.R's
# baseline-vs-IID comparison (project doc section 16, next-steps item c) --
# same fixed effects, f(idx, model="iid") in place of the BYM2 spatial term.
# Fit on every one of the M imputed datasets below (not just once on a
# complete-case sample) so the spatial-vs-non-spatial DIC/WAIC test is done
# under the same missing-data handling used for the headline pooled results.
formula_iid <- deaths ~ 1 + cei + well_water_fraction + poverty_rate + median_household_income_10k + lst_day_c_10 +
  f(idx, model = "iid") +
  f(year_idx, model = "iid")

fits <- vector("list", M)
fits_iid <- vector("list", M)
for (m in seq_len(M)) {
  cat(sprintf("Fitting imputation %d / %d (BYM2 + non-spatial IID comparator)...\n", m, M))
  fits[[m]] <- inla(
    formula_bym2,
    family = "poisson",
    data = imputed_datasets[[m]],
    offset = log_pop,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  fits_iid[[m]] <- inla(
    formula_iid,
    family = "poisson",
    data = imputed_datasets[[m]],
    offset = log_pop,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
}
saveRDS(fits, "outputs/model_mi/bym2_mi_fits.rds")
saveRDS(fits_iid, "outputs/model_mi/iid_mi_fits.rds")

# ---------------------------------------------------------------------------
# 5. Pool fixed effects with Rubin's rules; report hyperparameters/fit
#    stats as a range across imputations (not pooled -- see header note)
# ---------------------------------------------------------------------------
rubin_pool <- function(estimates, variances, M) {
  # estimates, variances: numeric vectors of length M for ONE term
  qbar <- mean(estimates)
  ubar <- mean(variances)               # within-imputation variance
  b <- var(estimates)                   # between-imputation variance
  total_var <- ubar + (1 + 1 / M) * b
  se <- sqrt(total_var)
  data.frame(
    pooled_mean = qbar,
    pooled_se = se,
    ci_lower95 = qbar - 1.96 * se,
    ci_upper95 = qbar + 1.96 * se,
    frac_missing_info = ((1 + 1 / M) * b) / total_var
  )
}

fixed_terms <- rownames(fits[[1]]$summary.fixed)
pooled_list <- lapply(fixed_terms, function(term) {
  ests <- sapply(fits, function(f) f$summary.fixed[term, "mean"])
  vars <- sapply(fits, function(f) f$summary.fixed[term, "sd"]^2)
  cbind(term = term, rubin_pool(ests, vars, M))
})
pooled_fixed <- do.call(rbind, pooled_list)
pooled_fixed$RR <- exp(pooled_fixed$pooled_mean)
pooled_fixed$RR_lower95 <- exp(pooled_fixed$ci_lower95)
pooled_fixed$RR_upper95 <- exp(pooled_fixed$ci_upper95)

fit_stats <- data.frame(
  imputation = seq_len(M),
  dic = sapply(fits, function(f) f$dic$dic),
  waic = sapply(fits, function(f) f$waic$waic),
  dic_iid = sapply(fits_iid, function(f) f$dic$dic),
  waic_iid = sapply(fits_iid, function(f) f$waic$waic),
  phi_mean = sapply(fits, function(f) f$summary.hyperpar["Phi for idx", "mean"]),
  sigma_mean = sapply(fits, function(f) 1 / sqrt(f$summary.hyperpar["Precision for idx", "mean"]))
)
write.csv(fit_stats, "outputs/model_mi/per_imputation_hyperpar_and_fit.csv", row.names = FALSE)

sink("outputs/model_mi/pooled_fixed_effects.txt")
cat("=== BYM2 with covariate-informed multiple imputation (M =", M, ") ===\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("Suppressed rows imputed: %d of %d total panel rows (%.1f%%)\n\n",
            sum(panel$deaths_suppressed), nrow(panel), 100 * mean(panel$deaths_suppressed)))
cat("--- Pooled fixed effects (Rubin's rules) ---\n")
print(pooled_fixed, row.names = FALSE)
cat("\n--- Fit statistics across imputations (NOT Rubin-pooled -- report as a range) ---\n")
cat(sprintf("BYM2 DIC:   %.1f to %.1f (mean %.1f)\n", min(fit_stats$dic), max(fit_stats$dic), mean(fit_stats$dic)))
cat(sprintf("BYM2 WAIC:  %.1f to %.1f (mean %.1f)\n", min(fit_stats$waic), max(fit_stats$waic), mean(fit_stats$waic)))
cat(sprintf("IID  DIC:   %.1f to %.1f (mean %.1f)\n", min(fit_stats$dic_iid), max(fit_stats$dic_iid), mean(fit_stats$dic_iid)))
cat(sprintf("IID  WAIC:  %.1f to %.1f (mean %.1f)\n", min(fit_stats$waic_iid), max(fit_stats$waic_iid), mean(fit_stats$waic_iid)))
cat(sprintf("Phi:   %.3f to %.3f (mean %.3f)\n", min(fit_stats$phi_mean), max(fit_stats$phi_mean), mean(fit_stats$phi_mean)))
cat(sprintf("Sigma: %.3f to %.3f (mean %.3f)\n", min(fit_stats$sigma_mean), max(fit_stats$sigma_mean), mean(fit_stats$sigma_mean)))
cat("\n--- Spatial vs. non-spatial, per-imputation DIC/WAIC delta (IID minus BYM2; positive = BYM2 fits better) ---\n")
cat(sprintf("Mean DIC delta:  %.1f\n", mean(fit_stats$dic_iid - fit_stats$dic)))
cat(sprintf("Mean WAIC delta: %.1f\n", mean(fit_stats$waic_iid - fit_stats$waic)))
sink()

cat("\nDone. See outputs/model_mi/pooled_fixed_effects.txt for the headline numbers.\n")
