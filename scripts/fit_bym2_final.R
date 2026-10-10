# ============================================================================
# fit_bym2_final.R  --  final analysis (replaces fit_bym2_model_mi.R for the paper)
#
# Changes vs. the MI model (all in response to the 2026-10-02 critique):
#  1. EXPOSURE rebuilt from correctly-masked CyAN pixels (see build_alt_cei.py).
#     PRE-SPECIFIED primary exposure (chosen 2026-10-02, BEFORE seeing results):
#     px_share_20k = annual mean share of lake pixels >= 20,000 cells/mL, per IQR.
#     Secondary: corrected cei, wk_frac_5pct_20k, log10_cells_mean. All reported.
#  2. SUPPRESSION handled by a censored likelihood (INLA 'cenpoisson', interval
#     [1,9]) -- no imputation, no circularity, no Rubin pooling needed.
#  3. WATER PATHWAY: public-supply surface-water fraction (+ well fraction kept
#     as a covariate), and exposure x surface-supply interaction.
#  4. Added covariates: % 65+, % Black, % Hispanic, (optional) RUCC rurality.
#     NOTE: % 65+ is a covariate adjustment, not true age standardization.
#  5. All-county model with a lake-county indicator (non-lake reference group).
#  6. Effects reported as RR per IQR (exposure) / per 10 percentage points
#     (shares) with 95% credible intervals + P(RR in 0.95-1.05) as an
#     equivalence-style summary of the null.
#
# Run from project root:  Rscript scripts/fit_bym2_final.R
# Writes outputs/model_final/{results_table.csv, summary.txt, fits.rds}
# Runtime: ~8-10 BYM2 fits (each similar to one fit of the old MI loop).
# ============================================================================
suppressWarnings({ library(sf); library(spdep); library(dplyr); library(INLA) })
inla.setOption(num.threads = "1:1")
dir.create("outputs/model_final", showWarnings = FALSE, recursive = TRUE)
CENS <- c(1, 9)

# ---- 0. self-test: does cenpoisson recover a known coefficient on YOUR INLA build? --------
set.seed(1); n <- 6000; x <- rnorm(n); y <- rpois(n, exp(2.0 + 0.30 * x))
cens <- y >= CENS[1] & y <= CENS[2]; ys <- y; ys[cens] <- CENS[1]
st <- inla(ys ~ x, family = "cenpoisson", data = data.frame(ys, x),
           control.family = list(cenpoisson.I = CENS))
b <- st$summary.fixed["x", "mean"]
cat(sprintf("SELF-TEST: true beta 0.30, cenpoisson estimate %.3f (%.0f%% of rows censored)\n", b, 100 * mean(cens)))
if (abs(b - 0.30) > 0.06) stop("cenpoisson self-test failed -- do not trust results; send me this output.")

# ---- 1. data -------------------------------------------------------------------------
panel <- read.csv("outputs/county_year_panel_v2.csv", colClasses = c(county_fips = "character"))
panel$deaths_suppressed <- as.logical(panel$deaths_suppressed)   # CSV stores True/False text
counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
counties <- counties[!(counties$STATEFP %in% c("02","15","72","78","66","60","69")), ]
county_ids <- sort(intersect(unique(panel$county_fips), as.character(counties$GEOID)))
g <- inla.read.graph("outputs/model_mi/county_graph.adj")
if (g$n != length(county_ids)) stop(sprintf("graph has %d nodes but %d counties -- re-run fit_bym2_model_mi.R's graph step", g$n, length(county_ids)))
idx_map <- data.frame(county_fips = county_ids, idx = seq_along(county_ids))

has_rucc <- "rucc" %in% names(panel)
panel <- panel %>% filter(county_fips %in% county_ids, !is.na(population), population > 0) %>%
  inner_join(idx_map, by = "county_fips") %>%
  mutate(year_idx = year - min(year) + 1, log_pop = log(as.numeric(population)),
         y = ifelse(deaths_suppressed, CENS[1], deaths),
         ww10 = well_water_fraction * 10, sw10 = surface_supply_fraction * 10,
         pov10 = poverty_rate * 10, inc10k = median_household_income / 1e4, lst10 = lst_day_c / 10,
         age5 = pct_65plus * 20, blk10 = pct_black * 10, hsp10 = pct_hispanic * 10,
         lake_county = as.integer(!is.na(cei)))
if (has_rucc) panel$rucc_f <- factor(panel$rucc)
stopifnot(!any(is.na(panel$y)))

covs <- c("ww10", "sw10", "pov10", "inc10k", "age5", "blk10", "hsp10")
if (has_rucc) covs <- c(covs, "rucc_f")
cat(sprintf("Rows %d | suppressed (censored) %.1f%% | rurality covariate: %s\n", nrow(panel), 100 * mean(panel$deaths_suppressed), has_rucc))

pc <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
fit_one <- function(dat, rhs, spatial = TRUE) {
  sp <- if (spatial) "f(idx, model='bym2', graph=g, scale.model=TRUE, hyper=pc)" else "f(idx, model='iid')"
  fm <- as.formula(paste("y ~ 1 +", rhs, "+", sp, "+ f(year_idx, model='iid')"))
  inla(fm, family = "cenpoisson", data = dat, offset = log_pop,
       control.family = list(cenpoisson.I = CENS),
       control.compute = list(dic = TRUE, waic = TRUE))
}
pull <- function(fit, term, label, scale_note) {
  s <- fit$summary.fixed[term, ]; m <- fit$marginals.fixed[[term]]
  p0 <- inla.pmarginal(log(1.05), m) - inla.pmarginal(log(0.95), m)
  data.frame(model = label, term = term, unit = scale_note, RR = exp(s$mean),
             lo95 = exp(s[["0.025quant"]]), hi95 = exp(s[["0.975quant"]]),
             P_RR_in_0.95_1.05 = p0, DIC = fit$dic$dic, WAIC = fit$waic$waic, row.names = NULL)
}

# ---- 2. lake-matched analytic sample ---------------------------------------------------
lake <- panel %>% filter(lake_county == 1, !is.na(lst_day_c), !is.na(px_share_20k), !is.na(surface_supply_fraction),
                         if_all(c(pov10, inc10k, age5, blk10, hsp10, ww10), ~ !is.na(.)))
if (has_rucc) lake <- lake %>% filter(!is.na(rucc_f))
cat(sprintf("Lake-matched analytic sample: %d rows, %d counties\n", nrow(lake), n_distinct(lake$county_fips)))

covs_lake <- c(covs, "lst10")   # LST only exists for lake-matched rows
exposures <- c(px_share_20k = "PRIMARY: share of lake pixels >= 20k cells/mL (per IQR)",
               cei = "corrected CEI: fraction of weeks county-mean >= 20k (per IQR)",
               wk_frac_5pct_20k = "fraction of weeks with >=5% of pixels >= 20k (per IQR)",
               log10_cells_mean = "mean log10 cells/mL (per IQR)")
res <- list(); fits <- list()
for (ex in names(exposures)) {
  d <- lake %>% filter(!is.na(.data[[ex]]))
  iqr <- IQR(d[[ex]]); if (!is.finite(iqr) || iqr == 0) { cat(sprintf("%s: IQR is 0 -- skipped (still saturated/degenerate)\n", ex)); next }
  d$expo_z <- (d[[ex]] - mean(d[[ex]])) / iqr
  cat(sprintf("Fitting %s (IQR = %.3f, n = %d)...\n", ex, iqr, nrow(d)))
  f <- fit_one(d, paste(c("expo_z", covs_lake), collapse = " + "))
  fits[[ex]] <- f
  res[[length(res) + 1]] <- pull(f, "expo_z", paste0("BYM2: ", ex), "RR per IQR")
  if (ex == "px_share_20k") {
    res[[length(res) + 1]] <- pull(f, "ww10", "BYM2 primary", "RR per 10 pp well users")
    res[[length(res) + 1]] <- pull(f, "sw10", "BYM2 primary", "RR per 10 pp surface-water supply")
    res[[length(res) + 1]] <- pull(f, "age5", "BYM2 primary", "RR per 5 pp age 65+")
    res[[length(res) + 1]] <- pull(f, "pov10", "BYM2 primary", "RR per 10 pp poverty")
    res[[length(res) + 1]] <- pull(f, "inc10k", "BYM2 primary", "RR per $10k income")
    res[[length(res) + 1]] <- pull(f, "lst10", "BYM2 primary", "RR per 10 C (land-surface proxy; weak covariate)")
    # non-spatial comparator and interaction on the primary exposure
    fi <- fit_one(d, paste(c("expo_z", covs_lake), collapse = " + "), spatial = FALSE)
    fits[["primary_iid"]] <- fi
    res[[length(res) + 1]] <- pull(fi, "expo_z", "Non-spatial IID: primary", "RR per IQR")
    fx <- fit_one(d, paste(c("expo_z", "expo_z:sw10", covs_lake), collapse = " + "))
    fits[["primary_interaction"]] <- fx
    res[[length(res) + 1]] <- pull(fx, "expo_z:sw10", "BYM2 interaction: exposure x surface-supply (per 10 pp)", "ratio of RRs")
  }
}

# ---- 3. all-county model with lake indicator (non-lake reference group) ----------------
ex <- "px_share_20k"
da <- panel %>% filter(if_all(c(pov10, inc10k, age5, blk10, hsp10, ww10, sw10), ~ !is.na(.)))
if (has_rucc) da <- da %>% filter(!is.na(rucc_f))
il <- IQR(da[[ex]], na.rm = TRUE); mu <- mean(da[[ex]], na.rm = TRUE)
da$expo_z <- ifelse(da$lake_county == 1 & !is.na(da[[ex]]), (da[[ex]] - mu) / il, 0)
da <- da %>% filter(lake_county == 0 | !is.na(da[[ex]]))
cat(sprintf("Fitting all-county model (%d rows, %d lake-matched)...\n", nrow(da), sum(da$lake_county)))
fa <- fit_one(da, paste(c("expo_z", "lake_county", covs), collapse = " + "))
fits[["all_county"]] <- fa
res[[length(res) + 1]] <- pull(fa, "expo_z", "BYM2 ALL counties (+lake indicator)", "RR per IQR (among lake counties)")
res[[length(res) + 1]] <- pull(fa, "lake_county", "BYM2 ALL counties (+lake indicator)", "RR lake-matched vs. non-lake")

out <- do.call(rbind, res)
write.csv(out, "outputs/model_final/results_table.csv", row.names = FALSE)
saveRDS(fits, "outputs/model_final/fits.rds")
sink("outputs/model_final/summary.txt")
cat("FINAL MODEL (censored Poisson, BYM2), generated", format(Sys.time()), "\n")
cat("Primary exposure pre-specified 2026-10-02: px_share_20k per IQR\n\n")
print(out, row.names = FALSE, digits = 3)
sink()
print(out, row.names = FALSE, digits = 3)
cat("\nWrote outputs/model_final/results_table.csv and summary.txt\n")
