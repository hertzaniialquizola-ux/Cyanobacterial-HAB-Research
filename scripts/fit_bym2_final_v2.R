# ============================================================================
# fit_bym2_final_v2.R  --  v2: adds overdispersion, AI/AN + PLACES covariates, lag/cumulative,
#   within/between, lead falsification, excl. 2020-21, zero-death rows retained (population fix)
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
# Run from project root:  Rscript scripts/fit_bym2_final_v2.R
# Writes outputs/model_final_v2/{results_table.csv, summary.txt, fits.rds}
# Runtime: ~8-10 BYM2 fits (each similar to one fit of the old MI loop).
# ============================================================================
suppressWarnings({ library(sf); library(spdep); library(dplyr); library(INLA) })
inla.setOption(num.threads = "1:1")
NO_LST <- Sys.getenv("NO_LST", "0") == "1"   # NO_LST=1: drop the (doubtful, sample-shrinking) MODIS lake-temperature covariate
OUT2 <- if (NO_LST) "outputs/model_final_v2_nolst" else "outputs/model_final_v2"
dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)
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
has_aian_col <- "pct_aian" %in% names(panel)
panel <- panel %>% filter(county_fips %in% county_ids, !is.na(population), population > 0) %>%
  inner_join(idx_map, by = "county_fips") %>%
  mutate(year_idx = year - min(year) + 1, log_pop = log(as.numeric(population)),
         y = ifelse(deaths_suppressed, CENS[1], deaths),
         ww10 = well_water_fraction * 10, sw10 = surface_supply_fraction * 10,
         pov10 = poverty_rate * 10, inc10k = median_household_income / 1e4, lst10 = lst_day_c / 10,
         age5 = pct_65plus * 20, blk10 = pct_black * 10, hsp10 = pct_hispanic * 10,
         aian10 = if (has_aian_col) pct_aian * 10 else NA_real_,
         lake_county = as.integer(!is.na(cei)))
if (has_rucc) panel$rucc_f <- factor(panel$rucc)
stopifnot(!any(is.na(panel$y)))
panel$obs_id <- seq_len(nrow(panel))
has_aian <- "aian10" %in% names(panel) && any(!is.na(panel$aian10))
has_places <- all(c("diabetes","obesity","binge","uninsured") %in% names(panel)) && any(!is.na(panel$diabetes))
if (has_places) { panel$diab <- panel$diabetes/10; panel$obes <- panel$obesity/10; panel$bing <- panel$binge/10; panel$unins <- panel$uninsured/10 }
cat(sprintf("Zero-death rows kept: %d | AI/AN available: %s | PLACES available: %s\n", sum(panel$deaths == 0, na.rm = TRUE), has_aian, has_places))

covs <- c("ww10", "sw10", "pov10", "inc10k", "age5", "blk10", "hsp10")
if (has_rucc) covs <- c(covs, "rucc_f")
if (has_aian) covs <- c(covs, "aian10")
cat(sprintf("Rows %d | suppressed (censored) %.1f%% | rurality covariate: %s\n", nrow(panel), 100 * mean(panel$deaths_suppressed), has_rucc))

pc <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
pc_obs <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
fit_one <- function(dat, rhs, spatial = TRUE, overdisp = TRUE) {
  sp <- if (spatial) "f(idx, model='bym2', graph=g, scale.model=TRUE, hyper=pc)" else "f(idx, model='iid')"
  fm <- as.formula(paste("y ~ 1 +", rhs, "+", sp, "+ f(year_idx, model='iid')",
                         if (overdisp) "+ f(obs_id, model='iid', hyper=pc_obs)" else ""))
  inla(fm, family = "cenpoisson", data = dat, offset = log_pop,
       control.family = list(cenpoisson.I = CENS),
       control.compute = list(dic = TRUE, waic = TRUE))
}
obs_sd <- function(fit) {
  h <- fit$marginals.hyperpar[["Precision for obs_id"]]
  if (is.null(h)) return(NA_real_)
  inla.emarginal(function(p) 1/sqrt(p), h)
}
pull <- function(fit, term, label, scale_note) {
  s <- fit$summary.fixed[term, ]; m <- fit$marginals.fixed[[term]]
  p0 <- inla.pmarginal(log(1.05), m) - inla.pmarginal(log(0.95), m)
  data.frame(model = label, term = term, unit = scale_note, RR = exp(s$mean),
             lo95 = exp(s[["0.025quant"]]), hi95 = exp(s[["0.975quant"]]),
             P_RR_in_0.95_1.05 = p0, obs_iid_sd = obs_sd(fit), DIC = fit$dic$dic, WAIC = fit$waic$waic, row.names = NULL)
}

# ---- 2. lake-matched analytic sample ---------------------------------------------------
lake <- panel %>% filter(lake_county == 1, (NO_LST | !is.na(lst_day_c)), !is.na(px_share_20k), !is.na(surface_supply_fraction),
                         if_all(c(pov10, inc10k, age5, blk10, hsp10, ww10), ~ !is.na(.)))
if (has_rucc) lake <- lake %>% filter(!is.na(rucc_f))
cat(sprintf("Lake-matched analytic sample: %d rows, %d counties\n", nrow(lake), n_distinct(lake$county_fips)))

covs_lake <- c(covs, if (!NO_LST) "lst10")
covs_places <- if (has_places) c(covs_lake, "diab", "obes", "bing", "unins") else NULL   # LST only exists for lake-matched rows
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
    if (!NO_LST) res[[length(res) + 1]] <- pull(f, "lst10", "BYM2 primary", "RR per 10 C (land-surface proxy; weak covariate)")
    # non-spatial comparator and interaction on the primary exposure
    fi <- fit_one(d, paste(c("expo_z", covs_lake), collapse = " + "), spatial = FALSE)
    fits[["primary_iid"]] <- fi
    res[[length(res) + 1]] <- pull(fi, "expo_z", "Non-spatial IID: primary", "RR per IQR")
    fx <- fit_one(d, paste(c("expo_z", "expo_z:sw10", covs_lake), collapse = " + "))
    fits[["primary_interaction"]] <- fx
    res[[length(res) + 1]] <- pull(fx, "expo_z:sw10", "BYM2 interaction: exposure x surface-supply (per 10 pp)", "ratio of RRs")
  }
}


# ---- 2b. v2 additions on the primary exposure ----------------------------------------------
d0 <- lake %>% filter(!is.na(px_share_20k))
iqr0 <- IQR(d0$px_share_20k); mu0 <- mean(d0$px_share_20k)
d0$expo_z <- (d0$px_share_20k - mu0) / iqr0
rhs0 <- paste(c("expo_z", covs_lake), collapse = " + ")
# (a) no overdispersion term (to show what it changes)
fits[["primary_noOD"]] <- fit_one(d0, rhs0, overdisp = FALSE)
res[[length(res)+1]] <- pull(fits[["primary_noOD"]], "expo_z", "BYM2 primary WITHOUT overdispersion term", "RR per IQR")
# (b) exclude 2020-2021 (COVID)
dc <- d0 %>% filter(!(year %in% c(2020, 2021)))
fits[["excl_covid"]] <- fit_one(dc, rhs0)
res[[length(res)+1]] <- pull(fits[["excl_covid"]], "expo_z", "BYM2 primary, excluding 2020-2021", "RR per IQR")
# (c) PLACES lifestyle/health covariates
if (has_places) {
  dp <- d0 %>% filter(if_all(c(diab, obes, bing, unins), ~ !is.na(.)))
  fits[["places"]] <- fit_one(dp, paste(c("expo_z", covs_places), collapse = " + "))
  res[[length(res)+1]] <- pull(fits[["places"]], "expo_z", "BYM2 primary + PLACES (diabetes, obesity, binge, uninsured)", "RR per IQR")
  for (tm in c("diab","obes","bing","unins")) res[[length(res)+1]] <- pull(fits[["places"]], tm, "BYM2 primary + PLACES", "RR per 10 pp")
}
# (d) lagged / cumulative exposure (uses same IQR scaling of the lagged variable)
for (lx in c("px_lag1", "px_cum3")) if (lx %in% names(lake)) {
  dl <- lake %>% filter(!is.na(.data[[lx]]), !is.na(px_share_20k))
  il <- IQR(dl[[lx]]); dl$expo_z <- (dl[[lx]] - mean(dl[[lx]])) / il
  fits[[lx]] <- fit_one(dl, paste(c("expo_z", covs_lake), collapse = " + "))
  res[[length(res)+1]] <- pull(fits[[lx]], "expo_z", paste0("BYM2 lagged exposure: ", lx), "RR per IQR")
}
# (e) lead falsification: next-year exposure should NOT predict this year's deaths once current exposure is in
if ("px_lead1" %in% names(lake)) {
  dd <- lake %>% filter(!is.na(px_lead1), !is.na(px_share_20k))
  dd$expo_z <- (dd$px_share_20k - mean(dd$px_share_20k)) / IQR(dd$px_share_20k)
  dd$lead_z <- (dd$px_lead1 - mean(dd$px_lead1)) / IQR(dd$px_lead1)
  fits[["lead"]] <- fit_one(dd, paste(c("expo_z", "lead_z", covs_lake), collapse = " + "))
  res[[length(res)+1]] <- pull(fits[["lead"]], "expo_z", "Falsification model: current exposure", "RR per IQR")
  res[[length(res)+1]] <- pull(fits[["lead"]], "lead_z", "Falsification model: NEXT-year exposure (should be ~1)", "RR per IQR")
}
# (f) within/between decomposition
if (all(c("px_cty_mean", "px_within") %in% names(lake))) {
  dw <- lake %>% filter(!is.na(px_cty_mean), !is.na(px_within))
  dw$btw_z <- (dw$px_cty_mean - mean(dw$px_cty_mean)) / IQR(dw$px_cty_mean)
  dw$wth_z <- dw$px_within / IQR(dw$px_within)
  fits[["within_between"]] <- fit_one(dw, paste(c("btw_z", "wth_z", covs_lake), collapse = " + "))
  res[[length(res)+1]] <- pull(fits[["within_between"]], "btw_z", "Within/between: BETWEEN-county (county mean)", "RR per IQR")
  res[[length(res)+1]] <- pull(fits[["within_between"]], "wth_z", "Within/between: WITHIN-county (year deviation)", "RR per IQR")
}
# (g) collinearity check poverty vs income
cat(sprintf("cor(poverty, income) in lake sample = %.2f\n", cor(d0$pov10, d0$inc10k, use = "complete.obs")))

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
write.csv(out, paste0(OUT2, "/results_table.csv"), row.names = FALSE)
saveRDS(fits, paste0(OUT2, "/fits.rds"))
sink(paste0(OUT2, "/summary.txt"))
cat("FINAL MODEL (censored Poisson, BYM2), generated", format(Sys.time()), "\n")
cat("Primary exposure designated 2026-10-02 (before results were seen, not formally registered): px_share_20k per IQR\n\n")
print(out, row.names = FALSE, digits = 3)
sink()
print(out, row.names = FALSE, digits = 3)
cat("\nWrote outputs/model_final_v2/results_table.csv and summary.txt\n")
