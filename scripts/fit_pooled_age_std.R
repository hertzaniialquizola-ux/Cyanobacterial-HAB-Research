# ============================================================================
# fit_pooled_age_std.R -- 7-year pooled, INDIRECTLY AGE-STANDARDIZED county models
# Run after: get_acs_age_groups.py, WONDER exports, build_pooled_panel.py
#   Rscript scripts/fit_pooled_age_std.R
# One row per county. Outcome = pooled 2018-2024 deaths (CDC-suppressed 1-9 -> censored interval),
# offset = log(expected deaths from national age-specific rates x county ACS age structure).
# Exposure = county mean (2018-2024) of the primary exposure, per IQR. Outcomes modeled SEPARATELY:
# kidney (N18), liver (K74+K76), and (if exported) K70 and a negative-control cause.
# Writes outputs/model_pooled/{results_pooled.csv, coef_<outcome>.csv}
# ============================================================================
suppressWarnings({ library(sf); library(spdep); library(dplyr); library(INLA) })
inla.setOption(num.threads = "1:1")
NO_LST <- Sys.getenv("NO_LST", "0") == "1"   # NO_LST=1: drop the (doubtful, sample-shrinking) MODIS lake-temperature covariate
OUTP <- if (NO_LST) "outputs/model_pooled_nolst" else "outputs/model_pooled"
dir.create(OUTP, showWarnings = FALSE, recursive = TRUE)
CENS <- c(1, 9)
dat <- read.csv("outputs/county_pooled_panel.csv", colClasses = c(county_fips = "character"))
counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
counties <- counties[!(counties$STATEFP %in% c("02","15","72","78","66","60","69")), ]
county_ids <- sort(intersect(unique(read.csv("outputs/county_year_panel_v2.csv", colClasses = c(county_fips = "character"))$county_fips),
                             as.character(counties$GEOID)))
g <- inla.read.graph("outputs/model_mi/county_graph.adj")
stopifnot(g$n == length(county_ids))      # same node order as the annual models
idx_map <- data.frame(county_fips = county_ids, idx = seq_along(county_ids))
dat <- dat %>% inner_join(idx_map, by = "county_fips") %>%
  mutate(ww10 = well * 10, sw10 = surf * 10, pov10 = pov * 10, inc10k = inc / 1e4, blk10 = blk * 10,
         hsp10 = hsp * 10, aian10 = aian * 10, rucc_f = factor(round(rucc)), lst10 = lst / 10)
lake <- dat %>% filter(!is.na(px_mean), n_years_exposure >= 3, (NO_LST | !is.na(lst10)), !is.na(rucc_f),
                       if_all(c(ww10, sw10, pov10, inc10k, blk10, hsp10, aian10), ~ !is.na(.)))
cat(sprintf("Pooled lake-matched counties: %d\n", nrow(lake)))
lake$expo_z <- (lake$px_mean - mean(lake$px_mean)) / IQR(lake$px_mean)
covs <- c("ww10", "sw10", "pov10", "inc10k", "blk10", "hsp10", "aian10", "rucc_f", if (!NO_LST) "lst10")   # NO % aged 65+: age is in the offset
pc <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
res <- list()
for (oc in c("kidney", "liver", "k70", "negctrl")) {
  dcol <- paste0(oc, "_deaths"); scol <- paste0(oc, "_suppressed"); ecol <- paste0(oc, "_expected")
  if (!(dcol %in% names(lake))) { cat(sprintf("[skip] %s (not in pooled panel)\n", oc)); next }
  lake[[scol]] <- as.logical(lake[[scol]])
  d <- lake %>% filter(!is.na(.data[[ecol]]), .data[[ecol]] > 0, !is.na(.data[[dcol]]) | .data[[scol]])
  d$y <- ifelse(d[[scol]], CENS[1], d[[dcol]]); d$off <- log(d[[ecol]])
  cat(sprintf("%s: %d counties, %.1f%% censored, total observed ~%.0f vs expected %.0f\n", oc, nrow(d),
              100 * mean(d[[scol]]), sum(ifelse(d[[scol]], 4.5, d$y)), sum(d[[ecol]])))
  fm <- as.formula(paste("y ~ 1 + expo_z +", paste(covs, collapse = " + "), "+ f(idx, model='bym2', graph=g, scale.model=TRUE, hyper=pc)"))
  f <- inla(fm, family = "cenpoisson", data = d, offset = off, control.family = list(cenpoisson.I = CENS),
            control.compute = list(dic = TRUE))
  s <- f$summary.fixed["expo_z", ]
  res[[length(res) + 1]] <- data.frame(outcome = oc, n_counties = nrow(d), pct_censored = 100 * mean(d[[scol]]),
        RR_per_IQR = exp(s$mean), lo95 = exp(s[["0.025quant"]]), hi95 = exp(s[["0.975quant"]]))
  cf <- f$summary.fixed
  write.csv(data.frame(term = rownames(cf), RR = exp(cf$mean), lo95 = exp(cf[["0.025quant"]]), hi95 = exp(cf[["0.975quant"]])),
            sprintf(paste0(OUTP, "/coef_%s.csv"), oc), row.names = FALSE)
}
out <- do.call(rbind, res)
write.csv(out, paste0(OUTP, "/results_pooled.csv"), row.names = FALSE)
print(out, row.names = FALSE, digits = 3)
cat("\nWrote outputs/model_pooled/\n")
