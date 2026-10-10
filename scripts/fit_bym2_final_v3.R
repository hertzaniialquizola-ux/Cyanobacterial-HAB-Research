# ============================================================================
# fit_bym2_final_v3.R -- answers the 2026-10-03 critique (no new data needed)
#
#  A. Region-by-year fixed effects and state-by-year random effects, re-testing the
#     FALSIFICATION (next-year exposure) model. If the 0.969 was residual
#     spatiotemporal confounding, it should shrink toward 1 here.
#  B. "Spatial+" (Dupont et al.): replace exposure by its residual after removing a
#     smooth spatial component, to check the BYM2 term is not distorting the estimate.
#  C. Full coefficient table for the primary model (race/ethnicity, AI/AN, RUCC...).
#  D. E-values (what unmeasured confounding could hide an effect).
#  E. OPTIONAL power / positive-control simulation (RUN_SIM=1): inject known effects
#     (RR 1.00/1.05/1.10/1.20 per IQR), degrade the exposure with classical
#     measurement error at several reliabilities, apply the CDC [1,9] censoring,
#     refit, and report detection probability. Uses a NON-spatial model for speed.
#
# Run from the project root:
#   Rscript scripts/fit_bym2_final_v3.R                 (A-D, ~10 fits)
#   RUN_SIM=1 N_SIM=30 Rscript scripts/fit_bym2_final_v3.R   (adds E; slow -- hours possible,
#                                                          start with N_SIM=5 to time it)
# Needs outputs/county_year_panel_v2.csv and outputs/model_mi/county_graph.adj (as before).
# Writes outputs/model_final_v3/{results_v3.csv, coef_primary.csv, evalues.csv, sim_power.csv}
# ============================================================================
suppressWarnings({ library(sf); library(spdep); library(dplyr); library(INLA) })
inla.setOption(num.threads = "1:1")
NO_LST <- Sys.getenv("NO_LST", "0") == "1"   # NO_LST=1: drop the (doubtful, sample-shrinking) MODIS lake-temperature covariate
OUT3 <- if (NO_LST) "outputs/model_final_v3_nolst" else "outputs/model_final_v3"
dir.create(OUT3, showWarnings = FALSE, recursive = TRUE)
CENS <- c(1, 9)
RUN_SIM <- Sys.getenv("RUN_SIM", "0") == "1"
N_SIM <- as.integer(Sys.getenv("N_SIM", "30"))

# ---- data (same construction as fit_bym2_final_v2.R) ------------------------------------
panel <- read.csv("outputs/county_year_panel_v2.csv", colClasses = c(county_fips = "character"))
panel$deaths_suppressed <- as.logical(panel$deaths_suppressed)
counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
counties <- counties[!(counties$STATEFP %in% c("02","15","72","78","66","60","69")), ]
county_ids <- sort(intersect(unique(panel$county_fips), as.character(counties$GEOID)))
g <- inla.read.graph("outputs/model_mi/county_graph.adj")
stopifnot(g$n == length(county_ids))
idx_map <- data.frame(county_fips = county_ids, idx = seq_along(county_ids))

ne <- c("09","23","25","33","44","50","34","36","42")
mw <- c("17","18","26","39","55","19","20","27","29","31","38","46")
so <- c("10","11","12","13","24","37","45","51","54","01","21","28","47","05","22","40","48")
region_of <- function(f) ifelse(f %in% ne, "NE", ifelse(f %in% mw, "MW", ifelse(f %in% so, "S", "W")))

panel <- panel %>% filter(county_fips %in% county_ids, !is.na(population), population > 0) %>%
  inner_join(idx_map, by = "county_fips") %>%
  mutate(year_idx = year - min(year) + 1, log_pop = log(as.numeric(population)),
         y = ifelse(deaths_suppressed, CENS[1], deaths),
         ww10 = well_water_fraction * 10, sw10 = surface_supply_fraction * 10,
         pov10 = poverty_rate * 10, inc10k = median_household_income / 1e4, lst10 = lst_day_c / 10,
         age5 = pct_65plus * 20, blk10 = pct_black * 10, hsp10 = pct_hispanic * 10,
         aian10 = pct_aian * 10, rucc_f = factor(rucc),
         lake_county = as.integer(!is.na(cei)),
         st = substr(county_fips, 1, 2), region = region_of(substr(county_fips, 1, 2)))
panel$obs_id <- seq_len(nrow(panel))
panel$region_year <- factor(paste(panel$region, panel$year))
panel$state_year_idx <- as.integer(factor(paste(panel$st, panel$year)))
covs <- c("ww10", "sw10", "pov10", "inc10k", "age5", "blk10", "hsp10", "aian10", "rucc_f", if (!NO_LST) "lst10")

lake <- panel %>% filter(lake_county == 1, (NO_LST | !is.na(lst_day_c)), !is.na(px_share_20k), !is.na(surface_supply_fraction),
                         if_all(c(pov10, inc10k, age5, blk10, hsp10, ww10, aian10), ~ !is.na(.)), !is.na(rucc_f))
cat(sprintf("Lake-matched analytic sample: %d rows, %d counties (v2 had 5498 / 845)\n", nrow(lake), n_distinct(lake$county_fips)))

pc <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
pc_obs <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
fit_one <- function(dat, rhs, spatial = TRUE, extra = "", pred = FALSE) {
  sp <- if (spatial) "f(idx, model='bym2', graph=g, scale.model=TRUE, hyper=pc)" else "f(idx, model='iid')"
  fm <- as.formula(paste("y ~ 1 +", rhs, "+", sp, "+ f(year_idx, model='iid') + f(obs_id, model='iid', hyper=pc_obs)", extra))
  inla(fm, family = "cenpoisson", data = dat, offset = log_pop, control.family = list(cenpoisson.I = CENS),
       control.compute = list(dic = TRUE, waic = TRUE, config = FALSE), control.predictor = list(compute = pred))
}
pull <- function(fit, term, label) {
  s <- fit$summary.fixed[term, ]
  data.frame(model = label, term = term, RR = exp(s$mean), lo95 = exp(s[["0.025quant"]]), hi95 = exp(s[["0.975quant"]]),
             DIC = fit$dic$dic, row.names = NULL)
}
std <- function(x) (x - mean(x)) / IQR(x)

d <- lake
d$expo_z <- std(d$px_share_20k)
d$lead_z <- NA_real_
dl <- d %>% filter(!is.na(px_lead1)); dl$expo_z <- std(dl$px_share_20k); dl$lead_z <- std(dl$px_lead1)
rhs <- paste(c("expo_z", covs), collapse = " + ")
rhs_lead <- paste(c("expo_z", "lead_z", covs), collapse = " + ")
res <- list()

# ---- baseline (refit so every comparison uses the identical specification) -------------
cat("Baseline fits...\n")
f0 <- fit_one(d, rhs);        res[[length(res)+1]] <- pull(f0, "expo_z", "Baseline BYM2 primary")
f0l <- fit_one(dl, rhs_lead); res[[length(res)+1]] <- pull(f0l, "expo_z", "Baseline falsification: current")
res[[length(res)+1]] <- pull(f0l, "lead_z", "Baseline falsification: NEXT-year")

# ---- A. region x year fixed effects ----------------------------------------------------
cat("A. region-by-year fixed effects...\n")
fa <- fit_one(d, paste(rhs, "+ region_year"));            res[[length(res)+1]] <- pull(fa, "expo_z", "Region-by-year FE: primary")
fal <- fit_one(dl, paste(rhs_lead, "+ region_year"))
res[[length(res)+1]] <- pull(fal, "expo_z", "Region-by-year FE: falsification current")
res[[length(res)+1]] <- pull(fal, "lead_z", "Region-by-year FE: falsification NEXT-year")
# ---- A2. state x year iid random effects -----------------------------------------------
cat("A2. state-by-year random effects...\n")
ex_sy <- "+ f(state_year_idx, model='iid', hyper=pc_obs)"
fs <- fit_one(d, rhs, extra = ex_sy);        res[[length(res)+1]] <- pull(fs, "expo_z", "State-by-year RE: primary")
fsl <- fit_one(dl, rhs_lead, extra = ex_sy)
res[[length(res)+1]] <- pull(fsl, "expo_z", "State-by-year RE: falsification current")
res[[length(res)+1]] <- pull(fsl, "lead_z", "State-by-year RE: falsification NEXT-year")

# ---- B. Spatial+ (corrected) ------------------------------------------------------------
# Remove only the SMOOTH (structured, ICAR/besag) spatial component of the county-mean exposure.
# (The first version used a BYM2 term, whose unstructured part can absorb ALL variation when there is
#  one observation per county; that gave a meaningless "99.7% spatial" and a near-zero residual.)
cat("B. Spatial+ ...\n")
cm <- d %>% group_by(idx) %>% summarise(x = mean(px_share_20k), .groups = "drop")
cm_full <- data.frame(idx = seq_len(g$n)) %>% left_join(cm, by = "idx")     # counties without data = NA (smoothed over)
fx <- inla(x ~ 1 + f(idx, model = "besag", graph = g, scale.model = TRUE,
                     hyper = list(prec = list(prior = "pc.prec", param = c(0.3, 0.01)))),
           data = cm_full, family = "gaussian")
cm_full$x_spatial <- fx$summary.random$idx$mean[match(cm_full$idx, fx$summary.random$idx$ID)]
sp_var <- var(cm_full$x_spatial[!is.na(cm_full$x)]); tot_var <- var(cm_full$x, na.rm = TRUE)
cat(sprintf("  share of between-county exposure variance in the smooth spatial component: %.1f%%\n", 100 * sp_var / tot_var))
dsp <- d %>% left_join(cm_full[, c("idx", "x_spatial")], by = "idx")
dsp$expo_sp <- dsp$px_share_20k - dsp$x_spatial
cat(sprintf("  IQR of residual exposure = %.3f (original IQR = %.3f)\n", IQR(dsp$expo_sp), IQR(d$px_share_20k)))
dsp$expo_z <- dsp$expo_sp / IQR(d$px_share_20k)            # scale by the ORIGINAL IQR so the RR is comparable
fp <- fit_one(dsp, rhs); res[[length(res)+1]] <- pull(fp, "expo_z", "Spatial+: residual exposure (per ORIGINAL IQR)")

out <- do.call(rbind, res)
write.csv(out, paste0(OUT3, "/results_v3.csv"), row.names = FALSE)
print(out, row.names = FALSE, digits = 3)

# ---- C. full coefficient table for the primary model ---------------------------------
cf <- f0$summary.fixed
coefs <- data.frame(term = rownames(cf), RR = exp(cf$mean), lo95 = exp(cf[["0.025quant"]]), hi95 = exp(cf[["0.975quant"]]))
write.csv(coefs, paste0(OUT3, "/coef_primary.csv"), row.names = FALSE)
cat("\nPRIMARY MODEL, all coefficients (RR):\n"); print(coefs, row.names = FALSE, digits = 3)
hp <- f0$summary.hyperpar
cat("\nHyperparameters:\n"); print(hp[, c("mean", "0.025quant", "0.975quant")], digits = 3)

# ---- D. E-values ---------------------------------------------------------------------
evalue <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ev <- data.frame(scenario = c("true RR 1.05 per IQR hidden by confounding (observed ~0.993)", "true RR 1.10 hidden", "true RR 1.20 hidden"),
                 ratio_needed = c(1.05 / 0.993, 1.10 / 0.993, 1.20 / 0.993))
ev$E_value <- evalue(ev$ratio_needed)
write.csv(ev, paste0(OUT3, "/evalues.csv"), row.names = FALSE)
cat("\nE-values (confounder-exposure and confounder-outcome RR both >= E needed):\n"); print(ev, row.names = FALSE, digits = 3)

# ---- E. power / positive-control simulation ------------------------------------------
if (RUN_SIM) {
  cat(sprintf("\nSIMULATION: N_SIM=%d per scenario (non-spatial model for speed)\n", N_SIM))
  set.seed(20261003)
  f_ns <- fit_one(d, rhs, spatial = FALSE, pred = TRUE)
  b0 <- f_ns$summary.fixed["expo_z", "mean"]
  eta <- f_ns$summary.linear.predictor$mean[1:nrow(d)] - b0 * d$expo_z   # linear predictor WITHOUT the exposure term
  tot_obs <- sum(ifelse(d$deaths_suppressed, 4.5, d$y)); tot_pred <- sum(exp(eta + b0 * d$expo_z))
  if (tot_pred < tot_obs / 10) eta <- eta + d$log_pop          # offset was not included in the INLA linear predictor
  cat(sprintf("  check: observed deaths ~%.0f, simulated-truth expected deaths %.0f\n", tot_obs, sum(exp(eta + b0 * d$expo_z))))
  xt <- d$expo_z                                                         # treated as the "true" exposure, in IQR units
  sd_t <- sd(xt)
  scen <- expand.grid(RR = c(1.00, 1.05, 1.10, 1.20), reliability = c(1.0, 0.8, 0.5, 0.3))
  rows <- list()
  for (s in seq_len(nrow(scen))) {
    beta <- log(scen$RR[s]); lam <- scen$reliability[s]; est <- lo <- hi <- rep(NA_real_, N_SIM)
    for (r in seq_len(N_SIM)) {
      mu <- exp(eta + beta * xt)
      yy <- rpois(nrow(d), mu)
      cens <- yy >= CENS[1] & yy <= CENS[2]; ys <- yy; ys[cens] <- CENS[1]
      xo <- xt + rnorm(nrow(d), 0, sd_t * sqrt((1 - lam) / lam))        # classical error: reliability = lam
      ds <- d; ds$y <- ys; ds$expo_z <- xo
      ft <- try(fit_one(ds, rhs, spatial = FALSE), silent = TRUE)
      if (inherits(ft, "try-error")) next
      z <- ft$summary.fixed["expo_z", ]
      est[r] <- z$mean; lo[r] <- z[["0.025quant"]]; hi[r] <- z[["0.975quant"]]
    }
    ok <- !is.na(est)
    rows[[s]] <- data.frame(true_RR = scen$RR[s], reliability = lam, reps = sum(ok),
                            median_est_RR = exp(median(est[ok])),
                            P_detect_harm = mean(lo[ok] > 0), P_detect_any = mean(lo[ok] > 0 | hi[ok] < 0),
                            P_ci_upper_below_1.05 = mean(exp(hi[ok]) < 1.05))
    cat(sprintf("  RR %.2f, reliability %.1f: median est %.3f, detect(harm) %.2f\n", scen$RR[s], lam, rows[[s]]$median_est_RR, rows[[s]]$P_detect_harm))
    write.csv(do.call(rbind, rows), paste0(OUT3, "/sim_power.csv"), row.names = FALSE)
  }
}
cat("\nWrote outputs/model_final_v3/*\n")
