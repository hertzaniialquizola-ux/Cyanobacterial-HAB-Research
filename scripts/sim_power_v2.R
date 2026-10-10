# ============================================================================
# sim_power_v2.R -- corrected power / measurement-error simulation
#
# Fixes the flaw in the first simulation (fit_bym2_final_v3.R, section E): that version
# held every random effect at its fitted posterior mean, so the simulated data had only
# Poisson noise and the refit standard errors were about 2.6x too small (~0.005 vs the
# real fit's ~0.013). Here every replicate DRAWS new county, year and county-year random
# effects from the fitted variances, so the simulated data are as noisy as the real data.
#
#  * Calibration check: at true RR 1.00 and perfect exposure, the mean posterior SD of the
#    exposure coefficient should be close to the real non-spatial fit's SD (printed).
#  * Two exposure-error types: "classical" (independent each county-year) and
#    "persistent" (the same error every year within a county, like a lake-specific bias,
#    which averaging over years cannot remove).
#  * Writes a minimum-detectable-effect table (80% power) for each error scenario.
#  * Still a NON-spatial data-generating model and fit (for speed): treat power as an upper bound.
#
# Run from the project root (needs outputs/county_year_panel_v2.csv, scripts/fit_bym2_final_v3.R):
#   N_REP=4 CORES=2 Rscript scripts/sim_power_v2.R      # quick timing test
#   N_REP=100 CORES=4 Rscript scripts/sim_power_v2.R    # real run (20 scenarios x 100 fits)
# Writes outputs/model_final_v3/{sim_power_v2.csv, sim_mde.csv}
# ============================================================================
suppressWarnings({ library(sf); library(spdep); library(dplyr); library(INLA); library(parallel) })
src <- readLines("scripts/fit_bym2_final_v3.R")
k <- grep("^# ---- baseline", src)[1] - 1
eval(parse(text = src[1:k]))          # re-use the exact data prep and fit_one() from v3
inla.setOption(num.threads = "1:1")
N_REP <- as.integer(Sys.getenv("N_REP", "100")); CORES <- as.integer(Sys.getenv("CORES", "4"))
RNGkind("L'Ecuyer-CMRG"); set.seed(20261004)

cat("Fitting the non-spatial reference model on the real data...\n")
f_ns <- fit_one(d, rhs, spatial = FALSE, pred = TRUE)
b0 <- f_ns$summary.fixed["expo_z", "mean"]; se_real <- f_ns$summary.fixed["expo_z", "sd"]
n <- nrow(d)
lp <- f_ns$summary.linear.predictor$mean[1:n]
re_at <- function(nm, id) { r <- f_ns$summary.random[[nm]]; r$mean[match(id, r$ID)] }
re_hat <- re_at("idx", d$idx) + re_at("year_idx", d$year_idx) + re_at("obs_id", d$obs_id)
eta_fix <- lp - re_hat - b0 * d$expo_z            # fixed effects (+ offset) only, no random effects, no exposure
tot_obs <- sum(ifelse(d$deaths_suppressed, 4.5, d$y))
if (sum(exp(lp)) < tot_obs / 10) eta_fix <- eta_fix + d$log_pop    # offset was not in the INLA linear predictor
hp <- f_ns$summary.hyperpar
tau <- function(pat) hp$mean[grep(pat, rownames(hp))[1]]
tau_u <- tau("Precision for idx"); tau_v <- tau("Precision for year_idx"); tau_e <- tau("Precision for obs_id")
cat(sprintf("  real non-spatial exposure SD = %.4f (log scale); precisions: county %.1f, year %.1f, county-year %.1f\n",
            se_real, tau_u, tau_v, tau_e))
uc <- unique(d$idx); cty <- match(d$idx, uc); ny <- max(d$year_idx)
xt <- d$expo_z; sd_t <- sd(xt)
cat(sprintf("  check: observed ~%.0f deaths, simulated-truth expectation ~%.0f\n", tot_obs, sum(exp(eta_fix + sd(re_hat)^2 / 2))))

one_rep <- function(r, beta, lam, persist) {
  u <- rnorm(length(uc), 0, 1 / sqrt(tau_u))[cty]
  v <- rnorm(ny, 0, 1 / sqrt(tau_v))[d$year_idx]
  e <- rnorm(n, 0, 1 / sqrt(tau_e))
  yy <- rpois(n, exp(eta_fix + u + v + e + beta * xt))
  ys <- yy; ys[yy >= CENS[1] & yy <= CENS[2]] <- CENS[1]            # CDC suppression of 1-9
  sd_e <- sd_t * sqrt((1 - lam) / lam)                              # reliability = lam on the full exposure
  xo <- xt
  if (lam < 1) xo <- xt + sd_e * (sqrt(persist) * rnorm(length(uc))[cty] + sqrt(1 - persist) * rnorm(n))
  ds <- d; ds$y <- ys; ds$expo_z <- xo
  ft <- try(fit_one(ds, rhs, spatial = FALSE), silent = TRUE)
  if (inherits(ft, "try-error")) return(rep(NA_real_, 4))
  z <- ft$summary.fixed["expo_z", ]
  c(z$mean, z$sd, z[["0.025quant"]], z[["0.975quant"]])
}

scen <- expand.grid(RR = c(1.00, 1.05, 1.10, 1.20), rel = c(1.0, 0.5, 0.3), persist = c(0, 1))
scen <- scen[!(scen$rel == 1 & scen$persist == 1), ]
rows <- list()
for (s in seq_len(nrow(scen))) {
  beta <- log(scen$RR[s]); lam <- scen$rel[s]; pr <- scen$persist[s]
  out <- do.call(rbind, mclapply(seq_len(N_REP), one_rep, beta = beta, lam = lam, persist = pr,
                                 mc.cores = CORES, mc.set.seed = TRUE))
  out <- out[complete.cases(out), , drop = FALSE]
  rows[[s]] <- data.frame(true_RR = scen$RR[s], reliability = lam,
                          error_type = ifelse(lam == 1, "none", ifelse(pr == 1, "persistent", "classical")),
                          reps = nrow(out), median_est_RR = exp(median(out[, 1])),
                          mean_post_sd = mean(out[, 2]), empirical_sd = sd(out[, 1]),
                          power_lower_gt_1 = mean(out[, 3] > 0))
  cat(sprintf("  RR %.2f rel %.1f %-10s: median est %.3f, post SD %.4f (emp %.4f), power %.2f  [n=%d]\n",
              scen$RR[s], lam, rows[[s]]$error_type, rows[[s]]$median_est_RR, rows[[s]]$mean_post_sd,
              rows[[s]]$empirical_sd, rows[[s]]$power_lower_gt_1, nrow(out)))
  write.csv(do.call(rbind, rows), paste0(OUT3, "/sim_power_v2.csv"), row.names = FALSE)
}
sp <- do.call(rbind, rows)
cal <- sp[sp$true_RR == 1 & sp$reliability == 1, ]
cat(sprintf("\nCALIBRATION: simulated mean posterior SD %.4f vs real-data SD %.4f (should be close)\n", cal$mean_post_sd, se_real))

# minimum detectable true effect (per IQR) with 80% power, per error scenario
mde <- do.call(rbind, lapply(split(sp[sp$true_RR > 1, ], list(sp$reliability[sp$true_RR > 1], sp$error_type[sp$true_RR > 1]), drop = TRUE), function(g) {
  a <- median(log(g$median_est_RR) / log(g$true_RR))            # attenuation factor (observed log RR / true log RR)
  se <- mean(g$mean_post_sd)
  data.frame(reliability = g$reliability[1], error_type = g$error_type[1], attenuation = a,
             MDE_RR_80pct_power = exp((1.96 + 0.8416) * se / a))
}))
print(mde, row.names = FALSE, digits = 3)
write.csv(mde, paste0(OUT3, "/sim_mde.csv"), row.names = FALSE)
cat("\nWrote outputs/model_final_v3/sim_power_v2.csv and sim_mde.csv\n")
