# ============================================================================
# fit_dml_pooled.R -- "Hybrid B": double/debiased machine learning (DML) on the pooled
# 2018-2024 age-standardised county data, as a SENSITIVITY check on the pooled BYM2 results.
#
# What it does (partially linear model, Chernozhukov et al. 2018):
#   Y = log(observed deaths / expected deaths)   (indirectly age-standardised log SMR)
#   D = exposure per IQR (same scaling as the pooled BYM2 models)
#   X = the same covariates as fit_pooled_age_std.R (well, surface supply, poverty, income,
#       % Black, % Hispanic, % AIAN, rurality), optionally plus county-centroid x,y.
#   Random forests predict Y and D from X on held-out data (cross-fitting). The effect is the
#   weighted regression of the Y-residual on the D-residual. Folds are SPATIAL BLOCKS
#   (k-means clusters of county centroids), so neighbouring counties never train and test
#   against each other; standard errors are cluster-robust by block. Repeated N_REP times with
#   different fold splits; the median is reported.
#
# LIMITS (say these in the paper): uses only counties with UNsuppressed death counts (the
# interval-censoring cannot be handled), so the sample is smaller and favours larger counties;
# counties with 0 deaths are dropped (log 0); it is a different estimator and sample from the
# BYM2 model, so it is a sensitivity check, not a replacement. The linear WLS benchmark on the
# SAME sample is printed so the DML-vs-linear difference is not confused with the sample change.
#
# Run from the project root. Needs the R package ranger (install.packages("ranger")).
#   Rscript scripts/fit_dml_pooled.R
#   N_REP=10 TAU2=0.02 K_BLOCKS=20 Rscript scripts/fit_dml_pooled.R
# Writes outputs/model_pooled_nolst/dml_results.csv
# ============================================================================
suppressWarnings({ library(sf); library(dplyr) })
if (!requireNamespace("ranger", quietly = TRUE)) stop('Install first:  install.packages("ranger")')
library(ranger)
N_REP    <- as.integer(Sys.getenv("N_REP", "5"))
K_BLOCKS <- as.integer(Sys.getenv("K_BLOCKS", "20"))
N_FOLDS  <- as.integer(Sys.getenv("N_FOLDS", "5"))
MIN_NODE <- as.integer(Sys.getenv("MIN_NODE", "5"))      # larger = smoother forests (try 20-40 if r2_exposure_from_X is negative)
TAU2     <- as.numeric(Sys.getenv("TAU2", "0.02"))   # between-county variance of log SMR used in the weights
OUTP <- "outputs/model_pooled_nolst"; dir.create(OUTP, showWarnings = FALSE, recursive = TRUE)

# ---- same county universe and filters as fit_pooled_age_std.R (NO_LST=1 version) -----------
dat <- read.csv("outputs/county_pooled_panel.csv", colClasses = c(county_fips = "character"))
counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
counties <- counties[!(counties$STATEFP %in% c("02","15","72","78","66","60","69")), ]
county_ids <- sort(intersect(unique(read.csv("outputs/county_year_panel_v2.csv", colClasses = c(county_fips = "character"))$county_fips),
                             as.character(counties$GEOID)))
dat <- dat %>% filter(county_fips %in% county_ids) %>%
  mutate(ww10 = well * 10, sw10 = surf * 10, pov10 = pov * 10, inc10k = inc / 1e4, blk10 = blk * 10,
         hsp10 = hsp * 10, aian10 = aian * 10, rucc_f = factor(round(rucc)))
lake <- dat %>% filter(!is.na(px_mean), n_years_exposure >= 3, !is.na(rucc_f),
                       if_all(c(ww10, sw10, pov10, inc10k, blk10, hsp10, aian10), ~ !is.na(.)))
lake$D <- (lake$px_mean - mean(lake$px_mean)) / IQR(lake$px_mean)     # per IQR, as in the BYM2 runs
cat(sprintf("Pooled lake-matched counties: %d (should match the BYM2 negative-control n = 978)\n", nrow(lake)))

# ---- county centroids and spatial blocks -------------------------------------------------------
cen <- st_centroid(st_geometry(st_transform(counties, 5070)))
xy <- st_coordinates(cen)
coords <- data.frame(county_fips = as.character(counties$GEOID), x = xy[, 1] / 1e5, y = xy[, 2] / 1e5)
lake <- left_join(lake, coords, by = "county_fips")
set.seed(20261009)
lake$block <- kmeans(lake[, c("x", "y")], centers = K_BLOCKS, nstart = 10)$cluster
covs <- c("ww10", "sw10", "pov10", "inc10k", "blk10", "hsp10", "aian10", "rucc_f")

# ---- estimators ---------------------------------------------------------------------------------
run_dml <- function(d, xcols) {
  X <- d[, xcols, drop = FALSE]
  est <- se2 <- r2 <- numeric(N_REP)
  for (r in seq_len(N_REP)) {
    set.seed(1000 + r)
    fold <- sample(rep_len(seq_len(N_FOLDS), K_BLOCKS))[d$block]
    gY <- gD <- numeric(nrow(d))
    for (k in seq_len(N_FOLDS)) {
      tr <- fold != k; te <- !tr
      if (!any(te) || !any(tr)) next
      fy <- ranger(x = X[tr, , drop = FALSE], y = d$Y[tr], num.trees = 500, min.node.size = MIN_NODE, case.weights = d$w[tr])
      fd <- ranger(x = X[tr, , drop = FALSE], y = d$D[tr], num.trees = 500, min.node.size = MIN_NODE, case.weights = d$w[tr])
      gY[te] <- predict(fy, data = X[te, , drop = FALSE])$predictions
      gD[te] <- predict(fd, data = X[te, , drop = FALSE])$predictions
    }
    rY <- d$Y - gY; rD <- d$D - gD
    den <- sum(d$w * rD^2); th <- sum(d$w * rD * rY) / den
    psi <- d$w * rD * (rY - th * rD)
    bs <- tapply(psi, d$block, sum); B <- length(bs)
    est[r] <- th; se2[r] <- B / (B - 1) * sum(bs^2) / den^2
    r2[r] <- 1 - sum(d$w * rD^2) / sum(d$w * (d$D - weighted.mean(d$D, d$w))^2)
  }
  th <- median(est)
  c(theta = th, se = sqrt(median(se2 + (est - th)^2)), r2_exposure_from_X = median(r2))
}
run_wls <- function(d) {
  fit <- lm(as.formula(paste("Y ~ D +", paste(covs, collapse = " + "))), data = d, weights = w)
  Xm <- model.matrix(fit); e <- residuals(fit); W <- d$w
  bread <- solve(t(Xm) %*% (W * Xm))
  sc <- rowsum(Xm * (W * e), d$block)
  B <- nrow(sc); n <- nrow(Xm); p <- ncol(Xm)
  V <- bread %*% (t(sc) %*% sc) %*% bread * (B / (B - 1)) * ((n - 1) / (n - p))
  c(theta = unname(coef(fit)["D"]), se = sqrt(V["D", "D"]), r2_exposure_from_X = NA_real_)
}

res <- list()
for (oc in c("kidney", "liver", "negctrl")) {
  dcol <- paste0(oc, "_deaths"); scol <- paste0(oc, "_suppressed"); ecol <- paste0(oc, "_expected")
  if (!(dcol %in% names(lake))) { cat(sprintf("[skip] %s\n", oc)); next }
  d <- lake
  d[[scol]] <- as.logical(d[[scol]])
  d <- d[!is.na(d[[ecol]]) & d[[ecol]] > 0 & !is.na(d[[dcol]]) & !d[[scol]] & d[[dcol]] > 0, ]
  d$Y <- log(d[[dcol]] / d[[ecol]])
  d$w <- 1 / (1 / d[[dcol]] + TAU2)
  cat(sprintf("\n%s: %d uncensored counties (of %d lake-matched), %d spatial blocks\n", oc, nrow(d), nrow(lake), length(unique(d$block))))
  for (m in c("DML, covariates only", "DML, covariates + coordinates", "Linear WLS benchmark (same sample)")) {
    z <- switch(m, "DML, covariates only" = run_dml(d, covs),
                   "DML, covariates + coordinates" = run_dml(d, c(covs, "x", "y")),
                   run_wls(d))
    res[[length(res) + 1]] <- data.frame(outcome = oc, method = m, n_counties = nrow(d),
      RR_per_IQR = exp(z[["theta"]]), lo95 = exp(z[["theta"]] - 1.96 * z[["se"]]), hi95 = exp(z[["theta"]] + 1.96 * z[["se"]]),
      r2_exposure_from_X = z[["r2_exposure_from_X"]])
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUTP, "dml_results.csv"), row.names = FALSE)
cat("\n=== DML sensitivity check (per IQR of exposure) ===\n"); print(out, row.names = FALSE, digits = 3)
ref <- file.path(OUTP, "results_pooled.csv")
if (file.exists(ref)) { cat("\nReference: pooled BYM2 (censored Poisson) results\n"); print(read.csv(ref), row.names = FALSE, digits = 3) }
cat("\nWrote", file.path(OUTP, "dml_results.csv"), "\n")
