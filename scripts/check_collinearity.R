# ============================================================================
# check_collinearity.R
#
# Quick pairwise-correlation + VIF (variance inflation factor) check between
# the fixed effects used in fit_bym2_model_mi.R -- specifically to answer:
# is poverty_rate's persistent, counter-intuitive-direction effect (section
# 14 of the project doc, RR 0.312 even after the lake-adjacent sample fix) a
# collinearity artifact with median_household_income, or a real, separate
# effect?
#
# This does NOT require INLA or spatial packages -- it's a plain-R check on
# the same analytic sample (lake-adjacent counties, i.e. filter(!is.na(cei)))
# used by both BYM2 scripts, so the numbers here are directly relevant to
# those models' fixed effects. VIF is computed manually via auxiliary
# regressions (VIF_j = 1 / (1 - R^2_j), where R^2_j comes from regressing
# predictor j on all the other predictors) rather than pulling in the `car`
# package, so this has no new dependencies beyond what's already installed.
#
# Rule of thumb: VIF > 5 is worth a closer look, VIF > 10 is a real problem
# (the coefficient's standard error is inflated enough that its estimate --
# including its sign -- becomes unstable and hard to trust).
#
# Run from the project root:
#   cd ~/Projects/Cyanobacterial-HAB-Research
#   Rscript scripts/check_collinearity.R
#
# Writes: outputs/model_mi/collinearity_check.txt
# ============================================================================

suppressWarnings(library(dplyr))

panel <- read.csv("outputs/county_year_panel_full.csv", colClasses = c(county_fips = "character"))

required_cols <- c("poverty_rate", "median_household_income", "cei", "well_water_fraction", "population")
missing_cols <- setdiff(required_cols, names(panel))
if (length(missing_cols) > 0) {
  stop("county_year_panel_full.csv is missing: ", paste(missing_cols, collapse = ", "))
}

# Same scope as both BYM2 scripts' analytic sample (section 14 of the project
# doc): population > 0, and restricted to lake-adjacent counties (a real
# CyAN match, cei not NA). Deliberately NOT restricted to non-suppressed rows
# -- poverty_rate/median_household_income/cei/well_water_fraction are all
# covariates, present regardless of whether that row's death count happens
# to be suppressed, so using the full 7,602-row analytic frame (not just the
# 3,945 non-suppressed rows) gives the most statistically stable read on
# collinearity among the predictors themselves.
df <- panel %>%
  mutate(population = as.numeric(population), cei = as.numeric(cei)) %>%
  filter(!is.na(population), population > 0, !is.na(cei)) %>%
  mutate(median_household_income_10k = median_household_income / 10000)

predictors <- c("cei", "well_water_fraction", "poverty_rate", "median_household_income_10k")
df_check <- df[, predictors]
n_before <- nrow(df_check)
df_check <- df_check[complete.cases(df_check), ]
n_after <- nrow(df_check)

cat(sprintf("Rows in lake-adjacent analytic sample: %d (%d dropped for missing covariates)\n\n", n_after, n_before - n_after))

# --- Pairwise correlation matrix ---
cor_mat <- cor(df_check, use = "complete.obs")

# --- VIF via auxiliary regressions (no `car` dependency) ---
vif_manual <- function(data, predictors) {
  out <- data.frame(term = predictors, VIF = NA_real_, aux_R2 = NA_real_)
  for (i in seq_along(predictors)) {
    y <- predictors[i]
    x <- setdiff(predictors, y)
    f <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
    r2 <- summary(lm(f, data = data))$r.squared
    out$aux_R2[i] <- r2
    out$VIF[i] <- 1 / (1 - r2)
  }
  out
}
vif_tab <- vif_manual(df_check, predictors)

flag <- function(v) if (v > 10) "PROBLEM (>10)" else if (v > 5) "worth a look (>5)" else "fine"
vif_tab$flag <- sapply(vif_tab$VIF, flag)

out_lines <- c(
  "=== Collinearity check: fixed effects in fit_bym2_model_mi.R ===",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Sample: lake-adjacent counties only (filter(!is.na(cei))), N = %d", n_after),
  "",
  "--- Pairwise correlation matrix ---",
  capture.output(print(round(cor_mat, 3))),
  "",
  sprintf("--- poverty_rate vs. median_household_income_10k: r = %.3f ---", cor_mat["poverty_rate", "median_household_income_10k"]),
  "",
  "--- Variance Inflation Factors (VIF_j = 1 / (1 - R^2_j)) ---",
  capture.output(print(vif_tab, row.names = FALSE)),
  "",
  "Interpretation: VIF > 5 means this predictor's estimated coefficient (and",
  "its standard error / significance) is meaningfully entangled with the",
  "other predictors in the same model. VIF > 10 means it's unstable enough",
  "that the sign/significance shouldn't be trusted without further checks",
  "(e.g. dropping one of the correlated pair and seeing if the other's",
  "estimate holds up)."
)

dir.create("outputs/model_mi", showWarnings = FALSE, recursive = TRUE)
writeLines(out_lines, "outputs/model_mi/collinearity_check.txt")
cat(paste(out_lines, collapse = "\n"))
cat("\n\nWritten to outputs/model_mi/collinearity_check.txt\n")
