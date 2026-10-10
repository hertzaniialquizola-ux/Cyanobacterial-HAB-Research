# ============================================================================
# make_bym2_spatial_map.R
# Maps the BYM2 county-level spatial random effect (residual relative risk after
# adjustment) from the final model. Uses saved fits -- no re-fitting, fast.
#
# Run from project root:  Rscript scripts/make_bym2_spatial_map.R
# Needs: sf, dplyr, ggplot2 (install.packages("ggplot2") if missing)
# Writes: outputs/figures/fig5_bym2_spatial_effect.png (+ _county_effects.csv)
# ============================================================================
suppressWarnings({ library(sf); library(dplyr); library(ggplot2); library(INLA) })

# Uses the FINAL model (run fit_bym2_final_v2.R first). Single fit -> posterior mean/sd directly.
NO_LST <- Sys.getenv("NO_LST", "0") == "1"
fit <- readRDS(if (NO_LST) "outputs/model_final_v2_nolst/fits.rds" else "outputs/model_final_v2/fits.rds")[["px_share_20k"]]
panel <- read.csv("outputs/county_year_panel_v2.csv", colClasses = c(county_fips = "character"))
cty <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
cty <- cty[!(cty$STATEFP %in% c("02","15","72","78","66","60","69")), ]
county_ids <- sort(intersect(unique(panel$county_fips), as.character(cty$GEOID)))
n <- length(county_ids)
map_idx <- data.frame(county_fips = county_ids, idx = seq_len(n)) %>%
  filter(county_fips %in% unique(panel$county_fips[!is.na(panel$px_share_20k)]))  # only counties with exposure data
pooled_mean <- fit$summary.random$idx$mean[1:n]
pooled_sd <- fit$summary.random$idx$sd[1:n]

eff <- data.frame(idx = 1:n, log_RR = pooled_mean, RR = exp(pooled_mean), sd = pooled_sd)
eff$sig <- ifelse(pooled_mean - 1.96 * pooled_sd > 0, "higher",
           ifelse(pooled_mean + 1.96 * pooled_sd < 0, "lower", "not different"))
eff <- inner_join(map_idx, eff, by = "idx")   # drops graph nodes with no data
write.csv(eff, "outputs/figures/fig5_bym2_county_effects.csv", row.names = FALSE)

counties <- st_read("data/raw/census/county_boundaries/cb_2022_us_county_20m.shp", quiet = TRUE)
counties <- counties[!(counties$STATEFP %in% c("02","15","72","78","66","60","69")), ]
counties$GEOID <- as.character(counties$GEOID)
counties <- st_transform(counties, 5070)
states <- counties %>% group_by(STATEFP) %>% summarise(.groups = "drop")
m <- left_join(counties, eff, by = c("GEOID" = "county_fips"))

lim <- max(as.numeric(quantile(abs(eff$log_RR), 0.98)), 0.1)
m$log_RR_c <- pmax(pmin(m$log_RR, lim), -lim)   # clip extremes for the colour scale
p <- ggplot() +
  geom_sf(data = m, aes(fill = log_RR_c), color = NA) +
  geom_sf(data = states, fill = NA, color = "grey35", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0,
                       limits = c(-lim, lim),
                       breaks = c(-lim, 0, lim), labels = function(x) sprintf("%.2f", exp(x)),
                       na.value = "grey90", name = "Spatial RR\n(white = 1.0)") +
  labs(title = "Residual spatial structure in kidney/liver mortality (BYM2)",
       subtitle = "Posterior mean county spatial effect after adjustment; grey = county not in analytic sample") +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.key.width = unit(1.6, "cm"),
        plot.title = element_text(face = "bold"), plot.margin = margin(8, 8, 8, 8))
ggsave("outputs/figures/fig5_bym2_spatial_effect.png", p, width = 8, height = 5.6, dpi = 400, bg = "white")

cat(sprintf("Counties mapped: %d. Credibly higher: %d, credibly lower: %d.\n",
            nrow(eff), sum(eff$sig == "higher"), sum(eff$sig == "lower")))
cat("Wrote outputs/figures/fig5_bym2_spatial_effect.png\n")
