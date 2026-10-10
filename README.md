# Satellite-measured cyanobacterial bloom intensity and kidney and liver disease mortality across U.S. counties

Analysis code and derived county-year tables for a county-level Bayesian spatial study (2018-2024) of
whether satellite-measured cyanobacterial bloom intensity in a county's lakes is associated with kidney
(ICD-10 N18) and liver (K74, K76) mortality. Data are public; this repository holds code and small derived
tables only (see [data/README.md](data/README.md) for how to obtain the raw inputs).

## Main results (final sample: 6,675 county-years, 1,015 counties)

| Analysis | Rate ratio per IQR of exposure (95% interval) |
|---|---|
| Primary: share of lake pixels >= 20,000 cells/mL (censored-Poisson BYM2) | 0.997 (0.977-1.018) |
| Falsification: next-year exposure (should be null) | 0.981 (0.957-1.006) |
| Pooled 2018-2024, age-standardized kidney | 1.007 (0.969-1.047) |
| Pooled 2018-2024, age-standardized liver | 1.001 (0.969-1.033) |
| Negative-control outcome (accidental deaths) | 0.989 (0.967-1.010) |

No specification showed increased mortality. This is a **null result at county level with a noisy exposure**, not
evidence of safety: the simulation shows the design can detect a 3% increase only if exposure were measured without
error, and error can hide much larger effects (`outputs/model_final_v3_nolst/sensitivity_reliability.csv`).
Satellite and measured microcystin agree moderately lake by lake but weakly at county level (EPA NLA 2022).

## Methods in brief
- **Exposure:** CyAN CI_cyano 300 m weekly composites decoded to cells/mL; per lake-in-county piece, share of valid pixels >= 20,000 cells/mL; annual mean over valid weeks (`build_alt_cei_v2.py`).
- **Outcome:** CDC WONDER county deaths. Suppressed counts (1-9) are treated as interval-censored in a censored-Poisson likelihood (no imputation).
- **Model:** R-INLA BYM2 spatial model with year and county-year random effects, adjusted for age, race/ethnicity, poverty, income, well-water and surface-supply shares, and rurality.
- **Checks:** next-year-exposure falsification, region- and state-by-year adjustment, Spatial+, E-values, negative-control outcome, power and measurement-error simulation, validation against EPA NLA 2022 microcystin, double machine learning on the pooled data.

## Repository layout
```
scripts/     all analysis code (see table below)
outputs/     small derived tables: county_year_panel_v2.csv, county_pooled_panel.csv, model result tables,
             validation tables, county adjacency graph (model_mi/county_graph.adj), figure script
data/        README only; raw data are not committed
```

## Reproducing the results
1. **Environment.** Python 3 (`pip install -r requirements.txt`) and R with the packages in `r-packages.txt` (R-INLA installs from its own repository). Run every command from the project root.
2. **Fast path (recommended):** the derived panels are included, so the models can be refit without any raw data:
   ```
   NO_LST=1 Rscript scripts/fit_bym2_final_v3.R        # primary, robustness, falsification, E-values
   NO_LST=1 Rscript scripts/fit_pooled_age_std.R       # pooled kidney, liver, negative control
   Rscript scripts/fit_dml_pooled.R                     # double machine learning check (needs the R package ranger)
   Rscript scripts/sim_power_v2.R                       # N_REP=100 CORES=4 for the full run (slow)
   python3 scripts/sensitivity_reliability.py           # reliability sensitivity table from the simulation
   python3 scripts/make_paper_stats.py                  # tables/figure statistics
   Rscript scripts/make_bym2_spatial_map.R              # spatial-effect map
   ```
   `NO_LST=1` drops the MODIS lake-temperature covariate, which is how the reported results were produced.
3. **Full rebuild from raw data.** Obtain the raw data (`data/README.md`), then in order:
   `get_cyan_data_full.py`, `get_nhdplus_national.R`, `combine_nhdplus_states.R`, the three `get_census_*` / `get_acs_age_groups.py` scripts (with `CENSUS_API_KEY` set in your shell), `get_usgs_surface_supply.py`, `get_well_water_usgs.py`, `build_alt_cei_v2.py` (slow), `assemble_panel_v2.py`, `build_pooled_panel.py`, then the model scripts above.
   **Known gap:** `assemble_panel_v2.py` starts from an earlier panel (`outputs/county_year_panel_full.csv`, built by `build_county_year_panel.py` from CDC annual mortality, population, well water and LST) that is not included, and the county adjacency graph was written by the earlier `fit_bym2_model_mi.R`. The graph is included; rebuilding the earlier panel needs the CDC WONDER annual exports and the legacy script.
4. Validation against measured toxin: `python3 scripts/validate_nla2022.py` (needs EPA NLA 2022 files).

## Scripts
| Role | Scripts |
|---|---|
| Data acquisition | `get_cyan_data_full.py`, `get_nhdplus_national.R`, `combine_nhdplus_states.R`, `get_census_demographics.py`, `get_census_poverty.py`, `get_acs_age_groups.py`, `get_usgs_surface_supply.py`, `get_well_water_usgs.py`, `get_modis_lst.py` |
| Exposure and panels | `build_alt_cei_v2.py` (final lake-to-county overlay), `assemble_panel_v2.py`, `build_pooled_panel.py` |
| Models | `fit_bym2_final_v3.R` (primary and robustness), `fit_pooled_age_std.R`, `fit_dml_pooled.R` |
| Simulation and sensitivity | `sim_power_v2.R`, `sensitivity_reliability.py` |
| Validation and figures | `validate_nla2022.py`, `make_bym2_spatial_map.R`, `make_paper_stats.py`, `outputs/figures/make_figures.py` |
| Exploratory (not in the main results) | `sdwis_feasibility.py` (drinking-water-weighted exposure pilot, Ohio and Iowa 2021-22) |
| Earlier versions, kept for history | `fit_bym2_model.R`, `fit_bym2_model_mi.R`, `fit_bym2_final.R`, `fit_bym2_final_v2.R`, `build_alt_cei.py`, `build_county_year_panel.py` (still imported as a helper), `get_cyan_data.py`, `check_collinearity.R`, `build_suppression_map.py` and `shp_dbf.py` (written for a sandbox without geospatial libraries) |

## Corrections made during development
Four errors were found and fixed before the reported results, and no result from earlier versions is used:
the CyAN digital-number codes (254 land, 255 no data) were averaged into lake values; lakes were assigned to
the county of their centroid (dropping Great Lakes counties); the "well water" covariate had been built from a
Census heating-fuel table; and the chlorophyll threshold constant was an uncited placeholder.

## Limitations
Ecological design; satellite bloom index is a weak proxy for toxin exposure at county level; exposure is not
population- or intake-weighted; suppressed counts and post hoc analysis choices (the primary exposure definition
was fixed after an earlier flawed pipeline was seen) are described in the paper. Interpret rate ratios as
sensitivity to these limits, not as causal estimates.

## AI assistance
Claude (Anthropic) was used for coding help, debugging, and drafting, as disclosed in the accompanying paper and
ISEF paperwork. The author ran the analyses and made the research decisions.

## License and citation
Code: MIT (see `LICENSE`). Derived tables are built from public data; please also cite the original data
providers listed in `data/README.md`. Citation metadata: `CITATION.cff`.
