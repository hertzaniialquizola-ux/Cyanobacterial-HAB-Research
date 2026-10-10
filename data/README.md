# Raw data (not in this repository)

Raw inputs are large (about 18 GB, mostly MODIS and CyAN rasters) and some are subject to the data
provider's terms, so none is committed. Scripts expect them under `data/raw/` as follows. All sources are public.

| Folder under `data/raw/` | Source | How it is obtained |
|---|---|---|
| `cyan/` | EPA/NASA Cyanobacteria Assessment Network (CyAN), merged Sentinel-3 OLCI 7-day CI_cyano composites, 300 m, 2018-2024 | `scripts/get_cyan_data_full.py` (free NASA Earthdata login; `earthaccess`) |
| `nhdplus/` | USGS NHDPlus V2 lakes/reservoirs (>= 1 km2 kept) | `scripts/get_nhdplus_national.R`, then `scripts/combine_nhdplus_states.R` |
| `cdc-wonder/`, `cdc_pooled/` | CDC WONDER Multiple Cause of Death 2018-2024 (kidney N18, liver K74/K76, accidental-injury negative control) | Manual web export (group by County, years 2018-2024, show zero and suppressed values; underlying cause N18; K74 + K76; unintentional injury V01-X59, Y85-Y86 as the negative control; plus national deaths by ten-year age group for age standardization; see `scripts/build_pooled_panel.py`). WONDER's terms forbid publishing counts of 1-9, so raw exports are not shared; the derived tables here contain no such counts. |
| `census/` | ACS 5-year 2022 (age, race/ethnicity, poverty, income); county boundaries `cb_2022_us_county_20m` | `scripts/get_census_demographics.py`, `get_census_poverty.py`, `get_acs_age_groups.py` (need a free Census API key in the environment variable `CENSUS_API_KEY`; never commit a key) |
| `usgs/` | USGS county water use (2015): surface-water supply share; self-supplied (well) population | `scripts/get_usgs_surface_supply.py`, `scripts/get_well_water_usgs.py` (the latter needs the USGS file downloaded manually; use `--diagnose` first) |
| `usda/` | USDA ERS Rural-Urban Continuum Codes 2023 (`rucc_2023.csv`) | Manual download |
| `places/` | CDC PLACES county data (secondary analysis) | Manual download |
| `modis-lst/` | MODIS MOD11A2 v061 land-surface temperature (sensitivity analysis only) | `scripts/get_modis_lst.py` |
| `nla/` | EPA National Lakes Assessment 2022 (microcystin) for validation | Manual download from EPA |
| `epa_sdwis/` | EPA public water system service-area polygons and SDWIS system search (exploratory pilot only) | Manual download from EPA ECHO and the EPA service-area boundaries page |
