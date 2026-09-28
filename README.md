# Econometrics

Climate-response models of regional (DOSE, Kummu) and national (PWT, WB)
GDP-per-capita growth.

## Layout

- `load_functions.R`: attaches dplyr, fixest and ggplot2 and sources every
  script in `functions/`. Every analysis script and notebook starts with
  `source("load_functions.R")`.
- `functions/`: one script per function or family of functions:
  - `build_dat.R`: `panel_config()`, `build_dat()` (the estimation panel) and
    `add_weather_variables()`
  - `econ_panel.R`, `climate_panel.R`: the economic and climate halves of the
    panel, built from raw sources. Economic sources: `DOSE` and `KUMMU`
    (GADM1), `PWT` and `WB` (GADM0)
  - `gdp_levels.R`: total GDP and population of any source in a common GADM
    layout
  - `national_subnational.R`: `compare_national_subnational()` (national vs
    aggregated subnational growth, with exact regional attribution) and
    `rank_subnational_outliers()`
  - `outlier_sensitivity.R`: re-estimation after dropping ranked outliers
  - `kummu_grid.R`: builds and caches the Kummu 2025 GADM1 panel from the
    raw grids (`load_kummu2025_grid()`)
  - `national_subnational_report.R`: `national_subnational_growth_report()`
  - `climate_download.R`: Weighted Climate Dataset downloader (`wcd_get()`)
  - `climate_transforms.R`: Hamilton filter, standardized anomalies
  - `data_cache.R`: data paths and the `data/cache` panel cache
  - `panel_utils.R`: year-matched lags, power columns
  - `binning.R`: signed anomaly bins, quantile bins, tail pooling
  - `model_spec.R`: `make_panel_formula()`, `fit_panel_model()`
  - `inference.R`: coefficient tables, linear combinations, Wald rows, EB
    shrinkage, year-block bootstrap
  - `spatial_inference.R`: residual spatial correlation, Conley cutoff
  - `cross_validation.R`: leave-fold-out RMSE on demeaned data
  - `response_functions.R`: distributed lags, bin coefficients, response
    curves, and their plots
  - `heterogeneity.R`: Carleton-style and lagged-income adaptation responses
  - `climate_runs.R`: warming/cooling runs and climate velocity
  - `projection.R`: growth effects along climate trajectories
  - `correlation_tests.R`: pairwise and trend correlation tests
- `test_functions.Rmd`: tests of alternative climate-response functions
- `harmonize_dose_pwt.R`: DOSE-PWT harmonization

National-vs-subnational growth comparison, ranked subnational outliers and the
outlier sensitivity of the signed-bin model, for every pair of sources:

    source("load_functions.R")
    national_subnational_growth_report(c("DOSE", "KUMMU"), c("PWT", "WB"))

## Data

Raw inputs live in `data/`; derived panels are cached in `data/cache/` and
rebuilt automatically when the specification or builder code changes.
`data/kummu2025/` holds the Kummu et al. (2025) rasters, their Zenodo metadata
and the GADM 4.1 ADM1 polygons; the KUMMU panel is built from them on first
use (several minutes). The folder is git-ignored because the rasters exceed
GitHub's 100 MB file limit.

| Source | Level | Variable | Prices |
|---|---|---|---|
| DOSE v2.14 | GADM1 | `grp_pc_lcu_2015` | constant 2015 local currency |
| KUMMU 2025 | GADM1 | total GDP / population | PPP, constant 2021 int. $ |
| PWT 11.0 | GADM0 | `rgdpna / pop` | constant national prices, 2021 PPP US$ |
| WB WDI | GADM0 | `NY.GDP.PCAP.KD` | constant 2015 US$, market rates |

All four growth rates are real growth, so they are comparable. Levels are not
comparable across sources, nor across countries for DOSE (local currencies).
No dataset is interpolated or extrapolated.
DOSE's own `grp_pc_usd_2015` is not used: its growth includes exchange-rate
movements against the dollar (see `functions/econ_panel.R`).
