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
  - `spatial_inference.R`: polygon centroids, residual spatial correlation,
    Conley cutoff
  - `cross_validation.R`: leave-fold-out RMSE on demeaned data
  - `response_functions.R`: distributed lags, bin coefficients, response
    curves, and their plots
  - `heterogeneity.R`: Carleton-style and lagged-income adaptation responses
  - `climate_runs.R`: warming/cooling runs and climate velocity
  - `projection.R`: growth effects along climate trajectories
  - `correlation_tests.R`: pairwise and trend correlation tests
  - `specification_grid.R`: one specification over a grid of data choices:
    samples, coefficient and dispersion tables, and the comparison figures
  - `analysis_spec.R`: `main_spec()`, `build_main_dat()`, `fit_main()`: the
    main_analysis.Rmd specification shared by the analysis scripts
  - `climate_slope_bins.R`: free climate-bin slopes and linearity tests
  - `functional_forms.R`: hinged anomalies, cold/hot symmetry test,
    coefficient stability
- `test_functions.Rmd`: tests of alternative climate-response functions
- `compare_signed_bins_across_datasets.R`: signed-bin specification across
  DOSE, KUMMU, PWT and WB x ERA5, CRU, UDel x population and area weighting
- `compare_dose_gdp_definitions.R`: signed-bin specification across the five
  DOSE GDP-per-capita definitions x ERA5, CRU, UDel (population weighting)
- `compare_dose_gdp_definitions_lags.R`: the same grid for the distributed-lag
  model l(dTM, 0:10) + l(dTM, 0:10):mean_TM_all + the same for dRR
- `plot_bhm_climate_slope_bins.R`, `plot_bhm_climate_slope_national.R`: free
  temperature slope per bin of long-run climate vs the linear BHM interaction
  (DOSE; PWT with DOSE alongside)
- `plot_preferred_functional_form.R`: hinged-anomaly powers vs the signed bins,
  cold/hot symmetry, stability of the base coefficients
- `compare_dose_conley_se.R`: country, region and Conley standard errors for
  the preferred and quadratic DOSE specifications
- `compare_dose_climate_levels.R`: DOSE growth on the world, national-minus-
  world and regional-minus-national components of zTM and zRR (CRU TS,
  area-weighted), under year, no-year and country-year fixed effects
- `harmonize_dose_pwt.R`: DOSE-PWT harmonization

## Fast versus slow warming

Implementation of `FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md` on DOSE
(`grp_pc_lcu_2015`) x ERA5 (area weighted). The frozen design is
`config/fast_slow_design.yml`; every output carries its SHA-256. Run in order
from the project root:

    Rscript tests/test_fast_slow_warming.R      # synthetic fixtures (a)-(g)
    Rscript fast_slow_stage1_panel.R            # Stages 0-1: manifest, panel, QA
    Rscript fast_slow_stage2_support.R          # Stage 2: outcome-blind support, power
    Rscript fast_slow_stage3_7_models.R         # Stages 3-7: models, paths, decisions
    Rscript fast_slow_stage8_12_checks.R        # Stages 8-12: anomaly model, robustness, falsification
    Rscript fast_slow_report.R                  # results/fast_slow_warming/report.md

Functions: `functions/fast_slow_rates.R` (the trailing 5-year OLS slope in
C/year, hinges, anomalies, equal-endpoint paths and their contrasts),
`functions/fast_slow_inference.R` (two-way CGM, wild cluster bootstrap-t with
synchronized country/year Rademacher draws, Conley HAC with temporal lags, the
section 7 decision rules), `functions/fast_slow_support.R` (support audit,
fallback contrast, design-range power), `functions/fast_slow_pipeline.R`
(design hash, stage status, registry, fixed-effect designs). Environment
variables `FS_BOOT_REPS`, `FS_POWER_SIMS`, `FS_HUBER_REPS`, `FS_PERM_REPS`
override the replication counts for quick runs.

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
| DOSE v2.14 | GADM1 | `grp_pc_lcu2015_usd` | constant 2015 local prices, at the 2015 exchange rate (2015 US$) |
| KUMMU 2025 | GADM1 | total GDP / population | PPP, constant 2021 int. $ |
| PWT 11.0 | GADM0 | `rgdpna / pop` | constant national prices, 2021 PPP US$ |
| WB WDI | GADM0 | `NY.GDP.PCAP.KD` | constant 2015 local prices, at the 2015 official exchange rate (2015 US$) |

All four growth rates are real growth, so they are comparable. In levels, DOSE
and WB are in exchange-rate dollars and KUMMU and PWT in PPP dollars; levels are
comparable within each pair but not across them (exchange rates understate
poorer countries' incomes relative to PPP).
No dataset is interpolated or extrapolated.
DOSE's own `grp_pc_usd_2015` is not used: its growth includes exchange-rate
movements against the dollar (see `functions/econ_panel.R`).
