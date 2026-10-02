# The main specification, in one place, for the analysis scripts.
#
# The defaults mirror the shared-specification and build-data chunks of
# main_analysis.Rmd; change them here (or pass arguments) and every script that
# uses main_spec() follows.

main_spec <- function(
    outcome = "dlgrp_pc_usd",
    fixed_effects = "year + GID_1[year] + GID_1[year^2]",
    panel_id = c("GID_1", "year"),
    cluster = ~GID_1,
    base_climate = "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
    econ_year_min = 1950L,
    clim_history_years = 6L,
    dose_variable = "lcu2015_usd"
) {
  list(
    outcome = outcome,
    fixed_effects = fixed_effects,
    panel_id = panel_id,
    cluster = cluster,
    base_climate = base_climate,
    econ_year_min = as.integer(econ_year_min),
    clim_history_years = as.integer(clim_history_years),
    dose_variable = dose_variable
  )
}

# build_dat() with the specification's sample settings; the DOSE GDP definition
# applies only to DOSE. Further arguments (e.g. `baseline`, `deviation`) go to
# build_dat(); `config_args` go to panel_config() (e.g. the climate source).
build_main_dat <- function(econ_data, spec = main_spec(), config_args = list(),
                           ...) {
  config <- do.call(panel_config, c(
    list(econ_year_min = spec$econ_year_min,
         clim_history_years = spec$clim_history_years),
    config_args
  ))
  is_dose <- resolve_econ_source(econ_data) == "DOSE_V2_14"
  build_dat(
    econ_data = econ_data,
    config = config,
    econ_variable = if (is_dose) spec$dose_variable else NULL,
    ...
  )
}

# fit_panel_model() with the specification's outcome, fixed effects, panel and
# clustering. `terms` are the regressors only.
fit_main <- function(terms, model_data, spec = main_spec(), ...) {
  fit_panel_model(
    terms,
    model_data,
    outcome = spec$outcome,
    fixed_effects = spec$fixed_effects,
    panel_id = spec$panel_id,
    cluster = spec$cluster,
    ...
  )
}
