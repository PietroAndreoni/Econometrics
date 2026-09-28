# Estimation panel: economic outcome joined to the derived climate series.
#
# Both halves are built from raw sources and cached (see econ_panel.R and
# climate_panel.R), so a fresh checkout needs no preprocessing script: the
# economic series come from the DOSE csv or the PWT workbook, the climate series
# from the Weighted Climate Dataset.

# Every setting of the panel lives here rather than in globals, so analysis
# scripts can build variants side by side. CRU TS starts in 1901, so the 30-year
# trailing window is complete from 1930 rather than from 1970 as it is with
# ERA5's 1940 start. "concurrent population" weights each year by that year's
# population instead of by a single base year, so it carries no base year and
# `climate_weight_year` is ignored.
# `clim_history_years` climate years are kept before `econ_year_min` as lag
# carriers: they have no outcome but let lags and trailing windows resolve.
panel_config <- function(
    econ_year_min = 1950L,
    climate_source = "CRU TS",
    climate_weight = "concurrent population",
    climate_weight_year = NULL,
    climate_windows = c(5L, 10L, 20L, 30L),
    period_all = c(1990L, 2019L),
    period_pre = c(1960L, 1989L),
    hamilton_lags = 4L,
    hamilton_h = 2L,
    clim_history_years = 6L
) {
  list(
    econ_year_min = econ_year_min,
    climate_source = climate_source,
    climate_source_code = .wcd_match(climate_source, .WCD_SOURCES, "source"),
    climate_weight = climate_weight,
    climate_weight_year = climate_weight_year,
    climate_windows = as.integer(climate_windows),
    period_all = as.integer(period_all),
    period_pre = as.integer(period_pre),
    hamilton_lags = as.integer(hamilton_lags),
    hamilton_h = as.integer(hamilton_h),
    clim_history_years = as.integer(clim_history_years)
  )
}

# "DOSE" and "PWT" are version-agnostic aliases, so a version bump is a change
# to ECON_RAW_FILES and to the alias target rather than to every caller.
ECON_SOURCE_ALIASES <- c(
  PWT = "PWT110",
  PWT110 = "PWT110",
  WB = "WB",
  DOSE = "DOSE_V2_14",
  DOSE_V2_14 = "DOSE_V2_14",
  KUMMU = "KUMMU2025_GRID",
  KUMMU2025_GRID = "KUMMU2025_GRID"
)
GADM_LEVEL_BY_SOURCE <- c(
  PWT110 = "gadm0",
  WB = "gadm0",
  DOSE_V2_14 = "gadm1",
  KUMMU2025_GRID = "gadm1"
)

# Canonical stored-source name for an alias, with a helpful error.
resolve_econ_source <- function(econ_data) {
  if (length(econ_data) != 1L || !econ_data %in% names(ECON_SOURCE_ALIASES)) {
    stop(
      "Unknown economic source '", paste(econ_data, collapse = ", "),
      "'. Available: ", paste(names(ECON_SOURCE_ALIASES), collapse = ", ")
    )
  }
  unname(ECON_SOURCE_ALIASES[[econ_data]])
}

build_dat <- function(
    econ_data = c("DOSE", "DOSE_V2_14", "KUMMU", "KUMMU2025_GRID",
                  "PWT", "PWT110", "WB"),
    climate_series = NULL,
    config = panel_config(),
    baseline_window = 30L,
    econ_variable = NULL
) {
  econ_data <- match.arg(econ_data)

  stored_econ_source <- resolve_econ_source(econ_data)
  target_gadm_level <- unname(GADM_LEVEL_BY_SOURCE[[stored_econ_source]])

  # `econ_variable` picks a DOSE GDP definition (see DOSE_GDP_VARIABLES);
  # NULL uses the default.
  econ <- load_econ_panel(stored_econ_source, econ_variable) %>%
    select(
      year, GID_0, GID_1, grp_pc_usd, lgrp_pc_usd, dlgrp_pc_usd,
      econ_source, gadm_level
    ) %>%
    distinct() %>%
    filter(
      gadm_level == target_gadm_level,
      year >= config$econ_year_min,
      !is.na(dlgrp_pc_usd)
    )

  if (!nrow(econ)) {
    stop(
      "No ", target_gadm_level, " economic observations found for econ_data = '",
      econ_data, "' (stored source '", stored_econ_source, "')."
    )
  }

  # `climate_series` lets a caller pass an already derived panel, which bypasses
  # the cache entirely when experimenting with alternative moments.
  climate <- climate_series
  if (is.null(climate)) {
    climate <- load_climate_panel(target_gadm_level, config)
  }
  climate <- climate %>% filter(gadm_level == target_gadm_level)

  if (!nrow(climate)) {
    stop(
      "No ", target_gadm_level, " climate observations for econ_data = '",
      econ_data, "'. Check climate_source and climate_weight in the config."
    )
  }

  # Every rolling moment is lagged one year, so a realisation never enters its
  # own baseline: mean_TM_30 -> mean_TM_lag_30, and so on for each window.
  moment_cols <- as.vector(outer(
    c("mean_TM", "mean_RR", "sd_TM", "sd_RR"),
    config$climate_windows,
    paste,
    sep = "_"
  ))
  moment_cols <- intersect(moment_cols, names(climate))

  # The moments are already computed over the full series, so trimming the
  # history here only drops rows that carry no lag any model needs.
  climate %>%
    filter(year >= config$econ_year_min - config$clim_history_years) %>%
    left_join(
      econ,
      by = c("gadm_level", "GID_0", "GID_1", "year")
    ) %>%
    group_by(gadm_level, GID_1) %>%
    filter(any(!is.na(dlgrp_pc_usd))) %>%
    filter(year <= max(year[!is.na(dlgrp_pc_usd)])) %>%
    mutate(
      log_gdp_av = log(mean(grp_pc_usd, na.rm = TRUE)),
      mean_TM_sample = mean(TM[!is.na(dlgrp_pc_usd)], na.rm = TRUE),
      mean_RR_sample = mean(RR[!is.na(dlgrp_pc_usd)], na.rm = TRUE)
    ) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(across(
      all_of(moment_cols),
      ~ lag_by_year(.x, year),
      .names = "lag__{.col}"
    )) %>%
    ungroup() %>%
    rename_with(
      ~ sub("^lag__(mean|sd)_(TM|RR)_", "\\1_\\2_lag_", .x),
      starts_with("lag__")
    ) %>%
    add_weather_variables(baseline_window = baseline_window)
}

# Weather-response regressors built from the lagged trailing moments of window
# `baseline_window`: standardized anomalies, their one-sided magnitudes,
# piecewise-linear segments, and unsigned, signed and level bins. Also centres
# the cross-sectional adaptation covariates over the estimation sample, so the
# main bin coefficients describe a region at the sample means.
add_weather_variables <- function(panel, baseline_window = 30L,
                                  outcome = "dlgrp_pc_usd") {
  mean_TM_base <- panel[[paste0("mean_TM_lag_", baseline_window)]]
  sd_TM_base <- panel[[paste0("sd_TM_lag_", baseline_window)]]
  mean_RR_base <- panel[[paste0("mean_RR_lag_", baseline_window)]]
  sd_RR_base <- panel[[paste0("sd_RR_lag_", baseline_window)]]
  if (is.null(mean_TM_base) || is.null(mean_RR_base)) {
    stop("Lagged moments for a ", baseline_window, "-year window are missing.")
  }
  in_sample <- !is.na(panel[[outcome]])

  panel %>%
    mutate(
      log_gdp_av_c = log_gdp_av - mean(log_gdp_av[in_sample], na.rm = TRUE),
      mean_TM_sample_c = mean_TM_sample -
        mean(mean_TM_sample[in_sample], na.rm = TRUE),
      mean_RR_sample_c = mean_RR_sample -
        mean(mean_RR_sample[in_sample], na.rm = TRUE),
      mean_TM = mean_TM_base,
      mean_RR = mean_RR_base,
      zTM = (TM - mean_TM_base) / sd_TM_base,
      zRR = (RR - mean_RR_base) / sd_RR_base,
      abs_zTM = abs(zTM),
      abs_zRR = abs(zRR),
      # Magnitudes of positive (+) and negative (-) standardized deviations.
      abs_zTMp = pmax(zTM, 0),
      abs_zTMm = pmax(-zTM, 0),
      abs_zRRp = pmax(zRR, 0),
      abs_zRRm = pmax(-zRR, 0),

      # NOTE: the segment caps (1.5, 2.5) exceed the segment widths (1, 1), and
      # RR_15_25 subtracts 0.5 where TM_15_25 subtracts 1.5. Both are kept as in
      # the original notebook so results are unchanged; review before relying
      # on the piecewise-linear test.
      TM_0_05 = pmin(abs_zTM, 0.5),
      TM_05_15 = pmin(pmax(abs_zTM - 0.5, 0), 1.5),
      TM_15_25 = pmin(pmax(abs_zTM - 1.5, 0), 2.5),
      TM_25p = pmax(abs_zTM - 2.5, 0),

      RR_0_05 = pmin(abs_zRR, 0.5),
      RR_05_15 = pmin(pmax(abs_zRR - 0.5, 0), 1.5),
      RR_15_25 = pmin(pmax(abs_zRR - 0.5, 0), 2.5),
      RR_25p = pmax(abs_zRR - 2.5, 0),

      TM_bin = cut(
        abs_zTM,
        breaks = c(0, 1, 2, Inf),
        labels = c("0-1", "1-2", ">2"),
        include.lowest = TRUE,
        right = FALSE
      ),
      RR_bin = cut(
        abs_zRR,
        breaks = c(0, 1, 2, Inf),
        labels = c("0-1", "1-2", ">2"),
        include.lowest = TRUE,
        right = FALSE
      ),
      TM_bin_fine = cut(
        abs_zTM,
        breaks = c(0, 0.5, 1, 1.5, 2, 2.5, Inf),
        labels = c("0-.5", ".5-1", "1-1.5", "1.5-2", "2-2.5", ">2.5"),
        include.lowest = TRUE,
        right = FALSE
      ),
      RR_bin_fine = cut(
        abs_zRR,
        breaks = c(0, 0.5, 1, 1.5, 2, 2.5, Inf),
        labels = c("0-.5", ".5-1", "1-1.5", "1.5-2", "2-2.5", ">2.5"),
        include.lowest = TRUE,
        right = FALSE
      ),
      TM_bin_level = cut(
        TM,
        breaks = c(-Inf, 0, 5, 10, 15, 20, 25, Inf),
        labels = c("<0", "0-5", "5-10", "10-15", "15-20", "20-25", ">=25"),
        right = FALSE
      ),
      RR_bin_level = cut(
        RR,
        breaks = c(-Inf, 0.5, 1, 1.5, 2, 2.5, Inf),
        labels = c("<0.5", "0.5-1", "1-1.5", "1.5-2", "2-2.5", ">=2.5"),
        right = FALSE
      ),
      TM_bin_signed = signed_anomaly_bin(zTM),
      RR_bin_signed = signed_anomaly_bin(zRR)
    )
}
