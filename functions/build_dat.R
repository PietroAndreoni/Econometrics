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

# `baseline`, `deviation`, `bin_width` and `signed_quantile` choose how the
# anomalies and bins are defined; see add_weather_variables() below. The
# defaults (standardized anomalies against the lagged 30-year moments)
# reproduce the original panel.
build_dat <- function(
    econ_data = c("DOSE", "DOSE_V2_14", "KUMMU", "KUMMU2025_GRID",
                  "PWT", "PWT110", "WB"),
    climate_series = NULL,
    config = panel_config(),
    econ_variable = NULL,
    baseline = "lag_30",
    deviation = c("standardized", "absolute"),
    bin_width = NULL,
    signed_quantile = c("symmetric", "per_side")
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

  # The climate panel is the left-hand side of the join below. Record every
  # economic observation that this choice excludes rather than silently losing
  # it. The exact keys are attached to the returned panel; the message keeps
  # routine calls readable while making the exclusion visible.
  join_keys <- c("gadm_level", "GID_0", "GID_1", "year")
  unit_keys <- c("gadm_level", "GID_0", "GID_1")
  climate_keys <- climate %>% distinct(across(all_of(join_keys)))
  climate_units <- climate_keys %>%
    distinct(across(all_of(unit_keys))) %>%
    mutate(.climate_unit_available = TRUE)

  excluded_region_years <- econ %>%
    anti_join(climate_keys, by = join_keys) %>%
    left_join(climate_units, by = unit_keys) %>%
    transmute(
      econ_data,
      econ_source,
      gadm_level,
      GID_0,
      GID_1,
      year,
      exclusion_reason = if_else(
        dplyr::coalesce(.climate_unit_available, FALSE),
        "climate_year_unavailable",
        "climate_region_unavailable"
      )
    ) %>%
    arrange(GID_0, GID_1, year)

  if (nrow(excluded_region_years)) {
    exclusion_counts <- excluded_region_years %>%
      count(exclusion_reason, name = "n") %>%
      transmute(label = paste0(exclusion_reason, " = ", n)) %>%
      pull(label)
    message(
      "build_dat(", econ_data, "): excluded ",
      format(nrow(excluded_region_years), big.mark = ","), " of ",
      format(nrow(econ), big.mark = ","), " economic region-years across ",
      n_distinct(excluded_region_years$GID_1), " regions because the selected ",
      "climate panel has no matching key (",
      paste(exclusion_counts, collapse = "; "), "). Exact GID_0, GID_1 and ",
      "year values are available in attr(result, \"excluded_region_years\")."
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
  panel <- climate %>%
    filter(year >= config$econ_year_min - config$clim_history_years) %>%
    left_join(
      econ,
      by = join_keys
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
    add_weather_variables(
      baseline = baseline,
      deviation = deviation,
      bin_width = bin_width,
      signed_quantile = signed_quantile
    )

  attr(panel, "excluded_region_years") <- excluded_region_years
  panel
}

# Weather-response regressors: anomalies, their one-sided magnitudes, and
# unsigned, signed and level bins, each as equal-width and equal-frequency (_q)
# bins. Also centres the cross-sectional adaptation covariates over the
# estimation sample, so the main bin coefficients describe a region at the
# sample means.
#
# `baseline` is the reference climate the anomaly is measured against:
#   "lag_30", "lag_20", "lag_10", "lag_5"  trailing window ending the year
#                                          before (mean_TM_lag_30, ...)
#   "all"                                  fixed period_all (mean_TM_all)
#   "pre"                                  fixed period_pre (mean_TM_pre)
# `deviation`:
#   "standardized"  zTM = (TM - mean) / sd, in standard deviations of the
#                   baseline (the default and previous behaviour);
#   "absolute"      zTM = TM - mean, in the variable's units (degrees Celsius,
#                   metres of precipitation), with no SD scaling.
# The names zTM, zRR, abs_zTM, ..., *_bin* are the same in both cases, so a
# model formula works unchanged.
#
# Equal-width bins (TM_bin, TM_bin_fine, TM_bin_signed) are multiples of
# `bin_width` w: coarse |z| breaks 0, 2w, 4w; fine |z| breaks 0, w, ..., 5w;
# signed breaks +-w, ..., +-5w with [-w, w) as the reference bin coded 0. With
# standardized anomalies w = 0.5 SD by default, which reproduces the original
# bins and labels exactly. With absolute anomalies there is no natural unit, so
# the default is w = 0.5 x the pooled SD of the anomaly over the estimation
# sample (per variable); pass `bin_width` (a number, or c(TM = ., RR = .)) to
# set it in the variable's units instead.
#
# Equal-frequency bins have the same number of bins as their equal-width
# counterparts, with quantiles taken over the estimation sample:
#   TM_bin_q, TM_bin_fine_q  3 and 6 quantile bins of |z|;
#   TM_bin_signed_q          5 bins per side plus a reference bin containing
#                            zero, coded -5..5; `signed_quantile` =
#                            "symmetric" (default: breaks mirrored around
#                            zero) or "per_side" (equal frequency within each
#                            side); see
#                            signed_quantile_bin();
#   TM_bin_level_q           7 (TM) and 6 (RR) quantile bins of the level.
# All breaks, widths and options are returned in attr(panel, "weather_bins").
add_weather_variables <- function(panel, baseline = "lag_30",
                                  deviation = c("standardized", "absolute"),
                                  bin_width = NULL,
                                  signed_quantile = c("symmetric", "per_side"),
                                  outcome = "dlgrp_pc_usd") {
  deviation <- match.arg(deviation)
  signed_quantile <- match.arg(signed_quantile)
  if (length(baseline) != 1L || !grepl("^(lag_[0-9]+|all|pre)$", baseline)) {
    stop("`baseline` must be \"lag_<window>\", \"all\" or \"pre\".")
  }
  in_sample <- !is.na(panel[[outcome]])

  panel <- panel %>%
    mutate(
      log_gdp_av_c = log_gdp_av - mean(log_gdp_av[in_sample], na.rm = TRUE),
      mean_TM_sample_c = mean_TM_sample -
        mean(mean_TM_sample[in_sample], na.rm = TRUE),
      mean_RR_sample_c = mean_RR_sample -
        mean(mean_RR_sample[in_sample], na.rm = TRUE)
    )

  # Original labels, used when the bins are the original 0.5-SD bins.
  default_labels <- list(
    coarse = c("0-1", "1-2", ">2"),
    fine = c("0-.5", ".5-1", "1-1.5", "1.5-2", "2-2.5", ">2.5")
  )
  level_bins <- list(
    TM = list(
      breaks = c(-Inf, 0, 5, 10, 15, 20, 25, Inf),
      labels = c("<0", "0-5", "5-10", "10-15", "15-20", "20-25", ">=25")
    ),
    RR = list(
      breaks = c(-Inf, 0.5, 1, 1.5, 2, 2.5, Inf),
      labels = c("<0.5", "0.5-1", "1-1.5", "1.5-2", "2-2.5", ">=2.5")
    )
  )
  bin_info <- list(baseline = baseline, deviation = deviation,
                   signed_quantile = signed_quantile)

  for (v in c("TM", "RR")) {
    mean_col <- paste0("mean_", v, "_", baseline)
    sd_col <- paste0("sd_", v, "_", baseline)
    if (!all(c(mean_col, sd_col) %in% names(panel))) {
      stop("Baseline moments ", mean_col, " and ", sd_col, " are missing; ",
           "check `baseline` against config$climate_windows.")
    }
    anomaly <- panel[[v]] - panel[[mean_col]]
    if (deviation == "standardized") anomaly <- anomaly / panel[[sd_col]]
    magnitude <- abs(anomaly)

    width <- if (!is.null(bin_width)) {
      if (is.null(names(bin_width))) bin_width[[1]] else bin_width[[v]]
    } else if (deviation == "standardized") {
      0.5
    } else {
      0.5 * stats::sd(anomaly[in_sample], na.rm = TRUE)
    }
    original_bins <- deviation == "standardized" && isTRUE(all.equal(width, 0.5))
    coarse_breaks <- c(0, 2 * width, 4 * width, Inf)
    fine_breaks <- c(0, width * seq_len(5), Inf)
    signed_breaks <- width * c(-5:-1, 1:5)
    label_or <- function(default, breaks) {
      if (original_bins) default else interval_labels(breaks)
    }
    signed_q <- signed_quantile_bin(anomaly, n_side = 5L, sample = in_sample,
                                    method = signed_quantile)

    new_cols <- list(
      mean = panel[[mean_col]],
      z = anomaly,
      abs_z = magnitude,
      abs_zp = pmax(anomaly, 0),       # Magnitude of positive deviations.
      abs_zm = pmax(-anomaly, 0),      # Magnitude of negative deviations.
      bin = cut(magnitude, coarse_breaks,
                labels = label_or(default_labels$coarse, coarse_breaks),
                include.lowest = TRUE, right = FALSE),
      bin_fine = cut(magnitude, fine_breaks,
                     labels = label_or(default_labels$fine, fine_breaks),
                     include.lowest = TRUE, right = FALSE),
      bin_level = cut(panel[[v]], level_bins[[v]]$breaks,
                      labels = level_bins[[v]]$labels, right = FALSE),
      bin_signed = signed_anomaly_bin(anomaly, signed_breaks),
      bin_q = quantile_bin(magnitude, 3L, in_sample, lower = 0),
      bin_fine_q = quantile_bin(magnitude, 6L, in_sample, lower = 0),
      bin_level_q = quantile_bin(panel[[v]], length(level_bins[[v]]$labels),
                                 in_sample),
      bin_signed_q = as.vector(signed_q)
    )
    names(new_cols) <- c(
      paste0("mean_", v), paste0("z", v), paste0("abs_z", v),
      paste0("abs_z", v, "p"), paste0("abs_z", v, "m"),
      paste0(v, c("_bin", "_bin_fine", "_bin_level", "_bin_signed", "_bin_q",
                  "_bin_fine_q", "_bin_level_q", "_bin_signed_q"))
    )
    panel[names(new_cols)] <- new_cols

    bin_info[[v]] <- list(
      bin_width = width,
      coarse_breaks = coarse_breaks,
      fine_breaks = fine_breaks,
      signed_breaks = signed_breaks,
      signed_q_breaks = attr(signed_q, "breaks"),
      level_q_levels = levels(new_cols$bin_level_q),
      fine_q_levels = levels(new_cols$bin_fine_q)
    )
  }

  attr(panel, "weather_bins") <- bin_info
  panel
}
