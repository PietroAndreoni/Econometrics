# Climate panels derived from the Weighted Climate Dataset (Gortan, Testa,
# Fagiolo and Lamperti 2024) through wcd_get() in climate_download.R.
# Every setting comes from a panel_config() list (see build_dat.R).

# Annual series for one geographical resolution. Temperature is the mean of the
# monthly values and precipitation their sum, which is how the dashboard
# aggregates them. Precipitation is converted from millimetres to metres so the
# level bins stay in metres. Downloads are cached on disk by the helper, so
# repeated calls cost nothing after the first.
load_climate_series <- function(geo_resolution, config = panel_config()) {
  is_national <- identical(geo_resolution, "gadm0")
  id_col <- if (is_national) "GID_0" else "GID_1"
  # The dashboard names gadm1 units "AFG_1_1" while GADM and the economic
  # panels use "AFG.1_1", so the first underscore becomes a dot. gadm0 ids are
  # plain ISO3 codes and need no change. Units that do not match the expected
  # shape are dropped, as in prepare_climate_data.R.
  unit_pattern <- if (is_national) "^[A-Z0-9]{3}$" else "^[A-Z0-9]{3}\\..+"

  fetch <- function(variable, value_code, value_name) {
    out <- tibble::as_tibble(wcd_get(
      variable = variable,
      source = config$climate_source,
      geo_resolution = geo_resolution,
      weight = config$climate_weight,
      weight_year = config$climate_weight_year,
      time_frequency = "yearly",
      verbose = FALSE
    ))
    missing_cols <- setdiff(c("date", id_col, value_code), names(out))
    if (length(missing_cols)) {
      stop(
        "The downloader returned no ", paste(missing_cols, collapse = ", "),
        " column for ", variable, " at ", geo_resolution, "."
      )
    }
    out %>%
      transmute(
        year = as.integer(date),
        unit = .data[[id_col]],
        !!value_name := .data[[value_code]]
      )
  }

  inner_join(
    fetch("avg. temperature", "tmp", "TM"),
    fetch("precipitation", "pre", "RR"),
    by = c("year", "unit")
  ) %>%
    mutate(
      GID_1 = if (is_national) unit else sub("^([^_]+)_", "\\1.", unit)
    ) %>%
    filter(grepl(unit_pattern, GID_1)) %>%
    transmute(
      gadm_level = geo_resolution,
      climate_source = config$climate_source_code,
      GID_0 = substr(GID_1, 1, 3),
      GID_1,
      year,
      TM,
      RR = RR / 1000
    ) %>%
    filter(!is.na(TM) | !is.na(RR))
}

# Rolling moments are right-aligned and include the current year, exactly as in
# prepare_climate_data.R. The models lag them by one year so that a realisation
# never enters its own baseline.
add_rolling_moments <- function(unit_series, windows) {
  for (variable in c("TM", "RR")) {
    values <- unit_series[[variable]]
    for (window in windows) {
      unit_series[[paste0("mean_", variable, "_", window)]] <-
        zoo::rollmeanr(values, k = window, fill = NA_real_)
      unit_series[[paste0("sd_", variable, "_", window)]] <-
        zoo::rollapplyr(values, width = window, FUN = stats::sd, fill = NA_real_)
    }
  }
  unit_series
}

# Fixed-period moments and trends, plus the trailing rolling moments. The
# rolling operators assume a gap-free annual series within each unit, which is
# checked rather than assumed.
derive_climate_moments <- function(series, config = panel_config()) {
  gaps <- series %>%
    group_by(gadm_level, GID_1) %>%
    summarise(
      complete = n() == max(year) - min(year) + 1L,
      .groups = "drop"
    ) %>%
    filter(!complete)
  if (nrow(gaps)) {
    stop(
      "Gaps in the annual climate series for ", nrow(gaps),
      " unit(s), e.g. ", paste(utils::head(gaps$GID_1, 3), collapse = ", "),
      ". The rolling moments assume consecutive years."
    )
  }

  period_all <- config$period_all
  period_pre <- config$period_pre

  series %>%
    group_by(gadm_level, climate_source, GID_0, GID_1) %>%
    arrange(year, .by_group = TRUE) %>%
    group_modify(~add_rolling_moments(.x, config$climate_windows)) %>%
    mutate(
      dTM = TM - lag_by_year(TM, year),
      dRR = RR - lag_by_year(RR, year),
      mean_TM_all = mean(TM[in_period(year, period_all)], na.rm = TRUE),
      sd_TM_all = sd(TM[in_period(year, period_all)], na.rm = TRUE),
      mean_RR_all = mean(RR[in_period(year, period_all)], na.rm = TRUE),
      sd_RR_all = sd(RR[in_period(year, period_all)], na.rm = TRUE),
      mean_TM_pre = mean(TM[in_period(year, period_pre)], na.rm = TRUE),
      sd_TM_pre = sd(TM[in_period(year, period_pre)], na.rm = TRUE),
      mean_RR_pre = mean(RR[in_period(year, period_pre)], na.rm = TRUE),
      sd_RR_pre = sd(RR[in_period(year, period_pre)], na.rm = TRUE),
      trend_TM = hamilton_trend(
        TM, h = config$hamilton_h, p = config$hamilton_lags
      ),
      trend_RR = hamilton_trend(
        RR, h = config$hamilton_h, p = config$hamilton_lags
      )
    ) %>%
    ungroup()
}

load_climate_panel <- function(geo_resolution, config = panel_config()) {
  cached_panel(
    paste0(
      "climate_", geo_resolution, "_", config$climate_source_code, "_",
      gsub("[^a-z0-9]+", "_", tolower(config$climate_weight))
    ),
    list(
      geo_resolution = geo_resolution,
      climate_source = config$climate_source,
      weight = config$climate_weight,
      weight_year = if (is.null(config$climate_weight_year)) {
        "none"
      } else {
        config$climate_weight_year
      },
      windows = config$climate_windows,
      period_all = config$period_all,
      period_pre = config$period_pre,
      hamilton_h = config$hamilton_h,
      hamilton_p = config$hamilton_lags,
      code = code_fingerprint(
        load_climate_series, add_rolling_moments, derive_climate_moments
      )
    ),
    function() {
      derive_climate_moments(load_climate_series(geo_resolution, config), config)
    }
  )
}
