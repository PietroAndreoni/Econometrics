# Small panel helpers shared by the data builders and the analysis scripts.

# A row is estimable when it has an outcome; everything else is a lag carrier.
est_row <- function(d, outcome = "dlgrp_pc_usd") !is.na(d[[outcome]])

# Match on calendar year so that gaps produce NA rather than an incorrect lag.
lag_by_year <- function(x, year, k = 1) {
  x[match(year - k, year)]
}

in_period <- function(year, period) {
  year >= period[[1]] & year <= period[[2]]
}

# Year-matched lags of several columns at once, within each unit. Lag k of `x`
# is named `l<x>` for k = 1 and `l<k><x>` otherwise, matching the names the
# older prepare_climate_data.R produced. Unlike that script's positional lag, a
# gap in a unit's series yields NA instead of a lag silently spanning the gap.
add_year_lags <- function(data, vars, lags, id = "GID_1", year = "year",
                          prefix = "l") {
  data <- dplyr::group_by(data, dplyr::across(dplyr::all_of(id)))
  for (k in lags) {
    for (var in vars) {
      name <- if (k == 1L) paste0(prefix, var) else paste0(prefix, k, var)
      data <- dplyr::mutate(
        data,
        !!name := lag_by_year(.data[[var]], .data[[year]], k)
      )
    }
  }
  dplyr::ungroup(data)
}

# `x_<power>` columns, e.g. TM_2 for power = 2.
add_power_columns <- function(data, vars, power) {
  dplyr::mutate(
    data,
    dplyr::across(
      dplyr::all_of(vars),
      ~ .x^power,
      .names = paste0("{.col}_", power)
    )
  )
}
