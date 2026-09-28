# Warming and cooling runs, climate velocity, and run-aware lags.
#
# A warming run is a maximal sequence of at least two consecutive calendar years
# with strictly increasing annual temperature. Position 1 is the year
# immediately before the first positive change. Cooling runs are defined
# symmetrically from strictly decreasing temperature, so a local peak or trough
# can end one run and start the opposite one. See add_climate_runs().

# Slope of temperature on year within one run, in degrees per decade. Cooling
# slopes remain negative.
fit_run_rate <- function(year, temperature,
                         direction = c("warming", "cooling")) {
  direction <- match.arg(direction)
  monotone <- if (direction == "warming") {
    all(diff(temperature) > 0)
  } else {
    all(diff(temperature) < 0)
  }
  stopifnot(
    length(year) == length(temperature),
    length(year) >= 2L,
    all(is.finite(year)),
    all(is.finite(temperature)),
    all(diff(year) == 1L),
    monotone
  )

  unname(stats::lm.fit(
    x = cbind(1, year),
    y = temperature
  )$coefficients[[2]]) * 10
}

# Trailing `window_years` slope ending in each year, in degrees per decade. NA
# until a complete, gap-free window is available; uses no future observations.
fit_trailing_climate_velocity <- function(year, temperature, window_years) {
  stopifnot(
    length(year) == length(temperature),
    length(window_years) == 1L,
    window_years >= 3L
  )

  vapply(seq_along(year), function(i) {
    first <- i - window_years + 1L

    if (first < 1L) {
      return(NA_real_)
    }

    idx <- seq.int(first, i)
    year_window <- year[idx]
    temperature_window <- temperature[idx]

    if (any(!is.finite(year_window)) ||
        any(!is.finite(temperature_window)) ||
        !all(diff(year_window) == 1L)) {
      return(NA_real_)
    }

    unname(stats::lm.fit(
      x = cbind(1, year_window),
      y = temperature_window
    )$coefficients[[2]]) * 10
  }, numeric(1))
}

# Centred-window slope (t - radius ... t + radius), in degrees per decade. Uses
# future observations, so it describes rather than predicts; NA near the ends.
fit_local_warming_rate <- function(year, temperature, window_radius = 2) {
  n <- length(year)
  vapply(seq_len(n), function(i) {
    idx <- seq.int(max(1, i - window_radius), min(n, i + window_radius))
    if (length(idx) < (2 * window_radius + 1)) {
      return(NA_real_)
    }
    stats::coef(stats::lm(temperature[idx] ~ year[idx]))[["year[idx]"]] * 10
  }, numeric(1))
}

# Run identifiers for one direction, within a unit sorted by year. The year
# preceding the first change is attached to the same run.
.run_ids <- function(dTM_run, sign) {
  change <- !is.na(dTM_run) & sign * dTM_run > 0
  run_start <- change & !dplyr::lag(change, default = FALSE)
  change_run_id <- if_else(change, cumsum(run_start), NA_integer_)
  coalesce(
    change_run_id,
    if_else(
      dplyr::lead(change, default = FALSE),
      dplyr::lead(change_run_id),
      NA_integer_
    )
  )
}

# Run descriptors for one direction ("warming" or "cooling"). The fixed-baseline
# deviation holds the baseline mean and SD at their position-1 values
# throughout the run; outside a run it equals the moving-baseline deviation.
.add_run_block <- function(data, direction, id, baseline_mean, baseline_sd) {
  run_id <- paste0(direction, "_run_id")
  run <- paste0(direction, "_run")
  mean_f <- paste0("mean_TM_", direction, "_f")
  sd_f <- paste0("sd_TM_", direction, "_f")

  data %>%
    group_by(across(all_of(c(id, run_id)))) %>%
    mutate(
      !!run := !is.na(.data[[run_id]]) & n() >= 2L,
      !!paste0(direction, "_run_length") :=
        if_else(.data[[run]], row_number(), 1L),
      !!paste0(direction, "_run_total_length") :=
        if_else(.data[[run]], n(), 1L),
      # Sum of within-run dTM through this year; zero at position 1.
      !!paste0("cumulative_", direction) :=
        if_else(.data[[run]], TM - first(TM), 0),
      !!mean_f := if_else(
        .data[[run]], first(.data[[baseline_mean]]), .data[[baseline_mean]]
      ),
      !!sd_f := if_else(
        .data[[run]], first(.data[[baseline_sd]]), .data[[baseline_sd]]
      ),
      !!paste0("zTM_", direction, "_f") :=
        (TM - .data[[mean_f]]) / .data[[sd_f]],
      !!paste0(direction, "_run_rate") := if (!first(.data[[run]])) {
        NA_real_
      } else {
        fit_run_rate(year, TM, direction)
      },
      !!run_id := if_else(.data[[run]], .data[[run_id]], NA_integer_)
    ) %>%
    ungroup()
}

# Adds, for every unit: the warming and cooling run descriptors (id, position
# `*_run_length`, eventual `*_run_total_length`, `cumulative_*` change, run
# slope `*_run_rate`, fixed-baseline anomalies `zTM_*_f`), the trailing
# `climate_velocity` over `velocity_years`, its run-restricted versions, and the
# one-sided fixed-baseline magnitudes zTMp_f and zTMm_f.
add_climate_runs <- function(
    data,
    velocity_years,
    baseline_mean = "mean_TM_lag_30",
    baseline_sd = "sd_TM_lag_30",
    id = "GID_1"
) {
  data %>%
    group_by(across(all_of(id))) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      consecutive_year = year == dplyr::lag(year) + 1L,
      dTM_run = if_else(consecutive_year, TM - dplyr::lag(TM), NA_real_),
      warming_run_id = .run_ids(dTM_run, 1),
      cooling_run_id = .run_ids(dTM_run, -1),
      climate_velocity = fit_trailing_climate_velocity(
        year,
        TM,
        window_years = velocity_years
      )
    ) %>%
    .add_run_block("warming", id, baseline_mean, baseline_sd) %>%
    .add_run_block("cooling", id, baseline_mean, baseline_sd) %>%
    mutate(
      # Velocity is active only during the corresponding runs; observations
      # outside those runs are kept with zero.
      climate_velocity_warming = if_else(
        warming_run & climate_velocity > 0,
        climate_velocity,
        0
      ),
      climate_velocity_cooling = if_else(
        cooling_run & climate_velocity < 0,
        climate_velocity,
        0
      ),
      zTMp_f = pmax(zTM_warming_f, 0),
      zTMm_f = pmax(-zTM_cooling_f, 0),
      # Elapsed run position, with the fifth and later years grouped as 5+.
      warming_run_length_capped = pmin(warming_run_length, 5L),
      cooling_run_length_capped = pmin(cooling_run_length, 5L)
    )
}

# Invariants of add_climate_runs(); stops on the first violation.
check_climate_runs <- function(x) {
  stopifnot(all(x$warming_run_length >= 1L))
  stopifnot(all(x$cooling_run_length >= 1L))
  stopifnot(all(x$cumulative_warming >= 0))
  stopifnot(all(x$cumulative_cooling <= 0))
  stopifnot(isTRUE(all.equal(
    x$zTMp_f[x$warming_run_length == 1L],
    x$abs_zTMp[x$warming_run_length == 1L],
    check.attributes = FALSE
  )))
  stopifnot(isTRUE(all.equal(
    x$zTMm_f[x$cooling_run_length == 1L],
    x$abs_zTMm[x$cooling_run_length == 1L],
    check.attributes = FALSE
  )))
  stopifnot(all(x$cumulative_warming[x$warming_run_length == 1L] == 0))
  stopifnot(all(x$cumulative_cooling[x$cooling_run_length == 1L] == 0))
  # Years with exactly zero temperature change (repeated values in the source
  # series) and no change on either side belong to neither run.
  outside_runs <- sum(!(x$warming_run | x$cooling_run))
  if (outside_runs) {
    warning(
      outside_runs, " observation(s) belong to neither a warming nor a ",
      "cooling run (flat temperature)."
    )
  }
  stopifnot(all(!x$warming_run | x$warming_run_total_length >= 2L))
  stopifnot(all(
    !x$warming_run | x$warming_run_length <= x$warming_run_total_length
  ))
  stopifnot(all(!x$cooling_run | x$cooling_run_total_length >= 2L))
  stopifnot(all(
    !x$cooling_run | x$cooling_run_length <= x$cooling_run_total_length
  ))
  invisible(x)
}

# Lags 1..max_lag of `value_column` that are zero whenever looking back that far
# would leave the current run, named <prefix><k>_same.
add_same_run_lags <- function(
    data,
    value_column,
    run_length_column,
    prefix,
    max_lag,
    id = "GID_1"
) {
  data <- data %>%
    arrange(across(all_of(c(id, "year")))) %>%
    group_by(across(all_of(id)))

  for (lag_number in seq_len(max_lag)) {
    lag_column <- paste0(prefix, lag_number, "_same")
    data <- data %>%
      mutate(
        !!lag_column := if_else(
          .data[[run_length_column]] > lag_number,
          dplyr::lag(.data[[value_column]], n = lag_number),
          0
        )
      )
  }

  ungroup(data)
}

# Association between two binned indicators: Spearman's rho treats the bins as
# ordered, Cramer's V as categorical. scope = "Both non-zero" drops the shared
# "0" reference level before computing both.
summarise_bin_association <- function(
    data,
    cumulative_bin,
    velocity_bin,
    direction,
    scope = c("All complete", "Both non-zero")) {
  scope <- match.arg(scope)

  pairs <- data %>%
    transmute(
      cumulative = .data[[cumulative_bin]],
      velocity = .data[[velocity_bin]]
    ) %>%
    filter(!is.na(cumulative), !is.na(velocity))

  if (scope == "Both non-zero") {
    pairs <- pairs %>%
      filter(as.character(cumulative) != "0", as.character(velocity) != "0")
  }

  contingency <- table(
    droplevels(pairs$cumulative),
    droplevels(pairs$velocity)
  )
  enough_bins <- nrow(contingency) >= 2L && ncol(contingency) >= 2L

  spearman_rho <- if (enough_bins) {
    stats::cor(
      as.integer(pairs$cumulative),
      as.integer(pairs$velocity),
      method = "spearman"
    )
  } else {
    NA_real_
  }

  cramers_v <- if (enough_bins) {
    chi_squared <- suppressWarnings(stats::chisq.test(
      contingency,
      correct = FALSE
    ))$statistic
    sqrt(
      unname(chi_squared) /
        (sum(contingency) * min(nrow(contingency) - 1L, ncol(contingency) - 1L))
    )
  } else {
    NA_real_
  }

  tibble::tibble(
    direction = direction,
    scope = scope,
    observations = nrow(pairs),
    spearman_rho = spearman_rho,
    cramers_v = cramers_v
  )
}
