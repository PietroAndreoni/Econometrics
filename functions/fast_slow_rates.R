# Warming-speed variables and equal-endpoint temperature paths for the
# fast-versus-slow warming design (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md,
# sections 5 and 6).
#
# One function, fs_add_climate_terms(), builds every climate regressor from a
# unit's annual temperature series. It is applied to the observed climate panel
# (before any outcome is merged) and to the simulated fast and slow paths, so
# estimation and prediction share centring, scaling, knots and the rate
# calculation exactly (assertion 11).

# Trailing OLS slope of annual temperature on time through the k years ending
# in year t, in degrees C per year:
#   r_t = sum_{j=0}^{k-1} (j - jbar) T_{t-k+1+j} / sum_{j=0}^{k-1} (j - jbar)^2.
# Temperatures are matched on calendar year, so the slope is NA unless all k
# years t-k+1, ..., t are present; a gap is never bridged and no year after t is
# used. fit_trailing_climate_velocity() (climate_runs.R) is the same slope in
# degrees per decade; this is the design's definition and unit.
trailing_temperature_slope <- function(temp, year, k = 5L) {
  k <- as.integer(k)
  stopifnot(length(temp) == length(year), k >= 3L, !anyDuplicated(year))
  j <- seq_len(k) - 1L
  weights <- (j - mean(j)) / sum((j - mean(j))^2)
  slope <- numeric(length(temp))
  for (m in seq_along(j)) {
    # Offset j[m] is year t - (k - 1) + j[m].
    slope <- slope + weights[m] * temp[match(year - (k - 1L) + j[m], year)]
  }
  slope
}

# Slope through the k years t+1, ..., t+k: the future-only negative-control
# exposure of Stage 12. Same weights and unit as trailing_temperature_slope().
future_temperature_slope <- function(temp, year, k = 5L) {
  shifted <- trailing_temperature_slope(temp, year, k)
  shifted[match(year + k, year)]
}

# Acute one-year change T_t - T_{t-1}, NA unless year t-1 is present.
acute_temperature_change <- function(temp, year) {
  temp - temp[match(year - 1L, year)]
}

# Leave-current-year-out trailing-normal anomaly T_t - mean(T_{t-1..t-L}), NA
# unless all L previous years are present.
trailing_normal_anomaly <- function(temp, year, normal_years = 20L) {
  normal <- numeric(length(temp))
  for (lag in seq_len(normal_years)) {
    normal <- normal + temp[match(year - lag, year)]
  }
  temp - normal / normal_years
}

# Warming and cooling hinges and their squares, named <prefix>_pos,
# <prefix>_pos2, <prefix>_neg, <prefix>_neg2.
rate_hinges <- function(rate, prefix = "r") {
  pos <- pmax(rate, 0)
  neg <- pmax(-rate, 0)
  out <- list(pos, pos^2, neg, neg^2)
  names(out) <- paste0(prefix, c("_pos", "_pos2", "_neg", "_neg2"))
  tibble::as_tibble(out)
}

# Term names of each block, used to build formulas and to split path contrasts
# into their level and speed components.
fs_rate_terms <- function(prefix = "r") {
  paste0(prefix, c("_pos", "_pos2", "_neg", "_neg2"))
}

# Transformation constants frozen before any outcome model: the temperature
# reference T0 and the precipitation centring and scaling (pooled over the
# analytic region-years, without outcomes), and the fixed spline knots.
fs_terms_spec <- function(temp, precip, spline_knot_quantiles = NULL,
                          spline_boundary_quantiles = c(0.01, 0.99)) {
  spec <- list(
    T0 = mean(temp),
    P_mean = mean(precip),
    P_sd = stats::sd(precip)
  )
  if (!is.null(spline_knot_quantiles)) {
    spec$spline_knots <- unname(stats::quantile(temp, spline_knot_quantiles))
    spec$spline_boundary <- unname(
      stats::quantile(temp, spline_boundary_quantiles)
    )
  }
  spec
}

# Natural cubic spline basis of temperature with the frozen knots; columns
# Tns1, ..., Tnsm. predict() on the same basis object keeps estimation and path
# prediction identical.
fs_spline_basis <- function(temp, spec) {
  basis <- splines::ns(
    temp, knots = spec$spline_knots, Boundary.knots = spec$spline_boundary
  )
  out <- unclass(basis)
  attributes(out) <- list(dim = dim(out))
  colnames(out) <- paste0("Tns", seq_len(ncol(out)))
  tibble::as_tibble(out)
}

# Every climate regressor for a panel of annual series, computed within `id`
# from TM (degrees C) and RR (metres):
#   Tc, Tc2          T - T0 and its square
#   Pz, Pz2          (P - P_mean) / P_sd and its square
#   rate<k>          trailing k-year slope for each k in `rate_windows`
#   r<k>_pos, ...    hinges of each slope (r_* for the primary window)
#   d1, d1_*         acute change and its hinges
#   anom<L>, a<L>_*  trailing-normal anomalies and their hinges
#   rfut<k>, rf_*    future-only slope (negative control) and its hinges
#   Tns*             spline basis, when the spec carries knots
# Lags and leads are year-matched within id. The data must not contain an
# outcome when called on the observed panel (assertion 8).
fs_add_climate_terms <- function(data, spec, id = "GID_1",
                                 primary_window = 5L,
                                 rate_windows = c(3L, 5L, 10L),
                                 normal_years = c(15L, 20L, 30L),
                                 temp_lags = 0:10, lead_years = 1:2) {
  outcome_like <- grep("grp|gdp|dlg", names(data), value = TRUE)
  if (length(outcome_like)) {
    stop("fs_add_climate_terms() must run before outcomes are merged; found ",
         paste(outcome_like, collapse = ", "))
  }
  rate_windows <- sort(unique(c(primary_window, rate_windows)))

  add_unit <- function(d) {
    d <- d[order(d$year), , drop = FALSE]
    temp <- d$TM
    year <- d$year
    for (k in rate_windows) {
      d[[paste0("rate", k)]] <- trailing_temperature_slope(temp, year, k)
      d <- dplyr::bind_cols(d, rate_hinges(d[[paste0("rate", k)]],
                                           paste0("r", k)))
    }
    d$rfut <- future_temperature_slope(temp, year, primary_window)
    d <- dplyr::bind_cols(d, rate_hinges(d$rfut, "rf"))
    d$d1 <- acute_temperature_change(temp, year)
    d <- dplyr::bind_cols(d, rate_hinges(d$d1, "d1"))
    for (L in normal_years) {
      d[[paste0("anom", L)]] <- trailing_normal_anomaly(temp, year, L)
      d <- dplyr::bind_cols(d, rate_hinges(d[[paste0("anom", L)]],
                                           paste0("a", L)))
    }
    for (lag in setdiff(temp_lags, 0L)) {
      d[[paste0("TM_l", lag)]] <- temp[match(year - lag, year)]
    }
    for (lead in lead_years) {
      d[[paste0("TM_f", lead)]] <- temp[match(year + lead, year)]
    }
    d
  }

  out <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id))) %>%
    dplyr::group_modify(~ add_unit(.x)) %>%
    dplyr::ungroup()

  prim <- paste0("r", primary_window)
  out$r_pos <- out[[paste0(prim, "_pos")]]
  out$r_pos2 <- out[[paste0(prim, "_pos2")]]
  out$r_neg <- out[[paste0(prim, "_neg")]]
  out$r_neg2 <- out[[paste0(prim, "_neg2")]]
  out$rate <- out[[paste0("rate", primary_window)]]
  fs_apply_spec(out, spec)
}

# The columns that depend on the frozen constants: Tc, Tc2 (and their lags and
# leads), Pz, Pz2 and the spline basis. Rerunning it with the final spec
# replaces them without recomputing the rates.
fs_apply_spec <- function(out, spec) {
  out <- out[, !grepl("^(Tc|Pz|Tns)", names(out)), drop = FALSE]
  out$Tc <- out$TM - spec$T0
  out$Tc2 <- out$Tc^2
  for (col in grep("^TM_[lf][0-9]+$", names(out), value = TRUE)) {
    out[[sub("^TM_", "Tc_", col)]] <- out[[col]] - spec$T0
  }
  if (!is.null(out$RR)) {
    out$Pz <- (out$RR - spec$P_mean) / spec$P_sd
    out$Pz2 <- out$Pz^2
  }
  if (!is.null(spec$spline_knots)) {
    out <- dplyr::bind_cols(out, fs_spline_basis(out$TM, spec))
  }
  out
}

# Equal-endpoint temperature path for one baseline B, relative years
# s = -history, ..., horizon: flat at B through s = 0 (the flat history that
# initializes every lag and window), then T_s = B + M min(s / D, 1). `delay`
# shifts the start of the ramp (the equal-degree-year diagnostic).
fs_path_temperature <- function(baseline, total_warming, duration, horizon,
                                history = 30L, delay = 0) {
  s <- seq.int(-history, horizon)
  ramp <- pmin(pmax(s - delay, 0) / duration, 1)
  tibble::tibble(s = s, TM = baseline + total_warming * ramp)
}

# Climate regressors along the fast and slow paths of every region, from the
# same fs_add_climate_terms() used on the data. `regions` holds GID_1, the
# baseline temperature B (1991-2010 mean) and precipitation P (held at the
# same baseline mean in both paths). Returns one row per region, path and
# relative year s >= 1.
fs_path_terms <- function(regions, spec, total_warming, fast_years, slow_years,
                          horizon, history = 30L, fast_delay = 0,
                          primary_window = 5L) {
  paths <- dplyr::bind_rows(
    fast = tibble::tibble(duration = fast_years, delay = fast_delay),
    slow = tibble::tibble(duration = slow_years, delay = 0),
    .id = "path"
  )
  # Every rate, change and anomaly is a contrast of temperatures whose weights
  # sum to zero, so it is the same for every baseline B. The terms are built
  # once per path by fs_add_climate_terms() on the B = 0 template; only the
  # level columns (TM and its lags and leads) are then shifted by each B, and
  # the frozen constants reapplied (tests/test_fast_slow_warming.R checks this
  # against building every region separately).
  templates <- dplyr::bind_rows(lapply(seq_len(nrow(paths)), function(p) {
    fs_path_temperature(0, total_warming, paths$duration[p], horizon,
                        history, paths$delay[p]) %>%
      dplyr::mutate(unit = paths$path[p], year = s + 3000L, RR = 0)
  }))
  templates <- fs_add_climate_terms(templates, spec, id = "unit",
                                    primary_window = primary_window) %>%
    dplyr::rename(path = unit) %>%
    dplyr::filter(s >= 1L)
  level_cols <- grep("^TM(_[lf][0-9]+)?$", names(templates), value = TRUE)
  out <- tidyr::crossing(regions %>% dplyr::select(GID_1, B, P), templates)
  for (col in level_cols) out[[col]] <- out[[col]] + out$B
  out$RR <- out$P
  fs_apply_spec(out, spec)
}

# Contrast vectors (fast minus slow, cumulated over s = 1..h, averaged over
# regions with `weights`) for each horizon: one row per horizon, one column per
# term in `terms`. Because every model is linear in its parameters, the path
# gap in log points is contrast %*% beta. `blocks` names the term groups that
# are switched on (e.g. level only).
fs_path_contrast <- function(path_terms, terms, horizons, weights = NULL) {
  if (is.null(weights)) {
    weights <- tibble::tibble(GID_1 = unique(path_terms$GID_1), weight = 1)
  }
  wide <- path_terms %>%
    dplyr::inner_join(weights, by = "GID_1") %>%
    dplyr::mutate(weight = weight / sum(weight[path == "fast" & s == 1L]))
  missing_terms <- setdiff(terms, names(wide))
  if (length(missing_terms)) {
    stop("Path terms missing: ", paste(missing_terms, collapse = ", "))
  }
  sign <- ifelse(wide$path == "fast", 1, -1)
  X <- as.matrix(wide[, terms]) * (sign * wide$weight)
  if (anyNA(X)) stop("Path regressors contain NA; extend `history`.")
  by_s <- rowsum(X, wide$s)
  cumulative <- apply(by_s, 2, cumsum)
  if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, nrow = 1)
  s_values <- as.integer(rownames(by_s))
  out <- cumulative[match(horizons, s_values), , drop = FALSE]
  dimnames(out) <- list(paste0("h", horizons), terms)
  out
}
