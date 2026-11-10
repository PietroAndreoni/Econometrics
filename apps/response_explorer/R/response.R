# Responses and GDP paths of a fitted model of build_models.R for a country
# with given characteristics. Without the cap every quantity is a linear
# combination of the model's coefficients, so its standard error is exact
# (delta method with the clustered covariance) and its 95% CI normal. With the
# cap the response is not linear in the coefficients, so the CI comes from
# draws of the coefficients (percentiles), capped draw by draw.

# The years of the baselines (the level and the NRK filter have no baseline
# mean that moves: as the horizon's).
EXPLORER_WINDOWS <- c(level = Inf, all = Inf, lag_5 = 5, lag_10 = 10, lag_20 = 20, lag_30 = 30, nrk = Inf)

# The first differences of deviations x (series of `years` years stacked by
# path; NULL: one series), with no deviation before the first year of a path.
difference_path <- function(x, years = NULL) {
  at <- lag_rows(length(x), 1, years %||% length(x))
  x - ifelse(is.na(at), 0, x[at])
}
# The deviations a model enters for deviations x: their first differences
# with the style "first difference"; with "both", list(x = x, d = their first
# differences), x for the terms of the unit and d for those of its first
# difference (dev<v>_d_...); otherwise x.
model_deviations <- function(model, x, years = NULL) {
  if (isTRUE(model$both)) return(list(x = x, d = difference_path(x, years)))
  if (isTRUE(model$diff)) difference_path(x, years) else x
}
# The length of deviations x (a vector or such a list), the deviations of the
# terms of the first difference (`diff`) or of the unit, the deviations plus
# the offset `off` of the unit (not of its differences), and the offset alone.
n_inputs <- function(x) length(if (is.list(x)) x$x else x)
term_inputs <- function(x, diff) if (is.list(x)) x[[if (diff) "d" else "x"]] else x
plus_offset <- function(x, off) {
  if (!is.list(x)) return(x + off)
  x$x <- x$x + off
  x
}
offset_inputs <- function(x, off) {
  n <- n_inputs(x)
  if (is.list(x)) list(x = rep(off, n), d = rep(0, n)) else rep(off, n)
}
# The value of the regressor without any anomaly: zero, except for the level
# as it is (not demeaned or differenced: symmetric models and those without
# unit effects), the country's long-run climate over the deviation's scale.
# Effects are those of x + offset less those of the offset.
deviation_offset <- function(model, country, v) {
  if (!identical(model$window, "level") || isTRUE(model$diff) || !identical(of_variable(model$dev_fe, v), "none")) return(0)
  (if (v == "TM") country$mean_T else country$mean_P) / deviation_scale(model, country, v)
}
# The effects of the years so far on log GDP: the sum of the effects on growth
# (by column of a matrix), or with a level outcome the effects themselves.
accumulated <- function(model, x) {
  if (level_outcome(model)) return(x)
  if (is.matrix(x)) matrix(apply(x, 2, cumsum), nrow = nrow(x), dimnames = dimnames(x)) else cumsum(x)
}
N_DRAWS <- 2000

# Value of an arm term (column suffix) at deviation x (SD); "_sym" and "_sym2"
# are the unsplit deviation and its square (symmetric models).
arm_value <- function(suffix, x) {
  switch(suffix,
         "_pos" = pmax(x, 0), "_neg" = pmax(-x, 0),
         "_pos2" = pmax(x, 0)^2, "_neg2" = pmax(-x, 0)^2,
         "_sym" = x, "_sym2" = x^2,
         stop("Unknown arm ", suffix))
}

# The weights of the lags of the NRK filter of variable v for a country:
# phi_j = g_j + h_j vbar, vbar its long-run climate (C, m a year).
filter_lags <- function(model, country, v) {
  f <- model$filter[[v]]
  f$g + f$h * (if (v == "TM") country$mean_T else country$mean_P)
}
# The innovation of the rows `rows` of anomaly matrix A (columns: paths) under
# a filter with lag weights phi: A_t - sum_j phi_j A_t-j, the rows before the
# first anomalies zero.
filtered_rows <- function(A, rows, phi) {
  A <- rbind(matrix(0, length(phi), ncol(A)), A)
  rows <- rows + length(phi)
  out <- A[rows, , drop = FALSE]
  for (j in seq_along(phi)) out <- out - phi[[j]] * A[rows - j, , drop = FALSE]
  out
}
# The deviations (SD, or C and m with sd = 1) of the years from `from` that an
# anomaly path `a` (named by year; nothing before its first year) adds to the
# innovation of the NRK filter with lag weights phi: (a_t - sum_j phi_j
# a_t-j) / sd. The country's own climate has no innovation of its own.
innovation_path <- function(a, phi, sd, from = 0) {
  years <- as.integer(names(a))
  drop(filtered_rows(matrix(a), which(years >= from), phi)) / sd
}

# Restricted cubic spline basis of x with `knots` (Harrell), as rcs_nonlinear()
# in weather_arms_robustness.Rmd: the K - 2 non-linear terms, scaled by the
# squared range of the knots; linear beyond the outer knots.
rcs_nonlinear <- function(x, knots) {
  K <- length(knots)
  cube <- function(z) pmax(z, 0)^3
  out <- vapply(seq_len(K - 2), function(j) {
    (cube(x - knots[[j]]) -
       cube(x - knots[[K - 1]]) * (knots[[K]] - knots[[j]]) / (knots[[K]] - knots[[K - 1]]) +
       cube(x - knots[[K]]) * (knots[[K - 1]] - knots[[j]]) / (knots[[K]] - knots[[K - 1]])) /
      (knots[[K]] - knots[[1]])^2
  }, numeric(length(x)))
  matrix(out, nrow = length(x))
}

# The country's values of the modifiers the deviation terms are interacted
# with, centred as in the estimation: long-run climate (C, m a year) less the
# median of the panel's units, log income, share and log capital per capita
# (PWT) less their means, and the low-income dummy (income below the model's
# threshold); `clim_shift` (C, m) is added to the long-run climate.
# For a spline model, also the non-linear terms <z>_nl<j>, less their value at
# zero (as build_models.R).
country_modifiers <- function(model, country, clim_shift = c(TM = 0, RR = 0)) {
  out <- modifier_columns(model, country, as.list(clim_shift))
  stats::setNames(out[1, ], colnames(out))
}
# The same as a matrix, one row per element of the shifts (list TM, RR) and of
# the country's income, share and capital (scalars or vectors of one length).
# This year's and last year's value of T and P are the country's long-run
# climate plus the anomalies `cur` and `lag` (list TM, RR; C, m), with none by
# default.
modifier_columns <- function(model, country, clim_shift = list(TM = 0, RR = 0),
                             cur = list(TM = 0, RR = 0), lag = list(TM = 0, RR = 0)) {
  median <- model$clim_median %||% c(TM = 0, RR = 0)
  # (A list, so that the climates can be vectors: one per row.)
  level <- list(TM = country$mean_T - median[["TM"]], RR = country$mean_P - median[["RR"]])
  out <- cbind(clim_TM = level[["TM"]] + clim_shift$TM, clim_RR = level[["RR"]] + clim_shift$RR,
               lag_TM = level[["TM"]] + (lag$TM %||% 0), lag_RR = level[["RR"]] + (lag$RR %||% 0),
               cur_TM = level[["TM"]] + (cur$TM %||% 0), cur_RR = level[["RR"]] + (cur$RR %||% 0),
               inc = log(country$income) - model$inc_mean, agr = country$share - model$agr_mean,
               kap = if (is.null(model$kap_mean)) 0 else log(country$capital) - model$kap_mean,
               low = if (is.null(model$low_threshold)) 0 else as.numeric(country$income < model$low_threshold))
  if (!is.null(model$knots)) {
    for (z in names(model$knots)) {
      nl <- sweep(rcs_nonlinear(out[, z], model$knots[[z]]), 2, rcs_nonlinear(0, model$knots[[z]]))
      colnames(nl) <- paste0(z, "_nl", seq_len(ncol(nl)))
      out <- cbind(out, nl)
    }
  }
  # (Ybar / Y)^(1 - xi) of the scaled terms of T and P (Ybar: the geometric
  # mean income of the rows, exp(inc_mean)).
  xi <- model$xi %||% c(TM = 0, RR = 0)
  out <- cbind(out, idol_TM = exp(-(1 - xi[["TM"]]) * out[, "inc"]), idol_RR = exp(-(1 - xi[["RR"]]) * out[, "inc"]))
  # The exposure of the terms (R/exposure.R): the country's slopes, or the
  # smooth at its long-run climate plus the shift (so that it follows the
  # climate in the paths and projections).
  for (v in c("TM", "RR")) {
    if (!exposure_on(model$exposure, v)) next
    clim <- (if (v == "TM") country$mean_T else country$mean_P) + clim_shift[[v]]
    out <- cbind(out, exposure_values(model$exposure[[v]], v, model$exposure_table[[v]], country$GID_0,
                                      rep_len(clim, nrow(out))))
  }
  out
}

# The cumulative effect `x` (sum of the growth effects of the years so far, as
# returned by the paths) as % of the level without the shock: growth effects
# compound in the log level.
cumulative_pct <- function(model, x) 100 * (exp(x) - 1)

# The terms of a model, parsed once: for each coefficient, its arm (devTM_pos,
# ...), variable, arm suffix, the modifier it is interacted with (NA for the
# arm's own term), whether that modifier is the climate (or its spline), and
# whether the term is times (Ybar / Y)^(1 - xi) (a scaled term, ending in :idol).
# With lags, a term of lag k has every part ending in .l<k> (lag: k; arm and
# modifier without it).
structure_cache <- new.env()
model_structure <- function(model) {
  key <- model$key %||% paste(names(model$coef), collapse = ",")
  if (is.null(structure_cache[[key]])) {
    parts <- strsplit(names(model$coef), ":", fixed = TRUE)
    arm <- vapply(parts, function(p) grep("^dev", p, value = TRUE), "")
    lag <- rep(0L, length(arm))
    lagged <- grepl("\\.l[0-9]+$", arm)
    lag[lagged] <- as.integer(sub(".*\\.l", "", arm[lagged]))
    parts <- lapply(parts, function(p) sub("\\.l[0-9]+$", "", p))
    arm <- sub("\\.l[0-9]+$", "", arm)
    # The exposure channel of a term ("exh", "exc", "exw", "exd"; R/exposure.R; NA:
    # none), a factor like the income scaling, not a modifier.
    expo <- vapply(parts, function(p) {
      e <- intersect(p, EXPOSURE_COLUMNS)
      if (length(e)) e else NA_character_
    }, "")
    modifier <- vapply(parts, function(p) {
      other <- p[!grepl("^dev", p) & !p %in% c("idol", EXPOSURE_COLUMNS)]
      if (length(other)) other else NA_character_
    }, "")
    # The first-difference terms of style "both" (dev<v>_d_<arm>).
    suffix <- substring(arm, 6)
    diff <- startsWith(suffix, "_d_")
    structure_cache[[key]] <- data.frame(
      term = names(model$coef), arm = arm, variable = substr(arm, 4, 5), suffix = sub("^_d_", "_", suffix),
      diff = diff, lag = lag,
      modifier = modifier, climate = !is.na(modifier) & grepl("^(clim|lag|cur)_", modifier),
      scaled = vapply(parts, function(p) "idol" %in% p, TRUE), expo = expo,
      stringsAsFactors = FALSE)
  }
  structure_cache[[key]]
}

# The rows of the values k years before in series of `years` years stacked by
# path (NA before the first year of each path); `years` NULL: not a series
# (every lag at the same values, e.g. the response to one year's anomaly).
lag_rows <- function(n, k, years = NULL) {
  if (is.null(years) || k == 0) return(seq_len(n))
  year <- rep(seq_len(years), length.out = n)
  ifelse(year > k, seq_len(n) - k, NA_integer_)
}
# The deviations x and the modifiers (a matrix, one row per element of x) k
# years before: no deviation before the first year (its modifiers: those of
# the year itself).
lagged_inputs <- function(x, mods, k, years = NULL) {
  n <- n_inputs(x)
  at <- lag_rows(n, k, years)
  of <- function(z) ifelse(is.na(at), 0, z[at])
  list(x = if (is.list(x)) lapply(x, of) else of(x), mods = mods[ifelse(is.na(at), seq_len(n), at), , drop = FALSE])
}

# The modifiers as a matrix with one row per element of x: a vector repeated,
# or already a matrix.
as_modifier_matrix <- function(modifiers, n) {
  if (is.matrix(modifiers)) return(modifiers)
  matrix(modifiers, nrow = n, ncol = length(modifiers), byrow = TRUE,
         dimnames = list(NULL, names(modifiers)))
}

# The country's modifiers in each year of a path whose long-run climate of v is
# shifted by `shift`, whose anomalies of v this year and last year are `cur`
# and `lag`, and whose long-run income and share are `income` and `share` (one
# value per year; NULL: none, the country's), one row per year (see
# country_modifiers()).
modifier_matrix <- function(model, country, v, shift, income = NULL, share = NULL, cur = NULL, lag = NULL) {
  if (all(shift == 0) && is.null(income) && is.null(share) && is.null(cur) && is.null(lag)) {
    return(as_modifier_matrix(country_modifiers(model, country), length(shift)))
  }
  clim_shift <- list(TM = 0, RR = 0)
  clim_shift[[v]] <- shift
  anomaly <- function(a) {
    out <- list(TM = 0, RR = 0)
    if (!is.null(a)) out[[v]] <- a
    out
  }
  if (!is.null(income)) country$income <- income
  if (!is.null(share)) country$share <- share
  out <- modifier_columns(model, country, clim_shift, anomaly(cur), anomaly(lag))
  # One row per year even when the shift is a single value.
  if (nrow(out) < length(shift)) out <- out[rep(1, length(shift)), , drop = FALSE]
  out
}

# Weights on the model's coefficients of deviations x (SD, a vector) of
# variable v ("TM" or "RR"): one row per element of x. `modifiers`: a named
# vector, or a matrix with one row per element of x.
# With lags, the terms of lag k take the deviations and modifiers k years
# before in a series of `years` years (stacked by path), or with `years` NULL
# the same ones (the sum over the lags: the cumulative effect of one year).
deviation_weights <- function(model, v, x, modifiers, years = NULL) {
  st <- model_structure(model)
  mods <- as_modifier_matrix(modifiers, n_inputs(x))
  G <- matrix(0, nrow = n_inputs(x), ncol = nrow(st), dimnames = list(NULL, st$term))
  for (k in unique(st$lag[st$variable == v])) {
    at <- lagged_inputs(x, mods, k, years)
    for (j in which(st$variable == v & st$lag == k)) {
      G[, j] <- arm_value(st$suffix[[j]], term_inputs(at$x, st$diff[[j]])) *
        (if (is.na(st$modifier[[j]])) 1 else at$mods[, st$modifier[[j]]]) *
        (if (st$scaled[[j]]) at$mods[, paste0("idol_", v)] else 1) *
        (if (is.na(st$expo[[j]])) 1 else at$mods[, st$expo[[j]]])
    }
  }
  G
}

# The scale of model's deviations of variable v for a country: its SD (C, or m)
# for a model in SD, 1 for a model in absolute units (the deviation is the
# anomaly itself).
deviation_scale <- function(model, country, v) {
  if (identical(of_variable(model$units, v), "absolute")) return(1)
  if (v == "TM") country$sd_T else country$sd_P
}

linear_combination <- function(model, G) {
  estimate <- drop(G %*% model$coef)
  se <- sqrt(pmax(rowSums((G %*% model$vcov) * G), 0))
  data.frame(estimate = estimate, se = se,
             low = estimate - stats::qnorm(0.975) * se, high = estimate + stats::qnorm(0.975) * se)
}

# The capped responses ------------------------------------------------------
# Each arm a of variable v responds with
#   c_a + f_a(clim) + g_a(inc) + s_a(agr)     (clim, inc, agr centred)
# (linear, or splines that vanish at zero), where c_a is the response at the
# mean income and share without the climate part. With the cap, the
# non-climate part c_a + g_a + s_a is at most max(c_a, 0), whatever the income
# and share: income and share can bring a damaging c_a to zero but not past it,
# and cannot raise a beneficial c_a. The climate part f_a, at the country's
# historical climate and any warming beyond it, is added unchanged. With
# linear and quadratic arms, the parts are those of the side (warm, cold, wet,
# dry) at the deviation: the linear and the squared arm's summed before the cap.

# Draws of the coefficients from their estimated distribution (normal, with
# the clustered covariance), the same for a model every time they are asked.
draw_cache <- new.env()
coefficient_draws <- function(model) {
  key <- model$key %||% paste(names(model$coef), collapse = ",")
  if (is.null(draw_cache[[key]])) {
    e <- eigen(model$vcov, symmetric = TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = length(e$values))
    old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    set.seed(20261023)
    Z <- matrix(stats::rnorm(N_DRAWS * length(model$coef)), nrow = N_DRAWS)
    if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv())
    draws <- sweep(Z %*% t(A), 2, model$coef, "+")
    colnames(draws) <- names(model$coef)
    draw_cache[[key]] <- draws
  }
  draw_cache[[key]]
}

# Effect on growth of deviations x (SD) of variable v for each coefficient
# vector of `beta` (rows): a matrix with one row per element of x and one
# column per row of beta. `modifiers`: a named vector or one row per element of
# x. With `cap`, each arm's own and income and share part is capped at zero
# (see above), row by row. `paired`: beta has one row per element of x, each
# element at its own coefficients (a one-column matrix). With lags, as
# deviation_weights() (`years`), each lag capped on its own.
growth_matrix <- function(beta, model, v, x, modifiers, cap, paired = FALSE, years = NULL) {
  # The cap of variable v (one value, or c(TM, RR)).
  cap <- isTRUE(of_variable(cap, v))
  st <- model_structure(model)
  n <- n_inputs(x)
  all_mods <- as_modifier_matrix(modifiers, n)
  n_out <- if (paired) 1 else nrow(beta)
  out <- matrix(0, nrow = n, ncol = n_out)
  zero <- matrix(0, nrow = n, ncol = n_out)
  # The modifiers times the coefficients: every row of beta, or row by row.
  times <- function(m, b) if (paired) matrix(rowSums(m * b), ncol = 1) else m %*% t(b)
  # Each side's (warm/cold, wet/dry) plain terms and its terms times (Ybar /
  # Y)^(1 - xi) form a block each, capped on its own; the second block is then
  # times (Ybar / Y)^(1 - xi). With linear and quadratic arms, a block holds
  # both, and the cap applies to their sum at each x.
  # With style "both", the terms of the unit and those of its first difference
  # form blocks of their own.
  # With exposure (T), each channel's terms form blocks of their own too, then
  # times the channel's exposure.
  side <- paste0(sub("2$", "", st$suffix), ifelse(st$diff, "_d", ""), ifelse(is.na(st$expo), "", paste0(":", st$expo)))
  for (k in unique(st$lag[st$variable == v])) {
    at <- lagged_inputs(x, all_mods, k, years)
    mods <- at$mods
    of_k <- st$variable == v & st$lag == k
    for (s in unique(side[of_k])) for (scaled in unique(st$scaled[of_k & side == s])) {
      own_total <- nonclim_total <- clim_total <- zero
      block <- st[of_k & side == s & st$scaled == scaled, ]
      x_k <- term_inputs(at$x, block$diff[[1]])
      for (suffix in unique(block$suffix)) {
        terms <- block[block$suffix == suffix, ]
        a <- arm_value(suffix, x_k)
        own_term <- terms$term[is.na(terms$modifier)]
        own <- if (!length(own_term)) zero else if (paired) beta[, own_term, drop = FALSE] else {
          matrix(beta[, own_term], nrow = n, ncol = nrow(beta), byrow = TRUE)
        }
        adjustment <- terms[!is.na(terms$modifier) & !terms$climate, ]
        nonclim <- own + if (nrow(adjustment)) {
          times(mods[, adjustment$modifier, drop = FALSE], beta[, adjustment$term, drop = FALSE])
        } else 0
        climate <- terms[terms$climate, ]
        clim <- if (nrow(climate)) times(mods[, climate$modifier, drop = FALSE], beta[, climate$term, drop = FALSE]) else zero
        own_total <- own_total + a * own
        nonclim_total <- nonclim_total + a * nonclim
        clim_total <- clim_total + a * clim
      }
      if (cap) nonclim_total <- pmin(nonclim_total, pmax(own_total, 0))
      total <- nonclim_total + clim_total
      if (scaled) total <- total * mods[, paste0("idol_", v)]
      if (!is.na(block$expo[[1]])) total <- total * mods[, block$expo[[1]]]
      out <- out + total
    }
  }
  out
}

# Estimate (at the estimated coefficients) and 95% CI from the draws (rows of
# `draws`: one per element of the estimate).
capped_summary <- function(estimate, draws) {
  q <- apply(draws, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  data.frame(estimate = estimate, se = apply(draws, 1, stats::sd), low = q[1, ], high = q[2, ])
}

point_beta <- function(model) matrix(model$coef, nrow = 1, dimnames = list(NULL, names(model$coef)))

# Effect on growth (log points) of a year at each deviation x (SD) of variable v,
# whose value is the country's long-run climate plus the anomaly (x times the
# deviation's scale), last year's at the long-run climate.
response_curve <- function(model, v, x, country, cap = FALSE) {
  cap <- isTRUE(of_variable(cap, v))
  modifiers <- modifier_matrix(model, country, v, rep(0, length(x)), cur = x * deviation_scale(model, country, v))
  off <- deviation_offset(model, country, v)
  # With style "both", a year at x after a year at zero: x is also its first
  # difference.
  inputs <- if (isTRUE(model$both)) list(x = x, d = x) else x
  at_offset <- offset_inputs(inputs, off)
  if (!cap) {
    G <- deviation_weights(model, v, plus_offset(inputs, off), modifiers)
    if (off != 0) G <- G - deviation_weights(model, v, at_offset, modifiers)
    return(cbind(x = x, linear_combination(model, G)))
  }
  growth <- function(beta) {
    out <- growth_matrix(beta, model, v, plus_offset(inputs, off), modifiers, TRUE)
    if (off != 0) out <- out - growth_matrix(beta, model, v, at_offset, modifiers, TRUE)
    out
  }
  cbind(x = x, capped_summary(drop(growth(point_beta(model))), growth(coefficient_draws(model))))
}

# The anomaly path (C or m, against the country's long-run climate) of a shock
# of size `delta`: "one-off" (year 0 only), "permanent" (every year from 0),
# "gradual" (rising linearly from 0 in year 0 to delta in year `ramp_years`,
# then permanent) or "delayed" (nothing before year `ramp_years`, delta from it
# on), for the years -window..horizon.
anomaly_path <- function(type, delta, horizon, window, ramp_years = NULL) {
  pre <- if (is.finite(window)) window else 0
  years <- seq(-pre, horizon)
  a <- switch(type, "one-off" = ifelse(years == 0, delta, 0),
              "permanent" = ifelse(years >= 0, delta, 0),
              "gradual" = delta * pmin(pmax(years, 0) / ramp_years, 1),
              "delayed" = ifelse(years >= ramp_years, delta, 0))
  stats::setNames(a, years)
}

# Shift of the baseline mean in the years from `from` given the anomaly path `a`
# (named by year): none with the horizon baseline; with a trailing baseline of
# `window` years, the mean of the anomalies of the `window` years before.
baseline_shift <- function(a, window, from = 0) {
  years <- as.integer(names(a))
  out <- years[years >= from]
  if (!is.finite(window)) return(rep(0, length(out)))
  vapply(out, function(t) mean(a[years >= t - window & years <= t - 1]), 0)
}

# Deviations (SD) of the years from `from` given the anomaly path `a` (named by
# year). Horizon baseline: against the long-run mean with the long-run SD.
# Trailing baseline of `window` years: against the mean of the anomalies of the
# `window` years before, divided by the SD `sd` ("fixed") or by the trailing SD
# that the anomalies add to it ("moving": sqrt(sd^2 + variance of the anomalies
# in the window)).
deviation_path <- function(a, sd, window, sd_mode = "fixed", from = 0) {
  years <- as.integer(names(a))
  out <- years[years >= from]
  if (!is.finite(window)) return(a[years >= from] / sd)
  vapply(out, function(t) {
    past <- a[years >= t - window & years <= t - 1]
    s <- if (sd_mode == "moving") sqrt(sd^2 + mean((past - mean(past))^2)) else sd
    (a[years == t] - mean(past)) / s
  }, 0)
}

# Effect on log GDP of years 0..horizon (the sum of the growth effects up to
# each year) of an anomaly path of variable v; with its 95% CI, as % of the
# GDP without the shock. The long-run climate of the interactions is the
# country's ("fixed") or follows the baseline of the deviation ("baseline":
# the country's plus the shift of the trailing mean; the same with the horizon
# baseline). `cap`: the income and share interactions capped (see above).
# `ramp_years`: the year a "gradual" shift reaches `delta` and a "delayed" one
# starts.
gdp_path <- function(model, v, country, type, delta, horizon, sd_mode = "fixed",
                     clim_mode = "fixed", cap = FALSE, ramp_years = NULL) {
  cap <- isTRUE(of_variable(cap, v))
  window <- EXPLORER_WINDOWS[[model$window]]
  sd <- deviation_scale(model, country, v)
  if (identical(of_variable(model$units, v), "absolute")) sd_mode <- "fixed"
  a <- anomaly_path(type, delta, horizon, window, ramp_years)
  x <- if (identical(model$window, "nrk")) innovation_path(a, filter_lags(model, country, v), sd)
       else deviation_path(a, sd, window, sd_mode)
  x <- model_deviations(model, x)
  n <- n_inputs(x)
  shift <- if (clim_mode == "baseline") baseline_shift(a, window) else rep(0, n)
  # The modifiers of each year: the long-run climate of v shifted with the
  # baseline, this year's and last year's value with the anomalies (and their
  # spline terms recomputed).
  years <- as.integer(names(a))
  from <- years >= 0
  mods <- modifier_matrix(model, country, v, shift, cur = a[from], lag = c(0, a)[which(from)])
  off <- deviation_offset(model, country, v)
  at_offset <- offset_inputs(x, off)
  if (cap) {
    growth <- function(beta) {
      out <- growth_matrix(beta, model, v, plus_offset(x, off), mods, TRUE, years = n)
      if (off != 0) out <- out - growth_matrix(beta, model, v, at_offset, mods, TRUE, years = n)
      out
    }
    log_point <- accumulated(model, drop(growth(point_beta(model))))
    log_draws <- accumulated(model, growth(coefficient_draws(model)))
    lc <- capped_summary(log_point, matrix(log_draws, nrow = n))
  } else {
    G <- deviation_weights(model, v, plus_offset(x, off), mods, years = n)
    if (off != 0) G <- G - deviation_weights(model, v, at_offset, mods, years = n)
    G <- accumulated(model, matrix(G, nrow = n, dimnames = list(NULL, names(model$coef))))
    lc <- linear_combination(model, G)
  }
  data.frame(year = seq(0, horizon), deviation = term_inputs(x, FALSE), log_effect = lc$estimate,
             pct = cumulative_pct(model, lc$estimate), pct_low = cumulative_pct(model, lc$low),
             pct_high = cumulative_pct(model, lc$high))
}
