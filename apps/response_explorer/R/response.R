# Responses and GDP paths of a fitted model of build_models.R for a country
# with given characteristics. Without the cap every quantity is a linear
# combination of the model's coefficients, so its standard error is exact
# (delta method with the clustered covariance) and its 95% CI normal. With the
# cap the response is not linear in the coefficients, so the CI comes from
# draws of the coefficients (percentiles), capped draw by draw.

EXPLORER_WINDOWS <- c(all = Inf, lag_20 = 20, lag_30 = 30)
N_DRAWS <- 2000

# Value of an arm term (column suffix) at deviation x (SD).
arm_value <- function(suffix, x) {
  switch(suffix,
         "_pos" = pmax(x, 0), "_neg" = pmax(-x, 0),
         "_pos2" = pmax(x, 0)^2, "_neg2" = pmax(-x, 0)^2,
         stop("Unknown arm ", suffix))
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
# (PWT) less their means; `clim_shift` (C, m) is added to the long-run climate.
# For a spline model, also the non-linear terms <z>_nl<j>, less their value at
# zero (as build_models.R).
country_modifiers <- function(model, country, clim_shift = c(TM = 0, RR = 0)) {
  median <- model$clim_median %||% c(TM = 0, RR = 0)
  out <- c(clim_TM = country$mean_T + clim_shift[["TM"]] - median[["TM"]],
           clim_RR = country$mean_P + clim_shift[["RR"]] - median[["RR"]],
           inc = log(country$income) - model$inc_mean, agr = country$share - model$agr_mean,
           kap = if (is.null(model$kap_mean)) 0 else log(country$capital) - model$kap_mean)
  if (identical(model$shape, "splines")) {
    for (z in names(model$knots)) {
      nl <- rcs_nonlinear(out[[z]], model$knots[[z]]) - rcs_nonlinear(0, model$knots[[z]])
      out <- c(out, stats::setNames(nl[1, ], paste0(z, "_nl", seq_len(ncol(nl)))))
    }
  }
  # (Ybar / Y)^(1 - xi) of the scaled terms of T and P (Ybar: the geometric
  # mean income of the rows, exp(inc_mean)).
  xi <- model$xi %||% c(TM = 0, RR = 0)
  c(out, idol_TM = exp(-(1 - xi[["TM"]]) * out[["inc"]]), idol_RR = exp(-(1 - xi[["RR"]]) * out[["inc"]]))
}

# The cumulative effect `x` (sum of the growth effects of the years so far, as
# returned by the paths) as % of the level without the shock: growth effects
# compound in the log level.
cumulative_pct <- function(model, x) 100 * (exp(x) - 1)

# The terms of a model, parsed once: for each coefficient, its arm (devTM_pos,
# ...), variable, arm suffix, the modifier it is interacted with (NA for the
# arm's own term), whether that modifier is the climate (or its spline), and
# whether the term is times (Ybar / Y)^(1 - xi) (a scaled term, ending in :idol).
structure_cache <- new.env()
model_structure <- function(model) {
  key <- model$key %||% paste(names(model$coef), collapse = ",")
  if (is.null(structure_cache[[key]])) {
    parts <- strsplit(names(model$coef), ":", fixed = TRUE)
    arm <- vapply(parts, function(p) grep("^dev", p, value = TRUE), "")
    modifier <- vapply(parts, function(p) {
      other <- p[!grepl("^dev", p) & p != "idol"]
      if (length(other)) other else NA_character_
    }, "")
    structure_cache[[key]] <- data.frame(
      term = names(model$coef), arm = arm, variable = substr(arm, 4, 5), suffix = substring(arm, 6),
      modifier = modifier, climate = !is.na(modifier) & startsWith(modifier, "clim_"),
      scaled = vapply(parts, function(p) "idol" %in% p, TRUE),
      stringsAsFactors = FALSE)
  }
  structure_cache[[key]]
}

# The modifiers as a matrix with one row per element of x: a vector repeated,
# or already a matrix.
as_modifier_matrix <- function(modifiers, n) {
  if (is.matrix(modifiers)) return(modifiers)
  matrix(modifiers, nrow = n, ncol = length(modifiers), byrow = TRUE,
         dimnames = list(NULL, names(modifiers)))
}

# The country's modifiers in each year of a path whose long-run climate of v is
# shifted by `shift`, and whose long-run income and share are `income` and
# `share` (one value per year; NULL: the country's), one row per year (see
# country_modifiers()).
modifier_matrix <- function(model, country, v, shift, income = NULL, share = NULL) {
  if (all(shift == 0) && is.null(income) && is.null(share)) {
    return(as_modifier_matrix(country_modifiers(model, country), length(shift)))
  }
  do.call(rbind, lapply(seq_along(shift), function(i) {
    clim_shift <- c(TM = 0, RR = 0)
    clim_shift[[v]] <- shift[[i]]
    if (!is.null(income)) country$income <- income[[i]]
    if (!is.null(share)) country$share <- share[[i]]
    country_modifiers(model, country, clim_shift)
  }))
}

# Weights on the model's coefficients of deviations x (SD, a vector) of
# variable v ("TM" or "RR"): one row per element of x. `modifiers`: a named
# vector, or a matrix with one row per element of x.
deviation_weights <- function(model, v, x, modifiers) {
  st <- model_structure(model)
  mods <- as_modifier_matrix(modifiers, length(x))
  G <- matrix(0, nrow = length(x), ncol = nrow(st), dimnames = list(NULL, st$term))
  for (j in which(st$variable == v)) {
    G[, j] <- arm_value(st$suffix[[j]], x) * (if (is.na(st$modifier[[j]])) 1 else mods[, st$modifier[[j]]]) *
      (if (st$scaled[[j]]) mods[, paste0("idol_", v)] else 1)
  }
  G
}

# The scale of model's deviations of variable v for a country: its SD (C, or m)
# for a model in SD, 1 for a model in absolute units (the deviation is the
# anomaly itself).
deviation_scale <- function(model, country, v) {
  if (identical(model$units, "absolute")) return(1)
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
# historical climate and any warming beyond it, is added unchanged.

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
# (see above), row by row.
growth_matrix <- function(beta, model, v, x, modifiers, cap) {
  st <- model_structure(model)
  mods <- as_modifier_matrix(modifiers, length(x))
  out <- matrix(0, nrow = length(x), ncol = nrow(beta))
  # Each arm's plain terms and its terms times (Ybar / Y)^(1 - xi) form a block
  # each, capped on its own; the second block is then times (Ybar / Y)^(1 - xi).
  for (arm in unique(st$arm[st$variable == v])) for (scaled in unique(st$scaled[st$arm == arm])) {
    terms <- st[st$arm == arm & st$scaled == scaled, ]
    own_term <- terms$term[is.na(terms$modifier)]
    own <- if (length(own_term)) matrix(beta[, own_term], nrow = length(x), ncol = nrow(beta), byrow = TRUE) else 0
    adjustment <- terms[!is.na(terms$modifier) & !terms$climate, ]
    total <- own + if (nrow(adjustment)) {
      mods[, adjustment$modifier, drop = FALSE] %*% t(beta[, adjustment$term, drop = FALSE])
    } else 0
    if (cap) total <- pmin(total, pmax(own, 0))
    climate <- terms[terms$climate, ]
    if (nrow(climate)) {
      total <- total + mods[, climate$modifier, drop = FALSE] %*% t(beta[, climate$term, drop = FALSE])
    }
    if (scaled) total <- total * mods[, paste0("idol_", v)]
    out <- out + arm_value(terms$suffix[[1]], x) * total
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

# Effect on growth (log points) of a year at each deviation x (SD) of variable v.
response_curve <- function(model, v, x, country, cap = FALSE) {
  modifiers <- country_modifiers(model, country)
  if (!cap) return(cbind(x = x, linear_combination(model, deviation_weights(model, v, x, modifiers))))
  cbind(x = x, capped_summary(
    drop(growth_matrix(point_beta(model), model, v, x, modifiers, TRUE)),
    growth_matrix(coefficient_draws(model), model, v, x, modifiers, TRUE)))
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
  window <- EXPLORER_WINDOWS[[model$window]]
  sd <- deviation_scale(model, country, v)
  if (identical(model$units, "absolute")) sd_mode <- "fixed"
  a <- anomaly_path(type, delta, horizon, window, ramp_years)
  x <- deviation_path(a, sd, window, sd_mode)
  shift <- if (clim_mode == "baseline") baseline_shift(a, window) else rep(0, length(x))
  # The modifiers of each year: the long-run climate of v shifted with the
  # baseline (and its spline terms recomputed).
  mods <- modifier_matrix(model, country, v, shift)
  if (cap) {
    log_point <- cumsum(drop(growth_matrix(point_beta(model), model, v, x, mods, TRUE)))
    log_draws <- apply(growth_matrix(coefficient_draws(model), model, v, x, mods, TRUE), 2, cumsum)
    lc <- capped_summary(log_point, matrix(log_draws, nrow = length(x)))
  } else {
    G <- deviation_weights(model, v, x, mods)
    G <- matrix(apply(G, 2, cumsum), nrow = length(x), dimnames = list(NULL, names(model$coef)))
    lc <- linear_combination(model, G)
  }
  data.frame(year = seq(0, horizon), deviation = x, log_effect = lc$estimate,
             pct = cumulative_pct(model, lc$estimate), pct_low = cumulative_pct(model, lc$low),
             pct_high = cumulative_pct(model, lc$high))
}
