# GDP losses of the climate-change scenarios of build_scenarios.R for every
# country of a setting's panel, under the setting's model (R/response.R).
#
# Each country's anomaly path is the forced signal only (no year-to-year
# variability): the global temperature path (median of the FaIR ensemble) less
# its mean over the panel's horizon, times the country's CMIP6 pattern (T: the
# median slope; P: the median relative slope times the country's long-run P).
# The deviations, their baselines (horizon or trailing 20 / 30 years, with the
# history of the global path before 2025) and the SD follow the setting's
# options, as in the GDP paths of a shift. The long-run characteristics of the
# interactions are the country's (fixed) or, with a trailing baseline, follow
# it: the long-run climate shifts with the trailing mean of the anomalies, and
# income and the agricultural share are the trailing means of the country's
# own historical series continued with the SSP of the scenario (the SSP's
# relative change from the last historical year, chained onto it); for a model
# estimated with trailing income (and share), over the window it was estimated
# with, also with the horizon baseline.
# Losses accumulate from PROJECTION_FROM (the counterfactual: no deviation from
# then on), in % of the GDP without climate change; the global loss is the
# average of the countries' losses weighted by the SSP's GDP of each year. Its
# 95% CI comes from draws of the coefficients, the same draws for every country.

PROJECTION_FROM <- 2025L
PROJECTION_TO <- 2100L
# Draws of the coefficients for the CI of the global loss (the first of the
# model's draws of coefficient_draws()).
N_PROJECTION_DRAWS <- 1000L

# Global warming with perfect geoengineering (a climate with the same global
# temperature is the same with or without it): the background path `background`
# (named by year) until PROJECTION_FROM; from then on held at `held` (or the
# background where it is lower); from the termination year `termination` (NULL:
# never) the gap to the background closes with characteristic time `tau`:
#   G(t) = B(t) - max(B(T) - H, 0) exp(-(t - T) / tau).
termination_path <- function(background, held, termination = NULL, tau = 4) {
  years <- as.integer(names(background))
  out <- background
  geo <- years >= PROJECTION_FROM
  if (!is.null(termination)) geo <- geo & years < termination
  out[geo] <- pmin(background[geo], held)
  if (!is.null(termination)) {
    gap <- max(background[[as.character(termination)]] - held, 0)
    after <- years >= termination
    out[after] <- background[after] - gap * exp(-(years[after] - termination) / tau)
  }
  out
}

# Cumulative sums down the rows of a matrix (each column a path).
cumsum_rows <- function(m) {
  if (nrow(m) > 1) for (i in 2:nrow(m)) m[i, ] <- m[i - 1, ] + m[i, ]
  m
}

# A country's series in the model's definition: the historical values `history`
# (year, value) to their last year, then that value times the SSP's relative
# change (`ssp`: year, value), named by year.
chained_series <- function(history, ssp) {
  history <- history[is.finite(history$value), ]
  base <- max(history$year)
  ssp_base <- ssp$value[ssp$year == base]
  future <- ssp[ssp$year > base, ]
  if (!length(ssp_base) || !is.finite(ssp_base) || ssp_base <= 0) future <- future[0, ]
  stats::setNames(c(history$value, history$value[history$year == base] * future$value / ssp_base),
                  c(history$year, future$year))
}

# Mean of `series` (named by year) over the `window` years before each year of
# `years` (the years with a value).
trailing_mean <- function(series, years, window) {
  y <- as.integer(names(series))
  vapply(years, function(t) {
    v <- series[y >= t - window & y <= t - 1]
    if (length(v)) mean(v) else NA_real_
  }, 0)
}

# Cumulative effect on log GDP of years PROJECTION_FROM..PROJECTION_TO of the
# anomaly paths `anomalies` (list TM, RR, named by year with the history the
# trailing baselines need), with the long-run income and share of each year
# (`income`, `share`; NULL: the country's), at the estimated coefficients
# (point, and its parts of T and of P) and at each draw (years x draws).
project_country <- function(model, country, anomalies, sd_mode, clim_mode, cap, draws,
                            income = NULL, share = NULL) {
  window <- EXPLORER_WINDOWS[[model$window]]
  growth_point <- list()
  growth_draws <- 0
  for (v in c("TM", "RR")) {
    a <- anomalies[[v]]
    sd <- deviation_scale(model, country, v)
    if (identical(model$units, "absolute")) sd_mode <- "fixed"
    x <- deviation_path(a, sd, window, sd_mode, from = PROJECTION_FROM)
    shift <- if (clim_mode == "baseline") baseline_shift(a, window, PROJECTION_FROM) else rep(0, length(x))
    mods <- modifier_matrix(model, country, v, shift, income, share)
    growth_point[[v]] <- drop(growth_matrix(point_beta(model), model, v, x, mods, cap))
    growth_draws <- growth_draws + growth_matrix(draws, model, v, x, mods, cap)
  }
  list(point = cumsum(growth_point$TM + growth_point$RR), point_T = cumsum(growth_point$TM),
       point_P = cumsum(growth_point$RR), draws = cumsum_rows(growth_draws))
}

# The losses of the warming of forcing scenario `scenario` with the economic
# development of `ssp` (NULL: the scenario's own SSP) under setting `spec`
# (from app.R): the global loss by year (estimate and 95% CI, % of GDP), each
# country's loss by year (% of GDP, at the estimated coefficients; also the
# losses of the T and of the P deviations, each on its own), the share of the
# SSP's world GDP the countries cover, and how many follow the SSP in income
# and share.
project_setting <- function(spec, scen_data, scenario, ssp = NULL, gmt = NULL) {
  model <- spec$model
  ssp <- ssp %||% scen_data$ssp$of_scenario[[scenario]]
  # The global warming path: the scenario's median, or the path `gmt` (named by
  # year), e.g. one with geoengineering (termination_path()).
  if (is.null(gmt)) {
    g <- scen_data$gmt[scen_data$gmt$scenario == scenario, ]
    gmt <- stats::setNames(g$median, g$year)
  }
  horizon <- model$horizon %||% c(1970L, 2023L)
  reference <- mean(gmt[as.character(seq(horizon[[1]], horizon[[2]]))])
  history <- seq(PROJECTION_FROM - 31L, PROJECTION_TO)
  warming <- gmt[as.character(history)] - reference
  years <- seq(PROJECTION_FROM, PROJECTION_TO)
  window <- EXPLORER_WINDOWS[[model$window]]
  # Income and share follow the SSP as trailing means: over the windows they
  # were estimated with (trailing moderators), otherwise over the baseline's.
  inc_window <- model$inc_window %||% window
  share_window <- model$share_window %||% window
  # spec$follow_econ = FALSE: the climate follows the baseline, but income and
  # share stay at the country's values (default: they follow too).
  follow_econ <- !isFALSE(spec$follow_econ)
  follow_inc <- spec$clim_mode == "baseline" && follow_econ && is.finite(inc_window)
  follow_share <- spec$clim_mode == "baseline" && follow_econ && is.finite(share_window)
  follow <- follow_inc || follow_share

  ssp_gdp <- scen_data$ssp$gdp[scen_data$ssp$gdp$ssp == ssp & scen_data$ssp$gdp$year %in% years, ]
  gdp_wide <- tapply(ssp_gdp$gdp, list(ssp_gdp$iso3, ssp_gdp$year), sum)
  countries <- spec$units |>
    dplyr::inner_join(scen_data$pattern[[model$weight %||% "population"]], by = c("GID_0" = "iso3")) |>
    dplyr::filter(is.finite(t_slope), is.finite(p_rel_slope),
                  GID_0 %in% rownames(gdp_wide)[stats::complete.cases(gdp_wide)])
  # A model with capital needs the country's (PWT; held at its long-run value).
  if (!is.null(model$kap_mean)) countries <- dplyr::filter(countries, is.finite(capital))
  weights <- gdp_wide[countries$GID_0, as.character(years), drop = FALSE]
  coverage <- sum(weights[, 1]) / sum(gdp_wide[, 1], na.rm = TRUE)
  weights <- sweep(weights, 2, colSums(weights), "/")

  draws <- coefficient_draws(model)[seq_len(N_PROJECTION_DRAWS), , drop = FALSE]
  global_point <- numeric(length(years))
  global_draws <- matrix(0, nrow = length(years), ncol = nrow(draws))
  country_loss <- matrix(NA_real_, nrow = nrow(countries), ncol = length(years))
  country_loss_T <- country_loss_P <- country_loss
  followed <- c(income = 0, share = 0)
  ssp_gdppc <- scen_data$ssp$gdppc[scen_data$ssp$gdppc$ssp == ssp, ]
  ssp_share <- scen_data$ssp$share[scen_data$ssp$share$ssp == ssp, ]
  income_history <- scen_data$history$income[[model$gdp]]
  share_history <- scen_data$history$share
  for (i in seq_len(nrow(countries))) {
    cn <- countries[i, ]
    country <- list(mean_T = cn$mean_T, sd_T = cn$sd_T, mean_P = cn$mean_P, sd_P = cn$sd_P,
                    income = cn$income, share = cn$share, capital = cn$capital)
    anomalies <- list(TM = cn$t_slope * warming, RR = cn$p_rel_slope * cn$mean_P * warming)
    income <- share <- NULL
    if (follow_inc) {
      h <- income_history[income_history$iso3 == cn$GID_0, ]
      s <- ssp_gdppc[ssp_gdppc$iso3 == cn$GID_0, ]
      if (nrow(h) && nrow(s)) {
        path <- trailing_mean(chained_series(data.frame(year = h$year, value = h$gdppc),
                                             data.frame(year = s$year, value = s$gdppc)), years, inc_window)
        if (all(is.finite(path))) {
          income <- path
          followed[["income"]] <- followed[["income"]] + 1
        }
      }
    }
    if (follow_share) {
      h <- share_history[share_history$iso3 == cn$GID_0, ]
      s <- ssp_share[ssp_share$iso3 == cn$GID_0, ]
      if (nrow(h) && nrow(s)) {
        path <- trailing_mean(chained_series(data.frame(year = h$year, value = h$share),
                                             data.frame(year = s$year, value = s$share)), years, share_window)
        if (all(is.finite(path))) {
          share <- path
          followed[["share"]] <- followed[["share"]] + 1
        }
      }
    }
    r <- project_country(model, country, anomalies, spec$sd_mode, spec$clim_mode, spec$cap, draws,
                         income, share)
    country_loss[i, ] <- cumulative_pct(model, r$point)
    country_loss_T[i, ] <- cumulative_pct(model, r$point_T)
    country_loss_P[i, ] <- cumulative_pct(model, r$point_P)
    global_point <- global_point + weights[i, ] * country_loss[i, ]
    global_draws <- global_draws + weights[i, ] * cumulative_pct(model, r$draws)
  }
  q <- apply(global_draws, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  list(
    global = data.frame(year = years, estimate = global_point, low = q[1, ], high = q[2, ]),
    countries = data.frame(iso3 = rep(countries$GID_0, times = length(years)),
                           country = rep(countries$country, times = length(years)),
                           year = rep(years, each = nrow(countries)),
                           loss = as.vector(country_loss), loss_T = as.vector(country_loss_T),
                           loss_P = as.vector(country_loss_P)),
    coverage = coverage, n = nrow(countries), follow = follow, followed = followed, ssp = ssp
  )
}
