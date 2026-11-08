# GDP losses of the climate-change scenarios of build_scenarios.R for every
# country of a setting's panel, under the setting's model (R/response.R).
#
# Each country's anomaly path is the forced signal only (no year-to-year
# variability): the global temperature path (median of the FaIR ensemble) less
# its mean over the panel's horizon, times the country's CMIP6 pattern (T: the
# median slope; P: the median relative slope times the country's long-run P).
# The deviations, their baselines (horizon or trailing 20 / 30 years, with the
# history of the global path before 2025) and the SD follow the setting's
# options, as in the GDP paths of a shift; with the NRK filter the deviation is
# the innovation the anomalies add (lags from the history before 2025). The long-run characteristics of the
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
#
# Optionally the anomalies also carry natural variability: each country's
# detrended historical T and P (build_models.R), resampled by year, added to
# the forced signal from the first year of the history (PROJECTION_FROM - 31)
# and to the counterfactual (no warming: the noise alone), the same noise path
# in both; the loss is the difference. Years are resampled independently for
# each country, jointly for all countries (the same calendar year everywhere:
# the spatial correlation of T and P kept), or jointly in blocks of
# consecutive years (also their persistence). The estimate is the mean loss
# over the noise paths at the estimated coefficients; the 95% CI pairs each
# path with a draw of the coefficients (both uncertainties), and the
# variability range is that of the paths at the estimated coefficients.

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

# Natural variability ----------------------------------------------------------
VARIABILITY_CHOICES <- c("None (forced signal only)" = "none", "Years independent across countries" = "iid",
                         "Joint years (spatial correlation)" = "joint",
                         "Joint 5-year blocks (also persistence)" = "blocks")
VARIABILITY_BLOCK <- 5L
VARIABILITY_SEED <- 20261102L

# `expr` evaluated with the random seed `seed`, leaving the session's random
# state as it was.
with_seed <- function(seed, expr) {
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()))
  set.seed(seed)
  expr
}

# Indices into `n_pool` pool years of `n` noise paths of `h` years (h x n):
# single years, or blocks of `block` consecutive years.
noise_indices <- function(n, h, n_pool, block = 1L) {
  if (block == 1L) return(matrix(sample.int(n_pool, h * n, replace = TRUE), h, n))
  k <- ceiling(h / block)
  starts <- matrix(sample.int(n_pool - block + 1L, k * n, replace = TRUE), k, n)
  apply(starts, 2, function(s) as.vector(outer(seq_len(block) - 1L, s, "+"))[seq_len(h)])
}

# The noise of a country (list TM, RR: h x n matrices, C and m a year) of
# `variability` (from build_models.R): `joint`, the indices common to all
# countries (NULL: independent years, drawn with the country's own seed).
country_noise <- function(variability, iso3, h, n, joint = NULL) {
  idx <- joint %||% with_seed(VARIABILITY_SEED + sum(utf8ToInt(iso3) * c(1L, 128L, 16384L)),
                              noise_indices(n, h, length(variability$years)))
  list(TM = matrix(variability$TM[iso3, idx], h), RR = matrix(variability$RR[iso3, idx], h))
}

# Mean of the rows `window` before each row of `rows` of matrix A (columns:
# paths), by cumulative sums.
window_means <- function(A, rows, window) {
  cs <- rbind(0, cumsum_rows(A))
  (cs[rows, , drop = FALSE] - cs[rows - window, , drop = FALSE]) / window
}

# As project_country() with noise paths: the anomalies `anomalies` (forced
# signal, named by year) plus the noise (h x n matrices, the rows the years of
# the anomalies) against the noise alone. Effects on log GDP (years x paths) at
# the estimated coefficients (point, point_T, point_P) and with path k at draw
# k of `draws` (draws).
project_country_noise <- function(model, country, anomalies, noise, sd_mode, clim_mode, cap, draws,
                                  income = NULL, share = NULL) {
  window <- EXPLORER_WINDOWS[[model$window]]
  years <- as.integer(names(anomalies$TM))
  rows <- which(years >= PROJECTION_FROM)
  n <- ncol(noise$TM)
  # (Absolute deviations of a variable: no SD to move; see growth().)
  paired <- draws[rep(seq_len(n), each = length(rows)), , drop = FALSE]
  if (!is.null(income)) income <- rep(income, times = n)
  if (!is.null(share)) share <- rep(share, times = n)
  growth <- function(v, A) {
    sd <- deviation_scale(model, country, v)
    # Absolute deviations of v: no SD to move (local to this variable).
    if (identical(of_variable(model$units, v), "absolute")) sd_mode <- "fixed"
    # The baseline means M and deviations X of the rows r.
    deviations_at <- function(r) {
      if (identical(model$window, "nrk")) {
        # The innovations of the NRK filter (no baseline mean that moves).
        list(M = matrix(0, length(r), n), X = filtered_rows(A, r, filter_lags(model, country, v)) / sd)
      } else if (is.finite(window)) {
        M <- window_means(A, r, window)
        # The moving SD widens with the forced signal only: the noise is the
        # variability the country's SD already measures.
        s <- if (sd_mode == "moving") {
          a <- anomalies[[v]]
          vapply(r, function(k) sqrt(sd^2 + mean((a[k - window:1] - mean(a[k - window:1]))^2)), 0)
        } else sd
        list(M = M, X = (A[r, , drop = FALSE] - M) / s)
      } else {
        list(M = matrix(0, length(r), n), X = A[r, , drop = FALSE] / sd)
      }
    }
    at <- deviations_at(rows)
    M <- at$M
    # First differences: less the deviation of the year before (in the same
    # path, from the history before PROJECTION_FROM); the level as it is: plus
    # the country's own level (deviation_offset()), in both paths compared.
    # Style "both": the deviations and their first differences.
    X <- as.vector(at$X)
    if (isTRUE(model$diff) || isTRUE(model$both)) {
      D <- X - as.vector(deviations_at(rows - 1L)$X)
      X <- if (isTRUE(model$both)) list(x = X, d = D) else D
    }
    X <- plus_offset(X, deviation_offset(model, country, v))
    shift <- if (clim_mode == "baseline") as.vector(M) else rep(0, n_inputs(X))
    # This year's and last year's value: the anomalies with the noise.
    mods <- modifier_matrix(model, country, v, shift, income, share,
                            cur = as.vector(A[rows, , drop = FALSE]), lag = as.vector(A[rows - 1, , drop = FALSE]))
    list(point = matrix(growth_matrix(point_beta(model), model, v, X, mods, cap, years = length(rows)), length(rows)),
         draws = matrix(growth_matrix(paired, model, v, X, mods, cap, paired = TRUE, years = length(rows)), length(rows)))
  }
  out <- list()
  for (v in c("TM", "RR")) {
    with_signal <- growth(v, anomalies[[v]] + noise[[v]])
    alone <- growth(v, noise[[v]])
    out[[v]] <- list(point = with_signal$point - alone$point, draws = with_signal$draws - alone$draws)
  }
  list(point = accumulated(model, out$TM$point + out$RR$point), point_T = accumulated(model, out$TM$point),
       point_P = accumulated(model, out$RR$point), draws = accumulated(model, out$TM$draws + out$RR$draws))
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
    # Absolute deviations of v: no SD to move (this variable only).
    sd_mode_v <- if (identical(of_variable(model$units, v), "absolute")) "fixed" else sd_mode
    # First differences: from the year before PROJECTION_FROM (its deviation
    # from the history), each year less the year before.
    # (Style "both": the deviations and their first differences.)
    differenced <- isTRUE(model$diff) || isTRUE(model$both)
    from <- PROJECTION_FROM - differenced
    x <- if (identical(model$window, "nrk")) {
      innovation_path(a, filter_lags(model, country, v), sd, from = from)
    } else deviation_path(a, sd, window, sd_mode_v, from = from)
    if (isTRUE(model$both)) x <- list(x = x[-1], d = diff(x)) else if (differenced) x <- diff(x)
    n <- n_inputs(x)
    shift <- if (clim_mode == "baseline") baseline_shift(a, window, PROJECTION_FROM) else rep(0, n)
    rows <- which(as.integer(names(a)) >= PROJECTION_FROM)
    mods <- modifier_matrix(model, country, v, shift, income, share, cur = a[rows], lag = a[rows - 1])
    # The level as it is: the effects of the country's level plus the
    # deviations less those of its level alone.
    off <- deviation_offset(model, country, v)
    at_offset <- offset_inputs(x, off)
    growth <- function(beta) {
      out <- growth_matrix(beta, model, v, plus_offset(x, off), mods, cap, years = n)
      if (off != 0) out <- out - growth_matrix(beta, model, v, at_offset, mods, cap, years = n)
      out
    }
    growth_point[[v]] <- drop(growth(point_beta(model)))
    growth_draws <- growth_draws + growth(draws)
  }
  list(point = accumulated(model, growth_point$TM + growth_point$RR), point_T = accumulated(model, growth_point$TM),
       point_P = accumulated(model, growth_point$RR), draws = accumulated(model, growth_draws))
}

# Growth (% a year; NA in the first year) of world GDP `world_gdp` (by year)
# less the global loss `loss` (% of it, by year: a vector, or a matrix with a
# column per draw or path).
world_growth <- function(world_gdp, loss) {
  level <- world_gdp * (1 + as.matrix(loss) / 100)
  n <- nrow(level)
  out <- rbind(NA, 100 * (level[-1, , drop = FALSE] / level[-n, , drop = FALSE] - 1))
  if (is.matrix(loss)) out else drop(out)
}

# The losses of the warming of forcing scenario `scenario` with the economic
# development of `ssp` (NULL: the scenario's own SSP) under setting `spec`
# (from app.R): the global loss by year (estimate and 95% CI, % of GDP), each
# country's loss by year (% of GDP, at the estimated coefficients; also the
# losses of the T and of the P deviations, each on its own), the share of the
# SSP's world GDP the countries cover, and how many follow the SSP in income
# and share. `variability`: the natural variability of the setting's panel
# (build_models.R; NULL: none), resampled as `noise` (a value of
# VARIABILITY_CHOICES) in `n_paths` paths.
project_setting <- function(spec, scen_data, scenario, ssp = NULL, gmt = NULL,
                            variability = NULL, noise = "none", n_paths = 500L) {
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
  # With noise, the countries with the whole historical climate.
  with_noise <- noise != "none" && !is.null(variability)
  if (with_noise) countries <- dplyr::filter(countries, GID_0 %in% rownames(variability$TM))
  weights <- gdp_wide[countries$GID_0, as.character(years), drop = FALSE]
  coverage <- sum(weights[, 1]) / sum(gdp_wide[, 1], na.rm = TRUE)
  # The SSP's GDP of the countries (no climate change), for the growth rates.
  world_gdp <- colSums(weights)
  weights <- sweep(weights, 2, colSums(weights), "/")

  n_draws <- if (with_noise) n_paths else N_PROJECTION_DRAWS
  draws <- coefficient_draws(model)[seq_len(n_draws), , drop = FALSE]
  global_point <- numeric(length(years))
  global_draws <- matrix(0, nrow = length(years), ncol = nrow(draws))
  # The global loss of each noise path at the estimated coefficients.
  global_paths <- if (with_noise) matrix(0, nrow = length(years), ncol = n_paths)
  joint <- if (noise %in% c("joint", "blocks")) {
    with_seed(VARIABILITY_SEED, noise_indices(n_paths, length(history), length(variability$years),
                                              if (noise == "blocks") VARIABILITY_BLOCK else 1L))
  }
  country_loss <- matrix(NA_real_, nrow = nrow(countries), ncol = length(years))
  country_loss_T <- country_loss_P <- country_loss
  followed <- c(income = 0, share = 0)
  ssp_gdppc <- scen_data$ssp$gdppc[scen_data$ssp$gdppc$ssp == ssp, ]
  ssp_share <- scen_data$ssp$share[scen_data$ssp$share$ssp == ssp, ]
  income_history <- scen_data$history$income[[model$gdp]]
  share_history <- scen_data$history$share
  for (i in seq_len(nrow(countries))) {
    cn <- countries[i, ]
    country <- list(GID_0 = cn$GID_0, mean_T = cn$mean_T, sd_T = cn$sd_T, mean_P = cn$mean_P, sd_P = cn$sd_P,
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
    if (with_noise) {
      r <- project_country_noise(model, country, anomalies,
                                 country_noise(variability, cn$GID_0, length(history), n_paths, joint),
                                 spec$sd_mode, spec$clim_mode, spec$cap, draws, income, share)
      paths <- cumulative_pct(model, r$point)
      country_loss[i, ] <- rowMeans(paths)
      country_loss_T[i, ] <- rowMeans(cumulative_pct(model, r$point_T))
      country_loss_P[i, ] <- rowMeans(cumulative_pct(model, r$point_P))
      global_paths <- global_paths + weights[i, ] * paths
    } else {
      r <- project_country(model, country, anomalies, spec$sd_mode, spec$clim_mode, spec$cap, draws,
                           income, share)
      country_loss[i, ] <- cumulative_pct(model, r$point)
      country_loss_T[i, ] <- cumulative_pct(model, r$point_T)
      country_loss_P[i, ] <- cumulative_pct(model, r$point_P)
    }
    global_point <- global_point + weights[i, ] * country_loss[i, ]
    global_draws <- global_draws + weights[i, ] * cumulative_pct(model, r$draws)
  }
  q <- apply(global_draws, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  global <- data.frame(year = years, estimate = global_point, low = q[1, ], high = q[2, ])
  # Growth of the countries' GDP (% a year, from the second year): without
  # climate change (the SSP's), and with it, whose level is the SSP's times one
  # plus the global loss (the loss is the GDP-weighted mean of the countries'),
  # draw by draw for the CI.
  global$baseline_growth <- c(NA, 100 * (world_gdp[-1] / world_gdp[-length(years)] - 1))
  global$growth <- world_growth(world_gdp, global_point)
  qg <- apply(world_growth(world_gdp, global_draws), 1, stats::quantile, probs = c(0.025, 0.975),
              names = FALSE, na.rm = TRUE)
  global$growth_low <- qg[1, ]
  global$growth_high <- qg[2, ]
  if (with_noise) {
    qv <- apply(global_paths, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
    global$noise_low <- qv[1, ]
    global$noise_high <- qv[2, ]
    qv <- apply(world_growth(world_gdp, global_paths), 1, stats::quantile, probs = c(0.025, 0.975),
                names = FALSE, na.rm = TRUE)
    global$growth_noise_low <- qv[1, ]
    global$growth_noise_high <- qv[2, ]
  }
  list(
    global = global, noise = if (with_noise) noise else "none", n_paths = if (with_noise) n_paths,
    countries = data.frame(iso3 = rep(countries$GID_0, times = length(years)),
                           country = rep(countries$country, times = length(years)),
                           year = rep(years, each = nrow(countries)),
                           loss = as.vector(country_loss), loss_T = as.vector(country_loss_T),
                           loss_P = as.vector(country_loss_P)),
    coverage = coverage, n = nrow(countries), follow = follow, followed = followed, ssp = ssp
  )
}

# Exceedances of natural variability -------------------------------------------
# How far a warming path takes each country's annual T beyond its natural
# variability: the forced deviation of the country (its CMIP6 slope times the
# deviation of global warming from a reference climate: the trailing mean of
# the years before, or the mean over the horizon) against the SD of its
# detrended historical T (build_models.R).

# The share of each country's area (`world`: sf borders with iso3) in latitude
# bands of `width` degrees: data.frame(iso3, lat (the band's centre), share).
# Areas in squared degrees times the cosine of the band's centre.
latitude_shares <- function(world, width = 5) {
  w <- sf::st_transform(world, 4326)
  w <- w[grepl("^[A-Z]{3}$", w$iso3), ]
  w <- sf::st_make_valid(sf::st_set_crs(w, NA))
  edges <- seq(-90, 90, by = width)
  out <- do.call(rbind, lapply(seq_len(length(edges) - 1), function(b) {
    band <- sf::st_as_sfc(sf::st_bbox(c(xmin = -180, ymin = edges[b], xmax = 180, ymax = edges[b + 1])))
    a <- suppressWarnings(sf::st_intersection(w, band))
    if (!nrow(a)) return(NULL)
    lat <- edges[b] + width / 2
    data.frame(iso3 = a$iso3, lat = lat, area = as.numeric(sf::st_area(a)) * cos(lat * pi / 180))
  }))
  out <- stats::aggregate(area ~ iso3 + lat, out, sum)
  out <- out[out$area > 0, ]
  out$share <- out$area / stats::ave(out$area, out$iso3, FUN = sum)
  out[, c("iso3", "lat", "share")]
}

# Global warming `gmt` (named by year) in `years`, its reference climate (the
# mean of the `window` years before each year, or with `window` NA the mean
# over the years `horizon`) and its deviation from it.
gmt_deviation <- function(gmt, years, window, horizon) {
  reference <- if (is.na(window)) {
    rep(mean(gmt[as.character(seq(horizon[[1]], horizon[[2]]))]), length(years))
  } else {
    trailing_mean(gmt, years, window)
  }
  warming <- unname(gmt[as.character(years)])
  list(warming = warming, reference = reference, deviation = warming - reference)
}

# The natural variability of T of the countries of `variability` with a CMIP6
# pattern (`pattern`) and a population (`population`, named by iso3): each
# country's SD of detrended annual T (sd) and slope (slope), named by iso3, and
# the SD of the population-weighted mean of the countries' detrended T in
# degrees of global warming (gmt: the SD over the mean slope, so that global
# warming of gmt gives a population-weighted country warming of one SD).
natural_sigma <- function(variability, pattern, population) {
  slope <- stats::setNames(pattern$t_slope, pattern$iso3)
  iso <- rownames(variability$TM)
  iso <- iso[iso %in% names(slope)[is.finite(slope)] & iso %in% names(population)[is.finite(population) & population > 0]]
  w <- population[iso] / sum(population[iso])
  TM <- variability$TM[iso, , drop = FALSE]
  list(sd = apply(TM, 1, stats::sd), slope = slope[iso],
       gmt = stats::sd(colSums(w * TM)) / sum(w * slope[iso]))
}

# The share of the population (%) exposed to deviations of T beyond `k` SD by
# latitude band and year, for global warming path `gmt`: each country's
# probability of such a year, from a normal distribution with its natural SD
# around its forced deviation (`sign`: "warm", above +k SD, or "both", beyond
# +-k SD), weighted by its population (`population`: iso3 x years) spread over
# the latitude bands by area (`shares`: latitude_shares()).
exposure_by_latitude <- function(gmt, years, window, horizon, sigma, population, shares, k = 2.5,
                                 sign = "warm") {
  dev <- gmt_deviation(gmt, years, window, horizon)$deviation
  iso <- intersect(intersect(names(sigma$sd), rownames(population)), shares$iso3)
  z <- outer(sigma$slope[iso] / sigma$sd[iso], dev)
  p <- 1 - stats::pnorm(k - z)
  if (sign == "both") p <- p + stats::pnorm(-k - z)
  pop <- population[iso, as.character(years), drop = FALSE]
  s <- shares[shares$iso3 %in% iso, ]
  S <- stats::xtabs(share ~ iso3 + lat, s)[iso, , drop = FALSE]
  S <- matrix(S, nrow = length(iso), dimnames = list(iso, colnames(S)))
  exposed <- crossprod(S, pop * p)
  total <- crossprod(S, pop)
  data.frame(lat = rep(as.numeric(colnames(S)), times = length(years)), year = rep(years, each = ncol(S)),
             exposed = as.vector(100 * exposed / total), population = as.vector(total))
}
