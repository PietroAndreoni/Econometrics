# Builds the estimation data of the response explorer (apps/response_explorer/
# app.R), which fits its models on demand (R/fit.R), and saves it to
# results/response_explorer/panels.rds.
#
# The models are those of "Income or the agricultural share" in
# weather_arms_robustness.Rmd: growth on the demeaned deviations of T and P, cut
# into arms (linear or quadratic) and interacted with the unit's climate,
# income and agricultural share, with the controls and year and unit effects,
# clustered by country. The options:
#
#   GDP data      PWT, WB
#   climate data  ERA5, CRU TS, UDelaware
#   weights       population density (2000) or area (plain area average) in the
#                 climate of each country
#   units         deviations in SD of the country's variability, or in absolute
#                 units (C and metres, demeaned in the same way)
#   regression unit (the baseline): the level of T and P (in SD units over
#                 the country's horizon SD), the deviation from the horizon
#                 mean, from the trailing 20 or 30 year mean (mean and SD of
#                 the deviation), or the innovation of the Nath-Ramey-Klenow
#                 AR(3) filter (see nrk_filter())
#   style         the regression unit itself, or its first difference (this
#                 year's less last year's, in the same units: with the level,
#                 Kalkuhl and Wenz), columns ending in "_d"
#   fixed effects year and unit effects; plus unit linear trends; plus unit
#                 linear and quadratic trends (the notebook's). The deviations
#                 are demeaned in the same way before their arms are cut: each
#                 unit's mean, linear trend, or linear and quadratic trend
#                 removed, fitted from the start year minus MAX_LAGS to the last
#                 year of the panel; the annual T and P are kept for the
#                 deviations as they are (symmetric models), and the filter's
#                 innovations are never demeaned
#   moderators    fixed: the horizon mean of the climate, long-run income and
#                 long-run share; trailing: the climate of the deviation's own
#                 baseline (the trailing 20 or 30 year mean; the horizon mean
#                 with the horizon baseline) and the mean income of the
#                 TRAILING_INCOME_YEARS years before, with their level terms;
#                 with WB GDP, optionally also the mean share of the
#                 TRAILING_SHARE_YEARS years before
#   arms, interactions (for T and P, any subset of climate, last year's and
#   this year's value of the variable, income, share and, with PWT, capital)
#   and their shape (linear or splines)
#   outcome       growth of GDP per capita (log points), the log of GDP per
#                 capita itself (lgdp_pc, a level); with PWT also the
#                 growth of the capital stock per capita (dlk_pc: rnna / pop)
#                 and of GDP per worker (dly_worker: rgdpna / emp), in log
#                 points; NA where PWT has no capital or employment
#   capital       with PWT, the log capital stock per capita as a moderator:
#                 log of its mean over the horizon (kap_level, as long-run
#                 income) or of its mean over the TRAILING_INCOME_YEARS years
#                 before (kap_trailing); NA without capital (and with WB)
#   convergence   last year's log GDP per capita (gdp_l1), the trailing income
#                 (inc_trailing), last year's growth of GDP per capita
#                 (growth_l1) or its growth two years before (growth_l2,
#                 from t-3 to t-2), as a control of its own
#   variability   for the noise of the projections: each country's annual T
#                 and P less its quadratic trend over the horizon
#
#   controls      conflict (minor, war), post-Soviet and post-Yugoslav collapse
#                 years, each group optional (CONTROL_SETS)
#   outliers      the growth outliers whose cause does not involve weather, in
#                 the 1% or 5% (the notebook's) tails of the GDP source's own
#                 ERA5 panel (PWT or WB), of the change of growth from the
#                 year before (outlier_r1, outlier_r5; the notebook's rule) or
#                 of growth itself (outlier_1, outlier_5), flagged so that a
#                 model can drop any set or none
#
# Each panel's rows are those with the long-run agricultural share, the 30-year
# moments and every deviation, the growth outliers included; the trailing
# models use those of them with the trailing values (all the years of the
# window). The units' values and the support of the deviations are those of
# the rows without the 5% outliers (the notebook's sample).
#
# The panels and helpers are read from the notebook's knitr cache, so render
# weather_arms_robustness.Rmd first. From the project root:
#   Rscript apps/response_explorer/build_models.R

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
})
source("load_functions.R")

cache_dir <- file.path("results", "knitr_cache", "weather_arms_robustness")
rdx <- list.files(cache_dir, pattern = "[.]rdx$", full.names = TRUE)
if (!length(rdx)) stop("No knitr cache in ", cache_dir, ": render weather_arms_robustness.Rmd first.")
rdx <- rdx[order(file.info(rdx)$mtime)]
for (f in rdx) lazyLoad(sub("[.]rdx$", "", f), envir = globalenv())
needed <- c("national_panels", "panel_labels", "wdi_growth", "add_income", "MAX_LAGS", "CONTROLS",
            "OUTCOME", "SPLINE_KNOT_PROBS", "prepare_panel", "OUTLIER_CAUSES", "OUTLIERS_DROPPED",
            "OUTLIER_RULE", "outlier_tails")
missing <- needed[!vapply(needed, exists, TRUE)]
if (length(missing)) stop("Not in the knitr cache: ", paste(missing, collapse = ", "))

# The controls the app can include, by group (together, the notebook's).
CONTROL_SETS <- list(conflict = grep("^conflict_", CONTROLS, value = TRUE),
                     soviet = grep("fsu_collapse", CONTROLS, value = TRUE),
                     yugoslav = grep("yug_collapse", CONTROLS, value = TRUE))
stopifnot(setequal(unlist(CONTROL_SETS), CONTROLS))

# The growth outliers the app can drop, the notebook's rules with tails of
# either size and for each GDP source: the country-years in the lowest and the
# highest 1% (or 5%, the notebook's) of the measure of the rule in the source's
# ERA5 panel (PWT x ERA5, WB x ERA5) from the start year whose annotated cause
# does not involve weather, dropped from every panel of that source. The rules
# (the notebook's OUTLIER_RULES): "reversal", the change of growth from the
# year before (keys r1, r5; the default), and "level", growth itself (keys 1,
# 5). The causes are the notebook's (OUTLIER_CAUSES, every file; a cause is
# that of the country-year, whatever the source). The 1% tails are inside the
# 5% ones, so every country-year of them is annotated.
OUTLIER_TAILS <- c("1" = 0.01, "5" = 0.05)
OUTLIER_RULE_KEYS <- c(reversal = "r", level = "")
OUTLIER_PANELS <- c(PWT = "PWT x ERA5", WB = "WB x ERA5")
non_weather <- with(OUTLIER_CAUSES[OUTLIER_CAUSES$weather == "no", ], paste(GID_0, year))
outlier_years <- lapply(OUTLIER_PANELS, function(label) {
  main_all <- national_panels[[match(label, panel_labels)]]$all
  sets <- list()
  for (rule in names(OUTLIER_RULE_KEYS)) for (s in names(OUTLIER_TAILS)) {
    tails <- outlier_tails(main_all, rule, 2 * OUTLIER_TAILS[[s]])$rows
    missing <- setdiff(paste(tails$GID_0, tails$year), paste(OUTLIER_CAUSES$GID_0, OUTLIER_CAUSES$year))
    if (length(missing)) stop(label, ": ", rule, " outliers without an annotated cause: ", paste(missing, collapse = ", "))
    sets[[paste0(OUTLIER_RULE_KEYS[[rule]], s)]] <- intersect(paste(tails$GID_0, tails$year), non_weather)
  }
  sets
})
# The PWT 5% tails of the notebook's rule are the notebook's (the app's
# default).
DEFAULT_OUTLIERS <- paste0(OUTLIER_RULE_KEYS[[OUTLIER_RULE]], "5")
stopifnot(setequal(outlier_years$PWT[[DEFAULT_OUTLIERS]], paste(OUTLIERS_DROPPED$GID_0, OUTLIERS_DROPPED$year)))
for (g in names(outlier_years)) {
  message("Non-weather growth outliers of ", OUTLIER_PANELS[[g]], ": ",
          paste(names(outlier_years[[g]]), lengths(outlier_years[[g]]), sep = " = ", collapse = ", "))
}

EXPLORER_GDP <- c("PWT", "WB")
EXPLORER_CLIMATE <- c("ERA5", "CRU TS", "UDelaware")
# The climate weights of the panels (the notebook's are population density in
# 2000; area is the plain area average, which has no base year).
EXPLORER_WEIGHTS <- list(population = list(weight = "population density", year = 2000L),
                         area = list(weight = "area", year = NULL))
# The units of the deviations: SD of the variability (the notebook's), or
# absolute (the suffix of the unstandardized columns).
EXPLORER_UNITS <- c(SD = "", absolute = "_unstd")
EXPLORER_BASELINES <- c("Level" = "level", "Horizon" = "all", "5-year" = "lag_5", "10-year" = "lag_10",
                        "20-year" = "lag_20", "30-year" = "lag_30",
                        "NRK filter" = "nrk")
# The trailing windows shorter than the notebook's 20 and 30 years: optional
# for a row (the rows are those with the 20 and 30-year deviations).
SHORT_WINDOWS <- c("lag_5", "lag_10")
# The SD of the level is the country's horizon SD.
baseline_sd_window <- function(b) if (b == "level") "all" else b
# The unit slopes of each style of fixed effects, in the regression and in the
# demeaning of the deviations.
EXPLORER_FE_SLOPES <- list(unit = character(0), linear = "t", quadratic = c("t", "t2"))
# The lags of the Nath-Ramey-Klenow filter and the first year of the climate it
# is fitted on (shock_filters.R).
NRK_LAGS <- 3L
NRK_START <- 1950L
TRAILING_INCOME_YEARS <- 5L
TRAILING_SHARE_YEARS <- 5L
INCOME_UNITS <- c(PWT = "2021 PPP US$ (PWT rgdpna / pop)", WB = "2015 US$ (World Bank)")
CAPITAL_UNITS <- "2021 PPP US$ per person (PWT rnna / pop)"

# The panel `econ` x `source` with the climate weights `w` (an element of
# EXPLORER_WEIGHTS): the notebook's with population weights; otherwise built by
# the notebook's prepare_panel() with those weights (the same sample rules and
# growth-outlier exclusions).
explorer_panel <- function(econ, source, w) {
  if (identical(w$weight, CLIMATE_WEIGHT)) {
    return(national_panels[[match(paste(econ, "x", source), panel_labels)]])
  }
  CLIMATE_WEIGHT <<- w$weight
  CLIMATE_WEIGHT_YEAR <<- w$year
  on.exit({
    CLIMATE_WEIGHT <<- EXPLORER_WEIGHTS$population$weight
    CLIMATE_WEIGHT_YEAR <<- EXPLORER_WEIGHTS$population$year
  })
  prepare_panel(econ, source)
}

# The column of the demeaned deviations of variable v, baseline b (a value of
# EXPLORER_BASELINES), units suffix u and fixed effects fe (first differences:
# b ending in "_d"); R/fit.R builds the same.
deviation_column <- function(v, b, u, fe) paste0("dm_", v, "_", b, u, "_", fe)
# The deviations of v in the rows d as a model uses them (R/fit.R reads them the
# same way): the filter's innovations with the NRK baseline (never demeaned);
# otherwise demeaned as the fixed effects fe, or with fe "none" as they are (v
# less the baseline mean, none for the level, over its SD in SD units). With
# `diff`, their first differences (with fe "none": the column fd_<v>_<b><u>).
deviation_values <- function(d, v, b, u, fe, diff = FALSE) {
  if (b == "nrk") return(d[[paste0("nrk_", v, u, if (diff) "_d")]])
  if (fe != "none") return(d[[deviation_column(v, paste0(b, if (diff) "_d"), u, fe)]])
  if (diff) return(d[[paste0("fd_", v, "_", b, u)]])
  x <- if (b == "level") d[[v]] else d[[v]] - d[[paste0("mean_", v, "_", b)]]
  if (u == "") x / d[[paste0("sd_", v, "_", baseline_sd_window(b))]] else x
}

# The nonlinear autoregressive filter of Nath, Ramey and Klenow (Eq. 4), as
# option (c) of shock_filters.R: pooled over the countries `countries` of panel
# `p` (climate weights `w`, an element of EXPLORER_WEIGHTS),
#   v_it = sum_j g_j v_i,t-j + sum_j h_j v_i,t-j vbar_i + mu_i + tau_it
# (j = 1..NRK_LAGS; vbar_i the country's mean over the horizon; country
# effects, no year effects), fitted on the annual climate from NRK_START to the
# last year of the climate data (the lags from the whole series). The shock is
# tau over its SD in the country over the horizon (nrk_<v>), or tau itself in
# absolute units (nrk_<v>_unstd); sd_<v>_nrk is that SD. Returns the shocks by
# unit and year and the filter's coefficients g and h of T and P.
nrk_filter <- function(p, w, countries) {
  cfg <- panel_config(econ_year_min = START_YEAR, climate_source = p$source, climate_weight = w$weight,
                      climate_weight_year = w$year, clim_history_years = 11L)
  raw <- as.data.frame(load_climate_panel(p$level, cfg))
  raw <- raw[raw$gadm_level == p$level & raw$GID_1 %in% countries & raw$year >= NRK_START,
             c("GID_1", "year", "TM", "RR")]
  key <- function(id, year) paste(id, year)
  in_horizon <- raw$year >= p$start_year & raw$year <= p$end_year
  lags <- paste0("l", seq_len(NRK_LAGS))
  shocks <- raw[c("GID_1", "year")]
  coefs <- list()
  for (v in c("TM", "RR")) {
    df <- raw[c("GID_1", "year", v)]
    for (k in seq_len(NRK_LAGS)) df[[lags[k]]] <- raw[[v]][match(key(raw$GID_1, raw$year - k), key(raw$GID_1, raw$year))]
    df$state <- ave(ifelse(in_horizon, raw[[v]], NA), raw$GID_1, FUN = function(x) mean(x, na.rm = TRUE))
    f <- stats::as.formula(paste(v, "~", paste(c(lags, paste0(lags, ":state")), collapse = " + "), "| GID_1"))
    m <- fixest::feols(f, df, notes = FALSE)
    tau <- rep(NA_real_, nrow(df))
    tau[fixest::obs(m)] <- stats::resid(m)
    sd_i <- ave(ifelse(in_horizon, tau, NA), raw$GID_1, FUN = function(x) stats::sd(x, na.rm = TRUE))
    shocks[[paste0("nrk_", v)]] <- tau / sd_i
    shocks[[paste0("nrk_", v, "_unstd")]] <- tau
    # Their first differences (the style "first difference").
    tau_l1 <- tau[match(key(raw$GID_1, raw$year - 1), key(raw$GID_1, raw$year))]
    shocks[[paste0("nrk_", v, "_d")]] <- (tau - tau_l1) / sd_i
    shocks[[paste0("nrk_", v, "_unstd_d")]] <- tau - tau_l1
    shocks[[paste0("sd_", v, "_nrk")]] <- sd_i
    b <- stats::coef(m)
    coefs[[v]] <- list(g = unname(b[lags]), h = unname(b[paste0(lags, ":state")]))
    stopifnot(all(is.finite(unlist(coefs[[v]]))))
  }
  list(shocks = shocks, coefs = coefs)
}

# The demeaned deviations of panel `p` from every baseline, in SD and in
# absolute units, with each style of fixed effects: each unit's mean (and its
# linear, or linear and quadratic, trend) removed, fitted from the start year
# minus MAX_LAGS to the last year of the panel, as prepare_panel() does with the
# quadratic trend; on all the panel's rows (the outliers too) with all of them
# and positive 30-year SDs.
add_explorer_deviations <- function(p) {
  data <- p$data
  window <- data$year >= p$start_year - MAX_LAGS & data$year <= p$end_year
  last_row <- match(paste(data$GID_1, data$year - 1L), paste(data$GID_1, data$year))
  demean <- function(raw, name) {
    ok <- which(window & is.finite(raw))
    sub <- data.frame(GID_1 = data$GID_1[ok], t = data$t[ok], t2 = data$t2[ok], x = raw[ok])
    for (fe in names(EXPLORER_FE_SLOPES)) {
      slopes <- EXPLORER_FE_SLOPES[[fe]]
      unit_fe <- if (length(slopes)) paste0("GID_1[", paste(slopes, collapse = ", "), "]") else "GID_1"
      m <- fixest::feols(stats::as.formula(paste("x ~ 1 |", unit_fe)), sub, notes = FALSE)
      demeaned <- rep(NA_real_, nrow(data))
      demeaned[ok[fixest::obs(m)]] <- stats::resid(m)
      data[[deviation_column(v, name, u, fe)]] <<- demeaned
    }
    deviation_column(v, name, u, fe = names(EXPLORER_FE_SLOPES))
  }
  # cols: the deviations every row must have (those of before the level and
  # the first differences, which keep the rows as they were); extra: the level
  # and the first differences (fd_: as they are, for fe "none").
  cols <- extra <- character(0)
  for (v in c("TM", "RR")) for (b in setdiff(EXPLORER_BASELINES, "nrk")) for (u in EXPLORER_UNITS) {
    raw <- if (b == "level") data[[v]] else data[[v]] - data[[paste0("mean_", v, "_", b)]]
    if (u == "") raw <- raw / data[[paste0("sd_", v, "_", baseline_sd_window(b))]]
    # (A short window can have no variability, e.g. CRU TS fill values: those
    # deviations are NA, and the rows stay those of the longer windows.)
    raw[!is.finite(raw)] <- NA_real_
    if (b %in% c("level", SHORT_WINDOWS)) extra <- c(extra, demean(raw, b)) else cols <- c(cols, demean(raw, b))
    fd <- raw - raw[last_row]
    data[[paste0("fd_", v, "_", b, u)]] <- fd
    extra <- c(extra, demean(fd, paste0(b, "_d")), paste0("fd_", v, "_", b, u))
  }
  # The notebook's deviations are the quadratic-trend ones.
  stopifnot(isTRUE(all.equal(data$dm_TM_lag_20_quadratic[match(paste(p$clean$GID_1, p$clean$year),
                                                                paste(data$GID_1, data$year))],
                             p$clean$dm_TM_lag_20, check.attributes = FALSE)))
  rows <- p$all %>%
    inner_join(select(data, GID_1, year, all_of(c(cols, extra))), by = c("GID_1", "year")) %>%
    filter(if_all(all_of(cols), is.finite), sd_TM_lag_30 > 0, sd_RR_lag_30 > 0)
  # The unit-years without a row from MAX_LAGS years before the start (the lags
  # of the first years, and of the years after a gap in the rows), with the
  # deviations and the climate.
  history <- data %>%
    filter(year >= p$start_year - MAX_LAGS, year <= p$end_year, GID_1 %in% rows$GID_1,
           !paste(GID_1, year) %in% paste(rows$GID_1, rows$year)) %>%
    select(GID_0, GID_1, year, TM, RR, matches("^(mean|sd)_(TM|RR)_(all|lag_5|lag_10|lag_20|lag_30)$"), all_of(c(cols, extra)))
  list(rows = rows, history = history)
}

# Log of the mean GDP per capita of each unit over the `years` years before
# each year (NA unless all of them have GDP), from the whole history of the
# panel's GDP source (the panel itself starts at the start year). The series is
# the panel's: checked against it over the panel's years.
trailing_income <- function(p, years) {
  stored <- resolve_econ_source(p$econ)
  history <- load_econ_panel(stored) %>%
    filter(gadm_level == GADM_LEVEL_BY_SOURCE[[stored]], is.finite(grp_pc_usd)) %>%
    distinct(GID_1, year, grp_pc_usd)
  check <- inner_join(filter(p$data, is.finite(grp_pc_usd)) %>% distinct(GID_1, year, grp_pc_usd),
                      history, by = c("GID_1", "year"))
  stopifnot(nrow(check) > 0, max(abs(check$grp_pc_usd.x / check$grp_pc_usd.y - 1)) < 1e-8)
  history %>%
    group_by(GID_1) %>%
    tidyr::complete(year = seq(min(year), max(year))) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(window_mean = zoo::rollapplyr(grp_pc_usd, years, mean, fill = NA_real_),
           inc_trailing = log(lag_by_year(window_mean, year))) %>%
    ungroup() %>%
    select(GID_1, year, inc_trailing)
}

# With PWT: the growth of capital per capita and of GDP per worker (log
# points), the log of the mean capital per capita over the horizon years with
# capital (kap_level, as long-run income) and over the `years` years before
# each year (kap_trailing, NA unless all of them have it), from the whole
# history of PWT; the growth of GDP per capita of that history is checked
# against the panel's. With WB, all NA (WB has no capital or employment).
add_capital <- function(p, d, years) {
  cols <- c("dlk_pc", "dly_worker", "kap_level", "kap_trailing")
  history <- capital_history(p, years)
  if (is.null(history)) {
    d[cols] <- NA_real_
    return(d)
  }
  level <- history %>%
    filter(year >= p$start_year, year <= p$end_year, is.finite(k_pc_usd)) %>%
    group_by(GID_1) %>%
    summarise(kap_level = log(mean(k_pc_usd)), .groups = "drop")
  d <- d %>%
    left_join(select(history, GID_1, year, dly_check, dlk_pc, dly_worker, kap_trailing), by = c("GID_1", "year")) %>%
    left_join(level, by = "GID_1")
  chk <- is.finite(d$dly_check)
  stopifnot(sum(chk) > 0.9 * nrow(d), max(abs(d$dly_check[chk] - d[[OUTCOME]][chk])) < 1e-8)
  select(d, -dly_check)
}
# The whole PWT history of capital of panel `p` (NULL with WB): the growth of
# GDP per capita (the check), of capital per capita and of GDP per worker, and
# the trailing capital of `years` years.
capital_history <- function(p, years) {
  stored <- resolve_econ_source(p$econ)
  if (stored != "PWT110") return(NULL)
  load_econ_panel(stored) %>%
    filter(gadm_level == GADM_LEVEL_BY_SOURCE[[stored]]) %>%
    distinct(GID_1, year, grp_pc_usd, k_pc_usd, share_emp_pop) %>%
    group_by(GID_1) %>%
    tidyr::complete(year = seq(min(year), max(year))) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(dly_check = log(grp_pc_usd) - lag_by_year(log(grp_pc_usd), year),
           dlk_pc = log(k_pc_usd) - lag_by_year(log(k_pc_usd), year),
           dly_worker = log(grp_pc_usd / share_emp_pop) - lag_by_year(log(grp_pc_usd / share_emp_pop), year),
           kap_trailing = log(lag_by_year(zoo::rollapplyr(k_pc_usd, years, mean, fill = NA_real_), year))) %>%
    ungroup()
}

# Mean agricultural share (agriculture over the three sectors' value added,
# WDI) of each country over the `years` years before each year (NA unless all
# of them have it).
trailing_share <- function(years) {
  wdi_growth %>%
    mutate(share = va_agr_pc / va_total_pc) %>%
    select(GID_0, year, share) %>%
    group_by(GID_0) %>%
    tidyr::complete(year = seq(min(year), max(year))) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(window_mean = zoo::rollapplyr(share, years, mean, fill = NA_real_),
           agr_trailing = lag_by_year(window_mean, year)) %>%
    ungroup() %>%
    select(GID_0, year, agr_trailing)
}

# The rows of panel `p` for the explorer: the deviations, the long-run
# agricultural share (mean of the lagged share over the horizon years with data;
# the countries without any are dropped), long-run income, and the trailing
# income and share (NA where the window is incomplete).
explorer_rows <- function(p) {
  share <- wdi_growth %>%
    filter(year >= p$start_year, year <= p$end_year, is.finite(s_agr)) %>%
    group_by(GID_0) %>%
    summarise(agr_level = mean(s_agr), .groups = "drop")
  deviations <- add_explorer_deviations(p)
  d <- inner_join(deviations$rows, share, by = "GID_0")
  d <- add_income(p, d)
  income <- trailing_income(p, TRAILING_INCOME_YEARS)
  share_trailing <- trailing_share(TRAILING_SHARE_YEARS)
  d <- left_join(d, income, by = c("GID_1", "year"))
  # Last year's log GDP per capita (the trailing income of one year) and its
  # growth last year and two years before (the convergence controls).
  last_gdp <- trailing_income(p, 1L) %>%
    rename(gdp_l1 = inc_trailing) %>%
    group_by(GID_1) %>%
    mutate(growth_l1 = gdp_l1 - lag_by_year(gdp_l1, year),
           growth_l2 = lag_by_year(growth_l1, year)) %>%
    ungroup()
  d <- left_join(d, last_gdp, by = c("GID_1", "year"))
  # This year's log GDP per capita (the outcome "lgdp").
  this_gdp <- p$data %>% filter(is.finite(grp_pc_usd)) %>% distinct(GID_1, year, grp_pc_usd) %>%
    transmute(GID_1, year, lgdp_pc = log(grp_pc_usd))
  d <- left_join(d, this_gdp, by = c("GID_1", "year"))
  d <- left_join(d, share_trailing, by = c("GID_0", "year"))
  d <- add_capital(p, d, TRAILING_INCOME_YEARS)
  # Last year's T and P (the interactions with the lagged level).
  last_year <- p$data %>% distinct(GID_1, year, TM, RR) %>%
    transmute(GID_1, year = year + 1L, TM_l1 = TM, RR_l1 = RR)
  d <- left_join(d, last_year, by = c("GID_1", "year"))
  stopifnot(all(is.finite(d$TM_l1)), all(is.finite(d$RR_l1)))
  # The unit-years without a row, for the lags of the regressors: the climate
  # and the trailing income, share and capital of those years (NA where their
  # window is incomplete; the long-run values are the unit's, in the rows).
  capital <- capital_history(p, TRAILING_INCOME_YEARS)
  history <- deviations$history %>%
    filter(GID_1 %in% d$GID_1) %>%
    left_join(last_year, by = c("GID_1", "year")) %>%
    left_join(income, by = c("GID_1", "year")) %>%
    left_join(share_trailing, by = c("GID_0", "year"))
  history$kap_trailing <- if (is.null(capital)) NA_real_ else {
    capital$kap_trailing[match(paste(history$GID_1, history$year), paste(capital$GID_1, capital$year))]
  }
  # Whether the row is a non-weather growth outlier of the 1% or 5% tails of
  # either rule for the panel's GDP source (outlier_r1, outlier_r5, outlier_1,
  # outlier_5).
  own <- outlier_years[[p$econ]]
  for (k in names(own)) d[[paste0("outlier_", k)]] <- paste(d$GID_0, d$year) %in% own[[k]]
  control_columns <- c("conflict_low", "conflict_high", "fsu_collapse_year", "yug_collapse_year")
  list(rows = d %>%
         select(GID_0, GID_1, year, t, t2, TM, RR, TM_l1, RR_l1, all_of(OUTCOME), dlk_pc, dly_worker, all_of(control_columns),
                starts_with("outlier_"),
                matches("^(mean|sd)_(TM|RR)_(all|lag_5|lag_10|lag_20|lag_30)$"), starts_with("dm_") & ends_with(names(EXPLORER_FE_SLOPES)),
                starts_with("fd_"), inc_level, inc_trailing, lgdp_pc, gdp_l1, growth_l1, growth_l2, agr_level, agr_trailing, kap_level, kap_trailing),
       history = as.data.frame(history))
}

# The natural variability of the climate of panel `p`, for the noise of the
# projections (R/projection.R): each country's annual T and P (C, m a year)
# over the horizon years less its own quadratic trend over them, as a matrix
# of countries x years for each (the countries with every year).
natural_variability <- function(p) {
  years <- seq(p$start_year, p$end_year)
  d <- p$data[p$data$year %in% years, c("GID_0", "year", "TM", "RR")]
  complete <- names(which(tapply(is.finite(d$TM) & is.finite(d$RR), d$GID_0, sum) == length(years)))
  d <- d[d$GID_0 %in% complete, ]
  d$t <- d$year - years[[1]]
  d$t2 <- d$t^2
  out <- list(years = years)
  for (v in c("TM", "RR")) {
    m <- fixest::feols(stats::as.formula(paste(v, "~ 1 | GID_0[t, t2]")), d, notes = FALSE)
    stopifnot(length(fixest::obs(m)) == nrow(d))
    out[[v]] <- tapply(stats::resid(m), list(d$GID_0, d$year), identity)
  }
  out
}

panels <- list()
units <- list()
support <- list()
variability <- list()
filters <- list()
for (weight in names(EXPLORER_WEIGHTS)) for (gdp in EXPLORER_GDP) for (climate in EXPLORER_CLIMATE) {
  label <- paste(gdp, "x", climate)
  p <- explorer_panel(gdp, climate, EXPLORER_WEIGHTS[[weight]])
  built <- explorer_rows(p)
  rows <- built$rows
  stopifnot(!anyDuplicated(rows[, c("GID_1", "year")]))
  # The NRK filter on the panel's countries; every row (and year before the
  # start) has its innovations.
  f <- nrk_filter(p, EXPLORER_WEIGHTS[[weight]], unique(rows$GID_1))
  rows <- left_join(rows, f$shocks, by = c("GID_1", "year"))
  history <- left_join(built$history, f$shocks, by = c("GID_1", "year"))
  stopifnot(all(is.finite(as.matrix(rows[grep("nrk", names(rows))]))))
  filters[[paste(gdp, climate, weight, sep = " | ")]] <- f$coefs
  # The default rows (the notebook's: the 5% outliers of its rule dropped) for
  # the units' values and the support of the deviations.
  base <- rows[!rows[[paste0("outlier_", DEFAULT_OUTLIERS)]], ]
  message(label, ", ", weight, " weights: ", nrow(rows), " rows (non-weather outliers: ",
          paste(names(outlier_years[[gdp]]), vapply(names(outlier_years[[gdp]]), function(k) sum(rows[[paste0("outlier_", k)]]), 0L),
                sep = " = ", collapse = ", "), "), ", n_distinct(rows$GID_0),
          " countries; with trailing income ", sum(is.finite(rows$inc_trailing)),
          ", last year's GDP ", sum(is.finite(rows$gdp_l1)), ", growth ", sum(is.finite(rows$growth_l1)),
          " and growth two years before ", sum(is.finite(rows$growth_l2)),
          ", and trailing share ", sum(is.finite(rows$inc_trailing) & is.finite(rows$agr_trailing)),
          "; capital growth ", sum(is.finite(rows$dlk_pc)), ", growth per worker ", sum(is.finite(rows$dly_worker)),
          ", long-run capital ", sum(is.finite(rows$kap_level)), ", trailing capital ", sum(is.finite(rows$kap_trailing)))
  message("  unit-years without a row, for the lags: ", nrow(history), " (", sum(history$year < p$start_year),
          " before the start), ",
          sum(is.finite(history$dm_TM_lag_20_quadratic)), " with deviations, ",
          sum(is.finite(history$inc_trailing)), " with trailing income")
  panels[[paste(gdp, climate, weight, sep = " | ")]] <- list(
    rows = as.data.frame(rows), history = history, horizon = c(p$start_year, p$end_year))
  variability[[paste(gdp, climate, weight, sep = " | ")]] <- natural_variability(p)
  for (unit in names(EXPLORER_UNITS)) for (b in names(EXPLORER_BASELINES)) {
    window <- EXPLORER_BASELINES[[b]]
    panel_key <- paste(gdp, climate, b, weight, unit, sep = " | ")
    # The units: long-run climate (horizon mean), the SD of the baseline (for a
    # trailing baseline, the median over the unit's rows), income and share.
    units[[panel_key]] <- base %>%
      group_by(GID_0) %>%
      summarise(mean_T = first(mean_TM_all), mean_P = first(mean_RR_all),
                sd_T = stats::median(.data[[paste0("sd_TM_", baseline_sd_window(window))]]),
                sd_P = stats::median(.data[[paste0("sd_RR_", baseline_sd_window(window))]]),
                income = exp(first(inc_level)), share = first(agr_level),
                capital = exp(first(kap_level)), rows = n(),
                .groups = "drop") %>%
      mutate(country = countrycode::countrycode(GID_0, "iso3c", "country.name", warn = FALSE),
             country = ifelse(is.na(country), GID_0, country))
    # The 1st-99th percentile of the deviations in the data, by fixed effects
    # (and "none": the deviations as they are, of the symmetric models), and
    # of their first differences (key ending in " | diff"). The level as it is
    # enters the paths as the anomaly from the country's horizon mean: its
    # support is that of the horizon deviations.
    for (fe in c(names(EXPLORER_FE_SLOPES), "none")) for (diff in c(FALSE, TRUE)) {
      dev <- lapply(c(TM = "TM", RR = "RR"), function(v) {
        b <- if (window == "level" && fe == "none" && !diff) "all" else window
        deviation_values(base, v, b, EXPLORER_UNITS[[unit]], fe, diff)
      })
      support[[paste(panel_key, fe, if (diff) "diff", sep = " | ")]] <- tibble::tibble(
        variable = c("TM", "RR"),
        # (A first difference is NA where the year before has no climate.)
        q01 = c(stats::quantile(dev$TM, 0.01, na.rm = TRUE), stats::quantile(dev$RR, 0.01, na.rm = TRUE)),
        q99 = c(stats::quantile(dev$TM, 0.99, na.rm = TRUE), stats::quantile(dev$RR, 0.99, na.rm = TRUE)),
        share_beyond_3sd = c(mean(abs(dev$TM) > 3, na.rm = TRUE), mean(abs(dev$RR) > 3, na.rm = TRUE))
      )
    }
  }
}

# The most lags the app offers. The history stays from the start year minus
# MAX_LAGS (the deviations are demeaned from then on), so with more than
# MAX_LAGS lags a model's rows start lags - MAX_LAGS years after the start.
APP_MAX_LAGS <- 20L
out_dir <- file.path("results", "response_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(panels = panels, units = units, support = support, variability = variability, filters = filters,
             nrk_lags = NRK_LAGS, nrk_start = NRK_START, max_lags = APP_MAX_LAGS, history_lags = MAX_LAGS, income_units = INCOME_UNITS,
             capital_units = CAPITAL_UNITS, outcome_columns = c(capital = "dlk_pc", worker = "dly_worker", lgdp = "lgdp_pc"),
             baselines = EXPLORER_BASELINES, units_suffix = EXPLORER_UNITS,
             fe_slopes = EXPLORER_FE_SLOPES, outcome = OUTCOME, controls = CONTROLS, control_sets = CONTROL_SETS,
             outlier_tails = OUTLIER_TAILS,
             spline_knot_probs = SPLINE_KNOT_PROBS, trailing_income_years = TRAILING_INCOME_YEARS,
             trailing_share_years = TRAILING_SHARE_YEARS, built = Sys.time(), start_year = START_YEAR),
        file.path(out_dir, "panels.rds"))
message("Saved ", length(panels), " panels to ", file.path(out_dir, "panels.rds"))
