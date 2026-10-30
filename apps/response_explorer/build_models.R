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
#   baseline      horizon, trailing 20 years, trailing 30 years (mean and SD of
#                 the deviation)
#   fixed effects year and unit effects; plus unit linear trends; plus unit
#                 linear and quadratic trends (the notebook's). The deviations
#                 are demeaned in the same way before their arms are cut: each
#                 unit's mean, linear trend, or linear and quadratic trend
#                 removed, fitted from the start year minus MAX_LAGS to the last
#                 year of the panel
#   moderators    fixed: the horizon mean of the climate, long-run income and
#                 long-run share; trailing: the climate of the deviation's own
#                 baseline (the trailing 20 or 30 year mean; the horizon mean
#                 with the horizon baseline) and the mean income of the
#                 TRAILING_INCOME_YEARS years before, with their level terms;
#                 with WB GDP, optionally also the mean share of the
#                 TRAILING_SHARE_YEARS years before
#   arms, interactions (for T and P, any subset of climate, income, share and,
#   with PWT, capital) and their shape (linear or splines)
#   outcome       growth of GDP per capita (log points); with PWT also the
#                 growth of the capital stock per capita (dlk_pc: rnna / pop)
#                 and of GDP per worker (dly_worker: rgdpna / emp), in log
#                 points; NA where PWT has no capital or employment
#   capital       with PWT, the log capital stock per capita as a moderator:
#                 log of its mean over the horizon (kap_level, as long-run
#                 income) or of its mean over the TRAILING_INCOME_YEARS years
#                 before (kap_trailing); NA without capital (and with WB)
#
# Each panel's rows are those with the long-run agricultural share, the 30-year
# moments and every deviation; the trailing models use those of them with the
# trailing values (all the years of the window).
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
            "OUTCOME", "SPLINE_KNOT_PROBS", "prepare_panel")
missing <- needed[!vapply(needed, exists, TRUE)]
if (length(missing)) stop("Not in the knitr cache: ", paste(missing, collapse = ", "))

EXPLORER_GDP <- c("PWT", "WB")
EXPLORER_CLIMATE <- c("ERA5", "CRU TS", "UDelaware")
# The climate weights of the panels (the notebook's are population density in
# 2000; area is the plain area average, which has no base year).
EXPLORER_WEIGHTS <- list(population = list(weight = "population density", year = 2000L),
                         area = list(weight = "area", year = NULL))
# The units of the deviations: SD of the variability (the notebook's), or
# absolute (the suffix of the unstandardized columns).
EXPLORER_UNITS <- c(SD = "", absolute = "_unstd")
EXPLORER_BASELINES <- c("Horizon" = "all", "20-year" = "lag_20", "30-year" = "lag_30")
# The unit slopes of each style of fixed effects, in the regression and in the
# demeaning of the deviations.
EXPLORER_FE_SLOPES <- list(unit = character(0), linear = "t", quadratic = c("t", "t2"))
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
# EXPLORER_BASELINES), units suffix u and fixed effects fe; R/fit.R builds the same.
deviation_column <- function(v, b, u, fe) paste0("dm_", v, "_", b, u, "_", fe)

# The demeaned deviations of panel `p` from every baseline, in SD and in
# absolute units, with each style of fixed effects: each unit's mean (and its
# linear, or linear and quadratic, trend) removed, fitted from the start year
# minus MAX_LAGS to the last year of the panel, as prepare_panel() does with the
# quadratic trend; on the panel's clean rows with all of them and positive
# 30-year SDs.
add_explorer_deviations <- function(p) {
  data <- p$data
  window <- data$year >= p$start_year - MAX_LAGS & data$year <= p$end_year
  cols <- character(0)
  for (v in c("TM", "RR")) for (b in EXPLORER_BASELINES) for (u in EXPLORER_UNITS) {
    raw <- data[[v]] - data[[paste0("mean_", v, "_", b)]]
    if (u == "") raw <- raw / data[[paste0("sd_", v, "_", b)]]
    ok <- which(window & is.finite(raw))
    sub <- data.frame(GID_1 = data$GID_1[ok], t = data$t[ok], t2 = data$t2[ok], x = raw[ok])
    for (fe in names(EXPLORER_FE_SLOPES)) {
      slopes <- EXPLORER_FE_SLOPES[[fe]]
      unit_fe <- if (length(slopes)) paste0("GID_1[", paste(slopes, collapse = ", "), "]") else "GID_1"
      m <- fixest::feols(stats::as.formula(paste("x ~ 1 |", unit_fe)), sub, notes = FALSE)
      demeaned <- rep(NA_real_, nrow(data))
      demeaned[ok[fixest::obs(m)]] <- stats::resid(m)
      cols <- c(cols, deviation_column(v, b, u, fe))
      data[[cols[[length(cols)]]]] <- demeaned
    }
  }
  # The notebook's deviations are the quadratic-trend ones.
  stopifnot(isTRUE(all.equal(data$dm_TM_lag_20_quadratic[match(paste(p$clean$GID_1, p$clean$year),
                                                                paste(data$GID_1, data$year))],
                             p$clean$dm_TM_lag_20, check.attributes = FALSE)))
  p$clean %>%
    inner_join(select(data, GID_1, year, all_of(cols)), by = c("GID_1", "year")) %>%
    filter(if_all(all_of(cols), is.finite), sd_TM_lag_30 > 0, sd_RR_lag_30 > 0)
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
  stored <- resolve_econ_source(p$econ)
  if (stored != "PWT110") {
    d[cols] <- NA_real_
    return(d)
  }
  history <- load_econ_panel(stored) %>%
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
  d <- inner_join(add_explorer_deviations(p), share, by = "GID_0")
  d <- add_income(p, d)
  d <- left_join(d, trailing_income(p, TRAILING_INCOME_YEARS), by = c("GID_1", "year"))
  d <- left_join(d, trailing_share(TRAILING_SHARE_YEARS), by = c("GID_0", "year"))
  d <- add_capital(p, d, TRAILING_INCOME_YEARS)
  control_columns <- c("conflict_low", "conflict_high", "fsu_collapse_year", "yug_collapse_year")
  d %>%
    select(GID_0, GID_1, year, t, t2, all_of(OUTCOME), dlk_pc, dly_worker, all_of(control_columns),
           matches("^(mean|sd)_(TM|RR)_(all|lag_20|lag_30)$"), starts_with("dm_") & ends_with(names(EXPLORER_FE_SLOPES)),
           inc_level, inc_trailing, agr_level, agr_trailing, kap_level, kap_trailing)
}

panels <- list()
units <- list()
support <- list()
for (weight in names(EXPLORER_WEIGHTS)) for (gdp in EXPLORER_GDP) for (climate in EXPLORER_CLIMATE) {
  label <- paste(gdp, "x", climate)
  p <- explorer_panel(gdp, climate, EXPLORER_WEIGHTS[[weight]])
  rows <- explorer_rows(p)
  stopifnot(!anyDuplicated(rows[, c("GID_1", "year")]))
  message(label, ", ", weight, " weights: ", nrow(rows), " rows, ", n_distinct(rows$GID_0),
          " countries; with trailing income ", sum(is.finite(rows$inc_trailing)),
          ", and trailing share ", sum(is.finite(rows$inc_trailing) & is.finite(rows$agr_trailing)),
          "; capital growth ", sum(is.finite(rows$dlk_pc)), ", growth per worker ", sum(is.finite(rows$dly_worker)),
          ", long-run capital ", sum(is.finite(rows$kap_level)), ", trailing capital ", sum(is.finite(rows$kap_trailing)))
  panels[[paste(gdp, climate, weight, sep = " | ")]] <- list(
    rows = as.data.frame(rows), horizon = c(p$start_year, p$end_year))
  for (unit in names(EXPLORER_UNITS)) for (b in names(EXPLORER_BASELINES)) {
    window <- EXPLORER_BASELINES[[b]]
    panel_key <- paste(gdp, climate, b, weight, unit, sep = " | ")
    # The units: long-run climate (horizon mean), the SD of the baseline (for a
    # trailing baseline, the median over the unit's rows), income and share.
    units[[panel_key]] <- rows %>%
      group_by(GID_0) %>%
      summarise(mean_T = first(mean_TM_all), mean_P = first(mean_RR_all),
                sd_T = stats::median(.data[[paste0("sd_TM_", window)]]),
                sd_P = stats::median(.data[[paste0("sd_RR_", window)]]),
                income = exp(first(inc_level)), share = first(agr_level),
                capital = exp(first(kap_level)), rows = n(),
                .groups = "drop") %>%
      mutate(country = countrycode::countrycode(GID_0, "iso3c", "country.name", warn = FALSE),
             country = ifelse(is.na(country), GID_0, country))
    # The 1st-99th percentile of the deviations in the data, by fixed effects.
    for (fe in names(EXPLORER_FE_SLOPES)) {
      dev <- lapply(c(TM = "TM", RR = "RR"), function(v) {
        rows[[deviation_column(v, window, EXPLORER_UNITS[[unit]], fe)]]
      })
      support[[paste(panel_key, fe, sep = " | ")]] <- tibble::tibble(
        variable = c("TM", "RR"),
        q01 = c(stats::quantile(dev$TM, 0.01), stats::quantile(dev$RR, 0.01)),
        q99 = c(stats::quantile(dev$TM, 0.99), stats::quantile(dev$RR, 0.99)),
        share_beyond_3sd = c(mean(abs(dev$TM) > 3), mean(abs(dev$RR) > 3))
      )
    }
  }
}

out_dir <- file.path("results", "response_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(panels = panels, units = units, support = support, income_units = INCOME_UNITS,
             capital_units = CAPITAL_UNITS, outcome_columns = c(capital = "dlk_pc", worker = "dly_worker"),
             baselines = EXPLORER_BASELINES, units_suffix = EXPLORER_UNITS,
             fe_slopes = EXPLORER_FE_SLOPES, outcome = OUTCOME, controls = CONTROLS,
             spline_knot_probs = SPLINE_KNOT_PROBS, trailing_income_years = TRAILING_INCOME_YEARS,
             trailing_share_years = TRAILING_SHARE_YEARS, built = Sys.time(), start_year = START_YEAR),
        file.path(out_dir, "panels.rds"))
message("Saved ", length(panels), " panels to ", file.path(out_dir, "panels.rds"))
