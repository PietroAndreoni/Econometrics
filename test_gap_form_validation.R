# Validation of the deviation forms of test_gap_form.R (signed, linear arms,
# quadratic arms) on one national panel, by two checks.
#
# The deviation term d (run option DEVIATION), for temperature and
# precipitation alike, scaled by the SD over the entire climate series:
#   gap      the Hamilton trend minus the trailing 30-year mean, with the
#            squared Hamilton cycles as symmetric shocks (test_gap_form.R)
#   anomaly  the annual value minus the trailing 30-year mean (the notebook's
#            zTM, zRR); no shock terms, since d already contains the cycle
#   anomaly_all  as anomaly, but temperature against its mean over the entire
#            series (mean_TM_all); precipitation still against its trailing
#            30-year mean
# Forms, each term also interacted with the variable's long-run climate m:
#   signed d; linear arms max(d, 0), max(-d, 0); quadratic arms max(d, 0)^2,
#   max(-d, 0)^2.
#
# (a) Binned response, for temperature and for precipitation. The interaction
#     of the variable's terms with its long-run climate (mean_TM_all,
#     mean_RR_all) is replaced by quarters (N_GROUPS) of that long-run climate.
#     Within each quarter the response is estimated (i) free, with bins of d
#     centred on zero (the bin around zero is the reference; each tail is
#     pooled inward until it holds MIN_UNITS countries), (ii) by each of the
#     three forms with its own coefficients per quarter, and (iii) by each form
#     interacted with the long-run climate as in test_gap_form.R (the
#     "_interacted" models), its curve evaluated at the quarter's mean long-run
#     climate (climate_groups.csv). Shape test: each form (per quarter, and
#     interacted) plus the bins; the bins are jointly zero when the form leaves
#     no systematic shape unexplained. The other variable enters every model in
#     the encompassing form (signed + quadratic arms, each interacted with its
#     long-run climate), so only the evaluated variable's block differs across
#     models. Restriction tests within each quarter: symmetry of the linear arms
#     (the signed restriction) and, in the encompassing model, the quadratic
#     arms given the signed term and the signed term given the quadratic arms.
#
# (b) Out-of-sample extrapolation, for temperature and for precipitation. Each
#     form of the variable, interacted with its long-run climate (the other
#     variable in the encompassing form, as in (a)), is re-estimated without
#     the country-years of most extreme deviation of that variable, and
#     predicts growth in those years from the fixed effects and country trends
#     estimated on the remaining years. Held out: gap, |d| > 1 SD; anomaly and
#     anomaly_all, the HOLDOUT_SHARE of country-years with the largest |d|.
#     Squared prediction errors are compared across forms with a test of equal
#     mean loss, clustered by country (a panel Diebold-Mariano test); the model
#     without the variable's deviation terms is the benchmark.
#     With HOLDOUT = "global_years", whole calendar years are held out instead:
#     the globally hottest and coldest years for temperature, the wettest and
#     driest for precipitation (highest and lowest mean deviation across the
#     sample's countries; or the years in holdout_years_tm / holdout_years_rr).
#     A held-out year has no estimated year effect, so predicted and actual
#     growth are both taken relative to that year's cross-country mean: the
#     test is whether the forms predict which countries grow more or less in
#     that year. Errors are also reported without the country-years coded as
#     in conflict.
#
# Hamilton trend over the entire climate series. Fixed effects year + country +
# country linear trends; errors clustered by country.
#
# Controls (run option CONTROLS, comma-separated; none by default) enter every
# model of (a) and (b), the benchmarks included:
#   conflict  the UCDP/PRIO armed-conflict indicator (add_conflict() in
#             functions/conflict.R)
#   soviet    former-Soviet-republic (and Albania) x year dummies for their
#             collapse years (add_soviet_collapse() in
#             functions/transition_collapse.R)
#   yugoslavia  former-Yugoslav-republic x year dummies for the 1991-1993
#             collapse (add_yugoslav_collapse(), same file)
# With controls, the bins models are also fitted without them on the same
# sample, and the bin estimates are compared (a_bins_control_comparison).
#
# Outliers (run option OUTLIERS; "none" by default):
#   pinatubo  drops the remaining growth outliers of the cool years after the
#             Pinatubo eruption (OUTLIER_YEARS, 1991-1993): country-years
#             whose growth, net of the fixed effects and of the controls (no
#             climate terms), is among the OUTLIER_SHARE largest in absolute
#             value (flag_growth_outliers() in functions/outlier_sensitivity.R).
#             Every model then uses the remaining sample; the dropped
#             country-years are listed (dropped_outliers).
#
# Temperature scale (run option TEMPERATURE_SCALE; "sd" by default):
#   sd       the temperature deviation in SDs of the country's series, as above
#   celsius  the absolute temperature deviation, in degrees C (no SD scaling;
#            the temperature bins are BIN_WIDTH / 2 degrees wide and the
#            holdout threshold is in degrees). Precipitation stays standardized.
#
# Construction of the climate variables, for both temperature and precipitation:
#   WINDOW        years of the trailing mean the deviation is measured from
#                 (5, 10, 20 or 30; default 30)
#   SD_BASELINE   the SD that scales the deviation: "all", over the entire
#                 series (default), or "trailing", over the same trailing
#                 WINDOW years ("lag_30" is the same with a 30-year window)
#   INTERACTION   the long-run climate the forms are interacted with (and whose
#                 quarters define the groups of (a)): "all", the mean over the
#                 entire series (default), or "lag_30", the trailing 30-year
#                 mean, which then also enters every model as a level term
# START_YEAR drops the economic years before it (e.g. 1970, the first year with
# an ERA5 trailing 30-year mean, to give every climate source the same years).
# PLOT_FORMS (comma-separated) limits the forms drawn in the figures; every form
# is still estimated and tested.
#
# Run arguments, named (defaults shown):
#   Rscript test_gap_form_validation.R econ=PWT source="CRU TS" weight_year=2000 \
#     deviation=gap controls= outliers=none temperature_scale=sd start_year= \
#     window=30 sd=all interaction=all plot_forms=signed,linear_arms,quadratic_arms
# The first seven may also be given by position, in that order.
# Output: results/gap_form_validation/
#   <econ>_<source>_pop<year>[_<deviation>][_<controls>][_dropped_<outliers>]
#   [_celsius][_from<start_year>][_win<window>][_sd<window>][_int30]/
#   (no deviation suffix for gap; no window suffix for 30 years)

source("load_functions.R")

# Run options: named (key=value) or, for the first seven, positional.
OPTION_NAMES <- c("econ", "source", "weight_year", "deviation", "controls", "outliers",
                  "temperature_scale", "start_year", "sd", "interaction", "window",
                  "plot_forms", "holdout", "holdout_years_tm", "holdout_years_rr")
args <- commandArgs(trailingOnly = TRUE)
named <- grepl("=", args, fixed = TRUE)
options_given <- c(
  stats::setNames(as.list(args[!named]), OPTION_NAMES[seq_len(sum(!named))]),
  stats::setNames(as.list(sub("^[^=]*=", "", args[named])), sub("=.*$", "", args[named]))
)
stopifnot(all(names(options_given) %in% OPTION_NAMES))
option <- function(name, default) {
  x <- options_given[[name]]
  if (is.null(x)) default else x
}
ECON <- option("econ", "PWT")
SOURCE <- option("source", "CRU TS")
WEIGHT <- "population density"
WEIGHT_YEAR <- as.integer(option("weight_year", 2000L))
DEVIATION <- option("deviation", "gap")
stopifnot(DEVIATION %in% c("gap", "anomaly", "anomaly_all"))
START_YEAR <- suppressWarnings(as.integer(option("start_year", NA)))
WINDOW <- as.integer(option("window", 30L))
SD_BASELINE <- option("sd", "all")
INTERACTION <- option("interaction", "all")
stopifnot(WINDOW %in% c(5L, 10L, 20L, 30L),
          SD_BASELINE %in% c("all", "trailing", "lag_30"), INTERACTION %in% c("all", "lag_30"))
if (SD_BASELINE == "lag_30") {
  stopifnot(WINDOW == 30L)
  SD_BASELINE <- "trailing"
}
# Column suffix of the SD that scales the deviation.
SD_COLUMN <- if (SD_BASELINE == "trailing") paste0("lag_", WINDOW) else "all"
PLOT_FORMS <- trimws(strsplit(option("plot_forms", "signed,linear_arms,quadratic_arms"), ",")[[1]])
HOLDOUT <- option("holdout", "extremes")
stopifnot(HOLDOUT %in% c("extremes", "global_years"))
IS_GLOBAL_HOLDOUT <- HOLDOUT == "global_years"
parse_years <- function(x) if (is.null(x) || !nzchar(x)) NULL else as.integer(strsplit(x, ",")[[1]])
HOLDOUT_YEARS_GIVEN <- list(TM = parse_years(option("holdout_years_tm", "")),
                            RR = parse_years(option("holdout_years_rr", "")))
IS_GAP <- DEVIATION == "gap"
# Each control: the function that adds it to the panel, its column and its
# formula term.
CONTROL_SPECS <- list(
  conflict = list(build = add_conflict, column = "conflict", term = "conflict"),
  soviet = list(build = add_soviet_collapse, column = "fsu_collapse_year",
                term = "i(fsu_collapse_year, ref = 'none')"),
  yugoslavia = list(build = add_yugoslav_collapse, column = "yug_collapse_year",
                    term = "i(yug_collapse_year, ref = 'none')")
)
CONTROLS <- trimws(strsplit(option("controls", ""), ",")[[1]])
CONTROLS <- CONTROLS[nzchar(CONTROLS)]
stopifnot(all(CONTROLS %in% names(CONTROL_SPECS)))
CONTROL_COLUMNS <- vapply(CONTROL_SPECS[CONTROLS], `[[`, character(1), "column")
CONTROL_TERMS <- vapply(CONTROL_SPECS[CONTROLS], `[[`, character(1), "term")
OUTLIERS <- option("outliers", "none")
TEMPERATURE_SCALE <- option("temperature_scale", "sd")
stopifnot(TEMPERATURE_SCALE %in% c("sd", "celsius"))
IS_CELSIUS <- TEMPERATURE_SCALE == "celsius"
# Units of each variable's deviation.
DEV_UNITS <- c(TM = if (IS_CELSIUS) "°C" else "SD", RR = "SD")
stopifnot(OUTLIERS %in% c("none", "pinatubo"))
OUTLIER_YEARS <- 1991:1993
OUTLIER_SHARE <- 0.025
# Baseline the deviation is measured from, per variable.
BASELINE <- c(
  TM = if (DEVIATION == "anomaly_all") "mean_TM_all" else paste0("mean_TM_lag_", WINDOW),
  RR = paste0("mean_RR_lag_", WINDOW)
)

FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"
ARM_BASELINE <- "lag_30"
ARM_SD <- "all"
CONF_LEVEL <- 0.95
# The anomaly spans about -3 to 4 SD against the gap's -2 to 2, hence wider bins.
BIN_WIDTH <- if (IS_GAP) 0.25 else 0.5
# Per variable. In degrees C the temperature bins are half as wide: the
# interannual SD is about 0.5 degrees in warm climates, so the resolution there
# matches the SD bins.
BIN_WIDTHS <- c(TM = if (IS_CELSIUS) BIN_WIDTH / 2 else BIN_WIDTH, RR = BIN_WIDTH)
MIN_UNITS <- 10L
N_GROUPS <- 4L
HOLDOUT_SD <- 1
HOLDOUT_SHARE <- 0.05
# The long-run climate the forms are interacted with; a trailing 30-year mean
# varies within countries, so it also enters every model as a level term.
CLIMATE <- c(TM = paste0("mean_TM_", INTERACTION), RR = paste0("mean_RR_", INTERACTION))
LEVEL_TERMS <- if (INTERACTION == "lag_30") unname(CLIMATE) else character(0)
VARIABLE_LABELS <- c(TM = "Temperature", RR = "Precipitation")
GROUP_NAMES <- list(TM = c("coldest", "cool", "warm", "hottest"),
                    RR = c("driest", "dry", "wet", "wettest"))
stopifnot(all(lengths(GROUP_NAMES) == N_GROUPS))
CLIMATE_UNITS <- c(TM = "°C", RR = "m")
CLIMATE_DIGITS <- c(TM = 1L, RR = 2L)
# Long-run temperatures (degrees C) at which the (b) curves are evaluated; for
# precipitation, the quartiles of long-run precipitation (set from the data).
AT_TM <- c(5, 15, 25)
DEV_NAME <- if (IS_GAP) "gap" else "deviation"
DEV_DEFINITION <- paste0(
  switch(
    DEVIATION,
    gap = sprintf("Hamilton trend minus trailing %d-year mean", WINDOW),
    anomaly = sprintf("annual value minus trailing %d-year mean", WINDOW),
    anomaly_all = sprintf(paste("annual value minus its mean over the entire series",
                                "(temperature) or its trailing %d-year mean (precipitation)"),
                          WINDOW)
  ),
  if (SD_BASELINE == "trailing") sprintf(", in trailing %d-year SDs", WINDOW) else ""
)
DEV_LABEL <- paste0(
  tools::toTitleCase(DEV_NAME),
  if (IS_CELSIUS) " (temperature in °C, precipitation in SD): " else " (SD): ",
  DEV_DEFINITION
)
BIN_TEXT <- if (IS_CELSIUS) {
  sprintf("%g °C for temperature and %g SD for precipitation", BIN_WIDTHS[["TM"]],
          BIN_WIDTHS[["RR"]])
} else {
  sprintf("%g SD", BIN_WIDTH)
}

# Deviation terms of each form, by suffix: the term's value at a deviation of
# d SD (as in test_gap_form.R). Columns are named g<variable><suffix> for
# either deviation.
TERM_VALUE <- list(
  s = function(d) d,
  p = function(d) pmax(d, 0),
  m = function(d) pmax(-d, 0),
  p2 = function(d) pmax(d, 0)^2,
  m2 = function(d) pmax(-d, 0)^2
)
FORMS <- list(
  signed = "s",
  linear_arms = c("p", "m"),
  quadratic_arms = c("p2", "m2")
)
FORM_LABELS <- c(signed = "Signed (BHM style)", linear_arms = "Linear arms",
                 quadratic_arms = "Quadratic arms")
# Signed + quadratic arms: the other variable's block in (a), and the model of
# the nested restriction tests.
ENCOMPASSING <- c("s", "p2", "m2")
SHOCKS <- if (IS_GAP) c("sTM2", "sRR2") else character(0)
# Climate terms in every model besides the deviation forms.
BASE_TERMS <- c(SHOCKS, LEVEL_TERMS)

OUTPUT_DIR <- file.path(
  "results", "gap_form_validation",
  paste0(ECON, "_", gsub("[^A-Za-z0-9]+", "", SOURCE), "_pop", WEIGHT_YEAR,
         if (IS_GAP) "" else paste0("_", DEVIATION),
         if (length(CONTROLS)) paste0("_", paste(CONTROLS, collapse = "_")) else "",
         if (OUTLIERS == "none") "" else paste0("_dropped_", OUTLIERS),
         if (IS_CELSIUS) "_celsius" else "",
         if (is.na(START_YEAR)) "" else paste0("_from", START_YEAR),
         if (IS_GLOBAL_HOLDOUT) {
           paste0("_holdoutyears",
                  if (!is.null(HOLDOUT_YEARS_GIVEN$TM)) paste0("_T", paste(HOLDOUT_YEARS_GIVEN$TM, collapse = "_")) else "",
                  if (!is.null(HOLDOUT_YEARS_GIVEN$RR)) paste0("_P", paste(HOLDOUT_YEARS_GIVEN$RR, collapse = "_")) else "")
         } else {
           ""
         },
         if (WINDOW != 30L) paste0("_win", WINDOW) else "",
         if (SD_BASELINE == "trailing") paste0("_sd", WINDOW) else "",
         if (INTERACTION == "lag_30") "_int30" else "")
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

fit <- function(terms, data) {
  fit_panel_model(terms, data, fixed_effects = FIXED_EFFECTS, cluster = ~GID_0)
}
# Right-hand side from terms, dropping empty blocks; "1" when nothing is left
# (fixed effects only).
join_terms <- function(...) {
  x <- unlist(list(...))
  x <- x[nzchar(x)]
  if (length(x)) paste(x, collapse = " + ") else "1"
}
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
gap_vars <- function(variable, suffixes) paste0("g", variable, suffixes)
interacted_terms <- function(variable, suffixes) {
  x <- gap_vars(variable, suffixes)
  c(x, paste0(x, ":", CLIMATE[[variable]]))
}
# test_gap_form.R's model: both variables' deviations in the form, interacted.
form_terms <- function(suffixes) {
  join_terms(interacted_terms("TM", suffixes), interacted_terms("RR", suffixes), BASE_TERMS,
             CONTROL_TERMS)
}

# Data ----------------------------------------------------------------------------------

message("Building ", ECON, " x ", SOURCE, ", ", WEIGHT, " ", WEIGHT_YEAR,
        "; deviation: ", DEVIATION)
cfg <- panel_config(
  econ_year_min = min(load_econ_panel(resolve_econ_source(ECON))$year),
  climate_source = SOURCE,
  climate_weight = WEIGHT,
  climate_weight_year = WEIGHT_YEAR,
  clim_history_years = 1L
)
d <- suppressMessages(build_dat(ECON, config = cfg, baseline = ARM_BASELINE,
                                sd_baseline = ARM_SD))
for (v in names(CLIMATE)) {
  # In degrees (no scaling) for an absolute temperature deviation.
  scale <- if (v == "TM" && IS_CELSIUS) 1 else d[[paste0("sd_", v, "_", SD_COLUMN)]]
  scale[!is.na(scale) & scale <= 0] <- NA_real_
  baseline <- d[[BASELINE[[v]]]]
  deviation <- if (IS_GAP) {
    (d[[paste0("trend_", v)]] - baseline) / scale
  } else {
    (d[[v]] - baseline) / scale
  }
  for (s in names(TERM_VALUE)) d[[paste0("g", v, s)]] <- TERM_VALUE[[s]](deviation)
  if (IS_GAP) d[[paste0("s", v, "2")]] <- ((d[[v]] - d[[paste0("trend_", v)]]) / scale)^2
}
# The estimation sample: every model uses the same rows (no lags).
for (control in CONTROLS) d <- CONTROL_SPECS[[control]]$build(d)
needed <- c("dlgrp_pc_usd", "gTMs", "gRRs", SHOCKS, CONTROL_COLUMNS, unname(CLIMATE))
d <- d[stats::complete.cases(as.data.frame(d)[, needed]), ]
if (!is.na(START_YEAR)) d <- d[d$year >= START_YEAR, ]

# Remaining growth outliers of the post-Pinatubo years, dropped from every model.
dropped_outliers <- NULL
if (OUTLIERS == "pinatubo") {
  outliers <- flag_growth_outliers(d, CONTROL_TERMS, FIXED_EFFECTS, OUTLIER_YEARS,
                                   OUTLIER_SHARE)
  dropped_outliers <- outliers$table %>%
    left_join(transmute(d, GID_0, year, temperature_deviation = gTMs), by = c("GID_0", "year"))
  d <- d[!outliers$flag, ]
}

# Quarters of each variable's long-run climate over the estimation sample, and
# bins of its deviation within each quarter: centred on zero, tails pooled to
# MIN_UNITS countries. Bin ids are "<quarter>_<index>"; the bin containing zero
# is the reference of its quarter.
group_info <- list()
bin_info <- list()
bin_terms <- character()
for (v in names(CLIMATE)) {
  climate <- d[[CLIMATE[[v]]]]
  group_col <- paste0(v, "_group")
  bin_col <- paste0("g", v, "_bin")
  breaks <- stats::quantile(climate, seq(0, 1, length.out = N_GROUPS + 1), names = FALSE)
  d[[group_col]] <- factor(as.integer(cut(climate, breaks, include.lowest = TRUE)))
  fmt <- function(x) formatC(x, format = "f", digits = CLIMATE_DIGITS[[v]])
  group_info[[v]] <- tibble::tibble(
    variable = v, group = as.character(seq_len(N_GROUPS)),
    lower = breaks[-length(breaks)], upper = breaks[-1],
    mean_climate = as.numeric(tapply(climate, d[[group_col]], mean)),
    label = sprintf("%s, %s quarter\n%s to %s %s", VARIABLE_LABELS[[v]], GROUP_NAMES[[v]],
                    fmt(lower), fmt(upper), CLIMATE_UNITS[[v]])
  )

  dev <- d[[paste0("g", v, "s")]]
  width <- BIN_WIDTHS[[v]]
  k_max <- ceiling((max(abs(dev)) + width / 2) / width)
  edges <- seq(-width / 2 - k_max * width, width / 2 + k_max * width, by = width)
  d[[bin_col]] <- NA_character_
  for (g in levels(d[[group_col]])) {
    rows <- d[[group_col]] == g
    group_breaks <- pool_tails(dev[rows], edges, d$GID_0[rows], MIN_UNITS)
    index <- as.integer(cut(dev[rows], group_breaks, include.lowest = TRUE))
    ref <- as.integer(cut(0, group_breaks))
    d[[bin_col]][rows] <- sprintf("%s_%02d", g, index)
    bin_info[[length(bin_info) + 1]] <- tibble::tibble(
      variable = v, group = g, index = index, deviation = dev[rows], GID_0 = d$GID_0[rows]
    ) %>%
      group_by(variable, group, index) %>%
      summarise(observations = n(), countries = n_distinct(GID_0),
                mean_deviation = mean(deviation), min_deviation = min(deviation),
                max_deviation = max(deviation), .groups = "drop") %>%
      mutate(bin = sprintf("%s_%02d", group, index), reference = index == ref)
  }
  references <- bind_rows(bin_info) %>% filter(variable == v, reference) %>% pull(bin)
  bin_terms[[v]] <- sprintf("i(%s, ref = c(%s))", bin_col,
                            paste0("'", references, "'", collapse = ", "))
}
group_info <- bind_rows(group_info)
bin_info <- bind_rows(bin_info)

# Held-out country-years of (b), per variable: |d| above the threshold, or
# whole calendar years (the highest and lowest mean deviation across countries).
SIDE_NAMES <- list(TM = c(positive = "warm", negative = "cold"),
                   RR = c(positive = "wet", negative = "dry"))
GLOBAL_YEAR_NAMES <- list(TM = c(max = "hottest", min = "coldest"),
                          RR = c(max = "wettest", min = "driest"))
yearly_means <- d %>%
  group_by(year) %>%
  summarise(countries = n(), TM = mean(gTMs), RR = mean(gRRs), .groups = "drop")
holdout_years <- lapply(stats::setNames(names(CLIMATE), names(CLIMATE)), function(v) {
  if (!is.null(HOLDOUT_YEARS_GIVEN[[v]])) return(HOLDOUT_YEARS_GIVEN[[v]])
  c(yearly_means$year[which.max(yearly_means[[v]])],
    yearly_means$year[which.min(yearly_means[[v]])])
})
year_label <- function(v, y) {
  m <- yearly_means[[v]][yearly_means$year == y]
  rank <- if (m == max(yearly_means[[v]])) {
    GLOBAL_YEAR_NAMES[[v]][["max"]]
  } else if (m == min(yearly_means[[v]])) {
    GLOBAL_YEAR_NAMES[[v]][["min"]]
  } else {
    r <- as.integer(rank(-yearly_means[[v]])[yearly_means$year == y])
    suffix <- if (r %% 100 %in% 11:13 || r %% 10 > 3 || r %% 10 == 0) "th" else c("st", "nd", "rd")[r %% 10]
    paste0(r, suffix, " ", GLOBAL_YEAR_NAMES[[v]][["max"]])
  }
  sprintf("%d (%s, mean %+.2f %s)", y, rank, m, DEV_UNITS[[v]])
}
holdout_threshold <- vapply(names(CLIMATE), function(v) {
  if (IS_GLOBAL_HOLDOUT) {
    NA_real_
  } else if (IS_GAP) {
    HOLDOUT_SD
  } else {
    stats::quantile(abs(d[[paste0("g", v, "s")]]), 1 - HOLDOUT_SHARE, names = FALSE)
  }
}, numeric(1))
holdout <- lapply(stats::setNames(names(CLIMATE), names(CLIMATE)), function(v) {
  if (IS_GLOBAL_HOLDOUT) {
    d$year %in% holdout_years[[v]]
  } else {
    abs(d[[paste0("g", v, "s")]]) > holdout_threshold[[v]]
  }
})
holdout_rule <- function(v) {
  if (IS_GLOBAL_HOLDOUT) {
    paste0("the years ", paste(vapply(holdout_years[[v]], year_label, "", v = v),
                               collapse = " and "))
  } else if (IS_GAP) {
    sprintf("|%s %s| > %g %s", tolower(VARIABLE_LABELS[[v]]), DEV_NAME, HOLDOUT_SD,
            DEV_UNITS[[v]])
  } else {
    sprintf("the %g%% largest |%s %s| (> %.2f %s)", 100 * HOLDOUT_SHARE,
            tolower(VARIABLE_LABELS[[v]]), DEV_NAME, holdout_threshold[[v]], DEV_UNITS[[v]])
  }
}

sample_summary <- tibble::tibble(
  econ = ECON, climate_source = SOURCE, weight = paste(WEIGHT, WEIGHT_YEAR),
  deviation = DEVIATION, temperature_scale = TEMPERATURE_SCALE,
  start_year = START_YEAR, window = WINDOW, sd_baseline = SD_BASELINE,
  interaction = INTERACTION,
  controls = paste(CONTROLS, collapse = ","), outliers = OUTLIERS,
  outliers_dropped = if (is.null(dropped_outliers)) 0L else nrow(dropped_outliers),
  observations = nrow(d), countries = n_distinct(d$GID_0),
  first_year = min(d$year), last_year = max(d$year),
  temperature_deviation_min = min(d$gTMs), temperature_deviation_max = max(d$gTMs),
  precipitation_deviation_min = min(d$gRRs), precipitation_deviation_max = max(d$gRRs),
  holdout_design = HOLDOUT,
  temperature_holdout_years = if (IS_GLOBAL_HOLDOUT) paste(holdout_years$TM, collapse = ",") else NA,
  precipitation_holdout_years = if (IS_GLOBAL_HOLDOUT) paste(holdout_years$RR, collapse = ",") else NA,
  temperature_holdout_threshold = holdout_threshold[["TM"]],
  temperature_holdout = sum(holdout$TM),
  temperature_holdout_countries = n_distinct(d$GID_0[holdout$TM]),
  temperature_holdout_warm = sum(holdout$TM & d$gTMs > 0),
  precipitation_holdout_threshold = holdout_threshold[["RR"]],
  precipitation_holdout = sum(holdout$RR),
  precipitation_holdout_countries = n_distinct(d$GID_0[holdout$RR]),
  precipitation_holdout_wet = sum(holdout$RR & d$gRRs > 0)
)
write_out(sample_summary, "sample_summary")
if (!is.null(dropped_outliers)) write_out(dropped_outliers, "dropped_outliers")
if (IS_GLOBAL_HOLDOUT) {
  write_out(yearly_means %>%
              rename(mean_temperature_deviation = TM, mean_precipitation_deviation = RR) %>%
              mutate(held_out_temperature = year %in% holdout_years$TM,
                     held_out_precipitation = year %in% holdout_years$RR),
            "holdout_years")
}
write_out(mutate(group_info, label = sub("\n", ", ", label, fixed = TRUE)), "climate_groups")
write_out(bin_info, "gap_bins")

# Deviations at which the fitted responses are evaluated, per variable.
GRID <- lapply(stats::setNames(names(CLIMATE), names(CLIMATE)), function(v) {
  dev <- d[[paste0("g", v, "s")]]
  if (IS_GAP) seq(-2, 2, by = 0.05) else seq(floor(2 * min(dev)) / 2, ceiling(2 * max(dev)) / 2,
                                             by = 0.05)
})
AT_CLIMATE <- list(
  TM = AT_TM,
  RR = unname(round(stats::quantile(d[[CLIMATE[["RR"]]]], c(0.25, 0.5, 0.75)), 2))
)

# (a) Binned response by quarter of long-run climate --------------------------------------

group_terms <- function(v, suffixes) {
  paste0("i(", v, "_group, g", v, suffixes, ")", collapse = " + ")
}
group_coef <- function(v, g, s) paste0(v, "_group::", g, ":g", v, s)
other_block <- function(v) interacted_terms(setdiff(names(CLIMATE), v), ENCOMPASSING)

# Response of a per-quarter form at deviations `devs`, relative to the mean
# response over the reference-bin observations, so it is on the bins' scale.
group_curve <- function(model, v, suffixes, g, devs, ref_devs) {
  terms <- names(stats::coef(model))
  G <- t(vapply(devs, function(x) {
    w <- stats::setNames(numeric(length(terms)), terms)
    for (s in suffixes) {
      w[[group_coef(v, g, s)]] <- TERM_VALUE[[s]](x) - mean(TERM_VALUE[[s]](ref_devs))
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(deviation = devs)
}

# The same for a form interacted with the long-run climate, evaluated at the
# long-run climate m.
interacted_curve <- function(model, v, suffixes, m, devs, ref_devs) {
  terms <- names(stats::coef(model))
  G <- t(vapply(devs, function(x) {
    w <- stats::setNames(numeric(length(terms)), terms)
    for (s in suffixes) {
      term <- paste0("g", v, s)
      value <- TERM_VALUE[[s]](x) - mean(TERM_VALUE[[s]](ref_devs))
      w[[term]] <- value
      w[[interaction_name(term, CLIMATE[[v]], terms)]] <- value * m
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(deviation = devs)
}

a_fit <- list()
a_tests <- list()
a_bins <- list()
a_curves <- list()
a_control_bins <- list()
a_control_coefs <- list()
for (v in names(CLIMATE)) {
  message("(a) binned response, ", VARIABLE_LABELS[[v]])
  group_col <- paste0(v, "_group")
  bin_prefix <- paste0("g", v, "_bin::")
  model_terms <- c(
    bins = join_terms(bin_terms[[v]], other_block(v), BASE_TERMS, CONTROL_TERMS),
    vapply(FORMS, function(f) join_terms(group_terms(v, f), other_block(v), BASE_TERMS,
                                         CONTROL_TERMS), character(1)),
    stats::setNames(
      vapply(FORMS, function(f) join_terms(interacted_terms(v, f), other_block(v), BASE_TERMS,
                                           CONTROL_TERMS), character(1)),
      paste0(names(FORMS), "_interacted")
    ),
    encompassing = join_terms(group_terms(v, ENCOMPASSING), other_block(v), BASE_TERMS,
                              CONTROL_TERMS)
  )
  models <- lapply(model_terms, fit, data = d)
  n <- vapply(models, nobs, numeric(1))
  stopifnot(length(unique(n)) == 1)
  compared <- c("bins", names(FORMS), paste0(names(FORMS), "_interacted"))
  a_fit[[v]] <- tibble::tibble(
    variable = v, model = compared, observations = n[compared],
    parameters = vapply(models[compared], function(m) length(stats::coef(m)), integer(1)),
    aic = vapply(models[compared], stats::AIC, numeric(1)),
    bic = vapply(models[compared], stats::BIC, numeric(1)),
    within_r2 = vapply(models[compared], function(m) fixest::fitstat(m, "wr2")[[1]], numeric(1))
  ) %>%
    mutate(aic_vs_signed = aic - aic[model == "signed"],
           bic_vs_signed = bic - bic[model == "signed"])

  # Shape test: the bins added to each form, per quarter and interacted.
  for (f in c(names(FORMS), paste0(names(FORMS), "_interacted"))) {
    shape <- fit(paste(model_terms[[f]], "+", bin_terms[[v]]), d)
    rows <- list(wald_row(shape, paste0("^", bin_prefix), label = "bins = 0, all quarters"))
    for (g in levels(d[[group_col]])) {
      rows[[length(rows) + 1]] <- wald_row(shape, paste0("^", bin_prefix, g, "_"),
                                           label = paste("bins = 0, quarter", g))
    }
    a_tests[[length(a_tests) + 1]] <- bind_rows(rows) %>%
      mutate(variable = v, form = f, test_type = "shape", .before = 1)
  }
  # Form restrictions within each quarter.
  for (g in levels(d[[group_col]])) {
    sym <- lincom(models$linear_arms,
                  stats::setNames(c(1, 1), c(group_coef(v, g, "p"), group_coef(v, g, "m"))),
                  label = paste("positive + negative linear arms = 0 (signed), quarter", g))
    a_tests[[length(a_tests) + 1]] <- bind_rows(
      tibble::tibble(test = sym$test, keep = NA_character_, F_stat = sym$statistic^2,
                     df1 = 1, df2 = NA_real_, p_value = sym$p_value),
      wald_row(models$encompassing, paste0("^", v, "_group::", g, ":g", v, "(p2|m2)$"),
               label = paste("quadratic arms = 0 given signed, quarter", g)),
      wald_row(models$encompassing, paste0("^", v, "_group::", g, ":g", v, "s$"),
               label = paste("signed = 0 given quadratic arms, quarter", g))
    ) %>%
      mutate(variable = v, form = NA_character_, test_type = "restriction", .before = 1)
  }

  bins <- tidy_fixest(models$bins) %>%
    filter(startsWith(term, bin_prefix)) %>%
    mutate(bin = sub(bin_prefix, "", term, fixed = TRUE))
  a_bins[[v]] <- bin_info %>%
    filter(variable == v) %>%
    left_join(select(bins, bin, estimate, std_error, conf_low, conf_high, p_value),
              by = "bin") %>%
    mutate(across(c(estimate, conf_low, conf_high), ~ ifelse(reference, 0, .x)))

  if (length(CONTROLS)) {
    # The same bins without the controls, on the same sample: how much of each
    # bin effect the controls absorb.
    plain <- fit(join_terms(bin_terms[[v]], other_block(v), BASE_TERMS), d)
    stopifnot(nobs(plain) == nobs(models$bins))
    plain_bins <- tidy_fixest(plain) %>%
      filter(startsWith(term, bin_prefix)) %>%
      transmute(bin = sub(bin_prefix, "", term, fixed = TRUE), estimate_without = estimate,
                std_error_without = std_error, conf_low_without = conf_low,
                conf_high_without = conf_high)
    a_control_bins[[v]] <- a_bins[[v]] %>%
      select(variable, group, bin, reference, observations, countries, mean_deviation,
             estimate_with = estimate, std_error_with = std_error,
             conf_low_with = conf_low, conf_high_with = conf_high) %>%
      left_join(plain_bins, by = "bin") %>%
      mutate(across(c(estimate_without, conf_low_without, conf_high_without),
                    ~ ifelse(reference, 0, .x)),
             absorbed = estimate_without - estimate_with,
             share_absorbed = ifelse(reference, NA_real_, absorbed / estimate_without))
    a_control_coefs[[v]] <- bind_rows(lapply(compared, function(m) {
      tidy_fixest(models[[m]]) %>%
        filter(Reduce(`|`, lapply(CONTROL_COLUMNS, function(x) startsWith(term, x)))) %>%
        mutate(variable = v, model = m, .before = 1)
    }))
  }

  dev <- d[[paste0("g", v, "s")]]
  references <- bin_info %>% filter(variable == v, reference) %>% pull(bin)
  for (g in levels(d[[group_col]])) {
    rows <- d[[group_col]] == g
    ref_devs <- dev[rows & d[[paste0("g", v, "_bin")]] %in% references]
    grid <- GRID[[v]][GRID[[v]] >= min(dev[rows]) & GRID[[v]] <= max(dev[rows])]
    m <- group_info$mean_climate[group_info$variable == v & group_info$group == g]
    for (f in names(FORMS)) {
      a_curves[[length(a_curves) + 1]] <- bind_rows(
        group_curve(models[[f]], v, FORMS[[f]], g, grid, ref_devs) %>%
          mutate(estimation = "quarter"),
        interacted_curve(models[[paste0(f, "_interacted")]], v, FORMS[[f]], m, grid, ref_devs) %>%
          mutate(estimation = "interacted")
      ) %>%
        mutate(variable = v, group = g, form = f, climate = m)
    }
  }
}
a_fit <- bind_rows(a_fit)
a_tests <- bind_rows(a_tests)
a_bins <- bind_rows(a_bins)
a_curves <- bind_rows(a_curves)
write_out(a_fit, "a_model_fit")
write_out(a_tests, "a_tests")
write_out(a_bins, "a_bin_estimates")
write_out(a_curves, "a_form_curves")
if (length(CONTROLS)) {
  a_control_bins <- bind_rows(a_control_bins)
  a_control_coefs <- bind_rows(a_control_coefs)
  write_out(a_control_bins, "a_bins_control_comparison")
  write_out(a_control_coefs, "a_control_coefficients")
}

# (b) Out-of-sample extrapolation to the most extreme deviations -------------------------

B_LABELS <- c(none = paste("No", DEV_NAME, "terms"), FORM_LABELS)
B_PAIRS <- list(
  c("quadratic_arms", "linear_arms"),
  c("quadratic_arms", "signed"),
  c("linear_arms", "signed"),
  c("signed", "none"),
  c("linear_arms", "none"),
  c("quadratic_arms", "none")
)

# Effect on growth of a deviation of `devs` (in DEV_UNITS) in variable v at
# long-run climate m.
form_curve <- function(model, v, suffixes, devs, m) {
  terms <- names(stats::coef(model))
  G <- t(vapply(devs, function(x) {
    w <- stats::setNames(numeric(length(terms)), terms)
    for (s in suffixes) {
      term <- paste0("g", v, s)
      w[[term]] <- TERM_VALUE[[s]](x)
      w[[interaction_name(term, CLIMATE[[v]], terms)]] <- TERM_VALUE[[s]](x) * m
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(deviation = devs)
}

b_perf <- list()
b_tests <- list()
b_predictions <- list()
b_curves <- list()
for (v in names(CLIMATE)) {
  message("(b) holdout, ", VARIABLE_LABELS[[v]])
  # The variable's form interacted with its long-run climate; the other
  # variable in the encompassing form, as in (a).
  b_models <- c(
    none = join_terms(other_block(v), BASE_TERMS, CONTROL_TERMS),
    vapply(FORMS, function(f) join_terms(interacted_terms(v, f), other_block(v), BASE_TERMS,
                                         CONTROL_TERMS), character(1))
  )
  train <- d[!holdout[[v]], ]
  test <- d[holdout[[v]], ]
  actual <- test$dlgrp_pc_usd
  if (IS_GLOBAL_HOLDOUT) {
    # A held-out year has no year effect. The year effect is fitted on a copy of
    # `year` (the country trends keep the actual year), the test rows are
    # predicted with any training year's effect, and actual and predicted growth
    # are both taken relative to each held-out year's mean, where that constant
    # cancels.
    train$year_fe <- train$year
    test$year_fe <- train$year[[1]]
    fe_train <- sub("^year \\+", "year_fe +", FIXED_EFFECTS)
    fits_train <- lapply(b_models, function(t) {
      fit_panel_model(t, train, fixed_effects = fe_train, cluster = ~GID_0)
    })
  } else {
    fits_train <- lapply(b_models, fit, data = train)
  }
  pred <- vapply(fits_train, function(m) as.numeric(stats::predict(m, newdata = test)),
                 numeric(nrow(test)))
  ok <- stats::complete.cases(pred)
  if (IS_GLOBAL_HOLDOUT) {
    within_year <- function(x) x - stats::ave(ifelse(ok, x, NA), test$year,
                                              FUN = function(z) mean(z, na.rm = TRUE))
    actual <- within_year(actual)
    pred <- apply(pred, 2, within_year)
  }
  loss <- (actual - pred)^2
  test_dev <- test[[paste0("g", v, "s")]]
  b_predictions[[v]] <- tibble::tibble(
    variable = v, GID_0 = test$GID_0, year = test$year, deviation = test_dev,
    climate = test[[CLIMATE[[v]]]],
    conflict = if ("conflict" %in% names(test)) test$conflict else NA,
    predicted = ok, actual = actual
  ) %>% bind_cols(tibble::as_tibble(pred))
  if (IS_GLOBAL_HOLDOUT) {
    year_keys <- paste0("year_", holdout_years[[v]])
    sides <- c(list(all = rep(TRUE, nrow(test))),
               stats::setNames(lapply(holdout_years[[v]], function(y) test$year == y), year_keys))
    side_names <- c(all = "all",
                    stats::setNames(as.character(holdout_years[[v]]), year_keys))
    if ("conflict" %in% names(test)) {
      sides$no_conflict <- test$conflict == 0
      side_names[["no_conflict"]] <- "all, excluding conflict"
    }
  } else {
    sides <- list(all = rep(TRUE, nrow(test)), positive = test_dev > 0, negative = test_dev < 0)
    side_names <- c(all = "all", SIDE_NAMES[[v]])
  }
  for (side in names(sides)) {
    on_side <- sides[[side]]
    side_name <- side_names[[side]]
    keep <- ok & on_side
    if (!any(keep)) next
    b_perf[[length(b_perf) + 1]] <- tibble::tibble(
      variable = v, side = side_name, model = names(b_models),
      observations = sum(keep), countries = n_distinct(test$GID_0[keep]),
      unpredicted = sum(!ok & on_side),
      mse = colMeans(loss[keep, , drop = FALSE]),
      rmse = sqrt(mse),
      mean_error = colMeans((actual - pred)[keep, , drop = FALSE])
    )
    for (pair in B_PAIRS) {
      diff_data <- tibble::tibble(
        diff = loss[keep, pair[[1]]] - loss[keep, pair[[2]]],
        GID_0 = test$GID_0[keep]
      )
      dm <- fixest::feols(diff ~ 1, data = diff_data, cluster = ~GID_0)
      ct <- fixest::coeftable(dm)
      critical <- stats::qt(1 - (1 - CONF_LEVEL) / 2, n_distinct(diff_data$GID_0) - 1)
      b_tests[[length(b_tests) + 1]] <- tibble::tibble(
        variable = v, side = side_name, model = pair[[1]], against = pair[[2]],
        observations = sum(keep), countries = n_distinct(diff_data$GID_0),
        mse_difference = ct[1, 1], std_error = ct[1, 2],
        conf_low = ct[1, 1] - critical * ct[1, 2], conf_high = ct[1, 1] + critical * ct[1, 2],
        p_value = ct[1, 4]
      )
    }
  }
  for (f in names(FORMS)) for (m in AT_CLIMATE[[v]]) {
    b_curves[[length(b_curves) + 1]] <- form_curve(fits_train[[f]], v, FORMS[[f]], GRID[[v]], m) %>%
      mutate(variable = v, form = f, climate = m)
  }
}
b_perf <- bind_rows(b_perf)
b_tests <- bind_rows(b_tests)
b_curves <- bind_rows(b_curves)
b_predictions <- bind_rows(b_predictions)
write_out(b_perf, "b_holdout_performance")
write_out(b_tests, "b_holdout_tests")
write_out(b_curves, "b_training_curves")
# Growth and each model's prediction in the held-out country-years (both relative
# to the held-out year's mean with global years).
write_out(b_predictions, "b_holdout_predictions")

# Console summary -----------------------------------------------------------------------

options(width = 220)
print(as.data.frame(sample_summary), row.names = FALSE)
cat("\n== (a) Fit by model (AIC/BIC relative to the signed form) ==\n")
print(as.data.frame(a_fit %>% mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)
cat("\n== (a) Tests ==\n")
print(as.data.frame(a_tests %>% select(variable, form, test, df1, F_stat, p_value) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
cat("\n== (b) Holdout performance ==\n")
print(as.data.frame(b_perf %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
cat("\n== (b) Equal-loss tests (model minus against; negative = model predicts better) ==\n")
print(as.data.frame(b_tests %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

# Plots ---------------------------------------------------------------------------------

# Colours follow every form, so a form keeps its colour when PLOT_FORMS drops others.
colours <- grid_colours(unname(FORM_LABELS))
form_levels <- unname(FORM_LABELS[intersect(names(FORM_LABELS), PLOT_FORMS)])
stopifnot(length(form_levels) > 0)
forms_title <- if (length(form_levels) == length(FORM_LABELS)) {
  "the three forms"
} else {
  paste(tolower(form_levels), collapse = " and ")
}
facet_levels <- group_info$label
facet_of <- function(variable, group) {
  factor(group_info$label[match(paste(variable, group), paste(group_info$variable,
                                                              group_info$group))],
         levels = facet_levels)
}
caption <- paste0(ECON, " x ", SOURCE, ", ", WEIGHT, " ", WEIGHT_YEAR, " weights; ",
                  if (IS_GAP) "squared Hamilton shocks included" else "no shock terms",
                  if (IS_CELSIUS) "; temperature deviation in °C (not SD-scaled)" else "",
                  if (length(CONTROLS)) paste0("; controls: ", paste(CONTROLS, collapse = ", ")) else "",
                  if (is.null(dropped_outliers)) "" else {
                    paste0("; ", nrow(dropped_outliers), " growth outliers of ",
                           min(OUTLIER_YEARS), "-", max(OUTLIER_YEARS), " dropped")
                  },
                  ". Fixed effects year + country + country linear trends; ",
                  "errors clustered by country.")
wrap_text <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")
caption <- wrap_text(caption, 175)
rows_note <- paste0("Top: temperature, by quarter of long-run temperature. Bottom: ",
                    "precipitation, by quarter of long-run precipitation.")

estimation_labels <- c(quarter = "Own coefficients per quarter",
                       interacted = "Interacted with long-run climate, at the quarter's mean")
curves_plot <- a_curves %>%
  filter(form %in% PLOT_FORMS) %>%
  mutate(form = factor(FORM_LABELS[form], levels = form_levels),
         estimation = factor(estimation_labels[estimation], levels = estimation_labels),
         facet = facet_of(variable, group),
         across(c(estimate, conf_low, conf_high), ~ 100 * .x))
bins_plot <- a_bins %>%
  mutate(facet = facet_of(variable, group),
         across(c(estimate, conf_low, conf_high), ~ 100 * .x))
p_a <- ggplot() +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_line(data = curves_plot,
            aes(deviation, estimate, colour = form, linetype = estimation,
                linewidth = estimation)) +
  geom_linerange(data = filter(bins_plot, !reference),
                 aes(mean_deviation, ymin = conf_low, ymax = conf_high),
                 colour = GRID_INK, linewidth = 0.5) +
  geom_point(data = bins_plot, aes(mean_deviation, estimate, shape = reference),
             colour = GRID_INK, fill = "white", size = 2.2) +
  facet_wrap(~facet, nrow = 2, scales = "free") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_linetype_manual(values = stats::setNames(c("solid", "22"), estimation_labels),
                        name = NULL) +
  scale_linewidth_manual(values = stats::setNames(c(0.9, 0.7), estimation_labels),
                         name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                     labels = c(`FALSE` = paste0(tools::toTitleCase(DEV_NAME), " bin (95% CI)"),
                                `TRUE` = "Reference bin"),
                     name = NULL) +
  guides(colour = guide_legend(nrow = 1, order = 1),
         linetype = guide_legend(nrow = 1, order = 2, override.aes = list(colour = GRID_INK)),
         linewidth = guide_legend(nrow = 1, order = 2),
         shape = guide_legend(order = 3)) +
  labs(
    x = DEV_LABEL,
    y = "Effect on growth that year (pp), relative to the reference bin",
    title = paste0(tools::toTitleCase(DEV_NAME), " response by quarter of long-run climate: ",
                   "free bins against ", forms_title),
    subtitle = wrap_text(paste0(rows_note, " Bins of ", BIN_TEXT, ", centred on zero ",
                                "(points at the bin's mean ", DEV_NAME, "; tails pooled to ",
                                MIN_UNITS, " countries). Forms with their own coefficients ",
                                "per quarter (solid) and interacted with the long-run climate, ",
                                "evaluated at the quarter's mean (dashed); all curves relative ",
                                "to the reference bin."), 175),
    caption = caption
  ) +
  theme_results()
ggsave(file.path(OUTPUT_DIR, "a_binned_response.png"), p_a, width = 15, height = 8.5, dpi = 300)

if (length(CONTROLS)) {
  control_label <- paste(CONTROLS, collapse = " + ")
  variant_levels <- c(paste("Without", control_label), paste("With", control_label))
  compare_plot <- bind_rows(
    a_control_bins %>% transmute(variable, group, mean_deviation, reference,
                                 variant = variant_levels[[1]], estimate = estimate_without,
                                 conf_low = conf_low_without, conf_high = conf_high_without),
    a_control_bins %>% transmute(variable, group, mean_deviation, reference,
                                 variant = variant_levels[[2]], estimate = estimate_with,
                                 conf_low = conf_low_with, conf_high = conf_high_with)
  ) %>%
    filter(!reference) %>%
    mutate(variant = factor(variant, levels = variant_levels),
           facet = facet_of(variable, group),
           across(c(estimate, conf_low, conf_high), ~ 100 * .x))
  dodge <- position_dodge(width = 0.2)
  p_compare <- ggplot(compare_plot, aes(mean_deviation, estimate, colour = variant,
                                        shape = variant)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_linerange(aes(ymin = conf_low, ymax = conf_high), position = dodge, linewidth = 0.5) +
    geom_point(position = dodge, size = 2.2) +
    facet_wrap(~facet, nrow = 2, scales = "free") +
    scale_colour_manual(values = grid_colours(variant_levels), name = NULL) +
    scale_shape_manual(values = grid_shapes(variant_levels), name = NULL) +
    labs(
      x = DEV_LABEL,
      y = "Effect on growth that year (pp), relative to the reference bin",
      title = paste0(tools::toTitleCase(DEV_NAME), " bins with and without ", control_label),
      subtitle = wrap_text(paste0(rows_note, " Same sample; bins of ", BIN_TEXT,
                                  ", centred on zero, points at the bin's mean ", DEV_NAME,
                                  "; 95% CI clustered by country."), 175),
      caption = caption
    ) +
    theme_results()
  ggsave(file.path(OUTPUT_DIR, "a_bins_control_comparison.png"), p_compare,
         width = 15, height = 8.5, dpi = 300)
}

side_label <- function(v, side) {
  if (IS_GLOBAL_HOLDOUT) {
    if (side == "all") {
      sprintf("%s: all held-out years (%s; n = %d)", VARIABLE_LABELS[[v]],
              paste(holdout_years[[v]], collapse = ", "), sum(holdout[[v]]))
    } else if (side == "all, excluding conflict") {
      sprintf("%s: all held-out years, excluding conflict", VARIABLE_LABELS[[v]])
    } else {
      sprintf("%s: %s", VARIABLE_LABELS[[v]], year_label(v, as.integer(side)))
    }
  } else if (side == "all") {
    sprintf("%s: all held out (|%s| > %.2f %s, n = %d)", VARIABLE_LABELS[[v]], DEV_NAME,
            holdout_threshold[[v]], DEV_UNITS[[v]], sum(holdout[[v]]))
  } else {
    sprintf("%s: %s side", VARIABLE_LABELS[[v]], side)
  }
}
panel_sides <- function(v) {
  if (IS_GLOBAL_HOLDOUT) {
    c("all", as.character(holdout_years[[v]]),
      if ("conflict" %in% names(d)) "all, excluding conflict")
  } else {
    c("all", SIDE_NAMES[[v]])
  }
}
panel_levels <- unlist(lapply(names(CLIMATE), function(v) {
  vapply(panel_sides(v), side_label, character(1), v = v)
}))
tests_plot <- b_tests %>%
  filter(against == "signed" | (model == "quadratic_arms" & against == "linear_arms")) %>%
  mutate(
    comparison = factor(paste(B_LABELS[model], "vs", B_LABELS[against]),
                        levels = unique(paste(B_LABELS[model], "vs", B_LABELS[against]))),
    panel = factor(mapply(side_label, variable, side), levels = panel_levels),
    across(c(mse_difference, conf_low, conf_high), ~ 1e4 * .x)
  )
p_b <- ggplot(tests_plot, aes(mse_difference, comparison)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high), colour = GRID_PALETTE[[1]],
                 linewidth = 0.6) +
  geom_point(colour = GRID_PALETTE[[1]], size = 2.2) +
  facet_wrap(~panel, nrow = 2, scales = "free_x",
             labeller = label_wrap_gen(if (IS_GLOBAL_HOLDOUT) 30 else 45)) +
  labs(
    x = "Difference in mean squared prediction error (pp²); negative = first model predicts better",
    y = NULL,
    title = if (IS_GLOBAL_HOLDOUT) {
      "Out-of-sample prediction of growth in the globally most extreme years"
    } else {
      paste0("Out-of-sample prediction of growth in years with an extreme ", DEV_NAME)
    },
    subtitle = wrap_text(paste0("For each variable, its forms (interacted with the long-run ",
                                "climate) re-estimated without ", holdout_rule("TM"), " (top) or ",
                                holdout_rule("RR"), " (bottom)",
                                if (IS_GLOBAL_HOLDOUT) {
                                  paste0("; growth predicted relative to each held-out ",
                                         "year's cross-country mean")
                                } else {
                                  ""
                                },
                                "; 95% CI clustered by country."), 140),
    caption = caption
  ) +
  theme_results()
ggsave(file.path(OUTPUT_DIR, "b_holdout_tests.png"), p_b, width = 12, height = 7, dpi = 300)

climate_label <- function(v, m) {
  paste0(VARIABLE_LABELS[[v]], ", long-run ", formatC(m, format = "f",
                                                      digits = if (v == "TM") 0L else 2L),
         " ", CLIMATE_UNITS[[v]])
}
curve_levels <- unlist(lapply(names(CLIMATE), function(v) {
  vapply(AT_CLIMATE[[v]], climate_label, character(1), v = v)
}))
train_plot <- b_curves %>%
  filter(form %in% PLOT_FORMS) %>%
  mutate(form = factor(FORM_LABELS[form], levels = form_levels),
         facet = factor(mapply(climate_label, variable, climate), levels = curve_levels),
         across(c(estimate, conf_low, conf_high), ~ 100 * .x))
# The held-out (extrapolated) region of each variable, shaded in its panels. Held-out
# global years cover the whole deviation range, so nothing is shaded then.
held_out <- if (IS_GLOBAL_HOLDOUT) NULL else bind_rows(lapply(names(CLIMATE), function(v) {
  facets <- vapply(AT_CLIMATE[[v]], climate_label, character(1), v = v)
  tibble::tibble(
    facet = factor(rep(facets, 2), levels = curve_levels),
    xmin = rep(c(-Inf, holdout_threshold[[v]]), each = length(facets)),
    xmax = rep(c(-holdout_threshold[[v]], Inf), each = length(facets))
  )
}))
p_c <- ggplot(train_plot, aes(deviation, estimate, colour = form, fill = form)) +
  (if (!is.null(held_out)) {
    geom_rect(data = held_out, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = GRID_MUTED, alpha = 0.08)
  }) +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~facet, nrow = 2, scales = "free") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  labs(
    x = DEV_LABEL,
    y = "Effect on growth that year (pp)",
    title = if (IS_GLOBAL_HOLDOUT) {
      "Forms estimated without each variable's globally most extreme years"
    } else {
      paste0("Forms estimated without each variable's extreme-", DEV_NAME,
             " years (shaded: held out, extrapolated)")
    },
    subtitle = wrap_text(paste0("Top: temperature terms at long-run temperatures of 5, 15 and ",
                                "25 °C. Bottom: precipitation terms at the quartiles of long-run ",
                                "precipitation. Bands: 95% CI."), 140),
    caption = caption
  ) +
  theme_results()
ggsave(file.path(OUTPUT_DIR, "b_training_curves.png"), p_c, width = 12, height = 8, dpi = 300)

message("Results written to ", normalizePath(OUTPUT_DIR))
