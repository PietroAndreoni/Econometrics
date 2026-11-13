# Weather response explorer: the growth response to temperature and
# precipitation deviations, and the GDP effect of a permanent shift, for up to
# four settings (a country's characteristics and a model each), under models
# fitted on demand (R/fit.R) on the panels of build_models.R.
#
#   Rscript apps/response_explorer/build_models.R   # once, after rendering the notebook
#   shiny::runApp("apps/response_explorer")         # or source() this file
#
# R/fit.R fits the models, R/response.R and R/projection.R hold the computations.

library(shiny)
library(bslib)
library(ggplot2)

# The app's folder: the file's own when it is source()d, otherwise the working
# directory (runApp() sets it to the app's folder) or apps/response_explorer
# under it (the project root).
SOURCED_FILE <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
APP_DIR <- local({
  candidates <- c(if (!is.null(SOURCED_FILE)) dirname(SOURCED_FILE), getwd(),
                  file.path(getwd(), "apps", "response_explorer"))
  found <- candidates[file.exists(file.path(candidates, "R", "response.R"))]
  if (!length(found)) stop("Cannot find the app's folder (apps/response_explorer).")
  normalizePath(found[[1]])
})
source(file.path(APP_DIR, "R", "response.R"), local = environment())
source(file.path(APP_DIR, "R", "projection.R"), local = environment())
source(file.path(APP_DIR, "R", "fit.R"), local = environment())
# The exposure-mediated terms, from the exposure explorer's histograms; the
# bin edges and percentiles of daily P they offer.
source(file.path(APP_DIR, "R", "exposure.R"), local = environment())
# The guide page.
source(file.path(APP_DIR, "R", "guide.R"), local = environment())
# The "By characteristic" page.
source(file.path(APP_DIR, "R", "moderation.R"), local = environment())
EXPOSURE_P_EDGES <- if (exposure_available("RR")) exposure_hist("pre")$lo else c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20)
EXPOSURE_P_PCT <- local({
  p <- if (exposure_available("RR")) exposure_hist("pre")$pct$wet$probs else c(0.9, 0.95, 0.99, 0.999)
  stats::setNames(p, paste0(format(100 * p, trim = TRUE), "th pct"))
})

PANELS_FILE <- normalizePath(file.path(APP_DIR, "..", "..", "results", "response_explorer", "panels.rds"),
                             mustWork = FALSE)
if (!file.exists(PANELS_FILE)) {
  stop("No ", PANELS_FILE, ": run Rscript apps/response_explorer/build_models.R from the project root.")
}
EXPLORER <- readRDS(PANELS_FILE)
# The default low-income thresholds of the GDP sources (R/fit.R).
EXPLORER$low_income <- low_income_thresholds(EXPLORER)
SCENARIOS_FILE <- sub("panels[.]rds$", "scenarios.rds", PANELS_FILE)
if (!file.exists(SCENARIOS_FILE)) {
  stop("No ", SCENARIOS_FILE, ": run Rscript apps/response_explorer/build_scenarios.R from the project root.")
}
SCEN <- readRDS(SCENARIOS_FILE)
# Regional patterns of geoengineering by injection latitude (CIDER,
# build_cider.R); without the file the geoengineering page has only the uniform
# pattern.
CIDER_FILE <- sub("panels[.]rds$", "cider.rds", PANELS_FILE)
CIDER <- if (file.exists(CIDER_FILE)) readRDS(CIDER_FILE)
# Historical growth of world GDP in PWT 11.0, for the growth charts of the
# scenarios and termination shocks: in each year the growth of the sum of
# rgdpna (constant national prices, millions of 2021 US$) over the countries
# with it in that year and the year before (no interpolation), from the
# estimation horizon's first year on; its median, interquartile range and
# range over those years.
HISTORY_FROM <- 1971L
PWT_GDP <- local({
  file <- file.path(APP_DIR, "..", "..", "data", "pwt110.xlsx")
  if (!file.exists(file)) return(NULL)
  d <- readxl::read_excel(file, sheet = "Data")
  d <- data.frame(GID_0 = d$countrycode, year = as.integer(d$year), gdp = suppressWarnings(as.numeric(d$rgdpna)))
  d <- d[is.finite(d$gdp) & d$gdp > 0, ]
  key <- paste(d$GID_0, d$year)
  d$before <- d$gdp[match(paste(d$GID_0, d$year - 1L), key)]
  d[is.finite(d$before), ]
})
# Historical growth (% a year, from HISTORY_FROM) of the summed GDP of the
# countries `region` (iso3; NULL: every country, the world) with GDP in both
# years, with its summary (range, quartiles, median); NULL without PWT or
# without data for the region.
pwt_growth <- function(region = NULL) {
  if (is.null(PWT_GDP)) return(NULL)
  d <- if (is.null(region)) PWT_GDP else PWT_GDP[PWT_GDP$GID_0 %in% region, ]
  if (!nrow(d)) return(NULL)
  g <- tapply(d$gdp, d$year, sum) / tapply(d$before, d$year, sum)
  out <- data.frame(year = as.integer(names(g)), growth = 100 * (as.numeric(g) - 1),
                    countries = as.integer(table(d$year)[names(g)]))
  out <- out[out$year >= HISTORY_FROM, ]
  if (!nrow(out)) return(NULL)
  q <- stats::quantile(out$growth, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  list(series = out, summary = data.frame(min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5]))
}
PWT_WORLD <- pwt_growth()
# The layers of the historical growth `h` (pwt_growth()) on a growth chart: the
# annual series up to its last year (left) and, right of the projections, its
# range (thin), interquartile range (thick) and median (point), labelled.
history_growth_layers <- function(h = PWT_WORLD) {
  if (is.null(h)) return(NULL)
  s <- h$summary
  x <- PROJECTION_TO + 8
  list(
    geom_line(data = h$series, mapping = aes(year, growth), inherit.aes = FALSE,
              colour = INK_MUTED, linewidth = 0.6),
    geom_segment(data = s, mapping = aes(x = x, xend = x, y = min, yend = max), inherit.aes = FALSE,
                 colour = INK_MUTED, linewidth = 0.5),
    geom_segment(data = s, mapping = aes(x = x, xend = x, y = q25, yend = q75), inherit.aes = FALSE,
                 colour = INK_MUTED, linewidth = 3),
    geom_point(data = s, mapping = aes(x = x, y = median), inherit.aes = FALSE,
               colour = SURFACE, fill = INK, shape = 21, size = 2.6, stroke = 0.8),
    geom_text(data = s, mapping = aes(x = x + 2.5, y = median, label = sprintf("%.1f", median)),
              inherit.aes = FALSE, colour = INK_MUTED, size = 3.2, hjust = 0)
  )
}
history_growth_note <- function(h = PWT_WORLD, name = "world") {
  if (is.null(h)) return(paste("; no historical growth of", name, "GDP in PWT 11.0"))
  s <- h$summary
  sprintf(paste0("; grey line: historical growth of %s GDP in PWT 11.0 (%d-%d); right: its median %.1f%%, ",
                 "interquartile range %.1f to %.1f%% (thick) and range %.1f to %.1f%% (thin)"),
          name, min(h$series$year), max(h$series$year), s$median, s$q25, s$q75, s$min, s$max)
}
# The countries offered for the regional aggregates of the global charts: those
# with SSP GDP, named as in the map borders (else by iso3), in name order.
COUNTRY_CHOICES <- local({
  iso <- sort(unique(SCEN$ssp$gdp$iso3))
  name <- SCEN$world$name[match(iso, SCEN$world$iso3)]
  name[is.na(name)] <- iso[is.na(name)]
  o <- order(name)
  stats::setNames(iso[o], paste0(name[o], " (", iso[o], ")"))
})
COUNTRY_NAMES <- stats::setNames(sub(" \\([A-Z]{3}\\)$", "", names(COUNTRY_CHOICES)), COUNTRY_CHOICES)
# The name of region `region` (iso3; NULL: the world) on the charts.
region_name <- function(region) {
  if (is.null(region)) return("world")
  if (length(region) <= 3) return(paste(COUNTRY_NAMES[region], collapse = " + "))
  sprintf("%d selected countries'", length(region))
}
# The panels of the global charts: one y axis for all, or each its own.
COMMON_Y_LABEL <- "Global chart: common y axis across the panels"
y_scales <- function(common) if (isTRUE(common)) "fixed" else "free_y"
REGION_LABEL <- "Global chart: countries (empty: the world; several: their aggregate GDP)"
region_input <- function(id) {
  selectizeInput(id, REGION_LABEL, COUNTRY_CHOICES, multiple = TRUE,
                 options = list(placeholder = "World", plugins = list("remove_button")))
}
# The maps draw the sf borders of SCEN$world.
loadNamespace("sf")
# Each country's area by 5-degree latitude band, for the exposure chart of the
# termination shocks (R/projection.R).
LATITUDE_SHARES <- latitude_shares(SCEN$world, width = 5)
# The reference climate of the exceedances of natural variability: a trailing
# mean (years) or the panel's horizon (NA).
EXCEEDANCE_REFERENCES <- c("Trailing 20 years" = "20", "Trailing 30 years" = "30",
                           "Horizon mean (fixed)" = "horizon")
# Each forcing scenario (code ssp119, ...) gives the global warming; its
# forcing in W/m2 from the code (ssp245: 4.5).
FORCING <- stats::setNames(sub("^ssp[0-9]([0-9])([0-9])$", "\\1.\\2", names(SCEN$scenarios)),
                           names(SCEN$scenarios))
FORCING_LABELS <- stats::setNames(paste0(FORCING, " W/m2"), names(FORCING))
SCENARIO_CHOICES <- stats::setNames(names(SCEN$scenarios),
                                    sprintf("%s (%s)", FORCING_LABELS, SCEN$scenarios))
# The economic development (income, population, agricultural share): each
# forcing's own SSP, or one SSP for all.
ECONOMY_CHOICES <- c("Follow SSP-RCP" = "follow", SSP1 = "SSP1", SSP2 = "SSP2", SSP3 = "SSP3",
                     SSP4 = "SSP4", SSP5 = "SSP5")
# The warming charts: global warming itself, or its rate (velocity), 10 times
# the change from the year before (C per decade).
WARMING_CHOICES <- c("Global warming (C)" = "level", "Rate of warming (C per decade)" = "rate")
RATE_AXIS <- "Rate of global warming (C per decade)"
RATE_TEXT <- "rate of warming, 10 x the change from the year before (C per decade)"
warming_rate <- function(years, g) {
  out <- 10 * (g - g[match(years - 1L, years)])
  out[is.na(out)] <- NA_real_
  out
}
# The option of the world aggregate without gains (R/projection.R no_gains).
NO_GAINS_LABEL <- paste("Ignore gains in the world aggregate: a country whose GDP gains counts at its",
                        "baseline GDP (zero effect; the maps and country results unchanged)")
# What the geoengineering page can compare side by side (each alternative takes
# the place of that option, the others as set).
TS_COMPARE_CHOICES <- c("Nothing" = "none", "Background forcings" = "forcing",
                        "Economic developments" = "economy", "Catch-up times" = "tau",
                        if (!is.null(CIDER)) c("Geoengineering patterns" = "sai"))
# The regional pattern of geoengineering: uniform (each country's CMIP6
# pattern of global warming, perfect geoengineering) or CIDER's by injection
# latitude, with the injection fractions as set ("mix") or all at one latitude.
CIDER_LATS <- c("60N", "30N", "15N", "0N", "15S", "30S", "60S")
CIDER_FRACTION_DEFAULTS <- c(`60N` = 0, `30N` = 0.25, `15N` = 0.25, `0N` = 0, `15S` = 0.25, `30S` = 0.25, `60S` = 0)
# The injection rate of CIDER's single-latitude training runs (Tg SO2/yr).
CIDER_TRAINING_TG <- 12
TS_SAI_CHOICES <- c("Uniform (CMIP6 pattern)" = "uniform", "CIDER, fractions as set" = "mix",
                    stats::setNames(CIDER_LATS, paste("CIDER,", CIDER_LATS, "only")))
# The injection fractions of pattern `sai` ("mix": `fractions`, as set; a
# latitude: all there).
sai_fractions <- function(sai, fractions) {
  if (sai == "mix") return(fractions)
  stats::setNames(as.numeric(CIDER_LATS == sai), CIDER_LATS)
}
sai_label <- function(sai, model, fractions) {
  if (sai == "uniform") return("uniform pattern")
  f <- sai_fractions(sai, fractions)
  f <- f / sum(f)
  where <- if (sai == "mix") {
    paste(sprintf("%s %.0f%%", names(f)[f > 0], 100 * f[f > 0]), collapse = " / ")
  } else paste(sai, "only")
  sprintf("CIDER %s, %s", sub("-.*$", "", CIDER$models[[model]]), where)
}
# The number of noise paths with natural variability (about 5 s a setting and
# forcing per 200 paths).
NOISE_PATH_CHOICES <- c("100", "200", "500", "1000")
# What the global charts of the scenarios show.
MEASURE_CHOICES <- c("Loss against no climate change (% of GDP)" = "loss",
                     "Growth rate of world GDP (% a year)" = "growth",
                     "Change of the growth rate by climate change (pp a year)" = "growth_change")
if (is.null(EXPLORER$variability)) {
  stop("No natural variability in ", PANELS_FILE, ": rerun Rscript apps/response_explorer/build_models.R.")
}
economy_of <- function(scenario, economy) {
  if (economy == "follow") SCEN$ssp$of_scenario[[scenario]] else economy
}
# The label of a forcing with an economy: "SSP2-4.5" when they are paired as in
# the scenario, otherwise "4.5 W/m2 x SSP3".
scenario_label <- function(scenario, economy) {
  ssp <- economy_of(scenario, economy)
  if (identical(ssp, SCEN$ssp$of_scenario[[scenario]])) SCEN$scenarios[[scenario]]
  else sprintf("%s x %s", FORCING_LABELS[[scenario]], ssp)
}
# One blue step per forcing, light to dark (ordinal ramp).
SCENARIO_COLOURS <- stats::setNames(c("#86b6ef", "#5598e7", "#3987e5", "#2a78d6", "#1c5cab",
                                      "#104281", "#0d366b"),
                                    names(SCEN$scenarios))
# Diverging map colours: losses red, gains blue, grey at zero; countries without data.
MAP_LOSS <- "#e34948"
MAP_GAIN <- "#2a78d6"
MAP_ZERO <- "#f0efec"
MAP_NA <- "#d9d8d4"

N_SETTINGS <- 4
# One fixed colour per setting (categorical slots 1-4 in their fixed order).
SETTING_COLOURS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
DEFAULT_COUNTRIES <- c("IND", "ITA", "USA", "BRA")
INK <- "#0b0b0b"
INK_MUTED <- "#52514e"
GRID <- "#e6e5e1"
SURFACE <- "#fcfcfb"
FORM_SHORT <- c("Linear arms" = "linear arms", "Quadratic arms" = "quadratic arms",
                "Linear + quadratic arms" = "linear + quadratic arms")
# The arms of variable v ("TM", "RR") of a model (one form: both).
form_of <- function(m, v) if (length(m$form) == 1) m$form else m$form[[v]]
# The short label of the arms of T and P (one when they are the same).
form_label <- function(form, symmetric = FALSE) {
  form <- c(TM = unname(form[[1]]), RR = unname(form[[length(form)]]))
  symmetric <- c(TM = isTRUE(of_variable(symmetric, "TM")), RR = isTRUE(of_variable(symmetric, "RR")))
  one <- function(v) paste0(FORM_SHORT[[form[[v]]]], if (symmetric[[v]]) " (symmetric)")
  if (form[["TM"]] == form[["RR"]] && symmetric[["TM"]] == symmetric[["RR"]]) one("TM")
  else paste0("T ", one("TM"), ", P ", one("RR"))
}
# The modifiers a variable's deviation terms can be interacted with (in the
# order of build_models.R), their names in the inputs and in the labels.
MODIFIERS <- c("Long-run climate" = "clim", "Last year's value (KW)" = "lag", "This year's value (BHM)" = "cur",
               "Income" = "inc", "Agricultural share" = "agr",
               "Capital per capita (PWT)" = "kap", "Low income (below the threshold)" = "low")
# Ticked by default, by variable: the preferred specification of
# handout_preferred_specification.Rmd (T: climate, scaled by income; P:
# climate and share).
DEFAULT_MODIFIERS <- list(inter_T = "clim", inter_P = c("clim", "agr"))
MODIFIER_SHORT <- c(clim = "climate", lag = "last year's value", cur = "this year's value", inc = "income",
                    agr = "share", kap = "capital", low = "low income")

# What the paths and losses of settings `sp` are of: GDP, capital per capita,
# GDP per worker, or "outcome" when the settings differ (each setting's label
# names its outcome).
quantity_of <- function(sp) {
  outcomes <- unique(vapply(sp, function(s) s$model$outcome %||% "log", ""))
  if (identical(outcomes, "log")) return("GDP")
  if (length(outcomes) > 1) return("outcome")
  c(lgdp = "GDP", capital = "capital per capita", worker = "GDP per worker")[[outcomes]]
}
first_upper <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))
# A regression unit (baseline) in the labels: "level", "20-year", "horizon",
# "NRK filter"; with first differences (style "diff", or "both" with the unit
# itself), e.g. "level, first difference (KW)".
baseline_label <- function(b, style = "simple") {
  out <- switch(b, "NRK filter" = b, tolower(b))
  switch(style, simple = out, diff = paste0(out, ", first difference", if (b == "Level") " (KW)"),
         both = paste0(out, " and its first difference"))
}
model_style <- function(m) if (isTRUE(m$diff)) "diff" else if (isTRUE(m$both)) "both" else "simple"

# Short description of the interactions: one list if T and P have the same.
interaction_label <- function(inter) {
  part <- function(x) if (length(x)) paste(MODIFIER_SHORT[x], collapse = " + ") else "none"
  if (identical(inter$TM, inter$RR)) {
    if (length(inter$TM)) part(inter$TM) else "no interactions"
  } else {
    paste0("T: ", part(inter$TM), ", P: ", part(inter$RR))
  }
}

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))
# Readable names of the model's terms: "T warm arm x long-run T (spline 1)",
# "level: log income".
ARM_LABELS <- c(TM_pos = "T warm arm", TM_neg = "T cold arm", RR_pos = "P wet arm", RR_neg = "P dry arm",
                TM_sym = "T deviation (symmetric)", RR_sym = "P deviation (symmetric)")
CHARACTERISTIC_LABELS <- c(clim_TM = "long-run T", clim_RR = "long-run P", lag_TM = "last year's T",
                           lag_RR = "last year's P", cur_TM = "this year's T", cur_RR = "this year's P",
                           inc = "log income",
                           agr = "agricultural share", kap = "log capital per capita", low = "low income (dummy)")
characteristic_label <- function(z) {
  base <- sub("_nl[0-9]+$", "", z)
  out <- CHARACTERISTIC_LABELS[[base]]
  if (base != z) out <- paste0(out, " (spline ", sub(".*_nl", "", z), ")")
  out
}
# A term's label, with ", lag k" for the terms (and levels) of lag k.
term_label <- function(term) {
  lag <- if (grepl("\\.l[0-9]+$", term)) sub(".*\\.l", "", term)
  term <- gsub("\\.l[0-9]+", "", term)
  label <- if (startsWith(term, "convergence:")) {
    paste("Convergence:", CONVERGENCE_SHORT[[sub("^convergence:", "", term)]])
  } else if (startsWith(term, "quad:")) {
    # "quad:<kind>:quad_TM" or "quad:<kind>:quad_TM_sq".
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    paste0("Control: ", QUAD_SHORT[[parts[[2]]]], " ", if (grepl("TM", parts[[3]])) "T" else "P",
           if (endsWith(parts[[3]], "_sq")) " (squared)")
  } else if (startsWith(term, "level:idol_")) {
    paste0("Level: (Ybar / Y)^(1 - xi_", if (endsWith(term, "TM")) "T" else "P", ")")
  } else if (startsWith(term, "level:")) {
    paste("Level:", characteristic_label(sub("^level:", "", term)))
  } else {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    scaled <- "idol" %in% parts
    # The exposure channel of a term (R/exposure.R).
    expo <- intersect(parts, EXPOSURE_COLUMNS)
    parts <- parts[!parts %in% c("idol", EXPOSURE_COLUMNS)]
    arm <- substring(parts[[1]], 4)
    # The first-difference terms of style "both".
    differenced <- grepl("^(TM|RR)_d_", arm)
    arm <- sub("_d_", "_", arm)
    out <- ARM_LABELS[[sub("2$", "", arm)]]
    if (differenced) out <- paste(out, "of the first difference")
    if (endsWith(arm, "2")) out <- paste(out, "(squared)")
    if (length(expo)) out <- paste(out, "×", EXPOSURE_COLUMN_LABELS[[expo]], "exposure (per 100)")
    if (length(parts) > 1) out <- paste(out, "×", characteristic_label(parts[[2]]))
    if (scaled) out <- paste0(out, " × (Ybar / Y)^(1 - xi_", if (startsWith(arm, "TM")) "T" else "P", ")")
    out
  }
  if (is.null(lag)) label else paste0(label, ", lag ", lag)
}

# The specification of a model in LaTeX (MathJax): the equation, and the
# definitions of its terms (HTML with inline maths).
spec_symbols <- function(m) {
  trailing <- identical(m$moderators, "trailing")
  clim_t <- trailing && trailing_window(m$window)
  list(
    clim = function(v) paste0("c_{", v, if (clim_t) ",it" else ",i", "}"),
    inc = if (trailing) "y_{it}" else "y_{i}",
    agr = if (isTRUE(m$share_trailing)) "s_{it}" else "s_{i}",
    kap = if (trailing) "k_{it}" else "k_{i}",
    low = if (trailing) "L_{it}" else "L_{i}",
    clim_t = clim_t, inc_t = trailing, agr_t = isTRUE(m$share_trailing)
  )
}
spec_equation <- function(m) {
  sym <- spec_symbols(m)
  # The powers of the arms of a variable: linear, squared, or both.
  powers_of <- function(column) {
    switch(form_of(m, column), "Quadratic arms" = "^{2}", "Linear + quadratic arms" = c("", "^{2}"), "")
  }
  # Per variable (v: "T" or "P"): symmetric terms and spline interactions; the
  # level terms have splines when any variable does.
  column_of <- function(v) if (v == "T") "TM" else "RR"
  symmetric_of <- function(v) isTRUE(of_variable(m$symmetric, column_of(v)))
  splines_of <- function(v) identical(of_variable(m$shape, column_of(v)), "splines")
  splines <- any_splines(m$shape)
  # The own-climate modifiers of v (T, P): long-run climate, last year's and
  # this year's value (centred).
  own <- list(clim = sym$clim, lag = function(v) paste0(v, "_{i,t-1}"), cur = function(v) paste0(v, "_{it}"))
  z_of <- function(v, z) if (z %in% names(own)) own[[z]](v) else sym[[z]]
  letter <- c(clim = "c", lag = "\\ell", cur = "v", inc = "y", agr = "s", kap = "k", low = "L")
  # The arm's coefficient and interactions; with a tilde, those of the terms
  # times Ybar / Y (equal dollars); with `sq`, those of the squared arm when
  # both enter (superscript a,2). Symmetric: one coefficient for both sides
  # (superscript 2 for the squared term when both enter).
  # With lags, every coefficient has its lag k (a left superscript).
  lags <- m$lags %||% 0L
  # With style "both", the coefficients of the first-difference terms (`d`)
  # wear a breve.
  inner_of <- function(v, tilde, sq = FALSE, d = FALSE) {
    t <- function(s) {
      if (d) s <- paste0("\\breve{", s, "}")
      paste0(if (lags > 0) "{}^{k}", if (tilde) paste0("\\tilde{", s, "}") else s)
    }
    a <- if (symmetric_of(v)) (if (sq) "2" else "") else if (sq) "a,2" else "a"
    sup <- if (nzchar(a)) paste0("^{", a, "}") else ""
    paste0(t("\\beta"), sup, "_{", v, "}", paste0(vapply(m$inter[[if (v == "T") "TM" else "RR"]], function(z) {
      if (splines_of(v)) paste0(" + ", t("g"), sup, "_{", v, ",", letter[[z]], "}\\!\\left(", z_of(v, z), "\\right)")
      else paste0(" + ", t("\\gamma"), sup, "_{", v, ",", letter[[z]], "}\\, ", z_of(v, z))
    }, ""), collapse = ""))
  }
  forms <- model_income_forms(m)
  ratio_of <- function(column) {
    power <- 1 - forms$xi[[column]]
    base <- paste0("\\frac{\\bar Y}{Y_{", if (identical(m$moderators, "trailing")) "it" else "i", "}}")
    if (power == 1) base else paste0("\\left(", base, "\\right)^{", format(power), "}")
  }
  # d: the block of the first differences (style "both": \Delta x, and the
  # arms D of \Delta x).
  block <- function(v, d = FALSE) {
    column <- if (v == "T") "TM" else "RR"
    ratio <- ratio_of(column)
    powers <- powers_of(column)
    delta <- if (d) "\\Delta " else ""
    one <- function(power) {
      sq <- length(powers) > 1 && power != ""
      inner <- switch(forms$form[[column]],
                      standard = paste0("\\Big(", inner_of(v, FALSE, sq, d), "\\Big)"),
                      dollars = paste0("\\Big(", inner_of(v, TRUE, sq, d), "\\Big) ", ratio),
                      both = paste0("\\Big[\\Big(", inner_of(v, FALSE, sq, d), "\\Big) + \\Big(", inner_of(v, TRUE, sq, d),
                                    "\\Big) ", ratio, "\\Big]"))
      paste0(inner, if (symmetric_of(v)) paste0(" \\, ", delta, "x_{", v, ",it}")
             else paste0(" \\big(D^{a}_{", delta, v, ",it}\\big)"), power)
    }
    terms <- vapply(powers, one, "")
    out <- if (symmetric_of(v)) paste(terms, collapse = " + ") else {
      paste0("\\sum_{a \\in \\{+,-\\}} ",
             if (length(terms) > 1) paste0("\\Big\\{", paste(terms, collapse = " + "), "\\Big\\}") else terms)
    }
    # Exposure: the terms through each channel's exposure slope, with
    # coefficients of their own.
    if (exposure_on(m$exposure, column)) {
      chans <- c(exh = "H", exc = "C", exw = "W", exd = "D")[exposure_columns(m$exposure, column)]
      out <- paste0("\\sum_{e \\in \\{", paste(chans, collapse = ", "), "\\}} \\frac{s^{e}_{i}}{100} \\Big[", out, "\\Big]^{e}")
    }
    over_lags(out)
  }
  # With lags: the sum over k = 0..L of the terms of year t - k.
  at_lag <- function(x) gsub("{i,t-1}", "{i,t-k-1}", gsub("{it}", "{i,t-k}", gsub(",it}", ",i,t-k}", x, fixed = TRUE),
                                                     fixed = TRUE), fixed = TRUE)
  over_lags <- function(x) {
    if (lags == 0 || !nzchar(x)) return(x)
    paste0("\\sum_{k=0}^{", lags, "} \\Big[", at_lag(x), "\\Big]")
  }
  levels <- unique(sub("_nl[0-9]+$", "", gsub("\\.l[0-9]+", "", m$levels$term)))
  # The levels with lags (the time-varying ones), and the constant ones (no
  # unit effects), outside the sum over the lags.
  lagged <- unique(sub("_nl[0-9]+$", "", sub("\\.l[0-9]+$", "", grep("\\.l[0-9]+$", m$levels$term, value = TRUE))))
  level_term <- function(z, k_of) {
    if (startsWith(z, "idol_")) {
      return(paste0(" + ", k_of, "\\phi_{r,", if (z == "idol_TM") "T" else "P", "}\\, ", ratio_of(sub("^idol_", "", z))))
    }
    if (grepl("^(clim|lag|cur)_", z)) {
      kind <- sub("_.*", "", z)
      v <- if (endsWith(z, "TM")) "T" else "P"
      x <- own[[kind]](v)
      name <- paste0(letter[[kind]], "_", v)
    } else {
      x <- switch(z, inc = sym$inc, agr = sym$agr, kap = sym$kap, low = sym$low)
      name <- switch(z, inc = "y", agr = "s", kap = "k", low = "L")
    }
    if (splines) paste0(" + ", k_of, "h_{", name, "}\\!\\left(", x, "\\right)") else paste0(" + ", k_of, "\\phi_{", name, "}\\, ", x)
  }
  with_lags <- if (lags > 0) intersect(levels, lagged) else levels
  level_terms <- vapply(with_lags, level_term, "", if (lags > 0) "{}^{k}" else "")
  constant_terms <- paste(vapply(setdiff(levels, with_lags), level_term, "", ""), collapse = "")
  effects <- m$effects %||% "both"
  unit_fe <- switch(m$fe, unit = "\\alpha_i", linear = "\\alpha_i + \\delta_i t",
                    quadratic = "\\alpha_i + \\delta_i t + \\eta_i t^{2}")
  fe <- switch(effects, both = paste("\\lambda_t +", unit_fe), unit = unit_fe, year = "\\lambda_t")
  convergence <- switch(m$convergence %||% "none", none = "",
                        lag1 = " + \\rho \\ln \\mathrm{GDPpc}_{i,t-1}",
                        growth1 = " + \\rho\\, \\Delta \\ln \\mathrm{GDPpc}_{i,t-1}",
                        growth2 = " + \\rho\\, \\Delta \\ln \\mathrm{GDPpc}_{i,t-2}",
                        mean5 = " + \\rho \\ln\\big(\\tfrac{1}{5}\\textstyle\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\big)")
  lhs <- switch(m$outcome %||% "log",
                log = "\\Delta \\ln \\mathrm{GDPpc}_{it}",
                lgdp = "\\ln \\mathrm{GDPpc}_{it}",
                capital = "\\Delta \\ln \\mathrm{Kpc}_{it}",
                worker = "\\Delta \\ln \\big(\\mathrm{GDP}_{it} / L_{it}\\big)")
  both <- isTRUE(m$both)
  paste0("\\begin{aligned} ", lhs, " = {} & ", block("T"), if (both) paste0(" \\\\ & + ", block("T", TRUE)),
         " \\\\ & + ", block("P"), if (both) paste0(" \\\\ & + ", block("P", TRUE)),
         " \\\\ & ", if (length(level_terms)) paste0(" + ", over_lags(sub("^ \\+ ", "", paste(level_terms, collapse = "")))),
         constant_terms, convergence,
         paste(vapply(names(m$quad %||% character(0))[(m$quad %||% character(0)) != "none"], function(v) {
           s <- if (v == "TM") "T" else "P"
           sprintf(" + \\theta_{%s,1} Z_{%s,it} + \\theta_{%s,2} Z_{%s,it}^{2}", s, s, s, s)
         }, ""), collapse = ""),
         if (length(m$controls %||% CONTROL_CHOICES)) " + \\kappa' C_{it}",
         " + ", fe, " + \\varepsilon_{it}",
         " \\end{aligned}")
}
spec_definitions <- function(m) {
  sym <- spec_symbols(m)
  used <- unique(unlist(m$inter))
  start <- m$horizon[[1]]
  end <- m$horizon[[2]]
  w <- if (!trailing_window(m$window)) NA else sub("lag_", "", m$window)
  # Per variable: the units of the deviations and symmetric terms.
  units_v <- c(TM = of_variable(m$units, "TM"), RR = of_variable(m$units, "RR"))
  absolute <- all(units_v == "absolute")
  nrk <- identical(m$window, "nrk")
  symmetric_v <- c(TM = isTRUE(of_variable(m$symmetric, "TM")), RR = isTRUE(of_variable(m$symmetric, "RR")))
  symmetric <- any(symmetric_v)
  level <- identical(m$window, "level")
  effects <- m$effects %||% "both"
  no_unit <- effects == "year"
  # How a constant characteristic enters on its own.
  constant_text <- if (no_unit) "with its level term (no country effects)" else "absorbed by \\(\\alpha_i\\) on its own"
  baseline <- if (is.na(w)) {
    sprintf(paste0("\\(\\mu_{v,i}\\) and \\(\\sigma_{v,i}\\): the mean and SD of \\(v_{it}\\) over %d-%d ",
                   "(horizon baseline)"), start, end)
  } else {
    sprintf(paste0("\\(\\mu_{v,it} = \\frac{1}{%s}\\sum_{k=1}^{%s} v_{i,t-k}\\) and \\(\\sigma_{v,it}\\) the SD ",
                   "over the same years (trailing %s-year baseline, without year \\(t\\))"), w, w, w)
  }
  mu <- if (is.na(w)) "\\mu_{v,i}" else "\\mu_{v,it}"
  sigma <- if (is.na(w)) "\\sigma_{v,i}" else "\\sigma_{v,it}"
  deviation <- if (level) {
    sprintf("Regression unit of \\(v \\in \\{T, P\\}\\): the level, %s.",
            if (absolute) "\\(x_{v,it} = v_{it}\\) (C, m)"
            else sprintf("\\(x_{v,it} = v_{it} / \\sigma_{v,i}\\), \\(\\sigma_{v,i}\\) the SD of \\(v_{it}\\) over %d-%d", start, end))
  } else if (nrk) {
    f <- m$filter
    sprintf(paste0(
      "Deviation of \\(v \\in \\{T, P\\}\\): the innovation \\(\\tau_{v,it}\\) of the nonlinear autoregressive ",
      "filter of Nath, Ramey and Klenow (Eq. 4), pooled over the countries, ",
      "\\(v_{it} = \\sum_{j=1}^{3} g_j v_{i,t-j} + \\sum_{j=1}^{3} h_j v_{i,t-j} \\bar v_i + \\mu_i + \\tau_{v,it}\\) ",
      "(\\(\\bar v_i\\) the country's mean over %d-%d; country effects, no year effects), fitted on the annual ",
      "climate from 1950; %s. The lags weigh \\(\\phi_j = g_j + h_j \\bar v_i\\): ",
      "T \\(g\\) = (%s), \\(h\\) = (%s); P \\(g\\) = (%s), \\(h\\) = (%s)."),
      start, end, if (absolute) "\\(x_{v,it} = \\tau_{v,it}\\) in absolute units (C, m)"
      else "\\(x_{v,it} = \\tau_{v,it} / \\sigma^{\\tau}_{v,i}\\), with \\(\\sigma^{\\tau}_{v,i}\\) its SD in the country over the horizon",
      paste(sprintf("%.3f", f$TM$g), collapse = ", "), paste(sprintf("%.4f", f$TM$h), collapse = ", "),
      paste(sprintf("%.3f", f$RR$g), collapse = ", "), paste(sprintf("%.3f", f$RR$h), collapse = ", "))
  } else if (absolute) {
    sprintf("Deviation of \\(v \\in \\{T, P\\}\\) in absolute units (C, m): \\(x_{v,it} = v_{it} - %s\\), with %s.",
            mu, baseline)
  } else {
    sprintf("Deviation of \\(v \\in \\{T, P\\}\\) in SD: \\(x_{v,it} = (v_{it} - %s) / %s\\), with %s.",
            mu, sigma, baseline)
  }
  trend <- switch(m$fe, unit = "\\hat a_i", linear = "\\hat a_i + \\hat b_i t",
                  quadratic = "\\hat a_i + \\hat b_i t + \\hat c_i t^{2}")
  outcome <- switch(m$outcome %||% "log",
    log = NULL,
    lgdp = paste0(
      "Outcome: the log of GDP per capita itself (a level, not its growth). The effect of a year is on that year's ",
      "log GDP: in the paths and scenarios the effects are not summed over the years."),
    capital = paste0(
      "Outcome: the growth of the capital stock per capita, \\(\\mathrm{Kpc}_{it}\\) = PWT rnna / pop (constant ",
      "national prices, as GDP). Responses and losses are of capital per capita."),
    worker = paste0(
      "Outcome: the growth of GDP per worker, \\(\\mathrm{GDP}_{it} / L_{it}\\) = PWT rgdpna / emp (persons ",
      "engaged): the growth of GDP per capita less that of employment per head. Responses and losses are of ",
      "GDP per worker."))
  # How the terms of T and P enter: one phrase when they are the same,
  # otherwise one for each, e.g. "; T: entering linearly; P: entering squared".
  per_form <- function(phrases) {
    p <- vapply(c("TM", "RR"), function(v) phrases[[form_of(m, v)]], "")
    if (p[["TM"]] == p[["RR"]]) p[["TM"]] else paste0("; T", sub("^[,:]", ":", p[["TM"]]), "; P", sub("^[,:]", ":", p[["RR"]]))
  }
  powers_text <- per_form(c("Linear arms" = ", entering linearly", "Quadratic arms" = ", entering squared",
                            "Linear + quadratic arms" = ", entering linearly and squared"))
  # How the deviations enter: unsplit (symmetric) or cut into arms.
  arms_text <- function(sym) {
    if (sym) {
      sprintf(paste0("Symmetric: the deviation \\(x_{v,it}\\) enters unsplit and %s (one coefficient for both ",
                     "sides%s)."),
              if (nrk) "as it is" else "not demeaned",
              per_form(c("Quadratic arms" = ": only its square, so years of the same size above and below the baseline have the same effect",
                         "Linear + quadratic arms" = ", linearly and squared",
                         "Linear arms" = ": a year \\(x\\) below the baseline has minus the effect of one \\(x\\) above")))
    } else if (nrk || no_unit) {
      sprintf(paste0("Not demeaned before the arms are cut (%s): arms ",
                     "\\(D^{+}_{v,it} = \\max(x_{v,it}, 0)\\), \\(D^{-}_{v,it} = \\max(-x_{v,it}, 0)\\)%s."),
              if (nrk) "the filter has country effects" else "no country effects", powers_text)
    } else {
      sprintf(paste0("Demeaned before the arms are cut: \\(d_{v,it} = x_{v,it} - (%s)\\), fitted by OLS for each ",
                     "country over %d-%d; arms \\(D^{+}_{v,it} = \\max(d_{v,it}, 0)\\), ",
                     "\\(D^{-}_{v,it} = \\max(-d_{v,it}, 0)\\)%s."),
              trend, start - 10L, end, powers_text)
    }
  }
  out <- c(
    outcome,
    deviation,
    if (units_v[["TM"]] != units_v[["RR"]]) {
      sprintf("Units by variable: temperature deviations %s, precipitation deviations %s (as above for each).",
              if (units_v[["TM"]] == "absolute") "in absolute units (C)" else "in SD",
              if (units_v[["RR"]] == "absolute") "in absolute units (m a year)" else "in SD")
    },
    if (isTRUE(m$diff)) {
      paste0("First order difference: below \\(x_{v,it}\\) stands for \\(x_{v,it} - x_{v,i,t-1}\\) of the unit above",
             if (level) " (with the level, as Kalkuhl and Wenz: this year's value less last year's)", ".")
    } else if (isTRUE(m$both)) {
      paste0("Both: the unit \\(x_{v,it}\\) and its first difference \\(\\Delta x_{v,it} = x_{v,it} - x_{v,i,t-1}\\) ",
             "enter, each with arms and interactions of its own (the coefficients of the first-difference terms ",
             "with a breve, \\(\\breve\\beta\\)); \\(\\Delta x\\) is demeaned and cut into arms \\(D^{a}_{\\Delta v,it}\\) ",
             "as \\(x\\) below. In the ",
             "responses a year at \\(x\\) follows a year at zero (\\(\\Delta x = x\\)); in the paths and scenarios both ",
             "follow the anomalies.",
             if (level && (symmetric || no_unit)) " The level as it is enters as the long-run level plus the anomaly, less the long-run level.")
    } else if (level && (symmetric || no_unit)) {
      paste0("The level as it is: in the responses, paths and scenarios the effect of an anomaly is that of the ",
             "country's long-run level plus the anomaly less that of the long-run level.")
    } else if (level) {
      "Demeaned like the deviations, the level is the same as the deviation from the horizon mean."
    },
    # Symmetric (unsplit) or arms, per variable: one sentence when T and P
    # agree, otherwise one for each.
    if (all(symmetric_v) || !any(symmetric_v)) arms_text(symmetric_v[[1]])
    else paste0(vapply(c("TM", "RR"), function(v) {
      paste0(if (v == "TM") "Temperature: " else "Precipitation: ", arms_text(symmetric_v[[v]]))
    }, ""), collapse = " ")
  )
  if (any(c("lag", "cur") %in% used)) {
    out <- c(out, sprintf(paste0(
      "%s: \\(v \\in \\{T, P\\}\\) (C, m a year) less the median over the countries of their %d-%d means; ",
      "time-varying, so with its level term%s. In the responses this year's value is the country's long-run ",
      "climate plus the anomaly (last year's at the long-run climate); in the paths and scenarios both follow ",
      "the anomalies."),
      paste(c(if ("lag" %in% used) "Last year's value \\(v_{i,t-1}\\) (as Kalkuhl and Wenz)",
              if ("cur" %in% used) "this year's value \\(v_{it}\\) (as Burke, Hsiang and Miguel)"), collapse = " and "),
      start, end, if (all(c("lag", "cur") %in% used)) "s" else ""))
  }
  if ("clim" %in% used) {
    out <- c(out, if (sym$clim_t) {
      sprintf(paste0("Climate: \\(c_{v,it} = \\mu_{v,it} - \\operatorname{med}_j \\bar v_j\\), %s (the deviation's ",
                     "own baseline) less the median over the countries of their %d-%d means."),
              sprintf("the trailing %s-year mean", w), start, end)
    } else {
      sprintf(paste0("Climate: \\(c_{v,i} = \\bar v_i - \\operatorname{med}_j \\bar v_j\\), the country's mean ",
                     "over %d-%d less the median over the countries (constant: %s)."), start, end, constant_text)
    })
  }
  if ("inc" %in% used) {
    out <- c(out, if (sym$inc_t) {
      paste0("Income: \\(y_{it} = \\ln\\big(\\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\big)\\) less ",
             "its mean over the rows (without year \\(t\\)).")
    } else {
      sprintf(paste0("Income: \\(y_{i} = \\ln\\big(\\overline{\\mathrm{GDPpc}}_i\\big)\\), the log of mean GDP ",
                     "per capita over %d-%d, less its mean over the rows."), start, end)
    })
  }
  if ("agr" %in% used) {
    out <- c(out, if (sym$agr_t) {
      paste0("Agricultural share: \\(s_{it} = \\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{share}_{i,t-k}\\) less its ",
             "mean over the rows (without year \\(t\\)).")
    } else {
      sprintf(paste0("Agricultural share: \\(s_{i}\\), the mean of the previous year's share over %d-%d ",
                     "(that is, %d-%d), less its mean over the rows."), start, end, start - 1L, end - 1L)
    })
  }
  if ("kap" %in% used) {
    out <- c(out, if (sym$inc_t) {
      paste0("Capital: \\(k_{it} = \\ln\\big(\\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{Kpc}_{i,t-k}\\big)\\) less ",
             "its mean over the rows (without year \\(t\\)), with \\(\\mathrm{Kpc}\\) = PWT rnna / pop.")
    } else {
      sprintf(paste0("Capital: \\(k_{i} = \\ln\\big(\\overline{\\mathrm{Kpc}}_i\\big)\\), the log of the mean ",
                     "capital stock per capita (PWT rnna / pop) over %d-%d, less its mean over the rows."),
              start, end)
    })
  }
  if ("low" %in% used) {
    y_low <- if (sym$inc_t) "Y_{it}" else "Y_{i}"
    out <- c(out, sprintf(paste0(
      "Low income: \\(%s = 1\\) when \\(%s\\) (the income of the interactions: %s) is below %s %s, otherwise 0. ",
      "The default threshold is the World Bank's low-income line (FY2017, US$1,025 of 2015 GNI per capita), ",
      "converted for PWT to 2021 PPP with the median PWT / WB income ratio of the country-years under it. ",
      "Entered as one more interaction of each arm%s."),
      sym$low, y_low, if (sym$inc_t) "the trailing mean of the 5 years before" else "the mean over the sample",
      format(round(m$low_threshold), big.mark = ","), EXPLORER$income_units[[m$gdp]],
      if (sym$inc_t) ", with its level term" else paste0(" (constant: ", constant_text, ")")))
  }
  forms <- model_income_forms(m)
  scaled_v <- c(TM = "T", RR = "P")[forms$form != "standard"]
  if (length(scaled_v)) {
    y <- if (sym$inc_t) "Y_{it}" else "Y_{i}"
    xi_text <- paste(vapply(names(scaled_v), function(k) {
      sprintf("\\(\\xi_{%s} = %s\\)", scaled_v[[k]], format(forms$xi[[k]]))
    }, ""), collapse = ", ")
    # The country's income of the ratio is that of the interactions (long-run
    # or trailing); Ybar is a constant, the geometric mean of it over the rows.
    income_text <- if (sym$inc_t) {
      paste0("\\(Y_{it} = \\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\), the country's trailing income")
    } else {
      sprintf("\\(Y_{i} = \\overline{\\mathrm{GDPpc}}_i\\), the country's mean GDP per capita over %d-%d", start, end)
    }
    out <- c(out, paste0(
      "Scaled terms of ", paste(scaled_v, collapse = " and "), ": times \\((\\bar Y / ", y, ")^{1-\\xi}\\), with ",
      income_text, " (the income of the interactions), and \\(\\bar Y = \\exp\\big(\\overline{\\ln ", y,
      "}\\big)\\) a constant, the geometric mean of ", if (sym$inc_t) "the trailing income" else "it",
      " over the rows (so \\((\\bar Y / ", y, ")^{1-\\xi} = e^{-(1-\\xi)", sym$inc, "}\\), with \\(", sym$inc,
      "\\) the log income less its mean over the rows); ", xi_text, ". A scaled growth effect is a dollar loss per ",
      "person proportional to \\(", y, "^{\\xi}\\) beyond the ticked interactions (\\(\\xi = 0\\): the same ",
      "dollars at every income, \\(\\tilde\\beta \\bar Y\\) at the median climate and mean share; \\(\\xi = 1\\): ",
      "the same %, the standard terms)",
      if (any(forms$form == "both")) paste0(
        "; with the plain terms (", paste(c(TM = "T", RR = "P")[forms$form == "both"], collapse = " and "),
        "), the model nests equal % (\\(\\tilde\\beta = \\tilde\\gamma = 0\\)) and the scaled form ",
        "(\\(\\beta = \\gamma = 0\\)), tested in the table"), "."))
  }
  if ((m$lags %||% 0L) > 0) {
    out <- c(out, sprintf(paste0(
      "Lags: every deviation term also enters at \\(t-k\\), \\(k = 1, \\dots, %d\\), with coefficients of its own ",
      "(\\({}^{k}\\beta\\), \\({}^{k}\\gamma\\)) and the characteristics of year \\(t-k\\) (the interactions lagged with ",
      "it; the time-varying characteristics with their lagged level terms). The lags of the first years and after ",
      "a gap come from the years without a row (from %d). Responses: the cumulative effect of one year's anomaly on ",
      "that year's growth and the next %d years' (the sum over \\(k\\)); paths and scenarios: each year's growth ",
      "sums the effects of the deviations of that year and the %d before."),
      m$lags, start - 10L, m$lags, m$lags))
  }
  if (any_splines(m$shape)) {
    out <- c(out, paste0("Splines:\\(g(z) = \\gamma z + \\sum_{j=1}^{2} \\gamma_j [r_j(z) - r_j(0)]\\), ",
                         "restricted cubic spline with knots at the 5, 35, 65 and 95th percentiles of the ",
                         "countries' values (linear beyond the outer knots); one set of coefficients per arm."))
  }
  if (nrow(m$levels)) {
    out <- c(out, paste0("Level terms (\\(\\phi\\)", if (any_splines(m$shape)) ", \\(h\\)", "): the ",
                         if (no_unit) "characteristics" else "time-varying characteristics",
                         " that enter an interaction also enter on their own",
                         if (no_unit) " (without country effects the constant ones too, unlagged)", "."))
  }
  if (!identical(m$convergence %||% "none", "none")) {
    growth <- m$convergence %in% c("growth1", "growth2")
    out <- c(out, sprintf(paste0(
      "Convergence: \\(\\rho\\) on %s (%s, the GDP source's), a control whatever the interactions; %s %s ",
      "In the GDP paths and scenarios the climate's effect on GDP feeds back through it: each year adds ",
      "\\(\\rho\\) times the change of the control by the effects of the years before (none before the first ",
      "year; with the 5-year mean, a flat baseline; not for capital), at each draw of \\(\\rho\\) jointly with the ",
      "other coefficients. The response curves are the year's own effect, without it."),
      switch(m$convergence, lag1 = "\\(\\ln \\mathrm{GDPpc}_{i,t-1}\\)",
             mean5 = "\\(\\ln\\big(\\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\big)\\)",
             growth1 = "\\(\\Delta \\ln \\mathrm{GDPpc}_{i,t-1}\\)",
             growth2 = "\\(\\Delta \\ln \\mathrm{GDPpc}_{i,t-2}\\)"),
      CONVERGENCE_SHORT[[m$convergence]],
      if (growth) paste0(if (m$convergence == "growth1") "the lagged dependent variable" else
                           "the dependent variable's second lag, which leaves out year t−1 and its shared shocks",
                         " (for the outcomes other than GDP per capita, still GDP per capita's growth): ",
                         "\\(\\rho\\) is the persistence of growth, and \\(\\rho < 0\\) means part of a shock to ",
                         "growth is reversed ", if (m$convergence == "growth1") "the next year." else "two years later.")
      else "\\(\\rho < 0\\): conditional convergence.",
      if (no_unit) {
        if (growth) "Without country effects, it also picks up differences in mean growth across countries."
        else "Without country effects, it compares income levels across countries."
      } else paste0("With country effects it is a within-country ", if (growth) "coefficient" else "rate",
                    ", with a Nickell bias of order 1/T (T about 50 years).")))
    rho <- m$convergence_coef[["estimate"]] %||% NA
    if (!is.finite(rho)) {
      out <- c(out, paste0("Here the control has no coefficient (collinear with the other terms, e.g. the trailing ",
                           "income's level term): no feedback in the paths and scenarios."))
    } else if (identical(m$outcome %||% "log", "capital")) {
      out <- c(out, "With capital as the outcome there is no feedback (its control is GDP's).")
    } else if (rho > 0 && m$convergence %in% c("lag1", "mean5") && !level_outcome(m)) {
      out <- c(out, sprintf(paste0(
        "Here \\(\\rho = %+.3f > 0\\): divergence, so in the paths and scenarios a loss compounds (by about ",
        "%.0f%% of itself a year) and long horizons explode."), rho, 100 * rho))
    }
    if (isTRUE(m$convergence_iv)) {
      out <- c(out, paste0(
        "Instrumented (2SLS) with its value one year earlier, ", CONVERGENCE_IV_SHORT[[m$convergence]],
        ", against the Nickell bias (as Harding et al. 2025 and Bond et al. 2010)",
        if (is.finite(m$convergence_coef[["first_stage_f"]] %||% NA)) {
          sprintf("; first-stage clustered Wald F = %.1f", m$convergence_coef[["first_stage_f"]])
        }, ". With unit effects the instrument is not strictly exogenous either (its bias is smaller, of the ",
        "order of the shocks it shares with the unit's mean)."))
    }
  }
  for (ev_v in c("TM", "RR")) {
    if (!exposure_on(m$exposure, ev_v)) next
    e <- m$exposure[[ev_v]]
    spec <- EXPOSURE_VARS[[ev_v]]
    chans <- exposure_columns(m$exposure, ev_v)
    sym <- c(exh = "H", exc = "C", exw = "W", exd = "D")
    name <- if (ev_v == "TM") "T" else "P"
    index_text <- function(k) {
      x <- e[[k]]
      above <- k == "hi"
      what <- if (x$measure == "days") "days" else if (ev_v == "TM") "degree days" else "mm"
      thr <- if (x$type == "pct") sprintf("the country's %gth percentile of its daily values (1970-2023, all days)", 100 * x$value)
             else paste(format(x$value), spec$unit)
      sprintf("\\(%s_{it}\\) = %s of the %s %s %s (%s)", sym[[spec$cols[[k]]]], what, EXPOSURE_SOURCE_SHORT[[x$src]],
              if (above) "above" else "below", thr,
              if (x$measure == "days") "count of the days" else if (above) "sum of the daily excess" else "sum of the daily shortfall")
    }
    ks <- names(spec$cols)[spec$cols %in% chans]
    out <- c(out, paste0(
      "Exposure of the ", name, " terms: every ", name, " deviation term times \\(s^{e}_{i}/100\\) for each channel ",
      paste(vapply(ks, function(k) sprintf("\\(e = %s\\) (%s)", sym[[spec$cols[[k]]]], spec$sides[[k]]), ""), collapse = " and "),
      ", with no plain term, so that a deviation acts through the change of the index it brings. Indices (ERA5 ",
      "daily, population weights 2015): ", paste(vapply(ks, index_text, ""), collapse = "; "),
      ". \\(s^{e}_{i}\\): the within-country OLS slope of the index's deviation from its 20-year trailing mean on ",
      "that of the ", spec$annual, ", 1970-2023 (climate data only), per ", if (ev_v == "TM") "C" else "m a year",
      if (e$mode == "smooth") paste0(", smoothed in climate: the loess (span 0.5) of the country slopes on the ",
                                     "countries' 1970-2023 ", spec$annual, ", at the country's long-run value (plus the ",
                                     "shift in the paths and scenarios, so that the exposure follows the climate)")
      else "; the country's own slope (the smooth at the country's climate for a country without one)",
      ". The deviations are in absolute units (C, m)."))
  }
  quad <- m$quad %||% c(TM = "none", RR = "none")
  for (v in names(quad)[quad != "none"]) {
    s <- if (v == "TM") "T" else "P"
    z <- switch(quad[[v]], cur = sprintf("%s_{it}", s), lag1 = sprintf("%s_{i,t-1}", s),
                mean20 = sprintf("\\frac{1}{20}\\sum_{k=1}^{20} %s_{i,t-k}", s),
                mean30 = sprintf("\\frac{1}{30}\\sum_{k=1}^{30} %s_{i,t-k}", s))
    out <- c(out, sprintf(paste0(
      "Quadratic control of %s: \\(Z_{%s,it} = %s - \\operatorname{med}_j \\bar %s_j\\) (%s, %s, less the ",
      "median over the countries of their %d-%d means) and its square, a control whatever the interactions, not ",
      "in the responses, paths or scenarios%s."),
      s, s, z, s, if (startsWith(quad[[v]], "mean")) paste0(QUAD_SHORT[[quad[[v]]]], ", without year \\(t\\)")
      else paste(QUAD_SHORT[[quad[[v]]]], "value"),
      if (v == "TM") "C" else "m a year", start, end,
      if (quad[[v]] == "cur") "; with this year's value it shares the variation of the deviation itself" else ""))
  }
  c(out, paste0(
    switch(effects, both = "\\(\\lambda_t\\): year effects; \\(\\alpha_i\\): country effects",
           unit = "No year effects; \\(\\alpha_i\\): country effects",
           year = "\\(\\lambda_t\\): year effects; no country effects or trends"),
    if (!no_unit) switch(m$fe, unit = "", linear = "; \\(\\delta_i t\\): country linear trends",
                         quadratic = "; \\(\\delta_i t + \\eta_i t^{2}\\): country linear and quadratic trends"),
    if (!no_unit && m$fe != "unit") " (\\(t\\) = year - 1990)", control_text(m$controls %||% CONTROL_CHOICES),
    "; \\(\\varepsilon_{it}\\) clustered by country. Sample: ", outlier_text(m$outliers %||% DEFAULT_OUTLIERS, m$gdp), "."))
}
# The controls and the outliers of a model, in words.
control_text <- function(controls) {
  parts <- c(conflict = "conflict (minor, war)", soviet = "post-Soviet collapse years",
             yugoslav = "post-Yugoslav collapse years")[controls]
  if (!length(parts)) return("; no controls")
  paste0("; \\(C_{it}\\): ", paste(parts, collapse = ", "), " dummies")
}
outlier_text <- function(outliers, gdp = NULL) {
  if (outliers == "none") return("every growth outlier kept")
  sprintf(paste("the growth outliers whose cause does not involve weather dropped (those in %s",
                "of %s x ERA5 from 1970)"), outlier_rule_words(outliers), gdp %||% "the GDP source's")
}

# Legend labels wrapped so that long setting labels stay inside the chart.
wrap_label <- function(x, width = 75) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "")
}

theme_explorer <- function() {
  theme_minimal(base_size = 13) +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA),
          panel.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = INK_MUTED),
          axis.title = element_text(colour = INK_MUTED, size = 11),
          plot.subtitle = element_text(colour = INK_MUTED, size = 10.5),
          legend.position = "top", legend.justification = "left", legend.direction = "vertical",
          legend.text = element_text(colour = INK_MUTED, size = 10.5))
}

# Inputs of setting n: its label, its
# model and its country.
setting_ui <- function(id, n) {
  ns <- NS(id)
  tagList(
    textInput(ns("label"), "Name", placeholder = paste("Setting", n)),
    if (n > 1) actionButton(ns("copy"), "Copy the model of setting 1",
                            class = "btn-sm btn-outline-secondary mb-2"),
    helpText("A setting is shown (charts, coefficients, scenarios) when a shown country on the Single country",
             "page uses it."),
    accordion(
      id = ns("model_sections"), open = "Common", multiple = TRUE,
      accordion_panel("Common", common_ui(ns)),
      accordion_panel("Temperature", variable_ui(ns, "T")),
      accordion_panel("Precipitation", variable_ui(ns, "P"))
    )
  )
}

# The options common to both variables: data, outcome, the regression design,
# sample, characteristics and the paths.
common_ui <- function(ns) {
  tagList(
    selectInput(ns("gdp"), "GDP data", c("PWT", "WB")),
    selectInput(ns("climate"), "Climate data", c("ERA5", "CRU TS", "UDelaware")),
    selectInput(ns("outcome"), "Outcome (capital and GDP per worker: PWT only)", OUTCOME_CHOICES, selected = "log"),
    radioButtons(ns("weight"), "Climate weights (estimation and CMIP6 downscaling)",
                 c("Population density" = "population", "Area (plain area average)" = "area"),
                 selected = "population", inline = TRUE),
    selectInput(ns("baseline"), "Regression unit", BASELINE_CHOICES, selected = "20-year"),
    conditionalPanel(
      "input.baseline == 'Level'", ns = ns,
      helpText("T and P themselves (in SD units over the country's horizon SD). Demeaned like the deviations,",
               "the level is the same as the horizon deviation: it differs as it is (symmetric, or without",
               "unit effects) and in first differences (Kalkuhl and Wenz).")
    ),
    conditionalPanel(
      "input.baseline == 'NRK filter'", ns = ns,
      helpText("The deviation is the innovation of the pooled AR(3) filter with lags times the country's",
               "mean climate (country effects, no year effects, fitted from 1950), not demeaned. The",
               "country SDs are those of the innovation.")
    ),
    radioButtons(ns("style"), "Regression style (first difference of the level: Kalkuhl and Wenz)",
                 STYLE_CHOICES, selected = "simple"),
    checkboxInput(ns("use_lags"), paste("Lags: every deviation term also at t-1, ..., t-L, each interacted with the",
                                        "characteristics of its own year"), value = FALSE),
    conditionalPanel(
      "input.use_lags", ns = ns,
      numericInput(ns("lags"), sprintf("Number of lags L (1-%d)", EXPLORER$max_lags %||% 20), 5,
                   min = 1, max = EXPLORER$max_lags %||% 20, step = 1)
    ),
    radioButtons(ns("effects"), "Fixed effects", EFFECTS_CHOICES, selected = "both"),
    conditionalPanel(
      "input.effects != 'year'", ns = ns,
      selectInput(ns("fe"), "Unit trends (the deviations are demeaned in the same way)",
                  FE_CHOICES, selected = "quadratic")
    ),
    radioButtons(ns("convergence"), "Convergence control (on its own, whatever the interactions)",
                 CONVERGENCE_CHOICES, selected = "none"),
    conditionalPanel(
      "input.convergence != 'none'", ns = ns,
      checkboxInput(ns("convergence_iv"), paste("Instrument it with its value one year earlier (2SLS; for last",
                                                "year's log GDP: log GDP two years before, as Harding et al. 2025)"),
                    value = FALSE)
    ),
    checkboxGroupInput(ns("controls"), "Controls", CONTROL_CHOICES, selected = CONTROL_CHOICES, inline = TRUE),
    checkboxInput(ns("drop_outliers"), paste("Drop the growth outliers whose cause is not weather (annotated",
                                             "country-years in the tails of the GDP source's growth reversals",
                                             "or growth from 1970)"),
                  value = TRUE),
    conditionalPanel(
      "input.drop_outliers", ns = ns,
      radioButtons(ns("outlier_tail"), paste("Rule and tails (reversal: the change of growth from the year",
                                             "before, flagged in the second year)"), OUTLIER_CHOICES,
                   selected = DEFAULT_OUTLIERS)
    ),
    radioButtons(ns("moderators"), "Characteristics in the interactions (estimation)",
                 c("Fixed: sample means" = "fixed",
                   "Trailing: climate of the deviation's baseline, income of the 5 years before, with level terms" =
                     "trailing"), selected = "trailing"),
    conditionalPanel(
      "input.moderators == 'trailing' && input.gdp == 'WB'", ns = ns,
      checkboxInput(ns("share_trailing"), "Trailing agricultural share too (5 years before)", value = FALSE),
      uiOutput(ns("share_warning"))
    ),
    # The low-income threshold, in the GDP source's units (default: the World
    # Bank's line, set for the source when it changes).
    conditionalPanel(
      "(input.inter_T || []).indexOf('low') > -1 || (input.inter_P || []).indexOf('low') > -1", ns = ns,
      numericInput(ns("low_threshold"), "Low-income threshold (GDP per capita)",
                   round(EXPLORER$low_income[["PWT"]]), min = 1, step = 50),
      uiOutput(ns("low_threshold_note"))
    ),
    conditionalPanel(
      "['5-year', '10-year', '20-year', '30-year'].indexOf(input.baseline) > -1 && (input.units_T == 'SD' || input.units_P == 'SD')", ns = ns,
      radioButtons(ns("sd_mode"), "SD of the trailing baseline in the GDP path",
                   c("Fixed at the country's SD" = "fixed",
                     "Moving (widened by the shift)" = "moving"), selected = "fixed")
    ),
    conditionalPanel(
      "['5-year', '10-year', '20-year', '30-year'].indexOf(input.baseline) > -1 || input.moderators == 'trailing'", ns = ns,
      radioButtons(ns("clim_mode"),
                   "Long-run characteristics of the interactions (GDP paths and climate scenarios)",
                   c("Fixed at the sample mean" = "fixed",
                     "Follow the trailing baseline (climate; in the scenarios also income and share along the SSP)" =
                       "baseline"), selected = "baseline")
    )
  )
}

# The options of the terms of one variable ("T" or "P"), each on its own.
variable_ui <- function(ns, v) {
  name <- if (v == "T") "temperature" else "precipitation"
  id <- function(x) ns(paste0(x, "_", v))
  tagList(
    radioButtons(id("units"), "Deviations in",
                 stats::setNames(c("SD", "absolute"), c("SD of the country's variability",
                                                         if (v == "T") "Absolute units (C)" else "Absolute units (m)")),
                 selected = "SD", inline = TRUE),
    selectInput(id("form"), "Arms", names(ARM_SUFFIXES), selected = "Quadratic arms"),
    checkboxInput(id("symmetric"), paste("Symmetric: one coefficient for both sides (the deviations unsplit:",
                                         "x, x², or both; not demeaned)"), value = FALSE),
    checkboxGroupInput(id("inter"), "Interactions", MODIFIERS,
                       selected = DEFAULT_MODIFIERS[[paste0("inter_", v)]], inline = TRUE),
    radioButtons(id("shape"), "Shape of the interactions",
                 c("Linear" = "linear",
                   "Restricted cubic splines (knots at the 5, 35, 65, 95th percentiles of the countries)" =
                     "splines"), selected = "linear"),
    checkboxInput(id("cap"), paste("Cap the income, share and capital interactions: they can bring the",
                                   "non-climate response to zero but never create a benefit"),
                  value = FALSE),
    # Default: T scaled with xi = 0.6 (the preferred specification), P standard.
    radioButtons(id("income_form"), "Income scaling",
                 INCOME_FORMS, selected = if (v == "T") "dollars" else "standard"),
    conditionalPanel(
      sprintf("input.income_form_%s != 'standard'", v), ns = ns,
      sliderInput(id("xi"), "Elasticity xi of the scaled terms",
                  min = XI_RANGE[1], max = XI_RANGE[2], value = if (v == "T") 0.6 else 0, step = 0.1),
      helpText("The dollar loss per person of a scaled term is proportional to Y^xi (xi = 0: equal dollars;",
               "xi = 1: equal %, the standard terms).")
    ),
    selectInput(id("quad"), "Quadratic climate control (Z + Z², on its own)", QUAD_CHOICES, selected = "none"),
    if (v == "T") exposure_T_ui(ns) else exposure_P_ui(ns)
  )
}
exposure_T_ui <- function(ns) {
  tagList(
    selectInput(ns("exposure_mode"), "Exposure", EXPOSURE_MODES, selected = "none"),
    if (!exposure_available("TM")) {
      helpText("No exposure histograms: run Rscript apps/exposure_explorer/build_exposure.R gadm0.")
    },
    conditionalPanel(
      "input.exposure_mode != 'none'", ns = ns,
      helpText(paste("Every temperature term (each arm and its ticked interactions, income scaling included)",
                     "times the country's exposure slope of each index (the change of the index per 1 C of annual",
                     "mean T), with no plain term; deviations in absolute units. Smooth: the slope as a function of",
                     "the long-run climate, which follows the warming in the paths and scenarios.")),
      radioButtons(ns("exposure_channels"), "Channels", EXPOSURE_CHANNELS, selected = "both", inline = TRUE),
      tags$div(class = "d-flex gap-2",
        selectInput(ns("exp_hot_src"), "Hot: daily series", EXPOSURE_SOURCES, selected = "tmpmax"),
        numericInput(ns("exp_hot_value"), "above (C)", 29, step = 1)),
      radioButtons(ns("exp_hot_measure"), NULL, EXPOSURE_MEASURES, selected = "dd", inline = TRUE),
      tags$div(class = "d-flex gap-2",
        selectInput(ns("exp_cold_src"), "Cold: daily series", EXPOSURE_SOURCES, selected = "tmpmin"),
        numericInput(ns("exp_cold_value"), "below (C)", 0, step = 1)),
      radioButtons(ns("exp_cold_measure"), NULL, EXPOSURE_MEASURES, selected = "dd", inline = TRUE)
    )
  )
}
exposure_P_ui <- function(ns) {
  tagList(
    selectInput(ns("exposure_mode_P"), "Exposure", EXPOSURE_MODES, selected = "none"),
    conditionalPanel(
      "input.exposure_mode_P != 'none'", ns = ns,
      helpText(paste("Every precipitation term times the country's exposure slope of the wet and dry indices",
                     "(their change per m of annual P). National daily P is a population-weighted average, which",
                     "rarely exceeds fixed thresholds in large countries: the default wet index is relative, the mm",
                     "beyond the country's 99th percentile of daily P.")),
      radioButtons(ns("exposure_channels_P"), "Channels", EXPOSURE_VARS$RR$channels, selected = "both", inline = TRUE),
      radioButtons(ns("exp_wet_type"), "Wet: threshold", c("Country percentile" = "pct", "Absolute (mm a day)" = "abs"),
                   selected = "pct", inline = TRUE),
      conditionalPanel("input.exp_wet_type == 'pct'", ns = ns,
        selectInput(ns("exp_wet_pct"), "daily P above the country's", EXPOSURE_P_PCT, selected = 0.99)),
      conditionalPanel("input.exp_wet_type == 'abs'", ns = ns,
        selectInput(ns("exp_wet_abs"), "daily P above (mm)", EXPOSURE_P_EDGES[EXPOSURE_P_EDGES >= 1], selected = 20)),
      radioButtons(ns("exp_wet_measure"), NULL, EXPOSURE_P_MEASURES, selected = "dd", inline = TRUE),
      selectInput(ns("exp_dry_value"), "Dry: daily P below (mm)", EXPOSURE_P_EDGES[EXPOSURE_P_EDGES > 0 & EXPOSURE_P_EDGES <= 10],
                  selected = 1),
      radioButtons(ns("exp_dry_measure"), NULL, c("Days" = "days", "mm short" = "dd"), selected = "days", inline = TRUE)
    )
  )
}

# Country k of the Single country page: whether it is shown, its name, the
# setting whose model it uses, and its characteristics.
country_ui <- function(id, k) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      checkboxInput(ns("active"), "Show this country", value = k == 1),
      textInput(ns("label"), "Name", placeholder = paste("Country", k)),
      selectInput(ns("setting"), "Model: setting", stats::setNames(seq_len(N_SETTINGS), paste("Setting", seq_len(N_SETTINGS))),
                  selected = 1),
      selectInput(ns("preset"), "Start from a country", choices = NULL)
    ),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      numericInput(ns("mean_T"), "Long-run mean temperature (°C)", 25, step = 0.5),
      numericInput(ns("sd_T"), "SD of annual temperature (°C)", 0.5, min = 0.01, step = 0.05),
      numericInput(ns("mean_P"), "Long-run mean precipitation (mm a year)", 1000, min = 0, step = 50),
      numericInput(ns("sd_P"), "SD of annual precipitation (mm)", 150, min = 1, step = 10)
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      tagList(numericInput(ns("income"), "GDP per capita", 5000, min = 100, step = 500), uiOutput(ns("income_units"))),
      numericInput(ns("share"), "Agricultural share of value added (%)", 15, min = 0, max = 100, step = 1),
      tagList(numericInput(ns("capital"), "Capital stock per capita (PWT)", 20000, min = 100, step = 1000),
              helpText(paste0("Long-run (sample mean) level, in ", EXPLORER$capital_units, "; used with the capital ",
                              "interaction.")))
    )
  )
}
country_tab <- function(k) {
  nav_panel(
    title = tags$span(tags$span("●", style = paste0("color:", SETTING_COLOURS[k], ";")), paste("Country", k)),
    value = paste0("c", k),
    country_ui(paste0("c", k), k)
  )
}

# The values a country preset fills in (rounded as shown in the inputs; capital
# NA without PWT capital).
preset_values <- function(row) {
  c(mean_T = round(row$mean_T, 2), sd_T = round(row$sd_T, 3), mean_P = round(1000 * row$mean_P),
    sd_P = round(1000 * row$sd_P), income = signif(row$income, 3), share = round(100 * row$share, 1),
    capital = signif(row$capital, 3))
}

# Setting n: its presets, the copy of setting 1's model, and its
# specification (label, model, country).
setting_server <- function(id, n, default_country, first = NULL) {
  moduleServer(id, function(input, output, session) {
    panel_key <- reactive({
      req(input$gdp, input$climate, input$baseline)
      paste(input$gdp, input$climate, input$baseline, input$weight %||% "population",
            input$units_T %||% "SD", sep = " | ")
    })
    # The countries of the setting's panel (the presets of the countries that
    # use it, on the Single country page).
    units <- reactive(EXPLORER$units[[panel_key()]])
    # The PWT outcomes and the capital interaction only with PWT GDP.
    observeEvent(input$gdp, {
      pwt <- input$gdp == "PWT"
      modifiers <- if (pwt) MODIFIERS else MODIFIERS[!MODIFIERS %in% PWT_MODIFIERS]
      # (The ticked ones kept; none ticked is NULL, and stays none.)
      for (k in c("inter_T", "inter_P")) {
        updateCheckboxGroupInput(session, k, choices = modifiers, inline = TRUE,
                                 selected = intersect(isolate(input[[k]]) %||% character(0), modifiers))
      }
      outcomes <- if (pwt) OUTCOME_CHOICES else OUTCOME_CHOICES[!OUTCOME_CHOICES %in% PWT_OUTCOMES]
      current <- isolate(input$outcome) %||% "log"
      updateSelectInput(session, "outcome", choices = outcomes,
                        selected = if (current %in% outcomes) current else "log")
      # The low-income threshold of the source (its units differ).
      updateNumericInput(session, "low_threshold", value = round(EXPLORER$low_income[[input$gdp]]))
    })
    output$low_threshold_note <- renderUI({
      req(input$gdp)
      helpText(sprintf(paste0("In %s. Default: the World Bank's low-income line, US$%s (2015 US$)%s. The dummy ",
                              "uses the income of the interactions (trailing or long-run)."),
                       EXPLORER$income_units[[input$gdp]], format(LOW_INCOME_LINE, big.mark = ","),
                       if (input$gdp == "PWT") sprintf(", times the median PWT / WB ratio under it: %s PPP US$",
                                                       format(round(EXPLORER$low_income[["PWT"]]), big.mark = ","))
                       else ""))
    })
    # Load a model (a list as model_inputs()) into the panels, and optionally a
    # name. The inputs some others reset when they change (the low-income
    # threshold with the GDP source, the long-run characteristics with the
    # moderators, the units with the exposure) are set again twice, after
    # those observers have run.
    pending <- reactiveVal(NULL)
    load_model <- function(m, label = NULL) {
      apply_model_inputs(session, m)
      if (!is.null(label)) updateTextInput(session, "label", value = label)
      pending(list(m = m, rounds = 2L, at = Sys.time() + 0.6))
    }
    observe({
      p <- pending()
      req(p)
      wait <- as.numeric(difftime(p$at, Sys.time(), units = "secs"))
      if (wait > 0) {
        invalidateLater(ceiling(1000 * wait) + 10)
        return()
      }
      apply_model_inputs(session, p$m)
      pending(if (p$rounds > 1) list(m = p$m, rounds = p$rounds - 1L, at = Sys.time() + 0.6))
    })
    if (!is.null(first)) {
      observeEvent(input$copy, {
        load_model(first$model_inputs())
      })
    }
    # Exposure terms are in the index per C or m: absolute deviations.
    observeEvent(input$exposure_mode, {
      if (input$exposure_mode != "none") updateRadioButtons(session, "units_T", selected = "absolute")
    }, ignoreInit = TRUE)
    observeEvent(input$exposure_mode_P, {
      if (input$exposure_mode_P != "none") updateRadioButtons(session, "units_P", selected = "absolute")
    }, ignoreInit = TRUE)
    # The characteristics estimated as trailing means are those the country has
    # as its climate changes: by default they follow the baseline in the paths.
    observeEvent(input$moderators, {
      updateRadioButtons(session, "clim_mode",
                         selected = if (input$moderators == "trailing") "baseline" else "fixed")
    }, ignoreInit = TRUE)
    # The sample of the trailing share against the trailing income alone.
    output$share_warning <- renderUI({
      req(input$gdp == "WB")
      r <- trailing_rows(input$gdp, input$climate, input$weight %||% "population", EXPLORER,
                         if (isTRUE(input$drop_outliers)) input$outlier_tail %||% DEFAULT_OUTLIERS else "none")
      div(class = "alert alert-warning small py-1 px-2",
          sprintf(paste("The trailing share needs all 5 years of WDI value added before each year:",
                        "%s rows and %d countries, against %s rows and %d countries without it (-%.0f%%),",
                        "mostly before 1990."),
                  format(r["share", "rows"], big.mark = ","), r["share", "countries"],
                  format(r["trailing", "rows"], big.mark = ","), r["trailing", "countries"],
                  100 * (1 - r["share", "rows"] / r["trailing", "rows"])))
    })

    model_inputs <- reactive(list(gdp = input$gdp, climate = input$climate, outcome = input$outcome,
                                  baseline = input$baseline, style = input$style,
                                  weight = input$weight, units_T = input$units_T, units_P = input$units_P, fe = input$fe,
                                  effects = input$effects, convergence = input$convergence,
                                  convergence_iv = isTRUE(input$convergence_iv),
                                  quad_T = input$quad_T, quad_P = input$quad_P,
                                  form_T = input$form_T, form_P = input$form_P,
                                  symmetric_T = input$symmetric_T, symmetric_P = input$symmetric_P,
                                  use_lags = input$use_lags, lags = input$lags,
                                  inter_T = input$inter_T, inter_P = input$inter_P,
                                  low_threshold = input$low_threshold,
                                  moderators = input$moderators, share_trailing = input$share_trailing,
                                  income_form_T = input$income_form_T, income_form_P = input$income_form_P,
                                  xi_T = input$xi_T, xi_P = input$xi_P,
                                  exposure_mode = input$exposure_mode, exposure_channels = input$exposure_channels,
                                  exp_hot_src = input$exp_hot_src, exp_hot_value = input$exp_hot_value,
                                  exp_hot_measure = input$exp_hot_measure, exp_cold_src = input$exp_cold_src,
                                  exp_cold_value = input$exp_cold_value, exp_cold_measure = input$exp_cold_measure,
                                  exposure_mode_P = input$exposure_mode_P, exposure_channels_P = input$exposure_channels_P,
                                  exp_wet_type = input$exp_wet_type, exp_wet_pct = input$exp_wet_pct,
                                  exp_wet_abs = input$exp_wet_abs, exp_wet_measure = input$exp_wet_measure,
                                  exp_dry_value = input$exp_dry_value, exp_dry_measure = input$exp_dry_measure,
                                  sd_mode = input$sd_mode, clim_mode = input$clim_mode,
                                  cap_T = input$cap_T, cap_P = input$cap_P, shape_T = input$shape_T, shape_P = input$shape_P,
                                  controls = input$controls,
                                  drop_outliers = input$drop_outliers, outlier_tail = input$outlier_tail))
    spec <- reactive({
      req(input$gdp, input$climate, input$baseline, input$form_T, input$form_P)
      form <- c(TM = input$form_T, RR = input$form_P)
      # The PWT outcomes and capital only with PWT (model_options()); an
      # unticked group is NULL: no interaction for that variable.
      checked <- model_options(list(gdp = input$gdp, outcome = input$outcome %||% "log",
                                    inter = list(TM = intersect(MODIFIERS, input$inter_T %||% character(0)),
                                                 RR = intersect(MODIFIERS, input$inter_P %||% character(0))),
                                    income_form = c(TM = input$income_form_T %||% "standard",
                                                    RR = input$income_form_P %||% "standard"),
                                    xi = c(TM = input$xi_T %||% 0, RR = input$xi_P %||% 0)))
      outcome <- checked$outcome
      inter <- checked$inter
      income <- checked
      uses_kap <- "kap" %in% unlist(inter)
      uses_low <- "low" %in% unlist(inter)
      req(!uses_low || isTRUE(input$low_threshold > 0))
      low_threshold <- if (uses_low) input$low_threshold
      effects <- input$effects %||% "both"
      # No unit trends without unit effects.
      fe <- if (effects == "year") "unit" else input$fe %||% "quadratic"
      convergence <- input$convergence %||% "none"
      convergence_iv <- convergence != "none" && isTRUE(input$convergence_iv)
      quad <- c(TM = input$quad_T %||% "none", RR = input$quad_P %||% "none")
      u <- units()
      trailing <- input$baseline %in% TRAILING_BASELINES
      # Per variable: symmetric terms.
      symmetric <- c(TM = isTRUE(input$symmetric_T), RR = isTRUE(input$symmetric_P))
      # The lags of the regressors (none unless ticked; 1 to the build's maximum).
      lags <- if (isTRUE(input$use_lags)) {
        req(is.finite(input$lags))
        as.integer(min(max(round(input$lags), 1), EXPLORER$max_lags %||% 20))
      } else 0L
      weight <- input$weight %||% "population"
      # Exposure of the T terms (R/exposure.R): thresholds on the histogram's
      # bin edges (whole degrees within the source's range).
      exposure_T <- if (!identical(input$exposure_mode %||% "none", "none")) {
        validate(need(exposure_available("TM"), "No exposure histograms: run Rscript apps/exposure_explorer/build_exposure.R gadm0"))
        index <- function(src, value, measure) {
          lo <- exposure_hist(src)$lo
          validate(need(is.finite(value) && value == round(value) && value %in% lo[-1],
                        sprintf("The %s threshold must be a whole degree from %g to %g C.",
                                EXPOSURE_SOURCE_SHORT[[src]], lo[2], max(lo))))
          list(src = src, type = "abs", value = value, measure = measure)
        }
        list(mode = input$exposure_mode, channels = input$exposure_channels %||% "both",
             hi = index(input$exp_hot_src %||% "tmpmax", input$exp_hot_value, input$exp_hot_measure %||% "dd"),
             lo = index(input$exp_cold_src %||% "tmpmin", input$exp_cold_value, input$exp_cold_measure %||% "dd"))
      }
      exposure_P <- if (!identical(input$exposure_mode_P %||% "none", "none")) {
        validate(need(exposure_available("RR"), "No exposure histograms: run Rscript apps/exposure_explorer/build_exposure.R gadm0"))
        wet_type <- input$exp_wet_type %||% "pct"
        list(mode = input$exposure_mode_P, channels = input$exposure_channels_P %||% "both",
             hi = list(src = "pre", type = wet_type, measure = input$exp_wet_measure %||% "dd",
                       value = as.numeric(if (wet_type == "pct") input$exp_wet_pct %||% 0.99 else input$exp_wet_abs %||% 20)),
             lo = list(src = "pre", type = "abs", value = as.numeric(input$exp_dry_value %||% 1),
                       measure = input$exp_dry_measure %||% "days"))
      }
      exposure <- if (!is.null(exposure_T) || !is.null(exposure_P)) {
        list(TM = exposure_T %||% list(mode = "none"), RR = exposure_P %||% list(mode = "none"))
      }
      # Per variable: the units of the deviations (absolute with exposure).
      units_v <- c(TM = if (!is.null(exposure_T)) "absolute" else input$units_T %||% "SD",
                   RR = if (!is.null(exposure_P)) "absolute" else input$units_P %||% "SD")
      absolute <- units_v == "absolute"
      moderators <- input$moderators %||% "fixed"
      share_trailing <- moderators == "trailing" && input$gdp == "WB" && isTRUE(input$share_trailing)
      # Absolute deviations have no SD (the paths keep it fixed for them).
      sd_mode <- if (trailing && !all(absolute)) input$sd_mode %||% "fixed" else "fixed"
      # With the horizon baseline the climate stays at the horizon mean, but
      # trailing income and share can still follow the SSP.
      clim_mode <- if (trailing || moderators == "trailing") input$clim_mode %||% "fixed" else "fixed"
      # Per variable: the cap (only with an income, share or capital
      # interaction) and the shape (splines only with an interaction).
      capped <- function(v, ticked) isTRUE(ticked) && any(c("inc", "agr", "kap", "low") %in% inter[[v]])
      cap <- c(TM = capped("TM", input$cap_T), RR = capped("RR", input$cap_P))
      shape <- c(TM = if (length(inter$TM)) input$shape_T %||% "linear" else "linear",
                 RR = if (length(inter$RR)) input$shape_P %||% "linear" else "linear")
      # A per-variable option in the label: "x" for both, "x (T)" for one.
      which_vars <- function(on, text) {
        if (all(on)) paste0(" · ", text) else if (any(on)) paste0(" · ", text, " (", c(TM = "T", RR = "P")[on], ")")
      }
      # Controls (none ticked: none) and the growth outliers dropped.
      controls <- intersect(CONTROL_CHOICES, input$controls %||% character(0))
      outliers <- if (isTRUE(input$drop_outliers)) input$outlier_tail %||% DEFAULT_OUTLIERS else "none"
      # The label names the options only where they differ from the default.
      auto <- paste0(input$gdp, " × ", input$climate, " · ", baseline_label(input$baseline, input$style %||% "simple"),
                     " · ", form_label(form, symmetric), if (lags > 0) sprintf(" · %d lag%s", lags, if (lags > 1) "s" else ""),
                     " · ", interaction_label(inter),
                     if (uses_low) paste0(" (below ", format(round(low_threshold), big.mark = ","), ")"),
                     if (outcome != "log") paste0(" · outcome: ", OUTCOME_SHORT[[outcome]]),
                     if (!is.null(income_form_label(income$income_form, income$xi))) {
                       paste0(" · ", income_form_label(income$income_form, income$xi))
                     },
                     if (!is.null(exposure)) paste0(" · ", exposure_label(exposure_options(exposure))),
                     if (effects != "both") paste0(" · ", EFFECTS_SHORT[[effects]]),
                     if (effects != "year" && fe != "quadratic") paste0(" · ", FE_SHORT[[fe]]),
                     if (convergence != "none") paste0(" · convergence: ", CONVERGENCE_SHORT[[convergence]],
                                                       if (convergence_iv) " (IV)"),
                     if (any(quad != "none")) paste0(" · quadratic control: ", paste(QUAD_SHORT[quad[quad != "none"]],
                                                    c("T", "P")[quad != "none"], collapse = ", ")),
                     if (moderators == "trailing") " · trailing characteristics",
                     if (share_trailing) " (with share)",
                     if (weight == "area") " · area weights",
                     which_vars(absolute, "absolute units"),
                     which_vars(shape == "splines", "splines"),
                     which_vars(cap, "capped"),
                     if (sd_mode == "moving") " · moving SD",
                     if (clim_mode == "baseline") " · characteristics follow baseline",
                     if (!is.null(sample_label(controls, outliers))) paste0(" · ", sample_label(controls, outliers)))
      # The user's name of the setting (empty: none; the charts of the
      # countries then say "Setting n").
      user_label <- trimws(input$label %||% "")
      list(n = n, label = if (nzchar(user_label)) user_label else auto, user_label = user_label, uses_kap = uses_kap,
           model = explorer_model(list(gdp = input$gdp, climate = input$climate, weight = weight,
                                       units = units_v,
                                       baseline = input$baseline, exposure = exposure,
                                       style = input$style %||% "simple",
                                       form = form, symmetric = symmetric, lags = lags, inter = inter, shape = shape, fe = fe,
                                       effects = effects, convergence = convergence,
                                       convergence_iv = convergence_iv, quad = quad,
                                       moderators = moderators, share_trailing = share_trailing,
                                       outcome = outcome, income_form = income$income_form, xi = income$xi,
                                       controls = controls, outliers = outliers, low_threshold = low_threshold),
                                  EXPLORER),
           sd_mode = sd_mode, clim_mode = clim_mode, cap = cap,
           # The support of the deviations the model uses, each variable's from
           # its own units and demeaning (symmetric or without unit effects:
           # as they are).
           support = do.call(rbind, lapply(c("TM", "RR"), function(v) {
             key <- paste(input$gdp, input$climate, input$baseline, weight, units_v[[v]],
                          deviation_fe(fe, symmetric[[v]], effects),
                          if (identical(input$style, "diff")) "diff", sep = " | ")
             s <- EXPLORER$support[[key]]
             s[s$variable == v, , drop = FALSE]
           })),
           units = u, gdp = input$gdp)
    })
    list(spec = spec, model_inputs = model_inputs, units = units, gdp = reactive(input$gdp),
         load = load_model, label = reactive(trimws(input$label %||% "")))
  })
}

# Country k of the Single country page: its characteristics, the setting
# whose model it uses, and its name. Presets from that setting's panel.
country_server <- function(id, k, default_country, settings) {
  moduleServer(id, function(input, output, session) {
    setting_n <- reactive(as.integer(input$setting %||% k))
    setting <- reactive(settings[[setting_n()]])
    units <- reactive(setting()$units())
    observeEvent(units(), {
      u <- units()
      current <- isolate(input$preset)
      selected <- if (is.null(current)) default_country else current
      if (!selected %in% c("", u$GID_0)) selected <- ""
      choices <- c("Custom" = "", stats::setNames(u$GID_0, u$country)[order(u$country)])
      updateSelectInput(session, "preset", choices = choices, selected = selected)
    })
    observeEvent(list(input$preset, units()), {
      req(input$preset)
      row <- units()[units()$GID_0 == input$preset, ]
      req(nrow(row) == 1)
      v <- preset_values(row)
      for (z in names(v)[is.finite(v)]) updateNumericInput(session, z, value = v[[z]])
    }, ignoreInit = TRUE)
    output$income_units <- renderUI({
      helpText(paste0("Long-run (sample mean) level, in ", EXPLORER$income_units[[setting()$gdp()]], "."))
    })
    # (Country 1 shown until its box is read.)
    active <- reactive(isTRUE(input$active %||% (k == 1)))
    # The country: its characteristics (P in m a year, share as a fraction),
    # its preset's code (the exposure slopes; NULL: custom) and its name.
    country <- reactive({
      req(is.finite(input$mean_T), isTRUE(input$sd_T > 0), is.finite(input$mean_P),
          isTRUE(input$sd_P > 0), isTRUE(input$income > 0), isTRUE(input$share >= 0))
      u <- units()
      row <- u[u$GID_0 == (input$preset %||% ""), ]
      values <- c(mean_T = input$mean_T, sd_T = input$sd_T, mean_P = input$mean_P,
                  sd_P = input$sd_P, income = input$income, share = input$share)
      name <- if (nrow(row) == 1) {
        if (isTRUE(all.equal(unname(values), unname(preset_values(row)[names(values)])))) row$country
        else paste(row$country, "(edited)")
      } else "custom"
      list(GID_0 = if (nrow(row) == 1) row$GID_0, mean_T = input$mean_T, sd_T = input$sd_T,
           mean_P = input$mean_P / 1000, sd_P = input$sd_P / 1000, income = input$income,
           share = input$share / 100, capital = if (isTRUE(input$capital > 0)) input$capital,
           name = name, user_label = trimws(input$label %||% ""))
    })
    list(active = active, setting_n = setting_n, country = country)
  })
}

setting_tab <- function(n) {
  nav_panel(
    title = tags$span(tags$span("●", style = paste0("color:", SETTING_COLOURS[n], ";")),
                      paste("Setting", n)),
    value = paste0("s", n),
    setting_ui(paste0("s", n), n)
  )
}

# Snapshots of a setting's model and the literature presets.
source(file.path(APP_DIR, "R", "snapshots.R"), local = environment())

ui <- page_navbar(
  title = "Weather response explorer",
  id = "page",
  # Scroll rather than squeeze the charts into the window's height (a short
  # window would leave them too small to draw).
  fillable = FALSE,
  theme = bs_theme(version = 5, bg = SURFACE, fg = INK, primary = SETTING_COLOURS[1]),
  # The settings (model and country), common to every page.
  sidebar = sidebar(
    width = 380,
    downloadButton("save_page", "Save the page as HTML", class = "btn-sm btn-outline-secondary"),
    helpText("The shown settings, charts, tables and notes; the climate scenarios and the termination",
             "shocks only if run (their last run)."),
    tags$hr(),
    tags$div(
      tags$strong("Snapshots of a model"),
      selectInput("snap_setting", "Setting", stats::setNames(seq_len(N_SETTINGS), paste("Setting", seq_len(N_SETTINGS))),
                  selected = 1),
      tags$div(class = "d-flex gap-2 align-items-end",
               textInput("snap_name", "Name", placeholder = "e.g. my baseline", width = "100%"),
               actionButton("snap_save", "Save", class = "btn-sm btn-outline-secondary mb-3")),
      selectInput("snap_pick", "Load (literature presets and saved snapshots)", snapshot_choices(list())),
      uiOutput("snap_note"),
      tags$div(class = "d-flex gap-2",
               actionButton("snap_load", "Load into the setting", class = "btn-sm btn-outline-primary"),
               actionButton("snap_delete", "Delete", class = "btn-sm btn-outline-danger")),
      helpText("A snapshot keeps every model option of the setting (not its countries); saved ones are kept in",
               "results/response_explorer/snapshots.rds.")
    ),
    tags$hr(),
    tags$div(
      tags$strong("Snapshots of the app"),
      tags$div(class = "d-flex gap-2 align-items-end",
               textInput("app_snap_name", "Name", placeholder = "e.g. draft figures", width = "100%"),
               actionButton("app_snap_save", "Save", class = "btn-sm btn-outline-secondary mb-3")),
      selectInput("app_snap_pick", "Load a saved app snapshot", character(0)),
      tags$div(class = "d-flex gap-2",
               actionButton("app_snap_load", "Load", class = "btn-sm btn-outline-primary"),
               actionButton("app_snap_delete", "Delete", class = "btn-sm btn-outline-danger")),
      helpText("Every setting (model and name), every country, the options of every page, the fitted models",
               "and the projections, and whether the climate scenarios and termination shocks were run (they",
               "are run again on loading, from the saved projections). Kept in",
               "results/response_explorer/app_snapshots/.")
    ),
    tags$hr(),
    do.call(navset_pill, c(lapply(seq_len(N_SETTINGS), setting_tab), list(id = "setting_tabs")))
  ),
  # The warnings of the settings, above every page.
  header = uiOutput("warnings"),
  nav_panel(
    "Single country",
    card(
      card_header("Shift (common to all settings)"),
      layout_columns(
        col_widths = c(2, 2, 3, 2, 3),
        numericInput("shock_T", "Temperature (°C)", 1, step = 0.5),
        numericInput("shock_P", "Precipitation (mm a year)", 100, step = 25),
        sliderInput("horizon", "Years shown in the GDP paths", 10, 100, 50, step = 5),
        numericInput("ramp_years", "Year X", 50, min = 1, step = 5),
        checkboxGroupInput("gdp_lines", "Shifts shown in the GDP paths",
                           c("Immediate (solid)" = "permanent", "Gradual until year X (dashed)" = "gradual",
                             "Sudden in year X (dotted)" = "delayed",
                             "Transitory: year 0 only (dot-dash; the response over the lags)" = "one-off"),
                           selected = "gradual")
      ),
      helpText("Year X: the gradual shift reaches its full size there (dashed) and the sudden late one starts",
               "(dotted).")
    ),
    card(
      card_header("Countries: each with its characteristics and the setting whose model it uses"),
      do.call(navset_pill, c(lapply(seq_len(N_SETTINGS), country_tab), list(id = "country_tabs")))
    ),
    checkboxInput("response_first_lag",
                  "Distributed lags: show the year-0 coefficient (the effect in the year of the anomaly) instead of the sum over the lags",
                  FALSE),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Temperature: growth response to an annual anomaly"),
           plotOutput("response_T", height = 380, hover = hoverOpts("hover_response_T", delay = 80)),
           uiOutput("readout_response_T")),
      card(card_header("Precipitation: growth response to an annual anomaly"),
           plotOutput("response_P", height = 380, hover = hoverOpts("hover_response_P", delay = 80)),
           uiOutput("readout_response_P"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header(textOutput("gdp_header_T", inline = TRUE)),
           plotOutput("gdp_T", height = 430, hover = hoverOpts("hover_gdp_T", delay = 80)),
           uiOutput("readout_gdp_T")),
      card(card_header(textOutput("gdp_header_P", inline = TRUE)),
           plotOutput("gdp_P", height = 430, hover = hoverOpts("hover_gdp_P", delay = 80)),
           uiOutput("readout_gdp_P"))
    )
  ),
  nav_panel(
    "By characteristic",
    card(
      card_header("Effect of a year k SD above and below the baseline, against a characteristic of the country"),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        selectInput("mod_x_T", "Temperature: characteristic on the horizontal axis", MODERATION_CHOICES,
                    selected = "mean_T"),
        selectInput("mod_x_P", "Precipitation: characteristic on the horizontal axis", MODERATION_CHOICES,
                    selected = "mean_P"),
        selectInput("mod_ref", "Other characteristics at (and the SD of)", choices = c("Country 1" = "1")),
        numericInput("mod_k", "Size of the anomaly, k (SD of the reference country)", 1, min = 0.1, step = 0.5)
      ),
      plotOutput("mod_plot_T", height = 380),
      plotOutput("mod_plot_P", height = 380),
      uiOutput("mod_notes")
    ),
    card(
      card_header("Heatmap: the effect over two characteristics at once"),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        radioButtons("heat_var", "Variable", c("Temperature" = "TM", "Precipitation" = "RR"), selected = "TM",
                     inline = TRUE),
        radioButtons("heat_sign", "Anomaly", c("+k SD (warmer / wetter)" = "1", "-k SD (cooler / drier)" = "-1"),
                     selected = "1"),
        selectInput("heat_x", "Horizontal axis", MODERATION_CHOICES, selected = "mean_T"),
        selectInput("heat_y", "Vertical axis", MODERATION_CHOICES, selected = "income")
      ),
      plotOutput("heat_plot", height = 520),
      helpText("Colour: the effect on growth of one year at the anomaly (red: a loss, blue: a gain), over the range of",
               "the setting's countries on both axes, the other characteristics (and k, the SD) at the reference",
               "country above (circle). Dots: the setting's countries, where the data are; away from them the map",
               "extrapolates. Point estimates (no CI); with exposure, the smooth at each climate.")
    ),
    card(
      card_header("Physical exposure: how the indices change with the annual value, by the country's climate"),
      plotOutput("exposure_plot", height = 460),
      uiOutput("exposure_notes")
    )
  ),
  nav_panel(
    "Coefficients",
    card(card_header("Each setting's model: specification, coefficients, fit and joint tests"),
         withMathJax(uiOutput("spec_math")),
         tableOutput("coef_table"),
         uiOutput("coef_notes"))
  ),
  nav_panel(
    "Climate scenarios",
    card(
      card_header("Climate scenarios: GDP losses under each setting's model, 2025-2100"),
      layout_columns(
        col_widths = c(6, 6),
        checkboxGroupInput("scen_sel", "Forcing: the global warming path", SCENARIO_CHOICES,
                           selected = SCENARIO_CHOICES, inline = TRUE),
        radioButtons("econ_ssp", paste("Economic development: income, population and agricultural",
                                       "share (Follow SSP-RCP: each forcing with its own SSP)"),
                     ECONOMY_CHOICES, selected = "follow", inline = TRUE)
      ),
      layout_columns(
        col_widths = c(6, 3),
        selectInput("scen_noise", "Natural variability (each country's detrended historical T and P)",
                    VARIABILITY_CHOICES, selected = "none"),
        selectInput("scen_noise_n", "Noise paths", NOISE_PATH_CHOICES, selected = "200"),
        checkboxInput("scen_no_gains", NO_GAINS_LABEL, value = FALSE)
      ),
      layout_columns(
        col_widths = c(4, 3, 5),
        selectInput("map_setting", "Maps: setting", choices = NULL),
        selectInput("map_scenario", "Maps: forcing", choices = NULL),
        radioButtons("map_component", "GDP maps: impact of",
                     c("T and P" = "loss", "T only" = "loss_T", "P only" = "loss_P"),
                     selected = "loss", inline = TRUE)
      ),
      div(class = "d-flex align-items-center gap-3 mb-2",
          actionButton("run_scenarios", "Run the climate scenarios", class = "btn-primary"),
          uiOutput("scenario_status")),
      radioButtons("scen_gmt_measure", "Warming chart", WARMING_CHOICES, selected = "level", inline = TRUE),
      radioButtons("global_measure", "Global chart", MEASURE_CHOICES, selected = "loss", inline = TRUE),
      region_input("scen_region"),
      checkboxInput("scen_common_y", COMMON_Y_LABEL, value = FALSE),
      conditionalPanel("input.global_measure == 'growth'",
        checkboxInput("scen_history", "Show the historical growth of world GDP, or of the selected countries' (PWT 11.0, 1971-2023: series, median, IQR, range)", value = TRUE)),
      layout_columns(
        col_widths = c(4, 8),
        plotOutput("gmt_plot", height = 380),
        plotOutput("global_plot", height = 380, hover = hoverOpts("hover_global", delay = 80))
      ),
      uiOutput("readout_global"),
      layout_columns(
        col_widths = c(6, 6),
        plotOutput("map_2050", height = 340),
        plotOutput("map_2100", height = 340)
      ),
      tags$h6("Physical drivers: the change of each country's climate in the scenario", class = "mt-3"),
      layout_columns(
        col_widths = c(6, 6),
        plotOutput("driver_T_2050", height = 320),
        plotOutput("driver_T_2100", height = 320)
      ),
      layout_columns(
        col_widths = c(6, 6),
        plotOutput("driver_P_2050", height = 320),
        plotOutput("driver_P_2100", height = 320)
      ),
      tableOutput("global_table"),
      uiOutput("scenario_notes")
    )
  ),
  nav_panel(
    "Geoengineering",
    card(
      card_header("Geoengineering and termination shocks, 2025-2100"),
      layout_columns(
        col_widths = c(3, 3, 6),
        selectInput("ts_forcing", "Background forcing (the warming without geoengineering)",
                    SCENARIO_CHOICES, selected = "ssp245"),
        tagList(
          radioButtons("ts_mode", "Geoengineering limits", c("Global warming (C above 1850-1900)" = "level",
                                                              "Rate of warming (C per decade)" = "rate"),
                       selected = "level", inline = TRUE),
          conditionalPanel("input.ts_mode != 'rate'",
            numericInput("ts_held", "Held at (C; 0: preindustrial)", 1.5, min = 0, step = 0.1)),
          conditionalPanel("input.ts_mode == 'rate'",
            numericInput("ts_rate", "At most (C per decade; 0: warming stops)", 0.1, min = 0, step = 0.05))
        ),
        checkboxGroupInput("ts_years", "Termination years (none: only sustained geoengineering)",
                           as.character(seq(2040, 2090, by = 10)), selected = c("2040", "2060", "2080"), inline = TRUE)
      ),
      layout_columns(
        col_widths = c(3, 9),
        numericInput("ts_tau", "Catch-up time after termination, tau (years)", 4, min = 0.5, step = 0.5),
        radioButtons("ts_econ", "Economic development (Follow SSP-RCP: the forcing's own SSP)",
                     ECONOMY_CHOICES, selected = "follow", inline = TRUE)
      ),
      if (!is.null(CIDER)) tagList(
        layout_columns(
          col_widths = c(6, 3),
          radioButtons("ts_sai", "Regional pattern of geoengineering",
                       c("Uniform: each country's CMIP6 pattern of global warming (perfect geoengineering)" = "uniform",
                         "CIDER: stratospheric aerosol injection at the latitudes below" = "cider"),
                       selected = "uniform"),
          conditionalPanel(
            "input.ts_sai == 'cider' || input.ts_compare == 'sai'",
            selectInput("ts_cider_model", "CIDER emulated model",
                        stats::setNames(names(CIDER$models), unname(CIDER$models))))
        ),
        conditionalPanel(
          "input.ts_sai == 'cider' || input.ts_compare == 'sai'",
          div(class = "d-flex flex-wrap gap-2 align-items-end",
              tags$span(class = "small me-2", "Injection fractions by latitude (non-negative, scaled to sum to one):"),
              lapply(CIDER_LATS, function(l) {
                div(style = "width: 6.5em;",
                    numericInput(paste0("ts_f_", l), l, CIDER_FRACTION_DEFAULTS[[l]], min = 0, step = 0.05))
              })),
          uiOutput("ts_sai_info")
        )
      ),
      layout_columns(
        col_widths = c(3, 9),
        radioButtons("ts_compare", "Compare side by side", TS_COMPARE_CHOICES, selected = "none"),
        tagList(
          conditionalPanel(
            "input.ts_compare == 'sai'",
            checkboxGroupInput("ts_cmp_sai", "Geoengineering patterns (in place of the one above; CIDER with the model and fractions above)",
                               TS_SAI_CHOICES, selected = c("uniform", "mix"), inline = TRUE)
          ),
          conditionalPanel(
            "input.ts_compare == 'forcing'",
            checkboxGroupInput("ts_cmp_forcings", "Background forcings (in place of the one above)", SCENARIO_CHOICES,
                               selected = intersect(c("ssp126", "ssp245", "ssp585"), SCENARIO_CHOICES), inline = TRUE)
          ),
          conditionalPanel(
            "input.ts_compare == 'economy'",
            checkboxGroupInput("ts_cmp_economies", "Economic developments (in place of the one above)",
                               ECONOMY_CHOICES[-1], selected = c("SSP1", "SSP2", "SSP3"), inline = TRUE)
          ),
          conditionalPanel(
            "input.ts_compare == 'tau'",
            textInput("ts_cmp_taus", "Catch-up times in years, comma separated (in place of tau above)", "2, 4, 8")
          ),
          conditionalPanel(
            "input.ts_compare != 'none'",
            helpText("The warming and global charts show one column per alternative; the map, the latitude chart",
                     "and the exceedances show the one chosen below them.")
          )
        )
      ),
      layout_columns(
        col_widths = c(6, 3),
        selectInput("ts_noise", "Natural variability (each country's detrended historical T and P)",
                    VARIABILITY_CHOICES, selected = "none"),
        selectInput("ts_noise_n", "Noise paths", NOISE_PATH_CHOICES, selected = "200"),
        checkboxInput("ts_no_gains", NO_GAINS_LABEL, value = FALSE)
      ),
      div(class = "d-flex align-items-center gap-3 mb-2",
          actionButton("run_termination", "Run the geoengineering scenarios", class = "btn-primary"),
          uiOutput("ts_status")),
      radioButtons("ts_gmt_measure", "Warming chart", WARMING_CHOICES, selected = "level", inline = TRUE),
      radioButtons("ts_measure", "Global chart", MEASURE_CHOICES, selected = "loss", inline = TRUE),
      region_input("ts_region"),
      checkboxInput("ts_common_y", COMMON_Y_LABEL, value = FALSE),
      conditionalPanel("input.ts_measure == 'growth'",
        checkboxInput("ts_history", "Show the historical growth of world GDP, or of the selected countries' (PWT 11.0, 1971-2023: series, median, IQR, range)", value = TRUE)),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        selectInput("ts_ref", "Exceedances: reference climate", EXCEEDANCE_REFERENCES, selected = "20"),
        numericInput("ts_k_gmt", "Warming chart: band above k SD, k =", 3, min = 0.5, step = 0.5),
        numericInput("ts_k_lat", "Latitude chart: events beyond k SD, k =", 2.5, min = 0.5, step = 0.5),
        radioButtons("ts_sign", "Latitude chart: events", c("Warm (above +k SD)" = "warm",
                                                            "Either sign (beyond +-k SD)" = "both"),
                     selected = "warm")
      ),
      uiOutput("ts_charts"),
      uiOutput("readout_ts"),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("ts_map_setting", "Map and exceedances: setting (its panel's climate data)", choices = NULL),
        selectInput("ts_map_year", "Map: termination year", choices = NULL),
        selectInput("ts_variant", "Map, latitude chart and exceedances: alternative", choices = NULL)
      ),
      plotOutput("ts_map", height = 380),
      if (!is.null(CIDER)) tagList(
        radioButtons("ts_sai_var", "Map of the regional pattern: climate left by geoengineering that offsets 1 C of global warming",
                     c("Temperature (C)" = "T", "Precipitation (%)" = "P"), selected = "T", inline = TRUE),
        plotOutput("ts_sai_map", height = 380)
      ),
      selectInput("ts_lat_paths", "Latitude chart: warming paths", choices = NULL, multiple = TRUE,
                  width = "100%"),
      plotOutput("ts_lat_plot", height = 420),
      tableOutput("ts_table"),
      uiOutput("ts_notes")
    )
  ),
  nav_panel("Models and assumptions", card(card_header("Models and assumptions"), uiOutput("notes"))),
  nav_panel("Guide", guide_ui())
)

server <- function(input, output, session) {
  settings <- list(setting_server("s1", 1, DEFAULT_COUNTRIES[1]))
  settings <- c(settings, lapply(2:N_SETTINGS, function(n) {
    setting_server(paste0("s", n), n, DEFAULT_COUNTRIES[n], first = settings[[1]])
  }))

  # Snapshots: save the model of a setting under a name, load a preset or a
  # saved snapshot into a setting. The setting follows the open settings tab.
  saved <- reactiveVal(read_snapshots())
  observeEvent(input$setting_tabs, {
    updateSelectInput(session, "snap_setting", selected = sub("^s", "", input$setting_tabs))
  })
  refresh_snapshots <- function(selected = NULL) {
    updateSelectInput(session, "snap_pick", choices = snapshot_choices(saved()), selected = selected)
  }
  observeEvent(TRUE, refresh_snapshots(), once = TRUE)
  observeEvent(input$snap_save, {
    s <- settings[[as.integer(input$snap_setting)]]
    name <- trimws(input$snap_name %||% "")
    if (!nzchar(name)) name <- s$label()
    if (!nzchar(name)) name <- format(Sys.time(), "Snapshot %Y-%m-%d %H:%M")
    all <- saved()
    replaced <- name %in% names(all)
    all[[name]] <- list(name = name, note = sprintf("Saved %s from setting %s.", format(Sys.time(), "%Y-%m-%d %H:%M"),
                                                    input$snap_setting),
                        inputs = s$model_inputs(), saved = Sys.time())
    write_snapshots(all)
    saved(all)
    refresh_snapshots(paste0("own:", name))
    showNotification(sprintf("%s snapshot \"%s\".", if (replaced) "Replaced" else "Saved", name), type = "message")
  })
  output$snap_note <- renderUI({
    x <- snapshot_of(input$snap_pick, saved())
    req(x)
    helpText(x$note)
  })
  observeEvent(input$snap_load, {
    x <- snapshot_of(input$snap_pick, saved())
    req(x)
    n <- as.integer(input$snap_setting)
    settings[[n]]$load(x$inputs, label = x$name)
    showNotification(sprintf("Loaded \"%s\" into setting %d.", x$name, n), type = "message")
  })
  observeEvent(input$snap_delete, {
    if (!startsWith(input$snap_pick %||% "", "own:")) {
      showNotification("Only saved snapshots can be deleted (not the literature presets).", type = "warning")
      return()
    }
    all <- saved()
    all[[sub("^own:", "", input$snap_pick)]] <- NULL
    write_snapshots(all)
    saved(all)
    refresh_snapshots()
  })

  # Snapshots of the app: every input, the settings' fitted models and the
  # cached projections; on loading the inputs are set three times (some
  # observers reset others when they change: the presets of a country, the
  # threshold with the GDP source, ...), then the runs that had been made are
  # made again (from the saved projections).
  scen_trigger <- reactiveVal(0)
  ts_trigger <- reactiveVal(0)
  refresh_app_snapshots <- function(selected = NULL) {
    updateSelectInput(session, "app_snap_pick", choices = app_snapshot_names(), selected = selected)
  }
  observeEvent(TRUE, refresh_app_snapshots(), once = TRUE)
  observeEvent(input$app_snap_save, {
    name <- trimws(input$app_snap_name %||% "")
    if (!nzchar(name)) name <- format(Sys.time(), "App %Y-%m-%d %H%M")
    ids <- setdiff(names(input), grep(APP_SNAPSHOT_SKIP, names(input), value = TRUE))
    values <- lapply(stats::setNames(ids, ids), function(id) isolate(input[[id]]))
    models <- list()
    # (The models of the shown settings: those a shown country uses.)
    for (s in settings[isolate(used_settings())]) {
      m <- tryCatch(isolate(s$spec())$model, error = function(e) NULL)
      if (!is.null(m$key)) models[[m$key]] <- m
    }
    x <- list(name = gsub("[^A-Za-z0-9 _.-]", "_", name), saved = Sys.time(), inputs = values,
              settings = lapply(settings, function(s) list(inputs = isolate(s$model_inputs()), label = isolate(s$label()))),
              runs = c(scenarios = isolate(scen_runs()) > 0, termination = isolate(ts_runs()) > 0),
              models = models, projections = as.list(projection_cache))
    replaced <- file.exists(app_snapshot_file(x$name))
    write_app_snapshot(x)
    refresh_app_snapshots(x$name)
    showNotification(sprintf("%s app snapshot \"%s\" (%s MB).", if (replaced) "Replaced" else "Saved", x$name,
                             format(round(file.size(app_snapshot_file(x$name)) / 2^20, 1))), type = "message")
  })
  app_pending <- reactiveVal(NULL)
  set_inputs <- function(values) {
    for (id in names(values)) {
      v <- values[[id]]
      if (is.null(v)) v <- character(0)
      if (is.numeric(v) && length(v) == 1 && !is.finite(v)) next
      session$sendInputMessage(id, list(value = v))
    }
  }
  observeEvent(input$app_snap_load, {
    x <- read_app_snapshot(input$app_snap_pick %||% "")
    if (is.null(x)) {
      showNotification("No such app snapshot (or the file cannot be read).", type = "error")
      return()
    }
    for (k in names(x$models)) model_cache[[k]] <- x$models[[k]]
    for (k in names(x$projections)) projection_cache[[k]] <- x$projections[[k]]
    for (n in seq_along(x$settings)) if (n <= N_SETTINGS) {
      settings[[n]]$load(x$settings[[n]]$inputs, label = x$settings[[n]]$label)
    }
    set_inputs(x$inputs)
    load_step <<- 0L
    app_pending(list(x = x, at = Sys.time()))
    showNotification(sprintf("Loading app snapshot \"%s\" (saved %s)...", x$name, format(x$saved, "%Y-%m-%d %H:%M")),
                     type = "message")
  })
  # Every 0.8 s after loading: the inputs again (twice), then the runs.
  load_step <- 0L
  observe({
    p <- app_pending()
    req(p)
    invalidateLater(800)
    load_step <<- load_step + 1L
    if (load_step == 1L) return()
    if (load_step <= 3L) {
      set_inputs(p$x$inputs)
      return()
    }
    if (isTRUE(p$x$runs[["scenarios"]])) scen_trigger(isolate(scen_trigger()) + 1)
    if (isTRUE(p$x$runs[["termination"]])) ts_trigger(isolate(ts_trigger()) + 1)
    app_pending(NULL)
  })
  observeEvent(input$app_snap_delete, {
    f <- app_snapshot_file(input$app_snap_pick %||% "")
    if (!nzchar(input$app_snap_pick %||% "") || !file.exists(f)) return()
    unlink(f)
    refresh_app_snapshots()
    showNotification(sprintf("Deleted app snapshot \"%s\".", input$app_snap_pick), type = "message")
  })

  # The countries of the Single country page, each with the setting it uses.
  countries <- lapply(seq_len(N_SETTINGS), function(k) {
    country_server(paste0("c", k), k, DEFAULT_COUNTRIES[k], settings)
  })
  shown_k <- reactive({
    k <- which(vapply(countries, function(cn) cn$active(), TRUE))
    validate(need(length(k) > 0, "Tick \"Show this country\" for at least one country (Single country page, Countries card)."))
    k
  })
  # The shown settings: those a shown country uses (in setting order), with
  # unique labels.
  used_settings <- reactive(sort(unique(vapply(countries[shown_k()], function(cn) cn$setting_n(), 0L))))
  specs <- reactive({
    out <- lapply(settings[used_settings()], function(s) s$spec())
    labels <- vapply(out, `[[`, "", "label")
    dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    for (i in which(dup)) out[[i]]$label <- paste0(labels[i], " (setting ", out[[i]]$n, ")")
    out
  })
  # The pairs of the single-country charts: each shown country with the
  # specification of its setting, labelled "Country k - Setting n" with the
  # user's names where given, and coloured by the country.
  pairs <- reactive({
    out <- lapply(shown_k(), function(k) {
      cn <- countries[[k]]
      s <- settings[[cn$setting_n()]]$spec()
      country <- cn$country()
      validate(need(!s$uses_kap || !is.null(country$capital),
                    sprintf("Country %d: its setting has the capital interaction; enter its capital per capita.", k)))
      country_name <- if (nzchar(country$user_label)) country$user_label else paste("Country", k)
      setting_name <- if (nzchar(s$user_label)) s$user_label else paste("Setting", s$n)
      s$country <- country
      s$setting_label <- s$label
      s$label <- paste(country_name, "-", setting_name)
      s$colour_n <- k
      s
    })
    # Unique labels (the user's names may repeat).
    labels <- vapply(out, `[[`, "", "label")
    dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    for (i in which(dup)) out[[i]]$label <- paste0(labels[i], " (country ", out[[i]]$colour_n, ")")
    out
  })
  # The colour of each shown country's curves by its label (in country order).
  colours <- reactive({
    pr <- pairs()
    stats::setNames(SETTING_COLOURS[vapply(pr, `[[`, 0, "colour_n")], vapply(pr, `[[`, "", "label"))
  })
  as_setting <- function(x) factor(x, levels = names(colours()))
  shift <- function(v) if (v == "TM") input$shock_T else input$shock_P / 1000
  to_display <- function(v) if (v == "TM") 1 else 1000
  sd_of <- function(s, v) if (v == "TM") s$country$sd_T else s$country$sd_P

  # Responses over -3..3 SD of each setting's country (wider if the shift is
  # larger), in physical units, flagged where the deviation is within the 1st-
  # 99th percentile of those in the data.
  responses <- function(v) {
    reactive({
      delta <- shift(v)
      req(is.finite(delta))
      do.call(rbind, lapply(pairs(), function(s) {
        sd <- sd_of(s, v)
        reach <- max(3, 1.1 * abs(delta) / sd)
        x <- seq(-reach, reach, length.out = 241)
        # The deviations in the model's units (SD, or C and m).
        x_model <- x * sd / deviation_scale(s$model, s$country, v)
        q <- s$support[s$support$variable == v, ]
        d <- response_curve(s$model, v, x_model, s$country, s$cap, isTRUE(input$response_first_lag))
        d$x <- x
        d$anomaly <- x * sd * to_display(v)
        d$in_support <- x_model >= q$q01 & x_model <= q$q99
        d$setting <- s$label
        d
      }))
    })
  }
  responses_T <- responses("TM")
  responses_P <- responses("RR")
  shift_points <- function(v) {
    do.call(rbind, lapply(pairs(), function(s) {
      d <- response_curve(s$model, v, shift(v) / deviation_scale(s$model, s$country, v), s$country, s$cap,
                          isTRUE(input$response_first_lag))
      d$anomaly <- shift(v) * to_display(v)
      d$setting <- s$label
      d
    }))
  }

  # The three scenarios of the shift: from year 0 on (solid), rising linearly
  # to its full size in year X (dashed), or nothing until year X and the full
  # shift from it on (dotted); the names are the types of gdp_path().
  # Only the ticked ones are drawn.
  scenarios <- reactive({
    all <- c(permanent = "immediate",
             gradual = sprintf("gradual until year %d", as.integer(input$ramp_years)),
             delayed = sprintf("sudden in year %d", as.integer(input$ramp_years)),
             "one-off" = "transitory (year 0 only)")
    shown <- all[names(all) %in% input$gdp_lines]
    validate(need(length(shown) > 0, "Tick at least one shift (immediate, gradual, sudden or transitory) to draw."))
    shown
  })
  # The line type and band transparency of each scenario, by its type.
  SCENARIO_LINETYPE <- c(permanent = "solid", gradual = "42", delayed = "13", "one-off" = "4212")
  SCENARIO_ALPHA <- c(permanent = 0.15, gradual = 0.07, delayed = 0.07, "one-off" = 0.07)
  paths <- function(v) {
    reactive({
      delta <- shift(v)
      req(is.finite(delta), input$horizon, isTRUE(input$ramp_years >= 1))
      do.call(rbind, lapply(pairs(), function(s) {
        do.call(rbind, lapply(names(scenarios()), function(type) {
          cbind(setting = s$label, scenario = scenarios()[[type]],
                gdp_path(s$model, v, s$country, type, delta, input$horizon, s$sd_mode,
                         s$clim_mode, s$cap, ramp_years = input$ramp_years))
        }))
      }))
    })
  }
  paths_T <- paths("TM")
  paths_P <- paths("RR")

  response_plot <- function(d, v) {
    d$setting <- as_setting(d$setting)
    points <- shift_points(v)
    points$setting <- as_setting(points$setting)
    unit <- if (v == "TM") "°C" else "mm"
    ggplot(d, aes(anomaly, 100 * estimate, colour = setting, fill = setting)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_vline(xintercept = 0, colour = GRID, linewidth = 0.4) +
      geom_ribbon(aes(ymin = 100 * low, ymax = 100 * high), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.6, linetype = "22") +
      geom_line(aes(y = ifelse(in_support, 100 * estimate, NA)), linewidth = 0.9, na.rm = TRUE) +
      geom_point(data = points, size = 3.2, shape = 21, fill = SURFACE, stroke = 1.4) +
      scale_colour_manual(values = colours(), name = NULL, labels = wrap_label) +
      scale_fill_manual(values = colours(), guide = "none") +
      labs(x = paste0("Annual anomaly against the baseline (", unit, ")"),
           y = if (any(vapply(pairs(), function(s) (s$model$lags %||% 0L) > 0, TRUE))) {
             if (isTRUE(input$response_first_lag)) "Effect on growth that year (pp; year-0 coefficient)" else "Effect on growth (pp; lags summed)"
           } else "Effect on growth that year (pp)",
           subtitle = paste0("Bands: 95% CI. Dashed: beyond the 1st-99th percentile of the ",
                             "deviations in the data. Points: the shift.")) +
      theme_explorer()
  }
  output$response_T <- renderPlot(response_plot(responses_T(), "TM"), bg = SURFACE)
  output$response_P <- renderPlot(response_plot(responses_P(), "RR"), bg = SURFACE)

  readout_list <- function(lines) tags$ul(class = "small mb-0", lapply(lines, tags$li))
  response_readout <- function(d, hover, v) {
    if (is.null(hover)) return(helpText("Hover over the chart to read the values."))
    unit <- if (v == "TM") "°C" else "mm"
    readout_list(lapply(split(d, as_setting(d$setting)), function(s) {
      i <- which.min(abs(s$anomaly - hover$x))
      sprintf("%s: %+.2f %s (%+.2f SD): %+.2f pp [%+.2f, %+.2f]", s$setting[i], s$anomaly[i], unit,
              s$x[i], 100 * s$estimate[i], 100 * s$low[i], 100 * s$high[i])
    }))
  }
  output$readout_response_T <- renderUI(response_readout(responses_T(), input$hover_response_T, "TM"))
  output$readout_response_P <- renderUI(response_readout(responses_P(), input$hover_response_P, "RR"))

  gdp_header <- function(v) {
    sprintf("Effect on %s of a %+g %s shift in %s: immediate, gradual or sudden in year %d, or for one year",
            quantity_of(pairs()), if (v == "TM") input$shock_T else input$shock_P, if (v == "TM") "°C" else "mm",
            if (v == "TM") "temperature" else "precipitation", as.integer(input$ramp_years))
  }
  output$gdp_header_T <- renderText(gdp_header("TM"))
  output$gdp_header_P <- renderText(gdp_header("RR"))

  gdp_plot <- function(d) {
    d$setting <- as_setting(d$setting)
    d$scenario <- factor(d$scenario, levels = unname(scenarios()))
    x <- as.integer(input$ramp_years)
    captions <- c(permanent = "Solid: the full shift from year 0 on.",
                  gradual = sprintf("Dashed: rising linearly to the full shift in year %d.", x),
                  delayed = sprintf("Dotted: nothing until year %d, the full shift from then on.", x),
                  "one-off" = "Dot-dash: the shift in year 0 only (without lags a step; with lags it moves for L years).")
    ggplot(d, aes(year, pct, colour = setting, fill = setting, linetype = scenario,
                  group = interaction(setting, scenario))) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = pct_low, ymax = pct_high, alpha = scenario), colour = NA,
                  show.legend = FALSE) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = colours(), name = NULL, labels = wrap_label) +
      scale_fill_manual(values = colours(), guide = "none") +
      scale_linetype_manual(values = stats::setNames(SCENARIO_LINETYPE[names(scenarios())], scenarios()),
                            name = NULL) +
      scale_alpha_manual(values = stats::setNames(SCENARIO_ALPHA[names(scenarios())], scenarios()),
                         guide = "none") +
      guides(colour = guide_legend(order = 1),
             linetype = guide_legend(order = 2, direction = "horizontal", nrow = 1,
                                     override.aes = list(colour = INK_MUTED))) +
      labs(x = "Years after the shift starts (0 = first year)", y = paste(first_upper(quantity_of(pairs())), "against no shift (%)"),
           subtitle = paste(strwrap(paste(c(captions[names(scenarios())], "Bands: 95% CI."),
                                          collapse = " "), width = 110), collapse = "\n")) +
      theme_explorer() +
      theme(legend.box = "vertical", legend.box.just = "left")
  }
  output$gdp_T <- renderPlot(gdp_plot(paths_T()), bg = SURFACE)
  output$gdp_P <- renderPlot(gdp_plot(paths_P()), bg = SURFACE)

  gdp_readout <- function(d, hover) {
    if (is.null(hover)) return(helpText("Hover over the chart to read the values."))
    year <- max(0, min(round(hover$x), max(d$year)))
    rows <- d[d$year == year, ]
    rows <- rows[order(as_setting(rows$setting), match(rows$scenario, scenarios())), ]
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(sprintf("%s, %s: %+.2f%% [%+.2f, %+.2f]", rows$setting, rows$scenario,
                                 rows$pct, rows$pct_low, rows$pct_high)))
  }
  output$readout_gdp_T <- renderUI(gdp_readout(paths_T(), input$hover_gdp_T))
  output$readout_gdp_P <- renderUI(gdp_readout(paths_P(), input$hover_gdp_P))

  # "By characteristic": the reference country is one of the shown countries.
  observe({
    k <- shown_k()
    names(k) <- vapply(k, function(i) {
      cn <- tryCatch(countries[[i]]$country(), error = function(e) NULL)
      lab <- if (!is.null(cn) && nzchar(cn$user_label)) cn$user_label else paste("Country", i)
      if (!is.null(cn)) paste0(lab, " (", cn$name, ")") else lab
    }, "")
    current <- isolate(input$mod_ref)
    updateSelectInput(session, "mod_ref", choices = k, selected = if (!is.null(current) && current %in% k) current else k[[1]])
  })
  mod_reference <- reactive({
    i <- as.integer(input$mod_ref %||% shown_k()[[1]])
    req(i %in% shown_k())
    countries[[i]]$country()
  })
  # One chart per variable, each against its own characteristic.
  mod_results <- function(v) {
    reactive({
      x <- input[[paste0("mod_x_", if (v == "TM") "T" else "P")]]
      k <- input$mod_k
      req(is.finite(k), k > 0, x)
      ref <- mod_reference()
      sp <- specs()
      d <- do.call(rbind, lapply(sp, function(s) {
        r <- moderation_data(s, ref, x, k, vars = v)
        if (!is.null(r)) r$setting <- s$label
        r
      }))
      validate(need(!is.null(d) && nrow(d) > 0, "No setting to show for this characteristic (capital needs PWT GDP)."))
      rugs <- do.call(rbind, lapply(sp, function(s) {
        z <- s$units[[x]]
        data.frame(setting = s$label, value = z[is.finite(z)])
      }))
      list(d = d, rugs = rugs[rugs$setting %in% d$setting, ], x = x)
    })
  }
  mod_results_T <- mod_results("TM")
  mod_results_P <- mod_results("RR")
  mod_plot_of <- function(r) {
    k <- input$mod_k
    colours <- stats::setNames(c("#c0392b", "#2a78d6"),
                               c(sprintf("+%g SD (warmer / wetter)", k), sprintf("-%g SD (cooler / drier)", k)))
    moderation_plot(r$d, r$rugs, r$x, k, colours)
  }
  output$mod_plot_T <- renderPlot(mod_plot_of(mod_results_T()), bg = SURFACE)
  output$mod_plot_P <- renderPlot(mod_plot_of(mod_results_P()), bg = SURFACE)
  # The heatmap over two characteristics.
  heat_results <- reactive({
    x <- input$heat_x; y <- input$heat_y; k <- input$mod_k
    req(x, y, is.finite(k), k > 0)
    validate(need(x != y, "Choose two different characteristics for the axes."))
    ref <- mod_reference()
    sp <- specs()
    v <- input$heat_var %||% "TM"
    sign <- as.numeric(input$heat_sign %||% "1")
    d <- do.call(rbind, lapply(sp, function(s) heat_data(s, ref, v, sign, k, x, y)))
    validate(need(!is.null(d) && nrow(d) > 0, "No setting to show for these characteristics (capital needs PWT GDP)."))
    pts <- do.call(rbind, lapply(sp, function(s) {
      u <- s$units
      ok <- is.finite(u[[x]]) & is.finite(u[[y]])
      data.frame(setting = s$label, xv = u[[x]][ok], yv = u[[y]][ok])
    }))
    refs <- do.call(rbind, lapply(sp, function(s) {
      b <- moderation_base(s, ref)
      if (is.null(b[[x]]) || is.null(b[[y]])) return(NULL)
      data.frame(setting = s$label, xv = b[[x]], yv = b[[y]])
    }))
    keep <- unique(d$setting)
    title <- sprintf("%s: a year %s%g SD of the reference country", if (v == "TM") "Temperature" else "Precipitation",
                     if (sign > 0) "+" else "-", k)
    list(d = d, pts = pts[pts$setting %in% keep, ], refs = refs[refs$setting %in% keep, , drop = FALSE],
         x = x, y = y, title = title)
  })
  output$heat_plot <- renderPlot({
    r <- heat_results()
    heat_plot(r$d, r$pts, r$refs, r$x, r$y, r$title)
  }, bg = SURFACE)
  output$mod_notes <- renderUI(helpText(
    "Each panel: the effect on growth of one year at +k or -k SD of the reference country's annual T (or P) from its",
    "baseline, as the characteristic on the horizontal axis varies over the range of the setting's countries (rug),",
    "with the other characteristics at the reference country (circle: the reference country itself, with its own",
    "exposure slopes). Along the curve the country has no preset, so exposure is the smooth at its climate. With",
    "climate on the axis, the long-run climate of the interactions and the climate of this year's and last year's",
    "value move together. 95% bands; with the cap, from draws of the coefficients."))
  exposure_results <- reactive({
    d <- exposure_curve_data(specs())
    validate(need(!is.null(d), "No shown setting has exposure-mediated terms (Temperature or Precipitation > Exposure)."))
    ref <- mod_reference()
    marks <- unique(d[c("setting", "variable", "channel")])
    marks$clim <- ifelse(marks$variable == "TM", ref$mean_T, 1000 * ref$mean_P)
    list(d = d, marks = marks)
  })
  output$exposure_plot <- renderPlot({
    r <- exposure_results()
    exposure_curve_plot(r$d, r$marks)
  }, bg = SURFACE)
  output$exposure_notes <- renderUI(helpText(
    "For each setting with exposure and each channel: every country's slope (points; the change of the index per",
    "1 C of annual mean T, or per m of annual P, from ERA5 daily data 1970-2023) and the smooth in climate (line),",
    "which the smooth version uses and which follows the climate in the paths and scenarios. Dashed: the reference",
    "country's long-run climate."))

  # The equation of each shown setting's model and the definitions of its terms.
  # The parts of the page that are more than a call are functions, used by
  # their outputs and by the saved page (save_page below).
  spec_math <- function() {
    blocks <- lapply(specs(), function(s) {
      tags$div(
        class = "mb-3",
        tags$div(tags$span("●", style = paste0("color:", SETTING_COLOURS[s$n], ";")),
                 tags$strong(paste0(" Setting ", s$n, ": ")), s$label),
        tags$div(style = "overflow-x: auto;", paste0("$$", spec_equation(s$model), "$$")),
        tags$ul(class = "small", lapply(spec_definitions(s$model), function(x) tags$li(HTML(x))))
      )
    })
    tagList(blocks)
  }
  output$spec_math <- renderUI(withMathJax(spec_math()))

  # The coefficients of each shown setting's model (rows: every term of any of
  # them; deviation terms, then the level terms of trailing characteristics),
  # its fit, the verdict against setting 1, and the joint Wald tests.
  coef_table <- function() {
    sp <- specs()
    models <- lapply(sp, `[[`, "model")
    cell <- function(est, se, p) {
      sprintf("%.3g%s<br><span style='color:%s'>(%.3g)</span>", est, stars(p), INK_MUTED, se)
    }
    coef_cells <- lapply(models, function(m) {
      se <- sqrt(diag(m$vcov))
      out <- stats::setNames(cell(m$coef, se, m$p), names(m$coef))
      if (nrow(m$levels)) out <- c(out, stats::setNames(cell(m$levels$estimate, m$levels$se, m$levels$p),
                                                        paste0("level:", m$levels$term)))
      if (!is.null(m$convergence_coef)) {
        k <- m$convergence_coef
        out <- c(out, stats::setNames(cell(k[["estimate"]], k[["se"]], k[["p"]]), paste0("convergence:", m$convergence)))
      }
      # The quadratic climate controls ("quad:<kind>:<column>").
      if (!is.null(m$quad_coef)) {
        q <- m$quad_coef
        kind <- m$quad[ifelse(grepl("TM", q$term), "TM", "RR")]
        out <- c(out, stats::setNames(cell(q$estimate, q$se, q$p), paste0("quad:", kind, ":", q$term)))
      }
      out
    })
    outcome_label <- function(m) paste0("growth of ", OUTCOME_SHORT[[m$outcome %||% "log"]], " (log points)")
    terms <- unique(unlist(lapply(coef_cells, names)))
    extra <- startsWith(terms, "level:") | startsWith(terms, "convergence:") | startsWith(terms, "quad:")
    terms <- c(terms[!extra], terms[startsWith(terms, "level:")], terms[startsWith(terms, "convergence:")],
               terms[startsWith(terms, "quad:")])
    test_names <- unique(unlist(lapply(models, function(m) m$tests$test)))
    rows <- c("Setting", "Outcome", terms, "Observations", "Countries", "Within adj. R2", "AIC", "BIC",
              "Against setting 1 (AIC / BIC / R2)", test_names)
    columns <- lapply(seq_along(models), function(i) {
      m <- models[[i]]
      tests <- stats::setNames(sprintf("p = %.3g%s<br><span style='color:%s'>(F = %.2f, %d terms)</span>",
                                       m$tests$p, stars(m$tests$p), INK_MUTED, m$tests$stat,
                                       as.integer(m$tests$df1)), m$tests$test)
      verdict <- if (i == 1) "" else record_verdicts(m, models[[1]]) %||% "other rows or outcome: not comparable"
      values <- c(Setting = htmltools::htmlEscape(sp[[i]]$label), Outcome = outcome_label(m), coef_cells[[i]],
                  Observations = format(m$nobs, big.mark = ","), Countries = as.character(m$countries),
                  `Within adj. R2` = sprintf("%.4f", m$fit[["war2"]]), AIC = sprintf("%.1f", m$fit[["aic"]]),
                  BIC = sprintf("%.1f", m$fit[["bic"]]), `Against setting 1 (AIC / BIC / R2)` = verdict, tests)
      out <- values[rows]
      out[is.na(out)] <- ""
      unname(out)
    })
    names(columns) <- vapply(sp, function(s) {
      sprintf("<span style='color:%s'>●</span> Setting %d", SETTING_COLOURS[s$n], s$n)
    }, "")
    labels <- rows
    labels[rows %in% terms] <- vapply(rows[rows %in% terms], term_label, "")
    labels[rows %in% test_names] <- paste("Joint test,", rows[rows %in% test_names])
    data.frame(Term = labels, columns, check.names = FALSE)
  }
  output$coef_table <- renderTable(coef_table(), striped = TRUE, spacing = "s", align = "l",
                                   sanitize.text.function = function(x) x)
  coef_notes <- function() tags$ul(
    class = "small",
    tags$li(paste(
      "Coefficients: effect on the outcome (growth of GDP per capita, capital per capita or GDP per",
      "worker, in log points) per unit of the arm (an SD, or a C or m, beyond the baseline; squared",
      "for quadratic arms), and per unit of the characteristic for the interactions (long-run T in C",
      "and P in m a year less the median of the countries; log income, share and log capital per",
      "capita less their means over the rows; spline terms as in the notes below).",
      "Clustered SEs in brackets; *** p < 0.01, ** p < 0.05, * p < 0.1.")),
    tags$li(paste(
      "Level terms: the trailing characteristics on their own (controls; not in the responses",
      "and projections).")),
    tags$li(paste(
      "Joint tests: clustered Wald tests (F) that the terms are all zero: all the deviation terms",
      "of T (P), the terms interacted with each characteristic (both arms, with their spline",
      "terms), each pair of characteristics when there are three or more, and all the interactions;",
      "with linear + quadratic arms, the squared arms (zero: linear arms suffice) and the linear arms",
      "(zero: quadratic arms suffice), for T, P and both.")),
    tags$li(paste(
      "AIC and BIC compare only models on the same rows. Against setting 1: better, worse, or",
      "indifferent when the difference cannot be told apart from noise (Vuong test clustered by",
      "country, 5%), for each criterion; on other rows (another panel, or trailing values that",
      "drop rows) or with another outcome the models are not compared."))
  )
  output$coef_notes <- renderUI(coef_notes())

  # Warnings for each setting: characteristics outside the countries of its
  # sample, and a shift beyond the deviations in the data.
  warnings_ui <- function() {
    msgs <- character(0)
    for (s in pairs()) {
      cn <- s$country
      u <- s$units
      m <- s$model
      prefix <- paste0(s$label, ": ")
      check <- function(value, range, label, fmt) {
        if (value < min(range) || value > max(range)) {
          msgs <<- c(msgs, sprintf(paste0("%s%s (", fmt, ") is outside the range of the sample's ",
                                          "countries (", fmt, " to ", fmt, ")."),
                                   prefix, label, value, min(range), max(range)))
        }
      }
      used <- unique(unlist(m$inter))
      if (any(OWN_MODIFIERS %in% m$inter$TM)) check(cn$mean_T, u$mean_T, "mean temperature", "%.1f °C")
      if (any(OWN_MODIFIERS %in% m$inter$RR)) check(1000 * cn$mean_P, 1000 * u$mean_P, "mean precipitation", "%.0f mm")
      if (any(c("inc", "low") %in% used)) check(cn$income, u$income, "GDP per capita", "%.0f")
      if ("agr" %in% used) check(100 * cn$share, 100 * u$share, "agricultural share", "%.1f%%")
      if ("kap" %in% used) check(cn$capital, u$capital[is.finite(u$capital)], "capital per capita", "%.0f")
      for (v in c("TM", "RR")) {
        # The shift and the data's deviations in the model's units.
        x <- shift(v) / deviation_scale(m, cn, v)
        q <- s$support[s$support$variable == v, ]
        if (x < q$q01 || x > q$q99) {
          unit <- if (!identical(of_variable(m$units, v), "absolute")) "SD" else if (v == "TM") "C" else "m"
          msgs <- c(msgs, sprintf(paste(
            "%sthe %s shift is %+.2f %s, beyond the 1st-99th percentile of the deviations in the",
            "data (%+.2f to %+.2f %s): its response is an extrapolation."),
            prefix, if (v == "TM") "temperature" else "precipitation", x, unit, q$q01, q$q99, unit))
        }
      }
      # (First differences and a level outcome: the effect of a permanent
      # shift is bounded.)
      grows <- !isTRUE(m$diff) && !level_outcome(m)
      if (grows && m$baseline %in% c("Horizon", "Level")) {
        msgs <- c(msgs, paste0(prefix, "with the ", tolower(m$baseline), " as the regression unit the permanent ",
                               "shift stays a deviation every year, so its GDP effect grows without limit."))
      }
      if (grows && identical(m$window, "nrk")) {
        # A permanent shift a leaves the innovation (1 - sum phi) a every year.
        left <- vapply(c("TM", "RR"), function(v) 1 - sum(filter_lags(m, cn, v)), 0)
        msgs <- c(msgs, sprintf(paste0(
          "%swith the NRK filter a permanent shift leaves an innovation of %.2f (T) and %.2f (P) times its ",
          "size every year after the first three (one minus the filter's persistence at the country's climate), ",
          "so its GDP effect grows without limit unless that is zero."), prefix, left[["TM"]], left[["RR"]]))
      }
      if (isTRUE(m$share_trailing)) {
        r <- trailing_rows(m$gdp, m$climate, m$weight, EXPLORER, m$outliers)
        msgs <- c(msgs, sprintf(paste(
          "%sthe trailing agricultural share needs all 5 years of WDI value added before each year:",
          "the model runs on %s rows (%d countries), %.0f%% fewer than with the trailing income alone",
          "(%s rows), mostly before 1990, so its sample differs from the other models'."),
          prefix, format(m$nobs, big.mark = ","), m$countries,
          100 * (1 - m$nobs / r["trailing", "rows"]), format(r["trailing", "rows"], big.mark = ",")))
      }
    }
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-warning", style = "margin-bottom: 0.5rem;",
        tags$ul(style = "margin-bottom: 0;", lapply(msgs, tags$li)))
  }
  output$warnings <- renderUI(warnings_ui())

  # Climate scenarios --------------------------------------------------------
  # The projections of a setting depend on its model and options only (every
  # country enters with its own characteristics), so they are cached by those.
  projection_cache <- new.env()
  # The runs: the clicks on Run plus those a loaded app snapshot asks for.
  scen_runs <- reactive(input$run_scenarios + scen_trigger())
  ts_runs <- reactive(input$run_termination + ts_trigger())
  # The natural variability of a run: its type (VARIABILITY_CHOICES) and number
  # of noise paths.
  noise_of <- function(type, n) {
    list(type = type %||% "none", n = if ((type %||% "none") == "none") 0L else as.integer(n %||% 200))
  }
  noise_key <- function(noise) paste(noise$type, noise$n)
  # The bands of the global charts, their dotted lines with noise (the range of
  # the noise paths at the estimated coefficients), and the values of rows of a
  # global loss as text.
  band_label <- function(noise) {
    if (noise$type == "none") return("bands: 95% CI of the coefficients")
    "bands: 95% of the coefficient draws, each with its noise path; dotted: 95% of the noise paths alone"
  }
  noise_lines <- function(d) {
    if (!"noise_low" %in% names(d)) return(NULL)
    list(geom_line(aes(y = noise_low), linetype = "dotted", linewidth = 0.5),
         geom_line(aes(y = noise_high), linetype = "dotted", linewidth = 0.5))
  }
  loss_text <- function(rows, measure = "loss") {
    unit <- if (measure == "growth_change") " pp" else "%"
    out <- sprintf("%+.2f%s [%+.2f, %+.2f]", rows$estimate, unit, rows$low, rows$high)
    if ("noise_low" %in% names(rows)) {
      out <- paste0(out, sprintf(", variability [%+.2f, %+.2f]", rows$noise_low, rows$noise_high))
    }
    if (measure == "growth") out <- paste0(out, sprintf("; without climate change %.2f%%", rows$baseline_growth))
    out
  }
  # The global series of the chosen measure: the loss (% of GDP against no
  # climate change), the growth of world GDP (% a year, from the second year),
  # or its change by climate change (pp a year: with less without climate
  # change), the estimate and bands in the loss's columns.
  GROWTH_NOTE <- "Growth of world GDP (the GDP-weighted countries') needs a GDP outcome, not capital."
  as_measure <- function(d, measure) {
    if (measure != "loss") {
      d <- d[d$year > PROJECTION_FROM, ]
      less <- if (measure == "growth_change") d$baseline_growth else 0
      d$estimate <- d$growth - less
      d$low <- d$growth_low - less
      d$high <- d$growth_high - less
      if ("noise_low" %in% names(d)) {
        d$noise_low <- d$growth_noise_low - less
        d$noise_high <- d$growth_noise_high - less
      }
    }
    d
  }
  measure_axis <- function(measure, sp, region = NULL) {
    name <- region_name(region)
    switch(measure, growth = sprintf("Growth of %s GDP (%% a year)", name),
           growth_change = sprintf("Change of %s GDP growth by climate change (pp a year)", name),
           paste(first_upper(name), quantity_of(sp), "against no climate change (%)"))
  }
  # The countries of a global chart (input `id`): NULL for the world.
  region_of <- function(id) {
    x <- input[[id]]
    if (length(x)) sort(unique(x)) else NULL
  }
  # The history layers and note of a growth chart of `region`.
  region_history <- function(region) if (is.null(region)) PWT_WORLD else pwt_growth(region)
  # "; aggregate of 3 of the 4 selected countries (98% of their GDP)" when a
  # region is chosen: the countries the setting projects, and their GDP share.
  region_note <- function(region, results) {
    if (is.null(region)) return("")
    n <- min(vapply(results, function(r) r$n, 0))
    coverage <- min(vapply(results, function(r) r$coverage, 0))
    sprintf("; aggregate of %d of the %d selected countries (at least %.0f%% of their SSP GDP in %d)",
            n, length(region), 100 * coverage, PROJECTION_FROM)
  }
  growth_allowed <- function(sp) !any(vapply(sp, function(s) identical(s$model$outcome, "capital"), TRUE))
  # Growth of world GDP of global series `g` of setting n of `sp` for the
  # tables: the compound mean over the years after the first, without and with
  # climate change, or (`largest_drop`) with climate change and the largest
  # fall of the growth rate by climate change (pp) and its year.
  growth_text <- function(g, sp, n, largest_drop) {
    if (!growth_allowed(Filter(function(s) s$n == n, sp))) return("-")
    g <- g[g$year > PROJECTION_FROM, ]
    mean_of <- function(x) 100 * (prod(1 + x / 100)^(1 / length(x)) - 1)
    if (!largest_drop) return(sprintf("%.2f, %.2f", mean_of(g$baseline_growth), mean_of(g$growth)))
    change <- g$growth - g$baseline_growth
    i <- which.min(change)
    sprintf("%.2f; %+.2f (%d)", mean_of(g$growth), change[i], g$year[i])
  }
  noise_label <- function(noise) {
    if (noise$type == "none") "no natural variability"
    else sprintf("%s, %d noise paths", tolower(names(VARIABILITY_CHOICES)[VARIABILITY_CHOICES == noise$type]), noise$n)
  }
  # A setting's projection with the natural variability of its panel.
  # (`sai`: list(model, fractions), geoengineering with CIDER's pattern under
  # the setting's climate weights; NULL: uniform.)
  project_with <- function(s, scenario, ssp, noise, gmt = NULL, no_gains = FALSE, region = NULL, sai = NULL) {
    m <- s$model
    variability <- EXPLORER$variability[[paste(m$gdp, m$climate, m$weight %||% "population", sep = " | ")]]
    sai_arg <- if (!is.null(sai)) {
      list(pattern = cider_mix(CIDER, sai$model, sai$fractions, m$weight %||% "population")$pattern,
           background = ts_background(scenario))
    }
    out <- tryCatch(project_setting(s, SCEN, scenario, ssp, gmt = gmt, variability = variability, noise = noise$type,
                                    n_paths = noise$n, no_gains = no_gains, region = region, sai = sai_arg),
                    error = function(e) if (is.null(region)) stop(e) else NULL)
    validate(need(!is.null(out), paste0("None of the selected countries can be projected with setting ", s$label, ".")))
    out
  }
  projection_of <- function(s, scenario, ssp, noise, no_gains = FALSE, region = NULL) {
    key <- paste(s$model$key, s$sd_mode, s$clim_mode, paste(s$cap, collapse = "/"), scenario, ssp, noise_key(noise),
                 if (no_gains) "no gains", if (!is.null(region)) paste(region, collapse = "+"), sep = " | ")
    if (is.null(projection_cache[[key]])) {
      projection_cache[[key]] <- project_with(s, scenario, ssp, noise, no_gains = no_gains, region = region)
    }
    projection_cache[[key]]
  }
  # The projections run only when asked: the shown settings, the selected
  # forcings, the economic development and the variability at the moment of
  # the click.
  NOT_RUN <- "Press \"Run the climate scenarios\" to project the global losses and the maps."
  run <- eventReactive(scen_runs(), {
    req(scen_runs() > 0)
    validate(need(length(input$scen_sel) > 0, "Select at least one forcing."))
    list(specs = specs(), scenarios = input$scen_sel, economy = input$econ_ssp %||% "follow",
         noise = noise_of(input$scen_noise, input$scen_noise_n), no_gains = isTRUE(input$scen_no_gains))
  })
  run_key <- function(sp, scenarios, economy, noise, no_gains) {
    paste(c(vapply(sp, function(s) paste(s$label, s$model$key, s$sd_mode, s$clim_mode, paste(s$cap, collapse = "/")), ""),
            scenarios, economy, noise_key(noise), no_gains), collapse = " || ")
  }
  # The label of a forcing of the last run, with its economy.
  run_label <- function(scenario) scenario_label(scenario, run()$economy)
  scenario_status <- function() {
    if (scen_runs() == 0) return(helpText(class = "mb-0", "Not run yet."))
    r <- run()
    current <- run_key(specs(), input$scen_sel, input$econ_ssp %||% "follow",
                       noise_of(input$scen_noise, input$scen_noise_n), isTRUE(input$scen_no_gains))
    if (!identical(current, run_key(r$specs, r$scenarios, r$economy, r$noise, r$no_gains))) {
      return(helpText(class = "mb-0 text-warning-emphasis",
                      "Settings or scenarios changed since the last run: press Run to update."))
    }
    helpText(class = "mb-0", sprintf("Up to date: %d setting(s) x %d forcing(s), economy: %s; %s.",
                                     length(r$specs), length(r$scenarios),
                                     names(ECONOMY_CHOICES)[ECONOMY_CHOICES == r$economy],
                                     paste0(noise_label(r$noise), if (r$no_gains) "; gains ignored in the world aggregate")))
  }
  output$scenario_status <- renderUI(scenario_status())
  # The colours of the settings of the last run (in setting order).
  run_colours <- reactive({
    sp <- run()$specs
    stats::setNames(SETTING_COLOURS[vapply(sp, `[[`, 0, "n")], vapply(sp, `[[`, "", "label"))
  })
  as_run_setting <- function(x) factor(x, levels = names(run_colours()))
  # Every setting and scenario of the last run: list of (setting, scenario, result).
  # (`region`: the countries aggregated, NULL: the world.)
  projections_for <- function(region = NULL) {
    r <- run()
    jobs <- expand.grid(i = seq_along(r$specs), scenario = r$scenarios, stringsAsFactors = FALSE)
    withProgress(message = if (is.null(region)) "Projecting the scenarios" else "Projecting the selected countries",
                 value = 0, {
      lapply(seq_len(nrow(jobs)), function(j) {
        s <- r$specs[[jobs$i[j]]]
        scenario <- jobs$scenario[j]
        incProgress(1 / nrow(jobs), detail = paste(s$label, "-", scenario_label(scenario, r$economy)))
        list(setting = s$label, n = s$n, scenario = scenario,
             label = scenario_label(scenario, r$economy),
             result = projection_of(s, scenario, economy_of(scenario, r$economy), r$noise, r$no_gains, region))
      })
    })
  }
  projections <- reactive(projections_for())
  # The global chart's: the world's, or the selected countries' aggregate.
  scen_region <- reactive(region_of("scen_region")) |> debounce(600)
  global_projections <- reactive({
    region <- scen_region()
    if (is.null(region)) projections() else projections_for(region)
  })
  global_losses <- reactive({
    do.call(rbind, lapply(global_projections(), function(p) {
      cbind(setting = p$setting, scenario = p$scenario, label = p$label, p$result$global)
    }))
  })
  # The maps offer the settings and scenarios of the last run.
  observeEvent(run(), {
    r <- run()
    choices <- stats::setNames(as.character(vapply(r$specs, `[[`, 0, "n")),
                               vapply(r$specs, `[[`, "", "label"))
    selected <- isolate(input$map_setting)
    if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
    updateSelectInput(session, "map_setting", choices = choices, selected = selected)
    scenario_choices <- SCENARIO_CHOICES[SCENARIO_CHOICES %in% r$scenarios]
    names(scenario_choices) <- vapply(scenario_choices, scenario_label, "", r$economy)
    selected <- isolate(input$map_scenario)
    if (is.null(selected) || !selected %in% scenario_choices) {
      selected <- if ("ssp245" %in% scenario_choices) "ssp245" else scenario_choices[[1]]
    }
    updateSelectInput(session, "map_scenario", choices = scenario_choices, selected = selected)
  })

  gmt_plot <- function() {
    req(length(input$scen_sel) > 0)
    g <- SCEN$gmt[SCEN$gmt$scenario %in% input$scen_sel & SCEN$gmt$year >= 1970, ]
    g$scenario <- factor(g$scenario, levels = names(SCEN$scenarios))
    rate <- identical(input$scen_gmt_measure, "rate")
    if (rate) {
      # The rate of the median path (the 5-95% range of a rate is not that of
      # the levels: not shown).
      g <- g[order(g$scenario, g$year), ]
      g$median <- unlist(lapply(split(g, g$scenario, drop = TRUE), function(x) warming_rate(x$year, x$median)))
      g <- g[is.finite(g$median), ]
    }
    ggplot(g, aes(year, median, colour = scenario, fill = scenario)) +
      geom_vline(xintercept = PROJECTION_FROM, colour = GRID, linewidth = 0.5) +
      (if (rate) geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4)
       else geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.10, colour = NA)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = SCENARIO_COLOURS, labels = FORCING_LABELS, name = NULL,
                          guide = guide_legend(direction = "horizontal", nrow = 1)) +
      scale_fill_manual(values = SCENARIO_COLOURS, guide = "none") +
      labs(x = NULL, y = if (rate) RATE_AXIS else "Global warming above 1850-1900 (°C)",
           subtitle = if (rate) wrap_label(paste("FaIR (AR6-calibrated) median:", RATE_TEXT), 70)
                      else "FaIR (AR6-calibrated): median and 5-95% range") +
      theme_explorer()
  }
  output$gmt_plot <- renderPlot(gmt_plot(), bg = SURFACE)

  global_plot <- function() {
    validate(need(scen_runs() > 0, NOT_RUN))
    measure <- input$global_measure %||% "loss"
    validate(need(measure == "loss" || growth_allowed(run()$specs), GROWTH_NOTE))
    d <- as_measure(global_losses(), measure)
    d$setting <- as_run_setting(d$setting)
    d$scenario <- factor(d$scenario, levels = names(SCEN$scenarios))
    labels <- vapply(names(SCEN$scenarios), run_label, "")
    # Growth: each scenario's growth without climate change, dashed.
    baseline <- if (measure == "growth") {
      geom_line(aes(y = baseline_growth), linetype = "22", linewidth = 0.5)
    }
    ggplot(d, aes(year, estimate, colour = scenario, fill = scenario)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.10, colour = NA) +
      noise_lines(d) +
      baseline +
      (if (measure == "growth" && isTRUE(input$scen_history)) history_growth_layers(region_history(scen_region()))) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~setting, ncol = min(3, nlevels(d$setting)), scales = y_scales(input$scen_common_y),
                 labeller = labeller(setting = function(x) wrap_label(x, 50))) +
      scale_colour_manual(values = SCENARIO_COLOURS, labels = labels, name = NULL,
                          guide = guide_legend(direction = "horizontal", nrow = 1)) +
      scale_fill_manual(values = SCENARIO_COLOURS, guide = "none") +
      labs(x = NULL, y = measure_axis(measure, run()$specs, scen_region()),
           subtitle = wrap_label(paste0(
             "Best guess (median warming, CMIP6-median pattern, estimated coefficients",
             if (run()$noise$type != "none") ", mean over the noise paths", "); ",
             band_label(run()$noise),
             region_note(scen_region(), lapply(global_projections(), `[[`, "result")),
             if (measure == "growth") paste0("; dashed: the SSP's growth without climate change",
                                             if (isTRUE(input$scen_history)) history_growth_note(region_history(scen_region()), region_name(scen_region())))), 110)) +
      theme_explorer() +
      theme(strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0))
  }
  output$global_plot <- renderPlot(global_plot(), bg = SURFACE)

  output$readout_global <- renderUI({
    hover <- input$hover_global
    if (scen_runs() == 0) return(NULL)
    if (is.null(hover)) return(helpText("Hover over the global chart to read the values."))
    measure <- input$global_measure %||% "loss"
    if (measure != "loss" && !growth_allowed(run()$specs)) return(NULL)
    d <- as_measure(global_losses(), measure)
    year <- max(min(d$year), min(round(hover$x), PROJECTION_TO))
    rows <- d[d$year == year, ]
    rows <- rows[order(as_run_setting(rows$setting), match(rows$scenario, names(SCEN$scenarios))), ]
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(paste0(rows$setting, ", ", rows$label, ": ", loss_text(rows, measure))))
  })

  # Best-guess country losses of the map setting and scenario in `year`, with
  # one colour scale for the two maps.
  map_data_of <- reactive({
    validate(need(scen_runs() > 0, NOT_RUN))
    req(input$map_setting, input$map_scenario)
    s <- Filter(function(x) as.character(x$n) == input$map_setting, run()$specs)
    req(length(s) == 1)
    r <- projection_of(s[[1]], input$map_scenario, economy_of(input$map_scenario, run()$economy), run()$noise,
                       run()$no_gains)
    list(setting = s[[1]]$label, label = run_label(input$map_scenario),
         countries = r$countries[r$countries$year %in% c(2050, 2100), ])
  })
  theme_map <- function() {
    theme_void(base_size = 12) +
      theme(plot.background = element_rect(fill = SURFACE, colour = NA),
            plot.title = element_text(colour = INK, face = "bold"),
            plot.subtitle = element_text(colour = INK_MUTED, size = 9.5),
            legend.position = "bottom", legend.key.width = grid::unit(1.6, "cm"),
            legend.title = element_text(colour = INK_MUTED, size = 9.5),
            legend.text = element_text(colour = INK_MUTED, size = 9))
  }
  # The world map with `value` of each country (iso3), Antarctica left out.
  world_with <- function(iso3, value) {
    world <- SCEN$world[SCEN$world$iso3 != "ATA", ]
    world$value <- value[match(world$iso3, iso3)]
    world
  }
  COMPONENT_TITLES <- c(loss = "impact of T and P", loss_T = "impact of T only",
                        loss_P = "impact of P only")
  # GDP impact of the chosen component (T and P, T only, P only), with one
  # colour scale for the two years.
  country_map <- function(year) {
    m <- map_data_of()
    component <- input$map_component %||% "loss"
    values <- m$countries[m$countries$year == year, ]
    limit <- max(abs(m$countries[[component]]), na.rm = TRUE)
    ggplot(world_with(values$iso3, values[[component]])) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      scale_fill_gradient2(low = MAP_LOSS, mid = MAP_ZERO, high = MAP_GAIN, midpoint = 0,
                           limits = c(-limit, limit), na.value = MAP_NA,
                           name = paste0(first_upper(quantity_of(run()$specs)), " against no\nclimate change (%)")) +
      labs(title = sprintf("%d: %s, %s", year, m$label, COMPONENT_TITLES[[component]]),
           subtitle = wrap_label(m$setting)) +
      theme_map()
  }
  output$map_2050 <- renderPlot(country_map(2050), bg = SURFACE)
  output$map_2100 <- renderPlot(country_map(2100), bg = SURFACE)

  # Physical drivers: each country's T change (C) and P change (% of its
  # long-run P) in 2050 and 2100 against the mean of the map setting's horizon,
  # the CMIP6 pattern times the global warming (as in the projections).
  driver_data <- reactive({
    validate(need(scen_runs() > 0, NOT_RUN))
    req(input$map_setting, input$map_scenario)
    s <- Filter(function(x) as.character(x$n) == input$map_setting, run()$specs)
    req(length(s) == 1)
    g <- SCEN$gmt[SCEN$gmt$scenario == input$map_scenario, ]
    gmt <- stats::setNames(g$median, g$year)
    horizon <- s[[1]]$model$horizon %||% c(1970L, 2023L)
    reference <- mean(gmt[as.character(seq(horizon[[1]], horizon[[2]]))])
    weight <- s[[1]]$model$weight %||% "population"
    p <- SCEN$pattern[[weight]]
    do.call(rbind, lapply(c(2050, 2100), function(y) {
      warming <- gmt[[as.character(y)]] - reference
      data.frame(iso3 = p$iso3, year = y, warming = warming,
                 dT = p$t_slope * warming, dP = 100 * p$p_rel_slope * warming)
    })) |>
      transform(label = FORCING_LABELS[[input$map_scenario]], weight = weight,
                reference = sprintf("%d-%d", horizon[[1]], horizon[[2]]))
  })
  driver_map <- function(variable, year) {
    d <- driver_data()
    values <- d[d$year == year, ]
    title <- sprintf("%d: %s change, %s (global +%.2f C)", year,
                     if (variable == "dT") "temperature" else "precipitation",
                     values$label[[1]], values$warming[[1]])
    subtitle <- sprintf(paste("Against the mean of %s; CMIP6 pattern (median over the models, %s",
                              "weights) times the FaIR median warming"),
                        values$reference[[1]], values$weight[[1]])
    base <- ggplot(world_with(values$iso3, values[[variable]])) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      labs(title = title, subtitle = wrap_label(subtitle, 90)) +
      theme_map()
    if (variable == "dT") {
      base + scale_fill_gradientn(colours = SCEN$palettes$temp_seq,
                                  limits = c(min(0, d$dT, na.rm = TRUE), max(d$dT, na.rm = TRUE)),
                                  na.value = MAP_NA, name = "Temperature\nchange (C)")
    } else {
      # A few arid countries have very large relative changes: the colours stop
      # at the 95th percentile of the absolute changes (beyond: the end colours).
      limit <- stats::quantile(abs(d$dP), 0.95, na.rm = TRUE, names = FALSE)
      base + scale_fill_gradientn(colours = SCEN$palettes$prec_div, limits = c(-limit, limit),
                                  oob = scales::squish, na.value = MAP_NA,
                                  name = sprintf("Precipitation change (%%)\n(colours capped at +/-%.0f%%)",
                                                 limit))
    }
  }
  output$driver_T_2050 <- renderPlot(driver_map("dT", 2050), bg = SURFACE)
  output$driver_T_2100 <- renderPlot(driver_map("dT", 2100), bg = SURFACE)
  output$driver_P_2050 <- renderPlot(driver_map("dP", 2050), bg = SURFACE)
  output$driver_P_2100 <- renderPlot(driver_map("dP", 2100), bg = SURFACE)

  global_table <- function() {
    if (scen_runs() == 0) return(NULL)
    out <- do.call(rbind, lapply(projections(), function(p) {
      g <- p$result$global
      cell <- function(y) loss_text(g[g$year == y, ])
      followed <- if (p$result$follow) {
        sprintf("%d / %d", p$result$followed[["income"]], p$result$followed[["share"]])
      } else "fixed"
      data.frame(Setting = p$setting, Scenario = p$label,
                 `Forcing (warming of)` = sprintf("%s (%s)", FORCING_LABELS[[p$scenario]],
                                                  SCEN$scenarios[[p$scenario]]),
                 `Economy (income, population, share)` = p$result$ssp,
                 `World GDP 2050 (%)` = cell(2050), `World GDP 2100 (%)` = cell(2100),
                 Countries = p$result$n,
                 `Following the SSP (income / share)` = followed,
                 `Share of the SSP's world GDP covered` = sprintf("%.0f%%", 100 * p$result$coverage),
                 `Mean growth of world GDP 2026-2100 (% a year): without, with climate change` =
                   growth_text(g, run()$specs, p$n, FALSE),
                 check.names = FALSE)
    }))
    names(out) <- sub("^World GDP", paste("World", quantity_of(run()$specs)), names(out))
    out[order(match(out$Setting, names(run_colours())),
              match(out$`Forcing (warming of)`, sprintf("%s (%s)", FORCING_LABELS, SCEN$scenarios))), ]
  }
  output$global_table <- renderTable(global_table(), striped = TRUE, spacing = "s", align = "l")

  scenario_notes <- function() {
    tags$ul(
      class = "small",
      tags$li(paste0(
        "Global warming: ", SCEN$source, "; temperature above 1850-1900, median and 5-95% over the ",
        "ensemble. The projections use the median path.")),
      tags$li(paste(
        "Country climate: the forced signal, with or without natural variability. The forced anomaly",
        "is the global warming above its mean over the panel's horizon, times the country's CMIP6 pattern",
        "(median over the models of the best-fit linear relationship with global temperature; T:",
        "the slope, C per C; P: the slope relative to each model's own mean P, times the country's",
        "long-run P).")),
      tags$li(paste(
        "Natural variability (optional): each country's annual T and P of the panel's climate data",
        "over its horizon less the country's own quadratic trend over it, resampled by year and added",
        "to the forced anomaly from 1994 (the history the trailing baselines need) to 2100, and on its",
        "own to the counterfactual (no warming); the loss is the difference, so the noise enters only",
        "through the non-linearity of the response (the arms' kink, the squared arms, the trailing",
        "climate of the interactions). Years are drawn independently for each country, jointly (the",
        "same calendar year for every country: the spatial correlation of T and P kept), or jointly",
        "in blocks of 5 consecutive years (also their persistence). The noise paths are the same for",
        "every setting and forcing. The estimate is the mean over the paths at the estimated",
        "coefficients (the maps too); the band pairs each path with a draw of the coefficients, the",
        "dotted lines are the 2.5-97.5% range of the paths alone. The moving SD widens with the forced",
        "signal only (the noise is the variability the SD already measures). Only the countries with",
        "the whole historical climate enter.")),
      tags$li(paste(
        "Damages: each setting's model and options (baseline, arms, interactions, shape, cap, SD and",
        "long-run characteristics) applied to every country of its panel with the country's own",
        "long-run climate, SD, income and share; losses accumulate from 2025 (no deviation from then",
        "on in the counterfactual). The custom country of a setting does not enter.")),
      tags$li(paste0(
        "Scenarios combine a forcing, which gives the global warming, with an economic development ",
        "(SSP1-5), which gives income, population (the GDP weights) and the agricultural share. With ",
        "'Follow SSP-RCP' each forcing takes its own SSP: SSP1 (1.9, 2.6), SSP2 (4.5), SSP3 (7.0), ",
        "SSP4 (3.4, 6.0), SSP5 (8.5); otherwise the chosen SSP goes with every forcing (labels like ",
        "'4.5 W/m2 x SSP3'). The warming of a forcing is that of its own emissions scenario. ",
        SCEN$ssp_source, ". With the long-run characteristics following a trailing baseline, the ",
        "interactions use the trailing 20 (30) year means of the long-run climate, of income and of ",
        "the agricultural share (for a model estimated with trailing characteristics, the means of ",
        "income and, if estimated so, share over the 5 years before, also with the horizon baseline; ",
        "its fixed share stays at the sample mean). Income and share are the country's historical series in the model's ",
        "definitions (GDP per capita of the panel's GDP source; agriculture over the three sectors' ",
        "value added) continued with the SSP's relative change from their last year. SSP GDP per ",
        "capita and population (5-year points) and the SSP share (5-year points 2030-2050) are joined ",
        "by a monotone cubic spline of the log (smooth growth through every point, no overshoot); ",
        "after 2050 the share follows SSP income with the country's ",
        "2030-2050 elasticity. Fixed: the sample means, as in estimation. The table counts the ",
        "countries that follow the SSP in income and in share (the others keep their sample means).")),
      tags$li(paste(
        "Growth of world GDP (the chart's other option, and the table): the yearly change of the",
        "summed GDP of the countries covered, from 2026 (the SSP's GDP starts in 2025). Without",
        "climate change it is the SSP's; with it, the SSP's GDP times one plus the global loss",
        "(the GDP-weighted average of the countries' losses), and its band is that of the draws of",
        "the coefficients, draw by draw; the third option is their difference, the change of the",
        "growth rate by climate change. Not available for the capital outcome.")),
      tags$li(paste(
        "Ignore gains (option): each country's effect enters the global loss capped at zero, year by year,",
        "at every draw of the coefficients and every noise path (a gain counts as the baseline GDP); the",
        "countries keep their GDP weights, and the maps and country results show the effects as they are.",
        "Global loss: the average of the countries' losses weighted by the SSP's GDP of each year",
        "(PPP), over the countries with a pattern and SSP GDP; the share of the SSP's world GDP they",
        "cover in 2025 is in the table. 95% CI: 1,000 draws of the coefficients, the same for every",
        "country. Not included: the uncertainty of the warming and of the CMIP6 pattern."))
    )
  }
  output$scenario_notes <- renderUI(scenario_notes())

  # Geoengineering and termination shocks -------------------------------------
  # Global warming paths with geoengineering: none (the background),
  # sustained (limited from 2025), and terminated in each year (termination_path()).
  # They run only when asked, for the shown settings at the moment of the click.
  TS_NOT_RUN <- "Press \"Run the geoengineering scenarios\" to project them."
  # The catch-up times of the comparison (positive numbers of a list).
  ts_compare_taus <- function(text) {
    x <- suppressWarnings(as.numeric(trimws(strsplit(text %||% "", ",")[[1]])))
    unique(x[is.finite(x) & x > 0])
  }
  # The alternatives of a run: each a forcing, an economy and a catch-up time
  # (the options as set, the compared one taking each of its values), with its
  # label.
  # (and a geoengineering pattern, `sai`: "uniform", "mix" or a latitude; the
  # CIDER model and fractions are those of the run.)
  ts_variants_of <- function(compare, forcing, economy, tau, forcings, economies, taus, sai = "uniform",
                             sais = character(0), model = NULL, fractions = NULL) {
    one <- function(f, e, t, s = sai) list(forcing = f, economy = e, tau = t, sai = s)
    v <- switch(compare,
                forcing = lapply(forcings, function(f) one(f, economy, tau)),
                economy = lapply(economies, function(e) one(forcing, e, tau)),
                tau = lapply(taus, function(t) one(forcing, economy, t)),
                sai = lapply(sais, function(s) one(forcing, economy, tau, s)),
                list(one(forcing, economy, tau)))
    for (i in seq_along(v)) {
      v[[i]]$label <- switch(compare, tau = sprintf("%s, tau = %g years", scenario_label(v[[i]]$forcing, v[[i]]$economy), v[[i]]$tau),
                             scenario_label(v[[i]]$forcing, v[[i]]$economy))
      if (compare == "sai" || v[[i]]$sai != "uniform") {
        v[[i]]$label <- paste0(v[[i]]$label, "; ", sai_label(v[[i]]$sai, model, fractions))
      }
    }
    v
  }
  # The injection fractions as set (named by latitude).
  ts_fractions <- function() {
    stats::setNames(vapply(CIDER_LATS, function(l) as.numeric(input[[paste0("ts_f_", l)]] %||% NA_real_), 0), CIDER_LATS)
  }
  ts_inputs <- function() {
    compare <- input$ts_compare %||% "none"
    model <- input$ts_cider_model %||% if (!is.null(CIDER)) names(CIDER$models)[[1]] else ""
    fractions <- ts_fractions()
    list(compare = compare, model = model, fractions = fractions,
         variants = ts_variants_of(compare, input$ts_forcing, input$ts_econ %||% "follow", input$ts_tau,
                                   intersect(SCENARIO_CHOICES, input$ts_cmp_forcings), input$ts_cmp_economies,
                                   ts_compare_taus(input$ts_cmp_taus),
                                   sai = if (identical(input$ts_sai, "cider")) "mix" else "uniform",
                                   sais = intersect(TS_SAI_CHOICES, input$ts_cmp_sai), model = model,
                                   fractions = fractions))
  }
  # The background warming of a forcing (FaIR median, named by year).
  ts_background <- function(forcing) {
    g <- SCEN$gmt[SCEN$gmt$scenario == forcing, ]
    stats::setNames(g$median, g$year)
  }
  # The CIDER pattern of alternative v of run r: NULL (uniform) or list(model,
  # fractions).
  ts_sai_of <- function(r, v) {
    if (is.null(v$sai) || v$sai == "uniform") return(NULL)
    list(model = r$cider_model, fractions = sai_fractions(v$sai, r$fractions))
  }
  cider_key <- function(model, fractions) paste(model, paste(signif(fractions, 4), collapse = "/"))
  # The injection of alternative r (a run with one pattern) with CIDER's
  # pattern, linear in the steady-state cooling per Tg (no saturation of the
  # aerosol): the largest global cooling of the sustained path, its injection
  # and the latitudes beyond the rate of the training runs; "" if uniform.
  # (A new line, after the label of the alternative when there are several.)
  ts_injection_text <- function(r, many = FALSE) {
    sai <- ts_sai_of(r, r)
    if (is.null(sai)) return("")
    mix <- cider_mix(CIDER, sai$model, sai$fractions)
    paths <- ts_paths_of(r)
    years <- as.character(seq(PROJECTION_FROM, PROJECTION_TO))
    cooling <- paths[[1]][years] - paths[[2]][years]
    i <- which.max(cooling)
    total <- cooling[[i]] * mix$tg_per_c
    per <- total * mix$fractions
    over <- names(per)[per > CIDER_TRAINING_TG]
    paste0("\n", if (many) paste0(r$variant, ": "),
           sprintf("injection up to %.0f Tg SO2/yr (%s, %.2f C of cooling; linear in CIDER's steady-state cooling per Tg)%s",
                   total, years[[i]], cooling[[i]],
                   if (length(over)) sprintf(", beyond the %g Tg/yr of the training runs at %s", CIDER_TRAINING_TG,
                                             paste(over, collapse = ", ")) else ""))
  }
  # The limit of geoengineering as set: list(mode = "level" or "rate", value),
  # and its key.
  ts_target <- function() {
    mode <- if (identical(input$ts_mode, "rate")) "rate" else "level"
    list(mode = mode, value = as.numeric(if (mode == "rate") input$ts_rate else input$ts_held))
  }
  target_key <- function(mode, value) paste(mode, value)
  ts_run <- eventReactive(ts_runs(), {
    req(ts_runs() > 0)
    x <- ts_inputs()
    target <- ts_target()
    uses_cider <- any(vapply(x$variants, function(v) v$sai != "uniform", TRUE))
    validate(need(isTRUE(target$value >= 0), "The limit of geoengineering must be a number of at least 0."),
             need(length(x$variants) > 0, "Select at least one alternative to compare."),
             need(all(vapply(x$variants, function(v) isTRUE(v$tau > 0), TRUE)), "The catch-up time must be positive."),
             need(!uses_cider || !is.null(CIDER), "No CIDER patterns: run Rscript apps/response_explorer/build_cider.R."),
             need(!("mix" %in% vapply(x$variants, `[[`, "", "sai")) || isTRUE(sum(pmax(x$fractions, 0), na.rm = TRUE) > 0),
                  "The injection fractions must have a positive sum."))
    v1 <- x$variants[[1]]
    # (forcing, economy, tau, pattern: those of the first alternative.)
    list(specs = specs(), forcing = v1$forcing, mode = target$mode, held = target$value,
         years = sort(as.integer(input$ts_years)), tau = v1$tau, economy = v1$economy, sai = v1$sai,
         cider_model = x$model, fractions = x$fractions,
         compare = x$compare, variants = x$variants, noise = noise_of(input$ts_noise, input$ts_noise_n),
         no_gains = isTRUE(input$ts_no_gains))
  })
  ts_key <- function(sp, variants, held, years, noise, no_gains, cider) {
    paste(c(vapply(sp, function(s) paste(s$label, s$model$key, s$sd_mode, s$clim_mode, paste(s$cap, collapse = "/")), ""),
            vapply(variants, function(v) paste(v$forcing, v$economy, v$tau, v$sai), ""), held, years, noise_key(noise), no_gains,
            cider),
          collapse = " || ")
  }
  ts_status <- function() {
    if (ts_runs() == 0) return(helpText(class = "mb-0", "Not run yet."))
    r <- ts_run()
    x <- ts_inputs()
    target <- ts_target()
    current <- ts_key(specs(), x$variants, target_key(target$mode, target$value), sort(as.integer(input$ts_years)),
                      noise_of(input$ts_noise, input$ts_noise_n), isTRUE(input$ts_no_gains), cider_key(x$model, x$fractions))
    if (!identical(current, ts_key(r$specs, r$variants, target_key(r$mode, r$held), r$years, r$noise, r$no_gains,
                                   cider_key(r$cider_model, r$fractions)))) {
      return(helpText(class = "mb-0 text-warning-emphasis",
                      "Settings or options changed since the last run: press Run to update."))
    }
    helpText(class = "mb-0", sprintf("Up to date: %d setting(s), %s, %s, %d termination year(s); %s.",
                                     length(r$specs), paste(vapply(r$variants, `[[`, "", "label"), collapse = " | "),
                                     geo_target_text(r$mode, r$held), length(r$years),
                                     paste0(noise_label(r$noise), if (r$no_gains) "; gains ignored in the world aggregate")))
  }
  output$ts_status <- renderUI(ts_status())
  # The alternative of the map, latitude chart and exceedances: the run with
  # its forcing, economy and catch-up time.
  ts_one <- reactive({
    r <- ts_run()
    k <- match(input$ts_variant %||% "", vapply(r$variants, `[[`, "", "label"))
    ts_with(r, r$variants[[if (is.na(k)) 1L else k]])
  })
  # The paths of alternative r (a run with one forcing and catch-up time):
  # name -> global warming (named by year).
  ts_paths_of <- function(r) {
    g <- SCEN$gmt[SCEN$gmt$scenario == r$forcing, ]
    background <- stats::setNames(g$median, g$year)
    paths <- c(stats::setNames(list(background), "no geoengineering"),
               stats::setNames(list(termination_path(background, r$held, mode = r$mode)),
                               paste0("sustained, ", geo_target_text(r$mode, r$held))))
    for (y in r$years) {
      paths[[paste("terminated in", y)]] <- termination_path(background, r$held, y, r$tau, mode = r$mode)
    }
    paths
  }
  # The run with alternative v.
  ts_with <- function(r, v) {
    r$forcing <- v$forcing
    r$economy <- v$economy
    r$tau <- v$tau
    r$sai <- v$sai %||% "uniform"
    r$variant <- v$label
    r
  }
  ts_paths <- reactive(ts_paths_of(ts_one()))
  # (The alternatives in their order, for the facets.)
  ts_variant_factor <- function(x) factor(x, levels = vapply(ts_run()$variants, `[[`, "", "label"))
  # The colours of the paths: grey without geoengineering, blue sustained, an
  # orange ramp for the termination years (darker later).
  ts_colours <- reactive({
    p <- names(ts_paths())
    n_term <- length(p) - 2
    stats::setNames(c(INK_MUTED, SETTING_COLOURS[[1]],
                      grDevices::colorRampPalette(c("#f3a272", "#eb6834", "#8f3510"))(n_term)), p)
  })
  # Every setting and path of the last run: list of (setting, path, result).
  # Every alternative, setting and path of the last run: list of (variant,
  # setting, path, result).
  # (`region`: the countries aggregated, NULL: the world.)
  ts_projections_for <- function(region = NULL) {
    r <- ts_run()
    paths <- lapply(r$variants, function(v) ts_paths_of(ts_with(r, v)))
    jobs <- expand.grid(i = seq_along(r$specs), path = names(paths[[1]]), k = seq_along(r$variants),
                        stringsAsFactors = FALSE)
    withProgress(message = if (is.null(region)) "Projecting the geoengineering scenarios" else "Projecting the selected countries",
                 value = 0, {
      lapply(seq_len(nrow(jobs)), function(j) {
        s <- r$specs[[jobs$i[j]]]
        v <- r$variants[[jobs$k[j]]]
        ssp <- economy_of(v$forcing, v$economy)
        path <- jobs$path[j]
        incProgress(1 / nrow(jobs), detail = paste(s$label, "-", v$label, "-", path))
        # (Without geoengineering the pattern does not matter.)
        sai <- if (path != "no geoengineering") ts_sai_of(r, v)
        key <- paste(s$model$key, s$sd_mode, s$clim_mode, paste(s$cap, collapse = "/"), v$forcing, ssp,
                     target_key(r$mode, r$held), v$tau, path,
                     noise_key(r$noise), if (r$no_gains) "no gains", if (!is.null(region)) paste(region, collapse = "+"),
                     if (!is.null(sai)) paste("CIDER", cider_key(sai$model, sai$fractions)),
                     sep = " | ")
        if (is.null(projection_cache[[key]])) {
          projection_cache[[key]] <- project_with(s, v$forcing, ssp, r$noise, gmt = paths[[jobs$k[j]]][[path]],
                                                  no_gains = r$no_gains, region = region, sai = sai)
        }
        list(variant = v$label, setting = s$label, n = s$n, path = path, result = projection_cache[[key]])
      })
    })
  }
  ts_projections_all <- reactive(ts_projections_for())
  # Those of the chosen alternative (map, table).
  ts_projections <- reactive({
    Filter(function(p) p$variant == ts_one()$variant, ts_projections_all())
  })
  # The global chart's: the world's, or the selected countries' aggregate.
  ts_region <- reactive(region_of("ts_region")) |> debounce(600)
  ts_global_projections <- reactive({
    region <- ts_region()
    if (is.null(region)) ts_projections_all() else ts_projections_for(region)
  })
  ts_global <- reactive({
    do.call(rbind, lapply(ts_global_projections(), function(p) {
      cbind(variant = p$variant, setting = p$setting, path = p$path, p$result$global)
    }))
  })
  observeEvent(ts_run(), {
    r <- ts_run()
    choices <- stats::setNames(as.character(vapply(r$specs, `[[`, 0, "n")),
                               vapply(r$specs, `[[`, "", "label"))
    updateSelectInput(session, "ts_map_setting", choices = choices, selected = choices[[1]])
    updateSelectInput(session, "ts_map_year", choices = as.character(r$years),
                      selected = if (length(r$years)) as.character(r$years[[1]]) else character(0))
    p <- names(ts_paths_of(r))
    updateSelectInput(session, "ts_lat_paths", choices = p, selected = utils::head(p[-1], 2))
    labels <- vapply(r$variants, `[[`, "", "label")
    updateSelectInput(session, "ts_variant", choices = labels, selected = labels[[1]])
  })
  ts_run_colours <- function(x) factor(x, levels = names(ts_colours()))

  # Exceedances of natural variability: the setting of the map (its panel's
  # climate data and pattern weight), the SSP population of the run, each
  # country's natural SD of T and the global one (natural_sigma()).
  ts_sigma_setting <- reactive({
    sp <- ts_run()$specs
    s <- Filter(function(s) as.character(s$n) == input$ts_map_setting, sp)
    if (length(s)) s[[1]] else sp[[1]]
  })
  ts_population_of <- function(r) {
    ssp <- economy_of(r$forcing, r$economy)
    years <- seq(PROJECTION_FROM, PROJECTION_TO)
    gdp <- SCEN$ssp$gdp[SCEN$ssp$gdp$ssp == ssp & SCEN$ssp$gdp$year %in% years, ]
    pc <- SCEN$ssp$gdppc[SCEN$ssp$gdppc$ssp == ssp & SCEN$ssp$gdppc$year %in% years, ]
    gdp <- tapply(gdp$gdp, list(gdp$iso3, gdp$year), sum)
    pc <- tapply(pc$gdppc, list(pc$iso3, pc$year), sum)
    iso <- intersect(rownames(gdp), rownames(pc))
    pop <- gdp[iso, as.character(years), drop = FALSE] / pc[iso, as.character(years), drop = FALSE]
    pop[stats::complete.cases(pop), , drop = FALSE]
  }
  ts_population <- reactive(ts_population_of(ts_one()))
  # The natural SDs of alternative r (its SSP's population weights).
  ts_sigma_of <- function(r) {
    m <- ts_sigma_setting()$model
    weight <- m$weight %||% "population"
    variability <- EXPLORER$variability[[paste(m$gdp, m$climate, weight, sep = " | ")]]
    pop <- ts_population_of(r)
    natural_sigma(variability, SCEN$pattern[[weight]], pop[, as.character(PROJECTION_FROM)])
  }
  ts_sigma <- reactive(ts_sigma_of(ts_one()))
  ts_reference <- reactive({
    ref <- input$ts_ref %||% "20"
    list(window = if (ref == "horizon") NA_real_ else as.numeric(ref),
         horizon = ts_sigma_setting()$model$horizon %||% c(1970L, 2023L),
         label = if (ref == "horizon") {
           h <- ts_sigma_setting()$model$horizon %||% c(1970L, 2023L)
           sprintf("the %d-%d mean", h[[1]], h[[2]])
         } else sprintf("the trailing %s-year mean", ref))
  })
  # For each terminated path, by year: warming, its reference plus k SD, and
  # whether the warming exceeds it.
  ts_exceedance_of <- function(paths, sg) {
    k <- input$ts_k_gmt
    validate(need(isTRUE(k > 0), "k must be positive."))
    ref <- ts_reference()
    years <- seq(PROJECTION_FROM, PROJECTION_TO)
    out <- do.call(rbind, lapply(grep("^terminated in", names(paths), value = TRUE), function(p) {
      d <- gmt_deviation(paths[[p]], years, ref$window, ref$horizon)
      data.frame(path = p, termination = as.integer(sub("^terminated in ", "", p)), year = years,
                 warming = d$warming, threshold = d$reference + k * sg$gmt,
                 above = d$deviation > k * sg$gmt)
    }))
    # (No termination: no rows.)
    out %||% data.frame(path = character(0), termination = integer(0), year = integer(0), warming = numeric(0),
                        threshold = numeric(0), above = logical(0))
  }
  ts_exceedance <- reactive(ts_exceedance_of(ts_paths(), ts_sigma()))
  # "2040: 12 years (2041-2052)", the years from the termination on.
  exceedance_text <- function(ex) {
    vapply(split(ex, ex$termination), function(e) {
      y <- e$year[e$above & e$year >= e$termination]
      if (!length(y)) return(sprintf("%d: never", e$termination[1]))
      sprintf("%d: %d years (%d-%d)", e$termination[1], length(y), min(y), max(y))
    }, "")
  }

  # One panel per alternative of the run (a forcing and catch-up time; the
  # economy changes only the population weights of the natural SD).
  ts_gmt_plot <- function() {
    validate(need(ts_runs() > 0, TS_NOT_RUN))
    r <- ts_run()
    rate <- identical(input$ts_gmt_measure, "rate")
    parts <- lapply(r$variants, function(v) {
      rv <- ts_with(r, v)
      paths <- ts_paths_of(rv)
      sg <- ts_sigma_of(rv)
      d <- do.call(rbind, lapply(names(paths), function(k) {
        data.frame(path = k, year = as.integer(names(paths[[k]])), warming = unname(paths[[k]]))
      }))
      if (rate) {
        d <- d[order(d$path, d$year), ]
        d$warming <- unlist(lapply(split(d, factor(d$path, levels = unique(d$path))), function(x) warming_rate(x$year, x$warming)))
      }
      d <- d[d$year >= 2000, ]
      if (rate) {
        # The fastest warming of each terminated path, from its termination on.
        peaks <- vapply(grep("^terminated in", names(paths), value = TRUE), function(p) {
          x <- d[d$path == p & d$year >= as.integer(sub("^terminated in ", "", p)), ]
          i <- which.max(x$warming)
          sprintf("%s: %.2f C per decade in %d", sub("^terminated in ", "", p), x$warming[i], x$year[i])
        }, "")
        return(list(d = cbind(variant = v$label, d), ex = NULL, sd = NA_real_,
                    text = paste0(if (length(r$variants) > 1) paste0(v$label, ": "), paste(peaks, collapse = "; "),
                                  ts_injection_text(rv, length(r$variants) > 1))))
      }
      # The bands: for each terminated path, a vertical band from the first to
      # the last year its warming is above the reference climate plus k SD,
      # from the termination on.
      ex0 <- ts_exceedance_of(paths, sg)
      ex <- ex0[ex0$above & ex0$year >= ex0$termination, ]
      ex <- if (nrow(ex)) do.call(rbind, lapply(split(ex, ex$path), function(e) {
        data.frame(path = e$path[1], start = min(e$year) - 0.5, end = max(e$year) + 0.5)
      })) else ex[0, ]
      list(d = cbind(variant = v$label, d), ex = if (nrow(ex)) cbind(variant = v$label, ex),
           text = paste0(if (length(r$variants) > 1) paste0(v$label, ": "), paste(exceedance_text(ex0), collapse = "; "),
                         ts_injection_text(rv, length(r$variants) > 1)),
           sd = sg$gmt)
    })
    d <- do.call(rbind, lapply(parts, `[[`, "d"))
    ex <- do.call(rbind, lapply(parts, `[[`, "ex"))
    d$path <- ts_run_colours(d$path)
    d$variant <- ts_variant_factor(d$variant)
    bands <- if (!is.null(ex)) {
      ex$path <- ts_run_colours(ex$path)
      ex$variant <- ts_variant_factor(ex$variant)
      geom_rect(data = ex, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = path),
                inherit.aes = FALSE, alpha = 0.15, colour = NA)
    }
    sds <- unique(round(vapply(parts, `[[`, 0, "sd"), 2))
    many <- length(r$variants) > 1
    caption <- if (rate) {
      paste0(wrap_label(paste0(first_upper(RATE_TEXT), ". Fastest warming after each termination:"), if (many) 140 else 70),
             "\n", paste(vapply(parts, `[[`, "", "text"), collapse = "\n"))
    } else {
      paste0(wrap_label(sprintf(
        "Bands: from the first to the last year after the termination with warming above %s + %g SD of natural variability (1 SD = %s C of global warming). Years above it:",
        ts_reference()$label, input$ts_k_gmt, paste(sprintf("%.2f", sds), collapse = " / ")), if (many) 140 else 70), "\n",
        paste(vapply(parts, `[[`, "", "text"), collapse = "\n"))
    }
    ggplot(d, aes(year, warming, colour = path)) +
      geom_vline(xintercept = PROJECTION_FROM, colour = GRID, linewidth = 0.5) +
      (if (rate) geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4)) +
      bands +
      geom_line(aes(linetype = path == "no geoengineering"), linewidth = 0.9) +
      (if (many) facet_wrap(~variant, nrow = 1, labeller = labeller(variant = function(x) wrap_label(x, 40)))) +
      scale_colour_manual(values = ts_colours(), name = NULL,
                          guide = guide_legend(ncol = if (many) 4 else 2)) +
      scale_fill_manual(values = ts_colours(), guide = "none") +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
      labs(x = NULL, y = if (rate) RATE_AXIS else "Global warming above 1850-1900 (C)",
           subtitle = if (many) "Background: FaIR median of each forcing"
                      else sprintf("Background: %s (FaIR median); catch-up time %.1f years",
                                   FORCING_LABELS[[r$forcing]], r$tau),
           caption = caption) +
      theme_explorer() +
      theme(plot.caption = element_text(colour = INK_MUTED, size = 9, hjust = 0),
            strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0))
  }
  output$ts_gmt_plot <- renderPlot(ts_gmt_plot(), bg = SURFACE)
  # The heights of the charts: the global chart grows with the settings when
  # the alternatives are in columns.
  ts_many <- reactive(ts_runs() > 0 && length(ts_run()$variants) > 1)
  ts_global_height <- function() if (ts_many()) 160 + 230 * length(ts_run()$specs) else 460
  output$ts_charts <- renderUI({
    if (ts_many()) {
      tagList(plotOutput("ts_gmt_plot", height = 420),
              plotOutput("ts_global_plot", height = ts_global_height(), hover = hoverOpts("hover_ts", delay = 80)))
    } else {
      layout_columns(col_widths = c(4, 8), plotOutput("ts_gmt_plot", height = 460),
                     plotOutput("ts_global_plot", height = 460, hover = hoverOpts("hover_ts", delay = 80)))
    }
  })

  ts_global_plot <- function() {
    validate(need(ts_runs() > 0, TS_NOT_RUN))
    sp <- ts_run()$specs
    measure <- input$ts_measure %||% "loss"
    validate(need(measure == "loss" || growth_allowed(sp), GROWTH_NOTE))
    d <- as_measure(ts_global(), measure)
    d$setting <- factor(d$setting, levels = vapply(sp, `[[`, "", "label"))
    d$variant <- ts_variant_factor(d$variant)
    d$path <- ts_run_colours(d$path)
    many <- length(ts_run()$variants) > 1
    # Growth: the SSP's growth without climate change (the same for every
    # path), thin grey.
    baseline <- if (measure == "growth") {
      geom_line(data = d[d$path == levels(d$path)[1], ], mapping = aes(year, baseline_growth, group = interaction(setting, variant)), inherit.aes = FALSE,
                colour = INK_MUTED, linetype = "dotted", linewidth = 0.6)
    }
    ggplot(d, aes(year, estimate, colour = path, fill = path)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.06, colour = NA) +
      noise_lines(d) +
      baseline +
      (if (measure == "growth" && isTRUE(input$ts_history)) history_growth_layers(region_history(ts_region()))) +
      geom_line(aes(linetype = path == "no geoengineering"), linewidth = 0.9) +
      (if (many) {
        # Side by side: the alternatives in columns, the settings in rows (one y
        # scale per row).
        facet_grid(setting ~ variant, scales = y_scales(input$ts_common_y),
                   labeller = labeller(setting = function(x) wrap_label(x, 30), variant = function(x) wrap_label(x, 40)))
      } else {
        facet_wrap(~setting, ncol = min(3, nlevels(d$setting)), scales = y_scales(input$ts_common_y),
                   labeller = labeller(setting = function(x) wrap_label(x, 50)))
      }) +
      scale_colour_manual(values = ts_colours(), name = NULL, guide = guide_legend(ncol = 3)) +
      scale_fill_manual(values = ts_colours(), guide = "none") +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
      labs(x = NULL, y = measure_axis(measure, sp, ts_region()),
           subtitle = wrap_label(sprintf("%s; %s%s%s%s", if (many) first_upper(geo_target_text(ts_run()$mode, ts_run()$held)) else ts_run()$variants[[1]]$label,
                                         if (ts_run()$noise$type != "none") "mean over the noise paths; " else "",
                                         band_label(ts_run()$noise),
                                         region_note(ts_region(), lapply(ts_global_projections(), `[[`, "result")),
                                         if (measure == "growth") paste0("; grey dotted: the SSP's growth without climate change",
                                                                         if (isTRUE(input$ts_history)) history_growth_note(region_history(ts_region()), region_name(ts_region()))) else ""),
                                 110)) +
      theme_explorer() +
      theme(strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0),
            strip.text.y = element_text(colour = INK_MUTED, size = 9, angle = 0, hjust = 0))
  }
  output$ts_global_plot <- renderPlot(ts_global_plot(), bg = SURFACE)

  output$readout_ts <- renderUI({
    if (ts_runs() == 0) return(NULL)
    hover <- input$hover_ts
    if (is.null(hover)) return(helpText("Hover over the global chart to read the values."))
    measure <- input$ts_measure %||% "loss"
    if (measure != "loss" && !growth_allowed(ts_run()$specs)) return(NULL)
    d <- as_measure(ts_global(), measure)
    year <- max(min(d$year), min(round(hover$x), PROJECTION_TO))
    rows <- d[d$year == year, ]
    rows <- rows[order(match(rows$variant, vapply(ts_run()$variants, `[[`, "", "label")), match(rows$path, names(ts_colours()))), ]
    many <- length(ts_run()$variants) > 1
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(paste0(if (many) paste0(rows$variant, ", "), rows$setting, ", ", rows$path, ": ", loss_text(rows, measure))))
  })

  # Map: the extra loss from terminating in the chosen year against sustained
  # geoengineering, by country in 2100 (pp of GDP; negative: a larger loss);
  # without termination years, sustained geoengineering against none.
  ts_map <- function() {
    validate(need(ts_runs() > 0, TS_NOT_RUN))
    req(input$ts_map_setting)
    pr <- Filter(function(p) as.character(p$n) == input$ts_map_setting, ts_projections())
    sustained <- Filter(function(p) startsWith(p$path, "sustained"), pr)[[1]]$result$countries
    target <- geo_target_text(ts_run()$mode, ts_run()$held)
    if (length(ts_run()$years)) {
      req(input$ts_map_year)
      terminated <- Filter(function(p) p$path == paste("terminated in", input$ts_map_year), pr)
      req(length(terminated) == 1)
      a <- terminated[[1]]$result$countries
      b <- sustained
      legend <- "terminated minus\nsustained (pp)"
      title <- sprintf("%d: termination in %s against sustained geoengineering (%s)", PROJECTION_TO,
                       input$ts_map_year, target)
    } else {
      a <- sustained
      b <- Filter(function(p) p$path == "no geoengineering", pr)[[1]]$result$countries
      legend <- "sustained minus no\ngeoengineering (pp)"
      title <- sprintf("%d: sustained geoengineering (%s) against none", PROJECTION_TO, target)
    }
    a <- a[a$year == PROJECTION_TO, ]
    b <- b[b$year == PROJECTION_TO, ]
    difference <- a$loss - b$loss[match(a$iso3, b$iso3)]
    limit <- max(abs(difference), na.rm = TRUE)
    ggplot(world_with(a$iso3, difference)) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      scale_fill_gradient2(low = MAP_LOSS, mid = MAP_ZERO, high = MAP_GAIN, midpoint = 0,
                           limits = c(-limit, limit), na.value = MAP_NA,
                           name = paste0(first_upper(quantity_of(ts_run()$specs)), ", ", legend)) +
      labs(title = title,
           subtitle = wrap_label(paste0(pr[[1]]$setting, "; ", ts_one()$variant))) +
      theme_map()
  }
  output$ts_map <- renderPlot(ts_map(), bg = SURFACE)

  # Map of the regional pattern of the chosen alternative: the climate left
  # when geoengineering offsets 1 C of global warming (global temperature
  # unchanged), each country's CMIP6 change per C of warming plus CIDER's per C
  # of SAI cooling: T t_slope - t_sai (C), P 100 (p_rel_slope - p_sai) (%),
  # under the climate weights of the map's setting. Zero with the uniform
  # pattern.
  ts_sai_map <- function() {
    validate(need(ts_runs() > 0, TS_NOT_RUN))
    r <- ts_one()
    sai <- ts_sai_of(r, r)
    validate(need(!is.null(sai), paste("Uniform pattern: geoengineering that offsets global warming leaves every",
                                       "country's climate unchanged. Choose a CIDER alternative to see its pattern.")))
    weight <- ts_sigma_setting()$model$weight %||% "population"
    mix <- cider_mix(CIDER, sai$model, sai$fractions, weight)
    d <- merge(SCEN$pattern[[weight]], mix$pattern, by = "iso3")
    is_T <- !identical(input$ts_sai_var, "P")
    value <- if (is_T) d$t_slope - d$t_sai else 100 * (d$p_rel_slope - d$p_sai)
    limit <- max(abs(value), na.rm = TRUE)
    shares <- mix$shares[mix$shares > 0]
    ggplot(world_with(d$iso3, value)) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      (if (is_T) {
        scale_fill_gradient2(low = "#2166ac", mid = MAP_ZERO, high = "#b2182b", midpoint = 0,
                             limits = c(-limit, limit), na.value = MAP_NA, name = "T change (C)")
      } else {
        scale_fill_gradientn(colours = SCEN$palettes$prec_div, limits = c(-limit, limit), na.value = MAP_NA,
                             name = "P change (%)")
      }) +
      labs(title = sprintf("%s left by geoengineering that offsets 1 C of global warming",
                           if (is_T) "Temperature change" else "Precipitation change"),
           subtitle = wrap_label(sprintf(
             "%s. CMIP6 change per C of global warming plus %s's per C of SAI cooling; shares of the cooling: %s; %s weights",
             r$variant, CIDER$models[[sai$model]],
             paste(sprintf("%s %.0f%%", names(shares), 100 * shares), collapse = ", "), weight), 140)) +
      theme_map()
  }
  if (!is.null(CIDER)) output$ts_sai_map <- renderPlot(ts_sai_map(), bg = SURFACE)

  # What the injection fractions as set give (population weights).
  output$ts_sai_info <- renderUI({
    req(CIDER)
    model <- input$ts_cider_model %||% names(CIDER$models)[[1]]
    f <- ts_fractions()
    if (!isTRUE(sum(pmax(f, 0), na.rm = TRUE) > 0)) {
      return(helpText(class = "text-warning-emphasis", "The injection fractions must have a positive sum."))
    }
    mix <- cider_mix(CIDER, model, f)
    e <- CIDER$efficiency[CIDER$efficiency$model == model, ]
    shares <- mix$shares[mix$shares > 0]
    helpText(class = "small", sprintf(
      paste("%s, steady-state global cooling per Tg SO2/yr: %s. With these fractions: shares of the cooling %s;",
            "1 C of cooling takes about %.0f Tg SO2/yr (linear: no saturation of the aerosol; CIDER is trained",
            "on %g Tg/yr at one latitude)."),
      CIDER$models[[model]], paste(sprintf("%s %.3f C", e$lat, e$cooling_per_tg), collapse = ", "),
      paste(sprintf("%s %.0f%%", names(shares), 100 * shares), collapse = ", "), mix$tg_per_c, CIDER_TRAINING_TG))
  })

  # Latitude chart: the share of the population exposed to years beyond k SD of
  # natural variability, by 5-degree latitude band and year, for each chosen path.
  ts_lat_plot <- function() {
    validate(need(ts_runs() > 0, TS_NOT_RUN))
    validate(need(length(input$ts_lat_paths) > 0, "Select at least one warming path."),
             need(isTRUE(input$ts_k_lat > 0), "k must be positive."))
    paths <- ts_paths()
    chosen <- intersect(names(paths), input$ts_lat_paths)
    req(length(chosen) > 0)
    sg <- ts_sigma()
    ref <- ts_reference()
    sign <- input$ts_sign %||% "warm"
    years <- seq(PROJECTION_FROM, PROJECTION_TO)
    # With CIDER's pattern, geoengineering's own T pattern (the map setting's weights).
    r <- ts_one()
    sai <- ts_sai_of(r, r)
    t_sai <- if (!is.null(sai)) {
      mix <- cider_mix(CIDER, sai$model, sai$fractions, ts_sigma_setting()$model$weight %||% "population")
      stats::setNames(mix$pattern$t_sai, mix$pattern$iso3)
    }
    background <- if (!is.null(sai)) ts_background(r$forcing)
    d <- do.call(rbind, lapply(chosen, function(p) {
      cbind(path = p, exposure_by_latitude(paths[[p]], years, ref$window, ref$horizon, sg,
                                           ts_population(), LATITUDE_SHARES, input$ts_k_lat, sign,
                                           background = background, t_sai = t_sai))
    }))
    d$path <- factor(d$path, levels = chosen)
    baseline <- 100 * if (sign == "both") 2 * stats::pnorm(-input$ts_k_lat) else stats::pnorm(-input$ts_k_lat)
    events <- if (sign == "both") sprintf("beyond +-%g SD", input$ts_k_lat) else sprintf("above +%g SD", input$ts_k_lat)
    ggplot(d, aes(year, lat, fill = exposed)) +
      geom_tile(width = 1, height = 5) +
      facet_wrap(~path, ncol = min(3, length(chosen))) +
      scale_fill_gradientn(colours = SCEN$palettes$temp_seq, limits = c(0, 100),
                           name = "Population\nexposed (%)") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(-60, 90, by = 30), expand = c(0, 0)) +
      labs(x = NULL, y = "Latitude (degrees)",
           title = sprintf("Share of the population with a year of T %s of natural variability", events),
           subtitle = wrap_label(sprintf(
             "%s; deviation from %s; each country's normal with its SD of detrended historical T (%s), population (%s) spread over its area; %.1f%% with no forced deviation",
             ts_one()$variant, ref$label,
             ts_sigma_setting()$label, economy_of(ts_one()$forcing, ts_one()$economy), baseline), 140)) +
      theme_explorer() +
      theme(strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0), legend.position = "right")
  }
  output$ts_lat_plot <- renderPlot(ts_lat_plot(), bg = SURFACE)

  ts_table <- function() {
    if (ts_runs() == 0) return(NULL)
    pr <- ts_projections_all()
    many <- length(ts_run()$variants) > 1
    cell <- function(g, y) loss_text(g[g$year == y, ])
    out <- do.call(rbind, lapply(pr, function(p) {
      g <- p$result$global
      sustained <- Filter(function(q) q$variant == p$variant && q$setting == p$setting && startsWith(q$path, "sustained"), pr)[[1]]
      s100 <- sustained$result$global$estimate[sustained$result$global$year == PROJECTION_TO]
      data.frame(Alternative = p$variant, Setting = p$setting, Path = p$path,
                 `World GDP 2050 (%)` = cell(g, 2050), `World GDP 2100 (%)` = cell(g, PROJECTION_TO),
                 `2100 against sustained (pp)` = sprintf("%+.2f", g$estimate[g$year == PROJECTION_TO] - s100),
                 `Growth of world GDP 2026-2100 (% a year): mean; largest change by climate change (pp, year)` =
                   growth_text(g, ts_run()$specs, p$n, TRUE),
                 check.names = FALSE)
    }))
    names(out) <- sub("^World GDP", paste("World", quantity_of(ts_run()$specs)), names(out))
    out <- out[order(match(out$Alternative, vapply(ts_run()$variants, `[[`, "", "label")),
                     match(out$Setting, vapply(ts_run()$specs, `[[`, "", "label")),
                     match(out$Path, names(ts_colours()))), ]
    if (!many) out$Alternative <- NULL
    out
  }
  output$ts_table <- renderTable(ts_table(), striped = TRUE, spacing = "s", align = "l")

  ts_notes <- function() {
    tags$ul(
      class = "small",
      tags$li(paste(
        "Perfect geoengineering: a climate with the same global temperature is the same with or",
        "without geoengineering, so each country's climate is its CMIP6 pattern times global warming,",
        "as in the climate scenarios. Background: the FaIR median warming of the chosen forcing.")),
      if (!is.null(CIDER)) tags$li(paste(
        "CIDER pattern (stratospheric aerosol injection by latitude): the global warming paths are the",
        "same, but the geoengineering's change of global temperature, G(t) - B(t) (zero before 2025),",
        "has its own regional pattern: each country's T is t_slope (B - ref) + t_sai (G - B) and its P",
        "mean_P (p_rel_slope (B - ref) + p_sai (G - B)). The cooling is split across the injection",
        "latitudes by fraction times each latitude's steady-state cooling per Tg SO2/yr, and t_sai and",
        "p_sai mix the latitudes' patterns with those shares (P: relative to the model's own climatology,",
        "times its global P change per C of cooling). Patterns and efficiencies: the last 20 years of",
        "CIDER's 12 Tg/yr single-latitude runs of the chosen model, averaged over each country with the",
        "setting's climate weights (population 2000 or area, as its climate data).",
        CIDER$source, "No injection -> AOD -> temperature dynamics: termination follows the catch-up",
        "time tau as with the uniform pattern; the injection shown is linear (no aerosol saturation).",
        "Caveats: trained over an SSP2-4.5 background; precipitation skill weak relative to",
        "interannual variability; no AMOC, seasonal cycle or injection longitude; 60N and 60S are",
        "spring injections.")),
      tags$li(paste(
        "Sustained, from 2025: a level limit holds global warming at H C, S(t) = min(B(t), H) (a level",
        "below today's warming of about 1.4 C is an immediate cooling, itself a deviation for the damage",
        "model; 0 is the preindustrial climate); a rate limit lets it rise by at most r C per decade,",
        "S(t) = min(B(t), S(t-1) + r/10) (0 stops warming at the 2024 level). Terminated in year T:",
        "sustained until T, then the gap to the background closes exponentially,",
        "G(t) = B(t) - max(B(T) - S(T), 0) exp(-(t - T) / tau). With no termination year ticked, only",
        "the background and sustained geoengineering run.")),
      tags$li(paste(
        "Damages: each setting's model and options, every country of its panel, losses accumulating",
        "from 2025 against no climate change, and the global loss weighted by the SSP's GDP, as in the",
        "climate scenarios, with the same options of natural variability (the same noise paths for",
        "every warming path). The map and the last column compare termination with sustained",
        "geoengineering in 2100 (negative: termination costs more); without termination, the map",
        "compares sustained geoengineering with none.")),
      tags$li(paste(
        "Exceedances of natural variability (climate only, no damage model): the deviation of global",
        "warming from the reference climate (the trailing mean of the years before, or the horizon",
        "mean), times each country's CMIP6 slope, against the SD of the country's detrended",
        "historical T in the climate data of the map's setting. Warming chart: one SD of global",
        "warming is the SD of the population-weighted (SSP, 2025) mean of the countries' detrended T",
        "over their weighted mean slope; the bands shade the warming above the reference plus k SD,",
        "and the subtitle counts those years from each termination on. Latitude chart: each",
        "country's probability of a year beyond k SD, a normal with its SD around its forced",
        "deviation, weighted by its SSP population spread uniformly over its area in each 5-degree",
        "band (Natural Earth borders; countries without borders at 1:110m left out)."))
    )
  }
  output$ts_notes <- renderUI(ts_notes())

  notes_ui <- function() {
    model_lines <- lapply(specs(), function(s) {
      m <- s$model
      # The splines against the linear interactions (same rows): the Wald test
      # of the non-linear terms and the verdicts.
      test <- if (!is.null(m$spline_test)) {
        sprintf(paste(" Splines against linear interactions: non-linear terms jointly p = %.3g;",
                      "AIC %+.1f, BIC %+.1f; verdict by AIC / BIC / within adj. R2: %s."),
                m$spline_test$nonlinear_p, m$spline_test$aic, m$spline_test$bic,
                m$spline_test$verdict)
      } else ""
      # The level terms of the trailing characteristics (controls: not in the
      # responses or projections).
      levels <- if (nrow(m$levels)) {
        sprintf(" Level terms (controls, not in the responses and projections): %s.",
                paste(sprintf("%s %+.4f (SE %.4f, p = %.3g)", m$levels$term, m$levels$estimate,
                              m$levels$se, m$levels$p), collapse = "; "))
      } else ""
      moderators <- if (m$moderators == "trailing") {
        paste0("trailing characteristics (climate: ",
               if (!trailing_window(m$window)) "horizon mean"
               else paste0("trailing ", sub("lag_", "", m$window), "-year mean"),
               "; income: 5 years before; share: ", if (m$share_trailing) "5 years before" else "sample mean", ")")
      } else "fixed characteristics (sample means)"
      effects <- m$effects %||% "both"
      effects_text <- if (effects == "year") "year effects only (no unit effects or trends)" else {
        paste(if (effects == "unit") "unit effects only (no year effects)" else "year and unit effects", "with", FE_SHORT[[m$fe]])
      }
      convergence <- if (!is.null(m$convergence_coef)) {
        paste0(sprintf(" Convergence control, %s%s: %+.4f (SE %.4f, p = %.3g)", CONVERGENCE_SHORT[[m$convergence]],
                       if (isTRUE(m$convergence_iv)) paste0(", instrumented with ", CONVERGENCE_IV_SHORT[[m$convergence]]) else "",
                       m$convergence_coef[["estimate"]], m$convergence_coef[["se"]], m$convergence_coef[["p"]]),
               if (is.finite(m$convergence_coef[["first_stage_f"]] %||% NA)) {
                 sprintf(", first-stage F %.1f", m$convergence_coef[["first_stage_f"]])
               }, ".")
      } else ""
      if (!is.null(m$quad_coef)) {
        q <- m$quad_coef
        convergence <- paste0(convergence, sprintf(" Quadratic climate controls: %s.", paste(sprintf(
          "%s %+.4f (SE %.4f, p = %.3g)", vapply(paste0("quad:", m$quad[ifelse(grepl("TM", q$term), "TM", "RR")], ":", q$term),
                                                  term_label, ""), q$estimate, q$se, q$p), collapse = "; ")))
      }
      tags$li(tags$span("●", style = paste0("color:", SETTING_COLOURS[s$n], ";")),
              sprintf(paste("%s: %s x %s, regression unit: %s, %s, %s, interactions: %s",
                            "(%s) on %s; %s rows, %s countries from %d; controls: %s; %s.%s%s%s"), s$label, m$gdp, m$climate,
                      baseline_label(m$baseline, model_style(m)), form_label(m$form, m$symmetric),
                      effects_text, interaction_label(m$inter),
                      if (length(unique(m$shape %||% "linear")) == 1) (m$shape %||% "linear")[[1]] else paste0("T ", m$shape[["TM"]], ", P ", m$shape[["RR"]]), moderators, format(m$nobs, big.mark = ","), m$countries,
                      EXPLORER$start_year,
                      if (length(m$controls)) paste(CONTROL_SHORT[m$controls], collapse = ", ") else "none",
                      outlier_text(m$outliers, m$gdp), test, levels, convergence))
    })
    tags$div(
      tags$ul(model_lines),
      tags$p(paste(
        "Every model: year and unit effects (or only one of them), with no unit trends, unit linear trends, or unit linear",
        "and quadratic trends (the notebook's); the deviations are demeaned in the same way before",
        "their arms are cut (each unit's mean, mean and linear trend, or mean and linear and",
        "quadratic trend, fitted from 1960), except the innovations of the NRK filter and the",
        "deviations of symmetric models and of models without unit effects, which enter as they are;",
        "optionally a convergence control (last year's or the 5 years' log GDP per capita, or its",
        "growth last year or two years before, not in the responses or projections; optionally instrumented",
        "by 2SLS with its value one year earlier); the ticked controls (by default all: conflict, minor",
        "and war, and the post-Soviet, with Albania, and post-Yugoslav collapse years); errors",
        "clustered by country. Growth outliers: by default those whose annotated cause does not",
        "involve weather among the lowest and highest 5% of the growth reversals of the GDP source's",
        "own ERA5 panel from 1970 (the change of growth from the year before, g(t) - g(t-1), flagged",
        "in year t; PWT x ERA5: the notebook's sample) are dropped; optionally those of the 1% tails",
        "only, those of the tails of growth itself (the rule until November 2026), or none. The",
        "causes are annotated in data/growth_outlier_causes*.csv. They are events of a country and",
        "year, dropped from every panel of that source; the weather-driven ones (droughts,",
        "storms, floods, harvests) always stay. Characteristics in the interactions: fixed",
        "(the sample means: long-run climate, log income and share, constant within a country) or",
        "trailing (the mean climate of the deviation's own baseline, the trailing 5, 10, 20 or 30 years, or",
        "the horizon mean with the horizon baseline; log of the mean GDP per capita of the 5 years",
        "before, from the whole GDP history of the source; with WB GDP optionally the mean",
        "agricultural share of the 5 years before). Trailing characteristics in an interaction also",
        "enter on their own (level terms): without them the interaction of a non-negative arm",
        "would stand in for the omitted level. Trailing models run on the rows with every trailing",
        "value they use. The climate is centred at the median of the countries' long-run climate,",
        "income and share at their means over the model's rows. Climate weights: population density (2000)",
        "or area (plain area average) in each country's climate, in the estimation and in the CMIP6",
        "downscaling of the scenarios (its population-weighted or area-weighted fits); the same rows",
        "and growth-outlier exclusions. Deviations in SD of the country's variability, or in absolute",
        "units (C and metres, demeaned in the same way; the SD then plays no part in the model and",
        "the option of a moving SD does not apply). T and P terms are estimated together, each",
        "with its own interactions; a characteristic neither variable is interacted with is not",
        "used. Splines: each interacted characteristic enters as a restricted cubic spline (knots at",
        "the 5, 35, 65 and 95th percentiles of the panel's countries, linear beyond the outer knots),",
        "with its non-linear terms set to zero at the median climate and mean income and share.",
        "Verdicts: 'indifferent' when the difference cannot be told apart from noise (the Wald test",
        "of the non-linear terms, 5%).")),
      tags$ul(
        tags$li(paste(
          "Response: the effect on that year's growth of a year whose temperature (precipitation) is",
          "the given anomaly away from the baseline, in SDs of the country's variability, at its",
          "long-run climate, income and share. Estimated on deviations demeaned as the fixed effects.",
          "Symmetric models: the arms are replaced by the unsplit deviation (linear arms), its square",
          "(quadratic arms) or both, with one coefficient (and one set of interactions) for both sides,",
          "on the deviations as they are: with linear terms the fixed effects absorb what demeaning",
          "would remove, up to the interactions with trailing characteristics.")),
        tags$li(paste(
          "NRK filter (Nath, Ramey and Klenow, Eq. 4; option (c) of shock_filters.R): the deviation",
          "is the innovation of a pooled regression of T (P) on its three lags and the lags times the",
          "country's mean over the horizon, with country effects and no year effects, fitted on the",
          "annual climate of the panel's countries from 1950; in SDs of the innovation in the country",
          "over the horizon (or in C and m). The climate of the interactions is the horizon mean, also",
          "with trailing characteristics. In the GDP paths and the scenarios the deviation is the",
          "innovation the anomaly adds, a_t - sum_j phi_j a_t-j with phi_j = g_j + h_j times the",
          "country's long-run climate (the lags before the shift, or the history before 2025, included):",
          "a one-off year is a deviation that year only, while a permanent shift leaves (1 - sum_j phi_j)",
          "of itself every year.")),
        tags$li(paste(
          "Lags (optional, L = 5 by default when ticked, up to 20): every deviation term, with its interactions and",
          "income scaling, also enters at t-1, ..., t-L with coefficients of its own, each lag interacted with the",
          "characteristics of its own year (and the time-varying ones with their lagged level terms). The lags of",
          "the first years come from the climate of the 10 years before the start (deviations demeaned from 1960, as",
          "in the notebook), and those after a gap in the rows from the climate of the missing years; a model uses",
          "the rows with every lagged value it needs (trailing income before 1970 is often missing; with L above 10",
          "the rows start L - 10 years after the start, e.g. 1980 with 20 lags). The response",
          "is the cumulative effect of one year's anomaly over that year and the L after; in the paths and",
          "scenarios each year's growth sums the effects of its own deviation and of the L before (none before the",
          "shift, or before 2025 in the scenarios: the counterfactual shares that history). The cap applies lag by",
          "lag, and the table tests the lagged terms jointly.")),
        tags$li(paste(
          "Last year's value (as Kalkuhl and Wenz) and this year's value (as Burke, Hsiang and Miguel): the",
          "deviation terms of T (P) interacted with the level of T (P) the year before or the same year, in C",
          "(m a year) less the median of the countries' long-run means, with their level terms (they vary",
          "every year). This year's value makes the response to a deviation x depend on the level it brings",
          "(the country's climate plus the anomaly), so the response bends even with linear arms. In the GDP",
          "paths and scenarios both levels are the country's climate plus the anomaly of that year (and the",
          "noise, with natural variability); they are climate terms, never capped.")),
        tags$li(paste(
          "Outcome, in log points: the growth of GDP per capita; with PWT also the growth of the capital",
          "stock per capita (rnna / pop) or of GDP per worker (rgdpna / emp, persons engaged), on the rows",
          "with them. Responses, paths and losses are of the outcome: with capital or GDP per worker the",
          "GDP paths and scenario losses are of capital per capita or GDP per worker (weighted by the",
          "SSP's GDP in the global loss), as the axis labels say.")),
        tags$li(paste(
          "Capital interaction (PWT): log of the capital stock per capita (rnna / pop), long-run (mean",
          "over the horizon) or trailing (mean of the 5 years before, with its level term), less its",
          "mean over the rows. Its cross-country correlation with log income is about 0.92 (0.46 within",
          "countries after the quadratic trends). In the GDP paths and scenarios it stays at the country's",
          "value (there is no SSP capital path); the scenarios use the countries with PWT capital.")),
        tags$li(paste(
          "Income scaling, for the temperature and the precipitation terms separately:",
          "standard, or every deviation term (each arm and its ticked interactions) times (Ybar / Y)^(1 - xi),",
          "with Y the country's income of the interactions and Ybar a constant, the geometric mean of that",
          "income over the rows, so that the loss in dollars per person is",
          "proportional to Y^xi beyond the ticked interactions (xi = 0: equal dollars; xi = 1: equal %, the",
          "standard terms; xi = 2: the % loss rises with income), or both the plain terms and their scaled",
          "copies (nests the two; the table tests each). Y is the income of the interactions: long-run, or",
          "trailing (then with its level term, and along the SSP in the scenarios). The cap applies to the",
          "plain terms and to the scaled terms separately.")),
        tags$li(paste(
          "GDP path: growth effects accumulate in log GDP. The shift is immediate (solid: its full",
          "size from year 0 on), gradual (dashed: rising linearly from zero in year 0 to its full",
          "size in year X, then constant), sudden and late (dotted: nothing until year X, the full",
          "shift from then on) or transitory (dot-dash: the shift in year 0 only; without lags its GDP effect",
          "is a one-off step, with lags it moves for L years as the lagged effects come in, showing the",
          "dynamics of the distributed lags; with a trailing baseline the year also enters the trailing mean).",
          "Horizon baseline: the shift is a",
          "deviation every year, so growth is lower (or higher) every year. Trailing baselines: the",
          "shift enters the trailing mean, so an immediate shift fades out over 20 (30) years and a",
          "gradual one stays a deviation of about the warming of the last 10 (15) years while it",
          "rises, then fades out once it stops. The SD is fixed at the",
          "country's, or widened by the shift (moving). The long-run climate of the interactions is",
          "fixed at the country's sample mean, or follows the trailing baseline (the country's mean",
          "plus the shift as it enters the trailing mean: the response adapts as if the country's",
          "climate were that of the baseline); in the GDP paths income and share stay at the",
          "country's (they follow the SSP only in the climate scenarios). With fixed characteristics",
          "the models are estimated with the interactions on the sample (horizon) means, so following",
          "the baseline is an assumption of the paths; with trailing characteristics it is how they",
          "were estimated (the default then), and in the scenarios income and share follow the SSP",
          "as their 5-year trailing means.")),
        tags$li(paste(
          "Cap: each arm's response is its own coefficient c (the response at the sample's mean",
          "income and share), plus the income and share terms, plus the climate term. With the cap,",
          "the non-climate part (c plus the income and share terms) is at most max(c, 0), whatever",
          "the income and share, however high or low: income and share can bring a damaging response",
          "to zero but not past it, and cannot raise a beneficial one. The climate term, at the",
          "country's historical climate and any warming beyond it, is added on top, uncapped. The",
          "capped response is not linear in the coefficients, so its 95% CI comes from 2,000 draws",
          "of the coefficients (percentiles), capped draw by draw.")),
        tags$li(paste(
          "Held fixed: income and share of the interactions (no feedback from the GDP losses to",
          "income), the long-run climate unless it follows the baseline, and the absence of lags",
          "(no recovery).")),
        tags$li(paste0(
          "Country SDs: horizon SD of the annual series (also of the level), or for the trailing baselines the ",
          "median of the trailing SD over the country's years, or for the NRK filter the SD of its innovation ",
          "over the horizon; first differences are those of the unit in these units (with the level: Kalkuhl ",
          "and Wenz, this year's value less last year's, over the horizon SD). Agricultural share: World Bank agriculture value ",
          "added over that of the three sectors (constant 2015 US$), mean over the sample years. ",
          "Built ", format(EXPLORER$built, "%Y-%m-%d %H:%M"), "."))
      )
    )
  }
  output$notes <- renderUI(notes_ui())

  # The saved page ---------------------------------------------------------------
  # The page as it is shown, as one self-contained HTML file: the settings, the
  # charts (PNG, embedded), the tables and the notes. The climate scenarios and
  # the termination shocks enter only if they were run (their last run, with its
  # status). The equations need MathJax from the web to render.
  output$save_page <- downloadHandler(
    filename = function() paste0("response_explorer_", format(Sys.time(), "%Y%m%d_%H%M"), ".html"),
    content = function(file) {
      withProgress(message = "Saving the page", value = 0, {
        htmltools::save_html(saved_page(), file)
      })
    }
  )
  saved_page <- function() {
    # A chart as an embedded PNG drawn `width` px wide (about the app's on a
    # wide window) and `height` px high, shown at the width of its column; a
    # chart that cannot be drawn (e.g. growth with the capital outcome) is the
    # reason why.
    chart <- function(make, height, width = 760) {
      incProgress(1 / 24)
      p <- tryCatch(make(), error = function(e) conditionMessage(e))
      if (is.character(p)) return(tags$p(class = "muted", p))
      png <- tempfile(fileext = ".png")
      on.exit(unlink(png))
      ggsave(png, p, width = width / 96, height = height / 96, dpi = 192, bg = SURFACE)
      tags$img(src = knitr::image_uri(png), style = "width: 100%;")
    }
    # Two columns, of equal width or as `columns` (CSS grid columns).
    pair <- function(a, b, columns = "1fr 1fr") {
      tags$div(class = "pair", style = paste0("grid-template-columns: ", columns, ";"), tags$div(a), tags$div(b))
    }
    table_of <- function(df, html = FALSE) {
      if (is.null(df) || !nrow(df)) return(NULL)
      cell <- function(x) if (html) HTML(as.character(x)) else as.character(x)
      tags$table(class = "tbl",
                 tags$thead(tags$tr(lapply(names(df), function(x) tags$th(cell(x))))),
                 tags$tbody(lapply(seq_len(nrow(df)), function(i) {
                   tags$tr(lapply(df[i, ], function(x) tags$td(cell(x))))
                 })))
    }
    card <- function(title, ...) tags$section(class = "card", tags$h2(title), ...)
    # The common shift and, for each shown setting, its label.
    sp <- specs()
    settings <- tagList(
      tags$ul(lapply(sp, function(s) {
        tags$li(tags$span("●", style = paste0("color:", SETTING_COLOURS[s$n], ";")),
                tags$strong(paste0(" Setting ", s$n, ": ")), s$label)
      })),
      # The countries of the single-country charts and the setting each uses.
      tags$ul(lapply(pairs(), function(p) {
        cn <- p$country
        tags$li(tags$span("●", style = paste0("color:", SETTING_COLOURS[p$colour_n], ";")),
                tags$strong(paste0(" ", p$label, ": ")),
                sprintf("%s; mean T %.1f °C (SD %.2f), mean P %.0f mm (SD %.0f), GDP per capita %s, agricultural share %.1f%%",
                        cn$name, cn$mean_T, cn$sd_T, 1000 * cn$mean_P, 1000 * cn$sd_P,
                        format(round(cn$income), big.mark = ","), 100 * cn$share))
      })))
    shift_text <- sprintf(paste("Permanent shift: %+g °C in temperature, %+g mm a year in precipitation; GDP paths",
                                "over %d years; year X = %d; shifts shown: %s."),
                          input$shock_T, input$shock_P, as.integer(input$horizon), as.integer(input$ramp_years),
                          paste(unname(scenarios()), collapse = ", "))
    measure_text <- function(measure) names(MEASURE_CHOICES)[MEASURE_CHOICES == (measure %||% "loss")]
    scenarios_part <- if (scen_runs() > 0) {
      component <- input$map_component %||% "loss"
      card("Climate scenarios: GDP losses under each setting's model, 2025-2100",
           tags$div(class = "muted", scenario_status()),
           tags$p(class = "muted", paste0("Global chart: ", measure_text(input$global_measure), ". Maps: ",
                                          COMPONENT_TITLES[[component]], ".")),
           pair(chart(gmt_plot, 380, 520), chart(global_plot, 380, 1040), "1fr 2fr"),
           pair(chart(function() country_map(2050), 340), chart(function() country_map(2100), 340)),
           tags$h3("Physical drivers: the change of each country's climate in the scenario"),
           pair(chart(function() driver_map("dT", 2050), 320), chart(function() driver_map("dT", 2100), 320)),
           pair(chart(function() driver_map("dP", 2050), 320), chart(function() driver_map("dP", 2100), 320)),
           table_of(global_table()), scenario_notes())
    }
    termination_part <- if (ts_runs() > 0) {
      card("Geoengineering and termination shocks, 2025-2100",
           tags$div(class = "muted", ts_status()),
           tags$p(class = "muted", paste0("Global chart: ", measure_text(input$ts_measure), ".")),
           if (ts_many()) tagList(chart(ts_gmt_plot, 420, 1180), chart(ts_global_plot, ts_global_height(), 1180))
           else pair(chart(ts_gmt_plot, 460, 520), chart(ts_global_plot, 460, 1040), "1fr 2fr"),
           chart(ts_map, 380, 1180),
           if (length(input$ts_lat_paths)) chart(ts_lat_plot, 420, 1180),
           table_of(ts_table()), ts_notes())
    }
    tagList(
      tags$head(
        tags$title("Weather response explorer"),
        tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.9/MathJax.js?config=TeX-AMS-MML_HTMLorMML"),
        tags$style(HTML(sprintf(paste(
          "body { background: %s; color: %s; font-family: system-ui, sans-serif; max-width: 1220px;",
          "margin: 0 auto; padding: 16px; font-size: 14px; }",
          "h1 { font-size: 22px; } h2 { font-size: 17px; margin-top: 0; } h3 { font-size: 15px; }",
          ".card { border: 1px solid %s; border-radius: 8px; padding: 14px; margin: 14px 0; background: #fff; }",
          ".pair { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }",
          "@media (max-width: 800px) { .pair { grid-template-columns: 1fr !important; } }",
          ".muted { color: %s; } .small { font-size: 12.5px; }",
          ".alert { background: #fff6e0; border: 1px solid #eda100; border-radius: 6px; padding: 8px 12px; }",
          ".tbl { border-collapse: collapse; font-size: 12.5px; display: block; overflow-x: auto; margin: 10px 0; }",
          ".tbl th, .tbl td { border-bottom: 1px solid %s; padding: 4px 8px; text-align: left; vertical-align: top; }",
          ".tbl tbody tr:nth-child(odd) { background: #f5f4f1; }"), SURFACE, INK, GRID, INK_MUTED, GRID)))
      ),
      tags$h1("Weather response explorer"),
      tags$p(class = "muted", paste("Saved", format(Sys.time(), "%Y-%m-%d %H:%M"), "· panels built",
                                    format(EXPLORER$built, "%Y-%m-%d %H:%M"))),
      settings, tags$p(shift_text),
      warnings_ui(),
      pair(card("Temperature: growth response to an annual anomaly",
                chart(function() response_plot(responses_T(), "TM"), 380)),
           card("Precipitation: growth response to an annual anomaly",
                chart(function() response_plot(responses_P(), "RR"), 380))),
      pair(card(gdp_header("TM"), chart(function() gdp_plot(paths_T()), 430)),
           card(gdp_header("RR"), chart(function() gdp_plot(paths_P()), 430))),
      card("Each setting's model: specification, coefficients, fit and joint tests",
           spec_math(), table_of(coef_table(), html = TRUE), coef_notes()),
      scenarios_part, termination_part,
      card("Models and assumptions", notes_ui())
    )
  }
}

app <- shinyApp(ui, server)
# source()d in an interactive session: launch it; loaded by runApp(): return it.
if (interactive() && !is.null(SOURCED_FILE)) runApp(app) else app
