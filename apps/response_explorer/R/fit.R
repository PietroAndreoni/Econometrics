# Fits a model of the explorer on demand from the estimation panels of
# build_models.R (`panels`: the content of panels.rds), and keeps it for the
# session. The model is that of "Income or the agricultural share" in
# weather_arms_robustness.Rmd: growth on the arms of the demeaned deviations of
# T and P, each arm interacted with the chosen moderators, with the controls,
# year and unit effects (and unit trends), errors clustered by country.
#
# Fixed effects and demeaning go together: the deviations were demeaned by the
# unit mean ("unit"), mean and linear trend ("linear") or mean and linear and
# quadratic trend ("quadratic", the notebook's) before their arms are cut, and
# the regression has the same unit slopes. Two exceptions, never demeaned: the
# innovations of the Nath-Ramey-Klenow filter (baseline "NRK filter", window
# "nrk"; build_models.R), and the deviations of a symmetric model, which enter
# unsplit (x, x^2 or both in place of the linear, quadratic or both arms: one
# coefficient for the two sides).
#
# Regression unit (the baseline): the level of T and P (over the country's
# horizon SD in SD units; demeaned, the same as the horizon deviation), the
# deviation from the horizon or trailing 20 / 30 year mean, or the filter's
# innovation; style: the unit itself or its first difference (x_t - x_t-1 of
# the unit, demeaned in the same way; with the level, Kalkuhl and Wenz), or both
# (each with arms and interactions of its own; the first-difference terms
# are named dev<v>_d_<arm>).
#
# Moderators: "fixed", the horizon mean of the climate, long-run income,
# long-run share and (PWT) long-run capital per capita (constant within a unit,
# so they have no level terms of their own); "trailing", the climate of the
# deviation's own baseline (its trailing 20 or 30 year mean, the horizon mean
# with the horizon baseline), the log of the mean income and (PWT) of the mean
# capital per capita of the panels$trailing_income_years years before, and
# with `share_trailing` (WB GDP only) the mean share of the
# panels$trailing_share_years years before. Time-varying moderators enter with
# their level terms too (for those in at least one interaction): without them
# an interaction of a non-negative arm stands in for the omitted level. The
# levels are controls: the app's responses and projections use the deviation
# terms only. The trailing models run on the rows with all the trailing values
# they use (the full window); a model with capital (or a PWT outcome other than
# GDP per capita) on the rows with it.
#
# Last year's ("lag") and this year's ("cur") value of the variable, T or P (C,
# m a year), less the median of the units' long-run climate: time-varying,
# with their level terms, whatever the other moderators.
#
# Low income ("low"): a dummy, 1 when the income of the interactions (long-run
# or trailing) is below a threshold in the GDP source's units (default: the
# World Bank's low-income line, low_income_thresholds()); with trailing
# moderators it varies and has its level term too.
#
# Effects (EFFECTS_CHOICES): year and unit (the notebook's), unit only, or year
# only; optionally a quadratic climate control of T and of P (QUAD_CHOICES);
# and optionally a convergence control (CONVERGENCE_CHOICES), last year's
# or the 5 years' log GDP per capita, whatever the interactions.
#
# Outcome, in log points: the growth of GDP per capita ("log", the notebook's);
# with PWT also the growth of the capital stock per capita ("capital", rnna /
# pop) and of GDP per worker ("worker", rgdpna / emp).

library(fixest)

ARM_SUFFIXES <- list("Linear arms" = c("_pos", "_neg"), "Quadratic arms" = c("_pos2", "_neg2"),
                     "Linear + quadratic arms" = c("_pos", "_neg", "_pos2", "_neg2"))
# The unsplit terms of a symmetric model, by the arms they replace.
SYMMETRIC_SUFFIXES <- list("Linear arms" = "_sym", "Quadratic arms" = "_sym2",
                           "Linear + quadratic arms" = c("_sym", "_sym2"))
arm_suffixes <- function(form, symmetric = FALSE) {
  if (isTRUE(symmetric)) SYMMETRIC_SUFFIXES[[form]] else ARM_SUFFIXES[[form]]
}
# Whether a baseline window is a trailing one (its mean moves with the years).
trailing_window <- function(window) grepl("^lag_", window)
# The SD of the level ("level": T or P itself) is the country's horizon SD.
baseline_sd_window <- function(window) if (window == "level") "all" else window
# The deviations of v in the rows d: the filter's innovations with the NRK
# baseline; otherwise demeaned as the fixed effects fe, or with fe "none" as
# they are (v less the baseline mean, none for the level, over its SD in SD
# units). With `diff`, their first differences. As deviation_values() of
# build_models.R.
deviation_values <- function(d, v, window, u, fe, diff = FALSE) {
  if (window == "nrk") return(d[[paste0("nrk_", v, u, if (diff) "_d")]])
  if (fe != "none") return(d[[paste0("dm_", v, "_", window, if (diff) "_d", u, "_", fe)]])
  if (diff) return(d[[paste0("fd_", v, "_", window, u)]])
  x <- if (window == "level") d[[v]] else d[[v]] - d[[paste0("mean_", v, "_", window)]]
  if (u == "") x / d[[paste0("sd_", v, "_", baseline_sd_window(window))]] else x
}
# The regression unit (the baseline of the deviations; values: the keys of
# panels$baselines) and the style of the regression: the unit itself, or its
# first difference (with the level: Kalkuhl and Wenz).
BASELINE_CHOICES <- c("Level" = "Level", "Deviation (horizon)" = "Horizon", "Deviation 20 years" = "20-year",
                      "Deviation 30 years" = "30-year", "Filter (Nath-Ramey-Klenow AR(3) innovation)" = "NRK filter")
STYLE_CHOICES <- c("Simple (the variable only)" = "simple", "First order difference" = "diff",
                   "Both (the variable and its first difference)" = "both")
FE_CHOICES <- c("No unit trends" = "unit", "Unit linear trends" = "linear",
                "Unit linear and quadratic trends" = "quadratic")
FE_SHORT <- c(unit = "no unit trends", linear = "unit linear trends",
              quadratic = "unit linear and quadratic trends")
# The effects of the regression: year and unit (the notebook's), unit only (no
# year effects), or year only (no unit effects: then no unit trends either,
# the deviations enter as they are, not demeaned, and the constant
# characteristics in an interaction enter with their levels).
EFFECTS_CHOICES <- c("Year and unit effects" = "both", "Unit effects only (no year effects)" = "unit",
                     "Year effects only (no unit effects or trends; deviations not demeaned)" = "year")
EFFECTS_SHORT <- c(both = "year and unit effects", unit = "unit effects only", year = "year effects only")
# The convergence control, on its own (not in the responses or projections):
# last year's log GDP per capita of the GDP source, the log of its mean over
# the 5 years before (the trailing income), or its growth last year (the lagged
# dependent variable) or two years before (from t-3 to t-2, which leaves out
# the year t-1 that shares the shocks of t); none by default.
CONVERGENCE_CHOICES <- c("None" = "none", "Log GDP per capita last year" = "lag1",
                         "Log mean GDP per capita of the 5 years before" = "mean5",
                         "Growth of GDP per capita last year" = "growth1",
                         "Growth of GDP per capita two years before (t-3 to t-2)" = "growth2")
CONVERGENCE_COLUMNS <- c(lag1 = "gdp_l1", mean5 = "inc_trailing", growth1 = "growth_l1", growth2 = "growth_l2")
CONVERGENCE_SHORT <- c(lag1 = "last year's log GDP per capita", mean5 = "log mean GDP per capita of the 5 years before",
                       growth1 = "last year's growth of GDP per capita",
                       growth2 = "growth of GDP per capita two years before")
# The quadratic climate control of T and of P, on its own (Z + Z^2, not in the
# responses or projections): this year's value, last year's, or the trailing
# 20 or 30 year mean (years t-20..t-1 or t-30..t-1, as the baselines), in C and
# m a year less the median of the units' long-run climate; none by default.
QUAD_CHOICES <- c("None" = "none", "This year's value" = "cur", "Last year's value" = "lag1",
                  "Trailing 20-year mean (to t-1)" = "mean20", "Trailing 30-year mean (to t-1)" = "mean30")
QUAD_SHORT <- c(cur = "this year's", lag1 = "last year's", mean20 = "trailing 20-year mean",
                mean30 = "trailing 30-year mean")
quad_column <- function(kind, v) {
  switch(kind, cur = v, lag1 = paste0(v, "_l1"), mean20 = paste0("mean_", v, "_lag_20"),
         mean30 = paste0("mean_", v, "_lag_30"))
}
# The demeaning of the deviations of a model (as the fixed effects; "none":
# as they are, symmetric or without unit effects).
deviation_fe <- function(fe, symmetric, effects) if (symmetric || effects == "year") "none" else fe
SPLINE_COLUMNS <- c("clim_TM", "clim_RR", "lag_TM", "lag_RR", "cur_TM", "cur_RR", "inc", "agr", "kap")
# The modifiers of a variable's own climate: its long-run climate ("clim"),
# last year's value ("lag", as Kalkuhl and Wenz) and this year's value ("cur",
# as Burke, Hsiang and Miguel); their columns end in the variable (clim_TM).
OWN_MODIFIERS <- c("clim", "lag", "cur")
modifier_column <- function(z, v) if (z %in% OWN_MODIFIERS) paste0(z, "_", v) else z
OUTCOME_CHOICES <- c("Growth of GDP per capita (Δ log)" = "log",
                     "Log GDP per capita (level)" = "lgdp",
                     "Growth of the capital stock per capita (Δ log, PWT)" = "capital",
                     "Growth of GDP per worker (Δ log, PWT)" = "worker")
OUTCOME_SHORT <- c(log = "GDP per capita", lgdp = "log GDP per capita (level)", capital = "capital per capita",
                   worker = "GDP per worker")
# The outcomes that are a level, not a growth rate: their effects on log GDP
# are not summed over the years in the paths and projections.
LEVEL_OUTCOMES <- "lgdp"
level_outcome <- function(model) isTRUE(model$outcome %in% LEVEL_OUTCOMES)
# The outcomes and the moderator only PWT has.
PWT_OUTCOMES <- c("capital", "worker")
PWT_MODIFIERS <- "kap"
# The World Bank's low-income line, GNI per capita (Atlas method) of the FY2017
# classification, from 2015 GNI: in 2015 US$, the units of WB GDP.
LOW_INCOME_LINE <- 1025
# The default low-income thresholds of the GDP sources: the line for WB (2015
# US$); for PWT (2021 PPP US$) the line times the median ratio of PWT to WB
# trailing income over the country-years under it (ERA5 panels, population
# weights). Cached in panels$low_income when set.
low_income_thresholds <- function(panels) {
  if (!is.null(panels$low_income)) return(panels$low_income)
  pick <- function(gdp) {
    r <- panels$panels[[paste(gdp, "ERA5 | population", sep = " | ")]]$rows
    r[is.finite(r$inc_trailing), c("GID_0", "year", "inc_trailing")]
  }
  m <- merge(pick("PWT"), pick("WB"), by = c("GID_0", "year"), suffixes = c("_pwt", "_wb"))
  under <- exp(m$inc_trailing_wb) <= LOW_INCOME_LINE
  c(WB = LOW_INCOME_LINE, PWT = LOW_INCOME_LINE * stats::median(exp(m$inc_trailing_pwt[under] - m$inc_trailing_wb[under])))
}
# How income scales the deviation terms of T and of P, each on its own:
# "standard", each arm with its interactions; "dollars", every
# deviation term (the arm and its interactions) times idol_v = (Ybar / Y_i)^(1 -
# xi_v) (Ybar: the geometric mean income of the rows; Y_i the income of the
# interactions, long-run or trailing), so that the dollar loss per person is
# proportional to Y_i^xi_v beyond the ticked interactions (xi = 0: the same
# dollars at every income; xi = 1: the same %, the standard terms); "both", the
# terms and their copies times idol_v (the encompassing model). In the
# coefficients the scaled terms end in ":idol" (idol_TM or idol_RR by arm).
INCOME_FORMS <- c("Standard (income through its interaction only)" = "standard",
                  "Both: every deviation term also times (Ybar / Y)^(1 - xi)" = "both",
                  "Scaled: every deviation term times (Ybar / Y)^(1 - xi)" = "dollars")
XI_RANGE <- c(0, 2)
# The groups of controls (panels$control_sets; all of them: the notebook's)
# and the growth outliers dropped: the non-weather ones of the 1% or 5% tails
# of the change of growth from the year before (reversals, r1 and r5; r5 is
# the notebook's) or of growth itself (levels, 1 and 5), or none.
CONTROL_CHOICES <- c("Conflict (minor, war)" = "conflict", "Post-Soviet collapse years" = "soviet",
                     "Post-Yugoslav collapse years" = "yugoslav")
CONTROL_SHORT <- c(conflict = "conflict", soviet = "Soviet collapse", yugoslav = "Yugoslav collapse")
OUTLIER_CHOICES <- c("Growth reversals, 5% tails" = "r5", "Growth reversals, 1% tails" = "r1",
                     "Growth levels, 5% tails" = "5", "Growth levels, 1% tails" = "1")
DEFAULT_OUTLIERS <- "r5"
# The outliers of a key, in words: "the 5% tails of growth reversals".
outlier_rule_words <- function(outliers) {
  sprintf("the %s%% tails of growth %s", sub("^r", "", outliers),
          if (startsWith(outliers, "r")) "reversals (change from the year before)" else "levels")
}
# The income scaling of a fitted model: form and xi of T and P.
model_income_forms <- function(m) {
  o <- model_options(list(outcome = m$outcome, income_form = m$income_form, xi = m$xi))
  list(form = o$income_form, xi = o$xi)
}
# The short label of the income scaling of a model (NULL: standard for both).
income_form_label <- function(form, xi) {
  part <- function(v) {
    f <- form[[v]]
    if (f == "standard") return("unscaled")
    scaled <- if (xi[[v]] == 0) "equal dollars" else sprintf("scaled, xi = %g", xi[[v]])
    if (f == "both") paste("+", sub("^equal dollars$", "equal-dollar terms", sub("^scaled", "scaled terms", scaled)))
    else scaled
  }
  if (all(form == "standard")) return(NULL)
  if (form[["TM"]] == form[["RR"]] && (form[["TM"]] == "standard" || xi[["TM"]] == xi[["RR"]])) return(part("TM"))
  paste0("T ", part("TM"), ", P ", part("RR"))
}

fixed_effects_formula <- function(fe, slopes, effects = "both") {
  if (effects == "year") return("year")
  s <- slopes[[fe]]
  unit <- if (length(s)) paste0("GID_1[", paste(s, collapse = ", "), "]") else "GID_1"
  if (effects == "unit") unit else paste("year +", unit)
}

# Differences in fit against noise, as in weather_arms_robustness.Rmd ("Model
# terms and effects"): Vuong test clustered by country for non-nested models,
# the clustered Wald test of the added terms for nested ones (5%).
NOISE_LEVEL <- 0.05
CRITERIA <- c(aic = "AIC", bic = "BIC", r2 = "within adj. R2")
fit_record <- function(m) {
  e <- stats::resid(m)
  s2 <- mean(e^2)
  list(aic = stats::AIC(m), bic = stats::BIC(m), war2 = fixest::fitstat(m, "war2")[[1]],
       n = stats::nobs(m), ll = -0.5 * log(2 * pi * s2) - e^2 / (2 * s2), country = m$country)
}
criterion_value <- function(rec, criterion) {
  switch(criterion, aic = -rec$aic / 2, bic = -rec$bic / 2, r2 = -(rec$n / 2) * log(1 - rec$war2))
}
fit_verdict <- function(a, b, criterion, nested_p = NULL) {
  stopifnot(a$n == b$n, identical(a$country, b$country))
  difference <- criterion_value(a, criterion) - criterion_value(b, criterion)
  beyond_noise <- if (!is.null(nested_p)) {
    !is.na(nested_p) && nested_p < NOISE_LEVEL
  } else {
    D <- tapply(a$ll - b$ll, a$country, sum)
    G <- length(D)
    se <- sqrt(G / (G - 1) * sum((D - mean(D))^2))
    is.finite(se) && se > 0 && abs(difference / se) >= stats::qnorm(1 - NOISE_LEVEL / 2)
  }
  if (!beyond_noise) "indifferent" else if (difference > 0) "better" else "worse"
}
fit_verdicts <- function(a, b, nested_p = NULL) {
  a <- fit_record(a)
  b <- fit_record(b)
  paste(vapply(names(CRITERIA), function(cr) fit_verdict(a, b, cr, nested_p), ""), collapse = " / ")
}
added_terms_p <- function(larger, smaller) {
  added <- setdiff(names(stats::coef(larger)), names(stats::coef(smaller)))
  stopifnot(length(added) > 0, all(names(stats::coef(smaller)) %in% names(stats::coef(larger))))
  b <- stats::coef(larger)[added]
  V <- stats::vcov(larger)[added, added, drop = FALSE]
  stats::pchisq(as.numeric(t(b) %*% solve(V, b)), df = length(added), lower.tail = FALSE)
}

# The options of a model (o: gdp, climate, weight, units, baseline, form (of T
# and P, or one for both), symmetric, inter = list(TM, RR), shape, fe,
# moderators, share_trailing, outcome), completed:
# the trailing share only with WB GDP and trailing moderators; the PWT outcomes
# and the capital interaction only with PWT GDP.
model_options <- function(o) {
  o$outcome <- o$outcome %||% "log"
  o$fe <- o$fe %||% "quadratic"
  o$effects <- o$effects %||% "both"
  if (o$effects == "year") o$fe <- "unit"
  o$convergence <- o$convergence %||% "none"
  # The quadratic climate control of T and of P (a single value: both).
  quad <- unlist(o$quad %||% "none")
  if (length(quad) == 1) quad <- c(TM = quad, RR = quad)
  o$quad <- c(TM = unname(quad[["TM"]]), RR = unname(quad[["RR"]]))
  # Per variable (a single value: both): symmetric terms, the units of the
  # deviations (SD or absolute) and the shape of the interactions.
  o$symmetric <- vapply(per_variable_option(o$symmetric, FALSE), isTRUE, TRUE)
  o$units <- per_variable_option(o$units, "SD")
  o$shape <- per_variable_option(o$shape, "linear")
  # The style: the regression unit itself, or its first difference.
  o$diff <- identical(o$style, "diff") || isTRUE(o$diff)
  # "both": the terms of the unit and, with arms of their own, those of its
  # first difference (dev<v>_d_<arm>).
  o$both <- !o$diff && (identical(o$style, "both") || isTRUE(o$both))
  # The lags of the regressors (0: none; at most panels$max_lags, in
  # explorer_model()).
  o$lags <- max(0L, as.integer(o$lags %||% 0L))
  o$moderators <- o$moderators %||% "fixed"
  o$share_trailing <- isTRUE(o$share_trailing) && o$moderators == "trailing" && identical(o$gdp, "WB")
  # Controls (in CONTROL_CHOICES' order) and outliers dropped ("none" or one of
  # OUTLIER_CHOICES).
  o$controls <- intersect(CONTROL_CHOICES, o$controls %||% CONTROL_CHOICES)
  o$outliers <- o$outliers %||% DEFAULT_OUTLIERS
  if (!is.null(o$gdp) && o$gdp != "PWT") {
    if (o$outcome %in% PWT_OUTCOMES) o$outcome <- "log"
    if (!is.null(o$inter)) o$inter <- lapply(o$inter, function(z) setdiff(z, PWT_MODIFIERS))
  }
  # The low-income threshold only with the low-income interaction (NA: the
  # source's default, filled in by explorer_model()).
  o$low_threshold <- if ("low" %in% unlist(o$inter)) as.numeric(o$low_threshold %||% NA)
  # The income scaling of T and of P (a single value: both), and its xi (in
  # XI_RANGE). Scaled terms with xi = 1 are the standard ones.
  per_variable <- per_variable_option
  # The arms of T and of P (a single value: both).
  o$form <- per_variable(o$form, "Linear arms")
  form <- per_variable(o$income_form, "standard")
  xi <- pmin(pmax(as.numeric(per_variable(o$xi, 0)), XI_RANGE[1]), XI_RANGE[2])
  names(xi) <- c("TM", "RR")
  form[xi == 1] <- "standard"
  xi[form == "standard"] <- 1
  o$income_form <- form
  o$xi <- xi
  # The exposure-mediated T terms (R/exposure.R; none without it): they are in
  # degree days, so the deviations are in absolute units.
  o$exposure <- if (exists("exposure_options", mode = "function")) exposure_options(o$exposure) else list(mode = "none")
  for (v in c("TM", "RR")) if (exposure_on(o$exposure, v)) o$units[[v]] <- "absolute"
  o
}
# A per-variable option c(TM, RR) from one value (both) or two.
per_variable_option <- function(x, default) {
  x <- unlist(x %||% default)
  if (length(x) == 1) x <- c(TM = x, RR = x)
  c(TM = unname(x[["TM"]]), RR = unname(x[["RR"]]))
}
# A model's value of a per-variable option for variable v (models fitted
# before the option was per variable hold one value).
of_variable <- function(x, v) if (length(x) > 1) x[[v]] else x[[1]]
# Whether any variable of a model has spline interactions.
any_splines <- function(shape) any(shape == "splines")
# Whether a model's options (or a fitted model) have exposure-mediated terms:
# of variable v ("TM", "RR"), or of either. The columns of the exposure
# channels (R/exposure.R): hot, cold (T), wet, dry (P).
exposure_on <- function(e, v = c("TM", "RR")) {
  !is.null(e) && any(vapply(v, function(z) !identical(e[[z]]$mode %||% "none", "none"), TRUE))
}
EXPOSURE_COLUMNS <- c("exh", "exc", "exw", "exd")
# The key of a pair of subsets of interactions, e.g. "T: clim+inc, P: clim+agr"
# ("none": no interaction).
interaction_key <- function(inter) {
  part <- function(x) if (length(x)) paste(x, collapse = "+") else "none"
  paste0("T: ", part(inter$TM), ", P: ", part(inter$RR))
}
model_key <- function(o, panels) {
  moderators <- if (o$share_trailing) "trailing + share" else o$moderators
  per <- function(x) if (x[["TM"]] == x[["RR"]]) as.character(x[["TM"]]) else paste0("T ", x[["TM"]], ", P ", x[["RR"]])
  paste(o$gdp, o$climate, o$baseline, o$weight, per(o$units),
        paste0("T ", o$form[["TM"]], if (o$symmetric[["TM"]]) " (symmetric)",
               ", P ", o$form[["RR"]], if (o$symmetric[["RR"]]) " (symmetric)"),
        interaction_key(o$inter), per(o$shape),
        o$fe, moderators, o$outcome,
        sprintf("income T: %s %g, P: %s %g", o$income_form[["TM"]], o$xi[["TM"]], o$income_form[["RR"]], o$xi[["RR"]]),
        paste0("controls: ", if (length(o$controls)) paste(o$controls, collapse = "+") else "none"),
        paste0("outliers dropped: ", o$outliers),
        if (!is.null(o$low_threshold)) paste0("low income below ", format(o$low_threshold)),
        if (o$lags > 0) paste0("lags: ", o$lags),
        if (o$diff) "first differences", if (o$both) "unit and first differences",
        if (o$effects != "both") paste0("effects: ", o$effects),
        if (o$convergence != "none") paste0("convergence: ", o$convergence),
        if (any(o$quad != "none")) paste0("quadratic control T: ", o$quad[["TM"]], ", P: ", o$quad[["RR"]]),
        if (exposure_on(o$exposure)) exposure_key(o$exposure), sep = " | ")
}
# The short label of the controls and outliers of a model (NULL: the
# notebook's, all controls and the 5% reversal outliers dropped).
sample_label <- function(controls, outliers) {
  dropped <- setdiff(CONTROL_CHOICES, controls)
  out <- c(if (!length(controls)) "no controls"
           else if (length(dropped)) paste("no", paste(CONTROL_SHORT[dropped], collapse = ", "), "controls"),
           if (outliers == "none") "outliers kept"
           else if (outliers != DEFAULT_OUTLIERS) paste0("outliers of ", outlier_rule_words(outliers), " dropped"))
  if (length(out)) paste(out, collapse = " · ")
}

# The verdicts of model a against model b by AIC / BIC / within adjusted R2
# (their fit records), when they run on the same rows and outcome (same panel,
# row numbers and outcome); NULL otherwise.
record_verdicts <- function(a, b) {
  if (!identical(a$panel_rows, b$panel_rows)) return(NULL)
  paste(vapply(names(CRITERIA), function(cr) fit_verdict(a$record, b$record, cr), ""), collapse = " / ")
}

# Joint Wald tests (fixest::wald, clustered) of the deviation terms of model m
# with interactions `inter`: for T and for P, all their terms; the terms
# interacted with each characteristic (both arms, with their spline terms); with
# three interactions each pair of them; with two or more, all the interactions;
# with both income forms, the plain terms (zero: equal dollars) and the terms
# times Ybar / Y (zero: equal %).
MODIFIER_TEST_LABELS <- c(clim = "climate", lag = "last year's value", cur = "this year's value", inc = "income",
                          agr = "share", kap = "capital", low = "low income")
# The column of a deviation term ("devTM_pos:clim_TM:idol") in the fit, and back.
safe_term <- function(x) gsub(":", "_X_", x, fixed = TRUE)
term_name <- function(x) gsub("_X_", ":", x, fixed = TRUE)
joint_tests <- function(m, inter, xi = c(TM = 0, RR = 0)) {
  terms <- term_name(grep("^dev", names(stats::coef(m)), value = TRUE))
  sets <- list()
  for (v in c("TM", "RR")) {
    short <- if (v == "TM") "T" else "P"
    own <- terms[substr(terms, 4, 5) == v]
    if (!length(own)) next
    sets[[paste0(short, ": all terms")]] <- own
    # With linear and quadratic arms: either set zero.
    squared <- grepl("^dev(TM|RR)_(d_)?(pos|neg|sym)2", own)
    if (any(squared) && !all(squared)) {
      sets[[paste0(short, ": squared arms (zero: linear arms suffice)")]] <- own[squared]
      sets[[paste0(short, ": linear arms (zero: quadratic arms suffice)")]] <- own[!squared]
    }
    # Style "both": the first-difference terms zero, or the unit's.
    differenced <- grepl("^dev(TM|RR)_d_", own)
    if (any(differenced)) {
      sets[[paste0(short, ": first-difference terms (zero: the unit only)")]] <- own[differenced]
      sets[[paste0(short, ": terms of the unit (zero: first differences only)")]] <- own[!differenced]
    }
    # Exposure: the terms of each channel zero.
    for (col in EXPOSURE_COLUMNS) {
      chan <- own[grepl(paste0(":", col, "(\\.l[0-9]+)?(:|$)"), own)]
      if (length(chan) && length(chan) < length(own)) {
        sets[[paste0(short, ": ", c(exh = "hot", exc = "cold", exw = "wet", exd = "dry")[[col]], " exposure terms")]] <- chan
      }
    }
    # With lags: the lagged terms zero (no lags).
    lagged <- grepl("^dev[^:]*\\.l[0-9]+", own)
    if (any(lagged)) sets[[paste0(short, ": lagged terms (zero: no lags)")]] <- own[lagged]
    scaled <- grepl(":idol(\\.l[0-9]+)?$", own)
    if (any(scaled) && !all(scaled)) {
      zero <- if (xi[[v]] == 0) "equal dollars" else sprintf("dollar loss ~ Y^%g", xi[[v]])
      sets[[paste0(short, ": plain terms (zero: ", zero, ")")]] <- own[!scaled]
      sets[[paste0(short, ": terms times (Ybar / Y)^", 1 - xi[[v]], " (zero: equal %)")]] <- own[scaled]
    }
    of <- function(z) {
      column <- modifier_column(z, v)
      own[grepl(paste0(":", column, "(_nl[0-9]+)?(\\.l[0-9]+)?(:idol(\\.l[0-9]+)?)?$"), own)]
    }
    mods <- inter[[v]]
    for (z in mods) sets[[paste0(short, ": ", MODIFIER_TEST_LABELS[[z]], " interactions")]] <- of(z)
    if (length(mods) >= 3) {
      for (pair in utils::combn(mods, 2, simplify = FALSE)) {
        sets[[paste0(short, ": ", paste(MODIFIER_TEST_LABELS[pair], collapse = " + "), " interactions")]] <-
          unlist(lapply(pair, of))
      }
    }
    if (length(mods) >= 2) sets[[paste0(short, ": all interactions")]] <- unlist(lapply(mods, of))
  }
  # Both variables with linear and squared terms.
  squared <- grepl("^dev(TM|RR)_(d_)?(pos|neg|sym)2", terms)
  mixed <- vapply(c("TM", "RR"), function(v) {
    own <- substr(terms, 4, 5) == v
    any(squared[own]) && !all(squared[own])
  }, TRUE)
  if (all(mixed)) {
    sets[["T and P: squared arms (zero: linear arms suffice)"]] <- terms[squared]
    sets[["T and P: linear arms (zero: quadratic arms suffice)"]] <- terms[!squared]
  }
  sets <- Filter(length, sets)
  do.call(rbind, lapply(names(sets), function(k) {
    w <- fixest::wald(m, keep = paste0("^", safe_term(sets[[k]]), "$"), print = FALSE)
    data.frame(test = k, stat = w$stat, p = w$p, df1 = w$df1, stringsAsFactors = FALSE)
  }))
}

# The rows of a model: the panel's without the growth outliers it drops, those
# with the trailing values it uses (.row: the row number in the panel).
model_rows <- function(o, panels) {
  d <- panels$panels[[paste(o$gdp, o$climate, o$weight, sep = " | ")]]$rows
  d$.row <- seq_len(nrow(d))
  if (o$outliers != "none") d <- d[!d[[paste0("outlier_", o$outliers)]], ]
  if (o$moderators == "trailing") d <- d[is.finite(d$inc_trailing), ]
  if (o$share_trailing) d <- d[is.finite(d$agr_trailing), ]
  if (o$outcome != "log") d <- d[is.finite(d[[panels$outcome_columns[[o$outcome]]]]), ]
  if (o$convergence != "none") d <- d[is.finite(d[[CONVERGENCE_COLUMNS[[o$convergence]]]]), ]
  if ("kap" %in% unlist(o$inter)) {
    d <- d[is.finite(d[[if (o$moderators == "trailing") "kap_trailing" else "kap_level"]]), ]
  }
  d
}

model_cache <- new.env()
explorer_model <- function(o, panels) {
  o <- model_options(o)
  if (!is.null(o$low_threshold) && !is.finite(o$low_threshold)) o$low_threshold <- low_income_thresholds(panels)[[o$gdp]]
  key <- model_key(o, panels)
  if (!is.null(model_cache[[key]])) return(model_cache[[key]])
  panel <- panels$panels[[paste(o$gdp, o$climate, o$weight, sep = " | ")]]
  window <- panels$baselines[[o$baseline]]
  trailing <- o$moderators == "trailing"
  d <- model_rows(o, panels)
  # The demeaning of the deviations of each variable (symmetric ones as they are).
  dev_fe <- c(TM = deviation_fe(o$fe, o$symmetric[["TM"]], o$effects), RR = deviation_fe(o$fe, o$symmetric[["RR"]], o$effects))
  # Exposure: the slopes of each variable's indices (R/exposure.R), on the
  # rows of the countries that have them (in both versions, so that the two
  # run on the same rows); the smooth at the unit's long-run (horizon) climate.
  uses_exposure <- exposure_on(o$exposure)
  exposure_vars <- Filter(function(v) exposure_on(o$exposure, v), c("TM", "RR"))
  exposure_tab <- NULL
  if (uses_exposure) {
    exposure_tab <- lapply(stats::setNames(exposure_vars, exposure_vars), function(v) exposure_table(o$exposure[[v]], v))
    for (v in exposure_vars) d <- d[d$GID_0 %in% exposure_tab[[v]]$units$GID_0, ]
    unit_climate <- panel$rows[!duplicated(panel$rows$GID_1), c("GID_1", "GID_0", "mean_TM_all", "mean_RR_all")]
  }
  exposure_of <- function(x) {
    u <- unit_climate[match(x$GID_1, unit_climate$GID_1), ]
    do.call(cbind, lapply(exposure_vars, function(v) {
      exposure_values(o$exposure[[v]], v, exposure_tab[[v]], u$GID_0, u[[paste0("mean_", v, "_all")]])
    }))
  }
  # First differences: the rows whose year before has the climate.
  if (o$diff || o$both) {
    for (v in c("TM", "RR")) d <- d[is.finite(deviation_values(d, v, window, panels$units_suffix[[o$units[[v]]]], dev_fe[[v]], TRUE)), ]
  }

  # The outcome (see the header).
  d$.outcome <- as.numeric(d[[if (o$outcome == "log") panels$outcome else panels$outcome_columns[[o$outcome]]]])

  # The columns the deviation terms are built from (the levels of the
  # characteristics, and where they come from).
  uses_kap <- "kap" %in% unlist(o$inter)
  uses_low <- !is.null(o$low_threshold)
  clim_window <- if (trailing && trailing_window(window)) window else "all"
  inc_column <- if (trailing) "inc_trailing" else "inc_level"
  agr_column <- if (o$share_trailing) "agr_trailing" else "agr_level"
  kap_column <- if (trailing) "kap_trailing" else "kap_level"
  # With lags: each deviation term also enters at t-1, ..., t-lags, with the
  # characteristics of that year, from the panel's rows (the outliers too) and
  # the years before the start (panels$history); the rows of the model are
  # those with all of them.
  o$lags <- min(o$lags, panels$max_lags %||% 0L)
  pool <- NULL
  if (o$lags > 0) {
    unit_level <- panel$rows[!duplicated(panel$rows$GID_1), c("GID_1", "inc_level", "agr_level", "kap_level")]
    history <- panel$history
    history[c("inc_level", "agr_level", "kap_level")] <- unit_level[match(history$GID_1, unit_level$GID_1), -1]
    common <- intersect(names(panel$rows), names(history))
    pool <- rbind(panel$rows[common], history[common])
    # The inputs of the terms in the year: the deviations, the climate, and
    # income, share and capital where they are used.
    has_inputs <- function(x) {
      ok <- rep(TRUE, nrow(x))
      for (v in c("TM", "RR")) {
        ok <- ok & is.finite(deviation_values(x, v, window, panels$units_suffix[[o$units[[v]]]], dev_fe[[v]], o$diff))
        if (o$both) ok <- ok & is.finite(deviation_values(x, v, window, panels$units_suffix[[o$units[[v]]]], dev_fe[[v]], TRUE))
        if ("clim" %in% o$inter[[v]]) ok <- ok & is.finite(x[[paste0("mean_", v, "_", clim_window)]])
        if ("lag" %in% o$inter[[v]]) ok <- ok & is.finite(x[[paste0(v, "_l1")]])
      }
      if (any(c("inc", "low") %in% unlist(o$inter)) || any(o$income_form != "standard")) ok <- ok & is.finite(x[[inc_column]])
      if ("agr" %in% unlist(o$inter)) ok <- ok & is.finite(x[[agr_column]])
      if (uses_kap) ok <- ok & is.finite(x[[kap_column]])
      ok
    }
    pool_ok <- has_inputs(pool)
    pool_key <- paste(pool$GID_1, pool$year)
    keep <- rep(TRUE, nrow(d))
    for (k in seq_len(o$lags)) keep <- keep & pool_ok[match(paste(d$GID_1, d$year - k), pool_key)] %in% TRUE
    d <- d[keep, ]
  }

  # The moderators. The climate is centred at the median of the long-run
  # climate of the panel's units (fixed or trailing), so that the arms' own
  # coefficients are the response at the median climate and at the mean income
  # and share of the rows (the cap of the app is relative to it).
  units <- panels$units[[paste(o$gdp, o$climate, o$baseline, o$weight, o$units[["TM"]], sep = " | ")]]
  clim_median <- c(TM = stats::median(units$mean_T), RR = stats::median(units$mean_P))
  inc_level <- d[[inc_column]]
  agr_level <- d[[agr_column]]
  # Log capital per capita less its mean over the rows (PWT; only when it is in
  # an interaction, the rows have it).
  kap_level <- if (!uses_kap) 0 else d[[kap_column]]
  centre <- c(inc = mean(inc_level), agr = mean(agr_level), kap = mean(kap_level))
  # The arms of the deviations demeaned as the fixed effects (symmetric: the
  # unsplit deviations as they are) and the moderators of rows x, centred by
  # the model's rows.
  with_moderators <- function(x) {
    for (v in c("TM", "RR")) {
      dev <- deviation_values(x, v, window, panels$units_suffix[[o$units[[v]]]], dev_fe[[v]], o$diff)
      for (s in arm_suffixes(o$form[[v]], o$symmetric[[v]])) x[[paste0("dev", v, s)]] <- arm_value(s, dev)
      if (o$both) {
        dev_d <- deviation_values(x, v, window, panels$units_suffix[[o$units[[v]]]], dev_fe[[v]], TRUE)
        for (s in arm_suffixes(o$form[[v]], o$symmetric[[v]])) x[[paste0("dev", v, "_d", s)]] <- arm_value(s, dev_d)
      }
      x[[paste0("clim_", v)]] <- x[[paste0("mean_", v, "_", clim_window)]] - clim_median[[v]]
      # Last year's and this year's value, centred in the same way.
      x[[paste0("lag_", v)]] <- x[[paste0(v, "_l1")]] - clim_median[[v]]
      x[[paste0("cur_", v)]] <- x[[v]] - clim_median[[v]]
    }
    x$inc <- x[[inc_column]] - centre[["inc"]]
    x$agr <- x[[agr_column]] - centre[["agr"]]
    x$kap <- if (!uses_kap) 0 else x[[kap_column]] - centre[["kap"]]
    # The low-income dummy: the income of the interactions below the threshold.
    x$low <- if (uses_low) as.numeric(x[[inc_column]] < log(o$low_threshold)) else 0
    # (Ybar / Y_i)^(1 - xi_v), Ybar the geometric mean income of the rows.
    for (v in c("TM", "RR")) x[[paste0("idol_", v)]] <- exp(-(1 - o$xi[[v]]) * x$inc)
    # The exposure of the terms (per 100 of the index per unit of deviation).
    if (uses_exposure) {
      ev <- exposure_of(x)
      for (z in colnames(ev)) x[[z]] <- ev[, z]
    }
    x
  }
  d <- with_moderators(d)
  # The non-linear spline terms <z>_nl<j> (knots at the spline percentiles of
  # the units' mean values), less their value at zero.
  spline_columns <- setdiff(SPLINE_COLUMNS, if (!uses_kap) "kap")
  unit_values <- stats::aggregate(d[spline_columns], list(GID_1 = d$GID_1), mean)
  knots <- lapply(stats::setNames(spline_columns, spline_columns), function(z) {
    stats::quantile(unit_values[[z]], panels$spline_knot_probs, names = FALSE)
  })
  with_splines_of <- function(x) {
    for (z in spline_columns) {
      nl <- rcs_nonlinear(x[[z]], knots[[z]])
      at_zero <- rcs_nonlinear(0, knots[[z]])
      for (j in seq_len(ncol(nl))) x[[paste0(z, "_nl", j)]] <- nl[, j] - at_zero[1, j]
    }
    x
  }
  d <- with_splines_of(d)
  # The lagged columns <column>.l<k>, from the pool's year t-k.
  if (o$lags > 0) {
    pool <- with_splines_of(with_moderators(pool))
    lag_columns <- c(grep("^dev(TM|RR)_", names(d), value = TRUE),
                     intersect(c(spline_columns, grep("_nl[0-9]+$", names(d), value = TRUE)), names(pool)),
                     "low", "idol_TM", "idol_RR", if (uses_exposure) intersect(EXPOSURE_COLUMNS, names(d)))
    for (k in seq_len(o$lags)) {
      at <- match(paste(d$GID_1, d$year - k), pool_key)
      for (z in lag_columns) d[[paste0(z, ".l", k)]] <- pool[[z]][at]
    }
  }

  # The time-varying moderators in at least one interaction (their levels).
  varying <- c(if (trailing && trailing_window(window) && "clim" %in% o$inter$TM) "clim_TM",
               if (trailing && trailing_window(window) && "clim" %in% o$inter$RR) "clim_RR",
               if (trailing && "inc" %in% unlist(o$inter)) "inc",
               if (trailing && uses_kap) "kap",
               if (trailing && uses_low) "low",
               if (o$share_trailing && "agr" %in% unlist(o$inter)) "agr",
               # Last year's and this year's value always vary.
               unlist(lapply(c("TM", "RR"), function(v) {
                 vapply(intersect(c("lag", "cur"), o$inter[[v]]), modifier_column, "", v, USE.NAMES = FALSE)
               })))
  # (Ybar / Y)^(1 - xi) varies with the trailing income: its level too (once
  # when T and P are scaled by the same power).
  scaled_v <- names(o$income_form)[o$income_form != "standard"]
  idol_levels <- if (length(scaled_v)) paste0("idol_", scaled_v[!duplicated(o$xi[scaled_v])])
  varying_idol <- if (trailing) idol_levels
  # Without unit effects, the constant characteristics in an interaction (and
  # the constant scaling) enter with their levels too, once (not lagged).
  constant <- if (o$effects == "year") {
    setdiff(c(if ("clim" %in% o$inter$TM) "clim_TM", if ("clim" %in% o$inter$RR) "clim_RR",
              if ("inc" %in% unlist(o$inter)) "inc", if ("agr" %in% unlist(o$inter)) "agr",
              if (uses_kap) "kap", if (uses_low) "low"), varying)
  }
  constant_idol <- if (o$effects == "year" && !trailing) idol_levels
  # The convergence control (its column conv_gdp).
  if (o$convergence != "none") d$conv_gdp <- d[[CONVERGENCE_COLUMNS[[o$convergence]]]]
  # The quadratic climate controls (columns quad_<v> and quad_<v>_sq), centred
  # as the climate moderator.
  quad_terms <- character(0)
  for (v in names(o$quad)[o$quad != "none"]) {
    z <- d[[quad_column(o$quad[[v]], v)]] - clim_median[[v]]
    d[[paste0("quad_", v)]] <- z
    d[[paste0("quad_", v, "_sq")]] <- z^2
    quad_terms <- c(quad_terms, paste0("quad_", v), paste0("quad_", v, "_sq"))
  }
  # The spline terms of the continuous moderators (not of the low-income dummy).
  with_splines <- function(z, shape) {
    if (shape == "splines") c(z, unlist(lapply(intersect(z, spline_columns), function(x) {
      paste0(x, "_nl", seq_len(length(panels$spline_knot_probs) - 2))
    }))) else z
  }
  # A term, or a level, at lag k: each part ends in .l<k>.
  lag_of <- function(x, k) {
    vapply(strsplit(x, ":", fixed = TRUE), function(p) paste(paste0(p, ".l", k), collapse = ":"), "")
  }
  # The level terms have spline terms when any variable has spline interactions.
  varying_levels <- function(shape) c(with_splines(varying, if (any_splines(shape)) "splines" else "linear"), varying_idol)
  levels_of <- function(shape) c(varying_levels(shape), with_splines(constant, if (any_splines(shape)) "splines" else "linear"), constant_idol)
  terms_of <- function(shape) {
    deviation <- unlist(lapply(c("TM", "RR"), function(v) {
      arms <- paste0("dev", v, arm_suffixes(o$form[[v]], o$symmetric[[v]]))
      if (o$both) arms <- c(arms, paste0("dev", v, "_d", arm_suffixes(o$form[[v]], o$symmetric[[v]])))
      # Exposure: every arm through the exposure of each of the variable's
      # channels (T: hot, cold; P: wet, dry), no plain arm; the interactions
      # then multiply the exposed arms.
      if (v %in% exposure_vars) arms <- unlist(lapply(exposure_columns(o$exposure, v), function(e) paste0(arms, ":", e)))
      modifiers <- with_splines(vapply(o$inter[[v]], modifier_column, "", v, USE.NAMES = FALSE), of_variable(shape, v))
      plain <- c(arms, unlist(lapply(modifiers, function(z) paste0(arms, ":", z))))
      switch(o$income_form[[v]], standard = plain, dollars = paste0(plain, ":idol"),
             both = c(plain, paste0(plain, ":idol")))
    }))
    levels <- varying_levels(shape)
    lagged <- unlist(lapply(seq_len(o$lags), function(k) c(lag_of(deviation, k), if (length(levels)) paste0(levels, ".l", k))))
    c(deviation, levels_of(shape), lagged, if (o$convergence != "none") "conv_gdp", quad_terms,
      unlist(panels$control_sets[o$controls], use.names = FALSE))
  }
  # Each deviation term enters as a column of its own (the product of its
  # parts), named safe_term(): in a formula, R orders the variables of an
  # interaction by their first appearance, which would scramble the names. The
  # part "idol" (or "idol.l<k>") is idol_TM or idol_RR, by the term's variable.
  fit <- function(shape) {
    terms <- terms_of(shape)
    for (t in grep("^dev", terms, value = TRUE)) {
      d[[safe_term(t)]] <- Reduce(`*`, lapply(strsplit(t, ":", fixed = TRUE)[[1]], function(x) {
        d[[if (startsWith(x, "idol")) sub("^idol", paste0("idol_", substr(t, 4, 5)), x) else x]]
      }))
    }
    terms[startsWith(terms, "dev")] <- safe_term(terms[startsWith(terms, "dev")])
    f <- stats::as.formula(paste(".outcome ~", paste(terms, collapse = " + "), "|",
                                 fixed_effects_formula(o$fe, panels$fe_slopes, o$effects)))
    m <- fixest::feols(f, data = d, cluster = ~GID_0, panel.id = ~GID_1 + year, notes = FALSE)
    m$country <- d$GID_0[fixest::obs(m)]
    m
  }
  m <- fit(o$shape)
  # The splines against the linear interactions (same rows): the Wald test of
  # the non-linear terms and the verdict by AIC / BIC / within adjusted R2.
  spline_test <- NULL
  if (any_splines(o$shape) && any(grepl("_nl", names(stats::coef(m))))) {
    linear <- fit("linear")
    nl_p <- added_terms_p(m, linear)
    spline_test <- list(nonlinear_p = nl_p, verdict = fit_verdicts(m, linear, nested_p = nl_p),
                        aic = stats::AIC(m) - stats::AIC(linear), bic = stats::BIC(m) - stats::BIC(linear))
  }
  ct <- fixest::coeftable(m)
  rownames(ct) <- term_name(rownames(ct))
  keep <- grep("^dev", rownames(ct), value = TRUE)
  lagged_levels <- varying_levels(o$shape)
  level_terms <- intersect(c(levels_of(o$shape), unlist(lapply(seq_len(o$lags), function(k) {
    if (length(lagged_levels)) paste0(lagged_levels, ".l", k)
  }))), rownames(ct))
  coefs <- stats::setNames(stats::coef(m), term_name(names(stats::coef(m))))
  V <- stats::vcov(m)
  dimnames(V) <- list(term_name(rownames(V)), term_name(colnames(V)))
  model <- list(
    key = key, gdp = o$gdp, climate = o$climate, weight = o$weight, units = o$units,
    baseline = o$baseline, window = window, diff = o$diff, both = o$both, dev_fe = dev_fe,
    form = o$form, symmetric = o$symmetric, lags = o$lags,
    inter = o$inter, shape = o$shape,
    # The coefficients g and h of the NRK filter of T and P (its baseline only).
    filter = if (window == "nrk") panels$filters[[paste(o$gdp, o$climate, o$weight, sep = " | ")]],
    fe = o$fe, effects = o$effects, moderators = o$moderators, share_trailing = o$share_trailing,
    # The quadratic climate controls and their coefficients (NULL: none).
    quad = o$quad,
    quad_coef = if (length(intersect(quad_terms, rownames(ct)))) {
      k <- intersect(quad_terms, rownames(ct))
      data.frame(term = k, estimate = ct[k, 1], se = ct[k, 2], p = ct[k, 4], row.names = NULL)
    },
    # The convergence control and its coefficient (NULL: none).
    convergence = o$convergence,
    convergence_coef = if (o$convergence != "none" && "conv_gdp" %in% rownames(ct)) {
      c(estimate = ct["conv_gdp", 1], se = ct["conv_gdp", 2], p = ct["conv_gdp", 4])
    },
    outcome = o$outcome, income_form = o$income_form, xi = o$xi, controls = o$controls, outliers = o$outliers,
    # The exposure of the T terms and its table of slopes (NULL: none).
    exposure = o$exposure, exposure_table = exposure_tab,
    coef = coefs[keep], vcov = V[keep, keep, drop = FALSE],
    # The coefficients' p-values (clustered), the fit, the joint tests, and
    # what the noise tests of differences in fit need (its rows of the panel).
    p = ct[keep, 4],
    fit = c(war2 = unname(fixest::fitstat(m, "war2")[[1]]), aic = stats::AIC(m), bic = stats::BIC(m)),
    tests = joint_tests(m, o$inter, o$xi),
    record = fit_record(m),
    panel_rows = c(paste(o$gdp, o$climate, o$weight, o$outcome), d$.row[fixest::obs(m)]),
    levels = data.frame(term = level_terms, estimate = ct[level_terms, 1], se = ct[level_terms, 2],
                        p = ct[level_terms, 4], row.names = NULL),
    nobs = stats::nobs(m), countries = length(unique(m$country)),
    clim_median = clim_median, inc_mean = centre[["inc"]], agr_mean = centre[["agr"]],
    kap_mean = if (uses_kap) mean(kap_level), low_threshold = o$low_threshold,
    horizon = panel$horizon,
    # The windows of the trailing income and share (NULL: fixed), for the
    # projections that let them follow the SSP (capital has no SSP path: it
    # stays at the country's value).
    inc_window = if (trailing) panels$trailing_income_years,
    share_window = if (o$share_trailing) panels$trailing_share_years,
    knots = if (any_splines(o$shape)) knots, spline_test = spline_test
  )
  model_cache[[key]] <- model
  model
}

# The rows of the fixed and of the trailing models of a panel, and with the
# trailing share (for the warning on its sample), without the growth outliers
# `outliers` (as in model_options()).
trailing_rows <- function(gdp, climate, weight, panels, outliers = DEFAULT_OUTLIERS) {
  d <- panels$panels[[paste(gdp, climate, weight, sep = " | ")]]$rows
  if (outliers != "none") d <- d[!d[[paste0("outlier_", outliers)]], ]
  count <- function(x) c(rows = sum(x), countries = length(unique(d$GID_0[x])))
  rbind(fixed = count(rep(TRUE, nrow(d))), trailing = count(is.finite(d$inc_trailing)),
        share = count(is.finite(d$inc_trailing) & is.finite(d$agr_trailing)))
}
