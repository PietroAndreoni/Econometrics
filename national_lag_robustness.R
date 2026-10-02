# Robustness of the climate-response models on national (WB, PWT) and
# subnational (DOSE) data, under one of the specifications chosen by the run
# option SPECIFICATION:
#
#   "lags" (default; moved from the lags-national chunk of main_analysis.Rmd)
#     l(dTM, 0:10) + l(dTM, 0:10):mean_TM_all + l(dRR, 0:10) + l(dRR, 0:10):mean_RR_all
#       + l(zTM^2, 0:10) + l(zRR^2, 0:10)
#     dTM and dRR are year-on-year changes. Output: results/national_lag_robustness/
#
#   "bhmdev" (the contemporaneous BHM + deviations model of main_analysis.Rmd)
#     TM + TM:mean_TM_all + RR + RR:mean_RR_all + zTM^2 + zRR^2
#     Output: results/national_bhmdev_robustness/
#
#   "bhmdev_celsius" the "bhmdev" model with zTM in degrees Celsius
#     (zTM = TM - baseline mean, no SD scaling); zRR stays standardized.
#     Output: results/national_bhmdev_celsius_robustness/
#
#   "bhmsplit" the "bhmdev" model with each anomaly split into its positive
#     and negative part: abs_zTMp^2 + abs_zTMm^2 + abs_zRRp^2 + abs_zRRm^2 in
#     place of zTM^2 + zRR^2. Adds tests of warm/cold symmetry and of equal
#     anomaly responses in hot and cold units (long-run mean above or below
#     HOT_THRESHOLD). Output: results/national_bhmsplit_robustness/
#
#   "bhmsplitx" the "bhmsplit" model with each split anomaly interacted with
#     its variable's long-run climate (abs_zTMp^2:mean_TM_all, ...,
#     abs_zRRm^2:mean_RR_all), so warm, cold, wet and dry anomaly effects vary
#     with the climate; they are evaluated at AT_VALUES, with the climate
#     where each changes sign. Output: results/national_bhmsplitx_robustness/
#
#   "anomx" the anomaly terms of "bhmsplitx" alone, without temperature and
#     precipitation levels. Compared on the same sample with the levels-only
#     BHM and with the combined "bhmsplitx" model (test_levels_vs_anomalies).
#     Output: results/national_anomx_robustness/
#
#   "bhmquad" the "bhmdev" model with a quadratic in the current climate in
#     place of the interactions with the long-run climate:
#     TM + TM^2 + RR + RR^2 + zTM^2 + zRR^2. The temperature effect at T is
#     b1 + 2 * b2 * T, evaluated at the temperatures of AT_VALUES (and the same
#     for precipitation). Output: results/national_bhmquad_robustness/
#
# Run with  Rscript national_lag_robustness.R bhmdev  or set SPECIFICATION
# before sourcing the script.
#
# mean_TM_all and mean_RR_all are each unit's (country's or region's) mean over
# its entire climate series. The anomaly definitions are set in DEVIATIONS (see
# add_weather_variables() in functions/build_dat.R):
#   k-year MA       zTM = (TM - mean_TM_lag_k) / sd_TM_all
#   k-year SD       zTM = (TM - mean_TM_lag_k) / sd_TM_lag_k
#   Hamilton trend  zTM = (TM - trend_TM) / sd_TM_trend (SD of the cycle)
# and the same for RR.
#
# Grid: WB, PWT and DOSE x ERA5, CRU TS and UDelaware x area and population
# weights (CLIMATE_WEIGHTS) x the anomaly definitions. Each uses all
# available years: the economic sample
# starts at the first year of the economic series, and MAX_LAG + 1 climate
# years are kept before it so that every lag resolves. Fixed effects as in
# main_analysis.Rmd (year + GID_1 + GID_1[year]); errors clustered by country
# (GID_0) in every model.
#
# Claims tested in every model. "Effect" is the cumulative effect through lag
# 10 under "lags" and the contemporaneous effect under "bhmdev":
#   (a) cold countries benefit from warming, hot countries are harmed: the
#       temperature effect (dTM or TM) is positive at a long-run mean
#       temperature (mean_TM_all) of 5 degrees Celsius and negative at 25; the
#       (summed) interaction coefficient is the slope of that effect in
#       mean_TM_all. All interacted effects are evaluated at the long-run
#       climates of AT_VALUES (5, 15, 25 degrees Celsius; 0.5, 1.5, 2.5 metres
#       of precipitation per year), the same in every model.
#   (b) squared anomalies are harmful: the effects of zTM^2 and zRR^2 are
#       negative.
#   (c) anomalies add information: a joint clustered Wald test that all zTM^2
#       and zRR^2 coefficients are zero, and the fit against the model without
#       them on the identical sample.
#
# Output: tables as csv and analysis.rds, figures as png.

source("load_functions.R")

# Run option -------------------------------------------------------------------

if (!exists("SPECIFICATION")) {
  SPECIFICATION <- commandArgs(trailingOnly = TRUE)[1]
  if (is.na(SPECIFICATION)) SPECIFICATION <- "lags"
}
if (!SPECIFICATION %in% c("lags", "bhmdev", "bhmdev_celsius", "bhmquad",
                          "bhmsplit", "bhmsplitx", "anomx")) {
  stop("SPECIFICATION must be \"lags\", \"bhmdev\", \"bhmdev_celsius\", ",
       "\"bhmquad\", \"bhmsplit\", \"bhmsplitx\" or \"anomx\".")
}
QUADRATIC <- SPECIFICATION == "bhmquad"
# Split anomaly terms (the split options); SPLIT adds the symmetry and
# hot/cold tests, ANOMALY_X the interactions of the anomalies with climate,
# NO_LEVELS drops the temperature and precipitation levels ("anomx").
SPLIT_TERMS <- SPECIFICATION %in% c("bhmsplit", "bhmsplitx", "anomx")
SPLIT <- SPECIFICATION == "bhmsplit"
ANOMALY_X <- SPECIFICATION %in% c("bhmsplitx", "anomx")
NO_LEVELS <- SPECIFICATION == "anomx"

# Configuration ----------------------------------------------------------------

# WB and PWT are national; DOSE is subnational (GADM1 regions, default DOSE
# GDP definition). Errors are clustered by country (GID_0) in every model,
# which for the national panels is the panel unit itself.
ECON_DATASETS <- c(WB = "WB", PWT = "PWT", DOSE = "DOSE")
CLIMATE_SOURCES <- c(ERA5 = "ERA5", `CRU TS` = "CRU TS", UDelaware = "UDelaware")
# Climate aggregation weights: grid-cell area, or each year's population.
CLIMATE_WEIGHTS <- c(area = "area", pop = "concurrent population")
DEVIATIONS <- list(
  `10-year MA` = list(baseline = "lag_10", sd_baseline = "all"),
  `20-year MA` = list(baseline = "lag_20", sd_baseline = "all"),
  `30-year MA` = list(baseline = "lag_30", sd_baseline = "all"),
  `10-year SD` = list(baseline = "lag_10", sd_baseline = NULL),
  `20-year SD` = list(baseline = "lag_20", sd_baseline = NULL),
  `30-year SD` = list(baseline = "lag_30", sd_baseline = NULL),
  `Hamilton trend` = list(baseline = "trend", sd_baseline = NULL)
)

OUTCOME <- "dlgrp_pc_usd"
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"
PANEL_ID <- c("GID_1", "year")
CONF_LEVEL <- 0.95

# Long-run climates at which the interacted effects are evaluated, the same in
# every model: mean_TM_all in degrees Celsius and mean_RR_all in metres of
# precipitation per year, coded low/mid/high.
AT_CODES <- c("low", "mid", "high")
AT_VALUES <- list(
  mean_TM_all = stats::setNames(c(5, 15, 25), AT_CODES),
  mean_RR_all = stats::setNames(c(0.5, 1.5, 2.5), AT_CODES)
)
DEG_C <- "°C"
AT_LABELS <- stats::setNames(
  paste0(AT_VALUES$mean_TM_all, " ", DEG_C, " / ", AT_VALUES$mean_RR_all, " m"),
  AT_CODES
)

# Squared anomaly terms by variable. Under "bhmsplit" each anomaly is split
# into its positive (warm, wet) and negative (cold, dry) part, abs_zTMp =
# max(zTM, 0) and abs_zTMm = max(-zTM, 0), each entering squared.
SQUARED_TERMS <- if (SPLIT_TERMS) {
  list(TM = c("abs_zTMp^2", "abs_zTMm^2"), RR = c("abs_zRRp^2", "abs_zRRm^2"))
} else {
  list(TM = "zTM^2", RR = "zRR^2")
}
# Under "bhmsplitx" each squared anomaly is interacted with its variable's
# long-run climate, so its effect is b + g * mean_TM_all (or mean_RR_all).
ANOMALY_INTERACTIONS <- c(TM = "mean_TM_all", RR = "mean_RR_all")
# Long-run mean temperature (degrees Celsius) above which a unit counts as hot
# in the "bhmsplit" test of equal anomaly responses in hot and cold units.
HOT_THRESHOLD <- 15

# Terms of each specification. CLIMATE_VARS are the climate regressors
# interacted with the long-run climate; coef_names() gives the coefficient
# names of a term (one per lag).
if (SPECIFICATION == "lags") {
  MAX_LAG <- 10L
  CLIMATE_VARS <- c(TM = "dTM", RR = "dRR")
  RESTRICTED_TERMS <- sprintf(
    paste("l(dRR, 0:%1$d) + l(dTM, 0:%1$d)",
          "+ l(dRR, 0:%1$d):mean_RR_all + l(dTM, 0:%1$d):mean_TM_all"),
    MAX_LAG
  )
  DEVIATION_TERMS <- sprintf("l(zTM^2, 0:%1$d) + l(zRR^2, 0:%1$d)", MAX_LAG)
  coef_names <- function(term) sprintf("l(%s, %d)", term, 0:MAX_LAG)
  EFFECT_LABEL <- "Ten-year cumulative effect"
  MODEL_LABEL <- "lag model"
  OUTPUT_DIR <- file.path("results", "national_lag_robustness")
} else {
  MAX_LAG <- 0L
  CLIMATE_VARS <- c(TM = "TM", RR = "RR")
  RESTRICTED_TERMS <- if (QUADRATIC) {
    "TM + TM^2 + RR + RR^2"
  } else {
    "TM + TM:mean_TM_all + RR + RR:mean_RR_all"
  }
  DEVIATION_TERMS <- paste(unlist(SQUARED_TERMS), collapse = " + ")
  if (ANOMALY_X) {
    DEVIATION_TERMS <- paste(
      DEVIATION_TERMS, "+",
      paste(unlist(lapply(names(SQUARED_TERMS), function(v) {
        paste0(SQUARED_TERMS[[v]], ":", ANOMALY_INTERACTIONS[[v]])
      })), collapse = " + ")
    )
  }
  coef_names <- function(term) {
    if (grepl("\\^", term)) paste0("I(", term, ")") else term
  }
  EFFECT_LABEL <- "Contemporaneous effect"
  MODEL_LABEL <- if (QUADRATIC) {
    "quadratic BHM + deviations model"
  } else if (NO_LEVELS) {
    "split deviations x climate model (no levels)"
  } else if (ANOMALY_X) {
    "BHM + split deviations x climate model"
  } else if (SPLIT) {
    "BHM + split deviations model"
  } else {
    "BHM + deviations model"
  }
  OUTPUT_DIR <- file.path(
    "results",
    if (QUADRATIC) {
      "national_bhmquad_robustness"
    } else if (NO_LEVELS) {
      "national_anomx_robustness"
    } else if (ANOMALY_X) {
      "national_bhmsplitx_robustness"
    } else if (SPLIT) {
      "national_bhmsplit_robustness"
    } else {
      "national_bhmdev_robustness"
    }
  )
}
# The estimated model: the restricted terms plus the anomalies, or under
# "anomx" the anomalies alone. There the restricted (levels-only BHM) model is
# a non-nested comparison on the same sample, and a combined model with both
# (the "bhmsplitx" model) tests whether either block adds to the other.
FULL_TERMS <- if (NO_LEVELS) {
  DEVIATION_TERMS
} else {
  paste(RESTRICTED_TERMS, "+", DEVIATION_TERMS)
}
# Units of the anomalies (see `deviation` in add_weather_variables()): under
# "bhmdev_celsius" zTM is in degrees Celsius and only zRR is standardized, so
# the SD part of each DEVIATIONS entry applies to zRR alone.
DEVIATION_UNITS <- c(TM = "standardized", RR = "standardized")
if (SPECIFICATION == "bhmdev_celsius") {
  DEVIATION_UNITS[["TM"]] <- "absolute"
  MODEL_LABEL <- "BHM + deviations model (zTM in degrees Celsius)"
  OUTPUT_DIR <- file.path("results", "national_bhmdev_celsius_robustness")
}
CLIM_HISTORY_YEARS <- MAX_LAG + 1L
INTERACTIONS <- stats::setNames(c("mean_TM_all", "mean_RR_all"), CLIMATE_VARS)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

fit_terms <- function(terms, model_data) {
  fit_panel_model(terms, model_data, outcome = OUTCOME,
                  fixed_effects = FIXED_EFFECTS, panel_id = PANEL_ID,
                  cluster = ~GID_0)
}

sum_weights <- function(names) stats::setNames(rep(1, length(names)), names)
# Regular expression matching exactly the given coefficient names.
exact_keep <- function(names) {
  paste0("^(", paste(gsub("([][{}()^$.|*+?\\\\])", "\\\\\\1", names),
                     collapse = "|"), ")$")
}

# Clustered Wald test of the linear restrictions R b = 0 (one row of R per
# restriction, columns named by coefficient), chi-squared with nrow(R) df.
wald_linear <- function(model, R) {
  terms <- colnames(R)
  b <- stats::coef(model)[terms]
  V <- stats::vcov(model)[terms, terms, drop = FALSE]
  r <- R %*% b
  stat <- drop(t(r) %*% solve(R %*% V %*% t(R)) %*% r)
  stats::pchisq(stat, df = nrow(R), lower.tail = FALSE)
}

# Restriction rows: coefficient on `a` equals coefficient on `b`.
equality_row <- function(a, b, terms) {
  stats::setNames(as.numeric(terms == a) - as.numeric(terms == b), terms)
}

# Name of the coefficient on the interaction of `a` and `b`, in whichever
# order the model reports it.
interaction_name <- function(a, b, terms) {
  hit <- intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)
  if (!length(hit)) stop("No interaction coefficient for ", a, " x ", b)
  hit[[1]]
}

# Variable (TM or RR) of a squared anomaly term.
squared_variable <- function(term) {
  names(SQUARED_TERMS)[vapply(SQUARED_TERMS, function(x) term %in% x, logical(1))]
}

# Coefficient weights of the effect of climate variable `v` at climate
# `value`: b + g * value with the interaction on the long-run mean, or under
# "bhmquad" the derivative of the quadratic, b1 + 2 * b2 * value.
effect_weights <- function(v, value) {
  if (QUADRATIC) {
    stats::setNames(c(1, 2 * value), c(v, coef_names(paste0(v, "^2"))))
  } else {
    stats::setNames(c(1, value), c(v, paste0(v, ":", INTERACTIONS[[v]])))
  }
}

# Effects of the contemporaneous model, in the layout of
# distributed_lag_effects() with lag 0 as both horizons: the climate variables
# at each `at` value (effect_weights()), and the squared anomalies.
contemporaneous_effects <- function(model, at) {
  grid <- bind_rows(
    # No level effects without level terms ("anomx").
    if (!NO_LEVELS) {
      lapply(names(INTERACTIONS), function(v) {
        variable <- INTERACTIONS[[v]]
        values <- at[[variable]]
        tibble::tibble(expression = v, at = names(values),
                       interaction = variable, at_value = unname(values))
      })
    },
    if (ANOMALY_X) {
      # Each squared anomaly at each long-run climate of its variable.
      bind_rows(lapply(unlist(SQUARED_TERMS), function(term) {
        variable <- ANOMALY_INTERACTIONS[[squared_variable(term)]]
        values <- at[[variable]]
        tibble::tibble(expression = term, at = names(values),
                       interaction = variable, at_value = unname(values))
      }))
    } else {
      tibble::tibble(expression = unlist(SQUARED_TERMS), at = "NA",
                     interaction = NA_character_, at_value = NA_real_)
    }
  )
  terms <- names(stats::coef(model))
  G <- t(vapply(seq_len(nrow(grid)), function(i) {
    w <- stats::setNames(numeric(length(terms)), terms)
    e <- grid$expression[[i]]
    if (is.na(grid$interaction[[i]])) {
      w[[coef_names(e)]] <- 1
    } else if (e %in% unlist(SQUARED_TERMS)) {
      w[[coef_names(e)]] <- 1
      w[[interaction_name(coef_names(e), grid$interaction[[i]], terms)]] <-
        grid$at_value[[i]]
    } else {
      weights <- effect_weights(e, grid$at_value[[i]])
      w[names(weights)] <- weights
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  effects <- bind_cols(grid, linear_combination(model, G, CONF_LEVEL))
  bind_rows(mutate(effects, horizon = "Effect at lag"),
            mutate(effects, horizon = "Cumulative through lag")) %>%
    transmute(
      expression = factor(expression, levels = unique(grid$expression)),
      lag = 0L,
      horizon = factor(horizon, levels = c("Effect at lag",
                                           "Cumulative through lag")),
      at, interaction, at_value,
      effect = estimate, standard_error = std_error,
      lower = conf_low, upper = conf_high
    )
}

model_effects <- function(model, at) {
  if (MAX_LAG > 0L) {
    distributed_lag_effects(model, CONF_LEVEL, interactions = INTERACTIONS,
                            at = at)
  } else {
    contemporaneous_effects(model, at)
  }
}

# Panels and models ----------------------------------------------------------------

specs <- tibble::as_tibble(expand.grid(
  climate = names(CLIMATE_SOURCES),
  weight = names(CLIMATE_WEIGHTS),
  econ = names(ECON_DATASETS),
  deviation = names(DEVIATIONS),
  stringsAsFactors = FALSE
)) %>%
  mutate(id = paste(econ, climate, weight, deviation, sep = " | "))

first_year <- vapply(ECON_DATASETS, function(econ) {
  min(load_econ_panel(resolve_econ_source(econ))$year)
}, integer(1))

panels <- list()
fits_full <- list()
fits_restricted <- list()
hot_tests <- list()
level_tests <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  message("Fitting ", s$id)
  cfg <- panel_config(
    econ_year_min = first_year[[s$econ]],
    climate_source = CLIMATE_SOURCES[[s$climate]],
    climate_weight = CLIMATE_WEIGHTS[[s$weight]],
    clim_history_years = CLIM_HISTORY_YEARS
  )
  panel <- suppressMessages(build_dat(
    ECON_DATASETS[[s$econ]], config = cfg,
    baseline = DEVIATIONS[[s$deviation]]$baseline,
    sd_baseline = DEVIATIONS[[s$deviation]]$sd_baseline,
    deviation = DEVIATION_UNITS
  ))
  full <- fit_terms(FULL_TERMS, panel)

  # The model without anomalies on the identical sample: the outcome is
  # withheld outside the full model's sample, and the lag carriers are kept.
  restricted_data <- panel
  outside <- setdiff(seq_len(nrow(panel)), fixest::obs(full))
  restricted_data[[OUTCOME]][outside] <- NA_real_
  restricted <- fit_terms(RESTRICTED_TERMS, restricted_data)
  stopifnot(identical(sort(fixest::obs(restricted)), sort(fixest::obs(full))))

  # "anomx": the combined model with levels and anomalies on the same sample.
  # Wald tests of each block given the other, and the fit of the
  # anomalies-only model against the combined one.
  if (NO_LEVELS) {
    combined <- fit_terms(paste(RESTRICTED_TERMS, "+", DEVIATION_TERMS),
                          restricted_data)
    combined_terms <- names(stats::coef(combined))
    level_names <- intersect(combined_terms, names(stats::coef(restricted)))
    anomaly_names <- setdiff(combined_terms, level_names)
    level_tests[[s$id]] <- tibble::tibble(
      id = s$id,
      levels_add_p = wald_row(combined, exact_keep(level_names))$p_value,
      anomalies_add_p = wald_row(combined, exact_keep(anomaly_names))$p_value,
      within_r2_combined = unname(fixest::fitstat(combined, "wr2")[[1]]),
      delta_aic_vs_combined = stats::AIC(full) - stats::AIC(combined),
      delta_bic_vs_combined = stats::BIC(full) - stats::BIC(combined)
    )
  }

  # "bhmsplit": do hot and cold units respond alike to anomalies? The full
  # model plus each split anomaly term times a hot-unit indicator, on the same
  # sample; the interactions are the hot-minus-cold differences.
  if (SPLIT) {
    hot <- as.numeric(restricted_data$mean_TM_all > HOT_THRESHOLD)
    hot_cols <- c(TMp = "abs_zTMp", TMm = "abs_zTMm", RRp = "abs_zRRp", RRm = "abs_zRRm")
    for (k in names(hot_cols)) {
      restricted_data[[paste0("hot_", k)]] <- restricted_data[[hot_cols[[k]]]]^2 * hot
    }
    hot_terms <- paste0("hot_", names(hot_cols))
    hot_fit <- fit_terms(paste(RESTRICTED_TERMS, "+", DEVIATION_TERMS, "+",
                               paste(hot_terms, collapse = " + ")),
                         restricted_data)
    hot_tests[[s$id]] <- tibble::tibble(
      id = s$id,
      hot_p = wald_row(hot_fit, exact_keep(hot_terms))$p_value,
      hot_p_TM = wald_row(hot_fit, exact_keep(hot_terms[1:2]))$p_value,
      hot_p_RR = wald_row(hot_fit, exact_keep(hot_terms[3:4]))$p_value,
      zTMp2_cold = stats::coef(hot_fit)[["I(abs_zTMp^2)"]],
      zTMp2_hot = sum(stats::coef(hot_fit)[c("I(abs_zTMp^2)", "hot_TMp")]),
      zTMp2_hot_minus_cold_p = 2 * stats::pnorm(-abs(
        stats::coef(hot_fit)[["hot_TMp"]] / sqrt(stats::vcov(hot_fit)["hot_TMp", "hot_TMp"])
      )),
      hot_units = n_distinct(restricted_data$GID_1[fixest::obs(full)][hot[fixest::obs(full)] == 1])
    )
  }

  # Only the estimation sample's identifiers and anomalies are kept, so the
  # grid fits in memory.
  panels[[s$id]] <- panel[fixest::obs(full), c("GID_0", "GID_1", "year", "zTM", "zRR")]
  fits_full[[s$id]] <- full
  fits_restricted[[s$id]] <- restricted
}

used_rows <- function(id) panels[[id]]

# Part 1: analysis tables ----------------------------------------------------------

sample_summary <- bind_rows(lapply(specs$id, function(id) {
  used <- used_rows(id)
  tibble::tibble(
    id = id, observations = nrow(used), units = n_distinct(used$GID_1),
    countries = n_distinct(used$GID_0),
    first_year = min(used$year), last_year = max(used$year),
    max_abs_zTM = max(abs(used$zTM)), max_abs_zRR = max(abs(used$zRR))
  )
})) %>%
  left_join(specs, by = "id")
write_out(sample_summary, "sample_summary")

# Effects (lag paths under "lags") at the long-run climates of AT_VALUES.
lag_effects <- bind_rows(lapply(specs$id, function(id) {
  model_effects(fits_full[[id]], AT_VALUES) %>%
    mutate(id = id, .before = 1)
})) %>%
  left_join(specs, by = "id") %>%
  mutate(p_value = 2 * stats::pnorm(-abs(effect / standard_error)))
write_out(lag_effects, "lag_effects")

# (a) Temperature effect by long-run temperature.
test_a <- NULL
if (!NO_LEVELS) {
  temperature_var <- CLIMATE_VARS[["TM"]]
  test_a <- lag_effects %>%
    filter(expression == temperature_var, horizon == "Cumulative through lag",
           lag == MAX_LAG, at %in% c("low", "high")) %>%
    select(id, at, at_value, effect, standard_error, p_value) %>%
    tidyr::pivot_wider(names_from = at,
                       values_from = c(at_value, effect, standard_error, p_value)) %>%
    left_join(
      bind_rows(lapply(specs$id, function(id) {
        # Slope of the temperature effect in temperature: the (summed)
        # interaction with mean_TM_all, or 2 * b2 for the quadratic.
        slope_weights <- if (QUADRATIC) {
          stats::setNames(2, coef_names("TM^2"))
        } else {
          sum_weights(paste0(coef_names(temperature_var), ":mean_TM_all"))
        }
        lincom(fits_full[[id]], slope_weights, label = "slope") %>%
          transmute(id = id, slope = estimate, slope_se = std_error,
                    slope_p = p_value)
      })),
      by = "id"
    ) %>%
    mutate(holds = effect_low > 0 & effect_high < 0,
           holds_significant = holds & p_value_low < 0.05 & p_value_high < 0.05) %>%
    left_join(specs, by = "id")
  write_out(test_a, "test_a_temperature_by_climate")
}

# (b) Effect of the squared anomalies (summed over the lags).
test_b <- bind_rows(lapply(specs$id, function(id) {
  bind_rows(lapply(unlist(SQUARED_TERMS), function(term) {
    lincom(fits_full[[id]], sum_weights(coef_names(term)), term)
  })) %>%
    mutate(id = id, .before = 1)
})) %>%
  rename(term = test) %>%
  mutate(holds = estimate < 0, holds_significant = holds & p_value < 0.05) %>%
  left_join(specs, by = "id")
if (ANOMALY_X) {
  # With the climate interactions, the effect of each anomaly at each
  # long-run climate (b + g * m), its slope in the climate (g), and the
  # climate where the effect changes sign (-b / g, delta-method SE).
  anomaly_slopes <- bind_rows(lapply(specs$id, function(id) {
    full <- fits_full[[id]]
    terms <- names(stats::coef(full))
    bind_rows(lapply(unlist(SQUARED_TERMS), function(term) {
      b_name <- coef_names(term)
      g_name <- interaction_name(
        b_name, ANOMALY_INTERACTIONS[[squared_variable(term)]], terms
      )
      cf <- stats::coef(full)[c(b_name, g_name)]
      V <- stats::vcov(full)[c(b_name, g_name), c(b_name, g_name)]
      grad <- c(-1 / cf[[2]], cf[[1]] / cf[[2]]^2)
      tibble::tibble(
        id = id, term = term,
        slope = cf[[2]], slope_se = sqrt(V[2, 2]),
        slope_p = 2 * stats::pnorm(-abs(cf[[2]] / sqrt(V[2, 2]))),
        sign_change = -cf[[1]] / cf[[2]],
        sign_change_se = sqrt(drop(t(grad) %*% V %*% grad))
      )
    }))
  }))
  test_b <- lag_effects %>%
    filter(expression %in% unlist(SQUARED_TERMS),
           horizon == "Cumulative through lag") %>%
    transmute(id, term = as.character(expression), at, at_value,
              estimate = effect, std_error = standard_error, p_value) %>%
    left_join(anomaly_slopes, by = c("id", "term")) %>%
    mutate(holds = estimate < 0, holds_significant = holds & p_value < 0.05) %>%
    left_join(specs, by = "id")
}
write_out(test_b, "test_b_squared_anomalies")

# (c) Information added by the anomalies.
test_c <- bind_rows(lapply(specs$id, function(id) {
  full <- fits_full[[id]]
  restricted <- fits_restricted[[id]]
  squared_names <- lapply(SQUARED_TERMS, function(x) unlist(lapply(x, coef_names)))
  if (ANOMALY_X) {
    # The anomaly terms and their climate interactions.
    squared_names <- lapply(SQUARED_TERMS, function(x) {
      grep(paste(gsub("^2", "", x, fixed = TRUE), collapse = "|"),
           names(stats::coef(full)), value = TRUE)
    })
  }
  joint <- wald_row(full, exact_keep(unlist(squared_names)), label = "all")
  tibble::tibble(
    id = id,
    wald_F = joint$F_stat, wald_df1 = joint$df1, wald_p = joint$p_value,
    wald_p_zTM = wald_row(full, exact_keep(squared_names$TM))$p_value,
    wald_p_zRR = wald_row(full, exact_keep(squared_names$RR))$p_value,
    within_r2_full = unname(fixest::fitstat(full, "wr2")[[1]]),
    within_r2_restricted = unname(fixest::fitstat(restricted, "wr2")[[1]]),
    delta_aic = stats::AIC(full) - stats::AIC(restricted),
    delta_bic = stats::BIC(full) - stats::BIC(restricted)
  )
})) %>%
  mutate(holds = wald_p < 0.05) %>%
  left_join(specs, by = "id")
write_out(test_c, "test_c_information")

# "anomx": anomalies-only model against the levels-only BHM (test_c, where
# negative delta_aic / delta_bic favour the anomalies) and against the
# combined model (levels_add_p, anomalies_add_p, fit differences).
test_levels <- NULL
if (NO_LEVELS) {
  test_levels <- bind_rows(level_tests) %>%
    left_join(test_c %>% select(id, wald_p, within_r2_anomalies = within_r2_full,
                                within_r2_levels = within_r2_restricted,
                                delta_aic_vs_levels = delta_aic,
                                delta_bic_vs_levels = delta_bic), by = "id") %>%
    left_join(specs, by = "id")
  write_out(test_levels, "test_levels_vs_anomalies")
}

# "bhmsplit" tests. Symmetry: warm and cold (wet and dry) anomalies have equal
# squared coefficients, per variable and jointly. Average curvature of the
# precipitation anomaly: the mean of its two coefficients, which under
# symmetry is the zRR^2 coefficient. Hot versus cold units: see hot_tests.
test_split <- NULL
if (SPLIT) {
  test_split <- bind_rows(lapply(specs$id, function(id) {
    full <- fits_full[[id]]
    terms <- names(stats::coef(full))
    row_TM <- equality_row("I(abs_zTMp^2)", "I(abs_zTMm^2)", terms)
    row_RR <- equality_row("I(abs_zRRp^2)", "I(abs_zRRm^2)", terms)
    rr_mean <- lincom(full, c(`I(abs_zRRp^2)` = 0.5, `I(abs_zRRm^2)` = 0.5))
    tibble::tibble(
      id = id,
      symmetry_p = wald_linear(full, rbind(row_TM, row_RR)),
      symmetry_p_TM = wald_linear(full, rbind(row_TM)),
      symmetry_p_RR = wald_linear(full, rbind(row_RR)),
      zRR2_mean = rr_mean$estimate, zRR2_mean_p = rr_mean$p_value
    )
  })) %>%
    left_join(bind_rows(hot_tests), by = "id") %>%
    left_join(specs, by = "id")
  write_out(test_split, "test_split_symmetry_hot_cold")
}

# Peak cumulative effect of every lagged term (each climate variable at each
# long-run climate, zTM^2, zRR^2): the cumulative effect through the lag (0-10) with
# the largest absolute value, sign kept. It shows effects that build up and
# then fade, which the lag-10 total misses. Its interval is the pointwise one
# at that lag and does not account for choosing the lag. Not defined for the
# contemporaneous model.
peak_effects <- NULL
if (MAX_LAG > 0L) {
  peak_effects <- lag_effects %>%
    filter(horizon == "Cumulative through lag") %>%
    mutate(expression = as.character(expression)) %>%
    group_by(id, expression, at) %>%
    slice_max(abs(effect), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(id, expression, at, peak_lag = lag, estimate = effect,
              std_error = standard_error, p_value) %>%
    left_join(specs, by = "id")
  write_out(peak_effects, "peak_cumulative_effects")
}

summary_table <- specs %>%
  select(id, econ, climate, weight, deviation)
if (!is.null(test_a)) {
  summary_table <- summary_table %>%
    left_join(test_a %>% select(id, a_holds = holds,
                                a_significant = holds_significant), by = "id")
}
summary_table <- summary_table %>%
  left_join(test_b %>%
              group_by(id) %>%
              summarise(b_holds = all(holds),
                        b_significant = all(holds_significant)), by = "id") %>%
  left_join(test_c %>% select(id, c_holds = holds, delta_bic), by = "id")
write_out(summary_table, "summary")

# Model scores: the number of these conditions each model meets, all at 5%:
# a significant temperature gain at the low and loss at the high long-run
# temperature (a), significantly negative zTM^2 and zRR^2 (b), a rejected
# joint Wald test on the anomalies and a lower BIC with them (c). Models are
# ranked by score, then by the Wald p-value.
# The scores are not defined for "bhmsplitx", whose anomaly effects vary
# with the climate.
model_scores <- NULL
if (!ANOMALY_X) {
  model_scores <- test_a %>%
    select(id, effect_low, p_value_low, effect_high, p_value_high) %>%
    left_join(
      test_b %>%
        transmute(id, term = sub("^", "", term, fixed = TRUE), estimate, p_value) %>%
        tidyr::pivot_wider(names_from = term, values_from = c(estimate, p_value)),
      by = "id"
    ) %>%
    left_join(test_c %>% select(id, wald_p, delta_aic, delta_bic), by = "id") %>%
    mutate(
      a_low = effect_low > 0 & p_value_low < 0.05,
      a_high = effect_high < 0 & p_value_high < 0.05,
      c_wald = wald_p < 0.05,
      c_bic = delta_bic < 0
    )
  if (SPLIT) {
    # Under "bhmsplit" the anomaly conditions are: warm temperature anomalies
    # significantly harmful (zTMp^2), precipitation anomalies harmful on average
    # (mean of the two zRR^2 coefficients), warm/cold and wet/dry symmetry not
    # rejected, and equal responses in hot and cold units not rejected.
    model_scores <- model_scores %>%
      left_join(test_split %>%
                  select(id, symmetry_p, symmetry_p_TM, symmetry_p_RR, zRR2_mean,
                         zRR2_mean_p, hot_p, hot_p_TM, hot_p_RR, zTMp2_cold,
                         zTMp2_hot), by = "id") %>%
      mutate(
        b_TMp = estimate_abs_zTMp2 < 0 & p_value_abs_zTMp2 < 0.05,
        b_RR = zRR2_mean < 0 & zRR2_mean_p < 0.05,
        symmetric = symmetry_p >= 0.05,
        hot_cold_equal = hot_p >= 0.05,
        score = a_low + a_high + b_TMp + b_RR + c_wald + c_bic + symmetric +
          hot_cold_equal
      )
  } else {
    model_scores <- model_scores %>%
      mutate(
        b_TM = estimate_zTM2 < 0 & p_value_zTM2 < 0.05,
        b_RR = estimate_zRR2 < 0 & p_value_zRR2 < 0.05,
        score = a_low + a_high + b_TM + b_RR + c_wald + c_bic
      )
  }
  model_scores <- model_scores %>%
    left_join(specs, by = "id") %>%
    arrange(desc(score), wald_p)
  write_out(model_scores, "model_scores")
}

saveRDS(
  list(specification = SPECIFICATION, specs = specs,
       sample_summary = sample_summary, lag_effects = lag_effects,
       test_a = test_a, test_b = test_b, test_c = test_c,
       peak_effects = peak_effects, summary = summary_table,
       test_split = test_split, test_levels = test_levels,
       model_scores = model_scores),
  file.path(OUTPUT_DIR, "analysis.rds")
)

options(width = 200)
round_cols <- function(d) mutate(d, across(where(is.double), ~ signif(.x, 3)))
cat("\n== Specification:", SPECIFICATION, "==\n")
cat("\n== Samples ==\n")
print(round_cols(sample_summary %>%
  select(econ, climate, weight, deviation, observations, units, countries, first_year,
         last_year, max_abs_zTM, max_abs_zRR)), n = Inf)
if (!is.null(test_a)) {
  cat("\n== (a)", EFFECT_LABEL, "of", temperature_var, "in cold (",
      AT_VALUES$mean_TM_all[["low"]], DEG_C, ") and hot (",
      AT_VALUES$mean_TM_all[["high"]], DEG_C, ") countries ==\n")
  print(round_cols(test_a %>%
    select(econ, climate, weight, deviation, effect_low, p_value_low, effect_high,
           p_value_high, slope, slope_p, holds, holds_significant)), n = Inf)
}
cat("\n== (b)", EFFECT_LABEL, "of the squared anomalies ==\n")
print(round_cols(test_b %>%
  select(econ, climate, weight, deviation, term, estimate, std_error, p_value, holds,
         holds_significant)), n = Inf)
cat("\n== (c) Information added by the anomalies ==\n")
print(round_cols(test_c %>%
  select(econ, climate, weight, deviation, wald_F, wald_p, wald_p_zTM, wald_p_zRR,
         within_r2_restricted, within_r2_full, delta_aic, delta_bic, holds)),
  n = Inf)
if (!is.null(peak_effects)) {
  cat("\n== Peak cumulative effect over lags 0-10 (largest absolute value) ==\n")
  print(round_cols(peak_effects %>%
    select(econ, climate, weight, deviation, expression, at, peak_lag, estimate,
           std_error, p_value)), n = Inf)
}
cat("\n== Summary ==\n")
print(round_cols(summary_table %>% select(-id)), n = Inf)

# Part 2: plots ----------------------------------------------------------------------

Y_LABEL <- "Effect on log GDP per-capita growth"
# Figure sizes grow with the number of economic datasets.
n_econ <- length(ECON_DATASETS)
n_models <- n_econ * length(CLIMATE_SOURCES) * length(CLIMATE_WEIGHTS)

# Lag paths, for the distributed-lag specification only.
if (MAX_LAG > 0L) {
  plot_data <- lag_effects %>%
    mutate(
      climate = factor(climate, levels = names(CLIMATE_SOURCES)),
      econ = factor(econ, levels = names(ECON_DATASETS)),
      at = factor(at, levels = c(AT_CODES, "NA"),
                  labels = c(unname(AT_LABELS), "NA"))
    )

  for (dev in names(DEVIATIONS)) for (wt in names(CLIMATE_WEIGHTS)) {
    d <- plot_data %>% filter(deviation == dev, weight == wt)
    slug <- gsub("[^a-z0-9]+", "_", tolower(paste(dev, wt)))
    label <- paste0(dev, ", ", CLIMATE_WEIGHTS[[wt]], " weights")

    p_diff <- plot_distributed_lags(
      d %>% filter(expression %in% CLIMATE_VARS),
      CONF_LEVEL, Y_LABEL, colour = "climate", colour_name = "Climate dataset",
      horizons = "Cumulative through lag",
      facets = expression + econ ~ at,
      title = paste0("First differences, cumulative effect (", label, ")"),
      subtitle = paste("Rows: variable and economic dataset; columns: long-run",
                       "mean temperature / precipitation per year at which the",
                       "effect is evaluated")
    )
    ggsave(file.path(OUTPUT_DIR, paste0("first_differences_", slug, ".png")),
           p_diff, width = 12, height = 2 + 4 * n_econ, dpi = 300)

    p_sq <- plot_distributed_lags(
      d %>% filter(expression %in% c("zTM^2", "zRR^2")),
      CONF_LEVEL, Y_LABEL, colour = "climate", colour_name = "Climate dataset",
      facets = expression ~ horizon + econ,
      title = paste0("Squared anomalies (", label, ")"),
      subtitle = "Columns: horizon and economic dataset"
    )
    ggsave(file.path(OUTPUT_DIR, paste0("squared_anomalies_", slug, ".png")),
           p_sq, width = 3 + 4.5 * n_econ, height = 7, dpi = 300)
  }
}

# Robustness summaries, one figure per variable: the effect of the climate
# variable at the long-run climates of AT_VALUES and of the squared anomaly,
# one colour and shape per anomaly definition. Under
# "lags" the cumulative effect through lag 10 and at its peak over lags 0-10;
# under "bhmdev" the contemporaneous effect.
critical <- stats::qnorm(1 - (1 - CONF_LEVEL) / 2)
STATISTICS <- if (MAX_LAG > 0L) {
  c("Through lag 10", "Peak over lags 0-10")
} else {
  "Contemporaneous"
}
summary_quantities <- function(v, values, unit) {
  x <- CLIMATE_VARS[[v]]
  z <- SQUARED_TERMS[[v]]
  # Under "bhmsplitx" each anomaly is shown at every long-run climate.
  z_labels <- if (ANOMALY_X) {
    as.vector(outer(paste0(values, " ", unit), z, function(a, b) paste(b, "at", a)))
  } else {
    z
  }
  z_keys <- if (ANOMALY_X) {
    as.vector(outer(names(values), z, function(a, b) paste(b, a)))
  } else {
    paste(z, "NA")
  }
  stats::setNames(
    c(paste0(x, " at ", values, " ", unit), z_labels),
    c(paste(x, names(values)), z_keys)
  )
}
SUMMARY_VARIABLES <- list(
  TM = summary_quantities("TM", AT_VALUES$mean_TM_all, DEG_C),
  RR = summary_quantities("RR", AT_VALUES$mean_RR_all, "m")
)
SUMMARY_TITLES <- c(TM = "temperature", RR = "precipitation")

forest <- lag_effects %>%
  filter(horizon == "Cumulative through lag", lag == MAX_LAG) %>%
  transmute(id, key = paste(expression, at), estimate = effect,
            std_error = standard_error, statistic = STATISTICS[[1]])
if (!is.null(peak_effects)) {
  forest <- bind_rows(
    forest,
    peak_effects %>%
      transmute(id, key = paste(expression, at), estimate, std_error,
                statistic = STATISTICS[[2]])
  )
}
forest <- forest %>%
  left_join(specs, by = "id") %>%
  mutate(
    statistic = factor(statistic, levels = STATISTICS),
    model = factor(paste(econ, climate, weight, sep = ", "),
                   levels = rev(unique(paste(specs$econ, specs$climate, specs$weight,
                                             sep = ", ")))),
    deviation = factor(deviation, levels = names(DEVIATIONS)),
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error
  )
deviation_levels <- levels(forest$deviation)

summary_subtitle <- if (MAX_LAG > 0L) {
  paste0("Rows: cumulative effect through lag 10, and at its largest absolute ",
         "value over lags 0-10 (sign kept); ", round(100 * CONF_LEVEL),
         "% clustered CI")
} else {
  paste0("Contemporaneous effect; ", round(100 * CONF_LEVEL), "% clustered CI")
}
summary_caption <- paste(
  if (QUADRATIC) {
    paste("Marginal effects of TM and RR (b1 + 2 * b2 * level) at a temperature",
          "and at a precipitation in metres per year.")
  } else {
    paste("Effects of TM / dTM evaluated at a long-run mean temperature",
          "(mean_TM_all) and of RR / dRR at a long-run precipitation in metres",
          "per year (mean_RR_all).")
  },
  "Anomalies scaled by the SD over the entire climate series (MA),",
  "of the same trailing window (SD),\nor of the Hamilton cycle (Hamilton",
  "trend), as set in DEVIATIONS.",
  if (DEVIATION_UNITS[["TM"]] == "absolute") {
    paste("zTM is not scaled: it is in degrees Celsius for every definition,",
          "which then differ for zTM only in the baseline mean.")
  },
  if (MAX_LAG > 0L) {
    "Peak intervals are pointwise at the peak lag and do not account for choosing it."
  }
)

for (v in names(SUMMARY_VARIABLES)) {
  quantities <- SUMMARY_VARIABLES[[v]]
  d <- forest %>%
    filter(key %in% names(quantities)) %>%
    mutate(quantity = factor(unname(quantities[key]), levels = unname(quantities)))

  p_forest <- ggplot(d, aes(x = estimate, y = model, colour = deviation,
                            shape = deviation)) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0,
                  orientation = "y", linewidth = 0.5,
                  position = position_dodge(width = 0.75)) +
    geom_point(size = 2, position = position_dodge(width = 0.75)) +
    facet_grid(statistic ~ quantity, scales = "free_x") +
    scale_colour_manual(values = grid_colours(deviation_levels),
                        name = "Anomaly definition") +
    scale_shape_manual(values = grid_shapes(deviation_levels),
                       name = "Anomaly definition") +
    labs(x = paste(EFFECT_LABEL, "on log GDP per-capita growth"), y = NULL,
         title = paste0("Robustness of the ", MODEL_LABEL, ": ",
                        SUMMARY_TITLES[[v]]),
         subtitle = summary_subtitle, caption = summary_caption) +
    theme_results() +
    theme(legend.position = "bottom")
  ggsave(file.path(OUTPUT_DIR, paste0("robustness_summary_", v, ".png")),
         p_forest, width = 4 * length(quantities),
         height = 2 + n_models * if (MAX_LAG > 0L) 1.5 else 0.75,
         dpi = 300)
}

message("Results written to ", normalizePath(OUTPUT_DIR))
