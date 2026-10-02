# Linear or quadratic arms in the deviations-only model? With linear arms the
# total loss from a warming path depends only on the total warming (the
# anomalies against a trailing k-year mean of any path from 0 to dT sum to
# dT * (k + 1) / 2); with convex arms faster warming costs more. This script
# compares, on the identical sample, three forms of the warm, cold, wet and
# dry arms (abs_zTMp, abs_zTMm, abs_zRRp, abs_zRRm), each interacted with its
# long-run climate:
#
#   linear     a * z + a_m * z * m
#   quadratic  q * z^2 + q_m * z^2 * m
#   both       the linear and the quadratic terms together
#
# and reports fit (within R2, AIC, BIC; linear and quadratic have the same
# number of parameters), Wald tests in the combined model of the quadratic
# terms given the linear ones and the reverse (temperature and precipitation
# separately), and the effect of a 1- and 2-SD warm or cold year at long-run
# temperatures of 5, 15 and 25 degrees Celsius. In the combined model the
# ratio of the 2-SD to the 1-SD effect is 2 for a linear response and 4 for a
# quadratic one.
#
# Grid: WB and PWT x ERA5, CRU TS and UDelaware, area and population weights.
# Anomalies against the trailing 30-year mean, scaled by the SD over the entire
# climate series. Fixed effects as in main_analysis.Rmd, errors clustered by
# country.
#
# Output: results/deviation_arm_form/ (fit_tests.csv, arm_effects.csv,
# arm_curves_<weight>.png)

source("load_functions.R")

ECON_DATASETS <- c(WB = "WB", PWT = "PWT")
CLIMATE_SOURCES <- c(ERA5 = "ERA5", `CRU TS` = "CRU TS", UDelaware = "UDelaware")
CLIMATE_WEIGHTS <- c(area = "area", pop = "concurrent population")
BASELINE <- "lag_30"
SD_BASELINE <- "all"
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"
CONF_LEVEL <- 0.95
AT_TM <- c(5, 15, 25)
ANOMALIES <- seq(-2, 2, by = 0.05)
DEG_C <- "°C"

# Arm variables: linear (l) and squared (q), warm/cold and wet/dry.
ARMS <- c(TMp = "mean_TM_all", TMm = "mean_TM_all", RRp = "mean_RR_all", RRm = "mean_RR_all")
arm_terms <- function(prefix) {
  paste(c(paste0(prefix, names(ARMS)), paste0(prefix, names(ARMS), ":", ARMS)),
        collapse = " + ")
}
FORMS <- list(
  linear = arm_terms("l"),
  quadratic = arm_terms("q"),
  both = paste(arm_terms("l"), "+", arm_terms("q"))
)
FORM_LABELS <- c(linear = "Linear arms", quadratic = "Quadratic arms",
                 both = "Linear + quadratic arms")

OUTPUT_DIR <- file.path("results", "deviation_arm_form")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit <- function(terms, data) {
  fit_panel_model(terms, data, fixed_effects = FIXED_EFFECTS, cluster = ~GID_0)
}
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
# Wald p-value of the arms with the given prefix for one variable (TM or RR),
# main terms and interactions.
block_p <- function(model, prefix, variable) {
  terms <- names(stats::coef(model))
  keep <- grep(paste0("(^|:)", prefix, variable, "[pm]($|:)"), terms, value = TRUE)
  fixest::wald(model, keep = paste0("^(", paste(keep, collapse = "|"), ")$"),
               print = FALSE)$p
}

# Effect of a year d SD above (d > 0) or below (d < 0) normal temperature at
# long-run temperature m, for each d in `anomalies`.
arm_curve <- function(model, anomalies, m) {
  terms <- names(stats::coef(model))
  G <- t(vapply(anomalies, function(d) {
    arm <- if (d >= 0) "TMp" else "TMm"
    w <- stats::setNames(numeric(length(terms)), terms)
    powers <- c(l = 1, q = 2)
    for (prefix in names(powers)) {
      x <- paste0(prefix, arm)
      p <- powers[[prefix]]
      if (x %in% terms) {
        w[[x]] <- abs(d)^p
        w[[interaction_name(x, "mean_TM_all", terms)]] <- abs(d)^p * m
      }
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(anomaly = anomalies)
}

tests <- list()
curves <- list()
for (w in names(CLIMATE_WEIGHTS)) for (econ in names(ECON_DATASETS)) {
  for (source_name in names(CLIMATE_SOURCES)) {
    message("Fitting ", econ, " x ", source_name, ", ", w, " weights")
    cfg <- panel_config(
      econ_year_min = min(load_econ_panel(resolve_econ_source(ECON_DATASETS[[econ]]))$year),
      climate_source = CLIMATE_SOURCES[[source_name]],
      climate_weight = CLIMATE_WEIGHTS[[w]],
      clim_history_years = 1L
    )
    d <- suppressMessages(build_dat(ECON_DATASETS[[econ]], config = cfg,
                                    baseline = BASELINE, sd_baseline = SD_BASELINE)) %>%
      mutate(lTMp = abs_zTMp, lTMm = abs_zTMm, lRRp = abs_zRRp, lRRm = abs_zRRm,
             qTMp = abs_zTMp^2, qTMm = abs_zTMm^2, qRRp = abs_zRRp^2, qRRm = abs_zRRm^2)

    models <- lapply(FORMS, fit, data = d)
    stopifnot(nobs(models$linear) == nobs(models$both),
              nobs(models$quadratic) == nobs(models$both))
    both <- models$both
    tests[[length(tests) + 1]] <- tibble::tibble(
      weight = w, econ = econ, climate_source = source_name, observations = nobs(both),
      wr2_linear = fixest::fitstat(models$linear, "wr2")[[1]],
      wr2_quadratic = fixest::fitstat(models$quadratic, "wr2")[[1]],
      wr2_both = fixest::fitstat(both, "wr2")[[1]],
      aic_quadratic_minus_linear = stats::AIC(models$quadratic) - stats::AIC(models$linear),
      aic_both_minus_linear = stats::AIC(both) - stats::AIC(models$linear),
      aic_both_minus_quadratic = stats::AIC(both) - stats::AIC(models$quadratic),
      bic_both_minus_linear = stats::BIC(both) - stats::BIC(models$linear),
      bic_both_minus_quadratic = stats::BIC(both) - stats::BIC(models$quadratic),
      TM_quadratic_given_linear_p = block_p(both, "q", "TM"),
      TM_linear_given_quadratic_p = block_p(both, "l", "TM"),
      RR_quadratic_given_linear_p = block_p(both, "q", "RR"),
      RR_linear_given_quadratic_p = block_p(both, "l", "RR")
    )
    for (form in names(models)) for (m in AT_TM) {
      curves[[length(curves) + 1]] <- arm_curve(models[[form]], ANOMALIES, m) %>%
        mutate(weight = w, econ = econ, climate_source = source_name, form = form,
               climate = m)
    }
  }
}

tests <- bind_rows(tests)
curves <- bind_rows(curves)
utils::write.csv(tests, file.path(OUTPUT_DIR, "fit_tests.csv"), row.names = FALSE)

# Effects of a 1- and 2-SD warm or cold year, and the ratio of the 2-SD to the
# 1-SD effect (2 if linear, 4 if quadratic).
effects <- curves %>%
  filter(abs(anomaly) %in% c(1, 2)) %>%
  mutate(direction = ifelse(anomaly > 0, "warm", "cold"), size = abs(anomaly)) %>%
  select(weight, econ, climate_source, form, climate, direction, size, estimate,
         std_error, conf_low, conf_high)
utils::write.csv(effects, file.path(OUTPUT_DIR, "arm_effects.csv"), row.names = FALSE)

options(width = 220)
cat("\n== Fit and Wald tests (combined model) ==\n")
print(as.data.frame(tests %>%
  select(weight, econ, climate_source, wr2_linear, wr2_quadratic, wr2_both,
         aic_quadratic_minus_linear, aic_both_minus_linear, aic_both_minus_quadratic,
         bic_both_minus_quadratic, ends_with("_p")) %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
cat("\n== Warm and cold years in the combined model: 1 and 2 SD, ratio ==\n")
print(as.data.frame(effects %>%
  filter(form == "both") %>%
  select(weight, econ, climate_source, climate, direction, size, estimate, std_error) %>%
  tidyr::pivot_wider(names_from = size, values_from = c(estimate, std_error)) %>%
  mutate(ratio = estimate_2 / estimate_1,
         across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

# Plots ---------------------------------------------------------------------------------

form_levels <- unname(FORM_LABELS)
colours <- grid_colours(form_levels)
for (w in names(CLIMATE_WEIGHTS)) {
  d <- curves %>%
    filter(weight == w) %>%
    mutate(
      form = factor(FORM_LABELS[form], levels = form_levels),
      dataset = factor(paste(econ, "×", climate_source),
                       levels = as.vector(outer(names(ECON_DATASETS), names(CLIMATE_SOURCES),
                                                paste, sep = " × "))),
      climate = factor(paste0("Long-run temperature ", climate, " ", DEG_C),
                       levels = paste0("Long-run temperature ", AT_TM, " ", DEG_C)),
      across(c(estimate, conf_low, conf_high), ~ 100 * .x)
    )
  p <- ggplot(d, aes(x = anomaly, y = estimate, colour = form, fill = form)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(dataset ~ climate, scales = "free_y") +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    labs(
      x = "Temperature anomaly in the year (SD)",
      y = "Effect on growth that year (pp)",
      title = paste0("Linear or quadratic arms in the deviations-only model (",
                     w, " weights)"),
      subtitle = paste("Warm and cold arms interacted with the long-run temperature,",
                       "evaluated at 5, 15 and 25 °C. Bands: 95% CI."),
      caption = paste("Anomalies against the trailing 30-year mean, scaled by the",
                      "full-record SD; precipitation arms in the same form; fixed effects",
                      "year + country + country trends; errors clustered by country.")
    ) +
    theme_results()
  ggsave(file.path(OUTPUT_DIR, paste0("arm_curves_", w, ".png")), p,
         width = 12, height = 16, dpi = 300)
}
message("Results written to ", normalizePath(OUTPUT_DIR))
