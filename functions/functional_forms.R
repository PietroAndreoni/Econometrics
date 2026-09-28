# Hinged anomaly terms and helpers for comparing parametric weather-response
# forms against the nonparametric signed bins.
#
# A hinge at threshold c is h_c = max(0, |z| - c): zero for anomalies within
# c of the baseline and linear beyond it, symmetric in the sign of z. The
# notebook's zzTM is the temperature hinge at 1.5 (SD, with standardized
# anomalies). A threshold of 0 gives |z| itself.

hinge_col <- function(threshold, var = "TM") {
  sprintf("h%s_%02d", var, round(threshold * 10))
}

# Adds h<var>_<10c> for each threshold c and, if `sided`, the one-sided arms
# h<var>_<10c>_cold = max(0, -z - c) and _hot = max(0, z - c), whose sum is the
# symmetric hinge. `z` is the anomaly column (zTM by default).
add_hinge_columns <- function(panel, thresholds, var = "TM", z = paste0("z", var),
                              sided = TRUE) {
  for (threshold in thresholds) {
    stem <- hinge_col(threshold, var)
    panel[[stem]] <- pmax(0, abs(panel[[z]]) - threshold)
    if (sided) {
      panel[[paste0(stem, "_cold")]] <- pmax(0, -panel[[z]] - threshold)
      panel[[paste0(stem, "_hot")]] <- pmax(0, panel[[z]] - threshold)
    }
  }
  panel
}

# fixest reads `x^k` in a formula as the single power I(x^k) and names the
# coefficient that way; writing I(x^k) would be wrapped twice. So formulas use
# `^` and coefficient names are derived separately.
power_formula <- function(x, k) ifelse(k == 1L, x, sprintf("%s^%d", x, k))
power_coef <- function(x, k) ifelse(k == 1L, x, sprintf("I(%s^%d)", x, k))

# Equality of the cold and hot arms of a sided hinge (see add_hinge_columns()),
# a t test on their difference with the model's degrees of freedom.
hinge_symmetry_test <- function(model, threshold, var = "TM") {
  stem <- hinge_col(threshold, var)
  terms <- paste0(stem, c("_cold", "_hot"))
  b <- stats::coef(model)[terms]
  test <- lincom(model, stats::setNames(c(1, -1), terms),
                 df = fixest::degrees_freedom(model, "t"))
  tibble::tibble(
    threshold = threshold,
    cold_slope = unname(b[[1]]),
    hot_slope = unname(b[[2]]),
    difference = test$estimate,
    se = test$std_error,
    statistic = test$statistic,
    p_value = test$p_value
  )
}

# Spread of the estimates of `terms` across models (a tidy_fixest()-style table
# with columns model, term, estimate, std_error, p_value), also in units of the
# median standard error, which is the scale that matters for inference.
coefficient_stability <- function(coefs, terms, scope) {
  coefs %>%
    filter(term %in% terms) %>%
    group_by(term) %>%
    summarise(
      scope = scope,
      n_models = n(),
      min = min(estimate), max = max(estimate), mean = mean(estimate),
      sd = stats::sd(estimate), median_se = stats::median(std_error),
      spread = max - min,
      spread_in_se = spread / median_se,
      sd_in_se = sd / median_se,
      same_sign = n_distinct(sign(estimate)) == 1L,
      all_p05 = all(p_value < 0.05),
      .groups = "drop"
    ) %>%
    relocate(scope, .after = term)
}
