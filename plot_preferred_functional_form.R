# The temperature and precipitation response implied by the preferred
# specification against the nonparametric signed-anomaly bins.
#
# The preferred specification replaces the signed bins with a power of a hinged
# absolute temperature anomaly, h_c = max(0, |zTM| - c) (functions/
# functional_forms.R), plus a linear term in wet anomalies, abs_zRRp =
# max(zRR, 0). The notebook uses zzTM^2 = h_1.5^2.
#
# This script sweeps the hinge over 0, 0.5, 1.0 and 1.5 SD and takes the linear,
# quadratic and cubic power of each. A threshold of 0 is the unhinged |zTM|: a
# power of the raw signed anomaly would not do, because odd powers such as
# zTM^3 are antisymmetric and would make cold years beneficial whenever hot
# years are harmful. Every specification keeps the same precipitation term and
# the same estimation sample, so only the temperature response varies. Every
# hinge is symmetric in the anomaly, so the script also frees the cold and hot
# arms and tests their equality, and checks how stable the BHM base
# coefficients are across the weather forms.
#
# Specification: main_spec() (functions/analysis_spec.R), i.e. main_analysis.Rmd.
# Output: results/preferred_functional_form/.

source("load_functions.R")

SPEC <- main_spec()
HINGES <- c(0, 0.5, 1.0, 1.5)
NOTEBOOK_HINGE <- 1.5
NOTEBOOK_POWER <- 2L         # main_analysis.Rmd: zzTM^2 + abs_zRRp.
GRID_LIMIT <- 3.25           # Plotted anomaly range; matches the outer bin codes.
PP <- 100                    # Report in percentage points.
PRECIP_TERM <- "abs_zRRp"
BASE_TERMS <- c("TM", "TM:mean_TM_all", "RR", "RR:mean_RR_all")
OUT <- file.path("results", "preferred_functional_form")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# Panel ----------------------------------------------------------------------------

panel <- build_main_dat("DOSE", SPEC) %>%
  add_hinge_columns(HINGES) %>%
  mutate(zzTM = pmax(NOTEBOOK_HINGE, abs((TM - mean_TM_lag_30) / sd_TM_lag_30)) -
           NOTEBOOK_HINGE)
stopifnot(
  isTRUE(all.equal(panel$zzTM, panel[[hinge_col(NOTEBOOK_HINGE)]])),
  isTRUE(all.equal(panel$hTM_15_cold + panel$hTM_15_hot, panel$zzTM))
)

# One sample for every model, so the comparison is about functional form.
needed <- c(SPEC$outcome, "TM", "mean_TM_all", "RR", "mean_RR_all", "zTM", "zRR",
            PRECIP_TERM, "TM_bin_signed", "RR_bin_signed")
est <- panel %>% filter(if_all(all_of(needed), is.finite))
stopifnot(!anyDuplicated(est[c("GID_1", "year")]))
cat(sprintf("Estimation sample: %d region-years, %d regions, %d countries, %d-%d\n",
            nrow(est), n_distinct(est$GID_1), n_distinct(est$GID_0),
            min(est$year), max(est$year)))

fit <- function(terms) fit_main(paste(SPEC$base_climate, "+", terms), est, SPEC)

# Models -----------------------------------------------------------------------------

DEGREES <- c(`1` = "Linear", `2` = "Quadratic", `3` = "Cubic")
threshold_label <- function(threshold) {
  ifelse(threshold == 0, "No hinge: |zTM|", sprintf("Hinge at %.1f SD", threshold))
}
specs <- tidyr::expand_grid(threshold = HINGES, degree = as.integer(names(DEGREES))) %>%
  mutate(
    variable_name = hinge_col(threshold),
    id = paste0("hinge_", threshold, "_", degree),
    term_formula = power_formula(variable_name, degree),
    term_coef = power_coef(variable_name, degree),
    degree_label = factor(DEGREES[as.character(degree)], levels = unname(DEGREES)),
    panel_label = factor(threshold_label(threshold), levels = threshold_label(HINGES)),
    specification = paste0(panel_label, ", ", tolower(degree_label))
  )

models <- stats::setNames(
  lapply(paste(specs$term_formula, "+", PRECIP_TERM), fit),
  specs$id
)
m_bins <- fit("i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)")
m_base <- fit("zTM^2 + zRR^2")   # The notebook's quadratic benchmark.

# The notebook's preferred specification, verbatim, must coincide with its hinge.
m_pref <- fit(paste0("zzTM^", NOTEBOOK_POWER, " + abs_zRRp"))
stopifnot(isTRUE(all.equal(
  unname(stats::coef(m_pref)),
  unname(stats::coef(models[[paste0("hinge_", NOTEBOOK_HINGE, "_", NOTEBOOK_POWER)]]))
)))

# Hinges with free cold and hot arms.
asym_models <- stats::setNames(lapply(HINGES, function(threshold) {
  stem <- hinge_col(threshold)
  fit(sprintf("%s_cold + %s_hot + %s", stem, stem, PRECIP_TERM))
}), sprintf("%.1f", HINGES))

# Implied response curves -------------------------------------------------------------
#
# Each curve is re-centred on its average over the reference-bin observations
# (-0.5 to 0.5 SD), so it is comparable to the bin coefficients.

grid <- seq(-GRID_LIMIT, GRID_LIMIT, length.out = 601)
tm_ref_z <- est$zTM[est$TM_bin_signed == 0]
rr_ref_z <- est$zRR[est$RR_bin_signed == 0]
temp_basis <- function(z, threshold, degree) {
  matrix(pmax(0, abs(z) - threshold)^degree, ncol = 1)
}
rr_basis <- function(z) matrix(pmax(z, 0), ncol = 1)

curves <- bind_rows(lapply(seq_len(nrow(specs)), function(i) {
  s <- specs[i, ]
  model <- models[[s$id]]
  bind_rows(
    tibble::tibble(variable = "Temperature", anomaly = grid) %>%
      bind_cols(response_curve(model, s$term_coef,
                               temp_basis(grid, s$threshold, s$degree),
                               temp_basis(tm_ref_z, s$threshold, s$degree))),
    tibble::tibble(variable = "Precipitation", anomaly = grid) %>%
      bind_cols(response_curve(model, PRECIP_TERM, rr_basis(grid),
                               rr_basis(rr_ref_z)))
  ) %>%
    mutate(id = s$id, degree_label = s$degree_label, panel_label = s$panel_label)
}))

asym_curves <- bind_rows(lapply(HINGES, function(threshold) {
  stem <- hinge_col(threshold)
  arms <- function(z) cbind(pmax(0, -z - threshold), pmax(0, z - threshold))
  tibble::tibble(panel_label = factor(threshold_label(threshold),
                                      levels = levels(specs$panel_label)),
                 anomaly = grid) %>%
    bind_cols(response_curve(asym_models[[sprintf("%.1f", threshold)]],
                             paste0(stem, c("_cold", "_hot")),
                             arms(grid), arms(tm_ref_z)))
}))

bins <- bin_coefs(m_bins) %>%
  transmute(
    variable = if_else(bin_var == "TM_bin_signed", "Temperature", "Precipitation"),
    anomaly = as.numeric(bin_label),
    estimate, std_error, p_value, ci_low, ci_high,
    # The outermost codes stand for the open-ended tails beyond +-2.5 SD.
    open_ended = abs(anomaly) > 2.5
  )

# Tests --------------------------------------------------------------------------------

# The parametric forms are not nested in the bins, but the bins are nested in an
# encompassing model, so a joint Wald test on the bin indicators asks whether
# the nonparametric shape adds anything to the parametric response.
encompass <- function(label, terms) {
  bind_rows(
    wald_row(fit(paste(terms, "+ i(TM_bin_signed, ref = 0)")), "TM_bin_signed",
             label = label),
    wald_row(fit(paste(terms, "+ i(RR_bin_signed, ref = 0)")), "RR_bin_signed",
             label = label)
  )
}
encompassing_tests <- bind_rows(
  lapply(seq_len(nrow(specs)), function(i) {
    encompass(specs$specification[i], paste(specs$term_formula[i], "+", PRECIP_TERM))
  }),
  lapply(HINGES, function(threshold) {
    stem <- hinge_col(threshold)
    encompass(paste0(threshold_label(threshold), ", cold and hot arms free"),
              sprintf("%s_cold + %s_hot + %s", stem, stem, PRECIP_TERM))
  })
)
symmetry_tests <- bind_rows(lapply(HINGES, function(threshold) {
  hinge_symmetry_test(asym_models[[sprintf("%.1f", threshold)]], threshold)
}))

fit_row <- function(m, id) {
  tibble::tibble(
    model = id, n = stats::nobs(m), k = length(stats::coef(m)),
    within_r2 = as.numeric(fixest::fitstat(m, "wr2")[[1]]),
    r2 = as.numeric(fixest::fitstat(m, "r2")[[1]]),
    aic = stats::AIC(m), bic = stats::BIC(m)
  )
}
all_models <- c(
  models,
  stats::setNames(asym_models, paste0("asym_hinge_", names(asym_models))),
  list(signed_bins = m_bins, notebook_quadratic_benchmark = m_base)
)
model_fit <- bind_rows(Map(fit_row, all_models, names(all_models))) %>%
  left_join(specs %>% select(model = id, specification), by = "model") %>%
  arrange(desc(within_r2))

coef_table <- bind_rows(Map(tidy_fixest, all_models, names(all_models)))

# Stability of the BHM base coefficients across weather forms. The bins model
# replaces abs_zRRp with precipitation bins, so it is not like-for-like for the
# two RR controls.
base_coefs <- coef_table %>%
  filter(term %in% BASE_TERMS) %>%
  mutate(term = factor(term, levels = BASE_TERMS),
         comparable = model != "signed_bins")
base_stability <- bind_rows(
  coefficient_stability(filter(base_coefs, model %in% specs$id), BASE_TERMS,
                        "12 hinge powers"),
  coefficient_stability(filter(base_coefs, comparable), BASE_TERMS,
                        "all but signed bins"),
  coefficient_stability(base_coefs, BASE_TERMS, "all models")
) %>%
  arrange(term, scope)

# Do the base estimates track how well the weather term fits? A rank
# correlation near +-1 means the base estimate depends on the functional form.
base_vs_fit <- base_coefs %>%
  filter(model %in% specs$id) %>%
  left_join(model_fit %>% select(model, within_r2), by = "model") %>%
  group_by(term) %>%
  summarise(
    n_models = n(),
    spearman_estimate_vs_within_r2 = stats::cor(estimate, within_r2,
                                                method = "spearman"),
    estimate_at_best_fit = estimate[which.max(within_r2)],
    estimate_at_worst_fit = estimate[which.min(within_r2)],
    .groups = "drop"
  )

write_out(bins, "signed_bin_coefficients")
write_out(curves %>% select(-degree_label, -panel_label), "implied_response_curves")
write_out(model_fit, "model_fit_comparison")
write_out(encompassing_tests, "encompassing_tests")
write_out(symmetry_tests, "hinge_symmetry_tests")
write_out(coef_table, "coefficients")
write_out(base_coefs %>% select(-comparable), "base_coefficients")
write_out(base_stability, "base_coefficient_stability")
write_out(base_vs_fit, "base_coefficients_vs_fit")
writeLines(
  capture.output(fixest::etable(c(models, list(bins = m_bins)),
                                fitstat = ~n + r2 + wr2 + aic + bic)),
  file.path(OUT, "model_comparison.txt")
)

# Figures ------------------------------------------------------------------------------

variable_levels <- c("Temperature", "Precipitation")
panel_levels <- levels(specs$panel_label)
degree_colours <- stats::setNames(GRID_PALETTE[1:3], unname(DEGREES))
as_facets <- function(x) {
  if ("variable" %in% names(x)) x$variable <- factor(x$variable, levels = variable_levels)
  if ("panel_label" %in% names(x)) {
    x$panel_label <- factor(as.character(x$panel_label), levels = panel_levels)
  }
  x
}
# A shared y-range from the bins keeps the variables comparable; steep powers
# leave the frame in the tails, which is part of what the figure shows.
y_limits <- range(c(bins$ci_low, bins$ci_high), na.rm = TRUE) * PP * 1.05
hinge_marks <- specs %>%
  filter(threshold > 0) %>%
  distinct(panel_label, threshold) %>%
  mutate(panel_label = as.character(panel_label)) %>%
  tidyr::crossing(sign = c(-1, 1)) %>%
  mutate(xintercept = sign * threshold)

bin_layers <- function(bin_data, degree_scales = TRUE) {
  c(list(
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3),
    geom_errorbar(
      data = bin_data, aes(x = anomaly, ymin = ci_low * PP, ymax = ci_high * PP),
      width = 0.08, linewidth = 0.4, colour = GRID_INK, inherit.aes = FALSE
    ),
    geom_point(
      data = bin_data, aes(x = anomaly, y = estimate * PP, shape = open_ended),
      colour = GRID_INK, size = 1.9, fill = "white", inherit.aes = FALSE
    ),
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 21),
      labels = c(`FALSE` = "Interior bin", `TRUE` = "Open-ended tail bin"),
      name = "Nonparametric bins"
    ),
    coord_cartesian(ylim = y_limits),
    labs(x = "Standardized anomaly, within-region standard deviations",
         y = "Effect on growth (pp)"),
    theme_results(),
    theme(panel.spacing = unit(0.9, "lines"))
  ), if (degree_scales) list(
    scale_colour_manual(values = degree_colours, name = "Power of the anomaly"),
    scale_fill_manual(values = degree_colours, name = "Power of the anomaly",
                      guide = "none")
  ))
}

# Figure 1: every hinge threshold and power against the bins.
p_main <- ggplot() +
  bin_layers(as_facets(tidyr::crossing(bins, panel_label = panel_levels))) +
  geom_vline(data = as_facets(mutate(hinge_marks, variable = "Temperature")),
             aes(xintercept = xintercept), linetype = "dotted",
             colour = GRID_MUTED, linewidth = 0.35) +
  geom_line(data = as_facets(curves),
            aes(x = anomaly, y = estimate * PP, colour = degree_label),
            linewidth = 0.8) +
  facet_grid(panel_label ~ variable) +
  labs(
    title = "Preferred specification against nonparametric signed-anomaly bins",
    subtitle = paste(
      "Same estimation sample and climate controls throughout; bins against the",
      "-0.5 to 0.5 SD reference.\nDotted lines mark each hinge, below which the",
      "temperature term is zero. Only the temperature term varies, so the",
      "precipitation curves coincide."
    )
  )
ggsave(file.path(OUT, "signed_bins_vs_functional_form.png"), p_main,
       width = 11, height = 12, dpi = 200)

# Figure 2: temperature only, one panel per threshold and power, with bands.
temperature_curves <- as_facets(curves) %>% filter(variable == "Temperature")
p_bands <- ggplot() +
  bin_layers(as_facets(tidyr::crossing(
    bins %>% filter(variable == "Temperature") %>% select(-variable),
    specs %>% distinct(panel_label, degree_label) %>%
      mutate(panel_label = as.character(panel_label))
  ))) +
  geom_vline(data = as_facets(hinge_marks), aes(xintercept = xintercept),
             linetype = "dotted", colour = GRID_MUTED, linewidth = 0.35) +
  geom_ribbon(data = temperature_curves,
              aes(x = anomaly, ymin = conf_low * PP, ymax = conf_high * PP,
                  fill = degree_label), alpha = 0.18) +
  geom_line(data = temperature_curves,
            aes(x = anomaly, y = estimate * PP, colour = degree_label),
            linewidth = 0.8) +
  facet_grid(panel_label ~ degree_label) +
  labs(
    title = "Temperature response with 95% confidence bands",
    subtitle = paste("Delta-method intervals, clustered by region. Every term is a",
                     "power of a non-negative hinge, so each response is symmetric",
                     "in the anomaly and keeps one sign."),
    x = "Standardized temperature anomaly, within-region standard deviations"
  )
ggsave(file.path(OUT, "temperature_bands_by_hinge.png"), p_bands,
       width = 11, height = 11, dpi = 200)

# Figure 3: is symmetry, not the power, what the bins reject?
symmetric_linear <- curves %>%
  filter(variable == "Temperature", degree_label == "Linear") %>%
  select(panel_label, anomaly, estimate)
arm_colours <- c(`Symmetric in |anomaly|` = GRID_PALETTE[1],
                 `Cold and hot arms free` = GRID_PALETTE[5])
p_asym <- ggplot() +
  bin_layers(as_facets(tidyr::crossing(
    bins %>% filter(variable == "Temperature") %>% select(-variable),
    panel_label = panel_levels
  )), degree_scales = FALSE) +
  geom_vline(data = as_facets(hinge_marks), aes(xintercept = xintercept),
             linetype = "dotted", colour = GRID_MUTED, linewidth = 0.35) +
  geom_ribbon(data = asym_curves,
              aes(x = anomaly, ymin = conf_low * PP, ymax = conf_high * PP),
              fill = GRID_PALETTE[5], alpha = 0.18) +
  geom_line(data = symmetric_linear,
            aes(x = anomaly, y = estimate * PP, colour = "Symmetric in |anomaly|"),
            linewidth = 0.8) +
  geom_line(data = asym_curves,
            aes(x = anomaly, y = estimate * PP, colour = "Cold and hot arms free"),
            linewidth = 0.8) +
  scale_colour_manual(values = arm_colours, name = "Linear hinge") +
  facet_wrap(~panel_label, nrow = 1) +
  labs(
    title = "Diagnostic: is symmetry, not the power, what the bins reject?",
    subtitle = paste0(
      "Cold vs hot slope equality: p = ",
      paste(sprintf("%.2g at %.1f SD", symmetry_tests$p_value,
                    symmetry_tests$threshold), collapse = ", "), "."
    ),
    x = "Standardized temperature anomaly, within-region standard deviations"
  )
ggsave(file.path(OUT, "asymmetry_diagnostic.png"), p_asym,
       width = 15, height = 5.5, dpi = 200)

# Figure 4: base climate controls across weather specifications.
p_base <- ggplot(
  base_coefs %>% mutate(model = factor(model, levels = rev(sort(unique(model))))),
  aes(x = estimate * PP, y = model, colour = comparable)
) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf_low * PP, xmax = conf_high * PP),
                orientation = "y", width = 0.25, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_colour_manual(
    values = c(`TRUE` = GRID_PALETTE[1], `FALSE` = GRID_PALETTE[2]),
    labels = c(`TRUE` = "Same precipitation term",
               `FALSE` = "Precipitation bins instead"),
    name = NULL
  ) +
  facet_wrap(~term, scales = "free_x", nrow = 1) +
  labs(title = "Base climate controls across weather specifications",
       subtitle = "95% intervals clustered by region.",
       x = "Coefficient (percentage points)", y = NULL) +
  theme_results() +
  theme(panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "base_coefficient_stability.png"), p_base,
       width = 13, height = 6, dpi = 200)

# Console summary ------------------------------------------------------------------------

options(width = 180)
cat("\n=== Fit comparison, best within-R2 first ===\n")
print(model_fit, n = Inf)
cat("\n=== Do the bins add to the parametric form? ===\n")
print(encompassing_tests, n = Inf)
cat("\n=== Is each hinge symmetric between cold and hot tails? ===\n")
print(symmetry_tests)
cat("\n=== Base climate controls: spread across specifications ===\n")
print(base_stability, n = Inf)
cat("\n=== Base controls vs fit of the weather term (12 hinge powers) ===\n")
print(base_vs_fit)
cat("\nWrote figures and tables to ", OUT, "\n", sep = "")
