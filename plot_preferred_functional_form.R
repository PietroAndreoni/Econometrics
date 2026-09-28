# Overlay the temperature and precipitation response implied by the preferred
# specification on the nonparametric signed-anomaly bins.
#
# The preferred specification replaces the signed bins with a power of a hinged
# absolute standardized temperature anomaly,
#
#   h_c = max(0, |zTM| - c),   zTM = (TM - mean_TM_lag_30) / sd_TM_lag_30,
#
# plus a linear term in positive (wet) precipitation anomalies,
# abs_zRRp = max(zRR, 0). The notebook sets c = 1.5, where h_1.5 is `zzTM`.
#
# This script sweeps the hinge threshold over 0, 0.5, 1.0 and 1.5 standard
# deviations and takes the linear, quadratic and cubic power of each. A
# threshold of zero leaves the anomaly unhinged at |zTM|, which is the natural
# unhinged benchmark: a power of the raw signed anomaly would not do, because
# odd powers such as zTM^3 are antisymmetric and would make cold years
# beneficial whenever hot years are harmful. Every specification keeps the same
# precipitation term and the same estimation sample, so only the temperature
# response varies.
#
# Note on `^` inside fixest formulas: fixest translates `zzTM^3` into
# `I(zzTM^3)`, a single cubic term with no lower-order terms, rather than the
# formula-crossing operator that base R would apply. Every specification here is
# therefore a single power term, matching how the notebook writes it.
#
# Prerequisite: the notebook's parquet inputs under data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
})

source("rmd_chunks.R")

out <- "results/preferred_functional_form"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

# Hinge thresholds in standard deviations; 0 leaves the anomaly unhinged at
# |zTM|.
HINGES <- c(0, 0.5, 1.0, 1.5)
NOTEBOOK_HINGE <- 1.5        # The threshold the notebook's `zzTM` uses.
GRID_LIMIT <- 3.25           # Plotted anomaly range; matches the outer bin codes.
PP <- 100                    # Report coefficients in log-growth percentage points.
PRECIP_TERM <- "abs_zRRp"

# ---------------------------------------------------------------- panel -------

nb <- load_notebook_env("test_functions.Rmd")

cat("Building the DOSE panel...\n")
panel <- nb$build_dat(econ_data = "DOSE")

# One column per hinge threshold, plus the same hinge split by the sign of the
# anomaly for the symmetry test below.
hinge_col <- function(threshold) sprintf("hTM_%02d", round(threshold * 10))
for (threshold in HINGES) {
  stem <- hinge_col(threshold)
  panel[[stem]] <- pmax(0, abs(panel$zTM) - threshold)
  panel[[paste0(stem, "_cold")]] <- pmax(0, -panel$zTM - threshold)
  panel[[paste0(stem, "_hot")]] <- pmax(0, panel$zTM - threshold)
}

# The notebook's `zzTM` is the 1.5 SD hinge, and it is non-negative by
# construction, so abs(zzTM) is the identity.
panel <- panel %>%
  mutate(zzTM = pmax(NOTEBOOK_HINGE, abs((TM - mean_TM_lag_30) / sd_TM_lag_30)) -
           NOTEBOOK_HINGE)
stopifnot(
  all(panel$zzTM >= 0, na.rm = TRUE),
  isTRUE(all.equal(panel$zzTM, panel[[hinge_col(NOTEBOOK_HINGE)]])),
  isTRUE(all.equal(panel$hTM_15_cold + panel$hTM_15_hot, panel$zzTM))
)

# Every model must share one sample so that the overlay and the fit statistics
# compare functional forms rather than sample coverage.
needed <- c(
  "dlgrp_pc_usd", "TM", "mean_TM_all", "RR", "mean_RR_all",
  "zTM", "zRR", PRECIP_TERM, "TM_bin_signed", "RR_bin_signed"
)
est <- panel %>% filter(if_all(all_of(needed), is.finite))
stopifnot(!anyDuplicated(est[c("GID_1", "year")]))

cat(sprintf(
  "Estimation sample: %d region-years, %d regions, %d countries, %d-%d\n",
  nrow(est), n_distinct(est$GID_1), n_distinct(est$GID_0),
  min(est$year), max(est$year)
))

# ---------------------------------------------------------------- models ------

fit <- function(terms) {
  fixest::feols(
    nb$make_formula(paste(nb$base_climate, "+", terms)),
    data = est, panel.id = nb$pan_id, cluster = ~GID_1
  )
}

# fixest reads `x^k` in a formula as the power `I(x^k)` and names the coefficient
# accordingly. Writing `I(x^k)` instead would be wrapped a second time, so the
# formula uses the notebook's `^` syntax and the name is derived separately.
power_formula <- function(x, k) ifelse(k == 1L, x, sprintf("%s^%d", x, k))
power_coef <- function(x, k) ifelse(k == 1L, x, sprintf("I(%s^%d)", x, k))

degrees <- c(`1` = "Linear", `2` = "Quadratic", `3` = "Cubic")

# Every temperature term is a power of a non-negative hinge, so the implied
# response is symmetric in the anomaly and keeps one sign throughout. A power of
# the raw signed anomaly would not: odd powers such as zTM^3 are antisymmetric
# and would make cold years beneficial whenever hot years are harmful. The
# unhinged case is therefore |zTM|, which is the hinge at a threshold of zero.
threshold_label <- function(threshold) {
  ifelse(threshold == 0, "No hinge: |zTM|",
         sprintf("Hinge at %.1f SD", threshold))
}

specs <- tidyr::expand_grid(
  threshold = HINGES, degree = as.integer(names(degrees))
) %>%
  mutate(
    variable_name = hinge_col(threshold),
    id = paste0("hinge_", threshold, "_", degree),
    term_formula = power_formula(variable_name, degree),
    term_coef = power_coef(variable_name, degree),
    degree_label = factor(degrees[as.character(degree)], levels = unname(degrees)),
    panel_label = factor(threshold_label(threshold),
                         levels = threshold_label(HINGES))
  )

models <- setNames(
  lapply(paste(specs$term_formula, "+", PRECIP_TERM), fit),
  specs$id
)

m_bins <- fit("i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)")
m_base <- fit("zTM^2 + zRR^2")   # The notebook's quadratic benchmark.

# The notebook's current preferred specification, verbatim, to confirm that it
# coincides with the 1.5 SD hinge raised to the third power.
m_pref <- fit("zzTM^3 + abs_zRRp")
stopifnot(isTRUE(all.equal(unname(coef(m_pref)), unname(coef(models$hinge_1.5_3)))))

# ------------------------------------------------- implied response curves ----

# Temperature regressor implied by a specification, evaluated at anomaly z.
temp_basis <- function(z, threshold, degree) {
  matrix(pmax(0, abs(z) - threshold)^degree, ncol = 1)
}
rr_basis <- function(z) matrix(pmax(z, 0), ncol = 1)

# Delta-method curve for a single term, re-centred on its own average over the
# reference-bin observations so that it is comparable to the bin coefficients.
curve_from <- function(model, term, basis, ref_basis) {
  b <- coef(model)[term]
  v <- vcov(model)[term, term]
  stopifnot(!anyNA(b))
  a <- as.vector(basis) - mean(ref_basis, na.rm = TRUE)
  estimate <- a * b
  se <- abs(a) * sqrt(v)
  data.frame(estimate = estimate, se = se,
             ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)
}

grid <- seq(-GRID_LIMIT, GRID_LIMIT, length.out = 601)
tm_ref_z <- est$zTM[est$TM_bin_signed == 0]
rr_ref_z <- est$zRR[est$RR_bin_signed == 0]

curves <- purrr::pmap_dfr(
  list(specs$id, specs$threshold, specs$degree, specs$term_coef),
  function(id, threshold, degree, term_coef) {
    model <- models[[id]]
    bind_rows(
      cbind(
        variable = "Temperature", anomaly = grid,
        curve_from(
          model, term_coef,
          temp_basis(grid, threshold, degree),
          temp_basis(tm_ref_z, threshold, degree)
        )
      ),
      cbind(
        variable = "Precipitation", anomaly = grid,
        curve_from(model, PRECIP_TERM, rr_basis(grid), rr_basis(rr_ref_z))
      )
    ) %>%
      mutate(id = id, .before = 1)
  }
) %>%
  left_join(specs %>% select(id, degree_label, panel_label), by = "id")

# ------------------------------------------------------- bin coefficients -----

bin_points <- function(model) {
  ct <- as.data.frame(coeftable(model))
  ci <- confint(model)
  rn <- rownames(ct)
  sel <- grepl("^(TM|RR)_bin_signed::", rn)
  bind_rows(
    data.frame(
      variable = ifelse(startsWith(rn[sel], "TM"), "Temperature", "Precipitation"),
      anomaly = as.numeric(sub("^.*::", "", rn[sel])),
      estimate = ct[sel, 1], se = ct[sel, 2], p = ct[sel, 4],
      ci_low = ci[sel, 1], ci_high = ci[sel, 2]
    ),
    # The omitted bin is the -0.5 to 0.5 standard-deviation interval.
    data.frame(
      variable = c("Temperature", "Precipitation"), anomaly = 0,
      estimate = 0, se = NA_real_, p = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_
    )
  )
}

bins <- bin_points(m_bins) %>%
  # The outermost codes stand for open-ended tails beyond +/-2.5 SD.
  mutate(open_ended = abs(anomaly) > 2.5)

# ------------------------------------------------------------------ tests -----

# The parametric forms are not nested in the bins model, but the bins are nested
# in an encompassing model, so a joint Wald test on the bin indicators asks
# whether the nonparametric shape adds anything to the parametric response.
encompass <- function(label, terms, bin_terms, keep) {
  w <- wald(fit(paste(terms, "+", bin_terms)), keep = keep, print = FALSE)
  data.frame(specification = label, bins_tested = keep,
             F_stat = w$stat, df1 = w$df1, df2 = w$df2, p = w$p)
}

encompass_both <- function(label, terms) bind_rows(
  encompass(label, terms, "i(TM_bin_signed, ref = 0)", "TM_bin_signed"),
  encompass(label, terms, "i(RR_bin_signed, ref = 0)", "RR_bin_signed")
)

spec_label <- function(panel_label, degree_label) {
  paste0(panel_label, ", ", tolower(degree_label))
}

encompassing_tests <- bind_rows(
  purrr::pmap_dfr(
    list(specs$panel_label, specs$degree_label, specs$term_formula),
    function(panel_label, degree_label, term_formula) encompass_both(
      spec_label(panel_label, degree_label),
      paste(term_formula, "+", PRECIP_TERM)
    )
  ),
  purrr::map_dfr(HINGES, function(threshold) {
    stem <- hinge_col(threshold)
    encompass_both(
      paste(threshold_label(threshold), "cold and hot arms free", sep = ", "),
      sprintf("%s_cold + %s_hot + %s", stem, stem, PRECIP_TERM)
    )
  })
)

# Each hinge is symmetric in |zTM| by construction, so it cannot reproduce a
# response that differs between cold and hot tails. Splitting the hinge by sign
# tests that restriction without extra dependencies.
asym_models <- setNames(lapply(HINGES, function(threshold) {
  stem <- hinge_col(threshold)
  fit(sprintf("%s_cold + %s_hot + %s", stem, stem, PRECIP_TERM))
}), sprintf("%.1f", HINGES))

symmetry_tests <- purrr::map_dfr(HINGES, function(threshold) {
  model <- asym_models[[sprintf("%.1f", threshold)]]
  terms <- paste0(hinge_col(threshold), c("_cold", "_hot"))
  b <- coef(model)[terms]
  V <- vcov(model)[terms, terms]
  diff <- b[[1]] - b[[2]]
  var_diff <- V[1, 1] + V[2, 2] - 2 * V[1, 2]
  df2 <- fixest::degrees_freedom(model, "t")
  data.frame(
    threshold = threshold, cold_slope = b[[1]], hot_slope = b[[2]],
    difference = diff, se = sqrt(var_diff),
    F_stat = diff^2 / var_diff, df1 = 1, df2 = df2,
    p = pf(diff^2 / var_diff, 1, df2, lower.tail = FALSE)
  )
})

fit_row <- function(m, id) data.frame(
  model = id, n = nobs(m), k = length(coef(m)),
  within_r2 = as.numeric(fitstat(m, "wr2")[[1]]),
  r2 = as.numeric(fitstat(m, "r2")[[1]]),
  aic = AIC(m), bic = BIC(m)
)

model_fit <- purrr::imap_dfr(
  c(list(signed_bins = m_bins, notebook_quadratic_benchmark = m_base),
    setNames(asym_models, paste0("asym_hinge_", names(asym_models))),
    models),
  fit_row
) %>%
  left_join(
    specs %>% transmute(model = id, specification = spec_label(panel_label, degree_label)),
    by = "model"
  ) %>%
  arrange(desc(within_r2))

# ------------------------------------------------------------------- plots ----

variable_levels <- c("Temperature", "Precipitation")
degree_colors <- c(Linear = "#0072B2", Quadratic = "#D55E00", Cubic = "#009E73")
panel_levels <- levels(specs$panel_label)

# Restore the facet keys as factors so that every layer orders its panels alike.
prep <- function(x) {
  if ("variable" %in% names(x)) {
    x$variable <- factor(x$variable, levels = variable_levels)
  }
  if ("panel_label" %in% names(x)) {
    x$panel_label <- factor(x$panel_label, levels = panel_levels)
  }
  x
}

# A shared y-range keeps the two variables comparable; the raw cubic leaves the
# frame in the tails, which is itself part of what the figure shows.
y_limits <- range(c(bins$ci_low, bins$ci_high), na.rm = TRUE) * PP * 1.05

# `tidyr::crossing()` expands a factor over all of its levels rather than its
# observed values, so faceting keys are carried as character and `prep()`
# converts them back.
# A zero threshold has no kink to mark.
hinge_marks <- specs %>%
  filter(threshold > 0) %>%
  distinct(panel_label, threshold) %>%
  mutate(panel_label = as.character(panel_label)) %>%
  tidyr::crossing(variable = "Temperature", sign = c(-1, 1)) %>%
  mutate(xintercept = sign * threshold)

bin_layers <- function(bin_data) list(
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35),
  geom_errorbar(
    data = bin_data, aes(x = anomaly, ymin = ci_low * PP, ymax = ci_high * PP),
    width = 0.08, linewidth = 0.4, color = "grey25", inherit.aes = FALSE
  ),
  geom_point(
    data = bin_data,
    aes(x = anomaly, y = estimate * PP, shape = open_ended),
    color = "grey15", size = 1.9, fill = "white", inherit.aes = FALSE
  ),
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 21),
    labels = c(`FALSE` = "Interior bin", `TRUE` = "Open-ended tail bin"),
    name = "Nonparametric bins"
  ),
  scale_color_manual(values = degree_colors, name = "Power of the anomaly"),
  scale_fill_manual(values = degree_colors, name = "Power of the anomaly"),
  coord_cartesian(ylim = y_limits),
  labs(
    x = "Standardized anomaly, within-region standard deviations",
    y = "Log GDP-per-capita growth (percentage points)"
  ),
  theme_classic(base_size = 11),
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.9, "lines")
  )
)

# Figure 1: every hinge threshold and the unhinged benchmark, against the bins.
p_main <- ggplot() +
  bin_layers(prep(tidyr::crossing(bins, panel_label = panel_levels))) +
  geom_vline(
    data = prep(hinge_marks), aes(xintercept = xintercept),
    linetype = "dotted", color = "grey45", linewidth = 0.35
  ) +
  geom_line(
    data = prep(curves), aes(x = anomaly, y = estimate * PP, color = degree_label),
    linewidth = 0.8
  ) +
  facet_grid(panel_label ~ variable) +
  labs(
    title = "Preferred specification against nonparametric signed-anomaly bins",
    subtitle = paste(
      "Same estimation sample and climate controls throughout; bins measured",
      "against the -0.5 to 0.5 SD reference.\nDotted lines mark each hinge,",
      "below which the temperature term is zero by construction.",
      "\nOnly the temperature term varies, so the precipitation curves coincide."
    )
  )

ggsave(file.path(out, "signed_bins_vs_functional_form.png"),
       p_main, width = 11, height = 12, dpi = 180)

# Figure 2: temperature only, one panel per threshold and power, with bands.
p_bands <- ggplot() +
  bin_layers(prep(tidyr::crossing(
    bins %>% filter(variable == "Temperature") %>% select(-variable),
    specs %>% distinct(panel_label, degree_label) %>%
      mutate(panel_label = as.character(panel_label))
  ))) +
  geom_vline(
    data = prep(hinge_marks %>% select(-variable)),
    aes(xintercept = xintercept),
    linetype = "dotted", color = "grey45", linewidth = 0.35
  ) +
  geom_ribbon(
    data = prep(curves) %>% filter(variable == "Temperature"),
    aes(x = anomaly, ymin = ci_low * PP, ymax = ci_high * PP, fill = degree_label),
    alpha = 0.18
  ) +
  geom_line(
    data = prep(curves) %>% filter(variable == "Temperature"),
    aes(x = anomaly, y = estimate * PP, color = degree_label), linewidth = 0.8
  ) +
  facet_grid(panel_label ~ degree_label) +
  labs(
    title = "Temperature response with 95% confidence bands",
    subtitle = paste(
      "Delta-method intervals on the temperature term, clustered by GID_1.",
      "\nEvery term is a power of a non-negative hinge, so each response is",
      "symmetric in the anomaly and keeps one sign."
    ),
    x = "Standardized temperature anomaly, within-region standard deviations"
  )

ggsave(file.path(out, "temperature_bands_by_hinge.png"),
       p_bands, width = 11, height = 11, dpi = 180)

# Figure 3: diagnostic. Every form above is symmetric in the anomaly, so this
# frees the cold and hot arms of each hinge.
asym_curves <- purrr::map_dfr(HINGES, function(threshold) {
  model <- asym_models[[sprintf("%.1f", threshold)]]
  stem <- hinge_col(threshold)
  arm <- function(side, z) if (side == "cold") pmax(0, -z - threshold) else pmax(0, z - threshold)
  terms <- paste0(stem, c("_cold", "_hot"))
  basis <- cbind(arm("cold", grid), arm("hot", grid))
  ref <- cbind(arm("cold", tm_ref_z), arm("hot", tm_ref_z))
  b <- coef(model)[terms]
  V <- vcov(model)[terms, terms]
  A <- sweep(basis, 2, colMeans(ref, na.rm = TRUE), "-")
  estimate <- as.vector(A %*% b)
  se <- sqrt(pmax(0, rowSums((A %*% V) * A)))
  data.frame(
    panel_label = factor(threshold_label(threshold), levels = panel_levels),
    anomaly = grid, estimate = estimate,
    ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se
  )
})

symmetric_linear <- curves %>%
  filter(variable == "Temperature", degree_label == "Linear") %>%
  select(panel_label, anomaly, estimate)

p_asym <- ggplot() +
  bin_layers(prep(tidyr::crossing(
    bins %>% filter(variable == "Temperature") %>% select(-variable),
    panel_label = panel_levels
  ))) +
  geom_vline(
    data = prep(hinge_marks %>% select(-variable)),
    aes(xintercept = xintercept),
    linetype = "dotted", color = "grey45", linewidth = 0.35
  ) +
  geom_ribbon(
    data = asym_curves,
    aes(x = anomaly, ymin = ci_low * PP, ymax = ci_high * PP),
    fill = "#CC79A7", alpha = 0.18
  ) +
  geom_line(
    data = symmetric_linear,
    aes(x = anomaly, y = estimate * PP, color = "Symmetric in |anomaly|"),
    linewidth = 0.8
  ) +
  geom_line(
    data = asym_curves,
    aes(x = anomaly, y = estimate * PP, color = "Cold and hot arms free"),
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(`Symmetric in |anomaly|` = "#0072B2",
               `Cold and hot arms free` = "#CC79A7"),
    name = "Linear hinge"
  ) +
  facet_wrap(~panel_label, nrow = 1) +
  labs(
    title = "Diagnostic: is symmetry, not the power, what the bins reject?",
    subtitle = paste(
      "Every specification above is symmetric in the anomaly, so none can trace",
      "the deeper cold tail.",
      sprintf("\nCold vs hot slope equality: p = %s.",
              paste(sprintf("%.2g at %.1f SD", symmetry_tests$p,
                            symmetry_tests$threshold), collapse = ", "))
    ),
    x = "Standardized temperature anomaly, within-region standard deviations"
  )

ggsave(file.path(out, "asymmetry_diagnostic.png"), p_asym,
       width = 15, height = 5.5, dpi = 180)

# ------------------------------------------------------------------ output ----

write_out(bins, "signed_bin_coefficients")
write_out(curves %>% select(-degree_label, -panel_label), "implied_response_curves")
write_out(model_fit, "model_fit_comparison")
write_out(encompassing_tests, "encompassing_tests")
write_out(symmetry_tests, "hinge_symmetry_tests")

all_models <- c(
  models,
  setNames(asym_models, paste0("asym_hinge_", names(asym_models))),
  list(signed_bins = m_bins, notebook_quadratic_benchmark = m_base)
)

coef_table <- purrr::imap_dfr(all_models, function(m, id) {
  ct <- as.data.frame(coeftable(m))
  ci <- confint(m)
  data.frame(model = id, term = rownames(ct), estimate = ct[, 1],
             se = ct[, 2], p = ct[, 4], ci_low = ci[, 1], ci_high = ci[, 2])
})
write_out(coef_table, "coefficients")

# ---------------------------------------------- base-coefficient stability ----

# The BHM climate controls are common to every specification, so the spread of
# their estimates across weather forms shows whether the temperature term is
# taking variation those controls would otherwise claim. The spread is reported
# in units of the clustered standard error, which is the scale that matters for
# inference: a range well below one standard error is sampling noise.
BASE_TERMS <- c("TM", "TM:mean_TM_all", "RR", "RR:mean_RR_all")

base_coefs <- coef_table %>%
  filter(term %in% BASE_TERMS) %>%
  mutate(
    term = factor(term, levels = BASE_TERMS),
    # The bins model replaces abs_zRRp with precipitation bins, so it is not a
    # like-for-like comparison for the two RR controls.
    comparable = model != "signed_bins"
  )

summarise_base <- function(x, scope) {
  x %>%
    group_by(term) %>%
    summarise(
      scope = scope, n_models = n(),
      min = min(estimate), max = max(estimate), mean = mean(estimate),
      sd = sd(estimate), median_se = median(se),
      spread = max(estimate) - min(estimate),
      spread_in_se = (max(estimate) - min(estimate)) / median(se),
      sd_in_se = sd(estimate) / median(se),
      same_sign = n_distinct(sign(estimate)) == 1L,
      all_p05 = all(p < 0.05),
      .groups = "drop"
    ) %>%
    relocate(scope, .after = term)
}

base_stability <- bind_rows(
  summarise_base(base_coefs %>% filter(model %in% specs$id), "12 hinge powers"),
  summarise_base(base_coefs %>% filter(comparable), "all but signed bins"),
  summarise_base(base_coefs, "all models")
) %>%
  arrange(term, scope)

write_out(base_coefs %>% select(-comparable), "base_coefficients")
write_out(base_stability, "base_coefficient_stability")

# The temperature controls do not wander at random across the 12 hinge powers:
# they track how well the weather term fits. A deviation term that explains
# little leaves nonlinearity for TM and its long-run interaction to absorb, so a
# rank correlation near +1 or -1 means the base estimate is an artefact of the
# chosen functional form rather than an independent finding.
base_vs_fit <- base_coefs %>%
  filter(model %in% specs$id) %>%
  left_join(model_fit %>% select(model, within_r2), by = "model") %>%
  group_by(term) %>%
  summarise(
    n_models = n(),
    spearman_estimate_vs_within_r2 = cor(estimate, within_r2, method = "spearman"),
    estimate_at_best_fit = estimate[which.max(within_r2)],
    estimate_at_worst_fit = estimate[which.min(within_r2)],
    ratio_worst_to_best = estimate[which.min(within_r2)] /
      estimate[which.max(within_r2)],
    .groups = "drop"
  )
write_out(base_vs_fit, "base_coefficients_vs_fit")

p_base <- ggplot(
  base_coefs %>% mutate(model = factor(model, levels = rev(sort(unique(model))))),
  aes(x = estimate * PP, y = model, color = comparable)
) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low * PP, xmax = ci_high * PP),
                orientation = "y", width = 0.25, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"),
    labels = c(`TRUE` = "Same precipitation term",
               `FALSE` = "Precipitation bins instead"),
    name = NULL
  ) +
  facet_wrap(~term, scales = "free_x", nrow = 1) +
  labs(
    title = "Base climate controls across weather specifications",
    subtitle = paste(
      "The BHM controls are common to every model; only the temperature term",
      "and its hinge change.\n95% intervals clustered by GID_1."
    ),
    x = "Coefficient in log-growth percentage points", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out, "base_coefficient_stability.png"), p_base,
       width = 13, height = 6, dpi = 180)

etable_text <- capture.output(etable(
  c(models, list(bins = m_bins)),
  fitstat = ~n + r2 + wr2 + aic + bic
))
writeLines(etable_text, file.path(out, "model_comparison.txt"))

cat("\n=== Fit comparison, best within-R2 first ===\n")
print(model_fit, row.names = FALSE)
cat("\n=== Do the bins add to the parametric form? ===\n")
print(encompassing_tests, row.names = FALSE)
cat("\n=== Is each hinge symmetric between cold and hot tails? ===\n")
print(symmetry_tests, row.names = FALSE)
cat("\n=== Base climate controls: spread across specifications ===\n")
print(
  base_stability %>%
    mutate(across(c(min, max, mean, sd, median_se, spread), ~ signif(.x, 3)),
           across(c(spread_in_se, sd_in_se), ~ round(.x, 3))),
  row.names = FALSE
)
cat("\n=== Base controls vs how well the weather term fits (12 hinge powers) ===\n")
print(as.data.frame(base_vs_fit), row.names = FALSE, digits = 3)
cat("\nWrote figures and tables to ", out, "\n", sep = "")
