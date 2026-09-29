# Compare DOSE climate-response estimates with and without currency movements,
# and with and without country-by-year fixed effects.
#
# Outcomes in the clean currency-only experiment:
#   no_currency_fluctuation    Growth of GDP per capita at constant 2015 local
#                              prices.
#   with_currency_fluctuation  The same real-output growth plus the change in
#                              DOSE's country-year exchange rate.
#
# All models use the same cleaned observations, population-weighted CRU climate,
# year effects, region-specific linear trends and region-clustered standard
# errors. The second FE design replaces the year effects with country-by-year
# effects. If the currency component is common to all regions in a country-year,
# those effects must absorb it and make the slope estimates and residuals equal.

source("load_functions.R")

OUTPUT_DIR <- file.path("results", "dose_currency_country_year_fe")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTCOMES <- c(
  no_currency_fluctuation = "no_currency_fluctuation",
  with_currency_fluctuation = "with_currency_fluctuation"
)

FE_DESIGNS <- c(
  without_country_year_fe = "year + GID_1[year]",
  with_country_year_fe = "GID_0^year + GID_1[year]"
)

BASE_TERMS <- "TM + TM:mean_TM_all + RR + RR:mean_RR_all"
MODEL_TERMS <- c(
  base = BASE_TERMS,
  quadratic = paste(BASE_TERMS, "+ I(zTM^2) + I(zRR^2)"),
  symmetric_hinge = paste(BASE_TERMS, "+ I(hTM_15^2) + abs_zRRp"),
  asymmetric_hinge = paste(
    BASE_TERMS,
    "+ hTM_15_cold + hTM_15_hot + abs_zRRp"
  ),
  signed_bins = paste(
    BASE_TERMS,
    "+ i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)"
  )
)

KEYS <- c("GID_0", "GID_1", "year")
MODEL_VARS <- c(
  KEYS, "dlgrp_pc_usd", "TM", "RR", "mean_TM_all", "mean_RR_all",
  "zTM", "zRR", "abs_zRRp", "TM_bin_signed", "RR_bin_signed"
)

cfg <- panel_config(
  econ_year_min = 1950L,
  climate_source = "CRU TS",
  climate_weight = "concurrent population",
  clim_history_years = 6L
)

# Build the real-LCU outcome once, then add only the observed exchange-rate
# change. This holds the underlying volume series fixed and makes the treatment
# of currency the sole difference between the two outcomes. DOSE's `fx` is in
# local-currency units per US dollar, so USD growth subtracts dlog(fx).
fx_raw <- read_dose(data_file(ECON_RAW_FILES[["DOSE_V2_14"]])) %>%
  transmute(GID_0, year, fx) %>%
  filter(!is.na(fx), fx > 0)
fx_by_country_year <- fx_raw %>%
  group_by(GID_0, year) %>%
  summarise(
    regional_fx_values = n_distinct(fx),
    regional_fx_range = max(fx) - min(fx),
    fx = stats::median(fx),
    .groups = "drop"
  )
fx_source_identity <- fx_by_country_year %>%
  summarise(
    country_years = n(),
    country_years_with_multiple_fx_values = sum(regional_fx_values > 1),
    max_regional_fx_range = max(regional_fx_range)
  )

base_panel <- build_dat(
  econ_data = "DOSE",
  config = cfg,
  econ_variable = "lcu_2015",
  drop_break_years = TRUE
) %>%
  select(all_of(MODEL_VARS)) %>%
  left_join(fx_by_country_year %>% select(GID_0, year, fx),
            by = c("GID_0", "year")) %>%
  group_by(GID_1) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    log_fx = if_else(!is.na(fx) & fx > 0, log(fx), NA_real_),
    currency_component = -(log_fx - lag_by_year(log_fx, year)),
    growth_no_currency_fluctuation = dlgrp_pc_usd,
    growth_with_currency_fluctuation = dlgrp_pc_usd + currency_component
  ) %>%
  ungroup() %>%
  filter(
    !is.na(growth_no_currency_fluctuation),
    !is.na(growth_with_currency_fluctuation)
  ) %>%
  mutate(
    hTM_15 = pmax(1.5, abs(zTM)) - 1.5,
    hTM_15_cold = pmax(-zTM - 1.5, 0),
    hTM_15_hot = pmax(zTM - 1.5, 0)
  ) %>%
  arrange(across(all_of(KEYS)))

panels <- list(
  no_currency_fluctuation = base_panel %>%
    mutate(dlgrp_pc_usd = growth_no_currency_fluctuation),
  with_currency_fluctuation = base_panel %>%
    mutate(dlgrp_pc_usd = growth_with_currency_fluctuation)
)

common_keys <- Reduce(
  function(x, y) inner_join(x, y, by = KEYS),
  lapply(panels, function(x) distinct(x, across(all_of(KEYS))))
)
panels <- lapply(panels, function(x) semi_join(x, common_keys, by = KEYS))

# Climate regressors must be numerically identical before comparing outcomes.
regressor_cols <- setdiff(names(panels[[1]]), c(KEYS, "dlgrp_pc_usd"))
regressor_check <- inner_join(
  panels[[1]] %>% select(all_of(KEYS), all_of(regressor_cols)),
  panels[[2]] %>% select(all_of(KEYS), all_of(regressor_cols)),
  by = KEYS,
  suffix = c("_lcu", "_usd")
)
numeric_regressors <- regressor_cols[vapply(
  panels[[1]][regressor_cols], is.numeric, logical(1)
)]
factor_regressors <- setdiff(regressor_cols, numeric_regressors)

regressor_identity <- bind_rows(
  lapply(numeric_regressors, function(v) {
    difference <- regressor_check[[paste0(v, "_lcu")]] -
      regressor_check[[paste0(v, "_usd")]]
    tibble::tibble(
      variable = v,
      type = "numeric",
      mismatches = sum(abs(difference) > 1e-12, na.rm = TRUE),
      max_abs_difference = if (all(is.na(difference))) NA_real_ else
        max(abs(difference), na.rm = TRUE)
    )
  }),
  lapply(factor_regressors, function(v) {
    left <- as.character(regressor_check[[paste0(v, "_lcu")]])
    right <- as.character(regressor_check[[paste0(v, "_usd")]])
    different <- xor(is.na(left), is.na(right)) |
      (!is.na(left) & !is.na(right) & left != right)
    tibble::tibble(
      variable = v,
      type = "categorical",
      mismatches = sum(different),
      max_abs_difference = NA_real_
    )
  })
)

if (any(regressor_identity$mismatches > 0)) {
  stop("The two common-sample panels do not have identical regressors.")
}

# Quantify whether the outcome difference is exactly common within each
# country-year, which is the algebraic condition behind the FE identity.
outcomes_wide <- inner_join(
  panels$no_currency_fluctuation %>%
    select(all_of(KEYS), growth_lcu = dlgrp_pc_usd),
  panels$with_currency_fluctuation %>%
    select(all_of(KEYS), growth_usd = dlgrp_pc_usd),
  by = KEYS
) %>%
  mutate(currency_component = growth_usd - growth_lcu) %>%
  group_by(GID_0, year) %>%
  mutate(
    currency_component_country_year_mean = mean(currency_component),
    currency_component_within_country_year =
      currency_component - currency_component_country_year_mean
  ) %>%
  ungroup()

outcome_identity <- outcomes_wide %>%
  summarise(
    observations = n(),
    countries = n_distinct(GID_0),
    regions = n_distinct(GID_1),
    first_year = min(year),
    last_year = max(year),
    max_abs_within_country_year_currency_component =
      max(abs(currency_component_within_country_year)),
    rms_within_country_year_currency_component =
      sqrt(mean(currency_component_within_country_year^2)),
    correlation_raw_growth = cor(growth_lcu, growth_usd),
    sd_ratio_raw_growth = sd(growth_usd) / sd(growth_lcu)
  )

# Published `usd_2015` is retained as a diagnostic, but it is not the clean
# currency-only outcome: it can differ within country-year from `lcu_2015`.
published_usd <- build_dat(
  econ_data = "DOSE",
  config = cfg,
  econ_variable = "usd_2015",
  drop_break_years = TRUE
) %>%
  filter(!is.na(dlgrp_pc_usd)) %>%
  select(all_of(KEYS), growth_published_usd = dlgrp_pc_usd)

published_outcome_identity <- panels$no_currency_fluctuation %>%
  select(all_of(KEYS), growth_lcu = dlgrp_pc_usd) %>%
  inner_join(published_usd, by = KEYS) %>%
  mutate(published_difference = growth_published_usd - growth_lcu) %>%
  group_by(GID_0, year) %>%
  mutate(
    published_difference_within_country_year = published_difference -
      mean(published_difference)
  ) %>%
  ungroup() %>%
  summarise(
    observations = n(),
    max_abs_within_country_year_published_difference =
      max(abs(published_difference_within_country_year)),
    rms_within_country_year_published_difference =
      sqrt(mean(published_difference_within_country_year^2)),
    correlation_raw_growth = cor(growth_lcu, growth_published_usd),
    sd_ratio_raw_growth = sd(growth_published_usd) / sd(growth_lcu)
  )

# Estimate the complete 2 x 2 design for every existing functional form.
fits <- list()
for (form in names(MODEL_TERMS)) {
  for (outcome_name in names(OUTCOMES)) {
    for (fe_name in names(FE_DESIGNS)) {
      id <- paste(form, outcome_name, fe_name, sep = "__")
      message("Estimating ", id)
      fits[[id]] <- fit_panel_model(
        MODEL_TERMS[[form]],
        panels[[outcome_name]],
        fixed_effects = FE_DESIGNS[[fe_name]],
        panel_id = c("GID_1", "year"),
        cluster = ~GID_1,
        fixef.tol = 1e-10,
        fixef.iter = 100000
      )
    }
  }
}

fit_index <- expand.grid(
  form = names(MODEL_TERMS),
  outcome = names(OUTCOMES),
  fe_design = names(FE_DESIGNS),
  stringsAsFactors = FALSE
) %>%
  mutate(id = paste(form, outcome, fe_design, sep = "__"))

model_summary <- bind_rows(lapply(seq_len(nrow(fit_index)), function(i) {
  row <- fit_index[i, ]
  fit <- fits[[row$id]]
  tibble::tibble(
    form = row$form,
    outcome = row$outcome,
    currency_fluctuation = row$outcome == "with_currency_fluctuation",
    fe_design = row$fe_design,
    country_year_fe = row$fe_design == "with_country_year_fe",
    observations = stats::nobs(fit),
    coefficients = length(stats::coef(fit)),
    within_r2 = unname(fixest::fitstat(fit, "wr2")[[1]]),
    rmse = sqrt(mean(stats::residuals(fit)^2)),
    aic = stats::AIC(fit),
    bic = stats::BIC(fit)
  )
}))

coefficients <- bind_rows(lapply(seq_len(nrow(fit_index)), function(i) {
  row <- fit_index[i, ]
  fit <- fits[[row$id]]
  ci <- stats::confint(fit)
  tibble::tibble(
    form = row$form,
    outcome = row$outcome,
    currency_fluctuation = row$outcome == "with_currency_fluctuation",
    fe_design = row$fe_design,
    country_year_fe = row$fe_design == "with_country_year_fe",
    term = names(stats::coef(fit)),
    estimate = unname(stats::coef(fit)),
    std_error = unname(fixest::se(fit)),
    p_value = unname(fixest::pvalue(fit)),
    conf_low = ci[, 1],
    conf_high = ci[, 2]
  )
}))

# Compare the two outcome definitions within each FE design. With country-year
# effects, every reported difference should be floating-point noise only.
coefficient_comparison <- coefficients %>%
  select(form, fe_design, term, outcome, estimate, std_error, p_value) %>%
  tidyr::pivot_wider(
    names_from = outcome,
    values_from = c(estimate, std_error, p_value)
  ) %>%
  mutate(
    estimate_difference = estimate_with_currency_fluctuation -
      estimate_no_currency_fluctuation,
    std_error_difference = std_error_with_currency_fluctuation -
      std_error_no_currency_fluctuation,
    p_value_difference = p_value_with_currency_fluctuation -
      p_value_no_currency_fluctuation
  )

equality_checks <- bind_rows(lapply(names(MODEL_TERMS), function(form) {
  bind_rows(lapply(names(FE_DESIGNS), function(fe_name) {
    fit_lcu <- fits[[paste(form, "no_currency_fluctuation", fe_name, sep = "__")]]
    fit_usd <- fits[[paste(form, "with_currency_fluctuation", fe_name, sep = "__")]]
    coef_lcu <- stats::coef(fit_lcu)
    coef_usd <- stats::coef(fit_usd)
    common_terms <- intersect(names(coef_lcu), names(coef_usd))
    residual_lcu <- stats::residuals(fit_lcu)
    residual_usd <- stats::residuals(fit_usd)
    tibble::tibble(
      form = form,
      fe_design = fe_name,
      country_year_fe = fe_name == "with_country_year_fe",
      observations_lcu = stats::nobs(fit_lcu),
      observations_usd = stats::nobs(fit_usd),
      terms_lcu = length(coef_lcu),
      terms_usd = length(coef_usd),
      max_abs_coefficient_difference =
        max(abs(coef_usd[common_terms] - coef_lcu[common_terms])),
      max_abs_standard_error_difference = max(abs(
        fixest::se(fit_usd)[common_terms] - fixest::se(fit_lcu)[common_terms]
      )),
      max_abs_residual_difference = max(abs(residual_usd - residual_lcu)),
      rmse_difference =
        sqrt(mean(residual_usd^2)) - sqrt(mean(residual_lcu^2)),
      within_r2_difference =
        unname(fixest::fitstat(fit_usd, "wr2")[[1]]) -
        unname(fixest::fitstat(fit_lcu, "wr2")[[1]])
    )
  }))
})) %>%
  mutate(
    equal_at_1e_9 = country_year_fe &
      observations_lcu == observations_usd &
      terms_lcu == terms_usd &
      max_abs_coefficient_difference < 1e-9 &
      max_abs_standard_error_difference < 1e-9 &
      max_abs_residual_difference < 1e-9
  )

utils::write.csv(
  model_summary,
  file.path(OUTPUT_DIR, "model_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  coefficients,
  file.path(OUTPUT_DIR, "coefficients.csv"),
  row.names = FALSE
)
utils::write.csv(
  coefficient_comparison,
  file.path(OUTPUT_DIR, "coefficient_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  equality_checks,
  file.path(OUTPUT_DIR, "equality_checks.csv"),
  row.names = FALSE
)
utils::write.csv(
  outcome_identity,
  file.path(OUTPUT_DIR, "outcome_country_year_identity.csv"),
  row.names = FALSE
)
utils::write.csv(
  published_outcome_identity,
  file.path(OUTPUT_DIR, "published_usd_country_year_diagnostic.csv"),
  row.names = FALSE
)
utils::write.csv(
  regressor_identity,
  file.path(OUTPUT_DIR, "regressor_identity.csv"),
  row.names = FALSE
)
utils::write.csv(
  fx_source_identity,
  file.path(OUTPUT_DIR, "fx_source_identity.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    model_summary = model_summary,
    coefficients = coefficients,
    coefficient_comparison = coefficient_comparison,
    equality_checks = equality_checks,
    outcome_identity = outcome_identity,
    published_outcome_identity = published_outcome_identity,
    regressor_identity = regressor_identity,
    fx_source_identity = fx_source_identity
  ),
  file.path(OUTPUT_DIR, "analysis.rds")
)

options(width = 180)
cat("\n== Outcome difference within country-year ==\n")
print(outcome_identity)
cat("\n== Published usd_2015 diagnostic (not a currency-only contrast) ==\n")
print(published_outcome_identity)
cat("\n== Raw DOSE regional FX duplication diagnostic ==\n")
print(fx_source_identity)
cat("\n== Equality checks ==\n")
print(equality_checks, n = Inf)
cat("\n== Asymmetric-hinge coefficient comparison ==\n")
print(
  coefficient_comparison %>%
    filter(form == "asymmetric_hinge") %>%
    select(
      fe_design, term,
      estimate_no_currency_fluctuation,
      estimate_with_currency_fluctuation,
      estimate_difference
    ),
  n = Inf
)

if (!all(equality_checks$equal_at_1e_9[equality_checks$country_year_fe])) {
  stop("At least one country-year-FE specification failed the 1e-9 identity check.")
}

message("Results written to ", normalizePath(OUTPUT_DIR))
