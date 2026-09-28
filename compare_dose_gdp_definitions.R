# Signed-bin climate response across DOSE GDP-per-capita definitions and
# climate datasets.
#
# Specification and defaults as in compare_signed_bins_across_datasets.R (and
# `model_bins_signed` in test_functions.Rmd): BHM base terms + signed half-SD
# anomaly bins, fixed effects year + GID_1[year], errors clustered by GID_1,
# econ_year_min = 1950.
#
# Grid: the five DOSE definitions (DOSE_GDP_VARIABLES in functions/econ_panel.R)
# x climate source ERA5, CRU TS, UDelaware, all with concurrent-population
# weights = 15 panels.
#   lcu2015_usd  constant 2015 local prices at the 2015 exchange rate (default)
#   lcu_2015     constant 2015 local currency; identical growth to lcu2015_usd,
#                so its models must reproduce lcu2015_usd exactly
#   usd_2015     current US$ divided by the US GDP deflator
#   lcu          current local currency: nominal, includes inflation
#   usd          current US$: nominal, includes US inflation and exchange rates
#
# Samples:
#   max     each panel's own estimation sample;
#   common  region-years estimable in all 15 panels.
#
# Output: results/signed_bins_dose_gdp_definitions/. Part 1 estimates and
# writes the analysis tables; part 2 plots.

source("load_functions.R")

# Configuration ----------------------------------------------------------------

DEFINITIONS <- c(
  lcu2015_usd = "Const. LCU, 2015 FX",
  lcu_2015 = "Const. LCU",
  usd_2015 = "US$ / US deflator",
  lcu = "Nominal LCU",
  usd = "Nominal US$"
)
CLIMATE_SOURCES <- c(ERA5 = "ERA5", CRU = "CRU TS", UDEL = "UDelaware")
CLIMATE_WEIGHT <- "concurrent population"

# Notebook defaults (test_functions.Rmd, shared-specification chunk).
OUTCOME <- "dlgrp_pc_usd"
FIXED_EFFECTS <- "year + GID_1[year]"
PANEL_ID <- c("GID_1", "year")
ECON_YEAR_MIN <- 1950L
CLIMATE_VELOCITY_YEARS <- 5L
MODEL_TERMS <- paste(
  "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
  "+ i(TM_bin_signed, ref = 0)",
  "+ i(RR_bin_signed, ref = 0)"
)
KEEP_VARS <- c("GID_0", "GID_1", "year", OUTCOME, "TM", "RR", "mean_TM_all",
               "mean_RR_all", "TM_bin_signed", "RR_bin_signed", "zTM", "zRR")

OUTPUT_DIR <- file.path("results", "signed_bins_dose_gdp_definitions")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

fit_signed_bins <- function(model_data) {
  fit_panel_model(
    MODEL_TERMS,
    model_data,
    outcome = OUTCOME,
    fixed_effects = FIXED_EFFECTS,
    panel_id = PANEL_ID,
    cluster = ~GID_1
  )
}

# Panels -------------------------------------------------------------------------

specs <- expand.grid(
  definition_code = names(DEFINITIONS),
  climate = names(CLIMATE_SOURCES),
  stringsAsFactors = FALSE
) %>%
  mutate(
    definition = unname(DEFINITIONS[definition_code]),
    id = paste(definition_code, climate, sep = "_")
  )

panels <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  message("Building ", s$id)
  cfg <- panel_config(
    econ_year_min = ECON_YEAR_MIN,
    climate_source = CLIMATE_SOURCES[[s$climate]],
    climate_weight = CLIMATE_WEIGHT,
    clim_history_years = max(6L, CLIMATE_VELOCITY_YEARS - 1L)
  )
  panels[[s$id]] <- build_dat(
    econ_data = "DOSE",
    config = cfg,
    econ_variable = s$definition_code
  ) %>%
    select(all_of(KEEP_VARS))
}

# Samples ------------------------------------------------------------------------

fits_max <- lapply(panels, fit_signed_bins)
used_keys <- lapply(names(fits_max), function(id) {
  panels[[id]][fixest::obs(fits_max[[id]]), KEY_COLS]
})
common_keys <- intersect_keys(used_keys)
common_data <- lapply(panels, function(d) semi_join(d, common_keys, by = KEY_COLS))
fits_common <- lapply(common_data, fit_signed_bins)

# Part 1: analysis tables ----------------------------------------------------------
#
# Table builders are in functions/specification_grid.R.

FACTORS <- c("definition", "climate")
all_fits <- list(max = fits_max, common = fits_common)
all_data <- list(max = panels, common = common_data)

sample_summary <- grid_sample_summary(all_fits, all_data, specs)
write_out(sample_summary, "sample_summary")

coefficients <- grid_bin_coefficients(all_fits, specs)
write_out(coefficients, "bin_coefficients")

base_terms <- grid_base_terms(all_fits, specs)
write_out(base_terms, "base_coefficients")

joint_tests <- grid_joint_tests(all_fits, specs)
write_out(joint_tests, "joint_tests")

variance_shares <- grid_dispersion(coefficients, FACTORS)
write_out(variance_shares, "cross_specification_dispersion")

factor_means <- grid_factor_means(coefficients, FACTORS)
write_out(factor_means, "mean_estimate_by_factor")

data_quality <- grid_data_quality(all_fits, all_data, specs)
write_out(data_quality, "data_quality")

# lcu_2015 and lcu2015_usd differ by one constant per country, which the
# region fixed effects absorb: their estimates must coincide.
identity_check <- coefficients %>%
  filter(definition_code %in% c("lcu2015_usd", "lcu_2015")) %>%
  select(sample, climate, bin_var, bin_label, definition_code, estimate) %>%
  tidyr::pivot_wider(names_from = definition_code, values_from = estimate) %>%
  summarise(max_abs_difference = max(abs(lcu2015_usd - lcu_2015)))
write_out(identity_check, "identity_check_lcu_definitions")

# How different the outcome itself is across definitions on the common sample:
# correlation and relative spread of year-demeaned regional growth (year fixed
# effects remove anything common to all regions in a year, such as US
# inflation). The outcome does not depend on the climate data, so one climate
# source is enough.
growth_by_definition <- bind_rows(lapply(names(DEFINITIONS), function(code) {
  common_data[[paste(code, "ERA5", sep = "_")]] %>%
    select(GID_1, year, growth = all_of(OUTCOME)) %>%
    mutate(definition_code = code)
})) %>%
  group_by(definition_code, year) %>%
  mutate(growth_demeaned = growth - mean(growth)) %>%
  ungroup()

growth_wide <- growth_by_definition %>%
  select(GID_1, year, definition_code, growth_demeaned) %>%
  tidyr::pivot_wider(names_from = definition_code, values_from = growth_demeaned)
definition_pairs <- utils::combn(names(DEFINITIONS), 2, simplify = FALSE)
growth_agreement <- bind_rows(lapply(definition_pairs, function(p) {
  tibble::tibble(
    variant_1 = p[1],
    variant_2 = p[2],
    observations = sum(stats::complete.cases(growth_wide[p])),
    correlation = stats::cor(growth_wide[[p[1]]], growth_wide[[p[2]]],
                             use = "complete.obs"),
    sd_ratio = stats::sd(growth_wide[[p[2]]], na.rm = TRUE) /
      stats::sd(growth_wide[[p[1]]], na.rm = TRUE)
  )
}))
write_out(growth_agreement, "growth_agreement")

growth_summary <- growth_by_definition %>%
  group_by(definition_code) %>%
  summarise(
    observations = n(),
    mean_growth = mean(growth),
    sd_growth = stats::sd(growth),
    sd_growth_year_demeaned = stats::sd(growth_demeaned),
    p01 = stats::quantile(growth, 0.01),
    p99 = stats::quantile(growth, 0.99),
    abs_log_growth_above_0_5 = sum(abs(growth) > 0.5),
    .groups = "drop"
  )
write_out(growth_summary, "growth_summary")

saveRDS(
  list(specs = specs, sample_summary = sample_summary, coefficients = coefficients,
       base_terms = base_terms, joint_tests = joint_tests,
       variance_shares = variance_shares, factor_means = factor_means,
       data_quality = data_quality, identity_check = identity_check,
       growth_agreement = growth_agreement, growth_summary = growth_summary),
  file.path(OUTPUT_DIR, "analysis.rds")
)

options(width = 200)
cat("\n== Samples ==\n")
print(sample_summary %>% select(sample, definition, climate, observations, units,
                                countries, first_year, last_year, within_r2),
      n = Inf)
cat("\n== Outcome by definition (common sample) ==\n")
print(growth_summary, n = Inf)
cat("\n== Agreement of year-demeaned growth across definitions ==\n")
print(growth_agreement, n = Inf)
cat("\n== lcu_2015 vs lcu2015_usd estimates ==\n")
print(identity_check)
cat("\n== Joint tests (p-values) ==\n")
print(joint_tests %>% select(sample, definition, climate, test, p_value) %>%
        tidyr::pivot_wider(names_from = test, values_from = p_value), n = Inf)
cat("\n== Cross-specification dispersion per bin ==\n")
print(variance_shares, n = Inf)
cat("\n== Mean estimate by definition ==\n")
print(factor_means %>% filter(factor == "definition") %>%
        mutate(mean_estimate = round(100 * mean_estimate, 2)) %>%
        tidyr::pivot_wider(names_from = factor_level, values_from = mean_estimate),
      n = Inf)
cat("\n== Data quality of estimation rows ==\n")
print(data_quality %>% select(sample, definition, climate, rows, zero_precipitation,
                              share_extreme_TM_bins, share_extreme_RR_bins),
      n = Inf)

# Part 2: plots ----------------------------------------------------------------------
#
# Plot builders are in functions/specification_grid.R: one column per
# definition, one colour and marker per climate source.

SAMPLE_LABELS <- c(
  max = "Maximum sample",
  common = "Common sample (1970-2017, 75 countries)"
)

plot_data <- grid_plot_data(
  coefficients,
  list(definition = unname(DEFINITIONS), climate = names(CLIMATE_SOURCES))
)

# Figure 1: every specification, one figure per sample.
for (s in names(SAMPLE_LABELS)) {
  p <- plot_bin_grid(
    plot_data %>% filter(sample == s),
    facet = "definition",
    colour = "climate",
    colour_name = "Climate data (concurrent-population weights)",
    title = paste("Signed-bin response by DOSE GDP definition and climate data:",
                  tolower(SAMPLE_LABELS[[s]])),
    subtitle = paste("BHM base terms + signed half-SD anomaly bins;",
                     "reference bin -0.5 to 0.5 SD; 95% CI clustered by region"),
    caption = paste(
      "Const. LCU, 2015 FX and Const. LCU have identical growth, as do",
      "US$ / US deflator and Nominal US$ once year fixed effects are included."
    )
  )
  ggsave(file.path(OUTPUT_DIR, paste0("signed_bins_", s, "_sample.png")), p,
         width = 13, height = 8, dpi = 300)
}

# Figure 2: the range across all 15 specifications for each bin, by sample.
ggsave(
  file.path(OUTPUT_DIR, "signed_bins_cross_specification_range.png"),
  plot_bin_range(
    plot_data,
    SAMPLE_LABELS,
    title = "How much the signed-bin response depends on the DOSE GDP definition",
    subtitle = paste("Point: mean over the 15 definition x climate specifications;",
                     "thick bar: interquartile range; thin line: min-max")
  ),
  width = 11, height = 5.5, dpi = 300
)

# Figure 3: agreement of the outcome across definitions (common sample).
agreement_data <- growth_agreement %>%
  mutate(
    variant_1 = factor(DEFINITIONS[variant_1], levels = DEFINITIONS),
    variant_2 = factor(DEFINITIONS[variant_2], levels = DEFINITIONS),
    value = correlation
  )
ggsave(
  file.path(OUTPUT_DIR, "growth_agreement.png"),
  plot_agreement_heatmap(
    agreement_data,
    fill_name = "Correlation",
    title = "Five DOSE definitions, three distinct growth series",
    subtitle = paste("Correlation of year-demeaned regional growth,",
                     "common sample (1970-2017)"),
    limits = c(0, 1),
    label_format = function(x) sprintf("%.2f", x),
    fill_labels = waiver()
  ),
  width = 9, height = 6.5, dpi = 300
)

message("Figures written to ", normalizePath(OUTPUT_DIR))
