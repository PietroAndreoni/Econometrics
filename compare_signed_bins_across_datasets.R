# Signed-bin climate response across economic and climate datasets.
#
# Specification: the BHM base terms plus signed half-SD temperature and
# precipitation anomaly bins (`model_bins_signed` in test_functions.Rmd), with
# the notebook's defaults: outcome dlgrp_pc_usd, fixed effects year +
# GID_1[year], errors clustered by GID_1, panel_config() defaults with
# econ_year_min = 1950 and six years of climate lag carriers.
#
# Grid: economic data DOSE, KUMMU (GADM1) and PWT, WB (GADM0) x climate source
# ERA5, CRU TS, UDelaware x weighting concurrent population, area (WCD
# "unweighted", i.e. the plain average over grid cells) = 24 panels.
#
# Samples:
#   max     each panel's own estimation sample;
#   common  subnational region-years estimable in all 12 subnational panels
#           and national country-years estimable in all 12 national panels,
#           both cut to the country-years present at both levels, so every
#           model covers the same countries and years.
#
# Output: results/signed_bins_across_datasets/. Part 1 estimates and writes
# the analysis tables; part 2 plots.

source("load_functions.R")

# Configuration ----------------------------------------------------------------

ECON_DATASETS <- c(DOSE = "gadm1", KUMMU = "gadm1", PWT = "gadm0", WB = "gadm0")
CLIMATE_SOURCES <- c(ERA5 = "ERA5", CRU = "CRU TS", UDEL = "UDelaware")
CLIMATE_WEIGHTS <- c(pop_concurrent = "concurrent population", area = "unweighted")

# Notebook defaults (test_functions.Rmd, shared-specification chunk).
OUTCOME <- "dlgrp_pc_usd"
FIXED_EFFECTS <- "year + GID_1[year]"
PANEL_ID <- c("GID_1", "year")
ECON_YEAR_MIN <- 1950L
CLIMATE_VELOCITY_YEARS <- 5L
BASE_CLIMATE <- "TM + TM:mean_TM_all + RR + RR:mean_RR_all"
MODEL_TERMS <- paste(
  BASE_CLIMATE,
  "+ i(TM_bin_signed, ref = 0)",
  "+ i(RR_bin_signed, ref = 0)"
)
MODEL_VARS <- c(OUTCOME, "TM", "RR", "mean_TM_all", "mean_RR_all",
                "TM_bin_signed", "RR_bin_signed")
KEEP_VARS <- c("GID_0", "GID_1", "year", MODEL_VARS, "zTM", "zRR")

OUTPUT_DIR <- file.path("results", "signed_bins_across_datasets")
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
  econ = names(ECON_DATASETS),
  climate = names(CLIMATE_SOURCES),
  weight = names(CLIMATE_WEIGHTS),
  stringsAsFactors = FALSE
) %>%
  mutate(
    level = unname(ECON_DATASETS[econ]),
    id = paste(econ, climate, weight, sep = "_")
  )

panels <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  message("Building ", s$id)
  cfg <- panel_config(
    econ_year_min = ECON_YEAR_MIN,
    climate_source = CLIMATE_SOURCES[[s$climate]],
    climate_weight = CLIMATE_WEIGHTS[[s$weight]],
    clim_history_years = max(6L, CLIMATE_VELOCITY_YEARS - 1L)
  )
  panels[[s$id]] <- build_dat(econ_data = s$econ, config = cfg) %>%
    select(all_of(KEEP_VARS))
}

# Maximum samples ----------------------------------------------------------------

fits_max <- lapply(panels, fit_signed_bins)

# Rows actually used by each maximum-sample model.
used_keys <- lapply(names(fits_max), function(id) {
  panels[[id]][fixest::obs(fits_max[[id]]), c("GID_0", "GID_1", "year")]
})
names(used_keys) <- names(fits_max)

# Common sample ------------------------------------------------------------------

common_within_level <- function(level) {
  intersect_keys(used_keys[specs$id[specs$level == level]])
}
common_sub <- common_within_level("gadm1")
common_nat <- common_within_level("gadm0")
common_country_years <- inner_join(
  common_sub %>% distinct(GID_0, year),
  common_nat %>% distinct(GID_0, year),
  by = c("GID_0", "year")
)
common_keys <- list(
  gadm1 = semi_join(common_sub, common_country_years, by = c("GID_0", "year")),
  gadm0 = semi_join(common_nat, common_country_years, by = c("GID_0", "year"))
)

common_data <- lapply(specs$id, function(id) {
  level <- specs$level[specs$id == id]
  semi_join(panels[[id]], common_keys[[level]], by = c("GID_0", "GID_1", "year"))
})
names(common_data) <- specs$id
fits_common <- lapply(common_data, fit_signed_bins)

# Part 1: analysis tables ----------------------------------------------------------
#
# Table builders are in functions/specification_grid.R.

FACTORS <- c("econ", "climate", "weight", "level")
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

# How much of the spread in each bin comes from the economic data, the climate
# source, the weighting and the geographic level.
variance_shares <- grid_dispersion(coefficients, FACTORS)
write_out(variance_shares, "cross_specification_dispersion")

factor_means <- grid_factor_means(coefficients, c("econ", "climate", "weight"))
write_out(factor_means, "mean_estimate_by_factor")

# Agreement of the climate data themselves on the common sample: correlation of
# the standardized anomalies and the share of observations falling in the same
# signed bin, for every pair of climate variants within an economic dataset.
climate_agreement <- bind_rows(lapply(names(ECON_DATASETS), function(econ) {
  ids <- specs$id[specs$econ == econ]
  common <- lapply(common_data[ids], function(d) {
    select(d, GID_1, year, zTM, zRR, TM_bin_signed, RR_bin_signed)
  })
  pairs <- utils::combn(ids, 2, simplify = FALSE)
  bind_rows(lapply(pairs, function(p) {
    j <- inner_join(common[[p[1]]], common[[p[2]]], by = c("GID_1", "year"))
    tibble::tibble(
      econ = econ,
      variant_1 = sub(paste0("^", econ, "_"), "", p[1]),
      variant_2 = sub(paste0("^", econ, "_"), "", p[2]),
      observations = nrow(j),
      cor_zTM = stats::cor(j$zTM.x, j$zTM.y, use = "complete.obs"),
      cor_zRR = stats::cor(j$zRR.x, j$zRR.y, use = "complete.obs"),
      same_TM_bin = mean(j$TM_bin_signed.x == j$TM_bin_signed.y, na.rm = TRUE),
      same_RR_bin = mean(j$RR_bin_signed.x == j$RR_bin_signed.y, na.rm = TRUE)
    )
  }))
}))
write_out(climate_agreement, "climate_data_agreement")

# WCD's ERA5 concurrent-population files store precipitation as 0 (not
# missing) for units without data; those have no temperature and never enter a
# model, but a zero next to a valid temperature would (Fiji 1993-2022).
data_quality <- grid_data_quality(all_fits, all_data, specs)
write_out(data_quality, "data_quality")

saveRDS(
  list(specs = specs, sample_summary = sample_summary, coefficients = coefficients,
       data_quality = data_quality,
       base_terms = base_terms, joint_tests = joint_tests,
       variance_shares = variance_shares, factor_means = factor_means,
       climate_agreement = climate_agreement),
  file.path(OUTPUT_DIR, "analysis.rds")
)

options(width = 200)
cat("\n== Samples ==\n")
print(sample_summary %>% select(sample, econ, climate, weight, observations, units,
                                countries, first_year, last_year, within_r2),
      n = Inf)
cat("\n== Joint tests (p-values) ==\n")
print(joint_tests %>% select(sample, econ, climate, weight, test, p_value) %>%
        tidyr::pivot_wider(names_from = test, values_from = p_value), n = Inf)
cat("\n== Cross-specification dispersion per bin ==\n")
print(variance_shares, n = Inf)
cat("\n== Climate data agreement (common sample) ==\n")
print(climate_agreement, n = Inf)
cat("\n== Data quality of estimation rows ==\n")
print(data_quality %>% select(sample, econ, climate, weight, rows, zero_precipitation,
                              RR_min_m, share_extreme_TM_bins, share_extreme_RR_bins),
      n = Inf)

# Part 2: plots ----------------------------------------------------------------------
#
# Plot builders are in functions/specification_grid.R: one colour and marker
# per climate source, one line type per weighting.

SAMPLE_LABELS <- c(
  max = "Maximum sample",
  common = "Common sample (1991-2017, 75 countries)"
)

plot_data <- grid_plot_data(
  coefficients,
  list(econ = names(ECON_DATASETS), climate = names(CLIMATE_SOURCES),
       weight = names(CLIMATE_WEIGHTS))
)

# Figure 1: every specification, one figure per sample.
for (s in names(SAMPLE_LABELS)) {
  p <- plot_bin_grid(
    plot_data %>% filter(sample == s),
    facet = "econ",
    colour = "climate",
    colour_name = "Climate data",
    linetype = "weight",
    linetype_values = c(pop_concurrent = "solid", area = "22"),
    linetype_labels = c(pop_concurrent = "Concurrent population", area = "Area"),
    linetype_name = "Weighting",
    title = paste("Signed-bin response by economic and climate data:",
                  tolower(SAMPLE_LABELS[[s]])),
    subtitle = paste("BHM base terms + signed half-SD anomaly bins;",
                     "reference bin -0.5 to 0.5 SD; 95% CI clustered by unit"),
    caption = paste("DOSE and KUMMU: GADM1 regions; PWT and WB: countries.",
                    "Fixed effects: year + unit-specific linear trends.")
  )
  ggsave(file.path(OUTPUT_DIR, paste0("signed_bins_", s, "_sample.png")), p,
         width = 12, height = 8, dpi = 300)
}

# Figure 2: the range across all 24 specifications for each bin, by sample.
ggsave(
  file.path(OUTPUT_DIR, "signed_bins_cross_specification_range.png"),
  plot_bin_range(
    plot_data,
    SAMPLE_LABELS,
    title = "How much the signed-bin response depends on the data",
    subtitle = paste("Point: mean over the 24 dataset x climate x weighting",
                     "specifications; thick bar: interquartile range; thin line: min-max")
  ),
  width = 11, height = 5.5, dpi = 300
)

# Figure 3: agreement of the climate data themselves on the common sample.
variant_levels <- as.vector(outer(names(CLIMATE_SOURCES), names(CLIMATE_WEIGHTS),
                                  paste, sep = "_"))
variant_label <- function(x) {
  sub("_pop_concurrent$", " pop.", sub("_area$", " area", x))
}
agreement_data <- climate_agreement %>%
  filter(econ %in% c("DOSE", "PWT")) %>%
  mutate(level = if_else(econ == "DOSE", "Regions (GADM1)", "Countries (GADM0)")) %>%
  tidyr::pivot_longer(c(same_TM_bin, same_RR_bin), names_to = "variable",
                      values_to = "value") %>%
  mutate(
    variable = if_else(variable == "same_TM_bin", "Temperature", "Precipitation"),
    variable = factor(variable, levels = c("Temperature", "Precipitation")),
    level = factor(level, levels = c("Regions (GADM1)", "Countries (GADM0)")),
    variant_1 = factor(variant_label(variant_1),
                       levels = variant_label(variant_levels)),
    variant_2 = factor(variant_label(variant_2),
                       levels = variant_label(variant_levels))
  )
ggsave(
  file.path(OUTPUT_DIR, "climate_data_agreement.png"),
  plot_agreement_heatmap(
    agreement_data,
    facet_rows = "variable",
    facet_cols = "level",
    fill_name = "Same signed bin",
    title = "Climate datasets often put the same year in different anomaly bins",
    subtitle = paste("Share of common-sample observations assigned to the same",
                     "signed half-SD bin by two climate variants"),
    limits = c(0.3, 1)
  ),
  width = 10, height = 8, dpi = 300
)

message("Figures written to ", normalizePath(OUTPUT_DIR))
