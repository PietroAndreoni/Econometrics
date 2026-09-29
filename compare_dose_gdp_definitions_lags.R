# Distributed-lag climate response across DOSE GDP-per-capita definitions and
# climate datasets: the grid of compare_dose_gdp_definitions.R applied to
#
#   l(dTM, 0:10) + l(dTM, 0:10):mean_TM_all + l(dRR, 0:10) + l(dRR, 0:10):mean_RR_all
#
# dTM and dRR are the year-on-year changes in annual temperature and
# precipitation. The effect of a change k years ago is b_k + g_k * mean_TM_all,
# so every lag path is evaluated at the 10th, 50th and 90th percentile of the
# long-run climate across regions (distributed_lag_effects()); the cumulative
# path sums the lags and keeps their covariances.
#
# Settings as in compare_dose_gdp_definitions.R: fixed effects year + GID_1 +
# GID_1[year] + GID_0^year, errors clustered by GID_1, econ_year_min = 1950.
# Eleven years of
# climate are kept before 1950, so the ten lags of a change are available from
# the first estimation year.
#
# Grid: the five DOSE definitions (DOSE_GDP_VARIABLES) x ERA5, CRU TS,
# UDelaware, all with concurrent-population weights = 15 panels.
# Samples: max (each panel's own) and common (region-years used by all 15).
#
# Output: results/distributed_lags_dose_gdp_definitions/. Part 1 estimates and
# writes the analysis tables; part 2 plots, following the distributed-lag
# template of main_analysis.Rmd (plot_distributed_lags()).

source("load_functions.R")

# Configuration ----------------------------------------------------------------

DEFINITIONS <- c(
  lcu2015_usd = "Const. LCU, 2015 FX",
  lcu_2015 = "Const. LCU",
  usd_2015 = "US$ / US deflator",
  lcu = "Nominal LCU",
  usd = "Nominal US$"
)
DEFAULT_DEFINITION <- "Const. LCU, 2015 FX"
CLIMATE_SOURCES <- c(ERA5 = "ERA5", CRU = "CRU TS", UDEL = "UDelaware")
CLIMATE_WEIGHT <- "concurrent population"

OUTCOME <- "dlgrp_pc_usd"
FIXED_EFFECTS <- "year + GID_1 + GID_1[year] + GID_0^year"
PANEL_ID <- c("GID_1", "year")
ECON_YEAR_MIN <- 1950L
MAX_LAG <- 10L
CLIM_HISTORY_YEARS <- MAX_LAG + 1L
MODEL_TERMS <- sprintf(
  paste("l(dTM, 0:%1$d) + l(dTM, 0:%1$d):mean_TM_all",
        "+ l(dRR, 0:%1$d) + l(dRR, 0:%1$d):mean_RR_all"),
  MAX_LAG
)
INTERACTIONS <- c(dTM = "mean_TM_all", dRR = "mean_RR_all")
PERCENTILES <- c(0.1, 0.5, 0.9)
CONF_LEVEL <- 0.95
KEEP_VARS <- c("GID_0", "GID_1", "year", OUTCOME, "dTM", "dRR", "TM", "RR",
               "mean_TM_all", "mean_RR_all")

OUTPUT_DIR <- file.path("results", "distributed_lags_dose_gdp_definitions")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

fit_lags <- function(model_data) {
  fit_panel_model(MODEL_TERMS, model_data, outcome = OUTCOME,
                  fixed_effects = FIXED_EFFECTS, panel_id = PANEL_ID,
                  cluster = ~GID_1)
}

# Panels -------------------------------------------------------------------------

specs <- expand.grid(
  definition_code = names(DEFINITIONS),
  climate = names(CLIMATE_SOURCES),
  stringsAsFactors = FALSE
) %>%
  mutate(definition = unname(DEFINITIONS[definition_code]),
         id = paste(definition_code, climate, sep = "_"))

panels <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  message("Building ", s$id)
  cfg <- panel_config(econ_year_min = ECON_YEAR_MIN,
                      climate_source = CLIMATE_SOURCES[[s$climate]],
                      climate_weight = CLIMATE_WEIGHT,
                      clim_history_years = CLIM_HISTORY_YEARS)
  panels[[s$id]] <- suppressMessages(
    build_dat("DOSE", config = cfg, econ_variable = s$definition_code)
  ) %>%
    select(all_of(KEEP_VARS))
}

# Samples --------------------------------------------------------------------------
#
# The lag carriers (years without an outcome) must stay in the data for l() to
# resolve, so the common sample keeps every climate row and withholds only the
# outcome of region-years outside it.

fits_max <- lapply(panels, fit_lags)
used_keys <- lapply(names(fits_max), function(id) {
  panels[[id]][fixest::obs(fits_max[[id]]), KEY_COLS]
})
common_keys <- intersect_keys(used_keys) %>% mutate(.common = TRUE)
common_data <- lapply(panels, function(d) {
  d %>%
    left_join(common_keys, by = KEY_COLS) %>%
    mutate(!!OUTCOME := if_else(coalesce(.common, FALSE), .data[[OUTCOME]],
                                NA_real_)) %>%
    select(-.common)
})
fits_common <- lapply(common_data, fit_lags)

all_fits <- list(max = fits_max, common = fits_common)
all_data <- list(max = panels, common = common_data)

# Part 1: analysis tables ----------------------------------------------------------

sample_summary <- grid_sample_summary(all_fits, all_data, specs)
write_out(sample_summary, "sample_summary")

# Lag paths at the 10th/50th/90th percentile of each panel's long-run climate
# across the regions it estimates on.
lag_effects <- bind_rows(lapply(names(all_fits), function(sample) {
  bind_rows(lapply(specs$id, function(id) {
    used <- grid_used_rows(all_fits, all_data, sample, id)
    at <- list(
      mean_TM_all = region_quantiles(used, "mean_TM_all", PERCENTILES),
      mean_RR_all = region_quantiles(used, "mean_RR_all", PERCENTILES)
    )
    distributed_lag_effects(all_fits[[sample]][[id]], CONF_LEVEL,
                            interactions = INTERACTIONS, at = at) %>%
      mutate(sample = sample, id = id, .before = 1)
  }))
})) %>%
  left_join(specs, by = "id") %>%
  mutate(
    variable = if_else(expression == "dTM", "Temperature change (dTM)",
                       "Precipitation change (dRR)"),
    p_value = 2 * stats::pnorm(-abs(effect / standard_error))
  )
write_out(lag_effects, "lag_effects")

# Headline numbers: the contemporaneous effect and the cumulative effect after
# ten years, at the median climate and at the 10th/90th percentiles.
headline <- lag_effects %>%
  filter((horizon == "Effect at lag" & lag == 0) |
           (horizon == "Cumulative through lag" & lag == MAX_LAG)) %>%
  mutate(quantity = if_else(lag == 0, "effect_lag0", "cumulative_10")) %>%
  select(sample, definition, climate, variable, at, at_value, quantity, effect,
         standard_error, p_value)
write_out(headline, "headline_effects")

# Joint tests: all lags of a variable (and their interactions) are zero; and,
# on the cumulative scale, no persistence beyond the contemporaneous year.
joint_tests <- bind_rows(lapply(names(all_fits), function(sample) {
  bind_rows(lapply(specs$id, function(id) {
    fit <- all_fits[[sample]][[id]]
    bind_rows(
      wald_row(fit, "^l\\(dTM", label = "All dTM lags and interactions = 0"),
      wald_row(fit, "^l\\(dTM, ([1-9]|10)\\)",
               label = "dTM lags 1-10 and interactions = 0"),
      wald_row(fit, "^l\\(dRR", label = "All dRR lags and interactions = 0"),
      wald_row(fit, "^l\\(dRR, ([1-9]|10)\\)",
               label = "dRR lags 1-10 and interactions = 0")
    ) %>%
      mutate(sample = sample, id = id, .before = 1)
  }))
})) %>%
  left_join(specs, by = "id")
write_out(joint_tests, "joint_tests")

# Spread of the ten-year cumulative effect at the median climate across the 15
# specifications, and the share explained by the definition and by the climate
# source.
cumulative_dispersion <- headline %>%
  filter(quantity == "cumulative_10", at == "p50") %>%
  group_by(sample, variable) %>%
  group_modify(function(d, key) {
    total <- sum((d$effect - mean(d$effect))^2)
    share <- function(f) {
      sum((stats::ave(d$effect, d[[f]]) - mean(d$effect))^2) / total
    }
    tibble::tibble(
      mean = mean(d$effect), sd = stats::sd(d$effect),
      min = min(d$effect), max = max(d$effect),
      share_negative = mean(d$effect < 0),
      share_significant = mean(d$p_value < 0.05),
      share_definition = share("definition"),
      share_climate = share("climate")
    )
  }) %>%
  ungroup()
write_out(cumulative_dispersion, "cumulative_dispersion")

saveRDS(
  list(specs = specs, sample_summary = sample_summary, lag_effects = lag_effects,
       headline = headline, joint_tests = joint_tests,
       cumulative_dispersion = cumulative_dispersion),
  file.path(OUTPUT_DIR, "analysis.rds")
)

options(width = 200)
cat("\n== Samples ==\n")
print(sample_summary %>% select(sample, definition, climate, observations, units,
                                countries, first_year, last_year, within_r2),
      n = Inf)
cat("\n== Ten-year cumulative effect at the median climate ==\n")
print(headline %>%
        filter(quantity == "cumulative_10", at == "p50") %>%
        mutate(across(c(effect, standard_error), ~ round(.x, 4)),
               p_value = signif(p_value, 2)) %>%
        select(sample, variable, definition, climate, effect, standard_error, p_value),
      n = Inf)
cat("\n== Ten-year cumulative effect by long-run climate (default definition) ==\n")
print(headline %>%
        filter(quantity == "cumulative_10", definition == DEFAULT_DEFINITION) %>%
        mutate(across(c(at_value, effect, standard_error), ~ round(.x, 4)),
               p_value = signif(p_value, 2)) %>%
        select(sample, variable, climate, at, at_value, effect, standard_error,
               p_value),
      n = Inf)
cat("\n== Joint tests (p-values) ==\n")
print(joint_tests %>% select(sample, definition, climate, test, p_value) %>%
        tidyr::pivot_wider(names_from = test, values_from = p_value), n = Inf)
cat("\n== Dispersion of the ten-year cumulative effect ==\n")
print(cumulative_dispersion, n = Inf)

# Part 2: plots ----------------------------------------------------------------------
#
# Every figure follows the main_analysis.Rmd distributed-lag template
# (plot_distributed_lags()): lag paths with interval ribbons, faceted by
# variable and horizon or by a data choice.

SAMPLE_LABELS <- c(max = "maximum sample", common = "common sample")
plot_data <- lag_effects %>%
  mutate(
    expression = factor(variable, levels = c("Temperature change (dTM)",
                                             "Precipitation change (dRR)")),
    definition = factor(definition, levels = unname(DEFINITIONS)),
    climate = factor(climate, levels = names(CLIMATE_SOURCES)),
    at = factor(at, levels = paste0("p", round(100 * PERCENTILES)),
                labels = paste0(round(100 * PERCENTILES), "th percentile"))
  )
Y_LABEL <- "Effect on log GDP per-capita growth"

for (s in names(SAMPLE_LABELS)) {
  d <- plot_data %>% filter(sample == s)

  # Figure 1: the template (effect at lag and cumulative), default definition,
  # one colour per climate source, at the median long-run climate.
  p1 <- plot_distributed_lags(
    d %>% filter(definition == DEFAULT_DEFINITION, at == "50th percentile"),
    CONF_LEVEL, Y_LABEL, colour = "climate", colour_name = "Climate data",
    title = paste0("Distributed-lag paths, DOSE ", DEFAULT_DEFINITION, ": ",
                   SAMPLE_LABELS[[s]]),
    subtitle = "Evaluated at the median long-run climate across regions"
  )
  ggsave(file.path(OUTPUT_DIR, paste0("lag_paths_default_", s, ".png")), p1,
         width = 10, height = 6.5, dpi = 300)

  # Figure 2: cumulative paths across the DOSE definitions.
  p2 <- plot_distributed_lags(
    d %>% filter(at == "50th percentile"),
    CONF_LEVEL, Y_LABEL, colour = "climate", colour_name = "Climate data",
    horizons = "Cumulative through lag",
    facets = expression ~ definition,
    title = paste0("Cumulative lag effect by DOSE GDP definition: ",
                   SAMPLE_LABELS[[s]]),
    subtitle = paste("Cumulative effect through each lag, at the median long-run",
                     "climate; Const. LCU and Const. LCU, 2015 FX coincide, as do",
                     "the two US$ series")
  )
  ggsave(file.path(OUTPUT_DIR, paste0("cumulative_by_definition_", s, ".png")), p2,
         width = 14, height = 6.5, dpi = 300)

  # Figure 3: cumulative paths by long-run climate, default definition.
  p3 <- plot_distributed_lags(
    d %>% filter(definition == DEFAULT_DEFINITION),
    CONF_LEVEL, Y_LABEL, colour = "at", colour_name = "Long-run climate",
    horizons = "Cumulative through lag",
    facets = expression ~ climate,
    title = paste0("Cumulative lag effect by long-run climate, DOSE ",
                   DEFAULT_DEFINITION, ": ", SAMPLE_LABELS[[s]]),
    subtitle = paste("Percentiles across regions of mean_TM_all (temperature)",
                     "and mean_RR_all (precipitation), per climate dataset")
  )
  ggsave(file.path(OUTPUT_DIR, paste0("cumulative_by_climate_", s, ".png")), p3,
         width = 11, height = 6.5, dpi = 300)
}

message("Figures written to ", normalizePath(OUTPUT_DIR))
