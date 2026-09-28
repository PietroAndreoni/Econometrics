# Split the local temperature anomaly into a global component and an orthogonal
# local component.
#
# For each region i the anomaly-space pattern regression
#
#   zTM_it = c_i + phi_i * zG_t + u_it
#
# is estimated on the full ERA5 record, and the residual
#
#   zTM_orth_it = u_it = zTM_it - c_i - phi_i * zG_t
#
# is the part of the local temperature anomaly the global signal does not
# explain. `phi_i` is the anomaly-space analogue of a pattern-scaling
# coefficient: how strongly region i's weather anomaly loads on the global one.
#
# Both zTM and zG are built with the same operator - deviation from the
# *preceding* year's 20-year rolling mean, divided by the preceding 20-year
# rolling SD - so `zTM_orth` is on `zTM`'s own scale and drops into the existing
# model strings without rescaling.
#
# `phi_i` is estimated on the whole climate record rather than on the GDP
# estimation window, which sharpens it and makes it largely predetermined
# relative to the outcome (see the generated-regressor discussion in
# global_local_temperature_models.R).

suppressPackageStartupMessages(library(dplyr))

required_packages <- c("arrow", "ggplot2", "tibble", "tidyr", "tidyselect")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (file.exists("Econometrics.Rproj")) {
  project_root <- "."
} else if (file.exists(file.path("..", "Econometrics.Rproj"))) {
  project_root <- ".."
} else {
  stop("Run this script from the project root or the econometrics directory.")
}

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "global_local_temperature")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Configuration ----------------------------------------------------------------

# Matches build_dat()'s gadm1 selection in test_functions.Rmd: era5, population
# weights, base year 2015, and the 20-year baseline behind zTM.
SHOCK_CLIMATE_SOURCE <- "era5"
SHOCK_WEIGHT <- "pop"
SHOCK_WEIGHT_YEAR <- "2015"
SHOCK_WINDOW <- 20L
SHOCK_MIN_YEARS <- 20L

global_file <- file.path(data_dir, "data_gm_gadmworld_era5_pop-area_2015.parquet")
local_file <- file.path(
  data_dir, "data_tm_gadm1_era5_pop-area_2000-2015.parquet"
)
missing_inputs <- c(global_file, local_file)[
  !file.exists(c(global_file, local_file))
]
if (length(missing_inputs)) {
  stop(
    "Input files not found: ", paste(missing_inputs, collapse = ", "),
    ". Run prepare_global_climate_data.R first."
  )
}

# Local anomalies ---------------------------------------------------------------

# Calendar-year lag, as in test_functions.Rmd: a gap in the series yields NA
# rather than a silently wrong reference year.
lag_by_year <- function(x, year, k = 1) {
  x[match(year - k, year)]
}

local_anomalies <- arrow::read_parquet(
  local_file,
  col_select = tidyselect::all_of(c(
    "gadm_level", "climate_source", "window_tm", "weight", "weight_year",
    "GID_0", "GID_1", "year", "TM", "mean_TM", "sd_TM", "mean_TM_all"
  ))
) %>%
  tibble::as_tibble() %>%
  filter(
    window_tm == SHOCK_WINDOW,
    weight == SHOCK_WEIGHT,
    weight_year == SHOCK_WEIGHT_YEAR,
    climate_source == SHOCK_CLIMATE_SOURCE,
    !is.na(TM)
  ) %>%
  distinct() %>%
  group_by(GID_0, GID_1) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    mean_TM_lag = lag_by_year(mean_TM, year),
    sd_TM_lag = lag_by_year(sd_TM, year),
    zTM = if_else(
      is.finite(sd_TM_lag) & sd_TM_lag > 0,
      (TM - mean_TM_lag) / sd_TM_lag,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  select(GID_0, GID_1, year, TM, mean_TM_all, zTM)

# Global anomalies --------------------------------------------------------------

global_anomalies <- arrow::read_parquet(global_file) %>%
  tibble::as_tibble() %>%
  filter(window == SHOCK_WINDOW, climate_source == SHOCK_CLIMATE_SOURCE) %>%
  transmute(
    global_weight = weight,
    year = as.integer(year),
    zG,
    Gcyc,
    Gtrend,
    G
  )

# Pattern loading ---------------------------------------------------------------

# Closed-form OLS of zTM on zG within region: phi = cov / var, and the usual
# homoskedastic standard error. Regions with too few overlapping years are kept
# but flagged rather than dropped, so the estimation sample is not silently
# altered.
#
# `lambda` is the companion regression in *levels*, TM_it = a_i + lambda_i * G_t,
# i.e. classical pattern scaling: degrees of local warming per degree of global
# warming. It is a different object from `phi` and must not be confused with it.
# Because both sides of the phi regression are standardised by their own SD, phi
# is a correlation-like loading bounded near one; lambda is unbounded and
# averages above one over land. phi drives the orthogonalisation (it is what
# makes zTM_orth live on zTM's scale); lambda is the physical exposure measure
# used by the heterogeneous-exposure specification.
estimate_pattern_loading <- function(paired) {
  paired %>%
    filter(is.finite(zTM), is.finite(zG)) %>%
    group_by(global_weight, GID_0, GID_1) %>%
    summarise(
      n_years = n(),
      first_year = min(year),
      last_year = max(year),
      var_zG = stats::var(zG),
      var_zTM = stats::var(zTM),
      phi = stats::cov(zTM, zG) / var_zG,
      intercept = mean(zTM) - phi * mean(zG),
      r2 = stats::cor(zTM, zG)^2,
      lambda = stats::cov(TM, G) / stats::var(G),
      lambda_r2 = stats::cor(TM, G)^2,
      mean_TM_all = first(mean_TM_all),
      .groups = "drop"
    ) %>%
    mutate(
      se_phi = sqrt(
        pmax((1 - r2), 0) * var_zTM / (pmax(n_years - 2L, 1L) * var_zG)
      ),
      estimable = n_years >= SHOCK_MIN_YEARS &
        is.finite(phi) &
        is.finite(se_phi) &
        var_zG > 0
    )
}

# Empirical-Bayes shrinkage toward the precision-weighted mean (Morris 1983).
# Reported as a robustness variant, not the headline: with 60+ annual
# observations per region the shrinkage weight is close to one almost everywhere.
add_shrinkage <- function(loading) {
  loading %>%
    group_by(global_weight) %>%
    mutate(
      phi_bar = stats::weighted.mean(
        phi[estimable], w = 1 / se_phi[estimable]^2
      ),
      tau2 = pmax(
        0,
        stats::var(phi[estimable]) - mean(se_phi[estimable]^2)
      ),
      shrinkage_weight = tau2 / (tau2 + se_phi^2),
      phi_eb = shrinkage_weight * phi + (1 - shrinkage_weight) * phi_bar
    ) %>%
    ungroup()
}

paired <- local_anomalies %>%
  inner_join(global_anomalies, by = "year", relationship = "many-to-many")

pattern_loading <- paired %>%
  estimate_pattern_loading() %>%
  add_shrinkage()

# Orthogonal local shock --------------------------------------------------------

local_shocks <- paired %>%
  inner_join(
    pattern_loading %>%
      select(
        global_weight, GID_0, GID_1, phi, intercept, phi_eb, lambda, estimable
      ),
    by = c("global_weight", "GID_0", "GID_1")
  ) %>%
  mutate(
    zTM_orth = zTM - intercept - phi * zG,
    zTM_orth_eb = zTM - intercept - phi_eb * zG
  ) %>%
  select(
    global_weight, GID_0, GID_1, year,
    zTM, zG, Gcyc, Gtrend, G,
    zTM_orth, zTM_orth_eb, phi, phi_eb, lambda, estimable
  ) %>%
  arrange(global_weight, GID_0, GID_1, year)

arrow::write_parquet(
  pattern_loading,
  file.path(data_dir, "pattern_loading_gadm1_era5_pop_2015.parquet")
)
arrow::write_parquet(
  local_shocks,
  file.path(data_dir, "data_local_shock_gadm1_era5_pop_2015.parquet")
)

# Diagnostics -------------------------------------------------------------------

# Within a region, zTM_orth is orthogonal to zG by OLS construction, so the
# residual correlation must be zero to machine precision. It is legitimately
# non-zero for the shrunk variant, which is exactly why that one is a robustness
# row: a non-zero residual correlation means the global coefficient is no longer
# a clean total effect.
orthogonality_checks <- local_shocks %>%
  filter(estimable, is.finite(zTM_orth), is.finite(zG)) %>%
  group_by(global_weight, GID_1) %>%
  filter(n() >= 3L, stats::sd(zG) > 0) %>%
  summarise(
    cor_orth = stats::cor(zTM_orth, zG),
    cor_orth_eb = stats::cor(zTM_orth_eb, zG),
    cor_raw = stats::cor(zTM, zG),
    .groups = "drop"
  ) %>%
  group_by(global_weight) %>%
  summarise(
    regions = n(),
    max_abs_cor_orthogonal = max(abs(cor_orth)),
    max_abs_cor_orthogonal_eb = max(abs(cor_orth_eb)),
    mean_abs_cor_orthogonal_eb = mean(abs(cor_orth_eb)),
    mean_cor_raw_zTM_zG = mean(cor_raw),
    .groups = "drop"
  )

pattern_loading_summary <- pattern_loading %>%
  filter(estimable) %>%
  group_by(global_weight) %>%
  summarise(
    regions = n(),
    countries = n_distinct(GID_0),
    mean_phi = mean(phi),
    median_phi = stats::median(phi),
    sd_phi = stats::sd(phi),
    p05_phi = unname(stats::quantile(phi, 0.05)),
    p95_phi = unname(stats::quantile(phi, 0.95)),
    share_phi_above_1 = mean(phi > 1),
    median_r2 = stats::median(r2),
    share_r2_below_02 = mean(r2 < 0.2),
    cor_phi_mean_TM = stats::cor(phi, mean_TM_all, use = "complete.obs"),
    mean_shrinkage_weight = mean(shrinkage_weight),
    # Levels pattern scaling, reported alongside for interpretation.
    mean_lambda = mean(lambda),
    median_lambda = stats::median(lambda),
    share_lambda_above_1 = mean(lambda > 1),
    cor_lambda_mean_TM = stats::cor(lambda, mean_TM_all, use = "complete.obs"),
    cor_phi_lambda = stats::cor(phi, lambda),
    .groups = "drop"
  )

variance_shares <- local_shocks %>%
  filter(estimable, is.finite(zTM), is.finite(zTM_orth)) %>%
  group_by(global_weight) %>%
  summarise(
    sd_zTM = stats::sd(zTM),
    sd_zTM_orth = stats::sd(zTM_orth),
    variance_share_explained_by_global = 1 -
      stats::var(zTM_orth) / stats::var(zTM),
    .groups = "drop"
  )

write.csv(
  pattern_loading,
  file.path(output_dir, "pattern_loading.csv"),
  row.names = FALSE
)
write.csv(
  pattern_loading_summary,
  file.path(output_dir, "pattern_loading_summary.csv"),
  row.names = FALSE
)
write.csv(
  orthogonality_checks,
  file.path(output_dir, "orthogonality_checks.csv"),
  row.names = FALSE
)
write.csv(
  variance_shares,
  file.path(output_dir, "global_variance_share.csv"),
  row.names = FALSE
)

phi_plot_data <- pattern_loading %>%
  filter(estimable, global_weight == "pop")

phi_distribution_plot <- ggplot2::ggplot(
  phi_plot_data,
  ggplot2::aes(x = phi)
) +
  ggplot2::geom_histogram(bins = 60, fill = "#0072B2", colour = NA) +
  ggplot2::geom_vline(
    xintercept = 1,
    colour = "black",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  ggplot2::labs(
    title = "Anomaly-space pattern loading by region",
    subtitle = paste0(
      "phi_i from zTM_it = c_i + phi_i * zG_t, ERA5 population-weighted, ",
      nrow(phi_plot_data), " GADM1 regions. Dashed line: phi = 1."
    ),
    x = "phi (local anomaly SD per unit global anomaly SD)",
    y = "Regions"
  ) +
  ggplot2::theme_classic()

# The two loadings move in opposite directions with baseline climate, and the
# contrast is the point: in standardised-anomaly space warm regions track the
# global signal most closely (low idiosyncratic variability), while in levels it
# is cold regions that warm fastest per degree of global warming.
phi_climate_plot <- phi_plot_data %>%
  select(mean_TM_all, phi, lambda) %>%
  tidyr::pivot_longer(c(phi, lambda), names_to = "loading", values_to = "value") %>%
  mutate(loading = if_else(
    loading == "phi",
    "phi: zTM on zG (standardised anomalies)",
    "lambda: TM on G (levels, °C per °C)"
  )) %>%
  ggplot2::ggplot(ggplot2::aes(x = mean_TM_all, y = value)) +
  ggplot2::geom_hline(yintercept = 1, colour = "black", linewidth = 0.4) +
  ggplot2::geom_point(alpha = 0.2, size = 0.6, colour = "#0072B2") +
  ggplot2::geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    colour = "#E69F00",
    linewidth = 0.8
  ) +
  ggplot2::facet_wrap(~loading, scales = "free_y") +
  ggplot2::labs(
    title = "Two loadings on the global signal, and why they differ",
    subtitle = paste(
      "Each point is a GADM1 region; long-run mean temperature 1990-2019.",
      "Reference line at 1."
    ),
    x = "Mean temperature (°C)",
    y = NULL
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  file.path(output_dir, "phi_distribution.png"),
  phi_distribution_plot,
  width = 8,
  height = 4.5,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  file.path(output_dir, "phi_vs_mean_temperature.png"),
  phi_climate_plot,
  width = 10,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

cat("\nPattern loading (zTM on zG):\n")
print(as.data.frame(pattern_loading_summary), digits = 4)
cat("\nShare of local anomaly variance explained by the global anomaly:\n")
print(as.data.frame(variance_shares), digits = 4)
cat("\nOrthogonality of the local shock to the global shock:\n")
print(as.data.frame(orthogonality_checks), digits = 4)
cat(
  "\nHighest and lowest pattern loadings (population-weighted global series):\n"
)
print(
  phi_plot_data %>%
    arrange(desc(phi)) %>%
    select(GID_0, GID_1, phi, lambda, r2, mean_TM_all, n_years) %>%
    slice(c(1:8, (n() - 7):n())) %>%
    as.data.frame(),
  digits = 3
)
cat("\nHighest and lowest levels pattern scaling (lambda, °C per °C):\n")
print(
  phi_plot_data %>%
    arrange(desc(lambda)) %>%
    select(GID_0, GID_1, lambda, phi, lambda_r2, mean_TM_all) %>%
    slice(c(1:8, (n() - 7):n())) %>%
    as.data.frame(),
  digits = 3
)
cat("\nWritten to: ", normalizePath(output_dir), "\n", sep = "")
