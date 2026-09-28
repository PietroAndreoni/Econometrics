# Global temperature shocks and orthogonal local shocks in the DOSE subnational
# growth panel.
#
# The notebook specification in test_functions.Rmd is
#
#   dlgrp_pc_usd ~ TM + TM^2 + RR + RR^2 + zTM^2 + zRR^2
#                  | year + GID_1[year] + GID_0[year]
#
# The `year` fixed effect absorbs everything common to all regions in a year, so
# the effect of global temperature is set to zero by assumption rather than
# estimated. This script relaxes that, using the global series built by
# prepare_global_climate_data.R and the orthogonal local shock built by
# estimate_global_local_shocks.R, and compares the three ways of identifying a
# global coefficient in a panel of this shape:
#
#   Model 1  drop the year fixed effect; identify from time variation in the one
#            global series (Bilal & Kanzig 2024)
#   Model 2  keep the year fixed effect; identify from heterogeneous exposure by
#            interacting the global signal with the region's pattern scaling
#            (the shift-share device in Kotz, Levermann & Wenz 2024)
#   Model 3  two step: take the year effects out of the notebook specification
#            and regress them on the global shock as a short time series
#
# plus a fixed-effect ladder that quantifies how much identifying variation each
# trend control destroys.
#
# A specification point that decides what the global coefficient means. `TM` is
# the region's annual temperature *level*, so it carries the global warming
# signal itself. Conditioning on it makes the global coefficient a partial
# (pure spillover) effect. The primary local block therefore lets temperature
# enter only through the orthogonalised anomaly, which is what makes the global
# coefficient the *total* effect of warming - the object the orthogonalisation
# was constructed to deliver. The conditional block is reported alongside so the
# two estimands can be compared rather than conflated.

suppressPackageStartupMessages(library(dplyr))

required_packages <- c(
  "arrow", "fixest", "ggplot2", "lmtest", "purrr", "sandwich", "tibble",
  "tidyr", "tidyselect"
)
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

source(file.path(project_root, "rmd_chunks.R"))

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "global_local_temperature")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Configuration ----------------------------------------------------------------

GLOBAL_WEIGHT <- "pop"          # "pop" matches the panel's own weighting
GLOBAL_ECON_DATA <- "DOSE_V2_11"
CONFIDENCE_LEVEL <- 0.95
LP_HORIZONS <- 0:10
# Bootstrap sizes are overridable from the environment so the script can be smoke
# tested cheaply without editing it.
integer_option <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
  if (is.na(value)) default else value
}
BOOTSTRAP_REPLICATIONS <- integer_option("GLOBAL_LOCAL_REGION_BOOT", 200L)
YEAR_BOOTSTRAP_REPLICATIONS <- integer_option("GLOBAL_LOCAL_YEAR_BOOT", 999L)
BOOTSTRAP_SEED <- 20260907L

test_functions_rmd <- file.path(project_root, "test_functions.Rmd")
global_file <- file.path(data_dir, "data_gm_gadmworld_era5_pop-area_2015.parquet")
shock_file <- file.path(data_dir, "data_local_shock_gadm1_era5_pop_2015.parquet")
loading_file <- file.path(data_dir, "pattern_loading_gadm1_era5_pop_2015.parquet")

missing_inputs <- c(test_functions_rmd, global_file, shock_file, loading_file)[
  !file.exists(c(test_functions_rmd, global_file, shock_file, loading_file))
]
if (length(missing_inputs)) {
  stop(
    "Input files not found: ", paste(missing_inputs, collapse = ", "),
    ". Run prepare_global_climate_data.R and estimate_global_local_shocks.R first."
  )
}

# Panel ------------------------------------------------------------------------

notebook <- load_notebook_env(test_functions_rmd, econ_data = GLOBAL_ECON_DATA)

baseline_panel <- notebook$build_dat(
  econ_data = GLOBAL_ECON_DATA,
  econ_files_by_level = c(
    gadm0 = file.path(data_dir, "econ_processed_WDI-WB-PWT110.parquet"),
    gadm1 = file.path(
      data_dir,
      "econ_processed_DOSE_V2_11-KUMMU2018-KUMMU2025-WDI-WB.parquet"
    )
  ),
  climate_rr_files_by_level = c(
    gadm0 = file.path(data_dir, "data_rr_gadm0_era5_pop-area_2000-2015.parquet"),
    gadm1 = file.path(data_dir, "data_rr_gadm1_era5_pop-area_2000-2015.parquet")
  ),
  climate_tm_files_by_level = c(
    gadm0 = file.path(data_dir, "data_tm_gadm0_era5_pop-area_2000-2015.parquet"),
    gadm1 = file.path(data_dir, "data_tm_gadm1_era5_pop-area_2000-2015.parquet")
  )
) %>%
  mutate(year = as.integer(year))

local_shocks <- arrow::read_parquet(shock_file) %>%
  tibble::as_tibble() %>%
  filter(global_weight == GLOBAL_WEIGHT) %>%
  transmute(
    GID_0, GID_1,
    year = as.integer(year),
    zTM_check = zTM,
    zTM_orth, zTM_orth_eb, zG, Gcyc, Gtrend, G
  )

pattern_loading <- arrow::read_parquet(loading_file) %>%
  tibble::as_tibble() %>%
  filter(global_weight == GLOBAL_WEIGHT, estimable) %>%
  select(GID_0, GID_1, phi, lambda, r2, n_years)

panel <- baseline_panel %>%
  left_join(local_shocks, by = c("GID_0", "GID_1", "year")) %>%
  left_join(pattern_loading, by = c("GID_0", "GID_1")) %>%
  mutate(
    # Centred exposure measures. `lambda` is levels pattern scaling (degrees of
    # local warming per degree of global warming) and is the physical exposure
    # measure; `phi` is the standardised-anomaly loading that defines the
    # orthogonalisation. They are different objects and are only weakly related.
    lambda_c = lambda - mean(lambda[!is.na(dlgrp_pc_usd)], na.rm = TRUE),
    phi_c = phi - mean(phi[!is.na(dlgrp_pc_usd)], na.rm = TRUE)
  )

if (nrow(panel) != nrow(baseline_panel)) {
  stop(
    "Joining the global and local-shock tables changed the panel row count (",
    nrow(baseline_panel), " -> ", nrow(panel),
    "). Check for duplicate keys."
  )
}

estimation_panel <- panel %>% filter(!is.na(dlgrp_pc_usd))

# Gates ------------------------------------------------------------------------

# The zTM recomputed in estimate_global_local_shocks.R must equal build_dat()'s
# own zTM; if it does not, the orthogonalisation is being applied to a different
# anomaly than the notebook estimates on.
zTM_gate <- estimation_panel %>%
  filter(is.finite(zTM), is.finite(zTM_check))
max_zTM_gap <- max(abs(zTM_gate$zTM - zTM_gate$zTM_check))
if (!is.finite(max_zTM_gap) || max_zTM_gap > 1e-9) {
  stop(
    "Recomputed zTM does not match build_dat()'s zTM (max |gap| = ",
    format(max_zTM_gap), ")."
  )
}

# The notebook specification must be untouched by the joins.
model_baseline <- "TM + TM^2 + RR + RR^2 + zTM^2 + zRR^2"
m_baseline_before <- notebook$fit_climate_model(
  model_baseline,
  model_data = baseline_panel
)
m_baseline_after <- notebook$fit_climate_model(
  model_baseline,
  model_data = panel
)
if (!isTRUE(all.equal(
  coef(m_baseline_before), coef(m_baseline_after), tolerance = 1e-10
))) {
  stop("The joins changed the baseline coefficients.")
}
if (!identical(stats::nobs(m_baseline_before), stats::nobs(m_baseline_after))) {
  stop("The joins changed the baseline sample size.")
}

# Term blocks ------------------------------------------------------------------

# Primary: temperature enters only through the orthogonal anomaly, so the global
# coefficient is the total effect of warming.
local_orthogonal <- "RR + RR^2 + zRR^2 + zTM_orth + zTM_orth^2"
# Comparison: the notebook's local temperature level is retained, so the global
# coefficient is the partial effect holding local temperature fixed.
local_conditional <- "TM + TM^2 + RR + RR^2 + zRR^2 + zTM_orth + zTM_orth^2"
# The un-orthogonalised local anomaly, for the reparameterisation check.
local_raw <- "RR + RR^2 + zRR^2 + zTM + zTM^2"
global_block <- "zG + zG^2"

fe_full <- "year + GID_1[year] + GID_0[year]"
fe_no_year <- "GID_1 + GID_0[year]"
fe_region_only <- "GID_1"

fit_model <- function(terms, fe, data = estimation_panel, cluster = ~ GID_1 + year) {
  fixest::feols(
    stats::as.formula(paste("dlgrp_pc_usd ~", terms, "|", fe)),
    data = data,
    panel.id = c("GID_1", "year"),
    cluster = cluster
  )
}

tidy_fit <- function(fit, label, data = estimation_panel, note = NA_character_) {
  coefficients <- as.data.frame(fixest::coeftable(fit))
  names(coefficients) <- c("estimate", "std_error", "statistic", "p_value")
  intervals <- as.data.frame(
    stats::confint(fit, level = CONFIDENCE_LEVEL)
  )[rownames(coefficients), , drop = FALSE]
  # fixest::obs() gives the rows actually used, which is the only reliable way to
  # count clusters once NA removal and singleton dropping have run.
  used <- fixest::obs(fit)
  tibble::tibble(
    model = label,
    note = note,
    term = rownames(coefficients),
    estimate = coefficients$estimate,
    std_error = coefficients$std_error,
    statistic = coefficients$statistic,
    p_value = coefficients$p_value,
    conf_low = intervals[[1L]],
    conf_high = intervals[[2L]],
    observations = as.integer(stats::nobs(fit)),
    n_regions = n_distinct(data$GID_1[used]),
    n_countries = n_distinct(data$GID_0[used]),
    n_years = n_distinct(data$year[used])
  )
}

# Models -----------------------------------------------------------------------

models <- list()

models[["M1 total effect"]] <- fit_model(
  paste(local_orthogonal, "+", global_block), fe_no_year
)
models[["M1 partial effect"]] <- fit_model(
  paste(local_conditional, "+", global_block), fe_no_year
)
models[["M1 region FE only"]] <- fit_model(
  paste(local_orthogonal, "+", global_block), fe_region_only
)
models[["M1 cyclical global"]] <- fit_model(
  paste(local_orthogonal, "+ Gcyc + Gcyc^2"), fe_no_year
)
models[["M1 with cubic trend"]] <- fit_model(
  paste(local_orthogonal, "+", global_block, "+ poly(year, 3)"), fe_no_year
)
models[["M2 exposure"]] <- fit_model(
  paste(local_orthogonal, "+ lambda_c:Gcyc + lambda_c:Gtrend"), fe_full
)
models[["M2 exposure, region trends"]] <- fit_model(
  paste(local_orthogonal, "+ lambda_c:Gcyc"),
  "year + GID_1[year] + GID_0[year]"
)

model_comparison <- bind_rows(lapply(names(models), function(label) {
  tidy_fit(models[[label]], label)
}))

# Reparameterisation check: with a *pooled* loading and no squared terms, the
# (zG, zTM_orth) and (zG, zTM) parameterisations span the same column space, so
# they must have identical fit. With region-specific phi and squared terms they
# do not - which is the substantive point, not a bug.
pooled_phi <- with(
  estimation_panel %>% filter(is.finite(zTM), is.finite(zG)),
  stats::cov(zTM, zG) / stats::var(zG)
)
reparam_panel <- estimation_panel %>%
  mutate(zTM_orth_pooled = zTM - pooled_phi * zG)
m_reparam_orth <- fit_model(
  "RR + zRR^2 + zTM_orth_pooled + zG", fe_no_year, data = reparam_panel
)
m_reparam_raw <- fit_model(
  "RR + zRR^2 + zTM + zG", fe_no_year, data = reparam_panel
)
reparameterisation_check <- tibble::tibble(
  pooled_phi = pooled_phi,
  rss_orthogonal = stats::deviance(m_reparam_orth),
  rss_raw = stats::deviance(m_reparam_raw),
  rss_gap = abs(stats::deviance(m_reparam_orth) - stats::deviance(m_reparam_raw)),
  beta_local_orthogonal = unname(coef(m_reparam_orth)["zTM_orth_pooled"]),
  beta_local_raw = unname(coef(m_reparam_raw)["zTM"]),
  beta_global_orthogonal = unname(coef(m_reparam_orth)["zG"]),
  beta_global_raw = unname(coef(m_reparam_raw)["zG"]),
  implied_global_from_raw = unname(
    coef(m_reparam_raw)["zG"] + pooled_phi * coef(m_reparam_raw)["zTM"]
  )
)

# Fixed-effect ladder ----------------------------------------------------------

# How much of the global signal each trend control destroys, and what happens to
# the global coefficient as a result.
fe_ladder_specs <- tibble::tribble(
  ~rung, ~fixed_effects,
  "FE0  GID_1",                                   "GID_1",
  "FE1  GID_1 + GID_0[year]",                     "GID_1 + GID_0[year]",
  "FE2  GID_1 + GID_1[year]",                     "GID_1 + GID_1[year]",
  "FE3  GID_1 + GID_1[year] + GID_0[year]",       "GID_1 + GID_1[year] + GID_0[year]",
  "FE4  year + GID_1[year] + GID_0[year]",        "year + GID_1[year] + GID_0[year]"
)

retained_variation <- function(fe) {
  residualised <- fixest::feols(
    stats::as.formula(paste("zG ~ 1 |", fe)),
    data = estimation_panel %>% filter(is.finite(zG))
  )
  stats::sd(stats::resid(residualised)) /
    stats::sd(estimation_panel$zG[is.finite(estimation_panel$zG)])
}

fe_ladder <- bind_rows(lapply(seq_len(nrow(fe_ladder_specs)), function(i) {
  rung <- fe_ladder_specs$rung[[i]]
  fe <- fe_ladder_specs$fixed_effects[[i]]
  fit <- fit_model(paste(local_orthogonal, "+", global_block), fe)
  estimates <- as.data.frame(fixest::coeftable(fit))
  names(estimates) <- c("estimate", "std_error", "statistic", "p_value")
  has_global <- "zG" %in% rownames(estimates)
  tibble::tibble(
    rung = rung,
    fixed_effects = fe,
    zG_identified = has_global,
    retained_zG_variation = retained_variation(fe),
    estimate = if (has_global) estimates["zG", "estimate"] else NA_real_,
    std_error = if (has_global) estimates["zG", "std_error"] else NA_real_,
    p_value = if (has_global) estimates["zG", "p_value"] else NA_real_,
    observations = as.integer(stats::nobs(fit))
  )
}))

# Model 3: two-step year effects -------------------------------------------------

# The stage-1 fixed effects must NOT contain a year slope term such as
# `GID_0[year]`. With year dummies and unit-specific linear trends in the same
# model, the year effects are identified only up to an *affine* function of year:
# a linear trend can be moved between the year effects and the slopes without
# changing the fit. Regressing such year effects on a trending global shock would
# then report a normalisation, not an estimate. Restricting stage 1 to `GID_1 +
# year` keeps the year effects identified up to a constant, which stage 2's
# intercept absorbs.
fe_stage1 <- "GID_1 + year"
if (grepl("\\[[^]]*year", fe_stage1)) {
  stop(
    "Stage-1 fixed effects contain a year slope term; the extracted year ",
    "effects would be identified only up to an affine function of year."
  )
}

m_stage1 <- fixest::feols(
  stats::as.formula(paste("dlgrp_pc_usd ~", local_orthogonal, "|", fe_stage1)),
  data = estimation_panel,
  panel.id = c("GID_1", "year"),
  cluster = ~GID_1
)

year_effects <- fixest::fixef(m_stage1, sorted = FALSE)$year
year_effect_series <- tibble::tibble(
  year = as.integer(names(year_effects)),
  gamma = as.numeric(year_effects)
) %>%
  left_join(
    estimation_panel %>% count(year, name = "n_regions"),
    by = "year"
  ) %>%
  left_join(
    estimation_panel %>%
      distinct(year, zG, Gcyc, Gtrend, G),
    by = "year"
  ) %>%
  filter(is.finite(gamma), is.finite(zG)) %>%
  arrange(year)

newey_west_lag <- max(1L, floor(4 * (nrow(year_effect_series) / 100)^(2 / 9)))

fit_stage2 <- function(formula, label) {
  fit <- stats::lm(formula, data = year_effect_series, weights = n_regions)
  tested <- lmtest::coeftest(
    fit,
    vcov. = sandwich::NeweyWest(fit, lag = newey_west_lag, prewhite = FALSE)
  )
  tibble::tibble(
    model = label,
    term = rownames(tested),
    estimate = tested[, 1],
    std_error = tested[, 2],
    statistic = tested[, 3],
    p_value = tested[, 4],
    n_years = nrow(year_effect_series),
    newey_west_lag = newey_west_lag
  )
}

two_step_results <- bind_rows(
  fit_stage2(gamma ~ zG, "M3 year effects ~ zG"),
  fit_stage2(gamma ~ zG + I(zG^2), "M3 year effects ~ zG + zG^2"),
  fit_stage2(gamma ~ Gcyc, "M3 year effects ~ Gcyc"),
  fit_stage2(gamma ~ G, "M3 year effects ~ G (levels)")
)

# The fixef() normalisation is identified up to a constant, so re-centring the
# year effects must leave every slope unchanged.
year_effect_series_shifted <- year_effect_series %>% mutate(gamma = gamma + 10)
normalisation_gap <- max(abs(
  coef(stats::lm(gamma ~ zG, year_effect_series, weights = n_regions))[-1] -
    coef(stats::lm(gamma ~ zG, year_effect_series_shifted, weights = n_regions))[-1]
))
if (normalisation_gap > 1e-10) {
  stop("Stage-2 slopes depend on the year-effect normalisation.")
}

# Frisch-Waugh cross-check. Stage 1 removes the year effects that a one-step
# regression on the same controls would have used to identify zG, so the weighted
# stage-2 slope must reproduce the one-step coefficient. A large gap means the
# year effects were not cleanly extracted.
m_one_step_matched <- fit_model(
  paste(local_orthogonal, "+ zG"), fe_region_only, cluster = ~GID_1
)
two_step_fwl_check <- tibble::tibble(
  one_step_zG = unname(coef(m_one_step_matched)["zG"]),
  two_step_zG_weighted = unname(
    coef(stats::lm(gamma ~ zG, year_effect_series, weights = n_regions))["zG"]
  ),
  two_step_zG_unweighted = unname(
    coef(stats::lm(gamma ~ zG, year_effect_series))["zG"]
  )
) %>%
  mutate(
    absolute_gap = abs(one_step_zG - two_step_zG_weighted),
    relative_gap = absolute_gap / abs(one_step_zG)
  )
if (two_step_fwl_check$relative_gap > 0.15) {
  warning(
    "Two-step and one-step global coefficients differ by ",
    round(100 * two_step_fwl_check$relative_gap, 1),
    "%; the year effects may not be cleanly separable from the other controls."
  )
}

# Inference --------------------------------------------------------------------

set.seed(BOOTSTRAP_SEED)

# Year block bootstrap for the global coefficient. With one global series the
# effective sample for zG is the number of years, not the number of region-years;
# resampling years directly is the dependency-free analogue of a wild cluster
# bootstrap on the year dimension (fwildclusterboot is not installed here).
#
# Both bootstraps resample by row index rather than by filtering and rebinding,
# which is what makes a few hundred replications feasible on a 45k-row panel.
bootstrap_years <- function(model_terms, fe, term, replications) {
  years <- sort(unique(estimation_panel$year))
  rows_by_year <- split(seq_len(nrow(estimation_panel)), estimation_panel$year)
  rows_by_year <- rows_by_year[as.character(years)]
  formula <- stats::as.formula(paste("dlgrp_pc_usd ~", model_terms, "|", fe))

  estimates <- vapply(seq_len(replications), function(b) {
    drawn <- sample(seq_along(years), length(years), replace = TRUE)
    index <- unlist(rows_by_year[drawn], use.names = FALSE)
    resampled <- estimation_panel[index, , drop = FALSE]
    # Relabel so each drawn block occupies a distinct slot on the time axis.
    resampled$year <- rep(years, lengths(rows_by_year[drawn]))
    fit <- try(fixest::feols(formula, data = resampled), silent = TRUE)
    if (inherits(fit, "try-error") || !term %in% names(coef(fit))) {
      return(NA_real_)
    }
    unname(coef(fit)[term])
  }, numeric(1))
  estimates[is.finite(estimates)]
}

year_bootstrap <- bootstrap_years(
  paste(local_orthogonal, "+", global_block),
  fe_no_year,
  "zG",
  YEAR_BOOTSTRAP_REPLICATIONS
)

# Region block bootstrap that redoes *both* stages, so the uncertainty in the
# estimated pattern loading propagates into the local coefficients (Pagan 1984;
# Murphy & Topel 1985). This is honest for the local terms but not for zG: the
# single global path is held fixed across draws.
paired_for_bootstrap <- estimation_panel %>%
  filter(is.finite(zTM), is.finite(zG)) %>%
  select(GID_1, GID_0, year, dlgrp_pc_usd, TM, RR, zTM, zRR, zG, Gcyc, Gtrend)

rows_by_region <- split(
  seq_len(nrow(paired_for_bootstrap)), paired_for_bootstrap$GID_1
)
n_regions_boot <- length(rows_by_region)
bootstrap_formula <- stats::as.formula(
  paste("dlgrp_pc_usd ~", local_orthogonal, "+", global_block, "| GID_1")
)

bootstrap_two_stage <- function(replications) {
  bind_rows(lapply(seq_len(replications), function(b) {
    drawn <- sample(n_regions_boot, n_regions_boot, replace = TRUE)
    blocks <- rows_by_region[drawn]
    index <- unlist(blocks, use.names = FALSE)
    resampled <- paired_for_bootstrap[index, , drop = FALSE]
    # A region drawn twice must enter as two distinct units, otherwise the
    # cluster structure of the resample is wrong.
    resampled$GID_1 <- rep(seq_along(blocks), lengths(blocks))

    # Stage 1: re-estimate the pattern loading on the resampled regions, so the
    # first-stage uncertainty propagates into the second-stage coefficients.
    resampled <- resampled %>%
      group_by(GID_1) %>%
      mutate(
        phi_b = stats::cov(zTM, zG) / stats::var(zG),
        zTM_orth = zTM - (mean(zTM) - phi_b * mean(zG)) - phi_b * zG
      ) %>%
      ungroup()

    fit <- try(fixest::feols(bootstrap_formula, data = resampled), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    tibble::tibble(
      replication = b,
      term = names(coef(fit)),
      estimate = unname(coef(fit))
    )
  }))
}

region_bootstrap <- bootstrap_two_stage(BOOTSTRAP_REPLICATIONS)

bootstrap_inference <- bind_rows(
  tibble::tibble(
    bootstrap = "year block (global coefficient)",
    term = "zG",
    replications = length(year_bootstrap),
    estimate = mean(year_bootstrap),
    std_error = stats::sd(year_bootstrap),
    conf_low = unname(stats::quantile(year_bootstrap, 0.025)),
    conf_high = unname(stats::quantile(year_bootstrap, 0.975))
  ),
  region_bootstrap %>%
    group_by(term) %>%
    summarise(
      bootstrap = "region block, both stages",
      replications = n(),
      # Every statistic is taken from the draws before `estimate` is reassigned:
      # naming an output column after an input column shadows it for the rest of
      # the summarise, which would silently collapse the interval onto the mean.
      boot_mean = mean(estimate),
      boot_sd = stats::sd(estimate),
      boot_low = unname(stats::quantile(estimate, 0.025)),
      boot_high = unname(stats::quantile(estimate, 0.975)),
      .groups = "drop"
    ) %>%
    transmute(
      bootstrap, term, replications,
      estimate = boot_mean,
      std_error = boot_sd,
      conf_low = boot_low,
      conf_high = boot_high
    )
)

# Local projections --------------------------------------------------------------

lp_panel <- estimation_panel %>%
  group_by(GID_1) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(lgrp_lag1 = lgrp_pc_usd[match(year - 1L, year)]) %>%
  ungroup()

local_projections <- bind_rows(lapply(LP_HORIZONS, function(h) {
  horizon_panel <- lp_panel %>%
    group_by(GID_1) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(cumulative = lgrp_pc_usd[match(year + h, year)] - lgrp_lag1) %>%
    ungroup() %>%
    filter(is.finite(cumulative))

  fit <- fixest::feols(
    stats::as.formula(
      paste("cumulative ~", local_orthogonal, "+", global_block, "|", fe_no_year)
    ),
    data = horizon_panel,
    panel.id = c("GID_1", "year"),
    cluster = ~ GID_1 + year
  )
  estimates <- as.data.frame(fixest::coeftable(fit))
  names(estimates) <- c("estimate", "std_error", "statistic", "p_value")
  intervals <- as.data.frame(stats::confint(fit, level = CONFIDENCE_LEVEL))

  tibble::tibble(
    horizon = h,
    term = rownames(estimates),
    estimate = estimates$estimate,
    std_error = estimates$std_error,
    p_value = estimates$p_value,
    conf_low = intervals[rownames(estimates), 1],
    conf_high = intervals[rownames(estimates), 2],
    observations = as.integer(stats::nobs(fit))
  )
}))

# Outputs -----------------------------------------------------------------------

sample_sizes <- tibble::tibble(
  panel_rows = nrow(panel),
  estimation_rows = nrow(estimation_panel),
  regions = n_distinct(estimation_panel$GID_1),
  countries = n_distinct(estimation_panel$GID_0),
  years = n_distinct(estimation_panel$year),
  first_year = min(estimation_panel$year),
  last_year = max(estimation_panel$year),
  rows_with_global_shock = sum(is.finite(estimation_panel$zG)),
  rows_with_local_shock = sum(is.finite(estimation_panel$zTM_orth)),
  regions_with_loading = n_distinct(
    estimation_panel$GID_1[is.finite(estimation_panel$phi)]
  )
)

write.csv(model_comparison, file.path(output_dir, "model_comparison.csv"), row.names = FALSE)
write.csv(fe_ladder, file.path(output_dir, "fe_ladder.csv"), row.names = FALSE)
write.csv(two_step_results, file.path(output_dir, "two_step_year_effects.csv"), row.names = FALSE)
write.csv(two_step_fwl_check, file.path(output_dir, "two_step_fwl_check.csv"), row.names = FALSE)
write.csv(year_effect_series, file.path(output_dir, "year_effect_series.csv"), row.names = FALSE)
write.csv(bootstrap_inference, file.path(output_dir, "bootstrap_inference.csv"), row.names = FALSE)
write.csv(local_projections, file.path(output_dir, "local_projection_irf.csv"), row.names = FALSE)
write.csv(reparameterisation_check, file.path(output_dir, "reparameterisation_check.csv"), row.names = FALSE)
write.csv(sample_sizes, file.path(output_dir, "sample_sizes.csv"), row.names = FALSE)

# Plots -------------------------------------------------------------------------

fe_ladder_plot <- fe_ladder %>%
  mutate(rung = factor(rung, levels = rev(fe_ladder_specs$rung))) %>%
  ggplot2::ggplot(ggplot2::aes(x = retained_zG_variation, y = rung)) +
  ggplot2::geom_col(fill = "#0072B2", width = 0.6) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "How much of the global shock each fixed-effect set destroys",
    subtitle = paste(
      "Share of the standard deviation of zG surviving the fixed effects",
      "(Frisch-Waugh residualisation)."
    ),
    x = "Retained variation in zG",
    y = NULL
  ) +
  ggplot2::theme_classic()

lp_plot <- local_projections %>%
  filter(term %in% c("zG", "zTM_orth")) %>%
  mutate(term = if_else(
    term == "zG", "Global shock (zG)", "Orthogonal local shock (zTM_orth)"
  )) %>%
  ggplot2::ggplot(ggplot2::aes(x = horizon, y = estimate)) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = conf_low, ymax = conf_high),
    width = 0.06,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(linewidth = 0.6, colour = "#0072B2") +
  ggplot2::geom_point(size = 2, colour = "#0072B2") +
  ggplot2::facet_wrap(~term, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = LP_HORIZONS) +
  ggplot2::labs(
    title = "Cumulative growth response to global and orthogonal local shocks",
    subtitle = paste0(
      "Local projections, outcome log GDP p.c. at t+h minus t-1; ",
      100 * CONFIDENCE_LEVEL, "% CI, two-way clustered by region and year."
    ),
    x = "Horizon (years)",
    y = "Cumulative log GDP per capita response"
  ) +
  ggplot2::theme_classic()

year_effect_plot <- ggplot2::ggplot(
  year_effect_series,
  ggplot2::aes(x = zG, y = gamma)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(size = n_regions), alpha = 0.6, colour = "#0072B2") +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "#E69F00",
    linewidth = 0.8
  ) +
  ggplot2::scale_size_continuous(range = c(0.8, 4)) +
  ggplot2::labs(
    title = "Common year effects against the global temperature shock",
    subtitle = paste0(
      "Stage-2 of the two-step estimator: ", nrow(year_effect_series),
      " years, weighted by regions observed, Newey-West lag ", newey_west_lag, "."
    ),
    x = "Global shock zG",
    y = "Estimated year effect",
    size = "Regions"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(file.path(output_dir, "fe_ladder.png"), fe_ladder_plot,
                width = 9, height = 4, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(output_dir, "local_projection_irf.png"), lp_plot,
                width = 10, height = 4.5, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(output_dir, "two_step_year_effects.png"), year_effect_plot,
                width = 8, height = 4.5, dpi = 300, bg = "white")

# Console report ----------------------------------------------------------------

cat("\nSample:\n")
print(as.data.frame(sample_sizes))

cat("\nFixed-effect ladder (coefficient on zG):\n")
print(as.data.frame(fe_ladder), digits = 4)

cat("\nModel comparison (global and local shock terms):\n")
print(
  as.data.frame(
    model_comparison %>%
      filter(grepl("^(zG|zTM_orth|zTM|lambda_c|Gcyc)", term)) %>%
      select(model, term, estimate, std_error, p_value, observations)
  ),
  digits = 4
)

cat("\nTwo-step year effects (Newey-West):\n")
print(as.data.frame(two_step_results %>% filter(term != "(Intercept)")), digits = 4)

cat("\nTwo-step vs one-step (Frisch-Waugh cross-check on zG):\n")
print(as.data.frame(two_step_fwl_check), digits = 4)

cat("\nReparameterisation check (pooled loading, no squares):\n")
print(as.data.frame(reparameterisation_check), digits = 6)

cat("\nBootstrap inference:\n")
print(as.data.frame(bootstrap_inference), digits = 4)

cat("\nWritten to: ", normalizePath(output_dir), "\n", sep = "")
