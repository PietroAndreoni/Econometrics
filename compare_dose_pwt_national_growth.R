# Compare national GDP growth implied by raw DOSE subnational observations with
# national GDP growth in PWT 11.0.
#
# DOSE regional GDP is reconstructed as:
#   regional GDP = regional GDP per capita (2015 USD) * regional population
# and summed by country-year. PWT `rgdpna` is the national-accounts real-GDP
# series (millions of 2021 USD). The level units differ, but annual growth rates
# are comparable because each series is measured consistently through time.
#
# The final section uses the region-year discrepancy ranking produced here as an
# outlier definition, and traces how the subnational climate response
# coefficients move as those outliers are dropped from the DOSE estimation
# sample.

required_packages <- c("dplyr", "readxl", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(library(dplyr))

# Configuration ----------------------------------------------------------------

FLAG_PERCENTILE <- 0.99
TOP_N_TO_WRITE <- 50L
TOP_N_ATTRIBUTION_DETAIL <- 20L
ATTRIBUTION_TOLERANCE_PP <- 1e-8

if (file.exists("Econometrics.Rproj")) {
  project_root <- "."
} else if (file.exists(file.path("..", "Econometrics.Rproj"))) {
  project_root <- ".."
} else {
  stop("Run this script from the project root or the econometrics directory.")
}

dose_path <- file.path(project_root, "data", "DOSE_V2.11.csv")
pwt_path <- file.path(project_root, "data", "pwt110.xlsx")
output_dir <- file.path(
  project_root,
  "results",
  "dose_pwt_national_growth"
)

missing_inputs <- c(dose_path, pwt_path)[!file.exists(c(dose_path, pwt_path))]
if (length(missing_inputs)) {
  stop("Raw input files not found: ", paste(missing_inputs, collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helpers ----------------------------------------------------------------------

lag_if_consecutive <- function(x, year) {
  previous <- dplyr::lag(x)
  if_else(year == dplyr::lag(year) + 1L, previous, NA_real_)
}

annual_log_growth <- function(x, year) {
  previous <- lag_if_consecutive(x, year)
  if_else(
    is.finite(x) & x > 0 & is.finite(previous) & previous > 0,
    log(x) - log(previous),
    NA_real_
  )
}

# Raw DOSE aggregation ---------------------------------------------------------

dose_raw <- read.csv(
  dose_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    country = as.character(country),
    iso3 = toupper(trimws(coalesce(as.character(GID_0), ""))),
    region_name = trimws(coalesce(as.character(region), "")),
    region_id_file = trimws(coalesce(as.character(GID_1), "")),
    year = as.integer(year),
    regional_gdp_per_capita = as.numeric(grp_pc_usd_2015),
    regional_population = as.numeric(pop)
  ) %>%
  mutate(
    # Some DOSE rows have a region name but no GID_1. Dropping those rows would
    # silently bias the national aggregate, so give them a stable fallback key.
    # GID_0, rather than the GID_1 prefix, is the DOSE country identifier.
    region_id = if_else(
      nzchar(region_id_file),
      region_id_file,
      paste0(iso3, "::UNMAPPED::", region_name)
    ),
    region_id_prefix = substr(region_id_file, 1L, 3L),
    used_fallback_region_id = !nzchar(region_id_file),
    region_id_country_mismatch = nzchar(region_id_file) &
      region_id_prefix != iso3,
    valid_component = is.finite(regional_gdp_per_capita) &
      regional_gdp_per_capita > 0 &
      is.finite(regional_population) &
      regional_population > 0,
    regional_gdp = if_else(
      valid_component,
      regional_gdp_per_capita * regional_population,
      NA_real_
    )
  )

invalid_dose_identifiers <- dose_raw %>%
  filter(
    is.na(year) |
      !grepl("^[A-Z]{3}$", iso3) |
      (!nzchar(region_id_file) & !nzchar(region_name))
  )
if (nrow(invalid_dose_identifiers)) {
  stop(
    "DOSE contains missing years, invalid country codes, or rows with neither ",
    "a region ID nor a region name. Resolve these fields before aggregating."
  )
}

dose_identifier_issues <- dose_raw %>%
  filter(used_fallback_region_id | region_id_country_mismatch) %>%
  distinct(
    country,
    iso3,
    region_name,
    region_id_file,
    region_id,
    used_fallback_region_id,
    region_id_country_mismatch
  ) %>%
  arrange(iso3, region_id)

duplicate_dose_keys <- dose_raw %>%
  count(iso3, region_id, year, name = "rows") %>%
  filter(rows > 1L)
if (nrow(duplicate_dose_keys)) {
  write.csv(
    duplicate_dose_keys,
    file.path(output_dir, "duplicate_dose_region_year_keys.csv"),
    row.names = FALSE
  )
  stop(
    "DOSE has duplicate region-year keys. See ",
    file.path(output_dir, "duplicate_dose_region_year_keys.csv")
  )
}

dose_national <- dose_raw %>%
  group_by(iso3, year) %>%
  summarise(
    country_dose = first(country),
    n_regions_reported = n_distinct(region_id),
    n_regions_complete = n_distinct(region_id[valid_component]),
    n_incomplete_regions = n_regions_reported - n_regions_complete,
    region_signature = paste(sort(unique(region_id[valid_component])), collapse = "|"),
    population_dose = sum(regional_population[valid_component]),
    gdp_dose = sum(regional_gdp[valid_component]),
    .groups = "drop"
  ) %>%
  filter(n_regions_complete > 0L, population_dose > 0, gdp_dose > 0) %>%
  mutate(
    gdppc_dose = gdp_dose / population_dose,
    complete_regional_coverage = n_incomplete_regions == 0L
  ) %>%
  group_by(iso3) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    consecutive_dose_year = year == dplyr::lag(year) + 1L,
    gdp_dose_previous = lag_if_consecutive(gdp_dose, year),
    stable_region_coverage = consecutive_dose_year &
      region_signature == dplyr::lag(region_signature) &
      complete_regional_coverage &
      dplyr::lag(complete_regional_coverage, default = FALSE),
    gdp_growth_log_dose = annual_log_growth(gdp_dose, year),
    gdppc_growth_log_dose = annual_log_growth(gdppc_dose, year),
    population_growth_log_dose = annual_log_growth(population_dose, year)
  ) %>%
  ungroup()

# Raw PWT national series ------------------------------------------------------

pwt_national <- readxl::read_excel(
  pwt_path,
  sheet = "Data"
) %>%
  transmute(
    iso3 = toupper(trimws(coalesce(as.character(countrycode), ""))),
    country_pwt = as.character(country),
    year = as.integer(year),
    # rgdpna is in millions of constant 2021 USD; pop is in millions.
    gdp_pwt = as.numeric(rgdpna),
    population_pwt = as.numeric(pop)
  ) %>%
  filter(
    nzchar(iso3),
    !is.na(year),
    is.finite(gdp_pwt),
    gdp_pwt > 0,
    is.finite(population_pwt),
    population_pwt > 0
  ) %>%
  mutate(gdppc_pwt = gdp_pwt / population_pwt) %>%
  group_by(iso3) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    consecutive_pwt_year = year == dplyr::lag(year) + 1L,
    gdp_pwt_previous = lag_if_consecutive(gdp_pwt, year),
    gdp_growth_log_pwt = annual_log_growth(gdp_pwt, year),
    gdppc_growth_log_pwt = annual_log_growth(gdppc_pwt, year),
    population_growth_log_pwt = annual_log_growth(population_pwt, year)
  ) %>%
  ungroup()

duplicate_pwt_keys <- pwt_national %>%
  count(iso3, year, name = "rows") %>%
  filter(rows > 1L)
if (nrow(duplicate_pwt_keys)) {
  write.csv(
    duplicate_pwt_keys,
    file.path(output_dir, "duplicate_pwt_country_year_keys.csv"),
    row.names = FALSE
  )
  stop(
    "PWT has duplicate country-year keys. See ",
    file.path(output_dir, "duplicate_pwt_country_year_keys.csv")
  )
}

# Growth comparison and discrepancy flags ------------------------------------

growth_comparison <- dose_national %>%
  select(-region_signature) %>%
  inner_join(pwt_national, by = c("iso3", "year")) %>%
  mutate(
    eligible_for_comparison = consecutive_dose_year &
      consecutive_pwt_year &
      stable_region_coverage &
      is.finite(gdp_growth_log_dose) &
      is.finite(gdp_growth_log_pwt),
    gdp_growth_pct_dose = 100 * expm1(gdp_growth_log_dose),
    gdp_growth_pct_pwt = 100 * expm1(gdp_growth_log_pwt),
    gdp_growth_discrepancy_pp = gdp_growth_pct_dose - gdp_growth_pct_pwt,
    abs_gdp_growth_discrepancy_pp = abs(gdp_growth_discrepancy_pp),
    gdppc_growth_pct_dose = 100 * expm1(gdppc_growth_log_dose),
    gdppc_growth_pct_pwt = 100 * expm1(gdppc_growth_log_pwt),
    gdppc_growth_discrepancy_pp = gdppc_growth_pct_dose - gdppc_growth_pct_pwt,
    population_growth_pct_dose = 100 * expm1(population_growth_log_dose),
    population_growth_pct_pwt = 100 * expm1(population_growth_log_pwt),
    population_growth_discrepancy_pp =
      population_growth_pct_dose - population_growth_pct_pwt
  )

eligible_discrepancies <- growth_comparison$abs_gdp_growth_discrepancy_pp[
  growth_comparison$eligible_for_comparison
]
if (!length(eligible_discrepancies)) {
  stop("No eligible overlapping country-year growth observations were found.")
}

flag_threshold_pp <- unname(stats::quantile(
  eligible_discrepancies,
  probs = FLAG_PERCENTILE,
  na.rm = TRUE
))

growth_comparison <- growth_comparison %>%
  arrange(desc(eligible_for_comparison), desc(abs_gdp_growth_discrepancy_pp)) %>%
  mutate(
    discrepancy_rank = if_else(
      eligible_for_comparison,
      as.integer(rank(
        if_else(
          eligible_for_comparison,
          -abs_gdp_growth_discrepancy_pp,
          NA_real_
        ),
        ties.method = "min",
        na.last = "keep"
      )),
      NA_integer_
    ),
    flag_high_discrepancy = eligible_for_comparison &
      abs_gdp_growth_discrepancy_pp >= flag_threshold_pp,
    flag_top_n = eligible_for_comparison &
      discrepancy_rank <= TOP_N_TO_WRITE
  )

flagged_country_years <- growth_comparison %>%
  filter(flag_high_discrepancy) %>%
  arrange(discrepancy_rank)

top_country_years <- growth_comparison %>%
  filter(flag_top_n) %>%
  arrange(discrepancy_rank)

# Subregional attribution -----------------------------------------------------

# This is an exact shift-share decomposition. For country c, year t, let
# Y_it be DOSE GDP in region i, Y_t = sum_i Y_it, and q_PWT be PWT's national
# GDP growth factor. Region i's contribution in percentage points is
#
#   100 * (Y_it - q_PWT * Y_i,t-1) / Y_t-1
# = lagged DOSE GDP share_i * (regional growth % - PWT growth %).
#
# Summing over regions gives exactly
# 100 * [(Y_t / Y_t-1) - q_PWT], the reported DOSE-PWT growth discrepancy.
# Positive contributions in the direction of the national discrepancy widen
# it; contributions with the opposite sign offset it. A leave-one-region-out
# calculation is included as a non-additive sensitivity check.

dose_regional_growth <- dose_raw %>%
  filter(valid_component) %>%
  group_by(iso3, region_id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    regional_gdp_previous = lag_if_consecutive(regional_gdp, year),
    regional_growth_pct = if_else(
      is.finite(regional_gdp_previous) & regional_gdp_previous > 0,
      100 * (regional_gdp / regional_gdp_previous - 1),
      NA_real_
    )
  ) %>%
  ungroup()

attribution_targets <- growth_comparison %>%
  filter(eligible_for_comparison) %>%
  mutate(previous_year = year - 1L) %>%
  select(
    iso3,
    year,
    previous_year,
    country_dose,
    country_pwt,
    discrepancy_rank,
    n_regions_complete,
    gdp_dose,
    gdp_dose_previous,
    gdp_growth_pct_dose,
    gdp_growth_pct_pwt,
    gdp_growth_discrepancy_pp,
    abs_gdp_growth_discrepancy_pp
  )

subregion_attribution <- attribution_targets %>%
  inner_join(
    dose_regional_growth %>%
      select(
        iso3,
        year,
        region_id,
        region_name,
        used_fallback_region_id,
        regional_gdp,
        regional_gdp_previous,
        regional_growth_pct
      ),
    by = c("iso3", "year")
  ) %>%
  mutate(
    lagged_gdp_share = regional_gdp_previous / gdp_dose_previous,
    additive_contribution_pp = lagged_gdp_share *
      (regional_growth_pct - gdp_growth_pct_pwt),
    abs_additive_contribution_pp = abs(additive_contribution_pp),
    discrepancy_aligned_contribution_pp = sign(gdp_growth_discrepancy_pp) *
      additive_contribution_pp,
    contribution_role = case_when(
      additive_contribution_pp * gdp_growth_discrepancy_pp > 0 ~ "widens",
      additive_contribution_pp * gdp_growth_discrepancy_pp < 0 ~ "offsets",
      TRUE ~ "neutral"
    ),
    share_of_national_discrepancy = if_else(
      abs(gdp_growth_discrepancy_pp) > ATTRIBUTION_TOLERANCE_PP,
      additive_contribution_pp / gdp_growth_discrepancy_pp,
      NA_real_
    ),
    gdp_dose_without_region = gdp_dose - regional_gdp,
    gdp_dose_previous_without_region =
      gdp_dose_previous - regional_gdp_previous,
    leave_one_out_growth_pct_dose = if_else(
      gdp_dose_without_region > 0 & gdp_dose_previous_without_region > 0,
      100 * (
        gdp_dose_without_region / gdp_dose_previous_without_region - 1
      ),
      NA_real_
    ),
    leave_one_out_discrepancy_pp =
      leave_one_out_growth_pct_dose - gdp_growth_pct_pwt,
    leave_one_out_abs_discrepancy_reduction_pp =
      abs_gdp_growth_discrepancy_pp - abs(leave_one_out_discrepancy_pp),
    leave_one_out_improves_fit =
      leave_one_out_abs_discrepancy_reduction_pp > 0
  ) %>%
  group_by(iso3, year) %>%
  arrange(
    desc(discrepancy_aligned_contribution_pp),
    desc(abs_additive_contribution_pp),
    region_id,
    .by_group = TRUE
  ) %>%
  mutate(
    discrepancy_driver_rank = row_number(),
    absolute_contribution_rank = row_number(
      desc(abs_additive_contribution_pp)
    ),
    leave_one_out_reduction_rank = row_number(
      desc(leave_one_out_abs_discrepancy_reduction_pp)
    )
  ) %>%
  ungroup()

attribution_sums <- subregion_attribution %>%
  group_by(iso3, year) %>%
  summarise(
    attributed_regions = n(),
    sum_additive_contribution_pp = sum(additive_contribution_pp),
    lagged_share_sum = sum(lagged_gdp_share),
    .groups = "drop"
  )

attribution_validation <- attribution_targets %>%
  transmute(
    iso3,
    year,
    discrepancy_rank,
    expected_regions = n_regions_complete,
    national_discrepancy_pp = gdp_growth_discrepancy_pp
  ) %>%
  left_join(attribution_sums, by = c("iso3", "year")) %>%
  mutate(
    attributed_regions = coalesce(attributed_regions, 0L),
    additive_closure_error_pp =
      sum_additive_contribution_pp - national_discrepancy_pp
  )

invalid_attributions <- attribution_validation %>%
  filter(
    attributed_regions != expected_regions |
      !is.finite(additive_closure_error_pp) |
      abs(additive_closure_error_pp) > ATTRIBUTION_TOLERANCE_PP |
      !is.finite(lagged_share_sum) |
      abs(lagged_share_sum - 1) > ATTRIBUTION_TOLERANCE_PP
  )
if (nrow(invalid_attributions)) {
  write.csv(
    invalid_attributions,
    file.path(output_dir, "invalid_subregion_attributions.csv"),
    row.names = FALSE
  )
  stop(
    "Subregional attribution did not reconcile with the national comparison. ",
    "See ", file.path(output_dir, "invalid_subregion_attributions.csv")
  )
}

primary_additive_driver <- subregion_attribution %>%
  filter(discrepancy_driver_rank == 1L) %>%
  transmute(
    iso3,
    year,
    primary_driver_region_id = region_id,
    primary_driver_region_name = region_name,
    primary_driver_contribution_pp = additive_contribution_pp,
    primary_driver_share_of_discrepancy = share_of_national_discrepancy,
    primary_driver_regional_growth_pct = regional_growth_pct
  )

largest_absolute_contributor <- subregion_attribution %>%
  filter(absolute_contribution_rank == 1L) %>%
  transmute(
    iso3,
    year,
    largest_absolute_region_id = region_id,
    largest_absolute_region_name = region_name,
    largest_absolute_contribution_pp = additive_contribution_pp,
    largest_absolute_contribution_role = contribution_role
  )

largest_leave_one_out_reduction <- subregion_attribution %>%
  filter(leave_one_out_reduction_rank == 1L) %>%
  transmute(
    iso3,
    year,
    leave_one_out_driver_region_id = region_id,
    leave_one_out_driver_region_name = region_name,
    leave_one_out_abs_discrepancy_reduction_pp,
    leave_one_out_discrepancy_pp,
    leave_one_out_improves_fit
  )

country_year_subregion_driver <- attribution_targets %>%
  left_join(primary_additive_driver, by = c("iso3", "year")) %>%
  left_join(largest_absolute_contributor, by = c("iso3", "year")) %>%
  left_join(largest_leave_one_out_reduction, by = c("iso3", "year")) %>%
  left_join(
    attribution_validation %>%
      select(
        iso3,
        year,
        sum_additive_contribution_pp,
        additive_closure_error_pp
      ),
    by = c("iso3", "year")
  ) %>%
  arrange(discrepancy_rank)

top_country_year_subregion_attribution <- subregion_attribution %>%
  filter(discrepancy_rank <= TOP_N_ATTRIBUTION_DETAIL) %>%
  arrange(discrepancy_rank, discrepancy_driver_rank)

country_summary <- growth_comparison %>%
  filter(eligible_for_comparison) %>%
  group_by(iso3) %>%
  summarise(
    country_dose = first(country_dose),
    country_pwt = first(country_pwt),
    observations = n(),
    mean_signed_discrepancy_pp = mean(gdp_growth_discrepancy_pp),
    mean_absolute_discrepancy_pp = mean(abs_gdp_growth_discrepancy_pp),
    rmse_discrepancy_pp = sqrt(mean(gdp_growth_discrepancy_pp^2)),
    maximum_absolute_discrepancy_pp = max(abs_gdp_growth_discrepancy_pp),
    flagged_country_years = sum(flag_high_discrepancy),
    .groups = "drop"
  ) %>%
  arrange(desc(rmse_discrepancy_pp))

source_coverage <- full_join(
  dose_national %>% distinct(iso3) %>% mutate(in_dose = TRUE),
  pwt_national %>% distinct(iso3) %>% mutate(in_pwt = TRUE),
  by = "iso3"
) %>%
  mutate(
    in_dose = coalesce(in_dose, FALSE),
    in_pwt = coalesce(in_pwt, FALSE)
  ) %>%
  arrange(iso3)

# Outputs ----------------------------------------------------------------------

write.csv(
  growth_comparison,
  file.path(output_dir, "country_year_growth_comparison.csv"),
  row.names = FALSE
)
write.csv(
  flagged_country_years,
  file.path(output_dir, "flagged_high_discrepancy_country_years.csv"),
  row.names = FALSE
)
write.csv(
  top_country_years,
  file.path(output_dir, paste0("top_", TOP_N_TO_WRITE, "_country_years.csv")),
  row.names = FALSE
)
write.csv(
  country_summary,
  file.path(output_dir, "country_discrepancy_summary.csv"),
  row.names = FALSE
)
write.csv(
  source_coverage,
  file.path(output_dir, "source_country_coverage.csv"),
  row.names = FALSE
)
write.csv(
  dose_identifier_issues,
  file.path(output_dir, "dose_identifier_issues.csv"),
  row.names = FALSE
)
write.csv(
  subregion_attribution,
  file.path(output_dir, "subregion_discrepancy_attribution.csv"),
  row.names = FALSE
)
write.csv(
  country_year_subregion_driver,
  file.path(output_dir, "country_year_subregion_driver.csv"),
  row.names = FALSE
)
write.csv(
  top_country_year_subregion_attribution,
  file.path(
    output_dir,
    paste0(
      "top_", TOP_N_ATTRIBUTION_DETAIL,
      "_country_year_subregion_attribution.csv"
    )
  ),
  row.names = FALSE
)
write.csv(
  attribution_validation,
  file.path(output_dir, "subregion_attribution_validation.csv"),
  row.names = FALSE
)

cat(
  "Compared ", sum(growth_comparison$eligible_for_comparison),
  " eligible country-year growth pairs across ",
  n_distinct(growth_comparison$iso3[growth_comparison$eligible_for_comparison]),
  " countries.\n",
  sep = ""
)
cat(
  "High-discrepancy threshold (", 100 * FLAG_PERCENTILE,
  "th percentile): ", round(flag_threshold_pp, 3),
  " percentage points.\n",
  sep = ""
)
cat("\nLargest discrepancies:\n")
print(
  top_country_years %>%
    select(
      discrepancy_rank,
      iso3,
      year,
      gdp_growth_pct_dose,
      gdp_growth_pct_pwt,
      gdp_growth_discrepancy_pp,
      n_regions_complete
    ) %>%
    head(20L),
  n = 20L
)
cat("\nPrimary subregional drivers for the largest discrepancies:\n")
print(
  country_year_subregion_driver %>%
    select(
      discrepancy_rank,
      iso3,
      year,
      gdp_growth_discrepancy_pp,
      primary_driver_region_id,
      primary_driver_region_name,
      primary_driver_contribution_pp,
      primary_driver_share_of_discrepancy,
      leave_one_out_driver_region_name,
      leave_one_out_abs_discrepancy_reduction_pp
    ) %>%
    head(TOP_N_ATTRIBUTION_DETAIL),
  n = TOP_N_ATTRIBUTION_DETAIL
)
cat(
  "\nSubregional contributions reconciled for ",
  nrow(attribution_validation),
  " country-year pairs; maximum absolute closure error: ",
  format(
    max(abs(attribution_validation$additive_closure_error_pp)),
    scientific = TRUE,
    digits = 3
  ),
  " percentage points.\n",
  sep = ""
)
cat("\nResults written to: ", normalizePath(output_dir), "\n", sep = "")

# Outlier-removal sensitivity of the subnational climate model ------------------

# The subregional attribution above ranks every eligible region-year by how much
# it contributes to its country's DOSE-PWT growth discrepancy. Those region-years
# are the outliers tested here: the worst-contributing ones are dropped from the
# DOSE subnational estimation sample in increasing numbers, the climate response
# signed-bin model is re-estimated on each reduced sample, and the temperature
# and precipitation response curves are overlaid. Estimation reuses `build_dat()` and
# `fit_climate_model()` from `test_functions.Rmd`, so sample construction, fixed
# effects, and clustered standard errors match the notebook exactly.

OUTLIER_REMOVAL_COUNTS <- seq(1000L, 5000L, by = 1000L)

# Region-year outlier ranking. "abs_additive_contribution_pp" is the magnitude of
# the region's exact shift-share contribution to its country's growth gap;
# "leave_one_out_abs_discrepancy_reduction_pp" instead ranks by how much dropping
# the region shrinks that gap.
OUTLIER_RANK_METRIC <- "abs_additive_contribution_pp"

# Joint signed-deviation-bin response from the notebook. The central -0.5 to 0.5
# standard-deviation interval is the omitted category for both weather variables.
SENSITIVITY_MODEL_TERMS <- paste(
  "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
  "+ i(TM_bin_signed, ref = 0)",
  "+ i(RR_bin_signed, ref = 0)"
)
SENSITIVITY_CONFIDENCE_LEVEL <- 0.95
SENSITIVITY_ECON_DATA <- "DOSE_V2_11"

test_functions_rmd <- file.path(project_root, "test_functions.Rmd")
data_dir <- file.path(project_root, "data")

sensitivity_packages <- c(
  "arrow", "fixest", "ggplot2", "purrr", "tibble", "tidyselect"
)
missing_sensitivity_packages <- sensitivity_packages[
  !vapply(sensitivity_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_sensitivity_packages)) {
  stop(
    "Install missing packages before running the outlier-removal section: ",
    paste(missing_sensitivity_packages, collapse = ", ")
  )
}
if (!file.exists(test_functions_rmd)) {
  stop("test_functions.Rmd not found at: ", test_functions_rmd)
}

# Reuse the notebook's data builder ---------------------------------------------

# The notebook chunks are evaluated in a private environment rather than sourced
# wholesale, so only the specification constants and the data builders are
# imported. `load_notebook_env()` suppresses the notebook's own `data` objects;
# the estimation sample is built below with explicit file paths instead.

source(file.path(project_root, "rmd_chunks.R"))

notebook <- load_notebook_env(
  test_functions_rmd,
  econ_data = SENSITIVITY_ECON_DATA
)

# DOSE subnational estimation sample -------------------------------------------

estimation_sample <- notebook$build_dat(
  econ_data = SENSITIVITY_ECON_DATA,
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
  filter(!is.na(dlgrp_pc_usd)) %>%
  mutate(year = as.integer(year))

# Restrict once to the exact complete-case/fixed-effect sample used by the
# untrimmed signed-bin model. Outlier counts below therefore count observations
# that would otherwise enter that model.
baseline_sensitivity_fit <- notebook$fit_climate_model(
  SENSITIVITY_MODEL_TERMS,
  model_data = estimation_sample
)
estimation_sample <- estimation_sample %>%
  ungroup() %>%
  slice(fixest::obs(baseline_sensitivity_fit))

# Region-year outlier ranking ---------------------------------------------------

if (!OUTLIER_RANK_METRIC %in% names(subregion_attribution)) {
  stop(
    "OUTLIER_RANK_METRIC '", OUTLIER_RANK_METRIC,
    "' is not a column of subregion_attribution."
  )
}

ranked_region_year_outliers <- subregion_attribution %>%
  filter(is.finite(.data[[OUTLIER_RANK_METRIC]])) %>%
  arrange(desc(.data[[OUTLIER_RANK_METRIC]]), iso3, region_id, year) %>%
  transmute(
    outlier_rank = row_number(),
    GID_0 = iso3,
    GID_1 = region_id,
    year = as.integer(year),
    region_name,
    outlier_metric = .data[[OUTLIER_RANK_METRIC]],
    additive_contribution_pp,
    regional_growth_pct,
    gdp_growth_discrepancy_pp
  ) %>%
  # Only region-years present in the DOSE regression sample can actually be
  # dropped; the rest are ranked outliers the regression never used.
  left_join(
    estimation_sample %>%
      distinct(GID_0, GID_1, year) %>%
      mutate(in_estimation_sample = TRUE),
    by = c("GID_0", "GID_1", "year")
  ) %>%
  mutate(in_estimation_sample = coalesce(in_estimation_sample, FALSE))

# Restrict after producing the original attribution ranking. Thus each removal
# count corresponds to that many individual region-year observations actually
# available to the subnational regression--never an entire country.
ranked_estimation_region_year_outliers <- ranked_region_year_outliers %>%
  filter(in_estimation_sample) %>%
  arrange(outlier_rank) %>%
  mutate(removal_rank = row_number())

removal_counts <- OUTLIER_REMOVAL_COUNTS[
  OUTLIER_REMOVAL_COUNTS <= nrow(ranked_estimation_region_year_outliers)
]
if (length(removal_counts) < length(OUTLIER_REMOVAL_COUNTS)) {
  warning(
    "Only ", nrow(ranked_estimation_region_year_outliers),
    " ranked estimation-sample region-year outliers are available; larger ",
    "counts were skipped."
  )
}
removal_counts <- c(0L, removal_counts)

# Re-estimation across the reduced samples --------------------------------------

estimate_with_outliers_removed <- function(n_removed) {
  removed <- ranked_estimation_region_year_outliers[seq_len(n_removed), ]
  model_data <- anti_join(
    estimation_sample,
    removed %>% select(GID_0, GID_1, year),
    by = c("GID_0", "GID_1", "year")
  )

  fit <- notebook$fit_climate_model(
    SENSITIVITY_MODEL_TERMS,
    model_data = model_data
  )

  coefficients <- as.data.frame(fixest::coeftable(fit))
  names(coefficients) <- c("estimate", "std_error", "statistic", "p_value")
  intervals <- as.data.frame(
    stats::confint(fit, level = SENSITIVITY_CONFIDENCE_LEVEL)
  )[rownames(coefficients), , drop = FALSE]

  tibble::tibble(
    n_outliers_removed = as.integer(n_removed),
    outliers_in_estimation_sample = nrow(removed),
    rows_dropped = nrow(estimation_sample) - nrow(model_data),
    min_outlier_metric_removed = if (n_removed > 0L) {
      min(removed$outlier_metric)
    } else {
      NA_real_
    },
    observations = as.integer(stats::nobs(fit)),
    regions = n_distinct(model_data$GID_1),
    countries = n_distinct(model_data$GID_0),
    term = rownames(coefficients),
    estimate = coefficients$estimate,
    std_error = coefficients$std_error,
    statistic = coefficients$statistic,
    p_value = coefficients$p_value,
    conf_low = intervals[[1L]],
    conf_high = intervals[[2L]]
  )
}

outlier_removal_coefficients <- bind_rows(
  lapply(removal_counts, estimate_with_outliers_removed)
) %>%
  mutate(term = factor(term, levels = unique(term)))

outlier_removal_samples <- outlier_removal_coefficients %>%
  distinct(
    n_outliers_removed,
    outliers_in_estimation_sample,
    rows_dropped,
    min_outlier_metric_removed,
    observations,
    regions,
    countries
  )

# Signed-bin response paths -----------------------------------------------------

signed_bin_outlier_coefficients <- outlier_removal_coefficients %>%
  filter(grepl("^(TM|RR)_bin_signed::", as.character(term))) %>%
  mutate(
    weather = if_else(
      startsWith(as.character(term), "TM_bin_signed::"),
      "Temperature",
      "Precipitation"
    ),
    bin = suppressWarnings(as.numeric(sub(
      "^(TM|RR)_bin_signed::",
      "",
      as.character(term)
    )))
  )

if (!nrow(signed_bin_outlier_coefficients) ||
    anyNA(signed_bin_outlier_coefficients$bin)) {
  stop("Could not recover the signed temperature and precipitation bins.")
}

# Add the omitted central bin so every curve is displayed relative to zero.
signed_bin_reference_rows <- bind_rows(
  outlier_removal_samples %>%
    transmute(
      n_outliers_removed,
      weather = "Temperature",
      bin = 0,
      estimate = 0,
      std_error = 0,
      conf_low = 0,
      conf_high = 0
    ),
  outlier_removal_samples %>%
    transmute(
      n_outliers_removed,
      weather = "Precipitation",
      bin = 0,
      estimate = 0,
      std_error = 0,
      conf_low = 0,
      conf_high = 0
    )
)

signed_bin_outlier_plot_data <- bind_rows(
  signed_bin_outlier_coefficients,
  signed_bin_reference_rows
) %>%
  arrange(weather, n_outliers_removed, bin)

outlier_removal_plot <- ggplot2::ggplot(
  signed_bin_outlier_plot_data,
  ggplot2::aes(
    x = bin,
    y = estimate,
    colour = n_outliers_removed,
    group = n_outliers_removed
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "#52514e",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = conf_low,
      ymax = conf_high,
      fill = n_outliers_removed
    ),
    alpha = 0.08,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::facet_wrap(~weather, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(signed_bin_outlier_plot_data$bin))
  ) +
  ggplot2::scale_colour_viridis_c(
    option = "C",
    direction = -1,
    breaks = removal_counts,
    labels = format(removal_counts, big.mark = ",")
  ) +
  ggplot2::scale_fill_viridis_c(
    option = "C",
    direction = -1,
    breaks = removal_counts,
    labels = format(removal_counts, big.mark = ","),
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Signed-bin response after removing DOSE-PWT regional outliers",
    subtitle = paste(
      "Each colour removes additional individual region-years;",
      paste0(
        100 * SENSITIVITY_CONFIDENCE_LEVEL,
        "% confidence ribbons use standard errors clustered by region"
      )
    ),
    x = "Signed anomaly bin (standard deviations)",
    y = "Estimated effect relative to the central bin",
    colour = "Regional outliers removed",
    caption = paste0(
      "Sample: DOSE V2.11 subnational panel from build_dat(); outliers ranked by ",
      OUTLIER_RANK_METRIC, " from the regional attribution--not by country. ",
      "Observations fall from ",
      format(max(outlier_removal_samples$observations), big.mark = ","),
      " to ",
      format(min(outlier_removal_samples$observations), big.mark = ","),
      "."
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", colour = "#0b0b0b"),
    plot.subtitle = ggplot2::element_text(colour = "#52514e"),
    plot.caption = ggplot2::element_text(colour = "#52514e", hjust = 0),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", colour = "#0b0b0b"),
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey92",
      linewidth = 0.3
    ),
    axis.title = ggplot2::element_text(colour = "#52514e")
  )

# Outputs -----------------------------------------------------------------------

write.csv(
  signed_bin_outlier_plot_data,
  file.path(output_dir, "signed_bin_outlier_removal_coefficients.csv"),
  row.names = FALSE
)
write.csv(
  outlier_removal_samples,
  file.path(output_dir, "outlier_removal_sample_sizes.csv"),
  row.names = FALSE
)
write.csv(
  ranked_estimation_region_year_outliers %>%
    filter(removal_rank <= max(removal_counts)),
  file.path(output_dir, "ranked_region_year_outliers.csv"),
  row.names = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "signed_bin_outlier_removal_paths.png"),
  outlier_removal_plot,
  width = 11,
  height = 6.5,
  dpi = 300,
  bg = "#fcfcfb"
)

cat("\nOutlier-removal sensitivity (", SENSITIVITY_MODEL_TERMS, "):\n", sep = "")
print(outlier_removal_samples, n = nrow(outlier_removal_samples))
cat("\nSigned-bin coefficient paths:\n")
print(
  signed_bin_outlier_plot_data %>%
    select(
      n_outliers_removed,
      weather,
      bin,
      estimate,
      std_error,
      conf_low,
      conf_high,
      p_value
    ),
  n = nrow(outlier_removal_coefficients)
)
