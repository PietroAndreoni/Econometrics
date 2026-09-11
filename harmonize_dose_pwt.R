# Harmonize DOSE regional GDP with the PWT national real-GDP growth path.
#
# Two complete, row-preserving output datasets are produced:
#   1. proportional: the same country-year correction ratio for every region;
#   2. weighted_denton: reliability-weighted, movement-preserving correction
#      ratios estimated jointly over regions and years within each country.
#
# The original DOSE columns are never overwritten. Harmonized regional GDP and
# GDP per capita, benchmark metadata, and diagnostics are appended. Because
# DOSE and PWT levels use different price bases, PWT is used only as a growth
# path: for each country it is anchored to the complete DOSE national aggregate
# in 2015, or the closest overlapping complete year if 2015 is unavailable.

required_packages <- c("dplyr", "readxl", "Matrix", "quadprog")
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

preferred_anchor_year <- 2015L
fallback_identifier_reliability <- 0.5
denton_level_penalty <- 10
positive_ratio_floor <- 1e-10
validation_relative_tolerance <- 1e-9
validation_positive_tolerance <- positive_ratio_floor * 0.99

if (file.exists("RICE50x_PR.Rproj")) {
  project_root <- "."
} else if (file.exists(file.path("..", "RICE50x_PR.Rproj"))) {
  project_root <- ".."
} else {
  stop("Run this script from the project root or the econometrics directory.")
}

dose_path <- file.path(project_root, "econometrics", "data", "DOSE_V2.11.csv")
pwt_path <- file.path(project_root, "econometrics", "data", "pwt110.xlsx")
proportional_output_path <- file.path(
  project_root,
  "econometrics",
  "data",
  "DOSE_harmonized_proportional.csv"
)
denton_output_path <- file.path(
  project_root,
  "econometrics",
  "data",
  "DOSE_harmonized_weighted_denton.csv"
)
diagnostic_dir <- file.path(
  project_root,
  "econometrics",
  "results",
  "dose_pwt_harmonization"
)

missing_inputs <- c(dose_path, pwt_path)[!file.exists(c(dose_path, pwt_path))]
if (length(missing_inputs)) {
  stop("Raw input files not found: ", paste(missing_inputs, collapse = ", "))
}
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare DOSE and PWT ----------------------------------------------------------

dose_original <- read.csv(
  dose_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dose_work <- dose_original %>%
  transmute(
    dose_row_id = row_number(),
    country_dose = as.character(country),
    iso3 = toupper(trimws(coalesce(as.character(GID_0), ""))),
    region_name = trimws(coalesce(as.character(region), "")),
    region_id_file = trimws(coalesce(as.character(GID_1), "")),
    year = as.integer(year),
    regional_gdp_per_capita = as.numeric(grp_pc_usd_2015),
    regional_population = as.numeric(pop)
  ) %>%
  mutate(
    region_id = if_else(
      nzchar(region_id_file),
      region_id_file,
      paste0(iso3, "::UNMAPPED::", region_name)
    ),
    used_fallback_region_id = !nzchar(region_id_file),
    valid_component = is.finite(regional_gdp_per_capita) &
      regional_gdp_per_capita > 0 &
      is.finite(regional_population) &
      regional_population > 0,
    regional_gdp_raw = if_else(
      valid_component,
      regional_gdp_per_capita * regional_population,
      NA_real_
    ),
    # This is the only explicit row-level quality signal available in DOSE.
    # Lower weight means the observation may absorb more of the adjustment.
    reliability_weight = if_else(
      used_fallback_region_id,
      fallback_identifier_reliability,
      1
    ),
    reliability_weight_reason = if_else(
      used_fallback_region_id,
      "fallback_region_identifier",
      "default_equal_reliability"
    )
  )

invalid_dose_identifiers <- dose_work %>%
  filter(
    is.na(year) |
      !grepl("^[A-Z]{3}$", iso3) |
      (!nzchar(region_id_file) & !nzchar(region_name))
  )
if (nrow(invalid_dose_identifiers)) {
  stop(
    "DOSE contains missing years, invalid country codes, or rows with neither ",
    "a region ID nor a region name."
  )
}

duplicate_dose_keys <- dose_work %>%
  count(iso3, region_id, year, name = "rows") %>%
  filter(rows > 1L)
if (nrow(duplicate_dose_keys)) {
  write.csv(
    duplicate_dose_keys,
    file.path(diagnostic_dir, "duplicate_dose_region_year_keys.csv"),
    row.names = FALSE
  )
  stop("DOSE has duplicate region-year keys; see the diagnostics directory.")
}

dose_country_year <- dose_work %>%
  group_by(iso3, year) %>%
  summarise(
    n_regions_reported = n_distinct(region_id),
    n_regions_complete = n_distinct(region_id[valid_component]),
    complete_dose_country_year = all(valid_component),
    national_gdp_dose_raw = if_else(
      complete_dose_country_year,
      sum(regional_gdp_raw),
      NA_real_
    ),
    .groups = "drop"
  )

pwt_raw <- readxl::read_excel(pwt_path, sheet = "Data") %>%
  transmute(
    iso3 = toupper(trimws(coalesce(as.character(countrycode), ""))),
    country_pwt = as.character(country),
    year = as.integer(year),
    pwt_rgdpna = as.numeric(rgdpna)
  )

duplicate_pwt_keys <- pwt_raw %>%
  filter(nzchar(iso3), !is.na(year)) %>%
  count(iso3, year, name = "rows") %>%
  filter(rows > 1L)
if (nrow(duplicate_pwt_keys)) {
  write.csv(
    duplicate_pwt_keys,
    file.path(diagnostic_dir, "duplicate_pwt_country_year_keys.csv"),
    row.names = FALSE
  )
  stop("PWT has duplicate country-year keys; see the diagnostics directory.")
}

pwt_national <- pwt_raw %>%
  filter(
    nzchar(iso3),
    !is.na(year),
    is.finite(pwt_rgdpna),
    pwt_rgdpna > 0
  )

# Construct a PWT growth benchmark in the level units of DOSE. -----------------

anchor_candidates <- dose_country_year %>%
  filter(complete_dose_country_year, national_gdp_dose_raw > 0) %>%
  inner_join(pwt_national, by = c("iso3", "year")) %>%
  arrange(
    iso3,
    abs(year - preferred_anchor_year),
    year
  ) %>%
  group_by(iso3) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  transmute(
    iso3,
    harmonization_anchor_year = year,
    anchor_national_gdp_dose_raw = national_gdp_dose_raw,
    anchor_pwt_rgdpna = pwt_rgdpna
  )

country_year_benchmarks <- dose_country_year %>%
  left_join(pwt_national, by = c("iso3", "year")) %>%
  left_join(anchor_candidates, by = "iso3") %>%
  mutate(
    national_gdp_benchmark = if_else(
      complete_dose_country_year &
        is.finite(national_gdp_dose_raw) &
        national_gdp_dose_raw > 0 &
        is.finite(pwt_rgdpna) &
        pwt_rgdpna > 0 &
        is.finite(anchor_national_gdp_dose_raw) &
        anchor_national_gdp_dose_raw > 0 &
        is.finite(anchor_pwt_rgdpna) &
        anchor_pwt_rgdpna > 0,
      anchor_national_gdp_dose_raw * pwt_rgdpna / anchor_pwt_rgdpna,
      NA_real_
    ),
    national_benchmark_ratio = national_gdp_benchmark /
      national_gdp_dose_raw,
    harmonizable_country_year =
      is.finite(national_gdp_benchmark) & national_gdp_benchmark > 0,
    harmonization_status = case_when(
      !complete_dose_country_year ~ "incomplete_dose_country_year",
      !is.finite(pwt_rgdpna) | pwt_rgdpna <= 0 ~ "pwt_unavailable",
      !is.finite(anchor_pwt_rgdpna) | anchor_pwt_rgdpna <= 0 ~
        "no_valid_country_anchor",
      harmonizable_country_year ~ "harmonized",
      TRUE ~ "not_harmonizable"
    )
  )

dose_model <- dose_work %>%
  left_join(
    country_year_benchmarks,
    by = c("iso3", "year")
  ) %>%
  arrange(dose_row_id)

# Method 1: proportional scaling -----------------------------------------------

proportional_ratios <- dose_model %>%
  transmute(
    dose_row_id,
    harmonization_ratio = if_else(
      harmonizable_country_year,
      national_benchmark_ratio,
      NA_real_
    ),
    solver = if_else(
      harmonizable_country_year,
      "closed_form_proportional",
      NA_character_
    )
  )

# Method 2: reliability-weighted spatial-temporal proportional Denton ----------

# Let z_it be the harmonized-to-raw GDP ratio. Within each country, solve
#
#   min sum_edges w_it (z_it - z_i,t-1)^2
#       + lambda sum_it w_it (z_it - k_t)^2
#
# subject to sum_i share_it * z_it = benchmark_t / raw_total_t for every year.
# Here k_t is the proportional country-year ratio. It is a neutral allocation
# prior: departures from proportional scaling are allowed only when they make
# regional correction ratios smoother over time. The level penalty also makes
# the system strictly convex and handles regions observed for only one year.

solve_weighted_denton_country <- function(country_data) {
  country_data <- country_data %>%
    arrange(year, region_id)

  n_values <- nrow(country_data)
  years <- sort(unique(country_data$year))
  n_years <- length(years)
  value_index <- seq_len(n_values)
  year_index <- match(country_data$year, years)

  regional_share <- country_data$regional_gdp_raw /
    country_data$national_gdp_dose_raw
  target_ratio <- country_data %>%
    distinct(year, national_benchmark_ratio) %>%
    arrange(year) %>%
    pull(national_benchmark_ratio)
  proportional_center <- target_ratio[year_index]

  aggregation_matrix <- Matrix::sparseMatrix(
    i = year_index,
    j = value_index,
    x = regional_share,
    dims = c(n_years, n_values)
  )

  region_year_key <- paste(country_data$region_id, country_data$year, sep = "\r")
  previous_key <- paste(
    country_data$region_id,
    country_data$year - 1L,
    sep = "\r"
  )
  previous_index <- match(previous_key, region_year_key)
  current_index <- which(!is.na(previous_index))
  previous_index <- previous_index[current_index]

  penalty_matrix <- Matrix::Diagonal(
    n_values,
    x = denton_level_penalty * country_data$reliability_weight
  )
  if (length(current_index)) {
    difference_matrix <- Matrix::sparseMatrix(
      i = rep(seq_along(current_index), 2L),
      j = c(current_index, previous_index),
      x = c(
        rep(1, length(current_index)),
        rep(-1, length(current_index))
      ),
      dims = c(length(current_index), n_values)
    )
    edge_weight <- sqrt(
      country_data$reliability_weight[current_index] *
        country_data$reliability_weight[previous_index]
    )
    penalty_matrix <- penalty_matrix +
      Matrix::crossprod(
        difference_matrix,
        Matrix::Diagonal(x = edge_weight) %*% difference_matrix
      )
  }

  level_linear_term <- denton_level_penalty *
    country_data$reliability_weight * proportional_center

  # Equality-constrained sparse quadratic solution.
  unconstrained_ratio <- as.numeric(Matrix::solve(
    penalty_matrix,
    level_linear_term
  ))
  inverse_q_at <- Matrix::solve(
    penalty_matrix,
    Matrix::t(aggregation_matrix)
  )
  constraint_covariance <- as.matrix(
    aggregation_matrix %*% inverse_q_at
  )
  multiplier <- solve(
    constraint_covariance,
    target_ratio - as.numeric(aggregation_matrix %*% unconstrained_ratio)
  )
  ratio <- as.numeric(unconstrained_ratio + inverse_q_at %*% multiplier)
  solver <- "sparse_equality"

  # The linear Denton solution is normally positive. If it is not, enforce the
  # lower bound with a dense quadratic-program fallback for that country only.
  if (any(!is.finite(ratio)) || any(ratio < positive_ratio_floor)) {
    dense_penalty <- as.matrix(2 * penalty_matrix)
    equality_constraints <- t(as.matrix(aggregation_matrix))
    inequality_constraints <- diag(n_values)
    qp_result <- quadprog::solve.QP(
      Dmat = dense_penalty,
      dvec = as.numeric(2 * level_linear_term),
      Amat = cbind(equality_constraints, inequality_constraints),
      bvec = c(target_ratio, rep(positive_ratio_floor, n_values)),
      meq = n_years
    )
    ratio <- as.numeric(qp_result$solution)
    solver <- "quadprog_positive"
  }

  tibble::tibble(
    dose_row_id = country_data$dose_row_id,
    harmonization_ratio = ratio,
    solver = solver
  )
}

denton_input <- dose_model %>%
  filter(harmonizable_country_year, valid_component)

denton_ratios <- denton_input %>%
  group_split(iso3) %>%
  lapply(solve_weighted_denton_country) %>%
  bind_rows()

# Build row-preserving method datasets -----------------------------------------

build_harmonized_output <- function(method_name, ratios) {
  diagnostics <- dose_model %>%
    select(
      dose_row_id,
      iso3,
      year,
      region_id,
      used_fallback_region_id,
      valid_component,
      reliability_weight,
      reliability_weight_reason,
      n_regions_reported,
      n_regions_complete,
      complete_dose_country_year,
      country_pwt,
      pwt_rgdpna,
      harmonization_anchor_year,
      anchor_national_gdp_dose_raw,
      anchor_pwt_rgdpna,
      national_gdp_dose_raw,
      national_gdp_benchmark,
      national_benchmark_ratio,
      harmonizable_country_year,
      harmonization_status,
      regional_gdp_raw,
      regional_population
    ) %>%
    left_join(ratios, by = "dose_row_id") %>%
    arrange(dose_row_id) %>%
    mutate(
      harmonization_method = method_name,
      regional_gdp_harmonized = if_else(
        harmonizable_country_year & valid_component,
        regional_gdp_raw * harmonization_ratio,
        NA_real_
      ),
      grp_pc_usd_2015_harmonized = if_else(
        is.finite(regional_gdp_harmonized) &
          is.finite(regional_population) & regional_population > 0,
        regional_gdp_harmonized / regional_population,
        NA_real_
      ),
      regional_gdp_adjustment =
        regional_gdp_harmonized - regional_gdp_raw,
      regional_gdp_adjustment_pct =
        100 * (harmonization_ratio - 1),
      harmonization_ratio_relative_to_proportional =
        harmonization_ratio / national_benchmark_ratio,
      anchor_is_preferred_year =
        harmonization_anchor_year == preferred_anchor_year
    ) %>%
    transmute(
      dose_row_id,
      harmonization_method,
      harmonization_status,
      harmonizable_country_year,
      harmonization_anchor_year,
      anchor_is_preferred_year,
      pwt_country = country_pwt,
      pwt_rgdpna_million_2021_usd = pwt_rgdpna,
      anchor_national_gdp_dose_2015_usd = anchor_national_gdp_dose_raw,
      anchor_pwt_rgdpna_million_2021_usd = anchor_pwt_rgdpna,
      national_gdp_dose_2015_usd_raw = national_gdp_dose_raw,
      national_gdp_benchmark_dose_2015_usd = national_gdp_benchmark,
      national_benchmark_ratio,
      region_id_harmonization = region_id,
      used_fallback_region_id,
      valid_regional_gdp_component = valid_component,
      complete_dose_country_year,
      n_regions_reported,
      n_regions_complete,
      reliability_weight,
      reliability_weight_reason,
      harmonization_ratio,
      harmonization_ratio_relative_to_proportional,
      harmonization_solver = solver,
      denton_proportional_prior_penalty = if_else(
        method_name == "weighted_spatiotemporal_denton",
        denton_level_penalty,
        NA_real_
      ),
      regional_gdp_dose_2015_usd_raw = regional_gdp_raw,
      regional_gdp_dose_2015_usd_harmonized = regional_gdp_harmonized,
      grp_pc_usd_2015_harmonized,
      regional_gdp_adjustment_dose_2015_usd = regional_gdp_adjustment,
      regional_gdp_adjustment_pct
    )

  if (nrow(diagnostics) != nrow(dose_original) ||
      !identical(diagnostics$dose_row_id, seq_len(nrow(dose_original)))) {
    stop(method_name, " output no longer has a one-to-one DOSE row mapping.")
  }

  bind_cols(dose_original, diagnostics)
}

proportional_output <- build_harmonized_output(
  "proportional",
  proportional_ratios
)
denton_output <- build_harmonized_output(
  "weighted_spatiotemporal_denton",
  denton_ratios
)

# Validation -------------------------------------------------------------------

validate_method <- function(output_data, method_name) {
  validation <- output_data %>%
    filter(harmonizable_country_year) %>%
    group_by(GID_0, year) %>%
    summarise(
      method = first(harmonization_method),
      expected_regions = first(n_regions_complete),
      harmonized_regions = sum(is.finite(
        regional_gdp_dose_2015_usd_harmonized
      )),
      national_benchmark = first(
        national_gdp_benchmark_dose_2015_usd
      ),
      harmonized_region_sum = sum(
        regional_gdp_dose_2015_usd_harmonized
      ),
      benchmark_error = harmonized_region_sum - national_benchmark,
      relative_benchmark_error = benchmark_error / national_benchmark,
      minimum_harmonization_ratio = min(harmonization_ratio),
      maximum_harmonization_ratio = max(harmonization_ratio),
      minimum_ratio_relative_to_proportional = min(
        harmonization_ratio_relative_to_proportional
      ),
      maximum_ratio_relative_to_proportional = max(
        harmonization_ratio_relative_to_proportional
      ),
      .groups = "drop"
    )

  invalid <- validation %>%
    filter(
      harmonized_regions != expected_regions |
        !is.finite(relative_benchmark_error) |
        abs(relative_benchmark_error) > validation_relative_tolerance |
        !is.finite(minimum_harmonization_ratio) |
        minimum_harmonization_ratio < validation_positive_tolerance
    )
  invalid_path <- file.path(
    diagnostic_dir,
    paste0("invalid_", method_name, "_harmonization.csv")
  )
  if (nrow(invalid)) {
    write.csv(
      invalid,
      invalid_path,
      row.names = FALSE
    )
    stop(method_name, " harmonization failed its accounting checks.")
  }
  if (file.exists(invalid_path)) {
    unlink(invalid_path)
  }

  validation
}

proportional_validation <- validate_method(
  proportional_output,
  "proportional"
)
denton_validation <- validate_method(
  denton_output,
  "weighted_denton"
)
harmonization_validation <- bind_rows(
  proportional_validation,
  denton_validation
)

temporal_movement_validation <- bind_rows(
  proportional_output,
  denton_output
) %>%
  filter(harmonizable_country_year) %>%
  arrange(harmonization_method, GID_0, region_id_harmonization, year) %>%
  group_by(harmonization_method, GID_0, region_id_harmonization) %>%
  mutate(
    consecutive_year = year == lag(year) + 1L,
    correction_ratio_change = if_else(
      consecutive_year,
      harmonization_ratio - lag(harmonization_ratio),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  group_by(harmonization_method) %>%
  summarise(
    consecutive_region_year_edges = sum(is.finite(correction_ratio_change)),
    correction_ratio_roughness = sum(
      correction_ratio_change^2,
      na.rm = TRUE
    ),
    mean_absolute_correction_ratio_change = mean(
      abs(correction_ratio_change),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

roughness_proportional <- temporal_movement_validation %>%
  filter(harmonization_method == "proportional") %>%
  pull(correction_ratio_roughness)
roughness_denton <- temporal_movement_validation %>%
  filter(harmonization_method == "weighted_spatiotemporal_denton") %>%
  pull(correction_ratio_roughness)
if (length(roughness_proportional) != 1L || length(roughness_denton) != 1L ||
    !is.finite(roughness_proportional) || !is.finite(roughness_denton) ||
    roughness_denton > roughness_proportional * (1 + 1e-10)) {
  stop("Weighted Denton did not improve correction-ratio temporal smoothness.")
}

validate_national_growth <- function(output_data) {
  output_data %>%
    filter(harmonizable_country_year) %>%
    group_by(harmonization_method, GID_0, year) %>%
    summarise(
      harmonized_national_gdp = sum(
        regional_gdp_dose_2015_usd_harmonized
      ),
      pwt_rgdpna = first(pwt_rgdpna_million_2021_usd),
      .groups = "drop"
    ) %>%
    group_by(harmonization_method, GID_0) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      consecutive_year = year == lag(year) + 1L,
      harmonized_national_growth_pct = if_else(
        consecutive_year,
        100 * (harmonized_national_gdp / lag(harmonized_national_gdp) - 1),
        NA_real_
      ),
      pwt_national_growth_pct = if_else(
        consecutive_year,
        100 * (pwt_rgdpna / lag(pwt_rgdpna) - 1),
        NA_real_
      ),
      national_growth_error_pp =
        harmonized_national_growth_pct - pwt_national_growth_pct
    ) %>%
    ungroup() %>%
    filter(consecutive_year)
}

national_growth_validation <- bind_rows(
  validate_national_growth(proportional_output),
  validate_national_growth(denton_output)
)
if (!nrow(national_growth_validation) ||
    any(!is.finite(national_growth_validation$national_growth_error_pp)) ||
    max(abs(national_growth_validation$national_growth_error_pp)) > 1e-9) {
  stop("Harmonized national growth does not reproduce the PWT growth path.")
}

# Write outputs only after both methods pass all checks. ------------------------

write.csv(proportional_output, proportional_output_path, row.names = FALSE)
write.csv(denton_output, denton_output_path, row.names = FALSE)
write.csv(
  harmonization_validation,
  file.path(diagnostic_dir, "harmonization_validation.csv"),
  row.names = FALSE
)
write.csv(
  anchor_candidates,
  file.path(diagnostic_dir, "country_anchor_years.csv"),
  row.names = FALSE
)
write.csv(
  temporal_movement_validation,
  file.path(diagnostic_dir, "temporal_movement_validation.csv"),
  row.names = FALSE
)
write.csv(
  national_growth_validation,
  file.path(diagnostic_dir, "national_growth_validation.csv"),
  row.names = FALSE
)

method_summary <- harmonization_validation %>%
  group_by(method) %>%
  summarise(
    country_years = n(),
    countries = n_distinct(GID_0),
    maximum_absolute_relative_error = max(abs(relative_benchmark_error)),
    minimum_ratio = min(minimum_harmonization_ratio),
    maximum_ratio = max(maximum_harmonization_ratio),
    minimum_ratio_relative_to_proportional = min(
      minimum_ratio_relative_to_proportional
    ),
    maximum_ratio_relative_to_proportional = max(
      maximum_ratio_relative_to_proportional
    ),
    .groups = "drop"
  )

write.csv(
  method_summary,
  file.path(diagnostic_dir, "method_summary.csv"),
  row.names = FALSE
)

cat(
  "Harmonized ", nrow(proportional_validation),
  " complete country-year totals across ",
  n_distinct(proportional_validation$GID_0), " countries.\n",
  sep = ""
)
print(method_summary)
cat("\nWrote:\n")
cat("  ", normalizePath(proportional_output_path), "\n", sep = "")
cat("  ", normalizePath(denton_output_path), "\n", sep = "")
