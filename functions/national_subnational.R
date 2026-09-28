# Compare the national GDP growth implied by a subnational dataset with a
# national dataset, and attribute every country-year discrepancy to the
# subnational regions that produce it.
#
# Both inputs follow the read_gdp_levels() layout (GADM identifiers):
#   subnational: GID_0, GID_1, region_name, country, year, gdp, population
#   national:    GID_0, country, year, gdp, population
# Regional GDP is summed by country-year. The two sources may use different
# units and price bases; only annual growth rates are compared.
#
#   comparison <- compare_national_subnational(read_gdp_levels("DOSE"),
#                                              read_gdp_levels("PWT"))
#   comparison$ranked_outliers   # region-years, largest contribution first

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

.stop_on_duplicates <- function(data, keys, what) {
  duplicates <- data %>%
    count(across(all_of(keys)), name = "rows") %>%
    filter(rows > 1L)
  if (nrow(duplicates)) {
    examples <- utils::head(
      do.call(paste, c(duplicates[keys], sep = " ")),
      5L
    )
    stop(
      what, " has ", nrow(duplicates), " duplicate ",
      paste(keys, collapse = "-"), " keys, e.g. ",
      paste(examples, collapse = "; ")
    )
  }
}

# Regions with a name but no GID_1 are kept under a stable fallback key rather
# than dropped, which would silently bias the national aggregate. A component is
# valid only with positive GDP and population.
.standardize_subnational <- function(subnational) {
  subnational <- subnational %>%
    transmute(
      country = as.character(country),
      iso3 = toupper(trimws(coalesce(as.character(GID_0), ""))),
      region_name = trimws(coalesce(as.character(region_name), "")),
      region_id_file = trimws(coalesce(as.character(GID_1), "")),
      year = as.integer(year),
      regional_population = as.numeric(population),
      regional_gdp_raw = as.numeric(gdp)
    ) %>%
    mutate(
      region_id = if_else(
        nzchar(region_id_file),
        region_id_file,
        paste0(iso3, "::UNMAPPED::", region_name)
      ),
      region_id_prefix = substr(region_id_file, 1L, 3L),
      used_fallback_region_id = !nzchar(region_id_file),
      region_id_country_mismatch = nzchar(region_id_file) &
        region_id_prefix != iso3,
      valid_component = is.finite(regional_gdp_raw) &
        regional_gdp_raw > 0 &
        is.finite(regional_population) &
        regional_population > 0,
      regional_gdp = if_else(valid_component, regional_gdp_raw, NA_real_)
    ) %>%
    select(-regional_gdp_raw)

  invalid <- subnational %>%
    filter(
      is.na(year) |
        !grepl("^[A-Z0-9]{3}$", iso3) |
        (!nzchar(region_id_file) & !nzchar(region_name))
    )
  if (nrow(invalid)) {
    stop(
      "The subnational data contain ", nrow(invalid), " row(s) with a missing ",
      "year, an invalid GID_0, or neither a GID_1 nor a region name."
    )
  }
  .stop_on_duplicates(subnational, c("iso3", "region_id", "year"),
                      "The subnational data")
  subnational
}

.standardize_national <- function(national) {
  national <- national %>%
    transmute(
      iso3 = toupper(trimws(coalesce(as.character(GID_0), ""))),
      country_nat = as.character(country),
      year = as.integer(year),
      gdp_nat = as.numeric(gdp),
      population_nat = as.numeric(population)
    ) %>%
    filter(
      nzchar(iso3),
      !is.na(year),
      is.finite(gdp_nat),
      gdp_nat > 0,
      is.finite(population_nat),
      population_nat > 0
    )
  .stop_on_duplicates(national, c("iso3", "year"), "The national data")
  national
}

compare_national_subnational <- function(
    subnational,
    national,
    flag_percentile = 0.99,
    top_n = 50L,
    tolerance_pp = 1e-8,
    outlier_metric = "abs_additive_contribution_pp"
) {
  sub <- .standardize_subnational(subnational)
  nat_levels <- .standardize_national(national)

  identifier_issues <- sub %>%
    filter(used_fallback_region_id | region_id_country_mismatch) %>%
    distinct(
      country, iso3, region_name, region_id_file, region_id,
      used_fallback_region_id, region_id_country_mismatch
    ) %>%
    arrange(iso3, region_id)

  # National aggregate of the subnational data. A growth rate is eligible only
  # when the same complete set of regions is reported in both years.
  sub_national <- sub %>%
    group_by(iso3, year) %>%
    summarise(
      country_sub = first(country),
      n_regions_reported = n_distinct(region_id),
      n_regions_complete = n_distinct(region_id[valid_component]),
      n_incomplete_regions = n_regions_reported - n_regions_complete,
      region_signature = paste(
        sort(unique(region_id[valid_component])),
        collapse = "|"
      ),
      population_sub = sum(regional_population[valid_component]),
      gdp_sub = sum(regional_gdp[valid_component]),
      .groups = "drop"
    ) %>%
    filter(n_regions_complete > 0L, population_sub > 0, gdp_sub > 0) %>%
    mutate(
      gdppc_sub = gdp_sub / population_sub,
      complete_regional_coverage = n_incomplete_regions == 0L
    ) %>%
    group_by(iso3) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      consecutive_sub_year = year == dplyr::lag(year) + 1L,
      gdp_sub_previous = lag_if_consecutive(gdp_sub, year),
      stable_region_coverage = consecutive_sub_year &
        region_signature == dplyr::lag(region_signature) &
        complete_regional_coverage &
        dplyr::lag(complete_regional_coverage, default = FALSE),
      gdp_growth_log_sub = annual_log_growth(gdp_sub, year),
      gdppc_growth_log_sub = annual_log_growth(gdppc_sub, year),
      population_growth_log_sub = annual_log_growth(population_sub, year)
    ) %>%
    ungroup()

  nat <- nat_levels %>%
    mutate(gdppc_nat = gdp_nat / population_nat) %>%
    group_by(iso3) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      consecutive_nat_year = year == dplyr::lag(year) + 1L,
      gdp_nat_previous = lag_if_consecutive(gdp_nat, year),
      gdp_growth_log_nat = annual_log_growth(gdp_nat, year),
      gdppc_growth_log_nat = annual_log_growth(gdppc_nat, year),
      population_growth_log_nat = annual_log_growth(population_nat, year)
    ) %>%
    ungroup()

  # Country-year growth comparison and discrepancy flags ----------------------

  growth_comparison <- sub_national %>%
    select(-region_signature) %>%
    inner_join(nat, by = c("iso3", "year")) %>%
    mutate(
      eligible_for_comparison = consecutive_sub_year &
        consecutive_nat_year &
        stable_region_coverage &
        is.finite(gdp_growth_log_sub) &
        is.finite(gdp_growth_log_nat),
      gdp_growth_pct_sub = 100 * expm1(gdp_growth_log_sub),
      gdp_growth_pct_nat = 100 * expm1(gdp_growth_log_nat),
      gdp_growth_discrepancy_pp = gdp_growth_pct_sub - gdp_growth_pct_nat,
      abs_gdp_growth_discrepancy_pp = abs(gdp_growth_discrepancy_pp),
      gdppc_growth_pct_sub = 100 * expm1(gdppc_growth_log_sub),
      gdppc_growth_pct_nat = 100 * expm1(gdppc_growth_log_nat),
      gdppc_growth_discrepancy_pp = gdppc_growth_pct_sub - gdppc_growth_pct_nat,
      population_growth_pct_sub = 100 * expm1(population_growth_log_sub),
      population_growth_pct_nat = 100 * expm1(population_growth_log_nat),
      population_growth_discrepancy_pp =
        population_growth_pct_sub - population_growth_pct_nat
    )

  eligible_discrepancies <- growth_comparison$abs_gdp_growth_discrepancy_pp[
    growth_comparison$eligible_for_comparison
  ]
  if (!length(eligible_discrepancies)) {
    stop("No eligible overlapping country-year growth observations were found.")
  }
  flag_threshold_pp <- unname(stats::quantile(
    eligible_discrepancies,
    probs = flag_percentile,
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
      flag_top_n = eligible_for_comparison & discrepancy_rank <= top_n
    )

  # Subregional attribution ---------------------------------------------------
  #
  # Exact shift-share decomposition. For country c, year t, let Y_it be
  # subnational GDP in region i, Y_t = sum_i Y_it, and q_nat the national GDP
  # growth factor. Region i's contribution in percentage points is
  #
  #   100 * (Y_it - q_nat * Y_i,t-1) / Y_t-1
  # = lagged GDP share_i * (regional growth % - national growth %).
  #
  # Summing over regions gives exactly 100 * [(Y_t / Y_t-1) - q_nat], the
  # reported discrepancy. Contributions in the direction of the national
  # discrepancy widen it; those with the opposite sign offset it. A
  # leave-one-region-out calculation is a non-additive sensitivity check.

  regional_growth <- sub %>%
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
      iso3, year, previous_year, country_sub, country_nat, discrepancy_rank,
      n_regions_complete, gdp_sub, gdp_sub_previous, gdp_growth_pct_sub,
      gdp_growth_pct_nat, gdp_growth_discrepancy_pp,
      abs_gdp_growth_discrepancy_pp
    )

  subregion_attribution <- attribution_targets %>%
    inner_join(
      regional_growth %>%
        select(
          iso3, year, region_id, region_name, used_fallback_region_id,
          regional_gdp, regional_gdp_previous, regional_growth_pct
        ),
      by = c("iso3", "year")
    ) %>%
    mutate(
      lagged_gdp_share = regional_gdp_previous / gdp_sub_previous,
      additive_contribution_pp = lagged_gdp_share *
        (regional_growth_pct - gdp_growth_pct_nat),
      abs_additive_contribution_pp = abs(additive_contribution_pp),
      discrepancy_aligned_contribution_pp = sign(gdp_growth_discrepancy_pp) *
        additive_contribution_pp,
      contribution_role = case_when(
        additive_contribution_pp * gdp_growth_discrepancy_pp > 0 ~ "widens",
        additive_contribution_pp * gdp_growth_discrepancy_pp < 0 ~ "offsets",
        TRUE ~ "neutral"
      ),
      share_of_national_discrepancy = if_else(
        abs(gdp_growth_discrepancy_pp) > tolerance_pp,
        additive_contribution_pp / gdp_growth_discrepancy_pp,
        NA_real_
      ),
      gdp_sub_without_region = gdp_sub - regional_gdp,
      gdp_sub_previous_without_region =
        gdp_sub_previous - regional_gdp_previous,
      leave_one_out_growth_pct_sub = if_else(
        gdp_sub_without_region > 0 & gdp_sub_previous_without_region > 0,
        100 * (gdp_sub_without_region / gdp_sub_previous_without_region - 1),
        NA_real_
      ),
      leave_one_out_discrepancy_pp =
        leave_one_out_growth_pct_sub - gdp_growth_pct_nat,
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
      absolute_contribution_rank = row_number(desc(abs_additive_contribution_pp)),
      leave_one_out_reduction_rank = row_number(
        desc(leave_one_out_abs_discrepancy_reduction_pp)
      )
    ) %>%
    ungroup()

  # The contributions must add up to each national discrepancy and the lagged
  # shares to one; anything else means the aggregation is inconsistent.
  attribution_validation <- attribution_targets %>%
    transmute(
      iso3, year, discrepancy_rank,
      expected_regions = n_regions_complete,
      national_discrepancy_pp = gdp_growth_discrepancy_pp
    ) %>%
    left_join(
      subregion_attribution %>%
        group_by(iso3, year) %>%
        summarise(
          attributed_regions = n(),
          sum_additive_contribution_pp = sum(additive_contribution_pp),
          lagged_share_sum = sum(lagged_gdp_share),
          .groups = "drop"
        ),
      by = c("iso3", "year")
    ) %>%
    mutate(
      attributed_regions = coalesce(attributed_regions, 0L),
      additive_closure_error_pp =
        sum_additive_contribution_pp - national_discrepancy_pp
    )

  invalid_attributions <- attribution_validation %>%
    filter(
      attributed_regions != expected_regions |
        !is.finite(additive_closure_error_pp) |
        abs(additive_closure_error_pp) > tolerance_pp |
        !is.finite(lagged_share_sum) |
        abs(lagged_share_sum - 1) > tolerance_pp
    )
  if (nrow(invalid_attributions)) {
    stop(
      "Subregional attribution did not reconcile with the national comparison ",
      "for ", nrow(invalid_attributions), " country-year(s), e.g. ",
      paste(
        utils::head(paste(invalid_attributions$iso3, invalid_attributions$year), 5L),
        collapse = "; "
      )
    )
  }

  top_driver <- function(rank_col, ...) {
    subregion_attribution %>%
      filter(.data[[rank_col]] == 1L) %>%
      transmute(iso3, year, ...)
  }
  country_year_subregion_driver <- attribution_targets %>%
    left_join(
      top_driver(
        "discrepancy_driver_rank",
        primary_driver_region_id = region_id,
        primary_driver_region_name = region_name,
        primary_driver_contribution_pp = additive_contribution_pp,
        primary_driver_share_of_discrepancy = share_of_national_discrepancy,
        primary_driver_regional_growth_pct = regional_growth_pct
      ),
      by = c("iso3", "year")
    ) %>%
    left_join(
      top_driver(
        "absolute_contribution_rank",
        largest_absolute_region_id = region_id,
        largest_absolute_region_name = region_name,
        largest_absolute_contribution_pp = additive_contribution_pp,
        largest_absolute_contribution_role = contribution_role
      ),
      by = c("iso3", "year")
    ) %>%
    left_join(
      top_driver(
        "leave_one_out_reduction_rank",
        leave_one_out_driver_region_id = region_id,
        leave_one_out_driver_region_name = region_name,
        leave_one_out_abs_discrepancy_reduction_pp,
        leave_one_out_discrepancy_pp,
        leave_one_out_improves_fit
      ),
      by = c("iso3", "year")
    ) %>%
    left_join(
      attribution_validation %>%
        select(iso3, year, sum_additive_contribution_pp, additive_closure_error_pp),
      by = c("iso3", "year")
    ) %>%
    arrange(discrepancy_rank)

  country_summary <- growth_comparison %>%
    filter(eligible_for_comparison) %>%
    group_by(iso3) %>%
    summarise(
      country_sub = first(country_sub),
      country_nat = first(country_nat),
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
    sub_national %>% distinct(iso3) %>% mutate(in_subnational = TRUE),
    nat %>% distinct(iso3) %>% mutate(in_national = TRUE),
    by = "iso3"
  ) %>%
    mutate(
      in_subnational = coalesce(in_subnational, FALSE),
      in_national = coalesce(in_national, FALSE)
    ) %>%
    arrange(iso3)

  comparison <- list(
    growth_comparison = growth_comparison,
    flag_threshold_pp = flag_threshold_pp,
    subregion_attribution = subregion_attribution,
    country_year_subregion_driver = country_year_subregion_driver,
    attribution_validation = attribution_validation,
    country_summary = country_summary,
    source_coverage = source_coverage,
    identifier_issues = identifier_issues
  )
  comparison$ranked_outliers <- rank_subnational_outliers(
    comparison,
    metric = outlier_metric
  )
  comparison
}

# Every eligible region-year, ranked from the largest to the smallest value of
# `metric`:
#   "abs_additive_contribution_pp": magnitude of the region's exact
#     shift-share contribution to its country's growth discrepancy;
#   "leave_one_out_abs_discrepancy_reduction_pp": how much dropping the region
#     shrinks that discrepancy.
rank_subnational_outliers <- function(
    comparison,
    metric = "abs_additive_contribution_pp"
) {
  attribution <- comparison$subregion_attribution
  if (!metric %in% names(attribution)) {
    stop("metric '", metric, "' is not a column of the subregion attribution.")
  }

  attribution %>%
    filter(is.finite(.data[[metric]])) %>%
    arrange(desc(.data[[metric]]), iso3, region_id, year) %>%
    transmute(
      outlier_rank = row_number(),
      GID_0 = iso3,
      GID_1 = region_id,
      year = as.integer(year),
      region_name,
      outlier_metric = .data[[metric]],
      additive_contribution_pp,
      contribution_role,
      regional_growth_pct,
      gdp_growth_pct_sub,
      gdp_growth_pct_nat,
      gdp_growth_discrepancy_pp,
      leave_one_out_abs_discrepancy_reduction_pp
    )
}

# Writes every table of a comparison as csv files in `output_dir`.
write_national_subnational_comparison <- function(comparison, output_dir,
                                                  top_n = 50L,
                                                  top_n_detail = 20L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(x, name) {
    utils::write.csv(x, file.path(output_dir, name), row.names = FALSE)
  }
  gc <- comparison$growth_comparison

  write(gc, "country_year_growth_comparison.csv")
  write(
    gc %>% filter(flag_high_discrepancy) %>% arrange(discrepancy_rank),
    "flagged_high_discrepancy_country_years.csv"
  )
  write(
    gc %>% filter(flag_top_n) %>% arrange(discrepancy_rank),
    paste0("top_", top_n, "_country_years.csv")
  )
  write(comparison$country_summary, "country_discrepancy_summary.csv")
  write(comparison$source_coverage, "source_country_coverage.csv")
  write(comparison$identifier_issues, "subnational_identifier_issues.csv")
  write(comparison$subregion_attribution, "subregion_discrepancy_attribution.csv")
  write(comparison$country_year_subregion_driver, "country_year_subregion_driver.csv")
  write(
    comparison$subregion_attribution %>%
      filter(discrepancy_rank <= top_n_detail) %>%
      arrange(discrepancy_rank, discrepancy_driver_rank),
    paste0("top_", top_n_detail, "_country_year_subregion_attribution.csv")
  )
  write(comparison$attribution_validation, "subregion_attribution_validation.csv")
  write(comparison$ranked_outliers, "ranked_subnational_outliers.csv")
  invisible(output_dir)
}
