# Stages 0-1 of the fast-versus-slow warming design
# (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md): freeze the protocol, record the
# data manifest, build and validate the DOSE x ERA5 panel.
#
# Economic data: DOSE v2.14, grp_pc_lcu_2015 (constant 2015 local prices),
# g = 100 * (log y_t - log y_{t-1}) for adjacent years only. Climate: Weighted
# Climate Dataset, ERA5, area weighting (the dataset's "unweighted" option is a
# grid-cell-area-weighted mean). Every climate regressor is computed on the
# climate panel before any outcome is merged (fs_add_climate_terms()).
#
# Every exclusion is a named flag; nothing is interpolated or filled.
# Outputs: results/fast_slow_warming/{data,qa}/.

source("load_functions.R")

design <- fs_design()

# Stage 0: protocol gate ---------------------------------------------------------

stage0_reasons <- character()
if (is.null(design$hash)) stop("Design hash missing: no outcome model may run.")
if (!isTRUE(design$scenario$sesoi_confirmed)) {
  if (!isTRUE(design$allow_estimation_without_equivalence_claim)) {
    stop("SESOI not affirmed and estimation without equivalence claims is not ",
         "allowed: stopping before any outcome model.")
  }
  stage0_reasons <- c(stage0_reasons, paste(
    "SESOI (1 log point) not affirmed by the study owner:",
    "effects are estimated, practical-equivalence claims are disabled"
  ))
}
if (is.null(design$prior_related_results_seen)) {
  stage0_reasons <- c(stage0_reasons, paste(
    "prior_related_results_seen is unset (never inferred): the analysis is a",
    "frozen analysis plan, not a preregistration"
  ))
}
writeLines(utils::capture.output(utils::sessionInfo()),
           fs_path("qa", "session_info.txt"))
fs_stage_status("0_protocol",
                if (length(stage0_reasons)) "PASS_WITH_WARNING" else "PASS",
                stage0_reasons, inputs = FS_DESIGN_FILE, design = design)

# Data manifest -------------------------------------------------------------------

wcd_dir <- wcd_cache_dir(create = FALSE)
wcd_files <- file.path(wcd_dir, c("gadm1_era_tmp_un__monthly.parquet",
                                  "gadm1_era_pre_un__monthly.parquet"))
gadm_files <- list.files(file.path(kummu2025_dir(), "gadm41"),
                         pattern = "\\.json$", full.names = TRUE)
dose_raw <- read_dose(data_file("DOSE_V2.14.csv"))
kg_raw <- utils::read.csv(data_file("kg_data.csv"))

manifest_row <- function(source, release, doi, path, rows, variables,
                         coverage, level, license) {
  tibble::tibble(
    source, release, doi_or_url = doi,
    access_timestamp = format(file.info(path)$mtime, "%Y-%m-%dT%H:%M:%S"),
    raw_filename = basename(path), bytes = file.size(path),
    sha256 = digest::digest(file = path, algo = "sha256"),
    rows, variables, temporal_coverage = coverage, geographic_level = level,
    license
  )
}
monthly_coverage <- function(path) {
  dates <- arrow::read_parquet(path, col_select = "Date")$Date
  paste(range(substr(dates, 2, 7)), collapse = "-")
}
manifest <- dplyr::bind_rows(
  manifest_row("DOSE", "2.14", "https://doi.org/10.5281/zenodo.20035157",
               data_file("DOSE_V2.14.csv"), nrow(dose_raw),
               paste(names(dose_raw), collapse = ";"),
               paste(range(dose_raw$year), collapse = "-"), "GADM1 (DOSE)",
               "CC BY 4.0 (per Zenodo record)"),
  lapply(wcd_files, function(p) manifest_row(
    "Weighted Climate Dataset, ERA5, area (un) weighting",
    "figshare collection v1", "https://doi.org/10.6084/m9.figshare.c.6973998.v1",
    p, 1032L, "Date; one column per GADM 4.1 ADM1 unit",
    monthly_coverage(p), "GADM 4.1 ADM1", "CC BY 4.0 (per data descriptor)"
  )),
  manifest_row("Koppen-Geiger zones by region-year", "local file", "local",
               data_file("kg_data.csv"), nrow(kg_raw),
               paste(names(kg_raw), collapse = ";"),
               paste(range(kg_raw$year), collapse = "-"), "GADM1", "unknown"),
  tibble::tibble(
    source = "GADM 4.1 ADM1 polygons (centroids, names)", release = "4.1",
    doi_or_url = "https://gadm.org", access_timestamp = NA_character_,
    raw_filename = "gadm41/*.json", bytes = sum(file.size(gadm_files)),
    sha256 = digest::digest(fs_hash_files(gadm_files), algo = "sha256"),
    rows = length(gadm_files), variables = "GID_0;GID_1;NAME_1;geometry",
    temporal_coverage = NA_character_, geographic_level = "GADM 4.1 ADM1",
    license = "GADM licence (non-commercial)"
  ),
  manifest_row("Design configuration", design$design_version, "local",
               FS_DESIGN_FILE, NA_integer_, "", "", "", "")
)
fs_write(manifest, "data", "data_manifest.csv")
design$manifest_hash <- substr(digest::digest(manifest$sha256, algo = "sha256"),
                               1, 16)

# Climate: annual ERA5 series and every climate regressor (no outcome) ----------

climate_config <- panel_config(
  climate_source = design$climate$primary_product,
  climate_weight = "area"
)
dose_ids <- unique(dose_raw$GID_1[!dose_raw$GID_1 %in% c(" ", "")])
climate <- load_climate_panel("gadm1", climate_config) %>%
  filter(GID_1 %in% dose_ids) %>%
  select(GID_0, GID_1, year, TM, RR)

# Monthly completeness: annual means must come from complete months. The
# dashboard's yearly series averages whatever months exist, so count them.
tmp_monthly <- arrow::read_parquet(wcd_files[1])
unit_cols <- setdiff(names(tmp_monthly), "Date")
unit_ids <- sub("^([^_]+)_", "\\1.", unit_cols)
keep_units <- unit_cols[unit_ids %in% dose_ids]
month_counts <- tibble::tibble(
  year = as.integer(substr(tmp_monthly$Date, 2, 5))
) %>%
  bind_cols(as.data.frame(!is.na(tmp_monthly[, keep_units]))) %>%
  group_by(year) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  tidyr::pivot_longer(-year, names_to = "unit", values_to = "months") %>%
  mutate(GID_1 = sub("^([^_]+)_", "\\1.", unit)) %>%
  select(GID_1, year, months)
rm(tmp_monthly)
complete_share_required <- design$climate$annual_completeness %||% 0.95
climate <- climate %>%
  left_join(month_counts, by = c("GID_1", "year")) %>%
  mutate(TM = if_else(!is.na(months) & months / 12 >= complete_share_required,
                      TM, NA_real_))

n_incomplete_months <- sum(climate$months < 12 &
                             climate$year <= max(dose_raw$year), na.rm = TRUE)
spec0 <- list(T0 = 0, P_mean = 0, P_sd = 1)
climate_terms <- fs_add_climate_terms(
  climate %>% select(-months), spec0,
  primary_window = design$climate$rate_window,
  rate_windows = c(design$climate$alternative_slope_windows,
                   design$climate$rate_window),
  normal_years = c(design$climate$anomaly_normal_years,
                   design$climate$alternative_anomaly_normal_years),
  temp_lags = 0:design$model$distributed_lag_max
)

# Region-level climate covariates, outcome-free: 1991-2010 baseline climate
# (path origin), 1951-1980 predetermined climate (heterogeneity) and the modal
# Koppen-Geiger family (5 = polar merged into 4, as in main_analysis.Rmd).
region_climate <- climate_terms %>%
  group_by(GID_0, GID_1) %>%
  summarise(
    B = mean(TM[year >= 1991 & year <= 2010]),
    P = mean(RR[year >= 1991 & year <= 2010]),
    T_pre = mean(TM[year >= 1951 & year <= 1980]),
    n_baseline_years = sum(!is.na(TM[year >= 1991 & year <= 2010])),
    .groups = "drop"
  )
zones <- kg_raw %>%
  filter(!is.na(mean_kg)) %>%
  mutate(zone = if_else(as.integer(mean_kg) == 5L, 4L, as.integer(mean_kg))) %>%
  count(GID_1, zone) %>%
  group_by(GID_1) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(GID_1, zone)
region_climate <- region_climate %>% left_join(zones, by = "GID_1")

# Economic panel with named flags ---------------------------------------------------

econ <- load_econ_panel("DOSE_V2_14", design$dose$outcome_pc)
econ_usd <- load_econ_panel("DOSE_V2_14", design$dose$outcome_pc_fixed_usd)

raw <- read_econ_source("DOSE_V2_14", design$dose$outcome_pc) %>%
  filter(year >= design$dose$econ_year_min)
stopifnot(!anyDuplicated(raw[c("GID_1", "year")]))    # assertion 1
raw <- raw %>%
  group_by(GID_1) %>%
  mutate(
    level_ok = !is.na(grp_pc_usd) & grp_pc_usd > 0 & !is.na(pop) & pop > 0,
    prev_level_ok = level_ok[match(year - 1L, year)],
    prev_level_ok = dplyr::coalesce(prev_level_ok, FALSE)
  ) %>%
  ungroup() %>%
  left_join(econ %>% select(GID_1, year, lgrp_pc_usd, dlgrp_pc_usd),
            by = c("GID_1", "year")) %>%
  left_join(econ_usd %>% select(GID_1, year, dl_fixed_usd = dlgrp_pc_usd),
            by = c("GID_1", "year")) %>%
  mutate(
    g = 100 * dlgrp_pc_usd,
    struct_code = dplyr::coalesce(as.integer(struct_change), 0L),
    f_bad_level = !level_ok,
    f_no_adjacent_pair = level_ok & (!prev_level_ok | is.na(g)),
    f_struct_break = level_ok & !f_no_adjacent_pair &
      struct_code %in% design$dose$struct_change_codes_excluded
  )

# Assertion 5: fixed-2015-LCU and fixed-2015-USD growth agree.
usd_gap <- with(raw, max(abs(100 * (dlgrp_pc_usd - dl_fixed_usd)),
                         na.rm = TRUE))
usd_pairs <- sum(!is.na(raw$dlgrp_pc_usd) & !is.na(raw$dl_fixed_usd))
usd_only_one <- sum(xor(is.na(raw$dlgrp_pc_usd), is.na(raw$dl_fixed_usd)))

panel <- raw %>%
  select(GID_0, GID_1, year, pop, grp_pc = grp_pc_usd, lgrp_pc_usd,
         g, struct_code, f_bad_level, f_no_adjacent_pair, f_struct_break) %>%
  left_join(climate_terms %>% select(-GID_0), by = c("GID_1", "year")) %>%
  left_join(region_climate %>% select(-GID_0), by = "GID_1") %>%
  mutate(
    f_no_climate = !f_bad_level & !f_no_adjacent_pair & !f_struct_break &
      (is.na(TM) | is.na(RR) | is.na(rate)),
    valid_growth = !f_bad_level & !f_no_adjacent_pair & !f_struct_break &
      !f_no_climate
  ) %>%
  group_by(GID_1) %>%
  mutate(n_valid = sum(valid_growth)) %>%
  ungroup() %>%
  mutate(
    f_min_obs = valid_growth & n_valid < design$dose$min_adjacent_growth_obs,
    primary = valid_growth & !f_min_obs,
    year_c = year - 1990L,
    year_c2 = year_c^2,
    decade = 10L * (year %/% 10L),
    continent = countrycode::countrycode(GID_0, "iso3c", "continent",
                                         warn = FALSE)
  )
# Region names come from the raw DOSE file.
panel$region <- dose_raw$region[match(paste(panel$GID_1, panel$year),
                                      paste(dose_raw$GID_1, dose_raw$year))]

# Sensitivity samples: add back StructChange codes 1 and/or 3 (a new raw-data
# source, an extension relative to an earlier release). Code 2 (boundary
# change) is never added back: no constant-boundary reconstruction exists.
sensitivity_flag <- function(codes_back) {
  ok <- !panel$f_bad_level & !panel$f_no_adjacent_pair &
    (!panel$f_struct_break | panel$struct_code %in% codes_back) &
    !is.na(panel$TM) & !is.na(panel$RR) & !is.na(panel$rate)
  n <- stats::ave(as.numeric(ok), panel$GID_1, FUN = sum)
  ok & n >= design$dose$min_adjacent_growth_obs
}
panel$sample_sc1 <- sensitivity_flag(1L)
panel$sample_sc3 <- sensitivity_flag(3L)
panel$sample_sc13 <- sensitivity_flag(c(1L, 3L))

# Final transformation constants over the primary rows, without outcomes.
primary_rows <- panel$primary
spec <- fs_terms_spec(
  panel$TM[primary_rows], panel$RR[primary_rows],
  spline_knot_quantiles = unlist(design$model$spline_temp_knots_quantiles),
  spline_boundary_quantiles = unlist(design$model$spline_temp_boundary_quantiles)
)
climate_terms <- fs_apply_spec(climate_terms, spec)
panel <- fs_apply_spec(panel, spec)
saveRDS(spec, fs_path("data", "terms_spec.rds"))

# Population in the fixed baseline year for the person-year estimand; regions
# without a value in that year are left out of that estimand (not imputed).
pop_base <- raw %>%
  filter(year == design$model$baseline_population_year, level_ok) %>%
  select(GID_1, pop_base = pop)
panel <- panel %>% left_join(pop_base, by = "GID_1")

# Sample flow -------------------------------------------------------------------------

flow <- tibble::tibble(
  step = c(
    "DOSE region-years (year >= econ_year_min)",
    "minus missing or nonpositive GRP pc or population",
    "minus no adjacent previous year (first year or gap)",
    "minus growth into a StructChange year (codes 1, 2, 3)",
    "minus no complete ERA5 temperature, precipitation or 5-year slope",
    "minus regions with fewer than 15 valid growth observations"
  ),
  flag = c(NA, "f_bad_level", "f_no_adjacent_pair", "f_struct_break",
           "f_no_climate", "f_min_obs"),
  removed = c(NA, sum(panel$f_bad_level), sum(panel$f_no_adjacent_pair),
              sum(panel$f_struct_break), sum(panel$f_no_climate),
              sum(panel$f_min_obs))
) %>%
  mutate(
    remaining = nrow(panel) - cumsum(dplyr::coalesce(removed, 0L)),
    regions = c(
      n_distinct(panel$GID_1),
      n_distinct(panel$GID_1[!panel$f_bad_level]),
      n_distinct(panel$GID_1[!panel$f_bad_level & !panel$f_no_adjacent_pair]),
      n_distinct(panel$GID_1[!panel$f_bad_level & !panel$f_no_adjacent_pair &
                               !panel$f_struct_break]),
      n_distinct(panel$GID_1[panel$valid_growth]),
      n_distinct(panel$GID_1[panel$primary])
    ),
    countries = c(
      n_distinct(panel$GID_0),
      n_distinct(panel$GID_0[!panel$f_bad_level]),
      n_distinct(panel$GID_0[!panel$f_bad_level & !panel$f_no_adjacent_pair]),
      n_distinct(panel$GID_0[!panel$f_bad_level & !panel$f_no_adjacent_pair &
                               !panel$f_struct_break]),
      n_distinct(panel$GID_0[panel$valid_growth]),
      n_distinct(panel$GID_0[panel$primary])
    )
  )
stopifnot(tail(flow$remaining, 1) == sum(panel$primary))
fs_write(flow, "data", "sample_flow.csv", design = design)

# Geography audit -----------------------------------------------------------------------
# The climate series are GADM 4.1 ADM1 means; DOSE is largely GADM 3.6 plus
# custom regions. Without DOSE polygons the join cannot be verified spatially:
# only identifiers and names are compared, and the gate is recorded as not met.

old_s2 <- suppressMessages(sf::sf_use_s2(FALSE))
gadm <- read_gadm1_polygons(file.path(kummu2025_dir(), "gadm41"),
                            keep_ids = NULL)
gadm_names <- sf::st_drop_geometry(gadm) %>% select(GID_1, NAME_1)
normalize_name <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  gsub("[^a-z]", "", tolower(x))
}
crosswalk <- panel %>%
  filter(!f_bad_level) %>%
  distinct(GID_0, GID_1, region) %>%
  group_by(GID_0, GID_1) %>%
  summarise(dose_name = first(region), .groups = "drop") %>%
  left_join(gadm_names, by = "GID_1") %>%
  mutate(
    climate_available = GID_1 %in% climate$GID_1,
    gadm41_polygon = !is.na(NAME_1),
    name_match = normalize_name(dose_name) == normalize_name(NAME_1),
    name_contained = mapply(function(a, b) {
      !is.na(a) && !is.na(b) && nzchar(a) && nzchar(b) &&
        (grepl(a, b, fixed = TRUE) || grepl(b, a, fixed = TRUE))
    }, normalize_name(dose_name), normalize_name(NAME_1)),
    in_primary = GID_1 %in% panel$GID_1[panel$primary],
    mapping_type = "1:1 identifier join (polygon equivalence NOT verified)",
    target_coverage = NA_real_, weight_mass = NA_real_
  ) %>%
  rename(gadm41_name = NAME_1)
fs_write(crosswalk, "data", "geography_crosswalk.csv", design = design)

centroids <- gadm[gadm$GID_1 %in% crosswalk$GID_1, ] %>% polygon_centroids()
suppressMessages(sf::sf_use_s2(old_s2))
fs_write(centroids, "data", "region_centroids.csv")

join_checks <- list(
  dose_regions = nrow(crosswalk),
  with_era5_series = sum(crosswalk$climate_available),
  with_gadm41_polygon = sum(crosswalk$gadm41_polygon),
  primary_regions = sum(crosswalk$in_primary),
  primary_exact_name_match = sum(crosswalk$name_match & crosswalk$in_primary,
                                 na.rm = TRUE),
  primary_name_match_or_contained = sum(
    (crosswalk$name_match | crosswalk$name_contained) & crosswalk$in_primary,
    na.rm = TRUE),
  polygon_equivalence_verified = FALSE,
  coverage_rule_evaluable = FALSE,
  gate = paste("NOT MET: string-key join of DOSE GID_1 to GADM 4.1 climate",
               "means; target coverage and weight mass cannot be computed",
               "without DOSE polygons (section 4.3). Headline claims are",
               "downgraded (checklist item 2).")
)
fs_write(join_checks, "qa", "join_checks.json")

# Assertions (section 8.3) --------------------------------------------------------------

prim <- panel %>% filter(primary)
key_unique <- !anyDuplicated(prim[c("GID_1", "year")])
prev_log <- panel$lgrp_pc_usd[match(paste(prim$GID_1, prim$year - 1L),
                                    paste(panel$GID_1, panel$year))]
consecutive <- all(abs(100 * (prim$lgrp_pc_usd - prev_log) - prim$g) < 1e-8)
logs_ok <- all(is.finite(prim$lgrp_pc_usd)) && all(is.finite(prev_log)) &&
  all(prim$grp_pc > 0)
no_breaks <- all(prim$struct_code == 0L)
set.seed(fs_seed(design$hash, "rate-assertion"))
check_rows <- prim[sample(nrow(prim), 500), c("GID_1", "year", "rate")]
rate_exact <- all(vapply(seq_len(nrow(check_rows)), function(i) {
  cl <- climate[climate$GID_1 == check_rows$GID_1[i] &
                  climate$year %in% (check_rows$year[i] - 4:0), ]
  nrow(cl) == 5 && !anyNA(cl$TM) &&
    abs(unname(coef(lm(TM ~ year, cl))[2]) - check_rows$rate[i]) < 1e-10
}, logical(1)))
schema_checks <- list(
  assertion_1_unique_keys = key_unique,
  assertion_2_consecutive_years = consecutive,
  assertion_3_logs_positive_finite = logs_ok,
  assertion_4_no_structural_break_pairs = no_breaks,
  assertion_5_fixed_usd_growth_max_gap_log_points = usd_gap,
  assertion_5_pass = usd_gap <= design$dose$fixed_usd_growth_tolerance,
  assertion_5_pairs_compared = usd_pairs,
  assertion_5_pairs_in_only_one_series = usd_only_one,
  assertion_6_crosswalk_weights = "not evaluable (string join)",
  assertion_7_rate_windows_exact_500_rows = rate_exact,
  assertion_8_climate_before_outcomes = TRUE,
  annual_completeness_rule = complete_share_required,
  dose_region_years_incomplete_months = n_incomplete_months,
  T0 = spec$T0, P_mean = spec$P_mean, P_sd = spec$P_sd,
  spline_knots = spec$spline_knots, spline_boundary = spec$spline_boundary,
  primary_n = nrow(prim), primary_regions = n_distinct(prim$GID_1),
  primary_countries = n_distinct(prim$GID_0),
  primary_years = paste(range(prim$year), collapse = "-"),
  regions_missing_zone = n_distinct(prim$GID_1[is.na(prim$zone)]),
  regions_missing_centroid = n_distinct(
    prim$GID_1[!prim$GID_1 %in% centroids$GID_1]),
  regions_missing_pop_base = n_distinct(prim$GID_1[is.na(prim$pop_base)])
)
fs_write(schema_checks, "qa", "schema_checks.json")
hard_fail <- !all(key_unique, consecutive, logs_ok, no_breaks, rate_exact,
                  schema_checks$assertion_5_pass)

# Outputs --------------------------------------------------------------------------------

fs_write(panel, "data", "panel_analysis.parquet")
# Stage 2 reads only these two files: climate and sample keys, no outcome.
fs_write(climate_terms, "data", "climate_terms.parquet")
fs_write(panel %>% select(GID_0, GID_1, year, primary, zone, continent,
                          decade, year_c, year_c2),
         "data", "sample_key.parquet")
fs_write(region_climate, "data", "region_climate.csv")
saveRDS(list(manifest_hash = design$manifest_hash),
        fs_path("data", "manifest_hash.rds"))

fs_stage_status(
  "1_panel",
  if (hard_fail) "STOP" else "PASS_WITH_WARNING",
  c(if (hard_fail) "a required assertion failed; see qa/schema_checks.json",
    "geography gate not met: DOSE-to-GADM 4.1 join is by identifier only",
    sprintf("%d of %d primary regions have an exact name match",
            join_checks$primary_exact_name_match, join_checks$primary_regions)),
  inputs = c(data_file("DOSE_V2.14.csv"), wcd_files, FS_DESIGN_FILE),
  outputs = fs_path("data", c("panel_analysis.parquet", "sample_flow.csv",
                              "geography_crosswalk.csv")),
  design = design
)
if (hard_fail) stop("Stage 1 assertion failure; see qa/schema_checks.json")
