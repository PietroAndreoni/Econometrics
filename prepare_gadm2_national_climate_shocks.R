# Build national climate-shock indices from GADM2 ERA5 monthly data.
#
# The climate inputs come from the Weighted Climate Dataset (WCD).  Its
# population-density option weights ERA5 grid cells *within* each GADM2 unit;
# it does not provide the population total of each GADM2 unit.  National
# aggregation therefore uses separate 2015 GADM 4.1 population totals from the
# European Commission's GHS-WUP-DUC R2025A data.
#
# For region r and year t, the script computes
#
#   z[r,t] = (x[r,t] - mean(x[r,t-30], ..., x[r,t-1])) /
#            sd(x[r,t-30], ..., x[r,t-1])
#
# using R's sample standard deviation.  Temperature is the mean of twelve
# monthly values and precipitation is their sum.  A year or reference window
# containing any missing month/year is unavailable.  The national outputs are
# fixed-2015-population-weighted sums of z^2.  Positive and negative components
# are selected using the sign of z before it is squared, so total = positive +
# negative.  Where WCD has no usable climate history for a GADM2 unit, weights
# are renormalised over the covered population and `population_coverage` records
# the original national population share represented.  Zero-coverage results
# remain missing.

WCD_GADM2_WEIGHT_YEAR <- 2015L
SHOCK_WINDOW_YEARS <- 20L

GHSL_POPULATION_URL <- paste0(
  "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/",
  "GHS_WUP_DUC_GLOBE_R2025A/V1-0/",
  "GHS_WUP_DUC_GLOBE_R2025A_V1_0.zip"
)
GHSL_POPULATION_MEMBER <-
  "GHS_WUP_DUC_GLOBE_R2025A_V1_0_GADM41_2015_level2.csv"

required_packages <- c("arrow", "curl", "data.table", "jsonlite")

check_required_packages <- function(packages = required_packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Install missing package(s) before running this script: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

find_project_root <- function() {
  candidates <- c(".", "econometrics")
  is_root <- vapply(
    candidates,
    function(path) {
      file.exists(file.path(path, "download_weighted_climate_data.R")) &&
        dir.exists(file.path(path, "data"))
    },
    logical(1)
  )
  if (!any(is_root)) {
    stop(
      "Run this script from the Econometrics project root or its parent.",
      call. = FALSE
    )
  }
  normalizePath(candidates[which(is_root)[1]], winslash = "/", mustWork = TRUE)
}

parse_wcd_month <- function(x) {
  value <- sub("^X", "", as.character(x))
  if (anyNA(value) || any(!grepl("^[0-9]{6}$", value))) {
    stop("WCD monthly dates must have the form XYYYYMM.", call. = FALSE)
  }
  out <- as.Date(paste0(value, "01"), format = "%Y%m%d")
  if (anyNA(out)) {
    stop("At least one WCD monthly date could not be parsed.", call. = FALSE)
  }
  out
}

validate_complete_month_axis <- function(dates) {
  if (anyDuplicated(dates) || is.unsorted(dates, strictly = TRUE)) {
    stop("Monthly dates must be unique and strictly increasing.", call. = FALSE)
  }
  expected <- seq.Date(dates[1], dates[length(dates)], by = "month")
  if (!identical(as.character(dates), as.character(expected))) {
    stop("The WCD monthly time axis contains a gap.", call. = FALSE)
  }
  year <- as.integer(format(dates, "%Y"))
  month <- as.integer(format(dates, "%m"))
  month_by_year <- split(month, year)
  complete <- vapply(
    month_by_year,
    function(m) identical(as.integer(m), 1:12),
    logical(1)
  )
  if (!all(complete)) {
    stop(
      "Every retained calendar year must contain January through December.",
      call. = FALSE
    )
  }
  as.integer(names(month_by_year))
}

# Annualise one region after the common date axis has been validated.  Missing
# values invalidate the whole region-year instead of changing its seasonality.
annualize_monthly_values <- function(x, variable = c("temperature", "precipitation")) {
  variable <- match.arg(variable)
  if (length(x) %% 12L != 0L) {
    stop("Monthly values must contain complete 12-month blocks.", call. = FALSE)
  }
  monthly <- matrix(as.numeric(x), nrow = 12L)
  complete <- colSums(is.finite(monthly)) == 12L
  annual <- if (variable == "temperature") {
    colMeans(monthly)
  } else {
    colSums(monthly)
  }
  annual[!complete] <- NA_real_
  annual
}

# Read one 47k-column parquet at a time.  Keeping the monthly table only until
# its much smaller annual matrix has been constructed avoids a 43-million-row
# long pivot and keeps peak memory bounded.
read_wcd_annual <- function(path, variable = c("temperature", "precipitation")) {
  variable <- match.arg(variable)
  cache_path <- sub("_monthly\\.parquet$", "_yearly.parquet", path)
  if (file.exists(cache_path) && file.mtime(cache_path) >= file.mtime(path)) {
    message("Reading cached annual data from ", basename(cache_path), " ...")
    cached <- as.data.frame(arrow::read_parquet(cache_path))
    if (!"year" %in% names(cached)) {
      stop("Annual cache has no `year` column: ", cache_path, call. = FALSE)
    }
    years <- as.integer(cached$year)
    cached$year <- NULL
    annual <- as.matrix(cached)
    storage.mode(annual) <- "double"
    rownames(annual) <- as.character(years)
    return(list(year = years, GID_2 = colnames(annual), value = annual))
  }

  message("Reading and annualising ", basename(path), " ...")
  monthly <- as.data.frame(arrow::read_parquet(path))
  date_col <- if ("Date" %in% names(monthly)) "Date" else names(monthly)[1]
  dates <- parse_wcd_month(monthly[[date_col]])
  years <- validate_complete_month_axis(dates)
  monthly[[date_col]] <- NULL
  keep <- !is.na(names(monthly)) & !(names(monthly) %in% c("", "?", "NA"))
  monthly <- monthly[, keep, drop = FALSE]
  if (!ncol(monthly)) {
    stop("No GADM2 columns found in ", basename(path), ".", call. = FALSE)
  }

  annual <- vapply(
    monthly,
    annualize_monthly_values,
    numeric(length(years)),
    variable = variable,
    USE.NAMES = TRUE
  )
  if (is.null(dim(annual))) {
    annual <- matrix(annual, ncol = 1L)
    colnames(annual) <- names(monthly)
  }
  rownames(annual) <- as.character(years)
  storage.mode(annual) <- "double"
  rm(monthly)
  invisible(gc(verbose = FALSE))

  annual_cache <- data.frame(year = years, annual, check.names = FALSE)
  arrow::write_parquet(annual_cache, cache_path)
  rm(annual_cache)
  message("Wrote annual cache ", basename(cache_path))

  list(year = years, GID_2 = colnames(annual), value = annual)
}

# Return the previous-window mean, sample SD, and signed standardised anomaly.
# A rolling set of sums makes this linear in the number of region-years.
previous_window_moments <- function(x, years, window = SHOCK_WINDOW_YEARS) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  years <- as.integer(years)
  window <- as.integer(window)

  if (nrow(x) != length(years)) {
    stop("`x` must have one row per year.", call. = FALSE)
  }
  if (length(years) > 1L && any(diff(years) != 1L)) {
    stop("Rolling moments require consecutive calendar years.", call. = FALSE)
  }
  if (window < 2L) {
    stop("`window` must be at least two years.", call. = FALSE)
  }

  previous_mean <- previous_sd <- z <- matrix(
    NA_real_, nrow = nrow(x), ncol = ncol(x),
    dimnames = dimnames(x)
  )
  if (nrow(x) <= window) {
    return(list(mean = previous_mean, sd = previous_sd, z = z))
  }

  initial <- x[seq_len(window), , drop = FALSE]
  initial_ok <- is.finite(initial)
  running_n <- colSums(initial_ok)
  running_sum <- colSums(ifelse(initial_ok, initial, 0))
  running_sumsq <- colSums(ifelse(initial_ok, initial^2, 0))

  for (i in seq.int(window + 1L, nrow(x))) {
    good_history <- running_n == window
    mu <- running_sum / window
    variance <- (running_sumsq - running_sum^2 / window) / (window - 1L)
    # Round-off can produce tiny negative values for an exactly constant series.
    variance[variance < 0 & variance > -1e-12] <- 0
    sigma <- rep(NA_real_, length(variance))
    variance_ok <- good_history & is.finite(variance) & variance >= 0
    sigma[variance_ok] <- sqrt(variance[variance_ok])
    good <- good_history & is.finite(x[i, ]) & is.finite(sigma) & sigma > 0

    previous_mean[i, good] <- mu[good]
    previous_sd[i, good] <- sigma[good]
    z[i, good] <- (x[i, good] - mu[good]) / sigma[good]

    outgoing <- x[i - window, ]
    incoming <- x[i, ]
    outgoing_ok <- is.finite(outgoing)
    incoming_ok <- is.finite(incoming)
    running_n <- running_n - outgoing_ok + incoming_ok
    running_sum <- running_sum - ifelse(outgoing_ok, outgoing, 0) +
      ifelse(incoming_ok, incoming, 0)
    running_sumsq <- running_sumsq - ifelse(outgoing_ok, outgoing^2, 0) +
      ifelse(incoming_ok, incoming^2, 0)
  }

  list(mean = previous_mean, sd = previous_sd, z = z)
}

validate_population_table <- function(population, expected_ids = NULL) {
  needed <- c(
    "GID_0", "GID_2", "GADM_GID_2", "population_2015",
    "population_share_2015"
  )
  missing <- setdiff(needed, names(population))
  if (length(missing)) {
    stop(
      "Population table lacks column(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(population$GID_2)) {
    stop("Population GID_2 values are not unique.", call. = FALSE)
  }
  if (any(!is.finite(population$population_2015)) ||
      any(population$population_2015 < 0)) {
    stop("Population values must be finite and nonnegative.", call. = FALSE)
  }
  if (!is.null(expected_ids) && !setequal(expected_ids, population$GID_2)) {
    missing_pop <- setdiff(expected_ids, population$GID_2)
    extra_pop <- setdiff(population$GID_2, expected_ids)
    stop(
      "Population/climate GADM2 identifiers differ. Missing population IDs: ",
      length(missing_pop), "; extra population IDs: ", length(extra_pop), ".",
      call. = FALSE
    )
  }
  share_sum <- tapply(
    population$population_share_2015,
    population$GID_0,
    sum
  )
  if (any(!is.finite(share_sum)) || any(abs(share_sum - 1) > 1e-10)) {
    stop("Population shares do not sum to one within country.", call. = FALSE)
  }
  invisible(TRUE)
}

read_or_download_population <- function(data_dir, expected_ids = NULL) {
  derived_path <- file.path(data_dir, "gadm2_population_2015_ghsl.parquet")
  if (file.exists(derived_path)) {
    population <- as.data.frame(arrow::read_parquet(derived_path))
    validate_population_table(population, expected_ids)
    return(population)
  }

  local_csv <- file.path(data_dir, GHSL_POPULATION_MEMBER)
  cleanup_dir <- NULL
  if (!file.exists(local_csv)) {
    local_archive <- file.path(
      data_dir,
      "GHS_WUP_DUC_GLOBE_R2025A_V1_0.zip"
    )
    if (!file.exists(local_archive)) {
      local_archive <- tempfile(fileext = ".zip")
      message("Downloading GHS-WUP-DUC 2015 GADM2 population totals ...")
      curl::curl_download(GHSL_POPULATION_URL, local_archive, quiet = FALSE)
    }
    cleanup_dir <- tempfile("ghsl_population_")
    dir.create(cleanup_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(
      local_archive,
      files = GHSL_POPULATION_MEMBER,
      exdir = cleanup_dir
    )
    local_csv <- file.path(cleanup_dir, GHSL_POPULATION_MEMBER)
    on.exit(unlink(cleanup_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  raw <- data.table::fread(
    local_csv,
    select = c("GID_2", "Tot_Pop"),
    showProgress = FALSE
  )
  raw <- raw[raw$GID_2 != "?", , drop = FALSE]
  population <- data.frame(
    GID_0 = substr(raw$GID_2, 1L, 3L),
    GID_2 = gsub("\\.", "_", raw$GID_2),
    GADM_GID_2 = raw$GID_2,
    population_2015 = as.numeric(raw$Tot_Pop),
    stringsAsFactors = FALSE
  )
  country_population <- ave(
    population$population_2015,
    population$GID_0,
    FUN = sum
  )
  if (any(!is.finite(country_population)) || any(country_population <= 0)) {
    stop("Every country must have positive total population.", call. = FALSE)
  }
  population$population_share_2015 <-
    population$population_2015 / country_population
  population <- population[order(population$GID_0, population$GID_2), ]
  rownames(population) <- NULL
  validate_population_table(population, expected_ids)
  arrow::write_parquet(population, derived_path)
  message("Wrote ", normalizePath(derived_path, winslash = "/"))
  population
}

download_gadm2_climate <- function(project_root) {
  data_dir <- file.path(project_root, "data")
  wcd_env <- new.env(parent = globalenv())
  sys.source(
    file.path(project_root, "download_weighted_climate_data.R"),
    envir = wcd_env
  )

  selections <- list(
    temperature = list(
      variable = "avg. temperature",
      expected = "gadm2_era_tmp_pop_2015_monthly.parquet"
    ),
    precipitation = list(
      variable = "precipitation",
      expected = "gadm2_era_pre_pop_2015_monthly.parquet"
    )
  )

  paths <- lapply(selections, function(selection) {
    cached <- wcd_env$wcd_download(
      variable = selection$variable,
      source = "ERA5",
      geo_resolution = "gadm2",
      weight = "population density",
      weight_year = WCD_GADM2_WEIGHT_YEAR,
      time_frequency = "monthly",
      verbose = TRUE
    )
    if (!identical(basename(cached), selection$expected)) {
      stop(
        "WCD resolved ", basename(cached), " instead of ",
        selection$expected, ".",
        call. = FALSE
      )
    }
    local <- file.path(data_dir, selection$expected)
    if (!file.exists(local) || file.size(local) != file.size(cached)) {
      copied <- file.copy(cached, local, overwrite = TRUE)
      if (!copied) {
        stop("Could not copy ", selection$expected, " into data/.", call. = FALSE)
      }
    }
    normalizePath(local, winslash = "/", mustWork = TRUE)
  })
  names(paths) <- names(selections)
  paths
}

aggregate_national_z2 <- function(
    z,
    years,
    region_ids,
    population,
    variable,
    window = SHOCK_WINDOW_YEARS) {
  z <- as.matrix(z)
  years <- as.integer(years)
  if (nrow(z) != length(years) || ncol(z) != length(region_ids)) {
    stop("Shock matrix dimensions do not match years/regions.", call. = FALSE)
  }
  match_index <- match(region_ids, population$GID_2)
  if (anyNA(match_index)) {
    stop("At least one climate region has no population match.", call. = FALSE)
  }
  pop <- population[match_index, , drop = FALSE]
  if (!identical(as.character(pop$GID_2), as.character(region_ids))) {
    stop("Population alignment failed.", call. = FALSE)
  }

  report_rows <- seq.int(as.integer(window) + 1L, nrow(z))
  z <- z[report_rows, , drop = FALSE]
  report_years <- years[report_rows]
  country_regions <- split(seq_along(region_ids), pop$GID_0)

  result <- lapply(names(country_regions), function(country) {
    idx <- country_regions[[country]]
    region_population <- pop$population_2015[idx]
    positive_population <- region_population > 0
    weights <- region_population / sum(region_population)
    z_country <- z[, idx, drop = FALSE]
    valid <- is.finite(z_country)
    valid_positive_population <- valid[, positive_population, drop = FALSE]
    n_regions_with_shock <- rowSums(valid_positive_population)
    n_regions_positive_population <- sum(positive_population)
    complete <- n_regions_with_shock == n_regions_positive_population
    population_coverage <- as.numeric(valid %*% weights)

    q <- z_country^2
    q[!valid] <- 0
    total <- as.numeric(q %*% weights)
    # `valid & ...` turns the sign flag FALSE for missing shocks.  This matters
    # for zero-population regions: in matrix multiplication, NA * 0 is still NA.
    positive_indicator <- valid & z_country > 0
    negative_indicator <- valid & z_country < 0
    positive <- as.numeric((q * positive_indicator) %*% weights)
    negative <- as.numeric((q * negative_indicator) %*% weights)
    has_coverage <- is.finite(population_coverage) & population_coverage > 0
    # Conditional population weights sum to one over regions with usable shocks.
    # Coverage is retained so downstream analyses can impose their own threshold.
    total[has_coverage] <- total[has_coverage] / population_coverage[has_coverage]
    positive[has_coverage] <-
      positive[has_coverage] / population_coverage[has_coverage]
    negative[has_coverage] <-
      negative[has_coverage] / population_coverage[has_coverage]
    total[!has_coverage] <- NA_real_
    positive[!has_coverage] <- NA_real_
    negative[!has_coverage] <- NA_real_

    data.frame(
      GID_0 = country,
      year = report_years,
      variable = variable,
      population_weighted_z2 = total,
      population_weighted_z2_positive = positive,
      population_weighted_z2_negative = negative,
      population_coverage = population_coverage,
      complete_population_coverage = complete,
      population_coverage_ge_0_95 = population_coverage >= 0.95,
      n_regions_with_shock = n_regions_with_shock,
      n_regions_positive_population = n_regions_positive_population,
      n_regions_total = length(idx),
      population_2015 = sum(region_population),
      window_start = report_years - as.integer(window),
      window_end = report_years - 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

validate_national_output <- function(data) {
  key <- paste(data$GID_0, data$year, data$variable, sep = "|")
  if (anyDuplicated(key)) {
    stop("National output keys are not unique.", call. = FALSE)
  }
  value_cols <- c(
    "population_weighted_z2",
    "population_weighted_z2_positive",
    "population_weighted_z2_negative"
  )
  for (column in value_cols) {
    value <- data[[column]]
    if (any(value[is.finite(value)] < -1e-12)) {
      stop(column, " contains a negative value.", call. = FALSE)
    }
  }
  available <- stats::complete.cases(data[, value_cols])
  decomposition_error <- abs(
    data$population_weighted_z2 -
      data$population_weighted_z2_positive -
      data$population_weighted_z2_negative
  )
  if (any(decomposition_error[available] > 1e-9)) {
    stop("Total shocks do not equal positive plus negative shocks.", call. = FALSE)
  }
  if (any(data$population_coverage < -1e-12) ||
      any(data$population_coverage > 1 + 1e-12)) {
    stop("Population coverage lies outside [0, 1].", call. = FALSE)
  }
  if (any(data$complete_population_coverage &
          abs(data$population_coverage - 1) > 1e-10)) {
    stop("Complete rows do not have full population coverage.", call. = FALSE)
  }
  invisible(TRUE)
}

write_outputs <- function(national, population, project_root) {
  data_dir <- file.path(project_root, "data")
  result_dir <- file.path(project_root, "results", "gadm2_national_climate_shocks")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  national_path <- file.path(
    data_dir,
    "national_gadm2_era5_population_weighted_30y_shocks.parquet"
  )
  csv_path <- file.path(result_dir, "national_climate_shocks.csv")
  summary_path <- file.path(result_dir, "validation_summary.csv")

  arrow::write_parquet(national, national_path)
  utils::write.csv(national, csv_path, row.names = FALSE, na = "")

  summary <- do.call(
    rbind,
    lapply(split(national, national$variable), function(x) {
      finite <- is.finite(x$population_weighted_z2)
      data.frame(
        variable = x$variable[1],
        first_year = min(x$year),
        last_year = max(x$year),
        countries = length(unique(x$GID_0)),
        country_years = nrow(x),
        available_country_years = sum(finite),
        full_coverage_country_years = sum(x$complete_population_coverage),
        partial_coverage_country_years = sum(
          finite & !x$complete_population_coverage
        ),
        zero_coverage_country_years = sum(!finite),
        coverage_ge_0_95_country_years = sum(
          finite & x$population_coverage_ge_0_95
        ),
        min_population_coverage = min(x$population_coverage),
        max_decomposition_error = max(
          abs(
            x$population_weighted_z2 -
              x$population_weighted_z2_positive -
              x$population_weighted_z2_negative
          ),
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  utils::write.csv(summary, summary_path, row.names = FALSE)

  paths <- c(
    national_parquet = national_path,
    national_csv = csv_path,
    validation_summary = summary_path,
    population_weights = file.path(data_dir, "gadm2_population_2015_ghsl.parquet")
  )
  message(
    "Wrote:\n",
    paste(normalizePath(paths, winslash = "/", mustWork = TRUE), collapse = "\n")
  )
  invisible(paths)
}

main <- function() {
  check_required_packages()
  project_root <- find_project_root()
  climate_paths <- download_gadm2_climate(project_root)

  temperature <- read_wcd_annual(climate_paths$temperature, "temperature")
  precipitation <- read_wcd_annual(climate_paths$precipitation, "precipitation")
  if (!identical(temperature$year, precipitation$year)) {
    stop("Temperature and precipitation cover different years.", call. = FALSE)
  }
  if (!identical(temperature$GID_2, precipitation$GID_2)) {
    stop("Temperature and precipitation contain different GADM2 IDs/order.", call. = FALSE)
  }

  population <- read_or_download_population(
    file.path(project_root, "data"),
    expected_ids = temperature$GID_2
  )

  message("Computing strictly preceding 30-year regional moments ...")
  temperature_moments <- previous_window_moments(
    temperature$value,
    temperature$year
  )
  precipitation_moments <- previous_window_moments(
    precipitation$value,
    precipitation$year
  )

  national <- rbind(
    aggregate_national_z2(
      temperature_moments$z,
      temperature$year,
      temperature$GID_2,
      population,
      variable = "temperature"
    ),
    aggregate_national_z2(
      precipitation_moments$z,
      precipitation$year,
      precipitation$GID_2,
      population,
      variable = "precipitation"
    )
  )
  national$climate_source <- "ERA5"
  national$gadm_level <- "gadm2"
  national$wcd_spatial_weight <- "population density"
  national$wcd_weight_year <- WCD_GADM2_WEIGHT_YEAR
  national$population_source <- "GHS-WUP-DUC R2025A (GADM 4.1)"
  national$population_year <- 2015L
  national$window_years <- SHOCK_WINDOW_YEARS
  national <- national[
    order(national$GID_0, national$year, national$variable),
  ]
  rownames(national) <- NULL
  validate_national_output(national)
  write_outputs(national, population, project_root)
  invisible(national)
}

if (sys.nframe() == 0L) {
  main()
}
