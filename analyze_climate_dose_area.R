#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else
  "econometrics/data/data_climate_dose_area_1950_2023.csv"
output_dir <- if (length(args) >= 2L) args[[2L]] else
  "econometrics/output/climate_dose_area_analysis"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file, na.strings = c("", "NA", "NaN"))
required_id <- c("year", "GID_1")
if (!all(required_id %in% names(dat))) {
  stop("Input must contain columns: ", paste(required_id, collapse = ", "))
}
if (anyDuplicated(dat[, .(GID_1, year)])) {
  stop("GID_1-year keys are not unique.")
}

climate_vars <- setdiff(
  names(dat)[vapply(dat, is.numeric, logical(1L))],
  "year"
)
if (!length(climate_vars)) stop("No numeric climate variables found.")

# Non-finite values cannot enter means, standard deviations, or correlations.
for (v in climate_vars) {
  set(dat, which(!is.finite(dat[[v]])), v, NA_real_)
}

# Calculate each area's 1950-2023 long-run mean and sample standard deviation.
# CV is supplied as a complementary relative-variation measure; the standardized
# anomalies below are the requested deviations in units of area-specific SD.
stats_list <- vector("list", length(climate_vars))
names(stats_list) <- climate_vars
z_dat <- dat[, .(year, GID_1)]

for (v in climate_vars) {
  area_stats <- dat[, {
    x <- get(v)
    n_ok <- sum(!is.na(x))
    mu <- if (n_ok) mean(x, na.rm = TRUE) else NA_real_
    sig <- if (n_ok >= 2L) sd(x, na.rm = TRUE) else NA_real_
    list(n = n_ok, long_run_mean = mu, long_run_sd = sig)
  }, by = GID_1]
  area_stats[, variable := v]
  area_stats[, cv_abs_mean := fifelse(
    is.finite(long_run_mean) & long_run_mean != 0,
    long_run_sd / abs(long_run_mean),
    NA_real_
  )]
  setcolorder(
    area_stats,
    c("GID_1", "variable", "n", "long_run_mean", "long_run_sd", "cv_abs_mean")
  )
  stats_list[[v]] <- area_stats

  idx <- match(dat$GID_1, area_stats$GID_1)
  z <- (dat[[v]] - area_stats$long_run_mean[idx]) / area_stats$long_run_sd[idx]
  z[!is.finite(z)] <- NA_real_ # Includes areas with zero historical SD.
  set(z_dat, j = paste0(v, "_z"), value = z)
}

long_run_stats <- rbindlist(stats_list, use.names = TRUE)
setorder(long_run_stats, GID_1, variable)

variable_summary <- long_run_stats[, .(
  areas_with_data = sum(n > 0L),
  areas_with_positive_sd = sum(is.finite(long_run_sd) & long_run_sd > 0),
  median_area_mean = median(long_run_mean, na.rm = TRUE),
  median_area_sd = median(long_run_sd, na.rm = TRUE),
  median_area_cv = median(cv_abs_mean, na.rm = TRUE)
), by = variable]

pairwise_n <- function(x) {
  ok <- !is.na(x)
  crossprod(ok * 1L)
}

correlation_tests <- function(x, scale_name, method) {
  r <- suppressWarnings(cor(x, use = "pairwise.complete.obs", method = method))
  n_mat <- pairwise_n(x)
  pairs <- which(upper.tri(r), arr.ind = TRUE)
  out <- data.table(
    scale = scale_name,
    method = method,
    variable_1 = colnames(r)[pairs[, 1L]],
    variable_2 = colnames(r)[pairs[, 2L]],
    n = as.integer(n_mat[pairs]),
    estimate = r[pairs]
  )
  out[, df := n - 2L]
  out[, statistic := estimate * sqrt(df / pmax(1 - estimate^2, .Machine$double.eps))]
  out[, p_value := 2 * pt(abs(statistic), df = df, lower.tail = FALSE)]
  out[!is.finite(estimate) | n < 3L, `:=`(
    statistic = NA_real_, p_value = NA_real_
  )]
  out[, p_adjust_bh := p.adjust(p_value, method = "BH")]
  out[]
}

raw_x <- as.matrix(dat[, ..climate_vars])
z_vars <- paste0(climate_vars, "_z")
z_x <- as.matrix(z_dat[, ..z_vars])
colnames(z_x) <- climate_vars

cor_tests <- rbindlist(list(
  correlation_tests(raw_x, "raw", "pearson"),
  correlation_tests(raw_x, "raw", "spearman"),
  correlation_tests(z_x, "within_area_z", "pearson"),
  correlation_tests(z_x, "within_area_z", "spearman")
))

trend_tests <- function(x, scale_name, method, years) {
  ans <- rbindlist(lapply(seq_len(ncol(x)), function(j) {
    ok <- is.finite(x[, j]) & is.finite(years)
    n <- sum(ok)
    r <- if (n >= 3L) suppressWarnings(cor(years[ok], x[ok, j], method = method)) else NA_real_
    df <- n - 2L
    statistic <- if (is.finite(r))
      r * sqrt(df / max(1 - r^2, .Machine$double.eps)) else NA_real_
    p <- if (is.finite(statistic))
      2 * pt(abs(statistic), df = df, lower.tail = FALSE) else NA_real_
    data.table(
      scale = scale_name, method = method, variable = colnames(x)[j],
      n = n, estimate = r, df = df, statistic = statistic, p_value = p
    )
  }))
  ans[, p_adjust_bh := p.adjust(p_value, method = "BH")]
  ans[]
}

trend_tests_all <- rbindlist(list(
  trend_tests(raw_x, "raw", "pearson", dat$year),
  trend_tests(raw_x, "raw", "spearman", dat$year),
  trend_tests(z_x, "within_area_z", "pearson", dat$year),
  trend_tests(z_x, "within_area_z", "spearman", dat$year)
))

matrix_table <- function(x, method) {
  result <- as.data.table(
    suppressWarnings(cor(x, use = "pairwise.complete.obs", method = method)),
    keep.rownames = "variable"
  )
  result
}

fwrite(long_run_stats, file.path(output_dir, "long_run_stats_by_area.csv"))
fwrite(variable_summary, file.path(output_dir, "variable_summary.csv"))
fwrite(z_dat, file.path(output_dir, "standardized_anomalies_by_area_year.csv"))
fwrite(cor_tests, file.path(output_dir, "pairwise_correlation_tests.csv"))
fwrite(trend_tests_all, file.path(output_dir, "correlation_with_year_tests.csv"))
fwrite(matrix_table(raw_x, "pearson"), file.path(output_dir, "correlation_matrix_raw_pearson.csv"))
fwrite(matrix_table(z_x, "pearson"), file.path(output_dir, "correlation_matrix_within_area_z_pearson.csv"))

cat("Input:", input_file, "\n")
cat("Rows:", nrow(dat), " Areas:", uniqueN(dat$GID_1),
    " Years:", min(dat$year), "-", max(dat$year), "\n")
cat("Climate variables:", paste(climate_vars, collapse = ", "), "\n")
cat("Outputs written to:", output_dir, "\n")
