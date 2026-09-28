source("prepare_gadm2_national_climate_shocks.R")

# Monthly annualisation: equal-month temperature mean, precipitation sum, and a
# strict missing-month rule.
monthly_temperature <- rep(c(1:30, 31, 100), each = 12L)
monthly_precipitation <- monthly_temperature / 12
stopifnot(isTRUE(all.equal(
  annualize_monthly_values(monthly_temperature, "temperature"),
  c(1:30, 31, 100)
)))
stopifnot(isTRUE(all.equal(
  annualize_monthly_values(monthly_precipitation, "precipitation"),
  c(1:30, 31, 100)
)))
monthly_incomplete <- monthly_temperature
monthly_incomplete[5] <- NA_real_
stopifnot(is.na(annualize_monthly_values(
  monthly_incomplete,
  "temperature"
)[1]))

# Region A rises and region B falls.  The first available shock must use exactly
# 1950:1979 and exclude the 1980 realisation itself.
years <- 1950:1981
annual <- cbind(
  A = c(1:30, 31, 100),
  B = c(30:1, 0, -50)
)
moments <- previous_window_moments(annual, years, window = 30L)
stopifnot(all(is.na(moments$z[1:30, ])))
stopifnot(isTRUE(all.equal(unname(moments$mean[31, ]), c(15.5, 15.5))))
stopifnot(isTRUE(all.equal(unname(moments$sd[31, ]), rep(sqrt(77.5), 2))))
stopifnot(isTRUE(all.equal(unname(moments$z[31, ]^2), c(3.1, 3.1))))
stopifnot(isTRUE(all.equal(
  unname(moments$z[32, ]^2),
  c(89.9645161290323, 53.6806451612903),
  tolerance = 1e-12
)))

# Constant and incomplete histories are unavailable, never Inf/NaN.
constant <- previous_window_moments(
  matrix(c(rep(5, 30), 6), ncol = 1),
  1950:1980,
  window = 30L
)
stopifnot(is.na(constant$z[31, 1]))
with_missing <- annual
with_missing[10, 1] <- NA_real_
missing_moments <- previous_window_moments(with_missing, years, window = 30L)
stopifnot(is.na(missing_moments$z[31, 1]))
calendar_gap_error <- try(
  previous_window_moments(annual[-10, ], years[-10], window = 30L),
  silent = TRUE
)
stopifnot(inherits(calendar_gap_error, "try-error"))

# Fixed population shares A=.25 and B=.75.  Sign is classified before squaring.
population <- data.frame(
  GID_0 = c("AAA", "AAA"),
  GID_2 = c("AAA_A", "AAA_B"),
  GADM_GID_2 = c("AAA.A", "AAA.B"),
  population_2015 = c(25, 75),
  population_share_2015 = c(0.25, 0.75),
  stringsAsFactors = FALSE
)
national <- aggregate_national_z2(
  moments$z,
  years,
  population$GID_2,
  population,
  variable = "temperature",
  window = 30L
)
row_1980 <- national[national$year == 1980, ]
stopifnot(isTRUE(all.equal(row_1980$population_weighted_z2, 3.1)))
stopifnot(isTRUE(all.equal(
  row_1980$population_weighted_z2_positive,
  0.775
)))
stopifnot(isTRUE(all.equal(
  row_1980$population_weighted_z2_negative,
  2.325
)))
stopifnot(isTRUE(all.equal(
  row_1980$population_weighted_z2,
  row_1980$population_weighted_z2_positive +
    row_1980$population_weighted_z2_negative
)))
stopifnot(row_1980$complete_population_coverage)
stopifnot(isTRUE(all.equal(row_1980$population_coverage, 1)))

# If a positive-population region is unavailable, the statistic is normalized
# over observed population and the original coverage remains explicit.
partial_z <- moments$z
partial_z[31, 1] <- NA_real_
partial <- aggregate_national_z2(
  partial_z,
  years,
  population$GID_2,
  population,
  variable = "temperature",
  window = 30L
)
partial_1980 <- partial[partial$year == 1980, ]
stopifnot(isTRUE(all.equal(partial_1980$population_coverage, 0.75)))
stopifnot(isTRUE(all.equal(partial_1980$population_weighted_z2, 3.1)))
stopifnot(isTRUE(all.equal(
  partial_1980$population_weighted_z2_negative,
  3.1
)))
stopifnot(isTRUE(all.equal(
  partial_1980$population_weighted_z2_positive,
  0
)))
stopifnot(!partial_1980$complete_population_coverage)
validate_national_output(national)

cat("GADM2 climate-shock unit checks passed.\n")
