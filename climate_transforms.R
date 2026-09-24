# Shared climate-series transforms.
#
# These functions are used both by `prepare_climate_data.R`, which builds the
# regional TM/RR panels, and by `prepare_global_climate_data.R`, which builds the
# single global temperature series. Keeping one copy guarantees the global series
# is filtered and standardised by exactly the same operators as the local panels.

# Hamilton (2018) regression filter. For a single series y ordered by year, the
# trend is the fitted value of
#   y_t = b0 + b1*y_{t-h} + b2*y_{t-h-1} + ... + bp*y_{t-h-p+1} + v_t
# estimated by OLS; the residual v_t is the cyclical part. The horizon h is set
# to the same window used for the rolling means.
#
# `train` is an optional logical vector (same length as y) selecting the rows on
# which the OLS coefficients are estimated; the lagged predictors are always
# built from the full series, so training rows can reach back before the training
# window for their lags. When `train` is NULL the whole series is used and the
# trend is returned for every row with available predictors; otherwise the trend
# is returned only for the training rows (the trend "over" that period). Returns
# a same-length vector, NA where predictors are unavailable or too few rows are
# available to estimate the regression.
hamilton_trend <- function(y, h, p = 4L, train = NULL) {
  h <- as.integer(h)
  p <- as.integer(p)
  n <- length(y)
  trend <- rep(NA_real_, n)
  if (is.null(train)) {
    train <- rep(TRUE, n)
  }
  train <- train & !is.na(train)
  if (n <= h + p) {
    return(trend)
  }
  lag_cols <- lapply(seq.int(0L, p - 1L), function(j) dplyr::lag(y, h + j))
  X <- as.data.frame(stats::setNames(lag_cols, paste0("lag", seq.int(h, h + p - 1L))))
  df <- cbind(data.frame(y = y), X)
  has_predictors <- stats::complete.cases(X)
  fit_rows <- train & stats::complete.cases(df)
  if (sum(fit_rows) <= p + 1L) {
    return(trend)
  }
  fit <- stats::lm(y ~ ., data = df[fit_rows, , drop = FALSE])
  out_rows <- train & has_predictors
  trend[out_rows] <- as.numeric(
    stats::predict(fit, newdata = df[out_rows, , drop = FALSE])
  )
  trend
}

# Standardised anomaly, using the reference moments of the *preceding* year so the
# current realisation never enters its own baseline. This is the operator behind
# `zTM` in test_functions.Rmd:
#   zTM = (TM - mean_TM_lag_20) / sd_TM_lag_20
# where mean_TM_20 / sd_TM_20 are right-aligned rolling moments.
standardized_anomaly <- function(x, rolling_mean, rolling_sd) {
  previous_mean <- dplyr::lag(rolling_mean)
  previous_sd <- dplyr::lag(rolling_sd)
  dplyr::if_else(
    is.finite(previous_sd) & previous_sd > 0,
    (x - previous_mean) / previous_sd,
    NA_real_
  )
}
