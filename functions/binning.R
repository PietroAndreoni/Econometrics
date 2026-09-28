# Binning helpers for weather anomalies and other regressors.

# Signed anomaly bins coded by their midpoints. Intervals are left-closed,
# [b_i, b_{i+1}); the interval containing zero is the reference and is coded 0,
# and the two open tails are coded half a bin width beyond the outer breaks. With
# the default breaks, z in [-0.5, 0.5) -> 0, [0.5, 1) -> 0.75, ..., >= 2.5 -> 2.75.
signed_anomaly_bin <- function(
    z,
    breaks = c(-2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5)
) {
  k <- length(breaks)
  codes <- c(
    breaks[[1]] - (breaks[[2]] - breaks[[1]]) / 2,
    (breaks[-k] + breaks[-1]) / 2,
    breaks[[k]] + (breaks[[k]] - breaks[[k - 1]]) / 2
  )
  codes[findInterval(z, breaks) + 1L]
}

# Bin the strictly positive values of a non-negative magnitude into `n_bins`
# quantile bins, with exact zeros kept as their own "0" reference level.
bin_nonzero_magnitude <- function(x, n_bins = 5L) {
  if (any(x < 0, na.rm = TRUE)) {
    stop("Binned magnitudes must be non-negative.")
  }

  positive <- !is.na(x) & x > 0
  result <- rep(NA_character_, length(x))
  result[!is.na(x) & x == 0] <- "0"

  if (!any(positive)) {
    return(factor(result, levels = "0"))
  }

  positive_bins <- fixest::bin(
    x[positive],
    paste0("cut::", n_bins)
  )
  result[positive] <- as.character(positive_bins)

  factor(result, levels = c("0", levels(positive_bins)))
}

# Equal-width bins leave very thin tails, whose coefficients carry almost no
# information yet still enter joint Wald tests. Merge each tail inward until it
# holds at least `min_units` distinct units; interior breaks are unchanged.
pool_tails <- function(x, breaks, unit, min_units) {
  index <- as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
  k <- length(breaks) - 1L
  lo <- 1L
  while (lo < k && dplyr::n_distinct(unit[index <= lo]) < min_units) lo <- lo + 1L
  hi <- k
  while (hi > lo && dplyr::n_distinct(unit[index >= hi]) < min_units) hi <- hi - 1L
  c(-Inf, breaks[(lo + 1L):hi], Inf)
}
