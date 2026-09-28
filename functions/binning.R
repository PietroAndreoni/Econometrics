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

# Equal-frequency signed bins centred on zero, coded -n_side..n_side with 0
# the reference bin, which always contains zero. The break probabilities are
# (0.5, 1.5, ..., n_side - 0.5) / (n_side + 0.5): the reference bin takes half
# a bin's share and every other bin a full share. Two ways to apply them:
#   "symmetric"  (default) quantiles of |x|, mirrored: breaks at +-q, so the
#                reference bin and each pair of bins are symmetric around zero
#                and each pair +-k holds the same number of observations, split
#                by sign.
#   "per_side"   quantiles of each side of zero separately. Every bin on a side
#                holds the same number of observations, so the bins are
#                equal-frequency even when one sign dominates (e.g. warm
#                anomalies under warming); the reference bin is then not
#                symmetric in value.
# `sample` selects the observations the quantiles are computed on (e.g. the
# estimation rows); every value of `x` is then assigned. The breaks are attached
# as attr "breaks".
signed_quantile_bin <- function(x, n_side = 5L, sample = rep(TRUE, length(x)),
                                method = c("symmetric", "per_side")) {
  method <- match.arg(method)
  probs <- (seq_len(n_side) - 0.5) / (n_side + 0.5)
  side <- function(v) {
    unname(stats::quantile(v[is.finite(v)], probs, names = FALSE))
  }
  breaks <- if (method == "per_side") {
    c(-rev(side(-x[sample & x < 0])), side(x[sample & x > 0]))
  } else {
    q <- side(abs(x[sample]))
    c(-rev(q), q)
  }
  codes <- findInterval(x, breaks) - n_side
  attr(codes, "breaks") <- breaks
  codes
}

# Equal-frequency bins: `n_bins` quantile bins of `x` over `sample`, as a
# factor labelled by interval. `lower` is the lowest break (0 for magnitudes,
# -Inf for levels); the top bin is open. Tied quantiles are merged, so a
# variable with a mass point can end up with fewer bins.
quantile_bin <- function(x, n_bins, sample = rep(TRUE, length(x)),
                         lower = -Inf) {
  interior <- stats::quantile(
    x[sample & is.finite(x)],
    probs = seq_len(n_bins - 1L) / n_bins,
    names = FALSE
  )
  breaks <- unique(c(lower, interior, Inf))
  cut(x, breaks = breaks, labels = interval_labels(breaks), right = FALSE,
      include.lowest = TRUE)
}

# "a-b" labels for consecutive breaks, "<b" and ">=a" for open ends.
interval_labels <- function(breaks, digits = 3L) {
  f <- function(v) format(signif(v, digits), trim = TRUE, scientific = FALSE,
                          drop0trailing = TRUE)
  k <- length(breaks)
  lo <- breaks[-k]
  hi <- breaks[-1]
  ifelse(is.infinite(lo), paste0("<", f(hi)),
         ifelse(is.infinite(hi), paste0(">=", f(lo)), paste0(f(lo), "-", f(hi))))
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
