# Out-of-sample fit of a fixed-effects model on its demeaned data.
#
# The outcome and regressors are first demeaned by the model's own fixed
# effects (including varying slopes such as GID_1[year]) over the full sample,
# then OLS is refit leaving out one fold at a time. Demeaning uses the held-out
# rows too, so this ranks specifications by how well the weather terms predict
# within-unit variation; it is not a true forecast test.

cv_rmse <- function(actual, pred) {
  keep <- !is.na(pred) & !is.na(actual)
  sqrt(mean((actual[keep] - pred[keep])^2))
}

# `folds` assigns every observation used by `fit` to a fold, e.g. its year:
#   cv_demeaned(m, folds = data$year[fixest::obs(m)])
# A logical vector gives a single split, with TRUE rows held out, e.g. the last
# year or the most extreme year of every region. Returns the RMSE per fold and
# the pooled RMSE (root mean of the fold MSEs).
cv_demeaned <- function(fit, folds) {
  demeaned <- fixest::demean(fit)
  if (length(folds) != nrow(demeaned)) {
    stop("`folds` must have one entry per observation used by the model.")
  }
  y <- demeaned[, 1]
  X <- demeaned[, -1, drop = FALSE]

  held_out <- if (is.logical(folds)) TRUE else sort(unique(folds))
  per_fold <- vapply(held_out, function(fold) {
    test <- if (is.logical(folds)) folds else folds == fold
    ols <- stats::lm.fit(X[!test, , drop = FALSE], y[!test])
    beta <- ols$coefficients
    beta[is.na(beta)] <- 0
    cv_rmse(y[test], as.vector(X[test, , drop = FALSE] %*% beta))
  }, numeric(1))

  list(
    per_fold = tibble::tibble(fold = held_out, rmse = per_fold),
    rmse = sqrt(mean(per_fold^2, na.rm = TRUE))
  )
}
