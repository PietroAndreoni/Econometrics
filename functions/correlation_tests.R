# Pairwise correlation tests across the columns of a numeric matrix, and tests
# of each column's correlation with time. Both use pairwise-complete
# observations, a t test on the correlation, and Benjamini-Hochberg adjusted
# p-values. With method = "spearman" the t test is the usual large-sample
# approximation.

.correlation_t_test <- function(estimate, n) {
  df <- n - 2L
  statistic <- estimate * sqrt(df / pmax(1 - estimate^2, .Machine$double.eps))
  p_value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  invalid <- !is.finite(estimate) | n < 3L
  statistic[invalid] <- NA_real_
  p_value[invalid] <- NA_real_
  list(df = df, statistic = statistic, p_value = p_value)
}

correlation_tests <- function(x, method = c("pearson", "spearman"),
                              label = NA_character_) {
  method <- match.arg(method)
  r <- suppressWarnings(
    stats::cor(x, use = "pairwise.complete.obs", method = method)
  )
  complete <- !is.na(x)
  n_mat <- crossprod(complete * 1L)
  pairs <- which(upper.tri(r), arr.ind = TRUE)
  n <- as.integer(n_mat[pairs])
  test <- .correlation_t_test(r[pairs], n)

  tibble::tibble(
    label = label,
    method = method,
    variable_1 = colnames(r)[pairs[, 1L]],
    variable_2 = colnames(r)[pairs[, 2L]],
    n = n,
    estimate = r[pairs],
    df = test$df,
    statistic = test$statistic,
    p_value = test$p_value,
    p_adjust_bh = stats::p.adjust(test$p_value, method = "BH")
  )
}

trend_tests <- function(x, years, method = c("pearson", "spearman"),
                        label = NA_character_) {
  method <- match.arg(method)
  x <- as.matrix(x)
  rows <- lapply(seq_len(ncol(x)), function(j) {
    ok <- is.finite(x[, j]) & is.finite(years)
    n <- sum(ok)
    r <- if (n >= 3L) {
      suppressWarnings(stats::cor(years[ok], x[ok, j], method = method))
    } else {
      NA_real_
    }
    c(n = n, estimate = r)
  })
  n <- as.integer(vapply(rows, `[[`, numeric(1), "n"))
  estimate <- vapply(rows, `[[`, numeric(1), "estimate")
  test <- .correlation_t_test(estimate, n)

  tibble::tibble(
    label = label,
    method = method,
    variable = colnames(x),
    n = n,
    estimate = estimate,
    df = test$df,
    statistic = test$statistic,
    p_value = test$p_value,
    p_adjust_bh = stats::p.adjust(test$p_value, method = "BH")
  )
}
