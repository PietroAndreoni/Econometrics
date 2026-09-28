# Coefficient tables, linear combinations and resampling inference for fitted
# fixest models.

# Coefficient table with confidence intervals. `vcov` is passed straight to
# fixest (e.g. ~GID_0, or a precomputed matrix such as a Conley VCOV). When
# `data` is supplied, the regions, countries and years actually used by the fit
# are counted from fixest::obs(), the only reliable source once NA removal and
# singleton dropping have run.
tidy_fixest <- function(fit, label = NULL, data = NULL, conf_level = 0.95,
                        vcov = NULL) {
  coefficients <- as.data.frame(fixest::coeftable(fit, vcov = vcov))
  names(coefficients) <- c("estimate", "std_error", "statistic", "p_value")
  intervals <- as.data.frame(
    stats::confint(fit, vcov = vcov, level = conf_level)
  )[rownames(coefficients), , drop = FALSE]

  out <- tibble::tibble(
    model = if (is.null(label)) NA_character_ else label,
    term = rownames(coefficients),
    estimate = coefficients$estimate,
    std_error = coefficients$std_error,
    statistic = coefficients$statistic,
    p_value = coefficients$p_value,
    conf_low = intervals[[1L]],
    conf_high = intervals[[2L]],
    observations = as.integer(stats::nobs(fit))
  )

  if (!is.null(data)) {
    used <- fixest::obs(fit)
    out$n_regions <- dplyr::n_distinct(data$GID_1[used])
    out$n_countries <- dplyr::n_distinct(data$GID_0[used])
    out$n_years <- dplyr::n_distinct(data$year[used])
  }
  out
}

# Estimates and delta-method standard errors of the linear combinations G %*% b,
# one per row of G. G is a matrix (or a named vector for a single combination)
# whose column names are coefficient names; coefficients not named in G get
# weight zero. Covariances across terms are retained. `df = Inf` gives normal
# intervals; pass e.g. the number of clusters minus one for t intervals.
linear_combination <- function(model, G, conf_level = 0.95, df = Inf) {
  if (is.null(dim(G))) {
    G <- matrix(G, nrow = 1L, dimnames = list(NULL, names(G)))
  }
  terms <- colnames(G)
  missing_terms <- setdiff(terms, names(stats::coef(model)))
  if (length(missing_terms)) {
    stop(
      "Coefficients missing from the model: ",
      paste(missing_terms, collapse = ", ")
    )
  }
  b <- stats::coef(model)[terms]
  V <- stats::vcov(model)[terms, terms, drop = FALSE]
  critical <- stats::qt(1 - (1 - conf_level) / 2, df)

  estimate <- as.numeric(G %*% b)
  std_error <- sqrt(pmax(rowSums((G %*% V) * G), 0))
  tibble::tibble(
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error
  )
}

# Single linear hypothesis w'b = 0 with its test statistic, e.g.
#   lincom(m, c(TM_dev_pos = 1, TM_dev_neg = -1), "TM symmetry")
lincom <- function(model, weights, label = NA_character_, conf_level = 0.95,
                   df = Inf) {
  out <- linear_combination(model, weights, conf_level = conf_level, df = df)
  out$statistic <- out$estimate / out$std_error
  out$p_value <- 2 * stats::pt(-abs(out$statistic), df)
  tibble::add_column(out, test = label, .before = 1)
}

# Joint Wald test on the coefficients matching `keep`, as a one-row table.
wald_row <- function(model, keep, vcov = NULL, label = NA_character_) {
  w <- fixest::wald(model, keep = keep, vcov = vcov, print = FALSE)
  tibble::tibble(
    test = label, keep = keep,
    F_stat = w$stat, df1 = w$df1, df2 = w$df2, p_value = w$p
  )
}

stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

# Empirical-Bayes shrinkage of unit-level estimates toward their
# precision-weighted mean (Morris 1983). `use` selects the estimates that
# inform the prior.
eb_shrink <- function(estimate, std_error, use = is.finite(estimate) &
                        is.finite(std_error)) {
  prior_mean <- stats::weighted.mean(estimate[use], w = 1 / std_error[use]^2)
  tau2 <- max(0, stats::var(estimate[use]) - mean(std_error[use]^2))
  weight <- tau2 / (tau2 + std_error^2)
  tibble::tibble(
    shrunk = weight * estimate + (1 - weight) * prior_mean,
    shrinkage_weight = weight,
    prior_mean = prior_mean,
    tau2 = tau2
  )
}

# Year block bootstrap of one coefficient. Whole years are resampled with
# replacement and relabelled onto the original time axis, so each drawn block
# occupies a distinct slot. Useful when a regressor varies only over time (e.g.
# a global climate series) and the effective sample is the number of years.
bootstrap_year_blocks <- function(formula, data, term, replications,
                                  year = "year", ...) {
  years <- sort(unique(data[[year]]))
  rows_by_year <- split(seq_len(nrow(data)), data[[year]])
  rows_by_year <- rows_by_year[as.character(years)]

  estimates <- vapply(seq_len(replications), function(b) {
    drawn <- sample(seq_along(years), length(years), replace = TRUE)
    index <- unlist(rows_by_year[drawn], use.names = FALSE)
    resampled <- data[index, , drop = FALSE]
    resampled[[year]] <- rep(years, lengths(rows_by_year[drawn]))
    fit <- try(fixest::feols(formula, data = resampled, ...), silent = TRUE)
    if (inherits(fit, "try-error") || !term %in% names(stats::coef(fit))) {
      return(NA_real_)
    }
    unname(stats::coef(fit)[term])
  }, numeric(1))
  estimates[is.finite(estimates)]
}
