# Monte Carlo for the handout (handout_preferred_specification.Rmd, Section 13):
# the estimate of xi and the robust test when the temperature response changes
# sign with income (the linear income interaction of M2b: rich countries gain
# from warm years), and, as a check of the test's size, when the scaled model
# with xi = 0.6 is true. Saved to
# results/handout_preferred_specification/xi_montecarlo.rds.
#
# Samples: the fitted values of the data-generating model plus its residuals
# times one Rademacher sign per country (wild cluster), on the handout's rows;
# in each, xi by profile NLS and the robust test of xi = 0, 0.6, 1.
#
#   Rscript handout_xi_montecarlo.R

# The handout's data and estimators: handout_preferred_specification.R up to
# the profile over the fine grid.
lines <- readLines("handout_preferred_specification.R")
eval(parse(text = lines[seq_len(grep("^XI_FINE <- ", lines) - 1)]))
setFixest_nthreads(2)

ar_p_data <- function(xi0, data) {
  m <- fit_terms(c(scaled_terms, ar_terms), data = scaled_data(xi0, data))
  unname(zero_test(m, names(ar_terms))[["p"]])
}
simulate <- function(m, R, seed) {
  set.seed(seed)
  stopifnot(length(fitted(m)) == nrow(d))
  f <- fitted(m)
  e <- resid(m)
  countries <- unique(d$GID_0)
  t(vapply(seq_len(R), function(r) {
    s <- sample(c(-1, 1), length(countries), replace = TRUE)
    x <- d
    x$g <- f + e * s[match(d$GID_0, countries)]
    c(xi = xi_nls(x), p0 = ar_p_data(0, x), p06 = ar_p_data(0.6, x), p1 = ar_p_data(1, x))
  }, numeric(4)))
}
R <- 200
dgp <- list(M2b = models$M2b, scaled_0.6 = fit_scaled(0.6))
# Saved after each data-generating model; a rerun skips those done.
mc_file <- file.path("results", "handout_preferred_specification", "xi_montecarlo.rds")
res <- if (file.exists(mc_file)) readRDS(mc_file) else list()
for (k in setdiff(names(dgp), names(res))) {
  message("DGP ", k)
  res[[k]] <- simulate(dgp[[k]], R, 7)
  saveRDS(res, mc_file)
}
for (k in names(res)) {
  x <- res[[k]]
  cat(sprintf("DGP %s: xi hat median %.2f, IQR [%.2f, %.2f]; robust test rejects xi = 0 / 0.6 / 1 in %.0f / %.0f / %.0f%%\n",
              k, median(x[, "xi"]), quantile(x[, "xi"], .25), quantile(x[, "xi"], .75),
              100 * mean(x[, "p0"] < .05), 100 * mean(x[, "p06"] < .05), 100 * mean(x[, "p1"] < .05)))
}
