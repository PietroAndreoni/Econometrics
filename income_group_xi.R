# The income scaling of the T response, xi, in two stages: the preferred
# specification's T arms (handout_preferred_specification.R) estimated
# separately by income group, then a power law in income fitted to the group
# coefficients by minimum distance. A diagnostic of the one-stage NLS xi
# (estimates.rds$xi_nls): the unrestricted group responses show whether the
# income gradient looks like a power law, flattens or changes sign. Saved to
# results/income_group_xi/<GDP>_<grouping>_J<groups>/estimates.rds and .png.
#
# First stage: one pooled regression (common year effects, country effects and
# quadratic trends, controls) with the T arms and their climate interaction
# times group indicators; P arms with climate and share as in the preferred
# specification (common to the groups); the levels of the trailing climate, of
# the trailing income and the group indicators. Groups: quantiles of the
# trailing income (5 years before, so predetermined) over the rows, a country
# moving group as it grows (GROUPS_BY=trailing, default), or of the country's
# mean trailing income over the countries, fixed (GROUPS_BY=country).
#
# Second stage, minimum distance with the inverse of the clustered covariance
# of the group coefficients as the weight, xi profiled (given xi the model is
# linear in the rest), no sign restriction (a positive group estimate is
# allowed as noise):
#   simple:     warm arm at the median climate, b_j = beta (Y_j / Ybar)^(xi - 1),
#               Y_j the geometric mean trailing income of group j;
#   simple T:   warm and cold arm, common xi, a beta each;
#   projection: all T coefficients, b = Pi(xi) beta, with Pi(xi) the first-stage
#               regression of the one-stage model's scaled T terms (arms and
#               arms x climate, times (Ybar / Y_it)^(1 - xi)) on the first-stage
#               design: what the group coefficients estimate if the one-stage
#               model holds, so income variation within the groups, and the
#               climate of the groups, are accounted for exactly.
# Inference: Wald SE from the MD Jacobian; distance-metric test of xi = xi0
# (criterion difference, chi2(1)) inverted for a confidence set; the criterion
# at the estimate as the test of the power law (overidentification); country
# bootstrap through both stages (groups, first stage, MD).
#
#   Rscript income_group_xi.R          (GROUPS_J=5 GROUPS_BY=trailing; HANDOUT_GDP, HANDOUT_ROWS as the handout)

# The handout's panel, deviations and test helpers (up to its models), and its
# confidence-set helper.
src <- readLines("handout_preferred_specification.R")
eval(parse(text = src[1:(grep("^models <- list[(]", src) - 1)]))
s0 <- grep("^accept_set <- function", src)
eval(parse(text = src[s0:(s0 - 1 + which(src[s0:length(src)] == "}")[1])]))

J <- as.integer(Sys.getenv("GROUPS_J", "5"))
BY <- Sys.getenv("GROUPS_BY", "trailing")
stopifnot(J >= 3, BY %in% c("trailing", "country"))
GROUPS <- paste0("G", seq_len(J))

# Group indicators: cut points at the quantiles of the data at hand (so redrawn
# in each bootstrap replication).
with_groups <- function(data) {
  key <- if (BY == "trailing") data$y else ave(data$y, data$GID_0)
  cuts <- quantile(if (BY == "trailing") key else key[!duplicated(data$GID_0)], seq(0, 1, length.out = J + 1))
  data$group <- findInterval(key, cuts[-c(1, J + 1)]) + 1L
  for (j in seq_len(J)) data[[GROUPS[j]]] <- as.numeric(data$group == j)
  data
}
T_ARMS <- c(warm = "pTM", cold = "nTM")
group_terms <- list()
for (j in GROUPS) for (a in T_ARMS) {
  group_terms[[paste0(a, "_", j)]] <- c(a, j)
  group_terms[[paste0(a, "_X_clim_TM_", j)]] <- c(a, "clim_TM", j)
}
first_terms <- c(group_terms, arm_terms("RR", c(CLIM[["RR"]], "agr")))
# The one-stage model's T terms (before the scaling).
one_stage_t <- arm_terms("TM", CLIM[["TM"]])

first_stage <- function(data) {
  m <- fit_terms(first_terms, levels = GROUPS[-1], data = data)
  b <- coef(m)[names(group_terms)]
  stopifnot(!anyNA(b))
  list(m = m, b = b, V = vcov(m)[names(group_terms), names(group_terms)], data = data,
       y_group = vapply(seq_len(J), function(j) mean(data$y[data$group == j]), 0))
}
FE_VARS <- c("year", "GID_1")
# Pi(xi): the first-stage coefficients of the group T terms when the outcome is
# each scaled one-stage T term. The first-stage design demeaned once (fixest);
# the scaled terms demeaned for each xi.
projector <- function(fs) {
  data <- fs$data
  for (k in names(first_terms)) data[[k]] <- Reduce(`*`, lapply(first_terms[[k]], function(z) data[[z]]))
  f <- as.formula(paste("g ~", paste(c(names(first_terms), unique(c(LEVELS, GROUPS[-1])), CONTROLS), collapse = " + "), "|", FE))
  m <- feols(f, data = data, panel.id = ~GID_1 + year, demeaned = TRUE, notes = FALSE)
  stopifnot(isTRUE(all.equal(coef(m), coef(fs$m))))
  rows <- obs(m)
  qx <- qr(m$X_demeaned)
  base <- sapply(one_stage_t, function(z) Reduce(`*`, lapply(z, function(v) data[[v]][rows])))
  y <- data$y[rows]
  fv <- data[rows, FE_VARS]
  sv <- data[rows, c("t", "t2")]
  function(xi) {
    Z <- fixest::demean(base * exp(-(1 - xi) * y), f = fv, slope.vars = sv, slope.flag = c(0L, 2L))
    qr.coef(qx, Z)[names(group_terms), , drop = FALSE]
  }
}

# Minimum distance, b = X(xi) beta: beta by GLS given xi, xi by a coarse grid
# over the bounds then a continuous search (as the one-stage profile NLS).
md_fit <- function(b, V, X_of, bounds = XI_BOUNDS, step = 0.25, se = TRUE) {
  W <- solve(V)
  at <- function(xi) {
    X <- X_of(xi)
    beta <- solve(crossprod(X, W %*% X), crossprod(X, W %*% b))
    e <- b - X %*% beta
    list(Q = as.numeric(crossprod(e, W %*% e)), beta = setNames(as.numeric(beta), colnames(X)))
  }
  Q <- function(xi) at(xi)$Q
  xs <- seq(bounds[1], bounds[2], step)
  i <- which.min(vapply(xs, Q, 0))
  xi <- if (i %in% c(1, length(xs))) xs[i] else optimize(Q, xs[i] + c(-1, 1) * step, tol = 1e-5)$minimum
  fit <- at(xi)
  out <- list(xi = xi, beta = fit$beta, Q = fit$Q, df = length(b) - length(fit$beta) - 1, Q_of = Q, at = at)
  if (se && !xi %in% bounds) {
    h <- 1e-4
    G <- cbind(xi = as.vector(X_of(xi + h) %*% fit$beta - X_of(xi - h) %*% fit$beta) / (2 * h), X_of(xi))
    out$vcov <- solve(crossprod(G, W %*% G))
    out$se <- sqrt(diag(out$vcov))
  }
  out
}
XI_BOUNDS <- c(-4, 4)
# The three second stages from a first stage.
warm_terms <- paste0("pTM_", GROUPS)
cold_terms <- paste0("nTM_", GROUPS)
power_x <- function(fs, xi) exp(-(1 - xi) * fs$y_group)
second_stages <- function(fs, se = TRUE) {
  iw <- match(warm_terms, names(fs$b))
  it <- match(c(warm_terms, cold_terms), names(fs$b))
  Pi <- projector(fs)
  list(
    simple = md_fit(fs$b[iw], fs$V[iw, iw], function(xi) cbind(warm = power_x(fs, xi)), se = se),
    simple_T = md_fit(fs$b[it], fs$V[it, it], function(xi) {
      x <- power_x(fs, xi)
      cbind(warm = c(x, 0 * x), cold = c(0 * x, x))
    }, se = se),
    projection = md_fit(fs$b, fs$V, Pi, se = se)
  )
}

fs <- first_stage(with_groups(d))
md <- second_stages(fs)

# Unrestricted group responses at the median climate (pp of growth per +1 SD),
# the group incomes, and Wald tests on them: equal across groups (xi = 1 with
# no climate scaling) and all zero.
groups <- data.frame(
  group = GROUPS, rows = as.vector(table(fs$data$group)),
  countries = vapply(seq_len(J), function(j) length(unique(fs$data$GID_0[fs$data$group == j])), 0),
  income = exp(inc_mean + fs$y_group),
  income_lo = vapply(seq_len(J), function(j) exp(inc_mean + min(fs$data$y[fs$data$group == j])), 0),
  income_hi = vapply(seq_len(J), function(j) exp(inc_mean + max(fs$data$y[fs$data$group == j])), 0),
  warm = 100 * fs$b[warm_terms], warm_se = 100 * sqrt(diag(fs$V)[warm_terms]),
  cold = 100 * fs$b[cold_terms], cold_se = 100 * sqrt(diag(fs$V)[cold_terms]), row.names = NULL)
diff_R <- function(terms) {
  R <- matrix(0, length(terms) - 1, length(coef(fs$m)), dimnames = list(NULL, names(coef(fs$m))))
  for (i in seq_len(nrow(R))) R[i, terms[c(i, i + 1)]] <- c(1, -1)
  R
}
group_tests <- rbind(
  `warm equal across groups` = wald_R(fs$m, diff_R(warm_terms)),
  `cold equal across groups` = wald_R(fs$m, diff_R(cold_terms)),
  `warm zero in every group` = zero_test(fs$m, warm_terms),
  `cold zero in every group` = zero_test(fs$m, cold_terms),
  `T terms equal across groups` = wald_R(fs$m, do.call(rbind, lapply(c(T_ARMS, paste0(T_ARMS, "_X_clim_TM")), function(a)
    diff_R(paste0(a, "_", GROUPS))))))

# Distance-metric tests of xi0, and the confidence sets by inversion.
XI_TESTED <- c(0, 0.6, 1)
XI_FINE <- seq(XI_BOUNDS[1], XI_BOUNDS[2], 0.05)
dm_p <- function(r) function(xi0) pchisq(r$Q_of(xi0) - r$Q, 1, lower.tail = FALSE)
summ <- lapply(md, function(r) {
  p_fun <- dm_p(r)
  p <- vapply(XI_FINE, p_fun, 0)
  list(xi = r$xi, se = if (is.null(r$se)) NA_real_ else r$se[["xi"]], beta = r$beta,
       overid = c(Q = r$Q, df = r$df, p = pchisq(r$Q, r$df, lower.tail = FALSE)),
       p_dm = setNames(vapply(XI_TESTED, p_fun, 0), paste0("xi=", XI_TESTED)),
       dm_set = accept_set(XI_FINE, p, p_fun),
       profile = data.frame(xi = XI_FINE, Q = vapply(XI_FINE, r$Q_of, 0), p = p))
})

# Country bootstrap through both stages: whole countries drawn with replacement
# (each draw a separate country), groups, first stage and the three second
# stages redone. Cached: rerun only if the setup changes.
BOOT <- 999
SEED <- 20261102
tag <- paste0(GDP, if (ROWS == "common") "_common", "_", BY, "_J", J)
out_dir <- file.path("results", "income_group_xi", tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
boot_file <- file.path(out_dir, "bootstrap.rds")
boot_key <- list(rows = nrow(d), terms = names(first_terms), controls = CONTROLS, fe = FE, bounds = XI_BOUNDS,
                 B = BOOT, seed = SEED, xi = round(vapply(md, `[[`, 0, "xi"), 4))
boot <- if (file.exists(boot_file)) readRDS(boot_file)
if (!identical(boot$key, boot_key)) {
  set.seed(SEED)
  countries <- unique(d$GID_0)
  rows_of <- split(seq_len(nrow(d)), d$GID_0)
  t0 <- Sys.time()
  draws <- t(vapply(seq_len(BOOT), function(b) {
    draw <- sample(countries, replace = TRUE)
    db <- d[unlist(rows_of[draw], use.names = FALSE), ]
    db$GID_0 <- db$GID_1 <- rep(paste0(draw, "_", seq_along(draw)), lengths(rows_of[draw]))
    if (b %% 100 == 0) message("bootstrap ", b, "/", BOOT, ", ", format(Sys.time() - t0, digits = 3))
    tryCatch({
      fsb <- first_stage(with_groups(db))
      r <- second_stages(fsb, se = FALSE)
      c(vapply(r, `[[`, 0, "xi"), 100 * fsb$b[warm_terms], 100 * fsb$b[cold_terms])
    }, error = function(e) rep(NA_real_, 3 + 2 * J))
  }, numeric(3 + 2 * J)))
  colnames(draws) <- c(names(md), paste0("warm_", GROUPS), paste0("cold_", GROUPS))
  boot <- list(key = boot_key, draws = draws)
  saveRDS(boot, boot_file)
}
boot_summary <- t(vapply(names(md), function(k) {
  x <- boot$draws[, k]
  x <- x[is.finite(x)]
  bp <- function(x0) min(1, 2 * min(mean(x <= x0), mean(x >= x0)))
  c(reps = length(x), se = sd(x), lo = unname(quantile(x, 0.025)), hi = unname(quantile(x, 0.975)),
    at_bound = mean(x %in% XI_BOUNDS), setNames(vapply(XI_TESTED, bp, 0), paste0("p_xi=", XI_TESTED)))
}, numeric(5 + length(XI_TESTED))))

# The one-stage NLS estimate on the same rows, for comparison.
handout_dir <- file.path("results", "handout_preferred_specification",
                         if (ROWS == "common") paste0(GDP, "_common") else if (GDP != "PWT") GDP else "")
one_stage <- readRDS(file.path(handout_dir, "estimates.rds"))$xi_nls
stopifnot(nrow(one_stage$profile) > 0)

table_xi <- data.frame(
  method = c("one-stage NLS", "MD simple (warm)", "MD simple T (warm + cold)", "MD projection (all T terms)"),
  xi = c(one_stage$estimate, vapply(summ, `[[`, 0, "xi")),
  se_wald = c(one_stage$wald[["se"]], vapply(summ, `[[`, 0, "se")),
  se_boot = c(one_stage$boot[["se"]], boot_summary[, "se"]),
  boot_lo = c(one_stage$boot[["lo"]], boot_summary[, "lo"]),
  boot_hi = c(one_stage$boot[["hi"]], boot_summary[, "hi"]),
  at_bound = c(one_stage$boot[["at_bound"]], boot_summary[, "at_bound"]),
  overid_p = c(NA, vapply(summ, function(s) s$overid[["p"]], 0)),
  p_dm_0 = c(NA, vapply(summ, function(s) s$p_dm[["xi=0"]], 0)),
  p_dm_06 = c(NA, vapply(summ, function(s) s$p_dm[["xi=0.6"]], 0)),
  p_dm_1 = c(NA, vapply(summ, function(s) s$p_dm[["xi=1"]], 0)),
  p_boot_0 = c(one_stage$boot[["p0"]], boot_summary[, "p_xi=0"]),
  p_boot_1 = c(one_stage$boot[["p1"]], boot_summary[, "p_xi=1"]),
  row.names = NULL)

# Fitted warm response at the median climate across income (pp per +1 SD):
# the simple MD, MD with xi fixed at 0 and 1, and the one-stage NLS.
inc_grid <- exp(seq(log(min(groups$income_lo)), log(max(groups$income_hi)), length.out = 200))
y_grid <- log(inc_grid) - inc_mean
gls_at <- function(xi) md$simple$at(xi)$beta[["warm"]]
b1 <- one_stage$coefs
curves <- rbind(
  data.frame(income = inc_grid, curve = sprintf("MD fit, xi = %.2f", md$simple$xi),
             warm = 100 * md$simple$beta[["warm"]] * exp(-(1 - md$simple$xi) * y_grid)),
  data.frame(income = inc_grid, curve = "MD fit, xi = 0", warm = 100 * gls_at(0) * exp(-y_grid)),
  data.frame(income = inc_grid, curve = "MD fit, xi = 1", warm = 100 * gls_at(1)),
  data.frame(income = inc_grid, curve = sprintf("one-stage NLS, xi = %.2f", one_stage$estimate),
             warm = 100 * b1$estimate[b1$term == "pTM:r"] * exp(-(1 - one_stage$estimate) * y_grid)))

saveRDS(list(groups = groups, group_tests = group_tests, md = summ, xi = table_xi, boot = boot_summary,
             boot_draws = boot$draws, curves = curves, first_stage = coeftable(fs$m), J = J, by = BY,
             inc_mean = inc_mean, sample = c(rows = nrow(d), countries = length(unique(d$GID_0))),
             built = Sys.time()),
        file.path(out_dir, "estimates.rds"))
print(groups, digits = 3); print(group_tests, digits = 3)
print(table_xi, digits = 3)
for (k in names(summ)) {
  cat("\n", k, ": beta", round(summ[[k]]$beta, 5), " overid", round(summ[[k]]$overid, 3), "\n")
  print(summ[[k]]$dm_set)
}
