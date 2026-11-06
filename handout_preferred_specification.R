# Estimates for the handout on the preferred specification
# (handout_preferred_specification.Rmd), saved to
# results/handout_preferred_specification/estimates.rds.
#
# Base: PWT 11.0 x ERA5 (population weights), the response explorer's panel
# (results/response_explorer/panels.rds; run apps/response_explorer/
# build_models.R first). Deviations from the 20-year trailing mean in SD of the
# country's variability over the horizon (fixed SD), demeaned by the country's
# linear and quadratic trend; linear arms; climate of the interactions = the
# same 20-year trailing mean, income = log mean GDP per capita of the 5 years
# before; agricultural share = long-run WDI share. Year effects, country
# effects and country quadratic trends, controls, errors clustered by country.
# Every model runs on the same rows (those with the trailing income).
#
#   Rscript handout_preferred_specification.R

suppressPackageStartupMessages(library(fixest))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- file.path("apps", "response_explorer")
source(file.path(app_dir, "R", "response.R"))
source(file.path(app_dir, "R", "fit.R"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))

# The explorer's model data (moderators centred as in the app) without the fit.
lines <- readLines(file.path(app_dir, "R", "fit.R"))
s0 <- grep("^explorer_model <- function", lines)
s1 <- s0 - 1 + which(lines[s0:length(lines)] == "}")[1]
fn <- paste(lines[s0:s1], collapse = "\n")
fn <- sub("explorer_model <- function", "model_data <- function", fn)
fn <- sub("if (!is.null(model_cache[[key]])) return(model_cache[[key]])", "", fn, fixed = TRUE)
fn <- sub("m <- fit(o$shape)", "return(d)", fn, fixed = TRUE)
stopifnot(grepl("return(d)", fn, fixed = TRUE))
eval(parse(text = fn))

# GDP source: PWT (the handout's), or WB with HANDOUT_GDP=WB (everything else
# the same; results in a WB subfolder). HANDOUT_ROWS=common: only the rows both
# sources have (the weather, climate, share and controls are the same there;
# outcome and income are the source's), results in <GDP>_common.
GDP <- Sys.getenv("HANDOUT_GDP", "PWT")
ROWS <- Sys.getenv("HANDOUT_ROWS", "own")
# HANDOUT_OUTLIERS: the growth outliers dropped, the non-weather ones of the 5%
# (the handout's) or 1% tails of growth in the full panel, or ("r5", "r1")
# those of the tails of growth reversals (the notebook's rule since
# 2026-11-04). HANDOUT_MIN_INCOME: only the countries whose mean trailing
# income over their rows is above it (2021 PPP US$). Results in a subfolder
# named by both.
OUTLIERS <- Sys.getenv("HANDOUT_OUTLIERS", "5")
MIN_INCOME <- as.numeric(Sys.getenv("HANDOUT_MIN_INCOME", "0"))
stopifnot(GDP %in% c("PWT", "WB"), ROWS %in% c("own", "common"), OUTLIERS %in% c("1", "5", "r1", "r5"))
data_of <- function(gdp) {
  model_data(list(gdp = gdp, climate = "ERA5", weight = "population", units = "absolute",
                  baseline = "20-year", form = "Linear arms", shape = "linear", fe = "quadratic",
                  moderators = "trailing", outcome = "log", outliers = OUTLIERS,
                  inter = list(TM = c("clim", "inc", "agr"), RR = c("clim", "inc", "agr"))), P)
}
d <- data_of(GDP)
if (ROWS == "common") {
  other <- data_of(setdiff(c("PWT", "WB"), GDP))
  d <- d[paste(d$GID_0, d$year) %in% paste(other$GID_0, other$year), ]
}
if (MIN_INCOME > 0) {
  mean_income <- tapply(exp(d$inc_trailing), d$GID_0, mean)
  d <- d[d$GID_0 %in% names(mean_income)[mean_income > MIN_INCOME], ]
}
if (ROWS == "common" || MIN_INCOME > 0) {
  # Income and share less their means over these rows, as on the own rows.
  d$inc <- d$inc_trailing - mean(d$inc_trailing)
  d$agr <- d$agr_level - mean(d$agr_level)
}
# Deviations from the trailing mean in SD of the horizon (fixed SD): the
# demeaned absolute deviation over the country's horizon SD (a constant within
# the country, so dividing commutes with the demeaning).
stopifnot(all(tapply(d$sd_TM_all, d$GID_1, function(x) diff(range(x))) == 0))
for (v in c("TM", "RR")) {
  x <- d[[paste0("dm_", v, "_lag_20_unstd_quadratic")]] / d[[paste0("sd_", v, "_all")]]
  d[[paste0("d", v)]] <- x
  d[[paste0("p", v)]] <- pmax(x, 0)
  d[[paste0("n", v)]] <- pmax(-x, 0)
  # Quadratic arms: the arms squared.
  d[[paste0("p", v, "2")]] <- pmax(x, 0)^2
  d[[paste0("n", v, "2")]] <- pmax(-x, 0)^2
}
d$y <- d$inc                         # log trailing income less its mean over the rows
d$g <- d[[P$outcome]]
inc_mean <- mean(d$inc_trailing)
stopifnot(all(is.finite(d$dTM)), all(is.finite(d$dRR)), all(is.finite(d$y)))
CONTROLS <- P$controls
FE <- "year + GID_1[t, t2]"

# A model from a list of terms: name = product of columns. Level terms are
# columns of d. Returns the fit and its pieces.
# Every model controls for the levels of the trailing climate (T, P) and of the
# trailing income (LEVELS), so that the models differ only in the weather terms.
LEVELS <- c("clim_TM", "clim_RR", "y")
fit_terms <- function(terms, levels = character(0), data = d) {
  levels <- unique(c(LEVELS, levels))
  for (k in names(terms)) data[[k]] <- Reduce(`*`, lapply(terms[[k]], function(z) data[[z]]))
  f <- as.formula(paste("g ~", paste(c(names(terms), levels, CONTROLS), collapse = " + "), "|", FE))
  m <- feols(f, data = data, cluster = ~GID_0, panel.id = ~GID_1 + year, notes = FALSE)
  m$country <- data$GID_0[obs(m)]
  m
}
arm_terms <- function(v, modifiers = character(0), scale = NULL, symmetric = FALSE, power = "") {
  arms <- if (symmetric) paste0("d", v) else paste0(c("p", "n"), v, power)
  out <- list()
  for (a in arms) {
    out[[a]] <- c(a, scale)
    for (z in modifiers) out[[paste0(a, ":", z)]] <- c(a, z, scale)
  }
  if (!is.null(scale)) names(out) <- paste0(names(out), ":r")
  names(out) <- gsub(":", "_X_", names(out))
  out
}
CLIM <- c(TM = "clim_TM", RR = "clim_RR")

# Clustered Wald test of R b = 0 (rows of R named by coefficient).
wald_R <- function(m, R) {
  b <- coef(m)[colnames(R)]
  V <- vcov(m)[colnames(R), colnames(R)]
  Rb <- R %*% b
  stat <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R), Rb)) / nrow(R)
  G <- length(unique(m$country))
  c(F = stat, df = nrow(R), p = pf(stat, nrow(R), G - 1, lower.tail = FALSE))
}
zero_test <- function(m, terms) {
  R <- matrix(0, length(terms), length(coef(m)), dimnames = list(terms, names(coef(m))))
  for (k in terms) R[k, k] <- 1
  wald_R(m, R)
}

models <- list(
  M0 = fit_terms(c(arm_terms("TM", symmetric = TRUE), arm_terms("RR", symmetric = TRUE))),
  M1 = fit_terms(c(arm_terms("TM"), arm_terms("RR"))),
  M0c = fit_terms(c(arm_terms("TM", CLIM[["TM"]], symmetric = TRUE), arm_terms("RR", CLIM[["RR"]], symmetric = TRUE))),
  M1c = fit_terms(c(arm_terms("TM", CLIM[["TM"]]), arm_terms("RR", CLIM[["RR"]]))),
  M2 = fit_terms(c(arm_terms("TM", c(CLIM[["TM"]], "y", "agr")), arm_terms("RR", c(CLIM[["RR"]], "y", "agr")))),
  M2b = fit_terms(c(arm_terms("TM", c(CLIM[["TM"]], "y")), arm_terms("RR", c(CLIM[["RR"]], "agr"))))
)

# Income scaling of the T terms: times r = (Ybar / Y_it)^(1 - xi) = exp(-(1 - xi) y)
# (Ybar: the geometric mean of the trailing income over the rows). P: arms with
# climate and share. The level controls are those of every model (LEVELS).
scaled_data <- function(xi, data = d) {
  data$r <- exp(-(1 - xi) * data$y)
  data
}
scaled_terms <- c(arm_terms("TM", CLIM[["TM"]], scale = "r"), arm_terms("RR", c(CLIM[["RR"]], "agr")))
fit_scaled <- function(xi, data = d) fit_terms(scaled_terms, data = scaled_data(xi, data))
# The Gauss-Newton regression at xi0: the model at xi0 with h, the derivative
# of its fitted T part in xi (y times the scaled T terms times their
# coefficients). The coefficient of h is the Gauss-Newton step in xi; its
# clustered t is the C(alpha) score test of xi = xi0; at the NLS estimate (step
# 0) its clustered SE is that of xi hat in the full nonlinear model, which the
# conditional regression at xi hat leaves out.
gn_fit <- function(xi0) {
  dd <- scaled_data(xi0)
  m <- fit_scaled(xi0)
  b <- coef(m)
  t_terms <- grep("_X_r$", names(b), value = TRUE)
  for (k in names(scaled_terms)) dd[[k]] <- Reduce(`*`, lapply(scaled_terms[[k]], function(z) dd[[z]]))
  dd$h <- dd$y * as.numeric(as.matrix(dd[t_terms]) %*% b[t_terms])
  f <- as.formula(paste("g ~", paste(c(names(scaled_terms), LEVELS, "h", CONTROLS), collapse = " + "), "|", FE))
  m2 <- feols(f, data = dd, cluster = ~GID_0, panel.id = ~GID_1 + year, notes = FALSE)
  m2$country <- dd$GID_0[obs(m2)]
  m2
}
xi_score_p <- function(xi0) unname(pvalue(gn_fit(xi0))[["h"]])
XI_GRID <- seq(0, 1, 0.1)
scaled <- lapply(XI_GRID, fit_scaled)
names(scaled) <- sprintf("xi=%.1f", XI_GRID)

# Continuous xi by nonlinear least squares, profiled: given xi the model is
# linear in all other coefficients (fixed effects and trends included), so the
# NLS estimate minimises over xi the residual sum of squares of fit_scaled(xi).
# xi is not restricted to [0, 1]; a coarse grid over wide bounds guards against
# local minima, then a continuous search around the best grid point.
XI_BOUNDS <- c(-4, 4)
xi_nls <- function(data = d, step = 0.5) {
  rss <- function(x) sum(resid(fit_scaled(x, data))^2)
  xs <- seq(XI_BOUNDS[1], XI_BOUNDS[2], step)
  i <- which.min(vapply(xs, rss, 0))
  if (i %in% c(1, length(xs))) return(xs[i])
  optimize(rss, xs[i] + c(-1, 1) * step, tol = 1e-5)$minimum
}
xi_hat <- xi_nls()
gn <- gn_fit(xi_hat)
# No step left at an interior optimum (at a bound the Wald SE does not apply).
stopifnot(xi_hat %in% XI_BOUNDS || abs(coef(gn)[["h"]]) < 0.01 * se(gn)[["h"]])
G <- length(unique(d$GID_0))
xi_se <- unname(se(gn)[["h"]])
xi_wald <- c(estimate = xi_hat, se = xi_se,
             lo = xi_hat - qt(0.975, G - 1) * xi_se, hi = xi_hat + qt(0.975, G - 1) * xi_se,
             p0 = 2 * pt(-abs(xi_hat / xi_se), G - 1), p1 = 2 * pt(-abs((xi_hat - 1) / xi_se), G - 1))
# A test of xi = xi0 robust to weak identification (Anderson-Rubin type): the
# model at xi0 with every scaled T term also times y; under xi = xi0 the model
# is linear with these 4 terms zero, whatever the size of the T response, so
# their clustered Wald test has the right size even where the response is near
# zero (where the score test has none) and does not accept at maxima of the
# profile (as the score test does: it only tests the slope).
t_scaled <- scaled_terms[grep("^[pn]TM", names(scaled_terms))]
ar_terms <- stats::setNames(lapply(t_scaled, c, "y"), paste0(names(t_scaled), "_X_y"))
xi_ar_p <- function(xi0) {
  m <- fit_terms(c(scaled_terms, ar_terms), data = scaled_data(xi0))
  unname(zero_test(m, names(ar_terms))[["p"]])
}
# The profile objective (AIC: same rows and terms, so n log RSS plus a
# constant) and the two tests over a fine grid; the 95% confidence set by test
# inversion = the xi not rejected, as intervals with endpoints refined by root
# finding (an interval reaching a bound of the grid is open).
XI_FINE <- seq(XI_BOUNDS[1], XI_BOUNDS[2], 0.05)
profile <- data.frame(xi = XI_FINE, aic = vapply(XI_FINE, function(x) AIC(fit_scaled(x)), 0),
                      p_score = vapply(XI_FINE, xi_score_p, 0), p_ar = vapply(XI_FINE, xi_ar_p, 0))
accept_set <- function(xs, p, p_fun, level = 0.05) {
  runs <- rle(p >= level)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  edge <- function(i, j) uniroot(function(x) p_fun(x) - level, xs[c(i, j)], tol = 1e-4)$root
  do.call(rbind, lapply(which(runs$values), function(k) {
    a <- starts[k]; b <- ends[k]
    data.frame(lo = if (a == 1) -Inf else edge(a - 1, a), hi = if (b == length(xs)) Inf else edge(b, b + 1))
  }))
}
xi_score_set <- accept_set(profile$xi, profile$p_score, xi_score_p)
xi_ar_set <- accept_set(profile$xi, profile$p_ar, xi_ar_p)
xi_ar <- c(p0 = xi_ar_p(0), p1 = xi_ar_p(1), p_hat = xi_ar_p(xi_hat))
# On the common rows, the PWT run: the two parts of the change of source, the
# outcome and the income of the scaling, swapped for WB's one at a time (xi
# with a finer coarse grid), and the spread of WB log income against PWT's.
xi_swaps <- NULL
if (ROWS == "common" && GDP == "PWT") {
  i <- match(paste(d$GID_0, d$year), paste(other$GID_0, other$year))
  wb_g <- other[[P$outcome]][i]
  wb_y <- other$inc_trailing[i] - mean(other$inc_trailing[i])
  variant <- function(g, y) {
    x <- d
    x$g <- g
    x$y <- y
    x
  }
  xi_swaps <- list(
    xi = c(`PWT outcome, PWT income` = xi_nls(variant(d$g, d$y), step = 0.25),
           `WB outcome, PWT income` = xi_nls(variant(wb_g, d$y), step = 0.25),
           `PWT outcome, WB income` = xi_nls(variant(d$g, wb_y), step = 0.25),
           `WB outcome, WB income` = xi_nls(variant(wb_g, wb_y), step = 0.25)),
    income_spread = unname(coef(lm(wb_y ~ d$y))[2]))
}

# Country bootstrap: whole countries drawn with replacement (each draw a
# separate country, for its fixed effects and trends), xi re-estimated by
# profile NLS in each replication. Cached: rerun only if the setup changes.
XI_BOOT <- 999
XI_SEED <- 20261102
out_dir <- file.path("results", "handout_preferred_specification")
if (ROWS == "common") {
  out_dir <- file.path(out_dir, paste0(GDP, "_common"))
} else if (GDP != "PWT") {
  out_dir <- file.path(out_dir, GDP)
}
if (OUTLIERS != "5" || MIN_INCOME > 0) {
  out_dir <- file.path(out_dir, paste0("outliers", OUTLIERS, if (MIN_INCOME > 0) paste0("_income", MIN_INCOME)))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
boot_file <- file.path(out_dir, "xi_bootstrap.rds")
boot_key <- list(rows = nrow(d), terms = names(scaled_terms), controls = CONTROLS, fe = FE,
                 bounds = XI_BOUNDS, B = XI_BOOT, seed = XI_SEED, xi_hat = round(xi_hat, 4))
boot <- if (file.exists(boot_file)) readRDS(boot_file)
if (!identical(boot$key, boot_key)) {
  set.seed(XI_SEED)
  countries <- unique(d$GID_0)
  rows_of <- split(seq_len(nrow(d)), d$GID_0)
  t0 <- Sys.time()
  xi_b <- vapply(seq_len(XI_BOOT), function(b) {
    draw <- sample(countries, replace = TRUE)
    db <- d[unlist(rows_of[draw], use.names = FALSE), ]
    db$GID_0 <- db$GID_1 <- rep(paste0(draw, "_", seq_along(draw)), lengths(rows_of[draw]))
    if (b %% 100 == 0) message("bootstrap ", b, "/", XI_BOOT, ", ", format(Sys.time() - t0, digits = 3))
    tryCatch(xi_nls(db), error = function(e) NA_real_)
  }, 0)
  boot <- list(key = boot_key, xi = xi_b)
  saveRDS(boot, boot_file)
}
xi_b <- boot$xi[is.finite(boot$xi)]
boot_p <- function(x0) min(1, 2 * min(mean(xi_b <= x0), mean(xi_b >= x0)))
xi_boot <- c(reps = length(xi_b), failed = sum(!is.finite(boot$xi)), se = sd(xi_b),
             lo = unname(quantile(xi_b, 0.025)), hi = unname(quantile(xi_b, 0.975)),
             at_bound = mean(xi_b %in% XI_BOUNDS), p0 = boot_p(0), p1 = boot_p(1))
xi_nls_out <- list(estimate = xi_hat, wald = xi_wald, score_set = xi_score_set, ar_set = xi_ar_set,
                   ar = xi_ar, profile = profile,
                   boot = xi_boot, boot_draws = boot$xi, bounds = XI_BOUNDS, swaps = xi_swaps)
grid <- data.frame(
  xi = XI_GRID,
  aic = vapply(scaled, AIC, 0),
  bic = vapply(scaled, BIC, 0),
  war2 = vapply(scaled, function(m) fitstat(m, "war2")[[1]], 0),
  p_score = vapply(XI_GRID, xi_score_p, 0)
)
XI_PREFERRED <- 0
models$Pref <- fit_scaled(XI_PREFERRED)
# The preferred specification with quadratic arms (the arms squared, with the
# same interactions and scaling), and with both (the encompassing model).
quadratic_terms <- c(arm_terms("TM", CLIM[["TM"]], scale = "r", power = "2"),
                     arm_terms("RR", c(CLIM[["RR"]], "agr"), power = "2"))
models$PrefQ <- fit_terms(quadratic_terms, data = scaled_data(XI_PREFERRED))
models$PrefLQ <- fit_terms(c(scaled_terms, quadratic_terms), data = scaled_data(XI_PREFERRED))

# Fit verdicts against noise (clustered Vuong, 5%), as in the app.
verdict <- function(a, b) fit_verdict(fit_record(a), fit_record(b), "aic")
best <- names(scaled)[which.min(grid$aic)]
grid$verdict_vs_best <- vapply(scaled, function(m) if (identical(m, scaled[[best]])) "best" else verdict(m, scaled[[best]]), "")
grid$verdict_vs_M2b <- vapply(scaled, function(m) verdict(m, models$M2b), "")

# Symmetry: the warm and cold arm (wet and dry) with opposite coefficients, for
# the arm and for each interaction: b_p + b_n = 0.
symmetry_test <- function(m, v, modifiers = character(0)) {
  suffixes <- c("", if (length(modifiers)) paste0("_X_", modifiers))
  pairs <- lapply(suffixes, function(s) paste0(c("p", "n"), v, s))
  R <- matrix(0, length(pairs), length(coef(m)), dimnames = list(NULL, names(coef(m))))
  for (i in seq_along(pairs)) R[i, pairs[[i]]] <- 1
  wald_R(m, R)
}
tests <- list(
  symmetry = rbind(
    `T, no interaction` = symmetry_test(models$M1, "TM"), `P, no interaction` = symmetry_test(models$M1, "RR"),
    `T, climate interaction` = symmetry_test(models$M1c, "TM", CLIM[["TM"]]),
    `P, climate interaction` = symmetry_test(models$M1c, "RR", CLIM[["RR"]])),
  arms = rbind(
    `M1: T terms` = zero_test(models$M1, c("pTM", "nTM")), `M1: P terms` = zero_test(models$M1, c("pRR", "nRR")),
    `M1c: T terms` = zero_test(models$M1c, grep("TM", names(coef(models$M1c)), value = TRUE)[1:4]),
    `M1c: P terms` = zero_test(models$M1c, grep("RR", names(coef(models$M1c)), value = TRUE)[1:4])),
  interactions = do.call(rbind, lapply(c(TM = "T", RR = "P"), function(s) {
    v <- if (s == "T") "TM" else "RR"
    of <- function(z) paste0(c("p", "n"), v, "_X_", z)
    rbind(stats::setNames(zero_test(models$M2, of(CLIM[[v]])), c("F", "df", "p")),
          zero_test(models$M2, of("y")), zero_test(models$M2, of("agr")))
  }))
)
rownames(tests$interactions) <- paste(rep(c("T", "P"), each = 3), c("climate", "income", "share"))
tests$M2b_vs_M2 <- zero_test(models$M2, c(paste0(c("pTM", "nTM"), "_X_agr"), paste0(c("pRR", "nRR"), "_X_y")))
# Linear against quadratic arms in the encompassing model: the squared terms
# zero (linear arms suffice) or the linear terms zero (quadratic arms suffice).
of_var <- function(terms, v) grep(paste0("^[pn]", v), names(terms), value = TRUE)
tests$quadratic <- rbind(
  `Linear suffices: squared terms = 0, T and P` = zero_test(models$PrefLQ, names(quadratic_terms)),
  `Linear suffices: squared terms = 0, T` = zero_test(models$PrefLQ, of_var(quadratic_terms, "TM")),
  `Linear suffices: squared terms = 0, P` = zero_test(models$PrefLQ, of_var(quadratic_terms, "RR")),
  `Quadratic suffices: linear terms = 0, T and P` = zero_test(models$PrefLQ, names(scaled_terms)),
  `Quadratic suffices: linear terms = 0, T` = zero_test(models$PrefLQ, of_var(scaled_terms, "TM")),
  `Quadratic suffices: linear terms = 0, P` = zero_test(models$PrefLQ, of_var(scaled_terms, "RR")))

# Fit of every model.
fit_table <- do.call(rbind, lapply(models, function(m) {
  data.frame(rows = nobs(m), countries = length(unique(m$country)), terms = length(coef(m)),
             aic = AIC(m), bic = BIC(m), war2 = fitstat(m, "war2")[[1]])
}))
verdicts <- c(`M1 vs M0` = verdict(models$M1, models$M0), `M1c vs M0c` = verdict(models$M1c, models$M0c),
              `M2 vs M1c` = verdict(models$M2, models$M1c), `M2b vs M2` = verdict(models$M2b, models$M2),
              `Pref vs M2b` = verdict(models$Pref, models$M2b), `Pref vs M1c` = verdict(models$Pref, models$M1c),
              `PrefQ vs Pref` = verdict(models$PrefQ, models$Pref))

# Coefficients (all, with clustered SE and p), level terms included; controls dropped.
coef_of <- function(m) {
  ct <- coeftable(m)
  ct <- ct[!grepl("conflict|collapse", rownames(ct)), , drop = FALSE]
  data.frame(term = gsub("_X_", ":", rownames(ct)), estimate = ct[, 1], se = ct[, 2], p = ct[, 4], row.names = NULL)
}
coefs <- lapply(models, coef_of)
coefs_scaled <- lapply(scaled[c("xi=0.0", "xi=0.6", "xi=1.0")], coef_of)
# The model at the NLS estimate with SEs of the full nonlinear model (xi
# estimated): the Gauss-Newton regression at xi hat; h is the xi step.
xi_nls_out$coefs <- coef_of(gn)

# Responses (pp of growth) to +-1 SD in the countries' last sample year (their
# trailing climate and income, long-run share), with clustered SEs; M2b's income
# at which the warm (cold) arm's response changes sign at the median climate
# and mean share.
last <- do.call(rbind, lapply(split(d, d$GID_0), function(x) x[which.max(x$year), ]))
response <- function(m, v, arm, row, xi = NULL, power = "", size = 1) {
  b <- coef(m)
  a <- paste0(arm, v, power)
  w <- stats::setNames(numeric(length(b)), names(b))
  vals <- c(clim_TM = row$clim_TM, clim_RR = row$clim_RR, y = row$y, agr = row$agr)
  r <- if (!is.null(xi)) exp(-(1 - xi) * row$y) else 1
  for (k in names(b)) {
    parts <- strsplit(k, "_X_", fixed = TRUE)[[1]]
    if (parts[[1]] != a) next
    mods <- setdiff(parts[-1], "r")
    w[[k]] <- (if (length(mods)) prod(vals[mods]) else 1) * (if ("r" %in% parts) r else 1) *
      (if (power == "2") size^2 else size)
  }
  est <- sum(w * b)
  se <- sqrt(as.numeric(t(w) %*% vcov(m) %*% w))
  c(estimate = 100 * est, se = 100 * se)
}
COUNTRIES <- intersect(c("BDI", "ETH", "IND", "NGA", "IDN", "BRA", "CHN", "ITA", "USA"), last$GID_0)
responses <- do.call(rbind, lapply(COUNTRIES, function(iso) {
  row <- last[last$GID_0 == iso, ]
  do.call(rbind, lapply(list(M1c = list(models$M1c, NULL), M2b = list(models$M2b, NULL),
                             Pref = list(models$Pref, XI_PREFERRED)), function(z) {
    out <- c(warm = response(z[[1]], "TM", "p", row, z[[2]]), cold = response(z[[1]], "TM", "n", row, z[[2]]),
             wet = response(z[[1]], "RR", "p", row, z[[2]]), dry = response(z[[1]], "RR", "n", row, z[[2]]))
    data.frame(t(out))
  })) -> x
  data.frame(iso = iso, year = row$year, income = exp(row$inc_trailing), model = c("M1c", "M2b", "Pref"), x,
             row.names = NULL)
}))
# Linear against quadratic preferred specification: responses to 1 and 2 SD.
responses_q <- do.call(rbind, lapply(COUNTRIES, function(iso) {
  row <- last[last$GID_0 == iso, ]
  do.call(rbind, lapply(c(1, 2), function(size) {
    do.call(rbind, lapply(c(Pref = "", PrefQ = "2"), function(power) {
      m <- if (power == "") models$Pref else models$PrefQ
      out <- c(warm = response(m, "TM", "p", row, XI_PREFERRED, power, size),
               cold = response(m, "TM", "n", row, XI_PREFERRED, power, size),
               wet = response(m, "RR", "p", row, NULL, power, size),
               dry = response(m, "RR", "n", row, NULL, power, size))
      data.frame(iso = iso, size = size, model = if (power == "") "Pref" else "PrefQ", t(out))
    }))
  }))
}))
b2 <- coef(models$M2b)
threshold <- c(warm = exp(inc_mean - b2[["pTM"]] / b2[["pTM_X_y"]]), cold = exp(inc_mean - b2[["nTM"]] / b2[["nTM_X_y"]]))
above <- vapply(threshold, function(th) mean(exp(last$inc_trailing) > th), 0)

vcovs <-lapply(models, function(m) {
  V <- vcov(m)
  dimnames(V) <- lapply(dimnames(V), function(x) gsub("_X_", ":", x))
  V
})
# The preferred models in the response explorer's format, for its projection
# code (R/projection.R): the deviation terms only (named as the app's:
# devTM_pos:clim_TM:idol, ...), the centring of the characteristics, the
# scaling of T (xi) and the windows (20-year trailing climate, 5-year income;
# the share follows the SSP as a 20-year trailing mean).
units_20 <- P$units[[paste(GDP, "ERA5 | 20-year | population | absolute", sep = " | ")]]
clim_median <- c(TM = stats::median(units_20$mean_T), RR = stats::median(units_20$mean_P))
stopifnot(isTRUE(all.equal(d$clim_TM, d$mean_TM_lag_20 - clim_median[["TM"]])))
app_name <- function(x) {
  x <- gsub("_X_", ":", x)
  x <- sub("^p(TM|RR)(2?)", "dev\\1_pos\\2", x)
  x <- sub("^n(TM|RR)(2?)", "dev\\1_neg\\2", x)
  sub(":r$", ":idol", x)
}
to_app <- function(m, key, form) {
  keep <- grep("^[pn](TM|RR)", names(coef(m)), value = TRUE)
  V <- vcov(m)[keep, keep]
  dimnames(V) <- list(app_name(keep), app_name(keep))
  list(key = key, gdp = GDP, climate = "ERA5", weight = "population", units = "SD", baseline = "20-year",
       window = "lag_20", form = form, shape = "linear", fe = "quadratic", moderators = "trailing",
       inter = list(TM = "clim", RR = c("clim", "agr")), outcome = "log",
       income_form = c(TM = "dollars", RR = "standard"), xi = c(TM = XI_PREFERRED, RR = 1),
       coef = stats::setNames(coef(m)[keep], app_name(keep)), vcov = V,
       clim_median = clim_median, inc_mean = inc_mean, agr_mean = mean(d$agr_level),
       horizon = P$panels[[paste(GDP, "ERA5 | population", sep = " | ")]]$horizon, inc_window = P$trailing_income_years,
       share_window = 20L, nobs = nobs(m), countries = length(unique(m$country)))
}
app_models <- list(Pref = to_app(models$Pref, "handout | Pref | linear", "Linear arms"),
                   PrefQ = to_app(models$PrefQ, "handout | PrefQ | quadratic", "Quadratic arms"))
saveRDS(app_models, file.path(out_dir, "app_models.rds"))

saveRDS(list(coefs = coefs, coefs_scaled = coefs_scaled, vcov = vcovs, last = last[, c("GID_0", "year", "y", "inc_trailing")],
             responses_q = responses_q,
             fit = fit_table, verdicts = verdicts, tests = tests,
             grid = grid, xi_hat = xi_hat, xi_nls = xi_nls_out, xi_preferred = XI_PREFERRED, responses = responses,
             threshold = threshold, share_above = above, inc_mean = inc_mean,
             sample = c(rows = nrow(d), countries = length(unique(d$GID_0)), first = min(d$year), last = max(d$year)),
             horizon = P$panels[[paste(GDP, "ERA5 | population", sep = " | ")]]$horizon, built = Sys.time()),
        file.path(out_dir, "estimates.rds"))
print(fit_table); print(verdicts); print(tests); print(grid)
cat("xi hat (NLS):", xi_hat, "\n"); print(xi_wald); print(xi_score_set); print(xi_ar_set); print(xi_ar); print(xi_boot)
print(threshold); print(above)
print(responses[responses$iso %in% c("BDI", "IND", "BRA", "USA"), ])
