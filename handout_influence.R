# Influence of the countries on the warm and dry responses, and which moderator
# of the warm response, income or the long-run temperature, is the robust one.
# Same models, datasets and design as handout_inference.R (country effects and
# trends partialled out, year dummies; CR1 with fixest's constant, t on G - 1).
#   1. Leave one country out (exact): each country's change of the key
#      coefficients, in SE of the full sample, for
#        Pref (xi = 0): warm arm (pTM:r), dry arm (nRR), dry x share (nRR:agr);
#        M2b:  warm x income (pTM:y), warm x climate (pTM:clim_TM).
#   2. Minimal sets (greedy, exact): drop one at a time the country whose
#      removal moves the coefficient most towards the opposite sign, until it
#      flips; the number dropped when it is no longer significant (5%) and when
#      it flips. For the warm x climate gradient (not significant) the reverse:
#      the countries that move it most towards a loss in hotter climates, until
#      significant.
#   3. Horse race of the warm moderators: the warm arm with income only, with
#      climate only and with both (M2b; the cold arm likewise, P as in M2b),
#      the minimal sets of each gradient, and M2b without the countries that
#      separate income from climate (the 15 richest and 15 poorest for their
#      climate: residuals of mean log income on mean temperature).
#   Rscript handout_influence.R
# Results in results/handout_preferred_specification/influence.rds.

src <- readLines("handout_inference.R")
eval(parse(text = src[seq_len(grep("^inference <- function", src) - 1)]))
STEPS <- 40
N_SEPARATE <- 15

# OLS and CR1 on the countries in keep; the columns zero or collinear there
# (a country's own controls) dropped.
ols <- function(ds, keep) {
  rows <- ds$cluster %in% keep
  X <- ds$X[rows, , drop = FALSE]
  y <- ds$y[rows]
  cl <- ds$cluster[rows]
  X <- X[, colSums(X^2) > 1e-12, drop = FALSE]
  qx <- qr(X, tol = 1e-9)
  if (qx$rank < ncol(X)) X <- X[, -qx$pivot[(qx$rank + 1):ncol(X)], drop = FALSE]
  stopifnot(all(ds$terms %in% colnames(X)))
  A <- solve(crossprod(X))
  b <- drop(A %*% crossprod(X, y))
  u <- y - drop(X %*% b)
  idx <- split(seq_len(nrow(X)), cl)
  G <- length(idx); N <- nrow(X)
  S <- rowsum(X * u, cl)
  V <- G / (G - 1) * (N - 1) / (N - ncol(X) - 2 * G) * A %*% crossprod(S) %*% A
  list(X = X, u = u, A = A, b = b, V = V, idx = idx, G = G, keep = keep)
}
stat <- function(f, term) {
  est <- f$b[[term]]
  se <- sqrt(f$V[term, term])
  c(est = est, se = se, p = 2 * pt(-abs(est / se), f$G - 1))
}
# Leave-one-out change of a coefficient for every country of a fit: the
# closed form, or a refit where a column is the country's own.
loo <- function(ds, f, term) {
  j <- match(term, colnames(f$X))
  vapply(names(f$idx), function(g) {
    i <- f$idx[[g]]
    Xi <- f$X[i, , drop = FALSE]
    out <- tryCatch(-(f$A[j, , drop = FALSE] %*% crossprod(Xi, solve(diag(length(i)) - Xi %*% f$A %*% t(Xi), f$u[i])))[1, 1],
                    error = function(e) NA_real_)
    if (is.na(out)) out <- ols(ds, setdiff(f$keep, g))$b[[term]] - f$b[[term]]
    out
  }, 0)
}
# Greedy minimal set: towards = -1 pushes the coefficient down, +1 up; stops at
# the sign flip ("flip") or at significance in the direction pushed ("signif").
minimal_set <- function(ds, term, towards, until) {
  keep <- unique(ds$cluster)
  rows <- list()
  dropped <- ""
  for (k in 0:STEPS) {
    f <- ols(ds, keep)
    s <- stat(f, term)
    rows[[k + 1]] <- data.frame(k = k, dropped = dropped, t(s))
    done <- if (until == "flip") sign(s[["est"]]) == towards else (s[["p"]] < 0.05 && sign(s[["est"]]) == towards)
    if (done) break
    dlt <- loo(ds, f, term)
    dropped <- names(dlt)[which.max(towards * dlt)]
    keep <- setdiff(keep, dropped)
  }
  path <- do.call(rbind, rows)
  first <- function(ok) if (any(ok)) path$k[which(ok)[1]] else NA
  list(path = path,
       k_insignificant = if (until == "flip") first(path$p >= 0.05) else NA,
       k_flip = if (until == "flip") first(sign(path$est) == towards) else NA,
       k_significant = if (until == "signif") first(path$p < 0.05 & sign(path$est) == towards) else NA)
}

# Country traits over the rows: mean income (US$), temperature (C), share.
traits <- function(ds) {
  x <- ds$data
  data.frame(GID_0 = sort(unique(x$GID_0)),
             income = tapply(exp(x$inc_trailing), x$GID_0, mean)[sort(unique(x$GID_0))],
             temp = tapply(x$mean_TM_lag_20, x$GID_0, mean)[sort(unique(x$GID_0))],
             agr = tapply(x$agr_level, x$GID_0, mean)[sort(unique(x$GID_0))], row.names = NULL)
}
influence_table <- function(ds, term, tr) {
  f <- ols(ds, unique(ds$cluster))
  s <- stat(f, term)
  dlt <- loo(ds, f, term)
  out <- data.frame(GID_0 = names(dlt), change = dlt, change_se = dlt / s[["se"]], row.names = NULL)
  out <- merge(out, tr, by = "GID_0")
  out[order(-abs(out$change_se)), ]
}

run <- function(gdp, climate) {
  env <- new.env()
  env$CLIMATE <- climate
  Sys.setenv(HANDOUT_GDP = gdp, HANDOUT_ROWS = "own", HANDOUT_OUTLIERS = "5", HANDOUT_MIN_INCOME = "0")
  eval(parse(text = code), envir = env)
  clim <- env$CLIM
  p_terms <- env$arm_terms("RR", c(clim[["RR"]], "agr"))
  ds <- list(
    Pref = design(env, env$scaled_terms, env$scaled_data(0)),
    M2b = design(env, c(env$arm_terms("TM", c(clim[["TM"]], "y")), p_terms), env$d),
    `T income only` = design(env, c(env$arm_terms("TM", "y"), p_terms), env$d),
    `T climate only` = design(env, c(env$arm_terms("TM", clim[["TM"]]), p_terms), env$d))
  tr <- traits(ds$M2b)

  targets <- list(c("Pref", "pTM_X_r"), c("Pref", "nRR"), c("Pref", "nRR_X_agr"),
                  c("M2b", "pTM_X_y"), c("M2b", "pTM_X_clim_TM"))
  infl <- lapply(targets, function(z) influence_table(ds[[z[1]]], z[2], tr))
  names(infl) <- vapply(targets, paste, "", collapse = " ")

  # Minimal sets: significant coefficients towards zero; the warm x climate
  # gradient towards a loss in hotter climates (negative).
  sets <- list(
    `Pref warm arm` = c("Pref", "pTM_X_r", "flip"), `Pref dry arm` = c("Pref", "nRR", "flip"),
    `Pref dry x share` = c("Pref", "nRR_X_agr", "flip"),
    `M2b warm x income` = c("M2b", "pTM_X_y", "flip"), `M2b warm x climate` = c("M2b", "pTM_X_clim_TM", "signif"),
    `Income only: warm x income` = c("T income only", "pTM_X_y", "flip"),
    `Climate only: warm x climate` = c("T climate only", "pTM_X_clim_TM", "flip"))
  minimal <- lapply(sets, function(z) {
    b0 <- ols(ds[[z[1]]], unique(ds[[z[1]]]$cluster))$b[[z[2]]]
    towards <- if (z[3] == "flip") -sign(b0) else -1
    minimal_set(ds[[z[1]]], z[2], towards, z[3])
  })

  # The countries that separate income from climate.
  res <- resid(lm(log(income) ~ temp, tr))
  separate <- list(`rich for their climate` = tr$GID_0[order(-res)][1:N_SEPARATE],
                   `poor for their climate` = tr$GID_0[order(res)][1:N_SEPARATE])
  separate$both <- c(separate[[1]], separate[[2]])
  all_c <- unique(ds$M2b$cluster)
  race <- do.call(rbind, lapply(c(list(none = character(0)), separate), function(out) {
    f <- ols(ds$M2b, setdiff(all_c, out))
    data.frame(t(stat(f, "pTM_X_y")), t(stat(f, "pTM_X_clim_TM")))
  }))
  names(race) <- c("income_est", "income_se", "income_p", "climate_est", "climate_se", "climate_p")
  race$dropped <- c("", vapply(separate, paste, "", collapse = " "))
  list(influence = infl, minimal = minimal, race = race, separate = separate, traits = tr,
       cor_income_temp = cor(log(tr$income), tr$temp))
}

out <- list()
for (k in seq_len(nrow(DATASETS))) {
  key <- paste(DATASETS$gdp[k], sub("UDelaware", "UDel", sub("CRU TS", "CRU", DATASETS$climate[k])))
  t0 <- Sys.time()
  out[[key]] <- run(DATASETS$gdp[k], DATASETS$climate[k])
  message(key, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
saveRDS(list(results = out, steps = STEPS, n_separate = N_SEPARATE),
        file.path("results", "handout_preferred_specification", "influence.rds"))
