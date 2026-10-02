# Inference for the fast-versus-slow warming design (sections 6, 7 and 10 of
# FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md).
#
# Every quantity of interest (a coefficient, a path gap at one horizon, a local
# projection response) is a linear contrast C %*% beta of a fixed-effects OLS
# fit, so all inference works on the Frisch-Waugh-Lovell (FWL) form of the fit:
# the regressors demeaned by the model's own fixed effects and varying slopes,
# X, and the residuals u. With A = X'WX:
#   beta_hat - beta = A^-1 sum_i w_i x_i u_i.
#
# Wild cluster bootstrap (WCU, Rademacher weights): y* = X beta_hat + v u with
# v drawn at the level of one clustering dimension (country or year) and
# synchronized across every fit that shares those clusters. Each draw is
# studentized with the two-way (country, year) Cameron-Gelbach-Miller variance,
# so intervals are percentile-t and simultaneous bands use max-|t|. Draws are
# computed from country-year cell sums, so 9,999 draws take seconds. Following
# MacKinnon, Nielsen and Webb (2021, JBES), the bootstrap DGP clusters in one
# dimension while the statistic is two-way; both dimensions are run and the
# decision interval is the wider of the two.

# FWL pieces of a fixest fit. `clusters` names the columns of `data` (the data
# the model was fit on) used for clustering and cells.
fs_fwl <- function(fit, data, clusters = c("GID_0", "year")) {
  used <- fixest::obs(fit)
  dm <- fixest::demean(fit)
  X <- dm[, -1, drop = FALSE]
  X <- X[, names(stats::coef(fit)), drop = FALSE]
  w <- stats::weights(fit)
  if (is.null(w)) w <- rep(1, nrow(X))
  # beta and u are recomputed from the demeaned data so that X'Wu = 0 holds
  # exactly for the bootstrap; they must agree with fixest's own solution.
  A <- crossprod(X * w, X)
  beta <- solve(A, crossprod(X * w, dm[, 1]))[, 1]
  gap <- max(abs(beta - stats::coef(fit)) /
               pmax(abs(stats::coef(fit)), sqrt(diag(stats::vcov(fit)))))
  if (gap > 1e-3) {
    stop("FWL coefficients differ from fixest's (relative gap ",
         signif(gap, 3), "): tighten fixef.tol.")
  }
  out <- list(
    X = X, u = as.numeric(dm[, 1] - X %*% beta), w = w, coef = beta,
    obs = used, n = nrow(X), k = ncol(X), fixest_gap = gap
  )
  for (cl in clusters) out[[cl]] <- as.character(data[[cl]][used])
  out$A <- A
  out$Ainv <- solve(A)
  out
}

# Two-way CGM variance from the FWL pieces, with the G_min / (G_min - 1)
# small-cluster factor applied to each term (fixest's cluster.df = "min") and
# (n - 1) / (n - k) for the regressors (fixed effects are nested in the
# clusters). A non-positive-semidefinite result is repaired by truncating
# negative eigenvalues at zero; the repair is recorded as an attribute.
fs_cgm_vcov <- function(fwl, dims = c("GID_0", "year"), u = fwl$u,
                        adjust = TRUE) {
  scores <- fwl$X * (fwl$w * u)
  meat_one <- function(g) crossprod(rowsum(scores, g, reorder = FALSE))
  if (length(dims) == 1L) {
    meat <- meat_one(fwl[[dims]])
    G <- length(unique(fwl[[dims]]))
  } else {
    g1 <- fwl[[dims[1]]]
    g2 <- fwl[[dims[2]]]
    meat <- meat_one(g1) + meat_one(g2) - meat_one(paste(g1, g2))
    G <- min(length(unique(g1)), length(unique(g2)))
  }
  factor <- if (adjust) G / (G - 1) * (fwl$n - 1) / (fwl$n - fwl$k) else 1
  V <- factor * fwl$Ainv %*% meat %*% fwl$Ainv
  fs_psd_repair(V)
}

fs_psd_repair <- function(V) {
  dn <- dimnames(V)
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  repaired <- any(e$values < -1e-12 * max(abs(e$values)))
  if (repaired) {
    V <- e$vectors %*% diag(pmax(e$values, 0), nrow(V)) %*% t(e$vectors)
  }
  dimnames(V) <- dn
  attr(V, "psd_repaired") <- repaired
  V
}

# Deterministic seed from the design hash and a task name.
fs_seed <- function(design_hash, task) {
  hex <- substr(digest::digest(paste(design_hash, task), algo = "sha256"), 1, 7)
  strtoi(hex, 16L)
}

# Rademacher weights, one column per cluster (named), one row per draw.
fs_rademacher <- function(levels, reps, seed) {
  set.seed(seed)
  levels <- sort(unique(levels))
  V <- matrix(sample(c(-1, 1), length(levels) * reps, replace = TRUE),
              nrow = reps, dimnames = list(NULL, levels))
  V
}

# Country-year cell sums used by every draw.
fs_boot_prepare <- function(fwl, dims = c("GID_0", "year"), u = fwl$u) {
  cell <- paste(fwl[[dims[1]]], fwl[[dims[2]]], sep = "\r")
  cell_f <- factor(cell, levels = unique(cell))
  X <- fwl$X
  k <- ncol(X)
  S <- rowsum(X * (fwl$w * u), cell_f, reorder = FALSE)
  XX <- array(0, c(nrow(S), k, k))
  for (j in seq_len(k)) {
    XX[, j, ] <- rowsum(X * (fwl$w * X[, j]), cell_f, reorder = FALSE)
  }
  first <- match(levels(cell_f), cell)
  G <- min(length(unique(fwl[[dims[1]]])), length(unique(fwl[[dims[2]]])))
  list(
    S = S,
    XXm = matrix(XX, nrow(S) * k, k),
    k = k, n_cell = nrow(S),
    dim1 = fwl[[dims[1]]][first],
    dim2 = fwl[[dims[2]]][first],
    Ainv = fwl$Ainv,
    coef = fwl$coef,
    factor = G / (G - 1) * (fwl$n - 1) / (fwl$n - k),
    dims = dims
  )
}

# Draws of beta* - beta_hat and of the two-way variance of beta*, for weights
# `v` (reps x clusters, from fs_rademacher) applied in dimension `boot_dim`.
# Returns delta (reps x k) and Vflat (reps x k^2, column-major V*).
fs_boot_draws <- function(prep, v, boot_dim) {
  key <- if (boot_dim == prep$dims[1]) prep$dim1 else prep$dim2
  idx <- match(key, colnames(v))
  if (anyNA(idx)) stop("Bootstrap weights miss some ", boot_dim, " clusters.")
  reps <- nrow(v)
  k <- prep$k
  delta <- matrix(NA_real_, reps, k, dimnames = list(NULL, names(prep$coef)))
  Vflat <- matrix(NA_real_, reps, k * k)
  repaired <- logical(reps)
  d1 <- factor(prep$dim1)
  d2 <- factor(prep$dim2)
  for (b in seq_len(reps)) {
    vc <- v[b, idx]
    weighted <- prep$S * vc
    db <- as.numeric(prep$Ainv %*% colSums(weighted))
    cells <- weighted - matrix(prep$XXm %*% db, prep$n_cell, k)
    meat <- crossprod(rowsum(cells, d1, reorder = FALSE)) +
      crossprod(rowsum(cells, d2, reorder = FALSE)) -
      crossprod(cells)
    delta[b, ] <- db
    Vb <- prep$factor * (prep$Ainv %*% meat %*% prep$Ainv)
    # The same eigenvalue repair as the original two-way variance
    # (fs_cgm_vcov), so original and bootstrap studentization match.
    e <- eigen((Vb + t(Vb)) / 2, symmetric = TRUE)
    if (any(e$values < -1e-12 * max(abs(e$values)))) {
      Vb <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
      repaired[b] <- TRUE
    }
    Vflat[b, ] <- Vb
  }
  list(delta = delta, Vflat = Vflat, boot_dim = boot_dim, reps = reps,
       psd_repaired_share = mean(repaired))
}

# Estimates, two-way SEs and percentile-t intervals for the contrasts C (rows,
# with column names among the coefficients). `boots` is a named list of
# fs_boot_draws() results for the same fit. For each bootstrap dimension the
# symmetric interval is estimate +- q(|t*|) se; `sup_t` adds a simultaneous band
# from max_rows |t*|. The decision interval is the wider of the dimensions.
fs_boot_contrasts <- function(prep, boots, C, V_hat, levels = c(0.95, 0.90),
                              sup_t = FALSE) {
  if (is.null(dim(C))) C <- matrix(C, 1, dimnames = list(NULL, names(C)))
  full <- matrix(0, nrow(C), prep$k,
                 dimnames = list(rownames(C), names(prep$coef)))
  full[, colnames(C)] <- C
  estimate <- as.numeric(full %*% prep$coef)
  se <- sqrt(pmax(rowSums((full %*% V_hat) * full), 0))

  # vec(c c') for each contrast, to read var(c' beta*) off Vflat.
  outer_vec <- vapply(seq_len(nrow(full)), function(m) {
    as.numeric(outer(full[m, ], full[m, ]))
  }, numeric(prep$k^2))
  if (is.null(dim(outer_vec))) outer_vec <- matrix(outer_vec, ncol = 1)

  per_dim <- lapply(names(boots), function(dim) {
    bt <- boots[[dim]]
    theta <- bt$delta %*% t(full)
    se_star <- sqrt(pmax(bt$Vflat %*% outer_vec, 0))
    tstar <- abs(theta) / se_star
    tstar[!is.finite(tstar)] <- NA
    q <- sapply(levels, function(l) {
      apply(tstar, 2, stats::quantile, probs = l, na.rm = TRUE)
    })
    q <- matrix(q, nrow(full))
    colnames(q) <- paste0("q", round(100 * levels))
    out <- tibble::as_tibble(q)
    out$boot_dim <- dim
    out$draws_used <- colSums(is.finite(tstar))
    out$contrast <- rownames(full)
    if (sup_t) {
      sup <- apply(tstar, 1, max, na.rm = TRUE)
      out$q_sup95 <- stats::quantile(sup, 0.95, na.rm = TRUE)
    }
    out
  })
  per_dim <- dplyr::bind_rows(per_dim)
  crit <- per_dim %>%
    dplyr::group_by(contrast) %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("q"), max),
                     .groups = "drop")
  crit <- crit[match(rownames(full), crit$contrast), ]
  res <- tibble::tibble(
    contrast = rownames(full),
    estimate = estimate,
    se = se
  )
  for (l in levels) {
    qn <- paste0("q", round(100 * l))
    res[[paste0("ci", round(100 * l), "_low")]] <- estimate - crit[[qn]] * se
    res[[paste0("ci", round(100 * l), "_high")]] <- estimate + crit[[qn]] * se
  }
  if (sup_t) {
    res$band95_low <- estimate - crit$q_sup95 * se
    res$band95_high <- estimate + crit$q_sup95 * se
  }
  attr(res, "per_dimension") <- per_dim
  res
}

# Analytic intervals for the contrasts with a given variance (t with
# G_min - 1 degrees of freedom by default).
fs_analytic_contrasts <- function(coef, V, C, df = Inf,
                                  levels = c(0.95, 0.90), label = NULL) {
  if (is.null(dim(C))) C <- matrix(C, 1, dimnames = list(NULL, names(C)))
  full <- matrix(0, nrow(C), length(coef),
                 dimnames = list(rownames(C), names(coef)))
  full[, colnames(C)] <- C
  est <- as.numeric(full %*% coef)
  se <- sqrt(pmax(rowSums((full %*% V) * full), 0))
  out <- tibble::tibble(contrast = rownames(full), estimate = est, se = se)
  for (l in levels) {
    z <- stats::qt(1 - (1 - l) / 2, df)
    out[[paste0("ci", round(100 * l), "_low")]] <- est - z * se
    out[[paste0("ci", round(100 * l), "_high")]] <- est + z * se
  }
  if (!is.null(label)) out$inference <- label
  out
}

# Conley spatial HAC with a temporal Bartlett kernel (Hsiang 2010): uniform
# spatial kernel within `cutoff_km` across units in the same year, plus
# Bartlett-weighted autocovariances up to `lag_years` within each unit.
# `coords` has GID_1, longitude, latitude for every unit in the fit.
fs_conley_hac <- function(fwl, unit, year, coords, cutoff_km, lag_years) {
  scores <- fwl$X * (fwl$w * fwl$u)
  k <- ncol(scores)
  units <- sort(unique(unit))
  coords <- coords[match(units, coords$GID_1), ]
  if (anyNA(coords$longitude)) stop("Coordinates missing for some units.")
  pairs <- region_pair_distances(coords, max_km = cutoff_km)
  n_u <- length(units)
  K <- Matrix::sparseMatrix(
    i = c(match(pairs$region_i, units), match(pairs$region_j, units),
          seq_len(n_u)),
    j = c(match(pairs$region_j, units), match(pairs$region_i, units),
          seq_len(n_u)),
    x = 1, dims = c(n_u, n_u)
  )
  years <- sort(unique(year))
  ui <- match(unit, units)
  yi <- match(year, years)
  meat <- matrix(0, k, k)
  for (t in seq_along(years)) {
    rows <- which(yi == t)
    St <- matrix(0, n_u, k)
    St[ui[rows], ] <- scores[rows, , drop = FALSE]
    meat <- meat + as.matrix(crossprod(St, K %*% St))
  }
  if (lag_years > 0) {
    grid <- array(0, c(n_u, length(years), k))
    for (j in seq_len(k)) grid[cbind(ui, yi, j)] <- scores[, j]
    for (l in seq_len(lag_years)) {
      if (l >= length(years)) break
      w <- 1 - l / (lag_years + 1)
      now <- matrix(grid[, (l + 1):length(years), , drop = FALSE], ncol = k)
      past <- matrix(grid[, 1:(length(years) - l), , drop = FALSE], ncol = k)
      cross <- crossprod(now, past)
      meat <- meat + w * (cross + t(cross))
    }
  }
  V <- fwl$Ainv %*% meat %*% fwl$Ainv
  dimnames(V) <- list(names(fwl$coef), names(fwl$coef))
  fs_psd_repair(V)
}

# Section 7 decision rules for one contrast. Negative estimates mean the fast
# path is worse. `supported` is the identification gate; `power_ok` is the
# Stage 3 power check at the SESOI.
fs_classify <- function(ci95_low, ci95_high, ci90_low, ci90_high, sesoi,
                        sesoi_confirmed = FALSE, supported = TRUE,
                        power_ok = TRUE) {
  n <- length(ci95_low)
  direction <- dplyr::case_when(
    !supported ~ "not identified for this path",
    ci95_high < 0 ~ "fast worse",
    ci95_low > 0 ~ "fast less harmful",
    TRUE ~ "direction unresolved"
  )
  materiality <- dplyr::case_when(
    !supported ~ "not identified for this path",
    ci95_high < -sesoi ~ "materially fast worse",
    ci95_low > sesoi ~ "materially fast less harmful",
    ci90_low >= -sesoi & ci90_high <= sesoi ~ "practically equivalent",
    TRUE ~ "magnitude unresolved"
  )
  # Low power is imprecision: a null or failed-equivalence result stays
  # unresolved, while a clear directional interval or a valid equivalence
  # interval stands.
  materiality <- ifelse(
    !power_ok & materiality %in% c("magnitude unresolved"),
    "magnitude unresolved (low power)", materiality
  )
  not_materially_worse <- supported & ci90_low > -sesoi
  claim <- if (sesoi_confirmed) materiality else ifelse(
    materiality %in% c("practically equivalent", "materially fast worse",
                       "materially fast less harmful"),
    paste0("provisional: ", materiality, " (SESOI not affirmed)"),
    materiality
  )
  tibble::tibble(
    direction = rep_len(direction, n),
    materiality = rep_len(claim, n),
    not_materially_worse_one_sided = rep_len(not_materially_worse, n)
  )
}

# All inference for one fit: two-way CGM variance (own and fixest's) and the
# wild bootstrap in both dimensions, with the synchronized Rademacher weights
# `boot_weights` (list(GID_0 = ..., year = ...) from fs_rademacher()). `u`
# replaces the residuals in the bootstrap DGP (the Stage 3 residual-calibrated
# power pass uses the level-model residuals).
fs_infer <- function(fit, data, boot_weights, u = NULL) {
  fwl <- fs_fwl(fit, data)
  prep <- fs_boot_prepare(fwl, u = if (is.null(u)) fwl$u else u)
  boots <- lapply(names(boot_weights), function(dim) {
    fs_boot_draws(prep, boot_weights[[dim]], dim)
  })
  names(boots) <- names(boot_weights)
  G <- min(length(unique(fwl$GID_0)), length(unique(fwl$year)))
  list(fit = fit, fwl = fwl, V_hat = fs_cgm_vcov(fwl), prep = prep,
       boots = boots, df = G - 1,
       V_fixest = stats::vcov(fit, vcov = ~ GID_0 + year))
}

# Contrast table: estimate, two-way SE, wild-bootstrap percentile-t intervals
# (the decision intervals) and fixest's analytic two-way intervals.
fs_contrast_table <- function(inf, C, sup_t = FALSE) {
  boot <- fs_boot_contrasts(inf$prep, inf$boots, C, inf$V_hat, sup_t = sup_t)
  analytic <- fs_analytic_contrasts(inf$fit$coefficients, inf$V_fixest, C,
                                    df = inf$df)
  boot$analytic_se <- analytic$se
  boot$analytic_ci95_low <- analytic$ci95_low
  boot$analytic_ci95_high <- analytic$ci95_high
  boot$bootstrap_draws <- min(vapply(inf$boots, `[[`, numeric(1), "reps"))
  boot
}

# Estimate, SE and bootstrap |t*| draws (one column per dimension) of a single
# contrast, for simultaneous bands across separately estimated fits that share
# the bootstrap weights (local projections).
fs_boot_tstats <- function(inf, crow) {
  full <- stats::setNames(numeric(inf$prep$k), names(inf$prep$coef))
  full[names(crow)] <- crow
  cc <- as.numeric(outer(full, full))
  est <- sum(full * inf$prep$coef)
  se <- sqrt(max(sum(cc * as.numeric(inf$V_hat)), 0))
  tstar <- vapply(inf$boots, function(bt) {
    abs(as.numeric(bt$delta %*% full)) / sqrt(pmax(bt$Vflat %*% cc, 0))
  }, numeric(nrow(inf$boots[[1]]$delta)))
  tstar[!is.finite(tstar)] <- NA
  list(estimate = est, se = se, tstar = tstar)
}

# Wald test of R beta = r with a given variance (chi-square and F forms).
fs_wald_linear <- function(coef, V, R, r = rep(0, nrow(R)), df2 = Inf) {
  full <- matrix(0, nrow(R), length(coef), dimnames = list(NULL, names(coef)))
  full[, colnames(R)] <- R
  d <- full %*% coef - r
  stat <- as.numeric(t(d) %*% solve(full %*% V %*% t(full)) %*% d)
  q <- nrow(R)
  tibble::tibble(chi2 = stat, df1 = q, F_stat = stat / q, df2 = df2,
                 p_value = stats::pf(stat / q, q, df2, lower.tail = FALSE))
}
