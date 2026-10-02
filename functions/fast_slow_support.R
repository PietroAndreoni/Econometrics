# Outcome-blind support audit, fallback contrast and design-range power for the
# fast-versus-slow warming design (Stage 2 of
# FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md). Nothing here reads an outcome.

# Scenario support. `obs` holds the analytic cells (GID_0, GID_1, TM, rate, RR);
# `scen` the scenario-region-years (GID_1, path, s, TM, rate, RR) with GID_0.
# A scenario-region-year is supported when
#   (1) its rate lies inside the pooled `rate_q` quantile range of `obs$rate`;
#   (2) its temperature lies inside the region's own observed [min, max];
#   (3) its standardized distance to the k-th nearest observed cell in
#       (TM, rate, RR) space is at most the `dist_q` quantile of the
#       leave-one-out k-th-neighbour distances among observed cells, and
#       those k neighbours span at least `min_regions` regions and
#       `min_countries` countries.
fs_support_scenario <- function(obs, scen, k = 50L, dist_q = 0.99,
                                rate_q = c(0.01, 0.99), min_regions = 10L,
                                min_countries = 3L, loo = NULL) {
  vars <- c("TM", "rate", "RR")
  scale <- vapply(vars, function(v) stats::sd(obs[[v]]), numeric(1))
  Z <- sweep(as.matrix(obs[, vars]), 2, scale, "/")
  Zq <- sweep(as.matrix(scen[, vars]), 2, scale, "/")
  if (is.null(loo)) {
    # Column 1 is the cell itself.
    loo <- RANN::nn2(Z, Z, k = k + 1L)$nn.dists[, k + 1L]
  }
  threshold <- stats::quantile(loo, dist_q)
  nn <- RANN::nn2(Z, Zq, k = k)
  n_regions <- apply(nn$nn.idx, 1, function(i) length(unique(obs$GID_1[i])))
  n_countries <- apply(nn$nn.idx, 1, function(i) length(unique(obs$GID_0[i])))

  rate_bounds <- stats::quantile(obs$rate, rate_q)
  local <- obs %>%
    dplyr::group_by(GID_1) %>%
    dplyr::summarise(T_min = min(TM), T_max = max(TM), .groups = "drop")
  out <- scen %>%
    dplyr::left_join(local, by = "GID_1") %>%
    dplyr::mutate(
      rate_ok = rate >= rate_bounds[[1]] & rate <= rate_bounds[[2]],
      temp_ok = !is.na(T_min) & TM >= T_min & TM <= T_max,
      knn_distance = nn$nn.dists[, k],
      knn_regions = n_regions,
      knn_countries = n_countries,
      knn_ok = knn_distance <= threshold & knn_regions >= min_regions &
        knn_countries >= min_countries,
      supported = rate_ok & temp_ok & knn_ok
    )
  attr(out, "rate_bounds") <- rate_bounds
  attr(out, "distance_threshold") <- threshold
  attr(out, "loo") <- loo
  out
}

# Weighted share of each path-year that is supported, and the scenario verdict:
# every evaluated year of both paths must reach `min_share`.
fs_support_summary <- function(support, weights, horizon, min_share = 0.95) {
  by_year <- support %>%
    dplyr::filter(s <= horizon) %>%
    dplyr::inner_join(weights, by = "GID_1") %>%
    dplyr::group_by(path, s) %>%
    dplyr::summarise(
      supported_share = sum(weight * supported) / sum(weight),
      rate_ok_share = sum(weight * rate_ok) / sum(weight),
      temp_ok_share = sum(weight * temp_ok) / sum(weight),
      knn_ok_share = sum(weight * knn_ok) / sum(weight),
      .groups = "drop"
    )
  list(
    by_year = by_year,
    min_share = min(by_year$supported_share),
    passes = min(by_year$supported_share) >= min_share
  )
}

# Region-years, countries and decades in a positive-rate tail.
fs_tail_counts <- function(obs, rate_col, threshold) {
  tail <- obs[obs[[rate_col]] >= threshold, ]
  tibble::tibble(
    threshold = threshold,
    region_years = nrow(tail),
    regions = dplyr::n_distinct(tail$GID_1),
    countries = dplyr::n_distinct(tail$GID_0),
    decades = dplyr::n_distinct(10L * (tail$year %/% 10L))
  )
}

# Deterministic fallback pair (section 1): every (M, D_F, D_S) with D_F < D_S
# that passes support, ranked by standardized distance of its ramp rates
# M / D_F and M / D_S to the target quantiles of positive observed slopes; ties
# go to larger M, then shorter D_F, then shorter D_S.
fs_fallback_table <- function(candidates, positive_rates, target_q) {
  targets <- stats::quantile(positive_rates, target_q)
  scale <- stats::sd(positive_rates)
  candidates %>%
    dplyr::mutate(
      fast_rate = total_warming / fast_years,
      slow_rate = total_warming / slow_years,
      target_fast = targets[[1]],
      target_slow = targets[[2]],
      distance = sqrt(((fast_rate - target_fast) / scale)^2 +
                        ((slow_rate - target_slow) / scale)^2)
    ) %>%
    dplyr::arrange(dplyr::desc(passes), distance, dplyr::desc(total_warming),
                   fast_years, slow_years) %>%
    dplyr::mutate(rank = ifelse(passes, cumsum(passes), NA_integer_))
}

# Stationary AR(1) country-year shocks (countries x years x sims), marginal SD
# `sd`.
fs_country_ar1 <- function(n_country, n_year, n_sim, sd, rho) {
  out <- array(0, c(n_country, n_year, n_sim))
  innov_sd <- sd * sqrt(1 - rho^2)
  out[, 1, ] <- stats::rnorm(n_country * n_sim, sd = sd)
  if (n_year > 1) for (t in 2:n_year) {
    out[, t, ] <- rho * out[, t - 1, ] +
      stats::rnorm(n_country * n_sim, sd = innov_sd)
  }
  out
}

# Two-way CGM variance of the contrast with influence vector `a` (a_i =
# w_i x_i' A^-1 c), for many residual columns at once: returns one variance per
# column. `cells`, `cell_d1`, `cell_d2` index the country-year cells.
fs_contrast_var_cols <- function(a, resid, cell, cell_d1, cell_d2, factor) {
  C <- rowsum(a * resid, cell, reorder = FALSE)
  factor * (colSums(rowsum(C, cell_d1, reorder = FALSE)^2) +
              colSums(rowsum(C, cell_d2, reorder = FALSE)^2) -
              colSums(C^2))
}

# Design-range power (Stage 2) for the primary rate model. `X` holds the
# demeaned regressors of the primary sample (no outcome is used), `fe` the
# fixed-effect identifiers (GID_1, GID_0, year, year_c), `contrasts` the path
# contrast rows (named "speed" and "total"), `lag_signal` the demeaned and raw
# misspecification signal sum_k w_k T_{t-k}, `lag_total` its true total path
# gap per unit scale. Returns one row per noise cell, DGP and planted gap.
fs_power_design <- function(X, fe, contrasts, lag_signal_demeaned, lag_total,
                            grid, targets, sesoi, n_sim, seed, chunk = 250L) {
  Ainv <- solve(crossprod(X))
  k <- ncol(X)
  n <- nrow(X)
  cell <- factor(paste(fe$GID_0, fe$year))
  cell_first <- match(levels(cell), paste(fe$GID_0, fe$year))
  cell_d1 <- fe$GID_0[cell_first]
  cell_d2 <- fe$year[cell_first]
  G <- min(dplyr::n_distinct(fe$GID_0), dplyr::n_distinct(fe$year))
  factor <- G / (G - 1) * (n - 1) / (n - k)
  df <- G - 1
  crit95 <- stats::qt(0.975, df)
  crit90 <- stats::qt(0.95, df)

  A_vec <- X %*% Ainv %*% t(contrasts)          # n x contrasts
  est_proj <- as.numeric(contrasts %*% Ainv %*% crossprod(X, lag_signal_demeaned))
  names(est_proj) <- rownames(contrasts)
  resid_sig <- lag_signal_demeaned - X %*% (Ainv %*% crossprod(X, lag_signal_demeaned))
  speed_row <- "speed"
  total_row <- "total"
  c_rpos2 <- contrasts[speed_row, "r_pos2"]

  countries <- sort(unique(fe$GID_0))
  years <- sort(unique(fe$year))
  ci <- match(fe$GID_0, countries)
  yi <- match(fe$year, years)
  fe_df <- data.frame(GID_1 = fe$GID_1, year = fe$year, GID_0 = fe$GID_0)
  slope <- data.frame(year_c = fe$year_c)

  results <- list()
  for (g in seq_len(nrow(grid))) {
    set.seed(seed + g)
    est_noise <- list()
    var_noise <- list()
    sums_sig <- list()
    remaining <- n_sim
    while (remaining > 0) {
      m <- min(chunk, remaining)
      remaining <- remaining - m
      shocks <- fs_country_ar1(length(countries), length(years), m,
                               grid$country_sd[g], grid$country_ar1[g])
      E <- matrix(stats::rnorm(n * m, sd = grid$idiosyncratic_sd[g]), n, m)
      E <- E + matrix(shocks[cbind(rep(ci, m), rep(yi, m),
                                   rep(seq_len(m), each = n))], n, m)
      year_shock <- matrix(stats::rnorm(length(years) * m,
                                        sd = grid$year_sd[g]),
                           length(years), m)
      E <- E + year_shock[yi, , drop = FALSE]
      Et <- fixest::demean(E, f = fe_df, slope.vars = slope,
                           slope.flag = c(0L, 0L, -1L), notes = FALSE)
      delta <- Ainv %*% crossprod(X, Et)
      R <- Et - X %*% delta
      est_noise[[length(est_noise) + 1]] <- contrasts %*% delta
      # Cell sums of the influence-weighted residuals, kept per contrast so the
      # misspecification component can be added for any scale.
      sums <- lapply(rownames(contrasts), function(r) {
        rowsum(A_vec[, r] * R, cell, reorder = FALSE)
      })
      names(sums) <- rownames(contrasts)
      sums_sig[[length(sums_sig) + 1]] <- sums
    }
    est_noise <- do.call(cbind, est_noise)
    cell_sums <- lapply(rownames(contrasts), function(r) {
      do.call(cbind, lapply(sums_sig, `[[`, r))
    })
    names(cell_sums) <- rownames(contrasts)
    sig_cells <- lapply(rownames(contrasts), function(r) {
      as.numeric(rowsum(A_vec[, r] * resid_sig, cell, reorder = FALSE))
    })
    names(sig_cells) <- rownames(contrasts)

    variance <- function(Ccells) {
      factor * (colSums(rowsum(Ccells, cell_d1, reorder = FALSE)^2) +
                  colSums(rowsum(Ccells, cell_d2, reorder = FALSE)^2) -
                  colSums(Ccells^2))
    }

    for (dgp in c("planted_speed_curvature", "misspecified_lags")) {
      for (target in targets) {
        for (r in rownames(contrasts)) {
          if (dgp == "planted_speed_curvature") {
            truth <- target            # speed gap = total gap
            est <- target + est_noise[r, ]
            v <- variance(cell_sums[[r]])
          } else {
            kappa <- target / lag_total
            truth <- if (r == total_row) target else NA_real_
            est <- kappa * est_proj[[r]] + est_noise[r, ]
            v <- variance(cell_sums[[r]] + kappa * sig_cells[[r]])
          }
          se <- sqrt(pmax(v, 0))
          lo95 <- est - crit95 * se
          hi95 <- est + crit95 * se
          lo90 <- est - crit90 * se
          hi90 <- est + crit90 * se
          cls <- fs_classify(lo95, hi95, lo90, hi90, sesoi,
                             sesoi_confirmed = TRUE)
          correct <- if (is.na(truth) || abs(abs(truth) - sesoi) < 1e-12) {
            NA_real_
          } else if (abs(truth) < sesoi) {
            mean(cls$materiality == "practically equivalent")
          } else {
            mean(cls$materiality == "materially fast worse")
          }
          results[[length(results) + 1]] <- tibble::tibble(
            grid[g, ], dgp = dgp, planted_gap = target, contrast = r,
            truth = truth,
            mean_estimate = mean(est),
            bias = if (is.na(truth)) NA_real_ else mean(est) - truth,
            coverage95 = if (is.na(truth)) NA_real_ else
              mean(lo95 <= truth & hi95 >= truth),
            p_fast_worse = mean(cls$direction == "fast worse"),
            p_fast_less_harmful = mean(cls$direction == "fast less harmful"),
            p_materially_worse = mean(cls$materiality == "materially fast worse"),
            p_equivalent = mean(cls$materiality == "practically equivalent"),
            p_correct_sesoi_class = correct,
            mean_se = mean(se)
          )
        }
      }
    }
    message("power cell ", g, "/", nrow(grid))
  }
  attr_out <- dplyr::bind_rows(results)
  attr(attr_out, "rpos2_contrast") <- c_rpos2
  attr_out
}

# Permutation null for the incremental speed gap (Stage 12). Only the rate
# block is reassigned; the level and precipitation columns and the outcome stay
# fixed, and the model is refitted on its FWL form. `d` is the estimation data
# (GID_0, GID_1, year, year_c, zone, decade), `climate_terms` the complete
# climate panel. Donor schemes:
#   "region_within_zone_decade"  region labels permuted within climate zone x
#                                decade strata (a region's decade block moves
#                                as a whole); breaks spatial correlation of
#                                the exposure, so the null is too narrow when
#                                errors are spatially correlated;
#   "country_time_shift"         every region of a country takes its own rate
#                                history lagged by one random k in `shifts`
#                                years; preserves within-country spatial and
#                                serial structure (dependence-consistent).
fs_rate_permutation <- function(d, climate_terms, X_fixed, y_dm, speed_row,
                                scheme = c("region_within_zone_decade",
                                           "country_time_shift"),
                                n_perm, seed, shifts = 6:17) {
  scheme <- match.arg(scheme)
  rate_terms <- names(speed_row)
  fe_df <- data.frame(GID_1 = d$GID_1, year = d$year, GID_0 = d$GID_0)
  strata <- paste(ifelse(is.na(d$zone), "none", d$zone), d$decade)
  ckey <- paste(climate_terms$GID_1, climate_terms$year)
  countries <- unique(d$GID_0)
  set.seed(seed)
  t(vapply(seq_len(n_perm), function(b) {
    if (scheme == "region_within_zone_decade") {
      donor <- d$GID_1
      for (st in unique(strata)) {
        rows <- which(strata == st)
        ids <- unique(d$GID_1[rows])
        map <- stats::setNames(ids[sample.int(length(ids))], ids)
        donor[rows] <- map[d$GID_1[rows]]
      }
      idx <- match(paste(donor, d$year), ckey)
    } else {
      k <- stats::setNames(sample(shifts, length(countries), replace = TRUE),
                           countries)
      idx <- match(paste(d$GID_1, d$year - k[d$GID_0]), ckey)
    }
    R <- as.matrix(climate_terms[idx, rate_terms])
    if (anyNA(R)) stop("Permuted rate histories contain NA.")
    R_dm <- fixest::demean(R, f = fe_df,
                           slope.vars = data.frame(year_c = d$year_c),
                           slope.flag = c(0L, 0L, -1L), notes = FALSE)
    X <- cbind(X_fixed, R_dm)
    beta <- solve(crossprod(X), crossprod(X, y_dm))[, 1]
    c(speed = sum(speed_row * beta[rate_terms]),
      r_pos2 = unname(beta["r_pos2"]))
  }, numeric(2)))
}
