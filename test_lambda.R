require(fixest)

# Centering used in the deviation terms:
#   "window" -> rolling mean_TM / mean_RR, LAGGED one year (window-specific)
#   "all"    -> mean_TM_all / mean_RR_all (fixed full-period, window-invariant)
#
# mean_TM / mean_RR come out of prepare_climate_data.R as right-aligned rolling
# means, i.e. the window ENDING at t, which contains TM_t / RR_t itself. Using
# them raw would net part of the current shock out of its own deviation term, so
# the centring is the window ending at t-1 -- the same lag(mean_.) that
# prepare_climate_data.R applies when it builds dev_TM / dev_RR.
moments_grid <- c("window", "all")

moment_cols <- list(
  window = c(mean_TM = "mean_TM_lag", mean_RR = "mean_RR_lag"),
  all    = c(mean_TM = "mean_TM_all", mean_RR = "mean_RR_all")
)

# The scaling SD is always the fixed full-period one, for every specification
sd_cols <- c(sd_TM = "sd_TM_pre", sd_RR = "sd_RR_pre")

# Every column any specification needs, used to define the common sample
moment_source_cols <- unique(c(unlist(moment_cols, use.names = FALSE), sd_cols))


# First year of the ESTIMATION sample. Climate before it is still carried (see
# CLIM_HISTORY_YEARS), only the outcome is withheld.
ECON_YEAR_MIN <- 1990

# Years of climate kept BEFORE ECON_YEAR_MIN. Those rows have no outcome, so
# feols drops them from every fit, but fixest still builds panel lags through
# them: without the carriers, the one-year centring lag and l(., 0:5) would eat
# the first years of the sample instead (checked: 440 vs 240 obs on a 40x16 toy).
# Six covers the deepest lag used anywhere below, plus one for the centring.
CLIM_HISTORY_YEARS <- 6

# A row is estimable when it has an outcome; everything else is a lag carrier
est_row <- function(d) !is.na(d$dlgrp_pc_usd)

# Value of x in year t-k, matched on the year itself rather than on row order,
# so gaps in a region's panel produce NA instead of silently shifting the lag.
lag_by_year <- function(x, year, k = 1) x[match(year - k, year)]


build_dat <- function(w_rr, w_tm) {

  # Econ is restricted to the estimation window here; the climate panel below is
  # NOT, so the pre-ECON_YEAR_MIN years survive the join as outcome-less rows.
  econ <- econ_data %>%
    filter(econ_source == "DOSE_V2_11" &
             !is.na(dlgrp_pc_usd) &
             year >= ECON_YEAR_MIN)

  data_rr %>%
    filter(window_rr==w_rr &
             weight=="area" &
             climate_source=="era5" &
             (weight_year=="2015"|weight_year=="un") &
             !is.na(RR)) %>%
    select(-window_rr,-weight,-weight_year) %>%
    inner_join( data_tm %>%
                  filter(window_tm==w_tm &
                           weight=="area" &
                           climate_source=="era5" &
                           (weight_year=="2015"|weight_year=="un") &
                           !is.na(TM)) %>%
                  select(-window_tm,-weight,-weight_year) ) %>%
    filter(year >= ECON_YEAR_MIN - CLIM_HISTORY_YEARS) %>%
    left_join(econ) %>%
    group_by(GID_1) %>%
    # Regions with nothing to estimate are pure ballast, and climate running past
    # a region's last usable year carries no lag anyone reads
    filter(any(!is.na(dlgrp_pc_usd))) %>%
    filter(year <= max(year[!is.na(dlgrp_pc_usd)])) %>%
    mutate(log_gdp_av=log(mean(grp_pc_usd,na.rm=TRUE))) %>% 
    arrange(year, .by_group = TRUE) %>%
    # Centring of the deviation terms: the rolling window ending at t-1, so the
    # current year never enters its own reference mean. Built here, before the
    # common-sample restriction, so the warm-up year that loses its lag is
    # dropped by the completeness check below like any other missing moment.
    mutate(
      mean_TM_lag = lag_by_year(mean_TM, year),
      mean_RR_lag = lag_by_year(mean_RR, year)
    )
}


# Deviation terms (one per centering) and the common relative SDs. Must be
# called AFTER the sample restriction so the reference SDs describe the
# estimation sample.
add_moment_terms <- function(d) {

  # Reference SDs only improve coefficient interpretation. Taken over the
  # estimable rows only, so the lag carriers cannot shift the scaling.
  e <- est_row(d)
  sT_ref <- median(d[[sd_cols["sd_TM"]]][e], na.rm = TRUE)
  sR_ref <- median(d[[sd_cols["sd_RR"]]][e], na.rm = TRUE)

  d$sd_TM_rel <- d[[sd_cols["sd_TM"]]] / sT_ref
  d$sd_RR_rel <- d[[sd_cols["sd_RR"]]] / sR_ref

  for (mo in moments_grid) {
    cols <- moment_cols[[mo]]

    # Absolute deviation, so the exponent can take non-integer values
    d[[paste0("TM_absdev_", mo)]] <- abs(d$TM - d[[cols["mean_TM"]]])
    d[[paste0("RR_absdev_", mo)]] <- abs(d$RR - d[[cols["mean_RR"]]])
  }

  attr(d, "refs") <- c(sT_ref = sT_ref, sR_ref = sR_ref)

  d
}


# Regressor construction, split out from the fit so cross-validation can build
# the columns once per specification and then refit on subsets of the rows
dev_lambda_cols <- function(d, lambda_T, lambda_R, moments = "window",
                            p_T = 2, p_R = 2) {
  
  d$TM_dev_abs_lambda <-
    d[[paste0("TM_absdev_", moments)]]^p_T 
  
  d$RR_dev_abs_lambda <-
    d[[paste0("RR_absdev_", moments)]]^p_R 
  
  d$TM_dev_lambda <-
    d[[paste0("TM_absdev_", moments)]]^p_T / d$sd_TM_rel^p_T

  d$RR_dev_lambda <-
    d[[paste0("RR_absdev_", moments)]]^p_R / d$sd_RR_rel^p_R

  d
}

fit_on <- function(d) {
  feols(
    dlgrp_pc_usd ~
      (TM + TM_2 +
      RR + RR_2 +
      TM_dev_lambda +
      RR_dev_lambda) |
      year + GID_1 + GID_0[year,year^2] + GID_1[year,year^2],
    cluster = ~GID_1,
    data = d
  )
}

fit_lambdas <- function(lambda_T, lambda_R, d, moments = "window",
                        p_T = 2, p_R = 2) {

  fit_on(dev_lambda_cols(d, lambda_T, lambda_R, moments, p_T, p_R))
}

window_rr_grid <- c(5,10,20,30)
window_tm_grid <- c(5,10,20,30)

# Exponent on |X - Xmean| in the deviation terms
p_T_grid <- seq(0.5, 3, by = 0.5)
p_R_grid <- seq(0.5, 3, by = 0.5)

# The window only matters for moments == "window"; the fixed-baseline
# specifications are estimated once, on the reference dataset.
grid <- rbind(
  expand.grid(
    p_T = p_T_grid,
    p_R = p_R_grid,
    window_rr = window_rr_grid,
    window_tm = window_tm_grid,
    moments = "window",
    stringsAsFactors = FALSE
  ),
  expand.grid(
    p_T = p_T_grid,
    p_R = p_R_grid,
    window_rr = NA_real_,
    window_tm = NA_real_,
    moments = setdiff(moments_grid, "window"),
    stringsAsFactors = FALSE
  )
)

# Build each climate-window dataset once and reuse it across the lambda grid
dat_cache <- list()
for (w_rr in window_rr_grid) {
  for (w_tm in window_tm_grid) {
    dat_cache[[paste(w_rr, w_tm, sep = "_")]] <- build_dat(w_rr, w_tm)
  }
}

# Restrict every dataset to the common estimation sample, i.e. the rows where
# all moment sets are defined for all windows. This is the sample of the widest
# rolling window (the smallest one), and makes SSR/logLik comparable across the
# whole grid.
row_key <- function(d) paste(d$GID_1, d$year, sep = "_")

common_keys <- Reduce(
  intersect,
  lapply(dat_cache, function(d) {
    ok <- est_row(d) & complete.cases(
      as.data.frame(ungroup(d))[, moment_source_cols, drop = FALSE]
    )
    row_key(d)[ok]
  })
)

# Rows outside the common sample are not deleted -- they are demoted to lag
# carriers by blanking the outcome. That keeps the climate history available to
# l() while making it impossible for any cell to estimate on a row the other
# cells cannot, which is what makes SSR comparable across the grid.
dat_cache <- lapply(dat_cache, function(d) {
  d$dlgrp_pc_usd[!(row_key(d) %in% common_keys)] <- NA
  add_moment_terms(d[d$GID_1 %in% unique(d$GID_1[est_row(d)]), , drop = FALSE])
})

cat("Common sample:", length(common_keys), "region-year observations",
    "from", ECON_YEAR_MIN, "onwards;",
    "carrier rows kept for lags:",
    paste(sapply(dat_cache, function(d) sum(!est_row(d))), collapse = ", "), "\n")

# Rows are now identical across windows (only the rolling moments differ), so
# any cached dataset serves the window-invariant specifications
ref_key <- paste(window_rr_grid[1], window_tm_grid[1], sep = "_")

# Which cached dataset a grid cell is estimated on
cell_key <- function(w_rr, w_tm) {
  if (is.na(w_rr)) ref_key else paste(w_rr, w_tm, sep = "_")
}


## ---- scaffolding for the fit indicators ------------------------------------
## Every grid cell has the same n and the same K (six regressors, identical
## fixed effects), so RMSE, R2, adjusted R2, AIC and BIC are all monotone
## transformations of SSR: they would reproduce the SSR ranking exactly. The
## indicators below are the ones that can actually RE-RANK the grid, and they
## are only comparable if every cell is scored on the SAME rows -- hence the
## tail split, the trimmed sample, the CV folds and the reference SSR are all
## built once here, on the reference dataset, and reused for every cell.

ref_dat <- dat_cache[[ref_key]]

# Estimable rows of the reference dataset. Every quantile and every fold below
# is defined on these only: the lag carriers have no outcome and must never
# influence a threshold, a fold or a score.
ref_est <- est_row(ref_dat)

# The reuse above is only valid if the cached datasets are row-aligned
stopifnot(all(vapply(
  dat_cache,
  function(d) identical(row_key(d), row_key(ref_dat)),
  logical(1)
)))

# Tail = the decile of the sample with the largest climate anomaly, measured on
# the fixed full-period centring and the fixed SD. Spec-independent on purpose:
# the same rows are the tail in every cell, so RMSE_tail compares like with like
# instead of scoring each cell on a subsample of its own choosing.
clim_anom <- pmax(
  abs(ref_dat$TM - ref_dat$mean_TM_all) / ref_dat[[sd_cols["sd_TM"]]],
  abs(ref_dat$RR - ref_dat$mean_RR_all) / ref_dat[[sd_cols["sd_RR"]]]
)
is_tail <- ref_est & clim_anom >= quantile(clim_anom[ref_est], 0.9, na.rm = TRUE)

# Trimmed sample for the outlier-robust loss. Trimming on the OUTCOME (not on
# each cell's own residuals) keeps the trimmed rows fixed across the grid, so a
# cell cannot win SSR_trim by relabelling which observations count as outliers.
y_lim <- quantile(ref_dat$dlgrp_pc_usd[ref_est], c(0.01, 0.99), na.rm = TRUE)
is_trim <- ref_est &
  ref_dat$dlgrp_pc_usd >= y_lim[1] & ref_dat$dlgrp_pc_usd <= y_lim[2]

# 5-fold CV folds drawn WITHIN region, so every region keeps its intercept and
# its quadratic trend identified in the training fold and its held-out rows stay
# predictable. Leaving out whole regions (or whole years) would remove the very
# fixed effect the prediction needs. Drawn once, so all cells see identical
# folds and the CV comparison is paired. Fold 0 = lag carriers: never held out,
# so they stay in every training fold and keep the lags resolvable there too.
cv_k <- 5
set.seed(20260825)
cv_fold <- integer(nrow(ref_dat))
cv_fold[ref_est] <- ave(
  seq_len(sum(ref_est)), ref_dat$GID_1[ref_est],
  FUN = function(i) sample(rep_len(seq_len(cv_k), length(i)))
)

# Yardstick for the win rate: the plain quadratic, with no variability term at
# all. A cell "wins" a region when it fits that region's years better than this
# baseline does, which separates a broad-based gain from one carried by a
# handful of regions.
m_base <- feols(
  dlgrp_pc_usd ~ TM + TM_2 + RR + RR_2 |
    GID_1[year, I(year^2)] + year,
  cluster = ~GID_1,
  data = ref_dat
)
base_region_ssr <- drop(rowsum(resid(m_base)^2, ref_dat$GID_1[fixest::obs(m_base)]))

# Indicators recorded for every cell, all free (no refit):
#   SSR, logLik, nobs -- as before
#   K                 -- non-NA coefficients, catches a transform gone collinear
#   RMSE              -- readable scale only; monotone in SSR, not a criterion
#   MAE, SSR_trim     -- outlier-robust losses: does the SSR win survive?
#   RMSE_tail/_core   -- fit on the high-anomaly decile vs the rest, i.e. where
#                        the deviation term is supposed to earn its keep
#   win_rate          -- share of regions fitted better than the baseline
fit_indicators <- c("SSR", "logLik", "nobs", "K", "RMSE", "MAE",
                    "SSR_trim", "RMSE_tail", "RMSE_core", "win_rate")

score_fit <- function(m) {

  i <- fixest::obs(m)
  r <- resid(m)

  reg_ssr <- drop(rowsum(r^2, ref_dat$GID_1[i]))
  reg_base <- base_region_ssr[names(reg_ssr)]

  c(
    SSR       = sum(r^2),
    logLik    = as.numeric(logLik(m)),
    nobs      = length(r),
    K         = sum(!is.na(coef(m))),
    RMSE      = sqrt(mean(r^2)),
    MAE       = mean(abs(r)),
    SSR_trim  = sum(r[is_trim[i]]^2),
    RMSE_tail = sqrt(mean(r[is_tail[i]]^2)),
    RMSE_core = sqrt(mean(r[!is_tail[i]]^2)),
    win_rate  = mean(reg_ssr < reg_base, na.rm = TRUE)
  )
}

cat("Estimable rows:", sum(ref_est),
    " tail rows (top anomaly decile):", sum(is_tail),
    " trimmed out by outcome:", sum(ref_est) - sum(is_trim), "\n")

for (cl in fit_indicators) grid[[cl]] <- NA_real_

for (j in seq_len(nrow(grid))) {

  d <- dat_cache[[cell_key(grid$window_rr[j], grid$window_tm[j])]]

  m <- fit_lambdas(
    grid$p_T[j],
    grid$p_R[j],
    d,
    grid$moments[j],
    grid$p_T[j],
    grid$p_R[j]
  )

  grid[j, fit_indicators] <- score_fit(m)
}

grid$spec <- ifelse(
  grid$moments == "window",
  paste0("window (rr=", grid$window_rr, ", tm=", grid$window_tm, ")"),
  paste0("fixed: ", grid$moments)
)

# All specs share the same sample now, so this should be a single value
tapply(grid$nobs, grid$spec, range)

best <- grid[which.min(grid$SSR), ]

best


## ---- out-of-sample check on the leading cells ------------------------------
## In-sample SSR always improves with a more flexible deviation term, so the
## leaders are usually near-tied and the question is which of them generalises.
## 5-fold CV, folds drawn within region, answers that. It costs cv_k refits per
## cell, so it runs on the top cv_top cells by SSR rather than on the full grid.

cv_top <- 20
cv_indicators <- c("cv_rmse", "cv_mae", "cv_na")

cv_score <- function(g) {

  d <- dev_lambda_cols(
    dat_cache[[cell_key(g$window_rr, g$window_tm)]],
    g$p_T, g$p_R, g$moments, g$p_T, g$p_R
  )

  err <- rep(NA_real_, nrow(d))

  for (f in seq_len(cv_k)) {
    te <- cv_fold == f
    # The training data keeps the fold-0 carriers, so lags resolve there as well
    m_f <- fit_on(d[!te, , drop = FALSE])
    err[te] <- d$dlgrp_pc_usd[te] - predict(m_f, newdata = d[te, , drop = FALSE])
  }

  held_out <- cv_fold > 0

  c(
    cv_rmse = sqrt(mean(err[held_out]^2, na.rm = TRUE)),
    cv_mae  = mean(abs(err[held_out]), na.rm = TRUE),
    # held-out rows no training fold could predict; carriers are not counted
    cv_na   = sum(is.na(err[held_out]))
  )
}

for (cl in cv_indicators) grid[[cl]] <- NA_real_

cv_cells <- order(grid$SSR)[seq_len(min(cv_top, nrow(grid)))]

cat("Cross-validating the", length(cv_cells), "best cells by SSR,",
    cv_k, "folds each (", length(cv_cells) * cv_k, "fits )\n")

for (j in cv_cells) grid[j, cv_indicators] <- cv_score(grid[j, ])

cv_tab <- grid[cv_cells, c("p_T", "p_R", "moments",
                           "window_rr", "window_tm", "SSR", "SSR_trim", "MAE",
                           "RMSE_tail", "RMSE_core", "win_rate", "cv_rmse")]

cat("\n--- leading cells, ranked by out-of-sample RMSE ---\n")
print(cv_tab[order(cv_tab$cv_rmse), ], digits = 5, row.names = FALSE)

best_row    <- which.min(grid$SSR)
best_cv_row <- cv_cells[which.min(grid$cv_rmse[cv_cells])]

# `best` stays the in-sample (SSR) winner, so everything downstream is unchanged
# -- the CV column is a diagnostic, not a new selection rule. Disagreement
# between the two is the signal worth reading. Re-taken here only so that `best`
# carries the CV columns too.
best <- grid[best_row, ]
best_cv <- grid[best_cv_row, ]

if (best_row != best_cv_row) {
  cat("\nNOTE: the SSR winner is NOT the CV winner.\n")
  cmp <- grid[c(best_row, best_cv_row), names(cv_tab)]
  rownames(cmp) <- c("SSR winner", "CV winner")
  print(cmp, digits = 5)
} else {
  cat("\nSSR winner and CV winner agree.\n")
}

dat <- dat_cache[[
  if (is.na(best$window_rr)) ref_key else paste(best$window_rr, best$window_tm, sep = "_")
]]
best_moments <- best$moments
sT_ref <- unname(attr(dat, "refs")["sT_ref"])
sR_ref <- unname(attr(dat, "refs")["sR_ref"])


grid$delta_SSR <- grid$SSR - min(grid$SSR)

# lambda slice, held at the best exponents
ggplot(subset(grid, moments == best$moments),
       aes(x = p_T, y = p_R, fill = delta_SSR)) +
  geom_tile() +
  geom_point(data=best,
    aes(x = p_T, y = p_R),
    inherit.aes = FALSE,
    shape = 4,
    size = 4,
    stroke = 1.2
  ) +
  geom_point(data=best_cv,
             aes(x = p_T, y = p_R),
             inherit.aes = FALSE,
             shape = 4,
             size = 4,
             stroke = 1.2
  ) +
  labs(
    x = expression(p[T]),
    y = expression(p[R]),
    fill = expression(Delta*" RMSE")
  ) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  coord_equal() +
  facet_grid(window_rr ~ window_tm) +
  theme_classic()

## plots
library(ggplot2)
m <- fit_lambdas(
  best$p_T,
  best$p_R,
  dat,
  best_moments,
  best$p_T,
  best$p_R)

lambda_R <- best$p_R
p_R <- best$p_R
lambda_T <- best$p_T
p_T <- best$p_T

# Columns backing the deviation term of the selected specification
best_cols <- moment_cols[[best_moments]]

# Reference SDs used when constructing the regression variables. Must match
# add_moment_terms exactly -- estimable rows only, carriers excluded -- or the
# curves below are traced with a different scaling than the one estimated.
dat_est <- est_row(dat)
sR_ref <- median(dat[[sd_cols["sd_RR"]]][dat_est], na.rm = TRUE)
sT_ref <- median(dat[[sd_cols["sd_TM"]]][dat_est], na.rm = TRUE)

# Coefficients
b <- coef(m)

# Cluster-robust variance-covariance matrix
V <- vcov(m, cluster = ~GID_1)

# One entry per climate variable: which columns describe it, the exponent and
# scaling actually estimated, and the range over which to trace the response
curve_specs <- list(
  RR = list(
    var      = "RR",
    panel    = "Precipitation",
    mean_col = unname(best_cols["mean_RR"]),
    sd_col   = unname(sd_cols["sd_RR"]),
    s_ref    = sR_ref,
    lambda   = lambda_R,
    p        = p_R
  ),
  TM = list(
    var      = "TM",
    panel    = "Temperature",
    mean_col = unname(best_cols["mean_TM"]),
    sd_col   = unname(sd_cols["sd_TM"]),
    s_ref    = sT_ref,
    lambda   = lambda_T,
    p        = p_T
  )
)

# Representative regional means (colours) and variability (line types)
for (v in names(curve_specs)) {
  cfg <- curve_specs[[v]]

  curve_specs[[v]]$mu_values <- quantile(
    dat[[cfg$mean_col]][dat_est], c(.25, .50, .75), na.rm = TRUE
  )

  curve_specs[[v]]$sd_values <- quantile(
    dat[[cfg$sd_col]][dat_est], c(.25, .50, .75), na.rm = TRUE
  )

}

# Curves are traced in standard deviations away from the regional mean, so both
# panels share one x axis and the mean/variability levels stay comparable
z_seq <- seq(-4, 4, length.out = 300)


make_curve <- function(cfg, mu, sd_X, mean_label, var_label) {

  # This assumes your estimated regressor was:
  #
  # abs(X - mean_X)^p_X /
  # (sd_X_all / sX_ref)^lambda_X

  nm <- c(cfg$var, paste0(cfg$var, "_2"), paste0(cfg$var, "_dev_lambda"))

  sd_rel <- sd_X / cfg$s_ref
  scale_X <- sd_rel^cfg$lambda

  out <- lapply(z_seq, function(z) {

    # z standard deviations away from this curve's regional mean
    X <- mu + z * sd_X

    # Regressor values relative to X = mu
    grad <- c(
      X - mu,
      X^2 - mu^2,
      abs(X - mu)^cfg$p / scale_X
    )
    names(grad) <- nm

    # Predicted change in log GDP relative to X = mu
    effect <- unname(
      sum(b[nm] * grad)
    )

    # Delta-method variance
    Vsub <- V[nm, nm, drop = FALSE]

    se <- sqrt(
      as.numeric(
        t(grad) %*% Vsub %*% grad
      )
    )

    data.frame(
      variable = cfg$panel,
      z = z,
      x = X,
      effect = effect,
      se = se,
      lower = effect - 1.96 * se,
      upper = effect + 1.96 * se,
      mean_level = mean_label,
      variability = var_label,
      mu = unname(mu)
    )
  })

  do.call(rbind, out)
}


mean_labels <- c("Cool / dry (25th pct.)", "Median", "Warm / wet (75th pct.)")
var_labels  <- c("Low variability (25th pct.)", "Median variability",
                 "High variability (75th pct.)")

curve_dat <- do.call(rbind, lapply(curve_specs, function(cfg) {
  do.call(rbind, lapply(1:3, function(i) {
    do.call(rbind, lapply(1:3, function(k) {
      make_curve(
        cfg,
        unname(cfg$mu_values[i]),
        unname(cfg$sd_values[k]),
        mean_labels[i],
        var_labels[k]
      )
    }))
  }))
}))

# Precipitation on top, temperature below
curve_dat$variable <- factor(
  curve_dat$variable,
  levels = c("Precipitation", "Temperature")
)
curve_dat$mean_level <- factor(curve_dat$mean_level, levels = mean_labels)
curve_dat$variability <- factor(curve_dat$variability, levels = var_labels)


ggplot(
  curve_dat,
  aes(
    x = z,
    y = effect,
    colour = mean_level,
    linetype = variability
  )
) +
  # ribbon only for the median-variability curve, to keep 9 lines readable
  geom_ribbon(
    data = subset(curve_dat, variability == var_labels[2]),
    aes(
      ymin = lower,
      ymax = upper,
      fill = mean_level,
      group = interaction(mean_level, variability)
    ),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  # each curve is measured against its own regional mean
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(variable = c(
      Precipitation = "Precipitation",
      Temperature = "Temperature"
    ))
  ) +
  labs(
    x = "Deviation from regional mean (standard deviations)",
    y = expression(Delta * " log GDP per capita"),
    colour = "Regional mean climate",
    fill = "Regional mean climate",
    linetype = "Historical variability"
  ) +
  theme_classic()


## best-performing specification, written in the format of the model list in

dat$RR_dev_lambda <-
  (dat$RR - dat[[best_cols["mean_RR"]]])^best$p_R / dat$sd_RR_pre^best$p_R

dat$TM_dev_lambda <-
  (dat$TM - dat[[best_cols["mean_TM"]]])^best$p_T / dat$sd_TM_pre^best$p_T

m_dyn <- feols(
  dlgrp_pc_usd ~
    l(TM, 0:5) +
    l(TM_2, 0:5) +
    l(TM_dev_lambda, 0:5) +
    l(RR, 0:5) +
    l(RR_2, 0:5) +
    l(RR_dev_lambda, 0:5) |
    GID_1[year, I(year^2)] + year,
  panel.id = ~GID_1 + year,
  cluster = ~GID_1,
  data = dat
)

m_dyn <- feols(
  dlgrp_pc_usd ~
    l(TM, 0:5) +
    l(TM_2, 0:5) +
    l(RR, 0:5) +
    l(RR_2, 0:5)|
    GID_1[year, I(year^2)] + year,
  panel.id = ~GID_1 + year,
  cluster = ~GID_1,
  data = dat
)


m_dyn <- feols(
  dlgrp_pc_usd ~
    l(TM_dev_lambda, 0:5) +
    l(RR_dev_lambda, 0:5) |
    GID_1[year, I(year^2)] + year,
  panel.id = ~GID_1 + year,
  cluster = ~GID_1,
  data = dat
)

## ---- lag response of m_dyn -------------------------------------------------
## At each lag the response to a climate anomaly is spread over three terms
## (level, square, deviation), so the effect is a linear combination of their
## coefficients and its SE needs the full covariance block. The cumulative
## response additionally needs the cross-lag covariances.

# These must mirror how TM_dev_lambda / RR_dev_lambda were built for m_dyn
dyn_specs <- list(
  RR = list(var = "RR", panel = "Precipitation",
            mean_col = unname(best_cols["mean_RR"]), sd_col = "sd_RR_pre",
            p = best$p_R, lambda = best$lambda_R),
  TM = list(var = "TM", panel = "Temperature",
            mean_col = unname(best_cols["mean_TM"]), sd_col = "sd_TM_pre",
            p = best$p_T, lambda = best$lambda_T)
)

dyn_z <- 1                          # anomaly size, in regional SDs
dyn_percentiles <- c(0.05, 0.50, 0.95)
conf_level <- 0.95
z_crit <- qnorm(1 - (1 - conf_level) / 2)

b_dyn <- coef(m_dyn)
V_dyn <- vcov(m_dyn)                # already clustered on GID_1

# fixest names lag terms "l(var, k)" -- with a space, and lag 0 is NOT named
# plainly "var". Match on whitespace-stripped names so either form works.
b_dyn_norm <- gsub("[[:space:]]", "", names(b_dyn))

lag_term <- function(var, k) {
  hit <- names(b_dyn)[b_dyn_norm == paste0("l(", var, ",", k, ")")]
  if (length(hit) == 0 && k == 0) hit <- names(b_dyn)[b_dyn_norm == var]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Read the lag depth off the fitted model instead of hard-coding it
dyn_term_vars <- unlist(lapply(dyn_specs, function(cfg) {
  c(cfg$var, paste0(cfg$var, "_2"), paste0(cfg$var, "_dev_lambda"))
}), use.names = FALSE)

dyn_lags <- {
  pat <- paste0("^l\\((", paste(dyn_term_vars, collapse = "|"), "),([0-9]+)\\)$")
  ks <- as.integer(sub(pat, "\\2", grep(pat, b_dyn_norm, value = TRUE)))
  if (!length(ks)) stop("no lag terms of the expected form found in m_dyn")
  0:max(ks)
}

# Fail loudly rather than silently returning zero effects
dyn_missing <- unlist(lapply(dyn_term_vars, function(v) {
  k <- dyn_lags[is.na(vapply(dyn_lags, function(k) lag_term(v, k), character(1)))]
  if (length(k)) paste0(v, " @ lag ", paste(k, collapse = ",")) else NULL
}))

if (length(dyn_missing)) {
  stop("m_dyn has no coefficient for: ", paste(dyn_missing, collapse = "; "))
}

# Regional mean climate is a covariate of the response (through the square and
# the centring), so curves are traced at percentiles of it across regions. Each
# percentile keeps that region's own SD, so the pair stays internally consistent.
region_levels <- function(cfg) {
  rc <- dat %>%
    ungroup() %>%
    group_by(GID_1) %>%
    summarise(
      mu = median(.data[[cfg$mean_col]], na.rm = TRUE),
      sd = median(.data[[cfg$sd_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(mu) & !is.na(sd)) %>%
    arrange(mu)

  rc[pmax(1, ceiling(dyn_percentiles * nrow(rc))), ]
}

lag_response <- function(cfg, mu, sd_X, level_label) {

  # a +dyn_z SD anomaly, expressed in each term's own units
  x <- mu + dyn_z * sd_X
  w_base <- c(
    x - mu,
    x^2 - mu^2,
    (x - mu)^cfg$p / sd_X^cfg$lambda
  )
  vars <- c(cfg$var, paste0(cfg$var, "_2"), paste0(cfg$var, "_dev_lambda"))

  estimate <- function(nm, w) {
    if (length(nm) == 0) stop("no coefficients matched for ", cfg$var)
    eff <- sum(b_dyn[nm] * w)
    se <- sqrt(max(as.numeric(t(w) %*% V_dyn[nm, nm, drop = FALSE] %*% w), 0))
    c(effect = eff, se = se)
  }

  nm_cum <- character(0)
  w_cum <- numeric(0)

  out <- lapply(dyn_lags, function(k) {

    nm_k <- vapply(vars, lag_term, character(1), k = k)
    keep <- !is.na(nm_k)
    nm_k <- unname(nm_k[keep])
    w_k <- w_base[keep]

    # accumulate for the cumulative response through lag k
    nm_cum <<- c(nm_cum, nm_k)
    w_cum <<- c(w_cum, w_k)

    a <- estimate(nm_k, w_k)
    cu <- estimate(nm_cum, w_cum)

    data.frame(
      variable = cfg$panel,
      lag = k,
      horizon = c("Effect at lag", "Cumulative through lag"),
      level = level_label,
      effect = c(a["effect"], cu["effect"]),
      se = c(a["se"], cu["se"]),
      row.names = NULL
    )
  })

  do.call(rbind, out)
}

level_labels <- paste0(round(100 * dyn_percentiles), "th pct.")

lag_dat <- do.call(rbind, lapply(dyn_specs, function(cfg) {
  rl <- region_levels(cfg)
  do.call(rbind, lapply(seq_len(nrow(rl)), function(i) {
    lag_response(cfg, rl$mu[i], rl$sd[i], level_labels[i])
  }))
}))

lag_dat$lower <- lag_dat$effect - z_crit * lag_dat$se
lag_dat$upper <- lag_dat$effect + z_crit * lag_dat$se

lag_dat$variable <- factor(lag_dat$variable, levels = c("Precipitation", "Temperature"))
lag_dat$horizon <- factor(lag_dat$horizon, levels = c("Effect at lag", "Cumulative through lag"))
lag_dat$level <- factor(lag_dat$level, levels = level_labels)


ggplot(lag_dat, aes(x = lag, y = effect, colour = level, fill = level, group = level)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_grid(variable ~ horizon, scales = "free_y") +
  scale_x_continuous(breaks = dyn_lags) +
  labs(
    x = "Lag (years)",
    y = bquote(Delta ~ "log GDP per capita, per +" * .(dyn_z) * " SD anomaly"),
    colour = "Regional mean climate",
    fill = "Regional mean climate",
    caption = paste0("Ribbons = ", round(100 * conf_level),
                     "% CI, delta method over the level, square and deviation terms")
  ) +
  theme_classic()


## ---- asymmetry: positive vs negative deviations ----------------------------
## The deviation term is symmetric by construction -- abs(X - mean)^p carries a
## single coefficient, so a hot year and a cold year of equal size are forced to
## have the same effect. Split it at the regional mean and test whether the two
## halves actually share sign and magnitude. Uses the refined (round 2) optimum.

asym_dev <- function(var, mo, p, lambda) {
  dev <- dat[[paste0(var, "_absdev_", mo)]]^p /
    dat[[paste0("sd_", var, "_rel")]]^lambda
  above <- dat[[var]] >= dat[[moment_cols[[mo]][paste0("mean_", var)]]]
  list(pos = ifelse(above, dev, 0), neg = ifelse(above, 0, dev))
}

dev_T <- asym_dev("TM", best_moments, best$p_T, best$lambda_T)
dev_R <- asym_dev("RR", best_moments, best$p_R, best$lambda_R)

dat$TM_dev_pos <- dev_T$pos
dat$TM_dev_neg <- dev_T$neg
dat$RR_dev_pos <- dev_R$pos
dat$RR_dev_neg <- dev_R$neg

m_asym <- feols(
  dlgrp_pc_usd ~
    TM + TM_2 +
    RR + RR_2 +
    TM_dev_pos + TM_dev_neg +
    RR_dev_pos + RR_dev_neg |
    GID_1[year, I(year^2)] + year,
  cluster = ~GID_1,
  data = dat
)

# the symmetric model is nested in this one, so the samples must coincide
stopifnot(nobs(m_asym) == nobs(m))

asym_terms <- c("TM_dev_pos", "TM_dev_neg", "RR_dev_pos", "RR_dev_neg")

# cluster-robust inference: t with (number of clusters - 1) dof
dof_asym <- length(unique(dat$GID_1)) - 1

lincom <- function(mod, w, label) {
  nm <- names(w)
  est <- sum(coef(mod)[nm] * w)
  se <- sqrt(as.numeric(t(w) %*% vcov(mod)[nm, nm, drop = FALSE] %*% w))
  tstat <- est / se
  data.frame(
    test = label,
    estimate = est,
    se = se,
    lower = est - qt(0.975, dof_asym) * se,
    upper = est + qt(0.975, dof_asym) * se,
    t = tstat,
    p = 2 * pt(-abs(tstat), dof_asym),
    row.names = NULL
  )
}

# the two halves, each on its own
asym_coefs <- do.call(rbind, lapply(asym_terms, function(v) {
  w <- 1
  names(w) <- v
  lincom(m_asym, w, v)
}))
asym_coefs$sign <- ifelse(asym_coefs$estimate >= 0, "+", "-")
asym_coefs$signif <- ifelse(asym_coefs$p < 0.05, "*", "")

# pos - neg = 0  ->  the symmetric restriction the pooled model imposes
# pos + neg = 0  ->  equal magnitude but opposite sign
asym_tests <- rbind(
  lincom(m_asym, c(TM_dev_pos = 1, TM_dev_neg = -1), "TM: pos - neg (symmetry)"),
  lincom(m_asym, c(TM_dev_pos = 1, TM_dev_neg =  1), "TM: pos + neg (mirror)"),
  lincom(m_asym, c(RR_dev_pos = 1, RR_dev_neg = -1), "RR: pos - neg (symmetry)"),
  lincom(m_asym, c(RR_dev_pos = 1, RR_dev_neg =  1), "RR: pos + neg (mirror)")
)

# joint test that BOTH deviation terms are symmetric
R_sym <- rbind(
  c(1, -1, 0, 0),
  c(0, 0, 1, -1)
)
colnames(R_sym) <- asym_terms

Rb <- R_sym %*% coef(m_asym)[asym_terms]
W_sym <- as.numeric(
  t(Rb) %*% solve(R_sym %*% vcov(m_asym)[asym_terms, asym_terms] %*% t(R_sym)) %*% Rb
)

asym_joint <- data.frame(
  test = "joint: both deviation terms symmetric",
  chisq = W_sym,
  df = nrow(R_sym),
  p = pchisq(W_sym, nrow(R_sym), lower.tail = FALSE)
)

cat("\n--- deviation coefficients, split at the regional mean ---\n")
print(asym_coefs, digits = 4)
cat("\n--- linear hypotheses ---\n")
print(asym_tests, digits = 4)
cat("\n--- joint symmetry test ---\n")
print(asym_joint, digits = 4)
cat("\nSSR symmetric:", sum(resid(m)^2),
    " asymmetric:", sum(resid(m_asym)^2), "\n")


## outliers 
m <- fit_lambdas(
  2,
  2,
  # na.rm: the lag carriers have no outcome, and keeping them (outcome NA, so
  # never estimated on) is what lets the lags resolve in this refit too
  dat %>% filter(is.na(dlgrp_pc_usd) |
                   (dlgrp_pc_usd < quantile(dat$dlgrp_pc_usd, 0.95, na.rm = TRUE) &
                      dlgrp_pc_usd > quantile(dat$dlgrp_pc_usd, 0.05, na.rm = TRUE))),
  best_moments,
  2,
  2)



## ---- projection along a climate trajectory ---------------------------------
## Uses the best static model `m` (no lags, symmetric deviation term). Growth is
## measured against the counterfactual where climate stays frozen at (T0, P0):
## there the rolling mean equals T0/P0 and the deviation term is exactly zero,
## so that baseline contributes nothing and drops out of every gradient.

# Rolling mean over the `window` years ENDING AT t-1, with the pre-trajectory
# years held at x0. Lagged, matching the lagged centring used in estimation, so
# year t is excluded from its own reference mean here too.
roll_mean_pad <- function(x, x0, window) {
  w <- if (is.na(window) || window < 1) 1 else window
  padded <- c(rep(x0, w), head(x, -1))
  vapply(seq_along(x), function(t) mean(padded[t:(t + w - 1)]), numeric(1))
}

# convenience: constant trend starting from x0, so year 1 still equals x0
linear_traj <- function(x0, trend, n) x0 + trend * (seq_len(n) - 1)

project_climate <- function(T_traj, P_traj,
                            T0, P0, sd_T, sd_P,
                            years = seq_along(T_traj),
                            model = m,
                            moments = best_moments,
                            p_T = best$p_T, lambda_T = best$lambda_T,
                            p_R = best$p_R, lambda_R = best$lambda_R,
                            window_T = best$window_tm, window_R = best$window_rr,
                            s_ref_T = sT_ref, s_ref_R = sR_ref,
                            conf_level = 0.95) {

  stopifnot(
    length(T_traj) == length(P_traj),
    length(years) == length(T_traj),
    sd_T > 0, sd_P > 0
  )

  # centring of the deviation term, mirroring how the model was estimated
  if (identical(moments, "window")) {
    stopifnot(!is.na(window_T), !is.na(window_R))
    mean_T <- roll_mean_pad(T_traj, T0, window_T)
    mean_P <- roll_mean_pad(P_traj, P0, window_R)
  } else {
    # fixed historical centring: by assumption the pre-trajectory climate was
    # constant at T0 / P0, so that IS the historical mean
    mean_T <- rep(T0, length(T_traj))
    mean_P <- rep(P0, length(P_traj))
  }

  # SD is constant, so the scaling factor is the same in every year
  dev_T <- abs(T_traj - mean_T)^p_T / (sd_T / s_ref_T)^lambda_T
  dev_P <- abs(P_traj - mean_P)^p_R / (sd_P / s_ref_R)^lambda_R

  nm <- c("TM", "TM_2", "TM_dev_lambda", "RR", "RR_2", "RR_dev_lambda")
  stopifnot(all(nm %in% names(coef(model))))

  # one row per year: the regressor change relative to the frozen-climate baseline
  G <- cbind(
    TM            = T_traj - T0,
    TM_2          = T_traj^2 - T0^2,
    TM_dev_lambda = dev_T,
    RR            = P_traj - P0,
    RR_2          = P_traj^2 - P0^2,
    RR_dev_lambda = dev_P
  )[, nm, drop = FALSE]

  b <- coef(model)[nm]
  V <- vcov(model)[nm, nm, drop = FALSE]
  z_c <- qnorm(1 - (1 - conf_level) / 2)

  # annual growth effect, and its cumulative sum = total log GDP deviation.
  # Delta method: diag(G V G') for each, so the cumulative CI carries the
  # covariance between years rather than treating them as independent.
  G_cum <- apply(G, 2, cumsum)
  if (is.null(dim(G_cum))) G_cum <- matrix(G_cum, nrow = 1, dimnames = list(NULL, nm))

  growth <- as.numeric(G %*% b)
  cumul <- as.numeric(G_cum %*% b)

  out <- rbind(
    data.frame(year = years, quantity = "Annual growth effect",
               value = growth,
               se = sqrt(pmax(rowSums((G %*% V) * G), 0))),
    data.frame(year = years, quantity = "Cumulative effect",
               value = cumul,
               se = sqrt(pmax(rowSums((G_cum %*% V) * G_cum), 0)))
  )

  out$lower <- out$value - z_c * out$se
  out$upper <- out$value + z_c * out$se
  out$quantity <- factor(out$quantity,
                         levels = c("Annual growth effect", "Cumulative effect"))

  gg <- ggplot(out, aes(x = year, y = value)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "steelblue") +
    geom_line(linewidth = 0.9, colour = "steelblue4") +
    facet_wrap(~ quantity, ncol = 1, scales = "free_y") +
    labs(
      x = "Year",
      y = expression(Delta * " log GDP per capita"),
      caption = paste0("Relative to climate frozen at (T0, P0); ribbons = ",
                       round(100 * conf_level), "% CI")
    ) +
    theme_classic()

  list(
    plot = gg,
    data = out,
    climate = data.frame(year = years, TM = T_traj, RR = P_traj,
                         mean_TM = mean_T, mean_RR = mean_P,
                         dev_TM = dev_T, dev_RR = dev_P)
  )
}


## example: a median region warming 0.04 C/yr with slowly rising precipitation
proj_n <- 40
T0_ex <- 20
P0_ex <- median(dat[[best_cols["mean_RR"]]][dat_est], na.rm = TRUE)
sd_T_ex <- median(dat[[sd_cols["sd_TM"]]][dat_est], na.rm = TRUE)
sd_P_ex <- median(dat[[sd_cols["sd_RR"]]][dat_est], na.rm = TRUE)

pulse_T <- rep(T0_ex,proj_n)
pulse_T[2] <- pulse_T[2] + 3*sd_T_ex

proj <- project_climate(
  T_traj = linear_traj(T0_ex, 0.04, proj_n),
  P_traj = linear_traj(P0_ex, 1, proj_n),
  T0 = T0_ex, P0 = P0_ex,
  sd_T = sd_T_ex, sd_P = sd_P_ex
)


proj <- project_climate(
  T_traj = pulse_T,
  P_traj = rep(P0_ex,proj_n),
  T0 = T0_ex, P0 = P0_ex,
  sd_T = sd_T_ex, sd_P = sd_P_ex
)

proj$plot


## ---- several trajectories on one plot --------------------------------------
## Input is a long data frame, one row per trajectory-year:
##   name             label shown in the legend
##   T_traj, P_traj   the yearly temperature / precipitation values
##   T0, P0           initial conditions   (constant within name)
##   SDt, SDp         constant SDs         (constant within name)
## an optional `year` column sets the x axis, otherwise years run 1..n.
## Extra arguments are passed straight through to project_climate().

project_climate_multi <- function(traj_df, conf_level = 0.95, ...) {

  req <- c("name", "T_traj", "P_traj", "T0", "P0")
  miss <- setdiff(req, names(traj_df))
  if (length(miss)) stop("traj_df is missing: ", paste(miss, collapse = ", "))

  # accept either SDt/SDp or sd_T/sd_P
  if (!"SDt" %in% names(traj_df) && "sd_T" %in% names(traj_df)) traj_df$SDt <- traj_df$sd_T
  if (!"SDp" %in% names(traj_df) && "sd_P" %in% names(traj_df)) traj_df$SDp <- traj_df$sd_P
  if (!all(c("SDt", "SDp") %in% names(traj_df))) stop("traj_df needs SDt and SDp")

  traj_df$name <- as.character(traj_df$name)
  nms <- unique(traj_df$name)

  res <- lapply(nms, function(nm_i) {

    g <- traj_df[traj_df$name == nm_i, , drop = FALSE]

    # the four scalars describe the trajectory, so they must not vary within it
    for (v in c("T0", "P0", "SDt", "SDp")) {
      if (length(unique(g[[v]])) != 1L) {
        stop(v, " is not constant within trajectory '", nm_i, "'")
      }
    }

    yrs <- if ("year" %in% names(g)) g$year else seq_len(nrow(g))
    o <- order(yrs)

    r <- project_climate(
      T_traj = g$T_traj[o], P_traj = g$P_traj[o],
      T0 = g$T0[1], P0 = g$P0[1],
      sd_T = g$SDt[1], sd_P = g$SDp[1],
      years = yrs[o],
      conf_level = conf_level,
      ...
    )

    r$data$name <- nm_i
    r$climate$name <- nm_i
    r
  })

  dat_out <- do.call(rbind, lapply(res, `[[`, "data"))
  clim_out <- do.call(rbind, lapply(res, `[[`, "climate"))
  dat_out$name <- factor(dat_out$name, levels = nms)
  clim_out$name <- factor(clim_out$name, levels = nms)

  gg <- ggplot(dat_out, aes(x = year, y = value, colour = name, fill = name)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ quantity, ncol = 1, scales = "free_y") +
    labs(
      x = "Year",
      y = expression(Delta * " log GDP per capita"),
      colour = "Trajectory",
      fill = "Trajectory",
      caption = paste0("Relative to each trajectory's own frozen climate (T0, P0); ribbons = ",
                       round(100 * conf_level), "% CI")
    ) +
    theme_classic()

  list(plot = gg, data = dat_out, climate = clim_out)
}


## example: three warming rates for a median region
make_traj <- function(name, trend_T, trend_P, n = 100) {
  data.frame(
    name = name,
    year = seq_len(n),
    T_traj = linear_traj(T0_ex, trend_T, n),
    P_traj = linear_traj(P0_ex, trend_P, n),
    T0 = T0_ex, P0 = P0_ex,
    SDt = sd_T_ex, SDp = sd_P_ex,
    stringsAsFactors = FALSE
  )
}

traj_df <- rbind(
  make_traj("Low (0.02 C/yr)",  0.02, 0.000),
  make_traj("Mid (0.04 C/yr)",  0.04, 0.002),
  make_traj("High (0.06 C/yr)", 0.06, 0.004)
)

proj_multi <- project_climate_multi(traj_df)

proj_multi$plot


## ---- step vs ramp: the same final +3 SD shock, reached at different speeds --
## Both trajectories end at exactly the same climate, so any difference is pure
## timing: the deviation term depends on the gap between the trajectory and its
## own trailing mean, which a sudden step opens wide and a slow ramp does not.

shock_n <- 100        # years simulated
ramp_years <- 20     # years the gradual trajectory takes to get there
shock_sd <- 3        # size of the final shock, in SDs

step_delta <- function(n, total) rep(total, n)
ramp_delta <- function(n, total, ramp) total * pmin(seq_len(n) / ramp, 1)

make_shock <- function(name, dT, dP) {
  data.frame(
    name = name,
    year = seq_len(shock_n),
    T_traj = T0_ex + dT,
    P_traj = P0_ex + dP,
    T0 = T0_ex, P0 = P0_ex,
    SDt = sd_T_ex, SDp = sd_P_ex,
    stringsAsFactors = FALSE
  )
}

no_change <- rep(0, shock_n)

shock_df <- rbind(
  make_shock("T: step at t1",
             step_delta(shock_n, shock_sd * sd_T_ex), no_change),
  make_shock(paste0("T: ramp over ", ramp_years, " yr"),
             ramp_delta(shock_n, shock_sd * sd_T_ex, ramp_years), no_change),
  make_shock("P: step at t1",
             no_change, step_delta(shock_n, -shock_sd * sd_P_ex)),
  make_shock(paste0("P: ramp over ", ramp_years, " yr"),
             no_change, ramp_delta(shock_n,-shock_sd * sd_P_ex, ramp_years))
)

proj_shock <- project_climate_multi(shock_df)

proj_shock$plot

# what actually drives the difference: the anomaly against the trailing mean
shock_anom <- rbind(
  data.frame(proj_shock$climate[, c("name", "year")],
             variable = "Temperature",
             anomaly = with(proj_shock$climate, (TM - mean_TM) / sd_T_ex)),
  data.frame(proj_shock$climate[, c("name", "year")],
             variable = "Precipitation",
             anomaly = with(proj_shock$climate, (RR - mean_RR) / sd_P_ex))
)

# keep each variable's own shocks in its own panel
shock_anom <- subset(
  shock_anom,
  (variable == "Temperature"   & grepl("^T:", name)) |
  (variable == "Precipitation" & grepl("^P:", name))
)

ggplot(shock_anom, aes(x = year, y = anomaly, colour = name)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ variable, ncol = 1) +
  labs(
    x = "Year",
    y = "Anomaly vs trailing mean (SD)",
    colour = "Trajectory",
    caption = "The deviation term is driven by this gap, not by the level of the shock"
  ) +
  theme_classic()


## ---- residualised TM_2 vs dev_TM_all_nosd_2 --------------------------------
## Both are squared temperature terms: TM_2 = TM^2 and dev_TM_all_nosd_2 =
## (TM - mean_TM_all)^2. They differ only by a region-specific linear term and a
## constant, so the question is how much INDEPENDENT variation survives once the
## fixed effects are absorbed. Partialling out with feols(x ~ 1 | FE) and
## plotting one residual against the other answers that directly.

resid_vars <- c("TM_2", "dev_TM_all_nosd_2")

# same rows for both, so the residuals are comparable point by point
dat_res <- dat[complete.cases(as.data.frame(ungroup(dat))[, resid_vars, drop = FALSE]), ]

# feols with an empty RHS just absorbs the fixed effects; resid() is the
# within-transformed variable (fixest::demean() is the lower-level equivalent)
res_list <- lapply(resid_vars, function(v) {
  f <- as.formula(paste0(v, " ~ 1 | GID_1[year, I(year^2)] + year"))
  as.numeric(resid(feols(f, data = dat_res)))
})
names(res_list) <- resid_vars

res_dat <- data.frame(
  GID_1 = dat_res$GID_1[-1],
  year = dat_res$year[-1],
  TM_2_res = res_list[["TM_2"]],
  dev_res = res_list[["dev_TM_all_nosd_2"]]
)

res_cor <- cor(res_dat$TM_2_res, res_dat$dev_res)
res_fit <- lm(dev_res ~ TM_2_res, data = res_dat)

cat("\n--- residualised on GID_1[year, year^2] + year ---\n")
cat("n =", nrow(res_dat), "\n")
cat("correlation =", round(res_cor, 4), "  R2 =", round(summary(res_fit)$r.squared, 4), "\n")
cat("sd(TM_2 resid) =", round(sd(res_dat$TM_2_res), 4),
    " sd(dev resid) =", round(sd(res_dat$dev_res), 4), "\n")

ggplot(res_dat, aes(x = TM_2_res, y = dev_res)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(alpha = 0.06, size = 0.5, colour = "steelblue4") +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick", linewidth = 0.8) +
  labs(
    x = expression("residualised  TM"^2),
    y = expression("residualised  (TM - " * bar(TM)[all] * ")"^2),
    subtitle = paste0("r = ", round(res_cor, 3), ",  n = ", nrow(res_dat),
                      "  (fixed effects: GID_1[year, year^2] + year)")
  ) +
  theme_classic()

