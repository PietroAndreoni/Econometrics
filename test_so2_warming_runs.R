# Air pollution and the warming-run response of national growth.
#
# The warming/cooling-run model of main_analysis.Rmd (Test 10, m_bool), with
# changes in SO2 emissions as a proxy for changes in air pollution:
#
#   M0  TM + TM:mean_TM_all + RR + RR:mean_RR_all + zRR^2
#         + zTMp_f^2 + zTMp_f^2:i(warming_run_length_capped, ref = 1)
#         + i(warming_run_length_capped, ref = 1)
#         + zTMm_f^2 + zTMm_f^2:i(cooling_run_length_capped, ref = 1)
#         + i(cooling_run_length_capped, ref = 1)
#   M1  M0 + dlSO2_pos + dlSO2_neg
#   M2  M0 + lags 0-3 of dlSO2_pos and of dlSO2_neg
#
# dlSO2 is the log change in national SO2 emissions (CEDS via OWID; see
# load_so2_emissions() in functions/air_pollution.R), split into increases
# dlSO2_pos = max(dlSO2, 0), mostly tied to output growth, and decreases
# dlSO2_neg = min(dlSO2, 0), mostly air-pollution control. Runs and
# fixed-baseline anomalies are built by add_climate_runs() exactly as in the
# notebook (lagged 30-year baseline, run position capped at 6). Fixed effects
# year + GID_1 + GID_1[year], errors clustered by country (GID_1 = GID_0 at
# GADM0).
#
# All three models are fitted on one sample: the country-years M2 can use. The
# SO2 coefficients are not causal (output drives contemporaneous emissions);
# the question is whether controlling for pollution changes moves the
# warming-run response of the squared warm anomaly.
#
# Climate side: on the same sample, each climate outcome (annual temperature
# change, the anomaly zTM, and the warming- and cooling-run positions) is
# regressed on lags 0-3 of dlSO2_pos and dlSO2_neg with country and year fixed
# effects. Year effects absorb global temperature, including the effect of
# world aerosol emissions, so the coefficients relate a country's own SO2
# changes to its own temperature. Contemporaneous terms can also run the other
# way (cold years raise heating emissions).
#
# Cross-model tests: M0 and M1 (or M2) are stacked with model-specific fixed
# effects and coefficients, which reproduces each model's estimates exactly,
# and clustering by country gives the covariance between the two models'
# coefficients. Differences in the run-position effects are then tested with
# Wald tests on that stacked covariance.
#
# Run with  Rscript test_so2_warming_runs.R  (WB and PWT, ERA5 area weights).
# Output: results/so2_warming_runs/

source("load_functions.R")

OUTPUT_DIR <- file.path("results", "so2_warming_runs")
ECON_DATASETS <- c("WB", "PWT")
ECON_YEAR_MIN <- 1979L
CLIMATE_VELOCITY_YEARS <- 5L
RUN_LENGTH_CAP <- 6L
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

base_climate <- "TM + TM:mean_TM_all + RR + RR:mean_RR_all"
model_bool <- paste(
  base_climate,
  "+ zRR^2",
  "+ zTMp_f^2 + zTMp_f^2:i(warming_run_length_capped, ref = 1)",
  "+ i(warming_run_length_capped, ref = 1)",
  "+ zTMm_f^2 + zTMm_f^2:i(cooling_run_length_capped, ref = 1)",
  "+ i(cooling_run_length_capped, ref = 1)"
)
SO2_LAGS <- 0:3
so2_terms <- function(sign, lags) {
  v <- paste0("dlSO2_", sign)
  ifelse(lags == 0L, v, paste0(v, "_l", lags))
}
SO2_TERMS <- list(
  M1 = c(so2_terms("pos", 0L), so2_terms("neg", 0L)),
  M2 = c(so2_terms("pos", SO2_LAGS), so2_terms("neg", SO2_LAGS))
)
MODELS <- c(
  M0 = model_bool,
  M1 = paste(model_bool, "+", paste(SO2_TERMS$M1, collapse = " + ")),
  M2 = paste(model_bool, "+", paste(SO2_TERMS$M2, collapse = " + "))
)

CLIMATE_FIXED_EFFECTS <- "year + GID_1"
CLIMATE_OUTCOMES <- c(
  dTM_run = "Annual temperature change (deg C)",
  zTM = "Temperature anomaly zTM (SD)",
  warming_run_length_capped = "Warming-run position (years)",
  cooling_run_length_capped = "Cooling-run position (years)"
)

# Sums of the increase and decrease terms over `lags`, and the test that the
# two sums are equal.
so2_sign_sums <- function(model, lags, label) {
  w_pos <- stats::setNames(rep(1, length(lags)), so2_terms("pos", lags))
  w_neg <- stats::setNames(rep(1, length(lags)), so2_terms("neg", lags))
  bind_rows(
    lincom(model, w_pos, "increases"),
    lincom(model, w_neg, "decreases"),
    lincom(model, c(w_pos, -w_neg), "increases - decreases")
  ) %>%
    mutate(model = label, lags = paste(unique(range(lags)), collapse = "-"), .before = 1)
}

fit_model <- function(terms, model_data) {
  fit_panel_model(
    terms, model_data,
    fixed_effects = FIXED_EFFECTS,
    panel_id = c("GID_1", "year"),
    cluster = ~GID_1
  )
}

# Effect of a squared anomaly at each run position: the main coefficient at
# position 1, plus the position's interaction thereafter.
run_position_weights <- function(model, anomaly, run) {
  main <- paste0("I(", anomaly, "^2)")
  inter <- names(coef(model))
  inter <- inter[startsWith(inter, paste0(main, ":", run, "_run_length_capped::"))]
  positions <- c(1L, as.integer(sub(".*::", "", inter)))
  G <- matrix(0, length(positions), 1L + length(inter),
              dimnames = list(NULL, c(main, inter)))
  G[, main] <- 1
  for (k in seq_along(inter)) G[k + 1L, inter[[k]]] <- 1
  list(G = G, position = positions)
}

run_position_effects <- function(model, label) {
  bind_rows(lapply(
    list(c("zTMp_f", "warming"), c("zTMm_f", "cooling")),
    function(x) {
      w <- run_position_weights(model, x[[1]], x[[2]])
      linear_combination(model, w$G) %>%
        mutate(model = label, run = x[[2]], position = w$position, .before = 1)
    }
  ))
}

# M0 and an alternative stacked as one regression: every regressor and fixed
# effect is model-specific, so the coefficients equal the separate fits, and
# the country clusters span both copies of the sample.
fit_stacked <- function(fit_a, fit_b, sample) {
  blocks <- list(a = fit_a, b = fit_b)
  X <- lapply(blocks, function(f) model.matrix(f, type = "rhs"))
  stopifnot(nrow(X$a) == nrow(sample), nrow(X$b) == nrow(sample))
  rhs_names <- unlist(lapply(names(X), function(k) paste0(k, "__", colnames(X[[k]]))))
  safe_names <- paste0("x", seq_along(rhs_names))
  n <- nrow(sample)
  Z <- matrix(0, 2L * n, length(rhs_names))
  Z[seq_len(n), seq_len(ncol(X$a))] <- X$a
  Z[n + seq_len(n), ncol(X$a) + seq_len(ncol(X$b))] <- X$b
  colnames(Z) <- safe_names
  stacked <- bind_cols(
    tibble::tibble(
      y = rep(sample$dlgrp_pc_usd, 2L),
      GID_1 = rep(sample$GID_1, 2L),
      year = rep(sample$year, 2L),
      block = rep(c("a", "b"), each = n)
    ),
    tibble::as_tibble(Z)
  )
  fit <- fixest::feols(
    stats::as.formula(paste(
      "y ~", paste(safe_names, collapse = " + "),
      "| year^block + GID_1^block + GID_1^block[year]"
    )),
    data = stacked, cluster = ~GID_1,
    ssc = fixest::ssc(adj = FALSE)
  )
  b <- stats::setNames(coef(fit)[safe_names], rhs_names)
  V <- vcov(fit)[safe_names, safe_names]
  dimnames(V) <- list(rhs_names, rhs_names)
  separate <- c(stats::setNames(coef(fit_a), paste0("a__", names(coef(fit_a)))),
                stats::setNames(coef(fit_b), paste0("b__", names(coef(fit_b)))))
  max_gap <- max(abs(b[names(separate)] - separate))
  if (max_gap > 1e-6) stop("Stacked and separate coefficients differ by ", max_gap)
  list(coef = b, vcov = V, n_clusters = fixest::fitstat(fit, "g")[[1]])
}

# Differences b - a of the run-position effects, one per position, and the
# joint Wald test that all positions are unchanged.
compare_run_effects <- function(stacked, fit_a, run, anomaly, label) {
  w <- run_position_weights(fit_a, anomaly, run)
  G <- cbind(-w$G, w$G)
  colnames(G) <- c(paste0("a__", colnames(w$G)), paste0("b__", colnames(w$G)))
  b <- stacked$coef[colnames(G)]
  V <- stacked$vcov[colnames(G), colnames(G)]
  d <- as.numeric(G %*% b)
  Vd <- G %*% V %*% t(G)
  se <- sqrt(diag(Vd))
  stat <- as.numeric(t(d) %*% solve(Vd, d))
  q <- length(d)
  g <- stacked$n_clusters
  F_stat <- stat / q
  list(
    positions = tibble::tibble(
      comparison = label, run = run, position = w$position,
      difference = d, std_error = se, statistic = d / se,
      p_value = 2 * stats::pnorm(-abs(d / se))
    ),
    joint = tibble::tibble(
      comparison = label, run = run, test = "all run positions unchanged",
      F_stat = F_stat, df1 = q, df2 = g - 1L,
      p_value = stats::pf(F_stat, q, g - 1L, lower.tail = FALSE)
    )
  )
}

so2 <- load_so2_emissions(lags = 1:3)

panel_cfg <- panel_config(
  climate_source = "era5",
  climate_weight = "area",
  econ_year_min = ECON_YEAR_MIN,
  clim_history_years = max(6L, CLIMATE_VELOCITY_YEARS - 1L)
)

coefficients <- list()
run_effects <- list()
so2_tests <- list()
so2_sums <- list()
climate_coefs <- list()
climate_sums <- list()
diff_positions <- list()
diff_joint <- list()
samples <- list()
tables <- character()

for (econ in ECON_DATASETS) {
  message("Fitting ", econ)
  data <- build_dat(econ_data = econ, config = panel_cfg)
  runs <- add_climate_runs(
    data,
    velocity_years = CLIMATE_VELOCITY_YEARS,
    baseline_mean = "mean_TM_lag_30",
    baseline_sd = "sd_TM_lag_30"
  )
  check_climate_runs(runs)
  runs <- runs %>%
    mutate(
      warming_run_length_capped = pmin(RUN_LENGTH_CAP, warming_run_length),
      cooling_run_length_capped = pmin(RUN_LENGTH_CAP, cooling_run_length)
    ) %>%
    left_join(select(so2, GID_0, year, starts_with("dlSO2")), by = c("GID_0", "year"))

  # The common sample: every country-year M2 can use. The models have no
  # panel lags, so the rows outside it carry nothing the fits need.
  m2_full <- fit_model(MODELS[["M2"]], runs)
  sample <- runs[fixest::obs(m2_full), ]
  fits <- lapply(MODELS, fit_model, model_data = sample)
  n_obs <- vapply(fits, stats::nobs, integer(1))
  if (length(unique(n_obs)) != 1L) stop(econ, ": samples differ across models")

  samples[[econ]] <- tibble::tibble(
    econ = econ,
    observations = n_obs[[1]],
    countries = n_distinct(sample$GID_0),
    first_year = min(sample$year),
    last_year = max(sample$year),
    country_years_without_so2 = sum(est_row(runs) & is.na(runs$dlSO2)),
    dlSO2_sd = stats::sd(sample$dlSO2),
    share_increases = mean(sample$dlSO2 > 0),
    share_decreases = mean(sample$dlSO2 < 0),
    mean_increase = mean(sample$dlSO2[sample$dlSO2 > 0]),
    mean_decrease = mean(sample$dlSO2[sample$dlSO2 < 0])
  )

  coefficients[[econ]] <- bind_rows(lapply(names(fits), function(k) {
    tidy_fixest(fits[[k]], label = k) %>% mutate(econ = econ, .before = 1)
  }))
  run_effects[[econ]] <- bind_rows(lapply(names(fits), function(k) {
    run_position_effects(fits[[k]], k) %>% mutate(econ = econ, .before = 1)
  }))

  so2_tests[[econ]] <- bind_rows(
    wald_row(fits$M1, "^dlSO2_", label = "M1: all SO2 terms = 0"),
    wald_row(fits$M2, "^dlSO2_", label = "M2: all SO2 terms = 0"),
    wald_row(fits$M2, "^dlSO2_pos", label = "M2: increase terms = 0"),
    wald_row(fits$M2, "^dlSO2_neg", label = "M2: decrease terms = 0")
  ) %>%
    mutate(econ = econ, .before = 1)
  so2_sums[[econ]] <- bind_rows(
    so2_sign_sums(fits$M1, 0L, "M1"),
    so2_sign_sums(fits$M2, SO2_LAGS, "M2")
  ) %>%
    mutate(econ = econ, .before = 1)

  # Climate side: the same sample, country and year effects.
  for (y in names(CLIMATE_OUTCOMES)) {
    fit_c <- fixest::feols(
      stats::as.formula(paste(
        y, "~", paste(SO2_TERMS$M2, collapse = " + "), "|", CLIMATE_FIXED_EFFECTS
      )),
      data = sample, cluster = ~GID_1
    )
    climate_coefs[[length(climate_coefs) + 1L]] <-
      tidy_fixest(fit_c, label = y) %>%
      rename(outcome = model) %>%
      mutate(econ = econ, .before = 1)
    climate_sums[[length(climate_sums) + 1L]] <- bind_rows(
      so2_sign_sums(fit_c, SO2_LAGS, y),
      so2_sign_sums(fit_c, 1:3, y),
      wald_row(fit_c, "^dlSO2_", label = "all SO2 terms = 0") %>%
        transmute(model = y, lags = "0-3", test, statistic = F_stat, p_value)
    ) %>%
      rename(outcome = model) %>%
      mutate(
        econ = econ, outcome_sd = stats::sd(sample[[y]]),
        within_r2 = fixest::fitstat(fit_c, "wr2")[[1]], .before = 1
      )
  }

  for (alt in c("M1", "M2")) {
    stacked <- fit_stacked(fits$M0, fits[[alt]], sample)
    for (x in list(c("warming", "zTMp_f"), c("cooling", "zTMm_f"))) {
      cmp <- compare_run_effects(stacked, fits$M0, x[[1]], x[[2]],
                                 paste(alt, "- M0"))
      diff_positions[[length(diff_positions) + 1L]] <-
        mutate(cmp$positions, econ = econ, .before = 1)
      diff_joint[[length(diff_joint) + 1L]] <-
        mutate(cmp$joint, econ = econ, .before = 1)
    }
  }

  old_width <- options(width = 250)
  tables[[econ]] <- paste(
    utils::capture.output(print(fixest::etable(
      fits, fitstat = ~n + r2 + wr2 + aic + bic
    ))),
    collapse = "\n"
  )
  options(old_width)
}

sample_summary <- bind_rows(samples)
coefficients <- bind_rows(coefficients)
run_effects <- bind_rows(run_effects)
so2_tests <- bind_rows(so2_tests)
so2_sums <- bind_rows(so2_sums)
climate_coefs <- bind_rows(climate_coefs)
climate_sums <- bind_rows(climate_sums)
diff_positions <- bind_rows(diff_positions)
diff_joint <- bind_rows(diff_joint)

write_out(sample_summary, "sample_summary")
write_out(coefficients, "coefficients")
write_out(run_effects, "run_position_effects")
write_out(so2_tests, "so2_tests")
write_out(so2_sums, "so2_sign_sums")
write_out(climate_coefs, "climate_on_so2_coefficients")
write_out(climate_sums, "climate_on_so2_sums")
write_out(diff_positions, "run_effect_differences")
write_out(diff_joint, "run_effect_differences_joint")
writeLines(
  unlist(lapply(names(tables), function(k) c(paste("=====", k, "====="), tables[[k]], ""))),
  file.path(OUTPUT_DIR, "etable.txt")
)

model_levels <- c(
  M0 = "M0: m_bool",
  M1 = "M1: + SO2 increases, decreases",
  M2 = "M2: + SO2 increases, decreases, lags 0-3"
)
plot_data <- run_effects %>%
  mutate(
    model = factor(model_levels[model], levels = model_levels),
    run = factor(run, levels = c("warming", "cooling"),
                 labels = c("Warm anomaly in a warming run (zTMp_f^2)",
                            "Cold anomaly in a cooling run (zTMm_f^2)")),
    econ = factor(econ, levels = ECON_DATASETS)
  )

p <- ggplot(plot_data, aes(position, estimate, colour = model, shape = model)) +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(ymin = conf_low, ymax = conf_high),
                 position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.5), size = 2.4) +
  # Free y per panel: the thin cooling 6+ cell would otherwise flatten the
  # warming panels.
  facet_wrap(vars(econ, run), ncol = 2, scales = "free_y",
             labeller = label_wrap_gen(multi_line = FALSE)) +
  scale_x_continuous(breaks = seq_len(RUN_LENGTH_CAP),
                     labels = c(seq_len(RUN_LENGTH_CAP - 1L),
                                paste0(RUN_LENGTH_CAP, "+"))) +
  scale_colour_manual(values = grid_colours(model_levels), name = NULL) +
  scale_shape_manual(values = grid_shapes(model_levels), name = NULL) +
  labs(
    title = "Squared-anomaly effect by run position, with and without SO2 controls",
    subtitle = paste0(
      "National panels, ERA5 area weights; effect on growth of a 1 SD anomaly ",
      "(main + position interaction); 95% CI clustered by country"
    ),
    x = "Position in run (1 = run start or no run)",
    y = "Effect on GDP per capita growth"
  ) +
  theme_results()
ggsave(file.path(OUTPUT_DIR, "run_position_effects.png"), p,
       width = 10, height = 6.5, dpi = 200)



# Climate outcomes on SO2 changes, per lag and summed over lags 0-3, scaled to
# a 0.1 log-point (about 10%) increase or decrease in emissions.
step <- 0.1
sign_levels <- c(pos = "10% increase in SO2", neg = "10% decrease in SO2")
lag_levels <- c(paste("Lag", SO2_LAGS), "Sum 0-3")
scale_sign <- function(d, sign) {
  s <- if (sign == "pos") step else -step
  d %>% mutate(
    effect = s * estimate,
    low = pmin(s * conf_low, s * conf_high),
    high = pmax(s * conf_low, s * conf_high),
    sign = factor(sign_levels[[sign]], levels = sign_levels)
  )
}
climate_plot_data <- bind_rows(lapply(c("pos", "neg"), function(sg) {
  lag_rows <- climate_coefs %>%
    filter(term %in% so2_terms(sg, SO2_LAGS)) %>%
    mutate(lag = paste("Lag", SO2_LAGS[match(term, so2_terms(sg, SO2_LAGS))]))
  sum_rows <- climate_sums %>%
    filter(lags == "0-3",
           test == if (sg == "pos") "increases" else "decreases") %>%
    mutate(lag = "Sum 0-3")
  scale_sign(bind_rows(lag_rows, sum_rows), sg)
})) %>%
  mutate(
    lag = factor(lag, levels = lag_levels),
    outcome = factor(CLIMATE_OUTCOMES[outcome], levels = CLIMATE_OUTCOMES),
    econ = factor(econ, levels = ECON_DATASETS)
  )

p_climate <- ggplot(climate_plot_data,
                    aes(lag, effect, colour = sign, shape = sign)) +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_vline(xintercept = length(SO2_LAGS) + 0.5, colour = "#e6e5df",
             linewidth = 0.3) +
  geom_linerange(aes(ymin = low, ymax = high),
                 position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.5), size = 2.4) +
  facet_grid(outcome ~ econ, scales = "free_y",
             labeller = labeller(outcome = label_wrap_gen(18))) +
  scale_colour_manual(values = grid_colours(sign_levels), name = NULL) +
  scale_shape_manual(values = grid_shapes(sign_levels), name = NULL) +
  labs(
    title = "Climate outcomes and a country's own SO2 changes",
    subtitle = paste0(
      "Lags 0-3 of SO2 increases and decreases jointly; country and year ",
      "fixed effects; 95% CI clustered by country"
    ),
    x = NULL,
    y = "Change in outcome"
  ) +
  theme_results()
ggsave(file.path(OUTPUT_DIR, "climate_on_so2.png"), p_climate,
       width = 9, height = 9, dpi = 200)
message("Results written to ", OUTPUT_DIR)
