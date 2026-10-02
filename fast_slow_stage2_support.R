# Stage 2 of the fast-versus-slow warming design: outcome-blind support and
# power audit (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md, section 6).
#
# This script reads only the climate regressors and the sample keys written by
# fast_slow_stage1_panel.R. It never opens a GDP/GRP column: the inputs are
# checked for outcome-like names before anything else runs.
#
# Outputs: results/fast_slow_warming/support/ and qa/stage_status_2_support.json.

source("load_functions.R")

design <- fs_design()
climate_terms <- arrow::read_parquet(fs_path("data", "climate_terms.parquet"))
sample_key <- arrow::read_parquet(fs_path("data", "sample_key.parquet"))
region_climate <- utils::read.csv(fs_path("data", "region_climate.csv"))
spec <- readRDS(fs_path("data", "terms_spec.rds"))
outcome_like <- grep("grp|gdp|dlg|^g$", c(names(climate_terms),
                                           names(sample_key)), value = TRUE)
if (length(outcome_like)) stop("Stage 2 inputs carry outcomes: ", outcome_like)

obs <- sample_key %>%
  filter(primary) %>%
  inner_join(climate_terms %>% select(-GID_0), by = c("GID_1", "year"))
obs <- fs_add_fe_columns(obs)
stopifnot(nrow(obs) == sum(sample_key$primary))

cfg <- design$support
sc <- design$scenario
weights_n <- obs %>% count(GID_1, name = "weight")
regions <- region_climate %>%
  filter(GID_1 %in% obs$GID_1) %>%
  select(GID_0, GID_1, B, P)
stopifnot(!anyNA(regions$B), !anyNA(regions$P))

# Rate distributions -------------------------------------------------------------

rate_cols <- c(d1 = "d1", rate3 = "rate3", rate5 = "rate5", rate10 = "rate10")
probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
rate_quantiles <- bind_rows(lapply(names(rate_cols), function(r) {
  x <- obs[[rate_cols[[r]]]]
  x <- x[!is.na(x)]
  pos <- x[x > 0]
  bind_rows(
    tibble::tibble(rate = r, sign = "signed", q = probs,
                   value = unname(stats::quantile(x, probs)), n = length(x)),
    tibble::tibble(rate = r, sign = "positive", q = probs,
                   value = unname(stats::quantile(pos, probs)), n = length(pos))
  )
}))
fs_write(rate_quantiles, "support", "rate_quantiles.csv", design = design)

long_rates <- obs %>%
  select(GID_1, year, all_of(unname(rate_cols))) %>%
  tidyr::pivot_longer(-c(GID_1, year), names_to = "rate", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(rate = factor(rate, levels = unname(rate_cols),
                       labels = c("1-year change", "3-year slope",
                                  "5-year slope (primary)", "10-year slope")))
ramp_lines <- tibble::tibble(
  path = c("fast", "slow"),
  value = sc$total_warming_c / c(sc$fast_years, sc$slow_years)
)
p_dist <- ggplot(long_rates, aes(value)) +
  geom_histogram(bins = 120, fill = "grey40") +
  geom_vline(data = ramp_lines, aes(xintercept = value, colour = path),
             linetype = "dashed", show.legend = FALSE) +
  scale_colour_manual(values = c(fast = "firebrick", slow = "steelblue")) +
  facet_wrap(~rate, scales = "free") +
  labs(x = "degrees C per year", y = "Region-years",
       caption = paste("Dashed: canonical fast (red) and slow (blue) ramp",
                       "rates. ERA5, area weighted; primary sample.")) +
  theme_classic()
ggsave(fs_path("support", "rate_distributions.png"), p_dist, width = 9,
       height = 6, dpi = 150)

# Positive runs of the 5-year slope.
runs <- obs %>%
  arrange(GID_1, year) %>%
  group_by(GID_1) %>%
  mutate(pos = rate > 0,
         run_id = cumsum(pos != lag(pos, default = FALSE) |
                           year != lag(year, default = -1L) + 1L)) %>%
  filter(pos) %>%
  count(GID_1, run_id, name = "length") %>%
  ungroup()
fs_write(runs %>% count(length, name = "runs"), "support",
         "positive_run_lengths.csv", design = design)

# Tails: region-years at or above the canonical fast rate and above the 90th
# percentile of positive slopes, by country, continent, decade, baseline
# temperature.
p90_pos <- stats::quantile(obs$rate[obs$rate > 0], 0.90)
fast_rate <- sc$total_warming_c / sc$fast_years
slow_rate <- sc$total_warming_c / sc$slow_years
obs <- obs %>%
  left_join(regions %>% select(GID_1, B), by = "GID_1") %>%
  mutate(B_tercile = dplyr::ntile(B, 3))
tail_by <- function(threshold, label) {
  tail <- obs %>% filter(rate >= threshold)
  bind_rows(lapply(c("GID_0", "continent", "decade", "B_tercile"), function(v) {
    tail %>% count(level = as.character(.data[[v]]), name = "region_years") %>%
      mutate(dimension = v, tail = label)
  }))
}
tail_tables <- bind_rows(
  tail_by(fast_rate, "rate5 >= canonical fast rate"),
  tail_by(p90_pos, "rate5 >= p90 of positive slopes")
)
fs_write(tail_tables, "support", "tail_counts_by_group.csv", design = design)

# Joint support of temperature and rate.
p_joint <- ggplot(obs, aes(TM, rate)) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(trans = "log10") +
  geom_hline(data = ramp_lines, aes(yintercept = value, colour = path),
             linetype = "dashed", show.legend = FALSE) +
  scale_colour_manual(values = c(fast = "firebrick", slow = "steelblue")) +
  labs(x = "Annual temperature (C)", y = "5-year slope (C / year)",
       fill = "Region-years") +
  theme_classic()
ggsave(fs_path("support", "joint_support_temperature_rate.png"), p_joint,
       width = 7, height = 5, dpi = 150)

# Residualized (identifying) support ---------------------------------------------------
# The rate after each Stage 4 fixed-effect projection, from climate data only.

residual_rate <- lapply(names(FS_FIXED_EFFECTS), function(fe) {
  d <- obs %>% mutate(zone = ifelse(is.na(zone), -1L, zone)) %>%
    fs_add_fe_columns()
  fit <- fixest::feols(stats::as.formula(paste("rate ~ 1 |",
                                               FS_FIXED_EFFECTS[[fe]])),
                       d, notes = FALSE, fixef.rm = "none")
  res <- as.numeric(stats::residuals(fit))
  by_country <- tibble::tibble(GID_0 = d$GID_0, res = res) %>%
    group_by(GID_0) %>%
    summarise(ss = sum(res^2), sd = stats::sd(res), .groups = "drop") %>%
    mutate(share = ss / sum(ss))
  meaningful <- sum(by_country$sd >= design$model$meaningful_country_rate_sd_share *
                      stats::sd(res))
  list(
    residual = res,
    summary = tibble::tibble(
      fixed_effects = fe,
      raw_sd = stats::sd(obs$rate),
      residual_sd = stats::sd(res),
      residual_share_of_variance = stats::var(res) / stats::var(obs$rate),
      effective_countries = 1 / sum(by_country$share^2),
      meaningful_countries = meaningful,
      weakly_identified = meaningful < design$model$min_meaningful_countries
    )
  )
})
names(residual_rate) <- names(FS_FIXED_EFFECTS)
fe_variation <- bind_rows(lapply(residual_rate, `[[`, "summary"))
fs_write(fe_variation, "support", "rate_identifying_variation.csv",
         design = design)

obs$rate_resid <- residual_rate$broad_primary$residual
contrast_rate <- fast_rate - slow_rate
tail_support <- bind_rows(
  fs_tail_counts(obs, "rate", fast_rate) %>%
    mutate(measure = "raw 5-year slope >= fast ramp rate"),
  fs_tail_counts(obs, "rate", p90_pos) %>%
    mutate(measure = "raw 5-year slope >= p90 of positive slopes"),
  fs_tail_counts(obs, "rate_resid", contrast_rate) %>%
    mutate(measure = "residualized slope >= fast - slow ramp rate")
) %>%
  mutate(passes = region_years >= cfg$min_tail_region_years &
           countries >= cfg$min_tail_countries &
           decades >= cfg$min_tail_decades)
fs_write(tail_support, "support", "tail_support.csv", design = design)

clusters <- tibble::tibble(
  countries = n_distinct(obs$GID_0), years = n_distinct(obs$year),
  regions = n_distinct(obs$GID_1),
  effective_countries_raw_rate = 1 / sum((obs %>% group_by(GID_0) %>%
    summarise(s = sum((rate - mean(obs$rate))^2)) %>%
    mutate(s = s / sum(s)))$s^2)
)
fs_write(clusters, "support", "effective_clusters.csv", design = design)

# Scenario support ------------------------------------------------------------------------

obs_support <- obs %>% select(GID_0, GID_1, year, TM, rate, RR)
loo <- RANN::nn2(
  sweep(as.matrix(obs_support[, c("TM", "rate", "RR")]), 2,
        vapply(c("TM", "rate", "RR"), function(v) stats::sd(obs_support[[v]]),
               numeric(1)), "/"),
  k = cfg$neighbor_k + 1L
)$nn.dists[, cfg$neighbor_k + 1L]

scenario_support <- function(M, DF, DS, H, rate_q = cfg$rate_quantiles,
                             delay = 0) {
  pt <- fs_path_terms(regions %>% select(GID_1, B, P), spec, M, DF, DS,
                      horizon = H, fast_delay = delay)
  scen <- pt %>%
    left_join(regions %>% select(GID_1, GID_0), by = "GID_1") %>%
    select(GID_0, GID_1, path, s, TM, rate, RR)
  sup <- fs_support_scenario(obs_support, scen, k = cfg$neighbor_k,
                             dist_q = cfg$max_distance_quantile,
                             rate_q = unlist(rate_q),
                             min_regions = cfg$min_neighbor_regions,
                             min_countries = cfg$min_neighbor_countries,
                             loo = loo)
  summ <- fs_support_summary(sup, weights_n, H,
                             cfg$min_scenario_weight_supported)
  list(terms = pt, support = sup, summary = summ)
}

H0 <- sc$common_horizon
canonical <- scenario_support(sc$total_warming_c, sc$fast_years, sc$slow_years,
                              max(H0, 30L))
canonical_strict <- scenario_support(sc$total_warming_c, sc$fast_years,
                                     sc$slow_years, max(H0, 30L),
                                     rate_q = cfg$strict_rate_quantiles)
canonical_summary <- fs_support_summary(canonical$support, weights_n, H0,
                                        cfg$min_scenario_weight_supported)
strict_summary <- fs_support_summary(canonical_strict$support, weights_n, H0,
                                     cfg$min_scenario_weight_supported)
fs_write(canonical_summary$by_year %>% mutate(rule = "1-99"), "support",
         "canonical_support_by_year.csv", design = design)
fs_write(strict_summary$by_year %>% mutate(rule = "5-95"), "support",
         "canonical_support_by_year_strict.csv", design = design)
fs_write(canonical$support %>%
           select(GID_1, path, s, TM, rate, RR, rate_ok, temp_ok,
                  knn_distance, knn_regions, knn_countries, knn_ok, supported),
         "support", "canonical_support_region_years.parquet")

# Leverage of scenario points against the centred regressor cloud.
terms_rate <- c(FS_LEVEL_TERMS, FS_PRECIP_TERMS, fs_rate_terms())
Xc <- scale(as.matrix(obs[, terms_rate]), scale = FALSE)
XtX_inv <- solve(crossprod(Xc))
centre <- attr(Xc, "scaled:center")
leverage <- function(M) {
  D <- sweep(as.matrix(M[, terms_rate]), 2, centre)
  rowSums((D %*% XtX_inv) * D) + 1 / nrow(Xc)
}
obs_leverage <- leverage(obs)
scen_leverage <- canonical$terms %>%
  filter(s <= H0) %>%
  mutate(leverage = leverage(.)) %>%
  group_by(path, s) %>%
  summarise(max_leverage = max(leverage), median_leverage = median(leverage),
            .groups = "drop") %>%
  mutate(observed_p99_leverage = stats::quantile(obs_leverage, 0.99),
         observed_max_leverage = max(obs_leverage))
fs_write(scen_leverage, "support", "canonical_leverage.csv", design = design)

p_support <- ggplot(canonical_summary$by_year %>%
                      tidyr::pivot_longer(ends_with("_share")),
                    aes(s, value, colour = name)) +
  geom_line() +
  geom_hline(yintercept = cfg$min_scenario_weight_supported,
             linetype = "dashed") +
  facet_wrap(~path) +
  labs(x = "Years since the start of warming", y = "Weighted share of regions",
       colour = NULL) +
  theme_classic()
ggsave(fs_path("support", "canonical_support_share.png"), p_support,
       width = 8, height = 4, dpi = 150)

# Deterministic fallback grid (always tabulated; used only if the canonical
# pair fails).
grid <- tidyr::crossing(total_warming = unlist(cfg$fallback_total_warming_c),
                        fast_years = unlist(cfg$fallback_fast_years),
                        slow_years = unlist(cfg$fallback_slow_years)) %>%
  filter(fast_years < slow_years) %>%
  mutate(horizon = pmax(H0, slow_years))
grid$min_supported_share <- NA_real_
grid$passes <- NA
for (i in seq_len(nrow(grid))) {
  ss <- scenario_support(grid$total_warming[i], grid$fast_years[i],
                         grid$slow_years[i], grid$horizon[i])
  grid$min_supported_share[i] <- ss$summary$min_share
  grid$passes[i] <- ss$summary$passes
}
tail_ok <- all(tail_support$passes)
grid$passes <- grid$passes & tail_ok
fallback <- fs_fallback_table(
  grid, obs$rate[obs$rate > 0],
  unlist(cfg$fallback_targets_positive_rate_quantiles)
)
fs_write(fallback, "support", "fallback_candidates.csv", design = design)

canonical_passes <- canonical_summary$passes && tail_ok
chosen <- if (canonical_passes) {
  list(id = "canonical", total_warming = sc$total_warming_c,
       fast_years = sc$fast_years, slow_years = sc$slow_years, horizon = H0,
       supported = TRUE)
} else if (any(fallback$passes)) {
  f <- fallback[which(fallback$rank == 1L), ]
  list(id = "within_support_fallback", total_warming = f$total_warming,
       fast_years = f$fast_years, slow_years = f$slow_years,
       horizon = f$horizon, supported = TRUE)
} else {
  list(id = "none_supported", total_warming = NA, fast_years = NA,
       slow_years = NA, horizon = NA, supported = FALSE)
}
scenario_record <- list(
  canonical = list(total_warming = sc$total_warming_c,
                   fast_years = sc$fast_years, slow_years = sc$slow_years,
                   horizon = H0, min_supported_share = canonical_summary$min_share,
                   min_supported_share_strict = strict_summary$min_share,
                   passes = canonical_passes,
                   label = if (canonical_passes) "confirmatory" else
                     "extrapolative policy scenario"),
  confirmatory = chosen,
  tail_support_passes = tail_ok,
  rate_bounds_1_99 = unname(attr(canonical$support, "rate_bounds")),
  knn_distance_threshold = unname(attr(canonical$support, "distance_threshold"))
)
fs_write(scenario_record, "support", "scenario_choice.json")

# Design-range power --------------------------------------------------------------------
# All climate coefficients zero except beta_2+ (planted), or the lag DGP
# sum_k w_k T_{t-k} (misspecified), each scaled to the planted gaps. Noise:
# idiosyncratic + AR(1) country-year shocks + year shocks (the last are
# absorbed exactly by the year effects, so results cannot vary with them).

pw <- design$power
conf <- if (chosen$supported) chosen else scenario_record$canonical
paths_conf <- fs_path_terms(regions %>% select(GID_1, B, P), spec,
                            conf$total_warming, conf$fast_years,
                            conf$slow_years, conf$horizon)
lag_w <- unlist(pw$misspecified_lag_weights)
lag_signal <- function(d) {
  out <- lag_w[1] * d$TM
  for (k in seq_along(lag_w)[-1]) out <- out + lag_w[k] * d[[paste0("TM_l", k - 1)]]
  out
}
paths_conf$lagsig <- lag_signal(paths_conf)
obs$lagsig <- lag_signal(obs)
C_all <- fs_path_contrast(paths_conf, c(terms_rate, "lagsig"), conf$horizon,
                          weights_n)
contrasts <- rbind(
  speed = replace(C_all[1, terms_rate], !terms_rate %in% fs_rate_terms(), 0),
  total = C_all[1, terms_rate]
)
fe_df <- data.frame(GID_1 = obs$GID_1, year = obs$year, GID_0 = obs$GID_0)
X_dm <- fixest::demean(as.matrix(obs[, terms_rate]), f = fe_df,
                       slope.vars = data.frame(year_c = obs$year_c),
                       slope.flag = c(0L, 0L, -1L), notes = FALSE)
sig_dm <- as.numeric(fixest::demean(obs$lagsig, f = fe_df,
                                    slope.vars = data.frame(year_c = obs$year_c),
                                    slope.flag = c(0L, 0L, -1L), notes = FALSE))
noise_grid <- tidyr::crossing(idiosyncratic_sd = unlist(pw$idiosyncratic_sd),
                              country_sd = unlist(pw$country_sd),
                              year_sd = unlist(pw$year_sd),
                              country_ar1 = unlist(pw$country_ar1))
power <- fs_power_design(
  X = X_dm, fe = obs %>% select(GID_1, GID_0, year, year_c),
  contrasts = contrasts, lag_signal_demeaned = sig_dm,
  lag_total = C_all[1, "lagsig"], grid = noise_grid,
  targets = unlist(pw$planted_speed_gaps), sesoi = sc$sesoi_log_points,
  n_sim = as.integer(Sys.getenv("FS_POWER_SIMS",
                                 design$inference$power_simulations)),
  seed = fs_seed(design$hash, "power-design")
)
fs_write(power, "support", "power_design_range.csv", design = design)
saveRDS(list(contrasts = contrasts, lag_total = C_all[1, "lagsig"],
             beta_rpos2_per_log_point = 1 / contrasts["speed", "r_pos2"]),
        fs_path("support", "power_design_objects.rds"))

fs_stage_status(
  "2_support",
  if (!chosen$supported) "STOP" else if (!canonical_passes)
    "PASS_WITH_WARNING" else "PASS",
  c(sprintf("canonical +%.2fC in %d vs %d years: min supported share %.3f (%s)",
            sc$total_warming_c, sc$fast_years, sc$slow_years,
            canonical_summary$min_share,
            if (canonical_passes) "supported" else "outside support"),
    sprintf("confirmatory contrast: %s", chosen$id),
    if (!tail_ok) "tail support rule failed"),
  outputs = fs_path("support", c("scenario_choice.json", "tail_support.csv",
                                 "power_design_range.csv")),
  design = design
)
