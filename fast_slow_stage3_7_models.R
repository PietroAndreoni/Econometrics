# Stages 3-7 of the fast-versus-slow warming design
# (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md): level benchmark, confirmatory
# nonlinear rate model, dynamics, equal-endpoint path contrasts and the
# section 7 decisions, on the DOSE x ERA5 (area-weighted) primary sample.
#
# Primary model (region + year FE + country linear trends, two-way country and
# year clustering):
#   g = f(T) + q(P) + b1+ r+ + b2+ r+^2 + b1- r- + b2- r-^2
# with f, q centred quadratics and r the trailing 5-year OLS slope (C / year).
# Every headline quantity is an equal-endpoint path gap C %*% beta (log points,
# negative = fast worse); decision intervals are wild-cluster percentile-t
# (country and year Rademacher DGPs, the wider of the two).
#
# Inputs: fast_slow_stage1_panel.R and fast_slow_stage2_support.R outputs.
# Outputs: results/fast_slow_warming/{models,tables,figures}/.

source("load_functions.R")

design <- fs_design()
design$manifest_hash <- readRDS(fs_path("data", "manifest_hash.rds"))$manifest_hash
panel <- arrow::read_parquet(fs_path("data", "panel_analysis.parquet"))
climate_terms <- arrow::read_parquet(fs_path("data", "climate_terms.parquet"))
spec <- readRDS(fs_path("data", "terms_spec.rds"))
region_climate <- utils::read.csv(fs_path("data", "region_climate.csv"))
scenario <- jsonlite::read_json(fs_path("support", "scenario_choice.json"),
                                simplifyVector = TRUE)
fallback <- utils::read.csv(fs_path("support", "fallback_candidates.csv"))
support_rows <- arrow::read_parquet(
  fs_path("support", "canonical_support_region_years.parquet")
)
sc <- design$scenario
sesoi <- sc$sesoi_log_points
reps <- as.integer(Sys.getenv("FS_BOOT_REPS", design$inference$bootstrap_reps))

d <- panel %>% filter(primary) %>% fs_add_fe_columns()
LEVEL <- c(FS_LEVEL_TERMS, FS_PRECIP_TERMS)
RATE <- c(LEVEL, fs_rate_terms())
is_speed_term <- function(x) grepl("^(r|r[0-9]+|d1|a[0-9]+|rf)_(pos|neg)", x)

boot_w <- list(
  GID_0 = fs_rademacher(unique(d$GID_0), reps,
                        fs_seed(design$hash, "boot-country")),
  year = fs_rademacher(as.character(unique(d$year)), reps,
                       fs_seed(design$hash, "boot-year"))
)

registry <- list()
path_rows <- list()
register <- function(...) registry[[length(registry) + 1]] <<- fs_registry_row(
  ..., design = design
)

# Scenarios and path contrasts ---------------------------------------------------

regions <- region_climate %>%
  filter(GID_1 %in% d$GID_1) %>%
  select(GID_1, B, P)
weights_n <- d %>% count(GID_1, name = "weight")
conf <- scenario$confirmatory
canon <- c(scenario$canonical, list(id = "canonical"))
if (!isTRUE(conf$supported)) {
  message("No supported equal-endpoint contrast: results are model-dependent ",
          "extrapolations and are labelled so.")
  conf <- c(canon[c("total_warming", "fast_years", "slow_years", "horizon")],
            list(id = "canonical_unsupported", supported = FALSE))
}
H <- conf$horizon
K <- design$climate$rate_window
eval_h <- sort(unique(c(10L, conf$slow_years + K - 1L, H, H + K - 1L, 30L)))
profile_h <- seq_len(max(30L, H + K - 1L))

scenario_terms <- function(s, delay = 0, regions_used = regions) {
  fs_path_terms(regions_used, spec, s$total_warming, s$fast_years,
                s$slow_years, horizon = max(profile_h), fast_delay = delay)
}
paths_conf <- scenario_terms(conf)
paths_canon <- if (identical(conf$id, "canonical")) paths_conf else
  scenario_terms(canon)

# Contrast rows total / level / speed for each horizon.
component_contrasts <- function(pt, terms, horizons, weights = weights_n) {
  C <- fs_path_contrast(pt, terms, horizons, weights)
  speed_cols <- is_speed_term(colnames(C))
  level <- C
  level[, speed_cols] <- 0
  speed <- C
  speed[, !speed_cols] <- 0
  out <- rbind(C, level, speed)
  rownames(out) <- paste(rep(c("total", "level", "speed"), each = nrow(C)),
                         rep(horizons, 3), sep = "_h")
  out
}

support_status <- function(s, horizon) {
  if (isTRUE(s$id %in% c("canonical", "canonical_unsupported"))) {
    sh <- fs_support_summary(
      support_rows, weights_n, horizon,
      design$support$min_scenario_weight_supported
    )
    list(status = if (sh$passes && scenario$tail_support_passes)
      "supported" else "outside support", share = sh$min_share)
  } else {
    list(status = if (isTRUE(s$supported)) "supported" else "outside support",
         share = NA_real_)
  }
}

add_path_rows <- function(tab, model_id, scenario_id, s, extra = list()) {
  parts <- strsplit(tab$contrast, "_h")
  tab$component <- vapply(parts, `[`, "", 1)
  tab$horizon <- as.integer(vapply(parts, `[`, "", 2))
  st <- lapply(tab$horizon, function(h) support_status(s, min(h, H)))
  tab$support_status <- vapply(st, `[[`, "", "status")
  tab$min_support_share <- vapply(st, `[[`, numeric(1), "share")
  tab$model_id <- model_id
  tab$scenario_id <- scenario_id
  for (nm in names(extra)) tab[[nm]] <- extra[[nm]]
  path_rows[[length(path_rows) + 1]] <<- tab
  invisible(tab)
}

# Stage 3: level-only benchmark ----------------------------------------------------

level_fit <- fs_fit(LEVEL, d)
level_fit_cy <- fs_fit(LEVEL, d, fe = FS_FIXED_EFFECTS[["country_year"]])
rate_fit <- fs_fit(RATE, d)
stopifnot(identical(fixest::obs(level_fit), fixest::obs(rate_fit)),
          stats::nobs(rate_fit) == nrow(d))            # assertion 9
register("level_broad", "confirmatory", level_fit, d, "none",
         rate_function = "none")
register("level_country_year", "core_robustness", level_fit_cy, d, "none",
         rate_function = "none", fixed_effects = "country_year")

coef_tables <- bind_rows(
  tidy_fixest(level_fit, "level_broad", d, vcov = ~ GID_0 + year),
  tidy_fixest(level_fit_cy, "level_country_year", d, vcov = ~ GID_0 + year),
  tidy_fixest(rate_fit, "rate_broad", d, vcov = ~ GID_0 + year)
)

T_grid <- seq(stats::quantile(d$TM, 0.01), stats::quantile(d$TM, 0.99),
              length.out = 100)
marginal <- bind_rows(lapply(list(broad = level_fit, country_year = level_fit_cy),
                             function(m) {
  G <- cbind(Tc = 1, Tc2 = 2 * (T_grid - spec$T0))
  bind_cols(tibble::tibble(TM = T_grid),
            linear_combination(m, G, df = min(n_distinct(d$GID_0),
                                              n_distinct(d$year)) - 1,
                               vcov = stats::vcov(m, vcov = ~ GID_0 + year)))
}), .id = "fixed_effects")
fs_write(marginal, "tables", "level_marginal_response.csv", design = design)
p_level <- ggplot(marginal, aes(TM, estimate, colour = fixed_effects,
                                fill = fixed_effects)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15,
              colour = NA) +
  geom_line() +
  geom_rug(data = d %>% slice_sample(n = 3000), aes(x = TM),
           inherit.aes = FALSE, alpha = 0.2) +
  labs(x = "Annual temperature (C)",
       y = "Marginal effect on growth (log points per C)",
       colour = "Fixed effects", fill = "Fixed effects",
       caption = "Two-way (country, year) clustered 95% CI; rug: region-years") +
  theme_classic()
ggsave(fs_path("figures", "level_marginal_response.png"), p_level,
       width = 7, height = 4.5, dpi = 150)

# Stage 4: confirmatory rate model -------------------------------------------------

inf_rate <- fs_infer(rate_fit, d, boot_w)
saveRDS(list(delta = lapply(inf_rate$boots, `[[`, "delta"),
             Vflat = lapply(inf_rate$boots, `[[`, "Vflat"),
             coef = inf_rate$prep$coef, V_hat = inf_rate$V_hat),
        fs_path("models", "primary_bootstrap_draws.rds"))
register("rate_broad", "confirmatory", rate_fit, d, "5-year trailing OLS slope")

C_conf <- component_contrasts(paths_conf, RATE, eval_h)
tab_conf <- add_path_rows(fs_contrast_table(inf_rate, C_conf), "rate_broad",
                          conf$id, conf)

# Stage 3 (cont.): residual-calibrated power -------------------------------------
# Level-benchmark residuals are wild-resampled with the same country/year
# Rademacher scheme; the planted incremental speed gap is added through beta_2+
# only. Intervals use the critical values of the primary analysis.

# The planted DGP moves only beta_2+, so the true total and speed gaps are
# equal; power is computed for each component with its own noise and SE.
power_draws <- lapply(boot_w, function(v) v[seq_len(min(1000L, nrow(v))), ,
                                            drop = FALSE])
inf_cal <- fs_infer(rate_fit, d, power_draws, u = stats::residuals(level_fit))
power_cal <- bind_rows(lapply(c("total", "speed"), function(comp) {
  crow <- C_conf[paste0(comp, "_h", H), ]
  crit <- attr(fs_boot_contrasts(inf_rate$prep, inf_rate$boots,
                                 rbind(crow), inf_rate$V_hat),
               "per_dimension")
  q95 <- max(crit$q95)
  q90 <- max(crit$q90)
  cc <- as.numeric(outer(crow, crow))
  bind_rows(lapply(names(inf_cal$boots), function(dim) {
  bt <- inf_cal$boots[[dim]]
  noise <- as.numeric(bt$delta %*% crow)
  se <- sqrt(pmax(as.numeric(bt$Vflat %*% cc), 0))
  bind_rows(lapply(unlist(design$power$planted_speed_gaps), function(target) {
    est <- target + noise
    cls <- fs_classify(est - q95 * se, est + q95 * se, est - q90 * se,
                       est + q90 * se, sesoi, sesoi_confirmed = TRUE)
    tibble::tibble(
      component = comp, boot_dim = dim, planted_gap = target,
      coverage95 = mean(est - q95 * se <= target & est + q95 * se >= target),
      p_fast_worse = mean(cls$direction == "fast worse"),
      p_materially_worse = mean(cls$materiality == "materially fast worse"),
      p_equivalent = mean(cls$materiality == "practically equivalent"),
      mean_se = mean(se), draws = length(est)
    )
  }))
  }))
}))
fs_write(power_cal, "tables", "power_residual_calibrated.csv", design = design)
power_at_sesoi <- power_cal %>%
  filter(planted_gap == -sesoi) %>%
  group_by(component) %>%
  summarise(power = min(p_fast_worse), .groups = "drop")
power_by_component <- stats::setNames(power_at_sesoi$power,
                                      power_at_sesoi$component)
power_ok_by_component <- power_by_component >=
  design$decision$min_power_at_sesoi

# Fixed-effect variants on the identical sample.
fe_variants <- lapply(names(FS_FIXED_EFFECTS), function(fe) {
  data_fe <- if (fe == "zone_year") d %>% filter(!is.na(zone)) else d
  fit_r <- fs_fit(RATE, data_fe, fe = FS_FIXED_EFFECTS[[fe]])
  fit_l <- fs_fit(LEVEL, data_fe, fe = FS_FIXED_EFFECTS[[fe]])
  stopifnot(identical(fixest::obs(fit_r), fixest::obs(fit_l)))
  inf <- if (fe == "broad_primary") inf_rate else fs_infer(fit_r, data_fe, boot_w)
  X <- inf$fwl$X
  Xs <- scale(X)
  ev <- eigen(crossprod(Xs) / nrow(Xs), symmetric = TRUE, only.values = TRUE)$values
  diag_row <- tibble::tibble(
    fixed_effects = fe,
    sample = if (fe == "zone_year") "primary minus regions without a zone" else
      "primary",
    n_obs = stats::nobs(fit_r),
    partial_r2_rate_block = 1 - sum(stats::residuals(fit_r)^2) /
      sum(stats::residuals(fit_l)^2),
    max_vif = max(diag(solve(stats::cor(X)))),
    condition_index = sqrt(max(ev) / min(ev)),
    rate_sd_after_projection = stats::sd(X[, "r_pos"] - X[, "r_neg"]),
    psd_repaired = attr(inf$V_hat, "psd_repaired")
  )
  sample_id <- if (fe == "zone_year") "primary_zone" else "primary"
  tier <- if (fe == "broad_primary") "confirmatory" else "core_robustness"
  register(paste0("rate_", fe), tier, fit_r, data_fe,
           "5-year trailing OLS slope", fixed_effects = fe,
           sample_id = sample_id)
  coef_tab <- tidy_fixest(fit_r, paste0("rate_", fe), data_fe,
                          vcov = ~ GID_0 + year)
  vc <- as.data.frame(as.table(inf$V_hat)) %>%
    mutate(fixed_effects = fe)
  if (fe != "broad_primary") {
    add_path_rows(fs_contrast_table(inf, C_conf), paste0("rate_", fe),
                  conf$id, conf)
  }
  list(diag = diag_row, coef = coef_tab, vcov = vc)
})
fs_write(bind_rows(lapply(fe_variants, `[[`, "diag")) %>%
           left_join(utils::read.csv(fs_path("support",
                                             "rate_identifying_variation.csv")) %>%
                       select(fixed_effects, meaningful_countries,
                              weakly_identified),
                     by = "fixed_effects"),
         "tables", "fe_variant_diagnostics.csv", design = design)
coef_tables <- bind_rows(coef_tables, lapply(fe_variants, `[[`, "coef"))
fs_write(bind_rows(lapply(fe_variants, `[[`, "vcov")), "models",
         "rate_model_vcov_twoway.csv", design = design)

# Stage 6: equal-endpoint paths -----------------------------------------------------

# Canonical contrast (confirmatory if supported, else an extrapolative policy
# scenario).
if (!identical(conf$id, "canonical")) {
  C_canon <- component_contrasts(paths_canon, RATE, eval_h)
  add_path_rows(fs_contrast_table(inf_rate, C_canon), "rate_broad",
                "canonical_extrapolative", canon)
}

# Horizon profile with simultaneous (sup-t) bands, per component.
profile <- bind_rows(lapply(c("total", "level", "speed"), function(comp) {
  C <- component_contrasts(paths_conf, RATE, profile_h)
  C <- C[grepl(paste0("^", comp, "_"), rownames(C)), , drop = FALSE]
  fs_contrast_table(inf_rate, C, sup_t = TRUE) %>%
    mutate(component = comp, horizon = profile_h)
}))
fs_write(profile, "tables", "path_gap_profile.csv", design = design)
p_profile <- ggplot(profile, aes(horizon, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-sesoi, sesoi), linetype = "dotted",
             colour = "grey50") +
  geom_ribbon(aes(ymin = band95_low, ymax = band95_high), alpha = 0.15) +
  geom_ribbon(aes(ymin = ci95_low, ymax = ci95_high), alpha = 0.25) +
  geom_line() +
  geom_vline(xintercept = c(H, H + K - 1), linetype = "dashed",
             colour = "grey60") +
  facet_wrap(~component, nrow = 1) +
  labs(x = "Years since the start of warming",
       y = "Fast minus slow log GDP pc (log points)",
       caption = paste0(
         "+", conf$total_warming, "C in ", conf$fast_years, " vs ",
         conf$slow_years, " years from each region's 1991-2010 climate. ",
         "Dark: pointwise 95% wild-bootstrap-t; light: sup-t band. ",
         "Dotted: +/- provisional SESOI."
       )) +
  theme_classic()
ggsave(fs_path("figures", "path_gap_profile.png"), p_profile, width = 10,
       height = 4, dpi = 150)

# Estimand weights (section 3.2): equal-region and fixed-baseline-population.
d_eq <- d %>% group_by(GID_1) %>% mutate(w_row = 1 / n()) %>% ungroup() %>%
  mutate(w_row = w_row / mean(w_row))
fit_eq <- fs_fit(RATE, d_eq, weights = ~ w_row)
inf_eq <- fs_infer(fit_eq, d_eq, boot_w)
register("rate_broad_equal_region", "core_robustness", fit_eq, d_eq,
         "5-year trailing OLS slope", regression_weight = "1/n_i")
add_path_rows(
  fs_contrast_table(inf_eq, component_contrasts(
    paths_conf, RATE, eval_h,
    weights = tibble::tibble(GID_1 = unique(d$GID_1), weight = 1))),
  "rate_broad_equal_region", conf$id, conf,
  list(estimand = "equal region")
)
d_pop <- d %>% filter(!is.na(pop_base)) %>% mutate(w_row = pop_base / mean(pop_base))
fit_pop <- fs_fit(RATE, d_pop, weights = ~ w_row)
inf_pop <- fs_infer(fit_pop, d_pop, boot_w)
register("rate_broad_population", "core_robustness", fit_pop, d_pop,
         "5-year trailing OLS slope",
         regression_weight = paste0("population ",
                                    design$model$baseline_population_year),
         sample_id = "primary_with_base_population")
w_pop <- d_pop %>% group_by(GID_1) %>%
  summarise(weight = first(pop_base) * n(), .groups = "drop")
add_path_rows(
  fs_contrast_table(inf_pop, component_contrasts(paths_conf, RATE, eval_h,
                                                 weights = w_pop)),
  "rate_broad_population", conf$id, conf,
  list(estimand = "fixed-baseline-population person-years")
)

# Outcome-blind grid of (M, D_F, D_S) with a simultaneous band across cells;
# unsupported cells are shaded and excluded from confirmatory statements.
grid_C <- do.call(rbind, lapply(seq_len(nrow(fallback)), function(i) {
  f <- fallback[i, ]
  pt <- fs_path_terms(regions, spec, f$total_warming, f$fast_years,
                      f$slow_years, horizon = f$horizon)
  C <- component_contrasts(pt, RATE, f$horizon)
  rownames(C) <- paste0(rownames(C), "_M", f$total_warming, "_F",
                        f$fast_years, "_S", f$slow_years)
  C
}))
grid_tab <- bind_rows(lapply(c("total", "speed"), function(comp) {
  Ck <- grid_C[grepl(paste0("^", comp, "_"), rownames(grid_C)), , drop = FALSE]
  fs_contrast_table(inf_rate, Ck, sup_t = TRUE) %>% mutate(component = comp)
})) %>%
  mutate(
    total_warming = as.numeric(sub(".*_M([0-9.]+)_F.*", "\\1", contrast)),
    fast_years = as.integer(sub(".*_F([0-9]+)_S.*", "\\1", contrast)),
    slow_years = as.integer(sub(".*_S([0-9]+)$", "\\1", contrast))
  ) %>%
  left_join(fallback %>% select(total_warming, fast_years, slow_years,
                                horizon, passes, min_supported_share),
            by = c("total_warming", "fast_years", "slow_years"))
fs_write(grid_tab, "tables", "path_grid.csv", design = design)
p_grid <- ggplot(grid_tab, aes(factor(fast_years), estimate,
                               colour = factor(slow_years),
                               alpha = ifelse(passes, "supported",
                                              "outside support"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = band95_low, ymax = band95_high),
                  position = position_dodge(width = 0.6), size = 0.2) +
  scale_alpha_manual(values = c(supported = 1, "outside support" = 0.3)) +
  facet_grid(component ~ total_warming, scales = "free_y",
             labeller = labeller(total_warming = function(x) paste0("+", x, "C"))) +
  labs(x = "Fast duration (years)", y = "Gap at max(20, D_S) (log points)",
       colour = "Slow duration", alpha = NULL,
       caption = "Sup-t 95% bands across the grid (wild bootstrap-t)") +
  theme_classic()
ggsave(fs_path("figures", "path_grid.png"), p_grid, width = 11, height = 6,
       dpi = 150)

# Heterogeneity of the path gap by baseline temperature and predetermined
# income (aggregation weights only: the model is homogeneous).
first_income <- d %>% arrange(GID_1, year) %>% group_by(GID_1) %>%
  summarise(income0 = mean(lgrp_pc_usd[seq_len(min(5, n()))]), .groups = "drop")
groups <- regions %>%
  left_join(first_income, by = "GID_1") %>%
  mutate(baseline_T_tercile = dplyr::ntile(B, 3),
         income_tercile = dplyr::ntile(income0, 3))
for (v in c("baseline_T_tercile", "income_tercile")) {
  for (g in 1:3) {
    ids <- groups$GID_1[groups[[v]] == g]
    tab <- fs_contrast_table(inf_rate, component_contrasts(
      paths_conf, RATE, H, weights = weights_n %>% filter(GID_1 %in% ids)))
    add_path_rows(tab, "rate_broad", conf$id, conf,
                  list(estimand = paste0(v, " = ", g)))
  }
}

# Secondary: pair with (approximately) equal cumulative degree-years, by
# delaying the fast ramp by (D_S - D_F) / 2 years.
delay <- (conf$slow_years - conf$fast_years) / 2
paths_dy <- scenario_terms(conf, delay = delay)
add_path_rows(fs_contrast_table(inf_rate, component_contrasts(paths_dy, RATE,
                                                              eval_h)),
              "rate_broad", paste0(conf$id, "_equal_degree_years"), conf,
              list(estimand = paste0("fast ramp delayed ", delay, " years")))

# Stage 5: dynamics --------------------------------------------------------------------

# Almon cubic distributed lag of annual temperature, lags 0-10.
Lmax <- design$model$distributed_lag_max
deg <- design$model$distributed_lag_almon_degree
almon <- function(df) {
  lag_cols <- c("Tc", paste0("Tc_l", seq_len(Lmax)))
  L <- as.matrix(df[, lag_cols])
  for (j in 0:deg) df[[paste0("Z", j)]] <- as.numeric(L %*% ((0:Lmax)^j))
  df
}
d_dl <- almon(d)
paths_dl <- almon(paths_conf)
DL <- c(paste0("Z", 0:deg), FS_PRECIP_TERMS)
fit_dl <- fs_fit(DL, d_dl)
stopifnot(identical(fixest::obs(fit_dl), fixest::obs(rate_fit)))
inf_dl <- fs_infer(fit_dl, d_dl, boot_w)
register("dl_almon", "core_robustness", fit_dl, d_dl, "none",
         level_function = "Almon cubic lags 0-10 (linear T)",
         rate_function = "implied by lags")
lag_map <- outer(0:Lmax, 0:deg, `^`)
dimnames(lag_map) <- list(paste0("lag", 0:Lmax), paste0("Z", 0:deg))
cum_map <- apply(lag_map, 2, cumsum)
rownames(cum_map) <- paste0("cum", 0:Lmax)
dl_tab <- fs_contrast_table(inf_dl, rbind(lag_map, cum_map), sup_t = FALSE)
fs_write(dl_tab, "tables", "distributed_lag_almon.csv", design = design)
C_dl <- fs_path_contrast(paths_dl, DL, eval_h, weights_n)
rownames(C_dl) <- paste0("total_h", eval_h)
add_path_rows(fs_contrast_table(inf_dl, C_dl), "dl_almon", conf$id, conf)

# Restrictions implied by the linear 5-year slope on unrestricted lags 0-4:
# lag coefficients 1-4 proportional to (0.1, 0, -0.1, -0.2).
RL <- c(FS_LEVEL_TERMS, paste0("Tc_l", 1:4), FS_PRECIP_TERMS)
fit_rl <- fs_fit(RL, d)
R_slope <- rbind(c(Tc_l1 = 0, Tc_l2 = 1, Tc_l3 = 0, Tc_l4 = 0),
                 c(1, 0, 1, 0), c(2, 0, 0, 1))
colnames(R_slope) <- paste0("Tc_l", 1:4)
restriction_tests <- bind_rows(
  fs_wald_linear(coef(fit_rl), stats::vcov(fit_rl, vcov = ~ GID_0 + year),
                 R_slope, df2 = min(n_distinct(d$GID_0), n_distinct(d$year)) - 1) %>%
    mutate(test = "lags 1-4 proportional to the 5-year slope weights"),
  fs_wald_linear(coef(fit_rl), stats::vcov(fit_rl, vcov = ~ GID_0 + year),
                 diag(4) %>% `colnames<-`(paste0("Tc_l", 1:4)),
                 df2 = min(n_distinct(d$GID_0), n_distinct(d$year)) - 1) %>%
    mutate(test = "lags 1-4 jointly zero (no memory)")
)
fs_write(restriction_tests, "tables", "slope_restriction_tests.csv",
         design = design)
coef_tables <- bind_rows(coef_tables,
                         tidy_fixest(fit_rl, "unrestricted_lags_0_4", d,
                                     vcov = ~ GID_0 + year))

# Panel local projections: cumulative growth from t-1 to t+h on the climate at
# t, with two lags of growth, temperature, rate and precipitation (primary) and
# without them. Windows never bridge an invalid growth year.
key <- paste(panel$GID_1, panel$year)
g_valid <- ifelse(panel$valid_growth, panel$g, NA_real_)
lookup_g <- function(k) g_valid[match(paste(d$GID_1, d$year + k), key)]
ckey <- paste(climate_terms$GID_1, climate_terms$year)
lookup_c <- function(col, k) climate_terms[[col]][match(paste(d$GID_1,
                                                              d$year - k), ckey)]
d_lp <- d %>%
  mutate(g_l1 = lookup_g(-1), g_l2 = lookup_g(-2),
         rate_l1 = lookup_c("rate", 1), rate_l2 = lookup_c("rate", 2),
         Pz_l1 = lookup_c("Pz", 1), Pz_l2 = lookup_c("Pz", 2))
H_lp <- design$model$local_projection_horizon_sensitivity
cum <- matrix(NA_real_, nrow(d_lp), H_lp + 1)
running <- 0
for (h in 0:H_lp) {
  running <- running + lookup_g(h)
  cum[, h + 1] <- running
}
LP_LAGS <- c("g_l1", "g_l2", "Tc_l1", "Tc_l2", "rate_l1", "rate_l2", "Pz_l1",
             "Pz_l2")
rF <- conf$total_warming / conf$fast_years
rS <- conf$total_warming / conf$slow_years
lp_contrasts <- rbind(
  rate_fast_vs_slow = c(r_pos = rF - rS, r_pos2 = rF^2 - rS^2)
)
lp_rows <- list()
lp_t <- list()
for (version in c("with_lags", "no_lags")) {
  for (h in 0:H_lp) {
    dh <- d_lp %>% mutate(cum_h = cum[, h + 1])
    terms <- if (version == "with_lags") c(RATE, LP_LAGS) else RATE
    dh <- dh[stats::complete.cases(dh[, c("cum_h", terms)]), ]
    fit_h <- fs_fit(terms, dh, outcome = "cum_h")
    inf_h <- fs_infer(fit_h, dh, boot_w)
    for (cn in c("rate_fast_vs_slow", "temperature_at_T0")) {
      crow <- if (cn == "temperature_at_T0") c(Tc = 1) else
        lp_contrasts[cn, ]
      bt <- fs_boot_tstats(inf_h, crow)
      lp_t[[paste(version, cn, h)]] <- bt$tstar
      lp_rows[[length(lp_rows) + 1]] <- tibble::tibble(
        version, contrast = cn, horizon = h, estimate = bt$estimate,
        se = bt$se, q95 = max(apply(bt$tstar, 2, stats::quantile, 0.95,
                                    na.rm = TRUE)),
        n_obs = stats::nobs(fit_h), n_countries = n_distinct(dh$GID_0)
      )
    }
  }
}
lp <- bind_rows(lp_rows) %>%
  group_by(version, contrast) %>%
  mutate(in_primary_range = horizon <= design$model$local_projection_horizons)
sup_crit <- lp %>%
  distinct(version, contrast) %>%
  rowwise() %>%
  mutate(q_sup95 = {
    hs <- 0:design$model$local_projection_horizons
    mats <- lapply(hs, function(h) lp_t[[paste(version, contrast, h)]])
    max(vapply(seq_len(ncol(mats[[1]])), function(j) {
      stats::quantile(do.call(pmax, c(lapply(mats, function(m) m[, j]),
                                      na.rm = TRUE)), 0.95,
                      na.rm = TRUE)
    }, numeric(1)))
  }) %>%
  ungroup()
lp <- lp %>%
  ungroup() %>%
  left_join(sup_crit, by = c("version", "contrast")) %>%
  mutate(ci95_low = estimate - q95 * se, ci95_high = estimate + q95 * se,
         band95_low = ifelse(in_primary_range, estimate - q_sup95 * se, NA),
         band95_high = ifelse(in_primary_range, estimate + q_sup95 * se, NA))
fs_write(lp, "tables", "local_projections.csv", design = design)
p_lp <- ggplot(lp, aes(horizon, estimate, colour = version, fill = version)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = band95_low, ymax = band95_high), alpha = 0.12,
              colour = NA) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high), width = 0.2,
                position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) +
  facet_wrap(~contrast, scales = "free_y",
             labeller = labeller(contrast = c(
               rate_fast_vs_slow = "5-year slope at fast vs slow ramp rate",
               temperature_at_T0 = "Temperature at T0 (+1C)"))) +
  labs(x = "Horizon h (years)",
       y = "100 x (log y[t+h] - log y[t-1])",
       caption = paste("Horizon-specific conditional projections, not causal",
                       "impulse responses. Bars: pointwise 95% wild",
                       "bootstrap-t; ribbons: sup-t band over h = 0-10.")) +
  theme_classic()
ggsave(fs_path("figures", "local_projections.png"), p_lp, width = 10,
       height = 4.5, dpi = 150)

# Long differences on non-overlapping 10- and 20-year windows (adaptation
# benchmark; does not identify speed).
long_diff <- bind_rows(lapply(unlist(design$model$long_difference_windows),
                              function(L) {
  ld <- d %>%
    mutate(period = (year - 1960L) %/% L) %>%
    group_by(GID_0, GID_1, period) %>%
    filter(n() == L) %>%
    summarise(g = mean(g), Tc = mean(Tc), Tc2 = mean(Tc2), Pz = mean(Pz),
              Pz2 = mean(Pz2), .groups = "drop") %>%
    group_by(GID_1) %>%
    filter(n() >= 2) %>%
    ungroup()
  if (n_distinct(ld$GID_1) < 30) {
    return(tibble::tibble(model = paste0("long_difference_", L),
                          term = "not estimable", window = L,
                          regions = n_distinct(ld$GID_1)))
  }
  fit <- fixest::feols(g ~ Tc + Tc2 + Pz + Pz2 | GID_1 + period, ld,
                       cluster = ~ GID_0, notes = FALSE)
  tidy_fixest(fit, paste0("long_difference_", L)) %>%
    mutate(window = L, regions = n_distinct(ld$GID_1),
           countries = n_distinct(ld$GID_0))
}))
long_diff <- bind_rows(
  long_diff,
  tidy_fixest(level_fit, "annual_level_broad", d, vcov = ~ GID_0 + year) %>%
    mutate(window = 1L)
)
fs_write(long_diff, "tables", "long_differences.csv", design = design)

# Stage 7: climate-only episodes (disabled unless predeclared) ------------------------

if (!file.exists(file.path("config", "fast_slow_episodes.yml"))) {
  fs_stage_status("7_episodes", "SKIPPED_PREDECLARATION_MISSING",
                  "config/fast_slow_episodes.yml was not supplied before outcome access",
                  design = design)
}

# Section 7: decisions ----------------------------------------------------------------

path_contrasts <- bind_rows(path_rows) %>%
  mutate(
    supported = support_status == "supported",
    sesoi = sesoi
  )
cls <- fs_classify(path_contrasts$ci95_low, path_contrasts$ci95_high,
                   path_contrasts$ci90_low, path_contrasts$ci90_high, sesoi,
                   sesoi_confirmed = isTRUE(sc$sesoi_confirmed),
                   supported = path_contrasts$supported,
                   power_ok = dplyr::coalesce(
                     power_ok_by_component[path_contrasts$component],
                     FALSE))
path_contrasts <- bind_cols(path_contrasts, cls) %>%
  mutate(classification = paste(direction, materiality, sep = " | "),
         estimate_log_points = estimate, ci_level = 0.95,
         seed = fs_seed(design$hash, "boot-country"),
         design_hash = substr(design$hash, 1, 16))
fs_write(path_contrasts, "models", "path_contrasts.csv")

# SESOI sensitivity for the confirmatory contrasts.
headline <- path_contrasts %>%
  filter(model_id == "rate_broad", scenario_id == conf$id,
         is.na(estimand),
         component %in% c("total", "level", "speed"))
sesoi_sens <- bind_rows(lapply(c(sesoi, unlist(sc$sesoi_alternatives)),
                               function(s) {
  bind_cols(headline %>% select(component, horizon, estimate, ci95_low,
                                ci95_high, ci90_low, ci90_high, supported),
            fs_classify(headline$ci95_low, headline$ci95_high,
                        headline$ci90_low, headline$ci90_high, s,
                        sesoi_confirmed = FALSE,
                        supported = headline$supported,
                        power_ok = dplyr::coalesce(
                          power_ok_by_component[headline$component],
                          FALSE))) %>%
    mutate(sesoi = s)
}))
fs_write(sesoi_sens, "tables", "decision_sesoi_sensitivity.csv",
         design = design)

fs_write(bind_rows(registry), "models", "model_registry.csv")
fs_write(coef_tables, "tables", "coefficients.csv", design = design)
saveRDS(list(power_ok = power_ok_by_component,
             power_at_sesoi = power_by_component,
             conf = conf, H = H, eval_h = eval_h),
        fs_path("models", "stage3_7_summary.rds"))

fs_stage_status(
  "3_7_models",
  if (!isTRUE(conf$supported) || !all(power_ok_by_component))
    "PASS_WITH_WARNING" else "PASS",
  c(sprintf("residual-calibrated power at the SESOI: total %.2f, speed %.2f",
            power_by_component[["total"]], power_by_component[["speed"]]),
    sprintf("bootstrap draws per dimension: %d", reps),
    if (!isTRUE(conf$supported)) "no supported equal-endpoint contrast"),
  outputs = fs_path("models", c("path_contrasts.csv", "model_registry.csv")),
  design = design
)
