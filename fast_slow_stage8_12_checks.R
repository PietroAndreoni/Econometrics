# Stages 8-12 of the fast-versus-slow warming design
# (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md) on the preferred data (DOSE x
# ERA5, area weighted): the competing trailing-normal anomaly model, mechanisms
# and heterogeneity, the robustness universe within ERA5/DOSE, inference
# variants, and falsification. The confirmatory contrast is the one chosen in
# Stage 2 (support/scenario_choice.json).
#
# Stage 11 (PWT/WDI triangulation) and the climate-product and weighting parts
# of Stage 10 are outside the preferred data and are recorded as skipped.
#
# Outputs: results/fast_slow_warming/{tables,figures,models}/.

source("load_functions.R")

design <- fs_design()
design$manifest_hash <- readRDS(fs_path("data", "manifest_hash.rds"))$manifest_hash
panel <- arrow::read_parquet(fs_path("data", "panel_analysis.parquet"))
climate_terms <- arrow::read_parquet(fs_path("data", "climate_terms.parquet"))
spec <- readRDS(fs_path("data", "terms_spec.rds"))
region_climate <- utils::read.csv(fs_path("data", "region_climate.csv"))
centroids <- utils::read.csv(fs_path("data", "region_centroids.csv"))
s37 <- readRDS(fs_path("models", "stage3_7_summary.rds"))
conf <- s37$conf
H <- s37$H
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
regions <- region_climate %>% filter(GID_1 %in% d$GID_1) %>% select(GID_1, B, P)
weights_n <- d %>% count(GID_1, name = "weight")
paths <- fs_path_terms(regions, spec, conf$total_warming, conf$fast_years,
                       conf$slow_years, horizon = max(30L, H + 4L))
paths_canon <- fs_path_terms(regions, spec, sc$total_warming_c, sc$fast_years,
                             sc$slow_years, horizon = max(30L, H + 4L))

# Total and speed contrasts at H for any term set; `rename` maps model terms
# onto path columns (e.g. the future-rate placebo onto the current rate).
contrast_rows <- function(terms, pt = paths, h = H, weights = weights_n,
                          rename = NULL) {
  path_cols <- if (is.null(rename)) terms else
    ifelse(terms %in% names(rename), rename[terms], terms)
  C <- fs_path_contrast(pt, path_cols, h, weights)
  colnames(C) <- terms
  speed <- C
  speed[, !is_speed_term(path_cols)] <- 0
  out <- rbind(C, speed)
  rownames(out) <- paste0(rep(c("total", "speed"), each = nrow(C)), "_h",
                          rep(h, 2))
  out
}
registry <- list()
register <- function(...) registry[[length(registry) + 1]] <<- fs_registry_row(
  ..., design = design
)
spec_rows <- list()
spec_curve_add <- function(tab, spec_id, tier, group, n_obs) {
  spec_rows[[length(spec_rows) + 1]] <<- tab %>%
    mutate(spec_id = spec_id, tier = tier, group = group, n_obs = n_obs)
}

inf_rate <- fs_infer(fs_fit(RATE, d), d, boot_w)
base_tab <- fs_contrast_table(inf_rate, contrast_rows(RATE))
spec_curve_add(base_tab, "primary", "confirmatory", "primary", nrow(d))

# Stage 8: competing trailing-normal anomaly model ----------------------------------

anomaly_rows <- list()
for (L in c(design$climate$anomaly_normal_years,
            unlist(design$climate$alternative_anomaly_normal_years))) {
  pre <- paste0("a", L)
  d8 <- d %>% filter(!is.na(.data[[paste0("anom", L)]]))
  ANOM <- c(LEVEL, fs_rate_terms(pre))
  fits <- list(level = fs_fit(LEVEL, d8), rate = fs_fit(RATE, d8),
               anomaly = fs_fit(ANOM, d8))
  stopifnot(length(unique(vapply(fits, stats::nobs, 1L))) == 1L)
  p90 <- stats::quantile(d8$rate[d8$rate > 0], 0.90)
  for (m in names(fits)) {
    fit <- fits[[m]]
    cv_country <- cv_demeaned(fit, d8$GID_0[fixest::obs(fit)])
    cv_decade <- cv_demeaned(fit, d8$decade[fixest::obs(fit)])
    cv_tail <- cv_demeaned(fit, d8$rate[fixest::obs(fit)] > p90)
    terms <- switch(m, level = LEVEL, rate = RATE, anomaly = ANOM)
    inf <- fs_infer(fit, d8, boot_w)
    ct <- if (m == "level") {
      fs_contrast_table(inf, contrast_rows(terms)[1, , drop = FALSE])
    } else {
      fs_contrast_table(inf, contrast_rows(terms))
    }
    anomaly_rows[[length(anomaly_rows) + 1]] <- ct %>%
      mutate(model = m, normal_years = L, n_obs = stats::nobs(fit),
             within_adj_r2 = fixest::r2(fit, "war2"),
             rmse_leave_country_out = cv_country$rmse,
             rmse_leave_decade_out = cv_decade$rmse,
             rmse_tail_above_p90_positive_rate = cv_tail$rmse)
    if (L == design$climate$anomaly_normal_years && m == "anomaly") {
      register("anomaly_20", "core_robustness", fit, d8, "none",
               rate_function = "trailing-normal anomaly hinge (20y)",
               sample_id = "primary_with_20y_normal")
    }
  }
}
anomaly_tab <- bind_rows(anomaly_rows)
fs_write(anomaly_tab, "tables", "anomaly_model_comparison.csv", design = design)

# Warming versus cooling: the mirror-image cooling paths through the r- terms.
paths_cool <- fs_path_terms(regions, spec, -conf$total_warming,
                            conf$fast_years, conf$slow_years, horizon = H)
cooling <- fs_contrast_table(inf_rate, contrast_rows(RATE, pt = paths_cool)) %>%
  mutate(direction = "cooling (-M)")
fs_write(bind_rows(base_tab %>% mutate(direction = "warming (+M)"), cooling),
         "tables", "warming_vs_cooling.csv", design = design)
fs_stage_status("8_anomaly", "PASS", character(),
                outputs = fs_path("tables", "anomaly_model_comparison.csv"),
                design = design)

# Stage 9: mechanisms and heterogeneity (secondary; run after Stages 3-7) ----------

dose_raw <- read_dose(data_file("DOSE_V2.14.csv"))
sector_growth <- dose_raw %>%
  filter(!GID_1 %in% c(" ", "")) %>%
  select(GID_1, year, ag = ag_grp_pc_lcu_2015, man = man_grp_pc_lcu_2015,
         serv = serv_grp_pc_lcu_2015) %>%
  group_by(GID_1) %>%
  mutate(across(c(ag, man, serv), function(x) {
    lx <- ifelse(!is.na(x) & x > 0, log(x), NA_real_)
    100 * (lx - lx[match(year - 1L, year)])
  }, .names = "g_{.col}")) %>%
  ungroup() %>%
  select(GID_1, year, starts_with("g_"))
d9 <- d %>% left_join(sector_growth, by = c("GID_1", "year"))
sector_tab <- bind_rows(lapply(c("g_ag", "g_man", "g_serv"), function(y) {
  ds <- d9 %>% filter(!is.na(.data[[y]]))
  fit <- fs_fit(RATE, ds, outcome = y)
  inf <- fs_infer(fit, ds, boot_w)
  register(paste0("rate_sector_", y), "exploratory", fit, ds,
           "5-year trailing OLS slope",
           outcome = paste0("100*dlog ", sub("g_", "", y),
                            "_grp_pc_lcu_2015"),
           sample_id = paste0("primary_with_", y))
  fs_contrast_table(inf, contrast_rows(RATE,
                                       weights = ds %>% count(GID_1, name = "weight"))) %>%
    mutate(outcome = y, n_obs = nrow(ds), n_regions = n_distinct(ds$GID_1),
           n_countries = n_distinct(ds$GID_0))
}))
fs_write(sector_tab, "tables", "mechanisms_sectors.csv", design = design)

# Continuous heterogeneity of the rate response by predetermined moderators:
# 1951-1980 mean temperature, and mean log income in the region's first five
# sample years (predetermined only relative to later years).
first_income <- d %>% arrange(GID_1, year) %>% group_by(GID_1) %>%
  summarise(income0 = mean(lgrp_pc_usd[seq_len(min(5, n()))]), .groups = "drop")
moderators <- region_climate %>%
  select(GID_1, T_pre) %>%
  inner_join(first_income, by = "GID_1") %>%
  filter(GID_1 %in% d$GID_1)
het_rows <- list()
het_tests <- list()
for (mod in c("T_pre", "income0")) {
  dh <- d %>% select(-any_of(c("T_pre", "income0"))) %>%
    left_join(moderators, by = "GID_1") %>%
    filter(!is.na(.data[[mod]])) %>%
    mutate(m_c = .data[[mod]] - mean(moderators[[mod]]))
  inter <- paste0(fs_rate_terms(), "_x_m")
  for (j in seq_along(inter)) dh[[inter[j]]] <- dh[[fs_rate_terms()[j]]] * dh$m_c
  fit <- fs_fit(c(RATE, inter), dh)
  inf <- fs_infer(fit, dh, boot_w)
  het_tests[[mod]] <- tidy_fixest(fit, paste0("heterogeneity_", mod), dh,
                                  vcov = ~ GID_0 + year) %>%
    filter(term %in% inter) %>%
    mutate(moderator = mod)
  base <- contrast_rows(RATE)
  for (q in c(0.10, 0.50, 0.90)) {
    mval <- stats::quantile(moderators[[mod]], q) - mean(moderators[[mod]])
    C <- cbind(base, matrix(0, nrow(base), length(inter),
                            dimnames = list(NULL, inter)))
    C[, inter] <- base[, fs_rate_terms()] * mval
    het_rows[[length(het_rows) + 1]] <- fs_contrast_table(inf, C) %>%
      mutate(moderator = mod, quantile = q,
             value = stats::quantile(moderators[[mod]], q))
  }
}
het_tests <- bind_rows(het_tests) %>%
  mutate(p_bh = stats::p.adjust(p_value, "BH"))
fs_write(het_tests, "tables", "heterogeneity_interactions.csv", design = design)
fs_write(bind_rows(het_rows), "tables", "heterogeneity_path_gaps.csv",
         design = design)
fs_stage_status("9_mechanisms", "PASS",
                "secondary: explains, does not rescue, the primary test",
                design = design)

# Stage 10: robustness universe within ERA5 / DOSE ----------------------------------

robust_specs <- list(
  list(id = "rate_window_3", tier = "core_robustness", group = "rate window",
       terms = c(LEVEL, fs_rate_terms("r3"))),
  list(id = "rate_window_10", tier = "core_robustness", group = "rate window",
       terms = c(LEVEL, fs_rate_terms("r10"))),
  list(id = "acute_one_year_change", tier = "exploratory",
       group = "rate window", terms = c(LEVEL, fs_rate_terms("d1"))),
  list(id = "temperature_spline", tier = "core_robustness",
       group = "temperature response",
       terms = c(grep("^Tns", names(d), value = TRUE), FS_PRECIP_TERMS,
                 fs_rate_terms()))
)
for (s in robust_specs) {
  ds <- d[stats::complete.cases(d[, s$terms]), ]
  fit <- fs_fit(s$terms, ds)
  inf <- fs_infer(fit, ds, boot_w)
  register(s$id, s$tier, fit, ds,
           switch(s$id, rate_window_3 = "3-year slope",
                  rate_window_10 = "10-year slope",
                  acute_one_year_change = "one-year change",
                  "5-year trailing OLS slope"),
           level_function = if (s$id == "temperature_spline")
             "natural spline, fixed knots" else "centered_quadratic")
  spec_curve_add(fs_contrast_table(inf, contrast_rows(s$terms)), s$id, s$tier,
                 s$group, nrow(ds))
}

# Fixed-effect designs and estimand weights from Stages 4 and 6.
pc <- utils::read.csv(fs_path("models", "path_contrasts.csv"))
from_stage6 <- pc %>%
  filter(scenario_id == conf$id, horizon == H,
         component %in% c("total", "speed"),
         model_id != "rate_broad" | !is.na(estimand),
         model_id != "dl_almon",
         is.na(estimand) | grepl("equal region|population", estimand)) %>%
  mutate(spec_id = ifelse(is.na(estimand), model_id,
                          paste(model_id, estimand)),
         tier = "core_robustness",
         group = ifelse(is.na(estimand), "fixed effects", "estimand weights"),
         n_obs = NA_integer_,
         contrast = paste0(component, "_h", horizon)) %>%
  select(contrast, estimate, se, ci95_low, ci95_high, ci90_low, ci90_high,
         spec_id, tier, group, n_obs)

# Outcome treatment: 1/99 winsorized outcome and FWL-Huber (k = 1.345) with a
# country-cluster pairs bootstrap of the full projection-and-fit.
q <- stats::quantile(d$g, c(0.01, 0.99))
d_w <- d %>% mutate(g_w = pmin(pmax(g, q[[1]]), q[[2]]))
fit_w <- fs_fit(RATE, d_w, outcome = "g_w")
spec_curve_add(fs_contrast_table(fs_infer(fit_w, d_w, boot_w),
                                 contrast_rows(RATE)),
               "winsorized_1_99", "core_robustness", "outcome", nrow(d_w))
register("winsorized_1_99", "core_robustness", fit_w, d_w,
         "5-year trailing OLS slope", outcome = "100*dlog, winsorized 1/99")

huber_fit <- function(data) {
  fe_df <- data.frame(GID_1 = data$GID_1, year = data$year, GID_0 = data$GID_0)
  M <- fixest::demean(as.matrix(data[, c("g", RATE)]), f = fe_df,
                      slope.vars = data.frame(year_c = data$year_c),
                      slope.flag = c(0L, 0L, -1L), notes = FALSE)
  fit <- MASS::rlm(M[, -1], M[, 1], k = 1.345, maxit = 100)
  stats::coef(fit)
}
C_rows <- contrast_rows(RATE)
huber_est <- as.numeric(C_rows %*% huber_fit(d)[RATE])
countries <- unique(d$GID_0)
set.seed(fs_seed(design$hash, "huber-bootstrap"))
n_huber <- as.integer(Sys.getenv("FS_HUBER_REPS",
                                 design$inference$robust_bootstrap_reps))
huber_draws <- t(vapply(seq_len(n_huber), function(b) {
  drawn <- sample(countries, length(countries), replace = TRUE)
  db <- bind_rows(lapply(seq_along(drawn), function(i) {
    d[d$GID_0 == drawn[i], ] %>%
      mutate(GID_0 = paste0(GID_0, "_", i), GID_1 = paste0(GID_1, "_", i))
  }))
  as.numeric(C_rows %*% huber_fit(db)[RATE])
}, numeric(nrow(C_rows))))
huber_se <- apply(huber_draws, 2, stats::sd)
spec_curve_add(tibble::tibble(
  contrast = rownames(C_rows), estimate = huber_est, se = huber_se,
  ci95_low = apply(huber_draws, 2, stats::quantile, 0.025),
  ci95_high = apply(huber_draws, 2, stats::quantile, 0.975),
  ci90_low = apply(huber_draws, 2, stats::quantile, 0.05),
  ci90_high = apply(huber_draws, 2, stats::quantile, 0.95)
), "huber_fwl", "core_robustness", "outcome", nrow(d))

# StructChange codes 1 and/or 3 added back (named sensitivities).
for (smp in c("sample_sc1", "sample_sc3", "sample_sc13")) {
  ds <- panel %>% filter(.data[[smp]]) %>% fs_add_fe_columns()
  bw <- list(
    GID_0 = fs_rademacher(unique(ds$GID_0), reps,
                          fs_seed(design$hash, "boot-country")),
    year = fs_rademacher(as.character(unique(ds$year)), reps,
                         fs_seed(design$hash, "boot-year"))
  )
  fit <- fs_fit(RATE, ds)
  register(paste0("rate_", smp), "core_robustness", fit, ds,
           "5-year trailing OLS slope", sample_id = smp)
  spec_curve_add(fs_contrast_table(fs_infer(fit, ds, bw), contrast_rows(RATE,
    weights = ds %>% count(GID_1, name = "weight"))),
    paste0("struct_change_", sub("sample_sc", "codes_", smp), "_added_back"),
    "core_robustness", "sample", nrow(ds))
}

spec_curve <- bind_rows(bind_rows(spec_rows), from_stage6) %>%
  mutate(component = sub("_h.*", "", contrast))
fs_write(spec_curve, "tables", "specification_curve.csv", design = design)
p_spec <- ggplot(spec_curve %>% filter(component %in% c("total", "speed")),
                 aes(reorder(spec_id, estimate), estimate, colour = tier)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-sesoi, sesoi), linetype = "dotted",
             colour = "grey50") +
  geom_pointrange(aes(ymin = ci95_low, ymax = ci95_high), size = 0.2) +
  coord_flip() +
  facet_wrap(~component, scales = "free_x") +
  labs(x = NULL, y = paste0("Fast minus slow at year ", H, " (log points)"),
       colour = "Tier",
       caption = paste0("+", conf$total_warming, "C in ", conf$fast_years,
                        " vs ", conf$slow_years, " years. 95% intervals.")) +
  theme_classic()
ggsave(fs_path("figures", "specification_curve.png"), p_spec, width = 10,
       height = 6, dpi = 150)

# Inference variants for the primary contrasts.
fit_primary <- inf_rate$fit
fwl_sp <- fs_fwl(fit_primary, d, clusters = c("GID_0", "year", "GID_1"))
inference_tab <- bind_rows(
  fs_analytic_contrasts(fit_primary$coefficients,
                        stats::vcov(fit_primary, vcov = ~ GID_0 + year),
                        C_rows, df = inf_rate$df, label = "CGM two-way country + year"),
  fs_analytic_contrasts(fit_primary$coefficients,
                        stats::vcov(fit_primary, vcov = ~ GID_1 + year),
                        C_rows, df = inf_rate$df, label = "two-way region + year"),
  base_tab %>% select(contrast, estimate, se, ci95_low, ci95_high, ci90_low,
                      ci90_high) %>%
    mutate(inference = "wild bootstrap-t, country and year DGPs (decision)"),
  bind_rows(lapply(unlist(design$inference$conley_distance_km), function(km) {
    bind_rows(lapply(unlist(design$inference$conley_time_lags_years),
                     function(lag) {
      V <- fs_conley_hac(fwl_sp, fwl_sp$GID_1, as.integer(fwl_sp$year),
                         centroids, km, lag)
      fs_analytic_contrasts(fit_primary$coefficients, V, C_rows,
                            label = paste0("Conley ", km, " km, ", lag,
                                           "-year Bartlett"))
    }))
  })),
  tibble::tibble(
    contrast = rownames(C_rows), estimate = as.numeric(C_rows %*%
                                                       fit_primary$coefficients[RATE]),
    se = NA_real_,
    inference = "spatial block bootstrap: not run (see report)"
  )
)
fs_write(inference_tab, "tables", "inference_variants.csv", design = design)
fs_stage_status("10_robustness", "PASS_WITH_WARNING",
                c("climate products (CRU, UDel) and fixed-2000 population weights not run: outside the preferred data (ERA5, area)",
                  "spatial block bootstrap replaced by Conley HAC sensitivities"),
                design = design)
fs_stage_status("11_national", "SKIPPED",
                "PWT/WDI triangulation and the global module are outside the preferred data",
                design = design)

# Stage 12: falsification and negative controls ----------------------------------------

# Leads of temperature and future-only 5-year slope, jointly tested.
LEADS <- c("Tc_f1", "Tc_f2", "rf_pos", "rf_neg")
d_lead <- d[stats::complete.cases(d[, LEADS]), ]
fit_lead <- fs_fit(c(RATE, LEADS), d_lead)
lead_test <- fs_wald_linear(
  coef(fit_lead), stats::vcov(fit_lead, vcov = ~ GID_0 + year),
  diag(length(LEADS)) %>% `colnames<-`(LEADS),
  df2 = inf_rate$df
) %>% mutate(test = "leads jointly zero: T(t+1), T(t+2), future 5-year slope hinges")

# Negative control: the future-only slope in place of the trailing slope, with
# its coefficients mapped onto the same path rate columns.
PLACEBO <- c(LEVEL, fs_rate_terms("rf"))
fit_pl <- fs_fit(PLACEBO, d_lead)
placebo_tab <- fs_contrast_table(
  fs_infer(fit_pl, d_lead, boot_w),
  contrast_rows(PLACEBO, rename = stats::setNames(fs_rate_terms(),
                                                  fs_rate_terms("rf")))
) %>% mutate(check = "future-only 5-year slope as exposure")

# Permutation nulls for the speed gap (fs_rate_permutation()): the literal
# scheme (region histories within climate zone x decade strata) and the
# dependence-consistent country-block time shift; only the rate block moves.
n_perm <- as.integer(Sys.getenv("FS_PERM_REPS",
                                design$inference$permutation_reps))
X_fixed <- inf_rate$fwl$X[, LEVEL]
y_dm <- as.numeric(inf_rate$fwl$X %*% inf_rate$fwl$coef + inf_rate$fwl$u)
speed_row <- C_rows[paste0("speed_h", H), fs_rate_terms()]
obs_stats <- c(speed = sum(speed_row * inf_rate$fwl$coef[fs_rate_terms()]),
               r_pos2 = unname(inf_rate$fwl$coef["r_pos2"]))
perm_seeds <- c(region_within_zone_decade = "permutation",
                country_time_shift = "permutation-country-shift")
permutation_tab <- bind_rows(lapply(
  names(perm_seeds),
  function(scheme) {
    perm_stats <- fs_rate_permutation(d, climate_terms, X_fixed, y_dm,
                                      speed_row, scheme, n_perm,
                                      fs_seed(design$hash, perm_seeds[[scheme]]))
    tibble::tibble(
      scheme = scheme, statistic = names(obs_stats), observed = obs_stats,
      perm_mean = colMeans(perm_stats),
      perm_sd = apply(perm_stats, 2, stats::sd),
      p_two_sided = (1 + colSums(abs(perm_stats) >=
                                   matrix(abs(obs_stats), n_perm, 2,
                                          byrow = TRUE))) / (n_perm + 1),
      permutations = n_perm,
      dependence_consistent = scheme == "country_time_shift"
    )
  }
))

# Leave-one-out: country, continent, decade, year (exact refits), with the
# path contrasts at H.
loo <- function(var) {
  bind_rows(lapply(sort(unique(d[[var]])), function(v) {
    ds <- d[d[[var]] != v, ]
    fit <- fs_fit(RATE, ds)
    fs_analytic_contrasts(fit$coefficients,
                          stats::vcov(fit, vcov = ~ GID_0 + year),
                          contrast_rows(RATE, weights = ds %>% count(GID_1, name = "weight")),
                          df = inf_rate$df) %>%
      mutate(left_out_dimension = var, left_out = as.character(v),
             n_obs = nrow(ds))
  }))
}
loo_tab <- bind_rows(loo("GID_0"), loo("continent"), loo("decade"), loo("year"))
fs_write(loo_tab, "tables", "leave_one_out.csv", design = design)

# One-step (FWL) region influence on the speed and total contrasts.
fwl <- inf_rate$fwl
region_of <- d$GID_1[fwl$obs]
Cm <- C_rows[, colnames(fwl$X)]
influence <- bind_rows(lapply(split(seq_len(fwl$n), region_of), function(i) {
  Xi <- fwl$X[i, , drop = FALSE]
  delta <- -solve(fwl$A - crossprod(Xi), crossprod(Xi, fwl$u[i]))
  tibble::tibble(contrast = rownames(Cm), change = as.numeric(Cm %*% delta))
}), .id = "GID_1") %>%
  group_by(contrast) %>%
  arrange(desc(abs(change)), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup()
fs_write(influence, "tables", "influence_regions.csv", design = design)

falsification <- list(
  leads = lead_test,
  placebo = placebo_tab,
  permutation = permutation_tab,
  loo_summary = loo_tab %>%
    group_by(left_out_dimension, contrast) %>%
    summarise(min_estimate = min(estimate), max_estimate = max(estimate),
              sign_changes = sum(sign(estimate) != sign(base_tab$estimate[
                match(first(contrast), base_tab$contrast)])),
              most_influential = left_out[which.max(abs(estimate -
                base_tab$estimate[match(first(contrast), base_tab$contrast)]))],
              .groups = "drop")
)
fs_write(falsification$leads, "tables", "falsification_leads.csv", design = design)
fs_write(falsification$placebo, "tables", "falsification_future_rate.csv",
         design = design)
fs_write(falsification$permutation, "tables", "falsification_permutation.csv",
         design = design)
fs_write(falsification$loo_summary, "tables", "leave_one_out_summary.csv",
         design = design)
fs_write(bind_rows(registry), "models", "model_registry_stage8_12.csv")

fs_stage_status(
  "12_falsification", "PASS",
  c(sprintf("leads joint p = %.3f", lead_test$p_value),
    sprintf("permutation p (speed gap): region scheme %.3f, country time-shift %.3f",
            permutation_tab$p_two_sided[permutation_tab$statistic == "speed" &
                                          !permutation_tab$dependence_consistent],
            permutation_tab$p_two_sided[permutation_tab$statistic == "speed" &
                                          permutation_tab$dependence_consistent]),
    "no plausibly unaffected non-climate outcome in the preferred data",
    "episode placebo dates not applicable (Stage 7 skipped)"),
  design = design
)
