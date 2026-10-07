# The decomposed adaptation-shock model, PWT x UDelaware: is dropping the
# unit trends innocent, and do the Hamilton shocks need arms and climate
# interactions?
#
# Model (see test_adaptation_shock_model.R): the 30-year anomaly splits into
# a slow adaptation gap and a Hamilton cycle,
#   zTM = gapTM + cycTM,
#   gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
#   cycTM = (TM - trend_TM) / sd_TM_all,
# and the same for precipitation. Arms: the squared warm and cold gaps (wet and
# dry for precipitation), each interacted with the long-run climate.
#
# (a) Trends. With growth as the outcome, the unit effects absorb each unit's
# average growth and the unit trends (GID_1[year]) a linear change in it, such
# as the slowdown of a catching-up economy. The slow gaps trend too, so
# dropping the unit trends can load growth trends onto the arms. The script
#   - estimates each unit's growth trend (net of year effects) and relates it
#     to its income relative to the other units at entry (convergence), its
#     long-run temperature and the trend of its temperature gap;
#   - measures the within-country convergence speed (growth on lagged log
#     income, unit and year effects);
#   - measures how much of the within variation of each regressor the unit
#     trends absorb;
#   - refits the model under five trend treatments on an identical sample:
#       trends       year + unit effects + unit trends (main_analysis.Rmd)
#       no trends    year + unit effects
#       convergence  no trends + lagged log income per capita
#       income       no trends + year x relative income at entry
#       climate_trends  no trends + year x long-run temperature and precipitation
#
# (b) Shocks. Three forms of the Hamilton shock, each for temperature and
# precipitation:
#       symmetric          cycTM^2
#       symmetric_climate  cycTM^2 and cycTM^2 x long-run climate
#       arms_climate       warm and cold cycles squared (max(cycTM, 0)^2,
#                          max(-cycTM, 0)^2), each x long-run climate
#
# Effects of a 1-SD gap or cycle at long-run climates of 5, 15 and 25 degrees
# Celsius (0.5, 1.5 and 2.5 m), Wald tests of the blocks, and fit. Area and
# population weights, errors clustered by country.
#
# Output: results/gap_shock_trends/ (unit_trends.csv, convergence.csv,
# absorption.csv, model_tests.csv, model_effects.csv, convergence_<weight>.png,
# effects_<weight>.png)

source("load_functions.R")

ECON <- "PWT"
CLIMATE_SOURCE <- "UDelaware"
CLIMATE_WEIGHTS <- c(area = "area", pop = "concurrent population")
ARM_BASELINE <- "lag_30"
ARM_SD <- "all"
# Run option: "full" (default) fits the Hamilton trend over each unit's entire
# climate series; "realtime" fits each year's trend on the 30 years up to the
# forecast origin only (panel_config(hamilton_window = 30)), output in a
# "_realtime30" folder. Rscript test_gap_shock_trends.R realtime
HAMILTON <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(HAMILTON)) HAMILTON <- "full"
stopifnot(HAMILTON %in% c("full", "realtime"))
HAMILTON_WINDOW <- if (HAMILTON == "realtime") 30L else NULL
BASE_YEAR <- 1990
MIN_UNIT_YEARS <- 20
CONF_LEVEL <- 0.95
AT_VALUES <- list(TM = c(5, 15, 25), RR = c(0.5, 1.5, 2.5))
CLIMATE <- c(TM = "mean_TM_all", RR = "mean_RR_all")
DEG_C <- "°C"

TRENDS <- list(
  trends = list(label = "Unit trends", fe = "year + GID_1 + GID_1[year]", controls = NULL),
  notrends = list(label = "No trends", fe = "year + GID_1", controls = NULL),
  convergence = list(label = "No trends + lagged log income", fe = "year + GID_1",
                     controls = "lag_lgdp"),
  income = list(label = "No trends + year × initial income", fe = "year + GID_1",
                controls = "trend_income"),
  climate_trends = list(label = "No trends + year × climate", fe = "year + GID_1",
                 controls = c("trend_TM_climate", "trend_RR_climate"))
)

# Gap arms (quadratic, x climate): variable, climate variable, label.
GAP_ARMS <- tibble::tribble(
  ~var,    ~variable, ~label,
  "gTMp2", "TM",      "Warm gap",
  "gTMm2", "TM",      "Cold gap",
  "gRRp2", "RR",      "Wet gap",
  "gRRm2", "RR",      "Dry gap"
)
# Shock terms of each form: variable, climate interaction (TRUE/FALSE), and the
# shock labels each term contributes to.
SHOCK_FORMS <- list(
  symmetric = tibble::tribble(
    ~var,    ~variable, ~interacted, ~labels,
    "sTM2",  "TM",      FALSE,       c("Warm shock", "Cold shock"),
    "sRR2",  "RR",      FALSE,       c("Wet shock", "Dry shock")
  ),
  symmetric_climate = tibble::tribble(
    ~var,    ~variable, ~interacted, ~labels,
    "sTM2",  "TM",      TRUE,        c("Warm shock", "Cold shock"),
    "sRR2",  "RR",      TRUE,        c("Wet shock", "Dry shock")
  ),
  arms_climate = tibble::tribble(
    ~var,    ~variable, ~interacted, ~labels,
    "sTMp2", "TM",      TRUE,        "Warm shock",
    "sTMm2", "TM",      TRUE,        "Cold shock",
    "sRRp2", "RR",      TRUE,        "Wet shock",
    "sRRm2", "RR",      TRUE,        "Dry shock"
  )
)
SHOCK_LABELS <- c(symmetric = "Symmetric shock",
                  symmetric_climate = "Symmetric shock × climate",
                  arms_climate = "Shock arms × climate")
TERM_LEVELS <- c("Warm gap", "Cold gap", "Wet gap", "Dry gap",
                 "Warm shock", "Cold shock", "Wet shock", "Dry shock")

OUTPUT_DIR <- file.path("results", paste0("gap_shock_trends", if (HAMILTON == "realtime") "_realtime30" else ""))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit <- function(terms, data, fe) {
  fit_panel_model(terms, data, fixed_effects = fe, cluster = ~GID_0)
}
withhold <- function(data, keep_rows) {
  data$dlgrp_pc_usd[setdiff(seq_len(nrow(data)), keep_rows)] <- NA
  data
}
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
model_terms <- function(shocks) {
  gap <- c(GAP_ARMS$var, paste0(GAP_ARMS$var, ":", CLIMATE[GAP_ARMS$variable]))
  shock <- c(shocks$var, paste0(shocks$var, ":", CLIMATE[shocks$variable])[shocks$interacted])
  list(gap = gap, shock = shock)
}
# Wald test of the named coefficients (any order of interaction names).
block_p <- function(model, names) {
  terms <- names(stats::coef(model))
  keep <- unique(unlist(lapply(names, function(x) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    if (length(parts) == 2) interaction_name(parts[[1]], parts[[2]], terms) else x
  })))
  fixest::wald(model, keep = paste0("^(", paste(keep, collapse = "|"), ")$"),
               print = FALSE)$p
}
# Joint test that the warm and cold shock arms are equal (main term and
# climate interaction).
symmetry_p <- function(model, up, down, climate) {
  b <- stats::coef(model)
  V <- stats::vcov(model)
  terms <- names(b)
  R <- matrix(0, 2, length(b), dimnames = list(NULL, terms))
  R[1, c(up, down)] <- c(1, -1)
  R[2, c(interaction_name(up, climate, terms), interaction_name(down, climate, terms))] <- c(1, -1)
  r <- R %*% b
  stat <- as.numeric(t(r) %*% solve(R %*% V %*% t(R)) %*% r)
  stats::pchisq(stat, df = 2, lower.tail = FALSE)
}
# Effect of a 1-SD term at long-run climate m (the term alone if it is not
# interacted).
effect_at <- function(model, var, climate, m) {
  terms <- names(stats::coef(model))
  inter <- intersect(c(paste0(var, ":", climate), paste0(climate, ":", var)), terms)
  weights <- if (length(inter)) stats::setNames(c(1, m), c(var, inter[[1]])) else
    stats::setNames(1, var)
  lincom(model, weights, conf_level = CONF_LEVEL)
}
# Each unit's linear trend in x (net of year effects), per decade.
unit_slopes <- function(data, x) {
  data$x_net <- NA_real_
  m <- fixest::feols(as.formula(paste(x, "~ 1 | year")), data)
  data$x_net[fixest::obs(m)] <- stats::resid(m)
  data %>%
    filter(!is.na(x_net)) %>%
    group_by(GID_1) %>%
    filter(n() >= MIN_UNIT_YEARS) %>%
    summarise(slope = 10 * stats::coef(stats::lm(x_net ~ year))[["year"]], .groups = "drop")
}
# Within variation of x net of the fixed effects.
within_var <- function(x, data, fe) {
  stats::var(stats::resid(fixest::feols(as.formula(paste(x, "~ 1 |", fe)), data,
                                        panel.id = c("GID_1", "year"))))
}

tests <- list()
effects <- list()
trends_out <- list()
convergence <- list()
absorption <- list()
for (w in names(CLIMATE_WEIGHTS)) {
  message("Building the ", w, "-weighted panel")
  cfg <- panel_config(
    econ_year_min = min(load_econ_panel(resolve_econ_source(ECON))$year),
    climate_source = CLIMATE_SOURCE,
    climate_weight = CLIMATE_WEIGHTS[[w]],
    hamilton_window = HAMILTON_WINDOW,
    clim_history_years = 1L
  )
  d <- suppressMessages(build_dat(ECON, config = cfg, baseline = ARM_BASELINE,
                                  sd_baseline = ARM_SD)) %>%
    arrange(GID_1, year) %>%
    group_by(GID_1) %>%
    mutate(lag_lgdp = ifelse(dplyr::lag(year) == year - 1, dplyr::lag(lgrp_pc_usd), NA)) %>%
    ungroup() %>%
    mutate(
      gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
      cycTM = (TM - trend_TM) / sd_TM_all,
      gapRR = (trend_RR - mean_RR_lag_30) / sd_RR_all,
      cycRR = (RR - trend_RR) / sd_RR_all,
      gTMp2 = pmax(gapTM, 0)^2, gTMm2 = pmax(-gapTM, 0)^2,
      gRRp2 = pmax(gapRR, 0)^2, gRRm2 = pmax(-gapRR, 0)^2,
      sTM2 = cycTM^2, sTMp2 = pmax(cycTM, 0)^2, sTMm2 = pmax(-cycTM, 0)^2,
      sRR2 = cycRR^2, sRRp2 = pmax(cycRR, 0)^2, sRRm2 = pmax(-cycRR, 0)^2,
      trend_TM_climate = (year - BASE_YEAR) * mean_TM_all,
      trend_RR_climate = (year - BASE_YEAR) * mean_RR_all
    )

  # Identical sample for every model: the outcome is withheld outside the
  # intersection of the samples of the largest models.
  largest <- model_terms(SHOCK_FORMS$arms_climate)
  all_terms <- paste(c(largest$gap, largest$shock), collapse = " + ")
  common <- Reduce(intersect, lapply(TRENDS, function(t) {
    fixest::obs(fit(paste(c(all_terms, t$controls[t$controls != "trend_income"]),
                          collapse = " + "), d, t$fe))
  }))
  d <- withhold(d, common)

  # Income relative to the other units at entry: log income per capita in the
  # unit's first sample year net of that year's sample mean.
  entry <- d[common, ] %>%
    group_by(year) %>%
    mutate(rel_lgdp = lgrp_pc_usd - mean(lgrp_pc_usd)) %>%
    group_by(GID_1) %>%
    filter(year == min(year)) %>%
    ungroup() %>%
    select(GID_1, entry_year = year, init_rel_lgdp = rel_lgdp)
  d <- d %>%
    left_join(entry, by = "GID_1") %>%
    mutate(trend_income = (year - BASE_YEAR) * init_rel_lgdp)

  # (a) Unit growth trends and convergence --------------------------------------------
  used <- d[common, ]
  units <- unit_slopes(used, "dlgrp_pc_usd") %>%
    rename(growth_slope = slope) %>%
    left_join(unit_slopes(used, "gapTM") %>% rename(gapTM_slope = slope), by = "GID_1") %>%
    left_join(unit_slopes(used, "TM") %>% rename(TM_slope = slope), by = "GID_1") %>%
    left_join(used %>%
                group_by(GID_1) %>%
                summarise(mean_growth = mean(dlgrp_pc_usd), years = n(),
                          mean_TM_all = first(mean_TM_all), .groups = "drop"),
              by = "GID_1") %>%
    left_join(entry, by = "GID_1") %>%
    mutate(weight = w)
  trends_out[[w]] <- units

  cross <- list(
    `Average growth ~ initial income` = "mean_growth ~ init_rel_lgdp",
    `Growth trend ~ initial income` = "growth_slope ~ init_rel_lgdp",
    `Growth trend ~ long-run temperature` = "growth_slope ~ mean_TM_all",
    `Growth trend ~ temperature-gap trend` = "growth_slope ~ gapTM_slope",
    `Growth trend ~ all three` = "growth_slope ~ init_rel_lgdp + mean_TM_all + gapTM_slope",
    `Temperature-gap trend ~ long-run temperature` = "gapTM_slope ~ mean_TM_all",
    `Temperature-gap trend ~ initial income` = "gapTM_slope ~ init_rel_lgdp"
  )
  for (name in names(cross)) {
    m <- fixest::feols(as.formula(cross[[name]]), units, vcov = "hetero")
    ct <- fixest::coeftable(m)
    ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
    convergence[[length(convergence) + 1]] <- tibble::tibble(
      weight = w, regression = name, term = rownames(ct), estimate = ct[, 1],
      std_error = ct[, 2], p_value = ct[, 4], r2 = fixest::r2(m, "r2")[[1]],
      units = nobs(m)
    )
  }
  # Within-country convergence: growth on lagged log income.
  within_conv <- fit("lag_lgdp", d, "year + GID_1")
  convergence[[length(convergence) + 1]] <- tibble::tibble(
    weight = w, regression = "Growth ~ lagged log income | unit + year",
    term = "lag_lgdp", estimate = stats::coef(within_conv)[["lag_lgdp"]],
    std_error = fixest::se(within_conv)[["lag_lgdp"]],
    p_value = fixest::pvalue(within_conv)[["lag_lgdp"]],
    r2 = fixest::fitstat(within_conv, "wr2")[[1]], units = n_distinct(used$GID_1)
  )

  # Share of the within variation (net of year and unit effects) absorbed by
  # the unit trends.
  used <- d[common, ]
  for (x in c("dlgrp_pc_usd", "gapTM", "gapRR", GAP_ARMS$var, "cycTM", "cycRR",
              "sTM2", "sRR2", "lag_lgdp")) {
    absorption[[length(absorption) + 1]] <- tibble::tibble(
      weight = w, variable = x,
      share_absorbed_by_trends = 1 - within_var(x, used, TRENDS$trends$fe) /
        within_var(x, used, TRENDS$notrends$fe)
    )
  }

  # Models: trend treatments x shock forms -----------------------------------------------
  for (tr in names(TRENDS)) for (sf in names(SHOCK_FORMS)) {
    message("  ", w, ": ", tr, ", ", sf)
    shocks <- SHOCK_FORMS[[sf]]
    terms <- model_terms(shocks)
    controls <- TRENDS[[tr]]$controls
    rhs <- paste(c(terms$gap, terms$shock, controls), collapse = " + ")
    m <- fit(rhs, d, TRENDS[[tr]]$fe)
    gap_only <- fit(paste(c(terms$gap, controls), collapse = " + "), d, TRENDS[[tr]]$fe)
    stopifnot(nobs(m) == length(common), nobs(gap_only) == length(common))
    shock_TM <- terms$shock[grepl("^sTM", terms$shock)]
    shock_TM_climate <- shock_TM[grepl(":", shock_TM)]
    tests[[length(tests) + 1]] <- tibble::tibble(
      weight = w, trends = tr, shock_form = sf, observations = nobs(m),
      wr2 = fixest::fitstat(m, "wr2")[[1]],
      aic = stats::AIC(m), bic = stats::BIC(m),
      aic_minus_gap_only = stats::AIC(m) - stats::AIC(gap_only),
      gap_arms_p = block_p(m, terms$gap),
      gap_TM_p = block_p(m, terms$gap[grepl("^gTM", terms$gap)]),
      shocks_p = block_p(m, terms$shock),
      shock_TM_p = block_p(m, shock_TM),
      shock_RR_p = block_p(m, terms$shock[grepl("^sRR", terms$shock)]),
      shock_TM_climate_p = if (length(shock_TM_climate)) block_p(m, shock_TM_climate) else NA,
      shock_TM_symmetry_p = if (sf == "arms_climate") {
        symmetry_p(m, "sTMp2", "sTMm2", CLIMATE[["TM"]])
      } else NA,
      shock_RR_symmetry_p = if (sf == "arms_climate") {
        symmetry_p(m, "sRRp2", "sRRm2", CLIMATE[["RR"]])
      } else NA,
      control = if (length(controls)) paste(sprintf("%s = %.3g (p %.3f)", controls,
                                                    stats::coef(m)[controls],
                                                    fixest::pvalue(m)[controls]),
                                            collapse = "; ") else NA
    )
    # Effects of a 1-SD gap or cycle at each long-run climate.
    rows <- bind_rows(
      GAP_ARMS %>% mutate(labels = as.list(label)),
      shocks %>% select(var, variable, labels) %>% mutate(labels = as.list(labels))
    )
    for (i in seq_len(nrow(rows))) for (value in AT_VALUES[[rows$variable[[i]]]]) {
      e <- effect_at(m, rows$var[[i]], CLIMATE[[rows$variable[[i]]]], value)
      for (lab in rows$labels[[i]]) {
        effects[[length(effects) + 1]] <- tibble::tibble(
          weight = w, trends = tr, shock_form = sf, term = lab,
          variable = rows$variable[[i]], climate = value, effect = e$estimate,
          std_error = e$std_error, conf_low = e$conf_low, conf_high = e$conf_high,
          p_value = e$p_value
        )
      }
    }
  }
}

tests <- bind_rows(tests)
effects <- bind_rows(effects)
units_all <- bind_rows(trends_out)
convergence <- bind_rows(convergence)
absorption <- bind_rows(absorption)
utils::write.csv(tests, file.path(OUTPUT_DIR, "model_tests.csv"), row.names = FALSE)
utils::write.csv(effects, file.path(OUTPUT_DIR, "model_effects.csv"), row.names = FALSE)
utils::write.csv(units_all, file.path(OUTPUT_DIR, "unit_trends.csv"), row.names = FALSE)
utils::write.csv(convergence, file.path(OUTPUT_DIR, "convergence.csv"), row.names = FALSE)
utils::write.csv(absorption, file.path(OUTPUT_DIR, "absorption.csv"), row.names = FALSE)

options(width = 220)
cat("\n== Convergence and unit growth trends (cross-section of units, robust SE) ==\n")
print(as.data.frame(convergence %>% mutate(across(where(is.double), ~ signif(.x, 3)))),
      row.names = FALSE)
cat("\n== Share of within variation absorbed by unit trends ==\n")
print(as.data.frame(absorption %>%
  tidyr::pivot_wider(names_from = weight, values_from = share_absorbed_by_trends) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\n== Blocks and fit ==\n")
print(as.data.frame(tests %>%
  select(-control, -aic, -bic) %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
cat("\n== Controls ==\n")
print(as.data.frame(tests %>% filter(!is.na(control)) %>%
  select(weight, trends, shock_form, control)), row.names = FALSE)
cat("\n== Effects of a 1-SD gap or shock (pp), p-values ==\n")
print(as.data.frame(effects %>%
  mutate(value = sprintf("%+.2f (%.2f)", 100 * effect, p_value)) %>%
  select(weight, shock_form, term, climate, trends, value) %>%
  tidyr::pivot_wider(names_from = trends, values_from = value)), row.names = FALSE)

# Plots ---------------------------------------------------------------------------------

trend_levels <- vapply(TRENDS, `[[`, character(1), "label")
colours <- grid_colours(unname(trend_levels))
for (w in names(CLIMATE_WEIGHTS)) {
  # Unit growth trends against income at entry and against the trend of the
  # temperature gap.
  u <- units_all %>%
    filter(weight == w) %>%
    tidyr::pivot_longer(c(init_rel_lgdp, gapTM_slope), names_to = "x_var",
                        values_to = "x") %>%
    mutate(x_var = factor(x_var, levels = c("init_rel_lgdp", "gapTM_slope"),
                          labels = c("Log income per capita at entry, relative to that year's mean",
                                     "Trend of the temperature gap (SD per decade)")),
           growth_slope = 100 * growth_slope)
  p <- ggplot(u, aes(x = x, y = growth_slope)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_point(aes(colour = mean_TM_all), size = 2, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, colour = GRID_INK, linewidth = 0.7,
                fill = "grey80") +
    facet_wrap(~x_var, scales = "free_x") +
    scale_colour_gradient(low = "#3b6fb6", high = "#d0563b",
                          name = paste0("Long-run temperature (", DEG_C, ")")) +
    labs(
      x = NULL, y = "Trend of growth, net of year effects (pp per decade)",
      title = paste0("Unit growth trends absorbed by the unit trends: ", ECON, " × ",
                     CLIMATE_SOURCE, " (", w, " weights)"),
      subtitle = paste("Each point is a country: the linear trend of its growth rate over",
                       "its sample years. Line: OLS fit with 95% band."),
      caption = paste("Countries with at least", MIN_UNIT_YEARS, "sample years.",
                      "Temperature gap: Hamilton trend minus the trailing 30-year mean,",
                      "scaled by the full-record SD.")
    ) +
    theme_results()
  ggsave(file.path(OUTPUT_DIR, paste0("convergence_", w, ".png")), p,
         width = 12, height = 6, dpi = 300)

  # Effects of a 1-SD gap or shock across trend treatments and shock forms.
  e <- effects %>%
    filter(weight == w) %>%
    mutate(
      trends = factor(trend_levels[trends], levels = unname(trend_levels)),
      shock_form = factor(SHOCK_LABELS[shock_form], levels = unname(SHOCK_LABELS)),
      term = factor(term, levels = TERM_LEVELS),
      at = factor(ifelse(variable == "TM", paste(climate, DEG_C), paste(climate, "m")),
                  levels = c(paste(AT_VALUES$TM, DEG_C), paste(AT_VALUES$RR, "m"))),
      across(c(effect, conf_low, conf_high), ~ 100 * .x)
    )
  p <- ggplot(e, aes(x = at, y = effect, colour = trends)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high),
                    position = position_dodge(width = 0.7), size = 0.25, linewidth = 0.5) +
    facet_grid(term ~ shock_form, scales = "free") +
    scale_colour_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(
      x = "Long-run climate", y = "Effect of a 1-SD gap or shock on growth (pp)",
      title = paste0("Gap arms and Hamilton shocks across trend treatments: ", ECON, " × ",
                     CLIMATE_SOURCE, " (", w, " weights)"),
      subtitle = paste("Columns: form of the Hamilton shock. Gap arms quadratic, each",
                       "interacted with the long-run climate. Bars: 95% CI."),
      caption = paste("Identical sample in every model; errors clustered by country.",
                      "Initial income: log income per capita at entry relative to that",
                      "year's sample mean.")
    ) +
    theme_results()
  ggsave(file.path(OUTPUT_DIR, paste0("effects_", w, ".png")), p,
         width = 12, height = 16, dpi = 300)
}
message("Results written to ", normalizePath(OUTPUT_DIR))
