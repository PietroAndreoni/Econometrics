# Level response against two anomaly responses, for temperature and for
# precipitation, on one common sample: UDelaware, population weights, all
# available years of the economic dataset (run option, WB by default:
# Rscript plot_bhm_vs_anomaly_curves.R PWT).
#
#   Levels:      TM + TM:mean_TM_all + RR + RR:mean_RR_all
#   Deviations:  abs_zTMp^2 + abs_zTMm^2 + abs_zRRp^2 + abs_zRRm^2, the warm,
#                cold, wet and dry anomalies against the trailing 30-year mean
#                (scaled by the SD over the entire climate series), each
#                interacted with its long-run climate (mean_TM_all, mean_RR_all)
#   Adaptation + shock: the same 30-year anomaly split into a slow adaptation
#                gap and a short-term surprise (see test_adaptation_shock_model.R),
#                  gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
#                  cycTM = (TM - trend_TM) / sd_TM_all,
#                with quadratic arms on the gap (gapTMp^2, gapTMm^2, and wet and
#                dry for precipitation, each interacted with the long-run
#                climate) plus the squared Hamilton surprises cycTM^2, cycRR^2.
#
# Both anomaly models use the trailing 30-year window; the adaptation + shock
# model adds the Hamilton trend (lags 2 to 5). For the figures 1 SD is set to
# 1 degree Celsius for temperature and 100 mm for precipitation; the other
# variable is held at normal.
#
#   B. A year d SD above or below normal, in units with three long-run
#      climates: (b + g * m) * d * (units per SD) in the level model;
#      (b_up + g_up * m) * d^2 or (b_down + g_down * m) * d^2 in the deviations
#      model; c * d^2 in the adaptation + shock model, where a one-off anomaly
#      is a pure Hamilton surprise (the trend uses lags 2 to 5 only).
#   C. (and D for precipitation) A permanent 1-SD step at year 0: the
#      cumulative effect on log GDP. Deviations model: the anomaly against the
#      trailing 30-year mean is 1 - t / 30 in year t. Adaptation + shock: the
#      Hamilton trend rises by the lag coefficients as the step enters lags 2
#      to 5 (median coefficients across the sample's units), the surprise is 1
#      minus that rise, and the gap is the rise minus the share of the step in
#      the 30-year mean; each year adds the shock and the gap arm. The level
#      model, b + g * m every year, is left out of these panels (its scale
#      would hide the others) and its 30-year values are given in the subtitle.
#
# All three models are fitted on the same country-years. Bands are 95%
# intervals from draws of the coefficients from their clustered (by country)
# sampling distribution.
#
# Output: results/bhm_vs_anomaly_curves/bhm_vs_anomaly_curves_{temperature,
# precipitation}[_<dataset>].png and .csv

source("load_functions.R")

ECON <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(ECON)) ECON <- "WB"
CLIMATE_SOURCE <- "UDelaware"
CLIMATE_WEIGHT <- "concurrent population"
BASELINE <- "lag_30"
BASELINE_YEARS <- 30L
SD_BASELINE <- "all"
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"
LEVEL_TERMS <- "TM + TM:mean_TM_all + RR + RR:mean_RR_all"
ANOMALY_TERMS <- paste(
  "abs_zTMp^2 + abs_zTMm^2 + abs_zRRp^2 + abs_zRRm^2",
  "+ abs_zTMp^2:mean_TM_all + abs_zTMm^2:mean_TM_all",
  "+ abs_zRRp^2:mean_RR_all + abs_zRRm^2:mean_RR_all"
)
ADAPT_TERMS <- paste(
  "gapTMp^2 + gapTMm^2 + gapRRp^2 + gapRRm^2",
  "+ gapTMp^2:mean_TM_all + gapTMm^2:mean_TM_all",
  "+ gapRRp^2:mean_RR_all + gapRRm^2:mean_RR_all",
  "+ cycTM^2 + cycRR^2"
)
HAMILTON_LAGS <- 2:5
ANOMALIES <- seq(-2, 2, by = 0.05)
HORIZON <- 30L
DRAWS <- 2000L
SEED <- 20261001L
MODEL_LABELS <- c(levels = "Level model (BHM)",
                  anomalies = "Deviations-only model",
                  adapt = "Adaptation gap + Hamilton shock")
DEG_C <- "°C"

# Per variable: the long-run climates, the size of 1 SD in display units and
# in the level model's units, and the permanent steps to show (in SD).
VARIABLES <- list(
  TM = list(name = "temperature", climate_var = "mean_TM_all",
            climates = c(5, 15, 25), climate_unit = DEG_C,
            unit = DEG_C, sd_display = 1, sd_level = 1,
            up = "warm", down = "cold", up_cmp = "warmer", down_cmp = "colder",
            steps = c(warmer = 1)),
  RR = list(name = "precipitation", climate_var = "mean_RR_all",
            climates = c(0.5, 1.5, 2.5), climate_unit = "m",
            unit = "mm", sd_display = 100, sd_level = 0.1,
            up = "wet", down = "dry", up_cmp = "wetter", down_cmp = "drier",
            steps = c(wetter = 1, drier = -1))
)

OUTPUT_DIR <- file.path("results", "bhm_vs_anomaly_curves")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
file_stem <- function(variable) {
  stem <- paste0("bhm_vs_anomaly_curves_", VARIABLES[[variable]]$name)
  if (ECON == "WB") stem else paste0(stem, "_", ECON)
}

# Models ---------------------------------------------------------------------------

cfg <- panel_config(
  econ_year_min = min(load_econ_panel(resolve_econ_source(ECON))$year),
  climate_source = CLIMATE_SOURCE,
  climate_weight = CLIMATE_WEIGHT,
  clim_history_years = 1L
)
panel <- suppressMessages(build_dat(ECON, config = cfg, baseline = BASELINE,
                                    sd_baseline = SD_BASELINE)) %>%
  mutate(
    gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
    cycTM = (TM - trend_TM) / sd_TM_all,
    gapRR = (trend_RR - mean_RR_lag_30) / sd_RR_all,
    cycRR = (RR - trend_RR) / sd_RR_all,
    gapTMp = pmax(gapTM, 0), gapTMm = pmax(-gapTM, 0),
    gapRRp = pmax(gapRR, 0), gapRRm = pmax(-gapRR, 0)
  )

fit <- function(terms, data) {
  fit_panel_model(terms, data, fixed_effects = FIXED_EFFECTS, cluster = ~GID_0)
}
# The common sample: country-years every model can use. The outcome is
# withheld elsewhere, so lags and baselines still resolve.
common <- Reduce(intersect, lapply(c(ANOMALY_TERMS, ADAPT_TERMS, LEVEL_TERMS),
                                   function(terms) fixest::obs(fit(terms, panel))))
sample_data <- panel
sample_data$dlgrp_pc_usd[setdiff(seq_len(nrow(panel)), common)] <- NA
m_levels <- fit(LEVEL_TERMS, sample_data)
m_anomalies <- fit(ANOMALY_TERMS, sample_data)
m_adapt <- fit(ADAPT_TERMS, sample_data)
stopifnot(nobs(m_levels) == length(common), nobs(m_anomalies) == length(common),
          nobs(m_adapt) == length(common))

# Coefficient draws, first row the point estimates. Interaction names are
# looked up in either order.
draws <- function(model) {
  set.seed(SEED)
  b <- stats::coef(model)
  rbind(b, MASS::mvrnorm(DRAWS, b, stats::vcov(model)))
}
term <- function(B, a, b = NULL) {
  if (is.null(b)) return(B[, a])
  hit <- intersect(c(paste0(a, ":", b), paste0(b, ":", a)), colnames(B))
  B[, hit[[1]]]
}
B_levels <- draws(m_levels)
B_anomalies <- draws(m_anomalies)
B_adapt <- draws(m_adapt)

# Hamilton step response: median lag coefficients across the sample's units,
# from the same regression as hamilton_trend() (y_t on y_{t-2}, ..., y_{t-5}
# over each unit's entire climate series), for each variable.
units <- unique(panel$GID_1[common])
climate <- load_climate_panel(unique(panel$gadm_level), cfg) %>%
  filter(GID_1 %in% units) %>%
  arrange(GID_1, year)
hamilton_beta <- lapply(stats::setNames(names(VARIABLES), names(VARIABLES)), function(v) {
  climate %>%
    group_by(GID_1) %>%
    group_modify(function(u, key) {
      y <- u[[v]]
      X <- sapply(HAMILTON_LAGS, function(j) dplyr::lag(y, j))
      ok <- stats::complete.cases(cbind(y, X))
      b <- stats::coef(stats::lm(y[ok] ~ X[ok, ]))[-1]
      tibble::as_tibble(as.list(stats::setNames(b, paste0("lag", HAMILTON_LAGS))))
    }) %>%
    ungroup() %>%
    summarise(across(starts_with("lag"), stats::median)) %>%
    unlist()
})

# Point estimate (first draw) and 95% band over the remaining draws, in
# percentage points of growth or of log GDP.
summarise_draws <- function(values) {
  c(estimate = 100 * values[[1]],
    lower = 100 * stats::quantile(values[-1], 0.025, names = FALSE),
    upper = 100 * stats::quantile(values[-1], 0.975, names = FALSE))
}
curve_frame <- function(x, f, ...) {
  bind_rows(lapply(x, function(v) {
    tibble::as_tibble(as.list(summarise_draws(f(v)))) %>% mutate(x = v, ...)
  }))
}

# Curves ---------------------------------------------------------------------------

curves_for <- function(v) {
  spec <- VARIABLES[[v]]
  cv <- spec$climate_var
  lev_b <- term(B_levels, v)
  lev_g <- term(B_levels, v, cv)
  up_b <- term(B_anomalies, paste0("I(abs_z", v, "p^2)"))
  up_g <- term(B_anomalies, paste0("I(abs_z", v, "p^2)"), cv)
  down_b <- term(B_anomalies, paste0("I(abs_z", v, "m^2)"))
  down_g <- term(B_anomalies, paste0("I(abs_z", v, "m^2)"), cv)
  gap_up_b <- term(B_adapt, paste0("I(gap", v, "p^2)"))
  gap_up_g <- term(B_adapt, paste0("I(gap", v, "p^2)"), cv)
  gap_down_b <- term(B_adapt, paste0("I(gap", v, "m^2)"))
  gap_down_g <- term(B_adapt, paste0("I(gap", v, "m^2)"), cv)
  shock <- term(B_adapt, paste0("I(cyc", v, "^2)"))
  beta <- hamilton_beta[[v]]
  trend_rise <- function(t) sum(beta[HAMILTON_LAGS <= t])

  # Effects of an anomaly of d SD, and of one year after a permanent step of
  # s SD, for a unit with long-run climate m.
  level_effect <- function(d, m) (lev_b + lev_g * m) * spec$sd_level * d
  deviation_effect <- function(d, m) {
    if (d >= 0) (up_b + up_g * m) * d^2 else (down_b + down_g * m) * d^2
  }
  adapt_year <- function(t, m, s) {
    cycle <- s * (1 - trend_rise(t))
    gap <- s * (trend_rise(t) - min(t, BASELINE_YEARS) / BASELINE_YEARS)
    arm <- if (gap >= 0) gap_up_b + gap_up_g * m else gap_down_b + gap_down_g * m
    shock * cycle^2 + arm * gap^2
  }

  one_year <- bind_rows(lapply(spec$climates, function(m) {
    bind_rows(
      curve_frame(ANOMALIES, function(d) level_effect(d, m), model = "levels", climate = m),
      curve_frame(ANOMALIES, function(d) deviation_effect(d, m), model = "anomalies",
                  climate = m),
      curve_frame(ANOMALIES, function(d) shock * d^2, model = "adapt", climate = m)
    )
  }))
  steps <- bind_rows(lapply(names(spec$steps), function(direction) {
    s <- spec$steps[[direction]]
    bind_rows(lapply(spec$climates, function(m) {
      bind_rows(
        curve_frame(0:HORIZON, function(t) level_effect(s, m) * (t + 1),
                    model = "levels", climate = m),
        curve_frame(0:HORIZON, function(t) {
          Reduce(`+`, lapply(0:t, function(u) {
            deviation_effect(s * max(1 - u / BASELINE_YEARS, 0), m)
          }))
        }, model = "anomalies", climate = m),
        curve_frame(0:HORIZON, function(t) {
          Reduce(`+`, lapply(0:t, function(u) adapt_year(u, m, s)))
        }, model = "adapt", climate = m)
      )
    })) %>% mutate(direction = direction)
  }))
  list(one_year = one_year %>% mutate(x = x * spec$sd_display), steps = steps)
}

# Plots ----------------------------------------------------------------------------------

model_levels <- unname(MODEL_LABELS)
colours <- grid_colours(model_levels)
label_models <- function(d, spec) {
  climate_labels <- paste0("Long-run mean ", spec$climates, " ", spec$climate_unit)
  d %>% mutate(
    model = factor(MODEL_LABELS[model], levels = model_levels),
    climate = factor(paste0("Long-run mean ", climate, " ", spec$climate_unit),
                     levels = climate_labels)
  )
}
sd_text <- function(spec) paste0("1 SD = ", spec$sd_display, " ", spec$unit)

plot_variable <- function(v) {
  spec <- VARIABLES[[v]]
  curves <- curves_for(v)
  beta_sum <- sum(hamilton_beta[[v]])

  p_one <- ggplot(label_models(curves$one_year, spec),
                  aes(x = x, y = estimate, colour = model, fill = model)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~climate, nrow = 1) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, name = NULL) +
    labs(x = paste0(tools::toTitleCase(spec$name), " anomaly in the year (",
                    spec$unit, ", ", sd_text(spec), ")"),
         y = "Effect on growth that year (pp)",
         title = paste0("A. A year ", spec$up_cmp, " or ", spec$down_cmp, " than normal"),
         subtitle = paste("Level model: linear; deviations model: climate-dependent",
                          spec$up, "and", spec$down, "curvature; adaptation + shock",
                          "model: a one-off\nanomaly is a pure Hamilton surprise,",
                          "symmetric and squared")) +
    theme_results()

  step_panels <- lapply(seq_along(spec$steps), function(i) {
    direction <- names(spec$steps)[[i]]
    s <- spec$steps[[i]]
    d <- curves$steps %>% filter(direction == !!direction)
    level_30 <- d %>% filter(model == "levels", x == HORIZON) %>% arrange(climate)
    level_note <- paste0(
      "Level model (not shown) after ", HORIZON, " years: ",
      paste(sprintf("%+.1f pp at %g %s", level_30$estimate, level_30$climate,
                    spec$climate_unit), collapse = ", ")
    )
    step_size <- paste0(if (s > 0) "+" else "-", spec$sd_display, " ", spec$unit)
    ggplot(label_models(d %>% filter(model != "levels"), spec),
           aes(x = x, y = estimate, colour = model, fill = model)) +
      geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~climate, nrow = 1) +
      scale_colour_manual(values = colours, name = NULL) +
      scale_fill_manual(values = colours, name = NULL) +
      labs(x = paste0("Years after a permanent ", step_size, " change"),
           y = "Cumulative effect on log GDP (pp)",
           title = paste0(LETTERS[i + 1], ". A permanent ", step_size, " (", direction,
                          ") change"),
           subtitle = paste0(
             "Deviations model: the anomaly fades as the trailing ", BASELINE_YEARS,
             "-year baseline absorbs the change. Adaptation + shock: the Hamilton\n",
             "trend absorbs ", round(100 * beta_sum), "% of the step within 5 years ",
             "(surprise) while the ", BASELINE_YEARS, "-year baseline lags (gap).\n",
             level_note)) +
      theme_results()
  })

  caption <- paste0(
    ECON, " x ", CLIMATE_SOURCE, ", ", CLIMATE_WEIGHT, " weights; ", length(common),
    " country-years common to all models. Both anomaly models use the trailing ",
    BASELINE_YEARS, "-year mean;\nthe adaptation + shock model splits it into the ",
    "gap to the Hamilton trend (quadratic arms) and the squared Hamilton surprise. ",
    sd_text(spec), ";\nthe other climate variable held at normal. Bands: 95% ",
    "intervals from ", DRAWS, " draws of the country-clustered coefficients."
  )
  figure <- patchwork::wrap_plots(c(list(p_one), step_panels), ncol = 1) +
    patchwork::plot_annotation(
      title = paste0(tools::toTitleCase(spec$name), ": level versus anomaly responses"),
      caption = caption,
      theme = theme(plot.caption = element_text(colour = GRID_MUTED, hjust = 0),
                    plot.title = element_text(face = "bold", size = 15))
    )
  ggsave(file.path(OUTPUT_DIR, paste0(file_stem(v), ".png")), figure,
         width = 11, height = 1.5 + 4.6 * (1 + length(step_panels)), dpi = 300)
  write.csv(bind_rows(curves$one_year %>% mutate(panel = "one year"),
                      curves$steps %>% mutate(panel = "permanent step")),
            file.path(OUTPUT_DIR, paste0(file_stem(v), ".csv")), row.names = FALSE)
}

for (v in names(VARIABLES)) plot_variable(v)
message("Hamilton coefficients: ",
        paste(names(hamilton_beta), sapply(hamilton_beta, function(b) {
          paste0(paste(round(b, 2), collapse = ", "), " (sum ", round(sum(b), 2), ")")
        }), sep = " ", collapse = "; "),
        "; figures written to ", normalizePath(OUTPUT_DIR))
