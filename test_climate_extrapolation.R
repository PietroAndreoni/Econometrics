# Does the response of one climate group inform that of another? In both the
# level (BHM) model and the deviations-only model the effect at long-run
# climate m is b + g * m, a line in m fitted across all units, so the
# response in hot countries partly borrows from cold ones and the reverse.
# For each tercile of long-run temperature this script compares:
#
#   pooled         the line fitted on all units, at the tercile's mean climate;
#   in sample      the tercile's own coefficients (i(tercile, x)), without the
#                  line;
#   out of sample  the line fitted without the tercile (the outcome of its
#                  units withheld) and evaluated at its mean climate:
#                  extrapolation for the cold and hot terciles, interpolation
#                  for the temperate one.
#
#   Levels:      TM + TM:mean_TM_all + RR + RR:mean_RR_all; a year d degrees
#                warmer than normal changes growth by (b + g * m) * d.
#   Deviations:  the warm, cold, wet and dry squared anomalies against the
#                trailing WINDOW-year mean (scaled by the SD over the entire
#                climate series), each interacted with its long-run climate;
#                a year d SD warmer changes growth by (b_warm + g_warm * m) *
#                d^2, a year d SD colder by (b_cold + g_cold * m) * d^2.
#
# The figures show the response to a year up to 2 degrees warmer or colder
# than normal, with 1 SD set to 1 degree Celsius (as in
# plot_bhm_vs_anomaly_curves.R), so both models share one axis. Precipitation
# terms enter every model as in the pooled specification. Both models use the
# same country-years; the terciles split the sample's units by mean_TM_all.
#
# Grid: PWT x CRU TS, population weights, 30-year trailing window (set
# ECON_DATASETS, CLIMATE_SOURCES and WINDOW for others). Fixed effects as in
# main_analysis.Rmd, errors clustered by country.
#
# Output: results/climate_extrapolation/ (tercile_estimates_lag<WINDOW>.csv,
# anomaly_curves_lag<WINDOW>.csv,
# anomaly_curves_<dataset>_<climate source>_lag<WINDOW>.png)

source("load_functions.R")

ECON_DATASETS <- c(PWT = "PWT")
CLIMATE_SOURCES <- c(`CRU TS` = "CRU TS")
CLIMATE_WEIGHT <- "concurrent population"
WINDOW <- 30L
BASELINE <- paste0("lag_", WINDOW)
SUFFIX <- paste0("_lag", WINDOW)
SD_BASELINE <- "all"
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"
CONF_LEVEL <- 0.95
ANOMALIES <- seq(-2, 2, by = 0.05)
TERCILE_LABELS <- c("Cold", "Temperate", "Hot")
DEG_C <- "°C"

# Coefficient names: the level slope and its interaction, and the warm and
# cold squared anomalies and their interactions.
MODELS <- list(
  levels = list(
    label = "Level model (BHM)",
    pooled = "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
    free = "i(tercile, TM) + RR + RR:mean_RR_all",
    up = "TM", down = "TM", power = 1
  ),
  deviations = list(
    label = "Deviations-only model",
    pooled = paste("zTMp2 + zTMm2 + zTMp2:mean_TM_all + zTMm2:mean_TM_all",
                   "+ zRRp2 + zRRm2 + zRRp2:mean_RR_all + zRRm2:mean_RR_all"),
    free = paste("i(tercile, zTMp2) + i(tercile, zTMm2)",
                 "+ zRRp2 + zRRm2 + zRRp2:mean_RR_all + zRRm2:mean_RR_all"),
    up = "zTMp2", down = "zTMm2", power = 2
  )
)
ESTIMATE_LABELS <- c(pooled = "Pooled (all units)",
                     in_sample = "In sample (the tercile's own coefficients)",
                     out_of_sample = "Out of sample (fitted without the tercile)")

OUTPUT_DIR <- file.path("results", "climate_extrapolation")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit <- function(terms, data) {
  fit_panel_model(terms, data, fixed_effects = FIXED_EFFECTS, cluster = ~GID_0)
}
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
withhold <- function(data, rows) {
  data$dlgrp_pc_usd[rows] <- NA
  data
}

# Effect of a year d SD (= d degrees) above or below normal, for each d in
# `anomalies`: the line b + g * m of the up (d >= 0) or down (d < 0) term at
# climate m, or, for the in-sample model, the tercile's own coefficient.
anomaly_curve <- function(model, spec, anomalies, m = NULL, tercile = NULL) {
  terms <- names(stats::coef(model))
  G <- t(vapply(anomalies, function(d) {
    x <- if (d >= 0) spec$up else spec$down
    w <- stats::setNames(numeric(length(terms)), terms)
    scale <- sign(d)^(spec$power == 1) * abs(d)^spec$power
    if (is.null(tercile)) {
      w[[x]] <- scale
      w[[interaction_name(x, "mean_TM_all", terms)]] <- scale * m
    } else {
      w[[grep(paste0("tercile::", tercile, ":", x, "$"), terms, value = TRUE)]] <- scale
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(anomaly = anomalies)
}

curves <- list()
for (econ in names(ECON_DATASETS)) for (source_name in names(CLIMATE_SOURCES)) {
  message("Fitting ", econ, " x ", source_name)
  cfg <- panel_config(
    econ_year_min = min(load_econ_panel(resolve_econ_source(ECON_DATASETS[[econ]]))$year),
    climate_source = CLIMATE_SOURCES[[source_name]],
    climate_weight = CLIMATE_WEIGHT,
    clim_history_years = 1L
  )
  d <- suppressMessages(build_dat(ECON_DATASETS[[econ]], config = cfg,
                                  baseline = BASELINE, sd_baseline = SD_BASELINE)) %>%
    mutate(zTMp2 = abs_zTMp^2, zTMm2 = abs_zTMm^2,
           zRRp2 = abs_zRRp^2, zRRm2 = abs_zRRm^2)

  # Common sample of both models, and terciles of its units by long-run
  # temperature.
  common <- Reduce(intersect, lapply(MODELS, function(s) fixest::obs(fit(s$pooled, d))))
  d <- withhold(d, setdiff(seq_len(nrow(d)), common))
  units <- d[common, ] %>%
    distinct(GID_1, mean_TM_all) %>%
    mutate(tercile = ntile(mean_TM_all, 3))
  d <- d %>% left_join(units %>% select(GID_1, tercile), by = "GID_1")
  ranges <- units %>%
    group_by(tercile) %>%
    summarise(mean = mean(mean_TM_all), low = min(mean_TM_all),
              high = max(mean_TM_all), units = n(), .groups = "drop")

  for (model in names(MODELS)) {
    spec <- MODELS[[model]]
    pooled <- fit(spec$pooled, d)
    free <- fit(spec$free, d)
    for (k in 1:3) {
      m_k <- ranges$mean[ranges$tercile == k]
      without_k <- fit(spec$pooled, withhold(d, which(is.na(d$tercile) | d$tercile == k)))
      curves[[length(curves) + 1]] <- bind_rows(
        anomaly_curve(pooled, spec, ANOMALIES, m = m_k) %>% mutate(estimate_type = "pooled"),
        anomaly_curve(free, spec, ANOMALIES, tercile = k) %>% mutate(estimate_type = "in_sample"),
        anomaly_curve(without_k, spec, ANOMALIES, m = m_k) %>%
          mutate(estimate_type = "out_of_sample")
      ) %>%
        mutate(econ = econ, climate_source = source_name, model = model,
               tercile = TERCILE_LABELS[[k]], climate = m_k,
               climate_low = ranges$low[ranges$tercile == k],
               climate_high = ranges$high[ranges$tercile == k],
               units = ranges$units[ranges$tercile == k])
    }
  }
}

curves <- bind_rows(curves)
utils::write.csv(curves, file.path(OUTPUT_DIR, paste0("anomaly_curves", SUFFIX, ".csv")), row.names = FALSE)

# The effect of a year 1 degree warmer and 1 degree colder, by tercile.
estimates <- curves %>%
  filter(anomaly %in% c(-1, 1)) %>%
  mutate(direction = ifelse(anomaly > 0, "warm", "cold")) %>%
  select(econ, climate_source, model, tercile, climate, units, direction, estimate_type,
         estimate, std_error, conf_low, conf_high)
utils::write.csv(estimates, file.path(OUTPUT_DIR, paste0("tercile_estimates", SUFFIX, ".csv")), row.names = FALSE)

options(width = 220)
cat("\n== Effect of a year 1 degree (1 SD) warmer: pooled, in sample, out of sample ==\n")
print(as.data.frame(estimates %>%
  filter(direction == "warm") %>%
  mutate(value = sprintf("%+.4f (%.4f)", estimate, std_error)) %>%
  select(econ, climate_source, model, tercile, estimate_type, value) %>%
  tidyr::pivot_wider(names_from = estimate_type, values_from = value)), row.names = FALSE)

# Plots ---------------------------------------------------------------------------------

estimate_levels <- unname(ESTIMATE_LABELS)
colours <- grid_colours(estimate_levels)
for (econ in names(ECON_DATASETS)) for (source_name in names(CLIMATE_SOURCES)) {
  d <- curves %>%
    filter(econ == !!econ, climate_source == source_name) %>%
    mutate(
      estimate_type = factor(ESTIMATE_LABELS[estimate_type], levels = estimate_levels),
      model = factor(vapply(MODELS, `[[`, character(1), "label")[model],
                     levels = vapply(MODELS, `[[`, character(1), "label")),
      tercile = factor(
        sprintf("%s tercile: mean %.1f %s\n(%.1f to %.1f, %d units)", tercile, climate,
                DEG_C, climate_low, climate_high, units),
        levels = unique(sprintf("%s tercile: mean %.1f %s\n(%.1f to %.1f, %d units)",
                                tercile, climate, DEG_C, climate_low, climate_high,
                                units)[order(match(tercile, TERCILE_LABELS))])
      ),
      across(c(estimate, conf_low, conf_high), ~ 100 * .x)
    )
  p <- ggplot(d, aes(x = anomaly, y = estimate, colour = estimate_type,
                     fill = estimate_type)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(model ~ tercile) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    labs(
      x = paste0("Temperature anomaly in the year (", DEG_C, ", with SD = 1 ", DEG_C, ")"),
      y = "Effect on growth that year (pp)",
      title = paste0("A year warmer or colder than normal, by climate: ", econ, " × ",
                     source_name),
      subtitle = paste("Pooled: b + g · m fitted on all units at the tercile's mean",
                       "climate. In sample: the tercile's own coefficients. Out of sample:",
                       "b + g · m fitted\nwithout the tercile (extrapolation for the cold",
                       "and hot terciles, interpolation for the temperate one).",
                       "Bands: 95% CI."),
      caption = paste0("Population weights; deviations against the trailing ", WINDOW,
                       "-year mean, scaled by the full-record SD; fixed effects year +",
                       " country + country trends; errors clustered by country.")
    ) +
    theme_results()
  file <- paste0("anomaly_curves_", econ, "_", gsub("[^A-Za-z0-9]+", "", source_name), SUFFIX, ".png")
  ggsave(file.path(OUTPUT_DIR, file), p, width = 12, height = 8, dpi = 300)
}
message("Results written to ", normalizePath(OUTPUT_DIR))
