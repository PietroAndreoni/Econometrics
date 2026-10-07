# The slow climate component next to the gap and the shock, PWT x UDelaware.
#
# Temperature splits exactly into the trailing 30-year mean, the gap and the
# Hamilton shock:
#   TM = mean_TM_lag_30 + (trend_TM - mean_TM_lag_30) + (TM - trend_TM).
# The model enters the trailing mean in levels (degrees Celsius, metres for
# precipitation), so a non-zero coefficient is a growth effect of the adapted
# climate, next to the signed gap and the squared shock (both scaled by the
# SD over the entire climate series), each level and gap term interacted with
# the long-run climate:
#   gTM + gTM:mean_TM_all + mean_TM_lag_30 + mean_TM_lag_30:mean_TM_all + sTM2
#   (and the same for precipitation).
#
# Expectation: Bilal and Kanzig's Hamilton trend (h = 2, years t-2 and t-3),
# fitted in real time on the 30 years up to t-2 only, from each country's
# entire climate series. Fixed effects year + country (no unit trends), area
# and population weights, errors clustered by country.
#
# Output: results/gap_climate_level/ (coefficients.csv, effects.csv)

source("load_functions.R")

ECON <- "PWT"
CLIMATE_SOURCE <- "UDelaware"
CLIMATE_WEIGHTS <- c(pop = "concurrent population", area = "area")
HAMILTON_H <- 2L
FILTERS <- list(bk = list(lags = 2L, window = 30L))
FIXED_EFFECTS <- c(notrends = "year + GID_1")
TERMS <- paste(
  "gTM + gTM:mean_TM_all + gRR + gRR:mean_RR_all",
  "+ mean_TM_lag_30 + mean_TM_lag_30:mean_TM_all",
  "+ mean_RR_lag_30 + mean_RR_lag_30:mean_RR_all + sTM2 + sRR2"
)
# Effects reported: a 1-SD gap, a 1-degree (0.1-m) warmer (wetter) trailing
# mean, at long-run climates of 5, 15 and 25 degrees (0.5, 1.5, 2.5 m).
EFFECTS <- tibble::tribble(
  ~term,            ~climate,      ~size, ~label,
  "gTM",            "mean_TM_all", 1,     "Temperature gap, +1 SD",
  "mean_TM_lag_30", "mean_TM_all", 1,     "Trailing temperature mean, +1 °C",
  "gRR",            "mean_RR_all", 1,     "Precipitation gap, +1 SD",
  "mean_RR_lag_30", "mean_RR_all", 0.1,   "Trailing precipitation mean, +100 mm"
)
AT_VALUES <- list(mean_TM_all = c(5, 15, 25), mean_RR_all = c(0.5, 1.5, 2.5))

OUTPUT_DIR <- file.path("results", "gap_climate_level")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}

coefficients <- list()
effects <- list()
for (w in names(CLIMATE_WEIGHTS)) {
  cfg <- panel_config(
    econ_year_min = min(load_econ_panel(resolve_econ_source(ECON))$year),
    climate_source = CLIMATE_SOURCE,
    climate_weight = CLIMATE_WEIGHTS[[w]],
    clim_history_years = 1L
  )
  panel <- suppressMessages(build_dat(ECON, config = cfg, baseline = "lag_30",
                                      sd_baseline = "all"))
  # Hamilton trends of each filter from the entire climate series.
  climate <- suppressMessages(load_climate_panel(unique(panel$gadm_level), cfg)) %>%
    filter(GID_1 %in% unique(panel$GID_1)) %>%
    group_by(GID_1) %>%
    arrange(year, .by_group = TRUE)
  for (f in names(FILTERS)) {
    climate <- climate %>%
      mutate(!!paste0("trend_TM_", f) := hamilton_trend(TM, HAMILTON_H, FILTERS[[f]]$lags,
                                                        window = FILTERS[[f]]$window),
             !!paste0("trend_RR_", f) := hamilton_trend(RR, HAMILTON_H, FILTERS[[f]]$lags,
                                                        window = FILTERS[[f]]$window))
  }
  panel <- panel %>%
    left_join(climate %>% ungroup() %>% select(GID_1, year, starts_with("trend_TM_"),
                                               starts_with("trend_RR_")),
              by = c("GID_1", "year"))

  data_for <- function(f) {
    d <- panel
    for (v in c("TM", "RR")) {
      trend <- d[[paste0("trend_", v, "_", f)]]
      sd_all <- d[[paste0("sd_", v, "_all")]]
      d[[paste0("g", v)]] <- (trend - d[[paste0("mean_", v, "_lag_30")]]) / sd_all
      d[[paste0("s", v, "2")]] <- ((d[[v]] - trend) / sd_all)^2
    }
    d
  }
  # Identical sample for every filter and fixed-effects set (when several are run).
  common <- Reduce(intersect, lapply(names(FILTERS), function(f) {
    fixest::obs(fit_panel_model(TERMS, data_for(f), fixed_effects = FIXED_EFFECTS[["notrends"]],
                                cluster = ~GID_0))
  }))

  models <- list()
  for (f in names(FILTERS)) for (fe in names(FIXED_EFFECTS)) {
    d <- data_for(f)
    d$dlgrp_pc_usd[setdiff(seq_len(nrow(d)), common)] <- NA
    m <- fit_panel_model(TERMS, d, fixed_effects = FIXED_EFFECTS[[fe]], cluster = ~GID_0)
    key <- paste(f, fe)
    models[[key]] <- m
    ct <- fixest::coeftable(m)
    coefficients[[length(coefficients) + 1]] <- tibble::tibble(
      weight = w, filter = f, trends = fe, term = rownames(ct), estimate = ct[, 1],
      std_error = ct[, 2], p_value = ct[, 4], observations = nobs(m)
    )
    terms <- names(stats::coef(m))
    for (i in seq_len(nrow(EFFECTS))) for (value in AT_VALUES[[EFFECTS$climate[[i]]]]) {
      weights <- stats::setNames(EFFECTS$size[[i]] * c(1, value),
                                 c(EFFECTS$term[[i]],
                                   interaction_name(EFFECTS$term[[i]], EFFECTS$climate[[i]], terms)))
      e <- lincom(m, weights)
      effects[[length(effects) + 1]] <- tibble::tibble(
        weight = w, filter = f, trends = fe, effect = EFFECTS$label[[i]], climate = value,
        estimate = e$estimate, std_error = e$std_error, p_value = e$p_value
      )
    }
  }
  cat("\n==== ", ECON, " x ", CLIMATE_SOURCE, ", ", w, " weights ====\n", sep = "")
  print(etable(models, headers = names(models), digits = 4, fitstat = ~ n + wr2 + aic + bic))
}

coefficients <- bind_rows(coefficients)
effects <- bind_rows(effects)
utils::write.csv(coefficients, file.path(OUTPUT_DIR, "coefficients.csv"), row.names = FALSE)
utils::write.csv(effects, file.path(OUTPUT_DIR, "effects.csv"), row.names = FALSE)

options(width = 220)
cat("\n== Effects on growth (pp), with p-values ==\n")
print(as.data.frame(effects %>%
  mutate(value = sprintf("%+.2f (%.2f)", 100 * estimate, p_value),
         model = paste(filter, trends)) %>%
  select(weight, effect, climate, model, value) %>%
  tidyr::pivot_wider(names_from = model, values_from = value)), row.names = FALSE)
message("Results written to ", normalizePath(OUTPUT_DIR))
