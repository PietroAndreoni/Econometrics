# Decompose each DOSE region's weather anomaly into a global, a national and a
# subnational component, and regress regional growth on each component.
#
# For v in TM, RR, with the anomaly z defined at each level exactly as
# build_dat() defines zTM and zRR (baseline BASELINE, deviation DEVIATION):
#   z_v_glob     world anomaly
#   z_v_nat_res  national anomaly minus the world anomaly
#   z_v_reg_res  regional anomaly minus the national anomaly
# so that zv = z_v_glob + z_v_nat_res + z_v_reg_res holds identically.
# With standardized anomalies each level is scaled by its own SD, so the
# components are differences of standardized anomalies; with DEVIATION =
# "absolute" they are in degrees Celsius and metres of precipitation.
#
# Climate: CRU TS, area-weighted (the Weighted Climate Dataset's "unweighted",
# a plain average over grid cells), at GADM1, GADM0 and world level. The world
# series is not part of the cached climate panels, so it is downloaded and its
# moments derived here, with the same code (derive_climate_moments()).
#
# Year fixed effects absorb the world component entirely, and country-year
# fixed effects the national one, so the models differ in which components
# they can identify:
#   Total anomaly, year FE           zTM + zRR (reference)
#   Components, year FE              national + regional
#   Total anomaly, no year FE        zTM + zRR
#   Components, no year FE           world + national + regional
#   Components, country-year FE      regional only
# All models also carry region-specific quadratic trends and share one sample.
# The world component varies only by year, so region-clustered standard errors
# overstate its precision; every table also reports two-way clustering by
# country and year.
#
# Specification: main_spec() (functions/analysis_spec.R).
# Output: results/dose_climate_levels/.

source("load_functions.R")

SPEC <- main_spec()
BASELINE <- "lag_30"
DEVIATION <- "absolute"
CONFIG_ARGS <- list(climate_source = "CRU TS", climate_weight = "area")
CONFIG <- do.call(panel_config, c(
  list(econ_year_min = SPEC$econ_year_min,
       clim_history_years = SPEC$clim_history_years),
  CONFIG_ARGS
))
REGION_TRENDS <- "GID_1[year] + GID_1[year^2]"   # SPEC$fixed_effects without year
VCOVS <- list("Cluster: GID_1" = ~GID_1, "Two-way: GID_0 + year" = ~GID_0 + year)
COMPONENTS <- c(glob = "World", nat_res = "National - world",
                reg_res = "Regional - national")
OUT <- file.path("results", "dose_climate_levels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# (a) Regional, national and world climate ------------------------------------------

# World series in the layout of load_climate_series(), under the id "WLD".
load_world_series <- function(config) {
  fetch <- function(variable, value_code) {
    tibble::as_tibble(wcd_get(
      variable = variable,
      source = config$climate_source,
      geo_resolution = "gadm_world",
      weight = config$climate_weight,
      weight_year = config$climate_weight_year,
      time_frequency = "yearly",
      verbose = FALSE
    )) %>%
      transmute(year = as.integer(date), value = .data[[value_code]])
  }
  inner_join(
    fetch("avg. temperature", "tmp") %>% rename(TM = value),
    fetch("precipitation", "pre") %>% rename(RR = value),
    by = "year"
  ) %>%
    transmute(gadm_level = "gadm_world", climate_source = config$climate_source_code,
              GID_0 = "WLD", GID_1 = "WLD", year, TM, RR = RR / 1000)
}

# Anomaly of `v` in one unit's series, as add_weather_variables() defines it:
# trailing moments lagged one year, fixed-period moments as they are, and the
# Hamilton trend scaled by the SD of its cycle over the unit's whole series.
unit_anomaly <- function(d, v) {
  if (BASELINE == "trend") {
    centre <- d[[paste0("trend_", v)]]
    scale <- stats::sd(d[[v]] - centre, na.rm = TRUE)
  } else if (startsWith(BASELINE, "lag_")) {
    window <- sub("^lag_", "", BASELINE)
    centre <- lag_by_year(d[[paste0("mean_", v, "_", window)]], d$year)
    scale <- lag_by_year(d[[paste0("sd_", v, "_", window)]], d$year)
  } else {
    centre <- d[[paste0("mean_", v, "_", BASELINE)]]
    scale <- d[[paste0("sd_", v, "_", BASELINE)]]
  }
  anomaly <- d[[v]] - centre
  if (DEVIATION == "standardized") anomaly / scale else anomaly
}

# Levels and anomalies of one resolution, suffixed with `level`.
level_anomalies <- function(climate, level) {
  climate %>%
    group_by(GID_0, GID_1) %>%
    arrange(year, .by_group = TRUE) %>%
    group_modify(~ mutate(.x, zTM = unit_anomaly(.x, "TM"),
                          zRR = unit_anomaly(.x, "RR"))) %>%
    ungroup() %>%
    select(GID_0, GID_1, year, TM, RR, zTM, zRR) %>%
    rename_with(~ paste0(.x, "_", level), c(TM, RR, zTM, zRR))
}

climate_levels <- level_anomalies(load_climate_panel("gadm1", CONFIG), "reg") %>%
  left_join(
    level_anomalies(load_climate_panel("gadm0", CONFIG), "nat") %>%
      select(-GID_1),
    by = c("GID_0", "year")
  ) %>%
  left_join(
    level_anomalies(derive_climate_moments(load_world_series(CONFIG), CONFIG),
                    "glob") %>%
      select(-GID_0, -GID_1),
    by = "year"
  )

# (b) Components of the regional anomaly ---------------------------------------------

for (v in c("TM", "RR")) {
  z <- function(level) climate_levels[[paste0("z", v, "_", level)]]
  climate_levels[[paste0("z", v, "_glob")]] <- z("glob")
  climate_levels[[paste0("z", v, "_nat_res")]] <- z("nat") - z("glob")
  climate_levels[[paste0("z", v, "_reg_res")]] <- z("reg") - z("nat")
}
arrow::write_parquet(climate_levels, file.path(OUT, "climate_levels.parquet"))

countries_without_national <- climate_levels %>%
  group_by(GID_0) %>%
  summarise(national_available = any(!is.na(TM_nat)), .groups = "drop") %>%
  filter(!national_available)
if (nrow(countries_without_national)) {
  warning(nrow(countries_without_national), " GADM1 countries have no GADM0 ",
          "climate series: ", paste(countries_without_national$GID_0,
                                    collapse = ", "))
}

# DOSE panel with the components ------------------------------------------------------

dose <- build_main_dat("DOSE", SPEC, config_args = CONFIG_ARGS,
                       baseline = BASELINE, deviation = DEVIATION) %>%
  left_join(
    climate_levels %>%
      select(GID_0, GID_1, year, TM_nat, RR_nat, TM_glob, RR_glob,
             zTM_reg, zRR_reg, zTM_nat, zRR_nat, matches("_(glob|nat_res|reg_res)$")),
    by = c("GID_0", "GID_1", "year")
  )

# The regional anomaly recomputed here must be build_dat()'s zTM and zRR, or the
# components would not add up to the regressor of the main models.
anomaly_gap <- with(dose, max(abs(c(zTM - zTM_reg, zRR - zRR_reg)), na.rm = TRUE))
if (anomaly_gap > 1e-8) {
  stop("Regional anomalies differ from build_dat()'s by up to ", anomaly_gap, ".")
}

component_cols <- as.vector(outer(c("zTM", "zRR"), names(COMPONENTS), paste,
                                  sep = "_"))
dose_sample <- dose %>%
  filter(!is.na(.data[[SPEC$outcome]]),
         if_all(all_of(c("zTM", "zRR", component_cols)), is.finite))
dropped <- dose %>%
  filter(!is.na(.data[[SPEC$outcome]])) %>%
  anti_join(dose_sample, by = c("GID_1", "year"))
write_out(dropped %>% distinct(GID_0, GID_1, year), "dropped_region_years")
message("DOSE region-years with growth but an incomplete decomposition ",
        "(dropped): ", nrow(dropped))

# How much of the regional anomaly each component carries, on the sample.
component_summary <- bind_rows(lapply(c("TM", "RR"), function(v) {
  total <- dose_sample[[paste0("z", v)]]
  bind_rows(lapply(names(COMPONENTS), function(k) {
    x <- dose_sample[[paste0("z", v, "_", k)]]
    tibble::tibble(
      variable = v, component = COMPONENTS[[k]], column = paste0("z", v, "_", k),
      mean = mean(x), sd = stats::sd(x),
      correlation_with_total = stats::cor(x, total),
      # Covariance shares sum to one across the three components.
      covariance_share = stats::cov(x, total) / stats::var(total)
    )
  }))
}))
component_correlations <- bind_rows(lapply(c("TM", "RR"), function(v) {
  cols <- paste0("z", v, "_", names(COMPONENTS))
  as.data.frame(stats::cor(dose_sample[cols])) %>%
    tibble::rownames_to_column("component") %>%
    mutate(variable = v, .before = 1)
}))
write_out(component_summary, "component_summary")
write_out(component_correlations, "component_correlations")

# (c) Regressions -----------------------------------------------------------------------

component_terms <- function(components) {
  paste(as.vector(outer(c("zTM", "zRR"), components, paste, sep = "_")),
        collapse = " + ")
}
MODELS <- list(
  "Total anomaly, year FE" = list(
    terms = "zTM + zRR", fixed_effects = SPEC$fixed_effects),
  "Components, year FE" = list(
    terms = component_terms(c("nat_res", "reg_res")),
    fixed_effects = SPEC$fixed_effects),
  "Total anomaly, no year FE" = list(
    terms = "zTM + zRR", fixed_effects = REGION_TRENDS),
  "Components, no year FE" = list(
    terms = component_terms(names(COMPONENTS)), fixed_effects = REGION_TRENDS),
  "Components, country-year FE" = list(
    terms = component_terms("reg_res"),
    fixed_effects = paste("GID_0^year +", REGION_TRENDS))
)

models <- lapply(MODELS, function(m) {
  fit_main(m$terms, dose_sample,
           spec = modifyList(SPEC, list(fixed_effects = m$fixed_effects)))
})

term_labels <- function(term) {
  v <- sub("^z(TM|RR).*", "\\1", term)
  k <- sub("^z(TM|RR)_?", "", term)
  tibble::tibble(
    variable = c(TM = "Temperature", RR = "Precipitation")[v],
    component = if_else(nzchar(k), unname(COMPONENTS[k]), "Total")
  )
}
coefficients <- bind_rows(lapply(names(models), function(label) {
  bind_rows(lapply(names(VCOVS), function(vc) {
    tidy_fixest(models[[label]], label, data = dose_sample,
                vcov = VCOVS[[vc]]) %>%
      mutate(vcov = vc, fixed_effects = MODELS[[label]]$fixed_effects)
  }))
})) %>%
  bind_cols(term_labels(.$term))
write_out(coefficients, "coefficients")

# Does growth respond equally to each scale of variation? Pairwise equality of
# the component coefficients, per model and variance estimator.
equality_tests <- bind_rows(lapply(names(models), function(label) {
  b <- names(stats::coef(models[[label]]))
  bind_rows(lapply(c("TM", "RR"), function(v) {
    terms <- intersect(paste0("z", v, "_", names(COMPONENTS)), b)
    if (length(terms) < 2L) return(NULL)
    pairs <- utils::combn(terms, 2L, simplify = FALSE)
    bind_rows(lapply(names(VCOVS), function(vc) {
      m <- summary(models[[label]], vcov = VCOVS[[vc]])
      bind_rows(lapply(pairs, function(p) {
        lincom(m, stats::setNames(c(1, -1), p), paste(p, collapse = " = "))
      })) %>%
        mutate(model = label, variable = v, vcov = vc, .before = 1)
    }))
  }))
}))
write_out(equality_tests, "component_equality_tests")

# Figure --------------------------------------------------------------------------------

model_levels <- names(MODELS)
plot_data <- coefficients %>%
  filter(vcov == "Two-way: GID_0 + year") %>%
  mutate(
    variable = factor(variable, levels = c("Temperature", "Precipitation")),
    component = factor(component, levels = rev(c("Total", COMPONENTS))),
    model = factor(model, levels = model_levels),
    across(c(estimate, conf_low, conf_high), ~ 100 * .x)
  )
x_label <- if (DEVIATION == "standardized") {
  "Growth effect of a one-SD anomaly (percentage points)"
} else {
  "Growth effect of a one-unit anomaly (percentage points; degrees C or metres)"
}
# Reversed dodge, so models stack top to bottom in legend order.
dodge <- position_dodge(width = 0.7, reverse = TRUE)
p_components <- ggplot(plot_data, aes(estimate, component, colour = model,
                                      shape = model)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0,
                linewidth = 0.5, position = dodge) +
  geom_point(size = 2.4, position = dodge) +
  facet_wrap(~variable, scales = "free_x") +
  scale_colour_manual(values = grid_colours(model_levels), name = NULL) +
  scale_shape_manual(values = grid_shapes(model_levels), name = NULL) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  labs(
    title = "DOSE growth response to world, national and regional weather",
    subtitle = paste0("CRU TS, area-weighted; ", DEVIATION, " anomalies ",
                      "against the ", BASELINE, " baseline\n95% CI, two-way ",
                      "clustered by country and year"),
    x = x_label, y = NULL,
    caption = paste0("Region-years: ", format(nrow(dose_sample), big.mark = ","),
                     ". All models include region-specific quadratic trends.")
  ) +
  theme_results() +
  theme(panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "component_coefficients.png"), p_components,
       width = 9, height = 5.5, dpi = 200)

print(component_summary)
print(coefficients %>%
        filter(vcov == "Two-way: GID_0 + year") %>%
        select(model, term, estimate, std_error, p_value, observations))

# (d) Anomalies against levels, with and without country-year FE -------------------
#
# zTM, zRR against the levels TM, RR. The region fixed effects demean the
# levels, so TM gives the same coefficient as TM minus its regional mean or any
# fixed-period mean; the anomaly instead removes the lagged baseline mean, which
# the region trends largely absorb. Each pair is fitted under year FE and under
# country-year FE, which leaves only variation within a country-year. Under
# country-year FE, zTM gives the same coefficient as zTM_reg_res above, since
# the world and national components are absorbed.
FE_SETS <- c("Year FE" = SPEC$fixed_effects,
             "Country-year FE" = paste("GID_0^year +", REGION_TRENDS))
REGRESSORS <- c("Anomaly (zTM, zRR)" = "zTM + zRR", "Level (TM, RR)" = "TM + RR")

level_models <- list()
for (fe in names(FE_SETS)) {
  for (r in names(REGRESSORS)) {
    level_models[[paste0(r, ", ", fe)]] <- fit_main(
      REGRESSORS[[r]], dose_sample,
      spec = modifyList(SPEC, list(fixed_effects = FE_SETS[[fe]]))
    )
  }
}
level_coefficients <- bind_rows(lapply(names(level_models), function(label) {
  bind_rows(lapply(names(VCOVS), function(vc) {
    tidy_fixest(level_models[[label]], label, data = dose_sample,
                vcov = VCOVS[[vc]]) %>%
      mutate(vcov = vc)
  }))
})) %>%
  mutate(
    regressor = sub(", (Year|Country-year) FE$", "", model),
    fixed_effects = sub("^.*, ", "", model),
    variable = if_else(grepl("TM", term), "Temperature", "Precipitation")
  )
write_out(level_coefficients, "anomaly_vs_level_coefficients")

# Variation each regressor keeps once the fixed effects are partialled out, and
# how closely the demeaned anomaly tracks the demeaned level.
within_variation <- bind_rows(lapply(names(FE_SETS), function(fe) {
  demeaning <- fixest::feols(
    stats::as.formula(paste("c(zTM, zRR, TM, RR) ~ 1 |", FE_SETS[[fe]])),
    data = dose_sample
  )
  demeaned <- lapply(as.list(demeaning), stats::resid)
  names(demeaned) <- c("zTM", "zRR", "TM", "RR")
  used <- fixest::obs(demeaning[[1]])
  bind_rows(lapply(names(demeaned), function(x) {
    raw <- dose_sample[[x]][used]
    v <- sub("^z", "", x)
    tibble::tibble(
      fixed_effects = fe, regressor = x,
      sd_raw = stats::sd(raw),
      sd_demeaned = stats::sd(demeaned[[x]]),
      share_of_variance_kept = stats::var(demeaned[[x]]) / stats::var(raw),
      correlation_anomaly_level = stats::cor(demeaned[[paste0("z", v)]],
                                             demeaned[[v]]),
      observations = length(used)
    )
  }))
}))
write_out(within_variation, "anomaly_vs_level_within_variation")

regressor_levels <- names(REGRESSORS)
level_plot_data <- level_coefficients %>%
  filter(vcov == "Two-way: GID_0 + year") %>%
  mutate(
    variable = factor(variable, levels = c("Temperature", "Precipitation")),
    fixed_effects = factor(fixed_effects, levels = rev(names(FE_SETS))),
    regressor = factor(regressor, levels = regressor_levels),
    across(c(estimate, conf_low, conf_high), ~ 100 * .x)
  )
level_x_label <- if (DEVIATION == "standardized") {
  "Growth effect (percentage points; anomaly per SD, level per degree C or metre)"
} else {
  "Growth effect of one degree C or one metre (percentage points)"
}
p_levels <- ggplot(level_plot_data, aes(estimate, fixed_effects,
                                        colour = regressor, shape = regressor)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0,
                linewidth = 0.5, position = dodge) +
  geom_point(size = 2.4, position = dodge) +
  facet_wrap(~variable, scales = "free_x") +
  scale_colour_manual(values = grid_colours(regressor_levels), name = NULL) +
  scale_shape_manual(values = grid_shapes(regressor_levels), name = NULL) +
  labs(
    title = "Weather anomalies against weather levels, DOSE",
    subtitle = paste0("CRU TS, area-weighted; ", DEVIATION, " anomalies ",
                      "against the ", BASELINE, " baseline\n95% CI, two-way ",
                      "clustered by country and year"),
    x = level_x_label, y = NULL,
    caption = "All models include region-specific quadratic trends."
  ) +
  theme_results() +
  theme(panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "anomaly_vs_level_coefficients.png"), p_levels,
       width = 9, height = 4.2, dpi = 200)

print(within_variation)
print(level_coefficients %>%
        filter(vcov == "Two-way: GID_0 + year") %>%
        select(model, term, estimate, std_error, p_value, observations))
