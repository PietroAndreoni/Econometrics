# Do the slow adaptation gaps need warm and cold arms, or is a BHM-style
# signed term enough? Decomposed adaptation-shock model (see
# test_gap_shock_trends.R): the gap is the Hamilton trend minus the trailing
# 30-year mean, scaled by the SD over the entire climate series; the shocks
# are the squared Hamilton cycles (symmetric). Forms of the gap terms, the
# same for temperature and precipitation, each interacted with the long-run
# climate m:
#
#   signed           b * gap + g * gap * m (BHM style: a warm and a cold gap
#                    of the same size have opposite effects)
#   linear arms      warm and cold gaps, max(gap, 0) and max(-gap, 0), each
#                    with its own slope and interaction; nests the signed
#                    form (warm + cold coefficients = 0)
#   quadratic arms   the squared warm and cold gaps (the form used so far)
#   signed + quad    the signed and the quadratic-arm terms together, nesting
#                    both
#
# Tests, for temperature and precipitation separately: symmetry of the linear
# arms (the signed restriction, 2 restrictions), the quadratic arms given the
# signed terms and the signed terms given the quadratic arms (in the
# encompassing model); fit (AIC, BIC) of each form against the signed one.
# Response of growth to a gap of -2 to 2 SD at long-run temperatures of 5, 15
# and 25 degrees Celsius.
#
# Grid: WB and PWT x ERA5, CRU TS and UDelaware, area and population weights,
# with and without unit trends. Errors clustered by country.
#
# Output: results/gap_form/ (tests.csv, gap_curves.csv,
# gap_curves_<econ>_<source>_<weight>.png for PWT x UDelaware)

source("load_functions.R")

ECON_DATASETS <- c(WB = "WB", PWT = "PWT")
CLIMATE_SOURCES <- c(ERA5 = "ERA5", `CRU TS` = "CRU TS", UDelaware = "UDelaware")
CLIMATE_WEIGHTS <- c(area = "area", pop = "concurrent population")
TRENDS <- c(trends = "year + GID_1 + GID_1[year]", notrends = "year + GID_1")
ARM_BASELINE <- "lag_30"
ARM_SD <- "all"
# Run option: "full" (default) fits the Hamilton trend over each unit's entire
# climate series; "realtime" fits each year's trend on the 30 years up to the
# forecast origin only (panel_config(hamilton_window = 30)), output in a
# "_realtime30" folder. Rscript test_gap_form.R realtime
HAMILTON <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(HAMILTON)) HAMILTON <- "full"
stopifnot(HAMILTON %in% c("full", "realtime"))
HAMILTON_WINDOW <- if (HAMILTON == "realtime") 30L else NULL
CONF_LEVEL <- 0.95
AT_TM <- c(5, 15, 25)
GAPS <- seq(-2, 2, by = 0.05)
PLOT_CASES <- list(c(econ = "PWT", source = "UDelaware"))
CLIMATE <- c(TM = "mean_TM_all", RR = "mean_RR_all")
DEG_C <- "°C"

# Gap terms of each form, by suffix: the term's value at a gap of d SD.
TERM_VALUE <- list(
  s = function(d) d,
  p = function(d) pmax(d, 0),
  m = function(d) pmax(-d, 0),
  p2 = function(d) pmax(d, 0)^2,
  m2 = function(d) pmax(-d, 0)^2
)
FORMS <- list(
  signed = "s",
  linear_arms = c("p", "m"),
  quadratic_arms = c("p2", "m2"),
  signed_quadratic = c("s", "p2", "m2")
)
FORM_LABELS <- c(signed = "Signed (BHM style)", linear_arms = "Linear arms",
                 quadratic_arms = "Quadratic arms",
                 signed_quadratic = "Signed + quadratic arms")
SHOCKS <- "sTM2 + sRR2"

OUTPUT_DIR <- file.path("results", paste0("gap_form", if (HAMILTON == "realtime") "_realtime30" else ""))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit <- function(terms, data, fe) {
  fit_panel_model(terms, data, fixed_effects = fe, cluster = ~GID_0)
}
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
gap_vars <- function(variable, suffixes) paste0("g", variable, suffixes)
form_terms <- function(suffixes) {
  vars <- unlist(lapply(names(CLIMATE), function(v) {
    x <- gap_vars(v, suffixes)
    c(x, paste0(x, ":", CLIMATE[[v]]))
  }))
  paste(c(vars, SHOCKS), collapse = " + ")
}
# Wald test of the gap terms of one variable with the given suffixes (main
# terms and climate interactions).
block_p <- function(model, variable, suffixes) {
  terms <- names(stats::coef(model))
  x <- gap_vars(variable, suffixes)
  keep <- c(x, vapply(x, interaction_name, character(1), b = CLIMATE[[variable]],
                      terms = terms))
  fixest::wald(model, keep = paste0("^(", paste(keep, collapse = "|"), ")$"),
               print = FALSE)$p
}
# Wald test that the warm and cold linear arms sum to zero, main term and
# interaction: the signed (BHM-style) restriction.
symmetry_p <- function(model, variable) {
  b <- stats::coef(model)
  V <- stats::vcov(model)
  terms <- names(b)
  up <- paste0("g", variable, "p")
  down <- paste0("g", variable, "m")
  climate <- CLIMATE[[variable]]
  R <- matrix(0, 2, length(b), dimnames = list(NULL, terms))
  R[1, c(up, down)] <- 1
  R[2, c(interaction_name(up, climate, terms), interaction_name(down, climate, terms))] <- 1
  r <- R %*% b
  stat <- as.numeric(t(r) %*% solve(R %*% V %*% t(R)) %*% r)
  stats::pchisq(stat, df = 2, lower.tail = FALSE)
}
# Effect on growth of a temperature gap of d SD at long-run temperature m.
gap_curve <- function(model, suffixes, gaps, m) {
  terms <- names(stats::coef(model))
  G <- t(vapply(gaps, function(d) {
    w <- stats::setNames(numeric(length(terms)), terms)
    for (s in suffixes) {
      x <- paste0("gTM", s)
      w[[x]] <- TERM_VALUE[[s]](d)
      w[[interaction_name(x, CLIMATE[["TM"]], terms)]] <- TERM_VALUE[[s]](d) * m
    }
    w
  }, numeric(length(terms))))
  colnames(G) <- terms
  linear_combination(model, G, CONF_LEVEL) %>% mutate(gap = gaps)
}

tests <- list()
curves <- list()
for (w in names(CLIMATE_WEIGHTS)) for (econ in names(ECON_DATASETS)) {
  for (source_name in names(CLIMATE_SOURCES)) {
    message("Fitting ", econ, " x ", source_name, ", ", w, " weights")
    cfg <- panel_config(
      econ_year_min = min(load_econ_panel(resolve_econ_source(ECON_DATASETS[[econ]]))$year),
      climate_source = CLIMATE_SOURCES[[source_name]],
      climate_weight = CLIMATE_WEIGHTS[[w]],
      hamilton_window = HAMILTON_WINDOW,
      clim_history_years = 1L
    )
    d <- suppressMessages(build_dat(ECON_DATASETS[[econ]], config = cfg,
                                    baseline = ARM_BASELINE, sd_baseline = ARM_SD))
    for (v in names(CLIMATE)) {
      gap <- (d[[paste0("trend_", v)]] - d[[paste0("mean_", v, "_lag_30")]]) /
        d[[paste0("sd_", v, "_all")]]
      cycle <- (d[[v]] - d[[paste0("trend_", v)]]) / d[[paste0("sd_", v, "_all")]]
      for (s in names(TERM_VALUE)) d[[paste0("g", v, s)]] <- TERM_VALUE[[s]](gap)
      d[[paste0("s", v, "2")]] <- cycle^2
    }

    for (tr in names(TRENDS)) {
      models <- lapply(FORMS, function(f) fit(form_terms(f), d, TRENDS[[tr]]))
      n <- vapply(models, nobs, numeric(1))
      stopifnot(length(unique(n)) == 1)
      aic <- vapply(models, stats::AIC, numeric(1))
      bic <- vapply(models, stats::BIC, numeric(1))
      row <- tibble::tibble(
        weight = w, econ = econ, climate_source = source_name, trends = tr,
        observations = n[[1]],
        aic_linear_arms = aic[["linear_arms"]] - aic[["signed"]],
        aic_quadratic_arms = aic[["quadratic_arms"]] - aic[["signed"]],
        aic_signed_quadratic = aic[["signed_quadratic"]] - aic[["signed"]],
        bic_linear_arms = bic[["linear_arms"]] - bic[["signed"]],
        bic_quadratic_arms = bic[["quadratic_arms"]] - bic[["signed"]],
        bic_signed_quadratic = bic[["signed_quadratic"]] - bic[["signed"]]
      )
      for (v in names(CLIMATE)) {
        row[[paste0(v, "_signed_p")]] <- block_p(models$signed, v, "s")
        row[[paste0(v, "_symmetry_p")]] <- symmetry_p(models$linear_arms, v)
        row[[paste0(v, "_quad_given_signed_p")]] <- block_p(models$signed_quadratic, v,
                                                             c("p2", "m2"))
        row[[paste0(v, "_signed_given_quad_p")]] <- block_p(models$signed_quadratic, v, "s")
      }
      tests[[length(tests) + 1]] <- row
      for (f in names(FORMS)) for (m in AT_TM) {
        curves[[length(curves) + 1]] <- gap_curve(models[[f]], FORMS[[f]], GAPS, m) %>%
          mutate(weight = w, econ = econ, climate_source = source_name, trends = tr,
                 form = f, climate = m)
      }
    }
  }
}

tests <- bind_rows(tests)
curves <- bind_rows(curves)
utils::write.csv(tests, file.path(OUTPUT_DIR, "tests.csv"), row.names = FALSE)
utils::write.csv(curves, file.path(OUTPUT_DIR, "gap_curves.csv"), row.names = FALSE)

options(width = 220)
cat("\n== Fit against the signed form (negative = better than signed) ==\n")
print(as.data.frame(tests %>%
  select(trends, weight, econ, climate_source, starts_with("aic_"), starts_with("bic_")) %>%
  arrange(trends, weight, econ, climate_source) %>%
  mutate(across(where(is.double), ~ round(.x, 1)))), row.names = FALSE)
cat("\n== Wald tests: signed alone, symmetry of linear arms, encompassing model ==\n")
print(as.data.frame(tests %>%
  select(trends, weight, econ, climate_source, ends_with("_p")) %>%
  arrange(trends, weight, econ, climate_source) %>%
  mutate(across(where(is.double), ~ signif(.x, 2)))), row.names = FALSE)

# Plots ---------------------------------------------------------------------------------

form_levels <- unname(FORM_LABELS)
colours <- grid_colours(form_levels)
for (case in PLOT_CASES) for (w in names(CLIMATE_WEIGHTS)) {
  d <- curves %>%
    filter(econ == case[["econ"]], climate_source == case[["source"]], weight == w) %>%
    mutate(
      form = factor(FORM_LABELS[form], levels = form_levels),
      trends = factor(ifelse(trends == "trends", "With unit trends", "Without unit trends"),
                      levels = c("With unit trends", "Without unit trends")),
      climate = factor(paste0("Long-run temperature ", climate, " ", DEG_C),
                       levels = paste0("Long-run temperature ", AT_TM, " ", DEG_C)),
      across(c(estimate, conf_low, conf_high), ~ 100 * .x)
    )
  p <- ggplot(d, aes(x = gap, y = estimate, colour = form, fill = form)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(trends ~ climate) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    labs(
      x = "Temperature gap: Hamilton trend minus trailing 30-year mean (SD)",
      y = "Effect on growth that year (pp)",
      title = paste0("Form of the temperature gap: ", case[["econ"]], " × ", case[["source"]],
                     " (", w, " weights)"),
      subtitle = paste("Gap terms interacted with the long-run temperature, evaluated at",
                       "5, 15 and 25 °C; precipitation gaps in the same form; symmetric",
                       "Hamilton shocks. Bands: 95% CI."),
      caption = "Fixed effects year + country (+ country trends); errors clustered by country."
    ) +
    theme_results()
  ggsave(file.path(OUTPUT_DIR, paste0("gap_curves_", case[["econ"]], "_",
                                      gsub("[^A-Za-z0-9]+", "", case[["source"]]), "_", w,
                                      ".png")),
         p, width = 12, height = 8, dpi = 300)
}
message("Results written to ", normalizePath(OUTPUT_DIR))
