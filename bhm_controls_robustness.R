# Standard BHM regressions (Burke, Hsiang and Miguel 2015) with the conflict and
# transition-collapse controls and the post-Pinatubo outlier removal: does the
# concave temperature response survive?
#
# Growth of GDP per capita on temperature and precipitation levels, with year
# and country fixed effects and errors clustered by country:
#   A  TM + TM^2 + RR + RR^2, country quadratic trends (BHM's specification)
#   B  the same with country linear trends
#   C  the same without country trends
#   D  TM + TM:mean_TM_all + RR + RR:mean_RR_all, country linear trends (the
#      interacted form of main_analysis.Rmd)
# each in three versions:
#   baseline       no controls
#   controls       + the UCDP/PRIO conflict indicator (add_conflict() in
#                  functions/conflict.R) and the post-Soviet and post-Yugoslav
#                  collapse-year dummies (functions/transition_collapse.R)
#   controls+drop  the controls, without the growth outliers of 1991-1993, the
#                  cool years after the Pinatubo eruption (flag_growth_outliers()
#                  in functions/outlier_sensitivity.R: the 2.5% largest growth
#                  shocks net of the fixed effects with country linear trends
#                  and of the controls, as in test_gap_form_validation.R)
#
# Panels: PWT and WB x ERA5, CRU TS and UDelaware, population-density weights at
# base year 2000, every available year.
#
# Reported: the linear temperature term b1, the concavity term b2 (TM^2, or
# TM:mean_TM_all in D) with clustered standard errors and p-values, the joint
# Wald test of both, and the turning point: the temperature at which the
# marginal effect of temperature is zero, -b1 / (2 b2) in A-C and the long-run
# temperature -b1 / b2 in D.
#
# Output: results/bhm_controls_robustness/ (estimates.csv, dropped_outliers.csv,
# concavity.png)

source("load_functions.R")

OUTPUT_DIR <- file.path("results", "bhm_controls_robustness")
ECON_DATASETS <- c("PWT", "WB")
CLIMATE_SOURCES <- c(ERA5 = "ERA5", CRU = "CRU TS", UDel = "UDelaware")
WEIGHT <- "population density"
WEIGHT_YEAR <- 2000L
CONTROLS <- c("conflict", "i(fsu_collapse_year, ref = 'none')",
              "i(yug_collapse_year, ref = 'none')")
OUTLIER_YEARS <- 1991:1993
OUTLIER_SHARE <- 0.025
OUTLIER_FE <- "year + GID_1 + GID_1[t]"
SPECS <- list(
  `A quadratic trends` = list(terms = "TM + TM2 + RR + RR2", fe = "year + GID_1[t, t2]",
                              concavity = "TM2"),
  `B linear trends` = list(terms = "TM + TM2 + RR + RR2", fe = "year + GID_1[t]",
                           concavity = "TM2"),
  `C no trends` = list(terms = "TM + TM2 + RR + RR2", fe = "year + GID_1",
                       concavity = "TM2"),
  `D interacted` = list(terms = "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
                        fe = "year + GID_1[t]", concavity = "TM:mean_TM_all")
)
VARIANTS <- c("baseline", "controls", "controls+drop")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

# The concavity term as fixest names it (interactions may be stored either way).
term_name <- function(model, term) {
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  candidates <- unique(c(term, paste(rev(parts), collapse = ":")))
  intersect(candidates, names(stats::coef(model)))[[1]]
}

estimates <- list()
dropped <- list()
for (econ in ECON_DATASETS) for (source_name in names(CLIMATE_SOURCES)) {
  message("Fitting ", econ, " x ", source_name)
  cfg <- panel_config(
    econ_year_min = min(load_econ_panel(resolve_econ_source(econ))$year),
    climate_source = CLIMATE_SOURCES[[source_name]],
    climate_weight = WEIGHT,
    climate_weight_year = WEIGHT_YEAR,
    clim_history_years = 1L
  )
  d <- suppressMessages(build_dat(econ, config = cfg))
  d <- d[est_row(d) & !is.na(d$TM) & !is.na(d$RR) & !is.na(d$mean_TM_all), ]
  d <- add_yugoslav_collapse(add_soviet_collapse(add_conflict(d)))
  d$t <- d$year - 1990L
  d$t2 <- d$t^2
  d$TM2 <- d$TM^2
  d$RR2 <- d$RR^2

  outliers <- flag_growth_outliers(d, CONTROLS, OUTLIER_FE, OUTLIER_YEARS, OUTLIER_SHARE)
  dropped[[length(dropped) + 1]] <- mutate(outliers$table, econ = econ,
                                           climate_source = source_name, .before = 1)
  samples <- list(baseline = d, controls = d, `controls+drop` = d[!outliers$flag, ])

  for (spec in names(SPECS)) for (variant in VARIANTS) {
    s <- SPECS[[spec]]
    terms <- if (variant == "baseline") s$terms else paste(c(s$terms, CONTROLS), collapse = " + ")
    m <- fit_panel_model(terms, samples[[variant]], fixed_effects = s$fe, cluster = ~GID_0)
    ct <- fixest::coeftable(m)
    k <- term_name(m, s$concavity)
    b1 <- ct["TM", 1]
    b2 <- ct[k, 1]
    estimates[[length(estimates) + 1]] <- tibble::tibble(
      econ = econ, climate_source = source_name, spec = spec, variant = variant,
      observations = stats::nobs(m), countries = n_distinct(samples[[variant]]$GID_0),
      b1 = b1, b1_se = ct["TM", 2], b1_p = ct["TM", 4],
      b2 = b2, b2_se = ct[k, 2], b2_p = ct[k, 4],
      joint_p = fixest::wald(m, keep = paste0("^(TM|", k, ")$"), print = FALSE)$p,
      turning_point = if (s$concavity == "TM2") -b1 / (2 * b2) else -b1 / b2
    )
  }
}
estimates <- bind_rows(estimates)
dropped <- bind_rows(dropped)
write_out(estimates, "estimates")
write_out(dropped, "dropped_outliers")

options(width = 220)
cat("\n== Outliers dropped (1991-1993) ==\n")
print(as.data.frame(dropped %>% count(econ, climate_source, name = "dropped")), row.names = FALSE)
cat("\n== Temperature terms (pp per degree C; b2 per degree C squared, or per degree C",
    "of long-run temperature in D) ==\n")
print(as.data.frame(estimates %>%
  transmute(econ, climate_source, spec, variant, observations,
            b1 = round(100 * b1, 3), b1_p = signif(b1_p, 2),
            b2 = round(100 * b2, 4), b2_p = signif(b2_p, 2),
            joint_p = signif(joint_p, 2), turning_point = round(turning_point, 1))),
  row.names = FALSE)

# Plot: the concavity term by panel and version ---------------------------------------

variant_labels <- c(baseline = "No controls", controls = "Conflict and collapse controls",
                    `controls+drop` = "Controls, 1991-93 outliers dropped")
plot_data <- estimates %>%
  mutate(
    panel = factor(paste(econ, climate_source, sep = " x "),
                   levels = rev(unique(paste(econ, climate_source, sep = " x ")))),
    variant = factor(variant_labels[variant], levels = variant_labels),
    spec = factor(spec, levels = names(SPECS)),
    estimate = 100 * b2,
    conf_low = 100 * (b2 - stats::qnorm(0.975) * b2_se),
    conf_high = 100 * (b2 + stats::qnorm(0.975) * b2_se)
  )
dodge <- position_dodge(width = 0.6)
p <- ggplot(plot_data, aes(estimate, panel, colour = variant, shape = variant)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high), position = dodge, linewidth = 0.5) +
  geom_point(position = dodge, size = 2) +
  facet_wrap(~spec, nrow = 1, scales = "free_x") +
  scale_x_continuous(n.breaks = 4) +
  scale_colour_manual(values = grid_colours(variant_labels), name = NULL) +
  scale_shape_manual(values = grid_shapes(variant_labels), name = NULL) +
  labs(
    x = "Concavity term (TM^2 in A-C, TM x long-run temperature in D), pp; 95% CI",
    y = NULL,
    title = "BHM temperature concavity with conflict and collapse controls",
    subtitle = paste0("Controls: UCDP/PRIO conflict, post-Soviet and post-Yugoslav collapse-year ",
                      "dummies; dropped: the 2.5% largest growth shocks of 1991-1993."),
    caption = paste0("Population-density weights at base year ", WEIGHT_YEAR,
                     ". Fixed effects year + country (+ country trends); errors clustered ",
                     "by country.")
  ) +
  theme_results() +
  theme(panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.spacing.x = grid::unit(1.5, "lines"))
ggsave(file.path(OUTPUT_DIR, "concavity.png"), p, width = 13, height = 5.5, dpi = 300)

message("Results written to ", normalizePath(OUTPUT_DIR))
