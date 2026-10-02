# Country pieces for the BHM + deviation specification:
#
#   growth ~ TM + TM^2 + RR + RR^2 + zTM^2 + zRR^2
#            | year + unit FE + unit quadratic trends
#
# (main_analysis.Rmd's `model`, with main_spec() fixed effects; zTM, zRR are
# standardized anomalies against the lagged 30-year moments, build_dat()). The
# pooled model is estimated on the dataset's own units: countries for PWT,
# regions for DOSE. The pieces are always countries: on DOSE each piece pools a
# country's regions (its own region-by-year panel), residualised on the pooled
# model's region fixed effects and trends and year effects.
#
# Figure 1, deviation pieces: the pooled quadratic b * zTM^2 with its 95% CI,
# and each country's own fit of the zTM partial residual
# (country_deviation_curves(), functions/country_pieces.R), over the 2nd-98th
# percentile of the country's zTM: a full quadratic a + c1 zTM + c2 zTM^2
# (BASIS = "quadratic", the default) or a natural spline (BASIS = "spline").
# Faceted by Koppen-Geiger climate zone (climate_zones(): each country takes
# the zone of most of its regions; countries without Koppen data are
# "Unclassified"). In each zone the countries departing most from the pooled
# curve (largest Wald chi2 of their deviation from it, clustered by year) are
# labelled at their largest gap.
#
# Figure 2, temperature pieces: plot_bhm_country_pieces.R's figure with
# zTM^2 + zRR^2 among the controls, so the TM pieces are net of the deviation
# terms.
#
# Usage: Rscript plot_bhm_dev_country_pieces.R [PWT|DOSE] [quadratic|spline]
#        (defaults DOSE, quadratic)
# Output: results/bhm_dev_country_pieces/<dataset>/; the deviation figure and
# tables carry the basis in their names.

source("load_functions.R")

args <- commandArgs(trailingOnly = TRUE)
ECON <- if (length(args) >= 1) args[[1]] else "DOSE"
BASIS <- if (length(args) >= 2) args[[2]] else "quadratic"
BASIS <- match.arg(BASIS, c("quadratic", "spline"))
SPEC <- main_spec(econ_year_min = 1900L)
UNIT <- "GID_0"       # Pieces are countries on both panels.
DEV_TERMS <- "zTM^2 + zRR^2"
CONTROLS <- paste("RR + RR^2 +", DEV_TERMS)
SPLINE_DF <- 3
N_LABEL_PER_ZONE <- 3
Z_LIMIT <- 4          # Plotted zTM range, +/- SDs.
OUT <- file.path("results", "bhm_dev_country_pieces", tolower(ECON))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

needed <- c(SPEC$outcome, "TM", "RR", "zTM", "zRR")
panel <- build_main_dat(ECON, SPEC) %>%
  filter(if_all(all_of(needed), is.finite)) %>%
  add_climate_zone(level = UNIT)
stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
national <- all(panel$GID_1 == panel$GID_0)
fe_label <- if (national) {
  "country FE + quadratic trends + year FE"
} else {
  "region FE + quadratic trends + year FE (pieces pool each country's regions)"
}
cat(sprintf("%s: %d unit-years, %d regions, %d countries, %d-%d\n", ECON,
            nrow(panel), n_distinct(panel$GID_1), n_distinct(panel$GID_0),
            min(panel$year), max(panel$year)))
print(panel %>% distinct(GID_0, climate_zone) %>% count(climate_zone))

# Estimation: the BHM + dev model is country_pieces()'s quadratic model --------------

res <- country_pieces(panel, CONTROLS, SPEC, unit = UNIT)
m_dev <- res$m_quad
dev <- country_deviation_curves(panel, m_dev, "zTM", basis = BASIS,
                                df = SPLINE_DF, unit = UNIT)
basis_label <- if (BASIS == "quadratic") {
  "full quadratic a + c₁·zTM + c₂·zTM²"
} else {
  sprintf("natural spline (df = %d)", SPLINE_DF)
}

z_grid <- seq(-Z_LIMIT, Z_LIMIT, by = 0.05)
pooled_dev <- tibble::tibble(z = z_grid) %>%
  bind_cols(linear_combination(m_dev, cbind(`I(zTM^2)` = z_grid^2)))

labelled <- dev$stats %>%
  group_by(climate_zone) %>%
  slice_max(wald_chi2, n = N_LABEL_PER_ZONE) %>%
  ungroup()
ylim <- stats::quantile(dev$curves$estimate[abs(dev$curves$z) <= Z_LIMIT],
                        c(0.01, 0.99))
# Label at the largest gap that is inside the plotted window.
label_points <- dev$curves %>%
  semi_join(labelled, by = "unit") %>%
  filter(abs(z) <= Z_LIMIT, estimate >= ylim[[1]], estimate <= ylim[[2]]) %>%
  group_by(unit) %>%
  slice_max(abs(estimate - pooled), n = 1, with_ties = FALSE) %>%
  ungroup()
zone_notes <- dev$stats %>%
  group_by(climate_zone) %>%
  summarise(note = sprintf("%d countries (%d tested); %d%% depart at p < 0.05",
                           n(), sum(!is.na(p_value)),
                           round(100 * mean(p_value < 0.05, na.rm = TRUE))),
            .groups = "drop")

write_out(dev$stats, paste0("deviation_", BASIS, "_tests"))
write_out(dev$curves, paste0("deviation_", BASIS, "_curves"))
write_out(res$pieces, "temperature_pieces")

# Figure 1: deviation pieces by climate zone --------------------------------------------

p_dev <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_line(data = dev$curves, aes(z, estimate, group = unit),
            colour = "#C0392B", alpha = 0.25, linewidth = 0.35) +
  geom_line(data = semi_join(dev$curves, labelled, by = "unit"),
            aes(z, estimate, group = unit),
            colour = "#7B241C", alpha = 0.9, linewidth = 0.6) +
  geom_ribbon(data = pooled_dev, aes(z, ymin = conf_low, ymax = conf_high),
              fill = "grey20", alpha = 0.2) +
  geom_line(data = pooled_dev, aes(z, estimate), linewidth = 1) +
  ggrepel::geom_label_repel(
    data = label_points, aes(z, estimate, label = unit),
    size = 2.6, colour = "#7B241C", label.size = 0.15, label.padding = 0.12,
    min.segment.length = 0, seed = 1
  ) +
  facet_wrap(~climate_zone, ncol = 2, labeller = as_labeller(
    stats::setNames(paste0(zone_notes$climate_zone, "  —  ",
                           zone_notes$note),
                    zone_notes$climate_zone)
  )) +
  coord_cartesian(xlim = c(-Z_LIMIT, Z_LIMIT), ylim = ylim) +
  labs(
    x = "Standardized temperature deviation zTM (SDs of the lagged 30-year baseline)",
    y = "Growth effect (partial residual)",
    title = sprintf("Deviation response: pooled quadratic vs country %ss (%s)",
                    BASIS, ECON),
    subtitle = sprintf(
      paste0("BHM + dev: TM + TM² + RR + RR² + zTM² + ",
             "zRR², %s.\nBlack: pooled b·zTM² = ",
             "%.4f·zTM², 95%% CI. Red: each country's %s, over ",
             "its 2nd-98th percentile of zTM.\nLabelled: the %d countries ",
             "per zone departing most from the pooled curve (Wald ",
             "χ² of the deviation, clustered by year)."),
      fe_label, dev$b, basis_label, N_LABEL_PER_ZONE
    )
  ) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, paste0("deviation_pieces_by_zone_", BASIS, ".png")),
       p_dev,
       width = 11, height = 10, dpi = 200)

# Figure 2: temperature pieces net of the deviation terms --------------------------

p_tm <- plot_country_pieces(
  res,
  title = sprintf("BHM + dev quadratic and its country pieces (%s)", ECON),
  spec_label = paste0("Specification: TM + TM² + RR + RR² + ",
                      "zTM² + zRR², ", fe_label, ";\npieces net of ",
                      "the deviation terms.")
)
ggsave(file.path(OUT, "bhm_dev_country_pieces.png"), p_tm,
       width = 10, height = 9, dpi = 200)

# Console summary -----------------------------------------------------------------

etable(m_dev, fitstat = ~n + r2 + wr2)
cat(sprintf("Optimum: %.2f C; pooled linear TM slope from pieces: %.5f\n",
            res$t_opt, res$pooled_from_pieces))
cat(sprintf(paste0("\nTemperature pieces: %.0f%% of the %d tested countries ",
                   "depart from the tangent at p < 0.05 (weight-weighted: ",
                   "%.0f%%)\n"),
            100 * mean(res$pieces$gap_p < 0.05, na.rm = TRUE),
            sum(!is.na(res$pieces$gap_p)),
            100 * with(res$pieces, sum(weight[!is.na(gap_p) & gap_p < 0.05]))))
print(res$pieces %>% arrange(desc(weight)) %>%
        transmute(unit, climate_zone, regions, obs,
                  TM_mean = round(TM_mean, 1),
                  slope_pp = round(100 * slope, 2),
                  curve_slope_pp = round(100 * curve_slope, 2),
                  gap_z = round(gap_z, 1), weight = round(weight, 3)),
      n = 12)
cat(sprintf("\nDeviation %ss: share departing from the pooled curve at ",
            BASIS),
    "p < 0.05 (5% expected under the pooled curve)\n")
print(zone_notes)
if (BASIS == "quadratic") {
  # Own coefficients by zone: the pooled model has no linear zTM term, so a
  # systematic own_linear is asymmetry it cannot represent.
  own_significant <- function(est, se, years) {
    2 * stats::pt(-abs(est / se), years - 1) < 0.05
  }
  cat("\nCountry quadratics by zone (tested countries; medians and shares",
      "significant at 5%):\n")
  print(dev$stats %>%
          filter(!is.na(p_value)) %>%
          group_by(climate_zone) %>%
          summarise(
            countries = n(),
            linear_median = median(own_linear),
            linear_pos_sig = mean(own_significant(own_linear, own_linear_se,
                                                  years) & own_linear > 0),
            linear_neg_sig = mean(own_significant(own_linear, own_linear_se,
                                                  years) & own_linear < 0),
            quadratic_median = median(own_quadratic),
            .groups = "drop"
          ) %>%
          mutate(across(where(is.double), ~signif(.x, 3))))
}
print(labelled %>% arrange(climate_zone, desc(wald_chi2)) %>%
        mutate(across(c(z_sd, wald_chi2), ~round(.x, 2)),
               p_value = signif(p_value, 2)),
      n = Inf)
cat("\nWrote figures and tables to ", OUT, "\n", sep = "")
