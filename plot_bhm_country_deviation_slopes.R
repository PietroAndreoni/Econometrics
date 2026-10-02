# BHM curve shared across countries, deviation curvature country-specific:
#
#   growth ~ TM + TM^2 + RR + RR^2 + i(country, zTM^2) + i(country, zRR^2)
#            | year + unit FE + unit quadratic trends          (main_spec())
#
# One joint regression: every country has its own zTM^2 and zRR^2 coefficient,
# while TM, TM^2, RR, RR^2 are common, so all countries share the BHM optimum.
# Only the deviation curvature varies by country. A country-specific linear
# zTM is not included: within a country, after its fixed effect and trends,
# zTM is almost a rescaled copy of TM (median correlation 0.95 on PWT, 0.975
# on DOSE), so a country zTM slope absorbs the country's own temperature slope
# and the BHM curve is no longer identified (the optimum moves to 1.7 C on
# PWT and -83 C on DOSE).
#
# Inference: the shared BHM terms are clustered as in main_spec() (by GID_1).
# Each country coefficient is identified within its own country, which on the
# national panel is its own GID_1 cluster, so a GID_1-clustered SE for it is
# degenerate; country coefficients use SEs clustered by year instead (all
# units hit by one year's shock count once), as in the per-country tests of
# plot_bhm_dev_country_pieces.R. Departures from the pooled coefficient are
# t-tested against t(years - 1); countries with fewer than MIN_YEARS years get
# no test. eb_shrink() (functions/inference.R) gives the implied true
# cross-country SD tau of each coefficient and shrunk country estimates.
#
# Figures: the shared vs pooled BHM curve; country zTM^2 and zRR^2 curves
# c_i * z^2 by Koppen-Geiger zone against the pooled b * z^2; and a forest plot
# of the country coefficients (raw with 95% CI, and EB-shrunk).
#
# Usage: Rscript plot_bhm_country_deviation_slopes.R [PWT|DOSE] (default DOSE)
# Output: results/bhm_country_deviation_slopes/<dataset>/.

source("load_functions.R")

args <- commandArgs(trailingOnly = TRUE)
ECON <- if (length(args) >= 1) args[[1]] else "DOSE"
SPEC <- main_spec(econ_year_min = 1900L)
BASE <- "TM + TM^2 + RR + RR^2"
MIN_YEARS <- 10
N_LABEL_PER_ZONE <- 3
Z_LIMIT <- 4          # Plotted deviation range, +/- SDs.
OUT <- file.path("results", "bhm_country_deviation_slopes", tolower(ECON))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}
DEV_LABELS <- c(zTM = "Temperature deviation zTM",
                zRR = "Precipitation deviation zRR")

needed <- c(SPEC$outcome, "TM", "RR", "zTM", "zRR")
panel <- build_main_dat(ECON, SPEC) %>%
  filter(if_all(all_of(needed), is.finite)) %>%
  add_climate_zone(level = "GID_0") %>%
  mutate(zTM2 = zTM^2, zRR2 = zRR^2)
stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
national <- all(panel$GID_1 == panel$GID_0)
fe_label <- if (national) "country FE + quadratic trends + year FE" else
  "region FE + quadratic trends + year FE"
cat(sprintf("%s: %d unit-years, %d regions, %d countries, %d-%d\n", ECON,
            nrow(panel), n_distinct(panel$GID_1), n_distinct(panel$GID_0),
            min(panel$year), max(panel$year)))

# Estimation ----------------------------------------------------------------------

m_pooled <- fit_main(paste(BASE, "+ zTM^2 + zRR^2"), panel, SPEC)
m_country <- fit_main(paste(BASE, "+ i(GID_0, zTM2) + i(GID_0, zRR2)"),
                      panel, SPEC)
stopifnot(nobs(m_pooled) == nrow(panel), nobs(m_country) == nrow(panel))
if (length(m_country$collin.var)) {
  cat("Dropped as collinear:", m_country$collin.var, "\n")
}
# With ~250-360 coefficients and 59-178 clusters both clustered VCOVs are rank
# deficient; they are PSD by construction, so fixest's eigenvalue "fix" (which
# perturbs every entry) is switched off.
V_year <- stats::vcov(m_country, cluster = ~year, vcov_fix = FALSE)
V_bhm <- stats::vcov(m_country, vcov_fix = FALSE)
pooled_b <- c(zTM = coef(m_pooled)[["I(zTM^2)"]],
              zRR = coef(m_pooled)[["I(zRR^2)"]])

country_info <- panel %>%
  group_by(GID_0, climate_zone) %>%
  summarise(years = n_distinct(year), regions = n_distinct(GID_1),
            TM_mean = mean(TM),
            zTM_lo = stats::quantile(zTM, 0.02),
            zTM_hi = stats::quantile(zTM, 0.98),
            zRR_lo = stats::quantile(zRR, 0.02),
            zRR_hi = stats::quantile(zRR, 0.98), .groups = "drop")

coefs <- bind_rows(lapply(names(DEV_LABELS), function(term) {
  pattern <- sprintf("^GID_0::(.+):%s2$", term)
  names_i <- grep(pattern, names(coef(m_country)), value = TRUE)
  tibble::tibble(
    term = term,
    GID_0 = sub(pattern, "\\1", names_i),
    estimate = unname(coef(m_country)[names_i]),
    std_error = sqrt(diag(V_year)[names_i])
  )
})) %>%
  left_join(country_info, by = "GID_0") %>%
  mutate(
    pooled = pooled_b[term],
    tested = years >= MIN_YEARS,
    gap_t = (estimate - pooled) / std_error,
    gap_p = if_else(tested, 2 * stats::pt(-abs(gap_t), years - 1),
                    NA_real_),
    conf_low = estimate - stats::qt(0.975, years - 1) * std_error,
    conf_high = estimate + stats::qt(0.975, years - 1) * std_error
  ) %>%
  group_by(term) %>%
  group_modify(~ bind_cols(.x, eb_shrink(.x$estimate, .x$std_error,
                                         use = .x$tested))) %>%
  ungroup()

heterogeneity <- coefs %>%
  filter(tested) %>%
  group_by(term) %>%
  summarise(
    countries = n(),
    pooled = first(pooled),
    precision_weighted_mean = first(prior_mean),
    raw_sd = stats::sd(estimate),
    mean_se = sqrt(mean(std_error^2)),
    tau = sqrt(first(tau2)),
    share_depart = mean(gap_p < 0.05),
    .groups = "drop"
  )
zone_heterogeneity <- coefs %>%
  filter(tested) %>%
  group_by(term, climate_zone) %>%
  summarise(countries = n(), median = stats::median(estimate),
            shrunk_median = stats::median(shrunk),
            share_depart = mean(gap_p < 0.05), .groups = "drop")

write_out(coefs, "country_deviation_coefficients")
write_out(heterogeneity, "heterogeneity")
write_out(zone_heterogeneity, "heterogeneity_by_zone")

# Figure 1: shared vs pooled BHM curve ---------------------------------------------

t_grid <- seq(floor(min(panel$TM)), ceiling(max(panel$TM)), by = 0.25)
bhm_curve <- function(model, label, vcov = NULL) {
  b <- coef(model)
  t_opt <- -b[["TM"]] / (2 * b[["I(TM^2)"]])
  tibble::tibble(TM = t_grid) %>%
    bind_cols(linear_combination(
      model, cbind(TM = t_grid - t_opt, `I(TM^2)` = t_grid^2 - t_opt^2),
      vcov = vcov
    )) %>%
    mutate(model = sprintf("%s (optimum %.1f°C)", label, t_opt))
}
bhm <- bind_rows(
  bhm_curve(m_pooled, "Pooled deviation terms"),
  bhm_curve(m_country, "Country-specific deviation terms", V_bhm)
) %>%
  mutate(model = factor(model, levels = unique(model)))

p_bhm <- ggplot(bhm, aes(TM, estimate, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("grey20", "#C0392B"),
                      aesthetics = c("colour", "fill")) +
  labs(x = "Temperature (°C)", y = "Growth relative to optimum",
       colour = NULL, fill = NULL,
       title = sprintf("Shared BHM curve with country-specific deviation terms (%s)",
                       ECON),
       subtitle = sprintf(
         paste0("TM + TM² + RR + RR² + zTM² + zRR², %s. ",
                "95%% CI clustered by %s."),
         fe_label, if (national) "country" else "region")) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "bhm_shared_vs_pooled.png"), p_bhm,
       width = 9, height = 5.5, dpi = 200)

# Figure 2: country deviation curves by zone, one figure per term ---------------------

plot_deviation_curves <- function(term) {
  d <- filter(coefs, term == !!term)
  lo <- paste0(term, "_lo")
  hi <- paste0(term, "_hi")
  curves <- d %>%
    rowwise() %>%
    reframe(GID_0, climate_zone, estimate_c = estimate,
            z = seq(.data[[lo]], .data[[hi]], length.out = 40)) %>%
    mutate(estimate = estimate_c * z^2)
  z_grid <- seq(-Z_LIMIT, Z_LIMIT, by = 0.05)
  pooled_curve <- tibble::tibble(z = z_grid) %>%
    bind_cols(linear_combination(
      m_pooled, matrix(z_grid^2, ncol = 1,
                       dimnames = list(NULL, paste0("I(", term, "^2)")))
    ))
  labelled <- d %>%
    filter(tested) %>%
    group_by(climate_zone) %>%
    slice_max(abs(gap_t), n = N_LABEL_PER_ZONE) %>%
    ungroup()
  ylim <- stats::quantile(curves$estimate[abs(curves$z) <= Z_LIMIT],
                          c(0.01, 0.99))
  label_points <- curves %>%
    semi_join(labelled, by = "GID_0") %>%
    filter(abs(z) <= Z_LIMIT, estimate >= ylim[[1]], estimate <= ylim[[2]]) %>%
    group_by(GID_0) %>%
    slice_max(abs(z), n = 1, with_ties = FALSE) %>%
    ungroup()
  notes <- d %>%
    group_by(climate_zone) %>%
    summarise(note = sprintf(
      "%d countries (%d tested); %d%% depart at p < 0.05",
      n(), sum(tested), round(100 * mean(gap_p < 0.05, na.rm = TRUE))),
      .groups = "drop")
  h <- filter(heterogeneity, term == !!term)

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_line(data = curves, aes(z, estimate, group = GID_0),
              colour = "#C0392B", alpha = 0.25, linewidth = 0.35) +
    geom_line(data = semi_join(curves, labelled, by = "GID_0"),
              aes(z, estimate, group = GID_0),
              colour = "#7B241C", linewidth = 0.6) +
    geom_ribbon(data = pooled_curve,
                aes(z, ymin = conf_low, ymax = conf_high),
                fill = "grey20", alpha = 0.2) +
    geom_line(data = pooled_curve, aes(z, estimate), linewidth = 1) +
    ggrepel::geom_label_repel(
      data = label_points, aes(z, estimate, label = GID_0),
      size = 2.6, colour = "#7B241C", label.size = 0.15,
      label.padding = 0.12, min.segment.length = 0, seed = 1
    ) +
    facet_wrap(~climate_zone, ncol = 2, labeller = as_labeller(
      stats::setNames(paste0(notes$climate_zone, "  —  ", notes$note),
                      notes$climate_zone)
    )) +
    coord_cartesian(xlim = c(-Z_LIMIT, Z_LIMIT), ylim = ylim) +
    labs(
      x = paste(DEV_LABELS[[term]], "(SDs of the lagged 30-year baseline)"),
      y = "Growth effect",
      title = sprintf("Country-specific %s² response with a shared BHM curve (%s)",
                      term, ECON),
      subtitle = sprintf(
        paste0("Joint model: TM + TM² + RR + RR² + ",
               "i(country, zTM²) + i(country, zRR²), %s.\n",
               "Black: pooled-model b·%s² = %.4f·%s², 95%% CI. ",
               "Red: each country's c·%s² over its 2nd-98th percentile ",
               "of %s.\nImplied true cross-country SD of c: τ = %.4f ",
               "(raw SD %.4f, mean SE %.4f). Labelled: largest |c - b| / SE ",
               "per zone (SE clustered by year)."),
        fe_label, term, pooled_b[[term]], term, term, term,
        h$tau, h$raw_sd, h$mean_se
      )
    ) +
    theme_classic() +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}
for (term in names(DEV_LABELS)) {
  ggsave(file.path(OUT, paste0("country_", term, "2_curves_by_zone.png")),
         plot_deviation_curves(term), width = 11, height = 10, dpi = 200)
}

# Figure 3: forest plot of the country coefficients ------------------------------

forest <- coefs %>%
  filter(tested) %>%
  mutate(term = factor(DEV_LABELS[term], levels = DEV_LABELS)) %>%
  group_by(term) %>%
  arrange(climate_zone, estimate, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  # One y slot per term and rank, so each panel keeps its own sort order.
  mutate(row_id = factor(sprintf("%s_%03d", term, rank)))
row_labels <- stats::setNames(forest$GID_0, as.character(forest$row_id))
ref_lines <- tibble::tibble(term = factor(DEV_LABELS, levels = DEV_LABELS),
                            pooled = pooled_b)
p_forest <- ggplot(forest, aes(y = row_id)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(data = ref_lines, aes(xintercept = pooled), linewidth = 0.6) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high, colour = climate_zone),
                 alpha = 0.45, linewidth = 0.4) +
  geom_point(aes(x = estimate, colour = climate_zone), size = 0.9) +
  geom_point(aes(x = shrunk), shape = 4, size = 0.9, colour = "grey15") +
  facet_wrap(~term, scales = "free") +
  scale_y_discrete(labels = row_labels) +
  coord_cartesian(xlim = stats::quantile(c(forest$conf_low, forest$conf_high),
                                         c(0.02, 0.98))) +
  labs(x = "Country coefficient on the squared deviation", y = NULL,
       colour = "Climate zone",
       title = sprintf("Country deviation coefficients (%s)", ECON),
       subtitle = paste0("Dots and bars: estimate and 95% CI (clustered by ",
                         "year), sorted within zone; crosses: empirical-",
                         "Bayes shrunk estimate.\nVertical line: pooled ",
                         "coefficient. Clipped to the 2nd-98th percentile ",
                         "of the interval ends.")) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom",
        axis.text.y = element_text(size = 5))
ggsave(file.path(OUT, "country_coefficients_forest.png"), p_forest,
       width = 10, height = max(9, 2.5 + 0.075 * nrow(forest) / 2), dpi = 200,
       limitsize = FALSE)

# Console summary -----------------------------------------------------------------

options(width = 160)
bhm_terms <- c("TM", "I(TM^2)", "RR", "I(RR^2)")
print(tibble::tibble(
  term = bhm_terms,
  pooled = coef(m_pooled)[bhm_terms],
  pooled_se = sqrt(diag(vcov(m_pooled))[bhm_terms]),
  shared = coef(m_country)[bhm_terms],
  shared_se = sqrt(diag(V_bhm)[bhm_terms])
) %>% mutate(across(where(is.double), ~signif(.x, 3))))
cat(sprintf("Pooled deviation coefficients: zTM^2 %.5f (%.5f), zRR^2 %.5f (%.5f)\n",
            pooled_b[["zTM"]], sqrt(vcov(m_pooled)["I(zTM^2)", "I(zTM^2)"]),
            pooled_b[["zRR"]], sqrt(vcov(m_pooled)["I(zRR^2)", "I(zRR^2)"])))
opt <- function(m) -coef(m)[["TM"]] / (2 * coef(m)[["I(TM^2)"]])
cat(sprintf("\nOptimum: pooled %.2f C, shared with country deviations %.2f C\n",
            opt(m_pooled), opt(m_country)))
cat("\nCross-country heterogeneity of the deviation coefficients",
    "(tested countries; 5% departures expected under a common coefficient):\n")
print(heterogeneity %>% mutate(across(where(is.double), ~signif(.x, 3))))
print(zone_heterogeneity %>% mutate(across(where(is.double), ~signif(.x, 3))),
      n = Inf)
cat("\nLargest departures from the pooled coefficient:\n")
print(coefs %>% filter(tested) %>% group_by(term) %>%
        slice_max(abs(gap_t), n = 6) %>% ungroup() %>%
        transmute(term, GID_0, climate_zone, years, regions,
                  estimate = signif(estimate, 3), se = signif(std_error, 3),
                  shrunk = signif(shrunk, 3), gap_t = round(gap_t, 1)),
      n = Inf)
cat("\nWrote figures and tables to ", OUT, "\n", sep = "")
