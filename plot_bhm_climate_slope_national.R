# The binned BHM climate interaction, base function only, on national PWT data.
#
# This repeats the climate-bin exercise of plot_bhm_climate_slope_bins.R with
# two restrictions. First, only the base climate block is estimated: no
# temperature or precipitation deviation term at all, so the temperature slope
# is not competing with a separate anomaly response. Second, the primary panel
# is national PWT data, where GID_1 equals GID_0 and clustering is therefore at
# the country level by construction. The DOSE subnational panel is estimated
# the same way alongside it for comparison.
#
# The block being relaxed is
#
#   TM + TM:mean_TM_all + RR + RR:mean_RR_all
#
# in which the marginal effect of a within-unit temperature deviation is forced
# to vary linearly with the unit's long-run mean temperature. Replacing
# `TM:mean_TM_all` with a free slope per bin of `mean_TM_all` removes that
# restriction.
#
# Equal-count bins compress badly at the hot end of the national distribution,
# where many tropical countries sit within a degree or two of each other. Each
# bin is drawn across the range of `mean_TM_all` it covers so that the width is
# visible rather than implicit.
#
# Prerequisite: the notebook's parquet inputs under data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
})

source("rmd_chunks.R")

out <- "results/bhm_climate_slope_national"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

PP <- 100
BIN_PATTERN <- "^mean_TM_bin::[0-9]+:TM$"
BIN_WIDTH <- 2          # Degrees, for the equal-width scheme.
THIN_UNITS <- 10        # Bins with fewer units than this are flagged as thin.

# Equal-count bins put the same number of observations in each bin, which
# compresses them wherever units pile up; equal-width bins keep the degree
# interval constant and let the counts fall where they may.
#
# Equal-width bins leave very thin tails. Their coefficients carry almost no
# information yet still enter the joint Wald test, and a cluster-robust Wald
# test over many restrictions with few clusters over-rejects. The pooled variant
# merges each tail inward until it holds at least `THIN_UNITS` units, which
# keeps the interior grid equally spaced and makes the test interpretable.
pool_tails <- function(x, breaks, unit, min_units) {
  index <- as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
  k <- length(breaks) - 1L
  lo <- 1L
  while (lo < k && n_distinct(unit[index <= lo]) < min_units) lo <- lo + 1L
  hi <- k
  while (hi > lo && n_distinct(unit[index >= hi]) < min_units) hi <- hi - 1L
  c(-Inf, breaks[(lo + 1L):hi], Inf)
}

bin_schemes <- function(d) {
  x <- d$mean_TM_all
  width_breaks <- seq(
    floor(min(x, na.rm = TRUE) / BIN_WIDTH) * BIN_WIDTH,
    ceiling(max(x, na.rm = TRUE) / BIN_WIDTH) * BIN_WIDTH,
    by = BIN_WIDTH
  )
  list(
    `Equal count, 5 bins` = quantile(x, seq(0, 1, length.out = 6), na.rm = TRUE),
    `Equal count, 10 bins` = quantile(x, seq(0, 1, length.out = 11), na.rm = TRUE),
    `Equal width, 2 degrees` = width_breaks,
    `Equal width, 2 degrees, tails pooled` =
      pool_tails(x, width_breaks, d$GID_1, THIN_UNITS)
  )
}

# ---------------------------------------------------------------- panels ------

nb <- load_notebook_env("test_functions.Rmd")

# Only the base climate block is estimated, so no deviation variables are
# needed beyond what the block itself uses.
NEEDED <- c("dlgrp_pc_usd", "TM", "mean_TM_all", "RR", "mean_RR_all")

datasets <- c(
  `PWT national (GADM0)` = "PWT",
  `DOSE subnational (GADM1)` = "DOSE"
)

panels <- lapply(names(datasets), function(label) {
  cat("Building", label, "...\n")
  panel <- nb$build_dat(econ_data = datasets[[label]]) %>%
    filter(if_all(all_of(NEEDED), is.finite))
  stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
  cat(sprintf(
    "  %d unit-years, %d units, %d countries, %d-%d\n",
    nrow(panel), n_distinct(panel$GID_1), n_distinct(panel$GID_0),
    min(panel$year), max(panel$year)
  ))
  panel
})
names(panels) <- names(datasets)

# ---------------------------------------------------------------- models ------

climate_forms <- function(n_bins) c(
  linear = "TM + TM:mean_TM_all",
  bins = "i(mean_TM_bin, TM)",
  bins_ref = "TM + i(mean_TM_bin, TM, ref = 1)",
  encompassing = sprintf(
    "TM + TM:mean_TM_all + i(mean_TM_bin, TM, ref = %d)", ceiling(n_bins / 2)
  )
)

BHM_PRECIP <- "RR + RR:mean_RR_all"

run <- function(label, scheme) {
  d <- panels[[label]]
  breaks <- bin_schemes(d)[[scheme]]
  # Equal-width breaks can leave empty bins at the tails, so unused levels are
  # dropped and the survivors renumbered to keep the `ref` indices simple.
  d$mean_TM_bin <- droplevels(cut(
    d$mean_TM_all, breaks = breaks, include.lowest = TRUE
  ))
  stopifnot(!anyNA(d$mean_TM_bin))
  d$mean_TM_bin <- factor(as.integer(d$mean_TM_bin))
  n_bins <- nlevels(d$mean_TM_bin)

  fit <- function(terms) {
    model <- fixest::feols(
      as.formula(paste(nb$outcome, "~", terms, "|", nb$fixed_effects)),
      data = d, panel.id = nb$pan_id, cluster = ~GID_1
    )
    stopifnot(!anyNA(coef(model)))
    model
  }

  bin_summary <- d %>%
    group_by(mean_TM_bin) %>%
    summarise(
      obs = n(), units = n_distinct(GID_1), countries = n_distinct(GID_0),
      mean_TM_all_min = min(mean_TM_all), mean_TM_all_max = max(mean_TM_all),
      mean_TM_all_mean = mean(mean_TM_all),
      .groups = "drop"
    )

  forms <- climate_forms(n_bins)
  models <- setNames(
    lapply(forms, function(f) fit(paste(f, "+", BHM_PRECIP))), names(forms)
  )

  # Free slope per bin, with both clustering levels. On the national panel the
  # two coincide because GID_1 is the country.
  m <- models$bins
  ct <- as.data.frame(coeftable(m))
  sel <- grepl(BIN_PATTERN, rownames(ct))
  se_unit <- se(m, cluster = ~GID_1)[sel]
  se_country <- se(m, cluster = ~GID_0)[sel]
  b <- coef(m)[sel]
  df_unit <- n_distinct(d$GID_1) - 1L
  df_country <- n_distinct(d$GID_0) - 1L

  bin_slopes <- data.frame(
    dataset = label, scheme = scheme, n_bins = n_bins,
    mean_TM_bin = sub("^mean_TM_bin::([0-9]+):TM$", "\\1", rownames(ct)[sel]),
    estimate = unname(b), se = unname(se_unit), se_country = unname(se_country),
    p = 2 * pt(-abs(unname(b / se_unit)), df_unit),
    p_country = 2 * pt(-abs(unname(b / se_country)), df_country),
    ci_low = unname(b - 1.96 * se_unit),
    ci_high = unname(b + 1.96 * se_unit)
  ) %>%
    left_join(bin_summary, by = "mean_TM_bin") %>%
    mutate(thin = units < THIN_UNITS)

  climate_grid <- seq(min(d$mean_TM_all), max(d$mean_TM_all), length.out = 400)
  terms <- c("TM", "TM:mean_TM_all")
  bl <- coef(models$linear)[terms]
  Vl <- vcov(models$linear)[terms, terms]
  A <- cbind(1, climate_grid)
  est_l <- as.vector(A %*% bl)
  se_l <- sqrt(pmax(0, rowSums((A %*% Vl) * A)))
  linear_slopes <- data.frame(
    dataset = label, scheme = scheme, n_bins = n_bins,
    mean_TM_all = climate_grid,
    estimate = est_l, ci_low = est_l - 1.96 * se_l, ci_high = est_l + 1.96 * se_l
  )

  wald_row <- function(model, hypothesis, vcov_arg, cluster_label) {
    w <- wald(model, keep = BIN_PATTERN, vcov = vcov_arg, print = FALSE)
    data.frame(dataset = label, scheme = scheme, n_bins = n_bins,
               hypothesis = hypothesis, clustering = cluster_label,
               F_stat = w$stat, df1 = w$df1, p = w$p)
  }
  tests <- bind_rows(
    wald_row(models$bins_ref, "All climate-bin slopes equal", ~GID_1, "GID_1"),
    wald_row(models$bins_ref, "All climate-bin slopes equal", ~GID_0, "GID_0"),
    wald_row(models$encompassing, "Linear-in-climate interaction is adequate",
             ~GID_1, "GID_1"),
    wald_row(models$encompassing, "Linear-in-climate interaction is adequate",
             ~GID_0, "GID_0")
  )

  model_fit <- purrr::imap_dfr(models[c("linear", "bins")], function(mm, id) {
    data.frame(
      dataset = label, scheme = scheme, n_bins = n_bins, climate_form = id,
      n = nobs(mm), k = length(coef(mm)),
      within_r2 = as.numeric(fitstat(mm, "wr2")[[1]]),
      aic = AIC(mm), bic = BIC(mm)
    )
  })

  list(bin_summary = bin_summary %>%
         mutate(dataset = label, scheme = scheme, n_bins = n_bins, .before = 1),
       bin_slopes = bin_slopes, linear_slopes = linear_slopes,
       tests = tests, model_fit = model_fit)
}

SCHEMES <- names(bin_schemes(panels[[1]]))

passes <- unlist(
  lapply(names(panels), function(label) {
    lapply(SCHEMES, function(scheme) run(label, scheme))
  }),
  recursive = FALSE
)
gather <- function(name) bind_rows(lapply(passes, `[[`, name))

bin_summary <- gather("bin_summary")
bin_slopes <- gather("bin_slopes")
linear_slopes <- gather("linear_slopes")
tests <- gather("tests")
model_fit <- gather("model_fit")

# ------------------------------------------------------------------- plots ----

dataset_levels <- names(datasets)

prep <- function(x) {
  mutate(
    x,
    dataset = factor(dataset, levels = dataset_levels),
    scheme = factor(scheme, levels = SCHEMES)
  )
}

# Both clustering levels are shown, since they differ on the subnational panel
# and coincide on the national one.
test_labels <- tests %>%
  select(dataset, scheme, hypothesis, clustering, p) %>%
  tidyr::pivot_wider(names_from = hypothesis, values_from = p) %>%
  group_by(dataset, scheme) %>%
  summarise(
    label = sprintf(
      "linear adequate, p:\n  cluster unit  %.2g\n  cluster country  %.2g",
      .data[["Linear-in-climate interaction is adequate"]][clustering == "GID_1"],
      .data[["Linear-in-climate interaction is adequate"]][clustering == "GID_0"]
    ),
    .groups = "drop"
  )

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_ribbon(
    data = prep(linear_slopes),
    aes(x = mean_TM_all, ymin = ci_low * PP, ymax = ci_high * PP),
    fill = "#0072B2", alpha = 0.15
  ) +
  geom_line(
    data = prep(linear_slopes),
    aes(x = mean_TM_all, y = estimate * PP, color = "Linear in mean_TM_all"),
    linewidth = 0.8
  ) +
  geom_segment(
    data = prep(bin_slopes),
    aes(x = mean_TM_all_min, xend = mean_TM_all_max,
        y = estimate * PP, yend = estimate * PP,
        color = "Free slope per climate bin"),
    linewidth = 0.6
  ) +
  geom_errorbar(
    data = prep(bin_slopes),
    aes(x = mean_TM_all_mean, ymin = ci_low * PP, ymax = ci_high * PP,
        color = "Free slope per climate bin"),
    width = 0.6, linewidth = 0.4
  ) +
  geom_point(
    data = prep(bin_slopes),
    aes(x = mean_TM_all_mean, y = estimate * PP,
        color = "Free slope per climate bin", shape = thin),
    size = 1.6, fill = "white"
  ) +
  geom_text(
    data = prep(test_labels), aes(x = -Inf, y = -Inf, label = label),
    hjust = -0.05, vjust = -0.2, size = 2.6, color = "grey25", lineheight = 1.1
  ) +
  scale_color_manual(
    values = c(`Linear in mean_TM_all` = "#0072B2",
               `Free slope per climate bin` = "#D55E00"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 21),
    labels = c(`FALSE` = paste0(THIN_UNITS, "+ units"),
               `TRUE` = paste0("fewer than ", THIN_UNITS, " units")),
    name = NULL
  ) +
  facet_grid(scheme ~ dataset, scales = "free_y") +
  labs(
    title = "Base climate block only: temperature slope by long-run mean temperature",
    subtitle = paste(
      "No deviation term in any model. Horizontal bars span each bin's range of",
      "mean_TM_all, so narrow bars mark\ncompressed equal-count bins.",
      "Hollow points mark bins with too few units to interpret."
    ),
    x = "Long-run mean temperature, mean_TM_all (degrees Celsius)",
    y = "Effect of +1 degree on growth (percentage points)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.9, "lines"))

ggsave(file.path(out, "climate_slope_bins_national.png"), p,
       width = 12, height = 14, dpi = 180)

# ------------------------------------------------------------------ output ----

write_out(bin_summary, "climate_bin_definitions")
write_out(bin_slopes, "climate_bin_slopes")
write_out(tests, "linearity_tests")
write_out(model_fit, "model_fit_comparison")

cat("\n=== Tests ===\n")
print(as.data.frame(tests), row.names = FALSE, digits = 4)
cat("\n=== Fit: linear vs binned climate interaction ===\n")
print(as.data.frame(model_fit), row.names = FALSE, digits = 5)
cat("\n=== PWT national, equal-width 2-degree bins, tails pooled ===\n")
print(
  bin_slopes %>%
    filter(dataset == "PWT national (GADM0)",
           scheme == "Equal width, 2 degrees, tails pooled") %>%
    transmute(bin = mean_TM_bin,
              range = sprintf("%.1f to %.1f", mean_TM_all_min, mean_TM_all_max),
              countries, obs, slope_pp = round(estimate * PP, 2),
              se_pp = round(se * PP, 2), p = signif(p, 3), thin) %>%
    as.data.frame(),
  row.names = FALSE
)
cat("\nWrote figures and tables to ", out, "\n", sep = "")
