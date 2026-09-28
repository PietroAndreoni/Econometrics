# Replace the linear BHM climate interaction with a nonparametric one.
#
# The notebook's base climate block is
#
#   TM + TM:mean_TM_all + RR + RR:mean_RR_all
#
# Because `mean_TM_all` is time-invariant per region and the fixed effects
# `year + GID_1[year]` include region intercepts and region-specific trends,
# `mean_TM_all` enters only through its interaction. The block therefore does
# not estimate a response in temperature levels: it estimates the marginal
# effect of a within-region temperature deviation, restricted to vary linearly
# with the region's long-run mean temperature.
#
# This script relaxes that restriction by estimating a separate temperature
# slope within each bin of `mean_TM_all`, and tests the linear restriction
# against it. Binning `TM` itself would not work here: annual mean temperature
# rarely crosses a fixed level-bin boundary within a region, so region fixed
# effects would absorb the level bins almost entirely.
#
# Because the within-region variation in TM *is* the anomaly, the estimated
# climate profile depends on the deviation term that is also in the model. Four
# are reported, including the fully nonparametric signed-deviation bins, which
# is the benchmark that leaves the deviation response unrestricted. Bin counts
# of 5 and 10 are both reported so that features of the profile can be checked
# against the resolution used to find them.
#
# Prerequisite: the notebook's parquet inputs under data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
})

source("rmd_chunks.R")

out <- "results/bhm_climate_slope_bins"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

N_BINS_SET <- c(5L, 10L)   # Bin counts for long-run mean temperature.
PP <- 100                  # Slopes in log-growth percentage points per degree.

# ---------------------------------------------------------------- panel -------

nb <- load_notebook_env("test_functions.Rmd")

cat("Building the DOSE panel...\n")
panel <- nb$build_dat(econ_data = "DOSE") %>%
  mutate(
    hTM_10 = pmax(0, abs(zTM) - 1.0),   # Best-fitting hinge.
    hTM_15 = pmax(0, abs(zTM) - 1.5)    # The notebook's `zzTM`.
  )

# One sample for every specification, including the bin variables, so that the
# comparison is about functional form rather than coverage.
needed <- c(
  "dlgrp_pc_usd", "TM", "mean_TM_all", "RR", "mean_RR_all",
  "zTM", "hTM_10", "hTM_15", "abs_zRRp", "TM_bin_signed", "RR_bin_signed"
)
est <- panel %>% filter(if_all(all_of(needed), is.finite))
stopifnot(!anyDuplicated(est[c("GID_1", "year")]))

cat(sprintf(
  "Estimation sample: %d region-years, %d regions, %d countries, %d-%d\n",
  nrow(est), n_distinct(est$GID_1), n_distinct(est$GID_0),
  min(est$year), max(est$year)
))

# ---------------------------------------------------------------- models ------

fit <- function(terms) {
  model <- fixest::feols(
    as.formula(paste(nb$outcome, "~", terms, "|", nb$fixed_effects)),
    data = est, panel.id = nb$pan_id, cluster = ~GID_1
  )
  stopifnot(!anyNA(coef(model)))
  model
}

# The BHM precipitation block. Deviation terms live in `deviation_specs` so that
# each benchmark controls both its temperature and its precipitation deviation.
BHM_PRECIP <- "RR + RR:mean_RR_all"

deviation_specs <- c(
  "No deviation term" = "",
  "Linear hinge at 1.0 SD" = "hTM_10 + abs_zRRp",
  "Cubic hinge at 1.5 SD (notebook)" = "hTM_15^3 + abs_zRRp",
  "Signed deviation bins" =
    "i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)"
)

BIN_PATTERN <- "^mean_TM_bin::[0-9]+:TM$"

climate_forms <- function(n_bins) {
  c(
    linear = "TM + TM:mean_TM_all",
    bins = "i(mean_TM_bin, TM)",
    # Reference-coded, so a joint test on the interactions is a test that every
    # climate-bin slope is equal.
    bins_ref = "TM + i(mean_TM_bin, TM, ref = 1)",
    # Both readings at once, so a joint test on the interactions is a test of
    # the linear restriction against the free slopes.
    encompassing = sprintf(
      "TM + TM:mean_TM_all + i(mean_TM_bin, TM, ref = %d)",
      ceiling(n_bins / 2)
    )
  )
}

# Everything below is computed per bin count. `est$mean_TM_bin` is rebuilt at
# the start of each pass, and every model for that pass is fitted before the
# next one overwrites it.
run_bin_count <- function(n_bins) {
  breaks <- quantile(
    est$mean_TM_all, seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE
  )
  est$mean_TM_bin <<- cut(
    est$mean_TM_all, breaks = breaks,
    include.lowest = TRUE, labels = as.character(seq_len(n_bins))
  )
  stopifnot(!anyNA(est$mean_TM_bin))

  bin_summary <- est %>%
    group_by(mean_TM_bin) %>%
    summarise(
      obs = n(), regions = n_distinct(GID_1), countries = n_distinct(GID_0),
      mean_TM_all_min = min(mean_TM_all), mean_TM_all_max = max(mean_TM_all),
      mean_TM_all_mean = mean(mean_TM_all),
      .groups = "drop"
    ) %>%
    mutate(n_bins = n_bins, .before = 1)

  forms <- climate_forms(n_bins)
  grid <- tidyr::expand_grid(
    deviation = names(deviation_specs), climate = names(forms)
  ) %>%
    mutate(
      id = paste(climate, deviation, sep = " | "),
      terms = purrr::map2_chr(climate, deviation, function(climate, deviation) {
        paste(c(forms[[climate]], BHM_PRECIP,
                deviation_specs[[deviation]][nzchar(deviation_specs[[deviation]])]),
              collapse = " + ")
      })
    )
  models <- setNames(lapply(grid$terms, fit), grid$id)

  # Free slope in each climate bin, from the unreferenced parameterization.
  bin_slopes <- purrr::map_dfr(names(deviation_specs), function(deviation) {
    model <- models[[paste("bins", deviation, sep = " | ")]]
    ct <- as.data.frame(coeftable(model))
    ci <- confint(model)
    sel <- grepl(BIN_PATTERN, rownames(ct))
    data.frame(
      n_bins = n_bins, deviation = deviation,
      mean_TM_bin = sub("^mean_TM_bin::([0-9]+):TM$", "\\1", rownames(ct)[sel]),
      estimate = ct[sel, 1], se = ct[sel, 2], p = ct[sel, 4],
      ci_low = ci[sel, 1], ci_high = ci[sel, 2]
    )
  }) %>%
    left_join(bin_summary %>% select(-n_bins), by = "mean_TM_bin")

  # Implied slope under the linear restriction, with a delta-method band.
  climate_grid <- seq(min(est$mean_TM_all), max(est$mean_TM_all), length.out = 400)
  linear_slopes <- purrr::map_dfr(names(deviation_specs), function(deviation) {
    model <- models[[paste("linear", deviation, sep = " | ")]]
    terms <- c("TM", "TM:mean_TM_all")
    b <- coef(model)[terms]
    V <- vcov(model)[terms, terms]
    A <- cbind(1, climate_grid)
    estimate <- as.vector(A %*% b)
    se <- sqrt(pmax(0, rowSums((A %*% V) * A)))
    data.frame(
      n_bins = n_bins, deviation = deviation, mean_TM_all = climate_grid,
      estimate = estimate, ci_low = estimate - 1.96 * se,
      ci_high = estimate + 1.96 * se
    )
  })

  wald_row <- function(model, label, deviation) {
    w <- wald(model, keep = BIN_PATTERN, print = FALSE)
    data.frame(n_bins = n_bins, deviation = deviation, hypothesis = label,
               F_stat = w$stat, df1 = w$df1, df2 = w$df2, p = w$p)
  }
  tests <- purrr::map_dfr(names(deviation_specs), function(deviation) {
    bind_rows(
      wald_row(models[[paste("bins_ref", deviation, sep = " | ")]],
               "All climate-bin slopes equal", deviation),
      wald_row(models[[paste("encompassing", deviation, sep = " | ")]],
               "Linear-in-climate interaction is adequate", deviation)
    )
  })

  model_fit <- purrr::imap_dfr(models, function(m, id) {
    parts <- strsplit(id, " | ", fixed = TRUE)[[1]]
    data.frame(
      n_bins = n_bins, climate_form = parts[[1]], deviation = parts[[2]],
      n = nobs(m), k = length(coef(m)),
      within_r2 = as.numeric(fitstat(m, "wr2")[[1]]),
      aic = AIC(m), bic = BIC(m)
    )
  }) %>%
    filter(climate_form %in% c("linear", "bins"))

  coefs <- purrr::imap_dfr(models, function(m, id) {
    ct <- as.data.frame(coeftable(m))
    data.frame(n_bins = n_bins, model = id, term = rownames(ct),
               estimate = ct[, 1], se = ct[, 2], p = ct[, 4])
  })

  list(bin_summary = bin_summary, bin_slopes = bin_slopes,
       linear_slopes = linear_slopes, tests = tests,
       model_fit = model_fit, coefs = coefs)
}

passes <- lapply(N_BINS_SET, run_bin_count)
gather <- function(name) bind_rows(lapply(passes, `[[`, name))

bin_summary <- gather("bin_summary")
bin_slopes <- gather("bin_slopes")
linear_slopes <- gather("linear_slopes")
tests <- gather("tests")
model_fit <- gather("model_fit") %>% arrange(n_bins, deviation, desc(within_r2))
coef_table <- gather("coefs")

# ------------------------------------------------------------------- plots ----

deviation_levels <- names(deviation_specs)
bin_labels <- sprintf("%d climate bins", N_BINS_SET)

prep <- function(x) {
  mutate(
    x,
    deviation = factor(deviation, levels = deviation_levels),
    bin_label = factor(sprintf("%d climate bins", n_bins), levels = bin_labels)
  )
}

test_labels <- tests %>%
  select(n_bins, deviation, hypothesis, p) %>%
  tidyr::pivot_wider(names_from = hypothesis, values_from = p) %>%
  transmute(
    n_bins, deviation,
    label = sprintf(
      "equal slopes: p = %.2g\nlinear adequate: p = %.2g",
      .data[["All climate-bin slopes equal"]],
      .data[["Linear-in-climate interaction is adequate"]]
    )
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
  # Each bin is drawn across the range of long-run mean temperature it covers.
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
        color = "Free slope per climate bin"),
    size = 1.5
  ) +
  geom_text(
    data = prep(test_labels), aes(x = -Inf, y = -Inf, label = label),
    hjust = -0.06, vjust = -0.3, size = 2.7, color = "grey25", lineheight = 1.1
  ) +
  scale_color_manual(
    values = c(`Linear in mean_TM_all` = "#0072B2",
               `Free slope per climate bin` = "#D55E00"),
    name = NULL
  ) +
  facet_grid(bin_label ~ deviation) +
  labs(
    title = "Marginal effect of a temperature deviation, by long-run mean temperature",
    subtitle = paste(
      "The BHM block restricts this slope to be linear in mean_TM_all; the bins",
      "let it vary freely. Each column holds a different deviation term,",
      "\nincluding the nonparametric signed bins. Horizontal bars span each",
      "bin's range of mean_TM_all. 95% intervals clustered by GID_1."
    ),
    x = "Long-run mean temperature, mean_TM_all (degrees Celsius)",
    y = "Effect of +1 degree on growth (percentage points)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.9, "lines"))

ggsave(file.path(out, "climate_slope_bins.png"), p,
       width = 15, height = 8.5, dpi = 180)

# ------------------------------------------------------------------ output ----

write_out(bin_summary, "climate_bin_definitions")
write_out(bin_slopes, "climate_bin_slopes")
write_out(tests, "linearity_tests")
write_out(model_fit, "model_fit_comparison")
write_out(coef_table, "coefficients")

cat("\n=== Climate bin definitions ===\n")
print(as.data.frame(bin_summary), row.names = FALSE, digits = 4)
cat("\n=== Tests ===\n")
print(as.data.frame(tests), row.names = FALSE, digits = 4)
cat("\n=== Fit: linear vs binned climate interaction ===\n")
print(as.data.frame(model_fit), row.names = FALSE, digits = 5)
cat("\n=== Free temperature slope by climate bin, 10 bins (pp per degree) ===\n")
print(
  bin_slopes %>%
    filter(n_bins == 10) %>%
    transmute(deviation, bin = mean_TM_bin,
              range = sprintf("%.1f to %.1f", mean_TM_all_min, mean_TM_all_max),
              regions, slope_pp = round(estimate * PP, 3),
              se_pp = round(se * PP, 3), p = signif(p, 3)) %>%
    as.data.frame(),
  row.names = FALSE
)
cat("\nWrote figures and tables to ", out, "\n", sep = "")
