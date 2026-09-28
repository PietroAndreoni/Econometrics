# Relax the linear BHM climate interaction.
#
# In the base block  TM + TM:mean_TM_all + RR + RR:mean_RR_all  `mean_TM_all` is
# time-invariant, so with unit fixed effects and trends it enters only through
# its interaction: the block estimates the marginal effect of a within-unit
# temperature deviation, restricted to vary linearly with the unit's long-run
# mean temperature. These functions replace that restriction with a free slope
# per bin of `mean_TM_all` and test it. (Binning TM itself would not work:
# annual temperature rarely crosses a fixed level boundary within a unit, so
# the fixed effects would absorb level bins almost entirely.)

CLIMATE_BIN_PATTERN <- "^mean_TM_bin::[0-9]+:TM$"

# The four parameterisations of the temperature slope, for `n_bins` climate bins.
#   linear        the BHM restriction;
#   bins          a free slope per bin (unreferenced: each coefficient is a slope);
#   bins_ref      reference-coded, so a joint test on the interactions tests that
#                 every bin slope is equal;
#   encompassing  both at once, so a joint test on the interactions tests the
#                 linear restriction against the free slopes.
climate_bin_forms <- function(n_bins) {
  c(
    linear = "TM + TM:mean_TM_all",
    bins = "i(mean_TM_bin, TM)",
    bins_ref = "TM + i(mean_TM_bin, TM, ref = 1)",
    encompassing = sprintf(
      "TM + TM:mean_TM_all + i(mean_TM_bin, TM, ref = %d)",
      ceiling(n_bins / 2)
    )
  )
}

# Adds `mean_TM_bin`: bins of mean_TM_all by `breaks`, with empty bins dropped
# and the occupied ones renumbered 1..k so reference indices stay simple.
assign_climate_bins <- function(d, breaks) {
  bins <- droplevels(cut(d$mean_TM_all, breaks = breaks, include.lowest = TRUE))
  if (anyNA(bins)) stop("Some mean_TM_all values fall outside the bin breaks.")
  d$mean_TM_bin <- factor(as.integer(bins))
  d
}

# Estimates the four forms on `d` (after assign_climate_bins()), each with the
# precipitation block and an optional deviation term, and returns:
#   bin_summary    range and coverage of each climate bin;
#   bin_slopes     free slope per bin, clustered by unit (GID_1) and country;
#   linear_slopes  slope implied by the linear restriction over mean_TM_all;
#   tests          equal-slopes and linear-adequacy Wald tests, both clusterings;
#   model_fit      fit of the linear and binned forms;
#   coefficients   every coefficient of every model.
# `labels` (a named list) is prepended to every table to identify the run.
climate_slope_analysis <- function(
    d,
    deviation_terms = "",
    precip_terms = "RR + RR:mean_RR_all",
    spec = main_spec(),
    thin_units = 10L,
    labels = list()
) {
  n_bins <- nlevels(d$mean_TM_bin)
  forms <- climate_bin_forms(n_bins)
  models <- lapply(forms, function(form) {
    terms <- paste(c(form, precip_terms, deviation_terms[nzchar(deviation_terms)]),
                   collapse = " + ")
    model <- fit_main(terms, d, spec)
    if (anyNA(stats::coef(model))) stop("Collinear terms in: ", terms)
    model
  })
  tag <- function(x) {
    bind_cols(tibble::as_tibble(labels), tibble::tibble(n_bins = n_bins), x)
  }

  bin_summary <- d %>%
    group_by(mean_TM_bin) %>%
    summarise(
      obs = n(), units = n_distinct(GID_1), countries = n_distinct(GID_0),
      mean_TM_all_min = min(mean_TM_all), mean_TM_all_max = max(mean_TM_all),
      mean_TM_all_mean = mean(mean_TM_all),
      .groups = "drop"
    ) %>%
    mutate(mean_TM_bin = as.character(mean_TM_bin))

  slope_table <- function(vcov, suffix) {
    ct <- as.data.frame(fixest::coeftable(models$bins, vcov = vcov))
    ci <- stats::confint(models$bins, vcov = vcov)
    sel <- grepl(CLIMATE_BIN_PATTERN, rownames(ct))
    out <- tibble::tibble(
      mean_TM_bin = sub("^mean_TM_bin::([0-9]+):TM$", "\\1", rownames(ct)[sel]),
      estimate = ct[sel, 1], se = ct[sel, 2], p = ct[sel, 4],
      ci_low = ci[sel, 1], ci_high = ci[sel, 2]
    )
    if (nzchar(suffix)) {
      out <- out %>%
        select(mean_TM_bin, se, p) %>%
        rename_with(~ paste0(.x, suffix), c(se, p))
    }
    out
  }
  bin_slopes <- slope_table(NULL, "") %>%
    left_join(slope_table(~GID_0, "_country"), by = "mean_TM_bin") %>%
    left_join(bin_summary, by = "mean_TM_bin") %>%
    mutate(thin = units < thin_units) %>%
    tag()

  climate_grid <- seq(min(d$mean_TM_all), max(d$mean_TM_all), length.out = 400)
  linear_slopes <- bind_cols(
    tibble::tibble(mean_TM_all = climate_grid),
    linear_combination(
      models$linear,
      cbind(TM = 1, `TM:mean_TM_all` = climate_grid)
    )
  ) %>%
    rename(se = std_error, ci_low = conf_low, ci_high = conf_high) %>%
    tag()

  tests <- bind_rows(
    wald_row(models$bins_ref, CLIMATE_BIN_PATTERN,
             label = "All climate-bin slopes equal") %>%
      mutate(clustering = "GID_1"),
    wald_row(models$bins_ref, CLIMATE_BIN_PATTERN, vcov = ~GID_0,
             label = "All climate-bin slopes equal") %>%
      mutate(clustering = "GID_0"),
    wald_row(models$encompassing, CLIMATE_BIN_PATTERN,
             label = "Linear-in-climate interaction is adequate") %>%
      mutate(clustering = "GID_1"),
    wald_row(models$encompassing, CLIMATE_BIN_PATTERN, vcov = ~GID_0,
             label = "Linear-in-climate interaction is adequate") %>%
      mutate(clustering = "GID_0")
  ) %>%
    select(-keep) %>%
    tag()

  model_fit <- bind_rows(lapply(c("linear", "bins"), function(form) {
    m <- models[[form]]
    tibble::tibble(
      climate_form = form, n = stats::nobs(m), k = length(stats::coef(m)),
      within_r2 = as.numeric(fixest::fitstat(m, "wr2")[[1]]),
      aic = stats::AIC(m), bic = stats::BIC(m)
    )
  })) %>%
    tag()

  coefficients <- bind_rows(lapply(names(models), function(form) {
    tidy_fixest(models[[form]], label = form)
  })) %>%
    tag()

  list(
    bin_summary = tag(bin_summary),
    bin_slopes = bin_slopes,
    linear_slopes = linear_slopes,
    tests = tests,
    model_fit = model_fit,
    coefficients = coefficients
  )
}

# Free bin slopes (orange: a bar across each bin's range of mean_TM_all, a
# point and 95% interval at its mean) against the linear restriction (blue
# line and band), in percentage points per degree. `test_labels` has the facet
# columns plus `label`; hollow points mark thin bins.
plot_climate_slope_bins <- function(bin_slopes, linear_slopes, test_labels,
                                    facets, title, subtitle, thin_units = 10L,
                                    scales = "fixed") {
  pp <- 100
  colours <- c(`Linear in mean_TM_all` = GRID_PALETTE[1],
               `Free slope per climate bin` = GRID_PALETTE[2])
  ggplot() +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_ribbon(
      data = linear_slopes,
      aes(x = mean_TM_all, ymin = ci_low * pp, ymax = ci_high * pp),
      fill = GRID_PALETTE[1], alpha = 0.15
    ) +
    geom_line(
      data = linear_slopes,
      aes(x = mean_TM_all, y = estimate * pp, colour = "Linear in mean_TM_all"),
      linewidth = 0.8
    ) +
    geom_segment(
      data = bin_slopes,
      aes(x = mean_TM_all_min, xend = mean_TM_all_max, y = estimate * pp,
          yend = estimate * pp, colour = "Free slope per climate bin"),
      linewidth = 0.6
    ) +
    geom_errorbar(
      data = bin_slopes,
      aes(x = mean_TM_all_mean, ymin = ci_low * pp, ymax = ci_high * pp,
          colour = "Free slope per climate bin"),
      width = 0.6, linewidth = 0.4
    ) +
    geom_point(
      data = bin_slopes,
      aes(x = mean_TM_all_mean, y = estimate * pp,
          colour = "Free slope per climate bin", shape = thin),
      size = 1.8, fill = "white"
    ) +
    geom_text(
      data = test_labels, aes(x = -Inf, y = -Inf, label = label),
      hjust = -0.05, vjust = -0.25, size = 2.6, colour = GRID_INK,
      lineheight = 1.05
    ) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 21),
      labels = c(`FALSE` = paste0(thin_units, "+ units"),
                 `TRUE` = paste0("fewer than ", thin_units, " units")),
      name = NULL, drop = FALSE
    ) +
    facet_grid(facets, scales = scales) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Long-run mean temperature, mean_TM_all (degrees Celsius)",
      y = "Effect of +1 degree on growth (pp)"
    ) +
    theme_results() +
    theme(panel.grid.major.x = element_line(colour = "#f0efea", linewidth = 0.3),
          panel.spacing = unit(0.9, "lines"))
}

# "equal slopes / linear adequate" p-values per facet, for plot labels.
climate_slope_test_labels <- function(tests, facet_cols) {
  tests %>%
    select(all_of(facet_cols), test, clustering, p_value) %>%
    tidyr::pivot_wider(names_from = c(test, clustering), values_from = p_value) %>%
    transmute(
      across(all_of(facet_cols)),
      label = sprintf(
        "equal slopes: p = %.2g (unit) / %.2g (country)\nlinear adequate: p = %.2g (unit) / %.2g (country)",
        .data[["All climate-bin slopes equal_GID_1"]],
        .data[["All climate-bin slopes equal_GID_0"]],
        .data[["Linear-in-climate interaction is adequate_GID_1"]],
        .data[["Linear-in-climate interaction is adequate_GID_0"]]
      )
    )
}
