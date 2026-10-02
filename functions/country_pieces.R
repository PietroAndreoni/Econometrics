# Country-by-country "pieces" of a pooled climate regression, and Koppen-Geiger
# climate zones for grouping them.
#
# country_pieces(): the BHM temperature quadratic and the per-country OLS
# slopes that identify it. Growth and TM are residualised on `controls` and
# the specification's fixed effects, so each country keeps only its
# within-country, detrended temperature variation. The country slope
# b_i = sum(x y) / sum(x^2) is its own OLS marginal effect, and by
# Frisch-Waugh-Lovell the pooled linear coefficient is exactly sum_i w_i b_i
# with w_i = sum_t x_it^2 / sum x^2 (checked).
#
# country_deviation_curves(): per-country quadratic or natural-spline fits of
# the partial residual of a deviation term, to set against the pooled
# quadratic in that term. See its comment.

# Broad Koppen-Geiger family per unit, as in main_analysis.Rmd: source code 5
# (polar) is merged into 4. `mean_kg` is constant within GID_1; at GID_0 each
# country takes the family of most of its regions (ties: the lower code).
KG_ZONES <- c("Tropical (A)", "Arid (B)", "Temperate (C)",
              "Continental/Polar (D/E)")

climate_zones <- function(id = c("GID_1", "GID_0"),
                          path = file.path("data", "kg_data.csv")) {
  id <- match.arg(id)
  kg <- utils::read.csv(path, stringsAsFactors = FALSE) %>%
    filter(!is.na(mean_kg)) %>%
    distinct(GID_0, GID_1, mean_kg) %>%
    mutate(mean_kg = pmin(as.integer(mean_kg), 4L))
  stopifnot(!anyDuplicated(kg$GID_1))
  kg %>%
    count(across(all_of(c(id, "mean_kg")))) %>%
    arrange(across(all_of(id)), desc(n), mean_kg) %>%
    distinct(across(all_of(id)), .keep_all = TRUE) %>%
    mutate(climate_zone = factor(KG_ZONES[mean_kg], levels = KG_ZONES)) %>%
    select(all_of(id), climate_zone)
}

# Adds `climate_zone` at `level`: "GID_1" gives each region its own zone,
# "GID_0" gives every region its country's majority zone. The default is GID_1
# for subnational panels and GID_0 for national ones (where GID_1 = GID_0).
# Units without Koppen data are "Unclassified".
add_climate_zone <- function(data, level = NULL) {
  if (is.null(level)) {
    level <- if (all(data$GID_1 == data$GID_0)) "GID_0" else "GID_1"
  }
  data %>%
    left_join(climate_zones(level), by = level) %>%
    mutate(climate_zone = unclassified_zone(climate_zone))
}
unclassified_zone <- function(zone) {
  out <- as.character(zone)
  out[is.na(out)] <- "Unclassified"
  droplevels(factor(out, levels = c(KG_ZONES, "Unclassified")))
}

# `unit` is the grouping of the pieces: GID_1 (one per region, or per country
# on national panels) or GID_0 (one per country, pooling its regions).
# `climate_zone` must be constant within `unit`.
country_pieces <- function(panel, controls, spec = main_spec(),
                           segment_sd = 2, unit = "GID_1", min_years = 10) {
  m_quad <- fit_main(paste("TM + TM^2 +", controls), panel, spec)
  m_lin <- fit_main(paste("TM +", controls), panel, spec)
  stopifnot(nobs(m_lin) == nrow(panel), nobs(m_quad) == nrow(panel))
  b <- coef(m_quad)
  t_opt <- -b[["TM"]] / (2 * b[["I(TM^2)"]])

  t_grid <- seq(floor(min(panel$TM)), ceiling(max(panel$TM)), by = 0.25)
  curve <- tibble::tibble(TM = t_grid) %>%
    bind_cols(linear_combination(
      m_quad, cbind(TM = t_grid - t_opt, `I(TM^2)` = t_grid^2 - t_opt^2)
    ))
  marginal <- tibble::tibble(TM = t_grid) %>%
    bind_cols(linear_combination(
      m_quad, cbind(TM = 1, `I(TM^2)` = 2 * t_grid)
    ))

  residualise <- function(lhs) {
    fit <- fit_main(controls, panel, modifyList(spec, list(outcome = lhs)))
    stopifnot(nobs(fit) == nrow(panel))
    stats::resid(fit)
  }
  resid_panel <- panel %>%
    mutate(unit = .data[[unit]], x = residualise("TM"),
           y = residualise(spec$outcome)) %>%
    group_by(unit) %>%
    mutate(e = y - x * sum(x * y) / sum(x^2)) %>%
    ungroup()
  # SE of each country slope clustered by year, so regions of one country
  # hit by the same shock count once (plain HC on national panels), with the
  # G / (G - 1) finite-cluster factor; gap_p refers gap_z to t(G - 1), G the
  # number of years. Countries with fewer than `min_years` years get no test.
  slope_se <- resid_panel %>%
    group_by(unit, year) %>%
    summarise(score = sum(x * e), sxx_t = sum(x^2), .groups = "drop") %>%
    group_by(unit) %>%
    summarise(years = n(),
              slope_se = sqrt(sum(score^2) * years / (years - 1)) /
                sum(sxx_t),
              .groups = "drop")

  pieces <- resid_panel %>%
    group_by(unit, climate_zone) %>%
    summarise(
      obs = n(),
      regions = n_distinct(GID_1),
      first_year = min(year),
      last_year = max(year),
      TM_mean = mean(TM),
      sxx = sum(x^2),
      slope = sum(x * y) / sxx,
      x_sd = sqrt(sxx / n()),
      .groups = "drop"
    ) %>%
    left_join(slope_se, by = "unit") %>%
    mutate(
      weight = sxx / sum(sxx),
      curve_at_mean = b[["TM"]] * (TM_mean - t_opt) +
        b[["I(TM^2)"]] * (TM_mean^2 - t_opt^2),
      curve_slope = b[["TM"]] + 2 * b[["I(TM^2)"]] * TM_mean,
      # Departure from the curve's tangent in units of the slope's SE.
      gap_z = (slope - curve_slope) / slope_se,
      gap_p = if_else(years >= min_years,
                      2 * stats::pt(-abs(gap_z), years - 1), NA_real_),
      x0 = TM_mean - segment_sd * x_sd,
      x1 = TM_mean + segment_sd * x_sd,
      y0 = curve_at_mean - slope * segment_sd * x_sd,
      y1 = curve_at_mean + slope * segment_sd * x_sd
    )

  pooled_from_pieces <- sum(pieces$weight * pieces$slope)
  stopifnot(isTRUE(all.equal(pooled_from_pieces, coef(m_lin)[["TM"]])))
  list(m_quad = m_quad, m_lin = m_lin, t_opt = t_opt, curve = curve,
       marginal = marginal, pieces = pieces,
       pooled_from_pieces = pooled_from_pieces, segment_sd = segment_sd,
       panel = panel)
}

# The two-panel pieces figure; the `n_label` countries furthest from the
# curve's tangent (largest |gap_z|) are labelled with their ISO3 code.
plot_country_pieces <- function(res, title, spec_label, n_label = 10) {
  pieces <- res$pieces
  labelled <- slice_max(pieces, abs(gap_z), n = n_label)
  weight_scale <- list(
    scale_alpha_continuous(range = c(0.15, 0.95), guide = "none"),
    scale_linewidth_continuous(range = c(0.3, 1.6), guide = "none")
  )

  p_levels <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_segment(data = pieces,
                 aes(x = x0, y = y0, xend = x1, yend = y1, alpha = weight,
                     linewidth = weight),
                 colour = "#C0392B") +
    geom_ribbon(data = res$curve,
                aes(x = TM, ymin = conf_low, ymax = conf_high),
                fill = "grey20", alpha = 0.18) +
    geom_line(data = res$curve, aes(x = TM, y = estimate), linewidth = 1) +
    ggrepel::geom_text_repel(
      data = labelled, aes(x = x1, y = y1, label = unit),
      size = 2.8, colour = "#7B241C", min.segment.length = 0, seed = 1
    ) +
    weight_scale +
    labs(
      x = NULL, y = "Growth relative to optimum", title = title,
      subtitle = sprintf(
        paste0("%s\nBlack: pooled quadratic, 95%% CI, optimum %.1f°C. ",
               "Red: country OLS slopes at mean temperature, ±%g SD of ",
               "residual temperature;\nopacity and width = weight in the ",
               "pooled estimate. Labels: the %d countries furthest from the ",
               "tangent (in SEs). %d countries, %d-%d."),
        spec_label, res$t_opt, res$segment_sd, n_label, nrow(pieces),
        min(res$panel$year), max(res$panel$year)
      )
    ) +
    theme_classic()

  ylim <- stats::quantile(pieces$slope, c(0.02, 0.98))
  p_slopes <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(data = pieces,
               aes(x = TM_mean, y = slope, alpha = weight, size = weight),
               colour = "#C0392B") +
    geom_ribbon(data = res$marginal,
                aes(x = TM, ymin = conf_low, ymax = conf_high),
                fill = "grey20", alpha = 0.18) +
    geom_line(data = res$marginal, aes(x = TM, y = estimate), linewidth = 1) +
    ggrepel::geom_text_repel(
      data = filter(labelled, slope >= ylim[1], slope <= ylim[2]),
      aes(x = TM_mean, y = slope, label = unit),
      size = 2.8, colour = "#7B241C", min.segment.length = 0, seed = 1
    ) +
    scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") +
    scale_size_area(max_size = 5, guide = "none") +
    coord_cartesian(ylim = ylim) +
    labs(
      x = "Temperature (°C)", y = "Marginal effect (dg / dT)",
      caption = sprintf(
        paste0("Weighted mean of country slopes = pooled linear slope = ",
               "%.4f. Lower panel clipped to the 2nd-98th percentile of ",
               "slopes."),
        res$pooled_from_pieces
      )
    ) +
    theme_classic()

  p_levels / p_slopes + plot_layout(heights = c(3, 2))
}

# Per-country response to the deviation `term` (e.g. zTM), set against the
# pooled quadratic b * term^2 of `model`. The partial residual
# r = resid(model) + b * term^2 is growth net of the fixed effects and every
# other regressor. Each country's curve is the pooled b * term^2 plus its own
# intercept and deviation in `term`, fitted to r - b * term^2 and evaluated
# over its [probs] quantile range of `term`. The deviation `basis` is
#   "quadratic": c1 * term + c2 * term^2 (2 df), so the country's own curve is
#                a full quadratic a + c1 * term + (b + c2) * term^2; c1 picks
#                up asymmetry that the pooled pure square cannot. Its
#                coefficients (own_linear, own_quadratic) and year-clustered
#                SEs are returned.
#   "spline":    a natural spline with `df` degrees of freedom.
# `wald_chi2` tests that the deviation is zero, i.e. that the country follows
# the pooled curve, with errors clustered by year so that regions of one
# country hit by the same shock count once (plain heteroskedasticity-robust on
# national panels). With ~50 years per country the robust test over-rejects,
# so the covariance is the leverage-corrected HC3 (cluster jackknife) one and
# p_value refers wald_chi2 / k to F(k, years - 1), k the deviation's degrees
# of freedom. Countries with fewer than `min_years` years get no test. `unit`
# groups as in country_pieces(); with GID_0 on a regional panel one curve is
# fitted to all of a country's region-years (r has mean zero within each
# region, up to b * mean(term^2), so a single intercept per country suffices).
country_deviation_curves <- function(panel, model, term = "zTM",
                                     coef_name = paste0("I(", term, "^2)"),
                                     basis = c("quadratic", "spline"),
                                     df = 3, probs = c(0.02, 0.98),
                                     n_grid = 40, unit = "GID_1",
                                     min_years = 10) {
  basis <- match.arg(basis)
  stopifnot(nobs(model) == nrow(panel))
  b <- coef(model)[[coef_name]]
  d <- panel %>%
    mutate(unit = .data[[unit]], z = .data[[term]],
           gap = stats::resid(model))  # = r - b * z^2
  deviation_formula <- if (basis == "quadratic") {
    gap ~ z + I(z^2)
  } else {
    gap ~ splines::ns(z, df = df)
  }
  k <- if (basis == "quadratic") 2L else as.integer(df)
  slopes <- seq(2, k + 1)

  fits <- lapply(split(d, d$unit), function(u) {
    deviation <- stats::lm(deviation_formula, data = u)
    n_years <- n_distinct(u$year)
    wald <- NA_real_
    se <- rep(NA_real_, k)
    if (n_years >= min_years) {
      V <- sandwich::vcovCL(deviation, cluster = u$year, type = "HC3")
      beta <- stats::coef(deviation)[slopes]
      wald <- as.numeric(t(beta) %*% solve(V[slopes, slopes], beta))
      se <- sqrt(diag(V)[slopes])
    }
    range_z <- stats::quantile(u$z, probs)
    grid <- tibble::tibble(z = seq(range_z[[1]], range_z[[2]],
                                   length.out = n_grid))
    grid$pooled <- mean(u$gap) + b * grid$z^2
    grid$estimate <- stats::predict(deviation, newdata = grid) + b * grid$z^2
    stats_row <- tibble::tibble(
      unit = u$unit[1], climate_zone = u$climate_zone[1],
      obs = nrow(u), regions = n_distinct(u$GID_1),
      years = n_years, z_sd = stats::sd(u$z),
      wald_chi2 = wald,
      p_value = stats::pf(wald / k, k, n_years - 1, lower.tail = FALSE)
    )
    if (basis == "quadratic") {
      c_dev <- stats::coef(deviation)
      stats_row <- mutate(
        stats_row,
        own_linear = c_dev[["z"]], own_linear_se = se[1],
        own_quadratic = b + c_dev[["I(z^2)"]], own_quadratic_se = se[2]
      )
    }
    list(
      curve = mutate(grid, unit = u$unit[1],
                     climate_zone = u$climate_zone[1]),
      stats = stats_row
    )
  })
  list(
    curves = bind_rows(lapply(fits, `[[`, "curve")),
    stats = bind_rows(lapply(fits, `[[`, "stats")),
    b = b, basis = basis, df = k, probs = probs
  )
}
