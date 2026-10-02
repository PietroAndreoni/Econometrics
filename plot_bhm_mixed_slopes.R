# The Burke-Hsiang-Miguel (2015) quadratic with country-specific random slopes,
# on the national PWT panel, fitted with brms. Target model:
#
#   growth ~ 0 + factor(country) + factor(year) + temp + temp2 +
#            (0 + temp + temp2 | country)
#
# i.e. country and year fixed effects, a population quadratic, and correlated
# country deviations from both of its coefficients.
#
# The 178 country dummies make HMC hit its maximum tree depth, so the country
# fixed effects are replaced by a country random intercept plus the exact
# unbalanced-panel Mundlak terms: the country means of every time-varying
# regressor, i.e. of x and x2 and of the year dummies. The latter reduce to a
# full-rank basis of the coverage patterns (mundlak_year_basis()). With these
# terms the random-intercept estimator of the x, x2 coefficients equals the
# two-way within estimator exactly whatever the intercept variance, so the
# intercepts behave like fixed effects (checked against fixest below, without
# random slopes, via lme4). A plain random intercept does not: it lets
# cross-country variation back in. The intercept is correlated with the slopes,
# as unrestricted fixed-effect intercepts would be.
#
# brms defaults: flat priors on the fixed coefficients, half-Student-t on the
# SDs, LKJ(1) on the correlations. Temperature enters as
# x = (TM - T_CENTRE) / T_SCALE for sampler efficiency; with an unstructured
# slope covariance this is a reparametrisation, not a different model, and all
# reported curves are in degrees Celsius. The fit is cached in
# OUT/brms_fit.rds and refitted only when the model or data change.
#
# Compared with:
#   FE: the same fixed part without random slopes (country + year FE), fixest,
#       clustered by country.
# Heterogeneity: posterior of the implied cross-country SD of dg/dT at 5, 15
# and 25 degrees.
#
# Figure: (top) population curves relative to their optimum, with each
# country's own curve (population + posterior-mean country deviation) over its
# observed temperature range, anchored on the population curve at its mean
# temperature; (bottom) each country's marginal effect at its mean
# temperature, raw country OLS slope (residualised on country and year FE, as
# in plot_bhm_country_pieces.R) versus its shrunk mixed-model slope.
#
# Sample: plot_bhm_country_pieces.R, i.e. PWT from its first year.
# Output: results/bhm_mixed_slopes/.

source("load_functions.R")

SPEC <- main_spec(econ_year_min = 1900L)
T_CENTRE <- 15
T_SCALE <- 10
CHAINS <- 4L
ITER <- 1000L         # 500 warmup; ~75 min, the sampler hits max treedepth often.
OUT <- file.path("results", "bhm_mixed_slopes")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# Country means of the year dummies, as a full-rank basis orthogonal to the
# intercept; one column per coverage pattern beyond the first.
mundlak_year_basis <- function(year, id) {
  means <- apply(stats::model.matrix(~ 0 + factor(year)), 2,
                 function(d) stats::ave(d, id))
  q <- qr(cbind(1, means))
  basis <- qr.Q(q)[, seq(2, q$rank), drop = FALSE]
  colnames(basis) <- paste0("coverage_", seq_len(ncol(basis)))
  basis
}

panel <- build_main_dat("PWT", SPEC) %>%
  filter(if_all(all_of(c(SPEC$outcome, "TM")), is.finite)) %>%
  mutate(x = (TM - T_CENTRE) / T_SCALE, x2 = x^2) %>%
  group_by(GID_1) %>%
  mutate(x_country = mean(x), x2_country = mean(x2)) %>%
  ungroup()
coverage <- mundlak_year_basis(panel$year, panel$GID_1)
panel <- bind_cols(panel, as.data.frame(coverage))
stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
cat(sprintf("PWT: %d country-years, %d countries, %d-%d, %d coverage terms\n",
            nrow(panel), n_distinct(panel$GID_1), min(panel$year),
            max(panel$year), ncol(coverage)))

# Estimation ----------------------------------------------------------------------

mundlak_terms <- paste(c("x_country", "x2_country", colnames(coverage)),
                       collapse = " + ")
fixed_part <- paste(SPEC$outcome, "~ factor(year) + x + x2 +", mundlak_terms)

# The Mundlak terms make the random-intercept model reproduce the two-way FE
# estimates exactly; verify it without random slopes before sampling.
m_fe <- fit_panel_model("x + x2", panel, outcome = SPEC$outcome,
                        fixed_effects = "GID_1 + year",
                        panel_id = SPEC$panel_id, cluster = SPEC$cluster)
m_mundlak_check <- lme4::lmer(
  stats::as.formula(paste(fixed_part, "+ (1 | GID_1)")), data = panel,
  control = lme4::lmerControl(optimizer = "bobyqa")
)
stopifnot(isTRUE(all.equal(lme4::fixef(m_mundlak_check)[c("x", "x2")],
                           coef(m_fe)[c("x", "x2")], tolerance = 1e-6)))

m_mixed <- brms::brm(
  brms::bf(
    stats::as.formula(paste(fixed_part, "+ (1 + x + x2 | GID_1)")),
    decomp = "QR"
  ),
  data = panel,
  chains = CHAINS, cores = CHAINS, iter = ITER, seed = 20260929,
  file = file.path(OUT, "brms_fit"), file_refit = "on_change",
  refresh = 100
)
stopifnot(nobs(m_mixed) == nrow(panel), nobs(m_fe) == nrow(panel))

draws <- posterior::as_draws_df(m_mixed) %>%
  as.data.frame() %>%
  select(b_x, b_x2, sd_x = sd_GID_1__x, sd_x2 = sd_GID_1__x2,
         cor = cor_GID_1__x__x2, sigma)
convergence <- brms::rhat(m_mixed)[c("b_x", "b_x2", "sd_GID_1__x",
                                     "sd_GID_1__x2", "cor_GID_1__x__x2")]
beta <- c(x = mean(draws$b_x), x2 = mean(draws$b_x2))

# dg/dT per degree for coefficients on x and x2.
slope_mean_at <- function(b, t) {
  (b[["x"]] + 2 * b[["x2"]] * (t - T_CENTRE) / T_SCALE) / T_SCALE
}
# Implied cross-country SD of dg/dT at temperature t, per draw.
slope_sd_at <- function(t) {
  g1 <- 1 / T_SCALE
  g2 <- 2 * (t - T_CENTRE) / T_SCALE^2
  sqrt(g1^2 * draws$sd_x^2 + g2^2 * draws$sd_x2^2 +
         2 * g1 * g2 * draws$cor * draws$sd_x * draws$sd_x2)
}
summarise_draws <- function(v) {
  c(estimate = mean(v), conf_low = unname(stats::quantile(v, 0.025)),
    conf_high = unname(stats::quantile(v, 0.975)))
}
heterogeneity <- bind_rows(lapply(c(5, 15, 25), function(t) {
  bind_rows(
    c(TM = t, quantity = "population slope",
      summarise_draws(slope_mean_at(list(x = draws$b_x, x2 = draws$b_x2), t))),
    c(TM = t, quantity = "cross-country SD of slope",
      summarise_draws(slope_sd_at(t)))
  )
})) %>%
  mutate(across(-quantity, as.numeric))

# Curves ------------------------------------------------------------------------------

t_grid <- seq(floor(min(panel$TM)), ceiling(max(panel$TM)), by = 0.25)
x_grid <- (t_grid - T_CENTRE) / T_SCALE
optimum_x <- function(b) -b[["x"]] / (2 * b[["x2"]])
x_opt <- optimum_x(beta)

# Mixed population curve: pointwise posterior mean and 95% interval, relative
# to growth at the posterior-mean optimum.
mixed_curves <- bind_rows(
  tibble::tibble(TM = t_grid, quantity = "Growth relative to optimum") %>%
    bind_cols(t(vapply(x_grid, function(xx) summarise_draws(
      draws$b_x * (xx - x_opt) + draws$b_x2 * (xx^2 - x_opt^2)
    ), numeric(3)))),
  tibble::tibble(TM = t_grid, quantity = "Marginal effect (dg / dT)") %>%
    bind_cols(t(vapply(t_grid, function(t) summarise_draws(
      slope_mean_at(list(x = draws$b_x, x2 = draws$b_x2), t)
    ), numeric(3))))
) %>%
  mutate(model = sprintf("Mixed: population curve (optimum %.1f°C)",
                         T_CENTRE + T_SCALE * x_opt))

b_fe <- coef(m_fe)[c("x", "x2")]
x_opt_fe <- optimum_x(b_fe)
fe_curves <- bind_rows(
  tibble::tibble(TM = t_grid, quantity = "Growth relative to optimum") %>%
    bind_cols(linear_combination(
      m_fe, cbind(x = x_grid - x_opt_fe, x2 = x_grid^2 - x_opt_fe^2)
    )),
  tibble::tibble(TM = t_grid, quantity = "Marginal effect (dg / dT)") %>%
    bind_cols(linear_combination(m_fe, cbind(x = 1, x2 = 2 * x_grid) / T_SCALE))
) %>%
  select(-std_error) %>%
  mutate(model = sprintf("Fixed effects, no random slopes (optimum %.1f°C)",
                         T_CENTRE + T_SCALE * x_opt_fe))

curves <- bind_rows(mixed_curves, fe_curves) %>%
  mutate(model = factor(model, levels = unique(model)))

# Country curves and pieces -------------------------------------------------------

residualise <- function(lhs) {
  stats::resid(fit_panel_model("1", panel, outcome = lhs,
                               fixed_effects = "GID_1 + year",
                               panel_id = SPEC$panel_id))
}
country_dev <- brms::ranef(m_mixed)$GID_1[, "Estimate", c("x", "x2")]
country_dev <- tibble::tibble(GID_1 = rownames(country_dev),
                              u_x = country_dev[, "x"],
                              u_x2 = country_dev[, "x2"])
pop_level <- function(xx) {
  beta[["x"]] * (xx - x_opt) + beta[["x2"]] * (xx^2 - x_opt^2)
}

countries <- panel %>%
  mutate(tx = residualise("TM"), ty = residualise(SPEC$outcome)) %>%
  group_by(GID_1) %>%
  summarise(
    years = n(), TM_mean = mean(TM), TM_min = min(TM), TM_max = max(TM),
    ols_slope = sum(tx * ty) / sum(tx^2), sxx = sum(tx^2),
    .groups = "drop"
  ) %>%
  left_join(country_dev, by = "GID_1") %>%
  mutate(
    weight = sxx / sum(sxx),
    x_mean = (TM_mean - T_CENTRE) / T_SCALE,
    mixed_slope = slope_mean_at(list(x = beta[["x"]] + u_x,
                                     x2 = beta[["x2"]] + u_x2), TM_mean),
    population_slope = slope_mean_at(beta, TM_mean)
  )
stopifnot(!anyNA(countries$u_x))

country_curves <- countries %>%
  select(GID_1, TM_min, TM_max, x_mean, u_x, u_x2) %>%
  rowwise() %>%
  reframe(GID_1, TM = seq(TM_min, TM_max, length.out = 20),
          x_mean, u_x, u_x2) %>%
  mutate(
    xx = (TM - T_CENTRE) / T_SCALE,
    estimate = pop_level(x_mean) +
      (beta[["x"]] + u_x) * (xx - x_mean) +
      (beta[["x2"]] + u_x2) * (xx^2 - x_mean^2)
  )

write_out(countries, "country_slopes")
write_out(curves, "curves")
write_out(heterogeneity, "slope_heterogeneity")

# Figure --------------------------------------------------------------------------

colours <- c("#2C6FB7", "#C0392B")
names(colours) <- levels(curves$model)
curve_layers <- function(quantity_name) {
  d <- filter(curves, quantity == quantity_name)
  list(
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3),
    geom_ribbon(data = d, aes(TM, ymin = conf_low, ymax = conf_high,
                              fill = model), alpha = 0.15),
    geom_line(data = d, aes(TM, estimate, colour = model), linewidth = 1),
    scale_colour_manual(values = colours, aesthetics = c("colour", "fill")),
    theme_classic()
  )
}
sd_15 <- filter(heterogeneity, TM == 15,
                quantity == "cross-country SD of slope")

p_levels <- ggplot() +
  geom_line(data = country_curves, aes(TM, estimate, group = GID_1),
            colour = colours[[1]], alpha = 0.35, linewidth = 0.3) +
  curve_layers("Growth relative to optimum") +
  labs(
    x = NULL, y = "Growth relative to optimum", colour = NULL, fill = NULL,
    title = "BHM quadratic with country random slopes (PWT, brms)",
    subtitle = sprintf(
      paste0("Year FE, country intercepts with Mundlak terms (= country FE); ",
             "correlated random slopes on T and T\u00b2 by country.\nThin ",
             "lines: country curves over each country's observed ",
             "temperatures. Cross-country SD of dg/dT at 15\u00b0C:\n%.4f ",
             "[%.4f, %.4f]. %d countries, %d-%d."),
      sd_15$estimate, sd_15$conf_low, sd_15$conf_high,
      nrow(countries), min(panel$year), max(panel$year)
    )
  ) +
  theme(legend.position = "none")

p_slopes <- ggplot() +
  geom_segment(data = countries,
               aes(x = TM_mean, xend = TM_mean, y = ols_slope,
                   yend = mixed_slope),
               colour = "grey65", linewidth = 0.25) +
  geom_point(data = countries, aes(TM_mean, ols_slope, size = weight),
             shape = 1, colour = "grey45") +
  geom_point(data = countries, aes(TM_mean, mixed_slope, size = weight),
             colour = colours[[1]], alpha = 0.75) +
  curve_layers("Marginal effect (dg / dT)") +
  scale_size_area(max_size = 4, guide = "none") +
  coord_cartesian(ylim = stats::quantile(countries$ols_slope, c(0.02, 0.98))) +
  labs(
    x = "Temperature (°C)", y = "Marginal effect (dg / dT)",
    colour = NULL, fill = NULL,
    caption = paste0(
      "Lower panel, at each country's mean temperature: hollow = country OLS ",
      "slope, filled = mixed-model slope, line = shrinkage; size = ",
      "identifying variance.\nClipped to the 2nd-98th percentile of OLS ",
      "slopes. Mixed band: 95% posterior interval; FE band: 95% CI clustered ",
      "by country."
    )
  ) +
  theme(legend.position = "bottom")

p <- p_levels / p_slopes + plot_layout(heights = c(3, 2))
ggsave(file.path(OUT, "bhm_mixed_slopes.png"), p,
       width = 10, height = 9, dpi = 200)

# Console summary -----------------------------------------------------------------

options(width = 160)
cat("\nR-hat:\n")
print(round(convergence, 3))
cat("\nPosterior summary, x = (TM - ", T_CENTRE, ") / ", T_SCALE, ":\n",
    sep = "")
print(t(vapply(draws, summarise_draws, numeric(3))), digits = 3)
cat("\nFE comparison (clustered SE):\n")
print(fixest::coeftable(m_fe), digits = 3)
cat("\nPopulation slope and cross-country SD of dg/dT (per degree):\n")
print(heterogeneity, digits = 3)
cat(sprintf(
  "\nSD of country slopes at own mean T: OLS %.4f, mixed %.4f\n",
  sd(countries$ols_slope), sd(countries$mixed_slope)
))
cat("\nWrote figure and tables to ", OUT, "\n", sep = "")
