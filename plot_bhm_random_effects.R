# The Burke-Hsiang-Miguel (2015) quadratic on the national PWT panel, estimated
# by random effects and by fixed effects with the same two-way structure.
#
# Both models regress growth on TM + TM^2 + RR + RR^2 with year fixed effects
# and no country-specific trends:
#   FE: country fixed effects (within estimator), fixest.
#   RE: country random intercepts (Amemiya GLS), plm. Year effects enter as
#       dummies, so only the country effect is random. Amemiya takes both
#       variance components from within residuals; Swamy-Arora's between
#       regression is singular here, because the country means of the year
#       dummies take only as many values as there are coverage patterns.
# Standard errors are clustered by country in both (Arellano for plm).
#
# RE is consistent only if the country effect is uncorrelated with climate. The
# robust Hausman test is the Mundlak version: add the country means of the four
# climate regressors to the RE model and Wald-test them jointly with the
# clustered VCOV. (The classical Hausman test is not used: it assumes RE is
# efficient, which fails with clustered errors, and plm's version also
# contrasts the year dummies.)
#
# The figure draws each curve relative to its own optimum (to the sample-mean
# temperature if the curve is not concave), with 95% delta-method intervals,
# and the implied marginal effect b1 + 2 b2 T below.
#
# Sample: plot_bhm_country_pieces.R, i.e. PWT from its first year.
# Output: results/bhm_random_effects/.

source("load_functions.R")

SPEC <- main_spec(econ_year_min = 1900L)
FE_TWO_WAY <- "GID_1 + year"
TERMS <- "TM + TM^2 + RR + RR^2"
CLIMATE <- c("TM", "I(TM^2)", "RR", "I(RR^2)")
OUT <- file.path("results", "bhm_random_effects")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

needed <- c(SPEC$outcome, "TM", "RR")
panel <- build_main_dat("PWT", SPEC) %>%
  filter(if_all(all_of(needed), is.finite)) %>%
  group_by(GID_1) %>%
  mutate(
    TM_bar = mean(TM), TM2_bar = mean(TM^2),
    RR_bar = mean(RR), RR2_bar = mean(RR^2)
  ) %>%
  ungroup() %>%
  arrange(GID_1, year)  # pdata.frame order, so theta lines up by row.
stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
cat(sprintf("PWT: %d country-years, %d countries, %d-%d\n", nrow(panel),
            n_distinct(panel$GID_1), min(panel$year), max(panel$year)))

# Estimation ----------------------------------------------------------------------

spec_fe <- modifyList(SPEC, list(fixed_effects = FE_TWO_WAY))
m_fe <- fit_main(TERMS, panel, spec_fe)
m_bhm <- fit_main(TERMS, panel, SPEC)  # Main spec: adds quadratic trends.

pdata <- plm::pdata.frame(as.data.frame(panel), index = c("GID_1", "year"))
re_formula <- stats::as.formula(paste(
  SPEC$outcome, "~ TM + I(TM^2) + RR + I(RR^2) + factor(year)"
))
m_re <- plm::plm(re_formula, data = pdata, model = "random",
                 effect = "individual", random.method = "amemiya")
m_within <- plm::plm(re_formula, data = pdata, model = "within",
                     effect = "individual")
cluster_vcov <- function(m) {
  plm::vcovHC(m, method = "arellano", type = "HC1", cluster = "group")
}
V_re <- cluster_vcov(m_re)
stopifnot(nobs(m_fe) == nrow(panel), nobs(m_re) == nrow(panel))
stopifnot(isTRUE(all.equal(coef(m_within)[CLIMATE], coef(m_fe)[CLIMATE],
                           tolerance = 1e-6)))

theta <- plm::ercomp(m_re)$theta
variance_shares <- plm::ercomp(m_re)$sigma2
theta_obs <- if (length(theta) == 1L) rep(theta, nrow(panel)) else theta
stopifnot(length(theta_obs) == nrow(panel))

# RE GLS as OLS on data quasi-demeaned with m_re's theta, clustered by country.
# plm cannot refit Amemiya with the time-invariant Mundlak means, so the
# Mundlak model reuses m_re's variance components this way.
quasi_demeaned_gls <- function(formula) {
  X <- stats::model.matrix(formula, panel)
  y <- panel[[SPEC$outcome]]
  qd <- function(v) v - theta_obs * stats::ave(v, panel$GID_1)
  X_qd <- apply(X, 2, qd)
  fit <- stats::lm(qd(y) ~ X_qd - 1)
  names(fit$coefficients) <- colnames(X)
  list(b = coef(fit),
       V = sandwich::vcovCL(fit, cluster = panel$GID_1, type = "HC1") %>%
         `dimnames<-`(list(colnames(X), colnames(X))))
}
stopifnot(isTRUE(all.equal(quasi_demeaned_gls(re_formula)$b[CLIMATE],
                           coef(m_re)[CLIMATE], tolerance = 1e-6)))

# Mundlak: RE plus the country means; the means' coefficients are zero under RE.
means <- c("TM_bar", "TM2_bar", "RR_bar", "RR2_bar")
m_mundlak <- quasi_demeaned_gls(stats::update(re_formula, . ~ . + TM_bar +
                                                TM2_bar + RR_bar + RR2_bar))
b_means <- m_mundlak$b[means]
V_means <- m_mundlak$V[means, means]
mundlak_stat <- as.numeric(t(b_means) %*% solve(V_means, b_means))

# Each mean's coefficient is the between-minus-within gap for that regressor.
tests <- tibble::tibble(
  test = "Mundlak (cluster-robust)",
  statistic = mundlak_stat,
  df = length(means),
  p_value = stats::pchisq(mundlak_stat, length(means), lower.tail = FALSE)
)
mundlak_means <- tibble::tibble(
  term = means, estimate = unname(b_means),
  std_error = sqrt(diag(V_means)),
  p_value = 2 * stats::pnorm(-abs(estimate / std_error))
)

# Coefficients ----------------------------------------------------------------------

coef_rows <- function(b, V, model) {
  tibble::tibble(
    model = model, term = CLIMATE,
    estimate = unname(b[CLIMATE]),
    std_error = sqrt(diag(V)[CLIMATE])
  ) %>%
    mutate(p_value = 2 * stats::pnorm(-abs(estimate / std_error)))
}
coefficients <- bind_rows(
  coef_rows(coef(m_re), V_re, "RE: country RE + year FE"),
  coef_rows(coef(m_fe), vcov(m_fe), "FE: country + year FE"),
  coef_rows(coef(m_bhm), vcov(m_bhm),
            "FE: country + year FE + quadratic trends")
)
write_out(coefficients, "coefficients")
write_out(tests, "hausman_tests")
write_out(mundlak_means, "mundlak_means")

# Curves ------------------------------------------------------------------------------

t_grid <- seq(floor(min(panel$TM)), ceiling(max(panel$TM)), by = 0.25)

curves_for <- function(model, label, vcov = NULL) {
  b <- coef(model)
  concave <- b[["I(TM^2)"]] < 0
  t_ref <- if (concave) -b[["TM"]] / (2 * b[["I(TM^2)"]]) else mean(panel$TM)
  levels <- tibble::tibble(TM = t_grid) %>%
    bind_cols(linear_combination(
      model, cbind(TM = t_grid - t_ref, `I(TM^2)` = t_grid^2 - t_ref^2),
      vcov = vcov
    ))
  slopes <- tibble::tibble(TM = t_grid) %>%
    bind_cols(linear_combination(
      model, cbind(TM = 1, `I(TM^2)` = 2 * t_grid), vcov = vcov
    ))
  bind_rows(
    mutate(levels, quantity = "Growth relative to optimum"),
    mutate(slopes, quantity = "Marginal effect (dg / dT)")
  ) %>%
    mutate(model = sprintf("%s (optimum %s)", label,
                           if (concave) sprintf("%.1f°C", t_ref) else "none"))
}

curves <- bind_rows(
  curves_for(m_re, "Random effects", V_re),
  curves_for(m_fe, "Fixed effects")
) %>%
  mutate(model = factor(model, levels = unique(model)))
write_out(curves, "curves")

p <- ggplot(curves, aes(TM, estimate, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~quantity, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#2C6FB7", "#C0392B"), aesthetics = c("colour", "fill")) +
  labs(
    x = "Temperature (°C)", y = NULL, colour = NULL, fill = NULL,
    title = "BHM quadratic: random vs fixed country effects (PWT)",
    subtitle = sprintf(
      paste0("Both with year fixed effects, no country trends; 95%% CI ",
             "clustered by country. %d countries, %d-%d.\n",
             "Mundlak test of RE: chi2(%d) = %.1f, p = %.3g."),
      n_distinct(panel$GID_1), min(panel$year), max(panel$year),
      tests$df[1], tests$statistic[1], tests$p_value[1]
    )
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "bhm_random_effects.png"), p,
       width = 9, height = 8, dpi = 200)

# Console summary -----------------------------------------------------------------

options(width = 160)
print(coefficients %>%
        mutate(across(c(estimate, std_error), ~signif(.x, 3)),
               p_value = signif(p_value, 3)),
      n = Inf)
cat("\nRE quasi-demeaning theta (distribution across countries):\n")
print(summary(theta))
cat("Variance components (idios, id):\n")
print(variance_shares)
print(tests)
print(mundlak_means)
cat("\nWrote figure and tables to ", OUT, "\n", sep = "")
