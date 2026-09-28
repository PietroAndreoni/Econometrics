# Project the growth effect of a climate trajectory through a fitted model.
#
# The caller builds G, the change in every climate regressor along the
# trajectory relative to a counterfactual (typically climate frozen at its
# initial value), one row per year and one named column per coefficient. For a
# model with TM, TM_2 and a deviation term, for example:
#   G <- cbind(TM = T_traj - T0, TM_2 = T_traj^2 - T0^2,
#              dev = abs(T_traj - roll_mean_pad(T_traj, T0, 30)) / sd_T)
#   proj <- project_effects(m, G)
#   plot_projection(proj)

# Rolling mean over the `window` years ENDING AT t-1, with the pre-trajectory
# years held at x0. Lagged like the estimation baselines, so year t never
# enters its own reference mean.
roll_mean_pad <- function(x, x0, window) {
  w <- if (is.na(window) || window < 1) 1 else window
  padded <- c(rep(x0, w), utils::head(x, -1))
  vapply(seq_along(x), function(t) mean(padded[t:(t + w - 1)]), numeric(1))
}

# Constant trend starting from x0, so year 1 equals x0.
linear_traj <- function(x0, trend, n) x0 + trend * (seq_len(n) - 1)

# Changes that reach `total` immediately, or linearly over `ramp` years.
step_delta <- function(n, total) rep(total, n)
ramp_delta <- function(n, total, ramp) total * pmin(seq_len(n) / ramp, 1)

# Annual growth effect G %*% b and its running sum (the log GDP level
# deviation), with delta-method intervals. The cumulative intervals carry the
# covariance between years rather than treating them as independent.
project_effects <- function(model, G, years = seq_len(nrow(G)),
                            conf_level = 0.95, label = NULL) {
  G <- as.matrix(G)
  G_cum <- apply(G, 2, cumsum)
  if (is.null(dim(G_cum))) {
    G_cum <- matrix(G_cum, nrow = 1, dimnames = list(NULL, colnames(G)))
  }

  out <- bind_rows(
    tibble::tibble(year = years, quantity = "Annual growth effect") %>%
      bind_cols(linear_combination(model, G, conf_level)),
    tibble::tibble(year = years, quantity = "Cumulative effect") %>%
      bind_cols(linear_combination(model, G_cum, conf_level))
  ) %>%
    mutate(quantity = factor(
      quantity,
      levels = c("Annual growth effect", "Cumulative effect")
    ))
  if (!is.null(label)) out$trajectory <- label
  out
}

# One panel per quantity; several trajectories (bind_rows of project_effects()
# outputs with a `trajectory` label) are overlaid by colour.
plot_projection <- function(projection, conf_level = 0.95) {
  has_label <- "trajectory" %in% names(projection)
  mapping <- if (has_label) {
    aes(x = year, y = estimate, colour = trajectory, fill = trajectory)
  } else {
    aes(x = year, y = estimate)
  }

  ggplot(projection, mapping) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.12,
                colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~quantity, ncol = 1, scales = "free_y") +
    labs(
      x = "Year",
      y = expression(Delta * " log GDP per capita"),
      colour = "Trajectory",
      fill = "Trajectory",
      caption = paste0(
        "Relative to the counterfactual climate; ribbons = ",
        round(100 * conf_level), "% CI"
      )
    ) +
    theme_classic()
}
