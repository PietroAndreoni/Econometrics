# Heterogeneous bin responses: Carleton-style interactions with income and
# long-run climate, and adaptation through a spline in lagged income.

# ---- Carleton-style income x climate heterogeneity ---------------------------

# Nine evaluation points: regions are split into terciles of sample-average
# income and of long-run climate, and each cell is summarised by its mean
# centred income and centred climate.
make_carleton_scenarios <- function(
    model_data,
    climate_var,
    climate_centered_var,
    climate_group_labels,
    income_var = "log_gdp_av",
    income_centered_var = "log_gdp_av_c",
    income_group_labels = c("Low income", "Middle income", "High income"),
    outcome = "dlgrp_pc_usd",
    id = "GID_1"
) {
  model_data %>%
    ungroup() %>%
    filter(!is.na(.data[[outcome]])) %>%
    distinct(
      .data[[id]],
      .data[[income_var]],
      .data[[income_centered_var]],
      .data[[climate_var]],
      .data[[climate_centered_var]]
    ) %>%
    mutate(
      income_group = factor(
        income_group_labels[dplyr::ntile(.data[[income_var]], 3L)],
        levels = income_group_labels
      ),
      climate_group = factor(
        climate_group_labels[dplyr::ntile(.data[[climate_var]], 3L)],
        levels = climate_group_labels
      )
    ) %>%
    group_by(income_group, climate_group) %>%
    summarise(
      income_value = mean(.data[[income_centered_var]]),
      climate_value = mean(.data[[climate_centered_var]]),
      regions = n(),
      .groups = "drop"
    )
}

# Response of every bin of `bin_var` at each scenario, from a model containing
# i(bin_var), i(bin_var, income_var) and i(bin_var, climate_var).
carleton_bin_response <- function(
    fitted_model,
    model_data,
    bin_var,
    reference_bin,
    income_var,
    climate_var,
    scenarios
) {
  bin_levels <- levels(droplevels(model_data[[bin_var]]))

  purrr::map_dfr(seq_len(nrow(scenarios)), function(scenario_id) {
    scenario <- scenarios[scenario_id, ]

    purrr::map_dfr(bin_levels, function(bin_level) {
      if (identical(bin_level, reference_bin)) {
        return(tibble::tibble(bin = bin_level, estimate = 0, std_error = 0))
      }
      main_coefficient <- paste0(bin_var, "::", bin_level)
      weights <- stats::setNames(
        c(1, scenario$income_value[[1]], scenario$climate_value[[1]]),
        c(
          main_coefficient,
          paste0(main_coefficient, ":", income_var),
          paste0(main_coefficient, ":", climate_var)
        )
      )
      combination <- linear_combination(fitted_model, weights)
      tibble::tibble(
        bin = bin_level,
        estimate = combination$estimate,
        std_error = combination$std_error
      )
    }) %>%
      mutate(
        income_group = scenario$income_group[[1]],
        climate_group = scenario$climate_group[[1]],
        regions = scenario$regions[[1]]
      )
  }) %>%
    mutate(
      bin = factor(bin, levels = bin_levels),
      conf_low = estimate - stats::qnorm(0.975) * std_error,
      conf_high = estimate + stats::qnorm(0.975) * std_error
    )
}

plot_carleton_response <- function(plot_data, title, x_label) {
  ggplot(plot_data, aes(x = bin, y = estimate, group = 1)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 0.55, color = "#0072B2") +
    geom_errorbar(
      aes(ymin = conf_low, ymax = conf_high),
      width = 0.18,
      linewidth = 0.45,
      color = "#0072B2"
    ) +
    geom_point(size = 1.8, color = "#0072B2") +
    facet_grid(income_group ~ climate_group) +
    labs(
      title = title,
      subtitle = paste(
        "Columns: long-run climate tercile; rows: sample-average income tercile"
      ),
      x = x_label,
      y = "Estimated effect relative to reference bin"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey92", color = "grey70")
    )
}

# ---- Adaptation through lagged income ----------------------------------------

# Natural-spline basis of `var`, stored as scalar columns <prefix>1..<prefix>df
# so fixest can interact each one with i() bins (ns() returns a matrix, which
# i() cannot multiply directly). Rows where `var` is not finite get NA. Returns
# the augmented data, the fitted basis (for predict() at new values) and the
# column names.
add_spline_basis <- function(data, var, df, prefix = paste0(var, "_spline_")) {
  valid <- is.finite(data[[var]])
  basis <- matrix(NA_real_, nrow = nrow(data), ncol = df)
  basis_fit <- splines::ns(data[[var]][valid], df = df)
  basis[valid, ] <- basis_fit
  columns <- paste0(prefix, seq_len(df))
  colnames(basis) <- columns

  list(
    data = bind_cols(data, tibble::as_tibble(basis)),
    basis = basis_fit,
    columns = columns
  )
}

# Response of every signed bin of `bin_var` at each income scenario, from a
# model containing i(bin_var) and i(bin_var, <spline column>) for every column
# of `spline_values` (one row per scenario, from predict() on the fitted basis).
adaptation_response_data <- function(
    fitted_model,
    bin_var,
    model_data,
    spline_values,
    scenarios,
    spline_columns = colnames(spline_values)
) {
  bin_values <- sort(unique(model_data[[bin_var]][!is.na(model_data[[bin_var]])]))

  purrr::map_dfr(seq_len(nrow(scenarios)), function(income_id) {
    purrr::map_dfr(bin_values, function(bin_value) {
      if (bin_value == 0) {
        return(tibble::tibble(bin = bin_value, estimate = 0, std_error = 0))
      }
      main_coefficient <- paste0(bin_var, "::", bin_value)
      weights <- stats::setNames(
        c(1, spline_values[income_id, ]),
        c(main_coefficient, paste0(main_coefficient, ":", spline_columns))
      )
      combination <- linear_combination(fitted_model, weights)
      tibble::tibble(
        bin = bin_value,
        estimate = combination$estimate,
        std_error = combination$std_error
      )
    }) %>%
      mutate(
        income_label = scenarios$income_label[[income_id]],
        income = scenarios$income[[income_id]]
      )
  }) %>%
    mutate(
      conf_low = estimate - stats::qnorm(0.975) * std_error,
      conf_high = estimate + stats::qnorm(0.975) * std_error
    )
}

plot_adaptation_response <- function(plot_data, climate_variable) {
  response_dodge <- position_dodge(width = 0.12)

  ggplot(
    plot_data,
    aes(x = bin, y = estimate, color = income_label, group = income_label)
  ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(
      aes(ymin = conf_low, ymax = conf_high),
      width = 0.06,
      linewidth = 0.5,
      position = response_dodge
    ) +
    geom_point(size = 2, position = response_dodge) +
    scale_x_continuous(breaks = sort(unique(plot_data$bin))) +
    labs(
      title = paste(climate_variable, "response by lagged income"),
      x = "Signed anomaly bin (standard deviations)",
      y = "Estimated response",
      color = "Lagged GDP per capita"
    ) +
    theme_classic() +
    theme(legend.position = "right")
}
