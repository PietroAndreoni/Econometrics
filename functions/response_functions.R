# Extract and plot estimated climate-response functions from fitted models.

# ---- Distributed lags ---------------------------------------------------------

# Each lag coefficient and the cumulative response through that lag, for every
# l(term, lags) block in the model. The cumulative intervals retain the
# cross-lag covariances of the model's VCOV. Terms and lag ranges are discovered
# from the coefficient names, and specifications written as term + l(term, 1:k)
# are supported (fixest prints lag 0 without the l(., 0) wrapper).
distributed_lag_effects <- function(model, conf_level = 0.95) {
  z_critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  coefficients <- stats::coef(model)
  V <- stats::vcov(model)
  coefficient_names <- names(coefficients)
  names_normalized <- gsub("[[:space:]]", "", coefficient_names)

  lag_name_pattern <- "^l\\((.*),([0-9]+)\\)$"
  lagged_positions <- grep(lag_name_pattern, names_normalized)
  if (!length(lagged_positions)) {
    stop("The model contains no coefficients named l(term, lag).")
  }

  lag_term_index <- tibble::tibble(
    expression = sub(lag_name_pattern, "\\1", names_normalized[lagged_positions]),
    lag = as.integer(sub(lag_name_pattern, "\\2", names_normalized[lagged_positions])),
    term = coefficient_names[lagged_positions]
  )

  purrr::map_dfr(unique(lag_term_index$expression), function(lag_expression) {
    lag_terms <- lag_term_index %>%
      filter(expression == lag_expression) %>%
      arrange(lag)

    if (!0L %in% lag_terms$lag) {
      lag_zero_position <- match(lag_expression, names_normalized)
      if (!is.na(lag_zero_position)) {
        lag_terms <- bind_rows(
          tibble::tibble(
            expression = lag_expression,
            lag = 0L,
            term = coefficient_names[lag_zero_position]
          ),
          lag_terms
        ) %>%
          arrange(lag)
      }
    }

    if (anyDuplicated(lag_terms$lag)) {
      stop("Duplicate coefficients found for ", lag_expression, ".")
    }

    purrr::map_dfr(lag_terms$lag, function(current_lag) {
      current_term <- lag_terms$term[lag_terms$lag == current_lag]
      cumulative_terms <- lag_terms$term[lag_terms$lag <= current_lag]

      tibble::tibble(
        expression = lag_expression,
        lag = current_lag,
        horizon = c("Effect at lag", "Cumulative through lag"),
        effect = c(
          unname(coefficients[current_term]),
          sum(coefficients[cumulative_terms])
        ),
        standard_error = c(
          sqrt(max(as.numeric(V[current_term, current_term]), 0)),
          sqrt(max(sum(V[cumulative_terms, cumulative_terms, drop = FALSE]), 0))
        )
      )
    })
  }) %>%
    mutate(
      expression = factor(expression, levels = unique(lag_term_index$expression)),
      horizon = factor(
        horizon,
        levels = c("Effect at lag", "Cumulative through lag")
      ),
      lower = effect - z_critical * standard_error,
      upper = effect + z_critical * standard_error
    )
}

plot_distributed_lags <- function(effects, conf_level = 0.95,
                                  y_label = "Effect on log GDP per-capita growth") {
  ggplot(effects, aes(x = lag, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_grid(expression ~ horizon, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(effects$lag))) +
    labs(
      x = "Lag (years)",
      y = y_label,
      title = "Distributed-lag coefficient paths",
      caption = paste0(
        "Ribbons = ", round(100 * conf_level),
        "% clustered CI; cumulative intervals include cross-lag covariance"
      )
    ) +
    theme_classic()
}

# ---- Bin coefficients -------------------------------------------------------

# Coefficients of i(<bin_var>) terms, one row per bin, with the omitted
# reference bin added at zero. `bin_vars` are matched as <bin_var>::<level>, so
# the default recovers both signed anomaly bin sets of a joint model.
bin_coefs <- function(model, bin_vars = c("TM_bin_signed", "RR_bin_signed"),
                      reference_bin = "0", conf_level = 0.95) {
  table <- as.data.frame(fixest::coeftable(model))
  intervals <- stats::confint(model, level = conf_level)
  names_all <- rownames(table)

  bind_rows(lapply(bin_vars, function(bin_var) {
    selected <- startsWith(names_all, paste0(bin_var, "::")) &
      !grepl(":", substring(names_all, nchar(bin_var) + 3L), fixed = TRUE)
    if (!any(selected)) {
      stop("No ", bin_var, " coefficients found in the model.")
    }
    bind_rows(
      tibble::tibble(
        bin_var = bin_var,
        bin_label = substring(names_all[selected], nchar(bin_var) + 3L),
        estimate = table[selected, 1],
        std_error = table[selected, 2],
        p_value = table[selected, 4],
        ci_low = intervals[selected, 1],
        ci_high = intervals[selected, 2]
      ),
      tibble::tibble(
        bin_var = bin_var, bin_label = reference_bin, estimate = 0,
        std_error = NA_real_, p_value = NA_real_, ci_low = 0, ci_high = 0
      )
    )
  }))
}

# Coefficients of i(<bin_var>, <group_var>) interactions, one row per bin and
# group, with the reference bin added at zero for each of `group_levels`.
interacted_bin_coefs <- function(model, bin_var, group_var, group_levels,
                                 reference_bin = "0") {
  if (!inherits(model, "fixest")) {
    stop("model must be a fixest object.")
  }
  interaction_marker <- paste0(":", group_var, "::")
  model_coef <- stats::coef(model)
  coef_names <- names(model_coef)
  selected <- startsWith(coef_names, paste0(bin_var, "::")) &
    grepl(interaction_marker, coef_names, fixed = TRUE)

  if (!any(selected)) {
    stop("No ", bin_var, " by ", group_var, " coefficients found in the model.")
  }

  coef_names <- coef_names[selected]
  coefficient_tail <- substring(coef_names, nchar(bin_var) + 3L)
  model_ci <- stats::confint(model)

  bind_rows(
    data.frame(
      estimate = unname(model_coef[coef_names]),
      ci_low = model_ci[coef_names, 1],
      ci_high = model_ci[coef_names, 2],
      bin_label = sub(paste0(interaction_marker, ".*$"), "", coefficient_tail),
      group = sub(paste0("^.*", interaction_marker), "", coefficient_tail)
    ),
    data.frame(
      estimate = 0,
      ci_low = 0,
      ci_high = 0,
      bin_label = reference_bin,
      group = group_levels
    )
  )
}

# The i.select-th set of i() bins from every model of a split-sample
# fixest_multi, labelled by sample. As with fixest::iplot, the reference bin is
# included at zero.
split_bin_coefs <- function(models, i_select) {
  if (!inherits(models, "fixest_multi")) {
    stop("models must be a fixest_multi object.")
  }
  plot_parameters <- fixest::iplot(models, i.select = i_select, only.params = TRUE)
  model_tree <- attr(models, "tree")

  plot_parameters$prms %>%
    transmute(
      estimate,
      ci_low,
      ci_high,
      bin_label = as.character(estimate_names),
      group = model_tree$sample[match(id, model_tree$id)]
    )
}

# Order bin labels numerically, e.g. "-0.75" before "0" before "0.75".
numeric_bin_factor <- function(bin_label) {
  bin_numeric <- suppressWarnings(as.numeric(bin_label))
  if (anyNA(bin_numeric)) {
    stop("Bin labels must be numeric.")
  }
  factor(
    as.character(bin_numeric),
    levels = as.character(sort(unique(bin_numeric)))
  )
}

# Dodged point-and-interval plot of bin coefficients by group. Expects columns
# bin_label, estimate, ci_low, ci_high and `group`.
plot_group_bins <- function(
    plot_data,
    colors,
    shapes,
    title,
    group_label = "Group",
    x_label = "Signed anomaly bin (standard deviations)",
    facet = FALSE,
    ylim = NULL
) {
  bin_dodge <- position_dodge(width = 0.65)

  p <- ggplot(
    plot_data,
    aes(x = bin_label, y = estimate, color = group, shape = group)
  ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.06,
      linewidth = 0.5,
      position = bin_dodge
    ) +
    geom_point(size = 2, position = bin_dodge) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    labs(
      title = title,
      x = x_label,
      y = "Coefficient estimate",
      color = group_label,
      shape = group_label
    ) +
    theme_classic() +
    theme(legend.position = "right")

  if (facet) p <- p + facet_wrap(group ~ .)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ---- Parametric response curves ----------------------------------------------

# Delta-method response curve of one or more coefficients over a grid. `basis`
# is a matrix with one row per grid point and one column per term in `terms`.
# If `ref_basis` is given, the curve is re-centred on the average basis over
# those reference observations (e.g. the rows in the omitted bin), which makes
# it comparable to bin coefficients.
response_curve <- function(model, terms, basis, ref_basis = NULL,
                           conf_level = 0.95) {
  basis <- as.matrix(basis)
  if (!is.null(ref_basis)) {
    basis <- sweep(basis, 2, colMeans(as.matrix(ref_basis), na.rm = TRUE))
  }
  colnames(basis) <- terms
  linear_combination(model, basis, conf_level = conf_level)
}
