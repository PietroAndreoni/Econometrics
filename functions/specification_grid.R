# Estimate one specification over a grid of data choices (economic dataset,
# GDP definition, climate source, weighting, ...) on a maximum and a common
# sample, and summarise how much the bin coefficients depend on each choice.
#
# Conventions shared by every function here:
#   specs  one row per grid cell: an `id` plus one column per choice (factor);
#   fits   list(max = list(<id> = fixest), common = list(<id> = fixest));
#   data   the data each fit was estimated on, in the same nested layout.

KEY_COLS <- c("GID_0", "GID_1", "year")

# Rows a model actually used; fixest::obs() indexes the data it was fit on.
grid_used_rows <- function(fits, data, sample, id) {
  data[[sample]][[id]][fixest::obs(fits[[sample]][[id]]), ]
}

# Unit-years present in every element of `key_list` (data frames of KEY_COLS).
intersect_keys <- function(key_list) {
  Reduce(function(a, b) inner_join(a, b, by = KEY_COLS), key_list) %>%
    distinct()
}

# Applies `f(fit, sample, id)` to every model and binds the rows, with the
# grid choices joined on.
.map_grid <- function(fits, specs, f) {
  bind_rows(lapply(names(fits), function(sample) {
    bind_rows(lapply(specs$id, function(id) {
      f(fits[[sample]][[id]], sample, id) %>%
        mutate(sample = sample, id = id, .before = 1)
    }))
  })) %>%
    left_join(specs, by = "id")
}

grid_sample_summary <- function(fits, data, specs) {
  .map_grid(fits, specs, function(fit, sample, id) {
    used <- grid_used_rows(fits, data, sample, id)
    tibble::tibble(
      observations = as.integer(stats::nobs(fit)),
      units = n_distinct(used$GID_1),
      countries = n_distinct(used$GID_0),
      first_year = min(used$year),
      last_year = max(used$year),
      within_r2 = as.numeric(fixest::fitstat(fit, "wr2")[[1]])
    )
  })
}

# Signed-bin coefficients with the omitted central bin at zero.
grid_bin_coefficients <- function(fits, specs) {
  .map_grid(fits, specs, function(fit, sample, id) bin_coefs(fit)) %>%
    mutate(
      variable = if_else(bin_var == "TM_bin_signed", "Temperature", "Precipitation"),
      bin = as.numeric(bin_label),
      significant = !is.na(p_value) & p_value < 0.05
    )
}

# Every coefficient that is not a signed bin (the BHM base terms).
grid_base_terms <- function(fits, specs) {
  .map_grid(fits, specs, function(fit, sample, id) {
    tidy_fixest(fit, label = id) %>% filter(!grepl("_bin_signed::", term))
  })
}

# Joint Wald tests that each set of signed bins is zero.
grid_joint_tests <- function(fits, specs) {
  .map_grid(fits, specs, function(fit, sample, id) {
    bind_rows(
      wald_row(fit, "^TM_bin_signed::", label = "Temperature bins = 0"),
      wald_row(fit, "^RR_bin_signed::", label = "Precipitation bins = 0")
    )
  })
}

# Per bin: spread of the estimates across specifications, and the share of that
# spread explained by each factor's level means (share_<factor>).
grid_dispersion <- function(coefficients, factors) {
  coefficients %>%
    filter(bin != 0) %>%
    group_by(sample, variable, bin) %>%
    group_modify(function(d, key) {
      total <- sum((d$estimate - mean(d$estimate))^2)
      shares <- vapply(factors, function(f) {
        means <- stats::ave(d$estimate, d[[f]])
        sum((means - mean(d$estimate))^2) / total
      }, numeric(1))
      bind_cols(
        tibble::tibble(
          specifications = nrow(d),
          mean_estimate = mean(d$estimate),
          sd_estimate = stats::sd(d$estimate),
          min_estimate = min(d$estimate),
          max_estimate = max(d$estimate),
          share_positive = mean(d$estimate > 0),
          share_significant = mean(d$significant)
        ),
        tibble::as_tibble(as.list(stats::setNames(shares, paste0("share_", factors))))
      )
    }) %>%
    ungroup()
}

# Mean estimate per bin for each level of each factor.
grid_factor_means <- function(coefficients, factors) {
  bind_rows(lapply(factors, function(f) {
    coefficients %>%
      filter(bin != 0) %>%
      group_by(sample, variable, bin, factor_level = .data[[f]]) %>%
      summarise(mean_estimate = mean(estimate), .groups = "drop") %>%
      mutate(factor = f, .before = factor_level)
  }))
}

# Counts of values in the estimation rows that deserve a look. Source-specific
# missing-value normalization, including ERA5 concurrent-population zeros, is
# applied upstream in load_climate_series().
grid_data_quality <- function(fits, data, specs) {
  .map_grid(fits, specs, function(fit, sample, id) {
    used <- grid_used_rows(fits, data, sample, id)
    tibble::tibble(
      rows = nrow(used),
      zero_precipitation = sum(used$RR == 0),
      TM_min = min(used$TM),
      TM_max = max(used$TM),
      RR_min_m = min(used$RR),
      RR_max_m = max(used$RR),
      share_extreme_TM_bins = mean(abs(used$TM_bin_signed) == 2.75),
      share_extreme_RR_bins = mean(abs(used$RR_bin_signed) == 2.75)
    )
  })
}

# Plots -------------------------------------------------------------------------
#
# Reference categorical palette, in fixed order; every categorical encoding
# also gets a marker shape or line type so identity never rests on colour.

GRID_PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                  "#008300", "#4a3aa7", "#e34948")
GRID_SHAPES <- c(16, 17, 15, 18, 3, 4, 8, 1)
GRID_INK <- "#1f1f1e"
GRID_MUTED <- "#6b6a65"

# One y-range per variable, so a few thin extreme bins do not flatten every
# other curve; estimates beyond it are drawn at the edge and labelled.
SIGNED_BIN_Y_LIMITS <- list(Temperature = c(-4, 4), Precipitation = c(-4, 2))

theme_results <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text = element_text(colour = GRID_INK),
      axis.text = element_text(colour = GRID_MUTED),
      plot.subtitle = element_text(colour = GRID_MUTED),
      plot.caption = element_text(colour = GRID_MUTED, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#e6e5df", linewidth = 0.3),
      strip.text = element_text(face = "bold", colour = GRID_INK),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Named colours and shapes for the levels of a factor, in palette order.
grid_colours <- function(levels) stats::setNames(GRID_PALETTE[seq_along(levels)], levels)
grid_shapes <- function(levels) stats::setNames(GRID_SHAPES[seq_along(levels)], levels)

# Coefficients in percentage points with ordered factors, ready to plot.
grid_plot_data <- function(coefficients, factor_levels) {
  out <- coefficients %>%
    mutate(
      variable = factor(variable, levels = c("Temperature", "Precipitation")),
      across(c(estimate, ci_low, ci_high), ~ 100 * .x)
    )
  for (f in names(factor_levels)) {
    out[[f]] <- factor(out[[f]], levels = factor_levels[[f]])
  }
  out
}

# The continuous BHM terms that accompany the signed-deviation bins. Keep the
# two level effects next to their corresponding long-run-climate interactions.
BASE_TERM_LABELS <- c(
  TM = "TM",
  "TM:mean_TM_all" = "TM x mean(TM)",
  RR = "RR",
  "RR:mean_RR_all" = "RR x mean(RR)"
)

# Base-term coefficients in percentage points with ordered specification
# factors, ready for the robustness plots below.
grid_base_plot_data <- function(base_terms, factor_levels) {
  out <- base_terms %>%
    filter(term %in% names(BASE_TERM_LABELS)) %>%
    mutate(
      term_label = factor(
        unname(BASE_TERM_LABELS[term]),
        levels = unname(BASE_TERM_LABELS)
      ),
      across(c(estimate, conf_low, conf_high), ~100 * .x)
    )
  for (f in names(factor_levels)) {
    out[[f]] <- factor(out[[f]], levels = factor_levels[[f]])
  }
  out
}

# The four non-deviation coefficients for every specification. Climate source
# is on the horizontal axis and in colour; an optional second choice such as
# weighting is encoded by marker shape. Each coefficient has its own y-scale.
plot_base_term_grid <- function(plot_data, facet, colour, colour_name,
                                shape = NULL, shape_values = NULL,
                                shape_labels = waiver(), shape_name = NULL,
                                title = NULL, subtitle = NULL,
                                caption = NULL) {
  colour_levels <- levels(plot_data[[colour]])
  shape_factor <- if (is.null(shape)) colour else shape
  shape_levels <- levels(plot_data[[shape_factor]])
  if (is.null(shape_values)) shape_values <- grid_shapes(shape_levels)
  dodge <- position_dodge(width = if (is.null(shape)) 0 else 0.55)
  group_cols <- unique(c(colour, shape_factor))

  ggplot(
    plot_data,
    aes(
      x = .data[[colour]], y = estimate,
      colour = .data[[colour]], shape = .data[[shape_factor]],
      group = interaction(!!!rlang::syms(group_cols))
    )
  ) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_linerange(
      aes(ymin = conf_low, ymax = conf_high),
      position = dodge, linewidth = 0.45, alpha = 0.7
    ) +
    geom_point(position = dodge, size = 2) +
    facet_grid(
      stats::as.formula(paste("term_label ~", facet)),
      scales = "free_y"
    ) +
    scale_colour_manual(
      values = grid_colours(colour_levels), name = colour_name
    ) +
    scale_shape_manual(
      values = shape_values, labels = shape_labels,
      name = if (is.null(shape)) colour_name else shape_name
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste(
        c(caption, "Each coefficient row has an independent vertical scale."),
        collapse = " "
      ),
      x = "Climate data",
      y = "Coefficient estimate (pp)"
    ) +
    theme_results() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Mean, interquartile range and min-max of each non-deviation coefficient across
# specifications, shown separately for maximum and common samples.
plot_base_term_range <- function(plot_data, sample_labels, title, subtitle) {
  range_data <- plot_data %>%
    group_by(sample, term_label) %>%
    summarise(
      mean = mean(estimate),
      low = min(estimate),
      high = max(estimate),
      q25 = stats::quantile(estimate, 0.25),
      q75 = stats::quantile(estimate, 0.75),
      .groups = "drop"
    ) %>%
    mutate(sample = factor(sample, levels = names(sample_labels)))
  sample_colours <- grid_colours(names(sample_labels))

  ggplot(range_data, aes(sample, mean, colour = sample, shape = sample)) +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_linerange(aes(ymin = low, ymax = high), linewidth = 0.5) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 2) +
    geom_point(size = 3, colour = "white", show.legend = FALSE) +
    geom_point(size = 1.8) +
    facet_wrap(~term_label, scales = "free_y", ncol = 2) +
    scale_colour_manual(
      values = sample_colours, labels = sample_labels, name = NULL
    ) +
    scale_shape_manual(
      values = grid_shapes(names(sample_labels)),
      labels = sample_labels, name = NULL
    ) +
    scale_x_discrete(labels = sample_labels) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste(
        "Thin line: min-max; thick line: interquartile range.",
        "Each coefficient panel has an independent vertical scale."
      ),
      x = NULL,
      y = "Coefficient estimate (pp)"
    ) +
    theme_results() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# Response curves of one variable: one column per `facet`, one colour (and
# shape) per `colour`, optionally one line type per `linetype`.
.bin_panel <- function(d, var, facet, colour, colour_values, colour_name,
                       linetype, linetype_values, linetype_labels,
                       linetype_name, show_x, y_limits) {
  limits <- y_limits[[var]]
  d <- d %>%
    filter(variable == var) %>%
    mutate(
      clipped = estimate < limits[1] | estimate > limits[2],
      shown = pmin(pmax(estimate, limits[1]), limits[2]),
      # Stack the labels of clipped estimates that share a bin.
      offset = as.integer(.data[[if (is.null(linetype)) colour else linetype]]) - 1L +
        as.integer(factor(bin)) %% 2L
    )
  group_cols <- c(colour, linetype)
  dodge <- position_dodge(width = 0.38)

  p <- ggplot(
    d,
    aes(bin, shown, colour = .data[[colour]], shape = .data[[colour]],
        group = interaction(!!!rlang::syms(group_cols)))
  )
  if (!is.null(linetype)) p <- p + aes(linetype = .data[[linetype]])

  p <- p +
    geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = dodge,
                   linewidth = 0.35, alpha = 0.55, show.legend = FALSE) +
    geom_line(position = dodge, linewidth = 0.5) +
    geom_point(data = function(x) filter(x, !clipped), position = dodge,
               size = 1.6) +
    geom_point(data = function(x) filter(x, clipped), position = dodge,
               size = 2.2, fill = "white", stroke = 0.8, shape = 21,
               show.legend = FALSE) +
    geom_text(
      data = function(x) filter(x, clipped),
      aes(label = sprintf("%+.1f", estimate),
          vjust = if_else(estimate > 0, 1.9, -1) + 1.5 * offset * sign(estimate)),
      position = dodge, size = 2.4, colour = GRID_INK, show.legend = FALSE
    ) +
    facet_grid(stats::as.formula(paste("variable ~", facet))) +
    coord_cartesian(ylim = limits) +
    scale_x_continuous(breaks = sort(unique(d$bin))) +
    scale_colour_manual(values = colour_values, name = colour_name) +
    scale_shape_manual(values = grid_shapes(names(colour_values)),
                       name = colour_name) +
    labs(
      x = if (show_x) "Signed anomaly bin (SD of the lagged 30-year climate)" else NULL,
      y = "Effect on growth (pp)"
    ) +
    theme_results() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  if (!is.null(linetype)) {
    p <- p + scale_linetype_manual(values = linetype_values,
                                   labels = linetype_labels, name = linetype_name)
  }
  p
}

# Temperature above precipitation, sharing one legend.
plot_bin_grid <- function(plot_data, facet, colour, colour_name,
                          linetype = NULL, linetype_values = NULL,
                          linetype_labels = waiver(), linetype_name = NULL,
                          title = NULL, subtitle = NULL, caption = NULL,
                          y_limits = SIGNED_BIN_Y_LIMITS) {
  colour_values <- grid_colours(levels(plot_data[[colour]]))
  panel <- function(var, show_x) {
    .bin_panel(plot_data, var, facet, colour, colour_values, colour_name,
               linetype, linetype_values, linetype_labels, linetype_name,
               show_x, y_limits)
  }
  (panel("Temperature", FALSE) / panel("Precipitation", TRUE)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = paste(
        c(caption, "Hollow markers: estimate outside the axis range, drawn at",
          "the edge and labelled with its value."),
        collapse = " "
      ),
      theme = theme_results()
    ) &
    theme(legend.position = "bottom")
}

# Mean, interquartile range and min-max of each bin across specifications,
# one colour per sample.
plot_bin_range <- function(plot_data, sample_labels, title, subtitle,
                           y_limits = SIGNED_BIN_Y_LIMITS) {
  range_data <- plot_data %>%
    group_by(sample, variable, bin) %>%
    summarise(
      mean = mean(estimate),
      low = min(estimate),
      high = max(estimate),
      q25 = stats::quantile(estimate, 0.25),
      q75 = stats::quantile(estimate, 0.75),
      .groups = "drop"
    ) %>%
    mutate(sample = factor(sample, levels = names(sample_labels)))
  sample_colours <- grid_colours(names(sample_labels))

  panel <- function(var) {
    d <- filter(range_data, variable == var)
    dodge <- position_dodge(width = 0.3)
    ggplot(d, aes(bin, mean, colour = sample, shape = sample)) +
      geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
      geom_linerange(aes(ymin = low, ymax = high), position = dodge,
                     linewidth = 0.4) +
      geom_linerange(aes(ymin = q25, ymax = q75), position = dodge,
                     linewidth = 1.6) +
      geom_point(position = dodge, size = 2.4, colour = "white",
                 show.legend = FALSE) +
      geom_point(position = dodge, size = 1.5) +
      coord_cartesian(ylim = y_limits[[var]]) +
      scale_x_continuous(breaks = sort(unique(d$bin))) +
      scale_colour_manual(values = sample_colours, labels = sample_labels,
                          name = NULL) +
      scale_shape_manual(values = grid_shapes(names(sample_labels)),
                         labels = sample_labels, name = NULL) +
      labs(title = var, x = "Signed anomaly bin (SD)",
           y = if (var == "Temperature") "Effect on growth (pp)" else NULL) +
      theme_results()
  }
  (panel("Temperature") | panel("Precipitation")) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = "Min-max lines running to the panel edge are truncated at the axis range.",
      theme = theme_results()
    ) &
    theme(legend.position = "bottom")
}

# Heatmap of a pairwise agreement measure between specification variants:
# columns `variant_1`, `variant_2`, `value`, plus any facet columns.
plot_agreement_heatmap <- function(d, facet_rows = NULL, facet_cols = NULL,
                                   fill_name, title,
                                   subtitle, limits = c(0, 1),
                                   label_format = function(x) sprintf("%.0f%%", 100 * x),
                                   fill_labels = scales::percent) {
  threshold <- mean(limits)
  ggplot(d, aes(variant_1, variant_2, fill = value)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = label_format(value), colour = value > threshold),
              size = 2.8, show.legend = FALSE) +
    {
      if (is.null(facet_rows) && is.null(facet_cols)) NULL else
        facet_grid(stats::as.formula(paste(
          if (is.null(facet_rows)) "." else facet_rows, "~",
          if (is.null(facet_cols)) "." else facet_cols
        )))
    } +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = limits,
                        labels = fill_labels, name = fill_name) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = GRID_INK)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_results() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}
