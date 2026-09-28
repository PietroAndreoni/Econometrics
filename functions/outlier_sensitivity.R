# How a climate-response model moves as ranked subnational outliers (from
# rank_subnational_outliers()) are dropped from its estimation sample.

# `estimation_sample` should already be restricted to the rows the untrimmed
# model uses, so each removal count is a number of region-years that actually
# enter the regression. Outliers outside the sample are ranked but skipped.
# Returns the coefficient table per removal count, the sample sizes, and the
# in-sample outlier ranking.
outlier_removal_sensitivity <- function(
    estimation_sample,
    ranked_outliers,
    terms,
    removal_counts,
    conf_level = 0.95,
    fit_fun = fit_panel_model
) {
  in_sample <- ranked_outliers %>%
    semi_join(
      estimation_sample %>% distinct(GID_0, GID_1, year),
      by = c("GID_0", "GID_1", "year")
    ) %>%
    arrange(outlier_rank) %>%
    mutate(removal_rank = row_number())

  usable <- removal_counts[removal_counts <= nrow(in_sample)]
  if (length(usable) < length(removal_counts)) {
    warning(
      "Only ", nrow(in_sample), " ranked outliers are in the estimation ",
      "sample; larger removal counts were skipped."
    )
  }
  usable <- sort(unique(c(0L, as.integer(usable))))

  coefficients <- bind_rows(lapply(usable, function(n_removed) {
    removed <- in_sample[seq_len(n_removed), ]
    model_data <- anti_join(
      estimation_sample,
      removed %>% select(GID_0, GID_1, year),
      by = c("GID_0", "GID_1", "year")
    )
    fit <- fit_fun(terms, model_data = model_data)

    bind_cols(
      tibble::tibble(
        n_outliers_removed = as.integer(n_removed),
        rows_dropped = nrow(estimation_sample) - nrow(model_data),
        min_outlier_metric_removed = if (n_removed > 0L) {
          min(removed$outlier_metric)
        } else {
          NA_real_
        },
        observations = as.integer(stats::nobs(fit)),
        regions = n_distinct(model_data$GID_1),
        countries = n_distinct(model_data$GID_0)
      ),
      tidy_fixest(fit, conf_level = conf_level) %>%
        select(term, estimate, std_error, statistic, p_value, conf_low, conf_high)
    )
  })) %>%
    mutate(term = factor(term, levels = unique(term)))

  list(
    coefficients = coefficients,
    samples = coefficients %>%
      distinct(
        n_outliers_removed, rows_dropped, min_outlier_metric_removed,
        observations, regions, countries
      ),
    ranked_in_sample = in_sample,
    removal_counts = usable
  )
}

# Signed-bin coefficients per removal count, with the omitted central bin added
# at zero so every curve is shown relative to it.
signed_bin_paths <- function(sensitivity,
                             pattern = "^(TM|RR)_bin_signed::") {
  bins <- sensitivity$coefficients %>%
    filter(grepl(pattern, as.character(term))) %>%
    mutate(
      weather = if_else(
        startsWith(as.character(term), "TM_"),
        "Temperature",
        "Precipitation"
      ),
      bin = suppressWarnings(as.numeric(sub("^.*::", "", as.character(term))))
    )
  if (!nrow(bins) || anyNA(bins$bin)) {
    stop("Could not recover the signed temperature and precipitation bins.")
  }

  reference <- tidyr::expand_grid(
    n_outliers_removed = sensitivity$removal_counts,
    weather = unique(bins$weather)
  ) %>%
    mutate(bin = 0, estimate = 0, std_error = 0, conf_low = 0, conf_high = 0)

  bind_rows(bins, reference) %>%
    arrange(weather, n_outliers_removed, bin)
}

plot_signed_bin_paths <- function(paths, removal_counts, title,
                                  subtitle = NULL, caption = NULL) {
  ggplot(
    paths,
    aes(x = bin, y = estimate, colour = n_outliers_removed,
        group = n_outliers_removed)
  ) +
    geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3,
               linetype = "dashed") +
    geom_ribbon(
      aes(ymin = conf_low, ymax = conf_high, fill = n_outliers_removed),
      alpha = 0.08,
      colour = NA
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.7) +
    facet_wrap(~weather, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(paths$bin))) +
    scale_colour_viridis_c(
      option = "C", direction = -1, breaks = removal_counts,
      labels = format(removal_counts, big.mark = ",")
    ) +
    scale_fill_viridis_c(
      option = "C", direction = -1, breaks = removal_counts,
      labels = format(removal_counts, big.mark = ","), guide = "none"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Signed anomaly bin (standard deviations)",
      y = "Estimated effect relative to the central bin",
      colour = "Regional outliers removed",
      caption = caption
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", colour = "#0b0b0b"),
      plot.subtitle = element_text(colour = "#52514e"),
      plot.caption = element_text(colour = "#52514e", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "#0b0b0b"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      axis.title = element_text(colour = "#52514e")
    )
}
