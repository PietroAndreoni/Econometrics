# Full national-vs-subnational growth report for one or more source pairs:
# the comparison tables, the ranked subnational outliers and, optionally, how
# the signed-bin climate response moves as those outliers are dropped.
#
#   national_subnational_growth_report("DOSE", "PWT")
#   national_subnational_growth_report(c("DOSE", "KUMMU"), c("PWT", "WB"))
#
# Sources are those of read_gdp_levels(). Every subnational source is compared
# with every national one; results go to
# <output_root>/<subnational>_vs_<national>/. Returns, invisibly, a list with
# one element per pair holding `comparison` and `sensitivity` (NULL when
# run_sensitivity = FALSE).

# The joint signed-deviation-bin model of test_functions.Rmd. The central -0.5
# to 0.5 standard-deviation interval is the omitted category.
SIGNED_BIN_MODEL_TERMS <- paste(
  "TM + TM:mean_TM_all + RR + RR:mean_RR_all",
  "+ i(TM_bin_signed, ref = 0)",
  "+ i(RR_bin_signed, ref = 0)"
)

national_subnational_growth_report <- function(
    subnational = c("DOSE", "KUMMU"),
    national = c("PWT", "WB"),
    output_root = file.path("results", "national_subnational_growth"),
    flag_percentile = 0.99,
    top_n = 50L,
    top_n_detail = 20L,
    outlier_metric = "abs_additive_contribution_pp",
    run_sensitivity = TRUE,
    removal_counts = seq(1000L, 5000L, by = 1000L),
    sensitivity_terms = SIGNED_BIN_MODEL_TERMS,
    conf_level = 0.95,
    config = panel_config(),
    verbose = TRUE
) {
  say <- function(...) if (verbose) cat(..., sep = "")
  levels_cache <- list()
  levels_for <- function(source) {
    if (is.null(levels_cache[[source]])) {
      levels_cache[[source]] <<- read_gdp_levels(source)
    }
    levels_cache[[source]]
  }

  pairs <- expand.grid(
    subnational = subnational,
    national = national,
    stringsAsFactors = FALSE
  )
  results <- list()

  for (i in seq_len(nrow(pairs))) {
    sub_source <- pairs$subnational[[i]]
    nat_source <- pairs$national[[i]]
    label <- paste0(sub_source, "_vs_", nat_source)
    output_dir <- file.path(output_root, label)
    say("\n==== ", label, " ====\n")

    comparison <- compare_national_subnational(
      levels_for(sub_source),
      levels_for(nat_source),
      flag_percentile = flag_percentile,
      top_n = top_n,
      outlier_metric = outlier_metric
    )
    write_national_subnational_comparison(
      comparison, output_dir, top_n = top_n, top_n_detail = top_n_detail
    )

    gc <- comparison$growth_comparison
    say(
      "Compared ", sum(gc$eligible_for_comparison),
      " eligible country-year growth pairs across ",
      n_distinct(gc$iso3[gc$eligible_for_comparison]), " countries.\n",
      "High-discrepancy threshold (", 100 * flag_percentile, "th percentile): ",
      round(comparison$flag_threshold_pp, 3), " percentage points.\n"
    )
    if (verbose) {
      cat("\nLargest subnational outliers (", outlier_metric, "):\n", sep = "")
      print(head(comparison$ranked_outliers, 20L), n = 20L)
    }

    sensitivity <- NULL
    if (run_sensitivity) {
      sensitivity <- .report_outlier_sensitivity(
        comparison, sub_source, nat_source, output_dir, outlier_metric,
        removal_counts, sensitivity_terms, conf_level, config
      )
      if (verbose) {
        cat("\nOutlier-removal sample sizes:\n")
        print(sensitivity$samples, n = Inf)
      }
    }
    say("Results written to: ", normalizePath(output_dir), "\n")

    results[[label]] <- list(comparison = comparison, sensitivity = sensitivity)
  }
  invisible(results)
}

# Re-estimates the model on the subnational source's build_dat() panel as
# ranked outliers are dropped, and writes the paths, sample sizes and plot.
.report_outlier_sensitivity <- function(comparison, sub_source, nat_source,
                                        output_dir, outlier_metric,
                                        removal_counts, terms, conf_level,
                                        config) {
  # Restrict once to the exact sample of the untrimmed model, so the removal
  # counts count region-years that would otherwise enter that model.
  estimation_sample <- build_dat(econ_data = sub_source, config = config) %>%
    filter(!is.na(dlgrp_pc_usd))
  baseline_fit <- fit_panel_model(terms, model_data = estimation_sample)
  estimation_sample <- estimation_sample %>% slice(fixest::obs(baseline_fit))

  sensitivity <- outlier_removal_sensitivity(
    estimation_sample,
    comparison$ranked_outliers,
    terms = terms,
    removal_counts = removal_counts,
    conf_level = conf_level
  )
  paths <- signed_bin_paths(sensitivity)

  utils::write.csv(
    paths,
    file.path(output_dir, "signed_bin_outlier_removal_coefficients.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    sensitivity$samples,
    file.path(output_dir, "outlier_removal_sample_sizes.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    sensitivity$ranked_in_sample %>%
      filter(removal_rank <= max(sensitivity$removal_counts)),
    file.path(output_dir, "ranked_region_year_outliers_in_sample.csv"),
    row.names = FALSE
  )
  ggsave(
    file.path(output_dir, "signed_bin_outlier_removal_paths.png"),
    plot_signed_bin_paths(
      paths,
      sensitivity$removal_counts,
      title = paste0(
        "Signed-bin response after removing ", sub_source, "-", nat_source,
        " regional outliers"
      ),
      subtitle = paste0(
        "Each colour removes additional individual region-years; ",
        100 * conf_level,
        "% confidence ribbons use standard errors clustered by region"
      ),
      caption = paste0(
        "Sample: ", sub_source, " subnational panel from build_dat(); ",
        "outliers ranked by ", outlier_metric, ". Observations fall from ",
        format(max(sensitivity$samples$observations), big.mark = ","), " to ",
        format(min(sensitivity$samples$observations), big.mark = ","), "."
      )
    ),
    width = 11, height = 6.5, dpi = 300, bg = "#fcfcfb"
  )
  sensitivity
}
