# Run the preferred specification exactly from the relevant test_functions.Rmd
# chunks in a clean R process, then compare it with the GID_1-clustered entry
# exported by compare_dose_conley_se.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(ggplot2)
})

source("rmd_chunks.R")

notebook <- new.env(parent = globalenv())
notebook$params <- list(
  project_setup_script = NULL,
  climate_velocity_years = 5L,
  econ_data = "DOSE_V2_11"
)

eval_rmd_chunks(
  "test_functions.Rmd",
  labels = c(
    "shared-specification",
    "build-data",
    "test-absolute",
    "test-double deviation"
  ),
  envir = notebook
)

notebook_table <- as.data.frame(coeftable(notebook$m_pref))
notebook_interval <- as.data.frame(confint(notebook$m_pref, level = 0.95))
notebook_table$notebook_conf_low <- notebook_interval[
  rownames(notebook_table), 1
]
notebook_table$notebook_conf_high <- notebook_interval[
  rownames(notebook_table), 2
]
notebook_table$term <- rownames(notebook_table)
rownames(notebook_table) <- NULL
names(notebook_table)[1:4] <- c(
  "notebook_estimate",
  "notebook_std_error",
  "notebook_t_value",
  "notebook_p_value"
)
notebook_table <- as_tibble(notebook_table) %>%
  select(
    term,
    notebook_estimate,
    notebook_std_error,
    notebook_t_value,
    notebook_p_value,
    notebook_conf_low,
    notebook_conf_high
  )

notebook_sample_rows <- fixest::obs(notebook$m_pref)
notebook_sample <- notebook$data[notebook_sample_rows, ] %>%
  select(
    GID_0, GID_1, year, dlgrp_pc_usd,
    TM, RR, mean_TM_all, mean_RR_all,
    zTM, zRR, zzTM, abs_zRRp
  )
write.csv(
  notebook_sample,
  file.path("results", "dose_conley", "notebook_preferred_estimation_sample.csv"),
  row.names = FALSE
)

comparison_path <- file.path(
  "results", "dose_conley", "inference_comparison.csv"
)
if (!file.exists(comparison_path)) {
  stop("Run compare_dose_conley_se.R before this validation script.")
}

comparison_table <- read.csv(comparison_path, check.names = FALSE) %>%
  filter(
    model == "Preferred specification",
    vcov == "Cluster: GID_1"
  ) %>%
  transmute(
    term,
    comparison_estimate = estimate,
    comparison_std_error = std_error,
    comparison_t_value = t_value,
    comparison_p_value = p_value,
    comparison_conf_low = conf_low,
    comparison_conf_high = conf_high
  )

validation <- full_join(notebook_table, comparison_table, by = "term") %>%
  mutate(
    estimate_difference = comparison_estimate - notebook_estimate,
    std_error_difference = comparison_std_error - notebook_std_error,
    t_value_difference = comparison_t_value - notebook_t_value,
    p_value_difference = comparison_p_value - notebook_p_value,
    conf_low_difference = comparison_conf_low - notebook_conf_low,
    conf_high_difference = comparison_conf_high - notebook_conf_high
  )

comparison_sample_path <- file.path(
  "results", "dose_conley", "comparison_preferred_estimation_sample.csv"
)
if (file.exists(comparison_sample_path)) {
  comparison_sample <- read.csv(comparison_sample_path, check.names = FALSE)
  sample_validation <- full_join(
    notebook_sample %>% mutate(in_notebook = TRUE),
    comparison_sample %>% mutate(in_comparison = TRUE),
    by = c("GID_0", "GID_1", "year"),
    suffix = c("_notebook", "_comparison")
  ) %>%
    mutate(
      in_notebook = coalesce(in_notebook, FALSE),
      in_comparison = coalesce(in_comparison, FALSE)
    )
  write.csv(
    sample_validation,
    file.path("results", "dose_conley", "preferred_sample_validation.csv"),
    row.names = FALSE
  )
  cat("Rows only in notebook sample: ",
      sum(sample_validation$in_notebook & !sample_validation$in_comparison),
      "\n", sep = "")
  cat("Rows only in comparison sample: ",
      sum(!sample_validation$in_notebook & sample_validation$in_comparison),
      "\n", sep = "")
}

dir.create(
  file.path("results", "dose_conley"),
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(
  notebook_table,
  file.path("results", "dose_conley", "notebook_preferred_clean.csv"),
  row.names = FALSE
)
write.csv(
  validation,
  file.path(
    "results", "dose_conley", "notebook_vs_comparison_validation.csv"
  ),
  row.names = FALSE
)

numeric_differences <- validation %>%
  select(ends_with("_difference")) %>%
  unlist(use.names = FALSE)
if (any(!is.finite(numeric_differences)) ||
    max(abs(numeric_differences)) > 1e-12) {
  print(validation, n = Inf)
  stop("Preferred-specification results do not match within 1e-12.")
}

cat("Clean notebook preferred fit observations: ", nobs(notebook$m_pref), "\n", sep = "")
cat("Maximum absolute coefficient/SE/t/p/CI difference: ",
    format(max(abs(numeric_differences)), scientific = TRUE), "\n", sep = "")
print(validation, n = Inf)
