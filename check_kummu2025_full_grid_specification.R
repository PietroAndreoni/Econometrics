# Test the preferred test_functions.Rmd specification using KUMMU's complete
# gridded product aggregated to the Lamperti GADM 4.1 ADM1 geography.
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(fixest)
  library(ggplot2)
})
source("rmd_chunks.R")

out_dir <- "results/kummu2025_grid_gadm1"
grid_econ_path <- file.path(out_dir, "econ_kummu2025_grid_gadm41.parquet")
loader_econ_path <- file.path(out_dir, "econ_kummu2025_grid_gadm41_for_loader.parquet")
dose_econ_path <- "data/econ_processed_DOSE_V2_11-KUMMU2018-KUMMU2025-WDI-WB.parquet"
stopifnot(file.exists(grid_econ_path), file.exists(dose_econ_path))

# Keep the notebook's source label so its unchanged loader can read this file;
# provenance remains explicit in the filename and in the aggregation output.
grid_econ <- arrow::read_parquet(grid_econ_path) %>%
  mutate(econ_source = "KUMMU2025")
arrow::write_parquet(grid_econ, loader_econ_path)

nb <- load_notebook_env("test_functions.Rmd")
climate_cache <- new.env(parent = emptyenv())
nb$read_parquet_files <- function(files, columns) {
  stopifnot(length(files) == 1L)
  path <- files[[1]]
  if (grepl("data_(rr|tm)_", basename(path))) {
    key <- basename(path)
    if (!exists(key, climate_cache, inherits = FALSE)) {
      value <- arrow::open_dataset(path) %>%
        filter(
          gadm_level == "gadm1", climate_source == "era5",
          weight == "pop", weight_year == "2015"
        ) %>%
        select(all_of(columns)) %>%
        collect()
      assign(key, value, climate_cache)
    }
    return(get(key, climate_cache))
  }
  arrow::read_parquet(path, col_select = all_of(columns))
}

build_panel <- function(source, econ_path) {
  nb$build_dat(
    econ_data = source,
    econ_files_by_level = c(gadm1 = econ_path),
    climate_rr_files_by_level = c(
      gadm1 = "data/data_rr_gadm1_era5_pop-area_2000-2015.parquet"
    ),
    climate_tm_files_by_level = c(
      gadm1 = "data/data_tm_gadm1_era5_pop-area_2000-2015.parquet"
    )
  ) %>%
    mutate(zzTM = pmax(abs(zTM) - 1.5, 0))
}

message("Building full KUMMU-grid panel ...")
kummu <- build_panel("KUMMU2025", loader_econ_path)
message("Building DOSE comparison panel ...")
dose <- build_panel("DOSE_V2_11", dose_econ_path)

needed <- c(
  "dlgrp_pc_usd", "TM", "mean_TM_all", "RR", "mean_RR_all",
  "zzTM", "abs_zRRp", "zTM", "zRR"
)
eligible <- function(x) x %>% filter(if_all(all_of(needed), is.finite))
kummu <- eligible(kummu) %>% filter(year <= 2023)
dose <- eligible(dose)
stopifnot(
  !anyDuplicated(kummu[c("GID_1", "year")]),
  !anyDuplicated(dose[c("GID_1", "year")]),
  all(substr(kummu$GID_1, 1, 3) == kummu$GID_0)
)

# Import the preferred specification directly from its labelled notebook chunk.
nb$data <- kummu
nb$fit_climate_model <- function(terms, model_data = nb$data) nb$make_formula(terms)
for (expr in parse(text = extract_rmd_chunks(
  "test_functions.Rmd", "test-double deviation"
)[[1]])) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-"))) eval(expr, nb)
}
quadratic_formula <- nb$make_formula(paste(nb$base_climate, "+ zTM^2 + zRR^2"))
preferred_formula <- nb$make_formula("zTM^2 + zRR^2")
writeLines(deparse(preferred_formula), file.path(out_dir, "preferred_formula.txt"))

common_keys <- inner_join(
  dose %>% select(GID_1, year),
  kummu %>% select(GID_1, year),
  by = c("GID_1", "year")
)
samples <- list(
  KUMMU_grid_full = kummu,
  KUMMU_grid_DOSE_years = kummu %>% filter(year <= 2020),
  DOSE_full = dose,
  KUMMU_common = semi_join(kummu, common_keys, by = c("GID_1", "year")) %>% arrange(GID_1, year),
  DOSE_common = semi_join(dose, common_keys, by = c("GID_1", "year")) %>% arrange(GID_1, year)
)
stopifnot(
  nrow(samples$KUMMU_common) == nrow(samples$DOSE_common),
  all(samples$KUMMU_common$GID_1 == samples$DOSE_common$GID_1),
  all(samples$KUMMU_common$year == samples$DOSE_common$year)
)
climate_vars <- setdiff(needed, "dlgrp_pc_usd")
stopifnot(isTRUE(all.equal(
  samples$KUMMU_common[climate_vars], samples$DOSE_common[climate_vars],
  tolerance = 1e-12
)))

models <- list()
coefficient_rows <- list()
sample_rows <- list()
test_rows <- list()
for (sample_name in names(samples)) {
  model_data <- samples[[sample_name]]
  for (specification in c("preferred", "quadratic")) {
    f <- if (specification == "preferred") preferred_formula else quadratic_formula
    fit <- feols(f, data = model_data, panel.id = nb$pan_id, cluster = ~GID_1)
    model_name <- paste(sample_name, specification, sep = "__")
    models[[model_name]] <- fit
    used <- model_data[fixest::obs(fit), ]
    ct <- as.data.frame(coeftable(fit))
    ci <- confint(fit)
    coefficient_rows[[model_name]] <- data.frame(
      sample = sample_name, specification = specification,
      clustering = "GID_1", term = rownames(ct), estimate = ct[, 1],
      se = ct[, 2], p = ct[, 4], conf_low = ci[, 1], conf_high = ci[, 2]
    )
    sample_rows[[model_name]] <- data.frame(
      sample = sample_name, specification = specification, n = nobs(fit),
      regions = n_distinct(used$GID_1), countries = n_distinct(used$GID_0),
      first_year = min(used$year), last_year = max(used$year),
      wr2 = as.numeric(fitstat(fit, "wr2")[[1]])
    )
    if (specification == "preferred") {
      joint <- wald(fit, keep = "zTM|zRR", print = FALSE)
      joint_country <- wald(
        fit, keep = "zTM|zRR", vcov = ~GID_0, print = FALSE
      )
      test_rows[[sample_name]] <- bind_rows(
        data.frame(
          sample = sample_name, clustering = "GID_1",
          joint_F = joint$stat, joint_p = joint$p
        ),
        data.frame(
          sample = sample_name, clustering = "GID_0",
          joint_F = joint_country$stat, joint_p = joint_country$p
        )
      )
      country_ct <- as.data.frame(coeftable(fit, vcov = ~GID_0))
      country_ci <- confint(fit, vcov = ~GID_0)
      coefficient_rows[[paste0(model_name, "__country")]] <- data.frame(
        sample = sample_name, specification = specification,
        clustering = "GID_0", term = rownames(country_ct),
        estimate = country_ct[, 1], se = country_ct[, 2], p = country_ct[, 4],
        conf_low = country_ci[, 1], conf_high = country_ci[, 2]
      )
    }
  }
}

# The paired outcome-difference regression tests coefficient equality while
# retaining the covariance induced by identical region-year climate regressors.
delta <- samples$KUMMU_common
delta$dlgrp_pc_usd <- samples$KUMMU_common$dlgrp_pc_usd - samples$DOSE_common$dlgrp_pc_usd
difference_fit <- feols(
  preferred_formula, data = delta, panel.id = nb$pan_id, cluster = ~GID_1
)
stopifnot(max(abs(
  coef(difference_fit) -
    (coef(models$KUMMU_common__preferred) - coef(models$DOSE_common__preferred))
)) < 1e-8)
difference_ct <- as.data.frame(coeftable(difference_fit))
difference_ci <- confint(difference_fit)
difference_table <- data.frame(
  term = rownames(difference_ct), estimate = difference_ct[, 1],
  se = difference_ct[, 2], p = difference_ct[, 4],
  conf_low = difference_ci[, 1], conf_high = difference_ci[, 2]
)
models$KUMMU_minus_DOSE_common <- difference_fit

coefficients <- bind_rows(coefficient_rows)
sample_sizes <- bind_rows(sample_rows)
joint_tests <- bind_rows(test_rows)
write.csv(coefficients, file.path(out_dir, "coefficients.csv"), row.names = FALSE)
write.csv(sample_sizes, file.path(out_dir, "sample_sizes.csv"), row.names = FALSE)
write.csv(joint_tests, file.path(out_dir, "joint_tests.csv"), row.names = FALSE)
write.csv(difference_table, file.path(out_dir, "common_coefficient_differences.csv"), row.names = FALSE)
write.csv(data.frame(
  n = nrow(samples$KUMMU_common),
  growth_correlation = cor(
    samples$KUMMU_common$dlgrp_pc_usd,
    samples$DOSE_common$dlgrp_pc_usd
  ),
  growth_rmse = sqrt(mean(
    (samples$KUMMU_common$dlgrp_pc_usd - samples$DOSE_common$dlgrp_pc_usd)^2
  ))
), file.path(out_dir, "common_growth_comparison.csv"), row.names = FALSE)
saveRDS(models, file.path(out_dir, "models.rds"))

plot_data <- coefficients %>%
  filter(
    specification == "preferred", clustering == "GID_1"
  )
plot <- ggplot(plot_data, aes(x = estimate * 100, y = sample)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_errorbar(
    aes(xmin = conf_low * 100, xmax = conf_high * 100),
    orientation = "y", width = .2
  ) +
  geom_point() +
  facet_wrap(~term, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Coefficient in log-growth percentage points (95% CI)", y = NULL,
    title = "Preferred specification: full grid-aggregated KUMMU2025",
    subtitle = "GADM 4.1 ADM1; standard errors clustered by region"
  )
ggsave(file.path(out_dir, "preferred_coefficients.png"), plot, width = 10, height = 5, dpi = 170)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "model_session_info.txt"))

print(plot_data %>% select(sample, term, estimate, se, p))
print(sample_sizes)
print(joint_tests)
