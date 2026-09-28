# Compare country-clustered, region-clustered, and Conley standard errors for
# the preferred and quadratic DOSE specifications in test_functions.Rmd.
#
# The script deliberately imports the notebook's data builder and formula
# helpers so that the regressions use the same variables, fixed effects, and
# panel definition as the main analysis. Only the variance estimator changes.

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(ggplot2)
  library(sf)
  library(tidyr)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
notebook_path <- file.path(project_root, "test_functions.Rmd")
chunk_helper_path <- file.path(project_root, "rmd_chunks.R")
gadm_dir <- file.path(
  project_root, "results", "kummu2025_robustness", "sources", "gadm41"
)
output_dir <- file.path(project_root, "results", "dose_conley")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(notebook_path, chunk_helper_path, gadm_dir)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required inputs: ", paste(missing_paths, collapse = ", "))
}

source(chunk_helper_path)
notebook <- load_notebook_env(notebook_path, econ_data = "DOSE_V2_11")

# Keep the preferred transformation exactly as written in test_functions.Rmd.
dose_data <- notebook$build_dat(econ_data = "DOSE_V2_11") %>% 
  filter(year>=1990) |> 
  mutate(
    zzTM = pmax(
      1.5,
      abs((TM - mean_TM_lag_30) / sd_TM_lag_30)
    ) - 1.5
  )

preferred_terms <- paste(notebook$base_climate, "+ zzTM^2 + abs_zRRp")
quadratic_terms <- paste(notebook$base_climate, "+ zTM^2 + zRR^2")
model_terms <- c(
  "Preferred specification" = preferred_terms,
  "Quadratic specification" = quadratic_terms
)

# GADM 4.1 ADM1 polygons were previously downloaded as one GeoJSON per country.
# Reading only countries in the DOSE estimation candidate data avoids loading
# irrelevant geometry. Exact GID joins are audited below before estimation.
candidate_countries <- sort(unique(na.omit(dose_data$GID_0)))
gadm_files <- file.path(gadm_dir, paste0(candidate_countries, ".json"))
missing_gadm_files <- gadm_files[!file.exists(gadm_files)]
if (length(missing_gadm_files)) {
  stop(
    "Missing GADM 4.1 country files for: ",
    paste(tools::file_path_sans_ext(basename(missing_gadm_files)), collapse = ", ")
  )
}

gadm_adm1 <- bind_rows(lapply(gadm_files, function(path) {
  sf::st_read(path, quiet = TRUE) %>%
    select(GID_0, GID_1, NAME_1, geometry)
}))

gadm_id_audit <- gadm_adm1 %>%
  st_drop_geometry() %>%
  transmute(
    GID_0,
    GID_1,
    NAME_1,
    valid_country_prefix = !is.na(GID_1) &
      sub("\\..*$", "", GID_1) == GID_0
  )
write.csv(
  gadm_id_audit %>% filter(!valid_country_prefix),
  file.path(output_dir, "gadm41_geometry_id_issues.csv"),
  row.names = FALSE
)
gadm_adm1 <- gadm_adm1[gadm_id_audit$valid_country_prefix, ]

if (anyDuplicated(gadm_adm1$GID_1)) {
  stop("Valid GADM 4.1 geometry contains duplicate GID_1 identifiers.")
}

# Compute polygon centroids in a global equal-area projection, then return the
# coordinates to WGS84 for great-circle distances in vcov_conley().
gadm_centroids <- gadm_adm1 %>%
  st_make_valid() %>%
  st_transform(6933) %>%
  st_centroid() %>%
  st_transform(4326)

centroid_xy <- st_coordinates(gadm_centroids)
centroids <- gadm_centroids %>%
  st_drop_geometry() %>%
  transmute(
    GID_0,
    GID_1,
    NAME_1,
    longitude = centroid_xy[, "X"],
    latitude = centroid_xy[, "Y"]
  )

if (any(!is.finite(centroids$longitude)) ||
    any(!is.finite(centroids$latitude)) ||
    any(abs(centroids$longitude) > 180) ||
    any(abs(centroids$latitude) > 90)) {
  stop("Invalid longitude or latitude produced by the centroid calculation.")
}

candidate_ids <- dose_data %>%
  filter(!is.na(dlgrp_pc_usd)) %>%
  distinct(GID_0, GID_1)
centroid_audit <- candidate_ids %>%
  left_join(centroids, by = c("GID_0", "GID_1")) %>%
  mutate(
    country_prefix_valid = sub("\\..*$", "", GID_1) == GID_0,
    centroid_matched = is.finite(longitude) & is.finite(latitude)
  )

if (any(!centroid_audit$country_prefix_valid)) {
  stop("DOSE contains a GID_1 whose prefix does not agree with GID_0.")
}
if (any(!centroid_audit$centroid_matched)) {
  missing_ids <- centroid_audit$GID_1[!centroid_audit$centroid_matched]
  stop(
    "DOSE estimation candidates without an exact GADM 4.1 centroid match: ",
    paste(missing_ids, collapse = ", ")
  )
}

write.csv(
  centroids %>% semi_join(candidate_ids, by = c("GID_0", "GID_1")),
  file.path(output_dir, "gadm41_centroids.csv"),
  row.names = FALSE
)
write.csv(
  centroid_audit,
  file.path(output_dir, "gadm41_centroid_audit.csv"),
  row.names = FALSE
)

dose_spatial <- dose_data %>%
  left_join(
    centroids %>% select(GID_0, GID_1, longitude, latitude),
    by = c("GID_0", "GID_1")
  )

# Estimate without committing to a VCOV first. The fits are reused for the
# residual-distance diagnostic and all three inference methods.
fit_model <- function(terms) {
  feols(
    notebook$make_formula(terms),
    data = dose_spatial,
    panel.id = notebook$pan_id,
    vcov = "iid"
  )
}
models <- lapply(model_terms, fit_model)

preferred_sample_rows <- fixest::obs(models[["Preferred specification"]])
write.csv(
  dose_spatial[preferred_sample_rows, ] %>%
    select(
      GID_0, GID_1, year, dlgrp_pc_usd,
      TM, RR, mean_TM_all, mean_RR_all,
      zTM, zRR, zzTM, abs_zRRp
    ),
  file.path(output_dir, "comparison_preferred_estimation_sample.csv"),
  row.names = FALSE
)

model_sample <- bind_rows(lapply(names(models), function(label) {
  model <- models[[label]]
  rows <- fixest::obs(model)
  model_data <- dose_spatial[rows, , drop = FALSE]
  tibble(
    observations = nobs(model),
    model = label,
    regions = n_distinct(model_data$GID_1),
    countries = n_distinct(model_data$GID_0),
    first_year = min(model_data$year),
    last_year = max(model_data$year)
  )
}))
write.csv(
  model_sample,
  file.path(output_dir, "model_sample_summary.csv"),
  row.names = FALSE
)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi / 180
  dlon <- (lon2 - lon1) * to_rad
  dlat <- (lat2 - lat1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# Create a common region-pair distance lookup. Pairwise standardized residual
# products are then averaged within year and distance bin. This avoids treating
# repeated annual observations at the same centroid as distinct spatial units.
diagnostic_regions <- sort(unique(unlist(lapply(models, function(model) {
  dose_spatial$GID_1[fixest::obs(model)]
}))))
diagnostic_centroids <- centroids %>%
  filter(GID_1 %in% diagnostic_regions) %>%
  arrange(GID_1)
pair_index <- which(
  upper.tri(matrix(FALSE, nrow(diagnostic_centroids), nrow(diagnostic_centroids))),
  arr.ind = TRUE
)
region_pairs <- tibble(
  region_i = diagnostic_centroids$GID_1[pair_index[, 1]],
  region_j = diagnostic_centroids$GID_1[pair_index[, 2]],
  distance_km = haversine_km(
    diagnostic_centroids$longitude[pair_index[, 1]],
    diagnostic_centroids$latitude[pair_index[, 1]],
    diagnostic_centroids$longitude[pair_index[, 2]],
    diagnostic_centroids$latitude[pair_index[, 2]]
  )
)

DISTANCE_BIN_KM <- 100
DIAGNOSTIC_MAX_KM <- 3000
NEAR_ZERO_CORRELATION <- 0.02
CONSECUTIVE_ZERO_BINS <- 3L

residual_distance_by_model <- function(model, label) {
  rows <- fixest::obs(model)
  residual_data <- dose_spatial[rows, c("GID_1", "year"), drop = FALSE] %>%
    mutate(residual = as.numeric(residuals(model))) %>%
    group_by(year) %>%
    mutate(
      residual_sd = sd(residual),
      standardized_residual = if_else(
        is.finite(residual_sd) & residual_sd > 0,
        (residual - mean(residual)) / residual_sd,
        NA_real_
      )
    ) %>%
    ungroup()

  annual_bins <- bind_rows(lapply(sort(unique(residual_data$year)), function(y) {
    year_residuals <- residual_data %>%
      filter(year == y) %>%
      select(GID_1, standardized_residual)
    pair_products <- region_pairs %>%
      filter(distance_km <= DIAGNOSTIC_MAX_KM) %>%
      left_join(year_residuals, by = c("region_i" = "GID_1")) %>%
      rename(residual_i = standardized_residual) %>%
      left_join(year_residuals, by = c("region_j" = "GID_1")) %>%
      rename(residual_j = standardized_residual) %>%
      filter(is.finite(residual_i), is.finite(residual_j)) %>%
      mutate(
        distance_bin = floor(distance_km / DISTANCE_BIN_KM) * DISTANCE_BIN_KM,
        residual_product = residual_i * residual_j
      ) %>%
      group_by(distance_bin) %>%
      summarise(
        spatial_correlation = mean(residual_product),
        region_pairs = n(),
        .groups = "drop"
      ) %>%
      mutate(year = y)
    pair_products
  }))

  annual_bins %>%
    group_by(distance_bin) %>%
    summarise(
      distance_midpoint_km = first(distance_bin) + DISTANCE_BIN_KM / 2,
      mean_spatial_correlation = weighted.mean(
        spatial_correlation, region_pairs, na.rm = TRUE
      ),
      median_annual_correlation = median(spatial_correlation, na.rm = TRUE),
      annual_sd = sd(spatial_correlation, na.rm = TRUE),
      years = sum(is.finite(spatial_correlation)),
      region_year_pairs = sum(region_pairs),
      se_across_years = annual_sd / sqrt(years),
      ci_low = mean_spatial_correlation - qt(0.975, pmax(years - 1, 1)) *
        se_across_years,
      ci_high = mean_spatial_correlation + qt(0.975, pmax(years - 1, 1)) *
        se_across_years,
      .groups = "drop"
    ) %>%
    mutate(model = label, .before = 1)
}

residual_distance <- bind_rows(Map(
  residual_distance_by_model,
  models,
  names(models)
))

choose_cutoff <- function(diagnostic) {
  diagnostic <- diagnostic %>%
    arrange(distance_bin) %>%
    mutate(
      near_zero = abs(mean_spatial_correlation) <= NEAR_ZERO_CORRELATION &
        ci_low <= 0 & ci_high >= 0
    )
  run_ok <- vapply(seq_len(nrow(diagnostic)), function(i) {
    end <- i + CONSECUTIVE_ZERO_BINS - 1L
    end <= nrow(diagnostic) && all(diagnostic$near_zero[i:end])
  }, logical(1))
  eligible <- which(run_ok & diagnostic$distance_bin > 0)
  if (!length(eligible)) {
    warning(
      "No sustained near-zero residual-correlation run found; using ",
      DIAGNOSTIC_MAX_KM, " km."
    )
    return(DIAGNOSTIC_MAX_KM)
  }
  diagnostic$distance_bin[min(eligible)]
}

cutoff_by_model <- residual_distance %>%
  group_split(model) %>%
  lapply(function(x) {
    tibble(model = first(x$model), selected_cutoff_km = choose_cutoff(x))
  }) %>%
  bind_rows()

# A common cutoff makes the two models' Conley standard errors comparable. Use
# the more conservative model-specific choice when the diagnostics differ.
conley_cutoff_km <- max(cutoff_by_model$selected_cutoff_km)
cutoff_selection <- cutoff_by_model %>%
  mutate(
    common_conley_cutoff_km = conley_cutoff_km,
    distance_bin_km = DISTANCE_BIN_KM,
    near_zero_threshold = NEAR_ZERO_CORRELATION,
    consecutive_bins_required = CONSECUTIVE_ZERO_BINS
  )

write.csv(
  residual_distance,
  file.path(output_dir, "residual_spatial_correlation_by_distance.csv"),
  row.names = FALSE
)
write.csv(
  cutoff_selection,
  file.path(output_dir, "conley_cutoff_selection.csv"),
  row.names = FALSE
)

residual_plot <- ggplot(
  residual_distance,
  aes(distance_midpoint_km, mean_spatial_correlation)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.18) +
  geom_line(linewidth = 0.65, color = "#176B87") +
  geom_point(size = 1.1, color = "#176B87") +
  geom_vline(
    xintercept = conley_cutoff_km,
    linetype = "dashed",
    color = "#B03A2E"
  ) +
  facet_wrap(~model, ncol = 1) +
  labs(
    title = "DOSE residual spatial correlation by GADM 4.1 centroid distance",
    subtitle = paste0(
      "Points are pair-weighted annual means; bands are 95% intervals across years. ",
      "Dashed line: common Conley cutoff = ", conley_cutoff_km, " km."
    ),
    x = "Distance between ADM1 centroids (km)",
    y = "Mean standardized residual product"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(output_dir, "residual_spatial_correlation_by_distance.png"),
  residual_plot,
  width = 9,
  height = 8,
  dpi = 300
)

extract_inference <- function(fitted_model, model_label, method_label, vcov_spec) {
  table <- as.data.frame(coeftable(fitted_model, vcov = vcov_spec))
  interval <- as.data.frame(confint(fitted_model, vcov = vcov_spec, level = 0.95))
  table$conf_low <- interval[rownames(table), 1]
  table$conf_high <- interval[rownames(table), 2]
  table$term <- rownames(table)
  rownames(table) <- NULL
  names(table)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  as_tibble(table) %>%
    transmute(
      observations = nobs(fitted_model),
      model = model_label,
      vcov = method_label,
      term,
      estimate,
      std_error,
      t_value,
      p_value,
      conf_low,
      conf_high
    )
}

# Inspect the uncorrected Conley matrices explicitly. fixest's PSD repair is
# retained for reported inference, but the correction is recorded rather than
# left only as a console warning.
conley_vcovs <- lapply(models, function(model) {
  raw <- suppressWarnings(vcov_conley(
    model,
    lat = ~latitude,
    lon = ~longitude,
    cutoff = conley_cutoff_km,
    distance = "spherical",
    vcov_fix = FALSE
  ))
  fixed <- suppressWarnings(vcov_conley(
    model,
    lat = ~latitude,
    lon = ~longitude,
    cutoff = conley_cutoff_km,
    distance = "spherical",
    vcov_fix = TRUE
  ))
  list(raw = raw, fixed = fixed)
})

conley_vcov_diagnostics <- bind_rows(lapply(names(conley_vcovs), function(label) {
  raw_eigenvalues <- eigen(
    conley_vcovs[[label]]$raw,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  fixed_eigenvalues <- eigen(
    conley_vcovs[[label]]$fixed,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  tibble(
    model = label,
    cutoff_km = conley_cutoff_km,
    raw_minimum_eigenvalue = min(raw_eigenvalues),
    fixed_minimum_eigenvalue = min(fixed_eigenvalues),
    psd_correction_applied = min(raw_eigenvalues) < 0
  )
}))
write.csv(
  conley_vcov_diagnostics,
  file.path(output_dir, "conley_vcov_diagnostics.csv"),
  row.names = FALSE
)

inference <- bind_rows(lapply(names(models), function(label) {
  model <- models[[label]]
  bind_rows(
    extract_inference(model, label, "Cluster: GID_0", ~GID_0),
    extract_inference(model, label, "Cluster: GID_1", ~GID_1),
    extract_inference(
      model,
      label,
      paste0("Conley: ", conley_cutoff_km, " km"),
      conley_vcovs[[label]]$fixed
    )
  )
}))

write.csv(
  inference,
  file.path(output_dir, "inference_comparison.csv"),
  row.names = FALSE
)

inference_se_comparison <- inference %>%
  group_by(model, term, estimate, observations) %>%
  summarise(
    se_cluster_gid0 = std_error[vcov == "Cluster: GID_0"],
    se_cluster_gid1 = std_error[vcov == "Cluster: GID_1"],
    se_conley = std_error[vcov == paste0("Conley: ", conley_cutoff_km, " km")],
    p_cluster_gid0 = p_value[vcov == "Cluster: GID_0"],
    p_cluster_gid1 = p_value[vcov == "Cluster: GID_1"],
    p_conley = p_value[vcov == paste0("Conley: ", conley_cutoff_km, " km")],
    .groups = "drop"
  ) %>%
  mutate(
    conley_to_gid0_se_ratio = se_conley / se_cluster_gid0,
    conley_to_gid1_se_ratio = se_conley / se_cluster_gid1
  )
write.csv(
  inference_se_comparison,
  file.path(output_dir, "standard_error_comparison.csv"),
  row.names = FALSE
)

inference_plot <- inference %>%
  mutate(vcov = factor(
    vcov,
    levels = c(
      "Cluster: GID_0",
      "Cluster: GID_1",
      paste0("Conley: ", conley_cutoff_km, " km")
    )
  )) %>%
  ggplot(aes(estimate, term, color = vcov)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    position = position_dodge(width = 0.55),
    width = 0,
    orientation = "y"
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 1.6) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  labs(
    title = "DOSE coefficient inference under alternative spatial VCOVs",
    subtitle = "Point estimates and samples are fixed within specification; bars are 95% intervals.",
    x = "Coefficient estimate",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "inference_comparison.png"),
  inference_plot,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Wrote DOSE Conley analysis to: ", output_dir, "\n", sep = "")
cat("Common Conley cutoff: ", conley_cutoff_km, " km\n", sep = "")
print(model_sample)
print(cutoff_selection)
