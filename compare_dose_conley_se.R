# Compare country-clustered, region-clustered and Conley standard errors for
# the preferred and quadratic DOSE specifications of main_analysis.Rmd.
#
# Point estimates and samples are fixed within a specification; only the
# variance estimator changes. The Conley cutoff is chosen from the data: the
# first distance at which the spatial correlation of year-standardized
# residuals stays near zero for several consecutive bins (functions/
# spatial_inference.R). A common cutoff, the larger of the two models' choices,
# keeps the two specifications comparable.
#
# Region coordinates: GADM 4.1 ADM1 polygon centroids from
# data/kummu2025/gadm41 (the same polygons the Kummu panel uses), computed in an
# equal-area projection.
#
# Specification: main_spec() (functions/analysis_spec.R), i.e. main_analysis.Rmd.
# Output: results/dose_conley/.

source("load_functions.R")

SPEC <- main_spec()
DISTANCE_BIN_KM <- 100
DIAGNOSTIC_MAX_KM <- 3000
NEAR_ZERO_CORRELATION <- 0.02
CONSECUTIVE_ZERO_BINS <- 3L
OUT <- file.path("results", "dose_conley")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# Panel and specifications ---------------------------------------------------------

dose <- build_main_dat("DOSE", SPEC) %>%
  mutate(zzTM = pmax(1.5, abs((TM - mean_TM_lag_30) / sd_TM_lag_30)) - 1.5)

MODEL_TERMS <- c(
  "Preferred specification" = paste(SPEC$base_climate, "+ zzTM^2 + abs_zRRp"),
  "Quadratic specification" = paste(SPEC$base_climate, "+ zTM^2 + zRR^2")
)

# Region centroids -------------------------------------------------------------------

candidate_ids <- dose %>%
  filter(!is.na(.data[[SPEC$outcome]])) %>%
  distinct(GID_0, GID_1)
centroids <- read_gadm1_polygons(
  file.path(kummu2025_dir(), "gadm41"),
  keep_ids = candidate_ids$GID_1
) %>%
  polygon_centroids()

centroid_audit <- candidate_ids %>%
  left_join(centroids, by = c("GID_0", "GID_1")) %>%
  mutate(centroid_matched = is.finite(longitude) & is.finite(latitude))
write_out(centroid_audit, "gadm41_centroid_audit")
write_out(centroids, "gadm41_centroids")
unmatched <- centroid_audit %>% filter(!centroid_matched)
if (nrow(unmatched)) {
  warning(nrow(unmatched), " DOSE regions have no GADM 4.1 polygon and are ",
          "excluded from every model: ", paste(unmatched$GID_1, collapse = ", "))
}

dose_spatial <- dose %>%
  inner_join(centroids %>% select(GID_0, GID_1, longitude, latitude),
             by = c("GID_0", "GID_1"))

# Models (the variance estimator is chosen at inference time) --------------------------

models <- lapply(MODEL_TERMS, fit_main, model_data = dose_spatial, spec = SPEC)

model_sample <- bind_rows(lapply(names(models), function(label) {
  used <- dose_spatial[fixest::obs(models[[label]]), ]
  tibble::tibble(
    model = label, observations = stats::nobs(models[[label]]),
    regions = n_distinct(used$GID_1), countries = n_distinct(used$GID_0),
    first_year = min(used$year), last_year = max(used$year)
  )
}))
write_out(model_sample, "model_sample_summary")
write_out(
  dose_spatial[fixest::obs(models[["Preferred specification"]]), ] %>%
    select(GID_0, GID_1, year, all_of(SPEC$outcome), TM, RR, mean_TM_all,
           mean_RR_all, zTM, zRR, zzTM, abs_zRRp, longitude, latitude),
  "preferred_estimation_sample"
)

# Residual spatial correlation and the Conley cutoff ------------------------------------

pairs <- region_pair_distances(
  centroids %>% semi_join(candidate_ids, by = c("GID_0", "GID_1")),
  max_km = DIAGNOSTIC_MAX_KM
)
residual_distance <- bind_rows(lapply(names(models), function(label) {
  residual_spatial_correlation(models[[label]], dose_spatial, pairs,
                               bin_km = DISTANCE_BIN_KM,
                               max_km = DIAGNOSTIC_MAX_KM) %>%
    mutate(model = label, .before = 1)
}))
cutoff_by_model <- residual_distance %>%
  group_by(model) %>%
  group_modify(~ tibble::tibble(selected_cutoff_km = choose_conley_cutoff(
    .x, near_zero = NEAR_ZERO_CORRELATION,
    consecutive_bins = CONSECUTIVE_ZERO_BINS, fallback_km = DIAGNOSTIC_MAX_KM
  ))) %>%
  ungroup()
conley_cutoff_km <- max(cutoff_by_model$selected_cutoff_km)
cutoff_selection <- cutoff_by_model %>%
  mutate(common_conley_cutoff_km = conley_cutoff_km,
         distance_bin_km = DISTANCE_BIN_KM,
         near_zero_threshold = NEAR_ZERO_CORRELATION,
         consecutive_bins_required = CONSECUTIVE_ZERO_BINS)
write_out(residual_distance, "residual_spatial_correlation_by_distance")
write_out(cutoff_selection, "conley_cutoff_selection")

# Inference under the three variance estimators ------------------------------------------

conley <- lapply(models, conley_vcov_check, cutoff_km = conley_cutoff_km)
write_out(
  bind_rows(lapply(names(conley), function(label) {
    mutate(conley[[label]]$diagnostics, model = label, .before = 1)
  })),
  "conley_vcov_diagnostics"
)

conley_label <- paste0("Conley: ", conley_cutoff_km, " km")
VCOV_LEVELS <- c("Cluster: GID_0", "Cluster: GID_1", conley_label)
inference <- bind_rows(lapply(names(models), function(label) {
  m <- models[[label]]
  bind_rows(
    tidy_fixest(m, label, vcov = ~GID_0) %>% mutate(vcov = VCOV_LEVELS[1]),
    tidy_fixest(m, label, vcov = ~GID_1) %>% mutate(vcov = VCOV_LEVELS[2]),
    tidy_fixest(m, label, vcov = conley[[label]]$fixed) %>%
      mutate(vcov = VCOV_LEVELS[3])
  )
}))
write_out(inference, "inference_comparison")

se_comparison <- inference %>%
  select(model, term, estimate, observations, vcov, std_error, p_value) %>%
  mutate(vcov = c("gid0", "gid1", "conley")[match(vcov, VCOV_LEVELS)]) %>%
  tidyr::pivot_wider(names_from = vcov, values_from = c(std_error, p_value)) %>%
  mutate(conley_to_gid0_se_ratio = std_error_conley / std_error_gid0,
         conley_to_gid1_se_ratio = std_error_conley / std_error_gid1)
write_out(se_comparison, "standard_error_comparison")

# Figures ------------------------------------------------------------------------------------

p_residual <- ggplot(residual_distance,
                     aes(distance_midpoint_km, mean_spatial_correlation)) +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = GRID_PALETTE[1],
              alpha = 0.18) +
  geom_line(colour = GRID_PALETTE[1], linewidth = 0.65) +
  geom_point(colour = GRID_PALETTE[1], size = 1.1) +
  geom_vline(xintercept = conley_cutoff_km, linetype = "dashed",
             colour = GRID_PALETTE[2]) +
  facet_wrap(~model, ncol = 1) +
  labs(
    title = "DOSE residual spatial correlation by GADM 4.1 centroid distance",
    subtitle = paste0("Pair-weighted annual means; bands are 95% intervals across ",
                      "years. Dashed line: common Conley cutoff = ",
                      conley_cutoff_km, " km."),
    x = "Distance between ADM1 centroids (km)",
    y = "Mean standardized residual product"
  ) +
  theme_results()
ggsave(file.path(OUT, "residual_spatial_correlation_by_distance.png"), p_residual,
       width = 9, height = 8, dpi = 300)

vcov_colours <- stats::setNames(GRID_PALETTE[1:3], VCOV_LEVELS)
p_inference <- inference %>%
  mutate(vcov = factor(vcov, levels = VCOV_LEVELS)) %>%
  ggplot(aes(estimate, term, colour = vcov, shape = vcov)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high),
                position = position_dodge(width = 0.55), width = 0,
                orientation = "y") +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_colour_manual(values = vcov_colours, name = NULL) +
  scale_shape_manual(values = stats::setNames(GRID_SHAPES[1:3], VCOV_LEVELS),
                     name = NULL) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  labs(
    title = "DOSE coefficient inference under alternative variance estimators",
    subtitle = "Estimates and samples are fixed within a specification; bars are 95% intervals.",
    x = "Coefficient estimate", y = NULL
  ) +
  theme_results() +
  theme(panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
        panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "inference_comparison.png"), p_inference,
       width = 10, height = 8, dpi = 300)

# Console summary --------------------------------------------------------------------------

options(width = 180)
print(model_sample)
print(cutoff_selection)
cat("\nCommon Conley cutoff: ", conley_cutoff_km, " km\n", sep = "")
print(se_comparison %>% select(model, term, estimate, starts_with("std_error"),
                               starts_with("p_value"), ends_with("ratio")),
      n = Inf)
cat("Wrote DOSE Conley analysis to: ", OUT, "\n", sep = "")
