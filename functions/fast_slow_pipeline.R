# Bookkeeping for the fast-versus-slow warming pipeline (scripts
# fast_slow_stage*.R): the frozen design and its hash, output paths, stage
# status files, the model registry, and the fixed-effect designs of Stage 4.

FS_DESIGN_FILE <- file.path("config", "fast_slow_design.yml")
FS_OUT <- file.path("results", "fast_slow_warming")

fs_design <- function(path = FS_DESIGN_FILE) {
  design <- yaml::read_yaml(path)
  design$hash <- digest::digest(file = path, algo = "sha256")
  design
}

fs_path <- function(...) {
  path <- file.path(FS_OUT, ...)
  for (d in unique(dirname(path))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

fs_write <- function(x, ..., design = NULL) {
  path <- fs_path(...)
  if (!is.null(design) && is.data.frame(x) && nrow(x)) {
    x$design_hash <- substr(design$hash, 1, 16)
  }
  if (grepl("\\.parquet$", path)) {
    arrow::write_parquet(x, path)
  } else if (grepl("\\.json$", path)) {
    writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE,
                                digits = NA, null = "null", na = "null"),
               path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

fs_hash_files <- function(paths) {
  paths <- paths[file.exists(paths)]
  vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256")
}

# stage_status.json: PASS, PASS_WITH_WARNING or STOP, with the reasons and the
# hashes of the stage's inputs and outputs.
fs_stage_status <- function(stage, status, reasons, inputs = character(),
                            outputs = character(), design) {
  stopifnot(status %in% c("PASS", "PASS_WITH_WARNING", "STOP", "SKIPPED",
                          "SKIPPED_PREDECLARATION_MISSING"))
  record <- list(
    stage = stage, status = status, reasons = as.list(reasons),
    design_hash = design$hash,
    inputs = as.list(fs_hash_files(inputs)),
    outputs = as.list(fs_hash_files(outputs)),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  fs_write(record, "qa", paste0("stage_status_", stage, ".json"))
  message("Stage ", stage, ": ", status,
          if (length(reasons)) paste0(" (", paste(reasons, collapse = "; "), ")"))
  invisible(record)
}

# Fixed-effect designs of Stage 4, in fixest syntax. Country trends are slopes
# on the centred year (year_c) nested in the region effects.
FS_FIXED_EFFECTS <- c(
  broad_primary = "GID_1 + year + GID_0[[year_c]]",
  no_trend = "GID_1 + year",
  quadratic_trend = "GID_1 + year + GID_0[[year_c]] + GID_0[[year_c2]]",
  country_year = "GID_1 + cy",
  zone_year = "GID_1 + zy"
)

# Explicit country-year (cy) and climate-zone-year (zy) identifiers for the
# interacted fixed effects (fixest cannot re-demean a fit with `a^b`).
fs_add_fe_columns <- function(data) {
  data$cy <- paste(data$GID_0, data$year)
  data$zy <- ifelse(is.na(data$zone), NA_character_,
                    paste(data$zone, data$year))
  data
}

FS_LEVEL_TERMS <- c("Tc", "Tc2")
FS_PRECIP_TERMS <- c("Pz", "Pz2")

fs_formula <- function(terms, outcome = "g", fe = FS_FIXED_EFFECTS[["broad_primary"]]) {
  stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + "), "|", fe))
}

fs_fit <- function(terms, data, outcome = "g",
                   fe = FS_FIXED_EFFECTS[["broad_primary"]], weights = NULL) {
  fixest::feols(
    fs_formula(terms, outcome, fe), data = data,
    vcov = ~ GID_0 + year, weights = weights, fixef.rm = "none",
    fixef.tol = 1e-10, fixef.iter = 50000, notes = FALSE, warn = FALSE
  )
}

# One registry row per model (section 8.4).
fs_registry_row <- function(model_id, tier, fit, data, design, rate_definition,
                            level_function = "centered_quadratic",
                            rate_function = "signed_quadratic_hinge",
                            fixed_effects = "broad_primary",
                            outcome = "100*dlog grp_pc_lcu_2015",
                            estimator = "OLS", regression_weight = "none",
                            cluster_method = "twoway GID_0 + year",
                            sample_id = "primary", warnings = "") {
  used <- data[fixest::obs(fit), ]
  tibble::tibble(
    model_id, tier, outcome,
    climate_product = design$climate$primary_product,
    climate_weight = design$climate$primary_weight,
    rate_definition, level_function, rate_function,
    fixed_effects,
    trends = unname(c(broad_primary = "country linear", no_trend = "none",
                      quadratic_trend = "country quadratic",
                      country_year = "absorbed", zone_year = "none")[
                        fixed_effects]),
    controls = "Pz + Pz2", estimator, regression_weight, cluster_method,
    sample_id, n_obs = stats::nobs(fit),
    n_regions = dplyr::n_distinct(used$GID_1),
    n_countries = dplyr::n_distinct(used$GID_0),
    n_years = dplyr::n_distinct(used$year),
    design_hash = substr(design$hash, 1, 16),
    manifest_hash = design$manifest_hash %||% NA_character_,
    convergence = isTRUE(fit$convStatus) || is.null(fit$convStatus),
    warnings
  )
}
