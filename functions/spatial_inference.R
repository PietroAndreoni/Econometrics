# Residual spatial-correlation diagnostics and Conley cutoff selection.
#
# Typical use, given region centroids with columns GID_1, longitude, latitude:
#   pairs <- region_pair_distances(centroids)
#   diag  <- residual_spatial_correlation(model, data, pairs)
#   km    <- choose_conley_cutoff(diag)
#   conley_vcov_check(model, km)

# Centroids of polygons (an sf object with GID_0, GID_1), computed in a global
# equal-area projection and returned as WGS84 longitude and latitude, the
# coordinates vcov_conley() expects.
polygon_centroids <- function(polygons) {
  old_s2 <- sf::sf_use_s2(FALSE)
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  centres <- polygons %>%
    sf::st_make_valid() %>%
    sf::st_transform(6933) %>%
    sf::st_centroid() %>%
    sf::st_transform(4326)
  xy <- sf::st_coordinates(centres)
  out <- sf::st_drop_geometry(centres) %>%
    transmute(GID_0, GID_1, longitude = xy[, "X"], latitude = xy[, "Y"])
  if (any(!is.finite(out$longitude) | !is.finite(out$latitude) |
          abs(out$longitude) > 180 | abs(out$latitude) > 90)) {
    stop("Invalid longitude or latitude from the centroid calculation.")
  }
  out
}

# Great-circle distance in kilometres.
haversine_km <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi / 180
  dlon <- (lon2 - lon1) * to_rad
  dlat <- (lat2 - lat1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# All unordered region pairs within `max_km`.
region_pair_distances <- function(centroids, id = "GID_1", lon = "longitude",
                                  lat = "latitude", max_km = Inf) {
  centroids <- centroids[order(centroids[[id]]), , drop = FALSE]
  n <- nrow(centroids)
  pair_index <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
  tibble::tibble(
    region_i = centroids[[id]][pair_index[, 1]],
    region_j = centroids[[id]][pair_index[, 2]],
    distance_km = haversine_km(
      centroids[[lon]][pair_index[, 1]],
      centroids[[lat]][pair_index[, 1]],
      centroids[[lon]][pair_index[, 2]],
      centroids[[lat]][pair_index[, 2]]
    )
  ) %>%
    filter(distance_km <= max_km)
}

# Mean product of year-standardized residuals by distance bin. Products are
# averaged within year and bin first, so repeated annual observations at the
# same centroid are not treated as distinct spatial units; the band is a t
# interval across years.
residual_spatial_correlation <- function(model, data, pairs, bin_km = 100,
                                         max_km = 3000, id = "GID_1",
                                         year = "year") {
  rows <- fixest::obs(model)
  residual_data <- tibble::tibble(
    unit = data[[id]][rows],
    year = data[[year]][rows],
    residual = as.numeric(stats::residuals(model))
  ) %>%
    group_by(year) %>%
    mutate(
      residual_sd = stats::sd(residual),
      standardized_residual = if_else(
        is.finite(residual_sd) & residual_sd > 0,
        (residual - mean(residual)) / residual_sd,
        NA_real_
      )
    ) %>%
    ungroup()
  pairs <- pairs %>% filter(distance_km <= max_km)

  annual_bins <- bind_rows(lapply(sort(unique(residual_data$year)), function(y) {
    year_residuals <- residual_data %>%
      filter(year == y) %>%
      select(unit, standardized_residual)
    pairs %>%
      left_join(year_residuals, by = c("region_i" = "unit")) %>%
      rename(residual_i = standardized_residual) %>%
      left_join(year_residuals, by = c("region_j" = "unit")) %>%
      rename(residual_j = standardized_residual) %>%
      filter(is.finite(residual_i), is.finite(residual_j)) %>%
      mutate(
        distance_bin = floor(distance_km / bin_km) * bin_km,
        residual_product = residual_i * residual_j
      ) %>%
      group_by(distance_bin) %>%
      summarise(
        spatial_correlation = mean(residual_product),
        region_pairs = n(),
        .groups = "drop"
      ) %>%
      mutate(year = y)
  }))

  annual_bins %>%
    group_by(distance_bin) %>%
    summarise(
      distance_midpoint_km = first(distance_bin) + bin_km / 2,
      mean_spatial_correlation = stats::weighted.mean(
        spatial_correlation, region_pairs, na.rm = TRUE
      ),
      median_annual_correlation = stats::median(spatial_correlation, na.rm = TRUE),
      annual_sd = stats::sd(spatial_correlation, na.rm = TRUE),
      years = sum(is.finite(spatial_correlation)),
      region_year_pairs = sum(region_pairs),
      se_across_years = annual_sd / sqrt(years),
      ci_low = mean_spatial_correlation -
        stats::qt(0.975, pmax(years - 1, 1)) * se_across_years,
      ci_high = mean_spatial_correlation +
        stats::qt(0.975, pmax(years - 1, 1)) * se_across_years,
      .groups = "drop"
    )
}

# First positive distance bin that starts a run of `consecutive_bins` bins whose
# correlation is within `near_zero` of zero with an interval covering zero.
choose_conley_cutoff <- function(diagnostic, near_zero = 0.02,
                                 consecutive_bins = 3L, fallback_km = 3000) {
  diagnostic <- diagnostic %>%
    arrange(distance_bin) %>%
    mutate(
      is_near_zero = abs(mean_spatial_correlation) <= near_zero &
        ci_low <= 0 & ci_high >= 0
    )
  run_ok <- vapply(seq_len(nrow(diagnostic)), function(i) {
    end <- i + consecutive_bins - 1L
    end <= nrow(diagnostic) && all(diagnostic$is_near_zero[i:end])
  }, logical(1))
  eligible <- which(run_ok & diagnostic$distance_bin > 0)
  if (!length(eligible)) {
    warning(
      "No sustained near-zero residual-correlation run found; using ",
      fallback_km, " km."
    )
    return(fallback_km)
  }
  diagnostic$distance_bin[min(eligible)]
}

# Conley VCOV with and without fixest's positive-semidefinite repair, so the
# repair is recorded rather than left as a console warning. Report inference
# with `$fixed`.
conley_vcov_check <- function(model, cutoff_km, lat = ~latitude,
                              lon = ~longitude) {
  conley <- function(fix) {
    suppressWarnings(fixest::vcov_conley(
      model, lat = lat, lon = lon, cutoff = cutoff_km,
      distance = "spherical", vcov_fix = fix
    ))
  }
  raw <- conley(FALSE)
  fixed <- conley(TRUE)
  min_eigen <- function(V) {
    min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
  }
  list(
    raw = raw,
    fixed = fixed,
    diagnostics = tibble::tibble(
      cutoff_km = cutoff_km,
      raw_minimum_eigenvalue = min_eigen(raw),
      fixed_minimum_eigenvalue = min_eigen(fixed),
      psd_correction_applied = min_eigen(raw) < 0
    )
  )
}
