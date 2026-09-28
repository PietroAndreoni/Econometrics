# Aggregate the gridded KUMMU 2025 GDP products to the same GADM 4.1 ADM1
# geography used by the Lamperti/Weighted Climate Dataset inputs.
#
# GDP per capita is aggregated through the accounting identity
#
#   regional GDP pc = sum(cell total GDP) / sum(cell total GDP / cell GDP pc)
#
# so the weights are the population embedded in KUMMU itself. Polygon-edge
# cells are apportioned by their exact covered fraction.
local_library <- normalizePath(".rlib", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(exactextractr)
  library(sf)
  library(terra)
  library(tidyr)
})
sf::sf_use_s2(FALSE)

out_dir <- "results/kummu2025_grid_gadm1"
source_dir <- "results/kummu2025_robustness/sources"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gdp_path <- file.path(source_dir, "rast_gdpTot_1990_2024_5arcmin.tif")
pc_path <- file.path(source_dir, "rast_adm2_gdp_perCapita_1990_2024.tif")
metadata_path <- file.path(source_dir, "zenodo.json")
gadm_dir <- file.path(source_dir, "gadm41")

stopifnot(file.exists(gdp_path), file.exists(pc_path), file.exists(metadata_path))
metadata <- jsonlite::fromJSON(metadata_path)
check_zenodo_file <- function(path) {
  row <- metadata$files[metadata$files$key == basename(path), ]
  stopifnot(nrow(row) == 1L)
  observed <- paste0("md5:", unname(tools::md5sum(path)))
  stopifnot(identical(tolower(observed), tolower(row$checksum)))
}
check_zenodo_file(gdp_path)
check_zenodo_file(pc_path)
message("Check 1 passed: both KUMMU raster checksums match Zenodo record 18429133.")

gdp <- terra::rast(gdp_path)
gdp_pc <- terra::rast(pc_path)
years_gdp <- as.integer(sub("^gdp_tot_", "", names(gdp)))
years_pc <- as.integer(sub("^gdp_pc_", "", names(gdp_pc)))
stopifnot(
  terra::compareGeom(gdp, gdp_pc, stopOnError = FALSE),
  identical(years_gdp, 1990:2024),
  identical(years_pc, years_gdp),
  terra::nlyr(gdp) == 35L
)
message("Check 2 passed: GDP and GDP-per-capita grids have identical 5-arc-minute geometry and 1990-2024 layers.")

# The country files collectively include the special Z01-Z09 GADM territories,
# although those territories do not have standalone country downloads.
gadm_files <- list.files(gadm_dir, pattern = "^[A-Z0-9]{3}\\.json$", full.names = TRUE)
stopifnot(length(gadm_files) > 200L)
message("Reading GADM 4.1 ADM1 polygons ...")
gadm_list <- lapply(gadm_files, sf::st_read, quiet = TRUE)
gadm <- dplyr::bind_rows(gadm_list) %>%
  dplyr::select(GID_0, GID_1, COUNTRY, NAME_1, geometry) %>%
  arrange(GID_1) %>%
  distinct(GID_1, .keep_all = TRUE)
rm(gadm_list)
invisible(gc())

climate_ids <- arrow::open_dataset(
  "data/data_tm_gadm1_era5_pop-area_2000-2015.parquet"
) %>%
  filter(
    gadm_level == "gadm1",
    climate_source == "era5",
    weight == "pop",
    weight_year == "2015"
  ) %>%
  dplyr::select(GID_0, GID_1) %>%
  distinct() %>%
  collect() %>%
  arrange(GID_1)

missing_geometry <- anti_join(climate_ids, st_drop_geometry(gadm), by = c("GID_0", "GID_1"))
extra_geometry <- anti_join(st_drop_geometry(gadm), climate_ids, by = c("GID_0", "GID_1"))
write.csv(missing_geometry, file.path(out_dir, "climate_ids_without_geometry.csv"), row.names = FALSE)
write.csv(extra_geometry, file.path(out_dir, "geometry_without_climate_id.csv"), row.names = FALSE)
stopifnot(nrow(missing_geometry) == 0L)
gadm <- semi_join(gadm, climate_ids, by = c("GID_0", "GID_1")) %>% arrange(GID_1)
stopifnot(!anyDuplicated(gadm$GID_1), nrow(gadm) == nrow(climate_ids))
message("Check 3 passed: every Lamperti climate GADM1 identifier has exactly one GADM 4.1 polygon.")

population_path <- file.path(out_dir, "kummu_population_1990_2024_5arcmin.tif")
if (!file.exists(population_path)) {
  message("Inferring the annual population grid from KUMMU's accounting identity ...")
  population <- gdp / gdp_pc
  names(population) <- paste0("population_", years_gdp)
  terra::writeRaster(
    population, population_path, overwrite = TRUE,
    datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3")
  )
  rm(population)
  invisible(gc())
}

message("Extracting exact-fraction regional GDP totals ...")
gdp_sum <- exactextractr::exact_extract(
  raster::brick(gdp_path), gadm, fun = "sum", progress = TRUE,
  max_cells_in_memory = 3e8
)
message("Extracting exact-fraction regional population totals ...")
pop_sum <- exactextractr::exact_extract(
  raster::brick(population_path), gadm, fun = "sum", progress = TRUE,
  max_cells_in_memory = 3e8
)
gdp_sum <- as.data.frame(gdp_sum)
pop_sum <- as.data.frame(pop_sum)
stopifnot(nrow(gdp_sum) == nrow(gadm), nrow(pop_sum) == nrow(gadm),
          ncol(gdp_sum) == length(years_gdp), ncol(pop_sum) == length(years_gdp))
names(gdp_sum) <- paste0("gdp_", years_gdp)
names(pop_sum) <- paste0("population_", years_gdp)

gdp_long <- gdp_sum %>%
  mutate(GID_1 = gadm$GID_1) %>%
  pivot_longer(starts_with("gdp_"), names_to = "year", values_to = "total_gdp") %>%
  mutate(year = as.integer(sub("^gdp_", "", year)))
pop_long <- pop_sum %>%
  mutate(GID_1 = gadm$GID_1) %>%
  pivot_longer(starts_with("population_"), names_to = "year", values_to = "population") %>%
  mutate(year = as.integer(sub("^population_", "", year)))
econ_base <- inner_join(gdp_long, pop_long, by = c("GID_1", "year"))

econ <- econ_base %>%
  left_join(st_drop_geometry(gadm) %>% dplyr::select(GID_0, GID_1, COUNTRY, NAME_1), by = "GID_1") %>%
  mutate(
    grp_pc_usd = if_else(population > 0, total_gdp / population, NA_real_),
    econ_source = "KUMMU2025_GRID_GADM41",
    gadm_level = "gadm1"
  ) %>%
  group_by(GID_1) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    lgrp_pc_usd = if_else(grp_pc_usd > 0, log(grp_pc_usd), NA_real_),
    dlgrp_pc_usd = if_else(year == lag(year) + 1L, lgrp_pc_usd - lag(lgrp_pc_usd), NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::select(year, GID_0, GID_1, COUNTRY, NAME_1, grp_pc_usd, lgrp_pc_usd,
         dlgrp_pc_usd, total_gdp, population, econ_source, gadm_level)

stopifnot(
  !anyDuplicated(econ[c("GID_1", "year")]),
  all(is.finite(econ$total_gdp)),
  all(is.finite(econ$population)),
  all(econ$total_gdp >= 0),
  all(econ$population >= 0)
)

# External numerical check against KUMMU's tabulated national GDP-per-capita
# file. This is not used to build the regional series.
adm0_raw <- read.csv("data/tabulated_adm0_gdp_perCapita.csv", check.names = FALSE) %>%
  dplyr::select(iso3, matches("^[0-9]{4}$")) %>%
  pivot_longer(matches("^[0-9]{4}$"), names_to = "year", values_to = "official_gdp_pc") %>%
  mutate(year = as.integer(year))
adm0_grid <- econ %>%
  filter(!grepl("^Z[0-9]{2}$", GID_0)) %>%
  group_by(GID_0, year) %>%
  summarise(total_gdp = sum(total_gdp), population = sum(population), .groups = "drop") %>%
  mutate(grid_gdp_pc = total_gdp / population)
adm0_check <- inner_join(adm0_grid, adm0_raw, by = c("GID_0" = "iso3", "year")) %>%
  filter(is.finite(grid_gdp_pc), is.finite(official_gdp_pc), official_gdp_pc > 0) %>%
  mutate(relative_difference = grid_gdp_pc / official_gdp_pc - 1)
write.csv(adm0_check, file.path(out_dir, "adm0_validation.csv"), row.names = FALSE)
validation_summary <- adm0_check %>%
  summarise(
    comparisons = n(),
    countries = n_distinct(GID_0),
    median_abs_relative_difference = median(abs(relative_difference)),
    p95_abs_relative_difference = quantile(abs(relative_difference), .95),
    max_abs_relative_difference = max(abs(relative_difference)),
    correlation = cor(grid_gdp_pc, official_gdp_pc)
  )
write.csv(validation_summary, file.path(out_dir, "adm0_validation_summary.csv"), row.names = FALSE)

arrow::write_parquet(econ, file.path(out_dir, "econ_kummu2025_grid_gadm41.parquet"))
write.csv(
  econ %>% group_by(year) %>% summarise(
    regions = sum(is.finite(grp_pc_usd) & grp_pc_usd > 0),
    total_gdp = sum(total_gdp), population = sum(population), .groups = "drop"
  ),
  file.path(out_dir, "annual_coverage.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "aggregation_session_info.txt"))
message("Wrote full grid-aggregated GADM1 economic panel.")
print(validation_summary)
