# Aggregate the gridded Kummu et al. (2025) GDP products (Zenodo record
# 18429133) to GADM regions. Units: PPP, constant 2021 international dollars.
#
# GDP per capita is aggregated through the accounting identity
#
#   regional GDP pc = sum(cell total GDP) / sum(cell total GDP / cell GDP pc)
#
# so the weights are the population embedded in Kummu itself. Polygon-edge
# cells are apportioned by their exact covered fraction. Needs sf, terra,
# raster and exactextractr, but only when the panel has to be (re)built.
#
# Raw inputs, in data/kummu2025/:
#   rast_gdpTot_1990_2024_5arcmin.tif, rast_adm2_gdp_perCapita_1990_2024.tif
#   and zenodo.json from Zenodo record 18429133;
#   gadm41/<ISO>.json, GADM 4.1 ADM1 polygons, one file per country;
#   kummu_population_1990_2024_5arcmin.tif, the population grid inferred from
#   the two rasters (recreated automatically if deleted).
# load_kummu2025_grid() is what the economic readers call; it builds the GADM1
# panel on first use and caches it in data/cache.

kummu2025_dir <- function() file.path(data_root(), "kummu2025")

# Full build: checksum checks, polygons, extraction. `keep_ids = NULL` keeps
# every GADM 4.1 ADM1 region.
build_kummu2025_grid <- function(source_dir = kummu2025_dir(), keep_ids = NULL) {
  gdp_path <- file.path(source_dir, "rast_gdpTot_1990_2024_5arcmin.tif")
  pc_path <- file.path(source_dir, "rast_adm2_gdp_perCapita_1990_2024.tif")
  missing_inputs <- c(gdp_path, pc_path)[!file.exists(c(gdp_path, pc_path))]
  if (length(missing_inputs)) {
    stop("Kummu 2025 rasters not found: ", paste(missing_inputs, collapse = ", "))
  }

  metadata <- jsonlite::fromJSON(file.path(source_dir, "zenodo.json"))
  check_zenodo_file(gdp_path, metadata)
  check_zenodo_file(pc_path, metadata)

  old_s2 <- sf::sf_use_s2(FALSE)
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  regions <- read_gadm1_polygons(file.path(source_dir, "gadm41"), keep_ids)
  message("Aggregating Kummu 2025 grids to ", nrow(regions), " GADM 4.1 regions.")

  aggregate_kummu2025_grid(
    gdp_path,
    pc_path,
    regions,
    population_path = file.path(
      source_dir, "kummu_population_1990_2024_5arcmin.tif"
    )
  )
}

# Cached GADM1 panel; rebuilt when the inputs, `keep_ids` or the builder code
# change. A full rebuild takes several minutes.
load_kummu2025_grid <- function(keep_ids = NULL) {
  cached_panel(
    "econ_kummu2025_grid_gadm41",
    list(
      rasters = unname(tools::md5sum(file.path(
        kummu2025_dir(),
        c("rast_gdpTot_1990_2024_5arcmin.tif",
          "rast_adm2_gdp_perCapita_1990_2024.tif")
      ))),
      keep_ids = if (is.null(keep_ids)) "all" else sort(keep_ids),
      code = code_fingerprint(
        build_kummu2025_grid, read_gadm1_polygons, aggregate_kummu2025_grid
      )
    ),
    function() build_kummu2025_grid(keep_ids = keep_ids)
  )
}

# Stops unless `path` matches the md5 checksum listed in the Zenodo metadata.
check_zenodo_file <- function(path, metadata) {
  row <- metadata$files[metadata$files$key == basename(path), ]
  if (nrow(row) != 1L) {
    stop(basename(path), " is not listed in the Zenodo metadata.")
  }
  observed <- paste0("md5:", unname(tools::md5sum(path)))
  if (!identical(tolower(observed), tolower(row$checksum))) {
    stop("Checksum mismatch for ", basename(path), ".")
  }
  invisible(TRUE)
}

# GADM 4.1 ADM1 polygons from one GeoJSON per country (gadm41_<ISO>_1.json
# downloads renamed <ISO>.json), optionally restricted to `keep_ids`.
read_gadm1_polygons <- function(gadm_dir, keep_ids = NULL) {
  files <- list.files(gadm_dir, pattern = "^[A-Z0-9]{3}\\.json$",
                      full.names = TRUE)
  if (!length(files)) {
    stop("No <ISO>.json GADM files found in ", gadm_dir)
  }
  # Some country files carry empty features whose identifiers are the string
  # "NA" (GBR.json has one); they match nothing and are dropped.
  polygons <- dplyr::bind_rows(lapply(files, sf::st_read, quiet = TRUE)) %>%
    select(GID_0, GID_1, COUNTRY, NAME_1, geometry) %>%
    filter(grepl("^[A-Z0-9]{3}\\.", GID_1)) %>%
    arrange(GID_1) %>%
    distinct(GID_1, .keep_all = TRUE)
  if (!is.null(keep_ids)) {
    missing_ids <- setdiff(keep_ids, polygons$GID_1)
    if (length(missing_ids)) {
      warning(length(missing_ids), " requested GID_1 have no polygon, e.g. ",
              paste(utils::head(missing_ids, 5L), collapse = ", "))
    }
    polygons <- polygons[polygons$GID_1 %in% keep_ids, ]
  }
  polygons
}

# Regional panel of total GDP, population and GDP per capita. The inferred
# population grid (GDP / GDP pc) is written to `population_path` the first
# time and reused afterwards.
aggregate_kummu2025_grid <- function(gdp_path, pc_path, regions,
                                     population_path,
                                     econ_source = "KUMMU2025_GRID") {
  gdp <- terra::rast(gdp_path)
  gdp_pc <- terra::rast(pc_path)
  years <- as.integer(sub("^gdp_tot_", "", names(gdp)))
  years_pc <- as.integer(sub("^gdp_pc_", "", names(gdp_pc)))
  if (!terra::compareGeom(gdp, gdp_pc, stopOnError = FALSE) ||
      !identical(years, years_pc)) {
    stop("GDP and GDP-per-capita grids differ in geometry or years.")
  }

  if (!file.exists(population_path)) {
    message("Inferring the annual population grid from Kummu's identity ...")
    population <- gdp / gdp_pc
    names(population) <- paste0("population_", years)
    terra::writeRaster(
      population, population_path, overwrite = TRUE,
      datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3")
    )
    rm(population)
    invisible(gc())
  }

  extract_sum <- function(path, prefix) {
    message("Extracting exact-fraction regional sums from ", basename(path))
    out <- as.data.frame(exactextractr::exact_extract(
      raster::brick(path), regions, fun = "sum", progress = FALSE,
      max_cells_in_memory = 3e8
    ))
    names(out) <- paste0(prefix, years)
    out$GID_1 <- regions$GID_1
    tidyr::pivot_longer(
      out,
      starts_with(prefix),
      names_to = "year",
      values_to = sub("_$", "", prefix)
    ) %>%
      mutate(year = as.integer(sub(prefix, "", year, fixed = TRUE)))
  }

  panel <- inner_join(
    extract_sum(gdp_path, "total_gdp_"),
    extract_sum(population_path, "population_"),
    by = c("GID_1", "year")
  ) %>%
    left_join(
      sf::st_drop_geometry(regions) %>% select(GID_0, GID_1, COUNTRY, NAME_1),
      by = "GID_1"
    ) %>%
    mutate(
      grp_pc_usd = if_else(population > 0, total_gdp / population, NA_real_),
      econ_source = econ_source,
      gadm_level = "gadm1"
    ) %>%
    select(year, GID_0, GID_1, COUNTRY, NAME_1, grp_pc_usd, total_gdp,
           population, econ_source, gadm_level) %>%
    arrange(GID_1, year)

  stopifnot(
    !anyDuplicated(panel[c("GID_1", "year")]),
    all(is.finite(panel$total_gdp)),
    all(is.finite(panel$population)),
    all(panel$total_gdp >= 0),
    all(panel$population >= 0)
  )
  panel
}

# National GDP per capita re-aggregated from the regional panel against Kummu's
# own tabulated national file (not used to build the regional series).
validate_kummu_adm0 <- function(panel, adm0_csv) {
  official <- utils::read.csv(adm0_csv, check.names = FALSE) %>%
    select(iso3, matches("^[0-9]{4}$")) %>%
    tidyr::pivot_longer(matches("^[0-9]{4}$"), names_to = "year",
                        values_to = "official_gdp_pc") %>%
    mutate(year = as.integer(year))

  panel %>%
    filter(!grepl("^Z[0-9]{2}$", GID_0)) %>%
    group_by(GID_0, year) %>%
    summarise(total_gdp = sum(total_gdp), population = sum(population),
              .groups = "drop") %>%
    mutate(grid_gdp_pc = total_gdp / population) %>%
    inner_join(official, by = c("GID_0" = "iso3", "year")) %>%
    filter(is.finite(grid_gdp_pc), is.finite(official_gdp_pc),
           official_gdp_pc > 0) %>%
    mutate(relative_difference = grid_gdp_pc / official_gdp_pc - 1)
}
