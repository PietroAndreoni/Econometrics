# Validate KUMMU reporting units against the GADM 4.1 climate geography.
# Run from the project root. Source files are cached in the results directory.
suppressPackageStartupMessages({library(sf); library(dplyr)})
out <- "results/kummu2025_robustness"
src <- file.path(out, "sources")
raw <- read.csv("data/tabulated_adm1_gdp_perCapita.csv", check.names = FALSE)
polygons <- st_read(file.path(src, "kummu_adm1.gpkg"), quiet = TRUE)
stopifnot(!anyDuplicated(raw$GID_nmbr), !anyDuplicated(polygons$GID_nmbr))
ix <- match(raw$GID_nmbr, polygons$GID_nmbr)
stopifnot(!anyNA(ix), nrow(raw) == nrow(polygons),
          all(raw$iso3 == polygons$iso3[ix]))
# The official CSV contains legacy name bytes; the GeoPackage is UTF-8.
# Verify ASCII names directly and retain the authoritative polygon names.
ascii <- !grepl("[^ -~]", raw$Subnat)
name_check <- data.frame(GID_nmbr = raw$GID_nmbr, csv_name = raw$Subnat,
                        polygon_name = polygons$Subnat[ix])
write.csv(name_check, file.path(out, "source_name_check.csv"), row.names = FALSE)
stopifnot(all(raw$Subnat[ascii] == polygons$Subnat[ix][ascii], na.rm = TRUE))
years <- grep("^[0-9]{4}$", names(raw), value = TRUE)
stopifnot(isTRUE(all.equal(unname(as.matrix(raw[, years])),
                unname(as.matrix(st_drop_geometry(polygons)[ix, paste0("X", years)])))))
metadata <- jsonlite::fromJSON(file.path(src, "zenodo.json"))
source_file <- metadata$files[metadata$files$key == "tabulated_adm1_gdp_perCapita.csv", ]
stopifnot(paste0("md5:", unname(tools::md5sum("data/tabulated_adm1_gdp_perCapita.csv"))) ==
            source_file$checksum)
cat("Check 1 passed: official CSV checksum, region identity, and all annual polygon values.\n")
norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(x, to = "ASCII//TRANSLIT")))
countries <- sort(unique(raw$iso3[raw$GID_nmbr >= 1000000]))
crosswalk <- list()
gadm_names <- list()
for (iso in countries) {
  cat("Spatial audit:", iso, "\n")
  k <- polygons[polygons$iso3 == iso & polygons$GID_nmbr >= 1000000, ]
  g <- st_read(file.path(src, "gadm41", paste0(iso, ".json")), quiet = TRUE)
  gadm_names[[iso]] <- st_drop_geometry(g)
  # Equal-area intersections; both source geographies undergo identical repair
  # and transformation. Requiring near-complete overlap in BOTH directions
  # excludes nesting, aggregates, changed boundaries, and accidental code matches.
  kp <- st_make_valid(st_transform(k, 6933))
  gp <- st_make_valid(st_transform(g, 6933))
  ka <- as.numeric(st_area(kp)); ga <- as.numeric(st_area(gp))
  kp$ki <- seq_len(nrow(kp)); gp$gi <- seq_len(nrow(gp))
  inter <- suppressWarnings(st_intersection(kp[, "ki"], gp[, "gi"]))
  inter$area <- as.numeric(st_area(inter))
  pairs <- st_drop_geometry(inter) %>% group_by(ki, gi) %>%
    summarise(area = sum(area), .groups = "drop") %>%
    group_by(ki) %>% slice_max(area, n = 1, with_ties = FALSE) %>% ungroup()
  m <- match(seq_len(nrow(k)), pairs$ki)
  gi <- pairs$gi[m]
  cw <- data.frame(GID_nmbr = k$GID_nmbr, iso3 = iso, Subnat = k$Subnat,
    old_GID_1 = paste0(iso, ".", k$GID_nmbr %% 1000, "_1"),
    GID_1 = g$GID_1[gi], NAME_1 = g$NAME_1[gi],
    kummu_overlap = pairs$area[m] / ka,
    gadm_overlap = pairs$area[m] / ga[gi])
  cw$name_equal <- norm_name(cw$Subnat) == norm_name(cw$NAME_1)
  cw$accepted <- !is.na(cw$GID_1) & cw$kummu_overlap >= .98 & cw$gadm_overlap >= .98 &
    substr(cw$GID_1, 1, 3) == cw$iso3
  crosswalk[[iso]] <- cw
}
cw <- bind_rows(crosswalk)
# No multiple reporting units may silently become one climate region.
duplicates <- cw %>% filter(accepted) %>% count(GID_1) %>% filter(n > 1)
stopifnot(nrow(duplicates) == 0, !anyDuplicated(cw$GID_nmbr))
write.csv(cw, file.path(out, "regional_crosswalk.csv"), row.names = FALSE)
write.csv(bind_rows(gadm_names), file.path(out, "gadm41_region_names.csv"), row.names = FALSE)
write.csv(cw %>% filter(!accepted), file.path(out, "excluded_regions.csv"), row.names = FALSE)
write.csv(cw %>% filter(accepted, !name_equal), file.path(out, "name_differences.csv"), row.names = FALSE)
write.csv(cw %>% group_by(iso3) %>% summarise(regions = n(),
  changed_codes = sum(accepted & old_GID_1 != GID_1), accepted = sum(accepted), .groups = "drop"),
  file.path(out, "mapping_country_summary.csv"), row.names = FALSE)
cat("Check 2 passed: reciprocal >=98% polygon overlap and unique accepted matches.\n")
print(cw %>% summarise(regions = n(),
  changed_codes = sum(accepted & old_GID_1 != GID_1),
  name_equal_accepted = sum(accepted & name_equal, na.rm = TRUE), accepted = sum(accepted)))
writeLines(capture.output(sessionInfo()), file.path(out, "mapping_session_info.txt"))
