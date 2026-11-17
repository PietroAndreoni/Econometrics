# Country groups for the regional aggregates of the response explorer's global
# charts, saved to results/response_explorer/country_groups.rds:
#
#   OECD           the 38 members (2025)
#   EU + UK        the 27 EU members and the United Kingdom
#   BRICS          the five founding members (Brazil, Russia, India, China,
#                  South Africa)
#   BRICS+         with the members that joined in 2024-2025 (Egypt, Ethiopia,
#                  Iran, United Arab Emirates, Indonesia; Saudi Arabia, invited
#                  but not confirmed, left out)
#   latitude       by the population-weighted mean latitude of the country
#                  (Kummu et al. 2025, 2024 population, 5 arcmin, inside the
#                  GADM 4.1 borders; area where there is no population):
#                    equatorial         |lat| < 10
#                    northern tropics   10 <= lat < 23.44 (Tropic of Cancer)
#                    southern tropics   -23.44 < lat <= -10 (Tropic of Capricorn)
#   income         World Bank income groups as of the download (data/
#                  country_groups/wb_countries.json, api.worldbank.org/v2/
#                  country): low, lower-middle, upper-middle and high income
#
# From the project root:  Rscript apps/response_explorer/build_country_groups.R

suppressPackageStartupMessages(library(dplyr))
source("load_functions.R")

GROUP_DIR <- file.path("data", "country_groups")
WB_FILE <- file.path(GROUP_DIR, "wb_countries.json")
POPULATION_LAYER <- "population_2024"
TROPIC <- 23.44
EQUATOR_BAND <- 10

OECD <- c("AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC",
          "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "KOR", "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", "NOR",
          "POL", "PRT", "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA")
EU27 <- c("AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "IRL",
          "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE")
BRICS <- c("BRA", "RUS", "IND", "CHN", "ZAF")
BRICS_PLUS <- c(BRICS, "EGY", "ETH", "IRN", "ARE", "IDN")
stopifnot(length(OECD) == 38, length(EU27) == 27)

# World Bank income groups (countries only: the aggregates have none).
if (!file.exists(WB_FILE)) {
  dir.create(GROUP_DIR, showWarnings = FALSE, recursive = TRUE)
  utils::download.file("https://api.worldbank.org/v2/country?format=json&per_page=400", WB_FILE, mode = "wb")
}
wb <- jsonlite::fromJSON(WB_FILE)[[2]]
wb <- data.frame(iso3 = wb$id, income = wb$incomeLevel$value, stringsAsFactors = FALSE)
wb <- wb[wb$income %in% c("Low income", "Lower middle income", "Upper middle income", "High income"), ]
income <- lapply(c(low = "Low income", lower_middle = "Lower middle income", upper_middle = "Upper middle income",
                   high = "High income"), function(k) sort(wb$iso3[wb$income == k]))

# Population-weighted mean latitude of each country (area where it has none).
old_s2 <- sf::sf_use_s2(FALSE)
regions <- read_gadm1_polygons(file.path(data_root(), "kummu2025", "gadm41"))
population <- terra::rast(file.path(data_root(), "kummu2025", "kummu_population_1990_2024_5arcmin.tif"))
population <- population[[POPULATION_LAYER]]
parts <- exactextractr::exact_extract(
  population, regions, include_xy = TRUE, include_area = TRUE, include_cols = "GID_0", summarize_df = TRUE,
  progress = FALSE,
  fun = function(df) {
    v <- df$value
    v[!is.finite(v)] <- 0
    data.frame(iso3 = df$GID_0[[1]], pop = sum(v * df$coverage_fraction),
               pop_lat = sum(v * df$coverage_fraction * df$y),
               area = sum(df$area * df$coverage_fraction), area_lat = sum(df$area * df$coverage_fraction * df$y))
  })
invisible(sf::sf_use_s2(old_s2))
latitude <- bind_rows(parts) %>%
  group_by(iso3) %>%
  summarise(lat = if (sum(pop) > 0) sum(pop_lat) / sum(pop) else sum(area_lat) / sum(area), .groups = "drop")
band <- function(keep) sort(latitude$iso3[keep])

groups <- list(
  "OECD" = OECD,
  "EU + UK" = c(EU27, "GBR"),
  "BRICS (founding five)" = BRICS,
  "BRICS+ (2025 members)" = BRICS_PLUS,
  "Equatorial (|lat| < 10)" = band(abs(latitude$lat) < EQUATOR_BAND),
  "Northern tropics (10 to 23.4 N)" = band(latitude$lat >= EQUATOR_BAND & latitude$lat < TROPIC),
  "Southern tropics (10 to 23.4 S)" = band(latitude$lat <= -EQUATOR_BAND & latitude$lat > -TROPIC),
  "Low income (WB)" = income$low,
  "Lower-middle income (WB)" = income$lower_middle,
  "Upper-middle income (WB)" = income$upper_middle,
  "High income (WB)" = income$high
)
groups <- lapply(groups, function(x) sort(unique(x)))

out_dir <- file.path("results", "response_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(groups = groups, latitude = latitude,
             wb_downloaded = file.mtime(WB_FILE),
             source = paste("OECD and EU membership 2025; BRICS 2025; World Bank income groups (API, downloaded",
                            format(file.mtime(WB_FILE), "%d %B %Y"), "); population-weighted latitude from",
                            "Kummu et al. (2025) 2024 population and GADM 4.1"),
             built = Sys.time()),
        file.path(out_dir, "country_groups.rds"))
for (k in names(groups)) message(sprintf("%-34s %3d countries", k, length(groups[[k]])))
