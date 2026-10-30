# Download of World Development Indicators series through the World Bank API
# (https://api.worldbank.org/v2), written in the layout of a WDI "Data Bank"
# csv export (one row per series and country, one column per year "1960
# [YR1960]", ".." for missing, the database and its last update in the footer),
# so read_wdi_csv() reads it like a manual export. Values are written as the API
# returns them: nothing is converted, interpolated or extrapolated.
#
#   source("functions/wdi_download.R")
#   download_wdi(WDI_PPP_SECTOR_SERIES, "data/wb_ppp_sector_data.csv")

WDI_API <- "https://api.worldbank.org/v2"

# GDP per capita at constant prices (market exchange rates and PPP), sector
# value added at constant prices, total GDP and population.
WDI_PPP_SECTOR_SERIES <- c(
  "NY.GDP.PCAP.KD",     # GDP per capita (constant 2015 US$)
  "NY.GDP.PCAP.PP.KD",  # GDP per capita, PPP (constant 2021 international $)
  "NY.GDP.MKTP.KD",     # GDP (constant 2015 US$)
  "NV.AGR.TOTL.KD",     # Agriculture, forestry, and fishing, value added (constant 2015 US$)
  "NV.IND.TOTL.KD",     # Industry (including construction), value added (constant 2015 US$)
  "NV.SRV.TOTL.KD",     # Services, value added (constant 2015 US$)
  "SP.POP.TOTL"         # Population, total
)

# Every country and aggregate, every year, of one indicator (WDI, source 2).
wdi_indicator <- function(code, first_year = 1960L, last_year = as.integer(format(Sys.Date(), "%Y"))) {
  url <- sprintf("%s/country/all/indicator/%s?format=json&source=2&per_page=30000&date=%d:%d",
                 WDI_API, code, first_year, last_year)
  page <- jsonlite::fromJSON(url, simplifyVector = TRUE)
  meta <- page[[1]]
  if (is.null(meta$pages) || meta$pages != 1) stop(code, ": expected one page of results.")
  rows <- page[[2]]
  list(
    data = tibble::tibble(
      series_name = rows$indicator$value, series = rows$indicator$id,
      # The income groups have no ISO3 code in the API; their API id is kept.
      country = rows$country$value,
      GID_0 = ifelse(rows$countryiso3code == "", rows$country$id, rows$countryiso3code),
      year = as.integer(rows$date), value = rows$value
    ),
    last_updated = meta$lastupdated
  )
}

download_wdi <- function(codes, path, first_year = 1960L,
                         last_year = as.integer(format(Sys.Date(), "%Y"))) {
  parts <- lapply(codes, wdi_indicator, first_year = first_year, last_year = last_year)
  long <- dplyr::bind_rows(lapply(parts, `[[`, "data"))
  years <- sort(unique(long$year))
  wide <- long %>%
    dplyr::mutate(value = ifelse(is.na(value), "..", format(value, digits = 17, scientific = FALSE,
                                                            trim = TRUE)),
                  year = paste0(year, " [YR", year, "]")) %>%
    tidyr::pivot_wider(names_from = year, values_from = value) %>%
    dplyr::transmute(`Series Name` = series_name, `Series Code` = series, `Country Name` = country,
                     `Country Code` = GID_0, dplyr::across(dplyr::all_of(paste0(years, " [YR", years, "]"))))
  wide <- wide[order(match(wide$`Series Code`, codes), wide$`Country Name`), ]
  utils::write.csv(wide, path, row.names = FALSE, na = "..")
  updated <- unique(vapply(parts, `[[`, "", "last_updated"))
  cat("\n\n\nData from database: World Development Indicators\n",
      "Last Updated: ", paste(updated, collapse = ", "), "\n", sep = "", file = path, append = TRUE)
  invisible(wide)
}
