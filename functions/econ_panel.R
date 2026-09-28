# Economic panels built from their raw sources, so no preprocessing script has
# to be run first.
#
#   source          level  variable                 prices / units
#   DOSE_V2_14      GADM1  grp_pc_lcu_2015          constant 2015 local currency
#   KUMMU2025_GRID  GADM1  total GDP / population   PPP, constant 2021 int. $
#   PWT110          GADM0  rgdpna / pop             constant national prices,
#                                                   in 2021 PPP US$
#   WB              GADM0  NY.GDP.PCAP.KD           constant 2015 US$, 2015
#                                                   market exchange rates
#
# All four growth rates are real (volume) growth and comparable. Levels are
# not comparable across sources, nor across countries for DOSE (each country
# in its own currency).
#
# DOSE's grp_pc_usd_2015 is NOT used: it is current US$ divided by the US GDP
# deflator, so its growth carries every local-currency movement against the
# dollar (ARG 2002: -61% against -7% at constant local prices).

# DOSE as published; no series is converted, interpolated or extrapolated.
read_dose <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

ECON_RAW_FILES <- c(
  DOSE_V2_14 = "DOSE_V2.14.csv",
  KUMMU2025_GRID = "kummu2025",        # Raw inputs; see kummu_grid.R.
  PWT110 = "pwt110.xlsx",
  WB = "wb_gdp_data.csv"
)

# World Bank WDI "Data Bank" csv export with GDP per capita (NY.GDP.PCAP.KD)
# and population (SP.POP.TOTL): one row per country and series, one column per
# year, ".." for missing. Aggregates (WLD, EUU, ...) are kept; they simply never
# match a GADM country.
read_wdi_csv <- function(path) {
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  year_cols <- grep("^[0-9]{4}", names(raw), value = TRUE)

  raw %>%
    filter(`Series Code` %in% c("NY.GDP.PCAP.KD", "SP.POP.TOTL")) %>%
    mutate(across(all_of(year_cols), as.character)) %>%
    tidyr::pivot_longer(
      all_of(year_cols),
      names_to = "year",
      values_to = "value"
    ) %>%
    transmute(
      GID_0 = as.character(`Country Code`),
      country = as.character(`Country Name`),
      series = `Series Code`,
      year = as.integer(substr(year, 1L, 4L)),
      value = suppressWarnings(as.numeric(if_else(value == "..", NA, value)))
    ) %>%
    tidyr::pivot_wider(names_from = series, values_from = value) %>%
    transmute(
      GID_0,
      country,
      year,
      gdp_pc = NY.GDP.PCAP.KD,
      population = SP.POP.TOTL
    )
}

read_econ_source <- function(stored_source) {
  raw_file <- data_file(ECON_RAW_FILES[[stored_source]])

  if (identical(stored_source, "DOSE_V2_14")) {
    # `grp_pc_usd` keeps its name across sources; for DOSE it holds constant
    # 2015 local currency, as does the agricultural share.
    read_dose(raw_file) %>%
      transmute(
        year,
        GID_0,
        GID_1,
        grp_pc_usd = grp_pc_lcu_2015,
        pop,
        share_ag_gdp = ag_grp_pc_lcu_2015 / grp_pc_lcu_2015,
        econ_source = "DOSE_V2_14",
        gadm_level = "gadm1"
      ) %>%
      filter(!GID_1 %in% c(" ", ""))
  } else if (identical(stored_source, "KUMMU2025_GRID")) {
    # GDP per capita is regional GDP over the population embedded in Kummu's
    # own grids, so it aggregates consistently with total_gdp. Built from the
    # rasters in data/kummu2025 on first use.
    load_kummu2025_grid() %>%
      transmute(
        year = as.integer(year),
        GID_0 = as.character(GID_0),
        GID_1 = as.character(GID_1),
        grp_pc_usd,
        pop = population,
        share_ag_gdp = NA_real_,
        econ_source = "KUMMU2025_GRID",
        gadm_level = "gadm1"
      )
  } else if (identical(stored_source, "WB")) {
    read_wdi_csv(raw_file) %>%
      transmute(
        year,
        GID_0,
        GID_1 = GID_0,
        grp_pc_usd = gdp_pc,          # Constant 2015 US$, market rates.
        pop = population,
        share_ag_gdp = NA_real_,
        econ_source = "WB",
        gadm_level = "gadm0"
      )
  } else {
    # Penn World Table 11.0. `rgdpna` and `rnna` are the constant-national-prices
    # pair, so real GDP and the capital stock are deflated the same way and their
    # ratio is meaningful over time; the chained-PPP series are the ones for
    # cross-country LEVEL comparisons, which is not what this panel is used for.
    # Both are millions of 2021 US$ against millions of people, so the ratios are
    # 2021 PPP US$ per person -- a different basis from the other sources, which is
    # fine because the sources are only ever compared through growth rates.
    readxl::read_excel(raw_file, sheet = "Data") %>%
      mutate(across(c(rgdpna, rnna, pop, emp), as.numeric)) %>%
      transmute(
        year = as.integer(year),
        GID_0 = as.character(countrycode),
        GID_1 = as.character(countrycode),
        grp_pc_usd = rgdpna / pop,      # Real GDP per capita, 2021 US$.
        k_pc_usd = rnna / pop,          # Capital stock per capita, 2021 US$.
        share_emp_pop = emp / pop,      # Persons engaged per head.
        share_ag_gdp = NA_real_,
        econ_source = "PWT110",
        gadm_level = "gadm0"
      )
  }
}

# Growth rates are matched on calendar year, so a gap in a unit's series yields
# NA rather than a growth rate silently spanning the gap. prepare_econ_data.R
# used a positional lag here, which is the one substantive difference between
# this builder and that script.
build_econ_panel <- function(stored_source) {
  read_econ_source(stored_source) %>%
    filter(!is.na(grp_pc_usd), grp_pc_usd > 0) %>%
    group_by(GID_1, gadm_level, econ_source) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      lgrp_pc_usd = log(grp_pc_usd),
      lag_lgrp_pc_usd = lag_by_year(lgrp_pc_usd, year),
      dlgrp_pc_usd = lgrp_pc_usd - lag_lgrp_pc_usd,
      dg = grp_pc_usd / lag_by_year(grp_pc_usd, year) - 1,
      lag_share_ag_gdp = lag_by_year(share_ag_gdp, year)
    ) %>%
    mutate(lag_dlgrp_pc_usd = lag_by_year(dlgrp_pc_usd, year)) %>%
    ungroup() %>%
    mutate(GID_0 = dplyr::coalesce(GID_0, substr(GID_1, 1, 3)))
}

load_econ_panel <- function(stored_source) {
  cached_panel(
    paste0("econ_", tolower(stored_source)),
    list(
      econ_source = stored_source,
      raw_file = unname(ECON_RAW_FILES[[stored_source]]),
      code = code_fingerprint(read_econ_source, build_econ_panel)
    ),
    function() build_econ_panel(stored_source)
  )
}
