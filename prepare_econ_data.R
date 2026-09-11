suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(readxl))

# Economic preprocessing selectors
SELECT_ECON_SOURCE <- c("DOSE_V2_11", "KUMMU2018","KUMMU2025","WDI","WB","PWT110")

read_kummu2025_tabulated <- function(path, gadm_level, gid_col, gid_prefix = NULL) {
  data_raw <- read.csv(path, check.names = FALSE)
  year_cols <- names(data_raw)[grepl("^\\d{4}$", names(data_raw))]

  data_raw %>%
    mutate(
      GID_0 = as.character(iso3),
      GID_1 = as.character(.data[[gid_col]]),
      GID_1 = if (is.null(gid_prefix)) GID_1 else paste0(gid_prefix, GID_1),
      share_ag_gdp = NA_real_
    ) %>%
    filter(!is.na(GID_0), !is.na(GID_1), !GID_1 %in% c("", "NA")) %>%
    pivot_longer(
      all_of(year_cols),
      names_to = "year",
      values_to = "grp_pc_usd"
    ) %>%
    mutate(
      year = as.integer(year),
      grp_pc_usd = as.numeric(grp_pc_usd),
      econ_source = "KUMMU2025",
      gadm_level = gadm_level
    ) %>%
    select(year, GID_0, GID_1, grp_pc_usd, share_ag_gdp, econ_source, gadm_level)
}

read_kummu2025_adm1 <- function(path) {
  data_raw <- read.csv(path, check.names = FALSE)
  year_cols <- names(data_raw)[grepl("^\\d{4}$", names(data_raw))]

  data_raw %>%
    mutate(
      GID_nmbr = as.integer(GID_nmbr),
      adm_level_digit = GID_nmbr %/% 1000000L,
      adm1_index = GID_nmbr %% 1000L
    ) %>%
    # Rows below 1,000,000 are country-level records already supplied by adm0.
    filter(adm_level_digit == 1L) %>%
    mutate(
      GID_0 = as.character(iso3),
      GID_1 = paste0(GID_0, ".", adm1_index, "_1"),
      share_ag_gdp = NA_real_
    ) %>%
    pivot_longer(
      all_of(year_cols),
      names_to = "year",
      values_to = "grp_pc_usd"
    ) %>%
    mutate(
      year = as.integer(year),
      grp_pc_usd = as.numeric(grp_pc_usd),
      econ_source = "KUMMU2025",
      gadm_level = "gadm1"
    ) %>%
    select(year, GID_0, GID_1, grp_pc_usd, share_ag_gdp, econ_source, gadm_level)
}

econ_data <- read.csv("econometrics/data/DOSE_V2.11.csv") %>%
  select(year, GID_0, GID_1, grp_pc_usd_2015, ag_grp_pc_usd_2015, pop) %>%
  rename(grp_pc_usd = grp_pc_usd_2015) %>%
  mutate(share_ag_gdp = ag_grp_pc_usd_2015/grp_pc_usd) %>% select(-ag_grp_pc_usd_2015) %>% 
  filter(!GID_1 %in% c(" ", "")) %>%
  mutate(econ_source = "DOSE_V2_11",gadm_level="gadm1")


econ_data <- econ_data %>%
  bind_rows(
    read.csv("econometrics/data/kummu_aggregated_gid1_gdp.csv") %>%
      mutate(GID_0 = stringr::str_extract(GID_1, "^.{3}")) %>%
      select(year, GID_0, GID_1, gdp_pc) %>%
      rename(grp_pc_usd = gdp_pc) %>%
      mutate(econ_source = "KUMMU2018",gadm_level="gadm1")
  )

econ_data <- econ_data %>%
  bind_rows(
    read_kummu2025_tabulated(
      "econometrics/data/tabulated_adm0_gdp_perCapita.csv",
      gadm_level = "gadm0",
      gid_col = "iso3"
    ),
    read_kummu2025_adm1("econometrics/data/tabulated_adm1_gdp_perCapita.csv"),
    read_kummu2025_tabulated(
      "econometrics/data/tabulated_adm2_gdp_perCapita.csv",
      gadm_level = "gadm2",
      gid_col = "GID_2"
    )
  )

econ_data <- econ_data %>%
  bind_rows(
    read.csv("econometrics/data/DOSE_V2.csv") %>%
      select(year, GID_0, GID_1, grp_pc_usd_2015) %>%
      rename(grp_pc_usd = grp_pc_usd_2015) %>%
      filter(!GID_1 %in% c(" ", "")) %>%
      mutate(econ_source = "DOSE_V2",gadm_level="gadm1")
  )

econ_data <- econ_data %>%
  bind_rows(arrow::read_parquet("econometrics/data/data_gdp_gilli.parquet") %>% 
  select(year,iso3,gdppc) %>%
  rename(grp_pc_usd = gdppc, GID_1=iso3) %>% 
  mutate(econ_source = "WDI",gadm_level="gadm0") )
  
econ_data <- econ_data %>%
  bind_rows(read.csv("econometrics/data/wb_gdp_data.csv") %>% 
              pivot_longer(!c(Series.Name,Series.Code,Country.Code,Country.Name),names_to="year",values_to="grp_pc_usd") %>% 
              select(Country.Code,year,grp_pc_usd) %>%
              rename(GID_1=Country.Code) %>% 
              mutate(econ_source = "WB",gadm_level="gadm0",
                     year=as.integer(substr(gsub("^X", "", as.character(year)), 1, 4)),
            grp_pc_usd=as.numeric(ifelse(grp_pc_usd=="..",NA,grp_pc_usd))) )

# Penn World Table 11.0 (national accounts, country level, 1950-2023).
#
# rgdpna and rnna are the constant-national-prices pair, so real GDP and the
# capital stock are deflated the same way and their ratio is meaningful over
# time; the chained-PPP series (rgdpe/rgdpo) are the ones for cross-country
# LEVEL comparisons, which is not what this panel is used for.
#
# Units: both are millions of 2021 US$ and pop/emp are millions of people, so
# the ratios come out directly in 2021 US$ per person. NOTE this is a different
# base year from every other source here, which is 2015 US$ -- fine as long as
# sources are compared through growth rates or kept apart (the processing below
# groups by econ_source), but PWT levels are NOT on the same scale as the rest.
econ_data <- econ_data %>%
  bind_rows(
    read_excel("econometrics/data/pwt110.xlsx", sheet = "Data") %>%
      mutate(across(c(rgdpna, rnna, pop, emp), as.numeric)) %>%
      transmute(
        year = as.integer(year),
        GID_0 = as.character(countrycode),
        GID_1 = as.character(countrycode),
        grp_pc_usd = rgdpna / pop,      # real GDP per capita, 2021 US$
        k_pc_usd = rnna / pop,          # capital stock per capita, 2021 US$
        share_emp_pop = emp / pop,      # persons engaged / population
        econ_source = "PWT110",
        gadm_level = "gadm0"
      )
  )

available_econ_sources <- sort(unique(econ_data$econ_source))
econ_data <- econ_data %>%
  filter(econ_source %in% SELECT_ECON_SOURCE & !is.na(grp_pc_usd) & grp_pc_usd > 0)

if (nrow(econ_data) == 0) {
  stop(
    paste0(
      "No rows found for SELECT_ECON_SOURCE in c(",
      paste(sprintf("'%s'", SELECT_ECON_SOURCE), collapse = ", "),
      "). Available: ",
      paste(available_econ_sources, collapse = ", ")
    )
  )
}

econ_data_processed <- econ_data %>%
  group_by(GID_1, gadm_level, econ_source) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    dlgrp_pc_usd = log(grp_pc_usd) - log(lag(grp_pc_usd)),
    dg = (grp_pc_usd - lag(grp_pc_usd)) / lag(grp_pc_usd),
    lgrp_pc_usd = log(grp_pc_usd),
    lag_lgrp_pc_usd = log(lag(grp_pc_usd)),
    lag_share_ag_gdp = lag(share_ag_gdp)
  ) %>%
  mutate(lag_dlgrp_pc_usd = lag(dlgrp_pc_usd)) %>% 
  ungroup() %>% 
  mutate(GID_0 = dplyr::coalesce(GID_0, stringr::str_extract(GID_1, "^[A-Z0-9]{3}")))

deflator <- read.csv("econometrics/data/DOSE_V2.11.csv") %>%
  select(year, GID_0, GID_1, cpi_2015, deflator_2015, fx)

sanitize_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
join_tag <- function(x) paste(sanitize_tag(unique(as.character(x))), collapse = "-")

econ_tag <- join_tag(SELECT_ECON_SOURCE)
econ_data_parquet <- file.path("econometrics", "data", paste0("econ_processed_", econ_tag, ".parquet"))

arrow::write_parquet(econ_data_processed, econ_data_parquet)

cat("Wrote:\n", econ_data_parquet, "\n", sep = "")
