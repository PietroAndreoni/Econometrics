deflator <- read.csv("econometrics/data/DOSE_V2.11.csv") %>% 
  select(year,GID_0,GID_1,cpi_2015, deflator_2015, fx) 

econ_data_subnat <- read.csv("econometrics/data/DOSE_V2.11.csv") %>% 
  select(year,GID_0,GID_1,grp_pc_usd_2015) %>% 
  rename(grp_pc_usd = grp_pc_usd_2015) %>% 
  filter(!GID_1 %in% c(" ","")) %>% 
  mutate(econ_source="DOSE_V2.11")

econ_data_subnat <- econ_data_subnat %>% 
  bind_rows(read.csv("econometrics/data/kummu_aggregated_gid1_gdp.csv") %>% 
  mutate(GID_0=stringr::str_extract(GID_1, "^.{3}")) %>% 
  select(year,GID_0,GID_1,gdp_pc) %>% 
  rename(grp_pc_usd = gdp_pc) %>% 
  mutate(econ_source="KUMMU2018"))

# econ_data_subnat <- econ_data_subnat %>%
#   bind_rows(read.csv("econometrics/data/DOSE_V1.csv") %>%
#   mutate(GID_0=iso,GID_1=paste0(iso,".",id_1,"_1"),grp_pc_usd=gdp_pc_usd) %>% 
#   select(year,GID_0,GID_1,grp_pc_usd) %>%
#   filter(!GID_1 %in% c(" ","")) %>%
#   mutate(econ_source="DOSE_V1"))

econ_data_subnat <- econ_data_subnat %>% 
  bind_rows(read.csv("econometrics/data/DOSE_V2.csv") %>% 
  select(year,GID_0,GID_1,grp_pc_usd_2015) %>% 
  rename(grp_pc_usd = grp_pc_usd_2015) %>% 
  filter(!GID_1 %in% c(" ","")) %>% 
  mutate(econ_source="DOSE_V2"))

econ_data_subnat <- econ_data_subnat %>% 
  group_by(GID_1,econ_source) %>%
  arrange(year,.by_group=TRUE) %>% 
  ungroup()

# Slice-processing configuration to reduce memory footprint.
# You can pass one or multiple values for each selector.
SELECT_ECON_SOURCE <- "DOSE_V2.11"
SELECT_GADM_LEVEL <- "gadm1"              # any of: gadm0, gadm1
SELECT_CLIMATE_SOURCE <- c("cru","era5")  # any of: cru, dela, era5
SELECT_WEIGHT <- "pop"                     # any of: pop, area, cropland, lights
SELECT_WEIGHT_YEAR <- "2015"               # e.g. 2015, un, concurrent, unspecified
WINDOWS <- c(5,10,15,20)          # moving average windows in years (e.g., 10 for 10-year rolling averages)
# Optional strict intersection filter:
# keep only GID_1-year pairs that are present in all economic datasets
# (non-missing dlgrp_pc_usd) and all climate datasets (non-missing TM and RR).
REQUIRE_COMPLETE_YEAR_GID1 <- FALSE

# Load all subnational climate parquet files and aggregate monthly values to yearly:
# - temperature (tmp): yearly mean
# - precipitation (pre): yearly sum
parse_climate_file_metadata <- function(path) {
  f <- basename(path)
  m <- stringr::str_match(f, "^(gadm[01])_([a-z0-9]+)_(tmp|pre)_(.+?)_monthly\\.parquet$")
  if (is.na(m[1, 1])) {
    stop(paste("Unexpected climate parquet filename:", f))
  }

  gadm_level <- m[1, 2]
  source_raw <- m[1, 3]
  variable <- m[1, 4]
  weight_raw <- tolower(m[1, 5])

  # Parse optional weighting vintage/scheme suffix in the filename
  # (e.g. pop_2015, lights_concurrent, un).
  if (grepl("_(\\d{4}|un|concurrent)$", weight_raw)) {
    weight_year <- sub("^.*_(\\d{4}|un|concurrent)$", "\\1", weight_raw)
    weight_base <- sub("_(\\d{4}|un|concurrent)$", "", weight_raw)
  } else if (grepl("^un_*$", weight_raw)) {
    weight_year <- "un"
    weight_base <- "un"
  } else {
    weight_year <- "unspecified"
    weight_base <- weight_raw
  }

  weight <- dplyr::case_when(
    weight_base == "pop" ~ "pop",
    weight_base == "cropland" ~ "cropland",
    weight_base == "lights" ~ "lights",
    grepl("^un_*$", weight_base) ~ "area",
    TRUE ~ weight_base
  )

  climate_source <- dplyr::case_when(
    source_raw == "era" ~ "era5",
    TRUE ~ source_raw
  )

  list(
    gadm_level = gadm_level,
    variable = variable,
    weight = weight,
    weight_year = weight_year,
    climate_source = climate_source
  )
}

aggregate_climate_parquet <- function(path) {
  meta <- parse_climate_file_metadata(path)
  yearly_fn <- if (meta$variable == "tmp") mean else sum

  arrow::read_parquet(path) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(-Date, names_to = "GID_1", values_to = "value") %>%
    filter(!is.na(GID_1), !is.na(value)) %>%
    mutate(
      # Align identifier format with other data sources:
      # - gadm1: convert ISO3_XXX to ISO3.XXX (e.g., PHL_110_1 -> PHL.110_1)
      # - gadm0: keep ISO3 code as-is in GID_1
      GID_1 = ifelse(meta$gadm_level == "gadm1", sub("^([^_]+)_", "\\1.", GID_1), as.character(GID_1)),
      # Date in source files can be "XYYYYMM", "YYYY-MM-DD", Date, etc.
      year = as.integer(substr(gsub("^X", "", as.character(Date)), 1, 4))
    ) %>%
    # Some files (e.g., cru_tmp_un) contain NA climate values for some months/regions.
    # Drop missing values before yearly aggregation.
    group_by(GID_1, year) %>%
    summarise(
      value = yearly_fn(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      variable = meta$variable,
      weight = meta$weight,
      weight_year = meta$weight_year,
      climate_source = meta$climate_source,
      gadm_level = meta$gadm_level
    )
}

climate_parquet_files <- list.files(
  "econometrics/data",
  pattern = "^gadm[01]_.*_monthly\\.parquet$",
  full.names = TRUE
)

climate_file_meta <- tibble::tibble(path = climate_parquet_files) %>%
  mutate(
    meta = purrr::map(path, parse_climate_file_metadata),
    gadm_level = purrr::map_chr(meta, "gadm_level"),
    climate_source = purrr::map_chr(meta, "climate_source"),
    weight = purrr::map_chr(meta, "weight"),
    weight_year = purrr::map_chr(meta, "weight_year"),
    variable = purrr::map_chr(meta, "variable")
  ) %>%
  select(-meta)

available_econ_sources <- sort(unique(econ_data_subnat$econ_source))
econ_data_subnat <- econ_data_subnat %>%
  filter(econ_source %in% SELECT_ECON_SOURCE)

if (nrow(econ_data_subnat) == 0) {
  stop(
    paste0(
      "No rows found for SELECT_ECON_SOURCE in c(",
      paste(sprintf("'%s'", SELECT_ECON_SOURCE), collapse = ", "),
      "). ",
      "Available econ_source values: ",
      paste(available_econ_sources, collapse = ", ")
    )
  )
}

selected_climate_meta <- climate_file_meta %>%
  filter(
    gadm_level %in% SELECT_GADM_LEVEL,
    climate_source %in% SELECT_CLIMATE_SOURCE,
    weight %in% SELECT_WEIGHT,
    weight_year %in% SELECT_WEIGHT_YEAR
  )

if (nrow(selected_climate_meta) == 0) {
  stop(
    paste0(
      "No climate files found for selectors: ",
      "SELECT_GADM_LEVEL in c(",
      paste(sprintf("'%s'", SELECT_GADM_LEVEL), collapse = ", "),
      "), ",
      "SELECT_CLIMATE_SOURCE in c(",
      paste(sprintf("'%s'", SELECT_CLIMATE_SOURCE), collapse = ", "),
      "), SELECT_WEIGHT in c(",
      paste(sprintf("'%s'", SELECT_WEIGHT), collapse = ", "),
      "), SELECT_WEIGHT_YEAR in c(",
      paste(sprintf("'%s'", SELECT_WEIGHT_YEAR), collapse = ", "),
      ")."
    )
  )
}

selected_balance <- selected_climate_meta %>%
  count(gadm_level, climate_source, weight, weight_year, variable, name = "n_files") %>%
  tidyr::pivot_wider(names_from = variable, values_from = n_files, values_fill = 0)

selected_unbalanced <- selected_balance %>%
  filter(tmp == 0 | pre == 0)

if (nrow(selected_unbalanced) > 0) {
  print(selected_unbalanced)
  stop("Selected climate subset is unbalanced: some climate_source/weight/weight_year combinations are missing tmp or pre.")
}

print(
  paste0(
    "Loading selected climate files for: ",
    "level=", paste(SELECT_GADM_LEVEL, collapse = ","),
    " | ",
    "sources=", paste(SELECT_CLIMATE_SOURCE, collapse = ","),
    " | weights=", paste(SELECT_WEIGHT, collapse = ","),
    " | weight_year=", paste(SELECT_WEIGHT_YEAR, collapse = ","),
    " ..."
  )
)
climate_data_long <- purrr::map_dfr(selected_climate_meta$path, aggregate_climate_parquet)

climate_data_all <- climate_data_long %>%
  tidyr::pivot_wider(
    names_from = variable,
    values_from = value
  )

if (!all(c("tmp", "pre") %in% names(climate_data_all))) {
  stop("Pivoted climate data is missing 'tmp' and/or 'pre' columns.")
}

climate_data_all <- climate_data_all %>%
  rename(TM = tmp, RR = pre) %>%
  mutate(
    # Keep RR unit aligned with previous pipeline convention.
    RR = RR / 1000
  ) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  filter(any(dplyr::coalesce(TM, 0) != 0 | dplyr::coalesce(RR, 0) != 0)) %>%
  arrange(year, .by_group = TRUE) %>%
  ungroup()

# simple means and sd over the full dataset period
means_all <- climate_data_all %>%
  filter(year >= 1990 & year <= 2019) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  summarise(
    sd_TM_all = sd(TM, na.rm = TRUE),
    sd_RR_all = sd(RR, na.rm = TRUE),
    mean_TM_all = mean(TM, na.rm = TRUE),
    mean_RR_all = mean(RR, na.rm = TRUE),
    .groups = "drop"
  )

# simple means and sd pre-dataset
means_pre <- climate_data_all %>%
  filter(year >= 1960 & year <= 1989) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  summarise(
    sd_TM_pre = sd(TM, na.rm = TRUE),
    sd_RR_pre = sd(RR, na.rm = TRUE),
    mean_TM_pre = mean(TM, na.rm = TRUE),
    mean_RR_pre = mean(RR, na.rm = TRUE),
    .groups = "drop"
  )

# calculate moving averages with different windows
data <- data.frame()
for (window in WINDOWS) {
  print(paste("Calculating moving average with a window of",window,"years"))

  data <- bind_rows(
    data,
    climate_data_all %>%
      group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
      mutate(
        mean_TM = zoo::rollmean(TM, window, align = "right", fill = NA),
        mean_RR = zoo::rollmean(RR, window, align = "right", fill = NA),
        sd_TM = zoo::rollapply(TM, FUN = sd, width = window, align = "right", fill = NA),
        sd_RR = zoo::rollapply(RR, FUN = sd, width = window, align = "right", fill = NA),
        window = window
      ) %>%
      ungroup()
  )
}

climate_panel <- data %>% 
  full_join(means_all, by = c("GID_1", "gadm_level", "climate_source", "weight", "weight_year")) %>%
  full_join(means_pre, by = c("GID_1", "gadm_level", "climate_source", "weight", "weight_year")) 

data_with_econ <- econ_data_subnat %>%
  full_join(climate_panel, by = c("GID_1", "year")) %>%
  group_by(GID_1,econ_source,gadm_level,climate_source,window,weight,weight_year) %>%
  arrange(year,.by_group=TRUE) %>%
  group_by(GID_1,econ_source,gadm_level,climate_source,window,weight,weight_year) %>% 
  mutate(dlgrp_pc_usd = log(grp_pc_usd) - log(lag(grp_pc_usd)),
         dg = ( grp_pc_usd - lag(grp_pc_usd) ) / lag(grp_pc_usd), 
         lgrp_pc_usd=log(grp_pc_usd),
         lag_lgrp_pc_usd=log(lag(grp_pc_usd))) %>% 
  group_by(GID_1,econ_source,gadm_level,climate_source,window,weight,weight_year) %>%
  mutate(dev_TM_nosd = (TM - lag(mean_TM)),
         dev_RR_nosd = (RR - lag(mean_RR)),
         dev_TM_all_nosd = (TM - mean_TM_all),
         dev_RR_all_nosd = (RR - mean_RR_all),
         dev_TM = (TM - lag(mean_TM))/lag(sd_TM),
         dev_RR = (RR - lag(mean_RR))/lag(sd_RR),
         dev_TM_all = (TM - mean_TM_all)/sd_TM_all,
         dev_RR_all = (RR - mean_RR_all)/sd_RR_all,
         dTM = TM - lag(TM),
         dRR = RR - lag(RR)) %>%
  group_by(GID_1,econ_source,gadm_level,climate_source,window,weight,weight_year) %>%
  mutate(ldTM=lag(dTM),
         ldRR=lag(dRR),
         lTM=lag(TM),
         lRR=lag(RR)) %>% 
  ungroup() %>% 
  mutate(dev_TM_2 = dev_TM^2, 
         dev_RR_2 = dev_RR^2, 
         dev_TM_2_nosd = dev_TM_nosd^2, 
         dev_RR_2_nosd = dev_RR_nosd^2, 
         dRR_2 = dRR^2, 
         dTM_2=dTM^2,
         dev_TM_all_2 = dev_TM_all^2,
         dev_RR_all_2 = dev_RR_all^2,
         dev_TM_all_2_nosd = dev_TM_all_nosd^2,
         dev_RR_all_2_nosd = dev_RR_all_nosd^2,
         RR_2=RR^2,
         TM_2=TM^2) %>% 
  mutate(GID_0 = stringr::str_extract(GID_1, "^.{3}"))

if (REQUIRE_COMPLETE_YEAR_GID1) {
  expected_econ_sources <- 1
  expected_climate_datasets <- climate_panel %>%
    distinct(gadm_level, climate_source, weight, weight_year) %>%
    nrow()

  complete_gid1_year <- data_with_econ %>%
    group_by(GID_1, year) %>%
    summarise(
      n_econ_complete = n_distinct(econ_source[!is.na(dlgrp_pc_usd)]),
      n_climate_complete = n_distinct(
        paste(gadm_level, climate_source, weight, weight_year, sep = "||")[!is.na(TM) & !is.na(RR)]
      ),
      .groups = "drop"
    ) %>%
    filter(
      n_econ_complete == expected_econ_sources,
      n_climate_complete == expected_climate_datasets
    ) %>%
    select(GID_1, year)

  data_with_econ <- data_with_econ %>%
    inner_join(complete_gid1_year, by = c("GID_1", "year"))
}

# separate precipitation and temperature related variables 
data_rr <- data_with_econ %>% 
  select(econ_source,gadm_level,climate_source,window,weight,weight_year,GID_0,GID_1,year,grp_pc_usd,dlgrp_pc_usd,dg,RR,mean_RR,sd_RR,mean_RR_all,sd_RR_all,dRR,dev_RR,dev_RR_all,dev_RR_2,dev_RR_all_2,dev_RR_nosd,dev_RR_all_nosd,dev_RR_2_nosd,dev_RR_all_2_nosd,RR_2,dRR_2,lRR,ldRR) %>% 
  rename(window_rr=window)
data_tm <- data_with_econ %>% 
  select(econ_source,gadm_level,climate_source,window,weight,weight_year,GID_0,GID_1,year,grp_pc_usd,dlgrp_pc_usd,dg,TM,mean_TM,sd_TM,mean_TM_all,sd_TM_all,dTM,dev_TM,dev_TM_all,dev_TM_2,dev_TM_all_2,dev_TM_nosd,dev_TM_all_nosd,dev_TM_2_nosd,dev_TM_all_2_nosd,TM_2,dTM_2,lTM,ldTM) %>% 
  rename(window_tm=window)

sanitize_tag <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}
join_tag <- function(x) {
  paste(sanitize_tag(unique(as.character(x))), collapse = "-")
}
econ_tag <- join_tag(SELECT_ECON_SOURCE)
climate_tag <- paste(
  join_tag(SELECT_GADM_LEVEL),
  join_tag(SELECT_CLIMATE_SOURCE),
  join_tag(SELECT_WEIGHT),
  join_tag(SELECT_WEIGHT_YEAR),
  sep = "_"
)
data_rr_parquet <- file.path("econometrics", "data", paste0("data_rr_", econ_tag, "__", climate_tag, ".parquet"))
data_tm_parquet <- file.path("econometrics", "data", paste0("data_tm_", econ_tag, "__", climate_tag, ".parquet"))

arrow::write_parquet(data_rr, data_rr_parquet)
arrow::write_parquet(data_tm, data_tm_parquet)
