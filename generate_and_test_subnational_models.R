# flags
run_data <- FALSE
overwrite_existing <- FALSE
require(tidyverse)

# specify formulas and models
# fixed effects
i <- "year + GID_1[year]+ GID_1[year^2] " #"GID_1 + GID_0 + year + GID_1[year] + GID_1[year^2] + GID_0[year] + GID_0[year^2]"
# panel indexes
pan_id <- c("GID_1", "year")
# output variable 
o <- "dlgrp_pc_usd"
# specifications including moving averages
r1 <- c("bhmadasd"="TM + TM_2 + RR + RR_2 + ((TM - mean_TM)/(sd_TM))^2 + ((RR - mean_RR)/(sd_RR))^2", #burke+ada
        "bhmadasd2"="TM + TM_2 + RR + RR_2 + (TM - mean_TM)^2 + (RR - mean_RR)^2", #burke+ada
        "bhmadasd3"="I((TM - mean_TM_all)/(sd_TM_all)) + mean_TM_all : I((TM - mean_TM_all)/(sd_TM_all)) + ((TM - mean_TM)/(sd_TM))^2 + I((RR - mean_RR_all)/(sd_RR_all)) + mean_RR_all : I((RR - mean_RR_all)/(sd_RR_all)) + ((RR - mean_RR)/(sd_RR))^2", #burke+ada
        "bhmadasd4"="I((TM - mean_TM)/(sd_TM_all)) + mean_TM_all : I((TM - mean_TM)/(sd_TM_all)) + ((TM - mean_TM)/(sd_TM_all))^2 + I((RR - mean_RR)/(sd_RR_all)) + mean_RR_all : I((RR - mean_RR)/(sd_RR_all)) + ((RR - mean_RR)/(sd_RR_all))^2", #burke+ada
        "bhmadasd5"="I((TM - mean_TM)) + mean_TM : I((TM - mean_TM)) + ((TM - mean_TM)/(sd_TM))^2 + I((RR - mean_RR)) + mean_RR : I((RR - mean_RR)) + ((RR - mean_RR)/(sd_RR))^2", #burke+ada
        "bhmadaagr"="lag_lgrp_pc_usd + lag_share_ag_gdp*TM + lag_share_ag_gdp*RR + lag_share_ag_gdp*TM_2 + lag_share_ag_gdp*RR_2", #burke+ada
#        "bhmada"="TM + TM_2 + RR + RR_2 + dev_TM_2_nosd + dev_RR_2_nosd + dev_TM_nosd + dev_RR_nosd", #burke+ada
        "bhmadanosd"="TM + TM_2 +  dev_TM_2_nosd + dev_TM_2 + dev_TM_2_nosd:lag_ag_share + RR + RR_2 + dev_RR_2_nosd + dev_RR_2 + dev_RR_2_nosd:lag_ag_share", #burke+ada
        "adanosd"="dev_TM_2_nosd + dev_RR_2_nosd", #ada
        "adasd"="dev_TM_2 + dev_RR_2 + dev_TM + dev_RR",
#        "adasdfull"="dev_TM_2 + dev_RR_2 + dev_TM + dev_RR",
#        "ada"="dev_TM_2_nosd + dev_RR_2_nosd + dev_TM_nosd + dev_RR_nosd" #adanosd
        "kotz"="mean_TM:dTM  + mean_TM:ldTM + mean_TM:l2dTM + mean_TM:l3dTM + mean_TM:l4dTM + mean_TM:l5dTM + mean_TM:l6dTM + mean_TM:l7dTM + mean_TM:l8dTM + mean_TM:l9dTM + mean_TM:l10dTM + mean_RR:dRR + mean_RR:ldRR + mean_RR:l2dRR + mean_RR:l3dRR + mean_RR:l4dRR + mean_RR:l5dRR + mean_RR:l6dRR + mean_RR:l7dRR + mean_RR:l8dRR + mean_RR:l9dRR + mean_RR:l10dRR + ldRR + dRR + l2dRR + l3dRR + l4dRR + l5dRR + l6dRR + l7dRR + l8dRR + l9dRR + l10dRR + dTM + ldTM +l2dTM + l3dTM + l4dTM + l5dTM + l6dTM + l7dTM + l8dTM + l9dTM + l10dTM",#kotz
        "me"=paste(
          "mean_TM_all:dTMp + mean_TM_all:ldTMp + mean_TM_all:l2dTMp + mean_TM_all:l3dTMp + mean_TM_all:l4dTMp + mean_TM_all:l5dTMp + mean_TM_all:l6dTMp + mean_TM_all:l7dTMp + mean_TM_all:l8dTMp + mean_TM_all:l9dTMp + mean_TM_all:l10dTMp + dTMp + ldTMp + l2dTMp + l3dTMp + l4dTMp + l5dTMp + l6dTMp + l7dTMp + l8dTMp + l9dTMp + l10dTMp",
          "mean_TM_all:dTMm + mean_TM_all:ldTMm + mean_TM_all:l2dTMm + mean_TM_all:l3dTMm + mean_TM_all:l4dTMm + mean_TM_all:l5dTMm + mean_TM_all:l6dTMm + mean_TM_all:l7dTMm + mean_TM_all:l8dTMm + mean_TM_all:l9dTMm + mean_TM_all:l10dTMm + dTMm + ldTMm + l2dTMm + l3dTMm + l4dTMm + l5dTMm + l6dTMm + l7dTMm + l8dTMm + l9dTMm + l10dTMm",
          "mean_RR_all:dRRp + mean_RR_all:ldRRp + mean_RR_all:l2dRRp + mean_RR_all:l3dRRp + mean_RR_all:l4dRRp + mean_RR_all:l5dRRp + mean_RR_all:l6dRRp + mean_RR_all:l7dRRp + mean_RR_all:l8dRRp + mean_RR_all:l9dRRp + mean_RR_all:l10dRRp + dRRp + ldRRp + l2dRRp + l3dRRp + l4dRRp + l5dRRp + l6dRRp + l7dRRp + l8dRRp + l9dRRp + l10dRRp",
          "mean_RR_all:dRRm + mean_RR_all:ldRRm + mean_RR_all:l2dRRm + mean_RR_all:l3dRRm + mean_RR_all:l4dRRm + mean_RR_all:l5dRRm + mean_RR_all:l6dRRm + mean_RR_all:l7dRRm + mean_RR_all:l8dRRm + mean_RR_all:l9dRRm + mean_RR_all:l10dRRm + dRRm + ldRRm + l2dRRm + l3dRRm + l4dRRm + l5dRRm + l6dRRm + l7dRRm + l8dRRm + l9dRRm + l10dRRm",
          sep = " + "),#me
        "me2"=paste(
  "mean_TM_all:dTMp + mean_TM_all:ldTMp + mean_TM_all:l2dTMp + mean_TM_all:l3dTMp + mean_TM_all:l4dTMp + mean_TM_all:l5dTMp + mean_TM_all:l6dTMp + mean_TM_all:l7dTMp + mean_TM_all:l8dTMp + mean_TM_all:l9dTMp + mean_TM_all:l10dTMp + mean_TM_all:dTMp_2 + mean_TM_all:ldTMp_2 + mean_TM_all:l2dTMp_2 + mean_TM_all:l3dTMp_2 + mean_TM_all:l4dTMp_2 + mean_TM_all:l5dTMp_2 + mean_TM_all:l6dTMp_2 + mean_TM_all:l7dTMp_2 + mean_TM_all:l8dTMp_2 + mean_TM_all:l9dTMp_2 + mean_TM_all:l10dTMp_2 + dTMp + ldTMp + l2dTMp + l3dTMp + l4dTMp + l5dTMp + l6dTMp + l7dTMp + l8dTMp + l9dTMp + l10dTMp + dTMp_2 + ldTMp_2 + l2dTMp_2 + l3dTMp_2 + l4dTMp_2 + l5dTMp_2 + l6dTMp_2 + l7dTMp_2 + l8dTMp_2 + l9dTMp_2 + l10dTMp_2",
  "mean_TM_all:dTMm + mean_TM_all:ldTMm + mean_TM_all:l2dTMm + mean_TM_all:l3dTMm + mean_TM_all:l4dTMm + mean_TM_all:l5dTMm + mean_TM_all:l6dTMm + mean_TM_all:l7dTMm + mean_TM_all:l8dTMm + mean_TM_all:l9dTMm + mean_TM_all:l10dTMm + mean_TM_all:dTMm_2 + mean_TM_all:ldTMm_2 + mean_TM_all:l2dTMm_2 + mean_TM_all:l3dTMm_2 + mean_TM_all:l4dTMm_2 + mean_TM_all:l5dTMm_2 + mean_TM_all:l6dTMm_2 + mean_TM_all:l7dTMm_2 + mean_TM_all:l8dTMm_2 + mean_TM_all:l9dTMm_2 + mean_TM_all:l10dTMm_2 + dTMm + ldTMm + l2dTMm + l3dTMm + l4dTMm + l5dTMm + l6dTMm + l7dTMm + l8dTMm + l9dTMm + l10dTMm + dTMm_2 + ldTMm_2 + l2dTMm_2 + l3dTMm_2 + l4dTMm_2 + l5dTMm_2 + l6dTMm_2 + l7dTMm_2 + l8dTMm_2 + l9dTMm_2 + l10dTMm_2",
  "mean_RR_all:dRRp + mean_RR_all:ldRRp + mean_RR_all:l2dRRp + mean_RR_all:l3dRRp + mean_RR_all:l4dRRp + mean_RR_all:l5dRRp + mean_RR_all:l6dRRp + mean_RR_all:l7dRRp + mean_RR_all:l8dRRp + mean_RR_all:l9dRRp + mean_RR_all:l10dRRp + mean_RR_all:dRRp_2 + mean_RR_all:ldRRp_2 + mean_RR_all:l2dRRp_2 + mean_RR_all:l3dRRp_2 + mean_RR_all:l4dRRp_2 + mean_RR_all:l5dRRp_2 + mean_RR_all:l6dRRp_2 + mean_RR_all:l7dRRp_2 + mean_RR_all:l8dRRp_2 + mean_RR_all:l9dRRp_2 + mean_RR_all:l10dRRp_2 + dRRp + ldRRp + l2dRRp + l3dRRp + l4dRRp + l5dRRp + l6dRRp + l7dRRp + l8dRRp + l9dRRp + l10dRRp + dRRp_2 + ldRRp_2 + l2dRRp_2 + l3dRRp_2 + l4dRRp_2 + l5dRRp_2 + l6dRRp_2 + l7dRRp_2 + l8dRRp_2 + l9dRRp_2 + l10dRRp_2",
  "mean_RR_all:dRRm + mean_RR_all:ldRRm + mean_RR_all:l2dRRm + mean_RR_all:l3dRRm + mean_RR_all:l4dRRm + mean_RR_all:l5dRRm + mean_RR_all:l6dRRm + mean_RR_all:l7dRRm + mean_RR_all:l8dRRm + mean_RR_all:l9dRRm + mean_RR_all:l10dRRm + mean_RR_all:dRRm_2 + mean_RR_all:ldRRm_2 + mean_RR_all:l2dRRm_2 + mean_RR_all:l3dRRm_2 + mean_RR_all:l4dRRm_2 + mean_RR_all:l5dRRm_2 + mean_RR_all:l6dRRm_2 + mean_RR_all:l7dRRm_2 + mean_RR_all:l8dRRm_2 + mean_RR_all:l9dRRm_2 + mean_RR_all:l10dRRm_2 + dRRm + ldRRm + l2dRRm + l3dRRm + l4dRRm + l5dRRm + l6dRRm + l7dRRm + l8dRRm + l9dRRm + l10dRRm + dRRm_2 + ldRRm_2 + l2dRRm_2 + l3dRRm_2 + l4dRRm_2 + l5dRRm_2 + l6dRRm_2 + l7dRRm_2 + l8dRRm_2 + l9dRRm_2 + l10dRRm_2",
  sep = " + "),#me2
        "me3"=paste(
  "mean_TM_all:dTMp_2 + mean_TM_all:ldTMp_2 + mean_TM_all:l2dTMp_2 + mean_TM_all:l3dTMp_2 + mean_TM_all:l4dTMp_2 + mean_TM_all:l5dTMp_2 + mean_TM_all:l6dTMp_2 + mean_TM_all:l7dTMp_2 + mean_TM_all:l8dTMp_2 + mean_TM_all:l9dTMp_2 + mean_TM_all:l10dTMp_2 + dTMp_2 + ldTMp_2 + l2dTMp_2 + l3dTMp_2 + l4dTMp_2 + l5dTMp_2 + l6dTMp_2 + l7dTMp_2 + l8dTMp_2 + l9dTMp_2 + l10dTMp_2",
  "mean_TM_all:dTMm_2 + mean_TM_all:ldTMm_2 + mean_TM_all:l2dTMm_2 + mean_TM_all:l3dTMm_2 + mean_TM_all:l4dTMm_2 + mean_TM_all:l5dTMm_2 + mean_TM_all:l6dTMm_2 + mean_TM_all:l7dTMm_2 + mean_TM_all:l8dTMm_2 + mean_TM_all:l9dTMm_2 + mean_TM_all:l10dTMm_2 + dTMm_2 + ldTMm_2 + l2dTMm_2 + l3dTMm_2 + l4dTMm_2 + l5dTMm_2 + l6dTMm_2 + l7dTMm_2 + l8dTMm_2 + l9dTMm_2 + l10dTMm_2",
  "mean_RR_all:dRRp_2 + mean_RR_all:ldRRp_2 + mean_RR_all:l2dRRp_2 + mean_RR_all:l3dRRp_2 + mean_RR_all:l4dRRp_2 + mean_RR_all:l5dRRp_2 + mean_RR_all:l6dRRp_2 + mean_RR_all:l7dRRp_2 + mean_RR_all:l8dRRp_2 + mean_RR_all:l9dRRp_2 + mean_RR_all:l10dRRp_2 + dRRp_2 + ldRRp_2 + l2dRRp_2 + l3dRRp_2 + l4dRRp_2 + l5dRRp_2 + l6dRRp_2 + l7dRRp_2 + l8dRRp_2 + l9dRRp_2 + l10dRRp_2",
  "mean_RR_all:dRRm_2 + mean_RR_all:ldRRm_2 + mean_RR_all:l2dRRm_2 + mean_RR_all:l3dRRm_2 + mean_RR_all:l4dRRm_2 + mean_RR_all:l5dRRm_2 + mean_RR_all:l6dRRm_2 + mean_RR_all:l7dRRm_2 + mean_RR_all:l8dRRm_2 + mean_RR_all:l9dRRm_2 + mean_RR_all:l10dRRm_2 + dRRm_2 + ldRRm_2 + l2dRRm_2 + l3dRRm_2 + l4dRRm_2 + l5dRRm_2 + l6dRRm_2 + l7dRRm_2 + l8dRRm_2 + l9dRRm_2 + l10dRRm_2",
  sep = " + ")#me3
        #"all"="TM + TM_2 + dTM + ldTM + TM:dTM + TM:ldTM + dRR + ldRR + RR:dRR + RR:ldRR + RR + RR_2 + dev_TM_2_nosd + dev_RR_2_nosd",
        #"kwada"="lTM + dTM + ldTM + lTM:dTM + lTM:ldTM + dRR + ldRR + lRR:dRR + lRR:ldRR + lRR + dev_TM_2_nosd + dev_RR_2_nosd"
        )

#specifications excluding moving averages   
r2 <- c(
        "bhmspec"="TM + TM_2 + RR + RR_2 + ((TM - mean_TM_all))^2 + (RR - mean_RR_all)^2", #specialization+burke
        "bhmspec1"="TM + TM_2 + RR + RR_2 + ((TM - mean_TM_all)/(sd_TM_all))^2 + ((RR - mean_RR_all)/sd_RR_all)^2", #specialization+burke
        "bhmspec2"="TM + TM_2 + RR + RR_2 + (TM - mean_TM_all)^2 + ((RR - mean_RR_all)/sd_RR_all)^2", #specialization+burke
        "spec"="dev_TM_all_2_nosd + dev_TM_all_2 + dev_RR_all_2",
        "bhm"="TM + TM_2 + RR + RR_2" # burke 
        )

## select gdp source data among available datasets 
SELECT_ECON_SOURCE <- c("DOSE_V2_11","KUMMU25")
SELECT_GADM_LEVEL <- c("gadm1")      # gadm0, gadm1
SELECT_CLIMATE_SOURCE <- c("era5") # cru, dela, era5
SELECT_WEIGHT <- c("area","pop")
SELECT_WEIGHT_YEAR <-  c("2015","2000") #c("2000", "2015", "un", "concurrent", "unspecified")
years <- c(1960,1990)
windows <- c(10,20)

# generate or load data
sanitize_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
join_tag <- function(x) paste(sanitize_tag(unique(as.character(x))), collapse = "-")

if(run_data==T) {
  source("econometrics/prepare_econ_data.R")
  source("econometrics/prepare_climate_data.R")
}

econ_files <- list.files("econometrics/data", pattern = "^econ_processed_.*\\.parquet$", full.names = TRUE)
climate_rr_files <- list.files("econometrics/data", pattern = "^data_rr_.*\\.parquet$", full.names = TRUE)
climate_tm_files <- list.files("econometrics/data", pattern = "^data_tm_.*\\.parquet$", full.names = TRUE)

requested_pairs <- tidyr::expand_grid(
  climate_source = tolower(as.character(SELECT_CLIMATE_SOURCE)),
  weight = tolower(as.character(SELECT_WEIGHT)),
  weight_year = tolower(as.character(SELECT_WEIGHT_YEAR))
)

token_in_name <- function(name, token) {
  stringr::str_detect(
    tolower(name),
    paste0("(^|[^a-z0-9])", stringr::str_replace_all(token, "([^a-z0-9])", "\\\\\\1"), "([^a-z0-9]|$)")
  )
}

dash_count <- function(path) stringr::str_count(basename(path), fixed("-"))

file_matches_requested_gadm_level <- function(path) {
  n <- basename(path)
  any(purrr::map_lgl(tolower(as.character(SELECT_GADM_LEVEL)), ~ token_in_name(n, .x)))
}

pick_files_by_requested_pairs <- function(files) {
  if (length(files) == 0) return(character(0))
  files <- files[purrr::map_lgl(files, file_matches_requested_gadm_level)]
  if (length(files) == 0) return(character(0))

  pair_ids <- seq_len(nrow(requested_pairs))
  file_to_pairs <- purrr::map(
    files,
    function(path) {
      which(
        purrr::pmap_lgl(
          requested_pairs,
          function(climate_source, weight, weight_year) {
            token_in_name(basename(path), climate_source) &&
              token_in_name(basename(path), weight) &&
              token_in_name(basename(path), weight_year)
          }
        )
      )
    }
  )
  names(file_to_pairs) <- files

  uncovered <- pair_ids
  remaining <- files
  picked <- character(0)

  while (length(uncovered) > 0 && length(remaining) > 0) {
    coverage <- purrr::map_int(
      remaining,
      ~ length(intersect(file_to_pairs[[.x]], uncovered))
    )
    best_coverage <- max(coverage)
    if (best_coverage == 0) break

    candidates <- remaining[coverage == best_coverage]
    best <- candidates[order(dash_count(candidates), basename(candidates))][1]

    picked <- c(picked, best)
    uncovered <- setdiff(uncovered, file_to_pairs[[best]])
    remaining <- setdiff(remaining, best)
  }

  unique(picked)
}

climate_rr_files <- pick_files_by_requested_pairs(climate_rr_files)
climate_tm_files <- pick_files_by_requested_pairs(climate_tm_files)

if (length(econ_files) == 0 || length(climate_rr_files) == 0 || length(climate_tm_files) == 0) {
  stop(
    paste0(
      "Required processed parquet files not found.\n",
      "Expected at least one of each:\n",
      " - econometrics/data/econ_processed_*.parquet\n",
      " - econometrics/data/data_rr_*.parquet\n",
      " - econometrics/data/data_tm_*.parquet\n",
      "Also, climate files are pre-filtered by requested climate_source-weight_year-weight combinations."
    )
  )
}

econ_data <- purrr::map_dfr(econ_files, ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
  filter(econ_source %in% SELECT_ECON_SOURCE) %>%
  distinct()

data_rr <- purrr::map_dfr(climate_rr_files, ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
  filter(
    gadm_level %in% SELECT_GADM_LEVEL,
    climate_source %in% SELECT_CLIMATE_SOURCE,
    weight %in% SELECT_WEIGHT,
    weight_year %in% c(SELECT_WEIGHT_YEAR,"un")
  ) %>%
  distinct()

data_tm <- purrr::map_dfr(climate_tm_files, ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
  filter(
    gadm_level %in% SELECT_GADM_LEVEL,
    climate_source %in% SELECT_CLIMATE_SOURCE,
    weight %in% SELECT_WEIGHT,
    weight_year %in% c(SELECT_WEIGHT_YEAR,"un")
  ) %>%
  distinct()

if (nrow(data_rr) == 0 || nrow(data_tm) == 0) {
  stop("No rows in assembled modeling data after applying selectors.")
}

estimate_model = function(model="TM+TM_2+RR+RR_2",name="bhm",rr=10,tm=10,weight_p="pop",weight_t="pop",year_min=1960,year_max=2025,select_econ="DOSE_V2_11",select_climate="era5",select_weight_yr="2015",return_model=F,drop_outlier_vars=NULL,outlier_quantile=99.9,plot_outliers=FALSE,validate=FALSE,drop_countries=NULL,o="dlgrp_pc_usd",cluster_errors="GID_1") {

models_summary <- data.frame()
models_valid <- data.frame()

cluster_errors <- toupper(trimws(as.character(cluster_errors)[1]))
allowed_cluster_errors <- c("IID", "GID_1", "GID_0")
if (is.na(cluster_errors) || !cluster_errors %in% allowed_cluster_errors) {
  stop(paste0("cluster_errors must be one of: ", paste(allowed_cluster_errors, collapse = ", ")))
}
cluster_vcov <- if (cluster_errors == "IID") {
  "iid"
} else {
  stats::as.formula(paste0("~", cluster_errors))
}

if (!is.null(drop_countries) && length(drop_countries) > 0) {
  drop_countries <- unique(toupper(trimws(as.character(drop_countries))))
  drop_countries <- drop_countries[!is.na(drop_countries) & nzchar(drop_countries)]
}

outlier_prob <- suppressWarnings(as.numeric(outlier_quantile))
if (!is.null(drop_outlier_vars) && length(drop_outlier_vars) > 0) {
  drop_outlier_vars <- unique(as.character(drop_outlier_vars))
  if (is.na(outlier_prob)) stop("outlier_quantile must be numeric.")
  if (outlier_prob > 1) outlier_prob <- outlier_prob / 100
  if (outlier_prob <= 0.5 || outlier_prob >= 1) {
    stop("outlier_quantile must be in (0.5, 1) or (50, 100).")
  }
}
plot_outliers <- as.logical(plot_outliers)[1]
if (is.na(plot_outliers)) stop("plot_outliers must be TRUE or FALSE.")
if (plot_outliers && (is.null(drop_outlier_vars) || length(drop_outlier_vars) == 0)) {
  print("plot_outliers=TRUE ignored because drop_outlier_vars is empty.")
}

plot_outliers_density <- function(data_before, data_after, vars, model_id, quantile_prob) {
  plot_data <- dplyr::bind_rows(
    data_before %>% ungroup() %>% 
      dplyr::select(dplyr::all_of(vars)) %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "variable", values_to = "value") %>%
      dplyr::mutate(sample = "Before"),
    data_after %>% ungroup() %>% 
      dplyr::select(dplyr::all_of(vars)) %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "variable", values_to = "value") %>%
      dplyr::mutate(sample = "After")
  ) %>%
    dplyr::filter(!is.na(value))

  if (nrow(plot_data) == 0) {
    print("Outlier plot skipped: selected variables have no non-missing values.")
    return(invisible(NULL))
  }

  vars_with_variation <- plot_data %>%
    dplyr::group_by(sample, variable) %>%
    dplyr::summarise(n_unique = dplyr::n_distinct(value), .groups = "drop") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(has_variation = all(n_unique > 1), .groups = "drop") %>%
    dplyr::filter(has_variation) %>%
    dplyr::pull(variable)

  if (length(vars_with_variation) == 0) {
    print("Outlier plot skipped: no selected variable has enough variation before and after filtering.")
    return(invisible(NULL))
  }

  p <- plot_data %>%
    dplyr::filter(variable %in% vars_with_variation) %>%
    ggplot(aes(x = value, color = sample, fill = sample)) +
    geom_density(alpha = 0.2, linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~variable, scales = "free") +
    labs(
      title = paste0("Outlier Filtering Diagnostics: ", model_id),
      subtitle = paste0("Quantile threshold = ", round(quantile_prob * 100, 3), "% (symmetric tails)"),
      x = NULL,
      y = "Density"
    ) +
    theme_minimal()

  print(p)
  invisible(p)
}

plot_outlier_origin_histograms <- function(outlier_rows, model_id) {
  if (nrow(outlier_rows) == 0) {
    print("Outlier origin plot skipped: no outlier rows were identified.")
    return(invisible(NULL))
  }

  origin_counts <- dplyr::bind_rows(
    outlier_rows %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::count(variable, origin = as.character(year), name = "n") %>%
      dplyr::mutate(dimension = "Year"),
    outlier_rows %>%
      dplyr::filter(!is.na(GID_0)) %>%
      dplyr::count(variable, origin = as.character(GID_0), name = "n") %>%
      dplyr::mutate(dimension = "Country")
  ) %>%
    dplyr::group_by(variable, dimension) %>%
    dplyr::slice_max(order_by = n, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()

  if (nrow(origin_counts) == 0) {
    print("Outlier origin plot skipped: outlier rows have no year or GID_0 values.")
    return(invisible(NULL))
  }

  p <- origin_counts %>%
    dplyr::mutate(origin = stats::reorder(origin, n)) %>%
    ggplot(aes(x = origin, y = n, fill = dimension)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_grid(variable ~ dimension, scales = "free_y") +
    labs(
      title = paste0("Outlier Origins: ", model_id),
      subtitle = "Top 10 years and countries among rows flagged by each outlier variable",
      x = NULL,
      y = "Flagged rows"
    ) +
    theme_minimal()

  print(p)
  invisible(p)
}
  
print(paste("Checking model", name," with rr:",rr,",tm:",tm,weight_t,"weights for temperature,",weight_p,"weights for precipitations.","econ dataset", select_econ, "climate dataset", select_climate, "weighting year", select_weight_yr,  "Panel starting from",year_min, "standard errors", cluster_errors))
id <- paste0(
  "PREC",weight_p,
  "_TEMP",weight_t,
  "_MODEL",name,
  "_RR",rr,
  "_TM",tm,
  "_YM",year_min,
  "_ECON", select_econ,
  "_CLIM", select_climate,
  "_WY", select_weight_yr,
  "_VCOV", cluster_errors
)
if (!is.null(drop_outlier_vars) && length(drop_outlier_vars) > 0) {
  q_tag <- sub("\\.?0+$", "", sprintf("%.3f", outlier_prob * 100))
  id <- paste0(id, "_OUT", join_tag(drop_outlier_vars), "_Q", gsub("\\.", "p", q_tag))
}
if (!is.null(drop_countries) && length(drop_countries) > 0) {
  id <- paste0(id, "_DROPGID0", join_tag(drop_countries))
}
if(id %in% id_already_run)  {
  print("this model was already estimated. Skipping...") 
  return(NULL) }

data_m <- econ_data %>% 
  filter(year>=year_min & year<=year_max & econ_source==select_econ & !is.na(dlgrp_pc_usd)) %>% 
  inner_join(data_rr %>% 
  filter(year>=year_min & year<=year_max &
      window_rr==ifelse(is.numeric(rr),rr,unique(data_rr$window_rr)[1]) &
      weight==weight_p &
      climate_source==select_climate &
      (weight_year==select_weight_yr|weight_year=="un") & 
        !is.na(RR)) %>% 
  select(-window_rr,-weight,-weight_year)  ) %>%
  inner_join( data_tm %>% 
      filter(year>=year_min & year<=year_max &
          window_tm==ifelse(is.numeric(tm),tm,unique(data_tm$window_tm)[1]) &
          weight==weight_t &
          climate_source==select_climate &
          (weight_year==select_weight_yr|weight_year=="un") & 
            !is.na(TM)) %>% 
      select(-window_tm,-weight,-weight_year) ) %>% 
  group_by(GID_1) %>% 
  arrange(year, .by_group = TRUE) 
  
if(nrow(data_m)==0) {
  print("this climate-econ pair was empty (hint: maybe a gadm0-gadm1 pair?)") 
  return(NULL)}
if (!is.null(drop_countries) && length(drop_countries) > 0) {
  n_before_country_drop <- nrow(data_m)
  data_m <- data_m %>%
    dplyr::filter(!GID_0 %in% drop_countries)
  print(paste0(
    "Dropped ",
    n_before_country_drop - nrow(data_m),
    " observations from ",
    length(drop_countries),
    " GID_0 country entries: ",
    paste(drop_countries, collapse = ", "),
    "."
  ))
}
if(nrow(data_m)==0) {
  print("all rows were removed after country filtering")
  return(NULL)
}
if (!is.null(drop_outlier_vars) && length(drop_outlier_vars) > 0) {
  missing_outlier_vars <- setdiff(drop_outlier_vars, names(data_m))
  if (length(missing_outlier_vars) > 0) {
    stop(paste0("Outlier filtering variables not found in data: ", paste(missing_outlier_vars, collapse = ", ")))
  }
  
  non_numeric_outlier_vars <- drop_outlier_vars[!purrr::map_lgl(drop_outlier_vars, ~ is.numeric(data_m[[.x]]))]
  if (length(non_numeric_outlier_vars) > 0) {
    stop(paste0("Outlier filtering variables must be numeric: ", paste(non_numeric_outlier_vars, collapse = ", ")))
  }
  
  data_m_before_outlier_drop <- data_m
  n_before_outlier_drop <- nrow(data_m)
  lower_prob <- 1 - outlier_prob
  outlier_rows_by_var <- data.frame()
  
  for (var_name in drop_outlier_vars) {
    values_non_na <- data_m[[var_name]][!is.na(data_m[[var_name]])]
    if (length(values_non_na) == 0) next
    
    bounds <- stats::quantile(values_non_na, probs = c(lower_prob, outlier_prob), na.rm = TRUE, names = FALSE)
    keep_rows <- is.na(data_m[[var_name]]) | (data_m[[var_name]] >= bounds[1] & data_m[[var_name]] <= bounds[2])
    outlier_rows_by_var <- dplyr::bind_rows(
      outlier_rows_by_var,
      data_m[!keep_rows, , drop = FALSE] %>%
        ungroup() %>%
        dplyr::transmute(
          variable = var_name,
          year = year,
          GID_0 = GID_0,
          value = .data[[var_name]]
        )
    )
    data_m <- data_m[keep_rows, , drop = FALSE]
  }
  
  print(paste0(
    "Dropped ",
    n_before_outlier_drop - nrow(data_m),
    " outlier rows using ",
    paste(drop_outlier_vars, collapse = ", "),
    " at ",
    round(outlier_prob * 100, 3),
    "% quantile."
  ))
  
  if (plot_outliers) {
    plot_outliers_density(
      data_before = data_m_before_outlier_drop,
      data_after = data_m,
      vars = drop_outlier_vars,
      model_id = id,
      quantile_prob = outlier_prob
    )
    plot_outlier_origin_histograms(
      outlier_rows = outlier_rows_by_var,
      model_id = id
    )
  }
}

if(nrow(data_m)==0) {
  print("all rows were removed after outlier filtering")
  return(NULL)
}
f <- as.formula(paste( o, "~", model, "|" ,i ))
m <- fixest::feols(f, data_m, panel.id=pan_id, vcov=cluster_vcov)
ci_66 <- confint(m, vcov=cluster_vcov, level=0.66)
ci_95 <- confint(m, vcov=cluster_vcov)
p_values <- fixest::pvalue(m)

models_summary <- data.frame(id_model=id,
                                   cluster_errors=cluster_errors,
                                   coef=names(m$coefficients),
                                   best=unname(m$coefficients),
                                   mlo=unname(ci_66[[1]]),
                                   lo=unname(ci_95[[1]]),
                                   hi=unname(ci_95[[2]]),
                                   mhi=unname(ci_66[[2]]),
                                   se=unname(m$se),
                                   p=unname(p_values))

### cross validation

if(validate==T) {
  
### functions for performing panel cross-validation
cv_rmse <- function(actual, pred) {
  keep <- !is.na(pred)
  actual <- actual[keep]
  pred <- pred[keep]
  sqrt(sum((actual - pred)^2) / length(pred))
}

cross_validation_fixest_1year = function (data_dem, f, y){
  data_train<-data_dem[ data_dem$year <y, ]
  data_test<-data_dem[ data_dem$year ==y, ]
  pred<-predict_cv_ols(data_train, data_test, f)
  return (cv_rmse(data_test[[o]], pred))
  
}

# rows_id is the vector of the row numbers to be used as a test
cross_validation_fixest_pick = function (data_dem, f, rows_id){
  
  data_train<-data_dem[setdiff(seq(1,nrow(data_dem)),rows_id), ]
  data_test<-data_dem[rows_id, ]
  pred<-predict_cv_ols(data_train, data_test, f)
  return (cv_rmse(data_test[[o]], pred))
  
}

cross_validation_fixest_3year = function (data_dem, f, y){
  
  data_train<-data_dem[ data_dem$year <y, ]
  data_test<-data_dem[ data_dem$year ==y, ]
  pred<-predict_cv_ols(data_train, data_test, f)
  e1= cv_rmse(data_test[[o]], pred)
  
  data_train<-data_dem[ data_dem$year <(y-1), ]
  data_test<-data_dem[ data_dem$year >=(y-1), ]
  pred<-predict_cv_ols(data_train, data_test, f)
  e2= cv_rmse(data_test[[o]], pred)
  
  data_train<-data_dem[ data_dem$year <(y-2), ]
  data_test<-data_dem[ data_dem$year >=(y-2), ]
  pred<-predict_cv_ols(data_train, data_test, f)
  e3= cv_rmse(data_test[[o]], pred)
  
  return ((e1+e2+e3)/3)
}

cross_validation_fixest = function (data_dem, f){
  
  e<-list()
  years<-unique(data_dem$year)
  for (y in 1:length(years)){ #nmax useful if we have lags, corresponds to max n lags in our reg to avoid using first years 
    data_train<-data_dem[ data_dem$year !=years[y], ]
    data_test<-data_dem[ data_dem$year ==years[y], ]
    pred<-predict_cv_ols(data_train, data_test, f)
    e[[length(e)+1]]<-cv_rmse(data_test[[o]], pred)
    
  }
  err= sqrt(mean((unlist(e)^2), na.rm=TRUE))
  
  return(err)
}

fit_local_warming_rate <- function(year, temperature, window_radius = 2) {
  n <- length(year)
  vapply(seq_len(n), function(i) {
    idx <- seq.int(max(1, i - window_radius), min(n, i + window_radius))
    if (length(idx) < (2 * window_radius + 1)) {
      return(NA_real_)
    }
    coef(stats::lm(temperature[idx] ~ year[idx]))[["year[idx]"]] * 10
  }, numeric(1))
}

data_dem <- fixest::demean(as.formula(paste0(o , "+", model , "~", i )), data = data_m, na.rm=FALSE)
model_fcv <- setdiff(names(data_dem), o)
names(data_dem)[match(model_fcv, names(data_dem))] <- paste0("cv_x", seq_along(model_fcv))
fcv <- paste0("cv_x", seq_along(model_fcv))
predict_cv_ols <- function(data_train, data_test, x_vars) {
  keep_train <- stats::complete.cases(data_train[, c(o, x_vars), drop = FALSE])
  if (!any(keep_train)) return(rep(NA_real_, nrow(data_test)))
  fit <- stats::lm.fit(
    x = as.matrix(data_train[keep_train, x_vars, drop = FALSE]),
    y = data_train[[o]][keep_train]
  )
  as.vector(as.matrix(data_test[, x_vars, drop = FALSE]) %*% fit$coefficients)
}
data_dem$year<-data_m$year

climate_velocity_hist <- data_m %>% 
  select(GID_0,GID_1,year,TM,dTM) %>% 
  mutate(
    warming_run = dTM > 0,
    warming_run_id = cumsum(warming_run != lag(warming_run, default = FALSE))
  ) %>% ungroup() %>% 
  mutate(nrow=row_number())

year_3y <- climate_velocity_hist %>% 
  group_by(GID_1, warming_run_id) %>%
  filter(year==max(year) & warming_run==TRUE & n() >= 3)

year_4y <- climate_velocity_hist %>% 
  group_by(GID_1, warming_run_id) %>%
  filter(year==max(year) & warming_run==TRUE & n() >= 4) 

year_5y <- climate_velocity_hist %>% 
  group_by(GID_1, warming_run_id) %>%
  filter(year==max(year) & warming_run==TRUE & n() >= 5) 

year_6y <- climate_velocity_hist %>% 
  group_by(GID_1, warming_run_id) %>%
  filter(year==max(year) & warming_run==TRUE & n() >= 6) 

year_7y <- climate_velocity_hist %>% 
  group_by(GID_1, warming_run_id) %>%
  filter(year==max(year) & warming_run==TRUE & n() >= 7) 

year_6y <- anti_join(year_6y,year_7y)
year_5y <- anti_join(year_5y,year_6y)
year_4y <- anti_join(year_4y,year_5y)
year_3y <- anti_join(year_3y,year_4y)
temp_outliers <- data_m |> ungroup() |> mutate(nrow=row_number()) |> group_by(GID_1) |> filter(dev_TM_2==max(dev_TM_2))
prec_outliers <- data_m |> ungroup() |> mutate(nrow=row_number()) |> group_by(GID_1) |> filter(dev_RR_2==max(dev_RR_2))
}

models_valid <- data.frame(id_model=id,
                                 cluster_errors=cluster_errors,
                                 model=as.character(f)[3],
                                 r2=fixest::r2(m,"ar2"),
                                 wr2=fixest::r2(m,type='wr2'),
                                 loocv_all = ifelse(validate==T,cross_validation_fixest(data_dem=data_dem, f=fcv),NA),
                                 loocv_last = ifelse(validate==T,cross_validation_fixest_1year(data_dem=data_dem, f=fcv,y=2018),NA),
                                 loocv_3y = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=year_3y$nrow),NA),
                                 loocv_4y = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=year_4y$nrow),NA),
                                 loocv_5y = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem,f=fcv,rows_id=year_5y$nrow),NA),
                                 loocv_6y = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=year_6y$nrow),NA),
                                 loocv_7y = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=year_7y$nrow),NA),
                                 loocv_Tx = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=temp_outliers$nrow),NA),
                                 loocv_Px = ifelse(validate==T,cross_validation_fixest_pick(data_dem=data_dem, f=fcv,rows_id=prec_outliers$nrow),NA),
                                 aic=AIC(m),
                                 bic=BIC(m),
                                 nosign = paste0(nrow(models_summary %>% filter(id_model==id & p<=0.05)),"/",nrow(models_summary %>% filter(id_model==id))),
                                 topt = - m$coefficients["TM"]/(2*m$coefficients["TM_2"]) ) 
if(return_model==F) return(list(models_summary,models_valid)) else return(list(models_summary,models_valid,m,data_m))
}

plot_distributed_lag <- function(model,
                                 lag_vars = c("dTMp", "dTMm", "dRRp", "dRRm"),
                                 interaction_vars = c(dTMp = "mean_TM_all",
                                                      dTMm = "mean_TM_all",
                                                      dRRp = "mean_RR_all",
                                                      dRRm = "mean_RR_all"),
                                 data = NULL,
                                 region_id = "GID_1",
                                 conf_level = 0.95,
                                 percentiles = c(0.05, 0.5, 0.95),
                                 max_lag_search = 20,
                                 verbose = TRUE) {

  if (is.null(data)) {
    data <- tryCatch(eval(model$call$data), error = function(e) NULL)
    if (is.null(data)) stop("Could not retrieve model data; pass `data` explicitly.")
  }

  coefs <- stats::coef(model)
  V <- stats::vcov(model)
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # percentiles are computed across regions: one value per region_id (e.g. GID_1),
  # since interaction variables like mean_TM_all / mean_RR_all are region-level constants
  region_quantiles <- function(var) {
    if (!region_id %in% names(data)) {
      stop(paste0("region_id '", region_id, "' not found in data."))
    }
    df <- data %>% dplyr::ungroup() %>%
      dplyr::select(dplyr::all_of(c(region_id, var))) %>%
      dplyr::filter(!is.na(.data[[var]])) %>%
      dplyr::distinct(.data[[region_id]], .keep_all = TRUE)
    stats::quantile(df[[var]], probs = percentiles, na.rm = TRUE)
  }

  build_term_name <- function(base_var, lag) {
    if (lag == 0) base_var
    else if (lag == 1) paste0("l", base_var)
    else paste0("l", lag, base_var)
  }

  find_interaction_name <- function(tname, inter) {
    if (is.na(inter)) return(NA_character_)
    cand <- c(paste0(inter, ":", tname), paste0(tname, ":", inter))
    hit <- cand[cand %in% names(coefs)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  result <- purrr::map_dfr(lag_vars, function(v) {
    inter <- if (v %in% names(interaction_vars)) interaction_vars[[v]] else NA_character_
    if (is.null(inter) || is.na(inter) || !nzchar(inter)) inter <- NA_character_

    max_lag <- -1
    for (k in 0:max_lag_search) {
      tname <- build_term_name(v, k)
      iname <- find_interaction_name(tname, inter)
      if (tname %in% names(coefs) || !is.na(iname)) max_lag <- k
    }
    if (max_lag < 0) return(NULL)

    if (!is.na(inter) && inter %in% names(data)) {
      z_vals <- region_quantiles(inter)
      if (verbose) {
        message(sprintf("[%s] interaction = %s ; percentiles across %d regions: %s",
                        v, inter,
                        dplyr::n_distinct(data[[region_id]]),
                        paste(sprintf("%s=%.3f", names(z_vals), z_vals), collapse = ", ")))
      }
    } else {
      if (verbose) message(sprintf("[%s] no interaction variable found in data; using z = 0.", v))
      z_vals <- stats::setNames(rep(0, length(percentiles)),
                                paste0(round(100 * percentiles), "%"))
    }

    purrr::map_dfr(0:max_lag, function(k) {
      tname <- build_term_name(v, k)
      iname <- find_interaction_name(tname, inter)

      beta  <- if (tname %in% names(coefs)) unname(coefs[tname]) else 0
      gamma <- if (!is.na(iname)) unname(coefs[iname]) else 0
      var_b  <- if (tname %in% rownames(V)) V[tname, tname] else 0
      var_g  <- if (!is.na(iname) && iname %in% rownames(V)) V[iname, iname] else 0
      cov_bg <- if (tname %in% rownames(V) && !is.na(iname) && iname %in% rownames(V)) V[tname, iname] else 0

      # flip sign for negative-side variables (dTMm / dRRm — anything matching m after TM/RR)
      sign_flip <- if (grepl("(TM|RR)m(_|$)", v)) -1 else 1

      purrr::map_dfr(seq_along(z_vals), function(i) {
        z <- unname(z_vals[i])
        me <- sign_flip * (beta + gamma * z)
        se <- sqrt(pmax(var_b + z^2 * var_g + 2 * z * cov_bg, 0))
        data.frame(
          variable    = v,
          interaction = ifelse(is.na(inter), NA_character_, inter),
          lag         = k,
          percentile  = names(z_vals)[i],
          z           = z,
          me          = me,
          se          = se,
          lo          = me - z_crit * se,
          hi          = me + z_crit * se,
          stringsAsFactors = FALSE
        )
      })
    })
  })

  if (is.null(result) || nrow(result) == 0) stop("No matching lag coefficients found in the model.")

  # pretty facet labels: dTMp -> TM+, dTMm -> TM-, dRRp -> RR+, dRRm -> RR- (preserves _2 etc.)
  result$variable <- result$variable %>%
    sub("^d(TM|RR)p", "\\1+", .) %>%
    sub("^d(TM|RR)m", "\\1-", .)

  result$percentile <- factor(result$percentile, levels = unique(result$percentile))

  p <- ggplot(result, aes(x = lag, y = me, color = percentile, fill = percentile, group = percentile)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      x = "Lag (years)",
      y = "Marginal effect",
      title = "Distributed lag marginal effects",
      subtitle = paste0("Ribbons = ", round(100 * conf_level), "% CI; lines evaluated at percentiles of the interaction term"),
      color = "Interaction percentile",
      fill  = "Interaction percentile"
    ) +
    theme_minimal()

  list(plot = p, data = result)
}


models_summary_rds <- "econometrics/data/models_summary.rds"
models_valid_rds <- "econometrics/data/models_valid.rds"

if (file.exists(models_summary_rds) & overwrite_existing==FALSE) {
  models_summary <- readRDS(models_summary_rds)
} else {
  models_summary <- data.frame()
}

if (file.exists(models_valid_rds) & overwrite_existing==FALSE) {
  models_valid <- readRDS(models_valid_rds)
} else {
  models_valid <- data.frame()
}
if (ncol(models_summary) > 0 && !"cluster_errors" %in% names(models_summary)) {
  models_summary$cluster_errors <- NA_character_
}
if (ncol(models_valid) > 0 && !"cluster_errors" %in% names(models_valid)) {
  models_valid$cluster_errors <- NA_character_
}
id_already_run <- if ("id_model" %in% names(models_summary)) unique(models_summary$id_model) else character(0)

ndone <- length(id_already_run)
for (weight in SELECT_WEIGHT) {
for (tstart in years) {
for (econ in SELECT_ECON_SOURCE) {
for (cli in SELECT_CLIMATE_SOURCE) {
for (wt_yr in SELECT_WEIGHT_YEAR) {
  
for(name in names(r1)) {
for (window in windows) {
  model <- r1[[name]]
  out <- estimate_model(model,name,window,window,weight,weight,tstart,select_econ=econ,select_climate=cli,select_weight_yr=wt_yr)
  if (is.null(out)) {next} else {
    ndone<- ndone+1
    print(paste("estimated model",ndone)) }
  models_summary <- rbind(models_summary,out[[1]])
  models_valid <- rbind(models_valid,out[[2]])
  }}

for (name in names(r2)) {
  model <- r2[[name]]
  out <- estimate_model(model,name,NA,NA,weight,weight,tstart,select_econ=econ,select_climate=cli,select_weight_yr=wt_yr)
  if (is.null(out)) next else {
    ndone<- ndone+1
    print(paste("estimated model",ndone)) }
  models_summary <- rbind(models_summary,out[[1]])
  models_valid <- rbind(models_valid,out[[2]])
}
  
  }} } }}

saveRDS(models_summary, models_summary_rds)
saveRDS(models_valid, models_valid_rds)

id_expanded <- models_valid %>% 
  select(id_model,model) %>% 
  mutate(spec=stringr::str_extract(id_model,"(?<=_MODEL).+?(?=_)|(?<=_MODEL).*"),
         rr_weight=stringr::str_extract(id_model,"(?<=PREC).+?(?=_)"),
         tm_weight=stringr::str_extract(id_model,"(?<=_TEMP).+?(?=_)"),
         rr=stringr::str_extract(id_model,"(?<=_RR).+?(?=_)"),
         tm=stringr::str_extract(id_model,"(?<=_TM).+?(?=_)"),
         ym=stringr::str_extract(id_model,"(?<=_YM).+?(?=_)"),
         cl=stringr::str_extract(id_model,"(?<=_CLIM).+?(?=_WY)"),
         ec=stringr::str_extract(id_model,"(?<=_ECON).+?(?=_CLIM)"),
         wy=stringr::str_extract(id_model,"(?<=_WY).+?(?=(_VCOV|_OUT|_DROPGID0)|$)"),
         vcov=stringr::str_extract(id_model,"(?<=_VCOV).+?(?=(_OUT|_DROPGID0)|$)")) %>% 
  mutate(fix=ifelse(str_detect(spec,"ada"),"no","yes" ))

##### analyze model and test performance
ggplot(models_summary %>% 
         full_join(id_expanded),aes(x=coef)) + 
  geom_boxplot(data=. %>% pivot_longer(c(mlo,lo,hi,mhi,best),values_to="all"), aes(y=all), color="red") + 
  geom_boxplot(aes(y=best)) + 
  ggrepel::geom_text_repel(data=. %>% 
              group_by(coef,fix) %>% summarise(count=n(),value=median(best)),
            aes(y=value,label=count), 
            position=position_dodge(width=1.0)) + 
  facet_wrap(fix~.,)

## evaluate kotz coefficients
ggplot(models_summary %>% 
         full_join(id_expanded) %>% 
         filter(spec=="kotz") %>% 
         mutate(lag = str_extract(coef, "(?<=l)\\d+(?=d[A-Z]{2})"),
                coef=str_remove_all(coef, "(?<=l)\\d+(?=d[A-Z]{2})|_all|_nosd")) %>% 
         mutate(lag=case_when(is.na(lag) & str_detect(coef,"l")~"0",
                              is.na(lag) & !str_detect(coef,"l")~"1",
                              .default=lag), 
                coef=str_remove(coef,"l")) ,
       aes(x=ordered(lag,c("NA",seq(0,10)) )) ) + 
  geom_hline(yintercept=0) +
  geom_rect(data=. %>% pivot_longer(c(mlo_clust,lo_clust,hi_clust,mhi_clust,best),values_to="all") %>% 
              group_by(lag,coef,ec) %>% 
              summarise(min=min(all),max=max(all)),
            aes(width=0.5,ymin=min,ymax=max,linetype=ec), color="red",fill=NA,position="dodge") + 
  geom_rect(data=. %>% pivot_longer(c(mlo,lo,hi,mhi,best),values_to="all") %>% 
              group_by(lag,coef,ec) %>% 
              summarise(min=min(all),max=max(all)),
            aes(width=0.5,ymin=min,ymax=max,linetype=ec), color="blue",fill=NA,position="dodge") + 
  geom_boxplot(aes(y=best,linetype=ec),position="dodge") + 
  facet_wrap(coef~.,scales="free")

## evaluate bhm coefficients
ggplot(models_summary %>% 
         full_join(id_expanded) %>% 
         mutate(coef=str_remove_all(coef,"_all|_nosd")) %>% 
       filter(ec=="WB" & spec %in% c("bhmadasd","adasd","bhm")),
       aes(x=spec ) ) + 
  geom_rect(data=. %>% 
              pivot_longer(c(mlo_clust,lo_clust,hi_clust,mhi_clust,best),values_to="all") %>% 
              group_by(spec,coef,ec) %>% 
              summarise(min=min(all),max=max(all)),
            aes(width=0.5,ymin=min,ymax=max,linetype=ec), color="red",fill=NA,position="dodge") + 
  geom_rect(data=. %>% pivot_longer(c(mlo,lo,hi,mhi,best),values_to="all") %>% 
              group_by(spec,coef,ec) %>% 
              summarise(min=min(all),max=max(all)),
            aes(width=0.5,ymin=min,ymax=max,linetype=ec), color="blue",fill=NA,position="dodge") + 
  geom_boxplot(aes(y=best,linetype=ec),position="dodge") + 
  geom_hline(yintercept=0) +
  facet_wrap(.~coef,scales="free")


ggplot(models_summary %>% 
         full_join(id_expanded) %>% 
         mutate(coef=str_remove_all(coef,"_all|_nosd")) %>% 
         filter(spec %in% c("bhm","bhmadanosd","bhmspecnosd") )) + 
  geom_boxplot(data=. %>% pivot_longer(c(mlo,lo,hi,mhi,best),values_to="all") %>%
                 group_by(name,tm,rr,spec,rr_weight,tm_weight,ym,ec) %>%
                 summarise(optt=-all[coef=="TM"]/(2*all[coef=="TM_2"])) %>%
                 filter(name!="C"),
               aes(x=ec,y=optt,color=interaction(ordered(tm,c("NA","5","10","15","20","25")),spec) ),position="dodge",linetype=2) +
  geom_boxplot(data=. %>% group_by(tm,rr,spec,rr_weight,tm_weight,ym,ec) %>% 
                 summarise(optt=-best[coef=="TM"]/(2*best[coef=="TM_2"])), 
               aes(x=ec,y=optt,color=interaction(ordered(tm,c("NA","5","10","15","20","25")),spec) ),position="dodge",linetype=1) + 
  coord_cartesian(ylim=c(5,25)) 

ggplot(models_summary %>% filter(p<0.1) %>% 
         full_join(id_expanded) %>% 
         filter(spec %in% c("all"))) + 
  # geom_boxplot(data=. %>% pivot_longer(c(mlo,lo,hi,mhi,best),values_to="all") %>% 
  #                group_by(name,tm,rr,spec,rr_weight,tm_weight) %>% 
  #                summarise(optt=-all[coef=="TM"]/(2*all[coef=="TM_2"])) %>% 
  #                filter(name!="C"), 
  #              aes(x=ym,y=optt,color=tm),position="dodge",linetype=2) +
  geom_boxplot(aes(x=spec,y=best,color=ordered(tm,c("NA","5","10","15","20","25"))),position="dodge",linetype=1) + 
  geom_hline(yintercept=0) +
  facet_wrap(coef~.,scales="free")

id_short <- function(id_model) {
  id_model %>%
    str_replace("^PREC", "P") %>%
    str_replace("_TEMP", "_T") %>%
    str_replace("_MODEL", "_M") %>%
    str_replace("_ECON", "_EC") %>%
    str_replace("_CLIM", "_CL") %>%
    str_replace("_ECDOSE_V2_11", "_ECdose") %>%
    str_replace("_ECWB", "_ECwb") %>%
    str_replace("_ECKUMMU2018", "_ECkummu")
}

models_summary <- models_summary %>% unique()
#store all models into GDX for GAMS 
gdxtools::write.gdx("input/data/PNAS_models.gdx",
                    params= list(
                      model_coefs=models_summary %>% 
                        mutate(
                          id_model = id_short(id_model)
                        ) %>%
                        tidyr::pivot_longer(c(best,hi,mhi,mlo,lo),names_to="ci") %>%
                        select(id_model,coef,ci,value),
                      model_pvalues=models_summary %>% 
                        mutate(
                          id_model = id_short(id_model)
                        ) %>%
                        rename(value=p) %>%
                        select(id_model,coef,value) ) ) 

## save main specification table for paper
library(fixest)
library(modelsummary)
library(flextable)

# If you want the table as-is:
# modelsummary(lms, stars = TRUE, output = "table.docx") 

# if you want to customize the table with `flextable`
modelsummary(lms, stars = TRUE, output = "flextable") |> 
  autofit() |> 
  save_as_docx(path = "table.docx")

ada <- estimate_model(model=r1["adasd"],name="ADA",return_model=T,validate=T)
bhm <- estimate_model(model=r2["bhm"],name="BHM",rr=NA,tm=NA,return_model=T,validate=T)
bhmada <- estimate_model(model=r1["bhmadasd"],name="BHM+ADA",return_model=T,validate=T)
kotz <- estimate_model(model=r1["kotz"],name="KOTZ",return_model=T,validate=T)
models <- list("ADA"=summary(ada[[3]],vcov = ~GID_0),
  "BHM"=summary(bhm[[3]],vcov = ~GID_0),
  "BHMADA"=summary(bhmada[[3]],vcov = ~GID_0),
  "LAG"=summary(kotz[[3]],vcov = ~GID_0))
 
cross_valid <- data.frame(IMP=c("ADA","BHM","BHM+ADA","LAG"),
                          cv_all=c(ada[[2]]$loocv_all,bhm[[2]]$loocv_all,bhmada[[2]]$loocv_all,kotz[[2]]$loocv_all),
                          cv_2018=c(ada[[2]]$loocv_last,bhm[[2]]$loocv_last,bhmada[[2]]$loocv_last,kotz[[2]]$loocv_last),
                          cv_3yw=c(ada[[2]]$loocv_3y,bhm[[2]]$loocv_3y,bhmada[[2]]$loocv_3y,kotz[[2]]$loocv_3y),
                          cv_4yw=c(ada[[2]]$loocv_4y,bhm[[2]]$loocv_4y,bhmada[[2]]$loocv_4y,kotz[[2]]$loocv_4y),
                          cv_5yw=c(ada[[2]]$loocv_5y,bhm[[2]]$loocv_5y,bhmada[[2]]$loocv_5y,kotz[[2]]$loocv_5y),
                          cv_6yw=c(ada[[2]]$loocv_6y,bhm[[2]]$loocv_6y,bhmada[[2]]$loocv_6y,kotz[[2]]$loocv_6y),
                          cv_7yw=c(ada[[2]]$loocv_7y,bhm[[2]]$loocv_7y,bhmada[[2]]$loocv_7y,kotz[[2]]$loocv_7y))

ggplot(cross_valid %>% 
         pivot_longer(c(cv_3yw,cv_4yw,cv_5yw,cv_6yw,cv_7yw)) %>% mutate(name=str_remove_all(name,"cv_|w"),mean=mean(cv_all,na.rm=TRUE))) + 
  geom_bar(aes(x=name,y=100*(value-cv_all)/cv_all,fill=IMP),position="dodge",stat="identity") +
  xlab("year window") + ylab("% variation in cross validation performance rtm")


modelsummary(models, stars = TRUE, coef_omit = "RR", output = "flextable") |> 
  autofit() |> 
  save_as_docx(path = "table.docx")
