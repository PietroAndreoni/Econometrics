suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(zoo))

# Climate preprocessing selectors
SELECT_GADM_LEVEL <- c("gadm0")      # gadm0, gadm1
SELECT_CLIMATE_SOURCE <- c("era5") # cru, dela, era5
SELECT_WEIGHT <- c("pop", "area")
SELECT_WEIGHT_YEAR <-  c("2000","2015") #c("2000", "2015", "un", "concurrent", "unspecified")
WINDOWS <- c(5,10,20,30)
# Number of lagged regressors (p) in the Hamilton (2018) filter used for the
# temperature/precipitation trends. The forecast horizon h is set per-series to
# the same value as WINDOWS; p is Hamilton's suggested default.
HAMILTON_LAGS <- 4
# Reference periods (inclusive year bounds) used both for the fixed-period means
# (means_all / means_pre) and for the period-trained Hamilton trends
# (trend_*_all / trend_*_pre).
PERIOD_ALL <- c(1990, 2019)
PERIOD_PRE <- c(1960, 1989)
# Horizon (h) for the fixed-period Hamilton trends (trend_*_all / trend_*_pre).
# These are trained only on the reference periods above, so - unlike the rolling
# panel trend, whose horizon is the WINDOWS width - their horizon is a single
# fixed value and the trends do not depend on WINDOWS.
HAMILTON_H_FIXED <- 30

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

  data_long <- arrow::read_parquet(path) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(-Date, names_to = "GID_1", values_to = "value") %>%
    filter(!is.na(GID_1), !is.na(value)) %>%
    # Some CRU files include a placeholder "?" column: drop it early.
    filter(GID_1 != "?")

  # Use a scalar branch here; vectorized ifelse() with length-1 condition
  # would collapse GID_1 to one value.
  if (meta$gadm_level == "gadm1") {
    data_long <- data_long %>%
      mutate(GID_1 = sub("^([^_]+)_", "\\1.", GID_1)) %>%
      # Keep only expected GADM1-like identifiers (e.g., PHL.110_1).
      filter(grepl("^[A-Z0-9]{3}\\..+", GID_1))
  } else {
    data_long <- data_long %>%
      mutate(GID_1 = as.character(GID_1)) %>%
      # For GADM0, IDs are ISO3-like codes.
      filter(grepl("^[A-Z0-9]{3}$", GID_1))
  }

  data_long %>%
    mutate(year = as.integer(substr(gsub("^X", "", as.character(Date)), 1, 4))) %>%
    group_by(GID_1, year) %>%
    summarise(value = yearly_fn(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      variable = meta$variable,
      weight = meta$weight,
      weight_year = meta$weight_year,
      climate_source = meta$climate_source,
      gadm_level = meta$gadm_level
    )
}

add_lag_columns <- function(data, vars, n_lags, prefix = "l") {
  lag_grid <- purrr::map_dfr(
    seq_len(n_lags),
    ~ tibble::tibble(lag = .x, var = vars)
  )
  lag_names <- ifelse(
    lag_grid$lag == 1,
    paste0(prefix, lag_grid$var),
    paste0(prefix, lag_grid$lag, lag_grid$var)
  )
  lag_exprs <- purrr::map2(
    lag_grid$var,
    lag_grid$lag,
    ~ rlang::expr(dplyr::lag(!!rlang::sym(.x), !!.y))
  )

  dplyr::mutate(data, !!!rlang::set_names(lag_exprs, lag_names))
}

add_power_columns <- function(data, vars, power) {
  dplyr::mutate(
    data,
    dplyr::across(
      dplyr::all_of(vars),
      ~ .x^power,
      .names = paste0("{.col}_", power)
    )
  )
}

# Hamilton (2018) regression filter. For a single series y ordered by year, the
# trend is the fitted value of
#   y_t = b0 + b1*y_{t-h} + b2*y_{t-h-1} + ... + bp*y_{t-h-p+1} + v_t
# estimated by OLS; the residual v_t is the cyclical part. The horizon h is set
# to the same window used for the rolling means.
#
# `train` is an optional logical vector (same length as y) selecting the rows on
# which the OLS coefficients are estimated; the lagged predictors are always
# built from the full series, so training rows can reach back before the training
# window for their lags. When `train` is NULL the whole series is used and the
# trend is returned for every row with available predictors; otherwise the trend
# is returned only for the training rows (the trend "over" that period). Returns
# a same-length vector, NA where predictors are unavailable or too few rows are
# available to estimate the regression.
hamilton_trend <- function(y, h, p = 4L, train = NULL) {
  h <- as.integer(h)
  p <- as.integer(p)
  n <- length(y)
  trend <- rep(NA_real_, n)
  if (is.null(train)) {
    train <- rep(TRUE, n)
  }
  train <- train & !is.na(train)
  if (n <= h + p) {
    return(trend)
  }
  lag_cols <- lapply(seq.int(0L, p - 1L), function(j) dplyr::lag(y, h + j))
  X <- as.data.frame(stats::setNames(lag_cols, paste0("lag", seq.int(h, h + p - 1L))))
  df <- cbind(data.frame(y = y), X)
  has_predictors <- stats::complete.cases(X)
  fit_rows <- train & stats::complete.cases(df)
  if (sum(fit_rows) <= p + 1L) {
    return(trend)
  }
  fit <- stats::lm(y ~ ., data = df[fit_rows, , drop = FALSE])
  out_rows <- train & has_predictors
  trend[out_rows] <- as.numeric(
    stats::predict(fit, newdata = df[out_rows, , drop = FALSE])
  )
  trend
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

selected_climate_meta <- climate_file_meta %>%
  # Area-weighted files are encoded as "un" with no vintage year.
  # If area is requested, always include "un" in allowed weight_year values.
  {
    allowed_weight_year <- unique(c(
      as.character(SELECT_WEIGHT_YEAR),
      if ("area" %in% as.character(SELECT_WEIGHT)) "un"
    ))
    filter(
      .,
      gadm_level %in% SELECT_GADM_LEVEL,
      climate_source %in% SELECT_CLIMATE_SOURCE,
      weight %in% SELECT_WEIGHT,
      weight_year %in% allowed_weight_year
    )
  }

if (nrow(selected_climate_meta) == 0) {
  stop("No climate files match current selectors.")
}

selected_balance <- selected_climate_meta %>%
  count(gadm_level, climate_source, weight, weight_year, variable, name = "n_files") %>%
  tidyr::pivot_wider(names_from = variable, values_from = n_files, values_fill = 0)

selected_unbalanced <- selected_balance %>%
  filter(tmp == 0 | pre == 0)

if (nrow(selected_unbalanced) > 0) {
  print(selected_unbalanced)
  stop("Selected climate subset is unbalanced: tmp or pre missing for at least one tuple.")
}

climate_data_long <- purrr::map_dfr(selected_climate_meta$path, aggregate_climate_parquet)

climate_data_yearly <- climate_data_long %>%
  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
  rename(TM = tmp, RR = pre) %>%
  mutate(RR = RR / 1000) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  filter(any(dplyr::coalesce(TM, 0) != 0 | dplyr::coalesce(RR, 0) != 0)) %>%
  arrange(year, .by_group = TRUE) %>%
  ungroup()

means_all <- climate_data_yearly %>%
  filter(year >= PERIOD_ALL[1] & year <= PERIOD_ALL[2]) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  summarise(
    sd_TM_all = sd(TM, na.rm = TRUE),
    sd_RR_all = sd(RR, na.rm = TRUE),
    mean_TM_all = mean(TM, na.rm = TRUE),
    mean_RR_all = mean(RR, na.rm = TRUE),
    .groups = "drop"
  )

means_pre <- climate_data_yearly %>%
  filter(year >= PERIOD_PRE[1] & year <= PERIOD_PRE[2]) %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  summarise(
    sd_TM_pre = sd(TM, na.rm = TRUE),
    sd_RR_pre = sd(RR, na.rm = TRUE),
    mean_TM_pre = mean(TM, na.rm = TRUE),
    mean_RR_pre = mean(RR, na.rm = TRUE),
    .groups = "drop"
  )

climate_panel <- data.frame()
for (window in WINDOWS) {
  climate_panel <- bind_rows(
    climate_panel,
    climate_data_yearly %>%
      group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
      mutate(
        mean_TM = zoo::rollmean(TM, window, align = "right", fill = NA),
        mean_RR = zoo::rollmean(RR, window, align = "right", fill = NA),
        sd_TM = zoo::rollapply(TM, FUN = sd, width = window, align = "right", fill = NA),
        sd_RR = zoo::rollapply(RR, FUN = sd, width = window, align = "right", fill = NA),
        trend_TM = hamilton_trend(TM, h = window, p = HAMILTON_LAGS),
        trend_RR = hamilton_trend(RR, h = window, p = HAMILTON_LAGS),
        window = window) %>%
      ungroup()
  )
}

# Period-trained Hamilton trends. Unlike the rolling-window trend in the panel,
# these are estimated once per region on each fixed reference period (horizon
# HAMILTON_H_FIXED, independent of WINDOWS) and broadcast onto every window when
# joined onto the panel by region-year below.
trends_fixed <- climate_data_yearly %>%
  group_by(GID_1, gadm_level, climate_source, weight, weight_year) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    trend_TM_all = hamilton_trend(TM, h = HAMILTON_H_FIXED, p = HAMILTON_LAGS,
                                  train = year >= PERIOD_ALL[1] & year <= PERIOD_ALL[2]),
    trend_RR_all = hamilton_trend(RR, h = HAMILTON_H_FIXED, p = HAMILTON_LAGS,
                                  train = year >= PERIOD_ALL[1] & year <= PERIOD_ALL[2]),
    trend_TM_pre = hamilton_trend(TM, h = HAMILTON_H_FIXED, p = HAMILTON_LAGS,
                                  train = year >= PERIOD_PRE[1] & year <= PERIOD_PRE[2]),
    trend_RR_pre = hamilton_trend(RR, h = HAMILTON_H_FIXED, p = HAMILTON_LAGS,
                                  train = year >= PERIOD_PRE[1] & year <= PERIOD_PRE[2])
  ) %>%
  ungroup() %>%
  select(GID_1, gadm_level, climate_source, weight, weight_year, year,
         trend_TM_all, trend_RR_all, trend_TM_pre, trend_RR_pre)

climate_panel <- climate_panel %>%
  full_join(means_all, by = c("GID_1", "gadm_level", "climate_source", "weight", "weight_year")) %>%
  full_join(means_pre, by = c("GID_1", "gadm_level", "climate_source", "weight", "weight_year")) %>%
  left_join(trends_fixed, by = c("GID_1", "gadm_level", "climate_source", "weight", "weight_year", "year")) %>%
  group_by(GID_1,gadm_level,climate_source,window,weight,weight_year) %>%
  arrange(year,.by_group=TRUE) %>%
  group_by(GID_1,gadm_level,climate_source,window,weight,weight_year) %>%
  mutate(dev_TM_nosd = (TM - lag(mean_TM)),
         dev_RR_nosd = (RR - lag(mean_RR)),
         dev_TM_all_nosd = (TM - mean_TM_all),
         dev_RR_all_nosd = (RR - mean_RR_all),
         dev_TM = (TM - lag(mean_TM))/sd_TM_all,
         dev_RR = (RR - lag(mean_RR))/sd_RR_all,
         dev_TM_all = (TM - mean_TM_all)/(sd_TM_all),
         dev_RR_all = (RR - mean_RR_all)/(sd_RR_all),
         dTM = TM - lag(TM),
         dRR = RR - lag(RR),
         dTMp = pmax(TM - lag(TM), 0),
         dTMm = pmin(TM - lag(TM), 0),
         dRRp = pmax(RR - lag(RR), 0),
         dRRm = pmin(RR - lag(RR), 0) ) %>%
  ungroup() %>%  
  add_power_columns(
    c("dev_TM",
      "dev_RR",
      "dev_TM_nosd",
      "dev_RR_nosd",
      "dRR",
      "dTM",
      "dTMp",
      "dTMm",
      "dRRp",
      "dRRm",
      "dev_TM_all",
      "dev_RR_all",
      "dev_TM_all_nosd",
      "dev_RR_all_nosd",
      "RR",
      "TM"),2) %>%
  group_by(GID_1,gadm_level,climate_source,window,weight,weight_year) %>%
  add_lag_columns(c("TM", 
                    "RR",
                    "dTM", 
                    "dRR",
                    "dTMp",
                    "dTMm",
                    "dRRp",
                    "dRRm",
                    "dTMp_2",
                    "dTMm_2",
                    "dRRp_2",
                    "dRRm_2",
                    "dev_TM_all", 
                    "dev_RR_all",
                    "dev_TM_all_2",
                    "dev_RR_all_2"), 10) %>%
  ungroup() %>%  
  mutate(GID_0 = stringr::str_extract(GID_1, "^.{3}"))

sanitize_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
join_tag <- function(x) paste(sanitize_tag(unique(as.character(x))), collapse = "-")

climate_tag <- paste(
  join_tag(SELECT_GADM_LEVEL),
  join_tag(SELECT_CLIMATE_SOURCE),
  join_tag(SELECT_WEIGHT),
  join_tag(SELECT_WEIGHT_YEAR),
  sep = "_"
)

col_select <- c("gadm_level", "climate_source", "window", "weight", "weight_year", "GID_0", "GID_1", "year")
# separate precipitation and temperature related variables 
data_rr <- climate_panel %>% 
  select(all_of(col_select), contains("RR")) %>% 
  rename(window_rr=window)
data_tm <- climate_panel %>% 
  select(all_of(col_select), contains("TM")) %>% 
  rename(window_tm=window)

arrow::write_parquet(data_rr, file.path("econometrics", "data", paste0("data_rr_", climate_tag, ".parquet")))
arrow::write_parquet(data_tm, file.path("econometrics", "data", paste0("data_tm_", climate_tag, ".parquet")))

cat("Wrote:\n", file.path("econometrics", "data", paste0("data_rr_", climate_tag, ".parquet")), "\n", file.path("econometrics", "data", paste0("data_tm_", climate_tag, ".parquet")), "\n", sep = "")
