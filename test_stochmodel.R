
# A row is estimable when it has an outcome; everything else is a lag carrier
est_row <- function(d) !is.na(d$dlgrp_pc_usd)

# Value of x in year t-k, matched on the year itself rather than on row order,
# so gaps in a region's panel produce NA instead of silently shifting the lag.
lag_by_year <- function(x, year, k = 1) x[match(year - k, year)]

build_dat <- function(w_rr, w_tm) {
  
  econ_data <- purrr::map_dfr("econometrics/data/econ_processed_WDI-WB-PWT110.parquet", ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
    filter(econ_source %in% c("PWT110")) %>%
    distinct()
  
  data_rr <- purrr::map_dfr("econometrics/data/data_rr_gadm0_era5_pop-area_2000-2015.parquet", ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
    filter(
      gadm_level=="gadm0",
      climate_source=="era5",
      weight %in% c("area","pop"),
      weight_year %in% c("2000","un")
    ) %>%
    distinct()
  
  data_tm <- purrr::map_dfr("econometrics/data/data_tm_gadm0_era5_pop-area_2000-2015.parquet", ~ arrow::read_parquet(.x) %>% as_tibble()) %>%
    filter(
      gadm_level=="gadm0",
      climate_source=="era5",
      weight %in% c("area","pop"),
      weight_year %in% c("2000","un")
    ) %>%
    distinct()
  
  # Econ is restricted to the estimation window here; the climate panel below is
  # NOT, so the pre-ECON_YEAR_MIN years survive the join as outcome-less rows.
  econ <- econ_data %>%
    filter(econ_source == "PWT110" &
             !is.na(dlgrp_pc_usd))
  
  data_rr %>%
    filter(window_rr==w_rr &
             weight=="area" &
             climate_source=="era5" &
             (weight_year=="2015"|weight_year=="un") &
             !is.na(RR)) %>%
    select(-window_rr,-weight,-weight_year) %>%
    inner_join( data_tm %>%
                  filter(window_tm==w_tm &
                           weight=="area" &
                           climate_source=="era5" &
                           (weight_year=="2015"|weight_year=="un") &
                           !is.na(TM)) %>%
                  select(-window_tm,-weight,-weight_year) ) %>%
    left_join(econ) %>%
    group_by(GID_1) %>%
    # Regions with nothing to estimate are pure ballast, and climate running past
    # a region's last usable year carries no lag anyone reads
    filter(any(!is.na(dlgrp_pc_usd))) %>%
    filter(year <= max(year[!is.na(dlgrp_pc_usd)])) %>%
    mutate(log_gdp_av=log(mean(grp_pc_usd,na.rm=TRUE))) %>% 
    arrange(year, .by_group = TRUE) %>%
    # Centring of the deviation terms: the rolling window ending at t-1, so the
    # current year never enters its own reference mean. Built here, before the
    # common-sample restriction, so the warm-up year that loses its lag is
    # dropped by the completeness check below like any other missing moment.
    mutate(
      mean_TM_lag = lag_by_year(mean_TM, year),
      mean_RR_lag = lag_by_year(mean_RR, year)
    )
}

df <- build_dat(20,20)
  
df$zT <- (df$TM - df$mean_TM_lag)/df$sd_TM_all
df$zT2 <- ((df$TM - df$mean_TM_lag)/df$sd_TM_all)^2

df$hot  <- pmax(df$zT,  0)
df$cold <- pmax(-df$zT, 0)

df$zP <- (df$RR - df$mean_RR_lag)/df$sd_RR_all
df$zP2 <- ((df$RR - df$mean_RR_lag)/df$sd_RR_all)^2
df$wet  <- pmax(df$zP,  0)
df$dry  <- pmax(-df$zP, 0)
df$k <- log(df$k_pc_usd/df$share_emp_pop)
df$y <- log(df$grp_pc_usd/df$share_emp_pop)
vars <- c(
  "y", "k",
  "mean_TM_lag", "mean_RR_lag",
  "zT2", "zP2",
  "GID_0", "year"
)

df_sfa <- df %>%
  filter(complete.cases(across(all_of(vars)))) %>%
  filter(
    is.finite(y),
    is.finite(k),
    is.finite(mean_TM_lag),
    is.finite(mean_RR_lag),
    is.finite(zT2),
    is.finite(zP2)
  ) %>%
  mutate(
    country_fe = factor(GID_0),
    year_fe = factor(year)
  )

table(df_sfa$year_fe)
table(df_sfa$country_fe)

ols_check <- lm(
  y ~
    k +
    mean_TM_lag +
    I(mean_TM_lag^2) +
    mean_RR_lag +
    I(mean_RR_lag^2) +
    country_fe +
    year_fe,
  data = df_sfa
)

coef(ols_check)[is.na(coef(ols_check))]

require(frontier)
mod <- sfa(
  y ~
    k +
    mean_TM +
    I(mean_TM^2) +
    mean_RR +
    I(mean_RR^2) + 
    year_fe + country_fe
  |
    zT2 + zP2,
  data = df_sfa
)
summary(mod)

library(sfa)
library(plm)

df_sfa <- pdata.frame(
  df_sfa,
  index = c("GID_0", "year")
)

m_tfe <- sfa::psfm(
  y ~
    k +
    mean_TM +
    I(mean_TM^2) +
    mean_RR +
    I(mean_RR^2),
  data = df_sfa,
  individual = "GID_0",
  model_name = "TFE"
)


library(fixest)

# 1. very simple production relationship
o1 <- feols(
  y ~ k,
  data = df
)

# 2. climate frontier
o2 <- feols(
  y ~ k +
    mean_TM_lag + I(mean_TM_lag^2) +
    mean_RR_lag + I(mean_RR_lag^2),
  data = df
)

# 3. region FE
o3 <- feols(
  y ~ k +
    mean_TM_lag + I(mean_TM_lag^2) +
    mean_RR_lag + I(mean_RR_lag^2)
  | GID_1,
  data = df
)

# 4. region + year FE
o4 <- feols(
  y ~ k +
    mean_TM_lag + I(mean_TM_lag^2) +
    mean_RR_lag + I(mean_RR_lag^2)
  | GID_1 + year,
  data = df
)


library(e1071)

sapply(
  list(o1, o2, o3, o4),
  \(m) skewness(resid(m), na.rm = TRUE, type = 2)
)

# test 2
ols <- lm(
  y ~
    k +
    mean_TM_lag +
    I(mean_TM_lag^2) +
    mean_RR_lag +
    I(mean_RR_lag^2),
  data = df
)

hist(resid(ols), breaks = 100)
e1071::skewness(resid(ols), type = 2)

quantile(
  resid(ols),
  c(.001, .01, .05, .5, .95, .99, .999),
  na.rm = TRUE
)

ols_fe <- feols(
  y ~
    k +
    mean_TM_lag +
    I(mean_TM_lag^2) +
    mean_RR_lag +
    I(mean_RR_lag^2)
  | GID_1 + year,
  data = df
)

hist(resid(ols_fe), breaks = 100)
e1071::skewness(resid(ols_fe), type = 2)
