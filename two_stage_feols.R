
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
             weight=="pop" &
             climate_source=="era5" &
             (weight_year=="2000"|weight_year=="un") &
             !is.na(RR)) %>%
    select(-window_rr,-weight,-weight_year) %>%
    inner_join( data_tm %>%
                  filter(window_tm==w_tm &
                           weight=="pop" &
                           climate_source=="era5" &
                           (weight_year=="2000"|weight_year=="un") &
                           !is.na(TM)) %>%
                  select(-window_tm,-weight,-weight_year) ) %>%
    left_join(econ) %>%
    group_by(GID_0,GID_1) %>%
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

df <- build_dat(30,30)

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
df <- df %>%
  group_by(GID_0) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    y_lag  = lag_by_year(y, year),
    zT_lag = lag_by_year(zT, year),
    zP_lag = lag_by_year(zP, year),
    
    dy  = y  - y_lag,
    dzT = zT - zT_lag,
    dzP = zP - zP_lag
  ) %>%
  ungroup() %>%
  filter(
    complete.cases(
      GID_0, year,
      y, k,
      mean_TM,mean_RR,
      zT, zP
    )
  )

m_lr <- feols(
  y ~
    k + 
    mean_TM +
    mean_TM^2 +
    mean_RR +
    mean_RR^2 +
    mean_TM:mean_RR
  | GID_0 + year,
  data = df,
  panel.id = c("GID_0", "year"),
  cluster = ~GID_0
)

summary(m_lr)

df$V <- resid(m_lr)

m_ecm <- feols(
  dy ~ 
      zT2 +  
      zP2 + 
      zT2:zP2 + 
      I(V-l(V,1:1))
  | GID_0 + year,
  data = df,
  panel.id = c("GID_0", "year"),
  cluster = ~GID_0
)

summary(m_ecm)

etable(
  m_lr,
  m_ecm,
  headers = c(
    "Long-run level",
    "Tol ECM"
  )
)
