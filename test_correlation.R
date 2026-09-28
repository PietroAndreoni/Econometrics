library(dplyr)

region_dat <- dat %>%
  group_by(GID_1) %>%
  summarise(
    country = first(GID_0),
    
    sd_T = first(sd_TM_all),
    sd_R = first(sd_RR_all),
    
    mean_T = first(mean_TM_all),
    mean_R = first(mean_RR_all),
    
    # Example historical socioeconomic measures
    gdp_hist = mean(grp_pc_usd[year >= 1990 & year <= 2010],
                    na.rm = TRUE),
    
    agr_share_hist = mean(
      share_ag_gdp[year >= 1990 & year <= 2010],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

region_dat <- region_dat %>%
  mutate(
    log_sd_T = log(sd_T),
    log_sd_R = log(sd_R),
    log_gdp = log(gdp_hist)
  )

cor(
  region_dat[, c(
    "log_sd_T",
    "log_sd_R",
    "log_gdp",
    "agr_share_hist",
    "mean_T",
    "mean_R"
  )],
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor(
  region_dat[, c(
    "log_sd_T",
    "log_sd_R",
    "log_gdp",
    "agr_share_hist",
    "mean_T",
    "mean_R"
  )],
  use = "pairwise.complete.obs",
  method = "spearman"
)
library(ggplot2)

ggplot(region_dat, aes(log_gdp, log_sd_T)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(
    x = "Log historical GDP per capita",
    y = "Log historical temperature SD"
  )
library(fixest)

m_sdT_country <- feols(
  log_sd_T ~
    mean_T +
    I(mean_T^2) +
    log_gdp +
    agr_share_hist |
    country,
  data = region_dat
)

m_sdR_country <- feols(
  log_sd_R ~
    mean_R +
    I(mean_R^2) +
    log_gdp +
    agr_share_hist |
    country,
  data = region_dat
)


m_sdR_country <- feols(
  agr_share_hist ~
    mean_T +
    I(mean_T^2) +
    mean_R +
    I(mean_R^2) +
    log_gdp +
    agr_share_hist |
    country,
  data = region_dat
)
