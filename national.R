library(sf)
library(ggplot2)

window_rr <- 10
window_tm <- 10
################################################################################
econ_data_iso <- arrow::read_parquet("econometrics/data/data_gdp_gilli.parquet") %>%
  select(iso3,year,gdppc)%>% filter(!is.na(gdppc)) %>% rename(grp_pc_usd=gdppc)
climate_data_gilli <- arrow::read_parquet("econometrics/data/data_gdp_gilli.parquet") %>%
  select(iso3,year,temperature_mean,precip) %>% rename(RR=precip,TM=temperature_mean)
climate_data_country_pop <- read.csv("econometrics/data/data_climate_country_pop_weight_1950_2023_updated.csv")
climate_data_country_area <- read.csv("econometrics/data/data_climate_country_1950_2023_updated.csv")

data_national <- econ_data_iso %>%
  full_join(climate_data_country_area %>% select(year,iso3,TM)) %>% 
  full_join(climate_data_country_area %>% select(year,iso3,RR)) %>%
#  full_join(climate_data_gilli) %>%
  rename(GID_1=iso3) %>%
  group_by(GID_1) %>%
  arrange(year,.by_group=TRUE) %>%
  group_by(GID_1) %>% 
  mutate(dlgrp_pc_usd = log(grp_pc_usd) - log(lag(grp_pc_usd)),
         lgrp_pc_usd=log(grp_pc_usd),
         lag_lgrp_pc_usd=log(lag(grp_pc_usd)),
         RR=RR/1000) %>% ungroup() 

# simple means and sd over the full dataset period
means <- data_national %>% 
  filter(year>=1979 & year<=2019) %>%
  group_by(GID_1) %>%
  summarise(sd_TM_all=sd(TM,na.rm=TRUE),sd_RR_all=sd(RR,na.rm=TRUE),mean_TM_all=mean(TM,na.rm=TRUE),mean_RR_all=mean(RR,na.rm=TRUE))

data_national <- data_national %>% 
  inner_join(means) %>%
  group_by(GID_1) %>%
  mutate(dev_TM_all = (TM - mean_TM_all)/sd_TM_all,
         dev_RR_all = (RR - mean_RR_all)/sd_RR_all,
         dTM = TM - lag(TM),
         dRR = RR - lag(RR)) %>%
  ungroup() %>%
  mutate(dRR_2 = dRR^2, 
         dTM_2=dTM^2,
         dev_TM_all_2 = dev_TM_all^2,
         dev_RR_all_2 = dev_RR_all^2,
         RR_2=RR^2,
         TM_2=TM^2)

# start with estimates
i <- "GID_1 + year + GID_1[year]+GID_1[year^2]"
pan_id<-c("GID_1", "year")

data_national <- data_national %>% filter(!is.na(dlgrp_pc_usd))

# first, burke: precipitations linear economic good and temperature u shape with maximum at 9 °C
f <- as.formula(paste( "dlgrp_pc_usd", "~", "TM + TM_2 +
                                             RR + RR_2", "|" ,i ))

m_bhm <- fixest::feols(f, data_national, panel.id=pan_id)
#m_bhm2 <- fixest::feols(f, data_national[which(abs(m_bhm$residuals)<1),], panel.id=pan_id)
summary(m_bhm)

t_bhm  <- m_bhm$coeftable["TM","Estimate"]
t2_bhm <-  m_bhm$coeftable["TM_2","Estimate"]
p_bhm <-  m_bhm$coeftable["RR","Estimate"]
p2_bhm <-  m_bhm$coeftable["RR_2","Estimate"]

t_bhm/-(2*t2_bhm)
p_bhm/-(2*p2_bhm)

# kalhul and wenz
f <- as.formula(paste( "dlgrp_pc_usd", "~", "dTM + dTM:TM + TM + TM_2 +
                                             dRR + dRR:RR + RR + RR_2", "|" ,i ))

m_kw <- fixest::feols(f, data_national, panel.id=pan_id)
summary(m_kw)

# specialization (dev_all_TM2 and dev_all_RR2) 
f <- as.formula(paste( "dlgrp_pc_usd", "~", "TM + TM_2 + RR + RR_2 + dev_TM_all_2 + dev_RR_all_2", "|" ,i ))

m_spec <- fixest::feols(f, data_national, panel.id=pan_id)
summary(m_spec)
