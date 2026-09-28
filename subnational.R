library(tidyverse)
# specify formulas and models
# fixed effects
i <- "GID_1 + year + GID_1[year] + GID_1[year^2]"#"GID_1 + GID_0 + year + GID_1[year] + GID_1[year^2] + GID_0[year] + GID_0[year^2]"
# panel indexes
pan_id <- c("GID_1", "year")
# output variable 
o <- "dlgrp_pc_usd"
# model
model <- "TM + TM_2 + dev_TM_all_2 + RR + RR_2 + dev_RR_all_2"

sel_source <- "KUMMU2018"
weight_p <- "area"
weight_t <- "area"
rr <- 10
tm <- 10

data_m <- data_rr %>% filter(!is.na(dlgrp_pc_usd) & window_rr==rr & weight==weight_p & source==sel_source)  %>% select(-window_rr,-weight,-source) %>%
  full_join(data_tm %>% filter(!is.na(dlgrp_pc_usd) & window_tm==tm & weight==weight_t & source==sel_source) %>% select(-window_tm,-weight,-source))

f <- as.formula(paste( o, "~", model, "|" ,i ))
m <- fixest::feols(f, data_m, panel.id=pan_id)
m1 <- summary(m, vcov=~GID_0)


data_m <- data_m %>% 
  group_by(GID_1) %>% 
  mutate(temp_bin = round(mean_TM_all/10/100,2)*100*10) %>%
  mutate(prec_bin = round(mean_RR_all/2/100,2)*100*2) %>%
  rowwise() %>% mutate(temp_bin = max(0,min(30,temp_bin)),prec_bin=max(0,min(4,prec_bin)))
data_m <- data_m %>% 
  group_by(GID_0) %>%
  mutate(poor=ifelse(mean(grp_pc_usd,na.rm=TRUE) < 14000, "yes","no"))
data_m <- data_m %>% mutate(temp_bin=as.numeric(temp_bin),RR=RR)

climate_labels <- c(
  `1` = "Tropical",
  `2` = "Dry",
  `3` = "Temperate",
  `4` = "Continental",
  `5` = "Polar")

f <- as.formula(paste( "dlgrp_pc_usd", "~", "dev_TM + dev_TM_2 + dev_RR + dev_RR_2", "|" ,i ))
for (zone in names(climate_labels) ) {
  print(climate_labels[zone])
  data <- data_m %>% inner_join(kg_data) %>% filter(mean_kg == as.numeric(zone))
  m <- fixest::feols(f, data, panel.id=pan_id)
  print(summary(m))
  summary(m, vcov=~GID_0)
  m$coefficients
  vcov(m) }

for (bin in unique(data_subnational$temp_bin) ) {
  print(bin)
  data <- data_m %>% filter(temp_bin == bin)
  m <- fixest::feols(f, data, panel.id=pan_id)
  print(summary(m))
  summary(m, vcov=~GID_0)
  m$coefficients
  vcov(m) }

for (poorsel in c("yes","no") ) {
  print(paste0("Poor:",poorsel))
  data <- data_m %>% filter(poor == poorsel)
  m <- fixest::feols(f, data, panel.id=pan_id)
  print(summary(m))
  summary(m, vcov=~GID_0)
  m$coefficients
  vcov(m) }

data <- data_m %>% filter(year >= 1990 & year<=2012)
m <- fixest::feols(f, data, panel.id=pan_id)
print(summary(m))
summary(m, vcov=~GID_0)
m$coefficients
vcov(m)