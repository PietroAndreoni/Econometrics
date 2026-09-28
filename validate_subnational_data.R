econ_data_iso <- read.csv("econometrics/data/wb_gdp_data.csv") %>%
  pivot_longer(!c(Series.Name,Series.Code,Country.Name,Country.Code),names_to="year",values_to="gdppc") %>%
  filter(Series.Name=="GDP per capita (constant 2015 US$)") %>% 
  rename(iso3=Country.Code) %>% 
  mutate(year=as.numeric(str_extract(year,"(?<=..YR)(.*)(?=.)")),gdppc=ifelse(gdppc=="..",NA,gdppc)) %>% 
  mutate(gdppc=as.numeric(gdppc)) %>% 
  select(iso3,year,gdppc)%>% 
  filter(!is.na(gdppc)) %>% 
  rename(gdpc_wb=gdppc,GID_0=iso3) %>% 
  group_by(GID_0) %>%
  arrange(year,.by_group=TRUE) %>% 
  group_by(GID_0) %>% 
  mutate(dglog_wb = log(gdpc_wb) - log(lag(gdpc_wb)) )

econ_data_dose <- read.csv("econometrics/data/DOSE_V2.11.csv") %>% 
  select(year,GID_0,GID_1,grp_pc_usd_2015,pop) %>% 
  rename(grp_pc_usd = grp_pc_usd_2015) %>% 
  filter(!GID_1 %in% c(" ","")) %>% 
  mutate(GID_0 = stringr::str_extract(GID_1, "^.{3}")) %>% 
  group_by(GID_0,year) %>%  
  summarise(gdpc_dose=weighted.mean(grp_pc_usd,pop,na.rm=TRUE)) 

econ_data_dose <- econ_data_dose %>% 
  full_join(cross_join(data.frame(GID_0=unique(econ_data_dose$GID_0)),data.frame(year=seq(1960,2020,by=1) ) ) ) %>% 
  group_by(GID_0) %>%
  arrange(year,.by_group=TRUE) %>% 
  ungroup() %>% 
  group_by(GID_0) %>% 
  mutate(dglog_dose = log(gdpc_dose) - log(lag(gdpc_dose)) )

econ_data_subnat <- read.csv("econometrics/data/DOSE_V2.11.csv") %>% 
  select(year,GID_0,GID_1,grp_pc_usd_2015,pop) %>% 
  rename(gdpc_dose = grp_pc_usd_2015) %>% 
  filter(!GID_1 %in% c(" ","")) %>% 
  mutate(GID_0 = stringr::str_extract(GID_1, "^.{3}"))

econ_data_subnat <- econ_data_subnat %>% 
  full_join(cross_join(data.frame(GID_1=unique(econ_data_subnat$GID_1)),data.frame(year=seq(1960,2020,by=1) ) ) ) %>% 
  group_by(GID_1) %>%
  arrange(year,.by_group=TRUE) %>%  
  group_by(GID_1) %>% 
  mutate(dglog_dose = log(gdpc_dose) - log(lag(gdpc_dose)) )

econ_data_kummu <-read.csv("econometrics/data/kummu_aggregated_gid1_gdp.csv") %>% 
  mutate(GID_0=stringr::str_extract(GID_1, "^.{3}")) %>% 
  select(year,GID_0,GID_1,gdp_pc) %>% 
  rename(gdpc_kummu = gdp_pc) 

econ_data_kummu <- econ_data_kummu %>% 
  full_join(cross_join(data.frame(GID_1=unique(econ_data_kummu$GID_1)),data.frame(year=seq(1960,2020,by=1) ) ) ) %>% 
  group_by(GID_1) %>%
  arrange(year,.by_group=TRUE) %>%  
  group_by(GID_1) %>% 
  mutate(dglog_kummu = log(gdpc_kummu) - log(lag(gdpc_kummu)) )

remove_outliers <- inner_join(econ_data_dose,econ_data_iso) %>% 
  filter(!is.na(dglog_dose) & !is.na(dglog_wb)) %>% 
#  mutate(err=abs(dglog_dose-dglog_wb)) %>% 
  group_by(GID_0) %>%
  summarise(err = sum((dglog_dose-dglog_wb)^2)/n()) %>% 
  filter(abs(err)<=quantile(err,0.95)) %>% 
  select(GID_0,err)

ncheck <- "IND"
ggplot(full_join(econ_data_dose,econ_data_iso) %>% filter(GID_0==ncheck)) +
  geom_line(data=econ_data_subnat %>% filter(GID_0==ncheck),aes(x=year,y=dglog_dose,group=GID_1),color="grey") +
  geom_line(data=econ_data_kummu %>% filter(GID_0==ncheck),aes(x=year,y=dglog_kummu,group=GID_1),color="blue") +
geom_line(aes(x=year,y=dglog_dose)) +
  geom_line(aes(x=year,y=dglog_wb), color="red") +
  geom_point(aes(x=year,y=dglog_dose)) +
  geom_point(aes(x=year,y=dglog_wb), color="red") + ylim(c(-0.5,0.5))
