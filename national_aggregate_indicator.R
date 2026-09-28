
aggregate <- read_parquet("data/national_gadm2_era5_population_weighted_30y_shocks.parquet") |> 
  filter(variable=="temperature") |> 
  select(GID_0,year,tmz2 = population_weighted_z2, tmz2m =population_weighted_z2_negative, tmz2p =population_weighted_z2_positive) |> ungroup() |> 
  full_join(read_parquet("data/national_gadm2_era5_population_weighted_30y_shocks.parquet") |> 
              filter(variable=="precipitation") |> 
              select(GID_0,year,rrz2 =population_weighted_z2, rrz2m =population_weighted_z2_negative, rrz2p =population_weighted_z2_positive) |> ungroup() )  |> unique()

data_iso3 <- build_dat(econ_data = "PWT")

data_iso3 <- data_iso3 |> select(year,GID_0,lgrp_pc_usd,dlgrp_pc_usd,TM,RR) |> 
  inner_join(aggregate) |> 
  rename()

fixest::demean(data_iso3$dlgrp_pc_usd, data_iso3[,c("GID_0","year")] )
fixest::feols(
  as.formula(paste("dlgrp_pc_usd ~",
    "+ tmz2p + tmz2m + rrz2p + rrz2m + rrz2m:rrz2p + tmz2p:tmz2m",
    "| GID_0 + year"
  )),
  data = data_iso3,
  panel.id = c("GID_0","year"),
  vcov = ~GID_0
)
