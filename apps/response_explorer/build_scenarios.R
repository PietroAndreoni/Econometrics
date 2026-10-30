# Inputs of the climate-scenario projections of the response explorer, saved
# to results/response_explorer/scenarios.rds:
#
#   gmt      global surface air temperature above 1850-1900 (C) of SSP1-1.9,
#            SSP1-2.6, SSP4-3.4, SSP2-4.5, SSP4-6.0, SSP3-7.0 and SSP5-8.5,
#            1750-2100: median and 5th
#            and 95th percentiles over the 1,001 members of FaIR v2.1 calibrated
#            to IPCC AR6 (Smith, Zenodo 10.5281/zenodo.7838148, CC BY 4.0; the
#            netCDF files are in data/fair_ssp_zenodo7838148)
#   pattern  each country's response to global warming, the median over the
#            CMIP6 models of data/cmip6_downscaling_temp.parq and _prec.parq
#            (best-fit linear relationship of the country's T or P with global
#            temperature; the columns without suffix for population weights,
#            with _area for area weights): the T slope (C per
#            C of global warming) and the relative P slope (each model's slope
#            over its own intercept, so the change does not depend on the
#            model's bias in the level of P)
#   ssp      the socioeconomic pathway of each forcing level (SSP1 for 1.9 and
#            2.6, SSP2 for 4.5, SSP4 for 3.4 and 6.0, SSP3 for 7.0, SSP5 for 8.5)
#            and, by country and SSP,
#            annual series of
#              gdppc  GDP per capita, PPP 2017 US$ (OECD ENV-Growth 2025, SSP
#                     database 2025 release), 1980-2100
#              gdp    GDP (gdppc x population, IIASA-WiC POP 2025), 2025-2100:
#                     the weights of the global losses
#              share  agricultural value added, % of GDP ("Structural Change
#                     Model (2023)" of the SSP extensions: WDI to 2024, SSP to
#                     2050), 1998-2100
#            The 5-year points are joined by log-linear interpolation (constant
#            growth between them). After 2050 the share follows income with each
#            country's elasticity of the share to SSP GDP per capita over
#            2030-2050. Downloaded from IIASA's public ixmp4 platforms ("ssp",
#            "ssp-extensions") into data/ssp_iiasa_ixmp4.
#   history  the model's own historical series, which the SSP changes are
#            chained onto: GDP per capita of each GDP source (PWT: rgdpna / pop;
#            WB: constant 2015 US$) and the agricultural share of the three
#            sectors' value added (World Bank, constant 2015 US$, lagged as in
#            the models)
#   world    country borders (Natural Earth, 1:110m) for the maps
#   palettes IPCC AR6 colour maps of T change (temp_seq) and P change (prec_div)
#
# From the project root:  Rscript apps/response_explorer/build_scenarios.R

suppressPackageStartupMessages(library(dplyr))
source("load_functions.R")

SCENARIOS <- c(ssp119 = "SSP1-1.9", ssp126 = "SSP1-2.6", ssp434 = "SSP4-3.4",
               ssp245 = "SSP2-4.5", ssp460 = "SSP4-6.0", ssp370 = "SSP3-7.0",
               ssp585 = "SSP5-8.5")
SCENARIO_SSP <- c(ssp119 = "SSP1", ssp126 = "SSP1", ssp434 = "SSP4", ssp245 = "SSP2",
                  ssp460 = "SSP4", ssp370 = "SSP3", ssp585 = "SSP5")
FAIR_DIR <- file.path("data", "fair_ssp_zenodo7838148")
SSP_DIR <- file.path("data", "ssp_iiasa_ixmp4")

read_fair_gmt <- function(s) {
  nc <- ncdf4::nc_open(file.path(FAIR_DIR, paste0(s, ".nc")))
  on.exit(ncdf4::nc_close(nc))
  temperature <- ncdf4::ncvar_get(nc, "temperature")
  year <- as.integer(ncdf4::ncvar_get(nc, "year"))
  if (nrow(temperature) != length(year)) temperature <- t(temperature)
  keep <- year <= 2100
  tibble::tibble(scenario = s, label = SCENARIOS[[s]], year = year[keep],
                 median = apply(temperature[keep, ], 1, stats::median),
                 p05 = apply(temperature[keep, ], 1, stats::quantile, 0.05, names = FALSE),
                 p95 = apply(temperature[keep, ], 1, stats::quantile, 0.95, names = FALSE))
}
gmt <- bind_rows(lapply(names(SCENARIOS), read_fair_gmt))

temp <- as.data.frame(arrow::read_parquet(data_file("cmip6_downscaling_temp.parq")))
prec <- as.data.frame(arrow::read_parquet(data_file("cmip6_downscaling_prec.parq")))
# The median of the models' slopes (T) and relative slopes (P). The R2-weighted
# means are kept for reference only: R2 grows with the size of a model's
# response, so weighting by it selects the models with the largest responses
# (for P it nearly doubles the spread across countries).
# One pattern for each climate weighting of the models: the columns without
# suffix (population) or with _area (plain area average).
pattern_of <- function(suffix) {
  t_coef <- temp[[paste0("gmt_coef", suffix)]]
  t_r2 <- temp[[paste0("r2", suffix)]]
  p_coef <- prec[[paste0("gmt_coef", suffix)]]
  p_intercept <- prec[[paste0("Intercept_coef", suffix)]]
  p_r2 <- prec[[paste0("r2_pr", suffix)]]
  full_join(
    data.frame(iso3 = temp$iso3, coef = t_coef, r2 = t_r2) %>% filter(is.finite(coef)) %>%
      group_by(iso3) %>%
      summarise(t_slope = stats::median(coef), t_slope_r2_weighted = stats::weighted.mean(coef, r2),
                t_r2 = mean(r2), t_models = n(), .groups = "drop"),
    data.frame(iso3 = prec$iso3, ratio = p_coef / p_intercept, intercept = p_intercept, r2 = p_r2) %>%
      filter(intercept > 0, is.finite(ratio)) %>%
      group_by(iso3) %>%
      summarise(p_rel_slope = stats::median(ratio),
                p_rel_slope_r2_weighted = stats::weighted.mean(ratio, r2),
                p_r2 = mean(r2), p_models = n(), .groups = "drop"),
    by = "iso3"
  )
}
pattern <- list(population = pattern_of(""), area = pattern_of("_area"))

# SSP data --------------------------------------------------------------------
# IAMC datapoints of a public IIASA ixmp4 platform (all pages), read from the
# local copy in SSP_DIR when there is one.
ixmp4_datapoints <- function(platform, filter, file, page = 10000L) {
  path <- file.path(SSP_DIR, file)
  if (file.exists(path)) return(utils::read.csv(path, stringsAsFactors = FALSE))
  url <- sprintf(paste0("https://ixmp4.ece.iiasa.ac.at/v1/%s/iamc/datapoints/",
                        "?table=true&join_parameters=true&join_runs=true"), platform)
  out <- list()
  offset <- 0
  repeat {
    h <- curl::new_handle(customrequest = "PATCH",
                          postfields = jsonlite::toJSON(filter, auto_unbox = TRUE))
    curl::handle_setheaders(h, "Content-Type" = "application/json")
    res <- curl::curl_fetch_memory(paste0(url, "&limit=", page, "&offset=", offset), handle = h)
    if (res$status_code != 200) stop("ixmp4 ", platform, ": HTTP ", res$status_code)
    j <- jsonlite::fromJSON(rawToChar(res$content))
    d <- as.data.frame(j$results$data, stringsAsFactors = FALSE)
    names(d) <- j$results$columns
    out[[length(out) + 1]] <- d
    offset <- offset + nrow(d)
    if (nrow(d) < page || offset >= j$total) break
  }
  d <- do.call(rbind, out)[, c("model", "scenario", "region", "variable", "unit", "step_year", "value")]
  dir.create(SSP_DIR, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(d, path, row.names = FALSE)
  d
}
# Country names of the SSP database to ISO3 (aggregates have none).
with_iso3 <- function(d) {
  d$iso3 <- suppressWarnings(countrycode::countrycode(d$region, "country.name", "iso3c"))
  d$year <- as.integer(d$step_year)
  d$value <- as.numeric(d$value)
  d[!is.na(d$iso3) & is.finite(d$value), ]
}
# Annual values between the first and the last point of each series by
# log-linear interpolation (constant growth between the points); no
# extrapolation.
log_linear_annual <- function(d, keys) {
  d %>%
    group_by(across(all_of(keys))) %>%
    group_modify(function(s, k) {
      s <- s[s$value > 0, ]
      if (nrow(s) < 2) return(tibble::tibble(year = s$year, value = s$value))
      years <- seq(min(s$year), max(s$year))
      tibble::tibble(year = years, value = exp(stats::approx(s$year, log(s$value), xout = years)$y))
    }) %>%
    ungroup()
}

gdppc_raw <- with_iso3(ixmp4_datapoints(
  "ssp", list(variable = list(name = "GDP|PPP [per capita]"), model = list(name = "OECD ENV-Growth 2025")),
  "gdppc_ppp2017_oecd_envgrowth_2025.csv")) %>%
  filter(unit == "USD_2017/yr")
population_raw <- with_iso3(ixmp4_datapoints(
  "ssp", list(variable = list(name = "Population"), model = list(name = "IIASA-WiC POP 2025")),
  "population_iiasa_wic_2025.csv"))
share_raw <- with_iso3(ixmp4_datapoints(
  "ssp-extensions", list(variable = list(name = "Value Added|Agriculture [Share]")),
  "agriculture_share_structural_change_2023.csv"))

ssps <- sort(unique(SCENARIO_SSP))
# GDP per capita: the historical reference before 2025 shared by every SSP,
# then the SSP's own points.
gdppc <- bind_rows(lapply(ssps, function(s) {
  bind_rows(gdppc_raw %>% filter(scenario == "Historical Reference", year < 2025),
            gdppc_raw %>% filter(scenario == s)) %>%
    mutate(ssp = s)
})) %>%
  distinct(ssp, iso3, year, .keep_all = TRUE) %>%
  log_linear_annual(c("ssp", "iso3")) %>%
  rename(gdppc = value)
population <- population_raw %>%
  filter(scenario %in% ssps) %>%
  transmute(ssp = scenario, iso3, year, value) %>%
  distinct(ssp, iso3, year, .keep_all = TRUE) %>%
  log_linear_annual(c("ssp", "iso3")) %>%
  rename(population = value)
gdp <- inner_join(gdppc, population, by = c("ssp", "iso3", "year")) %>%
  filter(year >= 2025) %>%
  transmute(ssp, iso3, year, gdp = gdppc * population)
# Agricultural share: the SSP series (history included) to 2050, then each
# country's 2030-2050 elasticity to SSP GDP per capita.
share_to_2050 <- share_raw %>%
  filter(scenario %in% ssps) %>%
  transmute(ssp = scenario, iso3, year, value) %>%
  distinct(ssp, iso3, year, .keep_all = TRUE) %>%
  log_linear_annual(c("ssp", "iso3")) %>%
  rename(share = value)
share_after_2050 <- share_to_2050 %>%
  filter(year %in% c(2030, 2050)) %>%
  inner_join(gdppc, by = c("ssp", "iso3", "year")) %>%
  group_by(ssp, iso3) %>%
  filter(n() == 2) %>%
  summarise(share_2050 = share[year == 2050], gdppc_2050 = gdppc[year == 2050],
            elasticity = (log(share[year == 2050]) - log(share[year == 2030])) /
              (log(gdppc[year == 2050]) - log(gdppc[year == 2030])),
            .groups = "drop") %>%
  filter(is.finite(elasticity)) %>%
  inner_join(gdppc %>% filter(year > 2050), by = c("ssp", "iso3")) %>%
  transmute(ssp, iso3, year, share = share_2050 * (gdppc / gdppc_2050)^elasticity)
share <- bind_rows(share_to_2050, share_after_2050) %>% arrange(ssp, iso3, year)

# The model's own historical series --------------------------------------------
income_history <- lapply(c(PWT = "PWT110", WB = "WB"), function(src) {
  load_econ_panel(src) %>%
    filter(gadm_level == "gadm0", !is.na(grp_pc_usd), grp_pc_usd > 0) %>%
    transmute(iso3 = GID_0, year = as.integer(year), gdppc = grp_pc_usd) %>%
    distinct(iso3, year, .keep_all = TRUE)
})
share_history <- load_wdi_growth() %>%
  filter(is.finite(s_agr)) %>%
  transmute(iso3 = GID_0, year = as.integer(year), share = s_agr)

# IPCC AR6 colour maps (IPCC-WG1/colormaps, continuous, RGB 0-1) for the maps of
# the changes of T (sequential) and P (diverging), in data/ipcc_colormaps.
ipcc_palette <- function(name) {
  rgb <- utils::read.table(file.path("data", "ipcc_colormaps", paste0(name, ".txt")))
  grDevices::rgb(rgb[[1]], rgb[[2]], rgb[[3]])
}
palettes <- lapply(c(temp_seq = "temp_seq", prec_div = "prec_div"), ipcc_palette)

world <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
  transmute(iso3 = adm0_a3, name = name) %>%
  sf::st_transform("+proj=robin")

out_dir <- file.path("results", "response_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(scenarios = SCENARIOS, gmt = gmt, pattern = pattern,
             ssp = list(of_scenario = SCENARIO_SSP, gdppc = gdppc, gdp = gdp, share = share),
             history = list(income = income_history, share = share_history), world = world,
             palettes = palettes,
             source = paste("FaIR v2.1 (AR6-calibrated, 1,001 members), Smith (2023),",
                            "Zenodo 10.5281/zenodo.7838148"),
             ssp_source = paste("SSP database 2025 release (OECD ENV-Growth 2025, IIASA-WiC POP",
                                "2025) and SSP extensions (Structural Change Model 2023), IIASA"),
             built = Sys.time()),
        file.path(out_dir, "scenarios.rds"))
message("Saved ", file.path(out_dir, "scenarios.rds"), ": ", length(SCENARIOS), " scenarios, ",
        nrow(pattern$population), " countries with a pattern; SSP GDP for ", n_distinct(gdp$iso3),
        " countries, agricultural share for ", n_distinct(share$iso3))
