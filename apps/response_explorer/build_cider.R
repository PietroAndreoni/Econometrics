# Regional patterns of stratospheric aerosol injection (SAI) by injection
# latitude for the termination page of the response explorer, from CIDER, the
# Climate Intervention Dynamical EmulatoR (Farley, MacMartin, Visioni, Kravitz,
# Bednarz, Duffey, Henry and Akherati, Geosci. Model Dev. 19, 1809-1831, 2026,
# doi:10.5194/gmd-19-1809-2026; code and parameters v1.1.0, Zenodo
# 10.5281/zenodo.17553537, CC BY 4.0). Saved to
# results/response_explorer/cider.rds.
#
# CIDER is trained separately on CESM2-WACCM6 and UKESM1. For each model it
# gives, for SO2 injected at 60N, 30N, 15N, 0N, 15S, 30S or 60S (and for the
# non-SAI forcing, "GHG"):
#   - patterns of annual T and P change (288 x 192 grid, the last 20 years
#     of the 12 Tg SO2/yr step runs), each scaled to a global mean of one;
#   - the impulse responses from injection to global-mean AOD and from AOD to
#     global-mean T and P (monthly).
# Here they give, by model and injection latitude:
#   efficiency  cooling_per_tg: the global cooling per Tg SO2/yr (C), and
#               dp_per_dt: the global P change per C of global T change
#               (mm/d per C), both at the steady state of the training runs:
#               the emulator's response to 12 Tg/yr over 35 years, mean of its
#               last 20 (the years the patterns come from), divided by 12 and
#               by the T response;
#   pattern     by country (iso3) and weight (population, area): t, the
#               country's T change per C of global T change, and p_rel, its
#               relative P change per mm/d of global P change (the P change
#               of the pattern over the model's own P climatology, each
#               averaged over the country: the change is relative to the
#               model's climate, not its bias). The same for "GHG".
# The country means weight each model grid cell by the population (Kummu et
# al. 2025, 2000 as the app's climate data, 5 arcmin) or the area of the country in it
# (GADM 4.1 borders, data/kummu2025/gadm41): every 5-arcmin cell's population
# or covered area goes to the model cell whose centre is nearest.
#
# The app uses no injection -> AOD -> T dynamics: the SAI cooling of each year
# is the gap between the geoengineered and the background warming, split
# across the latitudes by injection fraction times cooling_per_tg.
#
# Inputs: data/CIDER/CESM_params.mat and UKESM_params.mat (CIDER_code/
# parameters of the Zenodo archive).
# From the project root:  Rscript apps/response_explorer/build_cider.R

suppressPackageStartupMessages(library(dplyr))
source("load_functions.R")

CIDER_DIR <- file.path("data", "CIDER")
CIDER_MODELS <- c(CESM = "CESM2-WACCM6", UKESM = "UKESM1")
# The order of the injection latitudes in the parameter files (checked from
# the hemispheric AOD of the patterns), then the non-SAI forcing.
CIDER_LATS <- c("60N", "30N", "15N", "0N", "15S", "30S", "60S")
CIDER_SOURCES <- c(CIDER_LATS, "GHG")
# The training step: 12 Tg SO2/yr (1 Tg a month) for 35 years; the patterns
# are the mean of the last 20.
STEP_TG_MONTH <- 1
STEP_MONTHS <- 35L * 12L
PATTERN_MONTHS <- 20L * 12L
# The year of the population weights: 2000, as the climate data of the app's
# population-weighted panels (population density 2000, build_models.R).
POPULATION_LAYER <- "population_2000"

# The emulator (CIDER_code/CIDER.py) --------------------------------------------
# Global-mean AOD from a monthly injection series q (Tg/month) at one latitude:
#   AOD(k) = sum_{j=0}^{k-1} beta exp(-(alpha + gamma q(k-j)) (j + 1)) q(k-j).
cider_aod <- function(p, q) {
  n <- length(q)
  vapply(seq_len(n), function(k) {
    j <- 0:(k - 1L)
    qq <- q[k - j]
    sum(p[[1]] * exp(-(p[[2]] + p[[3]] * qq) * (j + 1)) * qq)
  }, 0)
}
# The semi-infinite diffusion impulse response (months 1..n) of parameters
# p = (tau, mu), with the asymptotic form where exp(s) erfc(sqrt(s)) overflows.
cider_impulse <- function(p, n) {
  s <- seq_len(n) / p[[1]]
  e <- exp(s) * 2 * stats::pnorm(-sqrt(2 * s))
  bad <- !is.finite(e)
  e[bad] <- 2 / sqrt(pi) / (sqrt(s[bad]) + sqrt(s[bad] + 2))
  p[[2]] * (1 / sqrt(pi * s) - e)
}
# The response to forcing f: y(k) = sum_{j=1}^{k} h(j) f(k-j+1).
cider_response <- function(p, f) {
  h <- cider_impulse(p, length(f))
  vapply(seq_along(f), function(k) sum(h[seq_len(k)] * f[k:1]), 0)
}

read_cider <- function(model) {
  m <- R.matlab::readMat(file.path(CIDER_DIR, paste0(model, "_params.mat")))
  stopifnot(identical(dim(m$param.AOD.all), c(7L, 3L)), identical(dim(m$param.T.all), c(8L, 2L)),
            identical(dim(m$pattern.T.all), c(288L, 192L, 8L)), identical(dim(m$pattern.P.base), c(288L, 192L)))
  m
}

efficiency_of <- function(model, m) {
  q <- rep(STEP_TG_MONTH, STEP_MONTHS)
  last <- seq(STEP_MONTHS - PATTERN_MONTHS + 1L, STEP_MONTHS)
  bind_rows(lapply(seq_along(CIDER_LATS), function(i) {
    aod <- cider_aod(m$param.AOD.all[i, ], q)
    dT <- mean(cider_response(m$param.T.all[i, ], aod)[last])
    dP <- mean(cider_response(m$param.P.all[i, ], aod)[last])
    tibble::tibble(model = model, lat = CIDER_LATS[[i]], aod = mean(aod[last]),
                   cooling_per_tg = -dT / (12 * STEP_TG_MONTH), dp_per_dt = dP / dT)
  }))
}

# Country weights over the model grid ------------------------------------------
# The model grid: 288 longitudes from 0 by 1.25 degrees, 192 latitudes from -90
# to 90 (CESM2 f09; the UKESM1 patterns are on the same grid). Each 5-arcmin
# cell of every GADM 4.1 region goes to the nearest model cell, with its
# population (times the share of it inside the region) and its covered area.
NLON <- 288L
NLAT <- 192L
model_cell <- function(x, y) {
  i <- (round((x %% 360) / 1.25) %% NLON) + 1L
  j <- round((y + 90) / (180 / (NLAT - 1L))) + 1L
  (j - 1L) * NLON + i
}
country_weights <- function() {
  old_s2 <- sf::sf_use_s2(FALSE)
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  regions <- read_gadm1_polygons(file.path(data_root(), "kummu2025", "gadm41"))
  population <- terra::rast(file.path(data_root(), "kummu2025", "kummu_population_1990_2024_5arcmin.tif"))
  population <- population[[POPULATION_LAYER]]
  message("Weights of ", nrow(regions), " GADM 4.1 regions over the model grid.")
  parts <- exactextractr::exact_extract(
    population, regions, include_xy = TRUE, include_area = TRUE, include_cols = "GID_0",
    summarize_df = TRUE, progress = FALSE,
    fun = function(df) {
      v <- df$value
      v[!is.finite(v)] <- 0
      cell <- model_cell(df$x, df$y)
      out <- stats::aggregate(cbind(population = v * df$coverage_fraction, area = df$area * df$coverage_fraction),
                              list(cell = cell), sum)
      out$iso3 <- df$GID_0[[1]]
      out
    })
  bind_rows(parts) %>%
    group_by(iso3, cell) %>%
    summarise(population = sum(population), area = sum(area), .groups = "drop")
}

# The country means of the patterns of model m under weights w (population or
# area; a country without population takes its area).
country_patterns <- function(model, m, weights) {
  flat <- function(a) as.vector(a)  # [lon, lat] column-major: cell = (lat - 1) * NLON + lon
  base_P <- flat(m$pattern.P.base)
  bind_rows(lapply(c(population = "population", area = "area"), function(wname) {
    w <- weights
    w$w <- w[[wname]]
    if (wname == "population") {
      no_pop <- w %>% group_by(iso3) %>% summarise(total = sum(population)) %>% filter(!(total > 0))
      w$w[w$iso3 %in% no_pop$iso3] <- w$area[w$iso3 %in% no_pop$iso3]
    }
    bind_rows(lapply(seq_along(CIDER_SOURCES), function(k) {
      t_cell <- flat(m$pattern.T.all[, , k])
      p_cell <- flat(m$pattern.P.all[, , k])
      w %>%
        group_by(iso3) %>%
        summarise(t = sum(w * t_cell[cell]) / sum(w),
                  p_rel = sum(w * p_cell[cell]) / sum(w * base_P[cell]), .groups = "drop") %>%
        mutate(model = model, source = CIDER_SOURCES[[k]], weight = wname)
    }))
  }))
}

weights <- country_weights()
models <- lapply(stats::setNames(names(CIDER_MODELS), names(CIDER_MODELS)), read_cider)
efficiency <- bind_rows(lapply(names(models), function(k) efficiency_of(k, models[[k]])))
patterns <- bind_rows(lapply(names(models), function(k) country_patterns(k, models[[k]], weights)))
pattern <- lapply(c(population = "population", area = "area"), function(w) {
  patterns %>% filter(weight == w) %>% select(model, source, iso3, t, p_rel)
})

print(as.data.frame(efficiency), digits = 3)
# Check: CIDER's GHG pattern against the CMIP6 slopes of the app.
scen_file <- file.path("results", "response_explorer", "scenarios.rds")
if (file.exists(scen_file)) {
  cmip <- readRDS(scen_file)$pattern$population
  for (k in names(CIDER_MODELS)) {
    g <- pattern$population %>% filter(model == k, source == "GHG") %>% inner_join(cmip, by = "iso3")
    message(sprintf("%s GHG pattern vs CMIP6 median slope, %d countries: T correlation %.2f, mean %.2f vs %.2f",
                    k, nrow(g), stats::cor(g$t, g$t_slope), mean(g$t), mean(g$t_slope)))
  }
}

out_dir <- file.path("results", "response_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(models = CIDER_MODELS, lats = CIDER_LATS, efficiency = efficiency, pattern = pattern,
             source = paste("CIDER v1.1.0 (Farley et al., GMD 2026, doi:10.5194/gmd-19-1809-2026;",
                            "Zenodo 10.5281/zenodo.17553537, CC BY 4.0)"),
             built = Sys.time()),
        file.path(out_dir, "cider.rds"))
message("Saved ", file.path(out_dir, "cider.rds"), ": ", length(CIDER_MODELS), " models, ",
        n_distinct(patterns$iso3), " countries.")
