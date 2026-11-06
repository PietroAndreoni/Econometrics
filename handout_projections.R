# Projections for the handout (handout_preferred_specification.Rmd), with the
# response explorer's projection code (apps/response_explorer/R/projection.R)
# and scenario data (results/response_explorer/scenarios.rds), on the preferred
# models of handout_preferred_specification.R (run it first). Saved to
# results/handout_preferred_specification/projections.rds.
#
#   1. The linear preferred specification under SSP1-2.6 and SSP2-4.5 (each
#      with its own SSP), the climate of the interactions and the trailing
#      baselines following the warming, with income and agricultural share held
#      at the country's 2025 values ("fixed") or following the SSP ("SSP").
#   2. Termination of perfect geoengineering (warming held at 1.5 C from 2025)
#      under SSP2-4.5 in 2050 and 2070, with the quadratic specification, with
#      income and share fixed or following the SSP.
#   3. Natural variability: the ways of resampling it compared, then 1 and 2
#      again with each country's detrended historical T and P resampled in
#      joint 5-year blocks (the explorer's option; needs the variability of
#      build_models.R).
#
# Deviations use the country's long-run SD (sd_mode "fixed"), as estimated.
#
#   Rscript handout_projections.R

suppressPackageStartupMessages(library(dplyr))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- file.path("apps", "response_explorer")
source(file.path(app_dir, "R", "response.R"))
source(file.path(app_dir, "R", "projection.R"))
source(file.path(app_dir, "R", "fit.R"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
SCEN <- readRDS(file.path("results", "response_explorer", "scenarios.rds"))
out_dir <- file.path("results", "handout_preferred_specification")
models <- readRDS(file.path(out_dir, "app_models.rds"))

# The countries: long-run climate and its SD over the horizon (the SD of the
# deviations), with income and share at their 2025 values under the SSP: the
# trailing means the projection code uses (income over 5 years, share over 20)
# of the country's history continued with the SSP; the long-run values where
# the history or the SSP is missing.
units <- P$units[["PWT | ERA5 | Horizon | population | SD"]]
chars_2025 <- function(units, ssp, model) {
  inc_h <- SCEN$history$income$PWT
  inc_s <- SCEN$ssp$gdppc[SCEN$ssp$gdppc$ssp == ssp, ]
  sh_h <- SCEN$history$share
  sh_s <- SCEN$ssp$share[SCEN$ssp$share$ssp == ssp, ]
  at_2025 <- function(h, s, column, window) {
    if (!nrow(h) || !nrow(s)) return(NA_real_)
    trailing_mean(chained_series(data.frame(year = h$year, value = h[[column]]),
                                 data.frame(year = s$year, value = s[[column]])), PROJECTION_FROM, window)
  }
  for (i in seq_len(nrow(units))) {
    iso <- units$GID_0[i]
    inc <- at_2025(inc_h[inc_h$iso3 == iso, ], inc_s[inc_s$iso3 == iso, ], "gdppc", model$inc_window)
    sh <- at_2025(sh_h[sh_h$iso3 == iso, ], sh_s[sh_s$iso3 == iso, ], "share", model$share_window)
    if (is.finite(inc)) units$income[i] <- inc
    if (is.finite(sh)) units$share[i] <- sh
  }
  units
}
spec_of <- function(model, ssp, follow_econ) {
  list(model = model, units = chars_2025(units, ssp, model), sd_mode = "fixed", clim_mode = "baseline",
       cap = FALSE, follow_econ = follow_econ)
}
ECON <- c(fixed = FALSE, SSP = TRUE)
ECON_LABELS <- c(fixed = "Income and share fixed at 2025", SSP = "Income and share follow the SSP")

# 1. Linear specification, SSP1-2.6 and SSP2-4.5.
runs <- list()
for (scenario in c("ssp126", "ssp245")) for (econ in names(ECON)) {
  ssp <- SCEN$ssp$of_scenario[[scenario]]
  message("Projecting ", scenario, ", ", econ)
  r <- project_setting(spec_of(models$Pref, ssp, ECON[[econ]]), SCEN, scenario, ssp)
  runs[[paste(scenario, econ)]] <- list(scenario = scenario, label = SCEN$scenarios[[scenario]], econ = econ, result = r)
}
global <- do.call(rbind, lapply(runs, function(x) {
  cbind(scenario = x$label, econ = ECON_LABELS[[x$econ]], x$result$global)
}))
countries_2100 <- do.call(rbind, lapply(runs, function(x) {
  cbind(scenario = x$label, econ = ECON_LABELS[[x$econ]], x$result$countries[x$result$countries$year == 2100, ])
}))
coverage <- vapply(runs, function(x) x$result$coverage, 0)
n_countries <- vapply(runs, function(x) x$result$n, 0)
followed <- do.call(rbind, lapply(runs, function(x) x$result$followed))

# 2. Termination, quadratic specification, SSP2-4.5.
HELD <- 1.5
TAU <- 4
g <- SCEN$gmt[SCEN$gmt$scenario == "ssp245", ]
background <- stats::setNames(g$median, g$year)
paths <- list(`No geoengineering (SSP2-4.5)` = background,
              `Geoengineering sustained at 1.5 C` = termination_path(background, HELD),
              `Terminated in 2050` = termination_path(background, HELD, 2050L, TAU),
              `Terminated in 2070` = termination_path(background, HELD, 2070L, TAU))
# The quadratic specification (the figure), and the linear one for comparison.
termination <- do.call(rbind, lapply(c(Quadratic = "PrefQ", Linear = "Pref"), function(k_model) {
  do.call(rbind, lapply(names(ECON), function(econ) {
    spec <- spec_of(models[[k_model]], "SSP2", ECON[[econ]])
    do.call(rbind, lapply(names(paths), function(k) {
      message("Termination: ", k_model, ", ", econ, ", ", k)
      r <- project_setting(spec, SCEN, "ssp245", "SSP2", gmt = paths[[k]])
      cbind(model = if (k_model == "PrefQ") "Quadratic arms" else "Linear arms", econ = ECON_LABELS[[econ]],
            path = k, r$global)
    }))
  }))
}))
gmt_paths <- do.call(rbind, lapply(names(paths), function(k) {
  data.frame(path = k, year = as.integer(names(paths[[k]])), warming = unname(paths[[k]]))
}))

# 3. Natural variability (the explorer's option): each country's detrended
# historical T and P resampled in joint 5-year blocks, N_PATHS noise paths,
# each paired with a draw of the coefficients. First the three ways of
# resampling against none, both specifications, SSP2-4.5, income and share
# fixed; then 1 and 2 again with joint 5-year blocks.
N_PATHS <- 1000L
NOISE <- "blocks"
variability <- P$variability[["PWT | ERA5 | population"]]
with_noise <- function(spec, scenario, ssp, noise = NOISE, gmt = NULL) {
  project_setting(spec, SCEN, scenario, ssp, gmt = gmt, variability = variability, noise = noise, n_paths = N_PATHS)
}
noise_versions <- do.call(rbind, lapply(c(Linear = "Pref", Quadratic = "PrefQ"), function(k_model) {
  spec <- spec_of(models[[k_model]], "SSP2", FALSE)
  do.call(rbind, lapply(c("none", "iid", "joint", "blocks"), function(noise) {
    message("Noise versions: ", k_model, ", ", noise)
    r <- if (noise == "none") project_setting(spec, SCEN, "ssp245", "SSP2") else with_noise(spec, "ssp245", "SSP2", noise)
    g <- r$global[r$global$year %in% c(2050, 2100), ]
    country <- r$countries[r$countries$year == 2100, ]
    w <- SCEN$ssp$gdp[SCEN$ssp$gdp$ssp == "SSP2" & SCEN$ssp$gdp$year == 2100, ]
    w <- w$gdp[match(country$iso3, w$iso3)]
    data.frame(model = if (k_model == "PrefQ") "Quadratic arms" else "Linear arms", noise = noise, year = g$year,
               estimate = g$estimate, low = g$low, high = g$high,
               noise_low = g$noise_low %||% NA_real_, noise_high = g$noise_high %||% NA_real_,
               T = c(NA, sum(w * country$loss_T) / sum(w)), P = c(NA, sum(w * country$loss_P) / sum(w)))
  }))
}))
runs_noise <- list()
for (scenario in c("ssp126", "ssp245")) for (econ in names(ECON)) {
  ssp <- SCEN$ssp$of_scenario[[scenario]]
  message("Projecting with noise ", scenario, ", ", econ)
  runs_noise[[paste(scenario, econ)]] <- list(scenario = scenario, label = SCEN$scenarios[[scenario]], econ = econ,
                                              result = with_noise(spec_of(models$Pref, ssp, ECON[[econ]]), scenario, ssp))
}
global_noise <- do.call(rbind, lapply(runs_noise, function(x) {
  cbind(scenario = x$label, econ = ECON_LABELS[[x$econ]], x$result$global)
}))
countries_2100_noise <- do.call(rbind, lapply(runs_noise, function(x) {
  cbind(scenario = x$label, econ = ECON_LABELS[[x$econ]], x$result$countries[x$result$countries$year == 2100, ])
}))
termination_noise <- do.call(rbind, lapply(c(Quadratic = "PrefQ", Linear = "Pref"), function(k_model) {
  do.call(rbind, lapply(names(ECON), function(econ) {
    spec <- spec_of(models[[k_model]], "SSP2", ECON[[econ]])
    do.call(rbind, lapply(names(paths), function(k) {
      message("Termination with noise: ", k_model, ", ", econ, ", ", k)
      r <- with_noise(spec, "ssp245", "SSP2", gmt = paths[[k]])
      cbind(model = if (k_model == "PrefQ") "Quadratic arms" else "Linear arms", econ = ECON_LABELS[[econ]],
            path = k, r$global)
    }))
  }))
}))

saveRDS(list(global = global, countries_2100 = countries_2100, coverage = coverage, n = n_countries,
             followed = followed, termination = termination, gmt_paths = gmt_paths, held = HELD, tau = TAU,
             econ_labels = ECON_LABELS, noise_versions = noise_versions, global_noise = global_noise,
             countries_2100_noise = countries_2100_noise, termination_noise = termination_noise,
             n_noise = runs_noise[[1]]$result$n, n_paths = N_PATHS, built = Sys.time()),
        file.path(out_dir, "projections.rds"))
print(global[global$year %in% c(2050, 2100), ])
print(termination[termination$year %in% c(2050, 2070, 2100), ])
print(coverage); print(followed)
print(noise_versions)
print(global_noise[global_noise$year %in% c(2050, 2100), ])
print(termination_noise[termination_noise$year %in% c(2050, 2070, 2100), ])
