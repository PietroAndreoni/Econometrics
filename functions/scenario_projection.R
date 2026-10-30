# Projections of the deviation models of weather_arms_robustness.Rmd (saved to
# results/preferred_model/preferred_model.rds) on the SAI seed-ensemble
# scenarios of data/Scenarios: 4 RCPs x 9 SAI settings (no SAI, SAI never
# terminated, SAI terminated in 2040, 2050, ..., 2100) x 10 seeds of natural
# variability; 131 countries, the residual region ROW and WORLD, 1990-2150.
#
# In scenario s, country i and year t, with X the scenario temperature or
# precipitation, VAR its natural-variability component, and m the country's
# ERA5 mean over the estimation horizon (the long-run climate of every
# interaction):
#
#   horizon:  d = (X - m) / sd,          d_var = VAR / sd       (sd: ERA5, horizon)
#   20-year:  d = (X - mean20(X)) / s20,  d_var = (VAR - mean20(VAR)) / s20
#             (mean20: the mean of the scenario's own series over t-20..t-1;
#             that of VAR, which starts in 2015, over the years of the window
#             with it, at least VARIABILITY_MIN_YEARS). The baseline:
#             "fixed_sd"  s20 the country's ERA5 trailing 20-year SD, median
#                         over its estimation years (the default);
#             "moving"    s20 the SD of the series over t-20..t-1 (of VAR for
#                         d_var), as the mean;
#             "runs"      the moving mean and SD held through each rising
#                         (falling) run of the series, at the year before its
#                         first rise (fall), as run_fixed_deviation() of
#                         weather_arms_robustness.Rmd: the warm (wet) arm from
#                         the rising runs, the cold (dry) arm from the falling;
#             "sigma_spells"  the same, through the spells of
#                         fix_sigma_spells() (climate_runs.R): a spell ends
#                         with a reversal larger than the held SD or with the
#                         third consecutive smaller one.
#   unstandardized (the "..., unstandardized" definitions; C and m a year, no
#   SD, so the baseline does not apply):
#             horizon:  d = X - m,           d_var = VAR
#             20-year:  d = X - mean20(X),   d_var = VAR - mean20(VAR)
#   g(d)  = sum over the arms of (a + b m) f_arm(d)
#   dg    = g(d) - g(d_var)          growth effect (log points)
#   L_t   = sum of dg over start..t  effect on log GDP per capita, = c_t' beta
#   loss  = 1 - exp(L_t)             population unchanged, so also of GDP; for a
#                                    group of countries 1 - sum GDP exp(L) / sum GDP
#
# The per-capita growth of the scenario is that of its GDP and population
# series, g = dlog(GDP / population); with the climate, g + dg from `start`.
# Uncertainty: draws of the weather coefficients from N(beta, VCOV clustered by
# country), each applied to every seed (one pooled distribution).

SCENARIO_FILE <- "Scenarios/sai_seed_country_global_timeseries.parquet"
SCENARIO_METADATA_FILE <- "Scenarios/sai_seed_scenario_metadata.csv"
PREFERRED_MODEL_RDS <- file.path("results", "preferred_model", "preferred_model.rds")
# Trailing window of the 20-year deviations, and the fewest years of natural
# variability its moments are taken from (VAR starts in 2015, so the
# counterfactual of the 20-year deviations exists from 2025).
TRAILING_YEARS <- 20L
VARIABILITY_MIN_YEARS <- 10L
# Baseline of the 20-year deviations (see above): "fixed_sd", "moving" or
# "runs".
SCENARIO_BASELINE <- "fixed_sd"
SCENARIO_BASELINES <- c(fixed_sd = "trailing mean, fixed ERA5 SD",
                        moving = "moving trailing mean and SD",
                        runs = "trailing mean and SD held through runs",
                        sigma_spells = "trailing mean and SD held through one-SD spells")

# Short names of the columns of the scenario file.
SCENARIO_COLUMNS <- c(
  scenario = "scenario", rcp = "rcp", term = "term", seed = "seed", kind = "kind",
  termination_year = "termination_year", year = "year", region = "region",
  region_type = "region_type", T = "temperature_deg_c", P = "precipitation_lt_per_yr",
  TVAR = "DTEMP_REGION_VAR", PVAR = "DPRECIP_REGION_VAR", TATM_VAR = "TATM_VAR",
  gdp = "ykali_trillion_2005_usd_ppp", pop = "population_million"
)
SCENARIO_RCPS <- c(RCP3PD = "RCP2.6", RCP45 = "RCP4.5", RCP6 = "RCP6.0", RCP85 = "RCP8.5")
SCENARIO_TERMS <- c("BASE", "NOTERM", as.character(seq(2040L, 2100L, 10L)))

scenario_rcp_label <- function(rcp) unname(SCENARIO_RCPS[rcp])
scenario_term_label <- function(term) {
  ifelse(term == "BASE", "No SAI",
         ifelse(term == "NOTERM", "SAI, never terminated", paste("SAI, terminated", term)))
}
# Colours of the terms: no SAI dark, SAI never terminated blue, the
# terminations from yellow (earliest) to red (latest).
scenario_term_colours <- function() {
  terminations <- grDevices::colorRampPalette(c("#e6b800", "#e8702a", "#b2182b"))(
    length(SCENARIO_TERMS) - 2L)
  stats::setNames(c(GRID_INK, GRID_PALETTE[[1]], terminations),
                  scenario_term_label(SCENARIO_TERMS))
}

# Reading ------------------------------------------------------------------------

# The parquet's R metadata carries an external pointer that arrow discards with
# a warning; only that warning is muffled.
quiet_collect <- function(x) {
  withCallingHandlers(
    dplyr::collect(x),
    warning = function(w) {
      if (grepl("unsafe or invalid elements", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

scenario_dataset <- function(file = SCENARIO_FILE) arrow::open_dataset(data_file(file))

# The scenario matrix, one row per scenario.
scenario_matrix <- function(file = SCENARIO_METADATA_FILE) {
  utils::read.csv(data_file(file), stringsAsFactors = FALSE,
                  colClasses = c(term = "character")) %>%
    tibble::as_tibble() %>%
    select(scenario, rcp, term, seed, kind, termination_year)
}

# Rows of the scenarios with the given RCPs ("RCP3PD", "RCP45", "RCP6",
# "RCP85"), terms ("BASE", "NOTERM" or a termination year), seeds, regions
# (ISO3 codes, "ROW", "WORLD") and range of years; NULL keeps every value.
read_scenarios <- function(rcp = NULL, term = NULL, seed = NULL, regions = NULL,
                           years = NULL, file = SCENARIO_FILE) {
  ds <- scenario_dataset(file)
  if (!is.null(rcp)) ds <- filter(ds, rcp %in% !!as.character(rcp))
  if (!is.null(term)) ds <- filter(ds, term %in% !!as.character(term))
  if (!is.null(seed)) ds <- filter(ds, seed %in% !!as.integer(seed))
  if (!is.null(regions)) ds <- filter(ds, region %in% !!as.character(regions))
  if (!is.null(years)) ds <- filter(ds, year >= !!min(years), year <= !!max(years))
  out <- quiet_collect(select(ds, all_of(SCENARIO_COLUMNS)))
  if (!nrow(out)) stop("No scenario rows match the selection.")
  tibble::as_tibble(as.data.frame(out)) %>% arrange(scenario, region, year)
}

# Model --------------------------------------------------------------------------

# A model saved by weather_arms_robustness.Rmd: the deviation `definition` and
# arms `form` given (NULL: those of the preferred model), with the panel's
# horizon moments and its AIC against the best saved model. Interactions with
# the horizon mean only (the matched ones would need trailing climates).
load_preferred_model <- function(definition = NULL, form = NULL, file = PREFERRED_MODEL_RDS) {
  path <- if (file.exists(file)) file else file.path("..", file)
  if (!file.exists(path)) {
    stop("Knit weather_arms_robustness.Rmd first: ", file, " not found.")
  }
  saved <- readRDS(path)
  if (is.null(definition)) definition <- saved$preferred[["definition"]]
  if (is.null(form)) form <- saved$preferred[["form"]]
  model <- saved$models[[definition]][[form]]
  if (is.null(model)) stop("No saved model with definition '", definition, "' and form '", form, "'.")
  if (!all(model$suffix[c("clim_TM", "clim_RR")] == "all")) {
    stop("The projections support interactions with the horizon mean only.")
  }
  aic <- unlist(lapply(saved$models, function(x) vapply(x, `[[`, 0, "aic")))
  c(model, saved[c("panel", "observations", "horizon", "moments", "preferred")],
    list(aic_vs_best = model$aic - min(aic),
         label = paste0(definition, " deviations, ", form)))
}

# n draws of the weather coefficients from N(coef, vcov), one per column. The
# global random-number state is left as it was.
coefficient_draws <- function(model, n = 1000L, seed = 1L) {
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  set.seed(seed)
  z <- matrix(stats::rnorm(n * length(model$coef)), n, length(model$coef))
  draws <- t(z %*% chol(model$vcov)) + model$coef
  rownames(draws) <- names(model$coef)
  draws
}

# Value of an arm (column suffix): the warm (wet) side from deviations `up`,
# the cold (dry) side from `down`, each capped at `cap`.
arm_value <- function(suffix, up, down = up, cap = c(pos = Inf, neg = Inf)) {
  pos <- pmin(pmax(up, 0), cap[["pos"]])
  neg <- pmin(pmax(-down, 0), cap[["neg"]])
  if (suffix == "") return(pos - neg)
  switch(suffix, "_pos" = pos, "_neg" = neg, "_pos2" = pos^2, "_neg2" = neg^2,
         stop("Unknown arm: ", suffix))
}

# Cap of each side of variable v: its largest deviation in the estimation
# sample when `capped`, none otherwise.
arm_caps <- function(model, v, capped) {
  if (!capped) return(c(pos = Inf, neg = Inf))
  c(pos = max(model$support[[v]]), neg = -min(model$support[[v]]))
}

# Weather terms of `model` at deviations dev (list TM, RR, each a vector or a
# list(up, down) of the deviations of the warm and cold sides) and long-run
# climates clim (list TM, RR): one row per observation, one column per
# coefficient.
weather_design <- function(model, dev, clim, capped = FALSE) {
  terms <- names(model$coef)
  sides <- lapply(dev, function(x) if (is.list(x)) x else list(up = x, down = x))
  X <- matrix(0, length(sides$TM$up), length(terms), dimnames = list(NULL, terms))
  for (v in c("TM", "RR")) {
    cap <- arm_caps(model, v, capped)
    for (s in model$form_suffixes[[v]]) {
      term <- paste0("dev", v, s)
      interaction <- intersect(c(paste0(term, ":clim_", v), paste0("clim_", v, ":", term)), terms)
      stopifnot(term %in% terms, length(interaction) == 1)
      value <- arm_value(s, sides[[v]]$up, sides[[v]]$down, cap)
      X[, term] <- value
      X[, interaction] <- value * clim[[v]]
    }
  }
  X
}

# Effects ------------------------------------------------------------------------

# Mean and SD of x over the `window` years before each year (x of one unit, in
# consecutive years), from the years of the window with x; NA with fewer than
# `min_years` of them.
trailing_moments <- function(x, window = TRAILING_YEARS, min_years = window) {
  ok <- !is.na(x)
  z <- ifelse(ok, x, 0)
  t <- seq_along(x)
  window_sum <- function(v) {
    s <- c(0, cumsum(v))
    s[t] - s[pmax(t - window, 1L)]
  }
  n <- window_sum(ok)
  mean <- window_sum(z) / n
  sd <- sqrt(pmax((window_sum(z^2) - n * mean^2) / (n - 1), 0))
  mean[n < min_years] <- NA_real_
  sd[n < min_years] <- NA_real_
  list(mean = mean, sd = sd)
}

# trailing_moments() of each unit (`group`, rows in consecutive years).
grouped_trailing_moments <- function(x, group, min_years, window = TRAILING_YEARS) {
  mean <- sd <- rep(NA_real_, length(x))
  for (r in split(seq_along(x), group)) {
    m <- trailing_moments(x[r], window = window, min_years = min_years)
    mean[r] <- m$mean
    sd[r] <- m$sd
  }
  list(mean = mean, sd = sd)
}

# Deviation of x from its trailing moments (mean, sd), held through each
# rising (direction 1) or falling (-1) run of x at the run's first year, the
# year before the first rise (fall), as .run_ids() and run_fixed_deviation() of
# weather_arms_robustness.Rmd; outside runs, and in a run whose first year has
# no moments yet, the moving deviation. Rows of each unit (`group`) together,
# in consecutive years.
run_held_deviation <- function(x, mean, sd, group, direction) {
  n <- length(x)
  same_previous <- c(FALSE, group[-1] == group[-n])
  same_next <- c(same_previous[-1], FALSE)
  change <- ifelse(same_previous, x - c(NA_real_, x[-n]), NA_real_)
  rising <- !is.na(change) & direction * change > 0
  run <- ifelse(rising, cumsum(rising & !(c(FALSE, rising[-n]) & same_previous)), NA_integer_)
  # The year before the first rise belongs to the run.
  before <- is.na(run) & c(rising[-1], FALSE) & same_next
  run[before] <- c(run[-1], NA_integer_)[before]
  first <- match(run, run)
  held <- !is.na(run) & tabulate(run, nbins = max(c(run, 0L), na.rm = TRUE))[run] >= 2L
  held <- held & is.finite(mean[first]) & is.finite(sd[first])
  held[is.na(held)] <- FALSE
  (x - ifelse(held, mean[first], mean)) / ifelse(held, sd[first], sd)
}

# Deviations of variable v ("TM" or "RR") of the rows `d`, and of its
# natural-variability component, under the model's definition of v: each a
# list(up, down) of the deviations of the warm (wet) and cold (dry) sides,
# which differ only with baseline = "runs".
scenario_deviations <- function(d, model, v, group, baseline = SCENARIO_BASELINE) {
  both <- function(x) list(up = x, down = x)
  x <- d[[c(TM = "T", RR = "P")[[v]]]]
  variability <- d[[c(TM = "TVAR", RR = "PVAR")[[v]]]]
  if (model$suffix[[v]] == "all") {
    sd <- d[[paste0("sd_", v, "_all")]]
    return(list(dev = both((x - d[[paste0("mean_", v, "_all")]]) / sd),
                dev_var = both(variability / sd)))
  }
  # Unstandardized deviations (physical units): as estimated, whatever the
  # baseline.
  if (model$suffix[[v]] == "all_unstd") {
    return(list(dev = both(x - d[[paste0("mean_", v, "_all")]]), dev_var = both(variability)))
  }
  if (model$suffix[[v]] == paste0("lag_", TRAILING_YEARS, "_unstd")) {
    own <- grouped_trailing_moments(x, group, TRAILING_YEARS)
    counterfactual <- grouped_trailing_moments(variability, group, VARIABILITY_MIN_YEARS)
    return(list(dev = both(x - own$mean), dev_var = both(variability - counterfactual$mean)))
  }
  if (model$suffix[[v]] == paste0("lag_", TRAILING_YEARS, "_sd_all")) {
    # The 20-year mean in the horizon SD: as estimated, whatever the baseline.
    sd <- d[[paste0("sd_", v, "_all")]]
    own <- grouped_trailing_moments(x, group, TRAILING_YEARS)
    counterfactual <- grouped_trailing_moments(variability, group, VARIABILITY_MIN_YEARS)
    return(list(dev = both((x - own$mean) / sd),
                dev_var = both((variability - counterfactual$mean) / sd)))
  }
  stopifnot(model$suffix[[v]] == paste0("lag_", TRAILING_YEARS),
            baseline %in% names(SCENARIO_BASELINES))
  own <- grouped_trailing_moments(x, group, TRAILING_YEARS)
  counterfactual <- grouped_trailing_moments(variability, group, VARIABILITY_MIN_YEARS)
  if (baseline == "fixed_sd") {
    sd <- d[[paste0("sd_", v, "_lag_", TRAILING_YEARS)]]
    if (is.null(sd)) stop("Re-knit weather_arms_robustness.Rmd: the saved model has no trailing SDs.")
    return(list(dev = both((x - own$mean) / sd),
                dev_var = both((variability - counterfactual$mean) / sd)))
  }
  if (baseline == "moving") {
    return(list(dev = both((x - own$mean) / own$sd),
                dev_var = both((variability - counterfactual$mean) / counterfactual$sd)))
  }
  held <- if (baseline == "runs") {
    function(z, m) {
      list(up = run_held_deviation(z, m$mean, m$sd, group, 1),
           down = run_held_deviation(z, m$mean, m$sd, group, -1))
    }
  } else {
    function(z, m) {
      list(up = sigma_spell_held_deviation(z, m$mean, m$sd, d$year, group, 1),
           down = sigma_spell_held_deviation(z, m$mean, m$sd, d$year, group, -1))
    }
  }
  list(dev = held(x, own), dev_var = held(variability, counterfactual))
}

# fix_sigma_spells() deviation of each unit (`group`, rows in consecutive
# years).
sigma_spell_held_deviation <- function(x, mean, sd, year, group, direction) {
  out <- rep(NA_real_, length(x))
  for (r in split(seq_along(x), group)) {
    out[r] <- fix_sigma_spells(x[r], mean[r], sd[r], year[r], direction)$deviation
  }
  out
}

# Label of a model and baseline for titles.
projection_label <- function(model, baseline = SCENARIO_BASELINE) {
  paste0(model$label, if (any(model$suffix[c("TM", "RR")] == paste0("lag_", TRAILING_YEARS))) {
    paste0(" (", SCENARIO_BASELINES[[baseline]], ")")
  })
}

# Every country-year of `data` (read_scenarios(), with the years before
# `start` that the trailing moments need) with the deviations, the growth
# effect dg from the first year with both deviations, the effect on log GDP
# per capita L from `start` (point estimate and standard error for the seed),
# the GDP loss, the per-capita growth of the scenario with and without dg, and
# the cumulative design (matrix column `design`, L = design %*% coef) for the
# draws. Countries without moments (ROW) are left out.
scenario_effects <- function(data, model, start = 2025L, capped = FALSE,
                             baseline = SCENARIO_BASELINE) {
  d <- data %>%
    filter(region %in% model$moments$GID_0) %>%
    arrange(scenario, region, year) %>%
    group_by(scenario, region) %>%
    mutate(g_base = ifelse(year - lag(year) == 1L, log(gdp / pop) - lag(log(gdp / pop)),
                           NA_real_)) %>%
    ungroup() %>%
    inner_join(model$moments, by = c(region = "GID_0"))
  group <- paste(d$scenario, d$region)
  same <- group[-1] == group[-length(group)]
  stopifnot(all(diff(d$year)[same] == 1L))
  deviations <- list(dT = scenario_deviations(d, model, "TM", group, baseline),
                     dP = scenario_deviations(d, model, "RR", group, baseline))
  # Each side's deviation (_up, _down), and the signed deviation of the arms
  # (warm side minus cold side) for the summaries.
  for (k in names(deviations)) for (part in c("dev", "dev_var")) {
    name <- paste0(k, if (part == "dev_var") "_var")
    sides <- deviations[[k]][[part]]
    d[[paste0(name, "_up")]] <- sides$up
    d[[paste0(name, "_down")]] <- sides$down
    d[[name]] <- pmax(sides$up, 0) - pmax(-sides$down, 0)
  }
  finite <- Reduce(`&`, lapply(paste0(rep(c("dT", "dT_var", "dP", "dP_var"), each = 2),
                                      c("_up", "_down")), function(k) is.finite(d[[k]])))
  d <- d[finite, ]
  if (!nrow(d) || start < min(d$year)) {
    stop("The growth effect starts in ", if (nrow(d)) min(d$year) else "no year",
         "; `start` must not be earlier.")
  }
  group <- paste(d$scenario, d$region)
  clim <- list(TM = d$mean_TM_all, RR = d$mean_RR_all)
  side <- function(k) list(up = d[[paste0(k, "_up")]], down = d[[paste0(k, "_down")]])
  X <- weather_design(model, list(TM = side("dT"), RR = side("dP")), clim, capped) -
    weather_design(model, list(TM = side("dT_var"), RR = side("dP_var")), clim, capped)
  cumulative <- apply(X * (d$year >= start), 2, function(x) stats::ave(x, group, FUN = cumsum))
  cumulative <- matrix(cumulative, nrow(X), dimnames = dimnames(X))

  d$dg <- drop(X %*% model$coef)
  d$L <- drop(cumulative %*% model$coef)
  d$L_se <- sqrt(pmax(rowSums((cumulative %*% model$vcov) * cumulative), 0))
  d$loss <- 1 - exp(d$L)
  d$g_cc <- d$g_base + ifelse(d$year >= start, d$dg, 0)
  d$gdp_cc <- d$gdp * exp(d$L)
  d$capped <- capped
  d$start <- as.integer(start)
  d$baseline <- baseline
  d$design <- cumulative
  d
}

# GDP loss of the countries `regions` (NULL: every country of `effects`) from
# the start year: the point estimate of each scenario in every year (`seeds`),
# and its distribution over the draws (columns of `draws`) and the seeds of
# each RCP, term and cap pooled, in `draw_years` (NULL: every year; `summary`:
# mean and the quantiles `probs`, the middle one the central line). The draws
# cost one exp() per country, year and draw, hence the coarser years.
gdp_loss <- function(effects, model, draws, regions = NULL, probs = c(0.05, 0.5, 0.95),
                     draw_years = NULL) {
  e <- effects[effects$year >= effects$start, ]
  if (!is.null(regions)) e <- e[e$region %in% regions, ]
  if (!nrow(e)) stop("No country-years for the regions selected.")
  b <- model$coef[rownames(draws)]
  runs <- split(seq_len(nrow(e)), paste(e$scenario, e$capped))
  per_run <- lapply(runs, function(r) {
    year <- e$year[r]
    gdp <- e$gdp[r]
    total <- rowsum(gdp, year)
    point <- 1 - rowsum(gdp * exp(drop(e$design[r, , drop = FALSE] %*% b)), year) / total
    if (!is.null(draw_years)) r <- r[year %in% draw_years]
    drawn <- 1 - rowsum(e$gdp[r] * exp(e$design[r, , drop = FALSE] %*% draws), e$year[r]) /
      drop(rowsum(e$gdp[r], e$year[r]))
    list(key = e[r[1], c("rcp", "term", "seed", "capped")],
         years = as.integer(rownames(total)), point = drop(point),
         draw_years = as.integer(rownames(drawn)), drawn = drawn)
  })
  seeds <- bind_rows(lapply(per_run, function(x) {
    tibble::tibble(x$key, year = x$years, loss = x$point)
  }))
  pooled <- split(per_run, vapply(per_run, function(x) paste(x$key$rcp, x$key$term, x$key$capped), ""))
  summary <- bind_rows(lapply(pooled, function(runs) {
    years <- Reduce(intersect, lapply(runs, `[[`, "draw_years"))
    m <- do.call(cbind, lapply(runs, function(x) x$drawn[as.character(years), , drop = FALSE]))
    q <- t(apply(m, 1, stats::quantile, probs = probs, names = FALSE))
    tibble::tibble(runs[[1]]$key[, c("rcp", "term", "capped")], year = years,
                   seeds = length(runs), mean = rowMeans(m),
                   low = q[, 1], central = q[, ceiling(length(probs) / 2)],
                   high = q[, length(probs)])
  }))
  list(summary = summary, seeds = seeds)
}

# GDP loss of each country in `years`, its distribution pooled over the draws
# and the seeds of each RCP, term and cap, with the mean annual per-capita
# growth of the scenario over start..growth_end with and without the climate
# (mean over the seeds, point estimate).
country_loss <- function(effects, draws, years, probs = c(0.05, 0.5, 0.95),
                         growth_end = 2100L) {
  e <- effects[effects$year %in% years, ]
  drawn <- 1 - exp(e$design %*% draws)
  keys <- paste(e$rcp, e$term, e$capped, e$region, e$year)
  rows <- split(seq_len(nrow(e)), keys)
  out <- bind_rows(lapply(rows, function(r) {
    q <- stats::quantile(drawn[r, ], probs = probs, names = FALSE)
    tibble::tibble(e[r[1], c("rcp", "term", "capped", "region", "year")],
                   seeds = length(r), point = mean(e$loss[r]), mean = mean(drawn[r, ]),
                   low = q[[1]], central = q[[ceiling(length(probs) / 2)]],
                   high = q[[length(probs)]])
  }))
  growth <- effects %>%
    filter(year >= start, year <= growth_end) %>%
    group_by(rcp, term, capped, region, seed) %>%
    summarise(g_base = mean(g_base), g_cc = mean(g_cc), .groups = "drop_last") %>%
    summarise(g_base = mean(g_base), g_cc = mean(g_cc), .groups = "drop")
  left_join(out, growth, by = c("rcp", "term", "capped", "region"))
}

# Climate ------------------------------------------------------------------------

# Temperature (C) or precipitation (m a year) of the countries `regions` (NULL:
# every country; population-weighted mean), or for "WORLD" the global
# temperature of the model (TATM), in each scenario and year: the series, and
# the series without its natural-variability component (`forced`).
scenario_climate <- function(data, variable = c("temperature", "precipitation"),
                             regions = NULL) {
  variable <- match.arg(variable)
  v <- if (variable == "temperature") "T" else "P"
  if (identical(regions, "WORLD")) {
    if (v == "P") stop("WORLD has no precipitation; select countries (NULL: all).")
    return(data %>%
             filter(region == "WORLD") %>%
             transmute(scenario, rcp, term, seed, year, value = T,
                       forced = T - dplyr::coalesce(TATM_VAR, 0)))
  }
  variability <- if (v == "T") "TVAR" else "PVAR"
  d <- filter(data, region_type == "country")
  if (!is.null(regions)) d <- filter(d, region %in% regions)
  d %>%
    group_by(scenario, rcp, term, seed, year) %>%
    summarise(value = stats::weighted.mean(.data[[v]], pop),
              forced = stats::weighted.mean(.data[[v]] - dplyr::coalesce(.data[[variability]], 0), pop),
              .groups = "drop")
}

# Selection and plots ------------------------------------------------------------

# Label of a group of regions for titles.
regions_label <- function(regions, n_countries = NULL) {
  if (is.null(regions)) {
    return(sprintf("All %s countries (ROW left out)", if (is.null(n_countries)) "" else n_countries))
  }
  if (identical(regions, "WORLD")) return("World: the model's global temperature (TATM)")
  if (length(regions) <= 6) return(paste(regions, collapse = ", "))
  sprintf("%d countries (%s, ...)", length(regions), paste(utils::head(regions, 4), collapse = ", "))
}

# The scenarios with the given RCPs and terms (and seeds; NULL: all ten), and
# for the regions `regions` (NULL: every country) either
#   - "gdp_loss": the GDP loss from `start` (share of the GDP of the scenario),
#     its central value and band pooled over the draws and the seeds (every
#     `draw_step` years), and the point estimate of each seed (every year);
#   - "temperature" or "precipitation": the population-weighted series of each
#     seed, and the mean over the seeds of the series without natural
#     variability ("WORLD": the model's global temperature).
# `baseline` ("fixed_sd", "moving" or "runs") applies to the 20-year deviations.
# Returns list(summary, seeds, what, regions, n_countries, capped, start, probs,
# model).
select_scenarios <- function(rcp, term, what = c("gdp_loss", "temperature", "precipitation"),
                             regions = NULL, seed = NULL, capped = FALSE, start = 2025L,
                             years = c(1990L, 2150L), model = load_preferred_model(),
                             draws = coefficient_draws(model), probs = c(0.05, 0.5, 0.95),
                             baseline = SCENARIO_BASELINE, draw_step = 5L) {
  what <- match.arg(what)
  if (what == "gdp_loss") {
    # (From the first year: the trailing moments need the years before `start`.)
    data <- read_scenarios(rcp, term, seed, regions, years = c(1990L, max(years)))
    effects <- scenario_effects(data, model, start = start, capped = capped,
                                baseline = baseline)
    draw_years <- unique(c(seq(start, max(years), by = draw_step), max(years)))
    loss <- gdp_loss(effects, model, draws, probs = probs, draw_years = draw_years)
    out <- list(summary = loss$summary %>% rename(value = central),
                seeds = loss$seeds %>% rename(value = loss),
                n_countries = length(unique(effects$region)))
  } else {
    data <- read_scenarios(rcp, term, seed, regions, years = years)
    climate <- scenario_climate(data, what, regions)
    out <- list(
      summary = climate %>%
        group_by(rcp, term, year) %>%
        summarise(seeds = n(), value = mean(forced), low = NA_real_, high = NA_real_,
                  .groups = "drop"),
      seeds = climate %>% select(rcp, term, seed, year, value),
      n_countries = length(unique(data$region[data$region_type == "country"]))
    )
  }
  list(summary = label_scenarios(out$summary), seeds = label_scenarios(out$seeds), what = what,
       regions = regions, n_countries = out$n_countries, capped = capped, start = start,
       probs = probs, model = if (what == "gdp_loss") projection_label(model, baseline))
}

# RCP and term labels, as ordered factors.
label_scenarios <- function(x) {
  mutate(x, rcp_label = factor(scenario_rcp_label(rcp), levels = unname(SCENARIO_RCPS)),
         term_label = factor(scenario_term_label(term), levels = scenario_term_label(SCENARIO_TERMS)))
}

# The global GDP loss of project_all_scenarios() (`projections`) for the RCPs
# and terms given (NULL: all), uncapped or capped, as a selection for
# plot_scenarios().
global_selection <- function(projections, rcp = NULL, term = NULL, capped = FALSE,
                             start = 2025L, probs = c(0.05, 0.5, 0.95), n_countries = NULL,
                             model = NULL, baseline = SCENARIO_BASELINE) {
  keep <- function(x) {
    x[x$capped == capped & (is.null(rcp) | x$rcp %in% rcp) & (is.null(term) | x$term %in% term), ]
  }
  list(summary = label_scenarios(rename(keep(projections$global), value = central)),
       seeds = label_scenarios(rename(keep(projections$global_seeds), value = loss)),
       what = "gdp_loss", regions = NULL, n_countries = n_countries, capped = capped,
       start = start, probs = probs, model = if (!is.null(model)) projection_label(model, baseline))
}

# Plot of select_scenarios() (or, without `selection`, of its arguments `...`):
# one facet per RCP, one colour per term; the GDP loss with its pooled band
# and thin lines for the seeds (bands of the terms `bands` only; NULL: all).
# The y range of the loss (`ylim`, NULL: the lines, and the band down to -100%
# or the lowest line) keeps the gains of the extreme draws from flattening the
# rest.
plot_scenarios <- function(selection = NULL, ylim = NULL, bands = NULL, ...) {
  s <- if (is.null(selection)) select_scenarios(...) else selection
  band_data <- if (is.null(bands)) s$summary else s$summary[s$summary$term %in% bands, ]
  colours <- scenario_term_colours()
  colours <- colours[names(colours) %in% c(levels(droplevels(s$summary$term_label)))]
  y_label <- switch(s$what,
                    gdp_loss = "GDP loss (% of the scenario's GDP)",
                    temperature = if (identical(s$regions, "WORLD")) {
                      "Global temperature above pre-industrial (C)"
                    } else "Temperature (C)",
                    precipitation = "Precipitation (m a year)")
  scale <- if (s$what == "gdp_loss") 100 else 1
  band <- if (s$what == "gdp_loss") {
    lines <- 100 * c(s$summary$value, s$seeds$value)
    floor <- min(-100, lines)
    if (is.null(ylim)) {
      ylim <- c(max(100 * min(band_data$low), floor), max(100 * band_data$high, lines))
    }
    sprintf("Line: median; band: %g-%g%% of the %s seeds x coefficient draws pooled%s; thin lines: each seed (point estimate).",
            100 * s$probs[[1]], 100 * s$probs[[length(s$probs)]], max(s$summary$seeds),
            if (100 * min(band_data$low) < ylim[[1]]) sprintf(" (cut at %.0f%%)", ylim[[1]]) else "")
  } else {
    "Thin lines: each seed; thick line: mean over the seeds without the natural-variability component."
  }
  title <- switch(s$what, gdp_loss = "GDP loss", temperature = "Temperature",
                  precipitation = "Precipitation")
  subtitle <- paste0(
    regions_label(s$regions, s$n_countries),
    if (s$what == "gdp_loss") {
      paste0("; loss accumulated from ", s$start, if (!is.null(s$model)) paste0("; ", s$model),
             if (s$capped) "; deviations capped at the sample support" else "")
    } else if (is.null(s$regions) || length(s$regions) > 1 && !identical(s$regions, "WORLD")) {
      "; population-weighted mean"
    } else ""
  )
  ggplot(mapping = aes(year, scale * value, colour = term_label, fill = term_label)) +
    (if (s$what == "gdp_loss") geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3)) +
    (if (s$what == "gdp_loss") {
      geom_ribbon(data = band_data, aes(ymin = scale * low, ymax = scale * high), alpha = 0.12,
                  colour = NA)
    }) +
    geom_line(data = s$seeds, aes(group = interaction(term_label, seed)), linewidth = 0.2,
              alpha = 0.35) +
    geom_line(data = s$summary, linewidth = 0.9) +
    facet_wrap(~rcp_label) +
    coord_cartesian(ylim = ylim) +
    scale_colour_manual(values = colours, name = NULL, drop = TRUE) +
    # (One legend, by colour: the bands may cover only some of the terms.)
    scale_fill_manual(values = colours, name = NULL, drop = TRUE, guide = "none") +
    labs(x = NULL, y = y_label, title = title, subtitle = subtitle, caption = band) +
    theme_results()
}

# All scenarios ------------------------------------------------------------------

# Every RCP and term of the scenario matrix, uncapped and capped: the global GDP
# loss (gdp_loss() over every country; the draws every `draw_step` years), the
# country losses in `country_years`, the share of the country-years from
# `start` with a deviation beyond the sample support, and the deviations
# without natural variability (d - d_var) in `country_years`. One RCP and term
# (ten seeds) at a time, on `workers` R processes.
project_all_scenarios <- function(model, draws, start = 2025L, capped = c(FALSE, TRUE),
                                  country_years = c(2050L, 2100L, 2150L),
                                  probs = c(0.05, 0.5, 0.95), draw_step = 5L,
                                  scenarios = scenario_matrix(), workers = 1L,
                                  baseline = SCENARIO_BASELINE) {
  # (Forced here: the workers receive run_group() with its environment, and an
  # unevaluated argument would be looked up in their global environment.)
  invisible(list(model, draws, capped, country_years, probs, start, baseline))
  groups <- distinct(scenarios, rcp, term)
  draw_years <- unique(c(seq(start, 2150L, by = draw_step), 2150L))
  run_group <- function(i) {
    data <- read_scenarios(groups$rcp[i], groups$term[i])
    lapply(capped, function(cap) {
      e <- scenario_effects(data, model, start = start, capped = cap, baseline = baseline)
      after <- e$year >= start
      support <- tibble::tibble(
        rcp = groups$rcp[i], term = groups$term[i], capped = cap,
        hot = mean(e$dT_up[after] > max(model$support$TM)),
        cold = mean(e$dT_down[after] < min(model$support$TM)),
        wet = mean(e$dP_up[after] > max(model$support$RR)),
        dry = mean(e$dP_down[after] < min(model$support$RR))
      )
      deviations <- e %>%
        filter(year %in% country_years) %>%
        group_by(rcp, term, year) %>%
        summarise(T_median = stats::median(dT - dT_var), T_p90 = stats::quantile(dT - dT_var, 0.9),
                  T_max = max(dT - dT_var), P_p10 = stats::quantile(dP - dP_var, 0.1),
                  P_median = stats::median(dP - dP_var), P_p90 = stats::quantile(dP - dP_var, 0.9),
                  .groups = "drop")
      loss <- gdp_loss(e, model, draws, probs = probs, draw_years = draw_years)
      list(global = loss$summary, global_seeds = loss$seeds,
           country = country_loss(e, draws, country_years, probs), support = support,
           deviations = if (!cap) deviations)
    })
  }
  results <- if (workers > 1L) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster))
    parallel::clusterCall(cluster, function(dir) {
      setwd(dir)
      suppressMessages(source("load_functions.R"))
      NULL
    }, getwd())
    parallel::parLapplyLB(cluster, seq_len(nrow(groups)), run_group)
  } else {
    lapply(seq_len(nrow(groups)), run_group)
  }
  flat <- unlist(results, recursive = FALSE)
  parts <- c("global", "global_seeds", "country", "support", "deviations")
  lapply(stats::setNames(nm = parts), function(x) bind_rows(lapply(flat, `[[`, x)))
}

# A result saved as rds with the specification that produced it, rebuilt when
# the specification changes (as cached_panel()).
cached_result <- function(path, spec, build) {
  spec_path <- sub("\\.rds$", ".spec.rds", path)
  if (file.exists(path) && file.exists(spec_path) && identical(readRDS(spec_path), spec)) {
    message("Reusing ", path)
    return(readRDS(path))
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  result <- build()
  saveRDS(result, path)
  saveRDS(spec, spec_path)
  result
}
