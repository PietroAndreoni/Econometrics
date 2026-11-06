# The linear interaction model (M2b: warm/cold x climate x income, wet/dry x
# climate x share) with four definitions of the weather shocks, PWT x ERA5
# (the handout's rows and outliers), each with country quadratic trends,
# linear trends and none (year effects always):
#   (a) default: x = (v - mu20) / sigma, demeaned by the country's trend of
#       the regression (quadratic, linear or mean only), then split into arms;
#   (b) no detrending before the split: the same x, not demeaned (the arms'
#       trends only absorbed by the regression's effects and trends);
#   (c) the nonlinear autoregressive filter of Nath, Ramey and Klenow (Eq. 4):
#       pooled over countries, v_it = sum_j g_j v_i,t-j + sum_j h_j v_i,t-j vbar_i
#       + mu_i + tau_it (j = 1..3; vbar_i the country's mean over the horizon;
#       country effects, no year effects); the shock is tau over its SD in the
#       country over the horizon, split into arms;
#   (d) as (c), with growth filtered in the same way (its own lags times the
#       country's mean growth), the innovation as the outcome.
# All on the rows every option has (the lags of (d) need the three previous
# years of growth). Errors clustered by country.
#   Rscript shock_filters.R    -> results/handout_preferred_specification/shock_filters/

suppressPackageStartupMessages({ library(ggplot2) })
Sys.setenv(HANDOUT_GDP = "PWT", HANDOUT_ROWS = "own", HANDOUT_OUTLIERS = "5", HANDOUT_MIN_INCOME = "0")
main <- readLines("handout_preferred_specification.R")
eval(parse(text = main[seq_len(grep("models <- list(", main, fixed = TRUE)[1] - 1)]))
out_dir <- file.path("results", "handout_preferred_specification", "shock_filters")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
LAGS <- 3
FILTER_START <- 1950
TRENDS <- c(quadratic = "year + GID_1[t, t2]", linear = "year + GID_1[t]", none = "year + GID_1")
OPTIONS <- c(a = "(a) default", b = "(b) not detrended", c = "(c) NRK filter", d = "(d) NRK filter, growth too")
FIG_COUNTRIES <- c("MOZ", "ZWE", "MAR", "NER", "GNB", "SYR", "IND", "CHN", "USA", "BRA")

# All the panel's rows (outliers included, for the lags of growth) and the
# annual climate of the same product and weights.
panel <- P$panels[["PWT | ERA5 | population"]]$rows
horizon <- range(panel$year)
raw <- as.data.frame(arrow::read_parquet(file.path("data", "cache", "climate_gadm0_era_population_density_2000.parquet")))
raw <- raw[raw$GID_0 %in% panel$GID_0 & raw$year >= FILTER_START, c("GID_0", "year", "TM", "RR")]
key <- function(id, year) paste(id, year)
# Check: the cache's T less the panel's trailing mean is the panel's deviation
# up to a constant per country (the panel's deviation is demeaned).
chk <- merge(panel[c("GID_0", "year", "mean_TM_lag_20", "dm_TM_lag_20_unstd_unit")], raw, by = c("GID_0", "year"))
gap <- with(chk, (TM - mean_TM_lag_20) - dm_TM_lag_20_unstd_unit)
stopifnot(max(tapply(gap, chk$GID_0, function(x) diff(range(x)))) < 1e-6)

# The filter: pooled regression of v on its LAGS lags and the lags times the
# country's mean (state), with country effects; fitted where `use`, the lags
# from the whole series. Returns the residuals and the fit.
lagged <- function(df, v, k) df[[v]][match(key(df$GID_0, df$year - k), key(df$GID_0, df$year))]
nrk_filter <- function(df, v, state, use = rep(TRUE, nrow(df))) {
  lags <- paste0("l", seq_len(LAGS))
  for (k in seq_len(LAGS)) df[[lags[k]]] <- lagged(df, v, k)
  df$state <- state
  df$use <- use
  f <- as.formula(paste(v, "~", paste(c(lags, paste0(lags, ":state")), collapse = " + "), "| GID_0"))
  m <- feols(f, df, subset = ~use, notes = FALSE)
  res <- rep(NA_real_, nrow(df))
  res[obs(m)] <- resid(m)
  list(resid = res, fit = m)
}
in_horizon <- raw$year >= horizon[1] & raw$year <= horizon[2]
vbar <- function(v) ave(ifelse(in_horizon, raw[[v]], NA), raw$GID_0, FUN = function(x) mean(x, na.rm = TRUE))
filters <- list()
# The size of one SD of the innovation against the SD of the default
# deviations (sigma of v over the horizon), median over the countries.
sd_ratio <- numeric(0)
for (v in c("TM", "RR")) {
  fl <- nrk_filter(raw, v, vbar(v))
  filters[[v]] <- fl$fit
  sd_i <- ave(ifelse(in_horizon, fl$resid, NA), raw$GID_0, FUN = function(x) sd(x, na.rm = TRUE))
  raw[[paste0("tau_", v)]] <- fl$resid / sd_i
  sigma <- panel[[paste0("sd_", v, "_all")]][match(raw$GID_0, panel$GID_0)]
  sd_ratio[v] <- median(tapply(sd_i / sigma, raw$GID_0, mean), na.rm = TRUE)
}
# Growth: the same filter (state: the country's mean growth over its rows that
# are not outliers), fitted on those rows.
panel$g_all <- panel[[P$outcome]]
keep <- !panel$outlier_5 & is.finite(panel$g_all)
panel$gbar <- ave(ifelse(keep, panel$g_all, NA), panel$GID_0, FUN = function(x) mean(x, na.rm = TRUE))
fg <- nrk_filter(panel, "g_all", panel$gbar, keep)
filters$growth <- fg$fit
panel$g_innov <- fg$resid

# The shocks on the estimation rows.
i_raw <- match(key(d$GID_0, d$year), key(raw$GID_0, raw$year))
i_pan <- match(key(d$GID_0, d$year), key(panel$GID_0, panel$year))
for (v in c("TM", "RR")) {
  d[[paste0("x_", v)]] <- (raw[[v]][i_raw] - d[[paste0("mean_", v, "_lag_20")]]) / d[[paste0("sd_", v, "_all")]]
  d[[paste0("tau_", v)]] <- raw[[paste0("tau_", v)]][i_raw]
}
d$g_innov <- panel$g_innov[i_pan]
rows_before <- nrow(d)
d <- d[is.finite(d$x_TM) & is.finite(d$x_RR) & is.finite(d$tau_TM) & is.finite(d$tau_RR) & is.finite(d$g_innov), ]
message("common rows: ", nrow(d), " of ", rows_before, "; countries ", length(unique(d$GID_0)))

shock_of <- function(option, v, trend) {
  switch(option,
         a = d[[paste0("dm_", v, "_lag_20_unstd_", c(quadratic = "quadratic", linear = "linear", none = "unit")[[trend]])]] /
           d[[paste0("sd_", v, "_all")]],
         b = d[[paste0("x_", v)]],
         c = , d = d[[paste0("tau_", v)]])
}
m2b_terms <- c(arm_terms("TM", c(CLIM[["TM"]], "y")), arm_terms("RR", c(CLIM[["RR"]], "agr")))
fit_option <- function(option, trend) {
  data <- d
  for (v in c("TM", "RR")) {
    x <- shock_of(option, v, trend)
    data[[paste0("p", v)]] <- pmax(x, 0)
    data[[paste0("n", v)]] <- pmax(-x, 0)
  }
  data$outcome <- if (option == "d") data$g_innov else data$g
  for (k in names(m2b_terms)) data[[k]] <- Reduce(`*`, lapply(m2b_terms[[k]], function(z) data[[z]]))
  f <- as.formula(paste("outcome ~", paste(c(names(m2b_terms), LEVELS, CONTROLS), collapse = " + "), "|", TRENDS[[trend]]))
  m <- feols(f, data, cluster = ~GID_0, notes = FALSE)
  m$country <- data$GID_0[obs(m)]
  list(model = m, data = data)
}
fits <- list()
for (o in names(OPTIONS)) for (tr in names(TRENDS)) fits[[paste(o, tr)]] <- fit_option(o, tr)

coefs <- do.call(rbind, lapply(names(fits), function(k) {
  ct <- coeftable(fits[[k]]$model)[names(m2b_terms), ]
  data.frame(option = sub(" .*", "", k), trend = sub(".* ", "", k), term = gsub("_X_", ":", rownames(ct)),
             estimate = 100 * ct[, 1], se = 100 * ct[, 2], p = ct[, 4], row.names = NULL)
}))
of <- function(v, z) paste0(c("p", "n"), v, "_X_", z)
tests <- do.call(rbind, lapply(names(fits), function(k) {
  m <- fits[[k]]$model
  data.frame(option = sub(" .*", "", k), trend = sub(".* ", "", k),
             `T terms` = zero_test(m, grep("TM", names(m2b_terms), value = TRUE))[["p"]],
             `T x income` = zero_test(m, of("TM", "y"))[["p"]], `T x climate` = zero_test(m, of("TM", "clim_TM"))[["p"]],
             `P terms` = zero_test(m, grep("RR", names(m2b_terms), value = TRUE))[["p"]],
             `P x share` = zero_test(m, of("RR", "agr"))[["p"]], `P x climate` = zero_test(m, of("RR", "clim_RR"))[["p"]],
             rows = nobs(m), within_r2 = fitstat(m, "wr2")[[1]], check.names = FALSE)
}))
# Where the warm arm sits in time: mean warm and cold arm by decade (quadratic
# trends for (a)).
by_decade <- do.call(rbind, lapply(c("a", "b", "c"), function(o) {
  x <- shock_of(o, "TM", "quadratic")
  dec <- 10 * floor(d$year / 10)
  data.frame(option = o, decade = sort(unique(dec)), warm = tapply(pmax(x, 0), dec, mean), cold = tapply(pmax(-x, 0), dec, mean),
             share_warm = tapply(x > 0, dec, mean), row.names = NULL)
}))
filter_coefs <- lapply(filters, function(m) coeftable(m)[, 1:2])
# Persistence implied by the filter at the 10th, 50th and 90th percentile of
# the state: sum of the lag coefficients.
persistence <- t(sapply(names(filters), function(v) {
  b <- coef(filters[[v]])
  st <- if (v == "growth") panel$gbar[keep] else vbar(v)[in_horizon]
  q <- quantile(unique(st), c(0.1, 0.5, 0.9), na.rm = TRUE)
  sapply(q, function(s) sum(b[paste0("l", seq_len(LAGS))]) + s * sum(b[paste0("l", seq_len(LAGS), ":state")]))
}))
saveRDS(list(coefs = coefs, tests = tests, by_decade = by_decade, filters = filter_coefs, persistence = persistence, sd_ratio = sd_ratio,
             rows = nrow(d), countries = length(unique(d$GID_0)), options = OPTIONS, trends = TRENDS),
        file.path(out_dir, "shock_filters.rds"))

# Figures. Colours: one per option, in a fixed order (Okabe-Ito).
OPTION_COLOURS <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7"), OPTIONS)
cname <- function(iso) countrycode::countrycode(iso, "iso3c", "country.name")

# 1. The shock series of each option, T and P, for the figure's countries
# (all the panel's years; (a) with quadratic trends).
series <- function(v) {
  p <- panel[panel$GID_0 %in% FIG_COUNTRIES, ]
  j <- match(key(p$GID_0, p$year), key(raw$GID_0, raw$year))
  sdv <- p[[paste0("sd_", v, "_all")]]
  rbind(data.frame(country = p$GID_0, year = p$year, option = OPTIONS[["a"]], x = p[[paste0("dm_", v, "_lag_20_unstd_quadratic")]] / sdv),
        data.frame(country = p$GID_0, year = p$year, option = OPTIONS[["b"]], x = (raw[[v]][j] - p[[paste0("mean_", v, "_lag_20")]]) / sdv),
        data.frame(country = p$GID_0, year = p$year, option = OPTIONS[["c"]], x = raw[[paste0("tau_", v)]][j]))
}
# What the model sees: a column of the estimation rows less the year effects
# and the country quadratic trends of the main regression (all countries; the
# other regressors not partialled out).
partial_fe <- function(x) {
  m <- feols(x ~ 1 | year + GID_1[t, t2], data.frame(x = x, year = d$year, GID_1 = d$GID_1, t = d$t, t2 = d$t2), notes = FALSE)
  out <- rep(NA_real_, length(x))
  out[obs(m)] <- resid(m)
  out
}
# The model's view on the estimation rows: the shock of each option, or one of
# its arms (fun: identity, or the warm arm pmax(x, 0)), then partialled out.
series_model <- function(v, fun = identity) {
  i <- d$GID_0 %in% FIG_COUNTRIES
  do.call(rbind, lapply(c("a", "b", "c"), function(o) {
    data.frame(country = d$GID_0[i], year = d$year[i], option = OPTIONS[[o]], x = partial_fe(fun(shock_of(o, v, "quadratic")))[i])
  }))
}
# On a complete grid of years, so that the lines break at the rows the model
# does not have.
complete_years <- function(s, group) {
  grid <- unique(s[c("country", group)])
  grid <- grid[rep(seq_len(nrow(grid)), each = diff(horizon) + 1), ]
  grid$year <- rep(horizon[1]:horizon[2], length.out = nrow(grid))
  merge(grid, s, all.x = TRUE)
}
series_figure <- function(v, label, file, mode = c("raw", "model", "warm")) {
  mode <- match.arg(mode)
  s <- switch(mode, raw = series(v), model = series_model(v), warm = series_model(v, function(x) pmax(x, 0)))
  if (mode != "raw") s <- complete_years(s, "option")
  s$country <- factor(cname(s$country), cname(FIG_COUNTRIES))
  subtitle <- switch(mode,
    raw = paste("(a) deviation from the 20-year trailing mean, less the country's quadratic trend; (b) the same deviation",
                "not detrended;\n(c) innovation of the Nath-Ramey-Klenow AR(3) filter (lags x country mean, country effects)",
                "over its SD. Positive: warm (wet) arm."),
    model = paste("As the model sees them: each shock on the estimation rows less the year effects and the country quadratic trends",
                  "(fitted on all countries).\nThe line is the model's warm arm less its cold arm. (a) and (b) coincide here:",
                  "they differ only in how this splits between the arms (next figure).\n(c) Nath-Ramey-Klenow AR(3) innovation",
                  "over its SD. Gaps: rows not in the model (dropped outliers, missing lags)."),
    warm = paste("The warm arm max(x, 0) as the model sees it: on the estimation rows less the year effects and the country",
                 "quadratic trends (all countries).\n(a) splits the deviation already detrended; (b) splits it first and the",
                 "regression then detrends the arm; (c) the Nath-Ramey-Klenow innovation's\npositive part. Gaps: rows not in the model."))
  title <- switch(mode, raw = paste(label, "shocks under three definitions"),
                  model = paste(label, "shocks as the model sees them (less year effects and country quadratic trends)"),
                  warm = paste("Warm arm as the model sees it (less year effects and country quadratic trends)"))
  g <- ggplot(s, aes(year, x, colour = option, linetype = option)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    facet_wrap(~country, ncol = 2) +
    scale_colour_manual(values = OPTION_COLOURS[1:3], name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", if (mode == "model") "22" else "solid", "solid"), OPTIONS[1:3]), name = NULL) +
    labs(x = NULL, y = paste(label, if (mode == "warm") "warm arm (SD)" else "shock (SD)"), title = title, subtitle = subtitle) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey30"))
  ggsave(file.path(out_dir, file), g, width = 7.5, height = 8.5, dpi = 200, bg = "white")
}
series_figure("TM", "Temperature", "shocks_temperature.png")
series_figure("RR", "Precipitation", "shocks_precipitation.png")
series_figure("TM", "Temperature", "shocks_temperature_model.png", "model")
series_figure("RR", "Precipitation", "shocks_precipitation_model.png", "model")
series_figure("TM", "Temperature", "warm_arm_model.png", "warm")

# 2. Growth and its filtered innovation, (d): as is (all the panel's years),
# and as the model sees them (the outcome of the default model and of (d) on
# the estimation rows, less the year effects and the country quadratic trends).
growth_figure <- function(file, model = FALSE) {
  if (model) {
    i <- d$GID_0 %in% FIG_COUNTRIES
    labs_s <- c("Growth (outcome of (a)-(c))", "Growth innovation (outcome of (d))")
    gs <- rbind(data.frame(country = d$GID_0[i], year = d$year[i], series = labs_s[1], x = 100 * partial_fe(d$g)[i]),
                data.frame(country = d$GID_0[i], year = d$year[i], series = labs_s[2], x = 100 * partial_fe(d$g_innov)[i]))
    gs <- complete_years(gs, "series")
  } else {
    gp <- panel[panel$GID_0 %in% FIG_COUNTRIES, ]
    labs_s <- c("Growth less the country mean", "Growth innovation (NRK filter)")
    gs <- rbind(data.frame(country = gp$GID_0, year = gp$year, series = labs_s[1], x = 100 * ifelse(gp$outlier_5, NA, gp$g_all - gp$gbar)),
                data.frame(country = gp$GID_0, year = gp$year, series = labs_s[2], x = 100 * ifelse(gp$outlier_5, NA, gp$g_innov)))
  }
  gs$series <- factor(gs$series, labs_s)
  gs$country <- factor(cname(gs$country), cname(FIG_COUNTRIES))
  g <- ggplot(gs, aes(year, x, colour = series)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    facet_wrap(~country, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = setNames(c("grey45", OPTION_COLOURS[[4]]), labs_s), name = NULL) +
    labs(x = NULL, y = "pp",
         title = if (model) "Growth and its innovation as the model sees them (less year effects and country quadratic trends)" else
           "Growth and its innovation under the Nath-Ramey-Klenow filter, option (d)",
         subtitle = paste("Innovation: residual of growth on its 3 lags and the lags x the country's mean growth, with country effects.",
                          if (model) "\nOn the estimation rows, less the year effects and the country quadratic trends (fitted on all countries). Gaps: rows not in the model." else
                            "Dropped outlier years blank.")) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey30"))
  ggsave(file.path(out_dir, file), g, width = 7.5, height = 8.5, dpi = 200, bg = "white")
}
growth_figure("growth_filter.png")
growth_figure("growth_filter_model.png", model = TRUE)

# 3. The original series, not demeaned or detrended: annual T and P (with the
# 20-year trailing mean the default deviation is taken from) and growth (all
# the panel's years; the dropped outlier years marked).
original_figure <- function(v, label, unit, file) {
  p <- panel[panel$GID_0 %in% FIG_COUNTRIES, ]
  if (v == "growth") {
    s <- data.frame(country = p$GID_0, year = p$year, series = "Growth", x = 100 * p$g_all, outlier = p$outlier_5)
  } else {
    r <- raw[raw$GID_0 %in% FIG_COUNTRIES & raw$year >= horizon[1] & raw$year <= horizon[2], ]
    s <- rbind(data.frame(country = r$GID_0, year = r$year, series = label, x = r[[v]], outlier = FALSE),
               data.frame(country = p$GID_0, year = p$year, series = "20-year trailing mean", x = p[[paste0("mean_", v, "_lag_20")]],
                          outlier = FALSE))
  }
  s <- complete_years(s, "series")
  s$outlier[is.na(s$outlier)] <- FALSE
  s$country <- factor(cname(s$country), cname(FIG_COUNTRIES))
  cols <- setNames(c("grey20", "#0072B2"), c(if (v == "growth") "Growth" else label, "20-year trailing mean"))
  g <- ggplot(s, aes(year, x, colour = series)) +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    facet_wrap(~country, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = NULL, y = unit, title = paste(label, "as observed (not demeaned or detrended)"),
         subtitle = if (v == "growth") "Growth of real GDP per capita, PWT 11.0. Red points: non-weather outlier years dropped from the estimation." else
           "ERA5, 2000 population-density weights. Blue: the 20-year trailing mean, from which the default deviation is taken.") +
    theme_minimal(base_size = 8) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey30"))
  if (v == "growth") {
    g <- g + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
      geom_point(data = s[s$outlier, ], colour = "#D55E00", size = 1.2, na.rm = TRUE)
  }
  ggsave(file.path(out_dir, file), g, width = 7.5, height = 8.5, dpi = 200, bg = "white")
}
original_figure("TM", "Temperature", "°C", "original_temperature.png")
original_figure("RR", "Precipitation", "m per year", "original_precipitation.png")
original_figure("growth", "Growth", "%", "original_growth.png")

# 3. The coefficients of the options and trends, with 95% CIs.
key_terms <- c(pTM = "Warm arm", `pTM:y` = "Warm × income", `pTM:clim_TM` = "Warm × climate", nTM = "Cold arm",
               `nTM:y` = "Cold × income", pRR = "Wet arm", nRR = "Dry arm", `nRR:agr` = "Dry × share", `nRR:clim_RR` = "Dry × climate")
cf <- coefs[coefs$term %in% names(key_terms), ]
cf$term <- factor(key_terms[cf$term], key_terms)
cf$option <- factor(OPTIONS[cf$option], rev(OPTIONS))
cf$trend <- factor(cf$trend, names(TRENDS), c("Quadratic trends", "Linear trends", "No trends"))
g <- ggplot(cf, aes(estimate, option, colour = option, shape = trend)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), width = 0, linewidth = 0.4,
                orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
  facet_wrap(~term, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = OPTION_COLOURS, guide = "none") +
  scale_shape_manual(values = c(16, 17, 1), name = NULL) +
  labs(x = "pp of growth per SD (interactions: per log point of income, per °C, per unit of share)", y = NULL,
       title = "Linear interaction model (M2b) by shock definition and trends, PWT × ERA5",
       subtitle = sprintf("Same %s rows for all; 95%% CI clustered by country.", format(nrow(d), big.mark = ","))) +
  theme_minimal(base_size = 8) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey30"))
ggsave(file.path(out_dir, "coefficients.png"), g, width = 8, height = 6.5, dpi = 200, bg = "white")

options(width = 200)
cat("\nFilter coefficients (Eq. 4 analogue)\n"); print(lapply(filter_coefs, round, 4))
cat("\nPersistence (sum of lag coefficients) at the 10th/50th/90th percentile of the state\n"); print(round(persistence, 3))
cat("\nSD of the innovation / SD of the default deviation (median over countries)\n"); print(round(sd_ratio, 2))
cat("\nWarm arm by decade\n"); print(by_decade, digits = 3, row.names = FALSE)
cat("\nKey coefficients (pp; SE; p)\n")
w <- reshape(transform(cf[c("option", "trend", "term")], v = sprintf("%+.3f (%.3f)%s", cf$estimate, cf$se,
             ifelse(cf$p < 0.01, "***", ifelse(cf$p < 0.05, "**", ifelse(cf$p < 0.1, "*", ""))))),
             idvar = c("option", "trend"), timevar = "term", direction = "wide")
names(w) <- sub("^v\\.", "", names(w))
print(w[order(w$option, w$trend), ], row.names = FALSE)
cat("\nJoint tests (p)\n"); print(cbind(tests[1:2], round(tests[-(1:2)], 3)), row.names = FALSE)
