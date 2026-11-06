# Exposure-mediated temperature response: the physical exposure of each
# country is data driven, the economic response to exposure is pooled.
#
#   g_it = b_h (s_i^h dT_it) + b_c (s_i^c dT_it) [+ income interactions] + P terms + controls + FE
#
# dT = annual T - its mean over t-20..t-1 (ERA5 daily file, population weights
# 2015). s_i^h, s_i^c = within-country OLS slope of the trailing deviation of
# hot degree days (daily mean > 25 C) / cold degree days (< 5 C) on dT,
# 1970-2023 (climate data only, no GDP), in DD per C. So s_i dT is the
# within-country projection of the degree-day deviation on dT.
#   "country" s_i: each country's own slope
#   "smooth"  s_i: loess of the country slopes on the 1970-2023 climatology
#             (the climate interaction that reproduces the hockey stick)
# Compared with dT alone, dT + dT x climate (linear interaction, as the
# preferred spec), and the actual degree-day deviations (daily information).
#
# Panel: the explorer's PWT x ERA5 (population) rows, reversal outliers (r5),
# year effects + country quadratic trends, the panel's controls, the levels of
# the trailing climate (T, P) and trailing income, and the same P terms in
# every model (P deviation from its trailing mean over the country SD, and x
# trailing P climate). Errors clustered by country. All models on the same rows.
#
#   Rscript bhm_exposure_mediated.R

suppressPackageStartupMessages({ library(fixest); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("apps", "exposure_explorer", "R", "exposure.R"))
out_dir <- file.path("results", "bhm_microfoundation")
WINDOW <- 20
# EXPOSURE=mean (default): degree days of the daily mean above 25 C / below 5 C.
# EXPOSURE=maxmin: degree days of the daily max above 29 C / of the daily min
# below 0 C (results with the suffix _maxmin).
EXPOSURE <- Sys.getenv("EXPOSURE", "mean")
stopifnot(EXPOSURE %in% c("mean", "maxmin"))
SUFFIX <- if (EXPOSURE == "mean") "" else "_maxmin"
cat("Exposure:", if (EXPOSURE == "mean") "daily mean DD > 25 C / < 5 C" else "Tmax DD > 29 C / Tmin DD < 0 C", "\n")

# --- Exposure (exposure explorer histograms, gadm0) --------------------------
H <- readRDS(file.path("results", "exposure_explorer", "gadm0_tmp.rds"))
if (EXPOSURE == "mean") {
  hot <- index_matrix(H, list(side = "above", measure = "dd", type = "abs", value = 25))
  cold <- index_matrix(H, list(side = "below", measure = "dd", type = "abs", value = 5))
} else {
  source(file.path("functions", "climate_download.R"))
  # Units x years, aligned with H (a unit-year with a missing day is NA).
  dd_of <- function(v, f) {
    x <- as.data.frame(arrow::read_parquet(wcd_download(v, "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)))
    yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
    s <- yr %in% H$years
    m <- as.matrix(x[s, intersect(H$units$id, names(x))]); storage.mode(m) <- "double"
    out <- matrix(NA_real_, nrow(H$units), length(H$years), dimnames = list(H$units$id, H$years))
    a <- t(rowsum(f(m), yr[s]))
    out[rownames(a), colnames(a)] <- a
    out
  }
  hot <- dd_of("tmpmax", function(m) pmax(m - 29, 0))
  cold <- dd_of("tmpmin", function(m) pmax(0 - m, 0))
}
yrs <- 1970:2023
ok <- valid_units(H, seq(min(yrs) - WINDOW, max(yrs)))
need <- match(seq(min(yrs) - WINDOW, max(yrs)), H$years)
ok <- ok[rowSums(is.na(hot[ok, need, drop = FALSE]) | is.na(cold[ok, need, drop = FALSE])) == 0]
dev <- function(M) trailing_dev(M[ok, , drop = FALSE], H$years, yrs, WINDOW)
dT_m <- dev(H$annual); dH_m <- dev(hot); dC_m <- dev(cold)
ids <- H$units$id[ok]
slope <- function(dy) vapply(seq_along(ids), function(i) {
  x <- dT_m[i, ]; y <- dy[i, ]; if (stats::sd(y) == 0) 0 else sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
}, 0)
tclim <- rowMeans(H$annual[ok, match(yrs, H$years)])
expo <- data.frame(GID_0 = ids, Tclim = tclim, s_h = slope(dH_m), s_c = slope(dC_m))
lo_h <- loess(s_h ~ Tclim, expo, span = 0.5); lo_c <- loess(s_c ~ Tclim, expo, span = 0.5)
expo$s_h_smooth <- pmax(predict(lo_h), 0); expo$s_c_smooth <- pmin(predict(lo_c), 0)
long_dev <- function(M, name) data.frame(GID_0 = rep(ids, length(yrs)), year = rep(yrs, each = length(ids)),
                                         v = as.vector(M)) |> setNames(c("GID_0", "year", name))
devs <- Reduce(function(a, b) merge(a, b, by = c("GID_0", "year")),
               list(long_dev(dT_m, "dT"), long_dev(dH_m, "dHDD"), long_dev(dC_m, "dCDD")))

# --- Panel -----------------------------------------------------------------------
app_dir <- file.path("apps", "response_explorer")
source(file.path(app_dir, "R", "response.R")); source(file.path(app_dir, "R", "fit.R"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
lines <- readLines(file.path(app_dir, "R", "fit.R"))
s0 <- grep("^explorer_model <- function", lines)
s1 <- s0 - 1 + which(lines[s0:length(lines)] == "}")[1]
fn <- paste(lines[s0:s1], collapse = "\n")
fn <- sub("explorer_model <- function", "model_data <- function", fn)
fn <- sub("if (!is.null(model_cache[[key]])) return(model_cache[[key]])", "", fn, fixed = TRUE)
fn <- sub("m <- fit(o$shape)", "return(d)", fn, fixed = TRUE)
stopifnot(grepl("return(d)", fn, fixed = TRUE))
eval(parse(text = fn))
d <- model_data(list(gdp = "PWT", climate = "ERA5", weight = "population", units = "absolute",
                     baseline = "20-year", form = "Linear arms", shape = "linear", fe = "quadratic",
                     moderators = "trailing", outcome = "log", outliers = "r5",
                     inter = list(TM = c("clim", "inc"), RR = c("clim"))), P)
d$g <- d[[P$outcome]]
n0 <- nrow(d)
d <- merge(d, devs, by = c("GID_0", "year"))
d <- merge(d, expo, by = "GID_0")
cat(sprintf("Rows %d of %d, %d countries, %d-%d\n", nrow(d), n0, length(unique(d$GID_0)), min(d$year), max(d$year)))
# P: deviation from its trailing mean over the country SD, and x trailing P climate.
d$zP <- (d$RR - d$mean_RR_lag_20) / d$sd_RR_all
d$zP_clim <- d$zP * d$clim_RR
d$y <- d$inc
d$dT_clim <- d$dT * d$clim_TM
d$dT_inc <- d$dT * d$y
# Exposure-mediated regressors, per 100 DD.
for (v in c("", "_smooth")) {
  d[[paste0("eh", v)]] <- d$dT * d[[paste0("s_h", v)]] / 100
  d[[paste0("ec", v)]] <- d$dT * d[[paste0("s_c", v)]] / 100
  d[[paste0("eh", v, "_inc")]] <- d[[paste0("eh", v)]] * d$y
  d[[paste0("ec", v, "_inc")]] <- d[[paste0("ec", v)]] * d$y
}
d$hdd <- d$dHDD / 100; d$cdd <- d$dCDD / 100
d$hdd_inc <- d$hdd * d$y; d$cdd_inc <- d$cdd * d$y

COMMON <- c("zP", "zP_clim", "clim_TM", "clim_RR", "y", P$controls)
FE <- "year + GID_0[t, t2]"
fit <- function(terms) {
  f <- as.formula(paste("g ~", paste(c(terms, COMMON), collapse = " + "), "|", FE))
  m <- feols(f, data = d, cluster = ~GID_0, panel.id = ~GID_0 + year, notes = FALSE)
  m$country <- d$GID_0[obs(m)]
  m
}
SPECS <- list(
  "dT"                      = "dT",
  "dT + dT x climate"       = c("dT", "dT_clim"),
  "exposure: country s_i"   = c("eh", "ec"),
  "exposure: smooth s(clim)" = c("eh_smooth", "ec_smooth"),
  "exposure country + dT"   = c("eh", "ec", "dT"),
  "actual DD deviations"    = c("hdd", "cdd"),
  "actual DD + dT"          = c("hdd", "cdd", "dT"),
  "dT + dT x clim + income" = c("dT", "dT_clim", "dT_inc"),
  "exposure country + income" = c("eh", "ec", "eh_inc", "ec_inc"),
  "exposure smooth + income"  = c("eh_smooth", "ec_smooth", "eh_smooth_inc", "ec_smooth_inc"),
  "actual DD + income"        = c("hdd", "cdd", "hdd_inc", "cdd_inc"),
  # Horse races: the kink (exposure) and income as separate moderators of dT.
  "dT + dT x income"                    = c("dT", "dT_inc"),
  "dT + dT x income + exposure country" = c("dT", "dT_inc", "eh", "ec"),
  "dT + dT x income + exposure smooth"  = c("dT", "dT_inc", "eh_smooth", "ec_smooth"),
  "dT + dT x income + dT x clim + exposure smooth" = c("dT", "dT_inc", "dT_clim", "eh_smooth", "ec_smooth")
)
models <- lapply(SPECS, fit)
stopifnot(length(unique(vapply(models, nobs, 0))) == 1)
etable(models[1:7], digits = 3, fitstat = ~ n + wr2 + aic + bic,
       keep = c("^dT", "^eh", "^ec", "^hdd", "^cdd"), headers = names(SPECS)[1:7]) |> print()
etable(models[8:11], digits = 3, fitstat = ~ n + wr2 + aic + bic,
       keep = c("^dT", "^eh", "^ec", "^hdd", "^cdd"), headers = names(SPECS)[8:11]) |> print()

# --- Comparisons (noise rule) -----------------------------------------------------
G <- length(unique(d$GID_0))
loglik_i <- function(m) { e <- resid(m); s2 <- mean(e^2); -0.5 * log(2 * pi * s2) - e^2 / (2 * s2) }
vuong <- function(m1, m2, penalty = c("AIC", "BIC")) {
  penalty <- match.arg(penalty)
  k <- length(coef(m1)) - length(coef(m2)); n <- nobs(m1)
  pen <- if (penalty == "AIC") k else k * log(n) / 2
  D <- tapply(loglik_i(m1) - loglik_i(m2), m1$country, sum)
  z <- (sum(D) - pen) / sqrt(sum((D - mean(D))^2) * G / (G - 1))
  z
}
wald_zero <- function(m, terms) {
  b <- coef(m)[terms]; V <- vcov(m)[terms, terms, drop = FALSE]
  F <- as.numeric(t(b) %*% solve(V, b)) / length(terms)
  pf(F, length(terms), G - 1, lower.tail = FALSE)
}
verdict <- function(z) if (abs(z) < 1.96) "indifferent" else if (z > 0) "better" else "worse"
ref <- models[["dT + dT x climate"]]
cmp <- do.call(rbind, lapply(names(models), function(k) {
  m <- models[[k]]
  zA <- vuong(m, ref, "AIC"); zB <- vuong(m, ref, "BIC")
  data.frame(model = k, k = length(SPECS[[k]]), wr2 = fitstat(m, "wr2")[[1]],
             dAIC_vs_ref = AIC(m) - AIC(ref), dBIC_vs_ref = BIC(m) - BIC(ref),
             vuong_z_AIC = zA, AIC_verdict = if (identical(m, ref)) "-" else verdict(zA),
             vuong_z_BIC = zB, BIC_verdict = if (identical(m, ref)) "-" else verdict(zB))
}))
cat("\n== vs dT + dT x climate (clustered Vuong, positive = better than the reference) ==\n")
print(cmp, digits = 3, row.names = FALSE)
nested <- rbind(
  "dT adds to country exposure"            = wald_zero(models[["exposure country + dT"]], "dT"),
  "dT adds to actual DD"                   = wald_zero(models[["actual DD + dT"]], "dT"),
  "cold arm = 0 (country exposure)"        = wald_zero(models[["exposure: country s_i"]], "ec"),
  "cold arm = 0 (smooth exposure)"         = wald_zero(models[["exposure: smooth s(clim)"]], "ec_smooth"),
  "cold arm = 0 (actual DD)"               = wald_zero(models[["actual DD deviations"]], "cdd"),
  "income terms = 0 (country exposure)"    = wald_zero(models[["exposure country + income"]], c("eh_inc", "ec_inc")),
  "income terms = 0 (actual DD)"           = wald_zero(models[["actual DD + income"]], c("hdd_inc", "cdd_inc")),
  "income term = 0 (dT x clim)"            = wald_zero(models[["dT + dT x clim + income"]], "dT_inc"),
  "climate term = 0 (dT x clim + income)"  = wald_zero(models[["dT + dT x clim + income"]], "dT_clim"),
  "exposure country = 0 | dT, dT x inc"    = wald_zero(models[["dT + dT x income + exposure country"]], c("eh", "ec")),
  "hot exposure country = 0 | dT, dT x inc" = wald_zero(models[["dT + dT x income + exposure country"]], "eh"),
  "dT x inc = 0 | dT, exposure country"    = wald_zero(models[["dT + dT x income + exposure country"]], "dT_inc"),
  "exposure smooth = 0 | dT, dT x inc"     = wald_zero(models[["dT + dT x income + exposure smooth"]], c("eh_smooth", "ec_smooth")),
  "hot exposure smooth = 0 | dT, dT x inc" = wald_zero(models[["dT + dT x income + exposure smooth"]], "eh_smooth"),
  "dT x inc = 0 | dT, exposure smooth"     = wald_zero(models[["dT + dT x income + exposure smooth"]], "dT_inc"),
  "dT x clim = 0 | dT, dT x inc, exposure smooth" = wald_zero(models[["dT + dT x income + dT x clim + exposure smooth"]], "dT_clim"),
  "exposure smooth = 0 | dT, dT x inc, dT x clim" = wald_zero(models[["dT + dT x income + dT x clim + exposure smooth"]], c("eh_smooth", "ec_smooth"))
)
cat(sprintf("corr across countries: s_h vs mean trailing income %.2f, smooth s_h vs income %.2f, climate vs income %.2f\n",
            cor(tapply(d$s_h, d$GID_0, mean), tapply(d$y, d$GID_0, mean)),
            cor(tapply(d$s_h_smooth, d$GID_0, mean), tapply(d$y, d$GID_0, mean)),
            cor(tapply(d$Tclim, d$GID_0, mean), tapply(d$y, d$GID_0, mean))))
etable(models[12:15], digits = 3, fitstat = ~ wr2 + aic + bic,
       keep = c("^dT", "^eh", "^ec"), headers = names(SPECS)[12:15]) |> print()
colnames(nested) <- "p"
cat("\n== Nested tests (clustered Wald) ==\n"); print(round(nested, 4))

# --- Implied response to +1 C by country ------------------------------------------
cc <- expo[expo$GID_0 %in% d$GID_0, ]
clim_med <- median(tapply(d$mean_TM_lag_20, d$GID_0, mean))   # rough centring of the T climate for the line
me <- function(m, X) {
  b <- coef(m)[colnames(X)]; V <- vcov(m)[colnames(X), colnames(X)]
  data.frame(est = as.vector(X %*% b), se = sqrt(rowSums((X %*% V) * X)))
}
# The climate moderator of the linear model at each country's climatology
# (clim_TM is the trailing climate less the panel median).
cen <- mean(d$mean_TM_lag_20 - d$clim_TM)
Xc <- cbind(dT = 1, dT_clim = cc$Tclim - cen)
Xe <- cbind(eh = cc$s_h / 100, ec = cc$s_c / 100)
Xs <- cbind(eh_smooth = cc$s_h_smooth / 100, ec_smooth = cc$s_c_smooth / 100)
resp <- rbind(cbind(model = "dT + dT x climate", cc[c("GID_0", "Tclim")], me(ref, Xc)),
              cbind(model = "exposure: country s_i", cc[c("GID_0", "Tclim")], me(models[["exposure: country s_i"]], Xe)),
              cbind(model = "exposure: smooth s(clim)", cc[c("GID_0", "Tclim")], me(models[["exposure: smooth s(clim)"]], Xs)))
p <- ggplot(resp, aes(Tclim, est, colour = model)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(data = resp[resp$model == "exposure: country s_i", ], size = 1.3, alpha = .7) +
  geom_line(data = resp[resp$model != "exposure: country s_i", ], linewidth = 1) +
  geom_ribbon(data = resp[resp$model != "exposure: country s_i", ],
              aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se, fill = model), alpha = .12, colour = NA) +
  scale_colour_manual(values = c("firebrick", "black", "steelblue")) +
  scale_fill_manual(values = c("firebrick", "steelblue"), guide = "none") +
  labs(x = "Climate (C), 1970-2023", y = "growth response to a +1 C deviation",
       colour = NULL, title = "Temperature response: linear climate interaction vs exposure-mediated",
       subtitle = "points: pooled response per degree day x each country's exposure slope; lines: 95% bands") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("B5_exposure_mediated_response", SUFFIX, ".png")), p, width = 8.5, height = 5.5, dpi = 150)
saveRDS(list(models = models, comparison = cmp, nested = nested, exposure = expo, response = resp),
        file.path(out_dir, paste0("exposure_mediated", SUFFIX, ".rds")))
cat("Saved", file.path(out_dir, paste0("B5_exposure_mediated_response", SUFFIX, ".png")), "\n")

# --- B6: the exposure x income model, response to +1 C by climate and income ---
# dg/dT = (b_H + g_H y) s_i^H / 100 + (b_C + g_C y) s_i^C / 100, at y = -1, 0, +1
# (log income less its mean: about 0.37x, 1x and 2.7x the geometric-mean income).
cat("\n== Exposure x income: coefficients (per 100 DD) ==\n")
for (k in c("exposure country + income", "exposure smooth + income")) {
  m <- models[[k]]; terms <- SPECS[[k]]
  ct <- coeftable(m)[terms, ]
  cat(k, "\n"); print(round(ct, 5))
}
inc_lines <- do.call(rbind, lapply(c(-1, 0, 1), function(yv) {
  do.call(rbind, lapply(list(c("exposure country + income", ""), c("exposure smooth + income", "_smooth")), function(sp) {
    m <- models[[sp[1]]]; v <- sp[2]
    X <- cbind(cc[[paste0("s_h", v)]], cc[[paste0("s_h", v)]] * yv, cc[[paste0("s_c", v)]], cc[[paste0("s_c", v)]] * yv) / 100
    colnames(X) <- c(paste0("eh", v), paste0("eh", v, "_inc"), paste0("ec", v), paste0("ec", v, "_inc"))
    cbind(model = sp[1], income = sprintf("log income %+d (x%.2f)", yv, exp(yv)), cc[c("GID_0", "Tclim")], me(m, X))
  }))
}))
p6 <- ggplot(inc_lines, aes(Tclim, est, colour = income)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(data = inc_lines[inc_lines$model == "exposure country + income", ], size = 1, alpha = .6) +
  geom_line(data = inc_lines[inc_lines$model == "exposure smooth + income", ], linewidth = 1) +
  geom_ribbon(data = inc_lines[inc_lines$model == "exposure smooth + income", ],
              aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se, fill = income), alpha = .12, colour = NA) +
  scale_colour_manual(values = c("firebrick", "black", "steelblue")) +
  scale_fill_manual(values = c("firebrick", "black", "steelblue"), guide = "none") +
  labs(x = "Climate (C), 1970-2023", y = "growth response to a +1 C deviation", colour = NULL,
       title = paste("Exposure x income:", if (EXPOSURE == "mean") "daily mean DD > 25 C / < 5 C" else "Tmax DD > 29 C / Tmin DD < 0 C"),
       subtitle = "points: country s_i; lines: smooth s(climate) with 95% bands; income relative to the geometric mean") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("B6_exposure_income_response", SUFFIX, ".png")), p6, width = 8.5, height = 5.5, dpi = 150)
cat("Saved", file.path(out_dir, paste0("B6_exposure_income_response", SUFFIX, ".png")), "\n")
