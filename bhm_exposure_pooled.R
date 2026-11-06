# The exposure x income model on pooled subnational data (DOSE regions of the
# USA, China, India and Mexico), with year effects or country x year effects.
#
#   g_rt = (b_H + g_H^N yN_ct + g_H^R yR_rt) e_r^H dT_rt / 100
#        + (b_C + g_C^N yN_ct + g_C^R yR_rt) e_r^C dT_rt / 100
#        + P terms + trailing levels + region FE + region quadratic trends + [year | country x year] FE
#
# Exposure: degree days of the daily max above 29 C and of the daily min below
# 0 C, ERA5 daily, population weights 2015, GADM1; e_r = within-region slope
# of their trailing deviations on dT (as bhm_exposure_subnational.R), "country"
# version.
# Income, with no currency conversion of DOSE:
#   yR = relative income: log mean regional GDP per capita (DOSE, constant
#        local prices) of t-5..t-1, less the country's mean over its rows
#   yN = national income: log mean PWT GDP per capita of t-5..t-1 (the
#        national panel's inc_trailing), less its mean over the national
#        rows, as in bhm_exposure_mediated.R; so b_H and g_H^N compare with
#        the national estimates
# Errors clustered by region; two-way (region + year) as a check. Rows are
# unbalanced (USA, CHN 1970-2019; IND 1985-2019; MEX 2005-2019, after the
# income lags).
#
#   Rscript bhm_exposure_pooled.R
#   OIL=keep Rscript bhm_exposure_pooled.R     # keep Campeche and Tabasco

suppressPackageStartupMessages({ library(fixest); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
COUNTRIES <- c("USA", "CHN", "IND", "MEX")
OIL <- Sys.getenv("OIL", "drop"); stopifnot(OIL %in% c("drop", "keep"))
DROP <- if (OIL == "drop") c("Campeche", "Tabasco") else character(0)
TAG <- if (OIL == "keep") "_oil" else ""
WINDOW <- 20; CLIM <- 1970:2023; YEARS <- 1940:2023; END <- 2019

# --- Weather and exposure (GADM1) -----------------------------------------------------
read_cols <- function(v, cols) {
  path <- wcd_download(v, "ERA5", "gadm1", "pop", 2015, "daily", verbose = FALSE)
  x <- as.data.frame(arrow::read_parquet(path, col_select = tidyselect::all_of(c("Date", cols))))
  yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
  s <- yr %in% YEARS
  m <- as.matrix(x[s, cols, drop = FALSE]); storage.mode(m) <- "double"
  list(m = m, yr = yr[s])
}
all_cols <- setdiff(names(arrow::open_dataset(wcd_download("tmp", "ERA5", "gadm1", "pop", 2015, "daily", verbose = FALSE))), "Date")
cols <- all_cols[substr(all_cols, 1, 3) %in% COUNTRIES]
tm <- read_cols("tmp", cols); tx <- read_cols("tmpmax", cols); tn <- read_cols("tmpmin", cols); pr <- read_cols("pre", cols)
yr <- tm$yr
bad <- colSums(is.na(tm$m) | is.na(tx$m) | is.na(tn$m) | is.na(pr$m)) > 0
cols <- cols[!bad]
agg <- function(m) rowsum(m[, !bad, drop = FALSE], yr)
Tann <- agg(tm$m) / agg((!is.na(tm$m)) * 1); Pann <- agg(pmax(pr$m, 0))
hot <- agg(pmax(tx$m - 29, 0)); cold <- agg(pmax(0 - tn$m, 0))
years <- as.integer(rownames(Tann))
gid1 <- sub("^([A-Z]{3})_", "\\1.", cols)
tdev <- function(M, yy) {
  ci <- match(yy, years); out <- M[ci, , drop = FALSE]
  for (j in seq_along(ci)) out[j, ] <- M[ci[j], ] - colMeans(M[(ci[j] - WINDOW):(ci[j] - 1), , drop = FALSE])
  out
}
trail_mean <- function(M, yy) {
  ci <- match(yy, years)
  t(vapply(ci, function(k) colMeans(M[(k - WINDOW):(k - 1), , drop = FALSE]), numeric(ncol(M))))
}
dTc <- tdev(Tann, CLIM)
slope <- function(dy) vapply(seq_along(cols), function(i) {
  x <- dTc[, i] - mean(dTc[, i]); y <- dy[, i] - mean(dy[, i])
  if (sum(y^2) == 0) 0 else sum(x * y) / sum(x^2)
}, 0)
expo <- data.frame(GID_1 = gid1, Tclim = colMeans(Tann[match(CLIM, years), ]),
                   s_h = slope(tdev(hot, CLIM)), s_c = slope(tdev(cold, CLIM)))
yy <- 1970:END
long <- function(M, name) setNames(data.frame(rep(gid1, each = length(yy)), rep(yy, length(gid1)), as.vector(M)),
                                   c("GID_1", "year", name))
wx <- Reduce(function(a, b) merge(a, b, by = c("GID_1", "year")), list(
  long(tdev(Tann, yy), "dT"), long(tdev(Pann, yy), "dP_mm"),
  long(trail_mean(Tann, yy), "climT"), long(trail_mean(Pann, yy), "climP")))
sdP <- setNames(apply(Pann[match(CLIM, years), ], 2, sd), gid1)
wx$zP <- wx$dP_mm / sdP[wx$GID_1]

# --- National income (the national panel's trailing PWT income, centred as there) ----
app_dir <- file.path("apps", "response_explorer")
source(file.path(app_dir, "R", "response.R")); source(file.path(app_dir, "R", "fit.R"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
lines <- readLines(file.path(app_dir, "R", "fit.R"))
s0 <- grep("^explorer_model <- function", lines); s1 <- s0 - 1 + which(lines[s0:length(lines)] == "}")[1]
fn <- paste(lines[s0:s1], collapse = "\n")
fn <- sub("explorer_model <- function", "model_data <- function", fn)
fn <- sub("if (!is.null(model_cache[[key]])) return(model_cache[[key]])", "", fn, fixed = TRUE)
fn <- sub("m <- fit(o$shape)", "return(d)", fn, fixed = TRUE); eval(parse(text = fn))
nat <- model_data(list(gdp = "PWT", climate = "ERA5", weight = "population", units = "absolute",
                       baseline = "20-year", form = "Linear arms", shape = "linear", fe = "quadratic",
                       moderators = "trailing", outcome = "log", outliers = "r5",
                       inter = list(TM = c("clim", "inc"), RR = c("clim"))), P)
nat_fit <- readRDS(file.path(out_dir, "exposure_mediated_maxmin.rds"))$models[["exposure country + income"]]
stopifnot(nobs(nat_fit) == nrow(nat))
# yN for every country-year with the PWT trailing income (not only the
# national estimation rows, which drop outlier years), centred at the
# national rows' mean.
rows_all <- P$panels[["PWT | ERA5 | population"]]$rows
yN <- data.frame(GID_0 = rows_all$GID_0, year = rows_all$year, yN = rows_all$inc_trailing - mean(nat$inc_trailing))
yN <- yN[yN$GID_0 %in% COUNTRIES & is.finite(yN$yN), ]

# --- DOSE ---------------------------------------------------------------------------------
dose <- utils::read.csv(file.path("data", "DOSE_V2.14.csv"), stringsAsFactors = FALSE)
dose <- dose[dose$GID_0 %in% COUNTRIES & is.finite(dose$grp_pc_lcu_2015) & dose$grp_pc_lcu_2015 > 0,
             c("GID_0", "GID_1", "region", "year", "grp_pc_lcu_2015")]
lag_by <- function(v, id, t, k) v[match(paste(id, t - k), paste(id, t))]
dose$g <- log(dose$grp_pc_lcu_2015) - log(lag_by(dose$grp_pc_lcu_2015, dose$GID_1, dose$year, 1))
dose$inc_reg <- log(rowMeans(sapply(1:5, function(k) lag_by(dose$grp_pc_lcu_2015, dose$GID_1, dose$year, k))))
d <- merge(dose, wx, by = c("GID_1", "year"))
d <- merge(d, expo, by = "GID_1")
d <- merge(d, yN, by = c("GID_0", "year"))
d <- d[d$year <= END & is.finite(d$g) & is.finite(d$inc_reg) & !(d$region %in% DROP), ]
d$yR <- d$inc_reg - ave(d$inc_reg, d$GID_0)
d$clim_TM <- d$climT - median(tapply(d$climT, d$GID_1, mean))
d$clim_RR <- (d$climP - median(tapply(d$climP, d$GID_1, mean))) / 1000
d$zP_clim <- d$zP * d$clim_RR
d$t <- d$year - 1990; d$t2 <- d$t^2
d$eh <- d$dT * d$s_h / 100; d$ec <- d$dT * d$s_c / 100
for (z in c("yN", "yR")) {
  d[[paste0("eh_", z)]] <- d$eh * d[[z]]; d[[paste0("ec_", z)]] <- d$ec * d[[z]]
  d[[paste0("dT_", z)]] <- d$dT * d[[z]]
}
d$dT_clim <- d$dT * d$clim_TM
d$cy <- paste(d$GID_0, d$year)
cat("Rows by country:\n")
print(do.call(rbind, lapply(split(d, d$GID_0), function(x) data.frame(GID_0 = x$GID_0[1], rows = nrow(x), regions = length(unique(x$GID_1)),
  years = paste(range(x$year), collapse = "-"), yN = round(mean(x$yN), 2), yR_p10_p90 = round(diff(quantile(x$yR, c(.1, .9))), 2),
  max_s_h = round(max(x$s_h)), min_s_c = round(min(x$s_c))))), row.names = FALSE)
if (length(DROP)) cat("Left out:", paste(DROP, collapse = ", "), "\n")

# --- Models -------------------------------------------------------------------------------
SPECS <- list(
  "dT + dT x climate"           = c("dT", "dT_clim"),
  "dT + dT x income (N, R)"     = c("dT", "dT_yN", "dT_yR"),
  "exposure"                    = c("eh", "ec"),
  "exposure x income (R)"       = c("eh", "ec", "eh_yR", "ec_yR"),
  "exposure x income (N, R)"    = c("eh", "ec", "eh_yN", "ec_yN", "eh_yR", "ec_yR"),
  "hot exposure x income (N, R)" = c("eh", "eh_yN", "eh_yR")
)
FES <- c("year FE" = "year + GID_1[t, t2]", "country x year FE" = "cy + GID_1[t, t2]")
COMMON <- c("zP", "zP_clim", "clim_TM", "clim_RR", "yR", "yN")
fit <- function(terms, fe, data = d) {
  f <- as.formula(paste("g ~", paste(c(terms, COMMON), collapse = " + "), "|", fe))
  m <- feols(f, data = data, cluster = ~GID_1, notes = FALSE, warn = FALSE)
  m$region <- data$GID_1[obs(m)]
  m
}
models <- lapply(FES, function(fe) lapply(SPECS, fit, fe = fe))
G <- length(unique(d$GID_1))
for (k in names(FES)) {
  cat(sprintf("\n==== %s: %d rows, %d regions ====\n", k, nobs(models[[k]][[1]]), G))
  etable(models[[k]], keep = c("^dT", "^eh", "^ec"), digits = 3, fitstat = ~ wr2 + aic) |> print()
}

# Comparisons (clustered Vuong vs dT + dT x climate, within each FE).
ll <- function(m) { e <- resid(m); s2 <- mean(e^2); -0.5 * log(2 * pi * s2) - e^2 / (2 * s2) }
vuong <- function(m1, m2) {
  k <- length(coef(m1)) - length(coef(m2))
  D <- tapply(ll(m1) - ll(m2), m1$region, sum)
  (sum(D) - k) / sqrt(sum((D - mean(D))^2) * G / (G - 1))
}
cmp <- do.call(rbind, lapply(names(FES), function(k) do.call(rbind, lapply(names(SPECS), function(s) {
  m <- models[[k]][[s]]; ref <- models[[k]][["dT + dT x climate"]]
  z <- if (s == "dT + dT x climate") NA else vuong(m, ref)
  data.frame(fe = k, model = s, wr2 = fitstat(m, "wr2")[[1]], dAIC = AIC(m) - AIC(ref), vuong_z = z,
             verdict = if (is.na(z)) "-" else if (abs(z) < 1.96) "indifferent" else if (z > 0) "better" else "worse")
}))))
cat("\n== vs dT + dT x climate (same FE) ==\n"); print(cmp, digits = 3, row.names = FALSE)

# Key model: coefficients, two-way SEs, and the national benchmark.
key_terms <- SPECS[["exposure x income (N, R)"]]
bn <- coef(nat_fit)[c("eh", "ec", "eh_inc", "ec_inc")]; sn <- se(nat_fit)[c("eh", "ec", "eh_inc", "ec_inc")]
key <- do.call(rbind, lapply(names(FES), function(k) {
  m <- models[[k]][["exposure x income (N, R)"]]; m2 <- summary(m, cluster = ~GID_1 + year)
  nat_b <- c(eh = bn[["eh"]], ec = bn[["ec"]], eh_yN = bn[["eh_inc"]], ec_yN = bn[["ec_inc"]], eh_yR = NA, ec_yR = NA)
  nat_s <- c(eh = sn[["eh"]], ec = sn[["ec"]], eh_yN = sn[["eh_inc"]], ec_yN = sn[["ec_inc"]], eh_yR = NA, ec_yR = NA)
  data.frame(fe = k, term = key_terms, est = coef(m)[key_terms], se = se(m)[key_terms], p = pvalue(m)[key_terms],
             p_2way = pvalue(m2)[key_terms], national = nat_b[key_terms],
             z_vs_national = (coef(m)[key_terms] - nat_b[key_terms]) / sqrt(se(m)[key_terms]^2 + nat_s[key_terms]^2),
             row.names = NULL)
}))
cat("\n== exposure x income (N, R): pooled regions vs the national (176 countries) estimates ==\n")
print(key, digits = 3, row.names = FALSE)

# Leave one country out (key model).
loo <- do.call(rbind, lapply(names(FES), function(k) do.call(rbind, lapply(COUNTRIES, function(c0) {
  m <- fit(key_terms, FES[[k]], d[d$GID_0 != c0, ])
  data.frame(fe = k, without = c0, t(setNames(round(coef(m)[key_terms], 4), key_terms)),
             p_eh = round(pvalue(m)[["eh"]], 3), p_eh_yN = round(pvalue(m)[["eh_yN"]], 3))
}))))
cat("\n== Leave one country out (exposure x income (N, R)) ==\n"); print(loo, row.names = FALSE)

# Implied hot damage per 100 DD at each country's national income (relative income 0).
imp <- do.call(rbind, lapply(names(FES), function(k) {
  m <- models[[k]][["exposure x income (N, R)"]]
  V <- vcov(m)[c("eh", "eh_yN"), c("eh", "eh_yN")]
  do.call(rbind, lapply(COUNTRIES, function(c0) {
    x <- c(1, mean(d$yN[d$GID_0 == c0]))
    xn <- c(1, x[2]); Vn <- vcov(nat_fit)[c("eh", "eh_inc"), c("eh", "eh_inc")]
    data.frame(fe = k, country = c0, yN = x[2], pooled = sum(x * coef(m)[c("eh", "eh_yN")]),
               pooled_se = sqrt(t(x) %*% V %*% x)[1], national = sum(xn * bn[c("eh", "eh_inc")]),
               national_se = sqrt(t(xn) %*% Vn %*% xn)[1])
  }))
}))
cat("\n== Hot damage per 100 DD at each country's mean national income ==\n"); print(imp, digits = 3, row.names = FALSE)

p <- ggplot(imp, aes(yN)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_abline(intercept = bn[["eh"]], slope = bn[["eh_inc"]], colour = "grey40", linetype = 2) +
  geom_pointrange(aes(y = national, ymin = national - 1.96 * national_se, ymax = national + 1.96 * national_se),
                  colour = "grey40", shape = 1, position = position_nudge(x = -0.06)) +
  geom_pointrange(aes(y = pooled, ymin = pooled - 1.96 * pooled_se, ymax = pooled + 1.96 * pooled_se, colour = fe),
                  position = position_dodge2(width = 0.12)) +
  geom_text(data = imp[imp$fe == "year FE", ], aes(y = pooled + 1.96 * pooled_se, label = country), vjust = -0.6, size = 3.3) +
  scale_colour_manual(values = c("firebrick", "steelblue")) +
  labs(x = "national log income (PWT, centred as in the national model)", y = "growth effect of +100 hot degree days (Tmax > 29 C)",
       colour = "pooled regions", title = "Hot damage by national income: pooled regions vs countries",
       subtitle = "grey open: national model (176 countries), dashed its line; coloured: pooled regions; 95% CI") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("C3_pooled_hot_damage", TAG, ".png")), p, width = 8.5, height = 5.5, dpi = 150)
saveRDS(list(models = models, comparison = cmp, key = key, loo = loo, implied = imp),
        file.path(out_dir, paste0("exposure_pooled", TAG, ".rds")))
cat("Saved C3_pooled_hot_damage", TAG, ".png in ", out_dir, "\n", sep = "")
