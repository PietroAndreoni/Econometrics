# Test of the BHM (Burke, Hsiang & Miguel 2015) micro-foundation: the global
# quadratic in annual mean temperature as the aggregate of a daily response,
# because exposure to hot and cold days is non-linear in the annual mean.
#
# Data: ERA5 daily mean / max / min temperature, population weighted (2015
# weights, the only daily ones), countries (Weighted Climate Dataset,
# gadm0_era_*_pop_2015_daily). GDP: the response explorer's PWT x ERA5
# (population) panel, 1970+, reversal outliers (r5), its controls, year effects
# and country quadratic trends, errors clustered by country. The annual mean T
# of every model is the mean of the daily file (not the panel's monthly
# pop-2000 T), so that T and the day counts come from the same series.
#
# A. Exposure (climate only)
#   A1 cross-country: climatological exposure against climatological T.
#   A2 within country: exposure on d = T - Tbar_i, d x Tbar_i and d^2 (country
#      FE). A global quadratic E = a + bT + cT^2 means coef(d x Tbar) = 2c
#      and coef(d^2) = c; the ratio says how much of the global curvature is
#      local.
# B. GDP growth
#   B1 BHM: g ~ T + T^2 + P + P^2.
#   B2 the same split exactly: g ~ d + d x Tbar + d^2 (+ P, P^2). BHM imposes
#      coef(d^2) = coef(d x Tbar) / 2; tested, and coef(d^2) = 0 (local
#      linearity).
#   B3 daily bins (5 C, ref 10-15): sufficiency of the bins for T, T^2 and
#      vice versa.
#   B4 the annual-T slope each country gets from the bins model,
#      sum_b beta_b dN_b/dT (within-country exposure slopes), and its
#      curvature, against BHM's beta1 + 2 beta2 Tbar_i.
#
#   Rscript bhm_microfoundation.R

suppressPackageStartupMessages({ library(fixest); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Daily data -> country-year exposure ------------------------------------
read_daily <- function(v) {
  path <- wcd_download(v, "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)
  x <- as.data.frame(arrow::read_parquet(path))
  dates <- as.Date(sub("^X", "", x$Date), "%Y%m%d")
  x$Date <- NULL
  x <- x[, !is.na(names(x)) & !(names(x) %in% c("NA", "?", "")), drop = FALSE]
  m <- as.matrix(x); storage.mode(m) <- "double"
  list(year = as.integer(format(dates, "%Y")), m = m)
}
tm <- read_daily("tmp"); tx <- read_daily("tmpmax"); tn <- read_daily("tmpmin")
stopifnot(identical(tm$year, tx$year), identical(tm$year, tn$year))
yr <- tm$year
# Sums by year; a year with any missing day is missing.
by_year <- function(m) {
  s <- rowsum(m, yr, na.rm = FALSE)
  s
}
to_long <- function(s, name) {
  data.frame(year = as.integer(rownames(s))[row(s)], GID_0 = colnames(s)[col(s)],
             value = as.vector(s), stringsAsFactors = FALSE) |>
    setNames(c("year", "GID_0", name))
}
ndays <- by_year((!is.na(tm$m)) * 1)
BREAKS <- c(-Inf, 0, 5, 10, 15, 20, 25, 30, Inf)
BIN_NAMES <- c("b_lt0", "b0_5", "b5_10", "b10_15", "b15_20", "b20_25", "b25_30", "b_gt30")
REF_BIN <- "b10_15"
ex <- to_long(by_year(tm$m) / ndays, "T")
for (j in seq_along(BIN_NAMES)) {
  ind <- (tm$m > BREAKS[j] & tm$m <= BREAKS[j + 1]) * 1
  ex[[BIN_NAMES[j]]] <- as.vector(by_year(ind))
}
measures <- list(
  dd_gt25  = pmax(tm$m - 25, 0),    # hot degree days, daily mean above 25 C
  dd_gt20  = pmax(tm$m - 20, 0),
  dd_lt5   = pmax(5 - tm$m, 0),     # cold degree days, daily mean below 5 C
  dd_lt0   = pmax(0 - tm$m, 0),
  tx_gt35  = (tx$m > 35) * 1,       # days with max above 35 C
  ddx_gt30 = pmax(tx$m - 30, 0),    # degree days of the max above 30 C
  tn_lt0   = (tn$m < 0) * 1         # frost days
)
for (k in names(measures)) ex[[k]] <- as.vector(by_year(measures[[k]]))
ex$ndays <- as.vector(ndays)
ex <- ex[is.finite(ex$T) & ex$ndays >= 365, ]
EXPOSURES <- c(BIN_NAMES, names(measures))

# --- The GDP panel (the explorer's PWT x ERA5 population rows) --------------
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
OUTLIERS <- Sys.getenv("BHM_OUTLIERS", "r5")
d <- model_data(list(gdp = "PWT", climate = "ERA5", weight = "population", units = "absolute",
                     baseline = "20-year", form = "Linear arms", shape = "linear", fe = "quadratic",
                     moderators = "trailing", outcome = "log", outliers = OUTLIERS,
                     inter = list(TM = character(0), RR = character(0))), P)
d$g <- d[[P$outcome]]
n_panel <- nrow(d)
d <- merge(d, ex, by = c("GID_0", "year"))
cat(sprintf("Panel rows %d, with daily ERA5 %d (%d countries, %d-%d)\n",
            n_panel, nrow(d), length(unique(d$GID_0)), min(d$year), max(d$year)))
cat(sprintf("corr(daily-file T, panel T) = %.4f; mean diff %.3f C\n",
            cor(d$T, d$TM), mean(d$T - d$TM)))
d$P1 <- d$RR / 1000; d$P2 <- d$P1^2         # precipitation in m, as BHM
d$T2 <- d$T^2
T_MED <- median(tapply(d$T, d$GID_0, mean))
d$Tbar <- ave(d$T, d$GID_0) - T_MED         # country mean T on the rows, centred
d$dT <- d$T - ave(d$T, d$GID_0)
d$dT2 <- d$dT^2
d$dT_Tbar <- d$dT * d$Tbar
CONTROLS <- P$controls
FE <- "year + GID_0[t, t2]"
fit <- function(rhs, data = d) {
  f <- as.formula(paste("g ~", paste(c(rhs, "P1", "P2", CONTROLS), collapse = " + "), "|", FE))
  feols(f, data = data, cluster = ~GID_0, panel.id = ~GID_0 + year, notes = FALSE)
}
G <- length(unique(d$GID_0))
# Clustered Wald test of R b = r.
wald <- function(m, R, r = rep(0, nrow(R))) {
  b <- coef(m)[colnames(R)]; V <- vcov(m)[colnames(R), colnames(R), drop = FALSE]
  Rb <- R %*% b - r
  F <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R), Rb)) / nrow(R)
  c(F = F, df = nrow(R), p = pf(F, nrow(R), G - 1, lower.tail = FALSE))
}
zero <- function(m, terms) {
  R <- diag(length(terms)); dimnames(R) <- list(terms, terms); wald(m, R)
}

# --- A. Exposure mapping -----------------------------------------------------
countries <- unique(d$GID_0)
e <- ex[ex$GID_0 %in% countries & ex$year >= 1970 & ex$year <= max(d$year), ]
clim <- aggregate(e[c("T", EXPOSURES)], list(GID_0 = e$GID_0), mean)
e$Tbar <- ave(e$T, e$GID_0) - T_MED
e$dT <- e$T - ave(e$T, e$GID_0); e$dT2 <- e$dT^2
sd_T <- tapply(e$dT, e$GID_0, sd)
A <- do.call(rbind, lapply(EXPOSURES, function(k) {
  # Cross-country: the climatology's curvature (quadratic in climate T).
  cs <- lm(clim[[k]] ~ clim$T + I(clim$T^2))
  # Within country, with the global-quadratic split.
  m <- feols(as.formula(paste(k, "~ dT + dT:Tbar + dT2 | GID_0")), e, cluster = ~GID_0, notes = FALSE)
  b <- coef(m); se <- se(m)
  # Country by country: share of the within variance the linear term leaves to d^2.
  r2 <- sapply(split(e, e$GID_0), function(z) {
    if (var(z[[k]]) == 0) return(c(lin = NA, quad = NA))
    c(lin = summary(lm(z[[k]] ~ z$dT))$r.squared, quad = summary(lm(z[[k]] ~ z$dT + z$dT2))$r.squared)
  })
  data.frame(exposure = k,
             cs_curv = unname(coef(cs)[3]), cs_curv_se = sqrt(vcov(cs)[3, 3]),
             slope = b[["dT"]], slope_x_clim = b[["dT:Tbar"]], slope_x_clim_se = se[["dT:Tbar"]],
             local_curv = b[["dT2"]], local_curv_se = se[["dT2"]],
             local_over_global = b[["dT2"]] / (b[["dT:Tbar"]] / 2),
             median_r2_lin = median(r2["lin", ], na.rm = TRUE),
             median_r2_gain_quad = median(r2["quad", ] - r2["lin", ], na.rm = TRUE))
}))
cat("\n== A. Exposure vs annual mean T: global (cross-country / slope x climate) vs local curvature ==\n")
print(within(A, { for (z in names(A)[-1]) assign(z, signif(get(z), 3)) }), row.names = FALSE)
cat(sprintf("Median within-country SD of annual T: %.2f C (IQR %.2f-%.2f)\n",
            median(sd_T), quantile(sd_T, .25), quantile(sd_T, .75)))

# Country-specific exposure slopes and curvatures (OLS of each exposure on d, d^2).
S <- C <- matrix(NA_real_, length(countries), length(EXPOSURES), dimnames = list(countries, EXPOSURES))
for (i in countries) {
  z <- e[e$GID_0 == i, ]
  for (k in EXPOSURES) {
    if (var(z[[k]]) == 0) { S[i, k] <- 0; C[i, k] <- 0; next }
    S[i, k] <- coef(lm(z[[k]] ~ z$dT))[2]
    C[i, k] <- coef(lm(z[[k]] ~ z$dT + z$dT2))[3]
  }
}
Tclim <- tapply(e$T, e$GID_0, mean)[countries]

# Pure mean shift of each country's daily 1970+ series (the BHM thought
# experiment): days per year in each bin at shift delta.
shift_bins <- function(delta) {
  sel <- yr >= 1970 & yr <= max(d$year)
  m <- tm$m[sel, countries, drop = FALSE] + delta
  ny <- length(unique(yr[sel]))
  sapply(seq_along(BIN_NAMES), function(j) colSums(m > BREAKS[j] & m <= BREAKS[j + 1], na.rm = TRUE) / ny) |>
    `colnames<-`(BIN_NAMES)
}
DELTAS <- seq(-2, 2, by = 0.25)
shifted <- lapply(DELTAS, shift_bins)

# --- B. GDP -----------------------------------------------------------------
BINS <- setdiff(BIN_NAMES, REF_BIN)
m_bhm  <- fit(c("T", "T2"))
m_dec  <- fit(c("dT", "dT_Tbar", "dT2"))
m_loc  <- fit(c("dT", "dT_Tbar"))
m_bins <- fit(BINS)
m_both <- fit(c(BINS, "T", "T2"))
m_dd   <- fit(c("dd_gt25", "dd_lt5"))
m_ddx  <- fit(c("ddx_gt30", "tn_lt0"))
m_dd_both <- fit(c("dd_gt25", "dd_lt5", "T", "T2"))
cat("\n== B1-B3 ==\n")
etable(m_bhm, m_dec, m_loc, m_bins, m_both, m_dd, m_ddx, m_dd_both,
       keep = c("^T$", "^T2$", "dT", "^b", "dd", "tn_", "ddx"), digits = 3, fitstat = ~ n + wr2 + aic) |> print()
b <- coef(m_bhm)
cat(sprintf("BHM peak: %.1f C\n", -b[["T"]] / (2 * b[["T2"]])))
# B2 exact split check: with Tbar the row mean, BHM = dec with a2 = 2 a3.
R <- matrix(c(0, 1, -2), 1, dimnames = list("a2-2a3", c("dT", "dT_Tbar", "dT2")))
tests <- rbind(
  "B2 BHM restriction coef(d^2) = coef(d x Tbar)/2" = wald(m_dec, R),
  "B2 local linearity coef(d^2) = 0"                 = zero(m_dec, "dT2"),
  "B2 slope varies with climate coef(d x Tbar) = 0"  = zero(m_dec, "dT_Tbar"),
  "B3 T, T^2 = 0 given bins"                         = zero(m_both, c("T", "T2")),
  "B3 bins = 0 given T, T^2"                         = zero(m_both, BINS),
  "B3 bins = 0 alone"                                = zero(m_bins, BINS),
  "B3 T, T^2 = 0 given degree days (25/5)"           = zero(m_dd_both, c("T", "T2")),
  "B3 degree days = 0 given T, T^2"                  = zero(m_dd_both, c("dd_gt25", "dd_lt5"))
)
print(round(tests, 4))

# B4 implied annual-T slope by country from the bins model: ME_i = S_i beta.
implied <- function(m, Smat, terms) {
  bb <- coef(m)[terms]; V <- vcov(m)[terms, terms]
  X <- Smat[, terms, drop = FALSE]
  list(est = as.vector(X %*% bb), se = sqrt(rowSums((X %*% V) * X)), X = X, V = V, b = bb)
}
me_bins <- implied(m_bins, S, BINS)              # within-country regression slopes
cv_bins <- implied(m_bins, 2 * C, BINS)          # second derivative within country
# Mean-shift slope and curvature (finite differences at 0 over +-1 C).
k0 <- which(DELTAS == 0); kp <- which(DELTAS == 1); km <- which(DELTAS == -1)
S_shift <- (shifted[[kp]] - shifted[[km]]) / 2
C_shift <- (shifted[[kp]] - 2 * shifted[[k0]] + shifted[[km]])
me_shift <- implied(m_bins, S_shift, BINS)
cv_shift <- implied(m_bins, C_shift, BINS)
bhm_me <- b[["T"]] + 2 * b[["T2"]] * Tclim
b_dec <- coef(m_dec)
dec_me <- b_dec[["dT"]] + b_dec[["dT_Tbar"]] * (Tclim - T_MED)

# Is the bins-implied slope linear in climate (= a global quadratic)? Project
# ME on [1, Tbar, Tbar^2]; the coefficients are linear in beta (A beta).
proj_test <- function(im, Tc) {
  Z <- cbind(1, Tc - T_MED, (Tc - T_MED)^2)
  Amat <- solve(crossprod(Z), t(Z)) %*% im$X
  g <- as.vector(Amat %*% im$b); Vg <- Amat %*% im$V %*% t(Amat)
  data.frame(term = c("intercept", "slope on climate", "curvature in climate"),
             est = g, se = sqrt(diag(Vg)), p = 2 * pt(-abs(g / sqrt(diag(Vg))), G - 1))
}
cat("\n== B4 bins-implied annual-T slope across countries, projected on climate ==\n")
cat("(within-country exposure slopes)\n"); print(proj_test(me_bins, Tclim), digits = 3)
cat("(pure mean shift)\n");                print(proj_test(me_shift, Tclim), digits = 3)
cat(sprintf("BHM: slope %.4f + %.5f x climate (2 beta2); DEC: %.4f + %.5f x climate\n",
            b[["T"]] + 2 * b[["T2"]] * T_MED, 2 * b[["T2"]], b_dec[["dT"]], b_dec[["dT_Tbar"]]))
cat(sprintf("Local curvature (2nd derivative, per C^2): BHM 2b2 = %.5f; DEC 2a3 = %.5f (se %.5f)\n",
            2 * b[["T2"]], 2 * b_dec[["dT2"]], 2 * se(m_dec)[["dT2"]]))
cat(sprintf("bins-implied local curvature: median %.5f (IQR %.5f..%.5f) [regression]; median %.5f [shift]\n",
            median(cv_bins$est), quantile(cv_bins$est, .25), quantile(cv_bins$est, .75), median(cv_shift$est)))
# Relative size of local curvature: slope change over +-1 SD of the country's
# annual T against the slope.
loc_rel <- abs(cv_shift$est * 2 * sd_T[countries]) / pmax(abs(me_shift$est), 1e-6)
cat(sprintf("Slope change over the country's +-1 SD from the shift curvature: median %.4f per C (median |slope| %.4f)\n",
            median(abs(cv_shift$est * 2 * sd_T[countries])), median(abs(me_shift$est))))

country_tab <- data.frame(GID_0 = countries, Tclim = as.numeric(Tclim), sd_T = as.numeric(sd_T[countries]),
                          me_bins = me_bins$est, me_bins_se = me_bins$se,
                          me_shift = me_shift$est, me_shift_se = me_shift$se,
                          cv_bins = cv_bins$est, cv_shift = cv_shift$est, cv_shift_se = cv_shift$se,
                          bhm_me = as.numeric(bhm_me), dec_me = as.numeric(dec_me))

# Bins-implied level curve across climates (the BHM aggregation done with the
# data's daily distributions): sum_b beta_b Nbar_b,i, relative to the median
# country; and each country's local response to shifts.
bb <- coef(m_bins)[BINS]
clim$lev_bins <- as.vector(as.matrix(clim[BINS]) %*% bb)
clim$lev_bhm <- b[["T"]] * clim$T + b[["T2"]] * clim$T^2
clim[c("lev_bins", "lev_bhm")] <- lapply(clim[c("lev_bins", "lev_bhm")], function(v) v - approx(clim$T, v, T_MED, ties = mean)$y)
local_curves <- do.call(rbind, lapply(seq_along(DELTAS), function(j) {
  data.frame(GID_0 = countries, Tclim = as.numeric(Tclim), delta = DELTAS[j],
             resp = as.vector((shifted[[j]][, BINS] - shifted[[k0]][, BINS]) %*% bb))
}))

saveRDS(list(A = A, tests = tests, models = list(bhm = m_bhm, dec = m_dec, loc = m_loc, bins = m_bins,
                                                 both = m_both, dd = m_dd, ddx = m_ddx, dd_both = m_dd_both),
             country = country_tab, clim = clim, local_curves = local_curves, S = S, C = C,
             proj = list(bins = proj_test(me_bins, Tclim), shift = proj_test(me_shift, Tclim)),
             T_MED = T_MED, outliers = OUTLIERS),
        file.path(out_dir, "results.rds"))

# --- Figures ------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
# A1: exposure climatology against climate.
long <- do.call(rbind, lapply(c("dd_gt25", "tx_gt35", "dd_lt5", "tn_lt0", "b25_30", "b_gt30"), function(k)
  data.frame(T = clim$T, value = clim[[k]], exposure = k)))
p1 <- ggplot(long, aes(T, value)) + geom_point(size = 1, alpha = .6) +
  geom_smooth(method = "loess", se = FALSE, formula = y ~ x, colour = "firebrick") +
  facet_wrap(~exposure, scales = "free_y") +
  labs(x = "Climatological annual mean T (C), 1970-2023", y = "days or degree days per year",
       title = "Cross-country: exposure is non-linear in climate")
ggsave(file.path(out_dir, "A1_exposure_climatology.png"), p1, width = 9, height = 5.5, dpi = 150)
# A2: within-country exposure slopes against climate.
sl <- do.call(rbind, lapply(c("dd_gt25", "tx_gt35", "dd_lt5", "tn_lt0"), function(k)
  data.frame(T = as.numeric(Tclim), slope = S[, k], curv = 2 * C[, k], exposure = k)))
p2 <- ggplot(sl, aes(T, slope)) + geom_point(size = 1, alpha = .6) +
  geom_smooth(method = "loess", se = FALSE, formula = y ~ x, colour = "firebrick") +
  facet_wrap(~exposure, scales = "free_y") +
  labs(x = "Climate (C)", y = "within-country d exposure / d annual T",
       title = "Within-country exposure slopes vary with climate")
ggsave(file.path(out_dir, "A2_exposure_slopes.png"), p2, width = 8, height = 5.5, dpi = 150)
# B4: implied annual-T slope by country.
p3 <- ggplot(country_tab, aes(Tclim)) +
  geom_errorbar(aes(ymin = me_shift - 1.96 * me_shift_se, ymax = me_shift + 1.96 * me_shift_se), colour = "grey75", width = 0) +
  geom_point(aes(y = me_shift), size = 1.2) +
  geom_line(aes(y = bhm_me, colour = "BHM quadratic: b1 + 2 b2 T")) +
  geom_line(aes(y = dec_me, colour = "split: a1 + a2 T (d^2 free)")) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_colour_manual(values = c("firebrick", "steelblue"), name = NULL) +
  labs(x = "Climate (C)", y = "d growth / d annual T (per C)",
       title = "Annual-T slope implied by the daily-bin model (mean shift), by country") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "B4_implied_slopes.png"), p3, width = 8, height = 5, dpi = 150)
# B4: level curves and local shift responses.
p4 <- ggplot() +
  geom_line(data = local_curves, aes(Tclim + delta, resp + approx(clim$T, clim$lev_bins, Tclim, ties = mean)$y, group = GID_0),
            colour = "grey60", linewidth = .3) +
  geom_point(data = clim, aes(T, lev_bins), size = 1) +
  geom_line(data = clim[order(clim$T), ], aes(T, lev_bhm), colour = "firebrick", linewidth = 1) +
  labs(x = "Annual mean T (C)", y = "growth relative to the median climate",
       title = "Bins-implied response: across countries (points) and local +-2 C shifts (grey); BHM (red)")
ggsave(file.path(out_dir, "B4_level_curves.png"), p4, width = 8, height = 5, dpi = 150)
cat("\nSaved to", out_dir, "\n")
