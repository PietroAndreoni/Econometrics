# The exposure x income model on subnational data, one country at a time
# (default: US states and Chinese provinces; COUNTRIES= for others): does the national shape survive?
#
#   g_rt = (b_H + g_H y_rt) s_r^H dT_rt / 100 + (b_C + g_C y_rt) s_r^C dT_rt / 100
#          + P terms + trailing levels + region FE + region quadratic trends + year FE
#
# Exposure as in bhm_exposure_mediated.R. EXPOSURE=maxmin (default): degree days of
# the daily max above 29 C (hot) and of the daily min below 0 C (cold), ERA5
# daily, population weights 2015, GADM1 (WCD). dT = annual mean T less its
# mean over t-20..t-1; s_r = within-region OLS slope of the trailing deviation
# of the degree days on dT, 1970-2023 (climate data only); smooth s_r = loess
# (span 0.75) of the region slopes on the region climatology, within the
# country.
# GDP: DOSE V2.14 regional GDP per capita, constant local prices
# (grp_pc_lcu_2015); g = growth of its log; y = log mean of t-5..t-1, less its
# mean over the rows. No outlier rule (the national one is for countries);
# TRIM=1 drops the 1% tails of growth within the country. Year effects absorb
# the national shocks; errors clustered by region (few clusters: 51 / 31), and
# two-way (region + year) for the key terms.
#
#   Rscript bhm_exposure_subnational.R            # USA and CHN, 1970-2019
#   START=1980 Rscript bhm_exposure_subnational.R # later start (END= too)
#   COUNTRIES=IND,MEX Rscript bhm_exposure_subnational.R   # any DOSE countries
#   TRENDS=linear Rscript bhm_exposure_subnational.R       # linear region trends
#   DROP=Campeche,Tabasco COUNTRIES=MEX Rscript bhm_exposure_subnational.R   # leave regions out
#   TRIM=1 Rscript bhm_exposure_subnational.R
#   EXPOSURE=mean Rscript bhm_exposure_subnational.R   # daily mean DD > 25 / < 5 C (GADM1 tmp only)

suppressPackageStartupMessages({ library(fixest); library(ggplot2) })
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
COUNTRIES <- strsplit(Sys.getenv("COUNTRIES", "USA,CHN"), ",")[[1]]
START <- as.integer(Sys.getenv("START", "1970")); END <- as.integer(Sys.getenv("END", "2019"))
TRENDS <- Sys.getenv("TRENDS", "quadratic")                 # region trends: quadratic or linear
stopifnot(TRENDS %in% c("quadratic", "linear"))
TRIM <- as.numeric(Sys.getenv("TRIM", "0")) / 100
DROP <- strsplit(Sys.getenv("DROP", ""), ",")[[1]]           # regions left out (GID_1 or DOSE name)
WINDOW <- 20; CLIM <- 1970:2023; YEARS <- 1940:2023
EXPOSURE <- Sys.getenv("EXPOSURE", "maxmin")      # maxmin: Tmax DD > 29 / Tmin DD < 0; mean: daily mean DD > 25 / < 5
stopifnot(EXPOSURE %in% c("maxmin", "mean"))
ARMS <- if (EXPOSURE == "maxmin") c("hot: Tmax DD > 29 C", "cold: Tmin DD < 0 C") else c("hot: daily mean DD > 25 C", "cold: daily mean DD < 5 C")
TAG <- paste0(if (!identical(COUNTRIES, c("USA", "CHN"))) paste0("_", paste(COUNTRIES, collapse = "_")) else "",
              if (EXPOSURE == "mean") "_mean" else "", if (START != 1970) paste0("_from", START) else "",
              if (END != 2019) paste0("_to", END) else "", if (TRENDS == "linear") "_lintrend" else "",
              if (TRIM > 0) paste0("_trim", 100 * TRIM) else "", if (length(DROP)) "_drop" else "")

# --- Daily -> region-year ---------------------------------------------------------
read_cols <- function(v, cols) {
  path <- wcd_download(v, "ERA5", "gadm1", "pop", 2015, "daily", verbose = FALSE)
  x <- as.data.frame(arrow::read_parquet(path, col_select = c("Date", cols)))
  yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
  s <- yr %in% YEARS
  m <- as.matrix(x[s, cols, drop = FALSE]); storage.mode(m) <- "double"
  list(m = m, yr = yr[s])
}
all_cols <- setdiff(names(arrow::open_dataset(wcd_download("tmp", "ERA5", "gadm1", "pop", 2015, "daily", verbose = FALSE))), "Date")
cols <- all_cols[substr(all_cols, 1, 3) %in% COUNTRIES]
tm <- read_cols("tmp", cols); pr <- read_cols("pre", cols)
if (EXPOSURE == "maxmin") { tx <- read_cols("tmpmax", cols); tn <- read_cols("tmpmin", cols) } else { tx <- tn <- tm }
yr <- tm$yr
bad <- colSums(is.na(tm$m) | is.na(tx$m) | is.na(tn$m) | is.na(pr$m)) > 0
if (any(bad)) cat("Dropped (missing days):", paste(cols[bad], collapse = ", "), "\n")
cols <- cols[!bad]
agg <- function(m) rowsum(m[, !bad, drop = FALSE], yr)
Tann <- agg(tm$m) / agg((!is.na(tm$m)) * 1)
Pann <- agg(pmax(pr$m, 0))
hot <- if (EXPOSURE == "maxmin") agg(pmax(tx$m - 29, 0)) else agg(pmax(tm$m - 25, 0))
cold <- if (EXPOSURE == "maxmin") agg(pmax(0 - tn$m, 0)) else agg(pmax(5 - tm$m, 0))
years <- as.integer(rownames(Tann))
gid1 <- sub("^([A-Z]{3})_", "\\1.", cols); gid0 <- substr(cols, 1, 3)

tdev <- function(M, yy) {
  ci <- match(yy, years)
  out <- M[ci, , drop = FALSE]
  for (j in seq_along(ci)) out[j, ] <- M[ci[j], ] - colMeans(M[(ci[j] - WINDOW):(ci[j] - 1), , drop = FALSE])
  out
}
trail_mean <- function(M, yy) {
  ci <- match(yy, years)
  t(vapply(ci, function(k) colMeans(M[(k - WINDOW):(k - 1), , drop = FALSE]), numeric(ncol(M))))
}
dTc <- tdev(Tann, CLIM); dHc <- tdev(hot, CLIM); dCc <- tdev(cold, CLIM)
slope <- function(dy) vapply(seq_along(cols), function(i) {
  x <- dTc[, i] - mean(dTc[, i]); y <- dy[, i] - mean(dy[, i])
  if (sum(y^2) == 0) 0 else sum(x * y) / sum(x^2)
}, 0)
expo <- data.frame(GID_1 = gid1, GID_0 = gid0, Tclim = colMeans(Tann[match(CLIM, years), ]),
                   s_h = slope(dHc), s_c = slope(dCc))
for (c0 in COUNTRIES) {
  i <- expo$GID_0 == c0
  lh <- loess(s_h ~ Tclim, expo[i, ], span = 0.75); lc <- loess(s_c ~ Tclim, expo[i, ], span = 0.75)
  expo$s_h_smooth[i] <- pmax(predict(lh), 0); expo$s_c_smooth[i] <- pmin(predict(lc), 0)
}
# Region-year weather on the estimation years.
yy <- START:END
long <- function(M, name) setNames(data.frame(rep(gid1, each = length(yy)), rep(yy, length(gid1)), as.vector(M)),
                                   c("GID_1", "year", name))
wx <- Reduce(function(a, b) merge(a, b, by = c("GID_1", "year")), list(
  long(tdev(Tann, yy), "dT"), long(tdev(Pann, yy), "dP_mm"),
  long(trail_mean(Tann, yy), "climT"), long(trail_mean(Pann, yy), "climP"),
  long(tdev(hot, yy), "dHDD"), long(tdev(cold, yy), "dCDD")))
sdP <- setNames(apply(Pann[match(CLIM, years), ], 2, sd), gid1)
wx$zP <- wx$dP_mm / sdP[wx$GID_1]

# --- DOSE ------------------------------------------------------------------------------
dose <- utils::read.csv(file.path("data", "DOSE_V2.14.csv"), stringsAsFactors = FALSE)
dose <- dose[dose$GID_0 %in% COUNTRIES & is.finite(dose$grp_pc_lcu_2015) & dose$grp_pc_lcu_2015 > 0,
             c("GID_0", "GID_1", "region", "year", "grp_pc_lcu_2015")]
dose <- dose[order(dose$GID_1, dose$year), ]
lag_by <- function(v, id, t, k) v[match(paste(id, t - k), paste(id, t))]
dose$lgdp <- log(dose$grp_pc_lcu_2015)
dose$g <- dose$lgdp - lag_by(dose$lgdp, dose$GID_1, dose$year, 1)
inc_lags <- sapply(1:5, function(k) lag_by(dose$grp_pc_lcu_2015, dose$GID_1, dose$year, k))
dose$inc_trailing <- log(rowMeans(inc_lags))

# --- Models ----------------------------------------------------------------------------
SPECS <- list(
  "dT"                       = "dT",
  "dT + dT x climate"        = c("dT", "dT_clim"),
  "dT + dT x income"         = c("dT", "dT_inc"),
  "exposure country"         = c("eh", "ec"),
  "exposure smooth"          = c("eh_smooth", "ec_smooth"),
  "exposure x income"        = c("eh", "ec", "eh_inc", "ec_inc"),
  "exposure smooth x income" = c("eh_smooth", "ec_smooth", "eh_smooth_inc", "ec_smooth_inc"),
  "hot only x income"        = c("eh", "eh_inc"),
  "dT + dT x inc + exposure" = c("dT", "dT_inc", "eh", "ec")
)
run_country <- function(c0) {
  d <- merge(dose[dose$GID_0 == c0, ], wx, by = c("GID_1", "year"))
  d <- merge(d, expo[, setdiff(names(expo), "GID_0")], by = "GID_1")
  d <- d[d$year >= START & d$year <= END & is.finite(d$g) & is.finite(d$inc_trailing), ]
  if (length(DROP)) d <- d[!(d$GID_1 %in% DROP | d$region %in% DROP), ]
  if (TRIM > 0) {
    q <- quantile(d$g, c(TRIM, 1 - TRIM)); d <- d[d$g >= q[1] & d$g <= q[2], ]
  }
  d$y <- d$inc_trailing - mean(d$inc_trailing)
  clim_med <- median(tapply(d$climT, d$GID_1, mean))
  d$clim_TM <- d$climT - clim_med; d$clim_RR <- (d$climP - median(tapply(d$climP, d$GID_1, mean))) / 1000
  d$zP_clim <- d$zP * d$clim_RR
  d$dT_clim <- d$dT * d$clim_TM; d$dT_inc <- d$dT * d$y
  for (v in c("", "_smooth")) {
    d[[paste0("eh", v)]] <- d$dT * d[[paste0("s_h", v)]] / 100
    d[[paste0("ec", v)]] <- d$dT * d[[paste0("s_c", v)]] / 100
    d[[paste0("eh", v, "_inc")]] <- d[[paste0("eh", v)]] * d$y
    d[[paste0("ec", v, "_inc")]] <- d[[paste0("ec", v)]] * d$y
  }
  common <- c("zP", "zP_clim", "clim_TM", "clim_RR", "y")
  fit <- function(terms, cl = ~GID_1) {
    f <- as.formula(paste("g ~", paste(c(terms, common), collapse = " + "), if (TRENDS == "quadratic") "| year + GID_1[t, t2]" else "| year + GID_1[t]"))
    d$t <- d$year - 1990; d$t2 <- d$t^2
    m <- feols(f, data = d, cluster = cl, notes = FALSE)
    m$region <- d$GID_1[obs(m)]
    m
  }
  d$t <- d$year - 1990; d$t2 <- d$t^2
  models <- lapply(SPECS, fit)
  G <- length(unique(d$GID_1))
  ref <- models[["dT + dT x climate"]]
  ll <- function(m) { e <- resid(m); s2 <- mean(e^2); -0.5 * log(2 * pi * s2) - e^2 / (2 * s2) }
  vuong <- function(m1, m2) {
    k <- length(coef(m1)) - length(coef(m2))
    D <- tapply(ll(m1) - ll(m2), m1$region, sum)
    (sum(D) - k) / sqrt(sum((D - mean(D))^2) * G / (G - 1))
  }
  cmp <- do.call(rbind, lapply(names(models), function(k) {
    m <- models[[k]]; z <- if (identical(m, ref)) NA else vuong(m, ref)
    data.frame(country = c0, model = k, wr2 = fitstat(m, "wr2")[[1]], dAIC = AIC(m) - AIC(ref),
               vuong_z = z, verdict = if (is.na(z)) "-" else if (abs(z) < 1.96) "indifferent" else if (z > 0) "better" else "worse")
  }))
  wald <- function(m, terms) {
    b <- coef(m)[terms]; V <- vcov(m)[terms, terms, drop = FALSE]
    pf(as.numeric(t(b) %*% solve(V, b)) / length(terms), length(terms), G - 1, lower.tail = FALSE)
  }
  tests <- data.frame(country = c0, test = c("income terms = 0 | exposure x income",
                                               "exposure = 0 | dT, dT x inc",
                                               "dT x inc = 0 | dT, exposure",
                                               "cold arm terms = 0 | exposure x income"),
                      p = c(wald(models[["exposure x income"]], c("eh_inc", "ec_inc")),
                            wald(models[["dT + dT x inc + exposure"]], c("eh", "ec")),
                            wald(models[["dT + dT x inc + exposure"]], "dT_inc"),
                            wald(models[["exposure x income"]], c("ec", "ec_inc"))))
  key <- models[["exposure x income"]]
  key2 <- summary(key, cluster = ~GID_1 + year)
  kt <- intersect(SPECS[["exposure x income"]], names(coef(key)))
  # Leave one region out: the exposure x income coefficients without each
  # region, and the regions that move them most (in SE of the full estimate).
  f_key <- as.formula(paste("g ~", paste(c(SPECS[["exposure x income"]], common), collapse = " + "),
                            if (TRENDS == "quadratic") "| year + GID_1[t, t2]" else "| year + GID_1[t]"))
  names_of <- tapply(d$region, d$GID_1, function(v) v[1])
  loo <- do.call(rbind, lapply(sort(unique(d$GID_1)), function(r) {
    m <- feols(f_key, data = d[d$GID_1 != r, ], notes = FALSE, warn = FALSE)
    data.frame(GID_1 = r, region = names_of[[r]], t(setNames(coef(m)[SPECS[["exposure x income"]]], SPECS[["exposure x income"]])))
  }))
  shift <- sweep(sweep(as.matrix(loo[kt]), 2, coef(key)[kt]),
                 2, se(key)[kt], "/")
  loo$max_shift_se <- apply(abs(shift), 1, max)
  loo$term_moved <- colnames(shift)[apply(abs(shift), 1, which.max)]
  cat(sprintf("\n%s leave-one-region-out, exposure x income (full: %s); range of each coefficient:\n", c0,
              paste(sprintf("%s %.4f", kt, coef(key)[kt]), collapse = ", ")))
  print(apply(loo[kt], 2, range), digits = 3)
  cat("Most influential regions (|shift| in full-sample SE):\n")
  print(head(loo[order(-loo$max_shift_se), ], 5), digits = 3, row.names = FALSE)
  # Where the cold arm has exposure at all.
  ce <- expo[expo$GID_1 %in% d$GID_1, ]
  cat(sprintf("%s regions with cold exposure slope below -1 DD/C: %d of %d (%s)\n", c0, sum(ce$s_c < -1), nrow(ce),
              paste(names_of[ce$GID_1[ce$s_c < -1]], collapse = ", ")))
  coefs <- data.frame(country = c0, term = kt, est = coef(key)[kt],
                      se_region = se(key)[kt], p_region = pvalue(key)[kt],
                      se_2way = se(key2)[kt], p_2way = pvalue(key2)[kt], row.names = NULL)
  # Response to +1 C by region at the 10th, 50th, 90th percentile of income.
  reg <- expo[expo$GID_1 %in% d$GID_1, ]
  yq <- quantile(tapply(d$y, d$GID_1, mean), c(.1, .5, .9))
  resp <- do.call(rbind, lapply(seq_along(yq), function(j) {
    do.call(rbind, lapply(list(c("exposure x income", ""), c("exposure smooth x income", "_smooth")), function(sp) {
      m <- models[[sp[1]]]; v <- sp[2]
      X <- cbind(reg[[paste0("s_h", v)]], reg[[paste0("s_c", v)]], reg[[paste0("s_h", v)]] * yq[j], reg[[paste0("s_c", v)]] * yq[j]) / 100
      colnames(X) <- c(paste0("eh", v), paste0("ec", v), paste0("eh", v, "_inc"), paste0("ec", v, "_inc"))
      X <- X[, intersect(SPECS[[sp[1]]], names(coef(m))), drop = FALSE]
      b <- coef(m)[colnames(X)]; V <- vcov(m)[colnames(X), colnames(X)]
      data.frame(country = c0, model = sp[1], income = sprintf("income p%d", c(10, 50, 90)[j]), income_ratio = exp(yq[j]),
                 GID_1 = reg$GID_1, Tclim = reg$Tclim, est = as.vector(X %*% b), se = sqrt(rowSums((X %*% V) * X)))
    }))
  }))
  cat(sprintf("\n==== %s: %d regions, %d rows, %d-%d; climate %.1f to %.1f C; income range (log, p10-p90) %.2f ====\n",
              c0, G, nrow(d), min(d$year), max(d$year), min(reg$Tclim), max(reg$Tclim), diff(yq[c(1, 3)])))
  etable(models[c("dT + dT x climate", "dT + dT x income", "exposure country", "exposure x income",
                  "exposure smooth x income", "dT + dT x inc + exposure")],
         keep = c("^dT", "^eh", "^ec"), digits = 3, fitstat = ~ n + wr2 + aic) |> print()
  print(cmp, digits = 3, row.names = FALSE)
  print(tests, digits = 3, row.names = FALSE)
  cat("exposure x income, SE by region vs two-way (region + year):\n"); print(coefs, digits = 3, row.names = FALSE)
  list(models = models, cmp = cmp, tests = tests, coefs = coefs, resp = resp, expo = reg, data = d)
}
res <- lapply(setNames(COUNTRIES, COUNTRIES), run_country)

# --- Figures ----------------------------------------------------------------------------
ex <- do.call(rbind, lapply(res, `[[`, "expo"))
el <- rbind(data.frame(ex[c("GID_0", "Tclim")], arm = ARMS[1], s = ex$s_h, sm = ex$s_h_smooth),
            data.frame(ex[c("GID_0", "Tclim")], arm = ARMS[2], s = ex$s_c, sm = ex$s_c_smooth))
nat <- readRDS(file.path(out_dir, if (EXPOSURE == "maxmin") "exposure_mediated_maxmin.rds" else "exposure_mediated.rds"))$exposure
nl <- rbind(data.frame(Tclim = nat$Tclim, arm = ARMS[1], sm = nat$s_h_smooth),
            data.frame(Tclim = nat$Tclim, arm = ARMS[2], sm = nat$s_c_smooth))
p1 <- ggplot(el, aes(Tclim, s, colour = GID_0)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_line(data = nl[order(nl$Tclim), ], aes(Tclim, sm), colour = "grey50", linetype = 2, inherit.aes = FALSE) +
  geom_point(size = 1.5, alpha = .8) +
  geom_line(data = el[order(el$Tclim), ], aes(y = sm), linewidth = .9) +
  facet_wrap(~arm, scales = "free_y") +
  labs(x = "Region climate (C), 1970-2023", y = "exposure slope (DD per C of annual mean T)", colour = NULL,
       title = paste("Regional exposure slopes:", paste(COUNTRIES, collapse = ", ")),
       subtitle = "points: regions; solid: loess within the country; dashed grey: the national (176 countries) smooth") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("C1_subnational_exposure_slopes", TAG, ".png")), p1, width = 9, height = 4.8, dpi = 150)
rp <- do.call(rbind, lapply(res, `[[`, "resp"))
p2 <- ggplot(rp, aes(Tclim, est, colour = income)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(data = rp[rp$model == "exposure x income", ], size = 1.3, alpha = .7) +
  geom_line(data = rp[rp$model == "exposure smooth x income", ], linewidth = .9) +
  geom_ribbon(data = rp[rp$model == "exposure smooth x income" & rp$income == "income p50", ],
              aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), fill = "grey50", alpha = .15, colour = NA) +
  facet_wrap(~country, scales = "free") +
  scale_colour_manual(values = c("firebrick", "black", "steelblue")) +
  labs(x = "Region climate (C), 1970-2023", y = "growth response to a +1 C deviation", colour = NULL,
       title = paste("Exposure x income on subnational data:", sub("^hot: ", "", ARMS[1]), "/", sub("^cold: ", "", ARMS[2])),
       subtitle = sprintf("%d-%d; points: region s_r; lines: smooth s(climate); band: 95%% at median income", START, END)) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("C2_subnational_exposure_income", TAG, ".png")), p2, width = 9.5, height = 5, dpi = 150)
saveRDS(lapply(res, function(r) r[setdiff(names(r), "data")]), file.path(out_dir, paste0("exposure_subnational", TAG, ".rds")))
cat("\nSaved C1/C2 figures and exposure_subnational", TAG, ".rds in ", out_dir, "\n", sep = "")
