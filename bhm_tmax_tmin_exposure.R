# Exposure with daily max / min temperature: degree days of Tmin below 0 C
# and of Tmax above 29 C, their sum, and frost days (Tmin < 0 C) for
# comparison, against the annual mean T. ERA5 daily, population weights 2015, countries (WCD gadm0), the 176
# countries of the bhm_microfoundation.R panel.
#
# A8 climatologies 1970-2023 against climate: quadratic across countries,
#    loess, and each country's within-country slope over +-2 SD (as A4).
# A9 exposure slopes s_i: within-country OLS slope of the index's deviation
#    from its 20-year trailing mean on that of annual mean T (1970-2023), per
#    C, and their loess on climate (the smooth version of
#    bhm_exposure_mediated.R).
#
#   Rscript bhm_tmax_tmin_exposure.R

suppressPackageStartupMessages(library(ggplot2))
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
r <- readRDS(file.path(out_dir, "results.rds"))
countries <- r$country$GID_0
YEARS <- 1940:2023; CLIM <- 1970:2023; WINDOW <- 20
FROST <- 0; HOT <- 29

daily <- function(v) {
  x <- as.data.frame(arrow::read_parquet(wcd_download(v, "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)))
  yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
  s <- yr %in% YEARS
  m <- as.matrix(x[s, countries]); storage.mode(m) <- "double"
  list(m = m, yr = yr[s])
}
tm <- daily("tmp"); tx <- daily("tmpmax"); tn <- daily("tmpmin")
stopifnot(identical(tm$yr, tx$yr), identical(tm$yr, tn$yr))
ok <- colSums(is.na(tm$m) | is.na(tx$m) | is.na(tn$m)) == 0
if (any(!ok)) cat("Dropped (missing days):", paste(countries[!ok], collapse = ", "), "\n")
countries <- countries[ok]
yr <- tm$yr
agg <- function(m) rowsum(m[, ok, drop = FALSE], yr)        # years x countries
Tann <- agg(tm$m) / agg((!is.na(tm$m)) * 1)
# Degree days of Tmin below 0 C weigh a -50 C night more than a -3 C one;
# with degree days of Tmax above 29 C, both extremes are in degree days and
# their sum needs no rescaling. Frost days kept for comparison.
IDX <- list(
  frost = agg((tn$m < FROST) * 1),
  frost_dd = agg(pmax(FROST - tn$m, 0)),
  hot_dd = agg(pmax(tx$m - HOT, 0)))
IDX$sum_dd <- IDX$frost_dd + IDX$hot_dd
LAB <- c(frost = "frost days (Tmin < 0 C)", frost_dd = "degree days of Tmin below 0 C",
         hot_dd = "degree days of Tmax above 29 C", sum_dd = "sum: Tmin DD < 0 + Tmax DD > 29")
COLD_IDX <- c("frost", "frost_dd")
years <- as.integer(rownames(Tann)); ci <- match(CLIM, years)

# Trailing deviations (t-20..t-1) on the CLIM years.
tdev <- function(M) {
  out <- M[ci, , drop = FALSE]
  for (j in seq_along(ci)) out[j, ] <- M[ci[j], ] - colMeans(M[(ci[j] - WINDOW):(ci[j] - 1), , drop = FALSE])
  out
}
dT <- tdev(Tann)
Tclim <- colMeans(Tann[ci, ]); sdT <- apply(Tann[ci, ], 2, sd)
slope_of <- function(dy) vapply(seq_along(countries), function(i) {
  x <- dT[, i] - mean(dT[, i]); y <- dy[, i] - mean(dy[, i])
  if (sum(y^2) == 0) 0 else sum(x * y) / sum(x^2)
}, 0)
corr_of <- function(dy) vapply(seq_along(countries), function(i) {
  if (sd(dy[, i]) == 0) NA_real_ else cor(dT[, i], dy[, i])
}, 0)

pts <- segs <- curves <- slopes <- list(); summ <- list()
for (k in names(IDX)) {
  y <- IDX[[k]]; ybar <- colMeans(y[ci, ])
  dy <- tdev(y); s <- slope_of(dy); rc <- corr_of(dy)
  # Slope of the level on the annual mean (for the A8 segments, as in A4).
  s_lev <- vapply(seq_along(countries), function(i) {
    x <- Tann[ci, i] - mean(Tann[ci, i]); z <- y[ci, i] - mean(y[ci, i])
    if (sum(z^2) == 0) 0 else sum(x * z) / sum(x^2)
  }, 0)
  q <- lm(ybar ~ Tclim + I(Tclim^2)); q4 <- lm(ybar ~ poly(Tclim, 4, raw = TRUE))
  lo <- loess(s ~ Tclim, span = 0.5)
  sm <- predict(lo)
  # Sign restrictions only where the index has one sign (not the sum).
  sm <- if (k %in% COLD_IDX) pmin(sm, 0) else if (k == "sum_dd") sm else pmax(sm, 0)
  pts[[k]] <- data.frame(index = k, T = Tclim, y = ybar)
  segs[[k]] <- data.frame(index = k, x0 = Tclim - 2 * sdT, x1 = Tclim + 2 * sdT,
                          y0 = ybar - 2 * sdT * s_lev, y1 = ybar + 2 * sdT * s_lev)
  g <- data.frame(index = k, T = seq(min(Tclim), max(Tclim), length.out = 200)); g$y <- predict(q, data.frame(Tclim = g$T))
  curves[[k]] <- g
  slopes[[k]] <- data.frame(index = k, GID_0 = countries, T = Tclim, s = s, s_smooth = sm, r = rc)
  grp <- cut(Tclim, c(-Inf, 10, 20, Inf), labels = c("cold <10C", "temperate 10-20C", "hot >20C"))
  summ[[k]] <- data.frame(index = k, R2_quadratic = summary(q)$r.squared, R2_quartic = summary(q4)$r.squared,
                          p_quartic = anova(q, lm(ybar ~ poly(Tclim, 4, raw = TRUE)))[2, "Pr(>F)"],
                          vertex = -coef(q)[2] / (2 * coef(q)[3]),
                          slope_cold = median(s[grp == "cold <10C"]), slope_temperate = median(s[grp == "temperate 10-20C"]),
                          slope_hot = median(s[grp == "hot >20C"]),
                          r_cold = median(rc[grp == "cold <10C"], na.rm = TRUE),
                          r_temperate = median(rc[grp == "temperate 10-20C"], na.rm = TRUE),
                          r_hot = median(rc[grp == "hot >20C"], na.rm = TRUE),
                          zero_exposure_countries = sum(apply(dy, 2, sd) == 0), row.names = NULL)
}
summ <- do.call(rbind, summ)
num <- sapply(summ, is.numeric); summ[num] <- lapply(summ[num], function(z) signif(z, 3))
cat(sprintf("%d countries, climatology %d-%d, slopes on %d-year trailing deviations\n",
            length(countries), min(CLIM), max(CLIM), WINDOW))
cat("slope_* = median within-country exposure slope (per C of annual mean T) by climate group; r_* = median correlation of the trailing deviations\n")
print(summ, row.names = FALSE)

lab <- function(v) factor(v, names(LAB), LAB)
P <- do.call(rbind, pts); S <- do.call(rbind, segs); Cq <- do.call(rbind, curves); SL <- do.call(rbind, slopes)
p8 <- ggplot(P, aes(T, y)) +
  geom_segment(data = S, aes(x = x0, xend = x1, y = y0, yend = y1), colour = "steelblue", alpha = .6) +
  geom_point(size = 1.2, alpha = .7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey40", linetype = 2, linewidth = .6) +
  geom_line(data = Cq, colour = "firebrick", linewidth = 1) +
  facet_wrap(~lab(index), scales = "free_y") +
  labs(x = "Climatological annual mean T (C), 1970-2023", y = "per year (1970-2023 mean)",
       title = "Exposure from daily max / min temperature across countries",
       subtitle = "red: quadratic; dashed: loess; blue: within-country slope over +-2 SD of annual T") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "A8_tmax_tmin_exposure.png"), p8, width = 9.5, height = 6.5, dpi = 150)
p9 <- ggplot(SL, aes(T, s)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 1.2, alpha = .7) +
  geom_line(data = SL[order(SL$T), ], aes(y = s_smooth), colour = "steelblue", linewidth = 1) +
  facet_wrap(~lab(index), scales = "free_y") +
  labs(x = "Climate (C), 1970-2023", y = "exposure slope s_i (per C of annual mean T)",
       title = "Exposure slopes: country (points) and smooth in climate (blue)",
       subtitle = "within-country OLS of the index deviation on the annual-T deviation, both from the 20-year trailing mean") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "A9_tmax_tmin_exposure_slopes.png"), p9, width = 9.5, height = 6.5, dpi = 150)
saveRDS(list(summary = summ, slopes = SL, climatology = P), file.path(out_dir, "tmax_tmin_exposure.rds"))
cat("Saved A8_tmax_tmin_exposure.png and A9_tmax_tmin_exposure_slopes.png in", out_dir, "\n")
