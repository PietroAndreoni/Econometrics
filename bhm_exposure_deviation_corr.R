# Do year-to-year swings in the hot / cold degree-day indexes track the
# temperature deviation from the 20-year trailing mean?
#
# Country-year, 1970-2023, the 176 countries of the bhm_microfoundation.R
# panel. ERA5 daily mean T, population weights 2015 (WCD gadm0).
#   dT   = T - mean(T, t-20..t-1)                (the panel's lag_20 baseline)
#   zT   = dT / sd_TM_all                        (the panel's country SD of annual T)
#   dX   = X - mean(X, t-20..t-1) for X = hot (DD > 25 C), cold (DD < 5 C), sum
# Within-country correlations are scale-free (dT and zT give the same); the
# pooled ones differ, so they are shown for both, with X raw and X over its
# own country SD.
#
#   Rscript bhm_exposure_deviation_corr.R

source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
r <- readRDS(file.path(out_dir, "results.rds"))
countries <- r$country$GID_0
Tclim <- setNames(r$country$Tclim, countries)

path <- wcd_download("tmp", "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)
x <- as.data.frame(arrow::read_parquet(path))
yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
m <- as.matrix(x[, countries]); storage.mode(m) <- "double"
keep <- yr <= 2023
m <- m[keep, ]; yr <- yr[keep]
nd <- rowsum((!is.na(m)) * 1, yr)
Tm <- rowsum(m, yr) / nd
hot <- rowsum(pmax(m - 25, 0), yr)
cold <- rowsum(pmax(5 - m, 0), yr)
stopifnot(all(nd >= 365))
years <- as.integer(rownames(Tm))

# Deviation from the trailing 20-year mean (t-20..t-1).
trail_dev <- function(s) {
  out <- s * NA
  for (k in which(years >= min(years) + 20)) out[k, ] <- s[k, ] - colMeans(s[(k - 20):(k - 1), , drop = FALSE])
  out
}
sel <- years >= 1970
dev <- lapply(list(T = Tm, hot = hot, cold = cold, sum = hot + cold), function(s) trail_dev(s)[sel, , drop = FALSE])

# The panel's country SD of annual T (sd_TM_all).
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
rows <- P$panels[["PWT | ERA5 | population"]]$rows
sdT <- tapply(rows$sd_TM_all, rows$GID_0, function(v) v[1])[countries]
stopifnot(all(is.finite(sdT)))
cat(sprintf("Country SD of annual T (sd_TM_all): median %.2f C (IQR %.2f-%.2f)\n",
            median(sdT), quantile(sdT, .25), quantile(sdT, .75)))

group <- cut(Tclim, c(-Inf, 10, 20, Inf), labels = c("cold <10C", "temperate 10-20C", "hot >20C"))
summ <- function(idx) {
  dx <- dev[[idx]]; dT <- dev$T
  rc <- sapply(countries, function(i) {
    if (sd(dx[, i]) == 0) NA else cor(dx[, i], dT[, i])
  })
  zT <- sweep(dT, 2, sdT, "/")
  zx <- sweep(dx, 2, apply(dx, 2, sd), "/")
  ok <- apply(dx, 2, sd) > 0
  pooled <- function(a, b) cor(as.vector(a[, ok]), as.vector(b[, ok]))
  by_g <- tapply(rc, group, function(v) sprintf("%.2f [%d/%d]", median(v, na.rm = TRUE), sum(!is.na(v)), length(v)))
  data.frame(index = idx,
             countries_with_variation = sum(ok),
             within_median = median(rc, na.rm = TRUE),
             within_p10 = quantile(rc, .1, na.rm = TRUE), within_p90 = quantile(rc, .9, na.rm = TRUE),
             pooled_dT_rawX = pooled(dT, dx), pooled_zT_rawX = pooled(zT, dx),
             pooled_dT_zX = pooled(dT, zx), pooled_zT_zX = pooled(zT, zx),
             cold_countries = by_g[[1]], temperate = by_g[[2]], hot_countries = by_g[[3]],
             row.names = NULL)
}
res <- do.call(rbind, lapply(c("hot", "cold", "sum"), summ))
num <- sapply(res, is.numeric)
res[num] <- lapply(res[num], round, 2)
cat("\nCorrelation of the index deviation with the T deviation (20-year trailing baseline)\n")
cat("within_*: country-by-country correlation (scale-free); by group: median [countries with variation / all]\n")
print(res[, 1:9], row.names = FALSE)
print(res[, c(1, 10:12)], row.names = FALSE)

# Within-country correlation against climate, for the figure.
fig <- do.call(rbind, lapply(c("hot", "cold", "sum"), function(idx) {
  dx <- dev[[idx]]
  data.frame(index = idx, GID_0 = countries, Tclim = as.numeric(Tclim),
             r = sapply(countries, function(i) if (sd(dx[, i]) == 0) NA else cor(dx[, i], dev$T[, i])))
}))
suppressPackageStartupMessages(library(ggplot2))
p <- ggplot(fig[is.finite(fig$r), ], aes(Tclim, r)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 1.2, alpha = .7) +
  facet_wrap(~factor(index, c("hot", "cold", "sum"),
                     c("hot: DD > 25 C", "cold: DD < 5 C", "sum: hot + cold"))) +
  labs(x = "Climate (C)", y = "within-country correlation with the T deviation",
       title = "Index deviations vs annual-T deviations, both from the 20-year trailing mean",
       subtitle = "one point per country, 1970-2023; countries with no variation in the index left out") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "A5_deviation_correlations.png"), p, width = 9, height = 4, dpi = 150)
saveRDS(list(summary = res, by_country = fig), file.path(out_dir, "deviation_correlations.rds"))
cat("Saved", file.path(out_dir, "A5_deviation_correlations.png"), "\n")
