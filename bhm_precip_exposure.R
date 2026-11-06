# Precipitation version of A4 (combined extreme exposure against climate) and
# A5 (index deviations vs annual deviations from the 20-year trailing mean).
#
# ERA5 daily precipitation, population weights 2015 (WCD gadm0, mm/day), the
# 176 countries of bhm_microfoundation.R less those with missing days, 1970-2023.
#   wet = sum over days of max(P_d - 20 mm, 0)    heavy-rain excess (mm/yr; ETCCDI R20mm threshold)
#   dry = days with P_d < 1 mm                    dry days (days/yr; ETCCDI)
#   sum = wet / SD_wet + dry / SD_dry             SD = cross-country SD of the 1970-2023
#                                                 climatologies (units differ, equal weight)
# A4P: climatologies against climatological annual P, quadratic across
#      countries, loess, and each country's within-country slope (OLS of the
#      index on P - country mean) over +-2 SD of its annual P.
# A5P: deviation from the trailing mean (t-20..t-1) of each index against that
#      of annual P, in mm and over the panel's sd_RR_all.
#
#   Rscript bhm_precip_exposure.R

suppressPackageStartupMessages(library(ggplot2))
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
r <- readRDS(file.path(out_dir, "results.rds"))

path <- wcd_download("pre", "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)
x <- as.data.frame(arrow::read_parquet(path))
yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
countries <- r$country$GID_0
keep <- yr <= 2023
m <- as.matrix(x[keep, countries]); storage.mode(m) <- "double"; yr <- yr[keep]
bad <- colSums(is.na(m)) > 0
if (any(bad)) cat("Dropped (missing days):", paste(countries[bad], collapse = ", "), "\n")
m <- m[, !bad]; countries <- countries[!bad]
years <- sort(unique(yr))
Pann <- rowsum(m, yr)
wet <- rowsum(pmax(m - 20, 0), yr)
dry <- rowsum((m < 1) * 1, yr)
sel <- years >= 1970

clim <- data.frame(GID_0 = countries, P = colMeans(Pann[sel, ]),
                   wet = colMeans(wet[sel, ]), dry = colMeans(dry[sel, ]))
SD <- c(wet = sd(clim$wet), dry = sd(clim$dry))
cat(sprintf("Cross-country SD of climatologies: wet %.1f mm/yr, dry %.1f days/yr\n", SD[["wet"]], SD[["dry"]]))
clim$sum <- clim$wet / SD[["wet"]] + clim$dry / SD[["dry"]]
INDEX <- list(wet = wet, dry = dry, sum = wet / SD[["wet"]] + dry / SD[["dry"]])
LABEL <- c(wet = "wet: mm above 20 mm/day", dry = "dry: days below 1 mm",
           sum = "sum: wet/SD + dry/SD")

# Within-country slope of each index on annual P (per 100 mm), 1970-2023.
slope_of <- function(s) sapply(countries, function(i) {
  y <- s[sel, i]; p <- Pann[sel, i]
  if (sd(y) == 0) 0 else coef(lm(y ~ p))[2]
})
sdP_country <- apply(Pann[sel, ], 2, sd)

# --- A4P ---
fits <- list(); segs <- list(); curves <- list()
for (k in names(INDEX)) {
  cl <- clim; cl$x <- cl[[k]]
  q <- lm(x ~ P + I(P^2), data = cl)
  quart <- lm(x ~ poly(P, 4, raw = TRUE), data = cl)
  lo <- loess(x ~ P, data = cl)
  sl <- slope_of(INDEX[[k]])
  q_slope <- coef(q)[2] + 2 * coef(q)[3] * cl$P
  cat(sprintf("\n== A4P %s ==\nquadratic in P (mm): R2 %.3f, coef P^2 %.3g (se %.2g), vertex %.0f mm; quartic R2 %.3f, F p = %.2g; loess %s at %.0f mm\n",
              k, summary(q)$r.squared, coef(q)[3], sqrt(vcov(q)[3, 3]), -coef(q)[2] / (2 * coef(q)[3]),
              summary(quart)$r.squared, anova(q, quart)[2, "Pr(>F)"],
              if (k == "dry") "max" else "min", cl$P[if (k == "dry") which.max(fitted(lo)) else which.min(fitted(lo))]))
  cat(sprintf("within-country slopes vs the quadratic's slope line: corr %.2f\n", cor(sl, q_slope)))
  segs[[k]] <- data.frame(index = k, x0 = cl$P - 2 * sdP_country, x1 = cl$P + 2 * sdP_country,
                          y0 = cl$x - 2 * sdP_country * sl, y1 = cl$x + 2 * sdP_country * sl)
  g <- data.frame(index = k, P = seq(min(cl$P), max(cl$P), length.out = 200))
  g$fit <- predict(q, g)
  curves[[k]] <- g
  fits[[k]] <- data.frame(index = k, P = cl$P, x = cl$x)
}
lab <- function(v) factor(v, names(LABEL), LABEL)
pts <- do.call(rbind, fits); sg <- do.call(rbind, segs); cv <- do.call(rbind, curves)
p4 <- ggplot(pts, aes(P, x)) +
  geom_segment(data = sg, aes(x = x0, xend = x1, y = y0, yend = y1), colour = "steelblue", alpha = .6) +
  geom_point(size = 1.1, alpha = .7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey40", linetype = 2, linewidth = .6) +
  geom_line(data = cv, aes(P, fit), colour = "firebrick", linewidth = 1) +
  facet_wrap(~lab(index), scales = "free_y") +
  labs(x = "Climatological annual precipitation (mm), 1970-2023", y = NULL,
       title = "Extreme-precipitation exposure across countries",
       subtitle = "red: quadratic in P; dashed: loess; blue: within-country slope over +-2 SD of annual P") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "A4P_precip_exposure_quadratic.png"), p4, width = 11, height = 4.5, dpi = 150)

# --- A5P ---
trail_dev <- function(s) {
  out <- s * NA
  for (k in which(years >= min(years) + 20)) out[k, ] <- s[k, ] - colMeans(s[(k - 20):(k - 1), , drop = FALSE])
  out[sel, , drop = FALSE]
}
dP <- trail_dev(Pann)
P_panel <- readRDS(file.path("results", "response_explorer", "panels.rds"))
rows <- P_panel$panels[["PWT | ERA5 | population"]]$rows
sdP <- 1000 * tapply(rows$sd_RR_all, rows$GID_0, function(v) v[1])[countries]   # m -> mm
stopifnot(all(is.finite(sdP)))
cat(sprintf("\nCountry SD of annual P (sd_RR_all): median %.0f mm (IQR %.0f-%.0f); as a share of mean P: median %.2f\n",
            median(sdP), quantile(sdP, .25), quantile(sdP, .75), median(sdP / clim$P)))
group <- cut(clim$P, c(-Inf, 500, 1500, Inf), labels = c("arid <500mm", "500-1500mm", "wet >1500mm"))
cat("Countries by group:", paste(names(table(group)), table(group), collapse = ", "), "\n")
summ <- function(k) {
  dx <- trail_dev(INDEX[[k]])
  ok <- apply(dx, 2, sd) > 0
  rc <- sapply(seq_along(countries), function(i) if (!ok[i]) NA else cor(dx[, i], dP[, i]))
  zP <- sweep(dP, 2, sdP, "/"); zx <- sweep(dx, 2, ifelse(ok, apply(dx, 2, sd), 1), "/")
  pooled <- function(a, b) cor(as.vector(a[, ok]), as.vector(b[, ok]))
  by_g <- tapply(rc, group, function(v) sprintf("%.2f [%d/%d]", median(v, na.rm = TRUE), sum(!is.na(v)), length(v)))
  list(row = data.frame(index = k, countries_with_variation = sum(ok),
                        within_median = median(rc, na.rm = TRUE),
                        within_p10 = quantile(rc, .1, na.rm = TRUE), within_p90 = quantile(rc, .9, na.rm = TRUE),
                        pooled_dP_rawX = pooled(dP, dx), pooled_zP_rawX = pooled(zP, dx),
                        pooled_dP_zX = pooled(dP, zx), pooled_zP_zX = pooled(zP, zx),
                        arid = by_g[[1]], mid = by_g[[2]], wet = by_g[[3]], row.names = NULL),
       fig = data.frame(index = k, GID_0 = countries, P = clim$P, r = rc))
}
S5 <- lapply(names(INDEX), summ)
res <- do.call(rbind, lapply(S5, `[[`, "row"))
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], round, 2)
cat("\n== A5P: correlation of index deviation with annual-P deviation (20-year trailing) ==\n")
print(res[, 1:9], row.names = FALSE)
print(res[, c(1, 10:12)], row.names = FALSE)
fig <- do.call(rbind, lapply(S5, `[[`, "fig"))
p5 <- ggplot(fig[is.finite(fig$r), ], aes(P, r)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 1.2, alpha = .7) +
  facet_wrap(~lab(index)) +
  scale_x_log10() +
  labs(x = "Climatological annual precipitation (mm, log scale)", y = "within-country correlation with the P deviation",
       title = "Index deviations vs annual-P deviations, both from the 20-year trailing mean",
       subtitle = "one point per country, 1970-2023; countries with no variation in the index left out") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "A5P_precip_deviation_correlations.png"), p5, width = 10, height = 4, dpi = 150)
saveRDS(list(clim = clim, SD = SD, A5 = res, by_country = fig), file.path(out_dir, "precip_exposure.rds"))
cat("Saved A4P_precip_exposure_quadratic.png and A5P_precip_deviation_correlations.png in", out_dir, "\n")
