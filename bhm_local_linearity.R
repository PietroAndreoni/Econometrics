# Is exposure linear in the annual mean within countries? For T (hot = DD >
# 25 C, cold = DD < 5 C, sum) and P (wet = mm above 20 mm/day, dry = days
# below 1 mm, sum = each over its cross-country SD, as in
# bhm_precip_exposure.R). ERA5 daily, population weights 2015, the countries
# of bhm_microfoundation.R, 1970-2023.
#
# Country by country: index ~ x + x^2, x = annual value - country mean.
#   r2_gain    R2 of the quadratic less that of the linear fit
#   share_p05  share of countries whose x^2 term has p < 0.05 (5% if linear)
#   kappa      |slope at +1 SD - slope at -1 SD| / |slope at the mean|
#              = |4 c SD / b|: 0 linear, 1 the slope changes by its own size
#   convex     share of countries with c > 0
# Countries whose index never varies are left out; kappa only where the
# linear slope is significant (p < 0.05), else the ratio is noise.
# Figure: deciles of the within-country deviation (country SD units) against
# the mean index deviation (index SD units), by climate group.
#
#   Rscript bhm_local_linearity.R

suppressPackageStartupMessages(library(ggplot2))
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
r <- readRDS(file.path(out_dir, "results.rds"))
countries0 <- r$country$GID_0

daily <- function(v) {
  x <- as.data.frame(arrow::read_parquet(wcd_download(v, "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)))
  yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
  s <- yr >= 1970 & yr <= 2023
  m <- as.matrix(x[s, countries0]); storage.mode(m) <- "double"
  list(m = m, yr = yr[s])
}
tm <- daily("tmp"); pr <- daily("pre")
ok <- colSums(is.na(tm$m)) == 0 & colSums(is.na(pr$m)) == 0
countries <- countries0[ok]
cat("Countries:", length(countries), "(dropped:", paste(countries0[!ok], collapse = ", "), ")\n")
Tm <- rowsum(tm$m[, ok], tm$yr) / rowsum((!is.na(tm$m[, ok])) * 1, tm$yr)
Pa <- rowsum(pr$m[, ok], pr$yr)
hot <- rowsum(pmax(tm$m[, ok] - 25, 0), tm$yr); cold <- rowsum(pmax(5 - tm$m[, ok], 0), tm$yr)
wet <- rowsum(pmax(pr$m[, ok] - 20, 0), pr$yr); dry <- rowsum((pr$m[, ok] < 1) * 1, pr$yr)
sdw <- sd(colMeans(wet)); sdd <- sd(colMeans(dry))
IDX <- list(
  T = list(hot = hot, cold = cold, sum = hot + cold),
  P = list(wet = wet, dry = dry, sum = wet / sdw + dry / sdd))
X <- list(T = Tm, P = Pa)
clim <- list(T = colMeans(Tm), P = colMeans(Pa))
groups <- list(T = cut(clim$T, c(-Inf, 10, 20, Inf), labels = c("cold <10C", "temperate 10-20C", "hot >20C")),
               P = cut(clim$P, c(-Inf, 500, 1500, Inf), labels = c("arid <500mm", "500-1500mm", "wet >1500mm")))

country_fit <- function(y, x) {
  if (sd(y) == 0) return(c(r2_gain = NA, p2 = NA, kappa = NA, c = NA, p1 = NA, rel_sd = NA))
  xc <- x - mean(x); s <- sd(x)
  l <- lm(y ~ xc); q <- lm(y ~ xc + I(xc^2))
  b <- unname(coef(q)[2]); cc <- unname(coef(q)[3])
  c(r2_gain = summary(q)$r.squared - summary(l)$r.squared,
    p2 = summary(q)$coefficients[3, 4],
    kappa = abs(4 * cc * s / b), c = cc,
    p1 = summary(l)$coefficients[2, 4], rel_sd = s / abs(mean(x)))
}
tab <- list(); bins <- list(); bins_abs <- list()
for (v in names(IDX)) for (k in names(IDX[[v]])) {
  y <- IDX[[v]][[k]]; x <- X[[v]]
  f <- t(sapply(countries, function(i) country_fit(y[, i], x[, i])))
  g <- groups[[v]]
  row <- function(sel, lab) {
    ff <- f[sel & !is.na(f[, "p2"]), , drop = FALSE]
    kk <- ff[ff[, "p1"] < 0.05, "kappa"]
    data.frame(var = v, index = k, group = lab, n = nrow(ff),
               r2_gain = median(ff[, "r2_gain"]), share_p05 = mean(ff[, "p2"] < 0.05),
               kappa_median = median(kk), kappa_p75 = quantile(kk, .75, names = FALSE),
               n_kappa = length(kk), convex = mean(ff[, "c"] > 0))
  }
  tab[[paste(v, k)]] <- rbind(row(rep(TRUE, length(countries)), "all"),
                              do.call(rbind, lapply(levels(g), function(l) row(g == l, l))))
  # Binned within-country shape, in country SD units.
  zx <- sweep(sweep(x, 2, colMeans(x)), 2, apply(x, 2, sd), "/")
  sy <- apply(y, 2, sd); use <- sy > 0
  zy <- sweep(sweep(y, 2, colMeans(y)), 2, ifelse(use, sy, 1), "/")
  dd <- data.frame(z = as.vector(zx[, use]), y = as.vector(zy[, use]),
                   group = rep(g[use], each = nrow(x)))
  dd$bin <- cut(dd$z, quantile(dd$z, seq(0, 1, .1)), include.lowest = TRUE)
  b <- aggregate(cbind(z, y) ~ bin + group, dd, mean)
  b$var <- v; b$index <- k
  bins[[paste(v, k)]] <- b
  # Absolute version (A6b): the index deviation in its own units, every
  # country (an index that never moves is a real zero change in exposure).
  da <- data.frame(z = as.vector(zx), y = as.vector(sweep(y, 2, colMeans(y))),
                   group = rep(g, each = nrow(x)))
  da$bin <- cut(da$z, quantile(da$z, seq(0, 1, .1)), include.lowest = TRUE)
  ba <- aggregate(cbind(z, y) ~ bin + group, da, mean)
  ba$var <- v; ba$index <- k
  bins_abs[[paste(v, k)]] <- ba
}
tab <- do.call(rbind, tab); rownames(tab) <- NULL
num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], function(z) round(z, 3))
cat(sprintf("Relative SD of the annual value (SD / |mean|): T median %.3f, P median %.3f\n",
            median(apply(Tm, 2, sd) / abs(colMeans(Tm))), median(apply(Pa, 2, sd) / colMeans(Pa))))
cat("\nWithin-country curvature of exposure in the annual mean (see header)\n")
print(tab, row.names = FALSE)

bins <- do.call(rbind, bins)
bins$panel <- factor(paste(bins$var, bins$index), c("T hot", "T cold", "T sum", "P wet", "P dry", "P sum"))
p <- ggplot(bins, aes(z, y, colour = group)) +
  geom_abline(slope = 0, intercept = 0, linetype = 3) +
  geom_line() + geom_point(size = 1.3) +
  facet_wrap(~panel, nrow = 2) +
  labs(x = "within-country deviation of the annual mean (country SD), deciles",
       y = "mean index deviation (country SD of the index)", colour = NULL,
       title = "Within-country shape of exposure in the annual mean",
       subtitle = "a straight line through the deciles = locally linear; groups: T by climate (C), P by annual mm") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "A6_local_linearity.png"), p, width = 10, height = 6.5, dpi = 150)
bins_abs <- do.call(rbind, bins_abs)
UNITS <- c("T hot" = "T hot (degree days > 25 C)", "T cold" = "T cold (degree days < 5 C)",
           "T sum" = "T sum (degree days)", "P wet" = "P wet (mm above 20 mm/day)",
           "P dry" = "P dry (days < 1 mm)", "P sum" = "P sum (cross-country SD units)")
bins_abs$panel <- factor(paste(bins_abs$var, bins_abs$index), names(UNITS), UNITS)
pa <- ggplot(bins_abs, aes(z, y, colour = group)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_line() + geom_point(size = 1.3) +
  facet_wrap(~panel, nrow = 2, scales = "free_y") +
  labs(x = "within-country deviation of the annual mean (country SD), deciles",
       y = "mean index deviation from the country mean (index units)", colour = NULL,
       title = "Within-country shape of exposure in the annual mean: absolute exposure",
       subtitle = "all countries, including those whose index never moves; groups: T by climate (C), P by annual mm") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "A6b_local_linearity_absolute.png"), pa, width = 10, height = 6.5, dpi = 150)
# Slope of absolute exposure per country SD of the annual mean, below and
# above the mean (the arms), by group: what a +-1 SD year does.
arms <- do.call(rbind, lapply(split(bins_abs, list(bins_abs$panel, bins_abs$group), drop = TRUE), function(b) {
  lo <- b[b$z < 0, ]; hi <- b[b$z > 0, ]
  data.frame(panel = b$panel[1], group = b$group[1],
             slope_below = unname(coef(lm(y ~ z, lo))[2]), slope_above = unname(coef(lm(y ~ z, hi))[2]))
}))
arms$ratio_above_below <- arms$slope_above / arms$slope_below
num <- sapply(arms, is.numeric); arms[num] <- lapply(arms[num], function(z) signif(z, 3))
cat("\nAbsolute exposure change per country SD, deciles below vs above the mean\n")
print(arms, row.names = FALSE)
saveRDS(list(table = tab, bins = bins, bins_abs = bins_abs, arms = arms), file.path(out_dir, "local_linearity.rds"))
cat("Saved", file.path(out_dir, "A6_local_linearity.png"), "and A6b_local_linearity_absolute.png\n")
