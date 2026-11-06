# Combined extreme-temperature exposure, from hot degree days (daily mean
# above 25 C) and cold degree days (daily mean below 5 C), country
# climatologies 1970-2023, against climate: a quadratic fit across countries
# and each country's within-country slope (segments over +-2 SD of its annual
# T).
#   sum (A4): heat + cold, both extremes bad (the BHM argument)
#   net (A3): heat - cold
# Reads results/bhm_microfoundation/results.rds (run bhm_microfoundation.R).
#
#   Rscript bhm_net_exposure.R

suppressPackageStartupMessages(library(ggplot2))
out_dir <- file.path("results", "bhm_microfoundation")
r <- readRDS(file.path(out_dir, "results.rds"))
cl <- r$clim
ct <- r$country

plot_index <- function(sign, label, file) {
  cl$x <- cl$dd_gt25 + sign * cl$dd_lt5
  q <- lm(x ~ T + I(T^2), data = cl)
  cf <- coef(q); se <- sqrt(diag(vcov(q)))
  cat(sprintf("\n== %s ==\nx = %.1f %+.2f T %+.3f T^2  (se %.1f, %.2f, %.3f), R2 %.3f, n %d\n",
              label, cf[1], cf[2], cf[3], se[1], se[2], se[3], summary(q)$r.squared, nrow(cl)))
  cat(sprintf("vertex of the quadratic: %.1f C\n", -cf[2] / (2 * cf[3])))
  cubic <- lm(x ~ poly(T, 3, raw = TRUE), data = cl)
  quart <- lm(x ~ poly(T, 4, raw = TRUE), data = cl)
  cat(sprintf("R2 cubic %.3f, quartic %.3f; F test quadratic vs quartic p = %.2g\n",
              summary(cubic)$r.squared, summary(quart)$r.squared, anova(q, quart)[2, "Pr(>F)"]))
  lo <- loess(x ~ T, data = cl)
  cat(sprintf("loess minimum at %.1f C\n", cl$T[which.min(fitted(lo))]))

  # Within-country slopes (from the country regressions in results.rds).
  ct$slope <- r$S[ct$GID_0, "dd_gt25"] + sign * r$S[ct$GID_0, "dd_lt5"]
  ct$curv <- 2 * (r$C[ct$GID_0, "dd_gt25"] + sign * r$C[ct$GID_0, "dd_lt5"])
  ct$xc <- cl$x[match(ct$GID_0, cl$GID_0)]
  q_slope <- cf[2] + 2 * cf[3] * ct$Tclim
  cat(sprintf("within-country slopes vs the quadratic's 2c-line: corr %.2f; quadratic 2nd derivative %.2f, ",
              cor(ct$slope, q_slope), 2 * cf[3]))
  cat(sprintf("within-country 2nd derivative median %.2f (IQR %.2f..%.2f)\n",
              median(ct$curv), quantile(ct$curv, .25), quantile(ct$curv, .75)))
  seg <- transform(ct, x0 = Tclim - 2 * sd_T, x1 = Tclim + 2 * sd_T,
                   y0 = xc - 2 * sd_T * slope, y1 = xc + 2 * sd_T * slope)

  grid <- data.frame(T = seq(min(cl$T), max(cl$T), length.out = 200))
  grid$fit <- predict(q, grid)
  p <- ggplot(cl, aes(T, x)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_segment(data = seg, aes(x = x0, xend = x1, y = y0, yend = y1), colour = "steelblue", alpha = .6) +
    geom_point(size = 1.3, alpha = .7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey40", linetype = 2, linewidth = .6) +
    geom_line(data = grid, aes(T, fit), colour = "firebrick", linewidth = 1) +
    labs(x = "Climatological annual mean T (C), 1970-2023",
         y = paste("degree days per year:", label),
         title = "Combined extreme-temperature exposure across countries",
         subtitle = sprintf("red: quadratic %.0f %+.1f T %+.2f T^2 (R2 %.2f); dashed: loess; blue: within-country slope over +-2 SD",
                            cf[1], cf[2], cf[3], summary(q)$r.squared)) +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, file), p, width = 8.5, height = 5.5, dpi = 150)
  cat("Saved", file.path(out_dir, file), "\n")
}
plot_index(+1, "above 25 C plus below 5 C", "A4_sum_exposure_quadratic.png")
plot_index(-1, "above 25 C minus below 5 C", "A3_net_exposure_quadratic.png")
