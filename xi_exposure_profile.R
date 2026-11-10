# The income elasticity xi of the dollar loss, for T and for P, by the AIC
# profile, in the exposure specification of xi_exposure_setup.R (linear T
# arms through Tmax > 29 C and Tmin < 0 C exposure; quadratic P arms with
# climate and agricultural share; country and year effects, no trends;
# 20-year trailing deviations and characteristics; last year's log GDP).
#
# Given xi the model is linear, so each point of the profile is one OLS fit
# with the same regressors. The AIC then differs from the residual sum of
# squares only by constants, so its minimum is the NLS estimate, and
# AIC(xi0) - AIC(xi_hat) is the likelihood-ratio statistic of xi = xi0
# (chi2(1): 3.84 at 5%) under iid errors. It ignores the clustering by
# country, so it over-rejects: xi_exposure_bootstrap.R is the check.
#
# Profiles: xi_T with xi_P at its estimate, and xi_P with xi_T at its estimate
# (the joint NLS estimate by coordinates).
#
#   Rscript xi_exposure_profile.R
#   EXPOSURE_MODE=smooth Rscript xi_exposure_profile.R   # smooth exposure slopes

source("xi_exposure_setup.R")
aic_at <- function(xi_T, xi_P) stats::AIC(fit_xi(xi_T, xi_P))
t0 <- Sys.time()
hat <- nls_xi()
cat(sprintf("Exposure mode: %s. %d rows, %d countries. Joint NLS: xi_T = %.3f, xi_P = %.3f (%.0f s)\n",
            MODE, nrow(D), length(unique(D$GID_0)), hat[["TM"]], hat[["RR"]],
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
GRID <- seq(-3, 3, by = 0.05)
prof <- rbind(
  data.frame(variable = "TM", xi = GRID, other_xi = hat[["RR"]], aic = vapply(GRID, function(z) aic_at(z, hat[["RR"]]), 0)),
  data.frame(variable = "RR", xi = GRID, other_xi = hat[["TM"]], aic = vapply(GRID, function(z) aic_at(hat[["TM"]], z), 0)))
a_hat <- aic_at(hat[["TM"]], hat[["RR"]])
summ <- do.call(rbind, lapply(c("TM", "RR"), function(v) {
  lr <- function(x0) (if (v == "TM") aic_at(x0, hat[["RR"]]) else aic_at(hat[["TM"]], x0)) - a_hat
  p <- prof[prof$variable == v, ]
  inside <- p$xi[p$aic - a_hat <= stats::qchisq(0.95, 1)]
  data.frame(variable = v, xi_hat = hat[[v]], lr_xi0 = lr(0), p_xi0 = stats::pchisq(lr(0), 1, lower.tail = FALSE),
             lr_xi1 = lr(1), p_xi1 = stats::pchisq(lr(1), 1, lower.tail = FALSE),
             lr_ci95_lo = min(inside), lr_ci95_hi = max(inside),
             ci_open = min(inside) == min(GRID) || max(inside) == max(GRID))
}))
cat("\nAIC profile (LR = AIC(xi0) - AIC(xi_hat), chi2(1) under iid errors; the other xi at its estimate):\n")
print(summ, digits = 3, row.names = FALSE)

m_hat <- fit_xi(hat[["TM"]], hat[["RR"]], vcov = ~GID_0)
ct <- coeftable(m_hat)
rownames(ct) <- ifelse(rownames(ct) %in% names(m_hat$dev_names), m_hat$dev_names[rownames(ct)], rownames(ct))
cat("\nThe model at the estimate (clustered by country):\n")
print(ct[c(DEV, intersect(c("conv_gdp", "inc"), rownames(ct))), ], digits = 3)
saveRDS(list(mode = MODE, spec = SPEC, other = OTHER, grid = GRID, profiles = prof, summary = summ, xi_hat = hat,
             aic_hat = a_hat, coefs = ct),
        file.path(out_dir, paste0("profile_", MODE, ".rds")))

suppressPackageStartupMessages(library(ggplot2))
prof$panel <- ifelse(prof$variable == "TM", sprintf("Temperature: xi_T (xi_P = %.2f)", hat[["RR"]]),
                     sprintf("Precipitation: xi_P (xi_T = %.2f)", hat[["TM"]]))
prof$daic <- prof$aic - a_hat
p <- ggplot(prof, aes(xi, daic)) +
  geom_hline(yintercept = stats::qchisq(0.95, 1), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = c(0, 1), colour = "grey80") +
  geom_line(linewidth = 1, colour = "#2a78d6") +
  facet_wrap(~panel) + coord_cartesian(ylim = c(0, 20)) +
  labs(x = "xi (dollar loss per person proportional to Y^xi)", y = "AIC less its minimum",
       title = sprintf("AIC profile of xi, exposure specification (%s exposure slopes)", MODE),
       subtitle = "Dashed: 3.84, the 5% LR threshold under iid errors (ignores the clustering by country)") +
  theme_minimal(base_size = 12)
png(file.path(out_dir, paste0("profile_", MODE, ".png")), width = 1100, height = 480); print(p); invisible(dev.off())
cat("Saved", file.path(out_dir, paste0("profile_", MODE, ".{rds,png}")), "\n")
