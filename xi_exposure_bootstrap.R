# Country-block bootstrap of the income elasticities (xi_T, xi_P) of the
# exposure specification (xi_exposure_setup.R): countries resampled with
# replacement (each copy its own country, with its own fixed effects), and in
# each sample (xi_T, xi_P) re-estimated by NLS (nls_xi(): grids of step 0.5
# over XI_BOUNDS and continuous searches, by coordinates, from the full-sample
# estimate). Percentile 95% intervals and two-sided bootstrap p-values of
# xi = 0 (equal dollars) and xi = 1 (equal %), the share of draws at a bound,
# and the joint distribution.
#
#   Rscript xi_exposure_bootstrap.R            # B = 999, all cores but one
#   B=199 Rscript xi_exposure_bootstrap.R
#   EXPOSURE_MODE=smooth Rscript xi_exposure_bootstrap.R

source("xi_exposure_setup.R")
B <- as.integer(Sys.getenv("B", "999"))
SEED <- 20261109
t0 <- Sys.time()
hat <- nls_xi(D)
cat(sprintf("Exposure mode %s. Full sample: xi_T = %.3f, xi_P = %.3f (%.0f s)\n", MODE, hat[["TM"]], hat[["RR"]],
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

countries <- unique(D$GID_0)
rows_of <- split(seq_len(nrow(D)), D$GID_0)
set.seed(SEED)
draws <- lapply(seq_len(B), function(b) sample(countries, replace = TRUE))
one_draw <- function(draw) {
  db <- D[unlist(rows_of[draw], use.names = FALSE), ]
  id <- rep(paste0(draw, "_", seq_along(draw)), lengths(rows_of[draw]))
  db$GID_0 <- db$GID_1 <- id
  tryCatch(nls_xi(db, start = hat), error = function(e) c(TM = NA_real_, RR = NA_real_))
}
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(fixest)))
parallel::clusterExport(cl, c("D", "DEV", "OTHER", "FE", "XI_BOUNDS", "with_terms", "fit_xi", "rss_at", "nls_xi",
                              "rows_of", "hat"))
cat(sprintf("Bootstrap: %d draws on %d cores...\n", B, n_cores))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, draws, one_draw)
parallel::stopCluster(cl)
boot <- do.call(rbind, res)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ok <- stats::complete.cases(boot)
bp <- function(x, x0) min(1, 2 * min(mean(x <= x0), mean(x >= x0)))
summ <- do.call(rbind, lapply(c("TM", "RR"), function(v) {
  x <- boot[ok, v]
  data.frame(variable = v, estimate = hat[[v]], boot_median = stats::median(x), boot_se = stats::sd(x),
             lo95 = unname(stats::quantile(x, 0.025)), hi95 = unname(stats::quantile(x, 0.975)),
             at_bound = mean(x %in% XI_BOUNDS), p_xi0 = bp(x, 0), p_xi1 = bp(x, 1))
}))
cat(sprintf("\nCountry-block bootstrap, %d of %d draws (percentile 95%% CI, two-sided bootstrap p):\n", sum(ok), B))
print(summ, digits = 3, row.names = FALSE)
cat(sprintf("corr(xi_T, xi_P) across draws: %.2f\n", stats::cor(boot[ok, 1], boot[ok, 2])))
cat("Share of draws with xi_P > 0.75 (the second mode of the P profile near 1.2):",
    round(mean(boot[ok, "RR"] > 0.75), 3), "\n")
saveRDS(list(mode = MODE, spec = SPEC, other = OTHER, hat = hat, draws = boot, summary = summ, B = B, seed = SEED,
             bounds = XI_BOUNDS), file.path(out_dir, paste0("bootstrap_", MODE, ".rds")))

suppressPackageStartupMessages(library(ggplot2))
bd <- as.data.frame(boot[ok, ])
p <- ggplot(bd, aes(TM, RR)) +
  geom_point(alpha = 0.25, size = 1, colour = "#2a78d6") +
  geom_vline(xintercept = c(0, 1), colour = "grey70") + geom_hline(yintercept = c(0, 1), colour = "grey70") +
  annotate("point", x = hat[["TM"]], y = hat[["RR"]], colour = "#c0392b", size = 3) +
  labs(x = "xi_T", y = "xi_P", title = sprintf("Country-block bootstrap of (xi_T, xi_P), %d draws", sum(ok)),
       subtitle = "Red: the full-sample estimate; lines at 0 (equal dollars) and 1 (equal %)") +
  theme_minimal(base_size = 12)
png(file.path(out_dir, paste0("bootstrap_", MODE, ".png")), width = 700, height = 620); print(p); invisible(dev.off())
cat("Saved", file.path(out_dir, paste0("bootstrap_", MODE, ".{rds,png}")), "\n")
