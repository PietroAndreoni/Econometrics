# Single-country version of A6b: every year 1970-2023, absolute exposure
# deviation from the country mean against the deviation of the annual mean in
# the country's SD (interannual, 1970-2023). Top row T (hot DD > 25 C, cold
# DD < 5 C, sum), bottom row P (wet = mm above 20 mm/day, dry = days < 1 mm,
# sum = each over its cross-country SD of the climatologies, as in
# bhm_precip_exposure.R). Dashed: linear fit per country.
#
#   Rscript bhm_country_exposure.R            # India and the USA
#   Rscript bhm_country_exposure.R BRA DEU    # any ISO3 codes

suppressPackageStartupMessages(library(ggplot2))
source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "bhm_microfoundation")
args <- commandArgs(trailingOnly = TRUE)
SEL <- if (length(args)) toupper(args) else c("IND", "USA")
r <- readRDS(file.path(out_dir, "results.rds"))
all_c <- r$country$GID_0

daily <- function(v) {
  x <- as.data.frame(arrow::read_parquet(wcd_download(v, "ERA5", "gadm0", "pop", 2015, "daily", verbose = FALSE)))
  yr <- as.integer(substr(sub("^X", "", x$Date), 1, 4))
  s <- yr >= 1970 & yr <= 2023
  m <- as.matrix(x[s, intersect(all_c, names(x))]); storage.mode(m) <- "double"
  m <- m[, colSums(is.na(m)) == 0, drop = FALSE]
  list(m = m, yr = yr[s])
}
tm <- daily("tmp"); pr <- daily("pre")
stopifnot(all(SEL %in% colnames(tm$m)), all(SEL %in% colnames(pr$m)))
# Cross-country SDs of the P climatologies (the weights of the P sum).
wet_all <- rowsum(pmax(pr$m - 20, 0), pr$yr); dry_all <- rowsum((pr$m < 1) * 1, pr$yr)
sdw <- sd(colMeans(wet_all)); sdd <- sd(colMeans(dry_all))

rows <- list()
for (i in SEL) {
  t_d <- tm$m[, i]; p_d <- pr$m[, i]
  Ta <- tapply(t_d, tm$yr, mean); Pa <- tapply(p_d, pr$yr, sum)
  hot <- tapply(pmax(t_d - 25, 0), tm$yr, sum); cold <- tapply(pmax(5 - t_d, 0), tm$yr, sum)
  wet <- tapply(pmax(p_d - 20, 0), pr$yr, sum); dry <- tapply((p_d < 1) * 1, pr$yr, sum)
  zT <- (Ta - mean(Ta)) / sd(Ta); zP <- (Pa - mean(Pa)) / sd(Pa)
  idx <- list(
    "T hot (DD > 25 C)" = list(zT, hot), "T cold (DD < 5 C)" = list(zT, cold), "T sum (DD)" = list(zT, hot + cold),
    "P wet (mm > 20 mm/day)" = list(zP, wet), "P dry (days < 1 mm)" = list(zP, dry),
    "P sum (cross-country SD units)" = list(zP, wet / sdw + dry / sdd))
  for (k in names(idx)) {
    y <- idx[[k]][[2]]
    rows[[paste(i, k)]] <- data.frame(GID_0 = i, panel = k, year = as.integer(names(y)),
                                      z = as.numeric(idx[[k]][[1]]), y = as.numeric(y - mean(y)))
  }
  cat(sprintf("%s: mean T %.1f C (SD %.2f), mean P %.0f mm (SD %.0f); mean hot DD %.0f, cold DD %.0f, wet %.0f mm, dry %.0f days\n",
              i, mean(Ta), sd(Ta), mean(Pa), sd(Pa), mean(hot), mean(cold), mean(wet), mean(dry)))
}
d <- do.call(rbind, rows)
d$panel <- factor(d$panel, unique(d$panel))
p <- ggplot(d, aes(z, y, colour = GID_0)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 1.4, alpha = .8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = 2, linewidth = .6) +
  facet_wrap(~panel, nrow = 2, scales = "free_y") +
  labs(x = "within-country deviation of the annual mean (country SD; T top row, P bottom row)",
       y = "exposure deviation from the country mean (index units)", colour = NULL,
       title = paste("Within-country exposure, one point per year 1970-2023:", paste(SEL, collapse = " vs ")),
       subtitle = "ERA5 daily, population weights 2015; dashed: linear fit") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
f <- file.path(out_dir, paste0("A7_country_exposure_", paste(SEL, collapse = "_"), ".png"))
ggsave(f, p, width = 10, height = 6.5, dpi = 150)
cat("Saved", f, "\n")
