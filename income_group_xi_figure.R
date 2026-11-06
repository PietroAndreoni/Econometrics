# Figure of income_group_xi.R: the warm arm's response at the median climate by
# income group (95% CI; the horizontal bar spans the group's incomes), in
# growth and in dollars per person, with the fitted power laws.
#
#   Rscript income_group_xi_figure.R [tag]     (default PWT_trailing_J5)

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args)) args[1] else "PWT_trailing_J5"
dir <- file.path("results", "income_group_xi", tag)
e <- readRDS(file.path(dir, "estimates.rds"))

GROWTH <- "Growth response (pp per +1 SD)"
DOLLAR <- "Dollar loss per person (PPP $ per +1 SD)"
g <- e$groups
# Dollar loss per person per +1 SD: -Y_j b_j (b in pp).
pts <- rbind(
  data.frame(panel = GROWTH, income = g$income, lo_x = g$income_lo, hi_x = g$income_hi,
             est = g$warm, lo = g$warm - 1.96 * g$warm_se, hi = g$warm + 1.96 * g$warm_se),
  data.frame(panel = DOLLAR, income = g$income, lo_x = g$income_lo, hi_x = g$income_hi,
             est = -g$income * g$warm / 100, lo = -g$income * (g$warm + 1.96 * g$warm_se) / 100,
             hi = -g$income * (g$warm - 1.96 * g$warm_se) / 100))
cv <- rbind(transform(e$curves, panel = GROWTH, est = warm),
            transform(e$curves, panel = DOLLAR, est = -income * warm / 100))
# The curves at extreme xi explode outside the groups: cut them at each
# panel's range of the CIs, padded.
lim <- lapply(split(pts, pts$panel), function(x) range(x$lo, x$hi) + c(-0.1, 0.1) * diff(range(x$lo, x$hi)))
cv$est[cv$est < vapply(lim, `[`, 0, 1)[cv$panel] | cv$est > vapply(lim, `[`, 0, 2)[cv$panel]] <- NA
cv$curve <- factor(cv$curve, levels = unique(e$curves$curve))
pts$panel <- factor(pts$panel, c(GROWTH, DOLLAR))
cv$panel <- factor(cv$panel, c(GROWTH, DOLLAR))

ink <- "#3d3d3a"
p <- ggplot() +
  geom_hline(yintercept = 0, colour = "#b4b2a9", linewidth = 0.4) +
  geom_line(data = cv, aes(income, est, colour = curve, linetype = curve), linewidth = 0.7, na.rm = TRUE) +
  geom_segment(data = pts, aes(x = lo_x, xend = hi_x, y = est, yend = est), colour = "#b4b2a9", linewidth = 0.4) +
  geom_errorbar(data = pts, aes(income, ymin = lo, ymax = hi), width = 0, colour = ink, linewidth = 0.5) +
  geom_point(data = pts, aes(income, est), colour = ink, size = 2.2) +
  facet_wrap(~panel, scales = "free_y") +
  scale_x_log10(labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "longdash")) +
  labs(x = "Trailing GDP per capita (point: group geometric mean; grey bar: group range)", y = NULL,
       colour = NULL, linetype = NULL,
       title = sprintf("Warm arm by income group (%d groups, %s), at the median climate", e$J,
                       if (e$by == "trailing") "by trailing income" else "by country mean income"),
       subtitle = "Points: unrestricted first-stage estimates, 95% CI (clustered by country). Lines: b = beta (Y/Ybar)^(xi - 1).") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "#6b6a63"), strip.text = element_text(face = "bold", hjust = 0))
ggsave(file.path(dir, "warm_by_group.png"), p, width = 9, height = 4.6, dpi = 200, bg = "white")
