# The nonlinear least-squares estimate of xi (estimates.rds$xi_nls of
# handout_preferred_specification.R): its inference table and figure, shared by
# handout_preferred_specification.Rmd and the other GDP sources.
#
#   Rscript handout_xi_figure.R WB   # writes xi_nls.png and xi_nls.csv in
#                                    # results/handout_preferred_specification/WB

# The pieces of the inference of `e` (estimates.rds).
xi_helpers <- function(e) {
  xn <- e$xi_nls
  G <- e$sample[["countries"]]
  prof <- xn$profile
  draws <- xn$boot_draws[is.finite(xn$boot_draws)]
  list(
    xn = xn, G = G, prof = prof, draws = draws,
    at = function(x, col) prof[[col]][abs(prof$xi - x) < 1e-9],
    wald_p = function(x) 2 * pt(-abs((xn$estimate - x) / xn$wald[["se"]]), G - 1),
    boot_p = function(x) min(1, 2 * min(mean(draws <= x), mean(draws >= x))),
    fmt_set = function(s) paste(apply(s, 1, function(r) {
      if (!is.finite(r[["lo"]])) sprintf("ξ ≤ %.2f", r[["hi"]])
      else if (!is.finite(r[["hi"]])) sprintf("ξ ≥ %.2f", r[["lo"]])
      else sprintf("[%.2f, %.2f]", r[["lo"]], r[["hi"]])
    }), collapse = " ∪ "),
    at_low = mean(draws <= xn$bounds[1])
  )
}

# 95% confidence sets and p-values of xi = 0, 0.6, 1 by each inference;
# `fmt_p` formats the p-values.
xi_table <- function(e, fmt_p = function(p) sprintf("%.3f", p)) {
  h <- xi_helpers(e)
  xn <- h$xn
  p_row <- function(x) fmt_p(c(h$wald_p(x), h$at(x, "p_ar"), h$at(x, "p_score"), h$boot_p(x)))
  data.frame(
    Inference = c("Wald, full-model clustered SE", "Robust test, inverted", "Score test, inverted",
                  "Country bootstrap, percentile"),
    `95% confidence set` = c(sprintf("[%.2f, %.2f]", xn$wald[["lo"]], xn$wald[["hi"]]), h$fmt_set(xn$ar_set),
                             h$fmt_set(xn$score_set), sprintf("[%.2f, %.2f]", xn$boot[["lo"]], xn$boot[["hi"]])),
    `p: ξ = 0` = p_row(0), `p: ξ = 0.6` = p_row(0.6), `p: ξ = 1` = p_row(1), check.names = FALSE)
}

# Profile (with the Wald CI), p-values of the tests (with the robust-test set)
# and bootstrap estimates (with the percentile CI), against xi.
xi_figure <- function(e) {
  h <- xi_helpers(e)
  xn <- h$xn
  prof <- h$prof
  panels <- c(sprintf("Profile: AIC − minimum. Shaded: Wald 95%% CI [%.2f, %.2f]", xn$wald[["lo"]], xn$wald[["hi"]]),
              sprintf("p-value of the test of ξ. Shaded: robust-test 95%% set %s", h$fmt_set(xn$ar_set)),
              sprintf("Bootstrap estimates (share of draws). Shaded: percentile 95%% CI [%.2f, %.2f]",
                      xn$boot[["lo"]], xn$boot[["hi"]]))
  clip <- function(x) pmin(pmax(x, xn$bounds[1]), xn$bounds[2])
  sets_df <- rbind(
    data.frame(panel = panels[1], lo = clip(xn$wald[["lo"]]), hi = clip(xn$wald[["hi"]])),
    data.frame(panel = panels[2], lo = clip(xn$ar_set$lo), hi = clip(xn$ar_set$hi)),
    data.frame(panel = panels[3], lo = xn$boot[["lo"]], hi = xn$boot[["hi"]]))
  sets_df$panel <- factor(sets_df$panel, panels)
  lines_df <- rbind(
    data.frame(panel = panels[1], xi = prof$xi, value = prof$aic - min(prof$aic), series = "Profile"),
    data.frame(panel = panels[2], xi = prof$xi, value = prof$p_ar, series = "Robust test"),
    data.frame(panel = panels[2], xi = prof$xi, value = prof$p_score, series = "Score test"))
  lines_df$panel <- factor(lines_df$panel, panels)
  hist_df <- data.frame(panel = factor(panels[3], panels), xi = h$draws)
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = sets_df, ggplot2::aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
                       fill = "#2a78d6", alpha = 0.12) +
    ggplot2::geom_vline(xintercept = c(0, 0.6, 1), colour = "grey75", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = xn$estimate, colour = "grey20", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_hline(data = data.frame(panel = factor(panels[2], panels), y = 0.05), ggplot2::aes(yintercept = y),
                        colour = "grey40", linewidth = 0.3, linetype = "dotted") +
    ggplot2::geom_histogram(data = hist_df, ggplot2::aes(xi, ggplot2::after_stat(count / sum(count))),
                            binwidth = 0.1, boundary = 0, fill = "grey55", colour = "white", linewidth = 0.2) +
    ggplot2::geom_line(data = lines_df, ggplot2::aes(xi, value, colour = series, linetype = series), linewidth = 0.7) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(Profile = "grey20", `Robust test` = "#2a78d6", `Score test` = "#eb6834"),
                                 breaks = c("Robust test", "Score test"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(Profile = "solid", `Robust test` = "solid", `Score test` = "longdash"),
                                   breaks = c("Robust test", "Score test"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = -4:4) +
    ggplot2::labs(x = "ξ (income elasticity of the dollar loss)", y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"))
}

# Several runs (a named list of estimates.rds): rows and countries, the
# estimate, the 95% sets and the p-values of xi = 0, 0.6, 1 (Wald / robust /
# bootstrap), one row per run.
xi_compare_table <- function(runs, fmt_p = function(p) sprintf("%.3f", p)) {
  do.call(rbind, lapply(names(runs), function(k) {
    e <- runs[[k]]
    h <- xi_helpers(e)
    xn <- h$xn
    p3 <- function(x) paste(fmt_p(c(h$wald_p(x), h$at(x, "p_ar"), h$boot_p(x))), collapse = " / ")
    data.frame(Run = k, `Rows (countries)` = sprintf("%s (%d)", format(e$sample[["rows"]], big.mark = ","),
                                                     e$sample[["countries"]]),
               `ξ̂` = sprintf("%.2f", xn$estimate),
               `Wald 95% CI` = sprintf("[%.2f, %.2f]", xn$wald[["lo"]], xn$wald[["hi"]]),
               `Robust-test 95% set` = h$fmt_set(xn$ar_set),
               `Bootstrap 95% CI` = sprintf("[%.2f, %.2f]", xn$boot[["lo"]], xn$boot[["hi"]]),
               `p: ξ = 0` = p3(0), `p: ξ = 0.6` = p3(0.6), `p: ξ = 1` = p3(1), check.names = FALSE)
  }))
}

# Two runs on the same rows with the same bootstrap draws (`a`, `b`: estimates;
# `labels`: their names): profiles, bootstrap estimates, and the paired
# bootstrap difference b - a.
xi_compare_figure <- function(a, b, labels) {
  panels <- c("Profile: AIC − minimum", "Bootstrap estimates (share of draws, bins of 0.1)",
              sprintf("Paired bootstrap difference: %s − %s", labels[2], labels[1]))
  run_df <- function(e, k) {
    p <- e$xi_nls$profile
    rbind(data.frame(panel = panels[1], xi = p$xi, value = p$aic - min(p$aic), run = k),
          {
            br <- seq(-4.05, 4.05, 0.1)
            x <- e$xi_nls$boot_draws[is.finite(e$xi_nls$boot_draws)]
            data.frame(panel = panels[2], xi = head(br, -1) + 0.05,
                       value = as.numeric(table(cut(x, br))) / length(x), run = k)
          })
  }
  lines_df <- rbind(run_df(a, labels[1]), run_df(b, labels[2]))
  lines_df$panel <- factor(lines_df$panel, panels)
  lines_df$run <- factor(lines_df$run, labels)
  diff_df <- data.frame(panel = factor(panels[3], panels), xi = b$xi_nls$boot_draws - a$xi_nls$boot_draws)
  est_df <- data.frame(panel = factor(rep(panels[1:2], each = 2), panels),
                       xi = rep(c(a$xi_nls$estimate, b$xi_nls$estimate), 2), run = factor(rep(labels, 2), labels))
  ggplot2::ggplot() +
    ggplot2::geom_vline(data = data.frame(panel = factor(panels[3], panels), x = 0), ggplot2::aes(xintercept = x),
                        colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_histogram(data = diff_df, ggplot2::aes(xi, ggplot2::after_stat(count / sum(count))),
                            binwidth = 0.2, boundary = 0, fill = "grey55", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(data = est_df, ggplot2::aes(xintercept = xi, colour = run), linewidth = 0.4,
                        linetype = "dashed", show.legend = FALSE) +
    ggplot2::geom_line(data = lines_df, ggplot2::aes(xi, value, colour = run, linetype = run), linewidth = 0.7) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free") +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#2a78d6", "#eb6834"), labels), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "longdash"), labels), name = NULL) +
    ggplot2::labs(x = "ξ (bottom: difference in ξ)", y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"))
}

if (sys.nframe() == 0) {
  gdp <- commandArgs(TRUE)[1]
  dir <- file.path("results", "handout_preferred_specification")
  if (!is.na(gdp) && gdp != "PWT") dir <- file.path(dir, gdp)
  e <- readRDS(file.path(dir, "estimates.rds"))
  h <- xi_helpers(e)
  tab <- xi_table(e)
  utils::write.csv(tab, file.path(dir, "xi_nls.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  ggplot2::ggsave(file.path(dir, "xi_nls.png"), xi_figure(e), width = 6.5, height = 5.4, dpi = 200, bg = "white")
  cat(sprintf("xi hat = %.3f, Wald SE %.3f, bootstrap SE %.3f; %.0f%% of draws at a bound (%.0f%% at the lower)\n",
              h$xn$estimate, h$xn$wald[["se"]], h$xn$boot[["se"]], 100 * h$xn$boot[["at_bound"]], 100 * h$at_low))
  print(tab, right = FALSE)
}
