# Signed deviation bins cut on the anomaly distribution instead of on fixed
# standard-deviation thresholds, national against subnational.
#
# The notebook's `TM_bin_signed` uses fixed cutpoints at +/-0.5, 1, 1.5, 2 and
# 2.5 within-region standard deviations. Because anomalies are measured against
# a 30-year trailing mean in a warming period, the resulting cells are wildly
# unbalanced: on the DOSE panel the coldest bin holds 71 observations while the
# hottest holds 1262. Every tail conclusion then rests on the 71.
#
# Cutting instead at quantiles of the anomaly gives each bin the same number of
# observations, so precision is comparable across the response. The negative
# side gets fewer bins because it holds less mass, but each one is estimable.
#
# Only deviation terms are estimated here: no TM, no RR, and no BHM climate
# interaction. With region fixed effects and region-specific trends the
# deviation bins carry the whole weather response.
#
# Prerequisite: the notebook's parquet inputs under data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
})

source("rmd_chunks.R")

out <- "results/deviation_bins_distribution"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

# 11 quantile bins match the 11 categories of the fixed-threshold scheme, so
# the two have the same number of free coefficients and are comparable on AIC.
N_QBINS <- 11L
PP <- 100

# ---------------------------------------------------------------- panels ------

nb <- load_notebook_env("test_functions.Rmd")

NEEDED <- c("dlgrp_pc_usd", "zTM", "zRR", "TM_bin_signed", "RR_bin_signed")

datasets <- list(
  `PWT national (GADM0)` = "PWT",
  `DOSE subnational (GADM1)` = "DOSE"
)

# Quantile bins are renumbered to consecutive integers so that the reference
# level can be named by index, and the bin containing a zero anomaly is used as
# that reference to match the fixed-threshold scheme's baseline.
add_quantile_bins <- function(d) {
  for (v in c("zTM", "zRR")) {
    breaks <- quantile(d[[v]], seq(0, 1, length.out = N_QBINS + 1L), na.rm = TRUE)
    breaks <- unique(breaks)
    index <- as.integer(cut(d[[v]], breaks = breaks, include.lowest = TRUE))
    d[[paste0(v, "_qbin")]] <- factor(index)
    attr(d, paste0(v, "_ref")) <- as.integer(
      cut(0, breaks = breaks, include.lowest = TRUE)
    )
  }
  d
}

panels <- lapply(names(datasets), function(label) {
  cat("Building", label, "...\n")
  d <- nb$build_dat(econ_data = datasets[[label]]) %>%
    filter(if_all(all_of(NEEDED), is.finite)) %>%
    add_quantile_bins()
  stopifnot(!anyDuplicated(d[c("GID_1", "year")]))
  cat(sprintf(
    "  %d unit-years, %d units, %d countries, %d-%d | zTM ref bin %d, zRR ref bin %d\n",
    nrow(d), n_distinct(d$GID_1), n_distinct(d$GID_0),
    min(d$year), max(d$year), attr(d, "zTM_ref"), attr(d, "zRR_ref")
  ))
  d
})
names(panels) <- names(datasets)

# ---------------------------------------------------------------- models ------

schemes <- c("Fixed SD thresholds", "Equal-count quantiles")

scheme_terms <- function(scheme, d) {
  if (scheme == "Fixed SD thresholds") {
    list(
      terms = "i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)",
      tm_var = "TM_bin_signed", rr_var = "RR_bin_signed",
      tm_ref = "0", rr_ref = "0"
    )
  } else {
    list(
      terms = sprintf(
        "i(zTM_qbin, ref = %d) + i(zRR_qbin, ref = %d)",
        attr(d, "zTM_ref"), attr(d, "zRR_ref")
      ),
      tm_var = "zTM_qbin", rr_var = "zRR_qbin",
      tm_ref = as.character(attr(d, "zTM_ref")),
      rr_ref = as.character(attr(d, "zRR_ref"))
    )
  }
}

# Coefficients for one bin set, plus the omitted category at zero, joined to the
# observed anomaly range each bin covers.
extract_bins <- function(model, d, bin_var, ref_level, variable, z_var) {
  ct <- as.data.frame(coeftable(model))
  ci <- confint(model)
  rn <- rownames(ct)
  sel <- startsWith(rn, paste0(bin_var, "::"))
  levels_out <- sub(paste0("^", bin_var, "::"), "", rn[sel])

  est <- data.frame(
    variable = variable, bin = c(levels_out, ref_level),
    estimate = c(ct[sel, 1], 0), se = c(ct[sel, 2], NA_real_),
    p = c(ct[sel, 4], NA_real_),
    ci_low = c(ci[sel, 1], NA_real_), ci_high = c(ci[sel, 2], NA_real_)
  )

  spans <- d %>%
    mutate(bin = as.character(.data[[bin_var]]), z = .data[[z_var]]) %>%
    group_by(bin) %>%
    summarise(obs = n(), units = n_distinct(GID_1),
              z_min = min(z), z_max = max(z), z_mean = mean(z),
              .groups = "drop")

  left_join(est, spans, by = "bin")
}

run <- function(label, scheme) {
  d <- panels[[label]]
  spec <- scheme_terms(scheme, d)
  model <- fixest::feols(
    as.formula(paste(nb$outcome, "~", spec$terms, "|", nb$fixed_effects)),
    data = d, panel.id = nb$pan_id, cluster = ~GID_1
  )
  stopifnot(!anyNA(coef(model)))

  bins <- bind_rows(
    extract_bins(model, d, spec$tm_var, spec$tm_ref, "Temperature", "zTM"),
    extract_bins(model, d, spec$rr_var, spec$rr_ref, "Precipitation", "zRR")
  ) %>%
    mutate(dataset = label, scheme = scheme, .before = 1)

  joint <- purrr::map_dfr(
    list(Temperature = spec$tm_var, Precipitation = spec$rr_var),
    function(v) {
      w <- wald(model, keep = paste0("^", v, "::"), print = FALSE)
      data.frame(F_stat = w$stat, df1 = w$df1, p = w$p)
    }, .id = "variable"
  ) %>%
    mutate(dataset = label, scheme = scheme, .before = 1)

  fit <- data.frame(
    dataset = label, scheme = scheme, n = nobs(model), k = length(coef(model)),
    within_r2 = as.numeric(fitstat(model, "wr2")[[1]]),
    aic = AIC(model), bic = BIC(model)
  )

  # Cell balance is the point of the exercise, so it is reported explicitly.
  balance <- bins %>%
    filter(!is.na(obs)) %>%
    group_by(dataset, scheme, variable) %>%
    summarise(min_obs = min(obs), max_obs = max(obs),
              imbalance = max(obs) / min(obs), .groups = "drop")

  list(bins = bins, joint = joint, fit = fit, balance = balance)
}

passes <- unlist(
  lapply(names(panels), function(label) {
    lapply(schemes, function(scheme) run(label, scheme))
  }),
  recursive = FALSE
)
gather <- function(name) bind_rows(lapply(passes, `[[`, name))

bins <- gather("bins")
joint_tests <- gather("joint")
model_fit <- gather("fit")
balance <- gather("balance")

# ------------------------------------------------------------------- plots ----

variable_levels <- c("Temperature", "Precipitation")
prep <- function(x) {
  mutate(
    x,
    dataset = factor(dataset, levels = names(datasets)),
    variable = factor(variable, levels = variable_levels),
    scheme = factor(scheme, levels = schemes)
  )
}

scheme_colors <- c(`Fixed SD thresholds` = "#0072B2",
                   `Equal-count quantiles` = "#D55E00")

p <- ggplot(prep(bins), aes(color = scheme)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.35) +
  # Horizontal bars span the anomaly range each bin covers, which is where the
  # two schemes differ most: the quantile tails are wide, the fixed ones narrow.
  geom_segment(
    aes(x = z_min, xend = z_max, y = estimate * PP, yend = estimate * PP),
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(x = z_mean, ymin = ci_low * PP, ymax = ci_high * PP),
    width = 0.1, linewidth = 0.4
  ) +
  geom_point(aes(x = z_mean, y = estimate * PP), size = 1.7) +
  scale_color_manual(values = scheme_colors, name = "Bin cutpoints") +
  facet_grid(variable ~ dataset, scales = "free_y") +
  labs(
    title = "Signed deviation bins: fixed thresholds against the anomaly distribution",
    subtitle = paste(
      "Deviation terms only; no TM, RR, or BHM climate interaction.",
      "Horizontal bars span the anomaly range each bin covers,\nso wide bars",
      "mark the pooled tails that equal-count binning produces.",
      "95% intervals clustered by GID_1."
    ),
    x = "Standardized anomaly, within-unit standard deviations",
    y = "Log GDP-per-capita growth (percentage points)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.9, "lines"))

ggsave(file.path(out, "deviation_bins_distribution.png"), p,
       width = 12, height = 8, dpi = 180)

# ------------------------------------------------------------------ output ----

write_out(bins, "deviation_bin_coefficients")
write_out(joint_tests, "joint_tests")
write_out(model_fit, "model_fit_comparison")
write_out(balance, "cell_balance")

cat("\n=== Cell balance (observations per bin) ===\n")
print(as.data.frame(balance), row.names = FALSE, digits = 4)
cat("\n=== Joint significance of each bin set ===\n")
print(as.data.frame(joint_tests), row.names = FALSE, digits = 4)
cat("\n=== Fit ===\n")
print(as.data.frame(model_fit), row.names = FALSE, digits = 5)
cat("\n=== Temperature bins ===\n")
print(
  bins %>%
    filter(variable == "Temperature") %>%
    transmute(dataset = sub(" .*", "", dataset), scheme, bin,
              z_range = sprintf("%6.2f to %5.2f", z_min, z_max),
              obs, coef_pp = round(estimate * PP, 2),
              se_pp = round(se * PP, 2), p = signif(p, 2)) %>%
    arrange(dataset, scheme, z_range) %>%
    as.data.frame(),
  row.names = FALSE
)
cat("\nWrote figures and tables to ", out, "\n", sep = "")
