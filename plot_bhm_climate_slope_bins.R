# Replace the linear BHM climate interaction with a free temperature slope per
# bin of long-run mean temperature, on the DOSE subnational panel, and test the
# linear restriction against it (see functions/climate_slope_bins.R).
#
# Because the within-region variation in TM *is* the anomaly, the estimated
# climate profile depends on the deviation term that is also in the model. Four
# are reported, including the fully nonparametric signed-deviation bins, which
# leave the deviation response unrestricted. Equal-count bin counts of 5 and 10
# are both reported so that features of the profile can be checked against the
# resolution used to find them.
#
# Specification: main_spec() (functions/analysis_spec.R), i.e. main_analysis.Rmd.
# Output: results/bhm_climate_slope_bins/.

source("load_functions.R")

SPEC <- main_spec()
N_BINS_SET <- c(5L, 10L)
OUT <- file.path("results", "bhm_climate_slope_bins")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# Deviation terms; each controls both its temperature and its precipitation
# deviation. hTM_10 is the best-fitting linear hinge, hTM_15 the notebook's zzTM.
DEVIATION_SPECS <- c(
  "No deviation term" = "",
  "Linear hinge at 1.0 SD" = "hTM_10 + abs_zRRp",
  "Quadratic hinge at 1.5 SD (notebook)" = "hTM_15^2 + abs_zRRp",
  "Signed deviation bins" = "i(TM_bin_signed, ref = 0) + i(RR_bin_signed, ref = 0)"
)

# Panel: one estimation sample for every specification, bin variables included,
# so the comparison is about functional form rather than coverage.
panel <- build_main_dat("DOSE", SPEC) %>%
  add_hinge_columns(c(1.0, 1.5), sided = FALSE)
needed <- c(SPEC$outcome, "TM", "mean_TM_all", "RR", "mean_RR_all", "zTM",
            "hTM_10", "hTM_15", "abs_zRRp", "TM_bin_signed", "RR_bin_signed")
est <- panel %>% filter(if_all(all_of(needed), is.finite))
stopifnot(!anyDuplicated(est[c("GID_1", "year")]))
cat(sprintf("Estimation sample: %d region-years, %d regions, %d countries, %d-%d\n",
            nrow(est), n_distinct(est$GID_1), n_distinct(est$GID_0),
            min(est$year), max(est$year)))

# Estimation ---------------------------------------------------------------------

runs <- list()
for (n_bins in N_BINS_SET) {
  breaks <- stats::quantile(est$mean_TM_all, seq(0, 1, length.out = n_bins + 1L))
  binned <- assign_climate_bins(est, breaks)
  for (deviation in names(DEVIATION_SPECS)) {
    runs[[length(runs) + 1L]] <- climate_slope_analysis(
      binned,
      deviation_terms = DEVIATION_SPECS[[deviation]],
      spec = SPEC,
      labels = list(deviation = deviation)
    )
  }
}
gather <- function(name) bind_rows(lapply(runs, `[[`, name))

bin_summary <- gather("bin_summary") %>% distinct(n_bins, mean_TM_bin, .keep_all = TRUE) %>%
  select(-deviation)
bin_slopes <- gather("bin_slopes")
linear_slopes <- gather("linear_slopes")
tests <- gather("tests")
model_fit <- gather("model_fit") %>% arrange(n_bins, deviation, desc(within_r2))
coefficients <- gather("coefficients")

write_out(bin_summary, "climate_bin_definitions")
write_out(bin_slopes, "climate_bin_slopes")
write_out(tests, "linearity_tests")
write_out(model_fit, "model_fit_comparison")
write_out(coefficients, "coefficients")

# Figure ---------------------------------------------------------------------------

as_facets <- function(x) {
  mutate(x,
         deviation = factor(deviation, levels = names(DEVIATION_SPECS)),
         bin_label = factor(sprintf("%d climate bins", n_bins),
                            levels = sprintf("%d climate bins", N_BINS_SET)))
}
p <- plot_climate_slope_bins(
  as_facets(bin_slopes),
  as_facets(linear_slopes),
  as_facets(climate_slope_test_labels(tests, c("n_bins", "deviation"))),
  facets = bin_label ~ deviation,
  title = "Marginal effect of a temperature deviation, by long-run mean temperature (DOSE)",
  subtitle = paste(
    "The BHM block restricts this slope to be linear in mean_TM_all; the bins let",
    "it vary freely. Each column holds a different deviation term.\nHorizontal bars",
    "span each bin's range of mean_TM_all; 95% intervals clustered by region."
  )
)
ggsave(file.path(OUT, "climate_slope_bins.png"), p, width = 15, height = 8.5, dpi = 200)

# Console summary ------------------------------------------------------------------

options(width = 180)
cat("\n=== Tests ===\n")
print(tests %>% select(n_bins, deviation, test, clustering, F_stat, p_value), n = Inf)
cat("\n=== Fit: linear vs binned climate interaction ===\n")
print(model_fit, n = Inf)
cat("\n=== Free temperature slope by climate bin, 10 bins (pp per degree) ===\n")
print(
  bin_slopes %>%
    filter(n_bins == 10) %>%
    transmute(deviation, bin = mean_TM_bin,
              range = sprintf("%.1f to %.1f", mean_TM_all_min, mean_TM_all_max),
              units, slope_pp = round(100 * estimate, 3), se_pp = round(100 * se, 3),
              p = signif(p, 3), p_country = signif(p_country, 3)),
  n = Inf
)
cat("\nWrote figures and tables to ", OUT, "\n", sep = "")
