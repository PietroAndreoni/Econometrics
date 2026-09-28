# The binned BHM climate interaction with the base climate block only, on the
# national PWT panel (with DOSE alongside for comparison). See
# functions/climate_slope_bins.R for the method.
#
# No temperature or precipitation deviation term enters any model, so the
# temperature slope is not competing with a separate anomaly response. On PWT
# GID_1 equals GID_0, so unit and country clustering coincide by construction.
#
# Bin schemes: equal-count bins put the same number of observations in each bin
# and compress where units pile up (many tropical countries sit within a degree
# or two of each other); equal-width bins keep the degree interval constant.
# Equal-width tails are very thin, and a cluster-robust Wald test over many
# restrictions with few clusters over-rejects, so a variant pools each tail
# inward until it holds at least THIN_UNITS units (pool_tails()).
#
# Specification: main_spec() (functions/analysis_spec.R), i.e. main_analysis.Rmd.
# Output: results/bhm_climate_slope_national/.

source("load_functions.R")

SPEC <- main_spec()
BIN_WIDTH <- 2          # Degrees, for the equal-width schemes.
THIN_UNITS <- 10L       # Bins with fewer units are flagged, and pooled below.
DATASETS <- c(`PWT national (GADM0)` = "PWT", `DOSE subnational (GADM1)` = "DOSE")
OUT <- file.path("results", "bhm_climate_slope_national")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

bin_schemes <- function(d) {
  x <- d$mean_TM_all
  width_breaks <- seq(floor(min(x) / BIN_WIDTH) * BIN_WIDTH,
                      ceiling(max(x) / BIN_WIDTH) * BIN_WIDTH, by = BIN_WIDTH)
  list(
    `Equal count, 5 bins` = stats::quantile(x, seq(0, 1, length.out = 6)),
    `Equal count, 10 bins` = stats::quantile(x, seq(0, 1, length.out = 11)),
    `Equal width, 2 degrees` = width_breaks,
    `Equal width, 2 degrees, tails pooled` =
      pool_tails(x, width_breaks, d$GID_1, THIN_UNITS)
  )
}

# Panels: only the base block is estimated, so it needs only its own variables.
needed <- c(SPEC$outcome, "TM", "mean_TM_all", "RR", "mean_RR_all")
panels <- lapply(DATASETS, function(econ) {
  panel <- build_main_dat(econ, SPEC) %>% filter(if_all(all_of(needed), is.finite))
  stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
  cat(sprintf("%s: %d unit-years, %d units, %d countries, %d-%d\n", econ,
              nrow(panel), n_distinct(panel$GID_1), n_distinct(panel$GID_0),
              min(panel$year), max(panel$year)))
  panel
})

# Estimation ---------------------------------------------------------------------

SCHEMES <- names(bin_schemes(panels[[1]]))
runs <- list()
for (label in names(DATASETS)) {
  schemes <- bin_schemes(panels[[label]])
  for (scheme in SCHEMES) {
    runs[[length(runs) + 1L]] <- climate_slope_analysis(
      assign_climate_bins(panels[[label]], schemes[[scheme]]),
      spec = SPEC,
      thin_units = THIN_UNITS,
      labels = list(dataset = label, scheme = scheme)
    )
  }
}
gather <- function(name) bind_rows(lapply(runs, `[[`, name))

bin_summary <- gather("bin_summary")
bin_slopes <- gather("bin_slopes")
linear_slopes <- gather("linear_slopes")
tests <- gather("tests")
model_fit <- gather("model_fit")

write_out(bin_summary, "climate_bin_definitions")
write_out(bin_slopes, "climate_bin_slopes")
write_out(tests, "linearity_tests")
write_out(model_fit, "model_fit_comparison")
write_out(gather("coefficients"), "coefficients")

# Figure ---------------------------------------------------------------------------

as_facets <- function(x) {
  mutate(x, dataset = factor(dataset, levels = names(DATASETS)),
         scheme = factor(scheme, levels = SCHEMES))
}
p <- plot_climate_slope_bins(
  as_facets(bin_slopes),
  as_facets(linear_slopes),
  as_facets(climate_slope_test_labels(tests, c("dataset", "scheme"))),
  facets = scheme ~ dataset,
  scales = "free_y",
  thin_units = THIN_UNITS,
  title = "Base climate block only: temperature slope by long-run mean temperature",
  subtitle = paste(
    "No deviation term in any model. Horizontal bars span each bin's range of",
    "mean_TM_all, so narrow bars mark compressed equal-count bins.\nHollow",
    "points mark bins with too few units to interpret."
  )
)
ggsave(file.path(OUT, "climate_slope_bins_national.png"), p,
       width = 12, height = 14, dpi = 200)

# Console summary ------------------------------------------------------------------

options(width = 180)
cat("\n=== Tests ===\n")
print(tests %>% select(dataset, scheme, n_bins, test, clustering, F_stat, p_value),
      n = Inf)
cat("\n=== Fit: linear vs binned climate interaction ===\n")
print(model_fit, n = Inf)
cat("\n=== PWT national, equal-width 2-degree bins, tails pooled ===\n")
print(
  bin_slopes %>%
    filter(dataset == names(DATASETS)[1],
           scheme == "Equal width, 2 degrees, tails pooled") %>%
    transmute(bin = mean_TM_bin,
              range = sprintf("%.1f to %.1f", mean_TM_all_min, mean_TM_all_max),
              countries, obs, slope_pp = round(100 * estimate, 2),
              se_pp = round(100 * se, 2), p = signif(p, 3), thin),
  n = Inf
)
cat("\nWrote figures and tables to ", OUT, "\n", sep = "")
