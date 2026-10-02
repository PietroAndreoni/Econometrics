# The Burke-Hsiang-Miguel (2015) quadratic on the national PWT panel, drawn in
# temperature-growth space together with the country-by-country pieces that
# identify it.
#
# (a) Pooled BHM regression: growth on TM + TM^2 + RR + RR^2 with country fixed
#     effects, country-specific quadratic trends and year fixed effects,
#     clustered by country. The curve is f(T) - f(T*), zero at the estimated
#     optimum T*, with 95% delta-method intervals (T* held at its point
#     estimate).
# (b) Country pieces: country_pieces() in functions/country_pieces.R. Each
#     piece is a segment of the country's own OLS slope, centred on the curve
#     at its mean temperature and spanning +/- 2 SD of its identifying
#     variation; opacity and width scale with its weight in the pooled linear
#     estimate, which the weighted pieces rebuild exactly.
#
# The lower panel repeats the pieces in marginal-effect space: b_i against mean
# temperature, with the quadratic's implied slope b1 + 2 b2 T.
#
# The same figure for the specification with deviation terms is in
# plot_bhm_dev_country_pieces.R.
#
# Specification: main_spec() (functions/analysis_spec.R) fixed effects and
# cluster, on PWT from its first year.
# Output: results/bhm_country_pieces/.

source("load_functions.R")

SPEC <- main_spec(econ_year_min = 1900L)
CONTROLS <- "RR + RR^2"
OUT <- file.path("results", "bhm_country_pieces")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

needed <- c(SPEC$outcome, "TM", "RR")
panel <- build_main_dat("PWT", SPEC) %>%
  filter(if_all(all_of(needed), is.finite)) %>%
  add_climate_zone()
stopifnot(!anyDuplicated(panel[c("GID_1", "year")]))
cat(sprintf("PWT: %d country-years, %d countries, %d-%d\n", nrow(panel),
            n_distinct(panel$GID_1), min(panel$year), max(panel$year)))

res <- country_pieces(panel, CONTROLS, SPEC)
write_out(res$pieces, "country_pieces")
write_out(res$curve, "bhm_curve")

p <- plot_country_pieces(
  res,
  title = "BHM quadratic and the country pieces that identify it (PWT)",
  spec_label = "Specification: TM + TM² + RR + RR²"
)
ggsave(file.path(OUT, "bhm_country_pieces.png"), p,
       width = 10, height = 9, dpi = 200)

etable(res$m_quad, res$m_lin, fitstat = ~n + r2 + wr2)
cat(sprintf("Optimum: %.2f C; pooled linear slope from pieces: %.5f\n",
            res$t_opt, res$pooled_from_pieces))
print(res$pieces %>% arrange(desc(weight)) %>%
        transmute(unit, climate_zone, obs, TM_mean = round(TM_mean, 1),
                  slope_pp = round(100 * slope, 2),
                  curve_slope_pp = round(100 * curve_slope, 2),
                  gap_z = round(gap_z, 1), weight = round(weight, 4)),
      n = 15)
cat("\nWrote figure and tables to ", OUT, "\n", sep = "")
