# Private or national adaptive capacity: the preferred specification's weather
# arms (handout_preferred_specification.R) estimated on the growth of each
# income decile, with the arms interacted with the decile's own income, the
# country's income, or both. Saved to results/decile_adaptation/estimates.rds.
#
# Deciles: net income shares of Narayan et al. (2023), as assembled by Gilli et
# al. (2024, JEEM) (data/data_gdp_gilli.parquet, 1958-2015). Decile income per
# capita = PWT GDP per capita x share / 10, so decile growth = the panel's
# country growth + the growth of the decile's share. Moderators: country income
# = the specification's log mean GDP per capita of the 5 years before; decile
# income = the log mean decile income of the same 5 years. Same rows for every
# decile: the specification's rows with shares in the year, the year before and
# the 5 years before (so up to 2015).
#
# Per decile, three models, all controlling for the levels of the trailing
# climate (T, P) and of both incomes, so that they differ only in the weather
# terms: arms x (climate, country income), arms x (climate, decile income),
# arms x (climate, both); P arms also x agricultural share. Then all deciles
# stacked, with decile-specific year effects and decile-country trends, and
# with country-year effects (weather, country income and climate absorbed: only
# the decile income interactions are identified, from the deciles of the same
# country-year).
#
#   Rscript decile_adaptation.R

# The handout's panel, deviations and test helpers (up to its models).
src <- readLines("handout_preferred_specification.R")
eval(parse(text = src[1:(grep("^models <- list\\(", src) - 1)]))
source("load_functions.R")
# As in the handout: fit verdicts against noise (clustered Vuong, 5%) and the
# coefficients without the controls.
verdict <- function(a, b) fit_verdict(fit_record(a), fit_record(b), "aic")
coef_of <- function(m) {
  ct <- coeftable(m)
  ct <- ct[!grepl("conflict|collapse", rownames(ct)), , drop = FALSE]
  data.frame(term = gsub("_X_", ":", rownames(ct)), estimate = ct[, 1], se = ct[, 2], p = ct[, 4], row.names = NULL)
}

DECILES <- paste0("D", 1:10)
WINDOW <- P$trailing_income_years

# PWT GDP per capita over its whole history (as the panel's trailing income).
gdp <- load_econ_panel(resolve_econ_source("PWT")) %>%
  filter(gadm_level == GADM_LEVEL_BY_SOURCE[[resolve_econ_source("PWT")]], is.finite(grp_pc_usd)) %>%
  distinct(GID_1, year, grp_pc_usd)
shares <- as.data.frame(arrow::read_parquet(file.path("data", "data_gdp_gilli.parquet")))[, c("iso3", "year", DECILES)]
trailing <- function(x, year) log(lag_by_year(zoo::rollapplyr(x, WINDOW, mean, fill = NA_real_), year))
hist <- gdp %>%
  inner_join(shares, by = c(GID_1 = "iso3", "year")) %>%
  group_by(GID_1) %>%
  tidyr::complete(year = seq(min(year), max(year))) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(inc_check = trailing(grp_pc_usd, year),
         across(all_of(DECILES), list(dls = ~ log(.x) - log(lag_by_year(.x, year)),
                                      yq = ~ trailing(grp_pc_usd * .x / 10, year)))) %>%
  ungroup()

dd <- d %>% inner_join(select(hist, GID_1, year, inc_check, ends_with("_dls"), ends_with("_yq")), by = c("GID_1", "year"))
dd <- dd[stats::complete.cases(dd[, c(paste0(DECILES, "_dls"), paste0(DECILES, "_yq"))]), ]
# The trailing country income is the panel's.
stopifnot(max(abs(dd$inc_check - dd$inc_trailing)) < 1e-8)
dd$y <- dd$inc_trailing - mean(dd$inc_trailing)

# A model from a list of terms (name = product of columns), with the levels of
# the trailing climate and of both incomes, the controls and fixed effects `fe`.
fit_with <- function(terms, data, fe = FE, levels = c("clim_TM", "clim_RR", "y", "yq")) {
  for (k in names(terms)) data[[k]] <- Reduce(`*`, lapply(terms[[k]], function(z) data[[z]]))
  f <- as.formula(paste("g ~", paste(c(names(terms), levels, CONTROLS), collapse = " + "), "|", fe))
  m <- feols(f, data = data, cluster = ~GID_0, notes = FALSE)
  m$country <- data$GID_0[obs(m)]
  m
}
spec_terms <- function(income) {
  c(arm_terms("TM", c(CLIM[["TM"]], income)), arm_terms("RR", c(CLIM[["RR"]], "agr", income)))
}
SPECS <- list(country = spec_terms("y"), decile = spec_terms("yq"), both = spec_terms(c("y", "yq")))
of <- function(v, z) paste0(c("p", "n"), v, "_X_", z)
income_tests <- function(m) {
  rbind(`T: decile income` = zero_test(m, of("TM", "yq")), `T: country income` = zero_test(m, of("TM", "y")),
        `P: decile income` = zero_test(m, of("RR", "yq")), `P: country income` = zero_test(m, of("RR", "y")))
}
income_coefs <- function(m) {
  ct <- coeftable(m)
  ct <- ct[grepl("_X_(y|yq)$", rownames(ct)), , drop = FALSE]
  data.frame(term = gsub("_X_", ":", rownames(ct)), estimate = ct[, 1], se = ct[, 2], p = ct[, 4], row.names = NULL)
}

decile_data <- function(q) {
  x <- dd
  x$g <- x$g + x[[paste0(q, "_dls")]]
  x$yq <- x[[paste0(q, "_yq")]] - mean(x[[paste0(q, "_yq")]])
  x
}
per_decile <- lapply(stats::setNames(DECILES, DECILES), function(q) {
  x <- decile_data(q)
  ms <- lapply(SPECS, fit_with, data = x)
  list(models = ms, tests = income_tests(ms$both), coefs = income_coefs(ms$both),
       verdict = verdict(ms$decile, ms$country),
       aic = vapply(ms, AIC, 0), cor_income = cor(x$y, x$yq))
})
# The country's growth on the same rows (decile income: none, so only country).
country_model <- fit_with(spec_terms("y"), dd, levels = c("clim_TM", "clim_RR", "y"))

# All deciles stacked: units are decile-countries; clustered by country.
stack <- do.call(rbind, lapply(DECILES, function(q) {
  x <- decile_data(q)
  x$decile <- q
  x$unit <- paste(x$GID_1, q)
  x
}))
pooled <- list(
  decile_fe = fit_with(SPECS$both, stack, fe = "decile^year + unit[t, t2]"),
  country_year_fe = fit_with(SPECS$both, stack, fe = "GID_0^year + unit[t, t2]")
)
pooled_tests <- list(decile_fe = income_tests(pooled$decile_fe),
                     country_year_fe = rbind(`T: decile income` = zero_test(pooled$country_year_fe, of("TM", "yq")),
                                             `P: decile income` = zero_test(pooled$country_year_fe, of("RR", "yq"))))

summary_table <- do.call(rbind, lapply(DECILES, function(q) {
  r <- per_decile[[q]]
  b <- coef(r$models$both)
  data.frame(decile = q, cor_income = r$cor_income,
             aic_country = r$aic[["country"]], aic_decile = r$aic[["decile"]], decile_vs_country = r$verdict,
             warm_yq = b[["pTM_X_yq"]], warm_y = b[["pTM_X_y"]], cold_yq = b[["nTM_X_yq"]], cold_y = b[["nTM_X_y"]],
             p_T_decile = r$tests["T: decile income", "p"], p_T_country = r$tests["T: country income", "p"],
             p_P_decile = r$tests["P: decile income", "p"], p_P_country = r$tests["P: country income", "p"],
             row.names = NULL)
}))

out_dir <- file.path("results", "decile_adaptation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(summary = summary_table,
             per_decile = lapply(per_decile, function(r) r[c("tests", "coefs", "verdict", "aic", "cor_income")]),
             country = coef_of(country_model), pooled = lapply(pooled, income_coefs), pooled_tests = pooled_tests,
             share_growth = sapply(DECILES, function(q) stats::quantile(dd[[paste0(q, "_dls")]], c(0.01, 0.5, 0.99))),
             share_constant = mean(dd$D1_dls == 0),
             sample = c(rows = nrow(dd), countries = length(unique(dd$GID_0)), first = min(dd$year), last = max(dd$year)),
             built = Sys.time()),
        file.path(out_dir, "estimates.rds"))
print(summary_table, digits = 3)
print(lapply(pooled, income_coefs)); print(pooled_tests)
cat("rows", nrow(dd), "countries", length(unique(dd$GID_0)), "years", range(dd$year),
    "; share of years with D1 share unchanged:", mean(dd$D1_dls == 0), "\n")
