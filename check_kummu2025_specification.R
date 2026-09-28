# Reproduce the preferred specification with geographically validated KUMMU data.
# Prerequisite: Rscript audit_kummu2025_mapping.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(fixest); library(ggplot2)})
source("rmd_chunks.R")
out <- "results/kummu2025_robustness"
write_out <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
meta <- jsonlite::fromJSON(file.path(out, "sources", "zenodo.json"))
polygon_meta <- meta$files[meta$files$key == "polyg_adm1_gdp_perCapita_1990_2024.gpkg", ]
stopifnot(paste0("md5:", unname(tools::md5sum(file.path(out, "sources", "kummu_adm1.gpkg")))) == polygon_meta$checksum)
cw <- read.csv(file.path(out, "regional_crosswalk.csv"))
raw <- read.csv("data/tabulated_adm1_gdp_perCapita.csv", check.names = FALSE)
stopifnot(!anyDuplicated(raw$GID_nmbr), !anyDuplicated(cw$GID_nmbr),
          !anyDuplicated(cw$GID_1[cw$accepted]))
econ_path <- "data/econ_processed_DOSE_V2_11-KUMMU2018-KUMMU2025-WDI-WB.parquet"
stored <- arrow::read_parquet(econ_path) %>% filter(econ_source == "KUMMU2025", gadm_level == "gadm1")
raw_long <- raw %>% select(GID_nmbr, iso3, matches("^[0-9]{4}$")) %>%
  pivot_longer(matches("^[0-9]{4}$"), names_to = "year", values_to = "grp_pc_usd") %>%
  mutate(year = as.integer(year)) %>% group_by(GID_nmbr) %>% arrange(year, .by_group = TRUE) %>%
  mutate(lgrp_pc_usd = log(grp_pc_usd),
    dlgrp_pc_usd = if_else(year == lag(year) + 1L, lgrp_pc_usd - lag(lgrp_pc_usd), NA_real_)) %>% ungroup()
# Check 3a: independently reproduce the existing cached economic series, proving
# that its problem is the geographic identifier rather than the growth formula.
old <- raw_long %>% filter(GID_nmbr >= 1000000) %>%
  mutate(GID_1 = paste0(iso3, ".", GID_nmbr %% 1000, "_1"))
cache_check <- inner_join(old, stored, by = c("GID_1", "year"), suffix = c("_raw", "_stored"))
stopifnot(nrow(cache_check) == nrow(old), nrow(cache_check) == nrow(stored),
  max(abs(cache_check$grp_pc_usd_raw - cache_check$grp_pc_usd_stored), na.rm = TRUE) < 1e-10,
  max(abs(cache_check$dlgrp_pc_usd_raw - cache_check$dlgrp_pc_usd_stored), na.rm = TRUE) < 1e-12)
corrected <- raw_long %>% inner_join(cw %>% filter(accepted) %>% select(GID_nmbr, GID_1), by = "GID_nmbr") %>%
  transmute(year, GID_0 = iso3, GID_1, grp_pc_usd, lgrp_pc_usd, dlgrp_pc_usd,
            econ_source = "KUMMU2025", gadm_level = "gadm1")
stopifnot(!anyDuplicated(corrected[c("GID_1", "year")]),
          all(substr(corrected$GID_1, 1, 3) == corrected$GID_0))
verified_path <- file.path(out, "econ_kummu2025_verified.parquet")
arrow::write_parquet(corrected, verified_path)
nb <- load_notebook_env("test_functions.Rmd")
cache <- new.env(parent = emptyenv())
# Push down the notebook's unchanged climate filters before materializing its
# large multi-weight parquet files. Return only its requested columns.
nb$read_parquet_files <- function(files, columns) {
  stopifnot(length(files) == 1)
  path <- files[[1]]
  if (grepl("data_(rr|tm)_", basename(path))) {
    key <- basename(path)
    if (!exists(key, cache, inherits = FALSE)) {
      x <- arrow::open_dataset(path) %>%
        filter(gadm_level == "gadm1", climate_source == "era5", weight == "pop", weight_year == "2015") %>%
        select(all_of(columns)) %>% collect()
      assign(key, x, cache)
    }
    return(get(key, cache))
  }
  arrow::read_parquet(path, col_select = all_of(columns))
}
build <- function(source, path) nb$build_dat(econ_data = source,
  econ_files_by_level = c(gadm1 = path),
  climate_rr_files_by_level = c(gadm1 = "data/data_rr_gadm1_era5_pop-area_2000-2015.parquet"),
  climate_tm_files_by_level = c(gadm1 = "data/data_tm_gadm1_era5_pop-area_2000-2015.parquet")) %>%
  mutate(zzTM = pmax(1.5, abs((TM - mean_TM_lag_30) / sd_TM_lag_30)) - 1.5)
cat("Building DOSE panel...\n")
dose <- build("DOSE_V2_11", econ_path)
cat("Building verified KUMMU panel...\n")
kummu <- build("KUMMU2025", verified_path)
needed <- c("dlgrp_pc_usd", "TM", "mean_TM_all", "RR", "mean_RR_all", "zzTM", "abs_zRRp")
eligible <- function(x) x %>% filter(if_all(all_of(needed), is.finite))
d <- eligible(dose); k <- eligible(kummu)
stopifnot(!anyDuplicated(d[c("GID_1", "year")]), !anyDuplicated(k[c("GID_1", "year")]))
climate_ids <- unique(kummu$GID_1)
cw$climate_joined <- cw$accepted & cw$GID_1 %in% climate_ids
write_out(cw, "regional_crosswalk_with_climate")
unmatched <- corrected %>% filter(!is.na(dlgrp_pc_usd)) %>%
  anti_join(k %>% select(GID_1, year), by = c("GID_1", "year"))
write_out(unmatched %>% count(GID_0, GID_1, name = "unusable_region_years"), "climate_sample_exclusions")
common <- inner_join(d %>% select(GID_1, year), k %>% select(GID_1, year), by = c("GID_1", "year"))
dc <- semi_join(d, common, by = c("GID_1", "year")) %>% arrange(GID_1, year)
kc <- semi_join(k, common, by = c("GID_1", "year")) %>% arrange(GID_1, year)
climate_vars <- setdiff(needed, "dlgrp_pc_usd")
stopifnot(nrow(dc) == nrow(kc), all(dc$GID_1 == kc$GID_1), all(dc$year == kc$year),
          isTRUE(all.equal(dc[climate_vars], kc[climate_vars], tolerance = 1e-12)))
cat("Check 3 passed: raw/cache growth round-trip, unique climate joins, identical common-sample regressors.\n")
samples <- list(DOSE_full = d, KUMMU_verified = k, DOSE_common = dc, KUMMU_common = kc,
  KUMMU_strict995 = k %>% filter(GID_1 %in% cw$GID_1[cw$accepted & cw$kummu_overlap >= .995 & cw$gadm_overlap >= .995]),
  DOSE_common_years = d %>% filter(year >= min(k$year), year <= max(k$year)))
# Import the preferred formula directly from its notebook chunk, suppressing
# only the display table. This guards against testing a stale copied formula.
nb$data <- dose
nb$fit_climate_model <- function(terms, model_data = nb$data) nb$make_formula(terms)
for (expr in parse(text = extract_rmd_chunks("test_functions.Rmd", "test-double deviation")[[1]])) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-"))) eval(expr, nb)
}
preferred_formula <- nb$m_pref
models <- list(); coefs <- list(); summaries <- list(); tests <- list()
for (label in names(samples)) {
  x <- samples[[label]]
  for (spec in c("preferred", "quadratic")) {
    f <- if (spec == "preferred") preferred_formula else nb$make_formula(paste(nb$base_climate, "+ zTM^2 + zRR^2"))
    fit <- feols(f, data = x, panel.id = nb$pan_id, cluster = ~GID_1)
    id <- paste(label, spec, sep = "__"); models[[id]] <- fit
    used <- x[fixest::obs(fit), ]
    ct <- as.data.frame(coeftable(fit)); ci <- confint(fit)
    coefs[[id]] <- data.frame(sample = label, specification = spec, term = rownames(ct),
      estimate = ct[,1], se = ct[,2], p = ct[,4], conf_low = ci[,1], conf_high = ci[,2])
    summaries[[id]] <- data.frame(sample = label, specification = spec, n = nobs(fit),
      regions = n_distinct(used$GID_1), countries = n_distinct(used$GID_0),
      first_year = min(used$year), last_year = max(used$year), wr2 = as.numeric(fitstat(fit, "wr2")[[1]]))
    if (spec == "preferred") {
      w <- wald(fit, keep = "zzTM|abs_zRRp", print = FALSE)
      tests[[label]] <- data.frame(sample = label, joint_F = w$stat, joint_p = w$p)
      # Country clustering checks inference when KUMMU shares national movements.
      cc <- as.data.frame(coeftable(fit, vcov = ~GID_0))
      ci2 <- confint(fit, vcov = ~GID_0)
      coefs[[paste0(id, "_country")]] <- data.frame(sample = label, specification = "preferred_country_cluster",
        term = rownames(cc), estimate = cc[,1], se = cc[,2], p = cc[,4], conf_low = ci2[,1], conf_high = ci2[,2])
    }
  }
}
coefficients <- bind_rows(coefs)
# On identical regressors, regressing the outcome difference estimates the
# coefficient difference and its covariance without treating the datasets as
# independent samples.
delta <- dc
delta$dlgrp_pc_usd <- kc$dlgrp_pc_usd - dc$dlgrp_pc_usd
diff_fit <- feols(preferred_formula, data = delta, panel.id = nb$pan_id, cluster = ~GID_1)
stopifnot(max(abs(coef(diff_fit) -
  (coef(models$KUMMU_common__preferred) - coef(models$DOSE_common__preferred)))) < 1e-8)
diff_ct <- as.data.frame(coeftable(diff_fit)); diff_ci <- confint(diff_fit)
write_out(data.frame(term = rownames(diff_ct), estimate = diff_ct[,1], se = diff_ct[,2],
  p = diff_ct[,4], conf_low = diff_ci[,1], conf_high = diff_ci[,2]), "common_coefficient_differences")
models$KUMMU_minus_DOSE_common <- diff_fit
write_out(coefficients, "coefficients"); write_out(bind_rows(summaries), "sample_sizes")
write_out(bind_rows(tests), "joint_tests")
saveRDS(models, file.path(out, "models.rds"))
writeLines(deparse(preferred_formula), file.path(out, "preferred_formula.txt"))
write_out(data.frame(n = nrow(dc), growth_correlation = cor(dc$dlgrp_pc_usd, kc$dlgrp_pc_usd),
  growth_rmse = sqrt(mean((dc$dlgrp_pc_usd - kc$dlgrp_pc_usd)^2))), "common_growth_comparison")
plot_data <- coefficients %>% filter(specification == "preferred")
p <- ggplot(plot_data, aes(x = estimate * 100, y = sample)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_errorbar(aes(xmin = conf_low * 100, xmax = conf_high * 100), orientation = "y", width = .2) +
  geom_point() + facet_wrap(~term, scales = "free_x") + theme_minimal() +
  labs(x = "Coefficient in log-growth percentage points (95% CI)", y = NULL,
       title = "Preferred specification: DOSE and verified KUMMU2025", subtitle = "Standard errors clustered by region")
ggsave(file.path(out, "preferred_coefficients.png"), p, width = 10, height = 5, dpi = 170)
writeLines(capture.output(sessionInfo()), file.path(out, "model_session_info.txt"))
print(plot_data %>% select(sample, term, estimate, se, p))
print(bind_rows(summaries)); print(bind_rows(tests))
