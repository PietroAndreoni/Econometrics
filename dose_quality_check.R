# Quality of DOSE v2.14 regional GDP for China, India, Russia and Brazil: the
# national growth implied by the sum of the regions against World Bank GDP
# (NY.GDP.MKTP via NY.GDP.PCAP.KD x population: constant 2015 US$). DOSE is
# grp_pc_lcu2015_usd x pop (constant 2015 local prices at DOSE's 2015
# exchange rate), so growth rates are comparable; levels only up to the two
# 2015 exchange rates. A growth rate of the sum is compared only when the same
# complete set of regions reports in both years (compare_national_subnational).
# Saved to results/dose_quality_check/ (.log, .csv, .png).
#
#   Rscript dose_quality_check.R

source("load_functions.R")
suppressPackageStartupMessages(library(ggplot2))
COUNTRIES <- c(CHN = "China", IND = "India", RUS = "Russia", BRA = "Brazil")
out_dir <- file.path("results", "dose_quality_check")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "dose_quality_check.log"), split = TRUE)

sub <- read_gdp_levels("DOSE") %>% filter(GID_0 %in% names(COUNTRIES))
nat <- read_gdp_levels("WB") %>% filter(GID_0 %in% names(COUNTRIES))
comp <- compare_national_subnational(sub, nat)
gc <- comp$growth_comparison %>% mutate(country = COUNTRIES[iso3])

# Coverage: years with DOSE data, years whose growth can be compared.
coverage <- gc %>%
  group_by(country) %>%
  summarise(first = min(year), last = max(year), years = n(),
            regions = paste(range(n_regions_complete), collapse = "-"),
            compared = sum(eligible_for_comparison),
            not_compared = paste(year[!eligible_for_comparison & year > min(year)], collapse = " "),
            .groups = "drop")
options(width = 220)
cat("Coverage (DOSE years overlapping WB; growth compared where the same complete set of regions reports in both years)\n")
print(as.data.frame(coverage))

# Growth of total GDP and of GDP per capita: sum of the regions less WB (pp).
stats_of <- function(x, s, n) {
  e <- x[[s]] - x[[n]]
  c(mean = mean(e), mae = mean(abs(e)), rmse = sqrt(mean(e^2)), corr = cor(x[[s]], x[[n]]),
    sd_dose = sd(x[[s]]), sd_wb = sd(x[[n]]), max_abs = max(abs(e)), max_year = x$year[which.max(abs(e))])
}
ok <- gc %>% filter(eligible_for_comparison, is.finite(gdppc_growth_pct_sub), is.finite(gdppc_growth_pct_nat))
for (k in c(gdp = "GDP", gdppc = "GDP per capita", population = "population")) {
  s <- paste0(names(which(c(gdp = "GDP", gdppc = "GDP per capita", population = "population") == k)), "_growth_pct_sub")
  tab <- t(sapply(split(ok, ok$country), stats_of, s = s, n = sub("_sub$", "_nat", s)))
  cat(sprintf("\nGrowth of %s, %% a year: DOSE sum of regions less WB (mean, mean absolute, RMSE, max), correlation, SDs\n", k))
  print(round(tab, 2))
}

# Levels: the sum of the regions over WB (GDP in 2015 US$, population), by
# decade, to show drift between the two.
lev <- gc %>%
  mutate(gdp_ratio = gdp_sub / gdp_nat, pop_ratio = population_sub / population_nat,
         decade = paste0(10 * (year %/% 10), "s")) %>%
  filter(complete_regional_coverage)
cat("\nLevels, DOSE sum of regions / WB (complete-coverage years), mean by decade\n")
print(as.data.frame(lev %>% group_by(country, decade) %>%
                      summarise(gdp = round(mean(gdp_ratio), 3), population = round(mean(pop_ratio), 3), years = n(),
                                .groups = "drop")))

# The largest discrepancies in GDP growth per country, and the region that
# contributes most to each (exact shift-share attribution).
drivers <- comp$country_year_subregion_driver %>%
  mutate(country = COUNTRIES[iso3]) %>%
  select(country, year, primary_driver_region_name, primary_driver_contribution_pp,
         primary_driver_regional_growth_pct)
top <- gc %>% filter(eligible_for_comparison) %>%
  group_by(country) %>% slice_max(abs_gdp_growth_discrepancy_pp, n = 6) %>% ungroup() %>%
  select(country, year, dose = gdp_growth_pct_sub, wb = gdp_growth_pct_nat, gap = gdp_growth_discrepancy_pp) %>%
  left_join(drivers, by = c("country", "year")) %>%
  mutate(across(where(is.numeric) & !year, ~ round(.x, 2)))
cat("\nLargest gaps in GDP growth (pp), with the region contributing most to each\n")
print(as.data.frame(top))

# Matched regions: the growth of the sum of the regions with GDP and
# population in both years (consecutive), whatever the others report, against
# WB; with their share of the regions' GDP in the 2010 cross-section (the year
# Brazil has every state). Brazil has no year with every state but 2010: four
# states (Rio de Janeiro, Rio Grande do Sul, Rio Grande do Norte, Mato Grosso
# do Sul) have population only in 2010, and before 1985 it is quinquennial.
valid <- sub %>% filter(is.finite(gdp), gdp > 0, is.finite(population), population > 0) %>%
  select(GID_0, GID_1, year, gdp, population)
pairs <- valid %>% inner_join(valid %>% mutate(year = year + 1L), by = c("GID_0", "GID_1", "year"),
                              suffix = c("", "_prev"))
share_2010 <- sub %>% filter(year == 2010) %>% group_by(GID_0) %>%
  mutate(share = gdp_pc * population / sum(gdp_pc * population, na.rm = TRUE)) %>% ungroup() %>%
  select(GID_1, share)
matched <- pairs %>% left_join(share_2010, by = "GID_1") %>%
  group_by(GID_0, year) %>%
  summarise(regions = n(), share_2010 = sum(share, na.rm = TRUE),
            dose = 100 * (sum(gdp) / sum(gdp_prev) - 1),
            dose_pc = 100 * ((sum(gdp) / sum(population)) / (sum(gdp_prev) / sum(population_prev)) - 1),
            .groups = "drop") %>%
  inner_join(gc %>% select(GID_0 = iso3, year, wb = gdp_growth_pct_nat, wb_pc = gdppc_growth_pct_nat), by = c("GID_0", "year")) %>%
  filter(is.finite(wb)) %>%
  mutate(country = COUNTRIES[GID_0], gap = dose - wb)
cat("\nMatched regions (sum of the regions reporting in both years): growth of GDP, DOSE less WB (pp)\n")
print(as.data.frame(matched %>% group_by(country) %>%
  summarise(years = n(), first = min(year), last = max(year), regions = paste(range(regions), collapse = "-"),
            share_2010 = round(median(share_2010), 2), mean = round(mean(gap), 2), mae = round(mean(abs(gap)), 2),
            rmse = round(sqrt(mean(gap^2)), 2), corr = round(cor(dose, wb), 2),
            corr_pc = round(cor(dose_pc, wb_pc), 2), .groups = "drop")))
cat("\nBrazil, matched regions, every year (growth %)\n")
print(as.data.frame(matched %>% filter(GID_0 == "BRA") %>%
  transmute(year, regions, share_2010 = round(share_2010, 2), dose = round(dose, 2), wb = round(wb, 2), gap = round(gap, 2))))

write.csv(gc, file.path(out_dir, "country_year_growth_comparison.csv"), row.names = FALSE)
write.csv(matched, file.path(out_dir, "matched_regions_growth_comparison.csv"), row.names = FALSE)
write.csv(comp$subregion_attribution %>% filter(iso3 %in% names(COUNTRIES)),
          file.path(out_dir, "region_attribution.csv"), row.names = FALSE)

# Figures: GDP per capita growth (sum of the regions and WB), and the level
# ratios. Okabe-Ito: blue WB, vermillion DOSE.
cols <- c(`World Bank` = "#0072B2", `DOSE, sum of regions` = "#D55E00")
# The DOSE series is that of the matched regions (with every region, the
# complete sum), on a full grid of years so that gaps show as breaks.
g <- gc %>% transmute(country, year, `World Bank` = gdppc_growth_pct_nat) %>%
  left_join(matched %>% transmute(country, year, `DOSE, sum of regions` = dose_pc), by = c("country", "year")) %>%
  tidyr::pivot_longer(c(`World Bank`, `DOSE, sum of regions`), names_to = "series", values_to = "growth") %>%
  tidyr::complete(country, series, year = seq(min(year), max(year)))
p1 <- ggplot(g, aes(year, growth, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 0.9, na.rm = TRUE) +
  facet_wrap(~country, scales = "free", ncol = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = NULL, y = "Growth of GDP per capita (%)",
       title = "National growth of GDP per capita: DOSE regions aggregated vs World Bank",
       subtitle = "DOSE: the regions with GDP and population in both years (Brazil: 23 of 27 states; quinquennial before 1985)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "growth_dose_vs_wb.png"), p1, width = 10, height = 7, dpi = 150, bg = "white")

r <- lev %>% filter(country != "Brazil") %>% select(country, year, GDP = gdp_ratio, Population = pop_ratio) %>%
  tidyr::pivot_longer(c(GDP, Population), names_to = "series", values_to = "ratio")
p2 <- ggplot(r, aes(year, ratio, colour = series)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~country, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c(GDP = "#D55E00", Population = "#009E73"), name = NULL) +
  labs(x = NULL, y = "DOSE sum of regions / World Bank",
       title = "Levels: DOSE regions summed over World Bank national totals",
       subtitle = "GDP in constant 2015 US$ (each at its own 2015 exchange rate); years in which every region listed that year reports (India: 24-32 states); Brazil has one such year (2010)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "levels_dose_over_wb.png"), p2, width = 10, height = 7, dpi = 150, bg = "white")
sink()
