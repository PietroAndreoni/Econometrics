# PWT against World Bank GDP, for the handout (handout_preferred_specification.Rmd):
# the growth outliers dropped, and how the two sources' growth and income
# differ. Saved to results/handout_preferred_specification/sources.rds.
#
# Both growth rates are real growth at constant national prices (PWT rgdpna /
# pop, WB NY.GDP.PCAP.KD), so their gap comes from the national-accounts
# sources, vintages, splices and population series; PPP against market
# exchange rates only enters the income levels (PWT 2021 PPP US$, WB 2015 US$).
# Rows: the response explorer's ERA5 panels, population weights.
#
#   Rscript handout_sources.R

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
})
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
pwt <- P$panels[["PWT | ERA5 | population"]]$rows
wb <- P$panels[["WB | ERA5 | population"]]$rows

# Outliers: the 5% tails of PWT x ERA5 growth from the start year, their
# annotated causes, the non-weather ones dropped (weather_arms_robustness.Rmd;
# read from its knitr cache, as apps/response_explorer/build_models.R does).
cache <- new.env()
rdx <- list.files(file.path("results", "knitr_cache", "weather_arms_robustness"), pattern = "[.]rdx$",
                  full.names = TRUE)
for (f in rdx[order(file.info(rdx)$mtime)]) lazyLoad(sub("[.]rdx$", "", f), envir = cache)
tails <- cache$OUTLIERS
key <- function(x) paste(x$GID_0, x$year)
outliers <- list(
  tails = nrow(tails), weather_kept = sum(tails$weather != "no"), dropped = sum(tails$weather == "no"),
  share = cache$OUTLIER_SHARE, start_year = cache$START_YEAR,
  # Dropped in the countries and years of the estimation rows (trailing income).
  dropped_sample = sum(pwt$outlier_5 & is.finite(pwt$inc_trailing)),
  categories = sort(table(tails$category[tails$weather == "no" & key(tails) %in% key(pwt[is.finite(pwt$inc_trailing), ])]),
                    decreasing = TRUE))

# Growth: the gap PWT - WB (pp) on the rows both have.
keep <- c("GID_0", "GID_1", "year", "t", "t2", P$outcome, "inc_level", "inc_trailing", "conflict_low", "conflict_high",
          "outlier_5", "agr_level")
m <- inner_join(pwt[, keep], wb[, c("GID_0", "year", P$outcome, "inc_level", "inc_trailing")], by = c("GID_0", "year"),
                suffix = c("", "_wb")) %>%
  mutate(g_pwt = .data[[P$outcome]], g_wb = .data[[paste0(P$outcome, "_wb")]], gap = 100 * (g_pwt - g_wb),
         abs_gap = abs(gap), region = countrycode::countrycode(GID_0, "iso3c", "region"),
         decade = paste0(10 * (year %/% 10), "s"), conflict = conflict_low > 0 | conflict_high > 0) %>%
  arrange(GID_0, year) %>%
  group_by(GID_0) %>%
  mutate(gap_lag = ifelse(year - lag(year) == 1, lag(gap), NA)) %>%
  ungroup()
by_group <- function(x, ...) {
  x %>% group_by(...) %>%
    summarise(rows = n(), sd = sd(gap), mean_abs = mean(abs_gap), over_2pp = mean(abs_gap > 2), .groups = "drop")
}
r2_of <- function(f) unname(fixest::r2(feols(f, m, notes = FALSE), "r2"))
growth <- list(
  rows = nrow(m), countries = n_distinct(m$GID_0),
  cor = cor(m$g_pwt, m$g_wb), sd_pwt = 100 * sd(m$g_pwt), sd_gap = sd(m$gap), median_abs = median(m$abs_gap),
  over = c(`1` = mean(m$abs_gap > 1), `2` = mean(m$abs_gap > 2), `5` = mean(m$abs_gap > 5)),
  r2 = c(country = r2_of(gap ~ 1 | GID_0), year = r2_of(gap ~ 1 | year), both = r2_of(gap ~ 1 | GID_0 + year)),
  ar1 = cor(m$gap, m$gap_lag, use = "complete.obs"),
  decade = by_group(m, decade), region = by_group(m, region) %>% arrange(desc(sd)),
  quintile = by_group(m %>% mutate(quintile = ntile(inc_level_wb, 5)), quintile),
  conflict_outlier = by_group(m, conflict, outlier_5),
  abs_gap_model = coeftable(feols(abs_gap ~ inc_level_wb + conflict + outlier_5 + I(year >= 2000) + agr_level, m,
                                  cluster = ~GID_0)),
  drift = m %>% group_by(GID_0) %>%
    summarise(rows = n(), mean = mean(gap), t = mean(gap) / (sd(gap) / sqrt(n())), .groups = "drop") %>%
    filter(rows >= 20) %>% arrange(desc(abs(t))) %>% head(8),
  largest = m %>% arrange(desc(abs_gap)) %>% select(GID_0, year, g_pwt, g_wb, gap, outlier_5) %>% head(10))

# Does the gap follow the weather? The gap on the handout's weather arms (fixed
# SD deviations from the 20-year trailing mean), year effects and country
# quadratic trends, on the estimation rows (PWT outliers dropped, trailing income).
w <- pwt[!pwt$outlier_5 & is.finite(pwt$inc_trailing), ] %>%
  inner_join(wb[, c("GID_0", "year", P$outcome)], by = c("GID_0", "year"), suffix = c("", "_wb"))
for (v in c("TM", "RR")) {
  x <- w[[paste0("dm_", v, "_lag_20_unstd_quadratic")]] / w[[paste0("sd_", v, "_all")]]
  w[[paste0("p", v)]] <- pmax(x, 0)
  w[[paste0("n", v)]] <- pmax(-x, 0)
}
w$gap <- w[[P$outcome]] - w[[paste0(P$outcome, "_wb")]]
on_weather <- function(y) {
  feols(as.formula(paste(y, "~ pTM + nTM + pRR + nRR | year + GID_1[t, t2]")), w, cluster = ~GID_0, notes = FALSE)
}
fits <- list(PWT = on_weather(P$outcome), WB = on_weather(paste0(P$outcome, "_wb")), gap = on_weather("gap"))
weather <- list(rows = nrow(w), coef = 100 * sapply(fits, coef), se = 100 * sapply(fits, se),
                p_gap = fixest::wald(fits$gap, print = FALSE)$p)

# Income levels: log(PWT PPP / WB MER) of the long-run and trailing income.
l <- m %>% filter(is.finite(inc_trailing), is.finite(inc_trailing_wb)) %>% mutate(ratio = inc_trailing - inc_trailing_wb)
levels <- list(
  mean_ratio = exp(mean(l$ratio)),
  slope = unname(coef(lm(ratio ~ inc_trailing_wb, l))[2]), r2 = summary(lm(ratio ~ inc_trailing_wb, l))$r.squared,
  rank_cor = cor(l$inc_trailing, l$inc_trailing_wb, method = "spearman"),
  sd_within = sd(l$ratio - ave(l$ratio, l$GID_0)), sd_between = sd(tapply(l$ratio, l$GID_0, mean)),
  quintile = l %>% mutate(quintile = ntile(inc_trailing_wb, 5)) %>% group_by(quintile) %>%
    summarise(mer = exp(mean(inc_trailing_wb)), ppp = exp(mean(inc_trailing)), ratio = exp(mean(ratio)), .groups = "drop"))

out <- list(outliers = outliers, growth = growth, weather = weather, levels = levels, built = Sys.time())
saveRDS(out, file.path("results", "handout_preferred_specification", "sources.rds"))
str(out[c("outliers", "weather")], max.level = 1)
print(growth[c("cor", "sd_gap", "over", "r2", "ar1")]); print(weather$p_gap); print(levels[c("slope", "rank_cor")])
