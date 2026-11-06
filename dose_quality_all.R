# Quality of DOSE v2.14 regional GDP for EVERY country: the national total and
# growth implied by the sum of the regions against World Bank GDP
# (NY.GDP.PCAP.KD x SP.POP.TOTL, constant 2015 US$), and a pre-stated rule
# selecting a "robust" subnational dataset of whole countries whose regional
# statistics match the national ones in level and growth. Generalises
# dose_quality_check.R (CHN, IND, RUS, BRA) to all countries.
#
# DOSE GDP is grp_pc_lcu2015_usd x pop (constant 2015 local prices at DOSE's
# 2015 exchange rate): growth rates are comparable with WB, levels only up to a
# constant factor (the two 2015 exchange rates). Levels are therefore judged by
# the drift of the ratio DOSE / WB, not by its mean. Nothing is interpolated,
# extrapolated or converted.
#
# Two growth comparisons per country:
#   strict  : compare_national_subnational(), the same complete set of regions
#             in both years (complete_regional_coverage, stable signature);
#   matched : the sum of the regions with GDP and population in both
#             consecutive years, whatever the others report (needed when one
#             region is always missing, as in Brazil); its coverage is the
#             matched regions' population over WB population.
#
# Saved to results/dose_quality_all/ (.log, .csv, .png).
#
#   Rscript dose_quality_all.R

source("load_functions.R")
suppressPackageStartupMessages(library(ggplot2))
out_dir <- file.path("results", "dose_quality_all")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "dose_quality_all.log"), split = TRUE)
options(width = 250, dplyr.summarise.inform = FALSE)

# The rule, stated before looking at the results. A country passes when, over
# its compared years (strict comparison, whole period):
#   population coverage (regional sum / WB), MINIMUM over
#     compared years (the median lets partial-country years through,
#     e.g. West Germany before 1991)                     >= 0.90
#   correlation of total GDP growth with WB              >= 0.80
#   RMSE of the total GDP growth gap                     <= 2 pp
#   SD of log(GDP level ratio), complete-coverage years  <= 0.10
#   compared years                                       >= 15
# Near misses fail on the whole period but pass the same rule (a) on one
# contiguous sub-period (the longest, by compared years; at least 15 years,
# no single years dropped) or (b) with the matched-regions comparison, where
# the level drift is that of GDP per capita.
RULE <- list(min_cov = 0.90, min_corr = 0.80, max_rmse = 2, max_lsd = 0.10, min_years = 15)
BAD_PP <- 5   # a "bad year": |growth gap| > 5 pp
cat("Rule: population coverage >=", RULE$min_cov, "; growth corr >=", RULE$min_corr,
    "; growth RMSE <=", RULE$max_rmse, "pp; SD log level ratio <=", RULE$max_lsd,
    "; compared years >=", RULE$min_years, "\n")
cat("Bad year: |GDP growth gap| >", BAD_PP, "pp\n\n")

# Data -----------------------------------------------------------------------

sub <- read_gdp_levels("DOSE")
nat_all <- read_gdp_levels("WB")
nat <- nat_all %>% filter(GID_0 %in% unique(sub$GID_0))
comp <- compare_national_subnational(sub, nat)
gc <- comp$growth_comparison %>% arrange(iso3, year)
names_of <- sub %>% distinct(GID_0, country) %>% group_by(GID_0) %>% slice(1) %>% ungroup()
cname <- setNames(names_of$country, names_of$GID_0)

no_wb <- setdiff(unique(sub$GID_0), unique(gc$iso3))
cat("DOSE countries:", n_distinct(sub$GID_0), "; with no overlapping WB data or no valid DOSE GDP:",
    paste(no_wb, collapse = " "), "\n")

# DOSE break flags (StructChange; 0 none, nonzero codes 1-3 as in the file).
raw <- read_dose(data_file("DOSE_V2.14.csv"))
breaks <- raw %>%
  filter(!is.na(StructChange), StructChange != 0) %>%
  group_by(GID_0, year) %>%
  summarise(regions = n_distinct(GID_1), code = paste(sort(unique(StructChange)), collapse = "/")) %>%
  group_by(GID_0) %>%
  summarise(break_years = paste0(year, "(", code, ")", collapse = " "), n_break_years = n())

# Estimation panel rows (dose$clean) from the knitr cache.
cache_dir <- file.path("results", "knitr_cache", "weather_arms_robustness")
rdx <- list.files(cache_dir, pattern = "[.]rdx$", full.names = TRUE)
rdx <- rdx[order(file.info(rdx)$mtime)]
# Loaded into its own environment: the cache also holds objects named sub, nat, ...
cache_env <- new.env()
for (f in rdx) lazyLoad(sub("[.]rdx$", "", f), envir = cache_env)
panel <- cache_env$dose$clean %>% distinct(GID_0, GID_1, year)
panel_rows <- panel %>% count(GID_0, name = "panel_rows")
cat("dose$clean:", nrow(panel), "rows,", n_distinct(panel$GID_0), "countries,",
    paste(range(panel$year), collapse = "-"), "\n")

# Matched regions: regions with positive GDP and population in t and t-1.
valid <- sub %>% filter(is.finite(gdp), gdp > 0, is.finite(population), population > 0) %>%
  select(GID_0, GID_1, year, gdp, population)
matched <- valid %>%
  inner_join(valid %>% mutate(year = year + 1L), by = c("GID_0", "GID_1", "year"), suffix = c("", "_prev")) %>%
  group_by(GID_0, year) %>%
  summarise(m_regions = n(), m_pop = sum(population),
            m_g = 100 * (sum(gdp) / sum(gdp_prev) - 1),
            m_gpc = 100 * ((sum(gdp) / sum(population)) / (sum(gdp_prev) / sum(population_prev)) - 1)) %>%
  ungroup()

# Country-year table: every DOSE year overlapping WB.
cy <- gc %>%
  transmute(iso3, year, n_regions = n_regions_complete, complete = complete_regional_coverage,
            strict = eligible_for_comparison & is.finite(gdppc_growth_pct_sub) & is.finite(gdppc_growth_pct_nat),
            cov = population_sub / population_nat, ratio = gdp_sub / gdp_nat, ratio_pc = gdppc_sub / gdppc_nat,
            g_sub = gdp_growth_pct_sub, g_nat = gdp_growth_pct_nat,
            gpc_sub = gdppc_growth_pct_sub, gpc_nat = gdppc_growth_pct_nat,
            pop_nat = population_nat) %>%
  left_join(matched, by = c("iso3" = "GID_0", "year")) %>%
  mutate(m_cov = m_pop / pop_nat,
         matched_ok = is.finite(m_g) & is.finite(g_nat) & is.finite(gpc_nat),
         gap = if_else(strict, g_sub - g_nat, NA_real_),
         gap_pc = if_else(strict, gpc_sub - gpc_nat, NA_real_),
         m_gap = if_else(matched_ok, m_g - g_nat, NA_real_),
         m_gap_pc = if_else(matched_ok, m_gpc - gpc_nat, NA_real_))

# Statistics of one country (or one window of it) --------------------------
#
# view = "strict": growth on strict years, coverage = median population
#   coverage over them, level = total GDP ratio over complete-coverage years.
# view = "matched": growth on matched years, coverage = median matched
#   population coverage, level = GDP per capita ratio over all years.
gstats <- function(s, n, yr, prefix) {
  e <- s - n
  k <- length(e)
  w <- if (k) which.max(abs(e)) else integer(0)
  out <- c(n = k, mean = if (k) mean(e) else NA, mae = if (k) mean(abs(e)) else NA,
           rmse = if (k) sqrt(mean(e^2)) else NA,
           corr = if (k > 2 && sd(s) > 0 && sd(n) > 0) cor(s, n) else NA,
           worst_year = if (k) yr[w] else NA, worst_gap = if (k) e[w] else NA,
           bad = sum(abs(e) > BAD_PP))
  setNames(out, paste0(prefix, names(out)))
}
lstats <- function(r, yr) {
  ok <- is.finite(r) & r > 0
  r <- r[ok]; yr <- yr[ok]
  if (!length(r)) return(c(lev_years = 0, lev_median = NA, lev_lsd = NA, lev_min = NA, lev_max = NA, lev_step = NA, lev_step_year = NA))
  l <- log(r)
  d <- diff(l); dy <- diff(yr)
  st <- if (any(dy == 1)) which.max(abs(d) * (dy == 1)) else NA
  c(lev_years = length(r), lev_median = median(r), lev_lsd = if (length(r) > 1) sd(l) else NA,
    lev_min = min(r), lev_max = max(r),
    lev_step = if (is.na(st)) NA else d[st], lev_step_year = if (is.na(st)) NA else yr[st + 1])
}
country_stats <- function(d, view = c("strict", "matched")) {
  view <- match.arg(view)
  if (view == "strict") {
    x <- d[d$strict, ]
    cv <- if (nrow(x)) min(x$cov) else NA
    cvm <- if (nrow(x)) median(x$cov) else NA
    lv <- lstats(d$ratio[d$complete], d$year[d$complete])
    g <- gstats(x$g_sub, x$g_nat, x$year, "g_")
    gp <- gstats(x$gpc_sub, x$gpc_nat, x$year, "gpc_")
  } else {
    x <- d[d$matched_ok, ]
    cv <- if (nrow(x)) min(x$m_cov) else NA
    cvm <- if (nrow(x)) median(x$m_cov) else NA
    lv <- lstats(d$ratio_pc, d$year)
    g <- gstats(x$m_g, x$g_nat, x$year, "g_")
    gp <- gstats(x$m_gpc, x$gpc_nat, x$year, "gpc_")
  }
  c(cov = cv, cov_median = cvm, lv, g, gp)
}
passes <- function(s, rule = RULE, pc = FALSE) {
  ok <- is.finite(s["cov"]) && s["cov"] >= rule$min_cov &&
    is.finite(s["g_n"]) && s["g_n"] >= rule$min_years &&
    is.finite(s["g_corr"]) && s["g_corr"] >= rule$min_corr &&
    is.finite(s["g_rmse"]) && s["g_rmse"] <= rule$max_rmse &&
    is.finite(s["lev_lsd"]) && s["lev_lsd"] <= rule$max_lsd
  if (pc) ok <- ok && is.finite(s["gpc_corr"]) && s["gpc_corr"] >= rule$min_corr &&
    is.finite(s["gpc_rmse"]) && s["gpc_rmse"] <= rule$max_rmse
  unname(ok)
}
# Years a..b of one country: growth years t with a < t <= b (both levels in
# the window), level years in [a, b].
window_of <- function(d, a, b) {
  w <- d[d$year >= a & d$year <= b, ]
  w$strict <- w$strict & w$year > a
  w$matched_ok <- w$matched_ok & w$year > a
  w
}
# Longest contiguous window passing the rule (most compared years).
best_window <- function(d, view, rule = RULE, pc = FALSE) {
  yrs <- sort(unique(d$year))
  flag <- if (view == "strict") d$strict else d$matched_ok
  best <- NULL
  for (a in yrs) for (b in yrs[yrs >= a + rule$min_years]) {
    nyr <- sum(flag & d$year > a & d$year <= b)
    if (nyr < rule$min_years || (!is.null(best) && nyr <= best$n)) next
    w <- window_of(d, a, b)
    s <- country_stats(w, view)
    if (passes(s, rule, pc)) best <- list(a = a, b = b, n = nyr, s = s)
  }
  best
}

by_country <- split(cy, cy$iso3)

# Per-country table --------------------------------------------------------

basic <- cy %>% group_by(iso3) %>%
  summarise(first = min(year), last = max(year), dose_years = n(),
            regions_max = max(n_regions), regions_min = min(n_regions),
            regions_modal_share = mean(n_regions == as.integer(names(which.max(table(n_regions))))),
            complete_years = sum(complete), strict_years = sum(strict), matched_years = sum(matched_ok))
regions_ever <- sub %>% group_by(GID_0) %>% summarise(regions_listed = n_distinct(GID_1))

st <- do.call(rbind, lapply(by_country, country_stats, view = "strict"))
mt <- do.call(rbind, lapply(by_country, country_stats, view = "matched"))
colnames(mt) <- paste0("m_", colnames(mt))
tab <- basic %>%
  bind_cols(as_tibble(st), as_tibble(mt)) %>%
  left_join(regions_ever, by = c("iso3" = "GID_0")) %>%
  left_join(breaks, by = c("iso3" = "GID_0")) %>%
  left_join(panel_rows, by = c("iso3" = "GID_0")) %>%
  mutate(country = cname[iso3], panel_rows = coalesce(panel_rows, 0L),
         n_break_years = coalesce(n_break_years, 0L), .after = iso3)
tab$pass_strict <- sapply(by_country, function(d) passes(country_stats(d, "strict")))
tab$pass_matched <- sapply(by_country, function(d) passes(country_stats(d, "matched")))

# Clean sub-periods for those failing on the whole period.
win <- lapply(names(by_country), function(i) {
  d <- by_country[[i]]
  ws <- if (!tab$pass_strict[tab$iso3 == i]) best_window(d, "strict") else NULL
  wm <- if (!tab$pass_strict[tab$iso3 == i] && !tab$pass_matched[tab$iso3 == i]) best_window(d, "matched") else NULL
  tibble(iso3 = i,
         sub_strict = if (is.null(ws)) NA_character_ else paste0(ws$a, "-", ws$b),
         sub_strict_years = if (is.null(ws)) NA_real_ else ws$n,
         sub_strict_rmse = if (is.null(ws)) NA_real_ else ws$s["g_rmse"],
         sub_strict_corr = if (is.null(ws)) NA_real_ else ws$s["g_corr"],
         sub_strict_lsd = if (is.null(ws)) NA_real_ else ws$s["lev_lsd"],
         sub_strict_cov = if (is.null(ws)) NA_real_ else ws$s["cov"],
         sub_matched = if (is.null(wm)) NA_character_ else paste0(wm$a, "-", wm$b),
         sub_matched_years = if (is.null(wm)) NA_real_ else wm$n,
         sub_matched_rmse = if (is.null(wm)) NA_real_ else wm$s["g_rmse"],
         sub_matched_corr = if (is.null(wm)) NA_real_ else wm$s["g_corr"],
         sub_matched_cov = if (is.null(wm)) NA_real_ else wm$s["cov"])
}) %>% bind_rows()
tab <- tab %>% left_join(win, by = "iso3") %>%
  mutate(status = case_when(pass_strict ~ "pass",
                            pass_matched | !is.na(sub_strict) | !is.na(sub_matched) ~ "near miss",
                            TRUE ~ "fail"),
         # Regional sum equal to the national total (growth RMSE ~ 0): DOSE
         # has scaled the regions to the national accounts, so the check
         # cannot validate the regional figures themselves.
         scaled_to_national = is.finite(g_rmse) & g_rmse < 0.05)

# Countries failing everything, re-searched with population coverage relaxed
# to 0.75 (e.g. Brazil, four states always missing): reported, not counted.
RELAXED <- modifyList(RULE, list(min_cov = 0.75))
relaxed <- lapply(tab$iso3[tab$status == "fail"], function(i) {
  out <- lapply(c("strict", "matched"), function(v) {
    w <- best_window(by_country[[i]], v, RELAXED)
    if (is.null(w)) return(NULL)
    tibble(iso3 = i, view = v, span = paste0(w$a, "-", w$b), years = w$n, cov_min = w$s["cov"],
           cov_median = w$s["cov_median"], rmse = w$s["g_rmse"], corr = w$s["g_corr"], lev_lsd = w$s["lev_lsd"])
  })
  bind_rows(out)
}) %>% bind_rows()

# Failure reasons (strict, whole period).
tab <- tab %>% mutate(fails_on = purrr::pmap_chr(list(cov, g_n, g_corr, g_rmse, lev_lsd), function(cv, n, r, e, l) {
  x <- c(coverage = !(is.finite(cv) && cv >= RULE$min_cov),
         years = !(is.finite(n) && n >= RULE$min_years),
         corr = !(is.finite(r) && r >= RULE$min_corr),
         rmse = !(is.finite(e) && e <= RULE$max_rmse),
         level = !(is.finite(l) && l <= RULE$max_lsd))
  paste(names(x)[x], collapse = ",")
}))

# Panel rows inside the clean span (the whole DOSE span for passes; the
# sub-period for near misses) and in bad years.
span_of <- function(i) {
  r <- tab[tab$iso3 == i, ]
  if (r$status == "pass" || r$pass_matched) return(c(r$first, r$last))
  s <- if (!is.na(r$sub_strict)) r$sub_strict else r$sub_matched
  if (is.na(s)) return(c(NA, NA))
  as.integer(strsplit(s, "-")[[1]])
}
spans <- t(sapply(tab$iso3, span_of))
tab$span_first <- spans[, 1]; tab$span_last <- spans[, 2]
# Bad years: strict gap, else the matched gap when the matched regions cover
# the country (population coverage >= RULE$min_cov).
bad_cy <- cy %>% filter(abs(coalesce(gap, if_else(m_cov >= RULE$min_cov, m_gap, NA_real_))) > BAD_PP) %>% select(iso3, year)
tab <- tab %>%
  left_join(panel %>% inner_join(tab %>% select(iso3, span_first, span_last), by = c("GID_0" = "iso3")) %>%
              filter(year >= span_first, year <= span_last) %>% count(GID_0, name = "panel_rows_span"),
            by = c("iso3" = "GID_0")) %>%
  left_join(panel %>% semi_join(bad_cy, by = c("GID_0" = "iso3", "year")) %>% count(GID_0, name = "panel_rows_bad"),
            by = c("iso3" = "GID_0")) %>%
  mutate(panel_rows_span = coalesce(panel_rows_span, 0L), panel_rows_bad = coalesce(panel_rows_bad, 0L))
# Bad years the strict comparison cannot see (region set changes that year):
# matched gap > BAD_PP, matched coverage >= RULE$min_cov, inside the span.
hidden <- cy %>% inner_join(tab %>% select(iso3, span_first, span_last), by = "iso3") %>%
  filter(!strict, year > span_first, year <= span_last, m_cov >= RULE$min_cov, abs(m_gap) > BAD_PP) %>%
  group_by(iso3) %>% summarise(hidden_bad = paste0(year, "(", round(m_gap, 1), ")", collapse = " "))
tab <- tab %>% left_join(hidden, by = "iso3")

# Print ----------------------------------------------------------------------

show <- function(x) print(as.data.frame(x), row.names = FALSE)
r2 <- function(x) mutate(x, across(where(is.double) & !matches("year|first|last"), ~ round(.x, 2)))
cat("\n== Coverage and levels (strict view; level = sum of regions / WB, complete-coverage years) ==\n")
show(tab %>% arrange(status, iso3) %>%
       transmute(iso3, status, first, last, dose_years, regions = paste0(regions_min, "-", regions_max),
                 listed = regions_listed, modal_share = regions_modal_share, complete_years, strict_years,
                 matched_years, cov_min = cov, cov_median, m_cov_min = m_cov, lev_years, lev_median, lev_lsd, lev_min, lev_max,
                 lev_step, lev_step_year, pc_lev_lsd = m_lev_lsd, panel_rows) %>% r2())

cat("\n== Growth of total GDP and GDP per capita vs WB, strict (same complete regions both years), pp ==\n")
show(tab %>% arrange(status, iso3) %>%
       transmute(iso3, status, n = g_n, mean = g_mean, mae = g_mae, rmse = g_rmse, corr = g_corr,
                 worst = g_worst_year, worst_gap = g_worst_gap, bad = g_bad,
                 pc_rmse = gpc_rmse, pc_corr = gpc_corr, pc_bad = gpc_bad, fails_on) %>% r2())

cat("\n== Growth vs WB, matched regions (regions with GDP and population in both years), pp ==\n")
show(tab %>% arrange(status, iso3) %>%
       transmute(iso3, status, n = m_g_n, cov = m_cov, mean = m_g_mean, mae = m_g_mae, rmse = m_g_rmse,
                 corr = m_g_corr, worst = m_g_worst_year, worst_gap = m_g_worst_gap, bad = m_g_bad,
                 pc_rmse = m_gpc_rmse, pc_corr = m_gpc_corr, pass_matched) %>% r2())

cat("\n== DOSE break flags (StructChange != 0): year(code) ==\n")
show(tab %>% filter(n_break_years > 0) %>% transmute(iso3, status, n_break_years, break_years))

cat("\n== PASS (strict rule, whole DOSE period) ==\n")
show(tab %>% filter(status == "pass") %>%
       transmute(iso3, country, years = paste0(first, "-", last), regions = paste0(regions_min, "-", regions_max),
                 compared = g_n, cov_min = cov, cov_median, lev_median, lev_lsd, lev_range = paste0(round(lev_min, 2), "-", round(lev_max, 2)),
                 rmse = g_rmse, corr = g_corr, pc_rmse = gpc_rmse, pc_corr = gpc_corr, bad = g_bad,
                 breaks = n_break_years, scaled = scaled_to_national, panel_rows, panel_rows_bad, hidden_bad) %>% r2())
cat("scaled = DOSE regional sum equals the national total (growth RMSE < 0.05 pp): agreement by construction\n")
cat("hidden_bad = years the strict comparison skips (region set changes) with a matched-regions gap > ", BAD_PP, " pp (gap)
", sep = "")
cat("panel_rows_bad = dose$clean rows in country-years with |gap| >", BAD_PP, "pp (strict, else matched); those years:\n")
show(tab %>% filter(status == "pass", g_bad + m_g_bad > 0) %>% select(iso3) %>%
       semi_join(x = cy, y = ., by = "iso3") %>% semi_join(bad_cy, by = c("iso3", "year")) %>%
       transmute(iso3, year, strict_gap = gap, matched_gap = m_gap, m_cov) %>% r2())

cat("\n== NEAR MISSES (fail on the whole period; pass on a contiguous sub-period or with matched regions) ==\n")
show(tab %>% filter(status == "near miss") %>%
       transmute(iso3, country, years = paste0(first, "-", last), fails_on, pass_matched,
                 sub_strict, sub_strict_years, sub_strict_rmse, sub_strict_corr, sub_strict_lsd, sub_strict_cov,
                 sub_matched, sub_matched_years, sub_matched_rmse, sub_matched_corr, sub_matched_cov,
                 panel_rows, panel_rows_span, break_years, hidden_bad) %>% r2())

cat("\n== Failing countries with a clean sub-period if population coverage is relaxed to 0.75 (not counted) ==\n")
show(relaxed %>% r2())

cat("\n== FAIL ==\n")
show(tab %>% filter(status == "fail") %>%
       transmute(iso3, country, years = paste0(first, "-", last), fails_on, compared = g_n, cov_min = cov,
                 rmse = g_rmse, corr = g_corr, lev_lsd, m_n = m_g_n, m_cov, m_rmse = m_g_rmse,
                 m_corr = m_g_corr, panel_rows) %>% r2() %>% arrange(desc(panel_rows)))

# Shares of the panel and of world GDP ---------------------------------------

wld <- nat_all %>% filter(GID_0 == "WLD", year == 2015) %>% pull(gdp)
gdp15 <- nat %>% filter(year == 2015) %>% select(iso3 = GID_0, gdp15 = gdp)
tab <- tab %>% left_join(gdp15, by = "iso3") %>% mutate(world_share_2015 = gdp15 / wld)
share <- function(set, rows) c(countries = length(set),
                               panel_rows = sum(rows[tab$iso3 %in% set]),
                               panel_share = sum(rows[tab$iso3 %in% set]) / nrow(panel),
                               world_gdp_share_2015 = sum(tab$world_share_2015[tab$iso3 %in% set], na.rm = TRUE))
pass_set <- tab$iso3[tab$status == "pass"]
near_set <- tab$iso3[tab$status == "near miss"]
cat("\n== Shares: panel rows (dose$clean) and WB world GDP in 2015 ==\n")
print(round(rbind(
  `robust (pass)` = share(pass_set, tab$panel_rows),
  `pass + near misses, their clean spans` = share(c(pass_set, near_set), tab$panel_rows_span),
  `all DOSE countries compared with WB` = share(tab$iso3, tab$panel_rows)), 3))

# Sensitivity: one threshold at a time, loosened and tightened ---------------

variants <- list(
  baseline = RULE,
  `coverage >= 0.80` = modifyList(RULE, list(min_cov = 0.80)),
  `coverage >= 0.95` = modifyList(RULE, list(min_cov = 0.95)),
  `corr >= 0.70` = modifyList(RULE, list(min_corr = 0.70)),
  `corr >= 0.90` = modifyList(RULE, list(min_corr = 0.90)),
  `RMSE <= 3` = modifyList(RULE, list(max_rmse = 3)),
  `RMSE <= 1.5` = modifyList(RULE, list(max_rmse = 1.5)),
  `level SD <= 0.15` = modifyList(RULE, list(max_lsd = 0.15)),
  `level SD <= 0.05` = modifyList(RULE, list(max_lsd = 0.05)),
  `years >= 10` = modifyList(RULE, list(min_years = 10)),
  `years >= 20` = modifyList(RULE, list(min_years = 20)))
strict_stats <- lapply(by_country, country_stats, view = "strict")
panel_stats <- lapply(by_country, function(d) country_stats(window_of(d, 1970L, 9999L), "strict"))
matched_stats <- lapply(by_country, country_stats, view = "matched")
sens_row <- function(v, stats, rule = RULE, pc = FALSE) {
  p <- names(stats)[sapply(stats, passes, rule = rule, pc = pc)]
  tibble(variant = v, n_pass = length(p), added = paste(setdiff(p, pass_set), collapse = " "),
         dropped = paste(setdiff(pass_set, p), collapse = " "),
         panel_share = round(sum(tab$panel_rows[tab$iso3 %in% p]) / nrow(panel), 3))
}
sens <- bind_rows(
  lapply(names(variants), function(v) sens_row(v, strict_stats, variants[[v]])),
  sens_row("also GDP per capita growth", strict_stats, pc = TRUE),
  sens_row("panel period only (1970-)", panel_stats),
  sens_row("matched regions instead of strict", matched_stats))
cat("\n== Sensitivity of the passing list (strict, whole period), one threshold at a time ==\n")
show(sens)

# Save ---------------------------------------------------------------------

write.csv(tab, file.path(out_dir, "country_quality.csv"), row.names = FALSE)
write.csv(cy, file.path(out_dir, "country_year_comparison.csv"), row.names = FALSE)
write.csv(sens, file.path(out_dir, "rule_sensitivity.csv"), row.names = FALSE)
write.csv(comp$country_year_subregion_driver, file.path(out_dir, "country_year_region_driver.csv"), row.names = FALSE)

# Figures (Okabe-Ito) --------------------------------------------------------

status_cols <- c(pass = "#009E73", `near miss` = "#E69F00", fail = "#D55E00")
status_shapes <- c(pass = 16, `near miss` = 17, fail = 4)
fd <- tab %>% filter(is.finite(g_rmse), is.finite(g_corr)) %>%
  mutate(status = factor(status, names(status_cols)))
p1 <- ggplot(fd, aes(g_corr, g_rmse, colour = status, shape = status)) +
  annotate("rect", xmin = RULE$min_corr, xmax = Inf, ymin = -Inf, ymax = RULE$max_rmse,
           fill = "#009E73", alpha = 0.06) +
  geom_vline(xintercept = RULE$min_corr, colour = "grey60", linewidth = 0.3, linetype = 2) +
  geom_hline(yintercept = RULE$max_rmse, colour = "grey60", linewidth = 0.3, linetype = 2) +
  geom_point(aes(size = pmax(panel_rows, 1))) +
  ggrepel::geom_text_repel(aes(label = iso3), size = 2.6, colour = "grey25", show.legend = FALSE,
                           max.overlaps = 30, segment.colour = "grey70", segment.size = 0.2) +
  scale_colour_manual(values = status_cols, name = NULL) +
  scale_shape_manual(values = status_shapes, name = NULL) +
  scale_size_area(max_size = 6, name = "Panel rows", breaks = c(100, 500, 1500)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(-0.4, 1)) +
  labs(x = "Correlation of total GDP growth, DOSE regions summed vs World Bank",
       y = "RMSE of the growth gap (pp)",
       title = "DOSE v2.14: does the sum of the regions track national GDP growth?",
       subtitle = paste0("Strict comparison (same complete set of regions in both years), whole DOSE period.
",
                         "Status also uses population coverage, level drift and years compared; off-scale (corr < -0.4): ",
                         paste(fd$iso3[fd$g_corr < -0.4], collapse = " "))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "growth_fit_by_country.png"), p1, width = 10, height = 7, dpi = 150, bg = "white")

gap_cols <- c(`< 2 pp` = "#009E73", `2-5 pp` = "#E69F00", `> 5 pp` = "#D55E00",
              `matched only, < 2 pp` = "#56B4E9", `matched only, 2-5 pp` = "#F0E442",
              `matched only, > 5 pp` = "#CC79A7", `not compared` = "grey88")
ord <- tab %>% mutate(status = factor(status, names(status_cols))) %>% arrange(desc(status), desc(iso3))
hd <- cy %>% mutate(
  m_gap = if_else(m_cov >= RULE$min_cov, m_gap, NA_real_),
  a = coalesce(abs(gap), abs(m_gap)),
  cls = case_when(is.finite(gap) & a < 2 ~ "< 2 pp", is.finite(gap) & a <= BAD_PP ~ "2-5 pp",
                  is.finite(gap) ~ "> 5 pp",
                  is.finite(m_gap) & a < 2 ~ "matched only, < 2 pp",
                  is.finite(m_gap) & a <= BAD_PP ~ "matched only, 2-5 pp",
                  is.finite(m_gap) ~ "matched only, > 5 pp", TRUE ~ "not compared"),
  cls = factor(cls, names(gap_cols)),
  lab = paste0(iso3, " (", substr(as.character(ord$status[match(iso3, ord$iso3)]), 1, 4), ")"),
  lab = factor(lab, paste0(ord$iso3, " (", substr(as.character(ord$status), 1, 4), ")")))
p2 <- ggplot(hd, aes(year, lab, fill = cls)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = gap_cols, name = "|GDP growth gap|", drop = FALSE) +
  scale_x_continuous(breaks = seq(1960, 2020, 10), expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "DOSE regions summed vs World Bank: GDP growth gap by country and year",
       subtitle = "Strict comparison where the same complete set of regions reports in both years; otherwise matched regions covering >= 90% of WB population; countries ordered by status") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "right", axis.text.y = element_text(size = 6.5))
ggsave(file.path(out_dir, "growth_gap_heatmap.png"), p2, width = 11, height = 12, dpi = 150, bg = "white")
sink()
