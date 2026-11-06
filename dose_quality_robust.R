# Robust version of the DOSE v2.14 country quality screen (dose_quality_all.R):
# the regional sum of DOSE GDP against World Bank national GDP (both constant
# 2015 US$), restricted to the estimation period (growth years 1970 onward)
# and tolerant of a few isolated bad years.
#
# The question: do countries fail the screen because of one or two bad years
# (a splice, a one-off data error) that could simply be dropped from the
# regression panel, or because the regional series is genuinely unreliable?
#
# Rule (pre-stated; same thresholds as dose_quality_all.R where they carry
# over). Over a country's compared growth years t >= 1970:
#   bad year  : |growth gap| > BAD_PP (5 pp), or population coverage < 0.90
#               in a compared year, or (strict view) a year the strict
#               comparison skips because the region set changes but whose
#               matched-regions gap (coverage >= 0.90) exceeds BAD_PP;
#   allowance : at most K_BAD (2) bad years, listed explicitly;
#   on the remaining (clean) years:
#     compared years                                   >= 15
#     population coverage, minimum                     >= 0.90
#     correlation of total GDP growth with WB          >= 0.80
#     RMSE of the growth gap                           <= 2 pp
#     level drift net of steps: SD of the step-adjusted
#       log level ratio (cumulated log growth gaps of
#       the clean years; bad years contribute 0)       <= 0.10
# The old SD of log(DOSE / WB level) is reported but not used: a one-off step
# moves it for every later year although it only contaminates one growth year.
# Robust statistics over all compared years (median |gap|, share of years
# with |gap| <= 2 pp, Spearman correlation) are reported alongside.
#
# Views: "strict" (primary) = compare_national_subnational(), the same
# complete set of regions in both years; "matched" = regions with GDP and
# population in both consecutive years, coverage = their population / WB.
#
# Status: pass (strict, whole 1970+ period); pass (matched) = fails strict
# but passes matched; sub-period = passes on the longest contiguous window of
# growth years (>= 15 clean years, same allowance); fail.
#
# Bad years are listed with the primary driver region from the attribution in
# compare_national_subnational() (strict years) or the same shift-share on
# the matched regions (years only the matched comparison sees). They are the
# country-years to drop or flag in the regression panel; this script changes
# neither the panel nor the outlier rules. Nothing is interpolated,
# extrapolated or converted.
#
# Saved to results/dose_quality_robust/ (.log, .csv).
#
#   Rscript dose_quality_robust.R

source("load_functions.R")
out_dir <- file.path("results", "dose_quality_robust")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "dose_quality_robust.log"), split = TRUE)
options(width = 250, dplyr.summarise.inform = FALSE)

YEAR0 <- 1970L
PAR <- list(min_cov = 0.90, min_corr = 0.80, max_rmse = 2, max_drift = 0.10, max_lsd = 0.10,
            min_years = 15, bad_pp = 5, k = 2, k_share = NULL, bad_rule = "abs", z_cut = 5, z_floor = 2)
cat("Growth years >=", YEAR0, "\n")
cat("Bad year: |gap| >", PAR$bad_pp, "pp, or coverage <", PAR$min_cov,
    "in a compared year, or (strict) a skipped year with matched gap >", PAR$bad_pp, "pp\n")
cat("Allowance: <=", PAR$k, "bad years. Clean years: >=", PAR$min_years, "; coverage >=", PAR$min_cov,
    "; corr >=", PAR$min_corr, "; RMSE <=", PAR$max_rmse, "pp; SD step-adjusted log level ratio <=",
    PAR$max_drift, "\n\n")

# Data -----------------------------------------------------------------------

dose_sub <- read_gdp_levels("DOSE")
nat_all <- read_gdp_levels("WB")
nat <- nat_all %>% filter(GID_0 %in% unique(dose_sub$GID_0))
comp <- compare_national_subnational(dose_sub, nat)
gc <- comp$growth_comparison %>% arrange(iso3, year)
names_of <- dose_sub %>% distinct(GID_0, country) %>% group_by(GID_0) %>% slice(1) %>% ungroup()
cname <- setNames(names_of$country, names_of$GID_0)

# DOSE metadata by country-year: break flags (StructChange != 0), source
# vintage (`version`, modal over regions) and the DOSE exchange rate (LCU per
# US$, national), whose log change diagnoses currency-unit errors.
raw <- read_dose(data_file("DOSE_V2.14.csv"))
meta <- raw %>%
  group_by(GID_0, year) %>%
  summarise(sc_code = paste(sort(unique(StructChange[!is.na(StructChange) & StructChange != 0])), collapse = "/"),
            sc_regions = sum(!is.na(StructChange) & StructChange != 0),
            version = names(which.max(table(version))),
            fx = median(fx, na.rm = TRUE)) %>%
  group_by(GID_0) %>% arrange(year, .by_group = TRUE) %>%
  mutate(consec = year == lag(year) + 1L,
         version_change = consec & version != lag(version),
         fx_dlog = if_else(consec & fx > 0 & lag(fx) > 0, 100 * (log(fx) - log(lag(fx))), NA_real_)) %>%
  ungroup() %>%
  mutate(sc_code = na_if(sc_code, ""), year = as.integer(year))

# Estimation panel rows (dose$clean) from the knitr cache, in its own
# environment: the cache also holds objects named sub, nat, ...
cache_dir <- file.path("results", "knitr_cache", "weather_arms_robustness")
rdx <- list.files(cache_dir, pattern = "[.]rdx$", full.names = TRUE)
rdx <- rdx[order(file.info(rdx)$mtime)]
cache_env <- new.env()
for (f in rdx) lazyLoad(sub("[.]rdx$", "", f), envir = cache_env)
panel <- cache_env$dose$clean %>% distinct(GID_0, GID_1, year) %>% mutate(year = as.integer(year))
rm(cache_env)
panel_cy <- panel %>% count(GID_0, year, name = "rows")
cat("dose$clean:", nrow(panel), "rows,", n_distinct(panel$GID_0), "countries,",
    paste(range(panel$year), collapse = "-"), "\n")

# Matched regions: positive GDP and population in t and t-1.
valid <- dose_sub %>% filter(is.finite(gdp), gdp > 0, is.finite(population), population > 0) %>%
  select(GID_0, GID_1, region_name, year, gdp, population)
reg_growth <- valid %>%
  inner_join(valid %>% mutate(year = year + 1L) %>% select(GID_0, GID_1, year, gdp_prev = gdp, population_prev = population),
             by = c("GID_0", "GID_1", "year"))
matched <- reg_growth %>%
  group_by(GID_0, year) %>%
  summarise(m_regions = n(), m_pop = sum(population),
            m_g = 100 * (sum(gdp) / sum(gdp_prev) - 1),
            m_gpc = 100 * ((sum(gdp) / sum(population)) / (sum(gdp_prev) / sum(population_prev)) - 1)) %>%
  ungroup()

# Country-year table, growth years >= YEAR0.
cy <- gc %>%
  transmute(iso3, year, n_regions = n_regions_complete, complete = complete_regional_coverage,
            strict = eligible_for_comparison & is.finite(gdppc_growth_pct_sub) & is.finite(gdppc_growth_pct_nat),
            cov = population_sub / population_nat, ratio = gdp_sub / gdp_nat, ratio_pc = gdppc_sub / gdppc_nat,
            g_sub = gdp_growth_pct_sub, g_nat = gdp_growth_pct_nat, pop_nat = population_nat, gdp_nat) %>%
  left_join(matched, by = c("iso3" = "GID_0", "year")) %>%
  left_join(meta %>% select(GID_0, year, sc_code, sc_regions, version, version_change, fx_dlog),
            by = c("iso3" = "GID_0", "year")) %>%
  mutate(m_cov = m_pop / pop_nat,
         matched_ok = is.finite(m_g) & is.finite(g_nat),
         gap = if_else(strict, g_sub - g_nat, NA_real_),
         m_gap = if_else(matched_ok, m_g - g_nat, NA_real_)) %>%
  filter(year >= YEAR0)
by_country <- split(cy, cy$iso3)

# Statistics of one country (or one window of it) ---------------------------

is_bad <- function(e, P) {
  if (P$bad_rule == "abs") return(abs(e) > P$bad_pp)
  s <- 1.4826 * mad(e)
  z <- if (is.finite(s) && s > 0) abs(e - median(e)) / s else rep(0, length(e))
  z > P$z_cut & abs(e) > P$z_floor
}
k_allowed <- function(m, P) if (is.null(P$k_share)) P$k else floor(P$k_share * m)
fmt_gap <- function(y, e) if (length(y)) paste0(y, "(", sprintf("%+.1f", e), ")", collapse = " ") else NA_character_
safe_cor <- function(a, b, method = "pearson")
  if (length(a) > 2 && isTRUE(sd(a) > 0) && isTRUE(sd(b) > 0)) cor(a, b, method = method) else NA_real_

# Returns a list of statistics; `full = FALSE` skips descriptive extras (window search).
rstats <- function(d, view = "strict", P = PAR, full = TRUE) {
  if (view == "strict") {
    x <- d[d$strict, ]; s <- x$g_sub; cv <- x$cov
    h <- d[!d$strict & d$matched_ok & is.finite(d$m_cov) & d$m_cov >= P$min_cov, ]
    hbad <- abs(h$m_gap) > P$bad_pp
    lev_r <- d$ratio[d$complete]
  } else {
    x <- d[d$matched_ok, ]; s <- x$m_g; cv <- x$m_cov
    h <- x[0, ]; hbad <- logical(0)
    lev_r <- d$ratio_pc
  }
  n <- x$g_nat; e <- s - n; yr <- x$year; m <- length(e)
  gapbad <- if (m) is_bad(e, P) else logical(0)
  lowcov <- !(is.finite(cv) & cv >= P$min_cov)
  bad <- gapbad | lowcov
  good <- !bad
  kmax <- k_allowed(m, P)
  n_bad <- sum(bad) + sum(hbad)
  eg <- e[good]
  dl <- log1p(s[good] / 100) - log1p(n[good] / 100)
  L <- c(0, cumsum(dl))
  r <- list(
    n = m, n_clean = sum(good), n_gapbad = sum(gapbad & !lowcov), n_lowcov = sum(lowcov),
    n_hidden = sum(hbad), n_bad = n_bad, k_allowed = kmax,
    cov_min = if (any(good)) min(cv[good]) else NA_real_,
    corr = safe_cor(s[good], n[good]),
    rmse = if (any(good)) sqrt(mean(eg^2)) else NA_real_,
    drift_sd = if (length(L) > 2) sd(L) else NA_real_)
  crit <- c(years = is.finite(r$n_clean) && r$n_clean >= P$min_years,
            coverage = is.finite(r$cov_min) && r$cov_min >= P$min_cov,
            corr = is.finite(r$corr) && r$corr >= P$min_corr,
            rmse = is.finite(r$rmse) && r$rmse <= P$max_rmse,
            drift = is.finite(r$drift_sd) && r$drift_sd <= P$max_drift)
  r$other_ok <- all(crit)
  r$pass <- r$other_ok && n_bad <= kmax
  if (!full) return(r)
  lr <- lev_r[is.finite(lev_r) & lev_r > 0]
  old_lsd <- if (length(lr) > 1) sd(log(lr)) else NA_real_
  w <- if (m) which.max(abs(e)) else integer(0)
  c(r, list(
    fails_on = paste(c(if (n_bad > kmax) paste0("bad_years(", n_bad, ">", kmax, ")"),
                       names(crit)[!crit]), collapse = ","),
    k_needed = if (r$other_ok) n_bad else NA_integer_,
    mean_gap = if (any(good)) mean(eg) else NA_real_,
    mae = if (any(good)) mean(abs(eg)) else NA_real_,
    sd_dlog_clean = if (sum(good) > 1) 100 * sd(dl) else NA_real_,
    drift_total = sum(dl), drift_range = diff(range(L)),
    old_lsd = old_lsd,
    cov_min_all = if (m) min(cv) else NA_real_,
    med_abs_gap = if (m) median(abs(e)) else NA_real_,
    share_le2 = if (m) mean(abs(e) <= 2) else NA_real_,
    spearman = safe_cor(s, n, "spearman"),
    corr_all = safe_cor(s, n), rmse_all = if (m) sqrt(mean(e^2)) else NA_real_,
    worst_year = if (m) yr[w] else NA_real_, worst_gap = if (m) e[w] else NA_real_,
    max_z = if (m > 2 && mad(e) > 0) max(abs(e - median(e)) / (1.4826 * mad(e))) else NA_real_,
    old_pass = m >= P$min_years && is.finite(min(cv)) && min(cv) >= P$min_cov &&
      isTRUE(safe_cor(s, n) >= P$min_corr) && sqrt(mean(e^2)) <= P$max_rmse &&
      is.finite(old_lsd) && old_lsd <= P$max_lsd,
    bad_gap_years = fmt_gap(yr[gapbad & !lowcov], e[gapbad & !lowcov]),
    lowcov_years = if (any(lowcov)) paste0(yr[lowcov], "(", round(cv[lowcov], 2), ")", collapse = " ") else NA_character_,
    hidden_bad_years = fmt_gap(h$year[hbad], h$m_gap[hbad]),
    bad_year_list = paste(sort(c(yr[bad], h$year[hbad])), collapse = " ")))
}
# Correlation and slope of the growth gap on the DOSE exchange-rate change.
fx_diag <- function(d, view = "strict") {
  x <- if (view == "strict") d[d$strict, ] else d[d$matched_ok, ]
  e <- if (view == "strict") x$gap else x$m_gap
  ok <- is.finite(e) & is.finite(x$fx_dlog)
  if (sum(ok) < 5 || sd(x$fx_dlog[ok]) == 0) return(c(fx_corr = NA, fx_slope = NA, fx_corr_v1 = NA))
  v1 <- ok & x$version == "1"
  c(fx_corr = cor(e[ok], x$fx_dlog[ok]), fx_slope = unname(coef(lm(e[ok] ~ x$fx_dlog[ok]))[2]),
    fx_corr_v1 = if (sum(v1) >= 5 && sd(x$fx_dlog[v1]) > 0) cor(e[v1], x$fx_dlog[v1]) else NA)
}
to_row <- function(r) as_tibble(r[!vapply(r, is.null, logical(1))])

# Longest contiguous window of growth years [a, b] passing the robust rule.
best_window <- function(d, view, P = PAR) {
  flag <- if (view == "strict") d$strict else d$matched_ok
  yrs <- sort(unique(d$year[flag]))
  best <- NULL
  for (a in yrs) for (b in yrs[yrs >= a]) {
    nyr <- sum(flag & d$year >= a & d$year <= b)
    if (nyr < P$min_years || (!is.null(best) && nyr < best$n_clean)) next
    w <- d[d$year >= a & d$year <= b, ]
    r <- rstats(w, view, P, full = FALSE)
    # Most clean years; ties go to the window with fewer bad years.
    better <- is.null(best) || r$n_clean > best$n_clean || (r$n_clean == best$n_clean && r$n_bad < best$n_bad)
    if (r$pass && better) best <- list(a = a, b = b, n = nyr, n_clean = r$n_clean, n_bad = r$n_bad)
  }
  best
}

# Per-country table ------------------------------------------------------------

st <- lapply(by_country, rstats, view = "strict")
mt <- lapply(by_country, rstats, view = "matched")
tab <- bind_rows(lapply(st, to_row), .id = "iso3")
tab_m <- bind_rows(lapply(mt, to_row), .id = "iso3") %>% rename_with(~ paste0("m_", .x), -iso3)
fxd <- bind_rows(lapply(by_country, function(d) as_tibble(as.list(fx_diag(d)))), .id = "iso3")
basic <- cy %>% group_by(iso3) %>%
  summarise(first = min(year), last = max(year), dose_years = n(),
            regions = paste0(min(n_regions), "-", max(n_regions)), strict_years = sum(strict),
            matched_years = sum(matched_ok),
            break_years = {b <- !is.na(sc_code); if (any(b)) paste0(year[b], "(", sc_code[b], ")", collapse = " ") else NA_character_},
            version_change_years = if (any(version_change, na.rm = TRUE)) paste(year[which(version_change)], collapse = " ") else NA_character_)
tab <- basic %>% left_join(tab, by = "iso3") %>% left_join(tab_m, by = "iso3") %>%
  left_join(fxd, by = "iso3") %>%
  mutate(country = cname[iso3], .after = iso3)

win <- lapply(tab$iso3, function(i) {
  if (tab$pass[tab$iso3 == i] || tab$m_pass[tab$iso3 == i]) return(NULL)
  d <- by_country[[i]]
  ws <- best_window(d, "strict"); wm <- best_window(d, "matched")
  use <- if (!is.null(ws) && (is.null(wm) || ws$n_clean >= wm$n_clean)) list(w = ws, v = "strict") else
    if (!is.null(wm)) list(w = wm, v = "matched") else NULL
  if (is.null(use)) return(NULL)
  r <- rstats(d[d$year >= use$w$a & d$year <= use$w$b, ], use$v)
  tibble(iso3 = i, span_view = use$v, span = paste0(use$w$a, "-", use$w$b),
         span_first = use$w$a, span_last = use$w$b, span_clean = r$n_clean, span_bad = r$bad_year_list,
         span_rmse = r$rmse, span_corr = r$corr, span_drift = r$drift_sd, span_cov = r$cov_min)
}) %>% bind_rows()
tab <- tab %>% left_join(win, by = "iso3") %>%
  mutate(status = case_when(pass ~ "pass", m_pass ~ "pass (matched)", !is.na(span) ~ "sub-period", TRUE ~ "fail"),
         scaled_to_national = is.finite(rmse_all) & rmse_all < 0.05)

# Previous screen (dose_quality_all.R, whole DOSE period).
prev <- read.csv(file.path("results", "dose_quality_all", "country_quality.csv"), stringsAsFactors = FALSE) %>%
  select(iso3, old_status = status, old_fails_on = fails_on)
tab <- tab %>% left_join(prev, by = "iso3") %>%
  mutate(verdict = case_when(
    pass & old_status == "pass" ~ "pass (unchanged)",
    pass & n_bad > 0 ~ paste0("now passes: ", n_bad, " isolated bad year(s) set aside"),
    pass ~ "now passes: 1970+ window / step-robust drift",
    m_pass ~ paste0("passes with matched regions only (strict: ", fails_on, ")"),
    n < PAR$min_years & m_n < PAR$min_years ~ "too few years",
    TRUE ~ paste0(
      if_else(is.finite(k_needed), paste0("a few too many bad years (passes at k=", k_needed, ")"),
              fails_on),
      if_else(!is.na(span), paste0("; clean sub-period ", span, " (", span_view, ")"), ""))))

# Per-country failure classes (for the explanation of every non-pass).
classify <- function(r) {
  if (r$pass) return("pass")
  if (r$m_pass) return("matched-only pass")
  if (r$n < PAR$min_years && r$m_n < PAR$min_years) return("too few years")
  # Too few strict years (region set changes, e.g. CAN, BRA): judge the matched view.
  if (r$n < PAR$min_years) {
    m <- r[, grepl("^m_", names(r))]
    names(m) <- sub("^m_", "", names(m))
    r <- m
  }
  f <- strsplit(r$fails_on, ",")[[1]]
  out <- c()
  if (any(grepl("coverage", f)) || r$n_lowcov > r$k_allowed) out <- c(out, "low coverage")
  if (r$n_gapbad + r$n_hidden > r$k_allowed + 2 || "rmse" %in% f) out <- c(out, "systematic noise")
  else if (r$n_gapbad + r$n_hidden > r$k_allowed) out <- c(out, "a few too many bad years")
  if ("corr" %in% f) out <- c(out, "low correlation")
  if ("drift" %in% f) out <- c(out, "level drift")
  if ("years" %in% f) out <- c(out, "too few clean years")
  paste(unique(out), collapse = " + ")
}
tab$problem <- vapply(seq_len(nrow(tab)), function(j) classify(tab[j, ]), character(1))

# Bad years: table with drivers ------------------------------------------------

drv <- comp$country_year_subregion_driver %>%
  transmute(iso3, year = as.integer(year), primary_driver_region_id, primary_driver_region_name,
            primary_driver_contribution_pp, primary_driver_share_of_discrepancy, primary_driver_regional_growth_pct)
# Same shift-share on matched regions, for years only the matched view sees,
# and the GDP-weighted share of regions moving with the gap by at least half
# its size (near 1: the whole country jumps, i.e. a national-level splice,
# deflator or currency problem; small: a regional data error).
reg_contrib <- reg_growth %>%
  inner_join(cy %>% select(iso3, year, g_nat, m_gap, gap), by = c("GID_0" = "iso3", "year")) %>%
  group_by(GID_0, year) %>%
  mutate(share_prev = gdp_prev / sum(gdp_prev), g_i = 100 * (gdp / gdp_prev - 1),
         contrib = share_prev * (g_i - g_nat), e = coalesce(gap, m_gap)) %>%
  summarise(m_driver_id = GID_1[which.max(sign(first(e)) * contrib)],
            m_driver_name = region_name[which.max(sign(first(e)) * contrib)],
            m_driver_contribution_pp = max(sign(first(e)) * contrib) * sign(first(e)),
            m_driver_growth_pct = g_i[which.max(sign(first(e)) * contrib)],
            uniform_share = sum(share_prev[sign(first(e)) * (g_i - g_nat) > 0.5 * abs(first(e))]),
            regional_growth_median = median(g_i), regional_growth_min = min(g_i), regional_growth_max = max(g_i)) %>%
  ungroup()

bad_years_of <- function(i, view = "strict", P = PAR) {
  d <- by_country[[i]]
  if (view == "strict") {
    x <- d %>% filter(strict) %>% mutate(comparison = "strict", e = gap, s = g_sub, cv = cov)
    h <- d %>% filter(!strict, matched_ok, m_cov >= P$min_cov, abs(m_gap) > P$bad_pp) %>%
      mutate(comparison = "matched (region set changes)", e = m_gap, s = m_g, cv = m_cov, reason = "gap")
  } else {
    x <- d %>% filter(matched_ok) %>% mutate(comparison = "matched", e = m_gap, s = m_g, cv = m_cov)
    h <- x[0, ] %>% mutate(reason = character(0))
  }
  if (!nrow(x)) return(NULL)
  x <- x %>% mutate(gb = is_bad(e, P), lc = !(is.finite(cv) & cv >= P$min_cov)) %>%
    filter(gb | lc) %>% mutate(reason = if_else(lc, if_else(gb, "low coverage + gap", "low coverage"), "gap"))
  bind_rows(x, h) %>% select(iso3, year, comparison, reason, dose_growth = s, wb_growth = g_nat, gap = e,
                             coverage = cv, sc_code, version, version_change, fx_dlog)
}
bad_all <- bind_rows(lapply(tab$iso3, bad_years_of)) %>%
  left_join(drv, by = c("iso3", "year")) %>%
  left_join(reg_contrib, by = c("iso3" = "GID_0", "year")) %>%
  mutate(driver_source = if_else(is.na(primary_driver_region_id), "matched shift-share", "compare_national_subnational"),
         primary_driver_region_id = coalesce(primary_driver_region_id, m_driver_id),
         primary_driver_region_name = coalesce(primary_driver_region_name, m_driver_name),
         primary_driver_contribution_pp = coalesce(primary_driver_contribution_pp, m_driver_contribution_pp),
         primary_driver_regional_growth_pct = coalesce(primary_driver_regional_growth_pct, m_driver_growth_pct)) %>%
  select(-starts_with("m_driver")) %>%
  left_join(panel_cy, by = c("iso3" = "GID_0", "year")) %>%
  left_join(tab %>% select(iso3, status, span_first, span_last), by = "iso3") %>%
  mutate(panel_rows = coalesce(rows, 0L), .keep = "unused") %>%
  mutate(in_robust_set = status == "pass" |
           (status == "sub-period" & year >= span_first & year <= span_last)) %>%
  rename(GID_0 = iso3) %>% arrange(GID_0, year)

# Print ------------------------------------------------------------------------

show <- function(x) print(as.data.frame(x), row.names = FALSE)
r2 <- function(x) mutate(x, across(where(is.double) & !matches("year|first|last"), ~ round(.x, 2)))
cat("\n== All countries, strict view, growth years >=", YEAR0, "==\n")
show(tab %>% arrange(factor(status, c("pass", "pass (matched)", "sub-period", "fail")), iso3) %>%
       transmute(iso3, status, old = old_status, yrs = paste0(first, "-", last), n, clean = n_clean,
                 bad = n_bad, k = k_allowed, cov_min, corr, rmse, mean_gap, drift_sd, drift_total, old_lsd,
                 med_abs = med_abs_gap, le2 = share_le2, spear = spearman, corr_all, rmse_all, max_z,
                 m_n, m_bad = m_n_bad, m_corr, m_rmse, m_drift = m_drift_sd, fx_corr, fails_on) %>% r2())

cat("\n== PASS: robust set (strict, whole 1970+ period, <=", PAR$k, "bad years) ==\n")
show(tab %>% filter(status == "pass") %>%
       transmute(iso3, country, yrs = paste0(first, "-", last), regions, n, clean = n_clean, cov_min, corr, rmse,
                 mean_gap, drift_sd, old_lsd, med_abs = med_abs_gap, le2 = share_le2, spear = spearman,
                 scaled = scaled_to_national, old_status, bad_gap_years, hidden_bad_years, lowcov_years) %>% r2())
cat("scaled = regional sum equals the national total (RMSE < 0.05 pp over all years): agreement by construction\n")

cat("\n== PASS with matched regions only ==\n")
show(tab %>% filter(status == "pass (matched)") %>%
       transmute(iso3, country, strict_fails = fails_on, m_n, m_clean = m_n_clean, m_cov = m_cov_min, m_corr, m_rmse,
                 m_drift = m_drift_sd, m_bad_gap_years, old_status) %>% r2())

cat("\n== SUB-PERIOD: longest contiguous window passing the robust rule ==\n")
show(tab %>% filter(status == "sub-period") %>%
       transmute(iso3, country, whole = fails_on, span, span_view, span_clean, span_bad, span_cov, span_corr,
                 span_rmse, span_drift, old_status) %>% r2())

cat("\n== Every previously failing or near-miss country under the new screen ==\n")
show(tab %>% filter(old_status != "pass" | is.na(old_status)) %>%
       arrange(factor(status, c("pass", "pass (matched)", "sub-period", "fail")), iso3) %>%
       transmute(iso3, old_status, old_fails_on, status, problem, n, bad = n_bad, k_needed, clean = n_clean,
                 cov_min, corr, rmse, drift_sd, m_corr, m_rmse, span, bad_years = bad_year_list) %>% r2())

cat("\n== Previously passing countries under the new screen ==\n")
show(tab %>% filter(old_status == "pass") %>% transmute(iso3, status, fails_on, n_bad, bad_year_list))

# Large economies: bad years and drivers -----------------------------------------

big <- c("ITA", "GBR", "CAN", "DEU", "ESP", "AUS", "MEX", "IND", "IDN", "TUR", "RUS", "BRA", "CHN", "ZAF",
         "THA", "NOR", "AUT", "NZL")
cat("\n== Large economies: summary ==\n")
show(tab %>% filter(iso3 %in% big) %>% arrange(match(iso3, big)) %>%
       transmute(iso3, status, problem, n, bad = n_bad, clean = n_clean, cov_min_all, cov_min, corr, rmse,
                 drift_sd, old_lsd, med_abs = med_abs_gap, le2 = share_le2, spear = spearman, m_n, m_cov_min_all,
                 m_bad = m_n_bad, m_corr, m_rmse, fx_corr, fx_slope, fx_corr_v1, span, break_years, version_change_years) %>% r2())
cat("\n== Large economies: every bad year (strict view; matched for CAN and BRA, which have no strict years) ==\n")
big_bad <- bind_rows(lapply(big, function(i) {
  v <- if (sum(by_country[[i]]$strict) == 0) "matched" else "strict"
  bad_years_of(i, v)
})) %>%
  left_join(drv, by = c("iso3", "year")) %>% left_join(reg_contrib, by = c("iso3" = "GID_0", "year")) %>%
  mutate(driver = coalesce(primary_driver_region_name, m_driver_name),
         drv_pp = coalesce(primary_driver_contribution_pp, m_driver_contribution_pp),
         drv_g = coalesce(primary_driver_regional_growth_pct, m_driver_growth_pct))
show(big_bad %>% transmute(iso3, year, comparison = substr(comparison, 1, 7), reason, dose = dose_growth, wb = wb_growth,
                           gap, cov = coverage, driver, drv_pp, drv_g, uniform = uniform_share,
                           reg_med = regional_growth_median, sc = sc_code, vchg = version_change, fx = fx_dlog) %>% r2())
cat("uniform = lagged-GDP share of regions whose growth deviates from WB growth in the gap's direction by > half the gap\n")
cat("sc = DOSE StructChange codes that year; vchg = DOSE source version changes that year; fx = 100 x dlog DOSE LCU/US$\n")

# Italy in detail ----------------------------------------------------------------

cat("\n== ITALY: every year ==\n")
ita <- by_country$ITA
show(ita %>% transmute(year, strict, cov, ratio, dose = g_sub, wb = g_nat, gap, m_gap, fx = fx_dlog,
                       version, vchg = version_change, sc = sc_code, sc_regions) %>% r2())
it_nom <- raw %>% filter(GID_0 == "ITA") %>% group_by(year) %>%
  summarise(nominal_lcu_bn = sum(grp_lcu, na.rm = TRUE) / 1e9, fx = median(fx), deflator = median(deflator_2015),
            version = names(which.max(table(version))))
cat("\nItaly, DOSE nominal regional GDP summed (grp_lcu, bn) with the DOSE exchange rate: nominal / fx is the series in US$ terms\n")
show(it_nom %>% filter(year >= YEAR0) %>% mutate(nominal_over_fx = nominal_lcu_bn / fx,
                                                  g_nominal = 100 * (nominal_lcu_bn / lag(nominal_lcu_bn) - 1)) %>% r2())
cat("\nItaly: growth gap on 100 dlog(fx), strict years, by regime\n")
for (rg in list(c(1970, 1984), c(1985, 1995), c(1996, 2014), c(2016, 2019), c(1970, 2014))) {
  x <- ita %>% filter(strict, year >= rg[1], year <= rg[2], is.finite(fx_dlog))
  f <- lm(gap ~ fx_dlog, data = x)
  cat(sprintf("  %d-%d: n = %2d, slope = %5.2f (se %.2f), R2 = %.2f; RMSE gap = %5.2f; RMSE gap - fx = %5.2f\n",
              rg[1], rg[2], nrow(x), coef(f)[2], summary(f)$coefficients[2, 2], summary(f)$r.squared,
              sqrt(mean(x$gap^2)), sqrt(mean((x$gap - x$fx_dlog)^2))))
}
cat("Italy 2015 regional growth (all regions):\n")
show(reg_growth %>% filter(GID_0 == "ITA", year == 2015) %>% mutate(g = 100 * (gdp / gdp_prev - 1)) %>%
       summarise(regions = n(), min = min(g), median = median(g), max = max(g)) %>% r2())
cat("Italy without its bad years (strict):\n")
show(tab %>% filter(iso3 == "ITA") %>% transmute(n, n_bad, clean = n_clean, corr, rmse, mean_gap, drift_sd,
                                                 med_abs = med_abs_gap, le2 = share_le2, spear = spearman) %>% r2())
for (a in c(1970L, 1995L, 2015L)) {
  r <- rstats(ita[ita$year >= a, ], "strict")
  cat("  from", a, ": n =", r$n, "bad =", r$n_bad, "clean =", r$n_clean, "corr =", round(r$corr, 2),
      "rmse =", round(r$rmse, 2), "drift_sd =", round(r$drift_sd, 3), "pass =", r$pass, "\n")
}

# Panel coverage ---------------------------------------------------------------

wld <- nat_all %>% filter(GID_0 == "WLD") %>% select(year, wld = gdp)
wld15 <- wld$wld[wld$year == 2015]
gdp_cy <- nat %>% transmute(iso3 = GID_0, year = as.integer(year), gdp)
rows_c <- panel_cy %>% group_by(GID_0) %>% summarise(panel_rows = sum(rows))
bad_rows <- bad_all %>% group_by(GID_0) %>% summarise(panel_rows_bad = sum(panel_rows))
tab <- tab %>% left_join(rows_c, by = c("iso3" = "GID_0")) %>% left_join(bad_rows, by = c("iso3" = "GID_0")) %>%
  left_join(gdp_cy %>% filter(year == 2015) %>% select(iso3, gdp15 = gdp), by = "iso3") %>%
  mutate(panel_rows = coalesce(panel_rows, 0L), panel_rows_bad = coalesce(panel_rows_bad, 0L),
         world_share_2015 = gdp15 / wld15)
# Country-years of the panel with WB GDP; set membership by country (and span).
pcy <- panel_cy %>% rename(iso3 = GID_0) %>% left_join(gdp_cy, by = c("iso3", "year")) %>%
  left_join(bad_all %>% transmute(iso3 = GID_0, year, bad = TRUE), by = c("iso3", "year")) %>%
  mutate(bad = coalesce(bad, FALSE)) %>%
  left_join(tab %>% select(iso3, status, span_first, span_last), by = "iso3")
wld_years <- sum(wld$wld[wld$year >= min(panel$year) & wld$year <= max(panel$year)])
share_set <- function(in_set, label) {
  s <- pcy[which(in_set), ]
  tibble(set = label, countries = n_distinct(s$iso3),
         rows_kept_bad = nrow(panel %>% semi_join(s, by = c("GID_0" = "iso3", "year"))),
         rows_dropped_bad = nrow(panel %>% semi_join(s %>% filter(!bad), by = c("GID_0" = "iso3", "year"))),
         share_rows_kept = rows_kept_bad / nrow(panel), share_rows_dropped = rows_dropped_bad / nrow(panel),
         world_gdp_2015 = sum(tab$world_share_2015[tab$iso3 %in% s$iso3], na.rm = TRUE),
         world_gdp_years_kept = sum(s$gdp, na.rm = TRUE) / wld_years,
         world_gdp_years_dropped = sum(s$gdp[!s$bad], na.rm = TRUE) / wld_years)
}
in_span <- pcy$status == "sub-period" & pcy$year >= pcy$span_first & pcy$year <= pcy$span_last
shares <- bind_rows(
  share_set(pcy$status == "pass", "robust set (pass)"),
  share_set(pcy$status %in% c("pass", "pass (matched)"), "pass + matched-only pass"),
  share_set(pcy$status %in% c("pass", "pass (matched)") | in_span, "pass + matched + sub-period spans"),
  share_set(rep(TRUE, nrow(pcy)), "all panel countries"))
prev_pass <- prev$iso3[prev$old_status == "pass"]
shares <- bind_rows(shares, share_set(pcy$iso3 %in% prev_pass, "previous pass list (whole period rule)"))
cat("\n== Panel coverage (dose$clean): rows, and WB GDP of panel country-years over world GDP 1970-2019 ==\n")
cat("kept = bad years kept; dropped = the set's bad years (strict, incl. skipped years seen by matched) removed\n")
show(shares %>% mutate(across(where(is.double), ~ round(.x, 3))))
rs <- bad_all %>% filter(in_robust_set)
cat("\nBad years inside the robust set and sub-period spans:", nrow(rs), "; with no dose$clean rows (already dropped, e.g. DOSE StructChange years):",
    sum(rs$panel_rows == 0), "(of which flagged StructChange:", sum(rs$panel_rows == 0 & !is.na(rs$sc_code)), ")",
    "; still in the panel:", sum(rs$panel_rows > 0), "country-years,", sum(rs$panel_rows), "rows\n")
cat("All bad years (all countries):", nrow(bad_all), "; still in dose$clean:", sum(bad_all$panel_rows > 0), "country-years,",
    sum(bad_all$panel_rows), "rows; StructChange-flagged among them:", sum(!is.na(bad_all$sc_code)), "\n")
cat("\nRows per country in the robust set and sub-period spans:\n")
show(tab %>% filter(status %in% c("pass", "pass (matched)", "sub-period")) %>%
       mutate(rows_span = if_else(status == "sub-period",
                                  purrr::pmap_int(list(iso3, span_first, span_last), function(i, a, b)
                                    sum(panel_cy$rows[panel_cy$GID_0 == i & panel_cy$year >= a & panel_cy$year <= b])),
                                  panel_rows)) %>%
       left_join(rs %>% group_by(GID_0) %>% summarise(rows_bad_in_set = sum(panel_rows),
                                                      bad_years_in_set = paste(year, collapse = " ")),
                 by = c("iso3" = "GID_0")) %>%
       transmute(iso3, status, span, panel_rows, rows_span, rows_bad_in_set = coalesce(rows_bad_in_set, 0L),
                 bad_years_in_set, world_share_2015) %>% r2() %>%
       mutate(world_share_2015 = round(world_share_2015, 4)))

# Sensitivity -------------------------------------------------------------------

base_set <- tab$iso3[tab$status == "pass"]
sens_one <- function(label, view, P) {
  p <- names(by_country)[vapply(by_country, function(d) rstats(d, view, P, full = FALSE)$pass, logical(1))]
  s <- share_set(pcy$iso3 %in% p, label)
  # Rows dropped: that variant's bad years.
  bad_v <- bind_rows(lapply(p, bad_years_of, view = view, P = P))
  kept_rows <- panel %>% filter(GID_0 %in% p)
  dropped_rows <- if (nrow(bad_v)) kept_rows %>% anti_join(bad_v, by = c("GID_0" = "iso3", "year")) else kept_rows
  tibble(variant = label, view = view, n_pass = length(p),
         added = paste(setdiff(p, base_set), collapse = " "), dropped = paste(setdiff(base_set, p), collapse = " "),
         share_rows_kept = round(nrow(kept_rows) / nrow(panel), 3),
         share_rows_dropped = round(nrow(dropped_rows) / nrow(panel), 3),
         world_gdp_2015 = round(s$world_gdp_2015, 3), passing = paste(p, collapse = " "))
}
grid <- list()
for (v in c("strict", "matched")) for (tpp in c(3, 5, 8)) {
  for (k in 0:3) grid[[length(grid) + 1]] <- sens_one(paste0("T=", tpp, " k=", k), v, modifyList(PAR, list(bad_pp = tpp, k = k)))
  grid[[length(grid) + 1]] <- sens_one(paste0("T=", tpp, " k=10%"), v, modifyList(PAR, list(bad_pp = tpp, k_share = 0.10)))
}
grid[[length(grid) + 1]] <- sens_one("robust z > 5 and |gap| > 2, k=2", "strict", modifyList(PAR, list(bad_rule = "z")))
old_1970 <- names(by_country)[vapply(by_country, function(d) isTRUE(rstats(d, "strict")$old_pass), logical(1))]
sens <- bind_rows(grid) %>%
  bind_rows(tibble(variant = "old rule, 1970+ (no allowance, SD of log level)", view = "strict", n_pass = length(old_1970),
                   added = paste(setdiff(old_1970, base_set), collapse = " "),
                   dropped = paste(setdiff(base_set, old_1970), collapse = " "),
                   share_rows_kept = round(sum(tab$panel_rows[tab$iso3 %in% old_1970]) / nrow(panel), 3),
                   passing = paste(old_1970, collapse = " ")),
            tibble(variant = "old rule, whole DOSE period (dose_quality_all.R)", view = "strict", n_pass = length(prev_pass),
                   added = paste(setdiff(prev_pass, base_set), collapse = " "),
                   dropped = paste(setdiff(base_set, prev_pass), collapse = " "),
                   share_rows_kept = round(sum(tab$panel_rows[tab$iso3 %in% prev_pass]) / nrow(panel), 3),
                   passing = paste(prev_pass, collapse = " ")))
cat("\n== Sensitivity: bad-year threshold T (pp), allowance k, strict vs matched; added/dropped vs baseline (strict, T=5, k=2) ==\n")
show(sens %>% select(-passing))
cat("\nBaseline set:", paste(base_set, collapse = " "), "\n")

# Save ---------------------------------------------------------------------------

write.csv(tab, file.path(out_dir, "country_quality_robust.csv"), row.names = FALSE)
write.csv(bad_all %>% filter(in_robust_set) %>%
            select(GID_0, year, dose_growth, wb_growth, gap, primary_driver_region_id, primary_driver_region_name,
                   primary_driver_contribution_pp, primary_driver_share_of_discrepancy, primary_driver_regional_growth_pct,
                   reason, comparison, coverage, uniform_share, sc_code, version_change, fx_dlog, driver_source,
                   status, panel_rows),
          file.path(out_dir, "bad_years_robust_set.csv"), row.names = FALSE)
write.csv(bad_all, file.path(out_dir, "bad_years_all_countries.csv"), row.names = FALSE)
write.csv(cy, file.path(out_dir, "country_year_comparison_1970.csv"), row.names = FALSE)
write.csv(sens, file.path(out_dir, "rule_sensitivity_robust.csv"), row.names = FALSE)
write.csv(shares, file.path(out_dir, "panel_coverage.csv"), row.names = FALSE)
sink()
