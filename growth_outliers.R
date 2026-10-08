# Growth outliers of the national panel and the events they coincide with.
#
# Panel: PWT 11.0 x ERA5 (population-density weights, base year 2000), country-
# years from 1960 with growth of GDP per capita and both weather variables (the
# sample of weather_arms_robustness.Rmd before any outlier is dropped).
#
# Outliers: the country-years whose growth of GDP per capita is below the 5th or
# above the 95th percentile of all the country-years of the panel (10% in
# total).
#
# Each outlier is matched with the UCDP/PRIO armed conflicts located in the
# country that year and the year before (functions/conflict.R: the locations of
# UCDP, the successor states of the Soviet Union and Yugoslavia, and the
# hand-reviewed territories and parties), with the post-Soviet and
# post-Yugoslav collapse years, and with the 1991-1993 outlier rule used so far
# (flag_growth_outliers()).
#
# Each outlier also has a likely cause, read from data/growth_outlier_causes.csv
# (reviewed by hand from the conflict matches and from economic history): a
# specific event, a category, a confidence (high: well-documented event of that
# year; medium: plausible main cause; low: uncertain, unverified, or no event
# identified) and whether weather (drought, storm, flood, harvest) is part of the
# cause. Weather-driven outliers are the response studied, not noise.
#
# Output: results/growth_outliers/ (outliers.csv, thresholds.csv, categories.csv)

source("load_functions.R")

START_YEAR <- 1960L
SHARE <- 0.10
OUTPUT_DIR <- file.path("results", "growth_outliers")
CONTROLS <- c("conflict_low", "conflict_high", "i(fsu_collapse_year, ref = 'none')",
              "i(yug_collapse_year, ref = 'none')")
FIXED_EFFECTS <- "year + GID_1 + GID_1[year]"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

# Panel -----------------------------------------------------------------------------

cfg <- panel_config(
  econ_year_min = START_YEAR,
  climate_source = "ERA5",
  climate_weight = "population density",
  climate_weight_year = 2000L,
  clim_history_years = 1L
)
d <- suppressMessages(build_dat("PWT", config = cfg))
d <- add_yugoslav_collapse(add_soviet_collapse(add_conflict(d)))
d <- d[est_row(d) & d$year >= START_YEAR & !is.na(d$TM) & !is.na(d$RR) &
         !is.na(d$conflict), ]

# Outliers ---------------------------------------------------------------------------

cutoffs <- stats::quantile(d$dlgrp_pc_usd, c(SHARE / 2, 1 - SHARE / 2), names = FALSE)
d$tail <- ifelse(d$dlgrp_pc_usd < cutoffs[[1]], "low",
                 ifelse(d$dlgrp_pc_usd > cutoffs[[2]], "high", NA_character_))

# The 1991-1993 rule used so far: the 2.5% largest growth shocks net of the fixed
# effects and the controls, in 1991-1993.
old_rule <- flag_growth_outliers(d, CONTROLS, FIXED_EFFECTS, 1991:1993, 0.025)
d$old_rule <- old_rule$flag

thresholds <- tibble::tibble(
  country_years = nrow(d), countries = n_distinct(d$GID_0),
  first_year = min(d$year), last_year = max(d$year),
  low_cutoff_pct = 100 * cutoffs[[1]], high_cutoff_pct = 100 * cutoffs[[2]],
  low_outliers = sum(d$tail == "low", na.rm = TRUE),
  high_outliers = sum(d$tail == "high", na.rm = TRUE),
  growth_median_pct = 100 * stats::median(d$dlgrp_pc_usd)
)

# Conflicts by country and year, with their names ---------------------------------

# The UCDP locations as in load_conflict(), one row per conflict, year and
# country, plus the hand-reviewed territories and parties.
conflict_events <- function() {
  raw <- utils::read.csv(data_file("UcdpPrioConflict_v26_1.csv"), stringsAsFactors = FALSE)
  raw$label <- sprintf("%s vs %s (%s)", raw$side_a, raw$side_b,
                       c("minor", "war")[as.integer(raw$intensity_level)])
  locations <- strsplit(as.character(raw$gwno_loc), ",", fixed = TRUE)
  long <- tibble::tibble(
    row = rep(seq_len(nrow(raw)), lengths(locations)),
    gw = as.integer(trimws(unlist(locations)))
  ) %>%
    mutate(year = as.integer(raw$year[row]), intensity = as.integer(raw$intensity_level[row]),
           territory = raw$territory_name[row], label = raw$label[row],
           single_location = lengths(locations)[row] == 1L)
  long$GID_0 <- countrycode::countrycode(long$gw, "gwn", "iso3c", warn = FALSE)
  fixed <- as.character(long$gw) %in% names(GW_TO_ISO3)
  long$GID_0[fixed] <- GW_TO_ISO3[as.character(long$gw[fixed])]
  for (g in names(PREDECESSOR_TERRITORY)) {
    map <- PREDECESSOR_TERRITORY[[g]]
    hit <- long$gw == as.integer(g) & long$single_location & long$territory %in% names(map)
    long$GID_0[hit] <- map[long$territory[hit]]
  }
  located <- transmute(long, GID_0, year, intensity, label)

  by_id <- tibble::tibble(conflict_id = raw$conflict_id, year = as.integer(raw$year),
                          intensity = as.integer(raw$intensity_level), label = raw$label)
  territory <- inner_join(select(CONFLICT_TERRITORY, conflict_id, GID_0), by_id,
                          by = "conflict_id", relationship = "many-to-many")
  party <- inner_join(select(CONFLICT_PARTY, conflict_id, GID_0, first, last), by_id,
                      by = "conflict_id", relationship = "many-to-many") %>%
    filter(is.na(first) | (year >= first & year <= last))
  added <- bind_rows(territory, party) %>%
    transmute(GID_0, year, intensity, label = paste(label, "[added]"))

  bind_rows(located, added) %>%
    filter(!is.na(GID_0)) %>%
    group_by(GID_0, year) %>%
    summarise(conflicts = paste(unique(label), collapse = "; "), .groups = "drop")
}
events <- conflict_events()

# Outlier table ------------------------------------------------------------------------

outliers <- d %>%
  filter(!is.na(tail)) %>%
  transmute(
    GID_0, country = countrycode::countrycode(GID_0, "iso3c", "country.name", warn = FALSE),
    year, tail, growth_pct = round(100 * dlgrp_pc_usd, 1),
    conflict = c("none", "minor", "war")[conflict_intensity + 1],
    collapse = ifelse(fsu_collapse_year != "none", paste0("post-Soviet, year ", fsu_collapse_year),
                      ifelse(yug_collapse_year != "none",
                             paste0("post-Yugoslav, year ", yug_collapse_year), "")),
    dropped_by_1991_93_rule = old_rule
  ) %>%
  left_join(rename(events, conflicts_that_year = conflicts), by = c("GID_0", "year")) %>%
  left_join(transmute(events, GID_0, year = year + 1L, conflicts_year_before = conflicts),
            by = c("GID_0", "year")) %>%
  mutate(
    across(c(conflicts_that_year, conflicts_year_before), ~ coalesce(.x, "")),
    matched = case_when(
      nzchar(conflicts_that_year) ~ "conflict that year",
      nzchar(collapse) ~ "collapse",
      nzchar(conflicts_year_before) ~ "conflict the year before",
      TRUE ~ ""
    )
  ) %>%
  arrange(growth_pct)

# Causes ------------------------------------------------------------------------------

causes <- utils::read.csv(data_file("growth_outlier_causes.csv"), stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(causes[, c("GID_0", "year")]))
unannotated <- anti_join(outliers, causes, by = c("GID_0", "year"))
not_outliers <- anti_join(causes, outliers, by = c("GID_0", "year"))
if (nrow(unannotated) || nrow(not_outliers)) {
  stop("Outliers without a cause: ", paste(unannotated$GID_0, unannotated$year, collapse = ", "),
       "; causes that are not outliers: ",
       paste(not_outliers$GID_0, not_outliers$year, collapse = ", "))
}
outliers <- left_join(outliers, causes, by = c("GID_0", "year")) %>%
  relocate(cause, category, confidence, weather, .after = growth_pct)

categories <- outliers %>%
  count(category, tail) %>%
  tidyr::pivot_wider(names_from = tail, values_from = n, values_fill = 0) %>%
  mutate(total = low + high) %>%
  arrange(desc(total))

write_out(thresholds, "thresholds")
write_out(outliers, "outliers")
write_out(categories, "categories")

options(width = 200)
print(as.data.frame(thresholds), row.names = FALSE)
cat("\nOutliers by tail and conflict match:\n")
print(as.data.frame(outliers %>% count(tail, matched) %>%
  tidyr::pivot_wider(names_from = tail, values_from = n, values_fill = 0)), row.names = FALSE)
cat("\nOutliers by category of cause:\n")
print(as.data.frame(categories), row.names = FALSE)
cat("\nConfidence:\n")
print(table(outliers$confidence, outliers$tail))
cat("\nWeather part of the cause:\n")
print(table(outliers$weather, outliers$tail))
message("Results written to ", normalizePath(OUTPUT_DIR))
