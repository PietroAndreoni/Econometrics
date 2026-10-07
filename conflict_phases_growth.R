# GDP growth through the phases of wars and armed conflicts, with the growth of
# capital and employment (PWT 11.0 and UCDP/PRIO v26.1).
#
#   dlog Y = g1 onset + g2 continuation + g3 end + a dlog K + b dlog L
#            + country + year + country linear trend
#
# Y real GDP and K the capital stock at constant national prices (PWT rgdpna,
# rnna), L persons engaged (emp); dlog the change in logs from the previous
# calendar year. Phases (add_conflict_phases() in functions/conflict.R), each
# against the years in which the country was in war (conflict) neither that
# year nor the year before:
#   onset         the first year of war: war this year, not the year before
#   continuation  war this year and the year before
#   end           the first year without war: no war this year, war the year
#                 before
# War: a UCDP/PRIO armed conflict with 1,000 or more battle deaths that year on
# the country's territory (conflict_high); conflict: any armed conflict, 25 or
# more deaths, wars included (conflict).
#
# Main models, all on the same rows (every PWT country-year with the growth of
# Y, K and L):
#   War, Conflict            the phases alone (the total effect)
#   War + dK dL, ...         with the growth of capital and employment (the
#                            effect given the factors)
#   War & conflict + dK dL   both sets: the conflict phases are then the
#                            phases of conflicts that are not wars
# Robustness of the models with dK and dL: GDP, capital and employment per
# capita (pop); 1990-2023; the post-Soviet and post-Yugoslav collapse-year
# dummies (functions/transition_collapse.R).
#
# Errors clustered by country.
#
# Output: results/conflict_phases_growth/ (estimates.csv, phase_counts.csv)

source("load_functions.R")

OUTPUT_DIR <- file.path("results", "conflict_phases_growth")
FIXED_EFFECTS <- "year + GID_0 + GID_0[year]"
FACTORS <- c("dlk", "dll")
PHASES <- list(
  war = c("war_onset", "war_continuation", "war_end"),
  conflict = c("conflict_onset", "conflict_continuation", "conflict_end")
)
COLLAPSE <- c("i(fsu_collapse_year, ref = 'none')", "i(yug_collapse_year, ref = 'none')")
RECENT_YEARS <- c(1990L, 2023L)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

# Growth of GDP, capital and employment, in total and per capita -----------------

pwt <- readxl::read_excel(data_file(ECON_RAW_FILES[["PWT110"]]), sheet = "Data") %>%
  transmute(GID_0 = as.character(countrycode), year = as.integer(year),
            rgdpna = as.numeric(rgdpna), rnna = as.numeric(rnna),
            emp = as.numeric(emp), pop = as.numeric(pop)) %>%
  filter(rgdpna > 0, rnna > 0, emp > 0, pop > 0) %>%
  arrange(GID_0, year) %>%
  group_by(GID_0) %>%
  mutate(
    dly = log(rgdpna) - lag_by_year(log(rgdpna), year),
    dlk = log(rnna) - lag_by_year(log(rnna), year),
    dll = log(emp) - lag_by_year(log(emp), year),
    dln = log(pop) - lag_by_year(log(pop), year),
    dly_pc = dly - dln,
    dlk_pc = dlk - dln,
    dll_pc = dll - dln
  ) %>%
  ungroup() %>%
  filter(!is.na(dly), !is.na(dlk), !is.na(dll))

conflict <- load_conflict()
data <- pwt %>%
  add_conflict(conflict) %>%
  add_conflict_phases(conflict) %>%
  add_soviet_collapse() %>%
  add_yugoslav_collapse()
stopifnot(!anyNA(data$war_phase), !anyNA(data$conflict_phase))

phase_counts <- bind_rows(lapply(c("war", "conflict"), function(kind) {
  data %>%
    count(phase = .data[[paste0(kind, "_phase")]], name = "country_years") %>%
    left_join(data %>%
                group_by(phase = .data[[paste0(kind, "_phase")]]) %>%
                summarise(countries = n_distinct(GID_0),
                          mean_growth_pct = 100 * mean(dly), .groups = "drop"),
              by = "phase") %>%
    mutate(kind = kind, .before = 1)
}))
write_out(phase_counts, "phase_counts")

# Models -------------------------------------------------------------------------

fit <- function(terms, model_data, outcome = "dly") {
  fit_panel_model(paste(terms, collapse = " + "), model_data, outcome = outcome,
                  fixed_effects = FIXED_EFFECTS, panel_id = c("GID_0", "year"),
                  cluster = ~GID_0)
}

main <- list(
  "War" = fit(PHASES$war, data),
  "War + dK dL" = fit(c(PHASES$war, FACTORS), data),
  "Conflict" = fit(PHASES$conflict, data),
  "Conflict + dK dL" = fit(c(PHASES$conflict, FACTORS), data),
  "War & conflict + dK dL" = fit(c(PHASES$war, PHASES$conflict, FACTORS), data)
)
stopifnot(length(unique(vapply(main, stats::nobs, 0))) == 1)

recent <- filter(data, in_period(year, RECENT_YEARS))
robustness <- list()
for (kind in names(PHASES)) {
  label <- c(war = "War", conflict = "Conflict")[[kind]]
  robustness[[paste(label, "main")]] <- main[[paste(label, "+ dK dL")]]
  robustness[[paste(label, "per capita")]] <- fit(c(PHASES[[kind]], "dlk_pc", "dll_pc"),
                                                  data, outcome = "dly_pc")
  robustness[[paste(label, "1990-2023")]] <- fit(c(PHASES[[kind]], FACTORS), recent)
  robustness[[paste(label, "collapse")]] <- fit(c(PHASES[[kind]], FACTORS, COLLAPSE), data)
}

tidy <- function(models, set) {
  bind_rows(lapply(names(models), function(k) {
    ct <- fixest::coeftable(models[[k]])
    tibble::tibble(set = set, model = k, term = rownames(ct), estimate = ct[, 1],
                   std_error = ct[, 2], p_value = ct[, 4],
                   observations = stats::nobs(models[[k]]),
                   countries = length(unique(models[[k]]$fixef_id$GID_0)))
  }))
}
estimates <- bind_rows(tidy(main, "main"), tidy(robustness, "robustness")) %>%
  filter(!grepl("collapse_year", term))
write_out(estimates, "estimates")

# Report -------------------------------------------------------------------------

options(width = 200)
cat("\n== Country-years by phase (PWT rows with growth of GDP, capital and employment,",
    paste(range(data$year), collapse = "-"), ") ==\n")
print(as.data.frame(phase_counts %>% mutate(mean_growth_pct = round(mean_growth_pct, 2))),
      row.names = FALSE)

DICT <- c(dly = "dlog GDP", dlk = "dlog capital", dll = "dlog employment",
          war_onset = "War: onset", war_continuation = "War: continuation",
          war_end = "War: first year without",
          conflict_onset = "Conflict: onset", conflict_continuation = "Conflict: continuation",
          conflict_end = "Conflict: first year without")
cat("\n== Main models: growth of GDP (log points) ==\n")
print(etable(main, dict = DICT, digits = 3, fitstat = ~ n + r2 + wr2,
             signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1)))

cat("\n== Robustness of the models with dK and dL: phase coefficients, percentage points",
    "(clustered SE) ==\n")
stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
print(as.data.frame(estimates %>%
  filter(set == "robustness", grepl("_(onset|continuation|end)$", term)) %>%
  transmute(model, phase = sub("^(war|conflict)_", "", term),
            cell = sprintf("%.2f (%.2f)%s", 100 * estimate, 100 * std_error, stars(p_value)),
            observations) %>%
  tidyr::pivot_wider(names_from = phase, values_from = cell)), row.names = FALSE)

message("Results written to ", normalizePath(OUTPUT_DIR))
