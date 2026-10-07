# National armed conflict from the UCDP/PRIO Armed Conflict Dataset v26.1
# (data/UcdpPrioConflict_v26_1.csv: one row per conflict and year, 1946-2025).
#
# A country-year is in conflict (conflict = 1) when the country is a location
# (gwno_loc) of at least one state-based armed conflict that year, i.e. one
# with at least 25 battle-related deaths, of any type (extrasystemic,
# interstate, intrastate, internationalised intrastate). Extrasystemic
# (colonial) conflicts are located in the colony, not the colonial power;
# interstate conflicts list every primary party as a location (e.g. the United
# States for Grenada 1983). conflict_intensity is the highest intensity that
# year: 1 minor (25-999 deaths), 2 war (1,000 or more); conflict_low and
# conflict_high are its two indicators (highest intensity minor, or war).
#
# Gleditsch-Ward codes are mapped to ISO3 with countrycode. Codes without an
# ISO3 equivalent go to the present state of their territory (GW_TO_ISO3).
# Conflicts that UCDP locates in the Soviet Union or Yugoslavia but that were
# fought on the territory of a present-day successor state go to that state,
# by their territory_name (PREDECESSOR_TERRITORY); other conflicts in those two
# states stay with Russia and Serbia. conflict_ucdp is this location-based
# indicator.
#
# conflict adds two kinds of country-years that UCDP does not list as
# locations, each conflict reviewed by hand from its territory_name and
# warring parties:
#   CONFLICT_TERRITORY  conflicts fought on the territory of another
#                       present-day state (e.g. the Israel-Palestinian conflict,
#                       located in Israel, in Palestine), every year;
#   CONFLICT_PARTY      states that fought, as a primary or secondary party,
#                       over a contested territory next to them (e.g. Armenia
#                       in Nagorno-Karabakh), in the years given (NA: every
#                       year of the conflict).
# Not added: states contributing to distant coalitions and peace operations,
# and neighbours intervening in conflicts over government.
GW_TO_ISO3 <- c(
  `55` = "GRD",   # Grenada
  `345` = "SRB",  # Yugoslavia / Serbia
  `678` = "YEM",  # North Yemen
  `680` = "YEM",  # South Yemen
  `751` = "IND",  # Hyderabad
  `816` = "VNM"   # North Vietnam
)
PREDECESSOR_TERRITORY <- list(
  `365` = c(Azerbaijan = "AZE", `Nagorno-Karabakh` = "AZE", Estonia = "EST",
            Latvia = "LVA", Lithuania = "LTU", Ukraine = "UKR"),
  `345` = c(Croatia = "HRV", Slovenia = "SVN")
)
CONFLICT_TERRITORY <- tibble::tribble(
  ~conflict_id, ~GID_0, ~note,
  228L,   "PSE", "Palestine 1948-49, interstate",
  234L,   "PSE", "Israel-Palestinian conflict, located in Israel",
  301L,   "PSE", "West Bank 1967, located in Israel and Jordan",
  426L,   "LBN", "Hezbollah in Southern Lebanon, located in Israel",
  312L,   "BGD", "East Pakistan 1971, located in Pakistan",
  298L,   "NAM", "Namibia, located in South Africa",
  275L,   "ERI", "Eritrea, located in Ethiopia",
  282L,   "SSD", "Southern Sudan, located in Sudan",
  11344L, "SSD", "Abyei 2011, contested by Sudan and South Sudan",
  330L,   "TLS", "East Timor, located in Indonesia",
  331L,   "ESH", "Western Sahara, located in Morocco",
  428L,   "ESH", "Western Sahara, located in Mauritania",
  412L,   "XKO", "Kosovo 1998-99, located in Serbia",
  217L,   "TWN", "Taiwan 1947, located in China",
  238L,   "PRI", "Puerto Rico 1950, located in the United States"
)
CONFLICT_PARTY <- tibble::tribble(
  ~conflict_id, ~GID_0, ~first, ~last, ~note,
  376L,   "ARM", NA,    NA,    "Nagorno-Karabakh 1990-91, Armenia the opposing side",
  388L,   "ARM", NA,    NA,    "Nagorno-Karabakh (Artsakh), Armenia a secondary party",
  384L,   "SRB", NA,    NA,    "Slovenia 1991, Yugoslavia the government side",
  385L,   "SRB", NA,    NA,    "Croatia 1991-92, Yugoslavia the government side",
  389L,   "SRB", 1992L, 1993L, "Bosnian Serbs, Yugoslavia a secondary party",
  390L,   "SRB", 1992L, 1993L, "Croatian Serbs, Yugoslavia a secondary party",
  398L,   "HRV", 1993L, 1993L, "Bosnian Croats, Croatia a secondary party",
  393L,   "RUS", 2008L, 2008L, "South Ossetia, Russia a secondary party",
  13246L, "RUS", NA,    NA,    "Donetsk 2014, Russia a secondary party",
  13247L, "RUS", NA,    NA,    "Lugansk 2014, Russia a secondary party",
  13306L, "RUS", NA,    NA,    "Novorossiya 2014-22, Russia a secondary party",
  331L,   "DZA", 1976L, 1976L, "Western Sahara, Algeria a secondary party"
)
UCDP_YEARS <- c(1946L, 2025L)

load_conflict <- function() {
  raw <- utils::read.csv(data_file("UcdpPrioConflict_v26_1.csv"),
                         stringsAsFactors = FALSE)
  # One row per conflict, year and location.
  locations <- strsplit(as.character(raw$gwno_loc), ",", fixed = TRUE)
  long <- tibble::tibble(
    row = rep(seq_len(nrow(raw)), lengths(locations)),
    gw = as.integer(trimws(unlist(locations)))
  ) %>%
    mutate(
      year = as.integer(raw$year[row]),
      intensity = as.integer(raw$intensity_level[row]),
      territory = raw$territory_name[row],
      single_location = lengths(locations)[row] == 1L
    )
  if (anyNA(long$gw)) stop("Unparsed gwno_loc codes in the UCDP file.")

  long$GID_0 <- countrycode::countrycode(long$gw, "gwn", "iso3c", warn = FALSE)
  fixed <- as.character(long$gw) %in% names(GW_TO_ISO3)
  long$GID_0[fixed] <- GW_TO_ISO3[as.character(long$gw[fixed])]
  for (g in names(PREDECESSOR_TERRITORY)) {
    map <- PREDECESSOR_TERRITORY[[g]]
    hit <- long$gw == as.integer(g) & long$single_location & long$territory %in% names(map)
    long$GID_0[hit] <- map[long$territory[hit]]
  }
  if (anyNA(long$GID_0)) {
    stop("No ISO3 code for Gleditsch-Ward codes: ",
         paste(sort(unique(long$gw[is.na(long$GID_0)])), collapse = ", "))
  }
  located <- transmute(long, GID_0, year, intensity, ucdp = TRUE)

  # The hand-reviewed additions, with the intensity of the conflict that year.
  years <- tibble::tibble(conflict_id = raw$conflict_id, year = as.integer(raw$year),
                          intensity = as.integer(raw$intensity_level))
  missing_ids <- setdiff(c(CONFLICT_TERRITORY$conflict_id, CONFLICT_PARTY$conflict_id),
                         years$conflict_id)
  if (length(missing_ids)) {
    stop("Conflict ids not in the UCDP file: ", paste(missing_ids, collapse = ", "))
  }
  territory <- inner_join(select(CONFLICT_TERRITORY, conflict_id, GID_0), years,
                          by = "conflict_id", relationship = "many-to-many")
  party <- inner_join(select(CONFLICT_PARTY, conflict_id, GID_0, first, last), years,
                      by = "conflict_id", relationship = "many-to-many") %>%
    filter(is.na(first) | (year >= first & year <= last))
  added <- bind_rows(territory, party) %>%
    transmute(GID_0, year, intensity, ucdp = FALSE)

  bind_rows(located, added) %>%
    group_by(GID_0, year) %>%
    summarise(conflict_intensity = max(intensity), conflict_ucdp = as.integer(any(ucdp)),
              .groups = "drop") %>%
    mutate(conflict = 1L, .after = year) %>%
    mutate(conflict_low = as.integer(conflict_intensity == 1L),
           conflict_high = as.integer(conflict_intensity == 2L))
}

# Adds conflict, conflict_ucdp, conflict_intensity, conflict_low and
# conflict_high to a panel by GID_0 and year. UCDP records every state-based
# armed conflict, so a country-year within its coverage that it does not list
# had none (0); years outside the coverage are NA.
add_conflict <- function(panel, conflict = load_conflict()) {
  covered <- panel$year >= UCDP_YEARS[[1]] & panel$year <= UCDP_YEARS[[2]]
  out <- left_join(panel, conflict, by = c("GID_0", "year"))
  for (v in c("conflict", "conflict_ucdp", "conflict_intensity", "conflict_low",
              "conflict_high")) {
    out[[v]] <- ifelse(covered, dplyr::coalesce(out[[v]], 0L), NA_integer_)
  }
  out
}

# Phases of war and of armed conflict, from a country's status this year and
# the year before (both from the full UCDP/PRIO series, so gaps in the panel do
# not create false onsets):
#   onset         in war this year, not the year before (the first year)
#   continuation  in war this year and the year before
#   end           not in war this year, in war the year before (the first year
#                 without war)
#   none          in war neither year
# War is conflict_high (a conflict with 1,000 or more battle deaths that year);
# conflict is any armed conflict (25 or more, wars included), so a war that
# winds down to a minor conflict is a war end and a conflict continuation.
# Adds war_phase and conflict_phase (factors, "none" first) and the indicators
# war_onset, war_continuation, war_end and conflict_onset, conflict_continuation,
# conflict_end. NA when the year or the year before is outside UCDP_YEARS.
CONFLICT_PHASES <- c("none", "onset", "continuation", "end")

add_conflict_phases <- function(panel, conflict = load_conflict()) {
  status <- function(column, year) {
    listed <- paste(conflict$GID_0, conflict$year)[conflict[[column]] == 1L]
    ifelse(in_period(year, UCDP_YEARS),
           as.integer(paste(panel$GID_0, year) %in% listed), NA_integer_)
  }
  for (prefix in c("war", "conflict")) {
    column <- c(war = "conflict_high", conflict = "conflict")[[prefix]]
    now <- status(column, panel$year)
    before <- status(column, panel$year - 1L)
    phase <- dplyr::case_when(
      now == 1L & before == 0L ~ "onset",
      now == 1L & before == 1L ~ "continuation",
      now == 0L & before == 1L ~ "end",
      now == 0L & before == 0L ~ "none"
    )
    panel[[paste0(prefix, "_phase")]] <- factor(phase, levels = CONFLICT_PHASES)
    for (p in CONFLICT_PHASES[-1]) {
      panel[[paste0(prefix, "_", p)]] <- as.integer(phase == p)
    }
  }
  panel
}
