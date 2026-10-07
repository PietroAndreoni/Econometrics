# Expected climate over time for selected countries: the trailing 30-year mean
# (years t-30 to t-1, the mean the anomalies use) against Hamilton trends
# fitted in real time, each year's forecast from a regression trained only on
# the 30 years up to the forecast origin t-2 (see hamilton_trend()). Two
# Hamilton settings: year t from years t-2 to t-5 (h = 2, four lags, the
# panel default) and from years t-2 and t-3 (Bilal and Kanzig: their p = 1
# counts the lags beyond the first; hamilton_trend() counts all lags). The
# default countries include Peru, whose temperature swings with El Nino.
#
# Run options (all optional, in order): countries as comma-separated ISO3
# codes, variable (TM or RR), climate sources (comma-separated), weight.
#   Rscript plot_hamilton_trends.R
#   Rscript plot_hamilton_trends.R USA,IND,PER,IDN TM "ERA5,CRU TS"
#
# A second figure splits each year's departure from its 1951-1980 mean into
# the climate (trailing 30-year mean), the gap and the Hamilton shock (Bilal
# and Kanzig's trend), in SDs of the trailing 30 years.
#
# Output: results/hamilton_trends/ (hamilton_trends_<variable>_<countries>.png,
# hamilton_decomposition_<variable>_<countries>.png)

source("load_functions.R")

args <- commandArgs(trailingOnly = TRUE)
COUNTRIES <- if (length(args) >= 1) strsplit(toupper(args[[1]]), ",")[[1]] else
  c("USA", "IND", "PER")
VARIABLE <- if (length(args) >= 2) args[[2]] else "TM"
CLIMATE_SOURCES <- if (length(args) >= 3) trimws(strsplit(args[[3]], ",")[[1]]) else
  c("ERA5", "CRU TS", "UDelaware")
CLIMATE_WEIGHT <- if (length(args) >= 4) args[[4]] else "area"
HAMILTON_WINDOW <- 30L
HAMILTON_H <- 2L
# Number of lags of each Hamilton setting.
HAMILTON_LAGS <- c(hamilton_p4 = 4L, hamilton_bk = 2L)
stopifnot(VARIABLE %in% c("TM", "RR"))

UNITS <- c(TM = "°C", RR = "m")
VARIABLE_LABELS <- c(TM = "Temperature", RR = "Precipitation")
SERIES_LABELS <- c(
  annual = "Annual value",
  mean_lag_30 = "Trailing 30-year mean (t-30 to t-1)",
  hamilton_p4 = "Hamilton trend from years t-2 to t-5",
  hamilton_bk = "Hamilton trend from years t-2 and t-3 (Bilal & Känzig)"
)

OUTPUT_DIR <- file.path("results", "hamilton_trends")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

wide <- bind_rows(lapply(CLIMATE_SOURCES, function(source) {
  cfg <- panel_config(climate_source = source, climate_weight = CLIMATE_WEIGHT)
  series <- suppressMessages(load_climate_panel("gadm0", cfg)) %>%
    filter(GID_0 %in% COUNTRIES) %>%
    select(GID_0, year, annual = all_of(VARIABLE),
           mean_30 = all_of(paste0("mean_", VARIABLE, "_30")),
           sd_30 = all_of(paste0("sd_", VARIABLE, "_30")))
  missing <- setdiff(COUNTRIES, series$GID_0)
  if (length(missing)) stop("Not in the ", source, " panel: ", paste(missing, collapse = ", "))
  series %>%
    group_by(GID_0) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      mean_lag_30 = dplyr::lag(mean_30),
      sd_lag_30 = dplyr::lag(sd_30),
      hamilton_p4 = hamilton_trend(annual, h = HAMILTON_H, p = HAMILTON_LAGS[["hamilton_p4"]],
                                   window = HAMILTON_WINDOW),
      hamilton_bk = hamilton_trend(annual, h = HAMILTON_H, p = HAMILTON_LAGS[["hamilton_bk"]],
                                   window = HAMILTON_WINDOW)
    ) %>%
    ungroup() %>%
    select(GID_0, year, sd_lag_30, all_of(names(SERIES_LABELS))) %>%
    mutate(source = source)
}))

d <- wide %>%
  select(-sd_lag_30) %>%
  tidyr::pivot_longer(-c(GID_0, year, source), names_to = "series", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(series = factor(SERIES_LABELS[series], levels = unname(SERIES_LABELS)),
         GID_0 = factor(GID_0, levels = COUNTRIES),
         source = factor(source, levels = CLIMATE_SOURCES))

# Purple for the two-lag trend, distinct from the orange four-lag one.
colours <- stats::setNames(c(GRID_MUTED, GRID_PALETTE[c(1, 2, 7)]), unname(SERIES_LABELS))
widths <- stats::setNames(c(0.35, 0.9, 0.8, 0.8), unname(SERIES_LABELS))

# Rows: climate sources; columns: countries; each panel its own y scale.
p <- ggplot(d, aes(x = year, y = value, colour = series, linewidth = series)) +
  geom_line() +
  facet_wrap(vars(source, GID_0), ncol = length(COUNTRIES), scales = "free_y",
             labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_linewidth_manual(values = widths, name = NULL) +
  guides(colour = guide_legend(nrow = 2), linewidth = guide_legend(nrow = 2)) +
  labs(
    x = NULL, y = paste0(VARIABLE_LABELS[[VARIABLE]], " (", UNITS[[VARIABLE]], ")"),
    title = paste0("Expected ", tolower(VARIABLE_LABELS[[VARIABLE]]),
                   ": trailing mean and real-time Hamilton trends, by climate source"),
    subtitle = paste("Hamilton trends in real time: each year's forecast from a regression",
                     "fitted only on the", HAMILTON_WINDOW, "years up to t-2, so they start",
                     "once a full window\nexists (1976 with ERA5, which starts in 1940)."),
    caption = paste0("Country series weighted by ", CLIMATE_WEIGHT, ".")
  ) +
  theme_results()

file <- paste0("hamilton_trends_", VARIABLE, "_", paste(COUNTRIES, collapse = "_"), ".png")
ggsave(file.path(OUTPUT_DIR, file), p, width = 4.5 * length(COUNTRIES) + 1,
       height = 3.2 * length(CLIMATE_SOURCES) + 2, dpi = 300)
message("Written ", normalizePath(file.path(OUTPUT_DIR, file)))


# Second figure: each year's departure from the reference climate (the mean of
# REFERENCE_PERIOD) split into three parts, all signed and in SDs of the
# trailing 30 years (t-30 to t-1):
#   climate  the trailing 30-year mean minus the reference
#   gap      Bilal and Kanzig's real-time Hamilton trend minus the trailing mean
#   shock    the year minus that Hamilton trend
# so total = climate + gap + shock.
REFERENCE_PERIOD <- c(1951L, 1980L)
decomposition <- wide %>%
  group_by(source, GID_0) %>%
  mutate(reference = mean(annual[year >= REFERENCE_PERIOD[[1]] &
                                   year <= REFERENCE_PERIOD[[2]]])) %>%
  ungroup() %>%
  transmute(
    GID_0, year, source,
    total = (annual - reference) / sd_lag_30,
    climate = (mean_lag_30 - reference) / sd_lag_30,
    gap = (hamilton_bk - mean_lag_30) / sd_lag_30,
    shock = (annual - hamilton_bk) / sd_lag_30
  ) %>%
  filter(!is.na(total), !is.na(gap)) %>%
  mutate(GID_0 = factor(GID_0, levels = COUNTRIES),
         source = factor(source, levels = CLIMATE_SOURCES))
PART_LABELS <- c(
  climate = paste0("Climate: trailing 30-year mean minus the ", REFERENCE_PERIOD[[1]], "-",
                   REFERENCE_PERIOD[[2]], " mean"),
  gap = "Gap: Hamilton trend minus trailing mean",
  shock = "Hamilton shock: year minus Hamilton trend"
)
parts <- decomposition %>%
  tidyr::pivot_longer(c(climate, gap, shock), names_to = "part", values_to = "value") %>%
  mutate(part = factor(PART_LABELS[part], levels = unname(PART_LABELS)))

# Same years on every panel, so the climate sources line up.
p <- ggplot() +
  geom_hline(yintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_col(data = decomposition,
           aes(x = year, y = total,
               fill = paste0("Total: year minus the ", REFERENCE_PERIOD[[1]], "-",
                             REFERENCE_PERIOD[[2]], " mean")),
           width = 0.8) +
  geom_line(data = parts, aes(x = year, y = value, colour = part), linewidth = 0.7) +
  facet_wrap(vars(source, GID_0), ncol = length(COUNTRIES),
             labeller = labeller(.multi_line = FALSE)) +
  scale_fill_manual(values = "grey80", name = NULL) +
  scale_colour_manual(values = stats::setNames(GRID_PALETTE[c(3, 1, 2)], unname(PART_LABELS)),
                      name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(
    x = NULL, y = "SDs of the trailing 30 years",
    title = paste0(VARIABLE_LABELS[[VARIABLE]],
                   " relative to the reference climate, split into climate, gap and Hamilton shock"),
    subtitle = paste0("Hamilton trend in real time from years t-2 and t-3 (Bilal & Känzig), ",
                      "fitted on the ", HAMILTON_WINDOW, " years up to t-2. Total = climate + ",
                      "gap + shock; all divided\nby the SD of years t-30 to t-1."),
    caption = paste0("Country series weighted by ", CLIMATE_WEIGHT, ". Reference: each ",
                     "country's mean over ", REFERENCE_PERIOD[[1]], "-", REFERENCE_PERIOD[[2]],
                     " in that climate source.")
  ) +
  theme_results()

file <- paste0("hamilton_decomposition_", VARIABLE, "_", paste(COUNTRIES, collapse = "_"), ".png")
ggsave(file.path(OUTPUT_DIR, file), p, width = 4.5 * length(COUNTRIES) + 1,
       height = 3.2 * length(CLIMATE_SOURCES) + 2, dpi = 300)
message("Written ", normalizePath(file.path(OUTPUT_DIR, file)))
