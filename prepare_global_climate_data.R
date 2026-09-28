# Build the global temperature series that the global-shock regressions use.
#
# The series comes from the Weighted Climate Dataset (Gortan, Testa, Fagiolo &
# Lamperti 2024, doi:10.1038/s41597-024-03304-1) at its `gadm_world` resolution -
# a single world-aggregate unit produced by the same gridded weighting machinery
# as the GADM0/GADM1 files already in `data/`. Two weightings are downloaded:
#
#   pop  ERA5, population density, base year 2015. Identical source, weight and
#        base year to the GADM1 panel the regressions are estimated on, so the
#        global and local temperatures are strictly comparable.
#   area ERA5, unweighted. The area-weighted global land mean, i.e. the closest
#        WCD analogue of a physical global mean surface temperature.
#
# Monthly files are downloaded and aggregated to years here rather than asking
# the WCD client for `time_frequency = "yearly"`, because the client's yearly
# aggregator is `mean()` without `na.rm` whereas prepare_climate_data.R drops NAs
# first. Aggregating locally keeps the global and regional annual series on the
# same operator.

suppressPackageStartupMessages(library(dplyr))

required_packages <- c("arrow", "curl", "ggplot2", "jsonlite", "tibble", "tidyr", "zoo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (file.exists("Econometrics.Rproj")) {
  project_root <- "."
} else if (file.exists(file.path("..", "Econometrics.Rproj"))) {
  project_root <- ".."
} else {
  stop("Run this script from the project root or the econometrics directory.")
}

source(file.path(project_root, "download_weighted_climate_data.R"))
source(file.path(project_root, "climate_transforms.R"))

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "global_local_temperature")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Configuration ----------------------------------------------------------------

# Mirrors prepare_climate_data.R so the global moments match the regional ones.
GLOBAL_WINDOWS <- c(5L, 10L, 20L, 30L)
GLOBAL_HAMILTON_LAGS <- 4L
GLOBAL_PERIOD_ALL <- c(1990, 2019)
GLOBAL_BASELINE_WINDOW <- 20L

# The requested selections. `expected_file` is asserted against what the WCD
# client actually resolves: .wcd_resolve() silently substitutes the only stored
# weight base year for a selection, and ERA5 at gadm_world exists only at 2015.
GLOBAL_SELECTIONS <- tibble::tribble(
  ~weight_label, ~climate_source, ~wcd_source, ~wcd_weight,          ~wcd_weight_year, ~expected_file,
  "pop",         "era5",          "ERA5",      "population density", 2015,             "gadm_world_era_tmp_pop_2015_monthly.parquet",
  "area",        "era5",          "ERA5",      "unweighted",         NA,               "gadm_world_era_tmp_un__monthly.parquet"
)

# Download and aggregate --------------------------------------------------------

# Same operator as aggregate_climate_parquet() in prepare_climate_data.R: pivot
# long, drop missing values, take the calendar-year mean of the monthly values.
aggregate_world_parquet <- function(path) {
  arrow::read_parquet(path) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(-Date, names_to = "unit", values_to = "value") %>%
    filter(!is.na(unit), !is.na(value), unit != "?") %>%
    mutate(year = as.integer(substr(gsub("^X", "", as.character(Date)), 1, 4))) %>%
    group_by(unit, year) %>%
    summarise(
      G = mean(value, na.rm = TRUE),
      months = n(),
      .groups = "drop"
    )
}

download_world_series <- function(selection) {
  cached_path <- wcd_download(
    variable = "avg. temperature",
    source = selection$wcd_source,
    geo_resolution = "gadm_world",
    weight = selection$wcd_weight,
    weight_year = if (is.na(selection$wcd_weight_year)) {
      NULL
    } else {
      selection$wcd_weight_year
    },
    time_frequency = "monthly",
    verbose = TRUE
  )

  if (!identical(basename(cached_path), selection$expected_file)) {
    stop(
      "The WCD client resolved '", basename(cached_path),
      "' but this script asked for '", selection$expected_file,
      "'. Weight base years may have changed upstream; reconcile before using ",
      "the series."
    )
  }

  # Keep a copy in data/ so the panel is reproducible without network access,
  # matching the existing raw gadm0_*/gadm1_* files.
  local_copy <- file.path(data_dir, selection$expected_file)
  if (!file.exists(local_copy) ||
      file.size(local_copy) != file.size(cached_path)) {
    file.copy(cached_path, local_copy, overwrite = TRUE)
  }

  annual <- aggregate_world_parquet(local_copy)

  n_units <- dplyr::n_distinct(annual$unit)
  if (n_units != 1L) {
    stop(
      "Expected exactly one world unit in ", selection$expected_file,
      " but found ", n_units, "."
    )
  }

  annual %>%
    mutate(
      gadm_level = "gadm_world",
      climate_source = selection$climate_source,
      weight = selection$weight_label,
      weight_year = if (is.na(selection$wcd_weight_year)) {
        "un"
      } else {
        as.character(selection$wcd_weight_year)
      },
      source_file = selection$expected_file
    )
}

world_annual <- bind_rows(
  lapply(seq_len(nrow(GLOBAL_SELECTIONS)), function(i) {
    download_world_series(GLOBAL_SELECTIONS[i, ])
  })
)

# Incomplete calendar years would bias the annual mean through the seasonal
# cycle, so they are dropped rather than carried.
incomplete_years <- world_annual %>% filter(months != 12L)
if (nrow(incomplete_years)) {
  message(
    "Dropping ", nrow(incomplete_years),
    " world-year(s) with fewer than 12 monthly observations: ",
    paste(sort(unique(incomplete_years$year)), collapse = ", ")
  )
  world_annual <- world_annual %>% filter(months == 12L)
}

# Moments and shocks ------------------------------------------------------------

# Rolling moments are right-aligned, exactly as in prepare_climate_data.R, and
# `zG` uses the *preceding* year's moments so the current realisation never
# enters its own baseline - the same operator behind `zTM` in test_functions.Rmd.
global_windows <- bind_rows(lapply(GLOBAL_WINDOWS, function(w) {
  world_annual %>%
    group_by(climate_source, weight, weight_year) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      mean_G = zoo::rollmean(G, w, align = "right", fill = NA),
      sd_G = zoo::rollapply(G, FUN = sd, width = w, align = "right", fill = NA),
      trend_G = hamilton_trend(G, h = w, p = GLOBAL_HAMILTON_LAGS),
      zG = standardized_anomaly(G, mean_G, sd_G),
      Gtrend = trend_G,
      Gcyc = G - trend_G,
      dG = G - dplyr::lag(G),
      window = w
    ) %>%
    ungroup()
}))

global_fixed <- world_annual %>%
  filter(year >= GLOBAL_PERIOD_ALL[[1]], year <= GLOBAL_PERIOD_ALL[[2]]) %>%
  group_by(climate_source, weight, weight_year) %>%
  summarise(
    mean_G_all = mean(G, na.rm = TRUE),
    sd_G_all = sd(G, na.rm = TRUE),
    .groups = "drop"
  )

global_panel <- global_windows %>%
  left_join(global_fixed, by = c("climate_source", "weight", "weight_year")) %>%
  mutate(dev_G_all = (G - mean_G_all) / sd_G_all) %>%
  select(
    gadm_level, climate_source, weight, weight_year, window, year,
    G, mean_G, sd_G, trend_G, zG, Gtrend, Gcyc, dG,
    mean_G_all, sd_G_all, dev_G_all, source_file
  ) %>%
  arrange(climate_source, weight, window, year)

arrow::write_parquet(
  global_panel,
  file.path(data_dir, "data_gm_gadmworld_era5_pop-area_2015.parquet")
)

# Diagnostics -------------------------------------------------------------------

# `zG` is reported before anything is interpreted: under monotone warming the
# global series sits above its own trailing mean in nearly every recent year, so
# zG is strongly positive and largely a trend rather than a mean-zero shock. That
# fact drives the fixed-effect ladder in the estimation script.
baseline_panel <- global_panel %>% filter(window == GLOBAL_BASELINE_WINDOW)

global_series_summary <- baseline_panel %>%
  group_by(climate_source, weight) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n(),
    mean_G_overall = mean(G),
    warming_per_decade = 10 * unname(coef(lm(G ~ year))[["year"]]),
    recent_minus_baseline = mean(G[year >= 2014 & year <= 2023]) -
      mean(G[year >= 1951 & year <= 1980]),
    n_zG = sum(is.finite(zG)),
    mean_zG = mean(zG, na.rm = TRUE),
    sd_zG = sd(zG, na.rm = TRUE),
    cor_zG_year = cor(zG, year, use = "complete.obs"),
    cor_Gcyc_year = cor(Gcyc, year, use = "complete.obs"),
    .groups = "drop"
  )

write.csv(
  baseline_panel,
  file.path(output_dir, "global_series.csv"),
  row.names = FALSE
)
write.csv(
  global_series_summary,
  file.path(output_dir, "global_series_summary.csv"),
  row.names = FALSE
)

series_plot_data <- baseline_panel %>%
  mutate(weight_label = if_else(
    weight == "pop",
    "Population-weighted (ERA5, 2015)",
    "Area-weighted (ERA5, unweighted)"
  ))

global_series_plot <- ggplot2::ggplot(
  series_plot_data,
  ggplot2::aes(x = year, y = G, colour = weight_label)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_line(
    ggplot2::aes(y = Gtrend),
    linewidth = 0.5,
    linetype = "dashed",
    na.rm = TRUE
  ) +
  ggplot2::scale_colour_manual(values = c("#0072B2", "#E69F00")) +
  ggplot2::labs(
    title = "Global land temperature, Weighted Climate Dataset (gadm_world, ERA5)",
    subtitle = "Solid: annual mean. Dashed: Hamilton (2018) trend, h = 20, p = 4.",
    x = NULL,
    y = "Global mean temperature (°C)",
    colour = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  file.path(output_dir, "global_series.png"),
  global_series_plot,
  width = 9,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

shock_plot_data <- series_plot_data %>%
  filter(is.finite(zG)) %>%
  tidyr::pivot_longer(
    c(zG, Gcyc),
    names_to = "shock",
    values_to = "value"
  ) %>%
  mutate(shock = if_else(
    shock == "zG",
    "zG: (G - lagged 20y mean) / lagged 20y SD",
    "Gcyc: G - Hamilton trend (°C)"
  ))

global_shock_plot <- ggplot2::ggplot(
  shock_plot_data,
  ggplot2::aes(x = year, y = value, colour = weight_label)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::facet_wrap(~shock, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c("#0072B2", "#E69F00")) +
  ggplot2::labs(
    title = "Two definitions of the global temperature shock",
    subtitle = paste(
      "zG is the same operator as the local zTM, and is dominated by trend;",
      "Gcyc is the detrended component that survives unit-specific trends."
    ),
    x = NULL,
    y = NULL,
    colour = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  file.path(output_dir, "global_shock_definitions.png"),
  global_shock_plot,
  width = 10,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

cat("\nGlobal temperature series (window ", GLOBAL_BASELINE_WINDOW, "):\n", sep = "")
print(as.data.frame(global_series_summary), digits = 4)
cat(
  "\nWritten: ",
  normalizePath(file.path(data_dir, "data_gm_gadmworld_era5_pop-area_2015.parquet")),
  "\n         ", normalizePath(output_dir), "\n",
  sep = ""
)
