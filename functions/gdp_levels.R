# Total GDP and population by source, in a common GADM layout, for comparing
# growth across datasets (see national_subnational.R). Levels are in each
# source's own units and prices; only growth rates are comparable across
# sources.
#
# Columns: source, gadm_level, GID_0, GID_1 (= GID_0 for national sources),
# region_name, country, year, gdp, population, gdp_pc.

read_gdp_levels <- function(source) {
  stored <- resolve_econ_source(source)
  raw_file <- data_file(ECON_RAW_FILES[[stored]])

  out <- switch(
    stored,
    DOSE_V2_14 = read_dose(raw_file) %>%
      transmute(
        GID_0 = as.character(GID_0),
        GID_1 = as.character(GID_1),
        region_name = as.character(region),
        country = as.character(country),
        year = as.integer(year),
        # Constant 2015 local currency; see econ_panel.R.
        gdp_pc = as.numeric(grp_pc_lcu_2015),
        population = as.numeric(pop),
        gdp = gdp_pc * population
      ),
    # Uninhabited regions (zero population in every year, e.g. Kerguelen)
    # carry no GDP; kept, they would mark their country's coverage incomplete.
    KUMMU2025_GRID = load_kummu2025_grid() %>%
      group_by(GID_1) %>%
      filter(any(population > 0)) %>%
      ungroup() %>%
      transmute(
        GID_0 = as.character(GID_0),
        GID_1 = as.character(GID_1),
        region_name = as.character(NAME_1),
        country = as.character(COUNTRY),
        year = as.integer(year),
        gdp_pc = grp_pc_usd,
        population = population,
        gdp = total_gdp
      ),
    PWT110 = readxl::read_excel(raw_file, sheet = "Data") %>%
      transmute(
        GID_0 = as.character(countrycode),
        GID_1 = GID_0,
        region_name = NA_character_,
        country = as.character(country),
        year = as.integer(year),
        # rgdpna: millions, constant national prices in 2021 PPP US$;
        # pop: millions.
        gdp = as.numeric(rgdpna),
        population = as.numeric(pop),
        gdp_pc = gdp / population
      ),
    WB = read_wdi_csv(raw_file) %>%
      transmute(
        GID_0,
        GID_1 = GID_0,
        region_name = NA_character_,
        country,
        year,
        gdp = gdp_pc * population,          # Constant 2015 US$, market rates.
        population,
        gdp_pc
      )
  )

  out %>%
    mutate(
      source = stored,
      gadm_level = unname(GADM_LEVEL_BY_SOURCE[[stored]]),
      .before = 1
    )
}
