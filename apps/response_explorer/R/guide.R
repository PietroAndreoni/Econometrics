# The guide page of the response explorer: how to use the app and what its
# options do.

guide_section <- function(title, ...) {
  card(card_header(title), card_body(...))
}
guide_item <- function(term, ...) tags$li(tags$strong(term), " ", ...)

guide_ui <- function() {
  tagList(
    guide_section(
      "What the app does",
      tags$p("The explorer estimates the growth response of GDP per capita to annual deviations of temperature (T)",
             "and precipitation (P) from their baseline, with the choices below, and uses each model for a single",
             "country (response curves and GDP paths) and for the world (climate scenarios and termination shocks",
             "of geoengineering). Up to four settings, each a model and a country, are compared side by side, each",
             "with a fixed colour."),
      tags$p("Models are fitted on demand when a setting changes (a few seconds; each one is then kept for the",
             "session). The climate scenarios and termination shocks run only when their button is pressed.")
    ),
    guide_section(
      "Workflow",
      tags$ol(
        tags$li("In the sidebar, set the model of each setting, in three sections: Common (data and regression design),",
                "Temperature and Precipitation (the terms of each variable, each with its own options). Give it a name",
                "if you like; \"Copy the model of setting 1\" copies setting 1's model."),
        tags$li("On the Single country page, Countries 1-4: tick \"Show this country\", give it a name, choose the",
                "setting whose model it uses and its characteristics. The charts label each curve \"country - setting\"",
                "(\"Country 1 - Setting 2\", or the names you gave) in the country's colour. A setting is shown on the",
                "other pages (coefficients, scenarios, termination) when a shown country uses it."),
        tags$li("Single country: the response of growth to a year's anomaly and the GDP path of a lasting shift, with",
                "the shift's size and timing set on that page."),
        tags$li("Coefficients: the equation, the estimates and the tests of each setting's model."),
        tags$li("Climate scenarios and Geoengineering: choose the paths and press Run; the results stay until the next run."),
        tags$li("\"Save the page as HTML\" (sidebar) saves the shown settings, charts, tables and notes, with the",
                "scenarios and termination shocks of their last run."),
        tags$li("Snapshots (sidebar): \"Save\" keeps every model option of the chosen setting under a name",
                "(results/response_explorer/snapshots.rds); \"Load into the setting\" puts a saved snapshot or a",
                "literature preset (Burke-Hsiang-Miguel 2015, Dell-Jones-Olken 2012, Kalkuhl-Wenz 2020,",
                "Kotz-Levermann-Wenz 2024, Harding et al. 2025, the preferred and the xi exposure specifications) into it. The",
                "presets use the app's options; the note under the menu says where they only approximate the paper."),
        tags$li("Snapshots of the app (sidebar): \"Save\" keeps the whole app under a name: every setting and country,",
                "the options of every page, the fitted models of the shown settings and the projections computed so",
                "far (results/response_explorer/app_snapshots/, one file each). \"Load\" sets it all back, then runs",
                "the climate scenarios and termination shocks again if they had been run, from the saved projections",
                "(a few seconds). Projections are computed only for pages that have been open, so open the pages",
                "you want to keep before saving.")
      )
    ),
    guide_section(
      "The model of a setting (sidebar)",
      tags$p("Common: GDP and climate data, outcome, weights, regression unit, style, lags, fixed effects and trends,",
             "convergence control, controls, outliers, the characteristics of the interactions (fixed or trailing), the",
             "low-income threshold, and the SD and characteristics in the paths. Temperature and Precipitation, each on",
             "its own: deviations in SD or absolute units, arms, symmetric, interactions, shape, cap, income scaling,",
             "quadratic climate control and exposure. The regression unit and style stay common: they define the rows",
             "and the demeaning of both variables."),
      tags$ul(
        guide_item("GDP data, climate data, outcome:", "PWT 11.0 or World Bank GDP; ERA5, CRU TS or UDelaware climate;",
                   "growth of GDP per capita (default), log GDP per capita (a level), or with PWT the growth of capital",
                   "per capita and of GDP per worker."),
        guide_item("Climate weights:", "population (default) or area averages over the country, in estimation and in",
                   "the CMIP6 downscaling of the scenarios."),
        guide_item("Deviations in:", "SD of the country's variability, or absolute units (C, m a year). Exposure",
                   "(below) needs absolute units and sets them."),
        guide_item("Regression unit:", "the level of T and P, the deviation from the horizon mean or from the trailing",
                   "5, 10, 20 or 30-year mean (years t-k..t-1; default 20), or the innovation of the Nath-Ramey-Klenow AR(3)",
                   "filter. Lags add every deviation term at t-1..t-L."),
        guide_item("Regression style:", "the unit itself, its first difference (with the level: Kalkuhl and Wenz), or both."),
        guide_item("Fixed effects and trends:", "year and unit effects (default), unit only or year only; unit linear and",
                   "quadratic trends (default), linear only or none. The deviations are demeaned in the same way."),
        guide_item("Convergence and quadratic controls:", "an optional control for past GDP, and optional Z + Z² controls",
                   "of the climate; both on their own, not in the responses or projections."),
        guide_item("Controls and outliers:", "conflict and the post-Soviet and post-Yugoslav collapse years; growth",
                   "outliers whose cause is not weather dropped (default: the 5% tails of growth reversals)."),
        guide_item("Arms:", "the positive and negative deviations of T and P enter separately (warm and cold, wet and",
                   "dry), linear, quadratic or both. Symmetric: one coefficient for both sides."),
        guide_item("Interactions:", "each arm times the long-run climate, last year's or this year's value, income,",
                   "agricultural share, capital (PWT) or a low-income dummy (threshold in the GDP source's units).",
                   "Linear or restricted cubic splines."),
        guide_item("Income scaling:", "every T (or P) term times (Ybar / Y)^(1 - xi): the dollar loss per person is",
                   "proportional to Y^xi (xi = 0: equal dollars; xi = 1: equal %, the standard terms). \"Both\" keeps",
                   "the plain and the scaled terms (a test of the two)."),
        guide_item("Exposure of the T and P terms:", "every term of the variable times the country's exposure slope",
                   "of two indices, with no plain term, so that a deviation acts through the extreme days it brings.",
                   "T: hot (e.g. degree days of Tmax above 29 C) and cold (e.g. of Tmin below 0 C); P: wet (mm beyond",
                   "the country's 99th percentile of daily P) and dry (days below 1 mm). The slope is the change of the",
                   "index per C (m) of the annual value, from ERA5 daily data; \"Country\" uses each country's own slope,",
                   "\"Smooth\" the slope as a function of the long-run climate, which follows the climate in the paths",
                   "and scenarios. Thresholds, daily series (Tmax, mean, Tmin) and measure (degree days or days) can be",
                   "changed; the exposure explorer (apps/exposure_explorer) shows the indices themselves."),
        guide_item("Characteristics in the interactions:", "fixed (sample means) or trailing (the climate of the",
                   "baseline, the income of the 5 years before, with their level terms; with WB GDP optionally the",
                   "share too)."),
        guide_item("Cap:", "income, share and capital interactions can bring the non-climate response to zero but",
                   "never create a benefit."),
        guide_item("SD in the GDP path and long-run characteristics:", "with a trailing baseline, the SD of the",
                   "deviation fixed at the country's or widened by the shift; the characteristics of the interactions",
                   "fixed, or following the trailing baseline (climate; in the scenarios also income and share along",
                   "the SSP).")
      ),
      tags$p("The countries are on the Single country page: start from a preset and edit its long-run climate,",
             "variability, income, share and capital. The warnings above every page flag, for each shown country and",
             "its setting, characteristics outside the sample's countries and shifts beyond the deviations in the data.")
    ),
    guide_section(
      "Single country",
      tags$ul(
        guide_item("Countries:", "up to four, each shown or not, with a name, the setting whose model it uses (several",
                   "countries can use the same setting, or one country can be repeated with different settings), a",
                   "preset from that setting's panel and its characteristics."),
        guide_item("Shift:", "the size of the temperature and precipitation shifts (common to all settings), the years",
                   "shown, the year X and which shifts to draw: immediate, gradual until year X, sudden in year X, or",
                   "transitory (year 0 only)."),
        guide_item("Response curves:", "the effect on growth of a year at each anomaly, with its 95% CI; solid within",
                   "the 1st-99th percentile of the deviations in the data. Hover for the values."),
        guide_item("GDP paths:", "the effect of the shift on the level of GDP (the sum of the growth effects), % of",
                   "the GDP without it, with its 95% CI.")
      )
    ),
    guide_section(
      "By characteristic",
      tags$p("The effect on growth of one year k SD above (red) and below (blue) the baseline, in one chart for",
             "temperature and one for precipitation (columns: the shown settings), each against a characteristic",
             "of its own varying over the range of the setting's countries (rug): long-run temperature or",
             "precipitation, GDP per capita, agricultural share or capital. The other characteristics, and the SD that sets the size of the anomaly, are those",
             "of the reference country (one of the shown countries; circle). The heatmap shows the same effect over",
             "two characteristics at once (axes chosen freely), for T or P and the +k or -k SD year, with the",
             "setting's countries as dots (where the data are). With exposure-mediated terms, the",
             "second chart shows the physical exposure curves: each country's slope of the indices and their smooth",
             "in climate, with the reference country's climate dashed.")
    ),
    guide_section(
      "Coefficients",
      tags$p("Each setting's equation and definitions, its deviation terms with clustered SEs and stars (* 10%,",
             "** 5%, *** 1%), the fit (within adjusted R², AIC, BIC) and joint Wald tests (all terms of a variable,",
             "each interaction, each exposure channel, ...). Settings 2-4 get verdicts against setting 1 by AIC / BIC",
             "/ within R² when they run on the same rows: \"indifferent\" when the difference cannot be told from noise",
             "(clustered Vuong test for non-nested models, Wald test of the added terms for nested ones). For",
             "nested pairs read the Wald test.")
    ),
    guide_section(
      "Climate scenarios",
      tags$ul(
        guide_item("Forcing and economy:", "the global warming paths (median of the forcing's runs) and the SSP of",
                   "income, population and share (each forcing's own, or one SSP for all)."),
        guide_item("Natural variability:", "none (forced signal only), or each country's detrended historical T and P",
                   "added as noise: years independent across countries, joint years (spatial correlation) or joint",
                   "5-year blocks (also persistence), over the chosen number of noise paths (slower)."),
        guide_item("Ignore gains (also on the Geoengineering page):", "in the world aggregate a country whose GDP",
                   "gains counts at its baseline GDP (its effect capped at zero each year, at every draw of the",
                   "coefficients and noise path); it keeps its GDP weight, and the maps and country results are",
                   "unchanged."),
        guide_item("Outputs:", "the warming paths, the world loss (or the world growth rate and its change), maps of",
                   "the loss by country in 2050 and 2100 (setting, forcing, T and / or P), the physical drivers (each",
                   "country's change of T and P) and a table."),
        guide_item("Countries (also on the Geoengineering page):", "the global chart, its hover readout and the",
                   "historical growth show one country, or the GDP-weighted aggregate of the selected countries (their",
                   "SSP GDP, with its CI and noise bands), in place of the world; no new run needed. Maps and tables stay",
                   "global.")
      )
    ),
    guide_section(
      "Geoengineering",
      tags$p("Geoengineering limits global warming over the background forcing from 2025: to a level (any number of",
             "C above 1850-1900, 0 being preindustrial) or to a rate (at most so many C per decade, 0 stopping",
             "warming). In each ticked termination year it stops and warming catches up with the background with",
             "time constant tau; with none ticked, only the background and sustained geoengineering run. The page",
             "shows the warming paths with a band of exceedances of natural variability (k SD above the reference",
             "climate), the world loss, a map of the loss by country, the exposure by latitude and year (events",
             "beyond k SD), and a table."),
      tags$p("Regional pattern: uniform (each country's CMIP6 pattern of global warming: perfect geoengineering) or",
             "CIDER's, stratospheric aerosol injection at the latitudes and fractions set (CESM2-WACCM6 or UKESM1),",
             "which leaves a regional climate even at constant global temperature; its map shows that residual."),
      tags$p("\"Compare side by side\" runs several background forcings, economic developments (SSPs), catch-up",
             "times or geoengineering patterns at once, each in place of that option: the warming and world charts get one column per",
             "alternative (the settings in rows), the table a column naming it, and the map, latitude chart and",
             "exceedances show the alternative chosen above them. With the growth measure, both world charts can show",
             "the historical growth of world GDP in PWT 11.0 (1971-2023: the series, and its median, interquartile",
             "range and range on the right).")
    ),
    guide_section(
      "Caveats",
      tags$ul(
        tags$li("Clustered SEs treat the countries as independent; residual spatial correlation makes them too small."),
        tags$li("Models chosen after looking at fits overstate significance; prefer carrying several indifferent",
                "models through the projections."),
        tags$li("Projections extrapolate the estimated responses to climates and incomes beyond the sample: see the",
                "warnings and the support of the response curves."),
        tags$li("Exposure slopes come from ERA5 with 2015 population weights at country level (also with CRU or",
                "UDelaware panels); national daily values smooth local extremes, rain more than temperature.")
      )
    )
  )
}
