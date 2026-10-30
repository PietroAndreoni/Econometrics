# Weather response explorer: the growth response to temperature and
# precipitation deviations, and the GDP effect of a permanent shift, for up to
# three settings (a country's characteristics and a model each), under models
# fitted on demand (R/fit.R) on the panels of build_models.R.
#
#   Rscript apps/response_explorer/build_models.R   # once, after rendering the notebook
#   shiny::runApp("apps/response_explorer")         # or source() this file
#
# R/fit.R fits the models, R/response.R and R/projection.R hold the computations.

library(shiny)
library(bslib)
library(ggplot2)

# The app's folder: the file's own when it is source()d, otherwise the working
# directory (runApp() sets it to the app's folder) or apps/response_explorer
# under it (the project root).
SOURCED_FILE <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
APP_DIR <- local({
  candidates <- c(if (!is.null(SOURCED_FILE)) dirname(SOURCED_FILE), getwd(),
                  file.path(getwd(), "apps", "response_explorer"))
  found <- candidates[file.exists(file.path(candidates, "R", "response.R"))]
  if (!length(found)) stop("Cannot find the app's folder (apps/response_explorer).")
  normalizePath(found[[1]])
})
source(file.path(APP_DIR, "R", "response.R"), local = environment())
source(file.path(APP_DIR, "R", "projection.R"), local = environment())
source(file.path(APP_DIR, "R", "fit.R"), local = environment())

PANELS_FILE <- normalizePath(file.path(APP_DIR, "..", "..", "results", "response_explorer", "panels.rds"),
                             mustWork = FALSE)
if (!file.exists(PANELS_FILE)) {
  stop("No ", PANELS_FILE, ": run Rscript apps/response_explorer/build_models.R from the project root.")
}
EXPLORER <- readRDS(PANELS_FILE)
SCENARIOS_FILE <- sub("panels[.]rds$", "scenarios.rds", PANELS_FILE)
if (!file.exists(SCENARIOS_FILE)) {
  stop("No ", SCENARIOS_FILE, ": run Rscript apps/response_explorer/build_scenarios.R from the project root.")
}
SCEN <- readRDS(SCENARIOS_FILE)
# The maps draw the sf borders of SCEN$world.
loadNamespace("sf")
# Each forcing scenario (code ssp119, ...) gives the global warming; its
# forcing in W/m2 from the code (ssp245: 4.5).
FORCING <- stats::setNames(sub("^ssp[0-9]([0-9])([0-9])$", "\\1.\\2", names(SCEN$scenarios)),
                           names(SCEN$scenarios))
FORCING_LABELS <- stats::setNames(paste0(FORCING, " W/m2"), names(FORCING))
SCENARIO_CHOICES <- stats::setNames(names(SCEN$scenarios),
                                    sprintf("%s (%s)", FORCING_LABELS, SCEN$scenarios))
# The economic development (income, population, agricultural share): each
# forcing's own SSP, or one SSP for all.
ECONOMY_CHOICES <- c("Follow SSP-RCP" = "follow", SSP1 = "SSP1", SSP2 = "SSP2", SSP3 = "SSP3",
                     SSP4 = "SSP4", SSP5 = "SSP5")
economy_of <- function(scenario, economy) {
  if (economy == "follow") SCEN$ssp$of_scenario[[scenario]] else economy
}
# The label of a forcing with an economy: "SSP2-4.5" when they are paired as in
# the scenario, otherwise "4.5 W/m2 x SSP3".
scenario_label <- function(scenario, economy) {
  ssp <- economy_of(scenario, economy)
  if (identical(ssp, SCEN$ssp$of_scenario[[scenario]])) SCEN$scenarios[[scenario]]
  else sprintf("%s x %s", FORCING_LABELS[[scenario]], ssp)
}
# One blue step per forcing, light to dark (ordinal ramp).
SCENARIO_COLOURS <- stats::setNames(c("#86b6ef", "#5598e7", "#3987e5", "#2a78d6", "#1c5cab",
                                      "#104281", "#0d366b"),
                                    names(SCEN$scenarios))
# Diverging map colours: losses red, gains blue, grey at zero; countries without data.
MAP_LOSS <- "#e34948"
MAP_GAIN <- "#2a78d6"
MAP_ZERO <- "#f0efec"
MAP_NA <- "#d9d8d4"

N_SETTINGS <- 3
# One fixed colour per setting (categorical slots 1-3, validated for any pair).
SETTING_COLOURS <- c("#2a78d6", "#eb6834", "#1baf7a")
DEFAULT_COUNTRIES <- c("IND", "ITA", "USA")
INK <- "#0b0b0b"
INK_MUTED <- "#52514e"
GRID <- "#e6e5e1"
SURFACE <- "#fcfcfb"
FORM_SHORT <- c("Linear arms" = "linear arms", "Quadratic arms" = "quadratic arms")
# The modifiers a variable's deviation terms can be interacted with (in the
# order of build_models.R), their names in the inputs and in the labels.
MODIFIERS <- c("Long-run climate" = "clim", "Income" = "inc", "Agricultural share" = "agr",
               "Capital per capita (PWT)" = "kap")
# Ticked by default (capital is an extra).
DEFAULT_MODIFIERS <- c("clim", "inc", "agr")
MODIFIER_SHORT <- c(clim = "climate", inc = "income", agr = "share", kap = "capital")

# What the paths and losses of settings `sp` are of: GDP, capital per capita,
# GDP per worker, or "outcome" when the settings differ (each setting's label
# names its outcome).
quantity_of <- function(sp) {
  outcomes <- unique(vapply(sp, function(s) s$model$outcome %||% "log", ""))
  if (identical(outcomes, "log")) return("GDP")
  if (length(outcomes) > 1) return("outcome")
  c(capital = "capital per capita", worker = "GDP per worker")[[outcomes]]
}
first_upper <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))

# Short description of the interactions: one list if T and P have the same.
interaction_label <- function(inter) {
  part <- function(x) if (length(x)) paste(MODIFIER_SHORT[x], collapse = " + ") else "none"
  if (identical(inter$TM, inter$RR)) {
    if (length(inter$TM)) part(inter$TM) else "no interactions"
  } else {
    paste0("T: ", part(inter$TM), ", P: ", part(inter$RR))
  }
}

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))
# Readable names of the model's terms: "T warm arm x long-run T (spline 1)",
# "level: log income".
ARM_LABELS <- c(TM_pos = "T warm arm", TM_neg = "T cold arm", RR_pos = "P wet arm", RR_neg = "P dry arm")
CHARACTERISTIC_LABELS <- c(clim_TM = "long-run T", clim_RR = "long-run P", inc = "log income",
                           agr = "agricultural share", kap = "log capital per capita")
characteristic_label <- function(z) {
  base <- sub("_nl[0-9]+$", "", z)
  out <- CHARACTERISTIC_LABELS[[base]]
  if (base != z) out <- paste0(out, " (spline ", sub(".*_nl", "", z), ")")
  out
}
term_label <- function(term) {
  if (startsWith(term, "level:idol_")) {
    return(paste0("Level: (Ybar / Y)^(1 - xi_", if (endsWith(term, "TM")) "T" else "P", ")"))
  }
  if (startsWith(term, "level:")) return(paste("Level:", characteristic_label(sub("^level:", "", term))))
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  scaled <- "idol" %in% parts
  parts <- parts[parts != "idol"]
  arm <- substring(parts[[1]], 4)
  label <- ARM_LABELS[[sub("2$", "", arm)]]
  if (endsWith(arm, "2")) label <- paste(label, "(squared)")
  if (length(parts) > 1) label <- paste(label, "×", characteristic_label(parts[[2]]))
  if (scaled) label <- paste0(label, " × (Ybar / Y)^(1 - xi_", if (startsWith(arm, "TM")) "T" else "P", ")")
  label
}

# The specification of a model in LaTeX (MathJax): the equation, and the
# definitions of its terms (HTML with inline maths).
spec_symbols <- function(m) {
  trailing <- identical(m$moderators, "trailing")
  clim_t <- trailing && m$window != "all"
  list(
    clim = function(v) paste0("c_{", v, if (clim_t) ",it" else ",i", "}"),
    inc = if (trailing) "y_{it}" else "y_{i}",
    agr = if (isTRUE(m$share_trailing)) "s_{it}" else "s_{i}",
    kap = if (trailing) "k_{it}" else "k_{i}",
    clim_t = clim_t, inc_t = trailing, agr_t = isTRUE(m$share_trailing)
  )
}
spec_equation <- function(m) {
  sym <- spec_symbols(m)
  power <- if (m$form == "Quadratic arms") "^{2}" else ""
  splines <- identical(m$shape, "splines")
  z_of <- function(v, z) if (z == "clim") sym$clim(v) else sym[[z]]
  letter <- c(clim = "c", inc = "y", agr = "s", kap = "k")
  # The arm's coefficient and interactions; with a tilde, those of the terms
  # times Ybar / Y (equal dollars).
  inner_of <- function(v, tilde) {
    t <- function(s) if (tilde) paste0("\\tilde{", s, "}") else s
    paste0(t("\\beta"), "^{a}_{", v, "}", paste0(vapply(m$inter[[if (v == "T") "TM" else "RR"]], function(z) {
      if (splines) paste0(" + ", t("g"), "^{a}_{", v, ",", letter[[z]], "}\\!\\left(", z_of(v, z), "\\right)")
      else paste0(" + ", t("\\gamma"), "^{a}_{", v, ",", letter[[z]], "}\\, ", z_of(v, z))
    }, ""), collapse = ""))
  }
  forms <- model_income_forms(m)
  ratio_of <- function(column) {
    power <- 1 - forms$xi[[column]]
    base <- paste0("\\frac{\\bar Y}{Y_{", if (identical(m$moderators, "trailing")) "it" else "i", "}}")
    if (power == 1) base else paste0("\\left(", base, "\\right)^{", format(power), "}")
  }
  block <- function(v) {
    column <- if (v == "T") "TM" else "RR"
    ratio <- ratio_of(column)
    inner <- switch(forms$form[[column]],
                    standard = paste0("\\Big(", inner_of(v, FALSE), "\\Big)"),
                    dollars = paste0("\\Big(", inner_of(v, TRUE), "\\Big) ", ratio),
                    both = paste0("\\Big[\\Big(", inner_of(v, FALSE), "\\Big) + \\Big(", inner_of(v, TRUE),
                                  "\\Big) ", ratio, "\\Big]"))
    paste0("\\sum_{a \\in \\{+,-\\}} ", inner, " \\big(D^{a}_{", v, ",it}\\big)", power)
  }
  levels <- unique(sub("_nl[0-9]+$", "", m$levels$term))
  level_terms <- vapply(levels, function(z) {
    if (startsWith(z, "idol_")) {
      return(paste0(" + \\phi_{r,", if (z == "idol_TM") "T" else "P", "}\\, ", ratio_of(sub("^idol_", "", z))))
    }
    x <- switch(z, clim_TM = sym$clim("T"), clim_RR = sym$clim("P"), inc = sym$inc, agr = sym$agr, kap = sym$kap)
    name <- switch(z, clim_TM = "c_T", clim_RR = "c_P", inc = "y", agr = "s", kap = "k")
    if (splines) paste0(" + h_{", name, "}\\!\\left(", x, "\\right)") else paste0(" + \\phi_{", name, "}\\, ", x)
  }, "")
  fe <- switch(m$fe, unit = "\\lambda_t + \\alpha_i", linear = "\\lambda_t + \\alpha_i + \\delta_i t",
               quadratic = "\\lambda_t + \\alpha_i + \\delta_i t + \\eta_i t^{2}")
  lhs <- switch(m$outcome %||% "log",
                log = "\\Delta \\ln \\mathrm{GDPpc}_{it}",
                capital = "\\Delta \\ln \\mathrm{Kpc}_{it}",
                worker = "\\Delta \\ln \\big(\\mathrm{GDP}_{it} / L_{it}\\big)")
  paste0("\\begin{aligned} ", lhs, " = {} & ", block("T"), " \\\\ & + ", block("P"),
         " \\\\ & ", paste(level_terms, collapse = ""), " + \\kappa' C_{it} + ", fe, " + \\varepsilon_{it}",
         " \\end{aligned}")
}
spec_definitions <- function(m) {
  sym <- spec_symbols(m)
  used <- unique(unlist(m$inter))
  start <- m$horizon[[1]]
  end <- m$horizon[[2]]
  w <- if (m$window == "all") NA else sub("lag_", "", m$window)
  absolute <- identical(m$units, "absolute")
  baseline <- if (is.na(w)) {
    sprintf(paste0("\\(\\mu_{v,i}\\) and \\(\\sigma_{v,i}\\): the mean and SD of \\(v_{it}\\) over %d-%d ",
                   "(horizon baseline)"), start, end)
  } else {
    sprintf(paste0("\\(\\mu_{v,it} = \\frac{1}{%s}\\sum_{k=1}^{%s} v_{i,t-k}\\) and \\(\\sigma_{v,it}\\) the SD ",
                   "over the same years (trailing %s-year baseline, without year \\(t\\))"), w, w, w)
  }
  mu <- if (is.na(w)) "\\mu_{v,i}" else "\\mu_{v,it}"
  sigma <- if (is.na(w)) "\\sigma_{v,i}" else "\\sigma_{v,it}"
  deviation <- if (absolute) {
    sprintf("Deviation of \\(v \\in \\{T, P\\}\\) in absolute units (C, m): \\(x_{v,it} = v_{it} - %s\\), with %s.",
            mu, baseline)
  } else {
    sprintf("Deviation of \\(v \\in \\{T, P\\}\\) in SD: \\(x_{v,it} = (v_{it} - %s) / %s\\), with %s.",
            mu, sigma, baseline)
  }
  trend <- switch(m$fe, unit = "\\hat a_i", linear = "\\hat a_i + \\hat b_i t",
                  quadratic = "\\hat a_i + \\hat b_i t + \\hat c_i t^{2}")
  outcome <- switch(m$outcome %||% "log",
    log = NULL,
    capital = paste0(
      "Outcome: the growth of the capital stock per capita, \\(\\mathrm{Kpc}_{it}\\) = PWT rnna / pop (constant ",
      "national prices, as GDP). Responses and losses are of capital per capita."),
    worker = paste0(
      "Outcome: the growth of GDP per worker, \\(\\mathrm{GDP}_{it} / L_{it}\\) = PWT rgdpna / emp (persons ",
      "engaged): the growth of GDP per capita less that of employment per head. Responses and losses are of ",
      "GDP per worker."))
  out <- c(
    outcome,
    deviation,
    sprintf(paste0("Demeaned before the arms are cut: \\(d_{v,it} = x_{v,it} - (%s)\\), fitted by OLS for each ",
                   "country over %d-%d; arms \\(D^{+}_{v,it} = \\max(d_{v,it}, 0)\\), ",
                   "\\(D^{-}_{v,it} = \\max(-d_{v,it}, 0)\\)%s."),
            trend, start - 10L, end, if (m$form == "Quadratic arms") ", entering squared" else "")
  )
  if ("clim" %in% used) {
    out <- c(out, if (sym$clim_t) {
      sprintf(paste0("Climate: \\(c_{v,it} = \\mu_{v,it} - \\operatorname{med}_j \\bar v_j\\), the trailing ",
                     "%s-year mean (the deviation's own baseline) less the median over the countries of their ",
                     "%d-%d means."), w, start, end)
    } else {
      sprintf(paste0("Climate: \\(c_{v,i} = \\bar v_i - \\operatorname{med}_j \\bar v_j\\), the country's mean ",
                     "over %d-%d less the median over the countries (constant: absorbed by \\(\\alpha_i\\) on ",
                     "its own)."), start, end)
    })
  }
  if ("inc" %in% used) {
    out <- c(out, if (sym$inc_t) {
      paste0("Income: \\(y_{it} = \\ln\\big(\\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\big)\\) less ",
             "its mean over the rows (without year \\(t\\)).")
    } else {
      sprintf(paste0("Income: \\(y_{i} = \\ln\\big(\\overline{\\mathrm{GDPpc}}_i\\big)\\), the log of mean GDP ",
                     "per capita over %d-%d, less its mean over the rows."), start, end)
    })
  }
  if ("agr" %in% used) {
    out <- c(out, if (sym$agr_t) {
      paste0("Agricultural share: \\(s_{it} = \\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{share}_{i,t-k}\\) less its ",
             "mean over the rows (without year \\(t\\)).")
    } else {
      sprintf(paste0("Agricultural share: \\(s_{i}\\), the mean of the previous year's share over %d-%d ",
                     "(that is, %d-%d), less its mean over the rows."), start, end, start - 1L, end - 1L)
    })
  }
  if ("kap" %in% used) {
    out <- c(out, if (sym$inc_t) {
      paste0("Capital: \\(k_{it} = \\ln\\big(\\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{Kpc}_{i,t-k}\\big)\\) less ",
             "its mean over the rows (without year \\(t\\)), with \\(\\mathrm{Kpc}\\) = PWT rnna / pop.")
    } else {
      sprintf(paste0("Capital: \\(k_{i} = \\ln\\big(\\overline{\\mathrm{Kpc}}_i\\big)\\), the log of the mean ",
                     "capital stock per capita (PWT rnna / pop) over %d-%d, less its mean over the rows."),
              start, end)
    })
  }
  forms <- model_income_forms(m)
  scaled_v <- c(TM = "T", RR = "P")[forms$form != "standard"]
  if (length(scaled_v)) {
    y <- if (sym$inc_t) "Y_{it}" else "Y_{i}"
    xi_text <- paste(vapply(names(scaled_v), function(k) {
      sprintf("\\(\\xi_{%s} = %s\\)", scaled_v[[k]], format(forms$xi[[k]]))
    }, ""), collapse = ", ")
    # The country's income of the ratio is that of the interactions (long-run
    # or trailing); Ybar is a constant, the geometric mean of it over the rows.
    income_text <- if (sym$inc_t) {
      paste0("\\(Y_{it} = \\frac{1}{5}\\sum_{k=1}^{5} \\mathrm{GDPpc}_{i,t-k}\\), the country's trailing income")
    } else {
      sprintf("\\(Y_{i} = \\overline{\\mathrm{GDPpc}}_i\\), the country's mean GDP per capita over %d-%d", start, end)
    }
    out <- c(out, paste0(
      "Scaled terms of ", paste(scaled_v, collapse = " and "), ": times \\((\\bar Y / ", y, ")^{1-\\xi}\\), with ",
      income_text, " (the income of the interactions), and \\(\\bar Y = \\exp\\big(\\overline{\\ln ", y,
      "}\\big)\\) a constant, the geometric mean of ", if (sym$inc_t) "the trailing income" else "it",
      " over the rows (so \\((\\bar Y / ", y, ")^{1-\\xi} = e^{-(1-\\xi)", sym$inc, "}\\), with \\(", sym$inc,
      "\\) the log income less its mean over the rows); ", xi_text, ". A scaled growth effect is a dollar loss per ",
      "person proportional to \\(", y, "^{\\xi}\\) beyond the ticked interactions (\\(\\xi = 0\\): the same ",
      "dollars at every income, \\(\\tilde\\beta \\bar Y\\) at the median climate and mean share; \\(\\xi = 1\\): ",
      "the same %, the standard terms)",
      if (any(forms$form == "both")) paste0(
        "; with the plain terms (", paste(c(TM = "T", RR = "P")[forms$form == "both"], collapse = " and "),
        "), the model nests equal % (\\(\\tilde\\beta = \\tilde\\gamma = 0\\)) and the scaled form ",
        "(\\(\\beta = \\gamma = 0\\)), tested in the table"), "."))
  }
  if (identical(m$shape, "splines")) {
    out <- c(out, paste0("Splines: \\(g(z) = \\gamma z + \\sum_{j=1}^{2} \\gamma_j [r_j(z) - r_j(0)]\\), ",
                         "restricted cubic spline with knots at the 5, 35, 65 and 95th percentiles of the ",
                         "countries' values (linear beyond the outer knots); one set of coefficients per arm."))
  }
  if (nrow(m$levels)) {
    out <- c(out, paste0("Level terms (\\(\\phi\\)", if (identical(m$shape, "splines")) ", \\(h\\)", "): the ",
                         "time-varying characteristics that enter an interaction also enter on their own."))
  }
  c(out, paste0(
    "\\(\\lambda_t\\): year effects; \\(\\alpha_i\\): country effects",
    switch(m$fe, unit = "", linear = "; \\(\\delta_i t\\): country linear trends",
           quadratic = "; \\(\\delta_i t + \\eta_i t^{2}\\): country linear and quadratic trends"),
    if (m$fe != "unit") " (\\(t\\) = year - 1990)", "; \\(C_{it}\\): conflict (minor, war) and Soviet and ",
    "Yugoslav collapse-year dummies; \\(\\varepsilon_{it}\\) clustered by country."))
}

# Legend labels wrapped so that long setting labels stay inside the chart.
wrap_label <- function(x, width = 75) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "")
}

theme_explorer <- function() {
  theme_minimal(base_size = 13) +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA),
          panel.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = INK_MUTED),
          axis.title = element_text(colour = INK_MUTED, size = 11),
          plot.subtitle = element_text(colour = INK_MUTED, size = 10.5),
          legend.position = "top", legend.justification = "left", legend.direction = "vertical",
          legend.text = element_text(colour = INK_MUTED, size = 10.5))
}

# Inputs of setting n: whether it is shown (settings 2 and 3), its label, its
# model and its country.
setting_ui <- function(id, n) {
  ns <- NS(id)
  tagList(
    if (n > 1) checkboxInput(ns("active"), "Show this setting", value = FALSE),
    textInput(ns("label"), "Label", placeholder = "automatic"),
    tags$h6("Model", class = "mt-3"),
    if (n > 1) actionButton(ns("copy"), "Copy the model of setting 1",
                            class = "btn-sm btn-outline-secondary mb-2"),
    selectInput(ns("gdp"), "GDP data", c("PWT", "WB")),
    selectInput(ns("climate"), "Climate data", c("ERA5", "CRU TS", "UDelaware")),
    selectInput(ns("outcome"), "Outcome (capital and GDP per worker: PWT only)", OUTCOME_CHOICES, selected = "log"),
    radioButtons(ns("weight"), "Climate weights (estimation and CMIP6 downscaling)",
                 c("Population density" = "population", "Area (plain area average)" = "area"),
                 selected = "population", inline = TRUE),
    radioButtons(ns("units"), "Deviations in",
                 c("SD of the country's variability" = "SD", "Absolute units (C, m)" = "absolute"),
                 selected = "SD", inline = TRUE),
    selectInput(ns("baseline"), "Long-run baseline of the deviations",
                c("Horizon (sample period)" = "Horizon", "Trailing 20 years" = "20-year",
                  "Trailing 30 years" = "30-year")),
    selectInput(ns("fe"), "Fixed effects (the deviations are demeaned in the same way)",
                FE_CHOICES, selected = "quadratic"),
    selectInput(ns("form"), "Arms", c("Linear arms", "Quadratic arms")),
    checkboxGroupInput(ns("inter_T"), "Interactions of the temperature terms", MODIFIERS,
                       selected = DEFAULT_MODIFIERS, inline = TRUE),
    checkboxGroupInput(ns("inter_P"), "Interactions of the precipitation terms", MODIFIERS,
                       selected = DEFAULT_MODIFIERS, inline = TRUE),
    tagList(
      helpText("Income scaling of the deviation terms (each arm and its ticked interactions), by variable:",
               "the dollar loss per person of a scaled term is proportional to Y^xi (xi = 0: equal dollars;",
               "xi = 1: equal %, the standard terms)."),
      lapply(c(T = "temperature", P = "precipitation"), function(name) {
        v <- if (name == "temperature") "T" else "P"
        tagList(
          radioButtons(ns(paste0("income_form_", v)), paste("Income scaling of the", name, "terms"),
                       INCOME_FORMS, selected = "standard"),
          conditionalPanel(
            sprintf("input.income_form_%s != 'standard'", v), ns = ns,
            sliderInput(ns(paste0("xi_", v)), paste0("Elasticity xi of the scaled ", name, " terms"),
                        min = XI_RANGE[1], max = XI_RANGE[2], value = 0, step = 0.1)
          )
        )
      })
    ),
    radioButtons(ns("moderators"), "Characteristics in the interactions (estimation)",
                 c("Fixed: sample means" = "fixed",
                   "Trailing: climate of the deviation's baseline, income of the 5 years before, with level terms" =
                     "trailing"), selected = "fixed"),
    conditionalPanel(
      "input.moderators == 'trailing' && input.gdp == 'WB'", ns = ns,
      checkboxInput(ns("share_trailing"), "Trailing agricultural share too (5 years before)", value = FALSE),
      uiOutput(ns("share_warning"))
    ),
    radioButtons(ns("shape"), "Shape of the interactions",
                 c("Linear" = "linear",
                   "Restricted cubic splines (knots at the 5, 35, 65, 95th percentiles of the countries)" =
                     "splines"), selected = "linear"),
    checkboxInput(ns("cap"), paste("Cap the income, share and capital interactions: they can bring the",
                                   "non-climate response to zero but never create a benefit"),
                  value = FALSE),
    conditionalPanel(
      "input.baseline != 'Horizon' && input.units == 'SD'", ns = ns,
      radioButtons(ns("sd_mode"), "SD of the trailing baseline in the GDP path",
                   c("Fixed at the country's SD" = "fixed",
                     "Moving (widened by the shift)" = "moving"), selected = "fixed")
    ),
    conditionalPanel(
      "input.baseline != 'Horizon' || input.moderators == 'trailing'", ns = ns,
      radioButtons(ns("clim_mode"),
                   "Long-run characteristics of the interactions (GDP paths and climate scenarios)",
                   c("Fixed at the sample mean" = "fixed",
                     "Follow the trailing baseline (climate; in the scenarios also income and share along the SSP)" =
                       "baseline"), selected = "fixed")
    ),
    tags$h6("Country", class = "mt-3"),
    selectInput(ns("preset"), "Start from a country", choices = NULL),
    numericInput(ns("mean_T"), "Long-run mean temperature (°C)", 25, step = 0.5),
    numericInput(ns("sd_T"), "SD of annual temperature (°C)", 0.5, min = 0.01, step = 0.05),
    numericInput(ns("mean_P"), "Long-run mean precipitation (mm a year)", 1000, min = 0, step = 50),
    numericInput(ns("sd_P"), "SD of annual precipitation (mm)", 150, min = 1, step = 10),
    numericInput(ns("income"), "GDP per capita", 5000, min = 100, step = 500),
    uiOutput(ns("income_units")),
    numericInput(ns("share"), "Agricultural share of value added (%)", 15, min = 0, max = 100, step = 1),
    conditionalPanel(
      "input.gdp == 'PWT'", ns = ns,
      numericInput(ns("capital"), "Capital stock per capita", 20000, min = 100, step = 1000),
      helpText(paste0("Long-run (sample mean) level, in ", EXPLORER$capital_units, "; used with the capital ",
                      "interaction."))
    )
  )
}

# The values a country preset fills in (rounded as shown in the inputs; capital
# NA without PWT capital).
preset_values <- function(row) {
  c(mean_T = round(row$mean_T, 2), sd_T = round(row$sd_T, 3), mean_P = round(1000 * row$mean_P),
    sd_P = round(1000 * row$sd_P), income = signif(row$income, 3), share = round(100 * row$share, 1),
    capital = signif(row$capital, 3))
}

# Setting n: its presets, the copy of setting 1's model, and two reactives:
# whether it is shown, and its specification (label, model, country).
setting_server <- function(id, n, default_country, first = NULL) {
  moduleServer(id, function(input, output, session) {
    panel_key <- reactive({
      req(input$gdp, input$climate, input$baseline)
      paste(input$gdp, input$climate, input$baseline, input$weight %||% "population",
            input$units %||% "SD", sep = " | ")
    })
    units <- reactive(EXPLORER$units[[panel_key()]])

    observeEvent(panel_key(), {
      u <- units()
      current <- isolate(input$preset)
      selected <- if (is.null(current)) default_country else current
      if (!selected %in% c("", u$GID_0)) selected <- ""
      choices <- c("Custom" = "", stats::setNames(u$GID_0, u$country)[order(u$country)])
      updateSelectInput(session, "preset", choices = choices, selected = selected)
    })
    observeEvent(list(input$preset, panel_key()), {
      req(input$preset)
      row <- units()[units()$GID_0 == input$preset, ]
      req(nrow(row) == 1)
      v <- preset_values(row)
      for (k in names(v)[is.finite(v)]) updateNumericInput(session, k, value = v[[k]])
    }, ignoreInit = TRUE)
    # The PWT outcomes and the capital interaction only with PWT GDP.
    observeEvent(input$gdp, {
      pwt <- input$gdp == "PWT"
      modifiers <- if (pwt) MODIFIERS else MODIFIERS[!MODIFIERS %in% PWT_MODIFIERS]
      for (k in c("inter_T", "inter_P")) {
        updateCheckboxGroupInput(session, k, choices = modifiers, inline = TRUE,
                                 selected = intersect(isolate(input[[k]]) %||% DEFAULT_MODIFIERS, modifiers))
      }
      outcomes <- if (pwt) OUTCOME_CHOICES else OUTCOME_CHOICES[!OUTCOME_CHOICES %in% PWT_OUTCOMES]
      current <- isolate(input$outcome) %||% "log"
      updateSelectInput(session, "outcome", choices = outcomes,
                        selected = if (current %in% outcomes) current else "log")
    })
    if (!is.null(first)) {
      observeEvent(input$copy, {
        m <- first$model_inputs()
        for (k in c("gdp", "climate", "outcome", "baseline", "form", "fe")) {
          updateSelectInput(session, k, selected = m[[k]])
        }
        updateRadioButtons(session, "moderators", selected = m$moderators)
        for (v in c("T", "P")) {
          updateRadioButtons(session, paste0("income_form_", v), selected = m[[paste0("income_form_", v)]])
          updateSliderInput(session, paste0("xi_", v), value = m[[paste0("xi_", v)]])
        }
        updateCheckboxInput(session, "share_trailing", value = isTRUE(m$share_trailing))
        for (k in c("inter_T", "inter_P")) {
          updateCheckboxGroupInput(session, k, selected = m[[k]] %||% character(0))
        }
        updateRadioButtons(session, "weight", selected = m$weight)
        updateRadioButtons(session, "units", selected = m$units)
        updateRadioButtons(session, "sd_mode", selected = m$sd_mode)
        updateRadioButtons(session, "clim_mode", selected = m$clim_mode)
        updateCheckboxInput(session, "cap", value = isTRUE(m$cap))
        updateRadioButtons(session, "shape", selected = m$shape)
      })
    }
    output$income_units <- renderUI({
      helpText(paste0("Long-run (sample mean) level, in ", EXPLORER$income_units[[input$gdp]], "."))
    })
    # The characteristics estimated as trailing means are those the country has
    # as its climate changes: by default they follow the baseline in the paths.
    observeEvent(input$moderators, {
      updateRadioButtons(session, "clim_mode",
                         selected = if (input$moderators == "trailing") "baseline" else "fixed")
    }, ignoreInit = TRUE)
    # The sample of the trailing share against the trailing income alone.
    output$share_warning <- renderUI({
      req(input$gdp == "WB")
      r <- trailing_rows(input$gdp, input$climate, input$weight %||% "population", EXPLORER)
      div(class = "alert alert-warning small py-1 px-2",
          sprintf(paste("The trailing share needs all 5 years of WDI value added before each year:",
                        "%s rows and %d countries, against %s rows and %d countries without it (-%.0f%%),",
                        "mostly before 1990."),
                  format(r["share", "rows"], big.mark = ","), r["share", "countries"],
                  format(r["trailing", "rows"], big.mark = ","), r["trailing", "countries"],
                  100 * (1 - r["share", "rows"] / r["trailing", "rows"])))
    })

    model_inputs <- reactive(list(gdp = input$gdp, climate = input$climate, outcome = input$outcome,
                                  baseline = input$baseline,
                                  weight = input$weight, units = input$units, fe = input$fe,
                                  form = input$form, inter_T = input$inter_T, inter_P = input$inter_P,
                                  moderators = input$moderators, share_trailing = input$share_trailing,
                                  income_form_T = input$income_form_T, income_form_P = input$income_form_P,
                                  xi_T = input$xi_T, xi_P = input$xi_P,
                                  sd_mode = input$sd_mode, clim_mode = input$clim_mode,
                                  cap = input$cap, shape = input$shape))
    active <- reactive(n == 1 || isTRUE(input$active))
    spec <- reactive({
      req(input$gdp, input$climate, input$baseline, input$form)
      req(is.finite(input$mean_T), isTRUE(input$sd_T > 0), is.finite(input$mean_P),
          isTRUE(input$sd_P > 0), isTRUE(input$income > 0), isTRUE(input$share >= 0))
      # The PWT outcomes and capital only with PWT (model_options()); an
      # unticked group is NULL: no interaction for that variable.
      checked <- model_options(list(gdp = input$gdp, outcome = input$outcome %||% "log",
                                    inter = list(TM = intersect(MODIFIERS, input$inter_T %||% character(0)),
                                                 RR = intersect(MODIFIERS, input$inter_P %||% character(0))),
                                    income_form = c(TM = input$income_form_T %||% "standard",
                                                    RR = input$income_form_P %||% "standard"),
                                    xi = c(TM = input$xi_T %||% 0, RR = input$xi_P %||% 0)))
      outcome <- checked$outcome
      inter <- checked$inter
      income <- checked
      uses_kap <- "kap" %in% unlist(inter)
      req(!uses_kap || isTRUE(input$capital > 0))
      fe <- input$fe %||% "quadratic"
      u <- units()
      # The country's capital only with the capital interaction.
      values <- c(mean_T = input$mean_T, sd_T = input$sd_T, mean_P = input$mean_P,
                  sd_P = input$sd_P, income = input$income, share = input$share,
                  capital = if (uses_kap) input$capital)
      row <- u[u$GID_0 == (input$preset %||% ""), ]
      name <- if (nrow(row) == 1) {
        if (isTRUE(all.equal(unname(values), unname(preset_values(row)[names(values)])))) row$country
        else paste(row$country, "(edited)")
      } else "Custom country"
      trailing <- input$baseline != "Horizon"
      weight <- input$weight %||% "population"
      absolute <- identical(input$units, "absolute")
      moderators <- input$moderators %||% "fixed"
      share_trailing <- moderators == "trailing" && input$gdp == "WB" && isTRUE(input$share_trailing)
      # Absolute deviations have no SD.
      sd_mode <- if (trailing && !absolute) input$sd_mode else "fixed"
      # With the horizon baseline the climate stays at the horizon mean, but
      # trailing income and share can still follow the SSP.
      clim_mode <- if (trailing || moderators == "trailing") input$clim_mode %||% "fixed" else "fixed"
      # The cap matters only with an income, share or capital interaction.
      cap <- isTRUE(input$cap) && any(c("inc", "agr", "kap") %in% unlist(inter))
      # Splines matter only with an interaction.
      shape <- if (length(unlist(inter))) input$shape %||% "linear" else "linear"
      # The label names the options only where they differ from the default.
      auto <- paste0(name, " · ", input$gdp, " × ", input$climate, " · ", tolower(input$baseline),
                     " · ", FORM_SHORT[[input$form]], " · ", interaction_label(inter),
                     if (outcome != "log") paste0(" · outcome: ", OUTCOME_SHORT[[outcome]]),
                     if (!is.null(income_form_label(income$income_form, income$xi))) {
                       paste0(" · ", income_form_label(income$income_form, income$xi))
                     },
                     if (fe != "quadratic") paste0(" · ", FE_SHORT[[fe]]),
                     if (moderators == "trailing") " · trailing characteristics",
                     if (share_trailing) " (with share)",
                     if (weight == "area") " · area weights",
                     if (absolute) " · absolute units",
                     if (shape == "splines") " · splines",
                     if (cap) " · capped",
                     if (sd_mode == "moving") " · moving SD",
                     if (clim_mode == "baseline") " · characteristics follow baseline")
      list(n = n, label = if (nzchar(trimws(input$label))) trimws(input$label) else auto,
           model = explorer_model(list(gdp = input$gdp, climate = input$climate, weight = weight,
                                       units = input$units %||% "SD", baseline = input$baseline,
                                       form = input$form, inter = inter, shape = shape, fe = fe,
                                       moderators = moderators, share_trailing = share_trailing,
                                       outcome = outcome, income_form = income$income_form, xi = income$xi),
                                  EXPLORER),
           country = list(mean_T = input$mean_T, sd_T = input$sd_T, mean_P = input$mean_P / 1000,
                          sd_P = input$sd_P / 1000, income = input$income, share = input$share / 100,
                          capital = if (uses_kap) input$capital),
           sd_mode = sd_mode, clim_mode = clim_mode, cap = cap,
           support = EXPLORER$support[[paste(panel_key(), fe, sep = " | ")]], units = u, gdp = input$gdp)
    })
    list(active = active, spec = spec, model_inputs = model_inputs)
  })
}

setting_tab <- function(n) {
  nav_panel(
    title = tags$span(tags$span("●", style = paste0("color:", SETTING_COLOURS[n], ";")),
                      paste("Setting", n)),
    value = paste0("s", n),
    setting_ui(paste0("s", n), n)
  )
}

ui <- page_sidebar(
  title = "Weather response explorer",
  # Scroll rather than squeeze the charts into the window's height (a short
  # window would leave them too small to draw).
  fillable = FALSE,
  theme = bs_theme(version = 5, bg = SURFACE, fg = INK, primary = SETTING_COLOURS[1]),
  sidebar = sidebar(
    width = 360,
    tags$h6("Permanent shift (common to all settings)"),
    numericInput("shock_T", "Temperature (°C)", 1, step = 0.5),
    numericInput("shock_P", "Precipitation (mm a year)", 100, step = 25),
    sliderInput("horizon", "Years shown in the GDP paths", 10, 100, 50, step = 5),
    numericInput("ramp_years", paste("Year X: the gradual shift reaches its full size (dashed) and",
                                     "the sudden late one starts (dotted)"), 50, min = 1, step = 5),
    checkboxGroupInput("gdp_lines", "Shifts shown in the GDP paths",
                       c("Immediate (solid)" = "permanent", "Gradual until year X (dashed)" = "gradual",
                         "Sudden in year X (dotted)" = "delayed"),
                       selected = c("permanent", "gradual", "delayed")),
    tags$hr(),
    do.call(navset_pill, c(lapply(seq_len(N_SETTINGS), setting_tab), list(id = "setting_tabs")))
  ),
  uiOutput("warnings"),
  layout_columns(
    col_widths = c(6, 6),
    card(card_header("Temperature: growth response to an annual anomaly"),
         plotOutput("response_T", height = 380, hover = hoverOpts("hover_response_T", delay = 80)),
         uiOutput("readout_response_T")),
    card(card_header("Precipitation: growth response to an annual anomaly"),
         plotOutput("response_P", height = 380, hover = hoverOpts("hover_response_P", delay = 80)),
         uiOutput("readout_response_P"))
  ),
  layout_columns(
    col_widths = c(6, 6),
    card(card_header(textOutput("gdp_header_T", inline = TRUE)),
         plotOutput("gdp_T", height = 430, hover = hoverOpts("hover_gdp_T", delay = 80)),
         uiOutput("readout_gdp_T")),
    card(card_header(textOutput("gdp_header_P", inline = TRUE)),
         plotOutput("gdp_P", height = 430, hover = hoverOpts("hover_gdp_P", delay = 80)),
         uiOutput("readout_gdp_P"))
  ),
  card(card_header("Each setting's model: specification, coefficients, fit and joint tests"),
       withMathJax(uiOutput("spec_math")),
       tableOutput("coef_table"),
       uiOutput("coef_notes")),
  card(
    card_header("Climate scenarios: GDP losses under each setting's model, 2025-2100"),
    layout_columns(
      col_widths = c(6, 6),
      checkboxGroupInput("scen_sel", "Forcing: the global warming path", SCENARIO_CHOICES,
                         selected = SCENARIO_CHOICES, inline = TRUE),
      radioButtons("econ_ssp", paste("Economic development: income, population and agricultural",
                                     "share (Follow SSP-RCP: each forcing with its own SSP)"),
                   ECONOMY_CHOICES, selected = "follow", inline = TRUE)
    ),
    layout_columns(
      col_widths = c(4, 3, 5),
      selectInput("map_setting", "Maps: setting", choices = NULL),
      selectInput("map_scenario", "Maps: forcing", choices = NULL),
      radioButtons("map_component", "GDP maps: impact of",
                   c("T and P" = "loss", "T only" = "loss_T", "P only" = "loss_P"),
                   selected = "loss", inline = TRUE)
    ),
    div(class = "d-flex align-items-center gap-3 mb-2",
        actionButton("run_scenarios", "Run the climate scenarios", class = "btn-primary"),
        uiOutput("scenario_status")),
    layout_columns(
      col_widths = c(4, 8),
      plotOutput("gmt_plot", height = 380),
      plotOutput("global_plot", height = 380, hover = hoverOpts("hover_global", delay = 80))
    ),
    uiOutput("readout_global"),
    layout_columns(
      col_widths = c(6, 6),
      plotOutput("map_2050", height = 340),
      plotOutput("map_2100", height = 340)
    ),
    tags$h6("Physical drivers: the change of each country's climate in the scenario", class = "mt-3"),
    layout_columns(
      col_widths = c(6, 6),
      plotOutput("driver_T_2050", height = 320),
      plotOutput("driver_T_2100", height = 320)
    ),
    layout_columns(
      col_widths = c(6, 6),
      plotOutput("driver_P_2050", height = 320),
      plotOutput("driver_P_2100", height = 320)
    ),
    tableOutput("global_table"),
    uiOutput("scenario_notes")
  ),
  card(
    card_header("Termination shocks: perfect geoengineering, 2025-2100"),
    layout_columns(
      col_widths = c(3, 3, 6),
      selectInput("ts_forcing", "Background forcing (the warming without geoengineering)",
                  SCENARIO_CHOICES, selected = "ssp245"),
      radioButtons("ts_held", "Geoengineering holds global warming at",
                   c("1.5 C" = "1.5", "1.0 C" = "1.0"), selected = "1.5", inline = TRUE),
      checkboxGroupInput("ts_years", "Termination years", as.character(seq(2040, 2090, by = 10)),
                         selected = c("2040", "2060", "2080"), inline = TRUE)
    ),
    layout_columns(
      col_widths = c(3, 9),
      numericInput("ts_tau", "Catch-up time after termination, tau (years)", 4, min = 0.5, step = 0.5),
      radioButtons("ts_econ", "Economic development (Follow SSP-RCP: the forcing's own SSP)",
                   ECONOMY_CHOICES, selected = "follow", inline = TRUE)
    ),
    div(class = "d-flex align-items-center gap-3 mb-2",
        actionButton("run_termination", "Run the termination scenarios", class = "btn-primary"),
        uiOutput("ts_status")),
    layout_columns(
      col_widths = c(4, 8),
      plotOutput("ts_gmt_plot", height = 380),
      plotOutput("ts_global_plot", height = 380, hover = hoverOpts("hover_ts", delay = 80))
    ),
    uiOutput("readout_ts"),
    layout_columns(
      col_widths = c(6, 6),
      selectInput("ts_map_setting", "Map: setting", choices = NULL),
      selectInput("ts_map_year", "Map: termination year", choices = NULL)
    ),
    plotOutput("ts_map", height = 380),
    tableOutput("ts_table"),
    uiOutput("ts_notes")
  ),
  card(card_header("Models and assumptions"), uiOutput("notes"))
)

server <- function(input, output, session) {
  settings <- list(setting_server("s1", 1, DEFAULT_COUNTRIES[1]))
  settings <- c(settings, lapply(2:N_SETTINGS, function(n) {
    setting_server(paste0("s", n), n, DEFAULT_COUNTRIES[n], first = settings[[1]])
  }))

  # The shown settings, with unique labels.
  specs <- reactive({
    out <- lapply(Filter(function(s) s$active(), settings), function(s) s$spec())
    labels <- vapply(out, `[[`, "", "label")
    dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    for (i in which(dup)) out[[i]]$label <- paste0(labels[i], " (setting ", out[[i]]$n, ")")
    out
  })
  # The colour of each shown setting by its label (in setting order).
  colours <- reactive({
    sp <- specs()
    stats::setNames(SETTING_COLOURS[vapply(sp, `[[`, 0, "n")], vapply(sp, `[[`, "", "label"))
  })
  as_setting <- function(x) factor(x, levels = names(colours()))
  shift <- function(v) if (v == "TM") input$shock_T else input$shock_P / 1000
  to_display <- function(v) if (v == "TM") 1 else 1000
  sd_of <- function(s, v) if (v == "TM") s$country$sd_T else s$country$sd_P

  # Responses over -3..3 SD of each setting's country (wider if the shift is
  # larger), in physical units, flagged where the deviation is within the 1st-
  # 99th percentile of those in the data.
  responses <- function(v) {
    reactive({
      delta <- shift(v)
      req(is.finite(delta))
      do.call(rbind, lapply(specs(), function(s) {
        sd <- sd_of(s, v)
        reach <- max(3, 1.1 * abs(delta) / sd)
        x <- seq(-reach, reach, length.out = 241)
        # The deviations in the model's units (SD, or C and m).
        x_model <- x * sd / deviation_scale(s$model, s$country, v)
        q <- s$support[s$support$variable == v, ]
        d <- response_curve(s$model, v, x_model, s$country, s$cap)
        d$x <- x
        d$anomaly <- x * sd * to_display(v)
        d$in_support <- x_model >= q$q01 & x_model <= q$q99
        d$setting <- s$label
        d
      }))
    })
  }
  responses_T <- responses("TM")
  responses_P <- responses("RR")
  shift_points <- function(v) {
    do.call(rbind, lapply(specs(), function(s) {
      d <- response_curve(s$model, v, shift(v) / deviation_scale(s$model, s$country, v), s$country, s$cap)
      d$anomaly <- shift(v) * to_display(v)
      d$setting <- s$label
      d
    }))
  }

  # The three scenarios of the shift: from year 0 on (solid), rising linearly
  # to its full size in year X (dashed), or nothing until year X and the full
  # shift from it on (dotted); the names are the types of gdp_path().
  # Only the ticked ones are drawn.
  scenarios <- reactive({
    all <- c(permanent = "immediate",
             gradual = sprintf("gradual until year %d", as.integer(input$ramp_years)),
             delayed = sprintf("sudden in year %d", as.integer(input$ramp_years)))
    shown <- all[names(all) %in% input$gdp_lines]
    validate(need(length(shown) > 0, "Tick at least one shift (immediate, gradual or sudden) to draw."))
    shown
  })
  # The line type and band transparency of each scenario, by its type.
  SCENARIO_LINETYPE <- c(permanent = "solid", gradual = "42", delayed = "13")
  SCENARIO_ALPHA <- c(permanent = 0.15, gradual = 0.07, delayed = 0.07)
  paths <- function(v) {
    reactive({
      delta <- shift(v)
      req(is.finite(delta), input$horizon, isTRUE(input$ramp_years >= 1))
      do.call(rbind, lapply(specs(), function(s) {
        do.call(rbind, lapply(names(scenarios()), function(type) {
          cbind(setting = s$label, scenario = scenarios()[[type]],
                gdp_path(s$model, v, s$country, type, delta, input$horizon, s$sd_mode,
                         s$clim_mode, s$cap, ramp_years = input$ramp_years))
        }))
      }))
    })
  }
  paths_T <- paths("TM")
  paths_P <- paths("RR")

  response_plot <- function(d, v) {
    d$setting <- as_setting(d$setting)
    points <- shift_points(v)
    points$setting <- as_setting(points$setting)
    unit <- if (v == "TM") "°C" else "mm"
    ggplot(d, aes(anomaly, 100 * estimate, colour = setting, fill = setting)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_vline(xintercept = 0, colour = GRID, linewidth = 0.4) +
      geom_ribbon(aes(ymin = 100 * low, ymax = 100 * high), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.6, linetype = "22") +
      geom_line(aes(y = ifelse(in_support, 100 * estimate, NA)), linewidth = 0.9, na.rm = TRUE) +
      geom_point(data = points, size = 3.2, shape = 21, fill = SURFACE, stroke = 1.4) +
      scale_colour_manual(values = colours(), name = NULL, labels = wrap_label) +
      scale_fill_manual(values = colours(), guide = "none") +
      labs(x = paste0("Annual anomaly against the long-run baseline (", unit, ")"),
           y = "Effect on growth that year (pp)",
           subtitle = paste0("Bands: 95% CI. Dashed: beyond the 1st-99th percentile of the ",
                             "deviations in the data. Points: the shift.")) +
      theme_explorer()
  }
  output$response_T <- renderPlot(response_plot(responses_T(), "TM"), bg = SURFACE)
  output$response_P <- renderPlot(response_plot(responses_P(), "RR"), bg = SURFACE)

  readout_list <- function(lines) tags$ul(class = "small mb-0", lapply(lines, tags$li))
  response_readout <- function(d, hover, v) {
    if (is.null(hover)) return(helpText("Hover over the chart to read the values."))
    unit <- if (v == "TM") "°C" else "mm"
    readout_list(lapply(split(d, as_setting(d$setting)), function(s) {
      i <- which.min(abs(s$anomaly - hover$x))
      sprintf("%s: %+.2f %s (%+.2f SD): %+.2f pp [%+.2f, %+.2f]", s$setting[i], s$anomaly[i], unit,
              s$x[i], 100 * s$estimate[i], 100 * s$low[i], 100 * s$high[i])
    }))
  }
  output$readout_response_T <- renderUI(response_readout(responses_T(), input$hover_response_T, "TM"))
  output$readout_response_P <- renderUI(response_readout(responses_P(), input$hover_response_P, "RR"))

  output$gdp_header_T <- renderText(sprintf(
    "Effect on %s of a permanent %+g °C shift in temperature: immediate, gradual or sudden in year %d",
    quantity_of(specs()), input$shock_T, as.integer(input$ramp_years)))
  output$gdp_header_P <- renderText(sprintf(
    "Effect on %s of a permanent %+g mm shift in precipitation: immediate, gradual or sudden in year %d",
    quantity_of(specs()), input$shock_P, as.integer(input$ramp_years)))

  gdp_plot <- function(d) {
    d$setting <- as_setting(d$setting)
    d$scenario <- factor(d$scenario, levels = unname(scenarios()))
    x <- as.integer(input$ramp_years)
    captions <- c(permanent = "Solid: the full shift from year 0 on.",
                  gradual = sprintf("Dashed: rising linearly to the full shift in year %d.", x),
                  delayed = sprintf("Dotted: nothing until year %d, the full shift from then on.", x))
    ggplot(d, aes(year, pct, colour = setting, fill = setting, linetype = scenario,
                  group = interaction(setting, scenario))) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = pct_low, ymax = pct_high, alpha = scenario), colour = NA,
                  show.legend = FALSE) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = colours(), name = NULL, labels = wrap_label) +
      scale_fill_manual(values = colours(), guide = "none") +
      scale_linetype_manual(values = stats::setNames(SCENARIO_LINETYPE[names(scenarios())], scenarios()),
                            name = NULL) +
      scale_alpha_manual(values = stats::setNames(SCENARIO_ALPHA[names(scenarios())], scenarios()),
                         guide = "none") +
      guides(colour = guide_legend(order = 1),
             linetype = guide_legend(order = 2, direction = "horizontal", nrow = 1,
                                     override.aes = list(colour = INK_MUTED))) +
      labs(x = "Years after the shift starts (0 = first year)", y = paste(first_upper(quantity_of(specs())), "against no shift (%)"),
           subtitle = paste(strwrap(paste(c(captions[names(scenarios())], "Bands: 95% CI."),
                                          collapse = " "), width = 110), collapse = "\n")) +
      theme_explorer() +
      theme(legend.box = "vertical", legend.box.just = "left")
  }
  output$gdp_T <- renderPlot(gdp_plot(paths_T()), bg = SURFACE)
  output$gdp_P <- renderPlot(gdp_plot(paths_P()), bg = SURFACE)

  gdp_readout <- function(d, hover) {
    if (is.null(hover)) return(helpText("Hover over the chart to read the values."))
    year <- max(0, min(round(hover$x), max(d$year)))
    rows <- d[d$year == year, ]
    rows <- rows[order(as_setting(rows$setting), match(rows$scenario, scenarios())), ]
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(sprintf("%s, %s: %+.2f%% [%+.2f, %+.2f]", rows$setting, rows$scenario,
                                 rows$pct, rows$pct_low, rows$pct_high)))
  }
  output$readout_gdp_T <- renderUI(gdp_readout(paths_T(), input$hover_gdp_T))
  output$readout_gdp_P <- renderUI(gdp_readout(paths_P(), input$hover_gdp_P))

  # The equation of each shown setting's model and the definitions of its terms.
  output$spec_math <- renderUI({
    blocks <- lapply(specs(), function(s) {
      tags$div(
        class = "mb-3",
        tags$div(tags$span("●", style = paste0("color:", SETTING_COLOURS[s$n], ";")),
                 tags$strong(paste0(" Setting ", s$n, ": ")), s$label),
        tags$div(style = "overflow-x: auto;", paste0("$$", spec_equation(s$model), "$$")),
        tags$ul(class = "small", lapply(spec_definitions(s$model), function(x) tags$li(HTML(x))))
      )
    })
    withMathJax(tagList(blocks))
  })

  # The coefficients of each shown setting's model (rows: every term of any of
  # them; deviation terms, then the level terms of trailing characteristics),
  # its fit, the verdict against setting 1, and the joint Wald tests.
  output$coef_table <- renderTable({
    sp <- specs()
    models <- lapply(sp, `[[`, "model")
    cell <- function(est, se, p) {
      sprintf("%.3g%s<br><span style='color:%s'>(%.3g)</span>", est, stars(p), INK_MUTED, se)
    }
    coef_cells <- lapply(models, function(m) {
      se <- sqrt(diag(m$vcov))
      out <- stats::setNames(cell(m$coef, se, m$p), names(m$coef))
      if (nrow(m$levels)) out <- c(out, stats::setNames(cell(m$levels$estimate, m$levels$se, m$levels$p),
                                                        paste0("level:", m$levels$term)))
      out
    })
    outcome_label <- function(m) paste0("growth of ", OUTCOME_SHORT[[m$outcome %||% "log"]], " (log points)")
    terms <- unique(unlist(lapply(coef_cells, names)))
    terms <- c(terms[!startsWith(terms, "level:")], terms[startsWith(terms, "level:")])
    test_names <- unique(unlist(lapply(models, function(m) m$tests$test)))
    rows <- c("Setting", "Outcome", terms, "Observations", "Countries", "Within adj. R2", "AIC", "BIC",
              "Against setting 1 (AIC / BIC / R2)", test_names)
    columns <- lapply(seq_along(models), function(i) {
      m <- models[[i]]
      tests <- stats::setNames(sprintf("p = %.3g%s<br><span style='color:%s'>(F = %.2f, %d terms)</span>",
                                       m$tests$p, stars(m$tests$p), INK_MUTED, m$tests$stat,
                                       as.integer(m$tests$df1)), m$tests$test)
      verdict <- if (i == 1) "" else record_verdicts(m, models[[1]]) %||% "other rows or outcome: not comparable"
      values <- c(Setting = htmltools::htmlEscape(sp[[i]]$label), Outcome = outcome_label(m), coef_cells[[i]],
                  Observations = format(m$nobs, big.mark = ","), Countries = as.character(m$countries),
                  `Within adj. R2` = sprintf("%.4f", m$fit[["war2"]]), AIC = sprintf("%.1f", m$fit[["aic"]]),
                  BIC = sprintf("%.1f", m$fit[["bic"]]), `Against setting 1 (AIC / BIC / R2)` = verdict, tests)
      out <- values[rows]
      out[is.na(out)] <- ""
      unname(out)
    })
    names(columns) <- vapply(sp, function(s) {
      sprintf("<span style='color:%s'>●</span> Setting %d", SETTING_COLOURS[s$n], s$n)
    }, "")
    labels <- rows
    labels[rows %in% terms] <- vapply(rows[rows %in% terms], term_label, "")
    labels[rows %in% test_names] <- paste("Joint test,", rows[rows %in% test_names])
    data.frame(Term = labels, columns, check.names = FALSE)
  }, striped = TRUE, spacing = "s", align = "l", sanitize.text.function = function(x) x)
  output$coef_notes <- renderUI(tags$ul(
    class = "small",
    tags$li(paste(
      "Coefficients: effect on the outcome (growth of GDP per capita, capital per capita or GDP per",
      "worker, in log points) per unit of the arm (an SD, or a C or m, beyond the baseline; squared",
      "for quadratic arms), and per unit of the characteristic for the interactions (long-run T in C",
      "and P in m a year less the median of the countries; log income, share and log capital per",
      "capita less their means over the rows; spline terms as in the notes below).",
      "Clustered SEs in brackets; *** p < 0.01, ** p < 0.05, * p < 0.1.")),
    tags$li(paste(
      "Level terms: the trailing characteristics on their own (controls; not in the responses",
      "and projections).")),
    tags$li(paste(
      "Joint tests: clustered Wald tests (F) that the terms are all zero: all the deviation terms",
      "of T (P), the terms interacted with each characteristic (both arms, with their spline",
      "terms), each pair of characteristics when there are three, and all the interactions.")),
    tags$li(paste(
      "AIC and BIC compare only models on the same rows. Against setting 1: better, worse, or",
      "indifferent when the difference cannot be told apart from noise (Vuong test clustered by",
      "country, 5%), for each criterion; on other rows (another panel, or trailing values that",
      "drop rows) or with another outcome the models are not compared."))
  ))

  # Warnings for each setting: characteristics outside the countries of its
  # sample, and a shift beyond the deviations in the data.
  output$warnings <- renderUI({
    msgs <- character(0)
    for (s in specs()) {
      cn <- s$country
      u <- s$units
      m <- s$model
      prefix <- paste0(s$label, ": ")
      check <- function(value, range, label, fmt) {
        if (value < min(range) || value > max(range)) {
          msgs <<- c(msgs, sprintf(paste0("%s%s (", fmt, ") is outside the range of the sample's ",
                                          "countries (", fmt, " to ", fmt, ")."),
                                   prefix, label, value, min(range), max(range)))
        }
      }
      used <- unique(unlist(m$inter))
      if ("clim" %in% m$inter$TM) check(cn$mean_T, u$mean_T, "mean temperature", "%.1f °C")
      if ("clim" %in% m$inter$RR) check(1000 * cn$mean_P, 1000 * u$mean_P, "mean precipitation", "%.0f mm")
      if ("inc" %in% used) check(cn$income, u$income, "GDP per capita", "%.0f")
      if ("agr" %in% used) check(100 * cn$share, 100 * u$share, "agricultural share", "%.1f%%")
      if ("kap" %in% used) check(cn$capital, u$capital[is.finite(u$capital)], "capital per capita", "%.0f")
      for (v in c("TM", "RR")) {
        # The shift and the data's deviations in the model's units.
        x <- shift(v) / deviation_scale(m, cn, v)
        q <- s$support[s$support$variable == v, ]
        if (x < q$q01 || x > q$q99) {
          unit <- if (!identical(m$units, "absolute")) "SD" else if (v == "TM") "C" else "m"
          msgs <- c(msgs, sprintf(paste(
            "%sthe %s shift is %+.2f %s, beyond the 1st-99th percentile of the deviations in the",
            "data (%+.2f to %+.2f %s): its response is an extrapolation."),
            prefix, if (v == "TM") "temperature" else "precipitation", x, unit, q$q01, q$q99, unit))
        }
      }
      if (m$baseline == "Horizon") {
        msgs <- c(msgs, paste0(prefix, "with the horizon baseline the permanent shift stays a ",
                               "deviation every year, so its GDP effect grows without limit."))
      }
      if (isTRUE(m$share_trailing)) {
        r <- trailing_rows(m$gdp, m$climate, m$weight, EXPLORER)
        msgs <- c(msgs, sprintf(paste(
          "%sthe trailing agricultural share needs all 5 years of WDI value added before each year:",
          "the model runs on %s rows (%d countries), %.0f%% fewer than with the trailing income alone",
          "(%s rows), mostly before 1990, so its sample differs from the other models'."),
          prefix, format(m$nobs, big.mark = ","), m$countries,
          100 * (1 - m$nobs / r["trailing", "rows"]), format(r["trailing", "rows"], big.mark = ",")))
      }
    }
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-warning", style = "margin-bottom: 0.5rem;",
        tags$ul(style = "margin-bottom: 0;", lapply(msgs, tags$li)))
  })

  # Climate scenarios --------------------------------------------------------
  # The projections of a setting depend on its model and options only (every
  # country enters with its own characteristics), so they are cached by those.
  projection_cache <- new.env()
  projection_of <- function(s, scenario, ssp) {
    key <- paste(s$model$key, s$sd_mode, s$clim_mode, s$cap, scenario, ssp, sep = " | ")
    if (is.null(projection_cache[[key]])) {
      projection_cache[[key]] <- project_setting(s, SCEN, scenario, ssp)
    }
    projection_cache[[key]]
  }
  # The projections run only when asked: the shown settings, the selected
  # forcings and the economic development at the moment of the click.
  NOT_RUN <- "Press \"Run the climate scenarios\" to project the global losses and the maps."
  run <- eventReactive(input$run_scenarios, {
    validate(need(length(input$scen_sel) > 0, "Select at least one forcing."))
    list(specs = specs(), scenarios = input$scen_sel, economy = input$econ_ssp %||% "follow")
  })
  run_key <- function(sp, scenarios, economy) {
    paste(c(vapply(sp, function(s) paste(s$label, s$model$key, s$sd_mode, s$clim_mode, s$cap), ""),
            scenarios, economy), collapse = " || ")
  }
  # The label of a forcing of the last run, with its economy.
  run_label <- function(scenario) scenario_label(scenario, run()$economy)
  output$scenario_status <- renderUI({
    if (input$run_scenarios == 0) return(helpText(class = "mb-0", "Not run yet."))
    r <- run()
    current <- run_key(specs(), input$scen_sel, input$econ_ssp %||% "follow")
    if (!identical(current, run_key(r$specs, r$scenarios, r$economy))) {
      return(helpText(class = "mb-0 text-warning-emphasis",
                      "Settings or scenarios changed since the last run: press Run to update."))
    }
    helpText(class = "mb-0", sprintf("Up to date: %d setting(s) x %d forcing(s), economy: %s.",
                                     length(r$specs), length(r$scenarios),
                                     names(ECONOMY_CHOICES)[ECONOMY_CHOICES == r$economy]))
  })
  # The colours of the settings of the last run (in setting order).
  run_colours <- reactive({
    sp <- run()$specs
    stats::setNames(SETTING_COLOURS[vapply(sp, `[[`, 0, "n")], vapply(sp, `[[`, "", "label"))
  })
  as_run_setting <- function(x) factor(x, levels = names(run_colours()))
  # Every setting and scenario of the last run: list of (setting, scenario, result).
  projections <- reactive({
    r <- run()
    jobs <- expand.grid(i = seq_along(r$specs), scenario = r$scenarios, stringsAsFactors = FALSE)
    withProgress(message = "Projecting the scenarios", value = 0, {
      lapply(seq_len(nrow(jobs)), function(j) {
        s <- r$specs[[jobs$i[j]]]
        scenario <- jobs$scenario[j]
        incProgress(1 / nrow(jobs), detail = paste(s$label, "-", scenario_label(scenario, r$economy)))
        list(setting = s$label, n = s$n, scenario = scenario,
             label = scenario_label(scenario, r$economy),
             result = projection_of(s, scenario, economy_of(scenario, r$economy)))
      })
    })
  })
  global_losses <- reactive({
    do.call(rbind, lapply(projections(), function(p) {
      cbind(setting = p$setting, scenario = p$scenario, label = p$label, p$result$global)
    }))
  })
  # The maps offer the settings and scenarios of the last run.
  observeEvent(run(), {
    r <- run()
    choices <- stats::setNames(as.character(vapply(r$specs, `[[`, 0, "n")),
                               vapply(r$specs, `[[`, "", "label"))
    selected <- isolate(input$map_setting)
    if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
    updateSelectInput(session, "map_setting", choices = choices, selected = selected)
    scenario_choices <- SCENARIO_CHOICES[SCENARIO_CHOICES %in% r$scenarios]
    names(scenario_choices) <- vapply(scenario_choices, scenario_label, "", r$economy)
    selected <- isolate(input$map_scenario)
    if (is.null(selected) || !selected %in% scenario_choices) {
      selected <- if ("ssp245" %in% scenario_choices) "ssp245" else scenario_choices[[1]]
    }
    updateSelectInput(session, "map_scenario", choices = scenario_choices, selected = selected)
  })

  output$gmt_plot <- renderPlot({
    req(length(input$scen_sel) > 0)
    g <- SCEN$gmt[SCEN$gmt$scenario %in% input$scen_sel & SCEN$gmt$year >= 1970, ]
    g$scenario <- factor(g$scenario, levels = names(SCEN$scenarios))
    ggplot(g, aes(year, median, colour = scenario, fill = scenario)) +
      geom_vline(xintercept = PROJECTION_FROM, colour = GRID, linewidth = 0.5) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.10, colour = NA) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = SCENARIO_COLOURS, labels = FORCING_LABELS, name = NULL,
                          guide = guide_legend(direction = "horizontal", nrow = 1)) +
      scale_fill_manual(values = SCENARIO_COLOURS, guide = "none") +
      labs(x = NULL, y = "Global warming above 1850-1900 (°C)",
           subtitle = "FaIR (AR6-calibrated): median and 5-95% range") +
      theme_explorer()
  }, bg = SURFACE)

  output$global_plot <- renderPlot({
    validate(need(input$run_scenarios > 0, NOT_RUN))
    d <- global_losses()
    d$setting <- as_run_setting(d$setting)
    d$scenario <- factor(d$scenario, levels = names(SCEN$scenarios))
    labels <- vapply(names(SCEN$scenarios), run_label, "")
    ggplot(d, aes(year, estimate, colour = scenario, fill = scenario)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.10, colour = NA) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~setting, ncol = min(3, nlevels(d$setting)), scales = "free_y",
                 labeller = labeller(setting = function(x) wrap_label(x, 50))) +
      scale_colour_manual(values = SCENARIO_COLOURS, labels = labels, name = NULL,
                          guide = guide_legend(direction = "horizontal", nrow = 1)) +
      scale_fill_manual(values = SCENARIO_COLOURS, guide = "none") +
      labs(x = NULL, y = paste("World", quantity_of(run()$specs), "against no climate change (%)"),
           subtitle = paste("Best guess (median warming, CMIP6-median pattern, estimated coefficients);",
                            "bands: 95% CI of the coefficients")) +
      theme_explorer() +
      theme(strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0))
  }, bg = SURFACE)

  output$readout_global <- renderUI({
    hover <- input$hover_global
    if (input$run_scenarios == 0) return(NULL)
    if (is.null(hover)) return(helpText("Hover over the global chart to read the values."))
    d <- global_losses()
    year <- max(PROJECTION_FROM, min(round(hover$x), PROJECTION_TO))
    rows <- d[d$year == year, ]
    rows <- rows[order(as_run_setting(rows$setting), match(rows$scenario, names(SCEN$scenarios))), ]
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(sprintf("%s, %s: %+.2f%% [%+.2f, %+.2f]", rows$setting, rows$label,
                                 rows$estimate, rows$low, rows$high)))
  })

  # Best-guess country losses of the map setting and scenario in `year`, with
  # one colour scale for the two maps.
  map_data_of <- reactive({
    validate(need(input$run_scenarios > 0, NOT_RUN))
    req(input$map_setting, input$map_scenario)
    s <- Filter(function(x) as.character(x$n) == input$map_setting, run()$specs)
    req(length(s) == 1)
    r <- projection_of(s[[1]], input$map_scenario, economy_of(input$map_scenario, run()$economy))
    list(setting = s[[1]]$label, label = run_label(input$map_scenario),
         countries = r$countries[r$countries$year %in% c(2050, 2100), ])
  })
  theme_map <- function() {
    theme_void(base_size = 12) +
      theme(plot.background = element_rect(fill = SURFACE, colour = NA),
            plot.title = element_text(colour = INK, face = "bold"),
            plot.subtitle = element_text(colour = INK_MUTED, size = 9.5),
            legend.position = "bottom", legend.key.width = grid::unit(1.6, "cm"),
            legend.title = element_text(colour = INK_MUTED, size = 9.5),
            legend.text = element_text(colour = INK_MUTED, size = 9))
  }
  # The world map with `value` of each country (iso3), Antarctica left out.
  world_with <- function(iso3, value) {
    world <- SCEN$world[SCEN$world$iso3 != "ATA", ]
    world$value <- value[match(world$iso3, iso3)]
    world
  }
  COMPONENT_TITLES <- c(loss = "impact of T and P", loss_T = "impact of T only",
                        loss_P = "impact of P only")
  # GDP impact of the chosen component (T and P, T only, P only), with one
  # colour scale for the two years.
  country_map <- function(year) {
    m <- map_data_of()
    component <- input$map_component %||% "loss"
    values <- m$countries[m$countries$year == year, ]
    limit <- max(abs(m$countries[[component]]), na.rm = TRUE)
    ggplot(world_with(values$iso3, values[[component]])) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      scale_fill_gradient2(low = MAP_LOSS, mid = MAP_ZERO, high = MAP_GAIN, midpoint = 0,
                           limits = c(-limit, limit), na.value = MAP_NA,
                           name = paste0(first_upper(quantity_of(run()$specs)), " against no\nclimate change (%)")) +
      labs(title = sprintf("%d: %s, %s", year, m$label, COMPONENT_TITLES[[component]]),
           subtitle = wrap_label(m$setting)) +
      theme_map()
  }
  output$map_2050 <- renderPlot(country_map(2050), bg = SURFACE)
  output$map_2100 <- renderPlot(country_map(2100), bg = SURFACE)

  # Physical drivers: each country's T change (C) and P change (% of its
  # long-run P) in 2050 and 2100 against the mean of the map setting's horizon,
  # the CMIP6 pattern times the global warming (as in the projections).
  driver_data <- reactive({
    validate(need(input$run_scenarios > 0, NOT_RUN))
    req(input$map_setting, input$map_scenario)
    s <- Filter(function(x) as.character(x$n) == input$map_setting, run()$specs)
    req(length(s) == 1)
    g <- SCEN$gmt[SCEN$gmt$scenario == input$map_scenario, ]
    gmt <- stats::setNames(g$median, g$year)
    horizon <- s[[1]]$model$horizon %||% c(1970L, 2023L)
    reference <- mean(gmt[as.character(seq(horizon[[1]], horizon[[2]]))])
    weight <- s[[1]]$model$weight %||% "population"
    p <- SCEN$pattern[[weight]]
    do.call(rbind, lapply(c(2050, 2100), function(y) {
      warming <- gmt[[as.character(y)]] - reference
      data.frame(iso3 = p$iso3, year = y, warming = warming,
                 dT = p$t_slope * warming, dP = 100 * p$p_rel_slope * warming)
    })) |>
      transform(label = FORCING_LABELS[[input$map_scenario]], weight = weight,
                reference = sprintf("%d-%d", horizon[[1]], horizon[[2]]))
  })
  driver_map <- function(variable, year) {
    d <- driver_data()
    values <- d[d$year == year, ]
    title <- sprintf("%d: %s change, %s (global +%.2f C)", year,
                     if (variable == "dT") "temperature" else "precipitation",
                     values$label[[1]], values$warming[[1]])
    subtitle <- sprintf(paste("Against the mean of %s; CMIP6 pattern (median over the models, %s",
                              "weights) times the FaIR median warming"),
                        values$reference[[1]], values$weight[[1]])
    base <- ggplot(world_with(values$iso3, values[[variable]])) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      labs(title = title, subtitle = wrap_label(subtitle, 90)) +
      theme_map()
    if (variable == "dT") {
      base + scale_fill_gradientn(colours = SCEN$palettes$temp_seq,
                                  limits = c(min(0, d$dT, na.rm = TRUE), max(d$dT, na.rm = TRUE)),
                                  na.value = MAP_NA, name = "Temperature\nchange (C)")
    } else {
      # A few arid countries have very large relative changes: the colours stop
      # at the 95th percentile of the absolute changes (beyond: the end colours).
      limit <- stats::quantile(abs(d$dP), 0.95, na.rm = TRUE, names = FALSE)
      base + scale_fill_gradientn(colours = SCEN$palettes$prec_div, limits = c(-limit, limit),
                                  oob = scales::squish, na.value = MAP_NA,
                                  name = sprintf("Precipitation change (%%)\n(colours capped at +/-%.0f%%)",
                                                 limit))
    }
  }
  output$driver_T_2050 <- renderPlot(driver_map("dT", 2050), bg = SURFACE)
  output$driver_T_2100 <- renderPlot(driver_map("dT", 2100), bg = SURFACE)
  output$driver_P_2050 <- renderPlot(driver_map("dP", 2050), bg = SURFACE)
  output$driver_P_2100 <- renderPlot(driver_map("dP", 2100), bg = SURFACE)

  output$global_table <- renderTable({
    if (input$run_scenarios == 0) return(NULL)
    out <- do.call(rbind, lapply(projections(), function(p) {
      g <- p$result$global
      cell <- function(y) {
        r <- g[g$year == y, ]
        sprintf("%+.2f [%+.2f, %+.2f]", r$estimate, r$low, r$high)
      }
      followed <- if (p$result$follow) {
        sprintf("%d / %d", p$result$followed[["income"]], p$result$followed[["share"]])
      } else "fixed"
      data.frame(Setting = p$setting, Scenario = p$label,
                 `Forcing (warming of)` = sprintf("%s (%s)", FORCING_LABELS[[p$scenario]],
                                                  SCEN$scenarios[[p$scenario]]),
                 `Economy (income, population, share)` = p$result$ssp,
                 `World GDP 2050 (%)` = cell(2050), `World GDP 2100 (%)` = cell(2100),
                 Countries = p$result$n,
                 `Following the SSP (income / share)` = followed,
                 `Share of the SSP's world GDP covered` = sprintf("%.0f%%", 100 * p$result$coverage),
                 check.names = FALSE)
    }))
    names(out) <- sub("^World GDP", paste("World", quantity_of(run()$specs)), names(out))
    out[order(match(out$Setting, names(run_colours())),
              match(out$`Forcing (warming of)`, sprintf("%s (%s)", FORCING_LABELS, SCEN$scenarios))), ]
  }, striped = TRUE, spacing = "s", align = "l")

  output$scenario_notes <- renderUI({
    tags$ul(
      class = "small",
      tags$li(paste0(
        "Global warming: ", SCEN$source, "; temperature above 1850-1900, median and 5-95% over the ",
        "ensemble. The projections use the median path.")),
      tags$li(paste(
        "Country climate: the forced signal only, no year-to-year variability. The anomaly is the",
        "global warming above its mean over the panel's horizon, times the country's CMIP6 pattern",
        "(median over the models of the best-fit linear relationship with global temperature; T:",
        "the slope, C per C; P: the slope relative to each model's own mean P, times the country's",
        "long-run P).")),
      tags$li(paste(
        "Damages: each setting's model and options (baseline, arms, interactions, shape, cap, SD and",
        "long-run characteristics) applied to every country of its panel with the country's own",
        "long-run climate, SD, income and share; losses accumulate from 2025 (no deviation from then",
        "on in the counterfactual). The custom country of a setting does not enter.")),
      tags$li(paste0(
        "Scenarios combine a forcing, which gives the global warming, with an economic development ",
        "(SSP1-5), which gives income, population (the GDP weights) and the agricultural share. With ",
        "'Follow SSP-RCP' each forcing takes its own SSP: SSP1 (1.9, 2.6), SSP2 (4.5), SSP3 (7.0), ",
        "SSP4 (3.4, 6.0), SSP5 (8.5); otherwise the chosen SSP goes with every forcing (labels like ",
        "'4.5 W/m2 x SSP3'). The warming of a forcing is that of its own emissions scenario. ",
        SCEN$ssp_source, ". With the long-run characteristics following a trailing baseline, the ",
        "interactions use the trailing 20 (30) year means of the long-run climate, of income and of ",
        "the agricultural share (for a model estimated with trailing characteristics, the means of ",
        "income and, if estimated so, share over the 5 years before, also with the horizon baseline; ",
        "its fixed share stays at the sample mean). Income and share are the country's historical series in the model's ",
        "definitions (GDP per capita of the panel's GDP source; agriculture over the three sectors' ",
        "value added) continued with the SSP's relative change from their last year. SSP GDP per ",
        "capita and population (5-year points) and the SSP share (5-year points 2030-2050) are joined ",
        "by log-linear interpolation; after 2050 the share follows SSP income with the country's ",
        "2030-2050 elasticity. Fixed: the sample means, as in estimation. The table counts the ",
        "countries that follow the SSP in income and in share (the others keep their sample means).")),
      tags$li(paste(
        "Global loss: the average of the countries' losses weighted by the SSP's GDP of each year",
        "(PPP), over the countries with a pattern and SSP GDP; the share of the SSP's world GDP they",
        "cover in 2025 is in the table. 95% CI: 1,000 draws of the coefficients, the same for every",
        "country. Not included: the uncertainty of the warming and of the CMIP6 pattern."))
    )
  })

  # Termination shocks ---------------------------------------------------------
  # Global warming paths with perfect geoengineering: none (the background),
  # sustained (held from 2025), and terminated in each year (termination_path()).
  # They run only when asked, for the shown settings at the moment of the click.
  TS_NOT_RUN <- "Press \"Run the termination scenarios\" to project the termination shocks."
  ts_run <- eventReactive(input$run_termination, {
    validate(need(length(input$ts_years) > 0, "Select at least one termination year."),
             need(isTRUE(input$ts_tau > 0), "The catch-up time must be positive."))
    list(specs = specs(), forcing = input$ts_forcing, held = as.numeric(input$ts_held),
         years = sort(as.integer(input$ts_years)), tau = input$ts_tau,
         economy = input$ts_econ %||% "follow")
  })
  ts_key <- function(sp, forcing, held, years, tau, economy) {
    paste(c(vapply(sp, function(s) paste(s$label, s$model$key, s$sd_mode, s$clim_mode, s$cap), ""),
            forcing, held, years, tau, economy), collapse = " || ")
  }
  output$ts_status <- renderUI({
    if (input$run_termination == 0) return(helpText(class = "mb-0", "Not run yet."))
    r <- ts_run()
    current <- ts_key(specs(), input$ts_forcing, as.numeric(input$ts_held),
                      sort(as.integer(input$ts_years)), input$ts_tau, input$ts_econ %||% "follow")
    if (!identical(current, ts_key(r$specs, r$forcing, r$held, r$years, r$tau, r$economy))) {
      return(helpText(class = "mb-0 text-warning-emphasis",
                      "Settings or options changed since the last run: press Run to update."))
    }
    helpText(class = "mb-0", sprintf("Up to date: %d setting(s), %s, held at %.1f C, %d termination year(s).",
                                     length(r$specs), scenario_label(r$forcing, r$economy), r$held,
                                     length(r$years)))
  })
  # The paths of the last run: name -> global warming (named by year).
  ts_paths <- reactive({
    r <- ts_run()
    g <- SCEN$gmt[SCEN$gmt$scenario == r$forcing, ]
    background <- stats::setNames(g$median, g$year)
    held <- sprintf("%.1f C", r$held)
    paths <- c(stats::setNames(list(background), "no geoengineering"),
               stats::setNames(list(termination_path(background, r$held)),
                               paste("sustained at", held)))
    for (y in r$years) {
      paths[[paste("terminated in", y)]] <- termination_path(background, r$held, y, r$tau)
    }
    paths
  })
  # The colours of the paths: grey without geoengineering, blue sustained, an
  # orange ramp for the termination years (darker later).
  ts_colours <- reactive({
    p <- names(ts_paths())
    n_term <- length(p) - 2
    stats::setNames(c(INK_MUTED, SETTING_COLOURS[[1]],
                      grDevices::colorRampPalette(c("#f3a272", "#eb6834", "#8f3510"))(n_term)), p)
  })
  # Every setting and path of the last run: list of (setting, path, result).
  ts_projections <- reactive({
    r <- ts_run()
    paths <- ts_paths()
    ssp <- economy_of(r$forcing, r$economy)
    jobs <- expand.grid(i = seq_along(r$specs), path = names(paths), stringsAsFactors = FALSE)
    withProgress(message = "Projecting the termination scenarios", value = 0, {
      lapply(seq_len(nrow(jobs)), function(j) {
        s <- r$specs[[jobs$i[j]]]
        path <- jobs$path[j]
        incProgress(1 / nrow(jobs), detail = paste(s$label, "-", path))
        key <- paste(s$model$key, s$sd_mode, s$clim_mode, s$cap, r$forcing, ssp, r$held, r$tau, path,
                     sep = " | ")
        if (is.null(projection_cache[[key]])) {
          projection_cache[[key]] <- project_setting(s, SCEN, r$forcing, ssp, gmt = paths[[path]])
        }
        list(setting = s$label, n = s$n, path = path, result = projection_cache[[key]])
      })
    })
  })
  ts_global <- reactive({
    do.call(rbind, lapply(ts_projections(), function(p) {
      cbind(setting = p$setting, path = p$path, p$result$global)
    }))
  })
  observeEvent(ts_run(), {
    r <- ts_run()
    choices <- stats::setNames(as.character(vapply(r$specs, `[[`, 0, "n")),
                               vapply(r$specs, `[[`, "", "label"))
    updateSelectInput(session, "ts_map_setting", choices = choices, selected = choices[[1]])
    updateSelectInput(session, "ts_map_year", choices = as.character(r$years),
                      selected = as.character(r$years[[1]]))
  })
  ts_run_colours <- function(x) factor(x, levels = names(ts_colours()))

  output$ts_gmt_plot <- renderPlot({
    validate(need(input$run_termination > 0, TS_NOT_RUN))
    paths <- ts_paths()
    d <- do.call(rbind, lapply(names(paths), function(k) {
      data.frame(path = k, year = as.integer(names(paths[[k]])), warming = unname(paths[[k]]))
    }))
    d <- d[d$year >= 2000, ]
    d$path <- ts_run_colours(d$path)
    ggplot(d, aes(year, warming, colour = path)) +
      geom_vline(xintercept = PROJECTION_FROM, colour = GRID, linewidth = 0.5) +
      geom_line(aes(linetype = path == "no geoengineering"), linewidth = 0.9) +
      scale_colour_manual(values = ts_colours(), name = NULL,
                          guide = guide_legend(ncol = 2)) +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
      labs(x = NULL, y = "Global warming above 1850-1900 (C)",
           subtitle = sprintf("Background: %s (FaIR median); catch-up time %.1f years",
                              FORCING_LABELS[[ts_run()$forcing]], ts_run()$tau)) +
      theme_explorer()
  }, bg = SURFACE)

  output$ts_global_plot <- renderPlot({
    validate(need(input$run_termination > 0, TS_NOT_RUN))
    d <- ts_global()
    sp <- ts_run()$specs
    d$setting <- factor(d$setting, levels = vapply(sp, `[[`, "", "label"))
    d$path <- ts_run_colours(d$path)
    ggplot(d, aes(year, estimate, colour = path, fill = path)) +
      geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.4) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.06, colour = NA) +
      geom_line(aes(linetype = path == "no geoengineering"), linewidth = 0.9) +
      facet_wrap(~setting, ncol = min(3, nlevels(d$setting)), scales = "free_y",
                 labeller = labeller(setting = function(x) wrap_label(x, 50))) +
      scale_colour_manual(values = ts_colours(), name = NULL, guide = guide_legend(ncol = 3)) +
      scale_fill_manual(values = ts_colours(), guide = "none") +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
      labs(x = NULL, y = paste("World", quantity_of(sp), "against no climate change (%)"),
           subtitle = sprintf("%s; bands: 95%% CI of the coefficients",
                              scenario_label(ts_run()$forcing, ts_run()$economy))) +
      theme_explorer() +
      theme(strip.text = element_text(colour = INK_MUTED, size = 9.5, hjust = 0))
  }, bg = SURFACE)

  output$readout_ts <- renderUI({
    if (input$run_termination == 0) return(NULL)
    hover <- input$hover_ts
    if (is.null(hover)) return(helpText("Hover over the global chart to read the values."))
    d <- ts_global()
    year <- max(PROJECTION_FROM, min(round(hover$x), PROJECTION_TO))
    rows <- d[d$year == year, ]
    rows <- rows[order(match(rows$path, names(ts_colours()))), ]
    tagList(tags$div(class = "small", paste("Year", year)),
            readout_list(sprintf("%s, %s: %+.2f%% [%+.2f, %+.2f]", rows$setting, rows$path,
                                 rows$estimate, rows$low, rows$high)))
  })

  # Map: the extra loss from terminating in the chosen year against sustained
  # geoengineering, by country in 2100 (pp of GDP; negative: a larger loss).
  output$ts_map <- renderPlot({
    validate(need(input$run_termination > 0, TS_NOT_RUN))
    req(input$ts_map_setting, input$ts_map_year)
    pr <- Filter(function(p) as.character(p$n) == input$ts_map_setting, ts_projections())
    sustained <- Filter(function(p) startsWith(p$path, "sustained"), pr)[[1]]$result$countries
    terminated <- Filter(function(p) p$path == paste("terminated in", input$ts_map_year), pr)
    req(length(terminated) == 1)
    terminated <- terminated[[1]]$result$countries
    a <- terminated[terminated$year == PROJECTION_TO, ]
    b <- sustained[sustained$year == PROJECTION_TO, ]
    difference <- a$loss - b$loss[match(a$iso3, b$iso3)]
    limit <- max(abs(difference), na.rm = TRUE)
    ggplot(world_with(a$iso3, difference)) +
      geom_sf(aes(fill = value), colour = SURFACE, linewidth = 0.1) +
      scale_fill_gradient2(low = MAP_LOSS, mid = MAP_ZERO, high = MAP_GAIN, midpoint = 0,
                           limits = c(-limit, limit), na.value = MAP_NA,
                           name = paste0(first_upper(quantity_of(ts_run()$specs)), ", terminated minus\nsustained (pp)")) +
      labs(title = sprintf("%d: termination in %s against sustained geoengineering (%.1f C)",
                           PROJECTION_TO, input$ts_map_year, ts_run()$held),
           subtitle = wrap_label(paste0(pr[[1]]$setting, "; ",
                                        scenario_label(ts_run()$forcing, ts_run()$economy)))) +
      theme_map()
  }, bg = SURFACE)

  output$ts_table <- renderTable({
    if (input$run_termination == 0) return(NULL)
    pr <- ts_projections()
    cell <- function(g, y) {
      r <- g[g$year == y, ]
      sprintf("%+.2f [%+.2f, %+.2f]", r$estimate, r$low, r$high)
    }
    out <- do.call(rbind, lapply(pr, function(p) {
      g <- p$result$global
      sustained <- Filter(function(q) q$setting == p$setting && startsWith(q$path, "sustained"), pr)[[1]]
      s100 <- sustained$result$global$estimate[sustained$result$global$year == PROJECTION_TO]
      data.frame(Setting = p$setting, Path = p$path,
                 `World GDP 2050 (%)` = cell(g, 2050), `World GDP 2100 (%)` = cell(g, PROJECTION_TO),
                 `2100 against sustained (pp)` = sprintf("%+.2f", g$estimate[g$year == PROJECTION_TO] - s100),
                 check.names = FALSE)
    }))
    names(out) <- sub("^World GDP", paste("World", quantity_of(ts_run()$specs)), names(out))
    out[order(match(out$Setting, vapply(ts_run()$specs, `[[`, "", "label")),
              match(out$Path, names(ts_colours()))), ]
  }, striped = TRUE, spacing = "s", align = "l")

  output$ts_notes <- renderUI({
    tags$ul(
      class = "small",
      tags$li(paste(
        "Perfect geoengineering: a climate with the same global temperature is the same with or",
        "without geoengineering, so each country's climate is its CMIP6 pattern times global warming,",
        "as in the climate scenarios. Background: the FaIR median warming of the chosen forcing.")),
      tags$li(paste(
        "Sustained: from 2025 global warming is held at 1.5 C or 1.0 C (or the background where it is",
        "lower; at 1.0 C this is an immediate cooling from today's warming of about 1.4 C, itself a",
        "deviation for the damage model). Terminated in year T: held until T, then the gap to the",
        "background closes exponentially, G(t) = B(t) - max(B(T) - H, 0) exp(-(t - T) / tau).")),
      tags$li(paste(
        "Damages: each setting's model and options, every country of its panel, losses accumulating",
        "from 2025 against no climate change, and the global loss weighted by the SSP's GDP, as in the",
        "climate scenarios. The map and the last column compare termination with sustained",
        "geoengineering in 2100 (negative: termination costs more)."))
    )
  })

  output$notes <- renderUI({
    model_lines <- lapply(specs(), function(s) {
      m <- s$model
      # The splines against the linear interactions (same rows): the Wald test
      # of the non-linear terms and the verdicts.
      test <- if (!is.null(m$spline_test)) {
        sprintf(paste(" Splines against linear interactions: non-linear terms jointly p = %.3g;",
                      "AIC %+.1f, BIC %+.1f; verdict by AIC / BIC / within adj. R2: %s."),
                m$spline_test$nonlinear_p, m$spline_test$aic, m$spline_test$bic,
                m$spline_test$verdict)
      } else ""
      # The level terms of the trailing characteristics (controls: not in the
      # responses or projections).
      levels <- if (nrow(m$levels)) {
        sprintf(" Level terms (controls, not in the responses and projections): %s.",
                paste(sprintf("%s %+.4f (SE %.4f, p = %.3g)", m$levels$term, m$levels$estimate,
                              m$levels$se, m$levels$p), collapse = "; "))
      } else ""
      moderators <- if (m$moderators == "trailing") {
        paste0("trailing characteristics (climate: ",
               if (m$window == "all") "horizon mean" else paste0("trailing ", sub("lag_", "", m$window), "-year mean"),
               "; income: 5 years before; share: ", if (m$share_trailing) "5 years before" else "sample mean", ")")
      } else "fixed characteristics (sample means)"
      tags$li(tags$span("●", style = paste0("color:", SETTING_COLOURS[s$n], ";")),
              sprintf(paste("%s: %s x %s, %s baseline, %s, year and unit effects with %s, interactions: %s",
                            "(%s) on %s; %s rows, %s countries from %d.%s%s"), s$label, m$gdp, m$climate,
                      tolower(m$baseline), tolower(m$form), FE_SHORT[[m$fe]], interaction_label(m$inter),
                      m$shape %||% "linear", moderators, format(m$nobs, big.mark = ","), m$countries,
                      EXPLORER$start_year, test, levels))
    })
    tags$div(
      tags$ul(model_lines),
      tags$p(paste(
        "Every model: year and unit effects, with no unit trends, unit linear trends, or unit linear",
        "and quadratic trends (the notebook's); the deviations are demeaned in the same way before",
        "their arms are cut (each unit's mean, mean and linear trend, or mean and linear and",
        "quadratic trend, fitted from 1960); conflict and",
        "collapse controls; errors clustered by country. Characteristics in the interactions: fixed",
        "(the sample means: long-run climate, log income and share, constant within a country) or",
        "trailing (the mean climate of the deviation's own baseline, the trailing 20 or 30 years, or",
        "the horizon mean with the horizon baseline; log of the mean GDP per capita of the 5 years",
        "before, from the whole GDP history of the source; with WB GDP optionally the mean",
        "agricultural share of the 5 years before). Trailing characteristics in an interaction also",
        "enter on their own (level terms): without them the interaction of a non-negative arm",
        "would stand in for the omitted level. Trailing models run on the rows with every trailing",
        "value they use. The climate is centred at the median of the countries' long-run climate,",
        "income and share at their means over the model's rows. Climate weights: population density (2000)",
        "or area (plain area average) in each country's climate, in the estimation and in the CMIP6",
        "downscaling of the scenarios (its population-weighted or area-weighted fits); the same rows",
        "and growth-outlier exclusions. Deviations in SD of the country's variability, or in absolute",
        "units (C and metres, demeaned in the same way; the SD then plays no part in the model and",
        "the option of a moving SD does not apply). T and P terms are estimated together, each",
        "with its own interactions; a characteristic neither variable is interacted with is not",
        "used. Splines: each interacted characteristic enters as a restricted cubic spline (knots at",
        "the 5, 35, 65 and 95th percentiles of the panel's countries, linear beyond the outer knots),",
        "with its non-linear terms set to zero at the median climate and mean income and share.",
        "Verdicts: 'indifferent' when the difference cannot be told apart from noise (the Wald test",
        "of the non-linear terms, 5%).")),
      tags$ul(
        tags$li(paste(
          "Response: the effect on that year's growth of a year whose temperature (precipitation) is",
          "the given anomaly away from the baseline, in SDs of the country's variability, at its",
          "long-run climate, income and share. Estimated on deviations demeaned as the fixed effects.")),
        tags$li(paste(
          "Outcome, in log points: the growth of GDP per capita; with PWT also the growth of the capital",
          "stock per capita (rnna / pop) or of GDP per worker (rgdpna / emp, persons engaged), on the rows",
          "with them. Responses, paths and losses are of the outcome: with capital or GDP per worker the",
          "GDP paths and scenario losses are of capital per capita or GDP per worker (weighted by the",
          "SSP's GDP in the global loss), as the axis labels say.")),
        tags$li(paste(
          "Capital interaction (PWT): log of the capital stock per capita (rnna / pop), long-run (mean",
          "over the horizon) or trailing (mean of the 5 years before, with its level term), less its",
          "mean over the rows. Its cross-country correlation with log income is about 0.92 (0.46 within",
          "countries after the quadratic trends). In the GDP paths and scenarios it stays at the country's",
          "value (there is no SSP capital path); the scenarios use the countries with PWT capital.")),
        tags$li(paste(
          "Income scaling, for the temperature and the precipitation terms separately:",
          "standard, or every deviation term (each arm and its ticked interactions) times (Ybar / Y)^(1 - xi),",
          "with Y the country's income of the interactions and Ybar a constant, the geometric mean of that",
          "income over the rows, so that the loss in dollars per person is",
          "proportional to Y^xi beyond the ticked interactions (xi = 0: equal dollars; xi = 1: equal %, the",
          "standard terms; xi = 2: the % loss rises with income), or both the plain terms and their scaled",
          "copies (nests the two; the table tests each). Y is the income of the interactions: long-run, or",
          "trailing (then with its level term, and along the SSP in the scenarios). The cap applies to the",
          "plain terms and to the scaled terms separately.")),
        tags$li(paste(
          "GDP path: growth effects accumulate in log GDP. The shift is immediate (solid: its full",
          "size from year 0 on), gradual (dashed: rising linearly from zero in year 0 to its full",
          "size in year X, then constant) or sudden and late (dotted: nothing until year X, the full",
          "shift from then on). Horizon baseline: the shift is a",
          "deviation every year, so growth is lower (or higher) every year. Trailing baselines: the",
          "shift enters the trailing mean, so an immediate shift fades out over 20 (30) years and a",
          "gradual one stays a deviation of about the warming of the last 10 (15) years while it",
          "rises, then fades out once it stops. The SD is fixed at the",
          "country's, or widened by the shift (moving). The long-run climate of the interactions is",
          "fixed at the country's sample mean, or follows the trailing baseline (the country's mean",
          "plus the shift as it enters the trailing mean: the response adapts as if the country's",
          "climate were that of the baseline); in the GDP paths income and share stay at the",
          "country's (they follow the SSP only in the climate scenarios). With fixed characteristics",
          "the models are estimated with the interactions on the sample (horizon) means, so following",
          "the baseline is an assumption of the paths; with trailing characteristics it is how they",
          "were estimated (the default then), and in the scenarios income and share follow the SSP",
          "as their 5-year trailing means.")),
        tags$li(paste(
          "Cap: each arm's response is its own coefficient c (the response at the sample's mean",
          "income and share), plus the income and share terms, plus the climate term. With the cap,",
          "the non-climate part (c plus the income and share terms) is at most max(c, 0), whatever",
          "the income and share, however high or low: income and share can bring a damaging response",
          "to zero but not past it, and cannot raise a beneficial one. The climate term, at the",
          "country's historical climate and any warming beyond it, is added on top, uncapped. The",
          "capped response is not linear in the coefficients, so its 95% CI comes from 2,000 draws",
          "of the coefficients (percentiles), capped draw by draw.")),
        tags$li(paste(
          "Held fixed: income and share of the interactions (no feedback from the GDP losses to",
          "income), the long-run climate unless it follows the baseline, and the absence of lags",
          "(no recovery).")),
        tags$li(paste0(
          "Country SDs: horizon SD of the annual series, or for the trailing baselines the median of ",
          "the trailing SD over the country's years. Agricultural share: World Bank agriculture value ",
          "added over that of the three sectors (constant 2015 US$), mean over the sample years. ",
          "Built ", format(EXPLORER$built, "%Y-%m-%d %H:%M"), "."))
      )
    )
  })
}

app <- shinyApp(ui, server)
# source()d in an interactive session: launch it; loaded by runApp(): return it.
if (interactive() && !is.null(SOURCED_FILE)) runApp(app) else app
