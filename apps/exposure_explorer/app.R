# Exposure explorer: how exposure to extreme days (ERA5 daily, population
# weights 2015) maps into the annual mean temperature / annual precipitation,
# across units (climatologies) and within units (year-to-year), for countries
# (GADM0) or regions (GADM1). Physical side only, no GDP.
#
#   Rscript apps/exposure_explorer/build_exposure.R   # once (histograms)
#   shiny::runApp("apps/exposure_explorer")           # or source() this file
#
# R/exposure.R holds the computations.

library(shiny)
library(bslib)
library(ggplot2)

SOURCED_FILE <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
APP_DIR <- local({
  candidates <- c(if (!is.null(SOURCED_FILE)) dirname(SOURCED_FILE), getwd(),
                  file.path(getwd(), "apps", "exposure_explorer"))
  found <- candidates[file.exists(file.path(candidates, "R", "exposure.R"))]
  if (!length(found)) stop("Cannot find the app's folder (apps/exposure_explorer).")
  normalizePath(found[[1]])
})
source(file.path(APP_DIR, "R", "exposure.R"), local = environment())
DATA_DIR <- normalizePath(file.path(APP_DIR, "..", "..", "results", "exposure_explorer"), mustWork = FALSE)
RESOLUTIONS <- c("Countries (GADM0)" = "gadm0", "Regions (GADM1)" = "gadm1")
available <- RESOLUTIONS[file.exists(file.path(DATA_DIR, paste0(RESOLUTIONS, "_tmp.rds")))]
if (!length(available)) stop("No histograms in ", DATA_DIR, ": run Rscript apps/exposure_explorer/build_exposure.R")

# Histograms, read once per resolution x variable.
CACHE <- new.env()
hist_of <- function(res, var) {
  key <- paste(res, var)
  if (is.null(CACHE[[key]])) CACHE[[key]] <- readRDS(file.path(DATA_DIR, paste0(res, "_", var, ".rds")))
  CACHE[[key]]
}
H0 <- hist_of(available[[1]], "tmp")
YEAR_RANGE <- range(H0$years)
COUNTRY_CHOICES <- local({
  u <- H0$units[!duplicated(H0$units$GID_0), ]
  stats::setNames(u$GID_0, paste0(u$country, " (", u$GID_0, ")"))[order(u$country)]
})
P_EDGES <- hist_of(available[[1]], "pre")$lo
VARS <- list(
  tmp = list(label = "Temperature", annual = "annual mean T (C)", unit = "C",
             hi = "hot", lo = "cold", groups = "10, 20", combine = "raw"),
  pre = list(label = "Precipitation", annual = "annual precipitation (mm)", unit = "mm",
             hi = "wet", lo = "dry", groups = "500, 1500", combine = "sd"))
# Daily series of the temperature indices (the annual value stays the mean of
# the daily mean): daily max / min only where their histograms were built.
T_SOURCES_ALL <- c("Daily mean" = "tmp", "Daily max" = "tmpmax", "Daily min" = "tmpmin")
has_file <- function(res, var) file.exists(file.path(DATA_DIR, paste0(res, "_", var, ".rds")))
T_SOURCES <- T_SOURCES_ALL[vapply(T_SOURCES_ALL, function(v) any(vapply(available, has_file, TRUE, var = v)), TRUE)]
T_SHORT <- c(tmp = "daily mean", tmpmax = "Tmax", tmpmin = "Tmin")
# Slider ranges (the bin edges of each file, less the open first bin) and the
# default thresholds by source.
T_RANGE <- list(tmp = c(-29, 40), tmpmax = c(-34, 50), tmpmin = c(-44, 40))
T_DEFAULT <- list(hot = c(tmp = 25, tmpmax = 29, tmpmin = 20), cold = c(tmp = 5, tmpmax = 5, tmpmin = 0))
HOT_SRC0 <- if ("tmpmax" %in% T_SOURCES) "tmpmax" else "tmp"     # default sources: Tmax hot, Tmin cold
COLD_SRC0 <- if ("tmpmin" %in% T_SOURCES) "tmpmin" else "tmp"
theme_set(theme_minimal(base_size = 12))
PCT_LABEL <- function(p) {
  v <- 100 * p
  suf <- ifelse(v %% 1 != 0 | v %% 100 %in% 11:13, "th",
                c("th", "st", "nd", "rd", rep("th", 6))[(floor(v) %% 10) + 1])
  paste0(format(v, trim = TRUE), suf, " pct")
}

ui <- page_sidebar(
  title = "Exposure explorer: extreme days and the annual mean",
  sidebar = sidebar(
    width = 330,
    radioButtons("res", "Units", available, inline = TRUE),
    radioButtons("var", "Variable", c("Temperature" = "tmp", "Precipitation" = "pre"), inline = TRUE),
    sliderInput("years", "Years", min = YEAR_RANGE[1] + 10, max = YEAR_RANGE[2],
                value = c(1970, YEAR_RANGE[2]), step = 1, sep = ""),
    selectizeInput("countries", "Countries (empty = all)", COUNTRY_CHOICES, multiple = TRUE),
    conditionalPanel("input.var == 'tmp'",
      h6("Hot"),
      radioButtons("hot_src", NULL, T_SOURCES, selected = HOT_SRC0, inline = TRUE),
      radioButtons("hot_type", NULL, c("Absolute" = "abs", "Unit percentile" = "pct"), inline = TRUE),
      conditionalPanel("input.hot_type == 'abs'",
        sliderInput("hot_abs", "Above (C)", T_RANGE[[HOT_SRC0]][1], T_RANGE[[HOT_SRC0]][2], T_DEFAULT$hot[[HOT_SRC0]], step = 1)),
      conditionalPanel("input.hot_type == 'pct'",
        selectInput("hot_pct", "Above the unit's", stats::setNames(H0$pct$hot$probs, PCT_LABEL(H0$pct$hot$probs)), 0.99)),
      radioButtons("hot_measure", NULL, c("Degree days" = "dd", "Days" = "days"), inline = TRUE),
      h6("Cold"),
      radioButtons("cold_src", NULL, T_SOURCES, selected = COLD_SRC0, inline = TRUE),
      radioButtons("cold_type", NULL, c("Absolute" = "abs", "Unit percentile" = "pct"), inline = TRUE),
      conditionalPanel("input.cold_type == 'abs'",
        sliderInput("cold_abs", "Below (C)", T_RANGE[[COLD_SRC0]][1], T_RANGE[[COLD_SRC0]][2], T_DEFAULT$cold[[COLD_SRC0]], step = 1)),
      conditionalPanel("input.cold_type == 'pct'",
        selectInput("cold_pct", "Below the unit's", stats::setNames(H0$pct$cold$probs, PCT_LABEL(H0$pct$cold$probs)), 0.01)),
      radioButtons("cold_measure", NULL, c("Degree days" = "dd", "Days" = "days"), inline = TRUE)
    ),
    conditionalPanel("input.var == 'pre'",
      h6("Wet"),
      radioButtons("wet_type", NULL, c("Absolute" = "abs", "Unit percentile" = "pct"), inline = TRUE),
      conditionalPanel("input.wet_type == 'abs'",
        selectInput("wet_abs", "Daily P above (mm)", P_EDGES[P_EDGES >= 1], 20)),
      conditionalPanel("input.wet_type == 'pct'",
        selectInput("wet_pct", "Daily P above the unit's (all days)",
                    stats::setNames(hist_of(available[[1]], "pre")$pct$wet$probs,
                                    PCT_LABEL(hist_of(available[[1]], "pre")$pct$wet$probs)), 0.99)),
      radioButtons("wet_measure", NULL, c("mm beyond" = "dd", "Days" = "days"), inline = TRUE),
      h6("Dry"),
      selectInput("dry_abs", "Daily P below (mm)", P_EDGES[P_EDGES > 0 & P_EDGES <= 10], 1),
      radioButtons("dry_measure", NULL, c("Days" = "days", "mm short" = "dd"), inline = TRUE)
    ),
    radioButtons("combine", "Sum of the two", c("Raw units" = "raw", "Each over its cross-unit SD" = "sd")),
    textInput("groups", "Climate group breaks (comma separated)", "10, 20")
  ),
  navset_card_tab(
    nav_panel("Across units",
      plotOutput("across_plot", height = "480px"),
      checkboxInput("segments", "Within-unit slopes (over +-2 SD of the annual value)", TRUE),
      tableOutput("across_table")),
    nav_panel("Within units",
      radioButtons("within_scale", NULL, c("Absolute exposure" = "abs", "Standardized (unit SD)" = "std"), inline = TRUE),
      plotOutput("within_plot", height = "480px"),
      h6("Exposure change per unit SD of the annual value, deciles below vs above the mean"),
      tableOutput("arms_table"),
      h6("Within-unit quadratic fits (kappa: slope change over +-1 SD relative to the slope; units with a significant slope)"),
      tableOutput("within_table")),
    nav_panel("Single units",
      selectizeInput("units", "Units (up to 6)", NULL, multiple = TRUE, options = list(maxItems = 6)),
      plotOutput("units_plot", height = "520px"),
      tableOutput("units_table")),
    nav_panel("Trailing deviations",
      radioButtons("window", "Trailing window (years before t)", c(10, 20, 30), selected = 20, inline = TRUE),
      plotOutput("trail_plot", height = "420px"),
      tableOutput("trail_table")),
    nav_panel("About", uiOutput("about"))
  )
)

server <- function(input, output, session) {
  # Defaults that follow the variable.
  observeEvent(input$var, {
    updateTextInput(session, "groups", value = VARS[[input$var]]$groups)
    updateRadioButtons(session, "combine", selected = VARS[[input$var]]$combine)
  })
  # The threshold slider follows the daily series of the index (range and
  # default threshold of that source).
  observeEvent(input$hot_src, {
    r <- T_RANGE[[input$hot_src]]
    updateSliderInput(session, "hot_abs", min = r[1], max = r[2], value = T_DEFAULT$hot[[input$hot_src]])
  }, ignoreInit = TRUE)
  observeEvent(input$cold_src, {
    r <- T_RANGE[[input$cold_src]]
    updateSliderInput(session, "cold_abs", min = r[1], max = r[2], value = T_DEFAULT$cold[[input$cold_src]])
  }, ignoreInit = TRUE)
  H <- reactive(hist_of(input$res, input$var))
  V <- reactive(VARS[[input$var]])
  defs <- reactive({
    if (input$var == "tmp") {
      list(hi = list(side = "above", measure = input$hot_measure, type = input$hot_type, group = "hot", src = input$hot_src,
                     value = if (input$hot_type == "abs") input$hot_abs else input$hot_pct),
           lo = list(side = "below", measure = input$cold_measure, type = input$cold_type, group = "cold", src = input$cold_src,
                     value = if (input$cold_type == "abs") input$cold_abs else input$cold_pct))
    } else {
      list(hi = list(side = "above", measure = input$wet_measure, type = input$wet_type, group = "wet", src = "pre",
                     value = if (input$wet_type == "abs") input$wet_abs else input$wet_pct),
           lo = list(side = "below", measure = input$dry_measure, type = "abs", src = "pre", value = input$dry_abs))
    }
  })
  # The files of the indices' daily series (other than the annual value's).
  sources <- reactive({
    s <- unique(c(defs()$hi$src, defs()$lo$src))
    s <- setdiff(s, input$var)
    validate(need(all(vapply(s, has_file, TRUE, res = input$res)),
                  sprintf("No %s histograms for %s: run Rscript apps/exposure_explorer/build_exposure.R %s %s",
                          paste(T_SHORT[s], collapse = " / "), input$res, input$res, paste(s, collapse = " "))))
    lapply(stats::setNames(s, s), function(v) hist_of(input$res, v))
  })
  # The two indices on every unit and year of H (the annual value's file).
  full <- reactive({
    d <- defs(); src <- sources()
    one <- function(def) if (def$src == input$var) index_matrix(H(), def) else index_aligned(H(), src[[def$src]], def)
    list(hi = one(d$hi), lo = one(d$lo))
  })
  def_label <- function(d, name) {
    thr <- if (d$type == "pct") PCT_LABEL(as.numeric(d$value)) else paste(d$value, V()$unit)
    what <- if (d$measure == "days") "days" else if (input$var == "tmp") "degree days" else "mm"
    of <- if (input$var == "tmp") paste0(" of ", T_SHORT[[d$src]]) else ""
    sprintf("%s: %s%s %s %s", name, what, of, d$side, thr)
  }
  labels <- reactive({
    d <- defs()
    c(hi = def_label(d$hi, V()$hi), lo = def_label(d$lo, V()$lo),
      sum = if (input$combine == "sd") "sum: each over its cross-unit SD" else "sum: raw units")
  })
  yrs <- reactive(seq(input$years[1], input$years[2]))
  units_ok <- reactive(valid_units(H(), yrs(), input$countries, sources()))
  D <- reactive({
    u <- units_ok()
    validate(need(length(u) >= 5, "Fewer than 5 units with complete data for this selection."))
    view_data(H(), full(), u, yrs(), input$combine)
  })
  breaks <- reactive({
    b <- suppressWarnings(as.numeric(strsplit(input$groups, "[,; ]+")[[1]]))
    b[is.finite(b)]
  })
  groups <- reactive(climate_groups(rowMeans(D()$x), breaks(), V()$unit))
  fits <- reactive(lapply(D()$idx, function(y) within_fits(D()$x, y)))
  panel_factor <- function(v) factor(v, c("hi", "lo", "sum"), labels()[c("hi", "lo", "sum")])

  # --- Across units ---
  output$across_plot <- renderPlot({
    d <- D(); xb <- rowMeans(d$x)
    pts <- do.call(rbind, lapply(names(d$idx), function(k)
      data.frame(index = k, x = xb, y = rowMeans(d$idx[[k]]))))
    curves <- do.call(rbind, lapply(names(d$idx), function(k) {
      q <- lm(y ~ x + I(x^2), pts[pts$index == k, ])
      g <- data.frame(index = k, x = seq(min(xb), max(xb), length.out = 200))
      g$y <- predict(q, g); g
    }))
    p <- ggplot(pts, aes(x, y))
    if (input$segments) {
      seg <- do.call(rbind, lapply(names(d$idx), function(k) {
        f <- fits()[[k]]; yb <- rowMeans(d$idx[[k]]); s <- ifelse(is.finite(f$slope), f$slope, 0)
        data.frame(index = k, x0 = xb - 2 * f$sd_x, x1 = xb + 2 * f$sd_x,
                   y0 = yb - 2 * f$sd_x * s, y1 = yb + 2 * f$sd_x * s)
      }))
      p <- p + geom_segment(data = seg, aes(x = x0, xend = x1, y = y0, yend = y1),
                            colour = "steelblue", alpha = if (nrow(d$x) > 500) .25 else .6)
    }
    p + geom_point(size = if (nrow(d$x) > 500) .6 else 1.2, alpha = .7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey40", linetype = 2, linewidth = .6) +
      geom_line(data = curves, colour = "firebrick", linewidth = 1) +
      facet_wrap(~panel_factor(index), scales = "free_y") +
      labs(x = paste("climatology:", V()$annual), y = "mean exposure per year",
           subtitle = sprintf("%d units, %d-%d; red: quadratic; dashed: loess", nrow(d$x), min(d$years), max(d$years)))
  })
  output$across_table <- renderTable({
    d <- D(); xb <- rowMeans(d$x)
    do.call(rbind, lapply(names(d$idx), function(k) {
      y <- rowMeans(d$idx[[k]])
      q <- lm(y ~ xb + I(xb^2)); q4 <- lm(y ~ poly(xb, 4, raw = TRUE))
      sl <- fits()[[k]]$slope; ok <- is.finite(sl)
      data.frame(index = labels()[[k]], R2_quadratic = summary(q)$r.squared, R2_quartic = summary(q4)$r.squared,
                 p_quartic_vs_quadratic = anova(q, q4)[2, "Pr(>F)"],
                 vertex = -coef(q)[2] / (2 * coef(q)[3]),
                 corr_within_slopes_vs_quadratic = cor(sl[ok], (coef(q)[2] + 2 * coef(q)[3] * xb)[ok]))
    }))
  }, digits = 3)

  # --- Within units ---
  bins <- reactive({
    d <- D()
    do.call(rbind, lapply(names(d$idx), function(k) {
      b <- decile_bins(d$x, d$idx[[k]], groups(), input$within_scale == "std"); b$index <- k; b
    }))
  })
  output$within_plot <- renderPlot({
    ggplot(bins(), aes(z, y, colour = group)) +
      geom_hline(yintercept = 0, linetype = 3) + geom_line() + geom_point(size = 1.3) +
      facet_wrap(~panel_factor(index), scales = "free_y") +
      labs(x = "within-unit deviation of the annual value (unit SD), deciles",
           y = if (input$within_scale == "std") "mean index deviation (unit SD of the index)"
               else "mean index deviation from the unit mean (index units)", colour = NULL) +
      theme(legend.position = "bottom")
  })
  output$arms_table <- renderTable({
    b <- bins()
    out <- do.call(rbind, lapply(split(b, b$index), function(g) cbind(index = labels()[[g$index[1]]], arm_slopes(g))))
    out$ratio_above_below <- out$slope_above / out$slope_below
    out
  }, digits = 3)
  output$within_table <- renderTable({
    do.call(rbind, lapply(names(fits()), function(k) cbind(index = labels()[[k]], within_summary(fits()[[k]], groups()))))
  }, digits = 3)

  # --- Single units ---
  observe({
    u <- H()$units[units_ok(), ]
    ch <- stats::setNames(u$id, u$name)
    keep <- intersect(isolate(input$units), u$id)
    if (!length(keep)) keep <- intersect(c("IND", "USA"), u$id)
    updateSelectizeInput(session, "units", choices = ch, selected = keep, server = TRUE)
  })
  output$units_plot <- renderPlot({
    validate(need(length(input$units) > 0, "Pick one or more units."))
    d <- D(); i <- match(input$units, d$units$id); i <- i[!is.na(i)]
    validate(need(length(i) > 0, "The units picked are not in the current selection."))
    rows <- do.call(rbind, lapply(names(d$idx), function(k) do.call(rbind, lapply(i, function(r) {
      x <- d$x[r, ]; y <- d$idx[[k]][r, ]
      data.frame(index = k, unit = d$units$name[r], year = d$years,
                 z = (x - mean(x)) / stats::sd(x), y = y - mean(y))
    }))))
    ggplot(rows, aes(z, y, colour = unit)) +
      geom_hline(yintercept = 0, linetype = 3) + geom_point(size = 1.5, alpha = .8) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = 2, linewidth = .6) +
      facet_wrap(~panel_factor(index), scales = "free_y") +
      labs(x = "within-unit deviation of the annual value (unit SD), one point per year",
           y = "exposure deviation from the unit mean (index units)", colour = NULL) +
      theme(legend.position = "bottom")
  })
  output$units_table <- renderTable({
    d <- D(); i <- match(input$units, d$units$id); i <- i[!is.na(i)]
    req(length(i) > 0)
    do.call(rbind, lapply(i, function(r) {
      out <- data.frame(unit = d$units$name[r], mean = mean(d$x[r, ]), sd = stats::sd(d$x[r, ]))
      for (k in names(d$idx)) {
        f <- fits()[[k]][r, ]
        out[[paste0(k, "_mean")]] <- mean(d$idx[[k]][r, ])
        out[[paste0(k, "_slope_per_sd")]] <- f$slope * f$sd_x
        out[[paste0(k, "_kappa")]] <- f$kappa
      }
      out
    }))
  }, digits = 2)

  # --- Trailing deviations ---
  trail <- reactive({
    w <- as.integer(input$window)
    need_years <- seq(input$years[1] - w, input$years[2])
    validate(need(min(need_years) >= min(H()$years),
                  sprintf("The %d-year window needs data from %d; start later.", w, min(need_years))))
    u <- valid_units(H(), need_years, input$countries, sources())
    validate(need(length(u) >= 5, "Too few units."))
    d <- D()
    u <- intersect(u, units_ok())
    rc <- trailing_corr(H(), full(), u, yrs(), w, input$combine, d$sds)
    xb <- rowMeans(H()$annual[u, match(yrs(), H()$years), drop = FALSE])
    list(r = rc, xbar = xb, group = climate_groups(xb, breaks(), V()$unit))
  })
  output$trail_plot <- renderPlot({
    t <- trail()
    d <- do.call(rbind, lapply(names(t$r), function(k) data.frame(index = k, x = t$xbar, r = t$r[[k]])))
    ggplot(d[is.finite(d$r), ], aes(x, r)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_point(size = if (nrow(d) > 1500) .6 else 1.2, alpha = .7) +
      facet_wrap(~panel_factor(index)) +
      labs(x = paste("climatology:", V()$annual), y = "within-unit correlation with the deviation of the annual value",
           subtitle = "deviations of both from their trailing mean; units whose index never moves left out")
  })
  output$trail_table <- renderTable({
    t <- trail()
    do.call(rbind, lapply(names(t$r), function(k) {
      r <- t$r[[k]]
      by <- tapply(r, t$group, function(v) sprintf("%.2f [%d/%d]", stats::median(v, na.rm = TRUE), sum(!is.na(v)), length(v)))
      cbind(data.frame(index = labels()[[k]], median = stats::median(r, na.rm = TRUE),
                       p10 = stats::quantile(r, .1, na.rm = TRUE, names = FALSE),
                       p90 = stats::quantile(r, .9, na.rm = TRUE, names = FALSE)),
            as.data.frame(as.list(by), check.names = FALSE))
    }))
  }, digits = 2)

  output$about <- renderUI({
    h <- H()
    tagList(
      p(sprintf("Source: %s (ERA5 daily, Weighted Climate Dataset, population weights 2015). %d units, %d-%d.",
                h$source, nrow(h$units), min(h$years), max(h$years))),
      p("Daily values are the population-weighted mean over the unit: extremes are smoothed, more so the larger the unit (rain more than temperature)."),
      p("Absolute thresholds sit on the histogram's bin edges, so the indices are exact. Unit percentiles are those of the unit's daily values over ",
        sprintf("%d-%d", h$base[1], h$base[2]), " (all days; for rain, dry days included)."),
      p("'Above k' counts days with x >= k, 'below k' days with x < k (ERA5 values are rounded to 0.01, so ties occur). ",
        "Degree days above k are the sum of max(x - k, 0), below k of max(k - x, 0); for rain 'mm beyond' / 'mm short' the same in mm."),
      p("The annual value is the mean of the daily mean T or the total of the daily P; the hot and cold indices can use the daily mean, max or min (built from the WCD daily Tmax / Tmin files, whose units are matched to the daily-mean file by id). Unit SDs are the interannual SD over the selected years, not detrended."),
      p("Units are kept only with every day of every selected year.")
    )
  })
}

app <- shinyApp(ui, server)
# source()d in an interactive session: launch it; loaded by runApp(): return it.
if (interactive() && !is.null(SOURCED_FILE)) runApp(app) else app
