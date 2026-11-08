# Exposure-mediated deviation terms (bhm_exposure_mediated.R), for T and / or
# P: every deviation term of the variable times the country's exposure slope
# of its upper index and, as a copy, of its lower index, with no plain term:
#
#   g = sum_a [(b_a^hi + ...) s^hi + (b_a^lo + ...) s^lo] arm_a(x) / 100 + ...
#
#   T: hot (column "exh") and cold ("exc"): e.g. degree days of Tmax above
#      29 C, of Tmin below 0 C; slopes per C of annual mean T;
#   P: wet ("exw") and dry ("exd"): e.g. mm of daily P beyond the unit's 99th
#      percentile, days below 1 mm; slopes per m of annual P (the model's
#      units), from per mm.
# so a deviation acts through the change of the indices it brings. The slopes
# come from the exposure explorer's histograms (GADM0, results/exposure_explorer,
# apps/exposure_explorer/build_exposure.R), ERA5 daily, population weights 2015:
#   s_i = within-country OLS slope of the deviation of the index from its
#         20-year trailing mean (t-20..t-1) on that of the annual value (mean
#         T, total P), 1970-2023 (climate data only);
#   "country": each country's own slope;
#   "smooth":  the loess (span 0.5) of the slopes on the countries' 1970-2023
#              climatologies, evaluated at the country's long-run climate (+
#              any shift in the paths and projections), with the sign of the
#              index (upper >= 0, lower <= 0); flat beyond the range of the
#              climatologies.
# The terms are in the index per 100, so the deviations must be in absolute
# units (C, m): model_options() sets them.
#
# Options: list(TM = list(mode, channels, hi, lo), RR = ...), each index
# list(src, type ("abs" or "pct"), value (threshold, or the percentile's
# probability), measure ("dd": degree days / mm beyond; "days")). The older
# form list(mode, channels, hot, cold) is T's.

EXPOSURE_MODES <- c("None (the standard terms)" = "none", "Country exposure slopes" = "country",
                    "Smooth exposure (slope as a function of the climate)" = "smooth")
EXPOSURE_VARS <- list(
  TM = list(base = "tmp", cols = c(hi = "exh", lo = "exc"), sides = c(hi = "hot", lo = "cold"),
            groups = c(hi = "hot", lo = "cold"),
            # The index per unit of the model's deviation (C) from per unit of
            # the histogram's annual value (C); the long-run climate in the
            # histogram's units from the model's (C).
            per_unit = 1, clim_scale = 1, unit = "C", annual = "annual mean T",
            channels = c("Hot and cold" = "both", "Hot only" = "hi", "Cold only" = "lo"),
            default = list(hi = list(src = "tmpmax", type = "abs", value = 29, measure = "dd"),
                           lo = list(src = "tmpmin", type = "abs", value = 0, measure = "dd"))),
  RR = list(base = "pre", cols = c(hi = "exw", lo = "exd"), sides = c(hi = "wet", lo = "dry"),
            groups = c(hi = "wet", lo = NA),
            per_unit = 1000, clim_scale = 1000, unit = "mm", annual = "annual P",
            channels = c("Wet and dry" = "both", "Wet only" = "hi", "Dry only" = "lo"),
            default = list(hi = list(src = "pre", type = "pct", value = 0.99, measure = "dd"),
                           lo = list(src = "pre", type = "abs", value = 1, measure = "days"))))
EXPOSURE_CHANNELS <- EXPOSURE_VARS$TM$channels
EXPOSURE_SOURCES <- c("Daily max" = "tmpmax", "Daily mean" = "tmp", "Daily min" = "tmpmin")
EXPOSURE_SOURCE_SHORT <- c(tmpmax = "Tmax", tmp = "daily mean T", tmpmin = "Tmin", pre = "daily P")
EXPOSURE_MEASURES <- c("Degree days" = "dd", "Days" = "days")
EXPOSURE_P_MEASURES <- c("mm beyond" = "dd", "Days" = "days")
EXPOSURE_COLUMN_LABELS <- c(exh = "hot", exc = "cold", exw = "wet", exd = "dry")
EXPOSURE_WINDOW <- 20
EXPOSURE_YEARS <- 1970:2023
EXPOSURE_SPAN <- 0.5
# The exposure explorer's functions (index_matrix, aligned, trailing_dev), in
# an environment of their own; its histograms.
EXPOSURE_APP_DIR <- if (exists("APP_DIR")) APP_DIR else file.path(getwd(), "apps", "response_explorer")
EXPOSURE_DATA_DIR <- normalizePath(file.path(EXPOSURE_APP_DIR, "..", "..", "results", "exposure_explorer"),
                                   mustWork = FALSE)
XP <- new.env()
local({
  f <- file.path(EXPOSURE_APP_DIR, "..", "exposure_explorer", "R", "exposure.R")
  if (file.exists(f)) sys.source(f, envir = XP)
})
exposure_available <- function(v = "TM") {
  files <- if (v == "TM") c("tmp", "tmpmax", "tmpmin") else "pre"
  exists("index_matrix", envir = XP, inherits = FALSE) &&
    all(file.exists(file.path(EXPOSURE_DATA_DIR, paste0("gadm0_", files, ".rds"))))
}

# The options of the exposure terms, completed per variable (mode "none": no
# exposure of that variable).
exposure_options <- function(e) {
  e <- e %||% list()
  # The older form: T's options at the top.
  if (!is.null(e$mode)) {
    e <- list(TM = list(mode = e$mode, channels = c(both = "both", hot = "hi", cold = "lo", hi = "hi", lo = "lo")[[e$channels %||% "both"]],
                        hi = e$hot %||% e$hi, lo = e$cold %||% e$lo))
  }
  out <- list()
  for (v in c("TM", "RR")) {
    ev <- e[[v]] %||% list()
    spec <- EXPOSURE_VARS[[v]]
    if (identical(ev$mode %||% "none", "none")) {
      out[[v]] <- list(mode = "none")
      next
    }
    o <- list(mode = ev$mode, channels = ev$channels %||% "both")
    for (k in c("hi", "lo")) {
      x <- spec$default[[k]]
      for (z in c("src", "type", "value", "measure")) if (!is.null(ev[[k]][[z]])) x[[z]] <- ev[[k]][[z]]
      if (v == "RR") x$src <- "pre"
      if (is.na(spec$groups[[k]])) x$type <- "abs"            # no percentiles of the dry index
      x$value <- as.numeric(x$value)
      o[[k]] <- x
    }
    out[[v]] <- o
  }
  out
}
# The exposure columns of variable v (or of both), and their channels.
exposure_columns <- function(e, v = c("TM", "RR")) {
  unlist(lapply(v, function(z) {
    ez <- e[[z]]
    if (is.null(ez) || identical(ez$mode %||% "none", "none")) return(character(0))
    cols <- EXPOSURE_VARS[[z]]$cols
    unname(c(if (ez$channels %in% c("both", "hi")) cols[["hi"]], if (ez$channels %in% c("both", "lo")) cols[["lo"]]))
  }))
}
# An index in words: "Tmax DD > 29 C", "daily P mm > 99th pct".
exposure_index_label <- function(x, v, k) {
  side <- if (k == "hi") ">" else "<"
  thr <- if (x$type == "pct") sprintf("the unit's %gth pct", 100 * x$value) else paste0(format(x$value), " ", EXPOSURE_VARS[[v]]$unit)
  what <- if (x$measure == "days") "days" else if (v == "TM") "DD" else "mm"
  sprintf("%s %s %s %s", EXPOSURE_SOURCE_SHORT[[x$src]], what, side, thr)
}
exposure_label <- function(e) {
  parts <- unlist(lapply(c("TM", "RR"), function(v) {
    cols <- exposure_columns(e, v)
    if (!length(cols)) return(NULL)
    ks <- names(EXPOSURE_VARS[[v]]$cols)[EXPOSURE_VARS[[v]]$cols %in% cols]
    paste0(if (e[[v]]$mode == "smooth") "smooth " else "", if (v == "TM") "T" else "P", " exposure: ",
           paste(vapply(ks, function(k) exposure_index_label(e[[v]][[k]], v, k), ""), collapse = " / "))
  }))
  if (length(parts)) paste(parts, collapse = " · ")
}
exposure_key <- function(e) {
  parts <- unlist(lapply(c("TM", "RR"), function(v) {
    ev <- e[[v]]
    if (!length(exposure_columns(e, v))) return(NULL)
    idx <- function(x) paste(x$src, x$type, x$measure, x$value)
    sprintf("exposure %s %s %s: hi %s, lo %s", v, ev$mode, ev$channels, idx(ev$hi), idx(ev$lo))
  }))
  if (length(parts)) paste(parts, collapse = "; ")
}

# The slopes of a variable's pair of indices (hi, lo) for every country, and
# their smooth in climate: list(units = data.frame(GID_0, clim, s_hi, s_lo),
# grid = data.frame(clim, s_hi, s_lo)), per unit of the model's deviation;
# cached by index.
exposure_cache <- new.env()
exposure_hist <- function(var) {
  key <- paste0("hist_", var)
  if (is.null(exposure_cache[[key]])) exposure_cache[[key]] <- readRDS(file.path(EXPOSURE_DATA_DIR, paste0("gadm0_", var, ".rds")))
  exposure_cache[[key]]
}
exposure_table <- function(ev, v) {
  spec <- EXPOSURE_VARS[[v]]
  idx <- function(x) paste(x$src, x$type, x$measure, x$value)
  key <- paste(v, idx(ev$hi), idx(ev$lo))
  if (!is.null(exposure_cache[[key]])) return(exposure_cache[[key]])
  if (!exposure_available(v)) stop("No exposure histograms in ", EXPOSURE_DATA_DIR,
                                   ": run Rscript apps/exposure_explorer/build_exposure.R gadm0")
  H <- exposure_hist(spec$base)
  index_of <- function(x, k) {
    Hs <- exposure_hist(x$src)
    XP$aligned(H, Hs, XP$index_matrix(Hs, list(side = if (k == "hi") "above" else "below", measure = x$measure,
                                               type = x$type, value = x$value, group = spec$groups[[k]])))
  }
  hi <- index_of(ev$hi, "hi"); lo <- index_of(ev$lo, "lo")
  need <- match(seq(min(EXPOSURE_YEARS) - EXPOSURE_WINDOW, max(EXPOSURE_YEARS)), H$years)
  complete <- function(M) rowSums(is.na(M[, need, drop = FALSE])) == 0
  ok <- which(rowSums(H$ndays[, need, drop = FALSE] >= 365) == length(need) & complete(hi) & complete(lo))
  dx <- XP$trailing_dev(H$annual[ok, , drop = FALSE], H$years, EXPOSURE_YEARS, EXPOSURE_WINDOW)
  slope <- function(M) {
    dy <- XP$trailing_dev(M[ok, , drop = FALSE], H$years, EXPOSURE_YEARS, EXPOSURE_WINDOW)
    xc <- dx - rowMeans(dx); yc <- dy - rowMeans(dy)
    s <- rowSums(xc * yc) / rowSums(xc^2)
    spec$per_unit * ifelse(rowSums(yc^2) > 0, s, 0)
  }
  units <- data.frame(GID_0 = H$units$id[ok], clim = rowMeans(H$annual[ok, match(EXPOSURE_YEARS, H$years), drop = FALSE]),
                      s_hi = slope(hi), s_lo = slope(lo), stringsAsFactors = FALSE)
  step <- if (v == "TM") 0.1 else 10
  grid <- data.frame(clim = seq(floor(min(units$clim) / step) * step, ceiling(max(units$clim) / step) * step, by = step))
  sm <- function(s, sign) {
    lo <- stats::loess(s ~ clim, data = data.frame(s = s, clim = units$clim), span = EXPOSURE_SPAN)
    p <- stats::predict(lo, data.frame(clim = pmin(pmax(grid$clim, min(units$clim)), max(units$clim))))
    if (sign > 0) pmax(p, 0) else pmin(p, 0)
  }
  grid$s_hi <- sm(units$s_hi, 1); grid$s_lo <- sm(units$s_lo, -1)
  out <- list(units = units, grid = grid)
  exposure_cache[[key]] <- out
  out
}

# The exposure columns of variable v (exh / exc, or exw / exd; per 100 of the
# index per unit of the deviation) of countries `GID_0` (one or a vector) with
# long-run climate `clim` (the model's units: C, m; vector or scalar): the
# country's slope ("country", where it has one) or the smooth at the climate
# ("smooth"; also a country without a slope).
exposure_values <- function(ev, v, table, GID_0, clim) {
  spec <- EXPOSURE_VARS[[v]]
  n <- max(length(GID_0), length(clim))
  GID_0 <- rep_len(as.character(GID_0 %||% NA), n)
  clim <- rep_len(clim * spec$clim_scale, n)
  at <- match(GID_0, table$units$GID_0)
  active <- exposure_columns(stats::setNames(list(ev), v), v)
  out <- sapply(c("hi", "lo"), function(k) {
    if (!spec$cols[[k]] %in% active) return(rep(0, n))
    s <- stats::approx(table$grid$clim, table$grid[[paste0("s_", k)]], xout = clim, rule = 2)$y
    if (ev$mode == "country") s <- ifelse(is.na(at), s, table$units[[paste0("s_", k)]][at])
    s / 100
  })
  out <- matrix(out, nrow = n)
  colnames(out) <- unname(spec$cols)
  out
}
