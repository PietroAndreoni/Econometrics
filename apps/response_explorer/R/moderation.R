# The "By characteristic" page: the effect on growth of a year k SD above and
# below the baseline (T and P), against one characteristic of the country,
# the others held at a reference country; and the physical exposure curves of
# the settings with exposure-mediated terms (R/exposure.R).

MODERATION_CHOICES <- c("Long-run temperature (°C)" = "mean_T", "Long-run precipitation (mm a year)" = "mean_P",
                        "GDP per capita (log scale)" = "income", "Agricultural share of value added (%)" = "share",
                        "Capital stock per capita (PWT, log scale)" = "capital")
# How each characteristic is shown: scale for the axis and a log axis.
MODERATION_DISPLAY <- list(mean_T = list(scale = 1, log = FALSE), mean_P = list(scale = 1000, log = FALSE),
                           income = list(scale = 1, log = TRUE), share = list(scale = 100, log = FALSE),
                           capital = list(scale = 1, log = TRUE))
MODERATION_GRID <- 61

# The reference country of setting s: the country `ref` (a country list of
# the Single country page); where it is a preset of the setting's panel, its
# characteristics from that panel (the units of the setting's GDP source).
moderation_base <- function(s, ref) {
  u <- s$units
  if (!is.null(ref$GID_0) && ref$GID_0 %in% u$GID_0) {
    row <- u[u$GID_0 == ref$GID_0, ]
    ref$mean_T <- row$mean_T; ref$sd_T <- row$sd_T; ref$mean_P <- row$mean_P; ref$sd_P <- row$sd_P
    ref$income <- row$income; ref$share <- row$share
    if (is.finite(row$capital)) ref$capital <- row$capital
  }
  ref
}

# The effect (log points, with 95% CI) of a year at +k and -k SD of the
# reference country (T and P), along a grid of characteristic `x` spanning the
# setting's countries, the other characteristics at the reference country;
# `vars`: the variables (each can have a characteristic of its own).
# Along the grid the country is "custom" (no preset code), so exposure takes
# the smooth at its climate; the reference point (point = TRUE) keeps the
# reference country's own values and slopes.
moderation_data <- function(s, ref, x, k, vars = c("TM", "RR")) {
  u <- s$units
  base <- moderation_base(s, ref)
  if (x == "capital" && is.null(base$capital)) return(NULL)
  values <- u[[x]][is.finite(u[[x]])]
  if (!length(values)) return(NULL)
  grid <- if (MODERATION_DISPLAY[[x]]$log) exp(seq(log(min(values)), log(max(values)), length.out = MODERATION_GRID))
          else seq(min(values), max(values), length.out = MODERATION_GRID)
  one <- function(cn, value, point) {
    do.call(rbind, lapply(vars, function(v) {
      sd <- if (v == "TM") base$sd_T else base$sd_P
      shocks <- c(k, -k) * sd / deviation_scale(s$model, cn, v)
      r <- response_curve(s$model, v, shocks, cn, s$cap)
      data.frame(variable = v, sign = c("plus", "minus"), value = value, estimate = r$estimate,
                 low = r$low, high = r$high, point = point)
    }))
  }
  along <- do.call(rbind, lapply(grid, function(g) {
    cn <- base
    cn[[x]] <- g
    cn$GID_0 <- NULL
    one(cn, g, FALSE)
  }))
  rbind(along, one(base, base[[x]], TRUE))
}

# The chart: rows T and P, columns the settings; colour the sign of the
# anomaly; the rug the setting's countries, the marker the reference country.
moderation_plot <- function(d, rugs, x, k, colours) {
  disp <- MODERATION_DISPLAY[[x]]
  d$xv <- d$value * disp$scale
  rugs$xv <- rugs$value * disp$scale
  var_label <- c(TM = "Temperature", RR = "Precipitation")
  d$var <- factor(var_label[d$variable], var_label)
  # The rug in the rows of the variables shown.
  rugs <- merge(rugs, data.frame(var = factor(unique(as.character(d$var)), var_label)), by = NULL)
  sign_label <- c(plus = sprintf("+%g SD (warmer / wetter)", k), minus = sprintf("-%g SD (cooler / drier)", k))
  d$sign <- factor(sign_label[d$sign], sign_label)
  line <- d[!d$point, ]
  pts <- d[d$point, ]
  p <- ggplot(line, aes(xv, 100 * estimate, colour = sign, fill = sign)) +
    geom_hline(yintercept = 0, colour = GRID) +
    geom_ribbon(aes(ymin = 100 * low, ymax = 100 * high), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = rugs, aes(x = xv), inherit.aes = FALSE, sides = "b", colour = INK_MUTED, alpha = 0.4) +
    geom_point(data = pts, size = 2.6, shape = 21, stroke = 1, fill = SURFACE) +
    facet_grid(var ~ setting, scales = "free_y", labeller = labeller(setting = label_wrap_gen(40))) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(x = names(MODERATION_CHOICES)[MODERATION_CHOICES == x],
         y = "Effect on growth of a year at the anomaly (% points)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  if (disp$log) p <- p + scale_x_log10(labels = scales::label_comma())
  p
}

# The exposure curves of the settings with exposure: for each variable and
# channel, each country's slope (points) and the smooth in climate (line), per
# 1 C of annual mean T or 1 m a year of annual P; the reference country's
# climate marked.
exposure_curve_data <- function(sp) {
  out <- list()
  for (s in sp) {
    m <- s$model
    for (v in c("TM", "RR")) {
      if (!exposure_on(m$exposure, v)) next
      tab <- m$exposure_table[[v]]
      spec <- EXPOSURE_VARS[[v]]
      for (k in c("hi", "lo")) {
        col <- spec$cols[[k]]
        if (!col %in% exposure_columns(m$exposure, v)) next
        label <- sprintf("%s: %s, %s", if (v == "TM") "T" else "P", spec$sides[[k]],
                         exposure_index_label(m$exposure[[v]][[k]], v, k))
        out[[length(out) + 1]] <- rbind(
          data.frame(setting = s$label, variable = v, channel = label, clim = tab$units$clim,
                     slope = tab$units[[paste0("s_", k)]], kind = "country"),
          data.frame(setting = s$label, variable = v, channel = label, clim = tab$grid$clim,
                     slope = tab$grid[[paste0("s_", k)]], kind = "smooth"))
      }
    }
  }
  if (length(out)) do.call(rbind, out)
}
exposure_curve_plot <- function(d, marks) {
  pts <- d[d$kind == "country", ]
  line <- d[d$kind == "smooth", ]
  p <- ggplot(pts, aes(clim, slope)) +
    geom_hline(yintercept = 0, colour = GRID) +
    geom_point(colour = INK_MUTED, alpha = 0.5, size = 1.2) +
    geom_line(data = line, colour = SETTING_COLOURS[1], linewidth = 1) +
    facet_wrap(~ setting + channel, scales = "free", labeller = labeller(setting = label_wrap_gen(40))) +
    labs(x = "Long-run climate of the country (annual mean T in C, or annual P in mm)",
         y = "Exposure slope (index per C of T, or per m of P)") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  if (!is.null(marks) && nrow(marks)) {
    p <- p + geom_vline(data = marks, aes(xintercept = clim), colour = INK, linetype = "dashed")
  }
  p
}

# The heatmap: the effect (log points) of a year at sign * k SD of variable v
# of the reference country, over a grid of two characteristics (x, y), the
# others at the reference country; all cells at once (point estimates).
HEAT_GRID <- 50
moderation_grid <- function(u, x, n) {
  z <- u[[x]][is.finite(u[[x]])]
  if (!length(z)) return(NULL)
  if (MODERATION_DISPLAY[[x]]$log) exp(seq(log(min(z)), log(max(z)), length.out = n)) else seq(min(z), max(z), length.out = n)
}
heat_data <- function(s, ref, v, sign, k, x, y) {
  base <- moderation_base(s, ref)
  if ("capital" %in% c(x, y) && is.null(base$capital)) return(NULL)
  gx <- moderation_grid(s$units, x, HEAT_GRID)
  gy <- moderation_grid(s$units, y, HEAT_GRID)
  if (is.null(gx) || is.null(gy)) return(NULL)
  cells <- expand.grid(xv = gx, yv = gy)
  # The country of every cell (no preset code: exposure is the smooth).
  cn <- base
  cn$GID_0 <- NULL
  cn[[x]] <- cells$xv
  cn[[y]] <- cells$yv
  sd <- if (v == "TM") base$sd_T else base$sd_P
  scale <- deviation_scale(s$model, base, v)
  xs <- rep(sign * k * sd / scale, nrow(cells))
  cur <- list(TM = 0, RR = 0)
  cur[[v]] <- xs * scale
  mods <- modifier_columns(s$model, cn, list(TM = 0, RR = 0), cur = cur)
  inputs <- if (isTRUE(s$model$both)) list(x = xs, d = xs) else xs
  off <- rep_len(deviation_offset(s$model, cn, v), nrow(cells))
  beta <- point_beta(s$model)
  effect <- drop(growth_matrix(beta, s$model, v, plus_offset(inputs, off), mods, s$cap))
  if (any(off != 0)) effect <- effect - drop(growth_matrix(beta, s$model, v, offset_inputs(inputs, off), mods, s$cap))
  data.frame(setting = s$label, xv = cells$xv, yv = cells$yv, effect = effect)
}
heat_plot <- function(d, countries, reference, x, y, title) {
  dx <- MODERATION_DISPLAY[[x]]; dy <- MODERATION_DISPLAY[[y]]
  for (z in c("d", "countries", "reference")) {
    t <- get(z)
    t$xs <- t$xv * dx$scale
    t$ys <- t$yv * dy$scale
    assign(z, t)
  }
  lim <- max(abs(100 * d$effect), na.rm = TRUE)
  p <- ggplot(d, aes(xs, ys)) +
    geom_raster(aes(fill = 100 * effect), interpolate = TRUE) +
    geom_point(data = countries, colour = INK, alpha = 0.35, size = 0.9) +
    geom_point(data = reference, shape = 21, size = 3.2, stroke = 1.2, colour = INK, fill = SURFACE) +
    facet_wrap(~ setting, labeller = labeller(setting = label_wrap_gen(40))) +
    scale_fill_gradient2(low = "#c0392b", mid = "#f7f7f5", high = "#2a78d6", midpoint = 0, limits = c(-lim, lim),
                         name = "Effect on growth\n(% points)") +
    labs(x = names(MODERATION_CHOICES)[MODERATION_CHOICES == x], y = names(MODERATION_CHOICES)[MODERATION_CHOICES == y],
         title = title) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  p <- p + (if (dx$log) scale_x_log10(labels = scales::label_comma(), expand = c(0, 0)) else scale_x_continuous(expand = c(0, 0)))
  p + (if (dy$log) scale_y_log10(labels = scales::label_comma(), expand = c(0, 0)) else scale_y_continuous(expand = c(0, 0)))
}
