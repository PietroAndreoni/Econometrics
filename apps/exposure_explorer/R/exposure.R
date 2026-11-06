# Computations of the exposure explorer, on the histograms of build_exposure.R.
#
# H (one resolution x variable): units, years, lo (bin lower edges), scale,
# cnt / exc [unit, year, bin] (days; sum of x - lo in 1/scale units),
# ndays, annual [unit, year] (mean T / total P), pct (own-percentile
# thresholds: probs, thr [unit, prob], n and dd [unit, year, prob]).

# An index [unit, year] from a definition:
#   side    "above" or "below" the threshold
#   measure "dd" (degree days / mm beyond the threshold) or "days"
#   type    "abs" (threshold = a bin edge, value) or "pct" (the unit's own
#           percentile; group = "hot", "cold" or "wet", value = the prob)
index_matrix <- function(H, def) {
  if (def$type == "pct") {
    p <- H$pct[[def$group]]
    h <- match(as.numeric(def$value), p$probs)
    stopifnot(!is.na(h))
    return(if (def$measure == "days") p$n[, , h] * 1 else p$dd[, , h])
  }
  k <- as.numeric(def$value)
  lo <- H$lo
  stopifnot(k %in% lo, k > lo[1] || H$variable == "pre")
  sel <- if (def$side == "above") which(lo >= k) else which(lo < k)
  out <- matrix(0, dim(H$cnt)[1], dim(H$cnt)[2])
  for (j in sel) {
    n <- H$cnt[, , j]
    if (def$measure == "days") {
      out <- out + n
    } else if (def$side == "above") {
      out <- out + n * (lo[j] - k) + H$exc[, , j] / H$scale
    } else {
      out <- out + n * (k - lo[j]) - H$exc[, , j] / H$scale
    }
  }
  out
}

# An index from another file of the same resolution (daily max / min) on the
# units and years of H (the file of the annual value): rows matched by unit
# id, NA where that file lacks the unit; with its day counts likewise.
aligned <- function(H, Hs, M) {
  stopifnot(identical(H$years, Hs$years))
  M[match(H$units$id, Hs$units$id), , drop = FALSE]
}
index_aligned <- function(H, Hs, def) aligned(H, Hs, index_matrix(Hs, def))

# Units with every day of every year in `need` (years), optionally in some
# countries; `sources` = further files (daily max / min) whose days must be
# complete too.
valid_units <- function(H, need, countries = NULL, sources = list()) {
  yi <- match(need, H$years)
  ok <- rowSums(H$ndays[, yi, drop = FALSE] >= 365) == length(yi)
  for (Hs in sources) {
    nd <- aligned(H, Hs, Hs$ndays)[, yi, drop = FALSE]
    ok <- ok & rowSums(!is.na(nd) & nd >= 365) == length(yi)
  }
  if (length(countries)) ok <- ok & H$units$GID_0 %in% countries
  which(ok)
}

# The three indices of the view (hi, lo, sum) and the annual value, on the
# rows u and years yrs; full = list(hi, lo) of index matrices [unit, year] on
# the units and years of H. Sum: "raw" = hi + lo, "sd" = each over its
# cross-unit SD of the climatology (over these units and years).
view_data <- function(H, full, u, yrs, combine = "raw") {
  yi <- match(yrs, H$years)
  hi <- full$hi[u, yi, drop = FALSE]
  lo <- full$lo[u, yi, drop = FALSE]
  sds <- c(hi = stats::sd(rowMeans(hi)), lo = stats::sd(rowMeans(lo)))
  sum <- if (combine == "sd") hi / sds[["hi"]] + lo / sds[["lo"]] else hi + lo
  list(x = H$annual[u, yi, drop = FALSE], idx = list(hi = hi, lo = lo, sum = sum), sds = sds,
       units = H$units[u, ], years = yrs)
}

# Climate groups of the units from breaks of the climatology.
climate_groups <- function(xbar, breaks, unit_label) {
  b <- sort(unique(breaks))
  labs <- if (!length(b)) "all" else {
    edges <- c(-Inf, b, Inf)
    vapply(seq_len(length(edges) - 1), function(i) {
      if (i == 1) sprintf("< %s %s", format(b[1]), unit_label)
      else if (i == length(edges) - 1) sprintf("> %s %s", format(b[length(b)]), unit_label)
      else sprintf("%s-%s %s", format(edges[i]), format(edges[i + 1]), unit_label)
    }, "")
  }
  if (!length(b)) return(factor(rep("all", length(xbar))))
  cut(xbar, c(-Inf, b, Inf), labels = labs)
}

# Within-unit fits of y on x (rows = units): linear slope, quadratic
# curvature c (y ~ x + x^2, x centred), its p-value, R2 gain, kappa =
# |4 c sd(x) / b| (slope change over +-1 SD relative to the slope), the
# linear slope's p-value.
within_fits <- function(x, y) {
  n <- ncol(x)
  xc <- x - rowMeans(x); yc <- y - rowMeans(y)
  q <- xc^2; q <- q - rowMeans(q)
  s11 <- rowSums(xc^2); s12 <- rowSums(xc * q); s22 <- rowSums(q^2)
  y1 <- rowSums(xc * yc); y2 <- rowSums(q * yc); yy <- rowSums(yc^2)
  b_lin <- y1 / s11
  rss_lin <- yy - b_lin * y1
  det <- s11 * s22 - s12^2
  b <- (s22 * y1 - s12 * y2) / det
  c2 <- (s11 * y2 - s12 * y1) / det
  rss_q <- yy - b * y1 - c2 * y2
  se_c <- sqrt(pmax(rss_q, 0) / (n - 3) * s11 / det)
  se_b_lin <- sqrt(pmax(rss_lin, 0) / (n - 2) / s11)
  varies <- yy > 1e-12
  data.frame(slope = b_lin, p_slope = 2 * stats::pt(-abs(b_lin / se_b_lin), n - 2),
             curv = c2, p_curv = 2 * stats::pt(-abs(c2 / se_c), n - 3),
             r2_lin = 1 - rss_lin / yy, r2_gain = (rss_lin - rss_q) / yy,
             kappa = abs(4 * c2 * sqrt(s11 / (n - 1)) / b),
             sd_x = sqrt(s11 / (n - 1)), varies = varies)
}

# Summary of within_fits by group (kappa only where the slope is significant).
within_summary <- function(f, group) {
  do.call(rbind, lapply(c("all", levels(group)), function(g) {
    s <- if (g == "all") rep(TRUE, nrow(f)) else group == g
    ff <- f[s & f$varies, ]
    kk <- ff$kappa[ff$p_slope < 0.05]
    data.frame(group = g, units = sum(s), varying = nrow(ff),
               median_slope = stats::median(ff$slope),
               r2_lin = stats::median(ff$r2_lin), r2_gain = stats::median(ff$r2_gain),
               share_curv_p05 = mean(ff$p_curv < 0.05), convex = mean(ff$curv > 0),
               kappa_median = if (length(kk)) stats::median(kk) else NA_real_,
               kappa_p75 = if (length(kk)) stats::quantile(kk, .75, names = FALSE) else NA_real_)
  }))
}

# Deciles of the within-unit deviation (unit SD) against the mean index
# deviation, absolute (index units) or standardized (unit SD of the index),
# by group; plus the slopes below / above the mean.
decile_bins <- function(x, y, group, standardize = FALSE) {
  zx <- (x - rowMeans(x)) / apply(x, 1, stats::sd)
  dy <- y - rowMeans(y)
  if (standardize) {
    s <- apply(y, 1, stats::sd)
    keep <- s > 0
    dy <- dy[keep, , drop = FALSE] / s[keep]; zx <- zx[keep, , drop = FALSE]; group <- group[keep]
  }
  d <- data.frame(z = as.vector(zx), y = as.vector(dy), group = rep(group, times = ncol(x)))
  d <- d[is.finite(d$z) & is.finite(d$y), ]
  d$bin <- cut(d$z, unique(stats::quantile(d$z, seq(0, 1, .1))), include.lowest = TRUE)
  stats::aggregate(cbind(z, y) ~ bin + group, d, mean)
}
arm_slopes <- function(b) {
  do.call(rbind, lapply(split(b, b$group, drop = TRUE), function(g) {
    lo <- g[g$z < 0, ]; hi <- g[g$z > 0, ]
    sl <- function(d) if (nrow(d) >= 2) unname(stats::coef(stats::lm(y ~ z, d))[2]) else NA_real_
    data.frame(group = g$group[1], slope_below = sl(lo), slope_above = sl(hi))
  }))
}

# Within-unit correlation between the deviation of the index and that of
# the annual value from their trailing mean over the `window` years before.
trailing_dev <- function(M, years, yrs, window) {
  cs <- cbind(0, t(apply(M, 1, cumsum)))
  ti <- match(yrs, years)
  prev <- (cs[, ti, drop = FALSE] - cs[, ti - window, drop = FALSE]) / window   # years t-window..t-1
  M[, ti, drop = FALSE] - prev
}
trailing_corr <- function(H, full, u, yrs, window, combine, sds) {
  dx <- trailing_dev(H$annual[u, , drop = FALSE], H$years, yrs, window)
  hi <- full$hi[u, , drop = FALSE]; lo <- full$lo[u, , drop = FALSE]
  idx <- list(hi = hi, lo = lo, sum = if (combine == "sd") hi / sds[["hi"]] + lo / sds[["lo"]] else hi + lo)
  lapply(idx, function(M) {
    dy <- trailing_dev(M, H$years, yrs, window)
    vapply(seq_len(nrow(dy)), function(i) {
      if (stats::sd(dy[i, ]) == 0) NA_real_ else stats::cor(dx[i, ], dy[i, ])
    }, 0)
  })
}
