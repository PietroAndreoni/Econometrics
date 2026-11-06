# Builds the data of the exposure explorer (apps/exposure_explorer/app.R):
# per unit-year histograms of the ERA5 daily mean, max and min temperature and
# daily precipitation (Weighted Climate Dataset, population weights 2015, the
# only daily ones), for countries (gadm0) and regions (gadm1).
#
# Each histogram bin j = [lo_j, lo_{j+1}) stores the number of days and the
# sum of (x - lo_j) over them, so that with a threshold k on a bin edge
#   degree days above k = sum_{lo_j >= k} n_j (lo_j - k) + e_j
#   degree days below k = sum_{lo_j <  k} n_j (k - lo_j) - e_j
#   days above / below k = sum of n_j over the same bins
# are exact (the first bin is open below and the last open above, with their
# lo as reference, so the formulas hold for them too). The annual mean T and
# the annual total P are stored apart. Thresholds at each unit's own daily
# percentiles (base 1970-2023, all days) are stored as exact counts and
# degree days too.
#
#   Rscript apps/exposure_explorer/build_exposure.R                      # everything
#   Rscript apps/exposure_explorer/build_exposure.R gadm0                # one resolution
#   Rscript apps/exposure_explorer/build_exposure.R gadm1 tmpmax tmpmin  # some variables
#
# Output: results/exposure_explorer/<resolution>_<variable>.rds

source(file.path("functions", "climate_download.R"))
out_dir <- file.path("results", "exposure_explorer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
YEARS <- 1940:2023                      # 1940 = the first ERA5 year, for 30-year trailing means
BASE <- 1970:2023                       # base period of the unit percentiles
CHUNK <- 250                            # units read at a time
T_PCT <- list(hot = c(0.90, 0.95, 0.99), cold = c(0.01, 0.05, 0.10))
SPEC <- list(
  tmp = list(lo = c(-30:40), scale = 100, pct = T_PCT),
  # Daily max / min: wider edges (Tmax reaches ~50 C, Tmin ~-50 C).
  tmpmax = list(lo = c(-35:50), scale = 100, pct = T_PCT),
  tmpmin = list(lo = c(-45:40), scale = 100, pct = T_PCT),
  pre = list(lo = c(0, 0.1, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150),
             scale = 100, pct = list(wet = c(0.90, 0.95, 0.99, 0.999)))
)
# Arguments: resolutions (gadm0, gadm1) and / or variables (tmp, tmpmax,
# tmpmin, pre); none = all.
args <- commandArgs(trailingOnly = TRUE)
RESOLUTIONS <- intersect(args, c("gadm0", "gadm1")); if (!length(RESOLUTIONS)) RESOLUTIONS <- c("gadm0", "gadm1")
VARIABLES <- intersect(args, names(SPEC)); if (!length(VARIABLES)) VARIABLES <- names(SPEC)

# Unit names: GADM0 from the WCD country list, GADM1 from DOSE where it has them.
countries <- wcd_countries("gadm0")
country_name <- setNames(countries$COUNTRY, countries$GID_0)
dose <- utils::read.csv(file.path("data", "DOSE_V2.14.csv"), stringsAsFactors = FALSE)[c("GID_1", "region")]
region_name <- setNames(dose$region, dose$GID_1)[!duplicated(dose$GID_1)]

build <- function(res, var) {
  sp <- SPEC[[var]]
  path <- wcd_download(var, "ERA5", res, "pop", 2015, "daily", verbose = TRUE)
  pq <- arrow::open_dataset(path)
  ids <- setdiff(names(pq), c("Date", "NA", "?", ""))
  dates <- as.Date(sub("^X", "", as.data.frame(arrow::read_parquet(path, col_select = "Date"))$Date), "%Y%m%d")
  yr <- as.integer(format(dates, "%Y"))
  keep <- yr %in% YEARS
  yi <- match(yr[keep], YEARS)
  nu <- length(ids); ny <- length(YEARS); nb <- length(sp$lo)
  cnt <- array(0L, c(nu, ny, nb)); exc <- array(0L, c(nu, ny, nb))
  ndays <- matrix(0L, nu, ny); annual <- matrix(NA_real_, nu, ny)
  pct <- lapply(sp$pct, function(p) list(probs = p, thr = matrix(NA_real_, nu, length(p)),
                                         n = array(0L, c(nu, ny, length(p))),
                                         dd = array(NA_real_, c(nu, ny, length(p)))))
  base <- yr[keep] %in% BASE
  for (s in seq(1, nu, by = CHUNK)) {
    u <- s:min(nu, s + CHUNK - 1)
    m <- as.matrix(as.data.frame(arrow::read_parquet(path, col_select = ids[u])))[keep, , drop = FALSE]
    storage.mode(m) <- "double"
    if (var == "pre") m <- pmax(m, 0)
    ok <- !is.na(m)
    ndays[u, ] <- t(rowsum(ok * 1L, yi))
    annual[u, ] <- t(if (var != "pre") rowsum(m, yi) / rowsum(ok * 1, yi) else rowsum(m, yi))
    for (k in seq_along(u)) {
      x <- m[, k]; o <- ok[, k]
      if (!any(o)) next                 # a unit without data stays empty (ndays 0)
      j <- pmax(findInterval(x[o], sp$lo), 1L)
      key <- (yi[o] - 1L) * nb + j
      cnt[u[k], , ] <- matrix(tabulate(key, ny * nb), ny, nb, byrow = TRUE)
      e <- rowsum(x[o] - sp$lo[j], key)
      ee <- numeric(ny * nb); ee[as.integer(rownames(e))] <- e[, 1]
      exc[u[k], , ] <- matrix(as.integer(round(ee * sp$scale)), ny, nb, byrow = TRUE)
      # Own-percentile thresholds and their exact annual exposure.
      xb <- x[base & o]
      for (g in names(pct)) {
        q <- stats::quantile(xb, pct[[g]]$probs, names = FALSE, type = 7)
        pct[[g]]$thr[u[k], ] <- q
        for (h in seq_along(q)) {
          if (g == "cold") { hit <- x[o] < q[h]; amt <- pmax(q[h] - x[o], 0) }
          else { hit <- x[o] > q[h]; amt <- pmax(x[o] - q[h], 0) }
          pct[[g]]$n[u[k], , h] <- tabulate(yi[o][hit], ny)
          rs <- rowsum(amt, yi[o]); dd <- numeric(ny); dd[as.integer(rownames(rs))] <- rs[, 1]
          pct[[g]]$dd[u[k], , h] <- dd
        }
      }
    }
    cat(sprintf("  %s %s: units %d-%d of %d\n", res, var, min(u), max(u), nu))
  }
  gid0 <- sub("^([A-Z]{3}).*$", "\\1", ids)
  gid1 <- if (res == "gadm1") sub("^([A-Z]{3})_", "\\1.", ids) else ids
  name <- if (res == "gadm1") {
    rn <- region_name[gid1]
    ifelse(is.na(rn), gid1, paste0(rn, " (", gid1, ")"))
  } else {
    cn <- country_name[ids]
    ifelse(is.na(cn), ids, paste0(cn, " (", ids, ")"))
  }
  list(resolution = res, variable = var, years = YEARS, base = range(BASE),
       units = data.frame(id = ids, GID_0 = gid0, GID_1 = gid1, name = unname(name),
                          country = unname(ifelse(is.na(country_name[gid0]), gid0, country_name[gid0])),
                          stringsAsFactors = FALSE),
       lo = sp$lo, scale = sp$scale, cnt = cnt, exc = exc, ndays = ndays, annual = annual, pct = pct,
       source = basename(path), built = Sys.time())
}

for (res in RESOLUTIONS) for (var in VARIABLES) {
  f <- file.path(out_dir, paste0(res, "_", var, ".rds"))
  cat("Building", f, "\n")
  t0 <- Sys.time()
  b <- build(res, var)
  saveRDS(b, f, compress = "xz")
  cat(sprintf("  saved (%.0f MB on disk, %.1f min)\n", file.size(f) / 1e6,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
