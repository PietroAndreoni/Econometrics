# =============================================================================
# download_weighted_climate_data.R
#
# Programmatic (R) access to the Weighted Climate Dataset (WCD) served by
#   https://weightedclimatedata.streamlit.app/Download_Data
#
# Gortan, M., Testa, L., Fagiolo, G., Lamperti, F. (2024)
# "A unified dataset for pre-processed climate indicators weighted by gridded
#  economic activity", Scientific Data. doi:10.1038/s41597-024-03304-1
#
# HOW IT WORKS
# ------------
# The Streamlit dashboard itself reads the parquet files from a *private* Google
# Cloud Storage bucket (credentials live in the app's `st.secrets`), so the app
# URL cannot be queried directly. The very same parquet files are published by
# the authors in public, authentication-free Box folders -- one per geographical
# resolution -- advertised in the dashboard's "Interested in bulk downloads?"
# panel. This script talks to those folders: it lists them, resolves the file
# implied by your selection, downloads it (with an on-disk cache) and then
# reproduces the dashboard's post-processing (year/unit subsetting, yearly
# aggregation, threshold indicators, ERA5 precipitation rescaling, long format).
#
# Files are named exactly as the dashboard names them:
#   <geo_resolution>_<source>_<variable>_<weight>_<weight_year>_<frequency>.parquet
# e.g. gadm0_cru_tmp_pop_2015_monthly.parquet
#      gadm1_era_tmpmax_un__daily.parquet   (no weight year -> double underscore)
#
# QUICK START
# -----------
#   source("econometrics/download_weighted_climate_data.R")
#
#   # population-weighted CRU mean temperature, countries, yearly, 1960-2019
#   tmp <- wcd_get(variable = "avg. temperature", source = "CRU TS",
#                  geo_resolution = "gadm0", weight = "population density",
#                  weight_year = 2015, time_frequency = "yearly",
#                  years = c(1960, 2019))
#
#   # what exists for a given resolution?
#   wcd_catalog("gadm0")
#
# Requires: curl, jsonlite, arrow (tibble optional).
# =============================================================================

.wcd_need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing required package(s): ", paste(miss, collapse = ", "),
         "\nInstall with install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Public Box bulk-download folders, one per geographical resolution.
# (Source: the "Interested in bulk downloads?" table of the dashboard's
#  Download Data page, https://github.com/CoMoS-SA/climaterepo)
.WCD_BOX_SHARES <- c(
  gadm_world = "q566o1o4xjlin83jvgbupwv2tgbk73t4",
  gadm0      = "qaej7c5swi73fr24fkodtieq9qymh7bz",
  gadm1      = "qjflsyu33qcl1r9xw5ki4jtcmv2z9gx6",
  gadm2      = "1jlbcza7vrsm8r1x1sw80qo0k3q2io6b",
  nuts0      = "jgs8ib9ac372zoy7pgfokotoc2lbhwau",
  nuts1      = "gv2psi0ntn0tvd9huyhopo84zfsfttr4",
  nuts2      = "njpa3noriiajwlyzsy4hni2nkp7am13z",
  nuts3      = "bjz591wg261pw3u1zn2nzknszdilgihd"
)

.WCD_UA <- paste0("Mozilla/5.0 (Windows NT 10.0; Win64; x64) ",
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36")

# Dashboard label -> file-name code. Both sides are accepted as input.
.WCD_VARIABLES <- c(
  "avg. temperature" = "tmp",  "mean temperature" = "tmp",  "temperature" = "tmp",
  "min. temperature" = "tmpmin", "max. temperature" = "tmpmax",
  "precipitation"    = "pre",  "max. wind gust"   = "gust",
  "SPEI1" = "spei01", "SPEI12" = "spei12", "SPEI36" = "spei36",
  tmp = "tmp", tmpmin = "tmpmin", tmpmax = "tmpmax", pre = "pre", gust = "gust",
  spei01 = "spei01", spei1 = "spei01", spei12 = "spei12", spei36 = "spei36"
)

.WCD_SOURCES <- c(
  "CRU TS" = "cru", "CRU" = "cru", cru = "cru",
  "ERA5" = "era", era = "era", era5 = "era",
  "UDelaware" = "dela", "UDel" = "dela", "Delaware" = "dela", dela = "dela",
  "CSIC" = "spei", csic = "spei", spei = "spei"
)

.WCD_WEIGHTS <- c(
  "population density"    = "pop",
  "night lights"          = "lights", "nightlights" = "lights",
  "cropland use"          = "cropland",
  "concurrent population" = "concurrent",
  "unweighted"            = "un",
  pop = "pop", lights = "lights", cropland = "cropland",
  concurrent = "concurrent", un = "un", area = "un", none = "un"
)

# Weights that carry no base year (the dashboard blanks the field for these).
.WCD_YEARLESS_WEIGHTS <- c("un", "concurrent")

.WCD_UNITS <- c(tmp = "degrees Celsius", tmpmin = "degrees Celsius",
                tmpmax = "degrees Celsius", pre = "mm", gust = "m/s",
                spei01 = "unitless", spei12 = "unitless", spei36 = "unitless")

.WCD_FILE_RE <- paste0(
  "^(gadm_world|gadm[0-2]|nuts[0-3])_",   # geo resolution
  "(cru|era|dela|spei)_",                 # source
  "(tmpmin|tmpmax|tmp|pre|gust|spei01|spei12|spei36)_",
  "(pop|lights|cropland|concurrent|un)_", # weight
  "([0-9]{4})?_",                         # weight base year (may be empty)
  "(monthly|daily)\\.parquet$"
)

# In-session memo for folder listings, so repeated calls hit neither Box nor disk.
.wcd_memo <- new.env(parent = emptyenv())

.wcd_msg <- function(verbose, ...) if (isTRUE(verbose)) message(...)

.wcd_match <- function(x, table, what) {
  if (length(x) != 1L || is.na(x)) stop("`", what, "` must be a single value.", call. = FALSE)
  x <- as.character(x)
  hit <- table[match(tolower(trimws(x)), tolower(names(table)))]
  if (is.na(hit)) {
    stop("Unknown ", what, ": '", x, "'.\nAccepted values: ",
         paste(unique(names(table)), collapse = ", "), call. = FALSE)
  }
  unname(hit)
}

# ---------------------------------------------------------------------------
# Cache directory
# ---------------------------------------------------------------------------

#' Directory where downloaded parquet files and folder indexes are kept.
#' Override with options(wcd.cache_dir = "some/path") or the `cache_dir` argument.
wcd_cache_dir <- function(create = TRUE) {
  dir <- getOption("wcd.cache_dir",
                   tools::R_user_dir("weightedclimatedata", "cache"))
  if (create && !dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# ---------------------------------------------------------------------------
# Low-level: public Box shared-folder access
# ---------------------------------------------------------------------------

# Open a session on a public Box shared folder: fetch the share page to collect
# the anonymous session cookies plus the request token, and read off the folder id.
.wcd_box_session <- function(geo_resolution) {
  hash <- .WCD_BOX_SHARES[[geo_resolution]]
  h <- curl::new_handle(followlocation = TRUE, useragent = .WCD_UA)
  res <- curl::curl_fetch_memory(paste0("https://cmu.box.com/s/", hash), handle = h)
  if (res$status_code != 200L) {
    stop("Could not open the public Box folder for '", geo_resolution,
         "' (HTTP ", res$status_code, ").", call. = FALSE)
  }
  html <- rawToChar(res$content)
  Encoding(html) <- "UTF-8"

  token  <- regmatches(html, regexpr('requestToken":"[a-f0-9]+', html))
  folder <- regmatches(html, regexpr('"itemID":[0-9]+,"itemType":"folder"', html))
  if (!length(folder)) folder <- regmatches(html, regexpr('"currentFolderID":[0-9]+', html))
  if (!length(token) || !length(folder)) {
    stop("Box changed its shared-folder page layout: could not read the request ",
         "token / folder id for '", geo_resolution, "'.", call. = FALSE)
  }

  list(handle = h,
       hash   = hash,
       token  = sub('requestToken":"', "", token),
       folder = regmatches(folder, regexpr("[0-9]+", folder)))
}

# List every file in the shared folder. The endpoint pages with `page`/`pageSize`
# (pageSize is capped at 100; larger values are silently rejected).
.wcd_box_list <- function(geo_resolution, verbose = TRUE) {
  ses <- .wcd_box_session(geo_resolution)
  curl::handle_setheaders(ses$handle,
    "Request-Token"     = ses$token,
    "X-Request-Token"   = ses$token,
    "X-Box-EndUser-API" = paste0("sharedName=", ses$hash),
    "Accept"            = "application/json")

  page_size <- 100L
  out <- list()
  for (page in 1:100) {
    url <- sprintf(paste0("https://cmu.app.box.com/app-api/enduserapp/shared-folder",
                          "?folder_id=%s&sortColumn=name&sortDirection=asc",
                          "&pageSize=%d&page=%d"),
                   ses$folder, page_size, page)
    res <- curl::curl_fetch_memory(url, handle = ses$handle)
    if (res$status_code != 200L) {
      stop("Box folder listing failed for '", geo_resolution,
           "' (HTTP ", res$status_code, ").", call. = FALSE)
    }
    items <- jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE)$items
    files <- Filter(function(it) identical(it$type, "file"), items)
    if (length(files)) {
      out[[length(out) + 1L]] <- data.frame(
        file       = vapply(files, function(it) as.character(it$name), character(1)),
        file_id    = vapply(files, function(it) format(it$id, scientific = FALSE), character(1)),
        size_bytes = vapply(files, function(it) as.numeric(it$itemSize), numeric(1)),
        stringsAsFactors = FALSE)
    }
    .wcd_msg(verbose, "  ...listed ", length(items), " item(s) on page ", page)
    if (length(items) < page_size) break
  }
  if (!length(out)) stop("The Box folder for '", geo_resolution, "' looks empty.", call. = FALSE)
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

#' List every WCD file available for a geographical resolution.
#'
#' @param geo_resolution one of gadm_world, gadm0, gadm1, gadm2, nuts0..nuts3
#' @param refresh        TRUE to re-query Box instead of using the cached index
#' @param max_age_days   age above which the cached index is refreshed
#' @return data.frame with file, file_id, size_bytes and the parsed selectors
#'         (source, variable, weight, weight_year, frequency)
wcd_catalog <- function(geo_resolution = "gadm0",
                        refresh        = FALSE,
                        max_age_days   = 7,
                        cache_dir      = wcd_cache_dir(),
                        verbose        = TRUE) {
  .wcd_need(c("curl", "jsonlite"))
  geo_resolution <- match.arg(geo_resolution, names(.WCD_BOX_SHARES))

  if (!refresh && !is.null(.wcd_memo[[geo_resolution]])) return(.wcd_memo[[geo_resolution]])

  idx_file <- file.path(cache_dir, paste0("catalog_", geo_resolution, ".rds"))
  if (!refresh && file.exists(idx_file) &&
      difftime(Sys.time(), file.mtime(idx_file), units = "days") < max_age_days) {
    cat_df <- readRDS(idx_file)
  } else {
    .wcd_msg(verbose, "Indexing the public WCD Box folder for '", geo_resolution, "'...")
    cat_df <- .wcd_box_list(geo_resolution, verbose = verbose)
    m <- regmatches(cat_df$file, regexec(.WCD_FILE_RE, cat_df$file))
    keep <- lengths(m) == 7L
    cat_df <- cat_df[keep, , drop = FALSE]
    m <- do.call(rbind, lapply(m[keep], function(x) x[-1]))
    cat_df$geo_resolution <- m[, 1]
    cat_df$source         <- m[, 2]
    cat_df$variable       <- m[, 3]
    cat_df$weight         <- m[, 4]
    cat_df$weight_year    <- m[, 5]
    cat_df$frequency      <- m[, 6]
    cat_df <- cat_df[order(cat_df$file), ]
    rownames(cat_df) <- NULL
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(cat_df, idx_file)
    .wcd_msg(verbose, "Indexed ", nrow(cat_df), " file(s).")
  }
  .wcd_memo[[geo_resolution]] <- cat_df
  cat_df
}

#' Build the canonical WCD file name for a selection (dashboard naming scheme).
wcd_file_name <- function(geo_resolution, source, variable, weight, weight_year,
                          frequency = "monthly") {
  if (weight %in% .WCD_YEARLESS_WEIGHTS) weight_year <- ""
  paste0(geo_resolution, "_", source, "_", variable, "_", weight, "_",
         weight_year, "_", frequency, ".parquet")
}

# Turn a selection into a catalogue row, with a helpful error when it does not exist.
.wcd_resolve <- function(cat_df, geo_resolution, source, variable, weight,
                         weight_year, frequency, strict = FALSE, verbose = TRUE) {
  want <- wcd_file_name(geo_resolution, source, variable, weight, weight_year, frequency)
  hit <- cat_df[cat_df$file == want, , drop = FALSE]
  if (nrow(hit) == 1L) return(hit)

  # Same selection but any weight base year: the dashboard silently forces the
  # only year it has (e.g. 2015 for tmpmin/tmpmax/gust and for nuts*/gadm_world).
  alt <- cat_df[cat_df$source == source & cat_df$variable == variable &
                cat_df$weight == weight & cat_df$frequency == frequency, , drop = FALSE]
  if (!strict && nrow(alt) == 1L) {
    .wcd_msg(verbose, "No '", want, "'; using the only available weight base year ",
             "for this selection: '", alt$file, "'.")
    return(alt)
  }

  avail <- cat_df[cat_df$source == source & cat_df$variable == variable, , drop = FALSE]
  msg <- paste0("No WCD file matches this selection.\n  wanted: ", want)
  if (nrow(alt)) {
    msg <- paste0(msg, "\n  available weight base years: ",
                  paste(ifelse(alt$weight_year == "", "(none)", alt$weight_year), collapse = ", "))
  } else if (nrow(avail)) {
    msg <- paste0(msg, "\n  available for ", source, "/", variable, " at ", geo_resolution, ":\n    ",
                  paste(utils::head(avail$file, 25), collapse = "\n    "),
                  if (nrow(avail) > 25) paste0("\n    ... and ", nrow(avail) - 25, " more") else "")
  } else {
    msg <- paste0(msg, "\n  '", source, "' has no '", variable, "' at ", geo_resolution,
                  ".\n  available source/variable pairs: ",
                  paste(unique(paste0(cat_df$source, "/", cat_df$variable)), collapse = ", "))
  }
  stop(msg, "\nCall wcd_catalog('", geo_resolution, "') to browse everything on offer.",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

#' Download the raw parquet file behind a selection and return its local path.
#'
#' Same selectors as wcd_get(); no post-processing is applied.
wcd_download <- function(variable       = "avg. temperature",
                         source         = "CRU TS",
                         geo_resolution = "gadm0",
                         weight         = "population density",
                         weight_year    = 2015,
                         time_frequency = "monthly",
                         cache_dir      = wcd_cache_dir(),
                         refresh        = FALSE,
                         verbose        = TRUE) {
  .wcd_need(c("curl", "jsonlite"))
  geo_resolution <- match.arg(geo_resolution, names(.WCD_BOX_SHARES))
  variable <- .wcd_match(variable, .WCD_VARIABLES, "variable")
  source   <- .wcd_match(source,   .WCD_SOURCES,   "source")
  weight   <- .wcd_match(weight,   .WCD_WEIGHTS,   "weight")
  time_frequency <- match.arg(tolower(time_frequency), c("yearly", "monthly", "daily"))
  # `yearly` is not a stored frequency: the dashboard derives it from monthly files.
  file_freq <- if (time_frequency == "yearly") "monthly" else time_frequency
  weight_year <- if (is.null(weight_year) || is.na(weight_year)) "" else as.character(weight_year)

  cat_df <- wcd_catalog(geo_resolution, cache_dir = cache_dir, verbose = verbose)
  row <- .wcd_resolve(cat_df, geo_resolution, source, variable, weight,
                      weight_year, file_freq, verbose = verbose)

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_dir, row$file)
  pretty <- format(structure(row$size_bytes, class = "object_size"), units = "auto")
  if (file.exists(dest) && !refresh && file.size(dest) == row$size_bytes) {
    .wcd_msg(verbose, "Using cached ", row$file, " (", pretty, ").")
    return(dest)
  }

  url <- sprintf(paste0("https://cmu.box.com/index.php?rm=box_download_shared_file",
                        "&shared_name=%s&file_id=f_%s"),
                 .WCD_BOX_SHARES[[geo_resolution]], row$file_id)
  .wcd_msg(verbose, "Downloading ", row$file, " (", pretty, ")...")
  tmp <- paste0(dest, ".part")
  on.exit(unlink(tmp), add = TRUE)
  # curl's progress bar is only readable on an interactive console; in a script
  # it would otherwise flood the log with one line per chunk.
  curl::curl_download(url, tmp, quiet = !(verbose && interactive()),
                      handle = curl::new_handle(followlocation = TRUE, useragent = .WCD_UA))

  # Box answers with an HTML page rather than an HTTP error when a share dies.
  if (!identical(rawToChar(readBin(tmp, "raw", 4L)), "PAR1")) {
    stop("The download did not return a parquet file -- the public Box share for '",
         geo_resolution, "' may have moved. Check ",
         "https://weightedclimatedata.streamlit.app/Download_Data", call. = FALSE)
  }
  file.rename(tmp, dest)
  dest
}

# ---------------------------------------------------------------------------
# Country / unit lookup
# ---------------------------------------------------------------------------

#' GID_0 <-> country-name lookup used by the dashboard (from the project repo).
wcd_countries <- function(geo_resolution = "gadm0", cache_dir = wcd_cache_dir()) {
  fname <- if (grepl("^nuts", geo_resolution)) "country_list_nuts0.csv" else "country_list.csv"
  dest <- file.path(cache_dir, fname)
  if (!file.exists(dest)) {
    .wcd_need("curl")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    curl::curl_download(
      paste0("https://raw.githubusercontent.com/CoMoS-SA/climaterepo/main/poly/", fname),
      dest, quiet = TRUE)
  }
  utils::read.csv(dest, stringsAsFactors = FALSE)
}

# Columns of a WCD file are unit ids: ISO3 at gadm0 ("ITA"), GID_1/GID_2 with dots
# replaced by underscores at finer GADM levels ("ITA_1_1"), NUTS codes at nuts*
# ("ITC1"). Country filtering therefore reduces to a prefix test.
.wcd_country_filter <- function(cols, countries, geo_resolution, cache_dir) {
  lut <- try(wcd_countries(geo_resolution, cache_dir), silent = TRUE)
  codes <- character(0)
  for (ct in countries) {
    if (!inherits(lut, "try-error") && tolower(ct) %in% tolower(lut$COUNTRY)) {
      codes <- c(codes, lut$GID_0[match(tolower(ct), tolower(lut$COUNTRY))])
    } else {
      codes <- c(codes, ct)
    }
  }
  codes <- toupper(unique(codes))
  prefix <- toupper(sub("_.*$", "", cols))
  if (grepl("^nuts", geo_resolution)) {
    keep <- vapply(toupper(cols), function(cl) any(startsWith(cl, codes)), logical(1))
  } else {
    keep <- prefix %in% codes
    unknown <- setdiff(codes, prefix)
    if (length(unknown)) warning("No columns for: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  cols[keep]
}

# ---------------------------------------------------------------------------
# Post-processing (mirrors the dashboard)
# ---------------------------------------------------------------------------

# Stored dates are strings: "X190101" (monthly) or "X19010101" (daily).
.wcd_parse_dates <- function(x, frequency) {
  x <- sub("^X", "", as.character(x))
  if (frequency == "daily") {
    as.Date(x, format = "%Y%m%d")
  } else {
    as.Date(paste0(x, "01"), format = "%Y%m%d")
  }
}

# Yearly aggregators, exactly as the dashboard applies them.
.wcd_year_fun <- function(variable) {
  switch(variable,
    pre    = function(v) if (anyNA(v)) NA_real_ else sum(v),
    tmp    = function(v) mean(v),
    tmpmin = function(v) min(v),
    tmpmax = ,
    gust   = function(v) max(v),
    stop("The dashboard offers no yearly aggregation for '", variable,
         "' (SPEI is monthly only). Use time_frequency = \"monthly\".", call. = FALSE))
}

# Collapse rows within each period key, keeping the column layout.
.wcd_by_period <- function(mat, keys, fun) {
  u <- unique(keys)
  out <- vapply(u, function(k) apply(mat[keys == k, , drop = FALSE], 2, fun),
                numeric(ncol(mat)))
  out <- if (ncol(mat) == 1L) matrix(out, ncol = 1L) else t(out)
  dimnames(out) <- list(NULL, colnames(mat))
  list(values = out, keys = u)
}

.wcd_aggregate <- function(mat, dates, fun) {
  agg <- .wcd_by_period(mat, format(dates, "%Y"), fun)
  list(values = agg$values, dates = as.Date(paste0(agg$keys, "-12-31")))
}

# Threshold indicators: count (or cumulative excess) of days beyond a threshold,
# summed by month or year. Daily ERA5 files only, as in the dashboard.
.wcd_threshold <- function(mat, dates, kind, value, to) {
  lim <- if (kind == "percentile") {
    apply(mat, 2, stats::quantile, probs = value / 100, na.rm = TRUE, names = FALSE)
  } else {
    rep(value, ncol(mat))
  }
  over <- sweep(mat, 2, lim, ">") & !is.na(mat)
  cnt <- if (kind == "cumulative") {
    ex <- sweep(mat, 2, lim, "-")
    ex[!over] <- 0
    ex
  } else {
    over * 1
  }
  keys <- if (to == "yearly") format(dates, "%Y") else format(dates, "%Y-%m")
  agg <- .wcd_by_period(cnt, keys, function(v) sum(v, na.rm = TRUE))
  list(values = agg$values,
       dates  = if (to == "yearly") as.Date(paste0(agg$keys, "-12-31"))
                else as.Date(paste0(agg$keys, "-01")))
}

.wcd_id_col <- function(geo_resolution) {
  switch(geo_resolution,
    gadm_world = "unit", gadm0 = "GID_0", gadm1 = "GID_1", gadm2 = "GID_2",
    nuts0 = "GID_0", nuts1 = "GID_1", nuts2 = "GID_2", nuts3 = "GID_3", "unit")
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Download and assemble a Weighted Climate Dataset selection.
#'
#' Reproduces what the dashboard's Download Data page hands you, from R.
#'
#' @param variable       "avg. temperature", "min. temperature", "max. temperature",
#'                       "precipitation", "max. wind gust", "SPEI1", "SPEI12",
#'                       "SPEI36" (file codes tmp/tmpmin/tmpmax/pre/gust/spei* also work)
#' @param source         "CRU TS", "ERA5", "UDelaware" or "CSIC" (SPEI)
#' @param geo_resolution "gadm_world", "gadm0" (countries), "gadm1", "gadm2",
#'                       "nuts0", "nuts1", "nuts2", "nuts3"
#' @param weight         "population density", "night lights", "cropland use",
#'                       "concurrent population", "unweighted"
#' @param weight_year    2000, 2005, 2010 or 2015 (ignored for unweighted /
#'                       concurrent population, and for selections stored at a
#'                       single base year)
#' @param time_frequency "yearly", "monthly" or "daily" ("yearly" is aggregated
#'                       from the monthly file the way the dashboard does it;
#'                       "daily" exists for ERA5 only)
#' @param years          NULL, or c(start, end) inclusive, or a vector of years
#' @param countries      NULL, or ISO3 codes / country names to keep. At gadm1+
#'                       and nuts* this keeps every unit inside those countries.
#' @param units          NULL, or exact unit ids to keep (e.g. c("ITA_1_1"))
#' @param format         "long" (default) or "wide"
#' @param threshold      NULL, or a threshold value applied to the daily ERA5
#'                       file; the result counts days beyond it per period
#' @param threshold_kind "percentile" (default), "absolute" or "cumulative"
#' @param cache_dir      where parquet files are kept between sessions
#' @param refresh        TRUE to re-download even if cached
#' @param verbose        progress messages
#'
#' @return A data.frame (tibble when tibble is installed). Long format has
#'         columns date, <unit id>, <variable>; wide format has date plus one
#'         column per unit. attr(x, "wcd") records the resolved selection.
wcd_get <- function(variable       = "avg. temperature",
                    source         = "CRU TS",
                    geo_resolution = "gadm0",
                    weight         = "population density",
                    weight_year    = 2015,
                    time_frequency = "monthly",
                    years          = NULL,
                    countries      = NULL,
                    units          = NULL,
                    format         = c("long", "wide"),
                    threshold      = NULL,
                    threshold_kind = c("percentile", "absolute", "cumulative"),
                    cache_dir      = wcd_cache_dir(),
                    refresh        = FALSE,
                    verbose        = TRUE) {
  .wcd_need(c("curl", "jsonlite", "arrow"))
  format         <- match.arg(format)
  threshold_kind <- match.arg(threshold_kind)
  geo_resolution <- match.arg(geo_resolution, names(.WCD_BOX_SHARES))
  var_code <- .wcd_match(variable, .WCD_VARIABLES, "variable")
  src_code <- .wcd_match(source,   .WCD_SOURCES,   "source")
  wgt_code <- .wcd_match(weight,   .WCD_WEIGHTS,   "weight")
  time_frequency <- match.arg(tolower(time_frequency), c("yearly", "monthly", "daily"))

  # Threshold indicators are computed on daily data, as in the dashboard.
  use_threshold <- !is.null(threshold)
  if (use_threshold) {
    if (time_frequency == "daily") {
      stop("A threshold produces monthly or yearly counts; set time_frequency to one of those.",
           call. = FALSE)
    }
    if (src_code != "era") {
      stop("Threshold indicators need daily data, i.e. source = \"ERA5\".", call. = FALSE)
    }
  }
  file_freq <- if (use_threshold) "daily" else if (time_frequency == "yearly") "monthly" else time_frequency
  year_fun <- if (time_frequency == "yearly" && !use_threshold) .wcd_year_fun(var_code) else NULL

  path <- wcd_download(variable = var_code, source = src_code,
                       geo_resolution = geo_resolution, weight = wgt_code,
                       weight_year = weight_year, time_frequency = file_freq,
                       cache_dir = cache_dir, refresh = refresh, verbose = verbose)
  used <- regmatches(basename(path), regexec(.WCD_FILE_RE, basename(path)))[[1]]

  .wcd_msg(verbose, "Reading ", basename(path), "...")
  dat <- as.data.frame(arrow::read_parquet(path))

  date_col <- if ("Date" %in% names(dat)) "Date" else names(dat)[1]
  dates <- .wcd_parse_dates(dat[[date_col]], file_freq)
  dat[[date_col]] <- NULL
  # Placeholder columns present in some source files.
  dat <- dat[, !is.na(names(dat)) & !(names(dat) %in% c("NA", "?", "")), drop = FALSE]

  # ERA5 precipitation is stored in metres at these resolutions; the dashboard
  # divides by 1000 to report mm.
  if (src_code == "era" && var_code == "pre" &&
      geo_resolution %in% c("gadm_world", "nuts0", "nuts1", "nuts2", "nuts3")) {
    dat[] <- lapply(dat, function(v) v / 1000)
  }

  keep <- names(dat)
  if (!is.null(countries)) keep <- .wcd_country_filter(keep, countries, geo_resolution, cache_dir)
  if (!is.null(units))     keep <- intersect(keep, units)
  if (!length(keep)) stop("No geographical unit left after filtering.", call. = FALSE)
  dat <- dat[, keep, drop = FALSE]

  if (!is.null(years)) {
    yrs <- if (length(years) == 2L) seq(min(years), max(years)) else years
    sel <- as.integer(format(dates, "%Y")) %in% yrs
    if (!any(sel)) {
      stop("No observation in ", paste(range(yrs), collapse = "-"), "; this file covers ",
           paste(range(as.integer(format(dates, "%Y"))), collapse = "-"), ".", call. = FALSE)
    }
    dat <- dat[sel, , drop = FALSE]
    dates <- dates[sel]
  }

  mat <- as.matrix(dat)
  storage.mode(mat) <- "double"
  if (use_threshold) {
    agg <- .wcd_threshold(mat, dates, threshold_kind, threshold,
                          if (time_frequency == "yearly") "yearly" else "monthly")
    mat <- agg$values
    dates <- agg$dates
  } else if (time_frequency == "yearly") {
    agg <- .wcd_aggregate(mat, dates, year_fun)
    mat <- agg$values
    dates <- agg$dates
  }

  id_col  <- .wcd_id_col(geo_resolution)
  val_col <- if (use_threshold) paste0(var_code, "_over_threshold") else var_code

  if (format == "long") {
    out <- data.frame(
      date  = rep(dates, times = ncol(mat)),
      unit  = rep(colnames(mat), each = nrow(mat)),
      value = as.vector(mat),
      stringsAsFactors = FALSE)
    names(out) <- c("date", id_col, val_col)
  } else {
    out <- data.frame(date = dates, mat, check.names = FALSE, stringsAsFactors = FALSE)
  }
  # Yearly series carry no meaningful day: report the year itself.
  if (time_frequency == "yearly") out$date <- as.integer(format(out$date, "%Y"))

  attr(out, "wcd") <- list(
    file = basename(path), geo_resolution = geo_resolution, source = src_code,
    variable = var_code, weight = wgt_code, weight_year = used[6],
    stored_frequency = file_freq, time_frequency = time_frequency,
    units_of_measure = unname(.WCD_UNITS[var_code]),
    threshold = threshold, threshold_kind = if (use_threshold) threshold_kind else NULL,
    citation = paste("Gortan, Testa, Fagiolo & Lamperti (2024),",
                     "Scientific Data, doi:10.1038/s41597-024-03304-1"))

  .wcd_msg(verbose, "Returned ", nrow(out), " row(s) x ", ncol(mat), " unit(s), ",
           time_frequency, ", in ", unname(.WCD_UNITS[var_code]), ".")
  if (requireNamespace("tibble", quietly = TRUE)) {
    a <- attr(out, "wcd")
    out <- tibble::as_tibble(out)
    attr(out, "wcd") <- a
  }
  out
}
