# Data paths and the on-disk panel cache.

# rmarkdown may knit from either the project root or a subdirectory, so data
# paths are resolved against both.
data_root <- function() {
  if (dir.exists("data")) "data" else file.path("..", "data")
}

data_file <- function(name) {
  path <- file.path(data_root(), name)
  if (!file.exists(path)) {
    stop("Raw input not found: ", path)
  }
  path
}

# Derived panels are cached on disk so a session does not repeat the climate
# download, the rolling moments, or the Excel read. Each cache carries a JSON
# sidecar recording exactly the specification that produced it; when the current
# settings no longer match the sidecar the panel is rebuilt and both are
# rewritten. Deleting the cache directory is always safe.
# The sidecar records the settings that produced a cache, but settings alone do
# not pin the result: editing a builder while the constants stay put would
# otherwise leave a stale panel in place. Including the builders' own source in
# the specification makes any change to the derivation invalidate the cache.
code_fingerprint <- function(...) {
  vapply(
    list(...),
    function(f) paste(deparse(body(f)), collapse = "\n"),
    character(1)
  )
}

cached_panel <- function(tag, spec, build) {
  cache_dir <- file.path(data_root(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  panel_path <- file.path(cache_dir, paste0(tag, ".parquet"))
  spec_path <- file.path(cache_dir, paste0(tag, ".json"))
  spec_json <- as.character(
    jsonlite::toJSON(spec, auto_unbox = TRUE, digits = NA, null = "null")
  )

  if (file.exists(panel_path) && file.exists(spec_path)) {
    stored <- paste(readLines(spec_path, warn = FALSE), collapse = "")
    if (identical(stored, spec_json)) {
      message("Reusing cached ", basename(panel_path))
      return(tibble::as_tibble(arrow::read_parquet(panel_path)))
    }
    message("Specification changed, rebuilding ", basename(panel_path))
  }

  panel <- build()
  arrow::write_parquet(panel, panel_path)
  writeLines(spec_json, spec_path)
  panel
}
