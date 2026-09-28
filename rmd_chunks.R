# Reuse code chunks from an R Markdown notebook inside a plain script.
#
# `test_functions.Rmd` holds the shared specification (`fixed_effects`, `model`,
# `make_formula`, `fit_climate_model`) and the estimation-sample builder
# (`build_dat`). Analysis scripts need those definitions without knitting the
# whole notebook, so the named chunks are extracted and evaluated in a private
# environment. Assignments listed in `skip_assignments` are not evaluated, which
# lets a caller import `build_dat()` without triggering the notebook's own data
# builds.

extract_rmd_chunks <- function(rmd_path, labels) {
  lines <- readLines(rmd_path, warn = FALSE)
  chunk_opens <- grep("^```\\{r", lines)
  chunk_closes <- grep("^```\\s*$", lines)

  chunks <- list()
  for (open in chunk_opens) {
    close <- chunk_closes[chunk_closes > open]
    if (!length(close)) {
      stop("Unterminated code chunk at line ", open, " of ", rmd_path)
    }
    header <- trimws(sub("\\}\\s*$", "", sub("^```\\{r", "", lines[[open]])))
    label <- gsub("^[\"']|[\"']$", "", trimws(sub(",.*$", "", header)))
    if (nzchar(label) && close[[1]] > open + 1L) {
      chunks[[label]] <- lines[(open + 1L):(close[[1]] - 1L)]
    }
  }

  missing_chunks <- setdiff(labels, names(chunks))
  if (length(missing_chunks)) {
    stop(
      "Chunks not found in ", rmd_path, ": ",
      paste(missing_chunks, collapse = ", ")
    )
  }
  chunks[labels]
}

eval_rmd_chunks <- function(
    rmd_path,
    labels,
    envir,
    skip_assignments = character()
) {
  chunks <- extract_rmd_chunks(rmd_path, labels)
  for (label in labels) {
    for (expression in parse(text = chunks[[label]])) {
      assigns_skipped_object <- is.call(expression) &&
        identical(as.character(expression[[1L]]), "<-") &&
        is.name(expression[[2L]]) &&
        as.character(expression[[2L]]) %in% skip_assignments
      if (!assigns_skipped_object) {
        eval(expression, envir = envir)
      }
    }
  }
  invisible(envir)
}

# Convenience wrapper: return an environment holding the notebook's specification
# constants and data builders, with the notebook's own `data` / `data_iso3`
# assignments suppressed so the caller controls when and how the panel is built.
load_notebook_env <- function(
    rmd_path,
    econ_data = "DOSE_V2_11",
    labels = c("shared-specification", "build-data"),
    skip_assignments = c("data", "data_iso3"),
    climate_velocity_years = 5L
) {
  notebook <- new.env(parent = globalenv())
  notebook$params <- list(
    project_setup_script = NULL,
    climate_velocity_years = climate_velocity_years,
    econ_data = econ_data
  )
  eval_rmd_chunks(
    rmd_path,
    labels = labels,
    envir = notebook,
    skip_assignments = skip_assignments
  )
  notebook
}
