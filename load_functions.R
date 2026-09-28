# Attach the packages the functions/ scripts use unqualified, then source every
# script in functions/ into the global environment. Analysis scripts and
# notebooks start with:
#   source("load_functions.R")
#
# Other packages are called with `pkg::` and only need to be installed:
# arrow, curl, jsonlite, purrr, readxl, splines, tibble, zoo.

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(ggplot2)
})

local({
  functions_dir <- if (dir.exists("functions")) {
    "functions"
  } else {
    file.path("..", "functions")
  }
  if (!dir.exists(functions_dir)) {
    stop("functions/ not found; run from the project root.")
  }
  for (script in sort(list.files(functions_dir, pattern = "\\.R$",
                                 full.names = TRUE))) {
    source(script, local = globalenv())
  }
})
