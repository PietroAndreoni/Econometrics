# Shared setup of xi_exposure_profile.R and xi_exposure_bootstrap.R: the
# exposure specification, its data and terms from the response explorer, and
# the fit at given income elasticities (xi_T, xi_P).
#
#   PWT x ERA5 (population weights), deviations from the 20-year trailing mean,
#   trailing characteristics (20-year climate, income of the 5 years before),
#   country and year effects without trends, last year's log GDP per capita
#   as a control, the explorer's controls and outlier rule (r5).
#   T: linear arms (warm, cold), each through the country's exposure slopes of
#      degree days of Tmax above 29 C and of Tmin below 0 C (absolute units).
#   P: quadratic arms, each also interacted with the trailing climate of P and
#      the agricultural share (SD units).
#   Every T (P) term times (Ybar / Y_it)^(1 - xi_T) (xi_P), Y the trailing
#   income: the dollar loss per person proportional to Y^xi (xi = 1: equal %,
#   the unscaled terms; xi = 0: equal dollars).
#
# The controls are the same at every xi: the explorer's level terms of the
# trailing characteristics, with the trailing log income y in place of the
# scaling's own level (Ybar / Y)^(1 - xi), which the explorer adds with a
# scaled model. That level changes with xi and vanishes at xi = 1 (a
# constant), and next to last year's log GDP it acts as a momentum control
# (last year's level less the 5-year mean): a profile with it compares models
# with and without that control rather than values of xi. With y always in,
# the profile changes only xi (as in handout_preferred_specification.R).

suppressPackageStartupMessages(library(fixest))
`%||%` <- function(a, b) if (is.null(a)) b else a
APP_DIR <- normalizePath(file.path("apps", "response_explorer"))
for (f in c("response.R", "projection.R", "fit.R", "exposure.R")) source(file.path(APP_DIR, "R", f))
XI_RANGE <- c(-5, 5)
XI_BOUNDS <- c(-4, 4)
MODE <- Sys.getenv("EXPOSURE_MODE", "country")
stopifnot(MODE %in% c("country", "smooth"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
out_dir <- file.path("results", "xi_exposure")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SPEC <- list(
  gdp = "PWT", climate = "ERA5", weight = "population", baseline = "20-year", style = "simple",
  effects = "both", fe = "unit", convergence = "lag1", moderators = "trailing", outliers = "r5",
  form = c(TM = "Linear arms", RR = "Quadratic arms"), symmetric = FALSE,
  units = c(TM = "absolute", RR = "SD"),
  inter = list(TM = character(0), RR = c("clim", "agr")), shape = "linear",
  exposure = list(TM = list(mode = MODE, channels = "both",
                            hi = list(src = "tmpmax", type = "abs", value = 29, measure = "dd"),
                            lo = list(src = "tmpmin", type = "abs", value = 0, measure = "dd"))),
  # Both variables scaled (the value of xi here only sets the structure).
  income_form = c(TM = "dollars", RR = "dollars"), xi = c(TM = 0.5, RR = 0.5)
)

# explorer_model() spliced to return its data, terms and fixed effects.
lines <- readLines(file.path(APP_DIR, "R", "fit.R"))
s0 <- grep("^explorer_model <- function", lines)
s1 <- s0 - 1 + which(lines[s0:length(lines)] == "}")[1]
fn <- paste(lines[s0:s1], collapse = "\n")
fn <- sub("explorer_model <- function", "model_parts <- function", fn)
fn <- sub("if (!is.null(model_cache[[key]])) return(model_cache[[key]])", "", fn, fixed = TRUE)
fn <- sub("m <- fit(o$shape)",
          "return(list(d = d, terms = terms_of(o$shape), fe = fixed_effects_formula(o$fe, panels$fe_slopes, o$effects)))",
          fn, fixed = TRUE)
stopifnot(grepl("return(list(d = d, terms", fn, fixed = TRUE))
eval(parse(text = fn))
parts <- model_parts(SPEC, P)
FE <- parts$fe
DEV <- grep("^dev", parts$terms, value = TRUE)
# The controls: the explorer's, with the trailing log income in place of the
# scaling's level.
OTHER <- unique(c(setdiff(parts$terms[!startsWith(parts$terms, "dev")], c("idol_TM", "idol_RR")), "inc"))
# (The variables of the controls, e.g. i(fsu_collapse_year, ref = 'none').)
other_vars <- unique(unlist(lapply(OTHER, function(t) all.vars(str2lang(t)))))
need <- unique(c(".outcome", "GID_0", "GID_1", "year", "inc", other_vars, unlist(strsplit(DEV, ":", fixed = TRUE))))
need <- setdiff(need, "idol")
stopifnot(all(need %in% names(parts$d)))
D <- as.data.frame(parts$d)[, need]

# The scaled deviation terms at (xi_T, xi_P) as columns t_1.. of x.
with_terms <- function(xi_T, xi_P, x) {
  x$idol_TM <- exp(-(1 - xi_T) * x$inc)
  x$idol_RR <- exp(-(1 - xi_P) * x$inc)
  for (j in seq_along(DEV)) {
    v <- substr(DEV[j], 4, 5)
    p <- strsplit(DEV[j], ":", fixed = TRUE)[[1]]
    p[p == "idol"] <- paste0("idol_", v)
    x[[paste0("t_", j)]] <- Reduce(`*`, lapply(p, function(z) x[[z]]))
  }
  x
}
fit_xi <- function(xi_T, xi_P, x = D, vcov = "iid") {
  x <- with_terms(xi_T, xi_P, x)
  f <- stats::as.formula(paste(".outcome ~", paste(c(paste0("t_", seq_along(DEV)), OTHER), collapse = " + "), "|", FE))
  m <- fixest::feols(f, data = x, notes = FALSE, warn = FALSE, vcov = vcov)
  m$dev_names <- stats::setNames(DEV, paste0("t_", seq_along(DEV)))
  m
}
rss_at <- function(xi_T, xi_P, x = D) sum(stats::resid(fit_xi(xi_T, xi_P, x))^2)
# NLS by coordinates: a coarse grid of xi_T (step over XI_BOUNDS) at the
# current xi_P, a continuous search around its best point, the same for xi_P,
# then xi_T once more.
nls_xi <- function(x = D, start = c(TM = 1, RR = 1), step = 0.5) {
  one <- function(v, at) {
    f <- function(z) if (v == "TM") rss_at(z, at, x) else rss_at(at, z, x)
    xs <- seq(XI_BOUNDS[1], XI_BOUNDS[2], step)
    i <- which.min(vapply(xs, f, 0))
    if (i %in% c(1, length(xs))) return(xs[i])
    stats::optimize(f, xs[i] + c(-1, 1) * step, tol = 1e-3)$minimum
  }
  xT <- one("TM", start[["RR"]])
  xP <- one("RR", xT)
  xT <- one("TM", xP)
  c(TM = xT, RR = xP)
}
