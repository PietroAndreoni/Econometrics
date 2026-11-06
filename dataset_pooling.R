# One estimate from all the datasets: the preferred specification of
# handout_preferred_specification.R (linear and quadratic arms; T scaled by
# (Ybar / Y)^(1 - xi) with its climate interaction, P with climate and share;
# xi as there, read from its estimates.rds) on the six GDP x climate panels
# (PWT, WB x ERA5, CRU TS, UDelaware; population weights), each panel's rows
# and outliers as in the explorer. Saved to results/dataset_pooling/pooling.rds.
#
# The six panels measure the same country-years, so they are not six samples:
# they are fitted jointly as one stacked regression in which everything is
# panel-specific (fixed effects, trends, level terms, controls and the weather
# terms), which gives each panel's own estimates with their joint covariance,
# errors clustered by country across the panels. From it:
#   - clustered Wald tests of equal weather terms across the six panels, across
#     the GDP sources (same climate) and across the climate products (same GDP);
#   - the pooled weather terms: the equal-weight average of the six, the
#     efficient minimum-distance combination (weights from the joint
#     covariance; its J statistic is the Wald test of equality), and the
#     restricted stacked regression (weather terms shared, the rest
#     panel-specific).
# Each on the panels' own rows and on the country-years all six panels have
# (there the panels differ only in the measurements; income and share centred
# on those rows, as in the handout's common rows).
#
#   Rscript dataset_pooling.R

suppressPackageStartupMessages(library(fixest))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- file.path("apps", "response_explorer")
source(file.path(app_dir, "R", "response.R"))
source(file.path(app_dir, "R", "fit.R"))
P <- readRDS(file.path("results", "response_explorer", "panels.rds"))
XI <- readRDS(file.path("results", "handout_preferred_specification", "estimates.rds"))$xi_preferred

# The explorer's model data without the fit (as in the handout).
lines <- readLines(file.path(app_dir, "R", "fit.R"))
s0 <- grep("^explorer_model <- function", lines)
s1 <- s0 - 1 + which(lines[s0:length(lines)] == "}")[1]
fn <- paste(lines[s0:s1], collapse = "\n")
fn <- sub("explorer_model <- function", "model_data <- function", fn)
fn <- sub("if (!is.null(model_cache[[key]])) return(model_cache[[key]])", "", fn, fixed = TRUE)
fn <- sub("m <- fit(o$shape)", "return(d)", fn, fixed = TRUE)
stopifnot(grepl("return(d)", fn, fixed = TRUE))
eval(parse(text = fn))

GDPS <- c("PWT", "WB")
CLIMATES <- c(ERA5 = "ERA5", CRU = "CRU TS", UDel = "UDelaware")
PANELS <- expand.grid(climate = names(CLIMATES), gdp = GDPS, stringsAsFactors = FALSE)[, c("gdp", "climate")]
PANELS$id <- paste0(PANELS$gdp, "_", PANELS$climate)
K <- nrow(PANELS)

# A panel's data as in the handout: deviations from the 20-year trailing mean
# over the country's horizon SD, demeaned by its quadratic trend, their arms
# and squared arms; trailing climate, income and share.
data_of <- function(gdp, climate) {
  d <- model_data(list(gdp = gdp, climate = climate, weight = "population", units = "absolute",
                       baseline = "20-year", form = "Linear arms", shape = "linear", fe = "quadratic",
                       moderators = "trailing", outcome = "log",
                       inter = list(TM = c("clim", "inc", "agr"), RR = c("clim", "inc", "agr"))), P)
  for (v in c("TM", "RR")) {
    x <- d[[paste0("dm_", v, "_lag_20_unstd_quadratic")]] / d[[paste0("sd_", v, "_all")]]
    d[[paste0("p", v)]] <- pmax(x, 0)
    d[[paste0("n", v)]] <- pmax(-x, 0)
    d[[paste0("p", v, "2")]] <- pmax(x, 0)^2
    d[[paste0("n", v, "2")]] <- pmax(-x, 0)^2
  }
  d$g <- d[[P$outcome]]
  d
}
# Income and share less their means over the rows, and the scaling of T.
centre <- function(d) {
  d$y <- d$inc_trailing - mean(d$inc_trailing)
  d$agr <- d$agr_level - mean(d$agr_level)
  d$r <- exp(-(1 - XI) * d$y)
  d
}
raw <- lapply(seq_len(K), function(k) data_of(PANELS$gdp[k], CLIMATES[[PANELS$climate[k]]]))
names(raw) <- PANELS$id
common_keys <- Reduce(intersect, lapply(raw, function(d) paste(d$GID_0, d$year)))
SAMPLES <- list(own = lapply(raw, centre),
                common = lapply(raw, function(d) centre(d[paste(d$GID_0, d$year) %in% common_keys, ])))

# The terms of the preferred specification (names as in the handout).
arm_terms <- function(v, modifiers = character(0), scale = NULL, power = "") {
  out <- list()
  for (a in paste0(c("p", "n"), v, power)) {
    out[[a]] <- c(a, scale)
    for (z in modifiers) out[[paste0(a, ":", z)]] <- c(a, z, scale)
  }
  if (!is.null(scale)) names(out) <- paste0(names(out), ":r")
  names(out) <- gsub(":", "_X_", names(out))
  out
}
SPECS <- list(linear = c(arm_terms("TM", "clim_TM", scale = "r"), arm_terms("RR", c("clim_RR", "agr"))),
              quadratic = c(arm_terms("TM", "clim_TM", scale = "r", power = "2"),
                            arm_terms("RR", c("clim_RR", "agr"), power = "2")))
LEVELS <- c("clim_TM", "clim_RR", "y")
NUMERIC_CONTROLS <- grep("^i\\(", P$controls, value = TRUE, invert = TRUE)
FACTOR_CONTROLS <- sub("^i\\(([a-z_]+),.*$", "\\1", grep("^i\\(", P$controls, value = TRUE))

# The stacked data of a sample and specification: each panel's terms, levels
# and controls as columns of their own (zero on the other panels' rows), the
# collapse-year controls as panel-specific categories, fixed effects and
# trends by panel. The weather terms also as shared columns (for the
# restricted fit).
stack_of <- function(sample, terms) {
  do.call(rbind, lapply(seq_len(K), function(k) {
    d <- sample[[k]]
    for (t in names(terms)) d[[t]] <- Reduce(`*`, lapply(terms[[t]], function(z) d[[z]]))
    out <- data.frame(panel = PANELS$id[k], GID_0 = d$GID_0, year_k = paste(k, d$year),
                      unit_k = paste(k, d$GID_1), t = d$t, t2 = d$t2, g = d$g)
    for (t in names(terms)) out[[t]] <- d[[t]]
    for (j in seq_len(K)) for (z in c(names(terms), LEVELS, NUMERIC_CONTROLS)) {
      out[[paste0(z, "__", PANELS$id[j])]] <- if (j == k) d[[z]] else 0
    }
    for (z in FACTOR_CONTROLS) out[[paste0(z, "_k")]] <- ifelse(d[[z]] == "none", "none", paste(k, d[[z]]))
    out
  }))
}
FE <- "year_k + unit_k[t, t2]"
factor_terms <- paste0("i(", FACTOR_CONTROLS, "_k, ref = 'none')")
panel_terms <- function(z) as.vector(outer(z, PANELS$id, paste, sep = "__"))
fit_stack <- function(s, rhs) {
  f <- as.formula(paste("g ~", paste(c(rhs, factor_terms), collapse = " + "), "|", FE))
  feols(f, data = s, cluster = ~GID_0, notes = FALSE)
}

# Clustered Wald test of R b = 0 (F with G - 1 denominator df, as in the handout).
wald <- function(b, V, R, G) {
  Rb <- R %*% b
  stat <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R), Rb)) / nrow(R)
  c(F = stat, df = nrow(R), p = pf(stat, nrow(R), G - 1, lower.tail = FALSE))
}
# Restrictions: each term equal across the panels of each group (consecutive
# differences within the group).
equality_R <- function(terms, groups, columns) {
  rows <- list()
  for (grp in groups) for (t in terms) for (i in seq_len(length(grp) - 1)) {
    r <- stats::setNames(numeric(length(columns)), columns)
    r[paste0(t, "__", grp[i])] <- 1
    r[paste0(t, "__", grp[i + 1])] <- -1
    rows[[length(rows) + 1]] <- r
  }
  do.call(rbind, rows)
}
GROUPS <- list(
  all = list(PANELS$id),
  gdp_source = lapply(names(CLIMATES), function(c) PANELS$id[PANELS$climate == c]),
  climate_product = lapply(GDPS, function(g) PANELS$id[PANELS$gdp == g]))

results <- list()
for (sample_name in names(SAMPLES)) for (spec in names(SPECS)) {
  message("Fitting: ", sample_name, " rows, ", spec, " arms")
  terms <- SPECS[[spec]]
  weather <- names(terms)
  arms <- grep("^[pn](TM|RR)2?(_X_r)?$", weather, value = TRUE)
  s <- stack_of(SAMPLES[[sample_name]], terms)
  G <- length(unique(s$GID_0))

  # Unrestricted: each panel its own model.
  m <- fit_stack(s, panel_terms(c(weather, LEVELS, NUMERIC_CONTROLS)))
  cols <- panel_terms(weather)
  b <- coef(m)[cols]
  V <- vcov(m)[cols, cols]
  # The same as the panels fitted one by one (coefficients; the SEs differ by
  # the small-sample factor only).
  one <- vapply(PANELS$id, function(id) {
    x <- s[s$panel == id, ]
    own <- paste0(c(LEVELS, NUMERIC_CONTROLS), "__", id)
    mm <- feols(as.formula(paste("g ~", paste(c(weather, own, factor_terms), collapse = " + "),
                                 "|", FE)), data = x, cluster = ~GID_0, notes = FALSE)
    max(abs(coef(mm)[weather] - b[paste0(weather, "__", id)]))
  }, 0)
  stopifnot(all(one < 1e-8))

  tests <- do.call(rbind, lapply(names(GROUPS), function(grp) {
    rbind(data.frame(group = grp, terms = "all weather terms",
                     t(wald(b, V, equality_R(weather, GROUPS[[grp]], cols), G))),
          data.frame(group = grp, terms = "arms only",
                     t(wald(b, V, equality_R(arms, GROUPS[[grp]], cols), G))))
  }))
  # Term by term (with many restrictions at once the clustered Wald test
  # over-rejects; these have 5 or fewer).
  by_term <- do.call(rbind, lapply(names(GROUPS), function(grp) do.call(rbind, lapply(weather, function(t) {
    data.frame(group = grp, term = t, t(wald(b, V, equality_R(t, GROUPS[[grp]], cols), G)))
  }))))

  # Pooled: H maps the common terms to the K copies.
  H <- do.call(rbind, replicate(K, diag(length(weather)), simplify = FALSE))
  combine <- function(W) {
    A <- solve(t(H) %*% W %*% H, t(H) %*% W)
    list(b = drop(A %*% b), V = A %*% V %*% t(A))
  }
  equal <- combine(diag(length(b)))
  efficient <- combine(solve(V))
  restricted <- fit_stack(s, c(weather, panel_terms(c(LEVELS, NUMERIC_CONTROLS))))
  # The implied weight of each panel in the efficient combination (of each term).
  A_eff <- solve(t(H) %*% solve(V) %*% H, t(H) %*% solve(V))
  weights_eff <- sapply(seq_len(K), function(k) diag(A_eff[, (k - 1) * length(weather) + seq_along(weather)]))
  dimnames(weights_eff) <- list(weather, PANELS$id)

  per_panel <- data.frame(term = rep(weather, K), panel = rep(PANELS$id, each = length(weather)),
                          estimate = unname(b), se = sqrt(diag(V)))
  pooled <- rbind(
    data.frame(term = weather, method = "equal-weight average", estimate = equal$b, se = sqrt(diag(equal$V))),
    data.frame(term = weather, method = "efficient minimum distance", estimate = efficient$b, se = sqrt(diag(efficient$V))),
    data.frame(term = weather, method = "restricted stacked", estimate = unname(coef(restricted)[weather]),
               se = unname(se(restricted)[weather])))
  rownames(pooled) <- NULL
  results[[paste(sample_name, spec)]] <- list(
    sample = sample_name, spec = spec, tests = tests, by_term = by_term, per_panel = per_panel, pooled = pooled,
    weights_efficient = weights_eff, coef = b, vcov = V,
    rows = table(s$panel)[PANELS$id], countries = G)
}

out_dir <- file.path("results", "dataset_pooling")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(results = results, panels = PANELS, xi = XI, common_rows = length(common_keys), built = Sys.time()),
        file.path(out_dir, "pooling.rds"))
for (r in results) {
  cat("\n====", r$sample, "rows,", r$spec, "arms; countries", r$countries, "\n")
  print(r$rows)
  print(r$tests, digits = 3)
  print(r$by_term, digits = 3)
  w <- reshape(r$per_panel[, c("term", "panel", "estimate")], idvar = "term", timevar = "panel", direction = "wide")
  print(w, digits = 3)
  print(r$pooled, digits = 3)
  print(round(r$weights_efficient, 2))
}
