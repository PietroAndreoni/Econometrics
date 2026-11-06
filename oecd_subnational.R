# The handout's linear-arm models on OECD regions: DOSE v2.14 regional GDP per
# capita (GADM1; constant 2015 local prices at the 2015 exchange rate, i.e.
# constant 2015 US$) x ERA5 (population weights), the OECD members DOSE covers.
# Saved to results/oecd_subnational/ (oecd_subnational.rds, .log).
#
# As in the handout (handout_preferred_specification.R), with the region as
# the unit: deviations from the 20-year trailing mean in absolute units,
# demeaned by the region's linear and quadratic trend, over the region's SD
# over the horizon (fixed SD); linear arms; moderators that move: the same
# 20-year trailing mean of T or P (centred at the median horizon climate of
# the regions) and the log mean regional GDP per capita of the 5 years before
# (less its mean over the rows). Year effects, region effects and region
# quadratic trends, the controls, the level terms of every moderator in every
# model; the non-weather growth outliers of the 5% tails (country-years)
# dropped from every region of the country. No income scaling. Errors
# clustered by region, with country clusters alongside. Each model also with
# country x year effects in place of the year effects (identified from the
# regions' deviations from their country's in each year; the national
# controls are absorbed), and without one country at a time.
#
#   A  arms only
#   B  arms x trailing climate (T with T's, P with P's)
#   C  arms x trailing climate and trailing income (T and P)
#
# Benchmark: the same models on the national panel with World Bank GDP
# (constant 2015 US$) x ERA5, the same OECD countries, clustered by country.
#
# Needs the knitr cache of weather_arms_robustness.Rmd and
# results/response_explorer/panels.rds. From the project root:
#   Rscript oecd_subnational.R

out_dir <- file.path("results", "oecd_subnational")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "oecd_subnational.log"), split = TRUE)

# The explorer's panel helpers (and the notebook's cache, with the DOSE panel).
bm <- readLines(file.path("apps", "response_explorer", "build_models.R"))
eval(parse(text = bm[seq_len(grep("^panels <- list\\(\\)", bm) - 1)]))
stopifnot(exists("dose"), dose$label == "DOSE x ERA5")

# OECD members (38).
OECD <- c("AUS", "AUT", "BEL", "CAN", "CHE", "CHL", "COL", "CRI", "CZE", "DEU", "DNK", "ESP", "EST", "FIN",
          "FRA", "GBR", "GRC", "HUN", "IRL", "ISL", "ISR", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX",
          "NLD", "NOR", "NZL", "POL", "PRT", "SVK", "SVN", "SWE", "TUR", "USA")

# The rows: the DOSE panel's sample without the outliers, OECD regions, with
# the deviations and the trailing income.
rows <- add_explorer_deviations(dose)$rows
rows <- rows[paste(rows$GID_1, rows$year) %in% paste(dose$clean$GID_1, dose$clean$year) & rows$GID_0 %in% OECD, ]
rows <- left_join(rows, trailing_income(dose, TRAILING_INCOME_YEARS), by = c("GID_1", "year"))
d <- as.data.frame(rows[is.finite(rows$inc_trailing), ])
cat("OECD members not in DOSE:", paste(setdiff(OECD, unique(dose$clean$GID_0)), collapse = ", "), "\n")
cat("Rows", nrow(d), "of", nrow(rows), "with the trailing income;", n_distinct(d$GID_1), "regions,",
    n_distinct(d$GID_0), "countries,", min(d$year), "-", max(d$year), "\n")

stopifnot(all(tapply(d$sd_TM_all, d$GID_1, function(x) diff(range(x))) == 0))
for (v in c("TM", "RR")) {
  x <- d[[paste0("dm_", v, "_lag_20_unstd_quadratic")]] / d[[paste0("sd_", v, "_all")]]
  d[[paste0("p", v)]] <- pmax(x, 0)
  d[[paste0("n", v)]] <- pmax(-x, 0)
}
unit_mean <- function(z) tapply(d[[z]], d$GID_1, function(x) x[[1]])
clim_median <- c(TM = median(unit_mean("mean_TM_all")), RR = median(unit_mean("mean_RR_all")))
d$clim_TM <- d$mean_TM_lag_20 - clim_median[["TM"]]
d$clim_RR <- d$mean_RR_lag_20 - clim_median[["RR"]]
d$y <- d$inc_trailing - mean(d$inc_trailing)
d$g <- d[[OUTCOME]]
cat(sprintf("Median horizon climate of the regions: T %.1f C, P %.2f m; mean trailing income %.0f 2015 US$\n",
            clim_median[["TM"]], clim_median[["RR"]], exp(mean(d$inc_trailing))))

LEVELS <- c("clim_TM", "clim_RR", "y")
CLIM <- c(TM = "clim_TM", RR = "clim_RR")
arm_terms <- function(v, modifiers = character(0)) {
  out <- list()
  for (a in paste0(c("p", "n"), v)) {
    out[[a]] <- a
    for (z in modifiers) out[[paste0(a, "_X_", z)]] <- c(a, z)
  }
  out
}
SPECS <- list(
  A = c(arm_terms("TM"), arm_terms("RR")),
  B = c(arm_terms("TM", CLIM[["TM"]]), arm_terms("RR", CLIM[["RR"]])),
  C = c(arm_terms("TM", c(CLIM[["TM"]], "y")), arm_terms("RR", c(CLIM[["RR"]], "y")))
)
LABELS <- c(A = "A: arms", B = "B: + climate", C = "C: + climate, income")
# Fixed effects: year effects (the handout's), or country x year effects
# (these absorb the national controls, which are then left out).
FIXED_EFFECTS <- c(year = "year + GID_1[t, t2]", country_year = "GID_0^year + GID_1[t, t2]")
fit_terms <- function(terms, data, cluster, fe = "year") {
  for (k in names(terms)) data[[k]] <- Reduce(`*`, lapply(terms[[k]], function(z) data[[z]]))
  controls <- if (fe == "year") CONTROLS else character(0)
  f <- as.formula(paste("g ~", paste(c(names(terms), LEVELS, controls), collapse = " + "), "|", FIXED_EFFECTS[[fe]]))
  feols(f, data = data, cluster = cluster, panel.id = ~GID_1 + year, notes = FALSE)
}

# Clustered Wald test of the terms = 0 with the covariance V and G clusters.
wald <- function(m, V, G, terms) {
  terms <- intersect(terms, names(coef(m)))
  b <- coef(m)[terms]
  stat <- as.numeric(t(b) %*% solve(V[terms, terms], b)) / length(terms)
  pf(stat, length(terms), G - 1, lower.tail = FALSE)
}
star <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
of <- function(v, z) paste0(c("p", "n"), v, "_X_", z)
TESTS <- list(`T x climate` = of("TM", "clim_TM"), `P x climate` = of("RR", "clim_RR"),
              `T x income` = of("TM", "y"), `P x income` = of("RR", "y"))
KEY <- c(pTM = "Warm arm", pTM_X_clim_TM = "Warm x climate", pTM_X_y = "Warm x income",
         nTM = "Cold arm", nTM_X_clim_TM = "Cold x climate", nTM_X_y = "Cold x income",
         pRR = "Wet arm", pRR_X_clim_RR = "Wet x climate", pRR_X_y = "Wet x income",
         nRR = "Dry arm", nRR_X_clim_RR = "Dry x climate", nRR_X_y = "Dry x income")

# The three models on data `data` with the clusterings `clusters` (the first
# is the fit's; named by label), as tables.
run <- function(data, clusters, fe = "year") {
  fits <- lapply(SPECS, fit_terms, data = data, cluster = clusters[[1]], fe = fe)
  out <- list(fits = fits, coef = list(), tests = list())
  for (cl in names(clusters)) {
    G <- length(unique(data[[all.vars(clusters[[cl]])]][obs(fits$A)]))
    cells <- sapply(fits, function(m) {
      V <- vcov(m, cluster = clusters[[cl]])
      b <- coef(m)
      s <- sqrt(diag(V))
      p <- 2 * pt(-abs(b / s), G - 1)
      ifelse(names(KEY) %in% names(b), sprintf("%+.3f (%.3f)%s", 100 * b[names(KEY)], 100 * s[names(KEY)], star(p[names(KEY)])), "")
    })
    rownames(cells) <- KEY
    out$coef[[cl]] <- cells
    added <- sapply(names(fits), function(k) {
      prev <- c(A = NA, B = "A", C = "B")[[k]]
      if (is.na(prev)) return(NA)
      wald(fits[[k]], vcov(fits[[k]], cluster = clusters[[cl]]), G, setdiff(names(SPECS[[k]]), names(SPECS[[prev]])))
    })
    tests <- rbind(sapply(fits, function(m) sapply(TESTS, function(t) {
      if (all(t %in% names(coef(m)))) sprintf("%.3f", wald(m, vcov(m, cluster = clusters[[cl]]), G, t)) else ""
    })), `Added terms (vs previous)` = ifelse(is.na(added), "", sprintf("%.3f", added)))
    out$tests[[cl]] <- tests
    out$G[[cl]] <- G
  }
  out$fit <- rbind(`Within R2` = sapply(fits, function(m) sprintf("%.4f", fitstat(m, "wr2")[[1]])),
                   AIC = sapply(fits, function(m) sprintf("%.0f", AIC(m))),
                   BIC = sapply(fits, function(m) sprintf("%.0f", BIC(m))),
                   Rows = sapply(fits, nobs))
  out
}
show <- function(res, title) {
  cat("\n\n=====", title, "=====\n")
  for (cl in names(res$coef)) {
    tab <- res$coef[[cl]]
    colnames(tab) <- LABELS
    cat(sprintf("\nCoefficients, pp of growth per SD of the deviation (x 100); SE clustered by %s (%d clusters)\n", cl, res$G[[cl]]))
    print(noquote(tab))
    tests <- res$tests[[cl]]
    colnames(tests) <- LABELS
    cat(sprintf("Joint Wald tests, p (clustered by %s)\n", cl))
    print(noquote(tests))
  }
  fit <- res$fit
  colnames(fit) <- LABELS
  print(noquote(fit))
}
options(width = 200)

sub <- run(d, list(region = ~GID_1, country = ~GID_0))
show(sub, "OECD regions: DOSE x ERA5, year effects")
cy <- run(d, list(region = ~GID_1, country = ~GID_0), fe = "country_year")
show(cy, "OECD regions: DOSE x ERA5, country x year effects")

# How much of the variation of the deviations is left within the country-year:
# the variance of the deviation net of country x year and region effects over
# that net of year and region effects.
for (v in c("TM", "RR")) {
  dd <- data.frame(x = d[[paste0("p", v)]] - d[[paste0("n", v)]], GID_0 = d$GID_0, GID_1 = d$GID_1, year = d$year)
  within <- var(resid(feols(x ~ 1 | GID_0^year + GID_1, dd, notes = FALSE)))
  cat(sprintf("%s deviation: variance within the country-year %.2f of that with year effects\n", v,
              within / var(resid(feols(x ~ 1 | year + GID_1, dd, notes = FALSE)))))
}

# The response of each arm (pp per SD) at the 10th, 50th and 90th percentile
# of the regions' trailing income and climate, from model C (country clusters).
qs <- c(p10 = 0.1, p50 = 0.5, p90 = 0.9)
effects_of <- function(m) {
  V <- vcov(m, cluster = ~GID_0)
  out <- list()
  for (v in c("TM", "RR")) for (a in paste0(c("p", "n"), v)) for (zi in names(qs)) for (zc in names(qs)) {
    yv <- quantile(d$y, qs[[zi]])
    cv <- quantile(d[[CLIM[[v]]]], qs[[zc]])
    w <- setNames(c(1, cv, yv), c(a, paste0(a, "_X_", CLIM[[v]]), paste0(a, "_X_y")))
    est <- sum(w * coef(m)[names(w)])
    se <- sqrt(as.numeric(t(w) %*% V[names(w), names(w)] %*% w))
    p <- 2 * pt(-abs(est / se), length(unique(d$GID_0)) - 1)
    out[[length(out) + 1]] <- data.frame(arm = KEY[[a]], income = zi, climate = zc,
                                         effect = sprintf("%+.2f (%.2f)%s", 100 * est, 100 * se, star(p)))
  }
  do.call(rbind, out)
}
cat(sprintf("\nIncome percentiles: %s 2015 US$; T climate: %s C; P climate: %s m\n",
            paste(round(exp(mean(d$inc_trailing) + quantile(d$y, qs))), collapse = " / "),
            paste(round(clim_median[["TM"]] + quantile(d$clim_TM, qs), 1), collapse = " / "),
            paste(round(clim_median[["RR"]] + quantile(d$clim_RR, qs), 2), collapse = " / ")))
eff <- list(year = effects_of(sub$fits$C), country_year = effects_of(cy$fits$C))
for (fe in names(eff)) {
  cat("\nModel C,", fe, "effects: response (pp per SD, country-clustered SE) by percentile of the trailing income (rows) and climate (columns)\n")
  for (a in unique(eff[[fe]]$arm)) {
    w <- reshape(eff[[fe]][eff[[fe]]$arm == a, -1], idvar = "income", timevar = "climate", direction = "wide")
    names(w) <- gsub("effect.", "climate ", names(w), fixed = TRUE)
    cat("\n", a, "\n")
    print(w, row.names = FALSE)
  }
}

# Model C without one country at a time (coefficient x 100, t with country
# clusters), with each kind of fixed effects.
LOO_KEYS <- c("pTM", "pTM_X_clim_TM", "pTM_X_y", "nTM", "nTM_X_clim_TM", "pRR", "pRR_X_clim_RR", "nRR", "nRR_X_clim_RR")
loo <- lapply(setNames(names(FIXED_EFFECTS), names(FIXED_EFFECTS)), function(fe) {
  t(sapply(c("none", sort(unique(d$GID_0))), function(cc) {
    m <- fit_terms(SPECS$C, d[d$GID_0 != cc, ], ~GID_0, fe = fe)
    sprintf("%+.2f (%+.1f)", 100 * coef(m)[LOO_KEYS], coef(m)[LOO_KEYS] / se(m)[LOO_KEYS])
  }))
})
for (fe in names(loo)) {
  colnames(loo[[fe]]) <- KEY[LOO_KEYS]
  cat("\nModel C,", fe, "effects, without one country (coefficient, t with country clusters)\n")
  print(noquote(loo[[fe]]))
}

# Benchmark: the national panel (WB x ERA5), the same OECD countries.
nat <- new.env()
Sys.setenv(HANDOUT_GDP = "WB", HANDOUT_ROWS = "own", HANDOUT_OUTLIERS = "5", HANDOUT_MIN_INCOME = "0")
main <- readLines("handout_preferred_specification.R")
eval(parse(text = main[seq_len(grep("models <- list(", main, fixed = TRUE)[1] - 1)]), envir = nat)
dn <- nat$d[nat$d$GID_0 %in% unique(d$GID_0), ]
country_climate <- function(z) tapply(dn[[z]], dn$GID_0, function(x) x[[1]])
dn$clim_TM <- dn$mean_TM_lag_20 - median(country_climate("mean_TM_all"))
dn$clim_RR <- dn$mean_RR_lag_20 - median(country_climate("mean_RR_all"))
dn$y <- dn$inc_trailing - mean(dn$inc_trailing)
stopifnot(identical(nat$CONTROLS, CONTROLS))
national <- run(dn, list(country = ~GID_0))
show(national, sprintf("National benchmark: WB x ERA5, the same %d OECD countries (%d-%d)",
                       n_distinct(dn$GID_0), min(dn$year), max(dn$year)))

saveRDS(list(sub = sub[c("coef", "tests", "fit", "G")], country_year = cy[c("coef", "tests", "fit", "G")],
             effects = eff, loo = loo, national = national[c("coef", "tests", "fit", "G")],
             coef_sub = lapply(sub$fits, coef), vcov_sub = lapply(sub$fits, vcov),
             coef_cy = lapply(cy$fits, coef), vcov_cy = lapply(cy$fits, vcov, cluster = ~GID_0),
             countries = sort(unique(d$GID_0)), regions = n_distinct(d$GID_1), clim_median = clim_median),
        file.path(out_dir, "oecd_subnational.rds"))
sink()
