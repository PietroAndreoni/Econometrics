# ==============================================================================
# National stochastic frontier: climate in the production frontier,
# weather anomalies in the inefficiency term.
#
# Replicates the main approach of
#   Tol, R.S.J. (2021), "The economic impact of weather and climate",
#   arXiv:2102.13110v3 [econ.GN].
#
# A Cobb-Douglas frontier, concentrated in per-worker terms:
#
#   ln y_ct = b1*ln k_ct + f(Tbar_ct, Rbar_ct) + mu_c + tau*t + v_ct - u_ct   (3)
#   u_ct ~ Exp( d0 + d1*g(z(T_ct)) + d2*g(z(R_ct)) )                          (4)
#
# where Tbar/Rbar are the NORMAL_WINDOW-year climate normals over the years
# preceding t, z(.) is the weather anomaly standardised by the same window's
# standard deviation, and g(.) = |.| in the base specification. Climate shifts
# potential output; weather shifts the output gap.
#
# Tol estimates (3)-(4) as a Greene (2005) "true fixed effects" model with the
# Stata sfmodel command. No CRAN package fits a panel SFA with heteroscedastic
# inefficiency, so (3)-(4) is estimated here with sfaR::sfacross, whose `uhet`
# argument parameterises ln(sigma_u^2) the same way sfmodel does; the country
# effects mu_c are handled by PANEL_EFFECTS below (dummies = true fixed
# effects, or Mundlak means = correlated random effects). As a genuinely
# random-effects cross-check the script also fits the Battese & Coelli (1995)
# panel model with frontier::sfa, in which the same weather drivers enter the
# mean of a truncated-normal u_ct with a country-level random component.
#
# Data. Economic variables come from Penn World Table 11.0, staged by
# prepare_econ_data.R with "PWT110" in SELECT_ECON_SOURCE. Climate comes from
# the country-level monthly parquets in econometrics/data. The climate defaults
# below -- University of Delaware temperature and precipitation weighted by
# year-2000 population -- are Tol's own data, and with LABOUR_UNITS = "worker"
# the panel reproduces his Table 1 descriptives closely (7972 observations on
# 162 countries against his 7753 on 160; mean |z(T)| 0.958 against his 0.961).
#
# Run from the repository root:
#   Rscript econometrics/national_stochastic_frontier.R
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(zoo)
})

# ------------------------------------------------------------------------------
# Selectors
# ------------------------------------------------------------------------------
SELECT_CLIMATE_SOURCE <- "dela"       # dela (Tol's source), cru, era
SELECT_WEIGHT         <- "pop_2000"   # pop_2000 (Tol), un_ (area), pop_2015, ...
SELECT_ECON_SOURCE    <- "PWT110"     # econ_source tag written by prepare_econ_data.R

NORMAL_WINDOW <- 30                   # years defining "climate" (Tol: 30)
SAMPLE_YEARS  <- c(1950, 2014)        # Tol's sample
MIN_OBS       <- 5                    # drop countries with fewer usable years

# Tol concentrates the Cobb-Douglas by dividing output and capital by the labour
# force, so his y and k are per worker. That needs PWT's `emp` series, which is
# missing for a number of country-years. "capita" divides by population instead:
# more coverage, the same units as the rest of this repository, and only a level
# shift in ln(y) and ln(k) as long as the employment rate is close to constant
# within country -- which the country fixed effects then absorb.
LABOUR_UNITS  <- "capita"             # capita | worker (Tol)
DEP_LABEL <- switch(LABOUR_UNITS, capita = "ln(output per capita)",
                    worker = "ln(output per worker)")

UDIST         <- "exponential"        # Tol's base; "hnormal" is his robustness check
PANEL_EFFECTS <- "tfe"                # tfe | mundlak | pooled
TREND         <- "linear"             # none | linear | quadratic | cubic
POOR_DEF      <- "wb_highincome"      # wb_highincome (Tol's base) | q25_1990
HOT_QUANTILE  <- 0.75                 # "hot" = country mean T above this quantile
CLUSTER_SE    <- TRUE                 # country-clustered sandwich standard errors
ORTHOGONALISE <- TRUE                 # fit on the Q factor of the frontier design

RUN_BASELINE   <- TRUE                # Table 2, columns (1)-(6)
RUN_ROBUSTNESS <- TRUE                # Tables 3, 4 and A2
RUN_BC95       <- TRUE                # Battese-Coelli (1995) random-effects check
RUN_ECM        <- TRUE                # Equations (5)-(6), Tables 5-6
RUN_IMPACT     <- TRUE                # Section 5, stylised warming scenarios
ECM_LONGRUN    <- "sfa"               # sfa (Tol) | ols -- cointegrating-vector estimator
WARMING        <- c(1, 2, 3)          # uniform warming scenarios, degrees C

# Each fitted object carries its own design matrix and per-observation gradients
# -- around 25 MB per specification, and the orthogonalised design is dense, so
# it does not compress. "slim" drops those before saving, which keeps every
# coefficient, covariance and diagnostic but means efficiencies() and marginal()
# can no longer be re-run against the saved object. Use "full" when you want to.
SAVE_MODELS    <- "slim"              # slim | full

OUT_DIR <- "econometrics/output"

# World Bank high-income economies, classification vintage FY2021 -- the one Tol
# refers to. His baseline poverty dummy is "not classified high income by the
# World Bank"; POOR_DEF = "q25_1990" switches to his alternative definition
# (bottom quartile of output per worker in 1990), which is derived from the data
# and needs no external list.
WB_HIGH_INCOME <- c(
  "AND", "ARE", "ARG", "ATG", "AUS", "AUT", "ABW", "BEL", "BHR", "BHS",
  "BMU", "BRB", "BRN", "CAN", "CHE", "CHI", "CHL", "CUW", "CYM", "CYP",
  "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "FRO", "GBR", "GIB",
  "GRC", "GRL", "GUM", "HKG", "HRV", "HUN", "IMN", "IRL", "ISL", "ISR",
  "ITA", "JPN", "KNA", "KOR", "KWT", "LIE", "LTU", "LUX", "LVA", "MAC",
  "MAF", "MCO", "MLT", "MNP", "NCL", "NLD", "NOR", "NRU", "NZL", "OMN",
  "PAN", "PLW", "POL", "PRI", "PRT", "PYF", "QAT", "ROU", "SAU", "SGP",
  "SMR", "SVK", "SVN", "SWE", "SXM", "SYC", "TCA", "TTO", "TWN", "URY",
  "USA", "VGB", "VIR"
)

# ------------------------------------------------------------------------------
# Packages that are only needed once we actually estimate
# ------------------------------------------------------------------------------
need_pkg <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s. Install it with install.packages('%s').",
                 pkg, what, pkg), call. = FALSE)
  }
}
need_pkg("sfaR", "the stochastic frontier estimates")
# sfaR 1.0.1 needs Rcpp >= 1.1.1. On Windows the update silently fails while any
# other R session holds Rcpp.dll open, and the only symptom is a namespace error
# here, so say what to do about it.
tryCatch(suppressPackageStartupMessages(library(sfaR)),
         error = function(e) {
           stop(conditionMessage(e),
                "\n\nIf this is a version conflict on Rcpp: close every other R and ",
                "RStudio session (they lock Rcpp.dll on Windows), then run\n",
                "  install.packages(c(\"Rcpp\", \"sfaR\"))", call. = FALSE)
         })
if (CLUSTER_SE) need_pkg("sandwich", "country-clustered standard errors")
if (RUN_BC95) { need_pkg("frontier", "the Battese-Coelli (1995) cross-check"); need_pkg("plm", "panel data for frontier::sfa") }
if (RUN_ECM)  need_pkg("fixest", "the error-correction model")

sanitize_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
join_tag <- function(x) paste(sanitize_tag(unique(as.character(x))), collapse = "-")

RUN_TAG <- paste(
  join_tag(SELECT_ECON_SOURCE), join_tag(SELECT_CLIMATE_SOURCE),
  join_tag(SELECT_WEIGHT), paste0("w", NORMAL_WINDOW), join_tag(LABOUR_UNITS),
  join_tag(UDIST), join_tag(PANEL_EFFECTS), join_tag(POOR_DEF),
  sep = "_"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Data
# ==============================================================================

# Monthly country parquets are wide (one column per ISO3). Temperature is
# averaged over the twelve months; precipitation is averaged too and converted
# from mm to cm per month, which is Tol's unit (his Table 1 reports a mean of
# 9.375 cm/month).
aggregate_gadm0_monthly <- function(path, statistic = c("tmp", "pre")) {
  statistic <- match.arg(statistic)
  arrow::read_parquet(path) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(-Date, names_to = "GID_0", values_to = "value") %>%
    filter(!is.na(value), grepl("^[A-Z0-9]{3}$", GID_0)) %>%
    mutate(year = as.integer(substr(gsub("^X", "", as.character(Date)), 1, 4))) %>%
    group_by(GID_0, year) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    mutate(value = if (statistic == "pre") value / 10 else value)
}

climate_path <- function(statistic) {
  f <- file.path("econometrics", "data",
                 sprintf("gadm0_%s_%s_%s_monthly.parquet",
                         SELECT_CLIMATE_SOURCE, statistic, SELECT_WEIGHT))
  if (!file.exists(f)) {
    stop("Missing climate file: ", f, "\nAvailable: ",
         paste(basename(list.files("econometrics/data", "^gadm0_.*_monthly\\.parquet$")),
               collapse = ", "), call. = FALSE)
  }
  f
}

message("Reading climate (", SELECT_CLIMATE_SOURCE, ", ", SELECT_WEIGHT, ") ...")
climate <- inner_join(
  aggregate_gadm0_monthly(climate_path("tmp"), "tmp") %>% rename(TM = value),
  aggregate_gadm0_monthly(climate_path("pre"), "pre") %>% rename(RR = value),
  by = c("GID_0", "year")
) %>%
  group_by(GID_0) %>%
  # Countries outside the grid come through as all-zero columns.
  filter(any(TM != 0 | RR != 0)) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    # "the average temperature c.q. precipitation in country c in the thirty
    # years preceding year t": right-aligned rolling window, then lagged once so
    # that year t itself is excluded from its own normal.
    TMbar = lag(zoo::rollmeanr(TM, NORMAL_WINDOW, fill = NA)),
    RRbar = lag(zoo::rollmeanr(RR, NORMAL_WINDOW, fill = NA)),
    sdTM  = lag(zoo::rollapplyr(TM, NORMAL_WINDOW, sd, fill = NA)),
    sdRR  = lag(zoo::rollapplyr(RR, NORMAL_WINDOW, sd, fill = NA))
  ) %>%
  ungroup() %>%
  mutate(zTM = (TM - TMbar) / sdTM,
         zRR = (RR - RRbar) / sdRR)

# ---- Penn World Table -------------------------------------------------------
# prepare_econ_data.R stores rgdpna/pop, rnna/pop and emp/pop. Those are already
# the per-capita series; dividing the first two by the third instead recovers
# rgdpna/emp and rnna/emp, i.e. Tol's per-worker terms. See LABOUR_UNITS.
message("Reading Penn World Table (", SELECT_ECON_SOURCE, ") ...")
econ_files <- list.files("econometrics/data", pattern = "^econ_processed_.*\\.parquet$",
                         full.names = TRUE)
read_pwt <- function(files) {
  for (f in rev(files)) {
    d <- arrow::read_parquet(f) %>% tibble::as_tibble()
    if (!all(c("k_pc_usd", "share_emp_pop") %in% names(d))) next
    d <- d %>% filter(econ_source == SELECT_ECON_SOURCE)
    if (nrow(d) > 0) {
      message("  using ", basename(f))
      return(d)
    }
  }
  stop("No econ_processed_*.parquet holds '", SELECT_ECON_SOURCE,
       "' rows with k_pc_usd/share_emp_pop.\n",
       "Run prepare_econ_data.R with '", SELECT_ECON_SOURCE,
       "' in SELECT_ECON_SOURCE first.", call. = FALSE)
}
econ <- read_pwt(econ_files) %>%
  mutate(per_head = switch(LABOUR_UNITS,
                           capita = 1,
                           worker = share_emp_pop,
                           stop("unknown LABOUR_UNITS: ", LABOUR_UNITS))) %>%
  transmute(GID_0,
            year,
            y = grp_pc_usd / per_head,
            k = k_pc_usd   / per_head) %>%
  filter(is.finite(y), is.finite(k), y > 0, k > 0)

# ---- merge and derive -------------------------------------------------------
panel <- inner_join(econ, climate, by = c("GID_0", "year")) %>%
  filter(year >= SAMPLE_YEARS[1], year <= SAMPLE_YEARS[2],
         !is.na(TMbar), !is.na(RRbar), !is.na(zTM), !is.na(zRR)) %>%
  group_by(GID_0) %>% filter(n() >= MIN_OBS) %>% ungroup()

poor_dummy <- function(data, definition) {
  switch(
    definition,
    wb_highincome = data %>% distinct(GID_0) %>%
      mutate(P = as.integer(!GID_0 %in% WB_HIGH_INCOME)),
    q25_1990 = {
      ref <- data %>% filter(year == 1990)
      cut <- stats::quantile(ref$y, 0.25, na.rm = TRUE)
      ref %>% transmute(GID_0, P = as.integer(y < cut))
    },
    stop("unknown poverty definition: ", definition)
  )
}

hot_lookup <- panel %>%
  group_by(GID_0) %>% summarise(mean_TM = mean(TM), .groups = "drop")
hot_cut <- stats::quantile(hot_lookup$mean_TM, HOT_QUANTILE)
hot_lookup <- hot_lookup %>% mutate(H = as.integer(mean_TM > hot_cut))

panel <- panel %>%
  inner_join(poor_dummy(panel, POOR_DEF), by = "GID_0") %>%
  inner_join(hot_lookup, by = "GID_0") %>%
  mutate(
    lny = log(y),
    lnk = log(k),
    trend   = year - min(year),
    trend_2 = trend^2,
    trend_3 = trend^3,
    TMbar_2 = TMbar^2,
    RRbar_2 = RRbar^2,
    TMbar_RRbar = TMbar * RRbar,
    cty = factor(GID_0),
    yr  = factor(year),
    # anomaly transforms used by the alternative specifications
    absTM = abs(zTM),     absRR = abs(zRR),      # g(.) = |.|   (base)
    sqTM  = zTM^2,        sqRR  = zRR^2,         # g(.) = (.)^2
    linTM = zTM,          linRR = zRR,           # g(.) = .
    posTM = pmax(zTM, 0), negTM = pmax(-zTM, 0), # asymmetric shocks
    posRR = pmax(zRR, 0), negRR = pmax(-zRR, 0)
  ) %>%
  group_by(GID_0) %>%
  mutate(m_lnk = mean(lnk), m_TMbar = mean(TMbar), m_RRbar = mean(RRbar)) %>%
  ungroup()

message(sprintf("Estimation panel: %d observations, %d countries, %d-%d",
                nrow(panel), dplyr::n_distinct(panel$GID_0),
                min(panel$year), max(panel$year)))

# ---- Table 1: descriptive statistics ----------------------------------------
descriptives <- panel %>%
  transmute(`ln(y)` = lny, `ln(k)` = lnk, `T` = TMbar, `R` = RRbar,
            `|z(T)|` = absTM, `|z(R)|` = absRR, `P` = P, `H` = H) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(mean = mean(value), sd = sd(value),
            min = min(value), max = max(value), n = n(), .groups = "drop")
cat("\n== Table 1: descriptive statistics ==\n")
print(as.data.frame(descriptives), digits = 4, row.names = FALSE)

# ==============================================================================
# 2. Estimation machinery
# ==============================================================================

trend_terms <- function(kind = TREND) {
  switch(kind,
         none      = character(0),
         linear    = "trend",
         quadratic = c("trend", "trend_2"),
         cubic     = c("trend", "trend_2", "trend_3"),
         stop("unknown TREND: ", kind))
}

effect_terms <- function(kind = PANEL_EFFECTS) {
  switch(kind,
         tfe     = "cty",
         mundlak = c("m_lnk", "m_TMbar", "m_RRbar"),
         pooled  = character(0),
         stop("unknown PANEL_EFFECTS: ", kind))
}

# Frontier building blocks. Column (6) of Tol's Table 2 -- the preferred
# specification -- is the full quadratic in the climate normals plus their
# interaction, with poverty interacted only with the rainfall terms.
FR_CORE     <- c("lnk", "TMbar", "TMbar_2", "RRbar", "RRbar_2")
FR_POOR_ALL <- c("P:TMbar", "P:TMbar_2", "P:RRbar", "P:RRbar_2")
FR_POOR_R   <- c("P:RRbar", "P:RRbar_2")
FR_TXR      <- "TMbar_RRbar"
FR_TXR_POOR <- "P:TMbar_RRbar"
FR_PREF     <- c(FR_CORE, FR_TXR, FR_POOR_R, FR_TXR_POOR)

# Inefficiency drivers for an arbitrary pair of anomaly transforms.
ineff_terms <- function(tvar, rvar, poor = TRUE, hot = TRUE) {
  c(tvar, rvar,
    if (poor) paste0("P:", c(tvar, rvar)),
    if (hot)  paste0("H:", c(tvar, rvar)))
}
INEFF_PREF <- ineff_terms("absTM", "absRR")

# One fit. Rows with missing values on any model variable are dropped up front
# so that the clustering vector stays aligned with the estimation sample.
#
# In raw units the frontier design is badly conditioned: the climate quadratics
# reach 840 while the country dummies are 0/1, and with 160 dummies the
# information matrix comes back numerically singular (reciprocal condition
# number around 1e-21), so the standard errors are rounding noise however tight
# they look. ORTHOGONALISE fits on the Q factor of the QR decomposition of the
# same design instead. That is an exact reparameterisation -- X = QR, so
# b = R^-1 b_tilde and V = R^-1 V_tilde R^-T -- and it reaches the identical
# log-likelihood in a fifth of the iterations with a gradient norm five orders
# of magnitude smaller. The coefficients reported downstream are always in raw
# units, so the printed tables stay comparable with Tol's.
fit_sfa <- function(label, frontier, uhet, udist = UDIST, trend = TREND,
                    effects = PANEL_EFFECTS, data = panel) {
  rhs <- c(frontier, trend_terms(trend), effect_terms(effects))
  f_frontier <- as.formula(paste("lny ~", paste(rhs, collapse = " + ")))
  f_uhet <- if (length(uhet) > 0) as.formula(paste("~", paste(uhet, collapse = " + ")))
  used <- unique(unlist(lapply(c(rhs, uhet), function(s) all.vars(str2lang(s)))))
  d <- data[stats::complete.cases(data[, c("lny", used), drop = FALSE]), , drop = FALSE]

  args <- list(udist = udist, S = 1L, method = "bfgs", logDepVar = TRUE)
  if (!is.null(f_uhet)) args$uhet <- f_uhet

  if (ORTHOGONALISE) {
    X <- stats::model.matrix(f_frontier, d)
    qrX <- qr(X)
    if (qrX$rank < ncol(X)) {
      message("    ", label, ": dropping ", ncol(X) - qrX$rank, " aliased column(s)")
      X <- X[, sort(qrX$pivot[seq_len(qrX$rank)]), drop = FALSE]
      qrX <- qr(X)
    }
    Q <- qr.Q(qrX)
    Rf <- qr.R(qrX)
    qnames <- sprintf(".q%04d", seq_len(ncol(Q)))
    dq <- as.data.frame(Q)
    names(dq) <- qnames
    dq$lny <- d$lny
    for (v in setdiff(if (is.null(f_uhet)) character(0) else all.vars(f_uhet),
                      names(dq))) {
      dq[[v]] <- d[[v]]
    }
    args$formula <- as.formula(paste("lny ~ -1 +", paste(qnames, collapse = " + ")))
    args$data <- dq
  } else {
    args$formula <- f_frontier
    args$data <- d
  }

  t0 <- Sys.time()
  m <- tryCatch(suppressWarnings(do.call(sfaR::sfacross, args)),
                error = function(e) {
                  warning("fit failed for '", label, "': ", conditionMessage(e),
                          call. = FALSE)
                  NULL
                })
  if (is.null(m)) return(NULL)

  # Transform from the fitted parameterisation back to raw units. The frontier
  # betas come first in the sfacross parameter vector, followed by the Zu and Zv
  # blocks, which are never reparameterised.
  b <- coef(m)
  if (ORTHOGONALISE) {
    nq <- ncol(Q)
    stopifnot(m$nXvar == nq)
    Tmat <- diag(length(b))
    Tmat[seq_len(nq), seq_len(nq)] <- backsolve(Rf, diag(nq))
    b <- as.numeric(Tmat %*% b)
    names(b) <- c(colnames(X), names(coef(m))[-seq_len(nq)])
    dimnames(Tmat) <- list(names(b), names(coef(m)))
  } else {
    Tmat <- diag(length(b))
    dimnames(Tmat) <- list(names(b), names(b))
  }

  # as.formula() captures the calling frame, which at this point holds the
  # design matrix and its QR factors -- around 100 MB per fit that would
  # otherwise be carried through the session and serialised. do.call() likewise
  # inlines the entire data frame into the fitted object's call. Neither is
  # needed downstream: every formula below is only ever evaluated against a
  # data frame that supplies all of its variables.
  environment(f_frontier) <- globalenv()
  environment(m$formula) <- globalenv()
  m$call$data <- as.name("data")
  # do.call also stores the formula and uhet objects themselves in the call, and
  # each carries its own copy of that same environment reference.
  for (nm in c("formula", "uhet")) {
    if (inherits(m$call[[nm]], "formula")) environment(m$call[[nm]]) <- globalenv()
  }

  message(sprintf("  %-22s logL = %10.1f   n = %5d   %4.0f s", label,
                  as.numeric(logLik(m)), nrow(d),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(label = label, fit = m, data = d, formula = f_frontier,
       frontier = frontier, uhet = uhet, udist = udist,
       coef_raw = b, tmat = Tmat)
}

# Country-clustered sandwich covariance, mapped back to raw units; sfaR supplies
# bread() and estfun() methods, so vcovCL applies directly to the fitted object.
vcov_sfa <- function(res) {
  V <- NULL
  if (CLUSTER_SE) {
    V <- tryCatch(sandwich::vcovCL(res$fit, cluster = res$data$GID_0),
                  error = function(e) NULL)
    if (is.null(V)) {
      warning("clustered vcov failed for '", res$label, "'; using the ML Hessian",
              call. = FALSE)
    }
  }
  if (is.null(V)) V <- vcov(res$fit)
  V <- res$tmat %*% V %*% t(res$tmat)
  dimnames(V) <- list(names(res$coef_raw), names(res$coef_raw))
  V
}

# R orders the components of an interaction by where each variable first appears
# in the formula, so the same term comes back as "P:TMbar_RRbar" in one
# specification and "TMbar_RRbar:P" in another. Reordering the components -- the
# regressor first, the poverty and heat dummies last -- keeps the rows of the
# printed tables lined up across columns.
canonical_term <- function(x) {
  prefix <- sub("^((Zu_|Zv_)?).*$", "\\1", x)
  body <- sub("^(Zu_|Zv_)", "", x)
  body <- vapply(
    strsplit(body, ":", fixed = TRUE),
    function(p) paste(c(sort(setdiff(p, c("P", "H"))),
                        intersect(c("P", "H"), p)), collapse = ":"),
    character(1)
  )
  paste0(prefix, body)
}

tidy_sfa <- function(res) {
  if (is.null(res)) return(NULL)
  b <- res$coef_raw
  se <- sqrt(diag(vcov_sfa(res)))[names(b)]
  tibble(
    spec = res$label, term = canonical_term(names(b)), term_raw = names(b),
    estimate = as.numeric(b),
    std_error = as.numeric(se), z = as.numeric(b) / as.numeric(se)
  ) %>%
    mutate(p_value = 2 * stats::pnorm(-abs(z)),
           block = case_when(grepl("^Zu_", term) ~ "inefficiency",
                             grepl("^Zv_", term) ~ "noise",
                             TRUE ~ "frontier")) %>%
    filter(!grepl("^cty|^yr", term))
}

stars <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", ""))))

print_specs <- function(coefs, title) {
  if (is.null(coefs) || nrow(coefs) == 0) return(invisible(NULL))
  wide <- coefs %>%
    mutate(cell = sprintf("%.4g%s (%.2f)", estimate, stars(p_value), z)) %>%
    select(block, term, spec, cell) %>%
    pivot_wider(names_from = spec, values_from = cell, values_fill = "") %>%
    arrange(factor(block, levels = c("frontier", "inefficiency", "noise")))
  cat("\n== ", title, " ==\n", sep = "")
  print(as.data.frame(wide), row.names = FALSE, right = FALSE)
  cat("estimate<stars> (z);  * p<0.05, ** p<0.01, *** p<0.001",
      if (CLUSTER_SE) "; country-clustered SEs" else "", "\n", sep = "")
}

fit_stats <- function(results) {
  bind_rows(lapply(Filter(Negate(is.null), results), function(r) {
    # efficiencies() returns u (the conditional inefficiency, Jondrow et al.
    # 1982) alongside teBC, the Battese-Coelli (1988) conditional efficiency.
    # It does not cover every udist/uhet combination -- notably a frontier with
    # no inefficiency drivers -- so fall back to NA rather than aborting.
    e <- tryCatch(sfaR::efficiencies(r$fit), error = function(err) NULL)
    tibble(spec = r$label, udist = r$udist, logLik = as.numeric(logLik(r$fit)),
           n_par = length(r$coef_raw), n_obs = nrow(r$data),
           n_countries = dplyr::n_distinct(r$data$GID_0),
           mean_u = if (is.null(e)) NA_real_ else mean(e$u, na.rm = TRUE),
           mean_efficiency = if (is.null(e)) NA_real_ else mean(e$teBC, na.rm = TRUE),
           # Diagnostics worth reading before trusting a column: a gradient norm
           # that is not close to zero means the optimiser stopped short, and a
           # tiny reciprocal condition number means the information matrix is
           # near-singular, so the standard errors on that column are not
           # interpretable however tight they look.
           grad_norm = r$fit$gradientNorm,
           rcond_hessian = tryCatch(rcond(r$fit$invHessian),
                                    error = function(err) NA_real_),
           # maxLik pads its return message, so trim before comparing.
           opt_status = trimws(as.character(r$fit$optStatus)))
  }))
}

results <- list()

# ==============================================================================
# 3. Table 2 -- baseline specifications
# ==============================================================================
if (RUN_BASELINE) {
  cat("\n-- Table 2: baseline --\n")
  baseline_specs <- list(
    `(1)` = list(fr = FR_CORE,
                 uh = ineff_terms("absTM", "absRR", poor = FALSE, hot = FALSE)),
    `(2)` = list(fr = c(FR_CORE, FR_POOR_ALL),
                 uh = ineff_terms("absTM", "absRR", poor = FALSE, hot = FALSE)),
    `(3)` = list(fr = c(FR_CORE, FR_POOR_ALL),
                 uh = ineff_terms("absTM", "absRR", hot = FALSE)),
    `(4)` = list(fr = c(FR_CORE, FR_POOR_ALL),
                 uh = INEFF_PREF),
    `(5)` = list(fr = c(FR_CORE, FR_POOR_ALL, FR_TXR, FR_TXR_POOR),
                 uh = INEFF_PREF),
    `(6)` = list(fr = FR_PREF, uh = INEFF_PREF)
  )
  for (nm in names(baseline_specs)) {
    s <- baseline_specs[[nm]]
    results[[nm]] <- fit_sfa(nm, s$fr, s$uh)
  }
  print_specs(bind_rows(lapply(names(baseline_specs), function(nm) tidy_sfa(results[[nm]]))),
              paste0("Table 2: baseline results, dependent variable ", DEP_LABEL))
}

# The preferred specification, re-used by everything below.
preferred <- results[["(6)"]]
if (is.null(preferred)) {
  cat("\n-- preferred specification (column 6) --\n")
  preferred <- fit_sfa("preferred", FR_PREF, INEFF_PREF)
  results[["preferred"]] <- preferred
}

# ==============================================================================
# 4. Tables 3, 4 and A2 -- robustness
# ==============================================================================
if (RUN_ROBUSTNESS && !is.null(preferred)) {
  cat("\n-- Tables 3/4/A2: robustness --\n")

  alt_poor <- if (POOR_DEF == "wb_highincome") "q25_1990" else "wb_highincome"
  panel_altpoor <- panel %>%
    select(-P) %>%
    inner_join(poor_dummy(panel, alt_poor), by = "GID_0")

  rob <- list()
  rob[["base"]]        <- preferred
  rob[["alt_poor"]]    <- fit_sfa(paste0("poor=", alt_poor), FR_PREF, INEFF_PREF,
                                  data = panel_altpoor)
  rob[["sq_anom"]]     <- fit_sfa("sq_anomalies", FR_PREF, ineff_terms("sqTM", "sqRR"))
  rob[["lin_anom"]]    <- fit_sfa("lin_anomalies", FR_PREF, ineff_terms("linTM", "linRR"))
  rob[["asym_anom"]]   <- fit_sfa("asym_anomalies", FR_PREF,
                                  c(ineff_terms("posTM", "posRR"),
                                    ineff_terms("negTM", "negRR")))
  # Weather moved out of inefficiency and into the frontier -- the specification
  # most of the literature uses, and the one Tol rejects.
  rob[["wx_frontier"]] <- fit_sfa("weather_in_frontier",
                                  c(FR_PREF, INEFF_PREF), character(0))
  rob[["hnormal"]]     <- fit_sfa("half_normal", FR_PREF, INEFF_PREF, udist = "hnormal")
  # Capital as a substitute for climate (Table 4).
  rob[["k_climate"]]   <- fit_sfa("capital_x_climate",
                                  c(FR_PREF, "lnk:TMbar", "lnk:TMbar_2",
                                    "lnk:RRbar", "lnk:RRbar_2"), INEFF_PREF)
  # Trend variants (Table A2).
  rob[["no_trend"]]    <- fit_sfa("no_trend", FR_PREF, INEFF_PREF, trend = "none")
  rob[["quad_trend"]]  <- fit_sfa("quadratic_trend", FR_PREF, INEFF_PREF,
                                  trend = "quadratic")

  results <- c(results, rob[setdiff(names(rob), "base")])
  print_specs(bind_rows(lapply(rob, tidy_sfa)),
              paste0("Tables 3/4/A2: robustness, dependent variable ", DEP_LABEL))
}

# ==============================================================================
# 5. Marginal effects and efficiency
# ==============================================================================
marginal_effects <- NULL
efficiency_scores <- NULL
if (!is.null(preferred)) {
  me <- tryCatch(sfaR::marginal(preferred$fit), error = function(e) NULL)
  if (!is.null(me)) {
    marginal_effects <- tibble(driver = names(me),
                               mean_effect = sapply(me, mean),
                               sd_effect = sapply(me, sd))
    cat("\n== Marginal effects of the weather drivers ==\n")
    cat("(Eu_ = on E[u], Vu_ = on V[u]; u is the shortfall from the frontier in log points)\n")
    print(as.data.frame(marginal_effects), digits = 4, row.names = FALSE)
  }
  eff <- tryCatch(sfaR::efficiencies(preferred$fit), error = function(e) NULL)
  if (!is.null(eff)) {
    efficiency_scores <- preferred$data %>%
      select(GID_0, year, P, H, TMbar, RRbar, zTM, zRR) %>%
      mutate(u = eff$u, efficiency = eff$teBC)
    cat("\n== Technical efficiency, preferred specification ==\n")
    print(efficiency_scores %>%
            group_by(poor = P, hot = H) %>%
            summarise(mean_u = mean(u), mean_efficiency = mean(efficiency),
                      n = n(), .groups = "drop") %>%
            as.data.frame(), digits = 4, row.names = FALSE)
  }
}

# ==============================================================================
# 6. Battese-Coelli (1995): a genuinely random panel inefficiency term
# ==============================================================================
# sfacross treats every country-year as an independent draw. frontier::sfa on a
# panel instead gives u_ct a country-level random component and lets the same
# weather drivers shift the mean of a truncated normal. Country heterogeneity in
# the frontier is picked up by Mundlak means rather than dummies, which the
# estimator cannot carry.
bc95 <- NULL
if (RUN_BC95) {
  cat("\n-- Battese-Coelli (1995) random-effects cross-check --\n")
  fr_bc <- paste(c(FR_PREF, trend_terms(), effect_terms("mundlak")), collapse = " + ")
  uh_bc <- paste(INEFF_PREF, collapse = " + ")
  f_bc <- as.formula(paste("lny ~", fr_bc, "|", uh_bc))
  d_bc <- panel %>%
    select(GID_0, year, lny, lnk, TMbar, TMbar_2, RRbar, RRbar_2, TMbar_RRbar,
           P, H, absTM, absRR, trend, trend_2, trend_3, m_lnk, m_TMbar, m_RRbar) %>%
    stats::na.omit() %>%
    as.data.frame()
  bc95 <- tryCatch(
    suppressWarnings(frontier::sfa(f_bc,
                                   data = plm::pdata.frame(d_bc, c("GID_0", "year")),
                                   ineffDecrease = TRUE, truncNorm = TRUE)),
    error = function(e) {
      warning("Battese-Coelli fit failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  if (!is.null(bc95)) {
    cat("\n== Battese-Coelli (1995), truncated normal, country random effects ==\n")
    print(summary(bc95))
    cat("Note: with ineffDecrease = TRUE the z-coefficients are signed so that a",
        "positive value means LESS inefficiency -- the opposite of the sfacross",
        "Zu_ coefficients above. The weather drivers enter the MEAN of a",
        "truncated normal here, which is far less well identified than the",
        "exponential scale used above: expect large, offsetting Z_ estimates",
        "with wide standard errors. Read this as a sign check, not as a",
        "second set of point estimates.\n", sep = "\n")
  }
}

# ==============================================================================
# 7. Error-correction model (Equations 5-6)
# ==============================================================================
ecm_short <- NULL
V <- NULL
if (RUN_ECM) {
  cat("\n-- Error-correction model --\n")
  # Long run: the cointegrating vector, with year dummies standing in for the
  # linear trend. V is the deviation from potential output.
  if (ECM_LONGRUN == "sfa") {
    lr <- fit_sfa("cointegrating_vector", c(FR_PREF, "yr"), INEFF_PREF, trend = "none")
    if (!is.null(lr)) {
      results[["cointegrating_vector"]] <- lr
      V <- lr$data %>% mutate(V = as.numeric(residuals(lr$fit)))
      print_specs(tidy_sfa(lr), "Table 5: cointegrating vector")
    }
  } else {
    f_lr <- as.formula(paste("lny ~", paste(FR_PREF, collapse = " + "),
                             "| GID_0 + year"))
    lr <- fixest::feols(f_lr, panel)
    V <- panel[fixest::obs(lr), ] %>% mutate(V = as.numeric(resid(lr)))
    cat("\n== Table 5: cointegrating vector (OLS, country + year FE) ==\n")
    print(summary(lr))
  }

  if (!is.null(V)) {
    short <- V %>%
      arrange(GID_0, year) %>%
      group_by(GID_0) %>%
      mutate(dlny = lny - lag(lny), lag_V = lag(V), gap = year - lag(year)) %>%
      ungroup() %>%
      filter(gap == 1, !is.na(dlny), !is.na(lag_V))
    ecm_short <- fixest::feols(
      dlny ~ zTM + zRR + P:zTM + P:zRR + H:zTM + H:zRR + lag_V | GID_0 + year,
      data = short, cluster = ~GID_0)
    cat("\n== Table 6: short-run error correction ==\n")
    print(summary(ecm_short))
    cat("The coefficient on lag_V is the speed of adjustment back to potential",
        "output (Tol reports 0.06).\n", sep = "\n")
  }
}

# ==============================================================================
# 8. Section 5 -- stylised warming scenarios
# ==============================================================================
impact <- NULL
if (RUN_IMPACT && !is.null(preferred)) {
  # Climate channel: shift the 30-year normal by dT and rebuild the frontier
  # design matrix, so that whatever terms the preferred specification happens to
  # carry are all differenced consistently.
  frontier_shift <- function(res, dT) {
    d0 <- res$data
    d1 <- d0 %>%
      mutate(TMbar = TMbar + dT,
             TMbar_2 = TMbar^2,
             TMbar_RRbar = TMbar * RRbar)
    X0 <- stats::model.matrix(res$formula, d0)
    X1 <- stats::model.matrix(res$formula, d1)
    b <- res$coef_raw
    keep <- intersect(colnames(X0), names(b))
    as.numeric((X1[, keep, drop = FALSE] - X0[, keep, drop = FALSE]) %*% b[keep])
  }

  key <- paste(preferred$data$GID_0, preferred$data$year)
  last_year <- preferred$data %>%
    group_by(GID_0) %>% filter(year == max(year)) %>% ungroup()
  idx <- match(paste(last_year$GID_0, last_year$year), key)

  impact <- bind_rows(lapply(WARMING, function(dT) {
    shift <- frontier_shift(preferred, dT)
    last_year %>%
      transmute(GID_0, year, P, H, TMbar, RRbar,
                warming = dT, dlny_frontier = shift[idx])
  }))

  # Weather channel: under the exponential, E[u] = exp(0.5 * Zu'delta), so the
  # semi-elasticity of expected inefficiency to a one-standard-deviation larger
  # absolute anomaly is half the relevant sum of Zu_ coefficients.
  b <- preferred$coef_raw
  zu <- b[grepl("^Zu_", names(b))]
  zu_sum <- function(pat) {
    v <- zu[grepl(pat, names(zu))]
    if (length(v)) sum(v) else 0
  }
  weather_semi <- impact %>%
    distinct(GID_0, P, H) %>%
    mutate(
      dlnEu_dabsT = 0.5 * (zu_sum("^Zu_absTM$") +
                             P * zu_sum("absTM:P$|^Zu_P:absTM$") +
                             H * zu_sum("absTM:H$|^Zu_H:absTM$")),
      dlnEu_dabsR = 0.5 * (zu_sum("^Zu_absRR$") +
                             P * zu_sum("absRR:P$|^Zu_P:absRR$") +
                             H * zu_sum("absRR:H$|^Zu_H:absRR$"))
    )
  impact <- impact %>% left_join(weather_semi, by = c("GID_0", "P", "H"))

  cat("\n== Section 5: impact of uniform warming on potential output ==\n")
  cat("dlny_frontier = change in potential", DEP_LABEL, "\n")
  cat("dlnEu_dabs*   = semi-elasticity of expected inefficiency to a one-sd anomaly.\n")
  print(impact %>%
          group_by(warming, poor = P, hot = H) %>%
          summarise(mean_dlny = mean(dlny_frontier),
                    median_dlny = median(dlny_frontier),
                    mean_dlnEu_T = mean(dlnEu_dabsT),
                    n = n(), .groups = "drop") %>%
          as.data.frame(), digits = 3, row.names = FALSE)
}

# ==============================================================================
# 9. Save
# ==============================================================================
coef_table <- bind_rows(lapply(results, tidy_sfa))
stat_table <- fit_stats(results)

coef_csv <- file.path(OUT_DIR, paste0("sfa_coefficients_", RUN_TAG, ".csv"))
stat_csv <- file.path(OUT_DIR, paste0("sfa_fitstats_", RUN_TAG, ".csv"))
eff_csv  <- file.path(OUT_DIR, paste0("sfa_efficiency_", RUN_TAG, ".csv"))
imp_csv  <- file.path(OUT_DIR, paste0("sfa_impact_", RUN_TAG, ".csv"))
mod_rds  <- file.path(OUT_DIR, paste0("sfa_models_", RUN_TAG, ".rds"))

write.csv(coef_table, coef_csv, row.names = FALSE)
write.csv(stat_table, stat_csv, row.names = FALSE)
if (!is.null(efficiency_scores)) write.csv(efficiency_scores, eff_csv, row.names = FALSE)
if (!is.null(impact)) write.csv(impact, imp_csv, row.names = FALSE)
slim_result <- function(r) {
  if (is.null(r) || SAVE_MODELS == "full") return(r)
  r$vcov <- vcov_sfa(r)
  r$data <- NULL
  r$tmat <- NULL
  r$fit$dataTable <- NULL
  r$fit$gradL_OBS <- NULL
  r
}
saved_results <- lapply(results, slim_result)

saveRDS(list(results = saved_results, bc95 = bc95, ecm_short = ecm_short,
             descriptives = descriptives, marginal_effects = marginal_effects,
             fit_stats = stat_table,
             settings = list(climate_source = SELECT_CLIMATE_SOURCE,
                             weight = SELECT_WEIGHT, econ_source = SELECT_ECON_SOURCE,
                             normal_window = NORMAL_WINDOW,
                             sample_years = SAMPLE_YEARS,
                             labour_units = LABOUR_UNITS, udist = UDIST,
                             panel_effects = PANEL_EFFECTS, trend = TREND,
                             poor_def = POOR_DEF, orthogonalise = ORTHOGONALISE,
                             cluster_se = CLUSTER_SE, save_models = SAVE_MODELS)),
        mod_rds)

cat("\n== Fit statistics ==\n")
print(as.data.frame(stat_table), digits = 5, row.names = FALSE)

suspect <- stat_table %>%
  filter(!is.finite(grad_norm) | grad_norm > 1 |
           is.na(rcond_hessian) | rcond_hessian < 1e-10 |
           opt_status != "successful convergence")
if (nrow(suspect) > 0) {
  cat("\nCaution -- these columns did not produce a usable information matrix,",
      "so read their point estimates only as a sign check and ignore their",
      "standard errors entirely:\n", sep = "\n")
  print(as.data.frame(suspect %>% select(spec, udist, grad_norm, rcond_hessian,
                                         opt_status)),
        digits = 4, row.names = FALSE)
}
cat("\nWrote:\n", coef_csv, "\n", stat_csv, "\n",
    if (!is.null(efficiency_scores)) paste0(eff_csv, "\n") else "",
    if (!is.null(impact)) paste0(imp_csv, "\n") else "",
    mod_rds, "\n", sep = "")
