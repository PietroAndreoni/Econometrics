# Fits a model of the explorer on demand from the estimation panels of
# build_models.R (`panels`: the content of panels.rds), and keeps it for the
# session. The model is that of "Income or the agricultural share" in
# weather_arms_robustness.Rmd: growth on the arms of the demeaned deviations of
# T and P, each arm interacted with the chosen moderators, with the controls,
# year and unit effects (and unit trends), errors clustered by country.
#
# Fixed effects and demeaning go together: the deviations were demeaned by the
# unit mean ("unit"), mean and linear trend ("linear") or mean and linear and
# quadratic trend ("quadratic", the notebook's) before their arms are cut, and
# the regression has the same unit slopes.
#
# Moderators: "fixed", the horizon mean of the climate, long-run income,
# long-run share and (PWT) long-run capital per capita (constant within a unit,
# so they have no level terms of their own); "trailing", the climate of the
# deviation's own baseline (its trailing 20 or 30 year mean, the horizon mean
# with the horizon baseline), the log of the mean income and (PWT) of the mean
# capital per capita of the panels$trailing_income_years years before, and
# with `share_trailing` (WB GDP only) the mean share of the
# panels$trailing_share_years years before. Time-varying moderators enter with
# their level terms too (for those in at least one interaction): without them
# an interaction of a non-negative arm stands in for the omitted level. The
# levels are controls: the app's responses and projections use the deviation
# terms only. The trailing models run on the rows with all the trailing values
# they use (the full window); a model with capital (or a PWT outcome other than
# GDP per capita) on the rows with it.
#
# Outcome, in log points: the growth of GDP per capita ("log", the notebook's);
# with PWT also the growth of the capital stock per capita ("capital", rnna /
# pop) and of GDP per worker ("worker", rgdpna / emp).

library(fixest)

ARM_SUFFIXES <- list("Linear arms" = c("_pos", "_neg"), "Quadratic arms" = c("_pos2", "_neg2"))
FE_CHOICES <- c("Year and unit effects" = "unit", "+ unit linear trends" = "linear",
                "+ unit linear and quadratic trends" = "quadratic")
FE_SHORT <- c(unit = "no unit trends", linear = "unit linear trends",
              quadratic = "unit linear and quadratic trends")
SPLINE_COLUMNS <- c("clim_TM", "clim_RR", "inc", "agr", "kap")
OUTCOME_CHOICES <- c("Growth of GDP per capita (Δ log)" = "log",
                     "Growth of the capital stock per capita (Δ log, PWT)" = "capital",
                     "Growth of GDP per worker (Δ log, PWT)" = "worker")
OUTCOME_SHORT <- c(log = "GDP per capita", capital = "capital per capita", worker = "GDP per worker")
# The outcomes and the moderator only PWT has.
PWT_OUTCOMES <- c("capital", "worker")
PWT_MODIFIERS <- "kap"
# How income scales the deviation terms of T and of P, each on its own:
# "standard", each arm with its interactions; "dollars", every
# deviation term (the arm and its interactions) times idol_v = (Ybar / Y_i)^(1 -
# xi_v) (Ybar: the geometric mean income of the rows; Y_i the income of the
# interactions, long-run or trailing), so that the dollar loss per person is
# proportional to Y_i^xi_v beyond the ticked interactions (xi = 0: the same
# dollars at every income; xi = 1: the same %, the standard terms); "both", the
# terms and their copies times idol_v (the encompassing model). In the
# coefficients the scaled terms end in ":idol" (idol_TM or idol_RR by arm).
INCOME_FORMS <- c("Standard (income through its interaction only)" = "standard",
                  "Both: every deviation term also times (Ybar / Y)^(1 - xi)" = "both",
                  "Scaled: every deviation term times (Ybar / Y)^(1 - xi)" = "dollars")
XI_RANGE <- c(0, 2)
# The income scaling of a fitted model: form and xi of T and P.
model_income_forms <- function(m) {
  o <- model_options(list(outcome = m$outcome, income_form = m$income_form, xi = m$xi))
  list(form = o$income_form, xi = o$xi)
}
# The short label of the income scaling of a model (NULL: standard for both).
income_form_label <- function(form, xi) {
  part <- function(v) {
    f <- form[[v]]
    if (f == "standard") return("unscaled")
    scaled <- if (xi[[v]] == 0) "equal dollars" else sprintf("scaled, xi = %g", xi[[v]])
    if (f == "both") paste("+", sub("^equal dollars$", "equal-dollar terms", sub("^scaled", "scaled terms", scaled)))
    else scaled
  }
  if (all(form == "standard")) return(NULL)
  if (form[["TM"]] == form[["RR"]] && (form[["TM"]] == "standard" || xi[["TM"]] == xi[["RR"]])) return(part("TM"))
  paste0("T ", part("TM"), ", P ", part("RR"))
}

fixed_effects_formula <- function(fe, slopes) {
  s <- slopes[[fe]]
  paste("year +", if (length(s)) paste0("GID_1[", paste(s, collapse = ", "), "]") else "GID_1")
}

# Differences in fit against noise, as in weather_arms_robustness.Rmd ("Model
# terms and effects"): Vuong test clustered by country for non-nested models,
# the clustered Wald test of the added terms for nested ones (5%).
NOISE_LEVEL <- 0.05
CRITERIA <- c(aic = "AIC", bic = "BIC", r2 = "within adj. R2")
fit_record <- function(m) {
  e <- stats::resid(m)
  s2 <- mean(e^2)
  list(aic = stats::AIC(m), bic = stats::BIC(m), war2 = fixest::fitstat(m, "war2")[[1]],
       n = stats::nobs(m), ll = -0.5 * log(2 * pi * s2) - e^2 / (2 * s2), country = m$country)
}
criterion_value <- function(rec, criterion) {
  switch(criterion, aic = -rec$aic / 2, bic = -rec$bic / 2, r2 = -(rec$n / 2) * log(1 - rec$war2))
}
fit_verdict <- function(a, b, criterion, nested_p = NULL) {
  stopifnot(a$n == b$n, identical(a$country, b$country))
  difference <- criterion_value(a, criterion) - criterion_value(b, criterion)
  beyond_noise <- if (!is.null(nested_p)) {
    !is.na(nested_p) && nested_p < NOISE_LEVEL
  } else {
    D <- tapply(a$ll - b$ll, a$country, sum)
    G <- length(D)
    se <- sqrt(G / (G - 1) * sum((D - mean(D))^2))
    is.finite(se) && se > 0 && abs(difference / se) >= stats::qnorm(1 - NOISE_LEVEL / 2)
  }
  if (!beyond_noise) "indifferent" else if (difference > 0) "better" else "worse"
}
fit_verdicts <- function(a, b, nested_p = NULL) {
  a <- fit_record(a)
  b <- fit_record(b)
  paste(vapply(names(CRITERIA), function(cr) fit_verdict(a, b, cr, nested_p), ""), collapse = " / ")
}
added_terms_p <- function(larger, smaller) {
  added <- setdiff(names(stats::coef(larger)), names(stats::coef(smaller)))
  stopifnot(length(added) > 0, all(names(stats::coef(smaller)) %in% names(stats::coef(larger))))
  b <- stats::coef(larger)[added]
  V <- stats::vcov(larger)[added, added, drop = FALSE]
  stats::pchisq(as.numeric(t(b) %*% solve(V, b)), df = length(added), lower.tail = FALSE)
}

# The options of a model (o: gdp, climate, weight, units, baseline, form, inter
# = list(TM, RR), shape, fe, moderators, share_trailing, outcome), completed:
# the trailing share only with WB GDP and trailing moderators; the PWT outcomes
# and the capital interaction only with PWT GDP.
model_options <- function(o) {
  o$outcome <- o$outcome %||% "log"
  o$fe <- o$fe %||% "quadratic"
  o$moderators <- o$moderators %||% "fixed"
  o$share_trailing <- isTRUE(o$share_trailing) && o$moderators == "trailing" && identical(o$gdp, "WB")
  if (!is.null(o$gdp) && o$gdp != "PWT") {
    if (o$outcome %in% PWT_OUTCOMES) o$outcome <- "log"
    if (!is.null(o$inter)) o$inter <- lapply(o$inter, function(z) setdiff(z, PWT_MODIFIERS))
  }
  # The income scaling of T and of P (a single value: both), and its xi (in
  # XI_RANGE). Scaled terms with xi = 1 are the standard ones.
  per_variable <- function(x, default) {
    x <- unlist(x %||% default)
    if (length(x) == 1) x <- c(TM = x, RR = x)
    c(TM = unname(x[["TM"]]), RR = unname(x[["RR"]]))
  }
  form <- per_variable(o$income_form, "standard")
  xi <- pmin(pmax(as.numeric(per_variable(o$xi, 0)), XI_RANGE[1]), XI_RANGE[2])
  names(xi) <- c("TM", "RR")
  form[xi == 1] <- "standard"
  xi[form == "standard"] <- 1
  o$income_form <- form
  o$xi <- xi
  o
}
# The key of a pair of subsets of interactions, e.g. "T: clim+inc, P: clim+agr"
# ("none": no interaction).
interaction_key <- function(inter) {
  part <- function(x) if (length(x)) paste(x, collapse = "+") else "none"
  paste0("T: ", part(inter$TM), ", P: ", part(inter$RR))
}
model_key <- function(o, panels) {
  moderators <- if (o$share_trailing) "trailing + share" else o$moderators
  paste(o$gdp, o$climate, o$baseline, o$weight, o$units, o$form, interaction_key(o$inter), o$shape,
        o$fe, moderators, o$outcome,
        sprintf("income T: %s %g, P: %s %g", o$income_form[["TM"]], o$xi[["TM"]], o$income_form[["RR"]], o$xi[["RR"]]),
        sep = " | ")
}

# The verdicts of model a against model b by AIC / BIC / within adjusted R2
# (their fit records), when they run on the same rows and outcome (same panel,
# row numbers and outcome); NULL otherwise.
record_verdicts <- function(a, b) {
  if (!identical(a$panel_rows, b$panel_rows)) return(NULL)
  paste(vapply(names(CRITERIA), function(cr) fit_verdict(a$record, b$record, cr), ""), collapse = " / ")
}

# Joint Wald tests (fixest::wald, clustered) of the deviation terms of model m
# with interactions `inter`: for T and for P, all their terms; the terms
# interacted with each characteristic (both arms, with their spline terms); with
# three interactions each pair of them; with two or more, all the interactions;
# with both income forms, the plain terms (zero: equal dollars) and the terms
# times Ybar / Y (zero: equal %).
MODIFIER_TEST_LABELS <- c(clim = "climate", inc = "income", agr = "share", kap = "capital")
# The column of a deviation term ("devTM_pos:clim_TM:idol") in the fit, and back.
safe_term <- function(x) gsub(":", "_X_", x, fixed = TRUE)
term_name <- function(x) gsub("_X_", ":", x, fixed = TRUE)
joint_tests <- function(m, inter, xi = c(TM = 0, RR = 0)) {
  terms <- term_name(grep("^dev", names(stats::coef(m)), value = TRUE))
  sets <- list()
  for (v in c("TM", "RR")) {
    short <- if (v == "TM") "T" else "P"
    own <- terms[substr(terms, 4, 5) == v]
    if (!length(own)) next
    sets[[paste0(short, ": all terms")]] <- own
    scaled <- grepl(":idol$", own)
    if (any(scaled) && !all(scaled)) {
      zero <- if (xi[[v]] == 0) "equal dollars" else sprintf("dollar loss ~ Y^%g", xi[[v]])
      sets[[paste0(short, ": plain terms (zero: ", zero, ")")]] <- own[!scaled]
      sets[[paste0(short, ": terms times (Ybar / Y)^", 1 - xi[[v]], " (zero: equal %)")]] <- own[scaled]
    }
    of <- function(z) {
      column <- if (z == "clim") paste0("clim_", v) else z
      own[grepl(paste0(":", column, "(_nl[0-9]+)?(:idol)?$"), own)]
    }
    mods <- inter[[v]]
    for (z in mods) sets[[paste0(short, ": ", MODIFIER_TEST_LABELS[[z]], " interactions")]] <- of(z)
    if (length(mods) >= 3) {
      for (pair in utils::combn(mods, 2, simplify = FALSE)) {
        sets[[paste0(short, ": ", paste(MODIFIER_TEST_LABELS[pair], collapse = " + "), " interactions")]] <-
          unlist(lapply(pair, of))
      }
    }
    if (length(mods) >= 2) sets[[paste0(short, ": all interactions")]] <- unlist(lapply(mods, of))
  }
  sets <- Filter(length, sets)
  do.call(rbind, lapply(names(sets), function(k) {
    w <- fixest::wald(m, keep = paste0("^", safe_term(sets[[k]]), "$"), print = FALSE)
    data.frame(test = k, stat = w$stat, p = w$p, df1 = w$df1, stringsAsFactors = FALSE)
  }))
}

# The rows of a model: the panel's, those with the trailing values it uses
# (.row: the row number in the panel).
model_rows <- function(o, panels) {
  d <- panels$panels[[paste(o$gdp, o$climate, o$weight, sep = " | ")]]$rows
  d$.row <- seq_len(nrow(d))
  if (o$moderators == "trailing") d <- d[is.finite(d$inc_trailing), ]
  if (o$share_trailing) d <- d[is.finite(d$agr_trailing), ]
  if (o$outcome != "log") d <- d[is.finite(d[[panels$outcome_columns[[o$outcome]]]]), ]
  if ("kap" %in% unlist(o$inter)) {
    d <- d[is.finite(d[[if (o$moderators == "trailing") "kap_trailing" else "kap_level"]]), ]
  }
  d
}

model_cache <- new.env()
explorer_model <- function(o, panels) {
  o <- model_options(o)
  key <- model_key(o, panels)
  if (!is.null(model_cache[[key]])) return(model_cache[[key]])
  panel <- panels$panels[[paste(o$gdp, o$climate, o$weight, sep = " | ")]]
  window <- panels$baselines[[o$baseline]]
  trailing <- o$moderators == "trailing"
  d <- model_rows(o, panels)

  # The outcome (see the header).
  d$.outcome <- as.numeric(d[[if (o$outcome == "log") panels$outcome else panels$outcome_columns[[o$outcome]]]])

  # The arms of the deviations demeaned as the fixed effects.
  for (v in c("TM", "RR")) {
    dev <- d[[paste0("dm_", v, "_", window, panels$units_suffix[[o$units]], "_", o$fe)]]
    for (s in ARM_SUFFIXES[[o$form]]) d[[paste0("dev", v, s)]] <- arm_value(s, dev)
  }
  # The moderators. The climate is centred at the median of the long-run
  # climate of the panel's units (fixed or trailing), so that the arms' own
  # coefficients are the response at the median climate and at the mean income
  # and share of the rows (the cap of the app is relative to it).
  units <- panels$units[[paste(o$gdp, o$climate, o$baseline, o$weight, o$units, sep = " | ")]]
  clim_median <- c(TM = stats::median(units$mean_T), RR = stats::median(units$mean_P))
  clim_window <- if (trailing) window else "all"
  for (v in c("TM", "RR")) {
    d[[paste0("clim_", v)]] <- d[[paste0("mean_", v, "_", clim_window)]] - clim_median[[v]]
  }
  inc_level <- if (trailing) d$inc_trailing else d$inc_level
  agr_level <- if (o$share_trailing) d$agr_trailing else d$agr_level
  d$inc <- inc_level - mean(inc_level)
  d$agr <- agr_level - mean(agr_level)
  # Log capital per capita less its mean over the rows (PWT; only when it is in
  # an interaction, the rows have it).
  uses_kap <- "kap" %in% unlist(o$inter)
  kap_level <- if (!uses_kap) 0 else if (trailing) d$kap_trailing else d$kap_level
  d$kap <- kap_level - mean(kap_level)
  # (Ybar / Y_i)^(1 - xi_v), Ybar the geometric mean income of the rows.
  for (v in c("TM", "RR")) d[[paste0("idol_", v)]] <- exp(-(1 - o$xi[[v]]) * d$inc)
  # The non-linear spline terms <z>_nl<j> (knots at the spline percentiles of
  # the units' mean values), less their value at zero.
  spline_columns <- setdiff(SPLINE_COLUMNS, if (!uses_kap) "kap")
  unit_values <- stats::aggregate(d[spline_columns], list(GID_1 = d$GID_1), mean)
  knots <- lapply(stats::setNames(spline_columns, spline_columns), function(z) {
    stats::quantile(unit_values[[z]], panels$spline_knot_probs, names = FALSE)
  })
  for (z in spline_columns) {
    nl <- rcs_nonlinear(d[[z]], knots[[z]])
    at_zero <- rcs_nonlinear(0, knots[[z]])
    for (j in seq_len(ncol(nl))) d[[paste0(z, "_nl", j)]] <- nl[, j] - at_zero[1, j]
  }

  # The time-varying moderators in at least one interaction (their levels).
  varying <- c(if (trailing && window != "all" && "clim" %in% o$inter$TM) "clim_TM",
               if (trailing && window != "all" && "clim" %in% o$inter$RR) "clim_RR",
               if (trailing && "inc" %in% unlist(o$inter)) "inc",
               if (trailing && uses_kap) "kap",
               if (o$share_trailing && "agr" %in% unlist(o$inter)) "agr")
  # (Ybar / Y)^(1 - xi) varies with the trailing income: its level too (once
  # when T and P are scaled by the same power).
  scaled_v <- names(o$income_form)[o$income_form != "standard"]
  varying_idol <- if (trailing && length(scaled_v)) paste0("idol_", scaled_v[!duplicated(o$xi[scaled_v])])
  with_splines <- function(z, shape) {
    if (shape == "splines") c(z, unlist(lapply(z, function(x) {
      paste0(x, "_nl", seq_len(length(panels$spline_knot_probs) - 2))
    }))) else z
  }
  terms_of <- function(shape) {
    deviation <- unlist(lapply(c("TM", "RR"), function(v) {
      arms <- paste0("dev", v, ARM_SUFFIXES[[o$form]])
      modifiers <- with_splines(sub("^clim$", paste0("clim_", v), o$inter[[v]]), shape)
      plain <- c(arms, unlist(lapply(modifiers, function(z) paste0(arms, ":", z))))
      switch(o$income_form[[v]], standard = plain, dollars = paste0(plain, ":idol"),
             both = c(plain, paste0(plain, ":idol")))
    }))
    c(deviation, with_splines(varying, shape), varying_idol, panels$controls)
  }
  # Each deviation term enters as a column of its own (the product of its
  # parts), named safe_term(): in a formula, R orders the variables of an
  # interaction by their first appearance, which would scramble the names. The
  # part "idol" is idol_TM or idol_RR, by the term's variable.
  fit <- function(shape) {
    terms <- terms_of(shape)
    for (t in grep("^dev", terms, value = TRUE)) {
      d[[safe_term(t)]] <- Reduce(`*`, lapply(strsplit(t, ":", fixed = TRUE)[[1]], function(x) {
        d[[if (x == "idol") paste0("idol_", substr(t, 4, 5)) else x]]
      }))
    }
    terms[startsWith(terms, "dev")] <- safe_term(terms[startsWith(terms, "dev")])
    f <- stats::as.formula(paste(".outcome ~", paste(terms, collapse = " + "), "|",
                                 fixed_effects_formula(o$fe, panels$fe_slopes)))
    m <- fixest::feols(f, data = d, cluster = ~GID_0, panel.id = ~GID_1 + year, notes = FALSE)
    m$country <- d$GID_0[fixest::obs(m)]
    m
  }
  m <- fit(o$shape)
  # The splines against the linear interactions (same rows): the Wald test of
  # the non-linear terms and the verdict by AIC / BIC / within adjusted R2.
  spline_test <- NULL
  if (o$shape == "splines" && any(grepl("_nl", names(stats::coef(m))))) {
    linear <- fit("linear")
    nl_p <- added_terms_p(m, linear)
    spline_test <- list(nonlinear_p = nl_p, verdict = fit_verdicts(m, linear, nested_p = nl_p),
                        aic = stats::AIC(m) - stats::AIC(linear), bic = stats::BIC(m) - stats::BIC(linear))
  }
  ct <- fixest::coeftable(m)
  rownames(ct) <- term_name(rownames(ct))
  keep <- grep("^dev", rownames(ct), value = TRUE)
  level_terms <- intersect(c(with_splines(varying, o$shape), varying_idol), rownames(ct))
  coefs <- stats::setNames(stats::coef(m), term_name(names(stats::coef(m))))
  V <- stats::vcov(m)
  dimnames(V) <- list(term_name(rownames(V)), term_name(colnames(V)))
  model <- list(
    key = key, gdp = o$gdp, climate = o$climate, weight = o$weight, units = o$units,
    baseline = o$baseline, window = window, form = o$form, inter = o$inter, shape = o$shape,
    fe = o$fe, moderators = o$moderators, share_trailing = o$share_trailing,
    outcome = o$outcome, income_form = o$income_form, xi = o$xi,
    coef = coefs[keep], vcov = V[keep, keep, drop = FALSE],
    # The coefficients' p-values (clustered), the fit, the joint tests, and
    # what the noise tests of differences in fit need (its rows of the panel).
    p = ct[keep, 4],
    fit = c(war2 = unname(fixest::fitstat(m, "war2")[[1]]), aic = stats::AIC(m), bic = stats::BIC(m)),
    tests = joint_tests(m, o$inter, o$xi),
    record = fit_record(m),
    panel_rows = c(paste(o$gdp, o$climate, o$weight, o$outcome), d$.row[fixest::obs(m)]),
    levels = data.frame(term = level_terms, estimate = ct[level_terms, 1], se = ct[level_terms, 2],
                        p = ct[level_terms, 4], row.names = NULL),
    nobs = stats::nobs(m), countries = length(unique(m$country)),
    clim_median = clim_median, inc_mean = mean(inc_level), agr_mean = mean(agr_level),
    kap_mean = if (uses_kap) mean(kap_level),
    horizon = panel$horizon,
    # The windows of the trailing income and share (NULL: fixed), for the
    # projections that let them follow the SSP (capital has no SSP path: it
    # stays at the country's value).
    inc_window = if (trailing) panels$trailing_income_years,
    share_window = if (o$share_trailing) panels$trailing_share_years,
    knots = if (o$shape == "splines") knots, spline_test = spline_test
  )
  model_cache[[key]] <- model
  model
}

# The rows of the fixed and of the trailing models of a panel, and with the
# trailing share (for the warning on its sample).
trailing_rows <- function(gdp, climate, weight, panels) {
  d <- panels$panels[[paste(gdp, climate, weight, sep = " | ")]]$rows
  count <- function(x) c(rows = sum(x), countries = length(unique(d$GID_0[x])))
  rbind(fixed = count(rep(TRUE, nrow(d))), trailing = count(is.finite(d$inc_trailing)),
        share = count(is.finite(d$inc_trailing) & is.finite(d$agr_trailing)))
}
