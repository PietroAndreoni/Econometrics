# Adaptation arms versus short-term shocks, PWT x UDelaware.
#
# Two channels through which an unusual year can affect growth: the distance
# from the climate the economy is adapted to (adaptation arms, entering
# linearly as a warm and a cold arm, wet and dry for precipitation, each
# interacted with the long-run climate), and surprises against short-term
# expectations (shocks, entering squared). Two designs:
#
#   overlapping  arms on the anomaly against the trailing 30-year mean scaled
#                by the SD over the entire climate series (abs_zTMp, abs_zTMm);
#                shocks on the deviation from the Hamilton trend scaled by the
#                SD of the Hamilton cycle (hTM^2). The two measures share most
#                of their variation.
#   decomposed   the same 30-year anomaly split exactly into two additive parts,
#                  zTM = (TM - mean_TM_lag_30) / sd_TM_all = gapTM + cycTM,
#                  gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
#                  cycTM = (TM - trend_TM) / sd_TM_all,
#                arms on the slow adaptation gap (gapTMp = max(gapTM, 0),
#                gapTMm = max(-gapTM, 0)) and shocks on the squared cycle
#                (cycTM^2); the same for precipitation.
#
# The Hamilton trend is the OLS forecast of each unit's series from its values
# 2 to 5 years earlier, with coefficients estimated over the whole series (see
# hamilton_trend() in functions/climate_transforms.R).
#
# For each design, three models on the identical sample: arms only, shocks
# only, and both. Wald tests of each block given the other, and fit (within
# R2, AIC, BIC). Arm effects per SD at long-run climates of 5, 15 and 25
# degrees Celsius (0.5, 1.5 and 2.5 m of precipitation). Fixed effects as in
# main_analysis.Rmd (or without the unit trends, run option "notrends"),
# errors clustered by country, all available years, area and population
# weights.
#
# Output: results/adaptation_shock_pwt_udel/ (..._notrends/ without trends)

source("load_functions.R")

ECON <- "PWT"
CLIMATE_SOURCE <- "UDelaware"
CLIMATE_WEIGHTS <- c(area = "area", pop = "concurrent population")
ARM_BASELINE <- "lag_30"
ARM_SD <- "all"
# Run option: "trends" (default) keeps the unit-specific linear trends of
# main_analysis.Rmd; "notrends" drops them, so slow climate gaps are not
# absorbed by the trends. Rscript test_adaptation_shock_model.R notrends
TRENDS <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(TRENDS)) TRENDS <- "trends"
stopifnot(TRENDS %in% c("trends", "notrends"))
FIXED_EFFECTS <- if (TRENDS == "trends") "year + GID_1 + GID_1[year]" else "year + GID_1"
AT_VALUES <- list(mean_TM_all = c(5, 15, 25), mean_RR_all = c(0.5, 1.5, 2.5))
ARM_LABELS <- c("warm arm", "cold arm", "wet arm", "dry arm")
# Functional form of the arms: linear (distance from the adapted climate) or
# quadratic (squared distance).
ARM_FORMS <- c("linear", "quadratic")

# Arms (with the long-run climate each is interacted with) and shocks of
# each design, in the order of ARM_LABELS.
DESIGNS <- list(
  overlapping = list(
    arms = c(abs_zTMp = "mean_TM_all", abs_zTMm = "mean_TM_all",
             abs_zRRp = "mean_RR_all", abs_zRRm = "mean_RR_all"),
    shocks = c("hTM", "hRR")
  ),
  decomposed = list(
    arms = c(gapTMp = "mean_TM_all", gapTMm = "mean_TM_all",
             gapRRp = "mean_RR_all", gapRRm = "mean_RR_all"),
    shocks = c("cycTM", "cycRR")
  )
)

OUTPUT_DIR <- file.path(
  "results",
  if (TRENDS == "trends") "adaptation_shock_pwt_udel" else "adaptation_shock_pwt_udel_notrends"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit <- function(terms, data) {
  fit_panel_model(terms, data, fixed_effects = FIXED_EFFECTS, cluster = ~GID_0)
}
# Coefficient on the interaction of a and b, in whichever order it is named.
interaction_name <- function(a, b, terms) {
  intersect(c(paste0(a, ":", b), paste0(b, ":", a)), terms)[[1]]
}
block_p <- function(model, names) {
  fixest::wald(model, keep = paste0("^(", paste(gsub("([()^:])", "\\\\\\1", names),
                                               collapse = "|"), ")$"),
               print = FALSE)$p
}
# Variation net of the fixed effects, for within-unit correlations.
within <- function(x, data) {
  stats::resid(fixest::feols(as.formula(paste(x, "~ 1 |", FIXED_EFFECTS)), data,
                             panel.id = c("GID_1", "year")))
}

effects <- list()
tests <- list()
models <- list()
for (w in names(CLIMATE_WEIGHTS)) {
  cfg <- panel_config(
    econ_year_min = min(load_econ_panel(resolve_econ_source(ECON))$year),
    climate_source = CLIMATE_SOURCE,
    climate_weight = CLIMATE_WEIGHTS[[w]],
    clim_history_years = 1L
  )
  panel <- suppressMessages(build_dat(ECON, config = cfg, baseline = ARM_BASELINE,
                                      sd_baseline = ARM_SD))
  hamilton <- suppressMessages(build_dat(ECON, config = cfg, baseline = "trend")) %>%
    select(GID_1, year, hTM = zTM, hRR = zRR)
  d <- panel %>%
    left_join(hamilton, by = c("GID_1", "year")) %>%
    mutate(
      gapTM = (trend_TM - mean_TM_lag_30) / sd_TM_all,
      cycTM = (TM - trend_TM) / sd_TM_all,
      gapRR = (trend_RR - mean_RR_lag_30) / sd_RR_all,
      cycRR = (RR - trend_RR) / sd_RR_all,
      gapTMp = pmax(gapTM, 0), gapTMm = pmax(-gapTM, 0),
      gapRRp = pmax(gapRR, 0), gapRRm = pmax(-gapRR, 0)
    )
  # The decomposition is exact: gap + cycle is the 30-year anomaly.
  ok <- stats::complete.cases(d[, c("zTM", "gapTM", "cycTM", "zRR", "gapRR", "cycRR")])
  stopifnot(isTRUE(all.equal(d$gapTM[ok] + d$cycTM[ok], d$zTM[ok])),
            isTRUE(all.equal(d$gapRR[ok] + d$cycRR[ok], d$zRR[ok])))

  for (design in names(DESIGNS)) for (form in ARM_FORMS) {
    arms <- DESIGNS[[design]]$arms
    shocks <- DESIGNS[[design]]$shocks
    # Linear arms enter as abs_zTMp, quadratic arms as abs_zTMp^2 (coefficient
    # I(abs_zTMp^2)), each with its climate interaction.
    arm_vars <- if (form == "quadratic") paste0(names(arms), "^2") else names(arms)
    arm_coefs <- if (form == "quadratic") paste0("I(", arm_vars, ")") else arm_vars
    arm_terms <- paste(c(arm_vars, paste0(arm_vars, ":", arms)), collapse = " + ")
    shock_terms <- paste(paste0(shocks, "^2"), collapse = " + ")

    both <- fit(paste(arm_terms, "+", shock_terms), d)
    # The other models on the identical sample: the outcome is withheld elsewhere.
    same <- d
    same$dlgrp_pc_usd[setdiff(seq_len(nrow(d)), fixest::obs(both))] <- NA
    arms_only <- fit(arm_terms, same)
    shocks_only <- fit(shock_terms, same)
    stopifnot(nobs(arms_only) == nobs(both), nobs(shocks_only) == nobs(both))
    models[[w]][[paste(design, form)]] <- list(arms_only = arms_only,
                                               shocks_only = shocks_only, both = both)

    terms <- names(stats::coef(both))
    shock_names <- paste0("I(", shocks, "^2)")
    arm_names <- setdiff(terms, shock_names)
    used <- d[fixest::obs(both), ]
    arm_measure <- if (design == "decomposed") "gapTM" else "zTM"
    shock_measure <- if (design == "decomposed") "cycTM" else "hTM"
    tests[[length(tests) + 1]] <- tibble::tibble(
      weight = w, design = design, arm_form = form, observations = nobs(both),
      countries = n_distinct(used$GID_0),
      years = paste(range(used$year), collapse = "-"),
      arms_given_shocks_p = block_p(both, arm_names),
      shocks_given_arms_p = block_p(both, shock_names),
      shock_TM_alone = stats::coef(shocks_only)[[shock_names[[1]]]],
      shock_TM_alone_p = fixest::pvalue(shocks_only)[[shock_names[[1]]]],
      shock_RR_alone = stats::coef(shocks_only)[[shock_names[[2]]]],
      shock_RR_alone_p = fixest::pvalue(shocks_only)[[shock_names[[2]]]],
      shock_TM_both = stats::coef(both)[[shock_names[[1]]]],
      shock_TM_both_p = fixest::pvalue(both)[[shock_names[[1]]]],
      shock_RR_both = stats::coef(both)[[shock_names[[2]]]],
      shock_RR_both_p = fixest::pvalue(both)[[shock_names[[2]]]],
      wr2_arms = fixest::fitstat(arms_only, "wr2")[[1]],
      wr2_shocks = fixest::fitstat(shocks_only, "wr2")[[1]],
      wr2_both = fixest::fitstat(both, "wr2")[[1]],
      aic_both_minus_arms = stats::AIC(both) - stats::AIC(arms_only),
      bic_both_minus_arms = stats::BIC(both) - stats::BIC(arms_only),
      aic_both_minus_shocks = stats::AIC(both) - stats::AIC(shocks_only),
      # Overlap of the arm and shock measures within countries.
      within_cor_TM = stats::cor(within(arm_measure, used), within(shock_measure, used)),
      # Share of the within variance of the 30-year anomaly in each part.
      gap_share_TM = stats::var(within("gapTM", used)) / stats::var(within("zTM", used)),
      cycle_share_TM = stats::var(within("cycTM", used)) / stats::var(within("zTM", used))
    )

    # Arm effects of a 1-SD deviation at each long-run climate, in the
    # combined and the arms-only model (the same weights for both forms; a
    # 2-SD deviation has twice the linear and four times the quadratic effect).
    for (k in c("both", "arms_only")) {
      m <- models[[w]][[paste(design, form)]][[k]]
      m_terms <- names(stats::coef(m))
      for (i in seq_along(arms)) {
        arm <- arm_coefs[[i]]
        climate <- arms[[i]]
        for (value in AT_VALUES[[climate]]) {
          weights <- stats::setNames(c(1, value),
                                     c(arm, interaction_name(arm, climate, m_terms)))
          e <- lincom(m, weights)
          effects[[length(effects) + 1]] <- tibble::tibble(
            weight = w, design = design, arm_form = form, model = k,
            arm = ARM_LABELS[[i]], climate = value, effect = e$estimate,
            std_error = e$std_error, p_value = e$p_value
          )
        }
      }
    }
  }
}

tests <- bind_rows(tests)
effects <- bind_rows(effects)
utils::write.csv(tests, file.path(OUTPUT_DIR, "block_tests.csv"), row.names = FALSE)
utils::write.csv(effects, file.path(OUTPUT_DIR, "arm_effects.csv"), row.names = FALSE)

options(width = 220)
cat("\n== Blocks, shocks and fit ==\n")
print(as.data.frame(tests %>% mutate(across(where(is.double), ~ signif(.x, 3)))))
cat("\n== Arm effects per SD (combined model) ==\n")
print(as.data.frame(effects %>%
  filter(model == "both") %>%
  mutate(across(c(effect, std_error, p_value), ~ signif(.x, 3)))))
for (w in names(models)) for (key in names(models[[w]])) {
  cat("\n== Models,", key, "arms,", w, "weights ==\n")
  print(etable(models[[w]][[key]], headers = c("Arms only", "Shocks only", "Both"),
               digits = 4, fitstat = ~n + wr2 + aic + bic))
}
