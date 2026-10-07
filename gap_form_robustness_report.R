# Robustness of the pooled deviation forms (signed, linear arms, quadratic arms,
# each interacted with the long-run climate) across datasets and variable
# constructions, from the runs of test_gap_form_validation.R, and the criteria
# for choosing one form per variable.
#
# Runs (all: population-density weights 2000, deviation from the trailing
# 30-year mean unless stated, controls soviet + yugoslavia + conflict, growth
# outliers of 1991-1993 dropped, economic years from 1970):
#   datasets       PWT and WB x ERA5, CRU TS and UDelaware (whole-series SD and
#                  whole-series mean in the interaction)
#   construction   PWT x ERA5 with the deviation baseline (trailing 30-year or
#                  whole-series mean, the latter for temperature only), the SD
#                  (whole-series or trailing 30-year) and the interaction mean
#                  (whole-series or trailing 30-year), 2 x 2 x 2
# The main specification is PWT x ERA5, trailing baseline, whole-series SD and
# interaction mean.
#
# Criteria, per variable and run, for the three pooled forms:
#   aic, bic        information criteria of the (a) interacted models
#   shape           p-value of the bins added to the pooled form (all quarters)
#   holdout         mean squared prediction error on the 5% most extreme
#                   deviations of the variable, and equal-loss tests between forms
#   effects         implied effect of a +/-2 SD deviation in the coldest and
#                   hottest (driest and wettest) quarter, and its sign agreement
#                   with the main specification
#
# Run the validation runs first, e.g.
#   Rscript test_gap_form_validation.R econ=PWT source=ERA5 deviation=anomaly \
#     controls=soviet,yugoslavia,conflict outliers=pinatubo start_year=1970
# Output: results/gap_form_robustness/ (CSV tables and figures; REPORT.md)

source("load_functions.R")

VALIDATION_DIR <- file.path("results", "gap_form_validation")
OUTPUT_DIR <- file.path("results", "gap_form_robustness")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  utils::write.csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
}

COMMON <- "_soviet_yugoslavia_conflict_dropped_pinatubo_from1970"
run_dir <- function(econ, source, deviation = "anomaly", sd = "all", interaction = "all") {
  paste0(econ, "_", source, "_pop2000_", deviation, COMMON,
         if (sd == "lag_30") "_sd30" else "", if (interaction == "lag_30") "_int30" else "")
}
RUNS <- bind_rows(
  tidyr::expand_grid(econ = c("PWT", "WB"), source = c("ERA5", "CRUTS", "UDelaware")) %>%
    mutate(block = "datasets", deviation = "anomaly", sd = "all", interaction = "all"),
  tidyr::expand_grid(deviation = c("anomaly", "anomaly_all"), sd = c("all", "lag_30"),
                     interaction = c("all", "lag_30")) %>%
    mutate(block = "construction", econ = "PWT", source = "ERA5")
) %>%
  mutate(
    dir = mapply(run_dir, econ, source, deviation, sd, interaction),
    source_label = c(ERA5 = "ERA5", CRUTS = "CRU", UDelaware = "UDel")[source],
    construction = paste0(
      ifelse(deviation == "anomaly", "trailing", "fixed"), " baseline, ",
      ifelse(sd == "all", "fixed", "30-yr"), " SD, ",
      ifelse(interaction == "all", "fixed", "30-yr"), " interaction"
    ),
    run = ifelse(block == "datasets", paste(econ, "x", source_label), construction),
    main = econ == "PWT" & source == "ERA5" & deviation == "anomaly" & sd == "all" &
      interaction == "all"
  )
missing <- RUNS$dir[!dir.exists(file.path(VALIDATION_DIR, RUNS$dir))]
if (length(missing)) stop("Missing validation runs: ", paste(missing, collapse = ", "))

read_run <- function(file) {
  bind_rows(lapply(seq_len(nrow(RUNS)), function(i) {
    x <- utils::read.csv(file.path(VALIDATION_DIR, RUNS$dir[[i]], file), stringsAsFactors = FALSE)
    cbind(RUNS[i, c("block", "run", "main")], x)
  }))
}
FORMS <- c("signed", "linear_arms", "quadratic_arms")
FORM_LABELS <- c(signed = "Signed (BHM style)", linear_arms = "Linear arms",
                 quadratic_arms = "Quadratic arms")
VARIABLE_LABELS <- c(TM = "Temperature", RR = "Precipitation")
pooled <- function(model) sub("_interacted$", "", model)

# Samples ----------------------------------------------------------------------------------

samples <- read_run("sample_summary.csv") %>%
  select(block, run, observations, countries, first_year, last_year, outliers_dropped,
         temperature_holdout_threshold, temperature_holdout, precipitation_holdout_threshold,
         precipitation_holdout)
write_out(samples, "samples")

# Criteria 1-2: information criteria of the pooled forms ---------------------------------

fit <- read_run("a_model_fit.csv") %>%
  filter(model %in% paste0(FORMS, "_interacted")) %>%
  mutate(form = pooled(model)) %>%
  group_by(block, run, variable) %>%
  mutate(delta_aic = aic - min(aic), delta_bic = bic - min(bic),
         best_aic = delta_aic == 0, best_bic = delta_bic == 0) %>%
  ungroup() %>%
  select(block, run, main, variable, form, parameters, aic, bic, within_r2, delta_aic, delta_bic,
         best_aic, best_bic)
write_out(fit, "criteria_information")

# Criterion 3: shape adequacy (bins added to the pooled form) ----------------------------

shape <- read_run("a_tests.csv") %>%
  filter(test_type == "shape", form %in% paste0(FORMS, "_interacted")) %>%
  mutate(form = pooled(form)) %>%
  group_by(block, run, main, variable, form) %>%
  summarise(p_all_quarters = p_value[test == "bins = 0, all quarters"],
            p_min_quarter = min(p_value[test != "bins = 0, all quarters"]),
            .groups = "drop") %>%
  mutate(adequate = p_all_quarters > 0.05)
write_out(shape, "criteria_shape")

# Criterion 4: out-of-sample prediction of the extreme deviations -------------------------

holdout_perf <- read_run("b_holdout_performance.csv") %>%
  filter(side == "all", model %in% FORMS) %>%
  rename(form = model) %>%
  group_by(block, run, variable) %>%
  mutate(delta_mse_pp2 = 1e4 * (mse - min(mse)), best_holdout = delta_mse_pp2 == 0) %>%
  ungroup() %>%
  select(block, run, main, variable, form, observations, rmse, delta_mse_pp2, best_holdout)
holdout_tests <- read_run("b_holdout_tests.csv") %>%
  filter(side == "all") %>%
  transmute(block, run, main, variable, model, against,
            mse_difference_pp2 = 1e4 * mse_difference, conf_low_pp2 = 1e4 * conf_low,
            conf_high_pp2 = 1e4 * conf_high, p_value)
write_out(holdout_perf, "criteria_holdout")
write_out(holdout_tests, "criteria_holdout_tests")

# Criterion 5: implied effects at the climate extremes -----------------------------------

curves <- read_run("a_form_curves.csv") %>%
  filter(estimation == "interacted", group %in% c(1, 4))
effect_at <- function(x, at) {
  x %>%
    group_by(block, run, main, variable, form, group) %>%
    filter(min(deviation) <= at + 1e-9, max(deviation) >= at - 1e-9) %>%
    slice_min(abs(deviation - at), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(at = at)
}
effects <- bind_rows(effect_at(curves, -2), effect_at(curves, 2)) %>%
  mutate(quarter = ifelse(group == 1,
                          ifelse(variable == "TM", "coldest", "driest"),
                          ifelse(variable == "TM", "hottest", "wettest")),
         significant = conf_low > 0 | conf_high < 0) %>%
  select(block, run, main, variable, form, quarter, climate, at, estimate, conf_low, conf_high,
         significant)
# The main specification is in both blocks; its signs are the reference.
main_sign <- effects %>%
  filter(main, block == "datasets") %>%
  transmute(variable, form, quarter, at, main_sign = sign(estimate))
effects <- effects %>%
  left_join(main_sign, by = c("variable", "form", "quarter", "at")) %>%
  mutate(agrees_with_main = sign(estimate) == main_sign)
write_out(effects, "criteria_effects")

# Scorecard ----------------------------------------------------------------------------

# Per block (the main specification is in both) and variable: in how many runs
# each form is best, adequate, or significantly better or worse out of sample.
keys <- c("block", "variable", "form")
holdout_wins <- holdout_tests %>%
  filter(model %in% FORMS, against %in% FORMS, p_value < 0.05) %>%
  transmute(block, variable, winner = ifelse(mse_difference_pp2 < 0, model, against),
            loser = ifelse(mse_difference_pp2 < 0, against, model)) %>%
  tidyr::pivot_longer(c(winner, loser), names_to = "role", values_to = "form") %>%
  count(block, variable, form, role) %>%
  tidyr::pivot_wider(names_from = role, values_from = n, values_fill = 0,
                     names_prefix = "holdout_significantly_")
for (col in c("holdout_significantly_winner", "holdout_significantly_loser")) {
  if (!col %in% names(holdout_wins)) holdout_wins[[col]] <- 0L
}
scorecard <- fit %>%
  group_by(block, variable, form) %>%
  summarise(runs = n(), best_aic = sum(best_aic), mean_delta_aic = mean(delta_aic),
            best_bic = sum(best_bic), mean_delta_bic = mean(delta_bic), .groups = "drop") %>%
  left_join(shape %>% group_by(block, variable, form) %>%
              summarise(shape_adequate = sum(adequate), median_shape_p = median(p_all_quarters),
                        .groups = "drop"), by = keys) %>%
  left_join(holdout_perf %>% group_by(block, variable, form) %>%
              summarise(best_holdout = sum(best_holdout),
                        mean_delta_mse_pp2 = mean(delta_mse_pp2), .groups = "drop"),
            by = keys) %>%
  left_join(effects %>% group_by(block, variable, form) %>%
              summarise(effect_sign_agreement = mean(agrees_with_main, na.rm = TRUE),
                        effects_significant = mean(significant), .groups = "drop"),
            by = keys) %>%
  left_join(holdout_wins, by = keys) %>%
  mutate(across(starts_with("holdout_significantly_"), ~ dplyr::coalesce(as.integer(.x), 0L))) %>%
  arrange(block, variable, form)
write_out(scorecard, "scorecard")
options(width = 220)
print(as.data.frame(scorecard %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)

# Figures ------------------------------------------------------------------------------

form_levels <- unname(FORM_LABELS)
colours <- grid_colours(form_levels)
shapes <- grid_shapes(form_levels)
run_levels <- c(RUNS$run[RUNS$block == "datasets"], RUNS$run[RUNS$block == "construction"])
block_labels <- c(datasets = "Datasets (main construction)",
                  construction = "Constructions (PWT x ERA5)")
prepare <- function(x) {
  x %>%
    mutate(form = factor(FORM_LABELS[form], levels = form_levels),
           run = factor(run, levels = rev(unique(run_levels))),
           block = factor(block_labels[block], levels = block_labels),
           variable = factor(VARIABLE_LABELS[variable], levels = VARIABLE_LABELS))
}
caption <- paste0("All runs: population-density weights 2000, controls (post-Soviet and ",
                  "post-Yugoslav collapse years, conflict), 1991-93 growth outliers dropped, ",
                  "1970 onwards, country and year effects with country trends; errors ",
                  "clustered by country. Main specification: PWT x ERA5, trailing baseline, ",
                  "fixed SD and interaction.")
caption <- paste(strwrap(caption, 170), collapse = "\n")
dodge <- position_dodge(width = 0.6)
save_plot <- function(p, name, height = 8.5) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), p, width = 13, height = height, dpi = 300)
}

p_aic <- ggplot(prepare(fit), aes(delta_aic, run, colour = form, shape = form)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_point(position = dodge, size = 2.4) +
  facet_grid(block ~ variable, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  labs(x = "AIC minus the lowest AIC of the three pooled forms (0 = best)", y = NULL,
       title = "In-sample fit of the pooled forms across datasets and constructions",
       subtitle = paste0("Each form interacted with the long-run climate; the other variable ",
                         "in the flexible (signed + quadratic arms) form."),
       caption = caption) +
  theme_results()
save_plot(p_aic, "fig1_aic")

p_shape <- ggplot(prepare(shape), aes(p_all_quarters, run, colour = form, shape = form)) +
  geom_vline(xintercept = 0.05, colour = GRID_MUTED, linewidth = 0.3, linetype = "22") +
  geom_point(position = dodge, size = 2.4) +
  facet_grid(block ~ variable, scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.05, 0.2, 1), labels = c("0.001", "0.01", "0.05", "0.2", "1")) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  labs(x = "p-value, bins added to the pooled form (all quarters; left of the dashed line: rejected)",
       y = NULL, title = "Shape adequacy: do free bins reject the pooled form?",
       subtitle = "Joint Wald test of the quarter x deviation-bin dummies, clustered by country.",
       caption = caption) +
  theme_results()
save_plot(p_shape, "fig2_shape")

comparisons <- c("quadratic_arms vs linear_arms" = "Quadratic arms vs linear arms",
                 "linear_arms vs signed" = "Linear arms vs signed",
                 "quadratic_arms vs signed" = "Quadratic arms vs signed")
holdout_plot <- holdout_tests %>%
  mutate(comparison = paste(model, "vs", against)) %>%
  filter(comparison %in% names(comparisons)) %>%
  mutate(comparison = factor(comparisons[comparison], levels = comparisons),
         run = factor(run, levels = rev(unique(run_levels))),
         block = factor(block_labels[block], levels = block_labels),
         variable = factor(VARIABLE_LABELS[variable], levels = VARIABLE_LABELS))
p_holdout <- ggplot(holdout_plot, aes(mse_difference_pp2, run, colour = comparison,
                                      shape = comparison)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf_low_pp2, xmax = conf_high_pp2), position = dodge,
                 linewidth = 0.5) +
  geom_point(position = dodge, size = 2.2) +
  facet_grid(block ~ variable, scales = "free", space = "free_y") +
  scale_colour_manual(values = grid_colours(unname(comparisons)), name = NULL) +
  scale_shape_manual(values = grid_shapes(unname(comparisons)), name = NULL) +
  labs(x = "Difference in mean squared prediction error (pp²); negative = first form predicts better",
       y = NULL, title = "Out-of-sample prediction of each variable's 5% most extreme deviations",
       subtitle = "Pooled forms re-estimated without the held-out country-years; 95% CI clustered by country.",
       caption = caption) +
  theme_results()
save_plot(p_holdout, "fig3_holdout")

effects_plot <- effects %>%
  mutate(point = paste0(ifelse(at > 0, "+2 SD", "-2 SD"), ", ", quarter, " quarter"),
         across(c(estimate, conf_low, conf_high), ~ 100 * .x)) %>%
  prepare()
point_levels <- c("-2 SD, coldest quarter", "+2 SD, coldest quarter", "-2 SD, hottest quarter",
                  "+2 SD, hottest quarter", "-2 SD, driest quarter", "+2 SD, driest quarter",
                  "-2 SD, wettest quarter", "+2 SD, wettest quarter")
effects_plot$point <- factor(effects_plot$point, levels = point_levels)
p_effects <- ggplot(filter(effects_plot, block == block_labels[["datasets"]]),
                    aes(estimate, run, colour = form, shape = form)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high), position = dodge, linewidth = 0.5) +
  geom_point(position = dodge, size = 2.2) +
  facet_wrap(~point, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  labs(x = "Effect on growth that year (pp), relative to the reference bin; 95% CI", y = NULL,
       title = "Implied effects of the pooled forms at the climate extremes, across datasets",
       subtitle = paste0("Evaluated at the mean long-run climate of the coldest and hottest ",
                         "(driest and wettest) quarter of each run."),
       caption = caption) +
  theme_results()
save_plot(p_effects, "fig4_effects_datasets", height = 9)

p_effects_c <- ggplot(filter(effects_plot, block == block_labels[["construction"]]),
                      aes(estimate, run, colour = form, shape = form)) +
  geom_vline(xintercept = 0, colour = GRID_MUTED, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf_low, xmax = conf_high), position = dodge, linewidth = 0.5) +
  geom_point(position = dodge, size = 2.2) +
  facet_wrap(~point, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  labs(x = "Effect on growth that year (pp), relative to the reference bin; 95% CI", y = NULL,
       title = "Implied effects of the pooled forms at the climate extremes, across constructions (PWT x ERA5)",
       subtitle = paste0("Evaluated at the mean long-run climate of the coldest and hottest ",
                         "(driest and wettest) quarter. SD units differ between fixed and 30-year SD."),
       caption = caption) +
  theme_results()
save_plot(p_effects_c, "fig5_effects_constructions", height = 9)

message("Report tables and figures written to ", normalizePath(OUTPUT_DIR))
