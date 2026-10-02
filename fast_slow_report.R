# Report for the fast-versus-slow warming design
# (FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md, sections 7, 9 and 12), built
# only from the saved outputs of fast_slow_stage*.R; nothing is refitted.
# Output: results/fast_slow_warming/report.md and checklist.csv.

source("load_functions.R")

design <- fs_design()
rd <- function(...) utils::read.csv(fs_path(...))
rj <- function(...) jsonlite::read_json(fs_path(...), simplifyVector = TRUE)
kable <- function(x, digits = 3) {
  paste(knitr::kable(x, format = "pipe", digits = digits), collapse = "\n")
}
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

flow <- rd("data", "sample_flow.csv")
joins <- rj("qa", "join_checks.json")
schema <- rj("qa", "schema_checks.json")
scen <- rj("support", "scenario_choice.json")
tails <- rd("support", "tail_support.csv")
idvar <- rd("support", "rate_identifying_variation.csv")
power_design <- rd("support", "power_design_range.csv")
s37 <- readRDS(fs_path("models", "stage3_7_summary.rds"))
pc <- rd("models", "path_contrasts.csv")
power_cal <- rd("tables", "power_residual_calibrated.csv")
fe_diag <- rd("tables", "fe_variant_diagnostics.csv")
restr <- rd("tables", "slope_restriction_tests.csv")
lp <- rd("tables", "local_projections.csv")
dl <- rd("tables", "distributed_lag_almon.csv")
ld <- rd("tables", "long_differences.csv")
anom <- rd("tables", "anomaly_model_comparison.csv")
cool <- rd("tables", "warming_vs_cooling.csv")
spec_curve <- rd("tables", "specification_curve.csv")
infv <- rd("tables", "inference_variants.csv")
leads <- rd("tables", "falsification_leads.csv")
placebo <- rd("tables", "falsification_future_rate.csv")
perm <- rd("tables", "falsification_permutation.csv")
loo <- rd("tables", "leave_one_out_summary.csv")
sectors <- rd("tables", "mechanisms_sectors.csv")
het <- rd("tables", "heterogeneity_interactions.csv")
het_gaps <- rd("tables", "heterogeneity_path_gaps.csv")
coefs <- rd("tables", "coefficients.csv")

conf <- s37$conf
H <- s37$H
sesoi <- design$scenario$sesoi_log_points
K <- design$climate$rate_window
main <- pc %>% filter(model_id == "rate_broad", scenario_id == conf$id,
                      is.na(estimand))
canon <- pc %>% filter(model_id == "rate_broad",
                       scenario_id == "canonical_extrapolative")
at <- function(tab, comp, h = H) tab[tab$component == comp & tab$horizon == h, ]
tot <- at(main, "total")
spd <- at(main, "speed")
lvl <- at(main, "level")

# Robustness envelope (section 7): core specifications only.
core <- spec_curve %>% filter(tier %in% c("confirmatory", "core_robustness"),
                              component %in% c("total", "speed"))
envelope <- core %>%
  group_by(component) %>%
  summarise(specifications = n(),
            envelope_low = min(ci95_low, na.rm = TRUE),
            envelope_high = max(ci95_high, na.rm = TRUE),
            n_excluding_zero_negative = sum(ci95_high < 0, na.rm = TRUE),
            n_excluding_zero_positive = sum(ci95_low > 0, na.rm = TRUE),
            min_estimate = min(estimate), max_estimate = max(estimate),
            .groups = "drop")
cy <- pc %>% filter(model_id == "rate_country_year", scenario_id == conf$id,
                    horizon == H, component %in% c("total", "speed"))

# Checklist (section 12) ----------------------------------------------------------

checklist <- tibble::tribble(
  ~item, ~question, ~answer, ~note,
  1, "Data versions, raw hashes and API dates frozen?", "yes",
  "data/data_manifest.csv (SHA-256 of DOSE, the WCD monthly files, GADM 4.1, KG zones); WDI not used",
  2, "Geography reconciled spatially rather than by unverified names/codes?", "no",
  sprintf("DOSE GID_1 joined to GADM 4.1 climate means by identifier; %d of %d primary regions match exactly by name; coverage and weight mass not computable",
          joins$primary_exact_name_match, joins$primary_regions),
  3, "Climate rate variables created without access to outcomes?", "yes",
  "fs_add_climate_terms() refuses outcome columns; Stage 2 reads no GDP column",
  4, "Canonical path checked against joint historical support?", "yes",
  sprintf("min weighted supported share %.2f; %s", scen$canonical$min_supported_share,
          scen$canonical$label),
  5, "All primary growth observations adjacent-year differences?", "yes", "assertion 2",
  6, "Documented structural breaks excluded or flagged?", "yes",
  "StructChange 1, 2, 3 excluded; 1 and 3 added back in named sensitivities",
  7, "Nested models on an identical row set?", "yes", "assertion 9 (stopifnot on obs())",
  8, "Speed represented nonlinearly or dynamically?", "yes",
  "signed quadratic hinges; Almon lags; local projections",
  9, "Fast and slow paths equal at start and endpoint, common horizon?", "yes",
  "tests/test_fast_slow_warming.R",
  10, "Total and incremental speed contrasts both reported?", "yes", "",
  11, "Broad and within-country-local variation kept distinct?", "yes",
  "country-by-year model reported separately",
  12, "Spatial and temporal dependence reflected in uncertainty?", "yes",
  "two-way CGM, wild bootstrap (country and year DGPs), Conley 500-2000 km with 5/10-year Bartlett",
  13, "Climate-product, weighting and model uncertainties shown without selective reporting?", "no",
  "model uncertainty shown in full; CRU/UDel and fixed-2000 population weights not run (preferred data only)",
  14, "Equivalence assessed with an ex ante SESOI rather than a non-significant p-value?", "no",
  "SESOI (1 log point) not affirmed by the study owner; equivalence reported as provisional only",
  15, "Unsupported scenario cells explicitly labelled?", "yes", "",
  16, "PWT and WDI described as triangulation?", "not applicable",
  "national module not run (preferred data only)",
  17, "Mechanisms and heterogeneity kept secondary?", "yes", "",
  18, "Report permits and states an inconclusive or null result?", "yes", ""
)
fs_write(checklist, "checklist.csv")
failed_items <- checklist %>% filter(answer == "no")

# Text helpers ------------------------------------------------------------------------

gap_line <- function(row) {
  sprintf("%s log points (95%% CI [%s, %s]; 90%% CI [%s, %s])",
          f2(row$estimate), f2(row$ci95_low), f2(row$ci95_high),
          f2(row$ci90_low), f2(row$ci90_high))
}
source_of_variation <- "country-wide and within-country variation net of region effects, year effects and country linear trends"
support_word <- if (isTRUE(conf$supported)) "inside" else "outside"
power <- s37$power_at_sesoi

decision_table <- main %>%
  filter(component %in% c("total", "level", "speed")) %>%
  transmute(component, horizon, estimate, ci95_low, ci95_high, ci90_low,
            ci90_high, support_status, direction, materiality,
            not_materially_worse_one_sided)
canon_table <- canon %>%
  transmute(component, horizon, estimate, ci95_low, ci95_high, support_status,
            direction)

power_design_sesoi <- power_design %>%
  filter(dgp == "planted_speed_curvature", planted_gap == -sesoi) %>%
  group_by(contrast) %>%
  summarise(min_p_fast_worse = min(p_fast_worse),
            median_p_fast_worse = stats::median(p_fast_worse),
            max_p_fast_worse = max(p_fast_worse),
            min_coverage95 = min(coverage95), .groups = "drop")
misspec <- power_design %>%
  filter(dgp == "misspecified_lags", contrast == "speed") %>%
  group_by(planted_gap) %>%
  summarise(mean_estimated_speed_gap = mean(mean_estimate),
            p_fast_worse = mean(p_fast_worse), .groups = "drop")

lp_show <- lp %>%
  filter(version == "with_lags", horizon %in% c(0, 1, 2, 5, 10)) %>%
  select(contrast, horizon, estimate, ci95_low, ci95_high, band95_low,
         band95_high, n_obs)
rate_coefs <- coefs %>% filter(model == "rate_broad") %>%
  select(term, estimate, std_error, p_value, conf_low, conf_high)
level_coefs <- coefs %>% filter(model %in% c("level_broad", "level_country_year"),
                                term %in% c("Tc", "Tc2")) %>%
  select(model, term, estimate, std_error, p_value)

anom_show <- anom %>%
  select(model, normal_years, contrast, estimate, ci95_low, ci95_high,
         within_adj_r2, rmse_leave_country_out, rmse_leave_decade_out,
         rmse_tail_above_p90_positive_rate)
spec_show <- spec_curve %>%
  filter(component %in% c("total", "speed")) %>%
  select(spec_id, tier, group, component, estimate, ci95_low, ci95_high) %>%
  arrange(component, group, spec_id)
inf_show <- infv %>% filter(!is.na(se)) %>%
  select(inference, contrast, se, ci95_low, ci95_high)

# Mechanical interpretation rules (sections 7, 8 and 11) ----------------------------

rate_cv <- anom %>% filter(model == "rate", normal_years ==
                             design$climate$anomaly_normal_years) %>% slice(1)
anom_cv <- anom %>% filter(model == "anomaly", normal_years ==
                             design$climate$anomaly_normal_years)
anom_block <- anom_cv %>% filter(grepl("^speed", contrast))
anom_fits_better <- anom_cv$rmse_leave_country_out[1] <=
  rate_cv$rmse_leave_country_out & anom_cv$rmse_leave_decade_out[1] <=
  rate_cv$rmse_leave_decade_out
anom_block_harmful <- anom_block$ci95_high < 0
speed_harmful <- spd$direction == "fast worse"
cool_speed <- cool %>% filter(direction == "cooling (-M)",
                              grepl("^speed", contrast))
no_trend <- pc %>% filter(model_id == "rate_no_trend", scenario_id == conf$id,
                          horizon == H, component %in% c("total", "speed"))
cy_dir <- cy$direction
names(cy_dir) <- cy$component
opposite <- function(a, b) (a == "fast worse" && b == "fast less harmful") ||
  (a == "fast less harmful" && b == "fast worse")
local_verdict <- if (opposite(cy_dir[["total"]], tot$direction) ||
                     opposite(cy_dir[["speed"]], spd$direction)) {
  "the local design contradicts the broad design, so the causal interpretation is unresolved."
} else if (all(cy_dir == c(total = tot$direction, speed = spd$direction)[names(cy_dir)])) {
  "same directional classification as the broad design."
} else {
  "imprecise rather than contradictory; any broad-design result relies on country-wide variation."
}
far_spec <- core %>%
  filter(component == "speed") %>%
  mutate(distance = abs(estimate - spd$estimate)) %>%
  slice_max(distance, n = 1, with_ties = FALSE)
rule_lines <- c(
  sprintf("- *Total versus speed.* %s",
          if (tot$direction == "fast worse" && !speed_harmful)
            "The fast path is worse in total, but there is no identified incremental speed penalty: earlier arrival at higher temperatures, not speed itself, explains the total gap under the model."
          else if (speed_harmful)
            "The explicit rate/sequence block itself makes the fast path worse."
          else sprintf("Neither the total gap nor the incremental speed component is directionally resolved. The total gap is imprecise (power %.2f at the SESOI, below 0.80), so its null is inconclusive rather than evidence of no effect. The speed component is precisely estimated (power %.2f) and its point estimate is %s zero.",
                       power[["total"]], power[["speed"]],
                       if (spd$estimate > 0) "above (fast *less* harmful), though not significantly different from," else "below, though not significantly different from,")),
  sprintf("- *Equivalence.* The 90%% interval of the speed component [%s, %s] %s the provisional +/-%s bounds; %s",
          f2(spd$ci90_low), f2(spd$ci90_high),
          if (spd$ci90_low >= -sesoi && spd$ci90_high <= sesoi) "lies inside" else "is not inside",
          f2(sesoi, 0),
          "this would support practical equivalence of fast and slow warming in the incremental speed channel only if the study owner affirms the SESOI ex ante."),
  sprintf("- *Adaptation/mismatch (Stage 8).* The trailing-normal anomaly model %s the rate model out of sample (leave-country-out RMSE %s vs %s; leave-decade-out %s vs %s), and its anomaly block implies a fast-slow gap of %s (95%% CI [%s, %s]). Rule: %s",
          if (anom_fits_better) "fits at least as well as" else "does not beat",
          f2(anom_cv$rmse_leave_country_out[1], 3), f2(rate_cv$rmse_leave_country_out, 3),
          f2(anom_cv$rmse_leave_decade_out[1], 3), f2(rate_cv$rmse_leave_decade_out, 3),
          f2(anom_block$estimate), f2(anom_block$ci95_low), f2(anom_block$ci95_high),
          if (anom_fits_better && anom_block_harmful && !speed_harmful)
            "the data favour a surprise/adaptation-mismatch account over a generic rate penalty."
          else "the anomaly model does not succeed where speed curvature fails (its mismatch block is not harmful at 95%), so no adaptation-mismatch conclusion is drawn."),
  sprintf("- *Warming versus cooling.* Cooling-path speed component %s (95%% CI [%s, %s]); %s",
          f2(cool_speed$estimate), f2(cool_speed$ci95_low), f2(cool_speed$ci95_high),
          if (speed_harmful && cool_speed$ci95_high < 0)
            "both rapid warming and rapid cooling are harmful: a rapid-change penalty."
          else "no rapid-change penalty is identified in either direction."),
  sprintf("- *Broad versus local design.* Country-by-year effects (within-country variation only) give total %s [%s, %s] and speed %s [%s, %s]: %s",
          f2(cy$estimate[cy$component == "total"]), f2(cy$ci95_low[cy$component == "total"]),
          f2(cy$ci95_high[cy$component == "total"]),
          f2(cy$estimate[cy$component == "speed"]), f2(cy$ci95_low[cy$component == "speed"]),
          f2(cy$ci95_high[cy$component == "speed"]),
          local_verdict),
  sprintf("- *Trends.* Without country trends the gaps are total %s (%s) and speed %s (%s); %s",
          f2(no_trend$estimate[no_trend$component == "total"]),
          no_trend$direction[no_trend$component == "total"],
          f2(no_trend$estimate[no_trend$component == "speed"]),
          no_trend$direction[no_trend$component == "speed"],
          if (identical(no_trend$direction[no_trend$component == "total"], tot$direction) &&
              identical(no_trend$direction[no_trend$component == "speed"], spd$direction))
            "the directional conclusions do not depend on the trend restriction."
          else "the directional conclusions are trend-sensitive."),
  sprintf("- *Robustness envelope.* Across %d core specifications the speed estimates range from %s to %s and the total from %s to %s; %s",
          envelope$specifications[1],
          f2(envelope$min_estimate[envelope$component == "speed"]),
          f2(envelope$max_estimate[envelope$component == "speed"]),
          f2(envelope$min_estimate[envelope$component == "total"]),
          f2(envelope$max_estimate[envelope$component == "total"]),
          sprintf("%s The specification furthest from the primary speed estimate is `%s` (%s).",
                  if (any(c(tot$direction, spd$direction) != "direction unresolved"))
                    "the envelope rule applies to the resolved components."
                  else "no component is directionally resolved, so the envelope documents sensitivity only.",
                  far_spec$spec_id, f2(far_spec$estimate))),
  "- *Canonical +1 C contrast.* Not identified by these data: it lies outside the region-specific temperature support. Its model-based extrapolation is shown only under that label."
)

decision_sentence <- function(row, label) {
  sprintf("- **%s at year %d:** %s -- *%s*; *%s*.", label, row$horizon,
          gap_line(row), row$direction, row$materiality)
}

# Report -----------------------------------------------------------------------------------

report <- c(
  "# Fast versus slow warming: DOSE x ERA5 implementation",
  "",
  sprintf("Design `config/fast_slow_design.yml` v%s (SHA-256 `%s`), frozen %s. Protocol: `FAST_VS_SLOW_WARMING_EMPIRICAL_STRATEGY.md`.",
          design$design_version, substr(design$hash, 1, 16), design$freeze_date),
  "",
  "**Status.** This is a frozen analysis plan, not a preregistration: `prior_related_results_seen` is unset and must be filled by the study owner. The 1-log-point SESOI has not been affirmed, so practical-equivalence claims are disabled and materiality labels are provisional. Headline claims are downgraded for the failed checklist items listed in section 10.",
  "",
  "## 1. Headline",
  "",
  sprintf("For paths producing the same %s C warming from each region's 1991-2010 climate, the model-implied real GDP per capita gap at year %d for the %d-year path relative to the %d-year path was **%s log points** (95%% CI [%s, %s]). The component attributed within the model to the explicit rate/sequence function was **%s** (95%% CI [%s, %s]). The contrast was %s historical support and was estimated primarily from %s.",
          f2(conf$total_warming), H, conf$fast_years, conf$slow_years,
          f2(tot$estimate), f2(tot$ci95_low), f2(tot$ci95_high),
          f2(spd$estimate), f2(spd$ci95_low), f2(spd$ci95_high),
          support_word, source_of_variation),
  "",
  sprintf("This contrast is the deterministic within-support fallback of section 1. The canonical +%s C in %d vs %d years contrast is **outside historical support** (weighted supported share %.2f against the 0.95 rule: after the ramp, B + %s C exceeds the region's own observed temperature range for about %.0f%% of the weight). It is reported only as a labelled extrapolative policy scenario (section 4).",
          f2(design$scenario$total_warming_c, 0), design$scenario$fast_years,
          design$scenario$slow_years, scen$canonical$min_supported_share,
          f2(design$scenario$total_warming_c, 0),
          100 * (1 - scen$canonical$min_supported_share)),
  "",
  "Decisions (section 7; negative = fast worse; intervals are wild-cluster bootstrap-t, the wider of the country and year Rademacher DGPs, 9,999 draws each):",
  "",
  decision_sentence(tot, "Total fast-slow gap"),
  decision_sentence(lvl, "Level-timing component"),
  decision_sentence(spd, "Incremental speed penalty"),
  "",
  sprintf("Residual-calibrated power to detect a -%s log-point gap: total %.2f, speed %.2f (rule: 0.80).",
          f2(sesoi, 0), power[["total"]], power[["speed"]]),
  "",
  "Interpretation under the protocol's rules:",
  "",
  rule_lines,
  "",
  "## 2. Decision table",
  "",
  "Confirmatory contrast (`rate_broad`, primary sample, path weights proportional to region observations):",
  "",
  kable(decision_table),
  "",
  "Canonical +1 C in 5 vs 20 years (extrapolative policy scenario; not identified):",
  "",
  kable(canon_table),
  "",
  "SESOI sensitivity (0.5, 1, 2 log points): `tables/decision_sesoi_sensitivity.csv`.",
  "",
  "Robustness envelope over the confirmatory and core-robustness specifications at the confirmatory horizon:",
  "",
  kable(envelope),
  "",
  "Broad versus local (country-by-year) design at the confirmatory horizon:",
  "",
  kable(cy %>% select(component, estimate, ci95_low, ci95_high, direction)),
  "",
  "## 3. Sample and geography",
  "",
  kable(flow %>% select(step, removed, remaining, regions, countries)),
  "",
  sprintf("Primary sample: %s region-years, %d regions, %d countries, %s. Constant-2015-LCU and fixed-2015-USD growth agree to %.1e log points over %s pairs (assertion 5). T0 = %.2f C.",
          format(schema$primary_n, big.mark = ","), schema$primary_regions,
          schema$primary_countries, schema$primary_years,
          schema$assertion_5_fixed_usd_growth_max_gap_log_points,
          format(schema$assertion_5_pairs_compared, big.mark = ","),
          schema$T0),
  "",
  sprintf("**Geography gate not met.** The Weighted Climate Dataset series are GADM 4.1 ADM1 means; DOSE regions are joined by GID_1 string. Of %d primary regions, %d have an exact normalized name match with GADM 4.1 and %d a match or containment; polygon equivalence, target coverage and weight mass cannot be computed without the DOSE polygons (section 4.3). Mismatches are listed in `data/geography_crosswalk.csv`.",
          joins$primary_regions, joins$primary_exact_name_match,
          joins$primary_name_match_or_contained),
  "",
  "## 4. Outcome-blind support and power (Stage 2)",
  "",
  sprintf("Five-year slopes (C/year): 1st-99th percentiles [%s, %s]; 90th percentile of positive slopes %s. Canonical ramp rates %s (fast) and %s (slow) are inside the rate support; the binding constraint is temperature, not rate.",
          f2(scen$rate_bounds_1_99[1], 3), f2(scen$rate_bounds_1_99[2], 3),
          f2(tails$threshold[2], 3),
          f2(design$scenario$total_warming_c / design$scenario$fast_years, 3),
          f2(design$scenario$total_warming_c / design$scenario$slow_years, 3)),
  "",
  kable(tails %>% select(measure, threshold, region_years, regions, countries,
                         decades, passes)),
  "",
  "Identifying variation of the slope after each fixed-effect projection:",
  "",
  kable(idvar %>% select(fixed_effects, residual_sd,
                         residual_share_of_variance, effective_countries,
                         meaningful_countries, weakly_identified)),
  "",
  "Design-range power (1,000 simulations per cell, 36 noise cells; analytic two-way intervals), planted gap = -SESOI:",
  "",
  kable(power_design_sesoi),
  "",
  "Misspecified lag DGP (lags 0-4, weights 1, .75, .5, .25, 0; no speed effect by construction): estimated incremental speed gap by planted total gap:",
  "",
  kable(misspec),
  "",
  "Residual-calibrated power (Stage 3; level-benchmark residuals, wild country/year resampling, 1,000 draws):",
  "",
  kable(power_cal),
  "",
  "## 5. Level benchmark and rate model (Stages 3-4)",
  "",
  "![](figures/level_marginal_response.png)",
  "",
  kable(level_coefs),
  "",
  "Primary rate model coefficients (two-way country and year clustered):",
  "",
  kable(rate_coefs),
  "",
  "Fixed-effect variants (identical rows except the zone-year model, which drops regions without a Koppen-Geiger zone):",
  "",
  kable(fe_diag %>% select(fixed_effects, n_obs, partial_r2_rate_block,
                           max_vif, condition_index,
                           rate_sd_after_projection, meaningful_countries,
                           weakly_identified)),
  "",
  "## 6. Path contrasts (Stage 6)",
  "",
  "![](figures/path_gap_profile.png)",
  "",
  sprintf("The level-timing component stops changing once both paths reach the endpoint (year %d), and the speed component once the %d-year window has washed out (year %d); reporting years %s.",
          conf$slow_years, K, conf$slow_years + K - 1,
          paste(sort(unique(main$horizon)), collapse = ", ")),
  "",
  "Estimands and subgroups (fixed-baseline population uses 2000 population; regions without it are excluded from that estimand). The rate terms do not depend on the baseline temperature, so in the homogeneous model the subgroup speed gaps are identical by construction; only the level-timing part varies with the subgroup's climate. Heterogeneous responses are in section 9.",
  "",
  kable(pc %>% filter(scenario_id %in% c(conf$id, paste0(conf$id, "_equal_degree_years")),
                      horizon == H, component %in% c("total", "speed"),
                      !is.na(estimand) | grepl("equal_degree", scenario_id)) %>%
          transmute(model_id, estimand = dplyr::coalesce(estimand, scenario_id),
                    component, estimate, ci95_low, ci95_high, direction)),
  "",
  "![](figures/path_grid.png)",
  "",
  "Grid of (M, D_F, D_S): `tables/path_grid.csv` (sup-t bands across cells; unsupported cells shaded).",
  "",
  "## 7. Dynamics (Stage 5)",
  "",
  "Restrictions implied by the linear 5-year slope on unrestricted temperature lags 0-4 (two-way clustered Wald):",
  "",
  kable(restr %>% select(test, F_stat, df1, df2, p_value)),
  "",
  "Local projections (horizon-specific conditional projections, not causal impulse responses; sup-t band over h = 0-10):",
  "",
  "![](figures/local_projections.png)",
  "",
  kable(lp_show),
  "",
  "Almon cubic distributed lag (lags 0-10, linear temperature): cumulative temperature effect",
  "",
  kable(dl %>% filter(grepl("^cum", contrast)) %>%
          select(contrast, estimate, se, ci95_low, ci95_high)),
  "",
  "Long differences (non-overlapping windows; mean annual growth on mean climate; adaptation benchmark, does not identify speed):",
  "",
  kable(ld %>% filter(term %in% c("Tc", "Tc2", "not estimable")) %>%
          select(model, term, estimate, std_error, p_value, window, regions)),
  "",
  "## 8. Competing adaptation/mismatch model (Stage 8)",
  "",
  kable(anom_show),
  "",
  "Warming versus cooling (mirror-image paths through the cooling hinges):",
  "",
  kable(cool %>% select(direction, contrast, estimate, ci95_low, ci95_high)),
  "",
  "## 9. Robustness, inference and falsification (Stages 10 and 12)",
  "",
  "![](figures/specification_curve.png)",
  "",
  kable(spec_show),
  "",
  "The StructChange add-back samples readmit spliced-series jumps (e.g. Honduras 2001, a log jump of about +1.5 in every region); under the country-level Rademacher DGP that cluster fattens the bootstrap-t tail (95% critical value near 7), hence their wide intervals. The Huber row uses a country-cluster pairs bootstrap (499 draws) of the full projection-and-fit; its percentile interval is the only core interval excluding zero, and on the positive (fast less harmful) side.",
  "",
  "Inference variants for the confirmatory contrasts:",
  "",
  kable(inf_show),
  "",
  "Falsification:",
  "",
  kable(leads %>% select(test, F_stat, df1, df2, p_value)),
  "",
  kable(placebo %>% select(check, contrast, estimate, ci95_low, ci95_high)),
  "",
  kable(perm %>% select(-design_hash)),
  "",
  sprintf("The literal permutation (region histories within climate zone x decade) reassigns each region an unrelated donor, which destroys the within-country spatial correlation of the exposure while the growth errors stay spatially correlated; its null SD (%s) is about half the bootstrap SE (%s), so its p-value overstates the evidence. The dependence-consistent country-block time shift (each country's own rate history lagged 6-17 years) has a null SD of %s and p = %s, in line with the bootstrap intervals.",
          f2(perm$perm_sd[perm$statistic == "speed" & !perm$dependence_consistent], 3),
          f2(spd$se, 3),
          f2(perm$perm_sd[perm$statistic == "speed" & perm$dependence_consistent], 3),
          f2(perm$p_two_sided[perm$statistic == "speed" & perm$dependence_consistent], 3)),
  "",
  kable(loo),
  "",
  "Region-level one-step influence: `tables/influence_regions.csv`.",
  "",
  "Mechanisms (sectoral real value added per capita; secondary):",
  "",
  kable(sectors %>% select(outcome, contrast, estimate, ci95_low, ci95_high,
                           n_obs, n_countries)),
  "",
  "Heterogeneity of the rate response (continuous, centred moderators; BH-adjusted):",
  "",
  kable(het %>% select(moderator, term, estimate, std_error, p_value, p_bh)),
  "",
  kable(het_gaps %>% filter(grepl("speed", contrast)) %>%
          select(moderator, quantile, value, estimate, ci95_low, ci95_high)),
  "",
  "## 10. Checklist (section 12)",
  "",
  kable(checklist),
  "",
  sprintf("Headline claims are downgraded because of items %s: %s",
          paste(failed_items$item, collapse = ", "),
          paste(failed_items$note, collapse = "; ")),
  "",
  "## 11. Not run, and why",
  "",
  "- Climate products CRU TS and UDel, and fixed-2000 population weighting (Stage 10) and PWT/WDI triangulation plus the global module (Stage 11): outside the preferred data (ERA5, area weighting, DOSE).",
  "- Stage 7 episode comparison: disabled by protocol because `config/fast_slow_episodes.yml` was not predeclared before outcome access.",
  "- Spatial block bootstrap: replaced by Conley spatial HAC with temporal Bartlett lags (500/1,000/2,000 km; 5/10 years).",
  "- Polygon-level geography harmonization: requires the DOSE polygons (not in `data/`).",
  "- Non-climate negative-control outcomes and hot-day / drought channels: not available in the preferred data.",
  "- `targets`/`renv` orchestration: the pipeline is four scripts plus `tests/test_fast_slow_warming.R`; package versions are in `qa/session_info.txt`."
)
writeLines(report, fs_path("report.md"))
message("Wrote ", fs_path("report.md"))
