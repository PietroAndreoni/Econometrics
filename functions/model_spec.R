# Formula builder and panel estimator shared by every climate-response model.
# Defaults are the main subnational specification: region-specific linear
# trends plus year fixed effects, errors clustered by region.

make_panel_formula <- function(
    terms,
    outcome = "dlgrp_pc_usd",
    fixed_effects = "year + GID_1[year]"
) {
  stats::as.formula(paste(outcome, "~", terms, "|", fixed_effects))
}

fit_panel_model <- function(
    terms,
    model_data,
    outcome = "dlgrp_pc_usd",
    fixed_effects = "year + GID_1[year]",
    panel_id = c("GID_1", "year"),
    cluster = ~GID_1,
    ...
) {
  fixest::feols(
    make_panel_formula(terms, outcome, fixed_effects),
    data = model_data,
    panel.id = panel_id,
    cluster = cluster,
    ...
  )
}
