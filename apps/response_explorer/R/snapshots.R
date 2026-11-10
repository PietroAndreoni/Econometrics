# Snapshots of a setting's model: every input of the setting's model panels
# (the list setting_server()'s model_inputs() returns), saved with a name and
# loaded back into any setting. The user's snapshots are kept in
# results/response_explorer/snapshots.rds; the literature presets below are
# built from the app's own options (approximations where the app lacks an
# option of the paper, said in the note).

# The defaults of the model inputs, as the setting panels start.
default_model_inputs <- function(gdp = "PWT") {
  list(gdp = gdp, climate = "ERA5", outcome = "log", baseline = "20-year", style = "simple",
       weight = "population", units_T = "SD", units_P = "SD", fe = "quadratic", effects = "both",
       convergence = "none", quad_T = "none", quad_P = "none",
       form_T = "Quadratic arms", form_P = "Quadratic arms", symmetric_T = FALSE, symmetric_P = FALSE,
       use_lags = FALSE, lags = 5, inter_T = DEFAULT_MODIFIERS$inter_T, inter_P = DEFAULT_MODIFIERS$inter_P,
       low_threshold = round(EXPLORER$low_income[[gdp]]), moderators = "trailing", share_trailing = FALSE,
       income_form_T = "dollars", income_form_P = "standard", xi_T = 0.6, xi_P = 0,
       exposure_mode = "none", exposure_channels = "both", exp_hot_src = "tmpmax", exp_hot_value = 29,
       exp_hot_measure = "dd", exp_cold_src = "tmpmin", exp_cold_value = 0, exp_cold_measure = "dd",
       exposure_mode_P = "none", exposure_channels_P = "both", exp_wet_type = "pct", exp_wet_pct = "0.99",
       exp_wet_abs = "20", exp_wet_measure = "dd", exp_dry_value = "1", exp_dry_measure = "days",
       sd_mode = "fixed", clim_mode = "baseline", cap_T = FALSE, cap_P = FALSE,
       shape_T = "linear", shape_P = "linear", controls = unname(CONTROL_CHOICES),
       drop_outliers = TRUE, outlier_tail = DEFAULT_OUTLIERS)
}
# A preset: the defaults (of its GDP source) with the given inputs changed.
preset_inputs <- function(...) {
  x <- list(...)
  m <- default_model_inputs(x$gdp %||% "PWT")
  for (k in names(x)) m[k] <- list(x[[k]])
  m
}
# A plain specification of the literature: T and P themselves (absolute units)
# entering as they are, no outlier drops or controls, standard income terms.
literature_inputs <- function(...) {
  base <- list(baseline = "Level", units_T = "absolute", units_P = "absolute", symmetric_T = TRUE,
               symmetric_P = TRUE, inter_T = character(0), inter_P = character(0), income_form_T = "standard",
               income_form_P = "standard", xi_T = 0, controls = character(0), drop_outliers = FALSE,
               moderators = "fixed", clim_mode = "fixed", fe = "unit")
  x <- list(...)
  for (k in names(x)) base[k] <- list(x[[k]])
  do.call(preset_inputs, base)
}

LITERATURE_PRESETS <- list(
  list(name = "App defaults (preferred specification)",
       note = "The panels as the app starts.",
       inputs = preset_inputs()),
  list(name = "Burke, Hsiang and Miguel (2015)",
       note = paste("Growth on T, T², P and P² (levels), country and year effects, country quadratic trends;",
                    "WB GDP and UDelaware climate with population weights, as the paper."),
       inputs = literature_inputs(gdp = "WB", climate = "UDelaware", form_T = "Linear + quadratic arms",
                                  form_P = "Linear + quadratic arms", fe = "quadratic")),
  list(name = "Dell, Jones and Olken (2012), approx.",
       note = paste("Growth on T and P levels and their interaction with a poor-country dummy; country and year",
                    "effects (the paper: region x year), PWT GDP, UDelaware climate. Poor: below the app's",
                    "low-income threshold (the paper: below the median PPP GDP of the first year)."),
       inputs = literature_inputs(gdp = "PWT", climate = "UDelaware", form_T = "Linear arms", form_P = "Linear arms",
                                  inter_T = "low", inter_P = "low")),
  list(name = "Kalkuhl and Wenz (2020), national",
       note = paste("Growth on the first differences of T and P and their interaction with last year's level;",
                    "country and year effects. The paper is on subnational regions (here national panels, ERA5)."),
       inputs = literature_inputs(gdp = "WB", style = "diff", form_T = "Linear arms", form_P = "Linear arms",
                                  inter_T = "lag", inter_P = "lag")),
  list(name = "Kotz, Levermann and Wenz (2024), approx.",
       note = paste("First differences of T and P interacted with the long-run climate, with 10 lags; country",
                    "and year effects with country linear trends. The paper also has daily variability, wet days",
                    "and extreme precipitation (not in the app) and is on subnational regions."),
       inputs = literature_inputs(gdp = "WB", style = "diff", form_T = "Linear arms", form_P = "Linear arms",
                                  inter_T = "clim", inter_P = "clim", use_lags = TRUE, lags = min(10, EXPLORER$max_lags %||% 10),
                                  fe = "linear")),
  list(name = "Exposure specification of the xi analysis (2026-11-09)",
       note = paste("T: linear arms through Tmax > 29 C and Tmin < 0 C degree-day exposure (country slopes);",
                    "P: quadratic arms x climate and share; both scaled with xi = 0; country and year effects,",
                    "no trends, last year's log GDP; trailing 20-year characteristics (xi_exposure_*.R)."),
       inputs = preset_inputs(fe = "unit", convergence = "lag1", form_T = "Linear arms", units_T = "absolute",
                              inter_T = character(0), inter_P = c("clim", "agr"), income_form_T = "dollars",
                              income_form_P = "dollars", xi_T = 0, xi_P = 0, exposure_mode = "country"))
)
names(LITERATURE_PRESETS) <- vapply(LITERATURE_PRESETS, `[[`, "", "name")

SNAPSHOTS_FILE <- normalizePath(file.path(APP_DIR, "..", "..", "results", "response_explorer", "snapshots.rds"),
                                mustWork = FALSE)
read_snapshots <- function() {
  if (!file.exists(SNAPSHOTS_FILE)) return(list())
  tryCatch(readRDS(SNAPSHOTS_FILE), error = function(e) list())
}
write_snapshots <- function(s) {
  dir.create(dirname(SNAPSHOTS_FILE), showWarnings = FALSE, recursive = TRUE)
  saveRDS(s, SNAPSHOTS_FILE)
}
# The choices of the load menu: the literature presets and the saved ones.
snapshot_choices <- function(saved) {
  out <- list(Literature = stats::setNames(paste0("lit:", names(LITERATURE_PRESETS)), names(LITERATURE_PRESETS)))
  if (length(saved)) out$Saved <- stats::setNames(paste0("own:", names(saved)), names(saved))
  out
}
# The snapshot of a choice (NULL: none).
snapshot_of <- function(choice, saved) {
  if (is.null(choice) || !nzchar(choice)) return(NULL)
  key <- sub("^(lit|own):", "", choice)
  if (startsWith(choice, "lit:")) LITERATURE_PRESETS[[key]] else saved[[key]]
}

# Set the model inputs of a setting's panels (session: the setting's module)
# to m; inputs m lacks are left as they are.
apply_model_inputs <- function(session, m) {
  has <- function(k) k %in% names(m)
  for (k in c("gdp", "climate", "outcome", "baseline", "form_T", "form_P", "fe", "quad_T", "quad_P",
              "exposure_mode", "exp_hot_src", "exp_cold_src", "exposure_mode_P", "exp_wet_pct", "exp_wet_abs",
              "exp_dry_value")) {
    if (has(k) && !is.null(m[[k]])) updateSelectInput(session, k, selected = as.character(m[[k]]))
  }
  for (k in c("moderators", "effects", "convergence", "style", "units_T", "units_P", "shape_T", "shape_P",
              "income_form_T", "income_form_P", "exposure_channels", "exp_hot_measure", "exp_cold_measure",
              "exposure_channels_P", "exp_wet_type", "exp_wet_measure", "exp_dry_measure", "weight", "sd_mode",
              "clim_mode", "outlier_tail")) {
    if (has(k) && !is.null(m[[k]])) updateRadioButtons(session, k, selected = m[[k]])
  }
  for (k in c("symmetric_T", "symmetric_P", "cap_T", "cap_P", "use_lags", "share_trailing", "drop_outliers")) {
    if (has(k)) updateCheckboxInput(session, k, value = isTRUE(m[[k]]))
  }
  for (k in c("lags", "exp_hot_value", "exp_cold_value", "low_threshold")) {
    if (has(k) && is.finite(m[[k]] %||% NA)) updateNumericInput(session, k, value = m[[k]])
  }
  for (k in c("xi_T", "xi_P")) if (has(k) && is.finite(m[[k]] %||% NA)) updateSliderInput(session, k, value = m[[k]])
  for (k in c("inter_T", "inter_P", "controls")) {
    if (has(k)) updateCheckboxGroupInput(session, k, selected = as.character(m[[k]] %||% character(0)))
  }
}
