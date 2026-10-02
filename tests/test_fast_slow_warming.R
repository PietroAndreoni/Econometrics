# Synthetic fixtures for the fast-versus-slow warming pipeline (section 8.3,
# tests (a)-(g), plus checks of the rate definition and of the two-way variance
# against fixest). Run from the project root:
#   Rscript tests/test_fast_slow_warming.R

source("load_functions.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", label, call. = FALSE)
  message("ok  ", label)
}

spec0 <- list(T0 = 0, P_mean = 0, P_sd = 1)
regions <- tibble::tibble(GID_1 = "X.1_1", B = 10, P = 1)

# Rate definition -------------------------------------------------------------

set.seed(1)
yrs <- 1950:1980
tmp <- 10 + cumsum(rnorm(length(yrs), 0, 0.3))
r5 <- trailing_temperature_slope(tmp, yrs, 5)
lm_slope <- vapply(seq_along(yrs), function(i) {
  if (i < 5) return(NA_real_)
  unname(coef(lm(tmp[(i - 4):i] ~ yrs[(i - 4):i]))[2])
}, numeric(1))
check("5-year slope equals the OLS slope on the five years ending in t",
      isTRUE(all.equal(r5, lm_slope)))
check("slope is in degrees per year (10x fit_trailing_climate_velocity)",
      isTRUE(all.equal(10 * r5, fit_trailing_climate_velocity(yrs, tmp, 5))))
gap_years <- yrs[-10]
r_gap <- trailing_temperature_slope(tmp[-10], gap_years, 5)
# Dropping 1959 leaves windows ending 1960-1963 incomplete (the 1964 window
# starts in 1960), plus the four initial years.
check("a missing year makes every window containing it NA",
      sum(is.na(r_gap)) == 8 && all(is.na(r_gap[10:13])) &&
        !is.na(r_gap[9]) && !is.na(r_gap[14]))
w <- trailing_temperature_slope(c(0, 0, 0, 0, 1), 1:5, 5)[5]
check("weight on year t is 2/10 for k = 5", isTRUE(all.equal(w, 0.2)))

# (e) no temperature after year t ------------------------------------------------
tmp_future <- tmp
tmp_future[21:length(tmp)] <- tmp_future[21:length(tmp)] + 5
r5_future <- trailing_temperature_slope(tmp_future, yrs, 5)
check("(e) slopes never use temperature after year t",
      isTRUE(all.equal(r5[1:20], r5_future[1:20])))
check("(e) the future-only negative control uses only t+1..t+5",
      isTRUE(all.equal(future_temperature_slope(tmp, yrs, 5)[1:21],
                       trailing_temperature_slope(tmp, yrs, 5)[6:26])))

# (f) invariance to row order ----------------------------------------------------
panel <- tidyr::crossing(GID_1 = c("A.1_1", "A.2_1", "B.1_1"), year = yrs) %>%
  mutate(TM = 10 + rnorm(n()), RR = 1 + rnorm(n(), 0, 0.1))
a <- fs_add_climate_terms(panel, spec0) %>% arrange(GID_1, year)
b <- fs_add_climate_terms(panel[sample(nrow(panel)), ], spec0) %>%
  arrange(GID_1, year)
check("(f) climate terms are invariant to input row order",
      isTRUE(all.equal(as.data.frame(a), as.data.frame(b))))
check("(assertion 8) terms refuse data that carry an outcome",
      inherits(try(fs_add_climate_terms(mutate(panel, dlgrp_pc_usd = 0),
                                        spec0), silent = TRUE), "try-error"))

# Paths --------------------------------------------------------------------------
paths <- fs_path_terms(regions, spec0, total_warming = 1, fast_years = 5,
                       slow_years = 20, horizon = 40)
fast <- paths %>% filter(path == "fast")
slow <- paths %>% filter(path == "slow")
check("paths share start and endpoint",
      isTRUE(all.equal(fast$TM[fast$s >= 20], slow$TM[slow$s >= 20])) &&
        abs(fast$TM[fast$s == 5] - 11) < 1e-12)
check("the path slope is the production slope",
      isTRUE(all.equal(fast$rate[fast$s == 3],
                       trailing_temperature_slope(c(10, 10, 10.2, 10.4, 10.6),
                                                  1:5, 5)[5])))

# The template shortcut equals building every region's path separately.
two <- tibble::tibble(GID_1 = c("X.1_1", "Y.1_1"), B = c(10, 25),
                      P = c(1, 0.3))
spec_knots <- list(T0 = 14, P_mean = 1, P_sd = 0.8,
                   spline_knots = c(5, 10, 15, 20), spline_boundary = c(-5, 30))
fast_way <- fs_path_terms(two, spec_knots, 1, 5, 20, 30) %>%
  arrange(GID_1, path, s)
long <- bind_rows(lapply(c(fast = 5, slow = 20), function(D) {
  tidyr::crossing(two, fs_path_temperature(0, 1, D, 30, 30))
}), .id = "path") %>%
  mutate(TM = TM + B, RR = P, year = s + 3000L,
         unit = paste(GID_1, path))
slow_way <- fs_add_climate_terms(long %>% select(unit, GID_1, path, s, year,
                                                 TM, RR),
                                 spec_knots, id = "unit") %>%
  filter(s >= 1L) %>%
  arrange(GID_1, path, s)
shared <- intersect(setdiff(names(slow_way), c("unit", "B", "P")),
                    names(fast_way))
check("path template shortcut equals the per-region construction",
      isTRUE(all.equal(as.data.frame(fast_way[, shared]),
                       as.data.frame(slow_way[, shared]),
                       check.attributes = FALSE)))

C <- fs_path_contrast(paths, c("d1", "rate", "r_pos2", "Tc"),
                      horizons = c(20, 24, 30))
# (a) a linear annual-change effect: sum_s dT_s = M on both paths.
check("(a) a linear one-year-change effect gives no fast-slow gap",
      all(abs(C[, "d1"]) < 1e-12))

# (b) linear rolling slope: hand-calculated filter, zero after the washout.
D <- fast$TM - slow$TM
D_full <- c(rep(0, 10), D)            # s = -9..40, D = 0 before the ramp
hand <- vapply(c(20, 24, 30), function(h) {
  sum(vapply(1:h, function(s) {
    idx <- s + 10 - (4:0)            # D_{s-4}, ..., D_s
    sum(((0:4) - 2) / 10 * D_full[idx])
  }, numeric(1)))
}, numeric(1))
check("(b) linear rolling-slope gap matches the hand-calculated filter",
      isTRUE(all.equal(unname(C[, "rate"]), hand)))
check("(b) it is nonzero at the slow endpoint and zero after K-1 years",
      abs(C["h20", "rate"]) > 1e-3 && abs(C["h24", "rate"]) < 1e-12 &&
        abs(C["h30", "rate"]) < 1e-12)

# (c) negative curvature in positive rates makes the fast ramp worse.
check("(c) the fast ramp accumulates more squared warming speed",
      all(C[, "r_pos2"] > 0))
check("(c) so beta_2+ < 0 gives a negative (fast-worse) speed gap",
      all(C[, "r_pos2"] * -0.5 < 0))

# (d) a level-only DGP creates no incremental rate penalty -----------------------
set.seed(2)
sim <- tidyr::crossing(GID_1 = sprintf("C%02d.%d_1", 1:30 %/% 3, 1:30),
                       year = 1940:2000) %>%
  mutate(GID_0 = substr(GID_1, 1, 3),
         TM = 12 + rnorm(n(), 0, 1) + 0.02 * (year - 1940), RR = 1)
spec_sim <- list(T0 = 12, P_mean = 1, P_sd = 1)
sim <- fs_add_climate_terms(sim, spec_sim) %>%
  filter(!is.na(rate)) %>%
  mutate(g = -0.4 * Tc - 0.05 * Tc2, year_c = year - 1970)
fit_d <- fixest::feols(g ~ Tc + Tc2 + r_pos + r_pos2 + r_neg + r_neg2 |
                         GID_1 + year, sim, notes = FALSE)
paths_sim <- fs_path_terms(tibble::tibble(GID_1 = "C", B = 12, P = 1),
                           spec_sim, 1, 5, 20, 30)
C_sim <- fs_path_contrast(paths_sim, names(coef(fit_d)), 20)
speed <- sum(C_sim[, fs_rate_terms()] * coef(fit_d)[fs_rate_terms()])
check("(d) level-only DGP: incremental speed gap is zero",
      abs(speed) < 1e-8)
total <- sum(C_sim * coef(fit_d)[colnames(C_sim)])
check("(d) ... while the total gap is the level-timing gap (nonzero)",
      total < -1e-3)

# (g) decision-rule boundaries ----------------------------------------------------
cls <- fs_classify(
  ci95_low  = c(-2,  -2,     -3,   -3,  1.5, -0.5, -1.2, -1),
  ci95_high = c( 0,  -1e-9,  -1,   -1.01, 3,  0.5,  0.9,  -0.2),
  ci90_low  = c(-1.8, -1.8,  -2.8, -2.8, 1.6, -0.4, -1,  -0.9),
  ci90_high = c(-0.2, -0.1,  -1.2, -1.2, 2.9,  0.4,  1,  -0.3),
  sesoi = 1, sesoi_confirmed = TRUE
)
check("(g) upper 95% bound exactly 0 is unresolved", cls$direction[1] ==
        "direction unresolved")
check("(g) upper bound just below 0 is fast worse", cls$direction[2] ==
        "fast worse")
check("(g) upper bound exactly -SESOI is not material", cls$materiality[3] ==
        "magnitude unresolved")
check("(g) upper bound below -SESOI is materially worse",
      cls$materiality[4] == "materially fast worse")
check("(g) lower bound above +SESOI is materially less harmful",
      cls$direction[5] == "fast less harmful" &&
        cls$materiality[5] == "materially fast less harmful")
check("(g) 90% interval inside the bounds is equivalent",
      cls$materiality[6] == "practically equivalent")
check("(g) 90% interval exactly on the bounds is equivalent (closed)",
      cls$materiality[7] == "practically equivalent")
check("(g) directional but practically equivalent is allowed",
      cls$direction[8] == "fast worse" &&
        cls$materiality[8] == "practically equivalent")
cls_unconfirmed <- fs_classify(-0.5, 0.5, -0.4, 0.4, 1)
check("(g) without an affirmed SESOI, equivalence is only provisional",
      grepl("^provisional", cls_unconfirmed$materiality))
cls_unsupported <- fs_classify(-3, -2, -2.8, -2.2, 1, supported = FALSE)
check("(g) an unsupported path is not identified",
      cls_unsupported$direction == "not identified for this path")

# Two-way variance and bootstrap machinery ------------------------------------
set.seed(3)
sim$g_noise <- sim$g + rnorm(nrow(sim)) +
  rnorm(dplyr::n_distinct(sim$year))[match(sim$year, sort(unique(sim$year)))]
fit_v <- fixest::feols(g_noise ~ Tc + Tc2 + r_pos + r_pos2 | GID_1 + year +
                         GID_0[[year_c]], sim, vcov = ~ GID_0 + year,
                       notes = FALSE)
fwl <- fs_fwl(fit_v, sim)
V_mine <- fs_cgm_vcov(fwl, adjust = FALSE)
V_fixest <- vcov(fit_v, vcov = ~ GID_0 + year,
                 ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE))
check("two-way CGM variance matches fixest (no small-sample factors)",
      isTRUE(all.equal(unname(V_mine), unname(V_fixest), tolerance = 1e-6,
                       check.attributes = FALSE)))
prep <- fs_boot_prepare(fwl)
V_hat <- fs_cgm_vcov(fwl)
v_one <- matrix(1, 1, dplyr::n_distinct(sim$GID_0),
                dimnames = list(NULL, sort(unique(sim$GID_0))))
draw <- fs_boot_draws(prep, v_one, "GID_0")
check("a bootstrap draw with all weights +1 reproduces beta_hat and V_hat",
      max(abs(draw$delta)) < 1e-8 &&
        isTRUE(all.equal(as.numeric(V_hat), draw$Vflat[1, ], tolerance = 1e-6)))

message("All fast-slow warming tests passed.")
