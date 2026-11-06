# OECD regions (oecd_subnational.R's sample and construction): the country's
# income as a moderator on top of the region's. Year effects only (country x
# year effects would absorb the country's income). Saved to
# results/oecd_country_income/ (oecd_country_income.rds, .log).
#
# Income moderators, both moving (log mean GDP per capita of the 5 years
# before, less its mean over the rows):
#   y      the region's (DOSE, constant 2015 US$), as in oecd_subnational.R
#   y_cty  the country's (World Bank, constant 2015 US$; the handout's
#          trailing income of the national panel)
# Check: the exact split of the region's income into the country-year mean
# over its regions (y_btw) and the region's deviation from it (y_wth).
#
#   C     arms x (climate, y)               the region's income
#   Cc    arms x (climate, y_cty)           the country's income
#   D     arms x (climate, y, y_cty)        both
#   Dw    arms x (climate, y_wth, y_btw)    both, from the regional data
# Every model with the level terms of its moderators. Errors clustered by
# country (region clusters alongside).
#
# Needs what oecd_subnational.R needs. From the project root:
#   Rscript oecd_country_income.R

# oecd_subnational.R's setup: the rows, deviations, moderators and helpers
# (without its sink and its models).
src <- readLines("oecd_subnational.R")
eval(parse(text = src[grep("^bm <- readLines", src):(grep("^sub <- run\\(", src) - 1)]))
out_dir <- file.path("results", "oecd_country_income")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "oecd_country_income.log"), split = TRUE)

# The country's trailing income: the World Bank national panel's.
wb <- national_panels[[match("WB x ERA5", panel_labels)]]
cty <- trailing_income(wb, TRAILING_INCOME_YEARS)
print(head(as.data.frame(cty[is.finite(cty$inc_trailing), ]), 3))
cty <- transmute(cty, GID_0 = substr(GID_1, 1, 3), year, inc_cty = inc_trailing)
stopifnot(!anyDuplicated(cty[c("GID_0", "year")]))
n0 <- nrow(d)
d <- left_join(d, cty, by = c("GID_0", "year"))
d <- d[is.finite(d$inc_cty), ]
cat(sprintf("Rows with the country's trailing income: %d of %d; %d regions, %d countries\n",
            nrow(d), n0, n_distinct(d$GID_1), n_distinct(d$GID_0)))
d$y <- d$inc_trailing - mean(d$inc_trailing)
d$y_cty <- d$inc_cty - mean(d$inc_cty)
d$y_btw <- ave(d$y, d$GID_0, d$year)
d$y_wth <- d$y - d$y_btw

# The income moderators: how they relate.
cat(sprintf("Correlation of the region's and the country's income: %.3f (region's with the country-year mean of the regions': %.3f)\n",
            cor(d$y, d$y_cty), cor(d$y, d$y_btw)))
cat(sprintf("Variance of the region's income within the country-year: %.3f of the total (SD %.3f log points; between SD %.3f)\n",
            var(d$y_wth) / var(d$y), sd(d$y_wth), sd(d$y_btw)))
cat(sprintf("Region's income net of the country's (WB): residual SD %.3f, R2 %.3f\n",
            sd(resid(lm(y ~ y_cty, d))), summary(lm(y ~ y_cty, d))$r.squared))

# Cl: C with the level of the country's income too (D's level terms), so that
# D against Cl is the country-income interactions alone.
MODS <- list(C = "y", Cl = "y", Cc = "y_cty", D = c("y", "y_cty"), Dw = c("y_wth", "y_btw"))
LEVEL_MODS <- modifyList(MODS, list(Cl = c("y", "y_cty")))
SPECS2 <- lapply(MODS, function(z) c(arm_terms("TM", c(CLIM[["TM"]], z)), arm_terms("RR", c(CLIM[["RR"]], z))))
fit2 <- function(k, data, cluster = ~GID_0) {
  terms <- SPECS2[[k]]
  for (t in names(terms)) data[[t]] <- Reduce(`*`, lapply(terms[[t]], function(z) data[[z]]))
  levels <- c("clim_TM", "clim_RR", LEVEL_MODS[[k]])
  f <- as.formula(paste("g ~", paste(c(names(terms), levels, CONTROLS), collapse = " + "), "| year + GID_1[t, t2]"))
  feols(f, data = data, cluster = cluster, panel.id = ~GID_1 + year, notes = FALSE)
}
fits <- lapply(setNames(names(SPECS2), names(SPECS2)), fit2, data = d)
G <- n_distinct(d$GID_0)
ARMS <- c(pTM = "Warm", nTM = "Cold", pRR = "Wet", nRR = "Dry")
ROWS <- unlist(lapply(names(ARMS), function(a) c(a, paste0(a, "_X_", c(CLIM[[substr(a, 2, 3)]], "y", "y_cty", "y_wth", "y_btw")))))
row_label <- function(r) {
  a <- substr(r, 1, 3)
  z <- sub("^...(_X_)?", "", r)
  paste(ARMS[[a]], c(clim_TM = "x climate", clim_RR = "x climate", y = "x region income", y_cty = "x country income",
                     y_wth = "x region income, within", y_btw = "x region income, country-year mean")[z] %|% "arm")
}
`%|%` <- function(a, b) ifelse(is.na(a), b, a)

for (cl in c(country = "GID_0", region = "GID_1")) {
  Gc <- n_distinct(d[[cl]])
  cells <- sapply(fits, function(m) {
    V <- vcov(m, cluster = as.formula(paste("~", cl)))
    b <- coef(m)
    s <- sqrt(diag(V))
    p <- 2 * pt(-abs(b / s), Gc - 1)
    ifelse(ROWS %in% names(b), sprintf("%+.3f (%.3f)%s", 100 * b[ROWS], 100 * s[ROWS], star(p[ROWS])), "")
  })
  rownames(cells) <- vapply(ROWS, row_label, "")
  cells <- cells[rowSums(cells != "") > 0, ]
  options(width = 220)
  cat(sprintf("\nCoefficients, pp of growth per SD (x 100); per log point of income for the income terms; SE clustered by %s (%d)\n",
              names(which(c(country = "GID_0", region = "GID_1") == cl)), Gc))
  print(noquote(cells))
  tests <- sapply(names(fits), function(k) {
    m <- fits[[k]]
    V <- vcov(m, cluster = as.formula(paste("~", cl)))
    sel <- function(z) paste0(names(ARMS), "_X_", z)
    one <- function(z) if (all(sel(z) %in% names(coef(m)))) sprintf("%.3f", wald(m, V, Gc, sel(z))) else ""
    c(`Region income (4 arms)` = one("y"), `Country income (4 arms)` = one("y_cty"),
      `Within-country income` = one("y_wth"), `Between (country-year mean)` = one("y_btw"),
      `T x country income` = if ("pTM_X_y_cty" %in% names(coef(m))) sprintf("%.3f", wald(m, V, Gc, c("pTM_X_y_cty", "nTM_X_y_cty"))) else "",
      `T x region income` = if ("pTM_X_y" %in% names(coef(m))) sprintf("%.3f", wald(m, V, Gc, c("pTM_X_y", "nTM_X_y"))) else "")
  })
  cat("Joint Wald tests, p\n")
  print(noquote(tests))
}
cat("\nFit\n")
print(noquote(rbind(`Within R2` = sapply(fits, function(m) sprintf("%.5f", fitstat(m, "wr2")[[1]])),
                    AIC = sapply(fits, function(m) sprintf("%.1f", AIC(m))),
                    BIC = sapply(fits, function(m) sprintf("%.1f", BIC(m))),
                    Rows = sapply(fits, nobs))))

# How much of the variation of the response to a shock each moderator carries:
# for every row, the response of an arm (pp per SD) implied by the model,
# b0 + b_clim clim + b_y y + b_c y_cty; the SD over the rows of each part and
# the share of the variance of the response from each (each part's variance
# plus half the covariances with the others: the shares add up to 1).
heterogeneity <- function(m, mods) {
  out <- list()
  for (a in names(ARMS)) {
    zs <- c(CLIM[[substr(a, 2, 3)]], mods)
    parts <- sapply(zs, function(z) 100 * coef(m)[[paste0(a, "_X_", z)]] * d[[z]])
    total <- rowSums(parts)
    C <- cov(parts)
    share <- (diag(C) + (rowSums(C) - diag(C))) / var(total)
    r <- 100 * coef(m)[[a]] + total
    out[[a]] <- data.frame(arm = ARMS[[a]], moderator = zs, sd_part = round(apply(parts, 2, sd), 3),
                           share = round(share, 2), sd_response = round(sd(total), 3),
                           p10_p90 = sprintf("%+.2f / %+.2f", quantile(r, 0.1), quantile(r, 0.9)))
  }
  do.call(rbind, out)
}
het <- lapply(setNames(names(MODS), names(MODS)), function(k) heterogeneity(fits[[k]], MODS[[k]]))
for (k in names(het)) {
  cat(sprintf("\nModel %s: variation of the response over the rows (pp per SD): SD of each moderator's part, its share of the variance of the response; SD and p10/p90 of the response\n", k))
  print(het[[k]], row.names = FALSE)
}

# Model D without one country at a time (coefficient x 100, t with country
# clusters): the income terms.
LOO <- c("pTM_X_y", "pTM_X_y_cty", "nTM_X_y", "nTM_X_y_cty", "pRR_X_y", "pRR_X_y_cty", "nRR_X_y", "nRR_X_y_cty")
loo <- t(sapply(c("none", sort(unique(d$GID_0))), function(cc) {
  m <- fit2("D", d[d$GID_0 != cc, ])
  sprintf("%+.2f (%+.1f)", 100 * coef(m)[LOO], coef(m)[LOO] / se(m)[LOO])
}))
colnames(loo) <- vapply(LOO, row_label, "")
cat("\nModel D without one country (coefficient x 100, t with country clusters)\n")
print(noquote(loo))
LOO_W <- c("pTM_X_y_wth", "pTM_X_y_btw", "nTM_X_y_wth", "nTM_X_y_btw", "pRR_X_y_wth", "pRR_X_y_btw", "nRR_X_y_wth", "nRR_X_y_btw")
loo_w <- t(sapply(c("none", sort(unique(d$GID_0))), function(cc) {
  m <- fit2("Dw", d[d$GID_0 != cc, ])
  sprintf("%+.2f (%+.1f)", 100 * coef(m)[LOO_W], coef(m)[LOO_W] / se(m)[LOO_W])
}))
colnames(loo_w) <- sub("x region income, ", "", vapply(LOO_W, row_label, ""))
cat("\nModel Dw without one country (coefficient x 100, t with country clusters)\n")
print(noquote(loo_w))

saveRDS(list(coef = lapply(fits, coef), vcov = lapply(fits, vcov, cluster = ~GID_0), het = het, loo = loo, loo_w = loo_w,
             countries = sort(unique(d$GID_0)), rows = nrow(d)),
        file.path(out_dir, "oecd_country_income.rds"))
sink()
