# Inference for the main regression beyond the country-clustered Wald tests:
# the preferred specification (xi = 0) and M2b on the six GDP x climate
# datasets (the handout's sample: 5% non-weather outliers, each source's own
# flags), every weather coefficient and the joint tests of the arms, variables
# and interactions, with
#   CR1     the handout's clustered (country) Wald tests, F(q, G - 1);
#   WCR     wild cluster restricted bootstrap (country; Webb weights; the null
#           imposed; CR1 Wald statistic), 9999 draws;
#   CR3     jackknife variance: leave one country out, exact,
#           (G - 1) / G sum (b_-g - b)(b_-g - b)';
#   2-way   clustered by country and year (fixest), F(q, min(G, years) - 1);
#   DK      Driscoll-Kraay (fixest, its default lag), F(q, years - 1);
#   Conley  country clusters (any serial correlation) plus spatial correlation
#           between countries within the year (Conley-Hsiang): Bartlett kernel
#           in the distance between country centroids (largest polygon,
#           Natural Earth), cutoffs 1000, 2000 and 3000 km; F(q, G - 1).
# And per coefficient the effective number of clusters G*(0) (Carter, Schnepel
# and Steigerwald 2017) and the countries with the largest partial leverage.
# All computed on the model's design with the country effects and quadratic
# trends partialled out (nested in the clusters, so the leave-one-out and the
# bootstrap are exact) and the year effects as dummies.
#   Rscript handout_inference.R
# Results in results/handout_preferred_specification/inference.rds.

suppressPackageStartupMessages(library(fixest))
`%||%` <- function(a, b) if (is.null(a)) b else a
B <- as.numeric(Sys.getenv("INFERENCE_B", "9999"))
SEED <- 20261103
CUTOFFS <- c(1000, 2000, 3000)
DATASETS <- expand.grid(climate = c("ERA5", "CRU TS", "UDelaware"), gdp = c("PWT", "WB"), stringsAsFactors = FALSE)

# The handout's setup up to the scaled model, with the climate product free.
main <- readLines("handout_preferred_specification.R")
code <- main[seq_len(grep("^gn_fit <- function", main) - 1)]
hit <- grepl('climate = "ERA5"', code, fixed = TRUE)
stopifnot(sum(hit) == 1)
code[hit] <- sub('climate = "ERA5"', "climate = CLIMATE", code[hit], fixed = TRUE)

# Country centroids: the centroid of the largest polygon (overseas territories
# would move France's or the US's).
suppressPackageStartupMessages(library(sf))
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
iso <- if ("iso_a3_eh" %in% names(world)) world$iso_a3_eh else world$iso_a3
iso[iso == "-99"] <- world$adm0_a3[iso == "-99"]
world <- st_make_valid(st_sf(iso = iso, geometry = st_geometry(world)))
parts <- suppressWarnings(st_cast(st_cast(world, "MULTIPOLYGON"), "POLYGON"))
parts$area <- as.numeric(st_area(parts))
parts <- parts[order(-parts$area), ]
parts <- parts[!duplicated(parts$iso), ]
centroids <- st_coordinates(suppressWarnings(st_centroid(parts)))
rownames(centroids) <- parts$iso
haversine <- function(a, b) {
  r <- pi / 180
  h <- outer(a[, 2], b[, 2], function(x, y) sin((y - x) * r / 2)^2) +
    outer(cos(a[, 2] * r), cos(b[, 2] * r)) * outer(a[, 1], b[, 1], function(x, y) sin((y - x) * r / 2)^2)
  2 * 6371 * asin(pmin(sqrt(h), 1))
}

webb <- function(n) sample(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)), n, replace = TRUE)
wald <- function(b, V, df) {
  q <- length(b)
  stat <- if (q == 1) b^2 / V[1, 1] else as.numeric(t(b) %*% solve(V, b)) / q
  pf(stat, q, df, lower.tail = FALSE)
}

# The design of a model (terms: name = product of columns of data) with the
# country effects and trends partialled out and the year dummies; checked
# against the handout's fit.
design <- function(env, terms, data) {
  for (k in names(terms)) data[[k]] <- Reduce(`*`, lapply(terms[[k]], function(z) data[[z]]))
  rhs <- c(names(terms), env$LEVELS, env$CONTROLS)
  m <- feols(as.formula(paste("g ~", paste(rhs, collapse = " + "), "|", env$FE)), data = data,
             cluster = ~GID_0, panel.id = ~GID_1 + year, notes = FALSE)
  m1 <- feols(as.formula(paste("g ~", paste(c(rhs, "i(year)"), collapse = " + "), "| GID_1[t, t2]")),
              data = data, cluster = ~GID_0, demeaned = TRUE, notes = FALSE)
  stopifnot(nobs(m) == nobs(m1), isTRUE(all.equal(coef(m)[names(terms)], coef(m1)[names(terms)], tolerance = 1e-8)))
  X <- m1$X_demeaned[, names(coef(m1)), drop = FALSE]
  rows <- obs(m1)
  list(m = m, data = data[rows, ], X = X, y = as.numeric(m1$y_demeaned), cluster = data$GID_0[rows],
       year = data$year[rows], terms = names(terms), f = formula(m))
}

# Joint tests: each arm, each variable, each interaction of a variable.
joint_sets <- function(terms) {
  is_v <- function(v) grepl(paste0("^[pn]", v), terms)
  has <- function(z) grepl(paste0("_X_", z, "($|_X_)"), terms)
  sets <- list(`T terms` = terms[is_v("TM")], `P terms` = terms[is_v("RR")],
               `warm arm` = terms[grepl("^pTM", terms)], `cold arm` = terms[grepl("^nTM", terms)],
               `wet arm` = terms[grepl("^pRR", terms)], `dry arm` = terms[grepl("^nRR", terms)],
               `T x climate` = terms[is_v("TM") & has("clim_TM")], `T x income` = terms[is_v("TM") & has("y")],
               `P x climate` = terms[is_v("RR") & has("clim_RR")], `P x share` = terms[is_v("RR") & has("agr")])
  sets[lengths(sets) > 0]
}

inference <- function(env, terms, data) {
  ds <- design(env, terms, data)
  X <- ds$X; y <- ds$y; cl <- ds$cluster
  # Columns collinear to within fixest's tolerance but not to solve's (dropped; never weather terms).
  qx <- qr(X, tol = 1e-9)
  if (qx$rank < ncol(X)) {
    alias <- colnames(X)[qx$pivot[(qx$rank + 1):ncol(X)]]
    stopifnot(!any(ds$terms %in% alias))
    message("dropped as collinear: ", paste(alias, collapse = ", "))
    X <- X[, setdiff(colnames(X), alias), drop = FALSE]
  }
  N <- nrow(X); K <- ncol(X)
  idx <- split(seq_len(N), cl)
  G <- length(idx); years <- length(unique(ds$year))
  A <- solve(crossprod(X))
  b <- drop(A %*% crossprod(X, y))
  u <- y - drop(X %*% b)
  cc <- G / (G - 1) * (N - 1) / (N - K)
  S <- vapply(idx, function(i) drop(crossprod(X[i, , drop = FALSE], u[i])), numeric(K))
  V <- list(CR1 = cc * A %*% tcrossprod(S) %*% A)
  dimnames(V$CR1) <- list(colnames(X), colnames(X))
  # fixest counts the country trends in K: its constant, so that CR1 is the handout's.
  ratio <- diag(V$CR1)[ds$terms] / diag(vcov(ds$m))[ds$terms]
  stopifnot(diff(range(ratio)) < 1e-6)
  cc <- cc / ratio[[1]]
  V$CR1 <- V$CR1 / ratio[[1]]

  # CR3: b_-g - b = -A X_g' (I - X_g A X_g')^-1 u_g, checked by refitting
  # without the country with the most rows. A country whose controls no other
  # country has (the collapse dummies): refit without it and without those
  # columns (zero once it is out; their change is set to 0).
  Dg <- vapply(idx, function(i) {
    Xi <- X[i, , drop = FALSE]
    out <- tryCatch(-drop(A %*% crossprod(Xi, solve(diag(length(i)) - Xi %*% A %*% t(Xi), u[i]))), error = function(e) NULL)
    if (!is.null(out)) return(out)
    keep <- colSums(X[-i, , drop = FALSE]^2) > 1e-12
    stopifnot(!any(ds$terms %in% colnames(X)[!keep]))
    out <- numeric(K)
    out[keep] <- qr.coef(qr(X[-i, keep, drop = FALSE]), y[-i]) - b[keep]
    out
  }, numeric(K))
  big <- names(which.max(lengths(idx)))
  loo <- feols(ds$f, data = ds$data[ds$data$GID_0 != big, ], notes = FALSE)
  stopifnot(isTRUE(all.equal(unname(coef(loo)[ds$terms]), unname(b[ds$terms] + Dg[ds$terms, big]), tolerance = 1e-6)))
  V$CR3 <- (G - 1) / G * tcrossprod(Dg)
  dimnames(V$CR3) <- dimnames(V$CR1)

  V$`2-way` <- vcov(ds$m, vcov = ~GID_0 + year)
  V$DK <- vcov(ds$m, vcov = "DK")
  # Conley-Hsiang: scores A x_i u_i; the country clusters plus the pairs of
  # different countries in the same year, weighted by the kernel.
  s <- (X * u) %*% A
  within <- crossprod(rowsum(s, cl))
  miss <- setdiff(names(idx), rownames(centroids))
  if (length(miss)) stop("no centroid for ", paste(miss, collapse = ", "))
  dist <- haversine(centroids[names(idx), ], centroids[names(idx), ])
  dimnames(dist) <- list(names(idx), names(idx))
  by_year <- split(seq_len(N), ds$year)
  for (cut in CUTOFFS) {
    Vs <- within
    for (i in by_year) {
      k <- pmax(1 - dist[cl[i], cl[i], drop = FALSE] / cut, 0)
      diag(k) <- 0
      Vs <- Vs + crossprod(s[i, , drop = FALSE], k %*% s[i, , drop = FALSE])
    }
    V[[paste0("Conley ", cut, " km")]] <- cc * Vs
    dimnames(V[[paste0("Conley ", cut, " km")]]) <- dimnames(V$CR1)
  }
  df <- c(CR1 = G - 1, CR3 = G - 1, `2-way` = min(G, years) - 1, DK = years - 1,
          setNames(rep(G - 1, length(CUTOFFS)), paste0("Conley ", CUTOFFS, " km")))

  # The WCR bootstrap of R b = 0 for the coefficients J: the restricted fit
  # (J dropped), y* = X b_r + u_r v_g; b* = b_r + A S_r v and the clustered
  # scores of the bootstrap fit v_g s_g - X_g'X_g A S_r v, so that only q x G
  # pieces per country are needed for the q x q CR1 matrix of every draw.
  set.seed(SEED)
  v <- matrix(webb(G * B), G, B)
  wcr <- function(J) {
    q <- length(J)
    Xr <- X[, -J, drop = FALSE]
    ur <- y - drop(Xr %*% solve(crossprod(Xr), crossprod(Xr, y)))
    Sr <- vapply(idx, function(i) drop(crossprod(X[i, , drop = FALSE], ur[i])), numeric(K))
    M <- A[J, , drop = FALSE]
    a <- M %*% Sr
    AS <- A %*% Sr
    num <- a %*% v
    meat <- array(0, c(q, q, B))
    for (g in seq_len(G)) {
      Xi <- X[idx[[g]], , drop = FALSE]
      mg <- a[, g] %o% v[g, ] - (M %*% crossprod(Xi) %*% AS) %*% v
      for (p in seq_len(q)) for (r in p:q) meat[p, r, ] <- meat[p, r, ] + mg[p, ] * mg[r, ]
    }
    stat_b <- if (q == 1) num[1, ]^2 / (cc * meat[1, 1, ]) else vapply(seq_len(B), function(j) {
      w <- meat[, , j]; w[lower.tri(w)] <- t(w)[lower.tri(w)]
      as.numeric(t(num[, j]) %*% solve(cc * w, num[, j]))
    }, 0)
    stat <- as.numeric(t(b[J]) %*% solve(V$CR1[J, J, drop = FALSE], b[J]))
    mean(stat_b >= stat)
  }

  # Per coefficient: G*(0) from the partial leverage of the countries (with no
  # within-country correlation, gamma_g is the country's share of the
  # coefficient's sum of squared influence) and the most influential ones.
  W <- (X %*% A)[, ds$terms, drop = FALSE]
  lev <- rowsum(W^2, cl)
  lev <- sweep(lev, 2, colSums(lev), "/")
  gstar <- apply(lev, 2, function(l) G / (1 + mean((l - mean(l))^2) / mean(l)^2))
  top <- apply(lev, 2, function(l) {
    o <- order(-l)[1:5]
    paste(sprintf("%s %.0f%%", rownames(lev)[o], 100 * l[o]), collapse = ", ")
  })

  sets <- c(setNames(as.list(ds$terms), ds$terms), joint_sets(ds$terms))
  p <- t(vapply(sets, function(J) {
    j <- match(J, colnames(X))
    c(vapply(names(V), function(k) wald(b[J], V[[k]][J, J, drop = FALSE], df[[k]]), 0), WCR = wcr(j))
  }, numeric(length(V) + 1)))
  se <- vapply(V, function(x) sqrt(diag(x))[ds$terms], numeric(length(ds$terms)))
  list(coef = data.frame(term = gsub("_X_", ":", ds$terms), estimate = b[ds$terms], se, gstar = gstar,
                         top_leverage = top, row.names = NULL, check.names = FALSE),
       p = data.frame(test = gsub("_X_", ":", names(sets)), q = lengths(sets), p, row.names = NULL, check.names = FALSE),
       G = G, N = N, years = years)
}

run <- function(gdp, climate) {
  env <- new.env()
  env$CLIMATE <- climate
  Sys.setenv(HANDOUT_GDP = gdp, HANDOUT_ROWS = "own", HANDOUT_OUTLIERS = "5", HANDOUT_MIN_INCOME = "0")
  eval(parse(text = code), envir = env)
  m2b_terms <- c(env$arm_terms("TM", c(env$CLIM[["TM"]], "y")), env$arm_terms("RR", c(env$CLIM[["RR"]], "agr")))
  list(Pref = inference(env, env$scaled_terms, env$scaled_data(0)),
       M2b = inference(env, m2b_terms, env$d))
}

out <- list()
# INFERENCE_DATASETS: a subset of the rows of DATASETS (a quick check; not saved).
which_sets <- seq_len(nrow(DATASETS))
if (nzchar(Sys.getenv("INFERENCE_DATASETS"))) which_sets <- as.integer(strsplit(Sys.getenv("INFERENCE_DATASETS"), ",")[[1]])
for (k in which_sets) {
  key <- paste(DATASETS$gdp[k], sub("UDelaware", "UDel", sub("CRU TS", "CRU", DATASETS$climate[k])))
  t0 <- Sys.time()
  out[[key]] <- run(DATASETS$gdp[k], DATASETS$climate[k])
  message(key, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
out_dir <- file.path("results", "handout_preferred_specification")
if (length(which_sets) == nrow(DATASETS)) saveRDS(list(results = out, B = B, seed = SEED, cutoffs = CUTOFFS), file.path(out_dir, "inference.rds"))

options(width = 200)
for (key in names(out)) for (mod in names(out[[key]])) {
  r <- out[[key]][[mod]]
  cat(sprintf("\n==== %s, %s: %d countries, %d rows, %d years\n", key, mod, r$G, r$N, r$years))
  pr <- r$p
  pr[-(1:2)] <- lapply(pr[-(1:2)], function(x) sprintf("%.3f", x))
  print(pr, row.names = FALSE)
  cf <- r$coef
  cols <- setdiff(names(cf), c("term", "gstar", "top_leverage"))
  cf[cols] <- lapply(cf[cols], function(x) sprintf("%.3f", 100 * x))
  cf$gstar <- round(cf$gstar)
  print(cf, row.names = FALSE)
}
