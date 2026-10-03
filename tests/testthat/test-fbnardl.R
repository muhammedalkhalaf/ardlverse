# fbnardl(): values pinned to the auditor's examples (Stata module fbnardl
# conventions: common sample, delta-method long-run test, sum-of-lags short-run
# test, full-recursion multipliers)

fb_sim <- function(n, seed) {
  set.seed(seed)
  x <- cumsum(rnorm(n)); z <- cumsum(rnorm(n))
  dx <- c(0, diff(x)); xp <- cumsum(pmax(dx, 0)); xn <- cumsum(pmin(dx, 0))
  y <- numeric(n); y[1] <- 0
  for (t in 3:n) {
    ec <- -0.3 * (y[t - 1] - 1.0 * xp[t - 1] - 0.4 * xn[t - 1] - 0.5 * z[t - 1])
    y[t] <- y[t - 1] + ec + 0.3 * (y[t - 1] - y[t - 2]) + 0.5 * dx[t] * (dx[t] > 0) +
      0.2 * (x[t - 1] - x[t - 2]) + 0.5 * sin(2 * pi * 1.3 * t / n) * 0.1 + rnorm(1, sd = 0.5)
  }
  data.frame(y = y, x = x, z = z)
}

test_that("partial sums: x+ + x- = x - x1, start at zero, monotone", {
  d <- fb_sim(60, 2)
  W <- ardlverse:::.fb_levels(as.matrix(d[, c("x", "z")]),
                              list(dec = "x", ctrl = "z"))
  expect_equal(W[, "x_pos"] + W[, "x_neg"], d$x - d$x[1], tolerance = 1e-12)
  expect_equal(unname(W[1, 1:2]), c(0, 0))
  expect_true(all(diff(W[, "x_pos"]) >= 0) && all(diff(W[, "x_neg"]) <= 0))
  # identical to .decompose_asymmetric(threshold = 0)
  expect_equal(unname(W[, 1:2]), ardlverse:::.decompose_asymmetric(d$x, 0), tolerance = 1e-12)
})

test_that("fbnardl reproduces the auditor's t3 example", {
  d <- fb_sim(150, 7)
  m <- fbnardl(y ~ x + z, d, decompose = "x", maxlag = 2,
               lags = list(p = 2, q = 2, r = 1), kstar = 1.3, bands = 0)
  # coefficients equal a hand lm.fit on the common sample t = 4..150
  t <- 4:150; n <- 150
  xp <- cumsum(pmax(c(0, diff(d$x)), 0)); xn <- cumsum(pmin(c(0, diff(d$x)), 0))
  D <- function(v, l) v[t - l] - v[t - l - 1]
  Z <- cbind(1, D(d$y, 1), D(d$y, 2), D(xp, 0), D(xp, 1), D(xp, 2), D(xn, 0), D(xn, 1),
             D(xn, 2), D(d$z, 0), D(d$z, 1), d$y[t - 1], xp[t - 1], xn[t - 1], d$z[t - 1],
             sin(2 * pi * 1.3 * t / n), cos(2 * pi * 1.3 * t / n))
  f <- lm.fit(Z, d$y[t] - d$y[t - 1])
  expect_equal(unname(m$coefficients), unname(f$coefficients), tolerance = 1e-10)
  expect_equal(m$asymmetry_tests$x$lr_chi2, 631.897, tolerance = 1e-6)
  expect_equal(m$asymmetry_tests$x$sr_f, 2.941682, tolerance = 1e-6)
  expect_equal(round(m$dynamic_multipliers$x$H_pos[1:5], 5),
               c(0.52224, 0.98575, 1.31794, 1.32148, 1.20926))
  expect_equal(m$multipliers$x$lr_pos, 1.050439, tolerance = 1e-6)
  expect_equal(m$multipliers$x$lr_pos_se, 0.04888556, tolerance = 1e-6)
  bt <- m$bounds_test
  expect_equal(c(bt$Fov, bt$t, bt$Find), c(14.8975, -7.359752, 19.82841), tolerance = 1e-6)
  expect_equal(bt$k, 3)
  expect_equal(bt$decision, "NOT_AVAILABLE")
  # long-run SE equals the numerical delta method
  b <- m$coefficients; V <- m$vcov
  g <- sapply(seq_along(b), function(i) {
    e <- 1e-6; b1 <- b; b2 <- b; b1[i] <- b1[i] + e; b2[i] <- b2[i] - e
    ((-b1["L_x_pos"] / b1["L_y"]) - (-b2["L_x_pos"] / b2["L_y"])) / (2 * e) })
  expect_equal(m$multipliers$x$lr_pos_se, sqrt(drop(t(g) %*% V %*% g)), tolerance = 1e-6)
  # multipliers converge to the long-run coefficient
  expect_equal(m$dynamic_multipliers$x$H_pos[21], m$multipliers$x$lr_pos, tolerance = 1e-4)
})

test_that("common-sample lag selection: seed 1 selects (1, 0, 3)", {
  d <- fb_sim(80, 1)
  m <- fbnardl(y ~ x + z, d, decompose = "x", maxlag = 3, maxk = 1,
               kgrid = "fractional", bands = 0)
  expect_equal(m$kstar, 0.8)
  expect_equal(c(m$lags$p, m$lags$q, m$lags$r), c(1, 0, 3), ignore_attr = TRUE)
  expect_equal(m$nobs, 80 - 3 - 1)
  mi <- fbnardl(y ~ x + z, d, decompose = "x", maxlag = 3, maxk = 3, bands = 0)
  expect_true(mi$kstar %in% 1:3)
  expect_equal(mi$ssr_by_k$k, 1:3)
})

test_that("bounds without Fourier terms equal .ks_bounds with n and sr", {
  d <- fb_sim(100, 3)
  m <- fbnardl(y ~ x + z, d, decompose = "x", maxlag = 2, fourier = FALSE,
               lags = list(p = 1, q = 1, r = 0), bands = 0)
  bt <- m$bounds_test
  expect_equal(bt$n, 100 - 3)
  expect_equal(bt$sr, 1 + 2 * 2 + 1)
  ks <- ardlverse:::.ks_bounds("F", 3, 3, 97, 6, siglevels = c(10, 5, 2.5, 1))$cv
  expect_equal(bt$bounds_F, ks)
  expect_true(bt$bounds_valid)
  expect_true(bt$decision %in% c("COINTEGRATION", "NO_COINTEGRATION", "INCONCLUSIVE"))
})

test_that("fbnardl bootstrap: engine design, NA failures, local RNG", {
  d <- fb_sim(80, 5)
  set.seed(9); s0 <- .Random.seed
  m <- fbnardl(y ~ x, d, decompose = "x", type = "fbnardl", maxlag = 1, maxk = 2,
               reps = 9, bands = 10, seed = 1)
  expect_identical(.Random.seed, s0)
  eng <- m$bounds_test$bootstrap
  expect_lt(eng$design_check, 1e-10)
  expect_lt(max(eng$dgpcheck), 1e-8)
  expect_true(eng$settings$reselect)
  expect_equal(eng$settings$nulls, "separate")
  m2 <- fbnardl(y ~ x, d, decompose = "x", type = "fbnardl", maxlag = 1, kstar = 1,
                lags = list(p = 1, q = 0), reps = 9, bands = 0, seed = 1,
                bootstrap = "mcnown")
  expect_false(m2$bounds_test$bootstrap$settings$reselect)
  expect_equal(m2$bounds_test$bootstrap$settings$nulls, "joint")
})

test_that("fbnardl input checks", {
  d <- fb_sim(60, 6)
  expect_error(fbnardl(log(y + 100) ~ x, d, decompose = "x"), "transformed")
  d2 <- d; d2$x[10] <- NA
  expect_error(fbnardl(y ~ x, d2, decompose = "x"), "missing")
  expect_error(fbnardl(y ~ x, d, decompose = "w"), "decompose")
})
