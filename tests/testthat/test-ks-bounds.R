## Reference values from Stata 'ardlbounds' (package 'ardl' 1.0.6,
## Kripfganz and Schneider 2020)

test_that("Kripfganz-Schneider bounds reproduce ardlbounds", {
  r <- ardlverse:::.ks_bounds("F", 3, 2, 80, 4, value = 3.9)
  expect_equal(unname(r$cv["I0", ]), c(3.2043, 3.8812, 5.4116), tolerance = 1e-4)
  expect_equal(unname(r$cv["I1", ]), c(4.2082, 4.9866, 6.7118), tolerance = 1e-4)
  expect_equal(unname(r$pvalue), c(0.04903, 0.13032), tolerance = 1e-4)
  r <- ardlverse:::.ks_bounds("t", 1, 1, NULL, value = -2.2)
  expect_equal(unname(r$pvalue), c(0.02666, 0.11369), tolerance = 1e-4)
})

test_that("pss_critical_values returns the 5% Kripfganz-Schneider bounds", {
  cv <- pss_critical_values(2, 3, "5%", n = 80, sr = 4)
  expect_equal(cv$F_bounds$I0, 3.8812, tolerance = 1e-4)
  expect_equal(cv$F_bounds$I1, 4.9866, tolerance = 1e-4)
})

test_that("AARDL critical values come from the response surfaces", {
  cv <- ardlverse:::.aardl_critical_values(2, 3, 80, 4)
  expect_equal(unname(cv$F$I0), c(3.2043, 3.8812, 5.4116), tolerance = 1e-4)
  expect_named(cv$F$I1, c("90%", "95%", "99%"))
})

test_that("CUSUM bound for OLS residuals is the constant 1.358 sqrt(n)", {
  e <- stats::rnorm(50)
  r <- ardlverse:::.cusum_test(e - mean(e))
  expect_equal(r$upper_bound, rep(1.358 * sqrt(50), 50))
})

test_that("Fourier bounds test reports no invented critical values", {
  set.seed(1)
  n <- 100
  x <- cumsum(stats::rnorm(n))
  y <- 0.5 * x + stats::rnorm(n)
  d <- data.frame(y = y, x = x)
  f <- fourier_ardl(y ~ x, data = d, p = 1, q = 1, k = 1, case = 3)
  r <- suppressMessages(utils::capture.output(res <- fourier_bounds_test(f)))
  expect_true(is.na(res$cv_F))
  expect_true(all(is.finite(res$cv_pss_reference)))
})

test_that("Hausman test returns NA rather than a fallback variance", {
  h <- ardlverse:::.panel_hausman_test(
    list(long_run = c(1, 2), vcov_lr = diag(2)),
    list(long_run = c(1.1, 2.1), vcov_lr = diag(2) * 0.5))
  expect_true(is.na(h$statistic))
  h <- ardlverse:::.panel_hausman_test(
    list(long_run = c(1, 2), vcov_lr = diag(2) * 0.5),
    list(long_run = c(1.5, 2), vcov_lr = diag(2)))
  expect_equal(h$statistic, 0.5)
  expect_equal(h$df, 2)
})
