# Monte Carlo size of the recursive bootstrap (skipped on CRAN; about five
# minutes). y and x are independent random walks, so the nulls of Fov, t and
# Find all hold. Each 5% rejection rate must be below 0.05 + 3 Monte Carlo
# standard errors (0.134 with 60 samples). Limits: with 60 samples and 99
# replications the test detects gross failures such as the fixed-design
# bootstraps of version 2.0.3 (29% to 98% rejections in this design), not
# over-rejection of a few percentage points; the reference sizes are those of
# the larger Monte Carlo reported in NEWS.md (200 samples, B = 199). With
# Fourier terms and partial sums the Find test of the separate-null bootstrap
# over-rejects (about 11% to 16% in the larger Monte Carlo), so only Fov and t
# are checked for aardl type "fbnardl".

size_mc <- function(fun, MC = 60, n = 100) {
  rej <- matrix(NA, MC, 3)
  for (m in seq_len(MC)) {
    set.seed(5000 + m)
    d <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
    rej[m, ] <- suppressWarnings(fun(d)) < 0.05
  }
  colMeans(rej, na.rm = TRUE)
}
tol <- 0.05 + 3 * sqrt(0.05 * 0.95 / 60)

test_that("boot_ardl size", {
  skip_on_cran()
  r <- size_mc(function(d) {
    b <- boot_ardl(y ~ x, d, nboot = 99, seed = 1)
    c(b$p_value_F, b$p_value_t, b$p_value_Find)
  })
  expect_true(all(r <= tol))
})

test_that("aardl bootstrap size (bnardl and fbnardl)", {
  skip_on_cran()
  r <- size_mc(function(d) aardl(y ~ x, d, type = "bnardl", nboot = 99, seed = 1)$conclusion$p_values)
  expect_true(all(r <= tol))
  r <- size_mc(function(d) aardl(y ~ x, d, type = "fbnardl", nboot = 99, seed = 1)$conclusion$p_values)
  expect_true(all(r[1:2] <= tol))
})

test_that("mtnardl bootstrap size", {
  skip_on_cran()
  r <- size_mc(function(d) mtnardl(y ~ x, d, bootstrap = TRUE, nboot = 99, seed = 1)$bounds_test$p_values)
  expect_true(all(r <= tol))
})

test_that("fbnardl bootstrap size", {
  skip_on_cran()
  r <- size_mc(function(d) fbnardl(y ~ x, d, decompose = "x", type = "fbnardl", maxlag = 1,
                                   reps = 99, bands = 0, seed = 1)$bounds_test$pvalue_boot)
  expect_true(all(r <= tol))
})
