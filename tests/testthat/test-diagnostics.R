# Tests for diagnostic functions

test_that("ardl_diagnostics works with fourier_ardl", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 123)
  
  model <- fourier_ardl(
    gdp ~ investment,
    data = data,
    p = 1, q = 1,
    k = 1
  )
  
  diag <- ardl_diagnostics(model)
  
  expect_s3_class(diag, "ardl_diagnostics")
  expect_true(!is.null(diag$serial_corr))
  expect_true(!is.null(diag$hetero_bp))
  expect_true(!is.null(diag$arch))
  expect_true(!is.null(diag$normality))
  expect_true(!is.null(diag$cusum))
  expect_true(!is.null(diag$cusumsq))
})

test_that("ardl_diagnostics works with boot_ardl", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 456)
  
  model <- boot_ardl(
    gdp ~ investment,
    data = data,
    p = 1, q = 1,
    nboot = 50
  )
  
  diag <- ardl_diagnostics(model)
  
  expect_s3_class(diag, "ardl_diagnostics")
  expect_true(!is.null(diag$residuals))
})

test_that("CUSUM test returns correct structure", {
  resid <- rnorm(100)
  
  cusum <- ardlverse:::.cusum_test(resid)
  
  expect_true(!is.null(cusum$cusum))
  expect_true(!is.null(cusum$upper_bound))
  expect_true(!is.null(cusum$lower_bound))
  expect_true(is.logical(cusum$crosses_bounds))
})

test_that("Jarque-Bera test works", {
  # Normal residuals
  resid_normal <- rnorm(100)
  jb_normal <- ardlverse:::.jarque_bera_test(resid_normal)
  
  expect_true(jb_normal$statistic >= 0)
  expect_true(jb_normal$p.value >= 0 && jb_normal$p.value <= 1)
  
  # Non-normal residuals (uniform)
  resid_uniform <- runif(100)
  jb_uniform <- ardlverse:::.jarque_bera_test(resid_uniform)
  
  # Uniform should have lower p-value (more likely to reject normality)
  # This is a probabilistic test, so we don't make strict assertions
  expect_true(!is.null(jb_uniform$skewness))
  expect_true(!is.null(jb_uniform$kurtosis))
})
