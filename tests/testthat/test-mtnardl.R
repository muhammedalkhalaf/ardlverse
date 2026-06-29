# Tests for Multiple-Threshold NARDL (mtnardl)

test_that("mtnardl basic functionality works", {
  skip_on_cran()
  
  # Generate test data
  set.seed(123)
  n <- 200
  x1 <- cumsum(rnorm(n, 0, 0.5))
  y <- 2 + cumsum(pmax(diff(c(0, x1)), 0)) * 0.5 - 
       cumsum(pmin(diff(c(0, x1)), 0)) * 0.3 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  # Test with single threshold (standard NARDL)
  result <- mtnardl(y ~ x1, data = data, thresholds = c(0))
  
  expect_s3_class(result, "mtnardl")
  expect_equal(result$n_regimes, 2)
  expect_true(!is.null(result$long_run))
})

test_that("mtnardl multiple thresholds work", {
  skip_on_cran()
  
  set.seed(456)
  n <- 200
  x1 <- cumsum(rnorm(n, 0, 0.5))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  # Test with multiple thresholds
  result <- mtnardl(y ~ x1, data = data, thresholds = c(-0.1, 0, 0.1))
  
  expect_s3_class(result, "mtnardl")
  expect_equal(result$n_regimes, 4)
  expect_equal(length(result$thresholds), 3)
})

test_that("mtnardl asymmetry tests work", {
  skip_on_cran()
  
  set.seed(789)
  n <- 200
  x1 <- cumsum(rnorm(n))
  y <- 2 + cumsum(pmax(diff(c(0, x1)), 0)) * 0.8 - 
       cumsum(pmin(diff(c(0, x1)), 0)) * 0.2 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- mtnardl(y ~ x1, data = data)
  
  expect_true(!is.null(result$asymmetry_tests))
  expect_true("x1" %in% names(result$asymmetry_tests))
})

test_that("mtnardl print, summary, and plot methods work", {
  skip_on_cran()
  
  set.seed(111)
  n <- 150
  x1 <- cumsum(rnorm(n))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- mtnardl(y ~ x1, data = data)
  
  expect_output(print(result))
  expect_output(summary(result))
  expect_silent(plot(result, type = "asymmetry"))
})
