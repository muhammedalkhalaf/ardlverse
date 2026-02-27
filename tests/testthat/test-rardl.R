# Tests for Rolling/Recursive ARDL (rardl)

test_that("rardl rolling method works", {
  skip_on_cran()
  
  # Generate test data
  set.seed(123)
  n <- 200
  x1 <- cumsum(rnorm(n, 0, 0.5))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  # Test rolling method
  result <- rardl(y ~ x1, data = data, method = "rolling", window = 60)
  
  expect_s3_class(result, "rardl")
  expect_equal(result$method, "rolling")
  expect_true(length(result$F_stats) > 0)
  expect_true(length(result$ec_coefs) > 0)
})

test_that("rardl recursive method works", {
  skip_on_cran()
  
  set.seed(456)
  n <- 200
  x1 <- cumsum(rnorm(n))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  # Test recursive method
  result <- rardl(y ~ x1, data = data, method = "recursive", min_obs = 50)
  
  expect_s3_class(result, "rardl")
  expect_equal(result$method, "recursive")
})

test_that("rardl detects structural breaks", {
  skip_on_cran()
  
  set.seed(789)
  n <- 200
  x1 <- cumsum(rnorm(n))
  # Create structural break at midpoint
  y <- c(2 + 0.5 * x1[1:100], 5 + 0.8 * x1[101:200]) + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- rardl(y ~ x1, data = data, method = "rolling", window = 50)
  
  # Should detect some change in coefficients
  expect_true(sd(result$ec_coefs, na.rm = TRUE) > 0)
})

test_that("rardl print, summary, and plot methods work", {
  skip_on_cran()
  
  set.seed(111)
  n <- 150
  x1 <- cumsum(rnorm(n))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- rardl(y ~ x1, data = data, window = 50)
  
  expect_output(print(result))
  expect_output(summary(result))
  expect_silent(plot(result, type = "F"))
})

test_that("rardl validates window size", {
  set.seed(222)
  n <- 100
  data <- data.frame(y = rnorm(n), x1 = rnorm(n))
  
  # Window too small
  expect_error(rardl(y ~ x1, data = data, window = 5))
})
