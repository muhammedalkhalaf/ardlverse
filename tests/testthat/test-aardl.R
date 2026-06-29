# Tests for Augmented ARDL (aardl)

test_that("aardl basic functionality works", {
  skip_on_cran()
  
  # Generate test data
  set.seed(123)
  n <- 150
  x1 <- cumsum(rnorm(n, 0, 0.5))
  x2 <- cumsum(rnorm(n, 0, 0.3))
  y <- 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n, 0, 1)
  data <- data.frame(y = y, x1 = x1, x2 = x2)
  
  # Test basic estimation
  result <- aardl(y ~ x1 + x2, data = data, p = 1, q = 1, case = 3)
  
  expect_s3_class(result, "aardl")
  expect_true(!is.null(result$F_pss))
  expect_true(!is.null(result$t_dep))
  expect_true(!is.null(result$conclusion))
})

test_that("aardl NARDL type works", {
  skip_on_cran()
  
  set.seed(456)
  n <- 150
  x1 <- cumsum(rnorm(n))
  y <- 2 + cumsum(pmax(diff(c(0, x1)), 0)) * 0.5 - 
       cumsum(pmin(diff(c(0, x1)), 0)) * 0.3 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- aardl(y ~ x1, data = data, type = "nardl")
  
  expect_s3_class(result, "aardl")
  expect_equal(result$type, "nardl")
})

test_that("aardl Fourier type works", {
  skip_on_cran()
  
  set.seed(789)
  n <- 200
  t <- 1:n
  x1 <- cumsum(rnorm(n))
  y <- 2 + 0.5 * x1 + sin(2 * pi * t / n) + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- aardl(y ~ x1, data = data, type = "fourier", fourier_k = 2)
  
  expect_s3_class(result, "aardl")
  expect_equal(result$type, "fourier")
  expect_equal(result$fourier_k, 2)
})

test_that("aardl print and summary methods work", {
  skip_on_cran()
  
  set.seed(111)
  n <- 100
  x1 <- cumsum(rnorm(n))
  y <- 2 + 0.5 * x1 + rnorm(n)
  data <- data.frame(y = y, x1 = x1)
  
  result <- aardl(y ~ x1, data = data)
  
  expect_output(print(result))
  expect_output(summary(result))
})

test_that("aardl validates inputs correctly", {
  set.seed(222)
  n <- 100
  data <- data.frame(y = rnorm(n), x1 = rnorm(n))
  
  # Invalid case
  expect_error(aardl(y ~ x1, data = data, case = 6))
  
  # Invalid fourier_k
  expect_error(aardl(y ~ x1, data = data, type = "fourier", fourier_k = 5))
})
