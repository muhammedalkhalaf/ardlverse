# Tests for Panel NARDL (pnardl)

test_that("pnardl basic PMG estimation works", {
  skip_on_cran()
  
  # Generate panel test data
  set.seed(123)
  N <- 10
  T <- 50
  
  panel_data <- data.frame()
  for (i in 1:N) {
    x1 <- cumsum(rnorm(T, 0, 0.5))
    y <- 2 + cumsum(pmax(diff(c(0, x1)), 0)) * 0.5 - 
         cumsum(pmin(diff(c(0, x1)), 0)) * 0.3 + rnorm(T)
    panel_data <- rbind(panel_data, data.frame(
      id = i,
      time = 1:T,
      y = y,
      x1 = x1
    ))
  }
  
  # Test PMG estimator
  result <- pnardl(y ~ x1, data = panel_data, id = "id", time = "time",
                   estimator = "pmg")
  
  expect_s3_class(result, "pnardl")
  expect_equal(result$estimator, "pmg")
  expect_true(!is.null(result$long_run_pos))
  expect_true(!is.null(result$long_run_neg))
})

test_that("pnardl MG estimator works", {
  skip_on_cran()
  
  set.seed(456)
  N <- 8
  T <- 40
  
  panel_data <- data.frame()
  for (i in 1:N) {
    x1 <- cumsum(rnorm(T))
    y <- 2 + 0.5 * x1 + rnorm(T)
    panel_data <- rbind(panel_data, data.frame(
      id = i,
      time = 1:T,
      y = y,
      x1 = x1
    ))
  }
  
  result <- pnardl(y ~ x1, data = panel_data, id = "id", time = "time",
                   estimator = "mg")
  
  expect_s3_class(result, "pnardl")
  expect_equal(result$estimator, "mg")
})

test_that("pnardl asymmetry test works", {
  skip_on_cran()
  
  set.seed(789)
  N <- 8
  T <- 50
  
  panel_data <- data.frame()
  for (i in 1:N) {
    x1 <- cumsum(rnorm(T))
    # Create asymmetric effect
    y <- 2 + cumsum(pmax(diff(c(0, x1)), 0)) * 0.8 - 
         cumsum(pmin(diff(c(0, x1)), 0)) * 0.2 + rnorm(T)
    panel_data <- rbind(panel_data, data.frame(
      id = i,
      time = 1:T,
      y = y,
      x1 = x1
    ))
  }
  
  result <- pnardl(y ~ x1, data = panel_data, id = "id", time = "time")
  
  expect_true(!is.null(result$asymmetry_test))
  expect_true("x1" %in% names(result$asymmetry_test))
})

test_that("pnardl print and summary methods work", {
  skip_on_cran()
  
  set.seed(111)
  N <- 5
  T <- 30
  
  panel_data <- data.frame()
  for (i in 1:N) {
    x1 <- cumsum(rnorm(T))
    y <- 2 + 0.5 * x1 + rnorm(T)
    panel_data <- rbind(panel_data, data.frame(
      id = i,
      time = 1:T,
      y = y,
      x1 = x1
    ))
  }
  
  result <- pnardl(y ~ x1, data = panel_data, id = "id", time = "time")
  
  expect_output(print(result))
  expect_output(summary(result))
})
