# Tests for boot_ardl functions

test_that("generate_ts_data creates correct structure", {
  data <- generate_ts_data(n = 50)
  
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 50)
  expect_true(all(c("quarter", "gdp", "inflation", "investment", "trade") %in% names(data)))
})

test_that("boot_ardl estimation works", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 123)
  
  # Use fewer bootstrap reps for testing
  model <- boot_ardl(
    gdp ~ inflation + investment,
    data = data,
    p = 1, q = 1,
    case = 3,
    nboot = 100  # Reduced for testing
  )
  
  expect_s3_class(model, "boot_ardl")
  expect_true(!is.null(model$F_stat))
  expect_true(!is.null(model$t_stat))
  expect_true(length(model$boot_F) == 100)
  expect_true(!is.null(model$conclusion))
})

test_that("boot_ardl validates case parameter", {
  data <- generate_ts_data(n = 50)
  
  expect_error(
    boot_ardl(gdp ~ inflation, data = data, case = 6),
    "'case' must be"
  )
})

test_that("pss_critical_values returns correct structure", {
  cv <- pss_critical_values(k = 2, case = 3, level = "5%")
  
  expect_true(!is.null(cv$F_bounds$I0))
  expect_true(!is.null(cv$F_bounds$I1))
  expect_true(!is.null(cv$t_bounds$I0))
  expect_true(!is.null(cv$t_bounds$I1))
  expect_equal(cv$k, 2)
  expect_equal(cv$case, 3)
})
