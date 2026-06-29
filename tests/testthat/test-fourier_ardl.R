# Tests for fourier_ardl functions

test_that("fourier_ardl estimation works with fixed k", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 123)
  
  model <- fourier_ardl(
    gdp ~ investment + trade,
    data = data,
    p = 1, q = 1,
    k = 1,
    selection = "fixed"
  )
  
  expect_s3_class(model, "fourier_ardl")
  expect_equal(model$k, 1)
  expect_true(!is.null(model$fourier_coefs))
  expect_true(!is.null(model$bounds_test$F_stat))
})

test_that("fourier_ardl works with AIC selection", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 456)
  
  model <- fourier_ardl(
    gdp ~ investment,
    data = data,
    p = 1, q = 1,
    selection = "aic"
  )
  
  expect_s3_class(model, "fourier_ardl")
  expect_true(model$k %in% 1:3)
  expect_true(!is.null(model$ic_comparison))
})

test_that("fourier_ardl validates k parameter", {
  data <- generate_ts_data(n = 50)
  
  expect_error(
    fourier_ardl(gdp ~ investment, data = data, k = 5),
    "'k' must be between"
  )
})

test_that("fourier_bounds_test works", {
  skip_on_cran()
  
  data <- generate_ts_data(n = 80, seed = 789)
  
  model <- fourier_ardl(
    gdp ~ investment,
    data = data,
    p = 1, q = 1,
    k = 1
  )
  
  # Should not error and produce output
  expect_output(fourier_bounds_test(model))
})
