# Tests for panel_ardl functions

test_that("generate_panel_data creates correct structure", {
  data <- generate_panel_data(n_groups = 5, n_time = 20)
  
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 5 * 20)
  expect_true(all(c("country", "year", "gdp", "inflation", "investment", "trade") %in% names(data)))
})

test_that("panel_ardl PMG estimation works", {
  skip_on_cran()
  
  data <- generate_panel_data(n_groups = 5, n_time = 30, seed = 123)
  
  model <- panel_ardl(
    gdp ~ inflation + investment,
    data = data,
    id = "country",
    time = "year",
    p = 1, q = 1,
    estimator = "pmg"
  )
  
  expect_s3_class(model, "panel_ardl")
  expect_equal(model$estimator, "pmg")
  expect_true(!is.null(model$long_run))
  expect_true(model$ec_coef < 0)  # EC coefficient should be negative
})

test_that("panel_ardl MG estimation works", {
  skip_on_cran()
  
  data <- generate_panel_data(n_groups = 5, n_time = 30, seed = 123)
  
  model <- panel_ardl(
    gdp ~ inflation + investment,
    data = data,
    id = "country",
    time = "year",
    p = 1, q = 1,
    estimator = "mg"
  )
  
  expect_s3_class(model, "panel_ardl")
  expect_equal(model$estimator, "mg")
})

test_that("panel_ardl DFE estimation works", {
  skip_on_cran()
  
  data <- generate_panel_data(n_groups = 5, n_time = 30, seed = 123)
  
  model <- panel_ardl(
    gdp ~ inflation + investment,
    data = data,
    id = "country",
    time = "year",
    p = 1, q = 1,
    estimator = "dfe"
  )
  
  expect_s3_class(model, "panel_ardl")
  expect_equal(model$estimator, "dfe")
})

test_that("panel_ardl validates inputs", {
  data <- generate_panel_data(n_groups = 3, n_time = 10)
  
  # Missing id variable
  expect_error(
    panel_ardl(gdp ~ inflation, data = data, id = "wrong_id", time = "year"),
    "not found"
  )
  
  # Missing time variable
  expect_error(
    panel_ardl(gdp ~ inflation, data = data, id = "country", time = "wrong_time"),
    "not found"
  )
})

test_that("hausman_test works", {
  skip_on_cran()
  
  data <- generate_panel_data(n_groups = 5, n_time = 30, seed = 456)
  
  pmg <- panel_ardl(
    gdp ~ inflation + investment,
    data = data,
    id = "country",
    time = "year",
    estimator = "pmg"
  )
  
  mg <- panel_ardl(
    gdp ~ inflation + investment,
    data = data,
    id = "country",
    time = "year",
    estimator = "mg"
  )
  
  result <- hausman_test(pmg, mg)
  
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_true(result$df > 0)
})
