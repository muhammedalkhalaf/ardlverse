# Tests for qnardl functions

test_that("generate_oil_data creates correct structure", {
  data <- generate_oil_data(n = 100)
  
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 100)
  expect_true(all(c("week", "gasoline", "oil_price", "exchange_rate") %in% names(data)))
})

test_that("qnardl estimation works", {
  skip_on_cran()
  
  data <- generate_oil_data(n = 100, seed = 123)
  
  model <- qnardl(
    gasoline ~ oil_price,
    data = data,
    tau = c(0.25, 0.5, 0.75),
    p = 1, q = 1
  )
  
  expect_s3_class(model, "qnardl")
  expect_equal(model$tau, c(0.25, 0.5, 0.75))
  expect_true(!is.null(model$long_run_pos))
  expect_true(!is.null(model$long_run_neg))
  expect_true(!is.null(model$asymmetry_test))
})

test_that("qnardl validates decompose variables", {
  data <- generate_oil_data(n = 50)
  
  expect_error(
    qnardl(gasoline ~ oil_price, data = data, decompose = c("wrong_var")),
    "must be in the formula"
  )
})

test_that("asymmetry_test works", {
  skip_on_cran()
  
  data <- generate_oil_data(n = 100, seed = 456)
  
  model <- qnardl(
    gasoline ~ oil_price,
    data = data,
    tau = c(0.25, 0.5, 0.75),
    p = 1, q = 1
  )
  
  # Should not error
  expect_output(asymmetry_test(model, var = "oil_price"))
})
