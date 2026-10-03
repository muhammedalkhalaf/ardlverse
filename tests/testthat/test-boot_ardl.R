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

test_that("boot_ardl statistics equal hand restricted-RSS tests", {
  set.seed(21)
  n <- 80
  d <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
  b <- boot_ardl(y ~ x, d, p = 2, q = 1, nboot = 19, seed = 1)
  t <- 3:n
  dy <- d$y[t] - d$y[t - 1]; ly <- d$y[t - 1]; lx <- d$x[t - 1]
  dyl <- d$y[t - 1] - d$y[t - 2]; dx0 <- d$x[t] - d$x[t - 1]; dx1 <- d$x[t - 1] - d$x[t - 2]
  u <- lm(dy ~ ly + lx + dyl + dx0 + dx1)
  expect_equal(b$F_stat, anova(lm(dy ~ dyl + dx0 + dx1), u)$F[2], tolerance = 1e-10)
  expect_equal(b$Find_stat, anova(lm(dy ~ ly + dyl + dx0 + dx1), u)$F[2], tolerance = 1e-10)
  expect_equal(unname(b$t_stat), summary(u)$coefficients["ly", 3], tolerance = 1e-10)
  expect_equal(b$F_overall, b$F_stat)
  expect_lt(b$engine$design_check, 1e-10)
  expect_lt(max(b$engine$dgpcheck), 1e-8)
  expect_true(b$decision %in% c("COINTEGRATION", "NO_COINTEGRATION",
                                "DEGENERATE_1", "DEGENERATE_2"))
  bj <- boot_ardl(y ~ x, d, p = 2, q = 1, nboot = 19, seed = 1, nulls = "joint")
  expect_equal(bj$engine$settings$xmodel, "var")
  expect_equal(bj$engine$settings$recentre, "once")
})

test_that("boot_ardl rejects transformed terms and missing values", {
  d <- data.frame(y = cumsum(rnorm(40)) + 50, x = cumsum(rnorm(40)) + 50)
  expect_error(boot_ardl(log(y) ~ x, d, nboot = 9), "transformed")
  expect_error(boot_ardl(y ~ log(x), d, nboot = 9), "transformed")
  d$x[5] <- NA
  expect_error(boot_ardl(y ~ x, d, nboot = 9), "missing")
})
