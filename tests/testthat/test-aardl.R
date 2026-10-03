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

test_that("aardl statistics equal hand Wald tests and Fourier types give no decision", {
  set.seed(31)
  n <- 90
  d <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
  a <- aardl(y ~ x, d, p = 1, q = 1, case = 3)
  t <- 3:n
  dy <- d$y[t] - d$y[t - 1]; ly <- d$y[t - 1]; lx <- d$x[t - 1]
  dyl <- d$y[t - 1] - d$y[t - 2]; dx0 <- d$x[t] - d$x[t - 1]
  u <- lm(dy ~ ly + lx + dyl + dx0)
  expect_equal(a$F_pss, anova(lm(dy ~ dyl + dx0), u)$F[2], tolerance = 1e-10)
  expect_equal(a$F_ind, anova(lm(dy ~ ly + dyl + dx0), u)$F[2], tolerance = 1e-10)
  expect_equal(a$t_dep, summary(u)$coefficients["ly", 3], tolerance = 1e-10)
  for (ty in c("fourier", "fnardl")) {
    f <- aardl(y ~ x, d, type = ty)
    expect_equal(f$conclusion$decision, "NOT_AVAILABLE")
    expect_match(f$conclusion$message, "not valid with Fourier terms")
  }
})

test_that("aardl bootstrap types use the recursive engine with the model's design", {
  set.seed(32)
  n <- 80
  d <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
  for (ty in c("bootstrap", "bnardl", "fbootstrap", "fbnardl")) {
    a <- aardl(y ~ x, d, type = ty, nboot = 9, seed = 1)
    eng <- a$boot_results$engine
    expect_lt(eng$design_check, 1e-10)
    expect_lt(max(eng$dgpcheck), 1e-8)
    # independent lm()/anova() computation of the three statistics
    t <- 3:n
    if (ty %in% c("bnardl", "fbnardl")) {
      dd <- c(0, diff(d$x))
      W <- cbind(cumsum(dd * (dd > 0)), cumsum(dd * (dd < 0)))
    } else {
      W <- cbind(d$x)
    }
    dy <- d$y[t] - d$y[t - 1]; ly <- d$y[t - 1]; lw <- W[t - 1, , drop = FALSE]
    dyl <- d$y[t - 1] - d$y[t - 2]; dw0 <- W[t, , drop = FALSE] - W[t - 1, , drop = FALSE]
    if (ty %in% c("fbootstrap", "fbnardl")) {
      fo <- cbind(sin(2 * pi * t / n), cos(2 * pi * t / n))
      u <- lm(dy ~ ly + lw + dyl + dw0 + fo)
      r_ov <- lm(dy ~ dyl + dw0 + fo); r_in <- lm(dy ~ ly + dyl + dw0 + fo)
    } else {
      u <- lm(dy ~ ly + lw + dyl + dw0)
      r_ov <- lm(dy ~ dyl + dw0); r_in <- lm(dy ~ ly + dyl + dw0)
    }
    expect_equal(a$F_pss, anova(r_ov, u)$F[2], tolerance = 1e-10)
    expect_equal(a$F_ind, anova(r_in, u)$F[2], tolerance = 1e-10)
    expect_equal(a$t_dep, summary(u)$coefficients["ly", 3], tolerance = 1e-10)
    expect_equal(unname(eng$statistic), c(a$F_pss, a$t_dep, a$F_ind), tolerance = 1e-10)
    expect_true(all(a$conclusion$p_values >= 0 & a$conclusion$p_values <= 1))
    expect_equal(a$conclusion$method, "bootstrap")
  }
})
