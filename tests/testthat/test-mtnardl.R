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

# Hand reference (independent of the package): partial sums by loop, ECM by
# lm.fit, multiplier by simulating the estimated ECM
mt_ps <- function(x, th) {
  TT <- length(x); R <- length(th) + 1; lo <- c(-Inf, th); hi <- c(th, Inf)
  out <- matrix(0, TT, R)
  for (t in 2:TT) {
    dd <- x[t] - x[t - 1]; out[t, ] <- out[t - 1, ]
    for (r in 1:R) if (dd > lo[r] && dd <= hi[r]) out[t, r] <- out[t, r] + dd
  }
  out
}
mt_sim <- function(TT, seed, th = c(-2, 0, 2)) {
  set.seed(seed)
  beta <- c(0.15, 0.3, 0.5, 0.8); gam <- c(0.1, 0.2, 0.3, 0.4)
  x <- 50 + cumsum(c(0, rnorm(TT - 1, 0.05, 3)))
  ps <- mt_ps(x, th)
  y <- numeric(TT); y[1:2] <- 30
  e <- rnorm(TT, 0, 0.8)
  for (t in 3:TT) {
    ec <- -0.3 * (y[t - 1] - 30 - sum(beta * ps[t - 1, ]))
    y[t] <- y[t - 1] + ec + sum(gam * (ps[t, ] - ps[t - 1, ])) + 0.2 * (y[t - 1] - y[t - 2]) + e[t]
  }
  data.frame(y = y, x = x)
}
mt_hand <- function(y, X, p, q) {
  TT <- length(y); k <- ncol(X); rows <- (max(p + 1, q) + 1):TT
  Z <- cbind(y[rows - 1], X[rows - 1, ])
  for (i in 1:p) Z <- cbind(Z, y[rows - i] - y[rows - i - 1])
  for (j in 1:k) for (i in 0:(q - 1)) Z <- cbind(Z, X[rows - i, j] - X[rows - i - 1, j])
  Z <- cbind(Z, 1)
  f <- lm.fit(Z, y[rows] - y[rows - 1])
  V <- sum(f$residuals^2) / (length(rows) - ncol(Z)) * solve(crossprod(Z))
  list(b = f$coefficients, V = V, k = k, p = p, q = q)
}
mt_mult_sim <- function(h, j, H = 30) {
  b <- h$b; k <- h$k; p <- h$p; q <- h$q; pad <- 20; N <- pad + H + 1
  y <- rep(0, N); X <- matrix(0, N, k); X[(pad + 1):N, j] <- 1
  for (t in (pad + 1):N) {
    dy <- b[1] * y[t - 1] + sum(b[1 + 1:k] * X[t - 1, ])
    for (i in 1:p) dy <- dy + b[1 + k + i] * (y[t - i] - y[t - i - 1])
    for (m in 1:k) for (i in 0:(q - 1))
      dy <- dy + b[1 + k + p + (m - 1) * q + i + 1] * (X[t - i, m] - X[t - i - 1, m])
    y[t] <- y[t - 1] + dy
  }
  y[(pad + 1):N]
}

test_that("mtnardl matches the hand ECM: F, t, long run, delta-method SE, joint Wald", {
  d <- mt_sim(200, 11)
  th <- c(-2, 0, 2)
  m <- mtnardl(y ~ x, d, thresholds = th, p = 2, q = 3)
  h <- mt_hand(d$y, mt_ps(d$x, th), 2, 3)
  b <- h$b; V <- h$V
  lr <- -b[2:5] / b[1]
  expect_equal(unname(m$long_run), unname(lr), tolerance = 1e-10)
  R <- diag(length(b))[1:5, ]
  Fh <- drop(t(R %*% b) %*% solve(R %*% V %*% t(R), R %*% b)) / 5
  expect_equal(m$F_stat, Fh, tolerance = 1e-9)
  expect_equal(m$t_stat, unname(b[1] / sqrt(V[1, 1])), tolerance = 1e-9)
  # delta-method SE equals the numerical-gradient delta method
  ng <- function(f, b, e = 1e-6) sapply(seq_along(b), function(i) {
    b1 <- b; b2 <- b; b1[i] <- b1[i] + e; b2[i] <- b2[i] - e; (f(b1) - f(b2)) / (2 * e) })
  for (j in 1:4) {
    g <- ng(function(bb) -bb[1 + j] / bb[1], b)
    expect_equal(m$long_run_table$std_error[j], sqrt(drop(t(g) %*% V %*% g)), tolerance = 1e-6)
  }
  Rj <- matrix(0, 3, length(b)); for (r in 2:4) { Rj[r - 1, 1 + r] <- 1; Rj[r - 1, 2] <- -1 }
  Wj <- drop(t(Rj %*% b) %*% solve(Rj %*% V %*% t(Rj), Rj %*% b))
  expect_equal(m$asymmetry_tests$x$joint$wald, Wj, tolerance = 1e-9)
  expect_equal(m$asymmetry_tests$x$joint$df, 3)
})

test_that("mtnardl dynamic multipliers equal the simulated level response of the ECM", {
  d <- mt_sim(200, 11)
  th <- c(-2, 0, 2)
  m <- mtnardl(y ~ x, d, thresholds = th, p = 2, q = 3, horizon = 30)
  h <- mt_hand(d$y, mt_ps(d$x, th), 2, 3)
  for (j in 1:4) expect_lt(max(abs(m$multipliers$cumulative[, j] - mt_mult_sim(h, j))), 1e-10)
  expect_equal(m$multipliers$horizons, 0:30)
  # auditor's reference values (sim_mt(250, 1), p = 1, q = 1): h = 0, 4, 29
  d2 <- mt_sim(250, 1)
  m2 <- mtnardl(y ~ x, d2, thresholds = th, p = 1, q = 1)
  ref <- rbind(c(0.090375, 0.303465, 0.339966, 0.416636),
               c(0.130116, 0.316002, 0.551650, 0.729079),
               c(0.132587, 0.299847, 0.573519, 0.766601))
  expect_equal(unname(m2$multipliers$cumulative[c(1, 5, 30), ]), ref, tolerance = 1e-5)
})

test_that("mtnardl auto_select with one threshold searches the candidates", {
  d <- mt_sim(200, 11)
  m <- mtnardl(y ~ x, d, auto_select = TRUE, n_thresholds = 1)
  cand <- sort(unique(c(0, quantile(diff(d$x), seq(0.1, 0.9, by = 0.1)))))
  aics <- sapply(cand, function(th) mtnardl(y ~ x, d, thresholds = th)$fit$AIC)
  expect_equal(m$thresholds, unname(cand[which.min(aics)]))
})

test_that("mtnardl bounds k, transformed terms, decompose subset and quantile partition", {
  d <- mt_sim(150, 3)
  m1 <- mtnardl(y ~ x, d, thresholds = c(-2, 0, 2))
  expect_equal(m1$critical_values$bounds_k_value, 4)
  m2 <- mtnardl(y ~ x, d, thresholds = c(-2, 0, 2), bounds_k = "original")
  expect_equal(m2$critical_values$bounds_k_value, 1)
  expect_error(mtnardl(log(y) ~ x, d), "transformed")
  set.seed(5); d$z <- cumsum(rnorm(150))
  m3 <- mtnardl(y ~ x + z, d, thresholds = c(0.25, 0.5, 0.75), threshold_type = "quantile",
                decompose = "x")
  expect_equal(unname(m3$thresholds_by_var$x), unname(quantile(diff(d$x), c(0.25, 0.5, 0.75))))
  expect_true("z" %in% names(m3$long_run))
  expect_equal(length(m3$long_run), 5)
})

test_that("mtnardl bootstrap uses the recursive engine and reports t and Find", {
  d <- mt_sim(120, 4)
  m <- mtnardl(y ~ x, d, thresholds = 0, bootstrap = TRUE, nboot = 9, seed = 1)
  eng <- m$boot_results$engine
  expect_lt(eng$design_check, 1e-10)
  expect_lt(max(eng$dgpcheck), 1e-8)
  expect_true(all(c("F", "t", "Find") %in% names(m$bounds_test$p_values)))
})
