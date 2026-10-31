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
  
  result <- suppressWarnings(hausman_test(mg, pmg))

  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_length(result$eigenvalues, length(mg$long_run))
  if (result$n_nonpositive == 0L) {
    expect_true(result$df > 0)
    expect_true(result$p.value >= 0 && result$p.value <= 1)
  } else {
    expect_true(is.na(result$statistic))
  }
})

## Regression tests for 2.1.1: non-positive-definite difference matrix.
## Up to 2.1.0 the function returned statistic 0 and p-value 1 and printed
## "Cannot reject H0" when V(MG) - V(eff) had a non-positive eigenvalue.
fake_fit <- function(est, b, V, sigma = 1) {
  structure(list(estimator = est, long_run = b, long_run_vcov = V, sigma = sigma),
            class = "panel_ardl")
}

test_that("hausman_test returns NA with a warning when the difference is not PD", {
  b_mg  <- c(1.0, 0.5); b_pmg <- c(0.8, 0.6)
  V_mg  <- diag(c(0.10, 0.01))
  V_pmg <- diag(c(0.05, 0.02))          # second eigenvalue of the difference < 0
  mg <- fake_fit("mg", b_mg, V_mg); pmg <- fake_fit("pmg", b_pmg, V_pmg)
  expect_warning(out <- capture.output(res <- hausman_test(mg, pmg, sigmamore = FALSE)),
                 "not positive definite")
  expect_true(is.na(res$statistic))
  expect_true(is.na(res$p.value))
  expect_equal(res$n_nonpositive, 1L)
  expect_false(any(grepl("Cannot reject H0", out, fixed = TRUE)))
  expect_true(any(grepl("Test not defined", out, fixed = TRUE)))
})

test_that("hausman_test method = 'ginv' uses the positive eigen-subspace", {
  b_mg  <- c(1.0, 0.5); b_pmg <- c(0.8, 0.6)
  V_mg  <- diag(c(0.10, 0.01)); V_pmg <- diag(c(0.05, 0.02))
  mg <- fake_fit("mg", b_mg, V_mg); pmg <- fake_fit("pmg", b_pmg, V_pmg)
  expect_warning(out <- capture.output(res <- hausman_test(mg, pmg, sigmamore = FALSE, method = "ginv")),
                 "generalised inverse")
  expect_equal(res$df, 1L)
  expect_equal(res$statistic, (1.0 - 0.8)^2 / 0.05, tolerance = 1e-12)
  expect_equal(res$p.value, pchisq((1.0 - 0.8)^2 / 0.05, 1, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("hausman_test is unchanged when the difference is positive definite", {
  b_mg  <- c(1.0, 0.5); b_pmg <- c(0.8, 0.6)
  V_mg  <- diag(c(0.10, 0.04)); V_pmg <- diag(c(0.05, 0.02))
  mg <- fake_fit("mg", b_mg, V_mg); pmg <- fake_fit("pmg", b_pmg, V_pmg)
  expect_silent(out <- capture.output(res <- hausman_test(mg, pmg, sigmamore = FALSE)))
  W <- V_mg - V_pmg; d <- b_mg - b_pmg
  expect_equal(res$statistic, as.numeric(t(d) %*% solve(W) %*% d))
  expect_equal(res$df, 2L)
  expect_equal(res$n_nonpositive, 0L)
})

test_that("hausman_test also accepts a DFE efficient model", {
  b_mg <- c(1.0, 0.5); b_dfe <- c(0.9, 0.55)
  mg <- fake_fit("mg", b_mg, diag(c(0.10, 0.04))); dfe <- fake_fit("dfe", b_dfe, diag(c(0.05, 0.02)))
  expect_silent(out <- capture.output(res <- hausman_test(mg, dfe, sigmamore = FALSE)))
  expect_true(is.finite(res$statistic))
})

## Tests added in 2.1.2: compatibility check, quiet, hand quadratic form,
## sigmamore scale, NA when not positive definite with the rank message.

# Small panel ARDL(1,1) DGP with homogeneous long-run slopes (no set.seed
# inside; the test sets the seed).
make_small_panel <- function(N, T, k, theta = seq(1, 0.5, length.out = k), burn = 50) {
  Tt <- T + burn
  out <- vector("list", N)
  for (i in seq_len(N)) {
    phi_i <- runif(1, -0.6, -0.2); delta_i <- runif(k, -0.3, 0.3); mu_i <- rnorm(1, 0, 0.5)
    X <- apply(matrix(rnorm(Tt * k), Tt, k), 2, cumsum); e <- rnorm(Tt)
    y <- numeric(Tt); y[1] <- sum(theta * X[1, ]) + e[1]
    for (t in 2:Tt) {
      y[t] <- y[t - 1] + mu_i + phi_i * (y[t - 1] - sum(theta * X[t, ])) +
        sum(delta_i * (X[t, ] - X[t - 1, ])) + e[t]
    }
    keep <- (burn + 1):Tt
    d <- data.frame(id = i, time = seq_len(T), y = y[keep])
    for (j in seq_len(k)) d[[paste0("x", j)]] <- X[keep, j]
    out[[i]] <- d
  }
  do.call(rbind, out)
}

test_that("hausman_test stops when the two fits are not comparable", {
  skip_on_cran()
  set.seed(20261030)
  d <- make_small_panel(N = 8, T = 30, k = 2)
  mg  <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", estimator = "mg")
  pmg <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", estimator = "pmg")
  # different regressors of the same length
  pmg_x <- panel_ardl(y ~ x2 + x1, d, id = "id", time = "time", estimator = "pmg")
  expect_error(hausman_test(mg, pmg_x, quiet = TRUE), "regressors")
  # different lag order
  pmg_p2 <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", p = 2, q = 2, estimator = "pmg")
  expect_error(hausman_test(mg, pmg_p2, quiet = TRUE), "lag order")
  # different data: fewer periods, same groups (nobs differs)
  pmg_sub <- panel_ardl(y ~ x1 + x2, d[d$time <= 25, ], id = "id", time = "time", estimator = "pmg")
  expect_error(hausman_test(mg, pmg_sub, quiet = TRUE), "number of observations")
  # different group count
  pmg_g <- panel_ardl(y ~ x1 + x2, d[d$id <= 6, ], id = "id", time = "time", estimator = "pmg")
  expect_error(hausman_test(mg, pmg_g, quiet = TRUE), "number of groups")
  # comparable fits pass the check
  expect_no_error(suppressWarnings(hausman_test(mg, pmg, quiet = TRUE)))
})

test_that("hausman_test quiet = TRUE prints nothing and returns the same list", {
  b_mg  <- c(1.0, 0.5); b_pmg <- c(0.8, 0.6)
  V_mg  <- diag(c(0.10, 0.04)); V_pmg <- diag(c(0.05, 0.02))
  mg <- fake_fit("mg", b_mg, V_mg); pmg <- fake_fit("pmg", b_pmg, V_pmg)
  out <- capture.output(res_q <- hausman_test(mg, pmg, sigmamore = FALSE, quiet = TRUE))
  expect_identical(out, character(0))
  out2 <- capture.output(res <- hausman_test(mg, pmg, sigmamore = FALSE))
  expect_true(length(out2) > 0)
  expect_identical(res_q, res)
  expect_error(hausman_test(mg, pmg, quiet = NA), "quiet")
})

test_that("hausman_test equals a hand quadratic form on real fits at 1e-10", {
  skip_on_cran()
  set.seed(20261030)
  d <- make_small_panel(N = 12, T = 30, k = 2)
  mg  <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", estimator = "mg")
  pmg <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", estimator = "pmg")
  dfe <- panel_ardl(y ~ x1 + x2, d, id = "id", time = "time", estimator = "dfe")
  # MG inputs rebuilt independently with lm(): theta_i = -beta_i / phi_i,
  # V(MG) = sum (theta_i - mean)(theta_i - mean)' / (N (N - 1))
  th <- t(sapply(split(d, d$id), function(g) {
    g <- g[order(g$time), ]
    df <- data.frame(dy = diff(g$y), ylag = g$y[-nrow(g)],
                     x1 = g$x1[-1], x2 = g$x2[-1],
                     dx1 = diff(g$x1), dx2 = diff(g$x2))
    cf <- coef(lm(dy ~ ylag + x1 + x2 + dx1 + dx2, df))
    -cf[c("x1", "x2")] / cf["ylag"]
  }))
  N <- nrow(th); theta_mg <- colMeans(th)
  V_mg <- crossprod(sweep(th, 2, theta_mg)) / (N * (N - 1))
  expect_equal(unname(mg$long_run), unname(theta_mg), tolerance = 1e-10)
  expect_equal(unname(mg$long_run_vcov), unname(V_mg), tolerance = 1e-10)
  for (eff in list(pmg, dfe)) {
    for (sm in c(FALSE, TRUE)) {
      res <- suppressWarnings(hausman_test(mg, eff, sigmamore = sm, quiet = TRUE))
      q <- mg$long_run - eff$long_run
      sc <- if (sm) (eff$sigma / mg$sigma)^4 else 1
      W <- sc * mg$long_run_vcov - eff$long_run_vcov
      expect_equal(res$scale, sc, tolerance = 1e-12)
      if (all(eigen(W, symmetric = TRUE)$values > 0)) {
        H <- as.numeric(t(q) %*% solve(W) %*% q)
        expect_equal(res$statistic, H, tolerance = 1e-10)
        expect_equal(res$df, length(q))
        expect_equal(res$p.value, pchisq(H, length(q), lower.tail = FALSE), tolerance = 1e-10)
      } else {
        expect_true(is.na(res$statistic))
      }
    }
  }
})

test_that("hausman_test sigmamore applies the fourth power of the RMSE ratio", {
  b_mg <- c(1.0, 0.5); b_pmg <- c(0.8, 0.6)
  V_mg <- diag(c(0.10, 0.04)); V_pmg <- diag(c(0.05, 0.02))
  mg <- fake_fit("mg", b_mg, V_mg, sigma = 1.0); pmg <- fake_fit("pmg", b_pmg, V_pmg, sigma = 1.1)
  res <- hausman_test(mg, pmg, sigmamore = TRUE, quiet = TRUE)
  W <- 1.1^4 * V_mg - V_pmg; q <- b_mg - b_pmg
  expect_equal(res$scale, 1.1^4, tolerance = 1e-12)
  expect_equal(res$statistic, as.numeric(t(q) %*% solve(W) %*% q), tolerance = 1e-12)
  res0 <- hausman_test(mg, pmg, sigmamore = FALSE, quiet = TRUE)
  expect_equal(res0$scale, 1)
})

test_that("hausman_test returns NA and names the rank bound when N - 1 < k", {
  # N = 3 groups, k = 3 long-run coefficients: V(MG) has rank at most 2
  set.seed(1)
  G <- matrix(rnorm(9), 3, 3); Gc <- sweep(G, 2, colMeans(G))
  V_mg <- crossprod(Gc) / (3 * 2)             # rank 2
  V_pmg <- diag(1e-3, 3)
  mg <- fake_fit("mg", c(1, 0.5, 0.2), V_mg); pmg <- fake_fit("pmg", c(0.9, 0.6, 0.1), V_pmg)
  mg$ec_group <- c(-0.3, -0.4, -0.5); mg$ngroups <- 3L; pmg$ngroups <- 3L
  expect_warning(res <- hausman_test(mg, pmg, sigmamore = FALSE, quiet = TRUE),
                 "rank at most N - 1 = 2 < k = 3")
  expect_true(is.na(res$statistic)); expect_true(is.na(res$df)); expect_true(is.na(res$p.value))
  expect_true(res$n_nonpositive >= 1L)
  out <- capture.output(suppressWarnings(hausman_test(mg, pmg, sigmamore = FALSE)))
  expect_true(any(grepl("rank at most N - 1 = 2 < k = 3", out, fixed = TRUE)))
  expect_true(any(grepl("Test not defined", out, fixed = TRUE)))
})

test_that("rank-deficient warnings do not suggest options that cannot help", {
  set.seed(1)
  G <- matrix(rnorm(9), 3, 3); Gc <- sweep(G, 2, colMeans(G))
  V_mg <- crossprod(Gc) / (3 * 2)
  V_pmg <- diag(1e-3, 3)
  mg <- fake_fit("mg", c(1, 0.5, 0.2), V_mg); pmg <- fake_fit("pmg", c(0.9, 0.6, 0.1), V_pmg)
  mg$ec_group <- c(-0.3, -0.4, -0.5); mg$ngroups <- 3L; pmg$ngroups <- 3L
  w <- character(0)
  withCallingHandlers(capture.output(hausman_test(mg, pmg, sigmamore = FALSE, quiet = TRUE)),
                      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("rank at most N - 1 = 2", w, fixed = TRUE)))
  expect_false(any(grepl("Options:", w, fixed = TRUE)))
  # ginv: statistic reported on a subspace of at most N - 1 dimensions, and said so
  w2 <- character(0)
  res <- withCallingHandlers(hausman_test(mg, pmg, sigmamore = FALSE, method = "ginv", quiet = TRUE),
                             warning = function(cnd) { w2 <<- c(w2, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("subspace of at most N - 1 = 2 dimensions", w2, fixed = TRUE)))
  expect_true(res$df <= 2L)
  out <- capture.output(suppressWarnings(hausman_test(mg, pmg, sigmamore = FALSE, method = "ginv")))
  expect_true(any(grepl("reduced-rank statistic", out, fixed = TRUE)))
})
