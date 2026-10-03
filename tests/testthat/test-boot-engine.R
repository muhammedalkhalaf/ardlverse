# Tests of the shared recursive bootstrap engine (R/ardl_boot_engine.R)

ps_dec <- function(x) {
  x <- as.matrix(x)
  d <- rbind(0, diff(x))
  cbind(apply(d * (d > 0), 2, cumsum), apply(d * (d < 0), 2, cumsum))
}

test_that("engine recursion reproduces the data with the original residuals", {
  set.seed(1)
  n <- 80
  y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
  tt <- 1:n
  Fo <- cbind(sin(2 * pi * tt / n), cos(2 * pi * tt / n))
  for (cs in c(2, 3, 4)) {
    r <- suppressWarnings(ardlverse:::.ardl_boot_engine(
      y, x, p = 2, q = 1, case = cs, det = Fo, decompose = ps_dec, B = 5, seed = 1))
    expect_lt(max(r$dgpcheck), 1e-10)
    expect_lt(r$design_check, 1e-10)
  }
  r <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 0, nulls = "joint",
                                     xmodel = "var", recentre = "once",
                                     init = "observed", B = 5, seed = 1)
  expect_lt(max(r$dgpcheck), 1e-10)
})

test_that("engine statistics equal restricted-RSS F and t tests by lm()", {
  set.seed(2)
  n <- 70
  y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
  r <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 1, B = 3, seed = 1)
  t <- 3:n
  dy <- diff(y)[t - 1]; ly <- y[t - 1]; lx <- x[t - 1]
  dyl <- y[t - 1] - y[t - 2]; dx0 <- x[t] - x[t - 1]; dx1 <- x[t - 1] - x[t - 2]
  u <- lm(dy ~ ly + lx + dyl + dx0 + dx1)
  Fh <- anova(lm(dy ~ dyl + dx0 + dx1), u)$F[2]
  Fi <- anova(lm(dy ~ ly + dyl + dx0 + dx1), u)$F[2]
  expect_equal(unname(r$statistic["Fov"]), Fh, tolerance = 1e-10)
  expect_equal(unname(r$statistic["Find"]), Fi, tolerance = 1e-10)
  expect_equal(unname(r$statistic["t"]), summary(u)$coefficients["ly", 3], tolerance = 1e-10)
})

test_that("order-statistic critical values follow BVZ eqs. 24-25", {
  d <- sample(1:199)
  expect_equal(ardlverse:::.abe_cv(d, 0.05), 190)
  expect_equal(ardlverse:::.abe_cv(d, 0.05, lower = TRUE), 10)
  cvF <- ardlverse:::.abe_cv(d, 0.05)
  expect_true(sum(d > cvF) <= 0.05 * 199 && sum(d > cvF - 1) > 0.05 * 199)
  cvt <- ardlverse:::.abe_cv(d, 0.05, lower = TRUE)
  expect_true(sum(d < cvt) <= 0.05 * 199 && sum(d < cvt + 1) > 0.05 * 199)
})

test_that("combined decision is Fov AND t AND Find, degenerate cases labelled", {
  dec <- ardlverse:::.abe_decision
  expect_equal(dec(c(Fov = TRUE, t = TRUE, Find = TRUE))$decision, "COINTEGRATION")
  expect_equal(dec(c(Fov = FALSE, t = TRUE, Find = TRUE))$decision, "NO_COINTEGRATION")
  expect_equal(dec(c(Fov = TRUE, t = FALSE, Find = TRUE))$decision, "DEGENERATE_1")
  expect_equal(dec(c(Fov = TRUE, t = TRUE, Find = FALSE))$decision, "DEGENERATE_2")
  expect_equal(dec(c(Fov = TRUE, t = TRUE, Find = NA), use_find = FALSE)$decision, "COINTEGRATION")
  expect_equal(dec(c(Fov = TRUE, t = NA, Find = TRUE))$decision, "UNDETERMINED")
})

test_that("failed replications are NA, counted and reported", {
  set.seed(3)
  n <- 60
  y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
  cnt <- 0
  sf <- function(yy, xx, sp) {
    cnt <<- cnt + 1
    if (cnt > 1 && cnt %% 3 == 0) stop("boom")
    ardlverse:::.abe_stats(ardlverse:::.abe_design(yy, xx, 1, 0, 3, 3), 3)
  }
  expect_warning(r <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 0, B = 9,
                                                    stat_fun = sf, seed = 1),
                 "failed and were set to NA")
  expect_true(any(is.na(r$boot[, "Fov"])))
  expect_equal(unname(r$n_fail["Fov"]), sum(is.na(r$boot[, "Fov"])))
  expect_false(any(r$boot == 0, na.rm = TRUE))
})

test_that("the seed is local: the caller's RNG state is restored", {
  set.seed(10)
  s0 <- .Random.seed
  d <- data.frame(y = cumsum(rnorm(50)), x = cumsum(rnorm(50)))
  s1 <- .Random.seed
  b1 <- boot_ardl(y ~ x, d, nboot = 9, seed = 5)
  expect_identical(.Random.seed, s1)
  b2 <- boot_ardl(y ~ x, d, nboot = 9, seed = 5)
  expect_identical(b1$boot_F, b2$boot_F)
})

test_that("bootstrap partial sums are rebuilt from x* and are monotone", {
  set.seed(4)
  n <- 60
  y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
  dec <- list(full = ps_dec, step = function(d) c(d * (d > 0), d * (d < 0)))
  r <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 0, decompose = dec, B = 3, seed = 1)
  expect_lt(max(r$dgpcheck), 1e-10)
  # one draw through the generator
  got <- NULL
  sf <- function(yy, xx, sp) {
    got <<- list(y = yy, x = xx)
    ardlverse:::.abe_stats(ardlverse:::.abe_design(yy, ps_dec(xx), 1, c(0, 0), 3, 3), 3)
  }
  ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 0, decompose = dec, B = 1,
                                stat_fun = sf, seed = 2)
  w <- ps_dec(got$x)
  expect_true(all(diff(w[, 1]) >= 0) && all(diff(w[, 2]) <= 0))
  expect_false(isTRUE(all.equal(got$x, x)))
})

test_that("the engine file carries its version header", {
  expect_equal(ardlverse:::.abe_engine_version, "1.1.0")
})

# Independent minimal recursion (adapted from the reviewer's BVZ/MSG script):
# case 3, level regressors w = x (linear) or the partial sums of x (NARDL),
# with the data's levels of w as initial values of a random block.
ind_design <- function(y, w, p, q, t0) {
  w <- as.matrix(w); n <- length(y); k <- ncol(w)
  lagv <- function(v, j) if (j == 0) v else c(rep(NA, j), head(v, -j))
  dy <- c(NA, diff(y)); dw <- rbind(NA, diff(w)); r <- t0:n
  Z <- cbind(1, c(NA, y[-n])[r])
  for (j in 1:k) Z <- cbind(Z, c(NA, w[-n, j])[r])
  if (p > 0) for (i in 1:p) Z <- cbind(Z, lagv(dy, i)[r])
  for (j in 1:k) for (l in 0:q[j]) Z <- cbind(Z, lagv(dw[, j], l)[r])
  list(Y = dy[r], Z = unname(Z), ily = 2, ilw = 2 + 1:k)
}
ind_stats <- function(y, w, p, q, t0) {
  D <- ind_design(y, w, p, q, t0); Y <- D$Y; Z <- D$Z
  f <- lm.fit(Z, Y); R <- sum(f$residuals^2); df <- length(Y) - ncol(Z)
  Fr <- function(drop) ((sum(lm.fit(Z[, -drop, drop = FALSE], Y)$residuals^2) - R) /
                          length(drop)) / (R / df)
  zo <- lm.fit(Z[, -D$ily, drop = FALSE], Z[, D$ily])$residuals
  c(Fov = Fr(c(D$ily, D$ilw)), t = unname(f$coefficients[D$ily] / sqrt(R / df / sum(zo^2))),
    Find = Fr(D$ilw))
}
ind_boot <- function(y, x, p, q, t0, B, seed, dec, step, nulls = "separate",
                     xlags = p, xvar = FALSE, block = TRUE, draw = TRUE) {
  x <- as.matrix(x); n <- length(y); k <- ncol(x)
  W <- dec(x); kw <- ncol(W)
  D <- ind_design(y, W, p, q, t0); N <- length(D$Y)
  drops <- list(Fov = c(D$ily, D$ilw), t = D$ily, Find = D$ilw)
  if (nulls == "joint") drops <- drops["Fov"]
  dy <- c(NA, diff(y)); dx <- rbind(NA, diff(x)); r <- t0:n
  XZ <- cbind(1, if (xvar) y[r - 1], x[r - 1, , drop = FALSE])
  for (l in seq_len(xlags)) XZ <- cbind(XZ, dy[r - l], dx[r - l, , drop = FALSE])
  fx <- lm.fit(XZ, dx[r, , drop = FALSE])
  Cx <- as.matrix(fx$coefficients); ex <- as.matrix(fx$residuals)
  set.seed(seed)
  out <- matrix(NA, B, 3, dimnames = list(NULL, c("Fov", "t", "Find")))
  for (s in names(drops)) {
    keep <- setdiff(seq_len(ncol(D$Z)), drops[[s]])
    fr <- lm.fit(D$Z[, keep, drop = FALSE], D$Y)
    by <- numeric(ncol(D$Z)); by[keep] <- fr$coefficients; e <- fr$residuals
    for (b in 1:B) {
      j <- sample.int(N, N, replace = TRUE)
      if (draw) {
        u <- e[j] - mean(e[j])
        v <- sweep(ex[j, , drop = FALSE], 2, colMeans(ex[j, , drop = FALSE]))
      } else {
        u <- (e - mean(e))[j]; v <- sweep(ex, 2, colMeans(ex))[j, , drop = FALSE]
      }
      rr <- if (block) sample.int(n - t0 + 2, 1) else 1
      h <- 1:(t0 - 1)
      ys <- numeric(n); xs <- matrix(0, n, k); ws <- matrix(0, n, kw)
      ys[h] <- y[rr + h - 1]; xs[h, ] <- x[rr + h - 1, ]; ws[h, ] <- W[rr + h - 1, ]
      for (t in t0:n) {
        i <- t - t0 + 1
        z <- c(1, if (xvar) ys[t - 1], xs[t - 1, ])
        for (l in seq_len(xlags)) z <- c(z, ys[t - l] - ys[t - l - 1], xs[t - l, ] - xs[t - l - 1, ])
        xs[t, ] <- xs[t - 1, ] + drop(z %*% Cx) + v[i, ]
        ws[t, ] <- ws[t - 1, ] + step(xs[t, ] - xs[t - 1, ])
        zy <- c(1, ys[t - 1], ws[t - 1, ])
        if (p > 0) for (l in 1:p) zy <- c(zy, ys[t - l] - ys[t - l - 1])
        for (jj in 1:kw) for (l in 0:q[jj]) zy <- c(zy, ws[t - l, jj] - ws[t - l - 1, jj])
        ys[t] <- ys[t - 1] + sum(by * zy) + u[i]
      }
      st <- ind_stats(ys, dec(xs), p, q, t0)
      if (nulls == "joint") out[b, ] <- st else out[b, s] <- st[s]
    }
  }
  out
}

test_that("bootstrap draws equal an independent recursion (BVZ, NARDL, block start)", {
  set.seed(41)
  n <- 70
  x <- cumsum(rnorm(n))
  y <- numeric(n)
  for (t in 2:n) y[t] <- 0.7 * y[t - 1] + 0.5 * x[t - 1] + rnorm(1)
  step <- function(d) c(d * (d > 0), d * (d < 0))
  dec <- list(full = ps_dec, step = step)
  e <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = c(1, 1), decompose = dec,
                                     B = 6, seed = 7)
  ib <- ind_boot(y, x, 1, c(1, 1), 3, 6, 7, ps_dec, step)
  expect_equal(unname(e$boot), unname(ib), tolerance = 1e-10)
})

test_that("bootstrap draws equal an independent recursion (MSG joint null, linear)", {
  set.seed(42)
  n <- 60
  y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
  id <- function(x) as.matrix(x)
  e <- ardlverse:::.ardl_boot_engine(y, x, p = 1, q = 1, nulls = "joint", xmodel = "var",
                                     init = "observed", recentre = "once", B = 5, seed = 3)
  ib <- ind_boot(y, x, 1, 1, 3, 5, 3, id, identity, nulls = "joint", xvar = TRUE,
                 block = FALSE, draw = FALSE)
  expect_equal(unname(e$boot), unname(ib), tolerance = 1e-10)
})

test_that("block initial values carry the data's partial-sum levels", {
  # with the y equation under the t null, w_{t-1} enters the recursion; if
  # the block restarted the partial sums at zero the first generated dy*
  # would differ from the restricted model evaluated at the data's levels
  set.seed(43)
  n <- 70
  x <- cumsum(rnorm(n)); y <- numeric(n)
  for (t in 2:n) y[t] <- 0.6 * y[t - 1] + 0.8 * max(x[t - 1] - x[1], 0) + rnorm(1)
  W <- ps_dec(x)
  m <- list(kx = 1, kw = 2, case = 3, det = NULL,
            iy = list(const = 1, trend = integer(0), det = integer(0), ly = 2,
                      lw = 3:4, dy = integer(0), dw = 5:6),
            ix = list(det0 = 1, ly = integer(0), lx = 2, dy = integer(0), dx = list()),
            lag_dw = c(0, 0), col_dw = 1:2, xlags = 0, decompose = ps_dec,
            step = function(d) c(d * (d > 0), d * (d < 0)))
  cy <- c(0.1, 0, 0.3, -0.2, 0.5, 0.4)   # t null: no y_{t-1}, w_{t-1} kept
  cx <- matrix(c(0, 0), 2, 1)
  r <- 30
  hh <- r + 0:1
  g <- ardlverse:::.abe_generate(n, 3, y[hh], matrix(x[hh]), cy, cx, rep(0, n - 2),
                                 matrix(0, n - 2, 1), m, w0 = W[hh, ])
  expect_equal(g$w[1:2, ], W[hh, ])
  # first step: x* stays at x[r + 1] (zero innovation), so dw*_3 = 0
  expect_equal(g$y[3] - g$y[2], 0.1 + sum(c(0.3, -0.2) * W[r + 1, ]))
  expect_true(abs(sum(c(0.3, -0.2) * W[r + 1, ])) > 0.1)
})
