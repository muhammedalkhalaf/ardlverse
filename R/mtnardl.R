#' @title Multiple-Threshold Nonlinear ARDL (MT-NARDL)
#' @description Extends NARDL to allow multiple threshold decomposition
#' for capturing complex asymmetric relationships.
#'
#' @details
#' Each decomposed regressor \eqn{x} is split into regimes according to the
#' size of its change \eqn{\Delta x_t}: with thresholds
#' \eqn{c_1 < \dots < c_{R-1}}, regime \eqn{r} collects \eqn{\Delta x_t} when
#' \eqn{c_{r-1} < \Delta x_t \le c_r} (\eqn{c_0 = -\infty},
#' \eqn{c_R = \infty}), and its partial sum is the cumulative sum of these
#' changes, starting at zero at the first observation. For example, with
#' thresholds \code{c(-0.02, 0, 0.02)} a variable is decomposed into large
#' decreases, small decreases, small increases and large increases. With
#' \code{threshold_type = "quantile"} the thresholds are probabilities and are
#' converted to the sample quantiles of \eqn{\Delta x_t} of each decomposed
#' variable (a quantile partition of the changes). Regressors not listed in
#' \code{decompose} enter linearly.
#'
#' The conditional error correction model is
#' \deqn{\Delta y_t = \rho y_{t-1} + \sum_m \theta_m w_{m,t-1} +
#'   \sum_{i=1}^{p} \phi_i \Delta y_{t-i} +
#'   \sum_m \sum_{l=0}^{q-1} \pi_{ml} \Delta w_{m,t-l} + d_t + u_t,}
#' where \eqn{w} are the regime partial sums and the linear regressors and
#' \eqn{d_t} holds the deterministic terms of the case.
#'
#' \strong{Bounds test.} The statistics are the F test on all lagged levels
#' (Fov), the t test on \eqn{y_{t-1}} and the F test on the lagged
#' regressors (Find). The Kripfganz and Schneider (2020) bounds are computed
#' with \code{k} equal to the number of level regressors (partial sums plus
#' linear regressors) by default (\code{bounds_k = "partial_sums"}), or the
#' number of original regressors (\code{bounds_k = "original"}). This is a
#' package choice pending verification against Pal and Mitra (2016) and Shin,
#' Yu and Greenwood-Nimmo (2014); in a Monte Carlo with independent random
#' walks the first choice gave sizes close to the nominal level. The bounds
#' decision requires both F and t to exceed their I(1) bounds. With
#' \code{bootstrap = TRUE} the recursive bootstrap of \code{\link{boot_ardl}}
#' is used (separate nulls of Bertelli, Vacca and Zoia (2022) by default):
#' x* is generated recursively, the regime partial sums are rebuilt from x*
#' with the same thresholds, and the decision requires Fov, t and Find to
#' reject (degenerate cases labelled).
#'
#' \strong{Long run.} The long-run coefficients \eqn{-\theta_m / \rho} have
#' delta-method standard errors computed with the full covariance matrix.
#' For each decomposed variable the joint Wald test of equal long-run
#' coefficients across regimes is the linear test \eqn{\theta_1 = \dots =
#' \theta_R} (chi-squared with R - 1 degrees of freedom); pairwise tests and
#' the delta-method differences of the long-run coefficients are also
#' reported.
#'
#' \strong{Dynamic multipliers.} \code{multipliers$cumulative} gives the
#' response of the level of y at horizons 0 to \code{horizon} to a permanent
#' unit increase in each partial sum (or linear regressor) at horizon 0,
#' computed by the full recursion of the estimated ECM; it converges to the
#' long-run coefficient. The definition of the multiplier follows the
#' cumulative multiplier of the NARDL literature and is pending verification
#' against Shin, Yu and Greenwood-Nimmo (2014).
#'
#' @param formula A formula such as \code{gasoline ~ oil_price}; only plain
#'   variable names are accepted (create transformed variables in
#'   \code{data}).
#' @param data A data frame containing the time series (no missing values).
#' @param thresholds Numeric vector of threshold values (default: c(0)), or
#'   probabilities when \code{threshold_type = "quantile"}.
#' @param p Integer. Number of lagged differences of y (default: 1)
#' @param q Integer or vector (one entry per regressor). Number of
#'   differences of each regressor, lags 0 to \code{q - 1} (default: 1)
#' @param case Integer from 1-5 specifying deterministic components (default: 3)
#' @param auto_select Logical. Automatically select thresholds (default:
#'   FALSE). With \code{n_thresholds = 1} the threshold minimising the AIC
#'   over 0 and the deciles of the changes of the first decomposed variable is
#'   chosen; with more thresholds, equally spaced quantiles of those changes
#'   are used.
#' @param n_thresholds Integer. Number of thresholds to select if auto_select = TRUE
#' @param bootstrap Logical. Use bootstrap inference (default: FALSE)
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param seed Random seed; used locally, the session's random number
#'   generator state is restored afterwards.
#' @param decompose Character vector of regressors to decompose (default:
#'   all regressors). The other regressors enter linearly.
#' @param threshold_type \code{"value"} (default) or \code{"quantile"}.
#' @param bounds_k \code{"partial_sums"} (default) or \code{"original"}; the
#'   k used for the bounds, see Details.
#' @param nulls Bootstrap nulls, \code{"separate"} (default) or
#'   \code{"joint"}; see \code{\link{boot_ardl}}.
#' @param horizon Largest horizon of the dynamic multipliers (default: 30).
#'
#' @return An object of class "mtnardl" containing:
#' \itemize{
#'   \item \code{model}: The estimated MT-NARDL model
#'   \item \code{F_stat}, \code{t_stat}, \code{Find_stat}: bounds test
#'     statistics
#'   \item \code{bounds_test}: decision, message, method (and bootstrap
#'     p-values)
#'   \item \code{critical_values}: Kripfganz and Schneider (2020) bounds,
#'     with \code{bounds_k_value} the k used
#'   \item \code{long_run}: Long-run coefficients for each regime;
#'     \code{long_run_table} adds delta-method standard errors
#'   \item \code{thresholds}: Threshold values used (a list per decomposed
#'     variable in \code{thresholds_by_var})
#'   \item \code{asymmetry_tests}: for each decomposed variable, pairwise
#'     Wald tests and the joint test (element \code{joint})
#'   \item \code{lr_differences}: delta-method differences of the long-run
#'     coefficients between regimes
#'   \item \code{multipliers}: long-run multipliers and cumulative dynamic
#'     multipliers (horizons 0 to \code{horizon})
#' }
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' Pal, D. and Mitra, S. K. (2016). Asymmetric oil product pricing in
#' India: Evidence from a multiple threshold nonlinear ARDL model.
#' \emph{Economic Modelling}, 59, 314-328.
#' \doi{10.1016/j.econmod.2016.08.003}
#'
#' Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In R. C. Sickles and W. C. Horrace (Eds.), \emph{Festschrift in Honor of
#' Peter Schmidt} (pp. 281-314). Springer. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' @examples
#' data <- generate_oil_data(n = 120)
#' result1 <- mtnardl(gasoline ~ oil_price, data = data)
#' result1
#' \donttest{
#' result2 <- mtnardl(gasoline ~ oil_price, data = data,
#'                    thresholds = c(-1, 0, 1))
#' summary(result2)
#' result3 <- mtnardl(gasoline ~ oil_price, data = data,
#'                    thresholds = c(1/3, 2/3), threshold_type = "quantile",
#'                    bootstrap = TRUE, nboot = 49, seed = 1)
#' result3$bounds_test
#' }
#'
#' @export
mtnardl <- function(formula, data, thresholds = c(0), p = 1, q = 1, case = 3,
                    auto_select = FALSE, n_thresholds = 2,
                    bootstrap = FALSE, nboot = 2000, seed = NULL,
                    decompose = NULL,
                    threshold_type = c("value", "quantile"),
                    bounds_k = c("partial_sums", "original"),
                    nulls = c("separate", "joint"), horizon = 30) {

  threshold_type <- match.arg(threshold_type)
  bounds_k <- match.arg(bounds_k)
  nulls <- match.arg(nulls)
  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  if (p < 1) stop("'p' must be at least 1")

  fv <- .ardl_formula_vars(formula, data)
  y_var <- fv$y_var
  x_vars <- fv$x_vars
  k <- length(x_vars)
  if (is.null(decompose)) decompose <- x_vars
  if (!all(decompose %in% x_vars))
    stop("all variables in 'decompose' must be regressors in the formula")
  decompose <- x_vars[x_vars %in% decompose]
  if (length(q) == 1) q <- rep(q, k)
  if (length(q) != k) stop("'q' must have length 1 or one entry per regressor")
  if (any(q < 1)) stop("'q' must be at least 1 (lags 0 to q - 1 of the differences)")
  names(q) <- x_vars

  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  n <- length(y)

  if (auto_select) {
    thresholds <- .select_optimal_thresholds(y, X, n_thresholds, p, q, case,
                                             x_vars, decompose)
    threshold_type <- "value"
  }
  thr <- .mt_threshold_list(X[, decompose, drop = FALSE], thresholds, threshold_type)
  spec <- list(x_vars = x_vars, decompose = decompose, thr = thr, q = q,
               p = p, case = case)
  fit <- .mtnardl_fit(y, X, spec)
  design <- fit$design
  regime_names <- fit$level_names
  k_total <- length(regime_names)
  n_regimes <- length(thr[[1]]) + 1

  model <- stats::lm(fit$dy ~ design - 1)
  coefs <- stats::coef(model)
  names(coefs) <- colnames(design)
  vcov_mat <- stats::vcov(model)
  dimnames(vcov_mat) <- list(colnames(design), colnames(design))
  ec_coef <- coefs[1]

  st <- .aardl_stats(fit$dy, design, k_total, case)
  F_stat <- unname(st["Fov"])
  t_stat <- unname(st["t"])
  Find_stat <- unname(st["Find"])

  # Kripfganz and Schneider (2020) bounds
  k_b <- if (bounds_k == "partial_sums") k_total else k
  sr <- ncol(design) - 1 - k_total - (case >= 2) - (case >= 4)
  cv <- pss_critical_values(k_b, case, n = nrow(design), sr = sr)
  cv$bounds_k_value <- k_b

  boot_results <- NULL
  if (bootstrap) {
    stat_fun <- function(yy, xx, sp) {
      ff <- .mtnardl_fit(yy, xx, spec)
      .aardl_stats(ff$dy, ff$design, k_total, case)
    }
    eng <- .ardl_boot_engine(y, X, p = p, q = fit$q_levels - 1, case = case,
                             t0 = max(p + 1, max(q)) + 1,
                             decompose = .mt_decomposer(spec), nulls = nulls,
                             xmodel = if (nulls == "separate") "vecm" else "var",
                             B = nboot,
                             init = if (nulls == "separate") "block" else "observed",
                             recentre = if (nulls == "separate") "draw" else "once",
                             stat_fun = stat_fun, use_find = TRUE, seed = seed)
    boot_results <- list(
      F_dist = eng$boot[, "Fov"],
      t_dist = eng$boot[, "t"],
      Find_dist = eng$boot[, "Find"],
      cv_F = stats::setNames(eng$cv[, "Fov"], c("90%", "95%", "97.5%", "99%")),
      cv_t = stats::setNames(eng$cv[, "t"], c("10%", "5%", "2.5%", "1%")),
      cv_Find = stats::setNames(eng$cv[, "Find"], c("90%", "95%", "97.5%", "99%")),
      engine = eng
    )
  }
  bounds_test <- .mtnardl_bounds_conclusion(F_stat, t_stat, cv, boot_results)

  # Long-run coefficients with delta-method standard errors
  ilev <- 1 + seq_len(k_total)
  if (abs(ec_coef) > 1e-10) {
    ld <- .lr_delta(coefs, vcov_mat, 1, ilev)
    lr_coefs <- stats::setNames(ld$lr, regime_names)
    zz <- ld$lr / ld$se
    lr_table <- data.frame(estimate = ld$lr, std_error = ld$se, z = zz,
                           p_value = 2 * stats::pnorm(-abs(zz)),
                           row.names = regime_names)
    lr_by_var <- list()
    for (v in x_vars) lr_by_var[[v]] <- lr_coefs[fit$var_of_col == v]
    VL <- ld$V
  } else {
    lr_coefs <- stats::setNames(rep(NA_real_, k_total), regime_names)
    lr_table <- NULL
    lr_by_var <- NULL
    VL <- NULL
  }

  asym <- .test_regime_asymmetry(coefs, vcov_mat, lr_coefs, VL, fit, decompose)

  multipliers <- .compute_mt_multipliers(coefs, fit, regime_names, horizon)

  fit_stats <- list(
    R2 = summary(model)$r.squared,
    adj_R2 = summary(model)$adj.r.squared,
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    sigma = summary(model)$sigma,
    df = model$df.residual
  )

  result <- list(
    model = model,
    bounds_test = bounds_test,
    F_stat = F_stat,
    t_stat = t_stat,
    Find_stat = Find_stat,
    critical_values = cv,
    boot_results = boot_results,
    long_run = lr_coefs,
    long_run_table = lr_table,
    long_run_by_var = lr_by_var,
    short_run = coefs,
    thresholds = thr[[1]],
    thresholds_by_var = thr,
    n_regimes = n_regimes,
    regime_names = regime_names,
    decompose = decompose,
    asymmetry_tests = asym$tests,
    lr_differences = asym$lr_differences,
    multipliers = multipliers,
    fit = fit_stats,
    call = match.call(),
    case = case,
    n = length(fit$dy),
    k = k,
    p = p,
    q = unname(q)
  )

  class(result) <- "mtnardl"
  return(result)
}


# Thresholds per decomposed variable (values, or quantiles of the changes)
.mt_threshold_list <- function(Xd, thresholds, type) {
  thresholds <- sort(unique(as.numeric(thresholds)))
  if (length(thresholds) == 0) stop("at least one threshold is needed")
  out <- list()
  for (v in colnames(Xd)) {
    if (type == "quantile") {
      if (any(thresholds <= 0 | thresholds >= 1))
        stop("with threshold_type = 'quantile' the thresholds must be probabilities in (0, 1)")
      out[[v]] <- unname(stats::quantile(diff(Xd[, v]), thresholds, type = 7))
    } else {
      out[[v]] <- thresholds
    }
  }
  out
}


#' @title Multiple-Threshold Decomposition
#' @description Regime partial sums: regime r of a column collects its change
#'   when it lies between thresholds r - 1 (excluded) and r (included), with the change at the first observation
#'   set to zero. \code{thresholds} is a numeric vector used for every
#'   column, or a list with one vector per column.
#' @keywords internal
.mt_decompose <- function(X, thresholds) {
  X <- as.matrix(X)
  k <- ncol(X)
  n <- nrow(X)
  if (!is.list(thresholds)) thresholds <- rep(list(thresholds), k)
  n_regimes <- length(thresholds[[1]]) + 1
  result <- matrix(0, n, k * n_regimes)
  for (j in 1:k) {
    dx <- c(0, diff(X[, j]))
    lo <- c(-Inf, thresholds[[j]])
    hi <- c(thresholds[[j]], Inf)
    for (r in 1:n_regimes) {
      result[, (j - 1) * n_regimes + r] <- cumsum(dx * (dx > lo[r] & dx <= hi[r]))
    }
  }
  result
}


# Level regressors of mtnardl(): regime partial sums of the decomposed
# variables and the linear regressors, in the order of x_vars
.mt_levels <- function(X, spec) {
  W <- NULL
  nm <- character(0)
  var_of_col <- character(0)
  qq <- integer(0)
  for (v in spec$x_vars) {
    if (v %in% spec$decompose) {
      Wv <- .mt_decompose(X[, v, drop = FALSE], list(spec$thr[[v]]))
      W <- cbind(W, Wv)
      nm <- c(nm, .make_regime_names(v, spec$thr[[v]]))
      var_of_col <- c(var_of_col, rep(v, ncol(Wv)))
      qq <- c(qq, rep(spec$q[[v]], ncol(Wv)))
    } else {
      W <- cbind(W, X[, v])
      nm <- c(nm, v)
      var_of_col <- c(var_of_col, v)
      qq <- c(qq, spec$q[[v]])
    }
  }
  colnames(W) <- nm
  list(W = W, names = nm, var_of_col = var_of_col, q = qq)
}

.mt_decomposer <- function(spec) {
  steps <- lapply(spec$x_vars, function(v) {
    if (v %in% spec$decompose) {
      lo <- c(-Inf, spec$thr[[v]])
      hi <- c(spec$thr[[v]], Inf)
      function(d) d * (d > lo & d <= hi)
    } else {
      function(d) d
    }
  })
  list(full = function(x) {
         x <- as.matrix(x)
         colnames(x) <- spec$x_vars
         .mt_levels(x, spec)$W
       },
       step = function(d) unlist(lapply(seq_along(steps), function(j) steps[[j]](d[j]))))
}


# Design of mtnardl(): dy_t on y_{t-1}, w_{t-1}, p lagged dy, dw lags 0..q-1
# and the deterministic terms; sample t = max(p + 1, max(q)) + 1, ..., n
.mtnardl_fit <- function(y, X, spec) {
  colnames(X) <- spec$x_vars
  L <- .mt_levels(X, spec)
  W <- L$W
  p <- spec$p
  case <- spec$case
  n <- length(y)
  max_lag <- max(p + 1, max(spec$q))
  valid_idx <- (max_lag + 1):n
  n_valid <- length(valid_idx)
  if (n_valid <= 5) stop("Insufficient observations for specified lag structure")
  dyv <- diff(y)
  dy <- dyv[(max_lag):(n - 1)]
  y_lag <- y[valid_idx - 1]
  dy_lags <- sapply(1:p, function(i) dyv[(max_lag - i):(n - 1 - i)])
  dy_lags <- matrix(dy_lags, n_valid, p)
  colnames(dy_lags) <- paste0("dy_l", 1:p)
  x_levels <- W[valid_idx - 1, , drop = FALSE]
  x_diffs <- NULL
  dn <- character(0)
  for (j in seq_len(ncol(W))) {
    dw <- diff(W[, j])
    for (i in 0:(L$q[j] - 1)) {
      x_diffs <- cbind(x_diffs, dw[(max_lag - i):(n - 1 - i)])
      dn <- c(dn, paste0("d_", L$names[j], "_l", i))
    }
  }
  colnames(x_diffs) <- dn
  design <- cbind(y_lag = y_lag, x_levels, dy_lags, x_diffs)
  if (case >= 2) design <- cbind(design, intercept = 1)
  if (case >= 4) design <- cbind(design, trend = 1:n_valid)
  list(dy = dy, design = design, level_names = L$names,
       var_of_col = L$var_of_col, q_levels = L$q)
}


#' @title Make Regime Names
#' @keywords internal
.make_regime_names <- function(x_vars, thresholds) {
  n_regimes <- length(thresholds) + 1
  names_list <- c()

  for (v in x_vars) {
    for (r in 1:n_regimes) {
      if (r == 1) {
        name <- paste0(v, "_r1_le", round(thresholds[1], 3))
      } else if (r == n_regimes) {
        name <- paste0(v, "_r", r, "_gt", round(thresholds[length(thresholds)], 3))
      } else {
        name <- paste0(v, "_r", r, "_", round(thresholds[r-1], 3), "_to_",
                       round(thresholds[r], 3))
      }
      names_list <- c(names_list, name)
    }
  }

  return(names_list)
}


#' @title Select Optimal Thresholds
#' @description With one threshold, the value minimising the AIC of the
#'   MT-NARDL model over 0 and the deciles of the changes of the first
#'   decomposed variable; with more, equally spaced quantiles of those
#'   changes.
#' @keywords internal
.select_optimal_thresholds <- function(y, X, n_thresholds, p, q, case,
                                       x_vars, decompose) {
  dx <- diff(X[, decompose[1]])
  candidates <- stats::quantile(dx, seq(0.1, 0.9, by = 0.1))

  if (n_thresholds == 1) {
    candidates <- sort(unique(c(0, candidates)))
    best_aic <- Inf
    best_threshold <- NA_real_
    for (th in candidates) {
      spec <- list(x_vars = x_vars, decompose = decompose,
                   thr = stats::setNames(rep(list(th), length(decompose)), decompose),
                   q = q, p = p, case = case)
      aic <- tryCatch({
        ff <- .mtnardl_fit(y, X, spec)
        m <- stats::lm.fit(ff$design, ff$dy)
        if (m$rank < ncol(ff$design)) NA_real_ else {
          nn <- length(ff$dy)
          rss <- sum(m$residuals^2)
          nn * (log(2 * pi) + log(rss / nn) + 1) + 2 * (ncol(ff$design) + 1)
        }
      }, error = function(e) NA_real_)
      if (is.finite(aic) && aic < best_aic) {
        best_aic <- aic
        best_threshold <- th
      }
    }
    if (is.na(best_threshold)) stop("threshold selection failed for every candidate")
    return(unname(best_threshold))
  }

  threshold_quantiles <- seq(1/(n_thresholds + 1),
                             n_thresholds/(n_thresholds + 1),
                             length.out = n_thresholds)
  as.numeric(stats::quantile(dx, threshold_quantiles))
}


#' @title MT-NARDL Bounds Test Conclusion
#' @description Bootstrap decision (Fov, t and Find must all reject) or
#'   Kripfganz and Schneider (2020) bounds decision on F and t.
#' @keywords internal
.mtnardl_bounds_conclusion <- function(F_stat, t_stat, cv, boot = NULL) {

  if (!is.null(boot)) {
    eng <- boot$engine
    return(list(decision = eng$decision,
                message = paste0(eng$label, " (bootstrap, ", 100 * eng$level, "% level)"),
                p_values = c(F = unname(eng$p_value["Fov"]), t = unname(eng$p_value["t"]),
                             Find = unname(eng$p_value["Find"])),
                method = "bootstrap"))
  }

  F_upper <- cv$F_bounds$I1
  F_lower <- cv$F_bounds$I0
  t_upper <- cv$t_bounds$I1
  t_lower <- cv$t_bounds$I0

  if (any(is.na(c(F_upper, F_lower, t_upper, t_lower, F_stat, t_stat)))) {
    return(list(decision = "INCONCLUSIVE",
                message   = "Critical values or statistics unavailable",
                method    = "asymptotic"))
  }

  if (F_stat > F_upper && t_stat < t_upper) {
    decision <- "COINTEGRATION"
    message <- "Cointegration: F > I(1) bound and t < I(1) bound"
  } else if (F_stat < F_lower || t_stat > t_lower) {
    decision <- "NO_COINTEGRATION"
    message <- "No cointegration: F < I(0) bound or t > I(0) bound"
  } else {
    decision <- "INCONCLUSIVE"
    message <- "Inconclusive: F or t between the I(0) and I(1) bounds"
  }

  list(decision = decision, message = message, method = "asymptotic")
}


#' @title Test Regime Asymmetry
#' @description Pairwise and joint Wald tests of equal level coefficients
#'   (equivalently equal long-run coefficients) across the regimes of each
#'   decomposed variable, and delta-method differences of the long-run
#'   coefficients.
#' @keywords internal
.test_regime_asymmetry <- function(coefs, vcov_mat, lr_coefs, VL, fit, decompose) {
  tests <- list()
  lrd <- NULL
  for (v in decompose) {
    idx <- which(fit$var_of_col == v)
    R_ <- length(idx)
    if (R_ < 2) next
    coef_idx <- idx + 1
    var_tests <- list()
    for (i in 1:(R_ - 1)) {
      for (j in (i + 1):R_) {
        Rm <- matrix(0, 1, length(coefs))
        Rm[1, coef_idx[i]] <- 1
        Rm[1, coef_idx[j]] <- -1
        wald <- .wald_lin(coefs, vcov_mat, Rm)
        p_value <- stats::pchisq(wald, 1, lower.tail = FALSE)
        var_tests[[paste0("regime", i, "_vs_regime", j)]] <- list(
          wald = wald, df = 1, p_value = p_value,
          diff = unname(coefs[coef_idx[i]] - coefs[coef_idx[j]]),
          significant = p_value < 0.05)
        if (!is.null(VL)) {
          cc <- numeric(length(lr_coefs))
          cc[idx[i]] <- 1
          cc[idx[j]] <- -1
          dl <- sum(cc * lr_coefs)
          se <- sqrt(drop(t(cc) %*% VL %*% cc))
          lrd <- rbind(lrd, data.frame(variable = v, comparison = paste0("regime", i, " - regime", j),
                                       difference = dl, std_error = se, wald = (dl / se)^2,
                                       p_value = stats::pchisq((dl / se)^2, 1, lower.tail = FALSE)))
        }
      }
    }
    Rj <- matrix(0, R_ - 1, length(coefs))
    for (r in 2:R_) {
      Rj[r - 1, coef_idx[r]] <- 1
      Rj[r - 1, coef_idx[1]] <- -1
    }
    wj <- .wald_lin(coefs, vcov_mat, Rj)
    pj <- stats::pchisq(wj, R_ - 1, lower.tail = FALSE)
    var_tests[["joint"]] <- list(wald = wj, df = R_ - 1, p_value = pj,
                                 diff = NA_real_, significant = pj < 0.05)
    tests[[v]] <- var_tests
  }
  list(tests = tests, lr_differences = lrd)
}


#' @title Compute Multiple-Threshold Dynamic Multipliers
#' @description Cumulative dynamic multipliers by the full recursion of the
#'   estimated ECM (see \code{.ecm_multiplier}).
#' @keywords internal
.compute_mt_multipliers <- function(coefs, fit, regime_names, horizons = 30) {
  k_total <- length(regime_names)
  rho <- coefs[1]
  lr_mult <- if (abs(rho) > 1e-10) -coefs[1 + seq_len(k_total)] / rho else rep(NA_real_, k_total)
  names(lr_mult) <- regime_names
  nm <- names(coefs)
  phi <- coefs[grep("^dy_l[0-9]+$", nm)]
  cum_mult <- matrix(NA_real_, horizons + 1, k_total, dimnames = list(NULL, regime_names))
  for (j in seq_len(k_total)) {
    pij <- coefs[paste0("d_", regime_names[j], "_l", 0:(fit$q_levels[j] - 1))]
    cum_mult[, j] <- .ecm_multiplier(rho, coefs[1 + j], phi, pij, horizons)
  }
  list(long_run = lr_mult, cumulative = cum_mult, horizons = 0:horizons)
}


#' @rdname mtnardl
#' @param x,object An object of class "mtnardl"
#' @param ... Not used
#' @export
print.mtnardl <- function(x, ...) {
  cat("\n")
  cat("Multiple-Threshold Nonlinear ARDL (MT-NARDL)\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Thresholds:", paste(round(x$thresholds, 4), collapse = ", "), "\n")
  cat("Regimes:", x$n_regimes, "\n")
  cat("Observations:", x$n, "\n\n")

  cat("Bounds Test:\n")
  cat(sprintf("  F-statistic: %.4f\n", x$F_stat))
  cat(sprintf("  t-statistic: %.4f\n", x$t_stat))
  if (!is.null(x$Find_stat)) cat(sprintf("  Find:        %.4f\n", x$Find_stat))
  cat("\nDecision:", x$bounds_test$decision, "\n")
  cat(x$bounds_test$message, "\n")

  invisible(x)
}


#' @rdname mtnardl
#' @export
summary.mtnardl <- function(object, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("    Multiple-Threshold NARDL Estimation Results\n")
  cat("=======================================================\n\n")

  cat("Model Specification:\n")
  cat("  Case:", object$case, "\n")
  for (v in names(object$thresholds_by_var))
    cat("  Thresholds (", v, "): ",
        paste(round(object$thresholds_by_var[[v]], 4), collapse = ", "), "\n", sep = "")
  cat("  Number of regimes:", object$n_regimes, "\n")
  cat("  Lags: p =", object$p, ", q =", paste(object$q, collapse = ","), "\n")
  cat("  Sample size:", object$n, "\n\n")

  cat("Bounds Test:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  F-statistic: %10.4f\n", object$F_stat))
  cat(sprintf("  t-statistic: %10.4f\n", object$t_stat))
  if (!is.null(object$Find_stat)) cat(sprintf("  Find:        %10.4f\n", object$Find_stat))
  if (identical(object$bounds_test$method, "bootstrap")) {
    pv <- object$bounds_test$p_values
    cat(sprintf("  Bootstrap p-values: F %.4f, t %.4f, Find %.4f\n", pv["F"], pv["t"], pv["Find"]))
  } else {
    cv <- object$critical_values
    cat(sprintf("  5%% bounds (k = %d): F [%.3f, %.3f], t [%.3f, %.3f]\n",
                cv$bounds_k_value, cv$F_bounds$I0, cv$F_bounds$I1,
                cv$t_bounds$I0, cv$t_bounds$I1))
  }
  cat("\n  Decision:", object$bounds_test$decision, "\n")
  cat(" ", object$bounds_test$message, "\n\n")

  cat("Long-Run Coefficients by Regime (delta-method s.e.):\n")
  cat("-------------------------------------------------------\n")
  if (!is.null(object$long_run_table)) print(round(object$long_run_table, 4))
  else print(round(object$long_run, 4))

  cat("\nAsymmetry Tests (Wald, level coefficients):\n")
  cat("-------------------------------------------------------\n")
  for (v in names(object$asymmetry_tests)) {
    cat("\n", v, ":\n")
    for (test_name in names(object$asymmetry_tests[[v]])) {
      test <- object$asymmetry_tests[[v]][[test_name]]
      sig <- if (isTRUE(test$significant)) "*" else ""
      cat(sprintf("  %-20s Wald = %7.3f, df = %d, p = %.4f %s\n",
                  test_name, test$wald, as.integer(test$df), test$p_value, sig))
    }
  }

  cat("\nModel Fit:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  R-squared:     %.4f\n", object$fit$R2))
  cat(sprintf("  Adj R-squared: %.4f\n", object$fit$adj_R2))
  cat(sprintf("  AIC:           %.2f\n", object$fit$AIC))
  cat(sprintf("  BIC:           %.2f\n", object$fit$BIC))

  cat("\n=======================================================\n\n")

  invisible(object)
}


#' @rdname mtnardl
#' @param type For \code{plot}: \code{"multipliers"} or \code{"asymmetry"}.
#' @export
plot.mtnardl <- function(x, type = c("multipliers", "asymmetry"), ...) {
  type <- match.arg(type)

  if (type == "multipliers") {
    mult <- x$multipliers$cumulative
    horizons <- x$multipliers$horizons

    n_vars <- ncol(mult)
    colors <- grDevices::rainbow(n_vars)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow = c(1, 1))
    graphics::matplot(horizons, mult, type = "l", lty = 1, col = colors,
                     xlab = "Horizon", ylab = "Cumulative Multiplier",
                     main = "Dynamic Multipliers by Regime")
    graphics::legend("topright", legend = colnames(mult), col = colors,
                    lty = 1, cex = 0.7)
    graphics::abline(h = 0, lty = 2, col = "gray")
  }

  if (type == "asymmetry") {
    lr <- x$long_run

    graphics::barplot(lr, las = 2, col = "steelblue",
                     main = "Long-Run Coefficients by Regime",
                     ylab = "Coefficient")
    graphics::abline(h = 0, lty = 2)
  }

  invisible(x)
}
