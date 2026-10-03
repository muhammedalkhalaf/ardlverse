#' @title Bootstrap ARDL Bounds Test
#' @description Bounds test for cointegration in an ARDL model with bootstrap
#'   critical values from a recursive bootstrap under the null hypothesis.
#'
#' @details
#' The model is the conditional error correction form of the ARDL(p, q)
#' model of Pesaran, Shin and Smith (2001):
#' \deqn{\Delta y_t = c + \rho y_{t-1} + \theta' x_{t-1} +
#'   \sum_{i=1}^{p-1}\phi_i \Delta y_{t-i} + \sum_{j=0}^{q}\pi_j' \Delta x_{t-j} + u_t.}
#' Three statistics are computed:
#' \itemize{
#'   \item \strong{Fov} (\code{F_stat}): F test that \eqn{\rho} and
#'     \eqn{\theta} are jointly zero (with the intercept in case 2 and the
#'     trend in case 4). \code{F_overall} is the same statistic and is kept
#'     for backward compatibility.
#'   \item \strong{t} (\code{t_stat}): t test that \eqn{\rho = 0}.
#'   \item \strong{Find} (\code{Find_stat}): F test that \eqn{\theta = 0}
#'     (McNown, Sam and Goh, 2018).
#' }
#'
#' The bootstrap is recursive: y* and x* are generated from the estimated
#' model, observation by observation, with resampled residual pairs, and the
#' three statistics are recomputed with the same design on each bootstrap
#' sample. With \code{nulls = "separate"} (default) the procedure follows
#' Bertelli, Vacca and Zoia (2022, Section 3): each statistic has its own
#' restricted model (eqs. 16-18), x* is generated from the marginal model of
#' \eqn{\Delta x_t} on \eqn{x_{t-1}} and lagged differences (eqs. 19-20),
#' the residuals are recentred after each draw (eqs. 21-22) and the initial
#' values are a random block of the data. With \code{nulls = "joint"} the
#' null of the Fov test generates the samples for all statistics (McNown,
#' Sam and Goh, 2018, Steps 1-8) applied to the conditional ECM above, which
#' is a package choice (MSG write the y equation without the contemporaneous
#' \eqn{\Delta x_t}); \eqn{\Delta x_t} follows the unrestricted equation
#' that also contains \eqn{y_{t-1}}, the residuals are recentred once (MSG
#' eq. 13 subtracts the residual mean; no rescaling) and the initial values
#' are the first observations. The x equation has \code{p - 1} lags of
#' \eqn{\Delta y} and \eqn{\Delta x}, as in BVZ eq. 19.
#' Bertelli, Vacca and Zoia (2022) treat cases 2 and 3; cases 1, 4 and 5
#' are handled by the same algorithm with the deterministic terms of the
#' case (the intercept or trend is restricted under the Fov null in cases 2
#' and 4), which is a package choice.
#'
#' Critical values are order statistics of the bootstrap distributions (BVZ
#' eqs. 24-25). Cointegration is concluded only when Fov, t and Find all
#' reject; Fov rejecting without t is reported as the degenerate case of the
#' first type and Fov and t rejecting without Find as the degenerate case of
#' the second type (Pesaran, Shin and Smith, 2001; McNown, Sam and Goh,
#' 2018).
#'
#' @param formula A formula such as \code{gdp ~ investment + trade}; only
#'   plain variable names are accepted.
#' @param data A data frame containing the time series (no missing values).
#' @param p Integer. ARDL order of the dependent variable; the model has
#'   \code{p - 1} lagged differences of y (default: 1)
#' @param q Integer or vector. Lags of the differenced regressors: lags
#'   0 to \code{q} enter the model (default: 1)
#' @param case Integer from 1-5 specifying deterministic components:
#'   \itemize{
#'     \item 1: No intercept, no trend
#'     \item 2: Restricted intercept, no trend
#'     \item 3: Unrestricted intercept, no trend (default)
#'     \item 4: Unrestricted intercept, restricted trend
#'     \item 5: Unrestricted intercept, unrestricted trend
#'   }
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param seed Random seed (default: NULL). The seed is used locally and the
#'   random number generator state of the session is restored afterwards.
#' @param parallel Deprecated and ignored (the bootstrap runs serially).
#' @param ncores Deprecated and ignored.
#' @param nulls \code{"separate"} (default; Bertelli, Vacca and Zoia, 2022) or
#'   \code{"joint"} (Fov null for all statistics, McNown, Sam and Goh,
#'   2018, applied to the conditional ECM); see Details.
#' @param xmodel Model generating \eqn{\Delta x^*}: \code{"vecm"} (default
#'   with separate nulls), \code{"var"} (default with the joint null) or
#'   \code{"rw"} (random walk with the deterministic terms).
#' @param level Significance level of the decision (default: 0.05).
#'
#' @return An object of class "boot_ardl" containing:
#' \itemize{
#'   \item \code{F_stat}, \code{t_stat}, \code{Find_stat}: the statistics;
#'     \code{F_overall} equals \code{F_stat}
#'   \item \code{boot_F}, \code{boot_t}, \code{boot_Find}: bootstrap
#'     distributions (failed replications are \code{NA})
#'   \item \code{cv_F}, \code{cv_t}, \code{cv_Find}: bootstrap critical values
#'     at the 10\%, 5\%, 2.5\% and 1\% levels
#'   \item \code{p_value_F}, \code{p_value_t}, \code{p_value_Find}: bootstrap
#'     p-values
#'   \item \code{decision}: one of \code{"COINTEGRATION"},
#'     \code{"NO_COINTEGRATION"}, \code{"DEGENERATE_1"},
#'     \code{"DEGENERATE_2"}, \code{"UNDETERMINED"}
#'   \item \code{conclusion}: a text version of the decision
#'   \item \code{model}: the estimated ARDL model
#'   \item \code{engine}: the full output of the bootstrap engine (including
#'     the number of valid replications and the reproduction check)
#' }
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' @examples
#' data(macro_data)
#' boot_test <- boot_ardl(gdp ~ inflation, data = macro_data[1:60, ],
#'                        p = 2, q = 1, nboot = 49, seed = 1)
#' boot_test
#' \donttest{
#' summary(boot_test)
#' plot(boot_test)
#' }
#'
#' @export
#' @importFrom stats lm coef residuals fitted var rnorm quantile
boot_ardl <- function(formula, data, p = 1, q = 1, case = 3,
                      nboot = 2000, seed = NULL,
                      parallel = FALSE, ncores = 2,
                      nulls = c("separate", "joint"), xmodel = NULL,
                      level = 0.05) {

  nulls <- match.arg(nulls)
  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  if (isTRUE(parallel))
    warning("'parallel' is deprecated and ignored; the bootstrap runs serially",
            call. = FALSE)
  if (is.null(xmodel)) xmodel <- if (nulls == "separate") "vecm" else "var"
  xmodel <- match.arg(xmodel, c("vecm", "var", "rw"))
  if (p < 1) stop("'p' must be at least 1")

  fv <- .ardl_formula_vars(formula, data)
  y_var <- fv$y_var
  x_vars <- fv$x_vars
  k <- length(x_vars)
  if (length(q) == 1) q <- rep(q, k)
  if (length(q) != k) stop("'q' must have length 1 or one entry per regressor")
  if (any(q < 0)) stop("'q' must be non-negative")

  ardl_data <- .prepare_ardl_ts_data(data, y_var, x_vars, p, q, case)
  if (is.null(ardl_data)) {
    stop("Insufficient observations for specified lag structure")
  }
  model_ur <- .estimate_ardl_unrestricted(ardl_data, case)
  model_r <- .estimate_ardl_restricted(ardl_data, case)
  F_stat <- .compute_F_stat(model_ur, model_r, k, case)
  t_stat <- .compute_t_stat(model_ur)
  F_overall <- F_stat

  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  t0 <- max(3, p + 1, max(q) + 2)
  stat_fun <- function(yy, xx, spec) .boot_ardl_stats(yy, xx, x_vars, p, q, case)
  eng <- .ardl_boot_engine(y, X, p = p - 1, q = q, case = case, t0 = t0,
                           nulls = nulls, xmodel = xmodel, B = nboot,
                           init = if (nulls == "separate") "block" else "observed",
                           recentre = if (nulls == "separate") "draw" else "once",
                           stat_fun = stat_fun, use_find = TRUE, level = level,
                           seed = seed)
  Find_stat <- unname(eng$statistic["Find"])
  cvn <- c("90%", "95%", "97.5%", "99%")
  cv_F <- stats::setNames(eng$cv[, "Fov"], cvn)
  cv_t <- stats::setNames(eng$cv[, "t"], c("10%", "5%", "2.5%", "1%"))
  cv_Find <- stats::setNames(eng$cv[, "Find"], cvn)

  result <- list(
    F_stat = F_stat,
    t_stat = t_stat,
    Find_stat = Find_stat,
    F_overall = F_overall,
    boot_F = eng$boot[, "Fov"],
    boot_t = eng$boot[, "t"],
    boot_Find = eng$boot[, "Find"],
    cv_F = cv_F,
    cv_t = cv_t,
    cv_Find = cv_Find,
    p_value_F = unname(eng$p_value["Fov"]),
    p_value_t = unname(eng$p_value["t"]),
    p_value_Find = unname(eng$p_value["Find"]),
    decision = eng$decision,
    conclusion = eng$label,
    level = level,
    model = model_ur,
    case = case,
    k = k,
    nboot = nboot,
    nulls = nulls,
    xmodel = xmodel,
    engine = eng,
    call = match.call(),
    formula = formula,
    y_var = y_var,
    x_vars = x_vars,
    p = p,
    q = q
  )

  class(result) <- c("boot_ardl", "list")
  return(result)
}


#' @title Statistics of boot_ardl on a (bootstrap) sample
#' @description Fov, t and Find with the design of boot_ardl().
#' @keywords internal
.boot_ardl_stats <- function(y, X, x_vars, p, q, case) {
  d <- data.frame(y, X)
  names(d) <- c(".y_dep", x_vars)
  ad <- .prepare_ardl_ts_data(d, ".y_dep", x_vars, p, q, case)
  Z <- as.matrix(ad[, -1, drop = FALSE])
  dy <- ad$dy
  lv <- attr(ad, "level_vars")
  f <- stats::lm.fit(Z, dy)
  if (f$rank < ncol(Z)) stop("rank-deficient design")
  rss <- sum(f$residuals^2)
  df <- length(dy) - ncol(Z)
  Ftest <- function(drop) {
    keep <- setdiff(colnames(Z), drop)
    rss_r <- if (length(keep)) sum(stats::lm.fit(Z[, keep, drop = FALSE], dy)$residuals^2) else sum(dy^2)
    ((rss_r - rss) / length(drop)) / (rss / df)
  }
  piv <- order(f$qr$pivot)
  XtXi <- chol2inv(f$qr$qr[seq_len(ncol(Z)), seq_len(ncol(Z)), drop = FALSE])[piv, piv]
  iy <- which(colnames(Z) == "y_lag1")
  c(Fov = Ftest(c(lv, if (case == 2) "const", if (case == 4) "trend")),
    t = unname(f$coefficients[iy] / sqrt(rss / df * XtXi[iy, iy])),
    Find = Ftest(setdiff(lv, "y_lag1")))
}


#' @title Prepare Time Series Data for ARDL
#' @keywords internal
.prepare_ardl_ts_data <- function(data, y_var, x_vars, p, q, case) {
  
  n <- nrow(data)
  max_lag <- max(p, max(q))
  
  if (n <= max_lag + 5) {
    return(NULL)
  }
  
  # Dependent variable
  y <- data[[y_var]]
  dy <- diff(y)
  y_lag1 <- y[-length(y)]
  
  # Lagged dy
  dy_lags <- NULL
  if (p > 1) {
    dy_lags <- sapply(1:(p-1), function(lag) {
      c(rep(NA, lag), dy[1:(length(dy) - lag)])
    })
    colnames(dy_lags) <- paste0("dy_L", 1:(p-1))
  }
  
  # X variables: levels and differences
  X_levels <- as.matrix(data[-nrow(data), x_vars, drop = FALSE])
  
  X_diff <- sapply(x_vars, function(v) diff(data[[v]]))
  if (is.vector(X_diff)) X_diff <- matrix(X_diff, ncol = 1)
  colnames(X_diff) <- paste0("d", x_vars)
  
  # Lagged X differences
  X_diff_lags <- NULL
  for (j in seq_along(x_vars)) {
    if (q[j] > 0) {
      for (lag in 1:q[j]) {
        dx <- diff(data[[x_vars[j]]])
        lagged <- c(rep(NA, lag), dx[1:(length(dx) - lag)])
        X_diff_lags <- cbind(X_diff_lags, lagged)
        colnames(X_diff_lags)[ncol(X_diff_lags)] <- paste0("d", x_vars[j], "_L", lag)
      }
    }
  }
  
  # Combine
  result <- data.frame(
    dy = dy[-1],
    y_lag1 = y_lag1[-1],
    X_levels[-1, , drop = FALSE],
    X_diff[-1, , drop = FALSE]
  )
  
  if (!is.null(dy_lags)) {
    result <- cbind(result, dy_lags[-1, , drop = FALSE])
  }
  
  if (!is.null(X_diff_lags)) {
    result <- cbind(result, X_diff_lags[-1, , drop = FALSE])
  }
  
  # Add deterministics based on case
  n_obs <- nrow(result)
  # Case 2 has a restricted intercept, case 3 an unrestricted one; both
  # need the constant in the model (it is restricted under H0 in case 2)
  if (case >= 2) {
    result$const <- 1
  }
  if (case >= 4) {
    result$trend <- 1:n_obs
  }
  
  # Remove NAs
  result <- na.omit(result)
  attr(result, "level_vars") <- c("y_lag1", names(result)[2 + seq_along(x_vars)])
  
  return(result)
}


#' @title Estimate Unrestricted ARDL Model
#' @keywords internal
.estimate_ardl_unrestricted <- function(ardl_data, case) {
  
  # All variables except dy
  xvars <- names(ardl_data)[-1]
  
  # The constant is an explicit column (cases 2-5), so lm() adds none
  formula_str <- paste("dy ~", paste(xvars, collapse = " + "), "- 1")
  
  model <- lm(as.formula(formula_str), data = ardl_data)
  
  return(model)
}


#' @title Estimate Restricted ARDL Model
#' @keywords internal
.estimate_ardl_restricted <- function(ardl_data, case) {

  # Under H0 the lagged levels are zero, together with the intercept
  # (case 2) or the trend (case 4)
  xvars <- names(ardl_data)[-1]
  level_vars <- attr(ardl_data, "level_vars")
  drop <- c(level_vars, if (case == 2) "const", if (case == 4) "trend")
  keep_vars <- setdiff(xvars, drop)
  if (length(keep_vars) == 0) {
    return(lm(dy ~ 0, data = ardl_data))
  }

  # The constant is an explicit column, so lm() adds none
  formula_str <- paste("dy ~", paste(keep_vars, collapse = " + "), "- 1")
  model <- lm(as.formula(formula_str), data = ardl_data)
  
  return(model)
}


#' @title Compute F-statistic for Bounds Test
#' @keywords internal
.compute_F_stat <- function(model_ur, model_r, k, case) {
  
  RSS_ur <- sum(residuals(model_ur)^2)
  RSS_r <- sum(residuals(model_r)^2)
  
  n <- length(residuals(model_ur))
  k_ur <- length(coef(model_ur))
  
  # Number of restrictions: the k + 1 lagged levels, plus the intercept
  # (case 2) or the trend (case 4)
  m <- k + 1 + (case %in% c(2, 4))
  
  F_stat <- ((RSS_r - RSS_ur) / m) / (RSS_ur / (n - k_ur))
  
  return(F_stat)
}


#' @title Compute t-statistic for EC Coefficient
#' @keywords internal
.compute_t_stat <- function(model_ur) {
  
  coefs <- summary(model_ur)$coefficients
  
  # Find y_lag1 coefficient
  y_lag1_idx <- which(rownames(coefs) == "y_lag1")
  
  if (length(y_lag1_idx) == 0) {
    return(NA)
  }
  
  t_stat <- coefs[y_lag1_idx, "t value"]
  
  return(t_stat)
}


#' @title Compute Overall F-statistic
#' @keywords internal
.compute_F_overall <- function(model_ur, model_r, k, case) {
  # Same as F_stat for standard bounds test
  return(.compute_F_stat(model_ur, model_r, k, case))
}


#' @title Critical Value Bounds for the PSS Bounds Test
#' @description Critical value bounds for the F and t statistics of the
#'   Pesaran, Shin and Smith (2001) bounds test, computed from the response
#'   surface regressions of Kripfganz and Schneider (2020). With \code{n}
#'   missing the asymptotic bounds are returned; with \code{n} (and
#'   \code{sr}) the finite-sample bounds.
#'
#' @param k Number of regressors in levels (excluding the lagged dependent
#'   variable)
#' @param case PSS case (1-5)
#' @param level Significance level: \code{"10\%"}, \code{"5\%"},
#'   \code{"2.5\%"} or \code{"1\%"}
#' @param n Number of observations (optional)
#' @param sr Number of short-run coefficients (regressors other than the
#'   deterministic terms and the lagged levels); used with \code{n}
#'
#' @return A list with \code{F_bounds} and \code{t_bounds} (each a list with
#'   \code{I0} and \code{I1}), \code{k}, \code{case} and \code{level}. The
#'   t statistic is not tabulated for cases 2 and 4; the values of cases 3
#'   and 5 are used, as in Kripfganz and Schneider (2020).
#'
#' @references
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' @examples
#' pss_critical_values(k = 2, case = 3)
#' pss_critical_values(k = 2, case = 3, n = 80, sr = 4)
#'
#' @export
pss_critical_values <- function(k, case = 3, level = "5%", n = NULL, sr = 0) {
  lev <- as.numeric(sub("%", "", level))
  Fb <- .ks_bounds("F", case, k, n, sr, siglevels = lev)
  tb <- .ks_bounds("t", case, k, n, sr, siglevels = lev)
  list(
    F_bounds = list(I0 = unname(Fb$cv["I0", 1]), I1 = unname(Fb$cv["I1", 1])),
    t_bounds = list(I0 = unname(tb$cv["I0", 1]), I1 = unname(tb$cv["I1", 1])),
    k = k, case = case, level = level
  )
}


#' @title PSS bounds test p-values (Kripfganz and Schneider, 2020)
#' @keywords internal
.pss_pvalues <- function(F_stat, t_stat, k, case, n = NULL, sr = 0) {
  list(F = .ks_bounds("F", case, k, n, sr, value = F_stat)$pvalue,
       t = .ks_bounds("t", case, k, n, sr, value = t_stat)$pvalue)
}


#' @title Summary method for boot_ardl
#' @param object An object of class "boot_ardl"
#' @param ... Not used
#' @export
summary.boot_ardl <- function(object, ...) {

  cat("\n")
  cat("====================================================================\n")
  cat("     Bootstrap ARDL Bounds Test for Cointegration\n")
  cat("====================================================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Model: ARDL(", object$p, ", ", paste(object$q, collapse = ", "), ")\n", sep = "")
  cat("Case:  ", object$case, " (",
      switch(object$case,
             "1" = "No intercept, no trend",
             "2" = "Restricted intercept, no trend",
             "3" = "Unrestricted intercept, no trend",
             "4" = "Unrestricted intercept, restricted trend",
             "5" = "Unrestricted intercept, unrestricted trend"),
      ")\n", sep = "")
  cat("Regressors (k):", object$k, "\n")
  nulls <- if (is.null(object$nulls)) "separate" else object$nulls
  cat("Bootstrap:", if (nulls == "separate")
    "separate nulls (Bertelli, Vacca and Zoia, 2022)" else
      "Fov null for all statistics (McNown, Sam and Goh, 2018, Steps 1-8)\n  applied to the conditional ECM (package choice)", "\n")
  cat("Replications:", object$nboot, " valid (Fov, t, Find):",
      paste(object$engine$n_valid, collapse = ", "), "\n\n")

  cat("--------------------------------------------------------------------\n")
  cat("            Statistic   p-value      10%       5%     2.5%       1%\n")
  cat("--------------------------------------------------------------------\n")
  row <- function(nm, st, pv, cv)
    cat(sprintf("%-10s %10.4f %9.4f %8.3f %8.3f %8.3f %8.3f\n",
                nm, st, pv, cv[1], cv[2], cv[3], cv[4]))
  row("Fov", object$F_stat, object$p_value_F, object$cv_F)
  row("t", object$t_stat, object$p_value_t, object$cv_t)
  if (!is.null(object$Find_stat))
    row("Find", object$Find_stat, object$p_value_Find, object$cv_Find)
  cat("Critical values are bootstrap order statistics; reject for Fov and\n")
  cat("Find above, and for t below, the critical value.\n\n")

  pss <- pss_critical_values(object$k, object$case, "5%")
  cat("--------------------------------------------------------------------\n")
  cat("  For reference: PSS asymptotic bounds, Kripfganz and Schneider\n")
  cat("  (2020), 5% level\n")
  cat("--------------------------------------------------------------------\n")
  cat("F-bounds: I(0) =", round(pss$F_bounds$I0, 2),
      ", I(1) =", round(pss$F_bounds$I1, 2), "\n")
  cat("t-bounds: I(0) =", round(pss$t_bounds$I0, 2),
      ", I(1) =", round(pss$t_bounds$I1, 2), "\n\n")

  cat("--------------------------------------------------------------------\n")
  cat("  Decision at the ", 100 * object$level,
      "% level (Fov, t and Find must all reject)\n", sep = "")
  cat("--------------------------------------------------------------------\n")
  cat(object$conclusion, "\n")
  cat("====================================================================\n")

  invisible(object)
}


#' @title Print method for boot_ardl
#' @param x An object of class "boot_ardl"
#' @param ... Not used
#' @export
print.boot_ardl <- function(x, ...) {

  cat("\nBootstrap ARDL Bounds Test\n")
  cat("Fov :", round(x$F_stat, 4), "(bootstrap p-value:", round(x$p_value_F, 4), ")\n")
  cat("t   :", round(x$t_stat, 4), "(bootstrap p-value:", round(x$p_value_t, 4), ")\n")
  if (!is.null(x$Find_stat))
    cat("Find:", round(x$Find_stat, 4), "(bootstrap p-value:",
        round(x$p_value_Find, 4), ")\n")
  cat("\n", x$conclusion, "\n")

  invisible(x)
}


#' @title Plot Bootstrap Distribution
#' @description Plot the bootstrap distribution of test statistics
#'
#' @param x An object of class "boot_ardl"
#' @param which Character. "F", "t", or "both" (default)
#' @param ... Additional arguments passed to plotting functions
#'
#' @export
plot.boot_ardl <- function(x, which = "both", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  
  plots <- list()
  
  if (which %in% c("F", "both")) {
    df_F <- data.frame(F_stat = x$boot_F[!is.na(x$boot_F)])
    
    p1 <- ggplot2::ggplot(df_F, ggplot2::aes(x = F_stat)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = 50, fill = "steelblue", alpha = 0.7) +
      ggplot2::geom_density(color = "darkblue", linewidth = 1) +
      ggplot2::geom_vline(xintercept = x$F_stat, color = "red", 
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = x$cv_F["95%"], color = "orange",
                          linewidth = 1, linetype = "dotted") +
      ggplot2::labs(
        title = "Bootstrap Distribution of F-statistic",
        subtitle = paste("Observed F =", round(x$F_stat, 3),
                        "| 5% CV =", round(x$cv_F["95%"], 3)),
        x = "F-statistic",
        y = "Density"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    plots$F <- p1
  }
  
  if (which %in% c("t", "both")) {
    df_t <- data.frame(t_stat = x$boot_t[!is.na(x$boot_t)])
    
    p2 <- ggplot2::ggplot(df_t, ggplot2::aes(x = t_stat)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = 50, fill = "darkgreen", alpha = 0.7) +
      ggplot2::geom_density(color = "forestgreen", linewidth = 1) +
      ggplot2::geom_vline(xintercept = x$t_stat, color = "red",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = x$cv_t["5%"], color = "orange",
                          linewidth = 1, linetype = "dotted") +
      ggplot2::labs(
        title = "Bootstrap Distribution of t-statistic",
        subtitle = paste("Observed t =", round(x$t_stat, 3),
                        "| 5% CV =", round(x$cv_t["5%"], 3)),
        x = "t-statistic",
        y = "Density"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    plots$t <- p2
  }
  
  if (which == "both" && requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(plots$F, plots$t, ncol = 2)
  } else if (which == "F") {
    print(plots$F)
  } else if (which == "t") {
    print(plots$t)
  } else if (which == "both") {
    print(plots$F)
    print(plots$t)
  }
  
  invisible(plots)
}
