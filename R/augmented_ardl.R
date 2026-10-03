#' @title Augmented ARDL Bounds Test (AARDL)
#' @description Augmented ARDL bounds test with the overall F test, the t
#'   test on the lagged dependent variable and the F test on the lagged
#'   independent variables, in linear, nonlinear (NARDL) and Fourier
#'   variants, with Kripfganz and Schneider (2020) bounds or bootstrap
#'   critical values.
#'
#' @details
#' Three statistics are computed from the conditional error correction model:
#' \itemize{
#'   \item \code{F_pss}: F test that the lagged levels of y and of the
#'     regressors are jointly zero (with the intercept in case 2 and the
#'     trend in case 4; Pesaran, Shin and Smith, 2001)
#'   \item \code{t_dep}: t test on the lagged dependent variable
#'   \item \code{F_ind}: F test on the lagged independent variables
#'     (McNown, Sam and Goh, 2018)
#' }
#'
#' The function supports 8 sub-models:
#' \enumerate{
#'   \item \code{"linear"}: ARDL, Kripfganz and Schneider (2020) bounds for
#'     F_pss and t_dep (F_ind has no tabulated bounds)
#'   \item \code{"nardl"}: NARDL with positive and negative partial sums;
#'     the bounds use k = number of partial sums
#'   \item \code{"fourier"}: ARDL with Fourier terms. No valid bounds exist
#'     for the bounds test with Fourier terms; no decision is reported and
#'     the bounds shown are those of the model without Fourier terms, for
#'     reference only (as in \code{\link{fourier_bounds_test}})
#'   \item \code{"fnardl"}: NARDL with Fourier terms; as \code{"fourier"}
#'   \item \code{"bootstrap"}: ARDL with bootstrap critical values
#'   \item \code{"bnardl"}: NARDL with bootstrap critical values
#'   \item \code{"fbootstrap"}: Fourier ARDL with bootstrap critical values
#'   \item \code{"fbnardl"}: Fourier NARDL with bootstrap critical values
#' }
#' The bootstrap types use the recursive bootstrap of
#' \code{\link{boot_ardl}} (separate nulls of Bertelli, Vacca and Zoia,
#' 2022, by default, or the Fov null for all statistics of McNown, Sam and
#' Goh, 2018, applied to the conditional ECM, a package choice): y* and
#' x* are generated recursively, the partial sums are rebuilt from x* with
#' the same decomposition as the data, and the Fourier terms at the chosen
#' \code{fourier_k} are held fixed. The decision requires F_pss, t_dep and
#' F_ind to reject; the degenerate cases are labelled. With Fourier terms the
#' bootstrap is valid for the given \code{fourier_k}; it does not account for
#' a data-based choice of \code{fourier_k}.
#' With Fourier terms and partial sums (type \code{"fbnardl"}) the F_ind test
#' of the separate-null bootstrap over-rejects: in a Monte Carlo with
#' independent random walks (n = 100, 200 samples, B = 199) it rejected in
#' 16.0\% (Monte Carlo s.e. 2.6\%) of the samples at the 5\% level, and an
#' independent review Monte Carlo gave 11.3\% (s.e. 1.8\%); the combined
#' decision rejected in 6.0\% (s.e. 1.7\%).

#' The name follows the augmented ARDL bounds test of Sam, McNown and Goh
#' (2019); that paper was not available for this release, so its bootstrap
#' details are not reproduced here (pending verification): the bootstrap is
#' the one of Bertelli, Vacca and Zoia (2022) or McNown, Sam and Goh (2018).
#'
#' @param formula A formula such as \code{gdp ~ investment + trade}; only
#'   plain variable names are accepted.
#' @param data A data frame containing the time series (no missing values).
#' @param p Integer. Number of lagged differences of the dependent variable
#'   (default: 1)
#' @param q Integer or vector. Number of differences of each regressor,
#'   lags 0 to \code{q - 1} (default: 1)
#' @param case Integer from 1-5 specifying deterministic components
#' @param type Character. Model type: "linear", "nardl", "fourier", "fnardl",
#'   "bootstrap", "bnardl", "fbootstrap", "fbnardl" (default: "linear")
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param fourier_k Integer. Number of Fourier frequencies (default: 1, max: 3)
#' @param threshold Numeric. Threshold value for NARDL decomposition (default: 0)
#' @param seed Random seed; used locally, the session's random number
#'   generator state is restored afterwards.
#' @param nulls Bootstrap nulls, \code{"separate"} (default) or
#'   \code{"joint"}; see \code{\link{boot_ardl}}.
#'
#' @return An object of class "aardl" containing:
#' \itemize{
#'   \item \code{F_pss}: PSS F-statistic for bounds test
#'   \item \code{t_dep}: t-statistic for the lagged dependent variable
#'   \item \code{F_ind}: F-statistic for the lagged independent variables
#'   \item \code{critical_values}: Kripfganz and Schneider (2020) bounds
#'     (for the Fourier types: bounds of the model without Fourier terms,
#'     not valid for the estimated model)
#'   \item \code{boot_results}: bootstrap distributions, critical values and
#'     the engine output (bootstrap types)
#'   \item \code{conclusion}: list with \code{decision}, \code{message},
#'     \code{method} (and \code{p_values} for the bootstrap types)
#'   \item \code{model}: The estimated ARDL model
#'   \item \code{long_run}: Long-run coefficients
#'   \item \code{short_run}: Short-run coefficients
#'   \item \code{diagnostics}: Model diagnostic tests
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
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Sam, C. Y., McNown, R. and Goh, S. K. (2019). An augmented autoregressive
#' distributed lag bounds test for cointegration. \emph{Economic Modelling},
#' 80, 130-141. \doi{10.1016/j.econmod.2018.11.001}
#'
#' @examples
#' data <- generate_ts_data(n = 80, seed = 1)
#' result <- aardl(gdp ~ investment + trade, data = data, p = 2, q = 2, case = 3)
#' result
#' boot <- aardl(gdp ~ investment, data = data, type = "bnardl", nboot = 29,
#'               seed = 1)
#' boot
#' \donttest{
#' summary(result)
#' result_fourier <- aardl(gdp ~ investment, data = data, type = "fbootstrap",
#'                         fourier_k = 1, nboot = 49, seed = 1)
#' summary(result_fourier)
#' }
#'
#' @export
aardl <- function(formula, data, p = 1, q = 1, case = 3,
                  type = c("linear", "nardl", "fourier", "fnardl",
                          "bootstrap", "bnardl", "fbootstrap", "fbnardl"),
                  nboot = 2000, fourier_k = 1, threshold = 0, seed = NULL,
                  nulls = c("separate", "joint")) {

  type <- match.arg(type)
  nulls <- match.arg(nulls)

  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  if (fourier_k < 1 || fourier_k > 3) {
    stop("'fourier_k' must be between 1 and 3")
  }
  if (p < 1) stop("'p' must be at least 1")

  fv <- .ardl_formula_vars(formula, data)
  y_var <- fv$y_var
  x_vars <- fv$x_vars
  k0 <- length(x_vars)
  if (length(q) == 1) q <- rep(q, k0)
  if (length(q) != k0) stop("'q' must have length 1 or one entry per regressor")
  if (any(q < 1)) stop("'q' must be at least 1 (lags 0 to q - 1 of the differences)")

  y <- data[[y_var]]
  X0 <- as.matrix(data[, x_vars, drop = FALSE])

  use_fourier <- type %in% c("fourier", "fnardl", "fbootstrap", "fbnardl")
  use_nardl <- type %in% c("nardl", "fnardl", "bnardl", "fbnardl")
  use_bootstrap <- type %in% c("bootstrap", "bnardl", "fbootstrap", "fbnardl")

  spec <- list(y_var = y_var, x_vars = x_vars, p = p, q = q, case = case,
               use_nardl = use_nardl, use_fourier = use_fourier,
               fourier_k = fourier_k, threshold = threshold)
  dd <- .aardl_design(y, X0, spec)
  design <- dd$design
  dy <- dd$dy
  k <- dd$k
  n_valid <- length(dy)
  q_used <- dd$q_levels

  model <- stats::lm(dy ~ design - 1)
  coefs <- stats::coef(model)
  ec_coef <- coefs[1]

  st <- .aardl_stats(dy, design, k, case)
  F_pss <- unname(st["Fov"])
  t_dep <- unname(st["t"])
  F_ind <- if (k > 0) unname(st["Find"]) else NA

  # Kripfganz and Schneider (2020) bounds; for the Fourier types these are
  # the bounds of the model without Fourier terms (reference only)
  n_level_vars <- 1 + k
  n_fourier <- if (use_fourier) 2 * fourier_k else 0
  sr <- ncol(design) - n_level_vars - (case >= 2) - (case >= 4) - n_fourier
  cv <- .aardl_critical_values(k, case, n_valid, sr)

  boot_results <- NULL
  if (use_bootstrap) {
    stat_fun <- function(yy, xx, sp) .aardl_stats_from_data(yy, xx, spec)
    dec <- if (use_nardl) .aardl_decomposer(threshold) else NULL
    eng <- .ardl_boot_engine(y, X0, p = p, q = q_used - 1, case = case,
                             t0 = max(p + 1, max(q_used)) + 1,
                             det = dd$fourier_terms, decompose = dec,
                             nulls = nulls,
                             xmodel = if (nulls == "separate") "vecm" else "var",
                             B = nboot,
                             init = if (nulls == "separate") "block" else "observed",
                             recentre = if (nulls == "separate") "draw" else "once",
                             stat_fun = stat_fun, use_find = k > 0, seed = seed)
    nm <- c("90%", "95%", "97.5%", "99%")
    boot_results <- list(
      F_dist = eng$boot[, "Fov"],
      t_dist = eng$boot[, "t"],
      F_ind_dist = eng$boot[, "Find"],
      cv_F = stats::setNames(eng$cv[, "Fov"], nm),
      cv_t = stats::setNames(eng$cv[, "t"], c("10%", "5%", "2.5%", "1%")),
      cv_F_ind = stats::setNames(eng$cv[, "Find"], nm),
      engine = eng
    )
  }

  conclusion <- .aardl_conclusion(F_pss, t_dep, F_ind, cv, boot_results,
                                  fourier = use_fourier)

  if (abs(ec_coef) > 1e-10) {
    lr_coefs <- -coefs[2:(1 + k)] / ec_coef
    names(lr_coefs) <- dd$level_names
  } else {
    lr_coefs <- rep(NA, k)
  }

  sr_start <- 1 + k + 1
  sr_end <- sr_start + p - 1
  sr_coefs <- coefs[sr_start:sr_end]

  fit_stats <- list(
    R2 = summary(model)$r.squared,
    adj_R2 = summary(model)$adj.r.squared,
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    sigma = summary(model)$sigma,
    df = model$df.residual
  )

  resid <- stats::residuals(model)
  diagnostics <- list(
    serial_corr = .breusch_godfrey_test(resid, 2),
    heteroskedasticity = .breusch_pagan_test(model),
    normality = stats::shapiro.test(resid)$p.value
  )

  result <- list(
    F_pss = F_pss,
    t_dep = t_dep,
    F_ind = F_ind,
    critical_values = cv,
    boot_results = boot_results,
    conclusion = conclusion,
    model = model,
    coefficients = coefs,
    long_run = lr_coefs,
    short_run = sr_coefs,
    fit = fit_stats,
    diagnostics = diagnostics,
    call = match.call(),
    type = type,
    case = case,
    n = n_valid,
    k = k,
    p = p,
    q = q_used,
    fourier_k = if (use_fourier) fourier_k else NULL
  )

  class(result) <- "aardl"
  return(result)
}


#' @title Design matrix of aardl()
#' @keywords internal
.aardl_design <- function(y, X0, spec) {
  y_var <- spec$y_var
  x_vars <- spec$x_vars
  p <- spec$p
  q <- spec$q
  case <- spec$case
  n <- length(y)
  X <- X0
  level_names <- x_vars
  if (spec$use_nardl) {
    X <- .decompose_asymmetric(X0, spec$threshold)
    level_names <- c(paste0(x_vars, "_pos"), paste0(x_vars, "_neg"))
    q <- rep(q, 2)
  }
  k <- ncol(X)
  fourier_terms <- if (spec$use_fourier) .create_fourier_terms(n, spec$fourier_k) else NULL

  # p lagged differences of y need p + 1 initial observations
  max_lag <- max(p + 1, max(q))
  valid_idx <- (max_lag + 1):n
  n_valid <- length(valid_idx)
  if (n_valid <= 5) stop("Insufficient observations for specified lag structure")

  dy <- diff(y)[(max_lag):(n - 1)]
  y_lag <- y[valid_idx - 1]
  dy_lags <- matrix(NA, n_valid, p)
  for (i in 1:p) {
    dy_lags[, i] <- diff(y)[(max_lag - i):(n - 1 - i)]
  }
  colnames(dy_lags) <- paste0("d.", y_var, ".l", 1:p)

  x_levels <- matrix(NA, n_valid, k)
  x_diff_list <- list()
  for (j in 1:k) {
    x_j <- X[, j]
    x_levels[, j] <- x_j[valid_idx - 1]
    dx_j <- diff(x_j)
    x_diff_j <- matrix(NA, n_valid, q[j])
    for (i in 0:(q[j] - 1)) {
      x_diff_j[, i + 1] <- dx_j[(max_lag - i):(n - 1 - i)]
    }
    x_diff_list[[j]] <- x_diff_j
  }
  colnames(x_levels) <- level_names
  x_diffs <- do.call(cbind, x_diff_list)
  design <- cbind(y_lag, x_levels, dy_lags, x_diffs)
  if (case >= 2) design <- cbind(design, intercept = 1)
  if (case >= 4) design <- cbind(design, trend = 1:n_valid)
  if (spec$use_fourier) design <- cbind(design, fourier_terms[valid_idx, , drop = FALSE])
  list(dy = dy, design = design, k = k, q = spec$q, q_levels = q,
       level_names = level_names, fourier_terms = fourier_terms)
}


#' @title Statistics of aardl(): F_pss (Fov), t_dep (t) and F_ind (Find)
#' @keywords internal
.aardl_stats <- function(dy, design, k, case) {
  f <- stats::lm.fit(design, dy)
  if (f$rank < ncol(design)) stop("rank-deficient design")
  dfr <- length(dy) - ncol(design)
  s2 <- sum(f$residuals^2) / dfr
  nc <- ncol(design)
  piv <- order(f$qr$pivot)
  V <- s2 * chol2inv(f$qr$qr[seq_len(nc), seq_len(nc), drop = FALSE])[piv, piv]
  b <- f$coefficients
  W <- function(i) drop(crossprod(b[i], solve(V[i, i, drop = FALSE], b[i]))) / length(i)
  dn <- colnames(design)
  idx_h0 <- seq_len(1 + k)
  if (case == 2) idx_h0 <- c(idx_h0, which(dn == "intercept"))
  if (case == 4) idx_h0 <- c(idx_h0, which(dn == "trend"))
  c(Fov = W(idx_h0), t = unname(b[1] / sqrt(V[1, 1])),
    Find = if (k > 0) W(1 + seq_len(k)) else NA_real_)
}

.aardl_stats_from_data <- function(y, X0, spec) {
  dd <- .aardl_design(y, X0, spec)
  .aardl_stats(dd$dy, dd$design, dd$k, spec$case)
}

# Decomposition passed to the bootstrap engine: the partial sums of
# .decompose_asymmetric() and their increments
.aardl_decomposer <- function(threshold) {
  list(full = function(x) .decompose_asymmetric(x, threshold),
       step = function(d) {
         a <- d - threshold
         b <- d + threshold
         c(a * (a > 0), b * (b < 0))
       })
}


#' @title Decompose Variables for Asymmetric Analysis
#' @keywords internal
.decompose_asymmetric <- function(X, threshold = 0) {
  X <- as.matrix(X)
  k <- ncol(X)
  n <- nrow(X)
  
  X_pos <- matrix(0, n, k)
  X_neg <- matrix(0, n, k)
  
  for (j in 1:k) {
    dx <- c(0, diff(X[, j]))
    X_pos[, j] <- cumsum(pmax(dx - threshold, 0))
    X_neg[, j] <- cumsum(pmin(dx + threshold, 0))
  }
  
  result <- cbind(X_pos, X_neg)
  return(result)
}


#' @title Create Fourier Terms
#' @keywords internal
.create_fourier_terms <- function(n, k) {
  fourier_mat <- matrix(NA, n, 2 * k)
  t_idx <- 1:n
  
  for (i in 1:k) {
    fourier_mat[, 2*i - 1] <- sin(2 * pi * i * t_idx / n)
    fourier_mat[, 2*i] <- cos(2 * pi * i * t_idx / n)
  }
  
  colnames(fourier_mat) <- c(rbind(paste0("sin", 1:k), paste0("cos", 1:k)))
  return(fourier_mat)
}


#' @title AARDL Critical Values
#' @description Bounds for the overall F and the t test from the response
#'   surfaces of Kripfganz and Schneider (2020). The F test on the lagged
#'   regressors has no tabulated bounds; the bootstrap types give critical
#'   values for it.
#' @keywords internal
.aardl_critical_values <- function(k, case, n, sr = 0) {
  lv <- c(10, 5, 1)
  Fb <- .ks_bounds("F", case, k, n, sr, siglevels = lv)$cv
  tb <- .ks_bounds("t", case, k, n, sr, siglevels = lv)$cv
  nm <- c("90%", "95%", "99%")
  list(
    F = list(I0 = stats::setNames(Fb["I0", ], nm), I1 = stats::setNames(Fb["I1", ], nm)),
    t = list(I0 = stats::setNames(tb["I0", ], nm), I1 = stats::setNames(tb["I1", ], nm))
  )
}


#' @title AARDL Conclusion
#' @description Bootstrap decision (F_pss, t_dep and F_ind must all reject;
#'   degenerate cases labelled), Kripfganz and Schneider (2020) bounds
#'   decision for the models without Fourier terms, and no decision for the
#'   Fourier models without bootstrap.
#' @keywords internal
.aardl_conclusion <- function(F_pss, t_dep, F_ind, cv, boot = NULL,
                              fourier = FALSE) {

  if (!is.null(boot)) {
    eng <- boot$engine
    return(list(
      decision = eng$decision,
      message = paste0(eng$label, " (bootstrap, ", 100 * eng$level, "% level)"),
      p_values = c(F_pss = unname(eng$p_value["Fov"]),
                   t_dep = unname(eng$p_value["t"]),
                   F_ind = unname(eng$p_value["Find"])),
      method = "bootstrap"
    ))
  }

  F_lower <- cv$F$I0["95%"]
  F_upper <- cv$F$I1["95%"]
  t_lower <- cv$t$I0["95%"]
  t_upper <- cv$t$I1["95%"]

  if (fourier) {
    return(list(
      decision = "NOT_AVAILABLE",
      message = paste("No valid bounds exist for the bounds test with Fourier",
                      "terms; no decision is reported. Use type = 'fbootstrap'",
                      "or 'fbnardl'. The bounds shown are for the model without",
                      "Fourier terms; not valid with Fourier terms."),
      bounds = list(F = c(F_lower, F_upper), t = c(t_lower, t_upper)),
      method = "none"
    ))
  }

  if (any(is.na(c(F_lower, F_upper, t_lower, t_upper, F_pss, t_dep)))) {
    conclusion <- "Bounds unavailable (too few degrees of freedom); use the bootstrap"
    decision <- "INCONCLUSIVE"
  } else if (F_pss > F_upper && t_dep < t_upper) {
    conclusion <- "Cointegration: F > I(1) bound and t < I(1) bound"
    decision <- "COINTEGRATION"
  } else if (F_pss < F_lower || t_dep > t_lower) {
    conclusion <- "No cointegration: Statistics within I(0) bounds"
    decision <- "NO_COINTEGRATION"
  } else {
    conclusion <- "Inconclusive: Statistics between I(0) and I(1) bounds"
    decision <- "INCONCLUSIVE"
  }

  list(
    decision = decision,
    message = conclusion,
    bounds = list(F = c(F_lower, F_upper), t = c(t_lower, t_upper)),
    method = "Kripfganz-Schneider bounds"
  )
}


#' @title Breusch-Godfrey Serial Correlation Test
#' @keywords internal
.breusch_godfrey_test <- function(resid, order = 2) {
  n <- length(resid)
  resid_lag <- stats::embed(resid, order + 1)
  y <- resid_lag[, 1]
  X <- resid_lag[, -1, drop = FALSE]
  
  aux_model <- stats::lm(y ~ X)
  r2 <- summary(aux_model)$r.squared
  
  LM <- (n - order) * r2
  p_value <- 1 - stats::pchisq(LM, order)
  
  list(statistic = LM, p.value = p_value, df = order)
}


#' @title Breusch-Pagan Heteroskedasticity Test
#' @keywords internal
.breusch_pagan_test <- function(model) {
  resid <- stats::residuals(model)
  fitted_vals <- stats::fitted(model)
  
  resid_sq <- resid^2
  aux_model <- stats::lm(resid_sq ~ fitted_vals)
  
  n <- length(resid)
  r2 <- summary(aux_model)$r.squared
  LM <- n * r2
  p_value <- 1 - stats::pchisq(LM, 1)
  
  list(statistic = LM, p.value = p_value, df = 1)
}


#' @rdname aardl
#' @param x,object An object of class "aardl"
#' @param ... Not used
#' @export
print.aardl <- function(x, ...) {
  cat("\n")
  cat("Augmented ARDL Bounds Test (", toupper(x$type), ")\n", sep = "")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  cat("Model: Case", x$case, "| p =", x$p, "| k =", x$k, "\n")
  cat("Observations:", x$n, "\n\n")
  
  cat("Test Statistics:\n")
  cat(sprintf("  F_pss (bounds test):    %8.4f\n", x$F_pss))
  cat(sprintf("  t_dep (EC coefficient): %8.4f\n", x$t_dep))
  if (!is.na(x$F_ind)) {
    cat(sprintf("  F_ind (indep. vars):    %8.4f\n", x$F_ind))
  }
  
  cat("\nConclusion:", x$conclusion$decision, "\n")
  cat(x$conclusion$message, "\n")
  
  invisible(x)
}


#' @rdname aardl
#' @export
summary.aardl <- function(object, ...) {
  cat("\n")
  cat("===============================================\n")
  cat("    Augmented ARDL Bounds Test Results\n")
  cat("===============================================\n\n")
  
  cat("Model Specification:\n")
  cat("  Type:", toupper(object$type), "\n")
  cat("  Case:", object$case, "\n")
  cat("  Lags: p =", object$p, ", q =", paste(object$q, collapse = ","), "\n")
  if (!is.null(object$fourier_k)) {
    cat("  Fourier frequencies:", object$fourier_k, "\n")
  }
  cat("  Sample size:", object$n, "\n\n")
  
  cat("Test Statistics:\n")
  cat("-----------------------------------------------\n")
  cat(sprintf("  %-25s %10.4f\n", "F_pss (PSS bounds test):", object$F_pss))
  cat(sprintf("  %-25s %10.4f\n", "t_dep (EC coefficient):", object$t_dep))
  if (!is.na(object$F_ind)) {
    cat(sprintf("  %-25s %10.4f\n", "F_ind (indep. variables):", object$F_ind))
  }
  
  if (object$conclusion$method == "bootstrap") {
    br <- object$boot_results
    eng <- br$engine
    cat("\nBootstrap (", if (eng$settings$nulls == "separate")
      "separate nulls, Bertelli, Vacca and Zoia 2022" else
        "Fov null for all statistics, McNown, Sam and Goh 2018, conditional ECM", "), ", eng$settings$B,
      " replications\n", sep = "")
    cat("-----------------------------------------------\n")
    cat(sprintf("  %-8s %9s %9s %9s\n", "", "p-value", "5% cv", "valid"))
    cat(sprintf("  %-8s %9.4f %9.4f %9d\n", "F_pss", eng$p_value["Fov"], br$cv_F["95%"], eng$n_valid["Fov"]))
    cat(sprintf("  %-8s %9.4f %9.4f %9d\n", "t_dep", eng$p_value["t"], br$cv_t["5%"], eng$n_valid["t"]))
    if (!is.na(object$F_ind))
      cat(sprintf("  %-8s %9.4f %9.4f %9d\n", "F_ind", eng$p_value["Find"], br$cv_F_ind["95%"], eng$n_valid["Find"]))
  } else {
    if (object$conclusion$method == "none") {
      cat("\nBounds for the model without Fourier terms; not valid with\n")
      cat("Fourier terms (reference only, 5% level):\n")
    } else {
      cat("\nKripfganz and Schneider (2020) bounds (5% level):\n")
    }
    cat("-----------------------------------------------\n")
    cat(sprintf("  F: I(0) = %.3f, I(1) = %.3f\n",
                object$critical_values$F$I0["95%"],
                object$critical_values$F$I1["95%"]))
    cat(sprintf("  t: I(0) = %.3f, I(1) = %.3f\n",
                object$critical_values$t$I0["95%"],
                object$critical_values$t$I1["95%"]))
  }

  cat("\nLong-Run Coefficients:\n")
  cat("-----------------------------------------------\n")
  print(round(object$long_run, 4))
  
  cat("\nModel Fit:\n")
  cat("-----------------------------------------------\n")
  cat(sprintf("  R-squared:     %.4f\n", object$fit$R2))
  cat(sprintf("  Adj R-squared: %.4f\n", object$fit$adj_R2))
  cat(sprintf("  AIC:           %.2f\n", object$fit$AIC))
  cat(sprintf("  BIC:           %.2f\n", object$fit$BIC))
  
  cat("\nDiagnostic Tests (p-values):\n")
  cat("-----------------------------------------------\n")
  cat(sprintf("  Serial correlation (BG): %.4f\n", object$diagnostics$serial_corr$p.value))
  cat(sprintf("  Heteroskedasticity (BP): %.4f\n", object$diagnostics$heteroskedasticity$p.value))
  cat(sprintf("  Normality (Shapiro):     %.4f\n", object$diagnostics$normality))
  
  cat("\n===============================================\n")
  cat("CONCLUSION:", object$conclusion$decision, "\n")
  cat(object$conclusion$message, "\n")
  cat("===============================================\n\n")
  
  invisible(object)
}
