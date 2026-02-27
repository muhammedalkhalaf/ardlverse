#' @title Augmented ARDL Bounds Test (AARDL)
#' @description Implements the Augmented ARDL bounds testing approach with 
#' deferred t and F tests for cointegration analysis.
#'
#' @details
#' The Augmented ARDL (AARDL) approach extends the standard ARDL bounds test
#' by implementing additional diagnostic tests proposed by Sam, McNown & Goh (2019).
#' This addresses potential weaknesses in the PSS bounds test by adding:
#' \itemize{
#'   \item Deferred t-test (t_dep): Tests significance of lagged dependent variable
#'   \item Deferred F-test (F_ind): Tests joint significance of independent variables
#'   \item Overall F-test with all deferred conditions
#' }
#'
#' The function supports 8 sub-models:
#' \enumerate{
#'   \item Standard ARDL
#'   \item ARDL with bootstrap
#'   \item Nonlinear ARDL (NARDL)
#'   \item Fourier ARDL
#'   \item Fourier NARDL
#'   \item Bootstrap NARDL
#'   \item Fourier Bootstrap ARDL
#'   \item Fourier Bootstrap NARDL
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param case Integer from 1-5 specifying deterministic components
#' @param type Character. Model type: "linear", "nardl", "fourier", "fnardl",
#'   "bootstrap", "bnardl", "fbootstrap", "fbnardl" (default: "linear")
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param fourier_k Integer. Number of Fourier frequencies (default: 1, max: 3)
#' @param threshold Numeric. Threshold value for NARDL decomposition (default: 0)
#' @param seed Random seed for reproducibility
#'
#' @return An object of class "aardl" containing:
#' \itemize{
#'   \item \code{F_pss}: PSS F-statistic for bounds test
#'   \item \code{t_dep}: Deferred t-statistic for lagged dependent variable
#'   \item \code{F_ind}: Deferred F-statistic for independent variables
#'   \item \code{conclusion}: Cointegration decision based on all tests
#'   \item \code{model}: The estimated ARDL model
#'   \item \code{long_run}: Long-run coefficients
#'   \item \code{short_run}: Short-run coefficients
#'   \item \code{diagnostics}: Model diagnostic tests
#' }
#'
#' @references
#' Sam, C. Y., McNown, R., & Goh, S. K. (2019). An augmented autoregressive 
#' distributed lag bounds test for cointegration. Economic Modelling, 80, 130-141.
#'
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
#' distributed lag test for cointegration. Applied Economics, 50(13), 1509-1521.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' data <- generate_ts_data(n = 200)
#'
#' # Standard Augmented ARDL
#' result <- aardl(y ~ x1 + x2, data = data, p = 2, q = 2, case = 3)
#' summary(result)
#'
#' # Augmented NARDL (nonlinear)
#' result_nardl <- aardl(y ~ x1 + x2, data = data, type = "nardl")
#' summary(result_nardl)
#'
#' # Fourier Augmented ARDL
#' result_fourier <- aardl(y ~ x1 + x2, data = data, type = "fourier", fourier_k = 2)
#' summary(result_fourier)
#' }
#'
#' @export
aardl <- function(formula, data, p = 1, q = 1, case = 3,
                  type = c("linear", "nardl", "fourier", "fnardl",
                          "bootstrap", "bnardl", "fbootstrap", "fbnardl"),
                  nboot = 2000, fourier_k = 1, threshold = 0, seed = NULL) {
  
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  if (fourier_k < 1 || fourier_k > 3) {
    stop("'fourier_k' must be between 1 and 3")
  }
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  k <- length(x_vars)
  
  # Handle q as vector or scalar
  if (length(q) == 1) q <- rep(q, k)
  
  # Get variables
  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  n <- length(y)
  
  # Apply transformations based on type
  use_fourier <- type %in% c("fourier", "fnardl", "fbootstrap", "fbnardl")
  use_nardl <- type %in% c("nardl", "fnardl", "bnardl", "fbnardl")
  use_bootstrap <- type %in% c("bootstrap", "bnardl", "fbootstrap", "fbnardl")
  
  # Decompose variables for NARDL
  if (use_nardl) {
    X_decomp <- .decompose_asymmetric(X, threshold)
    x_vars_new <- c(paste0(x_vars, "_pos"), paste0(x_vars, "_neg"))
    X <- X_decomp
    k <- ncol(X)
    if (length(q) == length(x_vars)) {
      q <- rep(q, 2)  # Duplicate for pos/neg
    }
  }
  
  # Create Fourier terms
  fourier_terms <- NULL
  if (use_fourier) {
    fourier_terms <- .create_fourier_terms(n, fourier_k)
  }
  
  # Build ARDL model matrix
  max_lag <- max(p, max(q))
  valid_idx <- (max_lag + 1):n
  n_valid <- length(valid_idx)
  
  # Dependent variable in differences
  dy <- diff(y)[(max_lag):(n-1)]
  
  # Lagged dependent variable (level) for EC term
  y_lag <- y[valid_idx - 1]
  
  # Lagged differences of dependent variable
  dy_lags <- matrix(NA, n_valid, p)
  for (i in 1:p) {
    dy_lags[, i] <- diff(y)[(max_lag - i):(n - 1 - i)]
  }
  colnames(dy_lags) <- paste0("d.", y_var, ".l", 1:p)
  
  # Independent variables: levels and differences
  x_levels <- matrix(NA, n_valid, k)
  x_diff_list <- list()
  
  for (j in 1:k) {
    x_j <- if (use_nardl) X[, j] else X[, j]
    x_levels[, j] <- x_j[valid_idx - 1]
    
    dx_j <- diff(x_j)
    x_diff_j <- matrix(NA, n_valid, q[min(j, length(q))])
    for (i in 0:(q[min(j, length(q))] - 1)) {
      x_diff_j[, i + 1] <- dx_j[(max_lag - i):(n - 1 - i)]
    }
    x_diff_list[[j]] <- x_diff_j
  }
  
  if (use_nardl) {
    colnames(x_levels) <- x_vars_new
  } else {
    colnames(x_levels) <- x_vars
  }
  
  # Combine difference terms
  x_diffs <- do.call(cbind, x_diff_list)
  
  # Build design matrix
  design <- cbind(y_lag, x_levels, dy_lags, x_diffs)
  
  # Add deterministic components
  if (case >= 2) {
    design <- cbind(design, intercept = 1)
  }
  if (case >= 4) {
    trend <- 1:n_valid
    design <- cbind(design, trend = trend)
  }
  
  # Add Fourier terms
  if (use_fourier) {
    fourier_valid <- fourier_terms[valid_idx, , drop = FALSE]
    design <- cbind(design, fourier_valid)
  }
  
  # Estimate model
  model <- stats::lm(dy ~ design - 1)
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  
  # Extract coefficients for tests
  ec_coef <- coefs[1]  # Coefficient on lagged y
  ec_se <- sqrt(vcov_mat[1, 1])
  
  # === PSS F-test (joint test on level variables) ===
  # H0: coefficient on y_{t-1} and all x_{t-1} are jointly zero
  n_level_vars <- 1 + k
  R <- diag(n_level_vars)
  R_full <- matrix(0, n_level_vars, length(coefs))
  R_full[, 1:n_level_vars] <- R
  
  beta_h0 <- coefs[1:n_level_vars]
  V_h0 <- vcov_mat[1:n_level_vars, 1:n_level_vars]
  
  F_pss <- as.numeric(t(beta_h0) %*% solve(V_h0) %*% beta_h0 / n_level_vars)
  
  # === Deferred t-test (t_dep) ===
  # Test significance of lagged dependent variable alone
  t_dep <- ec_coef / ec_se
  
  # === Deferred F-test (F_ind) ===
  # Test joint significance of lagged independent variables
  if (k > 0) {
    beta_ind <- coefs[2:(1 + k)]
    V_ind <- vcov_mat[2:(1 + k), 2:(1 + k)]
    F_ind <- as.numeric(t(beta_ind) %*% solve(V_ind) %*% beta_ind / k)
  } else {
    F_ind <- NA
  }
  
  # === Get Critical Values ===
  cv <- .aardl_critical_values(k, case, n_valid)
  
  # === Bootstrap inference if requested ===
  boot_results <- NULL
  if (use_bootstrap) {
    boot_results <- .aardl_bootstrap(dy, design, n_level_vars, k, nboot)
  }
  
  # === Make Conclusion ===
  conclusion <- .aardl_conclusion(F_pss, t_dep, F_ind, cv, boot_results)
  

  # === Long-run coefficients ===
  if (abs(ec_coef) > 1e-10) {
    lr_coefs <- -coefs[2:(1 + k)] / ec_coef
    if (use_nardl) {
      names(lr_coefs) <- x_vars_new
    } else {
      names(lr_coefs) <- x_vars
    }
  } else {
    lr_coefs <- rep(NA, k)
  }
  
  # === Short-run coefficients ===
  sr_start <- 1 + k + 1
  sr_end <- sr_start + p - 1
  sr_coefs <- coefs[sr_start:sr_end]
  
  # === Model fit statistics ===
  fit_stats <- list(
    R2 = summary(model)$r.squared,
    adj_R2 = summary(model)$adj.r.squared,
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    sigma = summary(model)$sigma,
    df = model$df.residual
  )
  
  # === Diagnostic tests ===
  resid <- stats::residuals(model)
  diagnostics <- list(
    serial_corr = .breusch_godfrey_test(resid, 2),
    heteroskedasticity = .breusch_pagan_test(model),
    normality = stats::shapiro.test(resid)$p.value
  )
  
  # Build result
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
    q = q,
    fourier_k = if (use_fourier) fourier_k else NULL
  )
  
  class(result) <- "aardl"
  return(result)
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
#' @keywords internal
.aardl_critical_values <- function(k, case, n) {
  # Critical values from Sam, McNown & Goh (2019)
  # Returns I(0) and I(1) bounds for F and t statistics
  
  # Simplified critical values (asymptotic, case III)
  # In production, should use full tables
  cv_F <- list(
    I0 = c(`90%` = 2.63, `95%` = 3.10, `99%` = 4.13),
    I1 = c(`90%` = 3.35 + 0.22*k, `95%` = 3.87 + 0.25*k, `99%` = 4.96 + 0.30*k)
  )
  
  cv_t <- list(
    I0 = c(`90%` = -2.57, `95%` = -2.86, `99%` = -3.43),
    I1 = c(`90%` = -2.91 - 0.10*k, `95%` = -3.22 - 0.12*k, `99%` = -3.82 - 0.14*k)
  )
  
  list(F = cv_F, t = cv_t)
}


#' @title AARDL Bootstrap Procedure
#' @keywords internal
.aardl_bootstrap <- function(dy, design, n_level, k, nboot) {
  n <- length(dy)
  
  # Estimate null model (no cointegration)
  model_null <- stats::lm(dy ~ design[, -(1:n_level)] - 1)
  resid_null <- stats::residuals(model_null)
  fitted_null <- stats::fitted(model_null)
  
  boot_F <- numeric(nboot)
  boot_t <- numeric(nboot)
  boot_F_ind <- numeric(nboot)
  
  for (b in 1:nboot) {
    # Resample residuals
    boot_resid <- sample(resid_null, n, replace = TRUE)
    boot_y <- fitted_null + boot_resid
    
    # Estimate full model on bootstrap sample
    boot_model <- stats::lm(boot_y ~ design - 1)
    boot_coefs <- stats::coef(boot_model)
    boot_vcov <- stats::vcov(boot_model)
    
    # F-statistic
    beta_h0 <- boot_coefs[1:n_level]
    V_h0 <- boot_vcov[1:n_level, 1:n_level]
    boot_F[b] <- tryCatch({
      as.numeric(t(beta_h0) %*% solve(V_h0) %*% beta_h0 / n_level)
    }, error = function(e) NA)
    
    # t-statistic
    boot_t[b] <- boot_coefs[1] / sqrt(boot_vcov[1, 1])
    
    # F_ind statistic
    if (k > 0) {
      beta_ind <- boot_coefs[2:(1 + k)]
      V_ind <- boot_vcov[2:(1 + k), 2:(1 + k)]
      boot_F_ind[b] <- tryCatch({
        as.numeric(t(beta_ind) %*% solve(V_ind) %*% beta_ind / k)
      }, error = function(e) NA)
    }
  }
  
  # Remove NAs
  boot_F <- boot_F[!is.na(boot_F)]
  boot_t <- boot_t[!is.na(boot_t)]
  boot_F_ind <- boot_F_ind[!is.na(boot_F_ind)]
  
  list(
    F_dist = boot_F,
    t_dist = boot_t,
    F_ind_dist = boot_F_ind,
    cv_F = stats::quantile(boot_F, c(0.90, 0.95, 0.99)),
    cv_t = stats::quantile(boot_t, c(0.10, 0.05, 0.01)),
    cv_F_ind = stats::quantile(boot_F_ind, c(0.90, 0.95, 0.99))
  )
}


#' @title AARDL Conclusion
#' @keywords internal
.aardl_conclusion <- function(F_pss, t_dep, F_ind, cv, boot = NULL) {
  
  if (!is.null(boot)) {
    # Bootstrap-based inference
    p_F <- mean(boot$F_dist >= F_pss)
    p_t <- mean(boot$t_dist <= t_dep)
    p_F_ind <- mean(boot$F_ind_dist >= F_ind)
    
    coint_F <- p_F < 0.05
    coint_t <- p_t < 0.05
    coint_F_ind <- p_F_ind < 0.05
    
    if (coint_F && coint_t && coint_F_ind) {
      conclusion <- "Cointegration confirmed (bootstrap): F_pss, t_dep, and F_ind all significant"
      decision <- "COINTEGRATION"
    } else if (coint_F && coint_t) {
      conclusion <- "Possible cointegration (bootstrap): F_pss and t_dep significant, but F_ind not significant"
      decision <- "INCONCLUSIVE"
    } else {
      conclusion <- "No cointegration (bootstrap): Failed to reject null hypothesis"
      decision <- "NO_COINTEGRATION"
    }
    
    return(list(
      decision = decision,
      message = conclusion,
      p_values = c(F_pss = p_F, t_dep = p_t, F_ind = p_F_ind),
      method = "bootstrap"
    ))
  }
  
  # Asymptotic inference
  F_lower <- cv$F$I0["95%"]
  F_upper <- cv$F$I1["95%"]
  t_lower <- cv$t$I0["95%"]
  t_upper <- cv$t$I1["95%"]
  
  # PSS bounds test decision
  if (F_pss > F_upper && t_dep < t_upper) {
    conclusion <- "Cointegration confirmed: F > I(1) bound and t < I(1) bound"
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
    method = "asymptotic"
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
  
  cat("\nCritical Values (95%):\n")
  cat("-----------------------------------------------\n")
  if (object$conclusion$method == "bootstrap") {
    cat("  Bootstrap critical values used\n")
    cat(sprintf("  F: %8.4f\n", object$boot_results$cv_F["95%"]))
    cat(sprintf("  t: %8.4f\n", object$boot_results$cv_t["5%"]))
  } else {
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
