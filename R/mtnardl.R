#' @title Multiple-Threshold Nonlinear ARDL (MT-NARDL)
#' @description Extends NARDL to allow multiple threshold decomposition
#' for capturing complex asymmetric relationships.
#'
#' @details
#' The Multiple-Threshold NARDL (MT-NARDL) model extends the standard NARDL
#' framework by allowing decomposition of variables into multiple regimes
#' based on user-specified thresholds. This captures more nuanced asymmetric
#' effects beyond simple positive/negative decomposition.
#'
#' For example, with thresholds c(-0.02, 0, 0.02), a variable is decomposed into:
#' \itemize{
#'   \item Large decreases (< -2\%)
#'   \item Small decreases (-2\% to 0)
#'   \item Small increases (0 to 2\%)
#'   \item Large increases (> 2\%)
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param thresholds Numeric vector of threshold values (default: c(0))
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param case Integer from 1-5 specifying deterministic components (default: 3)
#' @param auto_select Logical. Automatically select optimal thresholds (default: FALSE)
#' @param n_thresholds Integer. Number of thresholds to select if auto_select = TRUE
#' @param bootstrap Logical. Use bootstrap inference (default: FALSE)
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param seed Random seed for reproducibility
#'
#' @return An object of class "mtnardl" containing:
#' \itemize{
#'   \item \code{model}: The estimated MT-NARDL model
#'   \item \code{bounds_test}: Bounds test results
#'   \item \code{long_run}: Long-run coefficients for each regime
#'   \item \code{short_run}: Short-run coefficients for each regime
#'   \item \code{thresholds}: Threshold values used
#'   \item \code{asymmetry_tests}: Wald tests for asymmetry between regimes
#'   \item \code{multipliers}: Dynamic multipliers for each regime
#' }
#'
#' @references
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In Festschrift in Honor of Peter Schmidt (pp. 281-314). Springer.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' data <- generate_oil_data(n = 300)
#'
#' # Standard NARDL (single threshold at 0)
#' result1 <- mtnardl(consumption ~ oil_price, data = data)
#'
#' # Multiple thresholds for different shock sizes
#' result2 <- mtnardl(
#'   consumption ~ oil_price,
#'   data = data,
#'   thresholds = c(-0.05, 0, 0.05)
#' )
#' summary(result2)
#'
#' # Auto-select optimal thresholds
#' result3 <- mtnardl(
#'   consumption ~ oil_price,
#'   data = data,
#'   auto_select = TRUE,
#'   n_thresholds = 2
#' )
#' }
#'
#' @export
mtnardl <- function(formula, data, thresholds = c(0), p = 1, q = 1, case = 3,
                    auto_select = FALSE, n_thresholds = 2,
                    bootstrap = FALSE, nboot = 2000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  thresholds <- sort(unique(thresholds))
  
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
  
  # Auto-select thresholds if requested
  if (auto_select) {
    thresholds <- .select_optimal_thresholds(y, X, n_thresholds, p, q, case)
  }
  
  # Number of regimes = number of thresholds + 1
  n_regimes <- length(thresholds) + 1
  
  # Decompose variables into regimes
  X_decomp <- .mt_decompose(X, thresholds)
  regime_names <- .make_regime_names(x_vars, thresholds)
  colnames(X_decomp) <- regime_names
  
  k_total <- ncol(X_decomp)
  q_expanded <- rep(q, each = n_regimes)
  
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
  
  # Decomposed variables: levels and differences
  x_levels <- matrix(NA, n_valid, k_total)
  x_diff_list <- list()
  
  for (j in 1:k_total) {
    x_j <- X_decomp[, j]
    x_levels[, j] <- x_j[valid_idx - 1]
    
    dx_j <- diff(x_j)
    q_j <- q_expanded[min(j, length(q_expanded))]
    x_diff_j <- matrix(NA, n_valid, q_j)
    for (i in 0:(q_j - 1)) {
      x_diff_j[, i + 1] <- dx_j[(max_lag - i):(n - 1 - i)]
    }
    x_diff_list[[j]] <- x_diff_j
  }
  colnames(x_levels) <- regime_names
  
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
  
  # Estimate model
  model <- stats::lm(dy ~ design - 1)
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  
  # Extract EC coefficient
  ec_coef <- coefs[1]
  ec_se <- sqrt(vcov_mat[1, 1])
  
  # === Bounds Test ===
  n_level_vars <- 1 + k_total
  beta_h0 <- coefs[1:n_level_vars]
  V_h0 <- vcov_mat[1:n_level_vars, 1:n_level_vars]
  
  F_stat <- as.numeric(t(beta_h0) %*% solve(V_h0) %*% beta_h0 / n_level_vars)
  t_stat <- ec_coef / ec_se
  
  # Critical values (use effective k = original k for PSS)
  cv <- pss_critical_values(k, case)
  
  # Bootstrap if requested
  boot_results <- NULL
  if (bootstrap) {
    boot_results <- .mtnardl_bootstrap(dy, design, n_level_vars, nboot)
  }
  
  # Bounds test conclusion
  bounds_test <- .mtnardl_bounds_conclusion(F_stat, t_stat, cv, boot_results)
  
  # === Long-run coefficients by regime ===
  if (abs(ec_coef) > 1e-10) {
    lr_coefs <- -coefs[2:(1 + k_total)] / ec_coef
    names(lr_coefs) <- regime_names
    
    # Organize by original variable
    lr_by_var <- list()
    for (v in x_vars) {
      idx <- grep(paste0("^", v, "_"), regime_names)
      lr_by_var[[v]] <- lr_coefs[idx]
    }
  } else {
    lr_coefs <- rep(NA, k_total)
    lr_by_var <- NULL
  }
  
  # === Asymmetry Tests ===
  asymmetry_tests <- .test_regime_asymmetry(coefs, vcov_mat, x_vars, 
                                             regime_names, thresholds, k_total)
  
  # === Dynamic Multipliers ===
  multipliers <- .compute_mt_multipliers(coefs, ec_coef, regime_names, 
                                          k_total, n_regimes, horizons = 30)
  
  # === Model fit statistics ===
  fit_stats <- list(
    R2 = summary(model)$r.squared,
    adj_R2 = summary(model)$adj.r.squared,
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    sigma = summary(model)$sigma,
    df = model$df.residual
  )
  
  # Build result
  result <- list(
    model = model,
    bounds_test = bounds_test,
    F_stat = F_stat,
    t_stat = t_stat,
    critical_values = cv,
    boot_results = boot_results,
    long_run = lr_coefs,
    long_run_by_var = lr_by_var,
    short_run = coefs,
    thresholds = thresholds,
    n_regimes = n_regimes,
    regime_names = regime_names,
    asymmetry_tests = asymmetry_tests,
    multipliers = multipliers,
    fit = fit_stats,
    call = match.call(),
    case = case,
    n = n_valid,
    k = k,
    p = p,
    q = q
  )
  
  class(result) <- "mtnardl"
  return(result)
}


#' @title Multiple-Threshold Decomposition
#' @keywords internal
.mt_decompose <- function(X, thresholds) {
  X <- as.matrix(X)
  k <- ncol(X)
  n <- nrow(X)
  n_regimes <- length(thresholds) + 1
  
  result <- matrix(0, n, k * n_regimes)
  
  for (j in 1:k) {
    dx <- c(0, diff(X[, j]))
    
    for (r in 1:n_regimes) {
      col_idx <- (j - 1) * n_regimes + r
      
      if (r == 1) {
        # First regime: below lowest threshold
        regime_change <- ifelse(dx <= thresholds[1], dx, 0)
      } else if (r == n_regimes) {
        # Last regime: above highest threshold
        regime_change <- ifelse(dx > thresholds[length(thresholds)], dx, 0)
      } else {
        # Middle regimes: between thresholds
        lower <- thresholds[r - 1]
        upper <- thresholds[r]
        regime_change <- ifelse(dx > lower & dx <= upper, dx, 0)
      }
      
      result[, col_idx] <- cumsum(regime_change)
    }
  }
  
  return(result)
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
#' @keywords internal
.select_optimal_thresholds <- function(y, X, n_thresholds, p, q, case) {
  # Grid search over potential thresholds
  dx <- diff(X[, 1])  # Use first variable for threshold selection
  
  # Generate candidate thresholds (quantiles of changes)
  quantiles <- seq(0.1, 0.9, by = 0.1)
  candidates <- stats::quantile(dx, quantiles)
  
  if (n_thresholds == 1) {
    # Include 0 as candidate
    candidates <- sort(unique(c(0, candidates)))
    
    best_aic <- Inf
    best_threshold <- 0
    
    for (th in candidates) {
      tryCatch({
        test_model <- mtnardl(
          as.formula(paste(names(y), "~", paste(colnames(X), collapse = "+"))),
          data = data.frame(y, X),
          thresholds = th,
          p = p, q = q, case = case
        )
        if (test_model$fit$AIC < best_aic) {
          best_aic <- test_model$fit$AIC
          best_threshold <- th
        }
      }, error = function(e) NULL)
    }
    
    return(best_threshold)
  }
  
  # For multiple thresholds, use percentile-based selection
  threshold_quantiles <- seq(1/(n_thresholds + 1), 
                             n_thresholds/(n_thresholds + 1), 
                             length.out = n_thresholds)
  thresholds <- stats::quantile(dx, threshold_quantiles)
  
  return(as.numeric(thresholds))
}


#' @title MT-NARDL Bootstrap
#' @keywords internal
.mtnardl_bootstrap <- function(dy, design, n_level, nboot) {
  n <- length(dy)
  
  model_null <- stats::lm(dy ~ design[, -(1:n_level)] - 1)
  resid_null <- stats::residuals(model_null)
  fitted_null <- stats::fitted(model_null)
  
  boot_F <- numeric(nboot)
  boot_t <- numeric(nboot)
  
  for (b in 1:nboot) {
    boot_resid <- sample(resid_null, n, replace = TRUE)
    boot_y <- fitted_null + boot_resid
    
    boot_model <- stats::lm(boot_y ~ design - 1)
    boot_coefs <- stats::coef(boot_model)
    boot_vcov <- stats::vcov(boot_model)
    
    beta_h0 <- boot_coefs[1:n_level]
    V_h0 <- boot_vcov[1:n_level, 1:n_level]
    
    boot_F[b] <- tryCatch({
      as.numeric(t(beta_h0) %*% solve(V_h0) %*% beta_h0 / n_level)
    }, error = function(e) NA)
    
    boot_t[b] <- boot_coefs[1] / sqrt(boot_vcov[1, 1])
  }
  
  boot_F <- boot_F[!is.na(boot_F)]
  boot_t <- boot_t[!is.na(boot_t)]
  
  list(
    F_dist = boot_F,
    t_dist = boot_t,
    cv_F = stats::quantile(boot_F, c(0.90, 0.95, 0.99)),
    cv_t = stats::quantile(boot_t, c(0.10, 0.05, 0.01))
  )
}


#' @title MT-NARDL Bounds Test Conclusion
#' @keywords internal
.mtnardl_bounds_conclusion <- function(F_stat, t_stat, cv, boot = NULL) {
  
  if (!is.null(boot)) {
    p_F <- mean(boot$F_dist >= F_stat)
    p_t <- mean(boot$t_dist <= t_stat)
    
    if (p_F < 0.05 && p_t < 0.05) {
      decision <- "COINTEGRATION"
      message <- "Cointegration confirmed (bootstrap): Both F and t significant at 5%"
    } else {
      decision <- "NO_COINTEGRATION"
      message <- "No cointegration (bootstrap): Failed to reject null"
    }
    
    return(list(decision = decision, message = message, 
                p_values = c(F = p_F, t = p_t), method = "bootstrap"))
  }
  
  F_upper <- cv$F_I1["5%"]
  F_lower <- cv$F_I0["5%"]
  
  if (F_stat > F_upper) {
    decision <- "COINTEGRATION"
    message <- "Cointegration confirmed: F > I(1) upper bound"
  } else if (F_stat < F_lower) {
    decision <- "NO_COINTEGRATION"
    message <- "No cointegration: F < I(0) lower bound"
  } else {
    decision <- "INCONCLUSIVE"
    message <- "Inconclusive: F between I(0) and I(1) bounds"
  }
  
  list(decision = decision, message = message, method = "asymptotic")
}


#' @title Test Regime Asymmetry
#' @keywords internal
.test_regime_asymmetry <- function(coefs, vcov_mat, x_vars, regime_names, 
                                    thresholds, k_total) {
  n_regimes <- length(thresholds) + 1
  tests <- list()
  
  for (v in x_vars) {
    # Get indices for this variable's regimes
    idx <- grep(paste0("^", v, "_"), regime_names)
    if (length(idx) < 2) next
    
    # Indices in coefficient vector (offset by 1 for EC term)
    coef_idx <- idx + 1
    
    # Pairwise Wald tests
    var_tests <- list()
    for (i in 1:(n_regimes - 1)) {
      for (j in (i + 1):n_regimes) {
        idx_i <- coef_idx[i]
        idx_j <- coef_idx[j]
        
        diff_coef <- coefs[idx_i] - coefs[idx_j]
        var_diff <- vcov_mat[idx_i, idx_i] + vcov_mat[idx_j, idx_j] - 
                   2 * vcov_mat[idx_i, idx_j]
        
        wald <- diff_coef^2 / var_diff
        p_value <- 1 - stats::pchisq(wald, 1)
        
        test_name <- paste0("regime", i, "_vs_regime", j)
        var_tests[[test_name]] <- list(
          wald = wald,
          p_value = p_value,
          diff = diff_coef,
          significant = p_value < 0.05
        )
      }
    }
    
    tests[[v]] <- var_tests
  }
  
  return(tests)
}


#' @title Compute Multiple-Threshold Dynamic Multipliers
#' @keywords internal
.compute_mt_multipliers <- function(coefs, ec_coef, regime_names, k_total, 
                                     n_regimes, horizons = 30) {
  multipliers <- list()
  
  # Long-run multipliers by regime
  lr_mult <- if (abs(ec_coef) > 1e-10) {
    -coefs[2:(1 + k_total)] / ec_coef
  } else {
    rep(NA, k_total)
  }
  names(lr_mult) <- regime_names
  
  # Cumulative multipliers over horizons
  cum_mult <- matrix(NA, horizons, k_total)
  colnames(cum_mult) <- regime_names
  
  for (j in 1:k_total) {
    sr_coef <- coefs[1 + j]  # Short-run coefficient
    adjustment <- ec_coef
    
    for (h in 1:horizons) {
      if (h == 1) {
        cum_mult[h, j] <- sr_coef
      } else {
        cum_mult[h, j] <- cum_mult[h-1, j] + sr_coef * (1 + adjustment)^(h-1)
      }
    }
  }
  
  list(
    long_run = lr_mult,
    cumulative = cum_mult,
    horizons = 1:horizons
  )
}


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
  cat("\nDecision:", x$bounds_test$decision, "\n")
  
  invisible(x)
}


#' @export
summary.mtnardl <- function(object, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("    Multiple-Threshold NARDL Estimation Results\n")
  cat("=======================================================\n\n")
  
  cat("Model Specification:\n")
  cat("  Case:", object$case, "\n")
  cat("  Thresholds:", paste(round(object$thresholds, 4), collapse = ", "), "\n")
  cat("  Number of regimes:", object$n_regimes, "\n")
  cat("  Lags: p =", object$p, ", q =", paste(object$q, collapse = ","), "\n")
  cat("  Sample size:", object$n, "\n\n")
  
  cat("Bounds Test:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  F-statistic: %10.4f\n", object$F_stat))
  cat(sprintf("  t-statistic: %10.4f\n", object$t_stat))
  cat("\n  Decision:", object$bounds_test$decision, "\n")
  cat(" ", object$bounds_test$message, "\n\n")
  
  cat("Long-Run Coefficients by Regime:\n")
  cat("-------------------------------------------------------\n")
  print(round(object$long_run, 4))
  
  cat("\nAsymmetry Tests:\n")
  cat("-------------------------------------------------------\n")
  for (v in names(object$asymmetry_tests)) {
    cat("\n", v, ":\n")
    for (test_name in names(object$asymmetry_tests[[v]])) {
      test <- object$asymmetry_tests[[v]][[test_name]]
      sig <- if (test$significant) "***" else ""
      cat(sprintf("  %-20s Wald = %6.3f, p = %.4f %s\n", 
                  test_name, test$wald, test$p_value, sig))
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


#' @export
plot.mtnardl <- function(x, type = c("multipliers", "asymmetry"), ...) {
  type <- match.arg(type)
  
  if (type == "multipliers") {
    # Plot cumulative dynamic multipliers
    mult <- x$multipliers$cumulative
    horizons <- x$multipliers$horizons
    
    n_vars <- ncol(mult)
    colors <- grDevices::rainbow(n_vars)
    
    graphics::par(mfrow = c(1, 1))
    graphics::matplot(horizons, mult, type = "l", lty = 1, col = colors,
                     xlab = "Horizon", ylab = "Cumulative Multiplier",
                     main = "Dynamic Multipliers by Regime")
    graphics::legend("topright", legend = colnames(mult), col = colors, 
                    lty = 1, cex = 0.7)
    graphics::abline(h = 0, lty = 2, col = "gray")
  }
  
  if (type == "asymmetry") {
    # Plot long-run coefficients comparison
    lr <- x$long_run
    
    graphics::barplot(lr, las = 2, col = "steelblue",
                     main = "Long-Run Coefficients by Regime",
                     ylab = "Coefficient")
    graphics::abline(h = 0, lty = 2)
  }
  
  invisible(x)
}
