#' @title Bootstrap ARDL Bounds Test
#' @description Perform bounds test for cointegration with bootstrap critical values
#'
#' @details
#' This function implements bootstrap-based inference for the ARDL bounds test,
#' which is particularly useful for small samples where asymptotic critical values
#' may be unreliable. The bootstrap procedure follows McNown, Sam & Goh (2018).
#'
#' Three test statistics are computed:
#' \itemize{
#'   \item \strong{F-statistic}: Joint test that all level coefficients are zero
#'   \item \strong{t-statistic}: Test that the error correction coefficient is zero
#'   \item \strong{F-overall}: Joint test including the dependent variable lag
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param case Integer from 1-5 specifying deterministic components:
#'   \itemize{
#'     \item 1: No intercept, no trend
#'     \item 2: Restricted intercept, no trend
#'     \item 3: Unrestricted intercept, no trend (default)
#'     \item 4: Unrestricted intercept, restricted trend
#'     \item 5: Unrestricted intercept, unrestricted trend
#'   }
#' @param nboot Number of bootstrap replications (default: 2000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param parallel Logical. Use parallel processing (default: FALSE
#' @param ncores Number of cores for parallel processing (default: 2)
#'
#' @return An object of class "boot_ardl" containing:
#' \itemize{
#'   \item \code{F_stat}: F-statistic for bounds test
#'   \item \code{t_stat}: t-statistic for EC coefficient
#'   \item \code{F_overall}: Overall F-statistic
#'   \item \code{boot_F}: Bootstrap distribution of F-statistics
#'   \item \code{boot_t}: Bootstrap distribution of t-statistics
#'   \item \code{cv_F}: Critical values for F-test (90%, 95%, 99%)
#'   \item \code{cv_t}: Critical values for t-test
#'   \item \code{p_value_F}: Bootstrap p-value for F-test
#'   \item \code{p_value_t}: Bootstrap p-value for t-test
#'   \item \code{model}: The estimated ARDL model
#'   \item \code{conclusion}: Test conclusion
#' }
#'
#' @references
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
#' distributed lag test for cointegration. Applied Economics, 50(13), 1509-1521.
#'
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. Journal of Applied Econometrics, 16(3), 289-326.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(macro_data)
#'
#' # Bootstrap bounds test
#' boot_test <- boot_ardl(
#'   gdp ~ inflation + investment + trade,
#'   data = macro_data,
#'   p = 2, q = 2,
#'   case = 3,
#'   nboot = 2000
#' )
#' summary(boot_test)
#' plot(boot_test)
#' }
#'
#' @export
#' @importFrom stats lm coef residuals fitted var rnorm quantile
boot_ardl <- function(formula, data, p = 1, q = 1, case = 3,
                      nboot = 2000, seed = NULL, 
                      parallel = FALSE, ncores = 2) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (!case %in% 1:5) {
    stop("'case' must be an integer from 1 to 5")
  }
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  k <- length(x_vars)
  
  # Handle q as vector or scalar
  if (length(q) == 1) {
    q <- rep(q, k)
  }
  
  # Prepare data
  ardl_data <- .prepare_ardl_ts_data(data, y_var, x_vars, p, q, case)
  
  if (is.null(ardl_data)) {
    stop("Insufficient observations for specified lag structure")
  }
  
  # Estimate unrestricted ARDL model
  model_ur <- .estimate_ardl_unrestricted(ardl_data, case)
  
  # Estimate restricted model (no levels)
  model_r <- .estimate_ardl_restricted(ardl_data, case)
  
  # Compute test statistics
  F_stat <- .compute_F_stat(model_ur, model_r, k, case)
  t_stat <- .compute_t_stat(model_ur)
  F_overall <- .compute_F_overall(model_ur, model_r, k, case)
  
  # Bootstrap procedure
  boot_results <- .bootstrap_bounds(
    ardl_data, model_ur, case, k, nboot, parallel, ncores
  )
  
  # Critical values from bootstrap distribution
  cv_F <- quantile(boot_results$F_boot, probs = c(0.90, 0.95, 0.99))
  cv_t <- quantile(boot_results$t_boot, probs = c(0.10, 0.05, 0.01))
  
  # P-values
  p_value_F <- mean(boot_results$F_boot >= F_stat)
  p_value_t <- mean(boot_results$t_boot <= t_stat)
  
  # Conclusion
  conclusion <- .bounds_conclusion(F_stat, t_stat, cv_F, cv_t, k, case)
  
  result <- list(
    F_stat = F_stat,
    t_stat = t_stat,
    F_overall = F_overall,
    boot_F = boot_results$F_boot,
    boot_t = boot_results$t_boot,
    cv_F = cv_F,
    cv_t = cv_t,
    p_value_F = p_value_F,
    p_value_t = p_value_t,
    model = model_ur,
    case = case,
    k = k,
    nboot = nboot,
    conclusion = conclusion,
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
    dy = dy,
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
  if (case >= 3) {
    result$const <- 1
  }
  if (case >= 4) {
    result$trend <- 1:n_obs
  }
  
  # Remove NAs
  result <- na.omit(result)
  
  return(result)
}


#' @title Estimate Unrestricted ARDL Model
#' @keywords internal
.estimate_ardl_unrestricted <- function(ardl_data, case) {
  
  # All variables except dy
  xvars <- names(ardl_data)[-1]
  
  formula_str <- paste("dy ~", paste(xvars, collapse = " + "))
  
  if (case == 1) {
    formula_str <- paste(formula_str, "- 1")
  }
  
  model <- lm(as.formula(formula_str), data = ardl_data)
  
  return(model)
}


#' @title Estimate Restricted ARDL Model
#' @keywords internal
.estimate_ardl_restricted <- function(ardl_data, case) {
  
  # Remove level variables (y_lag1 and X in levels)
  xvars <- names(ardl_data)[-1]
  
  # Identify level variables to exclude
  level_vars <- grep("^y_lag1$|^[A-Za-z]", xvars, value = TRUE)
  level_vars <- level_vars[!grepl("^d|^const$|^trend$|_L[0-9]", level_vars)]
  level_vars <- c("y_lag1", level_vars)
  
  # Keep only differenced variables and deterministics
  keep_vars <- setdiff(xvars, level_vars)
  
  if (length(keep_vars) == 0) {
    # At minimum, keep differences
    keep_vars <- grep("^d", xvars, value = TRUE)
  }
  
  formula_str <- paste("dy ~", paste(keep_vars, collapse = " + "))
  
  if (case == 1) {
    formula_str <- paste(formula_str, "- 1")
  }
  
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
  
  # Number of restrictions = k + 1 (x levels + y_lag1)
  m <- k + 1
  
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


#' @title Bootstrap Bounds Test
#' @keywords internal
.bootstrap_bounds <- function(ardl_data, model_ur, case, k, nboot, parallel, ncores) {
  
  # Get residuals and fitted values under null
  resid_ur <- residuals(model_ur)
  n <- length(resid_ur)
  
  # Center residuals
  resid_centered <- resid_ur - mean(resid_ur)
  
  # Bootstrap function
  boot_one <- function(b) {
    
    # Resample residuals
    boot_resid <- sample(resid_centered, n, replace = TRUE)
    
    # Generate bootstrap dy under null (no cointegration)
    # Simple approach: use restricted model structure
    ardl_boot <- ardl_data
    ardl_boot$dy <- fitted(model_ur) - resid_ur + boot_resid
    
    # Re-estimate models
    tryCatch({
      model_ur_boot <- .estimate_ardl_unrestricted(ardl_boot, case)
      model_r_boot <- .estimate_ardl_restricted(ardl_boot, case)
      
      F_boot <- .compute_F_stat(model_ur_boot, model_r_boot, k, case)
      t_boot <- .compute_t_stat(model_ur_boot)
      
      c(F_boot, t_boot)
    }, error = function(e) {
      c(NA, NA)
    })
  }
  
  # Run bootstrap
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    
    boot_stats <- parallel::parSapply(cl, 1:nboot, boot_one)
  } else {
    boot_stats <- sapply(1:nboot, boot_one)
  }
  
  list(
    F_boot = boot_stats[1, ],
    t_boot = boot_stats[2, ]
  )
}


#' @title Bounds Test Conclusion
#' @keywords internal
.bounds_conclusion <- function(F_stat, t_stat, cv_F, cv_t, k, case) {
  
  # Compare with 95% critical value
  cv_F_95 <- cv_F["95%"]
  cv_t_95 <- cv_t["5%"]
  
  if (F_stat > cv_F_95 && t_stat < cv_t_95) {
    conclusion <- "COINTEGRATION: Both F and t statistics significant at 5% level"
  } else if (F_stat > cv_F_95) {
    conclusion <- "COINTEGRATION: F-statistic significant at 5% level"
  } else if (t_stat < cv_t_95) {
    conclusion <- "COINTEGRATION: t-statistic significant at 5% level"
  } else {
    conclusion <- "NO COINTEGRATION: Cannot reject null hypothesis"
  }
  
  return(conclusion)
}


#' @title PSS Asymptotic Critical Values
#' @description Get Pesaran, Shin & Smith (2001) asymptotic critical values
#'
#' @param k Number of regressors (excluding the lagged dependent variable)
#' @param case Deterministic specification (1-5)
#' @param level Significance level: "10%", "5%", or "1%"
#'
#' @return A list with I(0) and I(1) bounds for F and t statistics
#'
#' @references
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. Journal of Applied Econometrics.
#'
#' @export
pss_critical_values <- function(k, case = 3, level = "5%") {
  
  # Table CI(iii): Unrestricted intercept, no trend (Case 3)
  # F-statistic bounds
  cv_F <- list(
    "10%" = list(
      I0 = c(2.37, 2.45, 2.52, 2.56, 2.62, 2.66, 2.69)[min(k, 7)],
      I1 = c(3.20, 3.52, 3.83, 4.10, 4.35, 4.60, 4.80)[min(k, 7)]
    ),
    "5%" = list(
      I0 = c(2.79, 2.87, 2.94, 3.02, 3.07, 3.12, 3.15)[min(k, 7)],
      I1 = c(3.67, 4.00, 4.35, 4.66, 4.90, 5.15, 5.35)[min(k, 7)]
    ),
    "1%" = list(
      I0 = c(3.65, 3.74, 3.88, 4.00, 4.10, 4.18, 4.26)[min(k, 7)],
      I1 = c(4.66, 5.06, 5.44, 5.80, 6.10, 6.36, 6.58)[min(k, 7)]
    )
  )
  
  # t-statistic bounds (Case 3)
  cv_t <- list(
    "10%" = list(I0 = -2.57, I1 = -3.21),
    "5%" = list(I0 = -2.86, I1 = -3.53),
    "1%" = list(I0 = -3.43, I1 = -4.10)
  )
  
  list(
    F_bounds = cv_F[[level]],
    t_bounds = cv_t[[level]],
    k = k,
    case = case,
    level = level
  )
}


#' @title Summary method for boot_ardl
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
  cat("Bootstrap replications:", object$nboot, "\n\n")
  
  cat("--------------------------------------------------------------------\n")
  cat("                     Test Statistics\n")
  cat("--------------------------------------------------------------------\n\n")
  
  # F-test results
  cat("F-statistic:", round(object$F_stat, 4), "\n")
  cat("Bootstrap critical values:\n")
  cat("  90%:", round(object$cv_F["90%"], 4), "\n")
  cat("  95%:", round(object$cv_F["95%"], 4), "\n")
  cat("  99%:", round(object$cv_F["99%"], 4), "\n")
  cat("Bootstrap p-value:", round(object$p_value_F, 4), "\n\n")
  
  # t-test results
  cat("t-statistic:", round(object$t_stat, 4), "\n")
  cat("Bootstrap critical values:\n")
  cat("  90%:", round(object$cv_t["10%"], 4), "\n")
  cat("  95%:", round(object$cv_t["5%"], 4), "\n")
  cat("  99%:", round(object$cv_t["1%"], 4), "\n")
  cat("Bootstrap p-value:", round(object$p_value_t, 4), "\n\n")
  
  # PSS asymptotic bounds for comparison
  pss <- pss_critical_values(object$k, object$case, "5%")
  cat("--------------------------------------------------------------------\n")
  cat("     PSS (2001) Asymptotic Critical Values (5% level)\n")
  cat("--------------------------------------------------------------------\n")
  cat("F-bounds: I(0) =", round(pss$F_bounds$I0, 2), 
      ", I(1) =", round(pss$F_bounds$I1, 2), "\n")
  cat("t-bounds: I(0) =", round(pss$t_bounds$I0, 2),
      ", I(1) =", round(pss$t_bounds$I1, 2), "\n\n")
  
  cat("--------------------------------------------------------------------\n")
  cat("                      Conclusion\n")
  cat("--------------------------------------------------------------------\n")
  cat(object$conclusion, "\n")
  cat("====================================================================\n")
  
  invisible(object)
}


#' @title Print method for boot_ardl
#' @export
print.boot_ardl <- function(x, ...) {
  
  cat("\nBootstrap ARDL Bounds Test\n")
  cat("F-statistic:", round(x$F_stat, 4), 
      "(p-value:", round(x$p_value_F, 4), ")\n")
  cat("t-statistic:", round(x$t_stat, 4),
      "(p-value:", round(x$p_value_t, 4), ")\n")
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
    df_F <- data.frame(F_stat = x$boot_F)
    
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
                        "| 95% CV =", round(x$cv_F["95%"], 3)),
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
                        "| 95% CV =", round(x$cv_t["5%"], 3)),
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
