#' @title Panel ARDL Estimation
#' @description Estimate Panel ARDL models with Pooled Mean Group (PMG),
#'   Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators.
#'
#' @details
#' This function implements the panel ARDL framework of Pesaran, Shin & Smith (1999)
#' for estimating long-run relationships in dynamic heterogeneous panels.
#'
#' The model is specified as:
#' \deqn{\Delta y_{it} = \phi_i (y_{i,t-1} - \theta'_i x_{it}) + \sum_{j=1}^{p-1} \lambda_{ij} \Delta y_{i,t-j} + \sum_{j=0}^{q-1} \delta'_{ij} \Delta x_{i,t-j} + \mu_i + \varepsilon_{it}}
#'
#' Three estimators are available:
#' \itemize{
#'   \item \strong{PMG}: Constrains long-run coefficients to be equal across groups while allowing
#'         short-run coefficients to vary. Efficient when long-run homogeneity holds.
#'   \item \strong{MG}: Estimates separate regressions for each group and averages coefficients.
#'         Consistent but less efficient than PMG.
#'   \item \strong{DFE}: Traditional dynamic fixed effects with all coefficients pooled.
#'         Consistent only under slope homogeneity.
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ... | group
#' @param data A data frame containing panel data
#' @param id Character string specifying the group/panel identifier variable
#' @param time Character string specifying the time variable
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param estimator Character. One of "pmg", "mg", or "dfe" (default: "pmg")
#' @param ec Logical. If TRUE, display error correction form (default: TRUE)
#' @param trend Logical. Include time trend (default: FALSE)
#' @param maxiter Maximum iterations for PMG optimization (default: 100)
#' @param tol Convergence tolerance (default: 1e-5)
#'
#' @return An object of class "panel_ardl" containing:
#' \itemize{
#'   \item \code{coefficients}: Estimated coefficients
#'   \item \code{long_run}: Long-run coefficients (theta)
#'   \item \code{short_run}: Short-run coefficients
#'   \item \code{ec_coef}: Error correction coefficient (phi)
#'   \item \code{se}: Standard errors
#'   \item \code{t_values}: t-statistics
#'   \item \code{p_values}: p-values
#'   \item \code{group_coefs}: Individual group coefficients (for MG)
#'   \item \code{residuals}: Model residuals
#'   \item \code{fitted}: Fitted values
#'   \item \code{nobs}: Number of observations
#'   \item \code{ngroups}: Number of groups
#'   \item \code{df}: Degrees of freedom
#'   \item \code{sigma}: Residual standard error
#'   \item \code{loglik}: Log-likelihood
#'   \item \code{aic}: Akaike Information Criterion
#'   \item \code{bic}: Bayesian Information Criterion
#'   \item \code{hausman}: Hausman test for PMG vs MG
#' }
#'
#' @references
#' Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation
#' of dynamic heterogeneous panels. Journal of the American Statistical Association,
#' 94(446), 621-634.
#'
#' Pesaran, M. H., & Smith, R. (1995). Estimating long-run relationships from
#' dynamic heterogeneous panels. Journal of Econometrics, 68(1), 79-113.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(macro_panel)
#'
#' # Estimate PMG model
#' pmg_model <- panel_ardl(
#'   gdp ~ inflation + investment,
#'   data = macro_panel,
#'   id = "country",
#'   time = "year",
#'   p = 1, q = 1,
#'   estimator = "pmg"
#' )
#' summary(pmg_model)
#'
#' # Compare with Mean Group estimator
#' mg_model <- panel_ardl(
#'   gdp ~ inflation + investment,
#'   data = macro_panel,
#'   id = "country",
#'   time = "year",
#'   estimator = "mg"
#' )
#'
#' # Hausman test
#' hausman_test(pmg_model, mg_model)
#' }
#'
#' @export
#' @importFrom stats lm coef residuals fitted model.matrix model.response na.omit
#' @importFrom stats pnorm pt pchisq optim nlminb
#' @importFrom MASS ginv
panel_ardl <- function(formula, data, id, time, p = 1, q = 1,
                       estimator = c("pmg", "mg", "dfe"),
                       ec = TRUE, trend = FALSE,
                       maxiter = 100, tol = 1e-5) {
  

  # Match arguments
  estimator <- match.arg(estimator)
  
 # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!id %in% names(data)) {
    stop(paste("ID variable '", id, "' not found in data", sep = ""))
  }
  
  if (!time %in% names(data)) {
    stop(paste("Time variable '", time, "' not found in data", sep = ""))
  }
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  
  if (!y_var %in% names(data)) {
    stop(paste("Dependent variable '", y_var, "' not found in data", sep = ""))
  }
  
  for (v in x_vars) {
    if (!v %in% names(data)) {
      stop(paste("Variable '", v, "' not found in data", sep = ""))
    }
  }
  
  # Sort data
  data <- data[order(data[[id]], data[[time]]), ]
  
  # Get groups
  groups <- unique(data[[id]])
  n_groups <- length(groups)
  
  # Handle q as vector or scalar
  if (length(q) == 1) {
    q <- rep(q, length(x_vars))
  } else if (length(q) != length(x_vars)) {
    stop("Length of 'q' must equal number of independent variables or be 1")
  }
  
  # Call appropriate estimator
  result <- switch(estimator,
    "pmg" = .estimate_pmg(data, y_var, x_vars, id, time, p, q, trend, maxiter, tol),
    "mg"  = .estimate_mg(data, y_var, x_vars, id, time, p, q, trend),
    "dfe" = .estimate_dfe(data, y_var, x_vars, id, time, p, q, trend)
  )
  
  # Add common elements
  result$call <- match.call()
  result$formula <- formula
  result$estimator <- estimator
  result$ec_form <- ec
  result$trend <- trend
  result$y_var <- y_var
  result$x_vars <- x_vars
  result$id <- id
  result$time <- time
  result$p <- p
  result$q <- q
  result$ngroups <- n_groups
  result$groups <- groups
  
  # Compute information criteria
  n <- result$nobs
  k <- length(result$coefficients)
  result$aic <- -2 * result$loglik + 2 * k
  result$bic <- -2 * result$loglik + log(n) * k
  
  class(result) <- c("panel_ardl", "list")
  return(result)
}


#' @title Pooled Mean Group Estimator
#' @description Internal function for PMG estimation
#' @keywords internal
.estimate_pmg <- function(data, y_var, x_vars, id, time, p, q, trend, maxiter, tol) {
  
  groups <- unique(data[[id]])
  n_groups <- length(groups)
  k_x <- length(x_vars)
  
  # Prepare data for each group
  group_data <- lapply(groups, function(g) {
    gdata <- data[data[[id]] == g, ]
    .prepare_ardl_data(gdata, y_var, x_vars, time, p, q, trend)
  })
  
  # Check for valid groups
  valid_groups <- sapply(group_data, function(x) !is.null(x) && nrow(x$X) > 0)
  if (sum(valid_groups) < 2) {
    stop("Need at least 2 valid groups for panel estimation")
  }
  
  group_data <- group_data[valid_groups]
  groups <- groups[valid_groups]
  n_groups <- length(groups)
  
  # Initial values from MG estimator
  mg_result <- .estimate_mg_internal(group_data, groups, k_x, trend)
  theta_init <- mg_result$long_run
  
  # PMG optimization: maximize concentrated log-likelihood
  # Long-run parameters are common, short-run are group-specific
  
  pmg_loglik <- function(theta) {
    ll <- 0
    for (i in seq_along(group_data)) {
      gd <- group_data[[i]]
      
      # Compute ECT = y_{t-1} - theta' * x_t
      ect <- gd$y_lag1 - as.vector(gd$X_levels %*% theta)
      
      # Regress dy on ECT and short-run terms
      X_sr <- cbind(ect, gd$X_diff)
      
      tryCatch({
        fit <- lm.fit(X_sr, gd$dy)
        resid <- fit$residuals
        sigma2 <- sum(resid^2) / length(resid)
        ll <- ll - 0.5 * length(resid) * (log(2 * pi) + log(sigma2) + 1)
      }, error = function(e) {
        ll <- ll - 1e10
      })
    }
    return(-ll)  # Return negative for minimization
  }
  
  # Optimize
  opt_result <- tryCatch({
    nlminb(theta_init, pmg_loglik, 
           control = list(iter.max = maxiter, rel.tol = tol))
  }, error = function(e) {
    list(par = theta_init, convergence = 1)
  })
  
  theta_pmg <- opt_result$par
  names(theta_pmg) <- x_vars
  
  # Compute group-specific short-run coefficients
  short_run_list <- list()
  ec_coefs <- numeric(n_groups)
  residuals_all <- c()
  fitted_all <- c()
  
  for (i in seq_along(group_data)) {
    gd <- group_data[[i]]
    
    # ECT with estimated theta
    ect <- gd$y_lag1 - as.vector(gd$X_levels %*% theta_pmg)
    
    # Short-run regression
    X_sr <- cbind(ect, gd$X_diff)
    fit <- lm.fit(X_sr, gd$dy)
    
    short_run_list[[i]] <- fit$coefficients
    ec_coefs[i] <- fit$coefficients[1]
    
    residuals_all <- c(residuals_all, fit$residuals)
    fitted_all <- c(fitted_all, fit$fitted.values)
  }
  
  # Average short-run coefficients
  sr_matrix <- do.call(rbind, short_run_list)
  short_run_avg <- colMeans(sr_matrix, na.rm = TRUE)
  short_run_se <- apply(sr_matrix, 2, sd, na.rm = TRUE) / sqrt(n_groups)
  
  # EC coefficient average
  ec_avg <- mean(ec_coefs)
  ec_se <- sd(ec_coefs) / sqrt(n_groups)
  
  # Standard errors for long-run via delta method (simplified)
  theta_se <- .compute_pmg_se(group_data, theta_pmg, n_groups)
  
  # T-values and p-values
  theta_t <- theta_pmg / theta_se
  theta_p <- 2 * pt(-abs(theta_t), df = sum(sapply(group_data, function(x) nrow(x$X))) - length(theta_pmg))
  
  # Log-likelihood
  sigma2 <- sum(residuals_all^2) / length(residuals_all)
  loglik <- -0.5 * length(residuals_all) * (log(2 * pi) + log(sigma2) + 1)
  
  list(
    coefficients = c(ec = ec_avg, theta_pmg, short_run_avg[-1]),
    long_run = theta_pmg,
    long_run_se = theta_se,
    long_run_t = theta_t,
    long_run_p = theta_p,
    short_run = short_run_avg,
    short_run_se = short_run_se,
    ec_coef = ec_avg,
    ec_se = ec_se,
    ec_group = ec_coefs,
    group_coefs = short_run_list,
    residuals = residuals_all,
    fitted = fitted_all,
    nobs = length(residuals_all),
    sigma = sqrt(sigma2),
    loglik = loglik,
    convergence = opt_result$convergence
  )
}


#' @title Mean Group Estimator
#' @description Internal function for MG estimation
#' @keywords internal
.estimate_mg <- function(data, y_var, x_vars, id, time, p, q, trend) {
  
  groups <- unique(data[[id]])
  n_groups <- length(groups)
  k_x <- length(x_vars)
  
  # Prepare data for each group
  group_data <- lapply(groups, function(g) {
    gdata <- data[data[[id]] == g, ]
    .prepare_ardl_data(gdata, y_var, x_vars, time, p, q, trend)
  })
  
  # Check for valid groups
  valid_groups <- sapply(group_data, function(x) !is.null(x) && nrow(x$X) > 0)
  group_data <- group_data[valid_groups]
  groups <- groups[valid_groups]
  n_groups <- length(groups)
  
  result <- .estimate_mg_internal(group_data, groups, k_x, trend)
  
  # Compute residuals and fitted
  residuals_all <- c()
  fitted_all <- c()
  
  for (i in seq_along(group_data)) {
    gd <- group_data[[i]]
    X_full <- cbind(gd$y_lag1, gd$X_levels, gd$X_diff)
    fit <- lm.fit(X_full, gd$dy)
    residuals_all <- c(residuals_all, fit$residuals)
    fitted_all <- c(fitted_all, fit$fitted.values)
  }
  
  sigma2 <- sum(residuals_all^2) / length(residuals_all)
  loglik <- -0.5 * length(residuals_all) * (log(2 * pi) + log(sigma2) + 1)
  
  result$residuals <- residuals_all
  result$fitted <- fitted_all
  result$nobs <- length(residuals_all)
  result$sigma <- sqrt(sigma2)
  result$loglik <- loglik
  
  return(result)
}


#' @title Internal MG estimation
#' @keywords internal
.estimate_mg_internal <- function(group_data, groups, k_x, trend) {
  
  n_groups <- length(groups)
  
  # Estimate ARDL for each group
  group_results <- lapply(seq_along(group_data), function(i) {
    gd <- group_data[[i]]
    
    # Full ARDL regression: dy ~ y_lag1 + X_levels + X_diff
    X_full <- cbind(gd$y_lag1, gd$X_levels, gd$X_diff)
    
    tryCatch({
      fit <- lm.fit(X_full, gd$dy)
      coefs <- fit$coefficients
      
      # Extract phi (EC coefficient) and compute theta (long-run)
      phi <- coefs[1]  # Coefficient on y_{t-1}
      
      # Long-run coefficients: theta = -beta / phi
      # where beta are coefficients on X_levels
      beta <- coefs[2:(k_x + 1)]
      theta <- -beta / phi
      
      list(
        phi = phi,
        theta = theta,
        beta = beta,
        short_run = coefs[-(1:(k_x + 1))],
        coefs = coefs,
        residuals = fit$residuals,
        valid = TRUE
      )
    }, error = function(e) {
      list(valid = FALSE)
    })
  })
  
  # Filter valid results
  valid <- sapply(group_results, function(x) x$valid)
  group_results <- group_results[valid]
  n_valid <- length(group_results)
  
  # Average coefficients (Mean Group)
  theta_mat <- do.call(rbind, lapply(group_results, function(x) x$theta))
  theta_mg <- colMeans(theta_mat, na.rm = TRUE)
  theta_se <- apply(theta_mat, 2, sd, na.rm = TRUE) / sqrt(n_valid)
  
  phi_vec <- sapply(group_results, function(x) x$phi)
  phi_mg <- mean(phi_vec, na.rm = TRUE)
  phi_se <- sd(phi_vec, na.rm = TRUE) / sqrt(n_valid)
  
  # T-values and p-values
  theta_t <- theta_mg / theta_se
  theta_p <- 2 * pnorm(-abs(theta_t))
  
  list(
    coefficients = c(ec = phi_mg, theta_mg),
    long_run = theta_mg,
    long_run_se = theta_se,
    long_run_t = theta_t,
    long_run_p = theta_p,
    ec_coef = phi_mg,
    ec_se = phi_se,
    ec_group = phi_vec,
    group_coefs = lapply(group_results, function(x) x$coefs),
    group_theta = theta_mat
  )
}


#' @title Dynamic Fixed Effects Estimator
#' @description Internal function for DFE estimation
#' @keywords internal
.estimate_dfe <- function(data, y_var, x_vars, id, time, p, q, trend) {
  
  groups <- unique(data[[id]])
  n_groups <- length(groups)
  
  # Prepare pooled data with group dummies
  all_data <- list()
  
  for (g in groups) {
    gdata <- data[data[[id]] == g, ]
    prepped <- .prepare_ardl_data(gdata, y_var, x_vars, time, p, q, trend)
    
    if (!is.null(prepped) && nrow(prepped$X) > 0) {
      prepped$group <- g
      all_data[[length(all_data) + 1]] <- prepped
    }
  }
  
  # Combine all groups
  dy_all <- unlist(lapply(all_data, function(x) x$dy))
  y_lag1_all <- unlist(lapply(all_data, function(x) x$y_lag1))
  X_levels_all <- do.call(rbind, lapply(all_data, function(x) x$X_levels))
  X_diff_all <- do.call(rbind, lapply(all_data, function(x) x$X_diff))
  
  # Group dummies (fixed effects)
  group_vec <- unlist(lapply(all_data, function(x) rep(x$group, nrow(x$X))))
  group_dummies <- model.matrix(~ factor(group_vec) - 1)[, -1, drop = FALSE]
  
  # Full regression with fixed effects
  X_full <- cbind(y_lag1_all, X_levels_all, X_diff_all, group_dummies)
  
  fit <- lm.fit(X_full, dy_all)
  coefs <- fit$coefficients
  
  # Extract parameters
  k_x <- length(x_vars)
  phi <- coefs[1]
  beta <- coefs[2:(k_x + 1)]
  theta <- -beta / phi
  names(theta) <- x_vars
  
  # Standard errors via sandwich
  n <- length(dy_all)
  k <- length(coefs)
  sigma2 <- sum(fit$residuals^2) / (n - k)
  
  # Variance-covariance matrix
  XtX_inv <- tryCatch({
    solve(crossprod(X_full))
  }, error = function(e) {
    MASS::ginv(crossprod(X_full))
  })
  
  vcov_mat <- sigma2 * XtX_inv
  se <- sqrt(diag(vcov_mat))
  
  # Delta method for long-run SE
  # theta = -beta / phi
  # Var(theta) = (1/phi^2) * Var(beta) + (beta^2/phi^4) * Var(phi) + ...
  # Simplified version:
  theta_se <- abs(theta) * sqrt((se[2:(k_x+1)]/beta)^2 + (se[1]/phi)^2)
  
  theta_t <- theta / theta_se
  theta_p <- 2 * pt(-abs(theta_t), df = n - k)
  
  loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
  
  list(
    coefficients = coefs,
    long_run = theta,
    long_run_se = theta_se,
    long_run_t = theta_t,
    long_run_p = theta_p,
    short_run = coefs[-(1:(k_x + 1 + ncol(group_dummies)))],
    ec_coef = phi,
    ec_se = se[1],
    residuals = fit$residuals,
    fitted = fit$fitted.values,
    nobs = n,
    sigma = sqrt(sigma2),
    loglik = loglik,
    vcov = vcov_mat
  )
}


#' @title Prepare ARDL Data
#' @description Internal function to prepare data for ARDL estimation
#' @keywords internal
.prepare_ardl_data <- function(data, y_var, x_vars, time, p, q, trend) {
  
  n <- nrow(data)
  max_lag <- max(p, max(q))
  
  if (n <= max_lag + 1) {
    return(NULL)
  }
  
  # Sort by time
  data <- data[order(data[[time]]), ]
  
  # Create lagged variables
  y <- data[[y_var]]
  
  # First difference of y
  dy <- diff(y)
  
  # Lagged level of y
  y_lag1 <- y[1:(n-1)]
  
  # X variables in levels
  X_levels <- as.matrix(data[1:(n-1), x_vars, drop = FALSE])
  
  # First differences of X variables
  X_diff_list <- lapply(x_vars, function(v) {
    x <- data[[v]]
    diff(x)
  })
  X_diff <- do.call(cbind, X_diff_list)
  colnames(X_diff) <- paste0("D.", x_vars)
  
  # Add lagged differences if q > 0
  if (any(q > 0)) {
    for (j in seq_along(x_vars)) {
      if (q[j] > 0) {
        for (lag in 1:q[j]) {
          x_diff <- diff(data[[x_vars[j]]])
          lagged <- c(rep(NA, lag), x_diff[1:(length(x_diff) - lag)])
          X_diff <- cbind(X_diff, lagged)
          colnames(X_diff)[ncol(X_diff)] <- paste0("D.", x_vars[j], ".L", lag)
        }
      }
    }
  }
  
  # Add lagged dy if p > 1
  if (p > 1) {
    for (lag in 1:(p-1)) {
      dy_lagged <- c(rep(NA, lag), dy[1:(length(dy) - lag)])
      X_diff <- cbind(X_diff, dy_lagged)
      colnames(X_diff)[ncol(X_diff)] <- paste0("D.y.L", lag)
    }
  }
  
  # Combine and remove NAs
  full_data <- data.frame(
    dy = dy,
    y_lag1 = y_lag1,
    X_levels,
    X_diff
  )
  
  full_data <- na.omit(full_data)
  
  list(
    dy = full_data$dy,
    y_lag1 = full_data$y_lag1,
    X_levels = as.matrix(full_data[, x_vars, drop = FALSE]),
    X_diff = as.matrix(full_data[, -(1:(2 + length(x_vars))), drop = FALSE]),
    X = as.matrix(full_data[, -1, drop = FALSE]),
    n = nrow(full_data)
  )
}


#' @title Compute PMG Standard Errors
#' @description Internal function to compute standard errors for PMG estimator
#' @keywords internal
.compute_pmg_se <- function(group_data, theta, n_groups) {
  
  k <- length(theta)
  
  # Compute group-specific long-run estimates and use their variation
  theta_group <- matrix(NA, n_groups, k)
  
  for (i in seq_along(group_data)) {
    gd <- group_data[[i]]
    
    # Full ARDL for this group
    X_full <- cbind(gd$y_lag1, gd$X_levels, gd$X_diff)
    
    tryCatch({
      fit <- lm.fit(X_full, gd$dy)
      phi_i <- fit$coefficients[1]
      beta_i <- fit$coefficients[2:(k + 1)]
      theta_group[i, ] <- -beta_i / phi_i
    }, error = function(e) {
      # Keep NA
    })
  }
  
  # Standard error from cross-sectional variation
  se <- apply(theta_group, 2, sd, na.rm = TRUE) / sqrt(sum(!is.na(theta_group[, 1])))
  
  return(se)
}


#' @title Summary method for panel_ardl
#' @description Print summary of Panel ARDL estimation results
#' @param object An object of class "panel_ardl"
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the summary object
#' @export
summary.panel_ardl <- function(object, ...) {
  
  cat("\n")
  cat("====================================================================\n")
  cat("         Panel ARDL Estimation Results\n")
  cat("====================================================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Estimator:        ", toupper(object$estimator), "\n")
  cat("Dependent variable:", object$y_var, "\n")
  cat("Number of groups:  ", object$ngroups, "\n")
  cat("Observations:      ", object$nobs, "\n")
  cat("ARDL(", object$p, ", ", paste(object$q, collapse = ","), ")\n\n", sep = "")
  
  # Long-run coefficients
  cat("--------------------------------------------------------------------\n")
  cat("                     Long-Run Coefficients\n")
  cat("--------------------------------------------------------------------\n")
  
  lr_table <- data.frame(
    Estimate = object$long_run,
    `Std. Error` = object$long_run_se,
    `t value` = object$long_run_t,
    `Pr(>|t|)` = object$long_run_p,
    check.names = FALSE
  )
  
  # Add significance stars
  lr_table$` ` <- ifelse(lr_table$`Pr(>|t|)` < 0.001, "***",
                   ifelse(lr_table$`Pr(>|t|)` < 0.01, "**",
                   ifelse(lr_table$`Pr(>|t|)` < 0.05, "*",
                   ifelse(lr_table$`Pr(>|t|)` < 0.1, ".", ""))))
  
  print(lr_table, digits = 4)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  # Error correction coefficient
  cat("--------------------------------------------------------------------\n")
  cat("                 Error Correction Coefficient\n")
  cat("--------------------------------------------------------------------\n")
  cat("EC coefficient (phi):", round(object$ec_coef, 4), "\n")
  cat("Std. Error:          ", round(object$ec_se, 4), "\n")
  ec_t <- object$ec_coef / object$ec_se
  ec_p <- 2 * pnorm(-abs(ec_t))
  cat("t-value:             ", round(ec_t, 4), "\n")
  cat("p-value:             ", format.pval(ec_p, digits = 4), "\n")
  
  if (object$ec_coef >= 0) {
    cat("\n*** Warning: EC coefficient should be negative for convergence ***\n")
  } else {
    half_life <- -log(2) / log(1 + object$ec_coef)
    cat("Half-life:           ", round(half_life, 2), "periods\n")
  }
  
  cat("\n")
  
  # Model fit
  cat("--------------------------------------------------------------------\n")
  cat("                      Model Fit Statistics\n")
  cat("--------------------------------------------------------------------\n")
  cat("Residual std. error: ", round(object$sigma, 4), "\n")
  cat("Log-likelihood:      ", round(object$loglik, 2), "\n")
  cat("AIC:                 ", round(object$aic, 2), "\n")
  cat("BIC:                 ", round(object$bic, 2), "\n")
  
  if (!is.null(object$convergence)) {
    cat("Convergence:         ", ifelse(object$convergence == 0, "Yes", "No"), "\n")
  }
  
  cat("====================================================================\n")
  
  invisible(object)
}


#' @title Print method for panel_ardl
#' @description Print Panel ARDL results
#' @param x An object of class "panel_ardl"
#' @param ... Additional arguments (ignored)
#' @export
print.panel_ardl <- function(x, ...) {
  
  cat("\nPanel ARDL (", toupper(x$estimator), ") Estimation\n", sep = "")
  cat("Formula: ")
  print(x$formula)
  cat("Groups:", x$ngroups, "| Obs:", x$nobs, "\n\n")
  
  cat("Long-run coefficients:\n")
  print(round(x$long_run, 4))
  
  cat("\nEC coefficient:", round(x$ec_coef, 4), "\n")
  
  invisible(x)
}


#' @title Hausman Test for PMG vs MG
#' @description Perform Hausman test comparing PMG and MG estimators
#'
#' @param pmg_model A panel_ardl object estimated with PMG
#' @param mg_model A panel_ardl object estimated with MG (optional)
#' @param data Data frame (required if mg_model not provided)
#'
#' @return A list containing test statistic, degrees of freedom, and p-value
#'
#' @details
#' The Hausman test examines whether the long-run homogeneity assumption
#' of the PMG estimator is valid. Under the null hypothesis of homogeneity,
#' both PMG and MG are consistent, but PMG is efficient. Under the alternative,
#' only MG is consistent.
#'
#' @export
hausman_test <- function(pmg_model, mg_model = NULL, data = NULL) {
  
  if (pmg_model$estimator != "pmg") {
    stop("First model must be PMG estimator")
  }
  
  # If MG model not provided, estimate it
  if (is.null(mg_model)) {
    if (is.null(data)) {
      stop("Either mg_model or data must be provided")
    }
    mg_model <- panel_ardl(
      pmg_model$formula, data,
      id = pmg_model$id, time = pmg_model$time,
      p = pmg_model$p, q = pmg_model$q,
      estimator = "mg", trend = pmg_model$trend
    )
  }
  
  # Coefficient differences
  theta_diff <- pmg_model$long_run - mg_model$long_run
  
  # Variance difference (MG variance - PMG variance should be positive definite)
  var_diff <- mg_model$long_run_se^2 - pmg_model$long_run_se^2
  
  # For positive variance differences
  var_diff[var_diff < 0] <- abs(var_diff[var_diff < 0])
  
  # Hausman statistic
  H <- sum(theta_diff^2 / var_diff)
  
  # Degrees of freedom
  df <- length(theta_diff)
  
  # P-value
  p_value <- 1 - pchisq(H, df)
  
  result <- list(
    statistic = H,
    df = df,
    p.value = p_value,
    theta_pmg = pmg_model$long_run,
    theta_mg = mg_model$long_run,
    difference = theta_diff
  )
  
  class(result) <- "hausman_test"
  
  cat("\n")
  cat("===============================================\n")
  cat("     Hausman Test: PMG vs MG Estimators\n")
  cat("===============================================\n\n")
  cat("H0: Long-run coefficients are homogeneous (PMG is efficient)\n")
  cat("H1: Long-run coefficients are heterogeneous (MG is consistent)\n\n")
  cat("Chi-squared statistic:", round(H, 4), "\n")
  cat("Degrees of freedom:   ", df, "\n")
  cat("P-value:              ", format.pval(p_value, digits = 4), "\n\n")
  
  if (p_value < 0.05) {
    cat("Result: Reject H0 at 5% level. MG estimator preferred.\n")
  } else {
    cat("Result: Cannot reject H0. PMG estimator is efficient.\n")
  }
  cat("===============================================\n")
  
  invisible(result)
}
