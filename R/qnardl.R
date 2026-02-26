#' @title Quantile Nonlinear ARDL (QNARDL)
#' @description Estimate Quantile Nonlinear ARDL models combining distributional
#'   and asymmetric effects
#'
#' @details
#' This function implements the QNARDL model which extends the NARDL framework
#' of Shin, Yu & Greenwood-Nimmo (2014) to a quantile regression setting,
#' combining it with the QARDL approach of Cho, Kim & Shin (2015).
#'
#' The model allows for:
#' \itemize{
#'   \item \strong{Asymmetric effects}: Positive and negative changes in X
#'         can have different impacts on Y
#'   \item \strong{Distributional heterogeneity}: Effects can vary across
#'         quantiles of the conditional distribution
#' }
#'
#' The partial sum decomposition separates each regressor into positive and
#' negative components:
#' \deqn{x^+_t = \sum_{j=1}^{t} \max(\Delta x_j, 0)}
#' \deqn{x^-_t = \sum_{j=1}^{t} \min(\Delta x_j, 0)}
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param tau Numeric vector of quantiles to estimate (default: c(0.25, 0.5, 0.75))
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param decompose Character vector. Variables to decompose into +/- components.
#'   Default is all x variables.
#' @param trend Logical. Include time trend (default: FALSE)
#'
#' @return An object of class "qnardl" containing:
#' \itemize{
#'   \item \code{coefficients}: Estimated coefficients for each quantile
#'   \item \code{long_run_pos}: Positive long-run coefficients by quantile
#'   \item \code{long_run_neg}: Negative long-run coefficients by quantile
#'   \item \code{short_run}: Short-run coefficients by quantile
#'   \item \code{ec_coef}: Error correction coefficients by quantile
#'   \item \code{asymmetry_test}: Wald test for long-run asymmetry
#'   \item \code{tau}: Quantiles estimated
#'   \item \code{residuals}: Model residuals
#'   \item \code{fitted}: Fitted values
#' }
#'
#' @references
#' Cho, J. S., Kim, T. H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. Journal of Econometrics,
#' 188(1), 281-300.
#'
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In Festschrift in Honor of Peter Schmidt (pp. 281-314). Springer.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(oil_data)
#'
#' # Estimate QNARDL model
#' qnardl_model <- qnardl(
#'   gasoline ~ oil_price + exchange_rate,
#'   data = oil_data,
#'   tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
#'   p = 2, q = 2
#' )
#' summary(qnardl_model)
#' plot(qnardl_model)
#'
#' # Test for asymmetry
#' asymmetry_test(qnardl_model, var = "oil_price")
#' }
#'
#' @export
#' @importFrom quantreg rq
#' @importFrom stats pchisq model.matrix
qnardl <- function(formula, data, tau = c(0.25, 0.5, 0.75),
                   p = 1, q = 1, decompose = NULL, trend = FALSE) {
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  k <- length(x_vars)
  
  # Default: decompose all variables
  if (is.null(decompose)) {
    decompose <- x_vars
  }
  
  # Validate decompose variables
  if (!all(decompose %in% x_vars)) {
    stop("All 'decompose' variables must be in the formula")
  }
  
  # Handle q as vector or scalar
  if (length(q) == 1) {
    q <- rep(q, k)
  }
  
  # Prepare data with partial sum decomposition
  qnardl_data <- .prepare_qnardl_data(data, y_var, x_vars, decompose, p, q, trend)
  
  if (is.null(qnardl_data) || nrow(qnardl_data$model_data) < 10) {
    stop("Insufficient observations for specified model")
  }
  
  # Get variable names
  y_col <- "dy"
  x_cols <- setdiff(names(qnardl_data$model_data), "dy")
  
  # Estimate quantile regressions for each tau
  results_by_tau <- lapply(tau, function(t) {
    .estimate_qnardl_tau(qnardl_data$model_data, y_col, x_cols, t)
  })
  names(results_by_tau) <- paste0("tau_", tau)
  
  # Extract coefficients matrix
  coef_matrix <- do.call(cbind, lapply(results_by_tau, function(r) r$coefficients))
  colnames(coef_matrix) <- paste0("tau=", tau)
  
  # Extract long-run coefficients
  long_run_pos <- .extract_long_run(results_by_tau, qnardl_data$pos_vars, tau)
  long_run_neg <- .extract_long_run(results_by_tau, qnardl_data$neg_vars, tau)
  
  # EC coefficients
  ec_coefs <- sapply(results_by_tau, function(r) {
    idx <- grep("y_lag1", names(r$coefficients))
    if (length(idx) > 0) r$coefficients[idx] else NA
  })
  
  # Asymmetry tests
  asymmetry_tests <- lapply(decompose, function(v) {
    .test_asymmetry(results_by_tau, v, tau)
  })
  names(asymmetry_tests) <- decompose
  
  result <- list(
    coefficients = coef_matrix,
    long_run_pos = long_run_pos,
    long_run_neg = long_run_neg,
    ec_coef = ec_coefs,
    asymmetry_test = asymmetry_tests,
    tau = tau,
    results_by_tau = results_by_tau,
    decompose = decompose,
    model_data = qnardl_data$model_data,
    pos_vars = qnardl_data$pos_vars,
    neg_vars = qnardl_data$neg_vars,
    call = match.call(),
    formula = formula,
    y_var = y_var,
    x_vars = x_vars,
    p = p,
    q = q
  )
  
  class(result) <- c("qnardl", "list")
  return(result)
}


#' @title Prepare QNARDL Data with Partial Sum Decomposition
#' @keywords internal
.prepare_qnardl_data <- function(data, y_var, x_vars, decompose, p, q, trend) {
  
  n <- nrow(data)
  max_lag <- max(p, max(q))
  
  if (n <= max_lag + 5) {
    return(NULL)
  }
  
  # Dependent variable
  y <- data[[y_var]]
  dy <- diff(y)
  y_lag1 <- y[-length(y)]
  
  # Initialize result data frame
  result_df <- data.frame(
    dy = dy[-1],
    y_lag1 = y_lag1[-1]
  )
  
  # Track positive and negative variable names
  pos_vars <- c()
  neg_vars <- c()
  
  # Process each x variable
  for (j in seq_along(x_vars)) {
    v <- x_vars[j]
    x <- data[[v]]
    dx <- diff(x)
    
    if (v %in% decompose) {
      # Partial sum decomposition
      dx_pos <- pmax(dx, 0)
      dx_neg <- pmin(dx, 0)
      
      # Cumulative sums (partial sums)
      x_pos <- cumsum(dx_pos)
      x_neg <- cumsum(dx_neg)
      
      # Levels (lagged)
      result_df[[paste0(v, "_pos")]] <- x_pos[-1]
      result_df[[paste0(v, "_neg")]] <- x_neg[-1]
      
      pos_vars <- c(pos_vars, paste0(v, "_pos"))
      neg_vars <- c(neg_vars, paste0(v, "_neg"))
      
      # Differences
      result_df[[paste0("d_", v, "_pos")]] <- dx_pos[-1]
      result_df[[paste0("d_", v, "_neg")]] <- dx_neg[-1]
      
      # Lagged differences
      if (q[j] > 0) {
        for (lag in 1:q[j]) {
          result_df[[paste0("d_", v, "_pos_L", lag)]] <- c(rep(NA, lag), dx_pos[1:(length(dx_pos) - lag - 1)])
          result_df[[paste0("d_", v, "_neg_L", lag)]] <- c(rep(NA, lag), dx_neg[1:(length(dx_neg) - lag - 1)])
        }
      }
      
    } else {
      # No decomposition - standard ARDL
      result_df[[v]] <- x[-c(1, length(x))]
      result_df[[paste0("d_", v)]] <- dx[-1]
      
      if (q[j] > 0) {
        for (lag in 1:q[j]) {
          result_df[[paste0("d_", v, "_L", lag)]] <- c(rep(NA, lag), dx[1:(length(dx) - lag - 1)])
        }
      }
    }
  }
  
  # Lagged dy
  if (p > 1) {
    for (lag in 1:(p-1)) {
      result_df[[paste0("dy_L", lag)]] <- c(rep(NA, lag), dy[1:(length(dy) - lag - 1)])
    }
  }
  
  # Trend
  if (trend) {
    result_df$trend <- 1:nrow(result_df)
  }
  
  # Remove NAs
  result_df <- na.omit(result_df)
  
  list(
    model_data = result_df,
    pos_vars = pos_vars,
    neg_vars = neg_vars
  )
}


#' @title Estimate QNARDL for Single Quantile
#' @keywords internal
.estimate_qnardl_tau <- function(model_data, y_col, x_cols, tau) {
  
  # Build formula
  formula_str <- paste(y_col, "~", paste(x_cols, collapse = " + "))
  
  # Estimate quantile regression
  fit <- quantreg::rq(
    as.formula(formula_str),
    tau = tau,
    data = model_data,
    method = "br"
  )
  
  # Get summary for inference
  fit_summary <- tryCatch({
    summary(fit, se = "boot", R = 200)
  }, error = function(e) {
    summary(fit, se = "nid")
  })
  
  coef_table <- fit_summary$coefficients
  
  list(
    coefficients = coef_table[, 1],
    se = coef_table[, 2],
    t_values = coef_table[, 3],
    p_values = coef_table[, 4],
    residuals = fit$residuals,
    fitted = fit$fitted.values,
    tau = tau,
    model = fit
  )
}


#' @title Extract Long-Run Coefficients
#' @keywords internal
.extract_long_run <- function(results_by_tau, level_vars, tau) {
  
  lr_matrix <- matrix(NA, nrow = length(level_vars), ncol = length(tau))
  rownames(lr_matrix) <- level_vars
  colnames(lr_matrix) <- paste0("tau=", tau)
  
  for (i in seq_along(results_by_tau)) {
    res <- results_by_tau[[i]]
    phi <- res$coefficients["y_lag1"]
    
    for (j in seq_along(level_vars)) {
      v <- level_vars[j]
      if (v %in% names(res$coefficients)) {
        beta <- res$coefficients[v]
        # Long-run: theta = -beta / phi
        lr_matrix[j, i] <- -beta / phi
      }
    }
  }
  
  return(lr_matrix)
}


#' @title Test for Long-Run Asymmetry
#' @keywords internal
.test_asymmetry <- function(results_by_tau, var, tau) {
  
  pos_var <- paste0(var, "_pos")
  neg_var <- paste0(var, "_neg")
  
  test_results <- data.frame(
    tau = tau,
    theta_pos = NA,
    theta_neg = NA,
    diff = NA,
    wald_stat = NA,
    p_value = NA
  )
  
  for (i in seq_along(results_by_tau)) {
    res <- results_by_tau[[i]]
    phi <- res$coefficients["y_lag1"]
    
    if (pos_var %in% names(res$coefficients) && neg_var %in% names(res$coefficients)) {
      beta_pos <- res$coefficients[pos_var]
      beta_neg <- res$coefficients[neg_var]
      
      theta_pos <- -beta_pos / phi
      theta_neg <- -beta_neg / phi
      
      test_results$theta_pos[i] <- theta_pos
      test_results$theta_neg[i] <- theta_neg
      test_results$diff[i] <- theta_pos - theta_neg
      
      # Simplified Wald test (delta method approximation)
      se_pos <- res$se[pos_var]
      se_neg <- res$se[neg_var]
      se_phi <- res$se["y_lag1"]
      
      # Variance of difference (simplified)
      var_theta_pos <- (se_pos / phi)^2 + (beta_pos * se_phi / phi^2)^2
      var_theta_neg <- (se_neg / phi)^2 + (beta_neg * se_phi / phi^2)^2
      var_diff <- var_theta_pos + var_theta_neg
      
      wald <- (theta_pos - theta_neg)^2 / var_diff
      test_results$wald_stat[i] <- wald
      test_results$p_value[i] <- 1 - pchisq(wald, df = 1)
    }
  }
  
  return(test_results)
}


#' @title Test Asymmetry in QNARDL Model
#' @description Perform Wald test for long-run asymmetry
#'
#' @param object A qnardl object
#' @param var Variable name to test asymmetry for
#'
#' @return Data frame with test results by quantile
#' @export
asymmetry_test <- function(object, var) {
  
  if (!inherits(object, "qnardl")) {
    stop("Object must be of class 'qnardl'")
  }
  
  if (!var %in% object$decompose) {
    stop(paste("Variable", var, "was not decomposed. Cannot test asymmetry."))
  }
  
  result <- object$asymmetry_test[[var]]
  
  cat("\n")
  cat("============================================================\n")
  cat("     Wald Test for Long-Run Asymmetry:", var, "\n")
  cat("============================================================\n\n")
  cat("H0: theta+ = theta- (symmetric long-run effect)\n")
  cat("H1: theta+ != theta- (asymmetric long-run effect)\n\n")
  
  print(result, digits = 4, row.names = FALSE)
  
  cat("\n")
  sig_tau <- result$tau[result$p_value < 0.05]
  if (length(sig_tau) > 0) {
    cat("Significant asymmetry at 5% level for tau =", 
        paste(sig_tau, collapse = ", "), "\n")
  } else {
    cat("No significant asymmetry detected at 5% level\n")
  }
  cat("============================================================\n")
  
  invisible(result)
}


#' @title Summary method for qnardl
#' @export
summary.qnardl <- function(object, ...) {
  
  cat("\n")
  cat("====================================================================\n")
  cat("     Quantile Nonlinear ARDL (QNARDL) Estimation Results\n")
  cat("====================================================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Model: QNARDL(", object$p, ", ", paste(object$q, collapse = ", "), ")\n", sep = "")
  cat("Quantiles:", paste(object$tau, collapse = ", "), "\n")
  cat("Decomposed variables:", paste(object$decompose, collapse = ", "), "\n")
  cat("Observations:", nrow(object$model_data), "\n\n")
  
  # Long-run coefficients - Positive
  cat("--------------------------------------------------------------------\n")
  cat("     Long-Run Coefficients: Positive Changes (theta+)\n")
  cat("--------------------------------------------------------------------\n")
  print(round(object$long_run_pos, 4))
  cat("\n")
  
  # Long-run coefficients - Negative
  cat("--------------------------------------------------------------------\n")
  cat("     Long-Run Coefficients: Negative Changes (theta-)\n")
  cat("--------------------------------------------------------------------\n")
  print(round(object$long_run_neg, 4))
  cat("\n")
  
  # EC coefficients
  cat("--------------------------------------------------------------------\n")
  cat("     Error Correction Coefficients by Quantile\n")
  cat("--------------------------------------------------------------------\n")
  ec_df <- data.frame(
    tau = object$tau,
    phi = round(object$ec_coef, 4)
  )
  print(ec_df, row.names = FALSE)
  cat("\n")
  
  # Asymmetry tests
  cat("--------------------------------------------------------------------\n")
  cat("     Long-Run Asymmetry Tests (Wald)\n")
  cat("--------------------------------------------------------------------\n")
  for (v in object$decompose) {
    cat("\nVariable:", v, "\n")
    test_res <- object$asymmetry_test[[v]]
    test_summary <- data.frame(
      tau = test_res$tau,
      `theta+` = round(test_res$theta_pos, 4),
      `theta-` = round(test_res$theta_neg, 4),
      Wald = round(test_res$wald_stat, 3),
      `p-value` = round(test_res$p_value, 4),
      check.names = FALSE
    )
    test_summary$` ` <- ifelse(test_summary$`p-value` < 0.01, "***",
                         ifelse(test_summary$`p-value` < 0.05, "**",
                         ifelse(test_summary$`p-value` < 0.1, "*", "")))
    print(test_summary, row.names = FALSE)
  }
  cat("\n---\n")
  cat("Signif. codes:  0 '***' 0.01 '**' 0.05 '*' 0.1\n")
  cat("====================================================================\n")
  
  invisible(object)
}


#' @title Print method for qnardl
#' @export
print.qnardl <- function(x, ...) {
  
  cat("\nQuantile Nonlinear ARDL (QNARDL)\n")
  cat("Formula: ")
  print(x$formula)
  cat("Quantiles:", paste(x$tau, collapse = ", "), "\n")
  cat("Decomposed:", paste(x$decompose, collapse = ", "), "\n\n")
  
  cat("EC coefficients by quantile:\n")
  print(round(x$ec_coef, 4))
  
  invisible(x)
}


#' @title Plot QNARDL Coefficients Across Quantiles
#' @description Visualize how coefficients change across quantiles
#'
#' @param x A qnardl object
#' @param var Variable to plot (default: first decomposed variable)
#' @param type "long_run" or "asymmetry"
#' @param ... Additional arguments
#'
#' @export
plot.qnardl <- function(x, var = NULL, type = "long_run", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  
  if (is.null(var)) {
    var <- x$decompose[1]
  }
  
  pos_var <- paste0(var, "_pos")
  neg_var <- paste0(var, "_neg")
  
  # Extract coefficients
  theta_pos <- x$long_run_pos[pos_var, ]
  theta_neg <- x$long_run_neg[neg_var, ]
  
  df <- data.frame(
    tau = rep(x$tau, 2),
    theta = c(theta_pos, theta_neg),
    type = rep(c("Positive", "Negative"), each = length(x$tau))
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = tau, y = theta, color = type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c("Positive" = "#2E86AB", "Negative" = "#E94F37"),
      name = "Change Direction"
    ) +
    ggplot2::labs(
      title = paste("Long-Run Coefficients Across Quantiles:", var),
      subtitle = "Asymmetric effects of positive vs negative changes",
      x = expression(tau~" (Quantile)"),
      y = expression(theta~" (Long-run coefficient)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  print(p)
  invisible(p)
}


#' @title Dynamic Multipliers for QNARDL
#' @description Compute and plot cumulative dynamic multipliers
#'
#' @param object A qnardl object
#' @param var Variable name
#' @param tau Quantile to compute multipliers for
#' @param horizon Number of periods (default: 20)
#'
#' @return Data frame with multipliers
#' @export
dynamic_multipliers <- function(object, var, tau = 0.5, horizon = 20) {
  
  if (!inherits(object, "qnardl")) {
    stop("Object must be of class 'qnardl'")
  }
  
  tau_name <- paste0("tau_", tau)
  if (!tau_name %in% names(object$results_by_tau)) {
    stop(paste("Quantile", tau, "not estimated"))
  }
  
  res <- object$results_by_tau[[tau_name]]
  coefs <- res$coefficients
  
  # Get phi (EC coefficient)
  phi <- coefs["y_lag1"]
  
  # Get short-run coefficients
  pos_var <- paste0(var, "_pos")
  neg_var <- paste0(var, "_neg")
  d_pos_var <- paste0("d_", var, "_pos")
  d_neg_var <- paste0("d_", var, "_neg")
  
  # Initial impacts
  delta_pos_0 <- ifelse(d_pos_var %in% names(coefs), coefs[d_pos_var], 0)
  delta_neg_0 <- ifelse(d_neg_var %in% names(coefs), coefs[d_neg_var], 0)
  
  # Long-run multipliers
  beta_pos <- ifelse(pos_var %in% names(coefs), coefs[pos_var], 0)
  beta_neg <- ifelse(neg_var %in% names(coefs), coefs[neg_var], 0)
  theta_pos <- -beta_pos / phi
  theta_neg <- -beta_neg / phi
  
  # Cumulative multipliers
  h <- 0:horizon
  rho <- 1 + phi  # Persistence parameter
  
  m_pos <- theta_pos * (1 - rho^(h + 1)) / (1 - rho)
  m_neg <- theta_neg * (1 - rho^(h + 1)) / (1 - rho)
  
  df <- data.frame(
    horizon = rep(h, 2),
    multiplier = c(m_pos, m_neg),
    type = rep(c("Positive", "Negative"), each = length(h))
  )
  
  # Plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon, y = multiplier, color = type)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = c(theta_pos, theta_neg), 
                          linetype = "dashed", alpha = 0.5) +
      ggplot2::scale_color_manual(
        values = c("Positive" = "#2E86AB", "Negative" = "#E94F37")
      ) +
      ggplot2::labs(
        title = paste("Cumulative Dynamic Multipliers:", var),
        subtitle = paste("Quantile =", tau),
        x = "Horizon",
        y = "Cumulative Multiplier"
      ) +
      ggplot2::theme_minimal()
    
    print(p)
  }
  
  invisible(df)
}
