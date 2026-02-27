#' @title Rolling and Recursive ARDL (R-ARDL)
#' @description Time-varying ARDL bounds test using rolling or recursive windows.
#'
#' @details
#' The Rolling ARDL (R-ARDL) approach implements time-varying bounds testing
#' to detect structural changes in cointegration relationships. Two methods
#' are available:
#'
#' \itemize{
#'   \item \strong{Rolling}: Fixed-size window that moves through the sample
#'   \item \strong{Recursive}: Expanding window starting from initial observations
#' }
#'
#' This is useful for detecting:
#' \itemize{
#'   \item Time-varying cointegration relationships
#'   \item Structural breaks in long-run parameters
#'   \item Changes in error correction speed
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param method Character. Either "rolling" or "recursive" (default: "rolling")
#' @param window Integer. Window size for rolling method (default: 50)
#' @param min_obs Integer. Minimum observations for recursive method (default: 40)
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param case Integer from 1-5 specifying deterministic components (default: 3)
#' @param parallel Logical. Use parallel processing (default: FALSE)
#' @param ncores Integer. Number of cores for parallel processing
#'
#' @return An object of class "rardl" containing:
#' \itemize{
#'   \item \code{F_stats}: Time series of F-statistics
#'   \item \code{t_stats}: Time series of t-statistics (EC coefficient)
#'   \item \code{cointegration}: Time series of cointegration decisions
#'   \item \code{ec_coefs}: Time-varying error correction coefficients
#'   \item \code{lr_coefs}: Time-varying long-run coefficients
#'   \item \code{dates}: Index/dates for each window
#'   \item \code{breaks}: Detected structural break points
#' }
#'
#' @references
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. Journal of Applied Econometrics, 16(3), 289-326.
#'
#' @examples
#' \dontrun{
#' # Generate example data with structural break
#' n <- 300
#' data <- generate_ts_data(n = n)
#' # Add break at midpoint
#' data$y[151:n] <- data$y[151:n] + 5
#'
#' # Rolling ARDL
#' roll_result <- rardl(y ~ x1 + x2, data = data, method = "rolling", window = 60)
#' plot(roll_result)
#' summary(roll_result)
#'
#' # Recursive ARDL
#' rec_result <- rardl(y ~ x1 + x2, data = data, method = "recursive", min_obs = 50)
#' plot(rec_result)
#' }
#'
#' @export
rardl <- function(formula, data, method = c("rolling", "recursive"),
                  window = 50, min_obs = 40, p = 1, q = 1, case = 3,
                  parallel = FALSE, ncores = 2) {
  
  method <- match.arg(method)
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  k <- length(x_vars)
  
  # Handle q as vector or scalar
  if (length(q) == 1) q <- rep(q, k)
  
  # Get data
  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  n <- length(y)
  
  # Validate window size
  max_lag <- max(p, max(q))
  min_required <- max_lag + k + 5 + case  # Minimum observations needed
  
  if (method == "rolling" && window < min_required) {
    stop(paste("Window size must be at least", min_required, "observations"))
  }
  if (method == "recursive" && min_obs < min_required) {
    stop(paste("Minimum observations must be at least", min_required))
  }
  
  # Determine windows
  if (method == "rolling") {
    n_windows <- n - window + 1
    window_starts <- 1:n_windows
    window_ends <- window:n
  } else {
    n_windows <- n - min_obs + 1
    window_starts <- rep(1, n_windows)
    window_ends <- min_obs:n
  }
  
  # Initialize storage
  F_stats <- numeric(n_windows)
  t_stats <- numeric(n_windows)
  ec_coefs <- numeric(n_windows)
  ec_ses <- numeric(n_windows)
  lr_coefs <- matrix(NA, n_windows, k)
  colnames(lr_coefs) <- x_vars
  coint_decisions <- character(n_windows)
  
  # Get critical values
  cv <- pss_critical_values(k, case)
  
  # Run estimation for each window
  for (w in 1:n_windows) {
    start_idx <- window_starts[w]
    end_idx <- window_ends[w]
    
    # Subset data
    window_data <- data[start_idx:end_idx, , drop = FALSE]
    
    tryCatch({
      # Estimate ARDL model for this window
      result <- .estimate_ardl_window(window_data, y_var, x_vars, p, q, case)
      
      F_stats[w] <- result$F_stat
      t_stats[w] <- result$t_stat
      ec_coefs[w] <- result$ec_coef
      ec_ses[w] <- result$ec_se
      lr_coefs[w, ] <- result$lr_coefs
      
      # Make cointegration decision
      if (F_stats[w] > cv$F_I1["5%"]) {
        coint_decisions[w] <- "COINT"
      } else if (F_stats[w] < cv$F_I0["5%"]) {
        coint_decisions[w] <- "NO_COINT"
      } else {
        coint_decisions[w] <- "INCONC"
      }
      
    }, error = function(e) {
      F_stats[w] <<- NA
      t_stats[w] <<- NA
      ec_coefs[w] <<- NA
      ec_ses[w] <<- NA
      lr_coefs[w, ] <<- NA
      coint_decisions[w] <<- "ERROR"
    })
  }
  
  # Detect structural breaks
  breaks <- .detect_ardl_breaks(F_stats, ec_coefs, lr_coefs)
  
  # Create time index
  if ("Date" %in% names(data)) {
    dates <- data$Date[window_ends]
  } else if (!is.null(rownames(data))) {
    dates <- rownames(data)[window_ends]
  } else {
    dates <- window_ends
  }
  
  # Calculate summary statistics
  summary_stats <- list(
    n_windows = n_windows,
    n_coint = sum(coint_decisions == "COINT", na.rm = TRUE),
    n_no_coint = sum(coint_decisions == "NO_COINT", na.rm = TRUE),
    n_inconc = sum(coint_decisions == "INCONC", na.rm = TRUE),
    F_mean = mean(F_stats, na.rm = TRUE),
    F_sd = stats::sd(F_stats, na.rm = TRUE),
    ec_mean = mean(ec_coefs, na.rm = TRUE),
    ec_sd = stats::sd(ec_coefs, na.rm = TRUE)
  )
  
  result <- list(
    F_stats = F_stats,
    t_stats = t_stats,
    ec_coefs = ec_coefs,
    ec_ses = ec_ses,
    lr_coefs = lr_coefs,
    cointegration = coint_decisions,
    dates = dates,
    breaks = breaks,
    critical_values = cv,
    summary_stats = summary_stats,
    call = match.call(),
    method = method,
    window = if (method == "rolling") window else NULL,
    min_obs = if (method == "recursive") min_obs else NULL,
    n = n,
    k = k,
    p = p,
    q = q,
    case = case
  )
  
  class(result) <- "rardl"
  return(result)
}


#' @title Estimate ARDL for Single Window
#' @keywords internal
.estimate_ardl_window <- function(data, y_var, x_vars, p, q, case) {
  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  n <- length(y)
  k <- length(x_vars)
  
  if (length(q) == 1) q <- rep(q, k)
  
  max_lag <- max(p, max(q))
  valid_idx <- (max_lag + 1):n
  n_valid <- length(valid_idx)
  
  # Build model matrix
  dy <- diff(y)[(max_lag):(n-1)]
  y_lag <- y[valid_idx - 1]
  
  # Lagged differences
  dy_lags <- matrix(NA, n_valid, p)
  for (i in 1:p) {
    dy_lags[, i] <- diff(y)[(max_lag - i):(n - 1 - i)]
  }
  
  # Independent variables
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
  
  x_diffs <- do.call(cbind, x_diff_list)
  
  # Design matrix
  design <- cbind(y_lag, x_levels, dy_lags, x_diffs)
  
  if (case >= 2) design <- cbind(design, intercept = 1)
  if (case >= 4) design <- cbind(design, trend = 1:n_valid)
  
  # Estimate
  model <- stats::lm(dy ~ design - 1)
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  
  # Extract statistics
  ec_coef <- coefs[1]
  ec_se <- sqrt(vcov_mat[1, 1])
  
  n_level <- 1 + k
  beta_h0 <- coefs[1:n_level]
  V_h0 <- vcov_mat[1:n_level, 1:n_level]
  
  F_stat <- as.numeric(t(beta_h0) %*% solve(V_h0) %*% beta_h0 / n_level)
  t_stat <- ec_coef / ec_se
  
  # Long-run coefficients
  if (abs(ec_coef) > 1e-10) {
    lr_coefs <- -coefs[2:(1 + k)] / ec_coef
  } else {
    lr_coefs <- rep(NA, k)
  }
  
  list(
    F_stat = F_stat,
    t_stat = t_stat,
    ec_coef = ec_coef,
    ec_se = ec_se,
    lr_coefs = lr_coefs
  )
}


#' @title Detect Structural Breaks
#' @keywords internal
.detect_ardl_breaks <- function(F_stats, ec_coefs, lr_coefs) {
  n <- length(F_stats)
  
  # Remove NAs for analysis
  valid <- !is.na(F_stats) & !is.na(ec_coefs)
  
  breaks <- list()
  
  # Detect breaks in F-statistic series (CUSUM-like)
  if (sum(valid) > 20) {
    F_valid <- F_stats[valid]
    F_mean <- mean(F_valid)
    F_sd <- stats::sd(F_valid)
    
    # Find significant deviations
    z_scores <- (F_valid - F_mean) / F_sd
    break_candidates <- which(abs(z_scores) > 2)
    
    if (length(break_candidates) > 0) {
      # Cluster nearby breaks
      breaks$F_stat <- .cluster_breaks(break_candidates)
    }
  }
  
  # Detect breaks in EC coefficient
  if (sum(valid) > 20) {
    ec_valid <- ec_coefs[valid]
    
    # Check for sign changes
    sign_changes <- which(diff(sign(ec_valid)) != 0)
    if (length(sign_changes) > 0) {
      breaks$ec_sign <- sign_changes
    }
    
    # Check for magnitude changes
    ec_mean <- mean(ec_valid)
    ec_sd <- stats::sd(ec_valid)
    z_scores <- (ec_valid - ec_mean) / ec_sd
    break_candidates <- which(abs(z_scores) > 2)
    
    if (length(break_candidates) > 0) {
      breaks$ec_magnitude <- .cluster_breaks(break_candidates)
    }
  }
  
  # Detect breaks in long-run coefficients
  if (!is.null(lr_coefs) && sum(valid) > 20) {
    for (j in 1:ncol(lr_coefs)) {
      lr_j <- lr_coefs[valid, j]
      if (sum(!is.na(lr_j)) > 20) {
        lr_j <- lr_j[!is.na(lr_j)]
        lr_mean <- mean(lr_j)
        lr_sd <- stats::sd(lr_j)
        
        if (lr_sd > 0) {
          z_scores <- (lr_j - lr_mean) / lr_sd
          break_candidates <- which(abs(z_scores) > 2)
          
          if (length(break_candidates) > 0) {
            breaks[[paste0("lr_", colnames(lr_coefs)[j])]] <- .cluster_breaks(break_candidates)
          }
        }
      }
    }
  }
  
  return(breaks)
}


#' @title Cluster Nearby Breaks
#' @keywords internal
.cluster_breaks <- function(breaks, min_gap = 5) {
  if (length(breaks) <= 1) return(breaks)
  
  clusters <- list()
  current_cluster <- breaks[1]
  
  for (i in 2:length(breaks)) {
    if (breaks[i] - breaks[i-1] <= min_gap) {
      current_cluster <- c(current_cluster, breaks[i])
    } else {
      clusters <- c(clusters, list(current_cluster))
      current_cluster <- breaks[i]
    }
  }
  clusters <- c(clusters, list(current_cluster))
  
  # Return midpoint of each cluster
  sapply(clusters, function(x) round(mean(x)))
}


#' @export
print.rardl <- function(x, ...) {
  cat("\n")
  cat("Rolling/Recursive ARDL Bounds Test\n")
  cat(paste(rep("=", 45), collapse = ""), "\n\n")
  
  cat("Method:", x$method, "\n")
  if (x$method == "rolling") {
    cat("Window size:", x$window, "\n")
  } else {
    cat("Minimum obs:", x$min_obs, "\n")
  }
  cat("Number of windows:", x$summary_stats$n_windows, "\n\n")
  
  cat("Cointegration Results:\n")
  cat("  Cointegrated:   ", x$summary_stats$n_coint, 
      sprintf("(%.1f%%)", 100 * x$summary_stats$n_coint / x$summary_stats$n_windows), "\n")
  cat("  No cointegration:", x$summary_stats$n_no_coint,
      sprintf("(%.1f%%)", 100 * x$summary_stats$n_no_coint / x$summary_stats$n_windows), "\n")
  cat("  Inconclusive:   ", x$summary_stats$n_inconc,
      sprintf("(%.1f%%)", 100 * x$summary_stats$n_inconc / x$summary_stats$n_windows), "\n")
  
  if (length(x$breaks) > 0) {
    cat("\nDetected Breaks:\n")
    for (b in names(x$breaks)) {
      cat("  ", b, ": positions", paste(x$breaks[[b]], collapse = ", "), "\n")
    }
  }
  
  invisible(x)
}


#' @export
summary.rardl <- function(object, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("       Rolling/Recursive ARDL Analysis Results\n")
  cat("=======================================================\n\n")
  
  cat("Model Specification:\n")
  cat("  Method:", object$method, "\n")
  if (object$method == "rolling") {
    cat("  Window size:", object$window, "\n")
  } else {
    cat("  Minimum observations:", object$min_obs, "\n")
  }
  cat("  Case:", object$case, "\n")
  cat("  Lags: p =", object$p, ", q =", paste(object$q, collapse = ","), "\n")
  cat("  Total observations:", object$n, "\n")
  cat("  Number of windows:", object$summary_stats$n_windows, "\n\n")
  
  cat("F-Statistic Summary:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  Mean:   %8.4f\n", object$summary_stats$F_mean))
  cat(sprintf("  Std:    %8.4f\n", object$summary_stats$F_sd))
  cat(sprintf("  Min:    %8.4f\n", min(object$F_stats, na.rm = TRUE)))
  cat(sprintf("  Max:    %8.4f\n", max(object$F_stats, na.rm = TRUE)))
  
  cat("\nError Correction Coefficient Summary:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  Mean:   %8.4f\n", object$summary_stats$ec_mean))
  cat(sprintf("  Std:    %8.4f\n", object$summary_stats$ec_sd))
  cat(sprintf("  Min:    %8.4f\n", min(object$ec_coefs, na.rm = TRUE)))
  cat(sprintf("  Max:    %8.4f\n", max(object$ec_coefs, na.rm = TRUE)))
  
  cat("\nCointegration Decisions:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  Cointegrated:    %4d (%5.1f%%)\n", 
              object$summary_stats$n_coint,
              100 * object$summary_stats$n_coint / object$summary_stats$n_windows))
  cat(sprintf("  No cointegration: %4d (%5.1f%%)\n",
              object$summary_stats$n_no_coint,
              100 * object$summary_stats$n_no_coint / object$summary_stats$n_windows))
  cat(sprintf("  Inconclusive:    %4d (%5.1f%%)\n",
              object$summary_stats$n_inconc,
              100 * object$summary_stats$n_inconc / object$summary_stats$n_windows))
  
  cat("\nLong-Run Coefficients Summary:\n")
  cat("-------------------------------------------------------\n")
  lr_summary <- apply(object$lr_coefs, 2, function(x) {
    c(Mean = mean(x, na.rm = TRUE),
      SD = stats::sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE))
  })
  print(round(t(lr_summary), 4))
  
  if (length(object$breaks) > 0) {
    cat("\nDetected Structural Breaks:\n")
    cat("-------------------------------------------------------\n")
    for (b in names(object$breaks)) {
      cat("  ", b, ": positions", paste(object$breaks[[b]], collapse = ", "), "\n")
    }
  }
  
  cat("\n=======================================================\n\n")
  
  invisible(object)
}


#' @export
plot.rardl <- function(x, type = c("F", "ec", "lr", "all"), ...) {
  type <- match.arg(type)
  
  n <- length(x$F_stats)
  idx <- 1:n
  
  if (type == "all") {
    graphics::par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  } else {
    graphics::par(mfrow = c(1, 1))
  }
  
  if (type %in% c("F", "all")) {
    # Plot F-statistics
    graphics::plot(idx, x$F_stats, type = "l", col = "steelblue", lwd = 1.5,
                  xlab = "Window", ylab = "F-statistic",
                  main = "Rolling F-Statistic (Bounds Test)")
    
    # Add critical value bands
    graphics::abline(h = x$critical_values$F_I1["5%"], col = "red", lty = 2)
    graphics::abline(h = x$critical_values$F_I0["5%"], col = "red", lty = 2)
    
    # Shade regions
    graphics::polygon(c(idx, rev(idx)), 
                     c(rep(x$critical_values$F_I0["5%"], n), 
                       rep(x$critical_values$F_I1["5%"], n)),
                     col = grDevices::rgb(1, 0, 0, 0.1), border = NA)
    
    graphics::legend("topright", 
                    legend = c("F-stat", "I(0)/I(1) bounds"),
                    col = c("steelblue", "red"), lty = c(1, 2), cex = 0.8)
    
    # Mark breaks
    if (!is.null(x$breaks$F_stat)) {
      graphics::abline(v = x$breaks$F_stat, col = "darkgreen", lty = 3, lwd = 2)
    }
  }
  
  if (type %in% c("ec", "all")) {
    # Plot EC coefficients
    graphics::plot(idx, x$ec_coefs, type = "l", col = "darkblue", lwd = 1.5,
                  xlab = "Window", ylab = "EC Coefficient",
                  main = "Time-Varying Error Correction Coefficient")
    graphics::abline(h = 0, col = "gray", lty = 2)
    
    # Add confidence bands
    upper <- x$ec_coefs + 1.96 * x$ec_ses
    lower <- x$ec_coefs - 1.96 * x$ec_ses
    graphics::polygon(c(idx, rev(idx)), c(upper, rev(lower)),
                     col = grDevices::rgb(0, 0, 1, 0.2), border = NA)
    
    graphics::lines(idx, x$ec_coefs, col = "darkblue", lwd = 1.5)
  }
  
  if (type %in% c("lr", "all")) {
    # Plot long-run coefficients
    n_coefs <- ncol(x$lr_coefs)
    colors <- grDevices::rainbow(n_coefs)
    
    graphics::matplot(idx, x$lr_coefs, type = "l", lty = 1, col = colors, lwd = 1.5,
                     xlab = "Window", ylab = "Long-Run Coefficient",
                     main = "Time-Varying Long-Run Coefficients")
    graphics::abline(h = 0, col = "gray", lty = 2)
    graphics::legend("topright", legend = colnames(x$lr_coefs), 
                    col = colors, lty = 1, cex = 0.8)
  }
  
  invisible(x)
}
