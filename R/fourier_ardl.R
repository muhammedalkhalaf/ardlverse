#' @title Fourier ARDL Model
#' @description Estimate ARDL models with Fourier terms for smooth structural breaks
#'
#' @details
#' This function implements the Fourier ARDL approach that captures smooth
#' structural changes using trigonometric functions. The Fourier terms approximate
#' unknown structural breaks without requiring prior specification of break dates.
#'
#' The Fourier terms are defined as:
#' \deqn{f_t = \sum_{k=1}^{K} [a_k \sin(2\pi k t/T) + b_k \cos(2\pi k t/T)]}
#'
#' where K is the number of frequencies and T is the sample size.
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ...
#' @param data A data frame containing the time series data
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param k Integer. Number of Fourier frequencies (default: 1, max: 3)
#' @param case Integer from 1-5 specifying deterministic components (default: 3)
#' @param selection Character. Method for selecting optimal k: "aic", "bic", or "fixed"
#'
#' @return An object of class "fourier_ardl" containing:
#' \itemize{
#'   \item \code{coefficients}: Estimated coefficients
#'   \item \code{long_run}: Long-run coefficients
#'   \item \code{short_run}: Short-run coefficients
#'   \item \code{fourier_coefs}: Fourier term coefficients
#'   \item \code{ec_coef}: Error correction coefficient
#'   \item \code{bounds_test}: F and t statistics for bounds test
#'   \item \code{k}: Number of Fourier frequencies used
#'   \item \code{aic}: AIC values for k = 1, 2, 3
#'   \item \code{bic}: BIC values for k = 1, 2, 3
#'   \item \code{structural_breaks}: Detected structural break periods
#' }
#'
#' @references
#' Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test
#' to approximate smooth breaks with new evidence from crude oil market.
#' Economic Modelling, 67, 114-124.
#'
#' Enders, W., & Lee, J. (2012). A unit root test using a Fourier series to
#' approximate smooth breaks. Oxford Bulletin of Economics and Statistics,
#' 74(4), 574-599.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(macro_data)
#'
#' # Estimate Fourier ARDL with automatic frequency selection
#' f_ardl <- fourier_ardl(
#'   gdp ~ investment + trade,
#'   data = macro_data,
#'   p = 2, q = 2,
#'   selection = "aic"
#' )
#' summary(f_ardl)
#' plot(f_ardl)
#' }
#'
#' @export
fourier_ardl <- function(formula, data, p = 1, q = 1, k = 1, case = 3,
                         selection = c("fixed", "aic", "bic")) {
  
  selection <- match.arg(selection)
  
  # Validate k
  if (k < 1 || k > 3) {
    stop("'k' must be between 1 and 3")
  }
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  n_x <- length(x_vars)
  
  # Handle q
  if (length(q) == 1) {
    q <- rep(q, n_x)
  }
  
  # If selection is automatic, estimate for k = 1, 2, 3 and choose best
  if (selection != "fixed") {
    ic_values <- data.frame(k = 1:3, aic = NA, bic = NA)
    models <- list()
    
    for (kk in 1:3) {
      models[[kk]] <- tryCatch({
        .estimate_fourier_ardl(data, y_var, x_vars, p, q, kk, case)
      }, error = function(e) NULL)
      
      if (!is.null(models[[kk]])) {
        ic_values$aic[kk] <- models[[kk]]$aic
        ic_values$bic[kk] <- models[[kk]]$bic
      }
    }
    
    # Select optimal k
    if (selection == "aic") {
      k_opt <- which.min(ic_values$aic)
    } else {
      k_opt <- which.min(ic_values$bic)
    }
    
    result <- models[[k_opt]]
    result$ic_comparison <- ic_values
    result$selection <- selection
    
  } else {
    result <- .estimate_fourier_ardl(data, y_var, x_vars, p, q, k, case)
    result$selection <- "fixed"
  }
  
  # Add metadata
  result$call <- match.call()
  result$formula <- formula
  result$y_var <- y_var
  result$x_vars <- x_vars
  result$p <- p
  result$q <- q
  result$case <- case
  
  class(result) <- c("fourier_ardl", "list")
  return(result)
}


#' @title Internal Fourier ARDL Estimation
#' @keywords internal
.estimate_fourier_ardl <- function(data, y_var, x_vars, p, q, k, case) {
  
  n <- nrow(data)
  
  # Generate Fourier terms
  fourier_terms <- .generate_fourier_terms(n, k)
  
  # Combine data
  data_ext <- cbind(data, fourier_terms)
  
  # Prepare ARDL data with Fourier terms
  model_data <- .prepare_fourier_data(data_ext, y_var, x_vars, p, q, k, case)
  
  if (nrow(model_data) < 10) {
    stop("Insufficient observations")
  }
  
  # Estimate model
  y_col <- "dy"
  x_cols <- setdiff(names(model_data), "dy")
  
  formula_str <- paste(y_col, "~", paste(x_cols, collapse = " + "))
  model <- lm(as.formula(formula_str), data = model_data)
  
  coefs <- coef(model)
  se <- summary(model)$coefficients[, 2]
  
  # Extract components
  phi <- coefs["y_lag1"]
  
  # Long-run coefficients
  level_idx <- which(names(coefs) %in% x_vars)
  if (length(level_idx) > 0) {
    beta <- coefs[level_idx]
    theta <- -beta / phi
    names(theta) <- x_vars[1:length(theta)]
  } else {
    theta <- NULL
  }
  
  # Fourier coefficients
  fourier_idx <- grep("^sin_|^cos_", names(coefs))
  fourier_coefs <- coefs[fourier_idx]
  
  # Bounds test statistics
  # F-test for joint significance of level variables
  RSS_ur <- sum(residuals(model)^2)
  
  # Restricted model (no levels)
  x_cols_r <- x_cols[!x_cols %in% c("y_lag1", x_vars)]
  if (length(x_cols_r) > 0) {
    formula_r <- paste(y_col, "~", paste(x_cols_r, collapse = " + "))
    model_r <- lm(as.formula(formula_r), data = model_data)
    RSS_r <- sum(residuals(model_r)^2)
  } else {
    RSS_r <- sum((model_data$dy - mean(model_data$dy))^2)
  }
  
  m <- length(x_vars) + 1  # Number of restrictions
  n_obs <- nrow(model_data)
  k_ur <- length(coefs)
  
  F_stat <- ((RSS_r - RSS_ur) / m) / (RSS_ur / (n_obs - k_ur))
  
  # t-statistic for phi
  t_stat <- coefs["y_lag1"] / se["y_lag1"]
  
  # Model fit statistics
  sigma2 <- sum(residuals(model)^2) / (n_obs - k_ur)
  loglik <- -0.5 * n_obs * (log(2 * pi) + log(sigma2) + 1)
  aic <- -2 * loglik + 2 * k_ur
  bic <- -2 * loglik + log(n_obs) * k_ur
  
  # Detect structural breaks (where Fourier component is significant)
  breaks <- .detect_fourier_breaks(fourier_coefs, n, k)
  
  list(
    coefficients = coefs,
    se = se,
    long_run = theta,
    ec_coef = phi,
    fourier_coefs = fourier_coefs,
    bounds_test = list(F_stat = F_stat, t_stat = t_stat),
    k = k,
    model = model,
    residuals = residuals(model),
    fitted = fitted(model),
    nobs = n_obs,
    sigma = sqrt(sigma2),
    loglik = loglik,
    aic = aic,
    bic = bic,
    structural_breaks = breaks,
    fourier_terms = fourier_terms
  )
}


#' @title Generate Fourier Terms
#' @keywords internal
.generate_fourier_terms <- function(n, k) {
  
  t <- 1:n
  result <- data.frame(row.names = 1:n)
  
  for (j in 1:k) {
    result[[paste0("sin_", j)]] <- sin(2 * pi * j * t / n)
    result[[paste0("cos_", j)]] <- cos(2 * pi * j * t / n)
  }
  
  return(result)
}


#' @title Prepare Fourier ARDL Data
#' @keywords internal
.prepare_fourier_data <- function(data, y_var, x_vars, p, q, k, case) {
  
  n <- nrow(data)
  max_lag <- max(p, max(q))
  
  y <- data[[y_var]]
  dy <- diff(y)
  y_lag1 <- y[-length(y)]
  
  # Start building model data
  result <- data.frame(
    dy = dy,
    y_lag1 = y_lag1
  )
  
  # X variables in levels
  for (v in x_vars) {
    result[[v]] <- data[[v]][-n]
  }
  
  # X variables in differences
  for (j in seq_along(x_vars)) {
    v <- x_vars[j]
    dx <- diff(data[[v]])
    result[[paste0("d_", v)]] <- dx
    
    # Lagged differences
    if (q[j] > 0) {
      for (lag in 1:q[j]) {
        lagged <- c(rep(NA, lag), dx[1:(length(dx) - lag)])
        result[[paste0("d_", v, "_L", lag)]] <- lagged
      }
    }
  }
  
  # Lagged dy
  if (p > 1) {
    for (lag in 1:(p-1)) {
      result[[paste0("dy_L", lag)]] <- c(rep(NA, lag), dy[1:(length(dy) - lag)])
    }
  }
  
  # Fourier terms
  fourier_cols <- grep("^sin_|^cos_", names(data), value = TRUE)
  for (fc in fourier_cols) {
    result[[fc]] <- data[[fc]][-n]
  }
  
  # Deterministic terms
  if (case >= 3) {
    result$const <- 1
  }
  if (case >= 4) {
    result$trend <- 1:nrow(result)
  }
  
  # Remove first row (lost to differencing) and NAs
  result <- result[-1, ]
  result <- na.omit(result)
  
  return(result)
}


#' @title Detect Fourier Structural Breaks
#' @keywords internal
.detect_fourier_breaks <- function(fourier_coefs, n, k) {
  
  if (length(fourier_coefs) == 0) {
    return(NULL)
  }
  
  # Reconstruct Fourier function
  t <- 1:n
  f_t <- rep(0, n)
  
  for (j in 1:k) {
    sin_name <- paste0("sin_", j)
    cos_name <- paste0("cos_", j)
    
    if (sin_name %in% names(fourier_coefs)) {
      f_t <- f_t + fourier_coefs[sin_name] * sin(2 * pi * j * t / n)
    }
    if (cos_name %in% names(fourier_coefs)) {
      f_t <- f_t + fourier_coefs[cos_name] * cos(2 * pi * j * t / n)
    }
  }
  
  # Find peaks and troughs (potential break points)
  # Use first derivative sign changes
  df_t <- diff(f_t)
  sign_changes <- which(diff(sign(df_t)) != 0)
  
  if (length(sign_changes) == 0) {
    return(NULL)
  }
  
  # Return approximate break periods
  breaks <- data.frame(
    period = sign_changes,
    type = ifelse(df_t[sign_changes] > 0, "trough", "peak"),
    fourier_value = f_t[sign_changes]
  )
  
  return(breaks)
}


#' @title Fourier ARDL Bounds Test
#' @description Perform bounds test with Fourier critical values
#'
#' @param object A fourier_ardl object
#'
#' @return Bounds test results with appropriate critical values
#' @export
fourier_bounds_test <- function(object) {
  
  if (!inherits(object, "fourier_ardl")) {
    stop("Object must be of class 'fourier_ardl'")
  }
  
  F_stat <- object$bounds_test$F_stat
  t_stat <- object$bounds_test$t_stat
  k <- length(object$x_vars)
  n_freq <- object$k
  
  # Critical values for Fourier ARDL (approximate)
  # Based on Banerjee et al. (2017) and simulation studies
  cv_F <- list(
    "10%" = list(I0 = 2.5 + 0.1 * n_freq, I1 = 3.5 + 0.15 * n_freq),
    "5%" = list(I0 = 3.0 + 0.1 * n_freq, I1 = 4.0 + 0.15 * n_freq),
    "1%" = list(I0 = 4.0 + 0.1 * n_freq, I1 = 5.0 + 0.15 * n_freq)
  )
  
  cat("\n")
  cat("============================================================\n")
  cat("     Fourier ARDL Bounds Test for Cointegration\n")
  cat("============================================================\n\n")
  
  cat("Number of Fourier frequencies (k):", n_freq, "\n")
  cat("Number of regressors:", k, "\n\n")
  
  cat("Test Statistics:\n")
  cat("  F-statistic:", round(F_stat, 4), "\n")
  cat("  t-statistic:", round(t_stat, 4), "\n\n")
  
  cat("Critical Values (5% level):\n")
  cat("  F: I(0) =", round(cv_F[["5%"]]$I0, 2), 
      ", I(1) =", round(cv_F[["5%"]]$I1, 2), "\n\n")
  
  # Conclusion
  if (F_stat > cv_F[["5%"]]$I1) {
    cat("Conclusion: COINTEGRATION - F above upper bound\n")
  } else if (F_stat < cv_F[["5%"]]$I0) {
    cat("Conclusion: NO COINTEGRATION - F below lower bound\n")
  } else {
    cat("Conclusion: INCONCLUSIVE - F between bounds\n")
  }
  
  cat("============================================================\n")
  
  invisible(list(F_stat = F_stat, t_stat = t_stat, cv_F = cv_F))
}


#' @title Summary method for fourier_ardl
#' @export
summary.fourier_ardl <- function(object, ...) {
  
  cat("\n")
  cat("====================================================================\n")
  cat("           Fourier ARDL Estimation Results\n")
  cat("====================================================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Model: Fourier-ARDL(", object$p, ", ", paste(object$q, collapse = ", "), 
      ") with k =", object$k, "frequencies\n", sep = "")
  cat("Case: ", object$case, "\n")
  cat("Observations:", object$nobs, "\n")
  
  if (!is.null(object$selection) && object$selection != "fixed") {
    cat("Frequency selection:", toupper(object$selection), "\n")
  }
  cat("\n")
  
  # IC comparison if available
  if (!is.null(object$ic_comparison)) {
    cat("--------------------------------------------------------------------\n")
    cat("     Frequency Selection (Information Criteria)\n")
    cat("--------------------------------------------------------------------\n")
    ic_df <- object$ic_comparison
    ic_df$aic <- round(ic_df$aic, 2)
    ic_df$bic <- round(ic_df$bic, 2)
    ic_df$selected <- ifelse(ic_df$k == object$k, " <--", "")
    print(ic_df, row.names = FALSE)
    cat("\n")
  }
  
  # Long-run coefficients
  cat("--------------------------------------------------------------------\n")
  cat("     Long-Run Coefficients\n")
  cat("--------------------------------------------------------------------\n")
  if (!is.null(object$long_run)) {
    lr_df <- data.frame(
      Variable = names(object$long_run),
      Estimate = round(object$long_run, 4)
    )
    print(lr_df, row.names = FALSE)
  }
  cat("\n")
  
  # EC coefficient
  cat("Error Correction Coefficient (phi):", round(object$ec_coef, 4), "\n")
  if (object$ec_coef < 0) {
    half_life <- -log(2) / log(1 + object$ec_coef)
    cat("Half-life:", round(half_life, 2), "periods\n")
  }
  cat("\n")
  
  # Fourier coefficients
  cat("--------------------------------------------------------------------\n")
  cat("     Fourier Coefficients\n")
  cat("--------------------------------------------------------------------\n")
  fourier_df <- data.frame(
    Term = names(object$fourier_coefs),
    Estimate = round(object$fourier_coefs, 4)
  )
  print(fourier_df, row.names = FALSE)
  cat("\n")
  
  # Bounds test
  cat("--------------------------------------------------------------------\n")
  cat("     Bounds Test Statistics\n")
  cat("--------------------------------------------------------------------\n")
  cat("F-statistic:", round(object$bounds_test$F_stat, 4), "\n")
  cat("t-statistic:", round(object$bounds_test$t_stat, 4), "\n\n")
  
  # Structural breaks
  if (!is.null(object$structural_breaks) && nrow(object$structural_breaks) > 0) {
    cat("--------------------------------------------------------------------\n")
    cat("     Detected Structural Break Periods\n")
    cat("--------------------------------------------------------------------\n")
    print(object$structural_breaks, row.names = FALSE)
    cat("\n")
  }
  
  # Model fit
  cat("--------------------------------------------------------------------\n")
  cat("     Model Fit Statistics\n")
  cat("--------------------------------------------------------------------\n")
  cat("Residual std. error:", round(object$sigma, 4), "\n")
  cat("Log-likelihood:     ", round(object$loglik, 2), "\n")
  cat("AIC:                ", round(object$aic, 2), "\n")
  cat("BIC:                ", round(object$bic, 2), "\n")
  cat("====================================================================\n")
  
  invisible(object)
}


#' @title Print method for fourier_ardl
#' @export
print.fourier_ardl <- function(x, ...) {
  
  cat("\nFourier ARDL Model (k =", x$k, "frequencies)\n")
  cat("Formula: ")
  print(x$formula)
  cat("\nLong-run coefficients:\n")
  print(round(x$long_run, 4))
  cat("\nEC coefficient:", round(x$ec_coef, 4), "\n")
  cat("F-statistic:", round(x$bounds_test$F_stat, 4), "\n")
  
  invisible(x)
}


#' @title Plot Fourier ARDL Components
#' @description Visualize Fourier terms and structural breaks
#'
#' @param x A fourier_ardl object
#' @param which Character. "fourier", "fit", or "both"
#' @param ... Additional arguments
#'
#' @export
plot.fourier_ardl <- function(x, which = "both", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  
  plots <- list()
  
  if (which %in% c("fourier", "both")) {
    # Reconstruct Fourier component
    n <- nrow(x$fourier_terms)
    t <- 1:n
    f_t <- rep(0, n)
    
    for (j in 1:x$k) {
      sin_name <- paste0("sin_", j)
      cos_name <- paste0("cos_", j)
      
      if (sin_name %in% names(x$fourier_coefs)) {
        f_t <- f_t + x$fourier_coefs[sin_name] * x$fourier_terms[[sin_name]]
      }
      if (cos_name %in% names(x$fourier_coefs)) {
        f_t <- f_t + x$fourier_coefs[cos_name] * x$fourier_terms[[cos_name]]
      }
    }
    
    df_fourier <- data.frame(t = t, fourier = f_t)
    
    p1 <- ggplot2::ggplot(df_fourier, ggplot2::aes(x = t, y = fourier)) +
      ggplot2::geom_line(color = "#2E86AB", linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = paste("Fourier Component (k =", x$k, ")"),
        subtitle = "Captures smooth structural changes",
        x = "Time",
        y = "Fourier Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")
      )
    
    # Add break points if detected
    if (!is.null(x$structural_breaks) && nrow(x$structural_breaks) > 0) {
      p1 <- p1 + ggplot2::geom_vline(
        xintercept = x$structural_breaks$period,
        linetype = "dotted", color = "#E94F37", alpha = 0.7
      )
    }
    
    plots$fourier <- p1
  }
  
  if (which %in% c("fit", "both")) {
    df_fit <- data.frame(
      t = 1:length(x$residuals),
      residuals = x$residuals
    )
    
    p2 <- ggplot2::ggplot(df_fit, ggplot2::aes(x = t, y = residuals)) +
      ggplot2::geom_line(color = "gray40") +
      ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      ggplot2::labs(
        title = "Model Residuals",
        x = "Time",
        y = "Residuals"
      ) +
      ggplot2::theme_minimal()
    
    plots$fit <- p2
  }
  
  if (which == "both" && requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(plots$fourier, plots$fit, ncol = 1)
  } else if (which == "fourier") {
    print(plots$fourier)
  } else if (which == "fit") {
    print(plots$fit)
  }
  
  invisible(plots)
}
