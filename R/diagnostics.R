#' @title ARDL Model Diagnostics
#' @description Comprehensive diagnostic tests for ARDL models
#'
#' @details
#' This function performs a battery of diagnostic tests commonly used to
#' assess the validity of ARDL models:
#'
#' \itemize{
#'   \item \strong{Serial Correlation}: Breusch-Godfrey LM test
#'   \item \strong{Heteroskedasticity}: Breusch-Pagan and ARCH tests
#'   \item \strong{Normality}: Jarque-Bera test on residuals
#'   \item \strong{Functional Form}: RESET test
#'   \item \strong{Stability}: CUSUM and CUSUMSQ tests
#' }
#'
#' @param model An estimated model object (panel_ardl, boot_ardl, qnardl, or fourier_ardl)
#' @param lags Integer. Number of lags for serial correlation tests (default: 4)
#' @param arch_lags Integer. Number of lags for ARCH test (default: 4)
#'
#' @return An object of class "ardl_diagnostics" containing all test results
#'
#' @examples
#' \dontrun{
#' # Estimate model
#' model <- fourier_ardl(gdp ~ investment, data = macro_data, p = 2, q = 2)
#'
#' # Run diagnostics
#' diag <- ardl_diagnostics(model)
#' summary(diag)
#' plot(diag)
#' }
#'
#' @export
#' @importFrom stats residuals fitted lm pchisq pf shapiro.test Box.test
#' @importFrom lmtest bgtest bptest resettest
ardl_diagnostics <- function(model, lags = 4, arch_lags = 4) {
  
  # Extract residuals and model
  if (inherits(model, "panel_ardl")) {
    resid <- model$residuals
    fitted_vals <- model$fitted
    lm_model <- NULL
  } else if (inherits(model, "boot_ardl")) {
    resid <- residuals(model$model)
    fitted_vals <- fitted(model$model)
    lm_model <- model$model
  } else if (inherits(model, "fourier_ardl")) {
    resid <- model$residuals
    fitted_vals <- model$fitted
    lm_model <- model$model
  } else if (inherits(model, "qnardl")) {
    # Use median quantile results
    tau_50 <- which.min(abs(model$tau - 0.5))
    resid <- model$results_by_tau[[tau_50]]$residuals
    fitted_vals <- model$results_by_tau[[tau_50]]$fitted
    lm_model <- model$results_by_tau[[tau_50]]$model
  } else {
    stop("Model must be of class panel_ardl, boot_ardl, qnardl, or fourier_ardl")
  }
  
  n <- length(resid)
  results <- list()
  
  # 1. Serial Correlation Test (Breusch-Godfrey)
  if (!is.null(lm_model)) {
    results$serial_corr <- tryCatch({
      lmtest::bgtest(lm_model, order = lags)
    }, error = function(e) {
      .manual_bg_test(resid, lags)
    })
  } else {
    results$serial_corr <- .manual_bg_test(resid, lags)
  }
  
  # 2. Ljung-Box Test
  results$ljung_box <- Box.test(resid, lag = lags, type = "Ljung-Box")
  
  # 3. Heteroskedasticity Test (Breusch-Pagan)
  if (!is.null(lm_model)) {
    results$hetero_bp <- tryCatch({
      lmtest::bptest(lm_model)
    }, error = function(e) {
      .manual_bp_test(resid, fitted_vals)
    })
  } else {
    results$hetero_bp <- .manual_bp_test(resid, fitted_vals)
  }
  
  # 4. ARCH Test
  results$arch <- .arch_test(resid, arch_lags)
  
  # 5. Normality Test (Jarque-Bera)
  results$normality <- .jarque_bera_test(resid)
  
  # 6. Shapiro-Wilk Test (for smaller samples)
  if (n <= 5000) {
    results$shapiro <- shapiro.test(resid)
  }
  
  # 7. RESET Test (Functional Form)
  if (!is.null(lm_model)) {
    results$reset <- tryCatch({
      lmtest::resettest(lm_model, power = 2:3)
    }, error = function(e) NULL)
  }
  
  # 8. CUSUM Statistics
  results$cusum <- .cusum_test(resid)
  
  # 9. CUSUM of Squares
  results$cusumsq <- .cusumsq_test(resid)
  
  # Store residuals for plotting
  results$residuals <- resid
  results$fitted <- fitted_vals
  results$nobs <- n
  
  class(results) <- c("ardl_diagnostics", "list")
  return(results)
}


#' @title Manual Breusch-Godfrey Test
#' @keywords internal
.manual_bg_test <- function(resid, lags) {
  
  n <- length(resid)
  
  # Create lagged residuals
  resid_lags <- sapply(1:lags, function(l) {
    c(rep(NA, l), resid[1:(n-l)])
  })
  
  # Auxiliary regression
  df <- data.frame(resid = resid, resid_lags)
  df <- na.omit(df)
  
  aux_model <- lm(resid ~ ., data = df)
  r_sq <- summary(aux_model)$r.squared
  
  # LM statistic
  LM <- (n - lags) * r_sq
  p_value <- 1 - pchisq(LM, df = lags)
  
  list(
    statistic = LM,
    parameter = lags,
    p.value = p_value,
    method = "Breusch-Godfrey LM Test"
  )
}


#' @title Manual Breusch-Pagan Test
#' @keywords internal
.manual_bp_test <- function(resid, fitted) {
  
  n <- length(resid)
  resid_sq <- resid^2
  
  # Regression of squared residuals on fitted values
  aux_model <- lm(resid_sq ~ fitted)
  
  # Test statistic
  RSS <- sum((resid_sq - fitted(aux_model))^2)
  TSS <- sum((resid_sq - mean(resid_sq))^2)
  r_sq <- 1 - RSS/TSS
  
  LM <- n * r_sq
  p_value <- 1 - pchisq(LM, df = 1)
  
  list(
    statistic = LM,
    parameter = 1,
    p.value = p_value,
    method = "Breusch-Pagan Test"
  )
}


#' @title ARCH Test
#' @keywords internal
.arch_test <- function(resid, lags) {
  
  n <- length(resid)
  resid_sq <- resid^2
  
  # Create lagged squared residuals
  resid_sq_lags <- sapply(1:lags, function(l) {
    c(rep(NA, l), resid_sq[1:(n-l)])
  })
  
  df <- data.frame(resid_sq = resid_sq, resid_sq_lags)
  df <- na.omit(df)
  
  aux_model <- lm(resid_sq ~ ., data = df)
  r_sq <- summary(aux_model)$r.squared
  
  # LM statistic
  LM <- (n - lags) * r_sq
  p_value <- 1 - pchisq(LM, df = lags)
  
  list(
    statistic = LM,
    parameter = lags,
    p.value = p_value,
    method = "ARCH-LM Test"
  )
}


#' @title Jarque-Bera Test
#' @keywords internal
.jarque_bera_test <- function(resid) {
  
  n <- length(resid)
  m <- mean(resid)
  s <- sd(resid)
  
  # Standardized residuals
  z <- (resid - m) / s
  
  # Skewness and Kurtosis
  skewness <- mean(z^3)
  kurtosis <- mean(z^4)
  
  # Jarque-Bera statistic
  JB <- n * (skewness^2 / 6 + (kurtosis - 3)^2 / 24)
  p_value <- 1 - pchisq(JB, df = 2)
  
  list(
    statistic = JB,
    parameter = 2,
    p.value = p_value,
    skewness = skewness,
    kurtosis = kurtosis,
    method = "Jarque-Bera Normality Test"
  )
}


#' @title CUSUM Test
#' @keywords internal
.cusum_test <- function(resid) {
  
  n <- length(resid)
  sigma <- sd(resid)
  
  # Cumulative sum of standardized residuals
  cusum <- cumsum(resid) / sigma
  
  # Critical bounds (5% level)
  # Approximate: +/- 0.948 * sqrt(n)
  bound <- 0.948 * sqrt(n) * (1:n) / n
  
  # Check if CUSUM crosses bounds
  crosses <- any(abs(cusum) > bound)
  
  list(
    cusum = cusum,
    upper_bound = bound,
    lower_bound = -bound,
    crosses_bounds = crosses,
    method = "CUSUM Test for Parameter Stability"
  )
}


#' @title CUSUM of Squares Test
#' @keywords internal
.cusumsq_test <- function(resid) {
  
  n <- length(resid)
  resid_sq <- resid^2
  
  # Cumulative sum of squares
  cusumsq <- cumsum(resid_sq) / sum(resid_sq)
  
  # Expected value under null: t/n
  expected <- (1:n) / n
  
  # Critical bounds (5% level, approximate)
  # Based on Brownian bridge
  bound <- 1.358 / sqrt(n)
  
  upper <- expected + bound
  lower <- expected - bound
  
  # Check if CUSUMSQ crosses bounds
  crosses <- any(cusumsq > upper | cusumsq < lower)
  
  list(
    cusumsq = cusumsq,
    expected = expected,
    upper_bound = upper,
    lower_bound = lower,
    crosses_bounds = crosses,
    method = "CUSUM of Squares Test"
  )
}


#' @title Summary method for ardl_diagnostics
#' @export
summary.ardl_diagnostics <- function(object, ...) {
  
  cat("\n")
  cat("====================================================================\n")
  cat("              ARDL Model Diagnostic Tests\n")
  cat("====================================================================\n\n")
  
  cat("Number of observations:", object$nobs, "\n\n")
  
  # Serial Correlation
  cat("--------------------------------------------------------------------\n")
  cat("  SERIAL CORRELATION\n")
  cat("--------------------------------------------------------------------\n")
  
  if (!is.null(object$serial_corr)) {
    cat("Breusch-Godfrey LM Test:\n")
    cat("  Statistic:", round(object$serial_corr$statistic, 4), "\n")
    cat("  Lags:     ", object$serial_corr$parameter, "\n")
    cat("  p-value:  ", format.pval(object$serial_corr$p.value, digits = 4), 
        .signif_stars(object$serial_corr$p.value), "\n")
  }
  
  cat("\nLjung-Box Test:\n")
  cat("  Statistic:", round(object$ljung_box$statistic, 4), "\n")
  cat("  p-value:  ", format.pval(object$ljung_box$p.value, digits = 4),
      .signif_stars(object$ljung_box$p.value), "\n")
  cat("\n")
  
  # Heteroskedasticity
  cat("--------------------------------------------------------------------\n")
  cat("  HETEROSKEDASTICITY\n")
  cat("--------------------------------------------------------------------\n")
  
  if (!is.null(object$hetero_bp)) {
    cat("Breusch-Pagan Test:\n")
    cat("  Statistic:", round(object$hetero_bp$statistic, 4), "\n")
    cat("  p-value:  ", format.pval(object$hetero_bp$p.value, digits = 4),
        .signif_stars(object$hetero_bp$p.value), "\n")
  }
  
  cat("\nARCH-LM Test:\n")
  cat("  Statistic:", round(object$arch$statistic, 4), "\n")
  cat("  Lags:     ", object$arch$parameter, "\n")
  cat("  p-value:  ", format.pval(object$arch$p.value, digits = 4),
      .signif_stars(object$arch$p.value), "\n")
  cat("\n")
  
  # Normality
  cat("--------------------------------------------------------------------\n")
  cat("  NORMALITY\n")
  cat("--------------------------------------------------------------------\n")
  
  cat("Jarque-Bera Test:\n")
  cat("  Statistic:", round(object$normality$statistic, 4), "\n")
  cat("  Skewness: ", round(object$normality$skewness, 4), "\n")
  cat("  Kurtosis: ", round(object$normality$kurtosis, 4), "\n")
  cat("  p-value:  ", format.pval(object$normality$p.value, digits = 4),
      .signif_stars(object$normality$p.value), "\n")
  
  if (!is.null(object$shapiro)) {
    cat("\nShapiro-Wilk Test:\n")
    cat("  Statistic:", round(object$shapiro$statistic, 4), "\n")
    cat("  p-value:  ", format.pval(object$shapiro$p.value, digits = 4),
        .signif_stars(object$shapiro$p.value), "\n")
  }
  cat("\n")
  
  # Functional Form
  if (!is.null(object$reset)) {
    cat("--------------------------------------------------------------------\n")
    cat("  FUNCTIONAL FORM\n")
    cat("--------------------------------------------------------------------\n")
    cat("RESET Test:\n")
    cat("  Statistic:", round(object$reset$statistic, 4), "\n")
    cat("  p-value:  ", format.pval(object$reset$p.value, digits = 4),
        .signif_stars(object$reset$p.value), "\n")
    cat("\n")
  }
  
  # Stability
  cat("--------------------------------------------------------------------\n")
  cat("  PARAMETER STABILITY\n")
  cat("--------------------------------------------------------------------\n")
  
  cat("CUSUM Test:\n")
  cat("  Crosses bounds:", ifelse(object$cusum$crosses_bounds, "YES ***", "NO"), "\n")
  
  cat("\nCUSUM of Squares Test:\n")
  cat("  Crosses bounds:", ifelse(object$cusumsq$crosses_bounds, "YES ***", "NO"), "\n")
  
  cat("\n")
  cat("--------------------------------------------------------------------\n")
  cat("  SUMMARY\n")
  cat("--------------------------------------------------------------------\n")
  
  # Summary assessment
  issues <- c()
  
  if (!is.null(object$serial_corr) && object$serial_corr$p.value < 0.05) {
    issues <- c(issues, "Serial correlation detected")
  }
  if (!is.null(object$hetero_bp) && object$hetero_bp$p.value < 0.05) {
    issues <- c(issues, "Heteroskedasticity detected")
  }
  if (object$arch$p.value < 0.05) {
    issues <- c(issues, "ARCH effects detected")
  }
  if (object$normality$p.value < 0.05) {
    issues <- c(issues, "Non-normal residuals")
  }
  if (object$cusum$crosses_bounds || object$cusumsq$crosses_bounds) {
    issues <- c(issues, "Parameter instability detected")
  }
  
  if (length(issues) == 0) {
    cat("All diagnostic tests passed at 5% significance level.\n")
    cat("Model appears to be well-specified.\n")
  } else {
    cat("Potential issues detected:\n")
    for (issue in issues) {
      cat("  -", issue, "\n")
    }
  }
  
  cat("====================================================================\n")
  cat("Signif. codes: 0 '***' 0.01 '**' 0.05 '*' 0.1 '.' (reject null)\n")
  
  invisible(object)
}


#' @title Significance Stars Helper
#' @keywords internal
.signif_stars <- function(p) {
  if (p < 0.01) return(" ***")
  if (p < 0.05) return(" **")
  if (p < 0.1) return(" *")
  return("")
}


#' @title Plot Diagnostics
#' @export
plot.ardl_diagnostics <- function(x, which = 1:4, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  
  plots <- list()
  n <- length(x$residuals)
  t <- 1:n
  
  # 1. Residuals vs Time
  if (1 %in% which) {
    df1 <- data.frame(t = t, residuals = x$residuals)
    
    plots[[1]] <- ggplot2::ggplot(df1, ggplot2::aes(x = t, y = residuals)) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Residuals over Time", x = "Time", y = "Residuals") +
      ggplot2::theme_minimal()
  }
  
  # 2. Residuals vs Fitted
  if (2 %in% which) {
    df2 <- data.frame(fitted = x$fitted, residuals = x$residuals)
    
    plots[[2]] <- ggplot2::ggplot(df2, ggplot2::aes(x = fitted, y = residuals)) +
      ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "orange") +
      ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
      ggplot2::theme_minimal()
  }
  
  # 3. Q-Q Plot
  if (3 %in% which) {
    df3 <- data.frame(residuals = x$residuals)
    
    plots[[3]] <- ggplot2::ggplot(df3, ggplot2::aes(sample = residuals)) +
      ggplot2::stat_qq(color = "steelblue") +
      ggplot2::stat_qq_line(color = "red") +
      ggplot2::labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggplot2::theme_minimal()
  }
  
  # 4. CUSUM Plot
  if (4 %in% which) {
    df4 <- data.frame(
      t = t,
      cusum = x$cusum$cusum,
      upper = x$cusum$upper_bound,
      lower = x$cusum$lower_bound
    )
    
    plots[[4]] <- ggplot2::ggplot(df4, ggplot2::aes(x = t)) +
      ggplot2::geom_line(ggplot2::aes(y = cusum), color = "steelblue", linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(y = upper), color = "red", linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(y = lower), color = "red", linetype = "dashed") +
      ggplot2::labs(title = "CUSUM Test", x = "Time", y = "CUSUM") +
      ggplot2::theme_minimal()
  }
  
  # Arrange plots
  if (requireNamespace("gridExtra", quietly = TRUE) && length(plots) > 1) {
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    for (p in plots) {
      print(p)
    }
  }
  
  invisible(plots)
}
