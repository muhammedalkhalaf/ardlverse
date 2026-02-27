#' @title Panel Nonlinear ARDL (Panel NARDL)
#' @description Panel data estimation with asymmetric/nonlinear ARDL specification.
#'
#' @details
#' The Panel NARDL model extends the Shin et al. (2014) NARDL framework to
#' panel data settings. It allows for asymmetric short-run and long-run
#' effects while controlling for cross-sectional heterogeneity.
#'
#' Three estimators are available:
#' \itemize{
#'   \item \strong{PMG}: Pooled Mean Group - constrains long-run coefficients to be
#'     homogeneous while allowing short-run heterogeneity
#'   \item \strong{MG}: Mean Group - estimates separate models for each unit
#'     and averages coefficients
#'   \item \strong{DFE}: Dynamic Fixed Effects - assumes all coefficients
#'     are homogeneous except fixed effects
#' }
#'
#' @param formula A formula specifying the model: y ~ x1 + x2 + ... | group
#' @param data A data frame in long panel format
#' @param id Character. Name of the group/panel identifier variable
#' @param time Character. Name of the time variable
#' @param p Integer. Number of lags for dependent variable (default: 1)
#' @param q Integer or vector. Number of lags for independent variables (default: 1)
#' @param estimator Character. One of "pmg", "mg", or "dfe" (default: "pmg")
#' @param threshold Numeric. Threshold for asymmetric decomposition (default: 0)
#' @param effect Character. Type of effects: "individual", "time", "twoways" (default: "individual")
#' @param bootstrap Logical. Use bootstrap for inference (default: FALSE)
#' @param nboot Number of bootstrap replications (default: 500)
#'
#' @return An object of class "pnardl" containing:
#' \itemize{
#'   \item \code{coefficients}: Estimated coefficients
#'   \item \code{long_run_pos}: Long-run coefficients for positive changes
#'   \item \code{long_run_neg}: Long-run coefficients for negative changes
#'   \item \code{short_run}: Short-run coefficients
#'   \item \code{ec_coef}: Error correction coefficient
#'   \item \code{asymmetry_test}: Wald test for long-run asymmetry
#'   \item \code{unit_results}: Individual unit estimation results (for MG)
#'   \item \code{hausman}: Hausman test comparing PMG vs MG
#' }
#'
#' @references
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#'
#' Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation
#' of dynamic heterogeneous panels. Journal of the American Statistical Association.
#'
#' @examples
#' \dontrun{
#' # Generate panel data
#' data <- generate_panel_data(N = 20, T = 50)
#'
#' # Panel NARDL with PMG estimator
#' result <- pnardl(
#'   y ~ x1 + x2,
#'   data = data,
#'   id = "id",
#'   time = "time",
#'   estimator = "pmg"
#' )
#' summary(result)
#'
#' # Test for asymmetry
#' asymmetry_test(result)
#'
#' # Dynamic multipliers
#' plot(dynamic_multipliers(result))
#' }
#'
#' @export
pnardl <- function(formula, data, id, time, p = 1, q = 1,
                   estimator = c("pmg", "mg", "dfe"),
                   threshold = 0,
                   effect = c("individual", "time", "twoways"),
                   bootstrap = FALSE, nboot = 500) {
  
  estimator <- match.arg(estimator)
  effect <- match.arg(effect)
  
  # Parse formula
  formula_vars <- all.vars(formula)
  y_var <- formula_vars[1]
  x_vars <- formula_vars[-1]
  k <- length(x_vars)
  
  # Handle q
  if (length(q) == 1) q <- rep(q, k)
  
  # Get panel structure
  data <- data[order(data[[id]], data[[time]]), ]
  groups <- unique(data[[id]])
  N <- length(groups)
  T_total <- nrow(data) / N
  
  # Decompose variables into positive and negative components
  data_decomp <- .panel_nardl_decompose(data, x_vars, id, time, threshold)
  
  # New variable names
  x_vars_pos <- paste0(x_vars, "_pos")
  x_vars_neg <- paste0(x_vars, "_neg")
  x_vars_all <- c(x_vars_pos, x_vars_neg)
  k_total <- length(x_vars_all)
  
  # Estimate based on chosen estimator
  if (estimator == "pmg") {
    result <- .pnardl_pmg(data_decomp, y_var, x_vars_all, x_vars, id, time, 
                          p, q, effect, k)
  } else if (estimator == "mg") {
    result <- .pnardl_mg(data_decomp, y_var, x_vars_all, x_vars, id, time,
                         p, q, effect, k, groups)
  } else {
    result <- .pnardl_dfe(data_decomp, y_var, x_vars_all, x_vars, id, time,
                          p, q, effect, k)
  }
  
  # Bootstrap standard errors if requested
  if (bootstrap) {
    boot_se <- .pnardl_bootstrap(data_decomp, y_var, x_vars_all, id, time,
                                  p, q, estimator, effect, nboot, groups)
    result$boot_se <- boot_se
    result$ci_lower <- result$coefficients - 1.96 * boot_se
    result$ci_upper <- result$coefficients + 1.96 * boot_se
  }
  
  # Asymmetry test
  result$asymmetry_test <- .test_panel_asymmetry(result, x_vars, k)
  
  # Hausman test (PMG vs MG)
  if (estimator == "pmg") {
    mg_result <- .pnardl_mg(data_decomp, y_var, x_vars_all, x_vars, id, time,
                            p, q, effect, k, groups)
    result$hausman <- .panel_hausman_test(result, mg_result)
  }
  
  # Add metadata
  result$call <- match.call()
  result$estimator <- estimator
  result$effect <- effect
  result$threshold <- threshold
  result$N <- N
  result$T <- T_total
  result$k <- k
  result$p <- p
  result$q <- q
  result$x_vars <- x_vars
  result$y_var <- y_var
  result$groups <- groups
  
  class(result) <- "pnardl"
  return(result)
}


#' @title Panel NARDL Decomposition
#' @keywords internal
.panel_nardl_decompose <- function(data, x_vars, id, time, threshold) {
  groups <- unique(data[[id]])
  
  result_list <- list()
  
  for (g in groups) {
    group_data <- data[data[[id]] == g, ]
    group_data <- group_data[order(group_data[[time]]), ]
    
    for (v in x_vars) {
      x <- group_data[[v]]
      dx <- c(0, diff(x))
      
      pos_changes <- ifelse(dx > threshold, dx - threshold, 0)
      neg_changes <- ifelse(dx < -threshold, dx + threshold, 0)
      
      group_data[[paste0(v, "_pos")]] <- cumsum(pos_changes)
      group_data[[paste0(v, "_neg")]] <- cumsum(neg_changes)
    }
    
    result_list[[as.character(g)]] <- group_data
  }
  
  do.call(rbind, result_list)
}


#' @title Panel NARDL PMG Estimator
#' @keywords internal
.pnardl_pmg <- function(data, y_var, x_vars, orig_vars, id, time, p, q, effect, k) {
  groups <- unique(data[[id]])
  N <- length(groups)
  k_total <- length(x_vars)
  
  # Expand q for pos/neg variables
  if (length(q) == k) {
    q <- rep(q, 2)
  }
  
  # Initialize long-run coefficients (to be estimated pooled)
  lr_init <- rep(0, k_total)
  
  # Iterative PMG estimation
  max_iter <- 100
  tol <- 1e-6
  
  lr_coefs <- lr_init
  
  for (iter in 1:max_iter) {
    lr_old <- lr_coefs
    
    # Step 1: Estimate short-run for each group, given long-run
    sr_results <- list()
    ec_coefs <- numeric(N)
    
    for (i in seq_along(groups)) {
      g <- groups[i]
      group_data <- data[data[[id]] == g, ]
      group_data <- group_data[order(group_data[[time]]), ]
      
      # Build error correction term
      y <- group_data[[y_var]]
      X <- as.matrix(group_data[, x_vars, drop = FALSE])
      n <- length(y)
      
      # EC term: y - sum(lr * x)
      ec_term <- y - X %*% lr_coefs
      
      # Build ARDL model with EC term
      max_lag <- max(p, max(q))
      valid_idx <- (max_lag + 1):n
      
      dy <- diff(y)[(max_lag):(n-1)]
      ec_lag <- ec_term[valid_idx - 1]
      
      # Lagged differences
      design_list <- list(ec_lag = ec_lag)
      
      # Lagged dy
      for (j in 1:p) {
        design_list[[paste0("dy_l", j)]] <- diff(y)[(max_lag - j):(n - 1 - j)]
      }
      
      # Contemporaneous and lagged dx
      for (j in 1:k_total) {
        dx_j <- diff(X[, j])
        for (l in 0:(q[min(j, length(q))] - 1)) {
          design_list[[paste0("dx", j, "_l", l)]] <- dx_j[(max_lag - l):(n - 1 - l)]
        }
      }
      
      design <- do.call(cbind, design_list)
      
      # Add intercept
      design <- cbind(design, intercept = 1)
      
      # Estimate
      model_g <- stats::lm(dy ~ design - 1)
      coefs_g <- stats::coef(model_g)
      
      ec_coefs[i] <- coefs_g[1]
      sr_results[[i]] <- list(coefs = coefs_g, model = model_g)
    }
    
    # Step 2: Update long-run coefficients (concentrated likelihood)
    # Use average EC coefficient
    phi_avg <- mean(ec_coefs)
    
    # Stack data for pooled long-run estimation
    lr_data <- .stack_lr_data(data, y_var, x_vars, id, time, phi_avg, sr_results, groups)
    
    if (nrow(lr_data$Y) > 0) {
      lr_model <- stats::lm(lr_data$Y ~ lr_data$X - 1)
      lr_coefs <- stats::coef(lr_model)
    }
    
    # Check convergence
    if (max(abs(lr_coefs - lr_old)) < tol) break
  }
  
  # Compute standard errors
  vcov_lr <- stats::vcov(lr_model)
  se_lr <- sqrt(diag(vcov_lr))
  
  # Average short-run coefficients
  sr_avg <- colMeans(do.call(rbind, lapply(sr_results, function(x) x$coefs)))
  
  # Separate long-run pos and neg
  lr_pos <- lr_coefs[1:k]
  lr_neg <- lr_coefs[(k+1):(2*k)]
  
  list(
    coefficients = c(ec = mean(ec_coefs), lr_coefs, sr_avg[-1]),
    long_run = lr_coefs,
    long_run_pos = lr_pos,
    long_run_neg = lr_neg,
    short_run = sr_avg[-1],
    ec_coef = mean(ec_coefs),
    ec_by_unit = ec_coefs,
    se_lr = se_lr,
    convergence = list(iterations = iter, converged = iter < max_iter),
    unit_results = sr_results
  )
}


#' @title Panel NARDL MG Estimator
#' @keywords internal
.pnardl_mg <- function(data, y_var, x_vars, orig_vars, id, time, p, q, effect, k, groups) {
  N <- length(groups)
  k_total <- length(x_vars)
  
  if (length(q) == k) q <- rep(q, 2)
  
  # Estimate for each group
  unit_results <- list()
  all_lr <- matrix(NA, N, k_total)
  all_ec <- numeric(N)
  
  for (i in seq_along(groups)) {
    g <- groups[i]
    group_data <- data[data[[id]] == g, ]
    group_data <- group_data[order(group_data[[time]]), ]
    
    tryCatch({
      result_g <- .estimate_unit_nardl(group_data, y_var, x_vars, p, q)
      unit_results[[i]] <- result_g
      all_lr[i, ] <- result_g$long_run
      all_ec[i] <- result_g$ec_coef
    }, error = function(e) {
      unit_results[[i]] <<- NULL
      all_lr[i, ] <<- NA
      all_ec[i] <<- NA
    })
  }
  
  # Average coefficients
  lr_avg <- colMeans(all_lr, na.rm = TRUE)
  lr_se <- apply(all_lr, 2, function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  
  ec_avg <- mean(all_ec, na.rm = TRUE)
  ec_se <- stats::sd(all_ec, na.rm = TRUE) / sqrt(sum(!is.na(all_ec)))
  
  lr_pos <- lr_avg[1:k]
  lr_neg <- lr_avg[(k+1):(2*k)]
  
  list(
    coefficients = c(ec = ec_avg, lr_avg),
    long_run = lr_avg,
    long_run_pos = lr_pos,
    long_run_neg = lr_neg,
    ec_coef = ec_avg,
    ec_se = ec_se,
    se_lr = lr_se,
    unit_results = unit_results,
    unit_lr = all_lr,
    unit_ec = all_ec
  )
}


#' @title Panel NARDL DFE Estimator
#' @keywords internal
.pnardl_dfe <- function(data, y_var, x_vars, orig_vars, id, time, p, q, effect, k) {
  groups <- unique(data[[id]])
  N <- length(groups)
  k_total <- length(x_vars)
  
  if (length(q) == k) q <- rep(q, 2)
  
  # Create panel with fixed effects demeaning
  data_demeaned <- .demean_panel(data, y_var, x_vars, id)
  
  # Stack all data
  all_data <- list()
  for (g in groups) {
    group_data <- data_demeaned[data_demeaned[[id]] == g, ]
    group_data <- group_data[order(group_data[[time]]), ]
    
    y <- group_data[[y_var]]
    X <- as.matrix(group_data[, x_vars, drop = FALSE])
    n <- length(y)
    
    max_lag <- max(p, max(q))
    valid_idx <- (max_lag + 1):n
    
    dy <- diff(y)[(max_lag):(n-1)]
    y_lag <- y[valid_idx - 1]
    x_lag <- X[valid_idx - 1, , drop = FALSE]
    
    # Lagged dy
    dy_lags <- matrix(NA, length(valid_idx), p)
    for (j in 1:p) {
      dy_lags[, j] <- diff(y)[(max_lag - j):(n - 1 - j)]
    }
    
    # Lagged dx
    dx_list <- list()
    for (j in 1:k_total) {
      dx_j <- diff(X[, j])
      dx_lags_j <- matrix(NA, length(valid_idx), q[min(j, length(q))])
      for (l in 0:(q[min(j, length(q))] - 1)) {
        dx_lags_j[, l + 1] <- dx_j[(max_lag - l):(n - 1 - l)]
      }
      dx_list[[j]] <- dx_lags_j
    }
    dx_all <- do.call(cbind, dx_list)
    
    all_data[[as.character(g)]] <- list(
      dy = dy,
      y_lag = y_lag,
      x_lag = x_lag,
      dy_lags = dy_lags,
      dx_all = dx_all
    )
  }
  
  # Stack
  dy_stack <- unlist(lapply(all_data, function(x) x$dy))
  design_stack <- do.call(rbind, lapply(all_data, function(x) {
    cbind(x$y_lag, x$x_lag, x$dy_lags, x$dx_all)
  }))
  
  # Estimate
  model <- stats::lm(dy_stack ~ design_stack - 1)
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  se <- sqrt(diag(vcov_mat))
  
  ec_coef <- coefs[1]
  lr_coefs <- if (abs(ec_coef) > 1e-10) {
    -coefs[2:(1 + k_total)] / ec_coef
  } else {
    rep(NA, k_total)
  }
  
  lr_pos <- lr_coefs[1:k]
  lr_neg <- lr_coefs[(k+1):(2*k)]
  
  list(
    coefficients = coefs,
    long_run = lr_coefs,
    long_run_pos = lr_pos,
    long_run_neg = lr_neg,
    ec_coef = ec_coef,
    se = se,
    model = model,
    vcov = vcov_mat
  )
}


#' @title Estimate Unit NARDL
#' @keywords internal
.estimate_unit_nardl <- function(data, y_var, x_vars, p, q) {
  y <- data[[y_var]]
  X <- as.matrix(data[, x_vars, drop = FALSE])
  n <- length(y)
  k <- ncol(X)
  
  max_lag <- max(p, max(q))
  valid_idx <- (max_lag + 1):n
  
  dy <- diff(y)[(max_lag):(n-1)]
  y_lag <- y[valid_idx - 1]
  x_lag <- X[valid_idx - 1, , drop = FALSE]
  
  # Build design
  design <- cbind(y_lag, x_lag)
  
  # Add lagged differences
  for (j in 1:p) {
    design <- cbind(design, diff(y)[(max_lag - j):(n - 1 - j)])
  }
  
  for (j in 1:k) {
    dx_j <- diff(X[, j])
    for (l in 0:(q[min(j, length(q))] - 1)) {
      design <- cbind(design, dx_j[(max_lag - l):(n - 1 - l)])
    }
  }
  
  design <- cbind(design, 1)  # intercept
  
  model <- stats::lm(dy ~ design - 1)
  coefs <- stats::coef(model)
  
  ec_coef <- coefs[1]
  lr_coefs <- if (abs(ec_coef) > 1e-10) {
    -coefs[2:(1 + k)] / ec_coef
  } else {
    rep(NA, k)
  }
  
  list(
    coefficients = coefs,
    long_run = lr_coefs,
    ec_coef = ec_coef,
    model = model
  )
}


#' @title Stack Long-Run Data for PMG
#' @keywords internal
.stack_lr_data <- function(data, y_var, x_vars, id, time, phi, sr_results, groups) {
  Y_list <- list()
  X_list <- list()
  
  for (i in seq_along(groups)) {
    if (is.null(sr_results[[i]])) next
    
    g <- groups[i]
    group_data <- data[data[[id]] == g, ]
    group_data <- group_data[order(group_data[[time]]), ]
    
    y <- group_data[[y_var]]
    X <- as.matrix(group_data[, x_vars, drop = FALSE])
    n <- length(y)
    
    # Transform for long-run estimation
    y_transformed <- y[2:n] - (1 + phi) * y[1:(n-1)]
    X_transformed <- X[2:n, , drop = FALSE] - (1 + phi) * X[1:(n-1), , drop = FALSE]
    
    Y_list[[i]] <- y_transformed
    X_list[[i]] <- X_transformed
  }
  
  list(
    Y = matrix(unlist(Y_list), ncol = 1),
    X = do.call(rbind, X_list)
  )
}


#' @title Demean Panel Data
#' @keywords internal
.demean_panel <- function(data, y_var, x_vars, id) {
  groups <- unique(data[[id]])
  
  for (g in groups) {
    idx <- data[[id]] == g
    
    data[idx, y_var] <- data[idx, y_var] - mean(data[idx, y_var])
    
    for (v in x_vars) {
      data[idx, v] <- data[idx, v] - mean(data[idx, v])
    }
  }
  
  data
}


#' @title Test Panel Asymmetry
#' @keywords internal
.test_panel_asymmetry <- function(result, x_vars, k) {
  tests <- list()
  
  lr_pos <- result$long_run_pos
  lr_neg <- result$long_run_neg
  
  for (i in 1:k) {
    diff_lr <- lr_pos[i] - lr_neg[i]
    
    # Approximate Wald test
    if (!is.null(result$se_lr) && length(result$se_lr) >= 2*k) {
      se_diff <- sqrt(result$se_lr[i]^2 + result$se_lr[k + i]^2)
      wald <- (diff_lr / se_diff)^2
      p_value <- 1 - stats::pchisq(wald, 1)
    } else {
      wald <- NA
      p_value <- NA
    }
    
    tests[[x_vars[i]]] <- list(
      pos_coef = lr_pos[i],
      neg_coef = lr_neg[i],
      difference = diff_lr,
      wald = wald,
      p_value = p_value,
      asymmetric = if (!is.na(p_value)) p_value < 0.05 else NA
    )
  }
  
  tests
}


#' @title Panel Hausman Test
#' @keywords internal
.panel_hausman_test <- function(pmg_result, mg_result) {
  # Compare PMG (efficient under H0) vs MG (consistent)
  lr_pmg <- pmg_result$long_run
  lr_mg <- mg_result$long_run
  
  diff_coef <- lr_pmg - lr_mg
  
  # Approximate variance of difference
  var_mg <- if (!is.null(mg_result$se_lr)) mg_result$se_lr^2 else rep(0.01, length(lr_mg))
  var_pmg <- if (!is.null(pmg_result$se_lr)) pmg_result$se_lr^2 else rep(0.01, length(lr_pmg))
  
  var_diff <- var_mg - var_pmg
  var_diff[var_diff <= 0] <- 0.001  # Avoid negative variance
  
  # Hausman statistic
  H <- sum(diff_coef^2 / var_diff)
  df <- length(diff_coef)
  p_value <- 1 - stats::pchisq(H, df)
  
  list(
    statistic = H,
    df = df,
    p_value = p_value,
    conclusion = if (p_value > 0.05) "PMG preferred (cannot reject H0)" else "MG preferred (reject H0)"
  )
}


#' @title Panel NARDL Bootstrap
#' @keywords internal
.pnardl_bootstrap <- function(data, y_var, x_vars, id, time, p, q,
                               estimator, effect, nboot, groups) {
  N <- length(groups)
  n_coefs <- 1 + length(x_vars)  # EC + long-run
  
  boot_coefs <- matrix(NA, nboot, n_coefs)
  
  for (b in 1:nboot) {
    # Resample groups with replacement
    boot_groups <- sample(groups, N, replace = TRUE)
    
    boot_data_list <- list()
    for (i in seq_along(boot_groups)) {
      g <- boot_groups[i]
      gdata <- data[data[[id]] == g, ]
      gdata[[id]] <- paste0("boot_", i)  # Rename to avoid duplicates
      boot_data_list[[i]] <- gdata
    }
    boot_data <- do.call(rbind, boot_data_list)
    
    tryCatch({
      if (estimator == "pmg") {
        boot_result <- .pnardl_pmg(boot_data, y_var, x_vars, 
                                    unique(substr(x_vars, 1, nchar(x_vars) - 4)),
                                    id, time, p, q, effect, length(x_vars)/2)
      } else if (estimator == "mg") {
        boot_result <- .pnardl_mg(boot_data, y_var, x_vars,
                                   unique(substr(x_vars, 1, nchar(x_vars) - 4)),
                                   id, time, p, q, effect, length(x_vars)/2,
                                   unique(boot_data[[id]]))
      } else {
        boot_result <- .pnardl_dfe(boot_data, y_var, x_vars,
                                    unique(substr(x_vars, 1, nchar(x_vars) - 4)),
                                    id, time, p, q, effect, length(x_vars)/2)
      }
      
      boot_coefs[b, ] <- c(boot_result$ec_coef, boot_result$long_run)
    }, error = function(e) NULL)
  }
  
  apply(boot_coefs, 2, function(x) stats::sd(x, na.rm = TRUE))
}


#' @export
print.pnardl <- function(x, ...) {
  cat("\n")
  cat("Panel Nonlinear ARDL (Panel NARDL)\n")
  cat(paste(rep("=", 45), collapse = ""), "\n\n")
  
  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Panels:", x$N, "| Time periods:", x$T, "\n")
  cat("Threshold:", x$threshold, "\n\n")
  
  cat("Error Correction Coefficient:", round(x$ec_coef, 4), "\n\n")
  
  cat("Long-Run Coefficients:\n")
  cat("  Positive changes:", paste(round(x$long_run_pos, 4), collapse = ", "), "\n")
  cat("  Negative changes:", paste(round(x$long_run_neg, 4), collapse = ", "), "\n")
  
  invisible(x)
}


#' @export
summary.pnardl <- function(object, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("       Panel Nonlinear ARDL Estimation Results\n")
  cat("=======================================================\n\n")
  
  cat("Model Specification:\n")
  cat("  Estimator:", toupper(object$estimator), "\n")
  cat("  Effect:", object$effect, "\n")
  cat("  Panels (N):", object$N, "\n")
  cat("  Time periods (T):", object$T, "\n")
  cat("  Threshold:", object$threshold, "\n")
  cat("  Lags: p =", object$p, ", q =", paste(object$q, collapse = ","), "\n\n")
  
  cat("Error Correction Coefficient:\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  phi = %.4f", object$ec_coef))
  if (!is.null(object$ec_se)) {
    cat(sprintf(" (SE = %.4f)", object$ec_se))
  }
  cat("\n\n")
  
  cat("Long-Run Coefficients:\n")
  cat("-------------------------------------------------------\n")
  cat("  Variable        Positive    Negative    Asymmetric?\n")
  for (v in object$x_vars) {
    i <- which(object$x_vars == v)
    asym <- object$asymmetry_test[[v]]
    asym_mark <- if (!is.na(asym$asymmetric) && asym$asymmetric) "Yes ***" else "No"
    cat(sprintf("  %-14s %8.4f    %8.4f    %s\n",
                v, object$long_run_pos[i], object$long_run_neg[i], asym_mark))
  }
  
  cat("\nAsymmetry Tests (Wald):\n")
  cat("-------------------------------------------------------\n")
  for (v in names(object$asymmetry_test)) {
    test <- object$asymmetry_test[[v]]
    if (!is.na(test$wald)) {
      sig <- if (test$p_value < 0.01) "***" else if (test$p_value < 0.05) "**" else if (test$p_value < 0.1) "*" else ""
      cat(sprintf("  %-14s Wald = %6.3f, p = %.4f %s\n", v, test$wald, test$p_value, sig))
    }
  }
  
  if (!is.null(object$hausman)) {
    cat("\nHausman Test (PMG vs MG):\n")
    cat("-------------------------------------------------------\n")
    cat(sprintf("  H = %.3f, df = %d, p = %.4f\n",
                object$hausman$statistic, object$hausman$df, object$hausman$p_value))
    cat(" ", object$hausman$conclusion, "\n")
  }
  
  cat("\n=======================================================\n\n")
  
  invisible(object)
}
