
# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

#' @title Panel ARDL Estimation (PMG / MG / DFE)
#' @description Estimate Panel ARDL models with Pooled Mean Group (PMG),
#'   Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators following
#'   Pesaran, Shin & Smith (1999). This implementation replicates Stata's
#'   \code{xtpmg} command (Blackburne & Frank 2007).
#'
#' @details
#' The model is specified as:
#' \deqn{\Delta y_{it} = \phi_i (y_{i,t-1} - \theta'_i x_{it}) +
#'   \sum_{j=1}^{p-1} \lambda_{ij} \Delta y_{i,t-j} +
#'   \sum_{j=0}^{q-1} \delta'_{ij} \Delta x_{i,t-j} + \mu_i + \varepsilon_{it}}
#'
#' Three estimators are available:
#' \itemize{
#'   \item \strong{PMG}: Constrains long-run coefficients across groups while
#'         allowing short-run coefficients to vary. Efficient under long-run homogeneity.
#'   \item \strong{MG}: Estimates separate regressions per group and averages.
#'         Consistent but less efficient than PMG.
#'   \item \strong{DFE}: Traditional dynamic fixed effects with all coefficients pooled.
#'         Consistent only under slope homogeneity.
#' }
#'
#' @param formula A formula \code{y ~ x1 + x2 + ...}.
#' @param data A data frame containing panel data.
#' @param id Character. Name of the panel/group identifier column.
#' @param time Character. Name of the time variable column.
#' @param p Integer. Lags for dependent variable (default 1).
#' @param q Integer or numeric vector. Lags for independent variables (default 1).
#' @param estimator Character. One of \code{"pmg"}, \code{"mg"}, \code{"dfe"}.
#' @param ec Logical. If \code{TRUE}, display error-correction form (default \code{TRUE}).
#' @param trend Logical. Include a deterministic time trend (default \code{FALSE}).
#' @param start_time Optional. Restrict estimation to observations at or after this time.
#' @param cluster Logical. If \code{TRUE}, use cluster-robust standard errors where applicable (default \code{FALSE}).
#' @param maxiter Integer. Maximum iterations for PMG optimization (default 100).
#' @param tol Numeric. Convergence tolerance (default 1e-6).
#'
#' @return An object of class \code{"panel_ardl"} with components:
#' \itemize{
#'   \item \code{long_run}: long-run coefficients (theta).
#'   \item \code{short_run}: short-run coefficient summaries.
#'   \item \code{ec_coef}: error-correction coefficient (phi).
#'   \item \code{vcov}: variance-covariance matrix of long-run coefficients.
#'   \item \code{sigma}: residual standard error.
#'   \item \code{estimator}: the estimator used.
#' }
#'
#' @references
#' Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group
#' estimation of dynamic heterogeneous panels. \emph{Journal of the American
#' Statistical Association}, 94(446), 621-634.
#'
#' Blackburne, E. F., & Frank, M. W. (2007). Estimation of nonstationary
#' heterogeneous panels. \emph{Stata Journal}, 7(2), 197-208.
#'
#' @examples
#' \dontrun{
#' fit <- panel_ardl(y ~ x1 + x2, data = my_panel,
#'                   id = "country", time = "year",
#'                   p = 1, q = 1, estimator = "pmg")
#' summary(fit)
#' }
#'
#' @importFrom stats lm.fit pnorm pt pchisq nlminb na.omit model.matrix setNames
#' @importFrom MASS ginv
#' @export
panel_ardl <- function(formula, data, id, time, p = 1, q = 1,
                       estimator = c("pmg", "mg", "dfe"),
                       ec = TRUE, trend = FALSE,
                       start_time = NULL,
                       cluster = FALSE,
                       maxiter = 100, tol = 1e-6) {

  estimator <- match.arg(estimator)

  if (!is.data.frame(data))       stop("'data' must be a data frame")
  if (!id   %in% names(data))     stop(paste0("ID variable '",   id,   "' not found"))
  if (!time %in% names(data))     stop(paste0("Time variable '", time, "' not found"))

  formula_vars <- all.vars(formula)
  y_var  <- formula_vars[1]
  x_vars <- formula_vars[-1]

  if (!y_var %in% names(data)) stop(paste0("'", y_var, "' not found in data"))
  for (v in x_vars)
    if (!v %in% names(data)) stop(paste0("'", v, "' not found in data"))

  data   <- data[order(data[[id]], data[[time]]), ]
  groups <- unique(data[[id]])
  n_groups <- length(groups)

  if (length(q) == 1)
    q <- rep(q, length(x_vars))
  else if (length(q) != length(x_vars))
    stop("length(q) must equal number of x-vars or 1")

  result <- switch(estimator,
    pmg = .estimate_pmg(data, y_var, x_vars, id, time, p, q, trend, start_time, maxiter, tol),
    mg  = .estimate_mg (data, y_var, x_vars, id, time, p, q, trend, start_time),
    dfe = .estimate_dfe(data, y_var, x_vars, id, time, p, q, trend, start_time, cluster)
  )

  result$call      <- match.call()
  result$formula   <- formula
  result$estimator <- estimator
  result$trend     <- trend
  result$y_var     <- y_var
  result$x_vars    <- x_vars
  result$id        <- id
  result$time      <- time
  result$p         <- p
  result$q         <- q
  result$ngroups   <- n_groups
  result$groups    <- groups

  n <- result$nobs
  k <- length(result$long_run) + length(result$short_run)
  result$aic <- -2 * result$loglik + 2 * k
  result$bic <- -2 * result$loglik + log(n) * k

  class(result) <- c("panel_ardl", "list")
  result
}

# -----------------------------------------------------------------------------
# Data preparation  (single group, sorted by time)
#
# Returns list:
#   dy       - Dy_t
#   y_lag1   - y_{t-1}   (LRy in Stata sense: lagged level)
#   X_levels - X_t       (LRx: CURRENT-period levels for ECT)
#   X_diff   - [Dx_{jt}, ..., Dx_{j,t-qj+1}, Dy_{t-1}, ..., Dy_{t-p+1}]  (SR dynamics)
#   n        - number of usable obs after NA removal
# -----------------------------------------------------------------------------
.prepare_ardl_data <- function(data, y_var, x_vars, time, p, q, trend,
                               start_time = NULL) {

  data <- data[order(data[[time]]), ]
  n    <- nrow(data)

  if (n < max(p, max(q)) + 2) return(NULL)

  y      <- data[[y_var]]
  t_vals <- data[[time]]   # time values corresponding to each row

  # LHS: Dy_t  (length n-1, row k corresponds to period t_vals[k+1])
  dy <- diff(y)

  # ECT "y" component: y_{t-1}
  y_lag1 <- y[1:(n - 1)]

  # ECT "X" component: X at CURRENT period t -> rows 2:n
  X_levels <- as.matrix(data[2:n, x_vars, drop = FALSE])

  # Time vector aligned with dy (current period = t_vals[2:n])
  time_vec <- t_vals[2:n]

  # SR contemporaneous diffs  Dx_t
  X_diff <- do.call(cbind, lapply(x_vars, function(v) diff(data[[v]])))
  colnames(X_diff) <- paste0("D.", x_vars)

  # SR lagged diffs of x: lags 1 ... q[j]-1
  for (j in seq_along(x_vars)) {
    if (q[j] > 1) {
      xd <- diff(data[[x_vars[j]]])
      for (lag in seq_len(q[j] - 1)) {
        col <- c(rep(NA_real_, lag), xd[seq_len(length(xd) - lag)])
        X_diff <- cbind(X_diff, col)
        colnames(X_diff)[ncol(X_diff)] <- paste0("D.", x_vars[j], ".L", lag)
      }
    }
  }

  # SR lagged diffs of y: lags 1 ... p-1
  if (p > 1) {
    for (lag in seq_len(p - 1)) {
      col <- c(rep(NA_real_, lag), dy[seq_len(length(dy) - lag)])
      X_diff <- cbind(X_diff, col)
      colnames(X_diff)[ncol(X_diff)] <- paste0("D.", y_var, ".L", lag)
    }
  }

  # Combine and drop rows with any NA
  df <- data.frame(dy = dy, y_lag1 = y_lag1, .time = time_vec,
                   X_levels, X_diff,
                   check.names = FALSE)
  df <- na.omit(df)

  # Apply start_time filter AFTER lag computation (mirrors Stata's "if year>=X")
  if (!is.null(start_time)) {
    df <- df[df$.time >= start_time, , drop = FALSE]
  }
  df$.time <- NULL   # remove helper column

  if (nrow(df) < ncol(X_diff) + length(x_vars) + 3) return(NULL)

  xd_cols <- setdiff(names(df), c("dy", "y_lag1", x_vars))

  list(
    dy       = df[["dy"]],
    y_lag1   = df[["y_lag1"]],
    X_levels = as.matrix(df[, x_vars,   drop = FALSE]),
    X_diff   = as.matrix(df[, xd_cols,  drop = FALSE]),
    n        = nrow(df)
  )
}

# -----------------------------------------------------------------------------
# CalcMGE: mean and cross-sectional variance
#   V = B_demeaned' %*% B_demeaned / (N*(N-1))   (Stata CalcMGE formula)
# -----------------------------------------------------------------------------
.calc_mge <- function(B) {
  # B: matrix N x k  (each row = one group's coefficient vector)
  N    <- nrow(B)
  bbar <- colMeans(B, na.rm = TRUE)
  tmp  <- sweep(B, 2, bbar)
  V    <- crossprod(tmp) / (N * (N - 1))
  list(mean = bbar, vcov = V, se = sqrt(diag(V)))
}

# -----------------------------------------------------------------------------
# PMG estimator
# -----------------------------------------------------------------------------
.estimate_pmg <- function(data, y_var, x_vars, id, time, p, q, trend,
                          start_time = NULL, maxiter = 100, tol = 1e-6) {

  groups   <- unique(data[[id]])
  k_x      <- length(x_vars)

  # Prepare per-group data
  gdata_list <- lapply(groups, function(g) {
    .prepare_ardl_data(data[data[[id]] == g, ], y_var, x_vars, time, p, q, trend,
                       start_time)
  })
  ok <- !sapply(gdata_list, is.null)
  gdata_list <- gdata_list[ok]
  groups     <- groups[ok]
  N          <- length(groups)
  if (N < 2) stop("Need at least 2 valid groups")

  # -- Initial theta: pooled OLS of y_lag1 ~ X_levels (no constant) ----------
  # Mirrors Stata: regress $LRy $LRx, noconstant
  all_ylg  <- do.call(c,    lapply(gdata_list, `[[`, "y_lag1"))
  all_Xlv  <- do.call(rbind, lapply(gdata_list, `[[`, "X_levels"))
  theta0   <- drop(lm.fit(all_Xlv, all_ylg)$coefficients)
  names(theta0) <- x_vars

  # -- Concentrated log-likelihood --------------------------------------------
  # Given theta, for each group regress dy ~ ECT + X_diff + 1 (with constant)
  # and sum the MLE log-likelihoods.  Mirrors Stata's xtpmg_ml d0 evaluator.
  neg_ll <- function(theta) {
    ll <- 0
    for (gd in gdata_list) {
      ect  <- gd$y_lag1 - drop(gd$X_levels %*% theta)
      Xsr  <- cbind(ect, gd$X_diff, 1)
      fit  <- tryCatch(lm.fit(Xsr, gd$dy), error = function(e) NULL)
      if (is.null(fit)) { ll <- ll - 1e12; next }
      Ti   <- length(fit$residuals)
      s2i  <- sum(fit$residuals^2) / Ti          # MLE variance (/ T)
      if (!is.finite(s2i) || s2i <= 0) { ll <- ll - 1e12; next }
      ll   <- ll - 0.5 * Ti * (log(2 * pi) + log(s2i) + 1)
    }
    -ll   # return negative for minimisation
  }

  opt <- tryCatch(
    nlminb(theta0, neg_ll, control = list(iter.max = maxiter, rel.tol = tol)),
    error = function(e) list(par = theta0, convergence = 1)
  )
  theta_pmg <- setNames(opt$par, x_vars)

  # -- Group-specific SR regressions with estimated theta --------------------
  phi_vec    <- numeric(N)
  sig2_vec   <- numeric(N)
  XLR_list   <- vector("list", N)
  W_list     <- vector("list", N)   # [ECT, X_diff, 1] per group
  sr_list    <- vector("list", N)   # full SR coef vector per group
  res_list   <- vector("list", N)
  fit_list   <- vector("list", N)

  for (i in seq_len(N)) {
    gd       <- gdata_list[[i]]
    ect      <- gd$y_lag1 - drop(gd$X_levels %*% theta_pmg)
    Xsr      <- cbind(ect, gd$X_diff, 1)
    fit      <- lm.fit(Xsr, gd$dy)
    Ti       <- length(fit$residuals)

    phi_vec[i]  <- fit$coefficients[1]
    sig2_vec[i] <- sum(fit$residuals^2) / Ti
    XLR_list[[i]] <- gd$X_levels
    W_list[[i]]   <- Xsr
    sr_list[[i]]  <- fit$coefficients
    res_list[[i]] <- fit$residuals
    fit_list[[i]] <- fit$fitted.values
  }

  # -- PSS 1999 information matrix for theta SE ------------------------------
  # G = [ Gxx    Grow  ]    where Grow is kl * (N*ks)
  #     [ Grow'  G_SR  ]    G_SR block-diagonal  N * ks * ks blocks
  ks   <- ncol(W_list[[1]])   # = 1 (ec) + ncol(X_diff) + 1 (const)
  Gxx  <- matrix(0, k_x, k_x)
  Grow <- matrix(0, k_x, N * ks)
  Gsr  <- matrix(0, N * ks, N * ks)

  for (i in seq_len(N)) {
    phi_i <- phi_vec[i];  s2_i <- sig2_vec[i]
    XLRi  <- XLR_list[[i]];  Wi <- W_list[[i]]
    idx   <- seq((i - 1) * ks + 1, i * ks)

    Gxx          <- Gxx + (phi_i^2 / s2_i) * crossprod(XLRi)
    Grow[, idx]  <- -(phi_i / s2_i) * crossprod(XLRi, Wi)
    Gsr[idx,idx] <- (1 / s2_i) * crossprod(Wi)
  }

  G_full <- rbind(cbind(Gxx, Grow), cbind(t(Grow), Gsr))
  V_full <- tryCatch(solve(G_full),
                     error = function(e) MASS::ginv(G_full))
  V_theta  <- V_full[seq_len(k_x), seq_len(k_x), drop = FALSE]
  theta_se <- setNames(sqrt(diag(V_theta)), x_vars)
  theta_z  <- theta_pmg / theta_se
  theta_p  <- 2 * pnorm(-abs(theta_z))

  # -- Average SR coefficients (CalcMGE) -------------------------------------
  sr_mat  <- do.call(rbind, sr_list)
  sr_names <- c("ec", colnames(gdata_list[[1]]$X_diff), "_cons")
  colnames(sr_mat) <- sr_names
  mge     <- .calc_mge(sr_mat)
  sr_avg  <- mge$mean
  sr_se   <- mge$se

  # -- Log-likelihood (sum of group OLS MLE log-likelihoods) -----------------
  loglik <- -sum(sapply(seq_len(N), function(i) {
    Ti <- length(res_list[[i]])
    s2 <- sig2_vec[i]
    0.5 * Ti * (log(2 * pi) + log(s2) + 1)
  }))

  residuals_all <- do.call(c, res_list)
  fitted_all    <- do.call(c, fit_list)
  nobs          <- length(residuals_all)

  list(
    long_run      = theta_pmg,
    long_run_se   = theta_se,
    long_run_z    = theta_z,
    long_run_p    = theta_p,
    long_run_vcov = V_theta,
    short_run     = sr_avg,
    short_run_se  = sr_se,
    ec_coef       = unname(sr_avg["ec"]),
    ec_se         = unname(sr_se["ec"]),
    ec_group      = phi_vec,
    group_coefs   = sr_list,
    sigma2_group  = sig2_vec,
    residuals     = residuals_all,
    fitted        = fitted_all,
    nobs          = nobs,
    sigma         = sqrt(mean(residuals_all^2)),
    loglik        = loglik,
    convergence   = opt$convergence
  )
}

# -----------------------------------------------------------------------------
# MG estimator
# -----------------------------------------------------------------------------
.estimate_mg <- function(data, y_var, x_vars, id, time, p, q, trend,
                         start_time = NULL) {

  groups <- unique(data[[id]])
  k_x    <- length(x_vars)

  gdata_list <- lapply(groups, function(g) {
    .prepare_ardl_data(data[data[[id]] == g, ], y_var, x_vars, time, p, q, trend,
                       start_time)
  })
  ok         <- !sapply(gdata_list, is.null)
  gdata_list <- gdata_list[ok]
  groups     <- groups[ok]
  N          <- length(groups)
  if (N < 2) stop("Need at least 2 valid groups")

  # For each group: regress dy ~ y_lag1 + X_levels + X_diff + 1
  # Then transform: theta_j = -beta_j/phi,  phi = coef on y_lag1
  # Mirrors Stata: regress $SRy $LRy $LRx $SRx (with constant)
  group_results <- lapply(gdata_list, function(gd) {
    X_full <- cbind(gd$y_lag1, gd$X_levels, gd$X_diff, 1)
    tryCatch({
      fit   <- lm.fit(X_full, gd$dy)
      cf    <- fit$coefficients
      phi   <- cf[1]
      beta  <- cf[seq(2, k_x + 1)]
      theta <- -beta / phi
      sr    <- cf[seq(k_x + 2, length(cf))]   # [Dx cols, _cons]
      list(theta = theta, phi = phi, sr = sr,
           all = c(theta, phi, sr),
           residuals = fit$residuals, valid = TRUE)
    }, error = function(e) list(valid = FALSE))
  })

  valid          <- sapply(group_results, function(x) isTRUE(x$valid))
  group_results  <- group_results[valid]
  n_valid        <- length(group_results)

  # CalcMGE on full transformed vector [theta_1,...,theta_k, phi, SR...]
  B    <- do.call(rbind, lapply(group_results, `[[`, "all"))
  mge  <- .calc_mge(B)

  k_sr   <- length(group_results[[1]]$sr)
  all_names <- c(x_vars,
                 "ec",
                 colnames(gdata_list[[which(valid)[1]]]$X_diff),
                 "_cons")
  names(mge$mean) <- all_names
  names(mge$se)   <- all_names

  theta_mg <- mge$mean[seq_len(k_x)]
  theta_se <- mge$se  [seq_len(k_x)]
  theta_z  <- theta_mg / theta_se
  theta_p  <- 2 * pnorm(-abs(theta_z))

  phi_mg   <- mge$mean["ec"]
  phi_se   <- mge$se  ["ec"]

  sr_idx   <- seq(k_x + 1, k_x + 1 + k_sr)
  sr_avg   <- mge$mean[sr_idx]
  sr_se    <- mge$se  [sr_idx]

  # Residuals / fitted from full group regressions
  res_list <- lapply(seq_along(gdata_list)[valid], function(ii) {
    gd     <- gdata_list[[ii]]
    X_full <- cbind(gd$y_lag1, gd$X_levels, gd$X_diff, 1)
    fit    <- lm.fit(X_full, gd$dy)
    list(res = fit$residuals, fit = fit$fitted.values)
  })
  residuals_all <- do.call(c, lapply(res_list, `[[`, "res"))
  fitted_all    <- do.call(c, lapply(res_list, `[[`, "fit"))
  nobs          <- length(residuals_all)

  sig2   <- mean(residuals_all^2)
  loglik <- -0.5 * nobs * (log(2 * pi) + log(sig2) + 1)

  list(
    long_run      = theta_mg,
    long_run_se   = theta_se,
    long_run_z    = theta_z,
    long_run_p    = theta_p,
    long_run_vcov = mge$vcov[seq_len(k_x), seq_len(k_x)],
    short_run     = sr_avg,
    short_run_se  = sr_se,
    ec_coef       = unname(phi_mg),
    ec_se         = unname(phi_se),
    ec_group      = sapply(group_results, `[[`, "phi"),
    group_coefs   = lapply(group_results, `[[`, "all"),
    group_theta   = do.call(rbind, lapply(group_results, `[[`, "theta")),
    residuals     = residuals_all,
    fitted        = fitted_all,
    nobs          = nobs,
    sigma         = sqrt(sig2),
    loglik        = loglik
  )
}

# -----------------------------------------------------------------------------
# DFE estimator
# -----------------------------------------------------------------------------
.estimate_dfe <- function(data, y_var, x_vars, id, time, p, q, trend,
                          start_time = NULL, cluster = FALSE) {

  groups <- unique(data[[id]])
  k_x    <- length(x_vars)

  gdata_list <- lapply(groups, function(g) {
    gd <- .prepare_ardl_data(data[data[[id]] == g, ], y_var, x_vars, time, p, q, trend,
                             start_time)
    if (!is.null(gd)) gd$group <- g
    gd
  })
  ok         <- !sapply(gdata_list, is.null)
  gdata_list <- gdata_list[ok]
  N_valid    <- length(gdata_list)

  dy_all    <- do.call(c,    lapply(gdata_list, `[[`, "dy"))
  ylg_all   <- do.call(c,    lapply(gdata_list, `[[`, "y_lag1"))
  Xlv_all   <- do.call(rbind, lapply(gdata_list, `[[`, "X_levels"))
  Xdiff_all <- do.call(rbind, lapply(gdata_list, `[[`, "X_diff"))

  grp_vec   <- do.call(c, lapply(gdata_list, function(gd) rep(gd$group, gd$n)))
  grp_dum   <- model.matrix(~ factor(grp_vec) - 1)

  # Full regression: dy ~ y_lag1 + X_levels + X_diff + group_FE
  # Mirrors: xtreg $SRy $LRy $LRx $SRx, fe
  X_full <- cbind(ylg_all, Xlv_all, Xdiff_all, grp_dum)
  fit    <- lm.fit(X_full, dy_all)
  cf     <- fit$coefficients

  phi  <- cf[1]
  beta <- cf[seq(2, k_x + 1)]
  theta <- setNames(-beta / phi, x_vars)

  n     <- length(dy_all)
  n_xd  <- ncol(Xdiff_all)
  k_nfe <- 1L + k_x + n_xd          # phi + betas + deltas (no FE dummies)
  k_all <- ncol(X_full)
  G     <- N_valid

  XtXinv <- tryCatch(solve(crossprod(X_full)),
                     error = function(e) MASS::ginv(crossprod(X_full)))

  # Bread: top-left k_nfe * k_nfe block = (X_within'X_within)^{-1} by Frisch-Waugh
  B_slopes <- XtXinv[seq_len(k_nfe), seq_len(k_nfe)]
  X_slopes <- X_full[, seq_len(k_nfe), drop = FALSE]

  if (cluster) {
    # Cluster-robust sandwich (clustering by group = FE level)
    # Score for group g: s_g = X_slopes_g' * u_g  (within mean of u_g = 0)
    meat <- matrix(0, k_nfe, k_nfe)
    obs_start <- 1L
    for (gd in gdata_list) {
      Ti  <- gd$n
      idx <- seq_len(Ti) + obs_start - 1L
      sg  <- colSums(X_slopes[idx, , drop = FALSE] * fit$residuals[idx])
      meat <- meat + outer(sg, sg)
      obs_start <- obs_start + Ti
    }
    # Stata correction: G/(G-1) * (n-1)/(n-k_all) - k_all includes FE dummies
    correction  <- (G / (G - 1)) * ((n - 1) / (n - k_all))
    vcov_slopes <- B_slopes %*% meat %*% B_slopes * correction
  } else {
    s2          <- sum(fit$residuals^2) / (n - k_all)
    vcov_slopes <- s2 * B_slopes
  }

  # Full delta-method VCV for theta = -beta / phi
  # G_theta[j, 1]   = d(theta_j)/d(phi)   = beta[j] / phi^2
  # G_theta[j, 1+j] = d(theta_j)/d(beta_j) = -1 / phi
  G_theta <- matrix(0, k_x, k_nfe)
  for (j in seq_len(k_x)) {
    G_theta[j, 1]   <- beta[j] / phi^2
    G_theta[j, 1+j] <- -1 / phi
  }
  V_theta  <- G_theta %*% vcov_slopes %*% t(G_theta)
  theta_se <- setNames(sqrt(diag(V_theta)), x_vars)
  theta_z  <- theta / theta_se
  theta_p  <- 2 * pnorm(-abs(theta_z))

  # SR coefs and SEs
  ec_se_val    <- sqrt(vcov_slopes[1, 1])
  delta_idx    <- seq(k_x + 2L, k_nfe)
  sr_cf        <- cf[delta_idx]
  sr_se_diffs  <- sqrt(diag(vcov_slopes)[delta_idx])

  # _cons = grand-mean intercept (equivalent to Stata's xtreg fe _b[_cons])
  slope_means <- c(mean(ylg_all), colMeans(Xlv_all), colMeans(Xdiff_all))
  cons_val    <- mean(dy_all) - sum(cf[seq_len(k_nfe)] * slope_means)
  grad_cons   <- -slope_means
  se_cons     <- sqrt(as.numeric(t(grad_cons) %*% vcov_slopes %*% grad_cons))

  sr_names <- c("ec", colnames(Xdiff_all), "_cons")
  sr_avg   <- setNames(c(phi, sr_cf, cons_val), sr_names)
  sr_se_v  <- setNames(c(ec_se_val, sr_se_diffs, se_cons), sr_names)

  s2_mle <- sum(fit$residuals^2) / n
  loglik  <- -0.5 * n * (log(2 * pi) + log(s2_mle) + 1)

  list(
    long_run      = theta,
    long_run_se   = theta_se,
    long_run_z    = theta_z,
    long_run_p    = theta_p,
    long_run_vcov = V_theta,
    short_run     = sr_avg,
    short_run_se  = sr_se_v,
    ec_coef       = phi,
    ec_se         = ec_se_val,
    residuals     = fit$residuals,
    fitted        = fit$fitted.values,
    nobs          = n,
    sigma         = sqrt(s2_mle),
    loglik        = loglik
  )
}

# -----------------------------------------------------------------------------
# print / summary
# -----------------------------------------------------------------------------

#' @title Print Method for panel_ardl Objects
#' @description Concise display of a fitted \code{panel_ardl} model.
#' @param x A \code{panel_ardl} object.
#' @param ... Additional arguments (currently ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.panel_ardl <- function(x, ...) {
  cat("\nPanel ARDL (", toupper(x$estimator), ") - ARDL(",
      x$p, ",", paste(x$q, collapse = ","), ")\n", sep = "")
  cat("Groups:", x$ngroups, " | Obs:", x$nobs,
      " | Log-lik:", round(x$loglik, 3), "\n\n")
  cat("Long-run:\n"); print(round(x$long_run, 6))
  cat("\nError-correction (avg phi):", round(x$ec_coef, 6), "\n")
  invisible(x)
}

#' @title Summary Method for panel_ardl Objects
#' @description Detailed summary of a fitted \code{panel_ardl} model with
#'   long-run, short-run, and error-correction coefficients.
#' @param object A \code{panel_ardl} object.
#' @param digits Integer. Decimal digits to display (default 6).
#' @param ... Additional arguments (currently ignored).
#' @return Invisibly returns the summary as a list.
#' @export
summary.panel_ardl <- function(object, digits = 6, ...) {

  est_label <- switch(object$estimator,
    pmg = "Pooled Mean Group Regression",
    mg  = "Mean Group Estimation: Error Correction Form",
    dfe = "Dynamic Fixed Effects Regression: Error Correction Form"
  )

  cat("\n")
  cat(est_label, "\n")
  cat(rep("-", 78), "\n", sep = "")
  cat(sprintf("%-30s  Number of obs    = %10d\n",
              paste0("Panel var (i): ", object$id),   object$nobs))
  cat(sprintf("%-30s  Number of groups = %10d\n",
              paste0("Time  var (t): ", object$time),  object$ngroups))
  if (!is.null(object$loglik))
    cat(sprintf("%-30s  Log Likelihood   = %10.3f\n", "", object$loglik))
  cat(rep("-", 78), "\n", sep = "")

  # table helper
  .print_block <- function(label, coefs, se, zv, pv) {
    cat(sprintf("%-10s |  %12s  %10s  %7s  %6s\n",
                label, "Coef.", "Std. Err.", "z", "P>|z|"))
    cat(rep("-", 78), "\n", sep = "")
    for (nm in names(coefs)) {
      stars <- ifelse(pv[nm] < 0.01, "**",
               ifelse(pv[nm] < 0.05, "* ", "  "))
      cat(sprintf("  %-8s |  %12.7f  %10.7f  %7.2f  %6.4f %s\n",
                  nm, coefs[nm], se[nm], zv[nm], pv[nm], stars))
    }
    cat(rep("-", 78), "\n", sep = "")
  }

  # Long-run block (ec equation)
  lr_z <- object$long_run / object$long_run_se
  lr_p <- 2 * pnorm(-abs(lr_z))
  cat(sprintf("%-10s |\n", "ec"))
  .print_block("", object$long_run, object$long_run_se, lr_z, lr_p)

  # Short-run block
  sr_z <- object$short_run / object$short_run_se
  sr_p <- 2 * pnorm(-abs(sr_z))
  cat(sprintf("%-10s |\n", "SR"))
  .print_block("", object$short_run, object$short_run_se, sr_z, sr_p)

  cat("Signif.: ** p<0.01  * p<0.05\n")
  invisible(object)
}

# -----------------------------------------------------------------------------
# Hausman test
# -----------------------------------------------------------------------------

#' @title Hausman Test for Panel ARDL Estimators
#' @description Compare two \code{panel_ardl} estimators (typically MG vs PMG)
#'   to test for long-run slope homogeneity.
#'
#' @details
#' Under the null of long-run homogeneity, the PMG (efficient) estimator is
#' consistent and efficient; the MG (inefficient) estimator is consistent but
#' less efficient. A significant test statistic rejects PMG in favour of MG.
#'
#' Argument order in v2.0.0 follows Stata's convention: pass the
#' \strong{inefficient} estimator first (e.g. MG), then the \strong{efficient}
#' one (e.g. PMG). This is a breaking change from earlier versions.
#'
#' @param inefficient A \code{panel_ardl} object - the always-consistent
#'   estimator (typically MG).
#' @param efficient A \code{panel_ardl} object - the under-H0 efficient
#'   estimator (typically PMG).
#' @param sigmamore Logical. If \code{TRUE} (default, matches Stata's
#'   \code{sigmamore} option), rescale the inefficient variance matrix to
#'   address non-positive-definite difference matrices.
#'
#' @return An (invisible) list with \code{statistic}, \code{df},
#'   \code{p.value}, and \code{theta_diff}.
#'
#' @references
#' Hausman, J. A. (1978). Specification tests in econometrics.
#' \emph{Econometrica}, 46(6), 1251-1271.
#'
#' @examples
#' \dontrun{
#' fit_mg  <- panel_ardl(y ~ x, data = d, id = "i", time = "t", estimator = "mg")
#' fit_pmg <- panel_ardl(y ~ x, data = d, id = "i", time = "t", estimator = "pmg")
#' hausman_test(fit_mg, fit_pmg)
#' }
#'
#' @export
hausman_test <- function(inefficient, efficient, sigmamore = TRUE) {

  if (inefficient$estimator != "mg")
    stop("'inefficient' must be a MG model")
  if (!efficient$estimator %in% c("pmg", "dfe"))
    stop("'efficient' must be a PMG or DFE model")

  theta_diff <- inefficient$long_run - efficient$long_run
  df <- length(theta_diff)

  V_ineff <- inefficient$long_run_vcov
  V_eff   <- efficient$long_run_vcov

  if (sigmamore) {
    # xtpmg stores e(sigma) = mean(resid^2) = sigma2, while Stata's hausman
    # assumes e(sigma) = sigma (RMSE) and does scale = (e(sigma)_E/e(sigma)_b)^2.
    # Net effect: scale = (sigma2_eff / sigma2_ineff)^2 = (sigma_eff/sigma_ineff)^4.
    scale   <- (efficient$sigma / inefficient$sigma)^4
    V_ineff <- scale * V_ineff
  }

  W   <- V_ineff - V_eff
  W   <- (W + t(W)) / 2          # symmetrise numerically
  eig <- eigen(W, symmetric = TRUE, only.values = TRUE)$values

  if (any(eig <= 0)) {
    H  <- 0
    pv <- 1
  } else {
    H  <- as.numeric(t(theta_diff) %*% solve(W) %*% theta_diff)
    pv <- 1 - pchisq(H, df)
  }

  eff_label <- toupper(efficient$estimator)
  cat(sprintf("\n--- Hausman Test (sigmamore=%s): MG vs %s ---\n",
              sigmamore, eff_label))
  cat("H0: long-run homogeneity (efficient estimator preferred)\n")
  cat(sprintf("Chi2(%d) = %.4f   p-value = %.4f\n", df, H, pv))
  cat(if (pv < 0.05) "Reject H0: use MG\n" else "Cannot reject H0: efficient model preferred\n")

  invisible(list(statistic = H, df = df, p.value = pv, theta_diff = theta_diff))
}