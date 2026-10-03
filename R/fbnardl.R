#' @title Fourier (Bootstrap) Nonlinear ARDL
#' @description Nonlinear ARDL model with positive and negative partial sums
#'   (Shin, Yu and Greenwood-Nimmo, 2014), optional Fourier terms for smooth
#'   breaks, data-based selection of the Fourier frequency and of the lags,
#'   and cointegration tests with Kripfganz and Schneider (2020) bounds (no
#'   Fourier terms) or bootstrap critical values.
#'
#' @details
#' The model (case 3, unrestricted intercept) is
#' \deqn{\Delta y_t = c + \rho y_{t-1} + \sum_i (\beta_i^+ x^+_{i,t-1} +
#'   \beta_i^- x^-_{i,t-1}) + \sum_j \delta_j z_{j,t-1} +
#'   \sum_{l=1}^{p}\phi_l \Delta y_{t-l} +
#'   \sum_i \sum_{l=0}^{q_i} (\theta^+_{il}\Delta x^+_{i,t-l} +
#'   \theta^-_{il}\Delta x^-_{i,t-l}) +
#'   \sum_j \sum_{l=0}^{r_j}\pi_{jl}\Delta z_{j,t-l} + \gamma' e_t +
#'   a_1\sin(2\pi k^* t/T) + a_2\cos(2\pi k^* t/T) + u_t,}
#' where \eqn{x^+_t = \sum_{s \le t}\max(\Delta x_s, 0)} and
#' \eqn{x^-_t = \sum_{s \le t}\min(\Delta x_s, 0)} start at zero at the
#' first observation, \eqn{z} are the regressors that are not decomposed and
#' \eqn{e_t} the exogenous short-run variables (\code{exog}).
#'
#' \strong{Selection.} When \code{fourier = TRUE} and \code{kstar} is not
#' given, \eqn{k^*} minimises the sum of squared residuals of the model with
#' every lag at \code{maxlag}; the grid is \eqn{1, \dots,} \code{maxk}
#' (\code{kgrid = "integer"}; integer frequencies chosen by minimum SSR as in
#' Enders and Lee, 2012, Oxford Bulletin of Economics and Statistics) or \eqn{0.1, 0.2, \dots,} \code{maxk} (\code{kgrid =
#' "fractional"}). The default and the validity of the fractional grid are
#' package choices pending verification against Yilanci, Bozoklu and Gorus
#' (2020) and Omay (2015). Then p (1 to \code{maxlag}) and each \eqn{q_i},
#' \eqn{r_j} (0 to \code{maxlag}) minimise the AIC or BIC. Every candidate
#' and the final model use the same sample, \eqn{t > } \code{maxlag} + 1, so
#' the criteria are comparable. \code{lags} fixes the lags instead.
#'
#' \strong{Cointegration tests.} Fov is the F test that \eqn{\rho}, all
#' \eqn{\beta} and \eqn{\delta} are zero, t the t test on \eqn{\rho} and Find
#' the F test on the \eqn{\beta} and \eqn{\delta}. Without Fourier terms
#' (\code{fourier = FALSE}) and \code{type = "fnardl"}, the bounds of
#' Kripfganz and Schneider (2020) are used for Fov and t, with k = 2 times
#' the number of decomposed variables plus the number of other regressors,
#' n the number of observations used and the number of short-run
#' coefficients sr = p + sum 2(q_i + 1) + sum (r_j + 1) + number of
#' \code{exog}; the bounds decision requires both statistics beyond their
#' I(1) bounds. With Fourier terms no valid bounds exist: no decision is
#' given for \code{type = "fnardl"} and the bounds printed are those of the
#' model without Fourier terms, for reference only. Inference with Fourier
#' terms is by the bootstrap (\code{type = "fbnardl"}).
#'
#' \strong{Bootstrap.} \code{type = "fbnardl"} uses the recursive bootstrap
#' engine of \code{\link{boot_ardl}}: \code{bootstrap = "bvz"} has separate
#' nulls for Fov, t and Find (Bertelli, Vacca and Zoia, 2022), x* from the
#' marginal model of \eqn{\Delta x_t} on the intercept, \eqn{x_{t-1}}, p
#' lags of \eqn{\Delta y} and \eqn{\Delta x}, the Fourier terms and
#' \code{exog} (\code{xdgp = "vecm"}, default; BVZ eqs. 19-20), or from the
#' regression of \eqn{\Delta x_t} on the intercept, the Fourier terms and
#' \code{exog} only (\code{xdgp = "rw"}, a random walk with these
#' deterministic terms; a package choice that is not the x model of BVZ),
#' residuals recentred after each draw and a random block of initial values
#' (the partial sums continue the data's levels in the block);
#' \code{bootstrap = "mcnown"} uses the null of the Fov test for all
#' statistics (McNown, Sam and Goh, 2018, Steps 1-8) applied to the
#' conditional ECM, which is a package choice (MSG write the y equation
#' without the contemporaneous differences), with the unrestricted equation
#' for \eqn{\Delta x} that also contains \eqn{y_{t-1}}, residuals recentred
#' once and the observed initial values.
#' The partial sums are rebuilt from x*; the Fourier terms and \code{exog}
#' are held fixed. On every bootstrap sample \eqn{k^*} (unless fixed by
#' \code{kstar}) and the lags (unless fixed by \code{lags}) are selected
#' again, so the bootstrap distribution accounts for the selection
#' (package choice; BVZ step 5(c) re-estimates the unrestricted model but does not discuss re-selection). The
#' decision requires Fov, t and Find to reject; the degenerate cases are
#' labelled. Failed replications are \code{NA} and counted.
#' In a Monte Carlo with independent random walks (n = 100, 200 samples, B =
#' 199, maxlag = 1, integer grid) the 5\% rejection rates were 7.0\% (Monte
#' Carlo s.e. 1.8\%) for Fov, 11.5\% (2.3\%) for t and 11.5\% (2.3\%) for
#' Find, and 3.5\% (1.3\%) for the combined decision. Find tends to
#' over-reject with Fourier terms and partial sums (16.0\%, s.e. 2.6\%, for
#' \code{aardl(type = "fbnardl")} in the same design).
#'
#' \strong{Asymmetry and multipliers.} The long-run coefficients
#' \eqn{-\beta^\pm/\rho} have delta-method standard errors, and the long-run
#' symmetry test is the delta-method Wald test of \eqn{-\beta^+/\rho =
#' -\beta^-/\rho} (chi-squared, 1 df). The short-run symmetry test compares
#' the sums over lags 0 to \eqn{q_i}, \eqn{\sum_l\theta^+_{il} =
#' \sum_l\theta^-_{il}} (F(1, df)); testing the sums is a package choice
#' pending verification against Shin, Yu and Greenwood-Nimmo (2014). The
#' cumulative dynamic multipliers are the responses of the level of y to a
#' permanent unit increase in \eqn{x^+} or \eqn{x^-} at horizon 0, computed
#' by the full recursion of the estimated ECM (all lagged differences
#' included), with parametric confidence bands from \code{bands} draws of
#' the coefficients from their estimated normal distribution. The sign
#' convention of the multiplier of \eqn{x^-} (response to a unit increase in
#' the negative partial sum) is pending verification against Shin, Yu and
#' Greenwood-Nimmo (2014).
#'
#' The R code is a port written for this package, following the Stata
#' module fbnardl.
#'
#' @param formula A formula \code{y ~ x1 + x2 + ...}; only plain variable
#'   names are accepted.
#' @param data A data frame with consecutive observations (no missing
#'   values or gaps).
#' @param decompose Names of the regressors decomposed into positive and
#'   negative partial sums.
#' @param exog Optional names of exogenous variables in \code{data} that
#'   enter at time t in the short run only (for example dummies); held fixed
#'   in the bootstrap.
#' @param type \code{"fnardl"} (bounds) or \code{"fbnardl"} (bootstrap).
#' @param maxlag Largest lag considered (default 4).
#' @param maxk Largest Fourier frequency (default 3).
#' @param kgrid \code{"integer"} (default) or \code{"fractional"}.
#' @param fourier Logical; include Fourier terms (default \code{TRUE}).
#' @param kstar Optional Fourier frequency fixed by the user; it is then not
#'   selected, neither on the data nor in the bootstrap.
#' @param lags Optional list with elements \code{p}, \code{q} (one per
#'   decomposed variable) and \code{r} (one per other regressor) fixing the
#'   lags.
#' @param ic \code{"aic"} or \code{"bic"}.
#' @param bootstrap \code{"bvz"} (default) or \code{"mcnown"}.
#' @param xdgp \code{"vecm"} (default) or \code{"rw"}, the model of x* with
#'   \code{bootstrap = "bvz"}.
#' @param reps Bootstrap replications (default 999).
#' @param vcov Covariance matrix for the long-run standard errors,
#'   asymmetry tests and multiplier bands: \code{"ols"}, \code{"HC1"} or
#'   \code{"HAC"} (Newey-West with Bartlett weights and
#'   \eqn{\lfloor 4 (n/100)^{2/9} \rfloor} lags). The bounds and bootstrap
#'   statistics always use the OLS covariance.
#' @param horizon Largest horizon of the dynamic multipliers (default 20).
#' @param bands Number of parametric draws for the multiplier bands
#'   (default 500; 0 for none).
#' @param level Confidence level of the bands (default 0.95).
#' @param seed Optional seed for the bootstrap and the bands; the
#'   session's random number generator state is restored afterwards.
#' @param graph Logical; plot the multipliers (default \code{FALSE}).
#'
#' @return An object of class \code{"fbnardl"} with elements
#'   \code{coefficients}, \code{vcov}, \code{se}, \code{residuals},
#'   \code{fitted.values}, \code{lags}, \code{kstar}, \code{ic_value},
#'   \code{ssr_by_k}, \code{multipliers} (short- and long-run, with
#'   standard errors), \code{asymmetry_tests}, \code{bounds_test} (statistics,
#'   bounds or bootstrap results and decision), \code{dynamic_multipliers}
#'   (paths and bands), \code{diagnostics}, \code{nobs}, \code{df.residual},
#'   \code{r.squared}, \code{adj.r.squared}, \code{sigma}, \code{design},
#'   \code{call} and the settings.
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' Enders, W. and Lee, J. (2012). A unit root test using a Fourier series
#' to approximate smooth breaks. \emph{Oxford Bulletin of Economics and
#' Statistics}, 74(4), 574-599. \doi{10.1111/j.1468-0084.2011.00662.x}
#'
#' Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Omay, T. (2015). Fractional frequency flexible Fourier form to
#' approximate smooth breaks in unit root testing. \emph{Economics Letters},
#' 134, 123-126. \doi{10.1016/j.econlet.2015.07.010}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In R. C. Sickles and W. C. Horrace (Eds.), \emph{Festschrift in Honor of
#' Peter Schmidt} (pp. 281-314). Springer. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' Yilanci, V., Bozoklu, S. and Gorus, M. S. (2020). Are BRICS countries
#' pollution havens? Evidence from a bootstrap ARDL bounds testing approach
#' with a Fourier function. \emph{Sustainable Cities and Society}, 55,
#' 102035. \doi{10.1016/j.scs.2020.102035}
#'
#' @examples
#' d <- generate_oil_data(n = 100)
#' m <- fbnardl(gasoline ~ oil_price, data = d, decompose = "oil_price",
#'              maxlag = 2, fourier = FALSE)
#' m
#' \donttest{
#' mb <- fbnardl(gasoline ~ oil_price, data = d, decompose = "oil_price",
#'               type = "fbnardl", maxlag = 1, maxk = 2, reps = 49, bands = 50,
#'               seed = 1)
#' summary(mb)
#' }
#'
#' @export
fbnardl <- function(formula, data, decompose, exog = NULL,
                    type = c("fnardl", "fbnardl"), maxlag = 4, maxk = 3,
                    kgrid = c("integer", "fractional"), fourier = TRUE,
                    kstar = NULL, lags = NULL, ic = c("aic", "bic"),
                    bootstrap = c("bvz", "mcnown"), xdgp = c("vecm", "rw"),
                    reps = 999, vcov = c("ols", "HC1", "HAC"), horizon = 20,
                    bands = 500, level = 0.95, seed = NULL, graph = FALSE) {

  type <- match.arg(type)
  kgrid <- match.arg(kgrid)
  ic <- match.arg(ic)
  bootstrap <- match.arg(bootstrap)
  xdgp <- match.arg(xdgp)
  vcov_type <- match.arg(vcov)
  call <- match.call()

  fv <- .ardl_formula_vars(formula, data)
  depvar <- fv$y_var
  regs <- fv$x_vars
  if (missing(decompose) || length(decompose) == 0)
    stop("'decompose' must name at least one regressor")
  if (!all(decompose %in% regs))
    stop("all variables in 'decompose' must be regressors in the formula")
  decompose <- regs[regs %in% decompose]
  ctrl <- setdiff(regs, decompose)
  if (maxlag < 1) stop("'maxlag' must be at least 1")
  maxlag <- as.integer(maxlag)
  if (fourier && is.null(kstar)) {
    if (kgrid == "integer" && (maxk < 1 || maxk != round(maxk)))
      stop("'maxk' must be a positive integer with kgrid = 'integer'")
    if (kgrid == "fractional" && maxk < 0.1)
      stop("'maxk' must be at least 0.1 with kgrid = 'fractional'")
  }
  if (!is.null(kstar) && (!fourier || kstar <= 0))
    stop("'kstar' must be positive and requires fourier = TRUE")
  E <- NULL
  if (!is.null(exog)) {
    miss <- setdiff(exog, names(data))
    if (length(miss)) stop("exog variable(s) not found in 'data': ", paste(miss, collapse = ", "))
    if (any(exog %in% c(depvar, regs))) stop("'exog' variables must not appear in the formula")
    E <- as.matrix(data[, exog, drop = FALSE])
    if (!is.numeric(E) || anyNA(E)) stop("'exog' variables must be numeric without missing values")
  }
  y <- data[[depvar]]
  X <- as.matrix(data[, c(decompose, ctrl), drop = FALSE])
  n <- length(y)
  if (n - maxlag - 1 < 20) stop("too few observations for 'maxlag'")

  setup <- list(dec = decompose, ctrl = ctrl, E = E, maxlag = maxlag, n = n,
                fourier = fourier, ic = ic,
                kvalues = if (kgrid == "integer") seq_len(maxk) else round(seq(0.1, maxk, by = 0.1), 10))
  if (!is.null(lags)) lags <- .fb_check_lags(lags, setup)
  fix <- list(kstar = kstar, lags = lags)

  sel <- .fb_select(y, X, setup, fix)
  spec <- sel[c("p", "q", "r", "kstar")]
  B <- .fb_build(y, X, setup, spec)
  Z <- B$M
  dy <- B$Y
  N <- length(dy)
  K <- ncol(Z)
  qz <- qr(Z)
  if (qz$rank < K) stop("the selected model has a rank-deficient design")
  b <- qr.coef(qz, dy)
  e <- qr.resid(qz, dy)
  XtXi <- chol2inv(qr.R(qz))[order(qz$pivot), order(qz$pivot)]
  dimnames(XtXi) <- list(colnames(Z), colnames(Z))
  V <- .fb_vcov(Z, e, vcov_type, XtXi)
  V_ols <- sum(e^2) / (N - K) * XtXi
  dimnames(V) <- dimnames(V_ols) <- list(colnames(Z), colnames(Z))
  names(b) <- colnames(Z)

  # cointegration statistics (OLS)
  st <- .fb_stats(Z, dy, setup)
  kb <- 2 * length(decompose) + length(ctrl)
  sr <- spec$p + sum(2 * (spec$q + 1)) + sum(spec$r + 1) + (if (is.null(E)) 0 else ncol(E))
  bt <- list(Fov = unname(st["Fov"]), t = unname(st["t"]), Find = unname(st["Find"]),
             k = kb, n = N, sr = sr)
  ksF <- .ks_bounds("F", 3, kb, N, sr, siglevels = c(10, 5, 2.5, 1),
                    value = if (fourier) NULL else bt$Fov)
  kst <- .ks_bounds("t", 3, kb, N, sr, siglevels = c(10, 5, 2.5, 1),
                    value = if (fourier) NULL else bt$t)
  bt$bounds_F <- ksF$cv
  bt$bounds_t <- kst$cv
  bt$bounds_valid <- !fourier
  if (!fourier) {
    bt$pvalue_F <- ksF$pvalue
    bt$pvalue_t <- kst$pvalue
  }

  if (type == "fbnardl") {
    selfun <- NULL
    need_k <- fourier && is.null(kstar)
    need_l <- is.null(lags)
    if (need_k || need_l) {
      selfun <- function(yy, xx) .fb_select(yy, xx, setup, fix)[c("p", "q", "r", "kstar")]
    }
    stat_fun <- function(yy, xx, sp) {
      bb <- .fb_build(yy, xx, setup, sp)
      .fb_stats(bb$M, bb$Y, setup)
    }
    det <- cbind(E, if (spec$kstar > 0) .fb_fourier(n, spec$kstar))
    qc <- c(rep(spec$q, each = 2), spec$r)
    if (bootstrap == "bvz") {
      eng <- .ardl_boot_engine(y, X, p = spec$p, q = qc, case = 3,
                               t0 = maxlag + 2, det = det,
                               decompose = .fb_decomposer(setup),
                               nulls = "separate", xmodel = xdgp, B = reps,
                               init = "block", recentre = "draw",
                               stat_fun = stat_fun, select = selfun, spec = spec,
                               use_find = TRUE, seed = seed)
    } else {
      eng <- .ardl_boot_engine(y, X, p = spec$p, q = qc, case = 3,
                               t0 = maxlag + 2, det = det,
                               decompose = .fb_decomposer(setup),
                               nulls = "joint", xmodel = "var", B = reps,
                               init = "observed", recentre = "once",
                               stat_fun = stat_fun, select = selfun, spec = spec,
                               use_find = TRUE, seed = seed)
    }
    bt$bootstrap <- eng
    bt$pvalue_boot <- eng$p_value
    bt$cv_boot <- eng$cv
    bt$decision <- eng$decision
    bt$message <- paste0(eng$label, " (bootstrap, ", 100 * eng$level, "% level)")
    bt$method <- "bootstrap"
  } else if (fourier) {
    bt$decision <- "NOT_AVAILABLE"
    bt$message <- paste("No valid bounds exist with Fourier terms; no decision.",
                        "Use type = 'fbnardl'. The bounds shown are for the",
                        "model without Fourier terms; not valid with Fourier terms.")
    bt$method <- "none"
  } else {
    Fb <- ksF$cv[, "5%"]
    tb <- kst$cv[, "5%"]
    if (anyNA(c(Fb, tb))) {
      bt$decision <- "INCONCLUSIVE"
      bt$message <- "Bounds unavailable (too few degrees of freedom)"
    } else if (bt$Fov > Fb["I1"] && bt$t < tb["I1"]) {
      bt$decision <- "COINTEGRATION"
      bt$message <- "Cointegration: Fov > I(1) bound and t < I(1) bound (5%)"
    } else if (bt$Fov < Fb["I0"] || bt$t > tb["I0"]) {
      bt$decision <- "NO_COINTEGRATION"
      bt$message <- "No cointegration: Fov < I(0) bound or t > I(0) bound (5%)"
    } else {
      bt$decision <- "INCONCLUSIVE"
      bt$message <- "Inconclusive: Fov or t between the I(0) and I(1) bounds (5%)"
    }
    bt$method <- "Kripfganz-Schneider bounds"
  }

  mult <- .fb_multipliers(b, V, setup, spec, N - K)
  asym <- .fb_asymmetry(b, V, setup, spec, N - K)
  dmult <- .with_local_seed(seed, .fb_dynamic_multipliers(b, V, setup, spec, horizon, bands, level))
  diag <- .fb_diagnostics(Z, dy, e)

  rss <- sum(e^2)
  tss <- sum((dy - mean(dy))^2)
  result <- list(
    coefficients = b, vcov = V, vcov_ols = V_ols, se = sqrt(diag(V)),
    residuals = e, fitted.values = dy - e,
    lags = list(p = spec$p, q = stats::setNames(spec$q, decompose),
                r = stats::setNames(spec$r, ctrl)),
    kstar = spec$kstar, ic_value = sel$ic, ssr_by_k = sel$ssr_by_k,
    multipliers = mult, asymmetry_tests = asym, bounds_test = bt,
    dynamic_multipliers = dmult, diagnostics = diag,
    nobs = N, df.residual = N - K,
    r.squared = 1 - rss / tss,
    adj.r.squared = 1 - (rss / (N - K)) / (tss / (N - 1)),
    sigma = sqrt(rss / (N - K)), design = Z, response = dy,
    call = call, formula = formula, depvar = depvar,
    decomposed_vars = decompose, control_vars = ctrl, exog = exog,
    type = type, fourier = fourier, kgrid = kgrid, maxk = maxk,
    maxlag = maxlag, ic = ic, vcov_type = vcov_type, level = level,
    bootstrap_scheme = if (type == "fbnardl") bootstrap else NULL,
    xdgp = if (type == "fbnardl" && bootstrap == "bvz") xdgp else NULL,
    kstar_fixed = !is.null(kstar)
  )
  class(result) <- "fbnardl"
  if (graph) plot.fbnardl(result)
  result
}


.with_local_seed <- function(seed, expr) .abe_with_seed(seed, expr)

.fb_check_lags <- function(lags, setup) {
  if (!is.list(lags) || is.null(lags$p)) stop("'lags' must be a list with p, q and r")
  q <- if (is.null(lags$q)) integer(0) else as.integer(lags$q)
  r <- if (is.null(lags$r)) integer(0) else as.integer(lags$r)
  if (length(q) == 1 && length(setup$dec) > 1) q <- rep(q, length(setup$dec))
  if (length(r) == 1 && length(setup$ctrl) > 1) r <- rep(r, length(setup$ctrl))
  if (length(q) != length(setup$dec)) stop("'lags$q' needs one entry per decomposed variable")
  if (length(r) != length(setup$ctrl)) stop("'lags$r' needs one entry per other regressor")
  p <- as.integer(lags$p)
  if (p < 1 || any(c(p, q, r) > setup$maxlag) || any(c(q, r) < 0))
    stop("lags must satisfy 1 <= p <= maxlag and 0 <= q, r <= maxlag")
  list(p = p, q = q, r = r)
}

.fb_fourier <- function(n, k) {
  tt <- seq_len(n)
  cbind(fourier_sin = sin(2 * pi * k * tt / n), fourier_cos = cos(2 * pi * k * tt / n))
}

# Level regressors: x_pos, x_neg for each decomposed variable, then the others
.fb_levels <- function(X, setup) {
  W <- NULL
  nm <- character(0)
  for (j in seq_along(setup$dec)) {
    d <- c(0, diff(X[, j]))
    W <- cbind(W, cumsum(d * (d > 0)), cumsum(d * (d < 0)))
    nm <- c(nm, paste0(setup$dec[j], c("_pos", "_neg")))
  }
  for (j in seq_along(setup$ctrl)) {
    W <- cbind(W, X[, length(setup$dec) + j])
    nm <- c(nm, setup$ctrl[j])
  }
  colnames(W) <- nm
  W
}

.fb_decomposer <- function(setup) {
  nd <- length(setup$dec)
  nc <- length(setup$ctrl)
  list(full = function(x) .fb_levels(as.matrix(x), setup),
       step = function(d) {
         out <- numeric(2 * nd + nc)
         for (j in seq_len(nd)) out[2 * j - 1:0] <- c(d[j] * (d[j] > 0), d[j] * (d[j] < 0))
         if (nc) out[2 * nd + seq_len(nc)] <- d[nd + seq_len(nc)]
         out
       })
}

# All candidate columns up to maxlag on the common sample t > maxlag + 1
.fb_full <- function(y, X, setup) {
  n <- length(y)
  rows <- (setup$maxlag + 2):n
  W <- .fb_levels(X, setup)
  dy <- c(NA, diff(y))
  dW <- rbind(NA, diff(W))
  L <- function(v, j) if (j == 0) v else c(rep(NA, j), v[seq_len(n - j)])
  cols <- list(const = rep(1, n))
  for (j in seq_len(setup$maxlag)) cols[[paste0("dy_L", j)]] <- L(dy, j)
  for (m in colnames(W)) for (j in 0:setup$maxlag)
    cols[[if (j == 0) paste0("d_", m) else paste0("d_", m, "_L", j)]] <- L(dW[, m], j)
  cols[["L_y"]] <- L(y, 1)
  for (m in colnames(W)) cols[[paste0("L_", m)]] <- L(W[, m], 1)
  if (!is.null(setup$E)) for (m in colnames(setup$E)) cols[[m]] <- setup$E[, m]
  M <- do.call(cbind, cols)[rows, , drop = FALSE]
  list(M = M, Y = dy[rows], rows = rows, levels = colnames(W))
}

.fb_colnames <- function(spec, setup) {
  nm <- c("const", if (spec$p > 0) paste0("dy_L", seq_len(spec$p)))
  for (i in seq_along(setup$dec)) for (s in c("_pos", "_neg")) {
    m <- paste0(setup$dec[i], s)
    nm <- c(nm, paste0("d_", m), if (spec$q[i] > 0) paste0("d_", m, "_L", seq_len(spec$q[i])))
  }
  for (j in seq_along(setup$ctrl)) {
    m <- setup$ctrl[j]
    nm <- c(nm, paste0("d_", m), if (spec$r[j] > 0) paste0("d_", m, "_L", seq_len(spec$r[j])))
  }
  nm <- c(nm, "L_y", paste0("L_", c(rbind(paste0(setup$dec, "_pos"), paste0(setup$dec, "_neg")))),
          if (length(setup$ctrl)) paste0("L_", setup$ctrl),
          if (!is.null(setup$E)) colnames(setup$E))
  nm
}

.fb_build <- function(y, X, setup, spec) {
  Fu <- .fb_full(y, X, setup)
  M <- Fu$M[, .fb_colnames(spec, setup), drop = FALSE]
  if (spec$kstar > 0) M <- cbind(M, .fb_fourier(setup$n, spec$kstar)[Fu$rows, , drop = FALSE])
  list(M = M, Y = Fu$Y)
}

.fb_ic <- function(rss, N, K, ic) {
  if (ic == "aic") N * log(rss / N) + 2 * K else N * log(rss / N) + K * log(N)
}

# Two-step selection on the common sample: k* by minimum SSR of the model with
# every lag at maxlag, then the lags by AIC or BIC with k* fixed
.fb_select <- function(y, X, setup, fix) {
  Fu <- .fb_full(y, X, setup)
  Y <- Fu$Y
  N <- length(Y)
  nd <- length(setup$dec)
  nc <- length(setup$ctrl)
  ml <- setup$maxlag
  ssr_by_k <- NULL
  if (!setup$fourier) {
    kstar <- 0
  } else if (!is.null(fix$kstar)) {
    kstar <- fix$kstar
  } else {
    Mmax <- Fu$M[, .fb_colnames(list(p = ml, q = rep(ml, nd), r = rep(ml, nc)), setup), drop = FALSE]
    ssr <- vapply(setup$kvalues, function(k) {
      Xk <- cbind(Mmax, .fb_fourier(setup$n, k)[Fu$rows, , drop = FALSE])
      f <- stats::.lm.fit(Xk, Y)
      if (f$rank < ncol(Xk)) NA_real_ else sum(f$residuals^2)
    }, numeric(1))
    ssr_by_k <- data.frame(k = setup$kvalues, ssr = ssr)
    if (all(is.na(ssr))) stop("no Fourier frequency gives a full-rank model")
    kstar <- setup$kvalues[which.min(ssr)]
  }
  Fk <- if (kstar > 0) .fb_fourier(setup$n, kstar)[Fu$rows, , drop = FALSE] else NULL
  if (!is.null(fix$lags)) {
    sp <- c(fix$lags, kstar = kstar)
    Xs <- cbind(Fu$M[, .fb_colnames(sp, setup), drop = FALSE], Fk)
    f <- stats::.lm.fit(Xs, Y)
    return(list(p = fix$lags$p, q = fix$lags$q, r = fix$lags$r, kstar = kstar,
                ic = .fb_ic(sum(f$residuals^2), N, ncol(Xs), setup$ic),
                ssr_by_k = ssr_by_k, total_models = 1L))
  }
  grid <- expand.grid(c(list(p = seq_len(ml)), rep(list(0:ml), nd + nc)))
  best <- Inf
  bi <- 1
  for (g in seq_len(nrow(grid))) {
    gv <- as.integer(grid[g, ])
    sp <- list(p = gv[1], q = gv[1 + seq_len(nd)], r = gv[1 + nd + seq_len(nc)])
    Xs <- cbind(Fu$M[, .fb_colnames(sp, setup), drop = FALSE], Fk)
    f <- stats::.lm.fit(Xs, Y)
    if (f$rank < ncol(Xs)) next
    icv <- .fb_ic(sum(f$residuals^2), N, ncol(Xs), setup$ic)
    if (icv < best) {
      best <- icv
      bi <- g
    }
  }
  if (!is.finite(best)) stop("lag selection failed: every candidate is rank deficient")
  gv <- as.integer(grid[bi, ])
  list(p = gv[1], q = gv[1 + seq_len(nd)], r = gv[1 + nd + seq_len(nc)],
       kstar = kstar, ic = best, ssr_by_k = ssr_by_k, total_models = nrow(grid))
}

# Fov, t and Find (restricted residual sums of squares, OLS)
.fb_stats <- function(M, Y, setup) {
  f <- stats::.lm.fit(M, Y)
  if (f$rank < ncol(M)) stop("rank-deficient design")
  rss <- sum(f$residuals^2)
  df <- length(Y) - ncol(M)
  lev <- grep("^L_", colnames(M), value = TRUE)
  Ft <- function(drop) {
    fr <- stats::.lm.fit(M[, setdiff(colnames(M), drop), drop = FALSE], Y)
    ((sum(fr$residuals^2) - rss) / length(drop)) / (rss / df)
  }
  iy <- which(colnames(M) == "L_y")
  qx <- qr(M)
  XtXi <- chol2inv(qr.R(qx))[order(qx$pivot), order(qx$pivot)]
  bq <- qr.coef(qx, Y)
  c(Fov = Ft(lev), t = unname(bq[iy] / sqrt(rss / df * XtXi[iy, iy])),
    Find = Ft(setdiff(lev, "L_y")))
}

.fb_vcov <- function(X, e, type, XtXi) {
  n <- nrow(X)
  k <- ncol(X)
  if (type == "ols") return(sum(e^2) / (n - k) * XtXi)
  u <- X * e
  S <- crossprod(u)
  if (type == "HAC") {
    L <- floor(4 * (n / 100)^(2 / 9))
    for (l in seq_len(L)) {
      G <- crossprod(u[(l + 1):n, , drop = FALSE], u[1:(n - l), , drop = FALSE])
      S <- S + (1 - l / (L + 1)) * (G + t(G))
    }
  }
  n / (n - k) * XtXi %*% S %*% XtXi
}

.fb_multipliers <- function(b, V, setup, spec, dfr) {
  nm <- names(b)
  irho <- which(nm == "L_y")
  out <- list()
  lrrow <- function(lev) {
    i <- which(nm == lev)
    ld <- .lr_delta(b, V, irho, i)
    z <- ld$lr / ld$se
    c(lr = ld$lr, lr_se = ld$se, lr_t = z, lr_p = 2 * stats::pt(-abs(z), dfr))
  }
  for (i in seq_along(setup$dec)) {
    v <- setup$dec[i]
    sr <- function(s) {
      m <- paste0(v, s)
      sum(b[c(paste0("d_", m), if (spec$q[i] > 0) paste0("d_", m, "_L", seq_len(spec$q[i])))])
    }
    pp <- lrrow(paste0("L_", v, "_pos"))
    nn <- lrrow(paste0("L_", v, "_neg"))
    out[[v]] <- list(sr_pos = sr("_pos"), sr_neg = sr("_neg"),
                     lr_pos = unname(pp["lr"]), lr_neg = unname(nn["lr"]),
                     lr_pos_se = unname(pp["lr_se"]), lr_neg_se = unname(nn["lr_se"]),
                     lr_pos_t = unname(pp["lr_t"]), lr_neg_t = unname(nn["lr_t"]),
                     lr_pos_p = unname(pp["lr_p"]), lr_neg_p = unname(nn["lr_p"]))
  }
  for (j in seq_along(setup$ctrl)) {
    v <- setup$ctrl[j]
    s <- sum(b[c(paste0("d_", v), if (spec$r[j] > 0) paste0("d_", v, "_L", seq_len(spec$r[j])))])
    cc <- lrrow(paste0("L_", v))
    out[[v]] <- list(sr = s, lr = unname(cc["lr"]), lr_se = unname(cc["lr_se"]),
                     lr_t = unname(cc["lr_t"]), lr_p = unname(cc["lr_p"]))
  }
  out
}

.fb_asymmetry <- function(b, V, setup, spec, dfr) {
  nm <- names(b)
  out <- list()
  for (i in seq_along(setup$dec)) {
    v <- setup$dec[i]
    # long run: delta method on -b_pos/rho - (-b_neg/rho)
    ir <- which(nm == "L_y")
    ip <- which(nm == paste0("L_", v, "_pos"))
    ineg <- which(nm == paste0("L_", v, "_neg"))
    rho <- b[ir]
    d <- unname((b[ineg] - b[ip]) / rho)
    g <- numeric(length(b))
    g[ip] <- -1 / rho
    g[ineg] <- 1 / rho
    g[ir] <- (b[ip] - b[ineg]) / rho^2
    vd <- drop(t(g) %*% V %*% g)
    lr_w <- d^2 / vd
    # short run: sum over lags 0..q of the positive minus the negative terms
    R <- numeric(length(b))
    lg <- c("", if (spec$q[i] > 0) paste0("_L", seq_len(spec$q[i])))
    R[match(paste0("d_", v, "_pos", lg), nm)] <- 1
    R[match(paste0("d_", v, "_neg", lg), nm)] <- -1
    sr_w <- .wald_lin(b, V, matrix(R, 1))
    out[[v]] <- list(lr_diff = d, lr_chi2 = lr_w,
                     lr_p = stats::pchisq(lr_w, 1, lower.tail = FALSE),
                     sr_diff = sum(R * b), sr_f = sr_w,
                     sr_p = stats::pf(sr_w, 1, dfr, lower.tail = FALSE))
  }
  out
}

.fb_mult_path <- function(b, nm, v, s, p, q, H) {
  m <- paste0(v, s)
  pi <- b[match(c(paste0("d_", m), if (q > 0) paste0("d_", m, "_L", seq_len(q))), nm)]
  phi <- if (p > 0) b[match(paste0("dy_L", seq_len(p)), nm)] else numeric(0)
  .ecm_multiplier(b[match("L_y", nm)], b[match(paste0("L_", m), nm)], phi, pi, H)
}

.fb_dynamic_multipliers <- function(b, V, setup, spec, H, bands, level) {
  nm <- names(b)
  out <- list()
  draws <- NULL
  if (bands > 0) draws <- MASS::mvrnorm(bands, b, V)
  a <- (1 - level) / 2
  for (i in seq_along(setup$dec)) {
    v <- setup$dec[i]
    mp <- .fb_mult_path(b, nm, v, "_pos", spec$p, spec$q[i], H)
    mn <- .fb_mult_path(b, nm, v, "_neg", spec$p, spec$q[i], H)
    res <- list(horizon = 0:H, H_pos = mp, H_neg = mn, asymmetry = mp - mn)
    if (!is.null(draws)) {
      P <- t(apply(draws, 1, function(bb) .fb_mult_path(bb, nm, v, "_pos", spec$p, spec$q[i], H)))
      Nn <- t(apply(draws, 1, function(bb) .fb_mult_path(bb, nm, v, "_neg", spec$p, spec$q[i], H)))
      qf <- function(Mx) apply(Mx, 2, stats::quantile, probs = c(a, 1 - a), names = FALSE)
      res$pos_band <- t(qf(P))
      res$neg_band <- t(qf(Nn))
      res$asym_band <- t(qf(P - Nn))
      colnames(res$pos_band) <- colnames(res$neg_band) <- colnames(res$asym_band) <- c("lower", "upper")
    }
    out[[v]] <- res
  }
  out
}

# Diagnostics (package choices): Breusch-Godfrey LM with lagged residuals set
# to zero before the sample, ARCH LM, Breusch-Pagan on all regressors,
# Jarque-Bera and RESET with powers 2 and 3 of the fitted values
.fb_diagnostics <- function(Z, y, e) {
  n <- length(e)
  bg <- function(o) {
    El <- sapply(seq_len(o), function(j) c(rep(0, j), e[seq_len(n - j)]))
    f <- stats::.lm.fit(cbind(Z, El), e)
    r2 <- 1 - sum(f$residuals^2) / sum((e - mean(e))^2)
    s <- n * r2
    c(statistic = s, df = o, p_value = stats::pchisq(s, o, lower.tail = FALSE))
  }
  arch <- function(o) {
    e2 <- e^2
    Y <- e2[(o + 1):n]
    Xa <- cbind(1, sapply(seq_len(o), function(j) e2[(o + 1 - j):(n - j)]))
    f <- stats::.lm.fit(Xa, Y)
    r2 <- 1 - sum(f$residuals^2) / sum((Y - mean(Y))^2)
    s <- length(Y) * r2
    c(statistic = s, df = o, p_value = stats::pchisq(s, o, lower.tail = FALSE))
  }
  e2 <- e^2
  fb <- stats::.lm.fit(Z, e2)
  r2b <- 1 - sum(fb$residuals^2) / sum((e2 - mean(e2))^2)
  bp <- c(statistic = n * r2b, df = ncol(Z) - 1,
          p_value = stats::pchisq(n * r2b, ncol(Z) - 1, lower.tail = FALSE))
  s3 <- mean(e^3) / mean(e^2)^1.5
  k4 <- mean(e^4) / mean(e^2)^2
  jbs <- n / 6 * (s3^2 + (k4 - 3)^2 / 4)
  jb <- c(statistic = jbs, df = 2, p_value = stats::pchisq(jbs, 2, lower.tail = FALSE))
  yh <- y - e
  fr <- stats::.lm.fit(cbind(Z, yh^2, yh^3), y)
  rssu <- sum(fr$residuals^2)
  dfr <- n - ncol(Z) - 2
  Fr <- ((sum(e^2) - rssu) / 2) / (rssu / dfr)
  reset <- c(statistic = Fr, df1 = 2, df2 = dfr,
             p_value = stats::pf(Fr, 2, dfr, lower.tail = FALSE))
  list(bg1 = bg(1), bg4 = bg(4), arch1 = arch(1), arch4 = arch(4), bp = bp,
       jb = jb, reset = reset)
}


#' @rdname fbnardl
#' @param x,object An object of class "fbnardl"
#' @param ... Not used
#' @export
print.fbnardl <- function(x, ...) {
  cat("\n", if (x$type == "fbnardl") "Fourier bootstrap NARDL" else "Fourier NARDL",
      if (!x$fourier) " (no Fourier terms)", "\n", sep = "")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Decomposed:", paste(x$decomposed_vars, collapse = ", "), "\n")
  cat("Lags: p =", x$lags$p, "| q =", paste(x$lags$q, collapse = ","))
  if (length(x$control_vars)) cat(" | r =", paste(x$lags$r, collapse = ","))
  cat("\n")
  if (x$fourier) cat("k* =", format(x$kstar), "\n")
  cat("Observations:", x$nobs, "\n")
  bt <- x$bounds_test
  cat(sprintf("Fov = %.4f, t = %.4f, Find = %.4f\n", bt$Fov, bt$t, bt$Find))
  cat("Decision:", bt$decision, "\n")
  invisible(x)
}


#' @rdname fbnardl
#' @export
summary.fbnardl <- function(object, ...) {
  x <- object
  line <- function() cat(paste(rep("-", 72), collapse = ""), "\n")
  cat("\n")
  line()
  cat(" ", if (x$type == "fbnardl") "Fourier bootstrap NARDL" else "Fourier NARDL", "\n")
  line()
  cat("  Lags: p =", x$lags$p, "| q =", paste(x$lags$q, collapse = ","))
  if (length(x$control_vars)) cat(" | r =", paste(x$lags$r, collapse = ","))
  cat("\n")
  if (x$fourier) cat("  Fourier frequency k* =", format(x$kstar),
                     if (isTRUE(x$kstar_fixed)) "(fixed by the user)\n" else
                       sprintf("(selected on the %s grid)\n", x$kgrid))
  if (is.finite(x$ic_value)) cat("  ", toupper(x$ic), " = ", sprintf("%.4f", x$ic_value), "\n", sep = "")
  cat("  Observations:", x$nobs, " R-squared:", sprintf("%.4f", x$r.squared),
      " Adj. R-squared:", sprintf("%.4f", x$adj.r.squared), "\n")
  line()
  se <- x$se
  tt <- x$coefficients / se
  tab <- cbind(Estimate = x$coefficients, `Std. Error` = se, `t value` = tt,
               `Pr(>|t|)` = 2 * stats::pt(-abs(tt), x$df.residual))
  cat("  Coefficients (", x$vcov_type, " standard errors):\n", sep = "")
  print(round(tab, 4))
  line()
  cat("  Long-run coefficients (delta method):\n")
  for (v in names(x$multipliers)) {
    m <- x$multipliers[[v]]
    if (!is.null(m$lr_pos)) {
      cat(sprintf("  %-12s LR+ = %8.4f (%.4f)   LR- = %8.4f (%.4f)\n", v,
                  m$lr_pos, m$lr_pos_se, m$lr_neg, m$lr_neg_se))
    } else {
      cat(sprintf("  %-12s LR  = %8.4f (%.4f)\n", v, m$lr, m$lr_se))
    }
  }
  cat("  Symmetry tests:\n")
  for (v in names(x$asymmetry_tests)) {
    a <- x$asymmetry_tests[[v]]
    cat(sprintf("  %-12s long run: Wald = %.4f, p = %.4f; short run (sum of lags): F = %.4f, p = %.4f\n",
                v, a$lr_chi2, a$lr_p, a$sr_f, a$sr_p))
  }
  line()
  bt <- x$bounds_test
  cat(sprintf("  Cointegration: Fov = %.4f, t = %.4f, Find = %.4f (k = %d, n = %d, sr = %d)\n",
              bt$Fov, bt$t, bt$Find, bt$k, bt$n, bt$sr))
  if (identical(bt$method, "bootstrap")) {
    eng <- bt$bootstrap
    cat("  Bootstrap (", if (x$bootstrap_scheme == "bvz")
      paste0("separate nulls, Bertelli, Vacca and Zoia 2022; x* model: ",
             if (identical(x$xdgp, "vecm")) "marginal model with levels and lags (BVZ eqs. 19-20)"
             else "dx on intercept, Fourier terms and exog only (package choice, not BVZ)")
      else "Fov null for all statistics, McNown, Sam and Goh 2018 Steps 1-8, applied to the conditional ECM (package choice)",
        "), ", eng$settings$B,
        " replications, valid: ", paste(eng$n_valid, collapse = ", "), "\n", sep = "")
    if (isTRUE(eng$settings$reselect))
      cat("  k* and/or lags re-selected on each bootstrap sample (package choice;\n",
          "  BVZ step 5(c) re-estimates the unrestricted model but does not discuss re-selection)\n")
    cat(sprintf("  p-values: Fov %.4f, t %.4f, Find %.4f\n",
                eng$p_value["Fov"], eng$p_value["t"], eng$p_value["Find"]))
    cat(sprintf("  5%% critical values: Fov %.3f, t %.3f, Find %.3f\n",
                eng$cv["5%", "Fov"], eng$cv["5%", "t"], eng$cv["5%", "Find"]))
  } else {
    if (bt$bounds_valid) cat("  Kripfganz and Schneider (2020) bounds, 5%:\n")
    else cat("  Bounds for the model without Fourier terms; not valid with Fourier terms (5%):\n")
    cat(sprintf("  F: I(0) = %.3f, I(1) = %.3f   t: I(0) = %.3f, I(1) = %.3f\n",
                bt$bounds_F["I0", "5%"], bt$bounds_F["I1", "5%"],
                bt$bounds_t["I0", "5%"], bt$bounds_t["I1", "5%"]))
    if (bt$bounds_valid)
      cat(sprintf("  Approximate p-values (I(0), I(1)): F %.4f, %.4f; t %.4f, %.4f\n",
                  bt$pvalue_F[1], bt$pvalue_F[2], bt$pvalue_t[1], bt$pvalue_t[2]))
  }
  cat("  Decision:", bt$decision, "\n ", bt$message, "\n")
  line()
  dg <- x$diagnostics
  cat(sprintf("  Diagnostics (p-values): BG(1) %.4f, BG(4) %.4f, ARCH(1) %.4f, BP %.4f, JB %.4f, RESET %.4f\n",
              dg$bg1["p_value"], dg$bg4["p_value"], dg$arch1["p_value"],
              dg$bp["p_value"], dg$jb["p_value"], dg$reset["p_value"]))
  line()
  invisible(object)
}


#' Plot method for fbnardl objects
#'
#' @param x An object of class \code{"fbnardl"}
#' @param type \code{"multipliers"} (dynamic multipliers with bands) or
#'   \code{"kstar"} (sum of squared residuals by Fourier frequency).
#' @param ... Not used
#' @export
plot.fbnardl <- function(x, type = c("multipliers", "kstar"), ...) {
  type <- match.arg(type)
  if (type == "kstar") {
    if (is.null(x$ssr_by_k)) stop("no Fourier frequency search to plot")
    graphics::plot(x$ssr_by_k$k, x$ssr_by_k$ssr, type = "b", xlab = "k",
                   ylab = "SSR", main = "Fourier frequency selection")
    graphics::abline(v = x$kstar, lty = 2)
    return(invisible(x))
  }
  dm <- x$dynamic_multipliers
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, length(dm)))
  for (v in names(dm)) {
    d <- dm[[v]]
    rng <- range(c(d$H_pos, d$H_neg, d$asymmetry, d$pos_band, d$neg_band, d$asym_band), na.rm = TRUE)
    graphics::plot(d$horizon, d$H_pos, type = "l", ylim = rng, xlab = "Horizon",
                   ylab = "Cumulative multiplier", main = v)
    graphics::lines(d$horizon, d$H_neg, lty = 2)
    graphics::lines(d$horizon, d$asymmetry, col = "red")
    if (!is.null(d$asym_band)) {
      graphics::lines(d$horizon, d$asym_band[, 1], col = "red", lty = 3)
      graphics::lines(d$horizon, d$asym_band[, 2], col = "red", lty = 3)
    }
    graphics::abline(h = 0, col = "gray")
    graphics::legend("topleft", c("positive", "negative", "asymmetry"),
                     lty = c(1, 2, 1), col = c("black", "black", "red"), cex = 0.7)
  }
  invisible(x)
}


#' @rdname fbnardl
#' @importFrom stats coef vcov nobs
#' @export
coef.fbnardl <- function(object, ...) object$coefficients

#' @rdname fbnardl
#' @export
vcov.fbnardl <- function(object, ...) object$vcov

#' @rdname fbnardl
#' @export
nobs.fbnardl <- function(object, ...) object$nobs
