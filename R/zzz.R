#' @title Package Startup
#' @description Functions run on package load
#' @keywords internal

.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.ardlverse <- list(
    ardlverse.digits = 4,
    ardlverse.pvalue_threshold = 0.05,
    ardlverse.bootstrap_reps = 2000
  )
  
  toset <- !(names(op.ardlverse) %in% names(op))
  if (any(toset)) options(op.ardlverse[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "====================================================================\n",
    "                        ardlverse v1.0.0\n",
    "====================================================================\n",
    " Comprehensive ARDL Modeling Framework\n",
    "--------------------------------------------------------------------\n",
    " Components:\n",
    "   - panel_ardl()     : Panel ARDL (PMG, MG, DFE)\n",
    "   - boot_ardl()      : Bootstrap Bounds Test\n",
    "   - qnardl()         : Quantile Nonlinear ARDL\n",
    "   - fourier_ardl()   : Fourier ARDL\n",
    "   - ardl_diagnostics(): Model Diagnostics\n",
    "--------------------------------------------------------------------\n",
    " Author: Muhammad Alkhalaf <contact@rufyqelngeh.com>\n",
    " Documentation: https://github.com/muhammedalkhalaf/ardlverse\n",
    "===================================================================="
  )
}

# Global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  "dy", "y_lag1", "tau", "theta", "type", "multiplier", "horizon",
  "residuals", "fitted", "t", "cusum", "upper", "lower", "fourier",
  "F_stat", "t_stat", "density", "after_stat"
))
