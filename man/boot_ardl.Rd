\name{boot_ardl}
\alias{boot_ardl}
\alias{summary.boot_ardl}
\alias{print.boot_ardl}
\alias{plot.boot_ardl}
\title{Bootstrap ARDL Bounds Test}
\description{
Perform bounds test for cointegration with bootstrap critical values.
}
\usage{
boot_ardl(formula, data, p = 1, q = 1, case = 3,
          nboot = 2000, seed = NULL, 
          parallel = FALSE, ncores = 2)

\method{summary}{boot_ardl}(object, ...)

\method{print}{boot_ardl}(x, ...)

\method{plot}{boot_ardl}(x, which = "both", ...)
}
\arguments{
  \item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
  \item{data}{A data frame containing the time series data}
  \item{p}{Integer. Number of lags for dependent variable (default: 1)}
  \item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
  \item{case}{Integer from 1-5 specifying deterministic components (default: 3)}
  \item{nboot}{Number of bootstrap replications (default: 2000)}
  \item{seed}{Random seed for reproducibility (default: NULL)}
  \item{parallel}{Logical. Use parallel processing (default: FALSE)}
  \item{ncores}{Number of cores for parallel processing (default: 2)}
  \item{object, x}{An object of class "boot_ardl"}
  \item{which}{Character. "F", "t", or "both" for plotting}
  \item{...}{Additional arguments}
}
\details{
This function implements bootstrap-based inference for the ARDL bounds test,
which is particularly useful for small samples where asymptotic critical values
may be unreliable.

Case specifications:
\itemize{
  \item 1: No intercept, no trend
  \item 2: Restricted intercept, no trend
  \item 3: Unrestricted intercept, no trend (default)
  \item 4: Unrestricted intercept, restricted trend
  \item 5: Unrestricted intercept, unrestricted trend
}
}
\value{
An object of class "boot_ardl" containing:
\item{F_stat}{F-statistic for bounds test}
\item{t_stat}{t-statistic for EC coefficient}
\item{boot_F}{Bootstrap distribution of F-statistics}
\item{boot_t}{Bootstrap distribution of t-statistics}
\item{cv_F}{Critical values for F-test (90\%, 95\%, 99\%)}
\item{cv_t}{Critical values for t-test}
\item{p_value_F}{Bootstrap p-value for F-test}
\item{p_value_t}{Bootstrap p-value for t-test}
\item{conclusion}{Test conclusion}
}
\references{
McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
distributed lag test for cointegration. Applied Economics, 50(13), 1509-1521.

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
to the analysis of level relationships. Journal of Applied Econometrics, 16(3), 289-326.
}
\examples{
\dontrun{
# Generate example data
data <- generate_ts_data(n = 100)

# Bootstrap bounds test
boot_test <- boot_ardl(
  gdp ~ inflation + investment,
  data = data,
  p = 2, q = 2,
  nboot = 1000
)
summary(boot_test)
plot(boot_test)
}
}
\seealso{
\code{\link{pss_critical_values}}, \code{\link{generate_ts_data}}
}
\keyword{regression}
\keyword{ts}
