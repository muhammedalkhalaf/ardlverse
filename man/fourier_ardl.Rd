\name{fourier_ardl}
\alias{fourier_ardl}
\alias{summary.fourier_ardl}
\alias{print.fourier_ardl}
\alias{plot.fourier_ardl}
\title{Fourier ARDL Model}
\description{
Estimate ARDL models with Fourier terms for smooth structural breaks.
}
\usage{
fourier_ardl(formula, data, p = 1, q = 1, k = 1, case = 3,
             selection = c("fixed", "aic", "bic"))

\method{summary}{fourier_ardl}(object, ...)

\method{print}{fourier_ardl}(x, ...)

\method{plot}{fourier_ardl}(x, which = "both", ...)
}
\arguments{
  \item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
  \item{data}{A data frame containing the time series data}
  \item{p}{Integer. Number of lags for dependent variable (default: 1)}
  \item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
  \item{k}{Integer. Number of Fourier frequencies (default: 1, max: 3)}
  \item{case}{Integer from 1-5 specifying deterministic components (default: 3)}
  \item{selection}{Character. Method for selecting optimal k: "aic", "bic", or "fixed"}
  \item{object, x}{An object of class "fourier_ardl"}
  \item{which}{Character. "fourier", "fit", or "both" for plotting}
  \item{...}{Additional arguments}
}
\details{
This function implements the Fourier ARDL approach that captures smooth
structural changes using trigonometric functions. The Fourier terms approximate
unknown structural breaks without requiring prior specification of break dates.

The Fourier terms are defined as:
\deqn{f_t = \sum_{k=1}^{K} [a_k \sin(2\pi k t/T) + b_k \cos(2\pi k t/T)]}

where K is the number of frequencies and T is the sample size.
}
\value{
An object of class "fourier_ardl" containing:
\item{coefficients}{Estimated coefficients}
\item{long_run}{Long-run coefficients}
\item{fourier_coefs}{Fourier term coefficients}
\item{ec_coef}{Error correction coefficient}
\item{bounds_test}{F and t statistics for bounds test}
\item{k}{Number of Fourier frequencies used}
\item{aic}{AIC values for k = 1, 2, 3}
\item{bic}{BIC values for k = 1, 2, 3}
\item{structural_breaks}{Detected structural break periods}
}
\references{
Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test
to approximate smooth breaks with new evidence from crude oil market.
Economic Modelling, 67, 114-124.

Enders, W., & Lee, J. (2012). A unit root test using a Fourier series to
approximate smooth breaks. Oxford Bulletin of Economics and Statistics,
74(4), 574-599.
}
\examples{
\dontrun{
# Generate example data
data <- generate_ts_data(n = 100)

# Estimate Fourier ARDL with automatic selection
model <- fourier_ardl(
  gdp ~ investment + trade,
  data = data,
  selection = "aic"
)
summary(model)
plot(model)
}
}
\seealso{
\code{\link{fourier_bounds_test}}, \code{\link{generate_ts_data}}
}
\keyword{regression}
\keyword{ts}
