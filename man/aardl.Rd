\name{aardl}
\alias{aardl}
\alias{print.aardl}
\alias{summary.aardl}
\title{Augmented ARDL Bounds Test (AARDL)}
\description{
Implements the Augmented ARDL bounds testing approach with deferred t and F tests 
for cointegration analysis.
}
\usage{
aardl(formula, data, p = 1, q = 1, case = 3,
      type = c("linear", "nardl", "fourier", "fnardl",
               "bootstrap", "bnardl", "fbootstrap", "fbnardl"),
      nboot = 2000, fourier_k = 1, threshold = 0, seed = NULL)

\method{print}{aardl}(x, ...)
\method{summary}{aardl}(object, ...)
}
\arguments{
\item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
\item{data}{A data frame containing the time series data}
\item{p}{Integer. Number of lags for dependent variable (default: 1)}
\item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
\item{case}{Integer from 1-5 specifying deterministic components (default: 3)}
\item{type}{Character. Model type: "linear", "nardl", "fourier", "fnardl",
  "bootstrap", "bnardl", "fbootstrap", "fbnardl"}
\item{nboot}{Number of bootstrap replications (default: 2000)}
\item{fourier_k}{Integer. Number of Fourier frequencies (default: 1, max: 3)}
\item{threshold}{Numeric. Threshold value for NARDL decomposition (default: 0)}
\item{seed}{Random seed for reproducibility}
\item{x, object}{An object of class "aardl"}
\item{...}{Additional arguments passed to methods}
}
\details{
The Augmented ARDL (AARDL) approach extends the standard ARDL bounds test
by implementing additional diagnostic tests proposed by Sam, McNown & Goh (2019).
This addresses potential weaknesses in the PSS bounds test by adding:
\itemize{
  \item Deferred t-test (t_dep): Tests significance of lagged dependent variable
  \item Deferred F-test (F_ind): Tests joint significance of independent variables
  \item Overall F-test with all deferred conditions
}

The function supports 8 sub-models including standard ARDL, NARDL, Fourier ARDL, 
and their bootstrap variants.
}
\value{
An object of class "aardl" containing:
\item{F_pss}{PSS F-statistic for bounds test}
\item{t_dep}{Deferred t-statistic for lagged dependent variable}
\item{F_ind}{Deferred F-statistic for independent variables}
\item{conclusion}{Cointegration decision based on all tests}
\item{model}{The estimated ARDL model}
\item{long_run}{Long-run coefficients}
\item{short_run}{Short-run coefficients}
\item{diagnostics}{Model diagnostic tests}
}
\references{
Sam, C. Y., McNown, R., & Goh, S. K. (2019). An augmented autoregressive 
distributed lag bounds test for cointegration. Economic Modelling, 80, 130-141.

McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
distributed lag test for cointegration. Applied Economics, 50(13), 1509-1521.
}
\examples{
\dontrun{
# Generate example data
data <- generate_ts_data(n = 200)

# Standard Augmented ARDL
result <- aardl(y ~ x1 + x2, data = data, p = 2, q = 2, case = 3)
summary(result)

# Augmented NARDL (nonlinear)
result_nardl <- aardl(y ~ x1 + x2, data = data, type = "nardl")
summary(result_nardl)

# Fourier Augmented ARDL with bootstrap
result_fb <- aardl(y ~ x1, data = data, type = "fbootstrap", nboot = 1000)
}
}
\seealso{
\code{\link{boot_ardl}}, \code{\link{qnardl}}, \code{\link{fourier_ardl}}
}
\keyword{models}
\keyword{cointegration}
