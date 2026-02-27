\name{rardl}
\alias{rardl}
\alias{print.rardl}
\alias{summary.rardl}
\alias{plot.rardl}
\title{Rolling and Recursive ARDL (R-ARDL)}
\description{
Time-varying ARDL bounds test using rolling or recursive windows to detect
structural changes in cointegration relationships.
}
\usage{
rardl(formula, data, method = c("rolling", "recursive"),
      window = 50, min_obs = 40, p = 1, q = 1, case = 3,
      parallel = FALSE, ncores = 2)

\method{print}{rardl}(x, ...)
\method{summary}{rardl}(object, ...)
\method{plot}{rardl}(x, type = c("F", "ec", "lr", "all"), ...)
}
\arguments{
\item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
\item{data}{A data frame containing the time series data}
\item{method}{Character. Either "rolling" or "recursive" (default: "rolling")}
\item{window}{Integer. Window size for rolling method (default: 50)}
\item{min_obs}{Integer. Minimum observations for recursive method (default: 40)}
\item{p}{Integer. Number of lags for dependent variable (default: 1)}
\item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
\item{case}{Integer from 1-5 specifying deterministic components (default: 3)}
\item{parallel}{Logical. Use parallel processing (default: FALSE)}
\item{ncores}{Integer. Number of cores for parallel processing}
\item{x, object}{An object of class "rardl"}
\item{type}{Character. Plot type: "F", "ec", "lr", or "all"}
\item{...}{Additional arguments passed to methods}
}
\details{
The Rolling ARDL (R-ARDL) approach implements time-varying bounds testing
to detect structural changes in cointegration relationships. Two methods
are available:

\itemize{
  \item \strong{Rolling}: Fixed-size window that moves through the sample
  \item \strong{Recursive}: Expanding window starting from initial observations
}

This is useful for detecting:
\itemize{
  \item Time-varying cointegration relationships
  \item Structural breaks in long-run parameters
  \item Changes in error correction speed
}
}
\value{
An object of class "rardl" containing:
\item{F_stats}{Time series of F-statistics}
\item{t_stats}{Time series of t-statistics (EC coefficient)}
\item{cointegration}{Time series of cointegration decisions}
\item{ec_coefs}{Time-varying error correction coefficients}
\item{lr_coefs}{Time-varying long-run coefficients}
\item{dates}{Index/dates for each window}
\item{breaks}{Detected structural break points}
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
to the analysis of level relationships. Journal of Applied Econometrics, 16(3), 289-326.
}
\examples{
\dontrun{
# Generate example data
data <- generate_ts_data(n = 300)

# Rolling ARDL with 60-observation window
roll_result <- rardl(y ~ x1 + x2, data = data, method = "rolling", window = 60)
plot(roll_result, type = "all")
summary(roll_result)

# Recursive ARDL
rec_result <- rardl(y ~ x1 + x2, data = data, method = "recursive", min_obs = 50)
plot(rec_result, type = "F")
}
}
\seealso{
\code{\link{boot_ardl}}, \code{\link{aardl}}
}
\keyword{models}
\keyword{structural breaks}
