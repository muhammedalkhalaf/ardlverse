\name{mtnardl}
\alias{mtnardl}
\alias{print.mtnardl}
\alias{summary.mtnardl}
\alias{plot.mtnardl}
\title{Multiple-Threshold Nonlinear ARDL (MT-NARDL)}
\description{
Extends NARDL to allow multiple threshold decomposition for capturing 
complex asymmetric relationships.
}
\usage{
mtnardl(formula, data, thresholds = c(0), p = 1, q = 1, case = 3,
        auto_select = FALSE, n_thresholds = 2,
        bootstrap = FALSE, nboot = 2000, seed = NULL)

\method{print}{mtnardl}(x, ...)
\method{summary}{mtnardl}(object, ...)
\method{plot}{mtnardl}(x, type = c("multipliers", "asymmetry"), ...)
}
\arguments{
\item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
\item{data}{A data frame containing the time series data}
\item{thresholds}{Numeric vector of threshold values (default: c(0))}
\item{p}{Integer. Number of lags for dependent variable (default: 1)}
\item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
\item{case}{Integer from 1-5 specifying deterministic components (default: 3)}
\item{auto_select}{Logical. Automatically select optimal thresholds (default: FALSE)}
\item{n_thresholds}{Integer. Number of thresholds to select if auto_select = TRUE}
\item{bootstrap}{Logical. Use bootstrap inference (default: FALSE)}
\item{nboot}{Number of bootstrap replications (default: 2000)}
\item{seed}{Random seed for reproducibility}
\item{x, object}{An object of class "mtnardl"}
\item{type}{Character. Plot type: "multipliers" or "asymmetry"}
\item{...}{Additional arguments passed to methods}
}
\details{
The Multiple-Threshold NARDL (MT-NARDL) model extends the standard NARDL
framework by allowing decomposition of variables into multiple regimes
based on user-specified thresholds. This captures more nuanced asymmetric
effects beyond simple positive/negative decomposition.

For example, with thresholds c(-0.02, 0, 0.02), a variable is decomposed into:
\itemize{
  \item Large decreases (< -2\%)
  \item Small decreases (-2\% to 0)
  \item Small increases (0 to 2\%)
  \item Large increases (> 2\%)
}
}
\value{
An object of class "mtnardl" containing:
\item{model}{The estimated MT-NARDL model}
\item{bounds_test}{Bounds test results}
\item{long_run}{Long-run coefficients for each regime}
\item{short_run}{Short-run coefficients for each regime}
\item{thresholds}{Threshold values used}
\item{asymmetry_tests}{Wald tests for asymmetry between regimes}
\item{multipliers}{Dynamic multipliers for each regime}
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework.
In Festschrift in Honor of Peter Schmidt (pp. 281-314). Springer.
}
\examples{
\dontrun{
# Generate example data
data <- generate_oil_data(n = 300)

# Standard NARDL (single threshold at 0)
result1 <- mtnardl(consumption ~ oil_price, data = data)

# Multiple thresholds for different shock sizes
result2 <- mtnardl(
  consumption ~ oil_price,
  data = data,
  thresholds = c(-0.05, 0, 0.05)
)
summary(result2)
plot(result2, type = "multipliers")

# Auto-select optimal thresholds
result3 <- mtnardl(
  consumption ~ oil_price,
  data = data,
  auto_select = TRUE,
  n_thresholds = 2
)
}
}
\seealso{
\code{\link{qnardl}}, \code{\link{aardl}}
}
\keyword{models}
\keyword{asymmetry}
