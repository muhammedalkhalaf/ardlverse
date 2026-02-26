\name{qnardl}
\alias{qnardl}
\alias{summary.qnardl}
\alias{print.qnardl}
\alias{plot.qnardl}
\title{Quantile Nonlinear ARDL (QNARDL)}
\description{
Estimate Quantile Nonlinear ARDL models combining distributional
and asymmetric effects.
}
\usage{
qnardl(formula, data, tau = c(0.25, 0.5, 0.75),
       p = 1, q = 1, decompose = NULL, trend = FALSE)

\method{summary}{qnardl}(object, ...)

\method{print}{qnardl}(x, ...)

\method{plot}{qnardl}(x, var = NULL, type = "long_run", ...)
}
\arguments{
  \item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
  \item{data}{A data frame containing the time series data}
  \item{tau}{Numeric vector of quantiles to estimate (default: c(0.25, 0.5, 0.75))}
  \item{p}{Integer. Number of lags for dependent variable (default: 1)}
  \item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
  \item{decompose}{Character vector. Variables to decompose into +/- components.
    Default is all x variables.}
  \item{trend}{Logical. Include time trend (default: FALSE)}
  \item{object, x}{An object of class "qnardl"}
  \item{var}{Variable to plot (default: first decomposed variable)}
  \item{type}{"long_run" or "asymmetry"}
  \item{...}{Additional arguments}
}
\details{
This function implements the QNARDL model which extends the NARDL framework
of Shin, Yu & Greenwood-Nimmo (2014) to a quantile regression setting.

The model allows for:
\itemize{
  \item \strong{Asymmetric effects}: Positive and negative changes in X
        can have different impacts on Y
  \item \strong{Distributional heterogeneity}: Effects can vary across
        quantiles of the conditional distribution
}

The partial sum decomposition separates each regressor into positive and
negative components.
}
\value{
An object of class "qnardl" containing:
\item{coefficients}{Estimated coefficients for each quantile}
\item{long_run_pos}{Positive long-run coefficients by quantile}
\item{long_run_neg}{Negative long-run coefficients by quantile}
\item{short_run}{Short-run coefficients by quantile}
\item{ec_coef}{Error correction coefficients by quantile}
\item{asymmetry_test}{Wald test for long-run asymmetry}
\item{tau}{Quantiles estimated}
}
\references{
Cho, J. S., Kim, T. H., & Shin, Y. (2015). Quantile cointegration in the
autoregressive distributed-lag modeling framework. Journal of Econometrics,
188(1), 281-300.

Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework.
In Festschrift in Honor of Peter Schmidt (pp. 281-314). Springer.
}
\examples{
\dontrun{
# Generate oil price data
data <- generate_oil_data(n = 200)

# Estimate QNARDL
model <- qnardl(
  gasoline ~ oil_price + exchange_rate,
  data = data,
  tau = c(0.1, 0.25, 0.5, 0.75, 0.9)
)
summary(model)
plot(model, var = "oil_price")
}
}
\seealso{
\code{\link{asymmetry_test}}, \code{\link{dynamic_multipliers}}, 
\code{\link{generate_oil_data}}
}
\keyword{regression}
\keyword{nonlinear}
