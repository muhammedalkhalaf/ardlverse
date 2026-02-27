\name{pnardl}
\alias{pnardl}
\alias{print.pnardl}
\alias{summary.pnardl}
\title{Panel Nonlinear ARDL (Panel NARDL)}
\description{
Panel data estimation with asymmetric/nonlinear ARDL specification using
PMG, MG, or DFE estimators.
}
\usage{
pnardl(formula, data, id, time, p = 1, q = 1,
       estimator = c("pmg", "mg", "dfe"),
       threshold = 0,
       effect = c("individual", "time", "twoways"),
       bootstrap = FALSE, nboot = 500)

\method{print}{pnardl}(x, ...)
\method{summary}{pnardl}(object, ...)
}
\arguments{
\item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
\item{data}{A data frame in long panel format}
\item{id}{Character. Name of the group/panel identifier variable}
\item{time}{Character. Name of the time variable}
\item{p}{Integer. Number of lags for dependent variable (default: 1)}
\item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
\item{estimator}{Character. One of "pmg", "mg", or "dfe" (default: "pmg")}
\item{threshold}{Numeric. Threshold for asymmetric decomposition (default: 0)}
\item{effect}{Character. Type of effects: "individual", "time", "twoways" (default: "individual")}
\item{bootstrap}{Logical. Use bootstrap for inference (default: FALSE)}
\item{nboot}{Number of bootstrap replications (default: 500)}
\item{x, object}{An object of class "pnardl"}
\item{...}{Additional arguments passed to methods}
}
\details{
The Panel NARDL model extends the Shin et al. (2014) NARDL framework to
panel data settings. It allows for asymmetric short-run and long-run
effects while controlling for cross-sectional heterogeneity.

Three estimators are available:
\itemize{
  \item \strong{PMG}: Pooled Mean Group - constrains long-run coefficients to be
    homogeneous while allowing short-run heterogeneity
  \item \strong{MG}: Mean Group - estimates separate models for each unit
    and averages coefficients
  \item \strong{DFE}: Dynamic Fixed Effects - assumes all coefficients
    are homogeneous except fixed effects
}
}
\value{
An object of class "pnardl" containing:
\item{coefficients}{Estimated coefficients}
\item{long_run_pos}{Long-run coefficients for positive changes}
\item{long_run_neg}{Long-run coefficients for negative changes}
\item{short_run}{Short-run coefficients}
\item{ec_coef}{Error correction coefficient}
\item{asymmetry_test}{Wald test for long-run asymmetry}
\item{unit_results}{Individual unit estimation results (for MG)}
\item{hausman}{Hausman test comparing PMG vs MG}
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework.

Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation
of dynamic heterogeneous panels. Journal of the American Statistical Association.
}
\examples{
\dontrun{
# Generate panel data
data <- generate_panel_data(N = 20, T = 50)

# Panel NARDL with PMG estimator
result <- pnardl(
  y ~ x1 + x2,
  data = data,
  id = "id",
  time = "time",
  estimator = "pmg"
)
summary(result)

# Panel NARDL with MG estimator
result_mg <- pnardl(
  y ~ x1,
  data = data,
  id = "id",
  time = "time",
  estimator = "mg"
)

# Panel NARDL with bootstrap standard errors
result_boot <- pnardl(
  y ~ x1,
  data = data,
  id = "id",
  time = "time",
  bootstrap = TRUE,
  nboot = 200
)
}
}
\seealso{
\code{\link{panel_ardl}}, \code{\link{qnardl}}, \code{\link{mtnardl}}
}
\keyword{models}
\keyword{panel}
