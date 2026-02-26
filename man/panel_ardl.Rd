\name{panel_ardl}
\alias{panel_ardl}
\alias{summary.panel_ardl}
\alias{print.panel_ardl}
\title{Panel ARDL Estimation}
\description{
Estimate Panel ARDL models with Pooled Mean Group (PMG), Mean Group (MG), 
and Dynamic Fixed Effects (DFE) estimators.
}
\usage{
panel_ardl(formula, data, id, time, p = 1, q = 1,
           estimator = c("pmg", "mg", "dfe"),
           ec = TRUE, trend = FALSE,
           maxiter = 100, tol = 1e-5)

\method{summary}{panel_ardl}(object, ...)

\method{print}{panel_ardl}(x, ...)
}
\arguments{
  \item{formula}{A formula specifying the model: y ~ x1 + x2 + ...}
  \item{data}{A data frame containing panel data}
  \item{id}{Character string specifying the group/panel identifier variable}
  \item{time}{Character string specifying the time variable}
  \item{p}{Integer. Number of lags for dependent variable (default: 1)}
  \item{q}{Integer or vector. Number of lags for independent variables (default: 1)}
  \item{estimator}{Character. One of "pmg", "mg", or "dfe" (default: "pmg")}
  \item{ec}{Logical. If TRUE, display error correction form (default: TRUE)}
  \item{trend}{Logical. Include time trend (default: FALSE)}
  \item{maxiter}{Maximum iterations for PMG optimization (default: 100)}
  \item{tol}{Convergence tolerance (default: 1e-5)}
  \item{object, x}{An object of class "panel_ardl"}
  \item{...}{Additional arguments (ignored)}
}
\details{
This function implements the panel ARDL framework of Pesaran, Shin & Smith (1999)
for estimating long-run relationships in dynamic heterogeneous panels.

Three estimators are available:
\itemize{
  \item \strong{PMG}: Constrains long-run coefficients to be equal across groups 
        while allowing short-run coefficients to vary.
  \item \strong{MG}: Estimates separate regressions for each group and averages 
        coefficients.
  \item \strong{DFE}: Traditional dynamic fixed effects with all coefficients pooled.
}
}
\value{
An object of class "panel_ardl" containing:
\item{coefficients}{Estimated coefficients}
\item{long_run}{Long-run coefficients (theta)}
\item{short_run}{Short-run coefficients}
\item{ec_coef}{Error correction coefficient (phi)}
\item{se}{Standard errors}
\item{residuals}{Model residuals}
\item{fitted}{Fitted values}
\item{nobs}{Number of observations}
\item{ngroups}{Number of groups}
\item{aic}{Akaike Information Criterion}
\item{bic}{Bayesian Information Criterion}
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation
of dynamic heterogeneous panels. Journal of the American Statistical Association,
94(446), 621-634.
}
\examples{
\dontrun{
# Generate example data
data <- generate_panel_data(n_groups = 10, n_time = 50)

# Estimate PMG model
pmg_model <- panel_ardl(
  gdp ~ inflation + investment,
  data = data,
  id = "country",
  time = "year",
  estimator = "pmg"
)
summary(pmg_model)
}
}
\seealso{
\code{\link{hausman_test}}, \code{\link{generate_panel_data}}
}
\keyword{regression}
\keyword{panel}
