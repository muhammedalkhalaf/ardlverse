\name{hausman_test}
\alias{hausman_test}
\title{Hausman Test for PMG vs MG}
\description{
Perform Hausman test comparing PMG and MG estimators.
}
\usage{
hausman_test(pmg_model, mg_model = NULL, data = NULL)
}
\arguments{
  \item{pmg_model}{A panel_ardl object estimated with PMG}
  \item{mg_model}{A panel_ardl object estimated with MG (optional)}
  \item{data}{Data frame (required if mg_model not provided)}
}
\details{
The Hausman test examines whether the long-run homogeneity assumption
of the PMG estimator is valid. Under the null hypothesis of homogeneity,
both PMG and MG are consistent, but PMG is efficient. Under the alternative,
only MG is consistent.
}
\value{
A list containing:
\item{statistic}{Hausman chi-squared statistic}
\item{df}{Degrees of freedom}
\item{p.value}{P-value}
\item{theta_pmg}{PMG long-run coefficients}
\item{theta_mg}{MG long-run coefficients}
\item{difference}{Coefficient differences}
}
\examples{
\dontrun{
data <- generate_panel_data(n_groups = 10, n_time = 50)

pmg <- panel_ardl(gdp ~ inflation + investment, data = data,
                  id = "country", time = "year", estimator = "pmg")

mg <- panel_ardl(gdp ~ inflation + investment, data = data,
                 id = "country", time = "year", estimator = "mg")

hausman_test(pmg, mg)
}
}
\seealso{
\code{\link{panel_ardl}}
}
\keyword{htest}
