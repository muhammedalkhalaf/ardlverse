\name{asymmetry_test}
\alias{asymmetry_test}
\title{Test Asymmetry in QNARDL Model}
\description{
Perform Wald test for long-run asymmetry in QNARDL models.
}
\usage{
asymmetry_test(object, var)
}
\arguments{
  \item{object}{A qnardl object}
  \item{var}{Variable name to test asymmetry for}
}
\details{
Tests the null hypothesis that positive and negative long-run effects
are equal (H0: theta+ = theta-) against the alternative of asymmetric
effects (H1: theta+ != theta-).

The test is performed separately for each quantile.
}
\value{
A data frame with test results by quantile:
\item{tau}{Quantile}
\item{theta_pos}{Positive long-run coefficient}
\item{theta_neg}{Negative long-run coefficient}
\item{diff}{Difference (theta+ - theta-)}
\item{wald_stat}{Wald test statistic}
\item{p_value}{P-value}
}
\examples{
\dontrun{
data <- generate_oil_data(n = 200)

model <- qnardl(
  gasoline ~ oil_price,
  data = data,
  tau = c(0.25, 0.5, 0.75)
)

asymmetry_test(model, var = "oil_price")
}
}
\seealso{
\code{\link{qnardl}}, \code{\link{dynamic_multipliers}}
}
\keyword{htest}
