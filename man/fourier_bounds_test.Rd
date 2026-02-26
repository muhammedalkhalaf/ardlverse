\name{fourier_bounds_test}
\alias{fourier_bounds_test}
\title{Fourier ARDL Bounds Test}
\description{
Perform bounds test with Fourier-adjusted critical values.
}
\usage{
fourier_bounds_test(object)
}
\arguments{
  \item{object}{A fourier_ardl object}
}
\details{
Performs the bounds test for cointegration using critical values
adjusted for the presence of Fourier terms. The critical values
account for the additional parameters estimated.
}
\value{
A list containing:
\item{F_stat}{F-statistic}
\item{t_stat}{t-statistic}
\item{cv_F}{Critical values for F-test by significance level}
}
\examples{
\dontrun{
data <- generate_ts_data(n = 100)

model <- fourier_ardl(
  gdp ~ investment,
  data = data,
  k = 1
)

fourier_bounds_test(model)
}
}
\seealso{
\code{\link{fourier_ardl}}, \code{\link{pss_critical_values}}
}
\keyword{htest}
