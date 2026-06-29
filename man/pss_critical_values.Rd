\name{pss_critical_values}
\alias{pss_critical_values}
\title{PSS Asymptotic Critical Values}
\description{
Get Pesaran, Shin & Smith (2001) asymptotic critical values for bounds test.
}
\usage{
pss_critical_values(k, case = 3, level = "5\%")
}
\arguments{
  \item{k}{Number of regressors (excluding the lagged dependent variable)}
  \item{case}{Deterministic specification (1-5, default: 3)}
  \item{level}{Significance level: "10\%", "5\%", or "1\%"}
}
\details{
Returns the I(0) and I(1) critical value bounds from the original
Pesaran, Shin & Smith (2001) tables.

Case specifications:
\itemize{
  \item 1: No intercept, no trend
  \item 2: Restricted intercept, no trend
  \item 3: Unrestricted intercept, no trend (default)
  \item 4: Unrestricted intercept, restricted trend
  \item 5: Unrestricted intercept, unrestricted trend
}
}
\value{
A list with:
\item{F_bounds}{List with I0 and I1 bounds for F-test}
\item{t_bounds}{List with I0 and I1 bounds for t-test}
\item{k}{Number of regressors}
\item{case}{Case specification}
\item{level}{Significance level}
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
to the analysis of level relationships. Journal of Applied Econometrics,
16(3), 289-326.
}
\examples{
# Get 5\% critical values for 3 regressors, Case III
cv <- pss_critical_values(k = 3, case = 3, level = "5\%")
cv$F_bounds  # I(0) and I(1) bounds for F-test
cv$t_bounds  # I(0) and I(1) bounds for t-test
}
\seealso{
\code{\link{boot_ardl}}, \code{\link{fourier_bounds_test}}
}
\keyword{htest}
