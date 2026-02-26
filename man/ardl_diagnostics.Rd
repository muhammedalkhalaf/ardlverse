\name{ardl_diagnostics}
\alias{ardl_diagnostics}
\alias{summary.ardl_diagnostics}
\alias{plot.ardl_diagnostics}
\title{ARDL Model Diagnostics}
\description{
Comprehensive diagnostic tests for ARDL models.
}
\usage{
ardl_diagnostics(model, lags = 4, arch_lags = 4)

\method{summary}{ardl_diagnostics}(object, ...)

\method{plot}{ardl_diagnostics}(x, which = 1:4, ...)
}
\arguments{
  \item{model}{An estimated model object (panel_ardl, boot_ardl, qnardl, or fourier_ardl)}
  \item{lags}{Integer. Number of lags for serial correlation tests (default: 4)}
  \item{arch_lags}{Integer. Number of lags for ARCH test (default: 4)}
  \item{object, x}{An object of class "ardl_diagnostics"}
  \item{which}{Integer vector. Which plots to display (1:4)}
  \item{...}{Additional arguments}
}
\details{
This function performs a battery of diagnostic tests:
\itemize{
  \item \strong{Serial Correlation}: Breusch-Godfrey LM test and Ljung-Box test

  \item \strong{Heteroskedasticity}: Breusch-Pagan and ARCH tests
  \item \strong{Normality}: Jarque-Bera and Shapiro-Wilk tests
  \item \strong{Functional Form}: RESET test
  \item \strong{Stability}: CUSUM and CUSUMSQ tests
}
}
\value{
An object of class "ardl_diagnostics" containing:
\item{serial_corr}{Breusch-Godfrey test results}
\item{ljung_box}{Ljung-Box test results}
\item{hetero_bp}{Breusch-Pagan test results}
\item{arch}{ARCH-LM test results}
\item{normality}{Jarque-Bera test results}
\item{shapiro}{Shapiro-Wilk test results}
\item{reset}{RESET test results}
\item{cusum}{CUSUM test results}
\item{cusumsq}{CUSUM of squares test results}
\item{residuals}{Model residuals}
}
\examples{
\dontrun{
# Estimate a model
data <- generate_ts_data(n = 100)
model <- fourier_ardl(gdp ~ investment, data = data)

# Run diagnostics
diag <- ardl_diagnostics(model)
summary(diag)
plot(diag)
}
}
\seealso{
\code{\link{fourier_ardl}}, \code{\link{boot_ardl}}
}
\keyword{regression}
\keyword{ts}
