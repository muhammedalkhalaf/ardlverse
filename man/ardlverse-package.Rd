\name{ardlverse-package}
\alias{ardlverse-package}
\alias{ardlverse}
\docType{package}
\title{
Comprehensive ARDL Modeling Framework
}
\description{
A unified framework for Autoregressive Distributed Lag (ARDL) modeling 
and cointegration analysis. Implements Panel ARDL with Pooled Mean Group (PMG), 
Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators. Provides 
bootstrap-based bounds testing for small samples, Quantile Nonlinear ARDL 
(QNARDL) combining distributional and asymmetric effects, and Fourier ARDL 
for modeling smooth structural breaks.
}
\details{
The main functions are:
\itemize{
  \item \code{\link{panel_ardl}}: Panel ARDL with PMG, MG, DFE estimators
  \item \code{\link{boot_ardl}}: Bootstrap ARDL bounds test
  \item \code{\link{qnardl}}: Quantile Nonlinear ARDL
  \item \code{\link{fourier_ardl}}: Fourier ARDL for structural breaks
  \item \code{\link{ardl_diagnostics}}: Comprehensive model diagnostics
}

Supporting functions:
\itemize{
  \item \code{\link{hausman_test}}: PMG vs MG comparison
  \item \code{\link{asymmetry_test}}: Long-run asymmetry test
  \item \code{\link{dynamic_multipliers}}: Cumulative multipliers
  \item \code{\link{pss_critical_values}}: PSS (2001) critical values
}

Data generation for examples:
\itemize{
  \item \code{\link{generate_panel_data}}: Panel data
  \item \code{\link{generate_ts_data}}: Time series data
  \item \code{\link{generate_oil_data}}: Oil/gasoline price data
}
}
\author{
Muhammad Abdullah Alkhalaf <muhammedalkhalaf@gmail.com>

Maintainer: Muhammad Abdullah Alkhalaf <muhammedalkhalaf@gmail.com>
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation
of dynamic heterogeneous panels. Journal of the American Statistical Association.

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
to the analysis of level relationships. Journal of Applied Econometrics.

Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework.

Cho, J. S., Kim, T. H., & Shin, Y. (2015). Quantile cointegration in the
autoregressive distributed-lag modeling framework. Journal of Econometrics.

Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test.
Economic Modelling.
}
\keyword{package}
\examples{
\dontrun{
# Panel ARDL
data <- generate_panel_data()
model <- panel_ardl(gdp ~ inflation + investment, data = data,
                    id = "country", time = "year", estimator = "pmg")
summary(model)

# Bootstrap bounds test
ts_data <- generate_ts_data()
boot <- boot_ardl(gdp ~ investment, data = ts_data, nboot = 1000)
summary(boot)

# QNARDL
oil <- generate_oil_data()
qmodel <- qnardl(gasoline ~ oil_price, data = oil, tau = c(0.25, 0.5, 0.75))
asymmetry_test(qmodel, var = "oil_price")
}
}
