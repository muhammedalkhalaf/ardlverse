\name{generate_ts_data}
\alias{generate_ts_data}
\title{Generate Example Time Series Data}
\description{
Generate simulated time series data for examples and testing.
}
\usage{
generate_ts_data(n = 100, seed = 456)
}
\arguments{
  \item{n}{Number of observations (default: 100)}
  \item{seed}{Random seed for reproducibility (default: 456)}
}
\details{
Generates a simulated time series dataset with cointegrated variables
and a smooth structural break (modeled via Fourier component). Suitable
for demonstrating Fourier ARDL and bootstrap bounds testing methods.
}
\value{
A data frame with columns:
\item{quarter}{Quarter number (1 to n)}
\item{gdp}{Log real GDP (I(1) with structural break)}
\item{inflation}{Inflation rate (I(0))}
\item{investment}{Investment growth (I(1))}
\item{trade}{Trade balance (I(1))}
}
\examples{
# Generate 100 quarters of data
ts_data <- generate_ts_data(n = 100)
head(ts_data)

# Plot GDP series
plot(ts_data$quarter, ts_data$gdp, type = "l",
     main = "Simulated GDP with Structural Break")
}
\seealso{
\code{\link{fourier_ardl}}, \code{\link{boot_ardl}}, 
\code{\link{generate_panel_data}}
}
\keyword{datagen}
