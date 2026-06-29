\name{generate_panel_data}
\alias{generate_panel_data}
\title{Generate Example Panel Data}
\description{
Generate simulated panel data for examples and testing.
}
\usage{
generate_panel_data(n_groups = 10, n_time = 50, seed = 123)
}
\arguments{
  \item{n_groups}{Number of groups/countries (default: 10)}
  \item{n_time}{Number of time periods (default: 50)}
  \item{seed}{Random seed for reproducibility (default: 123)}
}
\details{
Generates a simulated panel dataset with cointegrated variables suitable
for demonstrating panel ARDL methods. The data generating process includes
country-specific fixed effects and error correction dynamics.
}
\value{
A data frame with columns:
\item{country}{Country identifier (1 to n_groups)}
\item{year}{Year (starting from 1970)}
\item{gdp}{Log GDP per capita (I(1))}
\item{inflation}{Inflation rate (I(0))}
\item{investment}{Investment as \% of GDP (I(1))}
\item{trade}{Trade openness (I(1))}
}
\examples{
# Generate panel data for 5 countries, 30 years
panel <- generate_panel_data(n_groups = 5, n_time = 30)
head(panel)

# Check dimensions
table(panel$country)
}
\seealso{
\code{\link{panel_ardl}}, \code{\link{generate_ts_data}}
}
\keyword{datagen}
