\name{generate_oil_data}
\alias{generate_oil_data}
\title{Generate Oil Price Data}
\description{
Generate simulated oil and gasoline price data with asymmetric effects.
}
\usage{
generate_oil_data(n = 200, seed = 789)
}
\arguments{
  \item{n}{Number of observations (default: 200)}
  \item{seed}{Random seed for reproducibility (default: 789)}
}
\details{
Generates a simulated dataset exhibiting the "rockets and feathers"
phenomenon where gasoline prices respond more quickly to oil price
increases than decreases. Suitable for demonstrating QNARDL methods
and asymmetry testing.
}
\value{
A data frame with columns:
\item{week}{Week number (1 to n)}
\item{gasoline}{Retail gasoline price (cents/gallon)}
\item{oil_price}{Crude oil price ($/barrel)}
\item{exchange_rate}{USD exchange rate index}
}
\examples{
# Generate oil price data
oil <- generate_oil_data(n = 200)
head(oil)

# Plot prices
plot(oil$week, oil$oil_price, type = "l", col = "blue",
     main = "Oil and Gasoline Prices")
lines(oil$week, oil$gasoline/4, col = "red")  # Scaled for comparison
}
\seealso{
\code{\link{qnardl}}, \code{\link{asymmetry_test}}
}
\keyword{datagen}
