\name{dynamic_multipliers}
\alias{dynamic_multipliers}
\title{Dynamic Multipliers for QNARDL}
\description{
Compute and plot cumulative dynamic multipliers for QNARDL models.
}
\usage{
dynamic_multipliers(object, var, tau = 0.5, horizon = 20)
}
\arguments{
  \item{object}{A qnardl object}
  \item{var}{Variable name}
  \item{tau}{Quantile to compute multipliers for (default: 0.5)}
  \item{horizon}{Number of periods (default: 20)}
}
\details{
Computes the cumulative dynamic multipliers showing how the effect of
a unit change in x accumulates over time towards the long-run equilibrium.

Separate multipliers are computed for positive and negative changes,
allowing visualization of asymmetric adjustment paths.
}
\value{
A data frame with:
\item{horizon}{Time horizon (0 to horizon)}
\item{multiplier}{Cumulative multiplier value}
\item{type}{Either "Positive" or "Negative"}
}
\examples{
\dontrun{
data <- generate_oil_data(n = 200)

model <- qnardl(
  gasoline ~ oil_price,
  data = data,
  tau = c(0.25, 0.5, 0.75)
)

# Plot dynamic multipliers at median
dynamic_multipliers(model, var = "oil_price", tau = 0.5, horizon = 30)
}
}
\seealso{
\code{\link{qnardl}}, \code{\link{asymmetry_test}}
}
\keyword{regression}
