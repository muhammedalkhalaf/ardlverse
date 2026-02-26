#' @title Simulated Macroeconomic Panel Data
#' @description A simulated panel dataset of macroeconomic variables
#'   for 10 countries over 50 years.
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{country}{Country identifier (1-10)}
#'   \item{year}{Year (1970-2019)}
#'   \item{gdp}{Log GDP per capita}
#'   \item{inflation}{Inflation rate}
#'   \item{investment}{Gross capital formation (% of GDP)}
#'   \item{trade}{Trade openness (exports + imports as % of GDP)}
#' }
#'
#' @examples
#' \dontrun{
#' data(macro_panel)
#' head(macro_panel)
#'
#' # Estimate Panel ARDL
#' model <- panel_ardl(
#'   gdp ~ inflation + investment,
#'   data = macro_panel,
#'   id = "country",
#'   time = "year",
#'   estimator = "pmg"
#' )
#' }
#'
#' @source Simulated data for demonstration purposes
"macro_panel"


#' @title Time Series Macroeconomic Data
#' @description A simulated time series dataset of macroeconomic variables
#'   for a single country over 100 quarters.
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{quarter}{Quarter number (1-100)}
#'   \item{gdp}{Log real GDP}
#'   \item{inflation}{Inflation rate}
#'   \item{investment}{Investment growth rate}
#'   \item{trade}{Trade balance (% of GDP)}
#' }
#'
#' @examples
#' \dontrun{
#' data(macro_data)
#'
#' # Bootstrap bounds test
#' test <- boot_ardl(
#'   gdp ~ inflation + investment,
#'   data = macro_data,
#'   p = 2, q = 2
#' )
#' }
#'
#' @source Simulated data for demonstration purposes
"macro_data"


#' @title Oil Price and Gasoline Data
#' @description Weekly data on oil prices and retail gasoline prices,
#'   suitable for demonstrating asymmetric price transmission.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{week}{Week number (1-200)}
#'   \item{gasoline}{Retail gasoline price (cents/gallon)}
#'   \item{oil_price}{Crude oil price ($/barrel)}
#'   \item{exchange_rate}{USD exchange rate index}
#' }
#'
#' @examples
#' \dontrun{
#' data(oil_data)
#'
#' # Test for asymmetric price transmission
#' model <- qnardl(
#'   gasoline ~ oil_price + exchange_rate,
#'   data = oil_data,
#'   tau = c(0.25, 0.5, 0.75)
#' )
#' asymmetry_test(model, var = "oil_price")
#' }
#'
#' @source Simulated data for demonstration purposes
"oil_data"


#' @title Generate Example Panel Data
#' @description Generate simulated panel data for examples
#'
#' @param n_groups Number of groups/countries
#' @param n_time Number of time periods
#' @param seed Random seed
#'
#' @return A data frame with panel structure
#'
#' @examples
#' # Generate 5 countries, 30 years
#' panel <- generate_panel_data(n_groups = 5, n_time = 30)
#' head(panel)
#'
#' @export
generate_panel_data <- function(n_groups = 10, n_time = 50, seed = 123) {
  
  set.seed(seed)
  
  # Generate data
  data_list <- lapply(1:n_groups, function(i) {
    
    # Country-specific parameters
    alpha_i <- rnorm(1, mean = 0, sd = 0.5)
    
    # Generate cointegrated series
    x1 <- cumsum(rnorm(n_time, sd = 0.5))  # Investment
    x2 <- cumsum(rnorm(n_time, sd = 0.3))  # Trade
    
    # Long-run relationship
    y_lr <- 5 + 0.3 * x1 + 0.2 * x2
    
    # Add error correction dynamics
    ec <- rep(0, n_time)
    y <- rep(0, n_time)
    y[1] <- y_lr[1] + rnorm(1, sd = 0.2)
    
    for (t in 2:n_time) {
      ec[t] <- y[t-1] - y_lr[t-1]
      dy <- -0.2 * ec[t] + 0.1 * (x1[t] - x1[t-1]) + 0.05 * (x2[t] - x2[t-1]) + 
            alpha_i + rnorm(1, sd = 0.1)
      y[t] <- y[t-1] + dy
    }
    
    # Inflation (I(0))
    inflation <- 2 + 0.5 * diff(c(0, y)) + rnorm(n_time, sd = 0.3)
    
    data.frame(
      country = i,
      year = 1970:(1970 + n_time - 1),
      gdp = y,
      inflation = inflation,
      investment = 20 + x1,
      trade = 50 + x2 * 10
    )
  })
  
  do.call(rbind, data_list)
}


#' @title Generate Example Time Series Data
#' @description Generate simulated time series data for examples
#'
#' @param n Number of observations
#' @param seed Random seed
#'
#' @return A data frame with time series
#'
#' @examples
#' ts_data <- generate_ts_data(n = 100)
#' head(ts_data)
#'
#' @export
generate_ts_data <- function(n = 100, seed = 456) {
  
  set.seed(seed)
  
  # Generate I(1) series
  x1 <- cumsum(rnorm(n, sd = 0.3))  # Investment
  x2 <- cumsum(rnorm(n, sd = 0.2))  # Trade
  
  # Cointegrated dependent variable with structural break
  # Fourier component for smooth break
  t <- 1:n
  fourier <- 0.5 * sin(2 * pi * t / n) + 0.3 * cos(2 * pi * t / n)
  
  y_lr <- 3 + 0.4 * x1 + 0.25 * x2 + fourier
  
  # ECM dynamics
  y <- rep(0, n)
  y[1] <- y_lr[1]
  
  for (i in 2:n) {
    ec <- y[i-1] - y_lr[i-1]
    dy <- -0.15 * ec + 0.2 * (x1[i] - x1[i-1]) + 0.1 * (x2[i] - x2[i-1]) + 
          rnorm(1, sd = 0.05)
    y[i] <- y[i-1] + dy
  }
  
  # Inflation (stationary)
  inflation <- 2 + 0.3 * diff(c(0, y)) + rnorm(n, sd = 0.2)
  
  data.frame(
    quarter = 1:n,
    gdp = y,
    inflation = inflation,
    investment = 5 + x1,
    trade = 2 + x2
  )
}


#' @title Generate Oil Price Data
#' @description Generate simulated oil and gasoline price data with asymmetric effects
#'
#' @param n Number of observations
#' @param seed Random seed
#'
#' @return A data frame
#'
#' @examples
#' oil <- generate_oil_data(n = 200)
#' head(oil)
#'
#' @export
generate_oil_data <- function(n = 200, seed = 789) {
  
  set.seed(seed)
  
  # Oil price (random walk with drift)
  d_oil <- rnorm(n, mean = 0.05, sd = 2)
  oil <- 50 + cumsum(d_oil)
  
  # Decompose into positive and negative changes
  d_oil_pos <- pmax(d_oil, 0)
  d_oil_neg <- pmin(d_oil, 0)
  
  oil_pos <- cumsum(d_oil_pos)
  oil_neg <- cumsum(d_oil_neg)
  
  # Exchange rate
  exchange <- 100 + cumsum(rnorm(n, sd = 0.5))
  
  # Gasoline with asymmetric response (rockets and feathers)
  # Positive oil changes have larger effect than negative
  gasoline <- rep(0, n)
  gasoline[1] <- 200 + 2 * oil[1]
  
  for (i in 2:n) {
    ec <- gasoline[i-1] - (150 + 1.5 * oil_pos[i-1] + 1.0 * oil_neg[i-1])
    
    # Asymmetric short-run: positive changes pass through faster
    dgas <- -0.1 * ec + 0.8 * d_oil_pos[i] + 0.3 * d_oil_neg[i] + 
            0.05 * (exchange[i] - exchange[i-1]) + rnorm(1, sd = 1)
    
    gasoline[i] <- gasoline[i-1] + dgas
  }
  
  data.frame(
    week = 1:n,
    gasoline = gasoline,
    oil_price = oil,
    exchange_rate = exchange
  )
}
