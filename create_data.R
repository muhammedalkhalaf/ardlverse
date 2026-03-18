# Create datasets for ardlverse
setwd("C:/Users/acad_/.openclaw/workspace/ardlverse")
set.seed(123)

# Generate panel data
n_groups <- 10
n_time <- 50
data_list <- lapply(1:n_groups, function(i) {
  alpha_i <- rnorm(1, mean = 0, sd = 0.5)
  x1 <- cumsum(rnorm(n_time, sd = 0.5))
  x2 <- cumsum(rnorm(n_time, sd = 0.3))
  y_lr <- 5 + 0.3 * x1 + 0.2 * x2
  y <- rep(0, n_time)
  y[1] <- y_lr[1] + rnorm(1, sd = 0.2)
  for (t in 2:n_time) {
    ec <- y[t-1] - y_lr[t-1]
    dy <- -0.2 * ec + 0.1 * (x1[t] - x1[t-1]) + 0.05 * (x2[t] - x2[t-1]) + alpha_i + rnorm(1, sd = 0.1)
    y[t] <- y[t-1] + dy
  }
  inflation <- 2 + 0.5 * diff(c(0, y)) + rnorm(n_time, sd = 0.3)
  data.frame(country = i, year = 1970:(1970 + n_time - 1), gdp = y, inflation = inflation, investment = 20 + x1, trade = 50 + x2 * 10)
})
macro_panel <- do.call(rbind, data_list)

# Generate time series data
n <- 100
set.seed(456)
x1 <- cumsum(rnorm(n, sd = 0.3))
x2 <- cumsum(rnorm(n, sd = 0.2))
t_idx <- 1:n
fourier <- 0.5 * sin(2 * pi * t_idx / n) + 0.3 * cos(2 * pi * t_idx / n)
y_lr <- 3 + 0.4 * x1 + 0.25 * x2 + fourier
y <- rep(0, n)
y[1] <- y_lr[1]
for (i in 2:n) {
  ec <- y[i-1] - y_lr[i-1]
  dy <- -0.15 * ec + 0.2 * (x1[i] - x1[i-1]) + 0.1 * (x2[i] - x2[i-1]) + rnorm(1, sd = 0.05)
  y[i] <- y[i-1] + dy
}
inflation <- 2 + 0.3 * diff(c(0, y)) + rnorm(n, sd = 0.2)
macro_data <- data.frame(quarter = 1:n, gdp = y, inflation = inflation, investment = 5 + x1, trade = 2 + x2)

# Generate oil data
n <- 200
set.seed(789)
d_oil <- rnorm(n, mean = 0.05, sd = 2)
oil <- 50 + cumsum(d_oil)
d_oil_pos <- pmax(d_oil, 0)
d_oil_neg <- pmin(d_oil, 0)
oil_pos <- cumsum(d_oil_pos)
oil_neg <- cumsum(d_oil_neg)
exchange <- 100 + cumsum(rnorm(n, sd = 0.5))
gasoline <- rep(0, n)
gasoline[1] <- 200 + 2 * oil[1]
for (i in 2:n) {
  ec <- gasoline[i-1] - (150 + 1.5 * oil_pos[i-1] + 1.0 * oil_neg[i-1])
  dgas <- -0.1 * ec + 0.8 * d_oil_pos[i] + 0.3 * d_oil_neg[i] + 0.05 * (exchange[i] - exchange[i-1]) + rnorm(1, sd = 1)
  gasoline[i] <- gasoline[i-1] + dgas
}
oil_data <- data.frame(week = 1:n, gasoline = gasoline, oil_price = oil, exchange_rate = exchange)

# Save
if (!dir.exists("data")) dir.create("data")
save(macro_panel, file = "data/macro_panel.rda", compress = "xz")
save(macro_data, file = "data/macro_data.rda", compress = "xz")
save(oil_data, file = "data/oil_data.rda", compress = "xz")
cat("Datasets created!\n")
