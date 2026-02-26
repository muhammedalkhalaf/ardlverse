# ardlverse <img src="man/figures/logo.png" align="right" height="139" />
  
<!-- badges: start -->
[![R-CMD-check](https://github.com/muhammedalkhalaf/ardlverse/workflows/R-CMD-check/badge.svg)](https://github.com/muhammedalkhalaf/ardlverse/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/ardlverse)](https://CRAN.R-project.org/package=ardlverse)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ardlverse)](https://CRAN.R-project.org/package=ardlverse)
<!-- badges: end -->

## Overview

**ardlverse** is a comprehensive R package for Autoregressive Distributed Lag (ARDL) modeling and cointegration analysis. It provides unified tools for:

- 📊 **Panel ARDL**: Pooled Mean Group (PMG), Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators
- 🎰 **Bootstrap ARDL**: Bootstrap-based bounds testing for small samples
- 📈 **Quantile NARDL**: Quantile Nonlinear ARDL combining distributional and asymmetric effects
- 🌊 **Fourier ARDL**: Smooth structural breaks using Fourier approximation
- 🔍 **Diagnostics**: Comprehensive model diagnostics and visualization

## Installation

```r
# Install from CRAN (once available)
install.packages("ardlverse")

# Or install development version from GitHub
# install.packages("devtools")
devtools::install_github("muhammedalkhalaf/ardlverse")
```

## Quick Start

### Panel ARDL with PMG Estimator

```r
library(ardlverse)

# Generate example data
data <- generate_panel_data(n_groups = 10, n_time = 50)

# Estimate PMG model
pmg_model <- panel_ardl(
  gdp ~ inflation + investment,
  data = data,
  id = "country",
  time = "year",
  p = 1, q = 1,
  estimator = "pmg"
)

summary(pmg_model)

# Hausman test: PMG vs MG
hausman_test(pmg_model)
```

### Bootstrap Bounds Test

```r
# Generate time series data
ts_data <- generate_ts_data(n = 100)

# Bootstrap bounds test
boot_test <- boot_ardl(
  gdp ~ inflation + investment,
  data = ts_data,
  p = 2, q = 2,
  case = 3,
  nboot = 2000
)

summary(boot_test)
plot(boot_test)
```

### Quantile Nonlinear ARDL

```r
# Generate oil price data
oil <- generate_oil_data(n = 200)

# Estimate QNARDL
qnardl_model <- qnardl(
  gasoline ~ oil_price + exchange_rate,
  data = oil,
  tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
  p = 2, q = 2
)

summary(qnardl_model)
plot(qnardl_model, var = "oil_price")

# Test for asymmetry
asymmetry_test(qnardl_model, var = "oil_price")

# Dynamic multipliers
dynamic_multipliers(qnardl_model, var = "oil_price", tau = 0.5)
```

### Fourier ARDL

```r
# Estimate Fourier ARDL with automatic frequency selection
f_model <- fourier_ardl(
  gdp ~ investment + trade,
  data = ts_data,
  p = 2, q = 2,
  selection = "aic"
)

summary(f_model)
plot(f_model)
fourier_bounds_test(f_model)
```

### Model Diagnostics

```r
# Run comprehensive diagnostics
diag <- ardl_diagnostics(f_model)
summary(diag)
plot(diag)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `panel_ardl()` | Panel ARDL with PMG, MG, DFE estimators |
| `boot_ardl()` | Bootstrap ARDL bounds test |
| `qnardl()` | Quantile Nonlinear ARDL |
| `fourier_ardl()` | Fourier ARDL for structural breaks |
| `ardl_diagnostics()` | Comprehensive model diagnostics |
| `hausman_test()` | Hausman test for PMG vs MG |
| `asymmetry_test()` | Test for long-run asymmetry |
| `dynamic_multipliers()` | Cumulative dynamic multipliers |
| `pss_critical_values()` | PSS (2001) critical value tables |

## Theoretical Background

### Panel ARDL (Pesaran, Shin & Smith, 1999)

The PMG estimator allows for heterogeneous short-run dynamics while constraining long-run coefficients to be equal across groups:

$$\Delta y_{it} = \phi_i (y_{i,t-1} - \theta' x_{it}) + \sum_{j=1}^{p-1} \lambda_{ij} \Delta y_{i,t-j} + \sum_{j=0}^{q-1} \delta'_{ij} \Delta x_{i,t-j} + \mu_i + \varepsilon_{it}$$

### QNARDL (Cho, Kim & Shin, 2015 + Shin, Yu & Greenwood-Nimmo, 2014)

Combines quantile regression with asymmetric decomposition:

$$x^+_t = \sum_{j=1}^{t} \max(\Delta x_j, 0), \quad x^-_t = \sum_{j=1}^{t} \min(\Delta x_j, 0)$$

### Fourier ARDL (Banerjee, Arcabic & Lee, 2017)

Captures smooth structural breaks using Fourier approximation:

$$f_t = \sum_{k=1}^{K} [a_k \sin(2\pi k t/T) + b_k \cos(2\pi k t/T)]$$

## References

- Pesaran, M. H., Shin, Y., & Smith, R. P. (1999). Pooled mean group estimation of dynamic heterogeneous panels. *Journal of the American Statistical Association*, 94(446), 621-634.

- Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. *Journal of Applied Econometrics*, 16(3), 289-326.

- Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework. In *Festschrift in Honor of Peter Schmidt* (pp. 281-314). Springer.

- Cho, J. S., Kim, T. H., & Shin, Y. (2015). Quantile cointegration in the autoregressive distributed-lag modeling framework. *Journal of Econometrics*, 188(1), 281-300.

- Banerjee, P., Arcabic, V., & Lee, H. (2017). Fourier ADL cointegration test to approximate smooth breaks with new evidence from crude oil market. *Economic Modelling*, 67, 114-124.

- McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive distributed lag test for cointegration. *Applied Economics*, 50(13), 1509-1521.

## Author

**Muhammad Alkhalaf**
- ORCID: [0009-0002-2677-9246](https://orcid.org/0009-0002-2677-9246)
- Email: contact@rufyqelngeh.com
- Website: [rufyqelngeh.com](https://www.rufyqelngeh.com)

## License

GPL-3
