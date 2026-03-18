---
title: 'ardlverse: An R Ecosystem for ARDL Modeling and Cointegration Analysis'
tags:
  - R
  - econometrics
  - ARDL
  - cointegration
  - nonlinear models
  - quantile regression
  - panel data
  - Fourier approximation
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 18 March 2026
bibliography: paper.bib
---

# Summary

The `ardlverse` is an ecosystem of R packages for Autoregressive Distributed Lag (ARDL) modeling and bounds testing, spanning: augmented ARDL [@Sam2019], rolling/recursive ARDL [@Shahbaz2023], Fourier NARDL [@Bertelli2022], bootstrap NARDL [@McNown2018], nonlinear ARDL [@Shin2014], multiple-threshold NARDL [@Pal2016], quantile ARDL [@Cho2015], Fourier quantile ARDL, panel ARDL/NARDL, panel quantile ARDL, and cross-sectionally augmented panel quantile ARDL [@Pesaran2006cce]. All packages are open-source under GPL-3.

# Statement of Need

The ARDL bounds test of @PSS2001 is among the most widely used cointegration frameworks, yet existing R implementations cover only a subset. The `ARDL` package [@Natsiopoulos2022] provides standard linear ARDL but lacks bootstrap bounds testing, nonlinear decompositions, Fourier approximations, quantile extensions, and panel variants. No existing R package implements rolling/recursive ARDL, panel quantile ARDL, or cross-sectionally augmented panel ARDL.

# Packages

## aardl (ardlverse)

Augmented ARDL of @Sam2019, resolving degenerate PSS bounds test cases via F-overall, t-dependent, and F-independent tests. Supports linear, bootstrap, Fourier, and NARDL variants.

```r
library(ardlverse)
result <- aardl(y ~ x1 + x2, data = ts_data, p = 4, q = 4,
                case = 3, type = "bootstrap", nboot = 2000)
summary(result)
```

## rardl (ardlverse)

Rolling-window and recursive ARDL bounds testing [@Shahbaz2023; @Khan2023] for detecting structural instability in cointegrating relationships.

```r
result <- rardl(y ~ x1 + x2, data = ts_data,
                method = "rolling", window = 60, case = 3)
plot(result)
```

## fbnardl

Fourier NARDL with @Kripfganz2020 critical values and Fourier Bootstrap NARDL with bootstrap cointegration tests [@Bertelli2022].

```r
library(fbnardl)
result <- fbnardl(y ~ x1, data = ts_data, decompose = "x1",
                  type = "fbnardl", maxlag = 4, maxk = 3, reps = 999)
summary(result)
```

## mtnardl (ardlverse)

Multiple-threshold NARDL of @Pal2016, decomposing regressors into regime-specific partial sums.

```r
result <- mtnardl(y ~ x1, data = ts_data,
                  thresholds = c(-0.5, 0, 0.5), case = 3)
summary(result)
```

## fqardl

Fourier Quantile ARDL combining Fourier terms with quantile regression in the ARDL framework.

```r
library(fqardl)
result <- fqardl(y = ts_data$y, x = ts_data$x1,
                 tau = c(0.25, 0.5, 0.75), type = "fqardl")
summary(result)
```

## qardlr

Quantile ARDL of @Cho2015 with quantile-specific long-run, short-run, and impact parameters, BIC-based lag selection, and Wald tests for parameter constancy.

```r
library(qardlr)
result <- qardl(y ~ x1 + x2, data = ts_data,
                tau = seq(0.1, 0.9, 0.1), p = 2, q = 2)
summary(result)
```

## pnardl (ardlverse)

Panel NARDL [@Shin2014] with PMG, MG, and DFE estimators, capturing asymmetric effects via partial sum decompositions.

```r
result <- pnardl(y ~ x1, data = panel_data,
                 id = "country", time = "year",
                 p = 1, q = 1)
summary(result)
```

## xtpqardl

Panel Quantile ARDL with PMG, MG, and DFE estimators across multiple quantiles [@Cho2015; @Bildirici2022].

```r
library(xtpqardl)
result <- xtpqardl(y ~ x1 + x2, data = panel_data,
                   id = "country", time = "year", lr = c("x1", "x2"),
                   tau = c(0.25, 0.5, 0.75), model = "pmg")
summary(result)
```

## xtcspqardl

Cross-sectionally augmented panel quantile ARDL with QCCEMG and CS-PQARDL estimators for panel data with cross-sectional dependence [@Pesaran2006cce; @Chudik2015].

```r
library(xtcspqardl)
result <- xtcspqardl(y ~ x1 + x2, data = panel_data,
                     id = "country", time = "year",
                     tau = 0.5, estimator = "qccemg")
summary(result)
```

# Acknowledgements

The author acknowledges Merwan Roudane for contributing the original Stata and Python implementations that informed these packages.

# References
