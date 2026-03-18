---
title: 'ardlverse: A Comprehensive R Ecosystem for ARDL Modeling and Cointegration Analysis'
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

We introduce the `ardlverse`, an integrated ecosystem of eleven R packages for Autoregressive Distributed Lag (ARDL) modeling and bounds testing for cointegration. The ecosystem spans the full methodological frontier of ARDL analysis: augmented ARDL [@Sam2019], rolling and recursive ARDL [@Shahbaz2023], Fourier ARDL [@Yilanci2020], bootstrap ARDL [@McNown2018], nonlinear ARDL [@Shin2014], multiple-threshold NARDL [@Pal2016], quantile ARDL [@Cho2015], panel ARDL and panel NARDL, panel quantile ARDL, and cross-sectionally augmented panel quantile ARDL [@Harding2018]. The meta-package `ardlverse` provides a unified interface to all methods. All packages are open-source under GPL-3.

# Statement of Need

The ARDL bounds testing approach of @PSS2001 has become one of the most widely used cointegration frameworks in applied economics, with over 20,000 citations. Its appeal lies in its flexibility: it accommodates mixtures of I(0) and I(1) variables without requiring pre-testing for integration order. However, the methodology has expanded considerably since 2001, and existing R implementations cover only a subset of available methods.

The `ARDL` package [@Natsiopoulos2022] provides standard linear ARDL estimation, but does not implement bootstrap bounds testing, nonlinear decompositions, Fourier approximations, quantile extensions, or panel variants. The `nardl` package offers basic nonlinear ARDL but lacks the augmented framework, multiple thresholds, and Fourier terms. No existing R package implements panel quantile ARDL, cross-sectionally augmented panel ARDL, or the rolling/recursive ARDL framework.

Our ecosystem addresses these gaps comprehensively, providing researchers with a single, coherent collection of packages covering the full range of ARDL methods currently available in the literature.

# Packages

## aardl: Augmented ARDL

Implements the augmented ARDL (A-ARDL) framework of @Sam2019, which resolves degenerate cases of the standard PSS bounds test through a three-test procedure: the F-overall test, the t-dependent variable test, and the F-independent variables test. Includes eight model variants combining standard/bootstrap/Fourier specifications with linear/nonlinear decompositions.

```r
library(aardl)
result <- aardl(y ~ x1 + x2, data = ts_data, case = 3,
                max_p = 4, max_q = 4, ic = "aic")
summary(result)
```

## rardl: Rolling and Recursive ARDL

Implements rolling-window and recursive ARDL bounds testing following @Shahbaz2023 and @Khan2023, enabling detection of structural instability in cointegrating relationships over time. Also provides recursive ADF unit root testing and recursive Granger causality.

```r
library(rardl)
result <- rardl(y ~ x1 + x2, data = ts_data, type = "rolling",
                window = 60, case = 3)
plot(result)
```

## fbnardl: Fourier Bootstrap Nonlinear ARDL

Provides the Fourier NARDL (FNARDL) with @Kripfganz2020 critical values and the Fourier Bootstrap NARDL (FBNARDL) with bootstrap cointegration tests following @Bertelli2022. Combines Fourier terms for smooth structural breaks with partial sum decompositions for asymmetric effects.

```r
library(fbnardl)
result <- fbnardl(y ~ x1, data = ts_data, case = 3,
                  max_freq = 3, nboot = 1000)
summary(result)
```

## mtnardl: Multiple-Threshold Nonlinear ARDL

Implements the MTNARDL model of @Pal2016, decomposing regressors into regime-specific partial sums based on quantile or custom cut-point partitions. Provides dynamic multipliers per regime and optional bootstrap bounds testing.

```r
library(mtnardl)
result <- mtnardl(y ~ x1, data = ts_data, thresholds = 3,
                  case = 3, nboot = 1000)
summary(result)
```

## fqardl: Fourier ARDL Methods

A comprehensive package implementing Fourier Quantile ARDL (FQARDL), Fourier Nonlinear ARDL (FNARDL), and Multi-Threshold NARDL (MTNARDL), along with Fourier unit root tests following @Enders2012 and @Becker2006. Features automatic lag and frequency selection with publication-ready visualizations.

```r
library(fqardl)
result <- fqardl(y ~ x1, data = ts_data,
                 taus = c(0.25, 0.5, 0.75), max_freq = 3)
summary(result)
```

## qardlr: Quantile ARDL

Implements the Quantile ARDL (QARDL) model of @Cho2015, estimating quantile-specific long-run, short-run autoregressive, and impact parameters. Features BIC-based lag selection, ECM parameterization, Wald tests for parameter constancy across quantiles, and rolling/recursive QARDL estimation.

```r
library(qardlr)
result <- qardl(y ~ x1 + x2, data = ts_data,
                taus = seq(0.1, 0.9, 0.1), p = 2, q = 2)
summary(result)
```

## pnardl: Panel Nonlinear ARDL

Implements the Panel NARDL model following @Shin2014, with Pooled Mean Group (PMG), Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators. Captures asymmetric long-run and short-run effects through partial sum decompositions. Includes asymmetry Wald tests, dynamic multipliers, and impulse response functions.

```r
library(pnardl)
result <- pnardl(y ~ x1 + x2, data = panel_data,
                 id = "country", time = "year", estimator = "pmg")
summary(result)
```

## xtpqardl: Panel Quantile ARDL

Estimates Panel Quantile ARDL (PQARDL) models combining panel ARDL methodology with quantile regression. Supports PMG, MG, and DFE estimators across multiple quantiles following @PSS1999, @Cho2015, and @Bildirici2022.

```r
library(xtpqardl)
result <- xtpqardl(y ~ x1 + x2, data = panel_data,
                   id = "country", time = "year",
                   taus = c(0.25, 0.5, 0.75), estimator = "pmg")
summary(result)
```

## xtcspqardl: Cross-Sectionally Augmented Panel Quantile ARDL

Implements the CS-PQARDL model and Quantile Common Correlated Effects Mean Group (QCCEMG) estimator for panel data with cross-sectional dependence. Handles unobserved common factors through cross-sectional averages following @Pesaran2006cce and @Chudik2015.

```r
library(xtcspqardl)
result <- cspqardl(y ~ x1 + x2, data = panel_data,
                   id = "country", time = "year",
                   taus = c(0.25, 0.5, 0.75))
summary(result)
```

## ardlverse: Meta-Package

The `ardlverse` package provides a unified framework that loads all ARDL packages and offers a common interface for model selection, estimation, and diagnostics. Includes comprehensive visualization tools and publication-ready output formatting.

```r
library(ardlverse)
# Access all ARDL methods through a unified interface
result <- ardl_estimate(y ~ x1 + x2, data = ts_data,
                        method = "bootstrap", case = 3)
summary(result)
```

# Acknowledgements

The author acknowledges Merwan Roudane for contributing the original Stata and Python implementations that informed the `fqardl` and `ardlverse` packages.

# References
