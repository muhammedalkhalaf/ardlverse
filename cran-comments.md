## ardlverse 2.0.3

This release corrects the computations below; if the 2.0.2 submission is still pending, it should be discarded in favour of this one.

* Critical values of the PSS bounds test: `pss_critical_values()` returned hand-typed bounds that did not match Pesaran, Shin and Smith (2001) for several cases. It now computes the bounds from the response surface regressions of Kripfganz and Schneider (2020), the same coefficients as the Stata command `ardlbounds`, with finite-sample values when `n` and `sr` are given. The new internal `.pss_pvalues()` gives approximate p-values. Reference values from `ardlbounds` are checked in the tests.
* `aardl()`: the critical values were produced by invented formulas. The F and t bounds now come from the Kripfganz and Schneider (2020) response surfaces for the sample size and number of short-run coefficients of the model; the F test on the lagged independent variables has no tabulated bounds and is left to the bootstrap.
* `boot_ardl()`, `aardl()`, `mtnardl()`, `rardl()` and `fourier_ardl()`: in case 2 (restricted intercept) and case 4 (restricted trend) the F statistic did not restrict the intercept or the trend, contrary to Pesaran, Shin and Smith (2001). The restriction is now included, and the bootstrap samples are generated from the same restricted model. `boot_ardl()` also omitted the intercept in case 2.
* `boot_ardl()`: the bootstrap samples were generated from the unrestricted model, so the bootstrap distribution was not taken under the null of no level relationship. It is now a fixed-regressor residual bootstrap under the restricted model.
* `rardl()`: the bounds are computed for each window's sample size.
* `fourier_bounds_test()`: the critical values were invented formulas. No bounds have been tabulated for the Fourier ARDL F test, so the function now reports no decision; the Kripfganz and Schneider (2020) bounds for the model without Fourier terms are returned for reference only. `fourier_ardl()` also added an implicit intercept to every case, including case 1, and treated case 2 as having no intercept; the deterministic terms are now explicit.
* `qnardl()`: the lagged differences were misaligned by one period and the partial-sum levels were contemporaneous rather than lagged. The long-run asymmetry Wald test ignored the covariances between the coefficients; it now uses the delta method with the full covariance matrix of the quantile regression.
* `pnardl()`: the long-run asymmetry test ignored the covariance between the positive and negative coefficients, and the Hausman test replaced missing or negative variances with 0.01 and 0.001. The asymmetry test now uses the full covariance matrix (mean group: covariance of the unit estimates divided by N; DFE: delta method), and the Hausman test uses a generalised inverse of the covariance difference, returning NA when that difference has no positive eigenvalue.
* `ardl_diagnostics()`: the CUSUM bound used the recursive-residual constant 0.948 with a linear boundary, although the test is computed from OLS residuals; it is now the constant 1.358 sqrt(n) of Ploberger and Kraemer (1992). The CUSUM of squares bound is scaled by the estimated kurtosis of the residuals (Deng and Perron, 2008).

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
