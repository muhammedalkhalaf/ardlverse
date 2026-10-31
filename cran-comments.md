## ardlverse 2.1.2

CRAN version: 2.1.1, published 2026-10-10 (the `Date/Publication` field of the CRAN-installed DESCRIPTION in the baseline source; the source `Date` field is 2026-10-16).

Changes in this version (all in `hausman_test()` and its documentation; `panel_ardl()` estimation is unchanged):

* `hausman_test()` stops with an error when the two fits differ in dependent variable, regressors, lag orders, trend, id or time variable, number of observations or number of groups.
* New argument `quiet` (default `FALSE`) suppresses printing.
* The returned list gains `scale` (factor applied to the MG covariance matrix).
* When N - 1 < k the warning and the printed output state that the MG covariance matrix has rank at most N - 1 so the test cannot be defined.
* Documentation of `sigmamore`: exact formula (sigma_eff / sigma_MG)^4 on root mean squared residuals. The fourth power is the Stata/xtpmg convention; Hausman (1978, eqs. 2.16-2.17, pp. 1258 and 1260) forms V(q) with one common sigma^2, which corresponds to the second power, so the fourth power is not in Hausman (1978). For MG versus DFE the factor is a package choice.
* Documentation of `method = "ginv"`: described as a package choice; the attribution to Stata's `hausman` removed.
* Documentation: statistic, degrees of freedom and p-value written out; Monte Carlo evidence on undefined tests and size at N = 5, 10 and 20 added.
* README example of `hausman_test()` corrected to the MG-first argument order.
* Tests added for the new checks and a hand quadratic form on `lm()`-rebuilt inputs at 1e-10.
* One ASCII arrow in a code comment replaced.

## Test environments

* Local: R 4.3.3 on Linux (x86_64), `R CMD check --as-cran --no-manual` with `_R_CHECK_SYSTEM_CLOCK_=0`. No R-devel was available in the build container, so the package was checked on R 4.3.3 only.

## R CMD check results

0 errors, 0 warnings, 1 NOTE.

* NOTE (CRAN incoming feasibility): "Found the following (possibly) invalid URLs" for https://github.com/muhammedalkhalaf/ardlverse, .../issues and .../actions with status 403. The build container reaches GitHub only through a proxy that answers 403 to the check's HEAD requests; the URLs are valid.

The DESCRIPTION fields that CRAN writes at publication (Author, Maintainer, Packaged, Repository, Date/Publication) were removed from the source so that `R CMD build` derives them from `Authors@R`.
