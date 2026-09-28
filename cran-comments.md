## ardlverse 2.0.1

* Bug fix in `rardl()`, `aardl()` and `mtnardl()`: the lagged differences of the dependent variable were misaligned by one period (lag i was built as lag i-1). With `p = 1` the first "lagged" difference was identical to the regressand, so every regression had a perfect fit; `rardl()` returned NA for all windows and its plot method failed. The sample now starts at `max(p + 1, max(q)) + 1` and lag i is built as lag i.
* Bug fix in `rardl()` and `mtnardl()`: the bounds decision read `cv$F_I1["5%"]` although `pss_critical_values()` returns `cv$F_bounds$I1`; the comparison therefore failed, which made every `rardl()` window an error and every asymptotic `mtnardl()` decision "INCONCLUSIVE". Both now use the returned bounds.
* Authors@R updated: the package is developed and maintained by Muhammad Alkhalaf; a former contributor entry was removed.

All DOIs in the package were verified against CrossRef before this submission.

The two code fixes were validated against an independent `lm()` fit of the same
conditional error-correction regression (identical coefficients and standard errors)
and the `\donttest` examples of `rardl()`, `aardl()` and `mtnardl()` now run cleanly.

## Test environments

* Ubuntu 24.04, R 4.3.3 (R CMD check --as-cran)
* CRAN check results for the previous version: OK on all platforms

## R CMD check results

0 errors | 0 warnings | 0 notes
