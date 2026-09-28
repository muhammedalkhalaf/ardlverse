# ardlverse 2.0.3

* Critical values of the PSS bounds test: `pss_critical_values()` returned hand-typed bounds that did not match Pesaran, Shin and Smith (2001) for several cases. It now computes the bounds from the response surface regressions of Kripfganz and Schneider (2020), the same coefficients as the Stata command `ardlbounds`, with finite-sample values when `n` and `sr` are given. The new internal `.pss_pvalues()` gives approximate p-values. Reference values from `ardlbounds` are checked in the tests.
* `aardl()`: the critical values were produced by invented formulas. The F and t bounds now come from the Kripfganz and Schneider (2020) response surfaces for the sample size and number of short-run coefficients of the model; the F test on the lagged independent variables has no tabulated bounds and is left to the bootstrap.
* `boot_ardl()`, `aardl()`, `mtnardl()`, `rardl()` and `fourier_ardl()`: in case 2 (restricted intercept) and case 4 (restricted trend) the F statistic did not restrict the intercept or the trend, contrary to Pesaran, Shin and Smith (2001). The restriction is now included, and the bootstrap samples are generated from the same restricted model. `boot_ardl()` also omitted the intercept in case 2.
* `boot_ardl()`: the bootstrap samples were generated from the unrestricted model, so the bootstrap distribution was not taken under the null of no level relationship. It is now a fixed-regressor residual bootstrap under the restricted model.
* `rardl()`: the bounds are computed for each window's sample size.
* `fourier_bounds_test()`: the critical values were invented formulas. No bounds have been tabulated for the Fourier ARDL F test, so the function now reports no decision; the Kripfganz and Schneider (2020) bounds for the model without Fourier terms are returned for reference only. `fourier_ardl()` also added an implicit intercept to every case, including case 1, and treated case 2 as having no intercept; the deterministic terms are now explicit.
* `qnardl()`: the lagged differences were misaligned by one period and the partial-sum levels were contemporaneous rather than lagged. The long-run asymmetry Wald test ignored the covariances between the coefficients; it now uses the delta method with the full covariance matrix of the quantile regression.
* `pnardl()`: the long-run asymmetry test ignored the covariance between the positive and negative coefficients, and the Hausman test replaced missing or negative variances with 0.01 and 0.001. The asymmetry test now uses the full covariance matrix (mean group: covariance of the unit estimates divided by N; DFE: delta method), and the Hausman test uses a generalised inverse of the covariance difference, returning NA when that difference has no positive eigenvalue.
* `ardl_diagnostics()`: the CUSUM bound used the recursive-residual constant 0.948 with a linear boundary, although the test is computed from OLS residuals; it is now the constant 1.358 sqrt(n) of Ploberger and Kraemer (1992). The CUSUM of squares bound is scaled by the estimated kurtosis of the residuals (Deng and Perron, 2008).

# ardlverse 2.0.2

* Bug fix in `pnardl()` (all three estimators): the same lag misalignment as in rardl(), aardl() and mtnardl() below; with `p = 1` the lagged difference of y was one observation short, `cbind()` recycled it with a warning on every unit and every bootstrap draw, and the design matrix was misaligned. The sample now starts at `max(p + 1, max(q)) + 1`.
* Bug fix in `pnardl(bootstrap = TRUE)`: the bootstrap standard errors (error-correction and long-run coefficients) were subtracted from the full short-run coefficient vector, so `ci_lower`/`ci_upper` were recycled across unrelated coefficients. The intervals are now attached to exactly the bootstrapped estimates, with names.
* The `pnardl()` example uses a smaller panel and `nboot = 20`; it ran for more than ten minutes in the CRAN incoming check.
* Bug fix in `rardl()`, `aardl()` and `mtnardl()`: the lagged differences of the dependent variable were misaligned by one period (lag i was built as lag i-1). With `p = 1` the first "lagged" difference was identical to the regressand, so every regression had a perfect fit; `rardl()` returned NA for all windows and its plot method failed. The sample now starts at `max(p + 1, max(q)) + 1` and lag i is built as lag i.
* Bug fix in `rardl()` and `mtnardl()`: the bounds decision read `cv$F_I1["5%"]` although `pss_critical_values()` returns `cv$F_bounds$I1`; the comparison therefore failed, which made every `rardl()` window an error and every asymptotic `mtnardl()` decision "INCONCLUSIVE". Both now use the returned bounds.
* Authors@R updated: the package is developed and maintained by Muhammad Alkhalaf; a former contributor entry was removed.

# ardlverse 2.0.0

## Major bug fixes in panel_ardl()

Thanks to Yeleazar (Lazar) Levchenko (Kyiv School of Economics), who
audited `panel_ardl()` against Stata's `xtpmg` (Blackburne & Frank 2007)
and contributed corrections that bring the implementation into strict
alignment with the original Pesaran, Shin & Smith (1999) framework.

Seven issues were identified and fixed:

1.  **Missing intercepts in short-run regressions.** The original code
    used `lm.fit()` for internal regressions; unlike `lm()`, `lm.fit()`
    does not append an intercept. All short-run regressions across PMG,
    MG, and DFE were forced through the origin. A column of 1s is now
    bound to the design matrices, and DFE reconstructs the grand-mean
    intercept to match standard fixed-effects output.

2.  **Misaligned error-correction term in `.prepare_ardl_data`.** The
    long-run matrix (`X_levels`) was constructed from rows 1 to (n-1),
    pairing the lagged dependent variable y(t-1) with lagged X(t-1).
    The standard ARDL error-correction term requires y(t-1) paired with
    contemporaneous X(t). Indexing corrected to rows 2 through n.

3.  **Statistically invalid Hausman test.** The previous test isolated
    only diagonal variances and used `abs()` to force-ignore negative
    variance differences, bypassing the covariance structure and
    invalidating the chi-squared statistic. The test is now built on the
    proper matrix quadratic form, with a new `sigmamore = TRUE` argument
    (matching Stata) that rescales the inefficient variance matrix when
    the difference matrix is non-positive-definite.

4.  **Incorrect PMG standard errors.** The previous code computed PMG
    SEs from the cross-sectional standard deviation of group-specific
    long-run estimates, contradicting PMG theory (long-run coefficients
    are constrained to be homogeneous). Replaced with the exact PSS
    (1999) Information Matrix formulation using the G-matrix blocks.

5.  **Simplified delta method for DFE standard errors.** The previous
    code assumed zero covariance between short-run coefficients and the
    error-correction parameter. The full multivariate delta method with
    the proper Jacobian is now used.

6.  **Incorrect MG standard errors.** Previously computed naively as
    SD / sqrt(N). Replaced with the exact cross-sectional
    variance-covariance formula used by Blackburne & Frank (2007).

7.  **Sub-optimal PMG initialization.** The previous code ran the full
    MG estimator to generate PMG starting values. PMG is now initialized
    from a simple pooled OLS of the lagged dependent variable on the
    levels of X — faster, avoids convergence risk if MG fails, and
    matches Stata's exact initialization.

## Breaking changes

*   **`hausman_test()` signature changed** to follow Stata's convention:

    *   Old: `hausman_test(pmg_model, mg_model, data)`
    *   New: `hausman_test(inefficient, efficient, sigmamore = TRUE)`

    Pass the inefficient (always consistent) estimator first (typically
    MG), then the efficient one (typically PMG). The previous third
    argument `data` is removed; required information is read from the
    model objects.

*   Internal helpers `.estimate_mg_internal()` and `.compute_pmg_se()`
    have been consolidated into `.estimate_mg()` and `.estimate_pmg()`
    respectively. Code that imported these internal functions (which is
    not supported usage) will need to be updated.

## Other changes

*   New optional arguments to `panel_ardl()`: `start_time` (restrict
    estimation to observations at or after a given time) and `cluster`
    (cluster-robust SEs where applicable).
*   Convergence tolerance tightened from 1e-5 to 1e-6 by default.
*   Reference replication script (`replicate_jasa.R`, validating
    against Blackburne & Frank 2007) added to the test suite.

# ardlverse 1.1.3

*   Initial CRAN release of comprehensive ARDL framework (Panel,
    Bootstrap, Fourier, Quantile, Augmented, NARDL, Rolling/Recursive).
