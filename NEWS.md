# ardlverse 2.1.0

## Model Specification & Methodological Improvements

Critical audit of the entire package led to significant fixes in model design and lag management:

1.  **Improved Model Output Readability**: Internal variable naming has been overhauled across all estimators (`aardl`, `mtnardl`, `pnardl`, etc.). Design matrix prefixes like `design` have been removed in favor of explicit names (e.g., `d.yvar`, `yvar.l1`, `d.xvar.l{lag}`).
2.  **Corrected Lag Alignment**: Fixed a major specification error where contemporaneous differences of the dependent variable ($\Delta y_t$) were accidentally included as regressors in multiple models (`aardl`, `mtnardl`, `rolling_ardl`). This previously caused spurious $R^2 = 1$ in some cases.
3.  **Unified Lag Definitions**: Synchronized the interpretation of $p$ (dependent lags) and $q$ (independent lags) across all model types (Linear, Nonlinear, Panel, Fourier, and Quantile).
4.  **Bootstrap ARDL Robustness**: Improved the internal data preparation for bootstrap procedures, ensuring that the restricted and unrestricted models are perfectly aligned in terms of sample size and lag structure.

## Enhanced Diagnostics

1.  **Corrected CUSUM/CUSUMSQ Bounds**: Re-implemented the critical bounds for stability tests following Brown, Durbin, and Evans (1975). The bounds are now correctly non-linear and fan out as expected.
2.  **Panel-Aware Serial Correlation**: The Breusch-Godfrey LM test (`ardl_diagnostics`) now supports panel data. It correctly computes lagged residuals within groups (avoiding "cross-contamination" at unit boundaries).
3.  **Augmented ARDL (AARDL) Diagnostics**: `aardl` models now feature a full diagnostic suite including Ramsey RESET, CUSUM/CUSUMSQ stability tests, and automated residuals visualization via `plot()`.
4.  **Improved Heteroskedasticity Tests**: "Manual" implementations of Breusch-Pagan and ARCH tests now correctly include all model regressors, which is essential for valid inference in dynamic ARDL models.
5.  **MT-NARDL Support**: `ardl_diagnostics` now fully supports Multiple-Threshold NARDL models.

## Panel Data Updates

1.  **Metadata Export**: All panel estimators (`pmg`, `mg`, `dfe`) for both linear and nonlinear specifications now export residuals, fitted values, and group identifiers. This enables the use of the standard diagnostic suite on panel models.
2.  **Harmonized Output**: Panel model objects now have a more consistent structure, facilitating easier post-estimation analysis.

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
