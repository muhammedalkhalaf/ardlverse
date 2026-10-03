# ardlverse 2.1.0

Monte Carlo rates below: y and x independent random walks, n = 100, 200 samples, B = 199, nominal level 5%; Monte Carlo standard errors in parentheses.

## Breaking changes

* `boot_ardl()`: the bootstrap is now recursive and its default follows Bertelli, Vacca and Zoia (2022, Section 3): one restricted model per statistic (Fov, t, Find), y* and x* generated observation by observation (x* from the marginal model of the differenced regressors with p - 1 lags, BVZ eq. 19), residuals recentred after each draw and a random block of initial values. `nulls = "joint"` uses the null of the Fov test for all statistics (McNown, Sam and Goh, 2018, Steps 1-8) applied to the conditional ECM, a package choice (MSG write the y equation without the contemporaneous differences). The 2.0.3 bootstrap regenerated only the dependent variable and kept the lagged levels and differences of the data fixed, with one null for F and t; it rejected in 29.0% (3.2%) (F) and 53.0% (3.5%) (t) of the samples. The new bootstrap rejects in 4.5% (1.5%) (Fov), 7.5% (1.9%) (t) and 10.0% (2.1%) (Find); with `nulls = "joint"` 7.0% (1.8%), 8.0% (1.9%) and 8.5% (2.0%). All bootstrap p-values and critical values change.
* `boot_ardl()`: the conclusion was "cointegration" when F OR t was significant. Cointegration is now concluded only when Fov, t and Find all reject; Fov rejecting without t is reported as the degenerate case of the first type, Fov and t rejecting without Find as the degenerate case of the second type. The new field `decision` holds the code (`"COINTEGRATION"`, `"NO_COINTEGRATION"`, `"DEGENERATE_1"`, `"DEGENERATE_2"`, `"UNDETERMINED"`). The rate of a "cointegration" conclusion under the null fell from 57.5% (3.5%) to 3.0% (1.2%).
* `boot_ardl()`: new statistic `Find_stat` (F test on the lagged regressors) with `boot_Find`, `cv_Find`, `p_value_Find`. Critical values are order statistics of the bootstrap distribution (BVZ eqs. 24-25) at 10%, 5%, 2.5% and 1% (names `"90%"`, `"95%"`, `"97.5%"`, `"99%"` for the F tests, `"10%"`, `"5%"`, `"2.5%"`, `"1%"` for t) instead of `quantile()`. `parallel` and `ncores` are deprecated and ignored.
* `boot_ardl(p = 0)` now gives an error. In 2.0.3 it ran and gave the same model and F statistic as `p = 1` (an ARDL(0, q) has no lagged level of y, so the bounds test is not defined); `p` must be at least 1.
* `aardl()`, types `"bootstrap"`, `"bnardl"`, `"fbootstrap"` and `"fbnardl"`: same recursive bootstrap (new argument `nulls`); the partial sums are rebuilt from x* with the same decomposition as the data and the Fourier terms at `fourier_k` are held fixed. The 2.0.3 bootstrap kept the design fixed, used one null and did not recentre the residuals; its rejection rates were 32.0% (3.3%) to 86.0% (2.5%) (F), 50.5% (3.5%) to 97.5% (1.1%) (t) and 18.0% (2.7%) to 49.5% (3.5%) (F_ind). New rates: F 4.0% (1.4%) to 9.0% (2.0%), t 6.5% (1.7%) to 12.0% (2.3%), F_ind 6.5% (1.7%) to 16.0% (2.6%), and 2.0% (1.0%) to 6.0% (1.7%) for a "cointegration" decision (14.5% (2.5%) to 47.0% (3.5%) in 2.0.3). F_ind over-rejects with Fourier terms and partial sums (type `"fbnardl"`: 16.0% (2.6%); an independent review Monte Carlo gave 11.3% (1.8%)). The decision codes are those of `boot_ardl()` (the former `"INCONCLUSIVE"` for F and t significant without F_ind is now `"DEGENERATE_2"`), and the bootstrap distributions keep failed replications as `NA`.
* `aardl()`, types `"fourier"` and `"fnardl"`: the bounds decision with Fourier terms is removed (no valid bounds exist for the bounds test with Fourier terms; under the null its rejection rate was 22.5% (3.0%) and 28.5% (3.2%)). The decision is now `"NOT_AVAILABLE"`, and the bounds shown are labelled as bounds for the model without Fourier terms, not valid with Fourier terms, consistent with `fourier_bounds_test()`; the short-run count of these bounds no longer includes the Fourier terms.
* `mtnardl()`: the dynamic multipliers were wrong (the level coefficient was used as a short-run coefficient and the lagged dynamics were ignored; for example regime 2 at horizon 0 gave 0.090 instead of 0.303). They are now computed by the full recursion of the estimated ECM, equal to the simulated response of the level of y to a permanent unit increase in each partial sum, for horizons 0 to `horizon` (new argument, default 30; `multipliers$cumulative` now has `horizon + 1` rows and `multipliers$horizons` starts at 0).
* `mtnardl()`: coefficient names changed. The coefficients were named after the design matrix (`designy_lag`, `designd.<y>.l1`, and the differenced partial sums had no names); they are now `y_lag`, the regime names, `dy_l1`, ... and `d_<regime>_l<i>` for the differenced partial sums, with the same names in `vcov()`.
* `mtnardl()`: the bounds used k = number of original regressors. The default is now k = number of level regressors (partial sums plus linear regressors), with `bounds_k = "original"` for the former choice; this is a package choice pending verification against Pal and Mitra (2016) and Shin, Yu and Greenwood-Nimmo (2014). The asymptotic decision now requires both the F and the t statistic beyond their I(1) bounds (it used F only).
* `mtnardl(bootstrap = TRUE)`: the recursive bootstrap replaces the fixed-design bootstrap (rejection rates 46.0% (3.5%) for F and 73.0% (3.1%) for t in 2.0.3; now 3.5% (1.3%), 7.0% (1.8%) and 9.5% (2.1%) for Fov, t and Find), the Find statistic is added (`Find_stat`) and the decision requires Fov, t and Find.
* `mtnardl(auto_select = TRUE, n_thresholds = 1)` always returned 0 because the internal call failed silently; it now returns the threshold that minimises the AIC over 0 and the deciles of the changes.
* `mtnardl()`: `asymmetry_tests[[v]]` gains an element `joint` (Wald test of equal level coefficients across all regimes, chi-squared with R - 1 degrees of freedom) and every test a field `df`.
* `boot_ardl()`, `aardl()` and `mtnardl()` used `all.vars()` on the formula, which silently turned `log(y) ~ log(x)` into `y ~ x`; transformed terms and interactions now give an informative error. Missing values give an error instead of being dropped (dropping rows closed gaps in the series).
* `boot_ardl()`, `aardl()` and `mtnardl()` called `set.seed()` on the session's generator; the seed is now used locally and the caller's random number state is restored.

## New features

* New exported function `fbnardl()`: Fourier (bootstrap) nonlinear ARDL with positive and negative partial sums, selection of the Fourier frequency (integer frequencies chosen by minimum SSR as in Enders and Lee (2012, Oxford Bulletin of Economics and Statistics) by default, or a fractional grid) and of the lags by AIC or BIC on a common sample, Kripfganz and Schneider (2020) bounds without Fourier terms, recursive bootstrap with the frequency and lags selected again on every bootstrap sample (package choice; BVZ step 5(c) re-estimates the unrestricted model but does not discuss re-selection; or fixed with `kstar` and `lags`), delta-method long-run coefficients and symmetry test, a short-run symmetry test on the sum over lags, cumulative dynamic multipliers by full recursion with parametric bands, and diagnostics. With `bootstrap = "bvz"` the default x model is the marginal model of BVZ eqs. 19-20 (`xdgp = "vecm"`); `xdgp = "rw"` regresses the differences on the intercept, Fourier terms and `exog` only (a package choice). `bootstrap = "mcnown"` uses the Fov null for all statistics applied to the conditional ECM (package choice). The R code is a port written for this package, following the Stata module fbnardl; it corrects the long-run Wald gradient (13.36 instead of 631.90 in the reference example), the fixed-design bootstrap, the hard-coded asymptotic case III bounds, the use of bounds with Fourier terms, the lag selection on different samples, the multipliers without lagged short-run terms (0.928 instead of 0.986 at horizon 1 in the reference example), and the storage of failed replications as 0 found in the stand-alone version 1.0.2. Rejection rates with the frequency and lags selected again in each bootstrap sample (package choice; BVZ step 5(c) re-estimates the unrestricted model but does not discuss re-selection): 7.0% (1.8%) (Fov), 11.5% (2.3%) (t), 11.5% (2.3%) (Find), 3.5% (1.3%) for a "cointegration" decision (stand-alone 1.0.2: 98.5% (0.9%), 99.0% (0.7%), 67.5% (3.3%)).
* `mtnardl()`: long-run coefficients with delta-method standard errors computed with the full covariance matrix (`long_run_table`), delta-method differences of the long-run coefficients between regimes (`lr_differences`), the Find statistic, `decompose` (a subset of regressors to decompose; the others enter linearly) and `threshold_type = "quantile"` (thresholds given as probabilities and converted to quantiles of the changes of each variable), which give a migration path from `fqardl::mtnardl()`.
* New internal file `R/ardl_boot_engine.R` (engine version 1.1.0): the recursive bootstrap engine shared by all bootstrap tests, with a reproduction check (the recursion with the original residuals reproduces the data), a design check (the bootstrap model equals the estimated model), failed replications stored as `NA` with a warning and counts, order-statistic critical values and the combined decision. With a random block of initial values the partial sums continue the data's levels in the block, so the generating model is the estimated restricted model. The x equation has p lags of the differences when the y equation has p lagged differences of y (BVZ eq. 19). The argument `j0 = 1` gives the unconditional y equation of MSG eq. 12 (no contemporaneous differences). Recentring once (MSG eq. 13) subtracts the residual mean: eq. 13 divides the sum of the residuals by n - q - 1, the number of residuals of their ARDL(1,1) model, and applies no rescaling factor.

## Other changes

* DESCRIPTION: "bootstrap-based bounds testing per Pesaran, Shin and Smith (2001)" replaced by a description of the bounds test with Kripfganz and Schneider (2020) critical values and recursive bootstrap critical values (Bertelli, Vacca and Zoia, 2022; McNown, Sam and Goh, 2018).
* `boot_ardl()` documentation: it no longer claims to follow McNown, Sam and Goh (2018) for a scheme that differed from it, documents that `F_overall` equals `F_stat`, and documents the decision rule. `aardl()` documentation no longer attributes its bootstrap to McNown, Sam and Goh (2018) and Sam, McNown and Goh (2019).
* The start-up message showed version 1.0.0; it now shows the installed version.
* Examples run with small numbers of bootstrap replications.
* Tests: statistics against independent `lm()`/`anova()` computations, bootstrap draws against an independent recursion (BVZ with partial sums and a random block; MSG joint null), and size Monte Carlo tests (skipped on CRAN) for `boot_ardl()`, `aardl()`, `mtnardl()` and `fbnardl()`.

## Pending verification

The following papers were not available for this release; the choices below are documented as package choices:

* Pal and Mitra (2016) and Shin, Yu and Greenwood-Nimmo (2014): k in the bounds of `mtnardl()`, regime formation, multiplier definition and sign convention for negative partial sums, the short-run symmetry test on the sum over lags in `fbnardl()`.
* Yilanci, Bozoklu and Gorus (2020) and Omay (2015): the fractional frequency grid of `fbnardl()` and its validity; the default is the integer grid.
* Sam, McNown and Goh (2019): bootstrap details of the augmented ARDL test in `aardl()`.

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
audited `panel_ardl()` against Stata's `xtpmg` (Blackburne and Frank 2007)
and contributed corrections that bring the implementation into strict
alignment with the original Pesaran, Shin and Smith (1999) framework.

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
    variance-covariance formula used by Blackburne and Frank (2007).

7.  **Sub-optimal PMG initialization.** The previous code ran the full
    MG estimator to generate PMG starting values. PMG is now initialized
    from a simple pooled OLS of the lagged dependent variable on the
    levels of X; faster, avoids convergence risk if MG fails, and
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
    against Blackburne and Frank 2007) added to the test suite.

# ardlverse 1.1.3

*   Initial CRAN release of comprehensive ARDL framework (Panel,
    Bootstrap, Fourier, Quantile, Augmented, NARDL, Rolling/Recursive).
