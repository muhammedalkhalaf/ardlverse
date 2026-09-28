# Critical values and approximate p-values for the Pesaran, Shin and Smith
# (2001) bounds test from the response surface regressions of Kripfganz and
# Schneider (2020). Port of 'ardlbounds' (Stata package 'ardl', version
# 1.0.6, by S. Kripfganz and D. C. Schneider); the coefficients in
# R/sysdata.rda ('ks_coefs') are the file 'ardl_surfreg_coefs.dta' of that
# package.
#
# stat      "F" or "t"
# case      PSS case 1 to 5 (the t statistic is tabulated for cases 1, 3, 5;
#           cases 2 and 4 use cases 3 and 5)
# k         number of long-run forcing variables
# n         number of observations; NULL gives asymptotic values
# sr        number of short-run coefficients (regressors other than the
#           deterministic terms and the k + 1 level terms)
# siglevels significance levels in percent
# value     observed statistic (optional) for approximate p-values
#
# Returns list(cv = 2 x length(siglevels) matrix (rows I0, I1),
#              pvalue = c(I0, I1) or NULL)
.ks_bounds <- function(stat = c("F", "t"), case, k, n = NULL, sr = 0,
                       siglevels = c(10, 5, 1), value = NULL) {
  stat <- match.arg(stat)
  case <- as.integer(case)
  if (!case %in% 1:5) stop("'case' must be 1 to 5.", call. = FALSE)
  case_t <- if (stat == "t" && case %in% c(2L, 4L)) case + 1L else case
  asym <- is.null(n) || is.na(n)
  df1 <- 1 + k + (case %in% c(2L, 4L))
  df2 <- if (asym) NA_real_ else n - (1 + k + sr + (case >= 2L) + (case >= 4L))
  if (!asym && df2 < n / 2) {
    na <- matrix(NA_real_, 2L, length(siglevels),
                 dimnames = list(c("I0", "I1"), paste0(siglevels, "%")))
    return(list(cv = na, pvalue = if (!is.null(value)) c(I0 = NA_real_, I1 = NA_real_)))
  }

  cvs <- matrix(NA_real_, 2L, length(siglevels),
                dimnames = list(c("I0", "I1"), paste0(siglevels, "%")))
  pv <- c(I0 = NA_real_, I1 = NA_real_)
  for (I in 0:1) {
    d <- ks_coefs[ks_coefs$stat == stat & ks_coefs$c == case_t & ks_coefs$I == I, ]
    d <- d[order(d$p), ]
    cv <- 0
    if (!asym)
      cv <- cv + (d$theta_0_2_0 + d$theta_0_2_1 * sr) / n^2 +
        (d$theta_0_3_0 + d$theta_0_3_1 * sr) / n^3
    for (j in 0:4) {
      cv <- cv + d[[sprintf("theta_%d_0_0", j)]] / (k + 1)^j
      if (!asym)
        cv <- cv + d[[sprintf("theta_%d_1_0", j)]] / ((k + 1)^j * n) +
          d[[sprintf("theta_%d_1_1", j)]] * sr / ((k + 1)^j * n)
    }
    for (s in seq_along(siglevels)) {
      hit <- abs(d$p - siglevels[s] * 100) < 1e-8
      if (any(hit)) cvs[I + 1L, s] <- mean(cv[hit])
    }
    if (!is.null(value) && is.finite(value)) {
      if (value > max(cv)) {
        pv[I + 1L] <- if (stat == "t") 1 else 0
      } else if (value < min(cv)) {
        pv[I + 1L] <- if (stat == "t") 0 else 1
      } else {
        ad <- abs(cv - value)
        mp <- mean(which(ad == min(ad)))
        use <- seq_along(cv) >= mp - 4 & seq_along(cv) <= mp + 4
        pr <- d$p[use] / 10000
        inv <- if (!asym) {
          if (stat == "F") stats::qf(pr, df1, df2, lower.tail = FALSE)
          else -stats::qt(pr, df2, lower.tail = FALSE)
        } else {
          if (stat == "F") stats::qchisq(pr, df1, lower.tail = FALSE)
          else stats::qnorm(pr)
        }
        z <- cv[use]
        b <- stats::lm.fit(cbind(1, z, z^2), inv)$coefficients
        crit <- b[1] + b[2] * value + b[3] * value^2
        pv[I + 1L] <- if (!asym) {
          if (stat == "F") stats::pf(crit, df1, df2, lower.tail = FALSE)
          else stats::pt(crit, df2)
        } else {
          if (stat == "F") stats::pchisq(crit, df1, lower.tail = FALSE)
          else stats::pnorm(crit)
        }
      }
    }
  }
  list(cv = cvs, pvalue = if (!is.null(value)) pv)
}
