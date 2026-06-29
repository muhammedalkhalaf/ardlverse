# Replication of Blackburne & Frank (2007, Stata Journal)
# Data: Pesaran, Shin & Smith (1999) OECD consumption panel
# Reference output from st0125.md
#
# Target PMG results:
#   LR:  pi = -0.4658438 (se 0.0567332), y = 0.9044336 (se 0.0086815)
#   SR:  ec = -0.1998761, D.pi = -0.0182588, D.y = 0.3269355, _cons = 0.1544507
#   Log-likelihood = 2327.075

source("panel_ardl.R")

# ── Load data ─────────────────────────────────────────────────────────────────
# Requires haven: install.packages("haven")
library(haven)
jasa2 <- read_dta("http://fmwww.bc.edu/repec/bocode/j/jasa2.dta")
# jasa2 <- read_dta("stata/jasa2.dta")

# Pass FULL dataset; use start_time=1962 to mirror Stata's "if year>=1962"
# STATA uses pre-sample data (1960-1961) to compute lags for the 1962 first obs.
df <- jasa2

cat("Data dimensions:", nrow(df), "rows,", length(unique(df$id)), "groups\n")
cat("Full years:", min(df$year), "to", max(df$year), "\n\n")

# ── PMG: ARDL(1,1,1)  ─────────────────────────────────────────────────────────
# Stata: xtpmg d.c d.pi d.y if year>=1962, lr(l.c pi y) ec(ec) replace pmg
# formula = c ~ pi + y   (levels; the function creates diffs/lags internally)
# p=1 (one lag of c), q=c(1,1) (one lag per x-var → ARDL(1,1,1))

cat("=== PMG ===\n")
pmg <- panel_ardl(
  formula    = c ~ pi + y,
  data       = df,
  id         = "id",
  time       = "year",
  p          = 1,
  q          = c(1, 1),
  estimator  = "pmg",
  start_time = 1962,
  maxiter    = 200,
  tol        = 1e-8
)
summary(pmg)

cat("\nLR coefficients (target: pi=-0.4658438, y=0.9044336):\n")
print(pmg$long_run)
cat("\nLR std errors  (target: pi=0.0567332,  y=0.0086815):\n")
print(pmg$long_run_se)
cat("\nSR coefficients (target: ec=-0.1998761, D.pi=-0.0182588, D.y=0.3269355, _cons=0.1544507):\n")
print(pmg$short_run)
cat("\nLog-likelihood (target: 2327.075):", pmg$loglik, "\n\n")

# ── MG ────────────────────────────────────────────────────────────────────────
# Stata: xtpmg d.c d.pi d.y if year>=1962, lr(l.c pi y) ec(ec) replace mg

cat("=== MG ===\n")
mg <- panel_ardl(
  formula    = c ~ pi + y,
  data       = df,
  id         = "id",
  time       = "year",
  p          = 1,
  q          = c(1, 1),
  estimator  = "mg",
  start_time = 1962
)
summary(mg)

# ── DFE ───────────────────────────────────────────────────────────────────────
# Stata: xtpmg d.c d.pi d.y if year>=1962, lr(l.c pi y) ec(ec) replace dfe cluster(id)

cat("=== DFE ===\n")
dfe <- panel_ardl(
  formula    = c ~ pi + y,
  data       = df,
  id         = "id",
  time       = "year",
  p          = 1,
  q          = c(1, 1),
  estimator  = "dfe",
  start_time = 1962,
  cluster    = TRUE      # cluster(id) robust SEs, as in the paper
)
summary(dfe)

cat("\nDFE targets (LR): pi=-0.266343 (SE=0.102506), y=0.9120574 (SE=0.0468008)\n")
cat("DFE targets (SR): ec=-0.1794146 (SE=0.0434584), D.pi=-0.0280826 (SE=0.0325622),",
    "D.y=0.3811944 (SE=0.070876), _cons=0.1257634 (SE=0.0805454)\n\n")

# ── Hausman tests ─────────────────────────────────────────────────────────────
# Stata: hausman mg pmg, sigmamore  → chi2(2)=1.06, p=0.5887
# Stata: hausman mg DFE, sigmamore  → chi2(2)=0.00, p=1.0000
cat("=== Hausman test (MG vs PMG) — target: chi2(2)=1.06, p=0.5887 ===\n")
hausman_test(mg, pmg)

cat("=== Hausman test (MG vs DFE) — target: chi2(2)=0.00, p=1.0000 ===\n")
hausman_test(mg, dfe)
