# Tools shared by mtnardl() and fbnardl().

# Cumulative dynamic multiplier of y to a permanent unit increase at h = 0 in
# one level regressor w of the ECM
#   dy_t = rho y_{t-1} + theta w_{t-1} + sum_i phi_i dy_{t-i}
#          + sum_{l=0}^{L} pi_l dw_{t-l} + ...
# m_h is the response of the level of y at horizon h (h = 0..H): with
# m_{-1} = 0, dm_h = rho m_{h-1} + theta 1(h >= 1) + sum_i phi_i dm_{h-i}
# + pi_h 1(h <= L), m_h = m_{h-1} + dm_h. This equals the response obtained
# by simulating the estimated ECM; m_h tends to -theta / rho.
.ecm_multiplier <- function(rho, theta, phi = numeric(0), pi = numeric(0), H = 20) {
  m <- dm <- numeric(H + 1)
  for (h in 0:H) {
    s <- if (h >= 1) rho * m[h] + theta else 0
    for (i in seq_along(phi)) if (h - i >= 0) s <- s + phi[i] * dm[h - i + 1]
    if (h < length(pi)) s <- s + pi[h + 1]
    dm[h + 1] <- s
    m[h + 1] <- (if (h >= 1) m[h] else 0) + s
  }
  m
}

# Long-run coefficients -theta_j / rho and their delta-method covariance with
# the full coefficient covariance matrix V. 'irho' and 'itheta' index b.
.lr_delta <- function(b, V, irho, itheta) {
  rho <- b[irho]
  th <- b[itheta]
  lr <- -th / rho
  G <- matrix(0, length(itheta), length(b))
  for (j in seq_along(itheta)) {
    G[j, irho] <- th[j] / rho^2
    G[j, itheta[j]] <- -1 / rho
  }
  VL <- G %*% V %*% t(G)
  list(lr = unname(lr), V = VL, se = sqrt(pmax(diag(VL), 0)), G = G)
}

# Wald statistic for R b = r
.wald_lin <- function(b, V, R, r = rep(0, nrow(R))) {
  d <- R %*% b - r
  drop(t(d) %*% solve(R %*% V %*% t(R), d))
}
