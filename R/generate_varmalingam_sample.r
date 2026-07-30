# =============================================================================
# Sample data generation for VARMA-LiNGAM
# =============================================================================


#' Generate sample data from a VARMA-LiNGAM model
#'
#' Generates a 3-variable time series following a VARMA(1,1)-LiNGAM model with
#' a strictly acyclic instantaneous structure B0, a lag-1 autoregressive matrix
#' Phi1, a lag-1 moving-average matrix Theta1, and non-Gaussian (uniform)
#' errors. The reduced-form recursion is
#' `x(t) = Phi1 x(t-1) + n(t) + Theta1 n(t-1)` with `n(t) = (I - B0)^{-1} e(t)`.
#'
#' @param n number of time points to return (after burn-in)
#' @param seed random seed (NULL allowed)
#' @return list with `data` (data frame, n x 3), the reduced-form true
#'   matrices `true_B0`, `true_phi1`, `true_theta1`, and the structural-form
#'   counterparts `true_psi1 = (I - B0) Phi1` and
#'   `true_omega1 = (I - B0) Theta1 (I - B0)^{-1}` for comparison against
#'   [lingam_varma()] estimates
#' @export
#' @examples
#' sample <- generate_varmalingam_sample(n = 500, seed = 1)
#' head(sample$data)
generate_varmalingam_sample <- function(n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- 3L

  # instantaneous structure: x0 -> x1 -> x2 (strictly acyclic, lower-triangular)
  B0 <- matrix(0, p, p)
  B0[2, 1] <- 0.6
  B0[3, 2] <- -0.5

  # lag-1 AR coefficients (spectral radius well below 1 => stationary)
  Phi1 <- diag(c(0.3, 0.25, 0.3))
  Phi1[1, 3] <- 0.2

  # lag-1 MA coefficients (eigenvalue moduli below 1 => invertible)
  Theta1 <- diag(c(0.25, 0.2, 0.3))
  Theta1[2, 1] <- 0.2

  burn_in <- 100L
  total <- n + burn_in
  # non-Gaussian errors (uniform)
  e <- matrix(stats::runif(total * p, min = -1, max = 1), nrow = total, ncol = p)

  ib0_inv <- solve(diag(p) - B0)
  # reduced-form disturbances n(t) = (I - B0)^{-1} e(t)
  N <- e %*% t(ib0_inv)

  X <- matrix(0, nrow = total, ncol = p)
  X[1, ] <- N[1, ]
  for (t in 2:total) {
    X[t, ] <- Phi1 %*% X[t - 1L, ] + N[t, ] + Theta1 %*% N[t - 1L, ]
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- paste0("x", seq_len(p) - 1L)

  list(
    data = as.data.frame(X),
    true_B0 = B0,
    true_phi1 = Phi1,
    true_theta1 = Theta1,
    true_psi1 = (diag(p) - B0) %*% Phi1,
    true_omega1 = (diag(p) - B0) %*% Theta1 %*% ib0_inv
  )
}
