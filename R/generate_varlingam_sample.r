# =============================================================================
# Sample data generation for VAR-LiNGAM
# =============================================================================


#' Generate sample data from a VAR-LiNGAM model
#'
#' Generates a 3-variable time series following a VAR-LiNGAM model with a
#' strictly acyclic instantaneous structure B0, a lag-1 coefficient matrix M1,
#' and non-Gaussian (uniform) errors.
#'
#' @param n number of time points to return (after burn-in)
#' @param seed random seed (NULL allowed)
#' @return list with `data` (data frame, n x 3), `true_B0` (instantaneous
#'   matrix), and `true_M1` (lag-1 coefficient matrix)
#' @export
#' @examples
#' sample <- generate_varlingam_sample(n = 500, seed = 1)
#' head(sample$data)
generate_varlingam_sample <- function(n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- 3L

  # instantaneous structure: x0 -> x1 -> x2 (strictly acyclic, lower-triangular)
  B0 <- matrix(0, p, p)
  B0[2, 1] <- 0.6
  B0[3, 2] <- -0.5

  # lag-1 coefficients
  M1 <- matrix(0, p, p)
  diag(M1) <- c(0.4, 0.3, 0.5)
  M1[1, 3] <- 0.3

  burn_in <- 100L
  total <- n + burn_in
  # non-Gaussian errors (uniform)
  e <- matrix(stats::runif(total * p, min = -1, max = 1), nrow = total, ncol = p)

  X <- matrix(0, nrow = total, ncol = p)
  ib0_inv <- solve(diag(p) - B0)
  for (t in 2:total) {
    X[t, ] <- ib0_inv %*% (M1 %*% X[t - 1L, ] + e[t, ])
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- paste0("x", seq_len(p) - 1L)

  list(
    data = as.data.frame(X),
    true_B0 = B0,
    true_M1 = M1
  )
}
