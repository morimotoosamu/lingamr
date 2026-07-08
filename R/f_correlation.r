# =============================================================================
# F-correlation (kernel canonical correlation, Bach & Jordan 2002)
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/utils/_f_correlation.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
#
# NOTE: this is a standalone port and intentionally does not reuse
# incomplete_cholesky_gauss() from search_causal_order.r. The two have
# different stopping rules (see incomplete_cholesky_fcorr() below), and
# reusing the DirectLiNGAM version would change the resulting rank and
# therefore the numerical value of f_correlation().
# =============================================================================


#' Pivoted incomplete Cholesky decomposition of a Gaussian kernel matrix,
#' as used by [f_correlation()]
#'
#' Same greedy pivoted algorithm as `incomplete_cholesky_gauss()`
#' (search_causal_order.r), but with the stopping rule used by the upstream
#' `_f_correlation.py`: continue while the *sum* of the remaining diagonal
#' residuals exceeds `tol` (no rank cap), rather than stopping once the
#' single largest residual drops below a fixed tolerance.
#'
#' @param x input vector (length n)
#' @param sigma width of the Gaussian kernel
#' @param tol stop once the sum of the remaining diagonal residuals is at or below this
#' @return n x d matrix `G` with `tcrossprod(G) ~= K`
#' @keywords internal
incomplete_cholesky_fcorr <- function(x, sigma, tol) {
  n <- length(x)
  d_vec <- rep(1, n) # diag(K) == 1 for the Gaussian kernel
  perm <- seq_len(n)
  G <- matrix(0, n, n)
  rank <- 0L

  for (k in seq_len(n)) {
    remaining <- k:n
    if (sum(d_vec[perm[remaining]]) <= tol) break

    pivot_rel <- which.max(d_vec[perm[remaining]])
    pivot_pos <- remaining[pivot_rel]

    if (pivot_pos != k) {
      perm[c(k, pivot_pos)] <- perm[c(pivot_pos, k)]
    }
    pk <- perm[k]
    if (d_vec[pk] <= .Machine$double.eps) break # remaining residual is numerically zero; further columns add no information
    G[pk, k] <- sqrt(d_vec[pk])
    rank <- k

    if (k < n) {
      rest <- perm[(k + 1):n]
      Kcol <- exp(-1 / (2 * sigma^2) * (x[rest] - x[pk])^2)
      if (k > 1) {
        Kcol <- Kcol - as.vector(G[rest, seq_len(k - 1), drop = FALSE] %*% G[pk, seq_len(k - 1)])
      }
      G[rest, k] <- Kcol / G[pk, k]
      d_vec[rest] <- pmax(0, d_vec[rest] - G[rest, k]^2)
    }
  }

  G[, seq_len(rank), drop = FALSE]
}


#' Low-rank SVD transform used inside [f_correlation()]
#'
#' Given the (already column-centered) incomplete-Cholesky factor `G`,
#' returns the orthonormalized basis `U` and the shrinkage vector `R`
#' used to assemble the block canonical-correlation matrix `R_kappa`.
#'
#' @param G n x d centered incomplete-Cholesky factor
#' @param kappa regularization parameter
#' @param n sample size
#' @return list(U = n x d' matrix, R = length d' vector); `R` has length 0
#'   if no eigenvalue of `crossprod(G)` clears the `kappa` threshold (i.e.
#'   `G` carries no retainable rank), which callers must handle explicitly.
#' @keywords internal
fcorr_svd_transform <- function(G, kappa, n) {
  A_mat <- crossprod(G) # t(G) %*% G
  eig <- eigen(A_mat, symmetric = TRUE)
  D <- eig$values
  V <- eig$vectors

  keep <- D >= n * kappa * 1e-2
  D <- D[keep]
  V <- V[, keep, drop = FALSE]

  if (length(D) == 0) {
    return(list(U = matrix(0, nrow = nrow(G), ncol = 0), R = numeric(0)))
  }

  ord <- order(D, decreasing = TRUE)
  D <- D[ord]
  V <- V[, ord, drop = FALSE]

  U <- G %*% V %*% diag(1 / sqrt(D), nrow = length(D))
  R <- D / (n * kappa / 2 + D)

  list(U = U, R = R)
}


#' F-correlation (kernel canonical correlation) between two variables
#'
#' Bach & Jordan (2002) kernel canonical correlation, as used by
#' BottomUpParceLiNGAM's `independence = "fcorr"` option. Returns a value in
#' (roughly) `[0, 1]`; larger means more dependent.
#'
#' A constant `x` or `y` (zero variance) carries no dependence information
#' and is treated as trivially independent (returns `0`) rather than
#' propagating a division-by-zero standardization into NaN/Inf.
#'
#' @param x numeric vector
#' @param y numeric vector (same length as x)
#' @return F-correlation value (scalar)
#' @keywords internal
f_correlation <- function(x, y) {
  if (length(x) != length(y)) {
    stop("f_correlation(): x and y must have the same length.", call. = FALSE)
  }
  if (anyNA(x) || anyNA(y)) {
    stop("f_correlation(): x and y must not contain NA/NaN values.", call. = FALSE)
  }
  n <- length(x)
  if (sd_pop(x) == 0 || sd_pop(y) == 0) {
    return(0)
  }
  x <- (x - mean(x)) / sd_pop(x)
  y <- (y - mean(y)) / sd_pop(y)

  if (n > 1000) {
    kappa <- 2e-3
    sigma <- 0.5
  } else {
    kappa <- 2e-2
    sigma <- 1.0
  }
  tol <- n * kappa * 1e-2

  G1 <- incomplete_cholesky_fcorr(x, sigma, tol)
  G2 <- incomplete_cholesky_fcorr(y, sigma, tol)
  G1 <- sweep(G1, 2, colMeans(G1))
  G2 <- sweep(G2, 2, colMeans(G2))

  t1 <- fcorr_svd_transform(G1, kappa, n)
  t2 <- fcorr_svd_transform(G2, kappa, n)
  d1 <- length(t1$R)
  d2 <- length(t2$R)
  if (d1 == 0 || d2 == 0) {
    return(0) # no retainable rank on either side => no canonical correlation structure
  }

  cross <- crossprod(t1$U, t2$U) # d1 x d2
  block12 <- cross * t1$R * rep(t2$R, each = d1)

  R_kappa <- rbind(
    cbind(diag(d1), block12),
    cbind(t(block12), diag(d2))
  )
  R_kappa <- (R_kappa + t(R_kappa)) / 2

  lambda_min <- min(eigen(R_kappa, symmetric = TRUE, only.values = TRUE)$values)
  1 - lambda_min
}
