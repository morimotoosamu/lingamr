# =============================================================================
# HSIC (Hilbert-Schmidt Independence Criterion) gamma approximation test
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/hsic.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
#
# Shared by lingam_parce.r (BottomUpParceLiNGAM) and potentially future
# HSIC-based ports (e.g. RCD).
# =============================================================================


#' Median-heuristic kernel width for HSIC
#'
#' Uses only the first 100 points (not a random subsample) to keep the
#' O(n^2) pairwise-distance computation cheap, matching the upstream
#' implementation exactly.
#'
#' @param x numeric vector
#' @return kernel width (scalar)
#' @keywords internal
hsic_kernel_width <- function(x) {
  n <- length(x)
  m <- min(n, 100L)
  xs <- x[seq_len(m)]
  D <- outer(xs, xs, function(a, b) (a - b)^2)
  d <- D[upper.tri(D)]
  d <- d[d > 0]
  if (length(d) == 0) {
    return(1.0)
  }
  sqrt(0.5 * stats::median(d))
}


#' Gaussian Gram matrix and its double-centered version
#'
#' @param x numeric vector
#' @param width kernel width from [hsic_kernel_width()]
#' @return list(K = Gram matrix, Kc = centered Gram matrix)
#' @keywords internal
hsic_gram_matrix <- function(x, width) {
  n <- length(x)
  H <- outer(x, x, function(a, b) (a - b)^2)
  K <- exp(-H / (2 * width^2))
  row_sums <- rowSums(K)
  col_sums <- colSums(K)
  total <- sum(K)
  Kc <- K - row_sums / n - rep(col_sums, each = n) / n + total / n^2
  list(K = K, Kc = Kc)
}


#' HSIC independence test with gamma approximation
#'
#' Faithful port of `hsic_test_gamma()` (hsic.py). O(n^2) in the sample size
#' because it forms the full n x n Gram matrices; not recommended for n in
#' the thousands.
#'
#' The gamma-approximation variance estimator is only defined for n >= 6
#' (its closed form divides by `(n-1)(n-2)(n-3)`); below that this errors
#' instead of silently returning a NaN p-value that would otherwise
#' propagate into `NA`-valued rejection decisions in callers. A constant
#' input (zero variance) carries no dependence information, so it is
#' treated as trivially independent (`p = 1`) rather than routed through
#' the degenerate kernel-width / centered-Gram-matrix computation.
#'
#' @param X numeric vector
#' @param Y numeric vector (same length as X)
#' @return list(stat = HSIC test statistic, p = gamma-approximated p-value)
#' @keywords internal
hsic_test_gamma <- function(X, Y) {
  if (length(X) != length(Y)) {
    stop("hsic_test_gamma(): X and Y must have the same length.", call. = FALSE)
  }
  n <- length(X)
  if (n < 6) {
    stop(
      "hsic_test_gamma(): the gamma-approximation test requires at least ",
      "6 observations (got n = ", n, ").",
      call. = FALSE
    )
  }
  if (anyNA(X) || anyNA(Y)) {
    stop("hsic_test_gamma(): X and Y must not contain NA/NaN values.", call. = FALSE)
  }
  if (stats::sd(X) == 0 || stats::sd(Y) == 0) {
    return(list(stat = 0, p = 1))
  }

  width_x <- hsic_kernel_width(X)
  width_y <- hsic_kernel_width(Y)
  gm_x <- hsic_gram_matrix(X, width_x)
  gm_y <- hsic_gram_matrix(Y, width_y)

  K <- gm_x$K
  Kc <- gm_x$Kc
  L <- gm_y$K
  Lc <- gm_y$Kc

  stat <- sum(Kc * Lc) / n

  var_mat <- (Kc * Lc / 6)^2
  var_ <- (sum(var_mat) - sum(diag(var_mat))) / n / (n - 1)
  var_ <- 72 * (n - 4) * (n - 5) / n / (n - 1) / (n - 2) / (n - 3) * var_

  diag(K) <- 0
  diag(L) <- 0
  mu_x <- sum(K) / n / (n - 1)
  mu_y <- sum(L) / n / (n - 1)
  mean_ <- (1 + mu_x * mu_y - mu_x - mu_y) / n

  alpha <- mean_^2 / var_
  beta <- var_ * n / mean_

  p <- stats::pgamma(stat, shape = alpha, scale = beta, lower.tail = FALSE)

  list(stat = stat, p = p)
}
