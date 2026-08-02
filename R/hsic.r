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
# Shared by lingam_parce.r (BottomUpParceLiNGAM), lingam_rcd.r (RCD) --
# both univariate callers -- and lingam_resit.r (RESIT), which uses the
# multivariate (matrix) path.
# =============================================================================


#' Normalize HSIC input to an (n, d) numeric matrix
#'
#' Vectors become one-column matrices so that all downstream computations
#' can work row-wise, matching the upstream `X.reshape(-1, 1)` behavior.
#'
#' @param x numeric vector or matrix
#' @return numeric matrix (n x d)
#' @keywords internal
as_hsic_matrix <- function(x) {
  if (is.matrix(x)) x else matrix(as.numeric(x), ncol = 1L)
}


#' Pairwise squared Euclidean distances between rows
#'
#' The single-column case keeps the original `outer()` formulation so that
#' existing univariate callers (Parce / RCD) get bit-identical results; the
#' multivariate case uses the rowSums/tcrossprod expansion, which is only
#' reached by matrix inputs (new with the RESIT port).
#'
#' @param X numeric matrix (n x d)
#' @return n x n matrix of squared distances
#' @keywords internal
hsic_sqdist <- function(X) {
  if (ncol(X) == 1L) {
    x <- X[, 1L]
    return(outer(x, x, function(a, b) (a - b)^2))
  }
  sq <- rowSums(X^2)
  # x^2 + y^2 - 2xy can deviate slightly from the true distance through
  # floating-point cancellation: clamp negatives and pin the self-distance
  # diagonal to exactly 0 (so Gram diagonals are exactly 1).
  D <- pmax(outer(sq, sq, "+") - 2 * tcrossprod(X), 0)
  diag(D) <- 0
  D
}


#' Median-heuristic kernel width for HSIC
#'
#' Uses only the first 100 rows (not a random subsample) to keep the
#' O(n^2) pairwise-distance computation cheap, matching the upstream
#' implementation exactly. Multivariate input is measured by the squared
#' Euclidean distance between rows, as in the upstream `get_kernel_width()`.
#'
#' @param x numeric vector or matrix (n x d)
#' @return kernel width (scalar)
#' @keywords internal
hsic_kernel_width <- function(x) {
  X <- as_hsic_matrix(x)
  m <- min(nrow(X), 100L)
  D <- hsic_sqdist(X[seq_len(m), , drop = FALSE])
  d <- D[upper.tri(D)]
  d <- d[d > 0]
  if (length(d) == 0) {
    return(1.0)
  }
  sqrt(0.5 * stats::median(d))
}


#' Gaussian Gram matrix and its double-centered version
#'
#' Multivariate input is combined into a single Gaussian kernel over the
#' row-wise squared Euclidean distances (upstream behavior), not treated
#' column by column.
#'
#' @param x numeric vector or matrix (n x d)
#' @param width kernel width from [hsic_kernel_width()]
#' @return list(K = Gram matrix, Kc = centered Gram matrix)
#' @keywords internal
hsic_gram_matrix <- function(x, width) {
  X <- as_hsic_matrix(x)
  n <- nrow(X)
  H <- hsic_sqdist(X)
  K <- exp(-H / (2 * width^2))
  row_sums <- rowSums(K)
  col_sums <- colSums(K)
  total <- sum(K)
  Kc <- K - row_sums / n - rep(col_sums, each = n) / n + total / n^2
  list(K = K, Kc = Kc)
}


#' Precompute the per-variable parts of the HSIC gamma test
#'
#' Computes everything [hsic_test_gamma()] derives from one argument alone
#' (kernel width, centered Gram matrix, and the mean term `mu`), so callers
#' that test many pairs sharing a variable can reuse the O(n^2) Gram
#' computation. Validation (n >= 6, no NA, constant input) matches
#' [hsic_test_gamma()] exactly.
#'
#' @param x numeric vector or matrix (n x d)
#' @return list(n, is_const, Kc = centered Gram matrix, mu = mean term);
#'   `Kc`/`mu` are NULL when `is_const` is TRUE
#' @keywords internal
hsic_precompute <- function(x) {
  X <- as_hsic_matrix(x)
  n <- nrow(X)
  if (n < 6) {
    stop(
      "hsic_test_gamma(): the gamma-approximation test requires at least ",
      "6 observations (got n = ", n, ").",
      call. = FALSE
    )
  }
  if (anyNA(X)) {
    stop("hsic_test_gamma(): X and Y must not contain NA/NaN values.", call. = FALSE)
  }
  if (all(apply(X, 2, stats::sd) == 0)) {
    return(list(n = n, is_const = TRUE, Kc = NULL, mu = NULL))
  }
  width <- hsic_kernel_width(X)
  gm <- hsic_gram_matrix(X, width)
  K <- gm$K
  diag(K) <- 0
  mu <- sum(K) / n / (n - 1)
  list(n = n, is_const = FALSE, Kc = gm$Kc, mu = mu)
}


#' HSIC gamma test from precomputed parts
#'
#' Bit-identical to `hsic_test_gamma(X, Y)` when `pre_x`/`pre_y` are
#' [hsic_precompute()] of the same arguments in the same roles. The roles
#' matter at the last-ulp level: the `mean_` expression is not symmetric in
#' floating point, so `pre_x` must correspond to the first argument of the
#' original call being replaced.
#'
#' @param pre_x precomputed first argument (from [hsic_precompute()])
#' @param pre_y precomputed second argument
#' @return list(stat = HSIC test statistic, p = gamma-approximated p-value)
#' @keywords internal
hsic_gamma_from_pre <- function(pre_x, pre_y) {
  if (pre_x$n != pre_y$n) {
    stop(
      "hsic_test_gamma(): X and Y must have the same length ",
      "(number of observations).",
      call. = FALSE
    )
  }
  if (pre_x$is_const || pre_y$is_const) {
    return(list(stat = 0, p = 1))
  }
  n <- pre_x$n
  Kc <- pre_x$Kc
  Lc <- pre_y$Kc

  stat <- sum(Kc * Lc) / n

  var_mat <- (Kc * Lc / 6)^2
  var_ <- (sum(var_mat) - sum(diag(var_mat))) / n / (n - 1)
  var_ <- 72 * (n - 4) * (n - 5) / n / (n - 1) / (n - 2) / (n - 3) * var_

  mu_x <- pre_x$mu
  mu_y <- pre_y$mu
  mean_ <- (1 + mu_x * mu_y - mu_x - mu_y) / n

  alpha <- mean_^2 / var_
  beta <- var_ * n / mean_

  p <- stats::pgamma(stat, shape = alpha, scale = beta, lower.tail = FALSE)

  list(stat = stat, p = p)
}


#' Lazy per-column cache of [hsic_precompute()] results
#'
#' Returns a closure `f(k)` that computes `hsic_precompute(X[, k])` on first
#' use and returns the cached object afterwards. For callers that test many
#' pairs of columns of a fixed matrix.
#'
#' @param X numeric matrix whose columns will be tested
#' @return function(k) returning the precompute object for column k
#' @keywords internal
hsic_pre_col_cache <- function(X) {
  cache <- vector("list", ncol(X))
  function(k) {
    if (is.null(cache[[k]])) cache[[k]] <<- hsic_precompute(X[, k])
    cache[[k]]
  }
}


#' HSIC independence test with gamma approximation
#'
#' Faithful port of `hsic_test_gamma()` (hsic.py). O(n^2) in the sample size
#' because it forms the full n x n Gram matrices; not recommended for n in
#' the thousands.
#'
#' Either argument may be a matrix (n x d); its columns are combined into a
#' single Gaussian kernel over row-wise Euclidean distances, exactly as in
#' the upstream multivariate implementation (used by RESIT, which tests a
#' residual against the joint distribution of several predictors). The
#' variables are not standardized beforehand (upstream behavior), so with
#' wildly different column scales the largest-scale column dominates the
#' distance.
#'
#' The gamma-approximation variance estimator is only defined for n >= 6
#' (its closed form divides by `(n-1)(n-2)(n-3)`); below that this errors
#' instead of silently returning a NaN p-value that would otherwise
#' propagate into `NA`-valued rejection decisions in callers. A constant
#' input (zero variance in every column) carries no dependence information,
#' so it is treated as trivially independent (`p = 1`) rather than routed
#' through the degenerate kernel-width / centered-Gram-matrix computation.
#'
#' @param X numeric vector or matrix (n x d)
#' @param Y numeric vector or matrix with the same number of rows as X
#' @return list(stat = HSIC test statistic, p = gamma-approximated p-value)
#' @keywords internal
hsic_test_gamma <- function(X, Y) {
  X <- as_hsic_matrix(X)
  Y <- as_hsic_matrix(Y)
  if (nrow(X) != nrow(Y)) {
    stop(
      "hsic_test_gamma(): X and Y must have the same length ",
      "(number of observations).",
      call. = FALSE
    )
  }
  n <- nrow(X)
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
  if (all(apply(X, 2, stats::sd) == 0) || all(apply(Y, 2, stats::sd) == 0)) {
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
