# =============================================================================
# High-Dimensional Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
#
# Reference:
#   Wang, Y. S. and Drton, M. (2020). High-dimensional causal discovery under
#   non-Gaussianity. Biometrika, 107(1), 41-59.
# =============================================================================


#' Moore-Penrose pseudo-inverse via SVD
#'
#' Internal replacement for `numpy.linalg.pinv`, used to solve the
#' conditional regression coefficients from Gram-matrix submatrices in
#' [lingam_high_dim()]. Implemented in base R to avoid a new dependency.
#'
#' @param A numeric matrix
#' @param tol singular-value cutoff, relative to the largest singular value
#' @return the pseudo-inverse of `A`
#' @keywords internal
pinv <- function(A, tol = max(dim(A)) * .Machine$double.eps) {
  s <- svd(A)
  pos <- s$d > tol * max(s$d)
  if (!any(pos)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}


#' Unconditional non-Gaussianity statistic (tau)
#'
#' @param k moment degree
#' @param pa candidate parent column
#' @param ch candidate child column
#' @return a single non-negative numeric statistic
#' @keywords internal
calc_tau <- function(k, pa, ch) {
  abs(mean((pa^(k - 1)) * ch) * mean(pa^2) - mean(pa^k) * mean(pa * ch))
}


#' Conditional non-Gaussianity statistics (tau), minimized over all
#' conditioning subsets
#'
#' Per Wang & Drton (2020), the pruning statistic for a candidate
#' variable is the minimum of the conditional tau statistic over every
#' size-appropriate conditioning subset. Upstream HighDimDirectLiNGAM
#' (cdt15/lingam) has a `return` mis-indented inside its loop over
#' conditioning sets, so it only ever evaluates the first subset; this R
#' port intentionally does NOT replicate that bug and instead evaluates
#' every subset in `cond_sets`, so results differ numerically from the
#' Python package (see `dev/high-dim-direct-lingam-implementation.md`).
#'
#' @param Y data matrix
#' @param yty cached Gram matrix `t(Y) %*% Y`
#' @param pa index of the candidate parent variable (scalar)
#' @param ch indices of candidate child variables
#' @param k moment degree
#' @param cond_sets list of conditioning sets (each an integer vector of
#'   1-based column indices, always including `last_root`)
#' @param an_sets list, parallel to `cond_sets`, of ancestor-candidate sets
#'   (conditioning-set variables excluded from the corresponding subset)
#' @return numeric vector of length `ncol(Y)`, minimum statistic per variable
#' @keywords internal
calc_taus <- function(Y, yty, pa, ch, k, cond_sets, an_sets) {
  n <- nrow(Y)
  p <- ncol(Y)
  ret <- rep(1e10, p)

  for (z in seq_along(cond_sets)) {
    cond <- cond_sets[[z]]

    b <- pinv(yty[cond, cond, drop = FALSE]) %*% yty[cond, pa, drop = FALSE]
    resid <- as.numeric(Y[, pa] - Y[, cond, drop = FALSE] %*% b)

    resid_k_1 <- resid^(k - 1)
    resid_var <- mean(resid^2)
    resid_k <- mean(resid^k)

    for (j in ch) {
      value <- (1 / n) * sum(resid_k_1 * Y[, j]) * resid_var -
        resid_k * (1 / n) * sum(resid * Y[, j])
      ret[j] <- min(ret[j], abs(value))
    }

    for (a in an_sets[[z]]) {
      value <- (1 / n) * sum(resid_k_1 * Y[, a]) * resid_var -
        resid_k * (1 / n) * sum(resid * Y[, a])
      ret[a] <- min(ret[a], abs(value))
    }
  }

  ret
}


#' Enumerate conditioning subsets and compute pruning statistics for one
#' candidate variable
#'
#' @param Y data matrix
#' @param yty cached Gram matrix `t(Y) %*% Y`
#' @param i candidate variable (scalar, 1-based index)
#' @param j current candidate set (`psi`; `i` is removed internally)
#' @param K moment degree
#' @param last_root most recently fixed causal-order variable, or `NULL` on
#'   the first iteration
#' @param condition_set conditioning-set variables (already includes
#'   `last_root` when non-empty)
#' @param J assumed largest in-degree
#' @return numeric vector of length `ncol(Y)`
#' @importFrom utils combn
#' @keywords internal
get_prune_stats <- function(Y, yty, i, j, K, last_root, condition_set, J) {
  p <- ncol(Y)
  j <- setdiff(j, i)
  prune_stat <- rep(1e5, p)

  if (is.null(last_root)) {
    prune_stat[j] <- vapply(j, function(j_) calc_tau(K, Y[, i], Y[, j_]), numeric(1))
    return(prune_stat)
  }

  size_of_set <- min(J, length(condition_set))
  rest <- setdiff(condition_set, last_root)

  if (length(rest) == 1) {
    # NOTE: mirrors the upstream single-element special case. Also required
    # in R because combn(x, m) treats a length-1 numeric x as "choose from
    # seq_len(x)" instead of treating x as the single candidate element.
    raw_subsets <- list(rest)
  } else {
    m <- size_of_set - 1
    raw_subsets <- combn(rest, m, simplify = FALSE)
  }

  an_sets <- lapply(raw_subsets, function(x) setdiff(condition_set, x))
  cond_sets <- lapply(raw_subsets, function(x) c(last_root, x))

  calc_taus(Y, yty, i, j, K, cond_sets, an_sets)
}


#' Causal-order search for HighDimDirectLiNGAM
#'
#' @param X data matrix (n_samples x n_features)
#' @param J assumed largest in-degree
#' @param K moment degree
#' @param alpha pruning cutoff coefficient
#' @return integer vector, 1-based causal order (upstream-most first)
#' @keywords internal
high_dim_causal_order <- function(X, J, K, alpha) {
  Y <- X
  p <- ncol(Y)
  yty <- t(Y) %*% Y

  cut_off <- 0
  theta <- integer(0)
  psi <- seq_len(p)

  prune_stats <- matrix(1e5, p, p)
  diag(prune_stats) <- 0

  while (length(psi) > 1) {
    # Zero-init is safe here (unlike prune_stats' 1e5 sentinel) only because
    # every row idx in seq_along(psi) is unconditionally overwritten by
    # get_prune_stats() in the loop below before new_stats is combined with
    # prune_stats via pmin().
    new_stats <- matrix(0, nrow = length(psi), ncol = p)
    for (idx in seq_along(psi)) {
      v <- psi[idx]
      cond_set <- intersect(theta, which(prune_stats[v, ] > cut_off))
      if (length(theta) > 0) cond_set <- union(cond_set, theta[length(theta)])
      last_root <- if (length(theta) > 0) theta[length(theta)] else NULL

      new_stats[idx, ] <- get_prune_stats(Y, yty, v, psi, K, last_root, cond_set, J)
    }

    prune_stats[psi, ] <- pmin(prune_stats[psi, , drop = FALSE], new_stats)
    diag(prune_stats) <- 0

    sub <- prune_stats[psi, psi, drop = FALSE]
    max_taus <- apply(sub, 1, max)
    r <- psi[which.min(max_taus)]
    cut_off <- max(cut_off, min(max_taus) * alpha)

    theta <- c(theta, r)
    psi <- setdiff(psi, r)
  }

  c(theta, psi)
}


#' Adaptive LASSO with CV-selected lambda (n <= p route)
#'
#' Replicates upstream `_predict_adaptive_lasso` (StandardScaler + OLS
#' weights + `LassoLarsCV`), substituting `glmnet::cv.glmnet(alpha = 1)`
#' for `LassoLarsCV` (see [lingam_high_dim()] Details for the rationale).
#'
#' @param X original-scale data matrix
#' @param predictors indices of predictor variables (1-based)
#' @param target index of the target variable (1-based)
#' @param gamma exponent of the adaptive weights (fixed at 1.0 upstream)
#' @return coefficient vector, same length and order as `predictors`
#' @keywords internal
predict_adaptive_lasso_cv <- function(X, predictors, target, gamma = 1.0) {
  # glmnet requires >= 2 predictor columns; with a single predictor there is
  # nothing to prune, so fall back to plain OLS (mirrors fit_lasso()/
  # fit_ridge_reg() in fit_regression.r).
  if (length(predictors) == 1) {
    return(fit_ols(X[, target], X[, predictors, drop = FALSE]))
  }

  check_glmnet_available("high_dim (n <= p)")

  x_means <- colMeans(X)
  x_sds <- apply(X, 2, sd_pop)
  x_sds[x_sds < 1e-10] <- 1e-10
  X_std <- sweep(sweep(X, 2, x_means, "-"), 2, x_sds, "/")

  y_std <- X_std[, target]
  Xp_std <- X_std[, predictors, drop = FALSE]

  ols_fit <- stats::lm.fit(x = cbind(1, Xp_std), y = y_std)
  ols_coef <- as.numeric(ols_fit$coefficients[-1])
  ols_coef[is.na(ols_coef)] <- 0

  weight <- abs(ols_coef)^gamma

  Xp_weighted <- sweep(Xp_std, 2, weight, "*")

  # nfolds defaults to 10, which errors with a cryptic message when there are
  # fewer than 10 samples (plausible in the n <= p regime this function
  # targets); no upstream caller guarantees n_samples >= 10, so clamp
  # explicitly, respecting glmnet's minimum of nfolds = 3.
  n_samples <- nrow(Xp_weighted)
  nfolds <- max(3, min(10, n_samples - 1))
  cv_fit <- glmnet::cv.glmnet(x = Xp_weighted, y = y_std, alpha = 1, nfolds = nfolds)
  lasso_coef <- as.numeric(stats::coef(cv_fit, s = "lambda.min"))[-1]

  pruned_idx <- abs(lasso_coef * weight) > 0

  coef_out <- rep(0, length(predictors))
  if (any(pruned_idx)) {
    y_orig <- X[, target]
    Xp_orig <- X[, predictors[pruned_idx], drop = FALSE]
    refit <- stats::lm.fit(x = cbind(1, Xp_orig), y = y_orig)
    coef_out[pruned_idx] <- as.numeric(refit$coefficients[-1])
  }

  coef_out
}


#' Estimate the adjacency matrix by causal order using the n <= p route
#'
#' @param X original-scale data matrix
#' @param causal_order integer vector, 1-based causal order
#' @return adjacency matrix B (n_features x n_features)
#' @keywords internal
estimate_adjacency_matrix_high_dim_np <- function(X, causal_order) {
  n_features <- ncol(X)
  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (idx in seq_along(causal_order)) {
    if (idx == 1) next
    target <- causal_order[idx]
    predictors <- causal_order[1:(idx - 1)]
    B[target, predictors] <- predict_adaptive_lasso_cv(X, predictors, target)
  }

  B
}


#' High-Dimensional Direct LiNGAM
#'
#' A variant of Direct LiNGAM for high-dimensional data (large `p`, or
#' `p > n`). Causal order search is based on moment statistics of
#' non-Gaussianity rather than pairwise independence measures, which is
#' considerably faster for many variables. Unlike [lingam_direct()], the
#' algorithm is deterministic (no random restarts).
#'
#' @param X Numeric matrix (n_samples x n_features), data frame or matrix
#' @param J Assumed largest in-degree (single integer, must be >= 3)
#' @param K Degree of the moment used to measure non-Gaussianity (single
#'   integer, must be >= 1)
#' @param alpha Cutoff coefficient for pruning away false parents (single
#'   numeric value in `[0, 1]`)
#' @param estimate_adj_mat Whether to estimate the adjacency matrix (single
#'   logical value). If `FALSE`, causal-order search still runs but
#'   `adjacency_matrix` is returned as an NA-filled matrix (not `NULL`, so
#'   downstream S3 methods keep working).
#' @return A `LingamResult` object (list), the same class returned by
#'   [lingam_direct()], containing:
#' * `adjacency_matrix`: adjacency matrix B (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i).** Zero elements indicate no causal relationship.
#'   All-`NA` when `estimate_adj_mat = FALSE`.
#' * `causal_order`: estimated causal order (integer vector of 1-based
#'   indices). Earlier elements are more upstream.
#'
#' @details
#' When `n_samples <= n_features`, the adjacency matrix cannot be estimated
#' with the usual BIC-based Adaptive LASSO, so this function falls back to a
#' cross-validated LASSO (`glmnet::cv.glmnet`) after emitting a warning.
#' The upstream Python implementation uses `LassoLarsCV` for this fallback;
#' `cv.glmnet` follows the same cross-validation design but is not
#' numerically identical (different solver: coordinate descent vs. LARS).
#'
#' The pruning statistic used during causal-order search (Wang & Drton 2020)
#' is the minimum of a conditional non-Gaussianity statistic over every
#' size-appropriate conditioning subset. The upstream Python implementation
#' (cdt15/lingam) has a `return` statement mis-indented inside the loop over
#' conditioning subsets, causing it to only ever evaluate the first subset.
#' This R implementation intentionally does not replicate that bug and
#' evaluates every subset, so `causal_order` and `adjacency_matrix` are not
#' numerically identical to the upstream Python package.
#'
#' @references
#' Wang, Y. S. and Drton, M. (2020). High-dimensional causal discovery under
#' non-Gaussianity. Biometrika, 107(1), 41-59.
#'
#' @export
#' @examples
#' sample <- generate_lingam_sample_6(n = 300, seed = 1)
#' result <- lingam_high_dim(sample$data)
#' result$causal_order
#' round(result$adjacency_matrix, 3)
#'
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   # n <= p: falls back to cross-validated LASSO with a warning
#'   wide_sample <- generate_lingam_large_sample(p = 30, n = 25, seed = 1)
#'   wide_result <- lingam_high_dim(wide_sample$data)
#'   wide_result$causal_order
#' }
#' }
lingam_high_dim <- function(X,
                            J = 3L,
                            K = 4L,
                            alpha = 0.5,
                            estimate_adj_mat = TRUE) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 2) stop("X must have at least 2 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names

  if (!is.numeric(J) || length(J) != 1 || is.na(J) || J != as.integer(J) || J <= 2) {
    stop("J must be a single integer >= 3.", call. = FALSE)
  }
  J <- as.integer(J)

  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K != as.integer(K) || K < 1) {
    stop("K must be a single integer >= 1.", call. = FALSE)
  }
  K <- as.integer(K)

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (!is.logical(estimate_adj_mat) || length(estimate_adj_mat) != 1 || is.na(estimate_adj_mat)) {
    stop("estimate_adj_mat must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  n_samples <- nrow(X)
  n_features <- ncol(X)

  causal_order <- high_dim_causal_order(X, J = J, K = K, alpha = alpha)

  if (!estimate_adj_mat) {
    B <- matrix(NA_real_, nrow = n_features, ncol = n_features)
    colnames(B) <- rownames(B) <- colnames(X)
    result <- list(adjacency_matrix = B, causal_order = causal_order)
    class(result) <- "LingamResult"
    return(result)
  }

  if (n_samples <= n_features) {
    warning(
      "Since n_samples <= n_features, the adjacency matrix is estimated with ",
      "cross-validated lasso (cv.glmnet) instead of BIC-based lambda selection.",
      call. = FALSE
    )
    B <- estimate_adjacency_matrix_high_dim_np(X, causal_order)
  } else {
    B <- estimate_adjacency_matrix(X, causal_order,
      method = "adaptive_lasso", lambda = "BIC", init_method = "ols"
    )
  }

  colnames(B) <- rownames(B) <- colnames(X)
  result <- list(adjacency_matrix = B, causal_order = causal_order)
  class(result) <- "LingamResult"
  result
}
