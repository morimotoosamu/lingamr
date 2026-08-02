# =============================================================================
# Direct LiNGAM - Causal order search (pwling / kernel) and prior knowledge handling
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
# =============================================================================


# Above this sample size the exact O(n^3) kernel path switches to the
# incomplete-Cholesky low-rank approximation. The kappa/sigma parameters are
# tied to this switch, so search_causal_order_kernel() and f_correlation()
# (in f_correlation.r) must both use this same constant.
KERNEL_LOWRANK_MAX_N <- 1000L


#' Extract partial orders from prior knowledge
#' @param pk Prior knowledge matrix (NaN = unknown)
#' @return matrix (n x 2), each row is a (from, to) partial order
#' @keywords internal
extract_partial_orders <- function(pk) {
  # Pairs with a path (pk == 1)
  path_idx <- which(pk == 1, arr.ind = TRUE)
  # Pairs with no path (pk == 0)
  no_path_idx <- which(pk == 0, arr.ind = TRUE)

  # --- Consistency check for path pairs ---
  if (nrow(path_idx) > 0) {
    check_pairs <- rbind(path_idx, path_idx[, 2:1, drop = FALSE])
    dup <- duplicated(check_pairs) | duplicated(check_pairs, fromLast = TRUE)
    if (any(dup)) {
      bad <- unique(check_pairs[dup, , drop = FALSE])
      stop(paste(
        "The prior knowledge contains inconsistencies at indices:",
        paste(apply(bad, 1, function(r) paste0("(", r[1], ",", r[2], ")")),
          collapse = ", "
        )
      ))
    }
  }

  # --- Deduplicate no-path pairs ---
  if (nrow(no_path_idx) > 0) {
    check_pairs2 <- rbind(no_path_idx, no_path_idx[, 2:1, drop = FALSE])
    # Find and exclude pairs that have 0 in both directions
    pair_key <- paste(check_pairs2[, 1], check_pairs2[, 2], sep = ",")
    tbl <- table(pair_key)
    dup_keys <- names(tbl[tbl > 1])
    # Exclude entries of no_path_idx that exist in both directions
    no_path_key <- paste(no_path_idx[, 1], no_path_idx[, 2], sep = ",")
    keep <- !(no_path_key %in% dup_keys)
    no_path_idx <- no_path_idx[keep, , drop = FALSE]
  }

  # Combine path_pairs with no_path_pairs[:, [1,0]]
  combined <- matrix(nrow = 0, ncol = 2)
  if (nrow(path_idx) > 0) {
    combined <- rbind(combined, path_idx)
  }
  if (nrow(no_path_idx) > 0) {
    combined <- rbind(combined, no_path_idx[, 2:1, drop = FALSE])
  }

  if (nrow(combined) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  combined <- unique(combined)
  # [to, from] -> [from, to]
  result <- combined[, 2:1, drop = FALSE]
  colnames(result) <- NULL
  return(result)
}


#' Residual (residual when xi is regressed on xj)
#' Compute the residual vector
#' @param xi Target variable vector
#' @param xj Explanatory variable vector
#' @return Residual vector after regression
#' @keywords internal
residual_vec <- function(xi, xj) {
  xi_c <- xi - mean(xi)
  xj_c <- xj - mean(xj)
  beta <- sum(xi_c * xj_c) / sum(xj_c * xj_c)
  xi - beta * xj
}


#' Maximum-entropy approximation of entropy
#' @param u Input vector
#' @return Approximate entropy value
#' @keywords internal
entropy_approx <- function(u) {
  n <- length(u)
  k1 <- 79.047
  k2 <- 7.4129
  gamma <- 0.37457
  (1 + log(2 * pi)) / 2 -
    k1 * (sum(log(cosh(u))) / n - gamma)^2 -
    k2 * (sum(u * exp(-u^2 / 2)) / n)^2
}


#' Search for candidate variables
#' @param U Set of currently undetermined variables
#' @param Aknw Prior knowledge matrix
#' @param apply_prior_knowledge_softly Whether to apply prior knowledge softly
#' @param partial_orders Extracted partial orders
#' @return list(Uc, Vj)
#' @keywords internal
search_candidate <- function(U, Aknw, apply_prior_knowledge_softly, partial_orders) {
  # No prior knowledge

  if (is.null(Aknw)) {
    return(list(Uc = U, Vj = integer(0)))
  }

  # --- Hard application ---
  if (!apply_prior_knowledge_softly) {
    if (!is.null(partial_orders) && nrow(partial_orders) > 0) {
      Uc <- setdiff(U, partial_orders[, 2])
      if (length(Uc) == 0) Uc <- U
      return(list(Uc = Uc, Vj = integer(0)))
    } else {
      return(list(Uc = U, Vj = integer(0)))
    }
  }

  # --- Soft application ---
  # Search for exogenous variables
  Uc <- integer(0)
  for (j in U) {
    index <- setdiff(U, j)
    # A row containing NA (unknown) yields NA from sum -> isTRUE() becomes FALSE and it is dropped from candidates
    # (same behavior as Python's NaN.sum() == 0 evaluating to False)
    if (isTRUE(sum(Aknw[j, index]) == 0)) {
      Uc <- c(Uc, j)
    }
  }

  # Search for endogenous variables -> narrow down candidates
  if (length(Uc) == 0) {
    U_end <- integer(0)
    for (j in U) {
      index <- setdiff(U, j)
      s <- sum(Aknw[j, index], na.rm = TRUE)
      if (!is.na(s) && s > 0) {
        U_end <- c(U_end, j)
      }
    }
    # Sink features
    for (i in U) {
      index <- setdiff(U, i)
      if (isTRUE(sum(Aknw[index, i]) == 0)) {
        U_end <- c(U_end, i)
      }
    }
    Uc <- setdiff(U, unique(U_end))
    if (length(Uc) == 0) Uc <- U
  }

  # Build V^(j)
  Vj <- integer(0)
  for (i in U) {
    if (i %in% Uc) next
    if (isTRUE(sum(Aknw[i, Uc]) == 0)) {
      Vj <- c(Vj, i)
    }
  }

  return(list(Uc = Uc, Vj = Vj))
}


#' Causal order search via pwling
#' @param X Data matrix
#' @param U Indices of all variables
#' @param Uc Indices of candidate variables
#' @param Vj Variable set based on prior knowledge
#' @return Index of the selected variable
#' @keywords internal
search_causal_order_pwling <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) return(Uc[1])
  n <- nrow(X)
  p <- ncol(X)

  # --- Standardize all columns at once and precompute each column's entropy (independent of pairs) ---
  X_std <- matrix(0, nrow = n, ncol = p)
  H <- numeric(p)
  for (k in U) {
    xk <- X[, k]
    xk <- xk - sum(xk) / n
    X_std[, k] <- xk / sqrt(sum(xk * xk) / n)
    H[k] <- entropy_approx(X_std[, k])
  }

  # Since the data is standardized, the correlation matrix is crossprod / n (computed at once via BLAS).
  # The regression coefficient beta and the population SD of the residual sqrt(1 - r^2) follow analytically.
  R <- crossprod(X_std[, U, drop = FALSE]) / n
  pos <- integer(p)
  pos[U] <- seq_along(U)

  in_Uc <- logical(p)
  in_Uc[Uc] <- TRUE
  in_Vj <- logical(p)
  in_Vj[Vj] <- TRUE

  M_acc <- numeric(p)
  for (i in Uc) {
    xi_std <- X_std[, i]
    for (j in U) {
      if (i == j) next
      # The pairwise mutual-information difference is antisymmetric
      # (dm_ji = -dm_ij), so for pairs where both are candidates, compute
      # once on the i < j side and add to both.
      if (in_Uc[j] && j < i) next
      xj_std <- X_std[, j]
      r_ij <- R[pos[i], pos[j]]
      sd_r <- sqrt(max(0, 1 - r_ij^2))

      H_ri_j <- if (in_Vj[i] && in_Uc[j]) {
        H[i]
      } else {
        entropy_approx((xi_std - r_ij * xj_std) / sd_r)
      }
      H_rj_i <- if (in_Vj[j] && in_Uc[i]) {
        H[j]
      } else {
        entropy_approx((xj_std - r_ij * xi_std) / sd_r)
      }

      dm <- (H[j] + H_ri_j) - (H[i] + H_rj_i)
      M_acc[i] <- M_acc[i] + min(0, dm)^2
      if (in_Uc[j]) M_acc[j] <- M_acc[j] + min(0, -dm)^2
    }
  }
  return(Uc[which.max(-M_acc[Uc])])
}


#' Kernel-based mutual information: precomputation for variable 1
#'
#' Computes the matrix `E1 = tmp1^-1 K1` (`tmp1 = K1 + n*kappa/2 * I`) used in
#' `kernel_mi_core()`. It only needs to be called once per candidate variable,
#' avoiding per-pair recomputation.
#'
#' @param x Vector of variable 1
#' @param kappa Regularization parameter
#' @param sigma Width of the Gaussian kernel
#' @return Matrix E1 (n x n)
#' @keywords internal
kernel_mi_prepare <- function(x, kappa, sigma) {
  n <- length(x)
  K <- exp(-1 / (2 * sigma^2) * outer(x, x, "-")^2)
  c0 <- n * kappa / 2
  tmp <- K
  diag(tmp) <- diag(tmp) + c0
  # tmp and K commute, so E = tmp^-1 K = I - c0 * tmp^-1 (symmetric)
  E <- -c0 * chol2inv(chol(tmp))
  diag(E) <- diag(E) + 1
  E
}


#' Kernel-based mutual information: core
#'
#' The target quantity is the difference of logdets of 2n x 2n matrices, but via
#' the block structure and the Schur complement it can be computed equivalently
#' using only an n x n Cholesky decomposition:
#' `MI = -1/2 * (logdet(tmp2^2 - K2 K1 tmp1^-2 K1 K2) - logdet(tmp2^2))`
#'
#' @param E1 Variable-1 matrix precomputed by `kernel_mi_prepare()`
#' @param x2 Vector of variable 2
#' @param kappa Regularization parameter
#' @param sigma Width of the Gaussian kernel
#' @return Mutual information
#' @keywords internal
kernel_mi_core <- function(E1, x2, kappa, sigma) {
  n <- length(x2)
  K2 <- exp(-1 / (2 * sigma^2) * outer(x2, x2, "-")^2)
  tmp2 <- K2
  diag(tmp2) <- diag(tmp2) + n * kappa / 2
  W <- E1 %*% K2                       # = tmp1^-1 K1 K2
  S <- crossprod(tmp2) - crossprod(W)  # Schur complement (tmp2 is symmetric, so crossprod = tmp2^2)
  logdet_S <- 2 * sum(log(diag(chol(S))))
  logdet_tmp2_sq <- 4 * sum(log(diag(chol(tmp2))))
  (-1 / 2) * (logdet_S - logdet_tmp2_sq)
}


#' Gaussian kernel via pivoted incomplete Cholesky
#'
#' Low-rank approximation `G` (n x d) with `tcrossprod(G) ~= K`, where
#' `K[i, j] = exp(-(x[i] - x[j])^2 / (2*sigma^2))`, via the greedy pivoted
#' incomplete Cholesky of Fine & Scheinberg (2001) / Bach & Jordan (2002,
#' kernel-ICA). Each step picks the index with the largest diagonal
#' residual and computes only the kernel column for that pivot, so the
#' full n x n Gram matrix is never formed. The Gaussian kernel's diagonal
#' is 1, so residuals start at 1 without evaluating `K`.
#' @param x Input vector (length n)
#' @param sigma Width of the Gaussian kernel
#' @param tol Stop once the largest remaining diagonal residual is at or below this
#' @param max_rank Upper bound on the approximation rank
#' @return n x d matrix `G` with `tcrossprod(G) ~= K`
#' @keywords internal
incomplete_cholesky_gauss <- function(x, sigma, tol = 1e-4, max_rank = min(length(x), 200)) {
  n <- length(x)
  max_rank <- min(max_rank, n)
  d_vec <- rep(1, n) # diag(K) == 1 for the Gaussian kernel
  perm <- seq_len(n)
  G <- matrix(0, n, max_rank)
  rank <- 0L

  for (k in seq_len(max_rank)) {
    remaining <- k:n
    pivot_rel <- which.max(d_vec[perm[remaining]])
    pivot_pos <- remaining[pivot_rel]
    if (d_vec[perm[pivot_pos]] <= tol) break

    if (pivot_pos != k) {
      perm[c(k, pivot_pos)] <- perm[c(pivot_pos, k)]
    }
    pk <- perm[k]
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

  if (rank == max_rank && max_rank < n) {
    unselected <- perm[(max_rank + 1):n]
    if (max(d_vec[unselected]) > tol) {
      warning(
        "incomplete_cholesky_gauss: max_rank reached before tol; ",
        "the low-rank kernel approximation may be inaccurate"
      )
    }
  }

  G[, seq_len(rank), drop = FALSE]
}


#' Kernel-based mutual information: low-rank precomputation for variable 1
#'
#' Low-rank counterpart of `kernel_mi_prepare()`. Via the Woodbury identity,
#' `E1 = tmp1^-1 K1` collapses to `G1 %*% M1^-1 %*% t(G1)` where
#' `M1 = c0*I + t(G1) %*% G1`, so the n x n matrix `E1` never needs to be
#' formed; only the d x d inverse `M1^-1` does.
#' @param x Vector of variable 1
#' @param kappa Regularization parameter
#' @param sigma Width of the Gaussian kernel
#' @return list(G, Minv, A, c0) describing `E1` in factored form
#' @keywords internal
kernel_mi_prepare_lowrank <- function(x, kappa, sigma) {
  n <- length(x)
  c0 <- n * kappa / 2
  G <- incomplete_cholesky_gauss(x, sigma)
  d <- ncol(G)
  A <- crossprod(G) # t(G) %*% G
  Minv <- solve(c0 * diag(d) + A)
  # Minv %*% A %*% Minv is independent of variable 2; precompute it here so
  # kernel_mi_core_lowrank() does a single d1 x d2 product per candidate.
  Q <- Minv %*% A %*% Minv
  list(G = G, Minv = Minv, A = A, c0 = c0, Q = Q)
}


#' Kernel-based mutual information: low-rank core
#'
#' Low-rank counterpart of `kernel_mi_core()`. With `K2 ~= G2 %*% t(G2)`,
#' both `S = tmp2^2 - t(W) %*% W` and `tmp2^2` reduce to `c0^2*I + G2 %*% C %*% t(G2)`
#' for a d2 x d2 matrix `C`, via the Woodbury identity on `W` and the matrix
#' determinant lemma on the resulting rank-d2 update. The `n*log(c0)` terms
#' in `logdet(S)` and `logdet(tmp2^2)` cancel in their difference, so only
#' d2 x d2 matrices remain.
#' @param prep1 Output of `kernel_mi_prepare_lowrank()` for variable 1
#' @param x2 Vector of variable 2
#' @param kappa Regularization parameter
#' @param sigma Width of the Gaussian kernel
#' @return Mutual information
#' @keywords internal
kernel_mi_core_lowrank <- function(prep1, x2, kappa, sigma) {
  c0 <- prep1$c0
  G2 <- incomplete_cholesky_gauss(x2, sigma)
  d2 <- ncol(G2)
  A2 <- crossprod(G2)
  P <- crossprod(prep1$G, G2) # t(G1) %*% G2, d1 x d2
  D <- crossprod(P, prep1$Q %*% P)
  Cprime <- 2 * c0 * diag(d2) + A2 - D

  # logdet(S) and logdet(tmp2^2) each carry a 2*n*log(c0) term that cancels
  # in the difference below, so only the d2 x d2 determinants remain.
  logdet_S <- as.numeric(determinant(diag(d2) + Cprime %*% A2 / c0^2, logarithm = TRUE)$modulus)
  logdet_tmp2 <- as.numeric(determinant(diag(d2) + A2 / c0, logarithm = TRUE)$modulus)
  (-1 / 2) * (logdet_S - 2 * logdet_tmp2)
}


#' Causal order search via the kernel method
#' @param X Data matrix
#' @param U All variables
#' @param Uc Candidate variables
#' @param Vj Prior knowledge set
#' @return Index of the selected variable
#' @keywords internal
search_causal_order_kernel <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) {
    return(Uc[1])
  }

  n <- nrow(X)
  # Above this n, the exact O(n^3) kernel computation dominates runtime, so
  # switch to the incomplete-Cholesky low-rank path
  # (kernel_mi_prepare_lowrank() / kernel_mi_core_lowrank()).
  use_lowrank <- n > KERNEL_LOWRANK_MAX_N
  if (use_lowrank) {
    param <- c(2e-3, 0.5)
  } else {
    param <- c(2e-2, 1.0)
  }

  kappa <- param[1]
  sigma <- param[2]
  Tkernels <- numeric(length(Uc))

  for (idx in seq_along(Uc)) {
    j <- Uc[idx]
    # The kernel matrix / inverse for candidate j is invariant across the inner loop, so compute it once
    E1 <- if (use_lowrank) {
      kernel_mi_prepare_lowrank(X[, j], kappa, sigma)
    } else {
      kernel_mi_prepare(X[, j], kappa, sigma)
    }
    Tkernel <- 0
    for (i in U) {
      if (i == j) next
      ri_j <- if (i %in% Vj && j %in% Uc) {
        X[, i]
      } else {
        residual_vec(X[, i], X[, j])
      }
      Tkernel <- Tkernel + if (use_lowrank) {
        kernel_mi_core_lowrank(E1, ri_j, kappa, sigma)
      } else {
        kernel_mi_core(E1, ri_j, kappa, sigma)
      }
    }
    Tkernels[idx] <- Tkernel
  }

  return(Uc[which.min(Tkernels)])
}
