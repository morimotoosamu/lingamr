# =============================================================================
# LiM (LiNGAM for Mixed data) - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam  (lingam/lim.py, lingam/utils/__init__.py)
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
#   Zeng Y, Shimizu S, Matsui H, Sun F. Causal discovery for linear mixed
#   data. In: Proceedings of the First Conference on Causal Learning and
#   Reasoning (CLeaR 2022). PMLR 177, pp. 994-1009, 2022.
#
# Internal orientation note: all computation below uses the *NOTEARS-style*
# orientation W[i, j] = i -> j (M = X %*% W), matching the Python source
# exactly. Only the final result returned to the user is transposed to the
# lingamr convention B[i, j] = j -> i (see lingam_lim() at the bottom).
# =============================================================================


#' Numerically stable log(1 + exp(M)) (elementwise)
#' @keywords internal
lim_log1pexp <- function(M) {
  out <- M
  small <- M <= 33
  out[small] <- log1p(exp(M[small]))
  out
}


#' Numerically stable log(cosh(r)) (elementwise)
#' @keywords internal
lim_logcosh <- function(r) {
  abs(r) + log1p(exp(-2 * abs(r))) - log(2)
}


#' Flatten a matrix in row-major (numpy "C") order
#' @keywords internal
lim_flatten_rowmajor <- function(M) {
  as.vector(t(M))
}


#' Convert the doubled parameter vector (length 2*d*d) into a d x d matrix W
#'
#' Mirrors Python's `w[:d*d].reshape([d, d])` (row-major reshape), so `w` and
#' the reverse operation `lim_flatten_rowmajor()` must stay consistent.
#' @keywords internal
lim_adj <- function(w, d) {
  matrix(w[seq_len(d * d)] - w[(d * d + 1):(2 * d * d)], nrow = d, ncol = d, byrow = TRUE)
}


#' Mixed loss (logistic for discrete columns, log-cosh/Laplace for continuous
#' columns) and its gradient. `con` is a length-d vector, 1 = continuous,
#' 0 = discrete. `W_dis_mask` / `W_con_mask` are precomputed d x d 0/1 masks.
#' @keywords internal
lim_loss_mixed <- function(W, X, con, W_dis_mask, W_con_mask) {
  n <- nrow(X)
  M <- X %*% W
  R <- X - M

  col_mask_dis <- 1 - con
  col_mask_con <- con

  loss_dis <- sum(sweep(lim_log1pexp(M) - X * M, 2, col_mask_dis, `*`))
  loss_con <- sum(sweep(-lim_logcosh(R), 2, col_mask_con, `*`))
  loss <- -(loss_dis + loss_con) / n

  G_dis <- (t(X) %*% (stats::plogis(M) - X)) * W_dis_mask
  G_con <- -(t(X) %*% tanh(R)) * W_con_mask
  G_loss <- (G_dis + G_con) / n

  list(loss = loss, G_loss = G_loss)
}


#' Acyclicity constraint h(W) and the matrix power term used in its gradient
#' (Yu et al. 2019 formulation, as used by the Python source)
#' @keywords internal
lim_h <- function(W, d) {
  M_h <- diag(d) + (W * W) / d
  E <- M_h
  if (d > 2) {
    for (k in seq_len(d - 2)) E <- E %*% M_h
  }
  h <- sum(t(E) * M_h) - d
  list(h = h, E = E)
}


#' Value and gradient of the augmented-Lagrangian objective for the doubled
#' parameter vector w (length 2*d*d)
#' @keywords internal
lim_obj_grad <- function(w, d, X, con, W_dis_mask, W_con_mask, rho, alpha, lambda1) {
  W <- lim_adj(w, d)
  lm <- lim_loss_mixed(W, X, con, W_dis_mask, W_con_mask)
  hh <- lim_h(W, d)

  obj <- lm$loss + 0.5 * rho * hh$h^2 + alpha * hh$h + lambda1 * sum(w)
  G_smooth <- lm$G_loss + (rho * hh$h + alpha) * t(hh$E) * W * 2
  g_flat <- lim_flatten_rowmajor(G_smooth)
  g_obj <- c(g_flat + lambda1, -g_flat + lambda1)

  list(obj = obj, grad = g_obj)
}


#' Global (NOTEARS-style) optimization phase: augmented-Lagrangian outer loop
#' with L-BFGS-B inner solves. Returns the thresholded W (i -> j orientation).
#' @keywords internal
lim_global_optimize <- function(X, con, lambda1, max_iter, h_tol, rho_max, w_threshold) {
  d <- ncol(X)

  W_dis_mask <- matrix(0, d, d)
  disc_idx <- which(con == 0)
  if (length(disc_idx) > 0) {
    W_dis_mask[disc_idx, ] <- 1
    W_dis_mask[, disc_idx] <- 1
  }
  W_con_mask <- matrix(0, d, d)
  cont_idx <- which(con == 1)
  if (length(cont_idx) > 0) {
    W_con_mask[cont_idx, ] <- 1
    W_con_mask[, cont_idx] <- 1
  }

  diag_idx <- (seq_len(d) - 1L) * d + seq_len(d)
  lower <- rep(0, 2 * d * d)
  upper <- rep(Inf, 2 * d * d)
  upper[diag_idx] <- 0
  upper[d * d + diag_idx] <- 0

  w_est <- stats::runif(2 * d * d)
  rho <- 1.0
  alpha <- 0.0
  h <- Inf

  for (iter in seq_len(max_iter)) {
    w_new <- NULL
    h_new <- NULL
    while (rho < rho_max) {
      fn <- function(w) lim_obj_grad(w, d, X, con, W_dis_mask, W_con_mask, rho, alpha, lambda1)$obj
      gr <- function(w) lim_obj_grad(w, d, X, con, W_dis_mask, W_con_mask, rho, alpha, lambda1)$grad
      sol <- stats::optim(w_est, fn = fn, gr = gr, method = "L-BFGS-B",
                           lower = lower, upper = upper)
      w_new <- sol$par
      h_new <- lim_h(lim_adj(w_new, d), d)$h
      if (h_new >= 0.25 * h) {
        rho <- rho * 10
      } else {
        break
      }
    }
    w_est <- w_new
    h <- h_new
    alpha <- alpha + rho * h
    if (h <= h_tol && h != 0) break
    if (rho >= rho_max * 1e-6 && h > 1e5) {
      # avoid the full graph
      w_est <- stats::runif(2 * d * d)
      rho <- 1.0
    } else if (rho >= rho_max) {
      break
    }
    if (sum(abs(lim_adj(w_est, d))) < 0.09) {
      # avoid the zero matrix
      w_est <- stats::runif(2 * d * d)
    }
  }

  W_est <- lim_adj(w_est, d)
  W_est[abs(W_est) < w_threshold] <- 0
  W_est
}


#' Local log-likelihood of continuous component i (super-Gaussian disturbance)
#'
#' Direct port of `likelihood_i`/`log_p_super_gaussian`/`variance_i` from
#' `lingam/utils/__init__.py` (lines 866-939 of the Python source).
#'
#' A (near-)zero residual variance means component i is fit deterministically
#' by its parents; `-n * log(sqrt(var_i))` would then diverge to `Inf` and
#' make this candidate model look artificially superior. Such candidates are
#' instead scored with a large finite penalty so they lose every BIC
#' comparison in [lim_bic_loss()]/[lim_local_search()] without introducing
#' `Inf`/`NaN` into the running score.
#' @keywords internal
lim_likelihood_i <- function(X, i, b_full, b0) {
  n <- nrow(X)
  fitted <- as.vector(X %*% b_full)
  e <- X[, i] - fitted
  var_i <- mean((e - mean(e))^2)
  if (!is.finite(var_i) || var_i < .Machine$double.eps) {
    return(-1e10)
  }
  resid_std <- (X[, i] - fitted - b0) / sqrt(var_i)
  sum(-sqrt(2) * abs(resid_std) - 0.35) - n * log(sqrt(var_i))
}


#' BIC-penalized DAG score (sign-flipped, i.e. a "loss") used during the local
#' search phase. `W` is interpreted as a 0/1 skeleton in the i -> j
#' orientation (parents of j are the nonzero rows of column j).
#' @keywords internal
lim_bic_loss <- function(W, X, is_continuous) {
  n <- nrow(X)
  d <- ncol(X)
  n_edges <- sum(W != 0)
  penalty <- log(n) / 2 * (n_edges + d)
  total_score <- -penalty

  for (i in seq_len(d)) {
    parents_i <- which(W[, i] != 0)
    xi <- X[, i]

    if (!is_continuous[i]) {
      if (length(parents_i) == 0) {
        tab <- table(xi)
        counts <- as.numeric(tab)
        total_score <- total_score + sum(counts * (log(counts) - log(sum(counts))))
      } else {
        df <- data.frame(y = xi, X[, parents_i, drop = FALSE])
        fit <- suppressWarnings(stats::glm(y ~ ., data = df, family = stats::binomial()))
        ll <- as.numeric(stats::logLik(fit))
        # quasi-separation (or a non-converged fit) can drive the logLik to
        # +-Inf/NaN; treat such a candidate as maximally bad rather than
        # letting Inf/NaN propagate into total_score.
        if (!is.finite(ll)) ll <- -1e10
        total_score <- total_score + ll
      }
    } else {
      b_full <- numeric(d)
      if (length(parents_i) > 0) {
        fit <- stats::lm(xi ~ X[, parents_i, drop = FALSE])
        coefs <- stats::coef(fit)
        if (anyNA(coefs)) {
          # rank-deficient design (collinear/constant parent columns); reject
          # this candidate instead of letting NA propagate through the score
          total_score <- total_score - 1e10
          next
        }
        b0 <- coefs[1]
        b_full[parents_i] <- coefs[-1]
      } else {
        b0 <- mean(xi)
      }
      total_score <- total_score + lim_likelihood_i(X, i, b_full, b0)
    }
  }

  -total_score
}


#' Local search phase: direction reversal, pruning, and edge addition
#' (Python source lines 324-397). Returns W in the i -> j orientation.
#' @keywords internal
lim_local_search <- function(W_est, X, con, is_continuous, d, h_tol) {
  aa <- lim_bic_loss(W_est, X, is_continuous)
  nz <- which(W_est != 0, arr.ind = TRUE)
  n_edges <- nrow(nz)

  if (n_edges > 15) {
    warning(
      "Local search skipped: too many edges (2^|E| direction combinations); returning the global solution.",
      call. = FALSE
    )
    return(W_est)
  }

  W_min_lss <- W_est

  # (a) direction reversal: exhaustively try flipping each edge's direction
  if (n_edges > 0) {
    combos <- as.matrix(expand.grid(rep(list(c(0, 1)), n_edges)))
    for (r in seq_len(nrow(combos))) {
      setting <- combos[r, ]
      if (all(setting == 0)) next  # the "unchanged" candidate is W_est itself
      W_tmp <- W_est
      for (k in seq_len(n_edges)) {
        if (setting[k] == 1) {
          W_tmp[nz[k, 1], nz[k, 2]] <- 0
          W_tmp[nz[k, 2], nz[k, 1]] <- 1
        }
      }
      lss <- lim_bic_loss(W_tmp, X, is_continuous)
      if (lss < aa && lim_h(W_tmp, d)$h < h_tol) {
        W_min_lss <- W_tmp
        aa <- lss
      }
    }
  }

  # (b) pruning: try deleting each remaining edge one at a time
  if (d > 2 && n_edges > (d - 1)) {
    W0 <- W_min_lss
    I_delete <- which(W0 != 0, arr.ind = TRUE)
    for (k in seq_len(nrow(I_delete))) {
      W_tmp <- W0
      W_tmp[I_delete[k, 1], I_delete[k, 2]] <- 0
      lss <- lim_bic_loss(W_tmp, X, is_continuous)
      if (lss < aa && lim_h(W_tmp, d)$h < h_tol) {
        W_min_lss <- W_tmp
        aa <- lss
      }
    }
    if (!identical(W_est != 0, W_min_lss != 0)) {
      W0b <- W_est
      for (k in seq_len(nrow(nz))) {
        W_tmp <- W0b
        W_tmp[nz[k, 1], nz[k, 2]] <- 0
        lss <- lim_bic_loss(W_tmp, X, is_continuous)
        if (lss < aa && lim_h(W_tmp, d)$h < h_tol) {
          W_min_lss <- W_tmp
          aa <- lss
        }
      }
    }
  }

  # (c) edge addition: try adding each currently-absent (undirected) edge
  if (d > 2 && n_edges < (d * (d - 1) / 2)) {
    W0 <- W_min_lss
    W_edges <- W0 + t(W0) + diag(d)
    I_add <- which(W_edges == 0, arr.ind = TRUE)
    for (k in seq_len(nrow(I_add))) {
      W_tmp <- W0
      W_tmp[I_add[k, 1], I_add[k, 2]] <- 1
      lss <- lim_bic_loss(W_tmp, X, is_continuous)
      if (lss < aa && lim_h(W_tmp, d)$h < h_tol) {
        W_min_lss <- W_tmp
        aa <- lss
      }
    }
    if (!identical(W_est != 0, W_min_lss != 0)) {
      W0b <- W_est
      W_edges2 <- W0b + t(W0b) + diag(d)
      I_add2 <- which(W_edges2 == 0, arr.ind = TRUE)
      for (k in seq_len(nrow(I_add2))) {
        W_tmp <- W0b
        W_tmp[I_add2[k, 1], I_add2[k, 2]] <- 1
        lss <- lim_bic_loss(W_tmp, X, is_continuous)
        if (lss < aa && lim_h(W_tmp, d)$h < h_tol) {
          W_min_lss <- W_tmp
          aa <- lss
        }
      }
    }
  }

  W_min_lss
}


#' Topological order (Kahn's algorithm) from a lingamr-convention adjacency
#' matrix B (`B[i, j]` = j -> i, i.e. row i's nonzero columns are its parents)
#' @keywords internal
lim_topological_order <- function(B) {
  d <- ncol(B)
  adj <- B != 0  # adj[i, j] TRUE means j -> i
  indeg <- rowSums(adj)
  removed <- rep(FALSE, d)
  topo_order <- integer(0)

  for (step in seq_len(d)) {
    candidates <- which(!removed & indeg == 0)
    if (length(candidates) == 0) {
      warning(
        "Failed to determine a valid causal order: the estimated adjacency matrix is not acyclic.",
        call. = FALSE
      )
      return(rep(NA_integer_, d))
    }
    i <- candidates[1]
    topo_order <- c(topo_order, i)
    removed[i] <- TRUE
    children <- which(adj[, i] & !removed)
    indeg[children] <- indeg[children] - 1L
  }

  topo_order
}


#' LiM: LiNGAM for Mixed Data
#'
#' Estimates a causal structure from data containing a mixture of continuous
#' and binary (0/1) discrete variables, following Zeng et al. (2022). The
#' method combines a NOTEARS-style continuous optimization (global phase)
#' with a combinatorial local search over edge directions, pruning, and
#' edge addition.
#'
#' @param X Numeric matrix (n_samples x n_features) or data frame
#' @param is_continuous Logical vector of length `ncol(X)`. `TRUE` marks a
#'   continuous variable, `FALSE` marks a discrete (binary 0/1) variable.
#' @param lambda1 L1 penalty parameter (default: 0.1)
#' @param max_iter Maximum number of dual ascent (outer loop) steps (default: 150)
#' @param h_tol Tolerance for the acyclicity constraint h(W) (default: 1e-8)
#' @param rho_max Maximum value of the augmented-Lagrangian penalty rho (default: 1e16)
#' @param w_threshold Edges with |weight| below this value are dropped after
#'   the global optimization phase (default: 0.1)
#' @param only_global If `TRUE`, skip the combinatorial local search phase and
#'   return the thresholded global-optimization result directly (default: FALSE)
#' @return A `LiMResult` object (list) containing the following elements:
#' * `adjacency_matrix`: adjacency matrix B (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i)**, i.e. the same convention as [lingam_direct()].
#'   Zero elements indicate no causal relationship.
#' * `causal_order`: estimated causal order (integer vector of 1-based
#'   indices), derived from `adjacency_matrix` via a topological sort. `NA`
#'   (with a warning) if the estimated matrix is not acyclic.
#' * `is_continuous`: the input `is_continuous` vector, stored for reference.
#'
#' @details
#' Only binary (0/1) discrete variables are supported; count/Poisson-type
#' discrete variables (the Python source's `loss_type = "poisson"` path) are
#' not implemented.
#'
#' The Python implementation's `adjacency_matrix_` uses the opposite
#' convention (`W[i, j]` = i -> j). This R implementation transposes the
#' internal result so that `adjacency_matrix` follows the lingamr convention
#' (`B[i, j]` = j -> i), consistent with [lingam_direct()].
#'
#' In the local search phase, edges that are reversed or newly added are
#' assigned a weight of exactly 1 rather than a re-estimated coefficient
#' (this matches the original Python implementation). As a result, non-zero
#' entries of `adjacency_matrix` are a mix of global-phase estimated
#' coefficients and local-phase placeholder weights of 1.
#'
#' The local search's BIC score for discrete variables with parents uses R's
#' `glm(family = binomial())`, an unregularized maximum-likelihood fit. The
#' Python implementation uses scikit-learn's `LogisticRegression`, which
#' applies L2 regularization by default; consequently, numeric results will
#' not exactly match the Python implementation.
#'
#' @references
#' Zeng Y, Shimizu S, Matsui H, Sun F. Causal discovery for linear mixed
#' data. In: Proceedings of the First Conference on Causal Learning and
#' Reasoning (CLeaR 2022). PMLR 177, pp. 994-1009, 2022.
#'
#' @export
#' @examples
#' # Reproducibility requires set.seed(), since the optimization starts from
#' # a random initial point.
#' set.seed(1)
#' dat <- generate_lim_sample(n = 300)
#' result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' print(result)
lingam_lim <- function(X,
                        is_continuous,
                        lambda1 = 0.1,
                        max_iter = 150L,
                        h_tol = 1e-8,
                        rho_max = 1e16,
                        w_threshold = 0.1,
                        only_global = FALSE) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 2) stop("X must have at least 2 observations (rows).", call. = FALSE)

  d <- ncol(X)

  if (!is.logical(is_continuous) || anyNA(is_continuous)) {
    stop("is_continuous must be a logical vector without NA.", call. = FALSE)
  }
  if (length(is_continuous) != d) {
    stop("is_continuous must have length equal to ncol(X).", call. = FALSE)
  }
  discrete_cols <- which(!is_continuous)
  for (j in discrete_cols) {
    uj <- unique(X[, j])
    # allow floating-point noise (e.g. 1 - 1e-15) around the two valid levels
    if (!all(abs(uj - round(uj)) < 1e-8 & round(uj) %in% c(0, 1))) {
      stop(sprintf(
        "Discrete columns (is_continuous = FALSE) must be binary (0/1); column %d is not.", j
      ), call. = FALSE)
    }
  }

  if (!is.numeric(lambda1) || length(lambda1) != 1 || lambda1 < 0) {
    stop("lambda1 must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) {
    stop("max_iter must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(h_tol) || length(h_tol) != 1 || h_tol <= 0) {
    stop("h_tol must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(rho_max) || length(rho_max) != 1 || rho_max <= 0) {
    stop("rho_max must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(w_threshold) || length(w_threshold) != 1 || w_threshold < 0) {
    stop("w_threshold must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.logical(only_global) || length(only_global) != 1 || is.na(only_global)) {
    stop("only_global must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  max_iter <- as.integer(max_iter)
  con <- as.numeric(is_continuous)

  W_est <- lim_global_optimize(X, con, lambda1, max_iter, h_tol, rho_max, w_threshold)

  W_min_lss <- if (only_global) {
    W_est
  } else {
    lim_local_search(W_est, X, con, is_continuous, d, h_tol)
  }

  B <- t(W_min_lss)
  var_names <- if (!is.null(col_names)) col_names else paste0("x", seq_len(d) - 1L)
  dimnames(B) <- list(var_names, var_names)

  causal_order <- lim_topological_order(B)

  result <- list(
    adjacency_matrix = B,
    causal_order = causal_order,
    is_continuous = is_continuous
  )
  class(result) <- "LiMResult"
  result
}


#' Print method for LiMResult
#'
#' @param x LiMResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' set.seed(1)
#' dat <- generate_lim_sample(n = 300)
#' result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' print(result)
print.LiMResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  var_names <- colnames(x$adjacency_matrix)
  order_labels <- if (anyNA(x$causal_order)) {
    "NA (cyclic result)"
  } else if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  var_types <- ifelse(x$is_continuous, "continuous", "discrete")

  cat("LiM Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Variable types: %s\n", paste(var_types, collapse = ", ")))
  cat(sprintf("  Causal order: %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nAdjacency matrix (row = to, col = from):\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}
