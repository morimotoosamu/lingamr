# =============================================================================
# MultiGroup Direct LiNGAM - R Implementation
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
# Joint estimation of Direct LiNGAM across multiple datasets ("groups") that
# share a common causal order (Shimizu 2012). Reuses the existing internal
# helpers from R/search_causal_order.r and R/fit_regression.r; only the
# causal-order search objective is new (search_causal_order_pwling_multi()
# below), which is the sample-size-weighted sum of the single-group pwling
# objective across groups.
# =============================================================================


#' Multi-Group Direct LiNGAM
#'
#' Jointly estimates a Direct LiNGAM model from multiple datasets ("groups")
#' that are assumed to share a common causal order but may have different
#' structural coefficients (Shimizu 2012).
#'
#' @param X_list A list of numeric matrices or data frames (length >= 2), one
#'   per group. Each element must have `n_d` rows (sample size may differ by
#'   group) and the same number of columns `p` across all groups.
#' @param prior_knowledge Prior knowledge matrix (n_features x n_features) or
#'   NULL. Applied identically to every group. Same encoding as
#'   [lingam_direct()]:
#'   0: no directed path from x_i to x_j
#'   1: directed path from x_i to x_j
#'  -1: unknown
#' @param apply_prior_knowledge_softly Whether to apply prior knowledge softly (logical)
#' @param reg_method Regression method for adjacency matrix estimation.
#' "ols": ordinary least squares,
#' "lasso": LASSO regression,
#' "adaptive_lasso": adaptive LASSO regression (default),
#' "ridge": Ridge regression (robust to multicollinearity; does not perform sparse estimation).
#' @param lambda LASSO penalty (lambda) selection. Same choices as [lingam_direct()].
#' @param init_method Method for estimating the initial weights of adaptive
#'   LASSO regression ("ols" (default) or "ridge").
#' @return A `MultiGroupLingamResult` object (list) containing:
#' * `adjacency_matrices`: a named list of adjacency matrices, one per group
#'   (name = group name). Each follows the same convention as
#'   [lingam_direct()]: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i).
#' * `causal_order`: the causal order shared by all groups (integer vector of
#'   1-based indices).
#' @details
#' Unlike [lingam_direct()], this function has no `measure` argument: the
#' multi-group causal-order search only supports the pwling (pairwise
#' likelihood / entropy approximation) objective, matching the upstream
#' Python `MultiGroupDirectLiNGAM`, which does not offer a kernel-based
#' multi-group search.
#'
#' For downstream analysis of a single group (total causal effects,
#' independence tests of residuals, plotting), extract that group as a plain
#' `LingamResult` with [get_group_result()] and pass it to the existing
#' single-group functions ([estimate_total_effect()],
#' [estimate_all_total_effects()], [get_error_independence_p_values()],
#' [plot_adjacency()], `autoplot()`, `tidy()`); no multi-group-specific
#' wrappers are provided for these.
#' @references
#' S. Shimizu. Joint estimation of linear non-Gaussian acyclic models.
#' Neurocomputing, 81: 104-107, 2012.
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' res <- lingam_multi_group(mg$data_list, reg_method = "ols")
#' print(res)
#'
#' # Analyze group 1 with the existing single-group tooling
#' g1 <- get_group_result(res, 1)
#' estimate_all_total_effects(mg$data_list[[1]], g1, method = "ols")
lingam_multi_group <- function(X_list,
                               prior_knowledge = NULL,
                               apply_prior_knowledge_softly = FALSE,
                               reg_method = "adaptive_lasso",
                               lambda = "BIC",
                               init_method = "ols") {
  if (!is.list(X_list) || is.data.frame(X_list)) {
    stop("X_list must be a list of numeric matrices or data frames, one per group.", call. = FALSE)
  }
  if (length(X_list) < 2) {
    stop("X_list must contain at least two items (groups).", call. = FALSE)
  }

  group_names <- names(X_list)
  if (is.null(group_names)) group_names <- character(length(X_list))
  missing_name <- !nzchar(group_names)
  if (any(missing_name)) {
    group_names[missing_name] <- paste0("group", seq_along(X_list))[missing_name]
  }

  # --- Per-group validation and coercion (mirrors lingam_direct()'s checks) ---
  n_features_list <- integer(length(X_list))
  X_mats <- vector("list", length(X_list))
  colname_list <- vector("list", length(X_list))
  for (d in seq_along(X_list)) {
    Xd <- X_list[[d]]
    colname_list[[d]] <- if (is.data.frame(Xd)) names(Xd) else colnames(Xd)
    Xd <- as.matrix(Xd)
    if (!is.numeric(Xd)) {
      stop(sprintf("X_list[[%d]] must be a numeric matrix or data frame.", d), call. = FALSE)
    }
    if (anyNA(Xd)) {
      stop(sprintf("X_list[[%d]] must not contain missing values (NA).", d), call. = FALSE)
    }
    if (nrow(Xd) < 2) {
      stop(sprintf("X_list[[%d]] must have at least 2 observations (rows).", d), call. = FALSE)
    }
    X_mats[[d]] <- Xd
    n_features_list[d] <- ncol(Xd)
  }
  if (length(unique(n_features_list)) > 1) {
    stop("All items in X_list must have the same number of columns (variables).", call. = FALSE)
  }
  n_features <- n_features_list[1]
  if (n_features < 2) stop("X_list items must have at least 2 variables (columns).", call. = FALSE)

  col_names <- colname_list[[1]]
  if (!is.null(col_names)) {
    for (d in seq_along(colname_list)[-1]) {
      cn_d <- colname_list[[d]]
      if (!is.null(cn_d) && !identical(cn_d, col_names)) {
        warning(sprintf(
          "Column names of X_list[[%d]] differ from X_list[[1]]; using X_list[[1]]'s names for all groups.",
          d
        ), call. = FALSE)
      }
    }
  }
  for (d in seq_along(X_mats)) colnames(X_mats[[d]]) <- col_names

  if (!is.logical(apply_prior_knowledge_softly) || length(apply_prior_knowledge_softly) != 1) {
    stop("apply_prior_knowledge_softly must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))
  if (reg_method %in% c("lasso", "ridge") && lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }

  # --- Prior knowledge preprocessing (shared across all groups) ---
  Aknw <- NULL
  partial_orders <- NULL
  if (!is.null(prior_knowledge)) {
    Aknw <- as.matrix(prior_knowledge)
    if (!all(dim(Aknw) == c(n_features, n_features))) {
      stop("The shape of prior knowledge must be (n_features, n_features)")
    }
    Aknw[Aknw < 0] <- NA
    if (!apply_prior_knowledge_softly) {
      partial_orders <- extract_partial_orders(Aknw)
    }
  }

  # --- Causal order search (common across groups, weighted by sample size) ---
  U <- seq_len(n_features)
  K <- integer(0)
  X_work <- X_mats # working copies, residualized in place per iteration
  for (iter in seq_len(n_features)) {
    cand <- search_candidate(U, Aknw, apply_prior_knowledge_softly, partial_orders)
    m <- search_causal_order_pwling_multi(X_work, U, cand$Uc, cand$Vj)
    for (d in seq_along(X_work)) {
      for (i in U) {
        if (i != m) X_work[[d]][, i] <- residual_vec(X_work[[d]][, i], X_work[[d]][, m])
      }
    }
    K <- c(K, m)
    U <- setdiff(U, m)
    if (!is.null(Aknw) && !apply_prior_knowledge_softly && !is.null(partial_orders)) {
      if (nrow(partial_orders) > 0) {
        partial_orders <- partial_orders[partial_orders[, 1] != m, , drop = FALSE]
      }
    }
  }

  # --- Adjacency matrix per group, estimated from the ORIGINAL (non-residualized) data ---
  B_list <- lapply(X_mats, function(Xd) {
    B <- estimate_adjacency_matrix(Xd, K, Aknw,
      method = reg_method,
      lambda = lambda,
      init_method = init_method
    )
    colnames(B) <- rownames(B) <- col_names
    B
  })
  names(B_list) <- group_names

  result <- list(adjacency_matrices = B_list, causal_order = K)
  class(result) <- "MultiGroupLingamResult"
  result
}


#' Print method for MultiGroupLingamResult
#'
#' @param x MultiGroupLingamResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' res <- lingam_multi_group(mg$data_list, reg_method = "ols")
#' print(res)
print.MultiGroupLingamResult <- function(x, digits = 3, ...) {
  n_groups <- length(x$adjacency_matrices)
  n <- length(x$causal_order)
  var_names <- colnames(x$adjacency_matrices[[1]])
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("Multi-Group Direct LiNGAM Result\n")
  cat(sprintf("  Groups      : %d (%s)\n", n_groups, paste(names(x$adjacency_matrices), collapse = ", ")))
  cat(sprintf("  Variables   : %d\n", n))
  cat(sprintf("  Causal order (common): %s\n", paste(order_labels, collapse = " -> ")))
  for (g in names(x$adjacency_matrices)) {
    cat(sprintf("\n[%s] Adjacency matrix (row = to, col = from):\n", g))
    print(round(x$adjacency_matrices[[g]], digits = digits))
  }
  invisible(x)
}


#' Extract a single group's result from a MultiGroupLingamResult
#'
#' Returns the adjacency matrix and (shared) causal order of one group as a
#' plain `LingamResult`, so that the existing single-group functions
#' ([estimate_total_effect()], [estimate_all_total_effects()],
#' [get_error_independence_p_values()], [plot_adjacency()], `autoplot()`,
#' `tidy()`) can be applied to it directly.
#'
#' @param x A `MultiGroupLingamResult`, as returned by [lingam_multi_group()].
#' @param group Group name (character) or 1-based group index (integer).
#' @return A `LingamResult` object (list) with `adjacency_matrix` and
#'   `causal_order`, identical in shape to the return value of
#'   [lingam_direct()].
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' res <- lingam_multi_group(mg$data_list, reg_method = "ols")
#' g1 <- get_group_result(res, 1)
#' class(g1)
get_group_result <- function(x, group) {
  if (!inherits(x, "MultiGroupLingamResult")) {
    stop("x must be the return value of lingam_multi_group().", call. = FALSE)
  }
  groups <- names(x$adjacency_matrices)
  if (is.character(group)) {
    if (length(group) != 1 || !(group %in% groups)) {
      stop(sprintf(
        "Group '%s' not found. Available groups: %s",
        paste(group, collapse = ", "), paste(groups, collapse = ", ")
      ), call. = FALSE)
    }
  } else if (is.numeric(group)) {
    group <- as.integer(group)
    if (length(group) != 1 || is.na(group) || group < 1 || group > length(groups)) {
      stop(sprintf("'group' must be an integer between 1 and %d.", length(groups)), call. = FALSE)
    }
  } else {
    stop("'group' must be a group name (character) or a 1-based index (integer).", call. = FALSE)
  }

  result <- list(
    adjacency_matrix = x$adjacency_matrices[[group]],
    causal_order = x$causal_order
  )
  class(result) <- "LingamResult"
  result
}


# =============================================================================
# Internal: multi-group causal order search
# =============================================================================

#' Causal order search via pwling, jointly across multiple groups
#'
#' Sample-size-weighted sum of the single-group pwling objective
#' ([search_causal_order_pwling()]) across all groups in `X_list`, following
#' the joint estimation objective of Shimizu (2012). Reuses the same
#' standardize-once / correlation-matrix / antisymmetry optimizations as
#' [search_causal_order_pwling()], applied independently within each group
#' before the weighted sum.
#'
#' @param X_list List of per-group data matrices (residualized so far),
#'   one per group. All must have the same number of columns.
#' @param U Indices of all currently undetermined variables (shared across groups)
#' @param Uc Indices of candidate variables (shared across groups)
#' @param Vj Variable set based on prior knowledge (shared across groups)
#' @return Index of the selected variable
#' @keywords internal
search_causal_order_pwling_multi <- function(X_list, U, Uc, Vj) {
  if (length(Uc) == 1) return(Uc[1])

  p <- ncol(X_list[[1]])
  n_list <- vapply(X_list, nrow, integer(1))
  total_size <- sum(n_list)

  pos <- integer(p)
  pos[U] <- seq_along(U)
  in_Uc <- logical(p)
  in_Uc[Uc] <- TRUE
  in_Vj <- logical(p)
  in_Vj[Vj] <- TRUE

  M_total <- numeric(p)
  for (d in seq_along(X_list)) {
    Xd <- X_list[[d]]
    n <- n_list[d]

    # --- Standardize all columns at once and precompute each column's entropy ---
    X_std <- matrix(0, nrow = n, ncol = p)
    H <- numeric(p)
    for (k in U) {
      xk <- Xd[, k]
      xk <- xk - sum(xk) / n
      X_std[, k] <- xk / sqrt(sum(xk * xk) / n)
      H[k] <- entropy_approx(X_std[, k])
    }
    R <- crossprod(X_std[, U, drop = FALSE]) / n

    M_acc <- numeric(p)
    for (i in Uc) {
      xi_std <- X_std[, i]
      for (j in U) {
        if (i == j) next
        # The pairwise mutual-information difference is antisymmetric, so for
        # pairs where both are candidates, compute once on the i < j side and
        # add to both.
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
    M_total <- M_total + (n / total_size) * M_acc
  }

  return(Uc[which.max(-M_total[Uc])])
}
