# =============================================================================
# Direct LiNGAM - R Implementation
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
# The internal implementation is split across the following files:
#   R/search_causal_order.r : causal order search (pwling / kernel) and prior knowledge handling
#   R/fit_regression.r      : adjacency matrix estimation and regression backends
#   R/paths.r               : DAG path enumeration and path effects
# =============================================================================


#' Direct LiNGAM
#'
#' @param X Numeric matrix (n_samples x n_features), data frame or matrix
#' @param prior_knowledge Prior knowledge matrix (n_features x n_features) or NULL.
#'   0: no directed path from x_i to x_j
#'   1: directed path from x_i to x_j
#'  -1: unknown
#' @param apply_prior_knowledge_softly Whether to apply prior knowledge softly (logical)
#' @param measure Independence evaluation measure ("pwling" or "kernel")
#' @param reg_method Regression method for adjacency matrix estimation.
#' "ols": ordinary least squares,
#' "lasso": LASSO regression,
#' "adaptive_lasso": adaptive LASSO regression (default),
#' "ridge": Ridge regression (robust to multicollinearity; does not perform sparse estimation).
#' @param init_method Method for estimating the initial weights of adaptive LASSO regression.
#' "ols": ordinary least squares (default),
#' "ridge": Ridge regression.
#' Ridge regression is recommended when multicollinearity is suspected.
#' @param lambda LASSO penalty (lambda) selection.
#' "lambda.min" : minimum CV prediction error, prioritizes prediction accuracy.
#' "lambda.1se" : CV 1SE rule, robust and less prone to overfitting.
#' "AIC": minimum AIC. Fast.
#' "BIC": minimum BIC. Fast, sparsest. Default.
#' "oracle" : adaptive LASSO regression only. Selects a lambda that guarantees the oracle property. Fast.
#' @return A `LingamResult` object (list) containing the following elements:
#' * `adjacency_matrix`: adjacency matrix B (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to variable i (j -> i).**
#'   Zero elements indicate no causal relationship.
#' * `causal_order`: estimated causal order (integer vector of 1-based indices).
#'   Earlier elements are more upstream (closer to exogenous variables).
#' @importFrom stats sd lm.fit cov median quantile
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # OLS (no additional packages required)
#' result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
#' round(result$adjacency_matrix, 3)
#'
#' \donttest{
#' # LASSO (requires glmnet)
#' result_lasso <- lingam_direct(LiNGAM_sample_1000$data)
#' round(result_lasso$adjacency_matrix, 3)
#' }
lingam_direct <- function(X,
                          prior_knowledge = NULL,
                          apply_prior_knowledge_softly = FALSE,
                          measure = "pwling",
                          reg_method = "adaptive_lasso",
                          lambda = "BIC",
                          init_method = "ols") {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 2) stop("X must have at least 2 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names

  measure <- match.arg(measure, c("pwling", "kernel"))
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  if (reg_method == "ridge" && lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }

  if (!is.logical(apply_prior_knowledge_softly) || length(apply_prior_knowledge_softly) != 1) {
    stop("apply_prior_knowledge_softly must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  n_samples <- nrow(X)
  n_features <- ncol(X)
  # --- Prior knowledge preprocessing ---
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
  U <- seq_len(n_features)
  K <- integer(0)
  X_ <- X
  if (measure == "kernel") {
    X_ <- apply(X_, 2, function(col) {
      pop_sd <- sqrt(mean((col - mean(col))^2))
      (col - mean(col)) / pop_sd
    })
  }
  # --- Causal order search ---
  for (iter in seq_len(n_features)) {
    cand <- search_candidate(U, Aknw, apply_prior_knowledge_softly, partial_orders)
    m <- if (measure == "kernel") {
      search_causal_order_kernel(X_, U, cand$Uc, cand$Vj)
    } else {
      search_causal_order_pwling(X_, U, cand$Uc, cand$Vj)
    }
    for (i in U) {
      if (i != m) X_[, i] <- residual_vec(X_[, i], X_[, m])
    }
    K <- c(K, m)
    U <- setdiff(U, m)
    if (!is.null(Aknw) && !apply_prior_knowledge_softly && !is.null(partial_orders)) {
      if (nrow(partial_orders) > 0) {
        partial_orders <- partial_orders[partial_orders[, 1] != m, , drop = FALSE]
      }
    }
  }
  # --- Adjacency matrix estimation (regression method is selectable) ---
  B <- estimate_adjacency_matrix(X, K, Aknw,
    method = reg_method,
    lambda = lambda,
    init_method = init_method
  )
  colnames(B) <- rownames(B) <- colnames(X)
  result <- list(adjacency_matrix = B, causal_order = K)
  class(result) <- "LingamResult"
  result
}


#' Print method for LingamResult
#'
#' @param x LingamResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @export
print.LingamResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  var_names <- colnames(x$adjacency_matrix)
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("Direct LiNGAM Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Causal order: %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nAdjacency matrix (row = to, col = from):\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}


# =============================================================================
# Shared internal utilities
# =============================================================================

#' Validate the return value of lingam_direct()
#' @keywords internal
validate_lingam_result <- function(x) {
  if (!inherits(x, "LingamResult")) {
    stop("lingam_result must be the return value of lingam_direct().", call. = FALSE)
  }
}

#' Population standard deviation (divided by n)
#' @keywords internal
sd_pop <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

#' Get variable names, falling back to x0, x1, ... when colnames is NULL
#' @keywords internal
get_var_names <- function(x) {
  nm <- colnames(x)
  if (is.null(nm)) paste0("x", seq_len(ncol(x)) - 1L) else nm
}
