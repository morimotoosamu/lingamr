# =============================================================================
# Bottom-Up ParceLiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/bottom_up_parce_lingam.py)
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
#   Tashiro, T., Shimizu, S., Hyvarinen, A., and Washio, T. (2014).
#   ParceLiNGAM: a causal ordering method robust against latent confounders.
#   Neural Computation, 26(1), 57-83.
#
# The internal implementation reuses:
#   R/hsic.r            : HSIC gamma-approximation independence test
#   R/f_correlation.r   : F-correlation independence measure
#   R/lingam_high_dim.r : pinv() (Moore-Penrose pseudo-inverse)
#   R/fit_regression.r  : per-target regression backends (fit_ols etc.)
#   R/search_causal_order.r : extract_partial_orders()
#   R/paths.r            : find_all_paths() / calculate_total_effect() (used by bootstrap)
# =============================================================================


#' Bottom-Up ParceLiNGAM
#'
#' A causal ordering method robust against latent confounders. Unlike
#' [lingam_direct()], which always returns a full causal order, this
#' algorithm searches from the sink (most downstream) side and stops as soon
#' as an independence test is rejected. Variables it could not order are
#' returned together as a single "unresolved block" (suspected to share a
#' latent confounder), and the corresponding entries of the adjacency matrix
#' are set to `NA` rather than estimated.
#'
#' @param X Numeric matrix (n_samples x n_features), data frame or matrix
#' @param alpha Significance level for the independence test. `alpha = 0`
#'   disables rejection entirely, so a full causal order is always returned
#'   (equivalent in spirit to [lingam_direct()], but using the ParceLiNGAM
#'   sink-search direction and regression). Must be non-negative.
#' @param prior_knowledge Prior knowledge matrix (n_features x n_features) or NULL.
#'   0: no directed path from x_i to x_j
#'   1: directed path from x_i to x_j
#'  -1: unknown
#' @param independence Independence measure used for the ordering search:
#'   "hsic" (default) uses the HSIC gamma-approximation test combined across
#'   explanatory variables via Fisher's method; "fcorr" uses the F-correlation
#'   (kernel canonical correlation) and rejects based on `ind_corr` instead
#'   of a p-value.
#' @param ind_corr Threshold on the F-correlation value used only when
#'   `independence = "fcorr"`: a candidate is rejected once the largest
#'   F-correlation against any explanatory variable is at or above this
#'   value. Must be non-negative. Ignored when `independence = "hsic"`.
#' @param reg_method Regression method for adjacency matrix estimation.
#' "ols": ordinary least squares,
#' "lasso": LASSO regression,
#' "adaptive_lasso": adaptive LASSO regression (default, matches the
#'  upstream Python implementation's `predict_adaptive_lasso`),
#' "ridge": Ridge regression.
#' @param lambda LASSO penalty (lambda) selection. Same options as
#'   [lingam_direct()]: "lambda.min", "lambda.1se", "AIC", "BIC" (default),
#'   "oracle" (adaptive LASSO only).
#' @param init_method Method for estimating the initial weights of adaptive
#'   LASSO regression ("ols" (default) or "ridge").
#' @return A `ParceLingamResult` object (list) containing:
#' * `adjacency_matrix`: adjacency matrix B (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i)**, same as [lingam_direct()]. Entries between two
#'   variables that ended up in the same unresolved block are `NA`.
#' * `causal_order`: a list of integer vectors. Elements of length 1 are
#'   variables with a fully resolved position; an element of length > 1 (at
#'   most one, always first) is the unresolved block. Earlier elements are
#'   more upstream.
#' * `p_values`: independence-test p-values (or F-correlation values, for
#'   `independence = "fcorr"`) for each step that successfully placed a
#'   variable, in the order variables were placed (diagnostic only).
#' * `independence`: the independence measure used.
#' @details
#' Because HSIC forms full n x n Gram matrices, it is O(n^2) per test; avoid
#' very large `n` (beyond a few thousand) with `independence = "hsic"`.
#'
#' `independence = "fcorr"` rejects based on the raw F-correlation value
#' (`ind_corr`), not a p-value, so it is not directly comparable to `alpha`.
#'
#' [get_error_independence_p_values_parce()] uses the HSIC test rather than
#' the correlation-based test used by [get_error_independence_p_values()]
#' for `LingamResult` objects.
#'
#' [lingam_parce_bootstrap()] treats `NA` (unresolved) edges as absent when
#' aggregating, and does not support [get_causal_order_stability()] (see its
#' documentation for details). This function does not expose a `regressor`
#' or `bw_method` argument, unlike the upstream Python implementation.
#' @references
#' Tashiro, T., Shimizu, S., Hyvarinen, A., and Washio, T. (2014).
#' ParceLiNGAM: a causal ordering method robust against latent confounders.
#' Neural Computation, 26(1), 57-83.
#' @export
#' @examples
#' confounded <- generate_parce_sample(n = 500, seed = 1)
#'
#' result <- lingam_parce(confounded$data, reg_method = "ols")
#' print(result)
#'
#' # The variable pair sharing the latent confounder is left unresolved (NA)
#' result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#'
#' # Total effect estimation warns and returns NA for confounded variables
#' estimate_total_effect_parce(confounded$data, result,
#'   from_index = confounded$confounded_pair[1], to_index = 1
#' )
lingam_parce <- function(X,
                         alpha = 0.1,
                         prior_knowledge = NULL,
                         independence = "hsic",
                         ind_corr = 0.5,
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

  independence <- match.arg(independence, c("hsic", "fcorr"))
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  if (reg_method %in% c("lasso", "ridge") && lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0) {
    stop("alpha must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(ind_corr) || length(ind_corr) != 1 || is.na(ind_corr) || ind_corr < 0) {
    stop("ind_corr must be a non-negative numeric scalar.", call. = FALSE)
  }

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
    partial_orders <- extract_partial_orders(Aknw)
  }

  # --- Centering (original 124 line) ---
  Xc <- scale(X, center = TRUE, scale = FALSE)
  attr(Xc, "scaled:center") <- NULL
  class(Xc) <- "matrix"

  # --- Bonferroni-corrected threshold p-value (fixed at p, not shrinking U) ---
  thresh_p <- alpha / (n_features - 1)

  search <- parce_search_causal_order(
    Xc, seq_len(n_features), partial_orders, independence, thresh_p, ind_corr
  )

  causal_order <- c(
    if (length(search$U_res) >= 2) list(sort(search$U_res)) else list(),
    lapply(search$K_bttm, function(k) k)
  )

  B <- estimate_adjacency_matrix_parce(
    X, causal_order, Aknw,
    method = reg_method, lambda = lambda, init_method = init_method
  )

  result <- list(
    adjacency_matrix = B,
    causal_order = causal_order,
    p_values = search$p_bttm,
    independence = independence
  )
  class(result) <- "ParceLingamResult"
  result
}


#' Print method for ParceLingamResult
#'
#' @param x ParceLingamResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' confounded <- generate_parce_sample(n = 300, seed = 42)
#' result <- lingam_parce(confounded$data, reg_method = "ols")
#' print(result)
print.ParceLingamResult <- function(x, digits = 3, ...) {
  var_names <- colnames(x$adjacency_matrix)
  label_for <- function(idx) {
    if (!is.null(var_names)) var_names[idx] else paste0("x", idx - 1L)
  }
  order_labels <- vapply(x$causal_order, function(blk) {
    if (length(blk) > 1) {
      paste0("(", paste(label_for(blk), collapse = ", "), ")")
    } else {
      label_for(blk)
    }
  }, character(1))

  cat("Bottom-Up ParceLiNGAM Result\n")
  cat(sprintf("  Variables : %d\n", ncol(x$adjacency_matrix)))
  cat(sprintf("  Independence measure: %s\n", x$independence))
  cat(sprintf("  Causal order: %s\n", paste(order_labels, collapse = " -> ")))
  cat("  (NA entries in the adjacency matrix = unresolved order / suspected latent confounding)\n")
  cat("\nAdjacency matrix (row = to, col = from):\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}


# =============================================================================
# Internal: bottom-up causal order search
# =============================================================================


#' Candidate sink variables under prior knowledge (ParceLiNGAM direction)
#'
#' Unlike [search_candidate()] (used by top-down DirectLiNGAM), ParceLiNGAM
#' searches from the sink side, so the filter is simply "exclude variables
#' that appear as the 'from' side of a known partial order" (they cannot be
#' a sink because something is known to cause a variable through them...
#' more precisely: a variable known to *cause* another remaining variable
#' cannot itself be the next sink).
#'
#' @param U set of currently undetermined variables
#' @param partial_orders matrix of (from, to) pairs from [extract_partial_orders()]
#' @return candidate index vector (subset of U, or U itself if no candidates remain)
#' @keywords internal
parce_search_candidate <- function(U, partial_orders) {
  if (is.null(partial_orders) || nrow(partial_orders) == 0) {
    return(U)
  }
  Uc <- setdiff(U, partial_orders[, 1])
  if (length(Uc) == 0) Uc <- U
  Uc
}


#' Residual of `X[, j]` regressed on `X[, xi_index]` via the pseudo-inverse
#' of the covariance matrix
#'
#' @param X data matrix
#' @param xi_index explanatory-variable indices (may be empty)
#' @param j target variable index
#' @param Cov precomputed `stats::cov(X)`, since `X` is invariant across all
#'   calls within a single [parce_search_causal_order()] search
#' @return residual vector
#' @keywords internal
parce_residual <- function(X, xi_index, j, Cov) {
  if (length(xi_index) == 0) {
    return(X[, j])
  }
  coef <- pinv(Cov[xi_index, xi_index, drop = FALSE]) %*% Cov[xi_index, j]
  X[, j] - as.vector(X[, xi_index, drop = FALSE] %*% coef)
}


#' Evaluate the independence of a residual against a fixed set of predictors
#'
#' Used for the `length(Uc) == 1` special case in [find_exo_vec()], where
#' there is nothing left to compare against.
#'
#' @param X data matrix
#' @param predictors predictor indices (may be empty)
#' @param R residual vector
#' @param independence "hsic" or "fcorr"
#' @return evaluation value (Fisher-combined p-value for hsic, max
#'   F-correlation for fcorr)
#' @keywords internal
parce_eval_independence <- function(X, predictors, R, independence) {
  if (length(predictors) == 0) {
    return(if (independence == "hsic") 1.0 else 0.0)
  }
  if (independence == "hsic") {
    if (length(predictors) == 1) {
      return(hsic_test_gamma(X[, predictors], R)$p)
    }
    fisher_stat <- 0
    for (k in predictors) {
      p_k <- hsic_test_gamma(X[, k], R)$p
      fisher_stat <- fisher_stat + if (p_k <= 0) Inf else -2 * log(p_k)
    }
    return(1 - stats::pchisq(fisher_stat, df = 2 * length(predictors)))
  }
  max(vapply(predictors, function(k) f_correlation(X[, k], R), numeric(1)))
}


#' Find the most sink-like candidate variable
#'
#' Faithful port of `_search_causal_order.find_exo_vec()` (original
#' 227-276 lines). For each candidate `j` in `Uc`, regresses `j` on the
#' other candidates (`setdiff(Uc, j)`, **not** `U`) and evaluates how
#' independent the residual is from those explanatory variables.
#'
#' @param X data matrix
#' @param Uc candidate variable indices
#' @param U all currently undetermined variable indices
#' @param independence "hsic" or "fcorr"
#' @param Cov precomputed `stats::cov(X)`, since `X` is invariant across all
#'   calls within a single [parce_search_causal_order()] search
#' @return list(m = selected variable index, eval = its evaluation value)
#' @keywords internal
find_exo_vec <- function(X, Uc, U, independence, Cov) {
  if (length(Uc) == 1) {
    j <- Uc[1]
    predictors <- setdiff(U, Uc)
    R <- parce_residual(X, predictors, j, Cov)
    ev <- parce_eval_independence(X, predictors, R, independence)
    return(list(m = j, eval = ev))
  }

  best_m <- NA_integer_
  best_eval <- NA_real_
  best_stat <- Inf

  for (j in Uc) {
    xi_index <- setdiff(Uc, j)
    R <- parce_residual(X, xi_index, j, Cov)

    if (independence == "hsic") {
      if (length(xi_index) == 1) {
        ht <- hsic_test_gamma(X[, xi_index], R)
        cand_stat <- ht$stat
        cand_eval <- ht$p
      } else {
        fisher_stat <- 0
        aborted <- FALSE
        for (k in xi_index) {
          p_k <- hsic_test_gamma(X[, k], R)$p
          fisher_stat <- fisher_stat + if (p_k <= 0) Inf else -2 * log(p_k)
          # Early break: statistics only grow as more terms are added, and we
          # are looking for the candidate with the *smallest* statistic, so a
          # partial sum already exceeding the current best can never win.
          if (fisher_stat > best_stat) {
            aborted <- TRUE
            break
          }
        }
        if (aborted) next
        cand_stat <- fisher_stat
        cand_eval <- 1 - stats::pchisq(fisher_stat, df = 2 * length(xi_index))
      }
    } else {
      fmax <- -Inf
      aborted <- FALSE
      for (k in xi_index) {
        fc <- f_correlation(X[, k], R)
        if (fc > fmax) fmax <- fc
        if (fmax > best_stat) {
          aborted <- TRUE
          break
        }
      }
      if (aborted) next
      cand_stat <- fmax
      cand_eval <- fmax
    }

    if (cand_stat < best_stat) {
      best_stat <- cand_stat
      best_eval <- cand_eval
      best_m <- j
    }
  }

  list(m = best_m, eval = best_eval)
}


#' Bottom-up causal order search
#'
#' Faithful port of `_search_causal_order()` (original 188-225 lines).
#' Repeatedly finds the most sink-like remaining variable and appends it to
#' the front of `K_bttm` (bottom-up, so more recently placed variables are
#' more upstream). Stops as soon as a candidate is rejected by the
#' independence test; the remaining undetermined variables are returned as
#' `U_res`.
#'
#' @param X (centered) data matrix
#' @param U all variable indices
#' @param partial_orders matrix of (from, to) pairs, or NULL
#' @param independence "hsic" or "fcorr"
#' @param thresh_p Bonferroni-corrected significance threshold (hsic only)
#' @param ind_corr F-correlation rejection threshold (fcorr only)
#' @return list(K_bttm = integer vector, p_bttm = numeric vector, U_res = integer vector)
#' @keywords internal
parce_search_causal_order <- function(X, U, partial_orders, independence, thresh_p, ind_corr) {
  K_bttm <- integer(0)
  p_bttm <- numeric(0)
  # X does not change during this search, so its covariance matrix is
  # computed once here rather than on every find_exo_vec()/parce_residual() call.
  Cov <- stats::cov(X)

  repeat {
    Uc <- parce_search_candidate(U, partial_orders)
    res <- find_exo_vec(X, Uc, U, independence, Cov)
    m <- res$m
    ev <- res$eval

    rejected <- if (independence == "hsic") ev < thresh_p else ev >= ind_corr

    if (rejected) break

    K_bttm <- c(m, K_bttm)
    p_bttm <- c(ev, p_bttm)
    U <- setdiff(U, m)

    if (!is.null(partial_orders) && nrow(partial_orders) > 0) {
      partial_orders <- partial_orders[partial_orders[, 2] != m, , drop = FALSE]
    }

    if (length(U) <= 1) {
      K_bttm <- c(U, K_bttm)
      p_bttm <- c(0.0, p_bttm)
      U <- integer(0)
      break
    }
  }

  list(K_bttm = K_bttm, p_bttm = p_bttm, U_res = U)
}


# =============================================================================
# Internal: adjacency matrix estimation
# =============================================================================


#' Estimate the adjacency matrix from a ParceLiNGAM causal order
#'
#' `causal_order` is a list whose first element may be an unresolved block
#' (length > 1); all remaining elements are length-1. Block members are
#' never regression targets (their parents cannot be identified), but they
#' are valid predictors for downstream variables. Pairs within the block are
#' set to `NA`.
#'
#' @param X original (uncentered) data
#' @param causal_order list as produced by [lingam_parce()]
#' @param prior_knowledge prior-knowledge matrix (NULL allowed)
#' @param method regression method
#' @param lambda lambda selection
#' @param init_method adaptive LASSO initial-weight method
#' @return adjacency matrix B (n_features x n_features)
#' @keywords internal
estimate_adjacency_matrix_parce <- function(X, causal_order, prior_knowledge,
                                            method, lambda, init_method) {
  n_features <- ncol(X)
  B <- matrix(0, n_features, n_features)

  if (length(causal_order) >= 2) {
    for (i in 2:length(causal_order)) {
      target <- causal_order[[i]][1]
      predictors <- unlist(causal_order[seq_len(i - 1)])

      if (!is.null(prior_knowledge)) {
        keep <- vapply(predictors, function(p) {
          val <- prior_knowledge[target, p]
          is.na(val) || val != 0
        }, logical(1))
        predictors <- predictors[keep]
      }

      if (length(predictors) == 0) next

      y <- X[, target]
      Xp <- X[, predictors, drop = FALSE]

      B[target, predictors] <- fit_coef_by_method(y, Xp, method, lambda, init_method)
    }
  }

  if (length(causal_order) >= 1 && length(causal_order[[1]]) > 1) {
    blk <- causal_order[[1]]
    B[blk, blk] <- NA
    diag(B)[blk] <- 0
  }

  colnames(B) <- rownames(B) <- colnames(X)
  B
}


# =============================================================================
# Total causal effect and error independence (Parce-specific)
# =============================================================================


#' Validate the return value of lingam_parce()
#' @keywords internal
validate_parce_result <- function(x) {
  if (!inherits(x, "ParceLingamResult")) {
    stop("parce_result must be the return value of lingam_parce().", call. = FALSE)
  }
}


#' Find the position (rank) of a variable within a ParceLiNGAM causal order
#'
#' All members of the unresolved block (if any) share the same rank (1).
#'
#' @param causal_order list as produced by [lingam_parce()]
#' @param idx 1-based variable index
#' @return integer rank, or NA if not found
#' @keywords internal
parce_order_rank <- function(causal_order, idx) {
  for (i in seq_along(causal_order)) {
    if (idx %in% causal_order[[i]]) return(i)
  }
  NA_integer_
}


#' Estimate the total causal effect between two variables (ParceLiNGAM)
#'
#' Analogous to [estimate_total_effect()], but for [lingam_parce()] results,
#' which may contain `NA` entries in the adjacency matrix.
#'
#' @param X Original data (matrix or data.frame)
#' @param parce_result Return value of [lingam_parce()]
#' @param from_index Cause variable (1-based index or variable name)
#' @param to_index Effect variable (1-based index or variable name)
#' @param method Regression method ("ols", "lasso", "adaptive_lasso", "ridge"). Default is adaptive_lasso
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle"). Default is BIC
#' @param init_method Method for estimating the initial weights of adaptive LASSO regression ("ols" or "ridge")
#' @return Estimated total causal effect, or `NA` (with a warning) if `from_index`
#'   is part of an unresolved block (its parents cannot be identified).
#' @export
#' @examples
#' confounded <- generate_parce_sample(n = 500, seed = 1)
#' result <- lingam_parce(confounded$data, reg_method = "ols")
#'
#' # A well-identified pair returns a numeric estimate
#' estimate_total_effect_parce(confounded$data, result, from_index = 1, to_index = 5)
estimate_total_effect_parce <- function(X, parce_result, from_index, to_index,
                                        method = "adaptive_lasso", lambda = "BIC",
                                        init_method = "ols") {
  validate_parce_result(parce_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  X <- as.matrix(X)
  col_names <- colnames(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  n_features <- ncol(X)
  B <- parce_result$adjacency_matrix
  if (ncol(B) != n_features) {
    stop(
      "X has ", n_features, " variables but parce_result was estimated from ",
      ncol(B), " variables.",
      call. = FALSE
    )
  }

  resolve_index <- function(idx, arg_name) {
    if (is.character(idx)) {
      if (is.null(col_names)) {
        stop(sprintf("'%s' was specified as a name, but X has no column names.", arg_name))
      }
      pos <- match(idx, col_names)
      if (is.na(pos)) {
        stop(sprintf(
          "Variable '%s' not found. Available: %s",
          idx, paste(col_names, collapse = ", ")
        ))
      }
      return(pos)
    } else if (is.numeric(idx)) {
      idx <- as.integer(idx)
      if (idx < 1 || idx > n_features) {
        stop(sprintf("'%s' must be between 1 and %d.", arg_name, n_features))
      }
      return(idx)
    }
    stop(sprintf("'%s' must be integer or character.", arg_name))
  }

  from_index <- resolve_index(from_index, "from_index")
  to_index <- resolve_index(to_index, "to_index")
  if (from_index == to_index) stop("from_index and to_index must differ.")

  causal_order <- parce_result$causal_order
  from_rank <- parce_order_rank(causal_order, from_index)
  to_rank <- parce_order_rank(causal_order, to_index)

  from_label <- if (!is.null(col_names)) col_names[from_index] else paste0("x", from_index)
  to_label <- if (!is.null(col_names)) col_names[to_index] else paste0("x", to_index)

  if (from_rank > to_rank) {
    warning(sprintf(
      "Causal order of %s (to) is earlier than %s (from). Result may be incorrect.",
      to_label, from_label
    ))
  }

  if (anyNA(B[from_index, ])) {
    warning(sprintf(
      "%s is part of an unresolved causal order (suspected latent confounding); total effect cannot be estimated.",
      from_label
    ))
    return(NA_real_)
  }

  parents <- which(!is.na(B[from_index, ]) & abs(B[from_index, ]) > 0)
  predictors <- unique(c(from_index, parents))
  from_pos <- which(predictors == from_index)

  y <- X[, to_index]
  Xp <- X[, predictors, drop = FALSE]

  coefs <- fit_coef_by_method(y, Xp, method, lambda, init_method)

  coefs[from_pos]
}


#' Compute p-values for the independence of ParceLiNGAM residuals (HSIC-based)
#'
#' Analogous to [get_error_independence_p_values()], but for [lingam_parce()]
#' results. Uses the HSIC gamma-approximation test ([hsic_test_gamma()])
#' rather than a correlation test, and returns `NA` for any pair involving a
#' variable whose row or column in the adjacency matrix contains `NA`
#' (residuals cannot be computed for those variables).
#'
#' @param X Original data (matrix or data.frame)
#' @param parce_result Return value of [lingam_parce()]
#' @return matrix of p-values (n_features x n_features)
#' @export
#' @examples
#' confounded <- generate_parce_sample(n = 500, seed = 1)
#' result <- lingam_parce(confounded$data, reg_method = "ols")
#' round(get_error_independence_p_values_parce(confounded$data, result), 3)
get_error_independence_p_values_parce <- function(X, parce_result) {
  validate_parce_result(parce_result)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  B <- parce_result$adjacency_matrix
  n_features <- ncol(B)
  if (ncol(X) != n_features) {
    stop(
      "X has ", ncol(X), " variables but parce_result was estimated from ",
      n_features, " variables.",
      call. = FALSE
    )
  }

  na_vars <- unique(c(
    which(apply(B, 1, anyNA)),
    which(apply(B, 2, anyNA))
  ))

  B0 <- B
  B0[is.na(B0)] <- 0
  E <- X - X %*% t(B0)

  p_values <- matrix(NA_real_, n_features, n_features)
  pairs <- which(upper.tri(matrix(TRUE, n_features, n_features)), arr.ind = TRUE)
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1]
    j <- pairs[r, 2]
    if (i %in% na_vars || j %in% na_vars) next
    p_val <- hsic_test_gamma(E[, i], E[, j])$p
    p_values[i, j] <- p_val
    p_values[j, i] <- p_val
  }

  colnames(p_values) <- rownames(p_values) <- colnames(X)
  p_values
}
