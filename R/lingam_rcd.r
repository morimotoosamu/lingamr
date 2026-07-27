# =============================================================================
# RCD (Repetitive Causal Discovery) - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/rcd.py)
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
#   Maeda, T. N. and Shimizu, S. (2020). RCD: Repetitive causal discovery of
#   linear non-Gaussian acyclic models with latent confounders. AISTATS 2020,
#   PMLR 108: 735-745.
#
# The internal implementation reuses:
#   R/hsic.r            : hsic_kernel_width() / hsic_gram_matrix() / hsic_test_gamma()
#   R/f_correlation.r   : f_correlation()
#   R/paths.r           : calculate_total_effect() (used by bootstrap)
#   R/get_error_independence_p_values.r : SHAPIRO_MAX_N constant (referenced,
#     not redefined; see rcd_is_non_gaussian_all())
# =============================================================================


# =============================================================================
# Internal: shared helpers
# =============================================================================


#' OLS residual and coefficients (with intercept)
#'
#' @param y response vector
#' @param X_exog explanatory-variable matrix (may have 0 columns)
#' @return list(resid = residual vector, coef = coefficient vector excluding intercept)
#' @keywords internal
rcd_ols_resid_coef <- function(y, X_exog) {
  if (is.null(X_exog) || ncol(X_exog) == 0) {
    return(list(resid = as.vector(y), coef = numeric(0)))
  }
  design <- cbind(1, X_exog)
  fit <- stats::lm.fit(design, y)
  list(resid = as.vector(fit$residuals), coef = fit$coefficients[-1])
}


#' Pearson correlation test used throughout RCD
#'
#' @param a numeric vector
#' @param b numeric vector
#' @param cor_alpha significance level
#' @return TRUE if correlated (p-value < cor_alpha)
#' @keywords internal
rcd_is_correlated <- function(a, b, cor_alpha) {
  stats::cor.test(a, b, method = "pearson")$p.value < cor_alpha
}


#' Independence judgment used throughout RCD (hsic or fcorr)
#'
#' @param a numeric vector
#' @param b numeric vector
#' @param independence "hsic" or "fcorr"
#' @param ind_alpha significance level (hsic only)
#' @param ind_corr rejection threshold (fcorr only)
#' @return TRUE if independent
#' @keywords internal
rcd_is_independent <- function(a, b, independence, ind_alpha, ind_corr) {
  if (independence == "hsic") {
    hsic_test_gamma(a, b)$p > ind_alpha
  } else {
    f_correlation(a, b) < ind_corr
  }
}


#' Non-Gaussianity judgment (Shapiro-Wilk) for a set of columns
#'
#' When `n > SHAPIRO_MAX_N`, `stats::shapiro.test()` cannot be used directly
#' (hard cap at 5000), so the deterministic evenly-spaced subsample from
#' `shapiro_subsample()` (defined in `R/get_error_independence_p_values.r`)
#' is tested instead, matching the behavior of `test_residual_normality()`.
#' The deterministic thinning keeps results reproducible across calls and
#' leaves the caller's RNG stream untouched.
#'
#' @param Y data matrix
#' @param cols column indices to test
#' @param shapiro_alpha significance level
#' @return TRUE if all columns reject normality (p <= shapiro_alpha)
#' @keywords internal
rcd_is_non_gaussian_all <- function(Y, cols, shapiro_alpha) {
  for (j in cols) {
    p <- stats::shapiro.test(shapiro_subsample(Y[, j]))$p.value
    if (p > shapiro_alpha) return(FALSE)
  }
  TRUE
}


# =============================================================================
# Internal: MLHSICR regression (original 186-224 lines)
# =============================================================================


#' HSIC-sum-minimizing ("MLHSICR") regression
#'
#' Regresses `Y[, xi]` on `Y[, xj_list]` (no intercept) by minimizing the sum
#' of the empirical HSIC statistics between the residual and each explanatory
#' variable, instead of ordinary least squares. Used as a fallback when the
#' OLS residual is not independent of the explanatory variables.
#'
#' The kernel width used to build the residual's Gram matrix is itself a
#' linear combination of the explanatory variables' kernel widths (faithful
#' to the upstream implementation; see original 207-208 lines). This is an
#' unusual design choice but is deliberately preserved as-is.
#'
#' @param Y data matrix (residualized already, if applicable)
#' @param xi target column index
#' @param xj_list explanatory column indices (length >= 2)
#' @return list(resid = residual vector, coef = coefficient vector)
#' @keywords internal
mlhsicr_regression <- function(Y, xi, xj_list) {
  n <- nrow(Y)
  y_target <- Y[, xi]
  X_exp <- Y[, xj_list, drop = FALSE]

  ols <- rcd_ols_resid_coef(y_target, X_exp)
  init_coef <- as.vector(ols$coef)

  width_list <- vapply(xj_list, function(j) hsic_kernel_width(Y[, j]), numeric(1))
  Lc_list <- Map(function(j, w) hsic_gram_matrix(Y[, j], w)$Kc, xj_list, width_list)
  width_xi <- hsic_kernel_width(y_target)

  objective <- function(coef) {
    resid <- y_target - as.vector(X_exp %*% coef)
    width <- width_xi - sum(coef * width_list)
    if (!is.finite(width) || width <= 0) return(1e10)
    Kc <- hsic_gram_matrix(resid, width)$Kc
    total <- 0
    for (Lc in Lc_list) {
      total <- total + sum(t(Kc) * Lc) / n
    }
    total
  }

  opt <- stats::optim(par = init_coef, fn = objective, method = "L-BFGS-B")

  if (!is.finite(opt$value) || opt$value >= 1e10 || opt$convergence != 0) {
    # optim() either never left the invalid-width plateau (objective stuck at
    # the 1e10 sentinel) or failed to converge; trust the OLS residual rather
    # than a coefficient vector that was never actually minimizing HSIC.
    return(list(resid = as.vector(ols$resid), coef = init_coef))
  }

  coef_final <- opt$par
  resid_final <- y_target - as.vector(X_exp %*% coef_final)

  list(resid = resid_final, coef = coef_final)
}


# =============================================================================
# Internal: extract_ancestors (original 254-316 lines)
# =============================================================================


#' Whether `xi` can be excluded from sink candidacy given the ancestor sets
#' known so far
#'
#' Faithful port of the exclusion rule in `rcd.py` (original 165-174 lines):
#' `xi` cannot be the sink of `U` if (a) `xi` is already known to be an
#' ancestor of some other member of `U`, or (b) every other member of `U` is
#' already known to be an ancestor of `xi`.
#'
#' @param M current ancestor-set list
#' @param U variable-set under consideration (unused except via `xj_list`,
#'   kept for signature parity with the upstream port)
#' @param xi candidate variable
#' @param xj_list `U` minus `xi`
#' @return TRUE if `xi` should be excluded from sink candidacy
#' @keywords internal
exists_ancestor_in_U <- function(M, U, xi, xj_list) {
  for (xj in xj_list) {
    if (xi %in% M[[xj]]) return(TRUE)
  }
  if (length(xj_list) > 0 && all(xj_list %in% M[[xi]])) return(TRUE)
  FALSE
}


#' Whether the OLS (or MLHSICR) residual of `xi` on `xj_list` is independent
#' of every `xj`
#'
#' Faithful port of `rcd.py`'s `is_independent()` helper (original 226-252
#' lines).
#'
#' @param Y residual matrix (from [extract_ancestors()])
#' @param xi candidate sink variable
#' @param xj_list explanatory-variable indices
#' @param MLHSICR whether to retry with [mlhsicr_regression()] on failure
#' @param independence "hsic" or "fcorr"
#' @param ind_alpha significance level (hsic only)
#' @param ind_corr rejection threshold (fcorr only)
#' @return TRUE if independent
#' @keywords internal
is_independent_of_resid <- function(Y, xi, xj_list, MLHSICR, independence,
                                    ind_alpha, ind_corr) {
  fit <- rcd_ols_resid_coef(Y[, xi], Y[, xj_list, drop = FALSE])
  resid <- fit$resid

  all_independent <- TRUE
  for (xj in xj_list) {
    if (!rcd_is_independent(resid, Y[, xj], independence, ind_alpha, ind_corr)) {
      all_independent <- FALSE
      break
    }
  }
  if (all_independent) return(TRUE)

  if (length(xj_list) == 1L || !MLHSICR) return(FALSE)

  mlh <- mlhsicr_regression(Y, xi, xj_list)
  resid2 <- mlh$resid
  for (xj in xj_list) {
    if (!rcd_is_independent(resid2, Y[, xj], independence, ind_alpha, ind_corr)) {
      return(FALSE)
    }
  }
  TRUE
}


#' Extract the ancestor sets M(x_i) for every variable
#'
#' Faithful port of `rcd.py`'s `extract_ancestors()` (original 254-316
#' lines). Repeatedly scans variable subsets of increasing size, growing each
#' variable's ancestor set whenever a unique "sink" is found within a subset.
#' See `dev/rcd-implementation.md` section 3.2 for the full algorithm
#' description and the reasoning behind the cache-update position.
#'
#' @param X (uncentered) data matrix
#' @param max_explanatory_num maximum subset size minus 1
#' @param cor_alpha correlation-test significance level
#' @param ind_alpha independence-test significance level (hsic)
#' @param shapiro_alpha non-Gaussianity significance level
#' @param MLHSICR whether to use MLHSICR regression as a fallback
#' @param independence "hsic" or "fcorr"
#' @param ind_corr F-correlation rejection threshold (fcorr)
#' @return list of length `ncol(X)`; each element is a sorted integer vector
#'   of ancestor indices (possibly empty)
#' @keywords internal
extract_ancestors <- function(X, max_explanatory_num, cor_alpha, ind_alpha,
                              shapiro_alpha, MLHSICR, independence, ind_corr) {
  p <- ncol(X)
  M <- vector("list", p)
  for (i in seq_len(p)) M[[i]] <- integer(0)

  l <- 1L
  hu_history <- new.env(parent = emptyenv())
  combn_cache <- vector("list", max_explanatory_num + 1L)
  iter_count <- 0L
  max_iter <- p * p + 10L

  repeat {
    iter_count <- iter_count + 1L
    if (iter_count > max_iter) {
      stop("extract_ancestors(): exceeded the maximum expected number of iterations; this indicates a bug in ancestor-set convergence.", call. = FALSE)
    }

    changed <- FALSE
    if (is.null(combn_cache[[l]])) {
      combn_cache[[l]] <- utils::combn(seq_len(p), l + 1L, simplify = FALSE)
    }
    U_list <- combn_cache[[l]]

    for (U in U_list) {
      H_U <- Reduce(intersect, M[U])
      key <- paste(U, collapse = ",")

      cached <- if (exists(key, envir = hu_history, inherits = FALSE)) {
        get(key, envir = hu_history, inherits = FALSE)
      } else {
        NULL
      }
      if (!is.null(cached) && identical(cached, H_U)) next

      if (length(H_U) == 0) {
        Y <- X
      } else {
        Y <- matrix(0, nrow(X), p)
        for (u in U) {
          Y[, u] <- rcd_ols_resid_coef(X[, u], X[, H_U, drop = FALSE])$resid
        }
      }

      if (!rcd_is_non_gaussian_all(Y, U, shapiro_alpha)) next

      all_correlated <- TRUE
      pairs <- utils::combn(U, 2, simplify = FALSE)
      for (pr in pairs) {
        if (!rcd_is_correlated(Y[, pr[1]], Y[, pr[2]], cor_alpha)) {
          all_correlated <- FALSE
          break
        }
      }
      if (!all_correlated) next

      sink_set <- integer(0)
      for (xi in U) {
        xj_list <- setdiff(U, xi)
        if (exists_ancestor_in_U(M, U, xi, xj_list)) next
        if (is_independent_of_resid(Y, xi, xj_list, MLHSICR, independence,
                                     ind_alpha, ind_corr)) {
          sink_set <- c(sink_set, xi)
        }
      }

      if (length(sink_set) == 1L) {
        xi <- sink_set[1]
        new_ancestors <- setdiff(setdiff(U, xi), M[[xi]])
        if (length(new_ancestors) > 0) {
          M[[xi]] <- sort(union(M[[xi]], new_ancestors))
          changed <- TRUE
        }
      }

      assign(key, H_U, envir = hu_history)
    }

    if (changed) {
      l <- 1L
    } else if (l < max_explanatory_num) {
      l <- l + 1L
    } else {
      break
    }
  }

  M
}


# =============================================================================
# Internal: extract_parents (original 318-343 lines)
# =============================================================================


#' Extract parents from the ancestor sets
#'
#' Faithful port of `rcd.py`'s `extract_parents()` (original 318-343 lines).
#'
#' @param X (uncentered) data matrix
#' @param M ancestor-set list from [extract_ancestors()]
#' @param cor_alpha correlation-test significance level
#' @return list of length `ncol(X)`; each element is a sorted integer vector
#'   of parent indices (subset of the corresponding ancestor set)
#' @keywords internal
extract_parents <- function(X, M, cor_alpha) {
  p <- ncol(X)
  P <- vector("list", p)
  for (i in seq_len(p)) P[[i]] <- integer(0)

  for (xi in seq_len(p)) {
    for (xj in M[[xi]]) {
      other_ancestors_i <- setdiff(M[[xi]], xj)
      zi <- if (length(other_ancestors_i) == 0) {
        X[, xi]
      } else {
        rcd_ols_resid_coef(X[, xi], X[, other_ancestors_i, drop = FALSE])$resid
      }

      common <- intersect(M[[xi]], M[[xj]])
      wj <- if (length(common) == 0) {
        X[, xj]
      } else {
        rcd_ols_resid_coef(X[, xj], X[, common, drop = FALSE])$resid
      }

      if (rcd_is_correlated(zi, wj, cor_alpha)) {
        P[[xi]] <- sort(union(P[[xi]], xj))
      }
    }
  }
  P
}


# =============================================================================
# Internal: extract_vars_sharing_confounders (original 352-366 lines)
# =============================================================================


#' Detect variable pairs sharing an unobserved latent confounder
#'
#' Faithful port of `rcd.py`'s `extract_vars_sharing_confounders()` (original
#' 352-366 lines). Only pairs with no parent-child relationship in either
#' direction are considered.
#'
#' @param X (uncentered) data matrix
#' @param P parent list from [extract_parents()]
#' @param cor_alpha correlation-test significance level
#' @return list of length `ncol(X)`; each element is a sorted integer vector
#'   of indices sharing a latent confounder with that variable (symmetric)
#' @keywords internal
extract_vars_sharing_confounders <- function(X, P, cor_alpha) {
  p <- ncol(X)
  C <- vector("list", p)
  for (i in seq_len(p)) C[[i]] <- integer(0)

  pairs <- utils::combn(seq_len(p), 2, simplify = FALSE)
  for (pr in pairs) {
    i <- pr[1]
    j <- pr[2]
    if (j %in% P[[i]] || i %in% P[[j]]) next

    ri <- if (length(P[[i]]) == 0) X[, i] else rcd_ols_resid_coef(X[, i], X[, P[[i]], drop = FALSE])$resid
    rj <- if (length(P[[j]]) == 0) X[, j] else rcd_ols_resid_coef(X[, j], X[, P[[j]], drop = FALSE])$resid

    if (rcd_is_correlated(ri, rj, cor_alpha)) {
      C[[i]] <- sort(union(C[[i]], j))
      C[[j]] <- sort(union(C[[j]], i))
    }
  }
  C
}


# =============================================================================
# Internal: adjacency matrix construction (original 368-408 lines)
# =============================================================================


#' Build the RCD adjacency matrix from parents and confounder pairs
#'
#' Faithful port of `rcd.py`'s adjacency-matrix construction (original
#' 368-408 lines). `B[i, j]` is the coefficient of `j -> i`, matching the
#' lingamr convention; no transpose is needed.
#'
#' @param X (uncentered) data matrix
#' @param P parent list from [extract_parents()]
#' @param C confounder-pair list from [extract_vars_sharing_confounders()]
#' @return adjacency matrix B (n_features x n_features), with `NA` entries
#'   for confounder pairs
#' @keywords internal
build_adjacency_matrix_rcd <- function(X, P, C) {
  p <- ncol(X)
  B <- matrix(0, p, p)

  for (xi in seq_len(p)) {
    parents <- sort(P[[xi]])
    if (length(parents) > 0) {
      coef <- rcd_ols_resid_coef(X[, xi], X[, parents, drop = FALSE])$coef
      B[xi, parents] <- coef
    }
  }

  for (xi in seq_len(p)) {
    if (length(C[[xi]]) > 0) {
      B[xi, C[[xi]]] <- NA
    }
  }

  colnames(B) <- rownames(B) <- colnames(X)
  B
}


# =============================================================================
# Public: lingam_rcd()
# =============================================================================


#' RCD (Repetitive Causal Discovery)
#'
#' A causal discovery method robust against latent confounders. Unlike
#' [lingam_direct()] or [lingam_parce()], RCD does not attempt to recover a
#' full or partial causal order. Instead, it repeatedly extracts each
#' variable's **ancestor set** by scanning variable subsets of increasing
#' size (`extract_ancestors()`), narrows each ancestor set down to direct
#' **parents** (`extract_parents()`), and finally tests remaining
#' parent-free pairs for a shared latent confounder
#' (`extract_vars_sharing_confounders()`). Pairs found to share a latent
#' confounder are marked `NA` in the adjacency matrix rather than estimated.
#'
#' @param X Numeric matrix (n_samples x n_features), data frame or matrix
#' @param max_explanatory_num Maximum number of explanatory variables
#'   considered when searching for ancestors (i.e. the search scans variable
#'   subsets of size up to `max_explanatory_num + 1`). Larger values increase
#'   statistical power but grow combinatorially in cost. Must be an integer
#'   of 1 or more.
#' @param cor_alpha Significance level for the Pearson correlation tests used
#'   throughout the algorithm (ancestor-subset screening, parent extraction,
#'   confounder-pair detection). Must be non-negative.
#' @param ind_alpha Significance level for the HSIC independence test (used
#'   when `independence = "hsic"`). Must be non-negative.
#' @param shapiro_alpha Significance level for the Shapiro-Wilk
#'   non-Gaussianity test used when screening candidate ancestor subsets.
#'   Must be non-negative.
#' @param MLHSICR If `TRUE`, falls back to HSIC-sum-minimizing regression
#'   (instead of OLS) when the OLS residual is not independent of the
#'   explanatory variables and more than one explanatory variable is
#'   present. Substantially increases computation time.
#' @param independence Independence measure used for the sink search:
#'   "hsic" (default) uses the HSIC gamma-approximation test; "fcorr" uses
#'   the F-correlation (kernel canonical correlation) and rejects based on
#'   `ind_corr` instead of a p-value.
#' @param ind_corr Threshold on the F-correlation value, used only when
#'   `independence = "fcorr"`. Must be non-negative. Ignored when
#'   `independence = "hsic"`.
#' @return An `RCDResult` object (list) containing:
#' * `adjacency_matrix`: adjacency matrix B (n_features x n_features).
#'   **Convention: `B[i, j]` is the causal coefficient from variable j to
#'   variable i (j -> i)**, same as [lingam_direct()]. Entries between two
#'   variables found to share a latent confounder are `NA`.
#' * `ancestors_list`: a list of length `n_features`; element `i` is the
#'   sorted integer vector of variables found to be ancestors of variable
#'   `i` (possibly empty). Unlike [lingam_parce()], there is no
#'   `causal_order`: RCD estimates ancestor relations directly rather than a
#'   total or partial order.
#' @details
#' The algorithm has three stages: (1) `extract_ancestors()` grows each
#' variable's ancestor set by repeatedly scanning variable subsets; (2)
#' `extract_parents()` narrows ancestor sets down to direct parents; (3)
#' `extract_vars_sharing_confounders()` tests remaining parent-free pairs for
#' a shared latent confounder. `NA` entries in `adjacency_matrix` mean the
#' corresponding pair is suspected to share a latent confounder, not that no
#' relationship was estimated.
#'
#' `max_explanatory_num` controls both statistical power and computational
#' cost: stage 1 scans `choose(n_features, k)` subsets for each subset size
#' `k` up to `max_explanatory_num + 1`, and each subset requires several
#' HSIC tests (each O(n^2) in the sample size when `independence = "hsic"`).
#' Cost grows quickly with both the number of variables and `n`.
#'
#' `MLHSICR = TRUE` replaces the OLS residual in the independence check with
#' a residual obtained by directly minimizing the sum of HSIC statistics
#' between the residual and each explanatory variable via numerical
#' optimization (`stats::optim(method = "L-BFGS-B")`). This can recover
#' independence in cases where OLS cannot, but requires re-optimizing for
#' every candidate subset where the OLS residual fails, and is therefore
#' substantially slower.
#'
#' The Shapiro-Wilk test (`stats::shapiro.test()`) used for the
#' non-Gaussianity check is limited to `n` between 3 and 5000. For `n` above
#' 5000, a deterministic evenly-spaced subsample of 5000 observations is
#' tested instead (same policy as [test_residual_normality()]), so results
#' remain reproducible without touching the RNG state. This subsampling has
#' no effect when `n <= 5000`.
#'
#' This function does not expose a `bw_method` argument (kernel widths are
#' always the median heuristic; see [hsic_kernel_width()]), unlike some
#' upstream implementations. [lingam_rcd_bootstrap()] does not support
#' [get_causal_order_stability()], since RCD has no causal order.
#' @references
#' Maeda, T. N. and Shimizu, S. (2020). RCD: Repetitive causal discovery of
#' linear non-Gaussian acyclic models with latent confounders. AISTATS 2020,
#' PMLR 108: 735-745.
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#'
#' result <- lingam_rcd(confounded$data)
#' print(result)
#'
#' # The variable pair sharing the latent confounder is left NA
#' result$adjacency_matrix[confounded$confounded_pair, confounded$confounded_pair]
#'
#' # Total effect estimation warns and returns NA for confounded variables
#' estimate_total_effect_rcd(confounded$data, result,
#'   from_index = confounded$confounded_pair[1], to_index = 1
#' )
lingam_rcd <- function(X,
                       max_explanatory_num = 2L,
                       cor_alpha = 0.01,
                       ind_alpha = 0.01,
                       shapiro_alpha = 0.01,
                       MLHSICR = FALSE,
                       independence = "hsic",
                       ind_corr = 0.5) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names

  max_explanatory_num <- suppressWarnings(as.integer(max_explanatory_num))
  if (length(max_explanatory_num) != 1 || is.na(max_explanatory_num) ||
        max_explanatory_num < 1) {
    stop("max_explanatory_num must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(cor_alpha) || length(cor_alpha) != 1 || is.na(cor_alpha) || cor_alpha < 0) {
    stop("cor_alpha must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(ind_alpha) || length(ind_alpha) != 1 || is.na(ind_alpha) || ind_alpha < 0) {
    stop("ind_alpha must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(shapiro_alpha) || length(shapiro_alpha) != 1 || is.na(shapiro_alpha) ||
        shapiro_alpha < 0) {
    stop("shapiro_alpha must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.logical(MLHSICR) || length(MLHSICR) != 1 || is.na(MLHSICR)) {
    stop("MLHSICR must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  independence <- match.arg(independence, c("hsic", "fcorr"))
  if (!is.numeric(ind_corr) || length(ind_corr) != 1 || is.na(ind_corr) || ind_corr < 0) {
    stop("ind_corr must be a non-negative numeric scalar.", call. = FALSE)
  }

  M <- extract_ancestors(X, max_explanatory_num, cor_alpha, ind_alpha,
                         shapiro_alpha, MLHSICR, independence, ind_corr)
  P <- extract_parents(X, M, cor_alpha)
  C <- extract_vars_sharing_confounders(X, P, cor_alpha)
  B <- build_adjacency_matrix_rcd(X, P, C)

  names(M) <- colnames(X)

  result <- list(adjacency_matrix = B, ancestors_list = M)
  class(result) <- "RCDResult"
  result
}


#' Print method for RCDResult
#'
#' @param x RCDResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' result <- lingam_rcd(confounded$data)
#' print(result)
print.RCDResult <- function(x, digits = 3, ...) {
  var_names <- colnames(x$adjacency_matrix)
  label_for <- function(idx) {
    if (!is.null(var_names)) var_names[idx] else paste0("x", idx - 1L)
  }

  cat("RCD Result\n")
  cat(sprintf("  Variables : %d\n", ncol(x$adjacency_matrix)))
  cat("\nAncestor sets:\n")
  for (i in seq_along(x$ancestors_list)) {
    anc <- x$ancestors_list[[i]]
    anc_str <- if (length(anc) == 0) {
      "{}"
    } else {
      paste0("{", paste(label_for(anc), collapse = ", "), "}")
    }
    cat(sprintf("  M(%s) = %s\n", label_for(i), anc_str))
  }
  cat("\n  (NA entries in the adjacency matrix = suspected shared latent confounder)\n")
  cat("\nAdjacency matrix (row = to, col = from):\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}


# =============================================================================
# Total causal effect and error independence (RCD-specific)
# =============================================================================


#' Validate the return value of lingam_rcd()
#' @keywords internal
validate_rcd_result <- function(x) {
  if (!inherits(x, "RCDResult")) {
    stop("rcd_result must be the return value of lingam_rcd().", call. = FALSE)
  }
}


#' Resolve a variable index or name to a 1-based integer index
#' @keywords internal
rcd_resolve_index <- function(idx, arg_name, n_features, col_names) {
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


#' Estimate the total causal effect between two variables (RCD)
#'
#' Analogous to [estimate_total_effect()], but for [lingam_rcd()] results,
#' which may contain `NA` entries in the adjacency matrix.
#'
#' @param X Original data (matrix or data.frame)
#' @param rcd_result Return value of [lingam_rcd()]
#' @param from_index Cause variable (1-based index or variable name)
#' @param to_index Effect variable (1-based index or variable name)
#' @param method Regression method ("ols", "lasso", "adaptive_lasso", "ridge"). Default is adaptive_lasso
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle"). Default is BIC
#' @param init_method Method for estimating the initial weights of adaptive LASSO regression ("ols" or "ridge")
#' @return Estimated total causal effect, or `NA` (with a warning) if
#'   `from_index` is part of a suspected latent confounder pair (its parents
#'   cannot be identified). Also warns (without altering the estimate) if
#'   `to_index` is an ancestor of `from_index` according to `ancestors_list`,
#'   since that is inconsistent with a `from -> to` effect.
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' result <- lingam_rcd(confounded$data)
#'
#' # A well-identified pair returns a numeric estimate
#' estimate_total_effect_rcd(confounded$data, result, from_index = 6, to_index = 1)
estimate_total_effect_rcd <- function(X, rcd_result, from_index, to_index,
                                      method = "adaptive_lasso", lambda = "BIC",
                                      init_method = "ols") {
  validate_rcd_result(rcd_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  X <- as.matrix(X)
  col_names <- colnames(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  n_features <- ncol(X)
  B <- rcd_result$adjacency_matrix
  M <- rcd_result$ancestors_list
  if (ncol(B) != n_features) {
    stop(
      "X has ", n_features, " variables but rcd_result was estimated from ",
      ncol(B), " variables.",
      call. = FALSE
    )
  }

  from_index <- rcd_resolve_index(from_index, "from_index", n_features, col_names)
  to_index <- rcd_resolve_index(to_index, "to_index", n_features, col_names)
  if (from_index == to_index) stop("from_index and to_index must differ.")

  from_label <- if (!is.null(col_names)) col_names[from_index] else paste0("x", from_index)
  to_label <- if (!is.null(col_names)) col_names[to_index] else paste0("x", to_index)

  if (to_index %in% M[[from_index]]) {
    warning(sprintf(
      "%s is an ancestor of %s according to ancestors_list; the requested direction (%s -> %s) is inconsistent with the estimated ancestor relations.",
      to_label, from_label, from_label, to_label
    ))
  }

  if (anyNA(B[from_index, ])) {
    warning(sprintf(
      "%s is part of a suspected latent confounder pair; ", from_label
    ), "total effect cannot be estimated.")
    return(NA_real_)
  }

  parents <- which(abs(B[from_index, ]) > 0)
  predictors <- unique(c(from_index, parents))
  from_pos <- which(predictors == from_index)

  y <- X[, to_index]
  Xp <- X[, predictors, drop = FALSE]

  coefs <- fit_coef_by_method(y, Xp, method, lambda, init_method)

  coefs[from_pos]
}


#' Compute p-values for the independence of RCD residuals (HSIC-based)
#'
#' Analogous to [get_error_independence_p_values_parce()], but for
#' [lingam_rcd()] results. Returns `NA` for any pair involving a variable
#' whose row or column in the adjacency matrix contains `NA` (residuals
#' cannot be computed for those variables).
#'
#' @param X Original data (matrix or data.frame)
#' @param rcd_result Return value of [lingam_rcd()]
#' @return matrix of p-values (n_features x n_features)
#' @export
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' result <- lingam_rcd(confounded$data)
#' round(get_error_independence_p_values_rcd(confounded$data, result), 3)
get_error_independence_p_values_rcd <- function(X, rcd_result) {
  validate_rcd_result(rcd_result)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)

  B <- rcd_result$adjacency_matrix
  n_features <- ncol(B)
  if (ncol(X) != n_features) {
    stop(
      "X has ", ncol(X), " variables but rcd_result was estimated from ",
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
