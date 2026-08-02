# =============================================================================
# RESIT (Regression with Subsequent Independence Test)
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/resit.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
#
# Reuses:
#   R/hsic.r                : hsic_test_gamma() (multivariate matrix path)
#   R/lingam_direct.r       : validate_no_degenerate_columns(),
#                             validate_prior_knowledge(), get_var_names()
#   R/search_causal_order.r : extract_partial_orders()
#   R/lingam_parce.r        : parce_search_candidate()
# =============================================================================


#' Check that mgcv is available
#'
#' Raise an informative error when the suggested package mgcv is missing.
#' Same pattern as [check_glmnet_available()].
#'
#' @keywords internal
check_mgcv_available <- function() {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required for regressor = \"gam\". Please install it.",
      call. = FALSE
    )
  }
}


#' Built-in GAM regressor for RESIT
#'
#' Fits `mgcv::gam(y ~ s(V1) + s(V2) + ...)` and returns the fitted values.
#' The basis dimension `k` is capped so that the model stays estimable on
#' small or resampled data: bootstrap resamples duplicate rows, which lowers
#' the number of unique covariate values, and mgcv errors out when a basis
#' has more dimensions than unique values (or more total coefficients than
#' observations). When even `k = 3` is not affordable the smooth terms
#' degrade to plain linear terms.
#'
#' @param X numeric matrix of predictors (n x d, d >= 1)
#' @param y numeric response vector (length n)
#' @return fitted values (length n)
#' @keywords internal
resit_gam_fitted <- function(X, y) {
  d <- ncol(X)
  n <- nrow(X)
  df <- as.data.frame(X)
  names(df) <- paste0("V", seq_len(d))
  df$.y <- y

  n_unique <- vapply(
    seq_len(d), function(j) length(unique(X[, j])), integer(1)
  )
  k <- min(10L, as.integer(floor((n - 1L) / d)), min(n_unique) - 1L)

  terms <- if (k >= 3L) {
    sprintf("s(V%d, k = %d)", seq_len(d), k)
  } else {
    paste0("V", seq_len(d))
  }
  f <- stats::as.formula(paste(".y ~", paste(terms, collapse = " + ")))
  as.vector(stats::fitted(mgcv::gam(f, data = df)))
}


#' Resolve the regressor argument to a fitting function and its label
#'
#' The label is normalized (`match.arg` accepts partial matches such as
#' `"g"`), so the `regressor` field of the result always reports the
#' canonical name of what actually ran.
#'
#' @param regressor either a string (currently only `"gam"`) or a function
#'   `function(X, y)` returning fitted values
#' @return `list(fn = <function (X, y) -> fitted values>, label = <string>)`
#' @keywords internal
resit_make_regressor <- function(regressor) {
  if (is.function(regressor)) {
    return(list(fn = regressor, label = "user function"))
  }
  if (is.character(regressor) && length(regressor) == 1 && !is.na(regressor)) {
    regressor <- match.arg(regressor, c("gam"))
    check_mgcv_available()
    return(list(fn = resit_gam_fitted, label = regressor))
  }
  stop(
    "regressor must be \"gam\" or a function(X, y) returning fitted values.",
    call. = FALSE
  )
}


#' Call a regressor function and validate its return value
#'
#' User-supplied regressors are called many times deep inside the order
#' search, so a malformed return value is caught here with a clear message
#' instead of surfacing as a cryptic arithmetic error.
#'
#' @param reg_fn function `(X, y) -> fitted values`
#' @param Xp predictor matrix (n x d)
#' @param y response vector (length n)
#' @return fitted values as a plain numeric vector
#' @keywords internal
resit_fitted <- function(reg_fn, Xp, y) {
  fv <- reg_fn(Xp, y)
  if (!is.numeric(fv) || length(fv) != length(y) || anyNA(fv)) {
    stop(
      "regressor must return a numeric vector of fitted values with ",
      "length nrow(X) and no NA.",
      call. = FALSE
    )
  }
  as.vector(fv)
}


#' RESIT step 1: estimate the topological (causal) order
#'
#' Faithful port of `RESIT._estimate_order()` (resit.py). Repeatedly finds
#' the most sink-like variable: the one whose regression residual on all
#' other remaining variables has the smallest HSIC dependence statistic.
#' Sinks are identified first and prepended, so the returned order runs
#' from source to sink.
#'
#' Deviation from the Python original: [parce_search_candidate()] falls
#' back to the full remaining set when prior knowledge excludes every
#' candidate (the Python code would crash on an empty argmin).
#'
#' @param X data matrix
#' @param reg_fn regressor function from [resit_make_regressor()]
#' @param partial_orders (from, to) matrix from [extract_partial_orders()],
#'   or NULL
#' @param Aknw validated prior-knowledge matrix (negative entries already
#'   NA), or NULL
#' @return list(pa = list of parent index vectors per variable,
#'   order = causal order, source first)
#' @keywords internal
resit_estimate_order <- function(X, reg_fn, partial_orders, Aknw) {
  p <- ncol(X)
  S <- seq_len(p)
  pa <- rep(list(integer(0)), p)
  pi_order <- integer(0)

  for (step in seq_len(p)) {
    Sc <- parce_search_candidate(S, partial_orders)
    if (length(Sc) == 1L) {
      k <- Sc[1L]
    } else {
      hsic_stats <- vapply(Sc, function(kk) {
        predictors <- setdiff(S, kk)
        fitted_vals <- resit_fitted(
          reg_fn, X[, predictors, drop = FALSE], X[, kk]
        )
        residual <- X[, kk] - fitted_vals
        # dependence between the residual and the JOINT distribution of the
        # predictors (multivariate kernel), as in the Python original
        hsic_test_gamma(residual, X[, predictors, drop = FALSE])$stat
      }, numeric(1))
      # ties resolve to the first candidate, like numpy's argmin; setdiff
      # preserves the order of S, so candidates enumerate in the same order
      # as the Python list comprehension
      k <- Sc[which.min(hsic_stats)]
    }

    S <- setdiff(S, k)
    pa[[k]] <- S # parents = remaining variables AFTER removing k (verbatim)
    pi_order <- c(k, pi_order) # pi.insert(0, k): prepend, sink found first
    if (!is.null(partial_orders) && nrow(partial_orders) > 0) {
      # drop partial orders pointing at the resolved variable (to == k)
      partial_orders <- partial_orders[
        partial_orders[, 2] != k, ,
        drop = FALSE
      ]
    }
  }

  # Prior knowledge prunes the parent sets: keep parents with pk != 0,
  # where NA (unknown) also counts as "keep" -- numpy's NaN != 0 is True,
  # but R's `!=` would return NA, hence the explicit is.na() branch.
  # (The Python original zeroes the pk diagonal first; pa[[k]] can never
  # contain k itself, so that step is a no-op here.)
  if (!is.null(Aknw)) {
    for (k in seq_len(p)) {
      v <- pa[[k]]
      keep <- vapply(v, function(par) {
        a <- Aknw[k, par]
        is.na(a) || a != 0
      }, logical(1))
      pa[[k]] <- v[keep]
    }
  }

  list(pa = pa, order = pi_order)
}


#' RESIT step 2: remove superfluous edges
#'
#' Faithful port of `RESIT._remove_edges()` (resit.py). For each variable in
#' causal order (skipping the source) and each of its candidate parents `l`,
#' the variable is regressed on its current parents excluding `l`; if the
#' residual is independent (HSIC p-value > alpha) of the variable's original
#' parent set, the edge from `l` is dropped.
#'
#' Note the two distinct parent sets, mirroring the original exactly:
#' `parents_snapshot` is the parent set frozen before the inner loop and is
#' always the second HSIC argument, while `predictors` is derived from the
#' *current* (already pruned) parent set minus `l`.
#'
#' @param X data matrix
#' @param pa list of parent index vectors from [resit_estimate_order()]
#' @param pi_order causal order (source first)
#' @param reg_fn regressor function
#' @param alpha significance level of the HSIC test
#' @return pruned `pa` list
#' @keywords internal
resit_remove_edges <- function(X, pa, pi_order, reg_fn, alpha) {
  p <- length(pi_order)
  # Python: range(1, n_features) -- the source pi_order[1] has no parents.
  # (2:p would degenerate for p = 1, but ncol >= 2 is enforced upstream.)
  for (k_pos in 2:p) {
    target <- pi_order[k_pos]
    parents_snapshot <- pa[[target]] # pa[pi[k]].copy() in the original
    if (length(parents_snapshot) == 0L) next

    # the second HSIC argument is the same joint parent matrix for every l;
    # precompute its Gram parts once per target
    pre_parents <- hsic_precompute(X[, parents_snapshot, drop = FALSE])

    for (l in parents_snapshot) {
      predictors <- setdiff(pa[[target]], l)
      if (length(predictors) >= 1L) {
        fitted_vals <- resit_fitted(
          reg_fn, X[, predictors, drop = FALSE], X[, target]
        )
        residual <- X[, target] - fitted_vals
      } else {
        residual <- X[, target]
      }
      hsic_p <- hsic_gamma_from_pre(hsic_precompute(residual), pre_parents)$p
      if (hsic_p > alpha) {
        pa[[target]] <- setdiff(pa[[target]], l)
      }
    }
  }
  pa
}


#' RESIT causal discovery for nonlinear additive noise models
#'
#' R port of the RESIT algorithm (Regression with Subsequent Independence
#' Test) from the Python `lingam` package. RESIT assumes a nonlinear
#' additive noise model `x_i = f_i(parents(x_i)) + e_i` and recovers the
#' causal structure in two phases: (1) the causal order is estimated by
#' repeatedly detaching the most sink-like variable, i.e. the variable
#' whose regression residual on all remaining variables is least dependent
#' on them (smallest HSIC statistic); (2) superfluous parents are pruned by
#' testing, for each candidate parent, whether the residual regressed on
#' the other parents is already independent of the parent set (HSIC
#' p-value greater than `alpha`).
#'
#' @details
#' Because the model is nonlinear, the returned `adjacency_matrix` contains
#' 0/1 edge indicators, **not** connection strengths, and total causal
#' effects are undefined: the Python implementation's
#' `estimate_total_effect()` (always 0) and
#' `get_error_independence_p_values()` (always a zero matrix) are
#' intentionally not ported.
#'
#' Note the direction of `alpha`: a parent is removed when the HSIC p-value
#' exceeds `alpha`, so *larger* values of `alpha` make the test stricter
#' about declaring independence and therefore *keep more edges*.
#'
#' The HSIC test in phase 1 measures dependence between a residual and the
#' joint (multivariate) kernel of up to `ncol(X) - 1` predictors. Following
#' the Python original, variables are not standardized beforehand; if the
#' variable scales differ wildly, the largest-scale variable dominates the
#' kernel distances. Each HSIC call builds n x n Gram matrices, and
#' `O(ncol(X)^2)` regressions and HSIC tests are performed overall, so the
#' method is not recommended for `nrow(X)` in the thousands.
#'
#' At least 6 observations are required (the lower bound of the HSIC
#' gamma-approximation test), which is stricter than the linear methods in
#' this package.
#'
#' Deviation from the Python original when `prior_knowledge` excludes every
#' remaining sink candidate: the candidate set falls back to all remaining
#' variables (the original would fail), matching [lingam_parce()].
#'
#' @param X numeric matrix or data frame of observed variables
#' @param regressor nonlinear regressor used for all internal regressions.
#'   Either the string `"gam"` (default; requires the suggested package
#'   mgcv, and fits a smoothing-spline GAM per regression) or a function
#'   `function(X, y)` that receives a predictor matrix and a response
#'   vector and returns the fitted values as a numeric vector of length
#'   `nrow(X)` (the model is only ever evaluated on its own training data,
#'   mirroring `regressor.fit(X, y); regressor.predict(X)` in Python)
#' @param alpha significance level of the HSIC independence test used for
#'   edge pruning (default: 0.01; must be non-negative)
#' @param prior_knowledge optional prior-knowledge matrix with elements 1
#'   (directed path exists), 0 (no directed path), and -1 or NA (unknown),
#'   with the same `[to, from]` orientation as the adjacency matrix
#' @return An object of class `ResitResult` with elements:
#' * `adjacency_matrix`: (p x p) 0/1 matrix; `B[i, j] = 1` means an edge
#'   `j -> i` (row = to, col = from). Entries are edge indicators, not
#'   coefficients.
#' * `causal_order`: estimated causal order (1-based column positions,
#'   source first).
#' * `regressor`: label of the regressor used (`"gam"` or
#'   `"user function"`).
#' @references
#' J. Peters, J. M. Mooij, D. Janzing, B. Schoelkopf. Causal discovery with
#' continuous additive noise models. Journal of Machine Learning Research,
#' 15: 2009-2053, 2014.
#' @examples
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   nonlinear <- generate_resit_sample(n = 300, seed = 1)
#'   result <- lingam_resit(nonlinear$data)
#'   print(result)
#' }
#' }
#' @export
lingam_resit <- function(X,
                         regressor = "gam",
                         alpha = 0.01,
                         prior_knowledge = NULL) {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 6) {
    stop(
      "X must have at least 6 observations (rows); the HSIC ",
      "gamma-approximation test is undefined below that.",
      call. = FALSE
    )
  }
  if (!is.null(col_names)) colnames(X) <- col_names
  validate_no_degenerate_columns(X)

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0) {
    stop("alpha must be a non-negative numeric scalar.", call. = FALSE)
  }

  reg <- resit_make_regressor(regressor)
  reg_fn <- reg$fn

  Aknw <- NULL
  partial_orders <- NULL
  if (!is.null(prior_knowledge)) {
    Aknw <- validate_prior_knowledge(prior_knowledge, ncol(X))
    partial_orders <- extract_partial_orders(Aknw)
  }

  est <- resit_estimate_order(X, reg_fn, partial_orders, Aknw)
  pa <- resit_remove_edges(X, est$pa, est$order, reg_fn, alpha)

  p <- ncol(X)
  var_names <- get_var_names(X)
  B <- matrix(0, p, p, dimnames = list(var_names, var_names))
  for (i in seq_len(p)) {
    B[i, pa[[i]]] <- 1
  }

  result <- list(
    adjacency_matrix = B,
    causal_order = est$order,
    regressor = reg$label
  )
  class(result) <- "ResitResult"
  result
}


#' Print method for ResitResult
#'
#' @param x ResitResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   nonlinear <- generate_resit_sample(n = 300, seed = 1)
#'   result <- lingam_resit(nonlinear$data)
#'   print(result)
#' }
#' }
print.ResitResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  var_names <- colnames(x$adjacency_matrix)
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("RESIT Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Regressor : %s\n", x$regressor))
  cat(sprintf("  Causal order: %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nAdjacency matrix (row = to, col = from):\n")
  cat("  (entries are 0/1 edge indicators, not coefficients)\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}
