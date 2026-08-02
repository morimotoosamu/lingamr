# =============================================================================
# CAM-UV (Causal Additive Models with Unobserved Variables)
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/camuv.py)
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
#   R/hsic.r          : hsic_test_gamma() (univariate and multivariate paths)
#   R/f_correlation.r : f_correlation()
#   R/lingam_resit.r  : resit_make_regressor() / resit_fitted() (the upstream
#                       hardcodes pygam's LinearGAM; this port uses the same
#                       pluggable regressor abstraction as lingam_resit())
#   R/lingam_direct.r : get_var_names(), validate_no_degenerate_columns()
# =============================================================================


#' Validate CAM-UV prior knowledge and build the forbidden-cause list
#'
#' Faithful port of `CAMUV._make_pk_dict()`. Unlike the matrix-based prior
#' knowledge of the other algorithms in this package, CAM-UV takes a list of
#' pairs `(i, j)` meaning "variable i cannot be a cause of variable j"
#' (**1-based** in this port; the Python original is 0-based).
#'
#' @param prior_knowledge NULL, a 2-column numeric matrix, or a list of
#'   length-2 numeric vectors
#' @param d number of variables
#' @return NULL, or a list of length `d`; element `j` holds the variables
#'   forbidden from being a cause of variable `j`
#' @keywords internal
camuv_make_pk_dict <- function(prior_knowledge, d) {
  if (is.null(prior_knowledge)) {
    return(NULL)
  }

  bad_format <- function() {
    stop(
      "prior_knowledge must be a 2-column numeric matrix or a list of ",
      "length-2 numeric vectors; each pair (i, j) means \"variable i ",
      "cannot be a cause of variable j\" (1-based).",
      call. = FALSE
    )
  }

  if (is.list(prior_knowledge) && !is.data.frame(prior_knowledge)) {
    ok <- vapply(
      prior_knowledge,
      function(p) is.numeric(p) && length(p) == 2 && !anyNA(p),
      logical(1)
    )
    if (length(prior_knowledge) == 0 || !all(ok)) bad_format()
    pk <- do.call(rbind, lapply(prior_knowledge, as.integer))
  } else if (is.matrix(prior_knowledge) && is.numeric(prior_knowledge) &&
               ncol(prior_knowledge) == 2 && !anyNA(prior_knowledge)) {
    if (nrow(prior_knowledge) == 0) bad_format()
    pk <- matrix(as.integer(prior_knowledge), ncol = 2)
  } else {
    bad_format()
  }

  if (any(pk < 1L) || any(pk > d)) {
    stop(
      "prior_knowledge indices must be between 1 and ", d,
      " (1-based column positions).",
      call. = FALSE
    )
  }
  if (any(pk[, 1] == pk[, 2])) {
    stop(
      "prior_knowledge pairs must have two distinct variables.",
      call. = FALSE
    )
  }

  forbidden <- rep(list(integer(0)), d)
  for (r in seq_len(nrow(pk))) {
    forbidden[[pk[r, 2]]] <- union(forbidden[[pk[r, 2]]], pk[r, 1])
  }
  forbidden
}


#' Residual of one variable regressed on a set of explanatory variables
#'
#' Faithful port of `CAMUV._get_residual()`. With no explanatory variables
#' the column is returned as is.
#'
#' @param X data matrix
#' @param explained_i column index of the explained variable
#' @param explanatory_ids column indices of the explanatory variables
#'   (possibly empty)
#' @param reg_fn regressor function from [resit_make_regressor()]
#' @return residual vector (length n)
#' @keywords internal
camuv_get_residual <- function(X, explained_i, explanatory_ids, reg_fn) {
  if (length(explanatory_ids) == 0) {
    return(X[, explained_i])
  }
  fitted_vals <- resit_fitted(
    reg_fn, X[, explanatory_ids, drop = FALSE], X[, explained_i]
  )
  X[, explained_i] - fitted_vals
}


#' Independence judgment with an explicit threshold, returning the value
#'
#' Faithful port of `CAMUV._is_independent_by()`. For `"hsic"` the value is
#' the gamma-approximated p-value and independence means `value > threshold`;
#' for `"fcorr"` the value is the F-correlation and independence means
#' `value < threshold`. Note the opposite directions.
#'
#' `f_correlation()` is univariate, so the fcorr path requires single-column
#' inputs; [lingam_camuv()] rejects `independence = "fcorr"` with
#' `num_explanatory_vals > 2` up front, which is the only way a multivariate
#' input could reach this point (the upstream implementation breaks on that
#' combination too).
#'
#' @param X numeric vector or single/multi-column matrix
#' @param Y numeric vector or single/multi-column matrix
#' @param threshold rejection threshold (see above)
#' @param independence "hsic" or "fcorr"
#' @return list(independent = logical, value = numeric)
#' @keywords internal
camuv_is_independent_by <- function(X, Y, threshold, independence) {
  if (independence == "hsic") {
    value <- hsic_test_gamma(X, Y)$p
    list(independent = value > threshold, value = value)
  } else {
    X <- as_hsic_matrix(X)
    Y <- as_hsic_matrix(Y)
    if (ncol(X) > 1L || ncol(Y) > 1L) {
      stop(
        "internal error: f_correlation() only supports univariate inputs.",
        call. = FALSE
      )
    }
    value <- f_correlation(X[, 1L], Y[, 1L])
    list(independent = value < threshold, value = value)
  }
}


#' Independence judgment at the configured threshold
#'
#' Faithful port of `CAMUV._is_independent()`: alpha for hsic, ind_corr for
#' fcorr.
#'
#' @inheritParams camuv_is_independent_by
#' @param alpha significance level (hsic only)
#' @param ind_corr rejection threshold (fcorr only)
#' @return TRUE if independent
#' @keywords internal
camuv_is_independent <- function(X, Y, independence, alpha, ind_corr) {
  threshold <- if (independence == "hsic") alpha else ind_corr
  camuv_is_independent_by(X, Y, threshold, independence)$independent
}


#' Pairwise dependence neighborhoods
#'
#' Faithful port of `CAMUV._get_neighborhoods()`: `N[[i]]` holds the
#' variables whose raw column is *dependent* on column i.
#'
#' @param X data matrix
#' @inheritParams camuv_is_independent
#' @return list of length `ncol(X)` of integer vectors
#' @keywords internal
camuv_get_neighborhoods <- function(X, independence, alpha, ind_corr) {
  d <- ncol(X)
  N <- rep(list(integer(0)), d)
  if (independence == "hsic") {
    # every column enters d - 1 pairwise tests; precompute each column's
    # O(n^2) Gram parts once instead of once per pair
    pres <- lapply(seq_len(d), function(k) hsic_precompute(X[, k]))
    for (i in seq_len(d - 1L)) {
      for (j in seq.int(i + 1L, d)) {
        if (!(hsic_gamma_from_pre(pres[[i]], pres[[j]])$p > alpha)) {
          N[[i]] <- c(N[[i]], j)
          N[[j]] <- c(N[[j]], i)
        }
      }
    }
    return(N)
  }
  for (i in seq_len(d - 1L)) {
    for (j in seq.int(i + 1L, d)) {
      if (!camuv_is_independent(X[, i], X[, j], independence, alpha, ind_corr)) {
        N[[i]] <- c(N[[i]], j)
        N[[j]] <- c(N[[j]], i)
      }
    }
  }
  N
}


#' Check that no pair within a variable set already has a parent relation
#'
#' Faithful port of `CAMUV._check_identified_causality()`.
#'
#' @param vars integer vector of variable indices
#' @param P current parent list
#' @return TRUE if no pair in `vars` has an identified causal relation yet
#' @keywords internal
camuv_check_identified_causality <- function(vars, P) {
  k <- length(vars)
  for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      i <- vars[a]
      j <- vars[b]
      if (j %in% P[[i]] || i %in% P[[j]]) {
        return(FALSE)
      }
    }
  }
  TRUE
}


#' Check that every candidate parent is in the child's neighborhood
#'
#' Faithful port of `CAMUV._check_correlation()`.
#'
#' @param child candidate child index
#' @param parents candidate parent indices
#' @param N neighborhood list from [camuv_get_neighborhoods()]
#' @return TRUE if all parents are dependent on the child
#' @keywords internal
camuv_check_correlation <- function(child, parents, N) {
  all(parents %in% N[[child]])
}


#' Check whether prior knowledge forbids a candidate parent set
#'
#' Faithful port of `CAMUV._check_prior_knowledge()`. Returns TRUE when the
#' combination is *blocked* (some candidate parent is forbidden from being a
#' cause of the child).
#'
#' @param pk_forbidden forbidden-cause list from [camuv_make_pk_dict()], or NULL
#' @param parents candidate parent indices
#' @param child candidate child index
#' @return TRUE if blocked by prior knowledge
#' @keywords internal
camuv_check_prior_knowledge <- function(pk_forbidden, parents, child) {
  if (is.null(pk_forbidden)) {
    return(FALSE)
  }
  any(parents %in% pk_forbidden[[child]])
}


#' Select the most sink-like child within a variable set
#'
#' Faithful port of `CAMUV._get_child()`. For each candidate child, the
#' child is regressed on the candidate parents plus its already-identified
#' parents, and the residual is tested against the *residual matrix* columns
#' of the candidate parents. The candidate whose residual is most
#' independent wins; `prev_independence` starts at 0 for hsic (p-value
#' scale, higher is more independent) and at 1 for fcorr (lower is more
#' independent), and each accepted candidate raises the bar for the next.
#'
#' @param X data matrix
#' @param vars variable subset (integer vector)
#' @param P current parent list
#' @param N neighborhood list
#' @param Y current residual matrix
#' @param pk_forbidden forbidden-cause list, or NULL
#' @param reg_fn regressor function
#' @param get_residual function(v, ids) returning the residual of `X[, v]`
#'   regressed on `X[, ids]`; defaults to an uncached [camuv_get_residual()]
#'   call. [camuv_find_parents()] passes a memoized version, since the same
#'   (child, predictor set) recurs across subset rescans.
#' @inheritParams camuv_is_independent
#' @return list(child = index or NA, independent = logical); `independent`
#'   reports whether the winning candidate's dependence value clears the
#'   configured threshold (alpha / ind_corr)
#' @keywords internal
camuv_get_child <- function(X, vars, P, N, Y, pk_forbidden,
                            independence, alpha, ind_corr, reg_fn,
                            get_residual = NULL) {
  if (is.null(get_residual)) {
    get_residual <- function(v, ids) camuv_get_residual(X, v, ids, reg_fn)
  }
  prev_independence <- if (independence == "hsic") 0.0 else 1.0
  max_independence_child <- NA_integer_

  for (child in vars) {
    parents <- setdiff(vars, child)

    if (camuv_check_prior_knowledge(pk_forbidden, parents, child)) next
    if (!camuv_check_correlation(child, parents, N)) next

    residual <- get_residual(child, union(parents, P[[child]]))
    res <- camuv_is_independent_by(
      matrix(residual, ncol = 1L), Y[, parents, drop = FALSE],
      prev_independence, independence
    )
    if (res$independent) {
      prev_independence <- res$value
      max_independence_child <- child
    }
  }

  independent <- if (independence == "hsic") {
    prev_independence > alpha
  } else {
    prev_independence < ind_corr
  }

  list(child = max_independence_child, independent = independent)
}


#' Check that the child's residual is dependent on each parent's residual
#'
#' Faithful port of `CAMUV._check_independence_withou_K()` (typo in the
#' upstream method name). If the child's current residual is already
#' independent of some candidate parent's residual *without* regressing on
#' the candidate set K, that parent adds no information and the candidate
#' set is rejected.
#'
#' @param parents candidate parent indices
#' @param child candidate child index
#' @param Y current residual matrix
#' @inheritParams camuv_is_independent
#' @return TRUE if the child's residual is dependent on every parent's residual
#' @keywords internal
camuv_check_independence_without_k <- function(parents, child, Y,
                                               independence, alpha, ind_corr) {
  for (parent in parents) {
    if (camuv_is_independent(Y[, child], Y[, parent],
                             independence, alpha, ind_corr)) {
      return(FALSE)
    }
  }
  TRUE
}


#' CAM-UV main search: identify each variable's parents
#'
#' Faithful port of `CAMUV._find_parents()`. Scans variable subsets of
#' increasing size `t` (starting at 2); whenever any parent is identified,
#' `t` resets to 2, otherwise `t` grows until it exceeds
#' `num_explanatory_vals`. The residual matrix `Y` is updated in place each
#' time a variable gains a parent. A final pruning pass removes parents
#' whose residuals turn out to be independent of the child's residual.
#'
#' @param X data matrix
#' @param maxnum_vals maximum subset size (`num_explanatory_vals`)
#' @param N neighborhood list from [camuv_get_neighborhoods()]
#' @param pk_forbidden forbidden-cause list, or NULL
#' @param reg_fn regressor function
#' @inheritParams camuv_is_independent
#' @return list of parent index vectors per variable
#' @keywords internal
camuv_find_parents <- function(X, maxnum_vals, N, pk_forbidden,
                               independence, alpha, ind_corr, reg_fn) {
  d <- ncol(X)
  P <- rep(list(integer(0)), d)
  t <- 2L
  Y <- X

  # The regression residual is fully determined by (variable, exact predictor
  # vector) and X; memoize it so subset rescans (t resets to 2 after every
  # change) and the final pruning pass do not refit identical GAMs.
  res_env <- new.env(parent = emptyenv())
  cached_residual <- function(v, ids) {
    key <- paste(v, paste(ids, collapse = ","), sep = "|")
    if (!exists(key, envir = res_env, inherits = FALSE)) {
      assign(key, camuv_get_residual(X, v, ids, reg_fn), envir = res_env)
    }
    get(key, envir = res_env, inherits = FALSE)
  }

  repeat {
    changed <- FALSE
    # itertools.combinations() yields nothing for t > d; combn() would error
    if (t <= d) {
      combos <- utils::combn(d, t)
      for (ci in seq_len(ncol(combos))) {
        vars <- combos[, ci]

        if (!camuv_check_identified_causality(vars, P)) next

        gc <- camuv_get_child(X, vars, P, N, Y, pk_forbidden,
                              independence, alpha, ind_corr, reg_fn,
                              get_residual = cached_residual)
        if (is.na(gc$child)) next
        if (!gc$independent) next

        child <- gc$child
        parents <- setdiff(vars, child)
        if (!camuv_check_independence_without_k(parents, child, Y,
                                                independence, alpha,
                                                ind_corr)) {
          next
        }

        for (parent in parents) {
          P[[child]] <- union(P[[child]], parent)
          changed <- TRUE
          Y[, child] <- cached_residual(child, P[[child]])
        }
      }
    }

    if (changed) {
      t <- 2L
    } else {
      t <- t + 1L
      if (t > maxnum_vals) break
    }
  }

  # Final pruning: drop parents whose residual is independent of the child's
  # residual computed without that parent. residual_j depends only on
  # (j, P[[j]]) and repeats across children i, so the memoization pays here too.
  for (i in seq_len(d)) {
    non_parents <- integer(0)
    for (j in P[[i]]) {
      residual_i <- cached_residual(i, setdiff(P[[i]], j))
      residual_j <- cached_residual(j, P[[j]])
      if (camuv_is_independent(residual_i, residual_j,
                               independence, alpha, ind_corr)) {
        non_parents <- c(non_parents, j)
      }
    }
    P[[i]] <- setdiff(P[[i]], non_parents)
  }

  P
}


#' Build the CAM-UV adjacency matrix
#'
#' Faithful port of `CAMUV._estimate_adjacency_matrix()`: `B[child, parent]
#' = 1` (edge indicators, not coefficients, since the causal functions are
#' nonlinear), and both entries of each confounded pair are `NA`.
#'
#' @param X data matrix (for dimnames)
#' @param P parent list from [camuv_find_parents()]
#' @param U list of confounded pairs (length-2 integer vectors)
#' @return adjacency matrix B (n_features x n_features)
#' @keywords internal
camuv_build_adjacency_matrix <- function(X, P, U) {
  d <- ncol(X)
  var_names <- get_var_names(X)
  B <- matrix(0, d, d, dimnames = list(var_names, var_names))
  for (i in seq_len(d)) {
    B[i, P[[i]]] <- 1
  }
  for (pair in U) {
    B[pair[1], pair[2]] <- NA
    B[pair[2], pair[1]] <- NA
  }
  B
}


# =============================================================================
# Public: lingam_camuv()
# =============================================================================


#' CAM-UV (Causal Additive Models with Unobserved Variables)
#'
#' A causal discovery method for **nonlinear additive** models that allows
#' unobserved variables. CAM-UV assumes each observed variable is a
#' generalized additive function of its parents plus independent noise, and
#' the causal structure is a DAG. It identifies each variable's direct
#' parents by scanning variable subsets and testing whether the GAM
#' regression residual is independent of the candidate parents
#' (`camuv_find_parents()`); variable pairs whose residuals remain dependent
#' without an identified edge are reported as connected by an **unobserved
#' causal path** (UCP: a directed path through an unobserved variable) or an
#' **unobserved backdoor path** (UBP: a common unobserved ancestor), and are
#' marked `NA` in the adjacency matrix rather than oriented.
#'
#' @param X Numeric matrix (n_samples x n_features), data frame or matrix
#' @param alpha Significance level of the HSIC independence test (used when
#'   `independence = "hsic"`). Must be non-negative.
#' @param num_explanatory_vals Maximum size of the variable subsets scanned
#'   when searching for parents (the maximum number of explanatory
#'   variables considered jointly is `num_explanatory_vals - 1`). Larger
#'   values increase statistical power but grow combinatorially in cost.
#'   Must be an integer of 1 or more.
#' @param independence Independence measure: "hsic" (default) uses the HSIC
#'   gamma-approximation test; "fcorr" uses the F-correlation (kernel
#'   canonical correlation) and rejects based on `ind_corr` instead of a
#'   p-value.
#' @param ind_corr Threshold on the F-correlation value, used only when
#'   `independence = "fcorr"` (independence is declared below this value).
#'   Must be non-negative. Ignored when `independence = "hsic"`.
#' @param prior_knowledge Optional prior knowledge as variable pairs: a
#'   2-column matrix or a list of length-2 vectors, where each pair
#'   `c(i, j)` means "variable i cannot be a cause of variable j".
#'   Indices are **1-based** column positions (the Python implementation
#'   uses 0-based pairs).
#' @param regressor Nonlinear regressor used for all internal regressions,
#'   same interface as [lingam_resit()]: either the string `"gam"`
#'   (default; requires the suggested package mgcv) or a function
#'   `function(X, y)` returning the fitted values as a numeric vector of
#'   length `nrow(X)`. The Python implementation hardcodes pygam's
#'   `LinearGAM`, whose spline basis differs from mgcv's, so numerical
#'   agreement with Python is not expected (structural results are).
#' @return An object of class `CAMUVResult` with elements:
#' * `adjacency_matrix`: (p x p) matrix. **Convention: `B[i, j] = 1` means
#'   an edge j -> i (row = to, col = from)**, same as [lingam_direct()].
#'   Entries are 0/1 edge indicators, not coefficients (the causal
#'   functions are nonlinear). Both entries of a variable pair suspected to
#'   be connected by a UCP or UBP are `NA`.
#' * `parents_list`: a list of length `n_features`; element `i` is the
#'   sorted integer vector of the identified direct parents of variable `i`
#'   (possibly empty). Like [lingam_rcd()], there is no `causal_order`.
#' * `confounded_pairs`: 2-column integer matrix of the variable pairs
#'   (1-based column positions) suspected to be connected by a UCP or UBP;
#'   0 rows if none.
#' * `regressor`: label of the regressor used (`"gam"` or
#'   `"user function"`).
#' @details
#' The three stages mirror the Python implementation: (1)
#' `camuv_get_neighborhoods()` records which raw variable pairs are
#' dependent; (2) `camuv_find_parents()` scans variable subsets of size 2 up
#' to `num_explanatory_vals`, identifying the most sink-like variable of
#' each subset and re-testing until no new parent is found; (3) remaining
#' dependent pairs without an identified edge are flagged as UCP/UBP pairs.
#'
#' Every subset test involves GAM regressions and an HSIC test on `n x n`
#' Gram matrices, so the method is not recommended for `nrow(X)` in the
#' thousands, and cost grows combinatorially with `num_explanatory_vals`.
#'
#' `independence = "fcorr"` is only supported with
#' `num_explanatory_vals = 2` (the default): larger subsets require testing
#' a residual against several variables jointly, and the F-correlation is
#' defined for univariate pairs only (the Python implementation breaks on
#' that combination as well). The HSIC path supports any
#' `num_explanatory_vals` via the multivariate kernel.
#'
#' At least 6 observations are required (the lower bound of the HSIC
#' gamma-approximation test).
#'
#' Unlike most other estimators in this package there is no bootstrap
#' variant, matching the Python implementation, and no total-effect
#' estimation (effects are nonlinear).
#' @references
#' Maeda, T. N. and Shimizu, S. (2021). Causal additive models with
#' unobserved variables. In Proc. Thirty-Seventh Conference on Uncertainty
#' in Artificial Intelligence (UAI), PMLR 161: 97-106.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   confounded <- generate_camuv_sample(n = 200, seed = 1)
#'   result <- lingam_camuv(confounded$data)
#'   print(result)
#'
#'   # Pairs connected through unobserved variables are left NA
#'   result$confounded_pairs
#' }
#' }
lingam_camuv <- function(X,
                         alpha = 0.01,
                         num_explanatory_vals = 2L,
                         independence = "hsic",
                         ind_corr = 0.5,
                         prior_knowledge = NULL,
                         regressor = "gam") {
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
  num_explanatory_vals <- suppressWarnings(as.integer(num_explanatory_vals))
  if (length(num_explanatory_vals) != 1 || is.na(num_explanatory_vals) ||
        num_explanatory_vals < 1) {
    stop("num_explanatory_vals must be an integer >= 1.", call. = FALSE)
  }
  independence <- match.arg(independence, c("hsic", "fcorr"))
  if (!is.numeric(ind_corr) || length(ind_corr) != 1 || is.na(ind_corr) ||
        ind_corr < 0) {
    stop("ind_corr must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (independence == "fcorr" && num_explanatory_vals > 2) {
    stop(
      "independence = \"fcorr\" only supports num_explanatory_vals = 2: ",
      "larger subsets require a multivariate independence test, and the ",
      "F-correlation is defined for univariate pairs only. Use ",
      "independence = \"hsic\" instead.",
      call. = FALSE
    )
  }

  pk_forbidden <- camuv_make_pk_dict(prior_knowledge, ncol(X))
  reg <- resit_make_regressor(regressor)
  reg_fn <- reg$fn

  d <- ncol(X)
  N <- camuv_get_neighborhoods(X, independence, alpha, ind_corr)
  P <- camuv_find_parents(X, num_explanatory_vals, N, pk_forbidden,
                          independence, alpha, ind_corr, reg_fn)

  # Remaining dependent pairs with no identified edge: UCP / UBP candidates.
  # Each variable's parent-adjusted residual is pair-independent, so compute
  # it (lazily) once instead of once per pair.
  U <- list()
  pair_res <- vector("list", d)
  pair_residual <- function(v) {
    if (is.null(pair_res[[v]])) {
      pair_res[[v]] <<- camuv_get_residual(X, v, P[[v]], reg_fn)
    }
    pair_res[[v]]
  }
  for (i in seq_len(d - 1L)) {
    for (j in seq.int(i + 1L, d)) {
      if (i %in% P[[j]] || j %in% P[[i]]) next
      if (!(j %in% N[[i]]) || !(i %in% N[[j]])) next

      if (!camuv_is_independent(pair_residual(i), pair_residual(j),
                                independence, alpha, ind_corr)) {
        U[[length(U) + 1L]] <- c(i, j)
      }
    }
  }

  B <- camuv_build_adjacency_matrix(X, P, U)

  P <- lapply(P, sort)
  names(P) <- colnames(B)

  U_mat <- if (length(U) > 0) {
    do.call(rbind, U)
  } else {
    matrix(integer(0), nrow = 0, ncol = 2)
  }
  colnames(U_mat) <- c("var1", "var2")

  result <- list(
    adjacency_matrix = B,
    parents_list = P,
    confounded_pairs = U_mat,
    regressor = reg$label
  )
  class(result) <- "CAMUVResult"
  result
}


#' Print method for CAMUVResult
#'
#' @param x CAMUVResult object
#' @param digits Number of digits to display
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   confounded <- generate_camuv_sample(n = 200, seed = 1)
#'   result <- lingam_camuv(confounded$data)
#'   print(result)
#' }
#' }
print.CAMUVResult <- function(x, digits = 3, ...) {
  var_names <- colnames(x$adjacency_matrix)
  label_for <- function(idx) {
    if (!is.null(var_names)) var_names[idx] else paste0("x", idx - 1L)
  }

  cat("CAM-UV Result\n")
  cat(sprintf("  Variables : %d\n", ncol(x$adjacency_matrix)))
  cat(sprintf("  Regressor : %s\n", x$regressor))
  cat("\nParent sets:\n")
  for (i in seq_along(x$parents_list)) {
    pa <- x$parents_list[[i]]
    pa_str <- if (length(pa) == 0) {
      "{}"
    } else {
      paste0("{", paste(label_for(pa), collapse = ", "), "}")
    }
    cat(sprintf("  P(%s) = %s\n", label_for(i), pa_str))
  }
  if (nrow(x$confounded_pairs) > 0) {
    cat("\nPairs with an unobserved causal/backdoor path (UCP/UBP):\n")
    for (r in seq_len(nrow(x$confounded_pairs))) {
      cat(sprintf(
        "  %s -- %s\n",
        label_for(x$confounded_pairs[r, 1]),
        label_for(x$confounded_pairs[r, 2])
      ))
    }
  } else {
    cat("\nNo pairs with an unobserved causal/backdoor path (UCP/UBP).\n")
  }
  cat("\nAdjacency matrix (row = to, col = from):\n")
  cat("  (entries are 0/1 edge indicators, not coefficients;\n")
  cat("   NA = pair connected through an unobserved variable)\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}
