# =============================================================================
# Bootstrap for VAR-LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam  (lingam/var_lingam.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Bootstrap for VAR-LiNGAM
#'
#' Evaluates the statistical reliability of the estimated time-series DAG by
#' resampling. Unlike the i.i.d. row resampling used for Direct LiNGAM, this
#' uses a **residual bootstrap**: the VAR is fitted once on the original data,
#' the residuals are resampled with replacement, and a new series is rebuilt by
#' the VAR recursion before re-estimating VAR-LiNGAM on it (this preserves the
#' autoregressive structure). Port of the Python reference `VARLiNGAM.bootstrap`.
#'
#' @param X numeric matrix or data frame (n_samples x n_features), rows ordered
#'   in time.
#' @param n_sampling number of bootstrap iterations (positive integer).
#' @param lags maximum lag order. When `criterion` is not NULL, the lag is
#'   selected once on the original data and then fixed across all iterations.
#' @param criterion lag-selection criterion ("bic", "aic", "hqic", "fpe") or
#'   NULL to use `lags` directly.
#' @param measure independence measure for [lingam_direct()] ("pwling"/"kernel").
#' @param reg_method regression method for the instantaneous matrix.
#' @param lambda penalty selection (see [lingam_direct()]).
#' @param init_method initial-weight method for adaptive LASSO.
#' @param prune logical; passed to [lingam_var()] on each iteration (default TRUE).
#' @param seed random seed (NULL allowed).
#' @param verbose whether to print progress (logical).
#' @param parallel whether to distribute iterations across cores (logical).
#' @param n_cores number of cores (integer or NULL; NULL caps at 2 for safety).
#' @return a `VARBootstrapResult` object.
#' @details
#' Reproducibility follows the same rules as [lingam_direct_bootstrap()]: with
#' `parallel = TRUE`, L'Ecuyer streams via `parallel::clusterSetRNGStream()` make
#' results reproducible for a given `seed` and `n_cores`, but they do not match
#' the sequential (`parallel = FALSE`) results.
#' @importFrom stats median
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 500, seed = 42)
#'
#' # Fast example: OLS instantaneous structure, no pruning (no glmnet needed)
#' bs <- lingam_var_bootstrap(s$data,
#'   n_sampling = 10L, lags = 1, criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_var_probabilities(bs)
lingam_var_bootstrap <- function(X,
                                 n_sampling,
                                 lags = 1L,
                                 criterion = "bic",
                                 measure = "pwling",
                                 reg_method = "adaptive_lasso",
                                 lambda = "BIC",
                                 init_method = "ols",
                                 prune = TRUE,
                                 seed = NULL,
                                 verbose = TRUE,
                                 parallel = FALSE,
                                 n_cores = NULL) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)

  # Validate up-front (otherwise errors would surface confusingly inside workers).
  measure <- match.arg(measure, c("pwling", "kernel"))
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))
  if (reg_method == "ridge" && lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }
  if (!is.logical(prune) || length(prune) != 1 || is.na(prune)) {
    stop("prune must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  lags <- suppressWarnings(as.integer(lags))
  if (length(lags) != 1 || is.na(lags) || lags < 1) {
    stop("lags must be a positive integer.", call. = FALSE)
  }

  # Select the lag order once on the original data, then fix it for all
  # iterations (mirrors the Python reference, which disables selection inside
  # the bootstrap loop).
  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("bic", "aic", "hqic", "fpe"))
    lags <- select_var_lag(X, max_lag = lags, criterion = criterion)
  }

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Pre-fit the VAR on the original data: M (ar coefs) and residuals drive the
  # residual bootstrap below.
  vf <- fit_var_ols(X, lags)
  M <- vf$coefs          # array (lags, n_features, n_features)
  residuals <- vf$residuals
  n_resid <- nrow(residuals)

  # One bootstrap iteration: residual resample -> VAR recursion -> re-estimate.
  run_one <- function(i) {
    # i.i.d. resample of residual rows, up to the original series length.
    ridx <- sample.int(n_resid, n_samples, replace = TRUE)
    sampled <- residuals[ridx, , drop = FALSE]

    resampled_X <- matrix(0, nrow = n_samples, ncol = n_features)
    for (j in seq_len(n_samples)) {
      if (j <= lags) {
        # seed the first `lags` rows with the resampled noise
        resampled_X[j, ] <- sampled[j, ]
      } else {
        ar <- numeric(n_features)
        for (k in seq_len(lags)) {
          ar <- ar + as.numeric(M[k, , ] %*% resampled_X[j - k, ])
        }
        resampled_X[j, ] <- ar + sampled[j, ]
      }
    }

    res <- lingam_var(resampled_X,
      lags = lags, criterion = NULL,
      measure = measure, reg_method = reg_method,
      lambda = lambda, init_method = init_method, prune = prune
    )
    am <- res$adjacency_matrices
    causal_order <- res$causal_order

    # joined adjacency cbind(B0, B1, ..., Bp): n_features x n_features*(1+lags)
    am_joined <- do.call(cbind, lapply(seq_len(lags + 1L), function(k) am[k, , ]))

    # total effects over the time-expanded graph (square-padded joined matrix)
    ncol_j <- ncol(am_joined)
    am_sq <- matrix(0, nrow = ncol_j, ncol = ncol_j)
    am_sq[seq_len(n_features), ] <- am_joined

    te <- matrix(0, nrow = n_features, ncol = n_features * (lags + 1L))
    for (ci in seq_len(n_features)) {
      to <- rev(causal_order)[ci]
      # contemporaneous sources: those preceding `to` in the causal order
      n_earlier <- n_features - ci
      if (n_earlier >= 1L) {
        for (from in causal_order[seq_len(n_earlier)]) {
          te[to, from] <- calculate_total_effect(am_sq, from, to)
        }
      }
      # lagged sources: all variables at each lag
      for (lag in seq_len(lags)) {
        for (from in seq_len(n_features)) {
          from_col <- from + n_features * lag
          te[to, from_col] <- calculate_total_effect(am_sq, from_col, to)
        }
      }
    }

    list(adjacency = am_joined, total_effects = te, idx = ridx, causal_order = causal_order)
  }

  # Resolve cores (default cap of 2, matching lingam_direct_bootstrap).
  if (parallel) {
    available <- parallel::detectCores()
    if (is.na(available)) available <- 1L
    if (is.null(n_cores)) {
      n_cores <- min(2L, available)
    } else {
      n_cores <- as.integer(n_cores)
      if (is.na(n_cores) || n_cores < 1L) stop("n_cores must be a positive integer.")
    }
    n_cores <- max(1L, min(n_cores, available, n_sampling))
    if (n_cores == 1L) parallel <- FALSE
  }

  if (verbose) {
    mode_str <- if (parallel) sprintf("parallel, %d cores", n_cores) else "sequential"
    message(sprintf(
      "VAR-LiNGAM bootstrap: %d iterations, lag=%d, method=%s (%s)",
      n_sampling, lags, reg_method, mode_str
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Make this package available on the workers (same approach as the Direct
    # LiNGAM bootstrap: attach when installed, otherwise export the namespace).
    pkg <- utils::packageName()
    ns <- environment(lingam_var)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    worker_ready <- FALSE
    if (!is.null(pkg) && nzchar(pkg)) {
      worker_ready <- all(unlist(parallel::clusterCall(cl, function(p) {
        suppressWarnings(suppressMessages(isTRUE(requireNamespace(p, quietly = TRUE))))
      }, pkg)))
      if (worker_ready) {
        parallel::clusterCall(cl, function(p) {
          suppressWarnings(suppressMessages(library(p, character.only = TRUE)))
          invisible(NULL)
        }, pkg)
      }
    }
    if (!worker_ready) {
      parallel::clusterExport(cl, varlist = ls(ns, all.names = TRUE), envir = ns)
    }

    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    res_list <- parallel::parLapply(cl, seq_len(n_sampling), run_one)
  } else {
    if (!is.null(seed)) set.seed(seed)
    res_list <- lapply(seq_len(n_sampling), function(i) {
      if (verbose && (i %% 10 == 0 || i == 1)) {
        message(sprintf("  iteration %d / %d", i, n_sampling))
      }
      run_one(i)
    })
  }

  # Aggregate.
  adjacency_matrices <- vector("list", n_sampling)
  total_effects <- array(0, dim = c(n_sampling, n_features, n_features * (lags + 1L)))
  causal_orders <- matrix(0L, nrow = n_sampling, ncol = n_features)
  resampled_indices <- vector("list", n_sampling)
  for (i in seq_len(n_sampling)) {
    adjacency_matrices[[i]] <- res_list[[i]]$adjacency
    total_effects[i, , ] <- res_list[[i]]$total_effects
    causal_orders[i, ] <- res_list[[i]]$causal_order
    resampled_indices[[i]] <- res_list[[i]]$idx
  }

  if (verbose) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    message(sprintf("Completed in %.1f seconds.", elapsed))
  }

  create_var_bootstrap_result(
    adjacency_matrices, total_effects, lags, resampled_indices, causal_orders
  )
}


# =============================================================================
# VARBootstrapResult object
# =============================================================================

#' Create a VARBootstrapResult
#'
#' @param adjacency_matrices list (length n_sampling); each element is a joined
#'   adjacency matrix (n_features x n_features*(1 + lags))
#' @param total_effects array (n_sampling x n_features x n_features*(1 + lags))
#' @param lags lag order used
#' @param resampled_indices list of residual-index vectors (NULL allowed)
#' @param causal_orders matrix (n_sampling x n_features) (NULL allowed)
#' @return a VARBootstrapResult (list with class attribute)
#' @keywords internal
create_var_bootstrap_result <- function(adjacency_matrices, total_effects, lags,
                                        resampled_indices = NULL, causal_orders = NULL) {
  obj <- list(
    adjacency_matrices = adjacency_matrices,
    total_effects      = total_effects,
    lags               = lags,
    resampled_indices  = resampled_indices,
    causal_orders      = causal_orders
  )
  class(obj) <- "VARBootstrapResult"
  obj
}


#' Print a VARBootstrapResult
#'
#' @param x a VARBootstrapResult object
#' @param ... additional arguments (unused)
#' @method print VARBootstrapResult
#' @export
print.VARBootstrapResult <- function(x, ...) {
  n_sampling <- length(x$adjacency_matrices)
  n_features <- nrow(x$adjacency_matrices[[1]])
  cat(sprintf(
    "VARBootstrapResult: %d samplings, %d features, lag order %d\n",
    n_sampling, n_features, x$lags
  ))
  invisible(x)
}


# =============================================================================
# VARBootstrapResult methods
# =============================================================================

#' Bootstrap probabilities for a VAR-LiNGAM model
#'
#' Returns, for each entry of the joined adjacency matrix, the fraction of
#' bootstrap samples in which that edge exceeded `min_causal_effect`.
#'
#' @param result a VARBootstrapResult object
#' @param min_causal_effect minimum |effect| threshold (NULL = 0)
#' @return probability matrix (n_features x n_features*(1 + lags)). Columns
#'   1..n_features are the instantaneous block; the next n_features are lag 1; etc.
#'   `P[i, j]` is the probability of the edge j -> i.
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 500, seed = 42)
#' bs <- lingam_var_bootstrap(s$data,
#'   n_sampling = 10L, criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_var_probabilities(bs)
get_var_probabilities <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "VARBootstrapResult"))
  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.", call. = FALSE)

  ams <- result$adjacency_matrices
  n_sampling <- length(ams)
  acc <- matrix(0, nrow = nrow(ams[[1]]), ncol = ncol(ams[[1]]))
  for (am in ams) {
    am[is.na(am)] <- 0
    acc <- acc + (abs(am) > min_causal_effect)
  }
  acc / n_sampling
}


#' Enumerate bootstrap paths between two variables in a VAR-LiNGAM model
#'
#' Builds the time-expanded graph for every bootstrap sample and enumerates all
#' directed paths from the source (at `from_lag`) to the destination (at
#' `to_lag`), reporting each path's bootstrap probability and median effect.
#' Port of the Python reference `VARBootstrapResult.get_paths`.
#'
#' Node indices in the returned `path` are 1-based positions in the time-expanded
#' graph: column j of block L (lag L) corresponds to index `n_features * L + j`.
#'
#' @param result a VARBootstrapResult object
#' @param from_index source variable (1-based)
#' @param to_index destination variable (1-based)
#' @param from_lag lag of the source (default 0)
#' @param to_lag lag of the destination (default 0); must satisfy `to_lag <= from_lag`
#' @param min_causal_effect minimum |effect| threshold (NULL = 0)
#' @return a data frame (path, effect, probability), one row per distinct path
#' @importFrom stats median
#' @export
#' @examples
#' s <- generate_varlingam_sample(n = 500, seed = 42)
#' bs <- lingam_var_bootstrap(s$data,
#'   n_sampling = 10L, criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_var_paths(bs, from_index = 1, to_index = 3)
get_var_paths <- function(result, from_index, to_index,
                          from_lag = 0, to_lag = 0, min_causal_effect = NULL) {
  stopifnot(inherits(result, "VARBootstrapResult"))
  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.", call. = FALSE)

  from_lag <- as.integer(from_lag)
  to_lag <- as.integer(to_lag)
  if (is.na(from_lag) || is.na(to_lag) || from_lag < 0 || to_lag < 0) {
    stop("from_lag and to_lag must be non-negative integers.", call. = FALSE)
  }
  if (to_lag > from_lag) stop("from_lag must be >= to_lag.", call. = FALSE)
  if (to_lag == from_lag && to_index == from_index) {
    stop("from_index and to_index refer to the same variable.", call. = FALSE)
  }

  ams <- result$adjacency_matrices
  n_sampling <- length(ams)
  nf <- nrow(ams[[1]])
  n_lags <- ncol(ams[[1]]) / nf - 1L

  paths_collector <- vector("list", n_sampling)
  effects_collector <- vector("list", n_sampling)
  for (s_idx in seq_along(ams)) {
    am <- ams[[s_idx]]
    dim_e <- ncol(am)
    # Time-expanded square graph: block (i, j) for j >= i holds B_{j-i}, so an
    # earlier-time node (larger lag block) points to a later-time node.
    expansion_m <- matrix(0, nrow = dim_e, ncol = dim_e)
    for (i in 0:n_lags) {
      for (j in i:n_lags) {
        row <- nf * i
        col <- nf * j
        lag <- col - row
        expansion_m[(row + 1L):(row + nf), (col + 1L):(col + nf)] <-
          am[, (lag + 1L):(lag + nf), drop = FALSE]
      }
    }
    fr <- nf * from_lag + from_index
    to <- nf * to_lag + to_index
    res <- find_all_paths(expansion_m, fr, to, min_causal_effect)
    if (length(res$paths) > 0) {
      paths_collector[[s_idx]] <- vapply(res$paths, paste, "", collapse = "_")
      effects_collector[[s_idx]] <- res$effects
    }
  }
  paths_all <- unlist(paths_collector, use.names = FALSE)
  effects_all <- unlist(effects_collector, use.names = FALSE)

  if (length(paths_all) == 0) {
    return(data.frame(path = character(0), effect = numeric(0), probability = numeric(0)))
  }

  tbl <- sort(table(paths_all), decreasing = TRUE)
  path_strs <- names(tbl)
  probs <- as.numeric(tbl) / n_sampling
  effects_median <- vapply(
    path_strs,
    function(ps) stats::median(effects_all[paths_all == ps]),
    numeric(1)
  )
  path_list <- lapply(path_strs, function(ps) as.integer(strsplit(ps, "_")[[1]]))

  data.frame(
    path        = I(path_list),
    effect      = as.numeric(effects_median),
    probability = probs,
    row.names   = NULL
  )
}
