# =============================================================================
# Bootstrap for VARMA-LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam  (lingam/varma_lingam.py)
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Bootstrap for VARMA-LiNGAM
#'
#' Evaluates the statistical reliability of the estimated time-series DAG by
#' resampling. Like [lingam_var_bootstrap()], this uses a **residual
#' bootstrap**: the VARMA model is fitted once on the original data, the
#' residuals are resampled with replacement, and a new series is rebuilt by
#' the VARMA recursion before re-estimating VARMA-LiNGAM on it. Port of the
#' Python reference `VARMALiNGAM.bootstrap`.
#'
#' @param X numeric matrix or data frame (n_samples x n_features), rows ordered
#'   in time.
#' @param n_sampling number of bootstrap iterations (positive integer).
#' @param order VARMA order `c(p, q)`. When `criterion` is not NULL, the order
#'   is selected once on the original data and then fixed across all iterations.
#' @param criterion order-selection criterion ("bic", "aic", "hqic") or NULL to
#'   use `order` directly.
#' @param measure independence measure for [lingam_direct()] ("pwling"/"kernel").
#' @param reg_method regression method for the instantaneous matrix.
#' @param lambda penalty selection (see [lingam_direct()]).
#' @param init_method initial-weight method for adaptive LASSO.
#' @param prune logical; passed to [lingam_varma()] on each iteration (default TRUE).
#' @param seed random seed (NULL allowed).
#' @param verbose whether to print progress (logical).
#' @param parallel whether to distribute iterations across cores (logical).
#' @param n_cores number of cores (integer or NULL; NULL caps at 2 for safety).
#' @return a `VARMABootstrapResult` object.
#' @details
#' Reproducibility follows the same rules as [lingam_direct_bootstrap()]: with
#' `parallel = TRUE`, L'Ecuyer streams via `parallel::clusterSetRNGStream()` make
#' results reproducible for a given `seed` and `n_cores`, but they do not match
#' the sequential (`parallel = FALSE`) results.
#'
#' **On iteration failures:** as in [lingam_direct_bootstrap()], each iteration
#' runs inside a `tryCatch()`; a failing iteration is reported as a warning and
#' excluded from the result instead of aborting the run. An error is raised
#' only if every iteration fails.
#'
#' As in the Python reference, the series regeneration omits the estimated
#' intercept, so resampled series are centered near zero even when the
#' original data are not; each refit re-estimates its own intercept, so the
#' resampled coefficient estimates are unaffected.
#'
#' Total effects are estimated by the back-door regression of
#' [estimate_varma_total_effect()] (the Python reference does the same) and
#' cover the instantaneous block and the AR lags 1..p; the MA (omega) blocks
#' describe effects of past disturbances, not of observed variables, and are
#' therefore excluded from `total_effects`.
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 300, seed = 42)
#'
#' # Fast example: OLS instantaneous structure, no pruning (no glmnet needed)
#' bs <- lingam_varma_bootstrap(s$data,
#'   n_sampling = 5L, order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_varma_probabilities(bs)
lingam_varma_bootstrap <- function(X,
                                   n_sampling,
                                   order = c(1L, 1L),
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
  reg_args <- validate_reg_args(reg_method, lambda, init_method)
  reg_method <- reg_args$reg_method
  lambda <- reg_args$lambda
  init_method <- reg_args$init_method
  if (!is.logical(prune) || length(prune) != 1 || is.na(prune)) {
    stop("prune must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  order <- suppressWarnings(as.integer(order))
  if (length(order) != 2 || anyNA(order) || any(order < 0) || sum(order) < 1) {
    stop("order must be two non-negative integers c(p, q), not both zero.", call. = FALSE)
  }

  # Select the order once on the original data, then fix it for all iterations
  # (mirrors the Python reference, which disables selection inside the loop).
  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("bic", "aic", "hqic"))
    order <- select_varma_order(X,
      max_p = order[1], max_q = order[2],
      criterion = criterion
    )
  }
  p_order <- order[1]
  q_order <- order[2]
  k0 <- max(p_order, q_order)

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Pre-fit the VARMA on the original data: Phi/Theta and the filtered
  # residuals drive the residual bootstrap below.
  hr <- fit_varma_hr(X, p_order, q_order)
  phis <- hr$phis
  thetas <- hr$thetas
  E_full <- filter_varma_residuals(X, phis, thetas, hr$const)
  residuals <- E_full[(k0 + 1L):n_samples, , drop = FALSE]
  n_resid <- nrow(residuals)
  # slice the coefficient arrays once; 3D indexing inside the recursion below
  # would copy phis[tau, , ] / thetas[w, , ] on every time step
  phi_list <- lapply(seq_len(p_order), function(tau) phis[tau, , ])
  theta_list <- lapply(seq_len(q_order), function(w) thetas[w, , ])

  # One bootstrap iteration: residual resample -> VARMA recursion -> re-estimate.
  # Wrapped in tryCatch (like lingam_direct_bootstrap) so that one pathological
  # resample does not abort the entire run; failed iterations are reported as
  # warnings and excluded from the aggregated result.
  run_one <- function(i) {
    tryCatch({
      # i.i.d. resample of residual rows, up to the original series length.
      ridx <- sample.int(n_resid, n_samples, replace = TRUE)
      sampled <- residuals[ridx, , drop = FALSE]

      resampled_X <- matrix(0, nrow = n_samples, ncol = n_features)
      for (j in seq_len(n_samples)) {
        if (j <= k0) {
          # seed the first max(p, q) rows with the resampled noise
          resampled_X[j, ] <- sampled[j, ]
        } else {
          pred <- numeric(n_features)
          for (tau in seq_len(p_order)) {
            pred <- pred + as.numeric(phi_list[[tau]] %*% resampled_X[j - tau, ])
          }
          for (w in seq_len(q_order)) {
            pred <- pred + as.numeric(theta_list[[w]] %*% sampled[j - w, ])
          }
          resampled_X[j, ] <- pred + sampled[j, ]
        }
      }

      res <- lingam_varma(resampled_X,
        order = order, criterion = NULL,
        measure = measure, reg_method = reg_method,
        lambda = lambda, init_method = init_method, prune = prune
      )
      psis <- res$adjacency_matrices$psis
      omegas <- res$adjacency_matrices$omegas
      causal_order <- res$causal_order

      # joined matrix cbind(psi_0..psi_p, omega_1..omega_q):
      # n_features x n_features*(1 + p + q)
      am_joined <- joined_varma_matrix(psis, omegas)

      # LiNGAM residuals of this fit, full length, for the total-effect core.
      E_boot <- matrix(0, n_samples, n_features)
      E_boot[(k0 + 1L):n_samples, ] <- res$residuals
      ee_full <- E_boot %*% t(diag(n_features) - psis[1, , ])

      # total effects by back-door regression (same as the Python reference);
      # instantaneous block plus AR lags only. The joined design depends only
      # on the source lag, so build it once per lag, not once per pair.
      designs <- lapply(0L:p_order, function(fl) {
        varma_joined_design(resampled_X, ee_full, order, fl)
      })
      te <- matrix(0, nrow = n_features, ncol = n_features * (1L + p_order))
      for (ci in seq_len(n_features)) {
        to <- rev(causal_order)[ci]
        # contemporaneous sources: those preceding `to` in the causal order
        n_earlier <- n_features - ci
        if (n_earlier >= 1L) {
          for (from in causal_order[seq_len(n_earlier)]) {
            te[to, from] <- varma_total_effect_core(
              resampled_X, ee_full, am_joined, order, from, to, 0L,
              X_joined = designs[[1L]]
            )
          }
        }
        # lagged sources: all variables at each AR lag
        for (lag in seq_len(p_order)) {
          for (from in seq_len(n_features)) {
            te[to, from + n_features * lag] <- varma_total_effect_core(
              resampled_X, ee_full, am_joined, order, from, to, lag,
              X_joined = designs[[lag + 1L]]
            )
          }
        }
      }

      list(ok = TRUE, adjacency = am_joined, total_effects = te, idx = ridx,
           causal_order = causal_order)
    }, error = function(e) {
      list(ok = FALSE, iteration = i, message = conditionMessage(e))
    })
  }

  cores <- resolve_bootstrap_cores(parallel, n_cores, n_sampling)
  parallel <- cores$parallel
  n_cores <- cores$n_cores

  if (verbose) {
    message(sprintf(
      "VARMA-LiNGAM bootstrap: %d iterations, order=(%d,%d), method=%s (%s)",
      n_sampling, p_order, q_order, reg_method,
      bootstrap_mode_string(parallel, n_cores)
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Make this package available on the workers (same approach as the Direct
    # LiNGAM bootstrap: attach when installed, otherwise export the namespace).
    setup_cluster_worker(cl, lingam_varma)

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

  # Report and drop any failed iterations before aggregating (same failure
  # policy as lingam_direct_bootstrap).
  res_list <- filter_bootstrap_failures(res_list)
  n_success <- length(res_list)

  # Aggregate.
  adjacency_matrices <- vector("list", n_success)
  total_effects <- array(0, dim = c(n_success, n_features, n_features * (1L + p_order)))
  causal_orders <- matrix(0L, nrow = n_success, ncol = n_features)
  resampled_indices <- vector("list", n_success)
  for (i in seq_len(n_success)) {
    adjacency_matrices[[i]] <- res_list[[i]]$adjacency
    total_effects[i, , ] <- res_list[[i]]$total_effects
    causal_orders[i, ] <- res_list[[i]]$causal_order
    resampled_indices[[i]] <- res_list[[i]]$idx
  }

  if (verbose) bootstrap_completion_message(t_start, n_success, n_sampling)

  create_varma_bootstrap_result(
    adjacency_matrices, total_effects, order, resampled_indices, causal_orders
  )
}


# =============================================================================
# VARMABootstrapResult object
# =============================================================================

#' Create a VARMABootstrapResult
#'
#' @param adjacency_matrices list (length n_sampling); each element is a joined
#'   matrix `cbind(psi_0..psi_p, omega_1..omega_q)`
#'   (n_features x n_features*(1 + p + q))
#' @param total_effects array (n_sampling x n_features x n_features*(1 + p))
#' @param order VARMA order c(p, q) used
#' @param resampled_indices list of residual-index vectors (NULL allowed)
#' @param causal_orders matrix (n_sampling x n_features) (NULL allowed)
#' @return a VARMABootstrapResult (list with class attribute)
#' @keywords internal
create_varma_bootstrap_result <- function(adjacency_matrices, total_effects, order,
                                          resampled_indices = NULL, causal_orders = NULL) {
  obj <- list(
    adjacency_matrices = adjacency_matrices,
    total_effects      = total_effects,
    order              = order,
    resampled_indices  = resampled_indices,
    causal_orders      = causal_orders
  )
  class(obj) <- "VARMABootstrapResult"
  obj
}


#' Print a VARMABootstrapResult
#'
#' @param x a VARMABootstrapResult object
#' @param ... additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @method print VARMABootstrapResult
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 300, seed = 42)
#' bs <- lingam_varma_bootstrap(s$data,
#'   n_sampling = 5L, order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' print(bs)
print.VARMABootstrapResult <- function(x, ...) {
  n_sampling <- length(x$adjacency_matrices)
  n_features <- nrow(x$adjacency_matrices[[1]])
  cat(sprintf(
    "VARMABootstrapResult: %d samplings, %d features, order (%d, %d)\n",
    n_sampling, n_features, x$order[1], x$order[2]
  ))
  invisible(x)
}


# =============================================================================
# VARMABootstrapResult methods
# =============================================================================

#' Bootstrap probabilities for a VARMA-LiNGAM model
#'
#' Returns, for each entry of the joined matrix, the fraction of bootstrap
#' samples in which that edge exceeded `min_causal_effect`.
#'
#' @param result a VARMABootstrapResult object
#' @param min_causal_effect minimum |effect| threshold (NULL = 0)
#' @return probability matrix (n_features x n_features*(1 + p + q)). Columns
#'   1..n_features are the instantaneous block, the next p blocks are the AR
#'   lags 1..p (psi), and the final q blocks are the MA terms 1..q (omega).
#'   `P[i, j]` is the probability of the edge j -> i.
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 300, seed = 42)
#' bs <- lingam_varma_bootstrap(s$data,
#'   n_sampling = 5L, order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_varma_probabilities(bs)
get_varma_probabilities <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "VARMABootstrapResult"))
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


#' Enumerate bootstrap paths between two variables in a VARMA-LiNGAM model
#'
#' Builds the time-expanded graph for every bootstrap sample and enumerates all
#' directed paths from the source (at `from_lag`) to the destination (at
#' `to_lag`), reporting each path's bootstrap probability and median effect.
#' Port of the Python reference `VARMABootstrapResult.get_paths`.
#'
#' Node indices in the returned `path` are 1-based positions in the time-expanded
#' graph: column j of block L (lag L) corresponds to index `n_features * L + j`.
#'
#' @param result a VARMABootstrapResult object
#' @param from_index source variable (1-based)
#' @param to_index destination variable (1-based)
#' @param from_lag lag of the source (default 0); must not exceed the AR order p
#' @param to_lag lag of the destination (default 0); must satisfy `to_lag <= from_lag`
#' @param min_causal_effect minimum |effect| threshold (NULL = 0)
#' @return a data frame (path, effect, probability), one row per distinct path
#' @details
#' Only the instantaneous and AR (psi) blocks enter the time-expanded graph;
#' the MA (omega) blocks describe effects of past unobserved disturbances,
#' which are not nodes of the variable graph (the Python reference does the
#' same).
#' @export
#' @examples
#' s <- generate_varmalingam_sample(n = 300, seed = 42)
#' bs <- lingam_varma_bootstrap(s$data,
#'   n_sampling = 5L, order = c(1, 1), criterion = NULL,
#'   reg_method = "ols", prune = FALSE, seed = 1, verbose = FALSE
#' )
#' get_varma_paths(bs, from_index = 1, to_index = 3)
get_varma_paths <- function(result, from_index, to_index,
                            from_lag = 0, to_lag = 0, min_causal_effect = NULL) {
  stopifnot(inherits(result, "VARMABootstrapResult"))
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
  nf <- nrow(ams[[1]])
  n_lags <- result$order[1]
  # Lags beyond the AR order would index past the time-expanded graph and
  # surface as a bare "subscript out of bounds" error below.
  if (from_lag > n_lags || to_lag > n_lags) {
    stop(sprintf(
      "from_lag and to_lag must not exceed the model's AR order (%d).", n_lags
    ), call. = FALSE)
  }

  # Keep only the instantaneous and AR (psi) blocks; drop the omega blocks.
  psi_ams <- lapply(ams, function(am) am[, seq_len(nf * (1L + n_lags)), drop = FALSE])
  collect_time_expanded_paths(
    psi_ams, nf, n_lags, from_index, to_index, from_lag, to_lag, min_causal_effect
  )
}
