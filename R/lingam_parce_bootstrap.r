# =============================================================================
# Bootstrap for Bottom-Up ParceLiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/bottom_up_parce_lingam.py, 572-618 lines)
#
# Structure (worker setup, RNG streams, failure handling) follows
# lingam_direct_bootstrap() (R/lingam_bootstrap.r).
# =============================================================================


#' Bootstrap for Bottom-Up ParceLiNGAM
#'
#' @param X Numeric matrix (n_samples x n_features)
#' @param n_sampling Number of bootstrap iterations
#' @param prior_knowledge Prior knowledge matrix (NULL allowed)
#' @param alpha Significance level, passed to [lingam_parce()]
#' @param independence Independence measure, passed to [lingam_parce()]
#' @param ind_corr F-correlation rejection threshold, passed to [lingam_parce()]
#' @param reg_method Regression method ("ols", "lasso", "adaptive_lasso", "ridge")
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")
#' @param init_method Method for estimating the initial weights of adaptive LASSO
#'   regression ("ols" or "ridge")
#' @param seed Random seed (NULL allowed)
#' @param verbose Whether to display progress (logical)
#' @param parallel Whether to use parallel processing (logical)
#' @param n_cores Number of cores to use (integer, NULL allowed)
#' @param compute_total_effects Whether to also estimate total causal effects
#'   for every variable pair on each bootstrap iteration (logical, default `TRUE`).
#' @return A `BootstrapResult` (list); see [lingam_direct_bootstrap()] for the
#'   query helpers that operate on it (`get_probabilities()`,
#'   `get_causal_direction_counts()`, `get_directed_acyclic_graph_counts()`,
#'   `get_total_causal_effects()`).
#' @details
#' **Total effects are path sums, not regression estimates.** Each
#' iteration's total-effect matrix is built from [calculate_total_effect()]
#' (summing products of adjacency-matrix coefficients along every directed
#' path), matching the upstream Python implementation's bootstrap method
#' (`estimate_total_effect2`). If a variable's row in the adjacency matrix
#' contains `NA` (it is part of an unresolved block), all of its outgoing
#' total effects are set to `NA` for that iteration, since its causal
#' parents cannot be identified.
#'
#' **`NA` (unresolved) edges are treated as absent when aggregating.** Both
#' the adjacency matrix and the total-effect matrix have `NA` replaced by
#' `0` before being stored in the returned `BootstrapResult`, matching the
#' numpy comparison semantics used by the upstream implementation (where
#' `np.abs(nan) > threshold` evaluates to `FALSE`). This means, for example,
#' `get_probabilities()` reports the confounded pair's edge probability as
#' the fraction of resamples in which the order happened to resolve, not as
#' `NA`.
#'
#' **`causal_orders` is not populated** (unlike [lingam_direct_bootstrap()]):
#' ParceLiNGAM's causal order can include an unresolved block, which does
#' not fit the fixed-length integer-vector format `causal_orders` requires.
#' As a result, [get_causal_order_stability()] cannot be used with a
#' `BootstrapResult` returned by this function.
#' @export
#' @examples
#' \donttest{
#' confounded <- generate_parce_sample(n = 500, seed = 1)
#'
#' bs <- lingam_parce_bootstrap(confounded$data,
#'   n_sampling = 10L,
#'   reg_method = "ols",
#'   seed = 42
#' )
#' get_probabilities(bs)
#' }
lingam_parce_bootstrap <- function(X,
                                   n_sampling,
                                   prior_knowledge = NULL,
                                   alpha = 0.1,
                                   independence = "hsic",
                                   ind_corr = 0.5,
                                   reg_method = "adaptive_lasso",
                                   lambda = "BIC",
                                   init_method = "ols",
                                   seed = NULL,
                                   verbose = TRUE,
                                   parallel = FALSE,
                                   n_cores = NULL,
                                   compute_total_effects = TRUE) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  if (!is.logical(compute_total_effects) || length(compute_total_effects) != 1 ||
        is.na(compute_total_effects)) {
    stop("compute_total_effects must be a single logical (TRUE or FALSE).", call. = FALSE)
  }

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
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  n_samples <- nrow(X)
  n_features <- ncol(X)

  run_one <- function(i) {
    tryCatch({
      idx <- sample(n_samples, replace = TRUE)
      resampled_X <- X[idx, , drop = FALSE]
      result <- lingam_parce(
        resampled_X,
        alpha = alpha,
        prior_knowledge = prior_knowledge,
        independence = independence,
        ind_corr = ind_corr,
        reg_method = reg_method,
        lambda = lambda,
        init_method = init_method
      )
      B <- result$adjacency_matrix

      te <- if (compute_total_effects) {
        has_na_row <- apply(B, 1, anyNA)
        TE <- matrix(0, n_features, n_features)
        for (from_idx in seq_len(n_features)) {
          for (to_idx in seq_len(n_features)) {
            if (from_idx == to_idx) next
            TE[to_idx, from_idx] <- if (has_na_row[from_idx]) {
              NA_real_
            } else {
              calculate_total_effect(B, from_idx, to_idx)
            }
          }
        }
        TE
      } else {
        NULL
      }

      # NA (unresolved) edges are treated as absent when aggregating (see @details)
      B[is.na(B)] <- 0
      if (!is.null(te)) te[is.na(te)] <- 0

      list(
        ok               = TRUE,
        idx              = idx,
        adjacency_matrix = B,
        total_effects    = te
      )
    }, error = function(e) {
      list(ok = FALSE, iteration = i, message = conditionMessage(e))
    })
  }

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
      "Bootstrap: %d iterations, method=%s (%s)",
      n_sampling, reg_method, mode_str
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    setup_cluster_worker(cl, lingam_parce)

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

  ok <- vapply(res_list, function(r) isTRUE(r$ok), logical(1))
  if (any(!ok)) {
    for (r in res_list[!ok]) {
      warning(sprintf(
        "Bootstrap iteration %d failed and was skipped: %s", r$iteration, r$message
      ), call. = FALSE)
    }
  }
  res_list <- res_list[ok]
  n_success <- length(res_list)
  if (n_success == 0) {
    stop("All bootstrap iterations failed; see warnings above for details.", call. = FALSE)
  }

  adjacency_matrices <- array(0, dim = c(n_success, n_features, n_features))
  total_effects <- if (compute_total_effects) {
    array(0, dim = c(n_success, n_features, n_features))
  } else {
    NULL
  }
  resampled_indices <- vector("list", n_success)
  for (i in seq_len(n_success)) {
    adjacency_matrices[i, , ] <- res_list[[i]]$adjacency_matrix
    if (compute_total_effects) total_effects[i, , ] <- res_list[[i]]$total_effects
    resampled_indices[[i]] <- res_list[[i]]$idx
  }

  if (verbose) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    if (n_success < n_sampling) {
      message(sprintf(
        "Completed in %.1f seconds (%d / %d iterations succeeded).",
        elapsed, n_success, n_sampling
      ))
    } else {
      message(sprintf("Completed in %.1f seconds.", elapsed))
    }
  }

  # causal_orders is intentionally NULL; see @details
  create_bootstrap_result(adjacency_matrices, total_effects, resampled_indices, causal_orders = NULL)
}
