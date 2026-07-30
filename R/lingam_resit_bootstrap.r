# =============================================================================
# Bootstrap for RESIT - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (BootstrapMixin in lingam/bootstrap.py,
# applied to lingam/resit.py)
#
# Structure (worker setup, RNG streams, failure handling) follows
# lingam_rcd_bootstrap() (R/lingam_rcd_bootstrap.r).
# =============================================================================


#' Bootstrap for RESIT
#'
#' @param X Numeric matrix (n_samples x n_features)
#' @param n_sampling Number of bootstrap iterations
#' @param regressor Nonlinear regressor, passed to [lingam_resit()]
#' @param alpha Significance level of the HSIC pruning test, passed to
#'   [lingam_resit()]
#' @param prior_knowledge Prior-knowledge matrix, passed to [lingam_resit()]
#' @param seed Random seed (NULL allowed)
#' @param verbose Whether to display progress (logical)
#' @param parallel Whether to use parallel processing (logical)
#' @param n_cores Number of cores to use (integer, NULL allowed)
#' @return A `BootstrapResult` (list); see [lingam_direct_bootstrap()] for the
#'   query helpers that operate on it (`get_probabilities()`,
#'   `get_causal_direction_counts()`, `get_directed_acyclic_graph_counts()`,
#'   `get_causal_order_stability()`).
#' @details
#' **`total_effects` is always `NULL`**: RESIT is a nonlinear method for
#' which total causal effects are undefined, so there is no
#' `compute_total_effects` argument and [get_total_causal_effects()] raises
#' its usual "no total effects" error. (The Python implementation instead
#' stores an all-zero total-effects array; storing nothing is deliberate, so
#' the zeros cannot be mistaken for estimated effects.)
#'
#' **`causal_orders` is populated** (RESIT estimates a full causal order),
#' so [get_causal_order_stability()] works with the returned object.
#'
#' Each iteration runs `O(ncol(X)^2)` nonlinear regressions plus HSIC tests,
#' so a bootstrap is substantially slower than the linear variants; keep
#' `n_sampling` modest. When `parallel = TRUE`, a user-supplied `regressor`
#' function is serialized to the PSOCK workers; any package it uses must be
#' referenced with the `pkg::fun` form inside the function body.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   nonlinear <- generate_resit_sample(n = 300, seed = 1)
#'
#'   bs <- lingam_resit_bootstrap(nonlinear$data,
#'     n_sampling = 3L,
#'     seed = 42
#'   )
#'   get_probabilities(bs)
#' }
#' }
lingam_resit_bootstrap <- function(X,
                                   n_sampling,
                                   regressor = "gam",
                                   alpha = 0.01,
                                   prior_knowledge = NULL,
                                   seed = NULL,
                                   verbose = TRUE,
                                   parallel = FALSE,
                                   n_cores = NULL) {
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
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0) {
    stop("alpha must be a non-negative numeric scalar.", call. = FALSE)
  }
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  # fail fast on an invalid regressor / missing mgcv before spawning workers
  resit_make_regressor(regressor)
  n_samples <- nrow(X)
  n_features <- ncol(X)

  run_one <- function(i) {
    tryCatch({
      idx <- sample(n_samples, replace = TRUE)
      resampled_X <- X[idx, , drop = FALSE]
      result <- lingam_resit(
        resampled_X,
        regressor = regressor,
        alpha = alpha,
        prior_knowledge = prior_knowledge
      )
      list(
        ok               = TRUE,
        idx              = idx,
        adjacency_matrix = result$adjacency_matrix,
        causal_order     = result$causal_order
      )
    }, error = function(e) {
      list(ok = FALSE, iteration = i, message = conditionMessage(e))
    })
  }

  cores <- resolve_bootstrap_cores(parallel, n_cores, n_sampling)
  parallel <- cores$parallel
  n_cores <- cores$n_cores

  if (verbose) {
    message(sprintf(
      "Bootstrap: %d iterations, RESIT (%s)",
      n_sampling, bootstrap_mode_string(parallel, n_cores)
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    setup_cluster_worker(cl, lingam_resit)

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

  res_list <- filter_bootstrap_failures(res_list)
  n_success <- length(res_list)

  adjacency_matrices <- array(0, dim = c(n_success, n_features, n_features))
  causal_orders <- matrix(0L, nrow = n_success, ncol = n_features)
  resampled_indices <- vector("list", n_success)
  for (i in seq_len(n_success)) {
    adjacency_matrices[i, , ] <- res_list[[i]]$adjacency_matrix
    causal_orders[i, ] <- res_list[[i]]$causal_order
    resampled_indices[[i]] <- res_list[[i]]$idx
  }

  if (verbose) bootstrap_completion_message(t_start, n_success, n_sampling)

  # total_effects is intentionally NULL; see @details
  create_bootstrap_result(
    adjacency_matrices,
    total_effects = NULL,
    resampled_indices = resampled_indices,
    causal_orders = causal_orders
  )
}
