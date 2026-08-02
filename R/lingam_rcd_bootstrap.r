# =============================================================================
# Bootstrap for RCD - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://github.com/cdt15/lingam (lingam/rcd.py, 523-561 lines)
#
# Structure (worker setup, RNG streams, failure handling) follows
# lingam_parce_bootstrap() (R/lingam_parce_bootstrap.r).
# =============================================================================


#' Bootstrap for RCD
#'
#' @param X Numeric matrix (n_samples x n_features)
#' @param n_sampling Number of bootstrap iterations
#' @param max_explanatory_num Maximum number of explanatory variables, passed to [lingam_rcd()]
#' @param cor_alpha Significance level for correlation tests, passed to [lingam_rcd()]
#' @param ind_alpha Significance level for the HSIC independence test, passed to [lingam_rcd()]
#' @param shapiro_alpha Significance level for the non-Gaussianity test, passed to [lingam_rcd()]
#' @param MLHSICR Whether to use MLHSICR regression, passed to [lingam_rcd()]
#' @param independence Independence measure, passed to [lingam_rcd()]
#' @param ind_corr F-correlation rejection threshold, passed to [lingam_rcd()]
#' @param seed Random seed (NULL allowed)
#' @param verbose Whether to display progress (logical)
#' @param parallel Whether to use parallel processing (logical)
#' @param n_cores Number of cores to use (integer, NULL allowed)
#' @param compute_total_effects Whether to also estimate total causal effects
#'   for every ancestor pair on each bootstrap iteration (logical, default `TRUE`).
#' @return A `BootstrapResult` (list); see [lingam_direct_bootstrap()] for the
#'   query helpers that operate on it (`get_probabilities()`,
#'   `get_causal_direction_counts()`, `get_directed_acyclic_graph_counts()`,
#'   `get_total_causal_effects()`).
#' @details
#' **Total effects are computed only for ancestor pairs**, driven by each
#' iteration's `ancestors_list` (unlike [lingam_direct_bootstrap()] and
#' [lingam_parce_bootstrap()], which loop over all variable pairs). For a
#' pair `(from, to)` with `from` in `to`'s ancestor set, the total effect is
#' obtained via [calculate_total_effect()] (summing products of
#' adjacency-matrix coefficients along every directed path). If `from`'s row
#' in the adjacency matrix contains `NA` (it shares a latent confounder with
#' some other variable), the effect is set to `NA` for that iteration.
#'
#' **`NA` (confounded) edges are treated as absent when aggregating.** Both
#' the adjacency matrix and the total-effect matrix have `NA` replaced by
#' `0` before being stored in the returned `BootstrapResult`, matching the
#' policy used by [lingam_parce_bootstrap()].
#'
#' **`causal_orders` is not populated**: RCD has no causal order (see
#' [lingam_rcd()]). As a result, [get_causal_order_stability()] cannot be
#' used with a `BootstrapResult` returned by this function.
#'
#' RCD's `fit` step is HSIC-heavy and can be slow per iteration, especially
#' with `MLHSICR = TRUE`; keep `n_sampling` modest in examples.
#' @export
#' @examples
#' \donttest{
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#'
#' bs <- lingam_rcd_bootstrap(confounded$data,
#'   n_sampling = 5L,
#'   seed = 42
#' )
#' get_probabilities(bs)
#' }
lingam_rcd_bootstrap <- function(X,
                                 n_sampling,
                                 max_explanatory_num = 2L,
                                 cor_alpha = 0.01,
                                 ind_alpha = 0.01,
                                 shapiro_alpha = 0.01,
                                 MLHSICR = FALSE,
                                 independence = "hsic",
                                 ind_corr = 0.5,
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
      result <- lingam_rcd(
        resampled_X,
        max_explanatory_num = max_explanatory_num,
        cor_alpha = cor_alpha,
        ind_alpha = ind_alpha,
        shapiro_alpha = shapiro_alpha,
        MLHSICR = MLHSICR,
        independence = independence,
        ind_corr = ind_corr
      )
      B <- result$adjacency_matrix
      M <- result$ancestors_list

      te <- if (compute_total_effects) {
        has_na_row <- apply(B, 1, anyNA)
        TE <- matrix(0, n_features, n_features)
        TE_all <- calculate_total_effects_all(B)
        for (to_idx in seq_len(n_features)) {
          from_idx <- M[[to_idx]]
          if (length(from_idx) > 0) {
            TE[to_idx, from_idx] <- ifelse(
              has_na_row[from_idx], NA_real_, TE_all[to_idx, from_idx]
            )
          }
        }
        TE
      } else {
        NULL
      }

      # NA (confounded) edges are treated as absent when aggregating (see @details)
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

  cores <- resolve_bootstrap_cores(parallel, n_cores, n_sampling)
  parallel <- cores$parallel
  n_cores <- cores$n_cores

  if (verbose) {
    message(sprintf(
      "Bootstrap: %d iterations, RCD (%s)",
      n_sampling, bootstrap_mode_string(parallel, n_cores)
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    setup_cluster_worker(cl, lingam_rcd)

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

  if (verbose) bootstrap_completion_message(t_start, n_success, n_sampling)

  # causal_orders is intentionally NULL; see @details
  create_bootstrap_result(adjacency_matrices, total_effects, resampled_indices, causal_orders = NULL)
}
