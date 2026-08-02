# =============================================================================
# Bootstrap for MultiGroup Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================


#' Compute a total-effect matrix from an adjacency matrix via path products
#'
#' Unlike [estimate_all_total_effects()] (regression-based), this sums path
#' products over the DAG defined by `B`, matching the upstream Python
#' MultiGroup bootstrap's `calculate_total_effect()`.
#'
#' @param B adjacency matrix (n_features x n_features), `B[i,j]` = j -> i
#' @param causal_order causal order (1-based indices)
#' @return total-effect matrix (n_features x n_features)
#' @keywords internal
multi_group_total_effect_matrix <- function(B, causal_order) {
  p <- ncol(B)
  TE <- matrix(0, nrow = p, ncol = p)
  n <- length(causal_order)
  if (n < 2) return(TE)
  TE_all <- calculate_total_effects_all(B)
  for (i in seq_len(n - 1)) {
    from_idx <- causal_order[i]
    to_idx <- causal_order[(i + 1):n]
    TE[to_idx, from_idx] <- TE_all[to_idx, from_idx]
  }
  TE
}


#' Bootstrap for Multi-Group Direct LiNGAM
#'
#' @param X_list A list of numeric matrices or data frames (length >= 2), one
#'   per group. Same requirements as [lingam_multi_group()].
#' @param n_sampling Number of bootstrap iterations
#' @param prior_knowledge Prior knowledge matrix (NULL allowed). Applied to
#'   every group, same as [lingam_multi_group()].
#' @param apply_prior_knowledge_softly Apply prior knowledge softly (logical)
#' @param reg_method Regression method ("ols", "lasso", "adaptive_lasso", "ridge")
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")
#' @param init_method Method for estimating the initial weights of adaptive LASSO
#'   regression ("ols" or "ridge")
#' @param seed Random seed (NULL allowed)
#' @param verbose Whether to display progress (logical)
#' @param parallel Whether to use parallel processing (logical)
#' @param n_cores Number of cores to use (integer, NULL allowed). When `NULL`,
#'   the number of cores is limited to a maximum of 2 for safety. Ignored when
#'   `parallel = FALSE`.
#' @param compute_total_effects Whether to also compute total causal effects
#'   for every variable pair on each bootstrap iteration (logical, default
#'   `TRUE`). Set to `FALSE` to skip it when only edge/order stability is
#'   needed.
#' @return A `MultiGroupBootstrapResult`: a named list (one element per group)
#'   of `BootstrapResult` objects (see [lingam_direct_bootstrap()]), with
#'   class `"MultiGroupBootstrapResult"`.
#' @details
#' Each element of the returned list is a regular `BootstrapResult`, so the
#' existing single-group bootstrap functions ([get_probabilities()],
#' [get_causal_direction_counts()], [get_directed_acyclic_graph_counts()],
#' [get_total_causal_effects()], [get_causal_order_stability()],
#' [plot_bootstrap_probabilities()], `tidy()`) all work by extracting a group
#' with `result[[group_name]]` or `result[[i]]`, mirroring the upstream Python
#' API (which returns a list of `BootstrapResult` per group).
#'
#' **Total effects use path products, not regression.** Each bootstrap
#' iteration's total-effect matrix is the sum of path-coefficient products
#' over the DAG defined by that iteration's adjacency matrix (matching the
#' upstream Python MultiGroup bootstrap), which is a different method from
#' [lingam_direct_bootstrap()]'s regression-based
#' [estimate_all_total_effects()].
#'
#' **On iteration failures:** as in [lingam_direct_bootstrap()], each
#' iteration is wrapped in `tryCatch()`; a failing iteration is skipped with a
#' warning, and only if every iteration fails is an error raised.
#'
#' **On reproducibility:** same policy as [lingam_direct_bootstrap()]. During
#' parallel execution, L'Ecuyer parallel random number streams via
#' `parallel::clusterSetRNGStream()` are used. Results are reproducible given
#' the same `seed` and same `n_cores`, but they do not numerically match the
#' results of sequential execution (`parallel = FALSE`). If you need results
#' that exactly match the sequential version, use `parallel = FALSE`.
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#'
#' bs <- lingam_multi_group_bootstrap(mg$data_list,
#'   n_sampling = 10L,
#'   reg_method = "ols",
#'   seed = 42
#' )
#' get_probabilities(bs[[1]])
#'
#' \donttest{
#' bs_par <- lingam_multi_group_bootstrap(mg$data_list,
#'   n_sampling = 30L,
#'   reg_method = "ols",
#'   seed = 42,
#'   parallel = TRUE,
#'   n_cores = 2L
#' )
#' }
lingam_multi_group_bootstrap <- function(X_list,
                                         n_sampling,
                                         prior_knowledge = NULL,
                                         apply_prior_knowledge_softly = FALSE,
                                         reg_method = "adaptive_lasso",
                                         lambda = "BIC",
                                         init_method = "ols",
                                         seed = NULL,
                                         verbose = TRUE,
                                         parallel = FALSE,
                                         n_cores = NULL,
                                         compute_total_effects = TRUE) {
  # --- Validate before starting the cluster (same policy as lingam_direct_bootstrap()) ---
  if (!is.list(X_list) || is.data.frame(X_list)) {
    stop("X_list must be a list of numeric matrices or data frames, one per group.", call. = FALSE)
  }
  if (length(X_list) < 2) {
    stop("X_list must contain at least two items (groups).", call. = FALSE)
  }

  group_names <- names(X_list)
  if (is.null(group_names)) group_names <- character(length(X_list))
  missing_name <- !nzchar(group_names)
  if (any(missing_name)) {
    group_names[missing_name] <- paste0("group", seq_along(X_list))[missing_name]
  }

  n_features_list <- integer(length(X_list))
  X_mats <- vector("list", length(X_list))
  for (d in seq_along(X_list)) {
    Xd <- as.matrix(X_list[[d]])
    if (!is.numeric(Xd)) {
      stop(sprintf("X_list[[%d]] must be a numeric matrix or data frame.", d), call. = FALSE)
    }
    if (anyNA(Xd)) {
      stop(sprintf("X_list[[%d]] must not contain missing values (NA).", d), call. = FALSE)
    }
    if (nrow(Xd) < 3) {
      stop(sprintf("X_list[[%d]] must have at least 3 observations (rows).", d), call. = FALSE)
    }
    X_mats[[d]] <- Xd
    n_features_list[d] <- ncol(Xd)
  }
  if (length(unique(n_features_list)) > 1) {
    stop("All items in X_list must have the same number of columns (variables).", call. = FALSE)
  }
  n_features <- n_features_list[1]
  if (n_features < 2) stop("X_list items must have at least 2 variables (columns).", call. = FALSE)
  names(X_mats) <- group_names

  if (!is.logical(compute_total_effects) || length(compute_total_effects) != 1 ||
        is.na(compute_total_effects)) {
    stop("compute_total_effects must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  reg_args <- validate_reg_args(reg_method, lambda, init_method)
  reg_method <- reg_args$reg_method
  lambda <- reg_args$lambda
  init_method <- reg_args$init_method
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }

  n_groups <- length(X_mats)
  n_list <- vapply(X_mats, nrow, integer(1))

  # Processing for one iteration: resample every group independently ->
  # joint fit -> per-group total effects (path-product based). Wrapped in
  # tryCatch so a single pathological resample does not abort the whole run.
  run_one <- function(i) {
    tryCatch({
      idx_list <- vector("list", n_groups)
      resampled <- vector("list", n_groups)
      for (d in seq_len(n_groups)) {
        idx <- sample(n_list[d], replace = TRUE)
        idx_list[[d]] <- idx
        resampled[[d]] <- X_mats[[d]][idx, , drop = FALSE]
      }
      names(resampled) <- group_names

      result <- lingam_multi_group(
        resampled,
        prior_knowledge = prior_knowledge,
        apply_prior_knowledge_softly = apply_prior_knowledge_softly,
        reg_method = reg_method,
        lambda = lambda,
        init_method = init_method
      )

      total_effects <- if (compute_total_effects) {
        lapply(result$adjacency_matrices, multi_group_total_effect_matrix, causal_order = result$causal_order)
      } else {
        NULL
      }

      list(
        ok = TRUE,
        idx_list = idx_list,
        adjacency_matrices = result$adjacency_matrices,
        total_effects = total_effects,
        causal_order = result$causal_order
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
      "Multi-group bootstrap: %d iterations, %d groups, method=%s (%s)",
      n_sampling, n_groups, reg_method, bootstrap_mode_string(parallel, n_cores)
    ))
    t_start <- proc.time()
  }

  if (parallel) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    setup_cluster_worker(cl, lingam_multi_group)

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

  # Report and drop any failed iterations before aggregating.
  res_list <- filter_bootstrap_failures(res_list)
  n_success <- length(res_list)

  # Aggregate per group into BootstrapResult-shaped arrays. The causal order
  # is shared across groups (jointly estimated), so the same causal_orders
  # matrix is reused for every group's BootstrapResult.
  causal_orders <- matrix(0L, nrow = n_success, ncol = n_features)
  for (i in seq_len(n_success)) causal_orders[i, ] <- res_list[[i]]$causal_order

  res <- lapply(seq_len(n_groups), function(d) {
    adjacency_matrices <- array(0, dim = c(n_success, n_features, n_features))
    total_effects <- if (compute_total_effects) {
      array(0, dim = c(n_success, n_features, n_features))
    } else {
      NULL
    }
    resampled_indices <- vector("list", n_success)
    for (i in seq_len(n_success)) {
      adjacency_matrices[i, , ] <- res_list[[i]]$adjacency_matrices[[d]]
      if (compute_total_effects) total_effects[i, , ] <- res_list[[i]]$total_effects[[d]]
      resampled_indices[[i]] <- res_list[[i]]$idx_list[[d]]
    }
    create_bootstrap_result(adjacency_matrices, total_effects, resampled_indices, causal_orders)
  })
  names(res) <- group_names
  class(res) <- "MultiGroupBootstrapResult"

  if (verbose) bootstrap_completion_message(t_start, n_success, n_sampling)

  res
}


#' Print method for MultiGroupBootstrapResult
#'
#' @param x MultiGroupBootstrapResult object
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' mg <- generate_multi_group_sample()
#' bs <- lingam_multi_group_bootstrap(mg$data_list,
#'   n_sampling = 10L, reg_method = "ols", seed = 42
#' )
#' print(bs)
print.MultiGroupBootstrapResult <- function(x, ...) {
  cat(sprintf("MultiGroupBootstrapResult: %d groups\n", length(x)))
  for (g in names(x)) {
    n_sampling <- dim(x[[g]]$adjacency_matrices)[1]
    n_features <- dim(x[[g]]$adjacency_matrices)[2]
    cat(sprintf("  [%s] %d samplings, %d features\n", g, n_sampling, n_features))
  }
  invisible(x)
}
