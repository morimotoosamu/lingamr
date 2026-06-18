# =============================================================================
# Bootstrap for Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================

# =============================================================================
# Bootstrap execution function
# =============================================================================

#' Bootstrap for Direct LiNGAM
#'
#' @param X Numeric matrix (n_samples x n_features)
#' @param n_sampling Number of bootstrap iterations
#' @param prior_knowledge Prior knowledge matrix (NULL allowed)
#' @param apply_prior_knowledge_softly Apply prior knowledge softly (logical)
#' @param measure Independence measure ("pwling" or "kernel")
#' @param reg_method Regression method ("ols", "lasso", "adaptive_lasso", "ridge")
#' @param lambda Lambda selection ("lambda.min", "lambda.1se", "AIC", "BIC","oracle")
#' @param init_method Method for estimating the initial weights of adaptive LASSO
#'   regression ("ols" or "ridge"). Same as the argument of the same name in
#'   `lingam_direct()`.
#' @param seed Random seed (NULL allowed)
#' @param verbose Whether to display progress (logical)
#' @param parallel Whether to use parallel processing (logical). When `TRUE`,
#'   each bootstrap iteration is distributed across multiple cores.
#' @param n_cores Number of cores to use (integer, NULL allowed). When `NULL`,
#'   the number of cores is limited to a maximum of 2 for safety. Ignored when
#'   `parallel = FALSE`.
#' @return BootstrapResult (list)
#' @details
#' When `parallel = TRUE` is specified, iterations are distributed across a
#' socket cluster created by `parallel::makePSOCKcluster()`. The cluster is
#' always released via `on.exit()`, whether the process finishes normally or
#' an error occurs.
#'
#' **On reproducibility:** During parallel execution, L'Ecuyer parallel random
#' number streams via `parallel::clusterSetRNGStream()` are used. Results are
#' reproducible given the same `seed` and same `n_cores`, but they do not
#' numerically match the results of sequential execution (`parallel = FALSE`).
#' If you need results that exactly match the sequential version, use
#' `parallel = FALSE`.
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # Fast example with OLS
#' bs <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
#'   n_sampling = 10L,
#'   reg_method = "ols",
#'   seed = 42
#' )
#' get_probabilities(bs)
#'
#' \donttest{
#' # With LASSO (requires glmnet)
#' bs_lasso <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
#'   n_sampling = 30L,
#'   seed = 42
#' )
#'
#' # Parallel execution on 2 cores
#' bs_par <- lingam_direct_bootstrap(LiNGAM_sample_1000$data,
#'   n_sampling = 30L,
#'   seed = 42,
#'   parallel = TRUE,
#'   n_cores = 2L
#' )
#' }
lingam_direct_bootstrap <- function(X,
                             n_sampling,
                             prior_knowledge = NULL,
                             apply_prior_knowledge_softly = FALSE,
                             measure = "pwling",
                             reg_method = "adaptive_lasso",
                             lambda = "BIC",
                             init_method = "ols",
                             seed = NULL,
                             verbose = TRUE,
                             parallel = FALSE,
                             n_cores = NULL) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(X)) stop("X must not contain missing values (NA).", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  # Invalid arguments would otherwise produce confusing errors inside the
  # iterations (within workers when parallel), so validate them here before
  # starting the cluster.
  measure <- match.arg(measure, c("pwling", "kernel"))
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols", "ridge"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  if (reg_method == "ridge" && lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }
  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Processing for one iteration: resampling -> estimation -> total effects
  run_one <- function(i) {
    idx <- sample(n_samples, replace = TRUE)
    resampled_X <- X[idx, , drop = FALSE]
    result <- lingam_direct(
      resampled_X,
      prior_knowledge = prior_knowledge,
      apply_prior_knowledge_softly = apply_prior_knowledge_softly,
      measure = measure,
      reg_method = reg_method,
      lambda = lambda,
      init_method = init_method
    )
    te <- estimate_all_total_effects(
      resampled_X, result,
      method = reg_method,
      lambda = lambda,
      init_method = init_method
    )
    list(
      idx              = idx,
      adjacency_matrix = result$adjacency_matrix,
      total_effects    = te,
      causal_order     = result$causal_order
    )
  }

  # Resolve the number of cores to use (defaults to a maximum of 2 cores)
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

    # Make this package's functions available on the workers.
    # requireNamespace only loads the namespace without adding it to the search
    # path, so attach the package with library() so that functions referenced
    # inside the closure can be resolved on the worker side.
    pkg <- utils::packageName()
    ns <- environment(lingam_direct)
    # Because the workers are fresh R sessions, set the same library paths as
    # the main session so that the same version of this package can be loaded.
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    worker_ready <- FALSE
    if (!is.null(pkg) && nzchar(pkg)) {
      worker_ready <- all(unlist(parallel::clusterCall(cl, function(p) {
        suppressWarnings(suppressMessages(
          isTRUE(requireNamespace(p, quietly = TRUE))
        ))
      }, pkg)))
      if (worker_ready) {
        parallel::clusterCall(cl, function(p) {
          suppressWarnings(suppressMessages(library(p, character.only = TRUE)))
          invisible(NULL)
        }, pkg)
      }
    }
    if (!worker_ready) {
      # Fallback when not installed (e.g. devtools::load_all):
      # send all objects of the namespace to the workers.
      parallel::clusterExport(cl, varlist = ls(ns, all.names = TRUE), envir = ns)
    }

    # Parallel-safe random number stream (ensures reproducibility; does not
    # match the sequential version)
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

  # Aggregate results into arrays
  adjacency_matrices <- array(0, dim = c(n_sampling, n_features, n_features))
  total_effects <- array(0, dim = c(n_sampling, n_features, n_features))
  causal_orders <- matrix(0L, nrow = n_sampling, ncol = n_features)
  resampled_indices <- vector("list", n_sampling)
  for (i in seq_len(n_sampling)) {
    adjacency_matrices[i, , ] <- res_list[[i]]$adjacency_matrix
    total_effects[i, , ] <- res_list[[i]]$total_effects
    causal_orders[i, ] <- res_list[[i]]$causal_order
    resampled_indices[[i]] <- res_list[[i]]$idx
  }

  if (verbose) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    message(sprintf("Completed in %.1f seconds.", elapsed))
  }
  create_bootstrap_result(adjacency_matrices, total_effects, resampled_indices, causal_orders)
}


# =============================================================================
# BootstrapResult object
# =============================================================================

#' Create a BootstrapResult
#' @param adjacency_matrices array (n_sampling x n_features x n_features)
#' @param total_effects array (n_sampling x n_features x n_features)
#' @param resampled_indices list of index vectors
#' @param causal_orders matrix (n_sampling x n_features). Each row is the causal
#'   order of one sample.
#' @return BootstrapResult (list with class attribute)
#' @keywords internal
create_bootstrap_result <- function(adjacency_matrices, total_effects, resampled_indices = NULL, causal_orders = NULL) {
  obj <- list(
    adjacency_matrices = adjacency_matrices,
    total_effects      = total_effects,
    resampled_indices  = resampled_indices,
    causal_orders      = causal_orders
  )
  class(obj) <- "BootstrapResult"
  return(obj)
}


#' Display the contents of a BootstrapResult
#'
#' @param x BootstrapResult object
#' @param ... Additional arguments (for S3 method compatibility)
#' @method print BootstrapResult
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#'
#' print(bs_model)
print.BootstrapResult <- function(x, ...) {
  n_sampling <- dim(x$adjacency_matrices)[1]
  n_features <- dim(x$adjacency_matrices)[2]
  cat(sprintf("BootstrapResult: %d samplings, %d features\n", n_sampling, n_features))
}


# =============================================================================
# BootstrapResult methods
# =============================================================================

#' Get counts, proportions, and causal effects of causal directions
#'
#' @param result BootstrapResult object
#' @param n_directions How many of the top entries to return (NULL = all)
#' @param min_causal_effect Minimum threshold for the causal effect (NULL = 0)
#' @param split_by_causal_effect_sign Whether to split by the sign of the causal effect
#' @param labels Vector of variable names (NULL allowed; if provided, adds from_name and to_name columns)
#' @return A data frame containing the following columns:
#' * `from`, `to`: 1-based indices of the causal (from) and effect (to) variables.
#' * `count`: Number of bootstrap samples in which this specific causal direction was identified.
#' * `proportion`: The frequency of the direction's occurrence (count / n_sampling), representing its bootstrap probability.
#' * `mean_effect`: The average value of the estimated causal effects across samples where this direction was identified.
#' * `median_effect`: The median value of the estimated causal effects, providing a robust estimate of the effect size.
#' * `sd_effect`: The standard deviation of the causal effect estimates, indicating the stability of the effect size.
#' * `ci_lower`, `ci_upper`: The lower (2.5%) and upper (97.5%) bounds of the bootstrap confidence interval for the causal effect.
#' * `sign` (optional): The sign of the causal effect (1 for positive, -1 for negative), included if `split_by_causal_effect_sign = TRUE`.
#' * `from_name`, `to_name` (optional): Character labels for the variables, included if `labels` were provided.
#' @importFrom stats sd median quantile
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#'
#' get_causal_direction_counts(bs_model, labels = names(LiNGAM_sample_1000$data))
get_causal_direction_counts <- function(result,
                                        n_directions = NULL,
                                        min_causal_effect = NULL,
                                        split_by_causal_effect_sign = FALSE,
                                        labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]
  n_features <- dim(am_array)[2]

  # Collect directions and effect sizes from all bootstrap samples
  directions_list <- list()
  for (s in seq_len(n_sampling)) {
    am <- am_array[s, , ]
    idx <- which(abs(am) > min_causal_effect, arr.ind = TRUE)
    if (nrow(idx) == 0) next

    effects <- sapply(seq_len(nrow(idx)), function(k) am[idx[k, 1], idx[k, 2]])

    if (split_by_causal_effect_sign) {
      directions_list[[length(directions_list) + 1]] <- data.frame(
        from   = idx[, 2],
        to     = idx[, 1],
        sign   = as.integer(sign(effects)),
        effect = effects
      )
    } else {
      directions_list[[length(directions_list) + 1]] <- data.frame(
        from   = idx[, 2],
        to     = idx[, 1],
        effect = effects
      )
    }
  }

  # When empty
  if (length(directions_list) == 0) {
    empty <- data.frame(
      from = integer(0), to = integer(0),
      count = integer(0), proportion = numeric(0),
      mean_effect = numeric(0), median_effect = numeric(0),
      sd_effect = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)
    )
    if (split_by_causal_effect_sign) empty$sign <- integer(0)
    if (!is.null(labels)) {
      empty$from_name <- character(0)
      empty$to_name <- character(0)
    }
    return(empty)
  }

  directions <- do.call(rbind, directions_list)

  # Build group keys
  if (split_by_causal_effect_sign) {
    group_key <- paste(directions$from, directions$to, directions$sign, sep = "_")
  } else {
    group_key <- paste(directions$from, directions$to, sep = "_")
  }

  # Aggregate by group
  unique_keys <- unique(group_key)
  results_list <- lapply(unique_keys, function(key) {
    mask <- group_key == key
    subset_df <- directions[mask, ]
    effects <- subset_df$effect

    row <- data.frame(
      from = subset_df$from[1],
      to = subset_df$to[1],
      count = length(effects),
      proportion = length(effects) / n_sampling,
      mean_effect = base::mean(effects),
      median_effect = stats::median(effects),
      sd_effect = if (length(effects) > 1) stats::sd(effects) else 0,
      ci_lower = stats::quantile(effects, 0.025, names = FALSE),
      ci_upper = stats::quantile(effects, 0.975, names = FALSE)
    )

    if (split_by_causal_effect_sign) {
      row$sign <- subset_df$sign[1]
    }

    row
  })

  agg <- do.call(rbind, results_list)

  # Sort in descending order
  agg <- agg[order(-agg$count), ]
  rownames(agg) <- NULL

  # Add variable names
  if (!is.null(labels)) {
    agg$from_name <- labels[agg$from]
    agg$to_name <- labels[agg$to]
  }

  # Top n_directions entries
  if (!is.null(n_directions)) {
    n_directions <- min(n_directions, nrow(agg))
    agg <- agg[seq_len(n_directions), ]
  }

  return(agg)
}


#' Get DAG counts
#'
#' @param result BootstrapResult object
#' @param n_dags How many of the top entries to return (NULL = all)
#' @param min_causal_effect Minimum threshold for the causal effect (NULL = 0)
#' @param split_by_causal_effect_sign Whether to split by the sign of the causal effect
#' @return list(dag = list of data.frames, count = integer vector)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#'
#' get_directed_acyclic_graph_counts(bs_model)
get_directed_acyclic_graph_counts <- function(result,
                                              n_dags = NULL,
                                              min_causal_effect = NULL,
                                              split_by_causal_effect_sign = FALSE) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]

  # Convert each DAG to a string key
  dag_keys <- character(n_sampling)
  dag_list <- vector("list", n_sampling)

  if (split_by_causal_effect_sign) {
    sign_array <- sign(am_array)
    sign_array[abs(am_array) <= min_causal_effect] <- 0
    dag_keys <- apply(sign_array, MARGIN = 1, FUN = paste, collapse = ",")
    dag_list <- lapply(seq_len(n_sampling), function(s) sign_array[s, , ])
  } else {
    bin_array <- abs(am_array) > min_causal_effect
    mode(bin_array) <- "integer" # convert to integer for string conversion
    dag_keys <- apply(bin_array, MARGIN = 1, FUN = paste, collapse = ",")
    dag_list <- lapply(seq_len(n_sampling), function(s) bin_array[s, , ])
  }

  # Count
  tbl <- sort(table(dag_keys), decreasing = TRUE)
  if (!is.null(n_dags)) {
    n_dags <- min(n_dags, length(tbl))
    tbl <- tbl[seq_len(n_dags)]
  }

  # Build the result
  dags_result <- list()
  counts_result <- as.integer(tbl)

  for (i in seq_along(tbl)) {
    key <- names(tbl)[i]
    match_idx <- which(dag_keys == key)[1]
    mat <- dag_list[[match_idx]]

    # Fix: unify both to use which(mat != 0)
    idx <- which(mat != 0, arr.ind = TRUE)

    if (nrow(idx) == 0) {
      if (split_by_causal_effect_sign) {
        dags_result[[i]] <- data.frame(
          from = integer(0), to = integer(0),
          sign = integer(0)
        )
      } else {
        dags_result[[i]] <- data.frame(from = integer(0), to = integer(0))
      }
    } else if (split_by_causal_effect_sign) {
      dags_result[[i]] <- data.frame(
        from = idx[, 2],
        to   = idx[, 1],
        sign = sapply(seq_len(nrow(idx)), function(k) mat[idx[k, 1], idx[k, 2]])
      )
    } else {
      dags_result[[i]] <- data.frame(
        from = idx[, 2],
        to   = idx[, 1]
      )
    }
  }

  return(list(dag = dags_result, count = counts_result))
}


#' Get bootstrap probabilities
#'
#' @param result BootstrapResult object
#' @param min_causal_effect Minimum threshold for the causal effect (NULL = 0)
#' @return Probability matrix (n_features x n_features)
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#'
#' get_probabilities(bs_model)
get_probabilities <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0

  bp <- apply(abs(am_array) > min_causal_effect, MARGIN = c(2, 3), FUN = mean)

  return(bp)
}


#' Get a list of total causal effects
#'
#' @param result BootstrapResult object
#' @param min_causal_effect Minimum threshold for the causal effect (NULL = 0)
#' @return data.frame (from, to, effect, probability)
#' @importFrom stats median
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#'
#' get_total_causal_effects(bs_model)
get_total_causal_effects <- function(result, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  te_array <- result$total_effects
  te_array[is.na(te_array)] <- 0

  probs <- apply(abs(te_array) > min_causal_effect, MARGIN = c(2, 3), FUN = mean)

  # Causal directions with probability > 0
  idx <- which(abs(probs) > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      effect = numeric(0), probability = numeric(0)
    ))
  }

  from_vec <- idx[, 2]
  to_vec <- idx[, 1]

  prob_vec <- probs[idx]

  median_mat <- apply(te_array, MARGIN = c(2, 3), FUN = function(x) {
    nonzero <- x[abs(x) > 0]
    if (length(nonzero) == 0) {
      return(0)
    }
    stats::median(nonzero)
  })

  # Extract the medians in one step using idx as well
  effect_vec <- median_mat[idx]

  # Sort by probability in descending order
  ord <- order(-prob_vec)
  data.frame(
    from        = from_vec[ord],
    to          = to_vec[ord],
    effect      = effect_vec[ord],
    probability = prob_vec[ord]
  )
}


#' Get all paths between two specified variables and their bootstrap probabilities
#'
#' @param result BootstrapResult object
#' @param from_index Start index (1-based)
#' @param to_index End index (1-based)
#' @param min_causal_effect Minimum threshold for the causal effect (NULL = 0)
#' @return data.frame (path, effect, probability)
#' @importFrom stats median
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#' get_paths(bs_model, 1, 6)
get_paths <- function(result, from_index, to_index, min_causal_effect = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")

  am_array <- result$adjacency_matrices
  n_sampling <- dim(am_array)[1]

  # Collect all paths
  # Accumulate into a list first, then unlist at the end
  paths_collector <- vector("list", n_sampling)
  effects_collector <- vector("list", n_sampling)

  for (s in seq_len(n_sampling)) {
    am <- am_array[s, , ]
    res <- find_all_paths(am, from_index, to_index, min_causal_effect)
    if (length(res$paths) > 0) {
      paths_collector[[s]] <- vapply(res$paths, paste, "", collapse = "_")
      effects_collector[[s]] <- res$effects
    }
  }
  paths_str_list <- unlist(paths_collector)
  effects_list <- unlist(effects_collector)

  if (length(paths_str_list) == 0) {
    return(data.frame(path = character(0), effect = numeric(0), probability = numeric(0)))
  }

  # Count
  tbl <- table(paths_str_list)
  tbl <- sort(tbl, decreasing = TRUE)
  probs <- as.numeric(tbl) / n_sampling

  # Median effect of each path
  path_strs <- names(tbl)
  effects_median <- sapply(path_strs, function(ps) {
    stats::median(effects_list[paths_str_list == ps])
  })

  # Convert path strings to lists
  path_list <- lapply(path_strs, function(ps) {
    as.integer(strsplit(ps, "_")[[1]])
  })

  data.frame(
    path        = I(path_list),
    effect      = as.numeric(effects_median),
    probability = probs,
    row.names   = NULL
  )
}


# =============================================================================
# Helpers for displaying and visualizing results
# =============================================================================

#' Draw bootstrap probabilities with DiagrammeR
#'
#' @param result BootstrapResult object
#' @param labels Vector of variable names (NULL allowed)
#' @param min_causal_effect Minimum causal effect to display
#' @param min_probability Minimum probability to display
#' @param rankdir Layout direction
#' @param shape Node shape
#' @return grViz object
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#' plot_bootstrap_probabilities(bs_model)
plot_bootstrap_probabilities <- function(result,
                                         labels = NULL,
                                         min_causal_effect = NULL,
                                         min_probability = 0.5,
                                         rankdir = "TB",
                                         shape = "circle") {
  stopifnot(inherits(result, "BootstrapResult"))

  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required. Please install it.", call. = FALSE)
  }

  bp <- get_probabilities(result, min_causal_effect)
  p <- ncol(bp)
  if (is.null(labels)) labels <- paste0("x", seq_len(p) - 1)

  # Generate edges
  edge_lines <- c()
  for (i in 1:p) {
    for (j in 1:p) {
      if (bp[i, j] >= min_probability) {
        edge_lines <- c(
          edge_lines,
          sprintf(
            "  %s -> %s [label = ' %.2f', penwidth = %.1f]",
            labels[j], labels[i], bp[i, j], bp[i, j] * 3 + 0.5
          )
        )
      }
    }
  }

  if (length(edge_lines) == 0) {
    message("No edges exceed the specified threshold. Please lower 'min_probability'.")
    return(invisible(NULL))
  }

  dot <- paste0(
    "digraph bootstrap_result {\n",
    sprintf("  graph [rankdir = %s, fontsize = 14,\n", rankdir),
    "         label = 'Bootstrap Probabilities',\n",
    "         labelloc = t, fontname = 'Helvetica-Bold']\n",
    sprintf("  node [shape = %s, style = filled, fillcolor = lightyellow,\n", shape),
    "        fontname = Helvetica, fontsize = 14, width = 0.6]\n",
    "  edge [fontname = Helvetica, fontsize = 10, fontcolor = blue, color = gray40]\n\n",
    paste(edge_lines, collapse = "\n"), "\n",
    "}\n"
  )

  DiagrammeR::grViz(dot)
}


#' Create an adjacency matrix of representative causal-effect values from bootstrap results
#'
#' @param result BootstrapResult object
#' @param stat Representative statistic ("mean" or "median")
#' @param min_causal_effect Minimum threshold for the causal effect (values at or
#'   below this are treated as zero) (NULL = 0)
#' @param min_probability Edges below this probability are set to zero (NULL = 0)
#' @param labels Vector of variable names (NULL allowed)
#' @return Adjacency matrix (n_features x n_features).
#'   **Rule: `B[i, j]` is the causal coefficient from variable j to variable i (j -> i).**
#'   Same rule as the `adjacency_matrix` of `lingam_direct()`.
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' bs_model <- lingam_direct_bootstrap(LiNGAM_sample_1000$data, n_sampling = 30L, seed = 42)
#' get_adjacency_matrix_summary(bs_model)
get_adjacency_matrix_summary <- function(result,
                                         stat = "median",
                                         min_causal_effect = NULL,
                                         min_probability = NULL,
                                         labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))
  if (!(stat %in% c("mean", "median"))) {
    stop("'stat' must be either 'mean' or 'median'.")
  }

  if (is.null(min_causal_effect)) min_causal_effect <- 0.0
  if (min_causal_effect < 0) stop("min_causal_effect must be >= 0.")
  if (is.null(min_probability)) min_probability <- 0.0

  am_array <- result$adjacency_matrices
  am_array[is.na(am_array)] <- 0
  n_sampling <- dim(am_array)[1]
  n_features <- dim(am_array)[2]

  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (i in 1:n_features) {
    for (j in 1:n_features) {
      if (i == j) next

      # Get the (i, j) element across all bootstrap samples
      vals <- am_array[, i, j]

      # Extract only values exceeding the threshold
      significant <- vals[abs(vals) > min_causal_effect]

      # Compute the probability
      prob <- length(significant) / n_sampling

      if (prob < min_probability || length(significant) == 0) {
        B[i, j] <- 0
      } else {
        B[i, j] <- if (stat == "mean") mean(significant) else median(significant)
      }
    }
  }

  # Add variable names
  if (!is.null(labels)) {
    rownames(B) <- labels
    colnames(B) <- labels
  }

  return(B)
}
