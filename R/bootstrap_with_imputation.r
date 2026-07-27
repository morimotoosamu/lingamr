# =============================================================================
# Bootstrap with Multiple Imputation - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam (lingam.tools.bootstrap_with_imputation)
#
# License: MIT + file LICENSE
#
# Causal discovery on data with missing values: each bootstrap resample is
# multiply imputed, and a common causal structure is jointly estimated across
# the imputed datasets via lingam_multi_group() (Shimizu 2012), treating the
# n_repeats imputed copies as "groups" that share one causal order.
# =============================================================================


#' Signal a hook contract violation (imputer / cd_fit return value)
#'
#' A plain error, but always called from outside the per-iteration
#' `tryCatch()` in [bootstrap_with_imputation()] (see `validate_imputer_output()`
#' / `validate_cd_fit_output()` call sites there), so that a programming error
#' in a user-supplied hook aborts the whole call immediately instead of being
#' swallowed as a per-iteration stochastic estimation failure.
#' @param msg Error message
#' @keywords internal
contract_violation <- function(msg) {
  stop(msg, call. = FALSE)
}


#' Check that the 'mice' package is available
#' @param context Short description used in the error message
#' @keywords internal
check_mice_available <- function(context) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop(sprintf(
      "Package 'mice' is required for %s. Please install it, or supply a custom `imputer` function.",
      context
    ), call. = FALSE)
  }
}


#' Default imputer: multiple imputation via mice::mice(method = "norm")
#'
#' The closest standard R equivalent of the upstream Python default
#' (`sklearn.impute.IterativeImputer(sample_posterior = TRUE)`, Bayesian
#' linear regression chained equations with posterior sampling) is
#' `mice::mice(method = "norm")` (Bayesian linear regression). Numeric values
#' will not match the Python implementation; the multiple-imputation design
#' is equivalent.
#'
#' @param X_boot Bootstrap-resampled data (matrix, may contain NA)
#' @param n_repeats Number of imputed datasets to generate
#' @return A list of `n_repeats` complete numeric matrices
#' @keywords internal
default_imputer <- function(X_boot, n_repeats) {
  check_mice_available("the default imputer")
  # mice emits informational warnings (e.g. near-collinearity among
  # predictors) via loggedEvents that are not actionable here; suppress them
  # so they do not surface as spurious R warnings on every bootstrap iteration.
  imp <- suppressWarnings(
    mice::mice(as.data.frame(X_boot), m = n_repeats, method = "norm", printFlag = FALSE)
  )
  lapply(seq_len(n_repeats), function(k) as.matrix(mice::complete(imp, k)))
}


#' Default causal-discovery fit: joint estimation via lingam_multi_group()
#'
#' Treats the imputed datasets as "groups" sharing a common causal order.
#'
#' @param X_list List of imputed datasets (one per repeat)
#' @param prior_knowledge Prior knowledge matrix (NULL allowed)
#' @param apply_prior_knowledge_softly Apply prior knowledge softly (logical)
#' @return `list(causal_order = <integer vector>, adjacency_matrices = <list>)`
#' @keywords internal
default_cd_fit <- function(X_list, prior_knowledge, apply_prior_knowledge_softly) {
  res <- lingam_multi_group(X_list,
    prior_knowledge = prior_knowledge,
    apply_prior_knowledge_softly = apply_prior_knowledge_softly
  )
  list(causal_order = res$causal_order, adjacency_matrices = res$adjacency_matrices)
}


#' Validate the return value of an `imputer` function
#' @param datasets Return value of `imputer(X_boot)`
#' @param X_boot The bootstrap-resampled data passed to `imputer`
#' @keywords internal
validate_imputer_output <- function(datasets, X_boot) {
  n <- nrow(X_boot)
  p <- ncol(X_boot)
  if (!is.list(datasets) || length(datasets) < 1) {
    contract_violation(paste0(
      "The return value of imputer violates its specification: ",
      "it must be a non-empty list of imputed matrices."
    ))
  }
  for (k in seq_along(datasets)) {
    Xk <- datasets[[k]]
    if (!is.numeric(Xk) || !is.matrix(Xk)) {
      Xk <- tryCatch(as.matrix(Xk), error = function(e) NULL)
    }
    if (is.null(Xk) || !is.numeric(Xk) || !all(dim(Xk) == c(n, p))) {
      contract_violation(sprintf(
        "The return value of imputer violates its specification: element %d must be a numeric matrix of dimension %d x %d.",
        k, n, p
      ))
    }
    if (anyNA(Xk)) {
      contract_violation(sprintf(
        "The return value of imputer violates its specification: element %d still contains NA (imputation is incomplete).",
        k
      ))
    }
    datasets[[k]] <- Xk
  }
  datasets
}


#' Validate the return value of a `cd_fit` function
#' @param cd_res Return value of `cd_fit(X_list)`
#' @param p Number of features
#' @param n_datasets Expected number of adjacency matrices (= number of imputed datasets)
#' @keywords internal
validate_cd_fit_output <- function(cd_res, p, n_datasets) {
  if (!is.list(cd_res) || !all(c("causal_order", "adjacency_matrices") %in% names(cd_res))) {
    contract_violation(paste0(
      "The return value of cd_fit violates its specification: ",
      "it must be a list with elements `causal_order` and `adjacency_matrices`."
    ))
  }
  if (!is.null(dim(cd_res$causal_order))) {
    contract_violation(paste0(
      "The return value of cd_fit violates its specification: ",
      "causal_order must be a plain vector, not a matrix/array."
    ))
  }
  co <- suppressWarnings(as.integer(cd_res$causal_order))
  if (length(co) != p || anyNA(co) || !identical(sort(co), seq_len(p))) {
    contract_violation(sprintf(
      "The return value of cd_fit violates its specification: causal_order must be a permutation of 1:%d.",
      p
    ))
  }
  ams <- cd_res$adjacency_matrices
  if (!is.list(ams) || length(ams) != n_datasets) {
    contract_violation(sprintf(
      "The return value of cd_fit violates its specification: adjacency_matrices must be a list of length %d (one per imputed dataset).",
      n_datasets
    ))
  }
  for (k in seq_along(ams)) {
    if (!is.numeric(ams[[k]]) || !all(dim(as.matrix(ams[[k]])) == c(p, p))) {
      contract_violation(sprintf(
        "The return value of cd_fit violates its specification: adjacency_matrices[[%d]] must be a %d x %d numeric matrix.",
        k, p, p
      ))
    }
  }
  invisible(NULL)
}


#' Bootstrap with Multiple Imputation for Direct LiNGAM
#'
#' Causal discovery on data containing missing values (NA). Each bootstrap
#' resample (drawn with replacement, missing values retained) is multiply
#' imputed into `n_repeats` complete datasets, and a common causal structure
#' is jointly estimated across those datasets with [lingam_multi_group()]
#' (Shimizu 2012), treating the imputed copies as "groups" sharing one causal
#' order. R port of the Python `lingam.tools.bootstrap_with_imputation()`.
#'
#' @param X A numeric matrix or data frame (n_samples x n_features). May
#'   contain `NA`. If `X` has no missing values, a warning suggests using
#'   [lingam_direct_bootstrap()] instead, and estimation proceeds anyway.
#' @param n_sampling Number of bootstrap iterations (positive integer)
#' @param n_repeats Number of imputed datasets generated per bootstrap sample
#'   (positive integer, default `10L`). Ignored when a custom `imputer` is
#'   supplied; the number of datasets it returns is used instead.
#' @param imputer `NULL`, or a `function(X_boot)` returning a list of
#'   complete (no-NA) numeric matrices, each with the same dimensions as
#'   `X_boot`. Defaults to multiple imputation via `mice::mice(method = "norm")`.
#' @param cd_fit `NULL`, or a `function(X_list)` returning
#'   `list(causal_order = <integer vector, 1-based permutation>,
#'   adjacency_matrices = <list of p x p matrices, one per element of X_list>)`.
#'   Defaults to joint estimation via [lingam_multi_group()].
#' @param prior_knowledge Prior knowledge matrix (NULL allowed). Only used
#'   when `cd_fit = NULL`; a warning is issued if supplied together with a
#'   custom `cd_fit`.
#' @param apply_prior_knowledge_softly Apply prior knowledge softly (logical).
#'   Same restriction as `prior_knowledge`.
#' @param seed Random seed (NULL allowed). Set once before the bootstrap loop;
#'   governs both the resampling and (via the global RNG) `mice`'s imputation.
#' @param verbose Whether to display progress (logical)
#' @return An `ImputationBootstrapResult` (list) containing:
#' * `causal_orders`: `n_sampling` x `p` integer matrix (1-based causal order per iteration).
#' * `adjacency_matrices`: `array(n_sampling, n_repeats, p, p)`; `[, , i, j]`
#'   follows the lingamr convention (`B[i, j]` = coefficient of j -> i).
#' * `resampled_indices`: `n_sampling` x `n` integer matrix of resampled row indices.
#' * `imputation_results`: `array(n_sampling, n_repeats, n, p)`; non-`NA` only
#'   at positions that were missing in that iteration's bootstrap resample.
#' @details
#' **Procedure:** for each of `n_sampling` iterations, (1) resample `X` with
#' replacement (missing values are retained), (2) impute the resample into
#' `n_repeats` complete datasets, (3) jointly estimate one causal structure
#' shared by all `n_repeats` datasets with [lingam_multi_group()]. This
#' assumes the same causal structure underlies every imputed copy.
#'
#' **Default imputer.** `mice::mice(method = "norm")` (Bayesian linear
#' regression, multiple imputation by chained equations) is the closest
#' standard R analogue of the upstream Python default
#' (`IterativeImputer(sample_posterior = TRUE)`). The two do not produce
#' numerically identical imputations.
#'
#' **Custom `imputer` / `cd_fit`.** Supply your own imputation or
#' causal-discovery routine by passing a function with the signature
#' described above; the return value is validated and a descriptive error is
#' raised on violation. This replaces the abstract base classes
#' (`BaseMultipleImputation`, `BaseMultiGroupCDModel`) of the Python original.
#'
#' **Downstream analysis.** The result's shape (an extra `n_repeats`
#' dimension for `adjacency_matrices` and `imputation_results`) differs from
#' [lingam_direct_bootstrap()]'s `BootstrapResult`, so it cannot be passed
#' directly to [get_probabilities()] etc. Use [as_bootstrap_result()] to
#' collapse the `n_repeats` dimension (aggregating by median or mean) into a
#' `BootstrapResult`.
#'
#' **On iteration failures:** each iteration is wrapped in `tryCatch()`; a
#' failing iteration (e.g. `mice` fails to converge on a particular resample)
#' is skipped with a warning, and only if every iteration fails is an error
#' raised, mirroring [lingam_direct_bootstrap()].
#'
#' **Sequential execution only.** Unlike [lingam_direct_bootstrap()], this
#' function does not support `parallel = TRUE`; the upstream Python
#' implementation is sequential as well. If needed in the future, it can be
#' parallelized following the `parallel::makePSOCKcluster()` pattern used by
#' [lingam_direct_bootstrap()].
#' @export
#' @examples
#' set.seed(1)
#' sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
#' X <- sample6$data
#' X$x5[sample.int(nrow(X), size = round(0.1 * nrow(X)))] <- NA # MCAR 10% on x5
#'
#' \donttest{
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   res <- bootstrap_with_imputation(X,
#'     n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
#'   )
#'   print(res)
#'
#'   # Collapse the n_repeats dimension to reuse the existing bootstrap tooling
#'   bs <- as_bootstrap_result(res, aggregate = "median")
#'   get_probabilities(bs)
#' }
#' }
bootstrap_with_imputation <- function(X,
                                      n_sampling,
                                      n_repeats = 10L,
                                      imputer = NULL,
                                      cd_fit = NULL,
                                      prior_knowledge = NULL,
                                      apply_prior_knowledge_softly = FALSE,
                                      seed = NULL,
                                      verbose = TRUE) {
  var_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 3) stop("X must have at least 3 observations (rows).", call. = FALSE)
  if (!anyNA(X)) {
    warning(
      "X has no missing values; consider using lingam_direct_bootstrap() instead. ",
      "Proceeding anyway.",
      call. = FALSE
    )
  }

  n_sampling <- suppressWarnings(as.integer(n_sampling))
  if (length(n_sampling) != 1 || is.na(n_sampling) || n_sampling <= 0) {
    stop("n_sampling must be a positive integer.", call. = FALSE)
  }
  n_repeats <- suppressWarnings(as.integer(n_repeats))
  if (length(n_repeats) != 1 || is.na(n_repeats) || n_repeats <= 0) {
    stop("n_repeats must be a positive integer.", call. = FALSE)
  }
  if (!is.null(imputer) && !is.function(imputer)) {
    stop("imputer must be a function(X_boot) or NULL.", call. = FALSE)
  }
  if (!is.null(imputer) && length(formals(imputer)) < 1) {
    stop("imputer must accept at least one argument (X_boot).", call. = FALSE)
  }
  if (!is.null(cd_fit) && !is.function(cd_fit)) {
    stop("cd_fit must be a function(X_list) or NULL.", call. = FALSE)
  }
  if (!is.null(cd_fit) && length(formals(cd_fit)) < 1) {
    stop("cd_fit must accept at least one argument (X_list).", call. = FALSE)
  }
  if (!is.null(cd_fit) && (!is.null(prior_knowledge) || isTRUE(apply_prior_knowledge_softly))) {
    warning(
      "prior_knowledge / apply_prior_knowledge_softly are ignored because a custom cd_fit was supplied.",
      call. = FALSE
    )
  }
  if (is.null(cd_fit) && !is.null(prior_knowledge)) {
    Aknw <- as.matrix(prior_knowledge)
    if (!all(dim(Aknw) == c(ncol(X), ncol(X)))) {
      stop("The shape of prior knowledge must be (n_features, n_features)", call. = FALSE)
    }
  }

  n <- nrow(X)
  p <- ncol(X)

  if (!is.null(seed)) set.seed(seed)

  if (verbose) {
    message(sprintf("Bootstrap with imputation: %d iterations, n_repeats=%s (sequential)",
      n_sampling, if (is.null(imputer)) as.character(n_repeats) else "determined by imputer"
    ))
    t_start <- proc.time()
  }

  results <- vector("list", n_sampling)
  effective_n_repeats <- NULL

  for (i in seq_len(n_sampling)) {
    if (verbose && (i %% 10 == 0 || i == 1)) {
      message(sprintf("  iteration %d / %d", i, n_sampling))
    }
    idx <- sample(n, replace = TRUE)
    X_boot <- X[idx, , drop = FALSE]

    # Only the estimation calls themselves (imputer / cd_fit fitting, e.g. a
    # resample on which mice fails to converge) are treated as recoverable,
    # per-iteration stochastic failures. Validating their return values
    # (below, outside this tryCatch) is a contract check: a violation means
    # the hook is broken and must abort the whole call, not just this
    # iteration.
    imputer_attempt <- tryCatch({
      ds <- if (is.null(imputer)) default_imputer(X_boot, n_repeats) else imputer(X_boot)
      list(ok = TRUE, datasets = ds)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))

    if (!imputer_attempt$ok) {
      results[[i]] <- list(ok = FALSE, iteration = i, message = imputer_attempt$message)
      next
    }
    datasets <- validate_imputer_output(imputer_attempt$datasets, X_boot)
    this_n_repeats <- length(datasets)
    if (!is.null(effective_n_repeats) && this_n_repeats != effective_n_repeats) {
      contract_violation(sprintf(
        "imputer returned %d datasets in this iteration, but %d in a previous iteration; the count must be consistent across bootstrap iterations.",
        this_n_repeats, effective_n_repeats
      ))
    }

    cd_fit_attempt <- tryCatch({
      res <- if (is.null(cd_fit)) {
        default_cd_fit(datasets, prior_knowledge, apply_prior_knowledge_softly)
      } else {
        cd_fit(datasets)
      }
      list(ok = TRUE, cd_res = res)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))

    if (!cd_fit_attempt$ok) {
      results[[i]] <- list(ok = FALSE, iteration = i, message = cd_fit_attempt$message)
      next
    }
    cd_res <- cd_fit_attempt$cd_res
    validate_cd_fit_output(cd_res, p, this_n_repeats)

    pos <- is.na(X_boot)
    imp_res <- array(NA_real_, dim = c(this_n_repeats, n, p))
    for (r in seq_len(this_n_repeats)) {
      slice <- matrix(NA_real_, nrow = n, ncol = p)
      slice[pos] <- datasets[[r]][pos]
      imp_res[r, , ] <- slice
    }

    results[[i]] <- list(
      ok = TRUE,
      idx = idx,
      causal_order = as.integer(cd_res$causal_order),
      adjacency_matrices = cd_res$adjacency_matrices,
      imputation_results = imp_res,
      n_repeats = this_n_repeats
    )
    if (is.null(effective_n_repeats)) effective_n_repeats <- this_n_repeats
  }

  results <- filter_bootstrap_failures(results)
  n_success <- length(results)
  n_repeats_final <- effective_n_repeats

  causal_orders <- matrix(0L, nrow = n_success, ncol = p)
  adjacency_matrices <- array(0, dim = c(n_success, n_repeats_final, p, p))
  resampled_indices <- matrix(0L, nrow = n_success, ncol = n)
  imputation_results <- array(NA_real_, dim = c(n_success, n_repeats_final, n, p))

  for (i in seq_len(n_success)) {
    r <- results[[i]]
    causal_orders[i, ] <- r$causal_order
    resampled_indices[i, ] <- r$idx
    for (k in seq_len(n_repeats_final)) {
      adjacency_matrices[i, k, , ] <- as.matrix(r$adjacency_matrices[[k]])
    }
    imputation_results[i, , , ] <- r$imputation_results
  }

  if (!is.null(var_names)) {
    dimnames(adjacency_matrices) <- list(NULL, NULL, var_names, var_names)
    dimnames(imputation_results) <- list(NULL, NULL, NULL, var_names)
    # causal_orders holds variable *indices* per order position (not variable
    # identities aligned to var_names), so it intentionally gets no dimnames.
  }

  if (verbose) bootstrap_completion_message(t_start, n_success, n_sampling)

  result <- list(
    causal_orders = causal_orders,
    adjacency_matrices = adjacency_matrices,
    resampled_indices = resampled_indices,
    imputation_results = imputation_results,
    n_missing = sum(is.na(X))
  )
  class(result) <- "ImputationBootstrapResult"
  result
}


#' Print method for ImputationBootstrapResult
#'
#' @param x ImputationBootstrapResult object
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
#' @examples
#' set.seed(1)
#' sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
#' X <- sample6$data
#' X$x5[sample.int(nrow(X), size = 30)] <- NA
#'
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   res <- bootstrap_with_imputation(X,
#'     n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
#'   )
#'   print(res)
#' }
print.ImputationBootstrapResult <- function(x, ...) {
  n_sampling <- dim(x$adjacency_matrices)[1]
  n_repeats <- dim(x$adjacency_matrices)[2]
  n_features <- dim(x$adjacency_matrices)[3]
  cat(sprintf(
    "ImputationBootstrapResult: %d samplings x %d repeats, %d features, %d missing cells (original data)\n",
    n_sampling, n_repeats, n_features, x$n_missing
  ))
  invisible(x)
}


#' Collapse an ImputationBootstrapResult into a BootstrapResult
#'
#' `bootstrap_with_imputation()`'s result carries an extra `n_repeats`
#' dimension (one causal-structure estimate per imputed dataset per bootstrap
#' iteration), which the existing bootstrap analysis functions
#' ([get_probabilities()], [get_causal_direction_counts()],
#' [get_directed_acyclic_graph_counts()], [get_causal_order_stability()], `tidy()`)
#' do not expect. This collapses that dimension by aggregating the
#' `n_repeats` adjacency matrices of each iteration into one, producing a
#' regular [lingam_direct_bootstrap()]-style `BootstrapResult`.
#'
#' @param x An `ImputationBootstrapResult`, as returned by
#'   [bootstrap_with_imputation()].
#' @param aggregate How to aggregate across the `n_repeats` dimension:
#'   `"median"` (default) or `"mean"`.
#' @return A `BootstrapResult` (see [lingam_direct_bootstrap()]) with
#'   `total_effects = NULL` (total effects are not computed by
#'   `bootstrap_with_imputation()`); calling [get_total_causal_effects()] on
#'   it raises the usual "not computed" error.
#' @export
#' @examples
#' set.seed(1)
#' sample6 <- generate_lingam_sample_6(n = 300, seed = 1)
#' X <- sample6$data
#' X$x5[sample.int(nrow(X), size = 30)] <- NA
#'
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   res <- bootstrap_with_imputation(X,
#'     n_sampling = 5L, n_repeats = 3L, seed = 42, verbose = FALSE
#'   )
#'   bs <- as_bootstrap_result(res, aggregate = "median")
#'   get_probabilities(bs)
#'
#'   # get_total_causal_effects() is not available: total effects were never computed
#'   tryCatch(get_total_causal_effects(bs), error = function(e) conditionMessage(e))
#' }
as_bootstrap_result <- function(x, aggregate = c("median", "mean")) {
  if (!inherits(x, "ImputationBootstrapResult")) {
    stop("x must be the return value of bootstrap_with_imputation().", call. = FALSE)
  }
  aggregate <- match.arg(aggregate)
  agg_fun <- if (aggregate == "median") stats::median else base::mean

  am <- x$adjacency_matrices
  adjacency_matrices <- apply(am, c(1, 3, 4), agg_fun)
  dn <- dimnames(am)
  if (!is.null(dn)) dimnames(adjacency_matrices) <- list(NULL, dn[[3]], dn[[4]])

  n_sampling <- nrow(x$resampled_indices)
  resampled_indices <- lapply(seq_len(n_sampling), function(i) x$resampled_indices[i, ])

  create_bootstrap_result(
    adjacency_matrices = adjacency_matrices,
    total_effects = NULL,
    resampled_indices = resampled_indices,
    causal_orders = x$causal_orders
  )
}
