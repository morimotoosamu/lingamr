# =============================================================================
# Causal order stability from bootstrap
# =============================================================================

#' Evaluate the stability of the causal order from bootstrap
#'
#' Aggregates the causal order (causal_order) estimated in each bootstrap
#' sample and quantifies how stable the order is. Returns the rank
#' distribution of each variable, the precedence probabilities for variable
#' pairs, and an overall stability score.
#'
#' @param result A BootstrapResult object (run with the current version)
#' @param labels A vector of variable names (if NULL, x0, x1, ... are
#'   generated automatically)
#' @return A list of class `causal_order_stability`, containing:
#' * `rank_summary`: A summary of the rank of each variable (variable,
#'   mean_rank, sd_rank, median_rank, mode_rank). Sorted in ascending order
#'   of mean_rank (from upstream). A rank of 1 is the most upstream.
#' * `precedence_matrix`: A precedence probability matrix. `P[i, j]` is the
#'   proportion of bootstrap samples in which variable i was located upstream
#'   of (before) variable j.
#' * `stability_score`: An overall stability score, from 0 (random order) to
#'   1 (order agrees across all samples). The closer the precedence
#'   probability of each variable pair is to 0/1, the higher the score.
#' * `n_sampling`: The number of bootstrap samples.
#' @importFrom stats sd median
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#' get_causal_order_stability(bs, labels = names(dat$data))
get_causal_order_stability <- function(result, labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))
  if (is.null(result$causal_orders)) {
    stop(
      "This BootstrapResult does not contain causal orders. ",
      "Please re-run lingam_direct_bootstrap() with the current version of the package.",
      call. = FALSE
    )
  }

  co <- result$causal_orders
  n_sampling <- nrow(co)
  n_features <- ncol(co)

  if (is.null(labels)) {
    labels <- paste0("x", seq_len(n_features) - 1)
  } else if (length(labels) != n_features) {
    stop(sprintf(
      "Length of 'labels' (%d) must equal the number of variables (%d).",
      length(labels), n_features
    ), call. = FALSE)
  }

  # Rank of each variable in each sample (1 = most upstream)
  rank_matrix <- matrix(0L, nrow = n_sampling, ncol = n_features)
  for (s in seq_len(n_sampling)) {
    rank_matrix[s, co[s, ]] <- seq_len(n_features)
  }

  mode_rank <- function(v) {
    tb <- table(v)
    as.integer(names(tb)[which.max(tb)])
  }

  rank_summary <- data.frame(
    variable    = labels,
    mean_rank   = colMeans(rank_matrix),
    sd_rank     = apply(rank_matrix, 2, stats::sd),
    median_rank = apply(rank_matrix, 2, stats::median),
    mode_rank   = apply(rank_matrix, 2, mode_rank)
  )
  rank_summary <- rank_summary[order(rank_summary$mean_rank), ]
  rownames(rank_summary) <- NULL

  # Precedence probability matrix P[i, j] = P(i precedes j)
  P <- matrix(0, n_features, n_features, dimnames = list(labels, labels))
  for (i in seq_len(n_features)) {
    P[i, ] <- colMeans(rank_matrix[, i] < rank_matrix)
  }
  diag(P) <- 0

  # Overall stability score: how close each pair's precedence probability is to 0/1
  ps <- P[upper.tri(P)]
  stability_score <- if (length(ps) == 0) NA_real_ else mean(2 * abs(ps - 0.5))

  out <- list(
    rank_summary      = rank_summary,
    precedence_matrix = P,
    stability_score   = stability_score,
    n_sampling        = n_sampling
  )
  class(out) <- "causal_order_stability"
  out
}


#' print method for causal_order_stability
#'
#' @param x A `causal_order_stability` object
#' @param ... Additional arguments (unused)
#' @return The input object `x`, invisibly.
#' @export
print.causal_order_stability <- function(x, ...) {
  cat("=== Causal Order Stability ===\n")
  cat(sprintf("Bootstrap samples:       %d\n", x$n_sampling))
  cat(sprintf("Overall stability score: %.3f  (0 = random, 1 = fully stable)\n\n",
              x$stability_score))

  cat("Rank summary (sorted by mean rank; 1 = most upstream):\n")
  rs <- x$rank_summary
  disp <- data.frame(
    variable    = rs$variable,
    mean_rank   = sprintf("%.2f", rs$mean_rank),
    sd_rank     = sprintf("%.2f", rs$sd_rank),
    median_rank = rs$median_rank,
    mode_rank   = rs$mode_rank
  )
  print(disp, row.names = FALSE)

  cat("\nPrecedence probability P[i, j] = P(variable i precedes j):\n")
  print(round(x$precedence_matrix, 2))
  invisible(x)
}
