# =============================================================================
# Sample data generator for RESIT (nonlinear additive noise model)
#
# License: MIT + file LICENSE
#
# Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Generate sample data from a nonlinear additive noise model (for RESIT)
#'
#' Generates a 4-variable nonlinear SEM whose structure can only be
#' recovered by a nonlinear method such as [lingam_resit()]; linear
#' LiNGAM variants are not expected to work on this data.
#'
#' The data-generating process (all error terms `e()` are
#' `runif(n, -0.5, 0.5)`):
#' \preformatted{
#' x0 ~ runif(-1, 1)
#' x1 = 3.0 * x0^2 + e()
#' x2 = 2.0 * tanh(x1) + 0.8 * x0^3 + e()
#' x3 = 1.5 * sin(x2) + e()
#' }
#'
#' @param n number of samples (default: 300)
#' @param seed random seed (default: NULL, i.e. do not reset the RNG state)
#' @return list with three elements:
#' * `data`: data.frame of the 4 observed variables (`x0`-`x3`).
#' * `adjacency_matrix`: the true 4x4 adjacency matrix following the
#'   `m[to, from]` convention. Entries are 0/1 edge indicators (not
#'   coefficients), directly comparable to the output of [lingam_resit()].
#' * `causal_order`: the true causal order (1-based column positions,
#'   source first).
#' @examples
#' nonlinear <- generate_resit_sample(n = 300, seed = 1)
#' head(nonlinear$data)
#' nonlinear$adjacency_matrix
#' @export
generate_resit_sample <- function(n = 300L, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 2) {
    stop("n must be an integer >= 2.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("seed must be a numeric scalar or NULL.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  e <- function() stats::runif(n, -0.5, 0.5)

  x0 <- stats::runif(n, -1, 1)
  x1 <- 3.0 * x0^2 + e()
  x2 <- 2.0 * tanh(x1) + 0.8 * x0^3 + e()
  x3 <- 1.5 * sin(x2) + e()

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x0", "x0", "x1", "x2"),
    to   = c("x1", "x2", "x2", "x3"),
    coef = c(1, 1, 1, 1)
  )

  list(
    data             = X,
    adjacency_matrix = m_true,
    causal_order     = 1:4
  )
}
