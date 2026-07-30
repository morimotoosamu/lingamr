# =============================================================================
# Sample data generator for CAM-UV (nonlinear model with unobserved variables)
#
# License: MIT + file LICENSE
#
# Copyright (c) 2026 O.Morimoto
# =============================================================================


#' Generate sample data with unobserved variables (for CAM-UV)
#'
#' Generates a 6-variable nonlinear additive model with two *unobserved*
#' variables, following the structure of the CAM-UV tutorial of the Python
#' `lingam` package: `u1` is an unobserved common cause of `x3` and `x4`
#' (an unobserved backdoor path, UBP), and `u2` is an unobserved
#' intermediate variable on the path from `x2` to `x5` (an unobserved
#' causal path, UCP). Only `x0`-`x5` are returned as observed data.
#'
#' The data-generating process (all error terms `e()` are
#' `runif(n, -2.5, 2.5)`; the tutorial's per-call random constants are
#' replaced by the fixed constants below for reproducibility; each column
#' is standardized before being used as a cause, as in the tutorial):
#' \preformatted{
#' u1 (latent) ~ e();  x0..x5 ~ e()
#' x3 = x3 + (u1 + 1.5)^2
#' x4 = x4 + (u1 + 1.2)^2
#' x1 = x1 + (x0 + 1.2)^2
#' x3 = x3 + (x0 - 1.5)^2
#' x4 = x4 + (x2 + 1.0)^2
#' u2 (latent) = (x2 - 1.2)^2 + e()
#' x5 = x5 + (u2 + 1.5)^2
#' }
#'
#' @param n number of samples (default: 500)
#' @param seed random seed (default: NULL, i.e. do not reset the RNG state)
#' @return list with three elements:
#' * `data`: data.frame of the 6 observed variables (`x0`-`x5`).
#' * `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
#'   variables, following the `m[to, from]` convention. Entries are 0/1
#'   edge indicators (the causal functions are nonlinear), directly
#'   comparable to the output of [lingam_camuv()]. The `x3`-`x4` (UBP) and
#'   `x2`-`x5` (UCP) entries are `NA`, matching the convention used by
#'   [lingam_camuv()].
#' * `confounded_pairs`: 2-column integer matrix of the variable pairs
#'   (1-based column positions) connected through an unobserved variable,
#'   usable as a test oracle for [lingam_camuv()]'s `confounded_pairs`.
#' @examples
#' confounded <- generate_camuv_sample(n = 200, seed = 1)
#' head(confounded$data)
#' confounded$adjacency_matrix
#' confounded$confounded_pairs
#' @export
generate_camuv_sample <- function(n = 500L, seed = NULL) {
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

  e <- function() stats::runif(n, -2.5, 2.5)
  std <- function(z) z / stats::sd(z)

  # External noise for every variable (observed and unobserved)
  x0 <- e(); x1 <- e(); x2 <- e(); x3 <- e(); x4 <- e(); x5 <- e()
  u1 <- std(e())

  # Unobserved common cause first (UBP: x3 <- u1 -> x4), as in the tutorial
  x3 <- x3 + (u1 + 1.5)^2
  x4 <- x4 + (u1 + 1.2)^2

  # Direct effects and the unobserved intermediate variable, in ascending
  # cause order; each cause column is standardized before use.
  x0 <- std(x0)
  x1 <- x1 + (x0 + 1.2)^2
  x3 <- x3 + (x0 - 1.5)^2

  x1 <- std(x1)
  x2 <- std(x2)
  x4 <- x4 + (x2 + 1.0)^2
  u2 <- std((x2 - 1.2)^2 + e()) # UCP: x2 -> u2 -> x5
  x5 <- x5 + (u2 + 1.5)^2

  x3 <- std(x3)
  x4 <- std(x4)
  x5 <- std(x5)

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x0", "x0", "x2"),
    to   = c("x1", "x3", "x4"),
    coef = c(1, 1, 1)
  )
  m_true["x2", "x5"] <- NA
  m_true["x5", "x2"] <- NA
  m_true["x3", "x4"] <- NA
  m_true["x4", "x3"] <- NA

  confounded_pairs <- rbind(c(3L, 6L), c(4L, 5L))
  colnames(confounded_pairs) <- c("var1", "var2")

  list(
    data             = X,
    adjacency_matrix = m_true,
    confounded_pairs = confounded_pairs
  )
}
