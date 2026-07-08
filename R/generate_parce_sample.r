#' Generate sample data with a latent confounder (for BottomUpParceLiNGAM)
#'
#' Generates the 7-variable model used in the ParceLiNGAM tutorial, where
#' `x6` is an unobserved (latent) common cause of `x2` and `x3`. Only
#' `x0`-`x5` are returned as observed data; `x6` is not included.
#'
#' The data-generating process (all error terms are `Uniform(0, 1)`):
#' \preformatted{
#' x6 (latent) ~ Uniform(0, 1)
#' x3 = 2.0 * x6 + e
#' x2 = 2.0 * x6 + e
#' x0 = 0.5 * x3 + e
#' x1 = 0.5 * x0 + 0.5 * x2 + e
#' x5 = 0.5 * x0 + e
#' x4 = 0.5 * x0 - 0.5 * x2 + e
#' }
#'
#' @param n number of samples (default: 1000)
#' @param seed random seed (default: NULL, i.e. do not reset the RNG state)
#' @return list with three elements:
#' * `data`: data.frame of the 6 observed variables (`x0`-`x5`).
#' * `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
#'   variables, following the `m[to, from]` convention. The `x2`-`x3` entries
#'   (which share the latent confounder `x6` and have no direct edge between
#'   them) are `NA`, matching the convention used by [lingam_parce()].
#' * `confounded_pair`: 1-based column positions of `x2` and `x3` (the pair
#'   sharing the latent confounder).
#' @examples
#' confounded <- generate_parce_sample(n = 500, seed = 42)
#' head(confounded$data)
#' confounded$adjacency_matrix
#' confounded$confounded_pair
#' @export
generate_parce_sample <- function(n = 1000L, seed = NULL) {
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

  x6 <- stats::runif(n)
  x3 <- 2.0 * x6 + stats::runif(n)
  x2 <- 2.0 * x6 + stats::runif(n)
  x0 <- 0.5 * x3 + stats::runif(n)
  x1 <- 0.5 * x0 + 0.5 * x2 + stats::runif(n)
  x5 <- 0.5 * x0 + stats::runif(n)
  x4 <- 0.5 * x0 - 0.5 * x2 + stats::runif(n)

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x3", "x0", "x0", "x2", "x0", "x2"),
    to   = c("x0", "x1", "x5", "x1", "x4", "x4"),
    coef = c(0.5, 0.5, 0.5, 0.5, 0.5, -0.5)
  )
  m_true["x2", "x3"] <- NA
  m_true["x3", "x2"] <- NA

  list(
    data             = X,
    adjacency_matrix = m_true,
    confounded_pair  = c(3L, 4L)
  )
}
