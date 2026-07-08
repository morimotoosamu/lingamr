#' Generate sample data with a latent confounder (for RCD)
#'
#' Generates the 7-variable model used in the RCD tutorial, where `x6` is an
#' unobserved (latent) common cause of `x2` and `x4`. Only `x0`-`x5` are
#' returned as observed data; `x6` is not included.
#'
#' The data-generating process (all error terms `e()` are super-Gaussian,
#' `rnorm(n, 0, 0.5)^3`):
#' \preformatted{
#' x5 ~ e();  x6 (latent) ~ e()
#' x1 = 0.6 * x5 + e()
#' x3 = 0.5 * x5 + e()
#' x0 = 1.0 * x1 + 1.0 * x3 + e()
#' x2 = 0.8 * x0 - 0.6 * x6 + e()
#' x4 = 1.0 * x0 - 0.5 * x6 + e()
#' }
#'
#' @param n number of samples (default: 300)
#' @param seed random seed (default: NULL, i.e. do not reset the RNG state)
#' @return list with four elements:
#' * `data`: data.frame of the 6 observed variables (`x0`-`x5`).
#' * `adjacency_matrix`: the true 6x6 adjacency matrix among the observed
#'   variables, following the `m[to, from]` convention. The `x2`-`x4` entries
#'   (which share the latent confounder `x6` and have no direct edge between
#'   them) are `NA`, matching the convention used by [lingam_rcd()].
#' * `ancestors_list`: the true ancestor sets (1-based column positions),
#'   usable as a test oracle for [lingam_rcd()]'s `ancestors_list`.
#' * `confounded_pair`: 1-based column positions of `x2` and `x4` (the pair
#'   sharing the latent confounder).
#' @examples
#' confounded <- generate_rcd_sample(n = 300, seed = 1)
#' head(confounded$data)
#' confounded$adjacency_matrix
#' confounded$ancestors_list
#' confounded$confounded_pair
#' @export
generate_rcd_sample <- function(n = 300L, seed = NULL) {
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

  e <- function() stats::rnorm(n, 0, 0.5)^3

  x5 <- e()
  x6 <- e()
  x1 <- 0.6 * x5 + e()
  x3 <- 0.5 * x5 + e()
  x0 <- 1.0 * x1 + 1.0 * x3 + e()
  x2 <- 0.8 * x0 - 0.6 * x6 + e()
  x4 <- 1.0 * x0 - 0.5 * x6 + e()

  X <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

  m_true <- build_true_adjacency(
    names(X),
    from = c("x1", "x3", "x5", "x0", "x0"),
    to   = c("x0", "x0", "x1", "x2", "x4"),
    coef = c(1.0, 1.0, 0.6, 0.8, 1.0)
  )
  m_true["x2", "x4"] <- NA
  m_true["x4", "x2"] <- NA

  ancestors_true <- list(
    x0 = c(2L, 4L, 6L),
    x1 = 6L,
    x2 = c(1L, 2L, 4L, 6L),
    x3 = 6L,
    x4 = c(1L, 2L, 4L, 6L),
    x5 = integer(0)
  )
  names(ancestors_true) <- names(X)

  list(
    data             = X,
    adjacency_matrix = m_true,
    ancestors_list   = ancestors_true,
    confounded_pair  = c(3L, 5L)
  )
}
