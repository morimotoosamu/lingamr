# =============================================================================
# Sample data generation for LiM (LiNGAM for Mixed data)
# =============================================================================


#' Generate sample data for LiM (3 mixed variables)
#'
#' Generates a small dataset with a known causal chain of continuous and
#' binary (0/1) discrete variables:
#' \code{x1 (continuous) -> x2 (discrete) -> x3 (continuous)}.
#' Continuous variables use Laplace-distributed noise (non-Gaussian, as
#' required by LiNGAM-family methods); the discrete variable is drawn from a
#' Bernoulli distribution whose logit is a linear function of its parent.
#'
#' @param n number of samples (default: 1000)
#' @param seed random seed. If `NULL` (default), no seed is set and results
#'   are not reproducible across calls.
#' @return A list with three elements:
#' * `data`: data frame with columns `x1`, `x2`, `x3`.
#' * `adjacency_matrix`: 3x3 true adjacency matrix, following the lingamr
#'   convention (`m[to, from]`, i.e. `adjacency_matrix["x2", "x1"]` is the
#'   coefficient of the x1 -> x2 edge).
#' * `is_continuous`: logical vector `c(TRUE, FALSE, TRUE)` marking which
#'   columns of `data` are continuous.
#'
#' @examples
#' dat <- generate_lim_sample(n = 500, seed = 1)
#' result <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
#' result$adjacency_matrix
#'
#' @export
generate_lim_sample <- function(n = 1000L, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1 || n < 2) {
    stop("n must be an integer >= 2.", call. = FALSE)
  }
  if (!is.null(seed) && !is.numeric(seed)) {
    stop("seed must be numeric or NULL.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.null(seed)) set.seed(as.integer(seed))

  coef_12 <- 1.3
  coef_23 <- 1.3

  # Laplace(0, 1) noise: difference of two independent Exponential(1) variables
  laplace_noise <- function(n) stats::rexp(n) - stats::rexp(n)

  x1 <- laplace_noise(n)
  eta2 <- coef_12 * x1
  x2 <- stats::rbinom(n, 1, stats::plogis(eta2))
  x3 <- coef_23 * x2 + laplace_noise(n)

  X <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  m_true <- matrix(0, 3, 3, dimnames = list(names(X), names(X)))
  m_true["x2", "x1"] <- coef_12
  m_true["x3", "x2"] <- coef_23

  list(
    data = X,
    adjacency_matrix = m_true,
    is_continuous = c(TRUE, FALSE, TRUE)
  )
}
