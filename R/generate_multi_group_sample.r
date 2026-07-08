# =============================================================================
# Sample data generation for Multi-Group Direct LiNGAM
# Reuses the noise / true-adjacency helpers from R/generate_lingam_sample.r
# =============================================================================


#' Generate sample data for Multi-Group Direct LiNGAM (2 groups, 6 variables)
#'
#' Generates two datasets that share the same causal structure as
#' [generate_lingam_sample_6()] (`x3 -> x0, x3 -> x2, x0 -> x1, x2 -> x1,
#' x0 -> x4, x2 -> x4, x0 -> x5`) but with different structural coefficients
#' per group, following the multi-dataset tutorial's setup.
#'
#' @param n Numeric vector of length 2: sample size per group (default `c(1000, 1000)`).
#' @param seed Random seed (default 42). Group 2 uses an internally offset
#'   seed so the two groups are independently drawn.
#' @return A list with:
#' * `data_list`: named list (`group1`, `group2`) of data frames, each with
#'   columns `x0`..`x5`.
#' * `adjacency_matrices`: named list of the true adjacency matrix per group
#'   (`m[to, from] = coef`, same convention as [lingam_direct()]).
#' * `causal_order`: the true causal order shared by both groups (1-based
#'   indices into `x0`..`x5`).
#' @examples
#' mg <- generate_multi_group_sample()
#' lapply(mg$data_list, head)
#' mg$adjacency_matrices$group1
#' @export
generate_multi_group_sample <- function(n = c(1000, 1000), seed = 42L) {
  if (!is.numeric(n) || length(n) != 2) {
    stop("n must be a numeric vector of length 2 (sample size per group).", call. = FALSE)
  }
  if (any(n < 2)) stop("Each element of n must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("seed must be a single non-missing numeric value.", call. = FALSE)
  }
  if (abs(seed) > .Machine$integer.max / 20000) {
    stop("seed is too large; it must satisfy abs(seed) <= .Machine$integer.max / 20000 to avoid integer overflow in the per-group seed offset.", call. = FALSE)
  }

  n <- as.integer(n)
  seed <- as.integer(seed)

  noise_fn <- make_noise_fn("uniform")
  var_names <- paste0("x", 0:5)

  # Group 1 / group 2 structural coefficients, in the order:
  # x3->x0, x3->x2, x0->x1, x2->x1, x0->x4, x2->x4, x0->x5
  coefs <- list(
    group1 = c(3.0, 6.0, 3.0, 2.0, 8.0, -1.0, 4.0),
    group2 = c(3.5, 6.5, 3.5, 2.5, 8.5, -1.5, 4.5)
  )

  data_list <- vector("list", 2)
  adjacency_matrices <- vector("list", 2)
  group_names <- names(coefs)

  for (g in seq_along(coefs)) {
    ng <- n[g]
    # Offset the base seed per group so the two groups are drawn independently.
    seed_g <- seed + 10000L * (g - 1L)
    E <- generate_noise_matrix(ng, 6L, seed_g, noise_fn)
    cf <- coefs[[g]]

    x3 <- E[, 4]
    x0 <- cf[1] * x3 + E[, 1]
    x2 <- cf[2] * x3 + E[, 3]
    x1 <- cf[3] * x0 + cf[4] * x2 + E[, 2]
    x4 <- cf[5] * x0 + cf[6] * x2 + E[, 5]
    x5 <- cf[7] * x0 + E[, 6]

    data_list[[g]] <- data.frame(x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

    adjacency_matrices[[g]] <- build_true_adjacency(
      var_names,
      from = c("x3", "x3", "x0", "x2", "x0", "x2", "x0"),
      to   = c("x0", "x2", "x1", "x1", "x4", "x4", "x5"),
      coef = cf
    )
  }
  names(data_list) <- group_names
  names(adjacency_matrices) <- group_names

  # A valid topological order: x3, x0, x2, x1, x4, x5 (1-based indices into x0..x5)
  causal_order <- c(4L, 1L, 3L, 2L, 5L, 6L)

  list(
    data_list = data_list,
    adjacency_matrices = adjacency_matrices,
    causal_order = causal_order
  )
}
