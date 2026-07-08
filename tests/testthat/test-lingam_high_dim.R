test_that("lingam_high_dim returns LingamResult with correct structure", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_high_dim(dat$data)

  expect_s3_class(res, "LingamResult")
  expect_named(res, c("adjacency_matrix", "causal_order"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_equal(rownames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_high_dim errors on invalid inputs", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

  expect_error(lingam_high_dim(matrix(letters[1:4], nrow = 2)), "numeric matrix")
  expect_error(
    lingam_high_dim(as.data.frame(matrix(1, nrow = 5, ncol = 1))),
    "at least 2 variables"
  )
  expect_error(
    lingam_high_dim(matrix(numeric(0), nrow = 0, ncol = 3)),
    "at least 2 observations"
  )
  mat_na <- as.matrix(dat$data)
  mat_na[1, 1] <- NA
  expect_error(lingam_high_dim(mat_na), "missing values")

  expect_error(lingam_high_dim(dat$data, J = 2), "J must be")
  expect_error(lingam_high_dim(dat$data, K = 0), "K must be")
  expect_error(lingam_high_dim(dat$data, alpha = 1.5), "alpha must be")
  expect_error(lingam_high_dim(dat$data, estimate_adj_mat = "yes"), "estimate_adj_mat must be")
})

test_that("lingam_high_dim recovers causal order and edge directions", {
  # True structure: x3 -> x0 (3.0), x3 -> x2 (6.0), x0 -> x1 (3.0),
  # x2 -> x1 (2.0), x0 -> x5 (4.0), x0 -> x4 (8.0), x2 -> x4 (-1.0)
  dat <- generate_lingam_sample_6(n = 1000, seed = 1, noise_dist = "uniform")
  res <- lingam_high_dim(dat$data)

  order_pos <- function(var_name) which(res$causal_order == which(names(dat$data) == var_name))

  # every true edge j -> i has j appearing before i in the causal order
  edges <- list(
    c("x3", "x0"), c("x3", "x2"), c("x0", "x1"),
    c("x2", "x1"), c("x0", "x5"), c("x0", "x4"), c("x2", "x4")
  )
  for (e in edges) {
    expect_lt(order_pos(e[1]), order_pos(e[2]))
  }

  # true edges are non-zero in the adjacency matrix with the correct sign
  B <- res$adjacency_matrix
  idx <- function(var_name) which(names(dat$data) == var_name)
  expect_gt(B[idx("x0"), idx("x3")], 0)
  expect_gt(B[idx("x2"), idx("x3")], 0)
  expect_gt(B[idx("x1"), idx("x0")], 0)
  expect_gt(B[idx("x1"), idx("x2")], 0)
  expect_gt(B[idx("x5"), idx("x0")], 0)
  expect_gt(B[idx("x4"), idx("x0")], 0)
  expect_lt(B[idx("x4"), idx("x2")], 0)

  # false edges (all off-diagonal positions not among the true edges above)
  # should carry negligible weight
  true_edge_linear_indices <- vapply(edges, function(e) {
    idx(e[2]) + (idx(e[1]) - 1L) * nrow(B)
  }, numeric(1))
  off_diag_positions <- which(row(B) != col(B))
  non_true_positions <- setdiff(off_diag_positions, true_edge_linear_indices)
  expect_lt(sum(abs(B[non_true_positions])), 1)
})

test_that("lingam_high_dim is deterministic for the n > p route", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res1 <- lingam_high_dim(dat$data)
  res2 <- lingam_high_dim(dat$data)

  expect_equal(res1, res2)
})

test_that("estimate_adj_mat = FALSE returns an NA adjacency matrix but the same causal order", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res_full <- lingam_high_dim(dat$data)
  res_order_only <- lingam_high_dim(dat$data, estimate_adj_mat = FALSE)

  expect_true(all(is.na(res_order_only$adjacency_matrix)))
  expect_equal(dim(res_order_only$adjacency_matrix), c(6L, 6L))
  expect_equal(colnames(res_order_only$adjacency_matrix), names(dat$data))
  expect_equal(res_order_only$causal_order, res_full$causal_order)
})

test_that("lingam_high_dim falls back to cross-validated lasso when n <= p", {
  skip_if_not_installed("glmnet")
  skip_on_cran()

  wide <- generate_lingam_large_sample(p = 45, n = 40, seed = 1)

  set.seed(1)
  expect_warning(
    res <- lingam_high_dim(wide$data),
    "n_samples <= n_features"
  )

  expect_s3_class(res, "LingamResult")
  expect_equal(dim(res$adjacency_matrix), c(45L, 45L))
  expect_equal(length(res$causal_order), 45L)
  expect_false(anyNA(res$adjacency_matrix))
})

test_that("calc_taus minimizes over every conditioning subset, not just the first", {
  # Regression test for a bug in the upstream Python implementation
  # (cdt15/lingam) where a mis-indented `return` inside the loop over
  # conditioning sets causes only the first subset to ever be evaluated.
  # This R port deliberately does NOT replicate that bug, so calc_taus()
  # must return the true minimum across ALL supplied conditioning sets.
  set.seed(30)
  n <- 200
  Y <- matrix(rnorm(n * 4), n, 4)
  yty <- t(Y) %*% Y

  # Two conditioning sets on the same (pa, ch) pair will generically produce
  # different tau values; the true minimum-over-both must differ from (and
  # be <=) the value obtained from evaluating only the first set.
  cond_sets <- list(2L, 3L)
  an_sets <- list(integer(0), integer(0))

  ret_both <- calc_taus(Y, yty, pa = 1L, ch = 4L, k = 3L, cond_sets = cond_sets, an_sets = an_sets)
  ret_first_only <- calc_taus(Y, yty, pa = 1L, ch = 4L, k = 3L, cond_sets = cond_sets[1], an_sets = an_sets[1])
  ret_second_only <- calc_taus(Y, yty, pa = 1L, ch = 4L, k = 3L, cond_sets = cond_sets[2], an_sets = an_sets[2])

  expect_equal(ret_both[4L], min(ret_first_only[4L], ret_second_only[4L]))
  # sanity: the two single-subset results actually differ, so this test
  # would fail under the old first-subset-only behavior
  expect_false(isTRUE(all.equal(ret_first_only[4L], ret_second_only[4L])))
})

test_that("glmnet absence errors on the n <= p route but not the causal-order-only route", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)
  wide <- generate_lingam_large_sample(p = 20, n = 15, seed = 1)

  local_mocked_bindings(
    check_glmnet_available = function(method) {
      stop(sprintf(
        "Package 'glmnet' is required for reg_method = \"%s\". Please install it.",
        method
      ), call. = FALSE)
    }
  )

  # causal-order search itself never calls glmnet
  expect_no_error(lingam_high_dim(dat$data, estimate_adj_mat = FALSE))

  # n > p with adjacency estimation still needs glmnet (adaptive_lasso backend)
  expect_error(lingam_high_dim(dat$data), "glmnet")

  # n <= p route needs glmnet (cv.glmnet backend)
  expect_error(suppressWarnings(lingam_high_dim(wide$data)), "glmnet")
})
