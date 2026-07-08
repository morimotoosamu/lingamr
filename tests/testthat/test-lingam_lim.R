test_that("lingam_lim returns LiMResult with correct structure", {
  dat <- generate_lim_sample(n = 300, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous)

  expect_s3_class(res, "LiMResult")
  expect_named(res, c("adjacency_matrix", "causal_order", "is_continuous"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(3L, 3L))
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_equal(rownames(res$adjacency_matrix), names(dat$data))
  expect_equal(length(res$causal_order), 3L)
  expect_equal(res$is_continuous, dat$is_continuous)
})

test_that("lingam_lim errors on invalid inputs", {
  dat <- generate_lim_sample(n = 100, seed = 1)

  expect_error(
    lingam_lim(matrix(letters[1:6], nrow = 2), is_continuous = c(TRUE, TRUE, TRUE)),
    "numeric matrix"
  )
  expect_error(
    lingam_lim(matrix(c(1, NA, 2, 3), nrow = 2), is_continuous = c(TRUE, TRUE)),
    "missing values"
  )
  expect_error(
    lingam_lim(as.data.frame(matrix(1, nrow = 5, ncol = 1)), is_continuous = TRUE),
    "at least 2 variables"
  )
  expect_error(
    lingam_lim(matrix(numeric(0), nrow = 0, ncol = 3), is_continuous = c(TRUE, TRUE, TRUE)),
    "at least 2 observations"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = c(TRUE, FALSE)),
    "length equal to ncol"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = c(TRUE, FALSE, NA)),
    "without NA"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = c(TRUE, 1, 0)),
    "logical vector"
  )
  # x1 is continuous in dat$data, so declaring it discrete should fail the
  # binary-value check (it is not 0/1-valued)
  expect_error(
    lingam_lim(dat$data, is_continuous = c(FALSE, FALSE, TRUE)),
    "binary \\(0/1\\)"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = dat$is_continuous, lambda1 = -1),
    "lambda1"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = dat$is_continuous, max_iter = 0),
    "max_iter"
  )
  expect_error(
    lingam_lim(dat$data, is_continuous = dat$is_continuous, only_global = "yes"),
    "only_global"
  )
})

test_that("print.LiMResult runs without error", {
  dat <- generate_lim_sample(n = 200, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous)

  expect_output(print(res), "LiM Result")
  expect_output(print(res), "Variable types")
  expect_output(print(res), "Causal order")
  expect_output(print(res), "Adjacency matrix")
})

test_that("lingam_lim is reproducible given the same seed", {
  dat <- generate_lim_sample(n = 300, seed = 1)

  set.seed(42)
  r1 <- lingam_lim(dat$data, is_continuous = dat$is_continuous)
  set.seed(42)
  r2 <- lingam_lim(dat$data, is_continuous = dat$is_continuous)

  expect_equal(r1, r2)
})

test_that("lingam_lim recovers the correct edge direction (transpose check)", {
  # Known structure: x1 (continuous) -> x2 (discrete) -> x3 (continuous).
  # This is the single most important test: a transpose bug in lingam_lim()
  # would silently swap adjacency_matrix["x2", "x1"] and
  # adjacency_matrix["x1", "x2"].
  dat <- generate_lim_sample(n = 2000, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous)

  B <- res$adjacency_matrix
  expect_true(B["x2", "x1"] != 0)
  expect_lt(abs(B["x1", "x2"]), 1e-6)
  expect_true(B["x3", "x2"] != 0)
  expect_lt(abs(B["x2", "x3"]), 1e-6)
})

test_that("only_global = TRUE returns a valid LiMResult", {
  dat <- generate_lim_sample(n = 300, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous, only_global = TRUE)

  expect_s3_class(res, "LiMResult")
  expect_equal(dim(res$adjacency_matrix), c(3L, 3L))
  expect_equal(length(res$causal_order), 3L)
})

test_that("an extreme w_threshold drops all edges", {
  # Local search can re-introduce genuine edges even when the global phase
  # is thresholded to zero (it evaluates BIC on candidate structures rather
  # than reusing the thresholded coefficients), so this check is scoped to
  # only_global = TRUE, which reflects the thresholded coefficients directly.
  dat <- generate_lim_sample(n = 2000, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous,
                     w_threshold = 10, only_global = TRUE)

  expect_equal(sum(res$adjacency_matrix != 0), 0L)
  expect_setequal(res$causal_order, 1:3)
})

test_that("causal_order is consistent with the adjacency matrix's non-zero pattern", {
  dat <- generate_lim_sample(n = 300, seed = 1)
  res <- lingam_lim(dat$data, is_continuous = dat$is_continuous)

  B <- res$adjacency_matrix
  order_pos <- match(seq_len(ncol(B)), res$causal_order)
  nz <- which(B != 0, arr.ind = TRUE)  # nz[, "row"] = to, nz[, "col"] = from
  if (nrow(nz) > 0) {
    # every parent ("from") must appear earlier in the causal order than its child ("to")
    expect_true(all(order_pos[nz[, "col"]] < order_pos[nz[, "row"]]))
  }
})

test_that("lingam_lim produces no console output on the normal path", {
  dat <- generate_lim_sample(n = 300, seed = 1)
  expect_silent(lingam_lim(dat$data, is_continuous = dat$is_continuous))
})

test_that("lim_likelihood_i returns a large finite penalty instead of Inf for a zero-variance residual", {
  set.seed(21)
  n <- 50
  X <- cbind(x1 = rnorm(n), x2 = rep(0, n)) # x2 perfectly explained by b_full below
  b_full <- c(0, 0)
  b0 <- 0
  # residual for column 2 is X[,2] - b0 == 0 everywhere -> var_i == 0
  res <- lim_likelihood_i(X, i = 2, b_full = b_full, b0 = b0)
  expect_true(is.finite(res))
  expect_lt(res, -1e9)
})

test_that("lim_bic_loss stays finite when a discrete component is quasi-separated", {
  set.seed(22)
  n <- 40
  parent <- rnorm(n)
  child <- as.numeric(parent > 0) # perfectly separated by the parent
  X <- cbind(x1 = parent, x2 = child)
  W <- matrix(0, 2, 2)
  W[1, 2] <- 1 # x1 -> x2

  loss <- lim_bic_loss(W, X, is_continuous = c(TRUE, FALSE))
  expect_true(is.finite(loss))
})

test_that("lim_bic_loss stays finite when a continuous component has a rank-deficient parent design", {
  set.seed(23)
  n <- 40
  x1 <- rnorm(n)
  X <- cbind(x1 = x1, x2 = x1, x3 = rnorm(n)) # x2 is a duplicate of x1
  W <- matrix(0, 3, 3)
  W[1, 3] <- 1
  W[2, 3] <- 1 # x3's parents (x1, x2) are collinear

  loss <- lim_bic_loss(W, X, is_continuous = c(TRUE, TRUE, TRUE))
  expect_true(is.finite(loss))
})
