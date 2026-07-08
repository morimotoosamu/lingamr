test_that("hsic_test_gamma gives a large p-value for independent data", {
  set.seed(1)
  n <- 300
  x <- rnorm(n)
  y <- rnorm(n)

  res <- hsic_test_gamma(x, y)
  expect_true(is.finite(res$stat))
  expect_true(is.finite(res$p))
  expect_gt(res$p, 0.05)
})

test_that("hsic_test_gamma gives a near-zero p-value for strongly dependent data", {
  set.seed(2)
  n <- 300
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.01)

  res <- hsic_test_gamma(x, y)
  expect_lt(res$p, 1e-6)
})

test_that("hsic_test_gamma is deterministic for the same input", {
  set.seed(3)
  n <- 200
  x <- rnorm(n)
  y <- 0.5 * x + rnorm(n)

  res1 <- hsic_test_gamma(x, y)
  res2 <- hsic_test_gamma(x, y)

  expect_identical(res1, res2)
})

test_that("hsic_kernel_width and hsic_gram_matrix produce a valid symmetric Gram matrix", {
  set.seed(4)
  x <- rnorm(150)
  width <- hsic_kernel_width(x)
  expect_gt(width, 0)

  gm <- hsic_gram_matrix(x, width)
  expect_equal(dim(gm$K), c(150L, 150L))
  expect_equal(gm$K, t(gm$K))
  expect_true(all(diag(gm$K) == 1))
  # centered matrix has (approximately) zero row/column sums
  expect_lt(max(abs(rowSums(gm$Kc))), 1e-8)
})

test_that("f_correlation is near 0 for independent data and large for dependent data", {
  set.seed(5)
  n <- 300
  x <- rnorm(n)
  y_ind <- rnorm(n)
  y_dep <- 0.9 * x + rnorm(n, sd = 0.1)

  fc_ind <- f_correlation(x, y_ind)
  fc_dep <- f_correlation(x, y_dep)

  expect_gt(fc_dep, fc_ind)
  expect_true(fc_ind >= 0 && fc_ind <= 1.01)
  expect_true(fc_dep >= 0 && fc_dep <= 1.01)
})

test_that("f_correlation is deterministic for the same input", {
  set.seed(6)
  n <- 150
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n)

  res1 <- f_correlation(x, y)
  res2 <- f_correlation(x, y)

  expect_identical(res1, res2)
})

test_that("hsic_test_gamma errors on mismatched lengths", {
  expect_error(hsic_test_gamma(1:10, 1:9), "same length")
})

test_that("hsic_test_gamma errors when n < 6", {
  set.seed(8)
  x <- rnorm(5)
  y <- rnorm(5)
  expect_error(hsic_test_gamma(x, y), "at least")
})

test_that("hsic_test_gamma errors on NA input", {
  set.seed(9)
  x <- rnorm(10)
  y <- rnorm(10)
  x[1] <- NA
  expect_error(hsic_test_gamma(x, y), "NA/NaN")
})

test_that("hsic_test_gamma treats a constant input as trivially independent", {
  set.seed(10)
  n <- 50
  x <- rep(1, n)
  y <- rnorm(n)

  res <- hsic_test_gamma(x, y)
  expect_equal(res$stat, 0)
  expect_equal(res$p, 1)
})

test_that("f_correlation errors on mismatched lengths", {
  expect_error(f_correlation(1:10, 1:9), "same length")
})

test_that("f_correlation errors on NA input", {
  x <- rnorm(20)
  y <- rnorm(20)
  y[3] <- NA
  expect_error(f_correlation(x, y), "NA/NaN")
})

test_that("f_correlation returns 0 for a constant input instead of NaN/Inf", {
  set.seed(11)
  n <- 50
  x <- rep(2, n)
  y <- rnorm(n)

  expect_equal(f_correlation(x, y), 0)
  expect_equal(f_correlation(y, x), 0)
  expect_equal(f_correlation(x, x), 0)
})

test_that("incomplete_cholesky_fcorr reconstructs the Gaussian Gram matrix reasonably well", {
  set.seed(7)
  n <- 200
  sigma <- 1.0
  x <- rnorm(n)
  kappa <- 2e-2
  tol <- n * kappa * 1e-2

  G <- incomplete_cholesky_fcorr(x, sigma, tol)
  K_true <- exp(-1 / (2 * sigma^2) * outer(x, x, "-")^2)
  K_approx <- tcrossprod(G)

  expect_lte(ncol(G), n)
  expect_lt(mean(abs(K_true - K_approx)), 0.05)
})
