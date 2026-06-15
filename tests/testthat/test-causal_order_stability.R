test_that("BootstrapResult stores causal_orders", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L,
                                 reg_method = "ols")

  expect_false(is.null(bs$causal_orders))
  expect_equal(dim(bs$causal_orders), c(15L, 6L))
  # Each row is a permutation of 1:6
  expect_true(all(apply(bs$causal_orders, 1, function(r) setequal(r, 1:6))))
})

test_that("get_causal_order_stability returns expected structure", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 20L, seed = 42L,
                                 reg_method = "ols")
  st  <- get_causal_order_stability(bs)

  expect_s3_class(st, "causal_order_stability")
  expect_true(is.matrix(st$precedence_matrix))
  expect_equal(dim(st$precedence_matrix), c(6L, 6L))

  off_diag <- st$precedence_matrix[row(st$precedence_matrix) != col(st$precedence_matrix)]
  expect_true(all(off_diag >= 0 & off_diag <= 1))
  expect_true(st$stability_score >= 0 && st$stability_score <= 1)
  expect_true(all(c("variable", "mean_rank", "sd_rank", "median_rank", "mode_rank")
                  %in% names(st$rank_summary)))
  expect_equal(nrow(st$rank_summary), 6L)
})

test_that("get_causal_order_stability: x3 is most upstream", {
  dat <- generate_lingam_sample_6(n = 1000, seed = 42)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 1L,
                                 reg_method = "ols")
  st  <- get_causal_order_stability(bs, labels = names(dat$data))

  # x3 is the root, so its mean rank is the smallest (first after sorting)
  expect_equal(st$rank_summary$variable[1], "x3")
})

test_that("precedence_matrix P[i,j] + P[j,i] == 1 for i != j", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, seed = 42L,
                                 reg_method = "ols")
  st  <- get_causal_order_stability(bs)
  P   <- st$precedence_matrix

  for (i in 1:5) {
    for (j in (i + 1):6) {
      expect_equal(P[i, j] + P[j, i], 1)
    }
  }
})

test_that("get_causal_order_stability errors on legacy result without causal_orders", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 1L,
                                 reg_method = "ols")
  bs$causal_orders <- NULL  # Simulate a result from an older version

  expect_error(get_causal_order_stability(bs), "causal order")
})

test_that("get_causal_order_stability validates labels length", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 1L,
                                 reg_method = "ols")

  expect_error(get_causal_order_stability(bs, labels = c("a", "b")))
})

test_that("print.causal_order_stability runs", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 10L, seed = 1L,
                                 reg_method = "ols")
  st  <- get_causal_order_stability(bs)

  expect_output(print(st), "Causal Order Stability")
  expect_output(print(st), "stability score")
})
