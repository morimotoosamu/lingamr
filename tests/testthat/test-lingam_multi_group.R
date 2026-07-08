test_that("lingam_multi_group returns MultiGroupLingamResult with expected structure", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  expect_s3_class(res, "MultiGroupLingamResult")
  expect_named(res, c("adjacency_matrices", "causal_order"))
  expect_named(res$adjacency_matrices, c("group1", "group2"))
  for (B in res$adjacency_matrices) {
    expect_equal(dim(B), c(6L, 6L))
    expect_false(is.null(dimnames(B)))
  }
  expect_length(res$causal_order, 6L)
})

test_that("lingam_multi_group validates X_list", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)$data

  expect_error(lingam_multi_group(dat), "X_list must be a list")
  expect_error(lingam_multi_group(list(dat)), "at least two items")
  expect_error(
    lingam_multi_group(list(dat, dat[, 1:3])),
    "same number of columns"
  )
  bad_na <- dat
  bad_na[1, 1] <- NA
  expect_error(lingam_multi_group(list(dat, bad_na)), "missing values")
  bad_char <- dat
  bad_char$x0 <- as.character(bad_char$x0)
  expect_error(lingam_multi_group(list(dat, bad_char)), "numeric matrix")
})

test_that("duplicating a single group matches lingam_direct's causal order", {
  dat <- generate_lingam_sample_6(n = 500, seed = 3)$data

  single <- lingam_direct(dat, reg_method = "ols")
  joint <- lingam_multi_group(list(dat, dat), reg_method = "ols")

  expect_equal(joint$causal_order, single$causal_order)
})

test_that("lingam_multi_group recovers the true structure", {
  mg <- generate_multi_group_sample(n = c(1000, 1000), seed = 42)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  # Every true edge j -> i has j preceding i in the estimated common order
  true_order_pos <- function(idx) which(res$causal_order == idx)
  for (g in names(mg$adjacency_matrices)) {
    true_B <- mg$adjacency_matrices[[g]]
    edges <- which(true_B != 0, arr.ind = TRUE)
    for (r in seq_len(nrow(edges))) {
      to_i <- edges[r, 1]
      from_j <- edges[r, 2]
      expect_lt(true_order_pos(from_j), true_order_pos(to_i))
    }
    est_B <- res$adjacency_matrices[[g]]
    for (r in seq_len(nrow(edges))) {
      to_i <- edges[r, 1]
      from_j <- edges[r, 2]
      expect_equal(est_B[to_i, from_j], true_B[to_i, from_j], tolerance = 0.5)
    }
  }
})

test_that("get_group_result extracts a usable LingamResult", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  g1 <- get_group_result(res, 1)
  expect_s3_class(g1, "LingamResult")
  expect_equal(g1$adjacency_matrix, res$adjacency_matrices$group1)

  g2 <- get_group_result(res, "group2")
  expect_equal(g2$adjacency_matrix, res$adjacency_matrices$group2)

  expect_no_error(estimate_all_total_effects(mg$data_list$group1, g1, method = "ols"))
  expect_no_error(get_error_independence_p_values(mg$data_list$group1, g1))

  expect_error(get_group_result(res, "no_such_group"), "not found")
  expect_error(get_group_result(res, 99), "between 1 and")
})

test_that("lingam_multi_group works with prior knowledge", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)
  pk <- make_prior_knowledge(6, exogenous_variables = 4) # x3 (1-based index 4) is exogenous

  res <- lingam_multi_group(mg$data_list, prior_knowledge = pk, reg_method = "ols")
  expect_s3_class(res, "MultiGroupLingamResult")
  expect_equal(res$causal_order[1], 4L)
})

test_that("print.MultiGroupLingamResult reports the expected header", {
  mg <- generate_multi_group_sample(n = c(200, 200), seed = 1)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  expect_output(print(res), "Multi-Group")
})

test_that("lingam_multi_group_bootstrap returns MultiGroupBootstrapResult", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)
  bs <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 10L, reg_method = "ols", seed = 42L, verbose = FALSE
  )

  expect_s3_class(bs, "MultiGroupBootstrapResult")
  expect_length(bs, 2L)
  expect_named(bs, c("group1", "group2"))
  for (b in bs) {
    expect_s3_class(b, "BootstrapResult")
    expect_equal(dim(b$adjacency_matrices), c(10L, 6L, 6L))
    p <- get_probabilities(b)
    expect_equal(dim(p), c(6L, 6L))
  }
})

test_that("lingam_multi_group_bootstrap is reproducible with the same seed", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)

  bs1 <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 10L, reg_method = "ols", seed = 99L, verbose = FALSE
  )
  bs2 <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 10L, reg_method = "ols", seed = 99L, verbose = FALSE
  )

  expect_equal(bs1$group1$adjacency_matrices, bs2$group1$adjacency_matrices)
  expect_equal(bs1$group2$adjacency_matrices, bs2$group2$adjacency_matrices)
})

test_that("lingam_multi_group_bootstrap compute_total_effects = FALSE skips total effects", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)

  bs <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 5L, reg_method = "ols", seed = 1L, verbose = FALSE,
    compute_total_effects = FALSE
  )

  expect_null(bs$group1$total_effects)
  expect_error(get_total_causal_effects(bs$group1), "compute_total_effects = FALSE")
})

test_that("parallel multi-group bootstrap is reproducible with same seed and n_cores", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2L, "requires >= 2 cores")

  mg <- generate_multi_group_sample(n = c(300, 300), seed = 1)

  bs1 <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 8L, reg_method = "ols", seed = 77L, verbose = FALSE,
    parallel = TRUE, n_cores = 2L
  )
  bs2 <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 8L, reg_method = "ols", seed = 77L, verbose = FALSE,
    parallel = TRUE, n_cores = 2L
  )

  expect_equal(bs1$group1$adjacency_matrices, bs2$group1$adjacency_matrices)
  expect_equal(bs1$group2$adjacency_matrices, bs2$group2$adjacency_matrices)
})

test_that("print.MultiGroupBootstrapResult reports per-group summaries", {
  mg <- generate_multi_group_sample(n = c(200, 200), seed = 1)
  bs <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 5L, reg_method = "ols", seed = 1L, verbose = FALSE
  )

  expect_output(print(bs), "MultiGroupBootstrapResult")
})
