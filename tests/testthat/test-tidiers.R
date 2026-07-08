test_that("tidy.LingamResult returns long data.frame of edges", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")
  td  <- tidy(res)

  expect_s3_class(td, "data.frame")
  expect_named(td, c("from", "to", "estimate"))
  expect_type(td$from, "character")
  expect_type(td$to, "character")
  expect_type(td$estimate, "double")
  # number of rows == number of non-zero edges
  expect_equal(nrow(td), sum(abs(res$adjacency_matrix) > 0))
})

test_that("tidy.LingamResult threshold filters edges", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  td_all  <- tidy(res, threshold = 0)
  td_high <- tidy(res, threshold = 100)  # exclude all edges

  expect_gte(nrow(td_all), nrow(td_high))
  expect_equal(nrow(td_high), 0L)
})

test_that("tidy.LingamResult from/to follow j -> i convention", {
  # pass the true structure directly to confirm the convention: B["x0","x3"] = 3 is x3 -> x0
  dat  <- generate_lingam_sample_6(n = 100, seed = 1)
  fake <- structure(
    list(adjacency_matrix = dat$true_adjacency, causal_order = 1:6),
    class = "LingamResult"
  )
  td <- tidy(fake)

  edge <- td[td$from == "x3" & td$to == "x0", ]
  expect_equal(nrow(edge), 1L)
  expect_equal(edge$estimate, 3.0)
})

test_that("tidy.LingamResult returns 0-row data.frame for empty graph", {
  fake <- structure(
    list(adjacency_matrix = matrix(0, 3, 3,
           dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
         causal_order = 1:3),
    class = "LingamResult"
  )
  td <- tidy(fake)

  expect_equal(nrow(td), 0L)
  expect_named(td, c("from", "to", "estimate"))
})

test_that("glance.LingamResult returns one-row summary", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")
  g   <- glance(res)

  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_true(all(c("n_variables", "n_edges", "causal_order") %in% names(g)))
  expect_equal(g$n_variables, 6L)
})

test_that("tidy.BootstrapResult returns direction counts", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  bs  <- lingam_direct_bootstrap(dat$data, n_sampling = 15L, reg_method = "ols", seed = 42L)
  td  <- tidy(bs)

  expect_s3_class(td, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion") %in% names(td)))
})

test_that("tidy/glance for ParceLingamResult keep NA entries", {
  B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  B["b", "a"] <- 1.5
  B["a", "c"] <- NA
  B["c", "a"] <- NA
  fake <- structure(
    list(adjacency_matrix = B, causal_order = list(2L, c(1L, 3L)),
         p_values = matrix(0, 3, 3), independence = "hsic"),
    class = "ParceLingamResult"
  )

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(sum(is.na(td$estimate)), 2L)
  expect_equal(td$estimate[td$from == "a" & td$to == "b"], 1.5)

  g <- glance(fake)
  expect_equal(nrow(g), 1L)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_na_entries, 2L)
  expect_match(g$causal_order, "(a, c)", fixed = TRUE)
})

test_that("tidy/glance for RCDResult count confounded pairs", {
  B <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  B["c", "b"] <- -0.8
  B["a", "b"] <- NA
  B["b", "a"] <- NA
  fake <- structure(
    list(adjacency_matrix = B,
         ancestors_list = list(integer(0), integer(0), 2L)),
    class = "RCDResult"
  )

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(sum(is.na(td$estimate)), 2L)

  g <- glance(fake)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_confounded_pairs, 1L)
  expect_false("causal_order" %in% names(g))
})

test_that("tidy/glance for LiMResult work", {
  B <- matrix(0, 3, 3, dimnames = list(paste0("x", 1:3), paste0("x", 1:3)))
  B["x2", "x1"] <- 1.2
  fake <- structure(
    list(adjacency_matrix = B, causal_order = 1:3,
         is_continuous = c(TRUE, FALSE, TRUE)),
    class = "LiMResult"
  )

  td <- tidy(fake)
  expect_named(td, c("from", "to", "estimate"))
  expect_equal(nrow(td), 1L)

  g <- glance(fake)
  expect_equal(g$n_edges, 1L)
  expect_equal(g$n_discrete, 1L)
  expect_equal(g$causal_order, "x1 -> x2 -> x3")
})

test_that("tidy/glance for MultiGroupLingamResult stack groups", {
  mg  <- generate_multi_group_sample(n = c(300, 300), seed = 42)
  res <- lingam_multi_group(mg$data_list, reg_method = "ols")

  td <- tidy(res)
  expect_named(td, c("group", "from", "to", "estimate"))
  expect_setequal(unique(td$group), names(res$adjacency_matrices))

  # each group's block equals tidy() of that group's LingamResult
  g1_name <- names(res$adjacency_matrices)[1]
  block <- td[td$group == g1_name, c("from", "to", "estimate")]
  expect_equal(block, tidy(get_group_result(res, 1)), ignore_attr = TRUE)

  gl <- glance(res)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_groups, length(res$adjacency_matrices))
  expect_true(all(c("n_variables", "causal_order") %in% names(gl)))
})

test_that("tidy.MultiGroupBootstrapResult adds a group column", {
  mg <- generate_multi_group_sample(n = c(300, 300), seed = 42)
  bs <- lingam_multi_group_bootstrap(mg$data_list,
    n_sampling = 10L, reg_method = "ols", seed = 42, verbose = FALSE
  )

  td <- tidy(bs)
  expect_true(all(c("group", "from", "to", "count", "proportion") %in% names(td)))
  expect_setequal(unique(td$group), names(bs))
})

test_that("tidy.ImputationBootstrapResult collapses imputations", {
  skip_if_not_installed("mice")

  set.seed(1)
  dat <- generate_lingam_sample_6(n = 200, seed = 1)$data
  dat$x5[sample.int(nrow(dat), 20)] <- NA
  bs <- bootstrap_with_imputation(dat,
    n_sampling = 5L, n_repeats = 2L, seed = 42, verbose = FALSE
  )

  td <- tidy(bs)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("from", "to", "count", "proportion") %in% names(td)))
})
