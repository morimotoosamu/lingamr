test_that("make_prior_knowledge returns correct dimensions and defaults", {
  pk <- make_prior_knowledge(4)

  expect_true(is.matrix(pk))
  expect_equal(dim(pk), c(4L, 4L))
  # 対角は -1
  expect_true(all(diag(pk) == -1L))
  # 非対角はすべて -1（制約なし）
  expect_true(all(pk[row(pk) != col(pk)] == -1L))
})

test_that("exogenous_variables sets row to 0", {
  pk <- make_prior_knowledge(4, exogenous_variables = 2)

  expect_true(all(pk[2, -2] == 0L))   # 行2（自分以外）はゼロ
  expect_equal(pk[2, 2], -1L)          # 対角は -1 のまま
})

test_that("sink_variables sets column to 0", {
  pk <- make_prior_knowledge(4, sink_variables = 3)

  expect_true(all(pk[-3, 3] == 0L))   # 列3（自分以外）はゼロ
  expect_equal(pk[3, 3], -1L)          # 対角は -1 のまま
})

test_that("paths sets specified edges to 1", {
  # (from=1, to=3) → pk[3, 1] = 1
  pk <- make_prior_knowledge(4, paths = list(c(1, 3)))

  expect_equal(pk[3, 1], 1L)
})

test_that("no_paths sets specified edges to 0", {
  # (from=2, to=4) → pk[4, 2] = 0
  pk <- make_prior_knowledge(4, no_paths = list(c(2, 4)))

  expect_equal(pk[4, 2], 0L)
})

test_that("labels are attached and name-based lookup works", {
  labs <- c("a", "b", "c", "d")
  pk <- make_prior_knowledge(4,
    exogenous_variables = "a",
    paths = list(c("a", "b")),
    labels = labs
  )

  expect_equal(rownames(pk), labs)
  expect_equal(colnames(pk), labs)
  expect_true(all(pk["a", c("b", "c", "d")] == 0L))  # exogenous 行
  expect_equal(pk["b", "a"], 1L)                       # paths
})

test_that("make_prior_knowledge errors on invalid inputs", {
  expect_error(make_prior_knowledge(1))          # n_variables < 2
  expect_error(make_prior_knowledge(4, exogenous_variables = 5))  # index out of range
  expect_error(make_prior_knowledge(4,
    exogenous_variables = "x",
    labels = NULL))                              # 名前指定なのに labels なし
  expect_error(make_prior_knowledge(4,
    labels = c("a", "b", "a", "c")))            # 重複ラベル
})
