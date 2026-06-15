test_that("lingam_direct returns LingamResult with correct structure", {
  dat <- generate_lingam_sample_6(n = 500, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
  expect_named(res, c("adjacency_matrix", "causal_order"))
  expect_true(is.matrix(res$adjacency_matrix))
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  # 列名が保持されている
  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct identifies x3 as root (first in causal order)", {
  # x3 は親を持たない外生変数なので因果順序の先頭になるはず
  dat <- generate_lingam_sample_6(n = 2000, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_equal(res$causal_order[1], which(names(dat$data) == "x3"))
})

test_that("lingam_direct accepts data.frame and preserves colnames", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_equal(colnames(res$adjacency_matrix), names(dat$data))
  expect_equal(rownames(res$adjacency_matrix), names(dat$data))
})

test_that("lingam_direct accepts matrix input", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  mat <- as.matrix(dat$data)
  res <- lingam_direct(mat, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})

test_that("lingam_direct errors on invalid inputs", {
  dat <- generate_lingam_sample_6(n = 100, seed = 1)

  expect_error(lingam_direct(dat$data, measure = "bad_measure"))
  expect_error(lingam_direct(dat$data, reg_method = "bad_method"))
  expect_error(lingam_direct(dat$data, lambda = "bad_lambda"))
  expect_error(lingam_direct(matrix(1:4, nrow = 2)))   # 非 numeric（整数はOKなので別ケース）
  expect_error(lingam_direct(as.data.frame(matrix(1, nrow = 5, ncol = 1))))  # 1変数
  expect_error(lingam_direct(matrix(numeric(0), nrow = 0, ncol = 3)))        # 0行
})

test_that("print.LingamResult runs without error", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  res <- lingam_direct(dat$data, reg_method = "ols")

  expect_output(print(res), "Direct LiNGAM Result")
  expect_output(print(res), "Causal order")
  expect_output(print(res), "Adjacency matrix")
})

test_that("lingam_direct with prior_knowledge runs without error", {
  dat <- generate_lingam_sample_6(n = 300, seed = 1)
  pk <- make_prior_knowledge(6,
    exogenous_variables = 4,
    labels = names(dat$data)
  )
  res <- lingam_direct(dat$data, prior_knowledge = pk, reg_method = "ols")

  expect_s3_class(res, "LingamResult")
})

test_that("soft prior knowledge with -1 (unknown) entries runs without error", {
  # 回帰テスト: -1 は前処理で NA になり、search_candidate() 内の
  # if (sum(...) == 0) が NA を返して落ちていた
  dat <- generate_lingam_sample_6(n = 200, seed = 1)
  pk <- matrix(-1L, 6, 6)
  pk[2, 1] <- 1L

  res <- lingam_direct(dat$data,
    prior_knowledge = pk,
    apply_prior_knowledge_softly = TRUE
  )

  expect_s3_class(res, "LingamResult")
  expect_setequal(res$causal_order, 1:6)
})

test_that("search_candidate handles NA entries like Python's NaN semantics", {
  # Python 版では NaN を含む行の sum() == 0 は False（エラーではない）
  U <- 1:3

  # 行1: すべて 0（既知の外生変数）→ Uc に入る
  # 行2: NA を含む → 外生候補にならない
  # 行3: 1 を含む → 外生候補にならない
  Aknw <- matrix(NA_real_, 3, 3)
  Aknw[1, 2:3] <- 0
  Aknw[2, 3] <- 0
  Aknw[3, 2] <- 1

  cand <- search_candidate(U, Aknw,
    apply_prior_knowledge_softly = TRUE,
    partial_orders = NULL
  )
  expect_equal(cand$Uc, 1L)

  # すべて NA（不明だらけ）でも落ちず、Uc は U 全体になる
  Aknw_all_na <- matrix(NA_real_, 3, 3)
  cand_na <- search_candidate(U, Aknw_all_na,
    apply_prior_knowledge_softly = TRUE,
    partial_orders = NULL
  )
  expect_equal(cand_na$Uc, U)
  expect_equal(cand_na$Vj, integer(0))
})

test_that("measure = 'kernel' returns a valid LingamResult", {
  dat <- generate_lingam_sample_6(n = 100, seed = 42)
  res <- lingam_direct(dat$data, measure = "kernel", reg_method = "ols")

  expect_s3_class(res, "LingamResult")
  expect_equal(sort(res$causal_order), 1:6)
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
})

test_that("mutual_information_kernel returns finite non-negative values", {
  set.seed(1)
  x <- rnorm(100)
  y <- 0.8 * x + runif(100)
  z <- rnorm(100)
  param <- c(2e-2, 1.0)

  mi_dep <- mutual_information_kernel(x, y, param)
  mi_ind <- mutual_information_kernel(x, z, param)

  expect_true(is.finite(mi_dep) && is.finite(mi_ind))
  # 依存ペアの MI は独立ペアより大きい
  expect_gt(mi_dep, mi_ind)
})

test_that("reg_method = 'ols' does not require glmnet", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)

  # glmnet が無い環境を模擬する
  local_mocked_bindings(
    check_glmnet_available = function(method) {
      stop(sprintf(
        "Package 'glmnet' is required for reg_method = \"%s\". Please install it.",
        method
      ), call. = FALSE)
    }
  )

  expect_no_error(lingam_direct(dat$data, reg_method = "ols"))
  expect_error(lingam_direct(dat$data, reg_method = "lasso"), "glmnet")
  expect_error(lingam_direct(dat$data, reg_method = "adaptive_lasso"), "glmnet")
  expect_error(lingam_direct(dat$data, reg_method = "ridge"), "glmnet")
})

test_that("reg_method = 'ridge' returns a valid LingamResult", {
  dat <- generate_lingam_sample_6(n = 300, seed = 42)
  res <- lingam_direct(dat$data, reg_method = "ridge")

  expect_s3_class(res, "LingamResult")
  expect_equal(dim(res$adjacency_matrix), c(6L, 6L))
  expect_equal(length(res$causal_order), 6L)
  # Ridge はスパースでないため非ゼロ係数が多い
  expect_gt(sum(res$adjacency_matrix != 0), 0L)
})

test_that("reg_method = 'ridge' with lambda = 'oracle' errors", {
  dat <- generate_lingam_sample_6(n = 200, seed = 1)

  expect_error(
    lingam_direct(dat$data, reg_method = "ridge", lambda = "oracle"),
    "oracle"
  )
})
