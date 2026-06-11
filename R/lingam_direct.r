# =============================================================================
# Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
#
# License: MIT + file LICENSE
#
# Original work:
#   Copyright (c) 2019 T.Ikeuchi, G.Haraoka, M.Ide, W.Kurebayashi, S.Shimizu
#
# Portions of this work:
#   Copyright (c) 2026 O.Morimoto
#
# 内部実装は以下のファイルに分割されている：
#   R/search_causal_order.r : 因果順序の探索（pwling / kernel）と事前知識の処理
#   R/fit_regression.r      : 隣接行列の推定と回帰バックエンド
#   R/paths.r               : DAG のパス列挙と経路効果
# =============================================================================


#' Direct LiNGAM
#'
#' @param X 数値行列 (n_samples x n_features), data frame or matrix
#' @param prior_knowledge 事前知識行列 (n_features x n_features) または NULL。
#'   0: x_i から x_j への有向パスなし
#'   1: x_i から x_j への有向パスあり
#'  -1: 不明
#' @param apply_prior_knowledge_softly 事前知識をソフトに適用するか (logical)
#' @param measure 独立性の評価尺度 ("pwling" または "kernel")
#' @param reg_method 隣接行列推定の回帰手法。
#' "ols": 最小二乗法、
#' "lasso": LASSO回帰、
#' "adaptive_lasso": 適応的LASSO回帰（デフォルト）。
#' @param init_method 適応的LASSO回帰の初期重みの推定手法。
#' "ols": 最小二乗法（デフォルト）、
#' "ridge": Ridge回帰。
#' 多重共線性が疑われる場合はRidge回帰がおすすめ。
#' @param lambda LASSO のペナルティ（ラムダ）選択。
#' "lambda.min" : CV予測誤差最小, 予測精度優先。
#' "lambda.1se" : CV 1SEルール、ロバストで過学習しにくい。
#' "AIC": AIC最小。高速。
#' "BIC": BIC最小。高速、最もスパース。デフォルト。
#' "oracle" ：適応的LASSO回帰のみ。オラクル性を担保したλを選択。高速。
#' @return `LingamResult` オブジェクト（リスト）。以下の要素を含む：
#' * `adjacency_matrix`: 隣接行列 B (n_features x n_features)。
#'   **規則: `B[i, j]` は変数 j から変数 i への因果係数（j → i）。**
#'   ゼロ要素は因果関係なしを意味する。
#' * `causal_order`: 推定された因果順序（1-based インデックスの整数ベクトル）。
#'   先頭ほど上流（外生変数に近い）。
#' @importFrom stats sd lm.fit cov median quantile
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' # OLS (no additional packages required)
#' result <- lingam_direct(LiNGAM_sample_1000$data, reg_method = "ols")
#' round(result$adjacency_matrix, 3)
#'
#' \donttest{
#' # LASSO (requires glmnet)
#' result_lasso <- lingam_direct(LiNGAM_sample_1000$data)
#' round(result_lasso$adjacency_matrix, 3)
#' }
lingam_direct <- function(X,
                          prior_knowledge = NULL,
                          apply_prior_knowledge_softly = FALSE,
                          measure = "pwling",
                          reg_method = "adaptive_lasso",
                          lambda = "BIC",
                          init_method = "ols") {
  col_names <- if (is.data.frame(X)) names(X) else colnames(X)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix or data frame.", call. = FALSE)
  if (ncol(X) < 2) stop("X must have at least 2 variables (columns).", call. = FALSE)
  if (nrow(X) < 2) stop("X must have at least 2 observations (rows).", call. = FALSE)
  if (!is.null(col_names)) colnames(X) <- col_names

  measure <- match.arg(measure, c("pwling", "kernel"))
  reg_method <- match.arg(reg_method, c("adaptive_lasso", "lasso", "ols"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  if (!is.logical(apply_prior_knowledge_softly) || length(apply_prior_knowledge_softly) != 1) {
    stop("apply_prior_knowledge_softly must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  n_samples <- nrow(X)
  n_features <- ncol(X)
  # --- 事前知識の前処理 ---
  Aknw <- NULL
  partial_orders <- NULL
  if (!is.null(prior_knowledge)) {
    Aknw <- as.matrix(prior_knowledge)
    if (!all(dim(Aknw) == c(n_features, n_features))) {
      stop("The shape of prior knowledge must be (n_features, n_features)")
    }
    Aknw[Aknw < 0] <- NA
    if (!apply_prior_knowledge_softly) {
      partial_orders <- extract_partial_orders(Aknw)
    }
  }
  U <- seq_len(n_features)
  K <- integer(0)
  X_ <- X
  if (measure == "kernel") {
    X_ <- apply(X_, 2, function(col) {
      pop_sd <- sqrt(mean((col - mean(col))^2))
      (col - mean(col)) / pop_sd
    })
  }
  # --- 因果順序の探索 ---
  for (iter in seq_len(n_features)) {
    cand <- search_candidate(U, Aknw, apply_prior_knowledge_softly, partial_orders)
    m <- if (measure == "kernel") {
      search_causal_order_kernel(X_, U, cand$Uc, cand$Vj)
    } else {
      search_causal_order_pwling(X_, U, cand$Uc, cand$Vj)
    }
    for (i in U) {
      if (i != m) X_[, i] <- residual_vec(X_[, i], X_[, m])
    }
    K <- c(K, m)
    U <- setdiff(U, m)
    if (!is.null(Aknw) && !apply_prior_knowledge_softly && !is.null(partial_orders)) {
      if (nrow(partial_orders) > 0) {
        partial_orders <- partial_orders[partial_orders[, 1] != m, , drop = FALSE]
      }
    }
  }
  # --- 隣接行列の推定（回帰手法を選択可能）---
  B <- estimate_adjacency_matrix(X, K, Aknw,
    method = reg_method,
    lambda = lambda,
    init_method = init_method
  )
  colnames(B) <- rownames(B) <- colnames(X)
  result <- list(adjacency_matrix = B, causal_order = K)
  class(result) <- "LingamResult"
  result
}


#' LingamResult の print メソッド
#'
#' @param x LingamResult オブジェクト
#' @param digits 表示桁数
#' @param ... 追加引数（未使用）
#' @export
print.LingamResult <- function(x, digits = 3, ...) {
  n <- length(x$causal_order)
  var_names <- colnames(x$adjacency_matrix)
  order_labels <- if (!is.null(var_names)) {
    var_names[x$causal_order]
  } else {
    paste0("x", x$causal_order - 1L)
  }
  cat("Direct LiNGAM Result\n")
  cat(sprintf("  Variables : %d\n", n))
  cat(sprintf("  Causal order: %s\n", paste(order_labels, collapse = " -> ")))
  cat("\nAdjacency matrix (row = to, col = from):\n")
  print(round(x$adjacency_matrix, digits = digits))
  invisible(x)
}


# =============================================================================
# 共有内部ユーティリティ
# =============================================================================

#' lingam_direct() の返り値を検証する
#' @keywords internal
validate_lingam_result <- function(x) {
  if (!inherits(x, "LingamResult")) {
    stop("lingam_result must be the return value of lingam_direct().", call. = FALSE)
  }
}

#' Population standard deviation (n割り)
#' @keywords internal
sd_pop <- function(x) {
  sqrt(mean((x - mean(x))^2))
}
