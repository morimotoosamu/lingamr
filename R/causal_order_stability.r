# =============================================================================
# Causal order stability from bootstrap
# =============================================================================

#' ブートストラップによる因果順序の安定性を評価
#'
#' 各ブートストラップ標本で推定された因果順序 (causal_order) を集計し、順序が
#' どれだけ安定しているかを数値化する。各変数の順位分布、変数ペアの先行確率、
#' および全体の安定性スコアを返す。
#'
#' @param result BootstrapResult オブジェクト（現行バージョンで実行したもの）
#' @param labels 変数名ベクトル (NULL の場合は x0, x1, ... を自動生成)
#' @return `causal_order_stability` クラスのリスト。以下を含む：
#' * `rank_summary`: 各変数の順位の要約 (variable, mean_rank, sd_rank,
#'   median_rank, mode_rank)。mean_rank 昇順（上流から）にソート済み。
#'   順位は 1 が最も上流。
#' * `precedence_matrix`: 先行確率行列。`P[i, j]` は変数 i が変数 j より上流
#'   （先）に位置したブートストラップ標本の割合。
#' * `stability_score`: 全体の安定性スコア。0（順序がランダム）〜
#'   1（全標本で順序が一致）。各変数ペアの先行確率が 0/1 に近いほど高い。
#' * `n_sampling`: ブートストラップ標本数。
#' @importFrom stats sd median
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#' get_causal_order_stability(bs, labels = names(dat$data))
get_causal_order_stability <- function(result, labels = NULL) {
  stopifnot(inherits(result, "BootstrapResult"))
  if (is.null(result$causal_orders)) {
    stop(
      "This BootstrapResult does not contain causal orders. ",
      "Please re-run lingam_direct_bootstrap() with the current version of the package.",
      call. = FALSE
    )
  }

  co <- result$causal_orders
  n_sampling <- nrow(co)
  n_features <- ncol(co)

  if (is.null(labels)) {
    labels <- paste0("x", seq_len(n_features) - 1)
  } else if (length(labels) != n_features) {
    stop(sprintf(
      "Length of 'labels' (%d) must equal the number of variables (%d).",
      length(labels), n_features
    ), call. = FALSE)
  }

  # 各標本での各変数の順位（1 = 最上流）
  rank_matrix <- matrix(0L, nrow = n_sampling, ncol = n_features)
  for (s in seq_len(n_sampling)) {
    rank_matrix[s, co[s, ]] <- seq_len(n_features)
  }

  mode_rank <- function(v) {
    tb <- table(v)
    as.integer(names(tb)[which.max(tb)])
  }

  rank_summary <- data.frame(
    variable    = labels,
    mean_rank   = colMeans(rank_matrix),
    sd_rank     = apply(rank_matrix, 2, stats::sd),
    median_rank = apply(rank_matrix, 2, stats::median),
    mode_rank   = apply(rank_matrix, 2, mode_rank),
    stringsAsFactors = FALSE
  )
  rank_summary <- rank_summary[order(rank_summary$mean_rank), ]
  rownames(rank_summary) <- NULL

  # 先行確率行列 P[i, j] = P(i precedes j)
  P <- matrix(0, n_features, n_features, dimnames = list(labels, labels))
  for (i in seq_len(n_features)) {
    for (j in seq_len(n_features)) {
      if (i == j) next
      P[i, j] <- mean(rank_matrix[, i] < rank_matrix[, j])
    }
  }

  # 全体安定性スコア: 各ペアの先行確率が 0/1 にどれだけ近いか
  ps <- P[upper.tri(P)]
  stability_score <- if (length(ps) == 0) NA_real_ else mean(2 * abs(ps - 0.5))

  out <- list(
    rank_summary      = rank_summary,
    precedence_matrix = P,
    stability_score   = stability_score,
    n_sampling        = n_sampling
  )
  class(out) <- "causal_order_stability"
  out
}


#' causal_order_stability の print メソッド
#'
#' @param x `causal_order_stability` オブジェクト
#' @param ... 追加引数（未使用）
#' @export
print.causal_order_stability <- function(x, ...) {
  cat("=== Causal Order Stability ===\n")
  cat(sprintf("Bootstrap samples:       %d\n", x$n_sampling))
  cat(sprintf("Overall stability score: %.3f  (0 = random, 1 = fully stable)\n\n",
              x$stability_score))

  cat("Rank summary (sorted by mean rank; 1 = most upstream):\n")
  rs <- x$rank_summary
  disp <- data.frame(
    variable    = rs$variable,
    mean_rank   = sprintf("%.2f", rs$mean_rank),
    sd_rank     = sprintf("%.2f", rs$sd_rank),
    median_rank = rs$median_rank,
    mode_rank   = rs$mode_rank,
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)

  cat("\nPrecedence probability P[i, j] = P(variable i precedes j):\n")
  print(round(x$precedence_matrix, 2))
  invisible(x)
}
