# =============================================================================
# broom (generics) tidiers for lingamr
# =============================================================================

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance


#' LingamResult を tidy な data.frame に変換
#'
#' 推定された隣接行列を、1 行が 1 エッジの long 形式 data.frame に変換する。
#' `B[i, j]`（j → i の係数）の規則に従い、`from` 列が原因、`to` 列が結果となる。
#' ggplot2 や ggraph での可視化、dplyr でのフィルタリングに便利。
#'
#' @param x [lingam_direct()] の返り値（`LingamResult` オブジェクト）
#' @param threshold この絶対値以下の係数はエッジとみなさない (default: 0)
#' @param ... 未使用
#' @return data.frame(from, to, estimate)。`from`/`to` は変数名（文字列）、
#'   `estimate` は因果係数。エッジが無ければ 0 行の data.frame。
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' model <- lingam_direct(dat$data, reg_method = "ols")
#' tidy(model)
tidy.LingamResult <- function(x, threshold = 0, ...) {
  B <- x$adjacency_matrix
  p <- ncol(B)
  var_names <- colnames(B)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(p) - 1)

  idx <- which(abs(B) > threshold, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0), to = character(0),
      estimate = numeric(0), stringsAsFactors = FALSE
    ))
  }

  # B[i, j] は j -> i。行 i が to、列 j が from。
  ord <- order(idx[, 2], idx[, 1])
  idx <- idx[ord, , drop = FALSE]
  data.frame(
    from     = var_names[idx[, 2]],
    to       = var_names[idx[, 1]],
    estimate = B[idx],
    stringsAsFactors = FALSE
  )
}


#' LingamResult の1行サマリを取得
#'
#' モデル全体を1行に要約する。残差を計算しないためデータ `X` は不要。
#' 残差ベースの診断が必要な場合は [summary_lingam()] を使用すること。
#'
#' @param x [lingam_direct()] の返り値（`LingamResult` オブジェクト）
#' @param ... 未使用
#' @return 1 行の data.frame(n_variables, n_edges, causal_order)
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' model <- lingam_direct(dat$data, reg_method = "ols")
#' glance(model)
glance.LingamResult <- function(x, ...) {
  B <- x$adjacency_matrix
  p <- ncol(B)
  var_names <- colnames(B)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(p) - 1)
  data.frame(
    n_variables  = p,
    n_edges      = sum(abs(B) > 0),
    causal_order = paste(var_names[x$causal_order], collapse = " -> "),
    stringsAsFactors = FALSE
  )
}


#' BootstrapResult を tidy な data.frame に変換
#'
#' 各因果方向の出現回数・割合・効果量の要約を返す。内部で
#' [get_causal_direction_counts()] を呼び出すため、同関数の引数を `...` で渡せる。
#'
#' @param x [lingam_direct_bootstrap()] の返り値（`BootstrapResult` オブジェクト）
#' @param ... [get_causal_direction_counts()] に渡す引数
#'   (`n_directions`, `min_causal_effect`, `split_by_causal_effect_sign`, `labels` など)
#' @return data.frame (from, to, count, proportion, ...)
#' @export
#' @examples
#' dat <- generate_lingam_sample_6()
#' bs <- lingam_direct_bootstrap(dat$data, n_sampling = 30L, seed = 42)
#' tidy(bs)
tidy.BootstrapResult <- function(x, ...) {
  get_causal_direction_counts(x, ...)
}
