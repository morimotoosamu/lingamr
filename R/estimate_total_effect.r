#' 指定した2変数間の総合因果効果を推定
#'
#' @param X 元データ (matrix or data.frame)
#' @param lingam_result lingam_direct() の返り値
#' @param from_index 原因変数 (1-based index or 変数名)
#' @param to_index 結果変数 (1-based index or 変数名)
#' @param method 回帰手法 ("ols", "lasso", "adaptive_lasso")デフォルトはadaptive_lasso
#' @param lambda ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")デフォルトはBIC
#' @param init_method 適応的LASSO回帰の初期重みの推定手法 ("ols" または "ridge")
#' @return 推定された総合因果効果
#' @importFrom stats cov
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' model <- LiNGAM_sample_1000$data |>
#'   lingam_direct()
#'
#' LiNGAM_sample_1000$data |>
#'   estimate_total_effect(model, 4, 1)
estimate_total_effect <- function(X, lingam_result, from_index, to_index,
                                  method = "adaptive_lasso", lambda = "BIC",
                                  init_method = "ols") {
  validate_lingam_result(lingam_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  if (is.data.frame(X)) {
    col_names <- colnames(X)
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
    col_names <- colnames(X)
  }
  if (!is.numeric(X)) stop("X must be a numeric matrix or data.frame.", call. = FALSE)

  n_features <- ncol(X)
  if (ncol(lingam_result$adjacency_matrix) != n_features) {
    stop(
      "X has ", n_features, " variables but lingam_result was estimated from ",
      ncol(lingam_result$adjacency_matrix), " variables.",
      call. = FALSE
    )
  }

  # --- 変数名 → インデックス変換 ---
  resolve_index <- function(idx, arg_name) {
    if (is.character(idx)) {
      if (is.null(col_names)) {
        stop(sprintf("'%s' was specified as a name, but X has no column names.", arg_name))
      }
      pos <- match(idx, col_names)
      if (is.na(pos)) {
        stop(sprintf(
          "Variable '%s' not found. Available: %s",
          idx, paste(col_names, collapse = ", ")
        ))
      }
      return(pos)
    } else if (is.numeric(idx)) {
      idx <- as.integer(idx)
      if (idx < 1 || idx > n_features) {
        stop(sprintf("'%s' must be between 1 and %d.", arg_name, n_features))
      }
      return(idx)
    }
    stop(sprintf("'%s' must be integer or character.", arg_name))
  }

  from_index <- resolve_index(from_index, "from_index")
  to_index <- resolve_index(to_index, "to_index")
  if (from_index == to_index) stop("from_index and to_index must differ.")

  adjacency_matrix <- lingam_result$adjacency_matrix
  causal_order <- lingam_result$causal_order

  # --- 因果順序チェック ---
  from_order <- which(causal_order == from_index)
  to_order <- which(causal_order == to_index)

  from_label <- if (!is.null(col_names)) col_names[from_index] else paste0("x", from_index)
  to_label <- if (!is.null(col_names)) col_names[to_index] else paste0("x", to_index)

  if (from_order > to_order) {
    warning(sprintf(
      "Causal order of %s (to) is earlier than %s (from). Result may be incorrect.",
      to_label, from_label
    ))
  }

  # --- 親変数の特定と回帰 ---
  parents <- which(abs(adjacency_matrix[from_index, ]) > 0)
  predictors <- unique(c(from_index, parents))
  from_pos <- which(predictors == from_index)

  y <- X[, to_index]
  Xp <- X[, predictors, drop = FALSE]

  coefs <- switch(method,
    "ols"            = fit_ols(y, Xp),
    "lasso"          = fit_lasso(y, Xp, lambda),
    "adaptive_lasso" = fit_adaptive_lasso(y, Xp, lambda, init_method = init_method)
  )

  return(coefs[from_pos])
}


#' 全変数間の総合因果効果を一括推定
#'
#' @param X 元データ (n_samples x n_features)
#' @param lingam_result lingam_direct() の返り値
#' @param method 回帰手法 ("ols", "lasso", "adaptive_lasso")
#' @param lambda ラムダ選択 ("lambda.min", "lambda.1se", "AIC", "BIC")
#' @param init_method 適応的LASSO回帰の初期重みの推定手法 ("ols" または "ridge")
#' @return 総合因果効果の行列 (n_features x n_features)。
#'   **規則: `TE[i, j]` は変数 j から変数 i への総合因果効果（j → i）。**
#'   隣接行列 `adjacency_matrix` と同じ添字規則。直接効果と間接効果の合計。
#' @importFrom stats cov
#' @export
#' @examples
#' LiNGAM_sample_1000 <- generate_lingam_sample_6()
#'
#' model <- LiNGAM_sample_1000$data |>
#'   lingam_direct()
#'
#' LiNGAM_sample_1000$data |>
#'   estimate_all_total_effects(model)
estimate_all_total_effects <- function(X,
                                       lingam_result,
                                       method = "adaptive_lasso",
                                       lambda = "BIC",
                                       init_method = "ols") {
  validate_lingam_result(lingam_result)
  method <- match.arg(method, c("adaptive_lasso", "lasso", "ols"))
  lambda <- match.arg(lambda, c("BIC", "AIC", "lambda.min", "lambda.1se", "oracle"))
  init_method <- match.arg(init_method, c("ols", "ridge"))

  X <- as.matrix(X)
  n_features <- ncol(X)
  if (ncol(lingam_result$adjacency_matrix) != n_features) {
    stop(
      "X has ", n_features, " variables but lingam_result was estimated from ",
      ncol(lingam_result$adjacency_matrix), " variables.",
      call. = FALSE
    )
  }

  causal_order <- lingam_result$causal_order
  adj_matrix <- lingam_result$adjacency_matrix

  TE <- matrix(0, nrow = n_features, ncol = n_features)
  if (!is.null(colnames(X))) {
    rownames(TE) <- colnames(X)
    colnames(TE) <- colnames(X)
  }

  # 共分散行列はループ不変なので1回だけ計算する（OLS 経路でのみ使用）
  cov_mat <- if (method == "ols") cov(X) else NULL

  for (i in 1:(n_features - 1)) {
    from_idx <- causal_order[i]

    parents <- which(abs(adj_matrix[from_idx, ]) > 0)
    predictors <- unique(c(from_idx, parents))
    from_pos <- which(predictors == from_idx)

    downstream <- causal_order[(i + 1):n_features]

    if (method == "ols") {
      # --- OLS: 共分散行列ベースで一括計算（最速）---
      cov_xx <- cov_mat[predictors, predictors, drop = FALSE]
      cov_xy <- cov_mat[predictors, downstream, drop = FALSE]
      beta_mat <- solve(cov_xx, cov_xy)
      TE[downstream, from_idx] <- beta_mat[from_pos, ]
    } else {
      # --- LASSO / Adaptive LASSO ---
      Xp <- X[, predictors, drop = FALSE]
      for (to_idx in downstream) {
        y <- X[, to_idx]
        coefs <- switch(method,
          "lasso"          = fit_lasso(y, Xp, lambda),
          "adaptive_lasso" = fit_adaptive_lasso(y, Xp, lambda,
                                                init_method = init_method)
        )
        TE[to_idx, from_idx] <- coefs[from_pos]
      }
    }
  }

  return(TE)
}
