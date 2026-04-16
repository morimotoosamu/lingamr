# =============================================================================
# Direct LiNGAM - R Implementation
# Based on the Python implementation from the LiNGAM Project
# https://sites.google.com/view/sshimizu06/lingam
# https://github.com/cdt15/lingam
# =============================================================================


#' Direct LiNGAM
#'
#' @param X 数値行列 (n_samples x n_features)
#' @param prior_knowledge 事前知識行列 (n_features x n_features) または NULL
#'   0: x_i から x_j への有向パスなし
#'   1: x_i から x_j への有向パスあり
#'  -1: 不明
#' @param apply_prior_knowledge_softly 事前知識をソフトに適用するか (logical)
#' @param measure 独立性の評価尺度 ("pwling" または "kernel")
#' @return list(adjacency_matrix, causal_order)
#' @importFrom stats sd lm.fit cov median quantile
#' @export
direct_lingam <- function(X,
                          prior_knowledge = NULL,
                          apply_prior_knowledge_softly = FALSE,
                          measure = "pwling") {
  # --- 入力チェック ---
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix.")
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

  # --- 因果探索メインループ ---
  U <- seq_len(n_features) # 1-based index
  K <- integer(0)
  X_ <- X # コピー

  if (measure == "kernel") {
    X_ <- scale(X_) # 標準化
  }

  for (iter in seq_len(n_features)) {
    # 候補探索
    cand <- search_candidate(U, Aknw, apply_prior_knowledge_softly, partial_orders)
    Uc <- cand$Uc
    Vj <- cand$Vj

    # 因果順序の探索
    if (measure == "kernel") {
      m <- search_causal_order_kernel(X_, U, Uc, Vj)
    } else {
      m <- search_causal_order_pwling(X_, U, Uc, Vj)
    }

    # 残差化
    for (i in U) {
      if (i != m) {
        X_[, i] <- residual_vec(X_[, i], X_[, m])
      }
    }

    K <- c(K, m)
    U <- setdiff(U, m)

    # 部分順序の更新
    if (!is.null(Aknw) && !apply_prior_knowledge_softly && !is.null(partial_orders)) {
      if (nrow(partial_orders) > 0) {
        partial_orders <- partial_orders[partial_orders[, 1] != m, , drop = FALSE]
      }
    }
  }

  # --- 隣接行列の推定 ---
  B <- estimate_adjacency_matrix(X, K, Aknw)

  return(list(
    adjacency_matrix = B,
    causal_order     = K
  ))
}


# =============================================================================
# 内部関数群
# =============================================================================

#' 事前知識から部分順序を抽出
#' @param pk 事前知識行列 (NaN = 不明)
#' @return matrix (n x 2), 各行は [from, to] の部分順序
#' @keywords internal
extract_partial_orders <- function(pk) {
  # パスがあるペア (pk == 1)
  path_idx <- which(pk == 1, arr.ind = TRUE)
  # パスがないペア (pk == 0)
  no_path_idx <- which(pk == 0, arr.ind = TRUE)

  # --- パスありペアの矛盾チェック ---
  if (nrow(path_idx) > 0) {
    check_pairs <- rbind(path_idx, path_idx[, 2:1, drop = FALSE])
    dup <- duplicated(check_pairs) | duplicated(check_pairs, fromLast = TRUE)
    if (any(dup)) {
      bad <- unique(check_pairs[dup, , drop = FALSE])
      stop(paste(
        "The prior knowledge contains inconsistencies at indices:",
        paste(apply(bad, 1, function(r) paste0("(", r[1], ",", r[2], ")")),
          collapse = ", "
        )
      ))
    }
  }

  # --- パスなしペアの重複除去 ---
  if (nrow(no_path_idx) > 0) {
    check_pairs2 <- rbind(no_path_idx, no_path_idx[, 2:1, drop = FALSE])
    # 双方向に 0 が入っているペアを見つけて除外
    pair_key <- paste(check_pairs2[, 1], check_pairs2[, 2], sep = ",")
    tbl <- table(pair_key)
    dup_keys <- names(tbl[tbl > 1])
    # no_path_idx のうち、双方向に存在するものを除外
    no_path_key <- paste(no_path_idx[, 1], no_path_idx[, 2], sep = ",")
    keep <- !(no_path_key %in% dup_keys)
    no_path_idx <- no_path_idx[keep, , drop = FALSE]
  }

  # path_pairs と no_path_pairs[:, [1,0]] を結合
  combined <- matrix(nrow = 0, ncol = 2)
  if (nrow(path_idx) > 0) {
    combined <- rbind(combined, path_idx)
  }
  if (nrow(no_path_idx) > 0) {
    combined <- rbind(combined, no_path_idx[, 2:1, drop = FALSE])
  }

  if (nrow(combined) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  combined <- unique(combined)
  # [to, from] -> [from, to]
  result <- combined[, 2:1, drop = FALSE]
  colnames(result) <- NULL
  return(result)
}


#' 残差 (xi を xj に回帰したときの残差)
#' 残差ベクトルの計算
#' @param xi 対象変数ベクトル
#' @param xj 説明変数ベクトル
#' @return 回帰後の残差ベクトル
#' @keywords internal
residual_vec <- function(xi, xj) {
  n <- length(xi)
  xi_c <- xi - mean(xi)
  xj_c <- xj - mean(xj)

  cov_val <- as.numeric(crossprod(xi_c, xj_c)) / n
  var_val <- as.numeric(crossprod(xj_c)) / n

  return(xi - (cov_val / var_val) * xj)
}


#' エントロピーの最大エントロピー近似
#' @param u 入力ベクトル
#' @return 近似エントロピー値
#' @keywords internal
entropy_approx <- function(u) {
  k1 <- 79.047
  k2 <- 7.4129
  gamma <- 0.37457
  (1 + log(2 * pi)) / 2 -
    k1 * (mean(log(cosh(u))) - gamma)^2 -
    k2 * (mean(u * exp(-u^2 / 2)))^2
}


#' 相互情報量の差
#' @param xi_std 標準化された xi
#' @param xj_std 標準化された xj
#' @param ri_j xi を xj で回帰した残差
#' @param rj_i xj を xi で回帰した残差
#' @return 相互情報量の差
#' @keywords internal
diff_mutual_info <- function(xi_std, xj_std, ri_j, rj_i) {
  (entropy_approx(xj_std) + entropy_approx(ri_j / stats::sd(ri_j))) -
    (entropy_approx(xi_std) + entropy_approx(rj_i / stats::sd(rj_i)))
}


#' 候補変数の探索
#' @param U 現在の未確定変数の集合
#' @param Aknw 事前知識行列
#' @param apply_prior_knowledge_softly ソフト適用の有無
#' @param partial_orders 抽出された部分順序
#' @return list(Uc, Vj)
#' @keywords internal
search_candidate <- function(U, Aknw, apply_prior_knowledge_softly, partial_orders) {
  # 事前知識なし

  if (is.null(Aknw)) {
    return(list(Uc = U, Vj = integer(0)))
  }

  # --- ハード適用 ---
  if (!apply_prior_knowledge_softly) {
    if (!is.null(partial_orders) && nrow(partial_orders) > 0) {
      Uc <- setdiff(U, partial_orders[, 2])
      if (length(Uc) == 0) Uc <- U
      return(list(Uc = Uc, Vj = integer(0)))
    } else {
      return(list(Uc = U, Vj = integer(0)))
    }
  }

  # --- ソフト適用 ---
  # 外生変数の探索
  Uc <- integer(0)
  for (j in U) {
    index <- setdiff(U, j)
    if (sum(Aknw[j, index], na.rm = FALSE) == 0) {
      Uc <- c(Uc, j)
    }
  }

  # 内生変数の探索 → 候補の絞り込み
  if (length(Uc) == 0) {
    U_end <- integer(0)
    for (j in U) {
      index <- setdiff(U, j)
      s <- sum(Aknw[j, index], na.rm = TRUE)
      if (!is.na(s) && s > 0) {
        U_end <- c(U_end, j)
      }
    }
    # シンク特徴量
    for (i in U) {
      index <- setdiff(U, i)
      if (sum(Aknw[index, i], na.rm = FALSE) == 0) {
        U_end <- c(U_end, i)
      }
    }
    Uc <- setdiff(U, unique(U_end))
    if (length(Uc) == 0) Uc <- U
  }

  # V^(j) の構築
  Vj <- integer(0)
  for (i in U) {
    if (i %in% Uc) next
    if (sum(Aknw[i, Uc], na.rm = FALSE) == 0) {
      Vj <- c(Vj, i)
    }
  }

  return(list(Uc = Uc, Vj = Vj))
}


#' pwling による因果順序の探索
#' @param X データ行列
#' @param U 全変数インデックス
#' @param Uc 候補変数のインデックス
#' @param Vj 事前知識に基づく変数集合
#' @return 選ばれた変数のインデックス
#' @keywords internal
search_causal_order_pwling <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) {
    return(Uc[1])
  }

  M_list <- numeric(length(Uc))

  X_std_matrix <- scale(X)

  for (idx in seq_along(Uc)) {
    i <- Uc[idx]
    M <- 0
    xi_std <- X_std_matrix[, i] # 事前計算済みのベクトルを抽出

    for (j in U) {
      if (i == j) next
      xj_std <- X_std_matrix[, j] # 事前計算済みのベクトルを抽出

      ri_j <- if (i %in% Vj && j %in% Uc) xi_std else residual_vec(xi_std, xj_std)
      rj_i <- if (j %in% Vj && i %in% Uc) xj_std else residual_vec(xj_std, xi_std)

      dm <- diff_mutual_info(xi_std, xj_std, ri_j, rj_i)
      M <- M + min(0, dm)^2
    }
    M_list[idx] <- -1.0 * M
  }

  return(Uc[which.max(M_list)])
}


#' カーネル法による相互情報量
#' @param x1 変数1
#' @param x2 変数2
#' @param param パラメータベクトル (kappa, sigma)
#' @return 相互情報量
#' @keywords internal
mutual_information_kernel <- function(x1, x2, param) {
  kappa <- param[1]
  sigma <- param[2]
  n <- length(x1)

  # カーネル行列の計算
  K1 <- exp(-1 / (2 * sigma^2) * outer(x1, x1, "-")^2)
  K2 <- exp(-1 / (2 * sigma^2) * outer(x2, x2, "-")^2)

  I_n <- diag(n)
  tmp1 <- K1 + n * kappa * I_n / 2
  tmp2 <- K2 + n * kappa * I_n / 2

  K_kappa <- rbind(
    cbind(tmp1 %*% tmp1, K1 %*% K2),
    cbind(K2 %*% K1, tmp2 %*% tmp2)
  )
  D_kappa <- rbind(
    cbind(tmp1 %*% tmp1, matrix(0, n, n)),
    cbind(matrix(0, n, n), tmp2 %*% tmp2)
  )

  sigma_K <- svd(K_kappa, nu = 0, nv = 0)$d
  sigma_D <- svd(D_kappa, nu = 0, nv = 0)$d

  # 数値安定性のため、非常に小さい特異値を除外
  sigma_K <- sigma_K[sigma_K > .Machine$double.eps]
  sigma_D <- sigma_D[sigma_D > .Machine$double.eps]

  return((-1 / 2) * (sum(log(sigma_K)) - sum(log(sigma_D))))
}


#' カーネル法による因果順序の探索
#' @param X データ行列
#' @param U 全変数
#' @param Uc 候補変数
#' @param Vj 事前知識集合
#' @return 選ばれた変数のインデックス
#' @keywords internal
search_causal_order_kernel <- function(X, U, Uc, Vj) {
  if (length(Uc) == 1) {
    return(Uc[1])
  }

  n <- nrow(X)
  if (n > 1000) {
    param <- c(2e-3, 0.5)
  } else {
    param <- c(2e-2, 1.0)
  }

  Tkernels <- numeric(length(Uc))

  for (idx in seq_along(Uc)) {
    j <- Uc[idx]
    Tkernel <- 0
    for (i in U) {
      if (i == j) next
      ri_j <- if (j %in% Vj && i %in% Uc) {
        X[, i]
      } else {
        residual_vec(X[, i], X[, j])
      }
      Tkernel <- Tkernel + mutual_information_kernel(X[, j], ri_j, param)
    }
    Tkernels[idx] <- Tkernel
  }

  return(Uc[which.min(Tkernels)])
}


#' 因果順序から隣接行列を推定
#' @param X 元データ
#' @param causal_order 因果順序 (1-based index のベクトル)
#' @param prior_knowledge 事前知識行列 (NULL可)
#' @return 隣接行列 B (n_features x n_features)
#' @keywords internal
estimate_adjacency_matrix <- function(X, causal_order, prior_knowledge = NULL) {
  n_features <- ncol(X)
  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (idx in seq_along(causal_order)) {
    target <- causal_order[idx]
    if (idx == 1) next # 最初の変数は親なし

    # 因果順序でこの変数より前の変数
    predictors <- causal_order[1:(idx - 1)]

    # 事前知識で制約
    if (!is.null(prior_knowledge)) {
      # prior_knowledge[from, to] == 0 のペアを除外
      # predictors -> target のパスが 0 なら除外
      keep <- sapply(predictors, function(p) {
        val <- prior_knowledge[p, target]
        is.na(val) || val != 0
      })
      predictors <- predictors[keep]
    }

    if (length(predictors) == 0) next

    # OLS 回帰
    y <- X[, target]
    Xp <- X[, predictors, drop = FALSE]
    fit <- stats::lm.fit(x = cbind(1, Xp), y = y)
    coefs <- fit$coefficients[-1] # 切片を除く

    B[target, predictors] <- coefs
  }

  return(B)
}
