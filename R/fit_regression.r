# =============================================================================
# Direct LiNGAM - 隣接行列の推定と回帰バックエンド (OLS / LASSO / Adaptive LASSO)
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
# =============================================================================


# glmnet に渡す lambda の探索範囲。
# デフォルトの自動生成だと早期打ち切りが起こりうるため明示的に指定する。
# fit_lasso() と fit_adaptive_lasso() で共通。
lasso_lambda_seq <- exp(seq(2, -7, length.out = 80))

# Ridge 用の lambda グリッド。LASSO より大きい値が最適になりうるため広めに設定。
ridge_lambda_seq <- exp(seq(6, -7, length.out = 100))


#' glmnet が利用可能か確認する
#'
#' 利用できない場合は、どの回帰手法で必要になったかを示すエラーを出す。
#'
#' @param method glmnet を必要とする回帰手法名（エラーメッセージ用）
#' @keywords internal
check_glmnet_available <- function(method) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(sprintf(
      "Package 'glmnet' is required for reg_method = \"%s\". Please install it.",
      method
    ), call. = FALSE)
  }
}


#' 因果順序から隣接行列を推定
#'
#' @param X 元データ
#' @param causal_order 因果順序 (1-based index のベクトル)
#' @param prior_knowledge 事前知識行列 (NULL可)
#' @param method 回帰手法
#'   "ols"           : 通常の最小二乗法（デフォルト）
#'   "lasso"         : LASSO回帰（glmnet）
#'   "adaptive_lasso": Adaptive LASSO（2段階）
#'   "ridge"         : Ridge回帰（glmnet）
#' @param init_method Adaptive LASSOの初期重みの推定手法
#'   "ols"   :最小二乗法（デフォルト）
#'   "ridge" :Ridge回帰
#' @param lambda LASSO のペナルティ (NULL = 交差検証で自動選択)
#'   "lambda.min" : 予測誤差最小
#'   "lambda.1se" : 1SE ルール（よりスパース）
#'   "AIC"       : AIC最小（CVなし、高速）
#'   "BIC"        : BIC最小（CVなし、高速、最もスパース）デフォルト
#'   "oracle"     : Adaptive LASSO 専用。Ridge では使用不可。
#' @return 隣接行列 B (n_features x n_features)
#' @keywords internal
estimate_adjacency_matrix <- function(X,
                                      causal_order,
                                      prior_knowledge = NULL,
                                      method = "adaptive_lasso",
                                      lambda = "BIC",
                                      init_method = "ols") {
  valid_methods <- c("ols", "lasso", "adaptive_lasso", "ridge")
  if (!(method %in% valid_methods)) {
    stop(sprintf(
      "'method' must be one of: %s.",
      paste(valid_methods, collapse = ", ")
    ))
  }

  # glmnet の確認（OLS は base R のみで動作するため不要）
  if (method != "ols") {
    check_glmnet_available(method)
  }

  n_features <- ncol(X)
  B <- matrix(0, nrow = n_features, ncol = n_features)

  for (idx in seq_along(causal_order)) {
    target <- causal_order[idx]
    if (idx == 1) next

    # 因果順序でこの変数より前の変数
    predictors <- causal_order[1:(idx - 1)]

    # 事前知識で制約
    if (!is.null(prior_knowledge)) {
      keep <- sapply(predictors, function(p) {
        val <- prior_knowledge[target, p]
        is.na(val) || val != 0
      })
      predictors <- predictors[keep]
    }

    if (length(predictors) == 0) next

    y <- X[, target]
    Xp <- X[, predictors, drop = FALSE]

    # --- 回帰手法の分岐 ---
    coefs <- switch(method,
      "ols"            = fit_ols(y, Xp),
      "lasso"          = fit_lasso(y, Xp, lambda = lambda),
      "adaptive_lasso" = fit_adaptive_lasso(y, Xp,
                                            lambda = lambda,
                                            init_method = init_method),
      "ridge"          = fit_ridge_reg(y, Xp, lambda = lambda)
    )

    B[target, predictors] <- coefs
  }

  return(B)
}


#' OLS regression
#'
#' @param y response variable (numeric vector)
#' @param Xp predictor matrix
#' @return coefficient vector (excluding intercept)
#' @keywords internal
fit_ols <- function(y, Xp) {
  fit <- stats::lm.fit(x = cbind(1, as.matrix(Xp)), y = y)
  fit$coefficients[-1]
}


#' Select lambda by information criterion
#'
#' @param glmnet_model a glmnet model object
#' @return list with lambda_AIC_best, lambda_BIC_best, idx_AIC_best,
#'   idx_BIC_best, ic_table
#' @keywords internal
ic_glmnet <- function(glmnet_model) {
  tLL <- glmnet_model$nulldev - deviance(glmnet_model)
  k <- glmnet_model$df
  n <- glmnet_model$nobs
  AIC <- -tLL + 2 * k + 2 * k * (k + 1) / pmax(n - k - 1, 1)
  BIC <- log(n) * k - tLL
  ic_table <- data.frame(
    lambda = glmnet_model$lambda,
    df     = k,
    AIC   = AIC,
    BIC    = BIC
  )
  idx_AIC_best <- which.min(ic_table$AIC)
  idx_BIC_best <- which.min(ic_table$BIC)
  list(
    lambda_AIC_best = ic_table$lambda[idx_AIC_best],
    lambda_BIC_best  = ic_table$lambda[idx_BIC_best],
    idx_AIC_best     = idx_AIC_best,
    idx_BIC_best     = idx_BIC_best,
    ic_table         = ic_table
  )
}


#' LASSO 回帰（情報量基準 or CV でラムダ選択）
#'
#' @param y 目的変数
#' @param Xp 説明変数行列
#' @param lambda ラムダ選択方法
#'   "lambda.min" : CV予測誤差最小
#'   "lambda.1se" : CV 1SEルール
#'   "AIC"       : AIC最小
#'   "BIC"        : BIC最小。デフォルト
#' @return 係数ベクトル
#' @keywords internal
fit_lasso <- function(y, Xp, lambda = "BIC") {
  if (ncol(Xp) == 1) {
    return(fit_ols(y, Xp))
  }

  check_glmnet_available("lasso")

  Xp_mat <- as.matrix(Xp)

  if (lambda %in% c("AIC", "BIC")) {
    fit <- glmnet::glmnet(
      x = Xp_mat, y = y,
      alpha = 1, intercept = TRUE, standardize = TRUE,
      lambda = lasso_lambda_seq
    )

    ic <- ic_glmnet(fit)
    # 選んだ lambda は fit$lambda 上の1点なので、coef(fit, s = ...) の
    # 補間を経由せず列インデックスで直接取り出す（結果は同一で高速）
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])

  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y,
      alpha = 1, intercept = TRUE, standardize = TRUE,
      lambda = lasso_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
    coef_vec <- as.numeric(stats::coef(cv_fit, s = lambda_val))[-1]
  }

  return(coef_vec)
}


#' Ridge 回帰（情報量基準 or CV でラムダ選択）
#'
#' @param y 目的変数
#' @param Xp 説明変数行列
#' @param lambda ラムダ選択方法
#'   "lambda.min" : CV予測誤差最小
#'   "lambda.1se" : CV 1SEルール
#'   "AIC"       : AIC最小
#'   "BIC"        : BIC最小。デフォルト
#'   "oracle"は使用不可（Adaptive LASSO 専用）。
#' @return 係数ベクトル
#' @keywords internal
fit_ridge_reg <- function(y, Xp, lambda = "BIC") {
  if (ncol(Xp) == 1) {
    return(fit_ols(y, Xp))
  }

  if (lambda == "oracle") {
    stop("lambda = \"oracle\" is only supported for reg_method = \"adaptive_lasso\".",
         call. = FALSE)
  }

  check_glmnet_available("ridge")

  Xp_mat <- as.matrix(Xp)

  if (lambda %in% c("AIC", "BIC")) {
    fit <- glmnet::glmnet(
      x = Xp_mat, y = y,
      alpha = 0, intercept = TRUE, standardize = TRUE,
      lambda = ridge_lambda_seq
    )

    ic <- ic_glmnet(fit)
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])

  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y,
      alpha = 0, intercept = TRUE, standardize = TRUE,
      lambda = ridge_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
    coef_vec <- as.numeric(stats::coef(cv_fit, s = lambda_val))[-1]
  }

  return(coef_vec)
}


#' Adaptive LASSO
#' @param y 目的変数
#' @param Xp 説明変数行列
#' @param lambda ラムダ選択方法 ("lambda.min", "lambda.1se", "AIC", "BIC", "oracle")
#' @param gamma_weight 重みの指数
#' @param init_method 初期重みの推定手法 ("ols" または "ridge")
#' @return 係数ベクトル
#' @keywords internal
fit_adaptive_lasso <- function(y, Xp,
                               lambda = "BIC",
                               gamma_weight = 1.0,
                               init_method = "ols") {
  if (ncol(Xp) == 1) return(fit_ols(y, Xp))
  check_glmnet_available("adaptive_lasso")

  Xp_mat <- as.matrix(Xp)
  n <- nrow(Xp_mat) # サンプルサイズ n を取得

  # --- Step 1: 初期推定量 (元のスケール) の計算 ---
  if (init_method == "ols") {
    init_fit <- stats::lm.fit(x = cbind(1, Xp_mat), y = y)
    init_coefs <- as.numeric(init_fit$coefficients[-1])
  } else {
    ridge_cv <- glmnet::cv.glmnet(
      x = Xp_mat, y = y, alpha = 0,
      intercept = TRUE, standardize = TRUE
    )
    init_coefs <- as.numeric(stats::coef(ridge_cv, s = "lambda.min"))[-1]
  }
  init_coefs[is.na(init_coefs)] <- 0

  # --- Step 2: penalty.factor の計算 ---
  x_sds <- apply(Xp_mat, 2, sd_pop)
  x_sds[x_sds < 1e-10] <- 1e-10

  init_coefs_std <- init_coefs * x_sds

  pf <- 1 / (abs(init_coefs_std)^gamma_weight)
  pf[is.infinite(pf) | is.na(pf)] <- 1e10

  # --- Step 3: Adaptive LASSO の実行 ---
  fit <- glmnet::glmnet(
    x = Xp_mat, y = y, alpha = 1,
    intercept = TRUE, standardize = TRUE,
    penalty.factor = pf,
    lambda = lasso_lambda_seq
  )

  if (lambda %in% c("AIC", "BIC")) {
    ic <- ic_glmnet(fit)
    # fit$lambda 上の1点なので補間せず列インデックスで直接取り出す
    k_best <- if (lambda == "AIC") ic$idx_AIC_best else ic$idx_BIC_best
    coef_vec <- as.numeric(fit$beta[, k_best])
    return(coef_vec)
  }

  if (lambda == "oracle") {
    # オラクル lambda は探索グリッド上にないため coef() で補間する
    lambda_val <- 5 / (n^(1.75))
  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = Xp_mat, y = y, alpha = 1,
      intercept = TRUE, standardize = TRUE,
      penalty.factor = pf,
      lambda = lasso_lambda_seq
    )

    lambda_val <- cv_fit[[lambda]]
  }

  coef_vec <- as.numeric(stats::coef(fit, s = lambda_val))[-1]

  return(coef_vec)
}
