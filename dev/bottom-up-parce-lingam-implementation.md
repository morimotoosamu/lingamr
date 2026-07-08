# BottomUpParceLiNGAM 実装指示書

lingamr パッケージに BottomUpParceLiNGAM（潜在交絡に頑健な因果順序探索）を追加する。
Python cdt15/lingam の `BottomUpParceLiNGAM`（`bottom_up_parce_lingam.py`、全618行）の R 移植。
**本指示書は原典618行と、依存する `hsic.py`（全165行）・`utils/_f_correlation.py`（全145行）・
`utils.predict_adaptive_lasso` を全文精読した上で書かれている。**

アルゴリズムの要点: 下流（sink）側から順に「残りの変数で回帰した残差が説明変数と独立か」を
検定しながら因果順序を決める。検定が棄却された時点で探索を打ち切り、**順序を決められなかった
変数集合（潜在交絡の疑い）を「順序不明ブロック」として返す**。隣接行列ではブロック内ペアが
NA になる。この「NA を含む隣接行列」と「ブロック入り causal_order」が既存クラスとの最大の違い。

## 0. 参照資料（実装前に必ず全文を読むこと）

- ソース: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/bottom_up_parce_lingam.py
- HSIC 検定（**lingamr に未実装。新規移植が必要**）:
  https://raw.githubusercontent.com/cdt15/lingam/master/lingam/hsic.py
- F-correlation: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/utils/_f_correlation.py
- `predict_adaptive_lasso`（utils/__init__.py 756-794行目付近）:
  https://raw.githubusercontent.com/cdt15/lingam/master/lingam/utils/__init__.py
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/bottom_up_parce.html
- 論文: T. Tashiro, S. Shimizu, A. Hyvärinen, T. Washio.
  ParceLiNGAM: a causal ordering method robust against latent confounders.
  Neural Computation, 26(1): 57-83, 2014.

着手前に必ず読む lingamr 側ファイル:

- `R/lingam_direct.r` — バリデーション・roxygen2 の様式
- `R/search_causal_order.r` — `extract_partial_orders`（再利用可）、
  `incomplete_cholesky_gauss`（**再利用不可**。理由は 3.6 参照）
- `R/fit_regression.r` — 隣接行列推定の回帰バックエンド（ターゲット単位で再利用）
- `R/lingam_bootstrap.r` — bootstrap の並列・再現性・失敗処理パターン（踏襲）と
  `create_bootstrap_result()`
- `R/paths.r` — パスベースの総合効果（bootstrap 内で再利用）
- `tests/testthat/test-lingam_direct.R`, `test-bootstrap.R` — テストパターン

## 1. スコープ（ユーザー確認済み: fcorr・bootstrap とも含める）

### 実装するもの

| 成果物 | ファイル |
|---|---|
| HSIC ガンマ近似検定（内部ユーティリティ） | `R/hsic.r`（新規） |
| F-correlation（内部ユーティリティ） | `R/f_correlation.r`（新規） |
| `lingam_parce()` + `ParceLingamResult` + print | `R/lingam_parce.r`（新規） |
| `estimate_total_effect_parce()` / `get_error_independence_p_values_parce()` | 同上 |
| `lingam_parce_bootstrap()` | `R/lingam_parce_bootstrap.r`（新規） |
| `generate_parce_sample()`（潜在交絡入りサンプル） | `R/generate_parce_sample.r`（新規） |
| testthat テスト | `tests/testthat/test-hsic.R`, `tests/testthat/test-lingam_parce.R`（新規） |
| man/NAMESPACE 再生成・NEWS.md・vignette・_pkgdown.yml | 既存様式で追記 |

### 実装しないもの（設計判断込み）

- `regressor` 引数（Python は sklearn 互換 regressor を差し替え可能。原典 288-302 行目の
  既定動作＝共分散行列の擬似逆行列による OLS 残差のみを実装し、引数は設けない。
  Details に記載）
- `random_state` 引数（Python 側でも fit では未使用。乱数は bootstrap の seed のみ）
- `hsic_test_gamma` の `bw_method="scott"/"silverman"`（statsmodels 依存。
  BottomUpParce からは常に既定の "mdbs"＝中央値距離で呼ばれるため mdbs のみ実装）
- 新規パッケージ依存なし（HSIC・F-correlation とも base R + stats で書ける）

## 2. 新規ユーティリティ

### 2.1 HSIC ガンマ近似検定（R/hsic.r、原典 hsic.py の忠実移植）

将来 RCD 等でも使うため独立ファイルに置く。全て `@keywords internal`。

```r
hsic_kernel_width(X)      # get_kernel_width: 中央値距離バンド幅
hsic_gram_matrix(X, width)  # get_gram_matrix: RBF グラム行列 K と中心化 Kc
hsic_test_gamma(X, Y)     # 検定統計量とガンマ近似 p 値を返す list(stat =, p =)
```

精読済みの数式（hsic.py 16-165 行目）:

- バンド幅: 先頭 **100 点まで**で点間平方距離を計算し、
  `sqrt(0.5 * median(正の距離))`（下三角を除外し dists > 0 のみの median。
  R では `dist()` か outer で可。**「先頭100点」の切り出しは X[1:100, ] であって
  ランダムサンプリングではない**点を原典どおりに）
- グラム行列: `K <- exp(-H / (2 * width^2))`（H は平方距離行列）。
  中心化 `Kc <- K - (colSums + rowSums)/n + sum(K)/n^2`（rank-1 補正の形。原典 74-77 行目）
- 統計量: `stat <- sum(t(Kc) * Lc) / n`
- ガンマ近似（原典 147-163 行目を式ごと忠実に）:
  ```
  var  <- (Kc * Lc / 6)^2
  var  <- (sum(var) - sum(diag(var))) / n / (n - 1)
  var  <- 72 * (n-4) * (n-5) / n / (n-1) / (n-2) / (n-3) * var
  K, L の対角を 0 にして mu_X <- sum(K) / n / (n-1)（mu_Y も同様）
  mean <- (1 + mu_X * mu_Y - mu_X - mu_Y) / n
  alpha <- mean^2 / var;  beta <- var * n / mean
  p <- pgamma(stat, shape = alpha, scale = beta, lower.tail = FALSE)   # gamma.sf 相当
  ```
- **性能注意**: n×n のグラム行列を毎回作る O(n²) 実装（原典どおり）。roxygen2 の
  Details に「n が数千を超えると遅い」と記載。最適化はしない（忠実性優先）

テスト（test-hsic.R）: 固定シードの独立データで p が大きく、強い非線形従属データ
（例: y = x^2 + 小ノイズ）で p が 0 近傍になること。同一入力で決定的であること。

### 2.2 F-correlation（R/f_correlation.r、原典 _f_correlation.py の忠実移植）

Bach & Jordan (2002) のカーネル正準相関。`f_correlation(x, y)` は [0, 1] 近傍の値を返し、
大きいほど従属。全て `@keywords internal`。

- 入力2変数を**母標準偏差**で標準化（原典 39-40 行目の `.std()` は母 SD = `sd_pop()`）
- `n > 1000` なら `kappa = 2e-3, sigma = 0.5`、以下なら `kappa = 2e-2, sigma = 1.0`
  （lingamr の kernel measure と同じ切替値）
- 不完全 Cholesky → 中心化 → SVD（固有分解）→ R_kappa 行列 → **最小固有値** λ を求め
  `1 - λ` を返す
- SVD 部（原典 131-145 行目）: `eigh(t(G) %*% G)` の固有値 D のうち `D >= n * kappa * 1e-2`
  を降順に採用し、`U <- G %*% A[, idx] %*% diag(1/sqrt(D))`、
  `R[j] <- D[j] / (n * kappa / 2 + D[j])`
- R_kappa: ブロック行列（対角 = 単位行列、非対角ブロック = `diag(R_i) (U_i' U_j) diag(R_j)`）。
  最小固有値は `eigen(R_kappa, symmetric = TRUE)` の最小値（scipy eigsh(SM) 相当）

**【重要】既存 `incomplete_cholesky_gauss` は再利用しないこと**
（search_causal_order.r:328-371）。停止条件が異なる:
既存版は「最大対角残差 <= tol（tol=1e-4）+ max_rank 上限」、原典 _f_correlation 版は
「**対角残差の総和** > tol（tol = n * kappa * 1e-2）の間続行、上限なし」（原典 91-128 行目）。
ランクが変わると数値が変わるため、f_correlation.r 内に原典準拠の
`incomplete_cholesky_fcorr(x, sigma, tol)` を別途実装する（既存関数の構造は参考にしてよい。
既存版同様、G の行を元のインデックス順で埋めれば原典の並べ替え（56行目）は不要になる —
ただし原典 57 行目の**列平均による中心化**を忘れないこと）。

## 3. 本体 `lingam_parce()`（R/lingam_parce.r）

### 3.1 API

```r
lingam_parce <- function(X,
                         alpha = 0.1,
                         prior_knowledge = NULL,
                         independence = "hsic",
                         ind_corr = 0.5,
                         reg_method = "adaptive_lasso",
                         lambda = "BIC",
                         init_method = "ols")
```

- `alpha`: 独立性検定の有意水準。**`alpha = 0` なら棄却が起きず完全な順序を返す**
  （原典 docstring 46 行目）。負値はエラー
- `independence`: `match.arg(c("hsic", "fcorr"))`
- `ind_corr`: fcorr の独立判定閾値（f_correlation 値がこれ**以上**で棄却。既定 0.5、負値エラー）
- `reg_method`/`lambda`/`init_method`: 隣接行列推定の回帰指定。lingam_direct と同じ選択肢・
  既定値（Python は adaptive lasso + BIC 固定 = `predict_adaptive_lasso`。R では既存
  fit_regression.r の資産があるため選択可能にする。既定が Python と同値であることを Details に記載）
- `apply_prior_knowledge_softly` は**設けない**（Python 版 Parce に無い）

バリデーションは lingam_direct.r:73-93 の様式。prior_knowledge の前処理
（`< 0 → NA`、`extract_partial_orders()` 再利用）は lingam_direct.r:98-109 と同じ。

### 3.2 fit の流れ（原典 92-142 行目）

```
1. X を列ごとに中心化（原典 124 行目）
2. Bonferroni 補正: thresh_p <- alpha / (p - 1)   （p = 全変数数で固定。原典 127 行目）
3. ボトムアップ探索（3.3）→ K_bttm（sink 側から確定した順序）と p_bttm（各段の検定 p 値）
4. U_res <- 確定できなかった変数集合
   causal_order <- 「U_res が2個以上なら先頭に順序不明ブロック」+ K_bttm
5. 隣接行列推定（3.5）
```

### 3.3 ボトムアップ探索 `_search_causal_order`（原典 188-225 行目）

```
U <- 1:p;  K_bttm <- integer(0);  p_bttm <- numeric(0)
repeat {
  Uc <- parce_search_candidate(U, partial_orders)   # 3.4
  (m, eval) <- find_exo_vec(X, Uc, U)               # 最も sink らしい候補と検定値
  if (棄却されない) {
    K_bttm <- c(m, K_bttm)          # 先頭に追加（bottom-up なので新しいものほど上流側に積む）
    p_bttm <- c(eval, p_bttm)
    U <- setdiff(U, m)
    partial_orders から to 列 == m の行を削除     # 原典 211-214 行目（[:, 1] != m）
    if (length(U) <= 1) { K_bttm <- c(U, K_bttm); p_bttm <- c(0.0, p_bttm); break }
  } else break                       # 棄却 → 残りは順序不明ブロック
}
```

棄却判定（原典 278-286 行目）: hsic は `eval < thresh_p` で棄却（eval = Fisher 合成 p 値）、
fcorr は `eval >= ind_corr` で棄却（eval = 最大 f_correlation）。

**prior knowledge の候補選択は DirectLiNGAM と逆向き**（原典 175-186 行目）:
sink を探すので「partial_orders の from 列（1列目）に現れる変数は候補から除く」だけの
単純なフィルタ。既存 `search_candidate()`（top-down 用・Vj あり）は使わず、
`parce_search_candidate <- function(U, partial_orders)` を新規内部関数として書く。
partial_orders 自体は既存 `extract_partial_orders()` の出力（[from, to] 形式）をそのまま使う。

### 3.4 sink 候補の探索 `find_exo_vec`（原典 227-276 行目）

各候補 j ∈ Uc について:

1. `xi_index <- setdiff(Uc, j)`（**U ではなく Uc** から除く。読み違えやすい）、
   `xj <- j`
2. 残差: `R <- X[, j] - X[, xi_index] %*% coef`。coef は**共分散行列の擬似逆行列**で計算
   （原典 288-302 行目）:
   `coef <- pinv(cov(X)[xi_index, xi_index]) %*% cov(X)[j, xi_index]`
   pinv は svd ベースで自前実装（HighDim 指示書 dev/high-dim-direct-lingam-implementation.md
   の 3.2 に同じものがある。同一の内部関数を実装してよい）。
   注: `np.cov` は不偏（n-1 割り）だが pinv を通すので分母は係数に影響しない
3. **hsic の場合**: 各説明変数 X[, i]（i ∈ xi_index）と R の HSIC p 値を Fisher 法で合成:
   `fisher_stat <- Σ -2 log(p_i)`（p_i = 0 なら Inf）、
   `fisher_p <- 1 - pchisq(fisher_stat, df = 2 * length(xi_index))`。
   説明変数が1個のときは合成せず HSIC の (stat, p) をそのまま使う（原典 309-310 行目）。
   **早期打ち切り**: 現在の最良 `max_p_stat` を超えたらループを break（原典 316-317 行目。
   統計量が小さい候補を探しているため）。
   最小の fisher_stat を持つ候補を選び、その fisher_p を eval とする
4. **fcorr の場合**: 各説明変数と R の f_correlation の**最大値**を取り、
   それが最小の候補を選ぶ。eval = その最大 f_correlation
5. `length(Uc) == 1` の特例（原典 233-243 行目): 候補は確定。
   `predictors <- setdiff(U, Uc)`（**ここは U ベース**）で残差を作り eval のみ計算

### 3.5 隣接行列推定（原典 349-394 行目。既存関数は流用できないので専用実装）

`causal_order` は「先頭に順序不明ブロック（整数ベクトル）+ 以降スカラー」のリスト構造。

```
B <- 0 行列 (p × p)
for (i in 2:length(causal_order_list)) {
  target <- causal_order_list[[i]]                    # スカラー
  predictors <- unlist(causal_order_list[1:(i-1)])    # ブロックは展開（原典 _flatten）
  prior knowledge があれば pk[target, predictors] == 0 の predictor を除外
      # 原典 365-367, 377-378 行目。対角は 0 に潰してから
  if (length(predictors) > 0)
    B[target, predictors] <- ターゲット単位の adaptive lasso 係数
}
ブロック内の全ペア (xi, xj) について B[xi, xj] <- B[xj, xi] <- NA   # 原典 384-391 行目
```

- ターゲット単位の回帰は **R/fit_regression.r の内部関数を再利用**する。
  `estimate_adjacency_matrix()` がターゲットごとにどう回帰backend
  （fit_adaptive_lasso 等）を呼んでいるかを読み、同じ呼び方をする
  （estimate_adjacency_matrix 自体は厳密な順序ベクトル前提なので丸ごとは使えない）。
  Python の `predict_adaptive_lasso`（標準化 → OLS 重み → LassoLarsIC(BIC) →
  元スケールで OLS 再フィット）と lingamr の adaptive_lasso + lambda="BIC" が
  同等の設計であることを確認して使う
- **ブロック内の変数はターゲットにならない**（親を推定しない）。ブロック変数は
  下流ターゲットの predictor としてのみ現れる。ブロック内相互は NA

### 3.6 返り値: S3 クラス `ParceLingamResult`

```r
result <- list(
  adjacency_matrix = B,            # p×p, lingamr 規約 B[i,j] = j -> i。順序不明ペアは NA
  causal_order     = K,            # リスト。各要素は整数ベクトル（長さ1 = 確定、長さ>1 = 順序不明ブロック）
  p_values         = p_bttm,       # 各確定ステップの検定 p 値（bottom-up で確定した順に対応。診断用）
  independence     = independence
)
class(result) <- "ParceLingamResult"
```

- **向き**: 原典は `B[target, predictors] <- coef`（382行目）なので lingamr 規約
  `B[i,j] = j→i` とそのまま一致。**転置不要**
- ブロックが空（全順序確定）のときも causal_order は「長さ1ベクトルのリスト」で統一する
  （構造が入力依存で変わると下流処理が書きにくいため）
- `print.ParceLingamResult()`: ヘッダ "Bottom-Up ParceLiNGAM Result"。causal_order は
  ブロックを括弧で表示（例: `(x2, x3) -> x0 -> x5 -> x1 -> x4`）。
  「NA = 順序不明（潜在交絡の疑い）」の注記を1行入れる。隣接行列表示は
  print.LingamResult の書式に倣う

### 3.7 総合効果と誤差独立性（原典 396-544 行目）

既存の `estimate_total_effect()` / `get_error_independence_p_values()` は
`validate_lingam_result()` で LingamResult 専用のため使えない。Parce 専用に新設する:

- `estimate_total_effect_parce(X, parce_result, from_index, to_index)`（原典 396-458 行目）:
  - from が to より因果順序で後なら warning（ブロックは同順位として扱う。原典 418-439 行目）
  - **from の行に NA がある（= from が交絡を受けている）場合は warning を出して NA を返す**
    （原典 442-447 行目）
  - それ以外は「from + from の親」を説明変数に to を OLS 回帰し、from の係数を返す
    （既存 estimate_total_effect.r の実装と同じ方式。内部を流用してよい）
- `get_error_independence_p_values_parce(X, parce_result)`（原典 510-544 行目）:
  - `E <- X - X %*% t(B)`（NA を含む行の残差は計算不能）
  - B に NA を含む行・列の変数が絡むペアの p 値は NA
  - それ以外のペアは **HSIC 検定**（2.1 の hsic_test_gamma）で p 値行列を埋める。
    既存 get_error_independence_p_values（相関ベース）と方式が違うことを Details に記載

## 4. bootstrap `lingam_parce_bootstrap()`（原典 572-618 行目）

`R/lingam_bootstrap.r`（87-284行目）の構造（引数検証 → run_one + tryCatch →
PSOCK + `.libPaths` 伝搬 + `clusterSetRNGStream` → 失敗スキップ → 集計）を踏襲する。
引数は `lingam_direct_bootstrap` に準じ、`measure` の代わりに
`alpha / independence / ind_corr` を取る。

Parce 固有の注意:

1. **総合効果はパス積和**（原典 606-616 行目は `estimate_total_effect2` = 
   `calculate_total_effect(B, from, to)`）。R/paths.r の既存パス列挙・効果計算を再利用。
   from が交絡行（NA あり）のときは NA（原典 500-505 行目）。
   ブロックが from 側のときはブロックの**各メンバー**を from として計算（原典 608-612 行目）
2. **NA の集計互換性**: numpy では `np.abs(nan) > threshold` が False になるため、
   Python の BootstrapResult 集計では NA 辺は「辺なし」として数えられる。
   R の既存 `get_probabilities()` 等は NA が混ざると NA を返しかねないので、
   **BootstrapResult に格納する隣接行列・総合効果は NA → 0 に置換してから格納する**
   （numpy の比較セマンティクスと同じ集計結果になる）。この仕様を roxygen2 の
   Details に明記する
3. `create_bootstrap_result()`（lingam_bootstrap.r:299-308）を再利用。ただし
   **causal_orders は渡さない（NULL）**。Parce の因果順序はブロック入りで
   n_features 長の整数行列に格納できないため。結果として
   `get_causal_order_stability()` はこの結果には使えない — Details に明記
4. 返り値は通常の `BootstrapResult`（Python も単一の BootstrapResult を返す。618行目）。
   `get_probabilities()` / `get_causal_direction_counts()` /
   `get_directed_acyclic_graph_counts()` / `get_total_causal_effects()` が動くことを
   テストで確認する

## 5. データ生成関数 `generate_parce_sample()`

チュートリアル準拠の**潜在交絡入り**7変数モデル（観測は6変数）:

```
x6 (潜在) ~ Uniform
x3 = 2.0*x6 + e;  x2 = 2.0*x6 + e          # x6 が x2, x3 の共通原因
x0 = 0.5*x3 + e
x1 = 0.5*x0 + 0.5*x2 + e
x5 = 0.5*x0 + e
x4 = 0.5*x0 - 0.5*x2 + e                    # e はすべて一様乱数
```

- シグネチャ例: `generate_parce_sample(n = 1000, seed = NULL)`
- 返り値: `list(data = <x0..x5 の観測 data.frame>, adjacency_matrix = <観測変数間の真の B。
  x2-x3 間は交絡ペアとして NA>, confounded_pair = c(3, 4)  # x2, x3 の 1-based 列位置)`
  （既存 generate_* の返り値形式を確認して合わせる）
- 期待される推定結果（チュートリアル、seed 依存）: causal_order = `[(x2, x3), x0, x5, x1, x4]`
  相当、B の x2-x3 間が NA

## 6. テスト（test-lingam_parce.R。HSIC 単体は test-hsic.R）

数値スナップショット不使用。n は実行時間と検定力のバランスで調整（目安 500-1000、
HSIC は O(n²) なので大きくしすぎない）。

1. 構造: `expect_s3_class(res, "ParceLingamResult")`, `expect_named()`,
   causal_order がリストで全変数を過不足なく含む、B の次元・dimnames
2. バリデーション: 非数値・NA 入り X、`alpha = -1`、`independence = "foo"`、
   `ind_corr = -1`、prior_knowledge の次元不正、それぞれ `expect_error()`
3. **潜在交絡の検出（最重要）**: `generate_parce_sample(n = 1000, seed 固定)` で
   - causal_order の先頭ブロックに x2, x3 が入る
   - `B[x2, x3]` と `B[x3, x2]` が NA
   - ブロック外の真の辺（x3→x0, x0→x5 等）が符号込みで回復される
   （検定ベースで seed に敏感なら、複数 seed を試して安定するものに固定してよい。
   その場合もテスト内の seed は1つに固定）
4. **alpha = 0 で棄却なし**: ブロックが生じず、全要素が長さ1の causal_order になる
5. **交絡なしデータとの整合**: `generate_lingam_sample_6()` のような交絡なしデータで、
   ブロックなしの順序が返り、真の構造とおおむね整合する
6. `independence = "fcorr"`: 同じ交絡データで valid な結果が返る（ブロック検出まで
   要求するかは検定力を見て判断。最低限エラーなく ParceLingamResult が返ること）
7. `estimate_total_effect_parce`: 交絡を受ける変数を from にすると warning + NA、
   健全なペアでは数値が返る
8. `get_error_independence_p_values_parce`: NA 絡みペアが NA、他が [0,1] の p 値
9. print: `expect_output(print(res), "ParceLiNGAM")` とブロックの括弧表示
10. bootstrap: `n_sampling = 5-10, seed = 42` で BootstrapResult が返り、
    `get_probabilities()` 等が NA を返さず動く。逐次再現性。並列（n_cores = 2）は
    `skip_on_cran()` 付きで再現性確認
11. prior_knowledge: `make_prior_knowledge()` の pk を渡して動作すること

実行: `devtools::test(filter = "parce")` と `devtools::test(filter = "hsic")` → 全体

## 7. ドキュメント

- roxygen2: 既存様式。`@references` は Tashiro et al. (2014) Neural Computation 26(1):57-83
- Details に記載する事項:
  1. 潜在交絡に対する頑健性の考え方と、順序不明ブロック・NA の意味
  2. `alpha = 0` の挙動（棄却なし＝完全順序）
  3. HSIC は O(n²)。n が数千超で遅くなる
  4. fcorr の判定は p 値でなく F-correlation 値の閾値（`ind_corr`）であること
  5. 誤差独立性検定が HSIC ベースで、既存 `get_error_independence_p_values()`
     （相関ベース）と方式が異なること
  6. bootstrap 集計では NA（順序不明）辺を辺なしとして扱うこと、
     `get_causal_order_stability()` が使えないこと
  7. Python 版の `regressor` / `bw_method` 引数を持たないこと
- `@examples`: `generate_parce_sample()` → fit → print → total effect の NA warning 実演。
  HSIC が重いので n は小さめ、遅ければ `\donttest{}`
- NEWS.md 追記、vignette 新セクション「Latent Confounders: Bottom-Up ParceLiNGAM」
  （動機 → 生成 → fit → ブロックと NA の読み方 → bootstrap）、_pkgdown.yml に
  `lingam_parce`, `lingam_parce_bootstrap`, `estimate_total_effect_parce`,
  `get_error_independence_p_values_parce`, `generate_parce_sample` を追加

## 8. 完了条件・検証手順

R は PATH 外の `C:\R\R-4.6.0` にある。実行例:

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] `devtools::document()` が警告なく通る
- [ ] `devtools::test()` 全緑（既存テストを壊さない）
- [ ] `devtools::check()` で新規 ERROR/WARNING/NOTE なし
- [ ] 潜在交絡検出テスト（セクション6-3）が通る
- [ ] alpha = 0 テスト（セクション6-4）が通る
- [ ] bootstrap の再現性・NA 非混入テストが通る
- [ ] vignette がビルドできる

## 9. 遵守事項

- **既存の R ファイルは一切変更しない**（extract_partial_orders・fit_regression の
  内部関数・create_bootstrap_result・paths.r は同一パッケージ内なのでそのまま呼べる）
- NEWS.md / vignettes/lingamr.Rmd / _pkgdown.yml の編集は追記のみ
- テストコードの削除・コメントアウト禁止
- commit / push はしない（ユーザーが行う）
- エラーは握りつぶさない。意味のあるメッセージ付きで stop/warning
- 範囲外のリファクタリング禁止
