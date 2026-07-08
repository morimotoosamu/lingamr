# LiM (LiNGAM for Mixed data) 実装指示書（v2: 原典精読済み）

lingamr パッケージに LiM アルゴリズムを追加する。Python cdt15/lingam の `LiM` クラス
（`lim.py`、全411行）の R 移植。**本指示書は原典411行と `lingam/utils/__init__.py` の
`likelihood_i`（866行目〜）を全文精読した上で書かれている。**

## 0. 参照資料（実装前に必ず全文を読むこと）

実装時も必ず原典を WebFetch で取得し、原典を正とすること。

- ソース: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/lim.py
- `likelihood_i` / `log_p_super_gaussian` / `variance_i`（連続変数の尤度。866-939行目付近）:
  https://raw.githubusercontent.com/cdt15/lingam/master/lingam/utils/__init__.py
  （注意: チュートリアルに登場する `simulate_linear_mixed_sem` / `simulate_dag` は
  この utils には**存在しない**。チュートリアル用のローカルスクリプト由来なので探さないこと。
  データ生成関数はセクション5のモデル定義から設計する）
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/lim.html
- 論文: Y. Zeng, S. Shimizu, H. Matsui, F. Sun. Causal discovery for linear mixed data.
  In Proc. First Conference on Causal Learning and Reasoning (CLeaR2022). PMLR 177, pp. 994-1009, 2022.

着手前に必ずリポジトリの CLAUDE.md と以下の既存実装を読み、規約を踏襲すること:

- `R/lingam_direct.r` — 関数シグネチャ、バリデーション、S3 クラス、print、roxygen2 の手本
- `R/generate_lingam_sample.r` — データ生成関数の手本
- `tests/testthat/test-lingam_direct.R` — テストパターンの手本

## 1. スコープ

### 実装するもの

| 成果物 | ファイル |
|---|---|
| `lingam_lim()` + `print.LiMResult()` | `R/lingam_lim.r`（新規） |
| `generate_lim_sample()`（テスト・例示用データ生成） | `R/generate_lim_sample.r`（新規） |
| testthat テスト | `tests/testthat/test-lingam_lim.R`（新規） |
| roxygen2 → man/NAMESPACE 再生成 | `devtools::document()` |
| NEWS.md 追記 | `NEWS.md` |
| vignette 追記（LiM セクション） | `vignettes/lingamr.Rmd` |
| pkgdown reference 追記 | `_pkgdown.yml` |

### 実装しないもの

- Poisson 経路（原典の `loss_type="poisson"` 分岐 102-108行目、`is_poisson=True` 分岐
  188-221行目、`_factorial` 271-277行目）。移植しない。離散変数は **2値（0/1）のみ**
  サポートし、roxygen2 にその旨を明記する
- `loss_type` 引数の公開（原典の実運用は "mixed" 固定。"logistic"/"laplace" 単独分岐も
  移植しない。local phase の "mixed_dag" は内部実装として持つ）
- bootstrap 対応（Python 版 LiM にも無い）
- 新規パッケージ依存。**base R + stats のみ**（`optim(method="L-BFGS-B")`、`glm`、`lm`）。
  networkx / pandas / sklearn 相当はすべて自前 or stats で代替（後述）

## 2. API 設計

```r
lingam_lim <- function(X,
                       is_continuous,
                       lambda1 = 0.1,
                       max_iter = 150L,
                       h_tol = 1e-8,
                       rho_max = 1e16,
                       w_threshold = 0.1,
                       only_global = FALSE)
```

- `X`: 数値行列または data.frame。列名を保存する
- `is_continuous`: 長さ ncol(X) の logical。TRUE=連続、FALSE=離散(2値)。
  Python の `dis_con`（1×p 行列, 1=連続, 0=離散）の置き換え
- 既定値はすべて Python と同値

### バリデーション（lingam_direct.r:73-93 のスタイルを踏襲）

- `X` の数値・NA・次元チェック（`ncol >= 2`, `nrow >= 2`）
- `is_continuous`: logical、長さ一致、NA なし
- 離散指定列（`!is_continuous`）の値が {0, 1} の2値であることをチェック
  （原典の logistic 損失 `logaddexp(0, M) - X*M` は X が 0/1 であることを前提とする）
- スカラー引数の型・範囲チェック。エラーは英語メッセージ + `stop(call. = FALSE)`

### 返り値: S3 クラス `LiMResult`

```r
result <- list(
  adjacency_matrix = B,        # p×p, lingamr 規約 B[i,j] = j -> i（転置済み。下記）
  causal_order     = K,        # B から導出したトポロジカル順序（1-based）
  is_continuous    = is_continuous
)
class(result) <- "LiMResult"
```

- `adjacency_matrix` の dimnames に X の列名
- `print.LiMResult()` は `print.LingamResult()` に倣う。ヘッダ "LiM Result"、
  変数型（continuous/discrete）を1行表示、"Adjacency matrix (row = to, col = from):"、
  `invisible(x)`
- causal_order は Kahn 法のトポロジカルソートで導出。巡回時（w_threshold や local phase の
  制約上ほぼ起きないが）は `NA` + warning

### 【確定】隣接行列の向き — 転置が必要

原典精読により確定。内部 W の損失は `M = X @ W`（96-98行目）であり、M の列 j =
変数 j の予測値。つまり **内部 W は `W[i,j] = i → j`（NOTEARS 規約）**。
`_bic_scoring` の `nx.DiGraph(W)` → `predecessors(i)` の使い方（147-156行目）も同じ規約で一貫。
Python の `adjacency_matrix_` はこの W をそのまま返す。

**R では最終結果を `B <- t(W_min_lss)` と転置して lingamr 規約（`B[i,j] = j→i`）に変換すること。**
roxygen2 の Details に「Python 版 `adjacency_matrix_` とは転置の関係」と明記。
単体テストでも検証する（セクション6）。

### 【重要な仕様注意】出力の非ゼロ値は純粋な係数ではない

local phase で辺を反転・追加するとき、原典は係数ではなく **1.0 を代入する**
（338行目 `W_tmp[I[1][iii], I[0][iii]] = 1`、374行目 `W_tmp[...] = 1`）。
したがって最終行列の非ゼロ値は「global 段階の推定係数」と「local 段階で置かれた 1」の混在。
原典仕様なので忠実に再現し、roxygen2 の Details に
「local phase で向きが反転・追加された辺の重みは 1 に設定される（原典仕様）」と注記する。

## 3. アルゴリズム（原典 93-399 行目の精読に基づく実装ガイド）

内部は Python と同じ向き（`W[i,j] = i→j`, `M = X %*% W`）で計算し、最後に転置する。
n = nrow(X), d = ncol(X)。

### 3.1 mixed 損失と勾配（原典 113-131 行目）

`dis_con` 相当として `con <- as.numeric(is_continuous)`（1=連続, 0=離散）を行ベクトル的に使う。

```
M <- X %*% W
R <- X - M
# 損失: 離散列は logistic、連続列は log-cosh（Laplace）。a1 = a2 = 1
loss_dis <- sum((log1pexp(M) - X * M) %*% diag(1 - con))   # 列マスク。実装は sweep か rep でよい
loss_con <- sum(-logcosh(R) %*% diag(con))
loss <- -(loss_dis + loss_con) / n
```

- `log1pexp(M)` = `np.logaddexp(0, M)` = log(1 + exp(M))。オーバーフロー対策込みで実装
  （M > 33 で M をそのまま返す等。R では `ifelse(M > 33, M, log1p(exp(M)))` 程度で可）
- `logcosh(r)` は `abs(r) + log1p(exp(-2*abs(r))) - log(2)` で安定化（数値同等）

勾配（原典 121-131 行目）。原典は離散/連続の**行・列マスク行列**を作る:

```
W_dis: 離散変数 ii について W_dis[ii, ] <- 1 かつ W_dis[, ii] <- 1（他は 0）
W_con: 連続変数 jj について同様
G_dis <- (t(X) %*% (plogis(M) - X)) * W_dis      # 要素ごとの積
G_con <- -(t(X) %*% tanh(R)) * W_con
G_loss <- (G_dis + G_con) / n
```

注意: 原典 123-128 行目の `for ii in np.where(...)` は numpy の tuple をイテレートするため
実質「全離散インデックスをまとめて1回で代入」している。R では
`W_dis[!is_continuous, ] <- 1; W_dis[, !is_continuous] <- 1` が等価。
離散と連続の両方に接する要素は G_dis と G_con の**和**になる（原典と同じ挙動にすること）。

### 3.2 非巡回制約 `h(W)`（原典 244-252 行目）

```
M_h <- diag(d) + W * W / d              # 要素ごとの積
E   <- matpow(M_h, d - 1)               # 行列べき乗（反復乗算で自前実装, d-1 乗）
h   <- sum(t(E) * M_h) - d
```

勾配は `_func` 内で `(rho * h + alpha) * t(E) * W * 2` として合成される（原典 265 行目）。
`_h` 自体は `E`（= G_h）を返す設計にする。

### 3.3 拡張ラグランジュの目的関数（原典 254-269 行目）

パラメータ w は長さ 2*d*d（非負分解 w = (w+, w-)、`W = matrix(w[1:d^2] - w[(d^2+1):(2*d^2)], d, d)`。
**numpy の reshape は行優先**なので R では `matrix(..., d, d, byrow = TRUE)` か、
全体を転置規約で統一するか、どちらかに決めて一貫させること。テストで向きを検証する）。

```
obj      <- loss + 0.5 * rho * h^2 + alpha * h + lambda1 * sum(w)
G_smooth <- G_loss + (rho * h + alpha) * t(E) * W * 2
g_obj    <- c(G_smooth + lambda1, -G_smooth + lambda1)    # w+ 側と w- 側
```

- `optim()` の境界: 全要素 `lower = 0`, `upper = Inf`、ただし **W の対角に対応する
  要素（w+ 側・w- 側の両方）は upper = 0**（原典 286-291 行目の `(0,0) if i==j`）
- scipy は `jac=True` で obj と勾配を同時計算する。R の `optim` は `fn`/`gr` を別々に
  呼ぶので、まず素直に2回計算する実装でよい（正しさ優先。遅ければ後でキャッシュ化）

### 3.4 外側ループ（原典 279-321 行目）

```
w_est <- runif(2 * d * d);  rho <- 1.0;  alpha <- 0.0;  h <- Inf
for (iter in 1:max_iter) {
  while (rho < rho_max) {
    sol   <- optim(w_est, fn, gr, method = "L-BFGS-B", lower = ..., upper = ...)
    w_new <- sol$par
    h_new <- h(adj(w_new))
    if (h_new >= 0.25 * h) rho <- rho * 10 else break
  }
  w_est <- w_new;  h <- h_new
  alpha <- alpha + rho * h
  if (h <= h_tol && h != 0) break
  if (rho >= rho_max * 1e-6 && h > 1e5) {        # full graph 回避
    w_est <- runif(2 * d * d);  rho <- 1.0
  } else if (rho >= rho_max) break
  if (sum(abs(adj(w_est))) < 0.09) w_est <- runif(2 * d * d)   # zero matrix 回避
}
W_est <- adj(w_est)
W_est[abs(W_est) < w_threshold] <- 0
```

- 乱数はユーザーの `set.seed()` に従う。内部でシードを固定しない
- **原典 322 行目に `print("W_est (without the 2nd phase) is: ...")` がある。
  R 版では出力しない**（cat/print/message いずれも入れない）

### 3.5 ローカル探索（原典 324-397 行目。`only_global = FALSE` のとき）

損失関数を「mixed_dag」= **DAG 骨格の BIC スコア（符号反転）** に切り替える。
原典はここで `self._loss_type` を破壊的に書き換えるが（325行目、再 fit で挙動が変わる
Python 側の潜在バグ）、R は関数内ローカルに loss 関数を切り替えるだけでよい。
**mixed_dag の勾配（原典 141-142 行目）はローカル探索で一切使われないので実装しない**
（loss 値のみ返す関数にする）。

ローカル探索は3つのサブフェーズからなる（要約ではなく原典の正確な構造）:

**(a) 向き反転の全列挙**（330-342行目）:
- `I <- which(W_est != 0, arr.ind = TRUE)`（非ゼロ位置 = 辺集合。|E| 本）
- 各辺について {0=そのまま, 1=反転} の全組合せ 2^|E| 通り（`expand.grid` 等）を、
  **恒等設定（全部0）を除いて**評価（原典はインデックス1から）
- 反転: `W_tmp[from, to] <- 0; W_tmp[to, from] <- 1`（**係数でなく 1 を置く**）
- `loss < aa`（現行最良）かつ `h(W_tmp) < h_tol` なら採用。aa の初期値は W_est の mixed_dag loss
- **計算量ガード（R 独自追加）**: |E| > 15 なら
  `warning("Local search skipped: too many edges (2^|E| direction combinations); returning the global solution.", call. = FALSE)`
  で global 解を返す（Python にはガードが無い）

**(b) 刈り込み**（344-366行目、`d > 2 && |E| > d-1` のとき）:
- W_min_lss の各非ゼロ要素を1つずつ 0 にして loss 改善 + 非巡回なら採用
- さらに、(a)(b) の結果の骨格が W_est と異なる場合、**W_est を基準にした削除も
  もう1周試す**（356-366行目。見落としやすい）

**(c) 辺の追加**（367-392行目、`d > 2 && |E| < d(d-1)/2` のとき）:
- `W_edges <- W0 + t(W0) + diag(d)` が 0 の位置（どちら向きの辺も無いペア）へ
  1 を置いて loss 改善 + 非巡回なら採用
- こちらも W_est 基準の追加をもう1周試す（380-392行目）

### 3.6 BIC スコアリング（原典 147-242 行目、Poisson 分岐は除く）

W（0/1 骨格として解釈。W[i,j]≠0 ⇒ 辺 i→j、親 pa(j) = {i : W[i,j]≠0}）に対し:

```
penalty     <- log(n) / 2 * (辺数 + d)
total_score <- -penalty + Σ_i (変数 i の局所対数尤度)
mixed_dag_loss <- -total_score
```

各変数 i の局所対数尤度:

- **離散・親なし**（160-168行目）: Bernoulli MLE。`tab <- table(x_i)` として
  `Σ_k tab_k * (log(tab_k) - log(n))`
- **離散・親あり**（170-184行目）: ロジスティック回帰の対数尤度。
  R では `fit <- glm(x_i ~ X[, pa], family = binomial)` → `as.numeric(logLik(fit))`。
  **既知の差異**: sklearn の `LogisticRegression()` は既定で L2 正則化（C=1.0）が
  かかっており厳密な MLE ではない。R の glm は無正則化 MLE なので数値は完全一致しない。
  ここは glm を採用し、roxygen2 の Details に差異を1行記載する
- **連続**（224-241行目）: `likelihood_i`（utils 866-895行目）を移植する。
  glm/lm の logLik では**ない**。以下の super-Gaussian 尤度:
  ```
  親あり: lm(x_i ~ X[, pa]) で係数 b と切片 b0 を推定
  親なし: b <- 0, b0 <- mean(x_i)
  resid_std <- (x_i - X %*% b_full - b0) / sqrt(var_i)
      # var_i = mean((e - mean(e))^2)  ただし e = x_i - X %*% b_full （切片は引かない。
      #         np.var は平均を引き n で割る。原典 variance_i 916-939行目）
  ll <- sum(-sqrt(2) * abs(resid_std) - 0.35) - n * log(sqrt(var_i))
  ```
  定数 -0.35 は正規化定数（utils 912行目）。b_full は長さ d のベクトルで親以外は 0

networkx の代替: 親集合の取得は W の列の非ゼロで直接得られる。DAG 判定は
`h(W) < h_tol` をそのまま使う（原典と同じ）。

## 4. `only_global = TRUE` のとき

W_est（閾値処理済み）をそのまま最終結果とする（原典 396-397 行目）。

## 5. データ生成関数 `generate_lim_sample()`

`R/generate_lingam_sample.r` のスタイルを踏襲した固定構造の生成関数。
チュートリアルの生成コードはパッケージ外のローカル utils 由来で入手不可のため、
LiM のモデル定義から設計する:

- 連続変数: 親の線形結合 + Laplace ノイズ（`rexp(n) - rexp(n)` で生成可。依存追加なし）
- 離散変数: 線形予測子 eta に対し `rbinom(n, 1, plogis(eta))`
- シグネチャ例: `generate_lim_sample(n = 1000, seed = NULL)`
- 構造は 3 変数程度の既知 DAG（例: x1(連続) → x2(離散) → x3(連続)、係数 ±1 前後で固定。
  チュートリアルは係数 1.3 の 2 変数で n=1000 → 回復可能な強さの係数にする）
- 返り値: `list(data, adjacency_matrix（真の B, lingamr 規約）, is_continuous)`
  （既存 generate_* の返り値形式を確認して合わせる）

## 6. テスト（tests/testthat/test-lingam_lim.R）

test-lingam_direct.R のパターンを踏襲。**数値スナップショットは使わない**（乱数初期化依存）。

1. 構造: `expect_s3_class(res, "LiMResult")`, `expect_named()`, 次元、dimnames
2. バリデーション: 非数値 X、NA、`is_continuous` 長さ不一致・非 logical、
   離散指定列が2値でない、それぞれ `expect_error()`
3. print: `expect_output(print(res), "LiM Result")`
4. 再現性: `set.seed(1)` で2回 fit して `expect_equal(r1, r2)`
5. **向きの検証（最重要）**: 既知構造 x_from → x_to（連続→離散など）で
   `adjacency_matrix[to, from] != 0` かつ `adjacency_matrix[from, to] == 0`。
   転置ミスを検出するテスト。n=1000、`set.seed` 固定。乱数初期化で不安定な場合は
   複数シードで多数決にするか `only_global = TRUE` で安定する設定を探す
6. `only_global = TRUE` でも valid な結果
7. `w_threshold = 10` など極端な値で辺が全て消えること（+ causal_order が全変数を含むこと）
8. causal_order と B の非ゼロパターンの整合（トポロジカル順序と矛盾しない）
9. **標準出力が無いこと**: `expect_silent(lingam_lim(...))`（原典の print 混入防止。
   warning が出るケースは除いて設計）

実行: `devtools::test(filter = "lingam_lim")` → 全体 `devtools::test()`

## 7. ドキュメント

- roxygen2: lingam_direct.r の書式に倣う。`@param`, `@return`, `@references`
  （Zeng et al. CLeaR2022）, `@examples`
- Details に必ず記載する4点:
  1. 離散変数は2値(0/1)のみ対応
  2. Python 版 `adjacency_matrix_` とは転置の関係（本実装は lingamr 規約 B[i,j]=j→i）
  3. local phase で反転・追加された辺の重みは 1 になる（原典仕様）
  4. 離散変数の尤度は R の glm（無正則化 MLE）を使うため、sklearn の L2 付き
     LogisticRegression を使う Python 版と数値が完全一致しない
- `@examples` は小さな n（300程度）+ 2〜3変数。遅ければ `\donttest{}`
- 乱数初期化に依存するため「結果の再現には `set.seed()` を使う」ことを examples で示す
- `devtools::document()` で man/NAMESPACE 再生成
- NEWS.md: 開発版セクションに追記
- vignette: 新セクション「LiNGAM for Mixed Data (LiM)」（モデル説明 → 生成 → 実行 → 解釈）
- _pkgdown.yml: reference に `lingam_lim`, `generate_lim_sample` を追加

## 8. 完了条件・検証手順

R は PATH 外の `C:\R\R-4.6.0` にある。実行例:

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] `devtools::document()` が警告なく通る
- [ ] `devtools::test()` 全緑（既存テストを壊さない）
- [ ] `devtools::check()` で新規 ERROR/WARNING/NOTE なし
- [ ] 向きの検証テスト（セクション6-5）が通る
- [ ] `expect_silent` テスト（セクション6-9）が通る
- [ ] vignette がビルドできる

## 9. 遵守事項

- 既存ファイル（NEWS.md, vignettes/lingamr.Rmd, _pkgdown.yml）の編集は追記のみ
- テストコードの削除・コメントアウト禁止
- commit / push はしない（ユーザーが行う）
- エラーは握りつぶさない。意味のあるメッセージ付きで stop/warning
- 範囲外のリファクタリング禁止
