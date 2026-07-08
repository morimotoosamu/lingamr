# HighDimDirectLiNGAM 実装指示書

lingamr パッケージに HighDimDirectLiNGAM（高次元 Direct LiNGAM）を追加する。
Python cdt15/lingam の `HighDimDirectLiNGAM` クラス（`high_dim_direct_lingam.py`、全242行）の R 移植。

アルゴリズムの本体は「モーメントベースの統計量 τ による因果順序探索」であり、
乱数を使わない決定的アルゴリズム。隣接行列推定は既存の lingamr 基盤を再利用できるため、
新規実装の中心は `fit()` の while ループと τ 統計量の計算である。

## 0. 参照資料（実装前に必ず全文を読むこと）

本指示書のアルゴリズム記述は原典の全文（242行）を読んで書かれているが、
**実装時も必ず原典を WebFetch で取得し、原典を正とする**こと。

- ソース: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/high_dim_direct_lingam.py
- 親クラス（n > p 時の隣接行列推定 `_estimate_adjacency_matrix` の確認用）:
  https://raw.githubusercontent.com/cdt15/lingam/master/lingam/direct_lingam.py
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/high_dim_direct_lingam.html
- 論文: Y. S. Wang and M. Drton. High-dimensional causal discovery under non-Gaussianity.
  Biometrika 107(1): 41-59, 2020.

また、着手前に必ずリポジトリの CLAUDE.md と以下の既存実装を読み、規約を踏襲すること:

- `R/lingam_direct.r` — 関数シグネチャ、バリデーション、S3 クラス（LingamResult）、roxygen2 の手本
- `R/fit_regression.r` — **`estimate_adjacency_matrix()` と `fit_adaptive_lasso()` を再利用する**（後述）
- `R/generate_lingam_sample.r` — テストデータ生成関数（既存を再利用。新規生成関数は作らない）
- `tests/testthat/test-lingam_direct.R` — テストパターンの手本

## 1. スコープ

### 実装するもの

| 成果物 | ファイル |
|---|---|
| `lingam_high_dim()` | `R/lingam_high_dim.r`（新規） |
| testthat テスト | `tests/testthat/test-lingam_high_dim.R`（新規） |
| roxygen2 → man/NAMESPACE 再生成 | `devtools::document()` |
| NEWS.md 追記 | `NEWS.md` |
| vignette 追記（HighDim セクション） | `vignettes/lingamr.Rmd` |
| pkgdown reference 追記 | `_pkgdown.yml` |

### 実装しないもの

- 新規の結果クラス。**返り値は既存の `LingamResult` をそのまま使う**
  （フィールドが `adjacency_matrix` + `causal_order` で完全一致し、
  print / tidy / plot_adjacency / estimate_total_effect 等の既存 S3 資産が全て使えるため）
- 新規のデータ生成関数（既存の `generate_lingam_sample_*` を再利用。
  p > n ケースのデータはテスト内でインラインに構築する）
- `random_state` 引数（Python 版でも保持するだけで未使用。アルゴリズムは決定的）
- 新規パッケージ依存。glmnet は既存の Suggests + `check_glmnet_available()` パターンで利用

## 2. API 設計

```r
lingam_high_dim <- function(X,
                            J = 3L,
                            K = 4L,
                            alpha = 0.5,
                            estimate_adj_mat = TRUE)
```

- `X`: 数値行列または data.frame。列名を保存する（lingam_direct と同様）
- `J`: 想定する最大入次数（Python: `check_scalar(J, min_val=2, include_boundaries="neither")`
  → **J >= 3 の整数**。J <= 2 はエラー）
- `K`: 非ガウス性を測るモーメントの次数（**K >= 1 の整数**）
- `alpha`: 偽の親を刈るカットオフ係数（**[0, 1] の数値**）
- `estimate_adj_mat`: FALSE なら隣接行列推定をスキップし因果順序のみ返す

### バリデーション（lingam_direct.r:73-93 のスタイルを踏襲）

- `X` の数値・NA・次元チェック（`ncol >= 2`, `nrow >= 2`）
- `J`, `K`, `alpha`, `estimate_adj_mat` の型・範囲チェック（上記の Python 準拠の範囲で
  `stop(call. = FALSE)`、メッセージは英語）

### 返り値: 既存クラス `LingamResult`

```r
result <- list(adjacency_matrix = B, causal_order = K_order)
class(result) <- "LingamResult"
```

- `adjacency_matrix` の dimnames に X の列名を付ける（既存実装の挙動に合わせる）
- `estimate_adj_mat = FALSE` のときは `adjacency_matrix` を **NA 埋めの p×p 行列**
  （dimnames は付ける）にする。NULL にしないこと（`print.LingamResult` や
  `validate_lingam_result()` を壊さないため）。roxygen2 に「FALSE 時は NA 行列」と明記

### 隣接行列の向き

Python 原典は `B[target, predictors] <- coef`（`_estimate_adjacency_matrix2`, 原典199行目）で、
**lingamr 規約 `B[i,j] = j → i`（行=to, 列=from）とそのまま一致する。転置は不要**。
念のため単体テストで既知構造により検証する（セクション5）。

## 3. アルゴリズム（原典を正とする。以下は精読済みの実装ガイド）

### 3.1 fit() のメインループ（原典 49-114 行目）

p = ncol(X), n = nrow(X)。**インデックスは Python が 0-based、R は 1-based** —
psi/theta/cond_set の集合演算を移植する際に最も間違えやすい点なので注意。

```
Y     <- X                          # check_array 相当のバリデーション後
yty   <- t(Y) %*% Y                 # Gram 行列をキャッシュ（毎回再計算しない）
cut_off <- 0
theta <- integer(0)                 # 因果順序（先頭=最上流）
psi   <- 1:p                        # 候補集合
prune_stats <- matrix(1e5, p, p); diag(prune_stats) <- 0

while (length(psi) > 1) {
  new_stats <- 各 v in psi について:
    cond_set  <- intersect(theta, which(prune_stats[v, ] > cut_off))
    cond_set  <- union(cond_set, tail(theta, 1))       # theta が空なら union しない
    last_root <- if (length(theta) > 0) tail(theta, 1) else NULL
    get_prune_stats(v, psi, K, last_root, cond_set, J)  # 長さ p のベクトル
  （new_stats は length(psi) × p 行列に rbind）

  prune_stats[psi, ] <- pmin(prune_stats[psi, ], new_stats)   # 要素ごとの min
  diag(prune_stats) <- 0

  max_taus <- apply(prune_stats[psi, psi, drop = FALSE], 1, max)
  r        <- psi[which.min(max_taus)]
  cut_off  <- max(cut_off, min(max_taus) * alpha)

  theta <- c(theta, r)
  psi   <- setdiff(psi, r)
}
causal_order <- c(theta, psi)       # 残った1変数を末尾に
```

- ルート選択は「行方向の最大 τ が最小の変数」（原典 93-95 行目）。which.min のタイ処理は
  R の「最初の最小値」で Python の argmin と一致する
- `prune_stats[psi, ] <- pmin(...)` の行対応: new_stats の行順は psi のループ順と同一

### 3.2 τ 統計量

**無条件版 `calc_tau(k, pa, ch)`**（原典 116-117 行目、theta が空の初回イテレーションで使用）:

```
tau = | mean(pa^(k-1) * ch) * mean(pa^2) - mean(pa^k) * mean(pa * ch) |
```

**条件付き版 `calc_taus(pa, ch, k, cond_sets, an_sets)`**（原典 119-144 行目）:

各条件集合 cond について
1. Gram 行列から回帰係数: `b <- pinv(yty[cond, cond]) %*% yty[cond, pa]`
2. 残差: `resid <- Y[, pa] - Y[, cond] %*% b`
3. `resid_k_1 <- resid^(k-1)`, `resid_var <- mean(resid^2)`, `resid_k <- mean(resid^k)`
4. 各 j in ch について
   `value = (1/n) * sum(resid_k_1 * Y[,j]) * resid_var - resid_k * (1/n) * sum(resid * Y[,j])`
   `ret[j] <- min(ret[j], abs(value))`（ret は 1e10 で初期化した長さ p のベクトル）
5. 各 an_sets[z, i]（条件集合から除外した祖先候補）についても同じ式で更新

**【上流バグ疑いに関する2026-07-06追記の方針転換】** 原典 144 行目の `return ret` は
`for z in range(len(cond_sets)):` ループの**内側**にインデントされており、
**最初の条件集合（z=0、Python の itertools.combinations の辞書順で先頭）だけを評価して
即 return する**。全条件集合を回すコードには見えるが、実際には回らない。

初回実装（2026-07-06）では「Python との数値一致を最優先し、この挙動を忠実に再現する」
方針を採ったが、レビュー指摘を受けてパッケージ作者が再検討した結果、
**理論的正確性（Wang & Drton 2020 のアルゴリズム通り、全条件集合を評価して最小値を取る）を
優先する方針に変更した**。`R/lingam_high_dim.r` の `calc_taus()` からは
ループ内の `return(ret)` を削除済み。この結果、`lingam_high_dim()` の
`causal_order`/`adjacency_matrix` は上流 Python 版（cdt15/lingam）とは
数値的に一致しなくなる。この点は `?lingam_high_dim` の `@details` に明記済み。

（旧方針時に残していた「忠実に再現する」旨のコメントは削除し、代わりに
上流のバグ内容と、あえて再現しない理由を roxygen コメントに明記した。）

**pinv（擬似逆行列）**: `np.linalg.pinv` 相当を base R で自前実装する（新規依存を増やさない。
MASS::ginv は使わない）。svd ベースで数行:

```r
pinv <- function(A, tol = max(dim(A)) * .Machine$double.eps) {
  s <- svd(A)
  pos <- s$d > tol * max(s$d)
  if (!any(pos)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}
```

（numpy の既定 rcond に相当する閾値。`@keywords internal`）

### 3.3 条件集合の列挙 `get_prune_stats(i, j, K, last_root, condition_set, J)`（原典 146-173 行目）

- `j <- setdiff(j, i)`（自分自身を除く）
- `last_root` が NULL（初回）: `prune_stat[j] <- calc_tau(K, Y[,i], Y[,j_])` を各 j_ に。
  それ以外の要素は 1e5 のまま返す
- それ以外:
  - `size_of_set <- min(J, length(condition_set))`
  - `rest <- setdiff(condition_set, last_root)`
  - `rest` の中から **size_of_set - 1 個**を選ぶ全組合せを列挙
    （原典は `itertools.combinations` — R では `utils::combn`。
    `combn(rest, m)` は rest がスカラーのとき `seq_len(rest)` に化ける罠があるので、
    `length(rest) == 1` と `m == 0`（空組合せ1個）のエッジケースを明示的に処理すること。
    原典 157-158 行目にも `len == 1` の特別扱いがある）
  - 各組合せ x について `an_set <- setdiff(condition_set, x)`（除外された変数 = 祖先候補）
  - 各組合せの先頭に `last_root` を連結して条件集合とする（原典 167-170 行目）
  - `calc_taus(i, j, K, condition_sub_set, an_sets)` を呼ぶ
- 原典 164-165 行目の `an_sets` の shape 補正（組合せが1つのときの転置）は numpy の
  次元潰れ対策。R では組合せを常に「リスト of integer vector」または
  `drop = FALSE` 付き行列で持ち、次元潰れ自体を起こさない設計にすること

### 3.4 隣接行列推定（原典 104-114 行目）

因果順序確定後、`estimate_adj_mat = TRUE` のとき:

- **n > p の場合**: 既存の **`estimate_adjacency_matrix()`（R/fit_regression.r:65）を再利用**する。
  Python の親クラス DirectLiNGAM の `_estimate_adjacency_matrix`（adaptive lasso）に相当する
  lingamr 側の実装が既にあるので、新規実装しない。呼び出しシグネチャ
  （reg_method, lambda 等のデフォルト）は fit_regression.r を読んで確認し、
  lingam_direct のデフォルト（adaptive_lasso / BIC）と同じ設定で呼ぶ
- **n <= p の場合**: `warning()` を出し（原典 108-111 行目のメッセージに準拠:
  "Since n_samples <= n_features, the adjacency matrix is estimated with
  cross-validated lasso instead of BIC-based lambda selection." 等）、
  原典 `_predict_adaptive_lasso`（204-242 行目）相当を実装する:
  1. X 全体を標準化（scale()。Python の StandardScaler は母標準偏差＝n で割る。
     lingamr には `sd_pop()` があるので合わせること）
  2. 標準化データで OLS: predictors → target
  3. 重み `w <- abs(coef_ols)^1`（gamma = 1.0 固定）
  4. `X_std[, predictors] * w`（列ごとのスケーリング）に対し **交差検証で λ を選ぶ lasso**
     を当てる。Python は LassoLarsCV — R では `glmnet::cv.glmnet(alpha = 1)` で代替する
     （LARS と座標降下の違いで数値は完全一致しないが、CV 選択という設計意図を再現する。
     この差異を roxygen2 の Details に一言記載する）
  5. `abs(coef * w) > 0` の変数だけ残し、**元スケールの X** で OLS を再フィット
  6. `B[target, predictors[pruned]] <- coef`
  - glmnet は `check_glmnet_available()`（R/fit_regression.r:34）で条件付きロード。
    n <= p かつ glmnet 不在なら意味のあるエラーで stop
  - cv.glmnet の CV 分割は乱数依存。**再現性テストでは set.seed を fit 呼び出し直前に置く**こと

## 4. 性能に関する注意

- `yty`（Gram 行列）は fit の最初に1回だけ計算し、`calc_taus` の回帰はすべて
  `yty` の部分行列から解く（原典と同じ。n が大きくても回帰コストが n に依存しない設計）
- `resid` の計算（原典 128 行目）だけは Y の実データを使う。ここは避けられない
- while ループは p-1 回、各回で |psi| 個の変数 × 条件集合の組合せ。
  cut_off による枝刈りで cond_set が絞られる設計なので、まず原典どおり素直に実装し、
  プロファイル前の最適化はしない
- `resid^k`（K=4 で4乗）はスケールの大きいデータでオーバーフローし得る。
  原典は対策していないので R でも対策不要だが、テストデータは標準的スケールにする

## 5. テスト（tests/testthat/test-lingam_high_dim.R）

test-lingam_direct.R のパターンを踏襲。数値スナップショットは使わない。

1. 構造: `expect_s3_class(res, "LingamResult")`, `expect_named()`, 行列次元、dimnames が列名
2. バリデーション: 非数値 X、NA 入り X、`J = 2`（境界外）、`K = 0`、`alpha = 1.5`、
   `estimate_adj_mat = "yes"` それぞれ `expect_error()`
3. **因果順序と向きの検証**: `generate_lingam_sample_6()` 等の既存生成関数
   （非ガウス誤差であることを確認して選ぶ）で n=1000 程度のデータを作り、
   - 真の DAG と整合する因果順序が返ること（完全一致ではなく「真の辺 j→i について
     j が i より先に並ぶ」ことの検証にすると頑健）
   - 真の辺位置 `adjacency_matrix[to, from]` が非ゼロで係数符号が真値と一致すること
4. 決定性: 同じ X で2回 fit して `expect_equal(r1, r2)`（n > p の経路。
   アルゴリズム本体は乱数を使わないことの確認）
5. `estimate_adj_mat = FALSE`: `adjacency_matrix` が全 NA の p×p 行列で、
   `causal_order` は TRUE のときと同一
6. **n <= p の経路**: 小さめの p > n データ（例: 真の DAG から n=40, p=45 を生成、
   または既存6変数構造を横に複製）で
   - `expect_warning(..., "n_samples <= n_features")` 相当の警告文を検証
   - 結果が valid な LingamResult であること
   - `skip_if_not_installed("glmnet")` を付ける
   - 実行が遅い場合は `skip_on_cran()` を付ける
7. glmnet 不在時: `local_mocked_bindings()` で `check_glmnet_available` を差し替え
   （test-lingam_direct.R:195-212 の既存パターンを流用）、n <= p でエラー、
   n > p の OLS 系経路は動くことを確認

実行: `devtools::test(filter = "lingam_high_dim")` → 全体 `devtools::test()`

## 6. ドキュメント

- roxygen2: lingam_direct.r の書式に倣う。`@param`, `@return`（箇条書き）,
  `@references`（Wang & Drton, Biometrika 2020）, `@examples`
  - Details に「高次元（p が大きい / p > n）向けの Direct LiNGAM 変種。因果順序探索が
    モーメント統計量ベースで高速」「n <= p では CV lasso（glmnet）で隣接行列を推定」
    「Python 版は LassoLarsCV、本実装は cv.glmnet のため数値は完全一致しない」を記載
  - `@examples` は n=300・6変数程度で直接実行、n<=p の例は `\donttest{}` +
    glmnet の存在チェック（`if (requireNamespace("glmnet", quietly = TRUE))`）で包む
- `@export` は `lingam_high_dim` のみ。pinv 等の内部関数は `@keywords internal`
- `devtools::document()` で man/NAMESPACE を再生成
- NEWS.md: 開発版セクションに「Added `lingam_high_dim()` ...」を追記
- vignette（vignettes/lingamr.Rmd）: 新セクション「High-Dimensional Direct LiNGAM」を追加。
  構成: モデル説明（どういう時に lingam_direct でなくこちらを使うか）→ 使用例 →
  結果解釈。返り値が LingamResult なので既存の plot_adjacency / estimate_total_effect が
  そのまま使えることを1文添える。チャンク実行時間に注意
- _pkgdown.yml: reference の lingam_direct と同じグループに `lingam_high_dim` を追加

## 7. 完了条件・検証手順

R は PATH 外の `C:\R\R-4.6.0` にある。実行例:

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] `devtools::document()` が警告なく通る
- [ ] `devtools::test()` 全緑（既存テストを壊さない）
- [ ] `devtools::check()` で新規 ERROR/WARNING/NOTE なし
- [ ] 因果順序・向きの検証テスト（セクション5-3）が通る
- [ ] n <= p 経路の警告テスト（セクション5-6）が通る
- [ ] vignette がビルドできる

## 8. 遵守事項

- 既存ファイル（NEWS.md, vignettes/lingamr.Rmd, _pkgdown.yml）の編集は追記のみ。
  既存記述を書き換えない。**R/fit_regression.r 等の既存関数は変更しない**（再利用のみ）
- テストコードの削除・コメントアウト禁止
- commit / push はしない（ユーザーが行う）
- エラーは握りつぶさない。意味のあるメッセージ付きで stop/warning
- 範囲外のリファクタリング禁止
