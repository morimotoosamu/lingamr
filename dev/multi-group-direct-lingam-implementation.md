# MultiGroupDirectLiNGAM 実装指示書

lingamr パッケージに MultiGroupDirectLiNGAM（複数データセットの同時推定）を追加する。
Python cdt15/lingam の `MultiGroupDirectLiNGAM` クラス（`multi_group_direct_lingam.py`、全329行）
の R 移植。**本指示書は原典329行と、再利用対象である lingamr 側の
`R/lingam_direct.r` / `R/search_causal_order.r` / `R/lingam_bootstrap.r` を精読した上で書かれている。**

このアルゴリズムは「複数のデータセット（群）が**共通の因果順序**を持つ」と仮定し、
因果順序探索の目的関数をサンプルサイズ加重で群横断に合算する Direct LiNGAM の拡張
（Shimizu 2012）。**新規に書くコードの大半は薄いラッパで、lingamr の既存内部関数を
最大限再利用すること**が本指示書の主旨。

## 0. 参照資料（実装前に必ず全文を読むこと）

- ソース: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/multi_group_direct_lingam.py
- 親クラス（継承メソッドの確認用）:
  https://raw.githubusercontent.com/cdt15/lingam/master/lingam/direct_lingam.py
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/multiple_dataset.html
- 論文: S. Shimizu. Joint estimation of linear non-Gaussian acyclic models.
  Neurocomputing, 81: 104-107, 2012.

着手前に必ず読む lingamr 側ファイル（本指示書はこれらの具体的な行を参照する）:

- `R/lingam_direct.r` — fit ループの手本（119-137行目）とバリデーション様式
- `R/search_causal_order.r` — `search_candidate` / `search_causal_order_pwling` /
  `residual_vec` / `entropy_approx` / `extract_partial_orders`（すべて再利用する）
- `R/fit_regression.r` — `estimate_adjacency_matrix()`（群ごとに再利用）
- `R/lingam_bootstrap.r` — bootstrap の並列・再現性・失敗処理パターン（踏襲する）
- `R/paths.r` — パスベースの総合効果計算（bootstrap 内で再利用）
- `R/estimate_total_effect.r`, `R/get_error_independence_p_values.r` — 群抽出後に
  そのまま使えることを確認する
- `tests/testthat/test-lingam_direct.R`, `test-bootstrap.R` — テストパターンの手本

## 1. スコープ

### 実装するもの

| 成果物 | ファイル |
|---|---|
| `lingam_multi_group()` + `MultiGroupLingamResult` + print | `R/lingam_multi_group.r`（新規） |
| `get_group_result()`（群→LingamResult 抽出ヘルパ） | 同上 |
| `lingam_multi_group_bootstrap()` + print | `R/lingam_multi_group_bootstrap.r`（新規） |
| `generate_multi_group_sample()` | `R/generate_multi_group_sample.r`（新規） |
| testthat テスト | `tests/testthat/test-lingam_multi_group.R`（新規） |
| roxygen2 → man/NAMESPACE 再生成 | `devtools::document()` |
| NEWS.md 追記 | `NEWS.md` |
| vignette 追記 | `vignettes/lingamr.Rmd` |
| pkgdown reference 追記 | `_pkgdown.yml` |

### 実装しないもの（設計判断込み）

- `measure = "kernel"`: Python の MultiGroup 探索はエントロピー近似（pwling 相当）のみ。
  `measure` 引数自体を設けない。roxygen2 の Details に記載
- Python の `estimate_total_effect` / `get_error_independence_p_values` に相当する
  **マルチグループ専用ラッパ関数は作らない**。原典（160-270行目）はどちらも
  「群ごとにループして単群処理を適用するだけ」なので、R では
  `get_group_result()` で群を `LingamResult` として取り出し、既存の
  `estimate_total_effect()` / `estimate_all_total_effects()` /
  `get_error_independence_p_values()` をそのまま使う設計とする
  （R らしい合成可能な API。vignette でこの使い方を必ず示す）
- 新規パッケージ依存なし

## 2. API 設計

### 2.1 本体

```r
lingam_multi_group <- function(X_list,
                               prior_knowledge = NULL,
                               apply_prior_knowledge_softly = FALSE,
                               reg_method = "adaptive_lasso",
                               lambda = "BIC",
                               init_method = "ols")
```

- `X_list`: 数値行列または data.frame の**リスト**（各要素 n_d × p。n_d は群ごとに
  異なってよい。p は全群で同一）
- `prior_knowledge` 以下は `lingam_direct()` と同じ意味・同じ既定値
  （prior knowledge は**全群共通**に適用される。原典 70-78 行目と同じ）

**バリデーション**（原典 `_check_X_list` 272-290行目 + lingam_direct.r:73-93 の様式）:

- `X_list` がリストでない → エラー "X_list must be a list ..."
- `length(X_list) < 2` → エラー（原典: "at least two items"）
- 各要素に lingam_direct と同じチェック（numeric・NA なし・nrow >= 2）
- 全要素の `ncol` が一致しない → エラー
- 列名: 最初の群の列名を全体の変数名として採用。他の群に列名があり食い違う場合は
  `warning`（R 独自の親切。原典はチェックしない）
- 群名: `names(X_list)` があれば使い、無ければ `"group1", "group2", ...` を付与

### 2.2 返り値: S3 クラス `MultiGroupLingamResult`

```r
result <- list(
  adjacency_matrices = B_list,   # 群名付きリスト。各要素 p×p、lingamr 規約 B[i,j] = j -> i
  causal_order       = K         # 全群共通の因果順序（1-based）
)
class(result) <- "MultiGroupLingamResult"
```

- 各 B に dimnames（変数名）を付ける
- `print.MultiGroupLingamResult()`: ヘッダ "Multi-Group Direct LiNGAM Result"、
  群数・変数数・共通 causal order（`print.LingamResult` の書式 lingam_direct.r:162-176
  に倣う）、各群の隣接行列を群名見出し付きで順に表示。`invisible(x)`

### 2.3 群抽出ヘルパ（エクスポート）

```r
get_group_result <- function(x, group) {
  # x: MultiGroupLingamResult, group: 群名 or インデックス
  # -> list(adjacency_matrix = ..., causal_order = ...) に class "LingamResult" を付けて返す
}
```

- これにより既存の `estimate_total_effect()`, `estimate_all_total_effects()`,
  `get_error_independence_p_values()`, `plot_adjacency()`, `autoplot()`, `tidy()` が
  群単位でそのまま動く。roxygen2 の `@examples` と vignette で実演すること
- `group` の範囲外・存在しない群名は意味のあるエラー

### 2.4 隣接行列の向き

原典は群ごとに親クラスの `_estimate_adjacency_matrix` を呼ぶだけ（101-104行目）。
R でも群ごとに既存 `estimate_adjacency_matrix()` を呼ぶため、
**lingamr 規約 `B[i,j] = j→i` が自動的に成立する。転置等の変換は一切不要**。

## 3. アルゴリズム（原典 52-105, 292-316 行目の精読に基づく）

### 3.1 fit のメインループ

`lingam_direct.r` の 97-137 行目とほぼ同型。差分は「残差化を全群に適用」と
「探索がマルチグループ版」の2点だけ:

```r
# prior knowledge 前処理: lingam_direct.r:98-109 をそのまま流用
U <- seq_len(p); K <- integer(0)
X_list_ <- lapply(X_list, identity)   # 作業コピー
for (iter in seq_len(p)) {
  cand <- search_candidate(U, Aknw, apply_prior_knowledge_softly, partial_orders)  # 既存関数
  m <- search_causal_order_pwling_multi(X_list_, U, cand$Uc, cand$Vj)              # 新規（3.2）
  for (d in seq_along(X_list_)) {
    for (i in U) {
      if (i != m) X_list_[[d]][, i] <- residual_vec(X_list_[[d]][, i], X_list_[[d]][, m])  # 既存関数
    }
  }
  K <- c(K, m)
  U <- setdiff(U, m)
  # partial_orders の更新も lingam_direct.r:132-136 と同一
}
```

隣接行列は**元の（残差化していない）各群データ**で推定（原典 101-104 行目）:

```r
B_list <- lapply(X_list, function(Xd) {
  estimate_adjacency_matrix(Xd, K, Aknw, method = reg_method,
                            lambda = lambda, init_method = init_method)
})
```

### 3.2 マルチグループ因果順序探索 `search_causal_order_pwling_multi()`

原典 `_search_causal_order`（292-316行目）の数式:

```
total_size = Σ_d n_d
候補 i ごとに:
  MG_i = Σ_d (n_d / total_size) * M_{i,d}
  M_{i,d} = Σ_{j ∈ U, j≠i} min(0, diff_mutual_info(xi_std, xj_std, ri_j, rj_i))^2
選択: r = argmax_i(-MG_i)   （= MG_i 最小の候補）
```

- 標準化は群ごと・反復ごと（残差化後のデータに対して母標準偏差で。原典 309-310 行目の
  `np.std` は母標準偏差 = lingamr の `sd_pop` と同じ）
- `ri_j` の prior knowledge 分岐（`i ∈ Vj かつ j ∈ Uc` なら残差化スキップ）は
  原典 311-312 行目のとおりで、単群版と同一
- **実装方針: 既存 `search_causal_order_pwling`（search_causal_order.r:205-260）の
  最適化構造をそのまま群ループで包む**。すなわち群 d ごとに:
  1. U の列を一括標準化し `entropy_approx` を事前計算（211-218行目と同じ）
  2. 相関行列 `R <- crossprod(X_std[, U]) / n` から回帰係数と残差 SD
     `sqrt(1 - r^2)` を解析的に得る（220-241行目と同じ）
  3. diff_mutual_info の反対称性を使い i<j 側だけ計算して両者に加算
     （236-256行目のトリックと同じ）
  で群ごとの `M_acc_d`（長さ p）を求め、`M_total <- M_total + (n_d / total_size) * M_acc_d`
  と加重合算。最後に `Uc[which.max(-M_total[Uc])]`
- この関数は新規ファイル `R/lingam_multi_group.r` 内に `@keywords internal` で置く
  （**search_causal_order.r は変更しない**）
- 数値検証: 同一データを2群に複製して渡した場合、加重合算は単群の M_acc と
  一致するため、`lingam_direct()`（pwling, 同 reg_method）と causal_order が
  完全一致するはず。これをテストにする（セクション6）

### 3.3 バリデーションとの整合

`length(Uc) == 1` の早期リターン（原典 295-296 行目）は既存関数と同じ位置に置く。

## 4. bootstrap（原典 107-158 行目）

### 4.1 API

```r
lingam_multi_group_bootstrap <- function(X_list,
                                         n_sampling,
                                         prior_knowledge = NULL,
                                         apply_prior_knowledge_softly = FALSE,
                                         reg_method = "adaptive_lasso",
                                         lambda = "BIC",
                                         init_method = "ols",
                                         seed = NULL,
                                         verbose = TRUE,
                                         parallel = FALSE,
                                         n_cores = NULL,
                                         compute_total_effects = TRUE)
```

`R/lingam_bootstrap.r`（87-284行目）の構造・流儀を**そのまま踏襲**する:

- 引数バリデーションを cluster 起動前に実施（lingam_bootstrap.r:100-124 と同様）
- 1反復 = `run_one(i)`: **各群を独立に復元抽出**（原典 141 行目
  `[resample(X) for X in X_list]`）→ `lingam_multi_group()` で同時 fit →
  群ごとの隣接行列・共通 causal_order・総合効果を返す。`tryCatch` で失敗反復を
  警告付きスキップ（lingam_bootstrap.r:128-166 と同じ設計。全滅時のみ stop）
- 並列: `parallel::makePSOCKcluster` + `on.exit(stopCluster)` +
  `.libPaths` 伝搬 + `library` フォールバック + `clusterSetRNGStream(cl, seed)`
  （lingam_bootstrap.r:191-228 のコードパターンを踏襲。
  「同じ seed + 同じ n_cores で再現、逐次とは一致しない」という既存の
  再現性ポリシーも roxygen2 に同文で記載）
- verbose の進捗表示・経過時間表示も既存と同様

### 4.2 総合効果の計算方法 —【重要】既存 bootstrap との違い

原典のマルチグループ bootstrap は総合効果を**回帰でなくパス積**で計算する
（148-152行目: `estimate_total_effect2` → `calculate_total_effect(am, from, to)`。
隣接行列の係数のパス積和）。lingamr の `lingam_direct_bootstrap` が使う
回帰ベースの `estimate_all_total_effects` とは方式が異なる。

- **原典に忠実にパス積方式を採用する**。`R/paths.r` の既存内部関数
  （`find_all_paths` / パス効果の総和。正確な関数名とシグネチャは paths.r を読んで確認）
  を再利用し、causal_order 上で from より後の to のペア全てについて
  各群の隣接行列から総合効果を計算する
- roxygen2 の Details に「総合効果はパス積和（Python 版準拠）であり、
  `lingam_direct_bootstrap()` の回帰ベース推定とは方式が異なる」と明記

### 4.3 返り値

**群ごとの `BootstrapResult` の名前付きリスト**に class `"MultiGroupBootstrapResult"` を付ける:

```r
res <- lapply(seq_len(n_groups), function(d) {
  create_bootstrap_result(adjacency_matrices_d,  # n_success × p × p
                          total_effects_d,       # 同上（compute_total_effects=FALSE なら NULL）
                          resampled_indices_d,   # 群 d の復元抽出インデックスのリスト
                          causal_orders)         # 共通（全群同じ行列を入れてよい）
})
names(res) <- group_names
class(res) <- "MultiGroupBootstrapResult"
```

- `create_bootstrap_result()`（lingam_bootstrap.r:299-308、internal）を再利用
- 各要素が正規の `BootstrapResult` なので、既存の `get_probabilities()`,
  `get_causal_direction_counts()`, `get_directed_acyclic_graph_counts()`,
  `get_total_causal_effects()`, `get_causal_order_stability()`,
  `plot_bootstrap_probabilities()`, `tidy()` が **`result[[1]]` のように群を
  取り出すだけで全部使える**。これが Python（BootstrapResult のリストを返す、
  154-158行目）と同型の設計であることを Details に記載
- `print.MultiGroupBootstrapResult()`: 群数・sampling 数・各群1行のサマリ。
  `[[` で取り出せば `print.BootstrapResult` が働くので詳細表示はそちらに任せる
- 注意: class を付けても `[[` の挙動は変わらないのでそのまま動くが、
  `lapply` 等でクラスが落ちても実害がない設計にする（クラスは表示専用）

## 5. データ生成関数 `generate_multi_group_sample()`

チュートリアル準拠: **共通の因果構造・群ごとに異なる係数**の2群を生成する。

- シグネチャ例: `generate_multi_group_sample(n = c(1000, 1000), seed = NULL)`
- 構造は既存 `generate_lingam_sample_6()`（R/generate_lingam_sample.r）と同じ
  6変数 DAG を使い、係数だけ群2で少しずらす（チュートリアルは
  群1: 3.0, 6.0, 3.0, 2.0, 8.0, -1.0, 4.0 / 群2: 3.5, 6.5, 3.5, 2.5, 8.5, -1.5, 4.5
  という対応。既存生成関数の構造・誤差分布（非ガウス）を確認して合わせる）
- 返り値: `list(data_list = <2要素の named list>, adjacency_matrices = <真の B のリスト>,
  causal_order = <真の因果順序>)`（既存 generate_* の返り値形式に揃える）

## 6. テスト（tests/testthat/test-lingam_multi_group.R）

test-lingam_direct.R / test-bootstrap.R のパターンを踏襲。数値スナップショット不使用。

1. 構造: `expect_s3_class(res, "MultiGroupLingamResult")`, `expect_named()`,
   `adjacency_matrices` が群数分、各 p×p、dimnames
2. バリデーション: リストでない X_list、要素1個、列数不一致、NA 入り、非数値、
   それぞれ `expect_error()`
3. **単群一致**: 同一データを `list(X, X)` と2群にして渡すと、
   `lingam_direct(X, reg_method = "ols")` と causal_order が完全一致
   （3.2 の数理的に成り立つはず。成り立たない場合は実装ミス）
4. **回復**: `generate_multi_group_sample()`（n=1000×2, seed 固定）で
   共通 causal_order が真の順序と整合（真の辺 j→i で j が先行）し、
   各群の係数が各群の真値に近い（`tolerance` は既存テストの相場に合わせる）
5. `get_group_result()`: 返り値が `LingamResult` で、既存
   `estimate_all_total_effects()` と `get_error_independence_p_values()` が
   エラーなく動く
6. prior knowledge: `make_prior_knowledge()` で作った pk を渡して動作すること
   （test-prior_knowledge.R の使い方を流用）
7. print: `expect_output(print(res), "Multi-Group")`
8. bootstrap:
   - `n_sampling = 10L, reg_method = "ols", seed = 42` で
     `MultiGroupBootstrapResult`（長さ = 群数）が返る
   - 各要素で `get_probabilities()` が動く
   - 逐次の再現性: 同じ seed で2回実行して `expect_equal`
   - `compute_total_effects = FALSE` で total_effects が NULL、
     `get_total_causal_effects()` が明確なエラー
   - 並列（`parallel = TRUE, n_cores = 2L`）: 同 seed で再現、`skip_on_cran()`
     （既存 test-bootstrap.R の並列テストの流儀・PSOCK まわりの注意を必ず確認。
     メモリにある PSOCK の .libPaths 問題は既存コードのパターン踏襲で回避される）

実行: `devtools::test(filter = "multi_group")` → 全体 `devtools::test()`

## 7. ドキュメント

- roxygen2: lingam_direct.r / lingam_bootstrap.r の書式に倣う。
  `@references` は Shimizu (2012) Neurocomputing 81: 104-107
- Details に記載する事項:
  1. 全群で共通の因果順序を仮定し、探索目的関数をサンプルサイズ加重で合算すること
  2. `measure` 引数が無い理由（Python 版マルチグループはエントロピー近似のみ）
  3. 群ごとの解析（総合効果・独立性検定・プロット）は `get_group_result()` +
     既存関数で行うこと（例を付ける）
  4. bootstrap の総合効果はパス積和方式（4.2）
  5. bootstrap の再現性ポリシー（既存と同文）
- `@examples`: `generate_multi_group_sample()` → fit → print →
  `get_group_result(res, 1)` → `plot_adjacency` 等。bootstrap は
  `n_sampling = 10L, reg_method = "ols"` で直接実行、重いものは `\donttest{}`
- NEWS.md: 開発版セクションに追記
- vignette: 新セクション「Multi-Group Direct LiNGAM」
  （動機: 多施設・多群データで構造は共通だが係数が異なる場合 → 生成 → fit →
  群抽出して既存機能を適用 → bootstrap）
- _pkgdown.yml: reference に `lingam_multi_group`, `lingam_multi_group_bootstrap`,
  `get_group_result`, `generate_multi_group_sample` を追加

## 8. 完了条件・検証手順

R は PATH 外の `C:\R\R-4.6.0` にある。実行例:

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] `devtools::document()` が警告なく通る
- [ ] `devtools::test()` 全緑（既存テストを壊さない）
- [ ] `devtools::check()` で新規 ERROR/WARNING/NOTE なし
- [ ] 単群一致テスト（セクション6-3）が通る
- [ ] bootstrap の逐次・並列再現性テストが通る
- [ ] vignette がビルドできる

## 9. 遵守事項

- **既存の R ファイル（search_causal_order.r, fit_regression.r, lingam_bootstrap.r,
  paths.r 等）は一切変更しない**。内部関数は同一パッケージ内なのでそのまま呼べる
- NEWS.md / vignettes/lingamr.Rmd / _pkgdown.yml の編集は追記のみ
- テストコードの削除・コメントアウト禁止
- commit / push はしない（ユーザーが行う）
- エラーは握りつぶさない。意味のあるメッセージ付きで stop/warning
- 範囲外のリファクタリング禁止
