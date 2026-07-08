# bootstrap_with_imputation 実装指示書

lingamr パッケージに `bootstrap_with_imputation()`（欠測データに対する多重代入つき bootstrap
因果探索）を追加する。Python cdt15/lingam の `lingam.tools.bootstrap_with_imputation`
（tools/__init__.py、全315行）の R 移植。**本指示書は原典315行を全文精読した上で書かれている。**

## 【前提条件 — 最初に確認すること】

このアルゴリズムは各 bootstrap 標本を多重代入して得た n_repeats 個のデータセットに
**共通の因果構造を仮定して MultiGroupDirectLiNGAM で同時推定する**（原典 110-137 行目）。
したがって **`lingam_multi_group()`（指示書: dev/multi-group-direct-lingam-implementation.md）
が実装済みであることが前提**。未実装ならこの指示書の作業は開始できない。
`R/lingam_multi_group.r` の存在と `devtools::test(filter = "multi_group")` の成功を
最初に確認すること。

## 0. 参照資料

- 原典: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/tools/__init__.py
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/bootstrap_with_imputation.html
  （欠測は x5 に MCAR 10%、n_sampling=30, n_repeats=5 の例）

着手前に読む lingamr 側ファイル:

- `R/lingam_multi_group.r`（実装済み前提）— fit の呼び出し方と返り値
- `R/lingam_bootstrap.r` — バリデーション様式（並列化は本機能では**行わない**。後述）
- `R/fit_regression.r` の `check_glmnet_available()` — オプション依存パターン（mice に流用）
- `DESCRIPTION` — Suggests の現状

## 1. スコープ

| 成果物 | ファイル |
|---|---|
| `bootstrap_with_imputation()` + 結果クラス + print | `R/bootstrap_with_imputation.r`（新規） |
| `as_bootstrap_result()`（repeats 集約ヘルパ、R 独自） | 同上 |
| testthat テスト | `tests/testthat/test-bootstrap_with_imputation.R`（新規） |
| **DESCRIPTION の Suggests に `mice` を追記** | `DESCRIPTION` |
| man/NAMESPACE・NEWS.md・vignette・_pkgdown.yml | 既存様式で追記 |

**新規依存**: mice（Suggests）。既定の多重代入に使う。
Python の `IterativeImputer(sample_posterior=True)`（ベイズ線形回帰の連鎖方程式 + 事後分布
サンプリング）に最も近い R の標準は `mice(method = "norm")`（ベイズ線形回帰）。
`requireNamespace` パターンで条件付きロードし、独自 `imputer` を渡された場合は
mice 不要で動くこと。

### 実装しないもの（設計判断込み）

- 並列実行（Python 原典も逐次のみ。将来必要になったら lingam_bootstrap.r の
  PSOCK パターンで拡張できる旨を Details に一言）
- Python の抽象クラス `BaseMultipleImputation` / `BaseMultiGroupCDModel`
  （R では**関数を渡す**イディオムに置き換える。下記 `imputer` / `cd_fit` 引数）
- `before_imputation` フック（原典 128, 310-311 行目。既定実装は何もしない pass。
  R の `cd_fit` は代入済みリストだけ受け取る設計とし、必要になったら拡張）

## 2. API 設計

```r
bootstrap_with_imputation <- function(X,
                                      n_sampling,
                                      n_repeats = 10L,
                                      imputer = NULL,
                                      cd_fit = NULL,
                                      prior_knowledge = NULL,
                                      apply_prior_knowledge_softly = FALSE,
                                      seed = NULL,
                                      verbose = TRUE)
```

- `X`: n×p の数値行列 / data.frame。**NA を含んでよい（それがこの関数の目的）**。
  逆に NA が1つも無い場合は「欠測が無いので `lingam_direct_bootstrap()` を使うべき」
  という警告を出して続行する（R 独自の親切）
- `n_sampling`: bootstrap 反復数（正の整数）
- `n_repeats`: 各 bootstrap 標本への代入回数（正の整数。`imputer` 指定時は無視 —
  原典 49-51 行目と同じ）
- `imputer`: `function(X_boot)` → **代入済み行列のリスト**を返す関数。NULL なら既定の
  mice 実装（3.2）。返り値はチェックする（3.4）
- `cd_fit`: `function(X_list)` → `list(causal_order = <整数ベクトル>,
  adjacency_matrices = <リスト or 3次元配列>)` を返す関数。NULL なら
  `lingam_multi_group()` を使う既定実装。返り値はチェックする（3.4）
- `prior_knowledge` / `apply_prior_knowledge_softly`: `cd_fit = NULL` のときだけ
  既定の lingam_multi_group に渡す（原典 58-69 行目と同じ）。`cd_fit` 指定時に
  これらも指定されたら warning
- `seed`: 再現用。冒頭で一度 `set.seed(seed)`（復元抽出と mice の乱数が両方これに従う。
  mice 呼び出しに seed 引数は渡さず、グローバル RNG に従わせる —
  `mice(..., seed = NA)` 既定で printFlag = FALSE を指定）

## 3. 実装（原典 86-157 行目）

### 3.1 メインループ

```
for (i in 1:n_sampling) {
  resampled_index <- sample(n, replace = TRUE)      # NA を含んだまま復元抽出（原典 124 行目）
  X_boot <- X[resampled_index, , drop = FALSE]
  datasets <- imputer(X_boot)                        # n_repeats 個の完全データ
  datasets のチェック（3.4）
  cd_res <- cd_fit(datasets)                         # 共通構造で同時推定
  cd_res のチェック（3.4）
  # 代入値の記録（原典 140-145 行目）: X_boot が NA だった位置の値だけ残し、他は NA
  imputation_result[r, pos] <- datasets[[r]][pos]    # pos = is.na(X_boot)
  結果を蓄積
}
```

- 1反復を `tryCatch` で包み、失敗反復は warning + スキップ、全滅なら stop
  （lingam_bootstrap.r:128-166 の設計を踏襲。原典には無い R 側の頑健化で、
  mice が特定の resample で収束しないケースに備える）
- verbose の進捗表示は lingam_bootstrap.r:231-236 の様式

### 3.2 既定 imputer（原典 279-292 行目の対応）

```r
default_imputer <- function(X_boot, n_repeats) {
  # mice は m 個の代入を一括生成できる（Python は n_repeats 回 fit_transform を回すが同義）
  imp <- mice::mice(as.data.frame(X_boot), m = n_repeats, method = "norm",
                    printFlag = FALSE)
  lapply(seq_len(n_repeats), function(k) as.matrix(mice::complete(imp, k)))
}
```

- `method = "norm"`（ベイズ線形回帰）が IterativeImputer(sample_posterior=TRUE) の対応。
  Details に「Python は scikit-learn IterativeImputer、本実装は mice::norm のため
  数値は一致しない（多重代入という設計は同じ）」と明記
- mice が欠測の無い列だけのデータや全行欠測などで失敗する場合のエラーは
  そのまま伝播させる（握りつぶさない）

### 3.3 既定 cd_fit

```r
default_cd_fit <- function(X_list, prior_knowledge, apply_prior_knowledge_softly) {
  res <- lingam_multi_group(X_list,
                            prior_knowledge = prior_knowledge,
                            apply_prior_knowledge_softly = apply_prior_knowledge_softly)
  list(causal_order = res$causal_order,
       adjacency_matrices = res$adjacency_matrices)
}
```

### 3.4 フック返り値のチェック（原典 187-198, 251-276 行目の対応）

Python は `_check_imputer_outout` / `_check_cd_output` で厳密検証している。R でも:

- imputer 出力: リストであること、各要素が `dim = c(nrow(X_boot), p)` の数値行列で
  **NA を含まない**こと。違反時 "The return value of imputer violates its
  specification: ..." で stop
- cd_fit 出力: 2要素（causal_order, adjacency_matrices）。causal_order は
  長さ p で `sort() == 1:p`（1-based の完全順列）であること（原典 265-266 行目の
  0-based 検査の 1-based 版）。adjacency_matrices は要素数 = 代入数、各 p×p。
  違反時も同様の明確なメッセージ

### 3.5 返り値: S3 クラス `ImputationBootstrapResult`

```r
result <- list(
  causal_orders      = <n_sampling × p の整数行列>,             # 1-based
  adjacency_matrices = <array(n_sampling, n_repeats, p, p)>,    # lingamr 規約 B[i,j]=j→i
  resampled_indices  = <n_sampling × n の整数行列>,
  imputation_results = <array(n_sampling, n_repeats, n, p)>     # 代入位置以外は NA
)
class(result) <- "ImputationBootstrapResult"
```

- Python の返り値タプル（原典 152-157 行目）と同じ4点セット。dimnames に変数名
- `print.ImputationBootstrapResult()`: n_sampling / n_repeats / 変数数 / 欠測セル数を表示
- **既存 `BootstrapResult` とは形が違う**（n_repeats 次元がある）ため
  そのままでは get_probabilities 等が使えない。そこで R 独自ヘルパを追加:

```r
as_bootstrap_result <- function(x, aggregate = c("median", "mean"))
```

  n_repeats 次元を `aggregate` で潰して（チュートリアルも median 集約を示している）
  既存の `create_bootstrap_result()`（lingam_bootstrap.r:299-308、internal）で
  `BootstrapResult` に変換する。total_effects は NULL（未計算）、causal_orders は
  そのまま渡す。これで `get_probabilities()` / `get_causal_direction_counts()` /
  `get_causal_order_stability()` が使える（`get_total_causal_effects()` は
  total_effects = NULL の明確なエラーになる — 既存仕様どおり）。
  roxygen2 でこの使い方を必ず実演する

## 4. テスト（test-bootstrap_with_imputation.R）

冒頭 `skip_if_not_installed("mice")`（mice 不在エラーメッセージのテストのみモックで）。
mice + 反復 fit は重いので n_sampling / n_repeats は小さく（例: 5 × 3）、
`reg_method` 既定が重ければデータを小さく（n = 300、変数4-6個）。

1. データ準備: `generate_lingam_sample_6()` のデータに MCAR で 10% の NA を
   1変数に注入（チュートリアル方式。seed 固定）
2. 構造: `ImputationBootstrapResult` のクラス・4フィールドの次元が仕様どおり
3. **imputation_results のマスク**: 非 NA 要素の位置が「元 bootstrap 標本の NA 位置」
   と一致する（原典 140-145 行目のセマンティクス）
4. 再現性: 同じ `seed` で2回実行して `expect_equal`
5. 欠測なしデータ → warning（2章の R 独自仕様）
6. カスタム `imputer`（例: 列平均代入を n_repeats 回返すだけの関数）で動く。
   不正な imputer（NA を残す・次元が違う）で仕様違反エラー
7. カスタム `cd_fit`（lingam_multi_group をそのまま包んだもの）で動く。
   不正な causal_order（順列でない）で仕様違反エラー
8. `as_bootstrap_result()`: BootstrapResult が返り、`get_probabilities()` が動く。
   `get_total_causal_effects()` は明確なエラー
9. バリデーション: n_sampling = 0、n_repeats = 0、prior_knowledge 次元不正で
   `expect_error()`
10. print: `expect_output(print(res), ...)`

実行: `devtools::test(filter = "bootstrap_with_imputation")` → 全体

## 5. ドキュメント

- roxygen2 Details に記載する事項:
  1. 手順（bootstrap → 多重代入 → 共通構造の同時推定）と、共通構造仮定の意味
  2. 前提として lingam_multi_group を内部で使うこと
  3. 既定 imputer は mice::norm。Python（IterativeImputer）と数値は一致しない
  4. `imputer` / `cd_fit` による差し替え方法（仕様と例）
  5. `as_bootstrap_result()` で既存の bootstrap 集計関数群につなげること
  6. 逐次実行のみ（現状）
- `@examples`: 小さな n で NA 注入 → 実行 → `as_bootstrap_result` → get_probabilities。
  mice 依存なので `\donttest{}` か requireNamespace ガード
- DESCRIPTION: Suggests に `mice` をアルファベット順で追記（**この1行の追記のみ許可**）
- NEWS.md・vignette（「Causal Discovery with Missing Data」小節）・_pkgdown.yml 追記

## 6. 完了条件

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] 前提確認: lingam_multi_group が実装済み・テスト緑
- [ ] document / test / check 全通過（既存テスト無破壊、新規 NOTE なし）
- [ ] imputation_results マスクのテスト（セクション4-3）が通る
- [ ] seed 再現性テストが通る

## 7. 遵守事項

- 既存 R ファイルは変更しない（`create_bootstrap_result` 等は同一パッケージ内なので
  そのまま呼べる）。DESCRIPTION は Suggests への1行追記のみ
- NEWS.md / vignette / _pkgdown.yml は追記のみ
- テストコードの削除・コメントアウト禁止、commit/push はユーザー
- エラーは握りつぶさない。範囲外リファクタリング禁止
