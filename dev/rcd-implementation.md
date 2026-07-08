# RCD (Repetitive Causal Discovery) 実装指示書

lingamr パッケージに RCD（潜在交絡下で「祖先集合 → 親 → 交絡ペア」を繰り返し推定する
因果探索）を追加する。Python cdt15/lingam の `RCD`（`rcd.py`、全561行）の R 移植。
**本指示書は原典561行とチュートリアルを精読した上で書かれている。**

アルゴリズムの3段構成（原典 fit、104-131 行目）:
1. `extract_ancestors`: 変数集合の組合せを走査し、各変数の**祖先集合 M** を反復的に拡大
2. `extract_parents`: 祖先集合から**親（直接原因）P** を絞り込む
3. `extract_vars_sharing_confounders`: 親子関係のないペアの残差相関から
   **同一潜在交絡を共有するペア C** を検出 → 隣接行列で NA

出力は Parce と同様「NA 入り隣接行列」だが、**因果順序は出力されない**
（代わりに祖先集合リストを持つ）。

## 0. 参照資料（実装前に必ず全文を読むこと）

- ソース: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/rcd.py
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/rcd.html
- 論文: T. N. Maeda and S. Shimizu. RCD: Repetitive causal discovery of linear
  non-Gaussian acyclic models with latent confounders. AISTATS 2020, PMLR 108: 735-745.

**着手前に必ず読む lingamr 側ファイル（再利用資産。すべて実装済み）:**

- `R/hsic.r` — `hsic_kernel_width(x)` / `hsic_gram_matrix(x, width)` /
  `hsic_test_gamma(X, Y)`（返り値 `list(stat =, p =)`）。**そのまま再利用**
- `R/f_correlation.r` — `f_correlation(x, y)`。**そのまま再利用**
- `R/lingam_parce.r` — `estimate_total_effect_parce()` /
  `get_error_independence_p_values_parce()` の構造（RCD 版はこれの変種）、
  NA 入り隣接行列の print の書き方
- `R/lingam_parce_bootstrap.r` — bootstrap の構造と **NA→0 置換ポリシー**（踏襲する）
- `R/paths.r` — `calculate_total_effect(adjacency_matrix, from, to)`（パス積和）
- `R/lingam_bootstrap.r` — `create_bootstrap_result()`（internal、再利用）
- `tests/testthat/test-lingam_parce.R` — NA 行列系のテストパターン

## 1. スコープ

| 成果物 | ファイル |
|---|---|
| `lingam_rcd()` + `RCDResult` + print | `R/lingam_rcd.r`（新規） |
| `estimate_total_effect_rcd()` / `get_error_independence_p_values_rcd()` | 同上 |
| `lingam_rcd_bootstrap()` | `R/lingam_rcd_bootstrap.r`（新規） |
| `generate_rcd_sample()`（チュートリアル準拠の潜在交絡データ） | `R/generate_rcd_sample.r`（新規） |
| testthat テスト | `tests/testthat/test-lingam_rcd.R`（新規） |
| man/NAMESPACE・NEWS.md・vignette・_pkgdown.yml・inst/WORDLIST | 既存様式で追記 |

### 実装しないもの（設計判断込み）

- `bw_method = "scott" / "silverman"`（statsmodels 依存。R/hsic.r は mdbs のみで実装済み。
  Parce 移植と同じ方針。引数自体を設けず Details に記載）
- prior knowledge（Python 版 RCD に無い）
- 新規パッケージ依存なし

### 含めるもの

- `MLHSICR = TRUE` オプション（HSIC 和最小化による頑健回帰。原典 186-224 行目。
  既存 hsic.r の部品 + `optim(method = "L-BFGS-B")` で実装可能）
- `independence = "fcorr"` オプション（f_correlation 実装済みのため低コスト）
- bootstrap（これまでの移植と同様に含める）

## 2. API 設計

```r
lingam_rcd <- function(X,
                       max_explanatory_num = 2L,
                       cor_alpha = 0.01,
                       ind_alpha = 0.01,
                       shapiro_alpha = 0.01,
                       MLHSICR = FALSE,
                       independence = "hsic",
                       ind_corr = 0.5)
```

既定値はすべて Python と同値。バリデーション（原典 70-90 行目 + lingamr 様式）:
`max_explanatory_num >= 1` の整数、各 alpha `>= 0`、`independence` は match.arg、
`ind_corr >= 0`、X の数値・NA・次元チェック（`nrow >= 3`。Shapiro-Wilk の下限）。

### 返り値: S3 クラス `RCDResult`

```r
result <- list(
  adjacency_matrix = B,      # p×p, lingamr 規約 B[i,j] = j -> i。交絡ペアは NA
  ancestors_list   = M       # 長さ p のリスト。各要素は祖先の整数ベクトル（1-based、昇順）
)
class(result) <- "RCDResult"
```

- **causal_order は無い**（Python 版にも無い。RCD は順序でなく祖先関係を推定する）
- dimnames に変数名。ancestors_list にも names を付ける
- `print.RCDResult()`: ヘッダ "RCD Result"。各変数の祖先集合を
  `M(x0) = {x1, x3, x5}` 形式で列挙し、続けて隣接行列
  （print.ParceLingamResult の NA 注記の書式に倣う）。`invisible(x)`

### 隣接行列の向き

原典は `B[xi, xj] <- coef`（xi を親 xj_list で回帰した係数。392-393 行目）なので
**lingamr 規約 `B[i,j] = j→i` とそのまま一致。転置不要**。

## 3. アルゴリズム（原典精読に基づく実装ガイド。1-based に読み替え）

### 3.1 共通部品

- **OLS 残差と係数** `get_resid_and_coef(X, endog, exog)`（原典 138-143 行目）:
  **切片あり**の OLS（sklearn LinearRegression 既定）。R では
  `fit <- lm.fit(cbind(1, X[, exog, drop = FALSE]), X[, endog])` として
  `resid = fit$residuals`, `coef = fit$coefficients[-1]`。
  【Parce 実装の教訓】`%*%` の結果は n×1 行列になるので、残差をベクトルとして
  使う前に `as.vector()` で剥がすこと（outer 等の非適合エラーの元）
- **相関判定**（原典 161-163 行目）: `cor.test(a, b, method = "pearson")$p.value < cor_alpha`
- **独立性判定**（原典 176-184 行目): hsic は `hsic_test_gamma(X, Y)$p > ind_alpha`、
  fcorr は `f_correlation(x, y) < ind_corr`
- **非ガウス性判定**（原典 154-159 行目）: U の**全変数**が
  `shapiro.test(Y[, xj])$p.value <= shapiro_alpha`（= 正規性が棄却される）なら TRUE。
  1つでも p > alpha なら FALSE で U をスキップ。
  **【R 固有の制約】`stats::shapiro.test` は n が 3〜5000 に制限される**
  （scipy は n > 5000 でも警告つきで計算する）。n > 5000 のときは**先頭 5000 行**で
  検定する（決定的なサブサンプル。乱数を使わない）方針で書いたが、実装
  （`R/lingam_rcd.r` の `rcd_is_non_gaussian_all()`）は既存の
  `test_residual_normality()`（`R/get_error_independence_p_values.r`）が
  採用しているランダムサブサンプリング（`sample(x, SHAPIRO_MAX_N)`）に
  意図的に合わせてある。パッケージ内で n > 5000 時の Shapiro-Wilk 挙動を
  統一する（＝この指示書の記述を優先する）ことを決めた場合は、両ファイルを
  同時に決定的サブサンプルへ移行すること。片方だけ変えると挙動が食い違う。
  この挙動（乱数消費・再現性は呼び出し側の RNG シード管理に依存する点）は
  `?lingam_rcd` の Details に既に記載済み。

### 3.2 祖先集合の抽出 `extract_ancestors`（原典 254-316 行目。アルゴリズムの心臓部）

```
M <- 各変数について空集合（長さ p のリスト）
l <- 1
hu_history <- 空（環境 or 名前付きリスト。キーは "i,j,k" 形式の文字列）
repeat {
  changed <- FALSE
  for (U in combn(1:p, l + 1) の各列) {                 # サイズ l+1 の変数集合
    H_U <- Reduce(intersect, M[U])                       # U の共通祖先
    key <- paste(U, collapse = ",")
    if (key が hu_history にあり、H_U が前回と同じ) next  # キャッシュ（原典 271-272 行目）

    Y <- 残差行列: H_U が空なら X そのもの。さもなくば U の各列を H_U で回帰した残差
         （U 以外の列は使われないので 0 のままでよい。原典 145-152 行目）

    if (!全 U が非ガウス) { hu_history 更新せず next }    # ←原典は続く2チェックの
    if (!全ペアが相関)    { 同上 next }                   #   失敗時 hu_history を更新しない
                                                          #   （307行目は sink 探索まで
                                                          #    到達した場合のみ実行される）
    sink_set <- {}
    for (xi in U) {
      xj_list <- setdiff(U, xi)
      if (exists_ancestor_in_U(M, U, xi, xj_list)) next   # 3.2.1
      if (is_independent_of_resid(Y, xi, xj_list)) sink_set に追加   # 3.2.2
    }
    if (length(sink_set) == 1) {
      xi <- sink_set[1]
      新規祖先 <- setdiff(U, xi) のうち M[[xi]] に無いもの
      あれば M[[xi]] に併合して changed <- TRUE
    }
    hu_history[key] <- H_U
  }
  if (changed) l <- 1
  else if (l < max_explanatory_num) l <- l + 1
  else break
}
```

**移植時の注意（読み違えやすい点）:**

1. `hu_history` の更新位置は原典 307 行目 = **U のループ本体の最後（sink 探索まで
   到達した場合）**。非ガウス性チェックや相関チェックで `continue` した場合は
   更新されない。この位置を変えると再走査の挙動が変わるので忠実に
2. `changed` が立つと **l を 1 にリセット**して小さい組合せから再走査する
3. `combn` の列挙順は Python `itertools.combinations` と同じ辞書順（両者とも昇順入力なら
   一致）。順序依存の結果になり得るので変えない
4. `H_U` の比較は集合として（`setequal`）

#### 3.2.1 `exists_ancestor_in_U`（原典 165-174 行目）

xi を sink 候補から外す条件: (a) いずれかの xj について xi ∈ M[[xj]]（xi が xj の祖先）、
または (b) xj_list ⊆ M[[xi]]（他の全変数が既に xi の祖先と判明済み）。

#### 3.2.2 `is_independent_of_resid`（原典 226-252 行目）

1. xi を xj_list で **OLS** 回帰した残差が、**全ての** xj と独立なら TRUE
2. 1つでも独立でない場合: `length(xj_list) == 1` または `MLHSICR = FALSE` なら FALSE
3. `MLHSICR = TRUE` かつ説明変数2個以上なら **MLHSICR 回帰**（3.3）で残差を作り直し、
   全 xj と独立なら TRUE

### 3.3 MLHSICR 回帰（原典 186-224 行目）

HSIC の和を最小化する回帰係数を数値最適化で求める:

- 各 xj のカーネル幅 `hsic_kernel_width()` と中心化グラム行列 Lc を**事前計算**
- 初期値は OLS 係数
- 目的関数（係数ベクトル coef を受け取る）:
  ```
  resid <- Y[, xi] - Σ_j coef[j] * Y[, xj]
  width <- width_xi - Σ_j coef[j] * width_list[j]     # 幅も線形結合（原典 207-208 行目。
                                                       #  奇妙に見えるが原典仕様。変えない）
  Kc <- hsic_gram_matrix(resid, width) の中心化行列
  objective <- Σ_j sum(t(Kc) * Lc_j) / n              # HSIC 統計量の和
  ```
  HSIC 統計量 `sum(t(Kc) * Lc) / n` は Python の `hsic_teststat` 相当。R/hsic.r には
  独立関数として存在しないので **lingam_rcd.r 内にインラインで書く**
  （**R/hsic.r は変更しない**）。`hsic_gram_matrix()` の返り値の形
  （K と Kc のリストのはず）は実装を読んで確認すること
- 最適化: `optim(par = initial_coef, fn = sum_empirical_hsic, method = "L-BFGS-B")`
  （勾配は渡さない = 数値微分。Python の `approx_grad=True` と同じ）
- 返り値: 最適化後の係数で作った残差と係数

### 3.4 親の抽出 `extract_parents`（原典 318-343 行目）

各 xi と各 xj ∈ M[[xi]] について:

- `zi` = xi を `M[[xi]] \ {xj}` で回帰した残差（空なら X[, xi] そのまま）
- `wj` = xj を `M[[xi]] ∩ M[[xj]]` で回帰した残差（空なら X[, xj] そのまま）
- `zi` と `wj` が相関（pearson, cor_alpha）していれば xj は xi の親 → P[[xi]] に追加

### 3.5 交絡ペアの検出（原典 352-366 行目）

全ペア (i, j)（どちらの方向にも親子関係が無いもの）について:

- 各変数をそれぞれの親 P で回帰した残差同士が相関していれば、
  同一潜在交絡を共有するペアとして C[[i]], C[[j]] に相互登録

### 3.6 隣接行列（原典 368-408 行目）

- 各 xi: 親 `sort(P[[xi]])` で OLS 回帰し `B[xi, 親] <- coef`（親なしはゼロ行のまま）
- 各 xi: `B[xi, C[[xi]]] <- NA`（C は対称に登録済みなので両方向とも NA になる。
  C ペアは親子関係が無いことが保証されているので係数の上書きは起きない）

### 3.7 総合効果・誤差独立性（原典 410-496 行目）

`R/lingam_parce.r` の `estimate_total_effect_parce()` /
`get_error_independence_p_values_parce()` を**読んでから**、その構造を流用した
RCD 版を lingam_rcd.r に実装する。Parce 版との差分は:

- `estimate_total_effect_rcd(X, rcd_result, from_index, to_index)`:
  順序の妥当性チェックが「causal_order の前後」でなく
  **「to_index ∈ ancestors_list[[from_index]] なら warning」**（原典 415-420 行目）。
  from が交絡行（NA あり）なら warning + NA（Parce と同じ）
- `get_error_independence_p_values_rcd(X, rcd_result)`: Parce 版と同一ロジック
  （E = X - B X'、NA 絡み変数のペアは NA、他は HSIC p 値）

## 4. bootstrap `lingam_rcd_bootstrap()`（原典 523-561 行目）

`R/lingam_parce_bootstrap.r` の構造を踏襲（引数バリデーション → run_one + tryCatch 失敗
スキップ → PSOCK 並列 + L'Ecuyer → 集計。verbose・seed・parallel・n_cores 引数も同様）。
RCD 固有の差分:

1. 総合効果のループが**祖先リスト駆動**（原典 557-559 行目）:
   `for (to in 1:p) for (from in ancestors_list[[to]])` →
   `calculate_total_effect(B, from, to)`（R/paths.r。パス積和）。
   from が交絡行なら NA（estimate_total_effect_rcd と同じ判定）
2. **NA→0 置換ポリシー**は lingam_parce_bootstrap.r と同一
   （BootstrapResult 格納時に NA を 0 に。numpy の nan 比較セマンティクスとの整合。
   実装済みコードの該当箇所を読んで同じ書き方をする）
3. `create_bootstrap_result()` を再利用。**causal_orders は NULL**
   （RCD に因果順序が無い）→ `get_causal_order_stability()` 不可を Details に明記
4. RCD の fit は HSIC 多用で 1 回が重い。roxygen2 の例は `n_sampling` を小さく

## 5. データ生成関数 `generate_rcd_sample()`

チュートリアル準拠（観測6変数 + 潜在1変数、n = 300、
ノイズは `rnorm(n, 0, 0.5)^3` の超ガウス分布）:

```
x5 <- e();  x6 <- e()                 # x6 は潜在
x1 <- 0.6*x5 + e();  x3 <- 0.5*x5 + e()
x0 <- 1.0*x1 + 1.0*x3 + e()
x2 <- 0.8*x0 - 0.6*x6 + e()
x4 <- 1.0*x0 - 0.5*x6 + e()           # e() = rnorm(n, 0, 0.5)^3
```

- シグネチャ例: `generate_rcd_sample(n = 300, seed = NULL)`
- 返り値: `list(data = <x0..x5 の data.frame>, adjacency_matrix = <真の B。x2-x4 間は NA>,
  ancestors_list = <真の祖先集合リスト>, confounded_pair = c(3, 5))`
  （既存 generate_* の形式に合わせる。1-based で x2 → 列3、x4 → 列5）
- 真の祖先集合（チュートリアルの期待値、1-based 列番号で）:
  M(x0)={x1,x3,x5}, M(x1)={x5}, M(x2)={x0,x1,x3,x5}, M(x3)={x5},
  M(x4)={x0,x1,x3,x5}, M(x5)=∅ — **テストオラクルとして使える**

## 6. テスト（test-lingam_rcd.R）

HSIC 多用で重いので n = 300（チュートリアルと同じ）を基本とし、
実行時間が問題になるテストのみ `skip_on_cran()`。
**【Parce 実装の教訓】seed は事前に複数スイープして安定するものに固定すること**
（検定ベースのアルゴリズムは seed により誤検出が起きる。教訓の実例:
generate_lingam_sample_6 の seed=42 で誤ブロック化 → seed=1 に変更で回避した）。

1. 構造: `expect_s3_class(res, "RCDResult")`, フィールド・次元・dimnames、
   ancestors_list が長さ p のリスト
2. バリデーション: 非数値・NA 入り X、`max_explanatory_num = 0`、負の alpha、
   `independence = "foo"` で `expect_error()`
3. **祖先集合の回復（最重要）**: `generate_rcd_sample(seed 固定)` で
   ancestors_list が真の集合と一致（または主要な関係を含む。検定力次第で
   「真の祖先 ⊆ 推定祖先」の包含チェックに緩めてよい。seed スイープで安定させる）
4. **交絡ペアの検出**: `B[x2 行, x4 列]` と `B[x4 行, x2 列]` が NA、
   非交絡の真の辺（x1→x0, x0→x2 等）の係数が符号込みで回復
5. `MLHSICR = TRUE` で valid な結果（値の厳密比較はしない。エラーなく動き
   構造が返ること）
6. `independence = "fcorr"` で valid な結果
7. `estimate_total_effect_rcd`: 交絡を受ける変数を from → warning + NA、
   健全ペア → 数値
8. `get_error_independence_p_values_rcd`: NA 絡みペアが NA、他が [0,1]
9. print: `expect_output(print(res), "RCD")` と祖先集合の表示
10. bootstrap: `n_sampling = 3-5, seed = 42` で BootstrapResult、
    `get_probabilities()` が NA を返さず動く、逐次再現性。
    並列は `skip_on_cran()` 付き
11. n > 5000 の Shapiro 制限: 6000 行程度のデータでエラーにならないこと
    （`skip_on_cran()` 可。HSIC が重いので変数2-3個の小さい構造でよい）

実行: `devtools::test(filter = "rcd")` → 全体 `devtools::test()`

## 7. ドキュメント

- roxygen2: 既存様式。`@references` は Maeda & Shimizu (AISTATS 2020, PMLR 108:735-745)
- Details に記載する事項:
  1. 3段構成（祖先 → 親 → 交絡ペア）と、causal_order でなく ancestors_list を返すこと
  2. NA = 同一潜在交絡を共有するペア
  3. `max_explanatory_num` と計算量（組合せ × HSIC O(n²)。変数が多いと急激に重い）
  4. MLHSICR の意味（OLS 残差が独立にならないときの HSIC 最小化回帰）と計算コスト
  5. Shapiro-Wilk の n ≤ 5000 制限と先頭5000行での検定
  6. bw_method を持たないこと（mdbs 固定）、bootstrap で
     `get_causal_order_stability()` が使えないこと
- `@examples`: `generate_rcd_sample()` → fit → print → ancestors_list の読み方 →
  総合効果の NA warning。重い部分は `\donttest{}`
- NEWS.md・vignette（「Latent Confounders: RCD」小節。Parce のセクションの隣に置き、
  Parce との使い分け — Parce は順序不明ブロック、RCD はペア単位の交絡検出 — を1-2文）・
  _pkgdown.yml 追記
- **`inst/WORDLIST`**: 新出の専門用語（RCD, MLHSICR, Maeda, AISTATS 等）を
  `spelling::update_wordlist()` で追加（Parce 実装時に NOTE の原因になった）

## 8. 完了条件・検証手順

R は PATH 外の `C:\R\R-4.6.0` にある（check には RSTUDIO_PANDOC 設定が必要 —
過去の実装セッションと同様）。実行例:

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] `devtools::document()` が警告なく通る
- [ ] `devtools::test()` 全緑（既存テストを壊さない。現状 605 PASS が基準）
- [ ] `devtools::check()` 0 errors / 0 warnings / 0 notes
- [ ] 祖先集合・交絡ペア検出テスト（セクション6-3, 6-4）が通る
- [ ] bootstrap 再現性テストが通る
- [ ] vignette がビルドできる

## 9. 遵守事項

- **既存の R ファイル（hsic.r, f_correlation.r, lingam_parce*.r, paths.r,
  lingam_bootstrap.r 等）は一切変更しない**。同一パッケージ内なのでそのまま呼べる。
  hsic_teststat 相当が必要なら lingam_rcd.r 内にインラインで書く
- NEWS.md / vignette / _pkgdown.yml / inst/WORDLIST は追記のみ
- テストコードの削除・コメントアウト禁止
- commit / push はしない（ユーザーが行う）
- エラーは握りつぶさない。意味のあるメッセージ付きで stop/warning
- 範囲外のリファクタリング禁止
