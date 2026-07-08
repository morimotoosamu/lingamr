# evaluate_model_fit 実装指示書

lingamr パッケージに `evaluate_model_fit()`（推定した隣接行列の SEM 適合度指標による評価）を
追加する。Python cdt15/lingam の `lingam.utils.evaluate_model_fit`
（utils/__init__.py 942-1018 行目、約77行）の R 移植。
**本指示書は原典を全文精読した上で書かれている。**

機能の要点: 隣接行列（因果グラフ）を SEM（構造方程式モデル）の仕様に変換して当てはめ、
CFI・RMSEA・AIC/BIC 等の適合度指標を返す。Python は semopy に丸投げしており、
**R では lavaan に丸投げする**のが対応方針。アルゴリズム的な移植は「隣接行列 →
モデル記述文字列」の変換部だけ。

## 0. 参照資料

- 原典: https://raw.githubusercontent.com/cdt15/lingam/master/lingam/utils/__init__.py
  の `evaluate_model_fit`（942-1018 行目）
- チュートリアル: https://lingam.readthedocs.io/en/latest/tutorial/evaluate_model_fit.html
- lavaan のモデル構文とfitMeasures: `?lavaan::sem`, `?lavaan::fitMeasures`

着手前に読む lingamr 側ファイル:

- `R/fit_regression.r` の `check_glmnet_available()`（34行目）— オプション依存の
  条件付きロードパターン。lavaan にも同じパターンを使う
- `R/lingam_direct.r` — バリデーション・roxygen2 の様式
- `DESCRIPTION` — Suggests の現状

## 1. スコープ

| 成果物 | ファイル |
|---|---|
| `evaluate_model_fit()` | `R/evaluate_model_fit.r`（新規） |
| testthat テスト | `tests/testthat/test-evaluate_model_fit.R`（新規） |
| **DESCRIPTION の Suggests に `lavaan` を追記** | `DESCRIPTION` |
| man/NAMESPACE・NEWS.md・vignette・_pkgdown.yml | 既存様式で追記 |

**新規依存**: lavaan（Suggests。Imports にしない）。
`requireNamespace("lavaan", quietly = TRUE)` で確認し、無ければ
「Package 'lavaan' is required for evaluate_model_fit(). Install it with
install.packages("lavaan").」の旨で `stop(call. = FALSE)`。

## 2. API 設計

```r
evaluate_model_fit <- function(adjacency_matrix, X, is_ordinal = NULL)
```

- `adjacency_matrix`: p×p 数値行列（**lingamr 規約 B[i,j] = j→i**。NA 可）。
  **R 独自の親切として `LingamResult` / `ParceLingamResult` / `LiMResult` 等の
  オブジェクトも受け付け、`$adjacency_matrix` を自動抽出する**
  （`inherits()` でなく「リストで adjacency_matrix 要素を持つか」で判定してよい）
- `X`: n×p の数値行列または data.frame（NA 不可）
- `is_ordinal`: 長さ p の logical または 0/1 ベクトル（NULL なら全て連続扱い）。
  TRUE の変数は順序カテゴリカルとして扱う
- 返り値: 適合度指標の 1 行 data.frame（セクション 3.3）

バリデーション（原典 962-976 行目 + lingamr 様式）:
正方行列であること、`ncol(X) == ncol(adj)`、X に NA が無いこと、
`length(is_ordinal) == ncol(adj)`。エラーは英語 + `call. = FALSE`。

## 3. 実装

### 3.1 隣接行列 → lavaan モデル構文（原典 978-1008 行目の変換ロジック）

変数名は X の列名を使う（原典は x0, x1, ... 固定だが、R では実列名を活かす。
列名が無ければ `x0, x1, ...` を生成 — 既存 `get_var_names()` を再利用）。
**lavaan 構文で使えない文字を含む列名への対策として、内部では `x0...` の機械名で
構文を組み、結果には影響させない実装でもよい**（どちらかに決めて一貫させる。
機械名方式が安全で推奨）。

各行 i（= 変数 i への流入）について:

1. 行が「NA なし かつ 全要素ゼロ」なら外生変数 → 式を作らない（原典 984-985 行目）
2. それ以外は回帰式 `xi ~ ...` を作る:
   - `B[i, j]` が非ゼロ（`!isTRUE(all.equal(elem, 0))` 相当。原典は `np.isclose`）
     → 右辺に `xj` を追加
   - `B[i, j]` が **NA**（ParceLiNGAM 等の「潜在交絡ペア」）
     → 右辺に潜在変数 `eta_i_j`（i < j になるよう正規化した名前。原典 990-994 行目）
       を追加。同じ eta は一度だけ登録

3. 潜在交絡の表現（**semopy との対応で最重要の設計点**）:
   原典は `DEFINE(latent) eta_i_j` を宣言し、`xi ~ eta_i_j` と `xj ~ eta_i_j` の
   回帰に入れる（= 2変数の潜在共通原因）。lavaan には DEFINE(latent) が無いので
   等価な表現に置き換える。**推奨: 残差共分散 `xi ~~ xj`**。
   2 指標の潜在共通原因（負荷1つ固定）と残差共分散はパラメータ数・適合度とも
   等価であり、lavaan では残差共分散が標準的なイディオム。
   roxygen2 の Details に「semopy 実装は潜在変数、本実装は等価な残差共分散で表現」と明記。
   注意: NA ペアの2変数がどちらも回帰式を持たない（他に親が無い）場合、
   外生変数間の共分散になるが lavaan ではそのまま `xi ~~ xj` で書ける
4. `is_ordinal` の指定変数は `lavaan::sem(..., ordered = <変数名ベクトル>)` に渡す
   （原典の DEFINE(ordinal) 相当。lavaan は自動的に WLSMV 系推定に切り替わる）

### 3.2 当てはめ

```r
fit <- lavaan::sem(model_str, data = as.data.frame(X_named), ordered = ordered_vars)
```

- 収束失敗・警告は握りつぶさない（lavaan の warning はそのまま伝播させる）
- モデル文字列が空（全変数が外生 = 辺なしグラフ）の場合の挙動を確認し、
  lavaan がエラーになるなら「the adjacency matrix has no edges」の明確なエラーを先に出す

### 3.3 返り値

`lavaan::fitMeasures(fit)` から主要指標を取り出し、**1 行の data.frame** で返す。
semopy の calc_stats が返す列（DoF, DoF Baseline, chi2, chi2 p-value, chi2 Baseline,
CFI, GFI, AGFI, NFI, TLI, RMSEA, AIC, BIC, LogLik）に対応する lavaan 指標:

| semopy 列 | lavaan fitMeasures 名 |
|---|---|
| DoF | df | 
| chi2 / chi2 p-value | chisq / pvalue |
| chi2 Baseline | baseline.chisq |
| CFI / GFI / AGFI / NFI / TLI | cfi / gfi / agfi / nfi / tli |
| RMSEA | rmsea |
| AIC / BIC / LogLik | aic / bic / logl |

- ordered 指定時は一部指標が定義されない（aic/bic/logl 等）— `NA` で埋める
  （`fitMeasures` に存在しない名前は NA にするユーティリティを内部に書く）
- **数値は semopy と完全一致しない**（推定エンジン・既定オプションの差）。
  Details に明記。指標の定義としては同じものを返す

## 4. テスト（test-evaluate_model_fit.R）

冒頭で `skip_if_not_installed("lavaan")`（glmnet の既存テストパターンに倣う。
lavaan 不在時のエラーメッセージのテストだけは `local_mocked_bindings()` で
チェック関数を差し替えて行う — test-lingam_direct.R:195-212 参照）。

1. **正しいモデルで良い適合**: `generate_lingam_sample_6()` → `lingam_direct()` →
   `evaluate_model_fit(result, X)` が 1 行 data.frame を返し、CFI が高い（> 0.95 目安）、
   RMSEA が小さいことを確認
2. **誤ったモデルで悪化**: 真の構造から辺の向きを逆転させた行列で CFI が低下する
   （厳密な閾値でなく「正しいモデルより悪い」比較にする）
3. **行列でなく結果オブジェクトを渡す**: LingamResult を直接渡しても同じ結果
4. **NA（潜在交絡）入り行列**: 2変数間を NA にした行列でエラーなく動き、
   指標が返る（ParceLingamResult の出力を意識）
5. バリデーション: 非正方行列、次元不一致、X に NA、is_ordinal 長さ不正で
   `expect_error()`
6. `is_ordinal`: 1変数を 0/1 に離散化して ordinal 指定し、エラーなく動く
   （値の厳密比較はしない）
7. 辺なし行列の明確なエラー（3.2 参照）

## 5. ドキュメント

- roxygen2: `@examples` は `skip` できないので
  `if (requireNamespace("lavaan", quietly = TRUE)) { ... }` で包むか `\donttest{}`
- Details: (1) lavaan 依存（Suggests）、(2) NA ペア = 潜在交絡の残差共分散表現、
  (3) semopy と数値は一致しない、(4) B[i,j]=j→i 規約
- `@references` は semopy でなく lavaan と、適合度指標の一般的な解釈に触れる程度でよい
- DESCRIPTION: Suggests に `lavaan` をアルファベット順で追記（**この1行の追記のみ許可**）
- NEWS.md・vignette（「Evaluating Model Fit」小節。lingam_direct の結果を評価する例）・
  _pkgdown.yml 追記

## 6. 完了条件

```powershell
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::document(); devtools::test()"
& "C:\R\R-4.6.0\bin\Rscript.exe" -e "devtools::check()"
```

- [ ] document / test / check 全通過（既存テスト無破壊、新規 NOTE なし）
- [ ] lavaan 未インストール環境相当（モック）で意味のあるエラー
- [ ] NA 入り行列（潜在交絡）のテストが通る

## 7. 遵守事項

- 既存 R ファイルは変更しない。DESCRIPTION は Suggests への1行追記のみ
- NEWS.md / vignette / _pkgdown.yml は追記のみ
- テストコードの削除・コメントアウト禁止、commit/push はユーザー
- エラーは握りつぶさない。範囲外リファクタリング禁止
