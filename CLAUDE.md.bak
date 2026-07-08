# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

```r
devtools::load_all()        # 開発時にパッケージを読み込む
devtools::document()        # roxygen2からRdファイルとNAMESPACEを生成
devtools::check()           # R CMD CHECK（コミット前に必ず実行）
devtools::test()            # テストを全実行
devtools::build_readme()    # README.Rmdからインデックスを生成

# 個別テストの実行
testthat::test_file("tests/testthat/test-lingam_direct.R")

# 特定関数に関連するテストだけ実行
devtools::test(filter = "lingam_direct")

# CRAN提出前確認
devtools::check_win_devel()
devtools::release()
```

## Architecture

### S3クラスと対応ファイル

| クラス | ソースファイル | 説明 |
|--------|--------------|------|
| `LingamResult` | `R/lingam_direct.r` | Direct LiNGAMの推定結果 |
| `BootstrapResult` | `R/lingam_bootstrap.r` | Bootstrap安定性評価 |
| `VARLiNGAMResult` | `R/lingam_var.r` | VAR-LiNGAMの推定結果 |
| `VARBootstrapResult` | `R/lingam_var_bootstrap.r` | VARのBootstrap結果 |
| `lingam_summary` | `R/summary_lingam.r` | 包括的な適合度サマリー |

**重要な規約:** 隣接行列 `B[i,j]` は変数 j → i の因果係数を表す。

### 主要アルゴリズム

**Direct LiNGAM** (`R/lingam_direct.r`):
- `measure`: "pwling"（ペアワイズ独立）or "kernel"（カーネルベース）
- `reg_method`: "ols" / "lasso" / "adaptive_lasso"（デフォルト） / "ridge"
- 回帰バックエンドは `R/fit_regression.r`（glmnetが必要な手法あり）
- 因果順序探索は `R/search_causal_order.r`

**VAR-LiNGAM** (`R/lingam_var.r`):
- VARモデルを当てはめ、残差にDirect LiNGAMを適用
- ラグ選択基準: "bic"（デフォルト）/ "aic" / "hqic" / "fpe"
- `adjacency_matrices[1,,]` が瞬時構造B0、以降がラグ行列

**Bootstrap並列実行** (`R/lingam_bootstrap.r`):
- `parallel::makePSOCKcluster()` + L'Ecuyer RNGストリームで再現性を保証
- `n_cores` が変わると数値結果が変わる（設計上の仕様）

### S3メソッドの在処

| 関数 | ファイル |
|------|---------|
| `tidy()`, `glance()` | `R/tidiers.r` |
| `autoplot.LingamResult()` | `R/autoplot.r` |
| `print.*()` | 各クラスのソースファイル |

### オプション依存パッケージ

- **glmnet**: lasso/ridge系の回帰が必須
- **DiagrammeR**: `plot_adjacency()` のインタラクティブDAG
- **ggplot2**: `autoplot()` と診断プロット
- **nortest / tseries**: 残差正規性・定常性検定

## Testing

テストは `tests/testthat/` に16ファイル。スナップショットは `tests/testthat/_snaps/` に保存。

入力バリデーション、S3クラス確認、並列実行再現性、スナップショット比較の4パターンが中心。

## CI/CD

`.github/workflows/pkgdown.yaml` のみ。main/masterへのpushまたはリリース時にpkgdownサイトをgh-pagesへデプロイする。R CMD checkのCIは設定されていない。
