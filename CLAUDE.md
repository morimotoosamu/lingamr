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
| `LingamResult` | `R/lingam_direct.r` | Direct LiNGAMの推定結果（`lingam_high_dim()` も同クラスを返す） |
| `BootstrapResult` | `R/lingam_bootstrap.r` | Bootstrap安定性評価 |
| `VARLiNGAMResult` | `R/lingam_var.r` | VAR-LiNGAMの推定結果 |
| `VARBootstrapResult` | `R/lingam_var_bootstrap.r` | VARのBootstrap結果 |
| `MultiGroupLingamResult` | `R/lingam_multi_group.r` | MultiGroup Direct LiNGAMの推定結果 |
| `MultiGroupBootstrapResult` | `R/lingam_multi_group_bootstrap.r` | MultiGroupのBootstrap結果 |
| `ParceLingamResult` | `R/lingam_parce.r` | BottomUpParceLiNGAMの推定結果 |
| `RCDResult` | `R/lingam_rcd.r` | RCD（潜在交絡ありの因果探索）の推定結果 |
| `LiMResult` | `R/lingam_lim.r` | LiM（連続・離散混合データ）の推定結果 |
| `ImputationBootstrapResult` | `R/bootstrap_with_imputation.r` | 多重代入つきBootstrap結果（`as_bootstrap_result()` で `BootstrapResult` に変換可） |
| `lingam_summary` | `R/summary_lingam.r` | 包括的な適合度サマリー |

**重要な規約:** 隣接行列 `B[i,j]` は変数 j → i の因果係数を表す。

### 主要アルゴリズム

**Direct LiNGAM** (`R/lingam_direct.r`):
- `measure`: "pwling"（ペアワイズ独立）or "kernel"（カーネルベース）
- `reg_method`: "ols" / "lasso" / "adaptive_lasso"（デフォルト） / "ridge"
- 回帰バックエンドは `R/fit_regression.r`（glmnetが必要な手法あり）
- 因果順序探索は `R/search_causal_order.r`
  - `measure = "kernel"` はn > 1000でincomplete Cholesky低ランク近似（`kernel_mi_prepare_lowrank`/`kernel_mi_core_lowrank`）に自動切替、n <= 1000は既存の正確計算（`kernel_mi_prepare`/`kernel_mi_core`）のまま

**VAR-LiNGAM** (`R/lingam_var.r`):
- VARモデルを当てはめ、残差にDirect LiNGAMを適用
- ラグ選択基準: "bic"（デフォルト）/ "aic" / "hqic" / "fpe"
- `adjacency_matrices[1,,]` が瞬時構造B0、以降がラグ行列

**Bootstrap並列実行** (`R/lingam_bootstrap.r`):
- `parallel::makePSOCKcluster()` + L'Ecuyer RNGストリームで再現性を保証
- `n_cores` が変わると数値結果が変わる（設計上の仕様）

**その他の移植アルゴリズム**（いずれも Python cdt15/lingam からの移植）:
- MultiGroup Direct LiNGAM (`R/lingam_multi_group.r`): 複数データセットの共通因果順序を同時推定
- BottomUpParceLiNGAM (`R/lingam_parce.r`): 潜在交絡に頑健、未解決ブロックは `NA`
- RCD (`R/lingam_rcd.r`): 潜在交絡ペアを `NA` で表現、MLHSICR 回帰オプションあり
- LiM (`R/lingam_lim.r`): 連続・離散混合データ対応
- HighDimDirectLiNGAM (`R/lingam_high_dim.r`): 高次元向け、`LingamResult` を返す
- 独立性検定の共通基盤: HSIC (`R/hsic.r`)・F-correlation (`R/f_correlation.r`)
- ユーティリティ: `evaluate_model_fit()`（lavaan による SEM 適合度）、`bootstrap_with_imputation()`（mice による多重代入つき Bootstrap）

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
- **lavaan**: `evaluate_model_fit()` のSEM適合度評価
- **mice**: `bootstrap_with_imputation()` の多重代入
- **igraph / pcalg**: 診断・vignetteでの比較用

## Testing

テストは `tests/testthat/` に23ファイル（testthat edition 3）。`_snaps/` は使わず、
決定的な golden-value 方式（期待値をテスト内にピン留め、例: `test-lingam_var_snapshot.R`）で
数値回帰を検出する。

入力バリデーション、S3クラス確認、並列実行再現性、golden-value 比較の4パターンが中心。
オプション依存（glmnet / lavaan / mice / DiagrammeR / ggplot2 等）を使うテストは
`skip_if_not_installed()` でガードする。

## CI/CD

- `.github/workflows/R-CMD-check.yaml`: Ubuntu (devel/release/oldrel-1)・Windows・macOS で R CMD check を実行
- `.github/workflows/pkgdown.yaml`: main/masterへのpushまたはリリース時にpkgdownサイトをgh-pagesへデプロイ
