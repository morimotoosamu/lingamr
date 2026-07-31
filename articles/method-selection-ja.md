# 因果探索手法の選び方

`lingamr`
には11個の因果探索推定器が実装されている。いずれも「観測データから有向の
因果構造を復元する」という同じ目標を持つが、それぞれが基本のLiNGAMモデルの異なる
仮定を緩和したものである。この vignette は、その使い分けの指針を与える。

基本のLiNGAMモデルは以下を仮定する。

1.  **線形**な因果関係
2.  **非巡回**な因果グラフ（DAG）
3.  **非ガウス**で互いに独立な誤差項
4.  **潜在交絡変数がない**（すべての共通原因が観測されている）
5.  観測が**i.i.d.**（時間構造もグループ構造もない）
6.  変数が連続で、データ行列が完全（`NA` がない）

6つすべてが成り立つなら
[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)
を使う。パッケージの他の推定器はすべて、
このうち1つ（または2つ）の仮定を緩和するために存在する。

## 手法選択の決定ガイド

上から順に質問に答え、最初に「はい」になったところの手法を使う。

1.  **観測は時系列か？**
    - はい、時間的依存は自己回帰的 →
      **VAR-LiNGAM**（[`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md)）
    - はい、しかもラグを増やしてもVAR残差に自己相関が残る（移動平均的な撹乱）→
      **VARMA-LiNGAM**（[`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)）
2.  **因果関係が非線形である可能性が高いか？**
    - はい、かつ主要な共通原因はすべて観測されている →
      **RESIT**（[`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)）
    - はい、かつ未観測の変数があるかもしれない →
      **CAM-UV**（[`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)）
3.  **潜在交絡変数が存在するかもしれないか（線形の場合）？**
    - *因果順序のどの部分なら*信頼できるかを知りたい →
      **BottomUpParceLiNGAM**（[`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)）
    - *どの特定のペアが*交絡されているかを知りたい →
      **RCD**（[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)）
4.  **連続変数と離散変数が混在しているか？** →
    **LiM**（[`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md)。既定は二値、`is_poisson = TRUE`
    でポアソンカウント）
5.  **因果構造を共有する複数のデータセットがあるか**（複数拠点・複数期間・複数
    コホートなど）？ → **MultiGroup Direct
    LiNGAM**（[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)）
6.  **データに欠測値（`NA`）が含まれるか？** →
    **[`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md)**（ブートストラップ +
    多重代入）
7.  **変数が多い（数十〜数百）、あるいは p \> n か？** →
    **HighDimDirectLiNGAM**（[`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md)）
8.  **上のどれにも当てはまらない** → **Direct
    LiNGAM**（[`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md)）

複数の複雑さが同時に当てはまる場合（たとえば非線形の時系列）、両方を扱える単一の
推定器はない。データにとって最も重大な仮定違反を優先するか、仮定違反が少なくなる
ようにデータを変換（たとえば非定常系列の差分化）することになる。

## 手法一覧

| 手法 | 関数 | 対応する状況 | 出力の注意点 | Bootstrap |
|----|----|----|----|:--:|
| Direct LiNGAM | [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md) | 基準となる手法: 線形・非巡回・非ガウス・i.i.d. | 因果順序 + 係数行列 | [`lingam_direct_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct_bootstrap.md) |
| HighDimDirectLiNGAM | [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md) | 変数が多い。p \> n | [`lingam_direct()`](https://morimotoosamu.github.io/lingamr/reference/lingam_direct.md) と同じオブジェクト | — |
| VAR-LiNGAM | [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md) | 定常時系列（AR） | 瞬時行列B0 + ラグ行列 | [`lingam_var_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var_bootstrap.md) |
| VARMA-LiNGAM | [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md) | MA誤差を持つ時系列 | AR（psi）+ MA（omega）行列 | [`lingam_varma_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma_bootstrap.md) |
| MultiGroup Direct LiNGAM | [`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md) | 複数データセット・共通の順序 | 共通順序 + グループ別行列 | [`lingam_multi_group_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group_bootstrap.md) |
| BottomUpParceLiNGAM | [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md) | 潜在交絡変数（線形） | 未解決ブロック。`NA` エントリ | [`lingam_parce_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce_bootstrap.md) |
| RCD | [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md) | 潜在交絡変数（線形） | 祖先集合。交絡ペアは `NA` | [`lingam_rcd_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd_bootstrap.md) |
| RESIT | [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md) | 非線形加法ノイズ | 0/1エッジ（係数なし） | [`lingam_resit_bootstrap()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit_bootstrap.md) |
| CAM-UV | [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md) | 非線形 + 未観測変数 | 親リスト。UCP/UBPペアは `NA` | — |
| LiM | [`lingam_lim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_lim.md) | 連続・離散の混合データ | 係数行列 | — |
| Bootstrap with imputation | [`bootstrap_with_imputation()`](https://morimotoosamu.github.io/lingamr/reference/bootstrap_with_imputation.md) | 欠測値（`NA`） | [`as_bootstrap_result()`](https://morimotoosamu.github.io/lingamr/reference/as_bootstrap_result.md) で集約 | （それ自体がbootstrap） |

計算コストに関する実務上の注意を2点。

- **HSICベースの手法**（ParceLiNGAM、RCD、RESIT、CAM-UV、および
  `lingam_direct(measure = "kernel")`）は検定ごとに $`n \times n`$
  のグラム行列を 構築する。$`n`$
  が数千のオーダーでは推奨されない。先にサブサンプリングすること。
- **Direct LiNGAM** の計算コストは変数数について $`O(p^3)`$
  である。$`p`$ が大きい 場合は
  [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md)
  に切り替える。

各手法の実行例はパッケージサイトにある。

- [Direct LiNGAM
  詳説](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.html)
- [ブートストラップと診断](https://morimotoosamu.github.io/lingamr/articles/bootstrap-diagnostics-ja.html)
- [時系列:
  VARとVARMA](https://morimotoosamu.github.io/lingamr/articles/time-series-ja.html)
- [潜在交絡変数:
  ParceLiNGAMとRCD](https://morimotoosamu.github.io/lingamr/articles/latent-confounders-ja.html)
- [非線形の手法:
  RESITとCAM-UV](https://morimotoosamu.github.io/lingamr/articles/nonlinear-ja.html)
- [特殊なデータ:
  LiM・MultiGroup・欠測値](https://morimotoosamu.github.io/lingamr/articles/special-data-ja.html)

## LiNGAMが使えない場合

仮定が満たされない場合、推定は失敗するか、あるいは誤った構造を系統的に復元して
しまう。

| 仮定 | 問題が生じる場合 | 対処法・代替手段 |
|----|----|----|
| **非ガウス誤差** | すべての誤差がガウス分布に従う場合、因果の方向は識別不能になる | どのLiNGAM変種でも解決できない。制約ベースの手法（`pcalg` のPCアルゴリズムなど。一意な方向ではなく同値類を返す）を検討 |
| **非巡回グラフ（DAG）** | フィードバックループ（x -\> y -\> x）が存在する場合 | Cyclic LiNGAM（Python `lingam` パッケージに実装あり）の使用を検討 |
| **潜在共通原因が存在しない** | 未観測の共通原因（隠れた交絡変数）が存在する場合 | [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)・[`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)（線形）、[`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)（非線形） |
| **線形な因果関係** | 変数間の関係が非線形である場合 | [`lingam_resit()`](https://morimotoosamu.github.io/lingamr/reference/lingam_resit.md)・[`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md) |
| **測定誤差がない（上流の変数）** | rootに近い変数に大きな測定誤差がある場合、方向が系統的に逆転する | [Direct LiNGAM記事](https://morimotoosamu.github.io/lingamr/articles/direct-lingam-ja.html)の測定誤差パラドックスを参照 |
| **独立同一分布（i.i.d.）** | 時系列データ、または異質なソースを混ぜたデータ | [`lingam_var()`](https://morimotoosamu.github.io/lingamr/reference/lingam_var.md) / [`lingam_varma()`](https://morimotoosamu.github.io/lingamr/reference/lingam_varma.md)（時系列）、[`lingam_multi_group()`](https://morimotoosamu.github.io/lingamr/reference/lingam_multi_group.md)（グループ化データ） |
| **十分なサンプルサイズ** | 変数数$`p`$に対して$`n`$が極端に小さい場合（目安: $`n < 10p`$）、推定は不安定になりやすい | 変数の数を減らす。`reg_method = "adaptive_lasso"` でスパース化する。$`p > n`$ では [`lingam_high_dim()`](https://morimotoosamu.github.io/lingamr/reference/lingam_high_dim.md) を使う |

## 事前に確認すべきチェックリスト

実際の分析を始める前に、以下を確認することを推奨する。

1.  **グラフの非巡回性**：ドメインの専門知識からフィードバックループを排除できるか
2.  **潜在変数の不在**：重要な観測変数はすべて揃っているか。揃っていなければ
    [`lingam_parce()`](https://morimotoosamu.github.io/lingamr/reference/lingam_parce.md)
    /
    [`lingam_rcd()`](https://morimotoosamu.github.io/lingamr/reference/lingam_rcd.md)
    /
    [`lingam_camuv()`](https://morimotoosamu.github.io/lingamr/reference/lingam_camuv.md)
    を優先する
3.  **誤差の非ガウス性**：[`test_residual_normality()`](https://morimotoosamu.github.io/lingamr/reference/test_residual_normality.md)
    で確認できる（ただしこれは
    推定後の診断である）。事前の簡易チェックとして、各変数のヒストグラムと歪度を
    目視で確認する
4.  **測定誤差の有無**：rootに近い変数に測定誤差はあるか。ある場合は解釈に注意する
5.  **サンプルサイズ**：$`n \geq 10p`$
    を目指す。それに満たない場合、結果を過度に 信頼しない

> **まとめ:**
> LiNGAMは、線形性・非巡回性・非ガウス性・潜在変数の不在・i.i.d.という
> 仮定がすべて成り立つ場合に強力であり、`lingamr`
> はそれぞれの仮定が破れる典型的な
> 状況に対応した専用の推定器を提供している。分析前にドメイン知識と残差診断を通じて
> 仮定を検証することが、信頼できる因果推論への第一歩である。
