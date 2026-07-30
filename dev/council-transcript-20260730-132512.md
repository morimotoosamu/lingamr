# LLM Council 議事録: lingamr 未実装アルゴリズムの実装優先順位

日時: 2026-07-30

## 原質問

未実装の手法について、実装の優先順位を評議会で検討してください。

## フレーミングされた質問

Rパッケージ lingamr（Pythonの因果探索パッケージ cdt15/lingam を単独作者が移植しているもの）について、未実装アルゴリズムの実装優先順位を検討してほしい。

**背景**: lingamrはPython版lingam(v1.13.0)の一部アルゴリズムを移植済み。実装済み(7つ): DirectLiNGAM, VARLiNGAM, MultiGroupDirectLiNGAM, BottomUpParceLiNGAM, RCD, LiM(Poisson型離散変数のみ未対応), HighDimDirectLiNGAM。これらの実装では、独立性測度としてHSICとF相関(fcorr)を共有基盤として使い回しており、bootstrap版・prune処理・テストスイートも各アルゴリズムに揃っている。

**検討対象(未実装、15手法)**:
- 基本モデル系: ICALiNGAM, GroupDirectLiNGAM, GroupLiNGAM, ABICLiNGAM, LEWIS
- 時系列系: VARMALiNGAM, LongitudinalLiNGAM, LongitudinalRESIT
- 非線形・潜在変数系: RESIT, CAMUV, LiNA, MDLiNA
- マルチグループ系: MultiGroupRCD, MultiGroupCAMUV, MultiGroupRESIT
- 欠損データ系: mLiNGAM

**評価観点**: 実装難易度・工数／依存関係／Rユーザーの実用需要／CRAN公開パッケージとしての差別化

---

## アドバイザー回答(匿名化: A-E)

### The Contrarian(=応答B)

# 実装優先順位(結論)

1. **RESIT** 2. **CAMUV** 3. **ICALiNGAM** 4. **GroupDirectLiNGAM** 5. **MultiGroupRESIT** — 以下、VARMALiNGAM、LongitudinalLiNGAM、LongitudinalRESIT、MultiGroupCAMUV、GroupLiNGAM、mLiNGAM、LiNA、MultiGroupRCD、MDLiNA、ABICLiNGAM、LEWIS

## この順位が壊れる前提を潰しておく

**RESITを1位に置く議論は「基盤再利用」を過大評価している。** RESITはHSIC非依存性検定を使う点でRCD/Parceと似ているが、両者の核心は「加法ノイズモデルの回帰残差」であり、RCDのHSIC実装は独立性検定モジュールとして流用できても、GAM/非線形回帰のノンパラ回帰エンジンをRに用意する工数を過小評価しがち。Rでは`mgcv::gam`があるので致命的ではないが、「既存コード流用で軽い」という前提はここで崩れる。工数は中〜重、単独作者の"少しずつ"には荷が重い。

**ICALiNGAMを軽視するのは誤り。** 最も基本的な手法が最後まで実装されていないのは、CRANパッケージとして「網羅性の欠如」という致命傷になる。lingam論文を読んだRユーザーが最初に探すのはICALiNGAMであり、DirectLiNGAMだけでは「元論文の主要手法が抜けている」と映る。工数は`fastICA`パッケージ一本で最軽量、依存関係もゼロ。これを1位に置かないロードマップは順序が転倒している。

**差別化の観点は全体的に楽観的すぎる。** 「Rエコシステムで先行」と言えるのは各アルゴリズムが単に存在することの価値であり、Python版に対する独自性ではない。実装優先度を差別化で語るのは的外れで、素直に「工数最小×需要最大」で決めるべき。ABICLiNGAM・LEWIS・MDLiNAは論文実装の難度も需要も不透明で、単独作者ロードマップの最下位が妥当。

### The First Principles Thinker(=応答D)

前提を疑う。「15手法を難易度×需要でスコアリングして並べる」という発想そのものが罠。単独作者のロードマップで意味を持つ指標はただ一つ——「1手法の実装が次に何を可能にするか」という連鎖の起点になり得るか。

第1位: ICALiNGAM。実装ではなく欠落の異常性が問題。DirectLiNGAM以下7手法を実装しながらFastICAベースの原型が無いのは、教科書の第1章を書かずに応用編だけ出版しているに等しい。R側にfastICAが既にあり、HSIC/pruneも流用可能で工数は最小。比較ベースライン・教育用途としての需要は他の応用系より底堅い。

第2位: RESIT。ここが本当の分岐点。MultiGroupRESIT・LongitudinalRESITの前提であり、非線形回帰＋HSIC独立性検定という既存基盤をそのまま転用できる。1つ実装すれば2つの未実装手法の工数が激減する「レバレッジ資産」。RESITを飛ばして他の非線形系(CAMUV/LiNA)に手を出すのは資産を使い捨てる愚策。

第3位: CAMUV。RCDのHSIC基盤をほぼそのまま使え、RESITと合わせて非線形因果探索の需要(医療・行動データ)を最速で満たす。

第4位: VARMALiNGAM。VARLiNGAM資産の直接延長で工数対効果が高い。

以降: LongitudinalLiNGAM→MultiGroupRESIT/CAMUV→GroupDirectLiNGAM/GroupLiNGAM→LongitudinalRESIT→MultiGroupRCD→LiNA/MDLiNA→LEWIS。

明確に後回しにすべきは2つ。mLiNGAMはbootstrap_with_imputation()と目的が重複し差別化を生まない。ABICLiNGAMは深層学習依存でPython橋渡し(reticulate等)が必須になり、「軽量・CRAN準拠」というRエコシステムでの差別化軸そのものを破壊する。むしろ「深層学習に手を出さない」ことこそがPython版に対する差別化戦略であり、実装しないという判断が最上位の意思決定であるべき。

### The Expansionist(=応答A)

15手法の実装優先順位。単独作者の逐次実装ロードマップとして提示する。

**上位ランク**

1. **RESIT** — 最優先。非線形加法ノイズモデルという「非線形因果探索」枠の入口であり、かつMultiGroupRESIT・LongitudinalRESITの土台になる。RCDで独立性検定(HSIC)・回帰残差評価の基盤はすでにある。1つ実装すれば後続2手法の限界コストが激減する、レバレッジ最大の一手。
2. **CAMUV** — RESITと並走可能。潜在共通原因への対処はマーケティング・医療データで「未観測交絡」への実務ニーズが直撃する。ParceのHSIC基盤を流用でき、Python版でもRESITと組コードなので同時着手が効率的。
3. **MultiGroupRESIT** — RESIT完成後は増分コスト小。MultiGroupDirectLiNGAMで多群比較の設計パターン(グループ間差分検定・共通構造推定)が確立済みなので、RESITの回帰残差評価にそれを差し込むだけ。医療(多施設データ)・A/Bテスト応用で差別化になる。
4. **ICALiNGAM** — 実装は軽い(FastICAはRのfastICAパッケージに直接依存できる)くせに、CRAN公開パッケージとしての「網羅性」の見栄えが強い。DirectLiNGAMより前の基本形が抜けているのは学術ユーザーへの説得力を欠く穴。工数対効果で見れば最速リターン案件。
5. **VARMALiNGAM** — VARLiNGAMがすでにあるので時系列系の次の一手として自然。マーケティングのMMM(マーケティングミックスモデリング)や経済時系列でVARでは不十分な場面の需要は厚い。

**下位に回すもの**: LongitudinalLiNGAM/RESIT、MultiGroupRCD/CAMUV、GroupDirectLiNGAM/GroupLiNGAM、LiNA/MDLiNA、ABICLiNGAM、LEWIS、mLiNGAMは、依存基盤ができてからの「後乗せ」で十分。特にABICLiNGAM(微分可能因果探索)はTorch的自動微分基盤がR側になくゼロから作る羽目になるため最後回し。mLiNGAMはbootstrap_with_imputation()で代替済みなので実装優先度は最低。

RESIT一本で後続3手法の実装コストが激減する構造こそが最大の狙い目。ここに全力集中すべき。

### The Outsider(=応答C)

# 優先順位評価（アウトサイダー視点）

## 最優先：RESIT → LongitudinalLiNGAM → MultiGroupRESIT
**理由**：ハブ・アンド・スポーク戦略。RESIT は非線形因果探索の基幹で、これ以降の Longitudinal/MultiGroup 拡張を解き放つ。R には縦断データの因果推定が実装された信頼できるツールが *ない*。医療・社会科学は「時系列 + 反復測定」に飢えている。Python にはあっても、R ユーザーは既存のパッケージ（lavaan の潜在増殖度は限界、bnlearn は構造仮定が強い）に妥協している。これは差別化点。

## 次点：VARMALiNGAM
既に VAR-LiNGAM があるなら、ARMA 拡張は低投資で経済時系列・金融の実用価値が劇的に上がる。

## 検討段階：ICALiNGAM, MultiGroupRCD
ICA は「基本」だが、直接法・多変量対応が既にあるなら後回し。MultiGroupRCD は既存 RCD を流用できるはず。ROI は低～中。

## 後回し：GroupDirectLiNGAM, ABICLiNGAM, LEWIS, MDLiNA, LiNA
- GroupDirectLiNGAM → MultiGroup 系で吸収される可能性
- ABIC → 微分最適化は保守負荷高で、VARLiNGAM で十分なことが多い
- LEWIS・MDLiNA・LiNA → 引用数・利用事例が明らかに少ない

**結論**：単独作者なら「縦断因果推定の R における唯一の信頼できる選択肢」に特化せよ。RESIT 系で市場を獲得する方が、15 手法全制覇より現実的で持続可能。

### The Executor(=応答E)

# 優先実装順位

**1位: ICALiNGAM**
最基本。既存 DirectLiNGAM と独立性測度を完全共有。工数 1-2 週。CRAN での「FastICA 系統」完全カバレッジが出来、Rユーザーにとって基本ツールとしての信頼性が高まる。差別化ではなく「必須」枠。

**2位: RESIT**
非線形への唯一の橋。MultiGroupRESIT・LongitudinalRESIT の上流。RCD・BottomUpParce と同じ独立性測度基盤を再利用できる相乗効果。時系列・マルチグループへの展開が見える。工数 3-4 週。応用広い。

**3位: GroupDirectLiNGAM**
既実装 MultiGroupDirectLiNGAM の親概念。段階的マルチグループ対応の進捗を示せる。工数 2-3 週。DirectLiNGAM の拡張で明快。

**4位: VARMALiNGAM**
VAR-LiNGAM 実装済みなので拡張工数中。時系列分析は Rユーザー（金融・マーケ）の実用需要が高い。工数 3-4 週。

**5位: MultiGroupRCD**
RCD 拡張。工数中。段階的に MultiGroup 系統を充実させる布石。

---

**避けるべき：** mLiNGAM（既に bootstrap_with_imputation で代替）、CAMUV・LiNA・MDLiNA（複雑度高い割にニッチ）、Longitudinal 系（多変量時系列の需要が限定的）。

**現実的 2 ヶ月ロード：** ICALiNGAM → RESIT → GroupDirectLiNGAM で「基本・非線形・マルチグループ」の軸足を固める。その後、時系列か潜在変数系かユースケースで判断。

---

## 匿名化マッピング

| 匿名ラベル | アドバイザー |
|---|---|
| A | The Expansionist |
| B | The Contrarian |
| C | The Outsider |
| D | The First Principles Thinker |
| E | The Executor |

---

## ピアレビュー(5件)

**レビュー1**: 最強はA。RESIT一択集中の論理(HSIC/RCD基盤流用→MultiGroupRESIT/LongitudinalRESIT連鎖)が明快で、mLiNGAM・ABICLiNGAMの後回し根拠も具体的。工数の楽観バイアスをBが指摘する分、AとBを合わせて読むと補完性が高い。最大の盲点はC。「Rにない」「lavaan/bnlearnは限界」という差別化主張に、実際の代替パッケージ(CausalImpact、tsDyn、networkDynamic等)との比較検証がなく、需要断定が根拠薄弱。中国語混入("纵断")も推敲不足の兆候。全員が見落とした点: LiM未対応のPoisson型離散変数の補完という「安い勝ち筋」に誰も触れていない。単独作者の燃え尽き/継続性リスクへの言及がBのみ。CRAN上でのメンテナンス負荷を定量的に検討した回答がない。Python版lingamでの各手法の実際のGitHub star/issue/citation数などの需要の一次データが全回答とも欠如し、印象論に留まっている。

**レビュー2**: 最も強いのはD。ICALiNGAM単独実装の妥当性を「原論文の第1章が抜けている」という具体的比喩で説得力高く示し、ABICLiNGAMがreticulate依存でR版の「軽量・CRAN準拠」という差別化軸自体を壊すという鋭い指摘を含む。最大の盲点はA・C。Aは「RESIT実装後の後続コスト激減」を主張するが、RESITの非線形回帰エンジン(GAMなど)の工数を過小評価。Cは縦断データの差別化価値を語るが、LongitudinalLiNGAMの実装難度が不透明で検証不足。全員が見落としたこと: テストスイート・診断ツールの実装コスト。既存手法との「テスト基盤の共有可能性」の詳細が未検証。手法数増加に伴うメンテナンス・ユーザーサポート負荷も定性的評価のみ。

**レビュー3**: 最強はB。Aの「基盤再利用で軽い」という前提を具体的に破壊(RESITはGAM/非線形回帰エンジンが必要で工数「中～重」)。ICALiNGAMの欠落が「教科書第1章を書かずに応用編」という本質的な問題と認識。最大のブラインドスポットはA。「既存コード流用で軽い」という楽観が全体を支配。RESITのGAM実装工数・各手法のテスト・CRANリリースサイクルを過小評価。全員が見落とした点: テストスイート・検証の累積工数、既実装7手法の保守・バグ対応(依存パッケージのメンテ終了リスクも未言及)、「実装しない」という戦略判断(工数対効果がマイナスの手法を明示的に「実装不要枠」として切る判断が抜けている)。

**レビュー4**: 最強はA。RESIT一位という結論の論拠(RCDのHSIC・回帰残差基盤を明示的に指差し、後続2手法への波及を定量的に語る)が具体的。CAMUV・MultiGroupRESITの接続もPython版の実装関係まで踏み込んでおり裏取りが厚い。最大の盲点はC。「LongitudinalLiNGAM」をRESITの直後・第2位に置くが、LongitudinalLiNGAMは非線形加法ノイズモデルでなくVAR-LiNGAM系の縦断拡張であり、RESIT基盤を再利用できる根拠が薄い。全員が見落としたのは検証コスト。シミュレーション検証・数値安定性確認・既存Python版とのユニットテスト突合の工数を見積もっていない。mgcv依存やfastICA依存が既存DESCRIPTIONのCRAN依存関係ポリシーに与える影響にも触れていない。

**レビュー5**: 最も強いのはD。ICALiNGAMの「欠落の異常性」という本質的指摘と、ABICLiNGAMは「CRAN軽量性という差別化軸を破壊する」という後ろ向き判断が秀逸。単独作者の持続可能性を中心に据えた思考。最大の盲点はA・C。Aは「RESIT実装後の後続コスト激減」を主張するがGAMなどの工数を過小評価。Cは縦断データの差別化価値を語るが実装難度の検証が不足。全員が見落としたこと: テストスイート・診断ツールの実装コスト、手法数増加に伴うメンテナンス・ユーザーサポート負荷も定性的評価のみ。

---

## 議長の統合見解

## Where the Council Agrees

- RESITを上位に置くことに5人中4人(A/B/D/Eの順位、Cも実質同格)が一致。理由も共通: HSIC独立性検定の資産(RCD/Parce由来)が使え、MultiGroupRESIT・LongitudinalRESITの前提になるという「レバレッジ資産」論。
- mLiNGAMは全員一致で最下位級。bootstrap_with_imputation()と目的が重複し、新規実装の差別化効果がないという評価に異論なし。
- ABICLiNGAM・LEWISは軒並み最下位圏。深層学習依存・微分最適化基盤・引用実績の乏しさが理由。Dは「実装しないこと自体が差別化戦略」と踏み込んで明言し、他の回答も暗黙にこれを支持している。
- CAMUVはRESITと並走できる、というRCD/Parce基盤の転用可能性についても異論がない。

## Where the Council Clashes

- **ICALiNGAMの順位が真っ二つ**。D・Bは「最も基本的な手法が抜けているのはCRANパッケージとして致命的」として1位に置く。A・Eは「工数最速回収」と自ら評価しながら4位・1位に分かれて置いており、特にAは工数対効果最速と言いながら4位に落とす論理的矛盾をピアレビューで指摘されている。ここは未解決のまま。
- **RESITの「軽さ」の評価が対立**。A/C/D/Eは既存のHSIC基盤流用を根拠に「レバレッジ大・軽い」と評価するが、Bだけが「RESITの核心は加法ノイズモデルの非線形回帰(GAM)であり、HSICが使えるのは独立性検定部分のみ。mgcv::gamの実装・検証工数は別途必要で工数は中〜重」と技術的に反論している。この反論は複数のピアレビューで「Aの最大の盲点」として支持されており、Bの指摘の方が信頼度が高い。
- Cの「RESIT→LongitudinalLiNGAM」という連鎖の主張は根拠薄弱。LongitudinalLiNGAMは線形VAR系の縦断拡張であり、RESIT(非線形加法ノイズモデル)の基盤が転用できる理由が説明されていない。複数のピアレビューが同じ点を指摘しており、Cの依存関係論は棄却してよい。
- 「差別化」を判断軸に使うべきかどうかも対立。A/Cは差別化(Rエコシステムでの先行)を重視、Bは「差別化ではなく工数×需要で決めるべき」と真っ向から否定。

## Blind Spots the Council Caught

- **LiMのPoisson型離散変数未対応という「安い勝ち筋」に5人とも触れていない**。既存7手法のうち1つが未完成なのに、15の新規手法着手を検討する前にこの穴を塞ぐ方が低コストで完成度も上がる。これは今回のピアレビューが唯一指摘した、最も重い見落とし。
- 単独作者の継続可能性(燃え尽きリスク)への言及がBに薄くあるのみで、他4人はほぼ無視。工数を「2-3週」「3-4週」と軽く積み上げているが、それが積み重なったときの持続可能性は誰も検証していない。
- CRAN上の保守負荷(mgcv・fastICA等の依存追加によるcheck負担、依存パッケージの終了リスク)は誰も定量的に触れていない。
- 全回答が「HSIC/fcorr基盤の再利用可能性」ばかりを語り、シミュレーション検証・数値安定性確認・Python版との出力突合というテスト工数を過小評価している。単独作者の実工数はコーディングよりデバッグ・検証に割かれるという構造的な事実が全員から抜けている。
- Python版lingamでの各手法のGitHub star/issue/citation数など、需要を裏付ける一次データが誰からも提示されていない。需要論はすべて印象論。

## The Recommendation

優先順位(現実的ロードマップ):

1. **LiMのPoisson型離散変数対応**(15手法の枠外だが最優先) — 既存実装の完成度を上げる方が新規着手よりROIが高い、全員が見落とした安い勝ち筋。
2. **RESIT** — ただしBの指摘通り「工数中〜重」と見積もり直した上で着手。GAM/非線形回帰エンジンの検証コストを事前に別途見積もること。MultiGroupRESIT/LongitudinalRESITへの波及を狙う本命。
3. **ICALiNGAM** — 網羅性の穴を塞ぐ実利が大きく、fastICA依存のみで工数最小。RESITと並行、または直後に着手して「基本形の欠落」を早期に解消する。
4. **CAMUV** — RESITと同じHSIC基盤で並走可能。潜在共通原因への対応は医療・マーケティング需要と直結。
5. **VARMALiNGAM** — VARLiNGAM資産の直接延長。工数対効果が高く、時系列需要も厚い。

明示的に実装不要枠とする: mLiNGAM(bootstrap_with_imputationと重複)、ABICLiNGAM(reticulate依存がCRAN軽量性という差別化軸を破壊)、LEWIS・MDLiNA(需要不透明)。Cの「縦断因果推定特化」戦略は依存関係の裏付けが弱く、単独の戦略としては採用しない。

## The One Thing to Do First

15手法のどれにも着手する前に、まず(1) LiMのPoisson対応を仕上げ、(2) Python版lingamの各手法についてGitHub star/issue数・引用数を1時間程度でよいので実際に調べること。5人の助言全員が需要判断を印象論で行っており、これがロードマップ全体の信頼性を損なっている最大の欠陥。データを見てからRESIT着手の是非を最終確定すべき。
