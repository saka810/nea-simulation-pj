# geosim — プログラム構成まとめ

幾何音響シミュレーション（音線法 + 虚音源バックトレース）の Fortran → Python 移植プロジェクト。
元コードは `fortran/` の 3 ファイル（`backtrace.f90` / `make_ipls_freq_monaural_fortran.f90` / `ipls2rt_fortran.f90`）。
最終更新: 2026-08-23（HBN240018）。**コードに手を入れたらこの文書も合わせて更新する。**

## 実行環境

**Python 3.10.11**（チーム方針）。依存は `requirements.txt`
（numpy / scipy / matplotlib / pyvista / **openpyxl**）。
構築手順は [README.md](README.md#環境構築) にある。2026-08-14 までは numpy すら
インストール前提にしない方針だったが、GUI 表示のためにこの日 venv を切って導入した。

## 処理パイプライン全体像

```
入力: 音源位置・受音点・DXF形状・受音球半径・反射回数・音線数・吸音率リスト

procedure.process()
 ├─ 1. read_dxffile.read_model()      … DXF → Meshオブジェクトのリスト
 ├─ 2. sound_ray.soundray_generator() … 球面上に均等な音線ベクトル群を生成
 ├─ 3. loop_reflectionmesh.loop()     … 音線追跡。受音した経路の反射面ID履歴を記録
 │                                       （元コード backtrace.f90 524行〜）
 ├─ 4. loop_deleteredundancy.delete() … 反射履歴の重複経路を削除（元コード 721行〜）
 ├─ 5. loop_noredundancy.loop()       … 虚音源法バックトレース。到来時刻・到来方向・
 │                                       バンド別エネルギーのパルス列（元コード 876行〜）
 ├─ 6. impulse.impulse_response_from_pulses()     … パルス列 → バンド合成でインパルス応答CSV
 │                                       （元コード make_ipls_freq_monaural_fortran.f90）
 ├─ 7. reverberation.reverberation_time() … 減衰曲線と残響指標 EDT / T20 / T30
 │                                       （元コード ipls2rt_fortran.f90）
 ├─ 7b. reverberation.statistical_reverberation() … 統計残響式の理論値
 │                                       （Sabine / Eyring / Eyring-Knudsen）
 └─ 8. sound_level.band_levels() / speech_transmission_index()
                                      … 帯域別の音圧レベルと STI（元コードにない追加）

2026-08-14 に 5〜7 を実装して全段つながった。8 は 2026-08-21 に追加。
```

**受音点が複数あるときの回し方**（`run_project.run()`。2026-08-21）

```
受音点に依らないもの（1 回だけ）              受音点ごと（点の数だけ）
 ├─ 1. DXF 読み込み                          ├─ 4. 重複経路の削除
 ├─ 2. 音線生成                              ├─ 5. バックトレース
 ├─ 3. 音線追跡（受音球を全点ぶん置く）      ├─ 6. インパルス応答
 └─ 7b. 統計残響式（1 点目の結果を配る）      ├─ 7. 残響指標
                                              └─ 8. 音圧レベル・STI
```

★**受音点に依らない計算を受音点ごとに回さない**（2026-08-21 ユーザー指摘）。
結果の置き場もそれに合わせてあり、受音点に依らないものは `結果/` 直下、
受音点ごとのものは `結果/recN/`。

★**経路の幾何（反射面の並びと入射角）は吸音に依らないので保存して使い回す**
（`path_cache.py`。F-9）。吸音材だけ変えた再計算は 5 の**エネルギー計算だけ**で済み、
1〜4 を丸ごと省ける。指紋（頂点・法線・パッチの分け方・音源・受音点・音線数・
反射回数・受音球）を突き合わせ、食い違えば理由を告げてやり直す。

**3 と 5 の役割分担**（この設計の要）

音線法（3）は「どの壁をどの順に反射する経路がありそうか」を**発見**するだけ。
受音球に入ったかどうかで判定するので、経路そのものには誤差がある。
虚音源法（5）は経路が分かっていれば到来時刻とエネルギーを**厳密**に出せるが、
経路の候補を自力で探すことはできない（反射回数が増えると組み合わせが爆発する）。
そこで「探索は音線法、計算は虚音源法」と分担している。
5 は受け取った経路を幾何的に検算し直すので、音線法が拾った偽の経路はここで却下される
（実測: 972 本 → 582 本が生き残り、390 本は却下）。

## ファイル一覧

| ファイル | 役割 | 対応する元コード | 状態 |
|---|---|---|---|
| `procedure.py` | 全体フローのオーケストレーション（受音点 1 点ぶん） | backtrace.f90 全体 | **実装済**（2026-08-14 に 1〜7 が全部つながった） |
| `mesh.py` | 三角形メッシュのデータクラス | vwpt/surf/vwnm/absper 配列群 | 実装済（2026-08-12 に vertexes の形状を修正） |
| `mesh_method.py` | 音線と面の交差判定の幾何計算 | 376〜447 行ほか | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `sound_ray.py` | 音線の生成・正規化・反射・エネルギー減衰 | 318〜326, 493〜503, 1091〜1110 行 | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `receiver_sphere.py` | 受音球の通過判定 | 649〜663 行 | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `read_dxffile.py` | DXF ファイルの読み込み | 132〜283 行 | **実装済**（2026-08-12。実 DXF で動作確認済） |
| `loop_reflectionmesh.py` | 音線追跡ループ本体 | 524〜717 行 | **実装済**（2026-08-12 に traceff の扱いを確定し全面修正・動作確認済） |
| `ray_recorder.py` | 可視化用の音線軌跡レコーダ | （移植元なし・新規） | 実装済（2026-08-12 新設） |
| `view_model.py` | モデルの 3D ビューア（HTML + WebGL を書き出す） | （移植元なし・新規） | 実装済（2026-08-14 新設） |
| `view_model_gui.py` | モデルの 3D ビューア（PyVista のネイティブウィンドウ） | （移植元なし・新規） | 実装済（2026-08-14 新設・描画確認済） |
| `view_rays.py` | 音線・音粒子の可視化（G-1 / G-2） | （移植元なし・新規） | **実装済**（2026-08-14 新設） |
| `view_images.py` | **虚音源の可視化**（G-31）。パルス列から復元するので計算しない | （移植元なし・新規） | **実装済**（2026-08-24 新設） |
| `source_placement.py` | **面の上に置いた音源**（半球放射・載っている面の 1 次反射なし） | （移植元なし・新規） | **実装済**（2026-08-24 新設） |
| `absorption_curve.py` | 細かい周波数の吸音率を**オクターブバンドに直す**（対数補間＋帯域平均） | （移植元なし・新規） | **実装済**（2026-08-24 新設） |
| `report_inverse_square.py` | **逆二乗の検討報告書**（Markdown ＋ PDF。サンプルの体裁に合わせる） | （移植元なし・新規） | **実装済**（2026-08-26 新設） |
| `mode_shape.py` | **音圧分布（モード形状）**を断面で塗る。直方体は虚音源を式で並べる | （移植元なし・新規） | **実装済**（2026-08-26 新設） |
| `frequency_response.py` | **伝達関数（周波数特性）と固有周波数**。干渉が見える形で出す | （移植元なし・新規） | **実装済**（2026-08-26 新設） |
| `inverse_square.py` | **逆二乗則からのずれ**（測線ごと。別ルートの出力） | （移植元なし・新規） | **実装済**（2026-08-24 新設） |
| `view_camera.py` | **画角（見る向き）のセーブ・ロード**。条件を変えても同じ向きで見比べるため | （移植元なし・新規） | **実装済**（2026-08-24 新設） |
| `atmosphere.py` | 温度・湿度・気圧 → 音速・空気吸収 | c0 定数 / mair 近似式 | **実装済**（2026-08-14 新設。ISO 9613-1） |
| `absorption.py` | 吸音材ライブラリと吸音率の種類の変換 | absper 配列 | **実装済**（2026-08-14 新設。Paris の式・GUI 対応） |
| `loop_deleteredundancy.py` | 重複経路の削除 | 721〜841 行 | 実装済（2026-08-12 に整理・動作確認済） |
| `loop_noredundancy.py` | バックトレースループ（虚音源法） | 876〜1134 行 | **実装済**（2026-08-14 全面書き直し・解析解と一致） |
| `reverberation.py` | 残響時間・減衰曲線 | ipls2rt_fortran.f90 | **実装済**（2026-08-14 新設・既知の減衰と一致） |
| `impulse.py` | インパルス応答の合成・CSV 出力 | make_ipls_freq_monaural_fortran.f90 | **実装済**（2026-08-14 全面書き直し・scipy 化） |
| `project.py` | プロジェクト（条件・法線指定・結果）の保存と読み込み | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `table.py` | 表の並べ方の共通ルール（**周波数は横**）と CSV の読み書き | （移植元なし・新規） | **実装済**（2026-08-17 新設、2026-08-21 に「区分付きの表」を追加） |
| `sound_level.py` | 帯域別の音圧レベル（絶対値・逆二乗の物差し）と **STI** | （移植元なし・新規） | **実装済**（2026-08-21 新設。IEC 60268-16） |
| `condition_table.py` | 条件表（レイヤー名 → 材料番号 ＋ 吸音率シートの xlsx） | （移植元なし・新規） | **実装済**（2026-08-21 新設） |
| `path_cache.py` | 経路の幾何を保存して**吸音材だけ変えた再計算**に使う（F-9） | （移植元なし・新規） | **実装済**（2026-08-21 新設） |
| `workbook.py` | 結果一式の Excel 出力（グラフ付き・雛形対応） | （移植元なし・新規） | **実装済**（2026-08-21 新設。openpyxl） |
| `summary.py` | 受音点をまたいだまとめ表 | （移植元なし・新規） | **実装済**（2026-08-21 新設） |
| `ray_filter.py` | 注目したい音線の絞り込み（近く・方向・1 本） | （移植元なし・新規） | **実装済**（2026-08-19 新設） |
| `plots.py` | 結果を PNG にする（G-3 / G-4） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `setup_window.py` | 計算条件の入力ウィンドウ（tkinter） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `face_editor.py` | 面の確認・修正ウィンドウ（法線・吸音材。G-8 / G-9） | （移植元なし・新規） | **実装済**（2026-08-15 新設、2026-08-19 に `normal_editor.py` から改名して吸音材を追加） |
| `run_project.py` | プロジェクトの条件で計算を回す薄い層 | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `app.py` | 入口（条件入力 → 法線確認 → 計算 → 可視化） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `progress_window.py` | 計算中の進捗ウィンドウ（別スレッド＋ログ表示） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `view_directions.py` | 音線がどの向きへ飛ぶかのプレビュー（室形状と無関係） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `__init__.py` | パッケージ化 | — | 空 |

備考: 各モジュールは `from mesh import Mesh` のような**ベア名インポート**なので、
現状は `geosim/` ディレクトリを cwd にした直接実行のみ想定（パッケージとして import するなら相対インポート化が必要）。

---

## 変数の命名対応（Fortran → Python）

| Fortran | Python | 意味 |
|---|---|---|
| `vray(3)` | `sound_ray` | 音線ベクトル（反射後も同名。必ずしも音源発ではない） |
| `vini(3)` | `soundray_comesfrom` | 音線の基点（未反射なら音源位置） |
| `node(3)` | `node` | 音線と面の交点 |
| `vwnm(j,:)` | `mesh[j].normal` | 面 j の単位法線 |
| `vwpt(surf(j,2..4),:)` | `mesh[j].vertexes` | 面 j の頂点座標 3 点 |
| `absper(id,1..6)` | `mesh[j].absorption_coefficient` | 吸音率（周波数バンド別） |
| `tractmp(0:nref+1)` | `reflectionmeshid_history` | 一時反射履歴（先頭は -1） |
| `traceff(count, 0:nref)` | `reflectionmeshid_history_2dim` | 受音した経路の反射履歴（2次元） |
| `tracred` | `reflection_history`（delete() の戻り値） | 重複除去後の反射履歴 |
| `isrc(0:nref, 3)` | `imaginary_sourcepoint` | 虚音源座標リスト（0 番目は音源） |
| `disttmp = 1.0d50` | `min_distance = np.inf` | 最寄り壁距離の初期値 |
| `enertmp(6)` | `soundray_energy` | 音線エネルギー（バンド別） |
| `rtime` | `reflection_timing` | 受音時刻（距離 / 音速） |

---

## 各ファイルの詳細

### procedure.py

全体フローの記述。エントリポイント。**1〜7 が全部つながっている**（2026-08-14）。

```python
process(soundsource_point, receiver_point, dxf_filename, sphere_radius, nref, soundray_number,
        absorption_csv=None, unit=None, orient_normals="cad",
        raylog_filename=None, raylog_max_rays=2000, sound_velocity=340.0,
        pulse_filename=None, impulse_filename=None,
        sampling_frequency=44100.0, max_time=1.0, compensate_filter_delay=False,
        reverberation_filename=None, decay_filename=None)
    -> {"model", "pulses", "impulse", "reverberation"}
```

| 引数 | 型 | 意味 |
|---|---|---|
| `soundsource_point` | (3,) \| None | 音源座標 [m]。None なら DXF の src レイヤから |
| `receiver_point` | (3,) \| None | 受音点座標 [m]。None なら DXF の rec レイヤから |
| `dxf_filename` | str | 形状 DXF ファイルパス |
| `sphere_radius` | float | 受音球半径 [m] |
| `nref` | int | 最大反射回数 |
| `soundray_number` | int | 音線数 |
| `sound_velocity` | float | 音速 [m/s]。既定 340（元コード c0）。**軌跡・バックトレース・空気吸収で共用** |
| `*_filename` | str \| None | 各段の出力先。None なら**その段の計算自体を行わない** |

- オクターブバンドは `impulse.OCTAVE_BAND_FREQUENCIES` = **125/250/500/1k/2k/4k の 6 バンド**。
  以前は 63〜8000 Hz の 8 バンドを定義していたが、吸音率もバックトレースも 6 バンドなので
  不整合だった。6 バンドである根拠は make_ipls の 6→32 展開表（94〜125 行）。
- `nref` に達した経路があると警告を出す（後部残響が切れている＝残響時間が信用できない合図）。
- `main()` があり、コマンドラインから通し実行できる。

```
cd geosim
python procedure.py ..\test.dxf --absorption ..\absorption.csv --out ..\結果 ^
       --rays 20000 --nref 8 --compensate-delay
```

### atmosphere.py

温度・湿度・気圧から**音速と空気吸収**を求める。移植元のない**新規モジュール**（2026-08-14）。

```python
Atmosphere(temperature=20.0, humidity=40.0, pressure=101.325)
    .sound_velocity                 # [m/s]
    .absorption_coefficient(f)      # [1/m]（E = E0 exp(-m d) の m）
    .absorption_db_per_metre(f)     # [dB/m]
    .replace(temperature=25.0)      # 一部だけ変えた新しい Atmosphere（GUI 用）
```

元コードは音速を `c0 = 340.0` の定数、空気吸収を `1.81e-8 * f^1.57` のべき乗近似で
固定していた。どちらも温度・湿度から計算する形にした。

**音速**は湿り空気の状態方程式 `c = sqrt(γRT/M)` から。
水蒸気（18.0 g/mol）は乾燥空気（28.96 g/mol）より軽いので、湿度が上がると
平均モル質量が下がって音速は**速くなる**。温度の効果のほうが大きく、20℃ 付近で 1℃ ≒ 0.6 m/s。

| 条件 | 音速 | 検証 |
|---|---|---|
| 0℃ / 0% | 331.39 m/s | 標準値 331.3〜331.5 ✓ |
| 20℃ / 0% | 343.31 m/s | 標準値 343.2 ✓ |
| **20℃ / 40%（基準）** | **343.81 m/s** | |
| 14℃ / 40% | 340.12 m/s | **元コードの 340.0 はおよそ 14℃ 相当** |

**空気吸収**は **ISO 9613-1**。古典吸収＋回転緩和（第 1 項）と、
酸素・窒素の**振動緩和**（第 2・3 項）の和で、後者は湿度に強く依存する。

| 20℃ 70%RH | 1 kHz | 4 kHz |
|---|---|---|
| 本実装 | 4.98 dB/km | 23.1 dB/km |
| ISO 9613-1 の代表値 | 約 4.8 | 約 26 |

元コードの近似は 1k〜8k Hz では ISO と ±20% 程度で合うが、**低音側は 1/3 程度しかない**。

`legacy_absorption_coefficient()` に元の近似式を比較用に残してある。

### absorption.py

吸音材のデータと、**吸音率の種類の変換**。移植元のない**新規モジュール**（2026-08-14）。

```python
# 吸音率の種類の変換（Paris の式）
statistical_absorption(z)              # z → 残響室法吸音率（数値積分）
statistical_absorption_closed_form(z)  # 同上（解析解。検算用）
impedance_from_statistical(alpha_s)    # 残響室法 → z（二分法。上限で丸める）
random_to_normal(alpha_s)              # 残響室法 → 垂直入射  ★実務で使う向き
normal_to_random(alpha_n)              # 垂直入射 → 残響室法

# 材料
Material(name, coefficients, kind='normal'|'random', note='')
MaterialLibrary.from_csv(path, kind=None) / .add() / .update() / .remove()
                .to_csv() / .to_json() / .absorption_table(assignment, band_number)
LayerAssignment({レイヤ名: 材料名}) / .assign() / .save() / .load()
```

**① 吸音率の種類を取り違えない**

| 種類 | 測り方 | 特徴 |
|---|---|---|
| 垂直入射吸音率 α_n | 音響管 | 正面から当てた場合。1 を超えない |
| 残響室法吸音率 α_s | 残響室 | **1 を超えることがある**（試料端部の回折等） |

反射計算（書籍 式2.64）が要求するのは**垂直入射**。残響室法の値をそのまま入れると
`sqrt(1 - α)` が NaN になる。Paris の式

$$\alpha_s = 2\int_0^{\pi/2}\alpha(\theta)\sin\theta\cos\theta\,d\theta
          = \frac{8}{z^2}\Big[(1+z) - 2\ln(1+z) - \frac{1}{1+z}\Big]$$

を逆に解いて変換する。**α_s は z について単調ではなく**、z ≈ 1.567 で最大値 **0.951** を取り、
両側で 0 に近づく。したがって 1 つの α_s に z は 2 つあるので、建築材料が乗る
**z ≥ z_max の枝**を選ぶ。0.951 を超える値は上限に丸めて警告する。

| 残響室法 α_s | z | 垂直入射 α_n |
|---|---|---|
| 0.20 | 32.6 | 0.116 |
| 0.50 | 9.66 | **0.340** |
| 0.80 | 3.88 | 0.652 |

**残響室法の値をそのまま垂直入射として使うと吸音を大幅に過大評価する**（0.50 → 0.34）。

**② 材料を CAD から切り離す（`LayerAssignment`）**

CAD のレイヤ名で材料を判別するが、条件を変えるたびに CAD を編集するのは大変。
そこで「レイヤ → 材料名」の対応を**外に持つ**。GUI からはこの対応表と
`MaterialLibrary` を編集すればよく、DXF は触らない。JSON で保存できるので条件を使い分けられる。

**③ Excel に無い材料を後から足せる（`MaterialLibrary.add()`）**

**バンド数**は可変。既定 8（63〜8k Hz）、6（125〜4k Hz）も可。
表と計算でバンド数が違う場合は端の値を複製／切り詰めして合わせ、注意を表示する。

### mesh.py

```python
class Mesh:
    def __init__(self, vertex_1, vertex_2, vertex_3, normal, material, absorption_coefficient)
```

| 属性 | 型（想定） | 意味 |
|---|---|---|
| `vertexes` | ndarray (3,3) | 頂点 3 点の座標。~~`np.array(([v1],[v2],[v3]))` で (3,1,3)~~ → 2026-08-12 に `np.array([v1,v2,v3])` に修正し (3,3) に |
| `normal` | ndarray (3,) | 単位法線ベクトル（CAD になければ別途計算予定） |
| `material` | str | 吸音材の名前 |
| `absorption_coefficient` | ndarray | バンド別吸音率 |

設計メモ（コード内コメントより）: 三角形メッシュに限定。データのみ保持しメソッドは `mesh_method.py` に分離する方針。

### mesh_method.py

音線と三角形面の交差判定。すべて純関数。

```python
collision_distance(sound_ray, soundray_comesfrom, normal, vertexes) -> (collision: bool, distance: float)
```
面との衝突判定と基点から交点までの距離を返す。内部で `parameter_t` → `node_renew` → `collision_detection` を順に呼ぶ。
法線との内積が負（面に向かう）かつ t > 0（前方）の場合のみ判定。

```python
collision_detection(node, vertexes) -> bool
```
交点が三角形の内側かを、基準頂点を変えた 2 回の `innerproduct_from3vertexes` の符号で判定
（両方 0 以下なら内側）。2026-08-12 に `< 0` → `<= 0` へ修正し元コード 625 行と一致させた。

```python
innerproduct_from3vertexes(node, vertex_origin, vertex_1, vertex_2) -> float
```
基準頂点から面を張る 2 辺ベクトルと交点ベクトルの外積 2 つを作り、その内積を返す（内外判定の素材）。
※2026-08-12 に作業配列の形状を `np.zeros((3, 2))` → `np.zeros((2, 3))` に修正（旧形状では代入時に ValueError）。
なお 2 回目の呼び出しは元コードと外積の引数順が入れ替わっているが、内積は可換なので結果は同値。

```python
parameter_d(normal, vertex) -> float          # 平面方程式 ax+by+cz+d=0 の d = -normal・vertex
parameter_t(sound_ray, soundray_comesfrom, normal, vertexes) -> float
                                              # 直線と平面の交点パラメータ t（元コード 384行）
node_renew(sound_ray, soundray_comesfrom, t) -> ndarray (3,)
                                              # 交点 node = 基点 + t * 音線（元コード 388行）
```

#### FaceArrays ― 交差判定のベクトル化（F-1 高速化、2026-08-14）

```python
faces = mm.FaceArrays(mesh)
hit_id, distance, node = faces.nearest_hit(origins, directions)   # (A,3) の束を一度に
```

scalar 版（`collision_distance` など）は「音線 1 本 × 面 1 枚」を 1 回ずつ扱う。
元コードの二重ループをそのまま写したもので読みやすいが、
音線 n 本 × 反射 k 回 × 面 m 枚ぶん Python のループが回るので実用速度に届かない。

`FaceArrays` は面をまとめて配列にしておき、**音線の束 × 全面**を一度の配列演算で処理する。

| 前もって持つもの | 内容 |
|---|---|
| `normal` (M,3) / `d` (M,) | 平面の方程式 n·x + d = 0 |
| `v0` / `v1` | 三角形内部判定の基準点 |
| `edge_from_v0_a/b`, `edge_from_v1_a/b` | 内部判定に使う辺ベクトル |

実装の要点:

- `denominator < 0`（表側から向かっている）と `t > 0` で候補を絞ったあと、
  **`np.nonzero` で (音線, 面) の組を 1 次元に潰してから**内部判定する。
  (A, M, 3) の配列を作らないのでメモリが小さい
- 最寄り面の選択は `np.lexsort((距離, 音線))` + 先頭抽出。
  lexsort は安定なので、距離が同じなら面インデックスの小さいほうが残る。
  scalar 版が `distance_j < min_distance`（狭義）で先勝ちにしているのと**同じ結果**
- (音線 × 面) の要素数が大きくなりすぎないよう音線側を自動で分割する

**scalar 版とビット単位で一致する**ことをテストで確認している（面 ID・距離・受音判定）。
scalar 版は**消していない**。読みやすい参照実装として、また一致確認の基準として残してある。

#### 両面判定 `two_sided`（2026-08-14）

```python
faces = mm.FaceArrays(mesh, two_sided=True)
hit_id, distance, node = faces.nearest_hit(origins, directions, ignore=last_face)
```

既定（`False`）は**法線の側から来た音線しか当たらない**（元コードと同じ）。
法線が逆を向いた面はすり抜けるので、モデルの誤りが結果に出て気づける。

`True` にすると裏からの入射も当てる。反射ベクトル `v - 2(v·n)n` も
エネルギー減衰の `cosθ = |v·n|` も**法線の向きに依らない**ので、
当たり判定だけ両面にすれば法線がまちまちな板モデルがそのまま計算できる。

⚠ **`ignore` が必須**。反射直後の基点は面の上に乗っているので、丸め誤差で `t` が
ごくわずかな正になり、**同じ面にもう一度当たる**（self-intersection）。片面判定では
反射後の音線が裏側を向くので自動的に捨てられており、両面にした瞬間に表に出る。
直前に当たった面を `ignore` に渡して外す。直線は同じ平面と 2 回は交わらないので取りこぼしは無い。
`JR\研修室.dxf` では、これを入れる前は音線 2000 本中 1952 本が室外へ逃げ、入れたら 3 本になった。

⚠ **音線追跡（`loop_reflectionmesh`）とバックトレース（`loop_noredundancy`）で
`two_sided` を必ず揃えること。** 食い違うと、追跡側が通した経路を
バックトレース側が全部却下する。`procedure.process(two_sided=...)` が両方へ渡す。

※ 法線を直せるなら `read_dxffile` の `orient_normals='inward'` のほうが望ましい
（ビューアの裏面（赤）表示も正しくなる）。両者は同じ結果になることを確認済み
（`JR\研修室.dxf` でパルス列 2298 本・バンド別エネルギー合計まで一致）。

### sound_ray.py

音源・虚音源・音線ベクトル操作。すべて純関数。

```python
soundray_generator(ray_number: int) -> ndarray (ray_number, 3)
```
Fibonacci スパイラルで球面上に均等分布する単位音線ベクトル群を生成（元コード 318〜326 行）。
※2026-08-12 に `np.zeros(ray_number, 3)` → `np.zeros((ray_number, 3))`、`sound_rays(i, 2)` → `sound_rays[i, 2]`（2 か所）を修正。
※添字は 0 始まりで、Fortran（1 始まり、最終音線が z>1 で NaN 化）と 1 本ずれるが数学的には Python 側が正しい形。
検証: nray=1000 で NaN なし・全ベクトルのノルム 1.0・重心が原点。

```python
normalized_soundray(sound_ray) -> ndarray (3,)          # L2ノルムで正規化（関数名は normalized の綴りミス）
reflection_generator(sound_ray, normal) -> ndarray (3,) # 鏡面反射ベクトル r = v - 2(v・n)n → 正規化
soundraycomesfrom_renew(node) -> ndarray (3,)           # 基点を交点に置換（恒等関数）
soundray_renew(imaginarysound_point, soundray_comesfrom) -> ndarray (3,)
                                                        # 虚音源 - 基点 で反射音線を再構成（元コード 1101行）
energy_decay(sound_ray, normal, absorption, initial_energy) -> float
```
`energy_decay` は斜入射を考慮した反射エネルギー減衰（元コード 1091〜1094 行）。
反射係数 R = |((1+√(1-α))cosθ − (1−√(1-α))) / ((1+√(1-α))cosθ + (1−√(1-α)))| として R² × 累積エネルギー を返す。

**元の数式の妥当性（2026-08-12 検証済み）**: これは局所反応性（locally reacting）壁面の斜入射圧力反射率で、導出は
法線入射吸音率 α₀ → |R₀| = √(1−α₀) → 規格化音響インピーダンス（実数と仮定）z = (1+|R₀|)/(1−|R₀|)
→ R(θ) = (z·cosθ − 1)/(z·cosθ + 1) の分母分子に (1−|R₀|) を掛けたもの。**式として正しい**。
θ→90° で |R|→1（掠め入射は吸音されない）、θ=0 で残存エネルギー = 1−α₀ となり整合。

※2026-08-12 修正: 式中の `(1−√(1-α))` が 2 か所とも `(1+√(1-α))` になっていた。
　旧実装は分母分子の共通因子 (1+√(1-α)) が約分され R = (cosθ−1)/(cosθ+1) となり、
　**吸音率に一切依存せず・垂直入射で常に R=0（＝どんな材料でも完全吸音）** という致命的な挙動だった。
※`abs(...)**2 * enertmp` の「2 乗の後に自身を掛けている＝3 乗では？」というコード内の疑問は誤読。
　2 乗するのは圧力反射率 R、最後に掛けるのは累積エネルギー `enertmp` で別物。3 乗ではない。
※想定は垂直入射（法線入射）吸音率。残響室法吸音率は 1 を超えることがあり √(1−α) が NaN になるため、
　使う場合は Paris の式を逆解きして α₀ / z を求める前処理が必要。
検証: α₀ = 0.0/0.2/0.5/0.9 で垂直入射の残存エネルギーが厳密に 1−α₀。掠め入射(89.9°)で 0.96。6 バンド一括計算も動作。

### receiver_sphere.py

```python
inside_sphere(sphere_radius, sound_ray, soundray_comesfrom, receiver_point, min_distance) -> bool
```
音線（線分）が受音球を通過したかの判定（元コード 649〜663 行）。
判定 3 条件: ①受音点から音線への垂線距離 ≤ 球半径、②足までの射影距離が最寄り壁より手前、③射影距離 ≥ 0（前方）。
※2026-08-12 修正: ②の比較に垂線距離を使っていたのを射影距離（`inner_product`）に変更（元コード 663 行 `distd .le. disttmp`）。
　旧実装では壁の向こう側にある受音球でも受音と誤判定しうる状態だった。
検証: 壁より手前の受音点 → True、壁の向こうの受音点 → False（旧実装は True）、後方の受音点 → False。

### read_dxffile.py

**2026-08-12 に本実装**（それまではスタブ）。作図ルールは `docs/DXFデータの作り方.md` に別途記載。

```python
read_model(file_name, unit=None, absorption_table=None, default_absorption=None,
           orient_normals="auto", band_number=6,
           source_layers=("src","source","音源"),
           receiver_layers=("rec","receiver","受音点"), verbose=True) -> DxfModel
read(file_name, ...) -> list[Mesh]          # mesh だけ返す簡易版（procedure.py から使う）
read_absorption_csv(file_name, band_number=6) -> dict[str, ndarray]
face_normal(v1, v2, v3) -> ndarray | None   # 外積→正規化。縮退面なら None
signed_volume(triangles) -> float           # 閉形状の法線の向き判定に使う
open_edge_count(triangles, tol=1e-9) -> int # 開いた辺の本数。0 なら閉じている
winding_is_consistent(triangles) -> bool    # 巻き順の一貫性
mesh_shells(triangles) -> list[list[int]]   # 辺の共有で連結成分に分割
analyse_shells(triangles) -> list[dict]     # シェルごとの 閉/開・体積・法線の向き・外殻か
triangulate_polygon(points) -> (tris, info) # 多角形→三角形。凹み対応・ねじれ検出・耳刈り法
quad_warp(points) -> float                  # 四角形のねじれ量（平面からのずれ÷代表辺長）
quad_warp_distance(points) -> float         # 四角形のねじれ量の実寸 [m]
polygon_area(points, normal) -> float       # 多角形の面積（分割の妥当性チェックに使う）
```

`DxfModel` の属性: `mesh` / `source_points` / `receiver_points` / `unit_scale` /
`unit_source` / `layer_counts` / `skipped` / `extents` / `is_closed` / `open_edges` /
`volume` / `shells` / `winding_consistent` / `polygon_notes` / `face_sources`、
および `summary()`。

**パーサは自前**（`ezdxf` などの外部依存なし）。ASCII DXF は「グループコード / 値」が
1 行ずつ交互に並ぶだけなので、2 行ずつ読んで `(code, value)` のタプル列にすれば済む。
エンコーディングは UTF-8 → CP932 → latin-1 の順に試す（日本語レイヤ名対応）。

読むエンティティ:

| エンティティ | 条件 | 扱い |
|---|---|---|
| `POLYLINE` | `70` に 64（ポリフェイスメッシュ） | 後続の `VERTEX` を集め `SEQEND` で終了 |
| └ `VERTEX` | `70 = 192` | 頂点座標（`10`/`20`/`30`） |
| └ `VERTEX` | `70 = 128`（64 が立っていない） | 面レコード。`71`〜`74` が頂点番号（1 始まり、符号は辺の可視性なので `abs()`） |
| `POLYLINE` | `70` に 1（閉じている）でポリフェイスメッシュでない | **閉じたポリライン＝面の輪郭**とみなして三角形分割（`test2.dxf` がこの形式。頂点数の制限なし） |
| `3DFACE` | — | `10`〜`13` 系の 3〜4 頂点。4 点目が 3 点目と同じなら三角形 |
| `LWPOLYLINE` | `70` に 1（閉じている） | 面の輪郭とみなす。頂点は `10`/`20` の繰り返し、高さは `38`。**押し出し方向（`210`〜`230`）が Z 以外は未対応**（読み飛ばして報告） |
| `POINT` | レイヤが `src` 系 | 音源座標 |
| `POINT` | レイヤが `rec` 系 | 受音点座標 |

面がどのエンティティ由来かは `DxfModel.face_sources`（Counter）に記録し `summary()` で報告。
読み飛ばした非対応エンティティも**種類別に Counter で数える**（原因が分かるように）。

**多角形 → 三角形の分割**（`triangulate_polygon()`）。CAD 側で三角形に割る必要はない。

- **三角形**: そのまま通す
- **四角形**: `_diagonal_is_inside()` で**内側を通る対角線を選ぶ**
  （他の 2 頂点が対角線の両側に分かれるか。素朴な扇状分割は凹んだ四角形で多角形の外に
  三角形を作ってしまう）。片方が対角線上（符号 0）の退化ケースも内側扱いにする
  （面積ゼロの三角形は後段で捨てられる）
- **ねじれた四角形**: `quad_warp()`（相対）と `quad_warp_distance()`（実寸 [m]）で
  平面からのずれを測り、`WARP_TOLERANCE` を超えたら警告。
  どの対角線で切るかで形が変わるため CAD 側で直してもらうが、**実寸を mm で示して
  影響の目安（1mm 未満なら無視できる／20mm 超は作図ミスの可能性）まで報告する**。
  作図方法は `docs/DXFデータの作り方.md` 1b 節（UCS + PLINE / 押し出し / 三角形）
- **5 角形以上**: `_ear_clip()` で耳刈り法。「耳」＝凸な頂点で、その両隣を結んだ三角形に
  他の頂点が入らないもの。これを 1 つずつ切り落とすので**凹んだ多角形でも正しく分割できる**。
  閉じたポリラインで面を描くと 5 角形以上が普通に来る（`test2.dxf` は 9 角形）
- **分割の妥当性チェック**: `polygon_area()`（靴ひも公式に相当）と三角形の合計面積を
  `_area_mismatch()` で照合する。自己交差した多角形などで破綻すると誤差が出るので検出できる
- 縮退面（`|n| < 1e-12`）は `face_normal()` が None を返すので捨てて件数を記録

結果は `DxfModel.polygon_notes`（ねじれた四角形の枚数・最大ねじれ量・対角線変更の枚数・
耳刈り法の枚数・最大面積誤差・分割失敗の枚数）に入る。

**単位換算**: ヘッダ `$INSUNITS` のコード（4=mm, 6=m など、`INSUNITS_TO_METER` に 18 種）から
m へ換算する。`unit='mm'` / `unit=0.001` のような明示指定が優先される
（`$INSUNITS=0` の unitless 対策）。

**法線の向き**: 面は片側だけ反射する（`v·n < 0` のときのみ衝突と判定）。したがって
**法線は「音が通る空気側」を向く**のが作図ルール（室の外殻は内向き、室内の物体は外向き、
両面反射させたいものは厚みを持たせる）。**この正解を持っているのは CAD モデル自身**なので、
既定は `'cad'`（そのまま使う）。
また面の片側性により**閉じた室でなくても計算できる**（一面反射の検証など）。

| 値 | 挙動 |
|---|---|
| `'cad'`（既定） | CAD の巻き順をそのまま使う |
| `'flip'` | 全反転。**元コード 276行 `ynnmrev='y'` に相当** |
| `'shells'` | シェル単位で空気側へ揃える（外殻は内向き、内側は外向き）。巻き順が一貫していなければ補正を中止 |
| `'inward'` | **面ごとにレイの偶奇で室内側へ揃える**（2026-08-14 追加。下記） |

⚠ かつてあった `'auto'` / `'toward'` / `'away'` / `'outward'` は**廃止**（`ValueError`）。
法線を音源方向や重心方向へ面ごとに向ける方式は、**凸凹の壁や宙に浮いた家具で破綻する**ため。

**`'inward'`（`orient_inward()`）** は 2026-08-14 に、実務の CAD データ（`JR\研修室.dxf`）で
**床は上向き・天井も上向き・壁はバラバラ**という状態に遭遇して追加した。
CAD で面を 1 枚ずつ描くと法線は巻き順と押し出し方向で決まるので、意識しないとこうなる。
`'shells'` は面のつながり（巻き順の一貫性）を要求するため、こういう「板の寄せ集め」には効かない。

判定は**レイの偶奇**による。閉じた面に囲まれた領域の内側から外へレイを飛ばすと境界を
必ず**奇数回**横切る。そこで面の重心を法線側へわずかに浮かせた点から飛ばし、
交差回数が奇数なら「法線側が室内」＝そのまま、偶数なら反転する。
辺や頂点をかすめたときの数え間違いを避けるため、方向を 9 本ばらまいて**多数決**を取る
（固定シードなので同じモデルなら毎回同じ結果）。交差判定は Möller–Trumbore で、
**平面の法線を使わないので CAD の巻き順に左右されない**。
判定が割れた面があれば `summary()` の「法線の向き」に本数が出る。

**OCS（オブジェクト座標系）**: LWPOLYLINE のような平面図形は、頂点を**その図形が乗る平面上の
2 次元座標**で持っている。押し出し方向（グループコード 210/220/230）が Z なら OCS = WCS だが、
**鉛直な壁は押し出し方向が水平**になるため変換が要る。`ocs_axes()` / `ocs_to_wcs()` が
DXF 仕様の **Arbitrary Axis Algorithm** で変換する（軸の取り方は仕様で決まっていて任意性はない）。
2026-08-14 まで未対応で、該当する面をまるごと読み飛ばしていた（`JR\研修室.dxf` で壁 19 枚が消えていた）。

**診断関数**（GUI から法線を調整するための土台。TODO G-8）:

| 関数 | 内容 |
|---|---|
| `open_edge_count()` | 開いた辺（三角形 1 枚にしか属さない辺）の本数。0 なら閉じている |
| `winding_is_consistent()` | 同じ有向辺が 2 回現れないかで巻き順の一貫性を判定 |
| `mesh_shells()` | 辺の共有で連結成分（シェル）に分割 |
| `analyse_shells()` | シェルごとに 閉/開・体積・法線の向き（inward/outward）・外殻かどうかを返す |
| `signed_volume()` | 閉じたシェルの法線の向き判定に使う（V>0 で外向き） |

結果は `DxfModel.shells` / `.winding_consistent` / `.is_closed` / `.open_edges` / `.volume`
に入り、`summary()` がシェルごとに「OK」「★要確認（outward が空気側）」まで表示する。

**吸音率**: レイヤ名を材料の識別子として `absorption_table`（dict または CSV パス）と突き合わせる。
引く順序は **① レイヤ名と完全一致 → ② レイヤ名の先頭の数字を材料 ID とみなす**
（`layer_number()`。`01__研修室_床` → ID `1`。先頭のゼロは落とす）。
②は実務の作図で「**レイヤ名の先頭 2 桁 ＝ 吸音率表の ID**」という運用があったため 2026-08-14 に追加。
どのレイヤがどの材料になったかは `DxfModel.layer_materials` に入り、`summary()` の
`レイヤ→材料` 行に出る（`AbsorptionTable.names` を持たせてあるので**ID で引いても材料名が言える**）。
未登録のレイヤは既定値（0.1）を使い、レイヤ名を列挙して警告する。
`read_absorption_csv()` は 2 形式を自動判別する。

| 形式 | 列 | 備考 |
|---|---|---|
| (A) 元コード付属 `absorption.csv` | `ID, 材料名, a1〜a6` | **ID と材料名の両方をキーに登録**するのでレイヤ名がどちらでも引ける。CP932 なので UTF-8 → CP932 → latin-1 の順にデコードを試す |
| (B) `data/absorption_sample.csv` | `材料名, a1〜a6` | 2 列目が数値かどうかで (A) と判別 |

`if __name__ == "__main__":` でリポジトリ直下の `test.dxf` を読むテストが走る。

### view_model.py

読み込んだモデルを 3D 表示するビューア。移植元のない**新規モジュール**。

```python
view(dxf_path, out_path=None, absorption=None, unit=None,
     orient_normals="cad", open_browser=True) -> (out_path, DxfModel)
export_html(model, out_path, title="", subtitle="") -> out_path
build_payload(model, normal_ratio=0.06) -> dict     # HTML に埋め込む JSON
```

コマンドラインからも使える。

```
cd geosim
python view_model.py ..\test2.dxf --absorption ..\absorption.csv
```

**自己完結した HTML（WebGL）を書き出してブラウザで開く方式**。作った当初この環境に
matplotlib / pyvista が入っておらず（`tkinter` のみ）、`pip install` を前提にしたく
なかったため。**2026-08-14 に依存ライブラリを入れた後も残してある**。理由は、
①相手に環境構築を求めずにモデルを共有できる、②将来 GUI を Web ベースにするなら
そのまま土台になる、から。ネイティブウィンドウ版は `view_model_gui.py`（後述）。

表示するもの:

| 要素 | 表現 |
|---|---|
| 三角形要素 | 面 + **辺を描く**ので分割が見える |
| 法線 | 重心から伸びる矢印（矢じり付き。長さはモデル対角の 6%） |
| **法線の向き** | **裏側（法線と反対から見ている側）を赤で塗る**。音線が通り抜ける側が一目で分かる |
| レイヤ | 色分け + 表示/非表示の切り替え |
| 音源・受音点 | 赤 / 青の点 |
| 読み込み結果 | `DxfModel.summary()` をパネルに表示 |

操作: ドラッグで回転、ホイールで拡大縮小、右ドラッグ（または Shift+ドラッグ）で平行移動、
視点プリセット（等角 / 上 / 正面 / 側面）。

実装上の要点:

- 表裏の判定は `gl_FrontFacing`（頂点の巻き順）ではなく**保持している法線と視線の内積**で行う。
  読み込み側は巻き順を変えずに法線だけ反転することがあるため
- **不透明度をレイヤごとに変えられる**（2026-08-14 追加、`add_opacity_control`）。
  左の縦スライダで設定、`Tab` で対象を切替（すべて ↔ 各レイヤ）、`m` でモデル表示 ON/OFF。
  起動時に決めるなら `--opacity`（全体）と `--layer-opacity "1=0.6,2=0.05"`（レイヤ別）。
  「壁だけ薄くして中の様子を見る」「床は残す」といった見方をするため。
  ※ 左上のチェックボックスは**表示のオンオフ**で、不透明度とは別物。
  ※ pyvista の `add_slider_widget` は**生成時にコールバックを 1 回呼ぶ**ので、
    そのまま作るとレイヤ別に指定した値が全部スライダの初期値で上書きされる。
    `ready` フラグでその 1 回を無視している。
- レイヤの表示切り替えは、**頂点をレイヤ順に並べて描画範囲を飛ばす**方式。
  シェーダで uniform 配列を動的インデックスする実装は GLSL ES の移植性に不安があるため避けた
- 出力 HTML（`*_view.html`）は `view_model.py` で再生成できるので `.gitignore` で除外

### view_model_gui.py

同じモデルを **PyVista（VTK）のネイティブウィンドウ**で表示するビューア。
移植元のない**新規モジュール**（2026-08-14）。

```python
view(dxf_path, absorption=None, unit=None, orient_normals="cad",
     screenshot=None, show_normals=True) -> DxfModel
build_plotter(model, title=..., off_screen=False, show_normals=True,
              normal_ratio=0.06, window_size=(1280, 860)) -> pv.Plotter
triangles_to_polydata(triangles) -> pv.PolyData
normal_arrows(poly, length) -> pv.PolyData
japanese_font() -> str | None
```

```
cd geosim
python view_model_gui.py ..\test.dxf --absorption ..\absorption.csv
python view_model_gui.py ..\test.dxf --screenshot shot.png   # ウィンドウを開かず画像だけ
```

表示内容は HTML 版と同じ（三角形要素・法線矢印・裏面を赤・レイヤ色分けと表示切り替え・
音源/受音点・`summary()`）。色は `view_model.LAYER_PALETTE` を import して共有しているので
2 つのビューアで見え方が揃う。

操作: ドラッグで回転、ホイールで拡大縮小、中ドラッグで平行移動、
`z`/`x`/`c`/`v` で上/正面/横/等角、`n` で法線矢印、`w`/`s` でワイヤフレーム/面（VTK 既定）、
`r` で視点リセット、`q` で終了。

実装上の要点:

- **巻き順を `t.normal` に合わせ直してから PolyData を作る**（`triangles_to_polydata`）。
  `read_dxffile` の `orient_normals='flip'` / `'shells'` は**法線だけを反転**して頂点順序を
  触らないため、そのまま渡すと VTK の表裏判定がモデルの法線と食い違う。揃えておけば
  VTK の `backface_params` に裏面を赤で塗らせるだけで「法線の向きの確認」が成立する
  （HTML 版が法線と視線の内積で自前判定しているのと目的は同じ、手段が違う）
- VTK の既定フォントは日本語を持たないので、`japanese_font()` で Windows 標準の
  日本語フォント（meiryo.ttc など）を探して `add_text(font_file=...)` に渡す。
  これをしないとレイヤ名（＝吸音材名）が豆腐になる
- チェックボックスウィジェットは interactor が要るので **`off_screen=True` では追加しない**
- `--screenshot` は off-screen レンダリング。**表示を自動で検証できる**ようにするためのもの
  （実際 `test.dxf` は法線内向きなので外から見ると全面が赤くなる、で確認した）

**どちらのビューアを使うか**

| | view_model.py（HTML） | view_model_gui.py（PyVista） |
|---|---|---|
| 依存 | なし | pyvista + vtk |
| 共有 | HTML を渡すだけ | 相手にも環境構築が要る |
| Python から操作 | しにくい | しやすい（Plotter を触れる） |
| 発展先 | Web ベース GUI | 音線・音粒子の重ね描き（G-1/G-2）、デスクトップ GUI（G-7） |

### view_rays.py

音線と音粒子の可視化（出力① G-1 / 出力② G-2）。移植元のない**新規モジュール**（2026-08-14）。
`ray_recorder.RayRecorder` が保存した npz を読み、`view_model_gui` のモデル表示に重ねる。

```python
RayLog(npz_path)                       # 軌跡データの読み込みと前処理
    .selection(received_only, max_rays)         # 描く音線を選ぶ
    .line_polydata(index, colour, max_reflection)   # ① 折れ線
    .positions_at(time, index)                  # ② 任意時刻の粒子位置（一括）
    .energy_at(time, index, band)
add_rays(plotter, raylog, ...)         # ① 音線を描く
ParticleAnimation / animate(...)       # ② アニメーション
save_movie(raylog, model, "out.gif")   # ② GIF に書き出す
view(dxf, raylog, mode="rays"|"particles", ...)
```

```
cd geosim
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --received-only --max-reflection 3
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --mode particles
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --mode particles --movie 広がり.gif
```

**① 音線（`--mode rays`）**

反射経路を折れ線で描く。点ごとにスカラーを持たせるので、1 本の線の中で色が変わる。
`--color` で `energy`（バンド平均の dB）/ `time` / `reflection` / `ray` を選ぶ。
受音した経路は太い黄色で重ね描きする（`--received-only` のときは全部が受音経路なので出さない）。

> **`--max-reflection` は「その回数までで折れ線を打ち切る」動作**。
> 「その回数以下の音線だけ残す」ではない。閉じた室では全音線が上限まで反射するので、
> 絞り込みでは何も残らない。初期反射だけを見たいときに使う。

> **`--max-rays` の既定は 80**。多すぎると線が重なって何も読めなくなる。
> 間引きは**等間隔**（ランダムに抜くと見た目の密度が変わるため）。

**② 音粒子（`--mode particles`）**

離散化時間ごとに粒子が飛ぶ様子を見る。色はエネルギー [dB]。
`スペース` 再生/停止、`←` `→` コマ送り、`Home` 先頭、下のスライダで時刻指定。

実装上の要点:

- **再生は VTK のタイマーではなく自前のループで回す**（`run_animation`）。
  最初は `Plotter.add_timer_event()` を使ったが**発火しなかった**
  （20 回 × 50 ms なら 1 秒で終わるはずが 60 秒経っても 0 回）。
  VTK の `CreateRepeatingTimer` は interactor の初期化後でないと効かないが、
  pyvista の API は `show()` より前に呼ぶ形になっているためと思われる。
  症状は「コマ送りは効くのに再生されない」。
  `show(interactive_update=True)` で非ブロッキング表示にして
  こちらから `update()` を呼ぶ方式に変えたら確実に動いた（45 コマ 1.5 秒）
- **任意時刻の粒子位置を全音線ぶんまとめて計算する**（`positions_at`）。
  可変長の軌跡を「行ごとに揃えた 2 次元配列（末尾を +inf で埋める）」に直しておくと、
  `(距離 <= d) の個数 - 1` で区間の添字が求まり、音線ごとの `searchsorted` が要らない。
  `RayTrajectory.position_at()` を毎フレーム本数分呼ぶより速く、結果は一致する（テスト済み）
- **点群は 1 つの actor にまとめる**。音線ごとに actor を作ると本数分の描画呼び出しになる
- 粒子数はフレームごとに変わるので、`points` だけでなく**頂点セル（`verts`）も作り直す**。
  座標だけ差し替えると古いセル（最初の 1 点）しか描かれない
- 動画は追加の依存を増やさず、コマを画像として集めて **Pillow で GIF** にする
  （Pillow は matplotlib の依存で既に入っている）
- カラーバーとスライダの見出しは VTK の既定フォントで描かれるため**日本語が出せない**。
  ここだけ英字にしてある（本文側は日本語フォントを指定している）

### loop_reflectionmesh.py

```python
loop(soundsource_point, receiver_point, soundray_list, nref, mesh, sphere_radius,
     recorder=None) -> list[list[int]]   # 受音した経路ごとの反射面ID履歴（先頭要素は -1）
```
音線追跡の本体（元コード 524〜717 行）。3 重ループ構造:

```
for 音線 i:                       # 元コード 527行 do i = 1, nray
    履歴初期化（先頭に番兵 -1）
    for 反射回数 k in range(nref+1):   # 元コード 545行 do k = 0, nref（k=0 が直接音の区間）
        音線を正規化
        collision / min_distance を初期化   ← 面ループの「外」
        for 面 j:                 # 元コード 576行 do j = 1, sfcount
            collision_distance() で衝突判定 → 最寄り面 mesh_nearestid を更新
        受音球判定 inside_sphere()
        if 受音: 履歴のコピーを 2 次元リストへ保存      # ★壁ID追記より前
        if 衝突なし: break
        履歴に壁 ID を追記                             # ★無条件
        交点を算出 → 基点を移動 → 反射ベクトルを生成
```

**2026-08-12 に traceff の扱いを確定し、全面的に書き直した**（詳細は `docs/技術説明書.md` 5.8 節）。

修正点:

| # | 修正前 | 修正後 |
|---|---|---|
| A-1 | 2 次元リストへの append が反射ループ毎周・同一リスト参照 | 受音時のみ `list()` でコピーして append |
| A-2 | `if inside:` の中でのみ壁 ID を追記 | `if collision:` で**無条件に**追記（受音とは独立） |
| A-3 | 面ループ内で `soundray_comesfrom` を上書き | 面ループは判定のみ。交点算出は最寄り面確定後に 1 回 |
| A-4 | `collision` を面ループ**内**で初期化（最後の面の結果しか残らない） | 面ループの**外**で初期化 |
| A-5 | 不要な `from numpy.ma.core import count` | 削除 |
| A-7 | `for k in range(nref)` | `range(nref + 1)`（元コード `do k = 0, nref` は nref+1 周） |
| A-8 | — | 可視化用の `recorder` フックを追加 |

⚠ 戻り値の各履歴の**先頭要素は番兵 `-1`** で面 ID ではない（元コード `tractmp(0) = -1`）。
反射回数 = `len(history) - 1`。下流で `mesh[history[k]]` とすると `-1` が最後の面を指すので注意
（TODO A-9 で扱いを再検討予定）。

### ray_recorder.py

可視化用の音線軌跡レコーダ。移植元のない**新規モジュール**（`docs/出力・可視化方針.md` 参照）。

```python
class RayRecorder(total_rays, max_rays=2000, sound_velocity=343.0,
                  band_number=6, record_energy=True)
class RayTrajectory   # 音線 1 本分の軌跡
```

`loop()` に `recorder=` で渡すと、間引いた音線について
**節点座標・累積距離・反射面 ID・バンド別エネルギー・受音イベント・終了理由**を記録する。
本線（受音経路の反射面 ID 列）とは別チャンネルで、`recorder=None` なら一切動かない。

- 受音しなかった音線も記録する（音の広がりを見るのが目的のため）
- 間引きは「絶対本数の上限 `max_rays` ＋ 等間隔ストライド」が既定
- `traj.times(c)` で到達時刻、`traj.position_at(t, c)` で任意時刻の音粒子位置（線形補間）
- `recorder.save_npz(path)` で保存（可変長なので連結＋オフセット方式）

### loop_deleteredundancy.py

```python
delete(reflectionhistory_redundancy: list[list[int]]) -> list[list[int]]
```
反射履歴の重複経路削除（元コード 721〜841 行）。

**2026-08-12 に元コードの意図を確定**（それまで「何をしているか分からない」とコメントされていた）:

- `traceffn(i)`（721〜741 行）は、固定長配列 `traceff` の各行を先頭から走査して最初の 0 の位置を探し、
  **その経路の有効反射回数（＝行の実質的な長さ）**を求めている。
  Python のリストは可変長なので `len()` がそのまま相当し、この処理自体が不要。
- その後の「反射回数でソート → 切替点探索 → 同一反射回数内で総当たり比較」は、
  **「反射回数が違えば絶対に重複しない」ことを使って比較対象を絞る高速化**。
  Python では tuple 化してハッシュで一括除去すれば同じ結果が O(n) で得られる。

実装は `list → tuple → dict.fromkeys → list`。**`set` ではなく `dict.fromkeys` を使うのは
入力順を保って結果を決定的にするため**（`set` だと実行のたびに順序が変わり、
結果の突き合わせやデバッグがしづらい）。毎回数百行を吐いていたデバッグ print も削除した。

### loop_noredundancy.py

虚音源法バックトレース（元コード 876〜1134 行）。**2026-08-14 に全面書き直し**。

```python
loop(soundsource_point, receiver_point, reflectionmeshid_history, mesh,
     sound_velocity=340.0, band_number=None, filename=None, verbose=True) -> PulseList
backtrace_path(soundsource_point, receiver_point, wall_ids, mesh,
               band_number, sound_velocity) -> dict | None
image_sources(soundsource_point, wall_ids, mesh) -> ndarray (n+1, 3)
```

処理:

```
for 非重複経路 i:
    反射履歴に沿って虚音源列 isrc(k) を鏡像で順次生成      # 元コード 900〜926行
    受音点から最終虚音源へ向けて逆向きに追跡:              # 元コード 948行 do k = ktmp, 0, -1
        各段で最寄り面を探索し、履歴の面と一致するか検証   # 違えば却下
        一致すれば energy_decay() でバンド別エネルギーを減衰、基点と音線を更新
        k=0 で遮蔽チェック（音源より手前に壁があれば却下）
    通れば「反射回数, 到来時刻, 到来方向, バンド別エネルギー」を記録
```

**`PulseList`**（元コード 1080 行の 11 列出力に対応）

| 属性 | 形状 | 元コード | 内容 |
|---|---|---|---|
| `reflection_count` | (n,) | ktmp | 反射回数 |
| `time` | (n,) | rtime | 到来時刻 [s] |
| `distance` | (n,) | — | 経路長 [m] |
| `direction` | (n,3) | -vtgt | 到来方向の**単位ベクトル**（受音点→虚音源） |
| `energy` | (n,6) | enertmp | バンド別エネルギー |

元コードの `-vtgt` は正規化されておらず長さが経路長そのものなので、「方向」と「距離」が
1 列に混ざっていた。分けておくと**出力⑤（伝搬方向の可視化）にそのまま使える**。
`direction * distance` が元コードの `-vtgt` に一致する。

**元コードとの相違点：鏡像の式に abs() を使わない**

元コード 911 行は面までの距離に `abs()` を掛けている。

```
temp = |n・p + d| / |n|          →  isrc(k) = isrc(k-1) - 2 * temp * n
```

鏡像の正しい式は符号付きで `p - 2 (n・p + d)/|n|^2 * n`。両者が一致するのは
p が面の**表側**（n・p + d > 0）にあるときだけで、裏側の虚音源に abs() を使うと
面を跨がずに**遠ざかる方向へ動いてしまう**。凸な部屋なら実害は出にくいが、
凹んだ形状や法線を反転させたモデルでは崩れるので、符号付きで実装した。

**検証**（2026-08-14、test.dxf = 2×3×1 m 直方体）

| 項目 | 結果 |
|---|---|
| 直接音の距離 | 解析解と**誤差 0.000e+00** |
| 1 次反射 6 本（6 面ぶんちょうど） | 全て解析解と**誤差 0.000e+00** |
| 2 次反射 18 本 | 全て虚音源の解析距離と一致 |
| 到来方向 | 受音点→音源の単位ベクトルと一致 |
| エネルギー | 入射角と吸音率から手計算した \|R(θ)\|² と**誤差 0.000e+00** |

### impulse.py

パルス列 → インパルス応答の合成（元コード make_ipls_freq_monaural_fortran.f90 全体）。
**2026-08-14 に全面書き直し**。

```python
impulse_response_from_pulses(filename, pulses, sound_velocity=340.0, sampling_frequency=44100.0,
                 max_time=1.0, compensate_filter_delay=False, verbose=True) -> (t, ir)
impulse_response(time, energy6, ...) -> (t, ir)      # ファイル出力なし
```

| 関数 | 対応元コード | 処理 |
|---|---|---|
| `third_octave_bands()` | 46 行 | 1/3 オクターブ 32 バンドの中心周波数（15.625 Hz〜20 kHz） |
| `airdamping_coefficient(mf)` | 47 行 | 空気吸収係数 `mair = 1.81e-8 * mf^1.57`（20℃・湿度 40%） |
| `expand_6_to_32(energy6)` | 94〜125 行 | **6 オクターブバンド → 32 バンドの展開表** |
| `apply_air_absorption(...)` | 127 行 | `E * exp(-mair * c0 * t)` |
| `transfer_function(...)` | 136〜143 行 | `Σ √E · e^(−j2πf·t) / (t·c0)` の重ね合わせ |
| `bandpass_edges(mf, fmax)` | 152〜154 行 | 各バンドの遮断周波数 `mf * 2^(±1/6)` |
| `filter_bandpass(numtaps, wmin, wmax)` | fir1_bandpass | ハミング窓付き FIR バンドパス |
| `extend_to_negative(hfp, nn)` | 173〜180 行 | 伝達関数を負の周波数へエルミート拡張 |
| `write_impulse_response(...)` | 232〜235 行 | 時間ベクトルを付けて CSV 出力 |

**周波数バンドは可変。既定 8 バンド（63〜8k Hz）**

63 Hz と 8 kHz を対象外にする運用もあるので 6 バンド（125〜4k Hz）でも動く。
オクターブ → 1/3 オクターブの割り当ては `band_mapping()` が
**対数軸で中心周波数が最も近いもの**を選ぶ規則で行う。
元コードは 6 バンド用の対応表を手書きしていた（94〜125 行）が、この規則は
**その表を完全に再現する**（テストで確認済み）。

**直したもの（E-3 / E-4 / E-7）**

| 箇所 | 旧実装 | 正 |
|---|---|---|
| ハミング窓 | `0.54 * 0.46 * cos(...)` | `0.54 + 0.46 * cos(...)` |
| 伝達関数の重ね合わせ | `transfer = transfer * ...`（積） | `+=`（和）。積だと初期値 0 に掛かり全部 0 になる |
| 空気吸収係数 | `(15.625*2)^(i/3)` かつスカラー代入 | `15.625 * 2^(i/3)` の配列 |
| 空気吸収の時刻 | `reflection_timing[i]`（i はバンド添字） | 全パルスの時刻ベクトル |
| バンド端 | `(mf*2)^(-1/6)`、指数に fmax が混入 | `mf * 2^(∓1/6) / fmax` |
| 負周波数の添字 | `transfer[j, nn-i+1]` | 0 始まりでは `hfp[nn - i]` |
| フィルタのタップ添字 | `n_shift = i - 1.0 - m`（1 始まり前提） | `k - m` |
| 逆 FFT | 呼ばずに合算していた | `np.fft.ifft` してから実部を合算 |
| `df` | `transfer_function` に渡っておらず未定義参照 | 引数で受け取る |

**信号処理は scipy に任せる（2026-08-14 の方針）**

元コードは Fortran で使えるライブラリが限られていたため FFT もフィルタも自作していた。
Python には scipy があるので、そこは既存ライブラリに任せる。

| 処理 | 元コード | 本実装 |
|---|---|---|
| FFT | 自作 Cooley-Tukey | `scipy.fft` の `rfft` / `irfft`（実数信号用） |
| バンドパス FIR | 自作 `fir1_bandpass` | `scipy.signal.firwin(..., scale=False)` |
| 負の周波数への拡張 | 手動でエルミート対称に詰める | **不要**（`irfft` が前提として持っている） |
| 畳み込み | 巡回畳み込み（長さ nn） | ゼロ詰めして**線形畳み込み** |
| フィルタ長 | nn = 131072（遅れ 1.49 秒） | 既定 8192（遅れ 93 ms、出力からは除去済み） |

`scale=False` が要点。scipy は既定で通過域ゲインを 1 に正規化するが、それをすると
**バンドを足し合わせたときに平坦にならない**。元コードの `fir1_bandpass` は正規化していない。
両者が**相対誤差 1e-16 で一致する**ことをテストで確認している
（元コードの移植版は `_fir1_bandpass_fortran()` に比較用として残してある）。

**結果として、元コードにあった 1.49 秒の遅れは無くなった。**
出力されるインパルス応答は 0 秒が音の始まりで、長さは `max_time` そのもの。

**検証**（2026-08-14）

| 項目 | 結果 |
|---|---|
| `filter_bandpass` vs `scipy.signal.firwin` | 相対誤差 1e-16 |
| 32 バンドの総和（40 Hz〜15 kHz） | 振幅 1.0000（過不足なく全帯域を覆う） |
| 6→32 展開表 | 各オクターブ中心で期待どおり |
| 単一パルスのピーク位置 | 期待位置と**ずれ 0 サンプル** |
| 遅れ補正後のピーク時刻 | 10.0000 ms（入力どおり） |
| 距離減衰 | `1/(t·c0)` と誤差 1e-12 未満 |

### reverberation.py

インパルス応答 → 減衰曲線・残響時間（元コード ipls2rt_fortran.f90）。
移植元はあるが Python 側は**新規モジュール**（2026-08-14）。

```python
reverberation_time(time, ir, rt_filename=None, decay_filename=None,
                   frequencies=None, db_max=-5.0, db_min=-35.0, numtaps=4096) -> dict
decay_curves(...) -> dict          # 保存なし
bandpass(signal, numtaps, lower, upper) -> ndarray
schroeder_integral(x) -> ndarray
```

**`decay_measures()` が EDT / T20 / T30 を一度に返す**（2026-08-14）。
60 dB 減を厳密に見ることは実務でほとんど無く、減衰の直線部分を測って外挿するため。

| 指標 | 評価区間 | 性質 |
|---|---|---|
| **EDT** | 0 〜 -10 dB | 初期減衰時間。**聴感上の響きの短さに近い**。初期反射の密度を反映する |
| **T20** | -5 〜 -25 dB | 暗騒音が高い実測でも取りやすい |
| **T30** | -5 〜 -35 dB | 最も一般的（ISO 3382 の標準） |

減衰曲線は 1 回だけ作って評価区間を変えて読み取るので、`decay_curves` を 3 回呼ぶより速い。
戻り値の dict は `frequencies` / `time` / `decay` (nf, n) [dB] /
`measures` {名前: (nf,)} / `curvature` (nf,) [%]。
`measures` 引数に `{名前: (開始dB, 終了dB)}` を渡せば評価区間を変えられる。

**元コードとの相違点：巡回畳み込みの自作 FIR → Butterworth（scipy）** ★重要

元コードは FFT の掛け算だけで済ませているので**巡回畳み込み**になる。
バンドパスは左右対称な直線位相 FIR なので入力の立ち上がりより**前の時刻**から応答が始まり
（プリリンギング）、巡回畳み込みではそれが**バッファ末尾に回り込む**。
しかもタップ数 = nn なので回り込みはバッファの半分に及ぶ。
結果、減衰曲線の後半がまるごとプリリンギングで埋まり、Schroeder 積分に -20 dB 程度の
**床**ができて -35 dB まで落ちなくなる。

減衰率が既知の合成インパルス応答（指数減衰させた白色雑音）で確かめた結果:

| 方式 | T60 = 0.3 s | 0.5 s | 1.0 s |
|---|---|---|---|
| 巡回畳み込み（元コード, タップ数 = nn） | 算出不可 or 5.89 s | 算出不可 or 5.85 s | 5.77 s |
| 線形畳み込み FIR 4096（`method='fir'`） | 0.29〜0.34 s | 0.48〜0.52 s | 0.98〜1.02 s |
| **Butterworth 6 次（既定）** | **0.29〜0.33 s** | **0.47〜0.51 s** | **0.98〜1.02 s** |

既定は **`scipy.signal.butter` + `sosfilt`**。**IEC 61260**（オクターブフィルタの規格。
ISO 3382 の残響測定が前提にしている）は Butterworth 系なので実務に沿ううえ、
因果 IIR なので回り込みが構造的に起きず、計算量も桁違いに少ない。
FIR で処理したい場合は `method='fir'` を指定できる（精度は同等）。

誤差は**中央値 1.0〜1.9 %**。125 Hz で減衰が短いときだけ 10〜14 % 程度になるが、
これは帯域幅×減衰時間（BT 積）が小さいことによる手法の限界で ISO 3382 でも既知。

**統計残響式（2026-08-14 追加）**

```python
statistical_reverberation(mesh, volume, frequencies=None, atmosphere=None,
                          convert_to_random=True, include_air_absorption=True) -> dict
statistical_reverberation_from_model(model, **kwargs) -> dict | None
surface_summary(mesh, convert_to_random=True) -> dict
triangle_areas(mesh) -> ndarray
```

音線を飛ばさず、**室容積と各面の面積・吸音率だけ**から残響時間を出す。
共通の形は `T = 24 ln(10) V / (c A)`（係数は c=343 のとき 0.161）で、
等価吸音面積 A の作り方が 3 通りある。

| 式 | A | 向いている場面 |
|---|---|---|
| Sabine | `S·ᾱ` | 吸音率が小さいとき（〜0.2） |
| Eyring | `-S·ln(1-ᾱ)` | 吸音率が大きいとき |
| Eyring-Knudsen | `-S·ln(1-ᾱ) + 4mV` | Eyring に**空気吸収**を足したもの |

- ★**名前**：アイリングの式は `-S ln(1-ᾱ)` までで、**空気吸収 `4mV` を足した形は
  ヌードセンの寄与**（ユーザー指摘 2026-08-17。それまで空気吸収込みなのに
  `Eyring` と表示していた）。**2 つを別の列にしてある**ので、
  差がそのまま空気吸収の効きになる。表示名は `STATISTICAL_LABELS`
- ★**Sabine と Eyring には空気吸収を入れない。** 素の古典的な形にしておかないと、
  Eyring → Eyring-Knudsen の差が「空気吸収ぶん」にならない
- ★**Millington は落とした**（2026-08-17 ユーザー判断）。ミリントン・セッテの式は
  面ごとに対数を取るので **α→1 の面が 1 枚でもあると T→0** になり（開口が 1 つで
  残響ゼロ）、実務では使われない。研修室では中高域で Eyring-Knudsen の約半分に出ていた
- **閉じた室が前提**。`statistical_reverberation_from_model()` は
  `model.is_closed` を見て、開いていれば None を返して警告する
- **`Mesh` が持つ垂直入射吸音率を Paris の式で乱入射に変換してから使う**
  （`convert_to_random=True` が既定）。変換しないと残響時間が長く出る
- 空気吸収の項 `4mV` も `atmosphere` から入れる
- `procedure.py` は計算結果の T30 と並べて表示する

test.dxf での突き合わせ（音線 20000 本 / 反射 120 回）:

| 周波数 | T30（計算） | Sabine | Eyring | T30/Eyring |
|---|---|---|---|---|
| 63 Hz | 0.088 | 0.205 | 0.182 | 0.48 |
| 125 Hz | 0.150 | 0.191 | 0.168 | 0.89 |
| 500 Hz | 0.122 | 0.128 | 0.105 | 1.17 |
| 2000 Hz | 0.104 | 0.106 | 0.083 | 1.26 |
| 8000 Hz | 0.081 | 0.100 | 0.079 | 1.03 |

125 Hz 以上は 0.89〜1.26 で概ね合う。63 Hz が外れるのは帯域幅×減衰時間が
小さすぎる（BT ≒ 4）ためで、手法の限界。

**曲率（curvature）**

評価区間を半分にした推定値との食い違いを % で返す（ISO 3382 の C）。
0 に近ければ減衰が直線。**10 % を超えたら結果を疑う**。
`nref` が足りずに後部残響が途中で切れているときの検知に使える。
`procedure.py` 側でも「nref に達した経路がある」場合に警告を出している。

---

## Fortran 側との対応・未移植の機能

| 元コードの機能 | Python 側の状況 |
|---|---|
| 吸音率 CSV（absorption.csv）読み込み | 実装済（`read_absorption_csv()`。サンプル `data/absorption_sample.csv`） |
| DXF 形状読み込み・法線計算 | 実装済（`read_dxffile.py`。単位換算・法線向き自動判定・音源受音点の読み込みは元コードにない追加機能） |
| 有効経路カウント用の事前ループ（backtrace.f90 330〜515 行） | 不要（Python はリストの動的伸長で代替） |
| 音線追跡ループ | **実装済**（`loop_reflectionmesh.py`。音源も受音点も複数取れる） |
| クイックソート + 重複削除 | set 方式で代替済み（`loop_deleteredundancy.py`） |
| バックトレース + 受音リスト出力 | **実装済**（`loop_noredundancy.py`。解析解と一致。経路の束で回す版が本番、1 本ずつの版は参照実装） |
| インパルス応答合成（make_ipls_...f90） | **実装済**（`impulse.py`。scipy 化。時間領域 → FFT が既定で、式(2.67) をそのまま解く版は参照実装） |
| 残響時間算出（ipls2rt_fortran.f90） | **実装済**（`reverberation.py`。EDT / T20 / T30。帯域分割は Butterworth（IEC 61260）に変えた） |
| OpenMP 並列化 | **していない。**代わりに NumPy のベクトル化で通し 10.7 倍にした（`docs/高速化の説明.md`）。numba は入れていない |

**元コードから意図的に変えたところ**（変更前に `docs/` の該当節を読むこと）

| どこ | 何を変えたか | なぜ |
|---|---|---|
| `loop_noredundancy.image_sources` | 鏡像の距離に `abs()` を使わず符号付きにした | **判断待ち**。`docs/議論_鏡像の式のabs問題.md` |
| `reverberation` | 巡回畳み込みの自作 FIR → Butterworth（IEC 61260） | 元のやり方だと回り込みで減衰曲線に床ができ T30 が測れない |
| `impulse` | 信号処理を scipy に（`rfft`/`irfft`/`firwin`） | フィルタ長を実用的にした。**元コードにあった 1.49 秒の遅れは無い** |
| `atmosphere` | 音速と空気吸収を温度・湿度・気圧から計算（ISO 9613-1） | 定数・近似式では条件を変えられない |
| 交差判定 | **同一平面パッチ**単位・配列演算 | 速度 2.3〜3.0 倍＋**受音経路が 2.2 倍**（`docs/技術説明書.md` 5.4.1 節） |
| 容積 | 巻き順ではなく**法線**から（発散定理） | 巻き順が一貫しないモデルでも正しい（同 3.4 節） |
| 吸音材の指定 | レイヤ名そのものではなく**条件表（xlsx）**から | CAD を触らずに材料を差し替えるため（同 9.5.0 節） |

---

## GUI とプロジェクト（2026-08-15 新設）

**1 つのウィンドウで完結する GUI が最終形**だが、どんな情報が要るかがまだ固まっていないので、
いまは「必要なウィンドウをそのつど開く」形にしてある。
各ウィンドウは独立して呼べるので、統合するときは `app.py` を差し替えれば済む。

```
app.py  ──┬─→ setup_window.py   条件入力（tkinter。依存を増やさないため）
          ├─→ face_editor.py    面の確認・修正（法線と吸音材。PyVista）
          ├─→ run_project.py ──→ procedure.process()
          │                      └→ 結果 CSV ＋ plots.py で図 PNG
          └─→ view_rays.py      可視化（音線 ↔ 音粒子を Tab で切替）
```

### project.py ― プロジェクトの保存と読み込み

一度作った条件と結果を、あとから開き直せるようにするためのもの。
条件を JSON、結果を CSV、図を PNG で**プロジェクトフォルダに全部置く**。
CSV にしてあるのは Excel でそのまま開けるようにするため。

```
プロジェクトフォルダ/
  project.json          条件（DXF・吸音率・音線数・受音球・温度湿度…）
  normals.json          法線の反転指定（face_editor.py が書く）
  materials.json        面ごとの吸音材の割り当て（同上。レイヤで分けられないモデル用）
  条件表.xlsx           レイヤー名 → 材料番号（**CAD を触らずに材料を変える入力**）
                        シート「吸音率」＝材料一覧 ／ その他のシート＝条件 1 つ
                        （**シート名が条件名**になり結果ファイル名に入る）
  視点.json             画角（見る向き）と左パネルの数値・タブ（`view_camera.py`）
                        **条件を変えても同じ向きで見比べる**ための入力。
                        画面左の「画角を保存」「画角を読込」で出し入れする。
                        **名前を付けて何本でも**入る（`views`。昔の形も読める）
  結果/                 ← **受音点に依らないもの**は直下
    研修室_条件A_まとめ_残響時間.csv   全受音点 ＋ 平均 ＋ 理論値（summary.py）
    研修室_条件A_まとめ_明瞭度.csv     全受音点 ＋ 平均
    研修室_条件A_吸音率と理論値.csv    材料別の吸音率 → 平均吸音率 → 残響時間理論値
    研修室_条件A_raylog.npz            音線軌跡（可視化用）
    研修室_条件A_まとめ_音圧レベル.csv  帯域別 Lp（合計・A 特性・直接音・反射音）
    研修室_まとめ_条件比較.csv          **条件を横に並べた比較**（一括計算のとき）
    研修室_条件A_まとめ_STI.csv         STI と帯域別 MTI
    研修室_条件A_結果一式.xlsx          体裁用（グラフ付き。`workbook.py`）
    研修室_測定点.csv                   音源・受音点の座標と正面方位（**条件にも依らない**）
    研修室_条件A_計算情報.csv           使った PC と所要時間（GUI から回したとき）
    研修室_条件A_計算ログ.txt           計算中の画面のログ（同上）
    rec1/  研修室_条件A_{pulses,ir,rt,decay,clarity,spl,sti}.csv  ← **受音点ごと**
           研修室_経路.npz             経路の幾何（**条件をまたいで使い回す**。F-9）
    rec2/ …
  図/
    rec1/  研修室_条件A_{impulse_response,decay,reverberation,clarity,
                         absorption,pulses,direction,modes,mode_buildup}.png
    rec2/ …
    研修室_points.png                   測定点の配置（平面＋立面 2 方向。**条件に依らない**）
    画面/  画面から手で撮った画像・動画（`clear_results` で消えない）
```

★★**受音点に依らない結果は「1 点目の掃除」だけが消す**（2026-08-23 に直した）。
`clear_results()` は受音点ごとに呼ばれるので、`結果/` 直下まで毎回消していると
**1 点目が書いた『吸音率と理論値.csv』が 2 点目の掃除で消える**
（統計残響式は受音点に依らないので 2 点目以降は書き直されない）。

★**受音点ごとのものは `結果/recN/`、受音点に依らないものは `結果/` 直下**
（2026-08-21 にこの形へ。それまでは 1 点目だけ `結果/` 直下、2 点目以降が
`rec2/結果/` という不揃いな置き方だった）。

★**ファイル名の頭に「対象室＋条件名」を付ける**（2026-08-21 ユーザー要望）。
プロジェクト名（`name`。既定はフォルダ名）を `Project.prefixed()` が付ける。
結果は報告書やメールでフォルダの外へ出るので、`rt.csv` のままだと
どの室・どの条件か分からなくなる。**図（PNG）にも同じ頭を付ける。**

| メソッド | 役割 |
|---|---|
| `result_path(key)` | **書き出し先**。頭を付けた今の名前 |
| `result_candidates(key)` | 読める名前を探す順に返す（今の名前 → 頭なし → 昔の名前） |
| `existing_result_path(key)` | **読む先**。上のうち実在するもの。無ければ今の名前 |
| `figure_path(name)` | 図のパス（頭を付ける） |

`clear_results()` は `result_candidates()` を全部消す。**昔の名前のファイルが
残っていると今回の結果と並んでしまう**ため。

#### 『吸音率と理論値.csv』（`write_room_csv` / `read_room_csv`）

元は `rt_statistical.csv`（統計残響式）と `surface.csv`（材料別の面積・吸音率）に
分かれていたが、**どちらも受音点に依らない室の性質**で、見るときは
「材料の吸音率 → 平均吸音率 → 理論値」と続けて追うので 1 枚にした
（2026-08-21 ユーザー要望。並び順もユーザー指定）。

```
区分,項目,面積_m2,63,125,…
材料別の吸音率,コンクリート,120.5,0.02,…     ← 材料ごとの α（乱入射に変換したあと）
材料別の吸音率,吸音板,45.0,0.15,…
平均吸音率,平均吸音率,340.2,0.153,…         ← ᾱ = Σ Sᵢαᵢ / S（面積で重み付け）
平均吸音率,等価吸音面積_m2,,52.1,…          ← A = S・ᾱ
平均吸音率,空気吸収_4mV_m2,,0.0,…
残響時間理論値,sabine_s,,1.91,…
残響時間理論値,eyring_s,,…
残響時間理論値,eyring_knudsen_s,,…
```

**周波数は横**（`table.py` の共通ルール）。1 行の意味が区分ごとに変わるので
`write_frequency_table` は使わず、1 列目に区分を立てて書く。
等価吸音面積を並べているのは、理論値が `T = 24 ln10 · V /(c·A)` の A から
出ていることを表の中で追えるようにするため。

読むのは `read_room_csv()`。**1 枚にまとめる前の 2 ファイルも読める**
（`_load_legacy_room`。作り直す前のプロジェクトを開き直せるように）。

- DXF と吸音率 CSV は**プロジェクトフォルダからの相対パスで持つ**
  （フォルダごと別の端末へ移してもそのまま開ける）。外にあるものは絶対パスのまま
- `DEFAULTS` に無いキーは保存されない。**新しい計算条件を足したらここにも足すこと**
- `normals.json` / `materials.json` は**書いたときの DXF と面数を控えてある**。
  食い違ったら使わずに知らせる（黙って間違った面を反転・差し替えしないため）
- `materials.json` は「**材料名 → 面番号の配列**」で持つ。逆向きより短く、
  人が開いたときに「どの材料をどこに貼ったか」が読み取れる

### 同一平面パッチ単位の交差判定（F-4。2026-08-19）

**三角形に割るのは「ねじれのない面」を保証するためであって、判定の都合ではない**
（ユーザー指摘）。同一平面に並んだ三角形は、まとめて 1 枚の多角形として判定してよい。

`mesh_method.PatchArrays` を新設し、`collision_arrays()` が既定でこれを返すようにした。
`FaceArrays` は残してある（一致確認の基準）。

#### パッチの作り方（`coplanar_patches`）

同一平面・辺で連結（`read_dxffile.coplanar_groups`）に加えて、
**同じ材料**・**同じ法線の向き**でも割る。材料で割らないと吸音率が引けなくなる
（視聴覚室では扉・窓が壁の帯と同じ平面で接しており、同一平面だけでまとめると
18 パッチで材料が混ざった）。

★**同一平面のしきい値は編集画面より厳しくする**（`PATCH_ANGLE_DEGREES` = 0.1°、
`PATCH_DISTANCE` = 0.1 mm）。編集画面用の 1°／1 mm は「壁 1 枚」を選びやすくするための
値で、判定に使うと**形が変わる**。実測（階段教室）：1°／1 mm だと代表面の平面から
最大 3.5 mm ずれた三角形まで同じパッチに入り、交点距離が 3.15 mm 動いた。
0.1°／0.1 mm ならずれは 0.7 µm で、パッチ数は 68 → 70 とほとんど変わらない。

#### 判定の中身

平面との交点を 1 回求め、**外周多角形の内側かを交差数の偶奇で数える**。
外周は輪に並べ替えず**辺の集まりのまま**扱うので、穴の開いた面もそのまま正しい。
`hit_id` は **三角形のインデックス**（そのパッチの代表面）を返すので、下流は
パッチを知らなくてよい。代表面を返してよいのは、パッチの中では法線も材料も同じだから。

★音線追跡とバックトレースで**必ず同じ入れ物**を使うこと。反射面の番号を突き合わせる
ので、片方だけパッチにすると経路が全部却下される。

#### 効果 ― 速度だけでなく、**取りこぼしていた経路が戻る**

| モデル | 三角形 → パッチ | 判定の速度 |
|---|---|---|
| ModelTest | 68 → 16 | 223 → 75 ms（3.0 倍） |
| 階段教室 | 272 → 70 | 546 → 213 ms（2.6 倍） |
| 視聴覚室 | 292 → 114 | 709 → 303 ms（2.3 倍） |

いずれも三角形版と**距離・交点がビット単位で一致**する。

通しで回すと、速度以上に効いたのは**受音できる経路が増えたこと**。

    視聴覚室 / 音線 30,000 本 / 最大反射 40 回
      三角形  62.4 秒  経路 14633 本 → 受音  314 本  エネルギー 193.94
      パッチ  31.2 秒  経路 14234 本 → 受音  703 本  エネルギー 253.59

バックトレースは「記録した反射面の並びを再現できるか」で経路を検証する。
三角形単位だと、**同じ壁の隣の三角形に当たっただけで却下**されていた。
鏡像は平面で決まるのでこれは取りこぼしで、パッチ単位にすると正しく通る。
経路が 2.2 倍になり、後期の経路が薄いという問題（T30 の曲率）にも効くはず。
候補経路がわずかに減る（14633 → 14234）のは、重複削除がパッチ単位で効くため。

### 三角形分割の質 ― 耳の選び方（2026-08-19）

**枚数は減らせない。** n 角形は必ず n-2 枚になり、これは単純多角形の三角形分割の
理論的な最小値。変えられるのは**形**だけ。

以前の `_ear_clip()` は「最初に見つかった耳」を切っていた。走査が毎回同じ側から
始まるので**扇状に分割**され、細長い三角形（スリバー）が並ぶ。
ユーザー指摘「1 面なのに細い三角形要素で 10 分割以上されて無駄に感じる」がこれ。

`triangle_min_angle()` を新設し、**候補の耳のうち最小角がいちばん大きいものを切る**
ようにした。四角形の対角線も、両方使えるときは形の良いほうを選ぶ。

| 形 | 扇状（旧） | 質で選ぶ（新） |
|---|---|---|
| 10:1 の細長い 12 角形 | 最小角 1.53° | **14.76°** |
| 階段教室 272 枚 | 中央値 9.2° / 1°未満 24 枚 | **中央値 11.6° / 1°未満 16 枚** |

**細長さは実利。** レイと三角形の交差判定は、面が潰れているほど浮動小数の丸めに
左右されやすい。形が良くなれば当たり外れが安定する。

★**残るスリバーはモデル側の形状。** 視聴覚室の反射板は 1 枚が
**14.375 × 0.024 m（縦横比 610:1）**で、これを三角形 2 枚に割れば必ず 0.09° になる。
分割の工夫では直せない（25 mm 厚の板をそう描いてあるため）。
どこまでがモデルの性質かは `DxfModel.triangle_quality`（最小角の中央値・最小・
1°未満と 5°未満の枚数）を見て判断する。

### T 字接合は問題になるか（2026-08-19）

壁を高さの帯や開口で分割したモデル、こちらの変換でスラブ分割した面などでは、
辺が 1 対 1 で合わず `open_edge_count()` が 0 にならない（T 字接合）。
**実測では問題にならなかった。**

- **本当の隙間は 0 本。** 「開いた辺」1 本ずつについて、同じ直線上の他の開いた辺で
  覆われているかを調べたところ、階段教室 16 本・視聴覚室 70 本すべてが覆われていた。
  面は連続していて、辺の分け方だけが違う
- **容積が AutoCAD の `MASSPROP` と一致する。** 隙間があれば発散定理が破れて値がずれる
  ので、これが watertight であることの裏付けになる（階段教室 3557.3 m³、ModelTest 770.33 m³）
- レイが抜けるのは「T 字の頂点をちょうど通る」場合だけで、幅ゼロ＝測度ゼロ

**効くのは別の 2 点。** どちらも対処済み。

1. `is_closed` が False になるので、以前は容積が出せなかった → 法線から出すようにした
2. `orient_normals='shells'` が使えない（巻き順の一貫性を要求する）
   → 既定の `'auto'`（面ごとにレイの偶奇で判定）は巻き順に依存しないので影響なし

### 総表面積と容積 ― 法線から出す（2026-08-19）

「音源・受音点を入れる前に、容積と総表面積だけ知りたい」という用途のために、
`read_model()` が読んだ時点で両方を持つようにした。**幾何だけで済む**ので
プロジェクトを作る前に使える。

    cd geosim
    python read_dxffile.py "C:\...\室.dxf"

出るもの：面数 / レイヤ別の枚数 / **総表面積** / **レイヤ別の面積** /
**容積と平均自由行程** / シェル内訳 / 作図チェック。

#### なぜ巻き順ではなく法線から出すのか

`signed_volume()` は `V = (1/6)Σ x1·(x2×x3)` で頂点の**巻き順**から体積を出す。
これは巻き順が一貫した閉曲面でしか成り立たない。ところが実際のモデルでは崩れる。

- **リージョン（ACIS）から輪郭をつないで作ったモデル**は面ごとの巻き順がまちまち
- **壁を高さの帯や開口で分割**すると辺が 1 対 1 で合わず（T 字接合）、
  `open_edge_count()` が 0 にならないので `signed_volume()` を使う条件を満たさない

実例（視聴覚室モデル）：反射板 1 枚（14.375 × 0.91 × 0.025 m）の体積が
**巻き順だと 0.1090 m³、正解は 0.3270 m³**。室全体の容積も `None` になっていた。

そこで `volume_from_normals(triangles, normals)` を新設した。
発散定理 `V = (1/3)∮ r·n dA` を**法線**で評価する。利点は 2 つ。

1. **巻き順に依存しない。** 法線が空気側を向いていれば値が合う。T 字接合でもよい
2. **家具・反射板の体積が自動で引かれる。** 宙に浮いた物体の法線は物体の外
   ＝空気側を向くので、その寄与が逆符号で入り、外殻の分から差し引かれる

★**面に穴が開いていると誤る。** 見分けるのが `uncovered_open_edges()` で、
「開いた辺」を 2 種類に分ける。

| 種類 | 例 | 面としては | 容積 |
|---|---|---|---|
| **T 字接合** | 壁を高さの帯で分割した継ぎ目。長い辺 1 本 対 短い辺 3 本 | 閉じている | **使える** |
| **本当の自由端** | 宙に浮いた片面の板の外周、面の抜け | 開いている | **目安** |

判定は「同じ直線上の他の開いた辺に覆われているか」。実測：階段教室 16 本・
視聴覚室 70 本の開いた辺はすべて覆われていた（＝T 字接合）。
研修室（板の寄せ集めモデル）は 136 本のうち 14 本・計 32.6 m が自由端で、
`volume_note` と `check_model` が「容積は目安」と知らせる。

そのうえで `read_model()` は
`encloses_point()` で囲まれているか確かめてから採用し、駄目なら `volume` を
`None` にして理由を `volume_note` に入れる。巻き順から出した値と 1% 以上
食い違ったときも `volume_note` で知らせる（**法線のほうを採用する**）。

シェルごとの体積も読み込み後に法線から出し直す（`analyse_shells` が返すのは
巻き順から出した値で、崩れていると合わないため。巻き順の値は `volume_winding` に残る）。

検算：視聴覚室モデルで **総表面積 1094.829 m² / 容積 1152.269 m³**。
「床面積 256.776 × 高さ 4.5 − 反射板 3.224」= 1152.268 m³ と一致した。
ModelTest では 770.328 m³ で、AutoCAD の `MASSPROP`（770.33 m³）と一致。

### face_editor.py ― 面の確認・修正（法線と吸音材。G-8 / G-9）

読み込み時の自動補正（`orient_normals='auto'`）やレイヤ→吸音材の対応に任せきりにせず、
**人が見て確認し、必要なら直せる**ようにするためのもの。
2026-08-19 に `normal_editor.py` から改名し、吸音材の割り当てを足した。

| 保存先 | 中身 | 渡し先 |
|---|---|---|
| `normals.json` | 法線を反転する面 | `read_model(flip_faces=...)` |
| `materials.json` | 面ごとの吸音材 | `read_model(face_materials=...)` |

#### 面グループ（同一平面パッチ）― この画面の要

3DSOLID を STL 経由で取り込むと、**設計者が描いた 1 枚の壁が三角形に割られる**。
そこで `read_dxffile.coplanar_groups()` で**同一平面かつ辺で連結した三角形をまとめ**、
それを選択・操作の単位にする。実例：ModelTest は三角形 68 枚 → **16 グループ**
（床・天井・壁 14 枚）で、元のソリッドが持っていた面と一致した。

まとめる条件は「辺を共有」かつ「法線が 1° 以内・平面からのずれが 1 mm 以内」。
法線の比較に**絶対値**を使うのは、**グループ分けを法線の向きから独立させる**ため
（向きは `flip_faces` で変わるが、グループは幾何の性質なので変わってはいけない）。
`y` で三角形単位にも切り替えられる。

#### 選ぶ → 適用する（二段階）

以前は枠で囲んだ瞬間に反転していたため、**何を選んだのか分からなかった**（ユーザー指摘）。
いまは選択が残り、明るい橙で塗られ、**外周が線で描かれる**（中の三角形の辺まで引くと
網目になって形が読めないので、1 回しか現れない辺だけを残す）。
左パネルに「1 グループ / 16 面 / 256.8 m²」のように枚数と面積が出る。

一括選択：`k` 同じ向きの面（床・天井・壁の一括）/ `l` 同じ吸音材の面 /
`j` 全選択 / `h` 選択反転 / 数字キーでレイヤ / `0` 解除。

適用：`i` 選択の法線を反転（選択が空なら全部）/ 左パネルの材料をクリックで吸音材を貼る。

| 色（法線モード） | 意味 |
|---|---|
| 緑 | 法線が室内（空気側）を向いている |
| 赤 | 法線が室外を向いている。**反転が要る** |
| 灰 | 判定できない（開いた形状など）。CAD の指定を尊重する |
| 橙 | いま選択している面（判定・材料の色より優先） |

`m` で表示を **法線の判定 → 吸音材（材料ごとの色）→ 容積の拾い方** と切り替える。

#### 「容積の拾い方」モード（2026-08-19）

**どの面を拾って、どの領域を容積として数えたか**を目で確かめるためのもの。

| 色 | 意味 |
|---|---|
| 青 | 外殻。**この面が囲む中身が容積になる** |
| 橙 | 内側の閉じた物体（家具・反射板）。**その体積は容積から引かれる** |
| 灰 | 内側の開いた板。体積には寄与しない |
| **赤い線** | **開いた辺**。穴なのか T 字接合なのかがここで分かる |

開いた辺を線で重ねるのが要点。視聴覚室モデルの「開いた辺 70 本」は、
赤い線が**壁の帯の継ぎ目と床際に沿って走る**ので、穴ではなく面の割り方が
違うだけの T 字接合だと一目で判断できた（`open_edge_segments()`）。
左パネルには総表面積・容積・シェルの内訳も出る。

面の色は cmap ではなく **RGB を直接持たせている**（モードで色の意味が変わるので、
cmap を差し替えるより配列を書き換えるほうが単純で確実）。

状態は `flipped`（**CAD の巻き順から反転する面のインデックス集合**）と
`assigned`（**{面インデックス: 材料名}**）の 2 つだけ。どちらも `read_model` に
そのまま渡せる形にしてある。自動判定の結果も `flipped` に畳み込んで保存するので、
**保存したものを読めば再現できる**。適用の順序は **自動判定 → 手動指定**。

面ごとに吸音材を割り当てると `Mesh.material` が材料名になるので、
**『吸音率と理論値.csv』が材料別の面積・吸音率の表になる**（レイヤ別ではなく）。
DXF の実際のレイヤ名は `DxfModel.face_layers` に別に持っている
（画面のレイヤ表示が材料名に化けないように）。

★キーを足すときは `view_model_gui.VTK_RESERVED_KEYS` を見ること。
`e` `q` `w` `s` `r` `f` `p` `u` `3` は VTK が使っている。
枠選択は VTK のラバーバンドで、**`r` で ON/OFF する**（撮影が `g` なのはこのため）。

### plots.py ― 結果を図にする（G-3 / G-4）

`図/` に PNG で書き出す。画面には出さず Agg バックエンドを使う
（GUI のイベントループと取り合わないため）。

**インパルス応答は最大値で正規化して保存する。**
絶対振幅の校正が未了なので、そのまま出しても値に意味が無いため（E-11）。

バンド別のエネルギー時間曲線は **5 ms の移動平均で包絡線にしてから**重ねる。
生のままだと反射音の干渉で激しく暴れて読めない。

`pulses.png` の下段は**経路数の時間分布**。拡散音場なら時刻の 2 乗で増えるはずなので、
減っていく形なら後期の経路を取りこぼしている（受音球が小さすぎる）と分かる。

### plots.pulse_spectrum() / plots.mode_buildup() ― 経路差から見たモード（G-6b）

パルス列を**位相込みで足す**。`H(f) = Σ a_n exp(-i 2π f t_n)`。
到来時刻の差がそのまま経路差なので、`exp` の中身は経路差ぶんの位相。

**`a_n = 1`（減衰なし）のとき `|H(f)|` は「その周波数で同位相に重なった経路の
本数」そのもの**になる。1 本なら 1、2 本重なれば 2。テストで確かめてある。

経路差 Δ の音は Δ が波長の整数倍のとき強め合うので `f = m·c/Δ`。
直方体の x 方向の往復は Δ = 2Lx となり `f = m·c/(2Lx)`、これは軸モード
`f = (c/2)(m/Lx)` に一致する。`room_modes()` の式そのもの。

重みを 3 通りにして比べる。**パルスのエネルギーには面の吸音しか入っていない**
（距離減衰 `1/(ct)` と空気吸収は `impulse.transfer_function` /
`apply_air_absorption` で掛ける約束）ので、この分離ができる。

| 重み | 意味 |
|---|---|
| `1` | 減衰なし。**室形状だけ**で決まる。吸音材を貼っても変わらない |
| `1/d` | 完全反射（距離減衰のみ） |
| `√E/d` | 設計の吸音あり |

後ろ 2 本の差が**吸音の効果**で、`mode_buildup.png` の下段で塗ってある。

**速さ**：素直に書くと「パルス本数 × 周波数点数」の行列積になるが、到来時刻を
細かい格子に投げ込んで（`np.bincount`）から FFT すれば O(n log n)。
要はインパルス応答を作るのと同じ手順。格子は既定 8 kHz（0.125 ms = 4.3 cm）で、
200 Hz での位相誤差は最大 4.5°。

**周波数の刻み**：FFT の刻みは応答長で決まる（3 秒なら 0.333 Hz）。細かすぎて線が
ぎざぎざするので、`frequency_step`（既定 **1 Hz**）でまとめ直す（`rebin_spectrum`）。
**間引くのではなく二乗平均する**。減衰なしの①は櫛の歯が鋭く、幅は応答長でしか
決まらないので、点を飛ばすと山を跨いで見落とす。二乗平均なら取りこぼさず、
**ランダム位相の目安 √N が帯域幅に依らず保たれる**ので基準線を引き直さずに済む。

**図を読むときの注意**
- 位相がばらばらなら和はランダムウォークで √N 程度にしかならない。
  図に √N の線を引いてあり、**超えている山だけが本当に強め合っている周波数**
- 0 Hz は全経路が同位相になるので左端に巨大な山が立つ。
  **いちばん低い固有周波数から下は描かない**
- 吸音率はオクターブバンドしか無いので、各周波数点には最も近いバンドの値を当てる
  （`_band_of`）。6 バンドだと 88 Hz 以下は全部 125 Hz の値になる
- パルス列は受音球が拾えた経路の**標本**。本数の絶対値は音線数と受音球半径で変わる

### run_project.redraw() ― 計算し直さずに図だけ作り直す

**音線追跡（重い部分）はやり直さない。** プロジェクトフォルダに残っている
`pulses.csv` と `ir.csv` を読み、そこから先だけを計算して `図/` を作り直す。
研修室（パルス 3901 本）で数秒、図 9 枚。

    cd geosim
    python run_project.py "C:\Users\...\JR" --redraw

図の描き方を直したあとに**過去のプロジェクトへ反映する**ための道。
GUI の「前回の結果を見る」も内部でこれを呼ぶ（以前は文字の要約だけで、
図は古い見た目のまま残っていた）。

**パルス列とインパルス応答はそのまま使う**（再合成しない）ので、前回の計算結果と
食い違わない。残響指標・明瞭度・統計残響式は**本番と同じ関数**で計算し直すため、
CSV の読み方をここに別途書かずに済む（実測で保存済み CSV と 1e-12 s 以内で一致）。
**結果 CSV は書き換えない**（テストで確認している）。

### table.py ― 表の並べ方の共通ルール（2026-08-17）

**周波数は「横」に並べる。** グラフにしたときの横軸が周波数なので、表もそれに揃える。
縦だと Excel でグラフを作るたびに行と列を選び直す手間が要る（ユーザー判断）。

    ○ 項目,125,250,500,1000,2000,4000
      EDT_s,0.945,0.571,0.459,0.295,0.267,0.284

    × frequency_hz,EDT_s,T20_s
      125,0.945,0.633

| ファイル | 1 行の意味 | 周波数の位置 |
|---|---|---|
| `rt.csv` / `clarity.csv` / `まとめ_*.csv` | 指標 | **列**（`write_frequency_table`） |
| `decay.csv` | 時刻 | 列（`decay_125Hz_db` …） |
| `吸音率と理論値.csv` | 材料／指標（1 列目が区分） | 列（`project.write_room_csv`） |
| `pulses.csv` | 経路 | 列（`energy_125Hz` …） |
| `ir.csv` | 時刻 | 周波数の軸を持たない |

- 1 行が指標なら `table.write_frequency_table()`。1 列目が `項目`、以降が周波数
- 1 行が周波数以外なら列名を `table.band_column()` で作る。
  **番号（`energy_1`）にしない**。バンド数 6 と 8 で意味が変わって読み違えるため
- 読むのは `table.read_frequency_table()`。**古い縦向きの CSV もそのまま読める**
  （向きを変える前に作ったプロジェクトを開いても壊れないように）
- Excel でそのまま開けるよう **BOM 付き UTF-8**。NaN は空欄で書く
- **画面に出す表も同じ向き**にしてある（`procedure` / `app` / `reverberation`）

**新しい出力を書くときもこの向きにすること。**
どうしても縦にしたい事情が出たら、勝手に決めずにユーザーに確認する。

### キーの割り当て ― VTK の予約キーを避ける（2026-08-17）

**`e` と `q` は VTK の終了キー。** ここに機能を割り当てると、その機能と一緒に
**終了処理（ExitEvent）まで走る**。

数値入力を `e` に割り当てていたため、**値を入れて OK を押した瞬間に
3D ウィンドウが閉じ、条件入力の画面に戻っていた**（ユーザー報告）。
`run_animation` が `ExitEvent` を見てループを抜けるので、
例外も出ずきれいに閉じてしまい、ログにも痕跡が残らなかった。

数値入力は **`t`** に変えた（`view_model_gui.VALUE_INPUT_KEY`）。

実測した予約キー（`InvokeEvent('CharEvent')` で再現できる）:

| キー | VTK の既定動作 |
|---|---|
| `e` / `q` | **終了（ExitEvent）** |
| `w` / `s` | ワイヤフレーム / 面 |
| `r` | 視点リセット |
| `f` | 注視点へ寄る |
| `p` | ピック（面の選択） |
| `u` | ユーザーイベント |
| `3` | ステレオ表示 |

一覧は `view_model_gui.VTK_RESERVED_KEYS`。`enable_value_input()` は
予約キーを渡されたら警告を出す。**テストでも押さえてある**（新しいキーを足すときは
`tests/test_geosim.py` の [22] にも足すこと）。

※ `face_editor` が `s`（保存して閉じる）を使っているのは承知のうえ。枠選択は VTK の `r`。
  `p` はピックそのもの、`s` は押しても面表示に戻るだけで実害がない。

### 絞り込みはスライダで指定する（G-16b、2026-08-19）

**指摘**：「スライダーの方では無く、モデル上で `p` を押して、という操作が難しい。
どの方向へ、というのも平面方向と縦方向のスライダーを用意するなどで
指定出来た方がやり易い」。

そのとおりに変えた。**クリックは任意**になり、スライダとキーだけで使える。

| 決めるもの | どう決めるか |
|---|---|
| 種類 | `k`（なし → 近くを通る → この方向に飛ぶ → 1 本だけ） |
| 基準点 | `j`（**受音点** → 音源 → 拾った点）。**既定は受音点** |
| 方向 | **「方位角 [°]」「仰角 [°]」スライダ** |
| 範囲 | 「近さ [m]」「方向の半角 [°]」スライダ |

- **方向の約束は `project.head_azimuth` と同じ**。方位角 0° = +X で反時計回り、
  仰角は 0° が水平・+90° が真上（`ray_filter.direction_from_angles`）。
  バラバラだと混乱するので揃えた
- **方向は音源から伸びる黄色い矢印で描く**。スライダを動かすと向きが変わるので、
  何を指定しているかが見える
- 方向の初期値は**音源から受音点を見る向き**（いちばん見たい向きのはず）
- `p` で点を拾ったときは、方向モードなら**スライダの値もその向きに合わせる**
  （`angles_from_direction`。拾う操作とスライダで状態が食い違わないように）
- そのために `ControlPanel.slider` の control に `show` を持たせた
  （外から値を入れ直したとき、**つまみだけ動いて数字が古いまま**にならないように）

#### パネルは「別画面」として背景を変える（2026-08-21）

> モデル表示部と左側で、背景を変えておいて欲しい。あくまでも別画面扱いで、
> 左側は設定画面とわかるような感じの背景色にして

| | 背景 | 見え方 |
|---|---|---|
| 左（設定・操作） | `PANEL_BG = #0d1220` 平ら | 濃紺の面。UI の側だと分かる |
| 右（モデル） | `#1c2027` → `#2e3540` のグラデーション | 上が明るい、絵の側 |

境目に細い線（`PANEL_EDGE`）を入れ、パネルの見出しも色を変える
（`screen_title()`。「設定・操作（左）／モデル（右）」の副題つき）。

★★**`plotter.set_background()` は既定で全レンダラに掛かる**
（`all_renderers=True`）。片方だけ変えるには `all_renderers=False` が要る。
これを忘れると後の呼び出しでパネル側が上書きされ、
**色を変えたつもりで何も変わらない**（実際に 1 回それで気づかず進めた。
画面を撮って画素の色を測って初めて分かった）。

#### パネルに入りきらない問題（2026-08-19 → 2026-08-21 に解決）

スライダが 9 本になり、レイヤ 8 つのモデルで**下が 458 px はみ出した**（操作説明ごと消える）。
2026-08-19 は「載せる内容を削る」で 268 px まで減らしたが、実案件
（研修室：レイヤ 9 種 ＋ 材料 19 種）では **1280×860 で 450 px、1920×1080 でも
230 px あふれた**（ユーザー指摘 2026-08-21）。削り切れないので**送れるようにした**。

**① ページ送り（`enable_scroll()`）**

`relayout()` が `self.scroll` のぶん下へずらして積み直す。PageUp / PageDown で
1 ページ（見えている高さの 8 割）ずつ送る。下端に固定した案内（`_hint`）に
送れることと位置（%）を出す。**案内は積まないので送られない。**

★2026-08-19 に取り下げた原因は「**置き直しを飛ばした**」こと。古い位置に
ウィジェットが残って文字に重なっていた。いまは**必ず全部を置き直し、
範囲外は `show(False)` で隠す**。文字はレンダラのビューポートで切られるが、
**スライダはウィンドウ全体の座標で描かれるので隠さないと 3D の上に残る**
（`PanelItem.show`）。

**② 操作の一覧を別ウィンドウ（`help_window()` / F1）**

パネルに全部書くと縦を食うので、要点だけ残して**全部は F1 の別ウィンドウ**で読む
（`enable_value_input` と同じく tkinter の窓。読むだけなので modal でよい）。
これで面の確認画面は 1282 → 1218 px になった。

**③ それでも入りきらないときはコンソールにも出す**（`hidden_height()` で判定）。
パネルの高さ＝ウィンドウの高さなので、縦に広げれば入る。

### 終了処理で落ちないようにする（2026-08-19）

**症状**：5 分ほど操作したあと、アプリの終了時に segfault（exit 139）。
`[app] 閉じました` の**あと**なので、計算結果と図は無事。
GL のエラー（`Could not set shader program` 多数、`Error binding ndCoords to VAO`）が
**すべて終了時刻に固まっている**＝ OpenGL の文脈が消えたあとに VTK が描こうとした形。

**再現できていない**（対話ループ・スライダ連打・絞り込みの切り替えを 3 通り試して全部正常終了）。
GL ドライバの状態にも依るらしい。以下は「落ちる余地を減らす」処置で、根治の確認は取れていない。

1. `run_animation` は**ウィンドウが無くなったら 1 コマも描かずに抜ける**
   （`_closed` が遅れることがあるので `plotter.iren` まで見る）
2. `view_model_gui.release_window()` … 閉じる前に**こちらが握っている actor の参照を外す**。
   VTK の GL 資源は文脈が生きているうちに解放する必要があるが、
   Python の後片付けは順序が決まっていない
3. `app.py` は終了直前に `pv.close_all()`（`finally` で必ず通す）

**また落ちたときの手がかり**：× で閉じたか `q` で閉じたか / 最大化したか /
数値入力（`t`）を使ったか。

### ray_filter.py ― 注目したい音線だけを残す（G-16、2026-08-19）

**きっかけ**：「この経路をもっと見たい」「この辺りの粒子をもっと見たい」
「この方向に飛ぶ音線・粒子に注目したい」（ユーザー要望）。
音線を 60 本描いても室内では線が重なって追えず、増やすともっと読めない。
**数を減らす方向ではなく、見たいものだけを残す方向**で絞る。

3 つとも「表示する音線の添字を絞る」話なので、既存の作りにそのまま乗った。

| 種類 | 関数 | 何が残るか |
|---|---|---|
| 近くを通る | `near_point()` | 折れ線が基準点から半径以内を通った音線 |
| この方向に飛ぶ | `in_direction()` | **音源から出たときの向き**が基準の向きに近い音線（円錐） |
| 1 本だけ | `nearest_ray()` | 基準点にいちばん近い 1 本（反射面の並びも出す） |

おまけで `through_face()`（この面で反射した経路）も引ける。

操作は共通で **`p` で 3D 上の点を拾う**（VTK 標準のピックキー）。壁でも音線でもよい。
`k` で種類を切り替え、`0` で解除。範囲は左パネルの「近さ [m]」「方向の半角 [°]」。
拾った点には**黄色い球**を置く（どこを拾ったか見えないと操作できない）。

- `ray_filter.py` … **画面に触らない純粋な計算**。テストは GL 不要
- `view_rays.RayFocus` … その結果を `RayDisplay.set_pool()` と
  `ParticleAnimation.set_focus()` に渡し直すだけ

★**方向は「出射方向」で見る**（`RayLog.directions`）。途中の向きで見ると反射のたびに
  向きが変わり、どの音線も条件に引っかかって絞り込みにならない。

★**「近くを通る」は描いている反射回数までで測る**（`max_reflection`）。
  50 回も反射を追うとどの音線も室内を回って基準点の近くを通ってしまう。
  研修室・半径 0.5 m で **2000 本中 741 本**が該当した。反射 4 回までなら **67 本**。

★**距離は線分まで測る**（節点までではない）。節点だけで見ると、長い区間の途中を
  通っている音線を拾い落とす。

★**粒子はトポロジを変えずに隠す**。`_shown` を書き換えてエネルギーを NaN にし
  `nan_opacity=0` で消す（フレームごとにジオメトリを組み直さないための既存の約束に従う）。

★**スライダの本数（wanted）と実際に描いた本数（count）を分ける。**
  分けないと、絞り込みで候補が 8 本に減ったあと**緩めても 8 本のまま戻らない**
  （実際にそうなった）。

### 画面の保存とウィンドウのタイトル（G-12、2026-08-17）

**`g` でいまの画面をそのまま画像に、`b`（音粒子のとき）で動画に**する。
角度も、絞った本数も、不透明度も、見えているままが残る。

置き場は **`図/画面/`**。`図/` の直下ではなく子フォルダにしてあるのは、
`clear_results()` が `図/` の PNG を消すため。同じ場所に置くと**計算し直すたびに
手で撮った画像が巻き添えで消える**（テストで確認している）。
ファイルは `法線_01.png` `音線_01.png` `音粒子_0123ms_01.png` のような連番で、
**撮るたびに増える**（上書きしない）。音粒子は止めた時刻を名前に入れる。

動画は GIF。追加の依存を増やさないため（Pillow は matplotlib の依存で既に入っている。
mp4 には ffmpeg が要る）。そのぶん幅 720 px・128 色に落とし、コマ数は 400 で頭打ちにして
等間隔に間引く（離散化時間を 0.2 ms にすると 3 秒で 15000 コマになる）。

`save_movie()`（off-screen で作り直す既存の関数）とは**別物**。
こちらは画面に出ているウィンドウをそのまま録る。

#### ウィンドウのタイトル

**タイトルバーには画面の種類だけを出す。物件名は入れない**（ユーザー判断 2026-08-17）。
知りたいのは「法線の確認なのか、音線の可視化なのか」だけで、物件名は
画面の中の見出しに出ている。題は `WINDOW_TITLES` に集めてあり、
画面が増えたら 1 行足すだけで済む。

| `screen` | タイトル | 化けたときの英字 |
|---|---|---|
| `normals` | 法線の確認 | `geosim - normals` |
| `rays` | 音線・音粒子 | `geosim - rays and particles` |
| `directions` | 音線の飛び方 | `geosim - ray directions` |
| `model` | モデルビューア | `geosim - model viewer` |

`make_plotter(title, ..., screen=...)` の `title` は**画面の中の見出し**
（物件名を含んでよい。フォント指定で描くので化けない）、
`screen` が**タイトルバー**を決める。分けてあるのはこのため。

#### タイトルが化ける件

VTK が Windows に作るウィンドウは **ANSI ウィンドウ**で、日本語のタイトルが
**文字によって**化ける。`研修室 — 法線の確認` が `研修室 ?E法線?E確?E` になった
（研・修・室・法・線・確は通るのに —・の・認 が化ける）。手は一通り試して全滅:

| やったこと | 結果 |
|---|---|
| `SetWindowTextW`（ワイド文字版） | **かえって全部 `?`** |
| `SendMessageW` で `WM_SETTEXT` | 同上 |
| `SetWindowTextA` に CP932 のバイト列 | やはり全部 `?` |

そこで `set_window_title()` は**設定したあと読み戻して検証**し、
一致しなければ `geosim - normals [JR]` のような英字の題に差し替える。
化けた文字列を出しておくよりましという判断。**画面の中の見出しは日本語のまま**
（フォントを指定して描いているので化けない）。

踏んだ落とし穴が 2 つある。

★**タイトルの設定を描画イベントの中でやってはいけない。** 最初はそうしていたが、
  `SetWindowName` は `WM_SETTEXT` を同期で送るため描画中にウィンドウプロシージャが
  再入し、**OpenGL のコンテキストごと落ちた**（segfault、以後そのマシンの新しい
  プロセスで `wglChoosePixelFormatARB` に失敗する状態になった）。
  いまは `finish_window()` を **`show()` の直前**に呼ぶ形にしてある。

★**決めた題は `plotter.title` にも書き戻す。** pyvista の `show()` は
  そのつど `self.title` をウィンドウ名に**貼り直す**ので、書き戻さないと
  差し替えた題が日本語（＝化ける方）に戻る。実際に音線の画面だけ戻っていた。

### sound_level.py ― 音圧レベルと STI（依頼 2026-08-21）

**どちらも「受音点で実際に受け取るエネルギー」から出る**ので同じモジュールに置いた。

パルス列（`loop_noredundancy.PulseList`）のエネルギーには**反射の吸音だけ**が
入っていて、距離減衰と空気吸収はインパルス応答を合成する段階で掛けている。
そのままでは「受音点でどれだけ受け取るか」が分からないので、ここで 1 か所にまとめる。

```
received_energy()   E_i = A_i · exp(-m·d_i) / (4π d_i²)
band_levels()       Lp = Lw + 10log10(Σ E_i) + 10log10(ρc/400)
modulation_transfer()  m(F) = |Σ E_i e^{-j2πF t_i}| / Σ E_i
speech_transmission_index()  m(F) → 実効SNR → TI → MTI → STI
```

★**逆二乗則が実装の物差しになる。**反射面が無ければパルスは直接音 1 本だけで
`Σ E_i = 1/(4πd²)` なので `Lp = Lw − 20log10 d − 11` に厳密に一致する
（教科書の「−11」は `10log10 4π = 10.99` の丸め）。
`freefield_level()` がその理論値を返し、`spl.csv` と `spl.png` に
**自由音場との差**を必ず並べる。無響室のモデルなら差はゼロ、
反射があればそのぶん上に外れ、外れ量が室の効きになる。

音源パワーレベル PWL が未入力なら **Lw = 0 dB（＝1 pW）として相対値**で出す。
帯域ごとの相対関係と距離依存性はそのまま読めるので逆二乗の確認には足りる。
`result['relative']` にどちらで計算したかが入り、出力にも書く。

STI は IEC 60268-16。**インパルス応答を帯域分割して積分するのではなく、
パルス列から直接 m(F) を出す**（フィルタの遅れや裾の影響を受けない。
幾何音響では普通この形）。背景騒音と聴覚マスキングは
**PWL と騒音レベルの両方がそろっているときだけ**効かせる（絶対値が要るため）。
6 バンド（125〜4k）計算では 8 kHz が無いので**重みを外して正規化し直し警告する**。

数式は `docs/技術説明書.md` 9.8・9.9 節。検算（解析解との一致）は
`tests/test_geosim.py` [31][32]。

### path_cache.py ― 経路の使い回し（F-9。依頼 2026-08-21）

> 音線ループからバックトレースまでは、とくに吸音の違いは関係なく、最後の虚像法による
> 計算の時にようやく吸音の効果が反映されるのかなと思っています。つまり、わざわざ
> 最初から計算し直さなくても、最後の部分だけデータを残しておけば、データの流用および、
> 再計算時に、そこから再開できませんか？

そのとおり。正確には吸音が効くのは**バックトレースのエネルギー計算だけ**。

| 段階 | 吸音に依存するか |
|---|---|
| ① 音線追跡 | **しない**（受音しても打ち切らない。反射は鏡面反射だけ） |
| ② 重複削除 | しない |
| ③ バックトレース・経路の検証 | **しない**（虚音源が成立するかは純粋に幾何） |
| ③' バックトレース・エネルギー | **する**（各反射で `|R(θ,α)|²` を掛ける） |

だから経路ごとに

- **反射面の並び**（`walls` (P,K)、-1 詰め）
- **各反射の入射角の余弦**（`cos_theta` (P,K)）
- 到来時刻・距離・到来方向（これも幾何）

を `結果/recN/<室>_経路.npz` に残しておけば、吸音材を変えた再計算は

    E = Π_k |R(cosθ_k, α_{面_k})|²        （`loop_noredundancy.energy_from_geometry`）

という**配列の掛け算だけ**で済む。実測（テスト室）で丸め誤差 1e-15 で一致する。

★**1 つだけ落とし穴がある。**交差判定は同一平面パッチ単位で、パッチは
「同一平面＋辺で連結＋**同じ材料**＋同じ法線」でまとめている。つまり
**材料の割り当て方を変えるとパッチの切れ目が動き、見つかる経路も変わる**。

- 材料の**値**（吸音率）を変えるだけ → パッチは同じ。**使い回せる**
- 材料の**分け方**が変わる（隣り合う同一平面の面が同じ材料になる／別になる）
  → パッチが変わる。**使い回せない**

そこで保存時に**指紋**（三角形の頂点と法線・パッチの分け方・音源・受音点・
音線数・最大反射回数・受音球・両面判定）を残し、読むときに突き合わせる。
食い違ったら理由を告げて①からやり直す。**吸音率の値は指紋に入れない。**

`clear_results()` はこの npz を消さない（`pj.KEEP_ON_CLEAR`）。
作り直すのに音線追跡が丸ごと要るうえ、古いかどうかは指紋で分かるため。

### 条件を切り替える／一括で回す（依頼 2026-08-21）

> 条件入力画面で、どの条件で計算するかファイル選択できる、という認識ですか？
> 同一フォルダ内に、複数条件のファイルを入れる予定です。
> おそらく、DXF のファイル名は物件名や部屋名、条件表のファイル名は条件名にすると
> 思うので、それらを保存する結果ファイルなどの名前に使ってください。
> また、複数条件やる場合、一括で回せると嬉しいです。

**ファイル名がそのまま名前になる。**

    結果ファイル名の頭 ＝ 対象室名 ＋ 条件名
                        ↑              ↑
                        DXF のファイル名  材料条件表のファイル名

| | 決め方 | 実装 |
|---|---|---|
| 対象室名 | 条件入力の「対象室名」欄。**空欄なら DXF のファイル名** | `Project.room_label` |
| 条件名 | **条件表のシート名**（CSV ならファイル名。既定名なら付けない） | `Project.condition_label` |

    研修室.dxf ＋ 条件表.xlsx のシート「吸音追加案」 → 研修室_吸音追加案_rt.csv

条件表は条件入力ウィンドウの「条件表（xlsx）」で選び、条件は
「条件（シート）」のコンボボックスで選ぶ。用意が無ければ「条件表を作成」で作る。

**一括計算**は「全条件を一括 ▶▶」ボタン、または

    python run_project.py <フォルダ> --conditions              … フォルダの条件表を全部
    python run_project.py <フォルダ> --conditions 現状.csv 案A.csv

`run_conditions()` が条件を順に回す。**1 件目で音線追跡まで済ませ、
2 件目以降は経路を使い回してエネルギーだけ計算し直す**（上記 F-9）。
最後に**条件を横に並べた比較表**（`<室>_まとめ_条件比較.csv`）を作り、
条件ごとの Excel に「条件比較」シートを入れる。

★**条件に依らないファイルは条件名を付けない**（`pj.ROOM_SCOPED_RESULTS`）。
経路の npz（使い回すため）と音線軌跡（形は吸音に依らないため）の 2 つ。

### condition_table.py ― 条件表（依頼 2026-08-21）

「レイヤー名 → **材料番号**」の対応と吸音率を **1 つの Excel** で持つ。

```
プロジェクトフォルダ/条件表.xlsx

  シート「吸音率」      ← PJ で使う材料の一覧。**1000 材料ぶんの枠**
    番号,材料名,63,125,…,8000,種類,備考
     1 ,スチールドア,0.15,0.29,…,残響室法,
    （番号は `=ROW()-1`。**行の位置がそのまま材料番号**。種類はドロップダウン）

  シート「現状」        ← **条件 1 つ ＝ シート 1 枚。シート名が条件名**。30 行の枠
    番号,区分,レイヤー名,材料番号,安全率,材料名（参考）,面数,面積_m2
     1 ,レイヤ,01__研修室_壁_扉,1,,=IFERROR(VLOOKUP($D2,吸音率!$A$2:$B$1001,2,FALSE),""),2,3.3

  シート「吸音追加案」  ← シートを複製して番号を書き換えるだけ
```

#### 枠の数と表示（ユーザー指定 2026-08-21）

| | 数 | 定数 |
|---|---|---|
| 吸音率の材料 | **1000** | `MATERIAL_SLOTS` |
| 条件シートのレイヤ | **30**（式も 30 行ぶん） | `LAYER_SLOTS` |

吸音率・安全率は**小数 2 桁**、面積は**小数 1 桁**、面数は整数で表示する
（`ALPHA_FORMAT` / `AREA_FORMAT` / `COUNT_FORMAT`）。

#### ★安全率（例 0.8 掛け）

「吸音率を見過ぎないように間引く」ための係数。**未入力なら掛けない。**
★**カタログ値（吸音率シートの値）に掛けてから垂直入射へ変換する**
（`absorption_table()`）。順番が逆だと意味が変わる。
実測（研修室・全レイヤ 0.8）：平均 α 1 kHz が 0.539 → 0.443、
理論 T30 が 0.198 → 0.261 s（安全側に動く）。

#### ★更新は「その場で直す」

既存の xlsx は**書き直さない**。利用者が作った体裁（列幅・書式・番号の式・
列の並び・増やした列）を壊さないため、**見出しの文字から列を探し**、
面数・面積・参考の式だけを書き換える（`_update_book` / `_columns_of`）。
足りない列（安全率・材料名）は**右端に足す**（間に挿し込むと既にある
`VLOOKUP($D2,…)` の参照がずれるため）。

#### ★入力は「材料番号」と「安全率」だけ（ユーザー指摘）

> 材料名だと、全角・半角とかのミスタイプを誘発する可能性が高いので。

そこで**番号を入力する列**を作り、材料名は隣に**参考として出す**。
参考の列は `VLOOKUP` で吸音率シートを引くので、Excel 上で番号を入れれば
その場で名前が出る。**読み込むときは番号の列しか見ない**（名前の列は人向け）。
番号は数値として書く（Excel で並べ替えられるように）。

#### ★吸音率は条件シートに載せない（ユーザー指摘）

> 条件表に吸音率のデータを載せるのは、ミスリードを誘発するのでやめて欲しいです。

条件シートには α を一切書かない。吸音率は「吸音率」シートにだけ置く。
そのぶん**PJ 固有の吸音データをこの 1 ファイルに閉じ込められる**
（`library_from_book()` が材料一覧として読むので `absorption.csv` は無くてもよい。
`run_project._library_for()` が「吸音率シート → 吸音率 CSV」の順で探す）。

#### 作る・選ぶ

| したいこと | やり方 |
|---|---|
| DXF のレイヤから新しく作る | 設定画面の**「条件表を作成」**（`create()`）。`python condition_table.py <フォルダ>` でも |
| 用意したものを使う | 設定画面の「条件表（xlsx）」で**ファイルを選ぶ**（`参照…`） |
| 条件を選ぶ | 設定画面の**「条件（シート）」**のコンボボックス |
| 条件を増やす | Excel でシートを複製して名前を変え、番号を書き換える |
| 全条件を回す | 「全条件を一括 ▶▶」／`--conditions`（**シートごとに 1 条件**） |

#### 更新の約束（★ここを崩さないこと）

**利用者が書いた「材料番号」は絶対に上書きしない。**
`update()` は既存の全シートを読み、番号はそのまま残して面数・面積と参考の式だけ
書き直す。新しいレイヤの行だけを足し、モデルから消えたレイヤは区分を
`（モデルに無し）` にして残す（消すと前の設定が分からなくなる）。

#### 区分（いまある値と、将来ありうる値）

| 区分 | 意味 | 入力か記録か |
|---|---|---|
| `レイヤ` | DXF のレイヤに材料番号を割り当てる | **入力** |
| `面ごとの指定` | `face_editor` で面を選んで貼った材料 | 記録 |
| `（モデルに無し）` | 前は書いてあったが、いまのモデルに無いレイヤ | 記録 |

**将来ありうるもの**（まだ実装していない。2026-08-21 にユーザーから
「区分は必要ですか？」と聞かれたので、可能性を書き残しておく）:

- `追加吸音` … 客席・人・什器のように**形を作らず面積だけで足したい吸音**
  （等価吸音面積 A [m²] を直接足す。実務でよく使う）
- `グループ` … 複数のレイヤをまとめて 1 つの材料にする行

★区分が要らないと決めたら、`面ごとの指定` と `（モデルに無し）` を別シートに
移してから列を落とす。**記録の置き場が無くなるので急に消さない。**

#### 昔の形も読める

`材料条件表.csv`（2026-08-21 の日中の形。3 列目が材料名）もそのまま読める。
ただし**既定名の CSV は xlsx があれば条件に数えない**
（移行後に同じ内容が 2 条件として並ぶのを防ぐ）。
名前を付けた CSV（`別案.csv` など）は条件として扱う。

### workbook.py ― 結果一式の Excel（依頼 2026-08-21）

> CSV の方が、逐次読み込む際に都合が良いですか？ .xlsx でまとめた方が、
> グラフを作ってもらうのも、テンプレに当てこむだけ、というのが出来るかなと。
> ソフトで読み込む用のデータと、体裁を整える用のファイルを分けても良いかも。

**分ける**のが答え。CSV はそのまま残し、その上に Excel を 1 枚かぶせる。

| | CSV | Excel |
|---|---|---|
| 誰が読むか | プログラム（`--redraw`・まとめ表・他ソフト） | 人（報告書・打ち合わせ） |
| 良いところ | 1 行ずつ読める・差分が見える・文字コードだけで開ける | 1 枚で全部見える・**グラフ付き**・体裁を作れる |

`結果/<室>_結果一式.xlsx` に 7 シート（概要・残響時間・明瞭度・音圧レベル・STI・
吸音率と理論値・材料条件表）。まとめ表の CSV を読んで並べ直すだけなので
**計算はやり直さない**。グラフは openpyxl のネイティブ図（Excel で編集できる）。

`--template 雛形.xlsx` を渡すと**雛形を開いて同じ名前のシートに値だけ書く**。
雛形側のグラフ・書式・ロゴは残る。**書き出す位置は必ず A1 から**にしてあるので
雛形のセル参照がずれない。雛形に無いシートは作るが、**雛形のシートは消さない**。

openpyxl が無い環境では Excel だけ作れずに済む（CSV は書ける）ようにしてある。

### setup_window.py ― 条件入力（tkinter）

**tkinter（標準ライブラリ）で書いている。依存を増やさないため。**
3D 表示が要るところは PyVista 側の別ウィンドウに任せる。

`ask()` が `(Project, 'run'|'normals'|None)` を返すだけなので、
将来 PySide6 などに載せ替えるときもこの契約を守れば `app.py` は変えずに済む。

**「対象室・条件名」欄がそのまま結果ファイル名の頭になる**
（2026-08-21。空欄ならフォルダ名）。容積の「モデルから見積もる」は
容積と**総表面積・平均自由行程 4V/S** も出す（参考値。ユーザー要望）。

### reverberation.clarity_measures() ― 明瞭度系の指標（G-4）

残響指標が「どれだけ長く響くか」を見るのに対し、こちらは
「**初期の音が後から来る音に対してどれだけ強いか**」を見る。会議室・教室で効く。

| 指標 | 定義 | 用途 |
|---|---|---|
| C50 | `10 log10( ∫₀^50ms p² / ∫_50ms^∞ p² )` [dB] | 音声 |
| C80 | 同上、境目 80 ms [dB] | 音楽 |
| D50 | `∫₀^50ms p² / ∫₀^∞ p²` [0〜1] | 音声（C50 の言い換え） |
| Ts | `∫ t·p² / ∫ p²` [s] | 重心時刻。小さいほど明瞭 |

**時刻の起点は直接音の到来時刻**にする。音源から受音点までの伝搬時間ぶん、
インパルス応答の先頭には無音があり、そこを含めると 50 ms の窓がずれるため。
各バンドでエネルギーが最大になる時刻を直接音とみなしている。

### read_dxffile.check_model() ― 作図ミスの自動チェック（B-10）

計算に入る前に「そもそもモデルとして成立しているか」を機械的に確かめる。
黙って変な結果を出すより、先に指摘して直してもらうほうが早いため。
`procedure.process()` が計算前に必ず呼ぶ。

面が読めているか／音源・受音点の有無と室内にあるか／吸音率が引けないレイヤ／
重なった面／開いた辺／巻き順／ねじれ／極端に小さい面／読み飛ばした種別を見て、
`{'level': 'error'|'warning'|'info', 'message': str}` のリストを返す。
