# geosim — プログラム構成まとめ

幾何音響シミュレーション（音線法 + 虚音源バックトレース）の Fortran → Python 移植プロジェクト。
元コードは `fortran/` の 3 ファイル（`backtrace.f90` / `make_ipls_freq_monaural_fortran.f90` / `ipls2rt_fortran.f90`）。
本書は 2026-07-17 時点の `geosim/` パッケージの内容整理（コードには手を加えていません）。

## 実行環境

**Python 3.10.11**（チーム方針）。依存は `requirements.txt`（numpy / scipy / matplotlib / pyvista）。
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
 ├─ 6. impulse.impulse_responce()     … パルス列 → バンド合成でインパルス応答CSV
 │                                       （元コード make_ipls_freq_monaural_fortran.f90）
 └─ 7. reverberation.reverberation_time() … 減衰曲線と残響時間 T30
                                         （元コード ipls2rt_fortran.f90）

2026-08-14 に 5〜7 を実装して全段つながった。
```

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
| `procedure.py` | 全体フローのオーケストレーション | backtrace.f90 全体 | 骨格のみ（後半はコメントアウト） |
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
| `atmosphere.py` | 温度・湿度・気圧 → 音速・空気吸収 | c0 定数 / mair 近似式 | **実装済**（2026-08-14 新設。ISO 9613-1） |
| `absorption.py` | 吸音材ライブラリと吸音率の種類の変換 | absper 配列 | **実装済**（2026-08-14 新設。Paris の式・GUI 対応） |
| `loop_deleteredundancy.py` | 重複経路の削除 | 721〜841 行 | 実装済（2026-08-12 に整理・動作確認済） |
| `loop_noredundancy.py` | バックトレースループ（虚音源法） | 876〜1134 行 | **実装済**（2026-08-14 全面書き直し・解析解と一致） |
| `reverberation.py` | 残響時間・減衰曲線 | ipls2rt_fortran.f90 | **実装済**（2026-08-14 新設・既知の減衰と一致） |
| `impulse.py` | インパルス応答の合成・CSV 出力 | make_ipls_freq_monaural_fortran.f90 | **実装済**（2026-08-14 全面書き直し・scipy 化） |
| `project.py` | プロジェクト（条件・法線指定・結果）の保存と読み込み | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `plots.py` | 結果を PNG にする（G-3 / G-4） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `setup_window.py` | 計算条件の入力ウィンドウ（tkinter） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `normal_editor.py` | 法線の確認・修正ウィンドウ（G-8） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `run_project.py` | プロジェクトの条件で計算を回す薄い層 | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
| `app.py` | 入口（条件入力 → 法線確認 → 計算 → 可視化） | （移植元なし・新規） | **実装済**（2026-08-15 新設） |
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
process(soundsource_point, reciever_point, dxf_filename, sphere_radius, nref, soundray_number,
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
| `reciever_point` | (3,) \| None | 受音点座標 [m]。None なら DXF の rec レイヤから |
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
noramlized_soundray(sound_ray) -> ndarray (3,)          # L2ノルムで正規化（関数名は normalized の綴りミス）
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
loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius,
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
loop(soundsource_point, reciever_point, reflectionmeshid_history, mesh,
     sound_velocity=340.0, band_number=None, filename=None, verbose=True) -> PulseList
backtrace_path(soundsource_point, reciever_point, wall_ids, mesh,
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
impulse_responce(filename, pulses, sound_velocity=340.0, sampling_frequency=44100.0,
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
| `write_impulseresponce(...)` | 232〜235 行 | 時間ベクトルを付けて CSV 出力 |

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
| Sabine | `S·ᾱ + 4mV` | 吸音率が小さいとき（〜0.2） |
| Eyring | `-S·ln(1-ᾱ) + 4mV` | 吸音率が大きいとき |
| Millington | `-Σ Sᵢ·ln(1-αᵢ) + 4mV` | 面ごとに吸音率が大きく違うとき |

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
| 音線追跡ループ | 実装中（loop_reflectionmesh.py） |
| クイックソート + 重複削除 | set 方式で代替済み（loop_deleteredundancy.py） |
| バックトレース + 受音リスト出力 | 下書き（loop_noredundancy.py） |
| インパルス応答合成（make_ipls_...f90） | 実装中（impulse.py） |
| 残響時間算出（ipls2rt_fortran.f90） | 未着手 |
| OpenMP 並列化 | 未検討（NumPy ベクトル化 / numba が候補） |

---

## GUI とプロジェクト（2026-08-15 新設）

**1 つのウィンドウで完結する GUI が最終形**だが、どんな情報が要るかがまだ固まっていないので、
いまは「必要なウィンドウをそのつど開く」形にしてある。
各ウィンドウは独立して呼べるので、統合するときは `app.py` を差し替えれば済む。

```
app.py  ──┬─→ setup_window.py   条件入力（tkinter。依存を増やさないため）
          ├─→ normal_editor.py  法線の確認・修正（PyVista）
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
  normals.json          法線の反転指定（normal_editor.py が書く）
  結果/  pulses.csv  ir.csv  rt.csv  rt_statistical.csv  clarity.csv
         decay.csv  surface.csv  raylog.npz
  図/    impulse_response.png  decay.png  reverberation.png
         clarity.png  absorption.png  pulses.png
  rec2/  受音点が複数ある場合、2 点目以降は同じ構造で枝分かれ（B-9）
```

- DXF と吸音率 CSV は**プロジェクトフォルダからの相対パスで持つ**
  （フォルダごと別の端末へ移してもそのまま開ける）。外にあるものは絶対パスのまま
- `DEFAULTS` に無いキーは保存されない。**新しい計算条件を足したらここにも足すこと**
- `normals.json` は**書いたときの DXF と面数を控えてある**。食い違ったら使わずに知らせる
  （黙って間違った面を反転させないため）

### normal_editor.py ― 法線の確認・修正（G-8）

読み込み時の自動補正（`orient_normals='auto'`）に任せきりにせず、
**人が見て確認し、必要なら反転できる**ようにするためのもの。

| 色 | 意味 |
|---|---|
| 緑 | 法線が室内（空気側）を向いている |
| 赤 | 法線が室外を向いている。**反転が要る** |
| 灰 | 判定できない（開いた形状など）。CAD の指定を尊重する |

操作は「枠で囲んだ面を反転 / 数字キーでレイヤごと反転 / `a` 自動 / `c` CAD に戻す /
`x` 全反転 / `s` 保存して閉じる」。

状態は `flipped`（**CAD の巻き順から反転する面のインデックス集合**）ただ 1 つで、
`read_model(flip_faces=...)` にそのまま渡せる形にしてある。
自動判定の結果もここに畳み込んで保存するので、**保存したものを読めば再現できる**。

適用の順序は **自動判定 → 手動指定** で、手動が最後の決定になる。

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

**図を読むときの注意**
- 位相がばらばらなら和はランダムウォークで √N 程度にしかならない。
  図に √N の線を引いてあり、**超えている山だけが本当に強め合っている周波数**
- 0 Hz は全経路が同位相になるので左端に巨大な山が立つ。
  **いちばん低い固有周波数から下は描かない**
- 吸音率はオクターブバンドしか無いので、各周波数点には最も近いバンドの値を当てる
  （`_band_of`）。6 バンドだと 88 Hz 以下は全部 125 Hz の値になる
- パルス列は受音球が拾えた経路の**標本**。本数の絶対値は音線数と受音球半径で変わる

### setup_window.py ― 条件入力（tkinter）

**tkinter（標準ライブラリ）で書いている。依存を増やさないため。**
3D 表示が要るところは PyVista 側の別ウィンドウに任せる。

`ask()` が `(Project, 'run'|'normals'|None)` を返すだけなので、
将来 PySide6 などに載せ替えるときもこの契約を守れば `app.py` は変えずに済む。

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
