# geosim — プログラム構成まとめ

幾何音響シミュレーション（音線法 + 虚音源バックトレース）の Fortran → Python 移植プロジェクト。
元コードは `fortran/` の 3 ファイル（`backtrace.f90` / `make_ipls_freq_monaural_fortran.f90` / `ipls2rt_fortran.f90`）。
本書は 2026-07-17 時点の `geosim/` パッケージの内容整理（コードには手を加えていません）。

## 処理パイプライン全体像

```
入力: 音源位置・受音点・DXF形状・受音球半径・反射回数・音線数・吸音率リスト

procedure.process()
 ├─ 1. read_dxffile.read()            … DXF → Meshオブジェクトのリスト
 ├─ 2. sound_ray.soundray_generator() … 球面上に均等な音線ベクトル群を生成
 ├─ 3. loop_reflectionmesh.loop()     … 音線追跡。受音した経路の反射面ID履歴を記録
 │                                       （元コード backtrace.f90 524行〜）
 ├─ 4. loop_deleteredundancy.delete() … 反射履歴の重複経路を削除（元コード 721行〜）
 ├─ 5. loop_noredundancy.loop()       … 虚音源法バックトレース。到来時間・エネルギーの
 │                                       パルス列を出力（元コード 876行〜）※未接続
 └─ 6. impulse.impulse_responce()     … パルス列 → バンド合成でインパルス応答CSV
                                         （元コード make_ipls_freq_monaural_fortran.f90）※未接続

未着手: 残響時間算出（ipls2rt_fortran.f90 相当）
```

## ファイル一覧

| ファイル | 役割 | 対応する元コード | 状態 |
|---|---|---|---|
| `procedure.py` | 全体フローのオーケストレーション | backtrace.f90 全体 | 骨格のみ（後半はコメントアウト） |
| `mesh.py` | 三角形メッシュのデータクラス | vwpt/surf/vwnm/absper 配列群 | 実装済（2026-08-12 に vertexes の形状を修正） |
| `mesh_method.py` | 音線と面の交差判定の幾何計算 | 376〜447 行ほか | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `sound_ray.py` | 音線の生成・正規化・反射・エネルギー減衰 | 318〜326, 493〜503, 1091〜1110 行 | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `receiver_sphere.py` | 受音球の通過判定 | 649〜663 行 | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| `read_dxffile.py` | DXF ファイルの読み込み | 132〜283 行 | スタブ（座標未取得） |
| `loop_reflectionmesh.py` | 音線追跡ループ本体 | 524〜717 行 | 実装中（構造の疑問点をコメントで整理中） |
| `loop_deleteredundancy.py` | 重複経路の削除 | 721〜841 行 | デモ実装済（set 方式） |
| `loop_noredundancy.py` | バックトレースループ | 876〜1134 行 | 下書き段階 |
| `impulse.py` | インパルス応答の合成・CSV 出力 | make_ipls_freq_monaural_fortran.f90 | 実装中（転記ミス多数） |
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

全体フローの記述。エントリポイント。

```python
process(soundsource_point, reciever_point, dxf_filename, sphere_radius, nref, soundray_number) -> None
```

| 引数 | 型（想定） | 意味 |
|---|---|---|
| `soundsource_point` | ndarray (3,) | 音源座標 [m] |
| `reciever_point` | ndarray (3,) | 受音点座標 [m] |
| `dxf_filename` | str | 形状 DXF ファイルパス |
| `sphere_radius` | float | 受音球半径 [m] |
| `nref` | int | 最大反射回数 |
| `soundray_number` | int | 音線数 |

- 冒頭でオクターブバンド中心周波数 `[63, 125, 250, 500, 1k, 2k, 4k, 8k]`（8 バンド）を定義。
  ※元コードの吸音率は 6 バンド（absorption.csv の a1〜a6）なので不整合あり。要調整。
- ステップ 5（バックトレース）以降は呼び出しがコメントアウト。吸音率リストの読み込み処理も未定。

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

```python
read(file_name: str) -> list[Mesh]
```
DXF から Mesh リストを作るのが目的。現状は `POLYLINE` 行の検出と直後数行の print までで、
座標・法線・吸音材の取得は未実装（戻り値は常に空リスト）。
元コード（132〜283 行）で行っている処理: グループコード 70 のフラグで頂点(192)/面(128)を判別、
レイヤ名から吸音材 ID 取得、mm→m 変換、縮退面（法線ゼロ）の除去、外積による法線計算と正規化。
`if __name__ == "__main__":` で単体テスト実行可（テスト用パスがハードコード）。

### loop_reflectionmesh.py

```python
loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius)
    -> list[list[int]]   # 受音した経路ごとの反射面ID履歴（先頭要素は -1）
```
音線追跡の本体（元コード 524〜717 行）。3 重ループ構造:

```
for 音線 i:                       # 元コード 527行 do i = 1, nray
    履歴初期化（先頭に -1）
    for 反射回数 k:               # 元コード 545行 do k = 0, nref
        音線を正規化
        for 面 j:                 # 元コード 576行 do j = 1, sfcount
            前方かつ面に向かう場合に衝突判定 → 最寄り面 mesh_nearestid を更新
        受音球判定 inside_sphere()
        if 衝突あり: 交点へ基点を移動し、反射ベクトルを生成して継続
        else: break（この音線は終了）
```

現状の未確定点（コード内の 12/04 メモ・疑問 1〜4 で整理中）:
- 元コードでは「衝突のたびに履歴へ壁番号追加」+「受音のたびに履歴のスナップショットを 2 次元リストへ保存」の 2 本立て。
  現状は `if inside` の中でのみ履歴追加する形になっており、置き場所を検討中。
- 履歴の 2 次元リストへの append が反射ループの毎周・同一リスト参照で行われている（コピーが必要）。
- 面ループ内で `soundray_comesfrom` を上書きしている箇所あり（元コードでは最寄り面確定後に更新）。
- 冒頭の `from numpy.ma.core import count` は不要な自動インポート。

### loop_deleteredundancy.py

```python
delete(reflectionhistory_redundancy: list[list[int]]) -> list[list[int]]
```
反射履歴の重複経路削除（元コード 721〜841 行）。
元コードの「反射回数でソート → 切替点探索 → 同一反射回数内で総当たり比較」を、
Python では `list → tuple → set → list` の一括重複除去で置き換え（結果は等価、順序は不定）。
`if __name__ == '__main__':` にデモ入力あり。

### loop_noredundancy.py

```python
loop(soundsource_point, reciever_point, absorption_list, reflectionmeshid_history, mesh) -> None
```
虚音源法バックトレース（元コード 876〜1134 行）。**下書き段階**。意図する処理:

```
for 非重複経路 i:
    反射履歴に沿って虚音源列 isrc(k) を鏡映で順次生成    # 元コード 900〜926行
    受音点から最終虚音源へ向けて逆向きに追跡:            # 元コード 948行 do k = ktmp, 0, -1
        各段で最寄り面を探索し、履歴の面と一致するか検証
        一致すれば energy_decay() でバンド別エネルギーを減衰、基点と音線を更新
        k=0 で遮蔽チェック → 通れば「反射回数, 到来時間, 方向ベクトル, バンド別エネルギー」を出力
```

出力形式（元コード 1080 行）: `ktmp, rtime, -vtgt(1:3), enertmp(1:6)` の 11 列テキスト。
現状は履歴を ndarray 前提（`.shape`）で扱う一方 delete() はリストを返す、
`mesh[mesh_id[i, k]]` のスカラー添字、`energy_decay` の引数不足、衝突判定ブロックの位置など、骨格から要整理。
ファイル出力は print のダミー。戻り値なし（ファイル保存が目的）。

### impulse.py

パルス列 → インパルス応答の合成（元コード make_ipls_freq_monaural_fortran.f90 全体）。

```python
impulse_responce(filename, sound_velocity, reflection_timing, soundsourceenergy_list,
                 frequency_number, count, nn, mf, fmax, dt) -> None
```
全体フロー関数。以下を順に呼ぶ:

| 関数 | 対応元コード | 処理 | 主な入出力 |
|---|---|---|---|
| `airdamping_coefficient()` | 41 行 | 空気吸収係数 mair（1/3 オクターブ 32 バンド、20℃湿度40%の近似） | → ndarray (32,) |
| `power_airdamping(sound_velocity, reflection_timing, mair, soundsourceenergy_list)` | 92 行 | 各パルスのエネルギーに exp(−mair·c·t) を乗算 | → ndarray (32, パルス数) |
| `transfer_function(reflection_timing, frequency_number, count, p, sound_velocity)` | 132 行 | パルス列を周波数領域の伝達関数に変換（√E · e^(−j2πf·t) / (t·c) の重ね合わせ） | → complex ndarray (32, nfreq) |
| `time_responce(nn, mf, fmax)` | 147 行 | バンドごとの時間領域バンドパスフィルタ生成 | → ndarray (32, nn) |
| `filter_bandpass(n, min, max)` | fir1_bandpass | Hamming 窓付き FIR バンドパス（sinc の差 × 窓） | → ndarray (n+1,) |
| `fft_timeresponce(ht, nn)` | 158 行 | フィルタの FFT | → ndarray (32, nn) |
| `fft_negativerange(transfer, nn)` | 169 行 | 伝達関数の負周波数側を複素共役で拡張 | → ndarray (32, nn) |
| `convolution_hfhtf(hf, htf, nn)` | 184 行 | 周波数領域の乗算（= 時間領域の畳み込み） | → ndarray (32, nn) |
| `inversefft_responce(hfcc, nn)` | 196〜216 行 | 逆 FFT 後に 32 バンドを合算して IR に | → ndarray (nn,) |
| `write_impulseresponce(filename, ir, dt, nn)` | 220, 230 行 | 時間ベクトルを付けて CSV 出力 | ファイル出力 |
| `fft_filter(x, n, sign)` | subroutine fft | 自作 Cooley-Tukey FFT（sign=1 順変換 / −1 逆変換+正規化） | → ndarray |

既知の課題（前回レビューで詳細指摘済み）:
- 転記ミス: Hamming 窓の `+`→`*`、伝達関数の `+`→`*`、空気吸収の括弧位置と配列代入、負周波数の添字、逆 FFT 呼び出しの欠落 ほか
- `df`（周波数離散化幅）が `transfer_function` に未受け渡し（コード内メモあり）
- Fortran 側にある 6 バンド吸音率 → 32 バンドへの展開表が未移植
- 自作 FFT は `np.fft`、バンドパスは `scipy.signal.firwin` への置き換えを検討中（コード内メモあり）

---

## Fortran 側との対応・未移植の機能

| 元コードの機能 | Python 側の状況 |
|---|---|
| 吸音率 CSV（absorption.csv）読み込み | 未実装（procedure.py にメモのみ） |
| DXF 形状読み込み・法線計算 | スタブ（read_dxffile.py） |
| 有効経路カウント用の事前ループ（backtrace.f90 330〜515 行） | 不要（Python はリストの動的伸長で代替） |
| 音線追跡ループ | 実装中（loop_reflectionmesh.py） |
| クイックソート + 重複削除 | set 方式で代替済み（loop_deleteredundancy.py） |
| バックトレース + 受音リスト出力 | 下書き（loop_noredundancy.py） |
| インパルス応答合成（make_ipls_...f90） | 実装中（impulse.py） |
| 残響時間算出（ipls2rt_fortran.f90） | 未着手 |
| OpenMP 並列化 | 未検討（NumPy ベクトル化 / numba が候補） |
