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
| `read_dxffile.py` | DXF ファイルの読み込み | 132〜283 行 | **実装済**（2026-08-12。実 DXF で動作確認済） |
| `loop_reflectionmesh.py` | 音線追跡ループ本体 | 524〜717 行 | **実装済**（2026-08-12 に traceff の扱いを確定し全面修正・動作確認済） |
| `ray_recorder.py` | 可視化用の音線軌跡レコーダ | （移植元なし・新規） | 実装済（2026-08-12 新設） |
| `loop_deleteredundancy.py` | 重複経路の削除 | 721〜841 行 | 実装済（2026-08-12 に整理・動作確認済） |
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
- **ねじれた四角形**: `quad_warp()` で平面からのずれを測り、`WARP_TOLERANCE` を超えたら警告。
  どの対角線で切るかで形が変わるため、CAD 側で直してもらう
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

⚠ かつてあった `'auto'` / `'toward'` / `'away'` / `'inward'` / `'outward'` は**廃止**（`ValueError`）。
法線を音源方向や重心方向へ面ごとに向ける方式は、**凸凹の壁や宙に浮いた家具で破綻する**ため。

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
未登録のレイヤは既定値（0.1）を使い、レイヤ名を列挙して警告する。
`read_absorption_csv()` は 2 形式を自動判別する。

| 形式 | 列 | 備考 |
|---|---|---|
| (A) 元コード付属 `absorption.csv` | `ID, 材料名, a1〜a6` | **ID と材料名の両方をキーに登録**するのでレイヤ名がどちらでも引ける。CP932 なので UTF-8 → CP932 → latin-1 の順にデコードを試す |
| (B) `data/absorption_sample.csv` | `材料名, a1〜a6` | 2 列目が数値かどうかで (A) と判別 |

`if __name__ == "__main__":` でリポジトリ直下の `test.dxf` を読むテストが走る。

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
| 吸音率 CSV（absorption.csv）読み込み | 実装済（`read_absorption_csv()`。サンプル `data/absorption_sample.csv`） |
| DXF 形状読み込み・法線計算 | 実装済（`read_dxffile.py`。単位換算・法線向き自動判定・音源受音点の読み込みは元コードにない追加機能） |
| 有効経路カウント用の事前ループ（backtrace.f90 330〜515 行） | 不要（Python はリストの動的伸長で代替） |
| 音線追跡ループ | 実装中（loop_reflectionmesh.py） |
| クイックソート + 重複削除 | set 方式で代替済み（loop_deleteredundancy.py） |
| バックトレース + 受音リスト出力 | 下書き（loop_noredundancy.py） |
| インパルス応答合成（make_ipls_...f90） | 実装中（impulse.py） |
| 残響時間算出（ipls2rt_fortran.f90） | 未着手 |
| OpenMP 並列化 | 未検討（NumPy ベクトル化 / numba が候補） |
