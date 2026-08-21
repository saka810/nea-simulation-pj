# nea-simulation-pj

NEA（日本環境アメニティ株式会社）のシミュレーション PJ。
幾何音響シミュレーション（音線法 ＋ 虚音源バックトレース）の **Fortran → Python 移植**。

- 各モジュールの役割・Fortran との対応・実装状況 … [PROGRAM_STRUCTURE.md](PROGRAM_STRUCTURE.md)
- 作業一覧 … [TODO.md](TODO.md)
- 数式とフローの解説 … [docs/技術説明書.md](docs/技術説明書.md)
- **要判断**：鏡像の式の `abs()` 問題 … [docs/議論_鏡像の式のabs問題.md](docs/議論_鏡像の式のabs問題.md)
- 出力・可視化の方針 … [docs/出力・可視化方針.md](docs/出力・可視化方針.md)
- CAD 側の作図ルール … [docs/DXFデータの作り方.md](docs/DXFデータの作り方.md)

## フォルダ構成

    geosim/                 Python 移植版のパッケージ（本体。27 ファイル）
    tests/test_geosim.py    数値検証。pytest 不要、素の Python で走る
    docs/                   数式・フロー・CAD の作図ルールの解説
    data/                   吸音率テーブルのサンプル
    Claude履歴/             セッション履歴（端末をまたいで作業を引き継ぐため）
    参考文献/               論文・書籍（PDF 本体はローカルのみ）
    fortran/                移植元の Fortran 3 ファイル（**Git 管理外**）

`geosim/` の中身は**パイプラインの順**に読むのが早い。

| 段階 | モジュール | 元コード |
|---|---|---|
| 土台 | `mesh.py`（面 1 枚）/ `mesh_method.py`（交差判定）/ `sound_ray.py`（音線）/ `receiver_sphere.py`（受音判定） | `backtrace.f90` 各所 |
| ① DXF 読込 | `read_dxffile.py` | 132〜283 行 |
| ②③ 音線生成・追跡 | `loop_reflectionmesh.py`（＋`ray_recorder.py` が可視化用の軌跡を別チャンネルで記録） | 318〜326 / 524〜717 行 |
| ④ 重複経路の削除 | `loop_deleteredundancy.py` | 721〜841 行 |
| ⑤ バックトレース | `loop_noredundancy.py` | 876〜1134 行 |
| ⑥ インパルス応答 | `impulse.py` | `make_ipls_freq_monaural_fortran.f90` |
| ⑦ 残響時間・音響指標 | `reverberation.py` | `ipls2rt_fortran.f90` |
| 材料・大気 | `absorption.py`（吸音率の種類の変換）/ `atmosphere.py`（音速・空気吸収） | — |
| 通し実行 | `procedure.py` | — |
| 出力 | `plots.py`（図）/ `table.py`（表の並べ方の共通ルール）/ `project.py`（保存・読込） | — |
| 絞り込み | `ray_filter.py`（注目したい音線を残す。近く・方向・1 本） | — |
| 画面 | `app.py`（入口）/ `setup_window.py`（条件入力）/ `progress_window.py`（進捗）/ `face_editor.py`（面の確認：法線・吸音材）/ `view_rays.py`（音線・音粒子）/ `view_directions.py`（音線の飛び方）/ `view_model_gui.py`・`view_model.py`（モデルビューア） | — |

各モジュールの詳細は [PROGRAM_STRUCTURE.md](PROGRAM_STRUCTURE.md)。

## 環境構築

**Python 3.10.11**（チーム方針。[.python-version](.python-version) を参照）。

端末によっては `python` が 3.10 以外を指すので、**インタプリタを明示して** venv を作る。

```powershell
# Windows PowerShell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

```bash
# Git Bash など
py -3.10 -m venv .venv
source .venv/Scripts/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

確認:

```powershell
python --version            # -> Python 3.10.11
python -c "import numpy, scipy, matplotlib, pyvista, openpyxl; print('ok')"
```

端末間で完全に同じ組み合わせにしたい場合は `requirements-lock.txt` を使う
（間接依存まで固定した、動作確認済みの一式）。

`.venv/` は Git 管理外。端末ごとに作り直す。

### 依存ライブラリ

| ライブラリ | 用途 |
|---|---|
| numpy | 幾何計算・配列演算（パッケージ全体の土台） |
| scipy | 帯域フィルタ・信号処理（E-6） |
| matplotlib | 2D グラフ（インパルス応答・音響指標・モード分布） |
| pyvista（+ vtk） | 3D 表示（`geosim/view_model_gui.py`） |
| openpyxl | 結果一式の Excel 出力（`geosim/workbook.py`）。CSV 出力には不要 |

## 実行

`geosim/` は `from mesh import Mesh` のような**ベア名インポート**なので、
`geosim` を作業ディレクトリにして直接実行する（C-2 で見直し予定）。

### GUI から使う（おすすめ）

```powershell
cd geosim
python app.py                            条件入力から始める
python app.py "C:\...\プロジェクト"        そのプロジェクトを開いて始める
python app.py "C:\...\プロジェクト" --run  入力ウィンドウを出さずにすぐ計算
```

流れは **条件入力 → 面の確認 → 計算 → 可視化**。

1. **条件入力**（`setup_window.py`）… モデル DXF・吸音率 CSV・保存先フォルダ・
   音線数・受音球の半径・温度湿度などを入力する。
   既存フォルダを選べば前回の条件を読み込む
2. **面の確認**（`face_editor.py`、「面を確認…（法線・吸音材）」ボタン）…
   **法線**は 緑=内向き（OK）／赤=外向き（要反転）／灰=判定できない で色分け。
   **吸音材**は `m` で切り替えて材料ごとの色で表示する。
   `r` で枠選択して面を選び（**同一平面の三角形はまとめて 1 枚として選ばれる**）、
   `i` で法線を反転、左パネルの材料をクリックで吸音材を貼る。`s` で保存。
   レイヤで吸音材を分けられないモデル（1 つの 3DSOLID で出来ているなど）は、
   ここで面ごとに貼る
3. **計算**（「計算する ▶」）… 結果 CSV と図 PNG がプロジェクトフォルダに保存される
4. **可視化**（`view_rays.py`）… 音線と音粒子を `Tab` で切り替えて見る

プロジェクトフォルダの中身:

```
project.json          条件（次に開いたときそのまま復元される）
normals.json          法線の反転指定
materials.json        面ごとの吸音材の割り当て（レイヤで分けられないモデル用）
条件表.xlsx           レイヤー名 → 材料番号（★CAD を触らずに材料を変える入力）
                      シート「吸音率」＝材料一覧／その他のシート＝条件 1 つ
結果/                 ← 受音点に依らないものは直下
  研修室_条件A_結果一式.xlsx          ★体裁用（全部入り・グラフ付き）
  研修室_条件A_まとめ_残響時間.csv    全受音点 ＋ 平均 ＋ 理論値
  研修室_条件A_まとめ_明瞭度.csv      全受音点 ＋ 平均
  研修室_条件A_まとめ_音圧レベル.csv  帯域別 Lp と自由音場（逆二乗）との差
  研修室_条件A_まとめ_STI.csv         STI と帯域別 MTI
  研修室_条件A_吸音率と理論値.csv     材料別の吸音率 → 平均吸音率 → 残響時間理論値
  研修室_条件A_raylog.npz             可視化用の音線軌跡
  rec1/ rec2/ …       受音点ごとに pulses / ir / rt / decay / clarity / spl / sti の CSV
図/
  rec1/ rec2/ …       受音点ごとの PNG
  画面/               画面から手で撮った画像・動画（計算し直しても消えない）
```

**CSV と Excel は役割が違う。** CSV はプログラムが読む（描き直し・まとめ表・他ソフト）、
Excel は人が読む（報告書・打ち合わせ）。Excel はまとめ表の CSV を並べ直したもので、
グラフが付いている。会社の雛形に流し込むこともできる:

```powershell
.\.venv\Scripts\python geosim\workbook.py "<プロジェクト>" --template 雛形.xlsx
```

**吸音材を変えるときは CAD を触らない。** 条件入力の**「条件表を作成」**を押すと、
DXF のレイヤを並べた `条件表.xlsx` ができる。

```
シート「吸音率」   番号・材料名・吸音率（★PJ 固有の吸音データをここに置ける）
シート「現状」     区分 / レイヤー名 / 材料番号 / 材料名（参考）/ 面数 / 面積
```

**書き換えるのは「材料番号」の列だけ**（材料名は VLOOKUP で自動表示。
名前を打つと全角・半角のミスタイプが起きるため）。**条件シートに吸音率は載せない**
（ミスリードを避けるため。吸音率は「吸音率」シートにだけ置く）。

**複数条件は一括で回せる。** シートを複製して名前を変え（シート名が条件名）、
条件入力の「全条件を一括 ▶▶」を押す（CLI なら
`python geosim\run_project.py "<プロジェクト>" --conditions`）。

- **結果ファイル名は「DXF のファイル名 ＋ 条件表のシート名」**になるので混ざらない
- **2 件目以降は速い。**経路（反射面の並びと入射角）は吸音に依らないので、
  1 件目で保存した `結果/recN/<室>_経路.npz` を使い回してエネルギーだけ計算し直す
- 最後に**条件を横に並べた比較表**（`<室>_まとめ_条件比較.csv`）と、
  Excel の「条件比較」シートができる

**ファイル名の頭には「対象室＋条件名」（プロジェクト名）が付く。**
報告書やメールでフォルダの外へ出したときに、どの室・どの条件のものか
分かるようにするため。図（PNG）にも同じ頭が付く。

**インパルス応答の図は最大値で正規化してある。** 絶対振幅の校正が未了なので、
そのままの値には意味が無いため（TODO E-11）。

### コマンドラインから使う

```powershell
cd geosim

# プロジェクトの条件で計算する（GUI を使わずに回す）
python run_project.py "C:\...\プロジェクト"

# 面の確認だけ開く（法線・吸音材）
python face_editor.py "C:\...\プロジェクト"

# モデルの幾何だけを確認する（容積・総表面積・レイヤ別面積。音源や受音点は不要）
python read_dxffile.py ..\test.dxf

# モデルビューア（ネイティブウィンドウ / PyVista）
python view_model_gui.py ..\test.dxf --absorption ..\absorption.csv

# モデルビューア（HTML + WebGL を書き出してブラウザで開く。依存ライブラリ不要）
python view_model.py ..\test.dxf --absorption ..\absorption.csv

# 音線と音粒子（既定は両方。Tab で切り替え。閉じずに見比べられる）
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --received-only --max-reflection 3 --color time

# 片方だけにする
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --mode particles
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --mode particles --movie 広がり.gif

# 壁の透過をレイヤごとに指定（実行中は左の縦スライダと Tab / m でも変えられる）
python view_rays.py ..\test.dxf ..\結果\test_raylog.npz --layer-opacity "1=0.6,2=0.05"

# シミュレーション本体（DXF → パルス列 → インパルス応答 → 残響時間）
python procedure.py ..\test.dxf --absorption ..\absorption.csv --absorption-kind normal ^
       --out ..\結果 --rays 20000 --nref 8 --bands 8 --temperature 20 --humidity 40
```

主なオプション

| オプション | 意味 |
|---|---|
| `--absorption-kind normal\|random` | 吸音率が**垂直入射**か**残響室法**か。★取り違え注意 |
| `--assignment 条件A.json` | レイヤ→材料の対応。**CAD を編集せずに材料を差し替える** |
| `--bands 8\|6` | 周波数バンド数（8 = 63〜8kHz / 6 = 125〜4kHz） |
| `--temperature` / `--humidity` / `--pressure` | 大気条件。**音速と空気吸収が連動する** |
| `--nref` | 最大反射回数。残響時間を出すには 35 dB 減衰するまで必要 |
| `--orient-normals auto` | **閉じていれば内向きに揃え、開いた形状なら CAD のまま**（おすすめ） |
| `--orient-normals inward` | 法線を面ごとに室内側へ揃える。CAD で法線を意識せずに描いたモデル用 |
| `--two-sided` | 面の裏からの入射も当てる。`inward` が使えない場合の代替 |
| `--radius` | 受音球の半径 [m]。**経路を見つけるための網**で、小さすぎると残響時間が短く出る（下記） |
| `--volume` | 統計残響式に使う室容積 [m³]。閉じていないモデルで統計式と比べたいとき |

**受音球の半径について**（2026-08-14 に実データで判明）

受音球はエネルギーを集める器ではなく**経路を見つけるための網**（エネルギーはバックトレースが
厳密に出し、成り立たない経路は却下する）。ある経路が見つかる確率はおよそ `N·r²/(4d²)`
（N=音線数, r=半径, d=経路長）で、`d = ct` だから**遅い時刻ほど取りこぼす**。
その結果、後期のエネルギーが不足して**残響時間が短く出る**。

目安は `t_max = r√N / (2c) > 求めたい減衰時間`。372 m³ の室・音線 20 万本では
既定の 0.15 m だと 0.098 秒ぶんしか拾えず、T30 が 3〜4 割短く出た。
**半径を変えて値が動かないことを確認してから数値を使うこと**（この例では r ≥ 1.0 m で収束）。

`--out` に以下が書き出される。

| ファイル | 内容 |
|---|---|
| `*_raylog.npz` | 可視化用の音線軌跡 |
| `*_pulses.csv` | パルス列（反射回数・到来時刻・到来方向・バンド別エネルギー） |
| `*_ir.csv` | インパルス応答 |
| `*_rt.csv` | 残響指標 **EDT / T20 / T30**（オクターブバンド別）と曲率 |
| `*_吸音率と理論値.csv` | 材料別の吸音率 → **平均吸音率** → **残響時間理論値**（Sabine / Eyring / Eyring-Knudsen）。閉じた室のみ |
| `*_decay.csv` | 減衰曲線 |
| `*_spl.csv` | **帯域別の音圧レベル**（合計・A 特性・直接音・反射音）。PWL 未入力なら相対値 |
| `*_sti.csv` | **STI**（音声伝送指数・IEC 60268-16）と帯域別 MTI・変調伝達関数 |

（`app.py` / `run_project.py` から回した場合は、プロジェクトフォルダの
`結果/recN/` に受音点ごとの CSV、`結果/` 直下に受音点に依らないものが入り、
あわせて `図/recN/` に PNG が出る。名前の頭にはプロジェクト名が付く）

## 検証

```powershell
.\.venv\Scripts\python tests\test_geosim.py
```

**464 項目**（2026-08-21 時点）。解析的に答えが分かる問題（直方体の虚音源距離、
減衰率が既知の応答など）で数式レベルの正しさを確かめる。
**数式に関わるコードを変更したら必ず走らせること。**

高速化した処理は**必ず scalar 版（参照実装）との一致を見る**方針にしてある
（交差判定・受音判定・反射）。scalar 版は消さないこと。

## 動作確認用の DXF

| ファイル | 内容 |
|---|---|
| `test.dxf` | 2×3×1 m の直方体（mm 単位、音源・受音点入り、法線内向き） |
| `test2.dxf` | 閉じたポリラインで描いた平面 9 角形 ＋ 立ち上げた壁の 2 面（開いた形状） |

## Git 管理外のもの

`fortran/`（移植元）、`absorption.csv`、`参考文献/**/*.pdf`、`.venv/`、`*_view.html`。
このリポジトリは public のため、著作物と大きな生成物は置いていない。
