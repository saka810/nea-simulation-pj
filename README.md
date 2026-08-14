# nea-simulation-pj

NEA（日本環境アメニティ株式会社）のシミュレーション PJ。
幾何音響シミュレーション（音線法 ＋ 虚音源バックトレース）の **Fortran → Python 移植**。

- 各モジュールの役割・Fortran との対応・実装状況 … [PROGRAM_STRUCTURE.md](PROGRAM_STRUCTURE.md)
- 作業一覧 … [TODO.md](TODO.md)
- 数式とフローの解説 … [docs/技術説明書.md](docs/技術説明書.md)
- **要判断**：鏡像の式の `abs()` 問題 … [docs/議論_鏡像の式のabs問題.md](docs/議論_鏡像の式のabs問題.md)
- 出力・可視化の方針 … [docs/出力・可視化方針.md](docs/出力・可視化方針.md)
- CAD 側の作図ルール … [docs/DXFデータの作り方.md](docs/DXFデータの作り方.md)

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
python -c "import numpy, scipy, matplotlib, pyvista; print('ok')"
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

流れは **条件入力 → 法線の確認 → 計算 → 可視化**。

1. **条件入力**（`setup_window.py`）… モデル DXF・吸音率 CSV・保存先フォルダ・
   音線数・受音球の半径・温度湿度などを入力する。
   既存フォルダを選べば前回の条件を読み込む
2. **法線の確認**（`normal_editor.py`、「法線を確認…」ボタン）…
   緑=内向き（OK）／赤=外向き（要反転）／灰=判定できない、で色分けして表示。
   枠で囲んだ面を反転、数字キーでレイヤごと反転、`s` で保存
3. **計算**（「計算する ▶」）… 結果 CSV と図 PNG がプロジェクトフォルダに保存される
4. **可視化**（`view_rays.py`）… 音線と音粒子を `Tab` で切り替えて見る

プロジェクトフォルダの中身:

```
project.json          条件（次に開いたときそのまま復元される）
normals.json          法線の反転指定
結果/  pulses.csv  ir.csv  rt.csv  rt_statistical.csv  clarity.csv
       decay.csv  surface.csv  raylog.npz
図/    impulse_response.png  decay.png  reverberation.png
       clarity.png  absorption.png  pulses.png
rec2/  受音点が複数ある場合、2 点目以降は同じ構造で枝分かれ
```

**インパルス応答の図は最大値で正規化してある。** 絶対振幅の校正が未了なので、
そのままの値には意味が無いため（TODO E-11）。

### コマンドラインから使う

```powershell
cd geosim

# プロジェクトの条件で計算する（GUI を使わずに回す）
python run_project.py "C:\...\プロジェクト"

# 法線の確認だけ開く
python normal_editor.py "C:\...\プロジェクト"

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
| `*_rt_statistical.csv` | **統計残響式**（Sabine / Eyring / Millington）。閉じた室のみ |
| `*_decay.csv` | 減衰曲線 |

（`app.py` / `run_project.py` から回した場合は、プロジェクトフォルダの `結果/` に
`pulses.csv` `ir.csv` … という素の名前で入り、あわせて `図/` に PNG が出る）

## 検証

```powershell
.\.venv\Scripts\python tests\test_geosim.py
```

解析的に答えが分かる問題（直方体の虚音源距離、減衰率が既知の応答など）で
数式レベルの正しさを確かめる。**数式に関わるコードを変更したら必ず走らせること。**

## 動作確認用の DXF

| ファイル | 内容 |
|---|---|
| `test.dxf` | 2×3×1 m の直方体（mm 単位、音源・受音点入り、法線内向き） |
| `test2.dxf` | 閉じたポリラインで描いた平面 9 角形 ＋ 立ち上げた壁の 2 面（開いた形状） |

## Git 管理外のもの

`fortran/`（移植元）、`absorption.csv`、`参考文献/**/*.pdf`、`.venv/`、`*_view.html`。
このリポジトリは public のため、著作物と大きな生成物は置いていない。
