# このリポジトリでの作業ルール

## 絶対ルール

- **セッション履歴の記録**：複数端末（PC）でこのリポジトリを clone して Claude Code を使うため、
  `Claude履歴\` フォルダに、やり取りを二段階（逐次ログ＋要約整理）で記録すること。
  - ファイル名は `YYYY-MM-DD_PC名.md`（PC名は環境変数 `%COMPUTERNAME%` / `$env:COMPUTERNAME` で取得）。
  - 同日・同PCで複数セッションがある場合は、同じファイルに追記していく。
  - **①逐次ログ**：明示的な依頼がなくても、やり取り（1ターン）ごとに「### 逐次ログ」節へ
    1〜数行程度の簡潔な記録（何を依頼され、何を行ったか）をその都度追記する。
    要約を頼み忘れても記録が失われないようにするための保険。
  - **②要約整理**：ユーザーから依頼があったタイミング、またはトピックが大きく切り替わる・
    セッションが一区切りつくタイミングで、その日蓄積された逐次ログを読み返し、
    冒頭の「本日のやり取り」節としてまとまった要約に書き直す（生ログは残したままでよい）。
  - 目的：別端末からアクセスした際に、他端末での作業内容・direction をすぐ把握できるようにするため。
  - このルールは新規セッション・別PCからのアクセスでも適用すること。

- **履歴の Git 同期**：OneDrive のような自動同期がないため、履歴の受け渡しは Git で行う。
  リモートは `origin` = https://github.com/saka810/nea-simulation-pj
  - **セッション開始時**：まず `git pull` して、他端末で書かれた履歴を取り込む。
    そのうえで `Claude履歴\` の直近ファイル（日付の新しいもの数件）を読み、前回までの流れを把握する。
  - **セッション終了時／一区切りごと**：履歴ファイルを commit して push する。
    - 履歴だけを独立させてコミットする（コード変更と混ぜない）：
      `git add Claude履歴 && git commit -m "履歴: YYYY-MM-DD PC名" && git push`
    - コードも変更した場合は、コード側を先に別コミットにしてから履歴コミットを積む。
  - **コンフリクトが起きたら**：履歴ファイルは追記主体なので、どちらの内容も消さずに
    両方を時系列で残す形でマージする（勝手に片方を捨てない）。
  - push が失敗した場合（認証切れ・ネットワーク等）はユーザーに報告する。黙って放置しない。

## リポジトリの目的

- 日本環境アメニティ株式会社（NEA）のシミュレーションPJ。
- 幾何音響シミュレーション（音線法＋虚音源バックトレース）の **Fortran → Python 移植**。
- 元コードは Fortran 3 ファイル：`backtrace.f90` / `make_ipls_freq_monaural_fortran.f90` /
  `ipls2rt_fortran.f90`。移植先は `geosim/` パッケージ。

## リポジトリ構成

- `geosim\` … Python 移植版のパッケージ（本体）
  - `procedure.py`（全体フロー）／`read_dxffile.py`（DXF読込）／`sound_ray.py`（音線生成・反射）／
    `mesh.py`・`mesh_method.py`（メッシュと交差判定）／`receiver_sphere.py`（受音判定）／
    `loop_reflectionmesh.py`（音線追跡）／`loop_deleteredundancy.py`（重複経路削除）／
    `loop_noredundancy.py`（バックトレース）／`impulse.py`（インパルス応答合成）
  - ビューアは 2 本並存：`view_model.py`（HTML+WebGL を書き出す。依存なし・共有向き）と
    `view_model_gui.py`（PyVista のネイティブウィンドウ。Python から操作しやすい）。
    **どちらか片方に寄せない**（用途が違う）。
- `requirements.txt` / `requirements-lock.txt` / `.python-version` … 実行環境の定義（上記参照）
- `fortran\` … 移植元の Fortran コード。**`.gitignore` で除外されており Git 管理外**。
  別端末には自動で配られないので、参照するコードの内容を履歴や `PROGRAM_STRUCTURE.md` 側に
  書き残しておくと他端末で追いやすい。
- `PROGRAM_STRUCTURE.md` … 各モジュールの役割・Fortran との変数対応・実装状況・既知のバグ一覧。
  **コードに手を入れたらこの文書も合わせて更新する。**
- `tests\test_geosim.py` … 数値検証（34 項目）。pytest 不要、素の Python で走る。
  **数式に関わるコードを変更したら必ず走らせる**：`.venv\Scripts\python tests\test_geosim.py`
- `TODO.md` … 作業一覧（A〜G にグルーピング）。着手・完了したらチェックボックスを更新する。
- `docs\技術説明書.md` … プログラムのフローと数式の意味の解説資料。
  数式の根拠は書籍『建築音響物理学』2.2 節。**数式に関わるコードを変更したらこの文書も更新する。**
- `docs\出力・可視化方針.md` … 最終的に作りたい 6 種の出力と、そのために蓄積すべきデータ。
  **新しい計算段階を書くときは、この文書のチェックリストを確認する。**
- `docs\DXFデータの作り方.md` … CAD 側の作図ルール（レイヤ＝吸音材、単位、法線の向き、src/rec）。
- `data\absorption_sample.csv` … 吸音率テーブルのサンプル（材料名 + a1〜a6 の形式）。
- `absorption.csv` … 元 Fortran パッケージ付属の吸音率テーブル（ID + 材料名 + a1〜a6、CP932）。
  **`fortran\` と同じくローカルのみ**（`.gitignore` で除外。public リポジトリのため保留中）。
  読み込みは両形式に対応済み。
- `test.dxf` … 動作確認用の DXF（2×3×1 m の直方体、mm 単位、音源・受音点入り）。
- `test2.dxf` … 動作確認用の DXF（閉じたポリラインで描いた平面9角形＋立ち上げた壁の2面、mm 単位）。
- `参考文献\` … 論文・書籍の PDF。**PDF 本体はローカルのみ**（`.gitignore` で除外。public リポジトリのため）。
  フォルダ構成と文献一覧（`参考文献\README.md`）だけを共有する。原本は OneDrive の `06_参考文献\`。
  `fortran\` と同じく他端末には配られないので、**参照した内容は README や技術説明書に書き残す**こと。
- `Claude履歴\` … セッション履歴（上記「絶対ルール」参照）

## 実行環境

- **Python は 3.10.11**（チーム方針）。`.python-version` に明記してある。
  端末によっては既定の `python` が 3.10 以外を指すので、**venv は必ずインタプリタを明示して作る**：
  `py -3.10 -m venv .venv`。
- 依存は `requirements.txt`（numpy / scipy / matplotlib / pyvista）。
  厳密に揃えたいときは `requirements-lock.txt`。手順は README.md「環境構築」。
- `.venv/` は `.gitignore` 済み。**端末ごとに作り直す**（中身は共有しない）。
- 新しいライブラリを入れたら **`requirements.txt` と `requirements-lock.txt` を必ず更新**する。

## 作業上の注意

- 現状の `geosim/` は `from mesh import Mesh` のような**ベア名インポート**のため、
  `geosim/` を cwd にした直接実行のみを想定している。パッケージとして import する形に変える場合は
  相対インポートへの一括変更が必要（勝手に変えず、必要になった時点で相談する）。
- パイプラインは **1〜7 が全部つながっている**（2026-08-14）。
  DXF 読込 → 音線生成 → 音線追跡 → 重複削除 → バックトレース → インパルス応答 → 残響時間。
  ただし**残響時間を実用的に出すには反射回数が足りない**（F-1 の高速化が前提。TODO E-10）。
- **元コードから意図的に変えた箇所**。変更するときは PROGRAM_STRUCTURE.md の該当節を読むこと。
  - `loop_noredundancy.image_sources` … 鏡像の距離に `abs()` を使わず符号付きにした
    → **判断待ち**。`docs/議論_鏡像の式のabs問題.md` に詳しく書いてある
  - `reverberation` … 巡回畳み込みの自作 FIR → Butterworth（IEC 61260 準拠）。
    元のやり方だと回り込みで減衰曲線に床ができ T30 が測れない
  - `impulse` … 信号処理を scipy に置き換え（`rfft`/`irfft`/`firwin`）。
    フィルタ長を実用的にしたので**元コードにあった 1.49 秒の遅れは無い**
  - `atmosphere` … 音速と空気吸収を定数・近似式でなく温度・湿度・気圧から計算（ISO 9613-1）
- **周波数バンドは既定 8（63〜8k Hz）**。6（125〜4k）も可。`band_number` で切り替える。
- **交差判定はベクトル化してある**（`mesh_method.FaceArrays`）。
  scalar 版（`collision_distance` など）は参照実装・一致確認の基準として残してあるので**消さない**。
  交差判定に手を入れたら、必ず `tests/test_geosim.py` の
  「ベクトル化した交差判定（scalar 版との一致）」が通ることを確認する。
- **表は「周波数を横」に並べる**（2026-08-17 ユーザー判断。`geosim/table.py` に共通ルール）。
  グラフにしたときの横軸が周波数なので、CSV も画面の表もそれに揃える。
  縦だと Excel でグラフを作るたびに行と列を選び直す手間が要る。
  - 1 行が指標のとき（`rt.csv` / `rt_statistical.csv` / `clarity.csv`）は
    `table.write_frequency_table()` を使う。1 列目が `項目`、2 列目以降が周波数
  - 1 行が周波数以外（時刻・材料・経路）のときは、周波数を**列名**にする。
    列名は `table.band_column()` で作る（`alpha_125Hz` / `energy_125Hz` / `decay_125Hz_db`）。
    **番号（`energy_1`）にしない**。バンド数 6 と 8 で意味が変わって読み違えるため
  - 読むときは `table.read_frequency_table()`（**古い縦向きも読める**）
  - **新しい出力を書くときもこの向きにすること。**
    どうしても縦にしたい事情が出たら、勝手に決めずにユーザーに確認する
- **残響指標は EDT / T20 / T30 を出す**（`reverberation.decay_measures`）。
  60 dB 減を厳密に見る運用ではない。
- **統計残響式は Sabine / Eyring / Eyring-Knudsen の 3 つ**（`reverberation.STATISTICAL_LABELS`）。
  アイリングの式は `-S ln(1-ᾱ)` までで、**空気吸収 `4mV` を足した形はヌードセンの寄与**。
  2 つを並べると差がそのまま空気吸収の効きになる。
  ミリントン（ミリントン・セッテ）は実務で使われないので**落とした**（2026-08-17）。
- **吸音率は「垂直入射」か「残響室法」かを必ず区別する**。取り違えると吸音を大きく誤る。
  CSV に `# kind: normal|random` を書くか `--absorption-kind` で指定する。
  残響室法なら Paris の式で自動変換される（`geosim/absorption.py`）。
- Fortran の行番号で議論することが多い（例：「backtrace.f90 524行〜」）。
  コミットメッセージや履歴にも、対応する元コードの行番号を書き添えると後から追いやすい。
