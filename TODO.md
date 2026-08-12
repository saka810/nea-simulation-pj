# 作業一覧（TODO）

最終更新: 2026-08-12（HBN240018）

Fortran → Python 移植の残作業リスト。着手・完了したら本ファイルを更新すること。
各モジュールの詳細は `PROGRAM_STRUCTURE.md`、経緯は `Claude履歴/` を参照。

## 現状サマリ

| # | 処理 | ファイル | 状態 |
|---|---|---|---|
| 1 | DXF読込 | `geosim/read_dxffile.py` | **スタブ**（座標未取得・常に空リスト） |
| 2 | 音線生成 | `geosim/sound_ray.py` | 実装済（2026-08-12 に転記ミス修正・動作確認済） |
| 3 | 音線追跡 | `geosim/loop_reflectionmesh.py` | **実装中**（traceff の扱いが未反映） |
| 4 | 重複削除 | `geosim/loop_deleteredundancy.py` | デモ実装済（set 方式） |
| 5 | バックトレース | `geosim/loop_noredundancy.py` | 下書き・**未接続** |
| 6 | インパルス応答 | `geosim/impulse.py` | 実装中・**未接続**（転記ミス多数） |
| 7 | 残響時間 | — | 未着手 |

ボトルネックは **1 がスタブ・3 が未確定**で、通しで一度も動かせていないこと。

---

## A. 最優先：音線追跡を完成させる（ステップ3）

元コード `backtrace.f90` 649〜677 行の traceff の読み解きは 2026-08-12 に完了。
結論は `Claude履歴/2026-08-12_HBN240018.md` に記載（`fortran/` は Git 管理外のため要参照）。

- [ ] A-1 `if inside:` の中で `history_2dim.append(list(history))` に修正
      （**コピー必須**。現状は同一リスト参照を反射ループ毎周 append している）
- [ ] A-2 `if collision:` の中で**無条件に** `history.append(mesh_nearestid)` に修正
      （現状は `if inside:` の時だけ追加していて逆。受音と衝突は独立した条件）
- [ ] A-3 面ループ内で `soundray_comesfrom` を上書きしているバグの修正（`loop_reflectionmesh.py:61`）
      元コードは `node` を一時変数として使い、最寄り面が確定してから基点を更新する
- [ ] A-4 `collision` フラグのバグ修正（`loop_reflectionmesh.py:55`）
      面ループ内で毎回 False に戻すため最後の面の結果しか残らない。
      元コードの `flag` は面ループの**前**で初期化し、当たったら立てっぱなし
- [ ] A-5 不要な `from numpy.ma.core import count` を削除
- [ ] A-6 12/04 メモの疑問コメント（疑問1〜4）を、確定した回答に書き換え

### traceff の結論（メモ）

- 受音判定（663行）が `tractmp(k+1) = jtmp`（677行）**より前**にあるのが鍵。
  受音時点の履歴は「そこまでに済んだ k 回の反射」で、これから当たる壁は含めない。
- 疑問1 → k ループの外には**出せない**。1 本の音線が別々の k で複数回受音しうるため。
- 疑問2 → nray×nref×nref には**ならない**。`do j = 0, nref` は受音時のみ（count 回）。
  固定長配列の行まるごとコピーなだけで、Python なら `list(history)` のコピー1回で等価。
- 疑問3 → `do j = 0, k` で**良い**。k+1 以降は必ず 0。
  0 埋めが要るのは Fortran が固定長＋`traceffn`（有効反射回数）で行長管理しているため。
- 疑問4 → inside と collision を**統合してはいけない**。両者は独立。

---

## B. 次点：DXF読み込みの実装（ステップ1）

- [ ] B-1 `read_dxffile.read()` の本実装（元コード 132〜283 行）
      グループコード 70 のフラグで頂点(192)/面(128)を判別、レイヤ名→吸音材ID、
      mm→m 変換、縮退面（法線ゼロ）の除去、外積による法線計算と正規化
- [ ] B-2 `ezdxf` ライブラリを使うか自前パーサかの方針決定
- [x] B-3 `mesh.py` の `vertexes` が (3,1,3) になっている形状バグの修正 → 2026-08-12 完了
- [ ] B-4 吸音率CSV（absorption.csv）読み込みの実装（`procedure.py` にメモのみ）

---

## C. 動作確認の足場

- [ ] C-1 直方体など単純形状での通し実行（1→2→3→4）と結果検証
- [ ] C-2 `geosim/` のベア名インポート問題（現状 cwd 依存）— 実行方法を決める
      ※ CLAUDE.md 記載どおり、相対インポートに変えるなら要相談

---

## D. 既知の転記ミス修正 → **2026-08-12 完了**

- [x] D-1 `sound_ray.py`：`np.zeros(ray_number, 3)` → `np.zeros((ray_number, 3))`
- [x] D-2 `sound_ray.py`：`sound_rays(i, 2)` → `sound_rays[i, 2]`（2 箇所）
- [x] D-3 `sound_ray.py`：`energy_decay` の反射係数、第2項 `(1+√(1-α))` → `(1−√(1-α))`（2 箇所）
      旧実装は約分により吸音率が効かず、垂直入射で常に完全吸音になっていた
- [x] D-4 `receiver_sphere.py`：条件②の比較対象を垂線距離 → 射影距離（元コード 663 行）
- [x] D-5 `mesh_method.py`：`innerproduct_from3vertexes` の作業配列 (3,2) → (2,3)
- [x] D-6 `mesh_method.py`：`collision_detection` の `< 0` → `<= 0`（元コード 625 行に一致）

元の数式（局所反応性壁面の斜入射反射率）は検証の結果**正しい式**と確認済み。
導出と検証結果は `PROGRAM_STRUCTURE.md` の `sound_ray.py` 節を参照。

### D から派生した新規課題

- [ ] D-7 残響室法吸音率への対応。α > 1 になりうるため `sqrt(1 - α)` が NaN になる。
      Paris の式を逆解きして法線入射吸音率 α₀ または規格化インピーダンス z を求める前処理が必要。
      NEA では残響室法吸音率を使うことが多いとのことなので、実運用では避けて通れない見込み

---

## E. 後半パイプライン（A が固まってから）

- [ ] E-1 `loop_noredundancy.py` の整理
      （リスト/ndarray の型不整合、`energy_decay` の引数不足、衝突判定ブロックの位置）
- [ ] E-2 `procedure.py` のステップ5・6のコメントアウト解除と接続
- [ ] E-3 `impulse.py` の転記ミス修正
      （Hamming 窓の `+`→`*`、伝達関数の `+`→`*`、空気吸収の括弧、負周波数の添字、逆FFT呼び出し欠落）
- [ ] E-4 `impulse.py`：`df`（周波数離散化幅）の `transfer_function` への受け渡し
- [ ] E-5 6バンド吸音率 → 32バンド展開表の移植（Fortran 側にあり未移植）
- [ ] E-6 自作FFT → `np.fft`、バンドパス → `scipy.signal.firwin` への置き換え判断
- [ ] E-7 周波数バンド数の不整合解消（`procedure.py` は 8 バンド、吸音率は 6 バンド）
- [ ] E-8 残響時間算出（`ipls2rt_fortran.f90` 相当）の移植

---

## F. その他（優先度低）

- [ ] F-1 OpenMP 並列化相当の検討（NumPy ベクトル化 / numba が候補）
