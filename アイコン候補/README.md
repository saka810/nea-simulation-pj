# NEA GeoSim（ネアシム）アプリアイコン候補

作業者：**Codex** ／ 作成日：2026-11-16

まず [一覧.png](一覧.png) を開いてください。構図・題材が異なる8案です。見本の「線画版」に合わせ、濃紺の太めの輪郭線、丸い線端、白い余白と少数の差し色で制作しました。ロゴ案は線をさらに絞っています。

## 保存内容

- `原寸PNG/NG01.png`〜`NG08.png`：組込みimage_genの生成結果をそのままコピー（各1254×1254px）。
- `PNG/NG01.png`〜`NG08.png`：各512×512px。
- `ICO/NG01.ico`〜`NG08.ico`：16・20・24・32・40・48・64・128・256pxの9サイズを各ファイルに格納。
- `一覧.png`：8案の大きな見本、16・32・48pxの実寸縮小、暗い背景に置いた同じ画像。100%表示で確認してください。
- `コード/export_icons.cjs`：縦横比を保つリサイズ、ICOの作成・検査、一覧の生成。
- `コード/デザイン.json`：各案の最終プロンプトと原寸の生成元パス。
- `コード/検証結果.json`：原寸寸法・SHA-256・ICOのサイズ検査結果。

見本と同様に**白い背景を含む不透明画像**です。透明化はしていません。暗色欄も白背景を含む同じPNGを置いており、透明背景版の比較ではありません。

## 案の一覧

| ID | 名前 | 狙い | 差し色（共通：濃紺＋白） |
|---|---|---|---|
| NG01 | 反射の道すじ | 室の断面を一筆の音線が反射し、受音球へ届く。 | ティール・コーラル |
| NG02 | 鏡の向こうの音源 | 壁の反対側の虚音源と反射経路を作図する。 | バイオレット・オレンジ・ティール |
| NG03 | 音線から応答へ | 放射する音線と減衰する応答の棒を組み合わせる。 | シアン・アンバー |
| NG04 | ホールをたどる | ホールの輪郭と、内部を走る一本の反射経路。 | ブルー・コーラル |
| NG05 | 響きを聴く | 耳に届く折れ線で、計算から可聴化までを表す。 | ティール・マスタード |
| NG06 | Nの反射 | Nの対角線を音線と見立てた、簡潔なロゴ案。 | ティール・オレンジ |
| NG07 | Gの受音点 | Gの輪郭に受音点を組み込み、幾何と音を象徴する。 | ブルー・ライム |
| NG08 | ネアシムのネ | 「ネ」の払いを音線に見立てた、日本語のロゴ案。 | コーラル・バイオレット |

## 小サイズの確認

- **NG01 反射の道すじ**：16pxでも室の枠と折れ線は残る。受音球と音源の違いは32px以上が見やすい。
- **NG02 鏡の向こうの音源**：16pxでは虚音源の点線と輪が弱い。作図の意味を伝える用途は32〜48px以上が向く。
- **NG03 音線から応答へ**：16pxでは応答の棒の本数と放射線の分離が読みにくい。32pxから題材が伝わりやすい。
- **NG04 ホールをたどる**：16pxではホールの奥行き線と反射経路が重なりやすい。32〜48px以上が向く。
- **NG05 響きを聴く**：耳の形は比較的残る。16pxでは音源と反射の細かな意味は読み取りにくい。
- **NG06 Nの反射**：8案中、小サイズで最も明確。16pxでもNと対角線を識別しやすい。
- **NG07 Gの受音点**：16pxでもGの輪郭は残る。中央の受音点と経路は32px以上が見やすい。
- **NG08 ネアシムのネ**：文字主体で識別しやすい。16pxでは受音球の内側と上の短い線が弱く、32px以上を推奨。

機能の説明を重視するならNG01、16pxでの識別を重視するならNG06、幾何と音のロゴとしてはNG07が比較しやすい候補です。いずれも候補段階で、アプリへの組込みは行っていません。反射の図は機能を示す象徴表現で、物理計算の厳密な作図ではありません。

## 作り方・再出力

1. 指定の見本フォルダの `README.md`、`線画版/README.md`、比較画像2点、全68案の最終プロンプトJSONを読み取り、作風と指示の構造を確認。
2. 組込み **image_gen** を1案ずつ使用。生成は新規画像として行い、見本画像の編集・切り抜きは行っていません。最終プロンプトは以下とJSONに保存。N案は対角線の向きを直して再生成し、正しいNの最終結果を採用。
3. 原寸を無加工で保存。SharpのLanczos3で全体を正方形に収め、余白と縦横比を保持して各サイズへ変換。切り抜き・縦横別の引き伸ばし・描き足し・色変更は行っていません。
4. ICOの9フレームはPNG圧縮で格納。全フレームを読み戻し、復号・寸法・サイズ数を検証。一覧を描画し、縮小時の見え方と文字の配置を目視確認。

Node.jsとSharpが利用できる環境では、以下で保存済み原寸から再出力できます（画像の再生成は行いません）。

```powershell
node .\アイコン候補\コード\export_icons.cjs
```

Sharpが通常のモジュール探索先にない場合は、Sharpパッケージの絶対パスを第1引数に渡してください。今回利用した既存環境の例：

```powershell
node .\アイコン候補\コード\export_icons.cjs "C:\Users\sakayoshi\AppData\Local\OpenAI\Codex\runtimes\cua_node\bfa21d187eb00b05\bin\node_modules\sharp"
```

一覧の日本語描画にはWindowsの `C:/Windows/Fonts/YuGothM.ttc` を使用します。別環境ではスクリプトのフォント指定を調整してください。出力先はこのフォルダ内に固定されています。原寸は上書きせず、PNG・ICO・一覧・検証結果だけ再出力します。

## 参照・作業範囲

見本：`C:\Users\sakayoshi\Desktop\Claude\資料\プログラム一覧・アイコン候補`（読み取りのみ）。

リポジトリのAGENTS.mdと直近履歴を確認。今回の個別指定「他のファイルは変更しない」「git commitしない」を優先し、履歴ファイルへの追記・Git同期・アプリ本体の変更は行っていません。シェル実行環境に制限があったため、ファイル操作と画像変換はNode.jsの実行環境で実施しました。APIキーやCLIによる画像生成は使用していません。

## 各案の最終プロンプト

英語は画像生成へ渡した原文です。各案の日本語の狙い・配色は上表に記載しています。

### NG01 反射の道すじ

- 狙い：室の断面を一筆の音線が反射し、受音球へ届く。
- 配色：濃紺・白、ティール・コーラル。
- [512px PNG](PNG/NG01.png) ／ [ICO](ICO/NG01.ico) ／ [原寸](原寸PNG/NG01.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: A square room cross-section as a simple navy open architectural outline, viewed straight-on, not a cube. Inside one continuous thick teal zigzag sound ray travels from a small coral source dot at lower left, bounces once off top wall then once off bottom wall, and ends in a large white receiver circle with navy border at right. Exactly one ray and two bounces, no grid or extra arrows. No letters.
```

### NG02 鏡の向こうの音源

- 狙い：壁の反対側の虚音源と反射経路を作図する。
- 配色：濃紺・白、バイオレット・オレンジ・ティール。
- [512px PNG](PNG/NG02.png) ／ [ICO](ICO/NG02.ico) ／ [原寸](原寸PNG/NG02.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: A minimal image-source construction diagram: one horizontal navy reflecting wall across the lower middle, orange solid source dot above left, violet hollow dashed ghost-source circle directly below the wall, teal receiver ring above right. One V-shaped orange specular reflected path touches the wall once; a sparse violet dashed line extends the segment behind the wall to the ghost source. No axes, labels, extra dots or letters. Spacious geometrical composition.
```

### NG03 音線から応答へ

- 狙い：放射する音線と減衰する応答の棒を組み合わせる。
- 配色：濃紺・白、シアン・アンバー。
- [512px PNG](PNG/NG03.png) ／ [ICO](ICO/NG03.ico) ／ [原寸](原寸PNG/NG03.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: On the left one small solid navy sound-source circle sends exactly three thick cyan straight rays fanning toward upper right. On the right a compact impulse-response emblem: exactly four round-capped vertical amber stems of decreasing height sitting on one short navy baseline. Integrate them as one balanced diagonal composition, source lower left and decay upper right, no extra axes or curves. No letters.
```

### NG04 ホールをたどる

- 狙い：ホールの輪郭と、内部を走る一本の反射経路。
- 配色：濃紺・白、ブルー・コーラル。
- [512px PNG](PNG/NG04.png) ／ [ICO](ICO/NG04.ico) ／ [原寸](原寸PNG/NG04.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: A simplified architectural concert-hall cutaway, a broad trapezoidal roof and two side edges with open front in a spare navy wireframe, flat 2D linework suggesting an axonometric volume with no shading. Exactly one bold coral angular sound ray travels from a blue source dot at front left to the roof and then down into a large outlined receiver ring at front right. No seats, grids, mesh subdivisions, decorative waves or letters.
```

### NG05 響きを聴く

- 狙い：耳に届く折れ線で、計算から可聴化までを表す。
- 配色：濃紺・白、ティール・マスタード。
- [512px PNG](PNG/NG05.png) ／ [ICO](ICO/NG05.ico) ／ [原寸](原寸PNG/NG05.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: One large recognizable human ear drawn using only two smooth navy contour paths on right. A single thick teal sound ray from a mustard source circle at left makes one angular upward bounce before entering the ear center. No wavy audio waveform: the connecting ray is geometric straight segments. Strong minimal silhouette for 16px, no architectural frame or letters.
```

### NG06 Nの反射

- 狙い：Nの対角線を音線と見立てた、簡潔なロゴ案。
- 配色：濃紺・白、ティール・オレンジ。
- [512px PNG](PNG/NG06.png) ／ [ICO](ICO/NG06.ico) ／ [原寸](原寸PNG/NG06.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: One large uppercase N constructed from exactly three bold rounded monoline strokes. Two navy vertical wall-like strokes and one teal diagonal sound-ray stroke descending from UPPER LEFT to LOWER RIGHT, with one small orange endpoint dot at the LOWER RIGHT. Integrate ray into the letter, no separate arrow or waveform, no enclosing frame. Extra sparse iconic monogram with generous whitespace, legible at 16px. Only letter N, no other text. The diagonal MUST connect the TOP of the LEFT vertical to the BOTTOM of the RIGHT vertical, like the standard Latin N; never mirror it or draw Cyrillic И.
```

### NG07 Gの受音点

- 狙い：Gの輪郭に受音点を組み込み、幾何と音を象徴する。
- 配色：濃紺・白、ブルー・ライム。
- [512px PNG](PNG/NG07.png) ／ [ICO](ICO/NG07.ico) ／ [原寸](原寸PNG/NG07.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: A single large navy open circular G contour with a short horizontal inward crossbar at right. A blue straight ray comes diagonally from the lower-left opening into the center and terminates at one large lime receiver dot. Treat G as an abstract geometric acoustic chamber; retain a clean recognizable uppercase G silhouette. Exactly one circle contour, one crossbar, one ray and one dot. Flat sparse rounded monoline, no other text.
```

### NG08 ネアシムのネ

- 狙い：「ネ」の払いを音線に見立てた、日本語のロゴ案。
- 配色：濃紺・白、コーラル・バイオレット。
- [512px PNG](PNG/NG08.png) ／ [ICO](ICO/NG08.ico) ／ [原寸](原寸PNG/NG08.png)

```text
Use case: logo-brand. Asset type: Windows engineering app icon for NEA GeoSim, a geometric room acoustics ray-tracing and image-source backtrace simulator. Create one clean flat monoline illustration in the approved style: consistent medium-bold dark navy #102C54 contours, rounded caps and joins, warm-white opaque square artboard #FFFCF8. White negative space dominates. Restrained flat accent colors, no enclosing colored tile or badge. Entirely two-dimensional: no shadows, gradients, texture, glow, highlights, bevels or realistic materials. Center a large memorable silhouette with 8 percent safe margins. Deliberately sparse paths and large gaps to survive 16px and 32px. Professional contemporary engineering software icon. No labels, numerals or watermarks. Subject: One large stylized Japanese katakana ネ as the sole motif, recognizable and typographically correct, built with a few navy round-capped monoline strokes. Make its descending right diagonal stroke coral like a sound ray, ending in a modest violet outlined receiver circle. Keep the top short tick and sweeping diagonal distinct with wide gaps. No frame, no waveform, no English words, no other text. Very sparse Japanese engineering app monogram, designed to stay readable at 16px.
```

