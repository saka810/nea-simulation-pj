"""条件表 ― 「レイヤー名 → 材料番号」の対応と吸音率を 1 つの Excel で持つ。

## 経緯（依頼 2026-08-21）

> 壁の吸音材条件を変える場合に、CAD データを変えて、それをまた読み込んでもらって、
> っていうのが大変な気がします。レイヤー名と壁材料の一覧表（条件表）を作って、
> それを読み込ませる形式に変えたい。

> DXF ファイルを読み込んでレイヤー名等を抽出して、条件シートを作る仕組みを作って欲しい。
> 最初の設定画面で、DXF 読み込んだうえで、条件シート作成 みたいなボタンを押すと、
> 条件表を作成してください。準備している場合は、そのボタンを押さずに
> そのファイルの場所を選択します。

> 条件表に吸音率のデータを載せるのは、ミスリードを誘発するのでやめて欲しいです。
> 材料名については一応載せておこうとは思いますが、材料番号入力列を作って、
> その番号をもってくるようにしてください。材料名だと、全角・半角とかの
> ミスタイプを誘発する可能性が高いので。

> 同ファイルで、条件違いをシートを分けて記載したいのと、「吸音率」シートを一緒に
> 同梱して、データの行き来のやり易さ、および、PJ 特有の吸音データなどの
> 反映がやりやすいかな、と思っています。

## かたち

    プロジェクトフォルダ/条件表.xlsx

      シート「吸音率」      ← PJ で使う材料の一覧（**番号・材料名・吸音率**）
        番号,材料名,63,125,250,500,1000,2000,4000,8000,種類,備考
        1,コンクリート,0.01,0.02,...,random,
        11,吸音板,0.05,0.09,...,random,

      シート「現状」        ← 条件 1 つ ＝ シート 1 枚。**シート名が条件名**
        区分,レイヤー名,材料番号,材料名（参考）,面数,面積_m2
        レイヤ,01__研修室_壁_扉,1,=VLOOKUP(...),2,3.255

      シート「吸音追加案」  ← 同じ形。番号を書き換えるだけで条件が変わる

### ★入力は「材料番号」だけ

材料名で書くと**全角・半角のミスタイプを誘発する**（ユーザー指摘）。
そこで**番号を入力する列**を作り、材料名は隣に**参考として出す**。
Excel 上では `VLOOKUP` で吸音率シートを引くので、番号を入れれば名前がその場で出る。
読み込むときは**番号の列しか見ない**（名前の列は人が読むためのもの）。

### ★吸音率は条件シートに載せない

以前は参考として α を並べていたが、**ミスリードを誘発する**のでやめた
（ユーザー指摘 2026-08-21）。吸音率は「吸音率」シートにだけ置く。
そのぶん、**PJ 固有の吸音データをこの 1 ファイルに閉じ込められる**
（`absorption.csv` を別に持ち回らなくてよい）。

### 2 種類の行

| 区分 | 意味 | 書き換えられるか |
|---|---|---|
| `レイヤ` | DXF のレイヤに材料番号を割り当てる（**これが入力**） | ○ |
| `面ごとの指定` | `face_editor` で面を選んで貼った材料の**記録** | ×（記録） |

## 更新の約束（★ここを崩さないこと）

**利用者が書いた「材料番号」は絶対に上書きしない。**
`update()` は既存のシートを読み、番号の列はそのまま残して、
面数・面積と参考の名前だけ書き直す。新しいレイヤの行だけを足し、
モデルから消えたレイヤは区分を `（モデルに無し）` にして残す
（消すと「前はどう設定していたか」が分からなくなる）。

## 古い CSV も読める

`材料条件表.csv`（2026-08-21 の日中に作った形）もそのまま読める。
そちらは 3 列目が材料名で、番号の列は無い。
"""

import csv
import os

import numpy as np

import absorption as ab

# 既定のファイル名
CONDITION_BOOK = "条件表.xlsx"
CONDITION_FILE = "材料条件表.csv"      # 昔の形（CSV）。読むためだけに残す

# 吸音率のシート名（**条件のシートではない**ので一括計算からは外す）
ABSORPTION_SHEET = "吸音率"
# 条件として扱わないシート名（説明や凡例を置けるように）
RESERVED_SHEETS = {ABSORPTION_SHEET, "説明", "凡例", "メモ", "readme", "README"}

# 区分の名前
SECTION_LAYER = "レイヤ"
SECTION_FACE = "面ごとの指定"
SECTION_GONE = "（モデルに無し）"

# 条件シートの見出し。★3 列目が入力（材料番号）、4 列目は参考（材料名）
HEADER = ["区分", "レイヤー名", "材料番号", "材料名（参考）", "面数", "面積_m2"]
# 昔の CSV の見出し（3 列目が材料名だった）
LEGACY_HEADER = ["区分", "レイヤー名", "材料名", "面数", "面積_m2"]

# 吸音率シートの見出し
ABSORPTION_HEADER_HEAD = ["番号", "材料名"]
ABSORPTION_HEADER_TAIL = ["種類", "備考"]

# 既定で作る条件シートの名前
FIRST_SHEET = "現状"


# ------------------------------------------------------------------------------
# 場所とシート
# ------------------------------------------------------------------------------

def path(project):
    """使う条件表のパス（`project.condition_csv`、無ければ既定名）。"""
    return project.condition_path


def sheet_of(project):
    """使う条件シートの名前。指定が無ければ最初の条件シート。"""
    if project.condition_sheet:
        return project.condition_sheet
    found = sheets(path(project))
    return found[0] if found else None


def exists(project):
    return os.path.exists(path(project))


def is_book(file_name):
    return bool(file_name) and str(file_name).lower().endswith((".xlsx", ".xlsm"))


def sheets(file_name):
    """条件表の中の**条件シートの名前**を並び順に返す（吸音率などは除く）。

    CSV なら [None]（シートの概念が無い）。読めなければ []。
    """
    if not file_name or not os.path.exists(file_name):
        return []
    if not is_book(file_name):
        return [None]
    try:
        from openpyxl import load_workbook
        book = load_workbook(file_name, read_only=True)
    except Exception as error:
        print(f"[条件表] {file_name} を開けませんでした: "
              f"{type(error).__name__}: {error}")
        return []
    try:
        return [name for name in book.sheetnames if name not in RESERVED_SHEETS]
    finally:
        book.close()


def discover(folder, verbose=False):
    """フォルダの中の条件を全部返す [(パス, シート名), …]。一括計算の入力。

    `条件表.xlsx` の**条件シートを 1 件ずつ**に展開する。
    昔の CSV（`材料条件表.csv` の形）も拾う（シート名は None）。
    """
    found, books, legacy = [], [], []
    for name in sorted(os.listdir(folder)):
        candidate = os.path.join(folder, name)
        if is_book(name) and not name.startswith("~$"):
            if _looks_like_condition_book(candidate):
                books.append(candidate)
                found.extend((candidate, sheet) for sheet in sheets(candidate))
        elif name.lower().endswith(".csv") and is_condition_table(candidate):
            # 既定名の CSV（`材料条件表.csv`）は**xlsx があれば無視する**。
            # 移行後に同じ内容が 2 条件として並ぶのを防ぐ（2026-08-21）
            stem = os.path.splitext(name)[0]
            if stem in ("条件表", "材料条件表"):
                legacy.append(candidate)
            else:
                found.append((candidate, None))
    if not books:
        found.extend((candidate, None) for candidate in legacy)
    elif legacy and verbose:
        print(f"[条件表] 昔の CSV は使いません（xlsx があるので）: "
              + " / ".join(os.path.basename(f) for f in legacy))
    if verbose:
        print(f"[条件表] {folder} に条件が {len(found)} 件: "
              + " / ".join(label_of(p, s) for p, s in found))
    return found


def label_of(file_name, sheet=None):
    """条件の名前。**シートがあればシート名**、無ければファイル名（拡張子なし）。"""
    if sheet:
        return str(sheet)
    stem = os.path.splitext(os.path.basename(str(file_name or "")))[0]
    return stem


def _looks_like_condition_book(file_name):
    """その xlsx が条件表かどうか（条件シートの見出しで判別する）。"""
    try:
        from openpyxl import load_workbook
        book = load_workbook(file_name, read_only=True)
    except Exception:
        return False
    try:
        for name in book.sheetnames:
            row = next(book[name].iter_rows(max_row=1, values_only=True), None)
            if not row:
                continue
            cells = [str(c).strip() if c is not None else "" for c in row]
            if cells[:2] == HEADER[:2]:
                return True
    except Exception:
        return False
    finally:
        book.close()
    return False


def is_condition_table(file_name):
    """その CSV が（昔の形の）条件表かどうか。見出し行で判別する。"""
    if not os.path.isfile(file_name):
        return False
    try:
        with open(file_name, encoding="utf-8-sig", newline="") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                cells = [c.strip() for c in line.split(",")]
                return cells[:2] == HEADER[:2]
    except (OSError, UnicodeDecodeError):
        return False
    return False


# ------------------------------------------------------------------------------
# 読む
# ------------------------------------------------------------------------------

def read(file_name, sheet=None):
    """条件表を読む。戻り値 (割り当て, 記録).

    割り当て … {レイヤー名: 材料番号}（**番号の列だけを見る**。空の行は入れない）
    記録     … [(区分, レイヤー名, 材料, 面数, 面積), …] そのままの並び
    """
    if not file_name or not os.path.exists(file_name):
        return {}, []
    rows = (_book_rows(file_name, sheet) if is_book(file_name)
            else _csv_rows(file_name))
    if not rows:
        return {}, []

    legacy = rows[0][:3] == LEGACY_HEADER[:3]      # 3 列目が材料名の昔の形
    if rows[0][:2] == HEADER[:2]:
        rows = rows[1:]

    assignment, records = {}, []
    for row in rows:
        row = list(row) + [""] * (len(HEADER) - len(row))
        section, layer = row[0].strip(), row[1].strip()
        material = row[2].strip()      # 昔の形では材料名、いまは材料番号
        records.append((section, layer, material, row[4].strip(), row[5].strip())
                       if not legacy else
                       (section, layer, material, row[3].strip(), row[4].strip()))
        if section == SECTION_LAYER and layer and material:
            assignment[layer] = material
    return assignment, records


def _csv_rows(file_name):
    with open(file_name, encoding="utf-8-sig", newline="") as f:
        return [[c for c in row] for row in csv.reader(f)
                if row and row[0].strip() and not row[0].startswith("#")]


def _book_rows(file_name, sheet=None):
    """xlsx の 1 シートを文字列の表で読む。**式ではなく入力値がほしい列だけ使う。**"""
    from openpyxl import load_workbook

    book = load_workbook(file_name, read_only=True)
    try:
        name = sheet or next((n for n in book.sheetnames
                              if n not in RESERVED_SHEETS), None)
        if name is None or name not in book.sheetnames:
            return []
        rows = []
        for values in book[name].iter_rows(values_only=True):
            if values is None or all(v is None or str(v).strip() == ""
                                     for v in values):
                continue
            rows.append([_text(v) for v in values])
        return rows
    finally:
        book.close()


def _text(value):
    """セルの値を文字列にする。**整数はそのまま**（`1.0` にしない）。"""
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def assignment_for(project, verbose=True):
    """条件表から {レイヤー名: 材料番号} を返す。無ければ project.json の指定。"""
    if not exists(project):
        if project.condition_csv:
            print(f"[条件表] 指定された条件表が見つかりません: "
                  f"{project.condition_csv}（レイヤー名で吸音材を引きます）")
        return project.assignment
    sheet = sheet_of(project)
    assignment, _ = read(path(project), sheet)
    if verbose:
        where = f"{os.path.basename(path(project))}"
        if sheet:
            where += f" のシート『{sheet}』"
        print(f"[条件表] {where} から {len(assignment)} レイヤの割り当てを読みました")
    if project.assignment and verbose:
        print("[条件表] 注意: project.json 側の割り当ては使いません（条件表を優先）")
    return assignment or None


def library_from_book(file_name, kind=None, verbose=True):
    """条件表の「吸音率」シートから材料一覧を作る。無ければ None。

    ★これがあると**吸音率のデータを条件表 1 ファイルに閉じ込められる**
    （`absorption.csv` を別に持ち回らなくてよい。ユーザー要望 2026-08-21）。
    番号は**別名として登録**するので、条件シートに番号を書けば引ける。
    """
    if not is_book(file_name) or not os.path.exists(file_name):
        return None
    try:
        from openpyxl import load_workbook
        book = load_workbook(file_name, read_only=True, data_only=True)
    except Exception as error:
        print(f"[条件表] {file_name} を開けませんでした: "
              f"{type(error).__name__}: {error}")
        return None
    try:
        if ABSORPTION_SHEET not in book.sheetnames:
            return None
        rows = [list(v) for v in book[ABSORPTION_SHEET].iter_rows(values_only=True)]
    finally:
        book.close()

    library = ab.MaterialLibrary()
    for row in rows[1:] if rows else []:
        cells = [_text(v) for v in row]
        if len(cells) < 3 or not cells[1]:
            continue
        number, name = cells[0], cells[1]
        values, material_kind = [], kind
        for cell in cells[2:]:
            if cell in ab.KINDS:
                material_kind = cell
                break
            try:
                values.append(float(cell))
            except ValueError:
                break
        if len(values) < 3:
            continue
        library.add(name, values, kind=material_kind or ab.KIND_NORMAL,
                    overwrite=True)
        if number:
            library.aliases[number] = name
    if not len(library):
        return None
    if verbose:
        print(f"[条件表] 「{ABSORPTION_SHEET}」シートから "
              f"{len(library)} 材料を読みました（{os.path.basename(file_name)}）")
    return library


# ------------------------------------------------------------------------------
# 作る・更新する
# ------------------------------------------------------------------------------

def resolve_material(layer, library, assignment):
    """そのレイヤに使われる材料の**キー**（番号か名前）を返す。引けなければ ""。

    引く順番は `read_dxffile._resolve_absorption` と同じ。
        ① 条件表／project.json の割り当て
        ② レイヤー名そのもの（材料名または別名）
        ③ レイヤー名の先頭の番号（元コードの材料 ID）
    """
    if assignment:
        name = (assignment.mapping if hasattr(assignment, "mapping")
                else assignment).get(layer)
        if name:
            return name
    if library is None:
        return ""
    if layer in library:
        return layer
    import read_dxffile as rd
    number = rd.layer_number(layer)
    if number is not None and number in library:
        return number
    return ""


def material_name(key, library):
    """材料番号（または名前）から材料名を引く。引けなければ ""。"""
    if not key or library is None:
        return ""
    if key in library.aliases:
        return library.aliases[key]
    return key if key in library.materials else ""


def create(project, model, library=None, sheet=FIRST_SHEET, verbose=True):
    """**DXF のレイヤから条件表（xlsx）を新しく作る。**

    設定画面の「条件表を作成」から呼ぶ（ユーザー要望 2026-08-21）。
    「吸音率」シートに材料一覧を、`sheet` の名前で条件シートを 1 枚作る。
    既にあるファイルは**上書きしない**（作ったものを壊さないため）。
    """
    file_name = path(project)
    if not is_book(file_name):
        # 既定は xlsx。CSV を指定されている場合はその隣に xlsx を作る
        file_name = os.path.join(os.path.dirname(file_name) or project.folder,
                                 CONDITION_BOOK)
    if os.path.exists(file_name):
        return update(project, model, library, verbose=verbose)
    _write_book(file_name, project, model, library, {sheet: {}}, verbose=verbose)
    if verbose:
        print(f"[条件表] 作りました（シート「{ABSORPTION_SHEET}」＋「{sheet}」）: "
              f"{file_name}")
    return file_name


def update(project, model, library=None, assignment=None, verbose=True):
    """条件表を更新する（面数・面積・参考の材料名を書き直す）。

    ★**利用者が書いた「材料番号」は上書きしない。**
    CSV が指定されている場合は CSV のまま更新する（昔の形を壊さない）。
    """
    # ★**拡張子で行き先を決める**。CSV を指定されているなら CSV のまま更新する
    #   （昔の形を勝手に xlsx へ移し替えない）
    file_name = path(project)
    if is_book(file_name):
        existing = {}
        for name in sheets(file_name) or [FIRST_SHEET]:
            existing[name or FIRST_SHEET] = read(file_name, name)[0]
        if assignment is not None and len(existing) == 1:
            name = next(iter(existing))
            existing[name] = existing[name] or _as_dict(assignment)
        _write_book(file_name, project, model, library, existing, verbose=verbose)
    else:
        previous, _ = read(file_name)
        if assignment is None:
            assignment = previous or project.assignment or {}
        _write_csv(file_name, project, model, library,
                   previous or _as_dict(assignment))
    if verbose:
        print(f"[条件表] 更新しました: {file_name}")
    return file_name


def _as_dict(assignment):
    if assignment is None:
        return {}
    return dict(assignment.mapping if hasattr(assignment, "mapping")
                else assignment)


def _rows_for(project, model, library, assignment):
    """条件シート 1 枚ぶんの行（区分・レイヤー名・材料番号・面数・面積）。"""
    layer_areas = getattr(model, "layer_areas", {}) or {}
    layer_counts = getattr(model, "layer_counts", {}) or {}
    rows = []
    for layer in sorted(set(layer_counts) | set(layer_areas)):
        # ★材料番号は「利用者が書いたもの」が最優先。無ければ自動で引けた番号
        key = assignment.get(layer) or resolve_material(layer, library, assignment)
        rows.append((SECTION_LAYER, layer, key, layer_counts.get(layer, 0),
                     layer_areas.get(layer, 0.0)))
    # モデルから消えたレイヤも残す（前の設定を見返せるように）
    for layer, key in sorted(assignment.items()):
        if layer not in layer_counts and layer not in layer_areas:
            rows.append((SECTION_GONE, layer, key, 0, 0.0))
    # 面ごとの指定（face_editor で貼ったもの）の記録
    for material, (count, area) in sorted(_face_records(project, model).items()):
        rows.append((SECTION_FACE, material, "", count, area))
    return rows


def _face_records(project, model):
    """面ごとの割り当ての集計 {材料名: (面数, 面積)}。"""
    face_materials = project.face_materials_for(len(model.mesh))
    if not face_materials:
        return {}
    import reverberation as rv
    areas = rv.triangle_areas(model.mesh)
    result = {}
    for index, material in face_materials.items():
        if not material or index >= len(areas):
            continue
        count, area = result.get(material, (0, 0.0))
        result[material] = (count + 1, area + float(areas[index]))
    return result


def _write_book(file_name, project, model, library, conditions, verbose=True):
    """条件表（xlsx）を書く。`conditions` は {シート名: 割り当て}。"""
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        raise ImportError("条件表（xlsx）には openpyxl が要ります。"
                          "`pip install -r requirements.txt` を実行してください")

    book = Workbook()
    book.remove(book.active)
    head_font = Font(bold=True)
    head_fill = PatternFill("solid", fgColor="FFEFF3F8")

    def decorate(sheet, widths):
        for cell in sheet[1]:
            cell.font = head_font
            cell.fill = head_fill
            cell.alignment = Alignment(horizontal="center")
        for i, width in enumerate(widths, start=1):
            sheet.column_dimensions[get_column_letter(i)].width = width
        sheet.freeze_panes = "A2"

    # ---- 吸音率シート（PJ 固有の吸音データをここに置く）----
    frequencies = ab.octave_bands(project.band_number)
    sheet = book.create_sheet(ABSORPTION_SHEET)
    sheet.append(ABSORPTION_HEADER_HEAD + [f"{v:.0f}" for v in frequencies]
                 + ABSORPTION_HEADER_TAIL)
    for number, name, material in _materials_for(library):
        sheet.append([_number_cell(number), name]
                     + [float(v) for v in material.resample(project.band_number)]
                     + [material.kind, material.note])
    decorate(sheet, [10, 34] + [8] * len(frequencies) + [10, 24])
    sheet["A1"].comment = None

    # ---- 条件シート（1 条件 = 1 シート。シート名が条件名）----
    last = len(frequencies) + 2
    for name, assignment in conditions.items():
        sheet = book.create_sheet(name)
        sheet.append(HEADER)
        for row, (section, layer, key, count, area) in enumerate(
                _rows_for(project, model, library, _as_dict(assignment)),
                start=2):
            # ★材料名は VLOOKUP で「吸音率」シートから引く。
            #   番号を書き換えたら Excel 上でその場で名前が変わる
            formula = (f'=IFERROR(VLOOKUP($C{row},{ABSORPTION_SHEET}!$A:$B,2,FALSE),"")'
                       if section != SECTION_FACE else "")
            sheet.append([section, layer, _number_cell(key), formula,
                          int(count), round(float(area), 3)])
        decorate(sheet, [14, 34, 12, 34, 8, 12])
    book.save(file_name)
    return file_name


def _materials_for(library):
    """吸音率シートに並べる材料 [(番号, 材料名, Material), …]。番号順。"""
    if library is None:
        return []
    numbers = {}
    for alias, target in library.aliases.items():
        if alias.strip().isdigit():
            numbers.setdefault(target, alias)
    rows = []
    for name in library.names():
        rows.append((numbers.get(name, ""), name, library.materials[name]))
    # 番号のあるものを番号順に、無いものは後ろに名前順で
    rows.sort(key=lambda r: (r[0] == "", int(r[0]) if r[0].isdigit() else 0, r[1]))
    return rows


def _number_cell(key):
    """材料番号のセル。数字なら**数値として**入れる（Excel で並べ替えられるように）。

    空なら None（空文字を入れると Excel 上で「見えない何か」が残るため）。
    """
    text = str(key or "").strip()
    if not text:
        return None
    return int(text) if text.isdigit() else text


def _write_csv(file_name, project, model, library, assignment):
    """昔の形（CSV）で更新する。**吸音率の列は書かない**（ミスリードを避ける）。"""
    with open(file_name, "w", encoding="utf-8-sig", newline="") as f:
        f.write("# 条件表 — レイヤー名と材料番号の対応\n")
        f.write(f"# 対象室: {project.room_label}\n")
        f.write(f"# モデル: {project.dxf}\n")
        f.write("# ★書き換えるのは「材料番号」の列だけ。"
                "空欄ならレイヤー名（または先頭の番号）で引きます\n")
        f.write("# 面数・面積は参考です（計算のたびに書き直されます）\n")
        writer = csv.writer(f)
        writer.writerow(HEADER)
        for section, layer, key, count, area in _rows_for(project, model, library,
                                                          assignment):
            writer.writerow([section, layer, key, material_name(key, library),
                             count, "%.3f" % area])
    return file_name


def main():
    import argparse

    p = argparse.ArgumentParser(description="条件表（xlsx）を作る／更新する")
    p.add_argument("folder", help="プロジェクトフォルダ")
    p.add_argument("--sheet", default=FIRST_SHEET, help="作る条件シートの名前")
    a = p.parse_args()

    import project as pj
    import read_dxffile as rd

    project = pj.Project.load(a.folder)
    library = (library_from_book(path(project), kind=project.absorption_kind)
               or (ab.MaterialLibrary.from_csv(project.absorption_path,
                                               kind=project.absorption_kind)
                   if project.absorption_path else None))
    assignment = assignment_for(project)
    table = (library.absorption_table(assignment, band_number=project.band_number)
             if library else None)
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          absorption_table=table,
                          band_number=project.band_number, verbose=False)
    print(create(project, model, library, sheet=a.sheet))


if __name__ == "__main__":
    main()
