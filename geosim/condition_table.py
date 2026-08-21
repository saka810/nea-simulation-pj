"""材料条件表 ― 「レイヤー名 → 吸音材」の対応を CSV で持つ（依頼 2026-08-21）。

> いまレイヤー名から番号判別してもらってますが、壁の吸音材条件を変える場合に、
> CAD データを変えて、それをまた読み込んでもらって、っていうのが大変な気がします。
> なので、レイヤー名と壁材料の一覧表（条件表）を作って、それを読み込ませる形式に
> 変えたいと思います。そうすれば、ソフト上で吸音材かえた場合に、その記録を
> 反映しておけば、見返しやすいとも思いました。

## 何が変わるか

これまでは **CAD のレイヤ名そのもの**（または先頭の番号）で吸音材を引いていた。
材料を変えるには CAD のレイヤ名を変えて DXF を書き出し直す必要があった。

これからは**プロジェクトフォルダの `材料条件表.csv` が対応表**になる。
CAD は触らずに、この 1 枚を書き換えれば材料が変わる。

    プロジェクトフォルダ/材料条件表.csv

        # 材料条件表 — レイヤー名と吸音材の対応
        # 対象室・条件名: 研修室_吸音追加案
        # モデル: 研修室_faces.dxf
        # 吸音率表: absorption.csv（残響室法）
        区分,レイヤー名,材料名,面数,面積_m2,63,125,250,500,1000,2000,4000,8000
        レイヤ,01__研修室_壁_扉,01_扉(木製),4,3.255,0.15,0.29,…
        レイヤ,09__研修室_床,09_タイルカーペット,64,123.964,0.01,0.02,…
        面ごとの指定,11_吸音板,11_吸音板,12,67.163,0.05,0.09,…

- **書き換えるのは「材料名」の列だけ。**空欄にすると従来どおり
  「レイヤー名（または先頭の番号）で引く」動きに戻る
- 面数・面積・吸音率の列は**参考として毎回書き直される**（手で直しても次で戻る）
- `#` で始まる行は読み飛ばす。条件を後から見返すための覚え書き

## 2 種類の行

| 区分 | 意味 | 書き換えられるか |
|---|---|---|
| `レイヤ` | DXF のレイヤに材料を割り当てる（**これが入力**） | ○ |
| `面ごとの指定` | `face_editor` で面を選んで貼った材料の**記録** | ×（記録） |

面ごとの割り当ては「材料名 → 面番号」で `materials.json` が持っており、
レイヤの表では表せない。**ただし見返せないと困る**ので、
どの材料が何枚・何 m² に貼られているかをこの表に書き出す。

## 更新の約束（★ここを崩さないこと）

**利用者が書いた「材料名」は絶対に上書きしない。**
`update()` は既存の表を読み、材料名の列はそのまま残して、
面数・面積・吸音率だけを書き直す。新しいレイヤの行だけを足し、
モデルから消えたレイヤは区分を `（モデルに無し）` にして残す
（消すと「前はどう設定していたか」が分からなくなる）。
"""

import csv
import os

import numpy as np

import absorption as ab

CONDITION_FILE = "材料条件表.csv"

# 区分の名前
SECTION_LAYER = "レイヤ"
SECTION_FACE = "面ごとの指定"
SECTION_GONE = "（モデルに無し）"

HEADER = ["区分", "レイヤー名", "材料名", "面数", "面積_m2"]


def path(project):
    """材料条件表のパス。**プロジェクトフォルダ直下**（`normals.json` と同じ扱い）。

    対象室・条件名の頭は付けない。これは結果ではなく**入力**で、
    プロジェクトを名前変更してもそのまま使えるほうがよいため。
    """
    return project.path(CONDITION_FILE)


def exists(project):
    return os.path.exists(path(project))


# ------------------------------------------------------------------------------
# 読む
# ------------------------------------------------------------------------------

def read(file_name):
    """材料条件表を読む。戻り値 (割り当て, 記録).

    割り当て … {レイヤー名: 材料名}（材料名が空の行は入れない）
    記録     … [(区分, レイヤー名, 材料名, 面数, 面積), …] そのままの並び
    """
    if not file_name or not os.path.exists(file_name):
        return {}, []
    with open(file_name, encoding="utf-8-sig", newline="") as f:
        rows = [row for row in csv.reader(f)
                if row and row[0].strip() and not row[0].startswith("#")]
    if not rows:
        return {}, []
    if rows[0][0].strip() == HEADER[0]:
        rows = rows[1:]

    assignment, records = {}, []
    for row in rows:
        row = row + [""] * (5 - len(row))
        section, layer, material = (row[0].strip(), row[1].strip(), row[2].strip())
        records.append((section, layer, material, row[3].strip(), row[4].strip()))
        if section == SECTION_LAYER and layer and material:
            assignment[layer] = material
    return assignment, records


def assignment_for(project, verbose=True):
    """プロジェクトの材料条件表から {レイヤー名: 材料名} を返す。

    表が無ければ `project.assignment`（project.json の指定）にそのまま戻る。
    **表があればそちらを優先する**（利用者が触るのは表のほうなので）。
    """
    if not exists(project):
        return project.assignment
    assignment, _ = read(path(project))
    if verbose:
        print(f"[条件表] {CONDITION_FILE} から {len(assignment)} レイヤの"
              f"材料割り当てを読みました")
    if project.assignment and verbose:
        print(f"[条件表] 注意: project.json 側の割り当ては使いません"
              f"（条件表を優先します）")
    return assignment or None


# ------------------------------------------------------------------------------
# 書く（作る・更新する）
# ------------------------------------------------------------------------------

def resolve_material(layer, library, assignment):
    """そのレイヤに実際に使われる材料名を返す。引けなければ ""。

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
        return library.aliases.get(layer, layer)
    number = _layer_number(layer)
    if number is not None and number in library:
        return library.aliases.get(number, number)
    return ""


def _layer_number(layer):
    import read_dxffile as rd
    return rd.layer_number(layer)


def update(project, model, library=None, assignment=None, verbose=True):
    """材料条件表を作る／更新する。書いたパスを返す。

    ★**利用者が書いた「材料名」は上書きしない。**
    面数・面積・吸音率の参考列だけ書き直し、新しいレイヤの行を足す。

    引数:
        model    `read_dxffile.read_model()` の結果（面数・面積・レイヤを取る）
        library  `absorption.MaterialLibrary`（吸音率の参考列に使う。無くてもよい）
        assignment 既に読んである割り当て（無ければ表と project.json から）
    """
    file_name = path(project)
    previous, _ = read(file_name)
    if assignment is None:
        assignment = previous or project.assignment or {}

    frequencies = ab.octave_bands(project.band_number)
    layer_areas = getattr(model, "layer_areas", {}) or {}
    layer_counts = getattr(model, "layer_counts", {}) or {}

    rows = []
    for layer in sorted(set(layer_counts) | set(layer_areas)):
        # ★材料名は「利用者が書いたもの」が最優先。無ければ自動で引けた名前を入れる
        material = previous.get(layer) or resolve_material(layer, library, assignment)
        rows.append((SECTION_LAYER, layer, material, layer_counts.get(layer, 0),
                     layer_areas.get(layer, 0.0),
                     _coefficients(material, library, project.band_number)))

    # モデルから消えたレイヤも残す（前の設定を見返せるように）
    for layer, material in sorted(previous.items()):
        if layer not in layer_counts and layer not in layer_areas:
            rows.append((SECTION_GONE, layer, material, 0, 0.0,
                         _coefficients(material, library, project.band_number)))

    # 面ごとの指定（face_editor で貼ったもの）の記録
    for material, (count, area) in sorted(_face_records(project, model).items()):
        rows.append((SECTION_FACE, material, material, count, area,
                     _coefficients(material, library, project.band_number)))

    _write(file_name, rows, frequencies, project, model)
    if verbose:
        layers = sum(1 for r in rows if r[0] == SECTION_LAYER)
        empty = sum(1 for r in rows if r[0] == SECTION_LAYER and not r[2])
        print(f"[条件表] {file_name}（レイヤ {layers} / 材料未設定 {empty}）")
    return file_name


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


def _coefficients(material, library, band_number):
    """参考として載せる吸音率（**計算に使うのと同じ垂直入射の値**）。"""
    if not material or library is None:
        return None
    entry = library.get(material)
    if entry is None:
        return None
    resampled = ab.Material(entry.name, entry.resample(band_number),
                            entry.kind, entry.note)
    return resampled.normal_incidence(warn=False)


def _write(file_name, rows, frequencies, project, model):
    kind = {"normal": "垂直入射", "random": "残響室法"}.get(project.absorption_kind,
                                                            "種類未指定")
    with open(file_name, "w", encoding="utf-8-sig", newline="") as f:
        f.write("# 材料条件表 — レイヤー名と吸音材の対応\n")
        f.write(f"# 対象室・条件名: {project.name}\n")
        f.write(f"# モデル: {project.dxf}\n")
        f.write(f"# 吸音率表: {project.absorption_csv}（{kind}）\n")
        f.write("# ★書き換えるのは「材料名」の列だけ。"
                "空欄ならレイヤー名（または先頭の番号）で引きます\n")
        f.write("# 面数・面積・吸音率は参考です（計算のたびに書き直されます）\n")
        writer = csv.writer(f)
        writer.writerow(HEADER + [f"{v:.0f}" for v in frequencies])
        for section, layer, material, count, area, alpha in rows:
            cells = [section, layer, material, count, "%.3f" % area]
            if alpha is None:
                cells += [""] * len(frequencies)
            else:
                cells += ["%.4g" % v for v in np.asarray(alpha).ravel()]
            writer.writerow(cells)
    return file_name


def main():
    import argparse

    p = argparse.ArgumentParser(
        description="材料条件表（レイヤー名 → 吸音材）を作る／更新する")
    p.add_argument("folder", help="プロジェクトフォルダ")
    a = p.parse_args()

    import project as pj
    import read_dxffile as rd

    project = pj.Project.load(a.folder)
    library = (ab.MaterialLibrary.from_csv(project.absorption_path,
                                           kind=project.absorption_kind)
               if project.absorption_path else None)
    assignment = assignment_for(project)
    table = (library.absorption_table(assignment, band_number=project.band_number)
             if library else None)
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          absorption_table=table,
                          band_number=project.band_number, verbose=False)
    print(update(project, model, library, assignment))


if __name__ == "__main__":
    main()
