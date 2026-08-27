# -*- coding: utf-8 -*-
"""**断面の決め方**（音圧分布・モード形状を切る面）。

★2026-08-27 ユーザーと相談して決めた（2026-08-26 の「断面の設定方法は後で
相談させてください」の続き）。

> 1 面で良い場合もあれば、複数面見たい場合もありますし、その複数面は
> XY 平面であったり、XZ、YZ 平面、斜め方向など、いろいろあり得ます。

決めごと:
  ・断面の並びは**プロジェクト直下の `断面.json`** に書く（1 枚 = 1 行）。
    ★**入力**なので対象室名・条件名の頭は付けない（`視点.json` と同じ扱い）。
    条件にも受音点にも依らない
  ・書き方は 2 通りだけ
      - **軸に平行**  … `{"plane": "z", "at": 1.2}`
      - **斜め（測線）** … `{"through": [[x1,y1],[x2,y2]]}`
        ＝**平面図に引いた線を含む鉛直面**。実務で欲しい「斜め」はほぼこれ
        （ホールの中心線に沿った縦断面など）。
        ★傾いた面（点＋法線）は入れない——ユーザー判断
  ・面の中の座標は **u（横軸）・v（縦軸）** で持つ。
    軸に平行なときは**従来どおり世界座標のまま**になるように基底を選ぶ
    （既存の図・CSV と数字が変わらない）
"""
import io
import json
import math
import os

import numpy as np

FILENAME = "断面.json"
FORMAT = 1

AXIS = {"x": 0, "y": 1, "z": 2}
AXIS_NAME = ("X", "Y", "Z")

KIND_AXIS = "axis"          # 軸に平行な面
KIND_VERTICAL = "vertical"  # 平面図の測線を含む鉛直面


class Section:
    """1 枚の断面。**基準点 + 面内の 2 本の基底（u, v）**で持つ。

    `u` が図の横軸、`v` が図の縦軸。法線は `u × v`。
    """

    def __init__(self, point, u, v, name=None, kind=KIND_AXIS,
                 axis=None, value=None, through=None):
        self.point = np.asarray(point, dtype=float)
        self.u = _unit(u)
        self.v = _unit(v)
        self.normal = _unit(np.cross(self.u, self.v))
        self.kind = kind
        self.axis = axis              # 軸に平行なとき 0/1/2
        self.value = value            # 同上：その座標 [m]
        self.through = through        # 鉛直面のとき [[x1,y1],[x2,y2]]
        self.name = name or self.default_name()

    # ---- 呼び名 --------------------------------------------------------

    def default_name(self):
        if self.kind == KIND_AXIS:
            return f"{AXIS_NAME[self.axis]}={self.value:.2f}m"
        (x1, y1), (x2, y2) = self.through
        return f"測線({x1:.1f},{y1:.1f})-({x2:.1f},{y2:.1f})"

    def label(self):
        """図の題に出す説明。"""
        if self.kind == KIND_AXIS:
            return f"{AXIS_NAME[self.axis]} = {self.value:.2f} m 断面"
        (x1, y1), (x2, y2) = self.through
        return (f"測線 ({x1:.2f}, {y1:.2f}) → ({x2:.2f}, {y2:.2f}) "
                f"の鉛直断面")

    def axis_labels(self):
        """(横軸のラベル, 縦軸のラベル)。"""
        if self.kind == KIND_AXIS:
            others = [k for k in range(3) if k != self.axis]
            return f"{AXIS_NAME[others[0]]} [m]", f"{AXIS_NAME[others[1]]} [m]"
        return "測線に沿った距離 [m]", "Z [m]"

    def slug(self):
        """ファイル名に使える短い名前。"""
        text = str(self.name)
        for bad in '\\/:*?"<>|,\r\n\t':
            text = text.replace(bad, "_")
        return text.strip() or "断面"

    # ---- 座標 ----------------------------------------------------------

    def coordinates(self, points):
        """世界座標 (n,3) → 面内の (a, b)。"""
        offset = np.atleast_2d(np.asarray(points, dtype=float)) - self.point
        return offset @ self.u, offset @ self.v

    def place(self, a, b):
        """面内の (a, b) → 世界座標。"""
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        return (self.point + a[..., None] * self.u + b[..., None] * self.v)

    def height(self, points):
        """面からの符号付き距離（交線を拾うときに使う）。"""
        offset = np.atleast_2d(np.asarray(points, dtype=float)) - self.point
        return offset @ self.normal

    # ---- 出し入れ ------------------------------------------------------

    def to_dict(self):
        if self.kind == KIND_AXIS:
            return {"name": self.name, "plane": "xyz"[self.axis],
                    "at": round(float(self.value), 6)}
        return {"name": self.name,
                "through": [[round(float(v), 6) for v in pair]
                            for pair in self.through]}


def _unit(vector):
    vector = np.asarray(vector, dtype=float)
    length = float(np.linalg.norm(vector))
    if length <= 0.0:
        raise ValueError("向きが決まりません（長さ 0 のベクトル）")
    return vector / length


# ---- 作り方 ----------------------------------------------------------------

def axis_section(plane="z", value=0.0, name=None):
    """軸に平行な断面。★面内の座標が**世界座標のまま**になる基底を選ぶ。"""
    if plane not in AXIS:
        raise ValueError(f"断面の向きは x/y/z のどれか（{plane}）")
    index = AXIS[plane]
    others = [k for k in range(3) if k != index]
    point = np.zeros(3)
    point[index] = float(value)
    u, v = np.zeros(3), np.zeros(3)
    u[others[0]] = 1.0
    v[others[1]] = 1.0
    return Section(point, u, v, name=name, kind=KIND_AXIS,
                   axis=index, value=float(value))


def vertical_section(first, second, name=None):
    """**平面図に引いた測線を含む鉛直面**。`first`/`second` は (x, y)。

    横軸は測線に沿った距離（`first` が 0）、縦軸は高さ Z。
    """
    first = [float(first[0]), float(first[1])]
    second = [float(second[0]), float(second[1])]
    direction = np.array([second[0] - first[0], second[1] - first[1], 0.0])
    if float(np.linalg.norm(direction)) <= 1.0e-9:
        raise ValueError("測線の 2 点が同じ位置です")
    return Section([first[0], first[1], 0.0], direction, [0.0, 0.0, 1.0],
                   name=name, kind=KIND_VERTICAL, through=[first, second])


def from_dict(entry):
    """`断面.json` の 1 行 → Section。"""
    name = entry.get("name") or None
    if "through" in entry:
        pair = entry["through"]
        if len(pair) != 2:
            raise ValueError("through は 2 点で書いてください")
        return vertical_section(pair[0], pair[1], name=name)
    plane = str(entry.get("plane", "z")).lower()
    value = entry.get("at")
    if value is None:
        raise ValueError(f"断面『{name or plane}』に at（位置）がありません")
    return axis_section(plane, float(value), name=name)


# ---- ファイル --------------------------------------------------------------

def default_path(project):
    return os.path.join(project.folder, FILENAME)


def load(project, verbose=False):
    """`断面.json` を読む。無ければ空の並び（黙って落とさない）。"""
    path = default_path(project)
    if not os.path.exists(path):
        if verbose:
            print(f"[断面] {FILENAME} がありません（`--plane`/`--at` で 1 枚だけ"
                  "指定するか、`python section.py <プロジェクト> --create` で"
                  "雛形を作れます）")
        return []
    with io.open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    entries = document.get("sections", []) if isinstance(document, dict) \
        else list(document)
    sections = []
    for index, entry in enumerate(entries):
        try:
            sections.append(from_dict(entry))
        except Exception as error:
            print(f"[断面] {index + 1} 枚目は読めません"
                  f"（{type(error).__name__}: {error}）。飛ばします")
    if verbose and sections:
        print(f"[断面] {FILENAME} から {len(sections)} 枚読みました"
              "（" + "・".join(s.name for s in sections) + "）")
    return sections


def save(project, sections):
    """`断面.json` に書く。"""
    path = default_path(project)
    # ★**1 枚 = 1 行**で書く（人が手で直すファイルなので、
    #   json.dumps(indent=2) だと測線の座標が縦に散らばって読みにくい）
    body = ",\n".join("    " + json.dumps(s.to_dict(), ensure_ascii=False)
                      for s in sections)
    with io.open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("{\n")
        handle.write('  "format": %d,\n' % FORMAT)
        handle.write('  "sections": [\n' + body + "\n  ]\n")
        handle.write("}\n")
    return path


def suggest(model, source=None):
    """モデルから**雛形**を作る（音源を通る 3 面 ＋ 床上 1.2 m）。"""
    low = np.asarray(model.extents[0], dtype=float)
    high = np.asarray(model.extents[1], dtype=float)
    middle = 0.5 * (low + high)
    if source is None:
        source = model.source_points[0] if model.source_points else middle
    source = np.asarray(source, dtype=float)

    sections = [axis_section("z", float(low[2]) + 1.2, name="床上1.2m"),
                axis_section("x", float(source[0]), name="音源を通るYZ断面"),
                axis_section("y", float(source[1]), name="音源を通るXZ断面"),
                vertical_section([float(low[0]), float(low[1])],
                                 [float(high[0]), float(high[1])],
                                 name="対角の鉛直断面")]
    return sections


def describe(sections):
    return "\n".join(f"  {index + 1}. {s.name}　（{s.label()}）"
                     for index, s in enumerate(sections))


def main(argv=None):
    import argparse

    import project as pj

    parser = argparse.ArgumentParser(description="断面の並びを見る・作る")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--create", action="store_true",
                        help="モデルから雛形の 断面.json を作る")
    args = parser.parse_args(argv)

    project = pj.Project.load(args.folder)
    if args.create:
        import read_dxffile as rd
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              band_number=project.band_number, verbose=False)
        sections = suggest(model, project.source)
        path = save(project, sections)
        print(f"[断面] 雛形を作りました: {path}")
    else:
        sections = load(project)
        if not sections:
            print(f"[断面] {FILENAME} がありません（--create で作れます）")
            return 0
    print(f"[断面] {len(sections)} 枚")
    print(describe(sections))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
