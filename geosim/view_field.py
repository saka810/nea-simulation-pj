# -*- coding: utf-8 -*-
"""**音圧分布（モード形状）を 3D の画面で見る**（音線・音粒子・虚音源と同じ Tab）。

★2026-08-26 ユーザー要望

> 出力は、切った断面を、境界面含め表示させた結果を図として出力してほしいのと、
> 音線とかの確認画面のタブ切替で、見れるようにして欲しい。

★2026-08-27 相談で決めた（断面の設定方法）

> 複数の断面はチェックで同時に出す。

決めごと:
  ・**計算はやり直さない**。`mode_shape.py` が書いた
    `結果/<室>_<条件>_音圧分布_<番号>_<断面名>.csv` を読んで面に貼るだけ
  ・断面は**何枚でも**読む（ファイルの数だけ）。★**チェックで複数同時**に出せる
  ・★平面の決め方は**ファイル名ではなく中身の `# section:` 行**から読む。
    斜め（測線の鉛直面）はファイル名で表せないため。
    昔のファイル（`…音圧分布_z1.20m.csv`）は名前から読む
  ・周波数は CSV の列。左の欄の数値で切り替える（全断面まとめて動く）
"""
import csv
import glob
import io
import json
import os
import re

import numpy as np
import pyvista as pv

import project as pj
import section as sec

# 色の幅 [dB]（最大を 0 dB とした下限）
DEFAULT_SPAN = 30.0

BAR_TITLE = "音圧レベル [dB]"

# 昔のファイル名から断面を読む（`…_音圧分布_z1.20m.csv`）
LEGACY_NAME = re.compile(r"音圧分布_([xyz])(-?\d+(?:\.\d+)?)m\.csv$")

# 新しいファイル名（`…_音圧分布_01_床上1.2m.csv`）
NAME = re.compile(r"音圧分布_(\d+)_(.+)\.csv$")


def load_sections(project, verbose=True):
    """`結果/` に置かれた音圧分布の CSV を読む。→ [断面の辞書, …]"""
    folder = os.path.join(project.folder, pj.RESULT_DIR)
    found = []
    for path in sorted(glob.glob(os.path.join(folder, "*音圧分布_*.csv"))):
        try:
            item = read_section(path)
        except Exception as error:
            print(f"[音圧分布] {os.path.basename(path)} は読めません"
                  f"（{type(error).__name__}: {error}）")
            continue
        if item is not None:
            found.append(item)
    # ★昔の名前（`…音圧分布_z1.20m.csv`）は、同じ面の新しいファイルがあれば落とす。
    #   作り直すと名前が変わるので、消さない限り**二重に並んでしまう**
    fresh = {_plane_key(item["section"]) for item in found
             if not item["legacy"]}
    # ★`list.remove()` は**中身で**比べるので使わない
    #   （numpy の配列が入っていると「どちらが真か決まらない」で落ちる）
    keep = []
    for item in found:
        if item["legacy"] and _plane_key(item["section"]) in fresh:
            if verbose:
                print(f"[音圧分布] 昔の名前の {os.path.basename(item['path'])} は"
                      "同じ面の新しいファイルがあるので飛ばしました")
            continue
        keep.append(item)
    found = keep
    if verbose:
        if found:
            print(f"[音圧分布] 断面 {len(found)} 枚を読みました"
                  + "（" + "・".join(item["name"] for item in found) + "）")
        else:
            print("[音圧分布] 断面の CSV がありません"
                  "（`python mode_shape.py …` で作れます）")
    return found


def _plane_key(section):
    """同じ面かどうかを見るための鍵。"""
    return (round(float(section.normal[0]), 6), round(float(section.normal[1]), 6),
            round(float(section.normal[2]), 6),
            round(float(np.dot(section.point, section.normal)), 6))


def _section_of(path, comments):
    """CSV の見出しかファイル名から Section を組み立てる。"""
    for line in comments:
        if line.startswith("# section:"):
            return sec.from_dict(json.loads(line[len("# section:"):].strip()))
    match = LEGACY_NAME.search(os.path.basename(path))
    if match is not None:            # 昔のファイル（軸に平行な面だけ）
        return sec.axis_section(match.group(1), float(match.group(2)))
    return None


def _is_legacy(path, comments):
    return (not any(line.startswith("# section:") for line in comments)
            and LEGACY_NAME.search(os.path.basename(path)) is not None)


def read_section(path):
    """1 枚ぶんの CSV を読む。→ {"points", "levels", "frequencies", …}"""
    with io.open(path, encoding="utf-8-sig", newline="") as handle:
        lines = handle.read().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    rows = [row for row in csv.reader(line for line in lines
                                      if line and not line.startswith("#"))
            if row]
    if len(rows) < 2:
        return None
    section = _section_of(path, comments)
    if section is None:
        return None

    header = rows[0]
    frequencies = []
    for name in header[3:]:
        try:
            frequencies.append(float(str(name).replace("Hz_dB", "")))
        except ValueError:
            frequencies.append(float("nan"))
    points, levels = [], []
    for row in rows[1:]:
        if len(row) < 4:
            continue
        points.append([float(v) for v in row[:3]])
        levels.append([float(v) if v not in ("", None) else np.nan
                       for v in row[3:3 + len(frequencies)]])
    if not points:
        return None

    match = NAME.search(os.path.basename(path))
    name = match.group(2) if match else section.name
    return {"points": np.asarray(points, dtype=float),
            "levels": np.asarray(levels, dtype=float),
            "frequencies": np.asarray(frequencies, dtype=float),
            "section": section, "name": name, "path": path,
            "legacy": _is_legacy(path, comments),
            # 後ろ向きの持ち物（昔の呼び出しが見ていた）
            "plane": "xyz"[section.axis] if section.axis is not None else "u",
            "value": float(section.value) if section.value is not None
            else float("nan")}


def _grid_of(item):
    """点の並びから、面に貼れる格子（StructuredGrid）を作る。

    ★**断面の面内の座標（u, v）で並べ直す**ので、斜めの鉛直面でも同じ。
    """
    points = item["points"]
    section = item["section"]
    along, up = section.coordinates(points)
    first = np.unique(np.round(along, 6))
    second = np.unique(np.round(up, 6))
    shape = (len(second), len(first))
    order = np.lexsort((np.round(along, 6), np.round(up, 6)))

    sorted_points = points[order]
    if len(sorted_points) != shape[0] * shape[1]:
        return None, None, None         # 格子が欠けている（塗れない）
    grid = pv.StructuredGrid(
        sorted_points[:, 0].reshape(shape)[:, :, None],
        sorted_points[:, 1].reshape(shape)[:, :, None],
        sorted_points[:, 2].reshape(shape)[:, :, None])
    return grid, order, shape


class FieldDisplay:
    """断面の音圧分布を面に貼って見せる。★**複数同時**に出せる。"""

    def __init__(self, plotter, sections, span=DEFAULT_SPAN, label=None):
        self.plotter = plotter
        self.sections = list(sections)
        self.span = float(span)
        self.label = label
        self.band = 0
        self.actors = []
        # ★どの断面を出すか（既定は 1 枚目だけ。重ねると奥が隠れるため）
        self.shown = {0} if self.sections else set()
        self.visible = False
        self._build()

    # ---- 中身 ----------------------------------------------------------

    def _build(self):
        for item in self.sections:
            grid, order, shape = _grid_of(item)
            item["_grid"], item["_order"] = grid, order
            if grid is None:
                print(f"[音圧分布] 断面『{item['name']}』は格子が欠けているので"
                      "画面には貼れません（図と CSV は使えます）")
                self.actors.append(None)
                continue
            grid[BAR_TITLE] = self._values(item)
            actor = self.plotter.add_mesh(
                grid, scalars=BAR_TITLE, cmap="turbo",
                clim=(-self.span, 0.0), show_scalar_bar=False,
                nan_opacity=0.0, opacity=0.92, lighting=False,
                reset_camera=False, name=f"field{len(self.actors)}")
            actor.SetVisibility(False)
            self.actors.append(actor)
        self._refresh_label()

    def _values(self, item):
        band = min(self.band, item["levels"].shape[1] - 1)
        return np.asarray(item["levels"][item["_order"], band], dtype=float)

    # ---- 操作 ----------------------------------------------------------

    def set_visible(self, flag, render=True):
        """タブごと出す／隠す。"""
        self.visible = bool(flag)
        self._apply(render=render)

    def set_section_visible(self, index, flag, render=True):
        """★チェック 1 つぶん（複数同時に出せる）。"""
        index = int(index)
        if flag:
            self.shown.add(index)
        else:
            self.shown.discard(index)
        self._apply(render=render)
        return index in self.shown

    def set_section(self, index, render=True):
        """1 枚だけにする（昔の呼び方・画角の読み込みから使う）。"""
        self.shown = {int(index) % max(1, len(self.sections))}
        self._apply(render=render)
        return min(self.shown) if self.shown else 0

    def _apply(self, render=True):
        for index, actor in enumerate(self.actors):
            if actor is None:
                continue
            actor.SetVisibility(self.visible and index in self.shown)
        bar = self._scalar_bar()
        if bar is not None:
            bar.SetVisibility(self.visible and bool(self.shown))
        self._refresh_label()
        if render:
            self.plotter.render()

    def _scalar_bar(self):
        try:
            return self.plotter.scalar_bars[BAR_TITLE]
        except (KeyError, AttributeError):
            return None

    def set_band(self, index, render=True):
        self.band = max(0, int(index))
        for item in self.sections:
            if item.get("_grid") is not None:
                item["_grid"][BAR_TITLE] = self._values(item)
        self._refresh_label()
        if render:
            self.plotter.render()
        return self.band

    def frequencies(self):
        if not self.sections:
            return np.zeros(0)
        return self.sections[0]["frequencies"]

    def _refresh_label(self):
        if self.label is None or not self.sections:
            return
        values = self.frequencies()
        band = min(self.band, len(values) - 1)
        names = [self.sections[i]["name"] for i in sorted(self.shown)
                 if i < len(self.sections)]
        import view_rays as vr
        vr.ParticleAnimation._set_text(
            self.label,
            f"断面 {len(names)}/{len(self.sections)} 枚"
            + ("　" + "・".join(names[:2])
               + ("…" if len(names) > 2 else "") if names else "（なし）")
            + f"\n周波数 {values[band]:.1f} Hz（{len(values)} 本）\n"
            f"色は最大を 0 dB とした相対値")


def add_controls(panel, display, font=None):
    """左の欄に「断面のチェック」と「周波数」を並べる。→ 操作の一覧に足す行"""
    if not display.sections:
        return []
    panel.heading("音圧分布")
    for index, item in enumerate(display.sections):
        panel.checkbox(item["name"], index in display.shown,
                       lambda flag, index=index:
                       display.set_section_visible(index, flag))
    count = len(display.frequencies())
    if count > 1:
        panel.slider("周波数の番号", [1, count], 1,
                     lambda v: display.set_band(int(round(v)) - 1), fmt="%.0f")
    return ["音圧分布はチェックで断面を選び、周波数は左の欄の数値で切り替える"]


def main(argv=None):
    """中身の確認（`python view_field.py <プロジェクト>`）。"""
    import argparse

    parser = argparse.ArgumentParser(description="音圧分布の断面を確かめる")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    args = parser.parse_args(argv)

    project = pj.Project.load(args.folder)
    found = load_sections(project)
    for item in found:
        section = item["section"]
        print(f"  {item['name']}　（{section.label()}）"
              f"　{len(item['points'])} 点 / 周波数 "
              + "・".join(f"{v:.1f}" for v in item["frequencies"]) + " Hz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
