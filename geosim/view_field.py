# -*- coding: utf-8 -*-
"""**音圧分布（モード形状）を 3D の画面で見る**（音線・音粒子・虚音源と同じ Tab）。

★2026-08-26 ユーザー要望

> 出力は、切った断面を、境界面含め表示させた結果を図として出力してほしいのと、
> 音線とかの確認画面のタブ切替で、見れるようにして欲しい。

決めごと:
  ・**計算はやり直さない**。`mode_shape.py` が書いた
    `結果/<室>_<条件>_音圧分布_<断面>.csv` を読んで面に貼るだけ
  ・断面は**何枚でも**読む（ファイルの数だけ）。切り替えは左の欄の数値で
  ・周波数も**何本でも**（CSV の列）。同じく数値で切り替える
  ・断面の位置・向きの決め方（1 枚／複数枚／斜め）は**これから相談**。
    いまは「作った CSV を並べて見せる」ところまで
"""
import csv
import glob
import io
import os
import re

import numpy as np
import pyvista as pv

import project as pj

# 色の幅 [dB]（最大を 0 dB とした下限）
DEFAULT_SPAN = 30.0

BAR_TITLE = "音圧レベル [dB]"

# ファイル名から断面を読む（`…_音圧分布_z1.20m.csv`）
NAME = re.compile(r"音圧分布_([xyz])(-?\d+(?:\.\d+)?)m\.csv$")


def load_sections(project, verbose=True):
    """`結果/` に置かれた音圧分布の CSV を読む。→ [断面の辞書, …]"""
    folder = os.path.join(project.folder, pj.RESULT_DIR)
    found = []
    for path in sorted(glob.glob(os.path.join(folder, "*音圧分布_*.csv"))):
        match = NAME.search(os.path.basename(path))
        if match is None:
            continue
        section = read_section(path)
        if section is not None:
            section["plane"] = match.group(1)
            section["value"] = float(match.group(2))
            section["path"] = path
            found.append(section)
    if verbose:
        if found:
            print(f"[音圧分布] 断面 {len(found)} 枚を読みました"
                  + "（" + "・".join(f"{s['plane']}={s['value']:.2f}m"
                                     for s in found) + "）")
        else:
            print("[音圧分布] 断面の CSV がありません"
                  "（`python mode_shape.py …` で作れます）")
    return found


def read_section(path):
    """1 枚ぶんの CSV を読む。→ {"points", "levels", "frequencies", …}"""
    with io.open(path, encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.reader(handle)
                if row and not row[0].startswith("#")]
    if len(rows) < 2:
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
    return {"points": np.asarray(points, dtype=float),
            "levels": np.asarray(levels, dtype=float),
            "frequencies": np.asarray(frequencies, dtype=float)}


def _grid_of(section):
    """点の並びから、面に貼れる格子（StructuredGrid）を作る。"""
    points = section["points"]
    axis = {"x": 0, "y": 1, "z": 2}[section["plane"]]
    others = [k for k in range(3) if k != axis]
    first = np.unique(np.round(points[:, others[0]], 6))
    second = np.unique(np.round(points[:, others[1]], 6))
    shape = (len(second), len(first))

    order = np.lexsort((np.round(points[:, others[0]], 6),
                        np.round(points[:, others[1]], 6)))
    coordinates = [None, None, None]
    mesh_a, mesh_b = np.meshgrid(first, second)
    coordinates[others[0]] = mesh_a
    coordinates[others[1]] = mesh_b
    coordinates[axis] = np.full(shape, float(section["value"]))
    grid = pv.StructuredGrid(coordinates[0][:, :, None],
                             coordinates[1][:, :, None],
                             coordinates[2][:, :, None])
    return grid, order, shape


class FieldDisplay:
    """断面の音圧分布を面に貼って見せる。"""

    def __init__(self, plotter, sections, span=DEFAULT_SPAN, label=None):
        self.plotter = plotter
        self.sections = list(sections)
        self.span = float(span)
        self.label = label
        self.section = 0
        self.band = 0
        self.actors = []
        self.visible = False
        self._build()

    # ---- 中身 ----------------------------------------------------------

    def _build(self):
        for section in self.sections:
            grid, order, shape = _grid_of(section)
            section["_grid"], section["_order"] = grid, order
            actor = self.plotter.add_mesh(
                grid, scalars=self._values(section), cmap="turbo",
                clim=(-self.span, 0.0), show_scalar_bar=False,
                nan_opacity=0.0, opacity=0.92, lighting=False,
                reset_camera=False, name=f"field{len(self.actors)}")
            actor.SetVisibility(False)
            self.actors.append(actor)
        self._refresh_label()

    def _values(self, section):
        band = min(self.band, section["levels"].shape[1] - 1)
        values = section["levels"][section["_order"], band]
        return np.asarray(values, dtype=float)

    # ---- 操作 ----------------------------------------------------------

    def set_visible(self, flag, render=True):
        self.visible = bool(flag)
        for index, actor in enumerate(self.actors):
            actor.SetVisibility(self.visible and index == self.section)
        bar = self._scalar_bar()
        if bar is not None:
            bar.SetVisibility(self.visible)
        if render:
            self.plotter.render()

    def _scalar_bar(self):
        try:
            return self.plotter.scalar_bars[BAR_TITLE]
        except (KeyError, AttributeError):
            return None

    def set_section(self, index, render=True):
        self.section = int(index) % max(1, len(self.sections))
        self.set_visible(self.visible, render=False)
        self._refresh_label()
        if render:
            self.plotter.render()
        return self.section

    def set_band(self, index, render=True):
        self.band = max(0, int(index))
        for section in self.sections:
            grid = section["_grid"]
            grid[BAR_TITLE] = self._values(section)
        self._refresh_label()
        if render:
            self.plotter.render()
        return self.band

    def frequencies(self):
        if not self.sections:
            return np.zeros(0)
        return self.sections[min(self.section, len(self.sections) - 1)]["frequencies"]

    def _refresh_label(self):
        if self.label is None or not self.sections:
            return
        section = self.sections[min(self.section, len(self.sections) - 1)]
        band = min(self.band, len(section["frequencies"]) - 1)
        import view_rays as vr
        vr.ParticleAnimation._set_text(
            self.label,
            f"断面 {self.section + 1}/{len(self.sections)}"
            f"　{section['plane'].upper()} = {section['value']:.2f} m\n"
            f"周波数 {section['frequencies'][band]:.1f} Hz"
            f"（{len(section['frequencies'])} 本）\n"
            f"色は最大を 0 dB とした相対値")


def add_controls(panel, display, font=None):
    """左の欄に「断面」「周波数」の数値を並べる。→ 操作の一覧に足す行"""
    if not display.sections:
        return []
    panel.heading("音圧分布")
    if len(display.sections) > 1:
        panel.slider("断面", [1, len(display.sections)], 1,
                     lambda v: display.set_section(int(round(v)) - 1),
                     fmt="%.0f")
    count = len(display.frequencies())
    if count > 1:
        panel.slider("周波数の番号", [1, count], 1,
                     lambda v: display.set_band(int(round(v)) - 1), fmt="%.0f")
    return ["音圧分布は左の欄で断面と周波数を選ぶ（`mode_shape.py` で作った断面）"]
