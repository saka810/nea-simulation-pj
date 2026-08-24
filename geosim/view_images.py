# -*- coding: utf-8 -*-
"""虚音源の可視化（G-31。2026-08-24 ユーザー要望）。

> 虚音源を可視化して見れるようにしたい。
> ・虚音源から指定した受音点（画面上で指定）を結ぶ線を描く。
> ・音線と同じように、N 回反射まで表示。N 回目反射のみを表示というのも欲しいので、
>   反射回数の開始と終了を設定できるようにすれば良いかと思います。

## 何を描いているか

虚像法（虚音源法）は、壁での反射を「壁の向こう側に音源を鏡像として置く」ことに
置き換える方法です（書籍 図2.18、`docs/技術説明書.md` 7.1 節）。
虚音源と受音点を結ぶ**直線**が、折れ曲がった実際の経路を**伸ばしたもの**で、
その長さがそのまま経路長（＝ 音速 × 到来時刻）になります。

    実際の経路            伸ばした経路（この画面で描くもの）
    音源 ─┐               虚音源 ────────── 受音点
          └─ 壁 ─→ 受音点        （直線。壁を貫いて見える）

**壁を貫いた線に見えるのが正しい姿**です。線が壁と交わる点が反射点にあたります。

## データはどこから来るか

★**計算をやり直しません。**パルス列（`結果/recN/<室>_pulses.csv`）から復元します。

パルス列は到来方向（受音点から虚音源を見込む単位ベクトル）と距離を持っているので、

    虚音源の位置 = 受音点 + 到来方向 × 距離

で虚音源そのものが出ます（`loop_noredundancy.PulseList.image_sources`）。
直接音（反射 0 回）の虚音源は音源そのものに一致します（検算に使える）。

## 操作

    ドラッグ 回転 / ホイール 拡大縮小 / `z` `x` `c` `v` 視点
    **`w` 受音点を切り替え**（左パネルのスライダでも選べる）
    **反射回数の「開始」「終了」スライダ**（同じ値にすれば「N 回目のみ」）
    `本数` スライダ … 多いと線で埋まるので、**エネルギーの大きい順**に絞る
    `m` 色（反射回数 / エネルギー / 到来時刻）   `l` 線の表示 ON/OFF
    `h` 室に寄る   `r` 全部が入るように引く（VTK のキー）
    `o` 不透明度の対象を切り替え   `t` スライダの値を数字で入れる
    `g` いまの画面を画像で保存   `F1` 操作の一覧   `q` 閉じる

## 資料用に図だけ書き出す

    python view_images.py <プロジェクト> --start 1 --end 1 --receiver 1                           --fit room --screenshot 1次反射.png

`--receiver`（1 始まり）と `--fit`（`all` = 虚音源まで全部入る / `room` = 室に寄る）は
画面の `w` と `h` に相当する。
"""

import os
import sys

import numpy as np

try:
    import pyvista as pv
except ImportError as error:      # pragma: no cover
    raise SystemExit("pyvista が要ります。`pip install -r requirements.txt`") from error

import read_dxffile as rd
import view_model_gui as vg

# 起動時の見え方。**反射回数を絞った状態で開く**（高次の虚音源は室から遠く離れた
# ところに並ぶので、全部描くと室が点になって何も読めない）
DEFAULT_ORDER_START = 0
DEFAULT_ORDER_END = 2
DEFAULT_COUNT = 40
POOL_SIZE = 4000            # 抱えておく候補の上限（スライダの上限になる）

COLOUR_MODES = ("reflection", "energy", "time")
COLOUR_LABEL = {"reflection": "反射回数", "energy": "エネルギー [dB]",
                "time": "到来時刻 [ms]"}
COLOUR_BAR_TITLE = {"reflection": "Reflection order", "energy": "Energy [dB]",
                    "time": "Arrival time [ms]"}
COLOUR_CMAP = {"reflection": "viridis", "energy": "plasma", "time": "cividis"}

IMAGE_COLOR = "#ffd166"         # 虚音源の点（色分けしないとき）
RECEIVER_COLOR = "#4cc9f0"
SOURCE_COLOR = "#ff5f5f"
DIM_COLOR = "#5a6270"           # 選んでいない受音点
OUTLINE_COLOR = "#8b93a3"       # 室の形（同一平面パッチの外周）

# ★壁の不透明度（2026-08-24 ユーザー指摘「壁が見えなくなっている」）。
#   音線の画面は 0.10 だが、この画面は虚音源との位置関係を見るものなので
#   **壁がはっきり見えるほうがよい**。あわせて外周の線も重ねる
DEFAULT_OPACITY = 0.30


# ------------------------------------------------------------------------------
# 虚音源の集まり（純計算。画面に依らない）
# ------------------------------------------------------------------------------

class ImageSourceSet:
    """1 つの受音点に対する虚音源の集まり。

    ★**画面に依らない純計算**にしてある（テストしやすいように）。
    描く側は `select()` が返す添字を使うだけ。
    """

    def __init__(self, receiver_point, pulses, name="rec"):
        self.name = name
        self.receiver = np.asarray(receiver_point, dtype=float).reshape(3)
        self.order = np.asarray(pulses.reflection_count, dtype=int)
        self.distance = np.asarray(pulses.distance, dtype=float)
        self.time = np.asarray(pulses.time, dtype=float)
        self.position = pulses.image_sources(self.receiver)
        energy = np.asarray(pulses.energy, dtype=float)
        # 帯域を合わせた値で強さを見る（どの虚音源が効いているか）
        self.energy = energy.sum(axis=1) if energy.ndim == 2 and energy.size \
            else np.ones(len(self.order))

    def __len__(self):
        return len(self.order)

    @property
    def order_limit(self):
        return int(self.order.max()) if len(self) else 0

    def select(self, start=0, end=None, count=None):
        """**反射回数の範囲**で絞り、強い順に `count` 本だけ残した添字を返す。

        ★「N 回目反射のみ」は `start == end == N` で表せる（ユーザー要望どおり
        開始と終了を別に持たせた）。範囲に何も無ければ空を返す。

        絞るときに**エネルギーの大きい順**にするのが要点。適当に間引くと
        「聞こえ方に効いている虚音源」が落ちて、弱いものばかり残ることがある。
        """
        end = self.order_limit if end is None else int(end)
        start = int(start)
        if start > end:
            start, end = end, start
        inside = np.nonzero((self.order >= start) & (self.order <= end))[0]
        if count is None or len(inside) <= count:
            return inside
        strongest = np.argsort(self.energy[inside])[::-1][:int(count)]
        # 添字の順は「並べ替えない」（時刻順のまま描くほうが見やすい）
        return np.sort(inside[strongest])

    def values(self, index, mode="reflection"):
        """色に使う値。`mode` は `COLOUR_MODES` のいずれか。"""
        index = np.asarray(index, dtype=int)
        if mode == "energy":
            return 10.0 * np.log10(np.maximum(self.energy[index], 1e-30))
        if mode == "time":
            return self.time[index] * 1000.0
        return self.order[index].astype(float)

    def line_polydata(self, index, mode="reflection"):
        """虚音源と受音点を結ぶ線（`index` の本数だけ）。"""
        index = np.asarray(index, dtype=int)
        if len(index) == 0:
            return None
        far = self.position[index]
        near = np.repeat(self.receiver.reshape(1, 3), len(index), axis=0)
        points = np.empty((2 * len(index), 3))
        points[0::2] = far
        points[1::2] = near
        lines = np.column_stack([np.full(len(index), 2),
                                 np.arange(0, 2 * len(index), 2),
                                 np.arange(1, 2 * len(index), 2)]).ravel()
        poly = pv.PolyData(points, lines=lines)
        value = self.values(index, mode)
        poly.point_data[mode] = np.repeat(value, 2)
        return poly

    def point_polydata(self, index, mode="reflection"):
        """虚音源そのもの（点）。"""
        index = np.asarray(index, dtype=int)
        if len(index) == 0:
            return None
        poly = pv.PolyData(self.position[index])
        poly.point_data[mode] = self.values(index, mode)
        return poly

    def summary(self, start=None, end=None):
        if not len(self):
            return f"{self.name}: 虚音源なし（受音した経路がありません）"
        text = (f"{self.name}: 虚音源 {len(self)} 個 / "
                f"反射 {self.order.min()}〜{self.order_limit} 回 / "
                f"経路長 {self.distance.min():.2f}〜{self.distance.max():.1f} m")
        if start is not None:
            picked = self.select(start, end)
            text += f" / いま {len(picked)} 個"
        return text


# ------------------------------------------------------------------------------
# 画面に描く側
# ------------------------------------------------------------------------------

class ImageSourceDisplay:
    """虚音源の点と、受音点を結ぶ線を描く。**あとから絞り込みを変えられる。**"""

    def __init__(self, plotter, sets, index=0, start=DEFAULT_ORDER_START,
                 end=DEFAULT_ORDER_END, count=DEFAULT_COUNT, mode="reflection",
                 line_width=1.6, point_size=11.0, opacity=0.85, label=None):
        self.plotter = plotter
        self.sets = list(sets)
        self.index = int(index)
        self.start = int(start)
        self.end = int(end)
        self.count = int(count)
        self.mode = mode
        self.line_width = line_width
        self.point_size = point_size
        self.opacity = opacity
        self.label = label
        self.show_lines = True
        self.visible = True
        self.actors = []
        self.receiver_actors = []
        self.rebuild(render=False)

    # ---- いまの状態 ---------------------------------------------------
    @property
    def current(self):
        return self.sets[self.index]

    @property
    def order_limit(self):
        return max(1, max(s.order_limit for s in self.sets))

    # ---- 操作 ---------------------------------------------------------
    def set_receiver(self, index, render=True):
        """★描く相手の受音点を変える（`w` キー／パネルのスライダ）。"""
        self.index = int(index) % len(self.sets)
        self.rebuild(render=render)

    def next_receiver(self):
        self.set_receiver(self.index + 1)

    def set_start(self, value, render=True):
        self.start = max(0, int(round(value)))
        self.rebuild(render=render)

    def set_end(self, value, render=True):
        self.end = max(0, int(round(value)))
        self.rebuild(render=render)

    def set_count(self, value, render=True):
        self.count = max(1, int(round(value)))
        self.rebuild(render=render)

    def next_mode(self):
        self.mode = COLOUR_MODES[(COLOUR_MODES.index(self.mode) + 1)
                                 % len(COLOUR_MODES)]
        self.rebuild()

    def toggle_lines(self):
        self.show_lines = not self.show_lines
        self.rebuild()

    # ---- 描き直し -----------------------------------------------------
    def rebuild(self, render=True):
        for actor in self.actors:
            self.plotter.remove_actor(actor, render=False)
        self.actors = []

        current = self.current
        index = current.select(self.start, self.end, self.count)
        bar = {"title": COLOUR_BAR_TITLE[self.mode], "color": vg.TEXT_COLOR,
               "n_labels": 5, "fmt": "%.1f", "position_x": 0.32,
               "position_y": 0.11, "width": 0.5, "height": 0.05}

        if len(index):
            # ★**値に幅が無いときは色分けしない**（2026-08-24）。
            #   「N 回目のみ」にすると反射回数は全部同じなので、カラーバーが
            #   「1.0 1.0 1.0 1.0」と並んで意味が無く、かえって紛らわしい。
            #   そのときは 1 色で塗る
            value = current.values(index, self.mode)
            flat = (len(value) < 2) or (float(value.max() - value.min()) <= 0.0)
            paint = ({"color": IMAGE_COLOR, "show_scalar_bar": False} if flat
                     else {"scalars": self.mode, "cmap": COLOUR_CMAP[self.mode]})

            if self.show_lines:
                lines = current.line_polydata(index, self.mode)
                # ★**カメラを動かさない**（`reset_camera=False`）。高次の虚音源は
                #   室から遠く離れるので、合わせに行くと室が点になって読めない
                self.actors.append(self.plotter.add_mesh(
                    lines, line_width=self.line_width, lighting=False,
                    opacity=self.opacity, reset_camera=False,
                    scalar_bar_args=None if flat else bar, **paint))
            points = current.point_polydata(index, self.mode)
            show_bar = (not flat) and (not self.show_lines)
            self.actors.append(self.plotter.add_mesh(
                points, point_size=self.point_size,
                render_points_as_spheres=True, lighting=False,
                reset_camera=False,
                **({**paint, "show_scalar_bar": show_bar,
                    "scalar_bar_args": bar if show_bar else None}
                   if not flat else paint)))

        self._draw_receivers()
        # 隠している間に作り直しても、隠れたままにする（Tab で戻したら見える）
        if not self.visible:
            for actor in self.actors + self.receiver_actors:
                actor.SetVisibility(False)
        self._update_label(len(index))
        if render:
            self.plotter.render()

    def set_visible(self, flag, render=True):
        """★**表示の ON/OFF**（音線の画面と Tab で切り替えるために要る。2026-08-24）。

        隠している間も持ち物は残しておく（作り直さずに戻せるように）。
        """
        self.visible = bool(flag)
        for actor in self.actors + self.receiver_actors:
            actor.SetVisibility(self.visible)
        for title in COLOUR_BAR_TITLE.values():
            try:
                self.plotter.scalar_bars[title].SetVisibility(self.visible)
            except (KeyError, AttributeError):
                pass
        if self.label is not None:
            self.label.SetVisibility(self.visible)
        if render:
            self.plotter.render()

    def _draw_receivers(self):
        """受音点を描く。**選んでいる点だけ濃く**して、どれが相手かを示す。

        ★大きさは変えない（2026-08-24 ユーザー指摘「受音点を大きくしないで良い」）。
        色だけで十分見分けられるので、大きさをいじると室の中で目障りになる。
        """
        for actor in self.receiver_actors:
            self.plotter.remove_actor(actor, render=False)
        self.receiver_actors = []
        span = self._span()
        for k, item in enumerate(self.sets):
            chosen = (k == self.index)
            self.receiver_actors.append(self.plotter.add_mesh(
                pv.Sphere(radius=span * 0.006, center=item.receiver),
                color=RECEIVER_COLOR if chosen else DIM_COLOR,
                lighting=False, reset_camera=False, show_scalar_bar=False))

    def _span(self):
        bounds = self.plotter.bounds
        size = max(bounds[1] - bounds[0], bounds[3] - bounds[2],
                   bounds[5] - bounds[4])
        return size if size > 0 else 1.0

    def _update_label(self, shown):
        if self.label is None:
            return
        current = self.current
        picked = current.select(self.start, self.end)
        band = (f"{self.start} 回目のみ" if self.start == self.end
                else f"{self.start}〜{self.end} 回")
        vg.set_actor_text(self.label,
                          f"{current.name}：反射 {band}\n"
                          f"該当 {len(picked)} 個 → {shown} 個を描画"
                          f"（全 {len(current)} 個）\n"
                          f"色 = {COLOUR_LABEL[self.mode]}"
                          + ("" if self.show_lines else "／線は非表示"))


# ------------------------------------------------------------------------------
# 読み込み（プロジェクトから）
# ------------------------------------------------------------------------------

def load_sets(project, verbose=True):
    """プロジェクトの受音点ごとに `ImageSourceSet` を作る。

    パルス列の CSV（`結果/recN/<室>_pulses.csv`）を読むだけなので**計算しない**。
    1 点でも読めれば、その点だけで開く（全部そろっている必要はない）。
    """
    import loop_noredundancy as ln
    import project as pj

    receivers = _receivers(project)
    if not len(receivers):
        raise ValueError("受音点が見つかりません（DXF の rec レイヤか project.json）")

    sets = []
    for k, point in enumerate(receivers):
        sub = pj.Project(project.folder,
                         **{key: getattr(project, key) for key in pj.DEFAULTS})
        sub.receiver_index = k + 1
        path = sub.existing_result_path("pulses")
        if not os.path.exists(path):
            if verbose:
                print(f"[虚音源] rec{k + 1} のパルス列がありません: {path}")
            continue
        try:
            pulses = ln.PulseList.from_csv(path)
        except Exception as error:
            print(f"[虚音源] rec{k + 1} を読めませんでした: "
                  f"{type(error).__name__}: {error}")
            continue
        item = ImageSourceSet(point, pulses, name=f"rec{k + 1}")
        sets.append(item)
        if verbose:
            print(f"[虚音源] {item.summary()}")
    if not sets:
        raise ValueError("パルス列（pulses.csv）が見つかりません。先に計算してください")
    return sets


def _receivers(project):
    """受音点の一覧（`run_project._receivers` と同じ決め方）。"""
    if project.receiver is not None:
        return [np.asarray(project.receiver, dtype=float)]
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          band_number=project.band_number, verbose=False)
    return [np.asarray(p, dtype=float) for p in model.receiver_points]


# ------------------------------------------------------------------------------
# 画面
# ------------------------------------------------------------------------------

def view(dxf_path, sets, absorption=None, unit=None, band_number=None,
         orient_normals="cad", source_points=None, opacity=DEFAULT_OPACITY,
         layer_opacity=None, start=DEFAULT_ORDER_START, end=DEFAULT_ORDER_END,
         count=DEFAULT_COUNT, mode="reflection", screenshot=None, save_dir=None,
         receiver=0, fit="all"):
    """虚音源をモデルの上に重ねて表示する。

    sets … `ImageSourceSet` のリスト（受音点ごと）
    """
    kwargs = {} if band_number is None else {"band_number": band_number}
    table = absorption
    if isinstance(absorption, str):
        import absorption as ab
        try:
            table = ab.MaterialLibrary.from_file(absorption).absorption_table(
                band_number=band_number, warn=False)
        except Exception as error:
            print(f"[虚音源] 吸音率を読めませんでした（{error}）。色は既定になります")
            table = None
    model = rd.read_model(dxf_path, unit=unit, absorption_table=table,
                          orient_normals=orient_normals, **kwargs)

    base = os.path.splitext(os.path.basename(dxf_path))[0]
    off_screen = screenshot is not None
    plotter = vg.build_plotter(model, title=f"{base} 虚音源", off_screen=off_screen,
                               show_normals=False, opacity=opacity,
                               layer_opacity=layer_opacity, show_summary=False,
                               panel=True, screen="images")
    panel = vg.control_panel(plotter)
    font = vg.japanese_font()

    # ★★**カメラと目盛りは「室」に合わせたまま保つ**（2026-08-24）。
    #   高次の虚音源は室の外へどこまでも離れていく（研修室で経路長 3.4 km）。
    #   VTK は描いてあるもの全部が入るように自動で引くので、放っておくと
    #   **室が画面の真ん中の点になって何も読めない**（実際にそうなった）。
    #   ここで室の範囲を覚えておき、虚音源を足したあとで戻す。
    room_bounds = tuple(plotter.bounds)

    # ★**室の形を線で重ねる**（2026-08-24 ユーザー指摘）。
    #   壁を薄く描くと形が読めなくなるので、同一平面パッチの外周を引く。
    #   これがあると不透明度を 0 まで下げても室の輪郭が残る
    _add_room_outline(plotter, model)

    # 音源（虚音源と見比べるため。直接音の虚音源はここに一致する）
    points = source_points if source_points is not None else model.source_points
    span = max(model.extents[1] - model.extents[0]) if model.extents is not None else 1.0
    for point in (points or []):
        plotter.add_mesh(pv.Sphere(radius=span * 0.012, center=np.asarray(point)),
                         color=SOURCE_COLOR, lighting=False, reset_camera=False,
                         show_scalar_bar=False)

    label = None
    if panel is not None:
        # ★題名は `build_plotter` がすでに出しているので、ここでは出さない
        panel.heading("いま描いているもの")
        label = panel.reserve_text(3, size=9)

    display = ImageSourceDisplay(plotter, sets, index=receiver, start=start,
                                 end=end, count=count, mode=mode, label=label)
    add_controls(plotter, display, sets, panel=panel, font=font,
                 room_bounds=room_bounds, save_dir=save_dir, fit=fit,
                 opacity_control=True)

    if off_screen:
        plotter.show(auto_close=False)
        plotter.screenshot(screenshot)
        print(f"[虚音源] 画像を保存しました: {screenshot}")
        vg.release_window(plotter, display)
        return display

    plotter.show()
    vg.release_window(plotter, display)
    return display


def add_controls(plotter, display, sets, panel=None, font=None,
                 room_bounds=None, save_dir=None, fit="all",
                 opacity_control=True, colour_key="m", help_window=True):
    """虚音源のスピンボックスとキーを組み立て、`(display, 操作の一覧)` を返す。

    ★**虚音源だけの画面（`view_images.view`）と、音線・音粒子と同居する画面
    （`view_rays.view`）の両方から呼ぶ**ので関数にしてある
    （2026-08-24 ユーザー指摘「虚音源を見るが最初の画面にあるのはなぜですか？
    音線確認画面と並列にある想定です」）。

    同居する画面では、すでにあるものと**取り合いにならないように**外から指定する。

    | 引数 | 同居する画面では | なぜ |
    |---|---|---|
    | `opacity_control` | False | 不透明度の欄は音線側が作っている |
    | `colour_key` | `"i"` | `m` は不透明度の表示 ON/OFF が使っている |
    | `help_window` | False | F1 は音線側が持っている（一覧に混ぜてもらう） |
    | `save_dir` | None | 画像保存の `g` も音線側が持っている |
    """
    if panel is not None:
        limit = display.order_limit
        if len(sets) > 1:
            panel.heading("受音点（w で切り替え）")
            panel.slider("受音点", [1.0, float(len(sets))], 1.0,
                         lambda v: display.set_receiver(int(round(v)) - 1),
                         fmt="%.0f")
            panel.text("濃い水色がいま選んでいる点です", size=8)
        else:
            panel.heading("受音点")
            panel.text(f"{sets[0].name} のみ", size=9)

        # ★開始と終了を別に持つ（「N 回目のみ」を出せるようにするため）
        # ★開始と終了を同じ値にすれば「N 回目のみ」になる（ユーザー要望）
        panel.heading("反射回数（開始＝終了 で N 回目のみ）")
        panel.slider("開始", [0.0, float(limit)], float(display.start),
                     display.set_start, fmt="%.0f")
        panel.slider("終了", [0.0, float(limit)], float(min(display.end, limit)),
                     display.set_end, fmt="%.0f")

        panel.text(f"上限は {limit} 回（計算時の最大反射回数）。"
                   f"つまみで足りないときは {vg.VALUE_INPUT_KEY} で数字を入れる",
                   size=8)

        panel.heading("表示する本数（強い順）")
        panel.slider("本数", [1.0, float(min(POOL_SIZE, max(len(s) for s in sets)))],
                     float(display.count), display.set_count, fmt="%.0f")

        if opacity_control:
            vg.add_opacity_control(plotter, font=font, panel=panel, target_key="o")

        panel.heading("虚音源の操作")
        panel.text(f"w 受音点を切り替え   {colour_key} 色を切り替え\n"
                   "l 線の表示 ON/OFF\n"
                   "h 室に寄る   r 全部が入る", size=8, color="#7f8794")
        keys = [
            "w  受音点を切り替え（濃い水色がいま選んでいる点）",
            f"{colour_key}  虚音源の色（反射回数 / エネルギー / 到来時刻）",
            "l  虚音源と受音点を結ぶ線の表示 ON/OFF",
            "h  室に寄る（虚音源は室の外へ遠く離れるので行き来する）",
            "r  いま描いているもの全部が入るように引く（VTK のキー）",
            f"{vg.VALUE_INPUT_KEY}  数値の欄をまとめて数字で入れる（欄ごとは … ボタン）",
        ]
        if opacity_control:
            keys.append("o  不透明度の対象を切り替え")
        if save_dir:
            keys.append("g  いまの画面を画像で保存")
        if help_window:
            keys += ["z / x / c / v  視点（上・正面・横・等角）",
                     "PageUp / PageDown  左の欄を送る", "q  閉じる"]
            panel.help_window(keys, title="虚音源 — 操作の一覧", note=IMAGE_NOTE)
        panel.enable_value_input()
        panel.enable_scroll()

    _register_keys(plotter, display, room_bounds=room_bounds, fit=fit,
                   colour_key=colour_key, save_dir=save_dir)
    return display, (keys if panel is not None else [])


# 虚音源の図を読むときの注意（操作の一覧の下に出す）
IMAGE_NOTE = ("虚音源と受音点を結ぶ直線は「折れた経路を伸ばしたもの」です。"
              "壁を貫いて見えるのが正しい姿で、"
              "線の長さが経路長（音速 × 到来時刻）になります。"
              "★描いているのは受音点に届いた経路の虚音源だけで、"
              "鏡像として作れる虚音源すべてではありません"
              "（届かないものは経路として成立していないので"
              "バックトレースが却下しています）。")


def _register_keys(plotter, display, room_bounds=None, fit="all",
                   colour_key="m", save_dir=None):
    """虚音源のキーを登録する。

    同居する画面（`view_rays`）と取り合いにならないよう、色のキーは外から渡す。
    """
    # ★目盛りは室の範囲で引き直す（虚音源に合わせると室の目盛りが読めない）。
    #   カメラは `fit="room"` のときだけ室に寄せる（既定は全部が入る形）
    if room_bounds is not None:
        _fit_to_room(plotter, room_bounds, camera=(fit == "room"))
        # ★★`h`（室に寄る）の登録が抜けていたことがある（2026-08-24）。
        #   操作の一覧には書いてあったのに**キーが効かなかった**。
        #   `r`（全部が入る）は VTK が元から持っているキーなのでそのまま使う
        plotter.add_key_event(
            "h", lambda: _fit_to_room(plotter, room_bounds, camera=True))
    plotter.add_key_event("w", display.next_receiver)
    plotter.add_key_event(colour_key, display.next_mode)
    plotter.add_key_event("l", display.toggle_lines)
    if save_dir:
        vg.add_screenshot_key(plotter, save_dir, "虚音源", key="g")
    return display


def _add_room_outline(plotter, model):
    """室の形（同一平面パッチの外周）を線で重ねる。

    ★壁を薄くしても形が分かるようにするためのもの（2026-08-24 ユーザー指摘
    「壁が見えなくなっている」）。三角形の辺を全部引くと網目になるので外周だけ。
    """
    try:
        segments = rd.patch_outline_segments(
            np.array([np.asarray(m.vertexes, dtype=float) for m in model.mesh]),
            np.array([np.asarray(m.normal, dtype=float) for m in model.mesh]))
    except Exception as error:
        print(f"[虚音源] 室の輪郭を描けませんでした: {type(error).__name__}: {error}")
        return None
    if not len(segments):
        return None
    points = segments.reshape(-1, 3)
    lines = np.column_stack([np.full(len(segments), 2),
                             np.arange(0, 2 * len(segments), 2),
                             np.arange(1, 2 * len(segments), 2)]).ravel()
    return plotter.add_mesh(pv.PolyData(points, lines=lines),
                            color=OUTLINE_COLOR, line_width=1.0, lighting=False,
                            opacity=0.75, reset_camera=False,
                            show_scalar_bar=False)


def _fit_to_room(plotter, bounds, camera=True):
    """カメラと目盛りを**室の範囲**に合わせ直す。

    虚音源は室の外へ遠く離れるので、VTK に任せると室が点になる（`view()` 参照）。
    目盛り（`show_bounds`）も描いてあるもの全部に合わせて伸びてしまうので、
    範囲を明示して引き直す。

    `camera=False` なら目盛りだけ引き直す。起動時は「全部が入る」ほうが
    虚音源と受音点を結ぶ線が見えるので、そちらを既定にしている。
    """
    if camera:
        try:
            plotter.reset_camera(bounds=bounds)
        except TypeError:           # 古い pyvista は bounds を取らない
            plotter.reset_camera()
        # ★**`show()` は最初の描画でカメラを引き直す**（描いてあるもの全部が
        #   入るように）。合わせたことを伝えておかないと、起動前にここで
        #   室に寄せても効かない（実際に `--fit room` が効かなかった）
        try:
            plotter.camera_set = True
        except Exception:
            pass
    try:
        plotter.show_bounds(bounds=bounds, grid="back", location="outer",
                            ticks="outside", font_size=9, color="#7f8794",
                            xtitle="X [m]", ytitle="Y [m]", ztitle="Z [m]")
    except Exception as error:      # 目盛りが出せなくても本体は見える
        print(f"[虚音源] 目盛りを引き直せませんでした: {type(error).__name__}: {error}")
    return bounds


def open_for_project(project, verbose=True, **kwargs):
    """プロジェクトを渡すだけで開く（`app.py` から呼ぶ入口）。"""
    sets = load_sets(project, verbose=verbose)
    return view(project.dxf_path, sets, absorption=project.absorption_path,
                unit=project.unit, band_number=project.band_number,
                orient_normals=project.orient_normals,
                source_points=None if project.source is None else [project.source],
                save_dir=project.screenshot_dir(), **kwargs)


def main():
    import argparse

    import project as pj

    p = argparse.ArgumentParser(
        description="虚音源を可視化する（パルス列から復元。計算はしない）")
    p.add_argument("folder", nargs="?", default=".",
                   help="プロジェクトフォルダ（project.json があるところ）")
    p.add_argument("--start", type=int, default=DEFAULT_ORDER_START,
                   help="描く反射回数の開始（既定 0 = 直接音から）")
    p.add_argument("--end", type=int, default=DEFAULT_ORDER_END,
                   help="描く反射回数の終了（開始と同じ値なら「N 回目のみ」）")
    p.add_argument("--count", type=int, default=DEFAULT_COUNT,
                   help="表示する本数（エネルギーの大きい順）")
    p.add_argument("--colour", default="reflection", choices=list(COLOUR_MODES),
                   help="色の付け方")
    p.add_argument("--opacity", type=float, default=DEFAULT_OPACITY,
                   help="壁の不透明度（既定 %(default)s。外周の線は常に出る）")
    p.add_argument("--receiver", type=int, default=1,
                   help="どの受音点を相手にするか（1 始まり。画面では w キー）")
    p.add_argument("--fit", default="all", choices=("all", "room"),
                   help="視野。all=虚音源まで全部入る / room=室に寄る（画面では h）")
    p.add_argument("--screenshot", default=None, help="画像に保存して閉じる")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf:
        raise SystemExit(f"{a.folder} に project.json がありません")
    open_for_project(project, start=a.start, end=a.end, count=a.count,
                     mode=a.colour, opacity=a.opacity, screenshot=a.screenshot,
                     receiver=max(0, a.receiver - 1), fit=a.fit)
    return 0


if __name__ == "__main__":
    try:
        code = main()
    finally:
        pass
    sys.exit(code)
