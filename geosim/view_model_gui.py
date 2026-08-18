"""読み込んだモデルを 3D 表示するビューア（PyVista / VTK のネイティブウィンドウ版）。

`view_model.py`（HTML + WebGL 版）と表示内容は同じで、実装だけが違う。

  view_model.py      … 依存ライブラリなし。HTML を書き出してブラウザで開く。
                       共有しやすい（相手に環境構築を求めない）。
  view_model_gui.py  … 本ファイル。PyVista で OS のウィンドウを開く。
                       Python 側からそのまま操作でき、音線・音粒子の重ね描き（G-1/G-2）や
                       将来の GUI 統合（G-7）に発展させやすい。

表示内容:
  ・三角形要素（辺を描くので分割が見える）
  ・法線ベクトル（矢印）
  ・**法線の裏側を赤で塗る** ← 向きの誤りが一目で分かる
  ・レイヤ別の色分け・チェックボックスで表示切り替え
  ・音源 / 受音点

操作:
  ドラッグ        回転
  ホイール        拡大縮小
  中ドラッグ      平行移動
  z / x / c / v   上 / 正面 / 横 / 等角 の視点
  n               法線矢印の表示切り替え
  g               いまの画面をそのまま画像で保存（`add_screenshot_key` を付けた画面のみ）
  t               スライダの値を数字で入力（`enable_value_input` を付けた画面のみ）
  w / s           ワイヤフレーム / 面（VTK の既定キー）
  r               視点リセット、q でウィンドウを閉じる

★**キーを割り当てる前に `VTK_RESERVED_KEYS` を見ること。**
  `e` と `q` は VTK の終了キー。`e` に数値入力を割り当てていたため、
  値を入れた瞬間にウィンドウが閉じていた（2026-08-17 に `t` へ変更）。

★**ウィンドウのタイトルに日本語を出すと化ける**（Windows）。
  `set_window_title()` が検証して、駄目なら英字の題に落とす。事情はその関数に書いてある。

使い方:
    cd geosim
    python view_model_gui.py ..\\test2.dxf
    python view_model_gui.py ..\\test.dxf --absorption ..\\absorption.csv
    python view_model_gui.py ..\\test.dxf --screenshot shot.png   # 画像だけ書き出す
"""

import argparse
import ctypes
import os
import sys

import numpy as np
import pyvista as pv

import read_dxffile as rd
from view_model import LAYER_PALETTE

# 法線の裏側の色。HTML 版と同じ赤にしてある
BACK_COLOR = "#C24540"

# 裏面の不透明度を表面の何倍にするか。
# **裏から見ている面はより透ける**ようにすると、室の外から中を覗いたときに
# 手前の壁（法線が内向きなので裏から見ることになる）が邪魔をしにくい。
# 0 にすると裏面が完全に消えて「向きの誤り」に気づけなくなるので、薄く残す。
BACKFACE_OPACITY_RATIO = 0.3
BG_BOTTOM = "#1c2027"
BG_TOP = "#2e3540"
TEXT_COLOR = "#d6dae2"

# VTK の既定フォントは日本語を持っていないので、レイヤ名（＝吸音材名）が
# 豆腐になる。Windows 標準の日本語フォントを順に探して使う。
_FONT_CANDIDATES = [
    r"C:\Windows\Fonts\meiryo.ttc",
    r"C:\Windows\Fonts\YuGothR.ttc",
    r"C:\Windows\Fonts\msgothic.ttc",
    r"C:\Windows\Fonts\NotoSansJP-Regular.ttf",
]


def japanese_font():
    """使える日本語フォントのパスを返す。見つからなければ None。"""
    for path in _FONT_CANDIDATES:
        if os.path.exists(path):
            return path
    return None


# 左の操作パネルが占める幅の割合。残りが 3D 表示になる
PANEL_RATIO = 0.25

# ★**VTK が既定で握っているキー。ここに機能を割り当ててはいけない。**
#
# `e` に数値入力を割り当てていたせいで、押すと入力ダイアログが出ると同時に
# **VTK の終了処理（ExitEvent）も走り**、値を入れて OK を押した瞬間に
# ウィンドウが閉じていた（ユーザー報告 2026-08-17）。`q` と同じ扱いだと気づいていなかった。
# 実測で確かめた結果が下記（`InvokeEvent('CharEvent')` で再現できる）。
VTK_RESERVED_KEYS = {
    "e": "終了（ExitEvent。q と同じ）",
    "q": "終了（ExitEvent）",
    "w": "ワイヤフレーム表示",
    "s": "面表示",
    "r": "視点リセット",
    "f": "注視点へ寄る",
    "p": "ピック（面の選択）",
    "u": "ユーザーイベント",
    "3": "ステレオ表示",
}

# 数値入力を開くキー。**予約キーを避けること**（上記）
VALUE_INPUT_KEY = "t"


class PanelItem:
    """パネルに積んだ要素 1 つ。**高さと、置き直し方**だけを持つ。

    ウィンドウの大きさが変わったら `place(y)` を呼び直すだけで並べ直せる。
    """

    __slots__ = ("height", "place")

    def __init__(self, height, place):
        self.height = height
        self.place = place


class ControlPanel:
    """ウィンドウ左端の操作パネル（2026-08-15）。

    以前はレイヤのチェックボックス・不透明度スライダ・操作説明を
    **3D 表示の上に直接重ねて**いたため、モデルが見えなくなるうえ、
    要素どうしが重なって読めなくなっていた。
    そこで**ウィンドウを左右に分け、左 1/4 を操作パネル専用**にした
    （`pv.Plotter(shape=(1,2), col_weights=[...])` の 2 レンダラ構成）。

    ★**大きさが変わったら並べ直す**（2026-08-15 の作り直し）。
    最初の版は「文字は絶対ピクセル・スライダは正規化座標」で置いていたので、
    **ウィンドウを最大化すると両者がずれて重なった**（ユーザー指摘）。
    いまは要素を `PanelItem`（高さ ＋ 置き直し方）として覚えておき、
    リサイズのたびに上から積み直す。位置を決めるのは `relayout()` 1 か所だけ。

    文字はレンダラのビューポート内の座標で置かれるので、パネル側のレンダラに
    追加すれば自然に左側へ収まる。ウィジェット（チェックボックス・スライダ）も
    **選ばれているレンダラのビューポート**が基準なので、置く前にパネル側を選ぶ。
    """

    CHECK = 18          # チェックボックスの一辺 [px]
    SLIDER = 34         # スライダのバーぶんの高さ [px]（見出しは別に 1 行取る）
    LABEL_FONT = 8      # レイヤ名など、横に長くなりがちな文字
    LINE = 18           # font_size 9 のときの 1 行 [px]（外から参照される既定値）

    def __init__(self, plotter, font=None, margin=14, width_ratio=PANEL_RATIO):
        self.plotter = plotter
        self.font = font
        self.margin = margin
        self.ratio = width_ratio
        self.items = []
        self._widgets = []                 # 参照を残さないと GC で消える
        self.controls = []                 # 数値入力ダイアログに出すスライダ
        self._measure()
        self._watch_resize()

    # ---- 数値を直接入力する ---------------------------------------------

    def enable_value_input(self, key=VALUE_INPUT_KEY):
        """`key` を押すと、スライダの値を**数字で打ち込める**ダイアログを出す。

        VTK のスライダはつまみを動かすことしかできず、
        「音線を 137 本にしたい」のような指定ができない（ユーザー指摘）。
        tkinter の小さな窓を出して、いまある全スライダをまとめて入力する。

        3D の描画とは別のイベントループになるが、**ダイアログを閉じるまで
        そちらに入りきり**なので取り合いにはならない。
        """
        if not self.controls:
            return
        if key in VTK_RESERVED_KEYS:
            # 黙って壊れるより、気づけるようにしておく
            print(f"[view] 警告: キー {key!r} は VTK が使っています"
                  f"（{VTK_RESERVED_KEYS[key]}）。別のキーにしてください")
        self.plotter.add_key_event(key, self.open_value_input)

    def open_value_input(self):
        import tkinter as tk
        from tkinter import ttk

        root = tk.Tk()
        root.title("値を入力")
        root.attributes("-topmost", True)
        frame = ttk.Frame(root, padding=12)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="数値を入れて OK を押してください",
                  foreground="#666").grid(row=0, column=0, columnspan=3,
                                          sticky="w", pady=(0, 8))

        entries = []
        for row, control in enumerate(self.controls, start=1):
            low, high = control["range"]
            ttk.Label(frame, text=control["label"]).grid(row=row, column=0,
                                                        sticky="w", pady=3)
            var = tk.StringVar(value=control["format"] % control["value"])
            ttk.Entry(frame, textvariable=var, width=12).grid(row=row, column=1,
                                                              padx=8)
            ttk.Label(frame, text=f"（{low:g} 〜 {high:g}）",
                      foreground="#666").grid(row=row, column=2, sticky="w")
            entries.append((control, var))

        def apply():
            for control, var in entries:
                text = var.get().strip()
                if not text:
                    continue
                try:
                    value = float(text)
                except ValueError:
                    continue
                low, high = control["range"]
                control["set"](float(min(max(value, low), high)))
            root.destroy()

        buttons = ttk.Frame(frame)
        buttons.grid(row=len(entries) + 1, column=0, columnspan=3,
                     sticky="e", pady=(12, 0))
        ttk.Button(buttons, text="キャンセル", command=root.destroy).pack(side="right",
                                                                        padx=4)
        ttk.Button(buttons, text="OK", command=apply).pack(side="right", padx=4)
        root.bind("<Return>", lambda _e: apply())
        root.mainloop()

    def _measure(self):
        width, height = self.plotter.window_size
        self.width = max(120, int(width * self.ratio))
        self.height = max(120, int(height))

    # ---- レンダラの選択 ------------------------------------------------

    def _panel(self):
        """パネル側のレンダラを選ぶ（文字やウィジェットを置く前に呼ぶ）。"""
        self.plotter.subplot(0, 0)

    def _model(self):
        """3D 側のレンダラに戻す。"""
        self.plotter.subplot(0, 1)

    # ---- 並べ直し ------------------------------------------------------

    def _add(self, height, place):
        item = PanelItem(height, place)
        self.items.append(item)
        return item

    def relayout(self, render=False):
        """いまのウィンドウの大きさに合わせて、上から順に積み直す。"""
        self._measure()
        y = self.height - self.margin
        for item in self.items:
            y -= item.height
            item.place(y)
        if render:
            self.plotter.render()

    def _watch_resize(self):
        """ウィンドウの大きさが変わったら並べ直す。

        VTK は `ConfigureEvent` でリサイズを知らせる。
        off_screen では interactor が無いことがあるので、失敗しても黙って諦める
        （その場合は作ったときの大きさのまま。画像書き出しなので問題ない）。
        """
        try:
            self.plotter.iren.add_observer("ConfigureEvent",
                                           lambda *a: self.relayout())
        except Exception:
            pass

    # ---- 幅に収める ----------------------------------------------------

    @staticmethod
    def line_height(size):
        """フォントサイズから 1 行の高さ [px] を見積もる。

        VTK は行間を含めておよそ文字サイズの 2 倍を使う。
        小さく見積もると**文字どうしが重なる**（実際に重なった）。
        """
        return int(round(size * 2.0))

    def fit(self, text, available, size=None):
        """パネルの幅に収まるように文字を切り詰める（末尾を … にする）。

        レイヤ名は日本語で長くなりがちで、そのままだとパネルからはみ出して
        **途中で切れて読めなくなる**。全角は文字サイズの約 1.9 倍、
        半角は約 1 倍として見積もる（正確な幅は VTK 側でしか測れないので概算）。
        """
        size = size or self.LABEL_FONT
        width = 0.0
        for i, ch in enumerate(text):
            width += size * (1.9 if ord(ch) > 0x2000 else 1.0)
            if width > available:
                return text[:max(1, i - 1)] + "…"
        return text

    # ---- 要素を積む ----------------------------------------------------

    def gap(self, pixels=10):
        self._add(pixels, lambda y: None)
        return self

    def text(self, message, size=9, color="#9aa2b1", indent=0):
        """複数行の文字を置く。

        ★VTK は**渡した位置がブロックの下端**で、そこから上へ書いていく。
        だから高さを引いた位置に置くのが正しい。
        """
        height = self.line_height(size)
        lines = str(message).count("\n") + 1
        self._panel()
        actor = self.plotter.add_text(" ", position=(self.margin + indent, 0),
                                      font_size=size, color=color,
                                      font_file=self.font)
        self._model()

        def place(y, actor=actor, message=message, size=size, indent=indent):
            # 幅は並べ直しのたびに測り直す（ウィンドウを広げたら切り詰めが緩む）
            available = self.width - 2 * self.margin - indent
            actor.SetInput("\n".join(self.fit(line, available, size)
                                     for line in str(message).split("\n")))
            actor.SetDisplayPosition(self.margin + indent, int(y))

        self._add(height * lines, place)
        return self

    def reserve_text(self, lines, size=9, color=TEXT_COLOR):
        """あとから書き換える文字のために場所を空け、その actor を返す。

        行数が変わりうる表示（判定の内訳など）に使う。
        **最大行数ぶんを確保**しておけば、行が減っても上に隙間ができるだけで重ならない。
        """
        self._panel()
        actor = self.plotter.add_text(" ", position=(self.margin, 0),
                                      font_size=size, color=color,
                                      font_file=self.font)
        self._model()
        self._add(self.line_height(size) * lines,
                  lambda y, a=actor: a.SetDisplayPosition(self.margin, int(y)))
        return actor

    def heading(self, message):
        self.gap(8)
        return self.text(message, size=10, color=TEXT_COLOR)

    def checkbox(self, label, value, callback, colour="#4cc9f0"):
        """チェックボックス 1 つと、その右のラベル。"""
        self._panel()
        widget = self.plotter.add_checkbox_button_widget(
            callback, value=value, position=(self.margin, 0),
            size=self.CHECK, border_size=2, color_on=colour,
            color_off="#454c58", background_color="#2b303a")
        text = self.plotter.add_text(" ", position=(0, 0),
                                     font_size=self.LABEL_FONT, color=TEXT_COLOR,
                                     font_file=self.font)
        self._model()
        self._widgets.append(widget)

        def place(y, widget=widget, text=text, label=label):
            left = self.margin + self.CHECK + 8
            representation = widget.GetRepresentation()
            representation.SetPlaceFactor(1.0)
            representation.PlaceWidget([self.margin, self.margin + self.CHECK,
                                        y, y + self.CHECK, 0.0, 0.0])
            text.SetInput(self.fit(label, self.width - left - self.margin))
            text.SetDisplayPosition(left, int(y) + 3)

        self._add(self.CHECK + 4, place)
        return widget

    def slider(self, title, value_range, value, callback, fmt="%.2f"):
        """横向きのスライダ。見出しと値は**自前の文字**で出す。

        VTK のスライダにも見出し・値の表示はあるが、こちらで位置を動かしたときに
        文字が追随せず**消えたように見えた**。自分で描けば位置も内容も確実で、
        ついでに**日本語の見出しが使える**（VTK の既定フォントは日本語を持たない）。
        """
        label = self.reserve_text(1, size=9)

        def show_value(v):
            set_actor_text(label, f"{title}  {fmt % v}")

        show_value(value)
        self._panel()

        # ★スライダは**生成時にコールバックを 1 回呼ぶ**。そのときはパネル側の
        #   レンダラが選ばれているので、コールバックの中で 3D の actor を作り直すと
        #   **パネルの中に 3D が描かれてしまう**（実際に音線がパネルに出た）。
        #   生成時の 1 回は無視し、以降は必ず 3D 側を選んでから呼ぶ。
        state = {"ready": False}

        def guarded(v):
            if not state["ready"]:
                return
            control["value"] = float(v)
            show_value(v)
            self._model()
            callback(v)

        widget = self.plotter.add_slider_widget(
            guarded, value_range, value=value, title=None,
            pointa=(0.0, 0.0), pointb=(1.0, 0.0), style="modern",
            color=TEXT_COLOR, tube_width=0.004, slider_width=0.018)
        # VTK 側の値表示は消す（自前のラベルと二重になるため）
        widget.GetRepresentation().ShowSliderLabelOff()
        state["ready"] = True
        self._model()
        self._widgets.append(widget)

        def set_value(v, widget=widget):
            """つまみを動かさずに値を入れる（数値入力ダイアログから使う）。"""
            widget.GetRepresentation().SetValue(float(v))
            control["value"] = float(v)
            show_value(v)
            self._model()
            callback(v)
            self.plotter.render()

        control = {"label": title, "range": tuple(value_range), "value": float(value),
                   "format": fmt, "set": set_value, "widget": widget}
        self.controls.append(control)

        def place(y, widget=widget):
            # ★位置は**ウィンドウ全体**の正規化座標（`Normalized Display`）。
            #   パネル側のレンダラに置いても、ビューポート基準にはならない。
            #   パネルの幅で割ると右端がウィンドウの外まで伸び、
            #   見えている部分だけが切り取られて**つまみの効く範囲がずれる**
            window_width, window_height = self.plotter.window_size
            representation = widget.GetRepresentation()
            representation.GetPoint1Coordinate().SetValue(
                self.margin / window_width, (y + self.SLIDER * 0.4) / window_height)
            representation.GetPoint2Coordinate().SetValue(
                (self.width - self.margin) / window_width,
                (y + self.SLIDER * 0.4) / window_height)
            # 座標を変えただけでは描き直されない（内部の計算結果を持っているため）
            representation.Modified()
            representation.BuildRepresentation()

        self._add(self.SLIDER, place)
        return widget


# ------------------------------------------------------------------------------
# ウィンドウのタイトル（日本語が化ける件の対処）
# ------------------------------------------------------------------------------

# 画面ごとのウィンドウタイトル。(日本語, 化けたときの英字) の組。
#
# **物件名は入れない**（ユーザー判断 2026-08-17）。タイトルバーに欲しいのは
# 「法線の確認なのか、音線の可視化なのか」だけ。物件名は画面の中の見出しに出ている。
# ここに集めてあるので、画面が増えたら 1 行足すだけで済む。
WINDOW_TITLES = {
    "normals":    ("法線の確認",     "geosim - normals"),
    "rays":       ("音線・音粒子",   "geosim - rays and particles"),
    "directions": ("音線の飛び方",   "geosim - ray directions"),
    "model":      ("モデルビューア", "geosim - model viewer"),
}


def ascii_tag(text):
    """ASCII だけで書かれていればそのまま返す。そうでなければ空。

    ウィンドウのタイトルに使えるかどうかの判定に使う（下の事情による）。
    """
    text = (text or "").strip()
    return text if text and all(ord(c) < 128 for c in text) else ""


def window_titles(screen, default=""):
    """画面の種類から (タイトル, 化けたときの英字) を引く。

    知らない種類なら `default` をそのまま使い、英字の予備は
    そこから ASCII だけ抜いて作る（無ければ `geosim`）。
    """
    if screen in WINDOW_TITLES:
        return WINDOW_TITLES[screen]
    return default, (ascii_tag(default) or "geosim")


def _window_handle(plotter):
    """Plotter から Windows のウィンドウハンドルを取り出す（無ければ None）。

    VTK は `_00000000001f026e_p_void` のような SWIG 形式の文字列で返してくる。
    """
    try:
        raw = str(plotter.ren_win.GetGenericWindowId())
    except Exception:
        return None
    try:
        return int(raw.split("_")[1], 16) if raw.startswith("_") else int(raw, 0)
    except (ValueError, IndexError):
        return None


def _shown_window_title(plotter):
    """いま実際にタイトルバーに出ている文字列を読む（Windows のみ）。"""
    handle = _window_handle(plotter)
    if handle is None:
        return None
    buffer = ctypes.create_unicode_buffer(512)
    ctypes.windll.user32.GetWindowTextW(ctypes.c_void_p(handle), buffer, 512)
    return buffer.value


def set_window_title(plotter, title, fallback=None):
    """ウィンドウのタイトルを設定する。**化けたら ASCII の題に落とす**。

    【なぜこんなことをするか】2026-08-16
    VTK が Windows に作るウィンドウは **ANSI ウィンドウ**（`IsWindowUnicode` が偽）で、
    日本語のタイトルが**文字によって**化ける。`研修室 — 法線の確認` が
    `研修室 ?E法線?E確?E` になった（研・修・室・法・線・確は通るのに
    —・の・認 が化ける）。ユーザーからの指摘で発覚。

    手は一通り試して、どれも駄目だった:

      ・`SetWindowTextW`（ワイド文字版）→ **かえって全部 `?` になる**
        （ANSI ウィンドウなので書き込み時に ANSI へ変換されてしまう）
      ・`SendMessageW` で `WM_SETTEXT` → 同上
      ・`SetWindowTextA` に CP932 のバイト列 → やはり全部 `?`
        （このウィンドウの ANSI↔Unicode 変換は CP932 ではない）

    そこで**設定したあと読み戻して検証**し、一致しなければ ASCII の題に差し替える。
    化けた文字列を出したままにするよりは、読める英字のほうがましという判断。
    **画面の中の見出しは日本語のまま**（フォントを指定して描いているので化けない）。

    Windows 以外では検証を飛ばす（化ける現象自体が Windows のもの）。

    ★**決まった題は `plotter.title` にも書き戻す。** pyvista の `show()` は
      そのつど `self.title` をウィンドウ名に**貼り直す**ので、書き戻しておかないと
      せっかく差し替えた題が日本語（＝化ける方）に戻ってしまう。実際に戻った。
    """
    def remember(text):
        try:                       # 版によっては書けないので黙って諦める
            plotter.title = text
        except Exception:
            pass
        return text

    try:
        plotter.ren_win.SetWindowName(title)
    except Exception:
        return title
    if sys.platform != "win32":
        return remember(title)

    shown = _shown_window_title(plotter)
    if shown is None or shown == title:
        return remember(title)

    safe = fallback or "geosim"
    try:
        plotter.ren_win.SetWindowName(safe)
    except Exception:
        pass
    return remember(safe)


def prepare_window(plotter, title, fallback=None):
    """**`show()` の直前に呼ぶ。** ウィンドウを先に作らせてタイトルを確定させる。

    OS のウィンドウは最初の描画で作られるので、それより前にタイトルを
    検証することはできない。かといって `show()` はブロックするので後からも呼べない。
    そこで**ここで 1 回だけ描いてウィンドウを作らせ**、そのうえで設定する。

    ★**描画イベントの中で設定してはいけない。** 最初はそうしていたが、
      `SetWindowName` は `WM_SETTEXT` を同期で送るため、描画の途中で
      ウィンドウプロシージャが再入し、**OpenGL のコンテキストごと落ちた**
      （segfault と、以後のプロセスで `wglChoosePixelFormatARB` に失敗する状態）。
      描画の外から呼べばこの問題は起きない。
    """
    try:
        plotter.ren_win.Render()
    except Exception:
        return title
    return set_window_title(plotter, title, fallback)


# ------------------------------------------------------------------------------
# 画面の保存
# ------------------------------------------------------------------------------

def next_free_path(folder, stem, suffix=".png", digits=2):
    """`stem_01.png` … まだ空いている連番のパスを返す。

    **撮るたびに増える**（上書きしない）。角度を変えて何枚も撮るための配慮。
    """
    os.makedirs(folder, exist_ok=True)
    number = 1
    while True:
        path = os.path.join(folder, f"{stem}_{number:0{digits}d}{suffix}")
        if not os.path.exists(path):
            return path
        number += 1


def add_screenshot_key(plotter, folder, stem, key="p"):
    """`p` で「いま画面に出ているとおり」を PNG に保存する。

    **左のパネルごと写す。** どのレイヤを出していたか・不透明度がいくつだったかも
    一緒に残るほうが、あとで見返したときに条件が分かるため。

    `stem` は文字列でも、呼ぶたびにファイル名の芯を返す関数でもよい
    （音粒子は時刻をファイル名に入れたいので関数を渡す）。
    """
    def save():
        name = stem() if callable(stem) else stem
        path = next_free_path(folder, name)
        try:
            plotter.screenshot(path)
        except Exception as e:
            print(f"[view] 画像を保存できませんでした: {type(e).__name__}: {e}")
            return None
        print(f"[view] 画像を保存しました: {path}")
        return path

    plotter.add_key_event(key, save)
    return save


def make_plotter(title, window_size, off_screen, panel=True, screen=None):
    """左に操作パネル、右に 3D を持つ Plotter を作る。

    `panel=False` なら従来どおり 1 画面（画像の書き出しなど、操作しないとき用）。

    `title` は**画面の中の見出し**（物件名を含んでよい。日本語も化けない）。
    `screen` は画面の種類（`'normals'` など）で、**ウィンドウのタイトル**を決める。
    分けてあるのは、タイトルバーには物件名が要らないという判断のため
    （`WINDOW_TITLES` 参照）。
    """
    bar, fallback = window_titles(screen, title)
    if not panel:
        plotter = pv.Plotter(window_size=window_size, title=bar,
                             off_screen=off_screen)
        plotter.set_background(BG_BOTTOM, top=BG_TOP)
    else:
        plotter = pv.Plotter(shape=(1, 2),
                             col_weights=[PANEL_RATIO, 1.0 - PANEL_RATIO],
                             window_size=window_size, title=bar,
                             off_screen=off_screen, border=False)
        plotter.subplot(0, 0)
        plotter.set_background(BG_BOTTOM)
        plotter.subplot(0, 1)
        plotter.set_background(BG_BOTTOM, top=BG_TOP)

    # タイトルの確定は `show()` の直前（`prepare_window`）。ここではまだ
    # OS のウィンドウが無いので検証できない
    _attach(plotter, "geosim_title", None if off_screen else (bar, fallback))
    if not panel:
        return plotter, None
    return plotter, ControlPanel(plotter, font=japanese_font(),
                                 width_ratio=PANEL_RATIO)


def finish_window(plotter):
    """`make_plotter` に渡した題で、`show()` の直前にタイトルを確定させる。

    各ビューアが `plotter.show()` を呼ぶ直前にこれを挟むだけで済むよう、
    題は `make_plotter` の時点で Plotter に持たせてある。
    """
    stored = (getattr(plotter, "geosim_title", None)
              or getattr(plotter, "_geosim_title", None))
    if not stored:
        return None
    return prepare_window(plotter, stored[0], stored[1])


def triangles_to_polydata(triangles):
    """Mesh のリストを pv.PolyData にする。

    巻き順を `t.normal` に合わせて並べ替えているのが要点。
    read_dxffile の `orient_normals='flip'` / `'shells'` は**法線だけを反転**して
    頂点の順序は触らないため、そのまま渡すと VTK 側の表裏（＝backface の判定）が
    モデルの持つ法線と食い違ってしまう。ここで揃えておけば、VTK に裏面を
    赤で塗らせるだけで「法線の向きの確認」がそのまま成立する。
    """
    n_tri = len(triangles)
    points = np.empty((n_tri * 3, 3), dtype=float)
    faces = np.empty((n_tri, 4), dtype=np.int64)
    normals = np.empty((n_tri, 3), dtype=float)

    for i, t in enumerate(triangles):
        v = np.asarray(t.vertexes, dtype=float)
        n = np.asarray(t.normal, dtype=float)
        points[3 * i:3 * i + 3] = v
        normals[i] = n

        geometric = np.cross(v[1] - v[0], v[2] - v[0])
        if float(np.dot(geometric, n)) >= 0.0:
            order = (0, 1, 2)
        else:
            order = (0, 2, 1)
        faces[i] = (3, 3 * i + order[0], 3 * i + order[1], 3 * i + order[2])

    poly = pv.PolyData(points, faces.ravel())
    poly.cell_data["normal"] = normals
    return poly


def normal_arrows(poly, length):
    """面の重心から法線方向に伸びる矢印を作る。"""
    centres = poly.cell_centers()
    centres.point_data["normal"] = poly.cell_data["normal"]
    return centres.glyph(orient="normal", scale=False, factor=length,
                         geom=pv.Arrow(tip_length=0.3, tip_radius=0.09,
                                       shaft_radius=0.03))


def build_plotter(model, title="モデルビューア", off_screen=False,
                  show_normals=True, normal_ratio=0.06, window_size=(1280, 860),
                  opacity=1.0, show_bounds=True, show_summary=True,
                  layer_opacity=None, panel=True, screen=None):
    """DxfModel から Plotter を組み立てて返す（show() はしない）。

    opacity … 面の不透明度。音線を重ねるときは 0.15 くらいにすると中が見える
    layer_opacity … {レイヤ名: 不透明度} でレイヤごとに指定する。
        指定の無いレイヤは opacity を使う。
        「床だけ残して壁を消す」といった見方ができる
    show_bounds … 目盛り付きの箱を描くか
    show_summary … 読み込み結果のサマリを左下に出すか（音線を重ねるときは邪魔）

    戻り値の Plotter には `geosim_layers` を付けてある。
    {レイヤ名: {'face', 'arrow', 'colour', 'opacity'}} で、
    あとから不透明度や表示を変えるのに使う（`add_opacity_control` が利用する）。
    """
    layer_opacity = dict(layer_opacity or {})
    mesh = model.mesh
    if not mesh:
        raise ValueError("表示できる三角形がありません")

    layers = sorted({t.material for t in mesh})
    font = japanese_font()

    lo, hi = model.extents
    diag = float(np.linalg.norm(np.asarray(hi) - np.asarray(lo))) or 1.0
    arrow_len = diag * normal_ratio

    plotter, panel = make_plotter(title, window_size, off_screen, panel=panel,
                                  screen=screen)

    face_actors = {}
    arrow_actors = {}
    for i, name in enumerate(layers):
        colour = LAYER_PALETTE[i % len(LAYER_PALETTE)]
        poly = triangles_to_polydata([t for t in mesh if t.material == name])
        alpha = float(layer_opacity.get(name, opacity))

        face_actors[name] = plotter.add_mesh(
            poly, color=colour, show_edges=True, edge_color=BG_BOTTOM,
            line_width=1, lighting=True, ambient=0.32, diffuse=0.70,
            specular=0.06, smooth_shading=False, opacity=alpha,
            backface_params={"color": BACK_COLOR, "ambient": 0.32,
                             "diffuse": 0.70,
                             "opacity": alpha * BACKFACE_OPACITY_RATIO},
        )
        arrows = normal_arrows(poly, arrow_len)
        arrow_actors[name] = plotter.add_mesh(arrows, color="#f2f4f8",
                                              lighting=False)
        arrow_actors[name].SetVisibility(show_normals)
        layer_opacity[name] = alpha

    marker_radius = diag * 0.012
    for point in model.source_points:
        plotter.add_mesh(pv.Sphere(radius=marker_radius, center=np.asarray(point)),
                         color="#ff5f5f", lighting=False)
    for point in model.receiver_points:
        plotter.add_mesh(pv.Sphere(radius=marker_radius, center=np.asarray(point)),
                         color="#4dd0a0", lighting=False)

    plotter.add_axes(color=TEXT_COLOR)
    if show_bounds:
        plotter.show_bounds(grid="back", location="outer", ticks="outside",
                            font_size=9, color="#7f8794", xtitle="X [m]",
                            ytitle="Y [m]", ztitle="Z [m]")

    if panel is None:
        # パネルが無いときだけ 3D の上に文字を重ねる（画像書き出しなど）
        plotter.add_text(f"{title}\n三角形 {len(mesh)} 枚 / レイヤ {len(layers)}",
                         position="upper_left", font_size=11,
                         color=TEXT_COLOR, font_file=font)
        if show_summary:
            plotter.add_text(model.summary(), position=(12, 12), font_size=8,
                             color="#9aa2b1", font_file=font)

    # ---- 視点プリセット（VTK 既定の w/s/r/q とぶつからないキーを選ぶ） ----
    plotter.add_key_event("z", plotter.view_xy)
    plotter.add_key_event("x", plotter.view_xz)
    plotter.add_key_event("c", plotter.view_yz)
    plotter.add_key_event("v", plotter.view_isometric)

    state = {"normals": show_normals}

    def toggle_normals():
        state["normals"] = not state["normals"]
        for actor in arrow_actors.values():
            actor.SetVisibility(state["normals"])
        plotter.render()

    plotter.add_key_event("n", toggle_normals)

    # ---- 左パネルにレイヤの表示切り替えを並べる ----
    if panel is not None:
        panel.text(f"{title}", size=11, color=TEXT_COLOR)
        panel.text(f"三角形 {len(mesh)} 枚 / レイヤ {len(layers)}", size=9)
        panel.heading("レイヤ表示")
        for i, name in enumerate(layers):
            count = model.layer_counts.get(name, 0)
            panel.checkbox(f"{name} ({count})", True,
                           _visibility_callback(plotter, face_actors[name],
                                                arrow_actors[name], state),
                           colour=LAYER_PALETTE[i % len(LAYER_PALETTE)])

    plotter.view_isometric()
    # あとから不透明度や表示を変えられるよう、レイヤごとの actor を Plotter に持たせる。
    # pyvista は新しい公開属性の追加を禁じているので、専用の API を使う
    # （無い版のために private 名へのフォールバックも用意しておく）
    registry = {name: {"face": face_actors[name], "arrow": arrow_actors[name],
                       "colour": LAYER_PALETTE[i % len(LAYER_PALETTE)],
                       "opacity": layer_opacity[name]}
                for i, name in enumerate(layers)}
    _attach(plotter, "geosim_layers", registry)
    _attach(plotter, "geosim_panel", panel)
    return plotter


def _attach(plotter, name, value):
    """Plotter に独自の情報を持たせる。

    pyvista は新しい公開属性の追加を禁じているので専用 API を使う
    （無い版のために private 名へのフォールバックも用意しておく）。
    """
    try:
        pv.set_new_attribute(plotter, name, value)
    except AttributeError:
        setattr(plotter, "_" + name, value)


def layer_actors(plotter):
    """`build_plotter` が登録したレイヤ情報を取り出す。"""
    return (getattr(plotter, "geosim_layers", None)
            or getattr(plotter, "_geosim_layers", None))


def control_panel(plotter):
    """`build_plotter` が作った左パネルを取り出す（無ければ None）。"""
    return (getattr(plotter, "geosim_panel", None)
            or getattr(plotter, "_geosim_panel", None))


def set_face_opacity(actor, value, ratio=BACKFACE_OPACITY_RATIO):
    """面の不透明度を変える。**裏面は表面より薄く**する。

    裏面は別の vtkProperty なので、表だけ変えると比率が崩れて
    「裏から見るほうが濃い」という妙な見え方になる。
    """
    actor.GetProperty().SetOpacity(value)
    back = actor.GetBackfaceProperty()
    if back is not None:
        back.SetOpacity(value * ratio)
    return actor


def set_actor_text(actor, text, corner=0):
    """`add_text()` が返す actor の文字を書き換える。

    位置を座標で渡すと `vtkTextActor`、`position='upper_right'` のような文字列だと
    `vtkCornerAnnotation` が返り、**差し替えの API が違う**。
    どちらが来ても書き換えられるようにここでまとめている。
    """
    if hasattr(actor, "SetInput"):
        actor.SetInput(text)
    else:
        actor.SetText(corner, text)
    return actor


def add_opacity_control(plotter, font=None, panel=None, target_key="o"):
    """モデルの不透明度を変えるスライダとキー操作を、左パネルに足す。

    - **スライダ** … 対象の不透明度を 0〜1 で設定する
    - **`o`**（`target_key`） … 対象を切り替える（すべて → 各レイヤ → すべて …）
    - **m** … モデル全体の表示 ON / OFF

    レイヤごとに変えられるようにしてあるのは、
    「壁だけ薄くして中の様子を見る」「床は残す」といった使い方をするため。

    ※ パネルのチェックボックスは**表示のオンオフ**（不透明度とは別）。
    """
    layers = layer_actors(plotter)
    if not layers:
        raise ValueError("build_plotter が作った Plotter を渡してください")
    if panel is None:
        panel = control_panel(plotter)
    if panel is None:
        raise ValueError("操作パネルがありません（build_plotter(panel=True) で作ること）")

    names = list(layers)
    # ready … スライダを作った直後と、対象切替で値を書き換えるときにコールバックが
    #         走ってしまうので、そのぶんを無視するための旗。
    #         これが無いと、生成時に「全レイヤの平均値」が全レイヤへ適用されてしまい、
    #         layer_opacity で個別に指定した値が消える
    state = {"target": 0, "visible": True, "ready": False}      # target 0 = すべて

    def target_name():
        return "すべて" if state["target"] == 0 else names[state["target"] - 1]

    def current_opacity():
        if state["target"] == 0:
            return float(np.mean([layers[n]["opacity"] for n in names]))
        return layers[target_name()]["opacity"]

    def refresh_label():
        set_actor_text(label, f"対象: {target_name()}\n"
                              f"{target_key} 対象切替  m 表示ON/OFF")

    def apply(value):
        if not state["ready"]:
            return
        value = float(np.clip(value, 0.0, 1.0))
        targets = names if state["target"] == 0 else [target_name()]
        for name in targets:
            layers[name]["opacity"] = value
            set_face_opacity(layers[name]["face"], value)
        refresh_label()
        plotter.render()

    panel.heading("透過設定")
    slider = panel.slider("不透明度", [0.0, 1.0], current_opacity(), apply)
    label = panel.reserve_text(2)
    state["ready"] = True

    def set_slider(value):
        state["ready"] = False          # 表示だけ更新し、適用はしない
        slider.GetRepresentation().SetValue(value)
        state["ready"] = True

    def next_target():
        state["target"] = (state["target"] + 1) % (len(names) + 1)
        set_slider(current_opacity())
        refresh_label()
        plotter.render()

    def toggle_model():
        state["visible"] = not state["visible"]
        for name in names:
            layers[name]["face"].SetVisibility(state["visible"])
        plotter.render()

    plotter.add_key_event(target_key, next_target)
    plotter.add_key_event("m", toggle_model)
    refresh_label()
    return slider


def _visibility_callback(plotter, face_actor, arrow_actor, state):
    """チェックボックス用のコールバックを作る（クロージャの取り違え防止）。"""
    def callback(flag):
        face_actor.SetVisibility(flag)
        arrow_actor.SetVisibility(flag and state["normals"])
        plotter.render()
    return callback


def view(dxf_path, absorption=None, unit=None, orient_normals="cad",
         screenshot=None, show_normals=True, opacity=1.0, layer_opacity=None):
    """DXF を読み込んで 3D ビューアのウィンドウを開く。"""
    model = rd.read_model(dxf_path, unit=unit, absorption_table=absorption,
                          orient_normals=orient_normals)
    base = os.path.splitext(os.path.basename(dxf_path))[0]
    plotter = build_plotter(model, title=f"{base} モデルビューア", screen="model",
                            off_screen=screenshot is not None,
                            show_normals=show_normals, opacity=opacity,
                            layer_opacity=layer_opacity)
    if screenshot is None:
        add_opacity_control(plotter, font=japanese_font())
        panel = control_panel(plotter)
        panel.heading("操作")
        panel.text("z/x/c/v 視点\nn 法線矢印\nw/s ワイヤ/面\nr リセット   q 終了",
                   color="#7f8794")
        panel.relayout()
    if screenshot is not None:
        plotter.screenshot(screenshot)
        plotter.close()
        print(f"\n[view_model_gui] 画像を書き出しました: {screenshot}")
    else:
        finish_window(plotter)
        plotter.show()
    return model


def main():
    p = argparse.ArgumentParser(
        description="DXF モデルの 3D ビューア（PyVista のネイティブウィンドウ）")
    p.add_argument("dxf", help="DXF ファイル")
    p.add_argument("--absorption", help="吸音率 CSV")
    p.add_argument("--unit", help="'mm' / 'm' など。省略すると $INSUNITS から自動判定")
    p.add_argument("--orient-normals", default="cad",
                   choices=["cad", "flip", "shells"])
    p.add_argument("--no-normals", action="store_true", help="法線矢印を最初は隠す")
    p.add_argument("--opacity", type=float, default=1.0,
                   help="面の不透明度（0=透明, 1=不透明）")
    p.add_argument("--layer-opacity",
                   help="レイヤごとの不透明度。例 \"1=0.6,2=0.05\"")
    p.add_argument("--screenshot", help="ウィンドウを開かず画像を書き出す（PNG）")
    a = p.parse_args()

    layer_opacity = None
    if a.layer_opacity:
        layer_opacity = {}
        for item in a.layer_opacity.split(","):
            name, value = item.rsplit("=", 1)
            layer_opacity[name.strip()] = float(value)

    view(a.dxf, absorption=a.absorption, unit=a.unit,
         orient_normals=a.orient_normals, screenshot=a.screenshot,
         show_normals=not a.no_normals, opacity=a.opacity,
         layer_opacity=layer_opacity)


if __name__ == "__main__":
    main()
