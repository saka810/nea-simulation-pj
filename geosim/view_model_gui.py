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
  数値の枠を押す  ★**押した桁から打ち込める**（Enter 確定 / Esc 取消）。別窓は開かない
                  ← → で桁を移動、BackSpace で前の 1 文字、Delete でその文字
  ▲▼             1 段ずつ増減（枠の右にくっついている。Excel と同じ形）
  t               数値の欄をまとめて入力（`enable_value_input` を付けた画面のみ）
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

# ★左の操作パネルは**別画面として扱う**ので、3D 側と背景を変える
# （2026-08-21 ユーザー要望）。3D 側は上が明るいグラデーション、
# パネル側は**平らで暗い**色にして「ここは設定の面」と分かるようにした。
# 明度差を付けたのは、境目がどの視点でもはっきり見えるようにするため
# （同系色の濃淡だと、グラデーションの暗い側と見分けが付かなくなる）。
PANEL_BG = "#0d1220"
# パネルと 3D の境目の線。細く明るい線を 1 本入れると「面が違う」ことが伝わる
PANEL_EDGE = "#44506a"
# パネルの見出し（画面の種類）の色
PANEL_TITLE_COLOR = "#8fb7ff"

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


def _caret_from_x(actor, text, left, x, renderer):
    """押した x が `text` の何文字目かを返す（文字の幅を測って決める）。

    ★数値の枠のどこを押したかで入力位置を決めるために使う（2026-08-24）。
    測れないときは末尾（`len(text)`）を返す——末尾なら必ず安全に編集できる。
    """
    try:
        keep = actor.GetInput()
        size = [0, 0]
        # ★**いちばん近い「文字の境目」**に付ける（「その文字より前／後」で
        #   決めると 1 桁ずれる。実際にずれた）
        best, gap = len(text), None
        for index in range(len(text) + 1):
            if index:
                actor.SetInput(text[:index])
                actor.GetSize(renderer, size)
                width = float(size[0])
            else:
                width = 0.0
            distance = abs(left + width - x)
            if gap is None or distance < gap:
                best, gap = index, distance
        actor.SetInput(keep)
        return best
    except Exception:
        return len(text)


def _rgb(colour):
    """`"#3a4150"` のような色を (r, g, b)（0〜1）にする。

    文字の背景・枠の色に使う（VTK は 0〜1 の実数で受ける）。
    """
    from matplotlib.colors import to_rgb

    return tuple(float(c) for c in to_rgb(colour))


# ★かつてここに `_colour_texture` と `_place_button`（テクスチャ版ボタンの
#   置き場を決めるもの）があったが、**2026-08-24 に消した**。
#   `vtkTexturedButtonRepresentation2D` を並べるとテクスチャが GPU の上限を
#   超え、実行ログに「Hardware does not support the number of textures
#   defined」が 391 回出た（シェーダも組めなくなる）。
#   いまは文字の背景と枠（`ControlPanel._label`）＋クリック座標の当たり判定
#   （`_hit_area`）で作ってあり、テクスチャを 1 枚も使わない。


def _default_step(fmt, low, high):
    """スピンボックスの 1 段。表示の桁と範囲の広さから決める。"""
    digits = 0
    if "." in fmt:
        try:
            digits = int(fmt.split(".")[1][0])
        except (IndexError, ValueError):
            digits = 0
    unit = 10.0 ** (-digits)                    # 表示できる最小の刻み
    span = abs(high - low)
    rough = span / 100.0 if span else unit
    if rough <= unit:
        return unit
    # 1 / 2 / 5 × 10^n の形に丸める（0.01 → 0.02 → 0.05 → 0.1 …）
    import math
    exponent = math.floor(math.log10(rough))
    base = rough / (10.0 ** exponent)
    nice = 1.0 if base < 1.5 else (2.0 if base < 3.5 else 5.0)
    return max(unit, nice * (10.0 ** exponent))


class PanelItem:
    """パネルに積んだ要素 1 つ。**高さ・置き直し方・隠し方**を持つ。

    ウィンドウの大きさが変わったら `place(y)` を呼び直すだけで並べ直せる。
    `show(False)` は**ページ送りで範囲外に出たとき**に呼ぶ（2026-08-21）。
    文字はレンダラのビューポートで切られるが、スライダは
    ウィンドウ全体の座標で描かれるので**隠さないと 3D の上に残る**。
    """

    # `group` … どのタブの欄か（None なら共通。2026-08-24）
    __slots__ = ("height", "place", "show", "group")

    def __init__(self, height, place, show=None, group=None):
        self.height = height
        self.place = place
        self.show = show
        self.group = group


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
    # スライダのバーぶんの高さ [px]（見出しは別に 1 行取る）。
    # 操作が増えてパネルに入りきらなくなったので 34 → 28 に詰めた（2026-08-19）
    SLIDER = 28
    SPIN = 22           # スピンボックス 1 行の高さ [px]（▲▼ を 2 段に置くため）
    FIELD = 84          # 数値の枠の幅 [px]
    ARROW = 15          # ▲▼ の幅 [px]（枠の右にくっつける）
    ARROW_FONT = 6      # ▲▼ の文字の大きさ（半分の高さに収める）
    BUTTON_FACE = "#3a4150"
    # ★数値の枠は**入力欄に見えるように**する（2026-08-24 ユーザー指摘
    #   「今どちらかと言うとボタンになってないですか？」）。
    #   パネルより暗い地＋明るい細枠＋左寄せの数字＝入力欄、
    #   ▲▼ はパネルより明るい地＝押すもの、と見た目を分ける
    FIELD_FACE = "#101318"      # 数値の枠（入力欄。ほぼ黒）
    FIELD_EDIT = "#18283f"      # 打ち込んでいる最中（青みを入れて分かるように）
    FIELD_TEXT = "#eef2f8"      # 枠の中の数字（明るく）
    FRAME_COLOR = "#7f8794"     # 枠線（入力欄らしく、はっきり見える明るさ）
    ARROW_FACE = "#4a5261"      # ▲▼ の地（パネルより明るくして押すものに見せる）
    FIELD_PAD = 7               # 枠の左の余白 [px]（数字の始まる位置）
    TAB_HEIGHT = 20             # タブの見出しの高さ [px]
    TAB_ACTIVE = "#2f6f9f"      # 選んでいるタブ
    TAB_IDLE = "#2b303a"        # 選んでいないタブ
    LABEL_FONT = 8      # レイヤ名など、横に長くなりがちな文字
    LINE = 18           # font_size 9 のときの 1 行 [px]（外から参照される既定値）
    WHEEL = 60          # ホイール 1 段で送る量 [px]（3 行ぶん）

    def __init__(self, plotter, font=None, margin=14, width_ratio=PANEL_RATIO):
        self.plotter = plotter
        self.font = font
        self.margin = margin
        self.ratio = width_ratio
        self.items = []
        self._widgets = []                 # 参照を残さないと GC で消える
        self.controls = []                 # 数値入力ダイアログに出すスライダ
        # タブ（2026-08-24）。`_group` は「いま作っている欄がどのタブのものか」、
        # `active_group` は「いま開いているタブ」。None なら全部出す
        self._group = None
        self.active_group = None
        # ★**枠の中で数字を打つ**ための状態（2026-08-24 ユーザー要望
        #   「任意入力する場合は別ウィンドウが開くのではなく、その中で完結させたい」）。
        #   `_editing` は編集中の欄（`controls` の 1 件）、`_buffer` は打った文字
        self._editing = None
        self._buffer = ""
        self._caret = 0            # 入力位置（`_buffer` の何文字目の前か）
        self._editor_tag = None
        # 押せる四角（テクスチャを使わないボタン。2026-08-24）
        self._hits = []
        self._click_tag = None
        self.scroll = 0.0                  # ページ送りの量 [px]
        self._help_lines = []              # 重ねて出す操作の一覧
        self._help_title = "操作の一覧"
        self._help_note = None
        self._help_actor = None
        self._hint = None                  # 下端に固定する案内（送られない）
        self._scroll_label = "PageUp/PageDown"
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

    # ---- 押せる四角（テクスチャを使わないボタン。2026-08-24）-----------
    #
    # ★`vtkTexturedButtonRepresentation2D` を並べると**テクスチャが GPU の上限を
    #   超える**（実行ログに「Hardware does not support the number of textures
    #   defined」が 391 回出た）。文字の背景と枠で見た目を作り、
    #   **クリックの座標を自分で見て**押されたことにする。

    def _install_clicker(self):
        """クリックの見張りを 1 回だけ仕込む。"""
        if self._click_tag is not None:
            return True
        interactor = getattr(getattr(self.plotter, "iren", None), "interactor", None)
        if interactor is None:
            return False

        def on_click(caller, _event):
            x, y = caller.GetEventPosition()
            if x > self.width:
                return                  # 3D 側のクリックは触らない（回転など）
            for hit in reversed(self._hits):
                if not hit["visible"]:
                    continue
                x0, y0, x1, y1 = hit["rect"]
                if x0 <= x <= x1 and y0 <= y <= y1:
                    # 押した位置も渡す（数値の枠は**どの桁を押したか**で使う）
                    hit["on_click"](x, y)
                    command = caller.GetCommand(self._click_tag)
                    if command is not None:
                        command.AbortFlagOn()   # 回転させない
                    return

        self._click_tag = interactor.AddObserver("LeftButtonPressEvent",
                                                 on_click, 10.0)
        return True

    def _hit_area(self, on_click):
        """押せる四角を 1 つ登録する。戻り値で位置と表示を差し替える。"""
        self._install_clicker()
        hit = {"rect": (0, 0, -1, -1), "on_click": on_click, "visible": True}
        self._hits.append(hit)
        return hit

    def _label(self, text="", size=None, background=None, frame=False,
               colour=TEXT_COLOR):
        """文字。`background` を渡すと**背景と枠**を付けてボタンらしくする。"""
        self._panel()
        actor = self.plotter.add_text(text or " ", position=(0, 0),
                                     font_size=size or self.LABEL_FONT,
                                     color=colour, font_file=self.font)
        if background is not None:
            prop = actor.GetTextProperty()
            prop.SetBackgroundColor(_rgb(background))
            prop.SetBackgroundOpacity(1.0)
            if frame:
                prop.SetFrame(True)
                prop.SetFrameColor(_rgb(self.FRAME_COLOR))
                prop.SetFrameWidth(1)
        self._model()
        return actor

    # ---- 枠の中で数字を打つ（2026-08-24）------------------------------
    #
    # VTK には文字入力の部品が無いので、**キー入力を横取りして自前で組む**。
    # `KeyPressEvent` を**高い優先度**で見張り、編集中だけ受け取って
    # `AbortFlagOn()` で先へ流さない（流すと `q` で閉じたり数字でレイヤが
    # 切り替わったりしてしまう）。編集していないときは素通しする。

    # 打ち込みに使える文字（数字・小数点・符号）
    EDIT_CHARS = "0123456789.-+"

    def _install_editor(self):
        """キー入力の見張りを 1 回だけ仕込む。"""
        if self._editor_tag is not None:
            return True
        interactor = getattr(getattr(self.plotter, "iren", None), "interactor", None)
        if interactor is None:
            return False                # off_screen などで interactor が無い

        def on_key(caller, _event):
            if self._editing is None:
                return                  # 編集していないので何もしない（素通し）
            keep = self._editor_key(caller.GetKeySym(), caller.GetKeyCode())
            if keep:
                # ★編集に使ったキーは**先へ流さない**
                command = caller.GetCommand(self._editor_tag)
                if command is not None:
                    command.AbortFlagOn()

        # 優先度を上げて pyvista のキー処理より先に受け取る
        self._editor_tag = interactor.AddObserver("KeyPressEvent", on_key, 10.0)
        return True

    def _editor_key(self, keysym, keycode):
        """編集中のキー 1 つを処理する。戻り値 True なら先へ流さない。

        ★**編集に関係ないキーが来たら編集をやめて素通しする**。
        こうしておかないと、打ち込んでいる最中に `q` を押しても閉じられない。
        """
        if keysym in ("Return", "KP_Enter"):
            self.commit_edit()
            return True
        if keysym == "Escape":
            self.cancel_edit()
            return True
        if keysym == "BackSpace":
            # キャレットの**前**の 1 文字を消す
            if self._caret > 0:
                self._buffer = (self._buffer[:self._caret - 1]
                                + self._buffer[self._caret:])
                self._caret -= 1
            self._draw_buffer()
            return True
        if keysym == "Delete":
            # キャレットの**位置**の 1 文字を消す
            self._buffer = (self._buffer[:self._caret]
                            + self._buffer[self._caret + 1:])
            self._draw_buffer()
            return True
        if keysym == "Left":
            self._caret = max(0, self._caret - 1)
            self._draw_buffer()
            return True
        if keysym == "Right":
            self._caret = min(len(self._buffer), self._caret + 1)
            self._draw_buffer()
            return True
        if keysym == "Home":
            self._caret = 0
            self._draw_buffer()
            return True
        if keysym == "End":
            self._caret = len(self._buffer)
            self._draw_buffer()
            return True
        if keycode and keycode in self.EDIT_CHARS:
            # ★キャレットの位置に**挿し込む**（全部を置き換えたりしない）
            self._buffer = (self._buffer[:self._caret] + keycode
                            + self._buffer[self._caret:])
            self._caret += 1
            self._draw_buffer()
            return True
        # 関係ないキー → 編集をやめて、そのキーは本来の役目に回す
        self.cancel_edit()
        return False

    def start_edit(self, control, caret=None):
        """★枠を押したときに呼ばれる。**その枠の中で**打ち込みを始める。

        `caret` は入力位置（文字数）。枠のどこを押したかから決めて渡す。
        """
        if self._editing is not None and self._editing is not control:
            self.commit_edit()          # 別の欄へ移ったら、それまでの分を確定
        if not self._install_editor():
            # 見張りを仕込めない環境では、従来どおり別窓で受ける（保険）
            return self.open_value_input(only=control)
        self._editing = control
        # ★★**押した時点では、いまの値をそのまま残す**（2026-08-24 ユーザー指摘
        #   「入力しようとしたら、数字が消えてしまいます」）
        self._buffer = (control["format"] % control["value"]).strip()
        # ★**押した桁にキャレットを立てる**（同日「0.45 の 4 の部分を消して
        #   3 に書き換えたい」）。位置が分からなければ末尾に置く
        self._caret = (len(self._buffer) if caret is None
                       else max(0, min(len(self._buffer), int(caret))))
        control["highlight"](True)
        self._draw_buffer()
        return control

    def commit_edit(self):
        """打ち込んだ値を確定する（Enter）。"""
        control, text = self._editing, self._buffer
        self._editing, self._buffer = None, ""
        if control is None:
            return None
        control["highlight"](False)
        try:
            value = float(text)
        except ValueError:
            control["show"]()           # 空・数字でない → もとの値に戻す
            self.plotter.render()
            return None
        control["set"](value)
        return value

    def cancel_edit(self):
        """打ち込みをやめる（Escape・関係ないキー・別の欄へ移ったとき）。"""
        control = self._editing
        self._editing, self._buffer = None, ""
        if control is None:
            return None
        control["highlight"](False)
        control["show"]()
        self.plotter.render()
        return control

    def _draw_buffer(self):
        """打っている途中の文字を枠に出す。**入力位置に `|` を挟む。**"""
        if self._editing is None:
            return
        caret = max(0, min(len(self._buffer), self._caret))
        self._editing["draw"](self._buffer[:caret] + "|" + self._buffer[caret:])
        self.plotter.render()

    def open_value_input(self, only=None):
        """数字を打ち込む窓を出す。`only` を渡すとその欄だけにする。

        ★欄ごとの `[123]` ボタンからは `only` 付きで呼ばれる（2026-08-24）。
        `t` キーからは全部まとめて出す（従来どおり）。
        """
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
        targets = [only] if only is not None else self.controls
        for row, control in enumerate(targets, start=1):
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

    def _pin(self, widget):
        """★ウィジェットを**パネル側のレンダラに固定する**。

        VTK のウィジェットは `SetEnabled(1)` のときに
        `FindPokedRenderer`（マウスの下のレンダラ）で置き場を決める。
        固定しないと、送って隠す・出すを繰り返すうちに**3D 側のレンダラにも
        登録され、同じものが 2 つ描かれる**（2026-08-21 ユーザー指摘。
        吸音材のチェックボックスと方向のスライダがモデル側に飛び出していた）。
        """
        try:
            panel = self.plotter.renderers[0]
            widget.SetDefaultRenderer(panel)
            widget.SetCurrentRenderer(panel)
        except Exception:
            pass
        return widget

    def _add(self, height, place, show=None):
        item = PanelItem(height, place, show)
        # ★いま開いているタブ（`begin_group`）のものとして覚える。
        #   `group` が None のものは**どのタブでも出る**（共通の欄）
        item.group = self._group
        self.items.append(item)
        return item

    # ---- タブ（2026-08-24 ユーザー要望）--------------------------------
    #
    # > 音線・音粒子・虚音源、すべて一画面で表示されていて、
    # > どれが何の設定かわからない
    #
    # 欄に「どのタブのものか」の印を付け、選んでいるタブのものだけ出す。
    # 印の無いものは共通（レイヤ表示・不透明度など）で、いつでも出る。

    def begin_group(self, name):
        """ここから下に作る欄を `name` のタブのものにする。"""
        self._group = name
        return name

    def end_group(self):
        """共通の欄に戻す。"""
        self._group = None

    def show_group(self, name, render=True):
        """`name` のタブを開く（他のタブの欄は隠す）。"""
        self.active_group = name
        self.scroll = 0.0
        self.relayout(render=render)
        return name

    def groups(self):
        """使われているタブの名前（作った順）。"""
        seen = []
        for item in self.items:
            if item.group and item.group not in seen:
                seen.append(item.group)
        return seen

    def relayout(self, render=False):
        """いまのウィンドウの大きさに合わせて、上から順に積み直す。

        ★**ページ送りの分だけ下へずらす**（`self.scroll`）。
        ずらして範囲外に出た要素は `show(False)` で隠す
        （スライダはウィンドウ全体の座標で描かれるので、隠さないと 3D に残る）。
        """
        self._measure()
        self.scroll = min(max(0.0, self.scroll), self.max_scroll())
        top, bottom = self.height - self.margin, self.margin + self._hint_height()
        y = top + self.scroll
        for item in self.items:
            # ★選んでいないタブの欄は隠し、**画面の外へ追い出す**（高さも取らない）。
            #   `Off()` だけだと稀に四角が残ることがあったので、
            #   置き場も範囲外にしておく（見えていた四角が消える）
            if not self._in_active_group(item):
                item.place(-10000.0)
                if item.show is not None:
                    item.show(False)
                continue
            y -= item.height
            item.place(y)
            if item.show is not None:
                item.show(y >= bottom - 1 and y + item.height <= top + 1)
        self._update_hint()
        if render:
            self.plotter.render()

    # ---- ページ送り（2026-08-21）----------------------------------------
    #
    # ★以前は「入りきらない分は載せる内容を削る」方針だったが、実案件では
    #   レイヤ 9 種 ＋ 材料 19 種で 1280 px でも 450 px、1920 px でも 230 px
    #   あふれた（ユーザー指摘）。**削り切れないので送れるようにした。**
    #   前に試して取り下げたのは「置き直しを飛ばす」やり方で、
    #   古い位置に残ったウィジェットが文字に重なっていた。いまは必ず置き直し、
    #   範囲外は隠す。

    def enable_scroll(self, back="Prior", forward="Next",
                      label="ホイール / PageUp・PageDown"):
        """入りきらないときにスクロールできるようにする。

        ★**普通のマウスホイールで送る**（2026-08-21 ユーザー要望）。
        パネルの上にカーソルがあるときだけパネルを送り、
        3D の上なら今までどおり拡大縮小になる。キーボードも残す。

        ホイールは 3D 側の「近づく・遠ざかる」にも割り当てられているが、
        **パネル側のレンダラには 2D の文字とウィジェットしか置いていない**ので、
        そちらのカメラが動いても見た目は変わらない（だから止めなくてよい）。
        """
        self._scroll_label = label
        self.plotter.add_key_event(back, lambda: self.scroll_by(-1))
        self.plotter.add_key_event(forward, lambda: self.scroll_by(+1))
        self._watch_wheel()
        self._ensure_hint()
        self.relayout()

    def _watch_wheel(self):
        """マウスホイールを見張る（パネルの上なら送る）。"""
        def wheel(step):
            def handler(*_args):
                if self.max_scroll() <= 0:
                    return
                try:
                    x = self.plotter.iren.interactor.GetEventPosition()[0]
                except Exception:
                    return
                if x <= self.width:        # パネルの上にカーソルがあるときだけ
                    self.scroll_pixels(step * self.WHEEL)
            return handler

        try:
            self.plotter.iren.add_observer("MouseWheelForwardEvent", wheel(-1))
            self.plotter.iren.add_observer("MouseWheelBackwardEvent", wheel(+1))
        except Exception:
            pass        # ホイールが取れない環境でもキーで送れる

    def help_window(self, lines, key="F1", title="操作の一覧", note=None):
        """操作の一覧を **3D の上に重ねて**出す（`key` で表示/非表示）。

        ★はじめは tkinter の別ウィンドウにしたが、**開いている間は元の画面を
        閉じられない**（tkinter が入力を握る）とユーザー指摘を受けて作り直した。
        いまは同じウィンドウの中に重ねるだけなので、取り合いが起きない。

        パネルには「F1 操作の一覧」の 1 行だけ置く。
        """
        self._help_lines = list(lines)
        self._help_title = title
        self._help_note = note
        self._model()
        head = [f"― {title} ―"] + ([note] if note else []) + [""]
        # 文字は小さめ（10 だと大きいというユーザー指摘 2026-08-21）
        self._help_actor = self.plotter.add_text(
            chr(10).join(head + list(lines)), position=(24, 24), font_size=8,
            color="#ffe9a8", font_file=self.font)
        self._help_actor.SetVisibility(False)
        self.plotter.add_key_event(key, self.toggle_help)
        self.text(f"{key} 操作の一覧（重ねて表示）", size=8, color="#ffd166")
        return self

    def toggle_help(self):
        """操作の一覧の表示を切り替える。"""
        if self._help_actor is None:
            return
        self._help_actor.SetVisibility(not self._help_actor.GetVisibility())
        self.plotter.render()

    def _in_active_group(self, item):
        """この欄をいま出すか（共通の欄はいつでも出す）。"""
        group = getattr(item, "group", None)
        return (group is None or self.active_group is None
                or group == self.active_group)

    def max_scroll(self):
        """送れる最大量 [px]（入りきっていれば 0）。"""
        return max(0.0, self.content_height()
                   - (self.height - 2 * self.margin - self._hint_height()))

    def scroll_by(self, pages):
        """`pages` 枚ぶん送る（1 枚 = 見えている高さの 8 割）。"""
        step = max(40.0, (self.height - 2 * self.margin) * 0.8)
        return self.scroll_pixels(pages * step)

    def scroll_pixels(self, pixels):
        """`pixels` だけ送る（ホイール 1 段は `WHEEL` px）。"""
        limit = self.max_scroll()
        scroll = min(max(0.0, self.scroll + pixels), limit)
        if scroll == self.scroll:
            return self.scroll
        self.scroll = scroll
        self.relayout(render=True)
        return self.scroll

    def _ensure_hint(self):
        """パネルの下端に固定する案内の文字（**積まないので送られない**）。"""
        if self._hint is not None:
            return self._hint
        self._panel()
        self._hint = self.plotter.add_text(" ", position=(self.margin, self.margin),
                                           font_size=8, color="#ffd166",
                                           font_file=self.font)
        self._model()
        return self._hint

    def _hint_height(self):
        return self.line_height(8) + 4 if self._hint is not None else 0

    def _update_hint(self):
        if self._hint is None:
            return
        limit = self.max_scroll()
        if limit <= 0:
            self._hint.SetInput(" ")
        else:
            page = int(self.scroll / max(1.0, limit) * 100.0)
            self._hint.SetInput(f"▲▼ {self._scroll_label} で送る（{page}%）")
        self._hint.SetDisplayPosition(self.margin, self.margin)

    def content_height(self):
        """並べた要素の高さの合計 [px]。**いま開いているタブのぶんだけ**数える。"""
        return sum(item.height for item in self.items
                   if self._in_active_group(item))

    def hidden_height(self):
        """入りきらずに下へはみ出している高さ [px]（入りきっていれば 0）。

        **ウィンドウを縦に広げれば減る**（パネルの高さはウィンドウの高さそのもの）。
        入りきらない場合は `enable_scroll()` でページ送りできる（2026-08-21）。
        """
        return self.max_scroll()

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

        self._add(height * lines, place,
                  lambda visible, a=actor: a.SetVisibility(bool(visible)))
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
                  lambda y, a=actor: a.SetDisplayPosition(self.margin, int(y)),
                  lambda visible, a=actor: a.SetVisibility(bool(visible)))
        return actor

    def screen_title(self, message, note="設定・操作"):
        """パネルのいちばん上に置く見出し。

        ★**ここが「設定の面」だと分かるように**、3D 側と違う色にする
        （2026-08-21 ユーザー要望。背景の色分けと合わせて別画面に見せる）。
        """
        self.text(message, size=11, color=PANEL_TITLE_COLOR)
        if note:
            self.text(note, size=8, color="#6b7385")
        return self

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
        self._pin(widget)
        self._widgets.append(widget)

        def place(y, widget=widget, text=text, label=label):
            left = self.margin + self.CHECK + 8
            representation = widget.GetRepresentation()
            representation.SetPlaceFactor(1.0)
            representation.PlaceWidget([self.margin, self.margin + self.CHECK,
                                        y, y + self.CHECK, 0.0, 0.0])
            # ★置き場所を変えただけでは描き直されない（内部に計算結果を持っている）。
            #   忘れると**送ったときに文字だけ動いて四角が取り残される**
            #   （2026-08-21 ユーザー指摘。スライダ側には元から入れてあった）
            representation.Modified()
            representation.BuildRepresentation()
            text.SetInput(self.fit(label, self.width - left - self.margin))
            text.SetDisplayPosition(left, int(y) + 3)

        # 範囲外に出たら**チェックボックスも文字も隠す**（ページ送りのため）
        shown = {"on": True}

        def show(visible, widget=widget, text=text):
            text.SetVisibility(bool(visible))
            if bool(visible) == shown["on"]:
                return          # 変わらないなら触らない（登録し直しを減らす）
            shown["on"] = bool(visible)
            widget.On() if visible else widget.Off()

        self._add(self.CHECK + 4, place, show)
        return widget

    def slider(self, title, value_range, value, callback, fmt="%.2f", step=None):
        """数値の設定。★**中身はスピンボックス**（2026-08-24 ユーザー要望）。

        > 数字の設定がスライダーになっていますが、
        > スピンボックス（手入力機能もあり）の方がよいです。

        呼んでいる場所が多いので**名前は `slider` のまま**にしてある
        （呼び出し側を触らずに全画面が切り替わる）。
        つまみのスライダが要るときは `slider_widget()` を直接呼ぶ。
        """
        return self.spinbox(title, value_range, value, callback, fmt=fmt,
                            step=step)

    def spinbox(self, title, value_range, value, callback, fmt="%.2f", step=None):
        """数値の欄。**Excel と同じ形**（枠の右に ▲▼、枠を押せばその場で手入力）。

        ★スライダをやめた理由（ユーザー要望 2026-08-24）：つまみは
        「だいたいの値」を選ぶのは速いが、**狙った値にできない**
        （反射回数を 7 回にしたい、音線を 137 本にしたい、など）。

        ★見た目（同日・画像つきの指摘）：「数値とボタンが離れて良くわからない」
        → **枠の右に ▲▼ をくっつけ、枠を押せばその場で打ち込める**。

            見出し            [   1.00 ][▲]
                                       [▼]

        ★★中身は**テクスチャを使わない**（同日、実行ログで
        「Hardware does not support the number of textures defined」が
        391 回出た）。見た目は文字の背景と枠、押した判定はクリック座標。

        `step` を省いたときの 1 段は、表示の桁と範囲から決める。
        """
        low, high = float(value_range[0]), float(value_range[1])
        if high < low:
            low, high = high, low
        step = float(step) if step else _default_step(fmt, low, high)
        state = {"value": float(min(max(value, low), high))}

        def show_value(v=None):
            if v is not None:
                state["value"] = float(v)
            set_actor_text(field, f" {fmt % state['value']} ")

        def apply(v, notify=True):
            v = float(min(max(v, low), high))
            if abs(v - state["value"]) < 1e-12 and notify:
                return                      # 端で押し続けたときに呼び直さない
            state["value"] = v
            control["value"] = v
            show_value()
            self._model()
            if notify:
                callback(v)
            self.plotter.render()

        name = self._label(size=self.LABEL_FONT)
        field = self._label(background=self.FIELD_FACE, frame=True,
                            colour=self.FIELD_TEXT)
        up = self._label("▲", size=self.ARROW_FONT, background=self.ARROW_FACE)
        down = self._label("▼", size=self.ARROW_FONT, background=self.ARROW_FACE)

        def begin_edit(x=None, y=None):
            """★**押した桁**から編集を始める（2026-08-24 ユーザー要望
            「0.45 の 4 の部分を消して 3 に書き換えたい」）。

            文字の幅を測って、押した x がどの文字の上かを出す。
            """
            caret = None
            if x is not None:
                text = (fmt % state["value"]).strip()
                left = field_hit["rect"][0] + self.FIELD_PAD
                caret = _caret_from_x(field, text, left, x, self.plotter.renderer)
            self.start_edit(control, caret=caret)

        # 押せる四角（枠＝手入力、▲▼＝1 段ずつ）
        field_hit = self._hit_area(begin_edit)
        up_hit = self._hit_area(lambda *_: apply(state["value"] + step))
        down_hit = self._hit_area(lambda *_: apply(state["value"] - step))

        def draw_text(text):
            """打ち込んでいる途中の文字をそのまま枠に出す。"""
            set_actor_text(field, f" {text} ")

        def highlight(flag):
            """打ち込み中は枠を明るくする（どこに入るか分かるように）。"""
            prop = field.GetTextProperty()
            prop.SetBackgroundColor(_rgb(self.FIELD_EDIT if flag
                                         else self.FIELD_FACE))

        control = {"label": title, "range": (low, high), "value": state["value"],
                   "format": fmt, "step": step,
                   "set": lambda v: apply(v), "show": show_value,
                   "draw": draw_text, "highlight": highlight,
                   "widget": None}
        self.controls.append(control)
        show_value()

        def place(y):
            right = self.width - self.margin
            arrow_x = right - self.ARROW
            field_x = max(self.margin + 40, arrow_x - self.FIELD)
            half = int(self.SPIN / 2)
            name.SetDisplayPosition(self.margin, int(y) + 5)
            name.SetInput(self.fit(title, field_x - self.margin - 6))
            field.SetDisplayPosition(int(field_x), int(y) + 4)
            up.SetDisplayPosition(int(arrow_x) + 3, int(y) + half + 1)
            down.SetDisplayPosition(int(arrow_x) + 3, int(y) + 1)
            field_hit["rect"] = (field_x, y, arrow_x - 1, y + self.SPIN)
            up_hit["rect"] = (arrow_x, y + half, right, y + self.SPIN)
            down_hit["rect"] = (arrow_x, y, right, y + half)

        def show(visible):
            for actor in (name, field, up, down):
                actor.SetVisibility(bool(visible))
            for hit in (field_hit, up_hit, down_hit):
                hit["visible"] = bool(visible)

        self._add(self.SPIN + 6, place, show)
        return control

    def tab_strip(self, labels, on_select, active=0):
        """★**タブの見出しを横に並べる**（2026-08-24 ユーザー要望）。

        > 音線・音粒子・虚音源、すべて一画面で表示されていて、どれが何の設定か
        > わからない。左の欄もタブで切り替える形式にしてはどうでしょう？
        > 画面上側に切替える用のタブを用意して…

        押されたタブを濃く、それ以外を薄く塗る。`on_select(タブ名)` が呼ばれる。
        戻り値は `select(タブ名)`。**Tab キーから切り替えたときも見出しを合わせる**
        ために使う（キーとタブで見た目がずれないように）。
        """
        labels = list(labels)
        state = {"active": labels[active] if labels else None}
        tabs = []
        for label in labels:
            actor = self._label(label, background=self.TAB_IDLE, frame=True)
            tabs.append({"label": label, "actor": actor,
                         "hit": self._hit_area(
                             lambda *_, name=label: choose(name))})

        def paint():
            for tab in tabs:
                colour = (self.TAB_ACTIVE if tab["label"] == state["active"]
                          else self.TAB_IDLE)
                tab["actor"].GetTextProperty().SetBackgroundColor(_rgb(colour))

        def choose(name):
            if name == state["active"]:
                return
            state["active"] = name
            paint()
            on_select(name)

        def place(y):
            if not tabs:
                return
            span = (self.width - 2 * self.margin) / len(tabs)
            for index, tab in enumerate(tabs):
                x0 = self.margin + index * span
                tab["actor"].SetInput(f" {self.fit(tab['label'], span - 14)} ")
                tab["actor"].SetDisplayPosition(int(x0), int(y) + 4)
                tab["hit"]["rect"] = (x0, y, x0 + span - 3, y + self.TAB_HEIGHT)

        def show(visible):
            for tab in tabs:
                tab["actor"].SetVisibility(bool(visible))
                tab["hit"]["visible"] = bool(visible)

        self._add(self.TAB_HEIGHT + 8, place, show)
        paint()

        def select(name):
            """外から切り替える（Tab キーなど）。見た目だけ合わせる。"""
            if name in labels:
                state["active"] = name
                paint()
            return name

        return select

    def button(self, label, on_click, colour=None):
        """★**押すだけのボタン**（2026-08-24 ユーザー指摘）。

        > 受音点の向きを一括で向けるボタンがチェックボックスで複数選択可に
        > なってます。チェックボックス・ラジオボタンでは無く、ボタンとして、
        > そのボタンを押せばその向きになる、というだけで良いです。
        """
        actor = self._label(label, background=colour or self.BUTTON_FACE,
                            frame=True)
        hit = self._hit_area(lambda *_: on_click())

        def place(y):
            width = self.width - 2 * self.margin
            actor.SetInput(f" {self.fit(label, width - 14)} ")
            actor.SetDisplayPosition(self.margin, int(y) + 4)
            hit["rect"] = (self.margin, y, self.margin + width, y + self.SPIN)

        def show(visible):
            actor.SetVisibility(bool(visible))
            hit["visible"] = bool(visible)

        self._add(self.SPIN + 6, place, show)
        return actor

    def slider_widget(self, title, value_range, value, callback, fmt="%.2f"):
        """横向きのスライダ（**参照実装として残してある**。いまは使っていない）。

        2026-08-24 に数値の設定はスピンボックスへ移した（`spinbox`）。
        つまみで連続的に動かしたい場面が出たらこちらを使う。

        見出しと値は**自前の文字**で出す。

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
        self._pin(widget)
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
                   "format": fmt, "set": set_value, "widget": widget,
                   # 外から値を入れ直したときに**見出しの数字も直せる**ようにしておく
                   # （つまみだけ動かして数字が古いままだと、どちらが本当か分からない）
                   "show": show_value}
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

        shown = {"on": True}

        def show(visible, widget=widget):
            # スライダは**ウィンドウ全体**の座標で描くので、隠さないと 3D に残る
            if bool(visible) == shown["on"]:
                return
            shown["on"] = bool(visible)
            widget.On() if visible else widget.Off()

        self._add(self.SLIDER, place, show)
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
    "images":     ("虚音源",         "geosim - image sources"),
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


# 画面に出す知らせの色。よい知らせ・待たせる知らせ・悪い知らせで分ける
NOTICE_COLORS = {"ok": "#7ee787", "busy": "#ffd166", "error": "#ff7b72"}
NOTICE_SECONDS = 3.0


def notice(plotter, message, kind="ok", seconds=NOTICE_SECONDS):
    """★**画面の真ん中下に一時的な知らせを出す**（2026-08-24 ユーザー要望）。

    > 画面の画像保存ボタンを押した時に、保存しました。みたいな表示をして欲しい
    > 動画保存時は、動画保存中みたいな表示が欲しい

    それまでは端末に `print` していただけで、画面を見ている人には
    保存できたのか分からなかった。

    `seconds` 秒たったら自分で消える（`add_timer_event`）。
    消えるのを待たずに次の知らせが来たら、そのまま置き換わる。
    `seconds=None` にすると消えない（`clear_notice()` で消す。動画の書き出し中など、
    終わるまで出しっぱなしにしたいとき）。
    """
    text = str(message)
    colour = NOTICE_COLORS.get(kind, NOTICE_COLORS["ok"])
    try:
        plotter.add_text(text, name="geosim_notice", position="lower_edge",
                         font_size=11, color=colour, font=japanese_font(),
                         shadow=True)
    except Exception:
        # 文字が出せなくても本体は動かす（フォントが無い環境など）
        print(f"[view] {text}")
        return None
    plotter.render()
    if seconds:
        _hide_notice_later(plotter, seconds)
    return text


def clear_notice(plotter):
    """出している知らせを消す。"""
    try:
        plotter.remove_actor("geosim_notice", render=False)
        plotter.render()
    except Exception:
        pass


def _hide_notice_later(plotter, seconds):
    """`seconds` 秒後に知らせを消す。

    ★VTK のタイマーを使う（`time.sleep` で待つと画面が固まる）。
    タイマーが使えない場面（off_screen など）では出したままにする
    ——消えないほうがましなので、黙って諦める。
    """
    state = {"count": 0}

    def tick(step=None):
        state["count"] += 1
        if state["count"] >= 1:
            clear_notice(plotter)

    try:
        plotter.add_timer_event(max_steps=1, duration=int(seconds * 1000),
                                callback=tick)
    except Exception:
        pass


def add_screenshot_key(plotter, folder, stem, key="p"):
    """`p` で「いま画面に出ているとおり」を PNG に保存する。

    **左のパネルごと写す。** どのレイヤを出していたか・不透明度がいくつだったかも
    一緒に残るほうが、あとで見返したときに条件が分かるため。

    `stem` は文字列でも、呼ぶたびにファイル名の芯を返す関数でもよい
    （音粒子は時刻をファイル名に入れたいので関数を渡す）。

    ★保存したら**画面にも知らせる**（2026-08-24 ユーザー要望）。
    それまでは端末に print するだけで、画面を見ている人には分からなかった。
    """
    def save():
        name = stem() if callable(stem) else stem
        path = next_free_path(folder, name)
        try:
            plotter.screenshot(path)
        except Exception as e:
            print(f"[view] 画像を保存できませんでした: {type(e).__name__}: {e}")
            notice(plotter, f"画像を保存できませんでした（{type(e).__name__}）",
                   kind="error")
            return None
        print(f"[view] 画像を保存しました: {path}")
        notice(plotter, f"画像を保存しました: {os.path.basename(path)}")
        return path

    plotter.add_key_event(key, save)
    return save


def release_window(plotter, *holders):
    """ウィンドウを閉じ、**VTK の持ち物を明示的に手放す**。

    【なぜ要るか】2026-08-19
    ウィンドウを閉じたあと、プロセス終了時に **segfault** した
    （`Could not set shader program` / `Error binding ndCoords to VAO` が
    大量に出たあとに落ちる）。VTK の OpenGL の資源は「文脈（コンテキスト）」が
    生きているうちに解放しないといけないが、Python の後片付けは順序が決まっていない。
    こちら（`RayDisplay.actors` / `ParticleAnimation.actor` / パネルのウィジェット）が
    参照を握ったままだと、**文脈が消えたあとに解放されて落ちる**。

    そこで閉じる前に、握っている参照をこちらから外してから閉じる。
    `holders` には actor を抱えているオブジェクトを渡す（属性を空にする）。

    ★**再現はできていない。** 実際に落ちたのは 5 分ほど操作したあとの終了時で、
      手元では同じ手順（対話ループ・スライダ操作・絞り込み）を繰り返しても
      再現しなかった。GL ドライバの状態にも依るらしい。
      これは「落ちる余地を減らす」ための処置で、根治の確認は取れていない。
    """
    for holder in holders:
        if holder is None:
            continue
        for name in ("actors", "actor", "marker", "label", "arrows", "_widgets",
                     "items", "controls"):
            if hasattr(holder, name):
                try:
                    setattr(holder, name, [] if name in ("actors", "_widgets",
                                                         "items", "controls")
                            else None)
                except Exception:
                    pass
    try:
        plotter.clear()
    except Exception:
        pass
    try:
        plotter.close()
    except Exception:
        pass
    try:
        plotter.deep_clean()
    except Exception:
        pass


def close_all():
    """開いているウィンドウを全部片付ける（アプリを終わる直前に呼ぶ）。"""
    try:
        pv.close_all()
    except Exception:
        pass


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
        # ★★`set_background` は既定で**全レンダラに掛かる**（`all_renderers=True`）。
        #   片方ずつ変えるには `all_renderers=False` が要る。
        #   これを忘れて後の呼び出しでパネル側が上書きされ、
        #   色を変えたつもりで**何も変わっていなかった**（2026-08-21 に気づいた）
        plotter.subplot(0, 0)
        # 設定の面。3D 側と**別画面に見えるように**平らで暗い色にする
        plotter.set_background(PANEL_BG, all_renderers=False)
        try:
            # 境目に細い線を入れる（枠なので 4 辺に付くが、暗い色なので目立たない）
            plotter.renderers[0].add_border(color=PANEL_EDGE, width=2)
        except Exception:
            pass        # 描けなくても中身は使えるので黙って諦める
        plotter.subplot(0, 1)
        plotter.set_background(BG_BOTTOM, top=BG_TOP, all_renderers=False)

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
        panel.screen_title(f"{title}")
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
    # ★スピンボックスになったので、戻り値は VTK のウィジェットではなく
    #   `controls` の 1 件（`{"set": …, "show": …}`）。2026-08-24
    control = panel.spinbox("不透明度", [0.0, 1.0], current_opacity(), apply,
                            step=0.05)
    label = panel.reserve_text(2)
    state["ready"] = True

    def set_slider(value):
        """対象を切り替えたときに、**表示だけ**その対象の値に合わせる。"""
        state["ready"] = False          # 表示だけ更新し、適用はしない
        control["show"](float(value))
        state["ready"] = True

    def next_target():
        state["target"] = (state["target"] + 1) % (len(names) + 1)
        set_slider(current_opacity())
        refresh_label()
        plotter.render()

    def toggle_model():
        state["visible"] = not state["visible"]
        for name in names:
            # そのレイヤの持ち物を全部（面・輪郭・矢印）まとめて消す
            for key in ("face", "edge", "arrow"):
                actor = layers[name].get(key)
                if actor is not None:
                    actor.SetVisibility(state["visible"])
        plotter.render()

    plotter.add_key_event(target_key, next_target)
    plotter.add_key_event("m", toggle_model)
    refresh_label()
    return control


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
