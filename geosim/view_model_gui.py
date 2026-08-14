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
  w / s           ワイヤフレーム / 面（VTK の既定キー）
  r               視点リセット、q でウィンドウを閉じる

使い方:
    cd geosim
    python view_model_gui.py ..\\test2.dxf
    python view_model_gui.py ..\\test.dxf --absorption ..\\absorption.csv
    python view_model_gui.py ..\\test.dxf --screenshot shot.png   # 画像だけ書き出す
"""

import argparse
import os

import numpy as np
import pyvista as pv

import read_dxffile as rd
from view_model import LAYER_PALETTE, _hex_to_rgb

# 法線の裏側の色。HTML 版と同じ赤にしてある
BACK_COLOR = "#C24540"
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
                  layer_opacity=None):
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

    plotter = pv.Plotter(window_size=window_size, title=title, off_screen=off_screen)
    plotter.set_background(BG_BOTTOM, top=BG_TOP)

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
                             "diffuse": 0.70},
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

    header = f"{title}\n三角形 {len(mesh)} 枚 / レイヤ {len(layers)}"
    plotter.add_text(header, position="upper_left", font_size=11,
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

    # ---- レイヤの表示切り替え（対話時のみ。off_screen では interactor が無い） ----
    if not off_screen:
        size = 18
        gap = 8
        # 要約テキストが左下を使うので、チェックボックスはヘッダの下から下向きに並べる
        top_y = window_size[1] - 78
        for i, name in enumerate(layers):
            y = top_y - i * (size + gap)
            plotter.add_checkbox_button_widget(
                _visibility_callback(plotter, face_actors[name], arrow_actors[name],
                                     state),
                value=True, position=(14, y), size=size, border_size=2,
                color_on=LAYER_PALETTE[i % len(LAYER_PALETTE)],
                color_off="#454c58", background_color="#2b303a",
            )
            count = model.layer_counts.get(name, 0)
            plotter.add_text(f"{name}  ({count})", position=(14 + size + 10, y + 2),
                             font_size=9, color=TEXT_COLOR, font_file=font)
        plotter.add_text("z/x/c/v 視点  n 法線  w/s 表示  r リセット  q 終了",
                         position=(14, top_y - len(layers) * (size + gap) - 12),
                         font_size=9, color="#7f8794", font_file=font)

    plotter.view_isometric()
    # あとから不透明度や表示を変えられるよう、レイヤごとの actor を Plotter に持たせる。
    # pyvista は新しい公開属性の追加を禁じているので、専用の API を使う
    # （無い版のために private 名へのフォールバックも用意しておく）
    registry = {name: {"face": face_actors[name], "arrow": arrow_actors[name],
                       "colour": LAYER_PALETTE[i % len(LAYER_PALETTE)],
                       "opacity": layer_opacity[name]}
                for i, name in enumerate(layers)}
    try:
        pv.set_new_attribute(plotter, "geosim_layers", registry)
    except AttributeError:
        plotter._geosim_layers = registry
    return plotter


def layer_actors(plotter):
    """`build_plotter` が登録したレイヤ情報を取り出す。"""
    return (getattr(plotter, "geosim_layers", None)
            or getattr(plotter, "_geosim_layers", None))


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


def add_opacity_control(plotter, font=None, pointa=(0.075, 0.30),
                        pointb=(0.075, 0.62), label_position=None,
                        target_key="Tab"):
    """モデルの不透明度を対話的に変えるスライダとキー操作を足す。

    - **Tab**（`target_key`） … 対象を切り替える（すべて → 各レイヤ → すべて …）
    - **スライダ**（左側の縦） … 対象の不透明度を 0〜1 で設定する
    - **m** … モデル全体の表示 ON / OFF

    `target_key` を変えられるようにしてあるのは、Tab を別の用途
    （view_rays の音線↔音粒子の切り替え）に使いたい場面があるため。

    レイヤごとに変えられるようにしてあるのは、
    「壁だけ薄くして中の様子を見る」「床は残す」といった使い方をするため。

    ※ 左上のチェックボックスは**表示のオンオフ**（不透明度とは別）。
    """
    layers = layer_actors(plotter)
    if not layers:
        raise ValueError("build_plotter が作った Plotter を渡してください")

    names = list(layers)
    # ready … スライダを作った直後と、対象切替で値を書き換えるときにコールバックが
    #         走ってしまうので、そのぶんを無視するための旗。
    #         これが無いと、生成時に「全レイヤの平均値」が全レイヤへ適用されてしまい、
    #         layer_opacity で個別に指定した値が消える
    state = {"target": 0, "visible": True, "ready": False}      # target 0 = すべて
    if label_position is None:
        # レイヤのチェックボックスと操作説明の下（左下の座標軸には被らない位置）
        label_position = (14, max(120, plotter.window_size[1] - 230))
    label = plotter.add_text(" ", position=label_position, font_size=9,
                             color=TEXT_COLOR, font_file=font)

    def target_name():
        return "すべて" if state["target"] == 0 else names[state["target"] - 1]

    def current_opacity():
        if state["target"] == 0:
            return float(np.mean([layers[n]["opacity"] for n in names]))
        return layers[target_name()]["opacity"]

    def refresh_label():
        text = (f"不透明度の対象: {target_name()}  ({current_opacity():.2f})\n"
                f"{target_key} 対象切替   m モデル表示 ON/OFF")
        if hasattr(label, "SetInput"):
            label.SetInput(text)
        else:
            label.SetText(0, text)

    def apply(value):
        if not state["ready"]:
            return
        value = float(np.clip(value, 0.0, 1.0))
        targets = names if state["target"] == 0 else [target_name()]
        for name in targets:
            layers[name]["opacity"] = value
            layers[name]["face"].GetProperty().SetOpacity(value)
        refresh_label()
        plotter.render()

    slider = plotter.add_slider_widget(
        apply, [0.0, 1.0], value=current_opacity(), title="opacity",
        pointa=pointa, pointb=pointb, style="modern", fmt="%.2f",
        color=TEXT_COLOR, title_color=TEXT_COLOR)
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
    plotter = build_plotter(model, title=f"{base} モデルビューア",
                            off_screen=screenshot is not None,
                            show_normals=show_normals, opacity=opacity,
                            layer_opacity=layer_opacity)
    if screenshot is None:
        add_opacity_control(plotter, font=japanese_font())
    if screenshot is not None:
        plotter.screenshot(screenshot)
        plotter.close()
        print(f"\n[view_model_gui] 画像を書き出しました: {screenshot}")
    else:
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
