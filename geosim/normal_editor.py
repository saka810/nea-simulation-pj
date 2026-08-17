"""法線の向きを目で確かめて直すウィンドウ（TODO G-8）。

読み込み時の自動補正（`orient_normals='auto'`）に任せきりにせず、
**人が見て確認し、必要なら反転できる**ようにするためのもの。
結果は `project.json` と同じフォルダの `normals.json` に保存され、
次からは自動補正のあとに重ねて適用される（`read_model(flip_faces=...)`）。

なぜ要るか：法線は「音が通る空気側」を向いていないと、その面をすり抜ける。
CAD で面を 1 枚ずつ描くと巻き順と押し出し方向で向きが決まってしまうので、
**モデル側が正解を持っているとは限らない**。自動判定も、隙間のある形状では外すことがある。

判定の色分け：

| 色 | 意味 |
|---|---|
| 緑 | 法線が室内（空気側）を向いている。そのままでよい |
| 赤 | 法線が室外を向いている。**反転が要る** |
| 灰 | 判定できない（開いた形状など）。CAD の指定を尊重する |

壁が不透明だと**奥の面の向きが見えない**ので、左パネルでレイヤごとに
表示 ON/OFF と不透明度を変えられるようにしてある（既定の不透明度は 0.55）。

操作:
    ドラッグ 回転 / ホイール 拡大縮小 / `z` `x` `c` `v` 視点 / `r` リセット
    `p` 面の選択モード（枠で囲むと、その面を反転）
    `1`〜`9`  そのレイヤをまとめて反転
    `a` 自動判定どおりに揃える   `d` CAD の巻き順に戻す   `i` 全反転
    `n` 法線の矢印 ON/OFF        `o` 不透明度の対象を切り替え
    `g` **いまの画面を画像で保存**（`図/画面/法線_01.png` … 連番）
    `s` **保存して閉じる**       `q` 保存せずに閉じる

`p` は VTK のピック（面の枠選択）に取られているので、撮影は `g`（grab）にしてある。
"""

import numpy as np
import pyvista as pv

import read_dxffile as rd
import view_model_gui as vg

BACKGROUND = "#12151c"
TEXT_COLOR = "#d6dae2"

# 判定の色。cell_data['verdict'] の値 0/1/2 に対応する
VERDICT_COLORS = ["#4cc38a",    # 0 = 内向き（OK）
                  "#e5484d",    # 1 = 外向き（要反転）
                  "#8b929e"]    # 2 = 判定できない
VERDICT_NAMES = ["内向き（OK）", "外向き（要反転）", "判定できない"]


def _visibility(plotter, actor):
    """チェックボックス用のコールバックを作る（クロージャの取り違え防止）。"""
    def callback(flag):
        actor.SetVisibility(flag)
        plotter.render()
    return callback


class NormalEditor:
    """モデルを表示して、面ごとの法線の向きを確認・反転する。

    `flipped` が唯一の状態。**CAD の巻き順から反転する面のインデックス集合**で、
    `read_model(flip_faces=...)` にそのまま渡せる形にしてある
    （自動判定の結果もここに畳み込んでおくので、保存したものを読めば再現できる）。
    """

    def __init__(self, model, flipped=None, title="法線の確認", head_azimuth=None,
                 save_dir=None):
        # 受音点に置く「人」の正面方向 [度]（真上から見て +X から反時計回り）。
        # G-5 の伝搬方向の図で「前・後ろ・左・右」を決めるのに使う。
        # CAD で表すのは難しいという判断で、ここ（3D が見えている画面）で決める
        self.head_azimuth = None if head_azimuth is None else float(head_azimuth)
        self.head_actor = None
        self.model = model
        self.mesh = model.mesh
        self.count = len(self.mesh)
        self.triangles = [tuple(f.vertexes) for f in self.mesh]

        # CAD の巻き順そのままの法線（反転を重ねる前の状態）を復元しておく。
        # model.flipped_faces には読み込み時に反転した面が入っている
        self.cad_normal = np.array([
            (-f.normal if j in model.flipped_faces else f.normal)
            for j, f in enumerate(self.mesh)], dtype=float)

        self.flipped = set(model.flipped_faces if flipped is None else flipped)
        self.saved = False
        self.title = title
        # `g` で撮った画像の置き場（`図/画面/`）。指定が無ければ撮影キーを出さない
        self.save_dir = save_dir

        # 自動判定（レイの偶奇）。ここでは**反転するか否かの判定にしか使わない**
        self.layers = sorted({f.material for f in self.mesh})
        self.layer_of = np.array([self.layers.index(f.material) for f in self.mesh])
        self.enclosure = rd.encloses_point(
            self.triangles,
            model.source_points[0] if model.source_points
            else np.mean([np.mean(t, axis=0) for t in self.triangles], axis=0))
        self.auto_flip, self.ambiguous = rd.orient_inward(self.triangles, self.cad_normal)
        self.reliable = self.enclosure >= rd.ENCLOSURE_THRESHOLD

        self.plotter = None
        self.panel = None
        self.surfaces = []      # [(面インデックス, PolyData)] をレイヤごとに
        self.arrows = None
        self.show_normals = False
        self.label = None

    # ---- 状態 ----------------------------------------------------------

    def normals(self):
        """いまの反転指定を適用した法線 (M,3)。"""
        sign = np.where(np.isin(np.arange(self.count), list(self.flipped)), -1.0, 1.0)
        return self.cad_normal * sign[:, None]

    def verdict(self):
        """面ごとの判定 0=内向き / 1=外向き / 2=判定できない。"""
        if not self.reliable:
            return np.full(self.count, 2, dtype=np.int64)
        # 自動判定は「CAD の巻き順から反転すべき面」。
        # いまの反転指定と一致していれば内向きになっている
        want = np.isin(np.arange(self.count), list(self.auto_flip))
        now = np.isin(np.arange(self.count), list(self.flipped))
        return np.where(want == now, 0, 1).astype(np.int64)

    def toggle(self, indices):
        for j in np.atleast_1d(indices):
            j = int(j)
            if j in self.flipped:
                self.flipped.discard(j)
            else:
                self.flipped.add(j)
        self.refresh()

    def set_auto(self):
        if not self.reliable:
            print("[normal_editor] 開いた形状なので自動判定は使えません（CAD のままにします）")
            self.flipped = set()
        else:
            self.flipped = set(self.auto_flip)
        self.refresh()

    def set_cad(self):
        self.flipped = set()
        self.refresh()

    def flip_all(self):
        self.flipped = set(range(self.count)) - self.flipped
        self.refresh()

    def toggle_layer(self, layer_index):
        if not 0 <= layer_index < len(self.layers):
            return
        self.toggle(np.nonzero(self.layer_of == layer_index)[0])

    # ---- 表示 ----------------------------------------------------------

    def _build_surface(self, faces):
        """指定した面だけの PolyData を作る（レイヤごとに 1 つ作る）。

        レイヤごとに分けるのは、**レイヤ単位で表示 ON/OFF と不透明度を変えたい**ため。
        壁が不透明なままだと中の面の向きが確認できない（ユーザー指摘）。
        """
        faces = np.asarray(faces, dtype=np.int64)
        points = np.concatenate([np.asarray(self.triangles[j]) for j in faces])
        cells = np.hstack([[3, 3 * i, 3 * i + 1, 3 * i + 2]
                           for i in range(len(faces))])
        surface = pv.PolyData(points, faces=cells.astype(np.int64))
        # 選択した面から**元の面番号**を引けるようにしておく（反転の対象を決めるのに要る）
        surface.cell_data["face_id"] = faces
        surface.cell_data["verdict"] = self.verdict()[faces]
        return surface

    def refresh(self, render=True):
        if not self.surfaces:
            return
        verdict = self.verdict()
        for faces, surface in self.surfaces:
            surface.cell_data["verdict"] = verdict[faces]
            surface.Modified()
        if self.arrows is not None:
            self.plotter.remove_actor(self.arrows, render=False)
            self.arrows = None
        if self.show_normals:
            self._add_arrows()
        self._refresh_label()
        if render:
            self.plotter.render()

    def _add_arrows(self):
        centres = np.array([np.mean(t, axis=0) for t in self.triangles])
        length = float(np.linalg.norm(self.model.extents[1] - self.model.extents[0])) * 0.04
        cloud = pv.PolyData(centres)
        cloud["vector"] = self.normals() * length
        cloud.point_data["verdict"] = self.verdict()
        glyph = cloud.glyph(orient="vector", scale=False, factor=1.0,
                            geom=pv.Arrow(tip_length=0.3, shaft_radius=0.02,
                                          tip_radius=0.07))
        self.arrows = self.plotter.add_mesh(glyph, scalars="verdict",
                                            cmap=VERDICT_COLORS, clim=(0, 2),
                                            show_scalar_bar=False, lighting=False)

    def _refresh_label(self):
        if self.label is None:
            return
        counts = np.bincount(self.verdict(), minlength=3)
        lines = [f"面 {self.count} 枚 / 反転中 {len(self.flipped)} 枚"]
        for value, name in enumerate(VERDICT_NAMES):
            if counts[value]:
                lines.append(f"  {name}: {counts[value]} 枚")
        if not self.reliable:
            lines.append(f"※開いた形状（囲まれ度 {self.enclosure:.2f}）なので")
            lines.append(f"  自動判定は使いません")
        elif self.ambiguous:
            lines.append(f"※{self.ambiguous} 枚は自動判定が割れました")
        vg.set_actor_text(self.label, "\n".join(lines))

    def show(self, off_screen=False, screenshot=None, window_size=(1280, 860),
             opacity=0.55, panel=None):
        """ウィンドウを開く。保存されたら True を返す。

        `opacity` の既定を 1.0 でなく 0.55 にしてあるのは、
        **不透明だと手前の壁に隠れて奥の面の向きが確認できない**ため。
        左パネルのスライダとレイヤのチェックボックスでさらに調整できる。
        """
        want_panel = (not off_screen) if panel is None else bool(panel)
        self.plotter, panel = vg.make_plotter(self.title, window_size, off_screen,
                                              panel=want_panel, screen="normals")
        self.panel = panel
        font = vg.japanese_font()

        # レイヤごとに面をまとめる（表示 ON/OFF と不透明度をレイヤ単位で効かせるため）
        self.surfaces = []
        registry = {}
        for k, name in enumerate(self.layers):
            faces = np.nonzero(self.layer_of == k)[0]
            surface = self._build_surface(faces)
            actor = self.plotter.add_mesh(
                surface, scalars="verdict", cmap=VERDICT_COLORS, clim=(0, 2),
                show_scalar_bar=False, show_edges=True, edge_color="#3a4150",
                line_width=1, opacity=opacity,
                # **裏から見ている面はより透ける**（音線ビューアと同じ扱い）。
                # 室の外から覗くと手前の壁は裏側なので、そこが薄くなって中が見える
                backface_params={"opacity": opacity * vg.BACKFACE_OPACITY_RATIO})
            self.surfaces.append((faces, surface))
            registry[name] = {"face": actor, "arrow": None,
                              "colour": "#8b929e", "opacity": opacity}
        vg._attach(self.plotter, "geosim_layers", registry)
        vg._attach(self.plotter, "geosim_panel", panel)

        # 音源・受音点は音線ビューアと同じ球で描く（見た目を揃える）
        lo, hi = self.model.extents
        radius = float(np.linalg.norm(np.asarray(hi) - np.asarray(lo))) * 0.012
        for points, colour in [(self.model.source_points, "#ff5f5f"),
                               (self.model.receiver_points, "#4dd0a0")]:
            for point in points:
                self.plotter.add_mesh(
                    pv.Sphere(radius=radius, center=np.asarray(point)),
                    color=colour, lighting=False)
        self.plotter.add_axes(color=TEXT_COLOR)
        self.plotter.show_bounds(grid="back", location="outer", ticks="outside",
                                 font_size=9, color="#7f8794", xtitle="X [m]",
                                 ytitle="Y [m]", ztitle="Z [m]")

        if panel is not None:
            panel.text(self.title, size=11, color=TEXT_COLOR)
            self.label = panel.reserve_text(5)      # 判定の内訳は最大 5 行

            panel.heading("レイヤ表示（数字キーで反転）")
            for k, name in enumerate(self.layers):
                faces = np.nonzero(self.layer_of == k)[0]
                label = f"{k + 1}: {name} ({len(faces)})" if k < 9 \
                    else f"{name} ({len(faces)})"
                panel.checkbox(label, True,
                               _visibility(self.plotter, registry[name]["face"]),
                               colour="#4cc9f0")

            vg.add_opacity_control(self.plotter, font=font, panel=panel,
                                   target_key="o")

            # 受音点に置く「人」の正面方向（G-5 の伝搬方向の図で使う）。
            # CAD で表すのは難しいので、3D が見えているここで決める
            if self.head_azimuth is not None and self.model.receiver_points:
                panel.heading("受音点の向き（伝搬方向の図で使う）")
                panel.slider("正面の方位 [°]", [0.0, 360.0], self.head_azimuth,
                             lambda v: self.set_head_azimuth(v), fmt="%.0f")
                panel.text("0°=+X / 90°=+Y（真上から見て反時計回り）\n"
                           "黄色い矢印が正面です", size=8)

            panel.heading("操作")
            panel.text("p 面を枠で選んで反転\n"
                       "数字キー レイヤごと反転\n"
                       "a 自動判定に揃える\n"
                       "d CAD の巻き順に戻す\n"
                       "i 全反転   n 法線矢印\n"
                       "z/x/c/v 視点   o 不透明度の対象\n"
                       + ("g いまの画面を画像で保存\n" if self.save_dir else "")
                       + "s 保存して閉じる\n"
                       "q 保存せず閉じる", color="#7f8794")
        else:
            self.label = self.plotter.add_text(" ", position=(14, window_size[1] - 110),
                                               font_size=10, color=TEXT_COLOR,
                                               font_file=font)

        for key, action in (("z", self.plotter.view_xy), ("x", self.plotter.view_xz),
                            ("c", self.plotter.view_yz),
                            ("v", self.plotter.view_isometric)):
            self.plotter.add_key_event(key, action)
        self.plotter.add_key_event("a", self.set_auto)
        self.plotter.add_key_event("d", self.set_cad)      # default（CAD の巻き順）
        self.plotter.add_key_event("i", self.flip_all)     # invert
        self.plotter.add_key_event("n", self._toggle_normals)
        self.plotter.add_key_event("s", self._save_and_close)
        # `p` は面の枠選択（VTK のピック）に取られているので、撮影は `g`（grab）
        if self.save_dir:
            vg.add_screenshot_key(self.plotter, self.save_dir, "法線", key="g")
        for k in range(9):
            self.plotter.add_key_event(str(k + 1),
                                       lambda k=k: self.toggle_layer(k))

        # 枠で囲んだ面を反転する。through=False で**手前に見えている面だけ**を拾う
        self.plotter.enable_cell_picking(callback=self._picked, through=False,
                                         show_message=False, color="#ffd166")

        if self.head_azimuth is not None:
            self.set_head_azimuth(self.head_azimuth, render=False)
        self.refresh(render=False)
        if panel is not None:
            panel.enable_value_input("e")
            panel.relayout()
        self.plotter.view_isometric()
        if off_screen:
            if screenshot:
                self.plotter.screenshot(screenshot)
            self.plotter.close()
            return self.saved
        vg.finish_window(self.plotter)
        self.plotter.show()
        return self.saved

    def _picked(self, picked):
        if picked is None or picked.n_cells == 0:
            return
        ids = picked.cell_data.get("face_id")
        if ids is None:
            print("[normal_editor] 選択から面番号を取れませんでした")
            return
        self.toggle(np.asarray(ids))

    def set_head_azimuth(self, degrees, render=True):
        """受音点に置く「人」の正面方向を変え、矢印を描き直す。

        上下の向きは扱わない（実務では水平面で足りるというユーザー判断）。
        矢印は**受音点から正面へ**伸ばす。長さは室の対角の 8% にしてあり、
        「どちらを向いているか」が分かればよい大きさ。
        """
        self.head_azimuth = float(degrees) % 360.0
        if self.plotter is None or not self.model.receiver_points:
            return
        if self.head_actor is not None:
            self.plotter.remove_actor(self.head_actor, render=False)
            self.head_actor = None

        lo, hi = self.model.extents
        length = float(np.linalg.norm(np.asarray(hi) - np.asarray(lo))) * 0.08
        angle = np.deg2rad(self.head_azimuth)
        direction = np.array([np.cos(angle), np.sin(angle), 0.0])
        start = np.asarray(self.model.receiver_points[0], dtype=float)
        arrow = pv.Arrow(start=start, direction=direction, scale=length,
                         tip_length=0.3, tip_radius=0.12, shaft_radius=0.04)
        self.head_actor = self.plotter.add_mesh(arrow, color="#ffd166",
                                                lighting=False)
        self._refresh_label()
        if render:
            self.plotter.render()

    def _toggle_normals(self):
        self.show_normals = not self.show_normals
        self.refresh()

    def _save_and_close(self):
        self.saved = True
        self.plotter.close()


def edit(project, model=None, off_screen=False, screenshot=None):
    """プロジェクトの法線指定と受音点の向きを確認・修正して保存する。

    法線は `normals.json`、受音点の向きは `project.json` に入る
    （前者はモデルの性質、後者は計算条件なので置き場所を分けている）。

    戻り値 (保存したか, 反転する面の集合)。
    """
    if model is None:
        model = load_model_for(project)
    editor = NormalEditor(model, title=f"{project.name} — 法線の確認",
                          head_azimuth=getattr(project, "head_azimuth", 0.0),
                          save_dir=project.screenshot_dir())
    saved = editor.show(off_screen=off_screen, screenshot=screenshot)
    if saved:
        path = project.save_flipped_faces(editor.flipped, editor.count,
                                          mode=model.orient_mode)
        print(f"[normal_editor] 法線の指定を保存しました: {path}"
              f"（反転 {len(editor.flipped)} / {editor.count} 枚）")
        if editor.head_azimuth is not None:
            project.head_azimuth = editor.head_azimuth
            project.save()
            print(f"[normal_editor] 受音点の向きを保存しました: "
                  f"正面 {editor.head_azimuth:.0f}°")
    else:
        print("[normal_editor] 保存せずに閉じました（前回の指定のままです）")
    return saved, editor.flipped


def load_model_for(project, verbose=True):
    """プロジェクトの設定で DXF を読む（法線の手動指定も適用する）。"""
    absorption = project.absorption_path
    table = None
    if absorption:
        import absorption as ab
        library = ab.MaterialLibrary.from_csv(absorption, kind=project.absorption_kind)
        table = library.absorption_table(project.assignment,
                                         band_number=project.band_number)
    # 面数が分からないと反転指定を照合できないので、まず一度読む
    first = rd.read_model(project.dxf_path, unit=project.unit, absorption_table=table,
                          orient_normals=project.orient_normals,
                          band_number=project.band_number, verbose=False)
    flipped = project.flipped_faces_for(len(first.mesh))
    if not flipped:
        if verbose:
            print("[read_dxffile] " + first.summary().replace("\n", "\n[read_dxffile] "))
        return first
    return rd.read_model(project.dxf_path, unit=project.unit, absorption_table=table,
                         orient_normals=project.orient_normals,
                         band_number=project.band_number,
                         flip_faces=flipped, verbose=verbose)


def main():
    import argparse
    import project as pj

    p = argparse.ArgumentParser(description="法線の向きの確認・修正（G-8）")
    p.add_argument("folder", help="プロジェクトフォルダ")
    p.add_argument("--screenshot", help="画像に書き出して終了（確認用）")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf_path:
        raise SystemExit(f"{a.folder} に project.json が無いか、DXF が設定されていません")
    edit(project, off_screen=a.screenshot is not None, screenshot=a.screenshot)


if __name__ == "__main__":
    main()
