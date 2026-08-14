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

操作:
    ドラッグ 回転 / ホイール 拡大縮小 / `r` 視点リセット
    `p` 面の選択モード（枠で囲むと、その面を反転）
    `1`〜`9`  そのレイヤをまとめて反転
    `a` 自動判定どおりに揃える   `c` CAD の巻き順に戻す   `x` 全反転
    `n` 法線の矢印 ON/OFF        `w` ワイヤフレーム ON/OFF
    `s` **保存して閉じる**       `q` 保存せずに閉じる
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


class NormalEditor:
    """モデルを表示して、面ごとの法線の向きを確認・反転する。

    `flipped` が唯一の状態。**CAD の巻き順から反転する面のインデックス集合**で、
    `read_model(flip_faces=...)` にそのまま渡せる形にしてある
    （自動判定の結果もここに畳み込んでおくので、保存したものを読めば再現できる）。
    """

    def __init__(self, model, flipped=None, title="法線の確認"):
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
        self.surface = None
        self.arrows = None
        self.show_normals = False

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

    def _build_surface(self):
        points = np.concatenate([np.asarray(t) for t in self.triangles])
        cells = np.hstack([[3, 3 * i, 3 * i + 1, 3 * i + 2] for i in range(self.count)])
        surface = pv.PolyData(points, faces=cells.astype(np.int64))
        # 選択した面から元の面番号を引けるようにしておく
        surface.cell_data["face_id"] = np.arange(self.count, dtype=np.int64)
        surface.cell_data["verdict"] = self.verdict()
        return surface

    def refresh(self, render=True):
        if self.surface is None:
            return
        self.surface.cell_data["verdict"] = self.verdict()
        self.surface.Modified()
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
        counts = np.bincount(self.verdict(), minlength=3)
        lines = [f"面 {self.count} 枚 / 反転中 {len(self.flipped)} 枚"]
        for value, name in enumerate(VERDICT_NAMES):
            if counts[value]:
                lines.append(f"  {name}: {counts[value]} 枚")
        if not self.reliable:
            lines.append(f"  ※開いた形状（囲まれ度 {self.enclosure:.2f}）なので")
            lines.append(f"    自動判定は使いません。CAD の指定を尊重します")
        elif self.ambiguous:
            lines.append(f"  ※{self.ambiguous} 枚は自動判定が割れました（目で確認してください）")
        vg.set_actor_text(self.label, "\n".join(lines))

    def _layer_label(self):
        lines = ["レイヤ（数字キーでまとめて反転）"]
        for k, name in enumerate(self.layers[:9]):
            n = int(np.count_nonzero(self.layer_of == k))
            flipped = int(np.count_nonzero(
                [j in self.flipped for j in np.nonzero(self.layer_of == k)[0]]))
            lines.append(f"  {k + 1}: {name}  {n} 枚"
                         + (f"（{flipped} 枚反転中）" if flipped else ""))
        if len(self.layers) > 9:
            lines.append(f"  … 他 {len(self.layers) - 9} レイヤ（数字キーは 9 まで）")
        return "\n".join(lines)

    def show(self, off_screen=False, screenshot=None, window_size=(1280, 860)):
        """ウィンドウを開く。保存されたら True を返す。"""
        font = vg.japanese_font()
        self.plotter = pv.Plotter(title=self.title, window_size=window_size,
                                  off_screen=off_screen)
        self.plotter.set_background(BACKGROUND)

        self.surface = self._build_surface()
        self.plotter.add_mesh(self.surface, scalars="verdict",
                              cmap=VERDICT_COLORS, clim=(0, 2),
                              show_scalar_bar=False, show_edges=True,
                              edge_color="#3a4150", line_width=1, opacity=1.0)

        for point, colour, name in [(self.model.source_points, "#ffd166", "音源"),
                                    (self.model.receiver_points, "#4cc9f0", "受音点")]:
            if point:
                self.plotter.add_mesh(pv.PolyData(np.array(point)), color=colour,
                                      point_size=16, render_points_as_spheres=True)

        self.label = self.plotter.add_text(" ", position=(14, window_size[1] - 110),
                                           font_size=10, color=TEXT_COLOR, font_file=font)
        self.plotter.add_text(self._layer_label(), position=(14, 150), font_size=9,
                              color="#9aa2b1", font_file=font)
        self.plotter.add_text(
            "p 面を選んで反転   a 自動   c CAD に戻す   x 全反転\n"
            "n 法線矢印   s 保存して閉じる   q 保存せず閉じる",
            position=(14, 14), font_size=10, color=TEXT_COLOR, font_file=font)

        self.plotter.add_key_event("a", self.set_auto)
        self.plotter.add_key_event("c", self.set_cad)
        self.plotter.add_key_event("x", self.flip_all)
        self.plotter.add_key_event("n", self._toggle_normals)
        self.plotter.add_key_event("s", self._save_and_close)
        for k in range(9):
            self.plotter.add_key_event(str(k + 1),
                                       lambda k=k: self.toggle_layer(k))

        # 枠で囲んだ面を反転する。through=False で**手前に見えている面だけ**を拾う
        self.plotter.enable_cell_picking(callback=self._picked, through=False,
                                         show_message=False, color="#ffd166")

        self.refresh(render=False)
        self.plotter.view_isometric()
        if off_screen:
            if screenshot:
                self.plotter.screenshot(screenshot)
            self.plotter.close()
            return self.saved
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

    def _toggle_normals(self):
        self.show_normals = not self.show_normals
        self.refresh()

    def _save_and_close(self):
        self.saved = True
        self.plotter.close()


def edit(project, model=None, off_screen=False, screenshot=None):
    """プロジェクトの法線指定を確認・修正して `normals.json` に保存する。

    戻り値 (保存したか, 反転する面の集合)。
    """
    if model is None:
        model = load_model_for(project)
    editor = NormalEditor(model, title=f"{project.name} — 法線の確認")
    saved = editor.show(off_screen=off_screen, screenshot=screenshot)
    if saved:
        path = project.save_flipped_faces(editor.flipped, editor.count,
                                          mode=model.orient_mode)
        print(f"[normal_editor] 法線の指定を保存しました: {path}"
              f"（反転 {len(editor.flipped)} / {editor.count} 枚）")
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
