"""面を目で確かめて直すウィンドウ ― **法線の向き**と**吸音材の割り当て**（TODO G-8 / G-9）。

読み込み時の自動補正（`orient_normals='auto'`）や、レイヤ→吸音材の対応に任せきりにせず、
**人が見て確認し、必要なら直せる**ようにするためのもの。結果は `project.json` と同じ
フォルダに保存され、次からは自動判定のあとに重ねて適用される。

    normals.json    法線を反転する面      → `read_model(flip_faces=...)`
    materials.json  面ごとの吸音材        → `read_model(face_materials=...)`

**なぜ法線が要るか**：法線は「音が通る空気側」を向いていないと、その面をすり抜ける。
CAD で面を 1 枚ずつ描くと巻き順と押し出し方向で向きが決まってしまうので、
**モデル側が正解を持っているとは限らない**。自動判定も、隙間のある形状では外すことがある。

**なぜ面ごとの吸音材が要るか**：吸音材は本来レイヤで分けるが、**1 つの 3DSOLID で
出来ていて面ごとのレイヤが無いモデル**（他ソフトの都合でその形式が要るとき）では
レイヤが使えない。そこで、この画面で面を選んで直接貼れるようにしてある。

---

## 面グループ ― この画面の要になる考え方

3DSOLID を STL 経由で取り込むと、設計者が描いた 1 枚の壁が三角形に割られてしまう。
そこで**同一平面で連結した三角形をひとまとめ**にして扱う（`read_dxffile.coplanar_groups`）。
実例：ModelTest は三角形 68 枚 → **16 グループ**（床・天井・壁 14 枚）で、
元のソリッドが持っていた面と一致した。**1 クリックで「壁 1 枚」が選べる。**
`y` で三角形単位にも切り替えられる（グループ分けが合わない形状のため）。

## 選ぶ → 適用する（二段階）

以前は枠で囲んだ瞬間に反転していたため、**何を選んだのか分からなかった**（ユーザー指摘）。
いまは選択が残り、明るい橙で塗られ、輪郭が描かれる。左パネルに枚数と面積が出る。
選び終えてから、法線の反転（`i`）や吸音材（左パネルの材料をクリック）を適用する。

判定の色分け（法線モード）:

| 色 | 意味 |
|---|---|
| 緑 | 法線が室内（空気側）を向いている。そのままでよい |
| 赤 | 法線が室外を向いている。**反転が要る** |
| 灰 | 判定できない（開いた形状など）。CAD の指定を尊重する |
| 橙 | いま選択している面（判定の色より優先して塗る） |

壁が不透明だと**奥の面の向きが見えない**ので、左パネルでレイヤごとに
表示 ON/OFF と不透明度を変えられるようにしてある（既定の不透明度は 0.55）。

操作:
    ドラッグ 回転 / ホイール 拡大縮小 / `z` `x` `c` `v` 視点
    `r` **枠選択のオン・オフ**（囲んだ面を選択に**追加**する。VTK のキー）
    `0` 選択を解除   `j` 全選択   `h` 選択を反転
    `k` 同じ向きの面をまとめて選ぶ（床・天井・壁の一括に使う）
    `l` 同じ吸音材の面をまとめて選ぶ
    `1`〜`9` そのレイヤの面を選択に足す
    `y` 面グループ ⇔ 三角形 の切り替え
    `m` 表示を 法線 ⇔ 吸音材 で切り替え
    `i` 選択した面の法線を反転（**選択が空なら全部**）
    `a` 自動判定どおりに揃える   `d` CAD の巻き順に戻す
    `n` 法線の矢印 ON/OFF        `o` 不透明度の対象を切り替え
    `t` 値を数字で入力（不透明度・正面の方位）
    `g` **いまの画面を画像で保存**（`図/画面/面_01.png` … 連番）
    `s` **保存して閉じる**       `q` 保存せずに閉じる

吸音材は**左パネルの材料をクリックすると、選択中の面に貼られる**。
「未設定（レイヤに戻す）」を選べば剥がせる。

**`e` は使えない**（VTK の終了キー。`view_model_gui.VTK_RESERVED_KEYS` 参照）。
"""

import numpy as np
import pyvista as pv

import read_dxffile as rd
import view_model_gui as vg
from view_model import LAYER_PALETTE

BACKGROUND = "#12151c"
TEXT_COLOR = "#d6dae2"

# 判定の色。法線モードのときの面の色
VERDICT_COLORS = ["#4cc38a",    # 0 = 内向き（OK）
                  "#e5484d",    # 1 = 外向き（要反転）
                  "#8b929e"]    # 2 = 判定できない
VERDICT_NAMES = ["内向き（OK）", "外向き（要反転）", "判定できない"]

# 選択中の面の色。判定・材料のどちらの色より優先して塗る
SELECTED_COLOR = "#ffd166"
# 吸音材が割り当てられていない面の色（吸音材モード）
UNASSIGNED_COLOR = "#4a5160"

# 表示モード
MODE_NORMALS = "normals"
MODE_MATERIALS = "materials"
MODE_NAMES = {MODE_NORMALS: "法線の向き", MODE_MATERIALS: "吸音材"}

# 「同じ向きの面」とみなす角度。床・天井・壁をまとめて選ぶのに使う。
# 面グループのしきい値（1°）より緩くしてあるのは、**別々の平面でも向きが揃っていれば
# まとめたい**ため（例：段差のある天井を一度に選ぶ）
SAME_NORMAL_DEGREES = 5.0

# 「未設定（レイヤに戻す）」を表す番兵。材料名として使えない文字を含めてある
UNASSIGNED = "\x00未設定"

# 左上の状態表示に確保する行数と、1 行に入れる目安の文字数。
# `Panel.reserve_text` は `Panel.text` と違って**切り詰めをしない**ので、
# 長い行を渡すとパネルからはみ出して右が読めなくなる。こちら側で短く保つ
LABEL_LINES = 7
LABEL_WIDTH = 26


def _wrap(parts, width=LABEL_WIDTH, separator=" / "):
    """短い断片を、1 行が `width` 文字を超えないように詰めて並べる。"""
    lines, current = [], ""
    for part in parts:
        candidate = part if not current else current + separator + part
        if len(candidate) > width and current:
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def _hex_to_rgb(code):
    code = code.lstrip("#")
    return np.array([int(code[i:i + 2], 16) for i in (0, 2, 4)], dtype=np.uint8)


def _visibility(plotter, actor):
    """チェックボックス用のコールバックを作る（クロージャの取り違え防止）。"""
    def callback(flag):
        actor.SetVisibility(flag)
        plotter.render()
    return callback


class FaceEditor:
    """モデルを表示して、面ごとの法線の向きと吸音材を確認・修正する。

    状態は 2 つだけ。

    `flipped`
        **CAD の巻き順から反転する面のインデックス集合**で、`read_model(flip_faces=...)`
        にそのまま渡せる形にしてある（自動判定の結果もここに畳み込んでおくので、
        保存したものを読めば再現できる）。
    `assigned`
        **{面インデックス: 材料名}**。`read_model(face_materials=...)` に渡す。
        入っていない面はレイヤから引く（従来どおり）。

    `selection` は作業用で、保存しない。
    """

    def __init__(self, model, flipped=None, face_materials=None, materials=None,
                 title="面の確認", head_azimuth=None, save_dir=None):
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
        self.assigned = {int(k): v for k, v in (face_materials or {}).items() if v}
        self.saved = False
        self.title = title
        # `g` で撮った画像の置き場（`図/画面/`）。指定が無ければ撮影キーを出さない
        self.save_dir = save_dir

        # 面グループ（同一平面パッチ）。**法線の向きに依らない幾何の性質**なので、
        # 反転しても組み替わらない（read_dxffile.coplanar_groups の説明を参照）
        self.group_of = rd.coplanar_groups(self.triangles, self.cad_normal)
        self.groups = [np.nonzero(self.group_of == g)[0]
                       for g in range(int(self.group_of.max()) + 1 if self.count else 0)]
        self.by_group = True        # `y` で三角形単位に切り替えられる

        self.area = np.array([_triangle_area(t) for t in self.triangles])
        self.selection = set()
        self.mode = MODE_NORMALS

        # 材料一覧（吸音率 CSV から読んだもの）。無ければ吸音材モードは使えない
        self.materials = list(materials or [])
        self.material_colour = {name: LAYER_PALETTE[i % len(LAYER_PALETTE)]
                                for i, name in enumerate(self.materials)}

        # レイヤは**DXF の実際のレイヤ**で持つ。`Mesh.material` は面ごとの割り当てで
        # 材料名に化けるので、そちらを使うとレイヤ表示が材料名になってしまう
        face_layers = (model.face_layers if len(model.face_layers) == self.count
                       else [f.material for f in self.mesh])
        self.layers = sorted(set(face_layers))
        self.layer_of = np.array([self.layers.index(name) for name in face_layers])

        # 自動判定（レイの偶奇）。ここでは**反転するか否かの判定にしか使わない**
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
        self.outline = None
        self.show_normals = False
        self.label = None
        self.status = None

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

    def face_colours(self):
        """面ごとの表示色 (M,3) uint8。表示モードと選択で決まる。"""
        colours = np.zeros((self.count, 3), dtype=np.uint8)
        if self.mode == MODE_NORMALS:
            palette = np.array([_hex_to_rgb(c) for c in VERDICT_COLORS])
            colours[:] = palette[self.verdict()]
        else:
            colours[:] = _hex_to_rgb(UNASSIGNED_COLOR)
            for j, name in self.assigned.items():
                colours[j] = _hex_to_rgb(self.material_colour.get(name, "#ffffff"))
        # 選択は判定・材料より優先して塗る（いま何を選んでいるかが最優先の情報）
        if self.selection:
            colours[sorted(self.selection)] = _hex_to_rgb(SELECTED_COLOR)
        return colours

    # ---- 選択 ----------------------------------------------------------

    def expand(self, indices):
        """面グループ単位なら、同じグループの面へ広げる。三角形単位ならそのまま。"""
        indices = np.atleast_1d(np.asarray(indices, dtype=np.int64))
        if not self.by_group or not len(indices):
            return indices
        wanted = set(int(g) for g in self.group_of[indices])
        return np.nonzero(np.isin(self.group_of, list(wanted)))[0]

    def select(self, indices, render=True):
        """選択に**足す**（枠選択を繰り返して集められるように）。"""
        self.selection.update(int(j) for j in self.expand(indices))
        self.refresh(render=render)

    def clear_selection(self):
        self.selection.clear()
        self.refresh()

    def select_all(self):
        self.selection = set(range(self.count))
        self.refresh()

    def invert_selection(self):
        self.selection = set(range(self.count)) - self.selection
        self.refresh()

    def select_same_normal(self):
        """選択中の面と**向きが揃う**面をまとめて選ぶ（床・天井・壁の一括）。"""
        if not self.selection:
            self._say("先に面を選んでください（r で枠選択）")
            return
        unit = self.normals()
        picked = unit[sorted(self.selection)]
        cosine = np.cos(np.deg2rad(SAME_NORMAL_DEGREES))
        # どれか 1 つでも向きが揃えば仲間とみなす
        match = (unit @ picked.T >= cosine).any(axis=1)
        self.selection = set(np.nonzero(match)[0].tolist())
        self.refresh()

    def select_same_material(self):
        """選択中の面と**同じ吸音材**が貼られている面をまとめて選ぶ。"""
        if not self.selection:
            self._say("先に面を選んでください（r で枠選択）")
            return
        names = {self.assigned.get(j) for j in self.selection}
        chosen = set()
        for j in range(self.count):
            if self.assigned.get(j) in names:
                chosen.add(j)
        self.selection = chosen
        self.refresh()

    def select_layer(self, layer_index):
        if not 0 <= layer_index < len(self.layers):
            return
        self.select(np.nonzero(self.layer_of == layer_index)[0])

    def toggle_unit(self):
        self.by_group = not self.by_group
        # 単位を変えたら選択も新しい単位に合わせて広げ直す（見た目と実態を合わせる）
        if self.by_group and self.selection:
            self.selection = set(int(j) for j in self.expand(sorted(self.selection)))
        self.refresh()

    # ---- 適用 ----------------------------------------------------------

    def flip_selection(self):
        """選択した面の法線を反転する。**選択が空なら全部**（従来の `i` と同じ）。"""
        target = sorted(self.selection) if self.selection else range(self.count)
        for j in target:
            if j in self.flipped:
                self.flipped.discard(j)
            else:
                self.flipped.add(j)
        self.refresh()

    def assign(self, name):
        """選択した面に吸音材を貼る。`UNASSIGNED` なら剥がしてレイヤに戻す。"""
        if not self.selection:
            self._say("先に面を選んでください（r で枠選択）")
            return
        for j in self.selection:
            if name is UNASSIGNED:
                self.assigned.pop(j, None)
            else:
                self.assigned[j] = name
        # 貼った結果が見えるように吸音材モードへ切り替える
        self.mode = MODE_MATERIALS
        self.refresh()

    def set_auto(self):
        if not self.reliable:
            self._say("開いた形状なので自動判定は使えません（CAD のままにします）")
            self.flipped = set()
        else:
            self.flipped = set(self.auto_flip)
        self.refresh()

    def set_cad(self):
        self.flipped = set()
        self.refresh()

    def toggle_mode(self):
        self.mode = MODE_MATERIALS if self.mode == MODE_NORMALS else MODE_NORMALS
        self.refresh()

    # ---- 表示 ----------------------------------------------------------

    def _build_surface(self, faces):
        """指定した面だけの PolyData を作る（レイヤごとに 1 つ作る）。

        レイヤごとに分けるのは、**レイヤ単位で表示 ON/OFF と不透明度を変えたい**ため。
        壁が不透明なままだと中の面の向きが確認できない（ユーザー指摘）。

        色は cmap ではなく **RGB を直接持たせる**。判定（3 色）と吸音材（材料ごと）で
        色の意味が変わるので、モードを切り替えるたびに cmap を差し替えるより、
        配列を書き換えるだけで済むこちらのほうが単純で確実。
        """
        faces = np.asarray(faces, dtype=np.int64)
        points = np.concatenate([np.asarray(self.triangles[j]) for j in faces])
        cells = np.hstack([[3, 3 * i, 3 * i + 1, 3 * i + 2]
                           for i in range(len(faces))])
        surface = pv.PolyData(points, faces=cells.astype(np.int64))
        # 選択した面から**元の面番号**を引けるようにしておく（適用の対象を決めるのに要る）
        surface.cell_data["face_id"] = faces
        surface.cell_data["rgb"] = self.face_colours()[faces]
        return surface

    def _selection_outline(self):
        """選択している面のかたまりの**外周**を線で返す（無ければ None）。

        中の三角形の辺まで引くと網目になって形が読めないので、
        **1 回しか現れない辺＝外周**だけを残す。面グループで選んだときに
        「壁 1 枚がまるごと選ばれている」ことが一目で分かる。
        """
        if not self.selection:
            return None
        seen = {}
        for j in self.selection:
            v = np.asarray(self.triangles[j])
            for a, b in ((0, 1), (1, 2), (2, 0)):
                ka = tuple(np.round(v[a], 9))
                kb = tuple(np.round(v[b], 9))
                edge = (ka, kb) if ka <= kb else (kb, ka)
                seen[edge] = seen.get(edge, 0) + 1
        border = [edge for edge, n in seen.items() if n == 1]
        if not border:
            return None
        points = np.array([p for edge in border for p in edge], dtype=float)
        lines = np.hstack([[2, 2 * i, 2 * i + 1] for i in range(len(border))])
        return pv.PolyData(points, lines=lines.astype(np.int64))

    def refresh(self, render=True):
        if not self.surfaces:
            return
        colours = self.face_colours()
        for faces, surface in self.surfaces:
            surface.cell_data["rgb"] = colours[faces]
            surface.Modified()
        if self.arrows is not None:
            self.plotter.remove_actor(self.arrows, render=False)
            self.arrows = None
        if self.show_normals:
            self._add_arrows()
        if self.outline is not None:
            self.plotter.remove_actor(self.outline, render=False)
            self.outline = None
        border = self._selection_outline()
        if border is not None:
            self.outline = self.plotter.add_mesh(border, color=SELECTED_COLOR,
                                                 line_width=4, lighting=False,
                                                 render_lines_as_tubes=True,
                                                 pickable=False)
        self._refresh_label()
        if render:
            self.plotter.render()

    def _add_arrows(self):
        centres = np.array([np.mean(t, axis=0) for t in self.triangles])
        length = float(np.linalg.norm(self.model.extents[1] - self.model.extents[0])) * 0.04
        cloud = pv.PolyData(centres)
        cloud["vector"] = self.normals() * length
        cloud.point_data["rgb"] = self.face_colours()
        glyph = cloud.glyph(orient="vector", scale=False, factor=1.0,
                            geom=pv.Arrow(tip_length=0.3, shaft_radius=0.02,
                                          tip_radius=0.07))
        self.arrows = self.plotter.add_mesh(glyph, scalars="rgb", rgb=True,
                                            show_scalar_bar=False, lighting=False,
                                            pickable=False)

    def _say(self, message):
        """パネル下部に一言出す（キーを押しても何も起きない理由を伝えるため）。"""
        print(f"[face_editor] {message}")
        if self.status is not None:
            vg.set_actor_text(self.status, message)
            if self.plotter is not None:
                self.plotter.render()

    def _refresh_label(self):
        """左上の状態表示を書き直す。

        `reserve_text` は**切り詰めをしない**ので、パネルからはみ出さないよう
        1 行を短く保つこと（長いと右が切れて読めなくなる）。
        """
        if self.label is None:
            return
        unit = "面グループ" if self.by_group else "三角形"
        lines = [f"表示:{MODE_NAMES[self.mode]} / 単位:{unit}"]
        if self.selection:
            groups = len({int(g) for g in self.group_of[sorted(self.selection)]})
            area = float(self.area[sorted(self.selection)].sum())
            lines.append(f"選択 {groups}グループ {len(self.selection)}面 {area:.1f}m2")
        else:
            lines.append(f"選択なし（{self.count}面 {len(self.groups)}グループ）")
        if self.mode == MODE_NORMALS:
            counts = np.bincount(self.verdict(), minlength=3)
            lines.append(f"反転中 {len(self.flipped)}面")
            lines += _wrap(f"{name} {counts[v]}"
                           for v, name in enumerate(VERDICT_NAMES) if counts[v])
            if not self.reliable:
                lines.append(f"※開いた形状（囲まれ度 {self.enclosure:.2f}）")
                lines.append("　自動判定は使いません")
            elif self.ambiguous:
                lines.append(f"※{self.ambiguous} 面は自動判定が割れました")
        else:
            done = len(self.assigned)
            lines.append(f"吸音材あり {done}面 / 未設定 {self.count - done}面")
            counts = {}
            for name in self.assigned.values():
                counts[name] = counts.get(name, 0) + 1
            lines += _wrap(f"{k} {v}" for k, v in sorted(counts.items()))
        vg.set_actor_text(self.label, "\n".join(lines[:LABEL_LINES]))

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
                surface, scalars="rgb", rgb=True,
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
                    color=colour, lighting=False, pickable=False)
        self.plotter.add_axes(color=TEXT_COLOR)
        self.plotter.show_bounds(grid="back", location="outer", ticks="outside",
                                 font_size=9, color="#7f8794", xtitle="X [m]",
                                 ytitle="Y [m]", ztitle="Z [m]")

        if panel is not None:
            panel.text(self.title, size=11, color=TEXT_COLOR)
            self.label = panel.reserve_text(LABEL_LINES)

            self._build_material_panel(panel)

            panel.heading("レイヤ表示（数字キーで選択）")
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
                panel.heading("受音点の向き（伝搬方向の図）")
                panel.slider("正面の方位 [°]", [0.0, 360.0], self.head_azimuth,
                             lambda v: self.set_head_azimuth(v), fmt="%.0f")
                panel.text("0°=+X / 90°=+Y（真上から見て反時計回り）\n"
                           "黄色い矢印が正面です", size=8)

            panel.heading("操作")
            # ★パネルは縦に伸びず、横も狭い（`Panel.text` は右を「…」で切り詰める）。
            #   材料が増えるほど下が押されるので、**小さめの字で短く**書く
            panel.text("r 枠選択  0 解除  j 全選択\n"
                       "h 選択反転  k 同じ向き\n"
                       "l 同じ吸音材  数字 レイヤ\n"
                       "y グループ⇔三角形  m 法線⇔吸音材\n"
                       "i 法線を反転（空なら全部）\n"
                       "a 自動判定  d CAD の巻き順\n"
                       "n 法線矢印  o 不透明度\n"
                       f"z/x/c/v 視点  {vg.VALUE_INPUT_KEY} 数値入力"
                       + ("  g 画像" if self.save_dir else "") + "\n"
                       "s 保存して閉じる  q 保存せず閉じる",
                       size=8, color="#7f8794")
            self.status = panel.reserve_text(2, size=9, color="#ffd166")
        else:
            self.label = self.plotter.add_text(" ", position=(14, window_size[1] - 130),
                                               font_size=10, color=TEXT_COLOR,
                                               font_file=font)

        for key, action in (("z", self.plotter.view_xy), ("x", self.plotter.view_xz),
                            ("c", self.plotter.view_yz),
                            ("v", self.plotter.view_isometric)):
            self.plotter.add_key_event(key, action)
        self.plotter.add_key_event("a", self.set_auto)
        self.plotter.add_key_event("d", self.set_cad)      # default（CAD の巻き順）
        self.plotter.add_key_event("i", self.flip_selection)   # invert
        self.plotter.add_key_event("n", self._toggle_normals)
        self.plotter.add_key_event("m", self.toggle_mode)      # mode
        self.plotter.add_key_event("y", self.toggle_unit)
        self.plotter.add_key_event("0", self.clear_selection)
        self.plotter.add_key_event("j", self.select_all)
        self.plotter.add_key_event("h", self.invert_selection)
        self.plotter.add_key_event("k", self.select_same_normal)
        self.plotter.add_key_event("l", self.select_same_material)
        self.plotter.add_key_event("s", self._save_and_close)
        # `r` は枠選択（VTK のラバーバンド）に取られているので、撮影は `g`（grab）
        if self.save_dir:
            vg.add_screenshot_key(self.plotter, self.save_dir, "面", key="g")
        for k in range(9):
            self.plotter.add_key_event(str(k + 1),
                                       lambda k=k: self.select_layer(k))

        # 枠で囲んだ面を**選択に足す**。through=False で手前に見えている面だけを拾う
        self.plotter.enable_cell_picking(callback=self._picked, through=False,
                                         show_message=False, color=SELECTED_COLOR)

        if self.head_azimuth is not None:
            self.set_head_azimuth(self.head_azimuth, render=False)
        self.refresh(render=False)
        if panel is not None:
            panel.enable_value_input()
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

    def _build_material_panel(self, panel):
        """吸音材の一覧。**クリックすると選択中の面に貼られる。**

        パネルが持っているのはチェックボックスとスライダだけなので、
        チェックボックスを**ボタンとして**使う（押したら貼って、すぐ戻す）。
        トグルとして残すと「チェックが付いている材料」と「実際に貼られている材料」が
        食い違いを起こすため、状態は持たせない。
        """
        panel.heading("吸音材（クリックで貼る）")
        if not self.materials:
            panel.text("吸音率 CSV が設定されていないので\n"
                       "材料一覧が空です（条件入力で指定してください）",
                       size=8, color="#e5a03d")
            return
        for name in self.materials + [UNASSIGNED]:
            label = "未設定（レイヤに戻す）" if name is UNASSIGNED else name
            colour = "#8b929e" if name is UNASSIGNED \
                else self.material_colour.get(name, "#ffffff")
            panel.checkbox(label, False, self._paint(name), colour=colour)

    def _paint(self, name):
        """材料ボタンのコールバックを作る（クロージャの取り違え防止）。"""
        state = {"busy": False}

        def callback(flag, name=name):
            # チェックを外す操作で自分自身が呼ばれるので、そのぶんは無視する
            if state["busy"] or not flag:
                return
            state["busy"] = True
            try:
                self.assign(name)
            finally:
                state["busy"] = False
        return callback

    def _picked(self, picked):
        """枠選択の結果を受け取る。

        レイヤごとに別の PolyData を出しているので、**複数レイヤにまたがって囲むと
        MultiBlock が返る**。1 つずつ取り出して面番号を集める。
        """
        if picked is None:
            return
        blocks = list(picked) if isinstance(picked, pv.MultiBlock) else [picked]
        ids = []
        for block in blocks:
            if block is None or block.n_cells == 0:
                continue
            found = block.cell_data.get("face_id")
            if found is not None:
                ids.append(np.asarray(found).ravel())
        if not ids:
            self._say("選択から面番号を取れませんでした")
            return
        self.select(np.concatenate(ids))

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
                                                lighting=False, pickable=False)
        self._refresh_label()
        if render:
            self.plotter.render()

    def _toggle_normals(self):
        self.show_normals = not self.show_normals
        self.refresh()

    def _save_and_close(self):
        self.saved = True
        self.plotter.close()


def _triangle_area(triangle):
    v = np.asarray(triangle, dtype=float)
    return 0.5 * float(np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0])))


def material_names(project):
    """プロジェクトの吸音率 CSV から材料名の一覧を返す（無ければ空）。"""
    if not project.absorption_path:
        return []
    import absorption as ab
    try:
        library = ab.MaterialLibrary.from_csv(project.absorption_path,
                                              kind=project.absorption_kind)
    except Exception as error:      # CSV が壊れていても画面は開けるようにする
        print(f"[face_editor] 吸音率 CSV を読めませんでした（{error}）。"
              f"材料一覧は空になります")
        return []
    return library.names()


def edit(project, model=None, off_screen=False, screenshot=None):
    """プロジェクトの法線指定・吸音材の割り当て・受音点の向きを確認して保存する。

    法線は `normals.json`、面ごとの吸音材は `materials.json`、受音点の向きは
    `project.json` に入る（前 2 つはモデルの性質、最後は計算条件なので分けている）。

    戻り値 (保存したか, 反転する面の集合)。
    """
    if model is None:
        model = load_model_for(project)
    editor = FaceEditor(model,
                        face_materials=project.face_materials_for(len(model.mesh)),
                        materials=material_names(project),
                        title=f"{project.name} — 面の確認（法線・吸音材）",
                        head_azimuth=getattr(project, "head_azimuth", 0.0),
                        save_dir=project.screenshot_dir())
    saved = editor.show(off_screen=off_screen, screenshot=screenshot)
    if saved:
        path = project.save_flipped_faces(editor.flipped, editor.count,
                                          mode=model.orient_mode)
        print(f"[face_editor] 法線の指定を保存しました: {path}"
              f"（反転 {len(editor.flipped)} / {editor.count} 枚）")
        path = project.save_face_materials(editor.assigned, editor.count)
        print(f"[face_editor] 吸音材の割り当てを保存しました: {path}"
              f"（{len(editor.assigned)} / {editor.count} 枚）")
        if editor.head_azimuth is not None:
            project.head_azimuth = editor.head_azimuth
            project.save()
            print(f"[face_editor] 受音点の向きを保存しました: "
                  f"正面 {editor.head_azimuth:.0f}°")
    else:
        print("[face_editor] 保存せずに閉じました（前回の指定のままです）")
    return saved, editor.flipped


def load_model_for(project, verbose=True):
    """プロジェクトの設定で DXF を読む（法線・吸音材の手動指定も適用する）。"""
    absorption = project.absorption_path
    table = None
    if absorption:
        import absorption as ab
        library = ab.MaterialLibrary.from_csv(absorption, kind=project.absorption_kind)
        table = library.absorption_table(project.assignment,
                                         band_number=project.band_number)
    # 面数が分からないと反転指定・材料の割り当てを照合できないので、まず一度読む
    first = rd.read_model(project.dxf_path, unit=project.unit, absorption_table=table,
                          orient_normals=project.orient_normals,
                          band_number=project.band_number, verbose=False)
    flipped = project.flipped_faces_for(len(first.mesh))
    materials = project.face_materials_for(len(first.mesh))
    if not flipped and not materials:
        if verbose:
            print("[read_dxffile] " + first.summary().replace("\n", "\n[read_dxffile] "))
        return first
    return rd.read_model(project.dxf_path, unit=project.unit, absorption_table=table,
                         orient_normals=project.orient_normals,
                         band_number=project.band_number,
                         flip_faces=flipped or None, face_materials=materials,
                         verbose=verbose)


def main():
    import argparse
    import project as pj

    p = argparse.ArgumentParser(description="面の確認・修正（法線と吸音材。G-8 / G-9）")
    p.add_argument("folder", help="プロジェクトフォルダ")
    p.add_argument("--screenshot", help="画像に書き出して終了（確認用）")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf_path:
        raise SystemExit(f"{a.folder} に project.json が無いか、DXF が設定されていません")
    edit(project, off_screen=a.screenshot is not None, screenshot=a.screenshot)


if __name__ == "__main__":
    main()
