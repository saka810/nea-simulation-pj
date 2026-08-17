"""音源から音線がどう飛ぶかだけを見る画面。

**室形状とは関係のない話**なので、モデルを読まずに音線ベクトルだけを描く。
条件入力の「音線の飛び方を見る」から開く。

見たいのは次の 2 点。

  ・**球状に、全方向へ均等に出ているか**
  ・本数を変えると密度がどう変わるか（音線数の決め方の目安になる）

音線は `sound_ray.soundray_generator()` が Fibonacci 螺旋で作る。
添字に対して高さ z が単調増加し、方位角が黄金角ずつ回る並びなので、
**番号で色分けすると螺旋がそのまま見える**。等立体角に並んでいることが
目で確かめられるので、色分けは番号（1 本目 → 最後）にしてある。

    cd geosim
    python view_directions.py            2000 本で開く
    python view_directions.py --rays 20000
"""

import numpy as np
import pyvista as pv

import sound_ray as sr
import view_model_gui as vg

TEXT_COLOR = "#d6dae2"

# 表示に耐える本数の上限。これを超えたら等間隔に間引いて描く
# （20 万本を線で描くと操作が重くなるだけで、密度は 2 万本でも十分わかる）
MAX_DRAWN = 20000


class DirectionPreview:
    """音線ベクトルを原点からの線分として描く。本数はスライダで変えられる。"""

    def __init__(self, plotter, total, drawn=None):
        self.plotter = plotter
        self.total = int(total)
        self.drawn = int(min(drawn or self.total, MAX_DRAWN, self.total))
        self.actor = None
        self.label = None
        # **全本数ぶんを一度だけ作っておく**。スライダを動かすたびに
        # 生成し直すと 20 万本では待たされる
        self.directions = sr.soundray_generator(self.total)
        self.rebuild(self.drawn, render=False)

    def _subset(self, count):
        """`count` 本を**等間隔に**抜く。

        先頭から取ってはいけない。Fibonacci 螺旋は添字に対して z が単調増加するので、
        先頭から取ると天頂付近の帽子状に偏る（球に見えなくなる）。
        """
        count = int(np.clip(count, 1, min(self.total, MAX_DRAWN)))
        index = np.unique(np.round(
            np.linspace(0, self.total - 1, count)).astype(int))
        return index

    def rebuild(self, count, render=True):
        index = self._subset(count)
        self.drawn = len(index)
        vectors = self.directions[index]

        # 原点 → 単位ベクトルの線分をまとめて 1 つの PolyData にする
        n = len(vectors)
        points = np.empty((2 * n, 3))
        points[0::2] = 0.0
        points[1::2] = vectors
        lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2),
                                 np.arange(1, 2 * n, 2)]).ravel()
        poly = pv.PolyData(points, lines=lines.astype(np.int64))
        # 色は音線の番号（1 本目 → 最後）。螺旋の並びがそのまま見える
        poly.point_data["ray"] = np.repeat(index.astype(float), 2)

        if self.actor is not None:
            self.plotter.remove_actor(self.actor, render=False)
        self.actor = self.plotter.add_mesh(
            poly, scalars="ray", cmap="hsv", line_width=1.0, lighting=False,
            scalar_bar_args={"title": "Ray number (1 - N)", "color": TEXT_COLOR,
                             "n_labels": 5, "fmt": "%.0f",
                             "position_x": 0.32, "position_y": 0.05,
                             "width": 0.5, "height": 0.05})
        self._refresh_label()
        if render:
            self.plotter.render()

    def _refresh_label(self):
        if self.label is None:
            return
        # 単位球の表面積 4π を本数で割ると、1 本あたりが受け持つ立体角になる。
        # 受音球に当たる確率の目安として使えるので添えておく
        solid_angle = 4.0 * np.pi / self.total
        text = (f"生成した音線 {self.total} 本\n"
                f"描いている本数 {self.drawn} 本\n"
                f"1 本の立体角 {solid_angle:.2e} sr")
        if self.drawn < self.total:
            text += f"\n（上限 {MAX_DRAWN} 本で間引き）"
        vg.set_actor_text(self.label, text)


def show(total=2000, drawn=None, title="音線の飛び方", off_screen=False,
         screenshot=None, window_size=(1100, 800), panel=None,
         save_dir=None):
    """音線ベクトルのプレビューを開く。"""
    want_panel = (not off_screen) if panel is None else bool(panel)
    plotter, panel = vg.make_plotter(title, window_size, off_screen,
                                     panel=want_panel, screen="directions")
    preview = DirectionPreview(plotter, total, drawn=drawn)

    # 音源の位置に球を置く（原点）。大きさは音線の長さ 1 に対する目安
    plotter.add_mesh(pv.Sphere(radius=0.03), color="#ff5f5f", lighting=False)
    plotter.add_axes(color=TEXT_COLOR)

    if panel is not None:
        panel.text(title, size=11, color=TEXT_COLOR)
        panel.text("音源から出る向きだけ。\n室形状は関係ありません", size=9)
        preview.label = panel.reserve_text(4)
        panel.heading("表示する本数")
        panel.slider("描く本数", [1, min(total, MAX_DRAWN)], preview.drawn,
                     lambda v: preview.rebuild(int(round(v))), fmt="%.0f")
        panel.heading("操作")
        if save_dir:
            vg.add_screenshot_key(plotter, save_dir, "音線の飛び方", key="g")
        panel.text("ドラッグ 回転 / ホイール 拡縮\n"
                   "z/x/c/v 視点   r リセット\n"
                   + ("g いまの画面を画像で保存\n" if save_dir else "")
                   + "q 閉じる", color="#7f8794")
        panel.relayout()
        preview._refresh_label()

    plotter.view_isometric()
    if off_screen:
        if screenshot:
            plotter.screenshot(screenshot)
        plotter.close()
        return preview
    vg.finish_window(plotter)
    plotter.show()
    return preview


def main():
    import argparse

    p = argparse.ArgumentParser(description="音線の飛び方を見る")
    p.add_argument("--rays", type=int, default=2000, help="生成する音線の本数")
    p.add_argument("--drawn", type=int, help="描く本数（既定は上限まで全部）")
    p.add_argument("--screenshot", help="画像に書き出して終了")
    a = p.parse_args()
    show(total=a.rays, drawn=a.drawn, off_screen=a.screenshot is not None,
         screenshot=a.screenshot)


if __name__ == "__main__":
    main()
