"""音線と音粒子の可視化（出力① G-1 / 出力② G-2）。

`ray_recorder.RayRecorder` が保存した npz を読んで、モデルに重ねて表示する。

  ① 音線の可視化   … 反射経路を折れ線で描く。**どの経路を通ってきたか**を見る
  ② 音粒子の可視化 … 離散化時間ごとに粒子が飛ぶ様子。**音がどう広がるか**を見る

同じデータから両方を作る。前者は経路の確認、後者は広がりの把握と、用途が違う
（docs/出力・可視化方針.md 参照）。

使い方:

    cd geosim
    # まず計算して軌跡を保存する
    python procedure.py ..\\test.dxf --absorption ..\\absorption.csv --out ..\\結果

    # ① 音線
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz
    # 受音した経路だけ
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --received-only
    # ② 音粒子（アニメーション）
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --mode particles
    # GIF に書き出す
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --mode particles --movie 広がり.gif

操作:
    共通      ドラッグ 回転 / ホイール 拡大縮小 / `r` 視点リセット / `q` 終了
              左上のチェックボックス レイヤの表示 ON/OFF
              **左の縦スライダ 不透明度** / `Tab` その対象を切り替え（すべて↔各レイヤ）
              `m` モデル表示の ON/OFF
    音粒子    `スペース` 再生・一時停止 / `←` `→` 1 コマ送り / `Home` 先頭へ /
              下の横スライダ 時刻を指定

壁の透過はレイヤごとに変えられる。起動時に決めておくなら
`--opacity`（全体）と `--layer-opacity "1=0.6,2=0.05"`（レイヤ別）。

※ 全体的な GUI の設計（G-7）は後回しにしている。ここは「見るための最小限」に留めてある。
"""

import argparse
import os

import numpy as np
import pyvista as pv

import read_dxffile as rd
import view_model_gui as vg

# 音線の色分けに使う指標
COLOR_MODES = ("energy", "time", "reflection", "ray")

TEXT_COLOR = "#d6dae2"
RECEIVED_COLOR = "#ffd166"      # 受音した経路の色


class RayLog:
    """`RayRecorder.save_npz()` が書いた軌跡データ。

    可変長データを「連結＋オフセット」で持っているので、そのままでは扱いにくい。
    ここで**行ごとに揃えた 2 次元配列**（足りないところは埋める）にしておくと、
    任意時刻の粒子位置を全音線ぶんまとめて計算できる。
    """

    def __init__(self, filename):
        data = np.load(filename)
        self.nodes = data["nodes"]                     # (N,3)
        self.distances = data["distances"]             # (N,)
        self.node_offsets = data["node_offsets"]       # (R+1,)
        self.mesh_ids = data["mesh_ids"]
        self.mesh_offsets = data["mesh_offsets"]
        self.receive_steps = data["receive_steps"]
        self.recv_offsets = data["recv_offsets"]
        self.ray_indexes = data["ray_indexes"]
        self.directions = data["directions"]
        self.terminations = data["terminations"]
        self.total_rays = int(data["total_rays"])
        self.stride = int(data["stride"])
        self.sound_velocity = float(data["sound_velocity"])
        self.energies = data["energies"] if "energies" in data.files else None

        self.ray_count = len(self.node_offsets) - 1
        self.node_counts = np.diff(self.node_offsets)
        self.reflection_counts = np.diff(self.mesh_offsets)
        self.received = np.diff(self.recv_offsets) > 0

        self._build_padded()

    def _build_padded(self):
        """行ごとに揃えた配列を作る（粒子位置の一括計算用）。"""
        width = int(self.node_counts.max())
        self.width = width
        rows = self.ray_count

        # 距離は末尾を +inf で埋める。こうすると「d 以下の要素数」を数えるだけで
        # 区間の添字が求まり、音線ごとの searchsorted が要らなくなる
        self.pad_distance = np.full((rows, width), np.inf)
        self.pad_nodes = np.zeros((rows, width, 3))
        for i in range(rows):
            start, stop = self.node_offsets[i], self.node_offsets[i + 1]
            n = stop - start
            self.pad_distance[i, :n] = self.distances[start:stop]
            self.pad_nodes[i, :n] = self.nodes[start:stop]
            self.pad_nodes[i, n:] = self.nodes[stop - 1]

        self.total_distance = self.pad_distance[
            np.arange(rows), self.node_counts - 1]
        self.max_time = float(self.total_distance.max()) / self.sound_velocity

    # ---- 抜き出し ----------------------------------------------------

    def selection(self, received_only=False, max_rays=None):
        """描画する音線の添字を選ぶ。

        max_rays を超える場合は**等間隔に間引く**（記録時と同じ考え方）。
        ランダムに抜くと見た目の密度が変わってしまうため。
        """
        keep = self.received.copy() if received_only else np.ones(self.ray_count, bool)
        index = np.nonzero(keep)[0]
        if max_rays is not None and len(index) > max_rays:
            index = index[:: max(1, len(index) // max_rays)][:max_rays]
        return index

    # ---- ① 音線（折れ線）--------------------------------------------

    def line_polydata(self, index=None, colour="energy", band=None,
                      max_reflection=None):
        """選んだ音線を折れ線の PolyData にする。

        点ごとにスカラー値を付けるので、1 本の線の中で色が変わる
        （音源側は明るく、反射を重ねるほど暗く、など）。

        max_reflection … **描く反射回数の上限**。指定するとそこで折れ線を打ち切る。
            残響の後半まで全部描くと線が重なって何も読めなくなるので、
            初期反射だけを見たいときに使う。
            「その回数以下の音線だけを残す」のではなく「途中まで描く」動作なのが要点。
        """
        if index is None:
            index = np.arange(self.ray_count)
        if len(index) == 0:
            raise ValueError("描画する音線がありません")

        points = []
        lines = []
        scalar = []
        cursor = 0
        for i in index:
            start, stop = self.node_offsets[i], self.node_offsets[i + 1]
            if max_reflection is not None:
                stop = min(stop, start + max_reflection + 1)
            n = stop - start
            if n < 2:
                continue
            points.append(self.nodes[start:stop])
            lines.append(np.concatenate([[n], np.arange(cursor, cursor + n)]))
            scalar.append(self._scalar(i, start, stop, colour, band))
            cursor += n
        if not points:
            raise ValueError("描画する区間がありません")

        poly = pv.PolyData(np.concatenate(points),
                           lines=np.concatenate(lines).astype(np.int64))
        poly.point_data[colour] = np.concatenate(scalar)
        return poly

    def _scalar(self, ray, start, stop, colour, band):
        if colour == "time":
            return self.distances[start:stop] / self.sound_velocity * 1000.0   # ms
        if colour == "reflection":
            return np.arange(stop - start, dtype=float)
        if colour == "ray":
            return np.full(stop - start, float(self.ray_indexes[ray]))
        if colour == "energy":
            if self.energies is None:
                return np.zeros(stop - start)
            energy = self.energies[start:stop]
            energy = energy[:, band] if band is not None else energy.mean(axis=1)
            with np.errstate(divide="ignore"):
                return 10.0 * np.log10(np.maximum(energy, 1e-12))
        raise ValueError(f"colour は {COLOR_MODES} のいずれかです: {colour!r}")

    # ---- ② 音粒子 ----------------------------------------------------

    def positions_at(self, time, index=None):
        """時刻 time [s] における音粒子の位置を**まとめて**求める。

        戻り値 (位置 (K,3), 元の音線の添字 (K,))。
        まだ出発していない／すでに消えた音線は含まれない。

        `RayTrajectory.position_at()` の一括版。毎フレーム音線ごとに呼ぶと
        Python のループが効いてくるので、行を揃えた配列で一度に計算する。
        """
        if index is None:
            index = np.arange(self.ray_count)
        distance = time * self.sound_velocity

        alive = distance <= self.total_distance[index]
        rows = index[alive]
        if len(rows) == 0:
            return np.zeros((0, 3)), rows

        pad = self.pad_distance[rows]
        # 「distance 以下の要素数 - 1」が区間の添字。inf 埋めなので末尾は数えられない
        k = np.count_nonzero(pad <= distance, axis=1) - 1
        k = np.clip(k, 0, self.node_counts[rows] - 2)

        order = np.arange(len(rows))
        d0 = pad[order, k]
        d1 = pad[order, k + 1]
        span = d1 - d0
        weight = np.where(span > 0.0, (distance - d0) / np.where(span > 0.0, span, 1.0), 0.0)

        n0 = self.pad_nodes[rows, k]
        n1 = self.pad_nodes[rows, k + 1]
        return n0 + weight[:, None] * (n1 - n0), rows

    def energy_at(self, time, index=None, band=None):
        """時刻 time における粒子のエネルギー（dB）。反射のたびに階段状に下がる。"""
        if self.energies is None:
            return None
        position, rows = self.positions_at(time, index)
        if len(rows) == 0:
            return np.zeros(0)
        distance = time * self.sound_velocity
        pad = self.pad_distance[rows]
        k = np.clip(np.count_nonzero(pad <= distance, axis=1) - 1,
                    0, self.node_counts[rows] - 1)
        flat = self.node_offsets[rows] + k
        energy = self.energies[flat]
        energy = energy[:, band] if band is not None else energy.mean(axis=1)
        with np.errstate(divide="ignore"):
            return 10.0 * np.log10(np.maximum(energy, 1e-12))

    def summary(self):
        return (f"音線 {self.ray_count} 本（元 {self.total_rays} 本を {self.stride} 本おき）"
                f" / 受音 {int(self.received.sum())} 本"
                f" / 反射 0〜{int(self.reflection_counts.max())} 回"
                f" / 最大時刻 {self.max_time * 1000:.1f} ms"
                f" / 音速 {self.sound_velocity:.1f} m/s")


# ------------------------------------------------------------------------------
# ① 音線の可視化
# ------------------------------------------------------------------------------

def add_rays(plotter, raylog, index=None, colour="energy", band=None,
             line_width=2.0, cmap="plasma", highlight_received=True, opacity=0.8,
             max_reflection=None):
    """折れ線として音線を描く。

    ※ カラーバーの見出しは VTK の既定フォントで描かれるため**日本語が出せない**。
      本文側は日本語フォントを指定しているが、ここだけは英字にしてある。
    """
    poly = raylog.line_polydata(index, colour=colour, band=band,
                                max_reflection=max_reflection)
    label = {"energy": "Energy [dB]", "time": "Time [ms]",
             "reflection": "Reflection", "ray": "Ray index"}[colour]
    actor = plotter.add_mesh(poly, scalars=colour, cmap=cmap, line_width=line_width,
                             lighting=False, opacity=opacity,
                             scalar_bar_args={"title": label, "color": TEXT_COLOR,
                                              "n_labels": 5, "fmt": "%.1f",
                                              "position_x": 0.32, "position_y": 0.02,
                                              "width": 0.5, "height": 0.05})
    # 受音した経路を目立たせる。ただし受音経路しか描いていないときは意味がないので出さない
    received = None
    chosen = np.arange(raylog.ray_count) if index is None else np.asarray(index)
    hit = chosen[raylog.received[chosen]]
    if highlight_received and len(hit) and len(hit) < len(chosen):
        received = plotter.add_mesh(
            raylog.line_polydata(hit, colour="reflection",
                                 max_reflection=max_reflection),
            color=RECEIVED_COLOR, line_width=line_width + 2.0, lighting=False,
            opacity=1.0, show_scalar_bar=False)
    return actor, received


# ------------------------------------------------------------------------------
# ② 音粒子の可視化
# ------------------------------------------------------------------------------

class ParticleAnimation:
    """音粒子アニメーションの状態を持つ。

    毎フレーム全音線の位置を計算し直して 1 つの点群を差し替える。
    音線ごとに actor を作ると本数分の描画呼び出しになるので、点群 1 つにまとめている。
    """

    def __init__(self, plotter, raylog, index=None, frames=240, band=None,
                 point_size=9.0, cmap="plasma", trail=0.0):
        self.plotter = plotter
        self.raylog = raylog
        self.index = np.arange(raylog.ray_count) if index is None else np.asarray(index)
        self.frames = int(frames)
        self.band = band
        self.trail = float(trail)          # 尾を引く長さ [s]
        self.step = 0
        self.playing = True

        self.times = np.linspace(0.0, raylog.max_time, self.frames)
        self.cloud = pv.PolyData(np.zeros((1, 3)))
        self.cloud.point_data["energy"] = np.zeros(1)
        self._dead = np.full(1, -120.0)

        energies = raylog.energy_at(self.times[len(self.times) // 3], self.index, band)
        low = float(np.min(energies)) if energies is not None and len(energies) else -60.0
        self.actor = plotter.add_mesh(
            self.cloud, scalars="energy", cmap=cmap, clim=(max(low, -60.0), 0.0),
            point_size=point_size, render_points_as_spheres=True, lighting=False,
            scalar_bar_args={"title": "Energy [dB]", "color": TEXT_COLOR,
                             "n_labels": 5, "fmt": "%.0f",
                             "position_x": 0.32, "position_y": 0.11,
                             "width": 0.5, "height": 0.05})
        # 位置を座標で渡すと vtkTextActor が返る（文字列だと CornerAnnotation になり
        # 差し替えの API が違うので、毎フレーム書き換えるここでは座標を使う）
        width, height = plotter.window_size
        self.label = plotter.add_text(" ", position=(width - 430, height - 150),
                                      font_size=11, color=TEXT_COLOR,
                                      font_file=vg.japanese_font())
        self.update(0)

    @staticmethod
    def _set_text(actor, text):
        if hasattr(actor, "SetInput"):
            actor.SetInput(text)
        else:                       # CornerAnnotation の場合（3 = 右上）
            actor.SetText(3, text)

    def update(self, step):
        self.step = int(step) % self.frames
        t = self.times[self.step]
        position, rows = self.raylog.positions_at(t, self.index)
        energy = self.raylog.energy_at(t, self.index, self.band)

        if len(position) == 0:
            position = np.zeros((1, 3))
            energy = np.full(1, -120.0)
        elif energy is None:
            energy = np.zeros(len(position))

        # ★点の数が変わるので、座標だけでなく**頂点セルも作り直す**。
        #   points を差し替えただけだと古いセル（＝最初の 1 点）しか描かれない
        count = len(position)
        self.cloud.points = position
        self.cloud.verts = np.column_stack([
            np.ones(count, dtype=np.int64),
            np.arange(count, dtype=np.int64)]).ravel()
        self.cloud.point_data["energy"] = energy
        self._set_text(self.label,
                       f"{t * 1000:7.2f} ms   粒子 {len(rows):5d} / {len(self.index)}"
                       f"   [{self.step + 1}/{self.frames}]")
        self.plotter.render()

    def advance(self):
        if self.playing:
            self.update(self.step + 1)

    def toggle(self):
        self.playing = not self.playing


def animate(plotter, raylog, index=None, frames=240, band=None, point_size=9.0):
    """音粒子アニメーションの部品（スライダ・キー操作）を組み立てる。

    実際にコマを進めるのは `run_animation()`。
    """
    animation = ParticleAnimation(plotter, raylog, index=index, frames=frames,
                                  band=band, point_size=point_size)

    plotter.add_key_event("space", animation.toggle)
    plotter.add_key_event("Right", lambda: animation.update(animation.step + 1))
    plotter.add_key_event("Left", lambda: animation.update(animation.step - 1))
    plotter.add_key_event("Home", lambda: animation.update(0))

    # 見出しは VTK の既定フォントで描かれるので英字にする（日本語は豆腐になる）
    animation.slider = plotter.add_slider_widget(
        lambda value: animation.update(int(round(value))),
        [0, frames - 1], value=0, title="frame", pointa=(0.32, 0.05),
        pointb=(0.90, 0.05), style="modern", fmt="%.0f",
        color=TEXT_COLOR, title_color=TEXT_COLOR)
    return animation


def run_animation(plotter, animation, interval=30):
    """ウィンドウを開いてアニメーションを回す（閉じられるまで戻らない）。

    【なぜ VTK のタイマーを使わないか】2026-08-14
    最初は `Plotter.add_timer_event()` で駒を進めていたが、**発火しなかった**
    （20 回 × 50 ms なら 1 秒で終わるはずが、60 秒経っても 0 回）。
    VTK の `CreateRepeatingTimer` は interactor の初期化後でないと効かないが、
    pyvista の API は `show()` より前に呼ぶ形になっているためだと思われる。
    コマ送り（キー操作）だけ効いて再生されない、という症状になっていた。

    そこで **`show(interactive_update=True)` で非ブロッキング表示にして、
    こちらのループから `update()` を呼ぶ**方式にした。
    pyvista のアニメーション例でも使われている確実な作り方で、
    回転などの操作は `update()` の中で処理される。
    """
    closed = {"flag": False}
    plotter.add_key_event("q", lambda: closed.__setitem__("flag", True))
    try:
        plotter.iren.add_observer("ExitEvent",
                                  lambda *a: closed.__setitem__("flag", True))
    except Exception:
        pass

    plotter.show(interactive_update=True, auto_close=False)
    while not closed["flag"]:
        if getattr(plotter, "_closed", False) or plotter.render_window is None:
            break
        if animation.playing:
            animation.update(animation.step + 1)
        try:
            plotter.update(interval)
        except (RuntimeError, AttributeError):
            break
    plotter.close()
    return animation


def save_movie(raylog, model, filename, index=None, frames=240, band=None,
               point_size=9.0, opacity=0.12, window_size=(960, 720), duration=40):
    """音粒子アニメーションを GIF に書き出す。

    追加の依存を増やさないよう、コマを画像として集めて Pillow で GIF にする
    （Pillow は matplotlib の依存で既に入っている）。
    """
    from PIL import Image

    plotter = vg.build_plotter(model, title="音粒子", off_screen=True,
                               show_normals=False, opacity=opacity,
                               window_size=window_size)
    animation = ParticleAnimation(plotter, raylog, index=index, frames=frames,
                                  band=band, point_size=point_size)
    plotter.view_isometric()

    images = []
    for step in range(frames):
        animation.update(step)
        images.append(Image.fromarray(plotter.screenshot(return_img=True)))
    plotter.close()

    images[0].save(filename, save_all=True, append_images=images[1:],
                   duration=duration, loop=0, optimize=True)
    return filename


# ------------------------------------------------------------------------------
# 入口
# ------------------------------------------------------------------------------

def view(dxf_path, raylog_path, mode="rays", absorption=None, unit=None,
         orient_normals="cad", received_only=False, max_rays=None,
         max_reflection=None, colour="energy", band=None, frames=240,
         opacity=0.12, layer_opacity=None, movie=None, point_size=9.0,
         screenshot=None, interval=30):
    """モデルの上に音線または音粒子を重ねて表示する。"""
    model = rd.read_model(dxf_path, unit=unit, absorption_table=absorption,
                          orient_normals=orient_normals)
    raylog = RayLog(raylog_path)
    print(f"[view_rays] {raylog.summary()}")

    index = raylog.selection(received_only=received_only, max_rays=max_rays)
    if len(index) == 0:
        raise ValueError("条件に合う音線がありません（--received-only を外してみてください）")
    print(f"[view_rays] 描画する音線 {len(index)} 本")

    if movie is not None:
        path = save_movie(raylog, model, movie, index=index, frames=frames,
                          band=band, point_size=point_size, opacity=opacity)
        print(f"[view_rays] 動画を書き出しました: {path}")
        return raylog

    base = os.path.splitext(os.path.basename(dxf_path))[0]
    title = f"{base} {'音粒子' if mode == 'particles' else '音線'}"
    off_screen = screenshot is not None
    plotter = vg.build_plotter(model, title=title, off_screen=off_screen,
                               show_normals=False, opacity=opacity,
                               layer_opacity=layer_opacity, show_summary=False)
    plotter.add_text(raylog.summary().replace(" / ", "\n"), position="upper_right",
                     font_size=9, color="#9aa2b1", font_file=vg.japanese_font())

    font = vg.japanese_font()
    animation = None
    if mode == "rays":
        add_rays(plotter, raylog, index=index, colour=colour, band=band,
                 max_reflection=max_reflection)
    elif mode == "particles":
        if off_screen:
            ParticleAnimation(plotter, raylog, index=index, frames=frames,
                              band=band, point_size=point_size).update(frames // 3)
        else:
            animation = animate(plotter, raylog, index=index, frames=frames,
                                band=band, point_size=point_size)
            plotter.add_text("スペース 再生/停止   ← → コマ送り   Home 先頭",
                             position=(14, 14), font_size=9, color="#7f8794",
                             font_file=font)
    else:
        raise ValueError(f"mode は 'rays' か 'particles' です: {mode!r}")

    # 壁の透過はレイヤごとに変えられる（Tab で対象切替、m で表示 ON/OFF）
    if not off_screen:
        vg.add_opacity_control(plotter, font=font)

    plotter.view_isometric()
    if off_screen:
        plotter.screenshot(screenshot)
        plotter.close()
        print(f"[view_rays] 画像を書き出しました: {screenshot}")
    elif animation is not None:
        run_animation(plotter, animation, interval=interval)
    else:
        plotter.show()
    return raylog


def parse_layer_opacity(text):
    """`"1=0.05,2=0.3"` のような指定を辞書にする。"""
    if not text:
        return None
    result = {}
    for item in text.split(","):
        if "=" not in item:
            raise ValueError(f"--layer-opacity の書式は レイヤ名=値 です: {item!r}")
        name, value = item.rsplit("=", 1)
        result[name.strip()] = float(value)
    return result


def main():
    p = argparse.ArgumentParser(description="音線・音粒子の可視化（G-1 / G-2）")
    p.add_argument("dxf", help="室形状の DXF")
    p.add_argument("raylog", help="RayRecorder が保存した npz")
    p.add_argument("--mode", default="rays", choices=["rays", "particles"],
                   help="rays=音線の折れ線 / particles=音粒子のアニメーション")
    p.add_argument("--absorption", help="吸音率 CSV（モデル表示用。省略可）")
    p.add_argument("--unit", help="'mm' / 'm' など")
    p.add_argument("--orient-normals", default="cad", choices=["cad", "flip", "shells"])
    p.add_argument("--received-only", action="store_true",
                   help="受音した経路だけを描く")
    p.add_argument("--max-rays", type=int, default=80,
                   help="描く音線の本数の上限（等間隔に間引く）。"
                        "多すぎると線が重なって読めなくなる。既定 80")
    p.add_argument("--max-reflection", type=int,
                   help="描く反射回数の上限。そこで折れ線を打ち切る（初期反射だけ見たいとき）")
    p.add_argument("--color", default="energy", choices=list(COLOR_MODES),
                   help="音線の色分けに使う指標")
    p.add_argument("--band", type=int, help="エネルギーで色分けするときのバンド番号（0 始まり）")
    p.add_argument("--frames", type=int, default=240, help="アニメーションのコマ数")
    p.add_argument("--point-size", type=float, default=9.0, help="音粒子の大きさ")
    p.add_argument("--opacity", type=float, default=0.12,
                   help="モデルの不透明度（0=透明, 1=不透明）。既定 0.12")
    p.add_argument("--layer-opacity",
                   help="レイヤごとの不透明度。例 \"1=0.6,2=0.05\"。"
                        "指定しないレイヤは --opacity を使う")
    p.add_argument("--interval", type=int, default=30,
                   help="アニメーションのコマ間隔 [ms]")
    p.add_argument("--movie", help="GIF に書き出す（ウィンドウを開かない）")
    p.add_argument("--screenshot", help="静止画に書き出す（ウィンドウを開かない）")
    a = p.parse_args()

    view(a.dxf, a.raylog, mode=a.mode, absorption=a.absorption, unit=a.unit,
         orient_normals=a.orient_normals, received_only=a.received_only,
         max_rays=a.max_rays, max_reflection=a.max_reflection, colour=a.color,
         band=a.band, frames=a.frames, opacity=a.opacity,
         layer_opacity=parse_layer_opacity(a.layer_opacity), movie=a.movie,
         point_size=a.point_size, screenshot=a.screenshot, interval=a.interval)


if __name__ == "__main__":
    main()
