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

    # ①② 両方を同居させる（既定。Tab で切り替え）
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz
    # 受音した経路だけ
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --received-only
    # 片方だけにする
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --mode particles
    # GIF に書き出す
    python view_rays.py ..\\test.dxf ..\\結果\\test_raylog.npz --mode particles --movie 広がり.gif

**画面は左右に分かれている**（2026-08-15）。左 1/4 が操作パネル、右 3/4 が 3D 表示。
以前は 3D の上に操作系を重ねていたため、モデルが見えなくなるうえ
要素どうしが重なって読めなかった（ユーザー指摘）。

左パネルでできること:
    レイヤの表示 ON/OFF（チェックボックス）
    不透明度（スライダ）。`o` で対象を すべて↔各レイヤ に切り替え
    **表示する音線の本数・音粒子の数**（スライダ）
    音線の情報・時刻・操作説明

操作:
    共通      ドラッグ 回転 / ホイール 拡大縮小 / `z` `x` `c` `v` 視点 / `r` リセット / `q` 終了
              **`Tab` 音線 ↔ 音粒子の切り替え**（ウィンドウを閉じずに見比べられる）
              `o` 不透明度の対象を切り替え / `m` モデル表示の ON/OFF
    音粒子    `スペース` 再生・一時停止 / `←` `→` 1 コマ送り / `Home` 先頭へ /
              下の横スライダ 時刻を指定

起動時の本数は `--max-rays`（既定 60）と `--max-particles`（既定は候補すべて）。
候補そのものの上限は `--pool`（既定 2000）＝スライダの上限になる。

★**本数を減らすときは等間隔に抜く**（`even_subset`）。音線は Fibonacci 螺旋で
作っており添字に対して z が単調増加するので、先頭から取ると天頂付近に偏り、
「片方向にしか飛んでいない」ように見えてしまう。

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

# カラーバーの見出し。VTK の既定フォントで描かれるので**日本語が出せない**（英字にしてある）。
# 音線と音粒子を同じウィンドウに同居させるので、**見出しは重ならないようにする**
# （pyvista はカラーバーを見出しで管理しているため、同名だと 1 つにまとめられてしまう）。
RAY_BAR_TITLE = {"energy": "Ray energy [dB]", "time": "Ray time [ms]",
                 "reflection": "Reflection", "ray": "Ray number (1 - N)"}

# 音線の色分けに使う配色。`ray`（音線の番号）は**1 本目から最後まで**を
# 一巡する色相にすると、全方向へ均等に散っているかが目で確かめられる。
# 明るさが単調に変わる plasma だと「どこが 1 本目か」が分かりにくい
RAY_CMAP = {"ray": "hsv"}
PARTICLE_BAR_TITLE = "Particle energy [dB]"


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

        order = np.arange(rows)
        self.total_distance = self.pad_distance[order, self.node_counts - 1]
        self.last_node = self.pad_nodes[order, self.node_counts - 1]
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

    def positions_masked(self, time, index=None):
        """時刻 time [s] の粒子位置を、`index` と同じ並びで返す。

        戻り値 (位置 (len(index),3), 生きているか (len(index),) の bool)。
        すでに消えた音線の位置は**最後の節点**に置く（描画側で隠す前提）。

        `RayTrajectory.position_at()` の一括版。毎フレーム音線ごとに呼ぶと
        Python のループが効いてくるので、行を揃えた配列で一度に計算する。
        """
        if index is None:
            index = np.arange(self.ray_count)
        distance = time * self.sound_velocity
        alive = distance <= self.total_distance[index]

        position = self.last_node[index].copy()
        rows = index[alive]
        if len(rows) == 0:
            return position, alive

        pad = self.pad_distance[rows]
        # 「distance 以下の要素数 - 1」が区間の添字。inf 埋めなので末尾は数えられない
        k = np.count_nonzero(pad <= distance, axis=1) - 1
        k = np.clip(k, 0, self.node_counts[rows] - 2)

        order = np.arange(len(rows))
        d0 = pad[order, k]
        d1 = pad[order, k + 1]
        span = d1 - d0
        weight = np.where(span > 0.0,
                          (distance - d0) / np.where(span > 0.0, span, 1.0), 0.0)

        n0 = self.pad_nodes[rows, k]
        n1 = self.pad_nodes[rows, k + 1]
        position[alive] = n0 + weight[:, None] * (n1 - n0)
        return position, alive

    def positions_at(self, time, index=None):
        """時刻 time における、**生きている粒子だけ**の位置。

        戻り値 (位置 (K,3), 元の音線の添字 (K,))。
        """
        if index is None:
            index = np.arange(self.ray_count)
        position, alive = self.positions_masked(time, index)
        return position[alive], index[alive]

    def energy_masked(self, time, index=None, band=None):
        """時刻 time の粒子エネルギー（dB）を `index` と同じ並びで返す。

        消えた粒子は **NaN**。描画側で `nan_opacity=0` にして隠すため。
        """
        if index is None:
            index = np.arange(self.ray_count)
        if self.energies is None:
            return np.zeros(len(index))

        distance = time * self.sound_velocity
        alive = distance <= self.total_distance[index]
        result = np.full(len(index), np.nan)
        rows = index[alive]
        if len(rows) == 0:
            return result

        pad = self.pad_distance[rows]
        k = np.clip(np.count_nonzero(pad <= distance, axis=1) - 1,
                    0, self.node_counts[rows] - 1)
        energy = self.energies[self.node_offsets[rows] + k]
        energy = energy[:, band] if band is not None else energy.mean(axis=1)
        with np.errstate(divide="ignore"):
            result[alive] = 10.0 * np.log10(np.maximum(energy, 1e-12))
        return result

    def energy_at(self, time, index=None, band=None):
        """時刻 time における、生きている粒子だけのエネルギー（dB）。"""
        if self.energies is None:
            return None
        if index is None:
            index = np.arange(self.ray_count)
        result = self.energy_masked(time, index, band)
        return result[~np.isnan(result)]

    def summary(self):
        return (f"音線 {self.ray_count} 本（元 {self.total_rays} 本を {self.stride} 本おき）"
                f" / 受音 {int(self.received.sum())} 本"
                f" / 反射 0〜{int(self.reflection_counts.max())} 回"
                f" / 最大時刻 {self.max_time * 1000:.1f} ms"
                f" / 音速 {self.sound_velocity:.1f} m/s")


# ------------------------------------------------------------------------------
# ① 音線の可視化
# ------------------------------------------------------------------------------

def even_subset(index, count):
    """`index` から `count` 本を**等間隔に**抜く。

    ★先頭から `count` 本を取ってはいけない。音線は Fibonacci 螺旋で作っており、
    **添字に対して z が単調増加**するので、先頭から取ると天頂付近の帽子状に偏る
    （「片方向にしか飛んでいない」ように見える原因になる）。
    等間隔に抜けば全方向に散ったままになる。
    """
    index = np.asarray(index)
    count = int(np.clip(count, 1, len(index)))
    if count >= len(index):
        return index
    picked = np.unique(np.round(np.linspace(0, len(index) - 1, count)).astype(int))
    return index[picked]


class RayDisplay:
    """音線の折れ線。**表示本数をあとから変えられる**ようにまとめたもの。

    本数を変えるたびに折れ線を作り直す（点群と違って線は本数で構造が変わるため）。
    数百本までなら作り直しても一瞬で終わる。
    """

    def __init__(self, plotter, raylog, index, colour="energy", band=None,
                 line_width=2.0, cmap="plasma", highlight_received=True,
                 opacity=0.8, max_reflection=None, count=None):
        self.plotter = plotter
        self.raylog = raylog
        self.pool = np.asarray(index)
        self.colour = colour
        self.band = band
        self.line_width = line_width
        self.cmap = RAY_CMAP.get(colour, cmap)
        self.highlight_received = highlight_received
        self.opacity = opacity
        self.max_reflection = max_reflection
        # スライダの上限。これ以上は「打ち切らない」と同じなので None にする
        self.reflection_limit = int(raylog.reflection_counts.max())
        self.count = int(count or len(self.pool))
        self.actors = []
        self.visible = True
        self.rebuild(self.count, render=False)

    def set_max_reflection(self, value, render=True):
        """描く反射回数の上限を変える。

        0（＝1 区間だけ）にすると**音源から最初に当たるまで**だけが描かれ、
        「どの向きへ音線を飛ばしたか」がそのまま見える。
        上げていくと反射のたびに折れ線が伸びる。
        全反射まで描くと室内が線で埋まって読めなくなるので、ここで加減する。
        """
        value = int(round(value))
        self.max_reflection = None if value >= self.reflection_limit else max(1, value)
        self.rebuild(self.count, render=render)

    def rebuild(self, count, render=True):
        self.count = int(np.clip(count, 1, len(self.pool)))
        for actor in self.actors:
            self.plotter.remove_actor(actor, render=False)
        self.actors = [a for a in add_rays(
            self.plotter, self.raylog, index=even_subset(self.pool, self.count),
            colour=self.colour, band=self.band, line_width=self.line_width,
            cmap=self.cmap, highlight_received=self.highlight_received,
            opacity=self.opacity, max_reflection=self.max_reflection) if a is not None]
        for actor in self.actors:
            actor.SetVisibility(self.visible)
        if render:
            self.plotter.render()

    def set_visible(self, flag):
        self.visible = bool(flag)
        for actor in self.actors:
            actor.SetVisibility(self.visible)


def add_rays(plotter, raylog, index=None, colour="energy", band=None,
             line_width=2.0, cmap="plasma", highlight_received=True, opacity=0.8,
             max_reflection=None):
    """折れ線として音線を描く。

    ※ カラーバーの見出しは VTK の既定フォントで描かれるため**日本語が出せない**。
      本文側は日本語フォントを指定しているが、ここだけは英字にしてある。
    """
    poly = raylog.line_polydata(index, colour=colour, band=band,
                                max_reflection=max_reflection)
    label = RAY_BAR_TITLE[colour]
    actor = plotter.add_mesh(poly, scalars=colour, cmap=cmap, line_width=line_width,
                             lighting=False, opacity=opacity,
                             scalar_bar_args={"title": label, "color": TEXT_COLOR,
                                              "n_labels": 5, "fmt": "%.1f",
                                              "position_x": 0.32, "position_y": 0.11,
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

    毎フレーム全音線の位置を計算し直して 1 つの点群を更新する。
    音線ごとに actor を作ると本数分の描画呼び出しになるので、点群 1 つにまとめている。

    【トポロジは変えない】2026-08-14
    最初は「生きている粒子だけ」を点群に入れていたので、フレームごとに点数が変わり
    **頂点セルを毎回作り直していた**。動きはするが、毎フレーム VTK 側で
    ジオメトリを組み直すことになり無駄が大きい。

    そこで **点の数を音線の本数で固定**し、毎フレーム書き換えるのは座標とスカラーだけにした。
    消えた粒子はスカラーを NaN にし、`nan_opacity=0` で見えなくしている。
    6000 フレーム連続で 60 fps を確認済み。
    """

    def __init__(self, plotter, raylog, index=None, frames=240, band=None,
                 point_size=9.0, cmap="plasma", label=None):
        self.plotter = plotter
        self.raylog = raylog
        self.index = np.arange(raylog.ray_count) if index is None else np.asarray(index)
        self.frames = int(frames)
        self.band = band
        self.step = 0
        self.playing = True
        # 表示する粒子の数。点の数（＝トポロジ）は変えず、**表示しないぶんは
        # エネルギーを NaN にして隠す**（`nan_opacity=0`）。作り直しが起きない
        self._shown = np.ones(len(self.index), dtype=bool)

        self.times = np.linspace(0.0, raylog.max_time, self.frames)

        # 点数は最初に決めたら変えない（上記の理由）
        count = len(self.index)
        self.cloud = pv.PolyData(raylog.last_node[self.index].copy())
        self.cloud.point_data["energy"] = np.full(count, np.nan)

        energies = raylog.energy_at(self.times[len(self.times) // 3], self.index, band)
        low = float(np.min(energies)) if energies is not None and len(energies) else -60.0
        self.actor = plotter.add_mesh(
            self.cloud, scalars="energy", cmap=cmap, clim=(max(low, -60.0), 0.0),
            nan_opacity=0.0,
            point_size=point_size, render_points_as_spheres=True, lighting=False,
            scalar_bar_args={"title": PARTICLE_BAR_TITLE, "color": TEXT_COLOR,
                             "n_labels": 5, "fmt": "%.0f",
                             "position_x": 0.32, "position_y": 0.11,
                             "width": 0.5, "height": 0.05})
        # 時刻と粒子数の表示。**左パネルに置くのが既定**（`label` で渡す）。
        # パネルが無いときだけ 3D の左下に出す。
        # ウィンドウ幅からの相対位置で右上に置いていたが、レンダラを左右に分けたら
        # ビューポートからはみ出して切れてしまった
        self.label = label or plotter.add_text(
            " ", position=(14, 14), font_size=11, color=TEXT_COLOR,
            font_file=vg.japanese_font())
        self.update(0)

    @staticmethod
    def _set_text(actor, text):
        if hasattr(actor, "SetInput"):
            actor.SetInput(text)
        else:                       # CornerAnnotation の場合（3 = 右上）
            actor.SetText(3, text)

    def update(self, step, render=True):
        self.step = int(step) % self.frames
        t = self.times[self.step]

        position, alive = self.raylog.positions_masked(t, self.index)
        energy = self.raylog.energy_masked(t, self.index, self.band)
        energy[~self._shown] = np.nan       # 表示数を絞ったぶんを隠す

        # ★座標とスカラーだけを書き換える（頂点セルはそのまま）。
        #   点数を変えると毎フレーム VTK 側でジオメトリを組み直すことになる
        self.cloud.points[:] = position
        self.cloud.point_data["energy"][:] = energy
        self.cloud.Modified()

        visible = int(np.count_nonzero(alive & self._shown))
        self._set_text(self.label,
                       f"{t * 1000:.1f} ms\n"
                       f"粒子 {visible}/{int(self._shown.sum())}"
                       f"  [{self.step + 1}/{self.frames}]")
        if render:
            self.plotter.render()

    def set_count(self, count, render=True):
        """表示する粒子の数を変える。**等間隔に抜く**ので分布の偏りは出ない。"""
        keep = even_subset(np.arange(len(self.index)), count)
        self._shown = np.zeros(len(self.index), dtype=bool)
        self._shown[keep] = True
        self.update(self.step, render=render)

    def advance(self):
        if self.playing:
            self.update(self.step + 1)

    def toggle(self):
        self.playing = not self.playing


def animate(plotter, raylog, index=None, frames=240, band=None, point_size=9.0,
            label=None):
    """音粒子アニメーションの部品（スライダ・キー操作）を組み立てる。

    実際にコマを進めるのは `run_animation()`。
    """
    animation = ParticleAnimation(plotter, raylog, index=index, frames=frames,
                                  band=band, point_size=point_size, label=label)

    plotter.add_key_event("space", animation.toggle)
    plotter.add_key_event("Right", lambda: animation.update(animation.step + 1))
    plotter.add_key_event("Left", lambda: animation.update(animation.step - 1))
    plotter.add_key_event("Home", lambda: animation.update(0))

    # 見出しは VTK の既定フォントで描かれるので英字にする（日本語は豆腐になる）
    animation.slider = plotter.add_slider_widget(
        lambda value: animation.update(int(round(value))),
        [0, frames - 1], value=0, title="frame", pointa=(0.32, 0.04),
        pointb=(0.90, 0.04), style="modern", fmt="%.0f",
        color=TEXT_COLOR, title_color=TEXT_COLOR)
    return animation


# ------------------------------------------------------------------------------
# ①② の切り替え表示
# ------------------------------------------------------------------------------

class RayParticleView:
    """音線と音粒子を**同じウィンドウに同居させ、Tab で切り替える**。

    以前は `--mode` でどちらかを選ぶ形だったので、見比べるにはウィンドウを
    一度閉じてコマンドを打ち直す必要があった。同じ軌跡データから作る 2 つの見せ方なので、
    両方を組み立てておいて**表示を切り替えるだけ**にした（作り直しは起きない）。

    切り替えるのは以下のひとまとまり。

    | | 音線 | 音粒子 |
    |---|---|---|
    | 本体 | 折れ線 ＋ 受音経路のハイライト | 点群 |
    | カラーバー | Ray … | Particle energy |
    | 下の横スライダ | なし | 時刻 |
    | 文字 | なし | 時刻・粒子数 |

    音粒子側を隠している間は `playing` を落として**コマを進めない**。
    進めたままだと、見えていないのに毎フレーム位置計算が走ることになる。
    """

    MODES = ("rays", "particles")

    def __init__(self, plotter, animation=None, rays=None, mode="rays", panel=None):
        self.plotter = plotter
        self.animation = animation
        self.rays = rays
        self.mode = mode if mode in self.MODES else "rays"
        self._was_playing = True
        panel.heading("表示の切り替え")
        self.label = panel.reserve_text(4)      # 音粒子のときは 4 行になる
        self.apply(render=False)

    def _scalar_bar(self, title):
        try:
            return self.plotter.scalar_bars[title]
        except (KeyError, AttributeError):
            return None

    def _set_visible(self, actor, flag):
        if actor is not None:
            actor.SetVisibility(bool(flag))

    def apply(self, render=True):
        rays = self.mode == "rays"

        if self.rays is not None:
            self.rays.set_visible(rays)
        for title in RAY_BAR_TITLE.values():
            self._set_visible(self._scalar_bar(title), rays)

        if self.animation is not None:
            self._set_visible(self.animation.actor, not rays)
            self._set_visible(self.animation.label, not rays)
            self._set_visible(self._scalar_bar(PARTICLE_BAR_TITLE), not rays)
            slider = getattr(self.animation, "slider", None)
            if slider is not None:
                slider.On() if not rays else slider.Off()
            if rays:
                # 隠している間は止めておく（見えないコマを進めても意味がない）
                self._was_playing = self.animation.playing
                self.animation.playing = False
            else:
                self.animation.playing = self._was_playing

        self._refresh_label()
        if render:
            self.plotter.render()

    def _refresh_label(self):
        if self.mode == "rays":
            text = "いま: 音線\nTab で音粒子へ"
        else:
            text = ("いま: 音粒子\nTab で音線へ\n"
                    "スペース 再生/停止\n← → コマ送り   Home 先頭")
        ParticleAnimation._set_text(self.label, text)

    def toggle(self):
        self.mode = "particles" if self.mode == "rays" else "rays"
        self.apply()


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
            # 描画は plotter.update() に任せる（ここで render すると二重に描くことになる）
            animation.update(animation.step + 1, render=False)
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

    # 動画はパネル無しで（3D だけを大きく写す）
    plotter = vg.build_plotter(model, title="音粒子", off_screen=True,
                               show_normals=False, opacity=opacity,
                               window_size=window_size, panel=False)
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

def view(dxf_path, raylog_path, mode="both", absorption=None, unit=None,
         band_number=None, orient_normals="cad", received_only=False, max_rays=60,
         max_particles=None, pool_size=2000,
         max_reflection=None, colour="energy", band=None, frames=240,
         opacity=0.12, layer_opacity=None, movie=None, point_size=9.0,
         screenshot=None, interval=30):
    """モデルの上に音線と音粒子を重ねて表示する。

    max_rays / max_particles
        **起動時に見せる本数**。実行中は左パネルのスライダで変えられる。
        音線は数十本にしないと線が重なって読めない。
        音粒子は点なので多いほど広がりが分かる（既定は候補すべて）。
    pool_size
        描く候補として抱えておく音線の本数。スライダの上限になる。
    """
    # バンド数は計算時と揃える。既定のままだと「列が足りない」と注意が出て紛らわしい
    kwargs = {} if band_number is None else {"band_number": band_number}
    model = rd.read_model(dxf_path, unit=unit, absorption_table=absorption,
                          orient_normals=orient_normals, **kwargs)
    raylog = RayLog(raylog_path)
    print(f"[view_rays] {raylog.summary()}")

    if mode not in ("rays", "particles", "both"):
        raise ValueError(f"mode は 'rays' / 'particles' / 'both' です: {mode!r}")

    # **描く候補**をまとめて選んでおき、実際に見せる本数はスライダで決める。
    # 音線と音粒子で見やすい本数がまるで違う（線は数十本、点は数千個）ので、
    # 候補は広めに取って表示側で絞る
    pool = raylog.selection(received_only=received_only, max_rays=pool_size)
    if len(pool) == 0:
        raise ValueError("条件に合う音線がありません（--received-only を外してみてください）")
    ray_count = min(max_rays or 60, len(pool))
    particle_count = min(max_particles or len(pool), len(pool))
    print(f"[view_rays] 候補 {len(pool)} 本 / 音線 {ray_count} 本・音粒子 "
          f"{particle_count} 個を表示（スライダで変えられます）")

    if movie is not None:
        path = save_movie(raylog, model, movie,
                          index=even_subset(pool, particle_count), frames=frames,
                          band=band, point_size=point_size, opacity=opacity)
        print(f"[view_rays] 動画を書き出しました: {path}")
        return raylog

    base = os.path.splitext(os.path.basename(dxf_path))[0]
    title = f"{base} 音線・音粒子"
    off_screen = screenshot is not None
    # 静止画（--screenshot）でもパネルごと写す。画面で見えるものと同じにするため
    plotter = vg.build_plotter(model, title=title, off_screen=off_screen,
                               show_normals=False, opacity=opacity,
                               layer_opacity=layer_opacity, show_summary=False,
                               panel=True)
    panel = vg.control_panel(plotter)
    font = vg.japanese_font()

    want_rays = mode in ("rays", "both")
    want_particles = mode in ("particles", "both")

    rays = None
    if want_rays:
        rays = RayDisplay(plotter, raylog, pool, colour=colour, band=band,
                          max_reflection=max_reflection, count=ray_count)

    animation = None
    if want_particles:
        time_label = None
        if panel is not None:
            panel.heading("時刻")
            time_label = panel.reserve_text(2, size=10)
        animation = animate(plotter, raylog, index=pool, frames=frames,
                            band=band, point_size=point_size, label=time_label)
        animation.set_count(particle_count, render=False)
        if off_screen:
            animation.update(frames // 3, render=False)

    if panel is not None:
        vg.add_opacity_control(plotter, font=font, panel=panel, target_key="o")

        panel.heading("表示する本数")
        if rays is not None:
            panel.slider("音線の本数", [1, len(pool)], ray_count,
                         lambda v: rays.rebuild(int(round(v))), fmt="%.0f")
            # 反射をどこまで描くか。0 に寄せると**音源から出た向き**がそのまま見える
            panel.slider("描く反射回数", [1, rays.reflection_limit],
                         rays.max_reflection or rays.reflection_limit,
                         lambda v: rays.set_max_reflection(v), fmt="%.0f")
        if animation is not None:
            panel.slider("音粒子の数", [1, len(pool)], particle_count,
                         lambda v: animation.set_count(int(round(v))), fmt="%.0f")

        switch = None
        if mode == "both":
            # 音線と音粒子を同居させて Tab で切り替える
            switch = RayParticleView(plotter, animation=animation, rays=rays,
                                     mode="rays", panel=panel)
            plotter.add_key_event("Tab", switch.toggle)

        panel.heading("音線の情報")
        panel.text(raylog.summary().replace(" / ", "\n"))
        panel.heading("操作")
        panel.text("z/x/c/v 視点   r リセット   q 終了", color="#7f8794")
        panel.relayout()

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
    p.add_argument("--mode", default="both", choices=["both", "rays", "particles"],
                   help="both=両方を同居させて Tab で切替（既定） / "
                        "rays=音線の折れ線だけ / particles=音粒子だけ")
    p.add_argument("--absorption", help="吸音率 CSV（モデル表示用。省略可）")
    p.add_argument("--unit", help="'mm' / 'm' など")
    p.add_argument("--orient-normals", default="auto",
                   choices=["auto", "cad", "flip", "shells", "inward"],
                   help="法線の扱い。auto=閉じていれば内向き・開いていれば CAD のまま（既定）")
    p.add_argument("--received-only", action="store_true",
                   help="受音した経路だけを描く")
    p.add_argument("--max-rays", type=int, default=60,
                   help="起動時に描く音線の本数（実行中はスライダで変えられる）。"
                        "多すぎると線が重なって読めなくなる。既定 60")
    p.add_argument("--max-particles", type=int, default=None,
                   help="起動時に描く音粒子の数（実行中はスライダで変えられる）。"
                        "既定は候補すべて。点なので多いほど広がりが分かる")
    p.add_argument("--pool", type=int, default=2000,
                   help="描く候補として抱える音線の本数＝スライダの上限。既定 2000")
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
         max_rays=a.max_rays, max_particles=a.max_particles, pool_size=a.pool,
         max_reflection=a.max_reflection, colour=a.color,
         band=a.band, frames=a.frames, opacity=a.opacity,
         layer_opacity=parse_layer_opacity(a.layer_opacity), movie=a.movie,
         point_size=a.point_size, screenshot=a.screenshot, interval=a.interval)


if __name__ == "__main__":
    main()
