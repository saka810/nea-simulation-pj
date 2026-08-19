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
    **注目する音線の絞り込み**（`ray_filter`。下記の「注目」を参照）
    時刻・操作説明（音線の情報は起動時にコンソールへ出す）

操作:
    共通      ドラッグ 回転 / ホイール 拡大縮小 / `z` `x` `c` `v` 視点 / `r` リセット / `q` 終了
              **`Tab` 音線 ↔ 音粒子の切り替え**（ウィンドウを閉じずに見比べられる）
              `o` 不透明度の対象を切り替え / `m` モデル表示の ON/OFF
    音粒子    `スペース` 再生・一時停止 / `←` `→` 1 コマ送り / `Home` 先頭へ /
              下の横スライダ 時刻を指定
    注目      **`k` 絞り込みの種類**（なし → 近くを通る → この方向に飛ぶ → 1 本だけ）
              **`j` 基準点**（受音点 → 音源 → 拾った点）。**既定は受音点**
              **`0` 絞り込みを解除**
              範囲は左パネルのスライダで決める：
                「近さ [m]」… 基準点からの半径
                **「方位角 [°]」「仰角 [°]」… 見たい方向**（平面方向と縦方向）
                「方向の半角 [°]」… その方向からの許容角
              `p` で 3D 上の点を拾って基準点にすることもできる（**任意**）
    入力      **`t` スライダの値を数字で入力**（`e` は VTK の終了キーなので使えない）
    保存      **`g` いまの画面をそのまま画像で保存**（角度も設定もそのまま）
              **`b` いまの視点で音粒子の動画（GIF）を保存**
              置き場はプロジェクトの `図/画面/`。撮るたびに連番が増える
              （`save_dir` を渡したときだけ。渡さなければキーも出ない）

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
        # ★**スライダで指定された本数（wanted）と、実際に描いた本数（count）を分ける。**
        #   絞り込みで候補が 8 本に減ると count も 8 に下がるが、
        #   wanted を持っていないと**緩めたときに戻らない**（実際に戻らなかった）
        self.wanted = int(count or len(self.pool))
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
        self.rebuild(self.wanted, render=render)

    def rebuild(self, count, render=True):
        """`count` は**スライダで指定された本数**。候補より多ければ候補を全部描く。"""
        self.wanted = max(1, int(count))
        self.count = int(np.clip(self.wanted, 1, len(self.pool)))
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

    def set_pool(self, pool, render=True):
        """**描く候補そのもの**を差し替える（絞り込み。`ray_filter` が作った添字）。

        本数スライダは差し替えたあとの候補に対して効く。
        絞り込んだ結果が本数より少なければ、そのまま全部描く。
        """
        self.pool = np.asarray(pool, dtype=int)
        self.rebuild(self.wanted, render=render)


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
        # 絞り込み（`set_focus`）。None なら絞っていない
        self._focus = None
        self._count = len(self.index)

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

    def set_time_step(self, milliseconds, render=True):
        """**離散化時間**（1 コマあたりの時間）を変える。

        粗くすれば同じ長さを少ないコマ数で見渡せて、細かくすれば
        初期反射のような速い動きが追える。
        全体の長さ（最大時刻）は変えずに、その中の刻み方だけを変える。

        コマ数が変わるので、時刻スライダの範囲も合わせて直す。
        """
        step = max(1e-4, float(milliseconds) / 1000.0)
        frames = int(round(self.raylog.max_time / step)) + 1
        self.frames = max(2, frames)
        self.times = np.arange(self.frames) * step
        self.step = min(self.step, self.frames - 1)

        slider = getattr(self, "slider", None)
        if slider is not None:
            representation = slider.GetRepresentation()
            representation.SetMinimumValue(0)
            representation.SetMaximumValue(self.frames - 1)
            representation.SetValue(self.step)
        self.update(self.step, render=render)

    @property
    def time_step_ms(self):
        """いまの離散化時間 [ms]。"""
        return float(self.times[1] - self.times[0]) * 1000.0 if len(self.times) > 1 \
            else 0.0

    def set_count(self, count, render=True):
        """表示する粒子の数を変える。**等間隔に抜く**ので分布の偏りは出ない。

        絞り込み（`set_focus`）が効いているときは、**その中から**抜く。
        """
        pool = np.nonzero(self._focus)[0] if self._focus is not None             else np.arange(len(self.index))
        keep = even_subset(pool, count)
        self._shown = np.zeros(len(self.index), dtype=bool)
        self._shown[keep] = True
        self._count = int(count)
        self.update(self.step, render=render)

    def set_focus(self, index=None, render=True):
        """**見せる粒子を音線の添字で絞る**（`ray_filter` が作った添字）。

        点の数（＝トポロジ）は変えない。`_shown` を書き換えて
        エネルギーを NaN にし、`nan_opacity=0` で隠すだけ
        （フレームごとにジオメトリを組み直さないための約束。このクラスの説明を参照）。

        `index=None` で解除。
        """
        if index is None:
            self._focus = None
        else:
            wanted = np.zeros(len(self.index), dtype=bool)
            # self.index は「描く候補の音線番号」。その中での位置に直す
            position = {int(v): k for k, v in enumerate(self.index)}
            for value in np.asarray(index, dtype=int):
                k = position.get(int(value))
                if k is not None:
                    wanted[k] = True
            self._focus = wanted
        self.set_count(self._count, render=render)

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


# ------------------------------------------------------------------------------
# 注目したい音線だけを残す（絞り込み）
# ------------------------------------------------------------------------------

# 絞り込みの種類。`k` で順に切り替える
FOCUS_MODES = ("off", "near", "direction", "single")
FOCUS_LABELS = {
    "off": "絞り込みなし",
    "near": "基準点の近くを通る",
    "direction": "この方向に飛ぶ",
    "single": "いちばん近い 1 本だけ",
}

# 基準点の選び方。`j` で順に切り替える。**既定は受音点**
ANCHOR_MODES = ("receiver", "source", "picked")
ANCHOR_LABELS = {"receiver": "受音点", "source": "音源", "picked": "拾った点"}


class RayFocus:
    """**注目したい音線だけを残す**（ユーザー要望 2026-08-19）。

    「この経路をもっと見たい」「この辺りの粒子を見たい」「この方向に飛ぶ音線を見たい」に
    対応する。`k` でやり方を切り替える。

    | 種類 | 何が残るか | 決め方 |
    |---|---|---|
    | 近くを通る | 折れ線が基準点から半径以内を通った音線 | 基準点（`j`）＋「近さ」スライダ |
    | この方向に飛ぶ | **音源から出たときの向き**が指定の向きに近い音線（円錐） | **方位角・仰角スライダ**＋「半角」 |
    | 1 本だけ | 基準点にいちばん近い音線 1 本（反射面の並びも出す） | 基準点（`j`） |

    ★**スライダで決められるようにしてある**（2026-08-19 の指摘）。
      最初はモデル上で `p` を押して点を拾う形だけにしていたが、
      **その操作が難しい**という指摘を受けた。いまは

      - 基準点 … **既定が受音点**。`j` で 受音点 → 音源 → 拾った点 と巡回
      - 方向 … **方位角（平面方向）と仰角（縦方向）のスライダ**

      で決められる。`p`（点を拾う）は**任意**で、使えば基準点になり、
      方向モードならスライダの値もその向きに合わせる（食い違わないように）。

    ★**方向の約束は `project.head_azimuth` と同じ**。
      方位角 0° = +X で反時計回り、仰角 0° が水平・+90° が真上。

    計算は `ray_filter`（画面に触らない純粋な関数）。
    ここは**その結果を `RayDisplay` と `ParticleAnimation` に渡し直すだけ**。

    ★**「近くを通る」は描いている反射回数までで測る**。50 回も反射を追うと
      どの音線も室内を回って基準点の近くを通ってしまい、絞ったつもりが絞れていない
      （研修室で半径 0.5 m に 2000 本中 741 本が該当した。4 回までなら 67 本）。
    """

    def __init__(self, plotter, raylog, pool, rays=None, animation=None,
                 panel=None, radius=0.5, half_angle=15.0, receiver=None,
                 marker_scale=0.02):
        import ray_filter as rfl

        self.plotter = plotter
        self.raylog = raylog
        self.pool = np.asarray(pool, dtype=int)
        self.rays = rays
        self.animation = animation
        self.panel = panel
        self.radius = float(radius)
        self.half_angle = float(half_angle)
        self.mode = "off"
        self.label = None
        self.marker = None
        self.arrow = None
        self.matched = len(self.pool)

        self.source = rfl.source_point(raylog)
        self.receiver = None if receiver is None else np.asarray(receiver, dtype=float)
        # 受音点が分からないモデルでは音源を基準にする
        self.anchor = "receiver" if self.receiver is not None else "source"
        self.picked = None

        # 方向の初期値は「音源から受音点を見る向き」。いちばん見たい向きのはず
        target = self.receiver if self.receiver is not None else None
        if target is not None and not np.allclose(target, self.source):
            self.azimuth, self.elevation = rfl.angles_from_direction(
                target - self.source)
        else:
            self.azimuth, self.elevation = 0.0, 0.0

        lo = raylog.pad_nodes.min(axis=(0, 1))
        hi = raylog.pad_nodes.max(axis=(0, 1))
        self.span = float(np.linalg.norm(hi - lo))
        self.marker_radius = self.span * marker_scale

    # ---- 基準点と方向 ---------------------------------------------------

    def anchor_point(self):
        """いまの基準点。`j` で選んだものを返す（拾っていなければ受音点か音源）。"""
        if self.anchor == "picked" and self.picked is not None:
            return self.picked
        if self.anchor == "receiver" and self.receiver is not None:
            return self.receiver
        return self.source

    def direction(self):
        """いまの方向（方位角・仰角から作る単位ベクトル）。"""
        import ray_filter as rfl
        return rfl.direction_from_angles(self.azimuth, self.elevation)

    # ---- 操作 ----------------------------------------------------------

    def next_mode(self):
        """`k` で種類を順に切り替える。"""
        self.mode = FOCUS_MODES[(FOCUS_MODES.index(self.mode) + 1) % len(FOCUS_MODES)]
        self.apply()

    def next_anchor(self):
        """`j` で基準点を順に切り替える（受音点 → 音源 → 拾った点）。"""
        candidates = [m for m in ANCHOR_MODES
                      if not (m == "receiver" and self.receiver is None)
                      and not (m == "picked" and self.picked is None)]
        here = candidates.index(self.anchor) if self.anchor in candidates else -1
        self.anchor = candidates[(here + 1) % len(candidates)]
        if self.mode == "off":
            self.mode = "near"
        self.apply()

    def set_point(self, point):
        """`p` で点を拾う（**任意**。スライダだけでも使える）。

        方向モードのときは**スライダの値もその向きに合わせる**
        （拾う操作とスライダで状態が食い違わないように）。
        """
        import ray_filter as rfl

        self.picked = np.asarray(point, dtype=float)
        self.anchor = "picked"
        if self.mode == "direction":
            try:
                self.azimuth, self.elevation = rfl.angles_from_direction(
                    self.picked - self.source)
                self._sync_sliders()
            except ValueError:
                pass
        elif self.mode == "off":
            self.mode = "near"
        self.apply()

    def set_radius(self, value, render=True):
        self.radius = float(value)
        if self.mode in ("near", "single"):
            self.apply(render=render)

    def set_half_angle(self, value, render=True):
        self.half_angle = float(value)
        if self.mode == "direction":
            self.apply(render=render)

    def set_azimuth(self, value, render=True):
        """平面方向（真上から見た向き）。0° = +X、反時計回り。"""
        self.azimuth = float(value) % 360.0
        if self.mode == "off":
            self.mode = "direction"
        if self.mode == "direction":
            self.apply(render=render)

    def set_elevation(self, value, render=True):
        """縦方向。0° が水平、+90° が真上。"""
        self.elevation = float(np.clip(value, -90.0, 90.0))
        if self.mode == "off":
            self.mode = "direction"
        if self.mode == "direction":
            self.apply(render=render)

    def reset(self):
        """絞り込みを解除して全部に戻す（`0`）。"""
        self.mode = "off"
        self.apply()

    def _sync_sliders(self):
        """スライダのつまみを内部の値に合わせる（拾って向きが変わったとき）。"""
        wanted = {"方位角 [°]": self.azimuth, "仰角 [°]": self.elevation}
        for control in (self.panel.controls if self.panel is not None else []):
            value = wanted.get(control["label"])
            if value is None:
                continue
            control["widget"].GetRepresentation().SetValue(value)
            control["value"] = value
            # 見出しの数字も直す（つまみだけ動くと、どちらが本当か分からない）
            if "show" in control:
                control["show"](value)

    # ---- 適用 ----------------------------------------------------------

    def _selected(self):
        """いまの条件に合う音線の添字と、パネルに出す補足。"""
        import ray_filter as rfl

        if self.mode == "off":
            return self.pool, ""
        # 描いている範囲に合わせる（このクラスの説明の★を参照）
        limit = self.rays.max_reflection if self.rays is not None else None
        if self.mode == "near":
            index = rfl.near_point(self.raylog, self.anchor_point(), self.radius,
                                   index=self.pool, max_reflection=limit)
            return index, f"半径 {self.radius:.2f} m"
        if self.mode == "direction":
            index = rfl.in_direction(self.raylog, self.direction(),
                                     self.half_angle, index=self.pool)
            return index, f"半角 {self.half_angle:.0f}°"
        one = rfl.nearest_ray(self.raylog, self.anchor_point(), index=self.pool,
                              max_reflection=limit)
        if one is None:
            return np.array([], dtype=int), ""
        return np.array([one], dtype=int), ""

    def apply(self, render=True):
        import ray_filter as rfl

        index, note = self._selected()
        self.matched = len(index)

        if self.matched == 0:
            # 空にすると折れ線が作れず落ちるので、絞り込みは効かせずに知らせる
            print("[view_rays] 条件に合う音線がありません（半径や半角を広げてください）")
            index = self.pool
        if self.rays is not None:
            self.rays.set_pool(index, render=False)
        if self.animation is not None:
            self.animation.set_focus(None if self.mode == "off" else index,
                                     render=False)
        self._draw_guides()
        self._refresh_label(note)
        if self.mode == "single" and self.matched == 1:
            text = rfl.describe_ray(self.raylog, int(index[0]))
            print("[view_rays] " + text.replace("\n", " / "))
        if render:
            self.plotter.render()

    def _draw_guides(self):
        """基準点の球と、方向の矢印を描く。**何を指定しているかが見えないと操作できない。**"""
        for name in ("marker", "arrow"):
            actor = getattr(self, name)
            if actor is not None:
                self.plotter.remove_actor(actor, render=False)
                setattr(self, name, None)
        if self.mode == "off":
            return
        if self.mode in ("near", "single"):
            self.marker = self.plotter.add_mesh(
                pv.Sphere(radius=self.marker_radius, center=self.anchor_point()),
                color="#ffd166", opacity=0.55, lighting=False)
        else:
            # 方向は音源から伸びる矢印で示す（スライダを動かすと向きが変わる）
            self.arrow = self.plotter.add_mesh(
                pv.Arrow(start=self.source, direction=self.direction(),
                         scale=self.span * 0.35, tip_length=0.18,
                         tip_radius=0.035, shaft_radius=0.012),
                color="#ffd166", opacity=0.8, lighting=False)

    def _refresh_label(self, note=""):
        if self.label is None:
            return
        lines = [FOCUS_LABELS[self.mode]]
        if self.mode == "off":
            lines.append(f"候補 {len(self.pool)} 本すべて")
            lines.append("k で種類を選ぶ")
        elif self.mode == "direction":
            lines.append(f"方位 {self.azimuth:.0f}° / 仰角 {self.elevation:+.0f}°"
                         + (f" / {note}" if note else ""))
            lines.append(f"該当 {self.matched} / {len(self.pool)} 本")
        else:
            point = self.anchor_point()
            lines.append(f"基準 {ANCHOR_LABELS[self.anchor]}"
                         f"（{' '.join(f'{v:.1f}' for v in point)}）"
                         + (f" / {note}" if note else ""))
            lines.append(f"該当 {self.matched} / {len(self.pool)} 本")
        vg.set_actor_text(self.label, "\n".join(lines))


def add_focus(plotter, raylog, pool, rays=None, animation=None, panel=None,
              radius=0.5, half_angle=15.0, receiver=None,
              mode_key="k", anchor_key="j", reset_key="0"):
    """絞り込みの部品（スライダとキー）を組み立てる。

    **スライダだけで使える**（`p` で点を拾うのは任意）。
    """
    focus = RayFocus(plotter, raylog, pool, rays=rays, animation=animation,
                     panel=panel, radius=radius, half_angle=half_angle,
                     receiver=receiver)

    if panel is not None:
        panel.heading("注目する音線（k で種類）")
        focus.label = panel.reserve_text(3, size=9)
        panel.slider("近さ [m]", [0.05, 3.0], focus.radius,
                     lambda v: focus.set_radius(v), fmt="%.2f")
        # **方向は平面方向と縦方向のスライダで指定する**（2026-08-19 の指摘）
        panel.slider("方位角 [°]", [0.0, 360.0], focus.azimuth,
                     lambda v: focus.set_azimuth(v), fmt="%.0f")
        panel.slider("仰角 [°]", [-90.0, 90.0], focus.elevation,
                     lambda v: focus.set_elevation(v), fmt="%.0f")
        panel.slider("方向の半角 [°]", [1.0, 90.0], focus.half_angle,
                     lambda v: focus.set_half_angle(v), fmt="%.0f")

    plotter.add_key_event(mode_key, focus.next_mode)
    plotter.add_key_event(anchor_key, focus.next_anchor)
    plotter.add_key_event(reset_key, focus.reset)

    def picked(point, *_):
        # pyvista の版によって引数の数が違うので可変長で受ける
        if point is None:
            return
        focus.set_point(np.asarray(point, dtype=float)[:3])

    try:
        plotter.enable_point_picking(callback=picked, show_message=False,
                                     show_point=False, tolerance=0.02,
                                     left_clicking=False)
    except Exception as e:
        print(f"[view_rays] 点のピックを有効にできませんでした: {type(e).__name__}: {e}")
    focus._refresh_label()
    return focus


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

    vg.finish_window(plotter)
    plotter.show(interactive_update=True, auto_close=False)
    while not closed["flag"]:
        # ★**ウィンドウが無くなったら 1 コマも描かずに抜ける。**
        #   閉じられたあとに描き続けると、GL の文脈が消えた状態で
        #   シェーダを組み直そうとして落ちる（2026-08-19 の segfault）。
        #   interactor まで見るのは、`_closed` が立つのが遅れることがあるため
        if (getattr(plotter, "_closed", False) or plotter.render_window is None
                or plotter.iren is None):
            break
        if animation.playing:
            # 描画は plotter.update() に任せる（ここで render すると二重に描くことになる）
            animation.update(animation.step + 1, render=False)
        try:
            plotter.update(interval)
        except (RuntimeError, AttributeError):
            break
    return animation


def add_movie_key(plotter, animation, folder, stem="音粒子", key="b",
                  max_width=720, max_frames=400, duration=40, colors=128):
    """`b` で、**いまの視点・いまの設定のまま**音粒子の動画（GIF）を保存する。

    `save_movie()` が別に off-screen で作り直すのと違い、こちらは
    **画面に出ているウィンドウをそのまま録る**。回して見つけた角度や、
    絞った粒子数・離散化時間がそのまま動画になる。

    GIF にしているのは追加の依存を増やさないため（Pillow は matplotlib の依存で
    既に入っている。mp4 には ffmpeg が要る）。

    そのぶんの制約:
      ・`max_width` に縮めて、コマごとに 128 色へ落とす（GIF の色数上限は 256）
      ・コマ数は `max_frames` で頭打ちにして等間隔に間引く
        （離散化時間を 0.2 ms にすると 3 秒で 15000 コマになり、
          そのまま溜めるとメモリが持たない）
    """
    def save():
        from PIL import Image

        path = vg.next_free_path(folder, stem, ".gif")
        # `frames` と `times` は `set_time_step()` が揃えているが、
        # 食い違っていても落ちないように短いほうに合わせる
        total = min(int(animation.frames), len(animation.times))
        steps = np.arange(total)
        if total > max_frames:
            steps = np.unique(np.linspace(0, total - 1,
                                          max_frames).round().astype(int))
            print(f"[view_rays] コマが多いので {total} → {len(steps)} に"
                  f"間引きます（動画の長さは変わりません）")

        playing, keep = animation.playing, animation.step
        animation.playing = False
        print(f"[view_rays] 動画を作っています（{len(steps)} コマ）…")

        images = []
        try:
            for step in steps:
                animation.update(int(step))
                image = Image.fromarray(plotter.screenshot(return_img=True))
                if image.width > max_width:
                    height = round(image.height * max_width / image.width)
                    image = image.resize((max_width, height), Image.LANCZOS)
                images.append(image.convert("P", palette=Image.ADAPTIVE,
                                            colors=colors))
            images[0].save(path, save_all=True, append_images=images[1:],
                           duration=duration, loop=0, optimize=True)
        except Exception as e:
            print(f"[view_rays] 動画を保存できませんでした: {type(e).__name__}: {e}")
            return None
        finally:
            animation.update(keep)
            animation.playing = playing

        size = os.path.getsize(path) / 1e6
        print(f"[view_rays] 動画を保存しました: {path}"
              f"（{len(images)} コマ / {size:.1f} MB）")
        return path

    plotter.add_key_event(key, save)
    return save


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
         screenshot=None, interval=30, save_dir=None):
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
                               panel=True, screen="rays")
    panel = vg.control_panel(plotter)
    font = vg.japanese_font()

    want_rays = mode in ("rays", "both")
    want_particles = mode in ("particles", "both")

    focus = None
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
            # 離散化時間（1 コマあたりの時間）。粗くすると全体を見渡せ、
            # 細かくすると初期反射のような速い動きが追える
            panel.slider("離散化時間 [ms]", [0.2, 50.0], animation.time_step_ms,
                         lambda v: animation.set_time_step(v), fmt="%.1f")

        # ---- 注目したい音線だけを残す（p で基準点、k で種類）----
        # 基準点の既定は受音点（クリックしなくても使えるように）
        focus = add_focus(plotter, raylog, pool, rays=rays, animation=animation,
                          panel=panel,
                          receiver=(model.receiver_points[0]
                                    if model.receiver_points else None))

        switch = None
        if mode == "both":
            # 音線と音粒子を同居させて Tab で切り替える
            switch = RayParticleView(plotter, animation=animation, rays=rays,
                                     mode="rays", panel=panel)
            plotter.add_key_event("Tab", switch.toggle)

        # ---- いまの画面をそのまま保存する（G-12）----
        help_lines = ["k 絞り込みの種類   j 基準点",
                      "0 絞り込みを解除   p 点を拾う（任意）",
                      f"{vg.VALUE_INPUT_KEY} 値を数字で入力",
                      "z/x/c/v 視点   r リセット   q 終了"]
        if save_dir:
            # ファイル名は**いま何を見ているか**で変える。音粒子は時刻も入れる
            # （止めた場面を何枚も撮ったとき、あとで見分けられるように）
            def stem():
                showing = switch.mode if switch is not None else \
                    ("particles" if animation is not None else "rays")
                if showing == "particles" and animation is not None:
                    return f"音粒子_{animation.times[animation.step] * 1000:.0f}ms"
                return "音線"

            vg.add_screenshot_key(plotter, save_dir, stem, key="g")
            help_lines.insert(0, "g いまの画面を画像で保存")
            if animation is not None:
                add_movie_key(plotter, animation, save_dir, key="b")
                help_lines.insert(1, "b いまの視点で動画（GIF）を保存")

        # 「音線の情報」は起動時にコンソールへ出しているのでパネルには載せない
        # （操作が増えて入りきらなくなったため。2026-08-19）
        panel.heading("操作")
        panel.text("\n".join(help_lines), color="#7f8794")
        panel.enable_value_input()
        panel.relayout()
        if panel.hidden_height() > 0:
            # 黙って下が切れると操作説明ごと消えるので知らせる。
            # パネルの高さ＝ウィンドウの高さなので、縦に広げれば入る
            print(f"[view_rays] 左パネルが {panel.hidden_height():.0f} px ぶん"
                  f"入りきりません（ウィンドウを縦に広げてください）。操作は以下:")
            for line in help_lines:
                print(f"[view_rays]   {line}")

    plotter.view_isometric()
    if off_screen:
        plotter.screenshot(screenshot)
        plotter.close()
        print(f"[view_rays] 画像を書き出しました: {screenshot}")
        return raylog

    if animation is not None:
        run_animation(plotter, animation, interval=interval)
    else:
        vg.finish_window(plotter)
        plotter.show()
    # **握っている actor をこちらから手放してから閉じる**（`release_window` 参照）。
    # 放っておくと、GL の文脈が消えたあとに解放されて落ちることがある
    vg.release_window(plotter, rays, animation, focus, panel)
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
