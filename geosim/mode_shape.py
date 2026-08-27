# -*- coding: utf-8 -*-
"""**音圧分布（モード形状）**を断面で見る。

★2026-08-26 ユーザー要望（TODO G-39）

> まずは剛な面として、モード分布を確認しながら、形状を検討します。
> …4（音圧分布）も進めましょう。

指定した**断面**（水平／立面）に受音点を格子で置き、指定周波数の音圧を
**位相ごと重ねて**（式(2) の複素和）求め、面で塗る。
モードの腹・節がそのまま見えるので、形状を変えたときの効き方が分かる。

決めごと:
  ・**音線追跡は 1 回**（格子の全点をまとめて受音判定する。F-6 と同じ考え）
  ・**剛な面なら位相ずれ 0 は厳密**（R = +1）なので、形状検討ではこの絵が使える
  ・★**反射回数が足りないとモードはぼやける**。剛に近い室ほど後部が効くので、
    `nref` を大きく取る（出力に使った値と経路長を書く）
  ・受音球の半径は**格子の間隔と同じ**を既定にする（面を隙間なく覆う）
  ・格子は「点」なので、面の中に入ってしまう点（室の外）は**外して塗らない**
"""
import io
import math
import os

import numpy as np

import project as pj

# 断面の向き（名前 → 座標の添字）。
# ★並びをそのまま添字に使わないこと（`z` を 0 番と取り違えると断面が変わる）
PLANES = ("x", "y", "z")
AXIS = {"x": 0, "y": 1, "z": 2}

# 既定の格子の間隔 [m]
DEFAULT_SPACING = 0.25


def grid_points(model, plane="z", value=None, spacing=DEFAULT_SPACING,
                margin=0.05):
    """断面の上に格子の点を作る。→ (点 (n,3), 横軸, 縦軸, 形 (nrow, ncol))

    `plane` は断面の法線の向き（`z` なら水平断面）、`value` はその座標。
    省略すると室の中央。`margin` は壁からどれだけ離すか [m]。
    """
    low, high = np.asarray(model.extents[0], dtype=float), \
        np.asarray(model.extents[1], dtype=float)
    axis = AXIS[plane]
    others = [k for k in range(3) if k != axis]
    if value is None:
        value = 0.5 * (low[axis] + high[axis])

    def line(index):
        start, stop = low[index] + margin, high[index] - margin
        count = max(2, int(round((stop - start) / spacing)) + 1)
        return np.linspace(start, stop, count)

    first, second = line(others[0]), line(others[1])
    mesh_a, mesh_b = np.meshgrid(first, second)      # (nrow, ncol)
    points = np.zeros((mesh_a.size, 3))
    points[:, others[0]] = mesh_a.ravel()
    points[:, others[1]] = mesh_b.ravel()
    points[:, axis] = float(value)
    return points, first, second, mesh_a.shape


# 室の中かどうかを見るときのレイの本数と、内側とみなすしきい値。
# ★`encloses_point` は「全方向へ飛ばして面に当たった割合」を返す。
#   格子の点ごとに呼ぶので本数は控えめにする（多くても判定は変わらない）
INSIDE_SAMPLES = 24
INSIDE_LEVEL = 0.9

# この点数を超えたら中の判定を省く（時間が読めなくなるため）
INSIDE_LIMIT = 4000


def inside_mask(model, points, samples=INSIDE_SAMPLES, verbose=True):
    """室の中にある点だけ True。**外の点は塗らない**ため。

    閉じていない形（一面だけの反射板など）では判定が効かないので、
    そのときは**全部残す**（黙って落とさない）。
    """
    import read_dxffile as rd

    keep = np.ones(len(points), dtype=bool)
    if len(points) > INSIDE_LIMIT:
        if verbose:
            print(f"[音圧分布] 点が多い（{len(points)}）ので"
                  f"「室の中か」の判定は省きます")
        return keep
    try:
        # `encloses_point` は**三角形の頂点の配列**を受ける（Mesh の並びではない）
        triangles = np.array([np.asarray(m.vertexes, dtype=float)
                              for m in model.mesh])
        for index, point in enumerate(points):
            keep[index] = (rd.encloses_point(triangles, point,
                                             samples=samples) >= INSIDE_LEVEL)
    except Exception as error:
        if verbose:
            print(f"[音圧分布] 室の中かを判定できません"
                  f"（{type(error).__name__}: {error}）。全部の点を使います")
        keep[:] = True
    if not np.any(keep):        # 開いた形なら全部外に見えることがある
        keep[:] = True
    return keep


def pressure_at(pulses, frequency, sound_velocity):
    """1 点の音圧（複素数）。式(2) の複素和。"""
    distance = np.asarray(pulses.distance, dtype=float)
    if not len(distance):
        return 0.0 + 0.0j
    energy = np.asarray(pulses.energy, dtype=float).sum(axis=1)
    amplitude = np.sqrt(np.maximum(energy, 0.0) / (4.0 * np.pi)) / distance
    wave = 2.0 * np.pi * float(frequency) / sound_velocity
    return complex(np.sum(amplitude * np.exp(-1j * wave * distance)))


def is_box(model, tolerance=0.02):
    """室が直方体とみなせるか（容積が寸法の積とほぼ同じか）。"""
    low, high = np.asarray(model.extents[0]), np.asarray(model.extents[1])
    size = high - low
    volume = float(np.prod(size))
    if not volume or model.volume is None:
        return False
    return abs(model.volume - volume) / volume <= tolerance


# 既定の吸音率。★**0（完全に剛）にしない**——減衰が無いと虚音源の級数は
# 収束せず、打ち切る位置で答えが変わってしまう（共鳴が無限に尖るため）。
# 0.02 はコンクリート程度で、モードの形はそのまま見える
DEFAULT_ALPHA = 0.02


def box_images(size, source, order=10, alpha=DEFAULT_ALPHA, distance=None):
    """直方体の虚音源。→ (位置 (n,3), 振幅の係数 (n,))

    `alpha` は 6 面共通の吸音率。振幅は圧力反射率の大きさ `β = √(1−α)` を
    反射回数ぶん掛けたもの。

    ★**室の中心からの距離**で打ち切る（`distance`。既定は `order` × 対角）。
      添字の範囲で切ると鏡像側が 1 殻ぶん足りず、**左右非対称になる**。
    """
    size = np.asarray(size, dtype=float)
    source = np.asarray(source, dtype=float)
    beta = math.sqrt(max(0.0, 1.0 - float(alpha)))
    limit = (float(distance) if distance
             else float(order) * float(np.linalg.norm(size)))

    # 距離で切るぶん、添字は余裕を持って回す
    reach = [int(math.ceil(limit / (2.0 * value))) + 1 for value in size]
    l, m, n = np.meshgrid(np.arange(-reach[0], reach[0] + 1),
                          np.arange(-reach[1], reach[1] + 1),
                          np.arange(-reach[2], reach[2] + 1), indexing="ij")
    l, m, n = l.ravel(), m.ravel(), n.ravel()

    centre = size / 2.0
    positions, weights = [], []
    for p_flag in (0, 1):
        for q_flag in (0, 1):
            for r_flag in (0, 1):
                x = (1 - 2 * p_flag) * source[0] + 2 * l * size[0]
                y = (1 - 2 * q_flag) * source[1] + 2 * m * size[1]
                z = (1 - 2 * r_flag) * source[2] + 2 * n * size[2]
                count = (np.abs(l - p_flag) + np.abs(l)
                         + np.abs(m - q_flag) + np.abs(m)
                         + np.abs(n - r_flag) + np.abs(n))
                place = np.column_stack([x, y, z])
                keep = np.linalg.norm(place - centre, axis=1) <= limit
                positions.append(place[keep])
                weights.append(beta ** count[keep])
    return np.vstack(positions), np.concatenate(weights)


def box_field(size, source, points, frequency, sound_velocity, order=10,
              alpha=DEFAULT_ALPHA, distance=None, chunk=4096):
    """直方体の虚音源から、格子の点の音圧（複素数）を出す。

        p(x) = Σ_i A_i · e^{-j k r_i} / r_i
    """
    images, weights = box_images(size, source, order=order, alpha=alpha,
                                 distance=distance)
    keep = weights > 1.0e-6                 # 効かない虚音源は落とす（速さ）
    images, weights = images[keep], weights[keep]
    wave = 2.0 * np.pi * float(frequency) / sound_velocity

    points = np.asarray(points, dtype=float)
    field = np.zeros(len(points), dtype=complex)
    for start in range(0, len(images), chunk):
        stop = min(start + chunk, len(images))
        delta = points[:, None, :] - images[None, start:stop, :]
        distance = np.sqrt(np.sum(delta ** 2, axis=2))
        distance = np.maximum(distance, 1.0e-6)
        field += np.sum(weights[None, start:stop]
                        * np.exp(-1j * wave * distance) / distance, axis=1)
    return field


def compute_box(project, frequencies, plane="z", value=None,
                spacing=DEFAULT_SPACING, order=10, alpha=DEFAULT_ALPHA,
                verbose=True):
    """**直方体の虚音源**で断面の音圧分布を出す（全点で同じ並びを使う）。"""
    import atmosphere as at
    import read_dxffile as rd
    import run_project as rp

    frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          band_number=project.band_number, verbose=False)
    low, high = np.asarray(model.extents[0]), np.asarray(model.extents[1])
    size = high - low

    source = project.source
    if source is None:
        if not model.source_points:
            raise ValueError("音源がありません（src レイヤか project.json で指定）")
        source = model.source_points[0]
    source = np.asarray(source, dtype=float) - low      # 室の隅を原点にする

    points, first, second, shape = grid_points(model, plane, value, spacing)
    air = at.Atmosphere(temperature=project.temperature,
                        humidity=project.humidity, pressure=project.pressure)
    images, _weights = box_images(size, source, order=order, alpha=alpha)
    if verbose:
        print(f"[音圧分布] 直方体の虚音源で計算します"
              f"（{len(images):,} 個・吸音率 {alpha:.2f}・"
              f"中心から {order * float(np.linalg.norm(size)):.0f} m まで）")
        if alpha <= 0.0:
            print("[音圧分布] ★注意: 完全に剛（吸音率 0）だと虚音源の級数が"
                  "収束しません（打ち切る位置で値が変わります）。"
                  "0.01〜0.05 程度を入れることをお勧めします")
        # ★打ち切りは「応答をどこで切るか」と同じこと。
        #   半径 R ＝ 時間 R/c まで見たことになり、
        #   **分解能 Δf ≈ c/R**、**切り落とした尾の高さ ≈ −60·(R/c)/T [dB]**
        limit = order * float(np.linalg.norm(size))
        window = limit / air.sound_velocity
        volume = float(np.prod(size))
        surface = 2.0 * (size[0] * size[1] + size[1] * size[2]
                         + size[2] * size[0])
        decay = (0.161 * volume / (surface * alpha)) if alpha > 0 else float("inf")
        tail = -60.0 * window / decay if np.isfinite(decay) else 0.0
        print(f"[音圧分布] 打ち切り {limit:.0f} m ＝ 応答 {window * 1000:.0f} ms"
              f"（分解能 {air.sound_velocity / limit:.2f} Hz／"
              f"残る尾 {tail:.0f} dB・残響 {decay:.2f} s の目安）")
        if tail > -20.0:
            print("[音圧分布] ★尾が十分に落ちていません。"
                  "次数（--order）を上げるか、吸音率を上げてください")
        print(f"[音圧分布] 断面 {plane} = {points[0][AXIS[plane]]:.2f} m / "
              f"格子 {shape[0]}×{shape[1]}（間隔 {spacing:.2f} m）")

    field = np.empty((len(points), len(frequencies)), dtype=complex)
    for band, frequency in enumerate(frequencies):
        field[:, band] = box_field(size, source, points - low, frequency,
                                   air.sound_velocity, order=order, alpha=alpha)
    return {"points": points, "field": field, "frequencies": frequencies,
            "shape": shape, "axes": (first, second), "plane": plane,
            "value": float(points[0][AXIS[plane]]),
            "inside": np.ones(len(points), dtype=bool),
            "source": np.asarray(project.source if project.source is not None
                                 else model.source_points[0], dtype=float),
            "spacing": spacing, "nref": order, "rays": 0,
            "radius": float("nan"), "path_length": float("nan"),
            "sound_velocity": air.sound_velocity, "model": model,
            "method": "box", "alpha": alpha}


def compute(project, frequencies, plane="z", value=None,
            spacing=DEFAULT_SPACING, rays=None, nref=None, radius=None,
            verbose=True):
    """断面の音圧分布を求める。→ 結果の辞書

    `frequencies` は 1 つでも並びでもよい（同じ音線追跡を使い回す）。
    """
    import atmosphere as at
    import loop_deleteredundancy as ld
    import loop_noredundancy as ln
    import loop_reflectionmesh as lr
    import run_project as rp
    import sound_ray as sr
    import source_placement as spl

    frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))
    import read_dxffile as rd
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          absorption_table=rp._absorption_table_for(project),
                          orient_normals=project.orient_normals,
                          band_number=project.band_number,
                          flip_faces=rp._flip_faces_for(project),
                          face_materials=rp._face_materials_for(project),
                          verbose=False)

    source = project.source
    if source is None:
        if not model.source_points:
            raise ValueError("音源がありません（src レイヤか project.json で指定）")
        source = model.source_points[0]
    source = np.asarray(source, dtype=float)

    placement = spl.detect(source, model.mesh,
                           tolerance=getattr(project, "source_surface_tolerance", 0.0),
                           direction=getattr(project, "source_direction", None),
                           enabled=bool(getattr(project, "source_on_surface", True)))

    points, first, second, shape = grid_points(model, plane, value, spacing)
    keep = inside_mask(model, points, verbose=verbose)
    if verbose:
        print(f"[音圧分布] 断面 {plane} = "
              f"{points[0][AXIS[plane]]:.2f} m / 格子 {shape[0]}×{shape[1]}"
              f"（間隔 {spacing:.2f} m）/ 室の中 {int(keep.sum())} 点")

    ray_count = int(rays or min(project.rays, 50000))
    reflections = int(nref if nref is not None else project.nref)
    sphere = float(radius if radius is not None else spacing)

    beam = sr.soundray_generator(ray_count)
    if placement.on_surface:
        beam = spl.hemisphere(beam, placement.normal)
        source = placement.point

    air = at.Atmosphere(temperature=project.temperature,
                        humidity=project.humidity, pressure=project.pressure)
    if verbose:
        print(f"[音圧分布] 音線 {ray_count:,} 本 / 最大反射 {reflections} 回 / "
              f"受音球 {sphere:.2f} m で追跡します"
              + ("（面上の音源：半球）" if placement.on_surface else ""))

    targets = points[keep]
    histories = lr.loop(source, targets, beam, reflections, model.mesh, sphere,
                        two_sided=project.two_sided,
                        progress=(lambda f: print(f"\r[音圧分布] 音線追跡 "
                                                  f"{f * 100:5.1f}%", end=""))
                        if verbose else None)
    if verbose:
        print()

    field = np.full((len(points), len(frequencies)), np.nan, dtype=complex)
    lengths = []
    for index, (point, history) in enumerate(zip(targets, histories)):
        cleaned = ld.delete(history)
        pulses = ln.loop(source, point, cleaned, model.mesh,
                         sound_velocity=air.sound_velocity,
                         band_number=project.band_number,
                         two_sided=project.two_sided, verbose=False)
        if len(pulses.distance):
            lengths.append(float(np.max(pulses.distance)))
        where = np.flatnonzero(keep)[index]
        for band, frequency in enumerate(frequencies):
            field[where, band] = pressure_at(pulses, frequency,
                                             air.sound_velocity)
        if verbose and (index + 1) % 50 == 0:
            print(f"\r[音圧分布] バックトレース {index + 1}/{len(targets)} 点",
                  end="")
    if verbose:
        print()

    return {"points": points, "field": field, "frequencies": frequencies,
            "shape": shape, "axes": (first, second), "plane": plane,
            "value": float(points[0][AXIS[plane]]), "inside": keep,
            "source": source, "spacing": spacing, "nref": reflections,
            "rays": ray_count, "radius": sphere,
            "path_length": float(np.max(lengths)) if lengths else float("nan"),
            "sound_velocity": air.sound_velocity, "model": model}


# ---- 出力 ------------------------------------------------------------------

def write_csv(project, result, verbose=True):
    """`結果/<室>_<条件>_音圧分布_<断面>.csv`（格子の点ごとの音圧レベル）。"""
    plane, value = result["plane"], result["value"]
    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(f"音圧分布_{plane}{value:.2f}m.csv"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    level = _levels(result)
    with io.open(path, "w", encoding="utf-8-sig", newline="") as handle:
        handle.write("# 断面 " + f"{plane} = {value:.3f} m / 格子 "
                     f"{result['shape'][0]}×{result['shape'][1]} / "
                     f"間隔 {result['spacing']:.3f} m / 音線 {result['rays']} 本 / "
                     f"最大反射 {result['nref']} 回\n")
        handle.write("X_m,Y_m,Z_m," + ",".join(
            f"{v:.1f}Hz_dB" for v in result["frequencies"]) + "\n")
        for index, point in enumerate(result["points"]):
            cells = ["" if not np.isfinite(level[index, band])
                     else f"{level[index, band]:.3f}"
                     for band in range(level.shape[1])]
            handle.write(f"{point[0]:.3f},{point[1]:.3f},{point[2]:.3f},"
                         + ",".join(cells) + "\n")
    if verbose:
        print(f"[音圧分布] 表を書き出しました: {path}")
    return path


def _levels(result):
    """音圧レベル（最大を 0 dB とした相対値）。"""
    magnitude = np.abs(result["field"])
    peak = np.nanmax(magnitude) or 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20.0 * np.log10(magnitude / peak)


def write_figures(project, result, span=30.0, verbose=True):
    """断面ごと・周波数ごとの図（音圧レベルを面で塗る）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    import plots as pl

    pl.use_japanese_font()
    level = _levels(result)
    first, second = result["axes"]
    axis_index = AXIS[result["plane"]]
    others = [k for k in range(3) if k != axis_index]
    names = ["X", "Y", "Z"]

    written = []
    for band, frequency in enumerate(result["frequencies"]):
        values = level[:, band].reshape(result["shape"])
        figure, axis = plt.subplots(figsize=(8.4, 6.4))
        mesh = axis.pcolormesh(first, second, values, shading="nearest",
                               cmap="turbo", vmin=-span, vmax=0.0)
        figure.colorbar(mesh, ax=axis, label="音圧レベル [dB]（最大を 0 dB）")
        source = result["source"]
        if abs(source[axis_index] - result["value"]) < 1.0:
            axis.plot(source[others[0]], source[others[1]], "*",
                      color="#ffffff", markersize=14, markeredgecolor="#000000",
                      markeredgewidth=0.8, label="音源")
            axis.legend(loc="upper right", fontsize=8)
        axis.set_aspect("equal")
        axis.set_xlabel(f"{names[others[0]]} [m]")
        axis.set_ylabel(f"{names[others[1]]} [m]")
        axis.set_title(f"{project.room_label}／{project.condition_label}\n"
                       f"{frequency:.1f} Hz の音圧分布"
                       f"（{names[axis_index]} = {result['value']:.2f} m 断面）"
                       + (f"／直方体の虚音源 {result['nref']} 次・"
                          f"吸音率 {result.get('alpha', 0.0):.2f}"
                          if result.get("method") == "box"
                          else f"／反射 {result['nref']} 回・"
                               f"音線 {result['rays']:,} 本"),
                       fontsize=10)
        figure.tight_layout()
        path = project.figure_path(
            f"音圧分布_{result['plane']}{result['value']:.2f}m_{frequency:.0f}Hz.png")
        figure.savefig(path, dpi=140)
        plt.close(figure)
        written.append(path)
    if verbose and written:
        print(f"[音圧分布] 図を {len(written)} 枚書き出しました"
              f"（例: {written[0]}）")
    return written


def run(project, frequencies, plane="z", value=None, spacing=DEFAULT_SPACING,
        rays=None, nref=None, radius=None, span=30.0, method="auto",
        order=10, alpha=DEFAULT_ALPHA, verbose=True):
    """`method` は `box`（直方体の虚音源）/ `rays`（音線から）/ `auto`。

    ★`auto` は**箱なら box**（モード形状がなめらかに出る）、
      箱でなければ `rays`（形は自由だが、点ごとに拾う経路が違うのでまだらになる）。
    """
    if method == "auto":
        import read_dxffile as rd
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              band_number=project.band_number, verbose=False)
        method = "box" if is_box(model) else "rays"
        if verbose:
            print(f"[音圧分布] やり方は自動で **{method}** にしました"
                  + ("（直方体とみなせる）" if method == "box"
                     else "（直方体ではないので音線から）"))
    if method == "box":
        result = compute_box(project, frequencies, plane=plane, value=value,
                             spacing=spacing, order=order, alpha=alpha,
                             verbose=verbose)
    else:
        result = compute(project, frequencies, plane=plane, value=value,
                         spacing=spacing, rays=rays, nref=nref, radius=radius,
                         verbose=verbose)
    if verbose and np.isfinite(result["path_length"]):
        print(f"[音圧分布] 最長経路 {result['path_length']:.1f} m"
              f"（分解能の目安 {result['sound_velocity'] / result['path_length']:.2f} Hz）")
    write_csv(project, result, verbose=verbose)
    write_figures(project, result, span=span, verbose=verbose)
    return result


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="断面の音圧分布（モード形状）を出す")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None, help="条件（シート名）")
    parser.add_argument("--frequency", type=float, action="append",
                        required=True, help="周波数 [Hz]（繰り返し指定可）")
    parser.add_argument("--plane", default="z", choices=PLANES,
                        help="断面の向き（既定 z ＝水平断面）")
    parser.add_argument("--at", type=float, default=None,
                        help="断面の位置 [m]（省略すると室の中央）")
    parser.add_argument("--spacing", type=float, default=DEFAULT_SPACING)
    parser.add_argument("--rays", type=int, default=None)
    parser.add_argument("--nref", type=int, default=None)
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--span", type=float, default=30.0,
                        help="色の幅 [dB]（既定 30）")
    parser.add_argument("--method", default="auto",
                        choices=("auto", "box", "rays"),
                        help="box=直方体の虚音源（なめらか）/ rays=音線から / auto")
    parser.add_argument("--order", type=int, default=10,
                        help="虚音源の次数（box のとき。既定 10）")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA,
                        help="6 面共通の吸音率（box のとき。"
                             "既定 0.02。★0 だと級数が収束しない）")
    args = parser.parse_args(argv)

    base = pj.Project.load(args.folder)
    project = pj.Project(base.folder,
                         **{k: getattr(base, k) for k in pj.DEFAULTS})
    if args.sheet:
        project.condition_sheet = args.sheet
    print(f"[音圧分布] 条件『{project.condition_label}』")
    run(project, args.frequency, plane=args.plane, value=args.at,
        spacing=args.spacing, rays=args.rays, nref=args.nref,
        radius=args.radius, span=args.span, method=args.method,
        order=args.order, alpha=args.alpha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
