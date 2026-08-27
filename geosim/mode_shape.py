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

★2026-08-27 相談で決めた（断面の設定方法）:
  ・**何枚でも一度に切れる**。並びは `断面.json`（`section.py` が読む）
  ・**斜めは「平面図の測線を含む鉛直面」**（傾いた面は入れない）
  ・★**色の基準は「周波数ごとに、全断面で共通」**。断面ごとに最大を取ると
    基準が動いて**枚どうしを比べられない**。周波数をまたいで共通にすると
    低域だけ真っ赤になるので、そこは分ける
  ・音線から出すときは、**全断面の点をまとめて 1 回で追跡する**（F-6 と同じ）
"""
import io
import math
import os

import numpy as np

import project as pj
import section as sec

# 断面の向き（名前 → 座標の添字）。
# ★並びをそのまま添字に使わないこと（`z` を 0 番と取り違えると断面が変わる）
PLANES = ("x", "y", "z")
AXIS = {"x": 0, "y": 1, "z": 2}

# 既定の格子の間隔 [m]
DEFAULT_SPACING = 0.25


def resolve_section(model, plane="z", value=None, section=None):
    """引数を Section にそろえる（`plane`/`value` の昔の呼び方も通す）。"""
    if isinstance(section, sec.Section):
        return section
    if isinstance(plane, sec.Section):
        return plane
    if value is None:
        low = np.asarray(model.extents[0], dtype=float)
        high = np.asarray(model.extents[1], dtype=float)
        value = 0.5 * (low[AXIS[plane]] + high[AXIS[plane]])
    return sec.axis_section(plane, float(value))


def grid_points(model, plane="z", value=None, spacing=DEFAULT_SPACING,
                margin=0.05, section=None):
    """断面の上に格子の点を作る。→ (点 (n,3), 横軸, 縦軸, 形 (nrow, ncol))

    断面は `section`（Section）で渡す。昔ながらに `plane`（`z` など）と
    `value` を渡してもよい（軸に平行な面として組み立てる）。

    ★**面内の 2 本の基底（u, v）に沿って**格子を張るので、斜めの鉛直面でも
      同じ書き方で通る。軸に平行なときは基底が座標軸そのものなので、
      横軸・縦軸は**世界座標のまま**になる（従来と数字が変わらない）。

    `margin` は壁からどれだけ離すか [m]。室の外に出た点は
    `inside_mask` が落とす。
    """
    section = resolve_section(model, plane, value, section)
    low = np.asarray(model.extents[0], dtype=float)
    high = np.asarray(model.extents[1], dtype=float)

    # 室の外接箱の 8 隅を面へ投影して、格子を張る範囲を決める
    corners = np.array([[x, y, z] for x in (low[0], high[0])
                        for y in (low[1], high[1])
                        for z in (low[2], high[2])], dtype=float)
    along, up = section.coordinates(corners)

    def line(values):
        start, stop = float(values.min()) + margin, float(values.max()) - margin
        if stop <= start:
            start, stop = float(values.min()), float(values.max())
        count = max(2, int(round((stop - start) / spacing)) + 1)
        return np.linspace(start, stop, count)

    first, second = line(along), line(up)
    mesh_a, mesh_b = np.meshgrid(first, second)      # (nrow, ncol)
    points = section.place(mesh_a.ravel(), mesh_b.ravel())
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


def section_segments(mesh, plane="z", value=0.0, section=None):
    """三角形と断面の**交線**を集める。→ (m, 2, 2) の線分（断面内の 2 次元座標）

    切った面の**境界**（壁・什器・反射板の切り口）を図に重ねるために使う。
    ★斜めの断面でも同じ（面からの符号付き距離で見るだけ）。
    """
    if isinstance(plane, sec.Section):
        section = plane
    if section is None:
        section = sec.axis_section(plane, float(value))
    segments = []
    for face in mesh:
        vertexes = np.asarray(face.vertexes, dtype=float)
        height = section.height(vertexes)
        if np.all(height > 0) or np.all(height < 0):
            continue                    # 断面と交わらない
        crossings = []
        for k in range(3):
            a, b = k, (k + 1) % 3
            if height[a] == 0.0:
                crossings.append(vertexes[a])
            if height[a] * height[b] < 0.0:
                t = height[a] / (height[a] - height[b])
                crossings.append(vertexes[a] + t * (vertexes[b] - vertexes[a]))
        if len(crossings) >= 2:
            pair = np.asarray(crossings[:2], dtype=float)
            if np.linalg.norm(pair[0] - pair[1]) > 1.0e-9:
                along, up = section.coordinates(pair)
                segments.append([[along[0], up[0]], [along[1], up[1]]])
    return np.asarray(segments) if segments else np.zeros((0, 2, 2))


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
                verbose=True, section=None, model=None):
    """**直方体の虚音源**で断面の音圧分布を出す（全点で同じ並びを使う）。"""
    import atmosphere as at
    import read_dxffile as rd

    frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))
    if model is None:
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              band_number=project.band_number, verbose=False)
    section = resolve_section(model, plane, value, section)
    low, high = np.asarray(model.extents[0]), np.asarray(model.extents[1])
    size = high - low

    source = project.source
    if source is None:
        if not model.source_points:
            raise ValueError("音源がありません（src レイヤか project.json で指定）")
        source = model.source_points[0]
    source = np.asarray(source, dtype=float) - low      # 室の隅を原点にする

    points, first, second, shape = grid_points(model, spacing=spacing,
                                               section=section)
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
        print(f"[音圧分布] 断面『{section.name}』（{section.label()}）/ "
              f"格子 {shape[0]}×{shape[1]}（間隔 {spacing:.2f} m）")

    field = np.empty((len(points), len(frequencies)), dtype=complex)
    for band, frequency in enumerate(frequencies):
        field[:, band] = box_field(size, source, points - low, frequency,
                                   air.sound_velocity, order=order, alpha=alpha)
    # ★**室の外は塗らない**。虚音源は箱の形で並べるが、断面（とくに測線の
    #   鉛直面）は外接箱いっぱいに張るので、実際の室の外へはみ出しうる
    keep = inside_mask(model, points, verbose=False)
    field[~keep, :] = np.nan
    return {"points": points, "field": field, "frequencies": frequencies,
            "outline": section_segments(model.mesh, section=section),
            "shape": shape, "axes": (first, second), "section": section,
            "plane": "xyz"[section.axis] if section.axis is not None else "u",
            "value": float(section.value) if section.value is not None
            else float("nan"),
            "inside": keep,
            "source": np.asarray(project.source if project.source is not None
                                 else model.source_points[0], dtype=float),
            "spacing": spacing, "nref": order, "rays": 0,
            "radius": float("nan"), "path_length": float("nan"),
            "sound_velocity": air.sound_velocity, "model": model,
            "method": "box", "alpha": alpha}


def compute(project, frequencies, plane="z", value=None,
            spacing=DEFAULT_SPACING, rays=None, nref=None, radius=None,
            verbose=True, section=None, model=None):
    """断面 1 枚の音圧分布を求める（音線から）。→ 結果の辞書"""
    return compute_rays(project, frequencies,
                        [resolve_section(model or _model_for(project),
                                         plane, value, section)],
                        spacing=spacing, rays=rays, nref=nref, radius=radius,
                        verbose=verbose, model=model)[0]


def _model_for(project):
    """断面を組み立てるためだけにモデルを読む（軽い読み方）。"""
    import read_dxffile as rd
    return rd.read_model(project.dxf_path, unit=project.unit,
                         band_number=project.band_number, verbose=False)


def compute_rays(project, frequencies, sections, spacing=DEFAULT_SPACING,
                 rays=None, nref=None, radius=None, verbose=True, model=None):
    """断面を**何枚でも**まとめて求める（音線から）。→ 結果の辞書の並び

    ★**音線追跡は 1 回**（全断面の格子点を一度に受音判定する）。
      追跡は受音点に依らないので、断面が増えても追跡の手間は変わらない（F-6）。
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
    sections = [resolve_section(model, section=s_) for s_ in sections]

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

    # ---- 断面ごとに格子を作り、**まとめて 1 本の並び**にする ----
    layout = []
    for section in sections:
        points, first, second, shape = grid_points(model, spacing=spacing,
                                                   section=section)
        keep = inside_mask(model, points, verbose=verbose)
        if verbose:
            print(f"[音圧分布] 断面『{section.name}』（{section.label()}）/ "
                  f"格子 {shape[0]}×{shape[1]}（間隔 {spacing:.2f} m）"
                  f"/ 室の中 {int(keep.sum())} 点")
        layout.append({"section": section, "points": points, "keep": keep,
                       "axes": (first, second), "shape": shape})

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

    targets = np.vstack([item["points"][item["keep"]] for item in layout]) \
        if layout else np.zeros((0, 3))
    histories = lr.loop(source, targets, beam, reflections, model.mesh, sphere,
                        two_sided=project.two_sided,
                        progress=(lambda f: print(f"\r[音圧分布] 音線追跡 "
                                                  f"{f * 100:5.1f}%", end=""))
                        if verbose else None)
    if verbose:
        print()

    for item in layout:
        item["field"] = np.full((len(item["points"]), len(frequencies)),
                                np.nan, dtype=complex)
        item["where"] = np.flatnonzero(item["keep"])

    # 追跡した順に断面へ振り分ける（並べた順と同じ）
    owner, slot = [], []
    for number, item in enumerate(layout):
        owner.extend([number] * len(item["where"]))
        slot.extend(range(len(item["where"])))

    lengths = []
    for index, (point, history) in enumerate(zip(targets, histories)):
        cleaned = ld.delete(history)
        pulses = ln.loop(source, point, cleaned, model.mesh,
                         sound_velocity=air.sound_velocity,
                         band_number=project.band_number,
                         two_sided=project.two_sided, verbose=False)
        if len(pulses.distance):
            lengths.append(float(np.max(pulses.distance)))
        item = layout[owner[index]]
        where = item["where"][slot[index]]
        for band, frequency in enumerate(frequencies):
            item["field"][where, band] = pressure_at(pulses, frequency,
                                                     air.sound_velocity)
        if verbose and (index + 1) % 50 == 0:
            print(f"\r[音圧分布] バックトレース {index + 1}/{len(targets)} 点",
                  end="")
    if verbose:
        print()

    length = float(np.max(lengths)) if lengths else float("nan")
    results = []
    for item in layout:
        section = item["section"]
        results.append({
            "points": item["points"], "field": item["field"],
            "frequencies": frequencies,
            "outline": section_segments(model.mesh, section=section),
            "shape": item["shape"], "axes": item["axes"], "section": section,
            "plane": "xyz"[section.axis] if section.axis is not None else "u",
            "value": float(section.value) if section.value is not None
            else float("nan"),
            "inside": item["keep"], "source": source, "spacing": spacing,
            "nref": reflections, "rays": ray_count, "radius": sphere,
            "path_length": length, "method": "rays",
            "sound_velocity": air.sound_velocity, "model": model})
    return results


# ---- 出力 ------------------------------------------------------------------

def result_name(result, index=None):
    """ファイルに使う短い名前（`01_床上1.2m` のような形）。"""
    section = result.get("section")
    number = result.get("index") if index is None else index
    slug = section.slug() if section is not None else \
        f"{result['plane']}{result['value']:.2f}m"
    return f"{int(number) + 1:02d}_{slug}" if number is not None else slug


def write_csv(project, result, verbose=True, peak=None):
    """`結果/<室>_<条件>_音圧分布_<番号>_<断面名>.csv`（点ごとの音圧レベル）。

    ★平面の決め方（軸に平行か／測線の鉛直面か）は**中身の見出しに書く**。
      斜めが入るとファイル名では表せないので、`view_field` は
      ファイル名ではなくこの行を読む。
    """
    import json

    section = result.get("section")
    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(f"音圧分布_{result_name(result)}.csv"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    level = _levels(result, peak=peak)
    with io.open(path, "w", encoding="utf-8-sig", newline="") as handle:
        if section is not None:
            handle.write("# section: "
                         + json.dumps(section.to_dict(), ensure_ascii=False)
                         + "\n")
        handle.write("# 断面 "
                     + (section.label() if section is not None
                        else f"{result['plane']} = {result['value']:.3f} m")
                     + f" / 格子 "
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


def _levels(result, peak=None):
    """音圧レベル（最大を 0 dB とした相対値）。

    ★`peak` は**周波数ごとの基準**（バンドの数だけの並び）。
      断面が複数あるときは `common_peak()` で全断面をまたいだ最大を渡す。
      断面ごとに最大を取ると基準が動いて**枚どうしを比べられない**。
    """
    magnitude = np.abs(result["field"])
    if peak is None:
        peak = np.nanmax(magnitude) if np.any(np.isfinite(magnitude)) else 1.0
    peak = np.asarray(peak, dtype=float)
    peak = np.where(np.isfinite(peak) & (peak > 0.0), peak, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20.0 * np.log10(magnitude / peak)


def common_peak(results):
    """**周波数ごとに、全断面で共通**の基準（最大の |p|）。→ (バンド数,)"""
    stack = []
    for result in results:
        magnitude = np.abs(result["field"])
        with np.errstate(invalid="ignore"):
            stack.append(np.where(np.isfinite(magnitude), magnitude,
                                  -np.inf).max(axis=0))
    if not stack:
        return None
    peak = np.max(np.vstack(stack), axis=0)
    return np.where(np.isfinite(peak) & (peak > 0.0), peak, 1.0)


def write_figures(project, result, span=30.0, verbose=True, peak=None):
    """断面ごと・周波数ごとの図（音圧レベルを面で塗る）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    import plots as pl

    pl.use_japanese_font()
    level = _levels(result, peak=peak)
    first, second = result["axes"]
    section = result.get("section")
    if section is None:
        section = sec.axis_section(result["plane"], result["value"])
    x_label, y_label = section.axis_labels()

    written = []
    for band, frequency in enumerate(result["frequencies"]):
        values = level[:, band].reshape(result["shape"])
        # ★図の大きさは**断面の形**から決める（縦横比を保つので、
        #   横長の鉛直断面を 8.4×6.4 で描くと上下が空いて題がはみ出す）
        wide = float(first.max() - first.min()) or 1.0
        tall = float(second.max() - second.min()) or 1.0
        height = min(7.0, max(3.2, 5.4 * min(1.0, tall / wide) + 1.6))
        width = min(16.0, max(6.4, (height - 1.6) * wide / tall + 2.8))
        figure, axis = plt.subplots(figsize=(width, height))
        mesh = axis.pcolormesh(first, second, values, shading="nearest",
                               cmap="turbo", vmin=-span, vmax=0.0)
        figure.colorbar(mesh, ax=axis,
                        label="音圧レベル [dB]（全断面の最大を 0 dB）"
                        if peak is not None else "音圧レベル [dB]（最大を 0 dB）")
        # ★切り口（境界面）を重ねる。どこに壁・什器があるかが分かる
        cut = result.get("outline")
        if cut is not None and len(cut):
            from matplotlib.collections import LineCollection
            axis.add_collection(LineCollection(cut, colors="#111111",
                                               linewidths=1.1, alpha=0.85,
                                               zorder=3))

        # 音源は**断面の近く**にあるときだけ描く（面からの距離で見る）
        source = np.asarray(result["source"], dtype=float)
        if abs(float(section.height([source])[0])) < 1.0:
            along, up = section.coordinates([source])
            axis.plot(along[0], up[0], "*",
                      color="#ffffff", markersize=14, markeredgecolor="#000000",
                      markeredgewidth=0.8, label="音源")
            axis.legend(loc="upper right", fontsize=8)
        axis.set_aspect("equal")
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        # 題は 3 行に分ける（断面の名前が長くなるので 1 行に詰めない）
        how = (f"直方体の虚音源 {result['nref']} 次・"
               f"吸音率 {result.get('alpha', 0.0):.2f}"
               if result.get("method") == "box"
               else f"反射 {result['nref']} 回・音線 {result['rays']:,} 本")
        axis.set_title(f"{project.room_label}／{project.condition_label}\n"
                       f"{frequency:.1f} Hz の音圧分布：{section.name}\n"
                       f"{section.label()}／{how}",
                       fontsize=9.5)
        figure.tight_layout()
        path = project.figure_path(
            f"音圧分布_{result_name(result)}_{frequency:.0f}Hz.png")
        figure.savefig(path, dpi=140)
        plt.close(figure)
        written.append(path)
    if verbose and written:
        print(f"[音圧分布] 図を {len(written)} 枚書き出しました"
              f"（例: {written[0]}）")
    return written


def sections_for(project, model=None, sections=None, plane="z", value=None,
                 verbose=True):
    """使う断面を決める。→ Section の並び

    優先は **① 引数で渡された並び → ② `断面.json` → ③ `--plane`/`--at`**。
    ★`--plane`/`--at` を明示したときは、`断面.json` があってもそちらを使う
      （その 1 枚だけ見たいという指示なので）。
    """
    if sections:
        return [s_ if isinstance(s_, sec.Section) else sec.from_dict(s_)
                for s_ in sections]
    if value is None:
        found = sec.load(project, verbose=verbose)
        if found:
            return found
    if model is None:
        model = _model_for(project)
    return [resolve_section(model, plane, value)]


def run(project, frequencies, plane="z", value=None, spacing=DEFAULT_SPACING,
        rays=None, nref=None, radius=None, span=30.0, method="auto",
        order=10, alpha=DEFAULT_ALPHA, verbose=True, sections=None):
    """断面を**何枚でも**計算して CSV と図を書く。→ 結果の並び

    `method` は `box`（直方体の虚音源）/ `rays`（音線から）/ `auto`。
    ★`auto` は**箱なら box**（モード形状がなめらかに出る）、
      箱でなければ `rays`（形は自由だが、点ごとに拾う経路が違うのでまだらになる）。
    """
    model = _model_for(project)
    sections = sections_for(project, model, sections, plane, value,
                            verbose=verbose)
    if verbose:
        print(f"[音圧分布] 断面 {len(sections)} 枚 × 周波数 "
              f"{len(np.atleast_1d(frequencies))} 本")
        print(sec.describe(sections))

    if method == "auto":
        method = "box" if is_box(model) else "rays"
        if verbose:
            print(f"[音圧分布] やり方は自動で **{method}** にしました"
                  + ("（直方体とみなせる）" if method == "box"
                     else "（直方体ではないので音線から）"))
    if method == "box":
        # ★打ち切りの説明（分解能・残る尾）は断面に依らないので**1 回だけ**出す
        results = [compute_box(project, frequencies, section=section,
                               spacing=spacing, order=order, alpha=alpha,
                               verbose=verbose and number == 0, model=model)
                   for number, section in enumerate(sections)]
    else:
        results = compute_rays(project, frequencies, sections,
                               spacing=spacing, rays=rays, nref=nref,
                               radius=radius, verbose=verbose)
    for number, result in enumerate(results):
        result["index"] = number

    if results and verbose and np.isfinite(results[0]["path_length"]):
        length = results[0]["path_length"]
        print(f"[音圧分布] 最長経路 {length:.1f} m"
              f"（分解能の目安 "
              f"{results[0]['sound_velocity'] / length:.2f} Hz）")

    # ★色の基準は**周波数ごとに全断面で共通**（断面ごとだと比べられない）
    peak = common_peak(results) if len(results) > 1 else None
    if verbose and peak is not None:
        print("[音圧分布] 色の基準は**周波数ごとに全断面で共通**にしました"
              "（断面どうしを比べられるように）")
    for result in results:
        write_csv(project, result, verbose=verbose, peak=peak)
        write_figures(project, result, span=span, verbose=verbose, peak=peak)
    return results


def parse_frequencies(values=None, sweep=None):
    """`--frequency` の列挙と `--sweep 20:100:5` の範囲を混ぜる。"""
    found = [float(v) for v in (values or [])]
    for text in (sweep or []):
        parts = str(text).split(":")
        if len(parts) not in (2, 3):
            raise ValueError(f"範囲は 低:高[:刻み] で書いてください（{text}）")
        low, high = float(parts[0]), float(parts[1])
        step = float(parts[2]) if len(parts) == 3 else 1.0
        if step <= 0 or high < low:
            raise ValueError(f"範囲がおかしいです（{text}）")
        count = int(round((high - low) / step)) + 1
        found.extend(float(v) for v in np.linspace(low, low + step * (count - 1),
                                                   count))
    # 同じ周波数は 1 本にまとめる（並びは小さい順）
    return sorted(set(round(v, 6) for v in found))


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="断面の音圧分布（モード形状）を出す")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None, help="条件（シート名）")
    parser.add_argument("--frequency", type=float, action="append",
                        help="周波数 [Hz]（繰り返し指定可）")
    parser.add_argument("--sweep", action="append",
                        help="周波数の範囲 低:高:刻み [Hz]（例 20:100:5）")
    parser.add_argument("--plane", default="z", choices=PLANES,
                        help="断面の向き（既定 z ＝水平断面）")
    parser.add_argument("--at", type=float, default=None,
                        help="断面の位置 [m]。★指定すると 断面.json より優先")
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

    frequencies = parse_frequencies(args.frequency, args.sweep)
    if not frequencies:
        parser.error("周波数を --frequency か --sweep で指定してください")

    base = pj.Project.load(args.folder)
    project = pj.Project(base.folder,
                         **{k: getattr(base, k) for k in pj.DEFAULTS})
    if args.sheet:
        project.condition_sheet = args.sheet
    print(f"[音圧分布] 条件『{project.condition_label}』")
    run(project, frequencies, plane=args.plane, value=args.at,
        spacing=args.spacing, rays=args.rays, nref=args.nref,
        radius=args.radius, span=args.span, method=args.method,
        order=args.order, alpha=args.alpha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
