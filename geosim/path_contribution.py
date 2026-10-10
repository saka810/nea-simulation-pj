# -*- coding: utf-8 -*-
"""**その周波数で「どの壁・どの方向」が効いているか**を出す（TODO G-43）。

★2026-10-02 ユーザー要望

> 〇Hz の定在波がどの方向からあるいはどの壁由来の伝搬経路が影響しているかを
> 確認できると、その壁を吸音すれば良い、とかがわかるかなと思っています。

決めごと:
  ・**計算はやり直さない**。パルス列（`…_pulses.csv`。条件ごとのエネルギー）と
    経路の幾何（`…_経路.npz`。通った壁の並び）を突き合わせるだけ。
    ★`pulses.csv` には通った壁が入っていない（`from_csv` は空のまま）ので、
    npz と**行を対応づける**（反射回数・経路長・到来方向で並べて照合し、
    1 本でも合わなければ止める）
  ・足し方は伝達関数と同じ**複素和**（査読論文 式(2)）。
    ★エネルギーで分けると干渉が消えて「定在波」の意味が変わる
  ・★★**物差しの本命は「その壁を外したら何 dB 変わるか」**（`除いたときの変化_dB`）
        20 log10 |H − H_壁| / |H|
    H_壁 は**その壁で 1 回でも反射した経路**の複素和。その壁で**反射が起きない
    （どの入射角でも R = 0。開口と同じ）としたときに消える経路そのもの**なので、
    ほかの経路は変わらない ＝ 近似ではなく、その壁を R = 0 にした計算と同じ値になる
    - 負なら「その壁の反射を減らすと、その周波数の山が下がる」
    - ★「垂直入射の吸音率 1」とは違う。このプログラムの反射率は入射角で変わる
      （`sound_ray.energy_decay`）ので、α_n = 1 でも斜めの反射は残る。
      ここで出すのは**効き目の上限**（実際の吸音材ではこれより小さい）
    - ★1 本の経路は何枚も壁を通るので、**壁ごとの値は足しても全体にならない**
      （TODO の注意 3 の (a)）。「外したらどうなるか」は重複があっても意味が崩れない
  ・**分け方は 5 通り**
    - レイヤ … 吸音材を決める単位（その壁を 1 回でも通った経路）
    - 面 … 同一平面パッチ（計算の単位。同じレイヤの向かい合う壁を分けて見たいとき）
    - 壁の組 … 経路が通った面の**組み合わせ**（例：東壁と西壁だけ＝x の軸モード）。
      ★これは**重なりのない分け方**なので、割合を足すと 1 になる
    - 到来方向 … 受音点から見てどの向きから来たか（±X・±Y・±Z のいちばん近い向き）
    - 反射回数 … 何回反射した経路か（G-42 と同じ見方）
  ・`同位相の割合` = Re(H_組 · conj H) / |H|²。重なりのない分け方なら**足すと 1**。
    負なら「その組は全体を打ち消す側にいる」
  ・`--reflections 1:3` で反射回数を絞れる（G-42）。★絞ると全体 H も絞ったものになる
  ・周波数を指定しなければ、（絞った）伝達関数の**山を自動で拾う**
  ・`--width` で周波数に幅を持たせられる（その幅で |·|² を平均）。
    山の位置が分解能（≈ c / 最長経路）の中で揺れるので、1 点より安定する

使い方:
    python path_contribution.py <プロジェクト> --frequency 42.5 --receiver 1
    python path_contribution.py <プロジェクト> --reflections 1:3 --peaks 5
"""
import csv
import math
import os

import numpy as np

import frequency_response as fr
import project as pj

FILE_NAME = "寄与.csv"

# 山を自動で拾うときの既定の数・範囲
DEFAULT_PEAKS = 5

# 図に並べる本数（多すぎると読めない）
FIGURE_ROWS = 12

# 到来方向の呼び名（受音点から見て、その向きから来た）
DIRECTION_NAMES = ("+X 側から", "−X 側から", "+Y 側から", "−Y 側から",
                   "+Z 側（上）から", "−Z 側（下）から")

DIRECT = "直接音"

# 経路の照合で同じとみなす差（pulses.csv は 12 桁で書いてある）
MATCH_DISTANCE = 1.0e-6      # [m]
MATCH_DIRECTION = 1.0e-6


# ---- 経路の幾何を突き合わせる ----------------------------------------------

def _order(reflection_count, distance, direction):
    """照合のための並び（反射回数 → 経路長 → 到来方向）。"""
    direction = np.asarray(direction, dtype=float)
    return np.lexsort((np.round(direction[:, 2], 6), np.round(direction[:, 1], 6),
                       np.round(direction[:, 0], 6), np.round(distance, 6),
                       np.asarray(reflection_count)))


def attach_geometry(pulses, paths):
    """`pulses`（エネルギーあり・壁なし）に `paths`（壁あり）の壁の並びを持たせる。

    どちらも同じ受音点・同じ経路の集合のはずなので、並べ替えて 1 対 1 に照合する。
    ★1 本でも合わなければ `ValueError`（別の計算の npz を読んでいる）。
    → 壁の並びを持った新しい PulseList（並びは `pulses` のまま）
    """
    import copy

    if len(pulses.distance) != len(paths.distance):
        raise ValueError(f"経路の本数が合いません（パルス列 {len(pulses.distance)} 本 / "
                         f"経路の npz {len(paths.distance)} 本）。"
                         f"別の計算の結果が残っているかもしれません")
    a = _order(pulses.reflection_count, pulses.distance, pulses.direction)
    b = _order(paths.reflection_count, paths.distance, paths.direction)
    same = (np.array_equal(np.asarray(pulses.reflection_count)[a],
                           np.asarray(paths.reflection_count)[b])
            and np.all(np.abs(np.asarray(pulses.distance)[a]
                              - np.asarray(paths.distance)[b]) <= MATCH_DISTANCE)
            and np.all(np.abs(np.asarray(pulses.direction)[a]
                              - np.asarray(paths.direction)[b]) <= MATCH_DIRECTION))
    if not same:
        raise ValueError("パルス列と経路の npz の中身が合いません"
                         "（反射回数・経路長・到来方向のどれかが違う）。"
                         "計算し直したあとに古い npz が残っているかもしれません")
    out = copy.copy(pulses)
    walls = np.full((len(pulses.distance), paths.walls.shape[1]), -1, dtype=np.int32)
    walls[a] = paths.walls[b]
    out.walls = walls
    cosines = np.zeros(walls.shape)
    cosines[a] = paths.cos_theta[b]
    out.cos_theta = cosines
    return out


def check_geometry(paths_file, mesh):
    """npz を保存したときのモデルと、いまのモデルが同じ形か。→ 理由（同じなら ""）"""
    import path_cache as pc

    data = np.load(paths_file, allow_pickle=False)
    if "fp_geometry" not in data:
        return ""
    vertices = np.array([np.asarray(m.vertexes, dtype=float) for m in mesh])
    normals = np.array([np.asarray(m.normal, dtype=float) for m in mesh])
    if str(data["fp_geometry"].item()) != pc._digest(vertices, normals):
        return ("モデル（DXF の形か法線の向き）が、経路を保存したときと違います。"
                "壁の番号が対応しないので寄与は出せません（計算し直してください）")
    return ""


# ---- 面とレイヤ ------------------------------------------------------------

class Surfaces:
    """面（同一平面パッチ）とレイヤの呼び名。`walls` の番号（代表の三角形）から引く。"""

    def __init__(self, model, two_sided=False):
        import mesh_method as mm

        mesh = model.mesh
        faces = mm.collision_arrays(mesh, two_sided=two_sided)
        self.patch_of_face = np.asarray(faces.patch_of_face, dtype=np.int64)
        count = int(self.patch_of_face.max()) + 1 if len(self.patch_of_face) else 0
        layers = list(getattr(model, "face_layers", []) or [])
        if len(layers) != len(mesh):
            layers = [str(m.material) for m in mesh]

        area = np.zeros(count)
        centre = np.zeros((count, 3))
        self.layer = [""] * count
        self.normal = np.zeros((count, 3))
        biggest = np.full(count, -1.0)
        for j, face in enumerate(mesh):
            p = int(self.patch_of_face[j])
            v = np.asarray(face.vertexes, dtype=float)
            a = 0.5 * float(np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0])))
            area[p] += a
            centre[p] += a * v.mean(axis=0)
            if a > biggest[p]:            # レイヤはいちばん大きい三角形のもの
                biggest[p] = a
                self.layer[p] = layers[j]
                self.normal[p] = np.asarray(face.normal, dtype=float)
        self.area = area
        self.centre = centre / np.maximum(area, 1e-300)[:, None]

        # 呼び名：レイヤ名＋向き＋（同じレイヤ・同じ向きが複数あれば）番号
        #   short は「壁の組」に並べるときの短い呼び名（`2[y=0.00]`）
        self.name = [""] * count
        self.short = [""] * count
        seen = {}
        for p in np.argsort(-area, kind="stable"):
            key = (self.layer[p], _facing(self.normal[p]))
            seen[key] = seen.get(key, 0) + 1
            number = f" #{seen[key]}" if seen[key] > 1 else ""
            # 向き（法線）だけだと向かい合う壁を取り違えやすいので位置も添える
            axis = int(np.argmax(np.abs(self.normal[p])))
            where = ("" if key[1] == "斜め"
                     else f"・{'xyz'[axis]}={self.centre[p][axis]:.2f} m")
            self.name[p] = f"{self.layer[p]}（法線 {key[1]}{where}）{number}"
            self.short[p] = (f"{self.layer[p]}[斜め]{number}" if key[1] == "斜め"
                             else f"{self.layer[p]}[{'xyz'[axis]}="
                                  f"{self.centre[p][axis]:.2f}]{number}")
        self.count = count

    def patches(self, walls):
        """(n,K) の壁の並び → (n,K) のパッチ番号（-1 はそのまま）。"""
        walls = np.asarray(walls)
        out = np.full(walls.shape, -1, dtype=np.int64)
        hit = walls >= 0
        out[hit] = self.patch_of_face[walls[hit]]
        return out


def _facing(normal):
    """法線のいちばん近い軸の向き（`+X` など。傾いていれば `斜め`）。"""
    normal = np.asarray(normal, dtype=float)
    axis = int(np.argmax(np.abs(normal)))
    if abs(normal[axis]) < math.cos(math.radians(30.0)):
        return "斜め"
    return ("+" if normal[axis] > 0 else "−") + "XYZ"[axis]


def arrival_direction(direction):
    """到来方向 (n,3) → 呼び名の番号 (n,)（DIRECTION_NAMES の並び）。"""
    direction = np.asarray(direction, dtype=float)
    axis = np.argmax(np.abs(direction), axis=1)
    sign = direction[np.arange(len(direction)), axis] < 0
    return axis * 2 + sign.astype(int)


# ---- 寄与 ------------------------------------------------------------------

def path_terms(pulses, frequencies, atmosphere, band=None, air_absorption=True):
    """経路ごとの複素振幅 (F, n)。足せば `frequency_response.response` と一致する。"""
    distance = np.asarray(pulses.distance, dtype=float)
    energy = np.asarray(pulses.energy, dtype=float)
    frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))
    weight = energy.sum(axis=1) if band is None else energy[:, int(band)]
    amplitude = np.sqrt(np.maximum(weight, 0.0) / (4.0 * np.pi)) / distance
    wave = 2.0 * np.pi * frequencies / atmosphere.sound_velocity
    terms = amplitude[None, :] * np.exp(-1j * wave[:, None] * distance[None, :])
    if air_absorption:
        coefficient = np.atleast_1d(atmosphere.absorption_coefficient(frequencies))
        terms *= np.exp(-0.5 * coefficient[:, None] * distance[None, :])
    return terms


def groups(pulses, surfaces):
    """分け方ごとの経路の集まり。→ {区分: [(項目, (n,) bool), …]}

    レイヤと面は**重なりあり**（1 本の経路が何枚も通る）、
    壁の組・到来方向・反射回数は**重なりなし**。
    """
    patches = surfaces.patches(pulses.walls)
    n = len(pulses.distance)
    out = {}

    # 面（パッチ）：その面で 1 回でも反射した経路
    touched = np.zeros((n, surfaces.count), dtype=bool)
    rows, cols = np.nonzero(patches >= 0)
    touched[rows, patches[rows, cols]] = True
    out["面"] = [(surfaces.name[p], touched[:, p])
                for p in range(surfaces.count) if touched[:, p].any()]

    # レイヤ：そのレイヤの面で 1 回でも反射した経路
    layer_names = sorted(set(surfaces.layer))
    out["レイヤ"] = []
    for name in layer_names:
        mask = np.zeros(n, dtype=bool)
        for p in range(surfaces.count):
            if surfaces.layer[p] == name:
                mask |= touched[:, p]
        if mask.any():
            out["レイヤ"].append((name, mask))

    # 壁の組：通った面の組み合わせ（重なりなし）
    keys = {}
    label = np.empty(n, dtype=object)
    for i in range(n):
        key = tuple(sorted(set(int(p) for p in patches[i] if p >= 0)))
        label[i] = key
        keys.setdefault(key, []).append(i)
    out["壁の組"] = []
    for key, members in keys.items():
        mask = np.zeros(n, dtype=bool)
        mask[members] = True
        name = (DIRECT if not key
                else " ＋ ".join(surfaces.short[p] for p in key))
        out["壁の組"].append((name, mask))

    # 到来方向（重なりなし）
    which = arrival_direction(pulses.direction)
    out["到来方向"] = [(DIRECTION_NAMES[k], which == k)
                    for k in range(len(DIRECTION_NAMES)) if np.any(which == k)]

    # 反射回数（重なりなし）
    counts = np.asarray(pulses.reflection_count)
    out["反射回数"] = [(f"{int(c)} 回", counts == c) for c in np.unique(counts)]
    return out


def measures(terms, mask):
    """周波数（幅の中の点）で平均した物差し。→ (除いたときの変化 dB, 大きさ dB, 同位相の割合)"""
    total = terms.sum(axis=1)
    part = terms[:, mask].sum(axis=1)
    power = float(np.mean(np.abs(total) ** 2)) or 1e-300
    removed = float(np.mean(np.abs(total - part) ** 2))
    own = float(np.mean(np.abs(part) ** 2))
    share = float(np.mean(np.real(part * np.conj(total)))) / power
    return (10.0 * math.log10(removed / power + 1e-30),
            10.0 * math.log10(own / power + 1e-30), share)


def pick_peaks(frequencies, magnitude, count=DEFAULT_PEAKS, spacing=0.0):
    """|H| の山（極大）を大きい順に `count` 個。`spacing` [Hz] より近い山は 1 つにする。"""
    magnitude = np.asarray(magnitude, dtype=float)
    if len(magnitude) < 3:
        return []
    inner = np.flatnonzero((magnitude[1:-1] >= magnitude[:-2])
                           & (magnitude[1:-1] > magnitude[2:])) + 1
    chosen = []
    for index in inner[np.argsort(-magnitude[inner], kind="stable")]:
        f = float(frequencies[index])
        if all(abs(f - g) >= spacing for g in chosen):
            chosen.append(f)
        if len(chosen) >= count:
            break
    return sorted(chosen)


def analyse(pulses, surfaces, frequencies, atmosphere, band=None, width=0.0,
            points=9):
    """周波数ごとの寄与。→ {"frequencies": [...], "rows": [(区分, 項目, 経路数, [(Δ, L, 割合), …]), …],
                             "level": [全体のレベル dB, …]}"""
    found = groups(pulses, surfaces)
    rows = {(section, name): (int(mask.sum()), mask, [])
            for section, items in found.items() for name, mask in items}
    level = []
    for frequency in frequencies:
        grid = (np.linspace(frequency - 0.5 * width, frequency + 0.5 * width, points)
                if width > 0 else np.array([frequency]))
        terms = path_terms(pulses, grid, atmosphere, band=band)
        level.append(10.0 * math.log10(
            float(np.mean(np.abs(terms.sum(axis=1)) ** 2)) + 1e-300))
        for key, (_count, mask, values) in rows.items():
            values.append(measures(terms, mask))
    order = ("レイヤ", "面", "壁の組", "到来方向", "反射回数")
    out = []
    for section in order:
        items = [(name, count, values) for (sec, name), (count, _m, values)
                 in rows.items() if sec == section]
        if section in ("レイヤ", "面", "壁の組"):
            # 効いている順（いちばん下がる周波数での「除いたときの変化」が小さい順）
            items.sort(key=lambda item: min(v[0] for v in item[2]))
        out.extend((section, name, count, values) for name, count, values in items)
    return {"frequencies": list(frequencies), "rows": out, "level": level}


# ---- プロジェクトから回す --------------------------------------------------

def _receiver_files(project, index):
    sub = fr._on_shelf(project, index)
    pulses = sub.existing_result_path("pulses")
    paths = sub.paths_cache()
    return sub, pulses, paths


def evaluate(project, frequencies=None, receivers=None, reflections=None,
             band=None, width=0.0, low=fr.DEFAULT_LOW, high=fr.DEFAULT_HIGH,
             step=fr.DEFAULT_STEP, peaks=DEFAULT_PEAKS, verbose=True):
    """受音点ごとの寄与。→ {"recN": 結果, …}"""
    import atmosphere as at
    import loop_noredundancy as ln
    import path_cache as pc
    import run_project as rp

    reflections = fr.parse_reflections(reflections)
    air = at.Atmosphere(temperature=project.temperature,
                        humidity=project.humidity, pressure=project.pressure)
    model = rp._model_for(project)
    surfaces = Surfaces(model, two_sided=project.two_sided)

    wanted = None
    if receivers:
        wanted = {int(str(value).lower().replace("rec", "")) for value in receivers}

    results = {}
    missing = 0
    for index in range(1, 1000):
        if wanted is not None and index not in wanted:
            if index > max(wanted):
                break
            continue
        sub, pulse_file, paths_file = _receiver_files(project, index)
        if pulse_file is None or not os.path.exists(pulse_file):
            missing += 1
            if missing > 3 and (results or wanted is None):
                break
            continue
        missing = 0
        name = f"rec{index}"
        if not os.path.exists(paths_file):
            print(f"[寄与] {name}: 経路の npz（{os.path.basename(paths_file)}）が"
                  f"ありません。通った壁が分からないので飛ばします")
            continue
        reason = check_geometry(paths_file, model.mesh)
        if reason:
            print(f"[寄与] {name}: {reason}")
            continue
        paths = pc.load(paths_file, sound_velocity=air.sound_velocity, verbose=False)
        everything = ln.PulseList.from_csv(pulse_file)
        try:
            everything = attach_geometry(everything, paths)
        except ValueError as error:
            print(f"[寄与] {name}: {error}")
            continue
        pulses = fr.select(everything, reflections)
        if not len(pulses.distance):
            print(f"[寄与] {name}: {fr.reflection_label(reflections)}の経路がありません")
            continue

        chosen = list(frequencies or [])
        if not chosen:
            grid = np.arange(float(low), float(high) + 0.5 * step, float(step))
            magnitude = np.abs(fr.response(pulses, grid, air, band=band))
            chosen = pick_peaks(grid, magnitude, count=peaks,
                                spacing=fr.resolution(pulses, air))
        result = analyse(pulses, surfaces, chosen, air, band=band, width=width)
        result.update({"pulses": len(pulses.distance),
                       "pulses_all": len(everything.distance),
                       "resolution": fr.resolution(pulses, air),
                       "reflections": reflections, "width": width,
                       "project": sub})
        results[name] = result
        if verbose:
            _report(name, result)
    return results


def _report(name, result):
    picked = fr.reflection_label(result["reflections"])
    print(f"[寄与] {name}: 経路 {result['pulses']} 本"
          + (f"（{picked}。全 {result['pulses_all']} 本のうち）" if picked else "")
          + f" / 分解能の目安 {result['resolution']:.2f} Hz")
    for k, frequency in enumerate(result["frequencies"]):
        layers = [(name_, values[k][0]) for section, name_, _c, values
                  in result["rows"] if section == "レイヤ"]
        layers.sort(key=lambda item: item[1])
        top = "、".join(f"{n_} {d:+.1f} dB" for n_, d in layers[:3])
        print(f"[寄与]   {frequency:.1f} Hz：外すと下がるレイヤ  {top}")


# ---- 出力 ------------------------------------------------------------------

def output_name(result):
    return FILE_NAME.replace(".csv", fr.reflection_suffix(result.get("reflections"))
                             + ".csv")


def write_csv(project, results, verbose=True):
    """`結果/[srcM/]recN/<室>_<条件>_寄与.csv`。周波数は**列名**に入れる（表の約束）。"""
    written = []
    for name, result in results.items():
        sub = result["project"]
        path = os.path.join(os.path.dirname(sub.result_path("pulses")),
                            sub.prefixed(output_name(result)))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        picked = fr.reflection_label(result["reflections"])
        with open(path, "w", encoding="utf-8-sig", newline="") as handle:
            handle.write("# 除いたときの変化_dB … その壁（組・向き）を通る経路を"
                         "すべて除いたときのレベルの変化＝その壁で反射が起きない"
                         "（開口と同じ）としたときの変化。負ならその壁の反射を減らすと"
                         "その周波数が下がる（吸音材で得られる効き目の上限）\n")
            handle.write("# 大きさ_dB … その経路だけを足した大きさ（全体を 0 dB）／"
                         "同位相の割合 … Re(H_組・conj H)/|H|²。壁の組・到来方向・"
                         "反射回数は足すと 1\n")
            handle.write("# レイヤ・面は 1 本の経路が何枚も通るので、足しても全体にならない\n")
            if picked:
                handle.write(f"# ★{picked}の経路だけで見ている（全体 H も絞ったもの）\n")
            if result["width"]:
                handle.write(f"# 周波数の幅 ±{result['width'] / 2:.2f} Hz で平均\n")
            writer = csv.writer(handle)
            header = ["区分", "項目", "経路数"]
            for frequency in result["frequencies"]:
                header += [f"除いたときの変化_{frequency:.1f}Hz_dB",
                           f"大きさ_{frequency:.1f}Hz_dB",
                           f"同位相の割合_{frequency:.1f}Hz"]
            writer.writerow(header)
            writer.writerow(["全体", "全体のレベル_dB（相対）", result["pulses"]]
                            + sum([[f"{value:.3f}", "", ""]
                                   for value in result["level"]], []))
            for section, item, count, values in result["rows"]:
                writer.writerow([section, item, count]
                                + sum([[f"{d:.3f}", f"{l:.3f}", f"{s:.4f}"]
                                       for d, l, s in values], []))
        written.append(path)
    if verbose and written:
        print(f"[寄与] 表を {len(written)} 点ぶん書き出しました（例: {written[0]}）")
    return written


def write_figures(project, results, verbose=True):
    """受音点・周波数ごとに 1 枚（レイヤ・面／壁の組・到来方向の「除いたときの変化」）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    import plots as pl

    pl.use_japanese_font()
    written = []
    for name, result in results.items():
        sub = result["project"]
        picked = fr.reflection_label(result["reflections"])
        for k, frequency in enumerate(result["frequencies"]):
            figure, axes = plt.subplots(2, 2, figsize=(15.0, 10.0))
            for axis, section in zip(axes.ravel(), ("レイヤ", "面", "壁の組", "到来方向")):
                items = [(item, values[k][0]) for sec, item, _c, values
                         in result["rows"] if sec == section]
                items.sort(key=lambda item: item[1])
                items = items[:FIGURE_ROWS]
                labels = [_short(item) for item, _v in items][::-1]
                values = [value for _i, value in items][::-1]
                colours = ["#2b6cb0" if value < 0 else "#c05621" for value in values]
                axis.barh(range(len(values)), values, color=colours)
                axis.set_yticks(range(len(values)))
                axis.set_yticklabels(labels, fontsize=8)
                axis.axvline(0.0, color="#333333", linewidth=0.8)
                axis.set_xlabel("除いたときの変化 [dB]（負＝反射を減らすと下がる）",
                                fontsize=8.5)
                axis.set_title(section + ("（重なりあり）" if section in ("レイヤ", "面")
                                          else "（重なりなし）"), fontsize=10)
                axis.grid(axis="x", alpha=0.25)
            figure.suptitle(
                f"{sub.room_label}／{sub.condition_label}　{name}　{frequency:.1f} Hz で"
                f"効いている経路（複素和）"
                + (f"　★{picked}のみ" if picked else "")
                + f"\n経路 {result['pulses']} 本・分解能の目安 {result['resolution']:.2f} Hz"
                + (f"・幅 ±{result['width'] / 2:.2f} Hz で平均" if result["width"] else "")
                + "　※その壁で反射が起きない（開口と同じ）としたときの変化"
                  "＝吸音で得られる効き目の上限", fontsize=10)
            figure.tight_layout()
            path = project.figure_path(
                f"寄与_{name}_{frequency:.1f}Hz"
                f"{fr.reflection_suffix(result['reflections'])}.png")
            figure.savefig(path, dpi=130)
            plt.close(figure)
            written.append(path)
    if verbose and written:
        print(f"[寄与] 図を {len(written)} 枚書き出しました（例: {written[0]}）")
    return written


def _short(text, limit=48):
    text = str(text)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def run(project, frequencies=None, receivers=None, reflections=None, band=None,
        width=0.0, peaks=DEFAULT_PEAKS, low=fr.DEFAULT_LOW, high=fr.DEFAULT_HIGH,
        step=fr.DEFAULT_STEP, verbose=True):
    results = evaluate(project, frequencies=frequencies, receivers=receivers,
                       reflections=reflections, band=band, width=width,
                       low=low, high=high, step=step, peaks=peaks, verbose=verbose)
    if not results:
        print("[寄与] 使える受音点がありません（パルス列と経路の npz が要ります）")
        return results
    write_csv(project, results, verbose=verbose)
    write_figures(project, results, verbose=verbose)
    return results


def main(argv=None):
    import argparse

    import source_mix as sx

    parser = argparse.ArgumentParser(
        description="その周波数でどの壁・どの方向が効いているかを出す（計算はやり直さない）")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None, help="条件（シート名）")
    parser.add_argument("--frequency", type=float, action="append", default=None,
                        help="周波数 [Hz]（繰り返し指定可）。省略すると伝達関数の山を拾う")
    parser.add_argument("--peaks", type=int, default=DEFAULT_PEAKS,
                        help=f"山を拾う数（既定 {DEFAULT_PEAKS}）")
    parser.add_argument("--low", type=float, default=fr.DEFAULT_LOW,
                        help="山を探す範囲の下限 [Hz]")
    parser.add_argument("--high", type=float, default=fr.DEFAULT_HIGH,
                        help="山を探す範囲の上限 [Hz]")
    parser.add_argument("--step", type=float, default=fr.DEFAULT_STEP)
    parser.add_argument("--width", type=float, default=0.0,
                        help="周波数の幅 [Hz]（その幅で平均。既定 0 ＝その周波数だけ）")
    parser.add_argument("--band", type=int, default=None,
                        help="使う帯域の番号（省略すると全帯域の合計）")
    parser.add_argument("--receiver", type=int, action="append", default=None,
                        help="受音点の番号（繰り返し指定可。省略すると全部）")
    parser.add_argument("--reflections", default=None,
                        help="反射回数の範囲 下限:上限（例 1:3）。省略すると全部")
    args = parser.parse_args(argv)

    base = pj.Project.load(args.folder)
    project = pj.Project(base.folder, **{k: getattr(base, k) for k in pj.DEFAULTS})
    if args.sheet:
        project.condition_sheet = args.sheet
    project = sx.default_shelf(project)
    print(f"[寄与] 条件『{project.condition_label}』")
    run(project, frequencies=args.frequency,
        receivers=[f"rec{n}" for n in args.receiver] if args.receiver else None,
        reflections=args.reflections, band=args.band, width=args.width,
        peaks=args.peaks, low=args.low, high=args.high, step=args.step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
