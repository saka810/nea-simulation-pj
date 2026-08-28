# -*- coding: utf-8 -*-
"""**閉じていない辺を確かめる**（どこが開いているのかを見る）。

★2026-08-28 ユーザー要望

> 閉じていない面（辺）を確認したい。なお、そもそも閉じてる必要のないモデルも
> あるので、閉じたモデルを想定する場合に限ります。

決めごと:
  ・開いた辺には**性質の違う 2 種類**が混ざっている（`read_dxffile` 参照）。
    ここでも必ず分けて出す
      - **自由端** … 他の辺に覆われていない。宙に浮いた板の外周・面の抜け。
        **閉じた室のつもりならここが作図ミス**
      - **T 字接合** … 覆われている。壁を帯で分割しただけで**面は閉じている**
  ・★**「閉じているべきか」はモデルの側では決まらない**（一面反射板の検討も
    ふつうにある）。`project.closed_model` で言ってもらい、
    **「閉じている」と言われたときだけ作図ミスとして扱う**
  ・辺は**かたまり（連なり）ごと**にまとめて出す。1 本ずつ並べても場所が分からない
  ・表（CSV）と図（平面＋立面 2 方向）と、3D の画面に重ねる線の 3 つ
"""
import io
import os

import numpy as np

import project as pj
import read_dxffile as rd

# 同じ点とみなす距離 [m]
JOIN_TOLERANCE = 1.0e-6

FREE = "自由端"
TEE = "T字接合"

# 3D に重ねるときの色
FREE_COLOR = "#ff4d4d"
TEE_COLOR = "#ffc14d"


def collect(model):
    """開いた辺を 2 種類に分けて集める。→ {"自由端": [線分…], "T字接合": [線分…]}"""
    triangles = np.array([np.asarray(m.vertexes, dtype=float)
                          for m in model.mesh])
    if not len(triangles):
        return {FREE: [], TEE: []}
    every = rd.open_edge_segments(triangles)
    free = rd.uncovered_open_edges(triangles)

    keys = {_key(a, b) for a, b in free}
    tee = [(a, b) for a, b in every if _key(a, b) not in keys]
    return {FREE: list(free), TEE: tee}


def _key(first, second):
    a = tuple(np.round(first, 6))
    b = tuple(np.round(second, 6))
    return (a, b) if a <= b else (b, a)


def chains(segments, tolerance=JOIN_TOLERANCE):
    """線分をつながりごとにまとめる。→ [点の並び, …]

    1 本ずつ並べても場所が分からないので、**連なりを 1 件**として扱う。
    """
    remaining = [(np.asarray(a, float), np.asarray(b, float))
                 for a, b in segments]
    found = []
    while remaining:
        chain = list(remaining.pop(0))
        moved = True
        while moved:
            moved = False
            for index, (first, second) in enumerate(remaining):
                for end, point in ((-1, chain[-1]), (0, chain[0])):
                    if np.allclose(first, point, atol=tolerance):
                        chain.append(second) if end == -1 \
                            else chain.insert(0, second)
                        remaining.pop(index)
                        moved = True
                        break
                    if np.allclose(second, point, atol=tolerance):
                        chain.append(first) if end == -1 \
                            else chain.insert(0, first)
                        remaining.pop(index)
                        moved = True
                        break
                if moved:
                    break
        found.append(chain)
    return found


def length_of(points):
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def nearby_layers(model, points, limit=0.05):
    """その辺のそばにある面のレイヤ名（どこの部位か分かるように）。"""
    points = np.asarray(points, dtype=float)
    names = {}
    for face in model.mesh:
        vertexes = np.asarray(face.vertexes, dtype=float)
        for point in points:
            if np.min(np.linalg.norm(vertexes - point, axis=1)) <= limit:
                layer = getattr(face, "material", None) or "?"
                names[layer] = names.get(layer, 0) + 1
                break
    return sorted(names, key=lambda k: -names[k])


def summarise(model, verbose=True):
    """まとめ。→ {"自由端": [連なり…], "T字接合": [連なり…], …}"""
    found = collect(model)
    result = {}
    for kind, segments in found.items():
        groups = chains(segments)
        groups.sort(key=length_of, reverse=True)
        result[kind] = groups
    result["segments"] = found
    if verbose:
        for kind in (FREE, TEE):
            total = sum(length_of(g) for g in result[kind])
            print(f"[開いた辺] {kind}: {len(found[kind])} 本 / "
                  f"{len(result[kind])} かたまり / 計 {total:.2f} m")
    return result


# ---- 表 ---------------------------------------------------------------------

def write_csv(project, model, result=None, verbose=True):
    """`結果/<室>_開いた辺.csv`（かたまりごとに 1 行）。"""
    result = result or summarise(model, verbose=False)
    # ★条件にも受音点にも依らない（形だけで決まる）ので `結果/` 直下・条件名なし
    path = project.result_path("open_edges")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with io.open(path, "w", encoding="utf-8-sig", newline="") as handle:
        handle.write("# 自由端＝他の辺に覆われていない開いた辺"
                     "（宙に浮いた板の外周・面の抜け）。"
                     "T字接合＝覆われている（面は閉じている。辺の分け方が違うだけ）\n")
        handle.write("区分,番号,辺の数,長さ_m,近いレイヤ,"
                     "X1_m,Y1_m,Z1_m,X2_m,Y2_m,Z2_m\n")
        for kind in (FREE, TEE):
            for number, chain in enumerate(result[kind], start=1):
                points = np.asarray(chain, dtype=float)
                layers = "・".join(nearby_layers(model, points)[:3]) or "?"
                handle.write(
                    f"{kind},{number},{len(points) - 1},{length_of(points):.3f},"
                    f"{layers},"
                    + ",".join(f"{v:.3f}" for v in points[0])
                    + "," + ",".join(f"{v:.3f}" for v in points[-1]) + "\n")
    if verbose:
        print(f"[開いた辺] 表を書き出しました: {path}")
    return path


# ---- 図 ---------------------------------------------------------------------

VIEWS = (("平面図（上から）", 0, 1, "X", "Y"),
         ("立面図（正面 X-Z）", 0, 2, "X", "Z"),
         ("立面図（側面 Y-Z）", 1, 2, "Y", "Z"))


def write_figure(project, model, result=None, verbose=True):
    """`図/<室>_開いた辺.png`（平面＋立面 2 方向に赤で重ねる）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
    except ImportError:
        return None
    import plots as pl

    pl.use_japanese_font()
    result = result or summarise(model, verbose=False)
    triangles = np.array([np.asarray(m.vertexes, dtype=float)
                          for m in model.mesh])
    normals = np.array([np.asarray(m.normal, dtype=float)
                        for m in model.mesh])
    outline = rd.patch_outline_segments(triangles, normals)

    # ★**図の大きさはモデルの縦横比から決める**。等倍で 3 枚並べるので、
    #   同じ幅にすると平面図だけ縦に伸びて立面 2 枚が潰れる
    low = np.min(triangles.reshape(-1, 3), axis=0)
    high = np.max(triangles.reshape(-1, 3), axis=0)
    span = np.maximum(high - low, 1.0e-6)
    ratios = [span[first] / span[second] for _t, first, second, _x, _y in VIEWS]
    tall = 4.6                                   # 各枚の高さ [inch]
    figure, axes = plt.subplots(
        1, 3, figsize=(min(20.0, tall * sum(ratios) + 2.4), tall + 1.4),
        gridspec_kw={"width_ratios": ratios})
    for axis, (title, first, second, xname, yname) in zip(axes, VIEWS):
        # 室の形は**同一平面パッチの外周**だけ（三角形の辺を全部引くと網目になる）
        if len(outline):
            axis.add_collection(LineCollection(
                outline[:, :, [first, second]], colors="#b8bec9",
                linewidths=0.7, zorder=1))
        for kind, colour, width in ((TEE, TEE_COLOR, 1.4),
                                    (FREE, FREE_COLOR, 2.4)):
            lines = []
            for chain in result[kind]:
                points = np.asarray(chain, dtype=float)
                lines.extend(zip(points[:-1][:, [first, second]],
                                 points[1:][:, [first, second]]))
            if lines:
                axis.add_collection(LineCollection(
                    lines, colors=colour, linewidths=width, zorder=3,
                    label=kind))
        axis.set_title(title, fontsize=10)
        axis.set_xlabel(f"{xname} [m]")
        axis.set_ylabel(f"{yname} [m]")
        axis.set_aspect("equal")
        axis.autoscale_view()
        axis.margins(0.05)
    handles = [plt.Line2D([], [], color=FREE_COLOR, lw=2.4, label="自由端"),
               plt.Line2D([], [], color=TEE_COLOR, lw=1.4, label="T字接合")]
    axes[0].legend(handles=handles, loc="upper right", fontsize=8)
    free_length = sum(length_of(g) for g in result[FREE])
    figure.suptitle(f"{project.room_label}／開いた辺"
                    f"（自由端 {len(result[FREE])} かたまり・{free_length:.1f} m）",
                    fontsize=11)
    figure.tight_layout()
    path = project.figure_path("開いた辺.png", shared=True)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    figure.savefig(path, dpi=140)
    plt.close(figure)
    if verbose:
        print(f"[開いた辺] 図を書き出しました: {path}")
    return path


# ---- 3D の画面に重ねる ------------------------------------------------------

def add_actors(plotter, model, result=None, visible=False):
    """開いた辺を 3D に重ねる。→ {"自由端": actor, "T字接合": actor}

    面の確認画面から呼ぶ。**チェックで出し入れする**。
    """
    import pyvista as pv

    result = result or summarise(model, verbose=False)
    actors = {}
    for kind, colour, width in ((TEE, TEE_COLOR, 3), (FREE, FREE_COLOR, 6)):
        segments = result["segments"][kind]
        if not segments:
            actors[kind] = None
            continue
        points = np.array([p for pair in segments for p in pair], dtype=float)
        cells = np.hstack([[2, 2 * k, 2 * k + 1]
                           for k in range(len(segments))])
        line = pv.PolyData(points, lines=cells)
        actor = plotter.add_mesh(line, color=colour, line_width=width,
                                 lighting=False, pickable=False,
                                 reset_camera=False, name=f"open_{kind}")
        actor.SetVisibility(bool(visible))
        actors[kind] = actor
    return actors


def add_controls(panel, actors, plotter, counts=None):
    """左の欄に「開いた辺を表示」のチェックを足す。→ 操作の一覧に足す行"""
    if not any(actors.values()):
        return []
    panel.heading("開いた辺（閉じているか確かめる）")
    for kind, colour in ((FREE, FREE_COLOR), (TEE, TEE_COLOR)):
        actor = actors.get(kind)
        if actor is None:
            continue
        number = (counts or {}).get(kind)
        label = f"{kind}" + (f"（{number} 本）" if number else "")

        def toggle(flag, actor=actor):
            actor.SetVisibility(bool(flag))
            plotter.render()

        panel.checkbox(label, bool(actor.GetVisibility()), toggle,
                       colour=colour)
    panel.text("自由端＝面の抜け・宙に浮いた板の外周\n"
               "T字接合＝面は閉じている（辺の分け方が違うだけ）", size=8)
    return ["開いた辺はチェックで重ねて表示（赤＝自由端／橙＝T字接合）"]


# ---- 入り口 -----------------------------------------------------------------

def run(project, model=None, verbose=True):
    """表と図を作る。→ 結果"""
    import run_project as rp

    if model is None:
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              absorption_table=rp._absorption_table_for(project),
                              orient_normals=project.orient_normals,
                              band_number=project.band_number,
                              flip_faces=rp._flip_faces_for(project),
                              face_materials=rp._face_materials_for(project),
                              verbose=False)
    result = summarise(model, verbose=verbose)
    if verbose:
        expect = getattr(project, "closed_model", pj.CLOSED_AUTO)
        print("[開いた辺] " + closure_note(result, expect))
    write_csv(project, model, result, verbose=verbose)
    write_figure(project, model, result, verbose=verbose)
    return result


def closure_note(result, expect=None):
    """閉じているかの言い方（★**閉じているべきかはモデルでは決まらない**）。"""
    free = result[FREE]
    total = sum(length_of(g) for g in free)
    if not free:
        return "自由端はありません（面は閉じています）"
    body = (f"自由端が {len(free)} かたまり・計 {total:.1f} m あります")
    if expect == pj.CLOSED_YES:
        return body + "。★閉じた室のつもりなら**作図ミス**です。" \
                      "図と表で場所を確かめてください"
    if expect == pj.CLOSED_NO:
        return body + "（閉じていないモデルとして扱う設定なので、そのままで結構です）"
    return body + "。閉じた室のつもりなら作図ミス、" \
                  "一面反射板などの開いた形なら問題ありません" \
                  "（条件入力の「閉じたモデル」で決められます）"


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="開いた辺（自由端）を確かめる")
    parser.add_argument("folder", help="プロジェクトのフォルダ、または DXF")
    parser.add_argument("--closed", action="store_true",
                        help="閉じた室を想定している（自由端は作図ミス）")
    parser.add_argument("--open", action="store_true",
                        help="閉じていないモデル（反射板だけなど）")
    parser.add_argument("--list", type=int, default=10,
                        help="端末に出すかたまりの数（既定 10）")
    parser.add_argument("--no-write", action="store_true",
                        help="表と図を書かない（端末に出すだけ）")
    args = parser.parse_args(argv)

    if args.folder.lower().endswith(".dxf"):
        # ★プロジェクトを作る前に確かめたいことが多いので、DXF 単体でも通す。
        #   置き場は DXF と同じフォルダ（`結果/` `図/` を作る）
        model = rd.read_model(args.folder, verbose=False)
        result = summarise(model)
        print("[開いた辺] " + closure_note(result, _expect(args)))
        for kind in (FREE, TEE):
            for number, chain in enumerate(result[kind][:args.list], start=1):
                points = np.asarray(chain, dtype=float)
                shape = "輪（閉じた縁）" if np.allclose(points[0], points[-1]) \
                    else "線（端から端）"
                print(f"  {kind} {number}: {length_of(points):6.2f} m  {shape}  "
                      f"{np.round(points[0], 2)} → {np.round(points[-1], 2)}"
                      f"  近いレイヤ "
                      + "・".join(nearby_layers(model, points)[:3]))
        if not args.no_write:
            folder = os.path.dirname(os.path.abspath(args.folder))
            project = pj.Project(folder, **dict(pj.DEFAULTS))
            project.dxf = os.path.basename(args.folder)
            write_csv(project, model, result)
            write_figure(project, model, result)
        return 0

    project = pj.Project.load(args.folder)
    run(project)
    return 0


def _expect(args):
    if args.closed:
        return pj.CLOSED_YES
    if args.open:
        return pj.CLOSED_NO
    return pj.CLOSED_AUTO


if __name__ == "__main__":
    raise SystemExit(main())
