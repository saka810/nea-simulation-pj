"""プロジェクトの設定で計算を回し、結果と図をプロジェクトフォルダに保存する。

`project.Project` の条件をそのまま `procedure.process()` に渡すだけの薄い層。
GUI（`app.py`）とコマンドラインの両方からここを呼ぶので、
**「どこに何を書くか」の決め方が 1 か所に集まる**ようにしてある。

    cd geosim
    python run_project.py "C:\\Users\\...\\JR"
"""

import os

import numpy as np

import plots
import procedure
import project as pj
from atmosphere import Atmosphere


def run(project, verbose=True, make_figures=True, progress=None):
    """プロジェクトの条件で計算し、結果 CSV と図を書き出す。

    **受音点が複数ある場合は 1 点ずつ回す**（TODO B-9）。
    受音球の位置が変わると音線追跡そのものをやり直す必要があるので、
    1 回の追跡で全受音点を賄うことはできない。
    2 点目以降は `結果/rec2/` `図/rec2/` のように**枝分かれしたフォルダ**に書く。
    """
    project.ensure_dirs()
    project.save()      # 実行した条件を必ず残す（あとで再現できるように）

    dxf = project.dxf_path
    if not dxf or not os.path.exists(dxf):
        raise FileNotFoundError(f"DXF が見つかりません: {project.dxf!r}")

    receivers = _receivers(project)
    if len(receivers) <= 1:
        return _run_one(project, receivers[0] if receivers else None,
                        verbose=verbose, make_figures=make_figures,
                        progress=progress)

    if verbose:
        print(f"[run] 受音点が {len(receivers)} 点あります。1 点ずつ計算します")
    results = []
    for k, point in enumerate(receivers):
        sub = _sub_project(project, k)
        if verbose:
            print(f"\n[run] ── 受音点 {k + 1}/{len(receivers)} "
                  f"{np.round(point, 3).tolist()} → {sub.folder}")
        # ★親（k=0）の project.json には受音点を書き戻さない。
        #   書き戻すと `receiver` が 1 点に固定され、**次回から 1 点目しか回らなくなる**
        results.append(_run_one(sub, point, verbose=verbose,
                                make_figures=make_figures, write_back=(k > 0),
                                progress=_prefixed(progress,
                                                   f"受音点{k + 1}/{len(receivers)} ")))
    return {"receivers": receivers, "results": results, **results[0]}


def _prefixed(progress, prefix):
    """受音点が複数あるとき、どの受音点の処理かを段階名に添える。"""
    if progress is None:
        return None
    return lambda stage, fraction=None: progress(prefix + stage, fraction)


def _receivers(project):
    """計算する受音点の一覧。project.receiver の指定が最優先、無ければ DXF から。"""
    import read_dxffile as rd
    if project.receiver is not None:
        return [np.asarray(project.receiver, dtype=float)]
    probe = rd.read_model(project.dxf_path, unit=project.unit,
                          band_number=project.band_number, verbose=False)
    return [np.asarray(p, dtype=float) for p in probe.receiver_points]


def _sub_project(project, index):
    """受音点 2 点目以降の書き出し先。条件は同じで**フォルダだけ分ける**。"""
    if index == 0:
        return project
    sub = pj.Project.load(project.folder)
    sub.__dict__.update({k: getattr(project, k) for k in pj.DEFAULTS})
    sub.folder = os.path.join(project.folder, f"rec{index + 1}")
    sub.name = f"{project.name} 受音点{index + 1}"
    # DXF と吸音率は親フォルダのものをそのまま使う（絶対パスにしておく）
    sub.dxf = project.dxf_path
    sub.absorption_csv = project.absorption_path or ""
    # ★ここではフォルダを作らない。実際に回す `_run_one` に任せる。
    #   先に作ると、受音点が 1 つしかないのに空の rec2/ が残って紛らわしい
    return sub


def _run_one(project, receiver, verbose=True, make_figures=True,
             write_back=True, progress=None):
    project.ensure_dirs()
    # 前回の結果を消してから回す。条件を変えたときに古いファイルが残っていると、
    # 今回の条件の値だと思って読んでしまう
    project.clear_results(verbose=verbose)
    dxf = project.dxf_path

    # 法線・吸音材の手動指定。面数が合わないときは project 側が警告して空を返す
    flip_faces = _flip_faces_for(project)
    face_materials = _face_materials_for(project)

    if verbose:
        print(f"[run] {project.summary()}")

    results = procedure.process(
        soundsource_point=project.source,
        receiver_point=receiver,
        dxf_filename=dxf,
        sphere_radius=project.radius,
        nref=project.nref,
        soundray_number=project.rays,
        absorption_csv=project.absorption_path,
        absorption_kind=project.absorption_kind,
        layer_assignment=project.assignment,
        band_number=project.band_number,
        unit=project.unit,
        orient_normals=project.orient_normals,
        two_sided=project.two_sided,
        volume=project.volume,
        flip_faces=flip_faces,
        face_materials=face_materials,
        atmosphere=Atmosphere(temperature=project.temperature,
                              humidity=project.humidity,
                              pressure=project.pressure),
        raylog_filename=project.result_path("raylog"),
        raylog_max_rays=project.raylog_max_rays,
        pulse_filename=project.result_path("pulses"),
        impulse_filename=project.result_path("ir"),
        max_time=project.max_time,
        reverberation_filename=project.result_path("rt"),
        decay_filename=project.result_path("decay"),
        statistical_filename=project.result_path("statistical"),
        surface_filename=project.result_path("surface"),
        clarity_filename=project.path(pj.RESULT_DIR, "clarity.csv"),
        statistical=project.statistical,
        progress=progress,
    )

    if make_figures:
        if progress is not None:
            progress("図を書き出し中", None)
        written = plots.save_all(project, results, verbose=verbose)
        if verbose:
            print(f"[run] 図を {len(written)} 枚書き出しました → {project.path(pj.FIGURE_DIR)}")

    # 実際に使った音源・受音点を project.json に残す（DXF から取った場合も分かるように）
    project.source = results["soundsource_point"].tolist()
    if write_back:
        project.receiver = results["receiver_point"].tolist()
    project.save()
    return results


def redraw(project, verbose=True):
    """**計算し直さずに**、保存済みの結果から図を一式描き直す。

    音線追跡（重い）はやり直さない。プロジェクトフォルダに残っている
    `pulses.csv` と `ir.csv` を読み、そこから先だけを計算して `図/` を作り直す。
    研修室（パルス 3901 本）で数秒。

    使いどころ：
    - 図の描き方を直したあと、**過去のプロジェクトに新しい図を反映する**
    - 図だけ消してしまった／新しい図（`mode_buildup.png` など）を後から足す

    ★パルス列とインパルス応答は**そのまま使う**（再合成しない）ので、
      前回の計算結果と食い違うことはない。残響指標・明瞭度・統計残響式は
      本番と同じ関数で計算し直すため、CSV の読み方を別に書かずに済む。
    """
    import read_dxffile as rd
    import reverberation as rv
    import absorption as ab
    import loop_noredundancy as ln

    saved = pj.load_results(project)
    if saved["pulses"] is None:
        raise FileNotFoundError(
            f"{project.result_path('pulses')} がありません。先に計算してください")

    atmosphere = Atmosphere(temperature=project.temperature,
                            humidity=project.humidity,
                            pressure=project.pressure)
    frequencies = ab.octave_bands(project.band_number)

    # ---- パルス列を PulseList に戻す ----
    rows = np.atleast_1d(saved["pulses"])
    names = [n for n in rows.dtype.names if n.startswith("energy_")]
    pulses = ln.PulseList(len(names), atmosphere.sound_velocity)
    pulses.reflection_count = rows["reflection_count"].astype(int)
    pulses.time = rows["time_s"].astype(float)
    pulses.distance = rows["distance_m"].astype(float)
    pulses.direction = np.column_stack([rows["dir_x"], rows["dir_y"], rows["dir_z"]])
    pulses.energy = np.column_stack([rows[n] for n in names])
    if verbose:
        print(f"[redraw] {pulses.summary()}")

    # ---- モデル（外形寸法・容積・レイヤ別面積に要る）----
    # 吸音率の作り方は procedure.process() と同じ手順に揃える
    # （残響室法なら Paris の式で垂直入射へ、レイヤ対応は assignment で差し替え）
    absorption_table = None
    if project.absorption_path:
        library = ab.MaterialLibrary.from_csv(project.absorption_path,
                                              kind=project.absorption_kind)
        absorption_table = library.absorption_table(project.assignment,
                                                    band_number=project.band_number)
    model = rd.read_model(project.dxf_path, band_number=project.band_number,
                          absorption_table=absorption_table, unit=project.unit,
                          orient_normals=project.orient_normals,
                          flip_faces=_flip_faces_for(project),
                          face_materials=_face_materials_for(project), verbose=False)

    results = {"model": model, "pulses": pulses, "frequencies": frequencies,
               "atmosphere": atmosphere, "impulse": None,
               "reverberation": None, "clarity": None, "statistical": None}

    # ---- インパルス応答から先を計算し直す ----
    impulse = saved["ir"]
    if impulse is not None:
        rows = np.atleast_1d(impulse)
        results["impulse"] = (rows["time_s"].astype(float), rows["ir"].astype(float))
        results["reverberation"] = rv.reverberation_time(
            results["impulse"][0], results["impulse"][1], frequencies=frequencies)
        results["clarity"] = rv.clarity_measures(
            results["impulse"][0], results["impulse"][1], frequencies=frequencies)

    # ---- 統計残響式（材料別の面積・吸音率の図に要る）----
    if project.statistical:
        volume = project.volume
        if volume is not None:
            results["statistical"] = rv.statistical_reverberation(
                model.mesh, volume, frequencies=frequencies, atmosphere=atmosphere)
        else:
            results["statistical"] = rv.statistical_reverberation_from_model(
                model, frequencies=frequencies, atmosphere=atmosphere)

    written = plots.save_all(project, results, verbose=verbose)
    if verbose:
        print(f"[redraw] 図を {len(written)} 枚書き出しました "
              f"→ {project.path(pj.FIGURE_DIR)}")
    return written


def _flip_faces_for(project):
    """法線の手動指定を読む。面数の照合のために DXF を軽く 1 回読む。

    受音点ごとの子フォルダには normals.json を置かないので、
    **親フォルダのものを探しに行く**（法線はモデルの性質で、受音点には依らない）。
    """
    owner = _owner_of(project, lambda p: p.load_flipped_faces()[1])
    if owner is None:
        return None
    return owner.flipped_faces_for(_face_count(owner)) or None


def _face_materials_for(project):
    """面ごとの吸音材の割り当てを読む。`_flip_faces_for` と同じ探し方をする
    （どちらもモデルの性質で、受音点には依らない）。"""
    owner = _owner_of(project, lambda p: p.load_face_materials()[1])
    if owner is None:
        return None
    return owner.face_materials_for(_face_count(owner)) or None


def _owner_of(project, load):
    """その指定を持っているプロジェクトを返す。自分に無ければ親フォルダを見る。"""
    if load(project):
        return project
    parent = pj.Project.load(os.path.dirname(project.folder))
    return parent if load(parent) else None


def _face_count(project):
    """面数の照合用に DXF を軽く 1 回読む。"""
    import read_dxffile as rd
    probe = rd.read_model(project.dxf_path, unit=project.unit,
                          band_number=project.band_number, verbose=False)
    return len(probe.mesh)


def main():
    import argparse

    p = argparse.ArgumentParser(description="プロジェクトの条件で計算を回す")
    p.add_argument("folder", help="プロジェクトフォルダ（project.json があるところ）")
    p.add_argument("--no-figures", action="store_true", help="図を書き出さない")
    p.add_argument("--redraw", action="store_true",
                   help="計算し直さず、保存済みの結果から図だけ作り直す")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf:
        raise SystemExit(f"{a.folder} に project.json が無いか、DXF が設定されていません。"
                         f"先に app.py で条件を入力してください")
    if a.redraw:
        redraw(project)
        return
    run(project, make_figures=not a.no_figures)


if __name__ == "__main__":
    main()
