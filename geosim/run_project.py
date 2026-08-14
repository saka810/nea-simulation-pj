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


def run(project, verbose=True, make_figures=True):
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
                        verbose=verbose, make_figures=make_figures)

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
                                make_figures=make_figures, write_back=(k > 0)))
    return {"receivers": receivers, "results": results, **results[0]}


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


def _run_one(project, receiver, verbose=True, make_figures=True, write_back=True):
    project.ensure_dirs()
    # 前回の結果を消してから回す。条件を変えたときに古いファイルが残っていると、
    # 今回の条件の値だと思って読んでしまう
    project.clear_results(verbose=verbose)
    dxf = project.dxf_path

    # 法線の手動指定。面数が合わないときは project 側が警告して空を返す
    flip_faces = _flip_faces_for(project)

    if verbose:
        print(f"[run] {project.summary()}")

    results = procedure.process(
        soundsource_point=project.source,
        reciever_point=receiver,
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
    )

    if make_figures:
        written = plots.save_all(project, results, verbose=verbose)
        if verbose:
            print(f"[run] 図を {len(written)} 枚書き出しました → {project.path(pj.FIGURE_DIR)}")

    # 実際に使った音源・受音点を project.json に残す（DXF から取った場合も分かるように）
    project.source = results["soundsource_point"].tolist()
    if write_back:
        project.receiver = results["reciever_point"].tolist()
    project.save()
    return results


def _flip_faces_for(project):
    """法線の手動指定を読む。面数の照合のために DXF を軽く 1 回読む。

    受音点ごとの子フォルダには normals.json を置かないので、
    **親フォルダのものを探しに行く**（法線はモデルの性質で、受音点には依らない）。
    """
    import read_dxffile as rd
    flipped, data = project.load_flipped_faces()
    if not data:
        parent = pj.Project.load(os.path.dirname(project.folder))
        flipped, data = parent.load_flipped_faces()
        if not data:
            return None
        project = parent
    probe = rd.read_model(project.dxf_path, unit=project.unit,
                          band_number=project.band_number, verbose=False)
    return project.flipped_faces_for(len(probe.mesh)) or None


def main():
    import argparse

    p = argparse.ArgumentParser(description="プロジェクトの条件で計算を回す")
    p.add_argument("folder", help="プロジェクトフォルダ（project.json があるところ）")
    p.add_argument("--no-figures", action="store_true", help="図を書き出さない")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf:
        raise SystemExit(f"{a.folder} に project.json が無いか、DXF が設定されていません。"
                         f"先に app.py で条件を入力してください")
    run(project, make_figures=not a.no_figures)


if __name__ == "__main__":
    main()
