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


def run(project, verbose=True, make_figures=True, progress=None,
        reuse_paths=True):
    """プロジェクトの条件で計算し、結果 CSV と図を書き出す。

    受音点が複数あるときは、**音線追跡を 1 回で済ませて**受音判定だけ点ごとに行う
    （F-6。追跡は受音点に依らない）。結果は受音点ごとに `結果/recN/`・`図/recN/`、
    受音点に依らないもの（室の吸音と理論値・音線軌跡）は `結果/` 直下へ書く。
    ファイル名の頭には対象室＋条件名（`project.name`）が付く。
    """
    project.ensure_dirs()
    project.save()      # 実行した条件を必ず残す（あとで再現できるように）

    dxf = project.dxf_path
    if not dxf or not os.path.exists(dxf):
        raise FileNotFoundError(f"DXF が見つかりません: {project.dxf!r}")

    receivers = _receivers(project)
    if len(receivers) <= 1:
        # 1 点でも `結果/rec1/` に入れる（点数によって置き場が変わらないように）
        result = _run_one(_sub_project(project, 0),
                          receivers[0] if receivers else None,
                          verbose=verbose, make_figures=make_figures,
                          write_back=False,
                          head_azimuth=project.head_azimuth_for(0),
                          reuse_paths=reuse_paths, progress=progress)
        _write_summaries(project, verbose=verbose)
        return result

    # ★保存した経路が全受音点ぶんそろっていれば、音線追跡そのものを省く（F-9）。
    #   吸音材だけ変えた計算はここで終わり（あとはエネルギーの掛け算だけ）
    if reuse_paths and _paths_ready(project, receivers, verbose=verbose):
        results = []
        for k, point in enumerate(receivers):
            sub = _sub_project(project, k)
            if verbose:
                print("")
                print(f"[run] ── 受音点 {k + 1}/{len(receivers)}"
                      f"（保存した経路から再開）")
            results.append(_run_one(sub, point, verbose=verbose,
                                    make_figures=make_figures, write_back=False,
                                    head_azimuth=project.head_azimuth_for(k),
                                    reuse_paths=True,
                                    progress=_prefixed(progress,
                                                       f"受音点{k + 1}/{len(receivers)} ")))
        _write_summaries(project, verbose=verbose)
        return {"receivers": receivers, "results": results, **results[0]}

    # ★音線追跡は**受音点に依らない**ので 1 回だけ回し、受音点ごとに配る（F-6）。
    #   受音しても音線は打ち切られないため、受音球をいくつ置いても追跡は同じ。
    #   受音点ごとに追い直していたときは、その回数だけ全部やり直していた
    if verbose:
        print(f"[run] 受音点が {len(receivers)} 点あります。"
              f"音線追跡は 1 回で済ませ、受音判定だけ {len(receivers)} 点ぶん行います")
    traced, recorder = _trace_once(project, receivers, verbose=verbose, progress=progress)

    results = []
    for k, point in enumerate(receivers):
        sub = _sub_project(project, k)
        if verbose:
            print(f"\n[run] ── 受音点 {k + 1}/{len(receivers)} "
                  f"{np.round(point, 3).tolist()} → {sub.folder}")
        # ★親（k=0）の project.json には受音点を書き戻さない。
        #   書き戻すと `receiver` が 1 点に固定され、**次回から 1 点目しか回らなくなる**
        # ★project.json は 1 つだけなので受音点は書き戻さない。
        #   書き戻すと `receiver` が 1 点に固定され、次回から 1 点しか回らなくなる
        results.append(_run_one(sub, point, verbose=verbose,
                                make_figures=make_figures, write_back=False,
                                head_azimuth=project.head_azimuth_for(k),
                                traced_history=None if traced is None else traced[k],
                                reuse_paths=False,
                                progress=_prefixed(progress,
                                                   f"受音点{k + 1}/{len(receivers)} ")))
    if recorder is not None:
        # 軌跡は受音点に依らないので `結果/` 直下に 1 つだけ置く。
        # `clear_results` のあとに置かないと消される
        recorder.save_npz(project.result_path("raylog"))
    _write_summaries(project, verbose=verbose)
    return {"receivers": receivers, "results": results, **results[0]}


def _paths_ready(project, receivers, verbose=True):
    """保存した経路が**全受音点ぶん**使えるかを、計算に入る前に確かめる。

    ★ここで確かめてから音線追跡を省く。1 点でも使えなければ
    従来どおり「1 回の追跡を全受音点で共有」する（F-6）ほうが速いので、
    **部分的な使い回しはしない**（点ごとに追跡し直すと共有の利点が消える）。

    指紋（モデルの形・法線・パッチの分け方・音源・受音点・音線数・
    最大反射回数・受音球）が全部合ったときだけ True。
    """
    import mesh_method as mm
    import path_cache as pc

    try:
        for index in range(len(receivers)):
            sub = _sub_project(project, index)
            if not os.path.exists(sub.paths_cache()):
                if verbose:
                    print(f"[run] 受音点 {index + 1} の経路が無いので"
                          f"音線追跡から回します")
                return False
        model = _model_for(project)
        source = _source_of(project, model)
        if source is None:
            return False
        faces = mm.collision_arrays(model.mesh, two_sided=project.two_sided)
        for index, point in enumerate(receivers):
            sub = _sub_project(project, index)
            mark = pc.fingerprint(model.mesh, faces, source, point, project.rays,
                                  project.nref, project.radius, project.two_sided)
            if pc.load(sub.paths_cache(), mark, verbose=False) is None:
                if verbose:
                    # 理由は `pc.compare` が出す。もう一度呼んで表示させる
                    pc.load(sub.paths_cache(), mark, verbose=True)
                return False
    except Exception as error:      # 判定に失敗したら安全側（追跡からやり直す）
        print(f"[run] 経路の使い回しを判定できませんでした: "
              f"{type(error).__name__}: {error}")
        return False
    if verbose:
        print(f"[run] ★保存した経路を使います（{len(receivers)} 点ぶん）。"
              f"音線追跡とバックトレースの幾何は省いて、吸音率だけ当て直します")
    return True


def _model_for(project, verbose=False):
    """プロジェクトの設定で DXF を読む（吸音率・法線・面ごとの材料まで反映）。"""
    import read_dxffile as rd

    table = _absorption_table_for(project, verbose=verbose)
    return rd.read_model(project.dxf_path, unit=project.unit,
                         absorption_table=table,
                         orient_normals=project.orient_normals,
                         band_number=project.band_number,
                         flip_faces=_flip_faces_for(project),
                         face_materials=_face_materials_for(project),
                         verbose=verbose)


def _source_of(project, model):
    if project.source is not None:
        return np.asarray(project.source, dtype=float)
    if model.source_points:
        return np.asarray(model.source_points[0], dtype=float)
    return None


def run_conditions(project, conditions=None, verbose=True, make_figures=True,
                   progress=None):
    """**複数の条件（材料条件表）をまとめて回す**（依頼 2026-08-21）。

    > 複数条件やる場合、一括で回せると嬉しいです。

    経路（反射面の並びと入射角）は吸音に依らないので、**1 つ目の条件で
    音線追跡まで済ませれば、2 つ目以降はエネルギーの掛け算だけ**で終わる（F-9）。
    実測（研修室・受音点 5 点）で 1 条件目 10 分 → 2 条件目以降 数十秒。

    引数:
        conditions : 条件のリスト。`(条件表のパス, シート名)` の組か、
                     パスだけ（その中の条件シートに展開する）。
                     None ならプロジェクトフォルダの条件を全部
                     （`condition_table.discover`。**xlsx はシートごとに 1 条件**）

    結果は条件ごとに別のファイル名で並ぶ（頭が「対象室名_条件名」になる）。
    最後に**条件を横に並べた比較表**を作る（`summary.write_condition_summary`）。
    """
    import condition_table as ct
    import summary as sm

    if conditions is None:
        conditions = ct.discover(project.folder, verbose=verbose)
    if not conditions:
        if verbose:
            print(f"[run] 条件表が見つかりません。1 条件として回します")
        return {"conditions": [], "results": [run(project, verbose=verbose,
                                                  make_figures=make_figures,
                                                  progress=progress)]}

    results, done = [], []
    for i, (file_name, sheet) in enumerate(conditions):
        sub = pj.Project(project.folder,
                         **{k: getattr(project, k) for k in pj.DEFAULTS})
        sub.condition_csv = file_name
        sub.condition_sheet = sheet or ""
        if verbose:
            print("")
            print("=" * 70)
            print(f"[run] 条件 {i + 1}/{len(conditions)}: "
                  f"{ct.label_of(file_name, sheet)} → 結果の頭 "
                  f"{sub.file_prefix!r}")
            print("=" * 70)
        stage = _prefixed(progress, f"条件{i + 1}/{len(conditions)} ")
        results.append(run(sub, verbose=verbose, make_figures=make_figures,
                           progress=stage))
        done.append((file_name, sheet))

    # 条件を横に並べた比較表。**全条件が終わってから**でないと作れない
    comparison = None
    try:
        comparison = sm.write_condition_summary(project, done, verbose=verbose)
    except Exception as error:
        print(f"[run] 条件の比較表を作れませんでした: "
              f"{type(error).__name__}: {error}")

    # 比較表ができたので、条件ごとの Excel を作り直して比較シートを入れる
    # （条件ごとの Excel は計算の途中で書いているので、まだ比較表が無かった）
    if comparison is not None:
        try:
            import workbook as wb
            for file_name, sheet in done:
                sub = pj.Project(project.folder,
                                 **{k: getattr(project, k) for k in pj.DEFAULTS})
                sub.condition_csv = file_name
                sub.condition_sheet = sheet or ""
                wb.write(sub, verbose=False)
            if verbose:
                print(f"[run] 条件ごとの Excel に比較シートを入れました"
                      f"（{len(done)} 件）")
        except Exception as error:
            print(f"[run] 結果一式（Excel）を作れませんでした: "
                  f"{type(error).__name__}: {error}")
    return {"conditions": done, "results": results, "comparison": comparison}


def _write_summaries(project, verbose=True):
    """受音点をまたいだまとめ表を作る（`結果/まとめ_*.csv`）。

    全測定点を 1 つのファイルで見たいという要望（2026-08-21）。
    残響時間には理論値（統計残響式）も同じ表に入れる。
    まとめだけ作り直したいときは `python summary.py <プロジェクト>`。
    """
    try:
        import summary as sm
        sm.write_all(project, verbose=verbose)
    except Exception as error:      # まとめが作れなくても本体の結果は残す
        print(f"[run] まとめ表を作れませんでした: {type(error).__name__}: {error}")
    # 体裁を整える用の Excel（CSV とは役割を分ける。`workbook.py` 冒頭参照）
    try:
        import workbook as wb
        wb.write(project, verbose=verbose)
    except Exception as error:      # Excel が作れなくても CSV は残る
        print(f"[run] 結果一式（Excel）を作れませんでした: "
              f"{type(error).__name__}: {error}")


def _trace_once(project, receivers, verbose=True, progress=None):
    """**全受音点ぶんの音線追跡を 1 回で**行い、受音点ごとの反射面 ID 履歴を返す。

    追跡そのものは受音点に依らない（受音しても音線は打ち切られない）。
    受音球だけ受音点の数だけ置いて判定すればよいので、ここで 1 回に畳む。
    可視化用の軌跡は 1 本ぶんしか作らないため、**「どれかの受音点に届いた」**が
    受音の印になる（受音点ごとの色分けが要るようになったら作り直す）。

    読み込みに失敗するなど何かあれば None を返し、呼び出し側は
    従来どおり受音点ごとに追跡する（安全側）。
    """
    import loop_reflectionmesh as lr
    import read_dxffile as rd
    import sound_ray as sr
    from ray_recorder import RayRecorder
    import absorption as ab

    try:
        table = _absorption_table_for(project)
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              absorption_table=table,
                              orient_normals=project.orient_normals,
                              band_number=project.band_number,
                              flip_faces=_flip_faces_for(project),
                              face_materials=_face_materials_for(project),
                              verbose=False)
        source = (project.source if project.source is not None
                  else (model.source_points[0] if model.source_points else None))
        if source is None:
            return None
        rays = sr.soundray_generator(project.rays)
        recorder = RayRecorder(total_rays=project.rays,
                               max_rays=project.raylog_max_rays,
                               sound_velocity=Atmosphere(
                                   temperature=project.temperature,
                                   humidity=project.humidity,
                                   pressure=project.pressure).sound_velocity,
                               band_number=project.band_number)
        if verbose:
            print("[run] 音線追跡（全受音点ぶんを 1 回で）")
        histories = lr.loop(np.asarray(source, dtype=float),
                            np.asarray(receivers, dtype=float), rays,
                            project.nref, model.mesh, project.radius,
                            recorder=recorder, two_sided=project.two_sided,
                            progress=(lambda f: progress("音線追跡（共有）", f))
                            if progress else None)
        print("音線軌跡:", recorder.summary())
        # ★軌跡の保存は**呼び出し側が `_run_one` のあとで**行う。
        #   `_run_one` の頭で `clear_results()` が走るので、先に置くと消される
        return histories, recorder
    except Exception as error:      # 何かあっても従来どおり受音点ごとに追跡できる
        print(f"[run] 音線追跡の共有に失敗したので受音点ごとに追跡します: "
              f"{type(error).__name__}: {error}")
        return None, None


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
    """受音点 `index` 番目（0 始まり）を扱う `Project` を返す。

    ★**フォルダは分けない。**`receiver_index` を立てるだけで、結果は
    `結果/recN/`・図は `図/recN/` に入る（`Project.result_dir` が振り分ける）。
    以前は 2 点目以降だけ `rec2/` という別フォルダを作っていたので、
    1 点目だけ `結果/` 直下という不揃いな置き方になっていた（2026-08-21 に直した）。

    顔の向きだけは受音点ごとに違うので、その点のぶんを取り出して入れ直す
    （`head_azimuth` は数値でもリストでもよい。`Project.head_azimuth_for` を参照）。
    """
    sub = pj.Project(project.folder,
                     **{k: getattr(project, k) for k in pj.DEFAULTS})
    sub.head_azimuth = project.head_azimuth_for(index)
    sub.receiver_index = index + 1
    # ★**名前は変えない。**`name` は結果ファイル名の頭に付く（対象室＋条件名）ので、
    #   受音点ごとに変えるとファイル名が受音点ごとに違ってしまう。
    #   何点目かは `receiver_index` が持っていて `Project.summary()` が表示する
    return sub


def _run_one(project, receiver, verbose=True, make_figures=True,
             write_back=True, head_azimuth=None, traced_history=None,
             reuse_paths=True, progress=None):
    project.ensure_dirs()
    # 前回の結果を消してから回す。条件を変えたときに古いファイルが残っていると、
    # 今回の条件の値だと思って読んでしまう。
    # ★経路を使い回すときは音線軌跡も作り直さないので消さない（`keep`）
    project.clear_results(verbose=verbose,
                          keep=("raylog",) if reuse_paths else ())
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
        material_library=_library_for(project),
        layer_assignment=_assignment_for(project),
        band_number=project.band_number,
        unit=project.unit,
        orient_normals=project.orient_normals,
        two_sided=project.two_sided,
        volume=project.volume,
        flip_faces=flip_faces,
        face_materials=face_materials,
        traced_history=traced_history,
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
        room_filename=project.result_path("room"),
        clarity_filename=project.clarity_path(),
        level_filename=project.result_path("spl"),
        sti_filename=project.result_path("sti"),
        paths_filename=project.paths_cache(),
        reuse_paths=reuse_paths,
        source_power_db=project.source_power_db,
        noise_level_db=project.noise_level_db,
        statistical=project.statistical,
        progress=progress,
    )

    # 実際に使った条件を材料条件表に書き戻す（受音点ごとに繰り返さない）
    if project.receiver_index in (None, 1) and results.get("model") is not None:
        _update_condition_table(project, results["model"], verbose=verbose)

    if make_figures:
        if progress is not None:
            progress("図を書き出し中", None)
        written = plots.save_all(project, results, verbose=verbose)
        if verbose:
            print(f"[run] 図を {len(written)} 枚書き出しました → {project.figure_dir()}")

    # 実際に使った音源・受音点を project.json に残す（DXF から取った場合も分かるように）
    # ★顔の向きは**結果に持たせる**（受音点ごとに違うため）。
    #   project に書き戻すと、複数受音点のときにリストが 1 点ぶんの数値に潰れる
    results["head_azimuth"] = (project.head_azimuth_for(0)
                               if head_azimuth is None else float(head_azimuth))
    project.source = results["soundsource_point"].tolist()
    if write_back:
        project.receiver = results["receiver_point"].tolist()
    project.save()
    return results


def redraw(project, verbose=True):
    """**計算し直さずに**、保存済みの結果から図を一式描き直す。

    受音点が複数あれば**全点ぶん**描き直す（`結果/recN/` を順に見る）。
    顔の向きを直したあと伝搬方向の図だけ作り直したいときの入口。

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
    # 受音点が複数あるときは 1 点ずつ描き直す（`receiver_index` を立てて再帰）
    if project.receiver_index is None:
        import summary as sm
        folders = [name for name, _ in sm.receiver_folders(project)
                   if name.startswith("rec")]
        indexes = [int(name[3:]) for name in folders
                   if os.path.isdir(project.path(pj.RESULT_DIR, name))]
        if indexes:
            written = []
            for k in indexes:
                sub = _sub_project(project, k - 1)
                written.extend(redraw(sub, verbose=verbose))
            _write_summaries(project, verbose=verbose)
            return written

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
    absorption_table = _absorption_table_for(project)
    model = rd.read_model(project.dxf_path, band_number=project.band_number,
                          absorption_table=absorption_table, unit=project.unit,
                          orient_normals=project.orient_normals,
                          flip_faces=_flip_faces_for(project),
                          face_materials=_face_materials_for(project), verbose=False)

    results = {"model": model, "pulses": pulses, "frequencies": frequencies,
               "atmosphere": atmosphere, "impulse": None,
               "reverberation": None, "clarity": None, "statistical": None,
               "level": None, "sti": None}

    # ---- 音圧レベルと STI（パルス列から出るので描き直しでも作れる）----
    if len(pulses):
        import sound_level as sl
        source = (np.asarray(project.source, dtype=float)
                  if project.source is not None else None)
        receiver = (np.asarray(project.receiver, dtype=float)
                    if project.receiver is not None else None)
        distance = (float(np.linalg.norm(receiver - source))
                    if source is not None and receiver is not None else None)
        results["level"] = sl.band_levels(
            pulses.time, pulses.energy, pulses.distance, atmosphere, frequencies,
            source_power_db=project.source_power_db, source_distance=distance,
            verbose=False)
        results["sti"] = sl.speech_transmission_index(
            pulses.time, pulses.energy, pulses.distance, atmosphere, frequencies,
            source_power_db=project.source_power_db,
            noise_level_db=project.noise_level_db, verbose=False)
        # ★**すでにある CSV は書き換えない**（描き直しは図だけのはずなので）。
        #   ただし**無いものは作る**。あとから足した指標（音圧レベル・STI）を
        #   過去に計算したプロジェクトへ反映するには、これが唯一の道になる
        for key, write in (("spl", sl.write_levels), ("sti", sl.write_sti)):
            target = project.result_path(key)
            if not os.path.exists(project.existing_result_path(key)):
                write(target, results["level" if key == "spl" else "sti"])
                if verbose:
                    print(f"[redraw] 無かったので作りました: {target}")

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
              f"→ {project.figure_dir()}")
    return written


def _assignment_for(project):
    """レイヤ → 材料番号の対応。**条件表が最優先**（依頼 2026-08-21）。

    CAD のレイヤ名を書き換えずに材料を差し替えられるようにするための仕組み。
    表が無ければ従来どおり `project.assignment`（project.json）を使う。
    受音点ごとの子フォルダには置かないので、法線と同じ探し方をする。
    """
    import condition_table as ct
    owner = _owner_of(project, lambda p: ct.exists(p) or None) or project
    return ct.assignment_for(owner, verbose=False)


def _library_for(project, verbose=False):
    """材料の一覧。**条件表の「吸音率」シートが最優先**（依頼 2026-08-21）。

    PJ 固有の吸音データを条件表 1 ファイルに閉じ込められるようにするため。
    シートが無ければ従来どおり吸音率 CSV（`absorption_csv`）を読む。
    """
    import absorption as ab
    import condition_table as ct

    library = ct.library_from_book(ct.path(project), kind=project.absorption_kind,
                                  verbose=verbose)
    if library is not None:
        return library
    if project.absorption_path:
        return ab.MaterialLibrary.from_csv(project.absorption_path,
                                           kind=project.absorption_kind)
    return None


def _absorption_table_for(project, verbose=False):
    """`read_dxffile.read_model(absorption_table=...)` に渡す辞書。無ければ None。"""
    library = _library_for(project, verbose=verbose)
    if library is None:
        return None
    return library.absorption_table(_assignment_for(project),
                                    band_number=project.band_number,
                                    warn=verbose)


def _update_condition_table(project, model, verbose=True):
    """材料条件表を作る／更新する（面数・面積・吸音率の参考列を書き直す）。

    ★**利用者が書いた「材料名」は上書きしない**（`condition_table.update`）。
    計算のたびに更新するので、**そのとき実際に使った条件が表に残る**。
    """
    try:
        import condition_table as ct
        return ct.update(project, model, _library_for(project),
                         _assignment_for(project), verbose=verbose)
    except Exception as error:     # 表が作れなくても計算結果は残す
        print(f"[run] 材料条件表を更新できませんでした: "
              f"{type(error).__name__}: {error}")
        return None


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
    p.add_argument("--conditions", nargs="*", default=None, metavar="CSV",
                   help="材料条件表を指定して**一括で回す**。"
                        "ファイル名を並べるか、値なしでフォルダ内の条件表を全部")
    p.add_argument("--no-reuse", action="store_true",
                   help="保存した経路を使わず、音線追跡からやり直す")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if not project.dxf:
        raise SystemExit(f"{a.folder} に project.json が無いか、DXF が設定されていません。"
                         f"先に app.py で条件を入力してください")
    if a.redraw:
        redraw(project)
        return
    if a.conditions is not None:
        run_conditions(project, a.conditions or None,
                       make_figures=not a.no_figures)
        return
    run(project, make_figures=not a.no_figures, reuse_paths=not a.no_reuse)


if __name__ == "__main__":
    main()
