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

# 計算に使った PC と所要時間を残すファイル（結果一式の「概要」に並ぶ）
RUN_INFO_FILE = "計算情報.csv"
from atmosphere import Atmosphere


def _freeze_condition_name(project, verbose=True):
    """**この実行で使う条件名を、計算を始める前に確定させる。**

    ★条件シートを選ばずに計算しても、材料は条件表の**先頭シート**のものが使われる
      （`condition_table.sheet_of`）。名前が実態と食い違わないよう、
      `Project.condition_label` はそのシート名へフォールバックする
      （2026-09-06 ユーザー判断・不具合報告 WIN240377 の ②）。

    ★★**ここで `condition_sheet` に焼き付けるのが肝**。条件表は計算の**途中**で
      作られる（`_update_condition_table` は面の情報が要るので `process` のあと）。
      焼き付けないと、同じ 1 回の実行の中で
      「表ができる前に書いた `室_吸音率と理論値.csv`」と
      「表ができた後に書いた `室_現状_まとめ_….csv`」が**混在する**（実際に踏んだ）。

    ★**`project.json` には書き戻さない**（`save()` のあとに立てる）。
      書き戻すと、利用者が選んでいないシート名が設定として残ってしまう。
      受音点ごとの子プロジェクトへは `DEFAULTS` 経由でそのまま伝わる。
    """
    import condition_table as ct

    if project.condition_sheet:
        return
    if _stem_is_named(project):
        return              # 条件表のファイル名が条件名（従来どおり）
    sheet = project._fallback_sheet()
    if not sheet:
        # ★条件表がまだ無い＝**この実行の中で作られる**
        #   （`_update_condition_table` → `ct.update`。面の情報が要るので後半）。
        #   そのとき付く名前（`ct.FIRST_SHEET`）を**先取りする**。
        #   先取りしないと「1 回目は名前なし・2 回目から名前あり」になり、
        #   描き直し（`--redraw`）やまとめ表の作り直しで名前が変わってしまう。
        #   作られるシートは**そのとき実際に使った材料**を記録したものなので、
        #   1 回目にこの名前を付けても実態と食い違わない
        if not ct.is_book(project.condition_path):
            return
        sheet = ct.FIRST_SHEET
    project.condition_sheet = sheet
    if verbose:
        print(f"[run] 条件シートの指定が無いので『{sheet}』として扱います"
              f"（結果ファイル名にもこの名前が付きます）")


def _stem_is_named(project):
    """条件表のファイル名そのものが条件名になっているか（`条件A.xlsx` など）。"""
    stem = os.path.splitext(os.path.basename(project.condition_path or ""))[0]
    return bool(stem) and stem not in pj.DEFAULT_CONDITION_STEMS


def run(project, verbose=True, make_figures=True, progress=None,
        reuse_paths=True, save_settings=True):
    """プロジェクトの条件で計算し、結果 CSV と図を書き出す。

    受音点が複数あるときは、**音線追跡を 1 回で済ませて**受音判定だけ点ごとに行う
    （F-6。追跡は受音点に依らない）。結果は受音点ごとに `結果/recN/`・`図/recN/`、
    受音点に依らないもの（室の吸音と理論値・音線軌跡）は `結果/` 直下へ書く。
    ファイル名の頭には対象室＋条件名（`project.name`）が付く。

    ★**音源が複数あるときは 1 点ずつ回す**（2026-09-15。不具合報告 ⑨）。
    結果は `結果/src1/recN/` `結果/src2/recN/` …に分かれ、そのあと
    `source_mix` が見方（平均・重ね合わせ）を作る。
    **音源が 1 点だけなら置き場は従来どおり**（`結果/recN/`）。
    """
    project.ensure_dirs()
    if save_settings:
        project.save()  # 実行した条件を必ず残す（あとで再現できるように）
    _freeze_condition_name(project, verbose=verbose)

    dxf = project.dxf_path
    if not dxf or not os.path.exists(dxf):
        raise FileNotFoundError(f"DXF が見つかりません: {project.dxf!r}")

    # ★★**音源が複数あれば 1 点ずつ回す**（2026-09-15 ユーザー要望。不具合報告 ⑨）。
    #   それまでは `model.source_points[0]` を黙って使い、2 点目以降を捨てていた。
    #   結果は `結果/srcM/recN/` に分かれ、そのあと `source_mix` が合成を作る
    sources = _sources(project)
    if len(sources) > 1:
        return _run_sources(project, sources, verbose=verbose,
                            make_figures=make_figures, progress=progress,
                            reuse_paths=reuse_paths)
    return _run_source(project, verbose=verbose, make_figures=make_figures,
                       progress=progress, reuse_paths=reuse_paths,
                       save_settings=save_settings)


def _sources(project, model=None):
    """計算する音源の一覧。`project.source` の指定が最優先、無ければ DXF から。

    ★**全部返す**のが肝（不具合報告 ⑨）。`_source_of()` は 1 点しか返さないので、
    音源が 2 点あっても 2 点目が黙って捨てられていた。
    """
    import read_dxffile as rd

    if project.source is not None:
        return [np.asarray(project.source, dtype=float)]
    if model is None:
        model = _ordered(project,
                         rd.read_model(project.dxf_path, unit=project.unit,
                                       band_number=project.band_number,
                                       verbose=False), verbose=False)
    return [np.asarray(p, dtype=float)
            for p in (getattr(model, "source_points", None) or [])]


def _ordered(project, model, verbose=True):
    """**測定点の並び**（`測定点順.json`）を当てる（2026-09-15 ユーザー要望）。

    ★DXF に出てきた順のままだと `結果/recN/` の N が CAD のラベルとずれる
    （実案件で `rec1` が R4 になっていた）。読んだ直後に 1 回だけ通す。
    """
    import point_order as po

    try:
        model = po.apply(project, model, verbose=verbose)
    except Exception as error:      # 並びが当たらなくても計算は続けられる
        print(f"[run] 測定点の並びを当てられませんでした: "
              f"{type(error).__name__}: {error}")
        return model
    # ★★**前回の計算と並びが変わっていたら知らせる**（2026-09-20。不具合報告 ⑯）。
    #   DXF を作り直すと `POINT` の出てくる順が変わりうるので、黙って
    #   `結果/src1/` `結果/rec1/` の中身が別の点にすり替わる。
    #   ここで気づけば 50 分の無駄打ち（経路の使い回しも効かなくなる）を防げる
    try:
        po.check_against_previous(project, model, applied=True, verbose=verbose)
    except Exception as error:      # 突き合わせに失敗しても計算は続けられる
        print(f"[run] 前回の測定点と突き合わせられませんでした: "
              f"{type(error).__name__}: {error}")
    return model


def _run_sources(project, sources, verbose=True, make_figures=True,
                 progress=None, reuse_paths=True):
    """**音源を 1 点ずつ**回し、そのあと合成（`source_mix`）を作る。

    ★**`project.json` に音源を書き戻さない**（`save_settings=False`）。
    書き戻すと `source` が 1 点に固定され、**次回から 1 点目しか回らなくなる**
    （2026-09-06 に一括実行で踏んだのと同じ落とし穴）。

    ★測定点の一覧（`結果/<室>_測定点.csv`）は**全部の音源を 1 枚に**並べたいので、
    音源ごとの `_run_source` では書かず、ここで最後に 1 回だけ書く。
    """
    import source_mix as sx

    if verbose:
        print(f"[run] ★音源が {len(sources)} 点あります。"
              f"**1 点ずつ**計算して `結果/src1` `結果/src2` … に分けます"
              f"（ISO 3382 は音源位置ごとに測る）")
    results = []
    for order, point in enumerate(sources):
        sub = _sub_source(project, order, point)
        if verbose:
            print("")
            print(f"[run] ══ 音源 {order + 1}/{len(sources)} "
                  f"{np.round(point, 3).tolist()} → 結果/{sub.source_folder}/")
        results.append(_run_source(
            sub, verbose=verbose, make_figures=make_figures,
            reuse_paths=reuse_paths, save_settings=False, write_points=False,
            progress=_prefixed(progress, f"音源{order + 1}/{len(sources)} ")))

    _write_points(project, model=results[0].get("model"), sources=sources,
                  verbose=verbose)
    try:
        sx.write_all(project, verbose=verbose)
    except Exception as error:      # 合成が作れなくても音源ごとの結果は残る
        print(f"[run] 合成（音源のまとめ方）を作れませんでした: "
              f"{type(error).__name__}: {error}")
    return {"sources": sources, "source_results": results, **results[0]}


def _sub_source(project, order, point):
    """音源 `order` 番目（0 始まり）を扱う `Project`。結果は `結果/srcM/` へ。"""
    sub = pj.Project(project.folder,
                     **{k: getattr(project, k) for k in pj.DEFAULTS})
    sub.source = np.asarray(point, dtype=float).tolist()
    sub.source_index = order + 1
    return sub


def _run_source(project, verbose=True, make_figures=True, progress=None,
                reuse_paths=True, save_settings=True, write_points=True):
    """**音源 1 点ぶん**の計算（受音点は全部）。従来の `run()` の中身。"""
    receivers = _receivers(project)
    if len(receivers) <= 1:
        # 1 点でも `結果/rec1/` に入れる（点数によって置き場が変わらないように）
        result = _run_one(_sub_project(project, 0),
                          receivers[0] if receivers else None,
                          verbose=verbose, make_figures=make_figures,
                          write_back=False,
                          head_azimuth=project.head_azimuth_for(0),
                          reuse_paths=reuse_paths, progress=progress,
                          save_settings=save_settings)
        if write_points:
            _write_points(project, receivers, result.get("model"), verbose=verbose)
        _write_summaries(project, verbose=verbose)
        return result

    # ★保存した経路が全受音点ぶんそろっていれば、音線追跡そのものを省く（F-9）。
    #   吸音材だけ変えた計算はここで終わり（あとはエネルギーの掛け算だけ）
    if reuse_paths and _paths_ready(project, receivers, verbose=verbose):
        results, shared_statistical = [], None
        for k, point in enumerate(receivers):
            sub = _sub_project(project, k)
            if verbose:
                print("")
                print(f"[run] ── 受音点 {k + 1}/{len(receivers)}"
                      f"（保存した経路から再開）")
            results.append(_run_one(sub, point, verbose=verbose,
                                    make_figures=make_figures, write_back=False,
                                    head_azimuth=project.head_azimuth_for(k),
                                    reuse_paths=True, save_settings=save_settings,
                                    statistical_result=shared_statistical,
                                    progress=_prefixed(progress,
                                                       f"受音点{k + 1}/{len(receivers)} ")))
            # 統計残響式は受音点に依らないので 1 点目の結果を配る（無駄なループを消す）
            shared_statistical = (shared_statistical
                                  or results[-1].get("statistical"))
        if write_points:
            _write_points(project, receivers, results[0].get("model"),
                          verbose=verbose)
        _write_summaries(project, verbose=verbose)
        return {"receivers": receivers, "results": results, **results[0]}

    # ★音線追跡は**受音点に依らない**ので 1 回だけ回し、受音点ごとに配る（F-6）。
    #   受音しても音線は打ち切られないため、受音球をいくつ置いても追跡は同じ。
    #   受音点ごとに追い直していたときは、その回数だけ全部やり直していた
    if verbose:
        print(f"[run] 受音点が {len(receivers)} 点あります。"
              f"音線追跡は 1 回で済ませ、受音判定だけ {len(receivers)} 点ぶん行います")
    traced, recorder = _trace_once(project, receivers, verbose=verbose, progress=progress)

    results, shared_statistical = [], None
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
                                reuse_paths=False, save_settings=save_settings,
                                statistical_result=shared_statistical,
                                progress=_prefixed(progress,
                                                   f"受音点{k + 1}/{len(receivers)} ")))
        shared_statistical = shared_statistical or results[-1].get("statistical")
    if recorder is not None:
        # 軌跡は受音点に依らないので `結果/` 直下に 1 つだけ置く。
        # `clear_results` のあとに置かないと消される
        recorder.save_npz(project.result_path("raylog"))
    if write_points:
        _write_points(project, receivers, results[0].get("model"), verbose=verbose)
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
        # ★**面の上に置いた音源は置き直した位置で指紋を取る**（2026-08-24）。
        #   `procedure` 側は置き直してから指紋を取るので、ここで合わせないと
        #   「経路が無い」と判断してしまい、条件を変えるたびに追跡し直しになる
        placement = placement_for(project, model, source, verbose=False)
        if placement.on_surface:
            source = placement.point
        for index, point in enumerate(receivers):
            sub = _sub_project(project, index)
            mark = pc.fingerprint(model.mesh, faces, source, point, project.rays,
                                  project.nref, project.radius, project.two_sided)
            # ★測定点を全部渡す。位置が合わないとき「並びが入れ替わっただけ」
            #   かどうかを見て対処法を添えてもらう（2026-09-20。不具合報告 ⑯）
            known = {"source": list(model.source_points or []),
                     "receiver": list(model.receiver_points or [])}
            if pc.load(sub.paths_cache(), mark, verbose=False,
                       points=known) is None:
                if verbose:
                    # 理由は `pc.compare` が出す。もう一度呼んで表示させる
                    pc.load(sub.paths_cache(), mark, verbose=True, points=known)
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
    model = rd.read_model(project.dxf_path, unit=project.unit,
                          absorption_table=table,
                          orient_normals=project.orient_normals,
                          band_number=project.band_number,
                          flip_faces=_flip_faces_for(project),
                          face_materials=_face_materials_for(project),
                          verbose=verbose)
    return _ordered(project, model, verbose=verbose)


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
        # ★★**条件ごとの複製から `project.json` を書かせない**
        #   （2026-09-06 不具合報告 WIN240377 の ③）。
        #   `sub` は `folder` が同じなので、`run()` の中の `save()` が
        #   **その条件の設定で project.json を上書きし、最後の条件が残っていた**。
        #   次に単発で回すと、選んでいたつもりの条件と違うものが使われる。
        #   設定はループのあとに**利用者が選んでいた条件のまま**1 回だけ書く
        results.append(run(sub, verbose=verbose, make_figures=make_figures,
                           progress=stage, save_settings=False))
        done.append((file_name, sheet))

    # 設定を 1 回だけ残す。★条件は**呼ばれたときのまま**（一括で回した最後の条件に
    #   すり替えない）。音源だけは DXF から取った値を書き戻しておく
    if results and results[0].get("soundsource_point") is not None:
        project.source = np.asarray(results[0]["soundsource_point"]).tolist()
    project.save()

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


def log_path(project):
    """計算中のログを残すテキストのパス（`結果/<室>_計算ログ.txt`）。

    画面のログは閉じると消えてしまうので、あとから追えるように残す
    （2026-08-21 ユーザー要望）。条件ごとに分けたいので条件名も付く。
    """
    return os.path.join(project.path(pj.RESULT_DIR),
                        project.prefixed("計算ログ.txt"))


def machine_info():
    """計算に使った PC の情報（2026-08-21 ユーザー要望）。

    ★**標準ライブラリだけで取れるもの**にしてある（依存を増やさないため）。
    メモリと CPU 名は Windows なら `wmic` / 環境変数から拾えるので試し、
    取れなければその項目だけ空にする。
    """
    import platform

    info = {
        "PC 名": platform.node(),
        "OS": f"{platform.system()} {platform.release()}（{platform.version()}）",
        "CPU": platform.processor() or os.environ.get("PROCESSOR_IDENTIFIER", ""),
        "論理コア数": os.cpu_count(),
        "Python": platform.python_version(),
    }
    try:
        import numpy as np
        info["numpy"] = np.__version__
    except Exception:
        pass
    # 物理メモリ。Windows は GlobalMemoryStatusEx、それ以外は os.sysconf
    try:
        if os.name == "nt":
            import ctypes

            class Status(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong),
                            ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong),
                            ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong),
                            ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong),
                            ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

            status = Status()
            status.dwLength = ctypes.sizeof(Status)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
            info["メモリ [GB]"] = round(status.ullTotalPhys / 1024 ** 3, 1)
        else:
            info["メモリ [GB]"] = round(
                os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                / 1024 ** 3, 1)
    except Exception:
        pass
    return info


def write_run_info(project, elapsed=None, conditions=None, verbose=True):
    """**使った PC と計算にかかった時間**を残す（`結果/<室>_計算情報.csv`）。

    結果一式（Excel）の「概要」シートにそのまま並ぶ（2026-08-21 ユーザー要望）。
    「同じ条件で誰の PC で何分かかったか」を後から見られるようにするためのもの。
    """
    import csv
    import time

    rows = [("項目", "値"), ("計算日時", time.strftime("%Y-%m-%d %H:%M:%S"))]
    if elapsed is not None:
        rows.append(("計算にかかった時間 [s]", round(float(elapsed), 1)))
        rows.append(("計算にかかった時間", _clock(elapsed)))
    if conditions:
        rows.append(("回した条件の数", conditions))
    rows += [("音線数", project.rays), ("最大反射回数", project.nref),
             ("受音球の半径 [m]", project.radius),
             ("応答の長さ [s]", project.max_time),
             ("周波数バンド", project.band_number)]
    rows += list(machine_info().items())

    path = os.path.join(project.path(pj.RESULT_DIR),
                        project.prefixed(RUN_INFO_FILE))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        csv.writer(f).writerows(rows)
    if verbose:
        print(f"[run] 計算情報（PC と所要時間）: {path}")
    return path


def _clock(seconds):
    seconds = int(max(0.0, float(seconds)))
    hours, rest = divmod(seconds, 3600)
    minutes, second = divmod(rest, 60)
    return (f"{hours}:{minutes:02d}:{second:02d}" if hours
            else f"{minutes}:{second:02d}")


def _closed_expected(project):
    """★**閉じた室を想定しているか**（2026-08-28 ユーザー要望）。

    > そもそも閉じてる必要のないモデルもあるので、閉じたモデルを想定する場合に限ります。

    `自動` は None（決めつけず、どちらの可能性も伝える）。
    """
    value = str(getattr(project, "closed_model", pj.CLOSED_AUTO) or "").strip()
    if value == pj.CLOSED_YES:
        return True
    if value == pj.CLOSED_NO:
        return False
    return None


def _receiver_groups(project, receivers):
    """受音点のレイヤ名（測線）。取れなければ空。"""
    try:
        points, names = _receivers(project, groups=True)
    except Exception:
        return []
    if len(points) != len(receivers):
        return []
    return names


def _write_points(project, receivers=None, model=None, verbose=True, sources=None):
    """**測定点の一覧（CSV）と配置図（平面＋立面 2 方向）**を書く。

    ★どれがどの点でどちらを向いているか、あとから分かるように
    （2026-08-23 ユーザー要望）。**受音点に依らない**ので `結果/` `図/` 直下、
    吸音材にも依らないので条件名は付けない（`ROOM_SCOPED_RESULTS`）。
    名前は `結果/recN/` に合わせて `rec1`… とする（表とフォルダが対応する）。
    """
    try:
        if receivers is None:
            receivers = _receivers(project)
        # ★音源は**全部**並べる（2026-09-15。不具合報告 ⑨）。
        #   呼び出し側が音源ごとに回しているときは、その一覧をそのまま渡してもらう
        if sources is None:
            sources = _sources(project, model)
        sources = [np.asarray(p, dtype=float) for p in sources]
        azimuths = [project.head_azimuth_for(k) for k in range(len(receivers))]

        shared = _sub_project(project, None) if False else project
        path = pj.Project(project.folder, **project.values()).result_path("points") \
            if False else None
        # `receiver_index` を外した状態で書く（`結果/` 直下・条件名なし）
        keep = project.receiver_index
        project.receiver_index = None
        try:
            written = [pj.write_points_csv(project.result_path("points"),
                                           sources, receivers, azimuths,
                                           groups=_receiver_groups(project,
                                                                   receivers))]
            if verbose:
                print(f"[run] 測定点の一覧: {written[0]}")
            try:
                import plots as pl
                figure = project.figure_path("points.png", shared=True)
                pl.measurement_points(model, figure, sources=sources,
                                      receivers=receivers, azimuths=azimuths)
                written.append(figure)
                if verbose:
                    print(f"[run] 測定点の配置図（平面＋立面 2 方向）: {figure}")
            except Exception as error:      # 図が作れなくても表は残す
                print(f"[run] 測定点の配置図を作れませんでした: "
                      f"{type(error).__name__}: {error}")
            return written
        finally:
            project.receiver_index = keep
    except Exception as error:              # 本体の結果は残す
        print(f"[run] 測定点の一覧を作れませんでした: "
              f"{type(error).__name__}: {error}")
        return []


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
        model = _ordered(project,
                         rd.read_model(project.dxf_path, unit=project.unit,
                                       absorption_table=table,
                                       orient_normals=project.orient_normals,
                                       band_number=project.band_number,
                                       flip_faces=_flip_faces_for(project),
                                       face_materials=_face_materials_for(project),
                                       verbose=False), verbose=False)
        source = (project.source if project.source is not None
                  else (model.source_points[0] if model.source_points else None))
        if source is None:
            return None
        rays = sr.soundray_generator(project.rays)
        # ★面の上に置いた音源なら半球に折り返し、音源も面の上へ置き直す
        placement = placement_for(project, model, source, verbose=verbose)
        if placement.on_surface:
            import source_placement as spl
            rays = spl.hemisphere(rays, placement.normal)
            source = placement.point
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


def placement_for(project, model, source, verbose=True):
    """音源の置かれ方（面の上か）を調べる。→ `source_placement.Placement`

    ★2026-08-24 ユーザー要望。設定は `project.json` の
    `source_on_surface` / `source_surface_tolerance` / `source_direction`。
    """
    import source_placement as spl

    place = spl.detect(np.asarray(source, dtype=float), model.mesh,
                       tolerance=getattr(project, "source_surface_tolerance",
                                         spl.DEFAULT_TOLERANCE),
                       direction=getattr(project, "source_direction", None),
                       enabled=bool(getattr(project, "source_on_surface", True)))
    if verbose:
        print(f"[run] 音源の置かれ方: {place.describe()}")
    return place


def _prefixed(progress, prefix):
    """受音点が複数あるとき、どの受音点の処理かを段階名に添える。"""
    if progress is None:
        return None
    return lambda stage, fraction=None: progress(prefix + stage, fraction)


def _receivers(project, groups=False):
    """計算する受音点の一覧。project.receiver の指定が最優先、無ければ DXF から。

    `groups=True` なら **(点, レイヤ名)** の 2 つを返す（受音点を方向別に
    レイヤ分けしたモデル用。2026-08-24）。
    ★**音源と重なる点は外す**（そこでは音圧が発散する）。半無響室のモデルは
      測線の起点（＝音源の位置）が受音点のレイヤに入っていた。
    """
    import read_dxffile as rd
    if project.receiver is not None:
        points = [np.asarray(project.receiver, dtype=float)]
        names = [""]
    else:
        probe = _ordered(project,
                         rd.read_model(project.dxf_path, unit=project.unit,
                                       band_number=project.band_number,
                                       verbose=False), verbose=False)
        points = [np.asarray(p, dtype=float) for p in probe.receiver_points]
        names = list(getattr(probe, "receiver_layer_names", []) or
                     ["" for _ in points])
        source = project.source
        if source is None and probe.source_points:
            source = probe.source_points[0]
        if source is not None:
            source = np.asarray(source, dtype=float)
            keep = [k for k, p in enumerate(points)
                    if float(np.linalg.norm(p - source)) > 1.0e-6]
            if len(keep) != len(points):
                print(f"[run] 音源と重なる受音点を "
                      f"{len(points) - len(keep)} 点はずしました"
                      f"（その位置では音圧が発散します）")
            points = [points[k] for k in keep]
            names = [names[k] for k in keep] if names else names
    return (points, names) if groups else points


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
    # ★音源の棚（`結果/srcM/`）は引き継ぐ（不具合報告 ⑨）。
    #   引き継がないと、音源ごとに回しているのに結果が 1 か所へ重なって書かれる
    sub.source_index = project.source_index
    sub.source_tag = project.source_tag
    # ★**名前は変えない。**`name` は結果ファイル名の頭に付く（対象室＋条件名）ので、
    #   受音点ごとに変えるとファイル名が受音点ごとに違ってしまう。
    #   何点目かは `receiver_index` が持っていて `Project.summary()` が表示する
    return sub


def _run_one(project, receiver, verbose=True, make_figures=True,
             write_back=True, head_azimuth=None, traced_history=None,
             reuse_paths=True, statistical_result=None, progress=None,
             save_settings=True):
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
        # ★★**条件表の安全率まで効かせた表をここで渡す**（2026-09-20。不具合報告 ⑱）。
        #   渡さないと `procedure` が material_library から組み立て直すので、
        #   **安全率が掛からないまま計算される**（危険側・警告なし）
        absorption_table=_absorption_table_for(project, verbose=verbose),
        band_number=project.band_number,
        # ★帯域の幅（1/1 か 1/3）と下端（2026-08-26）
        band_width=getattr(project, "band_width", "1/1"),
        band_start=getattr(project, "band_start", None),
        unit=project.unit,
        orient_normals=project.orient_normals,
        two_sided=project.two_sided,
        volume=project.volume,
        # ★面の上に置いた音源（半無響室の床置きなど。2026-08-24）
        level_method=getattr(project, "level_method", "both"),
        source_on_surface=getattr(project, "source_on_surface", True),
        source_surface_tolerance=getattr(project, "source_surface_tolerance", 0.0),
        source_direction=getattr(project, "source_direction", None),
        # 閉じた室を想定しているか（作図チェックの言い方が変わる）
        closed_expected=_closed_expected(project),
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
        # ★RTany（減衰曲線をどこで読むか。2026-09-15 ユーザー指示）
        rt_any_start_db=getattr(project, "rt_any_start_db", None),
        rt_any_end_db=getattr(project, "rt_any_end_db", None),
        decay_fit=getattr(project, "decay_fit", None),
        room_filename=project.result_path("room"),
        clarity_filename=project.clarity_path(),
        level_filename=project.result_path("spl"),
        sti_filename=project.result_path("sti"),
        paths_filename=project.paths_cache(),
        reuse_paths=reuse_paths,
        # 統計残響式は受音点に依らない。2 点目以降は 1 点目の結果を使う
        statistical_result=statistical_result,
        source_power_db=project.source_power_db,
        noise_level_db=project.noise_level_db,
        statistical=project.statistical,
        # ★合成のやり方（`fast` / `exact`）。波動解と突き合わせるときは exact
        impulse_method=getattr(project, "impulse_method", "fast"),
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
    if save_settings:
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
    # ★音源が複数あるときは**棚ごと**に描き直す（2026-09-15。不具合報告 ⑨）。
    #   `結果/src1/` `結果/src2/` …と合成の棚を順に見る。合成の「平均」には
    #   パルス列が無いので、描き直せない棚は知らせて飛ばす
    if not project.source_folder and project.receiver_index is None:
        shelves = project.source_folders()
        if shelves:
            written = []
            for tag in shelves:
                shelf = pj.Project(project.folder,
                                   **{k: getattr(project, k) for k in pj.DEFAULTS})
                shelf.source_tag = tag
                if verbose:
                    print(f"[run] ── 描き直し: 結果/{tag}/")
                try:
                    written.extend(redraw(shelf, verbose=verbose))
                except (ValueError, FileNotFoundError) as error:
                    import source_mix as sx
                    if tag == sx.FOLDERS[sx.MIX_AVERAGE]:
                        # ★「平均」は**指標だけ**の棚（波形が無いので図は作れない）
                        print(f"[run] 結果/{tag}/ は指標だけの棚なので図は作りません"
                              f"（音源ごとの結果を平均したもの）")
                    else:
                        print(f"[run] 結果/{tag}/ は描き直せません: {error}")
            return written

    # 受音点が複数あるときは 1 点ずつ描き直す（`receiver_index` を立てて再帰）
    if project.receiver_index is None:
        import summary as sm
        folders = [name for name, _ in sm.receiver_folders(project)
                   if name.startswith("rec")]
        root = sm.results_root(project)     # 音源が複数なら `結果/srcM/`
        indexes = [int(name[3:]) for name in folders
                   if os.path.isdir(os.path.join(root, name))]
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
    frequencies = ab.frequency_bands(project.band_number,
                                     getattr(project, "band_width", "1/1"),
                                     getattr(project, "band_start", None))

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
    if not len(model.mesh):
        # ★面 0 枚のまま進むと分かりにくい例外になる（2026-08-21）
        raise ValueError(
            f"{os.path.basename(project.dxf_path or '')} から面が 1 枚も"
            f"読めませんでした（REGION と 3DSOLID は読めません）。"
            f"モデルの指定を確かめてください")

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
            results["impulse"][0], results["impulse"][1], frequencies=frequencies,
            # ★描き直しでも読み方を合わせる（図と CSV が食い違わないように）
            measures=rv.measures_with_any(
                getattr(project, "rt_any_start_db", None),
                getattr(project, "rt_any_end_db", None)),
            fit=getattr(project, "decay_fit", None) or rv.DEFAULT_DECAY_FIT)
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
        # ★xlsx（「吸音率」シート）も CSV も選べる（2026-08-21 ユーザー要望）
        return ab.MaterialLibrary.from_file(project.absorption_path,
                                            kind=project.absorption_kind,
                                            verbose=verbose)
    return None


def _absorption_table_for(project, verbose=False):
    """`read_dxffile.read_model(absorption_table=...)` に渡す辞書。無ければ None。

    ★条件表の**安全率**（例 0.8 掛け）もここで効かせる。
    カタログ値に掛けてから垂直入射へ変換する（`condition_table.absorption_table`）。

    ★★**この表を本計算にも渡すこと**（`procedure.process(absorption_table=…)`）。
    2026-09-20 より前は `procedure` が material_library から組み立て直していたので、
    ここで掛けた安全率が計算に届いていなかった（不具合報告 ⑱）。
    """
    import condition_table as ct

    library = _library_for(project, verbose=verbose)
    if library is None:
        # ★安全率が書いてあるのに材料一覧が無いなら**黙って捨てない**。
        #   吸音率を見過ぎる（危険側）方向に外れるので必ず知らせる
        factors = ct.factors_for(project, verbose=False)
        if factors:
            print(f"[run] ★安全率が {len(factors)} レイヤに書かれていますが、"
                  f"材料一覧（条件表の「吸音率」シート／吸音率表）が読めないので"
                  f"**効きません**。材料一覧を用意してください")
        return None
    return ct.absorption_table(library, _assignment_for(project),
                               factors=ct.factors_for(project, verbose=verbose),
                               band_number=project.band_number, warn=verbose)


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
