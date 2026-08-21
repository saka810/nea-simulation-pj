"""アプリの入口。条件入力 → 法線の確認 → 計算 → 可視化 を順に開く。

    cd geosim
    python app.py                       条件入力から始める
    python app.py "C:\\...\\JR"          そのプロジェクトを開いた状態で始める
    python app.py "C:\\...\\JR" --run     入力ウィンドウを出さずにすぐ計算する

流れ:

    ┌ 条件入力（setup_window）─────────────────┐
    │  モデル・吸音率・保存先・音線数・受音球・温度湿度   │
    └──┬──────────────────┬────────────────┘
       │「法線を確認…」        │「計算する ▶」
       ▼                      ▼
    面の確認（face_editor）      計算（run_project）
       │ 保存すると normals.json  │  結果 CSV ＋ 図 PNG をプロジェクトフォルダへ
       └──→ 条件入力に戻る       ▼
                              可視化（view_rays: 音線 ↔ 音粒子を Tab で切替）

**1 つのウィンドウで完結する GUI が最終形**だが、どんな情報が要るかが
まだ固まっていないので、いまは必要なウィンドウをそのつど開く形にしてある。
ここを差し替えれば済むよう、各ウィンドウは独立して呼べるようにしてある。
"""

import argparse
import os
import sys
import traceback


def _visualise(project, results=None):
    """計算結果の可視化ウィンドウを開く（音線 ↔ 音粒子を Tab で切替）。"""
    import project as pj
    import view_rays as vr

    raylog = project.result_path("raylog")     # 受音点に依らないので `結果/` 直下
    if not os.path.exists(raylog):
        print("[app] 音線軌跡が無いので可視化は開きません（先に計算してください）")
        return
    try:
        # ここは**経路がどう通ったか**を見るための画面。
        # 受音した経路の初期反射だけを到来時刻で色分けするのがいちばん読みやすい。
        #
        # ※「音線がどう飛ぶか」（球状に出ている様子）は室形状と関係ないので、
        #   条件入力の「音線の飛び方を見る」（view_directions.py）で見る。
        #   こちらで全音線を描くと室内が線で埋まって読めなくなる
        vr.view(project.dxf_path, raylog, mode="both",
                absorption=project.absorption_path,
                unit=project.unit, band_number=project.band_number,
                orient_normals=project.orient_normals,
                received_only=True, max_rays=60, max_reflection=4,
                colour="time", opacity=0.10, frames=300, point_size=7,
                save_dir=project.screenshot_dir())
    except Exception:
        print("[app] 可視化ウィンドウでエラーが起きました:")
        traceback.print_exc()


def _run_with_progress(project):
    """進捗ウィンドウを出しながら計算する。

    条件入力を閉じてから結果が出るまで**何も出ない時間**があり、
    動いているのか止まっているのか分からなかった（ユーザー指摘）ので挟んでいる。
    計算は別スレッドで走らせる（tkinter のイベントループを止めないため）。
    """
    import progress_window
    import run_project

    return progress_window.run_with_progress(
        f"{project.display_name} を計算中",
        lambda progress: run_project.run(project, progress=progress),
        subtitle=(f"音線 {project.rays} 本 / 最大反射 {project.nref} 回 / "
                  f"受音球 {project.radius} m"))


def _run_all_with_progress(project):
    """全条件の一括計算を、進捗ウィンドウを出しながら回す。"""
    import progress_window
    import run_project

    return progress_window.run_with_progress(
        f"{project.display_name} の全条件を計算中",
        lambda progress: run_project.run_conditions(project, progress=progress),
        subtitle="材料条件表を順に当てます（2 件目以降は経路を使い回すので速い）")


def _report_conditions(project, outcome):
    """一括計算のあとの要約。**条件ごとの代表値を並べる**。"""
    import os

    print(chr(10) + "=" * 70)
    print(f"全条件の計算が終わりました → {project.folder}")
    print("=" * 70)
    conditions = (outcome or {}).get("conditions") or []
    for condition in conditions:
        print(f"  ・{os.path.basename(condition)}")
    comparison = (outcome or {}).get("comparison")
    if comparison:
        print(f"{chr(10)}  条件の比較表: {comparison}")
    print(f"  結果 CSV / Excel : {project.path('結果')}")


def _report_saved(project):
    """保存済みの CSV から要約を出す（計算し直さずに前回の結果を見るとき）。

    計算直後の `_report` と違い、手元にあるのは CSV だけなので読み直して並べる。
    **同じ見え方にする**ことで「いま計算したのか、前回のものか」で
    読み方が変わらないようにしている。
    """
    import numpy as np
    import project as pj

    print("\n" + "=" * 70)
    print(f"前回の結果を読み込みました → {project.folder}")
    print("=" * 70)

    saved = pj.load_results(project)
    rt, stat = saved.get("rt"), saved.get("statistical")
    if rt is None:
        print("  残響指標の CSV がありません")
    else:
        # CSV の列名はそのままだと読みにくいので短くする
        short = {"EDT_s": "EDT", "T20_s": "T20", "T30_s": "T30",
                 "curvature_percent": "曲率%", "sabine_s": "Sabine",
                 "eyring_s": "Eyring", "eyring_knudsen_s": "Eyring-Knudsen"}
        rows = dict(rt["rows"])
        if stat is not None:
            for key in ("sabine_s", "eyring_s", "eyring_knudsen_s"):
                if key in stat["rows"]:
                    rows[key] = stat["rows"][key]

        # **周波数は横**（table.py の共通ルール）。CSV・図と向きを揃える
        frequencies = rt["frequencies"]
        print("  " + " " * 16 + "".join(f"{f:>10.0f}" for f in frequencies))
        for name, values in rows.items():
            print(f"  {short.get(name, name):>16}"
                  + "".join("       ---" if np.isnan(v) else f"{v:10.3f}"
                            for v in values))

    print(f"\n  結果 CSV : {project.path(pj.RESULT_DIR)}")
    print(f"  図 PNG   : {project.path(pj.FIGURE_DIR)}")


def _report(project, results):
    """計算が終わったあとの要約を文字で出す。"""
    import project as pj

    print("\n" + "=" * 70)
    print(f"計算が終わりました → {project.folder}")
    print("=" * 70)

    import numpy as np
    import reverberation as rv

    rt = results.get("reverberation")
    stat = results.get("statistical")
    if rt is not None:
        # **周波数は横**（table.py の共通ルール）。CSV・図と向きを揃える
        rows = dict(rt["measures"])
        rows["曲率%"] = rt["curvature"]
        if stat is not None:
            for key, label in rv.STATISTICAL_LABELS.items():
                rows[label] = stat[key]
        print("  " + " " * 16 + "".join(f"{f:>10.0f}" for f in rt["frequencies"]))
        for name, values in rows.items():
            print(f"  {name:>16}"
                  + "".join("       ---" if np.isnan(v) else f"{v:10.3f}"
                            for v in values))

    # 音圧レベルと STI（帯域別の表は CSV に、ここでは代表値だけ）
    level = results.get("level")
    if level is not None:
        kind = "相対値" if level["relative"] else "絶対値"
        print("")
        print(f"  音圧レベル（{kind}）: 帯域合成 {level['overall']:.1f} dB / "
              f"A 特性 {level['overall_a']:.1f} dB(A)")
        print(f"    自由音場（逆二乗）との差: "
              + " / ".join(f"{f:.0f}Hz {d:+.1f}" for f, d
                           in zip(level["frequencies"], level["excess"])))
    sti = results.get("sti")
    if sti is not None:
        print(f"  STI: {sti['sti']:.3f}（{sti['rating']}）")

    if stat is None and project.statistical:
        print("\n  ※ 統計残響式（Sabine / Eyring-Knudsen）は計算できませんでした。")
        print("     このモデルは閉じていないので容積が自動で決まりません。")
        print("     条件入力の「室容積」に値を入れると比較できます"
              "（床面積 × 天井高で構いません）。")

    print(f"\n  結果 CSV : {project.path(pj.RESULT_DIR)}")
    print(f"  図 PNG   : {project.path(pj.FIGURE_DIR)}")


def main():
    p = argparse.ArgumentParser(description="幾何音響シミュレーション（GUI）")
    p.add_argument("folder", nargs="?", help="プロジェクトフォルダ")
    p.add_argument("--run", action="store_true",
                   help="条件入力ウィンドウを出さずにすぐ計算する")
    p.add_argument("--no-view", action="store_true", help="計算後に可視化を開かない")
    p.add_argument("--all-conditions", action="store_true",
                   help="フォルダ内の材料条件表を全部回す（--run と併用）")
    a = p.parse_args()

    import face_editor
    import project as pj
    import run_project
    import setup_window

    project = pj.Project.load(a.folder or os.getcwd())

    if a.run:
        if not project.dxf:
            raise SystemExit(f"{project.folder} に project.json がありません。"
                             f"先に条件を入力してください（--run を外す）")
        if a.all_conditions:
            outcome = run_project.run_conditions(project)
            _report_conditions(project, outcome)
            if not a.no_view:
                _visualise(project)
            return
        results = run_project.run(project)
        _report(project, results)
        if not a.no_view:
            _visualise(project, results)
        return

    # 条件入力 → 法線確認 / 結果表示 の往復。どちらも終わったら入力に戻る
    while True:
        project, action = setup_window.ask(project=project)
        if action is None:
            print("[app] 閉じました")
            return
        if action == "normals":
            try:
                face_editor.edit(project)
            except Exception:
                print("[app] 面の確認ウィンドウでエラーが起きました:")
                traceback.print_exc()
            continue
        if action == "view":
            # 計算し直さずに、保存済みの結果を見る。
            # **図も描き直す**（図の作り方を直したあと、古いプロジェクトを開いても
            # 前の見た目のままだった。音線追跡はやり直さないので数秒で済む）
            try:
                run_project.redraw(project)
            except Exception:
                print("[app] 図を描き直せませんでした（保存済みの図をそのまま使います）:")
                traceback.print_exc()
            _report_saved(project)
            if not a.no_view:
                _visualise(project)
            continue
        break

    # 全条件の一括計算（材料条件表を順に当てて回す。経路は使い回す）
    if action == "run_all":
        try:
            outcome = _run_all_with_progress(project)
        except Exception:
            print("[app] 一括計算でエラーが起きました:")
            traceback.print_exc()
            return
        _report_conditions(project, outcome)
        if not a.no_view:
            _visualise(project)
        return

    try:
        results = _run_with_progress(project)
    except Exception:
        print("[app] 計算でエラーが起きました:")
        traceback.print_exc()
        return
    _report(project, results)
    if not a.no_view:
        _visualise(project, results)


if __name__ == "__main__":
    try:
        code = main()
    finally:
        # **VTK の持ち物を文脈が生きているうちに片付ける。**
        # 放っておくとプロセス終了時に segfault することがある（2026-08-19）
        try:
            import view_model_gui as vg
            vg.close_all()
        except Exception:
            pass
    sys.exit(code)
