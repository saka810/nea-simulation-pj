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
    法線の確認（normal_editor）  計算（run_project）
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

    raylog = project.result_path("raylog")
    if not os.path.exists(raylog):
        print("[app] 音線軌跡が無いので可視化は開きません（先に計算してください）")
        return
    try:
        # ★受音した経路だけに絞らない。**音がどう広がるか**を見るのが目的なので、
        #   受音経路だけにすると（偏ってはいないが）本数が減って様子が分かりにくい。
        #   経路の確認をしたいときは view_rays.py に --received-only を付けて呼ぶ
        vr.view(project.dxf_path, raylog, mode="both",
                absorption=project.absorption_path,
                unit=project.unit, orient_normals=project.orient_normals,
                received_only=False, max_rays=60, max_reflection=4,
                colour="time", opacity=0.10, frames=300, point_size=7)
    except Exception:
        print("[app] 可視化ウィンドウでエラーが起きました:")
        traceback.print_exc()


def _report(project, results):
    """計算が終わったあとの要約を文字で出す。"""
    import project as pj

    print("\n" + "=" * 70)
    print(f"計算が終わりました → {project.folder}")
    print("=" * 70)

    rt = results.get("reverberation")
    stat = results.get("statistical")
    if rt is not None:
        names = list(rt["measures"])
        print("  " + "周波数".rjust(8) + "".join(f"{n:>9}" for n in names)
              + f"{'Sabine':>9}{'Eyring':>9}")
        for i, fc in enumerate(rt["frequencies"]):
            cells = "".join(f"{rt['measures'][n][i]:9.3f}" for n in names)
            extra = ""
            if stat is not None:
                extra = f"{stat['sabine'][i]:9.3f}{stat['eyring'][i]:9.3f}"
            print(f"  {fc:7.0f}Hz{cells}{extra}")

    if stat is None and project.statistical:
        print("\n  ※ 統計残響式（Sabine / Eyring）は計算できませんでした。")
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
    a = p.parse_args()

    import normal_editor
    import project as pj
    import run_project
    import setup_window

    project = pj.Project.load(a.folder or os.getcwd())

    if a.run:
        if not project.dxf:
            raise SystemExit(f"{project.folder} に project.json がありません。"
                             f"先に条件を入力してください（--run を外す）")
        results = run_project.run(project)
        _report(project, results)
        if not a.no_view:
            _visualise(project, results)
        return

    # 条件入力 → 法線確認 の往復。法線を保存したら入力に戻る
    while True:
        project, action = setup_window.ask(project=project)
        if action is None:
            print("[app] 閉じました")
            return
        if action == "normals":
            try:
                normal_editor.edit(project)
            except Exception:
                print("[app] 法線の確認ウィンドウでエラーが起きました:")
                traceback.print_exc()
            continue
        break

    try:
        results = run_project.run(project)
    except Exception:
        print("[app] 計算でエラーが起きました:")
        traceback.print_exc()
        return
    _report(project, results)
    if not a.no_view:
        _visualise(project, results)


if __name__ == "__main__":
    sys.exit(main())
