# -*- coding: utf-8 -*-
"""画角（カメラ）と表示の設定を残して、あとで読み直す。

★2026-08-24 ユーザー要望

> 条件違いで計算した結果を比較したい時に，音線の結果など画角などを合わせたい。
> 画角などの情報をセーブ・ロードできる機能とボタンを追加しておいて

条件（吸音材）を変えて計算し直すと画面は**開き直し**になるので、そのままでは
視点が毎回変わってしまい、2 つの条件の絵を並べても見比べられない。
そこで `プロジェクトフォルダ/視点.json` に画角を書き出し、別の条件の画面で
読み込めるようにした。

決めごと:
  ・置き場は**プロジェクト直下**。結果ではなく**入力（どう見るか）**なので
    対象室名・条件名の頭は付けない（`条件表.xlsx` と同じ扱い）
  ・**条件にも受音点にも依らない**（同じ室を見比べるためのもの）
  ・カメラだけでなく**左パネルの数値**（反射回数・表示本数・不透明度など）と
    **開いているタブ**も残す。「画角など」を合わせたいという要望なので、
    見え方を決めるものは一緒に残す
  ・読み込みは**あるものだけ**当てる。欄が違う画面（面の確認 ↔ 結果）で
    読んでも落ちないようにし、当てられた数を画面に知らせる

★カメラを当てたあとは `plotter.camera_set = True` を立てる。
  立てないと、そのあとの `show()` が視点をリセットしてしまう
  （虚音源の `--fit room` で踏んだのと同じ落とし穴）。
"""
import io
import json
import os

FILE_NAME = "視点.json"     # プロジェクト直下（結果ではないので頭は付けない）
FORMAT = 1                  # ファイルの版（読むときに見る）


# ---- 置き場 ----------------------------------------------------------------

def default_path(project):
    """プロジェクトの画角ファイルの場所。条件・受音点で分けない。"""
    return project.path(FILE_NAME)


# ---- カメラ ----------------------------------------------------------------

def camera_state(plotter):
    """いまの画角を、そのまま書き出せる形（辞書）で取り出す。"""
    camera = plotter.camera
    state = {
        "position": [float(v) for v in camera.GetPosition()],
        "focal_point": [float(v) for v in camera.GetFocalPoint()],
        "up": [float(v) for v in camera.GetViewUp()],
        "view_angle": float(camera.GetViewAngle()),
        "parallel_projection": bool(camera.GetParallelProjection()),
        "parallel_scale": float(camera.GetParallelScale()),
        "clipping_range": [float(v) for v in camera.GetClippingRange()],
    }
    try:                    # 記録だけ（当てはしない。ウィンドウの大きさは端末ごと）
        state["window"] = [int(v) for v in plotter.window_size]
    except Exception:
        pass
    return state


def apply_camera(plotter, state):
    """書き出した画角を当てる。**あるものだけ**当てる。"""
    if not state:
        return False
    camera = plotter.camera
    if "position" in state:
        camera.SetPosition(*[float(v) for v in state["position"]])
    if "focal_point" in state:
        camera.SetFocalPoint(*[float(v) for v in state["focal_point"]])
    if "up" in state:
        camera.SetViewUp(*[float(v) for v in state["up"]])
    if "view_angle" in state:
        camera.SetViewAngle(float(state["view_angle"]))
    if "parallel_projection" in state:
        camera.SetParallelProjection(bool(state["parallel_projection"]))
    if "parallel_scale" in state:
        camera.SetParallelScale(float(state["parallel_scale"]))
    if "clipping_range" in state:
        camera.SetClippingRange(*[float(v) for v in state["clipping_range"]])
    # ★これを立てないと、あとの `show()` が視点を戻してしまう
    try:
        plotter.camera_set = True
    except Exception:
        pass
    plotter.render()
    return True


# ---- 左パネルの数値 --------------------------------------------------------

def panel_state(panel):
    """左パネルの数値（反射回数・表示本数など）を見出しの名前ごとに書き出す。"""
    values = []
    for control in (getattr(panel, "controls", None) or []):
        label = control.get("label")
        if not label:
            continue
        try:
            values.append({"label": str(label), "value": float(control["value"])})
        except (KeyError, TypeError, ValueError):
            continue        # 数値でない欄は飛ばす（画角の話には関係ない）
    return values


def apply_panel(panel, saved):
    """書き出した数値を当てる。→ (当てた数, 当てられなかった見出し)

    ★同じ見出しが複数あることがある（タブごとに「反射回数」がある）ので、
    **見出しごとに出てきた順**で突き合わせる。
    """
    buckets = {}
    for control in (getattr(panel, "controls", None) or []):
        buckets.setdefault(str(control.get("label")), []).append(control)
    used, done, missed = {}, 0, []
    for entry in (saved or []):
        label = str(entry.get("label"))
        bucket = buckets.get(label) or []
        index = used.get(label, 0)
        used[label] = index + 1
        if index >= len(bucket):
            missed.append(label)        # この画面には無い欄（別の画面の設定）
            continue
        try:
            bucket[index]["set"](float(entry["value"]))
            done += 1
        except Exception:
            missed.append(label)
    return done, missed


# ---- 読み書き --------------------------------------------------------------

def save(path, plotter, panel=None, tab=None):
    """いまの画角（と表示の設定）を書き出す。→ 書いた場所"""
    if tab is None:
        tab = getattr(panel, "active_group", None) if panel is not None else None
    data = {"format": FORMAT,
            "camera": camera_state(plotter),
            "tab": tab,
            "controls": panel_state(panel) if panel is not None else []}
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with io.open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
    return path


def read(path):
    """書き出した画角を読む。無い・壊れているときは None。"""
    if not path or not os.path.exists(path):
        return None
    try:
        with io.open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (ValueError, OSError):
        return None
    return data if isinstance(data, dict) else None


def load(path, plotter, panel=None, on_tab=None):
    """書き出した画角を当てる。→ (当てられたか, 知らせる文)

    `on_tab` を渡すと、開いていたタブも合わせる（結果画面の音線／音粒子／虚音源）。
    """
    data = read(path)
    if data is None:
        return False, f"画角のファイルがありません（{os.path.basename(path)}）"
    done, missed = 0, []
    if panel is not None:
        # ★数値を先に当てる（反射回数などを変えると描き直しが走るので、
        #   そのあとにカメラを当てて上書きされないようにする）
        done, missed = apply_panel(panel, data.get("controls"))
    tab = data.get("tab")
    if tab and on_tab is not None:
        try:
            on_tab(tab)
        except Exception:
            pass
    ok = apply_camera(plotter, data.get("camera"))
    if panel is not None:
        try:
            panel.relayout(render=True)
        except Exception:
            pass
    message = "画角を読み込みました"
    if done:
        message += f"（設定 {done} 件も合わせました）"
    if missed:
        message += f"／この画面に無い設定 {len(missed)} 件は飛ばしました"
    return (True if ok else False), message


# ---- パネルのボタン --------------------------------------------------------

def add_controls(panel, plotter, path, on_tab=None, heading="画角（見る向き）"):
    """左パネルに「画角を保存」「画角を読込」を足す。

    ★**共通の欄**に置く（`panel.end_group()`）。どのタブを開いていても
      押せるようにするため（タブごとに置くと 3 つ並んで紛らわしい）。
    """
    import view_model_gui as vg

    panel.end_group()

    def do_save():
        try:
            save(path, plotter, panel=panel)
        except Exception as error:
            vg.notice(plotter, f"画角を保存できません"
                               f"（{type(error).__name__}: {error}）", kind="error")
            return
        print(f"[view_camera] 画角を保存しました: {path}")
        vg.notice(plotter, f"画角を保存しました: {os.path.basename(path)}")

    def do_load():
        try:
            ok, message = load(path, plotter, panel=panel, on_tab=on_tab)
        except Exception as error:
            ok = False
            message = f"画角を読み込めません（{type(error).__name__}: {error}）"
        print(f"[view_camera] {message}")
        vg.notice(plotter, message, kind="ok" if ok else "error")

    if heading:
        panel.heading(heading)
    panel.button("画角を保存", do_save, colour="#26402c")
    panel.button("画角を読込", do_load, colour="#26364a")
    return path


# ---- 中身の確認（端末から）--------------------------------------------------

def describe(path):
    """書き出した画角の中身を人が読める形にする（確認用）。"""
    data = read(path)
    if data is None:
        return f"画角のファイルがありません: {path}"
    camera = data.get("camera") or {}
    lines = [f"画角のファイル: {path}",
             f"  視点   {camera.get('position')}",
             f"  注視点 {camera.get('focal_point')}",
             f"  上向き {camera.get('up')}",
             f"  画角   {camera.get('view_angle')} 度"
             f"（平行投影 {camera.get('parallel_projection')}）",
             f"  タブ   {data.get('tab')}"]
    for entry in (data.get("controls") or []):
        lines.append(f"  設定   {entry.get('label')} = {entry.get('value')}")
    return "\n".join(lines)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(
        description="保存した画角の中身を見る（プロジェクトかファイルを指定）")
    parser.add_argument("target", help="プロジェクトのフォルダ、または視点.json")
    args = parser.parse_args(argv)
    path = args.target
    if os.path.isdir(path):
        path = os.path.join(path, FILE_NAME)
    print(describe(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
