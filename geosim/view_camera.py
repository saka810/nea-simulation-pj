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
FORMAT = 2                  # ファイルの版（1 = 1 ファイル 1 画角 / 2 = 名前つき）
DEFAULT_NAME = "既定"        # 名前を言われなかったときの 1 本


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


# ---- 読み書き（**名前つきで何本でも持てる**）----------------------------
#
# ★2026-08-24 ユーザー要望
# > 角度は色々保存したいので、保存ボタン押したらダイアログボックス出して
# > 名前決めれるように、読込も同じように。
#
# ファイルは 1 つ（`視点.json`）のまま、中に**名前 → 画角**を並べる。
# 昔の形（1 ファイル 1 画角）も読める（`既定` という名前の 1 本として扱う）。

def _document(path):
    """ファイルの中身を、いまの形（名前 → 画角）にそろえて返す。"""
    data = read(path) or {}
    views = data.get("views")
    if isinstance(views, dict):
        return {"format": FORMAT, "views": dict(views)}
    if "camera" in data:
        # 昔の形（2026-08-24 の最初の版）。1 本ぶんとして拾う
        one = {key: data[key] for key in ("camera", "tab", "controls")
               if key in data}
        return {"format": FORMAT, "views": {DEFAULT_NAME: one}}
    return {"format": FORMAT, "views": {}}


def _write(path, document):
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with io.open(path, "w", encoding="utf-8") as handle:
        json.dump(document, handle, ensure_ascii=False, indent=2)
    return path


def view_names(path):
    """保存されている画角の名前（保存した順）。"""
    return list(_document(path)["views"])


def next_name(path, stem="視点"):
    """まだ使っていない名前（`視点1`、`視点2`…）。ダイアログの初期値に使う。"""
    names = set(view_names(path))
    index = 1
    while f"{stem}{index}" in names:
        index += 1
    return f"{stem}{index}"


def save(path, plotter, panel=None, tab=None, name=DEFAULT_NAME):
    """いまの画角（と表示の設定）を `name` という名前で書き出す。

    同じ名前があれば**その 1 本だけ**を書き替える（他の画角は残す）。
    """
    if tab is None:
        tab = getattr(panel, "active_group", None) if panel is not None else None
    document = _document(path)
    document["views"][str(name)] = {
        "camera": camera_state(plotter),
        "tab": tab,
        "controls": panel_state(panel) if panel is not None else []}
    return _write(path, document)


def remove(path, name):
    """名前を指定して 1 本消す。→ 消せたか"""
    document = _document(path)
    if str(name) not in document["views"]:
        return False
    del document["views"][str(name)]
    _write(path, document)
    return True


def read(path):
    """書き出したファイルをそのまま読む。無い・壊れているときは None。"""
    if not path or not os.path.exists(path):
        return None
    try:
        with io.open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (ValueError, OSError):
        return None
    return data if isinstance(data, dict) else None


def _pick(views, name):
    """読み込む 1 本を決める。名前を言われなければ順に当てにいく。"""
    if name is not None:
        return str(name) if str(name) in views else None
    if len(views) == 1:
        return next(iter(views))
    if DEFAULT_NAME in views:
        return DEFAULT_NAME
    return next(iter(views), None)


def load(path, plotter, panel=None, on_tab=None, name=None):
    """保存した画角を当てる。→ (当てられたか, 知らせる文)

    `on_tab` を渡すと、開いていたタブも合わせる（結果画面の音線／音粒子／虚音源）。
    """
    views = _document(path)["views"]
    if not views:
        return False, f"保存された画角がありません（{os.path.basename(path)}）"
    chosen = _pick(views, name)
    if chosen is None:
        return False, f"「{name}」という画角はありません"
    state = views[chosen] or {}
    done, missed = 0, []
    if panel is not None:
        # ★数値を先に当てる（反射回数などを変えると描き直しが走るので、
        #   そのあとにカメラを当てて上書きされないようにする）
        done, missed = apply_panel(panel, state.get("controls"))
    tab = state.get("tab")
    if tab and on_tab is not None:
        try:
            on_tab(tab)
        except Exception:
            pass
    ok = apply_camera(plotter, state.get("camera"))
    if panel is not None:
        try:
            panel.relayout(render=True)
        except Exception:
            pass
    message = f"画角「{chosen}」を読み込みました"
    if done:
        message += f"（設定 {done} 件も合わせました）"
    if missed:
        message += f"／この画面に無い設定 {len(missed)} 件は飛ばしました"
    return (True if ok else False), message


# ---- 名前を決めるダイアログ（2026-08-24 ユーザー要望）----------------------
#
# ★開いている 3D の画面はそのままに、tkinter の小さな窓を出して閉じる
#   （数値のまとめ入力 `ControlPanel.open_value_input` と同じやり方）。
#   出しっぱなしにはしないので、画面の取り合いは起きない。

def ask_name(names, default=None, title="画角に名前を付けて保存"):
    """名前を決める窓。→ 名前（キャンセルなら None）"""
    import tkinter as tk
    from tkinter import ttk

    chosen = {"name": None}
    root = tk.Tk()
    root.title(title)
    root.attributes("-topmost", True)
    frame = ttk.Frame(root, padding=12)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="名前を付けて保存します（同じ名前にすると上書きします）",
              foreground="#666").pack(anchor="w", pady=(0, 8))
    text = tk.StringVar(value=default or DEFAULT_NAME)
    entry = ttk.Entry(frame, textvariable=text, width=32)
    entry.pack(fill="x")
    entry.selection_range(0, "end")
    entry.focus_set()

    if names:
        ttk.Label(frame, text="いま入っている画角（選ぶと名前が入ります）",
                  foreground="#666").pack(anchor="w", pady=(10, 2))
        listbox = tk.Listbox(frame, height=min(6, len(names)),
                             exportselection=False)
        for name in names:
            listbox.insert("end", name)
        listbox.pack(fill="x")

        def picked(_event=None):
            selection = listbox.curselection()
            if selection:
                text.set(listbox.get(selection[0]))

        listbox.bind("<<ListboxSelect>>", picked)

    def apply():
        name = text.get().strip()
        if name:
            chosen["name"] = name
        root.destroy()

    buttons = ttk.Frame(frame)
    buttons.pack(fill="x", pady=(12, 0))
    ttk.Button(buttons, text="キャンセル",
               command=root.destroy).pack(side="right", padx=4)
    ttk.Button(buttons, text="保存", command=apply).pack(side="right", padx=4)
    root.bind("<Return>", lambda _event: apply())
    root.bind("<Escape>", lambda _event: root.destroy())
    root.mainloop()
    return chosen["name"]


def ask_view(names, title="画角を読み込む"):
    """どの画角を読むか選ぶ窓。→ (名前, したいこと)。

    したいことは `"load"`（読み込む）か `"delete"`（消す）。
    キャンセルなら `(None, None)`。
    """
    import tkinter as tk
    from tkinter import ttk

    chosen = {"name": None, "action": None}
    root = tk.Tk()
    root.title(title)
    root.attributes("-topmost", True)
    frame = ttk.Frame(root, padding=12)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="読み込む画角を選んでください",
              foreground="#666").pack(anchor="w", pady=(0, 8))
    listbox = tk.Listbox(frame, height=min(10, max(3, len(names))),
                         exportselection=False)
    for name in names:
        listbox.insert("end", name)
    listbox.pack(fill="both", expand=True)
    listbox.selection_set(0)
    listbox.focus_set()

    def decide(action):
        selection = listbox.curselection()
        if selection:
            chosen["name"] = listbox.get(selection[0])
            chosen["action"] = action
        root.destroy()

    buttons = ttk.Frame(frame)
    buttons.pack(fill="x", pady=(12, 0))
    ttk.Button(buttons, text="キャンセル",
               command=root.destroy).pack(side="right", padx=4)
    ttk.Button(buttons, text="読み込む",
               command=lambda: decide("load")).pack(side="right", padx=4)
    # 名前を付けて何本も貯まるので、ここで消せるようにしておく
    ttk.Button(buttons, text="削除",
               command=lambda: decide("delete")).pack(side="left", padx=4)
    listbox.bind("<Double-Button-1>", lambda _event: decide("load"))
    root.bind("<Return>", lambda _event: decide("load"))
    root.bind("<Escape>", lambda _event: root.destroy())
    root.mainloop()
    return chosen["name"], chosen["action"]


# ---- パネルのボタン --------------------------------------------------------

def add_controls(panel, plotter, path, on_tab=None, heading="画角（見る向き）"):
    """左パネルに「画角を保存」「画角を読込」を足す。

    ★**共通の欄**に置く（`panel.end_group()`）。どのタブを開いていても
      押せるようにするため（タブごとに置くと 3 つ並んで紛らわしい）。
    """
    import view_model_gui as vg

    panel.end_group()

    # 直前に使った名前（続けて保存するときの初期値。上書きが楽になる）
    last = {"name": None}

    def do_save():
        try:
            names = view_names(path)
            name = ask_name(names, default=last["name"] or next_name(path))
            if not name:
                return                      # キャンセル
            save(path, plotter, panel=panel, name=name)
            last["name"] = name
        except Exception as error:
            vg.notice(plotter, f"画角を保存できません"
                               f"（{type(error).__name__}: {error}）", kind="error")
            return
        print(f"[view_camera] 画角「{name}」を保存しました: {path}")
        vg.notice(plotter, f"画角「{name}」を保存しました")

    def do_load():
        try:
            names = view_names(path)
            if not names:
                vg.notice(plotter, "保存された画角がありません"
                                   "（先に「画角を保存」を押してください）",
                          kind="error")
                return
            name, action = ask_view(names)
            if not name:
                return                      # キャンセル
            if action == "delete":
                remove(path, name)
                print(f"[view_camera] 画角「{name}」を削除しました")
                vg.notice(plotter, f"画角「{name}」を削除しました")
                return
            ok, message = load(path, plotter, panel=panel, on_tab=on_tab,
                               name=name)
            last["name"] = name
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
    if read(path) is None:
        return f"画角のファイルがありません: {path}"
    views = _document(path)["views"]
    lines = [f"画角のファイル: {path}（{len(views)} 本）"]
    for name, state in views.items():
        camera = (state or {}).get("camera") or {}
        lines += [f"■ {name}",
                  f"  視点   {camera.get('position')}",
                  f"  注視点 {camera.get('focal_point')}",
                  f"  上向き {camera.get('up')}",
                  f"  画角   {camera.get('view_angle')} 度"
                  f"（平行投影 {camera.get('parallel_projection')}）",
                  f"  タブ   {(state or {}).get('tab')}"]
        for entry in ((state or {}).get("controls") or []):
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
