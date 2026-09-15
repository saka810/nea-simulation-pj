"""**測定点（音源・受音点）の並び**を決める。

2026-09-15 ユーザー要望（不具合報告 A-3 の相談）:

> 基本は作成順で良いですが、GUI 上で修正できるようにできないかな？

`read_dxffile` は DXF に**出てきた順**に番号を振るので、`結果/rec1/` の 1 が
CAD 側のラベル（R1）と一致するとは限らない。実案件（階段教室・2026-09-15）で
**`rec1` が R4 になっていた**。気づかないと結果を取り違える。

CAD を描き直させるのが筋だが、図面をもらってから気づくことが多いので、
**プロジェクト側で並べ替えられる**ようにした。

    プロジェクトフォルダ/測定点順.json

★**入力**なので対象室名・条件名の頭は付けない（`視点.json` `断面.json` と同じ）。
★**座標は持たない**（並べ替えだけ）。DXF を直したら順番も見直すことになるので、
  `normals.json` `materials.json` と同じく**点の数を控えて、食い違ったら使わない**。

中身:

    {
      "dxf": "260914_階段教室.dxf",   控え（食い違っても使う。知らせるだけ）
      "source_count": 2,             控えた点の数。**合わなければ使わない**
      "receiver_count": 5,
      "sources": [0, 1],             DXF に出てきた順の番号を、使いたい順に並べる
      "receivers": [3, 0, 1, 4, 2]   ＝ rec1 は DXF の 4 点目、rec2 は 1 点目…
    }

使い方:

    python point_order.py <プロジェクト>          いまの並びを見る
    python point_order.py <プロジェクト> --edit    並べ替えの画面を開く
"""

import json
import os

ORDER_FILE = "測定点順.json"


def path(project):
    return project.path(ORDER_FILE)


def load(project, verbose=True):
    """`測定点順.json` を読む。無ければ None。"""
    target = path(project)
    if not os.path.exists(target):
        return None
    try:
        with open(target, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError) as error:
        if verbose:
            print(f"[測定点順] {ORDER_FILE} を読めませんでした: {error}")
        return None


def save(project, sources=None, receivers=None, model=None):
    """並びを保存する。→ 書いたパス

    `sources` / `receivers` は **DXF に出てきた順の番号**を並べたリスト。
    点の数も一緒に控える（DXF が変わったら使わないため）。
    """
    data = {
        "dxf": os.path.basename(project.dxf or ""),
        "source_count": len(sources) if sources is not None else
                        (len(model.source_points) if model is not None else 0),
        "receiver_count": len(receivers) if receivers is not None else
                          (len(model.receiver_points) if model is not None else 0),
        "sources": list(sources or []),
        "receivers": list(receivers or []),
    }
    with open(path(project), "w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
    return path(project)


def _valid(order, count, kind, verbose=True):
    """並びが使えるか。**数が合い、0〜count-1 が 1 回ずつ**でなければ使わない。"""
    if not order:
        return None
    order = [int(v) for v in order]
    if sorted(order) != list(range(count)):
        if verbose:
            print(f"[測定点順] {ORDER_FILE} の{kind}の並びが今のモデルと"
                  f"合わないので**使いません**（控え {len(order)} 点 / "
                  f"いま {count} 点）。並べ替えをやり直してください")
        return None
    return order


def orders_for(project, model, verbose=True):
    """このモデルに当てられる並びを返す → (音源の並び, 受音点の並び)。無ければ None。"""
    data = load(project, verbose=verbose)
    if not data:
        return None, None
    sources = _valid(data.get("sources"), len(model.source_points or []),
                     "音源", verbose=verbose)
    receivers = _valid(data.get("receivers"), len(model.receiver_points or []),
                       "受音点", verbose=verbose)
    return sources, receivers


def apply(project, model, verbose=True):
    """モデルの点を `測定点順.json` の並びに直す（その場で入れ替える）。→ model

    ★**読んだ直後に 1 回だけ通す**。ここを通ったあとは、どのモジュールから見ても
    「1 番目の受音点」は利用者が決めた 1 番目になる。
    """
    if model is None or project is None:
        return model
    sources, receivers = orders_for(project, model, verbose=verbose)
    if sources:
        model.source_points = [model.source_points[k] for k in sources]
    if receivers:
        model.receiver_points = [model.receiver_points[k] for k in receivers]
        names = getattr(model, "receiver_layer_names", None)
        if names and len(names) == len(receivers):
            model.receiver_layer_names = [names[k] for k in receivers]
    if verbose and (sources or receivers):
        print(f"[測定点順] {ORDER_FILE} の並びを当てました"
              f"（音源 {'○' if sources else '—'} / 受音点 {'○' if receivers else '—'}）")
    return model


def describe(project, model):
    """いまの並びを人が読む形で返す（確認用）。"""
    lines = []
    for kind, points, names in (
            ("音源", model.source_points or [], None),
            ("受音点", model.receiver_points or [],
             getattr(model, "receiver_layer_names", None))):
        prefix = "src" if kind == "音源" else "rec"
        for index, point in enumerate(points):
            layer = ""
            if names and index < len(names) and names[index]:
                layer = f"  レイヤ {names[index]}"
            lines.append(f"{prefix}{index + 1}  "
                         f"({point[0]:.3f}, {point[1]:.3f}, {point[2]:.3f}) m{layer}")
    return "\n".join(lines)


# ------------------------------------------------------------------------------
# 並べ替えの画面（tkinter。設定画面から開く）
# ------------------------------------------------------------------------------

def edit(project, model=None, parent=None):
    """並べ替えの窓を開く。→ 保存したら True

    ★**一覧の上から順に `rec1` `rec2` …** になる。上へ・下へで動かすだけの
    素直な作りにしてある（利用者がやりたいのは「R4 を 4 番目へ持っていく」だけ）。
    """
    import tkinter as tk
    from tkinter import messagebox, ttk

    if model is None:
        model = _read_model(project)
    if model is None:
        if parent is not None:
            messagebox.showerror("測定点の並び", "DXF を読めませんでした")
        return False

    # いまの並び（保存済みがあればそれを当てた状態）から始める
    apply(project, model, verbose=False)
    sources = list(model.source_points or [])
    receivers = list(model.receiver_points or [])
    names = list(getattr(model, "receiver_layer_names", None) or
                 ["" for _ in receivers])
    saved, _ = orders_for(project, model, verbose=False)
    # 画面では「いまの並び」を扱い、保存するときに DXF の順番へ翻訳する。
    # そのため**元の番号**を一緒に持ち歩く
    source_ids, receiver_ids = _current_ids(project, model)

    window = tk.Toplevel(parent) if parent is not None else tk.Tk()
    window.title("測定点の並び（上から src1 / rec1 …）")
    window.geometry("560x520")

    ttk.Label(window, justify="left", foreground="#333", text=(
        "上から順に src1 / rec1 … になります（結果フォルダ `結果/recN/` の N）。\n"
        "★CAD の作成順が既定です。ラベル（R1…）と合わないときだけ直してください。")
    ).pack(anchor="w", padx=12, pady=(10, 4))

    boxes = {}
    for kind, points, ids in (("音源", sources, source_ids),
                              ("受音点", receivers, receiver_ids)):
        if not points:
            continue
        frame = ttk.LabelFrame(window, text=kind)
        frame.pack(fill="both", expand=True, padx=12, pady=6)
        listbox = tk.Listbox(frame, height=min(8, max(3, len(points))),
                             exportselection=False)
        listbox.pack(side="left", fill="both", expand=True, padx=(8, 4), pady=8)
        rows = []
        for order, point in enumerate(points):
            layer = ""
            if kind == "受音点" and order < len(names) and names[order]:
                layer = f"   [{names[order]}]"
            rows.append((ids[order],
                         f"({point[0]:.3f}, {point[1]:.3f}, {point[2]:.3f}) m{layer}"))
        state = {"rows": rows, "points": list(points),
                 "names": list(names) if kind == "受音点" else []}
        boxes[kind] = (listbox, state)

        def redraw(listbox=listbox, state=state, kind=kind):
            keep = listbox.curselection()
            listbox.delete(0, tk.END)
            prefix = "src" if kind == "音源" else "rec"
            for order, (origin, text) in enumerate(state["rows"]):
                listbox.insert(tk.END,
                               f"{prefix}{order + 1}   {text}   "
                               f"（DXF の {origin + 1} 番目）")
            if keep:
                listbox.selection_set(keep[0])

        def move(step, listbox=listbox, state=state, redraw=redraw):
            picked = listbox.curselection()
            if not picked:
                return
            index = picked[0]
            target = index + step
            if not 0 <= target < len(state["rows"]):
                return
            rows = state["rows"]
            rows[index], rows[target] = rows[target], rows[index]
            redraw()
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(target)
            listbox.see(target)

        side = ttk.Frame(frame)
        side.pack(side="left", fill="y", padx=(0, 8), pady=8)
        ttk.Button(side, text="▲ 上へ", command=lambda m=move: m(-1)).pack(fill="x")
        ttk.Button(side, text="▼ 下へ", command=lambda m=move: m(1)).pack(fill="x",
                                                                        pady=(4, 0))
        ttk.Button(side, text="CAD の順に戻す",
                   command=lambda s=state, r=redraw: (_reset(s), r())).pack(
                       fill="x", pady=(12, 0))
        redraw()

    result = {"saved": False}

    def on_save():
        order = {}
        for kind, (_listbox, state) in boxes.items():
            order[kind] = [origin for origin, _text in state["rows"]]
        save(project,
             sources=order.get("音源", list(range(len(sources)))),
             receivers=order.get("受音点", list(range(len(receivers)))))
        result["saved"] = True
        messagebox.showinfo("測定点の並び",
                            f"保存しました:\n{path(project)}\n\n"
                            f"★次の計算から効きます。**計算済みの結果フォルダ"
                            f"（結果/recN/）の中身は並べ替えません**")
        window.destroy()

    row = ttk.Frame(window)
    row.pack(fill="x", padx=12, pady=(0, 12))
    ttk.Button(row, text="保存", command=on_save).pack(side="right")
    ttk.Button(row, text="取り消し", command=window.destroy).pack(side="right",
                                                              padx=(0, 8))
    window.transient(parent) if parent is not None else None
    window.grab_set()
    window.wait_window()
    return result["saved"]


def _reset(state):
    """CAD に出てきた順（元の番号の昇順）に戻す。"""
    state["rows"].sort(key=lambda row: row[0])


def _current_ids(project, model):
    """いまの並びの各点が「DXF の何番目か」を返す → (音源, 受音点)。"""
    sources, receivers = orders_for(project, model, verbose=False)
    if sources is None:
        sources = list(range(len(model.source_points or [])))
    if receivers is None:
        receivers = list(range(len(model.receiver_points or [])))
    return sources, receivers


def _read_model(project):
    import read_dxffile as rd
    try:
        return rd.read_model(project.dxf_path, unit=project.unit,
                             band_number=project.band_number, verbose=False)
    except Exception as error:
        print(f"[測定点順] DXF を読めませんでした: {type(error).__name__}: {error}")
        return None


def main():
    import argparse

    import project as pj

    p = argparse.ArgumentParser(description="測定点（音源・受音点）の並びを見る／直す")
    p.add_argument("folder", help="プロジェクトフォルダ")
    p.add_argument("--edit", action="store_true", help="並べ替えの画面を開く")
    p.add_argument("--reset", action="store_true",
                   help="CAD の作成順に戻す（`測定点順.json` を消す）")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if a.reset:
        target = path(project)
        if os.path.exists(target):
            os.remove(target)
            print(f"[測定点順] 消しました: {target}（CAD の作成順に戻ります）")
        else:
            print("[測定点順] もともと CAD の作成順です")
        return

    model = _read_model(project)
    if model is None:
        raise SystemExit("DXF を読めませんでした")
    if a.edit:
        edit(project, model)
        model = _read_model(project)
    apply(project, model)
    print(describe(project, model))


if __name__ == "__main__":
    main()
