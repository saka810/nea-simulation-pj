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
    python point_order.py <プロジェクト> --match-previous
                                                  前回の計算と同じ並びに戻す

★★**DXF を作り直すと並びが変わる**（2026-09-20。不具合報告 ⑯）。
実案件で**画層名を 1 つ変えただけ**の DXF に差し替えたら、三角形は 779 枚とも
完全に同じなのに `POINT` の出てくる順が変わり、**src1 が S2 に、rec3 が R5 に**
なっていた。黙って番号が振り直されるので、結果を見るまで気づけない。
しかも**経路キャッシュの指紋が合わなくなる**ので音線追跡からやり直しになる
（実案件で 50 分回してから気づいて中止）。

そこで `check_against_previous()` が、**前回の `結果/<室>_測定点.csv` と
座標を突き合わせて並びの変化を知らせる**（`run_project` がモデルを読んだ直後に
1 回だけ呼ぶ）。戻し方まで出すので、`--match-previous` で直せる。
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


# ------------------------------------------------------------------------------
# ★前回の計算と突き合わせる（2026-09-20。不具合報告 ⑯）
# ------------------------------------------------------------------------------

KINDS = (("音源", "sources", "src"), ("受音点", "receivers", "rec"))


def previous_points(project):
    """前回書いた `結果/<室>_測定点.csv` から座標を読む → {"音源": [...], "受音点": [...]}。

    ★**結果**なので、条件を変えても対象室が同じなら同じファイルを指す
    （`SOURCE_SHARED_RESULTS` ＋ `ROOM_SCOPED_RESULTS`）。
    """
    import csv

    try:
        target = project.existing_result_path("points")
    except Exception:
        target = None
    if not target or not os.path.exists(target):
        return {}
    found = {"音源": [], "受音点": []}
    try:
        with open(target, encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                kind = (row.get("区分") or "").strip()
                if kind not in found:
                    continue
                try:
                    found[kind].append([float(row["X_m"]), float(row["Y_m"]),
                                        float(row["Z_m"])])
                except (KeyError, TypeError, ValueError):
                    return {}           # 読めない列があるなら当てにしない
    except OSError:
        return {}
    return found


def _match_order(before, now, tolerance=1.0e-6):
    """前回の点が**いまの何番目**かを返す。集合として一致しなければ None。"""
    import numpy as np

    if not before or len(before) != len(now):
        return None
    now = [np.asarray(p, dtype=float).reshape(3) for p in now]
    order, used = [], set()
    for point in before:
        point = np.asarray(point, dtype=float).reshape(3)
        hit = [i for i, q in enumerate(now)
               if i not in used and np.allclose(q, point, atol=tolerance)]
        if not hit:
            return None                 # 位置そのものが変わっている（別のモデル）
        order.append(hit[0])
        used.add(hit[0])
    return order


def check_against_previous(project, model, applied=True, verbose=True):
    """前回の計算と**並びが変わっていないか**を見る → 直すための並び（無ければ {}）。

    ★★**黙って番号が振り直されるのを防ぐ**（不具合報告 ⑯）。座標の集合は同じで
    順番だけが違うときに、何がどう入れ替わったかと戻し方を知らせる。
    **位置そのものが変わっているときは何も言わない**（点を動かしたなら当然なので）。

    applied : bool
        ★渡されたモデルに `測定点順.json` が**当たっているか**。
        `run_project._ordered` は当てたあとに呼ぶので True、
        `point_order.main` は DXF に出てきた順のまま呼ぶので False。
        ★**戻す並びは必ず「DXF に出てきた順の番号」**（`測定点順.json` の約束）
        なので、既にある並びと**合成して**返す。
    """
    if project is None or model is None:
        return {}
    before = previous_points(project)
    if not before:
        return {}
    kept = dict(zip(("sources", "receivers"),
                    orders_for(project, model, verbose=False)))
    fix, notes = {}, []
    for label, key, prefix in KINDS:
        now = (model.source_points if key == "sources"
               else model.receiver_points) or []
        saved = kept.get(key)
        if applied or not saved:
            effective = list(now)
        else:
            effective = [now[k] for k in saved]
        order = _match_order(before.get(label), effective)
        if order is None or order == list(range(len(order))):
            continue
        # 「前回の i 番目 = いまの order[i] 番目」を DXF 順の番号に直す
        fix[key] = [saved[o] for o in order] if saved else list(order)
        for index, target in enumerate(order):
            if index != target:
                notes.append(f"  前回の {prefix}{index + 1} は"
                             f"**いまの {target + 1} 番目**です")
    if fix and verbose:
        print("[測定点順] ★★測定点の並びが前回の計算と変わっています"
              "（位置は同じで順番だけ違います）")
        for line in notes:
            print(f"[測定点順] {line}")
        print("[測定点順]   ★このまま回すと `結果/src1/` `結果/rec1/` に"
              "別の点の結果が入り、**経路の使い回しも効きません**")
        print(f"[測定点順]   前回と同じ並びに戻すには: "
              f'python point_order.py "{project.folder}" --match-previous')
        print(f"[測定点順]   画面で決めるなら --edit / このままでよければ無視して"
              f"構いません（`{ORDER_FILE}` は書き換えません）")
    return fix



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
    p.add_argument("--match-previous", action="store_true",
                   help="前回の計算（結果/<室>_測定点.csv）と同じ並びに戻す")
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
    if getattr(a, "match_previous", False):
        # ★**並べ替える前の（DXF に出てきた順の）モデル**で突き合わせる
        fix = check_against_previous(project, model, applied=False, verbose=True)
        if not fix:
            print("[測定点順] 前回と同じ並びです（直すところはありません）")
        else:
            save(project, sources=fix.get("sources"),
                 receivers=fix.get("receivers"), model=model)
            print(f"[測定点順] 前回と同じ並びにしました: {path(project)}")
            model = _read_model(project)
    else:
        check_against_previous(project, model, applied=False, verbose=True)
    apply(project, model)
    print(describe(project, model))


if __name__ == "__main__":
    main()
