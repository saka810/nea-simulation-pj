# -*- coding: utf-8 -*-
"""**ACIS（REGION / 3DSOLID）の DXF を、面の DXF に直す**（TODO B-22）。

★2026-08-28 ユーザー指摘「REGION で読めたはずでは無いですか？」がきっかけ。
読めたことは一度も無く（`docs/DXFデータの作り方.md` に「REGION は本ツールでは
読めない」と書いてある）、B-22 として残っていた宿題。

やり方は TODO B-22 に書いてあったとおり **`accoreconsole` を呼ぶ**。
AutoCAD の中で `EXPLODE` すると、REGION は**輪郭の LINE**に、
3DSOLID は**面の REGION**に分解される。それを拾って輪郭の順に並べ直し、
**閉じた 3D ポリライン**として DXF に書き出す（本ツールが読める形）。

決めごと:
  ・**元のファイルは触らない**。作業は控えに対して行い、`<名前>_faces.dxf` を作る
  ・**レイヤ名は保つ**（吸音材の割り当てに使うので、これが崩れると意味が無い）
  ・★**円弧・スプラインが混じったら黙って捨てず、数えて知らせる**
    （直線だけで囲まれた面しか作れない）
  ・★**穴のある面（輪郭が 2 つ以上）も知らせる**。外周だけ使うと開口が
    塞がってしまうので、勝手に決めない
  ・座標は**元の単位のまま**（`$INSUNITS` も引き継ぐ）
"""
import io
import math
import os
import re
import shutil
import subprocess
import tempfile

# 端点が同じとみなす距離（図面の単位。mm なら 0.01 mm）
JOIN_TOLERANCE = 1.0e-2

# 面が平面とみなせるかのしきい値（同上）
PLANE_TOLERANCE = 1.0

# accoreconsole を探す場所
ACCORE_GLOB = r"C:\Program Files\Autodesk\AutoCAD *\accoreconsole.exe"


def find_accoreconsole(path=None):
    """`accoreconsole.exe` を探す。→ 場所（見つからなければ None）"""
    import glob

    if path:
        return path if os.path.exists(path) else None
    found = sorted(glob.glob(ACCORE_GLOB))
    return found[-1] if found else None      # 新しい版を優先


# ---- AutoCAD 側のスクリプト -------------------------------------------------

# ★`MESHSMOOTH` は accoreconsole では**何も作らない**（画面が無いため）。
#   `CONVTOMESH` は**そもそも無い**。使えるのは `EXPLODE` だけだった。
# ★`EXPLODE` は選択を**空行で閉じる**こと（`(command "_.EXPLODE" o "")`）。
#   閉じないとコマンドが開いたままになり、次の行を選択として食べて止まる。
SCRIPT = """FILEDIA
0
CMDDIA
0
QAFLAGS
1
(setq fp (open "{dump}" "w"))
(defun pt3 (p) (strcat (rtos (car p) 2 6) "," (rtos (cadr p) 2 6) ","
                       (rtos (caddr p) 2 6)))
(defun burst (o / marker c out)
  (setq marker (entlast))
  (command "_.EXPLODE" o "")
  (setq c marker out (list))
  (while (setq c (entnext c)) (setq out (cons c out)))
  (reverse out)
)
(setq ss (ssget "X" (list (cons 0 "REGION,3DSOLID") (cons 410 "Model"))))
(setq queue (list) n 0)
(if ss (repeat (sslength ss)
  (setq queue (cons (ssname ss n) queue) n (1+ n))))
(setq queue (reverse queue) gid 0 guard 0)
(while (and queue (< guard 20000))
  (setq guard (1+ guard))
  (setq o (car queue) queue (cdr queue))
  (setq info (entget o))
  (setq kind (cdr (assoc 0 info)) lay (cdr (assoc 8 info)))
  (setq kids (burst o))
  (if (= kind "REGION")
    (progn
      (setq gid (1+ gid))
      (foreach c kids
        (setq d (entget c))
        (setq ty (cdr (assoc 0 d)))
        (cond
          ((= ty "LINE")
            (write-line (strcat "L," (itoa gid) "," lay ","
                                (pt3 (cdr (assoc 10 d))) ","
                                (pt3 (cdr (assoc 11 d)))) fp))
          ((member ty (list "REGION" "3DSOLID"))
            (setq queue (cons c queue)))
          (t (write-line (strcat "X," (itoa gid) "," lay "," ty) fp))
        )
      )
    )
    (foreach c kids
      (setq ty (cdr (assoc 0 (entget c))))
      (if (member ty (list "REGION" "3DSOLID"))
        (setq queue (cons c queue))
        (write-line (strcat "X,0," lay "," ty) fp)
      )
    )
  )
)
(write-line (strcat "END," (itoa gid)) fp)
(close fp)
(princ (strcat "|GROUPS=" (itoa gid) "|"))
QUIT
Y
"""


def _run_autocad(source, dump, accore, timeout=1800, verbose=True):
    """accoreconsole で ACIS を分解し、辺を `dump` に書き出させる。"""
    folder = os.path.dirname(dump)
    script = os.path.join(folder, "burst.scr")
    with io.open(script, "w", encoding="ascii", newline="\r\n") as handle:
        handle.write(SCRIPT.format(dump=dump.replace("\\", "/")))
    if verbose:
        print(f"[面に分解] AutoCAD を呼びます（{os.path.basename(accore)}）…")
    result = subprocess.run([accore, "/i", source, "/s", script, "/l", "en-US"],
                            capture_output=True, timeout=timeout)
    if not os.path.exists(dump):
        text = result.stdout.decode("utf-16-le", "replace")[-1200:]
        raise RuntimeError("AutoCAD が辺を書き出しませんでした。"
                           "最後の出力:\n" + text)
    return dump


# ---- 輪郭を組み立てる -------------------------------------------------------

def read_dump(path):
    """AutoCAD が書いた辺の一覧を読む。→ (グループ→辺, 落ちたもの)"""
    groups, dropped = {}, []
    with io.open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = [p.strip() for p in line.strip().split(",")]
            if not parts or parts[0] == "END":
                continue
            if parts[0] == "X":
                dropped.append((parts[2] if len(parts) > 2 else "",
                                parts[3] if len(parts) > 3 else "?"))
                continue
            if parts[0] != "L" or len(parts) < 9:
                continue
            gid, layer = int(parts[1]), parts[2]
            first = tuple(float(v) for v in parts[3:6])
            second = tuple(float(v) for v in parts[6:9])
            groups.setdefault(gid, {"layer": layer, "edges": []})
            groups[gid]["edges"].append((first, second))
    return groups, dropped


def _key(point, tolerance):
    return tuple(int(round(v / tolerance)) for v in point)


def loops_from_edges(edges, tolerance=JOIN_TOLERANCE):
    """辺の集まりを**閉じた輪**に組み直す。→ [輪（点の並び）, …]

    ★穴のある面は輪が 2 つ以上になる。**外周だけ使って黙って塞がない**ため、
      全部返して呼び側で数える。
    """
    remaining = list(edges)
    loops = []
    while remaining:
        chain = list(remaining.pop(0))
        while True:
            tail = chain[-1]
            for index, (first, second) in enumerate(remaining):
                if _key(first, tolerance) == _key(tail, tolerance):
                    chain.append(second)
                    remaining.pop(index)
                    break
                if _key(second, tolerance) == _key(tail, tolerance):
                    chain.append(first)
                    remaining.pop(index)
                    break
            else:
                break                       # つながる辺が無い
            if _key(chain[-1], tolerance) == _key(chain[0], tolerance):
                chain.pop()                 # 閉じた（最後は始点と同じ）
                loops.append(chain)
                chain = None
                break
        if chain is not None:
            loops.append(chain)             # 閉じなかった（開いた鎖）
    return loops


def plane_error(points):
    """輪が平面に乗っているか。→ 面からの最大のずれ"""
    if len(points) < 4:
        return 0.0
    centre = [sum(p[k] for p in points) / len(points) for k in range(3)]
    # 最も面積の大きい三角形から法線を作る（細長い並びで崩れないように）
    best, normal = 0.0, None
    first = points[0]
    for a in range(1, len(points) - 1):
        u = [points[a][k] - first[k] for k in range(3)]
        v = [points[a + 1][k] - first[k] for k in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1],
                 u[2] * v[0] - u[0] * v[2],
                 u[0] * v[1] - u[1] * v[0]]
        size = math.sqrt(sum(c * c for c in cross))
        if size > best:
            best, normal = size, [c / size for c in cross]
    if normal is None:
        return 0.0
    return max(abs(sum(normal[k] * (p[k] - centre[k]) for k in range(3)))
               for p in points)


# ---- DXF を書く -------------------------------------------------------------

def _header(insunits):
    return ("0\nSECTION\n2\nHEADER\n"
            f"9\n$INSUNITS\n70\n{int(insunits)}\n"
            "0\nENDSEC\n")


def _tables(layers):
    out = ["0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n70\n%d\n" % len(layers)]
    for name in layers:
        out.append("0\nLAYER\n2\n%s\n70\n0\n62\n7\n6\nCONTINUOUS\n" % name)
    out.append("0\nENDTAB\n0\nENDSEC\n")
    return "".join(out)


def write_faces_dxf(path, polygons, insunits=4):
    """輪を**閉じた 3D ポリライン**として書く（本ツールが読める形）。"""
    layers = sorted({layer for layer, _points in polygons}) or ["0"]
    with io.open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(_header(insunits))
        handle.write(_tables(layers))
        handle.write("0\nSECTION\n2\nENTITIES\n")
        for layer, points in polygons:
            # 70 = 1(閉じている) + 8(3D ポリライン)
            handle.write("0\nPOLYLINE\n8\n%s\n66\n1\n70\n9\n"
                         "10\n0.0\n20\n0.0\n30\n0.0\n" % layer)
            for point in points:
                handle.write("0\nVERTEX\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n%.6f\n"
                             "70\n32\n" % (layer, point[0], point[1], point[2]))
            handle.write("0\nSEQEND\n8\n%s\n" % layer)
        handle.write("0\nENDSEC\n0\nEOF\n")
    return path


def _insunits_of(dxf_path):
    """元の DXF の `$INSUNITS`（無ければ 4 ＝ mm とみなす）。"""
    try:
        with io.open(dxf_path, encoding="utf-8", errors="replace") as handle:
            text = handle.read(200000)
    except OSError:
        return 4
    match = re.search(r"\$INSUNITS\s*\n\s*70\s*\n\s*(\d+)", text)
    return int(match.group(1)) if match else 4


# ---- 入り口 -----------------------------------------------------------------

def convert(dxf_path, out_path=None, accore=None, verbose=True, keep=False):
    """ACIS の DXF → 面の DXF。→ 書き出した場所

    ★**元のファイルは触らない**（控えを作って、そちらを AutoCAD に開かせる）。
    """
    dxf_path = os.path.abspath(dxf_path)
    if out_path is None:
        base, ext = os.path.splitext(dxf_path)
        out_path = base + "_faces" + ext
    accore = find_accoreconsole(accore)
    if accore is None:
        raise RuntimeError(
            "accoreconsole.exe が見つかりません（AutoCAD が要ります）。"
            "CAD 側で面（3DFACE か閉じたポリライン）に分解して"
            "書き出してください")

    folder = tempfile.mkdtemp(prefix="geosim_faces_")
    try:
        # ★控えに対して作業する（元のファイルは開かせない）
        work = os.path.join(folder, "work.dxf")
        shutil.copy2(dxf_path, work)
        dump = os.path.join(folder, "edges.txt")
        _run_autocad(work, dump, accore, verbose=verbose)
        groups, dropped = read_dump(dump)

        polygons, holes, open_chains, small = [], 0, 0, 0
        twisted = []
        for gid in sorted(groups):
            layer = groups[gid]["layer"]
            loops = loops_from_edges(groups[gid]["edges"])
            closed = [loop for loop in loops if len(loop) >= 3]
            if not closed:
                open_chains += 1
                continue
            if len(closed) > 1:
                holes += 1
            # 輪が複数なら**いちばん長いもの**を面にする（残りは知らせる）
            closed.sort(key=len, reverse=True)
            points = closed[0]
            if len(points) < 3:
                small += 1
                continue
            error = plane_error(points)
            if error > PLANE_TOLERANCE:
                twisted.append((layer, error))
            polygons.append((layer, points))

        if verbose:
            print(f"[面に分解] ACIS {len(groups)} 面 → 輪郭 {len(polygons)} 枚")
            if holes:
                print(f"[面に分解] ★輪郭が 2 つ以上の面が {holes} 枚あります"
                      "（穴あき）。**外周だけ**を面にしたので、"
                      "開口が塞がっているかもしれません")
            if open_chains:
                print(f"[面に分解] ★閉じなかった輪郭が {open_chains} 枚"
                      "（辺が足りない）。この面は落としました")
            if twisted:
                worst = max(v for _l, v in twisted)
                print(f"[面に分解] ★平面に乗っていない輪郭が {len(twisted)} 枚"
                      f"（最大 {worst:.1f}）。読み込み時に三角形へ割られます")
            if dropped:
                kinds = {}
                for _layer, kind in dropped:
                    kinds[kind] = kinds.get(kind, 0) + 1
                print(f"[面に分解] ★直線でない辺・面がありました: {kinds}"
                      "（円弧やスプラインは面にできません）")

        write_faces_dxf(out_path, polygons, insunits=_insunits_of(dxf_path))
        if verbose:
            print(f"[面に分解] 書き出しました: {out_path}")
        return out_path
    finally:
        if keep:
            print(f"[面に分解] 作業フォルダを残しました: {folder}")
        else:
            shutil.rmtree(folder, ignore_errors=True)


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="ACIS（REGION / 3DSOLID）の DXF を面の DXF に直す")
    parser.add_argument("dxf", help="元の DXF")
    parser.add_argument("--out", default=None, help="書き出し先")
    parser.add_argument("--accore", default=None,
                        help="accoreconsole.exe の場所")
    parser.add_argument("--keep", action="store_true",
                        help="作業フォルダを残す（中身を確かめたいとき）")
    args = parser.parse_args(argv)

    path = convert(args.dxf, args.out, accore=args.accore, keep=args.keep)

    # そのまま読めるか確かめる
    import read_dxffile as rd
    model = rd.read_model(path, verbose=True)
    print(f"[面に分解] 読み込み確認: 三角形 {len(model.mesh)} 枚 / "
          f"レイヤ {len(model.layer_areas)} 種")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
