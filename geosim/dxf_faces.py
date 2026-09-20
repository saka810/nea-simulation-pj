# -*- coding: utf-8 -*-
"""**ACIS（REGION / 3DSOLID）の DXF を、面の DXF に直す**（TODO B-22）。

★2026-08-28 ユーザー指摘「REGION で読めたはずでは無いですか？」がきっかけ。
読めたことは一度も無く（`docs/DXFデータの作り方.md` に「REGION は本ツールでは
読めない」と書いてある）、B-22 として残っていた宿題。

やり方は TODO B-22 に書いてあったとおり **`accoreconsole` を呼ぶ**。
AutoCAD の中で `EXPLODE` すると、REGION は**輪郭の LINE**に、
3DSOLID は**面の REGION**に分解される。それを拾って輪郭の順に並べ直し、
**閉じた 3D ポリライン**として DXF に書き出す（本ツールが読める形）。
穴のある面だけは三角形に割って **3DFACE** で書く（下記）。

決めごと:
  ・**元のファイルは触らない**。作業は控えに対して行い、`<名前>_faces.dxf` を作る
  ・**レイヤ名は保つ**（吸音材の割り当てに使うので、これが崩れると意味が無い）
  ・★**円弧・スプラインが混じったら黙って捨てず、数えて知らせる**
    （直線だけで囲まれた面しか作れない）
  ・★★**穴のある面（輪郭が 2 つ以上）は穴を開けたまま三角形に割る**
    （2026-09-09。それまでは外周だけを面にしていて**開口が塞がっていた**）。
    穴を橋でつないで（キーホール法）耳刈りし、**3DFACE** で書く。
    ★面積が「外周 − 穴」と合わなければ穴を開けず、外周だけに戻して知らせる
  ・★★**元の図面の `POINT` を引き継ぐ**（2026-09-11。不具合報告 ⑧）。
    音源・受音点は `src` / `rec` 画層の POINT で渡す決めなのに面しか書いて
    いなかったので、**変換を通すたびに音源も受音点も消えていた**
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

# ★★**`accoreconsole` はシステムのコードページで書き出す**（日本語 Windows なら CP932）。
#   2026-09-06 の不具合報告 ⑤ で実案件を踏んだ：UTF-8 だけで読んでいたので
#   **日本語の画層名が全部 U+FFFD に置き換わり、そのまま出力 DXF に焼き付いていた**。
#   画層名は吸音材の割り当てに使うので、崩れると条件表が引けない。
#   さらに悪いのは**化け方が同じ画層が 1 つにまとめられる**ことで、
#   実案件（階段教室）では 22 画層が 20 に減った
#   （`PHP_階段裏`＋`PHP_階段下`、`開口_2F小`＋`開口_2F大` が統合された）。
#   ★`errors="replace"` は**最後の手段**にする（黙って壊すのを避けるため）。
DUMP_ENCODINGS = ("cp932", "utf-8")

# 元の DXF を読むときの順。こちらは UTF-8 が普通
# （`$DWGCODEPAGE ANSI_932` でも中身は UTF-8 のことが多い）。
# `read_dxffile` と同じ並びにしてある
DXF_ENCODINGS = ("utf-8", "cp932")


def read_text(path, encodings=DUMP_ENCODINGS, label=""):
    """テキストを、順に試して読めた文字コードで読む。

    どれでも読めなければ**最後の手段**として `errors="replace"` に落とし、
    ★黙って壊さずに理由を告げる（画層名が崩れると吸音材が引けなくなるため）。
    """
    with io.open(path, "rb") as handle:
        raw = handle.read()
    for encoding in encodings:
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    print(f"[面に分解] ★{label or os.path.basename(path)} の文字コードが分かりません"
          f"（{' / '.join(encodings)} のどれでもありません）。"
          f"読めない文字を置き換えて続けます。**画層名が崩れるかもしれません**")
    return raw.decode(encodings[0], errors="replace")


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
    # ★**CP932 を先に試す**（accoreconsole はシステムのコードページで書く）。
    #   UTF-8 決め打ちだと日本語の画層名が壊れる（不具合報告 ⑤。`DUMP_ENCODINGS`）
    text = read_text(path, DUMP_ENCODINGS, label="AutoCAD が書いた辺の一覧")
    for line in text.splitlines():
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


def write_faces_dxf(path, polygons, insunits=4, triangles=None, points=None):
    """輪を**閉じた 3D ポリライン**として書く（本ツールが読める形）。

    `triangles` は [(レイヤ, (点, 点, 点)), …]。**穴のある面**はポリライン 1 本では
    表せないので、三角形に割って **3DFACE** で書く（2026-09-09。`read_dxffile` は
    3DFACE も読めるし、読み込み側の**同一平面パッチ**が 1 枚にまとめ直す）。

    `points` は [(レイヤ, (x, y, z)), …]。元の図面の `POINT`（音源・受音点）を
    そのまま書き戻す（2026-09-11。不具合報告 ⑧）。
    """
    triangles = list(triangles or [])
    points = list(points or [])
    layers = sorted({layer for layer, _points in polygons}
                    | {layer for layer, _corners in triangles}
                    | {layer for layer, _xyz in points}) or ["0"]
    with io.open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(_header(insunits))
        handle.write(_tables(layers))
        handle.write("0\nSECTION\n2\nENTITIES\n")
        for layer, xyz in points:
            handle.write("0\nPOINT\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n%.6f\n"
                         % (layer, xyz[0], xyz[1], xyz[2]))
        for layer, corners in triangles:
            # 3DFACE は 4 隅を持つ。三角形なので 4 点目は 3 点目と同じにする
            first, second, third = corners
            handle.write("0\n3DFACE\n8\n%s\n" % layer)
            for index, point in enumerate((first, second, third, third)):
                handle.write("%d\n%.6f\n%d\n%.6f\n%d\n%.6f\n"
                             % (10 + index, point[0], 20 + index, point[1],
                                30 + index, point[2]))
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


# ---- 元の図面の POINT（音源・受音点）----------------------------------------
#
# ★★2026-09-11 ユーザー要望「POINT も引き継ぎたいですね」（不具合報告 ⑧）。
#   音源・受音点は `src` / `rec` 画層の POINT で渡す決め（`read_dxffile`）なのに、
#   ここは面しか書いていなかったので**変換を通すたびに両方とも消えていた**。
#   実案件（階段教室）では利用者が `Src` 画層に置いた音源 2 点が落ち、
#   「音源が DXF にありません」と言われる状態になっていた。
#
#   ★**accoreconsole は通さない。**読めないのは ACIS（REGION / 3DSOLID）だけで、
#     元の DXF は**テキストとして普通に読める**。POINT のグループコードは
#     そのまま取れるので、元のファイルから直接拾って書き足す。

def _tags(text):
    """DXF を (グループコード, 値) の並びにする。

    ASCII DXF は**コードと値が 1 行ずつ交互**に並ぶだけなので、2 行ずつ読む。
    数字でない行が来たら（壊れたファイル）その組は飛ばす。
    """
    lines = text.replace("\r\n", "\n").split("\n")
    for index in range(0, len(lines) - 1, 2):
        code = lines[index].strip()
        if not code.lstrip("-").isdigit():
            continue
        yield int(code), lines[index + 1].strip()


def read_points(dxf_path):
    """元の DXF の `POINT` を拾う。→ [(画層, (x, y, z)), …]

    ★**モデル空間（`ENTITIES`）だけ**を見る。`BLOCKS` の中の POINT は
    ブロック定義の座標系なので、そのまま置くと位置が合わない
    （`INSERT` の挿入点・尺度・回転を当てないと世界座標にならない）。

    ★**画層で絞らない**。`src` / `rec` かどうかは `read_dxffile` が
    画層名で選ぶので、ここで決め打ちすると命名の揺れ（`Src` / `音源` など）に
    追随できなくなる。余分に入っていても読み込み側が無視する。
    """
    try:
        text = read_text(dxf_path, DXF_ENCODINGS, label="元の DXF")
    except OSError:
        return []

    points = []
    section, kind = None, None
    layer, xyz = "0", [None, None, None]

    def keep():
        if kind == "POINT" and None not in xyz:
            points.append((layer, (xyz[0], xyz[1], xyz[2])))

    for code, value in _tags(text):
        if code == 0:
            keep()
            if value == "SECTION":
                section = "(名前待ち)"
            elif value == "ENDSEC":
                section = None
            kind = value
            layer, xyz = "0", [None, None, None]
            continue
        if section == "(名前待ち)" and code == 2:
            section = value
            continue
        if section != "ENTITIES" or kind != "POINT":
            continue
        if code == 8:
            layer = value
        elif code in (10, 20, 30):
            try:
                xyz[code // 10 - 1] = float(value)
            except ValueError:
                pass
    keep()
    return sort_points(points)


def sort_points(points):
    """`POINT` を**画層 → 座標**の順に並べ直す（2026-09-20。不具合報告 ⑯）。

    ★★**作り直しても番号が動かないようにする**ため。CAD は画層名を変えただけでも
    `ENTITIES` の並びを変えることがあり、実案件では**三角形が 779 枚とも同じなのに
    `src1` が S2 に、`rec3` が R5 になった**。`read_dxffile` は出てきた順に番号を
    振るので、黙って `結果/src1/` の中身が別の点にすり替わる
    （さらに経路キャッシュの指紋も合わなくなり、音線追跡からやり直しになる）。

    ★**元の図面の並びには依存しない**決め方にしておけば、何度変換しても同じ番号になる。
    画層名（`rec1` `rec2` …）でまとめてから座標で並べるので、
    測線ごとに作った点も測線の順に並ぶ。
    """
    def key(item):
        layer, xyz = item
        return (str(layer).lower(), round(xyz[0], 6), round(xyz[1], 6),
                round(xyz[2], 6))

    return sorted(points, key=key)


def _insunits_of(dxf_path):
    """元の DXF の `$INSUNITS`（無ければ 4 ＝ mm とみなす）。"""
    try:
        text = read_text(dxf_path, DXF_ENCODINGS, label="元の DXF")[:200000]
    except OSError:
        return 4
    match = re.search(r"\$INSUNITS\s*\n\s*70\s*\n\s*(\d+)", text)
    return int(match.group(1)) if match else 4


# ---- 穴のある面（輪が 2 つ以上）--------------------------------------------
#
# ★★2026-09-09 ユーザー指摘「元々ドーナッツ状（中が空いている）の床が、
#   一面の床になってしまっています」。実案件（階段教室）の `床_1F` が
#   **外周 367.48 m² ＋ 穴 102.27 m²** で、本来 265.22 m² のドーナツなのに
#   外周だけを面にしていたため**2 層吹き抜けが幻の床で塞がっていた**。
#
#   それまでは「輪はいちばん長いものだけ使う」と割り切って**数えて知らせる**
#   だけだった（書き出す形が「閉じた 3D ポリライン 1 本」で、1 本のポリラインでは
#   穴を表せないため）。ここで**穴を橋でつないで三角形に割り、3DFACE で書く**
#   ようにした（`read_dxffile` は 3DFACE も読める）。
#
#   ★三角形で書いても、読み込み側の**同一平面パッチ**が 1 枚にまとめ直すので、
#     面の確認画面では元どおり「1 枚の床」として見える。
#   ★★**面積が「外周 − 穴」と合わなければ穴を開けない**（外周だけに戻して知らせる）。
#     黙って変な形の床を作るほうが危ないため。

# 橋（キーホール）の線が他の辺と交わっていないとみなす許容（図面の単位）
BRIDGE_TOLERANCE = 1.0e-6

# 三角形に割ったあとの面積が合っているかのしきい値（相対）
AREA_TOLERANCE = 1.0e-6


def loop_normal(points):
    """輪の法線（ニューウェル法）。→ 単位ベクトル（決まらなければ None）"""
    normal = [0.0, 0.0, 0.0]
    count = len(points)
    for index in range(count):
        current, following = points[index], points[(index + 1) % count]
        normal[0] += (current[1] - following[1]) * (current[2] + following[2])
        normal[1] += (current[2] - following[2]) * (current[0] + following[0])
        normal[2] += (current[0] - following[0]) * (current[1] + following[1])
    size = math.sqrt(sum(c * c for c in normal))
    if size < 1.0e-12:
        return None
    return [c / size for c in normal]


def _basis(normal):
    """法線に垂直な 2 本の軸（面の中で 2 次元に落とすため）。"""
    # いちばん寝ている軸を種にすると、外積が縮退しない
    seed = [0.0, 0.0, 0.0]
    seed[min(range(3), key=lambda k: abs(normal[k]))] = 1.0
    u = [seed[1] * normal[2] - seed[2] * normal[1],
         seed[2] * normal[0] - seed[0] * normal[2],
         seed[0] * normal[1] - seed[1] * normal[0]]
    size = math.sqrt(sum(c * c for c in u))
    u = [c / size for c in u]
    v = [normal[1] * u[2] - normal[2] * u[1],
         normal[2] * u[0] - normal[0] * u[2],
         normal[0] * u[1] - normal[1] * u[0]]
    return u, v


def _flatten(points, u, v):
    """面の中の 2 次元座標にする。"""
    return [(sum(p[k] * u[k] for k in range(3)),
             sum(p[k] * v[k] for k in range(3))) for p in points]


def signed_area_2d(points):
    """2 次元の輪の符号付き面積（反時計回りなら正）。"""
    total = 0.0
    for index in range(len(points)):
        x1, y1 = points[index]
        x2, y2 = points[(index + 1) % len(points)]
        total += x1 * y2 - x2 * y1
    return total / 2.0


def _same2d(a, b, tolerance=JOIN_TOLERANCE):
    return abs(a[0] - b[0]) <= tolerance and abs(a[1] - b[1]) <= tolerance


def _crosses(a, b, c, d):
    """線分 ab と cd が**内部で**交わるか（端点で触れるだけなら False）。"""
    def side(p, q, r):
        return ((q[0] - p[0]) * (r[1] - p[1])
                - (q[1] - p[1]) * (r[0] - p[0]))

    d1, d2 = side(a, b, c), side(a, b, d)
    d3, d4 = side(c, d, a), side(c, d, b)
    straddle_first = ((d1 > BRIDGE_TOLERANCE and d2 < -BRIDGE_TOLERANCE)
                      or (d1 < -BRIDGE_TOLERANCE and d2 > BRIDGE_TOLERANCE))
    straddle_second = ((d3 > BRIDGE_TOLERANCE and d4 < -BRIDGE_TOLERANCE)
                       or (d3 < -BRIDGE_TOLERANCE and d4 > BRIDGE_TOLERANCE))
    return straddle_first and straddle_second


def _inside_2d(point, ring):
    """2 次元の点が輪の内側か（交差数。境界は数に入れない）。"""
    inside = False
    for index in range(len(ring)):
        x1, y1 = ring[index]
        x2, y2 = ring[(index + 1) % len(ring)]
        if (y1 > point[1]) != (y2 > point[1]):
            crossing = x1 + (point[1] - y1) * (x2 - x1) / (y2 - y1)
            if point[0] < crossing:
                inside = not inside
    return inside


def bridge_holes(outer, holes):
    """穴を**橋（キーホール）でつないで 1 つの輪にする**。→ 2 次元の点の並び

    引数・戻り値とも**面の中の 2 次元座標**。外周は反時計回り、穴は時計回りに
    揃えてから `outer[:i+1] + hole[j:] + hole[:j+1] + outer[i:]` の形でつなぐ
    （教科書どおりのキーホール法。橋の両端は 2 度使われる）。

    ★橋の選び方は**総当たりでいちばん短いもの**。実案件は外周 13 点・穴 4 点なので
      総当たりでも一瞬で、視線判定の取りこぼしが無い。
      橋が他の辺と交わるもの、リングの外へ出るものは弾く。
    """
    ring = list(outer)
    if signed_area_2d(ring) < 0.0:
        ring.reverse()
    loops = []
    for hole in holes:
        loop = list(hole)
        if signed_area_2d(loop) > 0.0:       # 穴は外周と逆向きに揃える
            loop.reverse()
        loops.append(loop)

    for loop in loops:
        edges = [(ring[k], ring[(k + 1) % len(ring)]) for k in range(len(ring))]
        for other in loops:
            edges += [(other[k], other[(k + 1) % len(other)])
                      for k in range(len(other))]
        best = None
        for i, outer_point in enumerate(ring):
            for j, hole_point in enumerate(loop):
                if _same2d(outer_point, hole_point):
                    continue
                length = ((outer_point[0] - hole_point[0]) ** 2
                          + (outer_point[1] - hole_point[1]) ** 2)
                if best is not None and length >= best[0]:
                    continue
                if any(_crosses(outer_point, hole_point, a, b) for a, b in edges):
                    continue
                # 橋の中点がリングの中（外周の内側・どの穴の外側）にあること
                middle = ((outer_point[0] + hole_point[0]) / 2.0,
                          (outer_point[1] + hole_point[1]) / 2.0)
                if not _inside_2d(middle, ring):
                    continue
                if any(_inside_2d(middle, other) for other in loops):
                    continue
                best = (length, i, j)
        if best is None:
            return None                     # つなげない（外周だけに戻して知らせる）
        _length, i, j = best
        ring = ring[:i + 1] + loop[j:] + loop[:j + 1] + ring[i:]
    return ring


def _convex_2d(a, b, c):
    return ((b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0])) > BRIDGE_TOLERANCE


def _in_triangle_2d(point, a, b, c):
    """点が三角形 abc の中か。★**角と重なる点は「中に無い」とみなす**。

    キーホール法は橋の両端を**2 度使う**（同じ座標の点が 2 つある）ので、
    素直に判定すると常に「中にある」ことになり、耳が 1 つも見つからなくなる。
    """
    if _same2d(point, a) or _same2d(point, b) or _same2d(point, c):
        return False
    d1 = (b[0] - a[0]) * (point[1] - a[1]) - (b[1] - a[1]) * (point[0] - a[0])
    d2 = (c[0] - b[0]) * (point[1] - b[1]) - (c[1] - b[1]) * (point[0] - b[0])
    d3 = (a[0] - c[0]) * (point[1] - c[1]) - (a[1] - c[1]) * (point[0] - c[0])
    return ((d1 >= -BRIDGE_TOLERANCE and d2 >= -BRIDGE_TOLERANCE
             and d3 >= -BRIDGE_TOLERANCE)
            or (d1 <= BRIDGE_TOLERANCE and d2 <= BRIDGE_TOLERANCE
                and d3 <= BRIDGE_TOLERANCE))


def ear_clip_2d(ring):
    """2 次元の輪を耳刈り法で三角形に割る。→ [(添字, 添字, 添字), …] | None

    反時計回りの輪を前提にする（`bridge_holes` がそう揃えて返す）。
    """
    order = list(range(len(ring)))
    triangles = []
    guard = 0
    while len(order) > 3 and guard <= len(ring) * len(ring) + 10:
        guard += 1
        cut = None
        for m in range(len(order)):
            i, j, k = order[m - 1], order[m], order[(m + 1) % len(order)]
            if not _convex_2d(ring[i], ring[j], ring[k]):
                continue
            if any(_in_triangle_2d(ring[q], ring[i], ring[j], ring[k])
                   for q in order if q not in (i, j, k)):
                continue
            cut = (m, i, j, k)
            break
        if cut is None:
            return None
        m, i, j, k = cut
        triangles.append((i, j, k))
        order.pop(m)
    if len(order) != 3:
        return None
    triangles.append(tuple(order))
    return triangles


def triangles_with_holes(loops):
    """穴のある面を三角形に割る。→ [(点, 点, 点), …] | None

    `loops` は**いちばん長いものを外周**とみなす（`convert` がそう並べて渡す）。
    ★面積が「外周 − 穴」と合わなければ `None` を返す
      （黙って変な形の面を作らないため）。
    """
    outer, holes = loops[0], loops[1:]
    normal = loop_normal(outer)
    if normal is None:
        return None
    u, v = _basis(normal)
    flat_outer = _flatten(outer, u, v)
    flat_holes = [_flatten(hole, u, v) for hole in holes]

    ring = bridge_holes(flat_outer, flat_holes)
    if ring is None:
        return None
    order = ear_clip_2d(ring)
    if order is None:
        return None

    # 2 次元に落とす前の座標に戻すための対応表
    lookup = {}
    for flat, points in [(flat_outer, outer)] + list(zip(flat_holes, holes)):
        for flat_point, point in zip(flat, points):
            lookup[(round(flat_point[0], 6), round(flat_point[1], 6))] = point

    want = abs(signed_area_2d(flat_outer))
    for hole in flat_holes:
        want -= abs(signed_area_2d(hole))

    triangles, got = [], 0.0
    for i, j, k in order:
        corners = []
        for index in (i, j, k):
            key = (round(ring[index][0], 6), round(ring[index][1], 6))
            if key not in lookup:
                return None
            corners.append(lookup[key])
        got += abs(signed_area_2d([ring[i], ring[j], ring[k]]))
        triangles.append(tuple(corners))
    if want <= 0.0 or abs(got - want) > AREA_TOLERANCE * max(want, 1.0):
        return None
    return triangles


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
        triangles = []          # 穴のある面（3DFACE で書く）
        holed, unfilled = 0, []
        for gid in sorted(groups):
            layer = groups[gid]["layer"]
            loops = loops_from_edges(groups[gid]["edges"])
            closed = [loop for loop in loops if len(loop) >= 3]
            if not closed:
                open_chains += 1
                continue
            # ★輪はいちばん長いものを外周とみなす（穴より外周のほうが点が多い）
            closed.sort(key=len, reverse=True)
            if len(closed) > 1:
                holes += 1
                # ★★**穴を開けたまま三角形に割る**（2026-09-09。それまでは
                #   外周だけを面にしていて、ドーナツ状の床が塞がっていた）
                parts = triangles_with_holes(closed)
                if parts:
                    holed += 1
                    for part in parts:
                        triangles.append((layer, part))
                    continue
                # 割れなかったら外周だけに戻す。★黙って落とさず知らせる
                unfilled.append(layer)
            points = closed[0]
            if len(points) < 3:
                small += 1
                continue
            error = plane_error(points)
            if error > PLANE_TOLERANCE:
                twisted.append((layer, error))
            polygons.append((layer, points))

        if verbose:
            print(f"[面に分解] ACIS {len(groups)} 面 → 輪郭 {len(polygons)} 枚"
                  + (f" ＋ 穴のある面 {holed} 枚（三角形 {len(triangles)} 枚）"
                     if holed else ""))
            if holed:
                print(f"[面に分解] ★穴のある面 {holed} 枚は**穴を開けたまま**"
                      f"三角形に割りました（開口は塞いでいません）")
            if unfilled:
                print(f"[面に分解] ★★穴を開けられなかった面が {len(unfilled)} 枚"
                      f"あります（{' / '.join(sorted(set(unfilled)))}）。"
                      f"**外周だけ**を面にしたので、開口が塞がっています。"
                      f"CAD 側で面を分けて書き出してください")
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

        # ★元の図面の POINT（音源・受音点）を引き継ぐ（2026-09-11。不具合報告 ⑧）
        points = read_points(dxf_path)
        if verbose:
            if points:
                counts = {}
                for layer, _xyz in points:
                    counts[layer] = counts.get(layer, 0) + 1
                detail = " / ".join(f"{name} {n}" for name, n
                                    in sorted(counts.items()))
                print(f"[面に分解] 点（POINT）を {len(points)} 個"
                      f"引き継ぎました（画層: {detail}）")
                # ★**並べ直したことを言う**（番号が変わりうるので黙って直さない）
                print("[面に分解] ★点は画層 → 座標の順に並べ直しました"
                      "（作り直しても `src1` `rec1` の番号が動かないように）。"
                      "前に計算したことがある室なら、結果を見る前に "
                      "`python point_order.py <プロジェクト>` で並びを確かめてください")
            else:
                print("[面に分解] 元の図面に点（POINT）はありませんでした。"
                      "★音源・受音点は `src` / `rec` 画層の POINT で渡します")

        write_faces_dxf(out_path, polygons, insunits=_insunits_of(dxf_path),
                        triangles=triangles, points=points)
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
