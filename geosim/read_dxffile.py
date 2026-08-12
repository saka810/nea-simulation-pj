import csv
import os

import numpy as np

from mesh import Mesh


# DXF ファイルの読み込み（元コード backtrace.f90 132〜283 行）
#
# ------------------------------------------------------------------------------
# 【DXF の構造】ASCII DXF は「グループコード」と「値」が 1 行ずつ交互に並ぶだけの形式。
#   0        ← グループコード（0 = エンティティの開始）
#   POLYLINE ← 値
#   8        ← グループコード（8 = レイヤ名）
#   1        ← 値
#   ...
# つまり 2 行ずつ読んで (code, value) のタプル列にすれば、あとは code を見るだけで済む。
#
# 【室形状として読むもの】
#   POLYLINE (70 & 64) … ポリフェイスメッシュ。直後に VERTEX が並び SEQEND で終わる
#     VERTEX (70 = 192) … 頂点の座標。10/20/30 が x/y/z
#     VERTEX (70 = 128) … 面レコード。71/72/73/74 が頂点番号（1 始まり、符号は辺の可視性）
#   3DFACE            … 3〜4 頂点の面を単体で表すエンティティ（10/20/30 〜 13/23/33）
#
# 【音源・受音点】
#   POINT … レイヤ名が src / rec のものを音源・受音点として拾う
#
# 【単位】
#   HEADER の $INSUNITS が単位コードを持つ（4 = mm、6 = m など）。
#   これを見て m に換算する。CAD 側で設定されていない（0 = unitless）場合は
#   read(unit="mm") のように明示的に指定する。
# ------------------------------------------------------------------------------

# $INSUNITS のコード → 1 単位が何メートルか
INSUNITS_TO_METER = {
    0: None,          # unitless（判定不能。呼び出し側で指定が必要）
    1: 0.0254,        # インチ
    2: 0.3048,        # フィート
    3: 1609.344,      # マイル
    4: 0.001,         # ミリメートル
    5: 0.01,          # センチメートル
    6: 1.0,           # メートル
    7: 1000.0,        # キロメートル
    8: 2.54e-8,       # マイクロインチ
    9: 2.54e-5,       # ミル
    10: 0.9144,       # ヤード
    11: 1.0e-10,      # オングストローム
    12: 1.0e-9,       # ナノメートル
    13: 1.0e-6,       # ミクロン
    14: 0.1,          # デシメートル
    15: 10.0,         # デカメートル
    16: 100.0,        # ヘクトメートル
    17: 1.0e9,        # ギガメートル
}

# 文字列で単位を指定するとき用
UNIT_ALIAS = {
    "mm": 0.001, "millimeter": 0.001, "ミリ": 0.001,
    "cm": 0.01, "centimeter": 0.01, "センチ": 0.01,
    "m": 1.0, "meter": 1.0, "メートル": 1.0,
    "km": 1000.0,
    "inch": 0.0254, "in": 0.0254,
    "feet": 0.3048, "ft": 0.3048,
}

# VERTEX の 70 フラグ
VERTEX_FLAG_POLYFACE_MESH = 64      # 3D ポリゴンメッシュ頂点
VERTEX_FLAG_FACE_RECORD = 128       # ポリフェイスメッシュ頂点

# POLYLINE の 70 フラグ
POLYLINE_FLAG_POLYFACE_MESH = 64

DEGENERATE_EPS = 1.0e-12

DEFAULT_SOURCE_LAYERS = ("src", "source", "音源")
DEFAULT_RECEIVER_LAYERS = ("rec", "receiver", "受音点")


class DxfModel:
    """DXF から読み取った内容の入れ物。

    属性:
        mesh            list[Mesh]  三角形メッシュ（座標は m に換算済み、法線は単位ベクトル）
        source_points   list[(3,)]  src レイヤの POINT（m）
        receiver_points list[(3,)]  rec レイヤの POINT（m）
        unit_scale      float       1 CAD 単位 = 何 m か
        unit_source     str         単位をどこから決めたか（'$INSUNITS' / '指定'）
        layer_counts    dict        レイヤ名 → 三角形の枚数
        skipped         dict        読み飛ばした要素の内訳
        extents         (2,3)       バウンディングボックス（m）
    """

    def __init__(self):
        self.mesh = []
        self.source_points = []
        self.receiver_points = []
        self.unit_scale = 1.0
        self.unit_source = ""
        self.layer_counts = {}
        self.skipped = {"縮退面": 0, "非対応エンティティ": 0, "面レコード不正": 0}
        self.extents = None

    def summary(self):
        lines = [
            f"三角形メッシュ: {len(self.mesh)} 枚",
            f"単位: 1 CAD単位 = {self.unit_scale} m（{self.unit_source}）",
            f"レイヤ別の枚数: {self.layer_counts}",
            f"音源: {[np.round(p, 4).tolist() for p in self.source_points]}",
            f"受音点: {[np.round(p, 4).tolist() for p in self.receiver_points]}",
        ]
        if self.extents is not None:
            size = self.extents[1] - self.extents[0]
            lines.append(f"寸法: {np.round(size, 4).tolist()} m "
                         f"(min={np.round(self.extents[0], 4).tolist()})")
        if any(self.skipped.values()):
            lines.append(f"読み飛ばし: {self.skipped}")
        return "\n".join(lines)


# ------------------------------------------------------------------------------
# DXF の低レベル読み込み
# ------------------------------------------------------------------------------

def _read_tags(file_name):
    """ASCII DXF を (グループコード, 値) のタプル列にする。

    DXF のエンコーディングは版によって異なる（R2007 以降は UTF-8、それ以前は CP932 など）。
    レイヤ名に日本語の吸音材名を使う可能性があるので、順に試してデコードする。
    """
    last_error = None
    for encoding in ("utf-8", "cp932", "latin-1"):
        try:
            with open(file_name, encoding=encoding) as f:
                lines = [line.strip() for line in f]
            break
        except UnicodeDecodeError as e:
            last_error = e
    else:
        raise last_error

    tags = []
    for i in range(0, len(lines) - 1, 2):
        code = lines[i]
        if not code.lstrip("-").isdigit():
            # コード行でないものが来たら壊れたファイル。そこで打ち切る
            break
        tags.append((int(code), lines[i + 1]))
    return tags


def _header_value(tags, name):
    """HEADER 変数（$INSUNITS など）の値を返す。見つからなければ None。"""
    for i, (code, value) in enumerate(tags):
        if code == 9 and value == name and i + 1 < len(tags):
            return tags[i + 1][1]
    return None


def _split_entities(tags):
    """ENTITIES セクションをエンティティ単位に切り分ける。

    戻り値: list[(type, list[(code, value)])]
    """
    start = None
    for i, (code, value) in enumerate(tags):
        if code == 2 and value == "ENTITIES":
            start = i + 1
            break
    if start is None:
        return []

    entities = []
    current = None
    for code, value in tags[start:]:
        if code == 0:
            if value == "ENDSEC":
                break
            current = (value, [])
            entities.append(current)
        elif current is not None:
            current[1].append((code, value))
    return entities


def _tag(entity_tags, code, default=None):
    for c, v in entity_tags:
        if c == code:
            return v
    return default


def _float(entity_tags, code, default=0.0):
    v = _tag(entity_tags, code)
    return default if v is None else float(v)


def _int(entity_tags, code, default=0):
    v = _tag(entity_tags, code)
    return default if v is None else int(float(v))


# ------------------------------------------------------------------------------
# 単位
# ------------------------------------------------------------------------------

def resolve_unit_scale(tags, unit=None):
    """1 CAD 単位が何メートルかを決める。

    unit を明示指定すればそれを使う（'mm' / 'm' などの文字列、または数値で m/単位）。
    指定がなければ HEADER の $INSUNITS を見る。
    どちらも決まらなければ 1.0（m 扱い）にして警告する。
    """
    if unit is not None:
        if isinstance(unit, str):
            key = unit.strip().lower()
            if key not in UNIT_ALIAS:
                raise ValueError(f"未知の単位指定: {unit!r}（使えるのは {sorted(UNIT_ALIAS)}）")
            return UNIT_ALIAS[key], f"指定 unit={unit!r}"
        return float(unit), f"指定 unit={unit}"

    raw = _header_value(tags, "$INSUNITS")
    if raw is not None:
        code = int(float(raw))
        scale = INSUNITS_TO_METER.get(code)
        if scale is not None:
            return scale, f"$INSUNITS={code}"
        print(f"[read_dxffile] 警告: $INSUNITS={code} は単位不定です。"
              f"unit='mm' のように明示指定してください。とりあえず m として扱います。")
        return 1.0, f"$INSUNITS={code}（不定→m と仮定）"

    print("[read_dxffile] 警告: $INSUNITS が見つかりません。m として扱います。")
    return 1.0, "$INSUNITS なし（m と仮定）"


# ------------------------------------------------------------------------------
# 吸音率
# ------------------------------------------------------------------------------

def read_absorption_csv(file_name, band_number=6):
    """吸音率 CSV を読んで {材料名: ndarray(band_number,)} を返す。

    想定フォーマット（1 行目がヘッダでもよい。数値に変換できなければヘッダとみなす）:
        材料名, a1, a2, a3, a4, a5, a6
        コンクリート, 0.02, 0.02, 0.03, 0.03, 0.04, 0.05
    """
    table = {}
    with open(file_name, encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f):
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            name = row[0].strip()
            try:
                values = [float(x) for x in row[1:1 + band_number]]
            except ValueError:
                continue  # ヘッダ行
            if len(values) < band_number:
                print(f"[read_dxffile] 警告: 吸音率の列が {len(values)} 個しかありません "
                      f"(材料 '{name}')。最後の値で埋めます。")
                values += [values[-1]] * (band_number - len(values)) if values else [0.0] * band_number
            table[name] = np.array(values[:band_number], dtype=float)
    return table


def _resolve_absorption(layer, absorption_table, default_absorption, band_number, unresolved):
    if absorption_table and layer in absorption_table:
        return absorption_table[layer]
    unresolved.add(layer)
    if default_absorption is None:
        return np.full(band_number, 0.1)
    return np.broadcast_to(np.asarray(default_absorption, dtype=float),
                           (band_number,)).copy()


# ------------------------------------------------------------------------------
# 法線
# ------------------------------------------------------------------------------

def face_normal(vertex_1, vertex_2, vertex_3):
    """三角形の単位法線を外積で求める。縮退面（面積ゼロ）なら None。

    n = (x2 - x1) × (x3 - x1) を正規化したもの。
    外積は 2 辺が張る平面に垂直なベクトルを返すので、それが面の法線になる。
    向きは頂点の並び順（右ねじ）で決まる。
    """
    v12 = np.asarray(vertex_2, dtype=float) - np.asarray(vertex_1, dtype=float)
    v13 = np.asarray(vertex_3, dtype=float) - np.asarray(vertex_1, dtype=float)
    n = np.cross(v12, v13)
    length = np.linalg.norm(n, ord=2)
    if length < DEGENERATE_EPS:
        return None
    return n / length


def signed_volume(triangles):
    """閉じたメッシュの符号付き体積 V = (1/6)Σ x1·(x2×x3)。

    頂点の巻き順が一貫している閉曲面なら、
      V > 0 … 法線（v12×v13）が外向き
      V < 0 … 法線が内向き
    凸でなくても成り立つのが利点（重心との比較より汎用）。
    """
    total = 0.0
    for x1, x2, x3 in triangles:
        total += float(np.dot(x1, np.cross(x2, x3)))
    return total / 6.0


# ------------------------------------------------------------------------------
# 本体
# ------------------------------------------------------------------------------

def read_model(file_name, unit=None, absorption_table=None, default_absorption=None,
               orient_normals="auto", band_number=6,
               source_layers=DEFAULT_SOURCE_LAYERS,
               receiver_layers=DEFAULT_RECEIVER_LAYERS,
               verbose=True):
    """DXF を読んで DxfModel を返す。

    引数:
        unit              None なら $INSUNITS から自動判定。'mm' / 'm' などで明示指定も可
        absorption_table  {レイヤ名: ndarray} または吸音率 CSV のパス
        default_absorption テーブルに無いレイヤに使う値（None なら 0.1 で警告）
        orient_normals    'auto'    … 符号付き体積で内向きに揃える（既定）
                          'cad'     … CAD の巻き順をそのまま信じる
                          'inward'  … 重心方向へ強制（凸形状のみ）
                          'outward' … 重心の逆へ強制（凸形状のみ）
        source_layers     音源として扱う POINT のレイヤ名
        receiver_layers   受音点として扱う POINT のレイヤ名
    """
    if isinstance(absorption_table, str):
        absorption_table = read_absorption_csv(absorption_table, band_number)

    tags = _read_tags(file_name)
    model = DxfModel()
    model.unit_scale, model.unit_source = resolve_unit_scale(tags, unit)
    scale = model.unit_scale

    entities = _split_entities(tags)

    # --- 面を集める（法線・吸音率はあとでまとめて付ける） ---
    faces = []          # (layer, x1, x2, x3)
    unresolved = set()

    i = 0
    while i < len(entities):
        etype, etags = entities[i]

        if etype == "POLYLINE":
            flag = _int(etags, 70)
            layer = _tag(etags, 8, "0")
            i += 1

            coords, records = [], []
            while i < len(entities) and entities[i][0] == "VERTEX":
                vtags = entities[i][1]
                vflag = _int(vtags, 70)
                if (vflag & VERTEX_FLAG_FACE_RECORD) and not (vflag & VERTEX_FLAG_POLYFACE_MESH):
                    # 面レコード。71〜74 が頂点番号（1 始まり、負値は辺が非表示という意味だけ）
                    records.append([_int(vtags, c) for c in (71, 72, 73, 74)])
                else:
                    coords.append(np.array([_float(vtags, 10), _float(vtags, 20),
                                            _float(vtags, 30)]) * scale)
                i += 1
            if i < len(entities) and entities[i][0] == "SEQEND":
                i += 1

            if not (flag & POLYLINE_FLAG_POLYFACE_MESH):
                model.skipped["非対応エンティティ"] += 1
                continue

            for record in records:
                idx = [abs(v) for v in record if v != 0]
                if len(idx) < 3 or max(idx) > len(coords):
                    model.skipped["面レコード不正"] += 1
                    continue
                pts = [coords[k - 1] for k in idx]
                # 四角形以上は扇状に三角形分割する
                for k in range(1, len(pts) - 1):
                    faces.append((layer, pts[0], pts[k], pts[k + 1]))
            continue

        if etype == "3DFACE":
            layer = _tag(etags, 8, "0")
            pts = [np.array([_float(etags, 10 + k), _float(etags, 20 + k),
                             _float(etags, 30 + k)]) * scale for k in range(4)]
            # 4 点目が 3 点目と同じなら三角形
            if np.allclose(pts[3], pts[2]):
                pts = pts[:3]
            for k in range(1, len(pts) - 1):
                faces.append((layer, pts[0], pts[k], pts[k + 1]))
            i += 1
            continue

        if etype == "POINT":
            layer = (_tag(etags, 8, "0") or "").strip()
            p = np.array([_float(etags, 10), _float(etags, 20), _float(etags, 30)]) * scale
            low = layer.lower()
            if low in [s.lower() for s in source_layers]:
                model.source_points.append(p)
            elif low in [s.lower() for s in receiver_layers]:
                model.receiver_points.append(p)
            else:
                model.skipped["非対応エンティティ"] += 1
            i += 1
            continue

        if etype not in ("VERTEX", "SEQEND"):
            model.skipped["非対応エンティティ"] += 1
        i += 1

    # --- 法線を計算し、縮退面を除く ---
    kept = []
    for layer, x1, x2, x3 in faces:
        n = face_normal(x1, x2, x3)
        if n is None:
            model.skipped["縮退面"] += 1
            continue
        kept.append((layer, x1, x2, x3, n))

    # --- 法線の向きを内向きに揃える ---
    flip = False
    note = ""
    if orient_normals == "auto" and kept:
        v = signed_volume([(f[1], f[2], f[3]) for f in kept])
        if v > 0.0:
            flip = True
            note = f"符号付き体積 {v:+.4f} m³ > 0 → 外向きだったので全反転"
        else:
            note = f"符号付き体積 {v:+.4f} m³ < 0 → すでに内向き"
        if abs(v) < 1.0e-9:
            note += "（体積がほぼ 0。閉じていないか巻き順が不揃いの可能性あり）"
    elif orient_normals in ("inward", "outward"):
        center = np.mean([f[1] for f in kept] + [f[2] for f in kept] + [f[3] for f in kept], axis=0)
        note = f"重心 {np.round(center, 3).tolist()} を基準に面ごとに {orient_normals} へ強制（凸形状前提）"

    for layer, x1, x2, x3, n in kept:
        if orient_normals in ("inward", "outward"):
            to_center = center - x1
            facing_in = np.dot(n, to_center) > 0
            want_in = (orient_normals == "inward")
            if facing_in != want_in:
                n = -n
        elif flip:
            n = -n
        absorption = _resolve_absorption(layer, absorption_table, default_absorption,
                                         band_number, unresolved)
        model.mesh.append(Mesh(x1, x2, x3, n, layer, absorption))
        model.layer_counts[layer] = model.layer_counts.get(layer, 0) + 1

    if model.mesh:
        allpts = np.concatenate([m.vertexes for m in model.mesh])
        model.extents = np.array([allpts.min(axis=0), allpts.max(axis=0)])

    if unresolved:
        used = 0.1 if default_absorption is None else default_absorption
        print(f"[read_dxffile] 警告: 吸音率が未指定のレイヤ {sorted(unresolved)} → {used} を使用")

    if verbose:
        print("[read_dxffile] " + model.summary().replace("\n", "\n[read_dxffile] "))
        if note:
            print(f"[read_dxffile] 法線の向き: {note}")

    return model


def read(file_name, unit=None, absorption_table=None, default_absorption=None,
         orient_normals="auto", band_number=6, verbose=True):
    """メッシュのリストだけを返す簡易版（procedure.py から使う）。

    音源・受音点も欲しい場合は read_model() を使う。
    """
    return read_model(file_name, unit=unit, absorption_table=absorption_table,
                      default_absorption=default_absorption,
                      orient_normals=orient_normals, band_number=band_number,
                      verbose=verbose).mesh


if __name__ == "__main__":
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    model = read_model(os.path.join(here, "test.dxf"))
    print()
    for j, m in enumerate(model.mesh[:4]):
        print(f"mesh[{j}] layer='{m.material}' normal={np.round(m.normal, 3).tolist()}")
        print(f"         vertexes={np.round(m.vertexes, 3).tolist()}")
