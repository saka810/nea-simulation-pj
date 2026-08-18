"""DXF を読んで室形状・吸音率・音源・受音点を取り出す（パイプライン①）。

元コード `backtrace.f90` 132〜283 行。**このパッケージでいちばん大きいモジュール。**
CAD 側の作図ルールは `docs/DXFデータの作り方.md` にまとめてある。

読めるもの:

| DXF のエンティティ | 扱い |
|---|---|
| `POLYLINE`（ポリフェイスメッシュ） | 面としてそのまま読む |
| `3DFACE` | 同上 |
| 閉じた `POLYLINE` / `LWPOLYLINE` | 輪郭を面として読み、三角形に分割する |
| `POINT`（`src` / `rec` レイヤ） | 音源 / 受音点 |

やっていること:

- **単位換算**（ヘッダ `$INSUNITS` を見て m に直す。`unit=` で上書きも可）
- **OCS→WCS 変換**（Arbitrary Axis Algorithm）。
  これが無いと**鉛直な壁が読めない**（押し出し方向が Z でないため。B-18）
- **多角形→三角形分割**。CAD 側で三角形に割る必要はない。
  凹んだ四角形・ねじれた四角形・5 角形以上（耳刈り法）に対応し、面積の照合で破綻を検出する
- **レイヤ名 → 吸音材**の突き合わせ。レイヤ名の先頭の数字（`01__研修室_床` → `1`）でも引ける
- **法線の向きの自動判定**（`orient_normals`）。詳細は `read_model()` の説明を参照
- **シェル診断**。辺の共有で連結成分に分け、閉じているか・体積・法線の向きを報告する
- **作図ミスの自動チェック**（`check_model()`）。計算に入る前に指摘する

法線は「音が通る空気側」を向くのが大原則。その正解を持っているのは CAD モデル自身なので、
既定は `'cad'`（モデルを信じる）。閉じた室なら `'auto'` で内向きへ揃えられる。
"""


import collections
import csv
import os
import re

import numpy as np

from mesh import Mesh

# 周波数バンド数の既定。63〜8k Hz の 8 バンド。
# 63 Hz と 8 kHz を対象外にする運用では 6（125〜4k）を渡す。
# 中心周波数そのものは absorption.octave_bands() が持つ。
DEFAULT_BAND_NUMBER = 8


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
POLYLINE_FLAG_CLOSED = 1            # 閉じたポリライン
POLYLINE_FLAG_3D = 8                # 3D ポリライン
POLYLINE_FLAG_POLYFACE_MESH = 64    # ポリフェイスメッシュ

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
        is_closed       bool        全体が閉じているか（開いた辺が 0 本か）
        open_edges      int         開いた辺（三角形1枚にしか属さない辺）の本数
        volume          float|None  全体が閉じている場合の体積 [m³]
        shells          list[dict]  シェル（連結成分）ごとの診断。analyse_shells() を参照
        winding_consistent bool     頂点の巻き順が一貫しているか
    """

    def __init__(self):
        self.mesh = []
        self.source_points = []
        self.receiver_points = []
        self.unit_scale = 1.0
        self.unit_source = ""
        self.layer_counts = {}
        self.layer_materials = {}   # レイヤ名 → 実際に引けた吸音率表のキー
        self.orient_mode = "cad"    # 実際に適用された法線モード（'auto' の結果もここに入る）
        self.flipped_faces = set()  # CAD の巻き順から反転した面のインデックス
        self.enclosure = None       # 囲まれ具合 0〜1（'auto' で判定したときだけ）
        self.skipped = {"縮退面": 0, "非対応エンティティ": collections.Counter(),
                        "面レコード不正": 0}
        self.face_sources = collections.Counter()
        self.extents = None
        self.is_closed = False
        self.open_edges = 0
        self.volume = None
        self.shells = []
        self.winding_consistent = True
        self.polygon_notes = {"ねじれた四角形": 0, "最大ねじれ量": 0.0,
                             "最大ねじれ実寸": 0.0,
                             "凹み対応で対角線を変更": 0,
                             "耳刈り法で分割": 0, "分割に失敗": 0,
                             "最大面積誤差": 0.0}

    def summary(self):
        lines = [
            f"三角形メッシュ: {len(self.mesh)} 枚",
            f"単位: 1 CAD単位 = {self.unit_scale} m（{self.unit_source}）",
            f"レイヤ別の枚数: {self.layer_counts}",
        ]
        if self.layer_materials:
            lines.append("レイヤ→材料: " + ", ".join(
                f"{layer}→{self.layer_materials[layer]}"
                for layer in sorted(self.layer_materials)))
        lines += [
            f"音源: {[np.round(p, 4).tolist() for p in self.source_points]}",
            f"受音点: {[np.round(p, 4).tolist() for p in self.receiver_points]}",
        ]
        if self.extents is not None:
            size = self.extents[1] - self.extents[0]
            lines.append(f"寸法: {np.round(size, 4).tolist()} m "
                         f"(min={np.round(self.extents[0], 4).tolist()})")

        # シェルごとの診断。法線の向きの正しさはここを見て判断する
        lines.append(f"シェル（連結した面のかたまり）: {len(self.shells)} 個"
                     + ("" if self.winding_consistent
                        else "  ★巻き順が一貫していません（隣り合う面で法線が反対を向いている）"))
        for k, s in enumerate(self.shells):
            tag = "外殻" if s["is_outer"] else "内側"
            size = np.round(s["bbox"][1] - s["bbox"][0], 3).tolist()
            if s["closed"]:
                want = "inward" if s["is_outer"] else "outward"
                mark = "OK" if s["normals"] == want else f"★要確認（{want} が空気側）"
                lines.append(f"  シェル{k} [{tag}] 面{len(s['faces'])}枚 閉 "
                             f"体積{s['volume']:.4f}m³ 法線={s['normals']} {mark} 寸法{size}")
            else:
                lines.append(f"  シェル{k} [{tag}] 面{len(s['faces'])}枚 開"
                             f"（開いた辺{s['open_edges']}本）法線=CAD のまま 寸法{size}")
        if not self.is_closed:
            lines.append("  ※開いた形状（一面反射など）も計算できます")

        if self.face_sources:
            lines.append(f"面の元になったエンティティ: {dict(self.face_sources)}")
        if self.polygon_notes["耳刈り法で分割"]:
            lines.append(f"三角形分割: 5角形以上を {self.polygon_notes['耳刈り法で分割']} 枚"
                         f"耳刈り法で分割（面積の相対誤差 最大 "
                         f"{self.polygon_notes['最大面積誤差']:.2e}）")
        if self.polygon_notes["分割に失敗"]:
            lines.append(f"★三角形分割に失敗した多角形が "
                         f"{self.polygon_notes['分割に失敗']} 枚あります"
                         f"（自己交差しているか、同一平面上にない可能性）")
        if self.polygon_notes["凹み対応で対角線を変更"]:
            lines.append(f"三角形分割: 凹んだ四角形 "
                         f"{self.polygon_notes['凹み対応で対角線を変更']} 枚を"
                         f"内側を通る対角線で分割しました")
        if self.polygon_notes["ねじれた四角形"]:
            mm = self.polygon_notes["最大ねじれ実寸"] * 1000.0
            impact = ("音響的には無視できる大きさです" if mm < 1.0
                      else "音響計算に影響しうる大きさです" if mm < 20.0
                      else "★形状が意図とずれている可能性があります")
            lines.append(
                f"三角形分割: ねじれた四角形 {self.polygon_notes['ねじれた四角形']} 枚"
                f"（4 点が同一平面上にない）"
                f"\n  最大ねじれ {mm:.3f} mm（相対 "
                f"{self.polygon_notes['最大ねじれ量']:.2e}）… {impact}"
                f"\n  切る対角線によって形が変わります。気になる場合は CAD 側で"
                f"平面に直すか三角形で作り直してください"
                f"（作図方法は docs/DXFデータの作り方.md 1節）")
        skipped = {k: (dict(v) if isinstance(v, collections.Counter) else v)
                   for k, v in self.skipped.items()
                   if (len(v) if isinstance(v, collections.Counter) else v)}
        if skipped:
            lines.append(f"読み飛ばし: {skipped}")
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


def _tags_all(entity_tags, code):
    """同じグループコードが複数回現れるもの（LWPOLYLINE の頂点など）を全部返す。"""
    return [v for c, v in entity_tags if c == code]


# ------------------------------------------------------------------------------
# OCS（オブジェクト座標系）→ WCS（ワールド座標系）
# ------------------------------------------------------------------------------

# 任意軸アルゴリズムの分岐しきい値（DXF 仕様で 1/64 と決まっている）
ARBITRARY_AXIS_EPS = 1.0 / 64.0


def ocs_axes(extrusion):
    """押し出し方向から OCS の 3 軸を作る（DXF の Arbitrary Axis Algorithm）。

    LWPOLYLINE や CIRCLE のような「平面図形」は、頂点を **その図形が乗る平面上の
    2 次元座標（OCS）** で持っている。押し出し方向 N（グループコード 210/220/230）が
    Z 軸なら OCS = WCS なのでそのまま使えるが、**壁のような鉛直な面は N が水平を向く**
    ため、この変換を通さないと座標がまるで違う場所になる。

    軸の取り方は DXF 仕様で決められていて、任意性はない。
    N が Z 軸に近いときだけ基準を Y 軸に切り替えるのは、
    Z 軸との外積が退化して軸が決まらなくなるのを避けるため。
    """
    normal = np.asarray(extrusion, dtype=float)
    length = np.linalg.norm(normal)
    if length == 0.0:
        return np.eye(3)
    normal = normal / length

    if abs(normal[0]) < ARBITRARY_AXIS_EPS and abs(normal[1]) < ARBITRARY_AXIS_EPS:
        reference = np.array([0.0, 1.0, 0.0])   # N が Z 軸に近い
    else:
        reference = np.array([0.0, 0.0, 1.0])
    axis_x = np.cross(reference, normal)
    axis_x /= np.linalg.norm(axis_x)
    axis_y = np.cross(normal, axis_x)
    return np.array([axis_x, axis_y, normal])


def ocs_to_wcs(point, extrusion):
    """OCS の点（x, y, 高度）をワールド座標に直す。"""
    axes = ocs_axes(extrusion)
    p = np.asarray(point, dtype=float)
    return p[0] * axes[0] + p[1] * axes[1] + p[2] * axes[2]


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

class AbsorptionTable(dict):
    """{キー: 吸音率 ndarray} に、キー → 材料名の対応を添えたもの。

    ID でも材料名でも引けるようにキーを 2 通り登録するので、
    「ID で引いたときに何の材料だったか」を後から言えるようにしておく。
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = {}


def read_absorption_csv(file_name, band_number=DEFAULT_BAND_NUMBER):
    """吸音率 CSV を読んで {キー: ndarray(band_number,)} を返す。

    2 つの形式を自動判別する。

    (A) 元コードの absorption.csv 形式（ID + 材料名 + 吸音率）:
            1,Concrete wall,0.01,0.02,0.02,0.02,0.03,0.04
        → **ID と材料名の両方をキーに登録する**ので、DXF のレイヤ名がどちらでも引ける
    (B) 材料名 + 吸音率:
            コンクリート,0.02,0.02,0.03,0.03,0.04,0.05

    1 行目がヘッダでもよい（数値に変換できない行は読み飛ばす）。
    行頭 # の行と空行も読み飛ばす。
    """
    # 元コード付属の absorption.csv は CP932（Shift-JIS）なので、順に試してデコードする
    last_error = None
    text = None
    for encoding in ("utf-8-sig", "cp932", "latin-1"):
        try:
            with open(file_name, encoding=encoding, newline="") as f:
                text = f.read()
            break
        except UnicodeDecodeError as e:
            last_error = e
    if text is None:
        raise last_error

    table = AbsorptionTable()
    padded = []          # 列数が足りずに補った材料（まとめて 1 回だけ知らせる）
    for row in csv.reader(text.splitlines()):
        if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
            continue
        cells = [c.strip() for c in row]

        # 2 列目が数値なら形式 (B)、数値でなければ形式 (A)（ID + 材料名）
        second_is_number = False
        if len(cells) >= 2:
            try:
                float(cells[1])
                second_is_number = True
            except ValueError:
                second_is_number = False

        if second_is_number:
            keys, values_raw = [cells[0]], cells[1:1 + band_number]
        elif len(cells) >= 3:
            keys, values_raw = [cells[0], cells[1]], cells[2:2 + band_number]
        else:
            continue

        try:
            values = [float(x) for x in values_raw]
        except ValueError:
            continue  # ヘッダ行

        if not values:
            continue
        if len(values) < band_number:
            padded.append((keys[-1], len(values)))
            values = values + [values[-1]] * (band_number - len(values))

        arr = np.array(values[:band_number], dtype=float)
        for key in keys:
            if key:
                table[key] = arr
                table.names[key] = keys[-1]   # ID で引いても材料名が分かるように

    if padded:
        columns = sorted({n for _, n in padded})
        print(f"[read_dxffile] 注意: 吸音率の列が {columns} 個しかない材料が "
              f"{len(padded)} 種あります（{band_number} バンドで計算）。"
              f"最後の値で埋めました。"
              f"例: {', '.join(name for name, _ in padded[:3])} …")
    return table


def _add_triangles(model, faces, layer, points):
    """多角形を三角形に分割して faces に追加し、ねじれ等を model 側に記録する。"""
    tris, info = triangulate_polygon(points)
    if not tris:
        model.skipped["縮退面"] += 1
        return
    if info["warp"] > WARP_TOLERANCE:
        model.polygon_notes["ねじれた四角形"] += 1
        model.polygon_notes["最大ねじれ量"] = max(model.polygon_notes["最大ねじれ量"],
                                                 info["warp"])
        model.polygon_notes["最大ねじれ実寸"] = max(model.polygon_notes["最大ねじれ実寸"],
                                                   info["warp_distance"])
    if info["diagonal_changed"]:
        model.polygon_notes["凹み対応で対角線を変更"] += 1
    if info["ear_clipped"]:
        model.polygon_notes["耳刈り法で分割"] += 1
        model.polygon_notes["最大面積誤差"] = max(model.polygon_notes["最大面積誤差"],
                                               info["area_error"])
    if info["failed"]:
        model.polygon_notes["分割に失敗"] += 1
    for t in tris:
        faces.append((layer, t[0], t[1], t[2]))


def layer_number(layer):
    """レイヤ名の先頭の数字を材料 ID として取り出す。無ければ None。

    CAD 側の運用で `01__研修室_床` のように **先頭 2 桁を吸音率表の ID** にして
    レイヤ名を付けることがある。名前そのものは材料名と一致しないので、
    この番号で引けるようにしておく。先頭のゼロは落とす（`01` → `1`）。
    """
    m = re.match(r"^\s*(\d+)", layer or "")
    if not m:
        return None
    return str(int(m.group(1)))


def _resolve_absorption(layer, absorption_table, default_absorption, band_number,
                        unresolved, resolved=None):
    names = getattr(absorption_table, "names", {})
    if absorption_table:
        if layer in absorption_table:
            if resolved is not None:
                resolved[layer] = names.get(layer, layer)
            return absorption_table[layer]
        # レイヤ名が材料名と一致しないとき、先頭の番号を材料 ID とみなして引く
        number = layer_number(layer)
        if number is not None and number in absorption_table:
            if resolved is not None:
                resolved[layer] = f"{number}:{names.get(number, number)}"
            return absorption_table[number]
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

    ★開いた形状（閉じていない形状）では意味を持たない。必ず open_edge_count() で
      閉じているかを確認してから使うこと。
    """
    total = 0.0
    for x1, x2, x3 in triangles:
        total += float(np.dot(x1, np.cross(x2, x3)))
    return total / 6.0


# ------------------------------------------------------------------------------
# 多角形 → 三角形への分割
# ------------------------------------------------------------------------------
#
# CAD 側で三角形を手作りする必要はない。四角形で描いてよく、ここで 2 枚に分割する。
#
# そもそも DXF から入ってくる面は最大 4 頂点しかない：
#   ・POLYLINE ポリフェイスメッシュの面レコードは頂点番号がグループコード 71〜74 の 4 つ
#   ・3DFACE も 4 隅まで
# したがって 5 角形以上は考えなくてよい（下のコードは保険として一般の n 角形も扱う）。
#
# 注意すべきは 2 つだけ：
#   ・ねじれた四角形（4 点が同一平面上にない）… どの対角線で切るかで形が変わる。警告を出す
#   ・凹んだ四角形                            … 単純な扇状分割だと多角形の外に三角形ができる。
#                                              内側を通る対角線を選んで対処する

WARP_TOLERANCE = 1.0e-6      # 平面からのずれ / 代表辺長 がこれを超えたら「ねじれ」とみなす


def _dedupe_consecutive(points, tol=1.0e-12):
    """連続する重複頂点を落とす（3DFACE が三角形を表すときに 4 点目を 3 点目と同じにする等）。"""
    out = []
    for p in points:
        if not out or not np.allclose(out[-1], p, atol=tol):
            out.append(p)
    if len(out) > 1 and np.allclose(out[0], out[-1], atol=tol):
        out.pop()
    return out


def quad_warp_distance(points):
    """四角形のねじれ量を**実寸**で返す [m]。

    最初の 3 点が作る平面から 4 点目までの距離。0 なら完全に平面。
    相対値より実寸のほうが「CAD で直すべきか」の判断がしやすい。
    """
    p0, p1, p2, p3 = [np.asarray(p, dtype=float) for p in points[:4]]
    n = np.cross(p1 - p0, p2 - p0)
    length = np.linalg.norm(n)
    if length < DEGENERATE_EPS:
        return 0.0
    return abs(float(np.dot(n / length, p3 - p0)))


def quad_warp(points):
    """四角形のねじれ具合（平面からのずれ ÷ 代表辺長）を返す。0 なら完全に平面。"""
    p0, p1, p2, p3 = [np.asarray(p, dtype=float) for p in points[:4]]
    scale = max(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p0),
                np.linalg.norm(p3 - p0), 1.0e-30)
    return quad_warp_distance(points) / scale


def _diagonal_is_inside(points, i, j, normal):
    """対角線 points[i]-points[j] が多角形の内側を通るか。

    他の 2 頂点が対角線の両側に分かれていれば内側を通る。
    どちらかが対角線上に乗っている（符号 0）場合も内側扱いにする
    ― その側の三角形が面積ゼロになるだけで、残る三角形が多角形を正しく覆う
    （面積ゼロの三角形は後段の face_normal() で縮退面として捨てられる）。
    """
    a, b = np.asarray(points[i], dtype=float), np.asarray(points[j], dtype=float)
    others = [k for k in range(len(points)) if k not in (i, j)]
    signs = []
    for k in others:
        c = np.asarray(points[k], dtype=float)
        signs.append(float(np.dot(np.cross(b - a, c - a), normal)))
    if all(abs(s) < DEGENERATE_EPS for s in signs):
        return False        # 全頂点が一直線＝退化した多角形
    return signs[0] * signs[1] <= 0.0


def triangulate_polygon(points):
    """多角形（同一平面・単純多角形を想定）を三角形のリストに分割する。

    戻り値: (三角形のリスト, 情報dict)
      情報dict のキー: 'warp'（ねじれ量。四角形のみ）, 'diagonal_changed'（対角線を変えたか）
    """
    pts = _dedupe_consecutive(list(points))
    info = {"warp": 0.0, "warp_distance": 0.0, "diagonal_changed": False,
            "ear_clipped": False, "area_error": 0.0, "failed": False}

    if len(pts) < 3:
        return [], info
    if len(pts) == 3:
        return [(pts[0], pts[1], pts[2])], info

    # 多角形の代表法線（面積重み付き。凹んでいても向きが安定する）
    normal = np.zeros(3)
    for k in range(len(pts)):
        normal += np.cross(np.asarray(pts[k], dtype=float),
                           np.asarray(pts[(k + 1) % len(pts)], dtype=float))
    if np.linalg.norm(normal) < DEGENERATE_EPS:
        return [], info
    normal = normal / np.linalg.norm(normal)

    if len(pts) == 4:
        info["warp"] = quad_warp(pts)
        info["warp_distance"] = quad_warp_distance(pts)
        # 内側を通る対角線を選ぶ。0-2 が使えなければ 1-3 を使う（凹んだ四角形への対処）
        if _diagonal_is_inside(pts, 0, 2, normal):
            tris = [(pts[0], pts[1], pts[2]), (pts[0], pts[2], pts[3])]
        else:
            info["diagonal_changed"] = True
            tris = [(pts[1], pts[2], pts[3]), (pts[1], pts[3], pts[0])]
        return tris, info

    # 5 角形以上 … 耳刈り法（ear clipping）
    # 「耳」＝ 凸な頂点で、その両隣を結んだ三角形の中に他の頂点が入らないもの。
    # 耳を 1 つずつ切り落としていけば、凹んだ多角形でも正しく分割できる。
    tris, ok = _ear_clip(pts, normal)
    info["ear_clipped"] = True
    info["area_error"] = _area_mismatch(pts, tris, normal)
    if not ok or info["area_error"] > 1.0e-6:
        info["failed"] = True
    return tris, info


def polygon_area(points, normal):
    """多角形の面積（法線方向に射影した符号なし面積）。"""
    total = np.zeros(3)
    n = len(points)
    for k in range(n):
        total = total + np.cross(np.asarray(points[k], dtype=float),
                                 np.asarray(points[(k + 1) % n], dtype=float))
    return abs(float(np.dot(total, normal))) / 2.0


def _area_mismatch(points, triangles, normal):
    """三角形分割の面積が元の多角形と合っているかの相対誤差。

    分割が正しければ 0 になる。自己交差した多角形などで破綻すると大きくなるので、
    分割の妥当性チェックとしてそのまま使える。
    """
    truth = polygon_area(points, normal)
    if truth < DEGENERATE_EPS:
        return 0.0
    got = sum(0.5 * float(np.linalg.norm(np.cross(np.asarray(t[1]) - np.asarray(t[0]),
                                                  np.asarray(t[2]) - np.asarray(t[0]))))
              for t in triangles)
    return abs(got - truth) / truth


def _is_convex_corner(prev_p, cur_p, next_p, normal):
    return float(np.dot(np.cross(np.asarray(cur_p) - np.asarray(prev_p),
                                 np.asarray(next_p) - np.asarray(cur_p)), normal)) > 0.0


def _point_in_triangle(p, a, b, c, normal, eps=1.0e-12):
    """点 p が三角形 abc の内部（辺上を含む）にあるか。すべて同一平面上を仮定。"""
    p, a, b, c = [np.asarray(x, dtype=float) for x in (p, a, b, c)]
    d1 = float(np.dot(np.cross(b - a, p - a), normal))
    d2 = float(np.dot(np.cross(c - b, p - b), normal))
    d3 = float(np.dot(np.cross(a - c, p - c), normal))
    return (d1 >= -eps and d2 >= -eps and d3 >= -eps) or \
           (d1 <= eps and d2 <= eps and d3 <= eps)


def _ear_clip(points, normal):
    """耳刈り法で多角形を三角形に分割する。戻り値: (三角形のリスト, 成功したか)。"""
    ring = list(range(len(points)))
    tris = []
    guard = 0
    while len(ring) > 3 and guard < len(points) * len(points) + 10:
        guard += 1
        for m in range(len(ring)):
            i = ring[m - 1]
            j = ring[m]
            k = ring[(m + 1) % len(ring)]
            if not _is_convex_corner(points[i], points[j], points[k], normal):
                continue
            # 他の頂点が三角形 ijk の中に入っていないこと
            if any(_point_in_triangle(points[q], points[i], points[j], points[k], normal)
                   for q in ring if q not in (i, j, k)):
                continue
            tris.append((points[i], points[j], points[k]))
            ring.pop(m)
            break
        else:
            return tris, False      # 耳が見つからない＝自己交差などで分割できない
    if len(ring) == 3:
        tris.append((points[ring[0]], points[ring[1]], points[ring[2]]))
        return tris, True
    return tris, False


def _point_key(point, digits):
    return tuple(np.round(np.asarray(point, dtype=float), digits))


def _edge_map(triangles, tol=1.0e-9):
    """無向辺 → その辺を持つ面インデックスのリスト、および有向辺の出現回数を返す。"""
    digits = int(round(-np.log10(tol)))
    undirected = {}
    directed = {}
    for j, (x1, x2, x3) in enumerate(triangles):
        for a, b in ((x1, x2), (x2, x3), (x3, x1)):
            ka, kb = _point_key(a, digits), _point_key(b, digits)
            directed[(ka, kb)] = directed.get((ka, kb), 0) + 1
            edge = (ka, kb) if ka <= kb else (kb, ka)
            undirected.setdefault(edge, []).append(j)
    return undirected, directed


def open_edge_count(triangles, tol=1.0e-9):
    """開いた辺（三角形 1 枚にしか属さない辺）の本数を返す。0 なら閉じている。

    閉じた多面体では、すべての辺がちょうど 2 枚の三角形に共有される。
    この本数を数えるだけで「閉じているか」が判定できる。
    一面だけの壁のような開いた形状では、外周の辺が 1 回しか現れないので 0 にならない。
    """
    undirected, _ = _edge_map(triangles, tol)
    return sum(1 for faces in undirected.values() if len(faces) != 2)


def winding_is_consistent(triangles, tol=1.0e-9):
    """頂点の巻き順が一貫しているかを返す。

    巻き順が一貫した面同士が辺を共有するとき、その辺は互いに逆向きに現れる
    （面 A が a→b なら面 B は b→a）。同じ向きで 2 回現れる辺があれば、
    その 2 面の巻き順は裏返っている＝法線が反対を向いている。
    """
    _, directed = _edge_map(triangles, tol)
    return not any(count > 1 for count in directed.values())


def mesh_shells(triangles, tol=1.0e-9):
    """辺の共有で連結成分（シェル）に分割し、面インデックスのリストを返す。

    室の外殻と、室内に置いた厚みのある家具などは別々のシェルになる。
    法線の望ましい向きがシェルごとに違う（外殻は内向き、家具は外向き＝空気側）ので、
    診断や自動補正はシェル単位で行う必要がある。
    """
    undirected, _ = _edge_map(triangles, tol)
    neighbours = {j: set() for j in range(len(triangles))}
    for faces in undirected.values():
        for a in faces:
            for b in faces:
                if a != b:
                    neighbours[a].add(b)

    seen = set()
    shells = []
    for start in range(len(triangles)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        group = []
        while stack:
            j = stack.pop()
            group.append(j)
            for k in neighbours[j]:
                if k not in seen:
                    seen.add(k)
                    stack.append(k)
        shells.append(sorted(group))
    return shells


def _count_crossings(origins, directions, triangles, eps=1.0e-9):
    """各レイが三角形群と何回交わるかを数える（両面。裏からの交差も数える）。

    origins (K,3) / directions (K,3) / triangles list[(v1,v2,v3)] → (K,) の交差回数。
    """
    origins = np.asarray(origins, dtype=float)
    directions = np.asarray(directions, dtype=float)
    v0 = np.array([t[0] for t in triangles], dtype=float)
    v1 = np.array([t[1] for t in triangles], dtype=float)
    v2 = np.array([t[2] for t in triangles], dtype=float)

    # Möller–Trumbore。平面の法線を使わないので、CAD の巻き順に左右されない
    edge1 = v1 - v0
    edge2 = v2 - v0
    pvec = np.cross(directions[:, None, :], edge2[None, :, :])      # (K,M,3)
    det = np.einsum("kmi,mi->km", pvec, edge1)
    parallel = np.abs(det) < eps
    inv_det = np.where(parallel, 1.0, det)

    tvec = origins[:, None, :] - v0[None, :, :]
    u = np.einsum("kmi,kmi->km", tvec, pvec) / inv_det
    qvec = np.cross(tvec, edge1[None, :, :])
    v = np.einsum("kmi,ki->km", qvec, directions) / inv_det
    t = np.einsum("kmi,mi->km", qvec, edge2) / inv_det

    hit = (~parallel) & (u >= 0.0) & (v >= 0.0) & (u + v <= 1.0) & (t > eps)
    return hit.sum(axis=1)


# 内向き判定に使う試行方向の本数（多数決を取るので奇数にする）
INWARD_PROBE_DIRECTIONS = 9


def orient_inward(triangles, normals, probes=INWARD_PROBE_DIRECTIONS, seed=0):
    """面ごとに「法線が室内側を向いているか」を判定し、反転すべき面の集合を返す。

    **CAD で法線を意識せずに描いたモデル**（床・壁・天井を 1 枚ずつ描いただけのもの）は
    巻き順と押し出し方向で法線が決まるため、床は上・天井も上、壁はバラバラ、
    といったことが普通に起きる。それだと衝突判定で壁がすり抜ける。

    判定は**レイの偶奇**による。閉じた面に囲まれた領域の内側から外へレイを飛ばすと、
    境界を必ず**奇数回**横切る（外側からなら偶数回）。そこで面の重心を法線側へわずかに
    浮かせた点から外へレイを飛ばし、交差回数が奇数なら「法線側が室内」＝そのまま、
    偶数なら反転する。

    1 方向だけだと、辺や頂点をかすめたときに数え間違える。ここでは方向を
    `probes` 本ばらまいて**多数決**を取る。方向は固定シードで生成するので、
    同じモデルなら毎回同じ結果になる（再現性のため）。

    戻り値: (反転する面インデックスの集合, 判定に迷った面の数)
    """
    triangles = list(triangles)
    normals = np.asarray(normals, dtype=float)
    centroids = np.array([(t[0] + t[1] + t[2]) / 3.0 for t in triangles])

    # 面の大きさに合わせて浮かせる量を決める（大きいモデルでも小さいモデルでも効くように）
    span = float(np.max(np.ptp(np.concatenate([np.asarray(t) for t in triangles]), axis=0)))
    offset = max(span * 1.0e-7, 1.0e-9)

    rng = np.random.default_rng(seed)
    flip = set()
    ambiguous = 0
    for j in range(len(triangles)):
        origin = centroids[j] + offset * normals[j]
        # 法線側の半球へ飛ばす。法線そのものと、それを軸にばらけさせた方向
        directions = rng.normal(size=(probes, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        # 法線と逆向きのものは反転して、必ず法線側の半球にする
        sign = np.sign(directions @ normals[j])
        sign[sign == 0.0] = 1.0
        directions *= sign[:, None]
        directions[0] = normals[j]

        counts = _count_crossings(np.repeat(origin[None, :], probes, axis=0),
                                  directions, triangles)
        odd = int(np.count_nonzero(counts % 2 == 1))
        if odd * 2 == probes:
            ambiguous += 1
        if odd * 2 < probes:      # 偶数が多数 ＝ 法線側は室外 → 反転
            flip.add(j)
    return flip, ambiguous


def encloses_point(triangles, point, samples=256, seed=0):
    """`point` がこの面群に**囲まれているか**の度合い（0〜1）を返す。

    点から全方向へレイを飛ばし、**面に当たった割合**を返す。
    完全に閉じた室の内側なら 1.0、一面だけの反射板なら 0.5 前後、外側なら 0 に近い。

    法線の自動補正（`orient_normals='auto'`）で「室として扱ってよいか」を判断するために使う。
    レイの偶奇による内向き判定は**囲まれた形状でしか意味がない**ので、
    一面反射板のような開いた形状に適用すると、反射させたい側を逆に向けてしまう。

    方向は固定シードで生成するので、同じモデルなら毎回同じ値になる。
    """
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(int(samples), 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    origins = np.repeat(np.asarray(point, dtype=float)[None, :], len(directions), axis=0)
    counts = _count_crossings(origins, directions, triangles)
    return float(np.count_nonzero(counts > 0)) / len(directions)


# これ以上なら「室として囲まれている」とみなす閾値。
# 完全な室なら 1.0 になるので、隙間や T 字接合ぶんの余裕を見て 0.95 にしてある。
ENCLOSURE_THRESHOLD = 0.95


def _bbox(points):
    p = np.asarray(points, dtype=float)
    return p.min(axis=0), p.max(axis=0)


def _bbox_contains(outer, inner, eps=1.0e-9):
    return bool(np.all(outer[0] - eps <= inner[0]) and np.all(inner[1] <= outer[1] + eps))


def analyse_shells(triangles, tol=1.0e-9):
    """シェルごとに閉/開・体積・法線の向き・外殻かどうかを調べる。

    戻り値: list[dict]  キーは
        faces      … 面インデックスのリスト
        closed     … そのシェルが閉じているか
        open_edges … 開いた辺の本数
        volume     … 閉じている場合の体積 [m³]（開いていれば None）
        normals    … 'inward' / 'outward' / 'unknown'（そのシェル自身から見た向き）
        is_outer   … 他のすべてのシェルを内包する外殻か
    """
    shells = mesh_shells(triangles, tol)
    result = []
    for faces in shells:
        tris = [triangles[j] for j in faces]
        oe = open_edge_count(tris, tol)
        closed = (oe == 0)
        vol = signed_volume(tris) if closed else None
        if vol is None:
            normals = "unknown"
        elif vol > 0.0:
            normals = "outward"
        else:
            normals = "inward"
        pts = np.concatenate([np.array(t) for t in tris])
        result.append({"faces": faces, "closed": closed, "open_edges": oe,
                       "volume": None if vol is None else abs(vol),
                       "normals": normals, "is_outer": False, "bbox": _bbox(pts)})

    # 他のすべてを内包するシェルを外殻とする
    for i, s in enumerate(result):
        if all(_bbox_contains(s["bbox"], o["bbox"]) for k, o in enumerate(result) if k != i):
            s["is_outer"] = True
            break
    else:
        if result:
            # 内包関係がはっきりしなければ、いちばん大きいものを外殻とみなす
            biggest = max(range(len(result)),
                          key=lambda k: float(np.prod(result[k]["bbox"][1] - result[k]["bbox"][0])))
            result[biggest]["is_outer"] = True
    return result


# ------------------------------------------------------------------------------
# 本体
# ------------------------------------------------------------------------------

def read_model(file_name, unit=None, absorption_table=None, default_absorption=None,
               orient_normals="cad", reference_point=None, band_number=DEFAULT_BAND_NUMBER,
               source_layers=DEFAULT_SOURCE_LAYERS,
               receiver_layers=DEFAULT_RECEIVER_LAYERS,
               flip_faces=None, verbose=True):
    """DXF を読んで DxfModel を返す。

    **閉じた室でも、一面だけの壁のような開いた形状でも読める**（音線追跡側も
    開いた形状に対応している。当たる壁がなくなった音線はそこで打ち切られる）。

    【法線の向きの大原則】法線は「音が通る空気側」を向く。この正解を持っているのは
    CAD モデル自身なので、既定は 'cad'（モデルを信じる）。
      ・両面で反射させたい物体（机など）は CAD 側で厚みを持たせる
        → 上面は上向き法線、下面は下向き法線となり、どちらも空気側を向く
      ・片面だけで良ければ、反射させたい側に法線を向けてモデルを作る
    法線を「音源方向へ向ける」ような補正はしない（凸凹の壁や宙に浮いた家具で破綻する）。

    引数:
        unit              None なら $INSUNITS から自動判定。'mm' / 'm' などで明示指定も可
        absorption_table  {レイヤ名: ndarray} または吸音率 CSV のパス
        default_absorption テーブルに無いレイヤに使う値（None なら 0.1 で警告）
        orient_normals    法線の向きの扱い
          'cad'（既定）… CAD の巻き順をそのまま使う
          'auto'       … **閉じている／室として囲まれていれば内向きに揃え、
                         開いた形状なら CAD のまま**。開いた形状に内向き補正をかけると
                         一面反射板の反射させたい側が逆になるので、
                         `encloses_point()` で実際に囲まれているか確かめてから決める
          'flip'       … 全反転（元コード 276行 ynnmrev='y' に相当）
          'shells'     … シェル（連結成分）単位で空気側へ揃える。
                         外殻は内向き、内側の物体は外向き。
                         巻き順が一貫していない場合は補正を中止して 'cad' と同じ挙動になる
          'inward'     … 面ごとにレイの偶奇で室内側へ揃える。面のつながりを要求しない
        flip_faces        **CAD の巻き順から反転する面インデックスの絶対集合**。
                          `normal_editor.py` が作り、`project.py` が normals.json に保存する。
                          渡すと `orient_normals` の判定を**丸ごと置き換える**（差分ではない）。
                          人が確認し終えた最終状態そのものなので、
                          保存したものを読めば必ず同じ法線になる
        reference_point   使わない（旧 'toward' 用。渡すと警告する）
        source_layers     音源として扱う POINT のレイヤ名
        receiver_layers   受音点として扱う POINT のレイヤ名

    どの向きになっているかは戻り値の `summary()` のシェル診断で確認できる。
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
                # ポリフェイスメッシュではない POLYLINE。
                # 「閉じたポリライン」は面の輪郭を描いたものとみなして三角形に分割する
                # （CAD で輪郭だけ描くのはよくある作り方。test2.dxf がこの形式）。
                if flag & POLYLINE_FLAG_CLOSED and len(coords) >= 3:
                    _add_triangles(model, faces, layer, coords)
                    model.face_sources["閉じたポリライン"] += 1
                else:
                    model.skipped["非対応エンティティ"][
                        "POLYLINE(開いた線)"] += 1
                continue

            model.face_sources["ポリフェイスメッシュ"] += 1
            for record in records:
                idx = [abs(v) for v in record if v != 0]
                if len(idx) < 3 or max(idx) > len(coords):
                    model.skipped["面レコード不正"] += 1
                    continue
                pts = [coords[k - 1] for k in idx]
                _add_triangles(model, faces, layer, pts)
            continue

        if etype == "3DFACE":
            layer = _tag(etags, 8, "0")
            pts = [np.array([_float(etags, 10 + k), _float(etags, 20 + k),
                             _float(etags, 30 + k)]) * scale for k in range(4)]
            # 三角形を表すときは 4 点目が 3 点目と同じになる（_add_triangles 側で重複を落とす）
            _add_triangles(model, faces, layer, pts)
            model.face_sources["3DFACE"] += 1
            i += 1
            continue

        if etype == "LWPOLYLINE":
            # 軽量ポリライン（AutoCAD の既定の 2D ポリライン）。閉じていれば面として読む。
            layer = _tag(etags, 8, "0")
            flag = _int(etags, 70)
            xs = _tags_all(etags, 10)
            ys = _tags_all(etags, 20)
            elevation = _float(etags, 38, 0.0)
            extrusion = np.array([_float(etags, 210, 0.0), _float(etags, 220, 0.0),
                                  _float(etags, 230, 1.0)])
            if not (flag & 1) or min(len(xs), len(ys)) < 3:
                model.skipped["非対応エンティティ"]["LWPOLYLINE(開いた線)"] += 1
            elif np.allclose(extrusion, [0.0, 0.0, 1.0]):
                pts = [np.array([float(x), float(y), elevation]) * scale
                       for x, y in zip(xs, ys)]
                _add_triangles(model, faces, layer, pts)
                model.face_sources["閉じたLWPOLYLINE"] += 1
            else:
                # 押し出し方向が Z 以外 ＝ 頂点が OCS で書かれている（鉛直な壁など）。
                # 変換しないと座標がまるで違う場所になるので、必ず通す
                pts = [ocs_to_wcs([float(x), float(y), elevation], extrusion) * scale
                       for x, y in zip(xs, ys)]
                _add_triangles(model, faces, layer, pts)
                model.face_sources["閉じたLWPOLYLINE(OCS)"] += 1
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
                model.skipped["非対応エンティティ"][f"POINT(レイヤ'{layer}')"] += 1
            i += 1
            continue

        if etype not in ("VERTEX", "SEQEND"):
            model.skipped["非対応エンティティ"][etype] += 1
        i += 1

    # --- 法線を計算し、縮退面を除く ---
    kept = []
    for layer, x1, x2, x3 in faces:
        n = face_normal(x1, x2, x3)
        if n is None:
            model.skipped["縮退面"] += 1
            continue
        kept.append((layer, x1, x2, x3, n))

    # --- 法線の向きを決める ---
    # ★大原則：法線は「音が通る空気側」を向く。この向きの正解を持っているのは CAD モデル自身。
    #   したがって既定は 'cad'（モデルを信じる）。
    #
    #   ・両面で反射させたい物体（机など）は CAD 側で厚みを持たせる。
    #     上面は上向き法線、下面は下向き法線となり、どちらも空気側を向く
    #   ・片面だけで良ければ、反射させたい側に法線を向けてモデルを作る
    #
    #   「面ごとに音源へ向ける」ようなことをしてはいけない。凸凹の壁や宙に浮いた家具で
    #   向きが破綻し、衝突判定がおかしくなる。
    triangles = [(f[1], f[2], f[3]) for f in kept]
    flip_all = False
    per_face_flip = None    # 'shells' モードで反転する面インデックスの集合
    note = ""

    if kept:
        model.open_edges = open_edge_count(triangles)
        model.is_closed = (model.open_edges == 0)
        model.winding_consistent = winding_is_consistent(triangles)
        model.shells = analyse_shells(triangles)
        if model.is_closed:
            model.volume = abs(signed_volume(triangles))

    mode = orient_normals
    auto_note = ""
    if mode == "auto":
        # **閉じている（または室として囲まれている）なら内向きに揃え、
        #   開いた形状なら CAD のまま**という自動判定。
        # 開いた形状に内向き補正をかけると、一面反射板の反射させたい側が逆になるので、
        # 「囲まれているか」を実際にレイを飛ばして確かめてから決める。
        if not kept:
            mode, auto_note = "cad", "自動判定: 面が無い / "
        elif model.is_closed and model.winding_consistent:
            mode = "shells"
            auto_note = "自動判定: 閉じていて巻き順も一貫 → シェル単位で空気側へ / "
        else:
            probe = (model.source_points[0] if model.source_points
                     else np.mean([np.mean(t, axis=0) for t in triangles], axis=0))
            model.enclosure = encloses_point(triangles, probe)
            if model.enclosure >= ENCLOSURE_THRESHOLD:
                mode = "inward"
                auto_note = (f"自動判定: 囲まれている"
                             f"（全方向の {model.enclosure * 100:.0f}% で面に当たる）"
                             f" → 内向きへ揃える / ")
            else:
                mode = "cad"
                auto_note = (f"自動判定: 開いた形状"
                             f"（全方向の {model.enclosure * 100:.0f}% しか面に当たらない）"
                             f" → CAD のまま（モデル側の指定を尊重する） / ")

    if mode == "cad":
        note = "CAD の巻き順をそのまま使う（既定。モデルが法線の正解を持っている前提）"
    elif mode == "flip":
        flip_all = True
        note = "全反転（元コード backtrace.f90 276行 ynnmrev='y' と同じ）"
    elif mode == "shells":
        # シェル（連結成分）単位で「空気側」に揃える。
        # 外殻は内向き、内側にある物体（家具など）は外向きが空気側になる。
        if not model.winding_consistent:
            note = ("巻き順が一貫していないため補正を中止し、CAD のまま使う"
                    "（同じ向きで 2 回現れる辺がある＝隣り合う面の法線が反対を向いている）")
        else:
            per_face_flip = set()
            details = []
            for k, shell in enumerate(model.shells):
                tag = "外殻" if shell["is_outer"] else "内側"
                if not shell["closed"]:
                    details.append(f"シェル{k}({tag}, 開): 補正せず")
                    continue
                want = "inward" if shell["is_outer"] else "outward"
                if shell["normals"] != want:
                    per_face_flip.update(shell["faces"])
                    details.append(f"シェル{k}({tag}): {shell['normals']}→{want} 反転")
                else:
                    details.append(f"シェル{k}({tag}): {shell['normals']} のまま")
            note = "シェル単位で空気側へ揃える / " + " , ".join(details)
    elif mode == "inward":
        # 面ごとにレイの偶奇で「室内側」を判定して揃える。
        # 'shells' と違って**面のつながり（巻き順の一貫性）を要求しない**ので、
        # 床・壁・天井を 1 枚ずつ描いた「板の寄せ集め」モデルでも効く。
        per_face_flip, ambiguous = orient_inward(triangles,
                                                 [f[4] for f in kept])
        note = (f"面ごとにレイの偶奇で室内側へ揃える "
                f"（{len(per_face_flip)} / {len(kept)} 枚を反転）")
        if ambiguous:
            note += (f"  ★{ambiguous} 枚は判定が割れました"
                     f"（隙間のある形状かもしれません。表示で確認してください）")
    else:
        raise ValueError(f"orient_normals に未知の値: {orient_normals!r}"
                         f"（使えるのは 'auto' / 'cad' / 'flip' / 'shells' / 'inward'）")
    note = auto_note + note
    model.orient_mode = mode

    if reference_point is not None:
        print("[read_dxffile] 警告: reference_point は使われません。"
              "法線を音源方向に向ける方式は凸凹の壁や宙に浮いた家具で破綻するため廃止しました。"
              "CAD 側で法線が空気側を向くようにモデルを作ってください。")

    # 自動判定のあとに、**人が目で見て直した指定**を重ねる（normal_editor.py / project.py）。
    # ★`flip_faces` は「**CAD の巻き順から反転する面**」の絶対集合であって、
    #   自動判定への差分ではない。**渡されたら自動判定を丸ごと置き換える。**
    #
    #   normal_editor が保存するのは「人が確認し終えた最終状態」で、
    #   自動判定の結果もそこに畳み込まれている。差分として上に重ねると、
    #   自動が反転した面を手動指定がもう一度反転して**元に戻ってしまう**
    #   （実際にそれで残響時間が半分近く変わった）。
    #   絶対集合として扱えば、保存したものを読めば必ず同じ法線になる。
    manual = None if flip_faces is None else set(int(i) for i in flip_faces)
    if manual is not None:
        note += (f" → 保存済みの指定で置き換え（{len(manual)} / {len(kept)} 枚を反転）")
    model.flipped_faces = set()

    for j, (layer, x1, x2, x3, n) in enumerate(kept):
        if manual is not None:
            flipped = j in manual
        else:
            flipped = flip_all or (per_face_flip is not None and j in per_face_flip)
        if flipped:
            n = -n
            model.flipped_faces.add(j)
        absorption = _resolve_absorption(layer, absorption_table, default_absorption,
                                         band_number, unresolved, model.layer_materials)
        model.mesh.append(Mesh(x1, x2, x3, n, layer, absorption))
        model.layer_counts[layer] = model.layer_counts.get(layer, 0) + 1

    if model.mesh:
        allpts = np.concatenate([m.vertexes for m in model.mesh])
        model.extents = np.array([allpts.min(axis=0), allpts.max(axis=0)])

    # verbose=False は「形状だけ知りたい下読み」なので警告も出さない
    # （面数の照合や受音点の取得で何度も読むため、そのたびに警告が出ると邪魔になる）
    if unresolved and verbose:
        used = 0.1 if default_absorption is None else default_absorption
        print(f"[read_dxffile] 警告: 吸音率が未指定のレイヤ {sorted(unresolved)} → {used} を使用")

    if verbose:
        print("[read_dxffile] " + model.summary().replace("\n", "\n[read_dxffile] "))
        if note:
            print(f"[read_dxffile] 法線の向き: {note}")

    return model


def check_model(model, absorption_table=None, verbose=True):
    """作図ミスを洗い出す（TODO B-10）。

    計算に入る前に「そもそもモデルとして成立しているか」を機械的に確かめる。
    黙って変な結果を出すより、**先に指摘して直してもらう**ほうが早いため。

    見るところ:

    | 項目 | まずい理由 |
    |---|---|
    | 面が 1 枚も無い | DXF のエンティティ種別が非対応かもしれない |
    | 音源／受音点が無い | src / rec レイヤの POINT を描き忘れ |
    | 音源／受音点がモデルの外 | 座標系の取り違え。直接音すら出ない |
    | 吸音率が引けないレイヤ | 既定値 0.1 で計算されるので、結果が意味を持たない |
    | 同じ位置に重なった面 | 音線がどちらに当たるか定まらず、隙間なく反射して見える |
    | 開いた辺が多い | 音線が室外へ逃げる。閉じた室のつもりなら作図ミス |
    | 巻き順が一貫していない | 隣り合う面で法線が反対を向いている |
    | ねじれた四角形 | 三角形への割り方で形が変わる |
    | 極端に小さい面 | 数値誤差で交差判定が不安定になる |

    戻り値: list[dict]  {'level': 'error'|'warning'|'info', 'message': str}
    """
    issues = []

    def add(level, message):
        issues.append({"level": level, "message": message})

    mesh = model.mesh
    if not mesh:
        add("error", "面が 1 枚も読めていません。"
                     f"読み飛ばした種別: {dict(model.skipped['非対応エンティティ']) or 'なし'}")
        return _report_issues(issues, verbose)

    lo, hi = model.extents
    margin = 0.001 * float(np.max(hi - lo))
    for name, points in (("音源", model.source_points), ("受音点", model.receiver_points)):
        if not points:
            add("warning", f"{name}が DXF にありません"
                           f"（src / rec レイヤの POINT。計算時に座標を直接指定すれば動きます）")
            continue
        for p in points:
            if np.any(p < lo - margin) or np.any(p > hi + margin):
                add("error", f"{name} {np.round(p, 3).tolist()} がモデルの外にあります"
                             f"（範囲 {np.round(lo, 3).tolist()}〜{np.round(hi, 3).tolist()}）")

    if absorption_table is not None:
        missing = [layer for layer in model.layer_counts
                   if layer not in model.layer_materials]
        if missing:
            add("warning", f"吸音率が引けないレイヤ {missing} → 既定値 0.1 で計算されます")

    # 同じ位置に重なった面。重心を丸めて突き合わせる
    centres = {}
    for j, face in enumerate(mesh):
        key = tuple(np.round(face.vertexes.mean(axis=0), 6))
        centres.setdefault(key, []).append(j)
    duplicated = [v for v in centres.values() if len(v) > 1]
    if duplicated:
        add("warning", f"同じ位置に重なった面が {len(duplicated)} 組あります"
                       f"（例: 面 {duplicated[0]}）。"
                       f"音線がどちらに当たるか定まりません")

    if not model.is_closed:
        level = "info" if model.open_edges < 4 else "warning"
        add(level, f"開いた辺が {model.open_edges} 本あります"
                   f"（閉じた室のつもりなら作図ミス。一面反射板などなら問題ありません）")
    if not model.winding_consistent:
        add("warning", "巻き順が一貫していません"
                       "（隣り合う面で法線が反対を向いている箇所があります）")

    if model.polygon_notes["ねじれた四角形"]:
        add("warning", f"ねじれた四角形が {model.polygon_notes['ねじれた四角形']} 枚"
                       f"（最大 {model.polygon_notes['最大ねじれ実寸'] * 1000:.1f} mm）。"
                       f"三角形への割り方で形が変わります")
    if model.polygon_notes["分割に失敗"]:
        add("error", f"三角形に分割できなかった多角形が "
                     f"{model.polygon_notes['分割に失敗']} 枚あります")

    areas = np.array([0.5 * np.linalg.norm(np.cross(f.vertexes[1] - f.vertexes[0],
                                                    f.vertexes[2] - f.vertexes[0]))
                      for f in mesh])
    tiny = int(np.count_nonzero(areas < 1.0e-6))
    if tiny:
        add("warning", f"面積が 1 mm² 未満の面が {tiny} 枚あります"
                       f"（交差判定が数値的に不安定になります）")

    skipped = model.skipped["非対応エンティティ"]
    if skipped:
        add("info", f"読み飛ばしたエンティティ: {dict(skipped)}")

    return _report_issues(issues, verbose)


def _report_issues(issues, verbose):
    if verbose:
        marks = {"error": "✗ エラー", "warning": "△ 注意", "info": "・"}
        if not issues:
            print("[check] 作図チェック: 問題は見つかりませんでした")
        else:
            print(f"[check] 作図チェック: {len(issues)} 件")
            for issue in issues:
                print(f"[check]   {marks[issue['level']]} {issue['message']}")
    return issues


def read(file_name, unit=None, absorption_table=None, default_absorption=None,
         orient_normals="cad", reference_point=None, band_number=DEFAULT_BAND_NUMBER, verbose=True):
    """メッシュのリストだけを返す簡易版（procedure.py から使う）。

    音源・受音点や、シェルごとの法線診断も欲しい場合は read_model() を使う。
    """
    return read_model(file_name, unit=unit, absorption_table=absorption_table,
                      default_absorption=default_absorption,
                      orient_normals=orient_normals, reference_point=reference_point,
                      band_number=band_number, verbose=verbose).mesh


if __name__ == "__main__":
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    model = read_model(os.path.join(here, "test.dxf"))
    print()
    for j, m in enumerate(model.mesh[:4]):
        print(f"mesh[{j}] layer='{m.material}' normal={np.round(m.normal, 3).tolist()}")
        print(f"         vertexes={np.round(m.vertexes, 3).tolist()}")
