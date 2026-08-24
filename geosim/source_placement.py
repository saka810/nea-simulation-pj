# -*- coding: utf-8 -*-
"""**面の上に置いた音源**（半無響室の床置き音源など）の扱い。

★2026-08-24 ユーザー要望

> 半無響室の検討を行う場合、例えば床面に音源を想定したい。
> CAD 上で平面上に置くことが想定されますが、そのままだとエラーを起こしてしまうかなと。
> 面上・辺上・頂点上に置いた場合、面上の場合は半球に広がっていく様をみたい。
> 逆に、面上に置いている場合、床面の一回目の反射は見ないようにしたい。
> 面上等に置いた場合は、それがどっち方向の音源なのか、等設定項目は必要かなと。
> また、微小にずれていた場合でも、面上とする設定が出来るようにしたい。
> 0.1 m くらいまでは、実際の音源の大きさ含めて、考慮する可能性がある。

**そのままでも落ちはしない**（実測：床ちょうどに置くと音線の半分が
面に当たって消え、受音経路が 113 → 76 本に減る）。減るだけで気づきにくいので、
ここで**面の上に置いたと分かったら半球に飛ばす**ようにする。

決めごと:
  ・**半球にするのは 1 面ぶん**（放射方向を決める面。「載っている面」が
    複数あっても 1/4・1/8 にはしない）。半無響室の壁は吸音材なので、
    壁の方へ飛んだ音線は**ふつうに壁で吸われる**のが正しい。
    1/8 にしてしまうと、その吸収を「最初から出さない」ことにすり替わる
  ・**載っている面の 1 次反射は起きない**。音源を面から `EPSILON`（1 mm）だけ
    浮かせ、その面から離れる向きの半球にしか飛ばさないので、
    その面へ戻るには**別の面で反射してから**しかない（要望どおり）
  ・**エネルギーは全部を半球に入れる**（音線 1 本あたり E/N のまま）。
    床置きの音源が持つ音響パワー W が 2π に出る、という素直な形。
    自由音場に比べて +3 dB になるのが正しい姿
  ・載っている面が複数（辺上・頂点上）のときは、**放射方向の面以外からは
    `EPSILON` だけ離す**。離しておかないと、その面に向かう音線が
    「室の外」へ出て消えてしまい、吸音率どおりに吸われない
"""
import numpy as np

# 面から浮かせる隙間 [m]。1 mm。実寸の音源も面にめり込んではいないので実害はなく、
# 「面のちょうど上」に置いたときの数値の揺れ（当たり判定が 0 距離になる）を避けられる
EPSILON = 1.0e-3

# 面の上とみなす既定の距離 [m]。0 なら「ぴったり載っているときだけ」
DEFAULT_TOLERANCE = 0.0

# 面に載っているかを見るときの、面の縁のはみ出し許容（面の内側判定を少し緩める）
EDGE_MARGIN = 1.0e-6

# 放射方向の指定（設定画面の選択肢）。None は自動
DIRECTIONS = {
    "自動": None,
    "+X": (1.0, 0.0, 0.0), "-X": (-1.0, 0.0, 0.0),
    "+Y": (0.0, 1.0, 0.0), "-Y": (0.0, -1.0, 0.0),
    "+Z": (0.0, 0.0, 1.0), "-Z": (0.0, 0.0, -1.0),
}

KIND_FREE = "自由空間"
KIND_FACE = "面上"
KIND_EDGE = "辺上"
KIND_VERTEX = "頂点上"


class Placement:
    """音源の置かれ方。

    属性:
        point     計算に使う音源位置（面から `EPSILON` だけ浮かせたもの）
        normal    半球の軸（面から室へ向かう単位ベクトル）。None なら全球
        kind      `面上` / `辺上` / `頂点上` / `自由空間`
        faces     載っていると判定した面の番号
        mount     放射方向を決めた面の番号（この面の 1 次反射は起こらない）
        distance  もとの点から `mount` の面までの距離 [m]
        planes    載っている面の数（向きの違う平面の数）
    """

    def __init__(self, point, normal=None, kind=KIND_FREE, faces=(), mount=None,
                 distance=0.0, planes=0, original=None):
        self.point = np.asarray(point, dtype=float)
        self.normal = None if normal is None else np.asarray(normal, dtype=float)
        self.kind = kind
        self.faces = list(faces)
        self.mount = mount
        self.distance = float(distance)
        self.planes = int(planes)
        self.original = (self.point if original is None
                         else np.asarray(original, dtype=float))

    @property
    def on_surface(self):
        return self.normal is not None

    @property
    def solid_angle(self):
        """音が出ていく立体角 [sr]（半球なら 2π）。出力の説明に使う。"""
        return 2.0 * np.pi if self.on_surface else 4.0 * np.pi

    def describe(self):
        if not self.on_surface:
            return "自由空間（全球に放射）"
        moved = float(np.linalg.norm(self.point - self.original))
        return (f"{self.kind}（半球に放射）／放射方向 "
                f"[{self.normal[0]:+.2f}, {self.normal[1]:+.2f}, {self.normal[2]:+.2f}]"
                f"／面までの距離 {self.distance * 1000:.1f} mm"
                f"／載っている平面 {self.planes} 枚"
                f"／音源を {moved * 1000:.1f} mm 動かしました")

    def summary_row(self):
        """測定点の表に足す 1 行ぶんの情報。"""
        return {"置かれ方": self.kind,
                "放射方向": ("" if self.normal is None
                             else f"[{self.normal[0]:+.3f}, {self.normal[1]:+.3f}, "
                                  f"{self.normal[2]:+.3f}]"),
                "面までの距離 [m]": round(self.distance, 4)}


def _plane_distance(point, triangle, normal):
    """点から三角形の**平面**までの符号なし距離と、面の内側に落ちるか。"""
    vertexes = np.asarray(triangle, dtype=float)
    normal = np.asarray(normal, dtype=float)
    offset = float(np.dot(point - vertexes[0], normal))
    projected = point - offset * normal
    # 面の内側か（重心座標）
    v0 = vertexes[1] - vertexes[0]
    v1 = vertexes[2] - vertexes[0]
    v2 = projected - vertexes[0]
    d00, d01, d11 = np.dot(v0, v0), np.dot(v0, v1), np.dot(v1, v1)
    d20, d21 = np.dot(v2, v0), np.dot(v2, v1)
    denominator = d00 * d11 - d01 * d01
    if abs(denominator) < 1.0e-18:
        return abs(offset), False
    v = (d11 * d20 - d01 * d21) / denominator
    w = (d00 * d21 - d01 * d20) / denominator
    inside = (v >= -EDGE_MARGIN and w >= -EDGE_MARGIN
              and v + w <= 1.0 + EDGE_MARGIN)
    return abs(offset), inside


def _inward(point, normal, mesh, faces):
    """面から**室の側**へ向かう単位ベクトル。

    モデルの法線が内向きに揃っていればそのまま使えるが、揃っていない
    モデルもあるので、**面の重心の平均**（室のだいたいの中心）を見て決める。
    """
    normal = np.asarray(normal, dtype=float)
    normal = normal / (np.linalg.norm(normal) or 1.0)
    centre = np.mean([np.mean(np.asarray(m.vertexes, dtype=float), axis=0)
                      for m in mesh], axis=0)
    if np.dot(centre - point, normal) < 0:
        return -normal
    return normal


def detect(point, mesh, tolerance=DEFAULT_TOLERANCE, direction=None,
           enabled=True):
    """音源が面の上に載っているかを調べる。→ `Placement`

    引数:
        point      CAD に描かれた音源の位置
        mesh       室の三角形（`Mesh`）
        tolerance  これ以下の距離なら「面の上」とみなす [m]。
                   ★実寸の音源の大きさを見込んで 0.1 m くらいまで使う想定
        direction  放射方向の指定（`DIRECTIONS` の名前、または単位ベクトル）。
                   `None` / `自動` なら、いちばん近い面から決める
        enabled    False なら何もしない（従来どおり全球に飛ばす）
    """
    point = np.asarray(point, dtype=float)
    if not enabled or not mesh:
        return Placement(point)

    limit = max(float(tolerance), EPSILON)
    found = []
    for index, face in enumerate(mesh):
        distance, inside = _plane_distance(point, face.vertexes, face.normal)
        if inside and distance <= limit:
            found.append((distance, index, np.asarray(face.normal, dtype=float)))
    if not found:
        return Placement(point)

    # 向きの違う平面をまとめる（同じ平面の三角形はひとかたまり）
    planes = []
    for distance, index, normal in found:
        unit = normal / (np.linalg.norm(normal) or 1.0)
        for plane in planes:
            if abs(abs(float(np.dot(unit, plane["normal"]))) - 1.0) < 1.0e-6:
                plane["faces"].append(index)
                plane["distance"] = min(plane["distance"], distance)
                break
        else:
            planes.append({"normal": unit, "faces": [index], "distance": distance})

    wanted = DIRECTIONS.get(direction, direction) if direction is not None else None
    if wanted is not None:
        wanted = np.asarray(wanted, dtype=float)
        wanted = wanted / (np.linalg.norm(wanted) or 1.0)

    def score(plane):
        """放射方向にいちばん合う面を選ぶ（指定が無ければ近い順）。"""
        inward = _inward(point, plane["normal"], mesh, plane["faces"])
        if wanted is not None:
            return (-float(np.dot(inward, wanted)), plane["distance"])
        # 自動：いちばん近い面。同じくらいなら**上向き**を優先（床置きが普通）
        return (round(plane["distance"], 6), -float(inward[2]))

    mount = min(planes, key=score)
    normal = _inward(point, mount["normal"], mesh, mount["faces"])

    # 放射方向の面には載せ、そのほかの面からは少しだけ離す
    vertex = np.asarray(mesh[mount["faces"][0]].vertexes, dtype=float)[0]
    offset = float(np.dot(point - vertex, normal))
    placed = point - offset * normal + EPSILON * normal
    for plane in planes:
        if plane is mount:
            continue
        other = _inward(point, plane["normal"], mesh, plane["faces"])
        vertex = np.asarray(mesh[plane["faces"][0]].vertexes, dtype=float)[0]
        gap = float(np.dot(placed - vertex, other))
        if gap < EPSILON:
            placed = placed + (EPSILON - gap) * other

    kind = {1: KIND_FACE, 2: KIND_EDGE}.get(len(planes), KIND_VERTEX)
    faces = [index for plane in planes for index in plane["faces"]]
    return Placement(placed, normal=normal, kind=kind, faces=faces,
                     mount=mount["faces"][0], distance=mount["distance"],
                     planes=len(planes), original=point)


def hemisphere(rays, normal):
    """全球の音線を**半球**に折り返す。

    `soundray_generator` が作る球面上の等分布を、面の裏へ行く分だけ
    **鏡映**して表側へ移す。鏡映は 1 対 1 なので、本数も分布の均さも保たれる
    （半球に N 本、密度は 2 倍）。
    """
    rays = np.asarray(rays, dtype=float)
    normal = np.asarray(normal, dtype=float)
    normal = normal / (np.linalg.norm(normal) or 1.0)
    projection = rays @ normal
    behind = projection < 0.0
    if not np.any(behind):
        return rays
    folded = rays.copy()
    folded[behind] -= 2.0 * projection[behind, None] * normal[None, :]
    return folded
