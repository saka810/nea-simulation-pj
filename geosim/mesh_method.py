"""音線と三角形の交差判定。**`Mesh` が値だけを持つので、計算はこちらに置く。**

元コード `backtrace.f90` 579〜636 行（音線ループの中の面との判定）。

同じ判定が **2 通り**入っている。**どちらも消さないこと。**

| 実装 | 何を扱うか | 用途 |
|---|---|---|
| `FaceArrays` | 音線の束 × 全面をまとめて配列演算 | **本番**（F-1 高速化。180〜250 倍） |
| `collision_distance` ほか | 音線 1 本 × 面 1 枚 | **参照実装**。元コードの二重ループをそのまま写したもの |

scalar 版は読みやすさと「ベクトル化版が正しいことの基準」のために残してある。
両者が一致することを `tests/test_geosim.py` の
「ベクトル化した交差判定（scalar 版との一致）」で確認している。

| 記号 | 意味 | 元コードの変数 |
|---|---|---|
| `t` | 交点パラメータ。基点から交点までの距離 | `t` |
| `d` | 平面の方程式 `ax+by+cz+d=0` の定数項 | `d` |
| `node` | 音線と面の交点 | `vnode` |
"""


import numpy as np


# ------------------------------------------------------------------------------
# ベクトル化した交差判定（F-1 高速化。2026-08-14 追加）
#
# 下の scalar 版（collision_distance など）は「音線 1 本 × 面 1 枚」を 1 回ずつ扱う。
# 元コードの二重ループをそのまま写したもので読みやすいが、
# 音線 n 本 × 反射 k 回 × 面 m 枚ぶん Python のループが回るので実用速度に届かない。
#
# FaceArrays は面をまとめて配列にしておき、**音線の束 × 全面**を一度の配列演算で処理する。
# 計算している中身は scalar 版とまったく同じ（同じ式・同じ比較・同じ同点の扱い）なので、
# 結果はビット単位で一致する（tests/test_geosim.py で確認）。
# ------------------------------------------------------------------------------

class FaceArrays:
    """メッシュを配列にまとめた入れ物。交差判定を配列演算で行うために使う。

    面ごとに毎回計算していた値（平面の d、三角形の辺ベクトル）を前もって持っておく。

    two_sided
        既定（False）は **面の表側（法線の側）から来た音線しか当たらない**。
        CAD 側で法線が空気側を向いていることを前提にした、元コードと同じ扱い。
        法線が逆を向いた面はすり抜けるので、モデルの誤りが結果に出て気づける。

        True にすると裏からの入射も当てる。CAD で面を 1 枚ずつ描いた「板の寄せ集め」
        モデルは、巻き順や押し出し方向で法線がまちまちになりがちで、
        そのままだと壁が抜ける。反射ベクトル v - 2(v·n)n も
        エネルギー減衰の cosθ = |v·n| も法線の向きに依らないので、
        当たり判定だけ両面にすれば板モデルがそのまま計算できる。
    """

    def __init__(self, mesh, two_sided=False):
        self.two_sided = bool(two_sided)
        self.count = len(mesh)
        if self.count == 0:
            raise ValueError("メッシュが空です")
        self.vertexes = np.array([np.asarray(m.vertexes, dtype=float) for m in mesh])
        self.normal = np.array([np.asarray(m.normal, dtype=float) for m in mesh])

        self.v0 = self.vertexes[:, 0, :]
        self.v1 = self.vertexes[:, 1, :]
        self.v2 = self.vertexes[:, 2, :]

        # 平面の方程式 n·x + d = 0 の d（scalar 版 parameter_d と同じ）
        self.d = -np.einsum("ij,ij->i", self.normal, self.v0)

        # 三角形内部の判定に使う辺ベクトル。
        # scalar 版は innerproduct_from3vertexes(node, 頂点0, 頂点1, 頂点2) と
        #             innerproduct_from3vertexes(node, 頂点1, 頂点2, 頂点0) を呼ぶので、
        # 基準点は頂点0と頂点1、辺はそれぞれ (v1-v0, v2-v0) と (v2-v1, v0-v1)
        self.edge_from_v0_a = self.v1 - self.v0
        self.edge_from_v0_b = self.v2 - self.v0
        self.edge_from_v1_a = self.v2 - self.v1
        self.edge_from_v1_b = self.v0 - self.v1

    def nearest_hit(self, origins, directions, chunk_elements=4_000_000, ignore=None):
        """音線の束について、最も手前で当たる面を求める。

        引数:
            origins    (A,3) 各音線の基点
            directions (A,3) 各音線の方向（単位ベクトル前提）
            chunk_elements   一度に扱う (音線×面) の要素数の上限。メモリ対策
            ignore     (A,) | None  音線ごとに「当てない面」の番号（-1 なら無し）。
                **反射直後は直前に当たった面を必ず渡すこと**（two_sided のとき必須）。
                基点が面の上に乗っているので、丸め誤差で t がごくわずかな正になり、
                同じ面にもう一度当たってしまう（いわゆる self-intersection）。
                片面判定なら反射後の音線は裏側を向くので裏面が捨てられて起きないが、
                両面判定では捨てられないため、その面を明示的に外す必要がある。
                直線は同じ平面と 2 回は交わらないので、外しても取りこぼしは無い。

        戻り値:
            hit_id   (A,)  面のインデックス。当たらなければ -1
            distance (A,)  基点から交点までの距離。当たらなければ inf
            node     (A,3) 交点座標。当たらなければ基点をそのまま返す

        同点（同じ距離の面が複数）のときは**面インデックスの小さいほうを採る**。
        scalar 版が `distance_j < min_distance`（狭義の不等号）で先勝ちにしているのと同じ。
        """
        origins = np.atleast_2d(np.asarray(origins, dtype=float))
        directions = np.atleast_2d(np.asarray(directions, dtype=float))
        n_ray = len(origins)
        if ignore is not None:
            ignore = np.asarray(ignore, dtype=np.int64).reshape(n_ray)

        hit_id = np.full(n_ray, -1, dtype=np.int64)
        distance = np.full(n_ray, np.inf)
        node = origins.copy()

        # (音線 × 面) の要素数が大きくなりすぎないよう音線側を分割する
        step = max(1, int(chunk_elements // self.count))
        for start in range(0, n_ray, step):
            stop = min(start + step, n_ray)
            block_id, block_distance, block_node = self._nearest_hit_block(
                origins[start:stop], directions[start:stop],
                None if ignore is None else ignore[start:stop])
            hit_id[start:stop] = block_id
            distance[start:stop] = block_distance
            node[start:stop] = block_node
        return hit_id, distance, node

    def _nearest_hit_block(self, origins, directions, ignore=None):
        """`nearest_hit` の 1 ブロック分。(hit_id, distance, node) を返す。"""
        n_ray = len(origins)
        hit_id = np.full(n_ray, -1, dtype=np.int64)
        distance = np.full(n_ray, np.inf)
        node = origins.copy()

        # 音線が面の表側から向かっているか（scalar 版 collision_distance の最初の判定）
        denominator = directions @ self.normal.T                   # (A,M)
        numerator = origins @ self.normal.T + self.d[None, :]      # (A,M)
        if self.two_sided:
            candidate = denominator != 0.0     # 面と平行でなければ表裏どちらでも当てる
        else:
            candidate = denominator < 0.0

        # 交点パラメータ t（scalar 版 parameter_t と同じ）。t > 0 の面だけ残す
        with np.errstate(divide="ignore", invalid="ignore"):
            t = -numerator / denominator
        candidate &= t > 0.0

        if ignore is not None:
            has_ignore = ignore >= 0
            if np.any(has_ignore):
                rows = np.nonzero(has_ignore)[0]
                candidate[rows, ignore[rows]] = False

        ray_index, face_index = np.nonzero(candidate)
        if len(ray_index) == 0:
            return hit_id, distance, node

        # ここから先は「候補になった (音線, 面) の組」だけを 1 次元で扱う。
        # (A, M, 3) の配列を作らずに済むのでメモリが小さくなる
        origin_candidate = origins[ray_index]
        node_candidate = (origin_candidate
                          + t[ray_index, face_index][:, None] * directions[ray_index])

        inside = _inside_triangle_batch(
            node_candidate, self.v0[face_index], self.v1[face_index],
            self.edge_from_v0_a[face_index], self.edge_from_v0_b[face_index],
            self.edge_from_v1_a[face_index], self.edge_from_v1_b[face_index])
        if not np.any(inside):
            return hit_id, distance, node

        ray_index = ray_index[inside]
        face_index = face_index[inside]
        node_candidate = node_candidate[inside]
        # scalar 版と同じく交点と基点の距離を取る（|t| ではなく実際のノルム）
        distance_candidate = np.linalg.norm(
            node_candidate - origin_candidate[inside], axis=1)

        # 音線ごとに最短のものを選ぶ。
        # lexsort は安定なので、距離が同じなら元の並び（面インデックス昇順）が保たれる。
        # scalar 版が「先に見つけた面を残す」のと同じ結果になる
        order = np.lexsort((distance_candidate, ray_index))
        ray_sorted = ray_index[order]
        first = np.empty(len(ray_sorted), dtype=bool)
        first[0] = True
        np.not_equal(ray_sorted[1:], ray_sorted[:-1], out=first[1:])

        winner = order[first]
        rays = ray_sorted[first]
        hit_id[rays] = face_index[winner]
        distance[rays] = distance_candidate[winner]
        node[rays] = node_candidate[winner]
        return hit_id, distance, node


def _inside_triangle_batch(node, v0, v1, edge0a, edge0b, edge1a, edge1b):
    """交点が三角形の内部にあるかの判定（ベクトル化）。

    scalar 版 `collision_detection` + `innerproduct_from3vertexes` と同じ計算。
    頂点 0 と頂点 1 のそれぞれを基準に「交点までのベクトル」と「2 辺」の外積を取り、
    その内積が**両方とも 0 以下**なら内部（元コード 625 行 `.le. 0`）。
    """
    to_v0 = node - v0
    cross_a = np.cross(to_v0, edge0a)
    cross_b = np.cross(to_v0, edge0b)
    inner_0 = np.einsum("ij,ij->i", cross_a, cross_b)

    to_v1 = node - v1
    cross_c = np.cross(to_v1, edge1a)
    cross_d = np.cross(to_v1, edge1b)
    inner_1 = np.einsum("ij,ij->i", cross_c, cross_d)

    return (inner_0 <= 0.0) & (inner_1 <= 0.0)


# 7/9打ち合わせ用


# 注意として　音線ベクトルsound_rayは必ずしも音源から出ていない。
# 壁面から反射した音についてsound_rayとしている場合がある（むしろ圧倒的にそっちが多い）
# sound_ray: 元コードvrayに相当
# soundray_comesfrom: 元コードviniに相当
# 一度も反射していない場合　soundray_comesfrom は音源位置と同じ


# 壁に音線が衝突をしたかを判断する
# 衝突判定collisionとその時の距離distanceを返す
# OK
def collision_distance(sound_ray, soundray_comesfrom, normal, vertexes):
    """音線 1 本と面 1 枚の交差判定。**参照実装**（本番は `FaceArrays`）。

    戻り値 `(当たったか, 基点から交点までの距離)`。当たらなければ `(False, 0.0)`。
    `v・n < 0`（法線に向かって進む）ときだけ当たりとするので、**面は片側だけ反射する**。
    """
    collision = False
    distance = 0.0

    # 音線ベクトル reflection_rayとメッシュ法線の内積を計算
    if np.dot(sound_ray, normal) < 0:
        # 平面の方程式ax + by + cz + d = 0のdを算出
        # 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
        # d = -np.dot(normal, vertexes[0])

        # 直線と壁面が交わるときのパラメータtを算出
        t = parameter_t(sound_ray, soundray_comesfrom, normal, vertexes)

        #  基点から音線方向を向いて正の方向にある面を検出
        if t > 0:
            # 交点nodeを計算
            # node = soundray_comesfrom + t * sound_ray
            node = node_renew(sound_ray, soundray_comesfrom, t)

            # 元コード　頂点2を基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
            # なぜ2?
            # 四角形か？
            # 三角形で　頂点　0 1 2があり　0を引き算の後ろと仮定して書いた場合
            # innerproduct_from3vertexes(node,vertex_origin,vertex_1,vertex_2)を使う
            # inner_product_0 = innerproduct_from3vertexes(node, vertexes[0], vertexes[1], vertexes[2])

            # 元コード　頂点3を基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
            # なぜ3?
            # 三角形で　頂点　0 1 2があり　1を引き算の後ろと仮定して書いた場合
            # inner_product_1 = innerproduct_from3vertexes(node, vertexes[1], vertexes[2], vertexes[0])

            # if inner_product_0<0 and inner_product_1 < 0:
            #     collision = True

            collision = collision_detection(node, vertexes)

            distance = np.linalg.norm(soundray_comesfrom - node, ord=2)

    return collision, distance


# 内積を２つ計算し音線と壁が衝突しているかを判定する
# OK
def collision_detection(node, vertexes):
    """交点が三角形の内側にあるか（書籍 式2.50）。**参照実装**。

    2 頂点を基準に外積の内積を取り、**どちらも負なら内側**。
    """
    collision = False
    inner_product_0 = innerproduct_from3vertexes(node, vertexes[0], vertexes[1], vertexes[2])
    inner_product_1 = innerproduct_from3vertexes(node, vertexes[1], vertexes[2], vertexes[0])

    # 2026-08-12 修正: 元コード625行は inpro(1) .le. 0 .and. inpro(2) .le. 0。
    # 旧実装は < 0 で、辺・頂点ちょうどに当たった音線（inpro = 0）を取りこぼしていた。
    if inner_product_0 <= 0 and inner_product_1 <= 0:
        collision = True

    return collision


# 外積を2つ計算し、その外積の内積を計算する
# OK
def innerproduct_from3vertexes(node, vertex_origin, vertex_1, vertex_2):
    """`vertex_origin` を基準に、面を張る 2 辺と交点までのベクトルの外積どうしの内積。

    内側判定の材料。負なら交点がその 2 辺の間にある。
    """
    # 頂点 vertex_originを基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
    # 三角形で　頂点　origin 1 2があり　vertex_originを引き算の後ろと仮定して書いた場合

    # vertexes_origin = np.array([np.zeros(3)], [np.zeros(3)])
    # 書き方変更
    # 2026-08-12 修正: (3, 2) → (2, 3)。3 次元ベクトル 2 本を格納するので行が 2・列が 3。
    # 旧実装のままだと vertexes_toorigin[0] = vertex_1 - vertex_origin で
    # 形状 (3,) を (2,) に代入することになり ValueError。
    vertexes_toorigin = np.zeros((2, 3))

    # node_origin = np.array(np.zeros(3))

    # vertex_originを基点に三角形を考える
    vertexes_toorigin[0] = vertex_1 - vertex_origin
    vertexes_toorigin[1] = vertex_2 - vertex_origin
    node_toorigin = node - vertex_origin

    # 外積2個を計算
    cross_product_0 = np.cross(node_toorigin, vertexes_toorigin[0])
    cross_product_1 = np.cross(node_toorigin, vertexes_toorigin[1])
    # 外積2個から 内積 １つ目 を計算
    inner_product = np.dot(cross_product_0, cross_product_1)

    return inner_product


# 平面の方程式ax + by + cz + d = 0のdを算出
# 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
# OK
def parameter_d(normal, vertex):
    """平面の方程式 `ax+by+cz+d=0` の定数項 `d = -n・x`。頂点はどれでもよい。"""
    # d = -np.dot(normal, vertexes[0])
    d = -np.dot(normal, vertex)
    return d


# 直線と壁面が交わるときのパラメータtを算出
# OK
def parameter_t(sound_ray, soundray_comesfrom, normal, vertexes):
    """音線と平面の交点パラメータ `t = -(n・p0 + d) / (n・v)`。

    `t > 0` なら基点から見て**前方**にある。距離そのもの（音線は単位ベクトル）。
    """
    # 平面の方程式ax + by + cz + d = 0のdを算出
    # 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
    # d = -np.dot(normal, vertexes[0])
    d = parameter_d(normal, vertexes[0])

    # 直線と壁面が交わるときのパラメータtを算出
    t = np.dot(normal, soundray_comesfrom) + d
    t = -t / np.dot(normal, sound_ray)
    return t

# 元コード 388行目
# OK
def node_renew(sound_ray, soundray_comesfrom, t):
    """交点の座標 `q = p0 + t v`。次の区間の基点になる。"""
    # new_node = np.array(np.zeros(3))
    new_node = soundray_comesfrom + t * sound_ray
    return new_node
