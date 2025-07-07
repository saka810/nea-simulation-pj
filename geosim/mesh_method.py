import numpy as np


# 7/9打ち合わせ用



# 注意として　音線ベクトルsound_rayは必ずしも音源から出ていない。
# 壁面から反射した音についてsound_rayとしている場合がある（むしろ圧倒的にそっちが多い）
# sound_ray: 元コードvrayに相当
# soundray_comesfrom: 元コードviniに相当
# 一度も反射していない場合　soundray_comesfrom は音源位置と同じ


# 壁に音線が衝突をしたかを判断する
# 衝突判定collisionとその時の距離distanceを返す
def collision_distance(sound_ray, soundray_comesfrom, normal, vertexes):
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
def collision_detection(node, vertexes):
    collision = False
    inner_product_0 = innerproduct_from3vertexes(node, vertexes[0], vertexes[1], vertexes[2])
    inner_product_1 = innerproduct_from3vertexes(node, vertexes[1], vertexes[2], vertexes[0])

    if inner_product_0 < 0 and inner_product_1 < 0:
        collision = True

    return collision


# 外積を2つ計算し、その外積の内積を計算する
def innerproduct_from3vertexes(node, vertex_origin, vertex_1, vertex_2):
    # 頂点 vertex_originを基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
    # 三角形で　頂点　origin 1 2があり　vertex_originを引き算の後ろと仮定して書いた場合

    # vertexes_origin = np.array([np.zeros(3)], [np.zeros(3)])
    # 書き方変更
    vertexes_toorigin = np.array(np.zeros((3, 2)))

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
def parameter_d(normal, vertex):
    # d = -np.dot(normal, vertexes[0])
    d = -np.dot(normal, vertex)
    return d


# 直線と壁面が交わるときのパラメータtを算出
def parameter_t(sound_ray, soundray_comesfrom, normal, vertexes):
    # 平面の方程式ax + by + cz + d = 0のdを算出
    # 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
    # d = -np.dot(normal, vertexes[0])
    d = parameter_d(normal, vertexes[0])

    # 直線と壁面が交わるときのパラメータtを算出
    t = sound_ray[0] * soundray_comesfrom[0] \
        + sound_ray[1] * soundray_comesfrom[1] \
        + sound_ray[2] * soundray_comesfrom[2] * d
    t = -t / np.dot(sound_ray, soundray_comesfrom)
    return t


def node_renew(sound_ray, soundray_comesfrom, t):
    new_node = np.array(np.zeros(3))
    new_node = soundray_comesfrom + t * sound_ray
    return new_node
