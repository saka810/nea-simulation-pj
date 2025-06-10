import numpy as np

###
###
#メッシュクラスでは値のみ保持でメソッドは別作成が良い気がする

# メッシュは三角形メッシュに限定
# 他を考慮すると一向にプロジェクトが進行しないため
class Mesh:
    # オクターブバンド中心周波数　meshオブジェクトに持たせるかは要検討
    frequencies = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

    # meshオブジェクトは以下を持つ

    # メッシュの識別番号
    # id_number: int = 0

    # 三角形メッシュの頂点座標 1行目からそれぞれ頂点1~3の(x,y,z)座標
    # vertexes = np.array(([np.zeros(3)], [np.zeros(3)], [np.zeros(3)]))

    # メッシュの法線
    # normal = np.array(np.zeros(3))

    # メッシュの吸音率情報
    # 吸音材の名前
    # material: str = "material"
    # 吸音率
    # absorption_coefficient = np.zeros(frequencies.size)

    def __init__(self, vertex_1, vertex_2, vertex_3, normal, material, absorption_coefficient):
        self.vertexes = np.array(([vertex_1], [vertex_2], [vertex_3]))
        self.normal = normal
        self.material = material
        self.absorption_coefficient = absorption_coefficient

    # 注意として　音線ベクトルsound_rayは必ずしも音源から出ていない。
    # 壁面から反射した音についてsound_rayとしている場合がある（むしろそっちが多い）
    # sound_ray: vrayに相当
    # soundray_comesfrom: viniに相当
    # 一度も反射していない場合　soundray_comesfromは音源位置

    def collision_distance(self, sound_ray, soundray_comesfrom):
        collision = False
        distance = 0.0

        # 音線ベクトル reflection_rayとメッシュ法線の内積を計算
        if np.dot(sound_ray, self.normal) < 0:
            # 平面の方程式ax + by + cz + d = 0のdを算出
            # 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
            # d = -np.dot(normal, vertexes[0])

            # 直線と壁面が交わるときのパラメータtを算出
            # t = sound_ray[0] * soundray_comesfrom[0] \
            #     + sound_ray[1] * soundray_comesfrom[1] \
            #     + sound_ray[2] * soundray_comesfrom[2] * d
            # t = -t / np.dot(sound_ray, soundray_comesfrom)

            t = self.parameter_t(sound_ray, soundray_comesfrom)

            #  基点から音線方向を向いて正の方向にある面を検出
            if t > 0:
                # 交点nodeを計算
                # node = soundray_comesfrom + t * sound_ray
                node = self.node_renew(sound_ray, soundray_comesfrom, t)

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

                collision = self.collision_detection(node, self.vertexes)

                distance = np.linalg.norm(soundray_comesfrom - node, ord=2)

        return collision, distance

    def collision_detection(self, node, vertexes):
        collision = False
        inner_product_0 = self.innerproduct_from3vertexes(node, vertexes[0], vertexes[1], vertexes[2])
        inner_product_1 = self.innerproduct_from3vertexes(node, vertexes[1], vertexes[2], vertexes[0])

        if inner_product_0 < 0 and inner_product_1 < 0:
            collision = True

        return collision

    def innerproduct_from3vertexes(self, node, vertex_origin, vertex_1, vertex_2):
        # 頂点 vertex_originを基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
        # 三角形で　頂点　origin 1 2があり　vertex_originを引き算の後ろと仮定して書いた場合

        vertexes_origin = np.array([np.zeros(3)], [np.zeros(3)])
        node_origin = np.array([np.zeros(3)])

        vertexes_origin[0] = vertex_1 - vertex_origin
        vertexes_origin[1] = vertex_2 - vertex_origin
        node_origin = node - vertex_origin

        # 外積2個を計算
        cross_product_0 = np.cross(node_origin, vertexes_origin[0])
        cross_product_1 = np.cross(node_origin, vertexes_origin[1])
        # 外積2個から 内積 １つ目 を計算
        inner_product = np.dot(cross_product_0, cross_product_1)

        return inner_product

    def parameter_d(self):
        d = -np.dot(self.normal, self.vertexes[0])
        return d

    def parameter_t(self, sound_ray, soundray_comesfrom):
        # 平面の方程式ax + by + cz + d = 0のdを算出
        # 頂点一つと法線の掛け算であっている。 頂点は任意で良い。
        # d = -np.dot(normal, vertexes[0])
        d = self.parameter_d()

        # 直線と壁面が交わるときのパラメータtを算出
        t = sound_ray[0] * soundray_comesfrom[0] \
            + sound_ray[1] * soundray_comesfrom[1] \
            + sound_ray[2] * soundray_comesfrom[2] * d
        t = -t / np.dot(sound_ray, soundray_comesfrom)
        return t

    def node_renew(self, sound_ray, soundray_comesfrom, t):
        new_node = np.array(np.zeros(3))
        new_node = soundray_comesfrom + t * sound_ray
        return new_node
