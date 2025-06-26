import numpy as np

# メッシュクラスでは値のみ保持でメソッドは別作成が良い気がする
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
