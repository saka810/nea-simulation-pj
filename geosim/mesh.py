"""三角形メッシュ 1 枚を表すクラス。**値だけを持ち、計算は `mesh_method.py` に置く。**

元コード `backtrace.f90` が配列で持っていた面の情報（頂点・法線・吸音率）を
1 つのオブジェクトにまとめたもの。面の並びは `read_dxffile.read_model()` が作る。

| このクラスの属性 | 意味 | 元コードの配列 |
|---|---|---|
| `vertexes` (3,3) | 三角形の頂点。行が頂点 1〜3、列が x/y/z | `pos(3,3,nmesh)` |
| `normal` (3,) | 面の法線（単位ベクトル）。**音が通る空気側を向く** | `vn(3,nmesh)` |
| `material` str | 吸音材の名前（＝DXF のレイヤ名） | — |
| `absorption_coefficient` (nband,) | バンド別の**垂直入射**吸音率 | `absper(6,nmesh)` |

**メッシュは三角形に限定**している（他を考慮するとプロジェクトが進まないため）。
4 角形以上の面は `read_dxffile` が三角形に分割してから渡す。
"""


import numpy as np


# 7/9打ち合わせ用


# メッシュクラスでは値のみ保持でメソッドは別作成が良い気がする
# メッシュは三角形メッシュに限定
# 他を考慮すると一向にプロジェクトが進行しないため
class Mesh:
    """三角形メッシュ 1 枚。**値だけを持つ**（計算は `mesh_method.py`）。

    属性は `vertexes` (3,3) / `normal` (3,) / `material` / `absorption_coefficient` (nband,)。
    意味と元コードとの対応はこのファイルの冒頭を参照。
    """
    # オクターブバンド中心周波数　meshオブジェクトに持たせるかは要検討
    # frequencies = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

    # meshオブジェクトは以下を持つ

    # メッシュの識別番号　多分いらない
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
        # 2026-08-12 修正: ([v1], [v2], [v3]) だと形状が (3, 1, 3) になっていた。
        # 意図は「頂点3点 × xyz」の (3, 3)。mesh_method 側は vertexes[0] が (3,) 前提。
        self.vertexes = np.array([vertex_1, vertex_2, vertex_3])
        self.normal = normal
        # normalはcadにデータがあればそのまま代入　なければ何か計算
        self.material = material
        self.absorption_coefficient = absorption_coefficient
