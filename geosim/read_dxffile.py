import numpy as np
from mesh import Mesh

# 目的 meshのリストを返すこと
# Meshオブジェクトに
# 各三角形メッシュの頂点座標(x,y,z)を合計3点分と、法線ベクトル
# 材料名と材料の吸音率を保持させる
# そのMeshオブジェクトをmeshにメッシュの個数分追加する
def read(file_name):
    # このflag、countは特定の部分のみ処理を行うためのもの
    # POLYLINEの直後数行のみ読み取る想定で記述しています
    flag = False
    count = 0
    mesh = []
    vertex_1 = 0.0
    vertex_2 = 0.0
    vertex_3 = 0.0

    with open(file_name, encoding='utf-8') as f:
    # with open(r'C:\Users\yatsunami\Documents\Pycharm\dxfread_test_20240919\Cuboid.dxf', encoding='utf-8') as f:
        for line in f:
            # if "AcDbPolyFaceMesh" in line:
            #     if "Vertex" in line:
            #         pass
            #     else:
            #         print(line)
            if not flag:
                # if "AcDbPolyFaceMesh" in line:
                #     if not "Vertex" in line:
                #         flag = True
                if "POLYLINE" in line:
                    flag = True
            if flag:
                # print(line.strip())
                # ここだけ読んだら頂点の番号だけ？　座標はどこに
                if count > 1 and count < 5:
                    print(line.strip())
                count += 1

            # 多分10進数変換がいる　座標点とかはどこに？？
            # 下の感じで情報を集めたあとにmeshに入れたい
            if count == 2:
                vertex_1 = line.strip()
            elif count == 3:
                vertex_2 = line.strip()
            elif count == 4:
                vertex_3 = line.strip()
            elif count > 5:
                count = 0
                flag = False

    return mesh


if __name__ == "__main__":
    read(r'C:\Users\yatsunami\Documents\Pycharm\dxfread_test_20240919\Cuboid.dxf')
