import numpy as np
from mesh import Mesh

nref = 1000  # 反射回数（ダミー）
count_noredanduncy = 1000  # 非重複経路の個数（ダミー）

mesh_count = 1000  # メッシュの個数（ダミー）
mesh = [Mesh] * mesh_count  # メッシュの個数分Meshオブジェクトを準備

# 壁面番号
mesh_id = np.zeros((count_noredanduncy, nref))

saoundray_energy = np.zeros(6)
imaginary_sourcepoint = np.zeros((nref, 3))
soundsource_point = np.zeros(3)


# 非重複経路ループ
def loop():
    # 初期条件　虚音源＝音源
    imaginary_sourcepoint[0, :] = soundsource_point

    # 非重複経路分ループ
    for i in range(count_noredanduncy):
        # 反射回数ループ
        for k in range(nref):
            d = mesh[mesh_id[i, k]].parameter_d(self=Mesh(mesh[mesh_id[i, k]]))
            t = np.dot(mesh[mesh_id[i, k]].normal, imaginary_sourcepoint) + d
            t = t / np.linalg.norm(mesh_id[i, k].normal, 2)

            imaginary_sourcepoint[k] = imaginary_sourcepoint[k - 1, :] - 2.0 * t * mesh_id[i, k].normal

    #バックトレースループがここから

    return

