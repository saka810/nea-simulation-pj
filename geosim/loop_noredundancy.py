import numpy as np
from mesh import Mesh
import sound_ray as sr
import mesh_method as mm

nref = 1000  # 反射回数（ダミー）
count_noredanduncy = 1000  # 非重複経路の個数（ダミー）

mesh_count = 1000  # メッシュの個数（ダミー）
mesh = []  # メッシュの個数分Meshオブジェクトを準備
for i in range(mesh_count):
    mesh.append(Mesh(np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3), "name", 0.1))

# 壁面番号
mesh_id_list = np.zeros((count_noredanduncy, nref))
mesh_id = 0

saoundray_energy = np.zeros(6)
imaginary_sourcepoint = np.zeros((nref, 3))
soundsource_point = np.zeros(3)

sound_ray = np.zeros(3)
soundray_comesfrom = np.zeros(3)
reciever_point = np.zeros(3)
node = np.zeros(3)


# collision = False

mesh_nearestid=0 #最も手前の交点を持つ面のインデックス
temp_distance=1000000000.0 #とりあえず大きい数字


# 非重複経路ループ
def loop():
    # 初期条件　虚音源＝音源
    imaginary_sourcepoint[0, :] = soundsource_point

    # 非重複経路分ループ
    for i in range(count_noredanduncy):
        # 反射回数ループ
        for k in range(nref):
            mesh_id = mesh_id_list[i, k]
            d = mm.parameter_d(mesh[mesh_id].normal, mesh[mesh_id].vertexes)
            t = np.dot(mesh[mesh_id].normal, imaginary_sourcepoint) + d
            t = t / np.linalg.norm(mesh[mesh_id].normal, 2)

            imaginary_sourcepoint[k] = imaginary_sourcepoint[k - 1, :] - 2.0 * t * mesh[mesh_id[i, k]].normal

        # バックトレースループがここから
        sound_ray = imaginary_sourcepoint - reciever_point
        sound_ray = sr.noramlized_soundray(sound_ray)
        soundray_comesfrom = reciever_point

        for k in range(nref, 0, -1):
            for j in range(mesh_count):
                collision = False
                t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[j].normal, mesh[j].vertexes)

                if t > 0:
                    node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                    collision, distance = mm.collision_distance(sound_ray, soundray_comesfrom, mesh[mesh_id].normal,
                                                                mesh[mesh_id].vertexes)

                    if distance<temp_distance:
                        temp_distance=distance
                        mesh_nearestid=mesh_id

                if collision:
                    t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[mesh_id].normal, mesh[mesh_id].vertexes)
                    node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                    # ここから元コード1086に移動

                    if mesh_nearestid==mesh_id:
                        # 　エネルギーを減衰させる式がここに入る 今は書いてない
                        soundray_comesfrom=sr.soundraycomesfrom_renew(node)
                        sound_ray=sr.soundray_renew(imaginary_sourcepoint[k-1],soundray_comesfrom)
                        sound_ray=sr.noramlized_soundray(sound_ray)

        # k==0となる部分　元コード1066はfor kの外に出す
        if distance > np.linalg.norm(soundsource_point - soundray_comesfrom, ord=2):
            # 受音リストの書き込み　csvファイルとかに1行ずつ書く　以下はダミー
            print("受音リストに書き込み")
            # 元コードは<でif にしてループを抜けていたが外に出すなら＞の条件だけで十分？









    return
